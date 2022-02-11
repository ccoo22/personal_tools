# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use Parallel::ForkManager;
$|=1;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $blastn      = "/home/genesky/software/blast+/2.9.0/bin/blastn";
my $lineages    = "/home/genesky/database/ncbi/nt_20191128/ncbi_lineages_2020-05-29.csv";
my $nt_db       = "/home/genesky/database/ncbi/nt/blast_idx/nt";
my $table2excel = "perl /home/pub/bin/NGS/chip/GATK4/tools/personal/table2excel.pl";

# 检测 -> 脚本输入
my ($fastq_dir, $output_dir, $sample_list, $reads_count, $thread, $word_size, $evalue, $perc_identity, $max_target_seqs, $max_hsps, $seqidlist, $if_help);
GetOptions(
    "fastq_dir|i=s"     => \$fastq_dir,
    "output_dir|o=s"    => \$output_dir,
    "sample_list|s=s"   => \$sample_list,
    "reads_count|c=s"   => \$reads_count,
    "thread|t:i"        => \$thread,
    "word_size|w:i"     => \$word_size,
    "evalue|e:s"        => \$evalue,
    "perc_identity|p:i" => \$perc_identity,
    "max_target_seqs|m:i" => \$max_target_seqs,
    "max_hsps:i"          => \$max_hsps,
    "seqidlist=s"       => \$seqidlist,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $fastq_dir or not defined $output_dir));
###################################################################### 主程序
$sample_list     = ""    if(not defined $sample_list);
$reads_count     = 10000 if(not defined $reads_count);
$thread          = 3     if(not defined $thread);
$word_size       = (not defined $word_size)       ? ""                     : "-word_size $word_size";
$evalue          = (not defined $evalue)          ? "-evalue 1e-5"         : "-evalue $evalue";
$perc_identity   = (not defined $perc_identity)   ? "-perc_identity 97"    : "-perc_identity $perc_identity";
$max_target_seqs = (not defined $max_target_seqs) ? "-max_target_seqs 200" : "-max_target_seqs $max_target_seqs";
$max_hsps        = (not defined $max_hsps)        ? "-max_hsps 1"          : "-max_hsps $max_hsps";
$seqidlist       = (not defined $seqidlist)       ? ""                     : "-seqidlist $seqidlist";

# 这个软件的功能是做物种筛查，故应该尽可能的提高物种序列数量，减少同一个物种的多次输出

#########
# (1) 结果目录创建
#########
print "1: create output dir\n";
make_dir($output_dir);

#########
# (2) 样本获取
#########
print "2: get sample name\n";
my @samples;
if($sample_list eq "") # 所有样本
{
    opendir FASTQ_DIR, $fastq_dir;
    while(my $file = readdir FASTQ_DIR)
    { 
        my ($sample) = $file =~ /(.*)_R1.fastq.gz/;
        push @samples, $sample if(defined $sample);
    }
    close FASTQ_DIR;
}
else # 自定义样本
{
    @samples = split /,/, $sample_list;
}
 
print "    process sample: " . (join ",", @samples) . "\n";

#########
# (3) blastn , anno species, summary
#########
print "3: blastn[10 minutes/sample] , anno species, summary\n";

# 临时结果目录
my $tmp_dir = "$output_dir/tmp";
make_dir($tmp_dir);

# 样本并行
my %hashTax2Species = read_lineage($lineages, 1, 2, 3, 4, 5, 6, 7, 12); # 读取物种注释 
my $pm = Parallel::ForkManager->new($thread);
   $pm->run_on_start(sub{my ($pid, $sample) = @_; process_bar_array($sample, \@samples)});# 进度条
foreach my $sample(@samples)
{
    $pm->start($sample) and next;
    my $fastq_r1   = "$fastq_dir/$sample\_R1.fastq.gz";
    my $fasta      = "$tmp_dir/$sample\_R1.$reads_count.fa";
    my $blast      = "$tmp_dir/$sample\_R1.$reads_count.blast";
    my $summary    = "$tmp_dir/$sample.summary.txt";
    
    fastq_to_fasta($fastq_r1, $fasta, $reads_count); # 提取指定数量的fastq序列
    system("$blastn -query $fasta -db $nt_db -out $blast -outfmt '7 qacc sseqid staxids pident qcovs length' $evalue $perc_identity $word_size $max_target_seqs $max_hsps -num_threads 20 $seqidlist   -dust no"); # 比对  qacc sseqid staxids pident qcovs length qstart qend sstart send evalue salltitles
    summary($blast, $fasta, $summary, \%hashTax2Species, $reads_count); # 注释species并汇总

    $pm->finish;    
}
$pm->wait_all_children;

#########
# (5) summary to excel
#########
print "5: summary to excel\n";

my $excel_summary        = "$output_dir/NT_Blast_Summary.xlsx";        # 每一条序列只取最优的一条比对结果
my $excel_sum_all        = "$output_dir/NT_Blast_Summary_SumAll.xlsx"; # 取所有比对结果，看每一个物种的比对情况。如果想充分利用这个结果，最好解除max_target_seqs参数限制。
my $excel_sum_decrease   = "$output_dir/NT_Blast_Summary_SumDecrease.xlsx"; # 在SumAll的基础上，依次挑选丰度最高物种，输出统计，并删除该物种包含的reads记录，剩下的reads再重新挑选丰度最高的物种...
my $excel_unmapped       = "$output_dir/NT_Blast_unmapped.xlsx";       # 比对失败的结果
my $summary_list       = "";
my $sum_all_list       = "";
my $sum_decrease_list  = "";
my $unmapped_list      = "";
my $sample_list2       = "";
foreach my $sample(@samples)
{
    my $summary         = "$tmp_dir/$sample.summary.txt";
    my $sum_all         = "$summary.sum.all.txt";
    my $sum_decrease    = "$summary.sum.decrease.txt";
    my $unmapped        = "$summary.unmapped.txt";
    
    $summary_list       .= "$summary,";
    $sum_all_list       .= "$sum_all,";
    $sum_decrease_list  .= "$sum_decrease,";
    $unmapped_list      .= "$unmapped,";
    $sample_list2       .= "$sample,";

}
system("$table2excel -i $summary_list       -s $sample_list2 -o $excel_summary");
system("$table2excel -i $sum_all_list       -s $sample_list2 -o $excel_sum_all");
system("$table2excel -i $sum_decrease_list  -s $sample_list2 -o $excel_sum_decrease");
system("$table2excel -i $unmapped_list      -s $sample_list2 -o $excel_unmapped");
print "final:           $excel_summary  \n";
print "accumlate stats: $excel_sum_all  \n";
print "decrease  stats: $excel_sum_decrease\n";
print "unmapped:        $excel_unmapped \n";



###################################################################### 子函数

# 注释比对结果
sub summary{
    my $blast           = shift @_; # 比对结果
    my $fasta           = shift @_; # 原始fasta序列
    my $summary         = shift @_; # 最终汇总结果
    my $hashTax2Species = shift @_; # tax -> species 映射
    my $reads_count     = shift @_; # 项目分析总reads数量
    
    my @need_annos = split /,/, "superkingdom,kingdom,phylum,class,order,family,genus";

    #####################
    # (1) 注释汇总
    #####################
    my %hashSUMALL;  # 记录一条序列的所有比对结果，汇总统计模式
    my %hashCount;  # 记录最优比对结果统计
    my %hashBlastReads;  # 记录比对成功的reads名称,物种
    my %hashAnno;  # 记录每一个物种的注释信息
    my $count_blast = 0; # 比对成功的reads数量
    my $last_reads  = "";  # 上一条reads的名称
    open BLAST, $blast;
    while(<BLAST>)
    {
        next if($_=~/^#/);
        $_=~s/[\r\n]//g;
        my ($reads_name, $gi_info, $tax_id, $tmp) = split /\t/, $_, 4;
        my ($gi, $gi_id, $gb, $gb_id) = split /[|]/, $gi_info;

        # 关键物种
        my $species      = (exists $hashTax2Species{$tax_id}{'species'}      and $hashTax2Species{$tax_id}{'species'}      =~ /\w/) ? $hashTax2Species{$tax_id}{'species'}      : "NO_RANK";

        # 每一个物种的比对reads数量
        $hashBlastReads{$reads_name}{$species}++;
        $hashSUMALL{$species}{$reads_name}++;
        $hashSUMALL{$species}{'UniqueCount'}++ if($hashSUMALL{$species}{$reads_name} == 1);
        foreach my $level(@need_annos)
        {
            my $anno_info = (exists $hashTax2Species{$tax_id}{$level} and $hashTax2Species{$tax_id}{$level} =~ /\w/) ? $hashTax2Species{$tax_id}{$level} : "NO_RANK";
            $hashAnno{$species}{$level} = $anno_info;
        }


        # 只看比对的第一个结果
        next if($reads_name eq $last_reads);  # 后面的物种注释，仅读取第一个比对结果
        $last_reads = $reads_name;
        

        # 种水平的注释记录
        $hashCount{'Species'}{$species}{'Count'}++;
        $hashCount{'Species'}{$species}{"GB"}{$gb_id}++;  # 记录物种GB信息

        # 物种其他注释
        foreach my $level(@need_annos)
        {
            my $anno_info = (exists $hashTax2Species{$tax_id}{$level} and $hashTax2Species{$tax_id}{$level} =~ /\w/) ? $hashTax2Species{$tax_id}{$level} : "NO_RANK";
            # $hashCount{'Species'}{$species}{$level}{$anno_info}++; 
            $hashCount{'Species'}{$species}{$level} = $anno_info; # 记录当前物种在每个水平上的注释信息。
            $hashCount{'Level'}{$level}{$anno_info}++; # 记录每一个水平下，每一种注释的数量  
        }

        $count_blast++;
    }
    close BLAST;
    $hashCount{'Species'}{'NO_BLAST'}{'Count'} = $reads_count - $count_blast if($count_blast < $reads_count); # 部分序列没有比对上

    #####################
    # (2) 汇总结果输出
    #####################
    open SUMMARY, ">$summary";
    print SUMMARY "Species\tReadsCount\tPerc\tGB\t" . (join "\t", @need_annos) . "\n";
    foreach my $species(sort {$hashCount{'Species'}{$b}{'Count'} <=> $hashCount{'Species'}{$a}{'Count'}} keys %{$hashCount{'Species'}})
    {
        my $count       = $hashCount{'Species'}{$species}{'Count'};
        my $perc        = sprintf "%0.4f", $count / $reads_count;
        my $gb_info     = join ",", sort {$hashCount{'Species'}{$species}{'GB'}{$b}     <=> $hashCount{'Species'}{$species}{'GB'}{$a}}     keys %{$hashCount{'Species'}{$species}{'GB'}}; # 物种的GB编号
        my @datas       = ($species, $count, $perc, $gb_info);

        foreach my $level(@need_annos)
        {   
            if($species eq "NO_BLAST")
            {
                push @datas, "";
            }else
            {
                my $anno_info = $hashCount{'Species'}{$species}{$level};
                my $count     = $hashCount{'Level'}{$level}{$anno_info};
                push @datas, "$anno_info($count)";
            }

        }
        print SUMMARY (join "\t", @datas) . "\n";
    }
    close SUMMARY;

    #####################
    # (3) 比对失败序列汇总
    #####################
    open UNMAPPED, ">$summary.unmapped.txt";
    print UNMAPPED "ReadsName\tSeq\n";

    my $IN = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
    while(my $inSeq = $IN->next_seq)
    {
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;
        next if(exists $hashBlastReads{$id});
        print UNMAPPED "$id\t$seq\n";
    }
    close UNMAPPED;

    #####################
    # (4) SUMALL汇总模式
    #####################
    open SUMALL, ">$summary.sum.all.txt";
    print SUMALL "Species\tReadsCount\tPerc\t" . (join "\t", @need_annos) . "\n";
    foreach my $species( sort {$hashSUMALL{$b}{'UniqueCount'} <=> $hashSUMALL{$a}{'UniqueCount'}} keys %hashSUMALL)
    {   
        my $count  = $hashSUMALL{$species}{'UniqueCount'};
        my $perc   = sprintf "%0.4f", $count / $reads_count;

        my @datas  = ($species, $count, $perc);
        foreach my $level(@need_annos)
        {   
            my $anno_info = exists $hashAnno{$species}{$level} ? $hashAnno{$species}{$level} : "";
            push @datas, $anno_info;
        }

        print SUMALL (join "\t", @datas) . "\n";
    }
    close SUMALL;

    #####################
    # (5) SUM Decrease汇总模式
    #####################
    open SUMDECREASE, ">$summary.sum.decrease.txt";
    print SUMDECREASE "Species\tReadsCount\tPerc\t" . (join "\t", @need_annos) . "\n";
    my $max_circle = 20;  # 最多显示前n个物种
    my $now_circle = 0;
    while(1)
    {   
        $now_circle++;


        # 记录当前的reads下，物种包含的reads名称与数量
        my @reads_all = keys %hashBlastReads;
        last if(scalar(@reads_all) == 0);  #  所有reads已分配完毕，跳出循环

        my %hashSpecies;
        foreach my $reads(@reads_all)
        {
            foreach my $species(keys %{$hashBlastReads{$reads}})
            {   
                $hashSpecies{$species}{$reads} = 1;  
                $hashSpecies{$species}{'UniqueCount'}++;
            }
        }
        
        # 获取数量最多的物种
        my ($max_species, @other_species) = sort {$hashSpecies{$b}{'UniqueCount'} <=> $hashSpecies{$a}{'UniqueCount'}} keys %hashSpecies;
        my $count  = $hashSpecies{$max_species}{'UniqueCount'};
        my $perc   = sprintf "%0.4f", $count / $reads_count;
        
        my @datas  = ($max_species, $count, $perc);
        foreach my $level(@need_annos)
        {   
            my $anno_info = exists $hashAnno{$max_species}{$level} ? $hashAnno{$max_species}{$level} : "";
            push @datas, $anno_info;
        }
        print SUMDECREASE (join "\t", @datas) . "\n";
 
        # 从数据库里删除已统计过的reads
        foreach my $reads(keys %{$hashSpecies{$max_species}})
        {   
            next if($reads eq "UniqueCount");
            delete $hashBlastReads{$reads};
        }
        
        # 
        last if($now_circle == $max_circle);
    }
    close SUMDECREASE;

}


# 读取lineage数据库: tax to species
sub read_lineage{
    my $lineage_csv = shift @_;
    my @need_cols   = @_; # 要保留的列

    # 数据split最多列
    my $max_col = (sort {$b<=>$a} @need_cols)[0];
    my $split_col = $max_col + 2;

    print "Read $lineage_csv\n";
    my %hashTaxInfo;
    open LINEAGE, $lineage_csv;
    my $line1 = <LINEAGE>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /,/, $line1, $split_col;
    while(<LINEAGE>)
    {
        $_=~s/[\r\n]//g;
        my @datas = split /,/, $_, $split_col;

        my $tax_id = $datas[0];
        foreach my $col(@need_cols)
        {
            $hashTaxInfo{$tax_id}{$heads[$col]} = $datas[$col];
        }
    }
    close LINEAGE;
    return %hashTaxInfo;
}

# fastq转fasta
sub fastq_to_fasta{
    my $fastq = shift @_;
    my $fasta = shift @_;
    my $limit = shift @_; # reads数量限制

    my $count = 0;
    # fastq句柄
    if($fastq =~ /\.gz$/) # 压缩
    {
        open FASTQ, "gzip -cd $fastq|";
    }
    else# 非压缩
    {
        open FASTQ, $fastq;
    }

    # fasta句柄
    if($fasta =~ /\.gz$/)
    {
        open FASTA, "|gzip > $fasta";
    }
    else
    {
        open FASTA, ">$fasta";
    }

    # 开始处理
    while( my $line1 = <FASTQ>)
    {
        my $line2 = <FASTQ>;
        my $line3 = <FASTQ>;
        my $line4 = <FASTQ>;

        $line1 =~ s/^@//;
        print FASTA ">$line1$line2";
        $count++;
        last if($count >= $limit);
    }
    close FASTQ;
    close FASTA;    
}

# 创建目录
sub make_dir{
    my $dir = shift @_;
    mkdir $dir if(not -e $dir);
}

# 进度条， 根据输入向量计算
sub process_bar_array{
    my $process_name   = shift @_; # 需要分析的对象
    my $process_arrays = shift @_; # 总列表
    my $process_all_count = @$process_arrays;# 需要分析的总列表数量
    my ($pos) = grep $$process_arrays[$_] eq $process_name, 0..$process_all_count-1;
    my $process_count = $pos+ 1; # 当前对象的位置
    my $process_perc = sprintf "%0.2f", 100 * $process_count/$process_all_count; # 进度
    my $windows = 100; # 窗口宽度
    my $finished = $windows * int($process_perc) / 100; # 完成
    my $unfinished = $windows - $finished; # 未完成
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% ";
    print "\n" if($process_count == $process_all_count);   
}

sub help{
    my $info = "
Program: blast+ blastn， 挑选fastq序列比对
Version: 2019-06-18
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
    必填：
        --fastq_dir/-i    fastq文件路径，例如：/home/pub/project/genetics/18B1125B
        --output_dir/-o   结果输出目录，例如：/home/wangly/work/exome/18B1125B/report/nt_test 。流程会自动创建结果目录。

    选填：
        --sample_list/-s     分析的样本，默认为fastq_dir下的所有样本，也可以指定分析样本，例如：A1,B1,C1
        --reads_count/-c     分析的序列数量，默认为：10000
        --thread/-t          并行样本数量，默认：3
        --word_size/-w       blastn word_size 参数，>= 4
        --evalue/-e          evalue 参数，默认：1e-5
        --perc_identity/-p   一致性比例参数，默认：97
        --max_target_seqs/-m 保留的数据库序列比对数量，如果一个数据库序列有多个同源序列，可能会生成n个这条序列的同源比对结果。值越大，生成的文件越大，运行时间也越长。默认：200
        --max_hsps           设置同一个query和同一个subject的匹配结果最多可以在结果文件中保留多少个。理论上，保留的是最佳比对结果。值越大，生成的文件越大，运行时间也越长。 默认：1
        --seqidlist          nt库比对到文件中指定的物种id中，如果不指定，则比对到nt库所有物种id。文件中只有一列数据，一行一个id编号，例如：U13103.1。默认：不过滤
                             参考示例（nt库中线粒体、叶绿体物种ID）：/home/genesky/database/ncbi/nt/mitochondrion_chloroplast_id.txt
        --help/-h            查看帮助文档

    # 输出文件简介
    NT_Blast_Summary.xlsx  仅统计reads的第一个比对结果，汇总物种注释
    NT_Blast_Summary_SumAll.xlsx  统计reads的所有比对结果，汇总注释。即：一条reads同时统计到多个物种上
    NT_Blast_Summary_SumDecrease.xlsx  累计减法统计reads的所有比对结果，汇总注释。过程：记录每一条reads的所有比对物种，统计丰度最高的物种并输出，然后删除所有比对到这个物种的reads记录（等同于这条reads没有参与比对过程），然后再次统计丰度最高的物种并输出，依次循环。
    NT_Blast_unmapped.xlsx  比对失败的reads序列

    示例：
    （1）常规用法
    perl nt_blast.pl -i /home/project/abc -o /home/test/nt_blast
    
    (2) 线粒体、叶绿体细胞器组装
    perl nt_blast.pl -i /home/project/abc -o /home/test/nt_blast -p 90 --seqidlist /home/genesky/database/ncbi/nt/mitochondrion_chloroplast_id.txt -c 100000
    \n";
    return $info;
}

# my $excel_summary        = "$output_dir/NT_Blast_Summary.xlsx";        # 每一条序列只取最优的一条比对结果
# my $excel_sum_all        = "$output_dir/NT_Blast_Summary_SumAll.xlsx"; # 取所有比对结果，看每一个物种的比对情况。如果想充分利用这个结果，最好解除max_target_seqs参数限制。
# my $excel_sum_decrease   = "$output_dir/NT_Blast_Summary_SumDecrease.xlsx"; # 在SumAll的基础上，依次挑选丰度最高物种，输出统计，并删除该物种包含的reads记录，剩下的reads再重新挑选丰度最高的物种...
# my $excel_unmapped       = "$output_dir/NT_Blast_unmapped.xlsx";       # 比对失败的结果
