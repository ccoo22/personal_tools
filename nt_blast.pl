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
my $blastn      = "/home/ganb/soft/ncbi-blast-2.7.1+/bin/blastn";
my $lineages    = "/home/genesky/database/ncbi/nt_20191128/lineages-2019-02-20.csv";
my $nt_db       = "/home/genesky/database/ncbi/nt/blast_idx/nt";
my $table2excel = "perl /home/pub/bin/NGS/chip/GATK4/tools/personal/table2excel.pl";

# 检测 -> 脚本输入
my ($fastq_dir, $output_dir, $sample_list, $reads_count, $thread, $word_size, $evalue, $perc_identity, $max_target_seqs, $if_help);
GetOptions(
    "fastq_dir|i=s"     => \$fastq_dir,
    "output_dir|o=s"    => \$output_dir,
    "sample_list|s=s"   => \$sample_list,
    "reads_count|c=s"   => \$reads_count,
    "thread|t:i"        => \$thread,
    "word_size|w:i"     => \$word_size,
    "evalue|e:i"        => \$evalue,
    "perc_identity|p:i" => \$perc_identity,
    "max_target_seqs|p:i" => \$max_target_seqs,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $fastq_dir or not defined $output_dir));
###################################################################### 主程序
$sample_list   = ""    if(not defined $sample_list);
$reads_count   = 10000 if(not defined $reads_count);
$thread        = 3     if(not defined $thread);
$word_size     = (not defined $word_size)     ? ""                  : "-word_size $word_size";
$evalue        = (not defined $evalue)        ? "-evalue 1e-5"      : "-evalue $evalue";
$perc_identity = (not defined $perc_identity) ? "-perc_identity 97" : "-perc_identity $perc_identity";
$max_target_seqs = (not defined $max_target_seqs) ? "-max_target_seqs 5" : "-max_target_seqs $max_target_seqs";
 
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
    while(readdir FASTQ_DIR)
    {
        my ($sample) = $_ =~ /(.*)_R1.fastq.gz/;
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
my %hashTax2Species = read_lineage($lineages, 7); # 读取物种注释
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
    system("$blastn -query $fasta -db $nt_db -out $blast -outfmt '7 qacc sseqid staxids pident qcovs length qstart qend sstart send evalue salltitles' $evalue $perc_identity $word_size $max_target_seqs -num_threads 20    -dust no"); # 比对  
    summary($blast, $fasta, $summary, \%hashTax2Species, $reads_count); # 注释species并汇总

    $pm->finish;    
}
$pm->wait_all_children;

#########
# (5) summary to excel
#########
print "5: summary to excel\n";

my $excel_summary        = "$output_dir/NT_Blast_Summary.xlsx";
my $excel_unmapped       = "$output_dir/NT_Blast_unmapped.xlsx";
my $summary_list  = "";
my $unmapped_list = "";
my $sample_list2  = "";
foreach my $sample(@samples)
{
    my $summary    = "$tmp_dir/$sample.summary.txt";
    my $unmapped   = "$summary.unmapped.txt";
    $summary_list  .= "$summary,";
    $unmapped_list .= "$unmapped,";
    $sample_list2  .= "$sample,";

}
system("$table2excel -i $summary_list  -s $sample_list2 -o $excel_summary");
system("$table2excel -i $unmapped_list -s $sample_list2 -o $excel_unmapped");
print "final: $excel_summary\n";



###################################################################### 子函数

# 注释比对结果
sub summary{
    my $blast           = shift @_; # 比对结果
    my $fasta           = shift @_; # 原始fasta序列
    my $summary         = shift @_; # 最终汇总结果
    my $hashTax2Species = shift @_; # tax -> species 映射
    my $reads_count     = shift @_; # 项目分析总reads数量
    
    # (1) 注释汇总
    my %hashCount;
    my %hashBlastReads;
    my $count_blast; # 比对成功的reads
    my $last_reads_name = "";
    open BLAST, $blast;
    while(<BLAST>)
    {
        next if($_=~/^#/);
        $_=~s/[\r\n]//g;
        my ($reads_name, $gi, $tax_id, $tmp) = split /\t/, $_, 4;
        my $species = (exists $hashTax2Species{$tax_id}{'species'} and $hashTax2Species{$tax_id}{'species'} =~ /\w/) ? $hashTax2Species{$tax_id}{'species'} : "NO_RANK";
        
        # 只看比对的第一个结果
        next if($reads_name eq $last_reads_name); 
        $last_reads_name = $reads_name;
        
        $hashBlastReads{$reads_name}++;
        $hashCount{$species}++;
        $count_blast++;
    }
    close BLAST;
    $hashCount{'NO_BLAST'} = $reads_count - $count_blast if($count_blast < $reads_count); # 部分序列没有比对上

    # (2) 汇总结果输出
    open SUMMARY, ">$summary";
    print SUMMARY "Species\tReadsCount\tPerc\n";
    foreach my $species(sort {$hashCount{$b} <=> $hashCount{$a}} keys %hashCount)
    {
        my $count = $hashCount{$species};
        my $perc  = sprintf "%0.4f", $count / $reads_count;
        print SUMMARY "$species\t$count\t$perc\n";
    }
    close SUMMARY;

    # (3) 比对失败序列汇总
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
        --sample_list/-s   分析的样本，默认为fastq_dir下的所有样本，也可以指定分析样本，例如：A1,B1,C1
        --reads_count/-c   分析的序列数量，默认为：10000
        --thread/-t        并行样本数量，默认：3
        --word_size/-w     blastn word_size 参数，>= 4
        --evalue/-e        evalue 参数，默认：1e-5
        --perc_identity/-p 一致性比例参数，默认：97
        --max_target_seqs/-m 保留的比对数量，默认：5
        --help/-h          查看帮助文档
    \n";
    return $info;
}

