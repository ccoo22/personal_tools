use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;
use Cwd 'abs_path';
$|=1;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];

my $DEFAULT_THREAD = 20;
my $blast_plus_dir = "/home/genesky/software/blast+/2.7.1/bin";
my $fastx_dir      = "/home/genesky/software/fastx_toolkit/0.0.14/bin";


# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($map_file, $fastq_dir, $output_dir, $thread, $primer_f, $primer_r, $if_help);
GetOptions(
    "map_file=s" => \$map_file,
    "fastq_dir=s"  => \$fastq_dir,
    "output_dir=s" => \$output_dir,
    "primer_f=s"   => \$primer_f,
    "primer_r=s"   => \$primer_r,
    "thread=i"   => \$thread,
    "help|h" => \$if_help,
);
die "
根据引物，提取当前样本的fastq数据，并分离出匹配失败的fastq，放入其他文件
[注意] 该脚本主要适用于微生物，只有一个引物的项目

Options: 必填

        --map_file            要分析的样本，以及其他样本映射关系，三列，有表头，分别是样本名、其他项目编号、其他样本名
                              样本名：需要用引物提取的样本。 最终输出到 output_dir/extract/sample_R1.fastq.gz  sample_R2.fastq.gz
                              其他项目编号：匹配失败的数据对应的项目编号。如果没有，可任意填写
                              其他样本名：匹配失败的数据的样本名。 匹配失败的数据存放到 output_dir/其他项目编号/其他样本名_R1.fastq.gz 其他样本名_R2.fastq.gz
        --fastq_dir           样本fastq路径
        --output_dir          结果输出目录
        --primer_f            F端引物
        --primer_r            R端引物

Options: 可选

        --thread                   使用线程数 (default: $DEFAULT_THREAD)
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $map_file or not defined $fastq_dir or not defined $output_dir or not defined $primer_f or not defined $primer_r);
$primer_f = uc($primer_f);
$primer_r = uc($primer_r);
$thread = $DEFAULT_THREAD if(not defined $thread);
$map_file   = abs_path($map_file);
$fastq_dir  = abs_path($fastq_dir);
$output_dir = abs_path($output_dir);
my %hashJianbing = get_jianbing();  # 获取简并碱基数据库

###################################################################### 主程序
my $tmp_dir = "$output_dir/tmp";
make_dir($output_dir, $tmp_dir);

my %hashMap = read_map_file($map_file);
my @samples = sort keys %hashMap;

# other项目目录创建
map{ make_dir("$output_dir/$hashMap{$_}{'other_project'}")  } keys %hashMap;
make_dir("$output_dir/extract");
##############
# （1）引物建立索引
##############
print "\n[process] 引物建立索引\n";
my $primer = "$tmp_dir/primer.fa";
build_primer_index($primer);
 
##############
# (2) 比对
##############
print "\n[process] 比对\n";
blast(\@samples, $primer, $tmp_dir);

##############
# (3) 比对结果解析
##############
print "\n[process] 比对结果解析 \n";
parse_blast(\@samples, $tmp_dir); 

##############
# (4) 数据拆分
##############
print "\n[process] 数据拆分\n";
parse_sample(\@samples, $tmp_dir, $output_dir); 



###################################################################### 子程序

sub parse_sample{
    my $samples = shift @_;
    my $tmp_dir = shift @_;
    my $output_dir = shift @_;

    my @sample_runs = @$samples;
    my $pm = Parallel::ForkManager->new($thread);
       $pm->run_on_start(sub{my ($pid, $sample) = @_; process_bar_array($sample, \@sample_runs)});# 进度条
    foreach my $sample(@sample_runs)
    {
        $pm->start($sample) and next;
        my $other_project = $hashMap{$sample}{'other_project'};
        my $other_sample  = $hashMap{$sample}{'other_sample'};
        
        # 原始fastq
        my $r1_raw = "$fastq_dir/$sample\_R1.fastq.gz";
        my $r2_raw = "$fastq_dir/$sample\_R2.fastq.gz";

        # 根据引物拆分出来的fastq
        my $r1_extract = "$output_dir/extract/$sample\_R1.fastq.gz";
        my $r2_extract = "$output_dir/extract/$sample\_R2.fastq.gz";

        # 根据引物拆分失败的fastq
        my $r1_other = "$output_dir/$other_project/$other_sample\_R1.fastq.gz";
        my $r2_other = "$output_dir/$other_project/$other_sample\_R2.fastq.gz";

        # 读入blast结果
        my %hashExtract;
        open RESULT, "$tmp_dir/$sample.result";
        while(<RESULT>)
        {   
            $_=~s/[\r\n]//g;
            my ($read_name, $status) = split /\t/, $_;
            $hashExtract{$read_name} = 1 if($status eq 'OK');
        }
        close RESULT;

        # 拆
        open R1_RAW, "gzip -cd $r1_raw|";
        open R2_RAW, "gzip -cd $r2_raw|";
        open R1_EXTRACT, "|gzip > $r1_extract";
        open R2_EXTRACT, "|gzip > $r2_extract";
        open R1_OTHER, "|gzip > $r1_other";
        open R2_OTHER, "|gzip > $r2_other";
        while(my $r1_line1 = <R1_RAW>)
        {
            my $r1_line2 = <R1_RAW>;
            my $r1_line3 = <R1_RAW>;
            my $r1_line4 = <R1_RAW>;

            my $r2_line1 = <R2_RAW>;
            my $r2_line2 = <R2_RAW>;
            my $r2_line3 = <R2_RAW>;
            my $r2_line4 = <R2_RAW>;

            my $read_name = $r1_line1;
               ($read_name) = split /\s+/, $read_name;
            $read_name=~s/^\@//;

            if(exists $hashExtract{$read_name})
            {
                print R1_EXTRACT "$r1_line1$r1_line2$r1_line3$r1_line4";
                print R2_EXTRACT "$r2_line1$r2_line2$r2_line3$r2_line4";
            }else{
                print R1_OTHER "$r1_line1$r1_line2$r1_line3$r1_line4";
                print R2_OTHER "$r2_line1$r2_line2$r2_line3$r2_line4"; 
            }
        }
        close R1_RAW;
        close R2_RAW;
        close R1_EXTRACT;
        close R2_EXTRACT;
        close R1_OTHER;
        close R2_OTHER;
        $pm->finish;    
    }
    $pm->wait_all_children;
}



sub parse_blast{
    my $samples = shift @_;
    my $tmp_dir = shift @_;

    my $summary_file = "$tmp_dir/summary.txt";
    open SUMMARY, ">$summary_file";
    print SUMMARY "sample\ttotal reads\tblast reads\tblast perc\n";
    my @sample_runs = @$samples;
    my $pm = Parallel::ForkManager->new($thread);
       $pm->run_on_start(sub{my ($pid, $sample) = @_; process_bar_array($sample, \@sample_runs)});# 进度条
    foreach my $sample(@sample_runs)
    {
        $pm->start($sample) and next;
        # 读入blast结果
        my $r1_fasta = "$tmp_dir/$sample.r1.fasta";
        my $r1_blast = "$tmp_dir/$sample.r1.blast";
        my $r2_blast = "$tmp_dir/$sample.r2.blast";
        my %hashBlast;
        read_blast($r1_blast, \%hashBlast);
        read_blast($r2_blast, \%hashBlast);

        # 总reads数量
        my $total = `wc -l $r1_fasta`;
           ($total) = split /\s+/, $total;
        $total = $total / 2;

        # 结果输出
        my $result = "$tmp_dir/$sample.result";
        open RESULT, ">$result";
        my $ok = 0;
        foreach my $reads_name(keys %hashBlast)
        {   
            if(exists $hashBlast{$reads_name}{'primerF'} and exists $hashBlast{$reads_name}{'primerR'})
            {
                print RESULT "$reads_name\tOK\n";
                $ok++;
            }else
            {
                print RESULT "$reads_name\tERROR\n";
            }
        }
        close RESULT;

        my $ok_perc = sprintf "%0.4f", $ok / $total;
        print SUMMARY "$sample\t$total\t$ok\t$ok_perc\n";

        $pm->finish;    
    }
    $pm->wait_all_children;
    close SUMMARY;
}
sub read_blast{
    my $blast = shift @_;
    my $hashBlast = shift @_;

    my %hashPrimerLen;
       $hashPrimerLen{'primerF'} = length($primer_f);
       $hashPrimerLen{'primerR'} = length($primer_r);

    open BLAST,$blast;
    while(<BLAST>)
    {
        $_=~ s/[\r\n]//g;
        my ($reads_name, $target_name, $identify_perc, $identify_len, $mismatch, $gap, $reads_start, $reads_end, $target_start, $target_end, $evalue, $score) = split /\t/, $_;
        ($target_name) = split /-/, $target_name;  # 去掉-后的所有字符，因为有简并碱基处理

        next if($reads_start > 5);  # 引物比对起始位置必须在开头5bp以内
        next if(($hashPrimerLen{$target_name} - $identify_len) > 4 ); # 引物只允许最多4个碱基匹配丢失
        $hashBlast->{$reads_name}{$target_name}++;
 
    }
    close BLAST;
}



sub blast{
    my $samples = shift @_;
    my $primer  = shift @_;
    my $tmp_dir = shift @_;

    my @sample_runs = @$samples;
    my $pm = Parallel::ForkManager->new($thread);
       $pm->run_on_start(sub{my ($pid, $sample) = @_; process_bar_array($sample, \@sample_runs)});# 进度条
    foreach my $sample(@sample_runs)
    {
        $pm->start($sample) and next;
        my $r1 = "$fastq_dir/$sample\_R1.fastq.gz";
        my $r2 = "$fastq_dir/$sample\_R2.fastq.gz";
        system "gzip -cd $r1 | $fastx_dir/fastq_to_fasta -Q 33 -n  -o $tmp_dir/$sample.r1.fasta";
        system "gzip -cd $r2 | $fastx_dir/fastq_to_fasta -Q 33 -n  -o $tmp_dir/$sample.r2.fasta";
        system "$blast_plus_dir/blastn -query $tmp_dir/$sample.r1.fasta -outfmt 6 -evalue 0.001 -db $primer -out  $tmp_dir/$sample.r1.blast -num_threads 4 -max_target_seqs 5 -word_size 4";
        system "$blast_plus_dir/blastn -query $tmp_dir/$sample.r2.fasta -outfmt 6 -evalue 0.001 -db $primer -out  $tmp_dir/$sample.r2.blast -num_threads 4 -max_target_seqs 5 -word_size 4";

        $pm->finish;    
    }
    $pm->wait_all_children;
}


sub read_map_file{
    my $map_file = shift @_;
    my %hashMap;
    
    open MAP, $map_file;
    <MAP>;
    while(<MAP>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        my ($sample, $other_project, $other_sample) = split /\t/, $_;
        my $r1 = "$fastq_dir/$sample\_R1.fastq.gz";
        my $r2 = "$fastq_dir/$sample\_R2.fastq.gz";

        die "[Error] map_file 必须包含三列，且不能为空\n" if(not defined $other_sample or not defined $other_project or $other_sample !~ /\w/ or $other_project !~ /\w/);
        die "[Error] fastq_dir 丢失 $sample 的R1/R2 fastq文件\n" if(is_file_ok($r1, $r2) == 0);

        $hashMap{$sample}{'other_project'} = $other_project;
        $hashMap{$sample}{'other_sample'} = $other_sample;
    }
    close MAP;
    return %hashMap;
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

sub build_primer_index{
    my $primer = shift @_;

    open PRIMER, ">$primer";
    my $count = 0;
    foreach my $seq(replace_jianbing_bases($primer_f))# 替换简并碱基/获取所有碱基组合后的新序列
    {
        $count++;
        print PRIMER ">primerF-$count\n$seq\n";
    }
    foreach my $seq(replace_jianbing_bases($primer_r))# 替换简并碱基/获取所有碱基组合后的新序列
    {
        $count++;
        print PRIMER ">primerR-$count\n$seq\n";
    }
    close PRIMER;
    system("$blast_plus_dir/makeblastdb -in $primer -dbtype nucl");
}

sub replace_jianbing_bases{
    my $seq = shift @_;
    
    my @seqs_replace;
    if($seq =~/$hashJianbing{'base_all'}/)
    {   
        my @seqs_replace_tmp = replace_first_jianbing($seq, $&);  # 替换碱基
        foreach my $seq_replace_tmp(@seqs_replace_tmp)
        {
            my @seqs_replace_final = replace_jianbing_bases($seq_replace_tmp);
            push @seqs_replace, @seqs_replace_final;
        }
    }else
    {
        push @seqs_replace, $seq;  # 迭代终止
    }
    return @seqs_replace;
}

# 替换第一个简并碱基
sub replace_first_jianbing{
    my $seq      = shift @_;
    my $jianbing = shift @_;
    
    my @seqs_replace_first;
    foreach my $base(keys %{$hashJianbing{'base'}{$jianbing}})
    {
        my $seq_replace_first = $seq;
        $seq_replace_first=~s/$jianbing/$base/;
        push @seqs_replace_first, $seq_replace_first;
    }
    return @seqs_replace_first;
}
# 简并碱基数据库
sub get_jianbing{
    my %hashJianbing;
    $hashJianbing{'base'}{'R'}{'A'} = 1;
    $hashJianbing{'base'}{'R'}{'G'} = 1;
    
    $hashJianbing{'base'}{'Y'}{'C'} = 1;
    $hashJianbing{'base'}{'Y'}{'T'} = 1;

    $hashJianbing{'base'}{'M'}{'A'} = 1;
    $hashJianbing{'base'}{'M'}{'C'} = 1;

    $hashJianbing{'base'}{'K'}{'G'} = 1;
    $hashJianbing{'base'}{'K'}{'T'} = 1;

    $hashJianbing{'base'}{'S'}{'G'} = 1;
    $hashJianbing{'base'}{'S'}{'C'} = 1;

    $hashJianbing{'base'}{'W'}{'A'} = 1;
    $hashJianbing{'base'}{'W'}{'T'} = 1;

    $hashJianbing{'base'}{'H'}{'A'} = 1;
    $hashJianbing{'base'}{'H'}{'T'} = 1;
    $hashJianbing{'base'}{'H'}{'C'} = 1;

    $hashJianbing{'base'}{'B'}{'G'} = 1;
    $hashJianbing{'base'}{'B'}{'T'} = 1;
    $hashJianbing{'base'}{'B'}{'C'} = 1;

    $hashJianbing{'base'}{'V'}{'G'} = 1;
    $hashJianbing{'base'}{'V'}{'A'} = 1;
    $hashJianbing{'base'}{'V'}{'C'} = 1;

    $hashJianbing{'base'}{'D'}{'G'} = 1;
    $hashJianbing{'base'}{'D'}{'A'} = 1;
    $hashJianbing{'base'}{'D'}{'T'} = 1;

    $hashJianbing{'base'}{'N'}{'A'} = 1;
    $hashJianbing{'base'}{'N'}{'T'} = 1;
    $hashJianbing{'base'}{'N'}{'C'} = 1;
    $hashJianbing{'base'}{'N'}{'G'} = 1;

    $hashJianbing{'base_all'} = join "|", sort keys %{$hashJianbing{'base'}};
    return %hashJianbing;
}

# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}
