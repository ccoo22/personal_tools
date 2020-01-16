use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
$|=1;

# 定义 -> 核心变量
my $TrimGalore = "/home/genesky/software/trim_galore/0.5.0/trim_galore";
my $CutAdapt   = "/home/genesky/software/python/2.7.15/bin/cutadapt";

# 定义 -> 默认参数
my $length = 50;
my $stringency = 5;
my $adapter    = "default";
my $run_sample = "default";
my $thread = 4;

# 检测 -> 脚本输入
my  ($input_dir,$output_dir,$help);
GetOptions(
     'input_dir|i=s'    => \$input_dir,
     'output_dir|o=s'	 => \$output_dir,
	 'length|l=s'        => \$length,
	 'stringency|s=s'    => \$stringency,
	 'adapter|a=s'       => \$adapter,
	 'runsample|r=s'     => \$run_sample,
	 'thread|m=s'        => \$thread,
     'help|h'            => \$help,
);
die help() if(defined $help or (not defined $input_dir or not defined $output_dir ));

###################################################################### 主程序

# 输入引物
my ($adapter1,$adapter2) = ();
if ($adapter ne "default"){
   ($adapter1,$adapter2) = split/,/, $adapter;
	die "adapter序列只能含ATCG,且adapter序列需成对\n" if ($adapter1 =~ /[^ATCG]/i or $adapter2 =~ /[^ATCG]/i or $adapter1 eq "" or $adapter2 eq "");
}

# 输入指定分析样本

my @samples;
   @samples = split /,/, $run_sample if ($run_sample ne "default");
   
if ($run_sample eq "default"){
    while (my $file = glob "$input_dir/*_R1.fastq.gz")
    { 
        $file =~ s/[\r\n]//g;
    	my ($sample) = $file =~ /$input_dir\/(.*)_R1.fastq.gz/;
        push @samples, $sample;
    }
}

# 检查样本原始序列
my $error = 0;
foreach my $sample(@samples){
    my $file1 = $input_dir."/".$sample."_R1.fastq.gz";
    my $file2 = $input_dir."/".$sample."_R2.fastq.gz";
	my @files = ($file1, $file2);
    if(is_file_ok(@files) == 0){
	    print "Lost $sample input_fastq\n";
		$error ++;
	}
}
die if ($error != 0);

# 运行TrimGalore

my $pm = Parallel::ForkManager->new($thread);
   $pm->run_on_start(sub{my ($pid, $sample) = @_; process_bar_array($sample, \@samples)});# 进度条

foreach my $sample(@samples)
{
    my $pid = $pm->start($sample) and next;
	my $Fastq1 = $input_dir."/".$sample."_R1.fastq.gz";
	my $Fastq2 = $input_dir."/".$sample."_R2.fastq.gz";
	if ($adapter ne 'default'){
        system ("$TrimGalore --phred33 --suppress_warn --no_report_file -q 20 --adapter $adapter1 --adapter2 $adapter2 --stringency $stringency --length $length --path_to_cutadapt $CutAdapt --paired $Fastq1 $Fastq2 -o $output_dir");
	}
    else {
	    system ("$TrimGalore --phred33 --suppress_warn --no_report_file -q 20 --stringency $stringency --length $length --path_to_cutadapt $CutAdapt --paired $Fastq1 $Fastq2 -o $output_dir");
	}
	# TrimGalore过滤结果压缩 重命名
	my $tmp1 = $output_dir."/".$sample."_R1_val_1.fq.gz";
	my $tmp2 = $output_dir."/".$sample."_R2_val_2.fq.gz";
	my $result1 = $output_dir."/".$sample."_R1.fastq.gz";
	my $result2 = $output_dir."/".$sample."_R2.fastq.gz";
	system ("mv $tmp1 $result1");	
	system ("mv $tmp2 $result2");
	
	$pm->finish;
}

$pm->wait_all_children;

###################################################################### 子函数


# 检验文件是否为空

sub is_file_ok{
    my @files = @_;
	my $isok = 1;
	foreach my $file(@files)
	{
	    $isok = 0 if (not -e $file or -s $file == 0);
	}
	return $isok;
} 

# 进度条，据输入向量计算

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
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc%";
    print "\n" if($process_count == $process_all_count);   
}

# 帮助文档

sub help{
    my $info = "
Program: TrimGalore 
Version: 2019-08-16

Usage:   Usage: perl $0  -i input_dir -o output_dir [options]

Options: 必填

         --input_dir/-i      原始fastq文件路径，例如：/home/pub/project/genetics/19B0628A
         --output_dir/-o     输出结果文件夹路径
Options: 可选
         --length/-l         保留的最短序列长度，默认：50
         --runsample/-r      输入指定分析样本名，逗号分隔。默认分析全部样本
         --stringency/-s     设定可接受的前后adapter重叠的碱基数。默认：5
         --adapter/-a        输入adapter序列，逗号分隔。不输入则软件默认自动寻找adapter
         --thread/-t         输入并行样本数。默认为4
         --help/-h           查看帮助文档
    \n";
    return $info;
}




