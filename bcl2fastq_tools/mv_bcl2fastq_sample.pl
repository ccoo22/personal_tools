# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Parallel::ForkManager;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量


# 检测 -> 脚本输入
my ($input_dir, $output_dir, $threshold, $if_help);
GetOptions(
    "input_dir|i=s"   => \$input_dir,
    "output_dir|o=s"  => \$output_dir,
    "threshold|t=s"   => \$threshold,
    "help|h"          => \$if_help,
);
die help() if(defined $if_help or (not defined $input_dir or not defined $output_dir));
$threshold = 1 if(not defined $threshold);
###################################################################### 主程序
 
# (1) 创建输出目录
print "1: create output directory\n";
mkdir $output_dir if(not -e $output_dir);
die "[Error] output directory cannot write\n" if(not -w $output_dir); # 输出目录不可写
die "[Error] input directory cannot write\n" if(not -w $input_dir); # 输入目录不可写

#(2) 获取输入项目目录下的样本编号ID/样本名，并检测文件读写权限
print "2: check sample file\n";
my %hashID;
my $is_exists = 0; # 样本是否已存在
opendir INPUT, $input_dir;
foreach my $sample_id(readdir INPUT)
{
    next if($sample_id !~ /\w/);
    die "[Warnings] file is not put into sample_id directory:  $sample_id\n" if(-f "$input_dir/$sample_id"); # 当样本ID=样本名时，拆分的数据直接裸露在项目目录下，此时需要手工按照数据放置规则进行存放
    
    # 获取样本名
    my $sample_name = "";
    my @files;
    opendir FILE, "$input_dir/$sample_id";
    foreach my $file(readdir FILE)
    {
        next if($file !~ /\w/);
        push @files, $file;
        ($sample_name) = $file =~ /(\S+)_S\d+_R\d/;
    }
    close FILE;    
    die "[Error] file count is not 2 ： $sample_id\n" if(scalar(@files) != 2);

    # 输入
    my $input_R1 = "$input_dir/$sample_id/$sample_name\_S$sample_id\_R1_001.fastq.gz";
    my $input_R2 = "$input_dir/$sample_id/$sample_name\_S$sample_id\_R2_001.fastq.gz";
    if(not -e $input_R1)
    {
        print "[Warnings] input file lost     : $input_R1\n";
        next;
    }
 
    # 输出
    my $output_R1 = "$output_dir/$sample_name\_R1.fastq.gz";
    my $output_R2 = "$output_dir/$sample_name\_R2.fastq.gz";
    $is_exists    = 1 if(-e $output_R1); # 输出文件已存在
    # bak
    my $bak_dir  = "$output_dir/bak";
    my $bak_R1   = "$bak_dir/$sample_name\_R1.fastq.gz";
    my $bak_R2   = "$bak_dir/$sample_name\_R2.fastq.gz";

    # 数据读写权限检测
    die "[Error] output file cannot write     : $output_R1\n" if(-e $output_R1 and not -w $output_R1); # 输出文件不可写
    die "[Error] output bak file cannot write : $bak_R1\n"    if(-e $bak_R1 and not -w $bak_R1); # bak输出文件不可写
    
    # 记录
    $hashID{$sample_id}{'sample_name'} = $sample_name;
    $hashID{$sample_id}{'R1'}          = $input_R1;
    $hashID{$sample_id}{'R2'}          = $input_R2;
    
}
close INPUT;

# 备份路径
my $bak_dir = "$output_dir/bak";
mkdir $bak_dir if($is_exists == 1 and not -e $bak_dir); # 需要合并
die "[Error] output bak directory cannot write\n" if($is_exists == 1 and not -w $bak_dir); # 输出目录不可写


#（3）输出
print "3: start output\n";
my @samples_id = sort {$a<=>$b} keys %hashID;
my $pm = Parallel::ForkManager->new($threshold);
   $pm->run_on_start(sub{my ($pid, $sample_id) = @_; process_bar_array($sample_id, \@samples_id)});# 进度条
foreach my $sample_id(@samples_id)
{   
    # 注：单线程不要用pm,会比较慢,多线程时才使用pm
    if($threshold == 1)
    {   
        process_bar_array($sample_id, \@samples_id); # 进度条
        process_sample($sample_id, \%hashID); # 处理
    }else
    {
        $pm->start($sample_id) and next;
        process_sample($sample_id, \%hashID);
        $pm->finish;          
    } 
}
$pm->wait_all_children;
 




###################################################################### 子函数

sub process_sample{
    my $sample_id   = shift @_;
    my $hashID      = shift @_;
    my $sample_name = $hashID->{$sample_id}{'sample_name'};
    
    # 输入
    my $input_R1 = $hashID->{$sample_id}{'R1'};
    my $input_R2 = $hashID->{$sample_id}{'R2'};
    # 输出
    my $output_R1 = "$output_dir/$sample_name\_R1.fastq.gz";
    my $output_R2 = "$output_dir/$sample_name\_R2.fastq.gz";
    # bak
    my $bak_dir  = "$output_dir/bak";
    my $bak_R1   = "$bak_dir/$sample_name\_R1.fastq.gz";
    my $bak_R2   = "$bak_dir/$sample_name\_R2.fastq.gz";

    # 结果不存在，直接移动
    if(not -e $output_R1)
    {
        print("mv $input_R1 $output_R1\n");
        print("mv $input_R2 $output_R2\n");

        system("mv $input_R1 $output_R1");
        system("mv $input_R2 $output_R2");
    }else{ # 结果已存在
        # 首先，上一次的结果放入备份中
        mkdir $bak_dir if(not -e $bak_dir); # 下机数据存在重名样本，流程第2步是检测不出来的，故需检测并创建bak目录
        system("mv $output_R1 $bak_dir/");
        system("mv $output_R2 $bak_dir/");

        # 然后，数据合并
        print("cat $input_R1 $bak_R1 > $output_R1\n");
        print("cat $input_R2 $bak_R2 > $output_R2\n");

        system("cat $input_R1 $bak_R1 > $output_R1");
        system("cat $input_R2 $bak_R2 > $output_R2");
    }

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
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% \n";
}

sub help{
    my $info = "
Program: move bcl2fastq to dir， 移动下机拆分数据到指定目录
Version: 2019-05-06
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_dir/-i    输入目录，例如：/home/lch/work/20190520_Hiseq_SH19GSL091/Raw_data/19B0506B
         --output_dir/-o   输出目录，例如：/home/pub/project/Transcriptome/19B0506B
         --threshold/-t    并行样本数量
         --help/-h         查看帮助文档
    \n";
    return $info;
}

