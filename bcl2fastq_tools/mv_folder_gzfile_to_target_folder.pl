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
print "2: check file\n";
my %hashFile;

my $is_exists = 0; # 样本是否已存在
opendir INPUT, $input_dir;
foreach my $file(readdir INPUT)
{
    next if($file !~ /\.fastq\.gz$/);
    my $input_file = "$input_dir/$file";
    my $output_file = "$output_dir/$file";

    $hashFile{$file}{"input"}  = $input_file;
    $hashFile{$file}{"output"} = $output_file;

    $is_exists = 1 if(-e $output_file);
}
closedir INPUT;

# 备份路径
my $bak_dir = "$output_dir/bak";
mkdir $bak_dir if($is_exists == 1 and not -e $bak_dir); # 需要合并
die "[Error] output bak directory cannot write\n" if($is_exists == 1 and not -w $bak_dir); # 输出目录不可写


#（3）输出
print "3: start output\n";
my @files = sort keys %hashFile;
my $pm = Parallel::ForkManager->new($threshold);
   $pm->run_on_start(sub{my ($pid, $file) = @_; process_bar_array($file, \@files)});# 进度条
foreach my $file(@files)
{   
    # 注：单线程不要用pm,会比较慢,多线程时才使用pm
    if($threshold == 1)
    {   
        process_bar_array($file, \@files); # 进度条
        process_file($file, \%hashFile); # 处理
    }else
    {
        $pm->start($file) and next;
        process_file($file, \%hashFile);
        $pm->finish;          
    } 
}
$pm->wait_all_children;

###################################################################### 子函数
sub process_file{
    my $file     = shift @_;
    my $hashFile = shift @_;
    
    # 输入
    my $input_file = $hashFile->{$file}{'input'};
    
    # 输出
    my $output_file = $hashFile->{$file}{'output'};
    
    # bak
    my $bak_dir  = "$output_dir/bak";
    my $bak_file = "$bak_dir/$file";

    # 结果不存在，直接移动
    if(not -e $output_file)
    {
        print("mv $input_file $output_file\n");

        system("mv $input_file $output_file");
    }else{ # 结果已存在
        # 首先，上一次的结果放入备份中
        system("mv $output_file $bak_dir/");

        # 然后，数据合并
        print("cat $input_file $bak_file > $output_file\n");

        system("cat $input_file $bak_file > $output_file");
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
Program: move folder .fastq.gz file to dir， 移动目录下的.fastq.gz文件到指定目录
Version: 2019-05-30
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_dir/-i    输入目录,把输入目录下的所有.gz结尾的文件移动到目标目录下，如果存在重名，则备份，并cat到一起
         --output_dir/-o   输出目录
         --threshold/-t    并行样本数量
         --help/-h         查看帮助文档
    \n";
    return $info;
}

