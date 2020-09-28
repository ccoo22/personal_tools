$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
检查目录下的所有.gz文件的完整性（是否损坏）
Version: v1.0 2020-09-28
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Parallel::ForkManager;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 变量
my $DEFAULT_THREAD          = 10;

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input, $output, $thread, $if_help);
GetOptions(
	"input|i=s"           => \$input,
	"output|o=s"          => \$output,

	"thread=s"            => \$thread,
	"help|h"              => \$if_help,
);
die "
Options: 必填
        --input/-i                输入数据路径
                                  自动检查该目录下的所有.gz文件
        --output/-o               输出文件，记录文件情况

Options: 可选
        --thread                  并行运行数量 (default: '$DEFAULT_THREAD') 
        --help/-h                 查看帮助文档

\n" if (defined $if_help or not defined $input );

$thread          = $DEFAULT_THREAD if (not defined $thread);
  
 
###################################################################### 主程序

opendir FILE, $input;
my @files = grep{ $_=~/\.gz$/ } readdir(FILE);
closedir FILE;

open OUTPUT, ">$output";
print OUTPUT "file\tcondition\n";
my $pm = Parallel::ForkManager->new($thread);
   $pm->run_on_start(sub{my ($pid, $file) = @_; process_bar_array($file, \@files)});# 进度条
foreach my $file(@files)
{
    $pm->start($file) and next;
    my $rtn = `gzip -t $input/$file 2>&1`;
    my $condition = ($rtn eq '') ? 'OK' : 'BAD';
    print OUTPUT "$file\t$condition\n";

    $pm->finish;    
}
$pm->wait_all_children;
close OUTPUT;

###################################################################### 子程序

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

