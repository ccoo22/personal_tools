# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
 

# 检测 -> 脚本输入
my ($input, $output, $count, $if_help);
GetOptions(
    "input|i=s"   => \$input,
    "output|o=s"  => \$output,
    "count|c=i"   => \$count,
    "help|h"      => \$if_help,
);
die help() if(defined $if_help or (not defined $input or not defined $output or not defined $count));
###################################################################### 主程序

open INPUT, "gzip -cd $input|";
open OUTPUT, "|gzip > $output";
my $reads = 0;
while(my $line1 = <INPUT>)
{
    $reads++;
    my $line2 = <INPUT>;
    my $line3 = <INPUT>;
    my $line4 = <INPUT>;
    print OUTPUT "$line1$line2$line3$line4";
    last if($reads >= $count);
}
close INPUT;
close OUTPUT;


###################################################################### 子函数

 

sub help{
    my $info = "
Program: fastq cut， 截取fastq文件前n条
         多用于文件损坏的场景
Version: 2019-09-06
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input/-i    fastq文件,示例： a_R1.fastq.gz
         --output/-o   输出文件，示例；b_R1.fastq.gz
         --count/-c    取前n行，示例：1000
         --help/-h           查看帮助文档
    \n";
    return $info;
}

