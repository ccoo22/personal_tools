# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $bowtie2      = "/home/genesky/software/bowtie2/2.3.4.3/bowtie2";


# 检测 -> 脚本输入
my ($genome, $r1, $r2, $output, $bowtie2_para, $if_help);
GetOptions(
    "genome|r=s"       => \$genome,
    "r1|1=s"           => \$r1,
    "r2|2=s"           => \$r2,
    "output|o=s"       => \$output,
    "bowtie2_para|p=s" => \$bowtie2_para,
    "help|h"           => \$if_help,
);
die help() if(defined $if_help or (not defined $genome or not defined $output or not defined $r1) );
$bowtie2_para = "" if(not defined $bowtie2_para);
my $output_para = ($output eq 'std') ? "" : "-S $output";
###################################################################### 主程序

my $fastq_para = (defined $r2) ? "-1 $r1 -2 $r2" : "-U $r1";

system("$bowtie2 -x $genome $fastq_para $output_para $bowtie2_para ");


###################################################################### 子函数


sub help{
    my $info = "
Program: bowtie2 
Version: 2019-03-25
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --genome/-r        [必填] 参考序列index文件, 示例：/home/pub/database/Human/hg19/bowtie2_db/hg19.fa
         --r1/-1            [必填] 输入fastq r1文件 
         --r2/-2            [选填] 输入fastq r2文件
         --output/-o        [必填] 输出文件sample.sam
         --bowtie2_para/-p  其他bowtie2参数，注意添加单引号，示例：'--np 0 --n-ceil L,0,0.5 -p 10 -a'
         --help/-h          查看帮助文档
    \n";
    return $info;
}

