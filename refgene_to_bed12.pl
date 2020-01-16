# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $genePredToBed = "/home/genesky/software/ucsc/bin/genePredToBed";

# 检测 -> 脚本输入
my ($refgene, $output, $if_help);
GetOptions(
    "refgene|i=s"       => \$refgene,
    "output|o=s"        => \$output,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $refgene or not defined $output ));

###################################################################### 主程序
system("awk '{\$1=\"\"; print \$0}' $refgene |sed  's/^\\s//'| $genePredToBed  /dev/stdin $output");
 

###################################################################### 子函数

sub help{
    my $info = "
Program: refGene to BED12, refGene(UCSC) 转 BED12(UCSC)
Version: 2019-08-01
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --refgene/-i     [必填] refGene文件, 示例：/home/genesky/database/ucsc/hg19/gene/hg19_refGene.txt
         --output/-o      [必填] 结果文件, 示例：hg19_refGene.bed12 
         --help/-h        查看帮助文档
    \n";
    return $info;
}

