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
my ($bed, $dict, $output, $if_help);
GetOptions(
    "bed|b=s"           => \$bed,
    "dict|d=s"          => \$dict,
    "output|o=s"        => \$output,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $bed or not defined $dict or not defined $output ));

###################################################################### 主程序
system("cp $dict $output");
open BED, $bed;
open OUT, ">>$output";
while(<BED>)
{
	$_=~s/[\r\n]//g;
	next if($_ =~ /^#/);
	my ($chr, $start, $end, $gene, $score, $strand) = split /\t/, $_;
	next if($chr =~ /_/);
	print OUT "$chr\t$start\t$end\t$strand\t$gene\n";
}
close OUT;
close BED;

###################################################################### 子函数

sub help{
    my $info = "
Program: abstract rrna intervars from UCSC， 从UCSC rRNA中提取rRNA区域
Version: 2019-04-20
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --bed/-b         [必填] bed文件, rRNA_UCSC_hg38.bed。示例：下载自UCSC TableBrowser。group:all tables, table:rmsk, and filter to 'repClass (does match) rRNA'。
         --dict/-d        [必填] 参考基因组头文件, 示例：/home/genesky/database/ucsc/hg38/genome/hg38.dict
         --output/-o      [必填] 结果文件, 示例：rRNA_UCSC_hg38.ribosomal_intervals
         --help/-h        查看帮助文档
    \n";
    return $info;
}

