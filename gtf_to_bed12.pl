# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $gtfToGenePred = "/home/genesky/software/ucsc/bin/gtfToGenePred";
my $genePredToBed = "/home/genesky/software/ucsc/bin/genePredToBed";

# 检测 -> 脚本输入
my ($gtf, $output, $if_help);
GetOptions(
    "gtf|i=s"           => \$gtf,
    "output|o=s"        => \$output,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $gtf or not defined $output ));

###################################################################### 主程序
system("grep -v rRNA $gtf | $gtfToGenePred /dev/stdin /dev/stdout| $genePredToBed /dev/stdin $output");


###################################################################### 子函数
sub parse_attrbutes{
	my $attributes = shift @_;
	my %hashAttributes;
	foreach my $attribute(split /;/, $attributes)
	{
		$attribute=~ s/^\s+//;
		$attribute=~ s/\s+$//;
		my ($name, $value) = split /\s+/, $attribute;
		$value =~ s/^\"//;
		$value =~ s/\"$//;
		$hashAttributes{$name} = $value;
	}
	return %hashAttributes;
}

sub help{
    my $info = "
Program: gtf to BED12, gtf 转 BED12(UCSC)
Version: 2019-02-21
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --gtf/-i         [必填] gtf文件, 示例：gencode.v30.annotation.gtf
         --output/-o      [必填] 结果文件, 示例：gencode.v30.annotation.bed12(注意：过程中，rRNA去掉了。另外：会产生大约208067个转录本)
         --help/-h        查看帮助文档
    \n";
    return $info;
}

