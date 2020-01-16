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
my ($gtf, $dict, $output, $if_help);
GetOptions(
    "gtf|i=s"           => \$gtf,
    "dict|d=s"           => \$dict,
    "output|o=s"        => \$output,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $gtf or not defined $dict or not defined $output ));

###################################################################### 主程序
system("cp $dict $output");
open GTF, $gtf;
open OUT, ">>$output";
while(<GTF>)
{
	$_=~s/[\r\n]//g;
	next if($_ =~ /^#/);
	my ($chr, $source, $region_type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
	next if($region_type ne 'gene');
	next if($attributes !~ /gene_type\s\"rRNA\";/ and $attributes !~ /gene_type\s\"rRNA_pseudogene\";/);
	my %hashAttributes = parse_attrbutes($attributes);
	my $gene_name = $hashAttributes{'gene_name'};
	print OUT "$chr\t$start\t$end\t$strand\t$gene_name\n";
}
close OUT;
close GTF;

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
Program: abstract rrna intervars， 从gtf中提取rRNA区域
Version: 2019-02-21
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --gtf/-i         [必填] gtf文件, 示例：gencode.v30.annotation.gtf
         --dict/-d        [必填] 参考基因组头文件, 示例：/home/genesky/database/ucsc/hg38/genome/hg38.dict
         --output/-o      [必填] 结果文件, 示例：gencode.v30.annotation.ribosomal_intervals
         --help/-h        查看帮助文档
    \n";
    return $info;
}

