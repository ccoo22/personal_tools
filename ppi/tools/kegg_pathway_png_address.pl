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
my ($kegg_enrichment, $color_file, $output_file, $if_help);
GetOptions(
    "kegg_enrichment|k=s" => \$kegg_enrichment,
    "color|c=s"           => \$color_file,
    "output_file|o=s"     => \$output_file,
    "help|h"              => \$if_help,
);
die help() if(defined $if_help or (not defined $kegg_enrichment or not defined $output_file));
$color_file = "" if(not defined $color_file);
###################################################################### 主程序
my %hashDiffPahtway = read_kegg_enrichment($kegg_enrichment, 0.05); # 读取富集分析结果,p值小于0.05
my %hashColor       = read_color_file($color_file); # 读取基因上下调状态

open ADD, ">$output_file";
foreach my $pathway_id(sort keys %hashDiffPahtway)
{
	my @colors;
	foreach my $gene_name(sort keys %{$hashDiffPahtway{$pathway_id}})
	{
		my $gene_color = 'red';
		   $gene_color = 'red'  if(exists $hashColor{$gene_name} and $hashColor{$gene_name} eq 'up');
		   $gene_color = 'blue' if(exists $hashColor{$gene_name} and $hashColor{$gene_name} eq 'down');
		push @colors, "$gene_name%09$gene_color"; # 基因设成红色
	}
	my $color = join "/", @colors;
	print ADD "http://www.genome.jp/kegg-bin/show_pathway?$pathway_id/$color\n";
}
close ADD;


###################################################################### 子函数
sub read_color_file{
	my $color_file = shift @_;
	my %hashColor;

	return %hashColor if($color_file eq '');

	open COLOR, $color_file;
	while(<COLOR>)
	{
		$_=~s/[\r\n]//g;
		my ($gene, $regular) = split /\t/, $_;
		$regular = 'up' if(not defined $regular);
		$hashColor{$gene} = lc($regular);
	}
	close COLOR;
	return %hashColor;

}
sub read_kegg_enrichment{
	my $kegg_enrichment = shift @_;
	my $pvalue_cutoff   = shift @_;

	my %hashDiffPahtway;
	open ENRICHMENT, $kegg_enrichment;
	my $line1 = <ENRICHMENT>;
	   $line1 =~ s/[\r\n]//g;
	my @heads = split /\t/, $line1;

	while(<ENRICHMENT>)
	{
		next if($_!~/\w/);
		$_=~s/[\r\n]//g;
		my @datas = split /\t/, $_;
		my %hashTmp = map{  ($heads[$_], $datas[$_])  } (0..$#heads);

		my $pathway_id     = $hashTmp{'ID'};
		my $pvalue         = $hashTmp{'pvalue'};
		my $gene_name_list = $hashTmp{'geneID'};

		next if($pvalue > $pvalue_cutoff);
		foreach my $gene_name(split /\//, $gene_name_list)
		{
			$hashDiffPahtway{$pathway_id}{$gene_name}++;
		}
	}
	close ENRICHMENT;
	return %hashDiffPahtway;
}

sub help{
    my $info = "
Program: get kegg pathway download address， 根据富集分析结果，获取下载地址,基因设成红色
Version: 2019-05-10
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --kegg_enrichment/-k [必填] kegg_enrichment.xls 结果文件。kegg通路富集分析原始结果
         --output_file/-o     [必填] 结果输出文件
         --color/-c           [选填] 基因颜色文件，第一列基因，第二列up/down, (up=read, down=blue)

         --help/-h         查看帮助文档
    \n";
    return $info;
}
