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
my ($kegg_enrichment, $color_file, $output_file, $add_address_to_excel, $if_help);
GetOptions(
    "kegg_enrichment|k=s"     => \$kegg_enrichment,
    "color|c=s"               => \$color_file,
    "output_file|o=s"         => \$output_file,
    "add_address_to_excel!" => \$add_address_to_excel,
    "help|h"                  => \$if_help,
);
die help() if(defined $if_help or (not defined $kegg_enrichment or not defined $output_file));
$color_file = "" if(not defined $color_file);
###################################################################### 主程序
my %hashDiffPahtway = read_kegg_enrichment($kegg_enrichment, 0.05); # 读取富集分析结果,p值小于0.05
my %hashColor       = read_color_file($color_file); # 读取基因上下调状态

# 制作下载地址
my %hashAdd;
open ADD, ">$output_file";
foreach my $pathway_id(sort keys %hashDiffPahtway)
{
	# my $add = "http://www.genome.jp/kegg-bin/show_pathway?$pathway_id"; # 下载地址
	my $add = "http://www.kegg.jp/kegg-bin/show_pathway?$pathway_id"; # 下载地址
	my @lost_genes;
	foreach my $gene_name(sort keys %{$hashDiffPahtway{$pathway_id}})
	{
		my $gene_color = 'red';
		   $gene_color = 'red'  if(exists $hashColor{$gene_name} and $hashColor{$gene_name} eq 'up');
		   $gene_color = 'blue' if(exists $hashColor{$gene_name} and $hashColor{$gene_name} eq 'down');
		my $add_tmp = "$add/$gene_name%09$gene_color";

		# kegg 提交的地址长度不能超过3040（大约），否则直接拒绝访问
		if(length($add_tmp) > 3040)
		{	
			push @lost_genes, $gene_name;
			next;
		}
		$add = $add_tmp;
 
	}
	print "[Warnings] 下载地址中，从通路 [$pathway_id] 里删除了后面的部分基因，原因：kegg的访问地址字符长度不能超过3040左右，否则直接拒绝访问\n通路[$pathway_id] 删除的基因：@lost_genes\n" if(@lost_genes > 0);
	print ADD "$pathway_id\t$add\n";
	$hashAdd{$pathway_id} = $add;
}
close ADD;

exit if(not defined $add_address_to_excel);
# 下载地址加入到kegg_enrichment文件中

open KEGG_ENRICHMENT, $kegg_enrichment;
open KEGG_ENRICHMENT_ADD, ">$kegg_enrichment.add_address.txt";
my $line1 = <KEGG_ENRICHMENT>;
   $line1 =~ s/[\r\n]//g;
my @heads = split /\t/, $line1;
print KEGG_ENRICHMENT_ADD "$line1\tkegg_link\n";  # 更新表头
while(<KEGG_ENRICHMENT>)
{
	next if($_!~/\w/);
	$_=~s/[\r\n]//g;
	my @datas = split /\t/, $_;
	my %hashTmp = map{  ($heads[$_], $datas[$_])  } (0..$#heads);

	my $pathway_id     = $hashTmp{'ID'};
	my $pathway_address = exists $hashAdd{$pathway_id} ? $hashAdd{$pathway_id} : " ";
	print KEGG_ENRICHMENT_ADD "$_\t$pathway_address\n";
}
close KEGG_ENRICHMENT;
close KEGG_ENRICHMENT_ADD;
system("mv $kegg_enrichment.add_address.txt $kegg_enrichment");

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

         --kegg_enrichment/-k   [必填] kegg_enrichment.xls 结果文件。kegg通路富集分析原始结果
         --output_file/-o       [必填] 结果输出文件
         --color/-c             [选填] 基因颜色文件，第一列基因，第二列up/down, (up=read, down=blue)
         --add_address_to_excel [选填] 把kegg图像地址加入kegg_enrichment文件的末尾，供客户随时查看

         --help/-h         查看帮助文档
    \n";
    return $info;
}
