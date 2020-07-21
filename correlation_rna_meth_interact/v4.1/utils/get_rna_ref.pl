#!/usr/bin/env perl
# use strict;
# use warnings;
use Getopt::Long;

my $ref = "/home/wangly/scripts/correlation/database/hg19/hg19_refGene.txt";

my ($diff, $out, $type, $help);

GetOptions(
	'diff|i=s'      => \$diff,
	'out|o=s'      => \$out,
	'type|t=s'      => \$type,
	"help|h!" 		=> \$help,
);
die help() if (defined $help or not defined $diff);

my ($id, %list, %hash, $gene, $nm);

open DIFF, $diff;
open OUT, qq{>$out};

my $line1 = <DIFF>;
   $line1 =~ s/[\r\n]//g;
my @heads = split /\t/, $line1;
my $delete = shift @heads;
my @title = qw/id	chr	start	strand/;
@heads = (@title, @heads);
print OUT (join "\t", @heads) . "\n"; 

while (<DIFF>){
	$_=~s/[\r\n]//g;
	my @array = split /\t/;
	$id = shift @array;
	$list{$id} = join("\t",@array);
}

open REF, $ref;
while (<REF>){
	$_=~s/[\r\n]//g;
	my ($chr, $strand, $start) = (split /\t/, $_)[2,3,4];
	$chr =~ s/chr//g;
	my @arr = split /\t/;
	$nm = $arr[1];
	$gene = $arr[12];
	
	if($type eq "gene" && not exists $hash{$gene}){
		$hash{$gene}{'chr'}	= $chr;
		$hash{$gene}{'strand'} = $strand;
		$hash{$gene}{'start'} = $start; 
	}
	elsif($type eq "gene" && exists $hash{$gene}){
		if($strand eq "+" && $start < $hash{$gene}{'start'}){
			$hash{$gene}{'start'} = $start;
			$hash{$gene}{'strand'} = $strand;
			$hash{$gene}{'chr'}	= $chr;
		}
		if($strand eq "-" && $start > $hash{$gene}{'start'}){
			$hash{$gene}{'start'} = $start;
			$hash{$gene}{'chr'}	= $chr;
			$hash{$gene}{'strand'} = $strand;
		}
	}
	
	if($type eq "transcript" && not exists $hash{$nm}){
		$hash{$nm}{'chr'} = $chr;
		$hash{$nm}{'strand'} = $strand;
		$hash{$nm}{'start'} = $start;
	}
	elsif($type eq "transcript" && exists $hash{$nm}){
		if($strand eq "+" && $start < $hash{$nm}{'start'}){
			$hash{$nm}{'start'} = $start;
			$hash{$nm}{'strand'} = $strand;
			$hash{$nm}{'chr'}	= $chr;
		}
		if($strand eq "-" && $start > $hash{$nm}{'start'}){
			$hash{$nm}{'start'} = $start;
			$hash{$nm}{'strand'} = $strand;
			$hash{$nm}{'chr'}	= $chr;
		}
	}
}
foreach my $gene(sort keys %hash){ 
	print OUT "$gene\t$hash{$gene}{'chr'}\t$hash{$gene}{'start'}\t$hash{$gene}{'strand'}\t$list{$gene}\n" if($type eq "gene" && exists $list{$gene});
}
foreach my $nm(sort keys %hash){ 
	print OUT "$nm\t$hash{$nm}{'chr'}\t$hash{$nm}{'start'}\t$hash{$nm}{'strand'}\t$list{$nm}\n" if($type eq "transcript" && exists $list{$nm});
}

close OUT;
close REF;
close DIFF;

sub help
{
my $info = "
Usage: perl $0  -i diff -o outfilename -t type
 Options: 必填
	--diff/-i		行为基因列为样本的归一化count文件
				example：/home/wangly/scripts/correlation/v2.0/example/diff_gene_expression.xls
				
	--outfilename/-o	输出文件名			
	--type/-t		gene or transcript
				
 Options: 可选
	--help/-h		查看帮助文档
	\n";
	return $info;
}