use strict;
use warnings;
die "Usage:perl $0 Inputfile\n" if(@ARGV!=1);
# Annovar建库必须保证输入数据是按照染色体->起始位置从小到大依次排序构成的。本代码用于检测
my $file=shift @ARGV;
open IN,$file;
my %hashchr;
my $lastchr="";
my $lastposition=0;
my $row = 0;
while(<IN>){
	$_=~s/[\r\n]//g;
        $row++;
	my ($chr,$position,$tmp)=split /\t/,$_,3;
	if($chr ne $lastchr){
		die "Chromosome Error:$chr,$position!\n" if(exists($hashchr{$chr}));
		$lastchr=$chr;
		$lastposition=$position;
		next;
	}
	die "Not Order:$chr,$position!(file row = $row)\n" if($position<$lastposition);
	$lastposition=$position;
}
