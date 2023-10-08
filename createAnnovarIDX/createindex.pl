#!/usr/bin perl;
use warnings;
use strict;

die "$0 <Annovar Database File> <BIN Size 100>" unless @ARGV == 2;
my $input_file = $ARGV[0];# 库名
my $bin_size = $ARGV[1];# 跨度
 
if (!-e $input_file) {
	die "$input_file not found\n";
}

my $file_size = -s $input_file;

my %index;
open(my $in, "<", $input_file) or die "Couldn't open $input_file for indexing\n";

my $previous_file_position = tell $in;

while (my $ln = <$in>) {
	next if($ln=~/^#/);
	#Check input file. Some are (chr,start,stop) and others are (id,chr,start,stop).
	#If you have the latter you'll need to change the next line to account for the id column   
	my ($chr,$start,$stop) = split "\t", $ln;
	my $bin_start = int($start/$bin_size) * $bin_size;
	my $current_file_position = tell $in;

	if (!exists $index{$chr}->{$bin_start}) {
		$index{$chr}->{$bin_start} = [$previous_file_position, $current_file_position];
	}
	else{
		$index{$chr}->{$bin_start}->[1] = $current_file_position;
	}
	
	$previous_file_position = $current_file_position;
}

close $in;
open OUT,">$input_file.idx";
print OUT "#BIN\t$bin_size\t$file_size\n";
foreach my $chr ((1,10..19,2,20,21,22,3..9,"MT","X","Y")){ #Ordered array to match other Annovar idx files
	foreach my $index_region (sort keys %{$index{$chr}}){
		my $start	= $index{$chr}->{$index_region}->[0];
		my $stop	= $index{$chr}->{$index_region}->[1];
		print OUT "$chr\t$index_region\t$start\t$stop\n";
	}
}