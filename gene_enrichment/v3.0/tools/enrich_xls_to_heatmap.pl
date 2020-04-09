#!/usr/bin/env perl

my ($enrich_xls) = @ARGV;


my @pathways = ();
my @genes    = ();
my %meta     = ();

open EOX, $enrich_xls or die "Can't oepn $enrich_xls!\n";
while (<EOX>) {
	chomp;
	next if /^ID/;
	my @arr = split /\t/;
	next if $arr[4] > 0.05;
	push @pathways, $arr[1] if not $arr[1] ~~ @pathways;
	foreach my $x (split /\//, $arr[7]) {
		push @genes, $x if not $x ~~ @genes;
		$meta{$arr[1]}{$x} = 1;
	}
}
close EOX;

my $title = join "\t", @genes;
print qq{ID\t$title\n};
foreach my $x (@pathways) {
	my @res = ();
	foreach my $y (@genes) {
		my $val = exists $meta{$x}{$y} ? 1 : 0;
		push @res, $val;
	}
	my $line = join "\t", @res;
	print qq{$x\t$line\n};
}