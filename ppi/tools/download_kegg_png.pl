
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;


my ($input, $out_dir,  $help);

GetOptions(
	'input|i=s'     => \$input,
	'out_dir|o=s'    => \$out_dir,
	'help|h!'       => \$help
);

if ($help or not $input or not $out_dir) {
	usage();
	exit();
}

sub usage
{
my $help  =<<EOF;
Usage: perl  -input input -out_dir output
	-input     -i  input file [required]
	-out_dir   -o  out_dir file [required]
	-help      -h  print help message
EOF

print $help;

}


my $max_threads = 10;
my $pm = Parallel::ForkManager->new($max_threads);

make_dir();

my @html_urls = prase_txt($input);

download();

system qq{cat $out_dir/html/*html | grep "tmp/mark_pathway" | sed 's/<img src="//' |   sed 's/" name="pathwayimage" name="pathwayimage" usemap="#mapdata" border="0" \\/>//' |   awk '{print "http://www.genome.jp"\$0}' > $out_dir/png.urls.txt};

my @png_urls = prase_png_txt(qq{$out_dir/png.urls.txt});

download_pathway();


sub make_dir
{
	system qq{mkdir -p $out_dir}      if not -d $out_dir;
	system qq{mkdir -p $out_dir/png}  if not -d qq{$out_dir/png};
	system qq{mkdir -p $out_dir/html} if not -d qq{$out_dir/html};
}



sub prase_txt
{
	my $txt  = shift;
	my @urls = ();

	open TXT, $txt or die "Can't open $txt!\n";
	while (<TXT>) {
		s/\s+$//g;
		my @arr = split /\s+/;
		if (scalar @arr > 2) {
			push @urls, $arr[2];
		} else {
			push @urls, $arr[0];
		}
		
	}
	close TXT;

	return @urls;
}


sub download
{
	my @rest_htmls = ();
	do {
		@rest_htmls = download_html(\@html_urls, $out_dir);
	} while (scalar @rest_htmls != 0)


}


sub download_html
{
	my $urls     = shift;
	my $out_dir  = shift;

	my $max_threads = 10;
	my $pm          = Parallel::ForkManager->new($max_threads);

	
	foreach my $x (@{$urls}) {
		my $pid = $pm->start and next;
		my ($name) = $x =~ /show_pathway\?(\w+)/;
		system qq{wget  -q -O $out_dir/html/$name.html "$x"};
		$pm->finish;
	}
	$pm->wait_all_children;

	my @rest     = ();
	foreach my $x (@{$urls}) {
		my ($name) = $x =~ /show_pathway\?(\w+)/;
		next if -e qq{$out_dir/html/$name.html} and end_html(qq{$out_dir/html/$name.html}) == 1;
		push @rest, $x;
	}

	return @rest;
}



sub prase_png_txt
{
	my $txt  = shift;
	my @urls = ();

	open TXT, $txt or die "Can't open $txt!\n";
	while (<TXT>) {
		s/\s+$//g;


		push @urls, $_;
	}
	close TXT;

	return @urls;
}


sub download_pathway
{
	my @rest_htmls = ();
	do {
		@rest_htmls = download_png(\@png_urls, $out_dir);
	} while (scalar @rest_htmls != 0)


}

sub download_png
{
	my $urls     = shift;
	my $out_dir  = shift;

	my $max_threads = 10;
	my $pm          = Parallel::ForkManager->new($max_threads);

	
	foreach my $x (@{$urls}) {
		my $pid = $pm->start and next;
		my ($name) = $x =~ /(\w+)\.png/;
		system qq{wget  -q -O $out_dir/png/$name.png "$x"};
		$pm->finish;
	}
	$pm->wait_all_children;

	my @rest     = ();
	foreach my $x (@{$urls}) {
		my ($name) = $x =~ /(\w+)\.png/;
		next if -e qq{$out_dir/png/$name.png} and end_png(qq{$out_dir/png/$name.png}) == 1;
		push @rest, $x;
	}

	return @rest;
}




sub end_html
{
	my $html = shift;
	local $/ = undef;
	open HTML, $html;
	my $context = <HTML>;
	close HTML;
	return $context =~ /<\/html>/ ? 1 : 0;
}




sub end_png
{
    my $png = shift;
    my $buff;
    my $str;
    open PNG, $png or die "Can't open $png";
    binmode(PNG);
    while(read(PNG, $buff, 1)){
        my $hex = unpack("H*", $buff);
        $str .= $hex;
    }
    close PNG;
    my $line = substr $str, -8;

    return $line  eq qq{ae426082} ? 1 : 0;
}

