package system_time;

use strict;
use warnings;
use Exporter;
use vars qw(@EXPORT_OK @ISA);

@ISA = 'Exporter';
@EXPORT_OK = qw(time_to_datetime);

sub time_to_datetime{
	my $time = shift;
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
	$year += 1900;
	$mon += 1;
	return "$year-$mon-$day $hour:$min:$sec";
}




1