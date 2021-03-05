package environment;

use strict;
use warnings;
use Exporter;
use vars qw(@EXPORT_OK @ISA);

@ISA = 'Exporter';
@EXPORT_OK = qw(get_perl_package get_r_library get_software get_database get_docker);

sub get_perl_package{
	my $file = shift;
	my @array;
	open FILE, $file;
	while(my $line = <FILE>){
		$line =~ s/#.+//g;
		$line =~ s/\s//g;
		push @array, $line if($line =~ /\w/);
	}
	close FILE;
	return @array;
}

sub get_r_library{
	my $file = shift;
	my @array;
	my ($r_bin, $r_lib);
	open FILE, $file;
	while(my $line = <FILE>){
		$line =~ s/#.+//g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		if($line =~ /::CLEAR/){
			($r_bin, $r_lib) = ();
		}elsif($line =~ /::R/){
			$r_bin = (split /\s+/, $line, 2)[1];
		}elsif($line =~ /::LIB/){
			$r_lib = (split /\s+/, $line, 2)[1];
		}else{
			push @array, "$line,$r_bin,$r_lib" if($line =~ /\w/);
		}
	}
	close FILE;
	return @array;
}

sub get_software{
	my $file = shift;
	my %hash;
	open FILE, $file;
	while(my $line = <FILE>){
		$line =~ s/#.+//g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		if($line =~ /(\S+)\s+argv=([^;]*);rtn=(\d+)\s+(.+)/){
			%{$hash{$1}} = ("argv"=>$2, "rtn"=>$3, "cmd"=>$4);
		}
	}
	close FILE;
	return %hash;
}

sub get_database{
	my $file = shift;
	my %hash;
	open FILE, $file;
	while(my $line = <FILE>){
		$line =~ s/#.+//g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		my ($a, $b) = split /\s+/, $line, 2;
		next if(not defined $b);
		$hash{$a} = $b;
	}
	close FILE;
	return %hash;
}

sub get_docker{
	my $file = shift;
	my %hash;
	open FILE, $file;
	while(my $line = <FILE>){
		$line =~ s/#.+//g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		my ($a, $b) = split /\s+/, $line, 2;
		next if(not defined $b);
		$hash{$a} = $b;
	}
	close FILE;
	return %hash;
}



1