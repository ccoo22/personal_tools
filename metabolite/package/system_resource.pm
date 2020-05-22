package system_resource;

use strict;
use warnings;
use Exporter;
use vars qw(@EXPORT_OK @ISA);

@ISA = 'Exporter';
@EXPORT_OK = qw(check_system_resource set_thread);

sub set_thread{
	my ($thread, $sample, $expect) = @_;
	my $t_sample = 1;
	while($t_sample < $sample and int($thread/($t_sample+1)) >= $expect){
		$t_sample++;
	}
	my $t_soft = 1;
	$t_soft = int($thread/$t_sample) if($thread > $t_sample);
	return ($t_soft, $t_sample);
}

sub get_cpu_info{
	my $test_time = shift;
	my $processor = `cat /proc/cpuinfo | grep "processor" | wc -l`;
	my ($max, $min, $avg);
	my $cpu_info = `sar -u 1 $test_time`;
	foreach(split /[\r\n]/, $cpu_info){
		next if($_ !~ m/^\d/);
		my ($idle) = $_ =~ /([\d\.]+)$/;
		next if(not defined $idle);
		$max = $idle if(not defined $max or $max < $idle);
		$min = $idle if(not defined $min or $min > $idle);
		$avg += $idle/$test_time;
	}
	my $cpu_max = $processor-int($processor*$min/100);
	my $cpu_min = $processor-int($processor*$max/100);
	my $cpu_avg = $processor-int($processor*$avg/100);
	my $cpu_idle = int($processor*$avg/100);
	return ($cpu_max, $cpu_min, $cpu_avg, $cpu_idle);
}

sub get_memory_info{
	my $memory_info = `free -g`;
	my $free = 0;
	foreach(split /[\r\n]/, $memory_info){
		next if($_ !~ m/^Mem:/);
		my @info = $_ =~ /(\d+)/g;
		$free = $info[0]-$info[1];
	}
	return $free;
}

sub check_system_resource{
	my $thread = shift;
	my $memory = shift;
	my ($cpu_max, $cpu_min, $cpu_avg, $cpu_idle) = get_cpu_info(2);
	my ($memory_free) = get_memory_info();
	my $err = "";
	$err .= "[WARNING] 空闲 CPU = $cpu_idle < $thread\n" if($cpu_idle < $thread);
	$err .= "[WARNING] 剩余内存 = $memory_free GB < $memory\n" if($memory_free < $memory);
	if($err ne ""){
		print "$err\n是否继续运行?(y/n)";
		chomp(my $type = <STDIN>);
		exit if($type !~ /y/i);
	}
}




1