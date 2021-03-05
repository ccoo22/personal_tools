$| = 1; print "\n"."#"x(30)."
系统环境测试
Version: 2019-10-15
Contact: 042 王立
\n";
# Perl 系统包
use warnings;
use strict;
use Exporter;
use File::Spec;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# Perl 自定义包
use lib SCRIPT_DIR."/../package";
use system_time qw(time_to_datetime);
use environment qw(get_perl_package get_r_library get_software get_database get_docker);

# 变量
my $PERL_PACKAGE = SCRIPT_DIR."/perl_package.txt";
my $R_LIBRARY = SCRIPT_DIR."/r_library.txt";
my $SOFTWARE = SCRIPT_DIR."/software.txt";
my $DATABASE = SCRIPT_DIR."/database.txt";
my $DOCKER = SCRIPT_DIR."/docker.txt";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;

###################################################################### 初始化
my @perl_package = get_perl_package($PERL_PACKAGE);
my @r_library = get_r_library($R_LIBRARY);
my %software = get_software($SOFTWARE);
my %database = get_database($DATABASE);
my %docker = get_docker($DOCKER);
my $DATA_TIME = time_to_datetime(time);
print "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
";

###################################################################### 主程序
my $isok = 1;
$isok *= test_perl_package(@perl_package);
$isok *= test_r_library(@r_library);
$isok *= test_software(\%software);
$isok *= test_database(\%database);
$isok *= test_docker(\%docker);
if($isok == 1){
	print "\n[OK] 系统环境测试通过！\n";
}else{
	print "\n[ERROR] 系统环境错误！\n";
}
###################################################################### 子程序
sub test_perl_package{
	my $isok = 1;
	foreach(@_){
		print "[RUNTIME] Test perl package '$_' ... ";
		my $test = `perl -e 'use $_;' 2>&1`;
		if($test ne ""){
			print "ERROR\n";
			$isok = 0;
		}else{
			print "OK\n";
		}
	}
	return $isok;
}

sub test_r_library{
	my $isok = 1;
	foreach(@_){
		my ($pak, $r_bin, $r_lib) = split /,/, $_;
		my $lib = ""; $lib = ".libPaths(\"$r_lib\");" if(defined $r_lib);
		print "[RUNTIME] Test R library '$_' ... ";
		my $test = system "$r_bin -e '$lib library($pak);' > /dev/null 2>&1";
		if($test ne 0){
			print "ERROR\n";
			$isok = 0;
		}else{
			print "OK\n";
		}
	}
	return $isok;
}

sub test_software{
	my $hash = shift;
	my $isok = 1;
	foreach my $key(sort{$a cmp $b} keys %{$hash}){
		print "[RUNTIME] Test Software '$key', '".$hash->{$key}{"cmd"}." ".$hash->{$key}{"argv"}."' ... ";
		my $test = system $hash->{$key}{"cmd"}." ".$hash->{$key}{"argv"}." > /dev/null 2>&1";
		if($test ne $hash->{$key}{"rtn"}){
			print "ERROR (rtn=$test)\n";
			$isok = 0;
		}else{
			print "OK\n";
		}
	}
	return $isok;
}

sub test_database{
	my $hash = shift;
	my $isok = 1;
	foreach my $key(sort{$a cmp $b} keys %{$hash}){
		print "[RUNTIME] Test Database '$key', '$hash->{$key}' ... ";
		if(not -e $hash->{$key}){
			print "ERROR\n";
			$isok = 0;
		}else{
			print "OK\n";
		}
	}
	return $isok;
}

sub test_docker{
	my $hash = shift;
	my $isok = 1;
	foreach my $key(sort{$a cmp $b} keys %{$hash}){
		print "[RUNTIME] Test Docker '$key', '$hash->{$key}' ... ";
		my $test = `docker images $hash->{$key} 2>&1`;
		my @split = split /[\r\n]/, $test;
		my $key = (split /:/, $hash->{$key})[0];
		if(not exists $split[1] or $split[1] !~ /^$key/){
			print "ERROR\n";
			$isok = 0;
		}else{
			print "OK\n";
		}
	}
	return $isok;
}

