$| = 1; print "\n"."#"x(30)."
metabolite 分析
Version: v1.0 2020-6-30
Contact: 129 甘斌
\n";
# Perl 系统包
use warnings;
use strict;
use Cwd 'abs_path'; 
use File::Spec;
use Getopt::Long;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];

# Perl 自定义包
use lib SCRIPT_DIR."/environment";
use lib SCRIPT_DIR."/package";
use system_resource qw(check_system_resource set_thread);
use system_time qw(time_to_datetime);
use environment qw(get_software get_database);

# 变量
my $SOFTWARE = SCRIPT_DIR."/environment/software.txt";
my $DATABASE = SCRIPT_DIR."/environment/database.txt";
my $DEFAULT_THREAD = 10;
my $DEFAULT_MEMORY = 1;
my $SOFT_METABOLITE_RUN = SCRIPT_DIR."/utils/metabolite_run.pl";
my $SOFT_REPORT = SCRIPT_DIR."/utils/report/report_metabolite.pl";
my $SOFT_HTML = "/home/wangly/scripts/wts/html_markdown/v3.2/pipeline.pl";
my $SOFT_HTML_TEMPLATE_ID = "t_wts_v2.2";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($config_ini, $thread, $force_cover, $if_help);
GetOptions(
	"config-ini|c=s" => \$config_ini,
	"force-cover:s" => \$force_cover,
	"thread=i" => \$thread,
	"help|h" => \$if_help,
);
die "
Options: 必填

        --config-ini/-c            流程配置文件

Options: 可选

        --force-cover [sample,]    无参数, 删除所有原始结果文件, 完全覆盖;
                                   指定sample, 删除特定样本结果文件, 多个样本间','分隔
        --thread                   使用线程数 (default: $DEFAULT_THREAD)
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $config_ini);
$thread = $DEFAULT_THREAD if (not defined $thread);
$DEFAULT_MEMORY *= $thread/$DEFAULT_THREAD;
my $Ref = "";
###################################################################### 初始化
my %software = get_software($SOFTWARE);
my %database = get_database($DATABASE);
my %config = read_config($config_ini);

# 运行参数
check_system_resource($thread, $DEFAULT_MEMORY);
my $DATA_TIME = time_to_datetime(time);
print "
---
Command: perl ".abs_path($0)." $ARGV_INFO
---
Start time: $DATA_TIME
[SET] CPU 总消耗 $thread, 内存总消耗 $DEFAULT_MEMORY
\n";
foreach my $key(sort{$config{$a}{"sort"} <=> $config{$b}{"sort"}} keys %config){
	next if(not defined $config{$key}{"value"});
	print "[SET] $key = ".$config{$key}{"value"}."\n";
}
my $set_fc = "";
if(defined $force_cover){
	$set_fc = "--force-cover $force_cover";
	print "[SET] 数据覆盖 = $set_fc\n";
}

my $set_normalized = ($config{'Already_Normalized'}{'value'} eq 'TRUE') ? "--already_normalized" : "";

###################################################################### 主程序
# 建立目录
my $Sample_Group = $config{"Sample_Group"}{"value"};
my $Metabolite = $config{"Metabolite"}{"value"};
my $Backup_Dir = $config{"Backup_Dir"}{"value"};
my $Report_Dir = $config{"Report_Dir"}{"value"};
foreach("metabolite", "log"
){
	system "mkdir -p $Backup_Dir/$_";
}

$DATA_TIME =~ s/\D+/\_/g;
my $log_file = "$Backup_Dir/log/$DATA_TIME.txt";

print "\nLog save to: $log_file\n";

# 软件 / 数据库 设置
my $config_SOFT_METABOLITE_RUN = "";
foreach("R_Lib", "compound_db", "syn_nms", "kegg_db"){ $config_SOFT_METABOLITE_RUN .= " --$_ '".$database{$_}."' "; }
foreach("Rscript"){$config_SOFT_METABOLITE_RUN .= " --$_ '".$software{$_}{"cmd"}."' ";}


my @commands;
prepare_commond("metabolite 代谢物分析",  "perl $SOFT_METABOLITE_RUN -c $config_ini -o $Backup_Dir/metabolite -m $Metabolite -s $Sample_Group  --no-check --thread $thread $config_SOFT_METABOLITE_RUN $set_fc $set_normalized", "$Backup_Dir/metabolite/metabolite_run.log.$DATA_TIME.txt", 
                                      'DATABASE', "R_Lib", "compound_db", "syn_nms", "kegg_db",
                                      "SOFTWARE", "Rscript");
prepare_commond("整理报告目录",  "perl $SOFT_REPORT -c $config_ini -b $Backup_Dir -r $Report_Dir $set_fc", "$Backup_Dir/log/log.report.$DATA_TIME.txt");
# prepare_commond("生成 HTML 报告", "perl $SOFT_HTML -r $Report_Dir -t $SOFT_HTML_TEMPLATE_ID", "$Backup_Dir/log/log.html.$DATA_TIME.txt");


 
open SAVE, ">$Report_Dir/project.txt";
map{ print SAVE "$_ = $config{$_}{'value'}\n"; 1 } ("id", "project", "department", "ref", "sequencing", "method");
close SAVE;



# 运行
my @steps = select_run_steps(\@commands);
foreach my $step(@steps){
    print "[RUNTIME] Step $step : $commands[$step]->[0] ... ";
    system "$commands[$step]->[1] 2>&1 | tee -a  $log_file >  $commands[$step]->[2]";
	my $last_word = `tail -n 1 $log_file`;
	   $last_word=~s/[\r\n]//g;

    if(is_log_ok($log_file)){
        print "$last_word\n";
    }else{
        print "ERROR, 请检查Log: $commands[$step]->[2] （或者检查汇总 Log 日志 $log_file）\n";
    }
}

print "\n[OK] 运行完成\n";
###################################################################### 子程序
sub read_config{
	my $file = shift;
	my $id = 0;
	my %hash;
	open FILE, $file or die "[ERROR] 无法读取配置文件, $file\n";
	while(my $line = <FILE>){
		$line =~ s/#.+//g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		my ($a, $b) = split /\s*=\s*/, $line, 2;
		next if(not defined $b);
		if(exists $hash{$a})
		{
			$hash{$a}{"value"} .= "\n                $b";
		}else
		{
			$hash{$a}{"value"} = $b;
			$hash{$a}{"sort"} = $id;
		}
		$id++;
	}
	close FILE;
	my $err = "";
	foreach("Backup_Dir", "Report_Dir", "Metabolite", "Sample_Group"){
		$err .= "[ERROR] 没有定义 $_ 或 路径不存在, ".$hash{$_}{"value"}."\n" if(not exists $hash{$_} or not -e $hash{$_}{"value"});
	}
	foreach("Already_Normalized", "Compare"){
		$err .= "[ERROR] 没有定义 $_\n" if(not exists $hash{$_});
	}
	die $err if($err ne "");
	return %hash;
}

sub select_run_steps{
	my $command = shift;
	print "\n---\n请选择需要运行的步骤\n---\n\n";
	for my $step(0..@$command-1){
		print " "x(4)."[$step] $command->[$step][0]\n";
	}
	print "\n"." "x(4)."注1：多选输入格式示例 \"7-3,1,8\" (自动排序)。注2：直接输入\"回车\"运行所有步骤。\n";
	print "\n"." "x(4)."请选择[]：";
	my %steps;
	my $input = <STDIN>; $input =~ s/\s//g;
	if($input eq ""){
		map{$steps{$_} = 1} 0..@$command-1;
	}else{
		my @split = split /,/, $input;
		foreach(@split){
			if($_ =~ /^\d+$/){
				$steps{$_} = 1;
			}elsif($_ =~ /^\d+-\d+$/){
				my ($a, $b) = sort{$a <=> $b} split /-/, $_;
				map{$steps{$_} = 1} $a..$b;
			}else{
				die "\n[ERROR] 输入格式无法识别, $_\n";
			}
		}
	}
	my @steps = sort{$a <=> $b} keys %steps;
	print "\n[SET] 运行步骤 : ".(join " ", @steps)."\n";
	foreach(@steps){
		die "\n[ERROR] 输入步骤不存在, $_\n" if($_ >= @$command);
	}
	return @steps;
}

sub is_log_ok{
	my $file = shift;
	my $line;
	if(-e $file){
		open FILE, $file;
		map{$line = $_} <FILE>;
		close FILE;
	}
	return 1 if(defined $line and $line =~ /\[OK\]/);
	return 0;
}


sub prepare_commond{
    my $name         = shift @_;
    my $command      = shift @_;
    my $log          = shift @_;
    my @config_infos = @_;

    # 数据库、带有物种前缀的数据库、带有物种前缀的非必须存在的数据库、软件
    my @key_types = ("DATABASE", "DATABASE_SPECIES", "DATABASE_SPECIES_NOT_ENSSENTIAL", "SOFTWARE");

    # 获取软件、数据库等信息
    my $type  = "";
    my $error = "";
    my $warning = "";
    foreach my $info(@config_infos)
    {
        die "[Error] pipline prepare_commond 错误， 第一个字符必须是 @key_types 中的一个\n " if($type eq "" and (not $info ~~ @key_types));
        if($info ~~ @key_types)
        {
            $type = $info;
            next;
        }

        # 数据库准备
        if($type eq "DATABASE")
        {   
            if(exists $database{$info})
            {
                $command .= " --$info '" . $database{$info} . "' ";
                next;
            }else
            {
                $error .= "[ERROR] LOST DATABASE $info,";
            }
        }

        # 数据库准备(物种)
        if($type eq "DATABASE_SPECIES" )
        {
            if(exists $database{"$Ref\_$info"})
            {
                $command .= " --$info '" . $database{"$Ref\_$info"} . "' ";
                next;
            }else
            {
                $error .= "[ERROR] LOST DATABASE_SPECIES $Ref\_$info,";
            }
        }

        # 数据库准备(物种，非必要)
        if($type eq "DATABASE_SPECIES_NOT_ENSSENTIAL" )
        {
            if(exists $database{"$Ref\_$info"})
            {
                $command .= " --$info '" . $database{"$Ref\_$info"} . "' ";
                next;
            }else
            {
                $warning .= "[WARNING] LOST DATABASE_SPECIES_NOT_ENSSENTIAL $Ref\_$info,";
            }
        }

        # 软件准备
        if($type eq "SOFTWARE")
        {
            if(exists $software{$info}{"cmd"})
            {
                $command .= " --$info '" . $software{$info}{"cmd"} . "' ";
                next;
            }else
            {
                $error .= "[ERROR] LOST SOFTWARE $info,";
            }
        }
    }
    push @commands, [$name, $command, $log, $error, $warning];
 
}
