$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
代谢物分析结果整理
Version: v1.0 2020-7-1
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
my $SCRIPT_DIR = SCRIPT_DIR;

# Perl 自定义包
use lib SCRIPT_DIR."../../package";
use system_time qw(time_to_datetime);



# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($config_ini, $backup_dir, $report_dir, $force_cover, $if_help);
GetOptions(
	"config-ini|c=s"   => \$config_ini,
	"backup-dir|b=s"   => \$backup_dir,
	"report-dir|r=s"   => \$report_dir,

	"force-cover:s"    => \$force_cover,
	"help|h" => \$if_help,
);
die "
Options: 必填

        --config-ini/-c            流程配置文件， 只需要 Compare 定义信息
        --backup-dir/-b            Backup_Dir 路径， 项目运行的中间结果路径
        --report-dir/-r            Report_Dir 路径， 输出路径，提供给客户

Options: 可选
        --force-cover              删除所有原始结果文件, 完全覆盖
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $config_ini  or not defined $backup_dir or not defined $report_dir);

###################################################################### 初始化
# 运行参数
my $DATA_TIME = time_to_datetime(time);
my $RUN_INFO = "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
\n";

# 日志
$DATA_TIME =~ s/\D+/\_/g;
$DATA_TIME =~ s/\D+/\_/g;
my $log_dir = "$backup_dir/log";
system("mkdir -p $log_dir") if(not -e $log_dir);

my $log = "$log_dir/log.report.$DATA_TIME.txt";  # 记录所有分析结果的情况，用于检查是否异常
open SAVE, ">>$log"; print SAVE $SCRIPT_INFO.$RUN_INFO; close SAVE;
print $RUN_INFO;


my $sample_group = "$backup_dir/metabolite/sample_group.txt";
my %config = read_config($config_ini);
my %sample_group = read_sample_group($sample_group);  # 读入分组文件中每一列下每一组的样本信息（第3列开始是分组信息）

# 分析分组
my %group;
my $id = 0;
foreach my $compare_info(sort{$config{$a} <=> $config{$b}} keys %config){
	my ($column_name, $group_list) = split /:/, $compare_info;
	die "\n[ERROR] --config 文件 Compare 格式错误, $compare_info\n" if(not defined $group_list);
	die "\n[ERROR] sample_group.txt 分组列名没有定义, $column_name\n" if(not exists $sample_group{$column_name});
	my @group_names = split /,/, $group_list;
	$group{$compare_info}{"compare_name"} = "$column_name-".(join "_vs_", @group_names);
	$group{$compare_info}{"sort_order"} = $id;
    my $group_order = 0;
	foreach my $group_name(@group_names){
		die "\n[ERROR] sample_group.txt 分组列 $column_name 中 分组没有定义, $group_name\n" if(not exists $sample_group{$column_name}{$group_name});
		my @samples = split /,/, $sample_group{$column_name}{$group_name};
		my $sample_count = scalar(@samples);

		$group{$compare_info}{"group"}{$group_name}        = $group_order;
		$group{$compare_info}{"sample"}{$group_name}       = $sample_group{$column_name}{$group_name};
		$group{$compare_info}{"sample_count"}{$group_name} = $sample_count;
		$group{$compare_info}{"group_list"}         .= "$group_name,";
		$group{$compare_info}{"group_count"}++;
		$group{$compare_info}{"sample_count_list"}  .= "$sample_count,";

        $group_order++;
	}
	$id++;
}

###################################################################### 主程序
# 丰度表整理
my @command = (
        ["metabolite", "aboundance", "$SCRIPT_DIR/script/metabolite_aboundance.pl"],
);
run_command(@command);

# 代谢物数量
my $compound_count = `wc -l $backup_dir/metabolite/metabolite_data.normalized.txt`;
   ($compound_count) = split /\s+/, $compound_count;
   $compound_count -= 1;

# 每个compare整理
foreach my $compare_info (sort { $group{$a}{"sort_order"} <=> $group{$b}{"sort_order"} } keys %group)
{		
	system "echo '######################################################################\n报告分组 $compare_info' | tee -a  $log ";

	my $group_list        = $group{$compare_info}{"group_list"};
	my $group_count       = $group{$compare_info}{"group_count"};
	my $sample_count_list = $group{$compare_info}{"sample_count_list"};
	my $compare_name      = $group{$compare_info}{"compare_name"};

	# 分组信息
	system("mkdir -p $report_dir/$compare_name");
	system("cp $backup_dir/metabolite/$compare_name/sample_group.all_samples.txt $report_dir/$compare_name/sample.group.txt");
	system("sed -i '1i sample\\tgroup' $report_dir/$compare_name/sample.group.txt ");

	# 常规QC
	my @command = (
        ["metabolite/$compare_name/pca",        "$compare_name/01_Quality_Control/pca",        "$SCRIPT_DIR/script/metabolite_pca.pl        -g '$group_list' -c '$sample_count_list' --compound_count $compound_count "],
        ["metabolite/$compare_name/plsda",      "$compare_name/01_Quality_Control/plsda",      "$SCRIPT_DIR/script/metabolite_plsda.pl      -g '$group_list' -c '$sample_count_list' --compound_count $compound_count "],
        ["metabolite/$compare_name/dendrogram", "$compare_name/01_Quality_Control/dendrogram", "$SCRIPT_DIR/script/metabolite_dendrogram.pl -g '$group_list' -c '$sample_count_list' --compound_count $compound_count "],
        ["metabolite/$compare_name/heatmap",    "$compare_name/01_Quality_Control/heatmap",    "$SCRIPT_DIR/script/metabolite_heatmap.pl    -g '$group_list' -c '$sample_count_list' --compound_count $compound_count "],
	);
	# 差异分析
	if($group_count > 1)
	{
		push @command, ["metabolite/$compare_name/diff",       "$compare_name/02_Different_Analysis/diff",            "$SCRIPT_DIR/script/metabolite_diff.pl                                                                   -g '$group_list' -c '$sample_count_list' "];
        push @command, ["metabolite/$compare_name/pca",        "$compare_name/02_Different_Analysis/pca_diff",        "$SCRIPT_DIR/script/metabolite_pca_diff.pl        --diff_dir '$backup_dir/metabolite/$compare_name/diff' -g '$group_list' -c '$sample_count_list' "];
        push @command, ["metabolite/$compare_name/plsda",      "$compare_name/02_Different_Analysis/plsda_diff",      "$SCRIPT_DIR/script/metabolite_plsda_diff.pl      --diff_dir '$backup_dir/metabolite/$compare_name/diff' -g '$group_list' -c '$sample_count_list' "];
        push @command, ["metabolite/$compare_name/dendrogram", "$compare_name/02_Different_Analysis/dendrogram_diff", "$SCRIPT_DIR/script/metabolite_dendrogram_diff.pl --diff_dir '$backup_dir/metabolite/$compare_name/diff' -g '$group_list' -c '$sample_count_list' "];
        push @command, ["metabolite/$compare_name/heatmap",    "$compare_name/02_Different_Analysis/heatmap_diff",    "$SCRIPT_DIR/script/metabolite_heatmap_diff.pl    --diff_dir '$backup_dir/metabolite/$compare_name/diff' -g '$group_list' -c '$sample_count_list' "];
        
		push @command, ["metabolite/$compare_name/kegg",       "$compare_name/03_Pathway_Analysis",                   "$SCRIPT_DIR/script/metabolite_pathway.pl         --diff_dir '$backup_dir/metabolite/$compare_name/diff' -g '$group_list' -c '$sample_count_list' "];

	}else
	{
		system("echo '[RUNTIME]   : 代谢物 差异分析等相关内容跳过，因为只有一组数据，  ' | tee -a  $log"); # 运行效果输出到log文件以及屏幕
	}
	run_command(@command);
}

my $error = `grep '\\[Error\\]' $log`;
if($error =~/\w/)
{
    print "\n[ERROR] 运行错误, $error , 请核查文件 $log\n";
}else{
    print "\n[OK] 运行完成\n";
}


###################################################################### 子程序
sub run_command{
	my @command = @_;
	foreach my $step(0..@command-1){
		my ($i, $o, $soft) = map{$command[$step]->[$_]} 0..2;
		system "rm -rf $report_dir/$o" if(defined $force_cover);
		my $rtn     = `perl $soft -i $backup_dir/$i -o $report_dir/$o 2>&1`;
		my @split   = split /[\r\n]/, $rtn;
		my $subsoft = (exists $split[2]) ? $split[2] : "N/A";  # 子软件名称
		my $print   = join ",", grep /\[\w+\]/, @split;  # 提取自软件的警告、错误输出信息
		   $print   = $split[$#split] if($print !~ /\w/);  # 没有的话，提取最后一个信息
		
		system("echo '[RUNTIME]   : $subsoft ... $print' | tee -a  $log"); # 运行效果输出到log文件以及屏幕
	}
}


sub read_config{
	my $file = shift;
	my %hash;
	my $sort = 0;
	open FILE, $file;
	while(my $line = <FILE>){
		$line =~ s/#.+//g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		my ($a, $b) = split /\s*=\s*/, $line, 2;
		next if(not defined $b);
		next if($a ne "Compare");
		$b =~ s/\s//g;
		$hash{$b} = $sort;
		$sort++;
	}
	close FILE;
	return %hash;
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

sub read_sample_group{
	my $file = shift;
	open FILE, $file or die "[ERROR] 无法读取分组文件, $file\n";
	my $title = <FILE>;
	$title =~ s/[\r\n]//g;
	my @split_title = split /\t/, $title;
	my %hash;
	my $sort = 0;
	while(my $line = <FILE>){
		$line =~ s/[\r\n]//g;
		my @split_line = split /\t/, $line;
		for my $col(2..@split_line-1){
			$hash{$split_title[$col]}{$split_line[$col]} .= "$split_line[1],";
		}
	}
	close FILE;
	return %hash;
}