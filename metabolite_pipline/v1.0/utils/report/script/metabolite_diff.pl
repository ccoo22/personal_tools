$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
代谢物 差异分析 结果整理
Version: v1.0 2020-7-1
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use Cwd 'abs_path'; 
use File::Spec;
use Getopt::Long;
use Excel::Writer::XLSX;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# Perl 自定义包
use lib SCRIPT_DIR. "../../../package";
use system_time qw(time_to_datetime);
use excel qw(set_excel_format txt_to_excel);
use utils qw(check_result_file parse_group);

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input, $output, $group_list, $sample_count_list, $if_help);
GetOptions(
	"input|i=s"               => \$input,
	"output|o=s"              => \$output,
	"group_list|g=s"          => \$group_list,
	"sample_count_list|c=s"   => \$sample_count_list,

	"help|h"    => \$if_help,
);
die "
Options: 必填

        --input/-i               代谢物diff分析输入路径， 例如： backup/metabolite/Group-A1_vs_A2_vs_A3/diff
        --output/-o              结果输出目录
        --group_list/-g          当前compare的分组名称列表。多个分组用逗号分隔， 例如： A,B,C  。 注意：务必保证顺序与config文件一致。 
        --sample_count_list/-c   与group_list一一对应的每一组包含的样本数量， 多个分组用逗号分隔， 例如: 1,2,3

    选填:
        --help/-h                查看帮助文档
\n" if (defined $if_help or not defined $input or not defined $output or not defined $group_list or not defined $sample_count_list);
system "mkdir -p $output";
my $DATA_TIME = time_to_datetime(time);
print "
---
Command: perl ".abs_path($0)." $ARGV_INFO
---
Start time: $DATA_TIME
";

$input       = abs_path($input);
$output      = abs_path($output);

my @group_names = split /,/, $group_list;
my @sample_counts = split /,/, $sample_count_list;


###################################################################### 初始化
my %hashGroup = parse_group(\@group_names, \@sample_counts);  # 分组文件拆分: (1) 所有样本 （2）两两分组的

my $is_ok = 1;
my $is_ok_compare = 0;  # 当前目录下是否有结果
foreach my $sub_compare_name(sort keys %hashGroup)
{
    my $output_sub_compare_dir = "$output/$sub_compare_name";
    system("mkdir -p $output_sub_compare_dir");

    my $diff_txt      = "$input/$sub_compare_name/diff.txt";
    my $ttest_volcano = "$input/$sub_compare_name/ttest.volcano.pdf";
    my $check = check_result_file('txt', $diff_txt); # 检查文件是否异常
    my $is_ok_sub    = 0;  # 当前子分组是否有结果
    if($check eq "")
    {
        my $excel = "$output_sub_compare_dir/diff.xlsx";
        my $workbook = Excel::Writer::XLSX->new($excel);
        my %format   = set_excel_format($workbook);

        txt_to_excel($workbook, \%format, $diff_txt, 'Diff');
        $workbook->close();

        system(" cp $ttest_volcano $output_sub_compare_dir/") if(-e $ttest_volcano);

        print "[Good] DIFF is OK [$group_list][$sub_compare_name]\n";
        $is_ok_sub++;
        $is_ok_compare++;
    }else
    {
        print "[Error] DIFF error [$group_list][$sub_compare_name]\n";
        $is_ok = 0;
    }

    # 当前子分组没有结果
    system("rm -r $output_sub_compare_dir") if($is_ok_sub == 0);
}

# 当前目录下没有任何结果
system("rm -r $output") if($is_ok_compare == 0);

if($is_ok == 1)
{
    print "[OK] 运行完成\n";
}
else
{
    print "[Error] DIFF 有问题，请仔细核查 \n";
}

###################################################################### 主程序

