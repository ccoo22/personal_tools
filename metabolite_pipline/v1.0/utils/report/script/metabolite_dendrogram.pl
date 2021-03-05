$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
代谢物 dendrogram 结果整理
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
use excel qw(set_excel_format);
use utils qw(check_result_file parse_group);

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input, $output, $group_list, $sample_count_list, $compound_count, $diff_model, $if_help);
GetOptions(
	"input|i=s"               => \$input,
	"output|o=s"              => \$output,
	"group_list|g=s"          => \$group_list,
	"sample_count_list|c=s"   => \$sample_count_list,
	"compound_count=s"        => \$compound_count,
	"diff_model"              => \$diff_model,

	"help|h"    => \$if_help,
);
die "
Options: 必填

        --input/-i               代谢物dendrogram分析输入路径， 例如： backup/metabolite/Group-A1_vs_A2_vs_A3/dendrogram
        --output/-o              结果输出目录
        --group_list/-g          当前compare的分组名称列表。多个分组用逗号分隔， 例如： A,B,C  。 注意：务必保证顺序与config文件一致。 
        --sample_count_list/-c   与group_list一一对应的每一组包含的样本数量， 多个分组用逗号分隔， 例如: 1,2,3
        --compound_count         代谢物的数量

    选填:
        --diff_model             差异代谢物模式。 当声明了这个参数，则处理diff_compound结果。
        --help/-h                查看帮助文档
\n" if (defined $if_help or not defined $input or not defined $output or not defined $group_list or not defined $sample_count_list or not defined $compound_count);
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

$diff_model = (defined $diff_model) ? "diff_model" : "";

###################################################################### 初始化
my %hashGroup = parse_group(\@group_names, \@sample_counts);  # 分组文件拆分: (1) 所有样本 （2）两两分组的

my $is_ok = 1;
my $is_ok_compare = 0;  # 当前目录下是否有结果
foreach my $sub_compare_name(sort keys %hashGroup)
{
    my $output_sub_compare_dir = "$output/$sub_compare_name";
    system("mkdir -p $output_sub_compare_dir");

    my $sample_count = $hashGroup{$sub_compare_name}{'SampleCount'};
    my $is_ok_sub    = 0;  # 当前子分组是否有结果
    if($sample_count > 1 and $compound_count > 1)
    {   
        my @dendrogram_results;
        # 距离计算方法
        foreach my $distance_method(qw(euclidean pearson spearman))  
        {   # 聚类计算方法
            foreach my $cluster_method(qw(ward.D single complete average))  
            {   
                my $dendrogram = "$input/$sub_compare_name/dendrogram_$distance_method\_$cluster_method.pdf";
                push @dendrogram_results, $dendrogram;
            }
        }
        my $check = check_result_file('pdf', @dendrogram_results);

        if($check eq "")
        {
            system("cp @dendrogram_results $output_sub_compare_dir/");
            print "[Good] dendrogram $diff_model is OK [$group_list][$sub_compare_name]\n";
            $is_ok_sub++;
            $is_ok_compare++;
        }else
        {
            print "[Error] dendrogram $diff_model error [$group_list][$sub_compare_name] $check\n";
            $is_ok = 0;
        }
    }else
    {
        print "[Skip] in dendrogram $diff_model, sample count is not over 1 or compound_count is not over 1 in [$group_list][$sub_compare_name]\n";
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
    print "[Error] dendrogram $diff_model 有问题，请仔细核查 \n";
}

###################################################################### 主程序
