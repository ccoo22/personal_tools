$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
代谢物丰度表格整理
Version: v1.0 2020-7-1
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
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
use utils qw(check_result_file is_file_ok);

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input, $output, $if_help);
GetOptions(
	"input|i=s"    => \$input,
	"output|o=s"   => \$output,

	"help|h"    => \$if_help,
);
die "
Options: 必填

        --input/-i             代谢物分析结果输入目录，例如： backup/metabolite. 该目录下至少包含：metabolite_data.normalized.txt、metabolite_name_map.format.txt
        --output/-o            结果输出目录
        --help/-h              查看帮助文档
\n" if (defined $if_help or not defined $input or not defined $output);
$input      = File::Spec->rel2abs($input);
$output      = File::Spec->rel2abs($output);
system "mkdir -p $output";

my $DATA_TIME = time_to_datetime(time);
print "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
";


###################################################################### 初始化
my $raw_file          = "$input/metabolite_data.raw.txt";
my $replace_miss_file = "$input/metabolite_data.replace.missing.txt";
my $normalized_file   = "$input/metabolite_data.normalized.txt";
my $map_file          = "$input/metabolite_name_map.format.txt";

die "[Error] 归一化 文件缺失, $normalized_file\n" if(is_file_ok($normalized_file) == 0);
die "[Error] 代谢物查找 文件缺失, $map_file\n"    if(is_file_ok($map_file) == 0);

my $excel = "$output/aboundance.xlsx";
my $workbook = Excel::Writer::XLSX->new($excel);
my %format   = set_excel_format($workbook);

txt_to_excel($workbook, \%format, $raw_file, 'Raw')                  if(is_file_ok($raw_file));
txt_to_excel($workbook, \%format, $replace_miss_file, 'ReplaceMiss') if(is_file_ok($replace_miss_file));
txt_to_excel($workbook, \%format, $normalized_file, 'Normalized')    if(is_file_ok($normalized_file));
txt_to_excel($workbook, \%format, $map_file, 'Map')                  if(is_file_ok($map_file));
$workbook->close();


if(is_file_ok($excel) == 1)
{
    print "[OK] 运行完成\n";
}
else
{
    print "[Error] excel丰度表格没有生成 : $excel\n";
}

###################################################################### 主程序


