# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Spreadsheet::XLSX;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# # 定义 -> 核心变量
# my $blastn      = "/home/ganb/soft/ncbi-blast-2.7.1+/bin/blastn";
# my $makeblastdb = "/home/ganb/soft/ncbi-blast-2.7.1+/bin/makeblastdb";
# my $nt_db       = "/home/pub/database/NCBI/NT/latest/nt";

# 检测 -> 脚本输入
my ($input_excel, $output_dir, $sheet_list, $if_help);
GetOptions(
    "input_excel|i=s"   => \$input_excel,
    "output_dir|o=s"    => \$output_dir,
    "sheet_list|s=s"    => \$sheet_list,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $input_excel or not defined $output_dir));
###################################################################### 主程序
 
# 需要的表格
my @need_sheets = (not defined $sheet_list) ? () : split /,/, $sheet_list;

# 开始拆解
my $workbook  = Spreadsheet::XLSX -> new ($input_excel); # 创建解析对象

foreach my $sheet (@{$workbook->{'Worksheet'}}) # 获取所有sheet对象
{   
    # 跳过不需要的表格
    next if(defined $sheet_list and not $sheet->{'Name'} ~~ @need_sheets); 
    # 输出文件
    my $txt_file = "$output_dir/$sheet->{'Name'}.txt";
    print "parse sheet : $sheet->{'Name'}, MaxRow=$sheet->{'MaxRow'}, MaxCol=$sheet->{'MaxCol'}\n";
    print "   txt file : $txt_file\n";

    open TXT, ">$txt_file";
    foreach my $row(0..$sheet->{'MaxRow'})
    {
        my @values;
        foreach my $col(0..$sheet->{'MaxCol'})
        {
            my $cell = $sheet->{'Cells'}[$row][$col]; # 获取对应行列的数值
            my $value = (defined $cell) ?  $cell -> {Val} : "";
            push @values, $value;
        }
        print TXT (join "\t", @values) . "\n";
    }
    close TXT;
}

# sheet 包含的键
# DefColWidth     默认列宽
# MinCol  sheet最小列数
# MaxCol  sheet最大列数
# MaxRow  sheet最大行数
# MinRow  sheet最小行数
# path    表格路径
# Name    sheet名称
# Cells   二元矩阵对象，数据存储在这里
###################################################################### 子函数

sub help{
    my $info = "
Program: parse excel to txt， 把excel表格转成txt文件
Version: 2019-09-20
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [必填]
         --input_excel/-i    输入excel文件
         --output_dir/-o     结果输出目录，文件以sheet名称命名
         
        [选填]
         --sheet_list/s      需要转换的sheet名称，多个名称用“,”逗号分隔，默认转换所有sheet
         --help/-h           查看帮助文档
    \n";
    return $info;
}

