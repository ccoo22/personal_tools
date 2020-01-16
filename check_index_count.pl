# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量

# 检测 -> 脚本输入
my ($index, $format, $if_help);
GetOptions(
    "index|i=s"   => \$index,
    "format|f=s"  => \$format,
    "help|h"      => \$if_help,
);
die help() if(defined $if_help or (not defined $index or not defined $format));
###################################################################### 主程序

open INDEX, $index;
my %hashSample;
my $lane = "";
while(<INDEX>)
{
    $_=~s/[\r\n]//g;
    my @infos = split /\t/, $_;
    if($infos[0] ne "")
    {
        $lane = $infos[0];
        next;
    }
    die "no lane code, please check input format\n" if($lane eq "");

    my $sample = $infos[1];
    my $project = ($format eq 'S') ? $infos[4] : $infos[6];

    die " project code error, please check input format\n" if(not defined $project or $project eq "");
    $hashSample{$project}{$lane}++;
}
close INDEX;

foreach my $project(sort keys %hashSample)
{
    foreach my $lane(sort keys %{$hashSample{$project}})
    {
        print "$project\t$lane\t$hashSample{$project}{$lane}\n";
    }
}





###################################################################### 子函数

sub get_seq_name{
    my $list_file = shift @_;
    my %hashSeqName;

    open LIST, $list_file;
    while(<LIST>)
    {
        next if($_!~/\w/);
        $_=~s/[\s\r\n]//g;
        $hashSeqName{$_} = 0;
    }
    close LIST;
    return %hashSeqName;
}

sub help{
    my $info = "
Program: check index count， 统计index文件每个项目的样本数量
Version: 2019-05-05
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --index/-i    index file . 注：输入文件第一行为表头，第一列表头为上机编号
         --format/-f   输入数据格式，S/D, S 表示single 单端index。Double 双端index
         --help/-h     查看帮助文档
    \n";
    return $info;
}

