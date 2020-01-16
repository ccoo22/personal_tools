# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 检测 -> 脚本输入
my ($dir, $output, $if_help);
GetOptions(
    "dir|d=s"     => \$dir,
    "output|o=s"  => \$output,
    "help|h"      => \$if_help,
);
die help() if(defined $if_help or (not defined $dir or not defined $output));
###################################################################### 主程序
 
# 目录检查
my @all_files;
check_dir($dir);

# 输出
open OUTPUT, ">$output";
print OUTPUT "File\tType\tSize\tDate\n";
print OUTPUT (join "\n", @all_files) . "\n";
close OUTPUT;






###################################################################### 子函数
sub check_dir{
    my $dir = shift @_;
    
    # 获取目录下所有文件
    opendir DIR, $dir;
    my @file_lists = readdir DIR;
    closedir DIR;

    # 检查每一个文件、目录
    foreach my $file_name(@file_lists)
    {   
        my $file = "$dir/$file_name";

        next if($file_name!~/\w/);

        if(-f $file)
        {   
            push @all_files, get_info($file, 'file');
        }else
        {   
            push @all_files, get_info($file, 'dir');
            check_dir($file);
        }
    }

}

# 提取文件、目录的信息
sub get_info{
    my $file = shift @_;
    my $type = shift @_;
    my @infos;

    push @infos, $file;
    push @infos, $type;
    push @infos, -s $file; #  文件大小

    my $mtime = (stat $file)[9];
    my @tmp = localtime $mtime; 
    my $date = sprintf "%02u-%02u-%02u %02u:%02u", $tmp[5] % 100,  $tmp[4] + 1 ,$tmp[3], $tmp[2], $tmp[1];

    push @infos, $date; #  文件日期

    return (join "\t", @infos);
}

 
sub help{
    my $info = "
Program: list_all_files， 列出目录下所有的文件信息
Version: 2019-10-18
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --dir/-d      要检测的目录
         --output/-o   输出文件
         --help/-h     查看帮助文档
    \n";
    return $info;
}

