# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 检测 -> 脚本输入
my ($input_dir, $output_file, $surffix_r1, $surffix_r2, $remove, $if_help);
GetOptions(
    "input_dir|i=s"   => \$input_dir,
    "output_file|o=s" => \$output_file,
    "surffix_r1=s"    => \$surffix_r1,
    "surffix_r2=s"    => \$surffix_r2,
    "remove=s"        => \$remove,
    "help|h"          => \$if_help,
);
die help() if(defined $if_help or (not defined $input_dir or not defined $surffix_r1 or not defined $surffix_r2));
###################################################################### 主程序
$output_file = "./fastq.detail.txt" if(not defined $output_file);
$remove      = "" if(not defined $remove);

# 读入样本名称
opendir SAMPLE, "$input_dir/";
my @samples = sort readdir SAMPLE;
closedir SAMPLE;

# 样本循环处理
open OUTPUT, ">$output_file";
print "output : $output_file\n";
print OUTPUT "Sample\tCondition\tR1_files\tR1_file_size\tR1_file_count\tR2_files\tR2_file_size\tR2_file_count\tSUM_size(GB)\n";
foreach my $sample(@samples)
{
    my $sample_dir = "$input_dir/$sample/";
    next if($sample !~ /\w/);
    next if(is_dir_ok($sample_dir) == 0); # 不是目录

    # 获取所有文件
    opendir FILE, "$sample_dir/";
    my @files = sort readdir FILE;
    closedir FILE;

    # R1/R2文件列表、大小，单位:GB
    my ($r1_file, $r1_size, $r1_count) = ('', '', 0);
    my ($r2_file, $r2_size, $r2_count) = ('', '', 0);
    my $sum_size = 0;
    foreach my $file_name(@files)
    {   
        next if($file_name !~ /\w/);
        next if($remove ne '' and $file_name =~ /$remove/); # 去掉指定字符文件

        my $file = "$sample_dir/$file_name";
        my $size = get_file_size($file, 'GB');

        # R1
        if($file_name =~ /(.*)($surffix_r1)$/)
        {   
            $r1_file .= "$file,";
            $r1_size .= "$size,";
            $r1_count++;
            $sum_size += $size;
        }
        # R2
        if($file_name =~ /(.*)($surffix_r2)$/)
        {
            $r2_file .= "$file,";
            $r2_size .= "$size,";
            $r2_count++;
            $sum_size += $size;
        }
    }
    # 样本状态检查
    my $condition = 'GOOD'; 
    $condition = 'R2 lost file' if($r1_count > $r2_count);
    $condition = 'R1 lost file' if($r1_count < $r2_count);
    $condition = 'lost file'    if($r1_count == $r2_count and $r1_count == 0);

    # 检查R1/R2样本名是否配对
    if($r1_count == $r2_count and $r1_count > 0)
    {
        my @r1_files = split /,/, $r1_file;
        my @r2_files = split /,/, $r2_file;
        foreach my $col(0..$#r1_files)
        {
            my $file1 = $r1_files[$col];
            my $file2 = $r2_files[$col];
            $file1 =~ s/($surffix_r1)$//;
            $file2 =~ s/($surffix_r2)$//;
            $condition = 'Not match' if($file1 ne $file2);
        }
    }

    # 输出
    print OUTPUT "$sample\t$condition\t$r1_file\t$r1_size\t$r1_count\t$r2_file\t$r2_size\t$r2_count\t$sum_size\n";    
}
close OUTPUT;




###################################################################### 子函数

# 获取文件大小
sub get_file_size{
    my $file  = shift @_;
    my $level = shift @_;

    my $size = -s $file;

    $size = $size / 1024       if($level eq 'KB');
    $size = $size / 1048576    if($level eq 'MB');
    $size = $size / 1073741824 if($level eq 'GB');

    $size = sprintf "%0.4f", $size;
    return $size;
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}

sub help{
    my $info = "
Program: 获取目录下，每个样本每条lane的文件信息（注：一个样本一个目录，且目录下有多条lanefastq）
Version: 2019-08-05
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [必填]
         --input_dir/-i      数据路径
         --surffix_r1        R1文件后缀，例如：_R1.fastq.gz,正则表达式匹配该字符，_1.fq.gz|_R1.fastq.gz 这样写也是可以的
         --surffix_r2        R2文件后缀，例如：_R2.fastq.gz
        
        [选填]
        --remove            要去掉包含该字符的文件，例如：Undetermined
        --output_file/-o    结果输出文件，默认：./fastq.detail.txt
        --help/-h           查看帮助文档
    \n";
    return $info;
}

