# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Cwd 'abs_path';

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
#$input_dir    = abs_path($input_dir);

###################################################################### 主程序
$output_file = "./fastq.detail.txt" if(not defined $output_file);
$remove      = "" if(not defined $remove);

open OUTPUT, ">$output_file";
print "output : $output_file\n";
print OUTPUT "Sample\tcomdition\tR1_files\tR1_file_size\tR2_files\tR2_file_size\tSUM_size(GB)\n";
check_dir_fastq($input_dir);
close OUTPUT;


###################################################################### 子函数

sub check_dir_fastq{
    my $dir = shift @_;
    print "check: $dir \n";
    if(not -r $dir)
    {
        print "[error] 权限不足，无法读取\n";
        return;
    }

    opendir DIR, "$dir/";
    my %hashFastq;
    my @sub_dirs;
    foreach my $file_name(readdir DIR)
    {    
        next if($file_name eq '.' or $file_name eq '..');
        next if($remove ne '' and $file_name =~ /$remove/); # 去掉指定字符文件
        my $file = "$dir/$file_name";
        my $size = (-d $file) ? 0 : get_file_size($file, 'GB');
         
        # 次级目录
        if(-d $file)
        {
            push @sub_dirs, $file;
            next;
        }
         
        # R1
        if($file_name =~ /^(.*)($surffix_r1)$/)
        {   
            $hashFastq{$1}{'R1-file'} .= "$file,";
            $hashFastq{$1}{'R1-size'} += $size;
            $hashFastq{$1}{'file-count'}++;
        }
        # R2
        if($file_name =~ /^(.*)($surffix_r2)$/)
        {
            $hashFastq{$1}{'R2-file'} .= "$file,";
            $hashFastq{$1}{'R2-size'} = $size;
            $hashFastq{$1}{'file-count'}++;
        }
    }
    closedir DIR;

    # 输出当前目录下的fastq信息
    my $count = 0;
    foreach my $sample(sort keys %hashFastq)
    {   
        my ($r1, $r2) = ("", "");
        my ($r1_size, $r2_size, $size) = (0, 0, 0);
        my $file_count = $hashFastq{$sample}{'file-count'};

        if(exists $hashFastq{$sample}{'R1-file'})
        {
            $r1 = $hashFastq{$sample}{'R1-file'};
            $r1_size = $hashFastq{$sample}{'R1-size'};
        }
        if(exists $hashFastq{$sample}{'R2-file'})
        {
            $r2 = $hashFastq{$sample}{'R2-file'};
            $r2_size = $hashFastq{$sample}{'R2-size'};
        }
        $size = $r1_size + $r2_size;


        my $condition = 'error';
           $condition = 'paired' if($r1 ne '' and $r2 ne '' and $file_count == 2);
           $condition = 'empty'  if($size == 0);
        
        print OUTPUT "$sample\t$condition\t$r1\t$r1_size\t$r2\t$r2_size\t$size\n";
        $count++;
    }
    print "        [find] sample count $count\n";

    # 检查子目录
    foreach my $sub_dir(@sub_dirs)
    {   
        check_dir_fastq2($sub_dir);
    }
}


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
Program: 获取目录下，所有的样本的fastq数据路径
         脚本会遍历目标路径下的所有子目录
Version: 2021-02-09
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

