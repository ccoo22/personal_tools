# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Parallel::ForkManager;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my %hashType = get_index_type();
my $db_dir   = SCRIPTDIR . "/bcl2fastq_ref";
my $example_dir = SCRIPTDIR . "/bcl2fastq_ref"; # 示例目录

# 检测 -> 脚本输入
my ($index_file, $output_dir, $index_type_list, $index_type, $if_help);
GetOptions(
    "index_file|i=s"        => \$index_file,
    "output_dir|o=s"        => \$output_dir,
    "index_type|t=s"        => \$index_type,
    "index_type_list|l:s"   => \$index_type_list,
    "help|h"                => \$if_help,
);
show_support_index(\%hashType) if(defined $index_type_list);
die help() if(defined $if_help or (not defined $index_type));
die "[Error] do not support index type: $index_type\n" if(not exists $hashType{$index_type});

$index_file = "$db_dir/input_$index_type.txt" if(not defined $index_file); # 如果没有指定预处理好的index文件，则默认使用数据库下的。作用：提示作用
$output_dir = PWD if(not defined $output_dir); # 如果没有指定输出目录，则默认是当前目录

my $output_file = "$output_dir/SampleSheet.csv"; # 输出文件
my $output_sh   = "$output_dir/work.sh"; # 数据拆分用脚本

print "Input:        $index_file\n";
print "Output:       $output_file\n";
print "Output shell: $output_sh\n";
print "Type:         $index_type\n";
###################################################################### 主程序

# sheet生成
open OUTPUT, ">$output_file";
print OUTPUT $hashType{$index_type}; # sheet表头

open INPUT, $index_file;
<INPUT>; # 第一行表头不看
while(<INPUT>)
{
    $_=~s/[\r\n]//g;
    next if($_!~/\w/);
    # hiseq数据，I7 index数据
    if($index_type eq 'hiseq')
    {
        my ($lane, $sample_id, $sample_name, $index_id_i7, $index_i7, $project, $tmp) = split /\t/, $_, 7;
        print OUTPUT "$lane,$sample_id,$sample_name,,,$index_id_i7,$index_i7,,,$project,\n";        
    }

    # hiseq数据，双index数据（I7、I5）
    if($index_type eq 'hiseq2index')
    {
        my ($lane, $sample_id, $sample_name, $index_id_i7, $index_i7, $index_id_i5, $index_i5, $project, $tmp) = split /\t/, $_, 9;
        print OUTPUT "$lane,$sample_id,$sample_name,,,$index_id_i7,$index_i7,$index_id_i5,$index_i5,$project,\n";    
    }

    # miseq数据，双index数据（I7、I5）。不需要Lane信息
    if($index_type eq 'hiseq250')
    {
        my ($sample_id, $sample_name, $index_id_i7, $index_i7, $index_id_i5, $index_i5, $project, $tmp) = split /\t/, $_, 8;
        print OUTPUT "$sample_id,$sample_name,,,$index_id_i7,$index_i7,$index_id_i5,$index_i5,$project,\n"; 
    }
}
close INPUT;
close OUTPUT;

# 生成数据拆分用脚本.注意：-R路径根据实际数据修改
my $bcl2fastq_cmd = "bcl2fastq -R /home/lch/media/disk2/190330_E00517_0431_BHW7TWCCXY/ -o Raw_data --sample-sheet SampleSheet.csv --barcode-mismatches 0 --no-lane-splitting -r 20  -p 20 -w 20 --minimum-trimmed-read-length 15 --mask-short-adapter-reads 15 > run.log 2>&1 &";
system("echo '$bcl2fastq_cmd' > $output_sh");
system("chmod 755 $output_sh");

###################################################################### 子函数

# 显示支持类型
sub show_support_index{
    my $hashType = shift @_;
    print "Support index type:\n";
    foreach my $type(sort keys %hashType)
    {
        print "    $type\n";
    }
    print "\nExample File dir: $example_dir\n\n";
    exit;
}

# 获取支持的index类型
sub get_index_type{
    my %hashType;
    my $id = `id`; 
    my ($user) = $id =~ /uid=\d+\((\w+)\)/;
    my $time   = `date +"%Y-%m-%d"`;
       $time   =~ s/[\r\n]//g;

    # sheet 头部配置参考
    # https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf
    # https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/sequencing-sheet-format-specifications-technical-note-970-2017-004.pdf
############
# hiseq头
############
    $hashType{'hiseq'} = <<"EOF";
[Header]
Date,$time
Investigator Name,$user
IEMFileVersion,4
Experiment Name,$user
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
151
151

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA+AGATCGGAAGAGCACACGTCTGAACTCCAGTCA+CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT+GATCGTCGGACTGTAGAACTCTGAACGTGT+CTGTCTCTTATACACATCTGACGCTGCCGACGA

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
EOF

############
# highseq2index头
############
    $hashType{'hiseq2index'} = <<"EOF";
[Header]
Date,$time
Investigator Name,$user
IEMFileVersion,4
Experiment Name,$user
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon


[Reads]
151
151

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA+AGATCGGAAGAGCACACGTCTGAACTCCAGTCA+CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT+GATCGTCGGACTGTAGAACTCTGAACGTGT+CTGTCTCTTATACACATCTGACGCTGCCGACGA

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
EOF

############
# highseq头 2X250
############
    $hashType{'hiseq250'} = <<"EOF";
[Header]
Date,$time
Investigator Name,$user
IEMFileVersion,4
Experiment Name,$user
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon


[Reads]
251
251

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA+CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT+CTGTCTCTTATACACATCTGACGCTGCCGACGA

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
EOF

    return %hashType;
# 注：bcl2fastq 主要识别的列信息有6列，Lane、Sample_ID、Sample_Name、Sample_Project、index、index2
# Lane When specified, the software generates FASTQ files only for the samples with the specified lane number.
# Sample_ID The sample ID.
# Sample_Name The sample name.
# Sample_Project The sample project name. The software creates a directory with the specified sample project name and puts the FASTQ files there. You can assign multiple samples to the same project.
# index The Index 1 (i7) index adapter sequence.
# index2 The Index 2 (i5) Index adapter sequence.

}

sub help{
    my $info = "
Program: bcl2fastq sheet generate， 生成下机数据拆分时使用的sheet文件
Version: 2019-05-30
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --index_file/-i            输入预处理好的index文件。注：可以为空，此时会根据index_type，从bcf2fastq_db目录下寻找输入文件
         --output_dir/-o            输出目录,默认当前目录
         --index_type/-t      [必填]index类型，例如：hiseq
         --index_type_list/-l       列出所有支持的index
         --help/-h                  查看帮助文档
    \n";
    return $info;
}


# 甘斌整理
# 接头列表

# truseq
# R1   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# R2   AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# NEB 
# R1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# R2 GATCGTCGGACTGTAGAACTCTGAACGTGT

# Nextra
# R1  CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
# R2  CTGTCTCTTATACACATCTGACGCTGCCGACGA

# Small RNA 接头
# R1  TGGAATTCTCGGGTGCCAAGG
# R2  GATCGTCGGACTGTAGAACTCTGAAC