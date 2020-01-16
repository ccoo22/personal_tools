# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $samtools      = "/home/genesky/software/samtools/1.9/samtools";
 
# 检测 -> 脚本输入
my ($input_bam, $output_bam, $mt_name_in_bam, $rename_mt_in_bam, $if_help);
GetOptions(
    "input_bam|i=s"      => \$input_bam,
    "output_bam|o=s"     => \$output_bam,
    "mt_name_in_bam=s"   => \$mt_name_in_bam,
    "rename_mt_in_bam=s" => \$rename_mt_in_bam,
    "help|h"             => \$if_help,
);
die help() if(defined $if_help or (not defined $input_bam or not defined $output_bam or not defined $mt_name_in_bam));
###################################################################### 主程序

##########
# 1 获取bam里记录的线粒体长度,并生成bed
#########
print "获取MT长度\n";
my $mt_length = get_mt_length_from_bam($input_bam, $samtools, $mt_name_in_bam);
my $mt_bed    = "$output_bam.mt.bed";
system("echo $mt_name_in_bam    1   $mt_length  +   MT > $mt_bed ");

##########
# 2 bam 转 sam，并处理比对信息
#########
print "MT bam 文件生成\n";
my $mt_sam = "$output_bam.sam";
open BAM, "$samtools view -h $input_bam -L $mt_bed|";
open SAM, "|$samtools view -bS - > $output_bam";
while(<BAM>)
{
    my @datas = split /\t/, $_;
    my $reads_name = $datas[0];
    next if($reads_name eq '@SQ' and $datas[1] ne "SN:$mt_name_in_bam");  # 去掉非MT的染色体信息

    # 如果需要修改染色体名称，需要处理一下
    if(defined $rename_mt_in_bam)
    {   
        # 染色体名称替换
        if($reads_name eq '@SQ' and $datas[1] eq "SN:$mt_name_in_bam")
        {
            $datas[1] = "SN:$rename_mt_in_bam";
        }

        # reads内部信息替换
        if($reads_name !~ /^@/)
        {   
            # 染色体名字替换
            $datas[2] = $rename_mt_in_bam;

            # mate序列比对信息调整
            if($datas[6] ne '=')
            {
                $datas[6] = '=';
                $datas[7] = $datas[3];
                $datas[8] = 0;
            }
        }
    }
    print SAM (join "\t", @datas);
}
close BAM;
close SAM;

# 建索引
system("samtools index -b $output_bam");


###################################################################### 子函数

sub get_mt_length_from_bam{
    my $input_bam      = shift @_;
    my $samtools       = shift @_;
    my $mt_name_in_bam = shift @_;

    # 获取bam头部
    my $mt_head = `$samtools view -H $input_bam|grep SN:$mt_name_in_bam`;

    my ($mt_info, $tmp) = ($mt_head =~ /LN:\d+/g);
    my ($mt_length) = $mt_info =~ /(\d+)/;
    die "[error] 表头里发现多个 SN:$mt_name_in_bam 字符，请仔细核查是否写错" if(defined $tmp);
    
    return $mt_length;
}

sub help{
    my $info = "
Program: 从bam文件里提取MT的比对结果，生成只包含MT信息的bam文件。
Version: 2019-12-26
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [必填]
         --input_bam/-i   输入bam文件，注：必须经过排序，有索引
         --output_bam/-o  输出bam文件
         --mt_name_in_bam bam文件中，线粒体的字符名名称，例如：MT/M

        [选填]
         --rename_mt_in_bam 把新bam里的线粒体名称修改为指定名称，例如：MT/M
         --help/-h        查看帮助文档
    \n";
    return $info;
}

