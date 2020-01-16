# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Cwd qw( abs_path );

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"}; 

# 定义 -> 核心变量
my $tools_dir                = SCRIPTDIR . "/tools";
my $kallisto_to_reads_count_script = "$tools_dir/kallisto_to_reads_count.r";                  # 格式转换工具

# 软件、环境设置
my $Rscript         = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $Rlib            = "/home/genesky/software/r/3.5.1/lib64/R/library";
$ENV{"R_LIBS"} = $Rlib; # R包路径
# 本流程需要的R包汇总
# library(tximport)
# library(EnsDb.Hsapiens.v86)
# library(ensembldb)
 
# 检测 -> 脚本输入
my ($file_list, $species, $output, $if_help);
GetOptions(
    "file_list|l=s"   => \$file_list,
    "species|s=s"     => \$species,
    "output|o=s"      => \$output,
    "help|h"          => \$if_help,
);
die help() if(defined $if_help or (not defined $file_list or not defined $output));

###################################################################### 主程序

system("$Rscript $kallisto_to_reads_count_script $file_list $output ");

###################################################################### 子函数



sub help{
    my $info = "
Program: kallisto_abundance_to_reads_count， kallisto定量结果格式转换（用于DESEQ2分析）
Version: 2019-05-16
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --file_list/-l    [必填] kallisto定量输出结果的文件列表,目前仅支持人。两列，第一列：样本名，第二列：文件路径，没有表头
         --output/-o       [必填] 输出文件
         --help/-h         查看帮助文档
         注意：GO分析要慢一些
    \n";
    return $info;
}
