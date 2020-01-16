# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
# my $gatk3      = "/home/pub/software/GATK3.8/GenomeAnalysisTK.jar";
my $gatk3      = "/home/pub/software/GATK3/GenomeAnalysisTK.jar"; # 3.5
my $java       = "/home/pub/software/jdk1.8.0_172/bin/java";
my $genome     = "/home/pub/database/Human/hg19/genome/hg19.fa";
my $tmp_dir    = "/home/tmp";


# 检测 -> 脚本输入
my ($bam, $region, $output_dir, $prefix, $gatk3_para, $if_help);
GetOptions(
    "bam|i=s"           => \$bam,
    "region|r=s"        => \$region,
    "output_dir|o=s"    => \$output_dir,
    "prefix|p=s"        => \$prefix,
    "gatk3_para|a=s"    => \$gatk3_para,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $bam or not defined $output_dir or not defined $prefix));
$gatk3_para = "" if(not defined $gatk3_para);
$region     = "" if(not defined $region);
###################################################################### 主程序

# (1) 确定区域
my $interval = "$output_dir/$prefix.interval.bed";
if($region ne '')
{
    open BED, ">$interval";
    foreach my $region1(split /,/, $region)
    {
        my ($chr, $start, $end) = split /[:-]/, $region1;
        die "Region error:$region1\tsupport format = chr:start-end\n" if(not defined $end);
        print BED "$chr\t$start\t$end\n";
    }
    close BED;   
}

# (2) GVCF分析
my $region_limit = ($region eq "") ? "" : " -L $interval"; # 区域限制
my $gVCF         = "$output_dir/$prefix.g.vcf.gz";
my $bam_out      = "$output_dir/$prefix.bamout.bam";
my $command_file = "$output_dir/$prefix.command.line.sh";
my $command      = "$java -Djava.io.tmpdir=$tmp_dir -jar -Xmx5g $gatk3 -T HaplotypeCaller -I $bam -R $genome --min_base_quality_score 20 --emitRefConfidence GVCF -o $gVCF $region_limit -bamout $bam_out  $gatk3_para --sample_ploidy 2 ";
system("$command");
system("echo '$command' > $command_file");
chmod 0755, $command_file;
 

print "program run parameter:\n";
print " bam        $bam\n";
print " gVCF       $gVCF\n";
print " Genome     $genome\n";
print " Region     $region\n";
print " bam_out    $bam_out\n";
print " gatk3_para $gatk3_para\n";

###################################################################### 子函数


sub help{
    my $info = "
Program: GATK4 GVCF， 人源， GVCF分析
Version: 2019-01-30
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --bam/-i         [必填] bam文件, 示例：/home/pub/output/Target/18B1025D/1557/1557_final.bam
         --region/-r      区域, 多个区域用‘,’分割, 示例：6:18139134-18139334,13:48611612-48611812
         --output_dir/-o  [必填] 输出目录, 结果输出目录, 示例：/home/ganb/work/tmp/18B1025D_check
         --prefix/-p      [必填] 输出文件前缀, 示例：1557_test
         --gatk3_para/-a  其他GATK3参数，注意添加单引号，示例：'--downsampling_type NONE --dontUseSoftClippedBases --kmerSize 30 --maxNumHaplotypesInPopulation 256  '
         --help/-h        查看帮助文档
    \n";
    return $info;
}

