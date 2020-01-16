# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
# my $gatk4      = "/home/genesky/software/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar";
my $gatk4      = "/home/genesky/software/gatk/4.1.2.0/gatk-package-4.1.2.0-local.jar";
my $java       = "/home/pub/software/jdk1.8.0_172/bin/java";
my $genome     = "/home/pub/database/Human/hg19/genome/hg19.fa";
my $tmp_dir    = "/home/tmp";


# 检测 -> 脚本输入
my ($bam, $region, $output_dir, $prefix, $gatk4_para, $if_help);
GetOptions(
    "bam|i=s"           => \$bam,
    "region|r=s"        => \$region,
    "output_dir|o=s"    => \$output_dir,
    "prefix|p=s"        => \$prefix,
    "gatk4_para|a=s"    => \$gatk4_para,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $bam or not defined $output_dir or not defined $prefix));
$gatk4_para = "" if(not defined $gatk4_para);
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
my $command      = "$java -jar -Xmx6g $gatk4 HaplotypeCaller -I $bam -O $gVCF -R $genome -ERC GVCF --min-base-quality-score 20  --sample-ploidy 2 $region_limit -bamout $bam_out $gatk4_para  --tmp-dir $tmp_dir ";
system("$command");
system("echo '$command' > $command_file");
chmod 0755, $command_file;
 

print "program run parameter:\n";
print " bam        $bam\n";
print " gVCF       $gVCF\n";
print " Genome     $genome\n";
print " Region     $region\n";
print " bam_out    $bam_out\n";
print " gatk4_para $gatk4_para\n";

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
         --gatk4_para/-a  其他GATK4参数，注意添加单引号，示例：pcr参数 '--max-reads-per-alignment-start 0 --dont-use-soft-clipped-bases true --kmer-size 10 --kmer-size 15 --kmer-size 20 --kmer-size 25 --kmer-size 30 --kmer-size 35 --kmer-size 40 --max-num-haplotypes-in-population 256 '
                          其他参数：--max-assembly-region-size 200  --bam-writer-type ALL_POSSIBLE_HAPLOTYPES -L test.bed
         --help/-h        查看帮助文档
    \n";
    return $info;
}

