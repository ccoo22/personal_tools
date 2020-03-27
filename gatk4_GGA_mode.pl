# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $gatk4      = "/home/genesky/software/gatk/4.1.3.0/gatk-package-4.1.3.0-local.jar";
my $java       = "/home/genesky/software/java/1.8.0_231/bin/java";
my $bcftools   = "/home/genesky/software/bcftools/1.9/bin/bcftools";
# my $genome     = "/home/pub/database/Human/hg19/genome/hg19.fa";
my $tmp_dir    = "/home/tmp";


# 检测 -> 脚本输入
my ($bam, $vcf, $genome, $output_prefix, $gatk4_para, $if_help);
GetOptions(
    "bam|i=s"           => \$bam,
    "vcf|v=s"           => \$vcf,
    "genome|g=s"        => \$genome,
    "output_prefix|p=s" => \$output_prefix,
    "gatk4_para|a=s"    => \$gatk4_para,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $bam or not defined $vcf or not defined $genome or not defined $output_prefix));
$gatk4_para = "--max-reads-per-alignment-start 0 --dont-use-soft-clipped-bases true --kmer-size 10 --kmer-size 15 --kmer-size 20 --kmer-size 25 --kmer-size 30 --kmer-size 35 --kmer-size 40 --max-num-haplotypes-in-population 256" if(not defined $gatk4_para);
###################################################################### 主程序

my $log = "$output_prefix.log";
# (1) 检查参考vcf文件是否创建索引
system("$java -jar $gatk4 IndexFeatureFile -F $vcf 2>&1 |tee -a $log") if(not -e "$vcf.idx");

# (2) VCF 转bed (交给gatk4，否则gatk4会读取全部bam数据)
my $bed = "$output_prefix.bed";
vcf2bed($vcf, $bed, 50);

# (2) 分型
my $vcf_raw    = "$output_prefix.raw.vcf";
my $vcf_filter = "$output_prefix.filter.vcf";
my $vcf_final  = "$output_prefix.final.vcf";
system("$java -jar $gatk4 HaplotypeCaller -R $genome --alleles $vcf -L $bed --min-base-quality-score 20 --genotyping-mode GENOTYPE_GIVEN_ALLELES --output-mode EMIT_ALL_SITES -I $bam -O $vcf_raw $gatk4_para 2>&1 |tee -a $log");
system("$java -jar $gatk4 VariantFiltration -R $genome -V $vcf_raw --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'CallError' -O $vcf_filter 2>&1 |tee -a $log");
system("$bcftools norm -m -both -f $genome -O z -o $vcf_final $vcf_filter 2>&1 |tee -a $log");

print "vcf raw   :\n    $vcf_raw\n\n";
print "vcf filter [对分型位点打分] :\n    $vcf_filter\n\n";
print "vcf final [拆分多态位点] :\n    $vcf_final\n\n";

 

###################################################################### 子函数

sub vcf2bed{
    my $vcf = shift @_;
    my $bed = shift @_;
    my $extend = shift @_;

    open VCF, $vcf;
    open BED, ">$bed";
    while(<VCF>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/ or $_=~/^#/);

        my ($chr, $pos, $tmp) = split /\t/, $_;
        my $start = $pos - $extend;
        my $end   = $pos + $extend;
        $start = 0 if($start < 0);

        print BED "$chr\t$start\t$end\t+\tNONE\n";
    }
    close BED;
    close VCF;
}


sub help{
    my $info = "
Program: GATK4 genotype vcf， 对指定的vcf位点进行分型
Version: 2020-03-19
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [必填]
         --bam/-i            bam文件, 示例：/home/pub/output/Target/18B1025D/1557/1557_final.bam
         --vcf/-v           [必填] 要分型的vcf文件
         --genome/-g        [必填] 参考基因组文件 例如： /home/pub/database/Human/hg19/genome/hg19.fa
         --output_prefix/-p [必填] 输出文件， 示例 ： /output/sample

        [选填]
         --gatk4_para/-a     其他GATK4参数，注意添加单引号，示例：pcr参数 '--max-reads-per-alignment-start 0 --dont-use-soft-clipped-bases true --kmer-size 10 --kmer-size 15 --kmer-size 20 --kmer-size 25 --kmer-size 30 --kmer-size 35 --kmer-size 40 --max-num-haplotypes-in-population 256 '
                             其他参数：--max-assembly-region-size 200  --bam-writer-type ALL_POSSIBLE_HAPLOTYPES -L test.bed
         --help/-h           查看帮助文档



        [注意事项 1] 输入vcf格式示例：前两列描述信息必须存在，目标位点只需要填写 CHROM POS REF ALT 共4个内容即可，其他内容随便写
##fileformat=VCFv4.2									
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	genotype
1	3800242	rs3205087	A	G	100	.	*	GT:AD:DP:GQ:PL	0/1:233,223:456:99:4829,0,5327
1	59147926	rs12139511	T	C	100	.	*	GT:AD:DP:GQ:PL	0/1:233,223:456:99:4829,0,5327
1	86557967	rs11161732	G	A	100	.	*	GT:AD:DP:GQ:PL	0/1:233,223:456:99:4829,0,5327
1	118644430	rs10754367	A	G	100	.	*	GT:AD:DP:GQ:PL	0/1:233,223:456:99:4829,0,5327
1	154744807	rs1051614	C	G	100	.	*	GT:AD:DP:GQ:PL	0/1:233,223:456:99:4829,0,5327
1	158612236	rs863931	A	G	100	.	*	GT:AD:DP:GQ:PL	0/1:233,223:456:99:4829,0,5327
        
        [注意事项 2] 所有样本的vcf合并命令示例
bcftools merge -O z sample1.final.vcf sample2.final.vcf sample3.final.vcf -o allsample.final.vcf.gz
bcftools index --tbi allsample.final.vcf.gz
    \n";
    return $info;
}

