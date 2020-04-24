# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};
 
# 定义 -> 核心变量
my $SOFT_BOWTIE_BUILD   = "/home/genesky/software/bowtie2/2.3.4.3/bowtie2-build";
my $SOFT_GFF3TOGENEPRED = "/home/genesky/software/ucsc/bin/gff3ToGenePred";
my $SOFT_GENEPREDTOBED  = "/home/genesky/software/ucsc/bin/genePredToBed";
my $SOFT_ANNOVAR        = "/home/genesky/software/annovar/2018Apr16";
my $SOFT_BWA            = "/home/genesky/software/bwa/0.7.17/bwa";
my $SOFT_SAMTOOLS       = "/home/genesky/software/samtools/1.10/samtools";

# 检测 -> 脚本输入
my ($fa, $gff, $output_dir, $prefix, $if_help);
GetOptions(
    "fa=s"           => \$fa,
    "gff=s"          => \$gff,
    "output_dir=s"   => \$output_dir,
    "prefix=s"       => \$prefix,
    "help|h"         => \$if_help,
);
die help() if(defined $if_help or (not defined $fa or not defined $gff or not defined $output_dir or not defined $prefix));
###################################################################### 主程序
my $bowtie2_index_dir = "$output_dir/bowtie2_index";
my $annotation_dir    = "$output_dir/annotation";
my $bwa_index_dir     = "$output_dir/bwa_index";
make_dir($output_dir, $bowtie2_index_dir, $annotation_dir, $bwa_index_dir);

# (1) 拷贝输入文件到结果目录
print "copy \n";
my $genome = "$output_dir/$prefix.fa";
my $gff3   = "$output_dir/$prefix.gff3";
system("cp $fa $genome");
system("cp $gff $gff3");

# (2) 创建bowtie2索引
print "build bowtie2 index \n";
chdir $bowtie2_index_dir;
system("ln -s ../$prefix.fa ./");
system("$SOFT_BOWTIE_BUILD --threads 10 $prefix.fa $prefix.fa");
chdir PWD;
# (3) 创建bwa索引
print "build bwa index \n";
my $fai            = "$bwa_index_dir/$prefix.fa.fai";
my $chrom_size     = "$bwa_index_dir/$prefix.chrom_size";
my $effective_size = "$bwa_index_dir/$prefix.effective_size";
chdir $bwa_index_dir;
system("ln -s ../$prefix.fa ./");
system("$SOFT_BWA index $prefix.fa  ");
system("$SOFT_SAMTOOLS faidx $prefix.fa");
chdir PWD;

system("cut -f 1,2 $fai > $chrom_size");
system("awk 'BEGIN{size = 0} { size = size + \$2} END{ print size}' $chrom_size| grep -v -P 'Mt|Pt'   > $effective_size");


# (4) refgene 注释数据库创建
print "buil gene database\n";
my $refgene   = "$annotation_dir/$prefix\_refGene.txt";
my $refgenefa = "$annotation_dir/$prefix\_refGeneMrna.fa";
my $bed12     = "$annotation_dir/$prefix\_refGene.bed12"; 
my $symbol    = "$annotation_dir/$prefix\_symbol.txt";
my $tss_bed   = "$annotation_dir/$prefix\_tss.bed.gz";
my $beta_db   = "$annotation_dir/$prefix\_beta.db";
system("$SOFT_GFF3TOGENEPRED $gff3 $refgene.tmp -useName");
system("sed 's/transcript://'  $refgene.tmp > $refgene");  # 去掉奇怪的字符
system("rm  $refgene.tmp  "); 
system("perl $SOFT_ANNOVAR/retrieve_seq_from_fasta.pl -format refGene --seqfile $genome $refgene --outfile  $refgenefa");  # annovar注释数据库
system("$SOFT_GENEPREDTOBED $refgene $bed12");  # BED12
system("cut -f 1,12 $refgene > $symbol");
system(" less $refgene|  awk '{end = \$4 + 1; if(\$2 != 'Mt' || \$2 != 'Pt') print \$2\"\\t\"\$4\"\\t\"end\"\\t\"\$1\"\\t0\\t\"\$3  }' | gzip > $tss_bed");
system("cut -f 1,2,3,4,5,12 $refgene > $beta_db.tmp");
system("sed '1i\\#name\\tchrom\\tstrand\\ttxStart\\ttxEnd\\tname2' $beta_db.tmp > $beta_db ");
system("rm $beta_db.tmp ");
 
###################################################################### 子函数
 
# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}

sub help{
    my $info = "
Program: ATAC 数据库创建
Version: 2020-04-14
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --fa               基因组fasta文件
         --gff              基因数据库文件，gff3格式
         --output_dir       结果输出目录
         --prefix           输出文件的前缀
         --help/-h          查看帮助文档
    \n";
    return $info;
}

