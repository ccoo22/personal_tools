# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $gtfToGenePred         = "/home/genesky/software/ucsc/bin/gtfToGenePred";
my $retrieve_seq_from_fasta    = "/home/genesky/software/annovar/2018Apr16/retrieve_seq_from_fasta.pl";

# 检测 -> 脚本输入
my ($gtf, $genome, $prefix, $output_dir, $if_help);
GetOptions(
	"gtf=s" => \$gtf,
	"genome|g=s" => \$genome,
	"prefix|p=s" => \$prefix,
	"output_dir|o=s" => \$output_dir,
	"help|h" => \$if_help,
);
die help() if (defined $if_help or not defined $gtf or not defined $prefix);

###################################################################### 主程序

my $refGene0 = "$output_dir/$prefix\_refGene0.txt";
my $refGeneMrna = "$output_dir/$prefix\_refGeneMrna.fa";
my $refGene = "$output_dir/$prefix\_refGene.txt";

system "$gtfToGenePred -genePredExt $gtf $refGene0";
system "sed -i s/transcript\://g $refGene0";
system "sed -i s/gene\://g $refGene0";
system "nl $refGene0 > $refGene";
system "rm $refGene0";

if (defined($genome)){

    if (substr("$genome", -2) eq 'gz'){
        my $tmp = substr("$genome", 0, -3);
        gunzip $genome => $tmp or die "gunzip failed: $GunzipError\n";

        system "perl $retrieve_seq_from_fasta --format refGene --seqfile $tmp $refGene --out $refGeneMrna";
        system "rm $tmp";
    }else{
        system "perl $retrieve_seq_from_fasta --format refGene --seqfile $genome $refGene --out $refGeneMrna";
    }
}



###################################################################### 子函数

sub help{
    my $info = "
Program: GTF 转 refGene
Version: 2019-11-08
Contact: 341 zhangqing

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
         --gtf               GTF文件，例：Arabidopsis_thaliana.TAIR10.45.gtf，必填
         --genome/-g         基因组文件，例：Arabidopsis_thaliana.TAIR10.dna.toplevel.fa，选填
         --prefix/-p         结果文件前缀（物种名称或其首字母缩写），必填
         --output_dir/-o     结果输出目录
         --help/-h           查看帮助文档
    \n";
    return $info;
}