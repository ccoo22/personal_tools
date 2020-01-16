# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
$|=1;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $bedtools      = "/home/genesky/software/bedtools/2.28.0/bin/bedtools";

# 检测 -> 脚本输入
my ($input, $gtf, $output_dir, $if_help);
GetOptions(
    "input|i=s"      => \$input,
    "gtf|g=s"        => \$gtf,
    "output_dir|o=s" => \$output_dir,
    "help|h"         => \$if_help,
);
die help() if(defined $if_help or (not defined $input or not defined $gtf or not defined $output_dir));
###################################################################### 主程序

# (1) 读入需要的基因信息
my %hashGene = read_gene($input); 

# (2) 读入GTF文件
my %hashGTF = read_gtf($gtf, \%hashGene);

# (3) 信息输出
my $bed            = "$output_dir/gene.mrna.region.exon.bed"; # 外显子区域BED
my $bed_sort       = "$output_dir/gene.mrna.region.exon.sort.bed";  # 外显子区域BED排序
my $transcript_txt = "$output_dir/gene.mrna.region.transcript.txt";  # 转录本区域信息
my $exon_txt       = "$output_dir/gene.mrna.region.exon.txt";  # 外显子区域信息

open BED, ">$bed"; # bed区域，用于从bam中提取外显子覆盖深度
open TRANSCRIPT, ">$transcript_txt";# Gene    Mrna    Chr     Start   End
open EXON, ">$exon_txt"; # Gene    Mrna    Exon    Chr     Start   End

print TRANSCRIPT "Gene\tMrna\tChr\tStart\tEnd\n";
print EXON "Gene\tMrna\tExon\tChr\tStart\tEnd\n";

foreach my $gene(sort keys %hashGene)
{
    foreach my $trans(sort keys %{$hashGene{$gene}})
    {   
        if(not exists $hashGTF{$gene}{$trans})
        {
            print "在GTF中，没有找到 [$gene,$trans] 的信息, 请仔细核查\n";
            next;
        }

        print TRANSCRIPT "$gene\t$trans\t$hashGTF{$gene}{$trans}{'Chr'}\t$hashGTF{$gene}{$trans}{'Start'}\t$hashGTF{$gene}{$trans}{'End'}\n";
        foreach my $exon_num(sort {$a <=> $b} keys %{$hashGTF{$gene}{$trans}{"Exon"}})
        {   
            my $chr   = $hashGTF{$gene}{$trans}{'Exon'}{$exon_num}{'Chr'};
            my $start = $hashGTF{$gene}{$trans}{'Exon'}{$exon_num}{'Start'};
            my $end   = $hashGTF{$gene}{$trans}{'Exon'}{$exon_num}{'End'};
            print EXON "$gene\t$trans\t$exon_num\t$chr\t$start\t$end\n";
            print BED "$chr\t$start\t$end\t+\t$gene|$trans|Exon$exon_num\n";
        }
    }
}
close TRANSCRIPT;
close EXON;
close BED;

# 排序
system("$bedtools  sort -i $bed |$bedtools merge -i - > $bed_sort");

print "Generate Files:\n";
print "外显子区域bed文件：       $bed \n";
print "外显子区域排序后的bed文件：$bed_sort\n";
print "转录本区域文件：          $transcript_txt\n";
print "外显子区域详细文件：       $exon_txt\n";

###################################################################### 子函数
sub read_gene{
    my $input = shift @_;
    my %hashGene;
    open INPUT, $input;
    while(<INPUT>)
    {
        $_ =~s /[\r\n]//g;
        next if($_!~/\w/);
        my ($gene, $trans) = split /\t/, $_;
        $hashGene{$gene}{$trans}++;
    }
    close INPUT;
    return %hashGene;
}


sub read_gtf{
    my $gtf      = shift @_;
    my $hashGene = shift @_;

    my %hashGTF;
    print "Reading $gtf ... ";
    open GTF,$gtf;
    while(<GTF>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/ or $_ =~ /^#/);

        my ($chr, $source, $RegionType, $start, $end, $tmp1, $strand, $tmp2, $anno) = split /\t/, $_;
        next if($RegionType ne 'exon'); # 只要exon信息
        my %hashAnno = split_anno($anno); # 拆分注释信息，便于后续使用
        my $transcript_id = (split /\./, $hashAnno{'transcript_id'})[0];
        my $gene = $hashAnno{'gene_name'};
        my $exon = $hashAnno{'exon_number'};

        next if(!defined $transcript_id or !defined $gene or !defined $exon);  # 必须都有定义 
        next if(not exists $hashGene->{$gene}{$transcript_id}); # 必须是我们需要的信息

        # 记录exon信息
        $hashGTF{$gene}{$transcript_id}{"Exon"}{$exon}{"Chr"}   = $chr;
        $hashGTF{$gene}{$transcript_id}{"Exon"}{$exon}{"Start"} = $start;
        $hashGTF{$gene}{$transcript_id}{"Exon"}{$exon}{"End"}   = $end;
        
        # 记录转录本区域信息
        $hashGTF{$gene}{$transcript_id}{"Chr"}    = $chr;
        $hashGTF{$gene}{$transcript_id}{"Start"}  = $start if(!exists $hashGTF{$gene}{$transcript_id}{"Start"} or $hashGTF{$gene}{$transcript_id}{"Start"} > $start);
        $hashGTF{$gene}{$transcript_id}{"End"}    = $end   if(!exists $hashGTF{$gene}{$transcript_id}{"End"} or $hashGTF{$gene}{$transcript_id}{"End"} < $end);
        $hashGTF{$gene}{$transcript_id}{"Length"} = $hashGTF{$gene}{$transcript_id}{"End"} - $hashGTF{$gene}{$transcript_id}{"Start"} + 1;
    }
    close GTF;
    print "OK\n";
    return %hashGTF;
 
}
sub split_anno{
    my $anno = shift @_;
    my %hashAnno;
    my @infos = split /;/, $anno;
    foreach my $info(@infos)
    {
        next if($info!~/\w/);
        $info=~s/^\s+//;
        $info=~s/\s+$//;
        my ($name, $value)=split /\s+/, $info;
        $value=~s/\"//g;
        $hashAnno{$name} = $value;
    }
    return %hashAnno;
}

sub help{
    my $info = "
Program: get_gene_transcript_structure， 获取基因转录本的CDS区域信息
Version: 2019-12-03
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input/-i          输入文件，不要表头，第一列是基因名，第二列是转录本名(注意：转录本名不要加版本号，即：.后面的版本号去掉)
         --gtf/-g            GTF数据库，例如：/home/ganb/work/research/mRNA_quantify/database/gencode.v30.annotation.gtf
         --output_dir/-o     结果输出目录
         --help/-h           查看帮助文档
    \n";
    return $info;
}

