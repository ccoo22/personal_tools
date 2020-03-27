use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Parallel::ForkManager;
$|=1;

my $gtfToGenePred = "/home/genesky/software/ucsc/bin/gtfToGenePred";
my $genePredToGtf = "/home/genesky/software/ucsc/bin/genePredToGtf";

my ($gtf, $outputdir, $genome, $keep, $help);
GetOptions(
	"gtf|g=s"        => \$gtf,
	"outputdir|o=s"  => \$outputdir,
	"genome=s"       => \$genome,
	"keep!"           => \$keep,
	"help|h!"	     => \$help,
);
die help() if (defined $help or (not defined $gtf or not defined $outputdir));


################################## 主流程 ################
mkdir $outputdir if(not -e $outputdir);
my $input_name = `basename $gtf`;
   $input_name =~s/[\r\n]//g;

###############
# （1） gtf -> refgene -> 保留编码转录本 -> gtf
###############
my $tmp_refgene      = "$outputdir/$input_name.refgene.txt";
my $tmp_refgene_mrna = "$outputdir/$input_name.refgene.mRNA.txt";
my $tmp_mrna_gtf     = "$outputdir/$input_name.refgene.mRNA.gtf";

print "gtf to refgene\n";
system("$gtfToGenePred  $gtf $tmp_refgene -genePredExt -geneNameAsName2");

print "select mRNA in refgene\n";
system("awk '{if(\$6 != \$7) print \$0}'  $tmp_refgene  > $tmp_refgene_mrna");

print "mRNA refgene to gtf with utr3 flag\n";
system("$genePredToGtf  file $tmp_refgene_mrna $tmp_mrna_gtf -utr");

###############
# (2) 读入基因组数据
###############
my %hashFasta;
print "Read $genome ...";
my $IN = Bio::SeqIO->new(-file => $genome, -format=>'Fasta') or die "Could not open up file $genome: $!";
while(my $inSeq = $IN->next_seq)
{
    my $id  = $inSeq->display_id;
    my $seq = $inSeq->seq;
    $hashFasta{$id}= $seq;
}
print "OK\n";

###############
# (3) 获取3UTR序列
###############
print "get utr3 seq ... ";
my %hash3utr;
open GTF, "$tmp_mrna_gtf" or die "Could not open $tmp_mrna_gtf!";
while(<GTF>){
	$_=~s/[\r\n]//g;
	my ($chr, $region, $start, $end, $strand, $describe) = (split /\t/, $_)[0,2,3,4,6,8];
	next if($region ne '3UTR');
	my $seq = substr($hashFasta{$chr}, $start-1, $end - $start + 1);
	if($strand eq '-')
	{
		$seq = reverse $seq;
		$seq=~tr/ATCGatcg/TAGCtagc/;
	}
	   
	my @des = split /;/, $describe;
	my ($gene_name, $trans_id);
	foreach my $tmp(@des){
		$gene_name = $tmp if($tmp =~ /gene_name/);
		$trans_id  = $tmp if($tmp =~ /transcript_id/);
	}
	$gene_name =~ s/gene_name//g;
	$gene_name =~ s/"//g;	
	$gene_name =~ s/\s//g;
	
	$trans_id =~ s/transcript_id//g;
	$trans_id =~ s/"//g;	
	$trans_id =~ s/\s//g;
	
	my $title = "$trans_id|$gene_name";
	$hash3utr{$title} .= $seq;
}
close GTF;

# 输出 3utr序列
my $utr3_fa = "$outputdir/$input_name.utr_3.fa";
open FASTA, ">$utr3_fa";
foreach my $title(sort keys %hash3utr)
{
	print FASTA ">$title\n$hash3utr{$title}\n";
}
close FASTA;

print "OK\n";

# 清理中间结果
system("rm -f $tmp_refgene $tmp_refgene_mrna $tmp_mrna_gtf") if(not defined $keep);

print "$utr3_fa\n";

#################################  子程序
sub help{
    my $info = "
Program: 从gtf文件提取UTR3序列
Version: 2020-03-27
Contact: 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        --gtf/-g          输入gtf文件 例如： /home/ganb/work/research/mRNA_quantify/database/gencode.v30.annotation.gtf
        --outputdir/-o    输出路径
        --genome          参考基因组文件,eg: /home/genesky/database/ucsc/hg38/genome/hg38.fa
        --keep            是否保留中间临时文件，默认： 删除
        --help/-h         查看帮助文档
    \n";
    return $info;  
}

