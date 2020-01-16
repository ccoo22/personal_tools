# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
$|=1;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};
 
# 检测 -> 脚本输入
my ($genome, $refgene, $if_help);
GetOptions(
    "genome|g=s"   => \$genome,
    "refgene|r=s"  => \$refgene,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $genome or not defined $refgene));
###################################################################### 主程序
my %hashGenome = read_fasta($genome);
my %hashCode   = get_homo_MT_code();
my %hashCodeCount; # 密码子计数

open REFGENE, $refgene;
while(<REFGENE>)
{
    $_ =~ s/[\r\n]//g;
    my ($tmp, $t_ID, $chr, $strand, $t_start, $t_end, $cd_start, $cd_end, $cd_count, $exon_start_list, $exon_end_list, $phase, $gene) = split /\t/, $_;
    next if($cd_start == $cd_end); # 非mRNA
    $cd_start++; # 向后移动一个位置，因为refGene表示区域时，起点是不要的
    my $ref_seq = ""; # 参考基因组序列
    if(exists $hashGenome{$chr})
    {
        $ref_seq = $hashGenome{$chr};
    }elsif(exists $hashGenome{"chr$chr"})
    {
        $ref_seq = $hashGenome{"chr$chr"};
    }else
    {
        $chr =~ s/chr//i;
        $ref_seq = $hashGenome{$chr} if(exists $hashGenome{$chr});
    }

    my $code_seq        = substr($ref_seq, $cd_start - 1, $cd_end - $cd_start + 1);
    my $code_seq_length = length($code_seq);
    my $code_count      = sprintf "%0.2f", $code_seq_length/3; # 密码子数量
    my $remain          = $code_seq_length % 3;

    if($strand eq '-')
    {
        $code_seq = reverse $code_seq;
        $code_seq =~ tr/acgtACGT/tgcaTGCA/;
    }

    

    my $protein = translate($code_seq, \%hashCode);
    count_code($code_seq, \%hashCodeCount);

    print "$gene\t$t_ID\t$strand\t$code_seq_length($code_count, remain = $remain)\t$code_seq\t$protein\n";

}
close REFGENE;

foreach my $code(sort keys %hashCodeCount)
{
    print "$code\t$hashCodeCount{$code}\n";
}



###################################################################### 子函数
sub count_code{
    my $code_seq      = shift @_;
    my $hashCodeCount = shift @_;
    while ($code_seq =~ /(...)/g) 
    {
        $hashCodeCount->{$1}++;
    }

}

sub translate{
    my $code_seq = shift @_;
    my $hashCode = shift @_;

    $code_seq = uc $code_seq;
    my $protein = "";
    my $count = 0;
    while ($code_seq =~ /(...)/g) 
    {
        $count++;
        if($count == 1)
        {
            my $amino_acid = exists $hashCode{'InitiationCodon'}{$1} ? $hashCode{'InitiationCodon'}{$1} : 'X';
            $protein   .= $amino_acid;
            next;
        }

        if (defined $hashCode{'Code'}{$1}) {
            $protein .= $hashCode{'Code'}{$1};
        } else {
            $protein .= "X";
        }
    }

    return $protein;
}

sub read_fasta{
    my $fasta_file = shift @_;

    my %hashFasta;
    my $fastaIn     = Bio::SeqIO->new(-file => $fasta_file, -format=>'Fasta') or die "Could not open up file $fasta_file: $!";
    while(my $inSeq = $fastaIn->next_seq)
    {
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;
        $hashFasta{$id} = $seq;
    }
    return %hashFasta;
}

sub get_homo_MT_code{
    my $code_map = "
TTT F Phe
TTC F Phe
TTA L Leu
TTG L Leu
CTT L Leu
CTC L Leu
CTA L Leu
CTG L Leu
ATT I Ile
ATC I Ile
ATA M Met
ATG M Met
GTT V Val
GTC V Val
GTA V Val
GTG V Val
TCT S Ser
TCC S Ser
TCA S Ser
TCG S Ser
CCT P Pro
CCC P Pro
CCA P Pro
CCG P Pro
ACT T Thr
ACC T Thr
ACA T Thr
ACG T Thr
GCT A Ala
GCC A Ala
GCA A Ala
GCG A Ala
TAT Y Tyr
TAC Y Tyr
TAA * Ter
TAG * Ter
CAT H His
CAC H His
CAA Q Gln
CAG Q Gln
AAT N Asn
AAC N Asn
AAA K Lys
AAG K Lys
GAT D Asp
GAC D Asp
GAA E Glu
GAG E Glu
TGT C Cys
TGC C Cys
TGA W Trp
TGG W Trp
CGT R Arg
CGC R Arg
CGA R Arg
CGG R Arg
AGT S Ser
AGC S Ser
AGA * Ter
AGG * Ter
GGT G Gly
GGC G Gly
GGA G Gly
GGG G Gly
    ";
    my %hashCode;
    foreach(split /\n/, $code_map)
    {
        next if($_ !~ /\w/);
        my ($code, $abbreviate, $amino_acid) = split /\s+/, $_;
        $hashCode{'Code'}{$code} = $abbreviate;
    }

    # 起始密码子候选
    $hashCode{'InitiationCodon'}{'ATT'} = 'M';
    $hashCode{'InitiationCodon'}{'ATG'} = 'M';
    $hashCode{'InitiationCodon'}{'ATA'} = 'M';

    return %hashCode;
}

sub help{
    my $info = "
Program: get refGene code count， 获取refGene密码子数量
Version: 2019-01-27
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --genome/-g   基因组fasta文件, example: /home/pub/database/Human/MT/mt.fa
         --refgene/-r  refGene文件, example: /home/pub/database/Human/MT/Annotation/hg38_MT_refGene.txt
         --help/-h     查看帮助文档
    \n";
    return $info;
}

