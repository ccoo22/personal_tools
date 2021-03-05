use strict;
use warnings;
use File::Spec;
use Getopt::Long;
$|=1;

my $DEFAULT_SOFT_SAMTOOLS = "/home/genesky/software/samtools/1.10/samtools";

# 检测 -> 脚本输入
my ($bam, $soft_samtools, $if_help);
GetOptions(
    "bam=s"       => \$bam,
    "samtools=s"  => \$soft_samtools,
    "help|h"      => \$if_help,
);
die help() if(defined $if_help or (not defined $bam));
$soft_samtools = $DEFAULT_SOFT_SAMTOOLS if(not defined $soft_samtools);

###################################################################### 主程序

my $mismatch_count = 0;
open BAM, "$soft_samtools view $bam|";
while(<BAM>)
{
    my ($reads_name, $flag, $chr, $pos, $quality, $cigar, $chr_mate, $pos_mate, $distance, $seq, $seq_quality, @tags) = split /\t/, $_;
    foreach my $tag(@tags)
    {
        if($tag =~/^MD/)
        {
            my $count = $tag =~ tr/ATCG/ATCG/; 
            $mismatch_count += $count;
            next;
        }
    }
}
close BAM;
print "$bam\t$mismatch_count\n";

#########################################################################
sub help{
    my $info = "
Program: 统计bam文件中有多少mismatch碱基（基于bam的MD标签）
Version: 2021-02-23
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [必填]
         --bam         bam文件，例如： ./abc.bam

        [选填]
         --samtools    samtools软件路径，默认： $DEFAULT_SOFT_SAMTOOLS
         --help/-h     查看帮助文档
    \n";
    return $info;
}
