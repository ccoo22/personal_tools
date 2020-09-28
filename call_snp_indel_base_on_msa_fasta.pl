$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
基于多重序列比对MSA文件，识别SNP/INDEL
Version: v1.0 2020-09-25
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use Cwd 'abs_path';

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input, $output, $if_help);
GetOptions(
    "input|i=s"            => \$input,
    "output|o=s"           => \$output,

    "help|h"               => \$if_help,
);
die "
Options: 必填
        --input/-i                多重比对fasta文件
                                  第一条序列作为参考序列

        --output/-o               突变结果输出

Options: 可选
        --help/-h                 查看帮助文档

\n" if (defined $if_help or not defined $input or not defined $output);

$input  = abs_path($input);
$output = abs_path($output);
###################################################################### 主程序
my %hashFasta = read_fasta($input, 'seq');
my @seq_ids = sort {$hashFasta{$a}{'sort_order'}<=>$hashFasta{$b}{'sort_order'}} keys %hashFasta;
my $reference_id  = $seq_ids[0];
my $reference_seq = $hashFasta{$reference_id}{'seq'};
my @ref_bases     = split //, $reference_seq;

my @heads = qw(SeqName type insertion deletion snp seq);
open OUTPUT, ">$output";
print OUTPUT (join "\t", @heads) . "\n";
foreach my $seq_id(@seq_ids)
{
    next if($seq_id eq $reference_id);
    my @seq_bases     = split //, $hashFasta{$seq_id}{'seq'};

    # SNP/INDEL识别
    my %hashInfo;
    $hashInfo{'SeqName'} = $seq_id;
    $hashInfo{'seq'}     = $hashFasta{$seq_id}{'seq'};

    my $pos_ref   = 0;  # 参考序列位置
    my %hashType;
    my %hashInsert;
    foreach my $pos(0..$#ref_bases)
    {
        my $ref_base = $ref_bases[$pos];  # 参考碱基
        my $seq_base = $seq_bases[$pos];  # 测序碱基
        $pos_ref++ if($ref_base =~ /\w/);  # 参考序列位置
        next if($seq_base eq $ref_base);
        
        if($ref_base=~/\w/ and $seq_base =~ /\w/)  
        {   # mismatch
            $hashInfo{'snp'} .= "$ref_base->$seq_base;";
            $hashType{$pos_ref}{'M'} = 1;
        }elsif($ref_base=~/\w/ and $seq_base !~ /\w/)  
        {   # delete
            $hashInfo{'deletion'} .= "$ref_base;";
            $hashType{$pos_ref}{'D'} = 1;
        }else{  
            # insert
            $hashInsert{$pos_ref} .= $seq_base;
            $hashType{$pos_ref}{'I'} = 1;
        }
    }
    # insertion信息汇总
    foreach my $pos_insert(sort {$a<=>$b} keys %hashInsert)
    {   
        $hashInfo{'insertion'} .= "$hashInsert{$pos_insert};";
    }
        
    # 突变位置信息汇总
    foreach my $pos_diff(sort {$a<=>$b} keys %hashType)
    {   
        foreach my $diff_type(sort keys %{$hashType{$pos_diff}} )
        {
            $hashInfo{'type'} .= $pos_diff.$diff_type.";";
        }
    }

    # 输出
    my @values = map{ my $value = exists $hashInfo{$_} ? $hashInfo{$_} : "-" } @heads;
    print OUTPUT (join "\t", @values) . "\n";
}
close OUTPUT;


###################################################################### 子程序

sub read_fasta{
    my $fasta = shift @_;
    my $type  = shift @_;
    my %hashFasta;
    my $IN = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
    my $count = 0;
    while(my $inSeq = $IN->next_seq)
    {   
        $count++;
        my $id  = $inSeq->display_id;
        my $seq = uc($inSeq->seq);
        if($type eq 'seq_count')
        {
            $hashFasta{$seq}++;
        }elsif($type eq 'seq'){
            $hashFasta{$id}{'seq'}        = $seq;
            $hashFasta{$id}{'sort_order'} = $count;
        }
        
    }
    return %hashFasta;
}
