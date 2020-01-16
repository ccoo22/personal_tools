# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};
 
# 检测 -> 脚本输入
my ($input_fasta, $region, $if_help);
GetOptions(
    "input_fasta|i=s"   => \$input_fasta,
    "region|r=s"        => \$region,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $input_fasta or not defined $region));
###################################################################### 主程序
my %hashFasta = read_fasta($input_fasta);
my $lost = "";
foreach my $region_one(split /,/, $region)
{
    my ($chr, $start, $end) = split /[:-]/, $region_one;
    my $seq = "";
       $seq = substr($hashFasta{$chr}, $start - 1, $end - $start + 1) if(exists $hashFasta{$chr});
    print ">$region_one\n$seq\n";
    $lost .= "$region_one," if($seq eq '');
}

print "Lost fasta : $lost\n " if($lost ne '');



###################################################################### 子函数

# 读取fasta文件
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


sub help{
    my $info = "
Program: get seq from fasta by region， 从fasta文件中获取指定序列
Version: 2019-03-13
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_fasta/-i    fasta文件
         --region/-r         要获取的序列区域信息，多个区域用逗号分割，示例：1:1-10,2:100-200
         --help/-h           查看帮助文档
    \n";
    return $info;
}

