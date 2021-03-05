# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use Cwd qw( abs_path );

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 检测 -> 脚本输入
my ($input, $output, $length, $if_help);
GetOptions(
    "input|i=s"   => \$input,
    "output|o=s"  => \$output,
    "length|l:i"  => \$length,
    "help|h"      => \$if_help,
);
die help() if(defined $if_help or (not defined $input or not defined $output or not defined $length));
 
$input  = Cwd::abs_path($input);
$output = Cwd::abs_path($output);

###################################################################### 主程序

my $IN = Bio::SeqIO->new(-file => $input, -format=>'Fasta') or die "Could not open up file $input: $!";
open OUT, ">$output";
my $count = 0;
my $ok    = 0;
while(my $inSeq = $IN->next_seq)
{   
    $count++;
    my $id  = $inSeq->display_id;
    my $seq = $inSeq->seq;
    my $seq_length = length($seq);
    if($seq_length >= $length)
    {
        print OUT ">$id\n$seq\n";
        $ok++;
    }
}
close OUT;

print "Total seq: $count\n";
print "Pass filter seq: $ok\n";
###################################################################### 子函数
sub help{
    my $info = "
Program: 从fasta文件中，提取长度大于n的序列
Version: 2021-01-22
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input/-i   fasta序列输入
         --output/-o  输出文件
         --length/-l  最小长度， 例如  10000
         --help/-h   查看帮助文档
    \n";
    return $info;
}

