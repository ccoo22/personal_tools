# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};
 
# 检测 -> 脚本输入
my ($input_fastq, $seq_name_list, $if_help);
GetOptions(
    "input_fastq|i=s"   => \$input_fastq,
    "seq_name_list|l=s" => \$seq_name_list,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $input_fastq or not defined $seq_name_list));
###################################################################### 主程序
my %hashSeqName = get_seq_name($seq_name_list);


open FASTQ, $input_fastq if($input_fastq !~ /.gz$/);
open FASTQ, "gzip -cd $input_fastq|" if($input_fastq =~ /.gz$/);
while( my $line1 = <FASTQ>)
{
    my $line2 = <FASTQ>;
    my $line3 = <FASTQ>;
    my $line4 = <FASTQ>;
    my ($seq_name) = split /\s+/, $line1;
        $seq_name =~ s/[\r\n]//g;
        $seq_name =~ s/^@//;
    next if(not exists $hashSeqName{$seq_name});

    print "$line1$line2$line3$line4";
    $hashSeqName{$seq_name}++;
}
close FASTQ;

my @losts = grep{ $hashSeqName{$_} == 0 } keys %hashSeqName;
print "Lost fasta : @losts\n " if(@losts > 0);



###################################################################### 子函数

sub get_seq_name{
    my $list_file = shift @_;
    my %hashSeqName;

    open LIST, $list_file;
    while(<LIST>)
    {
        next if($_!~/\w/);
        $_=~s/[\s\r\n]//g;
        $hashSeqName{$_} = 0;
    }
    close LIST;
    return %hashSeqName;
}

sub help{
    my $info = "
Program: get seq from fastq， 从fasta文件中获取指定序列
Version: 2019-02-26
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_fastq/-i    fastq文件
         --seq_name_list/-l  要获取的序列名称文件，不能有'\@'符号
         --help/-h           查看帮助文档
    \n";
    return $info;
}

