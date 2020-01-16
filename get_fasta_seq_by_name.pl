# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};
 
# 检测 -> 脚本输入
my ($input_fasta, $seq_name_list, $if_help);
GetOptions(
    "input_fasta|i=s"   => \$input_fasta,
    "seq_name_list|l=s" => \$seq_name_list,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $input_fasta or not defined $seq_name_list));
###################################################################### 主程序
my %hashSeqName = get_seq_name($seq_name_list);

my $is_find = 0;

open FASTA, $input_fasta;
while( my $info = <FASTA>)
{
    if($info =~ /^>/)
    {
        my ($seq_name) = split /\s+/, $info;
        $seq_name =~ s/^>//g;
        if(exists $hashSeqName{$seq_name})
        {
            $hashSeqName{$seq_name}++;
            $is_find = 1;
        }
        else
        {
            $is_find = 0;
        }
    }
    print $info if($is_find == 1); # 直接输出到屏幕
}
close FASTA;

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
Program: get seq from fasta， 从fasta文件中获取指定序列
Version: 2019-01-22
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_fasta/-i    fasta文件
         --seq_name_list/-l  要获取的序列名称文件，不能有'>'符号
         --help/-h           查看帮助文档
    \n";
    return $info;
}

