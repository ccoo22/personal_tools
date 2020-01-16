# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $blastn      = "/home/ganb/soft/ncbi-blast-2.7.1+/bin/blastn";
my $makeblastdb = "/home/ganb/soft/ncbi-blast-2.7.1+/bin/makeblastdb";
my $nt_db       = "/home/pub/database/NCBI/NT/latest/nt";

# 检测 -> 脚本输入
my ($input_fasta, $input_db, $format, $if_help);
GetOptions(
    "input_fasta|i=s"   => \$input_fasta,
    "db|d=s"            => \$input_db,
    "format|f=s"        => \$format,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $input_fasta or not defined $input_db or not defined $format));
###################################################################### 主程序
 
# (1) 建库
my $db = $nt_db;
if($input_db ne 'nt')
{
    $db = $input_db;
    system("$makeblastdb -in $input_db -dbtype nucl") if(not -e "$db.nhr"); 
}

# (2) 比对
my $format_code = 0;
   $format_code = 6 if($format eq 'tab');
system("$blastn -task blastn -query $input_fasta -evalue 0.00001 -db $db -outfmt $format_code -num_threads 4 -dust no");






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
Program: blast+ blastn， 序列比对
Version: 2019-01-22
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_fasta/-i    fasta文件
         --db/-d             数据库，自建数据库：.../database.fa  或 关键词 ：nt
         --format/-f         数据比对格式，tab/pairwise
         --help/-h           查看帮助文档
    \n";
    return $info;
}

