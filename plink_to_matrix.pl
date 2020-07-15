$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
plink to matrix
Version: v1.0 2020-05-28
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 变量
my $DEFAULT_SOFT_PLINK   = "/home/genesky/software/plink/1.07/plink";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input_plink, $output, $miss_blank, $SOFT_PLINK, $if_help);
GetOptions(
	"input_plink|i=s"     => \$input_plink,
	"output|o=s"          => \$output,

	"miss_blank|m"        => \$miss_blank,

	"plink=s"             => \$SOFT_PLINK,
	"help|h"              => \$if_help,
);
die "
Options: 必填

        --input_plink/-i     原始plink格式数据输入前缀， 例如sample.ped, sample.map文件， 则输入 -i sample
                             也可以是plink的二进制文件，例如 sample.bed,sample.bim,sample.fam, 则输入 -i sample
        --output/-o          输出文件，例如： out.txt

Options: 可选
        --miss_blank/-m      把plink的缺失值 0 0 替换为空值
        --plink              更改软件 plink 版本 (default: '$DEFAULT_SOFT_PLINK')
        --help/-h            查看帮助文档
\n" if (defined $if_help or not defined $input_plink or not defined $output);

$SOFT_PLINK   = $DEFAULT_SOFT_PLINK if (not defined $SOFT_PLINK);
###################################################################### 初始化
my $DATA_TIME = `date +\"\%Y-\%m-\%d \%H:\%M.\%S\"`;
my $RUN_INFO = "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
[SET] 软件 plink : $SOFT_PLINK
";

print $RUN_INFO;

###################################################################### 主程序

my $ped_file = "$input_plink.ped";
my $map_file = "$input_plink.map";

# 二进制文件转换为纯文本
my $is_binary = 0;
if(is_file_ok($ped_file, $map_file) == 0)
{   
    print "[message] 未检测到plink的纯文本文件，输入文件应该是二进制文件，尝试将输入文件转换为纯文本\n";
    system("$SOFT_PLINK --bfile $input_plink --recode --out $output.plink --noweb");
    $ped_file = "$output.plink.ped";
    $map_file = "$output.plink.map";
    $is_binary = 1;
}

# 数据提取
my %hashMap = read_map($map_file);
my %hashGeno = read_ped($ped_file, \%hashMap);

# 输出
print "output result\n";
my @samples = sort {$hashGeno{'sample'}{$a} <=> $hashGeno{'sample'}{$b}} keys %{$hashGeno{'sample'}};
my @snvids = sort {$hashMap{$a}{'sort_row'} <=> $hashMap{$b}{'sort_row'}} keys %hashMap;

open OUT, ">$output";
print OUT "SNV\tchr\tpos\t" . (join "\t", @samples) . "\n";
foreach my $snvid(@snvids)
{
    my @genos = map{ $hashGeno{'data'}{$_}{$snvid} } @samples;
    print OUT "$snvid\t$hashMap{$snvid}{'chr'}\t$hashMap{$snvid}{'pos'}\t" . (join "\t", @genos) . "\n";
}
close OUT;

print "result file : $output\n";

system("rm $ped_file $map_file") if($is_binary == 1);  # 删除中间文件
 
###################################################################### 子程序
sub read_ped{
    my $ped_file = shift @_;
    my $hashMap  = shift @_;

    print "read ped file\n";

    my %hashGeno;
    my @snvids = sort {$hashMap->{$a}{'sort_row'} <=> $hashMap->{$b}{'sort_row'}} keys %$hashMap;
    open PED, $ped_file;
    while(<PED>)
    {   
        $_=~s/[\r\n]//g;
        my ($fid, $iid, $pid, $mid, $sex, $pheno, @alleles) = split /\s+/, $_;
        foreach my $col(0..$#snvids)
        {
            my $allele1 = $alleles[2*$col];
            my $allele2 = $alleles[2*$col + 1];
            my $geno = "$allele1/$allele2";
               $geno = "" if(defined $miss_blank and $allele1 eq '0' and $allele2 eq '0');
            $hashGeno{'data'}{$iid}{$snvids[$col]} = $geno;
            $hashGeno{'sample'}{$iid} = $.;
        }
    }
    close PED;
    return %hashGeno;
}


sub read_map{
    my $map_file = shift @_;

    print "read map file\n";
    my %hashMap;
    open MAP, $map_file;
    while(<MAP>)
    {
        $_=~s/[\r\n]//g;
        my ($chr, $id, $tmp, $pos) = split /\s+/, $_;
        $hashMap{$id}{'chr'} = $chr;
        $hashMap{$id}{'pos'} = $pos;
        $hashMap{$id}{'sort_row'} = $.;
    }
    close MAP;
    return %hashMap;
}

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}