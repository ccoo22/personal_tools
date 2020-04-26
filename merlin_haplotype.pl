$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
基于merlin的 单倍型分析
Version: v1.0 2020-04-26
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
my $DEFAULT_SOFT_MERLIN  = "/home/genesky/software/merlin/1.1.2/executables";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input_plink, $snp_list, $output_dir, $SOFT_PLINK, $SOFT_MERLIN, $if_help);
GetOptions(
	"input_plink|i=s"     => \$input_plink,
	"snp_list|s=s"        => \$snp_list,
	"output_dir|o=s"      => \$output_dir,

	"plink=s"             => \$SOFT_PLINK,
	"merlin=s"            => \$SOFT_MERLIN,
	"help|h"              => \$if_help,
);
die "
Options: 必填

        --input_plink/-i                plink格式数据输入前缀， 例如sample.ped, sample.map文件， 则输入 -i sample
        --snp_list/-s                   snp名称列表，一列数据，没有表头。
                                        流程会从plink里提取这些snp数据，制作merlin格式数据，并根据遗传信息推断单倍型。
        --output_dir/-o                 结果输出路径

Options: 可选
        --plink                    更改软件 plink 版本 (default: '$DEFAULT_SOFT_PLINK')
        --merlin                   更改软件 merlin 版本 (default: '$DEFAULT_SOFT_MERLIN')
        --help/-h                  查看帮助文档

        注意：由于家系配置错误、样本量等问题，可能会分析失败。此时要仔细检查一下是否存在遗传异常的样本。
\n" if (defined $if_help or not defined $input_plink or not defined $snp_list or not defined $output_dir);

$SOFT_PLINK   = $DEFAULT_SOFT_PLINK if (not defined $SOFT_PLINK);
$SOFT_MERLIN  = $DEFAULT_SOFT_MERLIN if (not defined $SOFT_MERLIN);
 
$output_dir = File::Spec->rel2abs($output_dir);  # 转换为绝对路径
mkdir $output_dir if(not -e $output_dir);

###################################################################### 初始化


###################################################################### 主程序

print "(1) 提取SNP数据\n";
my $select_snp = "$output_dir/select_snp";
system("$SOFT_PLINK --file $input_plink --extract $snp_list --recode --out $select_snp --noweb");


print "(2) 生成merlin格式数据\n";
my $final = "$output_dir/final";
system("cp $select_snp.ped $final.ped.tmp");
system("sed 's/-9/0/' $final.ped.tmp > $final.ped && rm $final.ped.tmp");  # merlin不支持-9
system("cut -f1,2,3 $select_snp.map > $final.map.tmp ");
system("sed '1i\\CHROMOSOME\\tMARKER\\tPOSITION' $final.map.tmp > $final.map  && rm $final.map.tmp");
system(" awk '{print \"M\\t\"\$2}' $select_snp.map > $final.dat.tmp ");
system(" sed '1i\\A\\tX' $final.dat.tmp > $final.dat && rm  $final.dat.tmp");


print "(3) 单倍型分析  根据家系大小、snp数量，分析可能会很慢\n";
my $result_raw  = "$output_dir/haplotype";
system("$SOFT_MERLIN/merlin -d $final.dat -p $final.ped -m $final.map --best --bits 36 --prefix $result_raw --horizontal  ");

print "(4) 结果整理，方便使用阅读\n";
my $result = "$output_dir/haplotype.txt";
format_result("$final.ped", "$result_raw.chr", $result);


print "[final result] $result\n";
###################################################################### 子程序

sub format_result{
    my $ped = shift @_;
    my $raw = shift @_;
    my $result = shift @_;

    # 读入样本患病状态
    my %hashPheno;
    open PED, $ped;
    while(<PED>)
    {
        my ($fid, $iid, $pid, $mid, $sex, $pheno, $tmp) = split /\s+/, $_, 7;
        $hashPheno{$iid} = "$pid\t$mid\t$sex\t$pheno";
    }
    close PED;

    # 读入单倍型
    my %hashHap;
    open HAP, $raw;
    <HAP>;
    while(my $info1 = <HAP>)
    {   
        last if($info1 !~/\w/);

        my $info2 = <HAP>;

        $info1=~s/[\r\n]//g;
        $info1=~s/^\s+//;
        $info2=~s/[\r\n]//g;
        $info2=~s/^\s+//;

        my ($iid1, $tmp1, $hap1) = split /\s+/, $info1, 3;
        my ($iid2, $tmp2, $hap2) = split /\s+/, $info2, 3;
        $hap1=~s/\s//g;
        $hap2=~s/\s//g;
        $hashHap{$iid1}{'HAP1'} = $hap1;
        $hashHap{$iid1}{'HAP2'} = $hap2;
    }
    close HAP;

    # 输出
    open RESULT, ">$result";
    print RESULT "Sample\tHAP1\tHAP2\tfatherID\tmotherID\tsex\tpheno\n";
    foreach my $sample(sort keys %hashHap)
    {
        print RESULT "$sample\t$hashHap{$sample}{'HAP1'}\t$hashHap{$sample}{'HAP2'}\t$hashPheno{$sample}\n";
    }
    close RESULT;
}

 