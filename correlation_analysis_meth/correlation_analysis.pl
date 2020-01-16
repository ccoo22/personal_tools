# 导入 -> 系统 package   
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Statistics::R;
use Excel::Writer::XLSX;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $METHOD       = "pearson";
my $SOFT_RSCRIPT = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $ANALYSIS_R   = SCRIPTDIR."/correlation_analysis.r";
my $table2excel  = "perl /home/pub/bin/NGS/chip/GATK4/tools/personal/table2excel.pl";

# 检测 -> 脚本输入
my ($meth, $expr, $method, $thread, $prefix, $if_help);
GetOptions(
	"meth|m=s"              => \$meth,
	"expr|e=s"              => \$expr,
	"corrlation-method|c=s" => \$method,
	"thread|t=s"            => \$thread,
	"prefix|p=s"            => \$prefix,
	"help|h"                => \$if_help,
);
die help() if (defined $if_help or not defined $meth or not defined $expr or not defined $prefix);
$method   = $METHOD   if (not defined $method or ($method ne "pearson" and $method ne "spearman"));
$thread = 10        if(not defined $thread);

###################################################################### 主程序
my $sampleCheck = checkSample($meth, $expr); #检测两个文件的样本列是否完全一样
die "Error:Check your sample names in two files!\n" if($sampleCheck == 0);

# 调用r脚本计算相关性
my $output_file = "$prefix.correlation.result.txt";
system("$SOFT_RSCRIPT $ANALYSIS_R -c $method -p $thread $meth $expr $output_file");

# 输出excel
my $excel = "$prefix.correlation.xlsx";
print "final : $excel\n";
system("$table2excel -i $output_file -o $excel -s $method ");



###################################################################### 子函数
sub help{
	my $info = "
Program: Correlation analysis between two files
Version: 2019-6-26
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options: 必填

        --meth/-m              用于相关性分析的 “甲基化数据”,行为样本,列为甲基化位点,样本名顺序与expr需完全一致)
        --expr/-e              用于相关性分析的 “表型数据”,行为样本,列为表型,样本名顺序与meth需完全一致)       
        --prefix/-p            结果结果前缀，例如: /home/ganb/result
		
Options: 可选
        --corrlation-method/-c     相关性分析方法[pearson(default)/spearman]
        --thread/-t                并行线程数量，默认：10
        --help/-h                  查看帮助文档
	\n";
	return $info;
}

sub checkSample{
	my $file1 = shift @_;
	my $file2 = shift @_;
	my @samples1;
	my @samples2;
	open FILE1,$file1;
	<FILE1>;
	while(<FILE1>){
		$_=~s/[\r\n]//g;
		push(@samples1, (split/\t/,$_)[0]);
	}
    close FILE1;
    open FILE2,$file2;
	<FILE2>;
	while(<FILE2>){
		$_=~s/[\r\n]//g;
		push(@samples2, (split/\t/,$_)[0]);
	}
    close FILE2;

    # 是否一致
    my $is_same = 1;
       $is_same = 0 if(scalar(@samples1) != scalar(@samples2)); # 数量不一致
    foreach my $col(0..$#samples1)
    {
    	next if($samples1[$col] eq $samples2[$col]);
    	$is_same = 0; # 样本名不是一一对应
    }

    return $is_same;
}
