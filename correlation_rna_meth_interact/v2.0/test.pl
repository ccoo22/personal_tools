#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use Parallel::ForkManager;

#定义
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

my $scriptdir     = SCRIPTDIR;
my $bedtools      = "/home/pub/software/bedtools/bedtools2-2.27.1/bin/bedtools";
my $method        = "pearson";
my $r_script      = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $split_data_r  = "$scriptdir/utils/split_data.r";
my $correlation_r = "$scriptdir/utils/correlation.r";
my $table2excel   = "$scriptdir/utils/table2excel.pl";

my ($rna, $meth, $output_dir, $up_extend, $down_extend, $thread, $help);

GetOptions(
	'rna|r=s'               => \$rna,
	'meth|m=s'              => \$meth,
	'output_dir|o=s'        => \$output_dir,
	'corrlation-method|c=s' => \$method,
	'up_extend|up=s'		=> \$up_extend,
	'down_extend|down=s'	=> \$down_extend,
	'thread|T=s'            => \$thread,
	'help|h!' 		        => \$help,
);
die help() if (defined $help or (not defined $rna or not defined $meth or not defined $output_dir) );
$method   = $method   if (not defined $method or ($method ne "pearson" and $method ne "spearman"));
$thread   = 10 if(not defined $thread);
$up_extend   = 500000 if (not defined $up_extend );
$down_extend = 500000 if (not defined $down_extend );

###################################################################### 主程序

########## 
# (1) 目录准备
########## 
print "(1) 创建输出目录\n";
my $tmp_dir    = "$output_dir/tmp";
my $split_dir  = "$output_dir/split";
make_dir($output_dir, $tmp_dir, $split_dir);

########## 
# (2) 表达量、甲基化区域取交集
########## 
print "(2) 检测表达量与甲基化重叠区域\n";
my $rna_bed      = "$tmp_dir/rna.bed";
my $meth_bed     = "$tmp_dir/meth.bed";
my $interact_bed = "$tmp_dir/interact.bed";

get_rna_extend_bed($rna, $rna_bed, $up_extend, $down_extend);
system("awk -F '\\t' '{if(NR >1) print \$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$1}' $meth > $meth_bed");

system("$bedtools intersect -a $rna_bed  -b $meth_bed -wa -wb > $interact_bed");

die "[Warnings]  输入的数据中，转录本与甲基化没有区域重叠" if(is_file_ok($interact_bed) == 0);

###################################################################### 子函数

# 创建目录
sub make_dir{
    my @dir = @_;
    foreach(@dir){
        mkdir $_ if(not -e $_);
    }
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
# 进度条， 根据输入向量计算
sub process_bar_array{
    my $process_name   = shift @_; # 需要分析的对象
    my $process_arrays = shift @_; # 总列表
    my $process_all_count = @$process_arrays;# 需要分析的总列表数量
    my ($pos) = grep $$process_arrays[$_] eq $process_name, 0..$process_all_count-1;
    my $process_count = $pos+ 1; # 当前对象的位置
    my $process_perc = sprintf "%0.2f", 100 * $process_count/$process_all_count; # 进度
    my $windows = 100; # 窗口宽度
    my $finished = $windows * int($process_perc) / 100; # 完成
    my $unfinished = $windows - $finished; # 未完成
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% ";
    print "\n" if($process_count == $process_all_count);   
}

sub get_rna_extend_bed{
	my $input   	= shift @_;
	my $output   	= shift @_;
    my $up   		= shift @_;
    my $down 		= shift @_;
	print "Up_extend:$up\n";
	print "Down_extend:$down\n";
	
	open IN, $input;
	open OUT, qq{>$output};
	
	# my $line1 = <IN>;
	# $line1 =~ s/[\r\n]//g;
	# my @heads = split /\t/, $line1;
	# print OUT (join "\t", @heads) . "\n"; 
	readline IN;
	while (<IN>){
		$_=~s/[\r\n]//g;
		#next if /chr/;
		my @arr = split /\t/;
		my $id = shift @arr;
		my $chr = shift @arr;
		#   $chr =~ s/chr//g;
		my $start = shift @arr;
		$start = $start - $up;
		if($start < 0){$start = 0};
		my $end = shift @arr;
		$end = $end + $down;
		#next if($start < 0);
		next if($end < 0);
		
		print OUT "$chr\t$start\t$end\t$id\n";
	}
	close IN;
	close OUT;	
}


sub help
{
    my $info = "
Program: 表达量与甲基化相关性分析
         只对处于表达量定义染色体区域内的甲基化位点做相关性分析
Version: 2019-08-14
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
    必填：
        --rna/-r          表达量矩阵，每一行为基因、转录本，每一列为样本，第一列为ID,第2-4列为染色体、起始、终止位置。示例：/home/wangly/scripts/correlation/example/gene.xls
        --meth/-m         甲基化矩阵，每一行为位点ID，每一列为样本，第一列为ID,第2-4列为染色体、起始、终止位置。示例：/home/wangly/scripts/correlation/example/methyl.xls
        --sample/-s       样本名文件，包含一列，没有表头，样本名务必与rna/meth一致
        --output_dir/-o   结果输出目录。流程会自动创建结果目录。 
                          注意：最终结果只会保留相关系数为“负”的结果

    选填：
        --corrlation-method/-c   相关分析算法，pearson/spearman, 默认： pearson
        --thread/-t        并行线程数，默认：10
        --help/-h          查看帮助文档

    建议：
        分析前，先把RNA的起始、终止位置向两侧扩展500K，然后再分析，即：分析基因上下游500K范围内的甲基化与基因表达量的关系
        分析后，需要手工把excel表格中的RNA起始、终止位置调节回来    
    \n";
    return $info;

}