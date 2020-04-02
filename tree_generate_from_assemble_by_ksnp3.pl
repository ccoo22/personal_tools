# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Cwd 'abs_path';

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $ksnp3_dir_default = "/home/genesky/software/ksnp/3.1/kSNP3";
my $rscript           = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $plot_tree         = "/home/pub/bin/NGS/chip/GATK4/tools/personal/plot_tree.r";  # tree绘图工具

# 检测 -> 脚本输入
my ($genome_dir, $output_dir, $ksnp3_dir, $thread, $keep_tmp, $if_help);
GetOptions(
    "genome_dir|g=s"   => \$genome_dir,
    "output_dir|o=s"   => \$output_dir,
    "ksnp3_dir|k=s"    => \$ksnp3_dir,
    "thread|t=s"       => \$thread,
    "keep_tmp=s"       => \$keep_tmp,
    "help|h"           => \$if_help,
);
$ksnp3_dir = $ksnp3_dir_default if(not defined $ksnp3_dir);
$thread    = 10 if(not defined $thread);

$ENV{'PATH'} = "$ksnp3_dir:$ENV{'PATH'}";  # 环境加载

die help() if(defined $if_help or (not defined $genome_dir or not defined $output_dir or not defined $ksnp3_dir));

# 获取绝对路径，方便后面chdir处理
$genome_dir = abs_path($genome_dir);
$output_dir = abs_path($output_dir);
 
# 软件包下4个关键的软件检查
my $make_ksnp3_infile = "$ksnp3_dir/MakeKSNP3infile";
my $make_fasta        = "$ksnp3_dir/MakeFasta";
my $k_chooser         = "$ksnp3_dir/Kchooser";
my $ksnp3             = "$ksnp3_dir/kSNP3";
die "ksnp3 软件丢失，请仔细检查 : \n$make_ksnp3_infile \n$make_fasta \n$k_chooser \n$ksnp3\n" if(is_file_ok($make_ksnp3_infile, $make_fasta, $k_chooser, $ksnp3) == 0);
 
###################################################################### 主程序
my $tmp_dir = "$output_dir/tmp";
make_dir($output_dir, $tmp_dir);
 
##########
# (1) 创建ksnp3 基因组输入文件列表
##########
print "(1) make ksnp3 input\n";
chdir($genome_dir);
my $inlist = "$tmp_dir/inlist";
my $group_file = "$tmp_dir/group.txt";
system("$make_ksnp3_infile ./ $inlist A");  # 输入路径必须是相对路径，否则inlist不对，是ksnp3脚本的特性
system("awk '{print \$2\"\\tALL\"}' $inlist|sed '1i\\sample\\tgroup' > $group_file");  # 生成plot_tree.r需要的分组文件
chdir(PWD);  # 返回到当前目录

##########
# (2) 所有样本的fasta合并
##########
print "(2) combine all sample fasta together\n";
my $fastainput = "$tmp_dir/fastainput";
system("$make_fasta $inlist $fastainput");
 
##########
# (3) 确定最佳kmer长度
##########
print "(3) find the best kmer\n";
chdir($tmp_dir);
my $k_chooser_report = "$tmp_dir/Kchooser.report";
system("$k_chooser $fastainput");
my $k_best = get_best_k($k_chooser_report);
chdir(PWD);

##########
# (4) tree 分析
##########
print "(4) ksnp3 start \n";
my $log = "$output_dir/log.txt";
system("$ksnp3 -in $inlist -outdir $tmp_dir -k $k_best -ML -NJ -vcf -CPU $thread -core -min_frac 0.5 2>&1 |tee $log");

##########
# (5) tree 绘图，结果汇总
##########
print "(5) copy file and plot tree\n";

my @files = ("tree.core.tre",  # 简约法构建树FastTreeMP。基于核心SNP，核心SNP指的是在所有样本中都有基因型的SNP位点
             "tree.majority0.5.tre",  # 简约法构建树FastTreeMP。基于SNP分型比例在所有样本中 > majority fraction的SNP，位点比core要多一些
             "tree.parsimony.tre",  # 简约法构建树FastTreeMP。基于所有的SNP位点
             "tree.ML.tre",  # 最大释然法构建进化树FastTreeMP。基于所有的SNP位点
             "tree.NJ.tre",  # 邻接法构建进化树。基于所有的SNP位点
             "core_SNPs_matrix.fasta",  # 核心snp构成的多重比对序列
             "SNPs_in_majority0.5_matrix.fasta",  # majority 0.5 的SNP构成的多重比对序列
             "SNPs_all_matrix.fasta",  # 所有SNP构成的多重比对序列
             "SNPs_all"  # 所有SNP信息，kmer
               );
foreach my $file(@files)
{
    system("cp $tmp_dir/$file $output_dir/$file");  # 文件拷贝
    system("$rscript $plot_tree --tree_file $tmp_dir/$file      --group_file $group_file --output $output_dir/$file.pdf") if($file =~ /\.tre$/);  # 绘制进化树
}

print "Result: $output_dir\n";
print "    tree.core.tre                       简约法构建树FastTreeMP。基于核心SNP，核心SNP指的是在所有样本中都有基因型的SNP位点\n";
print "    tree.majority0.5.tre                简约法构建树FastTreeMP。基于SNP分型比例在所有样本中 > majority fraction的SNP，位点比core要多一些\n";
print "    tree.parsimony.tre                  简约法构建树FastTreeMP。基于所有的SNP位点\n";
print "    tree.ML.tre                         最大释然法构建进化树FastTreeMP。基于所有的SNP位点\n";
print "    tree.NJ.tre                         邻接法构建进化树。基于所有的SNP位点\n";
print "    core_SNPs_matrix.fasta              核心snp构成的多重比对序列\n";
print "    SNPs_in_majority0.5_matrix.fasta    majority 0.5 的SNP构成的多重比对序列\n";
print "    SNPs_all_matrix.fasta               所有SNP构成的多重比对序列\n";
print "    SNPs_all                            所有SNP位点，kmer形式展示\n";



system("rm -r $tmp_dir") if(not defined $keep_tmp);

###################################################################### 子函数

sub get_best_k{
    my $k_chooser_report = shift @_;
    my $k_best;
    open REPORT, $k_chooser_report;
    while(<REPORT>)
    {   
        ($k_best) = $_=~/(\d+)/ if($_=~/The optimum value of K is/)
    }
    close REPORT;
    return $k_best;
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

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}

# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}

sub help{
    my $info = "
Program: 基于细菌的基因组组装结果（可以是草图）， 制作进化树
Version: 2020-04-02
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [必填]
         --genome_dir/-g        基因组fasta文件目录。目录下保存了所有需要放在一起做进化树的组装fasta文件，一个样本一个文件，文件名格式：sample.fa
                                注意：sample名称里不能含有字符'.'等奇怪的符号。名字越简单越好，防止绘图时字符串过长。
         --output_dir/-o        结果输出目录

        [选填]
         --ksnp3_dir/-k         ksnp3软件bin路径，默认： /home/genesky/software/ksnp/3.1/kSNP3
         --thread/-t            并行线程数量，默认： 10
         --help/-h              查看帮助文档
         --keep_tmp             是否删除临时目录output_dir/tmp，默认删除

         注意：如果需要重新运行，请务必把上一次的结果删掉后再运行，否则会自动跳过一些步骤
    \n";
    return $info;
}

