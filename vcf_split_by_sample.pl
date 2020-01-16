# 导入 -> 系统 package
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
$|=1;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $bcftools = "/home/genesky/software/bcftools/1.9/bin/bcftools";

# 检测 -> 脚本输入
my ($input_vcf, $samples, $output_dir, $thread, $help);
GetOptions(
	'input_vcf|i=s'    => \$input_vcf,
	'samples|s=s'      => \$samples,
	'output_dir|o=s'   => \$output_dir,
	'thread|t=s'	   => \$thread,
	"help|h!"          => \$help,
);
die help() if (defined $help or not defined $input_vcf or not defined $samples or not defined $output_dir);
$thread   = 5   if (not defined $thread);

###################################################################### 主程序
mkdir $output_dir if(not -e $output_dir);

# (1) 获取样本
my @samples = split /,/, $samples;

# (2) 并行
my $pm = Parallel::ForkManager->new($thread);
   $pm->run_on_start(sub{my ($pid, $sample) = @_; process_bar_array($sample, \@samples)});# 进度条
foreach my $sample(@samples){
	my $pid = $pm->start($sample) and next;

	system("$bcftools view -s $sample $input_vcf | $bcftools filter --exclude 'GT=\"ref\"' -o $output_dir/$sample.vcf.gz -O z");  # 提取指定样本的所有突变
	system("$bcftools index $output_dir/$sample.vcf.gz");  # 建索引

	$pm->finish;
}
$pm->wait_all_children;

###################################################################### 子函数
 
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
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc%";
    print "\n" if($process_count == $process_all_count);   
}

sub help
{
my $info = "
Usage: perl $0  -i input_vcf -s samples -o output_dir
 Options: 必填
	--input_vcf/-i	sample.all.final.vcf.gz
	--samples/-s	需要的样本列表，示例： FHY_N,FHY_T,HM_N
	--output_dir/-o	输出目录, 流程自己可以创建
				
 Options: 可选
	--thread/-t	并行样本数量，默认: 5
	--help/-h	查看帮助文档
	\n";
	return $info;
}