$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
bam转fastq
Version: v1.0 2020-09-23
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Parallel::ForkManager;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 变量
my $DEFAULT_THREAD          = 10;
my $DEFAULT_SOFT_SAMTOOLS   = "/home/genesky/software/samtools/1.10/samtools";
my $DEFAULT_SOFT_MERLIN     = "/home/genesky/software/merlin/1.1.2/executables";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input, $SOFT_SAMTOOLS, $thread, $if_help);
GetOptions(
	"input|i=s"           => \$input,

	"samtools=s"          => \$SOFT_SAMTOOLS,
	"thread=s"            => \$thread,
	"help|h"              => \$if_help,
);
die "
Options: 必填
        --input/-i                输入文件，有表头，三列数据： bam文件绝对路径、样本名、输出目录(脚本自动创建)
                                  输出文件会使用样本名进行命名
                                  生成的文件有：
                                  sample_R1.fastq.gz  
                                  sample_R1.fastq.gz 
                                  sample.sort_name.bam

Options: 可选
        --samtools                更改软件 samtools 版本 (default: '$DEFAULT_SOFT_SAMTOOLS')
        --thread                  并行运行样本数量 (default: '$DEFAULT_THREAD') 
        --help/-h                 查看帮助文档

\n" if (defined $if_help or not defined $input );

$SOFT_SAMTOOLS   = $DEFAULT_SOFT_SAMTOOLS if (not defined $SOFT_SAMTOOLS);
$thread          = $DEFAULT_THREAD if (not defined $thread);
  
 
###################################################################### 主程序
my %hashInfo = read_matrix($input, 0);
my $sample_index     = $hashInfo{'HEAD'}->[1];
my $output_dir_index = $hashInfo{'HEAD'}->[2];

my @bams = sort {$hashInfo{'DATA'}{$a}{'sort_value'}<=>$hashInfo{'DATA'}{$b}{'sort_value'}} keys %{$hashInfo{'DATA'}}; 
my $pm = Parallel::ForkManager->new($thread);
   $pm->run_on_start(sub{my ($pid, $bam) = @_; process_bar_array($bam, \@bams)});# 进度条
foreach my $bam(@bams)
{
    $pm->start($bam) and next;
    my $sample       = $hashInfo{'DATA'}{$bam}{$sample_index};
    my $output_dir   = $hashInfo{'DATA'}{$bam}{$output_dir_index};
    mkdir $output_dir if(not -e $output_dir);
   
    print "[process] $sample 按照名字排序\n";
    my $bam_sort_name = "$output_dir/$sample.sort_name.bam";
    system("$SOFT_SAMTOOLS sort -n -o  $bam_sort_name $bam");

    print "[process] $sample bam 转 fastq\n";
    my $r1 = "$output_dir/$sample\_R1.fastq.gz";
    my $r2 = "$output_dir/$sample\_R2.fastq.gz";
 
    system("$SOFT_SAMTOOLS fastq -1 $r1  -2 $r2 -0 /dev/null -s /dev/null -N  $bam_sort_name ");

    $pm->finish;    
}
$pm->wait_all_children;


###################################################################### 子程序

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


# 读取带有表头的矩阵
sub read_matrix{
    my $file     = shift @_;
    my @key_cols = @_;

    print "Read $file\n";
    my %hashMatrix;
    open FILE, $file;
    my $line1 = <FILE>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;
    $hashMatrix{'HEAD'} = \@heads;

    while(<FILE>)
    {
        $_ =~ s/[\r\n]//g;
        next if($_!~/\w/);
        my @datas = split /\t/, $_;
        my @key_values = map{ $datas[$_] } @key_cols;
        my $key_value  = join "\t", @key_values;
        foreach my $col(0..$#heads)
        {
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $hashMatrix{'DATA'}{$key_value}{$heads[$col]} = $value;
        }
        $hashMatrix{'DATA'}{$key_value}{'sort_value'} = $.;
    }
    close FILE;
    return %hashMatrix;
}
