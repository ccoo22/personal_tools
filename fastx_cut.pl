# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Cwd qw( abs_path );
use Parallel::ForkManager;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $fastx_trimmer = "/home/genesky/software/fastx_toolkit/0.0.14/bin/fastx_trimmer";

# 检测 -> 脚本输入
my ($input_dir, $output_dir, $first_base, $last_base, $trim_tail, $parallel, $if_help);
GetOptions(
    "input_dir|i=s"   => \$input_dir,
    "output_dir|o=s"  => \$output_dir,
    "first_base|f:i"  => \$first_base,
    "last_base|l:i"   => \$last_base,
    "trim_tail|t:i"   => \$trim_tail,
    "parallel|p:i"    => \$parallel,
    "help|h"          => \$if_help,
);
die help() if(defined $if_help or (not defined $input_dir or not defined $output_dir));
my $first_base_para = (defined $first_base) ? "-f $first_base" : "";
my $last_base_para  = (defined $last_base)  ? "-l $last_base"  : "";
my $trim_tail_para  = (defined $trim_tail)  ? "-t $trim_tail"  : "";
$parallel   = 1 if(not defined $parallel);
$input_dir  = Cwd::abs_path($input_dir);
$output_dir = Cwd::abs_path($output_dir);
die "[Error] input_dir can not be the same as output_dir\n" if($input_dir eq $output_dir);
###################################################################### 主程序

# (1) 获取需要处理的文件
my %hashFile = get_seq_file($input_dir);

# (2) 处理
# 提示
my @file_lists = sort keys %hashFile;
my $file_count = scalar(@file_lists);
print "\nProcess file[$file_count]:\t" . (join ",", @file_lists) . "\n";
print "\nProcess Parameter: \n    first_base = $first_base_para\n    first_base = $last_base_para\n    trim_tail = $trim_tail_para\n    parallel = $parallel\n    input_dir = $input_dir\n    output_dir = $output_dir\n";
is_continue();

# 开始
my $pm = Parallel::ForkManager->new($parallel);
$pm->run_on_start(sub{my ($pid, $file_name) = @_; process_bar_array($file_name, \@file_lists)});# 进度条
foreach my $file_name(@file_lists)
{
    $pm->start($file_name) and next;
    my $input_file = "$input_dir/$file_name";
    my $output_file = "$output_dir/$file_name";
    my $is_gz = ($file_name =~ /\.gz$/) ? 1 : 0; # 有无压缩
    print "trim $input_file -> $output_file\n";

    system("                  $fastx_trimmer $first_base_para $last_base_para $trim_tail_para -i $input_file -o $output_file") if($is_gz == 0);
    system("zcat $input_file| $fastx_trimmer $first_base_para $last_base_para $trim_tail_para -z             -o $output_file") if($is_gz == 1);
    
    
    $pm->finish;
}
$pm->wait_all_children; 




###################################################################### 子函数
# 获取文件列表
sub get_seq_file{
    my $input_dir = shift @_;
    my %hashFile;

    opendir DIR, "$input_dir/";
    foreach my $file_name(readdir DIR)
    {   
        next if($file_name !~ /\w/ or not -f "$input_dir/$file_name");
        next if($file_name !~ /\.fq$/ and $file_name !~ /\.fastq$/ and $file_name !~ /\.fq\.gz$/ and $file_name !~ /\.fastq\.gz$/ and $file_name !~ /\.fa$/ and $file_name !~ /\.fasta$/ and $file_name !~ /\.fa\.gz$/ and $file_name !~ /\.fasta\.gz$/);       
        $hashFile{$file_name}++;
    }
    closedir DIR;
    return %hashFile;
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
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% \n";
}

# 是否继续
sub is_continue{
    print "\n[Option] Confirm and Start?[y/n]";
    my $input = get_input();
    exit if($input ne "y" and $input ne 'Y');
}

# 从屏幕获得输入
sub get_input{
    my $input=<STDIN>;
    $input=~ s/[\r\n]//g;
    return $input;
}
sub help{
    my $info = "
Program: fastx_trimmer， 序列截取
Version: 2019-05-22
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_dir/-i    输入序列目录，例如：/home/pub/project/genetics/18P1209A00A。会处理该目录下的所有以 .fq、.fastq、.fq.gz、.fastq.gz、.fa、.fasta、.fa.gz、.fasta.gz为结尾的文件
                           注意：核苷酸序列暂不支持小写字符。分析软件基于fastx_trimmer
         --output_dir/-o   输出结果目录
         --first_base/-f   序列截取起始位置，默认：1
         --last_base/-l    序列截取终止位置，默认：全长序列
         --trim_tail/-t    砍掉尾部Nbp序列， 默认：0，不砍
                           注意：trim_tail声明后，first_base,last_base将失去作用
         --parallel/-p     并行文件数量，默认：1
         --help/-h         查看帮助文档
    \n";
    return $info;
}

