# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Parallel::ForkManager;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $sra_prefetch   = "/home/genesky/software/sra_tools/2.11.3/bin/prefetch";
my $sra_fastq_dump = "/home/genesky/software/sra_tools/2.11.3/bin/fastq-dump";
my $ascp           = "/home/xudl/.aspera/connect/bin/ascp|/home/xudl/.aspera/connect/etc/asperaweb_id_dsa.openssh";

# 检测 -> 脚本输入
my ($input_sra_code, $output_dir, $convert_fastq, $if_help);
GetOptions(
    "input_sra_code|i=s"   => \$input_sra_code,
    "output_dir|o=s"       => \$output_dir,
    "convert_fastq|c:s"    => \$convert_fastq,
    "help|h"               => \$if_help,
);
die help() if(defined $if_help or (not defined $input_sra_code or not defined $output_dir));
###################################################################### 主程序
set_ncbi_config();


# (1) 下载
print "Download sra file \n";
mkdir $output_dir if(not -e $output_dir);
system("$sra_prefetch --option-file $input_sra_code   -O $output_dir --max-size 2000000000");

# (2) 转换fastq
exit if(not defined $convert_fastq);
print "covert sra file to fastq \n";

# 读取sra编号
my @sra_codes;
open CODE, $input_sra_code;
while(<CODE>)
{
	$_=~s/[\s]//g;
	next if($_!~/\w/);
	my $sra_code = $_;
	push @sra_codes, $sra_code;
	# 
}
close CODE;

# 并行运行进行转换
my $fastq_dir = "$output_dir/fastq";
mkdir $fastq_dir if(not -e $fastq_dir);
my $pm = Parallel::ForkManager->new(10);
   $pm->run_on_start(sub{my ($pid, $sra_code) = @_; process_bar_array($sra_code, \@sra_codes)});# 进度条
foreach my $sra_code(@sra_codes)
{   
    $pm->start($sra_code) and next;
    system("$sra_fastq_dump  $output_dir/$sra_code/$sra_code.sra -O $fastq_dir/ --split-files --gzip --defline-qual '+' --defline-seq '\@\$ac-\$si/\$ri' ");
    $pm->finish;          
}
$pm->wait_all_children;


###################################################################### 子函数
sub set_ncbi_config{
    my $user = $ENV{'USER'};
    my $ncbi_dir = "/home/$user/.ncbi";
    my $ncbi_set = "$ncbi_dir/user-settings.mkfg";
    return if(-e $ncbi_set);

    mkdir $ncbi_dir if(not -e $ncbi_dir);

    open SET, ">$ncbi_set";
    print SET "## auto-generated configuration file - DO NOT EDIT ##\n\n";
    print SET "/LIBS/GUID = \"2449952d-ecd1-44bc-9fbc-f68e8e62fb29\"\n";
    print SET "/config/default = \"false\"\n";
    print SET "/repository/user/ad/public/apps/file/volumes/flatAd = \".\"\n";
    print SET "/repository/user/ad/public/apps/refseq/volumes/refseqAd = \".\"\n";
    print SET "/repository/user/ad/public/apps/sra/volumes/sraAd = \".\"\n";
    print SET "/repository/user/ad/public/apps/sraPileup/volumes/ad = \".\"\n";
    print SET "/repository/user/ad/public/apps/sraRealign/volumes/ad = \".\"\n";
    print SET "/repository/user/ad/public/root = \".\"\n";
    print SET "/repository/user/default-path = \"/home/$user/ncbi\"\n";
    close SET;
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

sub help{
    my $info = "
Program: sra download， 根据sra编号下载数据
Version: 2019-05-16
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_sra_code/-i sra编号文件。注：只有一列数据，为SRA编号，没有表头
         --output_dir/-o     结果输出目录
         --convert_fastq/-c  是否转换为fastq，无需参数，只需声明。
         --help/-h           查看帮助文档

题外话：
    sra数据转fastq是非常慢的，特别是全基因组，可能需要一周的时间。
    故，可以的话，推荐直接从EBI-ENA上下载fastq数据。
    EBI-ENA与NCBI的数据类似，但是保留的是fastq数据，且双方有数据交换，即：sra文件对应的fastq文件，在EBI-ENA上都可以找到。
    例如：
        SRR8670768编号文件，对应的EBI-ENA展示页面为：https://www.ebi.ac.uk/ena/data/view/SRR8670768，可以在页面下方找到对应的fastq的ftp地址
        然后使用aspera进行下载。
        个人开发一键化下载工具：/home/pub/bin/NGS/chip/GATK4/tools/personal/ebi_ena_fastq_download.pl

    \n";
    return $info;
}

