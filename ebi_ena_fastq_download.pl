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
my $ascp           = "/home/xudl/.aspera/connect/bin/ascp";
my $ascp_openssh   = "/home/xudl/.aspera/connect/etc/asperaweb_id_dsa.openssh";

# 检测 -> 脚本输入
my ($address_file, $output_dir, $if_help);
GetOptions(
    "address_file|a=s"     => \$address_file,
    "output_dir|o=s"       => \$output_dir,
    "help|h"               => \$if_help,
);
die help() if(defined $if_help or (not defined $address_file or not defined $output_dir));
###################################################################### 主程序

# (1) 获取下载地址
my @download_addresses;
open ADDRESS, $address_file;
while(<ADDRESS>)
{
    $_=~s/[\r\n\s]//g;
    next if($_!~/\w/);
    die "[ERROR] address error : $_\n file download address has to be start with ftp://ftp.sra.ebi.ac.uk\n forexample : ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR867/008/SRR8670768/SRR8670768_1.fastq.gz\n" if($_ !~ /^ftp\:\/\/ftp\.sra\.ebi\.ac\.uk/);
    push @download_addresses, $_;
}
close ADDRESS;


# (2) 下载
# 因为是使用aspera进行了加速，故不需要进行文件并行了，一个一个下载即可
print "Download file \n";
mkdir $output_dir if(not -e $output_dir);
foreach my $download_address(@download_addresses)
{
    process_bar_array($download_address, \@download_addresses);
    print "download : $download_address\n";

    $download_address =~ s/^ftp\:\/\/ftp\.sra\.ebi\.ac\.uk//; # 去掉前缀
    system("$ascp -v -QT -l 400m -P33001 -k1 -i $ascp_openssh era-fasp\@fasp.sra.ebi.ac.uk:$download_address $output_dir"); # EBI-ENA 的ascp下载地址、用户名是固定的
}


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
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% \n";
}

sub help{
    my $info = "
Program: EBI-ENA fastq download， 从EBI-ENA上下载NCBI原始fastq数据
Version: 2019-06-27
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --address_file/-a   EBI-ENA中fastq文件地址列表文件，每一行只能有一个文件。
                             地址示例：ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR867/008/SRR8670768/SRR8670768_1.fastq.gz
         --output_dir/-o     结果输出目录
         --help/-h           查看帮助文档
    \n";
    return $info;
}

