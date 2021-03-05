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


# 检测 -> 脚本输入
my ($index_file, $output_dir, $undetermined_dir, $index_count, $if_help);
GetOptions(
    "index_file|i=s"         => \$index_file,
    "output_dir|o=s"         => \$output_dir,
    "undetermined_dir|u=s"   => \$undetermined_dir,
    "index_count|c=s"        => \$index_count,
    "help|h"                 => \$if_help,
);
die help() if(defined $if_help or (not defined $index_file or not defined $output_dir or not defined $undetermined_dir or not defined $index_count));
die "index_count can only be S/D\n" if($index_count ne 'S' and $index_count ne 'D');
###################################################################### 主程序

# (1) 输入fastq
print "\n1: check Undetermined fastq \n";
my $R1_raw = "$undetermined_dir/Undetermined_S0_R1_001.fastq.gz";
my $R2_raw = "$undetermined_dir/Undetermined_S0_R2_001.fastq.gz";
die "Lost $R1_raw, $R2_raw\n" if(not -e $R1_raw or not -e $R2_raw);
print "    R1 : $R1_raw\n";
print "    R2 : $R2_raw\n";


# (2) 读入index信息，创建句柄
print "\n2: creat ouput fastq handle \n";
my %hashHandle = cread_file_handle($index_file, $index_count, $output_dir);

#（3）提取
print "\n3: start output\n";
my $pm = Parallel::ForkManager->new(4);
foreach my $strand(("R1", "R2"))
{   
    $pm->start() and next;

    my $fastq_input = $strand eq "R1" ? $R1_raw : $R2_raw;
    abstract_fastq(\%hashHandle, $index_count, $strand, $fastq_input);

    $pm->finish;          
}
$pm->wait_all_children;
 
# (4) 关闭句柄
print "\n4: stop ouput fastq handle \n";
foreach my $index(keys %hashHandle)
{
    foreach my $strand(keys %{$hashHandle{$index}})
    {
        close $hashHandle{$index}{$strand};
    }
}



###################################################################### 子函数

# 提取fastq
sub abstract_fastq{
    my $hashHandle  = shift @_;
    my $index_count = shift @_;
    my $strand      = shift @_;
    my $fastq_input = shift @_;

    open FASTQ, "gzip -cd $fastq_input|";
    while(my $line1 = <FASTQ>)
    {
        my $line2 = <FASTQ>;
        my $line3 = <FASTQ>;
        my $line4 = <FASTQ>;
        my $head  = $line1;
           $head =~ s/[\r\n]//g;

        # 获取index信息
        my ($reads_name, $reads_anno) = split /\s+/, $head;
        my ($anno1, $anno2, $anno3, $anno4) = split /:/, $reads_anno; # 1:N:0:NTGAGTCC+NTACCGTG  ; 1:N:0:NATAAGAT
        my ($index_i7, $index_i5) = split /\+/, $anno4;

        my $index = $index_i7;
        if($index_count eq 'D')# 双index识别
        {
            $index = "$index_i7+$index_i5";
            die "[Error]  Undetermined数据是使用单index进行拆分的，只有 I7序列。你不能用双index模式从该文件中提取数据\n" if(not defined $index_i5);
        }

        # 输出到文件
        if(exists($hashHandle->{$index}))
        {
           my $file_handle = $hashHandle->{$index}{$strand}; # 获取文件句柄
           print $file_handle "$line1$line2$line3$line4";
        }
    }
    close FASTQ;    
}

# 创建文件句柄
sub cread_file_handle{
    my $index_file  = shift @_;
    my $index_count = shift @_;
    my $output_dir  = shift @_;

    my $sample_count = 0;
    my %hashHandle;
    open INDEX, $index_file;
    <INDEX>; 
    while(<INDEX>){
        $_ =~ s/[\r\n]//g;
        next if($_!~/\w/);
        my ($sample_name, $index_id_i7, $index_i7, $index_id_i5, $index_i5, $project, $tmp) = split /\t/, $_, 7;
        my $index = $index_i7;
           $index = "$index_i7+$index_i5" if($index_count eq 'D'); # 双index识别
 
        $sample_count++;
        # 打开文件句柄，放入哈希
        open $hashHandle{$index}{"R1"}, "| gzip > $output_dir/$sample_name\_R1.fastq.gz";
        open $hashHandle{$index}{"R2"}, "| gzip > $output_dir/$sample_name\_R2.fastq.gz";
    }
    close INDEX;
    print "    Process handle: $sample_count * 2\n"; 
    return %hashHandle;
}


sub help{
    my $info = "
Program: 根据index信息，从bcl2fastq的Undetermined文件里提取样本数据
Version: 2019-05-30
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --index_file/-i          输入index文件。任意表头（脚本会忽略第一行），
                                  双index格式示例（至少5列，重点看1、3、5列）：HXL_A1   G385    CATTATAC    508 GTCAGTAC    18AATH041_SZ16R0110A    18SZ0069
                                  单index格式示例（至少3列，重点看1、3列）   ：Ctrl13-Th-5w   G-2 CCGTGTCC    19B0318B-ATAC   ATAC
                                  例如：/home/pub/bin/NGS/chip/GATK4/tools/personal/bcl2fastq_tools/bcl2fastq_ref/infoI7.txt
         --output_dir/-o          输出目录
         --undetermined_dir/-u    Undetermined文件所在目录，例如：/home/lch/work/20190520_Hiseq_SH19GSL083_087-090/Raw_data
         --index_count/-c         index识别数目，S/D。S = Single，D = Double，即根据单/双index识别样本         
                                  注意：如果Undetermined是用单index拆分的数据，你的index_count只能是 Single模式。
         --help/-h                查看帮助文档
    \n";
    return $info;
}

