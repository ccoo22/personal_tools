# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Parallel::ForkManager;
$|=1;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 检测 -> 脚本输入
my ($input_file, $output_dir, $execute, $threshold, $if_help);
GetOptions(
    "input_file|i=s" => \$input_file,
    "output_dir|o=s" => \$output_dir,
    "execute|e!"     => \$execute,
    "threshold|t=s"  => \$threshold,
    "help|h"         => \$if_help,
);
die help() if(defined $if_help or (not defined $input_file or not defined $output_dir));
###################################################################### 主程序
$threshold = 5 if(not defined $threshold or $threshold !~ /^\d+$/);  # 并行数量

############
# 1 创建输出目录
############
mkdir $output_dir if(not -e $output_dir);

############
# 2 开始处理输入文件
############
my @bash_files;
open FILE, $input_file;
my $line1 = <FILE>;
   $line1 =~ s/[\r\n]//g;
my @heads = split /\t/, $line1;
my $row = 1;
while(<FILE>)
{   
    $row++;
    $_ =~ s/[\r\n]//g;
    next if($_!~/\w/);
    my @datas = split /\t/, $_;

    # 数据存入哈希，便于后续根据表头提取信息
    my %hashTmp;
    foreach my $col(0..$#heads)
    {
        my $value = exists $datas[$col] ? $datas[$col] : "";
        $hashTmp{$heads[$col]} = $value;
    }
    
    # 开始文件检测 
    my $new_name = $hashTmp{'new_name'};
    my $r1_file_list = $hashTmp{'R1_files'};
    my $r2_file_list = $hashTmp{'R2_files'};
    $new_name =~ s/^\s+//;
    $new_name =~ s/\s+$//;
    die "[Error] 在行$row 中，new_name  表头缺失\n" if(not defined $new_name or $new_name !~ /\w/);
    die "[Error] 在行$row R1_files     表头缺失\n" if(not defined $r1_file_list or $r1_file_list !~ /\w/);
    die "[Error] 在行$row R2_files     表头缺失\n" if(not defined $r2_file_list or $r2_file_list !~ /\w/);

    my @r1_files = split /,/, $r1_file_list;
    my @r2_files = split /,/, $r2_file_list;
    die "[Error] 在行$row 中，R1_files与R2_files中的文件数量不一致" if(scalar(@r1_files) != scalar(@r2_files));

    my ($file_condition, $error_file)= is_file_ok((@r1_files, @r2_files));  # 文件是否丢失
    die "[Error] 在行$row 中，部分fastq文件丢失：$error_file\n" if($file_condition == 0);

    # bash文件生成
    my $bash_file = "$output_dir/$new_name.sh";

    open BASH, ">$bash_file";
    if(scalar(@r1_files) > 1) # 2个以上的文件，cat合并
    {
        print BASH "cat @r1_files > $output_dir/$new_name\_R1.fastq.gz\n";
        print BASH "cat @r2_files > $output_dir/$new_name\_R2.fastq.gz\n";
        print BASH "sleep 10s\n";
        print BASH "rm @r1_files \n";
        print BASH "rm @r2_files \n";
    }
    else # 只有1个文件，直接mv
    {
        print BASH "mv @r1_files $output_dir/$new_name\_R1.fastq.gz\n";
        print BASH "mv @r2_files $output_dir/$new_name\_R2.fastq.gz\n";      
    }

    close BASH;
    system("chmod 755 $bash_file");
    
    push @bash_files, $bash_file;

}
close FILE;
print "bash 文件已生成\n";
exit 0 if(not defined $execute);

############
# 3 执行合并
############
print "开始合并文件\n";
my $pm = Parallel::ForkManager->new($threshold);
   $pm->run_on_start(sub{my ($pid, $bash_file) = @_; process_bar_array($bash_file, \@bash_files)});# 进度条
foreach my $bash_file(@bash_files)
{
    $pm->start($bash_file) and next;
    system(" bash $bash_file ");
    $pm->finish;    
}
$pm->wait_all_children;




###################################################################### 子函数

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    my $error_file = "";
    foreach my $file(@files)
    {
        if (not -e $file or not -f $file or -s $file == 0)
        {
            $isOK = 0;
            $error_file .= "$file,";
        }
    }    
    return ($isOK, $error_file);
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
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

sub help{
    my $info = "
Program: 根据get_dir_fastq.pl软件提供的R1_files、R2_files信息，以及new_name信息，对分lane的数据做合并处理。
         软件会先生成bash文件，作为记录与检查。
         如果确认无误，再执行合并。
Version: 2019-12-26
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [必填]
         --input_file/-i     数据信息文件，必须包含表头：new_name,R1_files,R2_files
                             new_name: 希望生成的文件的前缀。所有结果会被命名为new_name_R1.fastq.gz,new_name_R2.fastq.gz
                             R1_files: R1文件绝对路径列表，多个文件之间用逗号分隔，不允许空字符，文件配对顺序务必与R2_files一致，否则合并出错
                             R2_files: R2文件绝对路径列表，多个文件之间用逗号分隔，不允许空字符，文件配对顺序务必与R1_files一致，否则合并出错

         --output_dir/-o     输出目录
        [选填]
         --execute/-e        执行bash合并命令 默认：不执行
         --threshold/-t      并行样本数量 默认：5

        --help/-h           查看帮助文档
    \n";
    return $info;
}

