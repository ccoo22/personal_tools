use strict;
use warnings;
use Getopt::Long;
$|=1;

my  ($input_file,
     $output_file,
     $help,
);

GetOptions(
     'input|i=s'       => \$input_file,
     'output|o=s'      => \$output_file,
     'help|h!'         => \$help,
);

if ($help or not $input_file or not $output_file) {
	usage();
	exit();
}


open REF , $input_file or die "[ERROR] can not open $input_file";
# 统计总行数
my $totle_lines;
while (<REF>){
     $totle_lines = $.;
}
close REF;

my %allcds = ();
my $index  = 1;

print "Reading refGene...\n";

open REF , $input_file or die "[ERROR] can not open $input_file";
while (<REF>)
{
     $_ =~ s/[\r\n]//g;
     my @lines = split "\t", $_;
     my ($tmp1, $tran_name, $chrom, $strand, $tmp2, $tmp3, $cds_start, $cds_end, $tmp4, $exon_starts, $exon_ends, $tmp5, $gene_name) = split "\t", $_;

     # 显示refGene文件的读取进度
     process_bar_array($., $totle_lines);

     # 跳过非常规染色体以及无CDS区域的数据
     next if ($chrom =~ /_/);
     next if ($cds_start == $cds_end);

     # 提取exon区域与CDS起止位点信息
     my @exon_starts = split ",", $exon_starts;
     my @exon_ends   = split ",", $exon_ends;
     my @exon = (@exon_starts, @exon_ends, $cds_start, $cds_end);


     ######################################################################
     #
     #	根据cds起止位点与exon区域的相对位置关系，确定cds区域
     #
     ######################################################################


     @exon = sort {$a <=> $b} @exon;

     # 获取cds_start和cds_end在所有区域起止位点中的位置
     my ($cds_s_p) = grep {$exon[$_] == $cds_start} 0..$#exon;
     my ($cds_e_p) = grep {$exon[$_] == $cds_end} 0..$#exon;

     # 根据cds_end位置的奇偶情况，判断cds_end落在外显子区域或内含子区域，以此确定cds区域的终止位置，删除终止位置之后的位点
     $cds_e_p%2 == 0 ? splice (@exon, $cds_e_p+1) : splice (@exon, $cds_e_p);

     # 根据cds_start位置的奇偶情况，判断cds_start落在外显子区域或内含子区域，以此确定cds区域的起始位置，删除起始位置之前的位点
     @exon = $cds_s_p%2 == 0 ? splice (@exon, $cds_s_p+1) : splice (@exon, $cds_s_p);


     ######################################################################
     #
     #	获取片段来源信息，结合序列方向，确定cds片段对应的exon编号
     #
     ######################################################################


     my $form_info="$gene_name,$tran_name";

     # 确定总外显子数量，以及首个cds区域对应的外显子区域编号
     my $exon_num = @exon_starts;
     my $start_num = $cds_s_p%2 == 0 ? $cds_s_p/2 : ($cds_s_p-1)/2;

     foreach (0..($#exon-1)/2)
     {
          $allcds{$chrom}{$index}{'cdsstart'} = $exon[2*$_];
          $allcds{$chrom}{$index}{'cdsend'}   = $exon[2*$_+1];
          if ($strand eq "+")
          {
               $allcds{$chrom}{$index}{'info'} = "$form_info,exon".($start_num+1+$_);
          }
          else
          {
               $allcds{$chrom}{$index}{'info'} = "$form_info,exon".($exon_num-$start_num-$_);
          }
          $index++;
     }
}
close REF;

# 对染色体进行排序
my @chroms = keys %allcds;
@chroms = sort_chrom(\@chroms);

open BED , ">$output_file" or die "[ERROR] can not open $output_file";
print "Forming BED file...\n";
foreach (@chroms)
{
     # 调用unique子函数进行去重复处理
     my %chrcds = unique($allcds{$_}) if (exists $allcds{$_});
     # 按照bed文件格式进行输出
     print_bed(\%chrcds, $_);
}
close BED;
print "All Processes Completed!\n";


sub usage{
     my $help =<<EOF;
Program: cds2bed格式转换
Version: 2019-08-08
Contact: 312 余晨阳

Usage: perl $0  [ Options ]

[ Options ]
-input|-i  (Required)        refGene
-output|-o (Required)        BEDoutput
-help|-h                     print help message
EOF
     print $help;
     exit;
}


sub process_bar_array{
    # 当前对象的位置
    my $process_count   = shift @_;
    # 需要分析的数据总量
    my $process_all_count = shift @_;
    # 当前进度
    my $process_perc = sprintf "%0.2f", 100 * $process_count/$process_all_count;
    # 进度条总长度
    my $windows = 100;
    # 已完成部分长度
    my $finished = $windows * int($process_perc) / 100;
    # 待完成部分长度
    my $unfinished = $windows - $finished;
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% ";
    print "\nRead All refGene!\n" if($process_count == $process_all_count);   
}


sub sort_chrom{
     my $chroms = shift @_;
     my @chroms = @$chroms;

     # 判断是否有chr前缀，并去除前缀
     my $flag = $chroms[0] =~ /^chr/;
     map {$_ =~ s/^chr// } @chroms if $flag;

     # 对纯数字编号的染色体按编号大小进行排序
     my @chroms_num = grep {$_ !~ /\D/} @chroms;
     @chroms_num = sort {$a <=> $b} @chroms_num;

     # 对非纯数字编号的染色体按字符串方式进行排序
     my @chroms_str = grep {$_ =~ /\D/} @chroms;
     @chroms_str = sort  @chroms_str;

     @chroms = (@chroms_num, @chroms_str);
     # 加回前缀
     map {$_ = "chr".$_ } @chroms if $flag;

     return @chroms;
} 


sub unique{
     my $area = shift @_;
     my %area = %$area;
     my @dels  = ();

     # 按CDS起始位置排序 
     my @arr  = sort {$area{$a}{'cdsstart'} <=> $area{$b}{'cdsstart'}} keys %area;
     foreach my $order (1..$#arr)
     {
          # 存在区域重叠
          if ($area{$arr[$order]}{'cdsstart'} <= $area{$arr[$order-1]}{'cdsend'})
          {
               # 区域之间呈包含关系，第n个区域被包含在第n-1个区域内
               if ($area{$arr[$order]}{'cdsend'} <= $area{$arr[$order-1]}{'cdsend'})
               {
                    # 将第n个区域更新为第n-1个区域
                    $area{$arr[$order]}{'cdsstart'} =$area{$arr[$order-1]}{'cdsstart'};
                    $area{$arr[$order]}{'cdsend'} =$area{$arr[$order-1]}{'cdsend'};
               }
               else
               {
                    # 将第n个区域更新为两个区域的并集区域
                    $area{$arr[$order]}{'cdsstart'} =$area{$arr[$order-1]}{'cdsstart'};
               }
               # 合并第n-1个区域与第n个区域的来源信息
               $area{$arr[$order]}{'info'} = "$area{$arr[$order-1]}{'info'}|$area{$arr[$order]}{'info'}";
               # 将第n-1个区域记入CDS重复区域，等待删除
               push @dels, $order-1;
          }
     }

     # 删除CDS重复区域
     foreach my $del (@dels)
     {
          delete $area{$arr[$del]};
     }
     return %area;
}


sub print_bed{
     my $area  = shift @_;
     my $chrom = shift @_;
     my %area  = %$area;

     foreach my $index (sort {$area{$a}{'cdsstart'} <=> $area{$b}{'cdsstart'}} keys %area)
     {
          my $len = $area{$index}{'cdsend'}-$area{$index}{'cdsstart'};
          # 计算片段长度
          $area{$index}{'cdsstart'}++;
          # 输出bed文件
          print BED "$chrom\t$area{$index}{'cdsstart'}\t$area{$index}{'cdsend'}\t$len\t$area{$index}{'info'}\n";
     }
}
