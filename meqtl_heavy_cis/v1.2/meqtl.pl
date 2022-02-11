$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
QTL cis 分析
Version: v1.2 2020-07-30
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Excel::Writer::XLSX;
use Encode;
use Cwd 'abs_path';

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 变量
my $DEFAULT_CIS_DISTANCE   = "1e6";
my $DEFAULT_CIS_PVALUE     = "0.05";
my $DEFAULT_REPORT_PVALUE  = "0.001";
my $DEFAULT_ANNO_THREAD    = 10;
my $DEFAULT_BOXPLOT_TOP    = 300;  
my $DEFAULT_PLOT_Y_LAB     = "EXP";
my $DEFAULT_SOFT_RSCRIPT   = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $DEFAULT_SOFT_RLIB      = "/home/genesky/software/r/3.5.1/lib64/R/library";
my $DEFAULT_SOFT_MEQTL     = SCRIPT_DIR . "/meqtl_cis.r";
my $DEFAULT_SOFT_boxplot   = SCRIPT_DIR . "/meqtl_top_boxplot.r";
my $DEFAULT_READ_ME        = SCRIPT_DIR. "/readme.txt";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($snp_file, $exp_file, $cor_file, $sample_file, $output_dir, $cis_distance, $cis_pvalue, $report_pvalue, $snp_add_col, $exp_add_col, $anno_thread, $boxplot_top, $plot_y_lab, $SOFT_RSCRIPT, $SOFT_RLIB, $if_help);
GetOptions(
	"snp_file|s=s"     => \$snp_file,
	"exp_file|e=s"     => \$exp_file,
	"sample_file=s"    => \$sample_file,
	"output_dir|o=s"   => \$output_dir,

	"cor_file|c=s"     => \$cor_file,

	"cis_distance=s"   => \$cis_distance,
	"cis_pvalue=s"     => \$cis_pvalue,
	"report_pvalue=s"  => \$report_pvalue,

	"snp_add_col=s"  => \$snp_add_col,
	"exp_add_col=s"  => \$exp_add_col,

	"anno_thread=s"  => \$anno_thread,

	"boxplot_top=s"  => \$boxplot_top,
	"plot_y_lab=s"  => \$plot_y_lab,

    "rscript=s"      => \$SOFT_RSCRIPT,
    "rlib=s"        => \$SOFT_RLIB,

	"help|h"        => \$if_help,
);  
die "
Options: 必填

        --snp_file/-s     SNP分型文件，前三列数据固定，分别是：snp_id chr pos。后面的列可以是：注释列、样本列。
                          样本基因型格式要求： “A/T” 这种类型的数据，如果缺失，空白即可。注意不要放空格

        --exp_file/-e     表达量文件，前4列数据固定，分别是：id chr start end。后面的列可以是：注释列、样本列。
                          如果某个样本表达量缺失，空白即可。注意不要放空格。

        --sample_file     QTL分析要用到的样本列表，一列数据，有表头。务必保证snp_file/exp_file 含有这些样本。

        --output_dir/-o   结果输出目录       

Options: 可选
        --cor_file/-c     协变量矫正文件，qtl分析模型是线性回归，通过添加协变量，矫正环境因子。第一列是环境变量名称，第一行是样本名称。有表头。第二列开始是每一个样本的表型。
                          环境因子ID名称要保证只含有字母、数字、下划线。
        --cis_distance    cis分析时的距离限制 [默认： $DEFAULT_CIS_DISTANCE] 
        --cis_pvalue      cis分析时，pvalue过滤阈值， 大于该阈值的结果不会输出 [默认： $DEFAULT_CIS_PVALUE] 
        --report_pvalue   分析结果导入excel表格时，pvalue过滤阈值，大于该阈值的结果不会输出，因为结果可能会非常多，excel表格放不下 [默认： $DEFAULT_REPORT_PVALUE] 

        --snp_add_col     从snp_file中挑选指定的列数据，放入excel中，用于注释。 第一列的编号为0， 需要多列信息的话，用逗号分隔，例如： 1,2,3。 注意：snp_id chr pos默认会加入excel表格中，这里不必再指定。
                          另外，结果文件中，脚本会自动为这些添加的注释加上 ‘snp_’ 的前缀，防止冲突
        --exp_add_col     从exp_file中挑选指定的列数据，放入excel中，用于注释。 第一列的编号为0， 需要多列信息的话，用逗号分隔，例如： 1,2,3。 注意：id chr start end默认会加入excel表格中，这里不必再指定。
                          另外，结果文件中，脚本会自动为这些添加的注释加上 ‘gene_’ 的前缀，防止冲突

        --anno_thread     QTL原始结果添加注释时，使用的并行线程数量 [默认： $DEFAULT_ANNO_THREAD] 
                          QTL原始结果比较简陋，只有简单的pvalue等信息，所以我们自己对结果添加注释。但是pair结果较多，注释会很慢，引入了多线程功能。
        --boxplot_top     对最显著的前n个结果绘制boxplot图 [默认： $DEFAULT_BOXPLOT_TOP] 
        --plot_y_lab      boxplot绘图时，y轴名称设置 [默认： $DEFAULT_PLOT_Y_LAB]
        --rscript         修改软件 Rscript 路径 [默认： $DEFAULT_SOFT_RSCRIPT] 
        --rlib            修改软件 rlib 路径 [默认： $DEFAULT_SOFT_RLIB] 
        --help/-h            查看帮助文档
\n" if (defined $if_help or not defined $snp_file or not defined $exp_file or not defined $sample_file or not defined $output_dir);

$cis_distance  = $DEFAULT_CIS_DISTANCE  if(not defined $cis_distance);
$cis_pvalue    = $DEFAULT_CIS_PVALUE    if(not defined $cis_pvalue);
$report_pvalue = $DEFAULT_REPORT_PVALUE if(not defined $report_pvalue);
$snp_add_col   = "" if(not defined $snp_add_col);
$exp_add_col   = "" if(not defined $exp_add_col);

$anno_thread   = $DEFAULT_ANNO_THREAD if(not defined $anno_thread);
$boxplot_top   = $DEFAULT_BOXPLOT_TOP if(not defined $boxplot_top);
$plot_y_lab   = $DEFAULT_PLOT_Y_LAB if(not defined $plot_y_lab);

$SOFT_RSCRIPT = $DEFAULT_SOFT_RSCRIPT if(not defined $SOFT_RSCRIPT);
$SOFT_RLIB = $DEFAULT_SOFT_RLIB if(not defined $SOFT_RLIB);


# 切换绝对路径
$snp_file    = abs_path($snp_file);
$exp_file    = abs_path($exp_file);
$sample_file = abs_path($sample_file);
$output_dir  = abs_path($output_dir);
###################################################################### 初始化
my $DATA_TIME = `date +\"\%Y-\%m-\%d \%H:\%M.\%S\"`;
my $RUN_INFO = "
---
Command: perl ".abs_path($0)." $ARGV_INFO
---
Start time: $DATA_TIME
[SET] 参数 cis_distance : $cis_distance
[SET] 参数 cis_pvalue : $cis_pvalue
[SET] 参数 report_pvalue : $report_pvalue
[SET] 参数 snp_add_col : $snp_add_col
[SET] 参数 exp_add_col : $exp_add_col
[SET] 参数 anno_thread : $anno_thread
[SET] 参数 rscript : $SOFT_RSCRIPT
[SET] 参数 rlib : $SOFT_RLIB
";

print $RUN_INFO;

###################################################################### 主程序
my $tmp_dir = "$output_dir/tmp";
make_dir($output_dir);
make_dir($tmp_dir);


# (1) 检查样本是否异常
my %hashSample = read_matrix($sample_file, 0);
check_sample(\%hashSample, $snp_file, "SNP FILE");
check_sample(\%hashSample, $exp_file, "EXP FILE");
check_sample(\%hashSample, $cor_file, "COR FILE") if(defined $cor_file);

# (2) snp 数据准备
my $snp_cis_prefix = "$tmp_dir/snp";
prepare_snp($snp_file, \%hashSample, $snp_cis_prefix);

# (3) exp 数据准备
my $exp_cis_prefix = "$tmp_dir/exp";
prepare_exp($exp_file, \%hashSample, $exp_cis_prefix);

# (3.1) cor数据准备
my $cor_txt = "$tmp_dir/cor.txt";
my $cor_para = "--cor_file $cor_txt" if(defined $cor_file);
prepare_cor($cor_file, \%hashSample, $cor_txt) if(defined $cor_file);


# (4) QTL分析
print "[process] start QTL analysis\n";
my $result_cis = "$output_dir/result.cis.pvalue$cis_pvalue.txt";
my $snp_anno_para = ($snp_add_col eq "") ? "" : "--snp_anno $snp_cis_prefix.anno.txt";  # 两个额外注释文件
my $exp_anno_para = ($exp_add_col eq "") ? "" : "--exp_anno $exp_cis_prefix.anno.txt";
system("$SOFT_RSCRIPT $DEFAULT_SOFT_MEQTL --cis_p  $cis_pvalue --cisdist $cis_distance --thread $anno_thread --output $result_cis --exp_file  $exp_cis_prefix.raw.txt $cor_para --exp_loc_file $exp_cis_prefix.pos.txt --snv_file $snp_cis_prefix.code.txt --snv_loc_file $snp_cis_prefix.pos.txt  $snp_anno_para $exp_anno_para --rlib $SOFT_RLIB");

# (5) top300 boxplot绘图
print "[process] start boxplot top $boxplot_top \n";
my $top_info = "$tmp_dir/top_result.$boxplot_top.txt";
my $boxplot_pdf = "$output_dir/boxplot.top$boxplot_top.pdf";
$boxplot_top++;
system("head -n $boxplot_top $result_cis > $top_info");
$boxplot_top--;

system("$SOFT_RSCRIPT $DEFAULT_SOFT_boxplot --top_info $top_info --output $boxplot_pdf --exp_file  $exp_cis_prefix.raw.txt --snv_file $snp_cis_prefix.raw.txt --y_lab '$plot_y_lab' --rlib $SOFT_RLIB");


# (6) 结果输出
print "[process] excel output\n";
my $excel = "$output_dir/result.cis.pvalue$report_pvalue.xlsx";
text_to_excel($result_cis, $DEFAULT_READ_ME, $report_pvalue, $excel);

print "[result] 按照report_pvalue参数，写入excel表格的结果：    $excel\n";
print "[result] 按照cis_pvalue参数，生成的原始QTL的结果：       $result_cis\n";
print "[result] 以snp-gene距离为X轴，-log10(pvalue)为Y轴绘图： $output_dir/scatter.png\n";
print "[result] 选取每个exp_id最显著的pvalue绘制曼哈顿图：      $output_dir/gene.most_sig_pvalue.manhattan.jpg\n";
print "[result] 选取每个snp_id最显著的pvalue绘制曼哈顿图：      $output_dir/snp.most_sig_pvalue.manhattan.jpg\n";
print "[result] 选取最显著的 top $boxplot_top 绘制箱体图：     $boxplot_pdf\n";

############################################################ 子程序

sub text_to_excel{
	my $txt     = shift @_;
	my $read_me = shift @_;
    my $pvalue  = shift @_;
	my $xlsx    = shift @_;

    my $workbook = Excel::Writer::XLSX->new($xlsx);
    my %format   = format_run($workbook);
    my $sheet    = $workbook->add_worksheet("meQTL");
    my $row      = 0;
    my @titles   = ();
    open TXT, $txt;
    while(<TXT>){
        $_=~s/[\r\n]//g;
        my @datas = split/\t/, $_;

        if ($row==0) 
        {
            @titles = @datas;
            foreach my $col(0..$#datas) 
            {
                $sheet->write($row, $col ,$datas[$col], $format{'title'});
            }
            $row++;
            next;
        }

        # 过滤
        my %hashTmp;
        my @datas_full;
        foreach my $col(0..$#titles)
        {
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $hashTmp{$titles[$col]} = $value;
            push @datas_full, $value;
        }
        next if($hashTmp{'pvalue'} > $pvalue);  # 仅放入pvalue允许的结果

        # 输出
        foreach my $col(0..$#datas_full) 
        {
            $sheet->write($row, $col ,$datas_full[$col], $format{'normal'});
        }
        $row++;
        if($row >= 1048576)
        {
            print "[Warning] 结果数量超过 1048576， 无法全部放入excel表格\n";
            last;
        }
    }
    close TXT;
    ReadMe($workbook, \%format, $read_me);
}

sub ReadMe{
    my ($workbook,$format,$file) = @_;
    my $sheet = $workbook->add_worksheet("Read Me");
    $sheet->set_row(0, 65);
    $sheet->set_column('A:A', 35);
    $sheet->set_column('B:B', 25);
    $sheet->set_column('C:C', 110);
    my $row=0;

    open FILE,$file;
    while(<FILE>){
        $_=~ s/\^/\n/g;
        my @split_line=split /\t/,$_;
        my $col=0;
        foreach (@split_line) {
            my $text = decode("UTF-8",$_);
            if($row==0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme1'});}
            if($row==0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme2'});}
            if($row==0 and $col==2){$sheet->write($row,$col,$text,$format->{'readme2tmp'});}
            if($row>0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme3'});}
            if($row>0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme4'});}
            if($row>0 and $col==2){$sheet->write($row,$col,$text,$format->{'readme5'});}
            $col++;
        }
        $row++;
    }
    close FILE;
}

sub prepare_cor{
    my $cor_file       = shift @_;
    my $hashSample     = shift @_;
    my $cor_txt        = shift @_;
    my @samples        = sort {$hashSample{'DATA'}{$a}{'sort_value'} <=> $hashSample{'DATA'}{$b}{'sort_value'}}keys %{$hashSample{'DATA'}};

    # 
    print "[process] Prepare cor data\n";
    my %hashCOR = read_matrix($cor_file, 0);

    open COR, ">$cor_txt";
    print COR "id\t" . (join "\t", @samples) . "\n";
    foreach my $cor_id(sort {$hashCOR{'DATA'}{$a}{'sort_value'} <=> $hashCOR{'DATA'}{$a}{'sort_value'}} keys %{$hashCOR{'DATA'}})
    {   
        my @cor_values = map{ $hashCOR{'DATA'}{$cor_id}{$_} } @samples;
        print COR "$cor_id\t" . (join "\t", @cor_values) . "\n";
    }
    close COR;
}

sub prepare_exp{
    my $exp_file       = shift @_;
    my $hashSample     = shift @_;
    my $exp_cis_prefix = shift @_;
    my @samples        = sort {$hashSample{'DATA'}{$a}{'sort_value'} <=> $hashSample{'DATA'}{$b}{'sort_value'}}keys %{$hashSample{'DATA'}};

    # 
    print "[process] Prepare exp data\n";
    my %hashEXP = read_matrix($exp_file, 0);

    # 确定要添加的exp注释
    my @anno_heads = ();
    my @anno_heads_surffix = ();  # 添加前缀，防止与exp冲突
    if($exp_add_col ne "")
    {
        @anno_heads = map{ $hashEXP{'HEAD'}->[$_] } split /,/, $exp_add_col;
        @anno_heads_surffix = map{ "gene_$_" } @anno_heads;
    }

    # 文件输出
    my $exp_matrix  = "$exp_cis_prefix.raw.txt";
    my $exp_pos     = "$exp_cis_prefix.pos.txt";
    my $exp_anno    = "$exp_cis_prefix.anno.txt";
    print "output $exp_matrix\n";
    print "output $exp_pos\n";
    print "output $exp_anno\n";

    open RAW, ">$exp_matrix";
    open POS, ">$exp_pos";
    open ANNO, ">$exp_anno";

    print RAW "EXP_ID\t" . (join "\t", @samples) . "\n";
    print POS "EXP_ID\tChr\tStart\tEnd\n";
    print ANNO "EXP_ID\t" . (join "\t", @anno_heads_surffix) . "\n";

    foreach my $exp_id(sort {$hashEXP{'DATA'}{$a}{'sort_value'} <=> $hashEXP{'DATA'}{$a}{'sort_value'}} keys %{$hashEXP{'DATA'}})
    {
        my $chr   = $hashEXP{'DATA'}{$exp_id}{$hashEXP{'HEAD'}->[1]};  # 取出第一列 chr
        my $start = $hashEXP{'DATA'}{$exp_id}{$hashEXP{'HEAD'}->[2]};  # 取出第二列 start
        my $end   = $hashEXP{'DATA'}{$exp_id}{$hashEXP{'HEAD'}->[3]};  # 取出第三列 end

        my @exps = map{ $hashEXP{'DATA'}{$exp_id}{$_} } @samples;
        my @annos = map{ $hashEXP{'DATA'}{$exp_id}{$_} } @anno_heads;

        print RAW "$exp_id\t" . (join "\t", @exps) . "\n";
        print POS "$exp_id\t$chr\t$start\t$end\n";
        print ANNO "$exp_id\t" . (join "\t", @annos) . "\n";
    }
    close RAW;
    close POS;
    close ANNO;
}
sub prepare_snp{
    my $snp_file       = shift @_;
    my $hashSample     = shift @_;
    my $snp_cis_prefix = shift @_;
    my @samples        = sort {$hashSample{'DATA'}{$a}{'sort_value'} <=> $hashSample{'DATA'}{$b}{'sort_value'}}keys %{$hashSample{'DATA'}};

    # 任务:转成0、1、2 构成的矩阵
    print "[process] Prepare snp data\n";
    my %hashSNP = read_matrix($snp_file, 0);

    # 确定要添加的snp注释
    my @anno_heads = ();
    my @anno_heads_surffix = ();  # 添加前缀，防止与exp冲突
    if($snp_add_col ne "")
    {
        @anno_heads = map{ $hashSNP{'HEAD'}->[$_] } split /,/, $snp_add_col;
        @anno_heads_surffix = map{ "snp_$_" } @anno_heads;
    }
    
 
    
    # 统计每个位点的等位基因
    print "collect allele \n";
    my %hashAllele;
    foreach my $snp_id(keys %{$hashSNP{'DATA'}})
    {
        foreach my $sample(@samples)
        {
            my $geno = $hashSNP{'DATA'}{$snp_id}{$sample};
            next if($geno !~ /\//);
            my ($allele1, $allele2) = split /\//, $geno;
            $hashAllele{$snp_id}{$allele1}++;
            $hashAllele{$snp_id}{$allele2}++;
        }
    }

    # 确定 ref / alt
    print "define ref/alt by maf \n";
    my %hashCode;
    my $drop = 0;
    foreach my $snp_id(keys %hashAllele)
    {
        my @alleles = sort {$hashAllele{$snp_id}{$b} <=> $hashAllele{$snp_id}{$a}} keys %{$hashAllele{$snp_id}};
        if(scalar(@alleles) > 2)
        {
            die "[Error] snp位点 $snp_id 存在3个或以上的等位基因数量，不能分析。 仅允许2个等位基因存在。\n";
        }
        if(scalar(@alleles) == 1)
        {
            print "[Warning] snp位点 $snp_id 仅有一个等位基因，不能做qtl分析，从结果里删除。\n";
            $drop++;
            next;
        }
        
        my $ref = $alleles[0];
        my $alt = $alleles[1];
        $hashCode{$snp_id}{"$ref/$ref"} = 0;
        $hashCode{$snp_id}{"$ref/$alt"} = 1;
        $hashCode{$snp_id}{"$alt/$ref"} = 1;
        $hashCode{$snp_id}{"$alt/$alt"} = 2;
        $hashCode{$snp_id}{""} = "";  # 空字符
        $hashCode{$snp_id}{" "} = "";  # 空字符
    }
    print "[Warning] 有 $drop 个位点因为只有一个等位基因，从文件中排除\n" if($drop > 0);

    # 文件输出
    my $snp_matrix_raw   = "$snp_cis_prefix.raw.txt";
    my $snp_matrix_code  = "$snp_cis_prefix.code.txt";
    my $snp_pos          = "$snp_cis_prefix.pos.txt";
    my $snp_anno         = "$snp_cis_prefix.anno.txt";
    print "output $snp_matrix_raw\n";
    print "output $snp_matrix_code\n";
    print "output $snp_pos\n";
    print "output $snp_anno\n";

    open RAW, ">$snp_matrix_raw";
    open CODE, ">$snp_matrix_code";
    open POS, ">$snp_pos";
    open ANNO, ">$snp_anno";

    print RAW "SNP_ID\t" . (join "\t", @samples) . "\n";
    print CODE "SNP_ID\t" . (join "\t", @samples) . "\n";
    print POS "SNP_ID\tChr\tPos\n"; 
    print ANNO "SNP_ID\t" . (join "\t", @anno_heads_surffix) . "\n";

    foreach my $snp_id(sort {$hashSNP{'DATA'}{$a}{'sort_value'} <=> $hashSNP{'DATA'}{$a}{'sort_value'}} keys %{$hashSNP{'DATA'}})
    {
        next if(not exists $hashCode{$snp_id});
        my $chr = $hashSNP{'DATA'}{$snp_id}{$hashSNP{'HEAD'}->[1]};  # 取出第一列 chr
        my $pos = $hashSNP{'DATA'}{$snp_id}{$hashSNP{'HEAD'}->[2]};  # 去除第二列 pos

        my @raws = map{ $hashSNP{'DATA'}{$snp_id}{$_} } @samples;
        my @codes = map{ $hashCode{$snp_id}{$_} } @raws;
        my @annos = map{ $hashSNP{'DATA'}{$snp_id}{$_} } @anno_heads;

        print RAW "$snp_id\t" . (join "\t", @raws) . "\n";
        print CODE "$snp_id\t" . (join "\t", @codes) . "\n";
        print POS "$snp_id\t$chr\t$pos\n";
        print ANNO "$snp_id\t" . (join "\t", @annos) . "\n";
    }
    close RAW;
    close CODE;
    close POS;
    close ANNO;
}


sub check_sample{
    my $hashSample = shift @_;
    my $file       = shift @_;
    my $title      = shift @_;

    print "[process] Check sample in $title  ";

    open FILE, $file;
    my $line1=<FILE>;
       $line1=~s/[\r\n]//g;
    my %hashFileSample = map{ ($_, 1) } split /\t/, $line1;
    close FILE;

    my @errors = grep{ not exists $hashFileSample{$_} } keys %{$hashSample->{'DATA'}};

    if(scalar(@errors) > 0)
    {   
        print "\n[Error] sample lost in $title: @errors\n";
    }else{
        print "[OK]\n";
    } 

}


# 读取带有表头的矩阵
sub read_matrix{
    my $file = shift @_;
    my $key_col = shift @_;

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
        my $key_value = $datas[$key_col];
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

# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}





sub format_run{
	my ($workbook)=@_;
    my %format=();
    $format{'title'} = $workbook->add_format();
    $format{'title'} ->set_align('center');
    $format{'title'} ->set_align('vcenter');
    $format{'title'} ->set_size(12);
    $format{'title'} ->set_font("Times New Roman");
    $format{'title'} ->set_border();
    $format{'title'} ->set_bg_color("yellow");
    $format{'title'} ->set_color("black");

    $format{'normal'} = $workbook->add_format();
    $format{'normal'} ->set_align('center');
    $format{'normal'} ->set_align('vcenter');
    $format{'normal'} ->set_size(12);
    $format{'normal'} ->set_font("Times New Roman");
    $format{'normal'} ->set_border();

    $format{'orange'} = $workbook->add_format();
    $format{'orange'} ->set_align('vcenter');
    $format{'orange'} ->set_size(12);
    $format{'orange'} ->set_font("Times New Roman");
    $format{'orange'} ->set_bg_color("#fac090");
    $format{'orange'} ->set_border();

    $format{'readme1'} = $workbook->add_format();
    $format{'readme1'}->set_align('center');
    $format{'readme1'}->set_align('vcenter');
    $format{'readme1'}->set_bold();
    $format{'readme1'}->set_size(14);
    $format{'readme1'}->set_font("Times New Roman");
    $format{'readme1'}->set_border();

    $format{'readme2'} = $workbook->add_format();
    $format{'readme2'}->set_align('vcenter');
    $format{'readme2'}->set_bold();
    $format{'readme2'}->set_size(14);
    $format{'readme2'}->set_font("Times New Roman");

    $format{'readme2tmp'} = $workbook->add_format();
    $format{'readme2tmp'}->set_right();

    $format{'readme3'} = $workbook->add_format();
    $format{'readme3'}->set_align('center');
    $format{'readme3'}->set_align('vcenter');
    $format{'readme3'}->set_bold();
    $format{'readme3'}->set_size(11);
    $format{'readme3'}->set_font("Times New Roman");
    $format{'readme3'}->set_border();

    $format{'readme4'} = $workbook->add_format();
    $format{'readme4'}->set_align('vcenter');
    $format{'readme4'}->set_bold();
    $format{'readme4'}->set_size(11);
    $format{'readme4'}->set_font("Times New Roman");
    $format{'readme4'}->set_border();

    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    return %format;
}
