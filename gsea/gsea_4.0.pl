# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $gsea                = "/home/genesky/software/gsea/4.0.2/gsea-cli.sh";
my %hashDB              = database_list(); # 数据库。GSEA官网只支持human。其他物种可以在http://ge-lab.org/gskb/下载

# 检测 -> 脚本输入
my ($input_matrix, $sample_group, $gene_col_name, $database, $database_list, $min_gene_count, $output_dir, $if_help);
GetOptions(
    "input_matrix|i=s"   => \$input_matrix,
    "sample_group|s=s"   => \$sample_group,
    "gene_col_name|g=s"  => \$gene_col_name,
    "database|d=s"       => \$database,
    "database_list|l!"   => \$database_list,
    "min_gene_count|m=i" => \$min_gene_count,
    "output_dir|o=s"     => \$output_dir,
    "help|h"             => \$if_help,
);
 
die show_database_list(\%hashDB) if(defined $database_list); # 查看数据库
die help() if(defined $if_help or (not defined $input_matrix or not defined $sample_group or not defined $database or not defined $output_dir));
die "not find database\n" if(not exists $hashDB{$database} or not -e $hashDB{$database}); # 数据库不存在

$gene_col_name  = "" if(not defined $gene_col_name);
$min_gene_count = 15 if(not defined $min_gene_count);
###################################################################### 主程序
 
# (1) 数据读入
my %hashExpression = read_expression_data($input_matrix, $gene_col_name);
my %hashGroup      = read_sample_group($sample_group);
 
# (2) 生成GSEA输入文件
mkdir $output_dir if(not -e $output_dir);
my $gct_file = "$output_dir/data.gct";
my $cls_file = "$output_dir/data.cls";

my @groups       = sort keys %hashGroup; # 分组
my $group_count  = scalar(@groups); # 分组数量
my @samples      = map{ my @samples_tmp = sort keys %{$hashGroup{$_}}; @samples_tmp; } @groups; # 样本
my @sample_groups = map{ my $group = $_; my @tmps = map{ $group } (sort keys %{$hashGroup{$group}}); @tmps } @groups; # 每个样本所属分组
my $sample_count = scalar(@samples); # 样本数量

my @genes        = sort keys %hashExpression; # 基因名称
my $gene_count   = scalar(@genes); # 基因数量
check_is_sample_lost(\%hashExpression, $genes[0], \@samples); # 检查样本是否缺失
my $metric_para  = choose_metric_para(\%hashGroup); # GSEA 默认参数，要求每一组至少3个样本。如果不满足，需要更换-metric参数 # 详细信息参考https://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Metrics_for_Ranking


# 样本分组信息文件
print "Generate GSEA input file: $gct_file,$cls_file\n";
open CLS, ">$cls_file";
print CLS "$sample_count $group_count 1\n";
print CLS "# ". (join " ", @groups) . "\n";
print CLS (join " ", @sample_groups) . "\n";
close CLS;

# 样本表达量文件
open GCT, ">$gct_file";
print GCT "#1.2\n";
print GCT "$gene_count\t$sample_count\n";
print GCT "NAME\tDESCRIPTION\t" . (join "\t", @samples) . "\n";
foreach my $gene(@genes)
{   
    my @datas = ($gene, 'na');
       @datas = (uc($gene), 'na') if($database !~ /^Human/); # 人以外的物种数据库来自于GSKB，该数据库中所有基因都是大写
    foreach my $sample(@samples)
    {
        push @datas, $hashExpression{$gene}{$sample};
    }
    print GCT (join "\t", @datas) . "\n";
}
close GCT;

# (3) GSEA 分析
system(" $gsea GSEA -gmx $hashDB{$database} -collapse false -res $gct_file  -cls $cls_file -set_min  $min_gene_count -metric $metric_para  -out $output_dir -rpt_label GSEA");




###################################################################### 子函数
sub database_list{
    my %hashDB;

    $hashDB{'Human_KEGG'} = "/home/genesky/database/gsea/human/c2.cp.kegg.v7.0.symbols.gmt";

    $hashDB{'Mouse_GO'}        = "/home/genesky/database/gsea/mouse/MousePath_GO_gmt.gmt";
    $hashDB{'Mouse_Metabolic'} = "/home/genesky/database/gsea/mouse/MousePath_Metabolic_gmt.gmt";
    $hashDB{'Mouse_Pathway'}   = "/home/genesky/database/gsea/mouse/MousePath_Pathway_gmt.gmt";

    return %hashDB;
}

sub show_database_list{
    my $hashDB = shift @_;
    print "Support DB: \n";
    foreach my $db(sort keys %$hashDB)
    {
        print "\t$db\t$hashDB->{$db}\n";
    }
    return "\n";
}
sub check_is_sample_lost{
    my $hashExpression = shift @_;
    my $gene           = shift @_;
    my $samples        = shift @_;
    foreach my $sample(@$samples)
    {
        next if(exists $hashExpression{$gene}{$sample});
        die "Lost sample in Expression file: $sample\n";
    }
}
sub choose_metric_para{
    my $hashGroup = shift @_;

    my $is_over_3 = 1; # 是否每一组样本数都大于3
    foreach my $group(keys %$hashGroup)
    {
        my @samples = keys %{$hashGroup->{$group}};
        next if(scalar(@samples) >= 3);
        $is_over_3 = 0;
    }

    my $metric_para = "Signal2Noise"; # 默认参数
       $metric_para = "log2_Ratio_of_Classes" if($is_over_3 == 0);  

    return $metric_para;
}
sub read_sample_group{
    my $sample_group = shift @_;
    my %hashGroup;
    print "Read $sample_group\n";

    open GROUP, $sample_group;
    while(<GROUP>)
    {
        next if($_!~/\w/);
        $_=~s/[\r\n]//g;
        my ($sample, $group) = split /\t/, $_;
        $hashGroup{$group}{$sample}++;
    }
    close GROUP;
    return %hashGroup;    
}


sub read_expression_data{
    my $input_matrix  = shift @_;
    my $gene_col_name = shift @_;


    my %hashExpression;
    print "Read $input_matrix\n";

    open MATRIX, $input_matrix;
    my $line1 = <MATRIX>;
       $line1  =~s/[\r\n]//g;
    my @heads = split /\t/, $line1;
    while(<MATRIX>)
    {
        next if($_!~/\w/);
        $_=~s/[\r\n]//g;
        my @datas   = split /\t/, $_;
        my %hashTmp = map{ ($heads[$_], $datas[$_]) } (0..$#heads);
        my $gene    = ($gene_col_name eq "") ? $datas[0] : $hashTmp{$gene_col_name}; # 基因名称
        die "[Error] 在行$.中，基因名缺失，请仔细检查输入文件 : $input_matrix\n" if(not defined $gene or $gene eq "");

        foreach my $col(0..$#heads)
        {
            $hashExpression{$gene}{$heads[$col]} = $datas[$col];
        }
    }
    close MATRIX;
    return %hashExpression;
}

sub help{
    my $info = "
Program: GSEA， GSEA分析, 基于GSEA 4.0.2版本
Version: 2019-10-21
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

    [必填]
         --input_matrix/-i    样本表达矩阵，每一行为基因，每一列为样本，第一行为样本名，第一列为基因名。推荐FPKM
         --sample_group/-s    样本分组文件，两列，第一列样本名，第二列分组名,不能有表头。注：分组数量允许>2。（流程会根据样本名，从input_matrix中提取需要的样本）
         --database/-d        数据库名称，例如：Human_KEGG
         --output_dir/-o      输出目录

    [选填]
         --gene_col_name/-g   指定基因列名。如果输入矩阵第一列不是基因，可以在这里指定基因列名，例如：-g gene_name
         --min_gene_count/-m  数据库中的一条通路里，最低基因数量，低于阈值，则不分析。默认：15
         --database_list/-l   列出支持的数据库
         --help/-h            查看帮助文档
    \n";
    return $info;
}

# GSEA 输入文件有2个
# gct与cls文件
# 
# gct格式示例
# 第一行固定，'#1.2'
# 第二行，基因、样本数目
# 第三行，NAME,DESCRIPTION,样本1，样本2.。。
# 第四行，每一个基因的表达量
# #1.2
# 19611   48
# NAME    DESCRIPTION     H1      H10     H11     H12     H13     H14     H15     H16     H17     H18     H19     H2      H3      H4      H5      H6      H7      H9      L1      L10     L11     L12     L13     L2      L3      L4      L5
# A1BG    na      0
# 
# cls格式示例
# 第一行：样本数、分类数、1
# 第二行：可视化分组名称。注：组数可以>2
# 第三行：每个样本的分组名称，可以与第二行不同，但是顺序要与第二行一致
# 48 2 1
# # HL N
# HL HL HL HL HL HL HL 