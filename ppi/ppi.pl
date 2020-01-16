# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Cwd qw( abs_path );

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"}; 

# 定义 -> 核心变量
my $tools_dir                = SCRIPTDIR . "/tools";
my $stringdb_script          = "$tools_dir/stringdb.r";                  # ppi分析工具
my $clusterProfiler_script   = "$tools_dir/clusterProfiler.R";           # 富集分析工具
my $kegg_pathway_png_address = "$tools_dir/kegg_pathway_png_address.pl"; # 生成kegg图的地址
my $download_kegg_png        = "$tools_dir/download_kegg_png.pl";        # kegg下载
my $enrich_xls_to_heatmap    = "$tools_dir/enrich_xls_to_heatmap.pl";    # 通路热图文件准备
my $enrich_heatmap           = "$tools_dir/enrich_heatmap.R";            # 通路热图绘制
my $table2excel              = "$tools_dir/table2excel.pl";              # 表格转excel工具
my $readme_kegg              = "$tools_dir/readme_kegg.txt";
my $readme_go                = "$tools_dir/readme_go.txt";

# 软件、环境设置
my $Rscript         = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $Rlib            = "/home/genesky/software/r/3.5.1/lib64/R/library";
$ENV{"R_LIBS"} = $Rlib; # R包路径
# 本流程需要的R包汇总
# library(STRINGdb)
# library(igraph)
# library(clusterProfiler)
# library(topGO)
# library(rlist)
# library(docopt)
# library(pheatmap)
 
# 检测 -> 脚本输入
my ($gene_list, $species, $output_dir, $enrichment, $keep_tmp, $if_help);
GetOptions(
    "gene_list|g=s"   => \$gene_list,
    "species|s=s"     => \$species,
    "output_dir|o=s"  => \$output_dir,
    "enrichment|e!"   => \$enrichment,
    "keep_tmp|k!"     => \$keep_tmp,
    "help|h"          => \$if_help,
);
die help() if(defined $if_help or (not defined $gene_list or not defined $output_dir));
$species = "Homo_sapiens" if(not defined $species);
$output_dir = Cwd::abs_path($output_dir);
###################################################################### 主程序

# (1) 结果目录创建
mkdir $output_dir if(not -e $output_dir);

# （2） PPI分析
system("$Rscript $stringdb_script $gene_list $species $output_dir");

# （3）每个簇富集分析
exit if(not defined $enrichment);
print "##### start enrichment analysis #####\n";

# （3.1）富集分析物种ID检测
my %hashPahtwayDB = get_enrichment_db(); # 富集分析物种数据库
my $OrgDb_name    = $hashPahtwayDB{$species}{"OrgDb"};
if(not defined $OrgDb_name)
{
    print "[warnings] do not support [$species] to do enrichment analysis\n";
    exit;
}

# （3.2）读取每个簇的基因
my $cluster_file  = "$output_dir/cluster.xls";
my %hashCluster   = read_cluster($cluster_file);

# （3.3）富集分析汇总文件
my @go_files;
my @kegg_files;
my @sheet_names;

# 临时目录，分析完成后，可以之间删除
my $tmp_dir = "$output_dir/tmp_delete";
mkdir $tmp_dir if(not -e $tmp_dir);

foreach my $cluster_id(sort { $a<=>$b } keys %hashCluster)
{   
    # 少于5个基因，不做
    my @genes = sort keys %{$hashCluster{$cluster_id}};
    next if(scalar(@genes) < 5); 
    print "##### [cluster_$cluster_id] #####\n";

    # 簇富集分析目录
    my $cluster_dir     = "$output_dir/cluster_$cluster_id";
    my $cluster_tmp_dir = "$tmp_dir/cluster_$cluster_id";
    mkdir $cluster_dir     if(not -e $cluster_dir);
    mkdir $cluster_tmp_dir if(not -e $cluster_tmp_dir); # 临时
    mkdir "$cluster_dir/KEGG_Pathway_Illustrations" if(not -e "$cluster_dir/KEGG_Pathway_Illustrations"); # 通路图
    mkdir "$cluster_dir/Custom_Figures" if(not -e "$cluster_dir/Custom_Figures"); # 绘图

    # 簇基因列表
    my $gene_list = "$cluster_tmp_dir/gene_list.txt";
    open GENELIST, ">$gene_list";
    map{ print GENELIST "$_\n" } @genes;
    close GENELIST;

    # 富集分析
    system("$Rscript $clusterProfiler_script -s $OrgDb_name $gene_list $cluster_tmp_dir");
    push @go_files, "$cluster_tmp_dir/go_enrichment.xls";
    push @kegg_files, "$cluster_tmp_dir/kegg_enrichment.xls";
    push @sheet_names, "cluster_$cluster_id";

    # 获取通路图下载地址，并上色
    system("perl $kegg_pathway_png_address -k $cluster_tmp_dir/kegg_enrichment.xls -c $gene_list -o $cluster_tmp_dir/kegg.urls.txt");
    system("perl $download_kegg_png -i $cluster_tmp_dir/kegg.urls.txt -o $cluster_tmp_dir/");
    
    # 结果拷贝
    system("cp $cluster_tmp_dir/png/*png $cluster_dir/KEGG_Pathway_Illustrations"); # 通路图
    system("cp $cluster_tmp_dir/GO_barplot.pdf $cluster_tmp_dir/go.bp.pdf  $cluster_tmp_dir/go.cc.pdf $cluster_tmp_dir/go.mf.pdf $cluster_tmp_dir/kegg_dotplot.pdf $cluster_dir/Custom_Figures");
    
    # kegg  heatmap
    system(" perl $enrich_xls_to_heatmap $cluster_tmp_dir/kegg_enrichment.xls > $cluster_tmp_dir/data.txt ");
    system(" $Rscript $enrich_heatmap $cluster_tmp_dir/data.txt $cluster_dir/Custom_Figures/kegg_heatmap.pdf");     
}

# 富集分析汇总
my $excel_go      = "$output_dir/GO_Enrichment_Summary.xlsx";
my $excel_kegg    = "$output_dir/KEGG_Enrichment_Summary.xlsx";

my $go_file    = join ",", @go_files;
my $kegg_file  = join ",", @kegg_files;
my $sheet_name = join ",", @sheet_names;
system("perl $table2excel -i $go_file,$readme_go     -s $sheet_name,README -o $excel_go");
system("perl $table2excel -i $kegg_file,$readme_kegg -s $sheet_name,README -o $excel_kegg");

# 删除中间结果
system("rm -r $tmp_dir") if(not defined $keep_tmp);
###################################################################### 子函数


sub read_cluster{
    my $cluster_file = shift @_;
    my %hashCluster;

    open CLUSTER, $cluster_file;
    <CLUSTER>;
    while(<CLUSTER>)
    {
        next if($_!~/\w/);
        $_=~s/[\r\n]//g;
        my ($cluster_id, $gene) = split /\t/, $_;
        $hashCluster{$cluster_id}{$gene}++;
    }
    close CLUSTER;
    return %hashCluster;
}

sub get_enrichment_db{
    my %hashPahtwayDB;
    

    $hashPahtwayDB{"Homo_sapiens"}{'OrgDb'}                    = "hsa";
    
    $hashPahtwayDB{"Mus_musculus"}{'OrgDb'}                    = "mmu";
   
    $hashPahtwayDB{"Rattus_norvegicus"}{'OrgDb'}               = "rno";
 

    return %hashPahtwayDB;
}

sub help{
    my $info = "
Program: ppi， 蛋白互作网络分析
Version: 2019-05-10
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --gene_list/-g    [必填] gene_list文件，第一列数据是基因名，不能有表头。第二列允许字符up/down(不区分大小写)，用于在通路图中标记颜色默认为红色（up=red,down=blue）
         --species/-s      物种名称，默认为人。人=Homo_sapiens，小鼠=Mus_musculus，大鼠=Rattus_norvegicus，拟南芥=Arabidopsis_thaliana，杨树=Populus_trichocarpa，
         --output_dir/-o   [必填] 结果输出路径
         --enrichment/-e   是否对基因簇做富集分析，默认不做。支持物种人，小鼠，大鼠
         --keep_tmp/-k     是否保留临时文件目录，默认删除。临时目录：tmp_delete
         --help/-h         查看帮助文档
         注意：GO分析要慢一些
    \n";
    return $info;
}
