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
my $tools_dir                     = SCRIPTDIR . "/tools";
my $clusterProfiler_script        = "$tools_dir/clusterProfiler.R";           # 富集分析工具
my $enrich_xls_to_heatmap         = "$tools_dir/enrich_xls_to_heatmap.pl";    # 通路热图文件准备
my $enrich_heatmap                = "$tools_dir/enrich_heatmap.R";            # 通路热图绘制
my $table2excel                   = "$tools_dir/table2excel.pl";              # 表格转excel工具
my $readme_kegg                   = "$tools_dir/readme_kegg.txt";
my $readme_go                     = "$tools_dir/readme_go.txt";
my $readme_do                     = "$tools_dir/readme_do.txt";

# 软件、环境设置
my $Rscript         = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $Rlib            = "/home/genesky/software/r/3.5.1/lib64/R/library";
$ENV{"R_LIBS"} = $Rlib; # R包路径
# 本流程需要的R包汇总
# library(STRINGdb)
# library(igraph)
# library(Rgraphviz)
# library(clusterProfiler)
# library(topGO)
# library(rlist)
# library(docopt)
# library(pheatmap)
 
# 检测 -> 脚本输入
my ($gene_list, $species, $output_dir, $enrichment, $to_excel, $list, $if_help);
GetOptions(
    "gene_list|g=s"    => \$gene_list,
    "species|s=s"      => \$species,
    "output_dir|o=s"   => \$output_dir,
    "enrichment|e!"    => \$enrichment,
    "to_excel|e!"      => \$to_excel,
    "list|l!"          => \$list,
    "help|h"           => \$if_help,
);
show_list() if(defined $list);
die help() if(defined $if_help or (not defined $gene_list or not defined $output_dir));
$species = "Homo_sapiens" if(not defined $species);
$output_dir = Cwd::abs_path($output_dir);
 
# 富集分析物种ID检测
my %hashPahtwayDB = get_enrichment_db(); # 富集分析物种数据库
if(not exists $hashPahtwayDB{$species})
{
    print "[warnings] do not support [$species] to do enrichment analysis\n";
    exit;
}

###################################################################### 主程序

my $custom_figure_dir = "$output_dir/Custom_Figures";
mkdir $output_dir        if(not -e $output_dir);
mkdir $custom_figure_dir if(not -e $custom_figure_dir);# 个性化绘图结果



# 富集分析
my $go_raw   = "$output_dir/go_enrichment.xls";
my $kegg_raw = "$output_dir/kegg_enrichment.xls";
my $do_raw = "$output_dir/do_enrichment.xls";  # 人类疾病富集分析

system("$Rscript $clusterProfiler_script -i $gene_list -o $output_dir --org $hashPahtwayDB{$species}{'org'} --orgdb $hashPahtwayDB{$species}{'orgdb'} --orgdb_keggid $hashPahtwayDB{$species}{'orgdb_keggid'} --rlib $Rlib ");
 
# 结果整理
system("mv $output_dir/GO_barplot.pdf $output_dir/go.bp.pdf  $output_dir/go.cc.pdf $output_dir/go.mf.pdf $custom_figure_dir/");
system("mv $output_dir/kegg_dotplot.pdf $custom_figure_dir/");
system("mv $output_dir/do_dotplot.pdf $custom_figure_dir/") if($species eq "Homo_sapiens");

# kegg热图
system(" perl $enrich_xls_to_heatmap $kegg_raw > $output_dir/data.kegg.heatmap.txt ");
system(" $Rscript $enrich_heatmap $output_dir/data.kegg.heatmap.txt $custom_figure_dir/kegg_heatmap.pdf");  


# 删除中间文件
system("rm    $output_dir/data.kegg.heatmap.txt");

exit if(not defined $to_excel);

# 富集分析结果转换成excel
my $excel_go      = "$output_dir/go_enrichment.xlsx";
my $excel_kegg    = "$output_dir/kegg_enrichment.xlsx";
my $excel_do    = "$output_dir/do_enrichment.xlsx";
system("perl $table2excel -i $go_raw,$readme_go     -s GO,README -o $excel_go");
system("perl $table2excel -i $kegg_raw,$readme_kegg -s KEGG,README -o $excel_kegg");
system("perl $table2excel -i $do_raw,$readme_do -s DO,README -o $excel_do") if($species eq 'Homo_sapiens');
system("rm    $go_raw");
system("rm    $kegg_raw");
system("rm    $do_raw") if($species eq 'Homo_sapiens');

###################################################################### 子函数
sub show_list{
    my %hashPahtwayDB = get_enrichment_db(); # 富集分析物种数据库
    print "当前流程支持的物种：\n";
    foreach my $db(keys %hashPahtwayDB)
    {
        print "    $db\n";
    }
    die "\n";
}
 

sub get_enrichment_db{
    my %hashPahtwayDB;
    

    $hashPahtwayDB{"Homo_sapiens"}{'org'}               = "hsa";
    $hashPahtwayDB{"Homo_sapiens"}{'orgdb'}             = "/home/genesky/database_new/self_build_database/clusterprofiler/human.hsa.orgdb";
    $hashPahtwayDB{"Homo_sapiens"}{'orgdb_keggid'}      = "ENTREZID";
    
    $hashPahtwayDB{"Mus_musculus"}{'org'}               = "mmu";
    $hashPahtwayDB{"Mus_musculus"}{'orgdb'}             = "/home/genesky/database_new/self_build_database/clusterprofiler/mouse.mmu.orgdb";
    $hashPahtwayDB{"Mus_musculus"}{'orgdb_keggid'}      = "ENTREZID";

    $hashPahtwayDB{"Rattus_norvegicus"}{'org'}          = "rno";
    $hashPahtwayDB{"Rattus_norvegicus"}{'orgdb'}        = "/home/genesky/database_new/self_build_database/clusterprofiler/rat.rno.orgdb";
    $hashPahtwayDB{"Rattus_norvegicus"}{'orgdb_keggid'} = "ENTREZID";

    $hashPahtwayDB{"rabbit"}{'org'}                     = "ocu";
    $hashPahtwayDB{"rabbit"}{'orgdb'}                   = "/home/genesky/database_new/self_build_database/clusterprofiler/rabbit_oryctolagus_cuniculus.ocu.orgdb";
    $hashPahtwayDB{"rabbit"}{'orgdb_keggid'}            = "ENTREZID";

    $hashPahtwayDB{"zea_mays"}{'org'}                   = "zma";
    $hashPahtwayDB{"zea_mays"}{'orgdb'}                 = "/home/genesky/database_new/self_build_database/clusterprofiler/zea_mays.zma.orgdb";
    $hashPahtwayDB{"zea_mays"}{'orgdb_keggid'}          = "ENTREZID";

    $hashPahtwayDB{"escherichia_coli"}{'org'}           = "eco";
    $hashPahtwayDB{"escherichia_coli"}{'orgdb'}         = "/home/genesky/database_new/self_build_database/clusterprofiler/escherichia_coli.eco.orgdb";
    $hashPahtwayDB{"escherichia_coli"}{'orgdb_keggid'}  = "ENTREZID";

    $hashPahtwayDB{"honey_bee_apis_mellifera"}{'org'}               = "ame";
    $hashPahtwayDB{"honey_bee_apis_mellifera"}{'orgdb'}             = "/home/genesky/database_new/self_build_database/clusterprofiler/honey_bee_apis_mellifera.ame.orgdb";
    $hashPahtwayDB{"honey_bee_apis_mellifera"}{'orgdb_keggid'}      = "ENTREZID";

    $hashPahtwayDB{"zebrafish_danio_rerio"}{'org'}               = "dre";
    $hashPahtwayDB{"zebrafish_danio_rerio"}{'orgdb'}             = "/home/genesky/database_new/self_build_database/clusterprofiler/zebrafish_danio_rerio.dre.orgdb";
    $hashPahtwayDB{"zebrafish_danio_rerio"}{'orgdb_keggid'}      = "ENTREZID";

    return %hashPahtwayDB;
}

sub help{
    my $info = "
Program: gene_enrichment, 基因富集分析GO/KEGG
Version: 2020-12-04
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [必填]
         --gene_list/-g     gene_list文件，第一列数据是基因名，不能有表头。第二列允许字符up/down(区分大小写)，用于在通路图中标记颜色默认为红色（up=red,down=blue）
                            注意：基因严格区分大小写，务必与NCBI对应物种的基因名保持一致，否则可能匹配失败。
         --species/-s       物种名称，默认为人。人=Homo_sapiens，小鼠=Mus_musculus，大鼠=Rattus_norvegicus 
         --output_dir/-o    结果输出路径。流程自动创建

        [选填]
         --to_excel/-e     是否把富集分析结果转换成excel。默认不转换。
         --list/-l         列出支持的物种
         --help/-h         查看帮助文档
         
         注意：
             （1）服务器需要联网，kegg联网获取通路
             （2）GO分析速度相对要慢一些

    \n";
    return $info;
}
 
