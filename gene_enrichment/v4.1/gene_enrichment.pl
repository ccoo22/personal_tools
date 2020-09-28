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
my $self_database_dir             = SCRIPTDIR . "/self_database";  # 自建数据库路径
my $stringdb_script               = "$tools_dir/stringdb.r";                  # ppi分析工具
my $clusterProfiler_script        = "$tools_dir/clusterProfiler.R";           # 富集分析工具
my $kegg_pathway_png_address      = "$tools_dir/kegg_pathway_png_address.pl"; # 生成kegg图的地址
my $download_kegg_png             = "$tools_dir/download_kegg_png.pl";        # kegg下载
my $kegg_pathway_download         = "$tools_dir/kegg_pathway_download.pl";    # kegg下载，升级版
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
my ($gene_list, $species, $output_dir, $enrichment, $to_excel, $no_download_png, $list, $if_help);
GetOptions(
    "gene_list|g=s"    => \$gene_list,
    "species|s=s"      => \$species,
    "output_dir|o=s"   => \$output_dir,
    "enrichment|e!"    => \$enrichment,
    "to_excel|e!"      => \$to_excel,
    "no_download_png!" => \$no_download_png,
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

my $no_download_para = (defined $no_download_png) ? "--no_download_png" : "";
###################################################################### 主程序

my $custom_figure_dir = "$output_dir/Custom_Figures";
mkdir $output_dir        if(not -e $output_dir);
mkdir $custom_figure_dir if(not -e $custom_figure_dir);# 个性化绘图结果



# 富集分析
my $go_raw   = "$output_dir/go_enrichment.xls";
my $kegg_raw = "$output_dir/kegg_enrichment.xls";
my $do_raw = "$output_dir/do_enrichment.xls";  # 人类疾病富集分析

system("$Rscript $clusterProfiler_script -i $gene_list -o $output_dir --org $hashPahtwayDB{$species}{'org'} --orgdb $hashPahtwayDB{$species}{'orgdb'} --rlib $Rlib ");
 
 
# 获取P<0.05通路图下载地址，并上色
# system("perl $kegg_pathway_png_address -k $kegg_raw -c $gene_list --add_address_to_excel -o $output_dir/kegg.urls.txt");  #  --add_address_to_excel 把通路地址加入kegg_raw的最后一列，列名： kegg_link
# system("perl $download_kegg_png -i $output_dir/kegg.urls.txt -o $output_dir/") if(not defined $no_download_png);  # 图片下载：

# 获取P<0.05通路图下载地址，并上色, 升级版，应对IP被禁
system("perl $kegg_pathway_download -k $kegg_raw -c $gene_list --add_address_to_excel $no_download_para -o $output_dir/png");

# 结果整理
system("mv $output_dir/png $output_dir/KEGG_Pathway_Illustrations") if(not defined $no_download_png); # 通路图
system("mv $output_dir/GO_barplot.pdf $output_dir/go.bp.pdf  $output_dir/go.cc.pdf $output_dir/go.mf.pdf $custom_figure_dir");
system("mv $output_dir/kegg_dotplot.pdf $custom_figure_dir");

# kegg热图
system(" perl $enrich_xls_to_heatmap $kegg_raw > $output_dir/data.kegg.heatmap.txt ");
system(" $Rscript $enrich_heatmap $output_dir/data.kegg.heatmap.txt $custom_figure_dir/kegg_heatmap.pdf");  


# 删除中间文件
# system("rm -r $output_dir/html");
# system("rm    $output_dir/kegg.urls.txt");
# system("rm    $output_dir/png.urls.txt") if(not defined $no_download_png);
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
    
    $hashPahtwayDB{"Mus_musculus"}{'org'}               = "mmu";
    $hashPahtwayDB{"Mus_musculus"}{'orgdb'}             = "/home/genesky/database_new/self_build_database/clusterprofiler/mouse.mmu.orgdb";

    $hashPahtwayDB{"Rattus_norvegicus"}{'org'}          = "rno";
    $hashPahtwayDB{"Rattus_norvegicus"}{'orgdb'}        = "/home/genesky/database_new/self_build_database/clusterprofiler/rat.rno.orgdb";

    $hashPahtwayDB{"rabbit"}{'org'}                     = "ocu";
    $hashPahtwayDB{"rabbit"}{'orgdb'}                   = "/home/genesky/database_new/self_build_database/clusterprofiler/rabbit_oryctolagus_cuniculus.ocu.orgdb";
 
    $hashPahtwayDB{"zea_mays"}{'org'}                   = "zma";
    $hashPahtwayDB{"zea_mays"}{'orgdb'}                 = "/home/genesky/database_new/self_build_database/clusterprofiler/zea_mays.zma.orgdb";

    $hashPahtwayDB{"escherichia_coli"}{'org'}           = "eco";
    $hashPahtwayDB{"escherichia_coli"}{'orgdb'}         = "/home/genesky/database_new/self_build_database/clusterprofiler/escherichia_coli.eco.orgdb";

    $hashPahtwayDB{"honey_bee_apis_mellifera"}{'org'}   = "ame";
    $hashPahtwayDB{"honey_bee_apis_mellifera"}{'orgdb'} = "/home/genesky/database_new/self_build_database/clusterprofiler/honey_bee_apis_mellifera.ame.orgdb";

    $hashPahtwayDB{"zebrafish_danio_rerio"}{'org'}      = "dre";
    $hashPahtwayDB{"zebrafish_danio_rerio"}{'orgdb'}    = "/home/genesky/database_new/self_build_database/clusterprofiler/zebrafish_danio_rerio.dre.orgdb";
    return %hashPahtwayDB;
}

sub help{
    my $info = "
Program: gene_enrichment, 基因富集分析GO/KEGG
Version: 2019-05-15
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [必填]
         --gene_list/-g     gene_list文件，第一列数据是基因名，不能有表头。第二列允许字符up/down(不区分大小写)，用于在通路图中标记颜色默认为红色（up=red,down=blue）
                            注意：基因严格区分大小写，务必与NCBI对应物种的基因名保持一致，否则可能匹配失败。
         --species/-s       物种名称，默认为人。人=Homo_sapiens，小鼠=Mus_musculus，大鼠=Rattus_norvegicus 
         --output_dir/-o    结果输出路径。流程自动创建

        [选填]
         --to_excel/-e     是否把富集分析结果转换成excel。默认不转换。
         --no_download_png 是否下载kegg官方图，默认：下载 。注：如果服务器IP被禁止，可能会导致下载死循环。   
         --list/-l         列出支持的物种
         --help/-h         查看帮助文档
         注意：
         （1）GO分析速度相对要慢一些
         （2）下载通路图标记颜色时，基因不能过多，当字符长度超过3045个时，kegg会拒绝访问，故，在下载的通路图中，我们的流程会自动放弃标注后面的部分基因。

    \n";
    return $info;
}
 
