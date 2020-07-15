# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

my $kegg_api_url = "http://rest.kegg.jp";  # kegg api接口

# 检测 -> 脚本输入
my ($output_dir, $kegg_org, $keep_tmp, $if_help);
GetOptions(
    "output_dir|o=s" => \$output_dir,
    "kegg_org=s"     => \$kegg_org,
    "keep_tmp|k"   => \$keep_tmp,
    "help|h"         => \$if_help,
);

die help() if(defined $if_help or (not defined $output_dir or not defined $kegg_org ));

###################################################################### 主程序
make_dir($output_dir);



###################
# （1）通路列表下载
###################
print "(1) 下载通路列表\n";
my $pathway_desc = "$output_dir/$kegg_org.pathway.desc.txt";
do{
    system("wget $kegg_api_url/list/pathway/$kegg_org -O $pathway_desc -q");
} while(is_file_ok($pathway_desc) == 0);
parse_pathway_desc($pathway_desc);  # 去掉额外的标签

###################
# (2) 下载通路与基因id的关系
###################
print "(2) 下载通路与基因id的关系\n";
my $pathway_geneid_link = "$output_dir/$kegg_org.pathway.geneid.map.txt";
do{
    system("wget $kegg_api_url/link/$kegg_org/pathway -O $pathway_geneid_link -q");
} while(is_file_ok($pathway_geneid_link) == 0);
parse_pathway_geneid($pathway_geneid_link);  # 去掉额外的标签

###################
# (3) 下载物种基因名称 -> kegg基因id映射表
###################
print "(3) 下载物种基因名称 -》kegg基因id映射表\n";
my $geneid_symbol_map_list = "$output_dir/$kegg_org.geneid.symbol.map.txt";
do{
    system("wget $kegg_api_url/list/$kegg_org -O $geneid_symbol_map_list -q");
} while(is_file_ok($geneid_symbol_map_list) == 0);
parse_geneid_symbol($geneid_symbol_map_list);  # 去掉额外的标签


print "\n\n final result \n";
print "kegg 通路描述信息：                  $pathway_desc \n";
print "kegg 通路->基因映射关系：             $pathway_geneid_link \n";
print "kegg geneid -> gene symbol 对应关系：$geneid_symbol_map_list \n";

###################################################################### 子函数
sub parse_geneid_symbol{
    my $geneid_symbol_map_list = shift @_;
    my $tmp_file = "$geneid_symbol_map_list.tmp";
    system("mv $geneid_symbol_map_list $tmp_file");
    open INPUT, $tmp_file;
    open OUTPUT, ">$geneid_symbol_map_list";
    print OUTPUT "geneid\tgenesymbol\tcount\n";
    while(<INPUT>)
    {
        $_=~s/[\r\n]//g;
        my ($gene_id, $gene_symbol_info) = split /\t/, $_;
        $gene_id=~s/$kegg_org://;
        my ($gene_symbol_list, $gene_desc) = split /;/, $gene_symbol_info;
        my $count = 0;
        foreach my $gene_symbol(split /,/, $gene_symbol_list)
        {   
            $count++;
            # 去掉多余的字符，且转换为大写
            $gene_symbol =~s/^\s+//;
            $gene_symbol =~s/\s+$//;
            $gene_symbol = uc($gene_symbol);

            print OUTPUT "$gene_id\t$gene_symbol\t$count\n";
        }
    }
    close INPUT;
    close OUTPUT;
    system("rm $tmp_file") if(not defined $keep_tmp);   
}

sub parse_pathway_geneid{
    my $pathway_geneid_link = shift @_;
    my $tmp_file = "$pathway_geneid_link.tmp";
    system("mv $pathway_geneid_link $tmp_file");
    open INPUT, $tmp_file;
    open OUTPUT, ">$pathway_geneid_link";
    print OUTPUT "pathway\tgeneid\n";
    while(<INPUT>)
    {
        $_=~s/[\r\n]//g;
        my ($pathway_id, $gene_id) = split /\t/, $_;
        $pathway_id=~s/^path://;
        $gene_id=~s/$kegg_org://;

        print OUTPUT "$pathway_id\t$gene_id\n";
    }
    close INPUT;
    close OUTPUT;
    system("rm $tmp_file") if(not defined $keep_tmp);   
}

sub parse_pathway_desc{
    my $pathway_desc = shift @_;
    my $tmp_file = "$pathway_desc.tmp";
    system("mv $pathway_desc $tmp_file");
    open INPUT, $tmp_file;
    open OUTPUT, ">$pathway_desc";
    print OUTPUT "pathway\tpathway_desc\n";
    while(<INPUT>)
    {
        $_=~s/[\r\n]//g;
        my ($pathway_id, $desc) = split /\t/, $_;
        $pathway_id=~s/^path://;

        # 去掉通路描述部分末尾的物种信息
        my @tmps = split /-/, $desc;
        my $desc_simplify = join "", @tmps[0..($#tmps - 1)];
        $desc_simplify=~s/\s+$//;

        print OUTPUT "$pathway_id\t$desc_simplify\n";
    }
    close INPUT;
    close OUTPUT;
    system("rm $tmp_file") if(not defined $keep_tmp); 
}


# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}

# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}


sub help
{
    my $info = "
Program: 下载、整理kegg pathway数据库 
Version: 2020-4-8
Contact: 129、甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --output_dir/-o    输出目录 , 例如： ./kegg_hsa_20200408
         --kegg_org         kegg的物种缩写编号， https://www.genome.jp/kegg/catalog/org_list.html 在这里查看
                            人： hsa
                            小鼠：mmu
                            大鼠：rno
                            兔子：ocu
         --keep_tmp/-k      保留临时文件。默认删除
         --help/-h          查看帮助文档。
    \n";
    return $info;
} 
