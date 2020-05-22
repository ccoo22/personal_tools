# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

my $kegg_api_url = "http://rest.kegg.jp";  # kegg api接口
my $rdb_build    = SCRIPTDIR . "/update_kegg_metabolite_db.r";

# 检测 -> 脚本输入
my ($output_dir, $kegg_org, $if_help);
GetOptions(
    "output_dir|o=s" => \$output_dir,
    "kegg_org=s"     => \$kegg_org,
    "help|h"         => \$if_help,
);

die help() if(defined $if_help or (not defined $output_dir or not defined $kegg_org ));

###################################################################### 主程序
my $tmp1_dir = "$output_dir/tmp";   
my $tmp_dir  = "$output_dir/tmp/$kegg_org";  # 临时目录，按物种存放
make_dir($output_dir, $tmp1_dir, $tmp_dir);
 
# 最终需要的三个文件
my $compound_desc        = "$output_dir/kegg.$kegg_org.compound.desc.txt";  # 物种包含的代谢物
my $pathway_desc          = "$output_dir/kegg.$kegg_org.pathway.desc.txt";  # 物种的通路
my $pathway_compound_link = "$output_dir/kegg.$kegg_org.pathway.compound.map.txt";  # 物种的通路与代谢物映射关系
my $metabolite_db         = "$output_dir/kegg.$kegg_org.rda";  # 最终R格式数据库

###################
# （1）代谢物列表下载
###################
print "(1) 下载代谢物列表 [$kegg_org]\n";
my $compound_desc_tmp = "$tmp_dir/kegg.$kegg_org.compound.desc.txt";
do{
    system("wget $kegg_api_url/list/compound -O $compound_desc_tmp -q");  # compound 只能下载所有的，不能只下载对应物种的
} while(is_file_ok($compound_desc_tmp) == 0);


###################
# （2）通路列表下载
###################
print "(2) 下载通路列表 [$kegg_org]\n";
my $pathway_desc_tmp = "$tmp_dir/kegg.$kegg_org.pathway.desc.txt";
do{
    system("wget $kegg_api_url/list/pathway/$kegg_org -O $pathway_desc_tmp -q");  # 可以只下载对应物种的
} while(is_file_ok($pathway_desc_tmp) == 0);
 
 
###################
# （3）下载通路与代谢物的关系
###################
print "(3) 下载通路与代谢物的关系 [$kegg_org]\n";
my %hashPathway = read_pathway($pathway_desc_tmp);
my @pathways = sort keys %hashPathway;
foreach my $pathway(@pathways)
{
    process_bar_array($pathway, \@pathways);
    my $pathway_conf = "$tmp_dir/$pathway.conf";  # 记录通路元素构成
    my $pathway_kgml = "$tmp_dir/$pathway.kgml";  # 通路原始数据，记录方向关系。该文件用于分析拓扑结构
    do{
        system("wget $kegg_api_url/get/$pathway/conf -O $pathway_conf -q");  # 只能通过conf文件下载物种下的代谢物信息。link/cpd/pathway 的方法只能下载到所有map的，而不是对应物种的
    } while(is_file_ok($pathway_conf) == 0);

    do{
        system("wget $kegg_api_url/get/$pathway/kgml -O $pathway_kgml -q");  # 只能通过conf文件下载物种下的代谢物信息。link/cpd/pathway 的方法只能下载到所有map的，而不是对应物种的
    } while(is_file_ok($pathway_kgml) == 0);
}

###################
# （4）通路数据整理
###################
print "(4) 通路数据整理 [$kegg_org]\n";
my %hashCompound = read_compound($compound_desc_tmp);
my %hashLink = read_link(\@pathways, $tmp_dir);  # 从conf文件里，读入每个通路包含的代谢物

open PATHWAY,  ">$pathway_desc";
print PATHWAY "pathway\tpathway_desc\n";

open LINK,     ">$pathway_compound_link";
print LINK "pathway\tcompound\n";

open COMPOUND, ">$compound_desc";
print COMPOUND "compound\tcompound_desc\n";

my %hashCUnique;  # 需要输出的代谢物
foreach my $pathway(sort keys %hashLink)
{
    print PATHWAY "$pathway\t$hashPathway{$pathway}\n";
    foreach my $compound(sort keys %{$hashLink{$pathway}})
    {
        print LINK "$pathway\t$compound\n";
        $hashCUnique{$compound}++;
    }
}
foreach my $compound(sort keys %hashCUnique)
{   
    print COMPOUND "$compound\t$hashCompound{$compound}\n";
}

close PATHWAY;
close LINK;
close COMPOUND;


###################
# （5）R语言数据库整理
###################
print "(5) R 通路数据整理 [$kegg_org]\n";
system("$rdb_build --compound_desc $compound_desc --pathway_desc $pathway_desc --pathway_compound_map $pathway_compound_link --kgml_dir $tmp_dir --output_file $metabolite_db ");

######################
print "\n\nfinal result \n";
print "kegg 代谢物描述：          $compound_desc \n";
print "kegg 通路描述：            $pathway_desc \n";
print "kegg 通路->代谢物 映射关系：$pathway_compound_link \n";
print "kegg 通路->代谢物 R数据库： $metabolite_db \n";


###################################################################### 子函数

# 读入通路包含的代谢物
sub read_link{
    my $pathways = shift @_;
    my $tmp_dir  = shift @_;

    my %hashLink;
    foreach my $pathway(@$pathways)
    {
        my $pathway_conf = "$tmp_dir/$pathway.conf";
        open CONF, $pathway_conf;
        while(<CONF>)
        {
            $_=~s/[\r\n]//g;
            my ($shape, $add, $value)= split /\t/, $_;
            next if($shape !~ /^circ/);  # circ 表示的是代谢物C/多聚糖G
            next if($value !~/^C/); # 仅保留C
            my ($compound, $desc) = split /\s+/, $value;
            $hashLink{$pathway}{$compound}++;
        }
        close CONF;
    }
    return %hashLink;
}

# 读入代谢物描述
sub read_compound{
    my $compound_list_tmp = shift @_;

    my %hashCompound;
    open COMPOUND, $compound_list_tmp;
    while(<COMPOUND>)
    {
        $_=~s/[\r\n]//g;
        my ($compound, $desc) = split /\t/, $_;
        $compound =~s/cpd://;
        $hashCompound{$compound} = $desc;
    }
    close COMPOUND;
    return %hashCompound;
}

# 读入通路描述
sub read_pathway{
    my $pathway_desc_tmp = shift @_;

    my %hashPathway;
    open PATHWAY, $pathway_desc_tmp;
    while(<PATHWAY>)
    {
        $_=~s/[\r\n]//g;
        my ($pathway, $desc) = split /\t/, $_;
        $pathway =~s/path://;

        # 去掉通路描述信息尾部的物种信息
        my @tmps = split /-/, $desc;
        my $desc_clean = join "", @tmps[0..($#tmps-1)];

        $hashPathway{$pathway} = $desc_clean;
    }
    close PATHWAY;
    return %hashPathway;
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
        mkdir  $dir if(not -e $dir);
    }
}


sub help
{
    my $info = "
Program: 更新kegg compound、pathway数据库 
Version: 2020-4-8
Contact: 129、甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
         --kegg_org         kegg的物种缩写编号， https://www.genome.jp/kegg/catalog/org_list.html 在这里查看
                            人： hsa
                            小鼠：mmu
                            大鼠：rno
                            兔子：ocu

         --output_dir/-o    输出目录 , 例如： ./
                            会生成四个文件 + 1个中间结果目录
                            kegg.kegg_org.compound.desc.txt  # 物种相关的所有代谢物
                            kegg.kegg_org.pathway.desc.txt   # 物种相关的所有通路
                            kegg.kegg_org.pathway.compound.map.txt  # 物种所有通路-> 代谢物映射关系
                            kegg.kegg_org.rda  # 基于上面的结果 + 通路kgml文件整理的R数据库。便于后续富集分析统一调用。该数据结构参考了MetaboAnalyst软件
                                               # 需要注意的点：（1）参考官方网站的结果，他们把通路中没有连通的代谢物去掉了（网络关系图可使用KEGGgraph提取，参考我的R脚本），我们没有去除。所以富集分析时，背景有所不同，分析的p值也会不同（例如hsa00010中，实际含有代谢物31个，而只有26个代谢物连通。故官方的分析背景只有26个，我们的就是31个）
                            ./tmp/kegg_org  # 保存了下载的通路conf与kgml文件。conf用于提取同于与代谢物映射关系，kgml用于分析网络拓扑，计算impact。
         --help/-h          查看帮助文档。
    \n";
    return $info;
} 