$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
基于GB文件做进化树分析（基于基因聚类）
Version: v1.0 2020-11-24
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use Cwd 'abs_path';

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];

# 变量
# my $SOFT_ORTHOFINDER = "source /home/genesky/software/conda/4.8.3/bin/activate orthofinder && orthofinder";
my $SOFT_ORTHOFINDER = "source /home/genesky/software/conda/4.5.12/bin/activate /home/genesky/software/orthofinder/2.2.7 && orthofinder";
my $SOFT_R           = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $SOFT_R_LIB       = "/home/genesky/database/r/3.5.1";
my $DEFAULT_THREAD   = 10;
my $DEFAULT_METHOD   = 'msa';
my $SOFT_IQTREE      = "/home/genesky/software/iqtree/1.6.12/bin/iqtree";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($genbank_dir, $output_dir, $thread, $method, $if_help);
GetOptions(
    "genbank-dir|g=s" => \$genbank_dir,
    "output-dir|o=s" => \$output_dir,
    "method|m=s" => \$method,
    "thread|t=i" => \$thread,
    "help|h" => \$if_help,
);
die "
Options: 必填

        --genbank-dir/-g           GenBank 目录, 其中至少包含两个 *.gb 文件
                                   进化树样本命名，采用gb文件内部的 ORGANISM 所在行后面的信息，特殊字符会用_替换
        --output-dir/-o            结果输出路径，不用提前创建，流程自己生成

Options: 可选
        --method/-m                基因进化树创建方式， dendroblast or msa. 默认： $DEFAULT_METHOD
                                   注意： dendroblast 方法只生成进化树 newick文件
                                          msa         方法除了newick，还提供多重序列比对文件
        --thread/-t                线程数, 默认： $DEFAULT_THREAD
        --help/-h                  查看帮助文档
                                   基于软件： https://github.com/davidemms/OrthoFinder
\n" if (defined $if_help or not defined $genbank_dir or not defined $output_dir);
$thread = $DEFAULT_THREAD if(not defined $thread);
$method = $DEFAULT_METHOD if(not defined $method);

$output_dir = abs_path($output_dir);
mkdir $output_dir if(not -e $output_dir);

###################################################################### 初始化
my @genebanks = get_genbank($genbank_dir);
 
# 检查文件
foreach my $genebank(@genebanks){
    die "[ERROR] GeneBank 文件错误, $genebank\n" if(not is_genbank_ok($genebank));
}


###################################################################### 主程序
my $gb_amino_acid_dir = "$output_dir/gb_amino_acid";
if(-e $gb_amino_acid_dir){
    print "[warning] 已存在临时目录， 请务必保证这个目录是空的，否则可能出现异常： $gb_amino_acid_dir\n";
    print "是否删除该目录？[y/n] ";
    my $input = <STDIN>;
    $input=~s/\s//g;
    $input=lc($input);
    if($input eq 'y')
    {
        system("rm -rf $gb_amino_acid_dir");
        mkdir $gb_amino_acid_dir;
    }
}else{
    mkdir $gb_amino_acid_dir;
}

##########
# (1) genebank 文件提取基因的氨基酸序列，生成fasta文件
##########
print "[process] 提取gb文件中基因的氨基酸序列，制作fasta\n";
my $gb_count = @genebanks;
my $count = 0;
my %hashName;  # 检查ORGANISM是否有重名
foreach my $genebank(@genebanks){
    $count++;
    print "    [gb] $count/$gb_count $genebank ... ";
    genebank_2_faa($genebank, $gb_amino_acid_dir, \%hashName);
    print " [OK]\n";
}
 
##########
# (2) ORTHOFINDER 基因同源序列分析
##########
print "[process] ORTHOFINDER 基因同源序列分析, 进化树构建模型： $method\n";
my $log = "$gb_amino_acid_dir/log.txt";
system "$SOFT_ORTHOFINDER -f $gb_amino_acid_dir -t $thread -a $thread -M $method > $log";
my $tree = get_ortho_tree($log);  # 从log里确定进化树文件路径
die "[Error] 运行失败，请检查log日志： $log\n" if($tree eq '0');

my $tree_newick = "$output_dir/ortho_tree.newick";
my $tree_msa    = "$output_dir/ortho_tree.msa.fa";
system "cp $tree $tree_newick";
# 多重序列比对文件获取
if($method eq 'msa'){
    my $msa = $tree;
    $msa=~s/SpeciesTree_rooted.txt/Alignments\/SpeciesTreeAlignment.fa/;  
    system "cp $msa $tree_msa";

    # 重做进化树，补充bootstrap
    system "$SOFT_IQTREE -s $msa -m MFP -bb 1000  -bnni  -redo -nt 8"; 
    system "cp $msa.treefile $tree_newick";
}


##########
# (3) 进化树绘图
##########
my $tree_pdf = "$output_dir/ortho_tree.pdf";
plot_tree("$tree_newick", "$gb_amino_acid_dir/tree.r", "$tree_pdf");

print "[output] tree file:  $tree_newick\n";
print "[output] tree plot:  $tree_pdf\n";
print "[output] tree MSA:  $tree_msa\n" if($method eq 'msa');
print "\n[OK] 运行完成\n";
###################################################################### 子程序

sub get_genbank{
    my $dir = shift;
    my @array;
    opendir DIR, $dir;
    while(my $file = readdir DIR){
        next if($file !~ /\.gb$/ and $file !~ /\.gbf$/);
        push @array, abs_path("$dir/$file");
    }
    closedir DIR;
    return @array;
}

sub is_genbank_ok{
    my $file = shift;
    my $isok = 0;
    if(-e $file){
        open FILE, $file;
        while(my $line = <FILE>){
            $isok = 1 if($line =~ /^\/\//);
        }
        close FILE;
    }
    return $isok;
}

sub genebank_2_faa{
    my $file    = shift;
    my $outdir  = shift;
    my $hashName = shift;

    # 基于gb文件的ORGANISM内容设定样本名称
    my $species = "Unknow";  # 默认
    open FILE, $file;
    while(my $line = <FILE>){
        if($line =~ /^\s*ORGANISM\s+(.+)/){
            $species = $1;
            $species =~ s/^\s+//;
            $species =~ s/\s+$//;
            $species =~ s/\W/_/g;
        }
    }
    close FILE;
    print "$species ";

    if(exists $hashName->{$species})
    {
        die "[Error] GB 文件中 ORGANISM 名称与 $hashName->{$species} 重名，请修改GB文件中 ORGANISM信息 \n ";
    }else{  
        $hashName->{$species} = $file;
    }
    # 读入GB, 提取每个基因的氨基酸序列
    my $io_gb = Bio::SeqIO->new(-file => $file, -format => "genbank");
    my $io_seq = $io_gb->next_seq;
    my @feats = $io_seq->get_SeqFeatures;
    my $id = 1;
    open SAVE, ">$outdir/$species.faa";
    foreach my $feat(@feats){
        my $primary_tag = $feat->primary_tag;
        next if($primary_tag ne "CDS");
        next if($feat->has_tag("pseudo") or $feat->has_tag("pseudogene"));
        my $gene = ($feat->has_tag("gene")) ? ($feat->get_tag_values("gene"))[0] : 'NA';
        die "[ERROR] GenBank 缺失 translation 标签 $gene, $file\n" if(not $feat->has_tag("translation"));
        my $translation = ($feat->get_tag_values("translation"))[0];
        print SAVE ">$species|Gene_$id|$gene\n$translation\n"; $id++;
    }
    close SAVE;
}

sub get_ortho_tree{
    my $file = shift;
    open FILE, $file;
    while(my $line = <FILE>){
        if($line =~ /Rooted species tree:/){
            my $tree = <FILE>;
            $tree =~ s/^\s+//;
            $tree =~ s/\s+$//;
            return $tree;
        }
    }
    close FILE;
    return 0;
}

sub plot_tree{
    my $nwk = shift;
    my $r_script = shift;
    my $pdf = shift;
    my $script = "
.libPaths('$SOFT_R_LIB')
library(ggtree)

tree <- read.tree('$nwk')

max_width = max(ape::node.depth.edgelength(tree))  # 获取最深的边
max_height = max(ape::node.height(tree))

pdf_height = round(max_height*0.25, 0)
if(pdf_height < 6 ) pdf_height = 6

pdf('$pdf', height = pdf_height)
ggtree(tree) +
    geom_tippoint(aes(color=label), size=4, alpha=.4) +
    geom_tiplab() +
    xlim(NA, max_width * 1.4)  # 防止文字被覆盖
dev.off()
    ";
    open SAVE, ">$r_script"; print SAVE $script; close SAVE;
    system "$SOFT_R $r_script";
}


