# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use GD;
use Parallel::ForkManager;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"}; 

# 定义 -> 核心变量
 
# 检测 -> 脚本输入
my ($kegg_enrichment, $color_file, $output_dir, $keep_tmp, $add_address_to_excel, $no_download_png, $if_help);
GetOptions(
    "kegg_enrichment|k=s"     => \$kegg_enrichment,
    "output_dir|o=s"          => \$output_dir,
    "color|c=s"               => \$color_file,
    "keep_tmp!"               => \$keep_tmp,
    "add_address_to_excel!"   => \$add_address_to_excel,
    "no_download_png!"        => \$no_download_png,
    "help|h"                  => \$if_help,
);
die help() if(defined $if_help or (not defined $kegg_enrichment or not defined $output_dir));
$color_file = "" if(not defined $color_file);
###################################################################### 主程序
my $tmp_dir = "$output_dir/download_tmp";  # 临时目录，用后即删
mkdir $output_dir if(not -e $output_dir);  # 输出目录
mkdir $tmp_dir    if(not -e $tmp_dir);     # 临时目录

#############
# (1) 提取差异通路信息/颜色信息
#############
my %hashDiffPahtway = read_kegg_enrichment($kegg_enrichment, 0.05); # 读取富集分析结果,p值小于0.05。$hashDiffPahtway{$pathway_id}{$gene_name}
my %hashColor = read_color_file($color_file, \%hashDiffPahtway); # 读取基因上下调状态。 $hashColor{$gene} = lc($regular)

#############
# (2) 下载通路文件
#############
my @pathway_ids = sort keys %hashDiffPahtway;
download_pathway(\@pathway_ids, $tmp_dir, 'png')  if(not defined $no_download_png);
download_pathway(\@pathway_ids, $tmp_dir, 'conf') if(not defined $no_download_png);

#############
# (3) 标记颜色
#############
mark_color(\%hashColor, \@pathway_ids, $output_dir, $tmp_dir) if(not defined $no_download_png);

#############
# (4) 下载地址加入到excel表格里
#############
add_url_to_kegg_enrichment(\%hashDiffPahtway, \%hashColor, $kegg_enrichment) if(defined $add_address_to_excel);


#############
# (5) 清理临时目录
#############
system("rm -r $tmp_dir")    if(not defined $keep_tmp);
system("rm -r $output_dir") if(defined $no_download_png);

###################################################################### 子函数

sub add_url_to_kegg_enrichment{
    my $hashDiffPahtway = shift @_;
    my $hashColor       = shift @_;
    my $kegg_enrichment = shift @_;

    print "添加通路图下载地址到富集分析文件里\n";
    # 制作下载地址
    my %hashAdd;
    foreach my $pathway_id(sort keys %$hashDiffPahtway)
    {
        # my $add = "http://www.genome.jp/kegg-bin/show_pathway?$pathway_id"; # 下载地址
        my $add = "http://www.kegg.jp/kegg-bin/show_pathway?$pathway_id"; # 下载地址
        my @lost_genes;
        foreach my $gene_name(sort keys %{$hashDiffPahtway{$pathway_id}})
        {
            my $gene_color = 'red';
               $gene_color = 'red'  if(exists $hashColor{$gene_name} and $hashColor{$gene_name} eq 'up');
               $gene_color = 'blue' if(exists $hashColor{$gene_name} and $hashColor{$gene_name} eq 'down');
            my $add_tmp = "$add/$gene_name%09$gene_color";

            # kegg 提交的地址长度不能超过3040（大约），否则直接拒绝访问
            if(length($add_tmp) > 3040)
            {    
                push @lost_genes, $gene_name;
                next;
            }
            $add = $add_tmp;
 
        }
        # print "[Warnings] 下载地址中，从通路 [$pathway_id] 里删除了后面的部分基因，原因：kegg的访问地址字符长度不能超过3040左右，否则直接拒绝访问\n通路[$pathway_id] 删除的基因：@lost_genes\n" if(@lost_genes > 0);
        $hashAdd{$pathway_id} = $add;
    }
    # 下载地址加入到kegg_enrichment文件中

    open KEGG_ENRICHMENT, $kegg_enrichment;
    open KEGG_ENRICHMENT_ADD, ">$kegg_enrichment.add_address.txt";
    my $line1 = <KEGG_ENRICHMENT>;
          $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;
    print KEGG_ENRICHMENT_ADD "$line1\tkegg_link\n";  # 更新表头
    while(<KEGG_ENRICHMENT>)
    {
        next if($_!~/\w/);
        $_=~s/[\r\n]//g;
        my @datas = split /\t/, $_;
        my %hashTmp = map{  ($heads[$_], $datas[$_])  } (0..$#heads);

        my $pathway_id     = $hashTmp{'ID'};
        my $pathway_address = exists $hashAdd{$pathway_id} ? $hashAdd{$pathway_id} : " ";
        print KEGG_ENRICHMENT_ADD "$_\t$pathway_address\n";
    }
    close KEGG_ENRICHMENT;
    close KEGG_ENRICHMENT_ADD;
    system("mv $kegg_enrichment.add_address.txt $kegg_enrichment");
}

# 通路图标记颜色
sub mark_color{
    my $hashColor   = shift @_;
    my $pathway_ids = shift @_;
    my $output_dir  = shift @_;
    my $tmp_dir     = shift @_;

    print "开始标记颜色\n";
    foreach my $pathway_id(@$pathway_ids)
    {    
        # 原始下载文件
        my $png_raw   = "$tmp_dir/$pathway_id.png";
        my $conf_raw  = "$tmp_dir/$pathway_id.conf";

        # 标记颜色后的文件
        my $png_new   = "$output_dir/$pathway_id.png";
        my $html_new  = "$output_dir/$pathway_id.html";

        # 检查两个文件是否丢失
        if(not -e $conf_raw)
        {
            print "[Error] $pathway_id CONF文件丢失，无法标记颜色:$conf_raw\n";
            system("cp $png_raw $png_new") if(-e $png_raw);
            next;
        }
        if(not -e $png_raw)
        {
            print "[Error] $pathway_id png 文件丢失 : $png_raw\n";
            next;
        }

        # 开始标记
        my %hashConf = read_pathway_conf($conf_raw, $hashColor);
        my $image = GD::Image->newFromPng($png_raw);  # 读入png文件

        # 声明两种颜色
        my $red_color    = $image->colorAllocate(255, 0, 0);
        my $blue_color   = $image->colorAllocate(0, 0, 255);
        my $black_color  = $image->colorAllocate(0, 0, 0);

        # 创建html，供客户点击查看
        open HTML, ">$html_new";
        print HTML "<html>\n";
        print HTML "\t<head>\n";
        print HTML "\t\t<title>KEGG PATHWAY: $pathway_id</title>\n";
        print HTML "\t</head>\n";

        print HTML "\t<body>\n";
        print HTML "\t\t<img src='$pathway_id.png' name='pathwayimage' name='pathwayimage' usemap='#mapdata' border='0' />\n";
        print HTML "\t\t<map name='mapdata'>\n";

        # 开始处理通路图的每一个元素
        foreach my $row(sort {$a <=> $b} keys %hashConf)
        {
            my ($shape, $coords, $href, $title) = map{ $hashConf{$row}{$_} } ("shape", "coords", "href", "title");
            print HTML "\t\t<area shape=$shape coords=$coords href='$href' title='$title' target='_blank' />\n";  # 存储到HTML

            # 标记颜色
            if(exists $hashConf{$row}{'mark'})
            {	
                my ($x1, $y1, $x2, $y2) = split /,/, $coords;
                next if(not defined $y2);
                # 对目标区域每一个像素修改颜色
                foreach my $x($x1..$x2)
                {
                    foreach my $y($y1..$y2)
                    {
                        my ($R, $G, $B) = $image->rgb($image->getPixel($x, $y));
                        next if($R == 4 and $G == 2 and $B == 4); # 黑色不用替换，是基因名字
                        $image->setPixel($x, $y, $red_color)   if($hashConf{$row}{'mark'} eq 'up');   # 上调、红色
                        $image->setPixel($x, $y, $blue_color)  if($hashConf{$row}{'mark'} eq 'down');# 下调、蓝色
                        $image->setPixel($x, $y, $black_color) if($R == 252 and $G == 2 and $B == 4 and $hashConf{$row}{'mark'} eq 'up'); # 原始通路图上，部分基因名字是红色，故要把名字改成黑色
                    }
                }                
            }
        }
        print HTML "\t</body>\n";
        print HTML "</html>\n";
        close HTML;

        # PNG输出
        open PNG, ">$png_new";
        binmode PNG;
        print PNG $image->png;
        close PNG;
    }
    print "颜色标记完成\n";

}


# 读入通路conf文件
sub read_pathway_conf{
    my $file      = shift;
    my $hashColor = shift;
 
    open FILE, $file;
    my %hashConf;
    my $row = 0;
    while(my $line = <FILE>)
    {
        $line =~ s/[\r\n]//g;
        my @split_line = split /\t/, $line;
        # 形状
        my ($shape) = $split_line[0] =~ /^(\w+)/;
        # 坐标
        my @coords = $split_line[0] =~ /(\d+)/g;
        # 链接
        my $href = "https://www.kegg.jp/$split_line[1]";
        # 信息
        my $title = $split_line[2];
        %{$hashConf{$row}} = ("shape" => $shape, "coords" => (join ",", @coords), "href" => $href, "title" => $title);
        # 标记
        my @genes = $title =~ /\(([^\)]+)\)/g;
		
        my $color = "";
        foreach my $gene(@genes)
        {	
            $color = $hashColor->{$gene} if(exists $hashColor->{$gene} and $hashColor->{$gene} ne '');
        }
        $hashConf{$row}{"mark"} = $color if($shape eq "rect" and $color ne "");
        
        $row++;
    }
    close FILE;
    return %hashConf;
}


# 下载通路数据
sub download_pathway{
    my $pathway_ids = shift @_;
    my $tmp_dir     = shift @_;
    my $file_type   = shift @_;  # 文件类型

    # 下载状况记录
    my %hashCondition;
    map{ $hashCondition{'Error'}{$_} = 1;  $hashCondition{'ErrorCount'}++; } @$pathway_ids;
    
    # 下载（如果有文件下载失败，则重新循环）
    print "下载原始通路图的 $file_type 文件\n";
    my $time = 0;
    do
    {
        $time++;
        my $max_threads = 10;
        my $pm          = Parallel::ForkManager->new($max_threads);

        # 并行下载
        foreach my $pathway_id (sort keys %{$hashCondition{'Error'}}) 
        {
            my $pid = $pm->start and next;
            system("wget -q -O $tmp_dir/$pathway_id.png  http://rest.kegg.jp/get/$pathway_id/image")  if($file_type eq 'png');  # 通路png文件下载
            system("wget -q -O $tmp_dir/$pathway_id.conf http://rest.kegg.jp/get/$pathway_id/conf")   if($file_type eq 'conf'); # 通路conf文件下载
            $pm->finish;
        }
        $pm->wait_all_children;

        # 检查是否下载失败
        $hashCondition{'ErrorCount'} = 0;
        foreach my $pathway_id (sort keys %{$hashCondition{'Error'}}) 
        {    
            my $download_file = "";
               $download_file = "$tmp_dir/$pathway_id.png"  if($file_type eq 'png');
               $download_file = "$tmp_dir/$pathway_id.conf" if($file_type eq 'conf');
            # 检查文件是否完整
            if(($file_type eq 'png' and -e $download_file  and is_end_png($download_file) == 1 ) or ($file_type eq 'conf' and -e $download_file and -s $download_file > 0) )
            {
                delete $hashCondition{'Error'}{$pathway_id};
                next;
            }
            $hashCondition{'ErrorCount'}++;
        }

        # 提醒
        my @lost_hsa_ids = keys %{$hashCondition{'Error'}};
        print "[Warnings] kegg $file_type  文件下载循环 [第 $time 次]，有 [$hashCondition{'ErrorCount'}] 个通路下载失败，重新尝试下载失败通路(睡眠10s),如果反复失败，可手工协助下载或放弃该通路 \n失败通路列表：@lost_hsa_ids\n " if($hashCondition{'ErrorCount'} > 0);
        sleep 10;

    } while ($hashCondition{'ErrorCount'} != 0);     
    print "kegg $file_type  文件下载完成\n";
}

# png尾部是否有问题（存在png下载了一半，然后终止了的情况）
sub is_end_png
{
    my $png = shift;
    my $buff;
    my $str = "00000000";  # 补充前缀，防止下面提取数据失败，产生报错警告
    open PNG, $png or die "Can't open $png";
    binmode(PNG);
    while(read(PNG, $buff, 1)){
        my $hex = unpack("H*", $buff);
        $str .= $hex;
    }
    close PNG; 
    my $line = substr $str, -8;

    return $line  eq qq{ae426082} ? 1 : 0;
}

# 读入颜色标记信息
sub read_color_file{
    my $color_file      = shift @_;
    my $hashDiffPahtway = shift @_;
    my %hashColor;
    
    # 所有的基因默认设置为UP:红色
    foreach my $pahtway_id(keys %$hashDiffPahtway)
    {
        map{ $hashColor{$_} = 'up' } keys %{$hashDiffPahtway{$pahtway_id}};
    }

    return %hashColor if($color_file eq '');

    # 根据color文件设定颜色
    open COLOR, $color_file;
    while(<COLOR>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        my ($gene, $regular) = split /\t/, $_;
		my $color = (defined $regular and $regular ne '') ? $regular : "up";
        $hashColor{$gene} = lc($color);
    }
    close COLOR;
    return %hashColor;

}

# 读入通路信息
sub read_kegg_enrichment{
    my $kegg_enrichment = shift @_;
    my $pvalue_cutoff   = shift @_;

    print "提取要下载的通路信息\n";
    my %hashDiffPahtway;
    open ENRICHMENT, $kegg_enrichment;
    my $line1 = <ENRICHMENT>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;

    while(<ENRICHMENT>)
    {
        next if($_!~/\w/);
        $_=~s/[\r\n]//g;
        my @datas = split /\t/, $_;
        my %hashTmp = map{  ($heads[$_], $datas[$_])  } (0..$#heads);

        my $pathway_id     = $hashTmp{'ID'};
        my $pvalue         = $hashTmp{'pvalue'};
        my $gene_name_list = $hashTmp{'geneID'};

        next if($pvalue > $pvalue_cutoff);  # 只要满足阈值的
        foreach my $gene_name(split /\//, $gene_name_list)
        {
            $hashDiffPahtway{$pathway_id}{$gene_name}++;
        }
    }
    close ENRICHMENT;
    return %hashDiffPahtway;
}

sub help{
    my $info = "
Program: get kegg pathway download address， 根据富集分析结果，获取下载地址,基因设成红色
Version: 2019-05-10
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --kegg_enrichment/-k   [必填] kegg_enrichment.xls 结果文件。kegg通路富集分析原始结果
         --output_dir/-o        [必填] 输出目录

         --color/-c             [选填] 基因颜色文件，第一列基因，第二列up/down, (up=read, down=blue)
         --keep_tmp             [选填] 是否保留临时数据
         --add_address_to_excel [选填] 把通路图下载地址加入kegg_enrichment.xls文件里。
         --no_download_png      [选填] 是否下载kegg官方图，默认：下载 。注：如果服务器IP被禁止，可能会导致下载死循环。   
         --help/-h              查看帮助文档
    \n";
    return $info;
}
