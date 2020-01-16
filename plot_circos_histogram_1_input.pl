# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $circos      = "/home/ganb/soft/circos-0.69-5/bin/circos";
 

# 检测 -> 脚本输入
my ($input, $output_dir, $prefix, $karyotype, $color, $max, $if_help);
GetOptions(
    "input|i=s"         => \$input,
    "output_dir|o=s"    => \$output_dir,
    "prefix|p=s"        => \$prefix,
    "karyotype|k=s"     => \$karyotype,
    "max|m=s"           => \$max,
    "color|c=s"         => \$color,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $input or not defined $output_dir or not defined $prefix or not defined $karyotype));
###################################################################### 主程序
$color = 'vdred' if(not defined $color);
$max   = 'max'   if(not defined $max);


# 准备绘图
my $config  = "$output_dir/$prefix.config.ini";
open CONFIG,">$config";
print CONFIG "karyotype = $karyotype\n";# 基因组信息
    
# 染色体参数
print CONFIG "<ideogram> \n";
print CONFIG "\t<spacing> \n";
print CONFIG "\t\tdefault = 0.005r \n";# 设置圈图中染色体之间的空隙大小，以下设置为每个空隙大小为周长的 0.5% 
print CONFIG "\t</spacing> \n";
print CONFIG "\tradius            = 0.85r \n";# 设定 ideograms 的位置，以下设定 ideograms 在图离圆心的 85% 处 
print CONFIG "\tthickness         = 60p   \n";# 设定 ideograms 的厚度，可以使用 r（比例关系） 或 p（像素）作为单位 
print CONFIG "\tstroke_color      = black \n";# 设定 ideograms 轮廓的颜色及其厚度。如果没有该参数或设定其厚度为0，则表示没有轮廓
print CONFIG "\tstroke_thickness  = 2p    \n";# 设定 ideograms 轮廓的颜色及其厚度。如果没有该参数或设定其厚度为0，则表示没有轮廓
print CONFIG "\tshow_label        = yes   \n";# 设定是否显示 label 。 label 对应着 karyotype 文件的第 4 列
print CONFIG "\tlabel_font        = bold  \n";# label字体
print CONFIG "\tlabel_radius      = 1r + 130p \n";# 设定 label 的位置 
print CONFIG "\tlabel_size        = 60    \n";# 设定 label 的大小 
print CONFIG "\tlabel_parallel    = yes   \n";# 设定 label 的字体方向，yes 是易于浏览的方向。
print CONFIG "\tshow_bands        = yes   \n"; # 显示条带
print CONFIG "\tfill_bands        = yes   \n"; # 填充条带
print CONFIG "\tband_transparency = 0     \n"; # 条带透明度，0-100
print CONFIG "</ideogram>   \n"; 

# 设置染色体tick
print CONFIG "chromosomes_units = 1000000\n";# 单位，1u=chromosomes_units
print CONFIG "show_ticks        = yes \n";# 显示坐标刻度
print CONFIG "show_tick_labels  = yes \n";# 显示坐标刻度标签
print CONFIG "<ticks>\n";                 # 定义tick全局参数
print CONFIG "\tradius          = 1r \n"; # 刻度半径
print CONFIG "\tcolor           = black\n";# 刻度颜色
print CONFIG "\tthickness       = 2p \n";  # 刻度线条宽度
print CONFIG "\tmultiplier      = 1e-6 \n";# 
print CONFIG "\tformat          = \%0.3f \n"; # 数值整型
print CONFIG "\t<tick> \n";# 第二刻度
print CONFIG "\t\tspacing        = 30u \n";# 刻度间隔
print CONFIG "\t\tsize           = 15p \n";# 线条长度
print CONFIG "\t\tshow_label     = yes\n"; # 是否显示刻度所代表的大小
print CONFIG "\t\tlabel_font     = bold \n";# 刻度标签大小
print CONFIG "\t\tlabel_size     = 40p \n";# 刻度标签大小
print CONFIG "\t\tlabel_offset   = 10p \n";#  
print CONFIG "\t\tformat         = \%d \n";# 数值整型 
print CONFIG "\t</tick> \n"; 
print CONFIG "</ticks> \n\n"; 

# 开始绘图
print CONFIG "<plots>  \n";
print CONFIG "\t<plot>   \n";
print CONFIG "\t\ttype        = histogram    \n";
print CONFIG "\t\tfile        = $input   \n";
print CONFIG "\t\torientation = out   \n";# 方向向外
print CONFIG "\t\textend_bin  = no   \n"; # 区域不自动连接
print CONFIG "\t\tthickness   = 2p   \n"; # 边框宽度
print CONFIG "\t\tr0          = 0.65r   \n"; # 直方图绘制范围
print CONFIG "\t\tr1          = 0.90r   \n";
print CONFIG "\t\tfill_color  = $color   \n";# 填充色
print CONFIG "\t\tcolor       = $color   \n";# 边框色
print CONFIG "\t\tmin         = 0   \n";# 坐标最小值
print CONFIG "\t\tmax         = $max   \n" if($max ne 'max');# 最表最大值
print CONFIG "\t\t<axes>   \n"; # 显示梯度
print CONFIG "\t\t\tshow      = data  \n"; 
print CONFIG "\t\t\tthickness = 1  \n"; 
print CONFIG "\t\t\tcolor     = grey  \n"; 
print CONFIG "\t\t\t<axis>  \n"; 
print CONFIG "\t\t\t\tspacing = 0.2r   \n"; 
print CONFIG "\t\t\t</axis>    \n"; 
print CONFIG "\t\t</axes>   \n"; 
print CONFIG "\t</plot>    \n";# 最表最大值    
    
# plots end
print CONFIG "</plots>  \n";    
# 绘图
print CONFIG "<image> \n"; 
print CONFIG "\t<<include etc/image.conf>>  \n"; 
print CONFIG "</image>  \n"; 
print CONFIG "<<include etc/colors_fonts_patterns.conf>>  \n"; 
print CONFIG "<<include etc/housekeeping.conf>>  \n"; # 注意：有时候，输入start-end区域过长（超过了默认限制25000），会报错，此时应该把该文件拷贝过来修改，并指定。/home/ganb/soft/circos-0.69-5/etc/housekeeping.conf
close CONFIG;
system("$circos -conf $config -outputdir $output_dir -outputfile $prefix");  






###################################################################### 子函数

sub help{
    my $info = "
Program: circos 直方图绘制
Version: 2019-11-27
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
         必填：
         --input/-i         要绘制的直方图文件,要符合circos直方图文件格式
         --output_dir/-o    输出结果目录
         --prefix/-p        结果文件前缀
         --karyotype/-k     circos格式的染色体信息，例如：/home/ganb/soft/circos-0.69-5/data/karyotype/karyotype.human.hg38.txt
         
         选填：
         --color/-c         直方图填充颜色，circos颜色模版名称，默认：vdred
         --max/-m           显示时的最大值，必须是数字。默认：最大值
         --help/-h          查看帮助文档
    \n";
    return $info;
}

