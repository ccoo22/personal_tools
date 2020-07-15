$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
样本亲缘关系 分析
Version: v1.0 2020-04-24
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 变量
my $DEFAULT_SOFT_PLINK    = "/home/genesky/software/plink/1.07/plink";
my $DEFAULT_SOFT_KING     = "/home/genesky/software/king/2.2.5/king";
my $DEFAULT_SOFT_RSCRIPT  = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $DEFAULT_SOFT_RLIB     = "/home/genesky/software/r/3.5.1/lib64/R/library";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input_plink, $output_dir, $prefix, $SOFT_PLINK, $SOFT_KING, $SOFT_RSCRIPT, $SOFT_RLIB, $if_help);
GetOptions(
	"input_plink|i=s"     => \$input_plink,
	"output_dir|o=s"      => \$output_dir,
	"prefix|p=s"          => \$prefix,

	"plink=s"             => \$SOFT_PLINK,
	"king=s"              => \$SOFT_KING,
	"rscript=s"           => \$SOFT_RSCRIPT,
	"rlib=s"              => \$SOFT_RLIB,
	"help|h"              => \$if_help,
);
die "
Options: 必填

        --input_plink/-i                原始plink格式数据输入前缀， 例如sample.ped, sample.map文件， 则输入 -i sample
        --prefix/-p                     输出文件前缀， 例如 a
        --output_dir/-o                 结果输出路径

Options: 可选

        --plink                    更改软件 plink 版本 (default: $DEFAULT_SOFT_PLINK)
        --king                     更改软件 king 版本 (default: $DEFAULT_SOFT_KING)
        --rscript                  更改软件 rscript 版本 (default: $DEFAULT_SOFT_RSCRIPT)
        --rlib                     更改软件 rlib 版本 (default: '$DEFAULT_SOFT_RLIB')
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $input_plink or not defined $prefix or not defined $output_dir);


$SOFT_PLINK    = $DEFAULT_SOFT_PLINK if (not defined $SOFT_PLINK);
$SOFT_KING     = $DEFAULT_SOFT_KING if (not defined $SOFT_KING);
$SOFT_RSCRIPT  = $DEFAULT_SOFT_RSCRIPT if (not defined $SOFT_RSCRIPT);
$SOFT_RLIB     = $DEFAULT_SOFT_RLIB if (not defined $SOFT_RLIB);
 
$output_dir = File::Spec->rel2abs($output_dir);  # 转换为绝对路径
mkdir $output_dir if(not -e $output_dir);

###################################################################### 初始化
my $DATA_TIME = `date +\"\%Y-\%m-\%d \%H:\%M.\%S\"`;
my $RUN_INFO = "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
[SET] 软件 plink : $SOFT_PLINK
[SET] 软件 king : $SOFT_KING
[SET] 软件 rscript : $SOFT_RSCRIPT
[SET] 软件 rlib : $SOFT_RLIB
";
open SAVE, ">>$output_dir/kinship_run.info"; print SAVE $SCRIPT_INFO.$RUN_INFO; close SAVE;
print $RUN_INFO;

###################################################################### 主程序
my $ped_input = "$input_plink.ped";
my $map_input = "$input_plink.map";
die "[Error] lost $ped_input, $map_input file\n" if(is_file_ok($ped_input, $map_input) == 0);

# (1) 去掉家系信息
print "(1) 去掉家系信息\n";
my $clean_plink = "$output_dir/$prefix.clean";
clean_family_info($input_plink, $clean_plink); # 去掉家系信息

# （2）转二进制
print "(2) 转二进制\n";
my $clean_plink_b = "$output_dir/$prefix.clean.binary";
system("$SOFT_PLINK --file $clean_plink --make-bed --out $clean_plink_b --noweb"); # 二进制转化

# (3) 亲缘关系计算
print "(3) 亲缘关系计算\n";
my $relation = "$output_dir/$prefix.relationship";
system("$SOFT_KING -b $clean_plink_b.bed  --kinship --prefix $relation");# 关系计算 

# (4) 结果汇总、绘图
print "(4) 结果汇总、绘图\n";
my $relation_raw = "$relation.kin0";
my $relation_matrix = "$relation_raw.matrix";
my $relation_pdf    = "$output_dir/$prefix.relationship.pdf";
relation_format($relation_raw, $relation_matrix);
plot($relation_matrix, $relation_pdf);

print "[final] $relation_raw\n";
print "[final] $relation_pdf\n";

###################################################################### 子程序

sub plot{
    my $relation_matrix = shift @_;
    my $relation_pdf    = shift @_;

    # 计算样本数量
    my $sample_count_raw = `wc -l $relation_matrix`;
    my ($sample_count) = $sample_count_raw =~/^(\d+)/;
    $sample_count = $sample_count - 1;

    # pdf宽度设置
    my $pdf_width = ($sample_count > 10) ? $sample_count : 10;
 

    my $r_command = "";
    foreach my $command( (".libPaths('$SOFT_RLIB')",
                         "relation = read.table('$relation_matrix', sep = '\\t', header = T, row.name = 1)",
                         "pdf('$relation_pdf', width = $pdf_width, height = $pdf_width)",
                         "panel.cor <- function(x, y, digits=2, prefix='', cex.cor, ...){usr <- par('usr'); on.exit(par(usr)); par(usr = c(0, 1, 0, 1)); z=x[!is.na(y)]; txt=as.numeric( sprintf( '\%0.4f', z[length(z)] ) ); if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt);color=1;if(txt>=0.354) color=2;if(txt>=0.177 && txt<0.354) color=3;if(txt>=0.0884 && txt<0.177) color=4;if(txt>=0.0442 && txt<0.0884) color=5; text(0.5, 0.5, txt, cex = cex.cor,col=color);}",
                         "panel.text <- function(x, y, labels, cex, font, ...){cex.cor <- 0.8/strwidth(labels);text(0.5, 0.5, labels, cex = cex.cor);}",
                         "pairs(relation ,lower.panel = NULL ,upper.panel = panel.cor,text.panel=panel.text ,font.labels = 2, main='Sample Relationship (Based On King software)')",
                         "info = c('>0.354                = duplicate/MZ twin\n\n','[0.177, 0.354]     = 1st-degree\n\n','[0.0884, 0.177]   = 2nd-degree\n\n','[0.0442, 0.0884] = 3rd-degree\n\n')",
                         "textExpand = $pdf_width",
                         "if($sample_count < 10) textExpand = 10",
                         "mtext(info, side = 1, adj = 0, cex = 2*textExpand/10, line = c(-2*textExpand/10, 0, 2*textExpand/10, 4*textExpand/10), col = c(2,3,4,5))",
                         "dev.off()")
            )
    {
        $r_command .= "$command\n";
    }

    open SCRIPT, ">$relation_pdf.rscript";
    print SCRIPT $r_command;
    close SCRIPT;

    system("$SOFT_RSCRIPT $relation_pdf.rscript");
	
}

sub relation_format{
    my $input = shift @_;
    my $output = shift @_;

    my %hashRelation;  # 记录样本间亲缘关系
    my %hashSample;  # 记录样本ID

    open INPUT, $input;
    my $line1 = <INPUT>;
       $line1=~s/[\r\n]//g;
    my @heads = split /\t/, $line1;

    my $sample_count = 0;
	while(<INPUT>)
    {
		$_=~s/[\r\n]//g;
		next if($_!~/\w/);
		my @datas = split /\s+/, $_;

		my %hashTmp;
		foreach my $col(0..$#heads)
        {
			my $value = ( exists $datas[$col] ) ? $datas[$col] : "";
			$hashTmp{$heads[$col]} = $value;
		}

        my $sample1 = $hashTmp{'ID1'};
        my $sample2 = $hashTmp{'ID2'};
        my $kinship = $hashTmp{'Kinship'};

        # 数据保存
        $hashRelation{$sample1}{$sample2} = $kinship;

        $sample_count++;
        $hashSample{$sample1} = $sample_count if(not exists $hashSample{$sample1});
        $sample_count++;
        $hashSample{$sample2} = $sample_count if(not exists $hashSample{$sample2});
	}
	close INPUT;    

    # 构建新矩阵，方便R读取
    
    my @samples = sort {$hashSample{$a} <=> $hashSample{$b}} keys %hashSample;

    open OUTPUT, ">$output";
    print OUTPUT "Sample\t" . (join "\t", @samples) . "\n";

    my $row = 0;
    foreach my $sample_row(@samples)
    {   
        my @values = ($sample_row);

        $row++;
        my $col = 0;
        
        foreach my $sample_col (@samples)
        {   
            $col++;
            my $kinship = 1;
               $kinship = $hashRelation{$sample_row}{$sample_col} if(exists $hashRelation{$sample_row} and exists $hashRelation{$sample_row}{$sample_col});
               $kinship = $hashRelation{$sample_col}{$sample_row} if(exists $hashRelation{$sample_col} and exists $hashRelation{$sample_col}{$sample_row});
            $kinship = "" if($row > $col);

            push @values, $kinship;
        }

        print OUTPUT (join "\t", @values) . "\n";
    }
    close OUTPUT;
}


sub clean_family_info{
    my $input        = shift @_;
    my $output       = shift @_;

    system("cp $input.map $output.map");

    open INPUT, "$input.ped";
    open OUTPUT, ">$output.ped";

    while(<INPUT>)
    {
        my ($fid, $iid, $pid, $mid, $sex, $pheno, $tmp) = split /\s+/, $_, 7;
        print OUTPUT "$iid\t$iid\t0\t0\t$sex\t$pheno\t$tmp";
    }
    close INPUT;
    close OUTPUT;
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

