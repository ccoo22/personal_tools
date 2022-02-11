# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

my $DEFAULT_SOFT_RSCRIPT = "/home/genesky/software/r/4.0.3/bin/Rscript";
my $DEFAULT_DB_RLIB = "/home/genesky/software/r/4.0.3/lib64/R/library";

# 检测 -> 脚本输入
my $ARGV_INFO = join " ", @ARGV;
my ($input_fasta, $output_dir, $prefix, $reads_file, $rscript, $rlib, $taxon, $cpu,  $if_help);
GetOptions(
    "input_fasta|i=s"    => \$input_fasta,
    "output_dir|o=s"     => \$output_dir,
    "prefix|p=s"          => \$prefix,
    "reads_file|r=s"     => \$reads_file,
    "rscript=s"          => \$rscript,
    "rlib=s"             => \$rlib,
    "help|h"             => \$if_help,
);

die "
将基因的DNA序列翻译（务必保证每一条序列能够被3整除），并统计密码子的数量、频率等信息，同时绘图

Options: 必填

        --input_fasta/-i            需要翻译的DNA序列
                                       当前仅支持密码子表 1
        --output_dir/-o             结果输出目录， 流程自动创建，不需要自己提前创建
        --prefix/-p                 输出文件的前缀，生成5个文件
                                        prefix.amino_acid.freq.txt  # 每个氨基酸的出现频率，用于绘图
                                        prefix.amino_acid.rscu1.txt  # RSCU大于1的密码子，用于绘图
                                        prefix.amino_acid.statistic.pdf  # 绘图
                                        prefix.amino_acid.statistic.txt  # 汇总统计表
                                        prefix.plot.r  # 绘图R脚本
Options: 可选
        --reads_file                一列数据，reads名称，限制只翻译input_fasta中指定的序列。
        --rscript                   rscript 软件 (default: $DEFAULT_SOFT_RSCRIPT)
        --rlib                      rscript 包数据库 （ default: $DEFAULT_DB_RLIB ）
        --help/-h                   查看帮助文档
\n" if (defined $if_help or not defined $input_fasta or not defined $output_dir or not defined $prefix);

$rscript  = $DEFAULT_SOFT_RSCRIPT if(not defined $rscript);
$rlib     = $DEFAULT_DB_RLIB if(not defined $rlib);

system "mkdir -p $output_dir" if(not -e $output_dir);

# 运行参数
my $DATA_TIME = time_to_datetime(time);
my $RUN_INFO = "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
[SET] 需要分析的序列      : $input_fasta
[SET] 结果输出目录        ： $output_dir
[SET] 输出文件前缀        : $prefix
[SET] 软件 rscript      : $rscript
[SET] 软件 rscript包路径   : $rlib
\n";
 print $RUN_INFO;

##################  主流程  ##########################################

my %hashDNA      = read_fasta($input_fasta);
my %hashCode     = get_amino_code();
my %hashCodeAbbr = get_amino_code_abbreviation();

###### 指定需要分析的id
my @seq_ids;
if(defined $reads_file)
{
    open FILE, $reads_file;
    while(<FILE>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        my ($id) = split /\s+/, $_;
        push @seq_ids, $id;
    }
}else{
    @seq_ids = keys %hashDNA;
}

############# 计数
print "翻译、统计\n";
my %hashCount;
my $sum = 0;
foreach my $id(@seq_ids)
{   
    die "[error] not find seq $id\n" if(not exists $hashDNA{$id});
    my $seq_dna   = $hashDNA{$id};

    my $code_count = length($seq_dna) / 3;
    foreach my $i(0..($code_count - 1))
    {
        my $base3 = substr($seq_dna, $i * 3, 3);
        my $code   = $hashCode{$base3};
        my $code_abbr = ($code eq '_') ? 'TER' : $hashCodeAbbr{$code};  # 三字符表示的氨基酸名称
        $base3=~s/T/U/g;
        $hashCount{$code_abbr}{$base3}++;
        $sum++;
    }
}

#################### 输出
print "输出\n";

my  @code_abbrs = qw/Gly Ala Val Leu Ile Ser Thr Cys  Met Asp  Glu  Asn   Gln  Lys   Arg  Phe Tyr His  Trp    Pro TER   /;

open OUTPUT, ">$output_dir/$prefix.amino_acid.statistic.txt";
my @heads = ("Codon", "Number", "Amino acids", "Ratio of Codon", "RSCU", "Number of amino acid", "Ratio of amino acid");
print OUTPUT (join "\t", @heads) . "\n";

open AMINO, ">$output_dir/$prefix.amino_acid.freq.txt";  # 氨基酸频率
print AMINO "amino_acid\tFreq\n";

open RSCU1_CODE, ">$output_dir/$prefix.amino_acid.rscu1.txt";  # RSCU大于1的密码子
print RSCU1_CODE "nucleotide\trscu\tclass\n";
foreach my $code_abbr(@code_abbrs)
{   
    next if(not exists $hashCount{$code_abbr});
    my @base3_list   = sort keys %{$hashCount{$code_abbr}};  #当前氨基酸包括的密码子 
    my $this_code_sum = 0;  # 当前氨基酸总数
    map{ $this_code_sum += $hashCount{$code_abbr}{$_}} @base3_list;
    my $this_code_perc = sprintf "%.2f", 100 * $this_code_sum / $sum;
    print AMINO "$code_abbr\t$this_code_perc\n";

    foreach my $base3(@base3_list)
    {   
        my $perc = sprintf "%.2f", 100 * $hashCount{$code_abbr}{$base3} / $sum;
        my $rscu = sprintf "%.2f", scalar(@base3_list) * $hashCount{$code_abbr}{$base3} / $this_code_sum;
        print OUTPUT "$base3\t$hashCount{$code_abbr}{$base3}\t$code_abbr\t$perc%\t$rscu\t$this_code_sum\t$this_code_perc\n";
        if($rscu > 1)
        {
            my $base_last = substr($base3,2, 1);
            print RSCU1_CODE "$base3\t$rscu\t$base_last-ending\n";
        }
    }
}
close AMINO;
close RSCU1_CODE;
close OUTPUT;

################ 绘图
print "绘图\n";
open SCRIPT, ">$output_dir/$prefix.plot.r";
print SCRIPT "
.libPaths('$rlib');
library(data.table)
library(ggpubr)
pdf('$output_dir/$prefix.amino_acid.statistic.pdf', width = 10)
data = fread('$output_dir/$prefix.amino_acid.freq.txt')
p <- ggbarplot(data, 
          x='amino_acid', 
          y='Freq', 
          fill='#fe0000', 
          color = '#fe0000', 
          xlab='',
          ylab='',
          legend.title = '',
          width=0.5
          )  + 
          scale_y_continuous(expand = expansion(mult = c(0, 0.05), add = c(0, 0)))   
p
# 图2
# 颜色设定
color_design = c('#ff1919', '#1ab861', '#5b9bd4', '#ffff15', '#854555')
names(color_design) = c('A-ending','U-ending','C-ending','G-ending','total')

data = fread('$output_dir/$prefix.amino_acid.rscu1.txt')
count = as.data.frame(data[,.N,by=class])
count2 = rbind(count, c('total', nrow(data)))
count2\$N = as.integer(count2\$N)
p <- ggbarplot(count2, 
          x='class', 
          y='N', 
          fill='class', 
          color = 'class', 
          palette = color_design,
          label = TRUE, 
          lab.pos = 'out', 
          xlab='',
          ylab='',
          legend.title = '',
          width=0.5,
          legend='none'
          )  + 
          scale_y_discrete(expand = expansion(mult = c(0, 0.05), add = c(0, 0)))   
p
dev.off()

";
close SCRIPT;
system("chmod 755 $output_dir/$prefix.plot.r");
system("$rscript $output_dir/$prefix.plot.r");










################################################# 子函数

# 提取fasta序列
sub read_fasta{
    my $fasta = shift @_;
    my %hashFasta;
    print "read $fasta ... ";
    my $IN = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
    while(my $inSeq = $IN->next_seq)
    {
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;
        $hashFasta{$id}= $seq;
    }
    print "OK\n";
    return %hashFasta;
}

sub get_amino_code{
    my(%genetic_code) = ( 
    'TCA' => 'S',    # Serine 
    'TCC' => 'S',    # Serine 
    'TCG' => 'S',    # Serine 
    'TCT' => 'S',    # Serine 
    'TTC' => 'F',    # Phenylalanine 
    'TTT' => 'F',    # Phenylalanine 
    'TTA' => 'L',    # Leucine 
    'TTG' => 'L',    # Leucine 
    'TAC' => 'Y',    # Tyrosine  
    'TAT' => 'Y',    # Tyrosine 
    'TAA' => '_',    # Stop 
    'TAG' => '_',    # Stop 
    'TGC' => 'C',    # Cysteine 
    'TGT' => 'C',    # Cysteine 
    'TGA' => '_',    # Stop 
    'TGG' => 'W',    # Tryptophan 
    'CTA' => 'L',    # Leucine 
    'CTC' => 'L',    # Leucine 
    'CTG' => 'L',    # Leucine 
    'CTT' => 'L',    # Leucine 
    'CCA' => 'P',    # Proline 
    'CCC' => 'P',    # Proline 
    'CCG' => 'P',    # Proline 
    'CCT' => 'P',    # Proline 
    'CAC' => 'H',    # Histidine 
    'CAT' => 'H',    # Histidine 
    'CAA' => 'Q',    # Glutamine 
    'CAG' => 'Q',    # Glutamine 
    'CGA' => 'R',    # Arginine 
    'CGC' => 'R',    # Arginine 
    'CGG' => 'R',    # Arginine 
    'CGT' => 'R',    # Arginine 
    'ATA' => 'I',    # Isoleucine 
    'ATC' => 'I',    # Isoleucine 
    'ATT' => 'I',    # Isoleucine 
    'ATG' => 'M',    # Methionine 
    'ACA' => 'T',    # Threonine 
    'ACC' => 'T',    # Threonine 
    'ACG' => 'T',    # Threonine 
    'ACT' => 'T',    # Threonine 
    'AAC' => 'N',    # Asparagine 
    'AAT' => 'N',    # Asparagine 
    'AAA' => 'K',    # Lysine 
    'AAG' => 'K',    # Lysine 
    'AGC' => 'S',    # Serine 
    'AGT' => 'S',    # Serine 
    'AGA' => 'R',    # Arginine 
    'AGG' => 'R',    # Arginine 
    'GTA' => 'V',    # Valine 
    'GTC' => 'V',    # Valine 
    'GTG' => 'V',    # Valine 
    'GTT' => 'V',    # Valine 
    'GCA' => 'A',    # Alanine 
    'GCC' => 'A',    # Alanine 
    'GCG' => 'A',    # Alanine 
    'GCT' => 'A',    # Alanine     
    'GAC' => 'D',    # Aspartic Acid 
    'GAT' => 'D',    # Aspartic Acid 
    'GAA' => 'E',    # Glutamic Acid 
    'GAG' => 'E',    # Glutamic Acid 
    'GGA' => 'G',    # Glycine 
    'GGC' => 'G',    # Glycine 
    'GGG' => 'G',    # Glycine 
    'GGT' => 'G',    # Glycine 
    ); 
    return %genetic_code;
}
sub get_amino_code_abbreviation{
    my %genetic_code_abbr = ('A' => 'Ala', 
                'R' => 'Arg', 
                'N' => 'Asn', 
                'D' => 'Asp', 
                'C' => 'Cys', 
                'Q' => 'Gln', 
                'E' => 'Glu', 
                'G' => 'Gly', 
                'H' => 'His', 
                'I' => 'Ile', 
                'L' => 'Leu', 
                'K' => 'Lys', 
                'M' => 'Met', 
                'F' => 'Phe', 
                'P' => 'Pro', 
                'S' => 'Ser', 
                'T' => 'Thr', 
                'W' => 'Trp', 
                'Y' => 'Tyr', 
                'V' => 'Val'
               );
    return %genetic_code_abbr;
}



sub time_to_datetime{
    my $time = shift;
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
    $year += 1900;
    $mon += 1;
    return "$year-$mon-$day $hour:$min:$sec";
}
