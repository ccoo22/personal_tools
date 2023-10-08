# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
$|=1;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 检测 -> 脚本输入
my ($primer_file, $genome, $genome_dict, $output_dir, $is_with_primer, $if_help);
GetOptions(
    "primer|p=s"            => \$primer_file,
    "genome|ge=s"           => \$genome,
    "genome_dict|gt=s"      => \$genome_dict,
    "outdir|o=s"            => \$output_dir,
    "seq_with_primer|a!"    => \$is_with_primer,
    "help|h"                => \$if_help,
);
die help() if(defined $if_help or (not defined $primer_file or not defined $genome or not defined $genome_dict or not defined $output_dir or not defined $primer_file));

###################################################################### 主程序


my %hashPrimer = read_primer($primer_file);
my %hashGenome = read_fasta($genome, \%hashPrimer);

order_seq(\%hashPrimer);
print "\n[Notice] get fasta referrence with primer \n" if(defined $is_with_primer); # 参考序列添加引物

#输出
my $seq_fasta          = "$output_dir/seq.fa";
my $target_bed         = "$output_dir/target.bed";
my $target_realign_bed = "$output_dir/targetRealign.bed";

print "output : $seq_fasta\n";
print "output : $target_bed\n";
print "output : $target_realign_bed\n";

open FASTA, ">$seq_fasta";
open BED, ">$target_realign_bed";

foreach my $target(sort {$hashPrimer{$a}{"SortOrder"} <=> $hashPrimer{$b}{"SortOrder"}} keys %hashPrimer)
{
    my $chr  = $hashPrimer{$target}{'chr'};
    my $pos1 = $hashPrimer{$target}{'pos1'};
    my $pos2 = $hashPrimer{$target}{'pos2'};
    my $pos3 = $hashPrimer{$target}{'pos3'};
    my $pos4 = $hashPrimer{$target}{'pos4'};
    die "[error]参考基因组文件没有发现染色体编号 $chr，请确认引物文件是否写错了，或者选错参考基因组了!\n " if(not exists($hashGenome{$chr}));    

    my $start = $pos2;
    my $end   = $pos3;

    if(defined $is_with_primer)
    {   
        $start = $pos1;
        $end   = $pos4;
    }
    my $length = $end - $start + 1;

    my $seq = substr($hashGenome{$chr}, $start - 1, $length);
    my $seq_title = "$target|$chr:$start-$end";

    print FASTA ">$seq_title\n$seq\n";
    print BED "$chr\t$pos2\t$pos3\t+\t$target\n";

}
close FASTA;
close BED;
system "cat $genome_dict $target_realign_bed > $target_bed";

###############################################################子函数
# 排序
sub order_seq{
    my $hashPrimer = shift @_;

    my %hashTmp;
    foreach my $target(keys %$hashPrimer)
    {
        my ($chr, $pos1, $pos4) = ($hashPrimer->{$target}{'chr'}, $hashPrimer->{$target}{'pos1'}, $hashPrimer->{$target}{'pos4'});
        $hashTmp{$chr}{$pos1}{$pos4} = $target;
    }

    # 排序
    my $order = 0;
    foreach my $chr(sort keys %hashTmp)
    {
        foreach my $pos1(sort {$a <=> $b}keys %{$hashTmp{$chr}})
        {
            foreach my $pos4(sort {$a <=> $b}keys %{$hashTmp{$chr}{$pos1}})
            {
                $order++;
                my $target = $hashTmp{$chr}{$pos1}{$pos4};
                $hashPrimer{$target}{"SortOrder"} = $order;
            }            
        }
    }
}

# 读取引物信息
sub read_primer{
    my $primer_file = shift @_;

    my %hashPrimer;
    open PRIMER, $primer_file;
    my $line1 = <PRIMER>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;

    my @lost_heads = grep{ not $_ ~~ @heads } ('chr', 'pos1', 'pos2', 'pos3', 'pos4');
    die "[Error] primer 文件缺失表头： @lost_heads\n" if (scalar(@lost_heads) > 0);

    while(<PRIMER>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        my @datas = split /\t/, $_;
        my $title = $datas[0];
        die "[Error] duplicate $title\n" if(exists $hashPrimer{$title});
        foreach my $col(0..$#heads)
        {   
            my $value = $datas[$col];
            $value=~s/^\s+//;
            $value=~s/\s+$//;
            $value = "chr$value"  if($heads[$col] eq 'chr' and $genome=~/hg38_modify.fa$/ and $value!~/^chr/);  # hg38 添加chr符号
            $hashPrimer{$title}{$heads[$col]} = $value;
        }
        my @positions = map{$hashPrimer{$title}{$_}} ('pos1', 'pos2', 'pos3', 'pos4');
        my @positions_sort = sort {$a<=>$b} @positions;

        foreach my $col(0..$#positions)
        {
        die "[Error] $title 信息异常。输入的引物信息文件中，pos1/pos2/pos3/pos4  四列信息必须按照从小到大排列！\n" if($positions[$col] != $positions_sort[$col])
        }
    }
    close PRIMER;

    return %hashPrimer;
}

# 读取fasta文件
sub read_fasta{
    my $fasta_file = shift @_;
    my $hashPrimer = shift @_;

    # 需要的染色体
    my %hashChr;
    map{$hashChr{$hashPrimer->{$_}{'chr'}}++; } keys %$hashPrimer;

    # 读取
    print "Read $fasta_file ... ";
    my %hashFasta;
    my $fastaIn     = Bio::SeqIO->new(-file => $fasta_file, -format=>'Fasta') or die "Could not open up file $fasta_file: $!";
    while(my $inSeq = $fastaIn->next_seq)
    {
        my $id  = $inSeq->display_id;

        print "$id,";
        next if(not exists $hashChr{$id});

        my $seq = $inSeq->seq;
        $hashFasta{$id} = $seq;
    }
    print " OK\n";

    return %hashFasta;
}

sub help {
    my $info = "
Program: pcr_primer2ref， 根据PCR引物提取seq.fa文件信息
Version: 2019-06-06
Contact: 320 朱喜

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --primer/-p               [必填] 引物文件。使用在线引物匹配工具primer_position_PCR获取。 
                                          该文件第一行是表头，第一列必须是序列名称，其他列的顺序随意，但必须含有以下5列信息：chr pos1 pos2 pos3 pos4
         --genome/-ge              [必填] 样本物种对应的genome的绝对路径  
                                          示例：/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.fa 
                                               /home/genesky/database_new/self_build_database/ucsc/hg38_rm_alt_genome/hg38_modify.fa
         --genome_dict/-gt         [必填] 样本物种对应的genome_dict的绝对路径  
                                          示例：/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.dict  
                                               /home/genesky/database_new/self_build_database/ucsc/hg38_rm_alt_genome/hg38_modify.dict
         --outdir/-o               [必填] 输出目录  
         --seq_with_primer/-a      [选填] 参考序列是否添加引物
         --help/-h                  查看帮助文档。
    \n";
    return $info;
}
 
