$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
从一条序列里，提取中间特定区域的序列。同时分析序列的变异。
Version: v1.0 2020-09-25
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Bio::SearchIO;
use Cwd 'abs_path';

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 变量
my $DEFAULT_MODEL            = "blast";
my $DEFAULT_WORD_SIZE        = 11;
my $DEFAULT_NUM_THREADS      = 4;
my $DEFAULT_SOFT_CUTADAPT    = "/home/genesky/software/python/2.7.13/bin/cutadapt";
my $DEFAULT_SOFT_BLASTN      = "/home/genesky/software/blast+/2.9.0/bin/blastn";
my $DEFAULT_SOFT_MAKEBLASTDB = "/home/genesky/software/blast+/2.9.0/bin/makeblastdb";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($fasta, $seq_left, $seq_middle, $seq_right, $output_dir, $prefix, $model, $SOFT_CUTADAPT, $SOFT_BLASTN, $SOFT_MAKEBLASTDB, $word_size, $num_threads, $if_help);
GetOptions(
    "fasta|f=s"            => \$fasta,
    "seq_left|l=s"         => \$seq_left,
    "seq_middle|m=s"       => \$seq_middle,
    "seq_right|r=s"        => \$seq_right,
    "output_dir|o=s"       => \$output_dir,
    "prefix|p=s"           => \$prefix,

    "model=s"              => \$model,
    "cutadapt=s"           => \$SOFT_CUTADAPT,
    "blastn=s"             => \$SOFT_BLASTN,
    "makeblastdb=s"        => \$SOFT_MAKEBLASTDB,
    "word_size=s"          => \$word_size,
    "num_threads=s"        => \$num_threads,
    "help|h"               => \$if_help,
);
die "
Options: 必填
        --fasta/-f                输入fasta文件
                                  必须保证所有序列都是+向： 即与 seq_left.seq_middle.seq_right 构成的序列方向一致（--model blast 可以不考虑方向）。
                                  否则fasta序列中的seq_middle序列提取会失败

        --seq_left/-l             模版序列中，目标区域序列的左侧序列。 不能缺失
        --seq_middle/-m           模版序列中，目标区域序列。 不能缺失
        --seq_right/-r            模版序列中，目标区域序列的右侧序列。 不能缺失

        --output_dir/-o           结果输出目录, 脚本可以自己创建
        --prefix/-p               结果文件前缀， 通常是样本名

Options: 可选
        --model                   分析模型： blast/cutadapt   (default: '$DEFAULT_MODEL')
                                  blast模型：用 seq_left+seq_middle+seq_right 建库，然后blast， 提取中间N比对的序列
                                       优点：识别准确
                                       缺点：稍慢，识别过程稍麻烦。 
                                       注意：当三段序列总长度很小的时候，适当调节word_size参数。只有该模式才会做snv indel分析。


                                  cutadapt: 用cutadapt软件，以去接头的方式切掉头、尾的序列。
                                       有点：快速、方便、易于理解
                                       缺点：reads中的seq_left、seq_right序列靠近seq_middle的边缘序列不能有del碱基（3bp以内），否则会导致seq_middle序列边缘部分被多切掉一些，切掉的长度等于del的碱基数量。原因猜测：软件认为这部分边缘匹配失败是由于测序错误导致的，所以一起切掉。
                                            


        --cutadapt                更改软件 cutadapt 版本 (default: '$DEFAULT_SOFT_CUTADAPT')
        --blastn                  更改软件 blastn 版本 (default: '$DEFAULT_SOFT_BLASTN')
        --makeblastdb             更改软件 makeblastdb 版本 (default: '$DEFAULT_SOFT_MAKEBLASTDB')
        --word_size               调整blastn比对参数 word_size (default: '$DEFAULT_WORD_SIZE')
        --num_threads             调整blastn比对参数 num_threads 线程数 (default: '$DEFAULT_NUM_THREADS')
        --help/-h                 查看帮助文档

\n" if (defined $if_help or not defined $fasta or not defined $seq_left or not defined $seq_middle or not defined $seq_right or not defined $output_dir or not defined $prefix );
$model           = $DEFAULT_MODEL if (not defined $model);


$SOFT_CUTADAPT      = $DEFAULT_SOFT_CUTADAPT if (not defined $SOFT_CUTADAPT);
$SOFT_BLASTN        = $DEFAULT_SOFT_BLASTN if (not defined $SOFT_BLASTN);
$SOFT_MAKEBLASTDB   = $DEFAULT_SOFT_MAKEBLASTDB if (not defined $SOFT_MAKEBLASTDB);
$word_size          = $DEFAULT_WORD_SIZE if (not defined $word_size);
$num_threads        = $DEFAULT_NUM_THREADS if (not defined $num_threads);

$output_dir = abs_path($output_dir);
$fasta      = abs_path($fasta);
###################################################################### 主程序
make_dir($output_dir);


my $final_fa = "$output_dir/$prefix.final.fa";

# 使用trim adapter的方式去掉左右测序列
if($model eq 'cutadapt')
{
    print "[process] remove left/right seq\n";
    my $fasta_rm_right = "$output_dir/$prefix.rm_right.fa";
    my $fasta_rm_lr = "$output_dir/$prefix.rm_right_left.fa";
    system("$SOFT_CUTADAPT -f fasta  -a $seq_right -o $fasta_rm_right --times 1 --discard-untrimmed  $fasta > $fasta_rm_right.log");   # 没有接头的reads去掉
    system("$SOFT_CUTADAPT -f fasta  -g $seq_left  -o $fasta_rm_lr    --times 1 --discard-untrimmed  $fasta_rm_right > $fasta_rm_lr.log");   # 没有接头的reads去掉
    
    system("ln -s $fasta_rm_lr $final_fa");
}

# 使用blast的方式提取中间序列
if($model eq 'blast')
{
    # (1) 建库
    print "[process] build bastn db\n";
    my $blast_db = "$output_dir/$prefix.blastdb.fa";
    my $indelN   = 'N' x length($seq_middle);# N填充,放弃了，原因：当left/right边缘存在del时，会导致部分N区域的序列匹配到left/right区域
    open DB,">$blast_db";
    print DB ">db\n";
    print DB "$seq_left$seq_middle$seq_right\n";
    close DB;
    
    system "$SOFT_MAKEBLASTDB -in $blast_db -dbtype nucl";
    
    # (2) 比对
    print "[process] blast\n";
    my $blast_out = "$output_dir/$prefix.blast.txt";
    system "$SOFT_BLASTN -task blastn -query $fasta -evalue 0.00001 -db $blast_db -out $blast_out -num_threads $num_threads -gapopen 2 -gapextend 2 -word_size $word_size";  # 调整gap值，否则N太长，导致匹配断开。

    # (3) 序列识别提取
    print "[process] parse blast\n";
    my $blast_parse = "$output_dir/$prefix.blast.parse.fa";
    parse_blast($blast_out, $blast_parse);
 
    system("ln -sf $blast_parse $final_fa");
}
print "[process] final result: $final_fa\n";

# 数量统计
my $reads_count_raw = `grep '>' $fasta|wc -l`;
my $reads_count_ok  = `grep '>' $final_fa|wc -l`;
$reads_count_raw =~s/[\r\n]//g;
$reads_count_ok =~s/[\r\n]//g;
my $ratio = $reads_count_ok / $reads_count_raw;

open STATS, ">$output_dir/$prefix.final.stats";
print STATS "reads_raw\t$reads_count_raw\n";
print STATS "reads_find_middle\t$reads_count_ok\n";
print STATS "ratio\t$ratio\n";
close STATS;
###################################################################### 子程序


sub parse_blast{
    my $blast_out      = shift @_;
    my $blast_parse    = shift @_;
    
    # 目标起止位置，为seq_middle起始的前一个位置，终止的后一个位置
    my $start = length($seq_left);
    my $end   = $start + length($seq_middle) + 1;

    open PARSE, ">$blast_parse";
    open PARSEDETAIL, "|gzip >$blast_parse.detail.gz";
    open PARSESNVINDEL, "|gzip >$blast_parse.snv_indel.gz";
    print PARSESNVINDEL "SeqName\ttype\tinsertion\tdeletion\tsnp\trefseq\tcutseq\n";
    my $blastin = Bio::SearchIO->new(-format => 'blast',-file => $blast_out);
    my %hashOK;
    my %hashReads;
    while( my $r = $blastin->next_result ) {
        while( my $h = $r->next_hit ) {  
            while( my $hsp = $h->next_hsp ) {  
                my $hit_name   = $h->name; # 参照序列名称
                my $query_name = $r->query_name ; # fastq序列名称
                $hashReads{$query_name}++;
                last if(exists $hashOK{$query_name} );  # 识别过的reads不再使用

                my $hit_string = $hsp->hit_string;
                   $hit_string = uc($hit_string);# 参照序列
                my $qerry_string = $hsp-> query_string;
                   $qerry_string = uc($qerry_string);# fastq序列 

                my ($hit_start, $hit_end) = $hsp->range('hit') ;# 永远是从小到大的顺序，即使是Minus
                my $isok = 0;
                   $isok = 1 if($start >= $hit_start and $end <= $hit_end); # 比对区域包含中间部分
                next if($isok == 0);
                
                my $fangxiang = $hsp->strand('hit') ;# 获取参照序列的方向，查询序列永远是正向的
                if($fangxiang == -1){ # 如果参照序列是反向匹配的，则需要反转
                   $qerry_string =reverse $qerry_string;
                   $qerry_string =~tr/ATCGatcg/TAGCtagc/;
                   $hit_string   =reverse $hit_string; 
                   $hit_string   =~tr/ATCGatcg/TAGCtagc/;                                        
                }
                my @hit_strings   = split //,$hit_string;
                my @qerry_strings = split //,$qerry_string;
                my $position      = $hit_start;# 当前hit_string位置
                my $recordmark    = 0;# 是否可以提取
                my $markseq       = "";# 提取序列
                my $refseq        = "";# 参考序列
                foreach my $num(0..$#hit_strings)
                {
                    my $hit_base   = $hit_strings[$num];
                    my $query_base = $qerry_strings[$num];                    
                    $recordmark = ($position > $start and $position < $end) ? 1 : 0;
                    $position++ if($hit_base ne '-');# 下一个碱基的位置
                    # 目标区域
                    if($recordmark == 1)
                    {
                        $refseq.=$hit_base;
                        $markseq.=$query_base;
                    }
                    # 结尾前有插入缺失的情况，也算在middle中
                    if($hit_base eq '-' and $position eq $end)
                    {
                        $refseq.=$hit_base;
                        $markseq.=$query_base;
                    }
                }
                $hashOK{$query_name}++;

                # 用于后期debug, 检查是否识别错误
                print PARSEDETAIL "$query_name\tHitString\t$refseq\t$hit_string\n" ;       
                print PARSEDETAIL "$query_name\tQuerryString\t$markseq\t$qerry_string\n";   

                # snv_indel 识别
                my @snv_indels = call_snp_indel($refseq, $markseq);
                print PARSESNVINDEL "$query_name\t" . (join "\t", @snv_indels) . "\t$refseq\t$markseq\n";

                # fasta输出
                $markseq=~s/-//g;  # 清理
                print PARSE ">$query_name\n";
                print PARSE "$markseq\n";


        
            }
        }
    }
    close PARSE;
    close PARSEDETAIL;
    close PARSESNVINDEL;
    
    open PARSEFAILED, "| gzip >$blast_parse.failed.gz";
    foreach my $read_name(sort keys %hashReads)
    {
        next if(exists $hashOK{$read_name} );
        print PARSEFAILED "$read_name\n";
    }
    close PARSEFAILED;


} 



# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}

sub call_snp_indel{
    my $ref_seq = shift @_;
    my $alt_seq = shift @_;

    my @ref_bases = split //, $ref_seq;
    my @alt_bases = split //, $alt_seq;

    # SNP/INDEL识别
    my %hashInfo;

    my $pos_ref   = 0;  # 参考序列位置
    my %hashType;
    my %hashInsert;
    foreach my $pos(0..$#ref_bases)
    {
        my $ref_base = $ref_bases[$pos];  # 参考碱基
        my $alt_base = $alt_bases[$pos];  # 测序碱基
        $pos_ref++ if($ref_base =~ /\w/);  # 参考序列位置
        next if($alt_base eq $ref_base);
        
        if($ref_base=~/\w/ and $alt_base =~ /\w/)  
        {   # mismatch
            $hashInfo{'snp'} .= "$ref_base->$alt_base;";
            $hashType{$pos_ref}{'M'} = 1;
        }elsif($ref_base=~/\w/ and $alt_base !~ /\w/)  
        {   # delete
            $hashInfo{'deletion'} .= "$ref_base;";
            $hashType{$pos_ref}{'D'} = 1;
        }else{  
            # insert
            $hashInsert{$pos_ref} .= $alt_base;
            $hashType{$pos_ref}{'I'} = 1;
        }
    }
    # insertion信息汇总
    foreach my $pos_insert(sort {$a<=>$b} keys %hashInsert)
    {   
        $hashInfo{'insertion'} .= "$hashInsert{$pos_insert};";
    }
        
    # 突变位置信息汇总
    foreach my $pos_diff(sort {$a<=>$b} keys %hashType)
    {   
        foreach my $diff_type(sort keys %{$hashType{$pos_diff}} )
        {
            $hashInfo{'type'} .= $pos_diff.$diff_type.";";
        }
    }
    my @heads = qw(type insertion deletion snp);
    my @values = map{ my $value = exists $hashInfo{$_} ? $hashInfo{$_} : ""; $value } @heads;
    return @values;
}