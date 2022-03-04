use strict;
use warnings;
use Bio::SeqIO;
use Parallel::ForkManager;
$|=1;
# 当前序列在R1端
my @samples = qw(804_1 804_2 804_3 805_1 805_2 805_3 806_1 806_2 806_3 807_1 807_2 807_3 808_1 808_2 808_3 809_1 809_2 809_3 810_1 810_2 810_3 811_1 811_2 811_3 812_1 812_2 812_3);
my $fastq_dir = "/home/pub/project/genetics/21D1226A_ACC2027_PCR";
my $output_dir = "/home/pub/output2/Target/21D1226A/output_blue";
my $report_dir = "/home/pub/output2/Target/21D1226A/report_blue";
my $prefix = "ACC2027";
make_dir($output_dir);
make_dir($report_dir);

my $SOFT_BLASTN         = "/home/genesky/software/blast+/2.9.0/bin/blastn";
my $SOFT_MAKEBLASTDB    = "/home/genesky/software/blast+/2.9.0/bin/makeblastdb";
my $SOFT_FASTQ_TO_FASTA = "/home/genesky/software/fastx_toolkit/0.0.14/bin/fastq_to_fasta";
my $SOFT_MUSCLE         = "/home/genesky/software/muscle/3.8.31/muscle3.8.31_i86linux64";
my $SOFT_MAFFT          = "/home/genesky/software/mafft/7.487/core/mafft";
my $SOFT_MIDDLE_SEQ_SNV = "/home/pub/bin/NGS/chip/GATK4/tools/personal/fasta_cut_middle_seq_and_snv_indel_analysis.pl";
my $SOFT_TABLE2EXCEL    = "/home/pub/bin/NGS/chip/GATK4/tools/personal/table2excel.pl";

my $seq_left   = "TTCCTCGTGCTGGACAAGTGTGGTTTCCAGATTCTGCAACCAAGACTGCGCAGGCATTGCTGGACTTCAACCGTGAAGGATTACCTCTG";
my $seq_middle = "TTCATCCTCGCTAACTGGAG";
my $seq_right  = "AGGCTTCTCTGGTGGACAAAGAGATCTTTTTGAAGGAATTCTTCAGGCTGGCTCGACTATTGTTGAGAACCTTAGGACATACAATCAGCCTGCC";
my $seq_ref    = "AGACAGTGGTTACTGGCAGAGCAAAGCTTGGTGGAATTCCAGTGGGTGTGATAGCTGTGGAGACTCAGACCATGATGCAAACTATCCCTGCTGACCCTGGTCAGCTTGATTCCCGTGAGCAATCTGTTCCTCGTGCTGGACAAGTGTGGTTTCCAGATTCTGCAACCAAGACTGCGCAGGCATTGCTGGACTTCAACCGTGAAGGATTACCTCTGTTCATCCTCGCTAACTGGAGAGGCTTCTCTGGTGGACAAAGAGATCTTTTTGAAGGAATTCTTCAGGCTGGCTCGACTATTGTTGAGAACCTTAGGACATACAATCAGCCTGCCTTTGTCTACATTCCCATGGCTGCAGAGCTACGAGGAGGGGCTTGGGTTGTGGTTGATAGCAAGATAAACCCAGACCGCATTGAGTGCTATGCTGAGAGGAC";
my $strand     = 'R1'; # R1/R2  中间20bp序列在R1端还是R2端？

# 并行线程参考示例
my $stats = "$report_dir/$prefix.stats";
open STATS, ">$stats";
print STATS "sample\tRawReads\tCleanReads\tFinalReads\n";
my $pm = Parallel::ForkManager->new(10);
foreach my $sample(@samples)
{
    $pm->start() and next;
    
    my $r1_fastq = "$fastq_dir/$sample\_R1.fastq.gz";
    my $r2_fastq = "$fastq_dir/$sample\_R2.fastq.gz";
    my $sample_dir = "$output_dir/$sample";
    make_dir($sample_dir);
    ####################
    # (1) 基于参考序列，对reads做过滤： 90%比对成功， 前10bp比对成功
    ####################
    print "[process] $sample: filter\n";
    my $r1_clean = "$sample_dir/$sample\_R1.clean.fa";
    my $r2_clean = "$sample_dir/$sample\_R2.clean.fa";
    my $filter_stats = "$sample_dir/$sample.filter.stats";
    filter($r1_fastq, $r2_fastq, $sample_dir, $sample);

    ####################
    # (2) 提取中间目标序列
    ####################
    print "[process] $sample: abstract target seq\n";
    my $input_fa = $strand eq 'R1' ? $r1_clean : $r2_clean;  # 输入文件
    my $final_stats     = "$sample_dir/$sample.final.stats";  # 有效reads数量
    my $final_fa  = "$sample_dir/$sample.final.fa";  # 每条reads识别结果
    system("perl $SOFT_MIDDLE_SEQ_SNV -f $input_fa -l $seq_left -m $seq_middle -r $seq_right -o $sample_dir -p $sample");

   
    ####################
    # (3) 结果整理
    ####################
    my $snv_file = "$sample_dir/$sample.blast.parse.fa.snv_indel.gz";
    snv_count($snv_file, "$output_dir/$sample/$sample.tongj.txt");


    ####################
    # (4) reads 统计汇总
    ####################
    my %hashStats;
    open STATS_TMP, $filter_stats;
    my $reads_raw = <STATS_TMP>; 
       ($reads_raw) = $reads_raw=~/(\d+)/;
    my $reads_clean = <STATS_TMP>;
       ($reads_clean) = $reads_clean=~/(\d+)/;
    close STATS_TMP;
    open STATS_TMP, $final_stats;
    <STATS_TMP>;
    my $reads_middle = <STATS_TMP>;
      ($reads_middle) = $reads_middle=~/(\d+)/;
    close STATS_TMP;

    print STATS "$sample\t$reads_raw\t$reads_clean\t$reads_middle\n";

    $pm->finish;    
}
$pm->wait_all_children;

close STATS;


# EXCEL生成
my $input = "$stats,";
my $info  = "stats,";
foreach my $sample(@samples)
{
    $input .= "$output_dir/$sample/$sample.tongj.txt,";
    $info  .= "$sample,"
}

system("perl $SOFT_TABLE2EXCEL -i $input -s $info -o $report_dir/$prefix.xlsx");


###################################################### -----------子函数-------------------
sub filter{
    my $r1_fastq   = shift @_;
    my $r2_fastq   = shift @_;
    my $sample_dir = shift @_;
    my $sample     = shift @_;

    my $ref_db = "$sample_dir/ref.fa";
    open REF, ">$ref_db";
    print REF ">ref\n$seq_ref\n";
    close REF;
    system "$SOFT_MAKEBLASTDB -in $ref_db -dbtype nucl";

    my $r1_fa = "$sample_dir/$sample\_R1.fa";
    my $r2_fa = "$sample_dir/$sample\_R2.fa";
    my $r1_blast = "$sample_dir/$sample\_R1.blast";
    my $r2_blast = "$sample_dir/$sample\_R2.blast";
    my $r1_clean = "$sample_dir/$sample\_R1.clean.fa";
    my $r2_clean = "$sample_dir/$sample\_R2.clean.fa";
    system("zcat $r1_fastq | $SOFT_FASTQ_TO_FASTA - -o $r1_fa -n");
    system("zcat $r2_fastq | $SOFT_FASTQ_TO_FASTA - -o $r2_fa -n");
    system "$SOFT_BLASTN -task blastn -outfmt '6 qseqid sseqid qlen length qstart ' -query $r1_fa -evalue 0.00001 -db $ref_db -out $r1_blast";
    system "$SOFT_BLASTN -task blastn -outfmt '6 qseqid sseqid qlen length qstart ' -query $r2_fa -evalue 0.00001 -db $ref_db -out $r2_blast";
    
    # 过滤
    my %hashBlast;
    read_blast(\%hashBlast, $r1_blast, 'R1');
    read_blast(\%hashBlast, $r2_blast, 'R2');

    my $all = 0;
    my $clean = 0;
    open R1FA, $r1_fa;
    open R2FA, $r2_fa;
    open R1FACLEAN, ">$r1_clean";
    open R2FACLEAN, ">$r2_clean";
    while(my $r1_line1 = <R1FA>)
    {
        my $r1_line2 = <R1FA>;
        my $r2_line1 = <R2FA>;
        my $r2_line2 = <R2FA>;
        $all++;
        my $reads_name = (split /\s+/, $r1_line1)[0];
        $reads_name =~s/^>//;
        next if(not exists $hashBlast{$reads_name}{'R1'});
        next if(not exists $hashBlast{$reads_name}{'R2'});

        print R1FACLEAN "$r1_line1$r1_line2";
        print R2FACLEAN "$r2_line1$r2_line2";
        $clean++;
    }
    close R1FA;
    close R2FA;
    close R1FACLEAN;
    close R2FACLEAN;

    open STAT, ">$sample_dir/$sample.filter.stats";
    print STAT "raw reads\t$all\n";
    print STAT "clean reads\t$clean\n";
    close STAT;

}

sub read_blast{
    my $hashBlast = shift @_;
    my $blast = shift @_;
    my $type  = shift @_;

    open BLAST, $blast;
    while(<BLAST>)
    {
        $_=~s/[\r\n]//g;
        my ($reads_name, $ref_name, $reads_len, $align_len, $qstart) = split /\t/, $_;
        
        my $align_freq = $align_len / $reads_len;
        next if($align_freq < 0.9);  # 90% 比对成功
        next if($qstart > 10); # 前端是引物
        $hashBlast->{$reads_name}{$type}++;
 
    }
    close BLAST;
}

# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}

# 读取snv并输出统计结果
sub snv_count{
    my $snv_file = shift @_;
    my $tongji   = shift @_;
    my %hash;
    my $count = 0;
    open IN, "gzip -dc $snv_file |";
    open OUT, ">$tongji";
    print OUT "SeqName\tCount\tratio\ttype\tinsertion\tdeletion\tsnp\trefseq\tcutseq\n";
    <IN>;
    while(<IN>){
        $_ =~ s/[\r\n]//g;
        my($SeqName, $type, $insertion, $deletion, $snp, $refseq, $cutseq) = split/\t/,$_;
        my $seq = $cutseq;
        $seq =~ s/-//g;
        $hash{$seq}{'count'} ++;
        map{$_ = "-" if($_ !~ /\w/)}($type, $insertion, $deletion, $snp);
        $hash{$seq}{'print'}{"$type\t$insertion\t$deletion\t$snp\t$refseq\t$cutseq"} ++;
        $count ++;
    }
    close IN;

    my @seq = sort{$hash{$b}{'count'} <=> $hash{$a}{'count'}} keys(%hash);
    my $count_seq = 0;
    foreach my $seq_tmp(@seq){
        my $seq_count = (sort{$hash{$seq_tmp}{'print'}{$b} <=> $hash{$seq_tmp}{'print'}{$a}} keys %{$hash{$seq_tmp}{'print'}})[0];
        my $seq_sort = $hash{$seq_tmp}{'count'};
        my $ratio = $seq_sort/$count;
        $count_seq ++;
        print OUT "seq$count_seq\t$seq_sort\t$ratio\t$seq_count\n";
    }
    close OUT;
    print "[FILE] $tongji Finish\n";
}
