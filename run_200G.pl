use strict;
use warnings;
use Bio::SeqIO;
use Parallel::ForkManager;
$|=1;
# 当前序列在R1端
my @samples = qw(200-1-G 200-2-G 200-3-G);
my $fastq_dir = "/home/ganb/work/research/seq_abstract/20D0828A_PCR/fastq";
my $output_dir = "/home/ganb/work/research/seq_abstract/20D0828A_PCR/output";
my $report_dir = "/home/ganb/work/research/seq_abstract/20D0828A_PCR/report";
my $prefix = "200G";
make_dir($output_dir);
make_dir($report_dir);

my $SOFT_BLASTN         = "/home/genesky/software/blast+/2.9.0/bin/blastn";
my $SOFT_MAKEBLASTDB    = "/home/genesky/software/blast+/2.9.0/bin/makeblastdb";
my $SOFT_FASTQ_TO_FASTA = "/home/genesky/software/fastx_toolkit/0.0.14/bin/fastq_to_fasta";
my $SOFT_MUSCLE         = "/home/genesky/software/muscle/3.8.31/muscle3.8.31_i86linux64";
my $SOFT_MIDDLE_SEQ     = "perl /home/ganb/work/research/seq_abstract/20D0828A_PCR/fasta_cut_middle_seq.pl";
my $SOFT_SNP_MSA        = "perl /home/ganb/work/research/seq_abstract/20D0828A_PCR/call_snp_indel_base_on_msa_fasta.pl";
my $SOFT_TABLE2EXCEL    = "perl /home/pub/bin/NGS/chip/GATK4/tools/personal/table2excel.pl";

my $seq_left   = "TGGTCAGCTTGATTCCCGTGAGCAATCTG";
my $seq_middle = "TTCCTCGTGCTGGACAAGTGTGG";
my $seq_right  = "TTTCCAGATTCTGCAACCAAGACTGCGCAGGCATTGCTGGACTTCAACCGTGAAGGATTACCTCTGTTCATCCTCGCTAACTGGAGAGGCTTCTCTGGT";
my $seq_ref    = "TATCCGTGGTGTTGATGACAGCCAAGGGAAATGGTTAGGTGGTATGTTTGATAAAGACAGCTTTGTGGAAACATTTGAAGGTTGGGCTAAGACAGTGGTTACTGGCAGAGCAAAGCTTGGTGGAATTCCAGTGGGTGTGATAGCTGTGGAGACTCAGACCATGATGCAAACTATCCCTGCTGACCCTGGTCAGCTTGATTCCCGTGAGCAATCTGTTCCTCGTGCTGGACAAGTGTGGTTTCCAGATTCTGCAACCAAGACTGCGCAGGCATTGCTGGACTTCAACCGTGAAGGATTACCTCTGTTCATCCTCGCTAACTGGAGAGGCTTCTCTGGTGGACAAAGAGATCTTTTTGAAGGAATTCTTCAGGCTGGCTCGACTATTGTTGAGAACCTTAGGACATACAATCAGCCTGCCTTTGTCTACATTCCCATGGCTGCAGAGCTACGAGGAGGGGCTTGGGTTGTGGTTGATAGCAAGATAAACCCAGACCGCATTGAGTGCTATGCTGAGAGGACT";
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
    system("$SOFT_MIDDLE_SEQ -f $input_fa -l $seq_left -m $seq_middle -r $seq_right -o $sample_dir -p $sample");

    ####################
    # (3) 基于中间目标序列，做多重序列比对
    ####################
    print "[process] $sample: multiple sequence alignment\n";
    # MSA 输入文件制作
    my $msa_input = "$sample_dir/$sample.msa.input.fa";
    my %hashFasta = read_fasta($final_fa, 'seq_count');
    open MSA, ">$msa_input";
    print MSA ">reference\n$seq_middle\n";  # 参考序列
    my $seq_id = 0;
    foreach my $seq(sort {$hashFasta{$b}<=>$hashFasta{$a}} keys %hashFasta)
    {   
        $seq_id++;
        my $reads_count = $hashFasta{$seq};
        $seq= 'N' x length($seq_middle) if($seq eq '');  # 整个序列缺失，用N填充
        my $seq_name = "seq$seq_id\_$reads_count";  # 序列名称 = id_数量
        print MSA ">$seq_name\n$seq\n";
    }
    close MSA;
    
    # MSA 比对
    my $msa_output       = "$sample_dir/$sample.msa.output.fa";
    my $msa_output_order = "$sample_dir/$sample.msa.output.order.fa";
    system("$SOFT_MUSCLE -maxiters 16 -in $msa_input -out $msa_output ");
    order_seq($msa_input, $msa_output, $msa_output_order);  # 排序，第一个是reference

    ####################
    # (4) 基于多重序列比对结果，识别突变：第一条序列是reference
    ####################
    # MSA 突变识别
    my $msa_snp_indel = "$sample_dir/$sample.msa.snp_indel.txt";
    system("$SOFT_SNP_MSA -i $msa_output_order -o $msa_snp_indel");

    ####################
    # (5) 结果整理
    ####################
    my %hashMatrix = read_matrix($msa_snp_indel);
    open OUTPUT, ">$sample_dir/$sample.final.snp_indel.txt";
    my @heads = qw(SeqName Count ratio type insertion deletion snp seq);
    print OUTPUT (join "\t", @heads) . "\n";
    my $count_sum = $hashMatrix{'CountSum'};
    foreach my $seq_id(sort  {$hashMatrix{'DATA'}{$a}{'sort_value'}<=>$hashMatrix{'DATA'}{$b}{'sort_value'}} keys %{$hashMatrix{'DATA'}})
    {
        my @values;
        foreach my $head(@heads)
        {
            my $value = exists $hashMatrix{'DATA'}{$seq_id}{$head} ? $hashMatrix{'DATA'}{$seq_id}{$head} : "";
               $value = $seq_id if($head eq 'SeqName');
               $value = $hashMatrix{'DATA'}{$seq_id}{'Count'} / $count_sum if($head eq 'ratio');
            push @values, $value;
        }
        print OUTPUT (join "\t", @values) . "\n";
    }
    close OUTPUT;

    ####################
    # (6) reads 统计汇总
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
    $input .= "$output_dir/$sample/$sample.final.snp_indel.txt,";
    $info  .= "$sample,"
}

system("$SOFT_TABLE2EXCEL -i $input -s $info -o $report_dir/$prefix.xlsx");


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

 

sub order_seq{
    my $msa_input        = shift @_;
    my $msa_output       = shift @_;
    my $msa_output_order = shift @_;
    my %hashFastaIN = read_fasta($msa_input, 'seq');
    my %hashFastaOT = read_fasta($msa_output, 'seq');

    open ORDER, ">$msa_output_order";
    foreach my $seq_id(sort {$hashFastaIN{$a}{'sort_order'}<=>$hashFastaIN{$b}{'sort_order'}} keys %hashFastaIN)
    {
        print ORDER ">$seq_id\n$hashFastaOT{$seq_id}{'seq'}\n";
    }
    close ORDER;
}

sub read_fasta{
    my $fasta = shift @_;
    my $type  = shift @_;
    my %hashFasta;
    my $IN = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
    my $count = 0;
    while(my $inSeq = $IN->next_seq)
    {   
        $count++;
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;
        if($type eq 'seq_count')
        {
            $hashFasta{$seq}++;
        }elsif($type eq 'seq'){
            $hashFasta{$id}{'seq'}        = $seq;
            $hashFasta{$id}{'sort_order'} = $count;
        }
        
    }
    return %hashFasta;
}

# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}


# 读取带有表头的矩阵
sub read_matrix{
    my $file     = shift @_;

    print "Read $file\n";
    my %hashMatrix;
    open FILE, $file;
    my $line1 = <FILE>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;
    $hashMatrix{'HEAD'} = \@heads;

    while(<FILE>)
    {
        $_ =~ s/[\r\n]//g;
        next if($_!~/\w/);
        my @datas = split /\t/, $_;
        my $seq_name = $datas[0];
        my ($seq_id, $reads_count) = $seq_name=~/(seq\d+)_(\d+)/;
        foreach my $col(0..$#heads)
        {
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $hashMatrix{'DATA'}{$seq_id}{$heads[$col]} = $value;
        }
        $hashMatrix{'DATA'}{$seq_id}{'sort_value'} = $.;
        $hashMatrix{'DATA'}{$seq_id}{'Count'} = $reads_count;
        $hashMatrix{'CountSum'} += $reads_count;
    }
    close FILE;
    
    return %hashMatrix;
}


# # 读取blast结果
# sub read_blast{
#     my $blast           = shift @_;
#     my $type            = shift @_;###single表示对merge的reads进行比对。pair表示对没有merge上的R1或R2的reads进行比对。
    
#     # 读取blast比对信息，并进行基础过滤
#     my %hashBlast;
#     open BLAST,$blast;
#     while(<BLAST>){
#         $_=~ s/[\r\n]//g;
#         my ($reads_name, $target_name, $identify_perc, $identify_len, $mismatch, $gap, $reads_start, $reads_end, $target_start, $target_end, $evalue, $score) = split /\t/, $_;
#         my $reads_length  = length( $hashFastq->{$reads_name}{2} );
#         my $target_length = length( $hashTargetFasta->{$target_name} );
#         ($target_start, $target_end) = ($target_end, $target_start) if($target_end < $target_start);
 
#         my $judge_value = $identify_len; # 默认用匹配长度作为后续候选的筛选规则
#            $judge_value = $score if($target_name=~/^Amelo_[XY]\|[XY]/); # 如果是PCR的XY染色体判断片段的话，则采用score值作为筛选准则，因为这个片段是高度同源的，有6bp的Gap
#         #########
#         # 情况一：Reads覆盖整个目标区域
#         # 条件，覆盖超过95%的区域
#         # 条件，覆盖起始位置小于开头加10，结束位置大于末尾减10
#         #########
#         if($type eq 'Merge' and $identify_len > $target_length * 0.95 and $target_start < 10 and $target_end > $target_length - 10){
#             $hashBlast{$reads_name}{$judge_value} = "$reads_name\t$target_name\t$reads_start\t$reads_end\t$target_start\t$target_end\t$identify_len\n";
#         }
#         #########
#         # 情况二：Reads在目标区域的一端
#         # 条件，Reads长度减去一般引物长度（30bp）后，覆盖超过95%
#         # 条件，覆盖起始位置小于开头加10，或结束位置大于末尾减10
#         #########
#         if($type eq 'NoMerge' and $identify_len > ($reads_length - 30) * 0.95 and ($target_start < 10 or $target_end > $target_length - 10) ){
#             $hashBlast{$reads_name}{$judge_value} = "$reads_name\t$target_name\t$reads_start\t$reads_end\t$target_start\t$target_end\t$identify_len\n";
#         }
#     }
#     close BLAST;

#     # 进一步筛选
#     my %hashBest=();
#     if($type eq 'NoMerge')# Merge失败特殊处理
#     {
#         foreach my $title(keys %hashBlast){
#             foreach my $len(sort {$b <=> $a} keys %{$hashBlast{$title}}){
#                 my ($reads_name, $target_name, $tmp) = split /\t/, $hashBlast{$title}{$len}, 3;
#                 $hashBest{$title}{$target_name} = $hashBlast{$title}{$len} if(!exists($hashBest{$title}{$target_name}));# 保留reads所有比对结果，用于应对两个PCR片段重叠超过150的错配情况。用一端校正另一端
#             }
#         }    
#     }
#     else # 合并成功的，直接保留最佳匹配
#     {
#         foreach my $title(keys %hashBlast){
#             foreach my $len (sort {$b <=> $a} keys %{$hashBlast{$title}}){
#                 $hashBest{$title} = $hashBlast{$title}{$len};
#                 last;
#             }
#         }       
#     }

#     return %hashBest;
# }