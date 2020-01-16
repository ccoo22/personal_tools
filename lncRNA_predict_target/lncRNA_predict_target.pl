# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use Parallel::ForkManager;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $annotate_variation = "/home/genesky/software/annovar/2018Apr16/annotate_variation.pl";
my $gtf_to_fasta       = "/home/genesky/software/tophat/2.1.1/local_build/bin/gtf_to_fasta";
my $lnctar             = "/home/xudl/soft/LncTar/LncTar.pl"; # 注意，/home/tmp目录下会生成临时文件“$mRNA_base\_mRNA_temp.txt”与“$lncRNA_base\_lncRNA_temp.txt”。注意输入文件不要重复
my $table_to_excel     = SCRIPTDIR . "/tools/table2excel.pl";
my $readme             = SCRIPTDIR . "/tools/readme.txt";

# 检测 -> 脚本输入
my ($lncrna_gene, $gtf, $genome, $output_dir, $threads, $lncrna_extend, $if_help);
GetOptions(
    "lncrna_gene|l=s"   => \$lncrna_gene,
    "gtf|f=s"           => \$gtf,
    "genome|g=s"        => \$genome,
    "output_dir|o=s"    => \$output_dir,
    "threads|t=s"       => \$threads,
    "lncrna_extend|e=s" => \$lncrna_extend,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $lncrna_gene or not defined $gtf or not defined $genome or not defined $output_dir));
$threads = 10 if(not defined $threads);
$lncrna_extend = 10000 if(not defined $lncrna_extend); # 默认lncRNA上下扩展10K
###################################################################### 主程序
#
make_dir($output_dir);
my %hashGTFLncRNA = get_gtf_lncrna_biotype(); # # 根据GENCODE注释，获取LncRNA基因的biotype

# (1) 计算数据库中 需要的LncRNA上下10K范围内的mRNA
my $lncrna_10k_mrna_dir = "$output_dir/1.lncrna_10k_mrna"; # 临近10K计算目录
my $lncrna_10k_mrna     = "$lncrna_10k_mrna_dir/lncrna_10k_mrna.txt"; # 最终文件
make_dir($lncrna_10k_mrna_dir);
check_lncrna_10k_mrna($lncrna_gene, $gtf, $lncrna_10k_mrna_dir, $lncrna_10k_mrna);

# (2) gtf 转 fasta(把1中的RNA提取fasta)
my $gtf_to_fasta_dir = "$output_dir/2.gtf_to_fasta";
my $gtf_fasta        = "$gtf_to_fasta_dir/gtf.fasta";
make_dir($gtf_to_fasta_dir);
gtf_to_fasta_convert($lncrna_10k_mrna, $gtf_to_fasta_dir, $gtf_fasta);


# (3) 靶基因预测
my $predict_dir    = "$output_dir/3.predict";
my $lncrna_targets = "$output_dir/lncrna_targets.txt"; # 结果汇总文件
make_dir($predict_dir);
predict_target($predict_dir, $lncrna_targets, $lncrna_10k_mrna, $gtf_fasta);

# (4) 结果输出
my $excel = "$output_dir/lncRNA_Targets.xlsx";
print "\n\n[4] result : $excel\n";
system("perl $table_to_excel -i $lncrna_targets,$readme -s Target,ReadMe -o $excel");



###################################################################### 子函数
 

# 预测靶基因
sub predict_target{
    my $predict_dir     = shift @_;
    my $lncrna_targets  = shift @_;
    my $lncrna_10k_mrna = shift @_;
    my $gtf_fasta       = shift @_;
    print "\n\n[3]: Predict target\n";

    # (1) 读取所有fasta序列
    print "    Read fasta: $gtf_fasta\n";
    my %hashFasta;
    my $fastaIn   = Bio::SeqIO->new(-file => $gtf_fasta, -format=>'Fasta') or die "Could not open up file $gtf_fasta: $!";
    while(my $inSeq = $fastaIn->next_seq)
    {
        my $id   = $inSeq->display_id;
        my $desc = $inSeq->desc;
        my $seq  = $inSeq->seq;
        my ($tx) = split /\s+/, $desc;
        die "duplicate transcript : $tx\n" if(exists $hashFasta{$tx});
        $hashFasta{$tx} .= $seq;
    }    

    # (2) 读取所有的差异数据lncRNA -> mRNA配对情况,并生成fasta文件
    print "    output fasta :";
    my $random_word = get_random_word(10); # 项目随机长度字符，防止重复冲突
    my %hashMAP;# 记录要分析的lncRNA以及文件路径
    my $count = 0;
    open MAP, $lncrna_10k_mrna; # 10K范围内对应列表
    while(<MAP>)
    {
        $_=~s/[\r\n]//g;
        my ($lncRNA_info, $mRNA_info_list) = split /\t/, $_;
        my ($lncRNA_gene, $lncRNA_tx) = split /\|/, $lncRNA_info;
        $count++;

        # 当前lncRNA结果目录
        my $lncRNA_gene_tx = "$lncRNA_gene\_$lncRNA_tx";
        my $lncrna_dir     = "$predict_dir/$lncRNA_gene_tx";
        system("rm -r $lncrna_dir") if(-f $lncrna_dir); # 删除已存在的，如果项目重复运行的话，数据会因为$random_word字符，而积累很多
        make_dir($lncrna_dir);

        # 注意，lncTar会在/home/tmp目录下会生成临时文件“$mRNA_base\_mRNA_temp.txt”与“$lncRNA_base\_lncRNA_temp.txt”。注意不同的分析输入文件前缀不要重复，否则分析结果不对
        # lncRNA 序列
        my $lncRNA_fasta = "$lncrna_dir/$lncRNA_gene_tx.$random_word.lncRNA.fasta";
        open FASTA, ">$lncRNA_fasta";
        print FASTA ">$lncRNA_info\n";
        print FASTA "$hashFasta{$lncRNA_tx}\n";
        close FASTA;
        # mRNA 序列
        my $mRNA_fasta = "$lncrna_dir/$lncRNA_gene_tx.$random_word.mRNA.fasta";
        open FASTA, ">$mRNA_fasta";
        foreach my $mRNA_info(split /,/, $mRNA_info_list)
        {
            my ($mRNA_gene, $mRNA_tx) = split /\|/, $mRNA_info;
            print FASTA ">$mRNA_info\n";
            print FASTA "$hashFasta{$mRNA_tx}\n";
        }
        close FASTA;

        $hashMAP{$lncRNA_gene_tx}{'lncRNA_fasta'} = $lncRNA_fasta; 
        $hashMAP{$lncRNA_gene_tx}{'mRNA_fasta'}   = $mRNA_fasta; 
    }
    close MAP;

    print " total find lncRNA count: $count\n";
    if($count == 0)
    {
        print "[warnings] not find lncRNA\n";
        exit;
    }

    # (3) 靶基因预测
    print "    start predict\n";
    my @lncRNA_gene_txs = sort keys %hashMAP;
    my $pm = Parallel::ForkManager->new($threads);
    $pm->run_on_start(sub{my ($pid, $lncRNA_gene_tx) = @_; process_bar_array($lncRNA_gene_tx, \@lncRNA_gene_txs)});# 进度条
    foreach my $lncRNA_gene_tx(@lncRNA_gene_txs)
    {
        $pm->start($lncRNA_gene_tx) and next;

        my $lncrna_dir = "$predict_dir/$lncRNA_gene_tx";
        system("perl $lnctar -p 1 -l $hashMAP{$lncRNA_gene_tx}{'lncRNA_fasta'} -m $hashMAP{$lncRNA_gene_tx}{'mRNA_fasta'}  -d -0.1 -s F -o $lncrna_dir/target.txt &> $lncrna_dir/run.log\n");
        # 删除软件临时文件
        system("rm /home/tmp/$lncRNA_gene_tx.$random_word.mRNA\_mRNA_temp.txt"); 
        system("rm /home/tmp/$lncRNA_gene_tx.$random_word.lncRNA\_lncRNA_temp.txt");
        $pm->finish;
    }
    $pm->wait_all_children;

    # (4) 结果汇总
    print "    summary result ： $lncrna_targets\n";
    open SUMMARY, ">$lncrna_targets";
    my $count_result = 0;  # 读入的文件个数
    foreach my $lncRNA_gene_tx(@lncRNA_gene_txs)
    {
        my $result = "$predict_dir/$lncRNA_gene_tx/target.txt";
        next if(not -e $lncrna_targets);
        $count_result++;

        # 读入
        open RESULT, $result;
        my $line1 = <RESULT>;
           $line1 =~ s/[\r\n]//g;
        my @heads = split /\s+/, $line1; # 软件原始结果中，表头空字符很混乱
        if($count_result == 1) # 第一个文件读入，输出表头
        {
            print SUMMARY (join "\t", @heads) . "\n";
            next;
        }
        while(<RESULT>)
        {
            $_=~s/[\r\n]//g;
            my @datas = split /\t+/, $_;
            print SUMMARY (join "\t", @datas) . "\n";
        }
        close RESULT;

    }
    close SUMMARY;
}

# 获取随机长度字符
sub get_random_word{
    my $length   = shift @_;
    my @aword    = (0..9,'a'..'z','A'..'Z');
    my $password = join '', map { $aword[int rand @aword] } 0..($length-1);
    return $password;
}
sub gtf_to_fasta_convert{
    my $lncrna_10k_mrna  = shift @_;
    my $gtf_to_fasta_dir = shift @_; # 输出目录
    my $gtf_fasta        = shift @_; # 输出文件
    print "\n\n[2]: gtf_to_fasta_convert \n";

    # （1）获取需要的转录本编号
    print "    Read need transcript: $lncrna_10k_mrna\n";
    my %hashTX;
    open MAP, $lncrna_10k_mrna;
    while(<MAP>)
    {
        $_=~s/[\r\n]//g;
        $_=~s/\t/,/;
        foreach my $gene_info(split /,/, $_)
        {
            my ($gene, $tx) = split /\|/, $gene_info;
            $hashTX{$tx}++;
        }
    }
    close MAP;

    # （2）提取GTF
    my $select_gtf = "$gtf_to_fasta_dir/select.gtf";
    print "    Select gtf: $select_gtf\n";
    open GTF, $gtf;
    open SELECT, ">$select_gtf";
    while(<GTF>)
    {
        $_=~s/[\r\n]//g;
        next if($_ =~ /^#/);
        my ($transcript_id) = $_ =~ /transcript_id\s+\"(.+?)\"/;
        next if(not defined $transcript_id);
        next if(not exists $hashTX{$transcript_id}); # 按照转录本编号提取
        print SELECT "$_\n";
    }
    close GTF;    
    close SELECT;

    # （3）gtf to fasta
    print "    gtf to fasta: $gtf_fasta\n";
    system("$gtf_to_fasta $select_gtf $genome $gtf_fasta");
}

sub check_lncrna_10k_mrna{
    my $lncrna_gene         = shift @_;
    my $gtf                 = shift @_; # 输入总数据库
    my $lncrna_10k_mrna_dir = shift @_; # 结果目录
    my $lncrna_10k_mrna     = shift @_; # 结果

    print "\n\n[1]: check_lncrna_10k_mrna \n";    

    # (1) 获取需要的基因
    print "    Read need gene\n";
    my %hashNeedGene;
    open GENE, $lncrna_gene;
    while(<GENE>)
    {
        $_=~s/[\r\n\s]//g;
        next if($_!~/\w/);
        $hashNeedGene{$_}++;
    }
    close GENE;

    # （2） 获取RNA区域bed信息
    print "    Read gtf $gtf\n";
    my $bed_lncrna  = "$lncrna_10k_mrna_dir/gtf_lncRNA_bed.txt";
    my $bed_mrna    = "$lncrna_10k_mrna_dir/gtf_mRNA_bed.txt";

    open BED_LINCRNA, ">$bed_lncrna";
    open BED_MRNA, ">$bed_mrna";
    my $count_lincrna = 0;
    my $count_mrna    = 0;
    open GTF, $gtf;
    while(<GTF>)
    {
        $_=~s/[\r\n]//g;
        next if($_ =~ /^#/);
        my ($chr, $source, $region_type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
        next if($region_type ne 'transcript'); # 只需要看转录组
        next if($chr eq 'MT'); # 不要MT

        my %hashAttributes  = parse_attrbutes($attributes); # 拆分GTF注释
        my $gene_id            = $hashAttributes{'gene_id'};
        my $gene_name          = $hashAttributes{'gene_name'};
        my $gene_biotype       = $hashAttributes{'gene_biotype'}; # 注意：按照基因名，存在即是protein_coding基因、也是antisense的基因...。故按照gene_id或者transcript_id区分最佳。
        my $transcript_id      = $hashAttributes{'transcript_id'};
        my $transcript_biotype = $hashAttributes{'transcript_biotype'}; # 
        if(exists($hashGTFLncRNA{$transcript_biotype}) and exists $hashNeedGene{$gene_name}) # LncRNA,且是分析的基因
        {   
            $hashNeedGene{$gene_name}++;
            my $start_extend = $start - $lncrna_extend;
            my $end_extend   = $end + $lncrna_extend;
            print BED_LINCRNA "$chr\t$start_extend\t$end_extend\t0\t-\t$gene_name|$transcript_id|$chr|$start|$end\n";
            $count_lincrna++;
        }
        if($transcript_biotype eq 'protein_coding')
        {   
            print BED_MRNA "$count_mrna\t$chr\t$start\t$end\t$gene_name|$transcript_id|$chr|$start|$end\n";
            $count_mrna++;
        }        

    }
    close BED_LINCRNA;
    close BED_MRNA;
    close GTF;

    # (3) 注释(通过annovar的区域注释功能，注释上下10K内的mRNA)
    print "    annovar anno\n";
    my $region_anno = "$bed_lncrna.gtf_mRNA_bed";
    system("$annotate_variation -regionanno -dbtype mRNA_bed --buildver gtf $bed_lncrna $lncrna_10k_mrna_dir");

    # （4）汇总
    print "    summary : $lncrna_10k_mrna\n";
    open ANNO, $region_anno;
    open OUT, ">$lncrna_10k_mrna";
    while(<ANNO>)
    {
        $_=~s/[\r\n]//g;
        my ($title, $mrna_list, $chr, $start, $end, $ref, $alt, $lincRNA_info) = split /\t/, $_;
        my ($lincRNA_gene, $lincRNA_tx, $tmp) = split /\|/, $lincRNA_info, 3;
        
        my @near_mRNA;
        $mrna_list =~ s/^Name=//;
        foreach my $anno_value(split /,/, $mrna_list)
        {
            my ($mRNA_gene, $mRNA_tx, $tmp) = split /\|/, $anno_value;
            push @near_mRNA, "$mRNA_gene|$mRNA_tx";
        }
        print OUT "$lincRNA_gene|$lincRNA_tx\t" . (join ",", @near_mRNA) . "\n";
    }
    close ANNO;
    close OUT;

    # 检查是否有丢失的lncRNA基因
    my @lost_genes = grep{ $hashNeedGene{$_} == 1} sort keys %hashNeedGene;
    my $lost_gene_count = scalar(@lost_genes);
    if($lost_gene_count > 0)
    {   
        my $lost = join ",", @lost_genes;
        system "echo -e '\\033[41;37m[Warnings] $lost_gene_count lncRNA not find in GTF: $lost \\033[0m'"; # 错误，红底白字
    }
}

# 拆分GTF注释列
sub parse_attrbutes{
    my $attributes = shift @_;
    my %hashAttributes;
    foreach my $attribute(split /;/, $attributes)
    {
        $attribute=~ s/^\s+//;
        $attribute=~ s/\s+$//;
        my ($name, $value) = split /\s+/, $attribute;
        $value =~ s/^\"//;
        $value =~ s/\"$//;
        $hashAttributes{$name} = $value;
    }
    return %hashAttributes;
}

# 根据GENCODE注释，获取LncRNA基因的biotype
sub get_gtf_lncrna_biotype{
    my %hashGTFLncRNA;
    $hashGTFLncRNA{"3prime_overlapping_ncRNA"}++;
    $hashGTFLncRNA{"antisense"}++;
    $hashGTFLncRNA{"bidirectional_promoter_lncRNA"}++;
    $hashGTFLncRNA{"lincRNA"}++;
    $hashGTFLncRNA{"macro_lncRNA"}++;
    $hashGTFLncRNA{"non_coding"}++;
    $hashGTFLncRNA{"processed_transcript"}++;
    $hashGTFLncRNA{"sense_intronic"}++;
    $hashGTFLncRNA{"sense_overlapping"}++;
    $hashGTFLncRNA{"TEC"}++;
    return %hashGTFLncRNA;
}

sub make_dir{
    my $dir = shift @_;
    mkdir $dir if(not -e $dir);
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


sub help{
    my $info = "
Program: NCRNA 靶基因预测
Version: 2019-01-22
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --lncrna_gene/-l    【必填】要分析的lncRNA基因列表，一列。示例：/home/ganb/work/tmp/19C0515A/report/deseq2_diff_lincRNA/diff.rna.txt
         --gtf/-f            【必填】gtf文件。示例：/home/ganb/work/tmp/19C0515A/report/Homo_sapiens.GRCh38.86.chr.gtf
         --genome/-g         【必填】参考基因组，注意序列名称与gtf保持一致。示例：home/ganb/work/tmp/HG38/ENSEMBL/chromosose/Homo_sapiens.GRCh38.dna_sm.chromosome.fa
         --output_dir/-o     【必填】结果输出目录
         --threads/-t        并行线程数量。控制并行lncRNA基因的数量。默认：10
         --lncrna_extend/-e  lncRNA上下扩展的范围。默认：10000。靶基因预测lncRNA上下lncrna_extend范围内的mRNA。
         --help/-h           查看帮助文档
    \n";
    return $info;
}

 