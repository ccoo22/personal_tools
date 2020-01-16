# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $samtools      = "/home/pub/software/samtools/samtools-1.6/samtools";
my $tmp_dir    = "/home/tmp";


# 检测 -> 脚本输入
my ($bam, $output, $region, $if_help);
GetOptions(
    "bam|i=s"           => \$bam,
    "output|o=s"        => \$output,
    "region|r=s"        => \$region,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $bam or not defined $output ));
$region     = "" if(not defined $region);
###################################################################### 主程序

my $reads_count = 0;
open OUTPUT, ">$output";
print OUTPUT "reads name\tchr\tstart\tcigar\tMD tag\treads error\tgenome error\treads seq\n";
open BAM, "$samtools view $bam $region |";
while(<BAM>)
{   
    $reads_count++;
    my ($reads_name, $flag, $chr, $start, $quality, $cigar, $chr_mate, $start_mate, $distance, $reads_seq, $reads_qual, @tags) = split /\t/, $_;
    
    # (1) CIGAR 提取, CIGAR反应了序列的比对情况，包括缺失D,插入I
    my %hashCigar;
    my $cigar_count       = 0; # cigar数量
    my $match_position    = 0; # reads上的位置
    my $start_clip        = 0; # reads开头被softclip的长度
    my $insert_accumulate = 0; # 累计插入碱基长度
    my $genome_position   = $start - 1; # 比对基因组位置
    my @inserts; # 记录插入碱基的位置与长度（位置已校正）
    while($cigar =~ /(\d+)([A-Z])/g)
    {   
        my ($length, $type) = ($1, $2);
        $cigar_count++;
        $hashCigar{$cigar_count}{'Length'}   = $length;
        $hashCigar{$cigar_count}{'Type'}     = $type;
        $hashCigar{$cigar_count}{'Original'} = "$1$2";
        
        # 起始S处理
        if($cigar_count == 1 and $type eq 'S')
        {
            $start_clip = $length;
            $match_position += $start_clip;
            next;
        }
        if($type eq 'D')
        {
            $genome_position += $length;
            next;
        }
        next if($type eq 'H' or $type eq 'S');
        # 只需要考虑M与I
        if($type eq 'M')
        {
            $match_position  += $length;
            $genome_position += $length;
            next;
        }
        if($type eq 'I')
        {   
            my $insert_pos = $match_position + $insert_accumulate; # 当前插入的实际位置，需要考虑前面有多少插入
            push @inserts, "$insert_pos\t$length\t$chr:$genome_position"; # 记录插入的位置与长度,包括：在reads上的位置、基因组上的位置

            $insert_accumulate += $length;
        }
    }

    # (2) MD 提取, MD信息反应了Match序列的 错配与delete
    my ($MD_info) = grep{ $_ =~ /^MD/} @tags;
        $MD_info  = ":::" if(not defined $MD_info);
    my ($tag_name, $tag_vtype, $tag_value) = split /:/, $MD_info;
    my @MD_values;
    while($tag_value =~ /(\d+|[ATCGN^]+)/g)
    {
        push @MD_values, $1;
    }
    
    # (3) 分析reads错误情况，# MD中，长度只反应了CIGAR中M的部分
    my @errors;
    $genome_position = $start - 1; # 比对基因组位置
    $match_position  = $start_clip; # reads中匹配位置
    foreach my $MD_value(@MD_values)
    {   
        # match
        if($MD_value =~ /\d/)
        {   
            $match_position  += $MD_value;
            $genome_position += $MD_value;
            next;
        }
        # delete
        if($MD_value =~ /^\^/)
        {   
            $MD_value =~ s/^\^//;
            push @errors, "$match_position:D:$MD_value\t$chr:$genome_position"; # 与参考相比，丢掉的碱基
            $genome_position += length($MD_value);
            next;
        }

        # mismatch
        foreach my $ref_base(split //, $MD_value)
        {
            $match_position++;
            $genome_position++;
            push @errors, "$match_position:M\t$chr:$genome_position"; # 与参考相比，错配的碱基，错误碱基需要后面再提取
        }
    }

    # （4）位置校正（主要是插入碱基的校正）
    foreach my $insert_info(@inserts)
    {
        my ($insert_pos, $insert_length) = split /\t/, $insert_info;

        # 遍历所有错配信息，条件合适时进行校正
        foreach my $col(0..$#errors)
        {
            my ($match_position, $tmp) = split /:/, $errors[$col], 2;
            # 错配在当前插入位置之后，则需要校正
            if($match_position > $insert_pos) 
            {
                $match_position += $insert_length;
                $errors[$col] = "$match_position:$tmp";
            }
        }
    }    

    # (5) 提取错配、插入碱基
    my @finals_reads;
    my @finals_genome;
    foreach my $error_info(@errors)
    {   
        my ($error_info_reads, $genome_info) = split /\t/, $error_info;
        my ($match_position, $type, $value) = split /:/, $error_info_reads;
        if($type eq 'D')
        {
            push @finals_reads, $error_info_reads;
            push @finals_genome, "$genome_info:$type:$value";
            next;
        }
        if($type eq 'M')
        {
            my $mis_base = substr($reads_seq, $match_position - 1, 1);
            push @finals_reads, "$match_position:M:$mis_base";
            push @finals_genome, "$genome_info:M:$mis_base";
        }
    }
    foreach my $insert_info(@inserts)
    {
        my ($insert_pos, $insert_length, $genome_info) = split /\t/, $insert_info;
        my $insert_base = substr($reads_seq, $insert_pos - 1, $insert_length);
        push @finals_reads, "$insert_pos:I:$insert_base";
        push @finals_genome, "$genome_info:I:$insert_base";
    }
    
    # (6) 开头结尾的clip
    my $final_reads = join ",", @finals_reads;
       # 起始、终止的clip
       $final_reads = "Start:$hashCigar{'1'}{'Original'},$final_reads"          if(exists $hashCigar{'1'} and ($hashCigar{'1'}{'Type'} eq 'S' or $hashCigar{'1'}{'Type'} eq 'H'));
       $final_reads = "$final_reads,End:$hashCigar{$cigar_count}{'Original'}"   if(exists $hashCigar{$cigar_count} and ($hashCigar{$cigar_count}{'Type'} eq 'S' or $hashCigar{$cigar_count}{'Type'} eq 'H'));
       $final_reads = 'CompleteMatch' if($final_reads eq '' and $cigar =~ /\w/);
       $final_reads = "mapping error" if($cigar eq '*');
       $final_reads =~ s/^,|,$//g;

    my $final_genome = join ",", @finals_genome;

    # (7) 输出
    print OUTPUT "$reads_name\t$chr\t$start\t$cigar\t$tag_value\t$final_reads\t$final_genome\t$reads_seq\n";
    
    # 注意：输出结果中final中，M = mismatch; D = Deletion; I = Insertion。且D/I标注的位置为变异的前一个位置
    # die if($cigar =~ /I/);
}
close BAM;
close OUTPUT;

###################################################################### 子函数


sub help{
    my $info = "
Program: bam_reads_mark， 对每一个reads生成异常标志
Version: 2019-02-21
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --bam/-i         [必填] bam文件, 示例：/home/pub/output/Target/18B1208B/NC1003/NC1003_final.bam
         --output/-o      [必填] 结果文件, 示例：/home/ganb/work/tmp/bam_reads_mis/reads_mis.txt
         --region/-r      区域, 多个区域用‘,’分割, 示例：6:18139134-18139334,13:48611612-48611812
         --help/-h        查看帮助文档
    \n";
    return $info;
}

