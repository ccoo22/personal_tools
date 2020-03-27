use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
$|=1;

my $fasta_formatter = "/home/panrf/Softwares/fastx_toolkit/fasta_formatter";
my $miranda         = "/home/genesky/software/miranda/3.3/bin/miranda";
my $rnahybrid       = "/home/genesky/software/rnahybrid/2.1.2/bin/RNAhybrid";

my ($mirna, $utr3, $output_dir, $help);
GetOptions(
	"mirna|m=s"       => \$mirna,
	"utr3|u=s"        => \$utr3,
	"output_dir|o=s"  => \$output_dir,
	"help|h!"	      => \$help,
);
die help() if (defined $help or (not defined $mirna or not defined $utr3 or not defined $output_dir));
 
############################################ 主流程
mkdir $output_dir if(not -e $output_dir);

###########
# （1） miranda 预测
###########
print "Process : miranda miRNA target predict \n";
my $miranda_out     = "$output_dir/miranda_out.xls";
my $miranda_out_fmt = "$output_dir/miranda_out.fmt.xls";
system("$miranda $mirna $utr3 -sc 150 -en -15 -out $miranda_out");
fmt_miranda_out($miranda_out, $miranda_out_fmt);  # 对原始结果整理一下

###########
# （2） rnahybrid 预测
###########
print "Process : rnahybrid miRNA target predict\n";
my $max_utr3_length = call_max_length($utr3);  # 计算最长的utr3序列
my $rnahybrid_out     = "$output_dir/rnahybrid_out.xls";
my $rnahybrid_out_fmt = "$output_dir/rnahybrid_out.fmt.xls";
system("$rnahybrid -f 2,8 -c -v 3 -u 3 -e -15 -p 0.05 -m $max_utr3_length -q $mirna -t $utr3 -s 3utr_human > $rnahybrid_out");
fmt_rnahybrid_out($rnahybrid_out, $rnahybrid_out_fmt);  # 对原始结果整理一下

###########
# （3） 两个预测取交集、合并
###########
my $overlap_out = "$output_dir/rnahybrid_miranda.xls";
miranda_rnahybrid_overlap($miranda_out_fmt, $rnahybrid_out_fmt, $overlap_out);

print "Finished : $overlap_out\n";

##############子程序##############

sub miranda_rnahybrid_overlap{
	my $miranda_out_fmt   = shift @_; 
	my $rnahybrid_out_fmt = shift @_; 
	my $overlap_out       = shift @_; 
	
	my %meta;
	open RHB, $rnahybrid_out_fmt or die "Can't open $rnahybrid_out_fmt!\n";
	<RHB>;
	while(<RHB>)
	{
		$_ =~ s/[\r\n]//g;
		my ($mirna, $target, $mirna_length, $target_length, $energy, $pvalue) = split /\t/;
		$meta{"$mirna|$target"} = "$energy\t$pvalue";
	}
	close RHB;

	open OUT, ">$overlap_out" or die "Can't open $overlap_out!\n";
	print OUT "miRNA\tTranscript\tsymbol\tMiranda_Tot_score\tMiranda_Tot_energy\trnahybrid_energy\trnahybrid_Pvalue\n";
	open MRD, $miranda_out_fmt or die "Can't open $miranda_out_fmt!\n";
	<MRD>;
	while(<MRD>)
	{
		$_ =~ s/[\r\n]//g;
		my ($mirna, $target, $tot_score, $tot_energy, $max_score, $max_energy, $strand, $mirna_length, $target_length, $pos) = split /\t/, $_;
		my ($transcript, $gene) = split /[\|]/, $target;
		$gene = "" if(not defined $gene);

		next if (-s $rnahybrid_out_fmt > 0 and  not exists $meta{"$mirna\t$target"});  # 分析了rnahybrid，但是没有预测出靶向关系

		my $rnahybrid_value = $meta{"$mirna\t$target"};
		print OUT "$mirna\t$transcript\t$gene\t$tot_score\t$tot_energy\t$rnahybrid_value\n";
	}
	close MRD;
	close OUT;
}

sub fmt_rnahybrid_out{
	my $rnahybrid_out     = shift @_;
	my $rnahybrid_out_fmt = shift @_;
	
	open OUT, ">$rnahybrid_out_fmt" or die "Can't open $rnahybrid_out_fmt\n";
	open IN, $rnahybrid_out or die "Can't open $rnahybrid_out!\n";
	print OUT "miRNA\ttarget\tmirna_length\ttarget_length\trnahybrid_energy\trnahybrid_Pvalue\n";
	while(<IN>){
		$_ =~ s/[\r\n]//g;
		next if($_=~/target too long/) ;
		my @arr =  split /\s+/, $_, 2;
		my @res =  split /:/, $arr[0];
		print OUT "$res[2]\t$res[0]\t$res[3]\t$res[1]\t$res[4]\t$res[5]\n";
	}
	close IN;
	close OUT;
}


sub call_max_length{
	my $fasta = shift @_;
	my $IN = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
	my $max_length = 0;
	while(my $inSeq = $IN->next_seq)
	{
		my $length = length($inSeq->seq);
    	$max_length = $length if($length > $max_length);
	}
	return $max_length;
}

sub fmt_miranda_out{
	my $miranda_out = shift @_;
	my $miranda_out_fmt = shift @_;
	
	open OUT, ">$miranda_out_fmt" or die "Can't open $miranda_out_fmt\n";
	print OUT "miRNA\ttarget\tTot Score\tTot Energy\tMax Score\tMax Energy\tStrand\tmirna_length\ttarget_length\tPositions\n";
	my %unique;
	my $flag = 0;
	open IN, $miranda_out or die "Can't open $miranda_out\n";
	while(<IN>)
	{
		$_ =~ s/[\r\n]//g;
		$flag = 1 if($_ =~ /Score for this Scan:/);
		$flag = 0 if($_ =~ /^Complete/);
		next if $flag == 0;
		if ($_=~/>>/) 
		{
			$_=~s/^>>//;
			my @arr = split /\s+/, $_, 10;
			next if(exists $unique{$arr[0]}{$arr[1]});
			$arr[$#arr] =~ s/\s+/\|/g;
			$unique{$arr[0]}{$arr[1]} = "-";
			my $res = join "\t", @arr;
			print OUT "$res\n";
		}
	}
	close IN;
	close OUT;
}

sub help{
    my $info = "
Program: 基于miRNA序列和UTR3序列，预测靶向关系
Version: 2020-03-27
Contact: 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        --mirna/-m         成熟的miRNA序列fasta文件。,eg: /home/genesky/database/self_build_database/mirna/database/human/hsa_mature_dna.fa
        --utr3/-u          转录本UTR3序列 fasta文件。序列命名格式建议：转录本名称|基因名称
        --output_dir/-o    输出路径
        --help/-h          查看帮助文档
		
        注：UTR3序列的获取方式参考 get_3utr_from_gtf.pl 脚本
    \n";
    return $info;  
}

