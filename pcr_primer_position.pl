$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
基于引物序列查找参考序列
Version: v1.0 2020-06-29
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

my $DEFAULT_SOFT_MEGABLAST   = "/home/genesky/software/blast/20110130/bin/megablast";
my $DEFAULT_BLAST_WORD       = 12;
my $DEFAULT_RM_UNNORMAL_CHR  = 'yes';
my $DEFAULT_PRODUCT_SIZE_MIN = 150;
my $DEFAULT_PRODUCT_SIZE_MAX = 400;

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($primer, $genome, $output_dir, $SOFT_MEGABLAST, $BLAST_WORD, $RM_UNNORMAL_CHR, $PRODUCT_SIZE_MIN, $PRODUCT_SIZE_MAX, $if_help);
GetOptions(
	"primer|p=s"        => \$primer,
	"genome|g=s"        => \$genome,
	"output_dir|o=s"    => \$output_dir,

	"megablast=s"       => \$SOFT_MEGABLAST,
	"blast_word=s"      => \$BLAST_WORD,
	"rm_unnormal_chr=s" => \$RM_UNNORMAL_CHR,
	"product_size_min=s" => \$PRODUCT_SIZE_MIN,
	"product_size_max=s" => \$PRODUCT_SIZE_MAX,
	"help|h"            => \$if_help,
);
die "
Options: 必填

        --primer/-p                引物文件。一行对应一条引物，两列构成，第一列引物名称，第二列引物序列。
                                   ACAT1_14F GGCCTGCTGTAATCAGTGTGA
                                   ACAT1_14R AGACTCAGAATGCAAAAATGTATCAAA
                                   注：格式采用 ‘基因名’+‘_’+‘任意字符’+‘F/R’。
                                   上述示例中，基因名是PVT1，分割符：’_‘
        --genome/-g                参考基因组fasta文件.
                                   注意：需要有blast+索引文件， 例如： /home/genesky/database/ucsc/hg19/genome/blast+_idx/hg19.fa
        --output_dir/-o            结果输出路径

Options: 可选
        --megablast                更改软件 megablast 版本 (default: '$DEFAULT_SOFT_MEGABLAST')
        --blast_word               修改blast word 参数 (default: '$DEFAULT_BLAST_WORD')
        --rm_unnormal_chr          剔除类似 chr6_cox_hap2的结果 yes/no (default: '$DEFAULT_RM_UNNORMAL_CHR')
        --product_size_min         产物最小长度 (default: '$DEFAULT_PRODUCT_SIZE_MIN')
        --product_size_max         产物最大长度 (default: '$DEFAULT_PRODUCT_SIZE_MAX')
        --help/-h                  查看帮助文档

\n" if (defined $if_help or not defined $primer or not defined $genome or not defined $output_dir);

$SOFT_MEGABLAST   = $DEFAULT_SOFT_MEGABLAST if (not defined $SOFT_MEGABLAST); 
$BLAST_WORD       = $DEFAULT_BLAST_WORD if (not defined $BLAST_WORD); 
$RM_UNNORMAL_CHR  = $DEFAULT_RM_UNNORMAL_CHR if (not defined $RM_UNNORMAL_CHR); 
$PRODUCT_SIZE_MIN  = $DEFAULT_PRODUCT_SIZE_MIN if (not defined $PRODUCT_SIZE_MIN); 
$PRODUCT_SIZE_MAX  = $DEFAULT_PRODUCT_SIZE_MAX if (not defined $PRODUCT_SIZE_MAX); 

$output_dir = File::Spec->rel2abs($output_dir);  # 转换为绝对路径
mkdir $output_dir if(not -e $output_dir);

###################################################################### 初始化


###################################################################### 主程序

# （1）引物读入
print "[process] read primer\n";
my %hashPrimer=readPrimer($primer);

# （2）引物序列输出，准备比对
print "[process] blast to genome\n";
my $primerSeqFile  = "$output_dir/primerSeq.fa";# 引物Fasta序列
my $primerSeqBlast = "$output_dir/primerSeq.fa.blast";# 引物Fasta序列blast结果
outputPrimer(\%hashPrimer,$primerSeqFile); 
system("$SOFT_MEGABLAST -a 20 -d $genome -i $primerSeqFile -p 0.99 -F F -o $primerSeqBlast -D 3 -W $BLAST_WORD");
my %hashPrimerPos=readBlast($primerSeqBlast,\%hashPrimer);# 读取比对位置

# (3) 确定引物的产物序列
print "[process] defined product\n";
my %hashLost;# 引物没有比对上记录
my %hashWarnings;# 警告记录
my %hashProduct=primerProduct(\%hashPrimer,\%hashPrimerPos);

# (3) 结果输出
print "[process] output\n";
my $productFile = "$output_dir/product.txt";
my $warningFile = "$output_dir/warnings.txt";
OutputProduct(\%hashProduct,$productFile);# 产物输出
showOthers($warningFile);

print "\n";
print "[result] product file $productFile \n";
print "[result] warning file $warningFile \n";
###################################################################### 子程序


sub OutputProduct{
	my $hashProduct = shift @_;
	my $productFile = shift @_;
	open PRODUCT,">$productFile";
	print PRODUCT "title\tprimer_F\tprimer_R\tchr\tpos1\tpos2\tpos3\tpos4\tlen\tExtra\n";
	foreach my $group(keys %$hashProduct){
		my @mapPairs=keys %{$hashProduct->{$group}};
		my $mark=(@mapPairs>1) ? 'HOMOLOGY' : "";
		foreach my $mapPair(@mapPairs){
			print PRODUCT "$hashProduct->{$group}{$mapPair}\t$mark\n";
		}
	}
	close PRODUCT;
}
sub showOthers{
	my $warningFile = shift @_;
    open WARNING, ">$warningFile";

	my @losts    = sort keys %hashLost;
	my @warnings = sort keys %hashWarnings;	
	if(@losts>=1){
		my $lostCount=@losts;
	    print WARNING "### Primer Not Map $lostCount ###\n";
	    foreach my $lost(@losts){
	    	print WARNING "$lost\t$hashLost{$lost}\n";
	    }
	}
    if(@warnings>=1){
	    print "### Warnings ###\n";
	    foreach my $warning(@warnings){
	    	print WARNING "$warning\t$hashWarnings{$warning}\n";
	    }
    }
    close WARNING;
}
sub primerProduct{
	my $hashPrimer    = shift @_;
	my $hashPrimerPos = shift @_; 
	my %hashProduct;
	my %hashSizeBad;
	foreach my $group(keys %$hashPrimer){
		# 引物没有匹配到基因组
		if(!exists($hashPrimerPos{$group})){
			$hashLost{$group}="No Map";
			next;
		}
		# 引物只有一端匹配到基因组
		my @primerNames=keys %{$hashPrimerPos{$group}};
		if(@primerNames!=2){
			my @lostInfos;
			foreach my $primerName(@primerNames){
				my @positions=sort keys %{$hashPrimerPos->{$group}{$primerName}};
				my $positionInfo=join ",",@positions;
		        push @lostInfos,"$primerName = $positionInfo";
			}
			$hashLost{$group}="Single Map\t" .(join "\t",@lostInfos);
			next;
		}
		# 结果配对
		my ($primerFName,$primerRName)=getPrimerPair($hashPrimer->{$group});# 获取F/R端引物名称
		my $primerFSeq              = $hashPrimer->{$group}{$primerFName}{'Seq'};# 前端引物
		my $primerRSeq              = $hashPrimer->{$group}{$primerRName}{'Seq'};# 后端引物
		my @primerFPositions        = keys %{$hashPrimerPos->{$group}{$primerFName}};# 前端引物匹配的位置
		my @primerRPositions        = keys %{$hashPrimerPos->{$group}{$primerRName}};# 后端引物匹配的位置
		my $primerFPositionString   = join ",",@primerFPositions;
		my $primerRPositionString   = join ",",@primerRPositions;
		$hashWarnings{$primerFName} = "HOMOLOGY:\t$primerFPositionString" if(@primerFPositions>1);#同源警告
		$hashWarnings{$primerRName} = "HOMOLOGY:\t$primerRPositionString" if(@primerRPositions>1);

		foreach my $FPosition(@primerFPositions){
			my ($Fchr,$FStart,$FEnd) = split /\|/,$FPosition;
			foreach my $RPosition(@primerRPositions){
				my ($Rchr,$RStart,$REnd) = split /\|/,$RPosition;
				next if($Fchr ne $Rchr);
				my $ProductLength = ($FStart<$RStart) ? $REnd-$FStart+1 : $FEnd-$RStart+1;
				if($ProductLength<$PRODUCT_SIZE_MIN or $ProductLength>$PRODUCT_SIZE_MAX){# 不满足片段大小要求
					my $info = "Length=$ProductLength\t$primerFName=$FPosition,$primerRName=$RPosition";
					if(exists($hashSizeBad{$group})){
						$hashSizeBad{$group}.="\n$info";
					}else{
						$hashSizeBad{$group}=$info;
					}
					next;
				}
				my @positions      = ($FStart,$FEnd,$RStart,$REnd);
				   @positions      = sort {$a<=>$b} @positions;
				my $positionString = join "\t",@positions;
				$hashProduct{$group}{"$FPosition|$RPosition"}="$group\t$primerFSeq\t$primerRSeq\t$Fchr\t$positionString\t$ProductLength";
			}
		}
	}
	# 检查有比对信息，但无法构成合适成对引物的结果，并放入lost与warnings
	foreach my $group(keys %hashSizeBad){
		if(exists($hashProduct{$group})){# 某对引物找到了合适的产物，但是因为同源的关系，存在产物长度不合格的信息，这部分放入警告信息中 
			$hashWarnings{$group}="Size Bad\t".$hashSizeBad{$group};
		}else{
			$hashLost{$group}="Size Bad\t".$hashSizeBad{$group};
		}
	}
	# 检查是否有染色体不同导致的无序列匹配
	foreach my $group(keys %$hashPrimer){
		next if(exists($hashLost{$group}));# 已记录为丢失
		next if(exists($hashProduct{$group}));# 已找到序列
		my ($primerFName,$primerRName)=getPrimerPair($hashPrimer->{$group});# 获取F/R端引物名称
		my @primerFPositions        = keys %{$hashPrimerPos->{$group}{$primerFName}};# 前端引物匹配的位置
		my @primerRPositions        = keys %{$hashPrimerPos->{$group}{$primerRName}};# 后端引物匹配的位置
		my $primerFPositionString   = join ",",@primerFPositions;
		my $primerRPositionString   = join ",",@primerRPositions;
		$hashLost{$group}="Map Chr Error\t$primerFName\t$primerFPositionString\t$primerRName\t$primerRPositionString";
	}
	return %hashProduct;
}

sub readBlast{
	my $primerSeqBlast = shift @_;
	my $hashPrimer     = shift @_;
	my %hashPrimerPos;
	open BLAST,$primerSeqBlast or die "No Blast Result!\n";
	while(<BLAST>){
		$_=~s/[\r\n]//g;
		next if($_=~/^#/);
		next if($_!~/\w/);
		my ($Query_id,$Subject_id,$identity,$alignment_length,$mismatches,$gap_openings,$q_start,$q_end,$s_start,$s_end,$e_value,$bit_score)=split /\t/,$_;
		my ($group,$primerName)=split /\|/,$Query_id;
		my $chr=$Subject_id;
		$chr=~s/^lcl\|//;
		next if($RM_UNNORMAL_CHR eq 'yes' and $chr =~ /_/);
		my $alignP=($alignment_length)/length($hashPrimer->{$group}{$primerName}{'Seq'});
		if($alignP==1){# 完全匹配
			my $primerStart = ($s_start<$s_end) ? $s_start : $s_end;
			my $primerEnd   = ($s_start<$s_end) ? $s_end   : $s_start;
			$hashPrimerPos{$group}{$primerName}{"$chr|$primerStart|$primerEnd"}=1;
		}
	}
	close BLAST;
	return %hashPrimerPos;
}

# 输出引物序列用于比对
sub outputPrimer{
	my $hashPrimer    = shift @_;
	my $primerSeqFile = shift @_;
	open PRIMER,">$primerSeqFile";
	foreach my $group(sort keys %$hashPrimer){
		foreach my $primerName(sort keys %{$hashPrimer->{$group}}){
			print PRIMER ">$group|$primerName\n";
			print PRIMER "$hashPrimer->{$group}{$primerName}{'Seq'}\n";
		}
	}
	close PRIMER;
}

sub readPrimer{
	my $file       = shift @_;
	my %hashPrimer = ();
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		$line=~ s/^\s+//;
		next if($line!~/\w/);
		my @split_line=split /[\s\t]+/,$line;
		if(@split_line>1){
			my $lastFR=0;
			while($split_line[0]=~ /[FR]/ig){
				my $pos=pos($split_line[0]);
				$lastFR=$pos;
			}
			my $group=substr($split_line[0],0,$lastFR-1);
			my $primerType=substr($split_line[0],$lastFR-1,1);# 引物方向F/R
			if(exists($hashPrimer{$group}{$split_line[0]})){
				print "Duplicate Primer Name : $split_line[0] (Program will only keep  the last one)<br>";
			}
			$hashPrimer{$group}{$split_line[0]}{'Seq'}=$split_line[1];
			$hashPrimer{$group}{$split_line[0]}{'Strand'}=$primerType;
		}
	}
	close FILE;
	my %del=();
	foreach my $group(keys %hashPrimer){
		if(keys %{$hashPrimer{$group}}!=2){
			$del{$group}=1;
		}
		my ($primerFName,$primerRName)=getPrimerPair($hashPrimer{$group});# 要求F/R成对
		if($primerFName eq ''){
			$del{$group}=1;
		}
	}
	foreach my $group(keys %del){
		print "[WARN] Delete Group Primer (not pair), $group";
		delete $hashPrimer{$group};
	}
	if((keys %hashPrimer)==0){
		die "No Primer Input!\n";
	}
	return %hashPrimer;
}
sub getPrimerPair{
	my $hashPrimerGroup = shift @_;
	my @primerNames=keys %{$hashPrimerGroup};
	my %hashTmp;
	$hashTmp{$hashPrimerGroup->{$primerNames[0]}{'Strand'}}=$primerNames[0];
	$hashTmp{$hashPrimerGroup->{$primerNames[1]}{'Strand'}}=$primerNames[1];
	return ('','') if(!exists($hashTmp{'F'}) or !exists($hashTmp{'R'}));
	my $primerFName=$hashTmp{'F'};
	my $primerRName=$hashTmp{'R'};
	return ($primerFName,$primerRName);
}
