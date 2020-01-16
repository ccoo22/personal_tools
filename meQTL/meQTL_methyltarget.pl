# 导入 -> 系统 package   
use warnings;
use strict;
use Encode;
use File::Spec;
use Getopt::Long;
use Statistics::R;
use Excel::Writer::XLSX;

my $RBIN = "/home/genesky/software/r/3.5.1/bin/R";
my $RLIB = "/home/genesky/software/r/3.5.1/lib64/R/library/";

my ($methyl_file, $snp_file, $genotype_file, $group_file, $covariates_file, $output_dir, $if_help);
GetOptions(
	"methyl_file|m=s" => \$methyl_file,
	"snp_file|s=s" => \$snp_file,
	"genotype_file|geno=s" => \$genotype_file,
	"group_file|g=s" => \$group_file,
	"covariates_file|c=s" => \$covariates_file,
	"output_dir|o=s" => \$output_dir,
	"help|h" => \$if_help,
);
die help() if(defined $if_help or not defined $methyl_file or not defined $snp_file);
$genotype_file   = ""   if(not defined $genotype_file);
$group_file      = ""   if(not defined $group_file);
$covariates_file = ""   if(not defined $covariates_file);
$output_dir      = "./" if(not defined $output_dir);
######################################################################################## 主程序
if($group_file ne ""){
	my %hashGroup = readGroup($group_file);
	foreach my $group_name(sort keys %hashGroup){
		my $sub_output_dir = "$output_dir/$group_name";
		mkdir $sub_output_dir if(not -e $sub_output_dir);		
		my @samples = split/,/,$hashGroup{$group_name};
		my %hashMethyl = extract_methyl_snv_file(\@samples, $methyl_file, "$sub_output_dir/methyl.in"); # 提取指定分组样本的数据
		my %hashSNP    = extract_methyl_snv_file(\@samples, $snp_file, "$sub_output_dir/snp.in");
		my %hashCovar  = extract_methyl_snv_file(\@samples, $covariates_file, "$sub_output_dir/covariates.in") if($covariates_file ne "");

		check_sample_names("$sub_output_dir/methyl.in", "$sub_output_dir/snp.in"); # 检查输入文件样本是否对应
		check_sample_names("$sub_output_dir/methyl.in", "$sub_output_dir/covariates.in") if($covariates_file ne "");

		eqtl($sub_output_dir,"$sub_output_dir/methyl.in","$sub_output_dir/snp.in","$sub_output_dir/covariates.in");

		my %hashQTL = txt2excel("$sub_output_dir/eqtl.out", "$sub_output_dir/meQTL.xlsx"); #输出excel并返回分析结果

		if($genotype_file ne ""){
	        my %hashGeno = extract_methyl_snv_file(\@samples, $genotype_file, "$sub_output_dir/genotype.in"); # 提取基因型数据
	        plot_boxplot(\@samples, \%hashQTL, \%hashMethyl, \%hashGeno, $sub_output_dir);
        }
	}
}else{
	check_sample_names($methyl_file, $snp_file); # 检查输入文件样本是否对应
	check_sample_names($methyl_file, $covariates_file) if($covariates_file ne "");

    my @samples;
	my %hashMethyl = extract_methyl_snv_file("", $methyl_file, "", \@samples); # 提取数据

	eqtl($output_dir,$methyl_file, $snp_file, $covariates_file);

	my %hashQTL = txt2excel("$output_dir/eqtl.out", "$output_dir/meQTL.xlsx"); #输出excel

	if($genotype_file ne ""){		
	    my %hashGeno = extract_methyl_snv_file("", $genotype_file, ""); # 提取基因型数据
	    plot_boxplot(\@samples, \%hashQTL, \%hashMethyl, \%hashGeno, $output_dir);
    }
}

######################################################################################## 子函数
sub plot_boxplot{
	my $samples       = shift @_;
	my $hashQTL       = shift @_;
	my $hashMethyl    = shift @_;
	my $hashGeno      = shift @_;
	my $output_dir    = shift @_;
    my $plot_data_dir = "$output_dir/plot_data";
    mkdir $plot_data_dir if(not -e $plot_data_dir);    
    my $pdf  = "$output_dir/SNV_Methyl.boxplot.pdf";
    print "###Start Ploting $pdf...";
	my $R = Statistics::R->new(bin => $RBIN);
       $R->startR;
       $R->send(qq` .libPaths("$RLIB") \n`);
       $R->send(qq` pdf("$pdf") \n`) ;

    my $count = 0;
    foreach my $pair(sort { $hashQTL->{$a}{'p-value'} <=> $hashQTL->{$b}{'p-value'} } keys %$hashQTL){
        $count++;
        last if($hashQTL->{$pair}{'p-value'} < 0.05); # 默认只绘制p值<0.05的
        last if($count > 300);                        # 默认最多只绘制300个
        my ($snp_id, $methyl_pos) = split/,,,/,$pair;
        open PLOT, ">$plot_data_dir/$pair.data";
        print PLOT "Sample\tMethyl\tGeno\n";
        foreach my $sample(@$samples){
            my $methyl = $hashMethyl->{$methyl_pos}{$sample};
            my $geno   = $hashGeno->{$snp_id}{$sample};
            next if($methyl !~ /\d/ or $geno !~/\w/);
            print PLOT "$sample\t$methyl\t$geno\n";
        }
        close PLOT;
        
        my $pvalue = $hashQTL->{$pair}{'p-value'};
        $R->send(qq` data=read.table("$plot_data_dir/$pair.data", sep="\t", header=TRUE) \n`) ;    	
    	$R->send(qq` boxplot(Methyl~Geno, data, xlab = "$snp_id Genotype",ylab = "Methylation", main="$snp_id-$methyl_pos") \n`) ;
    	$R->send(qq` mtext("p.value = $pvalue", 1, line=2, adj=1, col="blue") \n`) ;
    }
    $R->send(qq` dev.off \n`) ;	
}

sub txt2excel{
	my $txt = shift @_;
	my $xlsx = shift @_;
	my %hashQTL;
	print "###Start Generating $xlsx...";
	my $workbook = Excel::Writer::XLSX->new($xlsx);
    my %format   = format_run($workbook);
    my $sheet    = $workbook->add_worksheet("meQTL");
    my $row      = 0;
    my @titles   = ();
	open TXT, $txt;
	while(<TXT>){
		$_=~s/[\r\n]//g;
		my @datas = split/\t/, $_;
		if($row==0){
			@titles = @datas;
            foreach my $col(0..$#datas){
            	$sheet->write($row, $col ,$datas[$col], $format{'title'});
            }
            $row++;
            next;
		}
		foreach my $col(0..$#datas){
			my $color = "normal";
			   $color = "orange" if($titles[$col] eq "p-value" and $datas[$col]<0.05);
            $sheet->write($row, $col ,$datas[$col], $format{$color});

            my $pair = "$datas[0],,,$datas[1]";
            $hashQTL{$pair}{$titles[$col]} = (defined $datas[$col]) ? $datas[$col] : "";
        }
		$row++;
	}
	close TXT;
    ReadMe($workbook, \%format, "readme.txt");
    print "OK\n";
    return %hashQTL;
}

sub eqtl{
	my $output_dir     = shift @_;
	my $methyfile      = shift @_;
	my $snvfile        = shift @_;
	my $covariatesfile = shift @_;
	   $covariatesfile = "" if($covariatesfile=~/\w/ and not -e $covariatesfile);
	my $output         = "$output_dir/eqtl.out";
	print "###Start generating $output...";
	my $R = Statistics::R->new(bin => $RBIN);
       $R->startR;
       $R->run( qq` .libPaths("$RLIB") ` );
       $R->run( qq` library("MatrixEQTL")   ` );
    $R->run( qq` SNP_file_name = "$snvfile"  `) ;
	$R->run( qq` expression_file_name = "$methyfile"  `) ;# 表达量测试数据
	$R->run( qq` covariates_file_name =  "$covariatesfile"  `) ;#协变量，如果没有协变量，需要把该变量字符设定为空
	$R->run( qq` output_file_name = "$output"  `) ;#输出文件
	$R->run( qq` pvOutputThreshold = 1  `) ;# p值阈值，只输出小于该值的结果，节省空间
	$R->run( qq` errorCovariance = numeric()  `) ;#默认即可
	# 数据读入，三种数据采用同样的规则读入
	$R->run( qq` snps = SlicedData\$new()  `) ;
	$R->run( qq` snps\$fileDelimiter = "\t"  `) ;     # the TAB character
	$R->run( qq` snps\$fileOmitCharacters = "NA"  `) ;# denote missing values;
	$R->run( qq` snps\$fileSkipRows = 1  `) ;          # one row of column labels
	$R->run( qq` snps\$fileSkipColumns = 1  `) ;      # one column of row labels
	$R->run( qq` snps\$fileSliceSize = 2000  `) ;      # read file in pieces of 2,000 rows
	$R->run( qq` snps\$LoadFile( SNP_file_name )  `) ;
	#
	$R->run( qq` gene = SlicedData\$new()  `) ;
	$R->run( qq` gene\$fileDelimiter = "\t"  `) ;     # the TAB character
	$R->run( qq` gene\$fileOmitCharacters = "NA"  `) ; # denote missing values;
	$R->run( qq` gene\$fileSkipRows = 1  `) ;          # one row of column labels
	$R->run( qq` gene\$fileSkipColumns = 1  `) ;       # one column of row labels
	$R->run( qq` gene\$fileSliceSize = 2000  `) ;     # read file in pieces of 2,000 rows
	$R->run( qq` gene\$LoadFile( expression_file_name )  `) ;
	#
	$R->run( qq` cvrt= SlicedData\$new()  `) ;
	$R->run( qq` cvrt\$fileDelimiter = "\t"  `) if($covariatesfile ne "");      # the TAB character
	$R->run( qq` cvrt\$fileOmitCharacters = "NA"  `) if($covariatesfile ne ""); # denote missing values;
	$R->run( qq` cvrt\$fileSkipRows = 1  `) if($covariatesfile ne "");         # one row of column labels
	$R->run( qq` cvrt\$fileSkipColumns = 1  `) if($covariatesfile ne "");      # one column of row labels
	$R->run( qq` cvrt\$fileSliceSize = 2000  `) if($covariatesfile ne "");     # read file in pieces of 2,000 rows
	$R->run( qq` cvrt\$LoadFile( covariates_file_name )  `) if($covariatesfile ne "");
	# 开始差异化分析
	$R->run( qq` me = Matrix_eQTL_engine(snps = snps,
		                                 gene = gene,
		                                 cvrt = cvrt,
		                                 output_file_name = output_file_name,
		                                 pvOutputThreshold = pvOutputThreshold,
                                         useModel = modelLINEAR, 
                                         errorCovariance = errorCovariance, 
                                         verbose = TRUE,
                                         pvalue.hist = TRUE,
                                         min.pv.by.genesnp = FALSE,
                                         noFDRsaveMemory = FALSE) `) ;
	print "OK\n";
}

sub check_sample_names{
	my $methyl_file = shift @_;
	my $snp_file    = shift @_;
	my $methyl_head = `cat $methyl_file| head -n1`;
	   $methyl_head =~s/[\r\n]//g;
	my $snp_head    = `cat $snp_file| head -n1`;
	   $snp_head    =~s/[\r\n]//g;
	my @methyl_samples = split/\t/,$methyl_head;
	my @snp_samples    = split/\t/,$snp_head;
	shift @methyl_samples;
	shift @snp_samples;
	if(!(@methyl_samples~~@snp_samples)){
		die "Error : methyl_file and snp_file/covariates_file's samples NOT match!\n";
	}
}
sub extract_methyl_snv_file{
    my $samples     = shift @_;
    my $input_file  = shift @_;
    my $output_file = shift @_;
    my $allsamples  = "";
       $allsamples  = shift @_ if(exists $_[0]);
    my $row = 0;
    my @titles = ();
    my %hashInfo = ();
    open IN,$input_file;
    while(<IN>){
    	$_=~s/[\r\n]//g;
        my @datas = split/\t/, $_;
        if($row==0){
            @titles = @datas;
            $row++;
            next;
        }
        foreach my $i(0..@titles-1){
            $hashInfo{$datas[0]}{$titles[$i]} = (defined $datas[$i])? $datas[$i] : "";
        }
        $row++;
    }
    close IN;

    if($allsamples ne ""){
    	@$allsamples = @titles;shift @$allsamples;
        return %hashInfo ; #不需要输出下列中间文件
    }    

    my $sample_string = join"\t", @$samples;
    open OUT,">$output_file";
    print OUT "title\t$sample_string\n";
    foreach my $name(sort keys %hashInfo){
    	my $outString = $name;
    	foreach my $sample(@$samples){
    		if(not exists $hashInfo{$name}{$sample}){
    			die "Error : $sample isn't in your $input_file\n";
    		}
    		my $value = $hashInfo{$name}{$sample};
    		$outString .= "\t$value";
    	}
    	print OUT "$outString\n";
    }
    close OUT;
    return %hashInfo;
}

sub readGroup{
	my $group_file = shift @_;
	my $row = 0;
	my @groups = ();
	my %hashGroup = ();
	open GROUP, $group_file;
	while(<GROUP>){
        $_=~s/[\r\n]//g;
        my @datas = split/\t/, $_;
        if($row==0){
            @groups = @datas;
            $row++;
            next;
        }
        foreach my $i(0..@groups-1){
            $hashGroup{$groups[$i]} .= "$datas[$i]," if(defined $datas[$i] and $datas[$i] =~/\w/);
        }
        $row++;
	}
	close GROUP;
	return %hashGroup;
}

sub help{
	my $info = "
Program: MethylTarget与SNP的meQTL相关性分析
Version: 2019-4-27
Contact: 230 钱妤雯

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
         
         --methyl_file/-m          甲基化文件(行为位点,列为样本,缺失值为空)
         --snp_file/-s             突变表格(行为位点,列为样本,值为0/1/2,高频碱基纯合设为0,低频碱基纯合设为2,杂合型设为1,缺失值为空)
         --genotype_file/-geno     原始突变基因型表格(行为位点,列为样本,值为基因型,缺失值为空)
         --group_file/-g           样本分组文件(每一列为一组样本,列名为组名)
         --covariates_file/-c      样本协变量文件(行为变量名,列为样本,缺失值为空)
         --output_dir/-o           结果输出目录
         --help/-h                 查看帮助文档
	\n";
	return $info;
}

sub format_run{
	my ($workbook)=@_;
    my %format=();
    $format{'title'} = $workbook->add_format();
    $format{'title'} ->set_align('center');
    $format{'title'} ->set_align('vcenter');
    $format{'title'} ->set_size(12);
    $format{'title'} ->set_font("Times New Roman");
    $format{'title'} ->set_border();
    $format{'title'} ->set_bg_color("yellow");
    $format{'title'} ->set_color("black");
    
    $format{'normal'} = $workbook->add_format();
    $format{'normal'} ->set_align('center');
    $format{'normal'} ->set_align('vcenter');
    $format{'normal'} ->set_size(12);
    $format{'normal'} ->set_font("Times New Roman");
    $format{'normal'} ->set_border();
    
    $format{'orange'} = $workbook->add_format();
    $format{'orange'} ->set_align('vcenter');
    $format{'orange'} ->set_size(12);
    $format{'orange'} ->set_font("Times New Roman");
    $format{'orange'} ->set_bg_color("#fac090");
    $format{'orange'} ->set_border();
   
    $format{'readme1'} = $workbook->add_format();
    $format{'readme1'}->set_align('center');
    $format{'readme1'}->set_align('vcenter');
    $format{'readme1'}->set_bold();
    $format{'readme1'}->set_size(14);
    $format{'readme1'}->set_font("Times New Roman");
    $format{'readme1'}->set_border();

    $format{'readme2'} = $workbook->add_format();
    $format{'readme2'}->set_align('vcenter');
    $format{'readme2'}->set_bold();
    $format{'readme2'}->set_size(14);
    $format{'readme2'}->set_font("Times New Roman");

    $format{'readme2tmp'} = $workbook->add_format();
    $format{'readme2tmp'}->set_right();

    $format{'readme3'} = $workbook->add_format();
    $format{'readme3'}->set_align('center');
    $format{'readme3'}->set_align('vcenter');
    $format{'readme3'}->set_bold();
    $format{'readme3'}->set_size(11);
    $format{'readme3'}->set_font("Times New Roman");
    $format{'readme3'}->set_border();

    $format{'readme4'} = $workbook->add_format();
    $format{'readme4'}->set_align('vcenter');
    $format{'readme4'}->set_bold();
    $format{'readme4'}->set_size(11);
    $format{'readme4'}->set_font("Times New Roman");
    $format{'readme4'}->set_border();

    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    return %format;
}

sub ReadMe{
	my ($workbook,$format,$file)=@_;
	my $sheet=$workbook->add_worksheet("Read Me");
	$sheet->set_row(0, 65);
	$sheet->set_column('A:A', 35);
	$sheet->set_column('B:B', 25);
	$sheet->set_column('C:C', 110);
	my $row=0;
	open FILE,$file;
	while(<FILE>){
		$_=~ s/\^/\n/g;
		my @split_line=split /\t/,$_;
		my $col=0;
		foreach(@split_line)
		{
			my $text=decode("UTF-8",$_);
			if($row==0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme1'});}
			if($row==0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme2'});}
			if($row==0 and $col==2){$sheet->write($row,$col,$text,$format->{'readme2tmp'});}
			if($row>0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme3'});}
			if($row>0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme4'});}
			if($row>0 and $col==2){$sheet->write($row,$col,$text,$format->{'readme5'});}
			$col++;
		}
		$row++;
	}
	close FILE;
}
