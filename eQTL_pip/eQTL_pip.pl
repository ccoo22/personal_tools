use strict;
use warnings;
use Statistics::R;
use Excel::Writer::XLSX;
use Parallel::ForkManager;

die "Usage:perl $0 config.ini\n" if(@ARGV!=1);
my $startTime  = localtime();print "The whole programme start at $startTime\n";
my $configFile = shift @ARGV;
my %config     = readConfig($configFile);
showConfig(\%config);
my @allsamples = split/,/,$config{"Sample"};
my $snvfile    = $config{"allsnv"};
my %hashSNV;
my %hashSNVinfo;
readsnv(\%config,\@allsamples,$snvfile,\%hashSNV,\%hashSNVinfo);
my @compares   = getCompair(\%config);
foreach my $compare(@compares){
	my($group1Name,$group2Name)=split/,/,$compare;
	my @cases       = split/,/,$config{'Group'}{$group1Name};
	my @controls    = split/,/,$config{'Group'}{$group2Name};
	my @samples     = (@cases,@controls);
	my $add_sig_snv = "/home/zhangyd/work/association_analysis/18B0719E0A/report_7groups/group_yang.VS.msdel/logistic/logistic_ALL_add.assoc.logistic" if($compare eq "yang,msdel");
	   $add_sig_snv = "/home/zhangyd/work/association_analysis/18B0719E0A/report_7groups/group_yang.VS.msadd/logistic/logistic_ALL_add.assoc.logistic" if($compare eq "yang,msadd");
	   $add_sig_snv = "/home/zhangyd/work/association_analysis/18B0719E0A/report_7groups/group_yang.VS.yin/logistic/logistic_ALL_add.assoc.logistic" if($compare eq "yang,yin");
	   $add_sig_snv = "/home/zhangyd/work/association_analysis/18B0719E0A/report_7groups/group_yin.VS.msadd/logistic/logistic_ALL_add.assoc.logistic" if($compare eq "yin,msadd");
	   $add_sig_snv = "/home/zhangyd/work/association_analysis/18B0719E0A/report_7groups/group_yin.VS.msdel/logistic/logistic_ALL_add.assoc.logistic" if($compare eq "yin,msdel");
	   $add_sig_snv = "/home/zhangyd/work/association_analysis/18B0719E0A/report_7groups/group_msadd.VS.msdel/logistic/logistic_ALL_add.assoc.logistic" if($compare eq "msadd,msdel");
	   $add_sig_snv = "/home/zhangyd/work/association_analysis/18B0719E0A/report_7groups/group_jun.VS.shai/logistic/logistic_ALL_add.assoc.logistic" if($compare eq "jun,mishai");
       my @Sigsnvs     = ();
	   @Sigsnvs     = get_add_sig_snv($add_sig_snv);
	my $methyfile   = $config{"methyl"}."/".$group1Name.".VS.".$group2Name."/DMP.txt.gz";
	my $reportDir   = $config{"report"};
	my $covariatesfile = (exists $config{"covariates"})?$config{"covariates"}:"";
	my %hashMethy;
	my %hashMethyinfo;
	my %hashCovar;
	readMethy(\%config,\@samples,$methyfile,\%hashMethy,\%hashMethyinfo);
	readCovar(\@samples,$covariatesfile,\%hashCovar) if(-e $covariatesfile and -s $covariatesfile>0);
	my @runs = ($group1Name,$group2Name,$group1Name."_".$group2Name);
	# my @runs = ($group1Name."_".$group2Name);
	foreach my $run(@runs){
    	subRun($compare,$run,\%config,\%hashSNV,\%hashMethy,\%hashCovar,\%hashSNVinfo,\%hashMethyinfo,\@Sigsnvs); 
	}
	my $endTime=localtime();print "$endTime\n";
	
}


sub subRun{
	my $compare       = shift @_;
	my $run           = shift @_;
	my $config        = shift @_;
	my $hashSNV       = shift @_;
	my $hashMethy     = shift @_;
	my $hashCovar     = shift @_;
	my $hashSNVinfo   = shift @_;
	my $hashMethyinfo = shift @_;
	my $Sigsnvs       = shift @_;
	my $Threshold     = $config->{"Process"};
	my $snv_probe     = $config->{"snv_probe"};
	my($group1Name,$group2Name)=split/,/,$compare;
	my $output        = $config->{"report"}."/".$group1Name.".VS.".$group2Name;
	my $outputDir     = $output."/$run";
	my $prepareDir    = $outputDir."/prepare";
	mkdir $output if(not -e $output);
	mkdir $outputDir if(not -e $outputDir);
	mkdir $prepareDir if(not -e $prepareDir);	
	my @cases      = split/,/,$config{'Group'}{$group1Name};
	my @controls   = split/,/,$config{'Group'}{$group2Name};
	my @samples;    
	   @samples    = @cases if($run eq $group1Name);
	   @samples    = @controls if($run eq $group2Name);
	   @samples    = (@cases,@controls) if($run eq $group1Name."_".$group2Name);
	outputsnv(\@samples,$hashSNV,"$prepareDir/$run.snv",$Sigsnvs) if(not -e "$prepareDir/$run.snv");
	outputmethy(\@samples,$hashMethy,"$prepareDir/$run.methyl") if(not -e "$prepareDir/$run.methyl");
	outputcovar(\@samples,$hashCovar,"$prepareDir/$run.covar") if(-e $config->{"covariates"} and -s $config->{"covariates"}>0 and not -e "$prepareDir/$run.covar");
	my $methyNum_cutoff = $config->{"methyNum_cutoff"};
	my $snvNum_cutoff   = $config->{"snvNum_cutoff"};
	my @snvfiles   = splitfile("$prepareDir/$run.snv",$outputDir,"snv",$snvNum_cutoff);
	my @methyfiles = splitfile("$prepareDir/$run.methyl",$outputDir,"methyl",$methyNum_cutoff);
	
	my $para_covar = "";
	   $para_covar = "$prepareDir/$run.covar" if(-e $config->{"covariates"} and -s $config->{"covariates"}>0);
	my $pm = Parallel::ForkManager->new($Threshold);
	foreach my $snvfile(@snvfiles){
		foreach my $methyfile(@methyfiles){
			$pm->start and next;
			eQTL($outputDir,"$outputDir/methyl_split/$methyfile","$outputDir/snv_split/$snvfile",$para_covar);
			$pm->finish;
		}		
	}
	$pm->wait_all_children;

	my $FinalResultFile = "$outputDir/DisIn1M.eQTL.result";
	my $SigResultFile = "$outputDir/DisOver1M_Sig.eQTL.result";
	my %hashDis1M ;
	%hashDis1M = filterDis($config,$hashSNV,$hashSNVinfo,$hashMethyinfo);
	generateFinalResult($config,$outputDir,\@snvfiles,\@methyfiles,$hashSNVinfo,$hashMethyinfo,\%hashDis1M,$FinalResultFile,$SigResultFile);
	outputExcel($config,$FinalResultFile,$SigResultFile);
	plotPoint($config,$outputDir,$FinalResultFile);
	plotCircos($outputDir,$SigResultFile);
	plotBoxplot($outputDir,$snv_probe,$hashSNV,$hashSNVinfo,$hashMethy,$hashMethyinfo,\@samples,\@snvfiles,\@methyfiles) if(defined $snv_probe and $snv_probe=~/\w+/);
}

sub get_add_sig_snv{
	my $file = shift @_;
	my @Sigsnvs;
	open IN,$file;<IN>;
	while(<IN>){
		$_=~s/[\r\n]//g;
		$_=~s/^\s+//;
		my @datas = split/\s+/,$_;
		push(@Sigsnvs, $datas[1]) if($datas[$#datas]=~/\d/ and $datas[$#datas]<0.05);
	}
	close IN;
	return @Sigsnvs;
}

sub splitfile{
	my $file      = shift @_;
	my $outputDir = shift @_;
	my $type      = shift @_;
	my $posnum    = shift @_;
	my $dir       = $outputDir."/$type"."_split";
	mkdir $dir if(not -e $dir);
	my $head = `head -n 1 $file`;chomp($head);
	system("split -l $posnum $file -d $dir/$type");
	opendir(DIR,$dir);
	my @dir = readdir DIR;
	my @files;
	foreach my $filet(@dir){
		push(@files,$filet) if($filet=~/$type/);
	}
	foreach my $i(1..$#files){
		system("sed -i '1i\\$head' $dir/$files[$i]");
	}
	my $filenum = @files;
	print "###Note: $file has been splitted into $filenum files (each has $posnum lines)\n";
	return @files;
}

sub filterDis{
	my $config        = shift @_;
	my $hashSNV       = shift @_;
	my $hashSNVinfo   = shift @_;
	my $hashMethyinfo = shift @_;
	my $snv_probe_dis = (exists $config->{"snv_probe_dis"} and $config->{"snv_probe_dis"}=~/\d+/)?$config->{"snv_probe_dis"}:1;
	my $startTime     = localtime();
	print "###Start filtering SNP and Probe  (distance between $snv_probe_dis"."M) at $startTime...";
	my %hashDis;
	my $count=0;
	my @chrs = (1..22,"X","Y");
	foreach my $chr(@chrs){
		print "###Start dealing with snvs and probes in chr$chr\n";
		my @snvs;my @probes;
		foreach my $snv(sort keys %$hashSNV){
			push(@snvs,$snv) if($hashSNVinfo->{$snv}{'snv_chr'} eq $chr);
		}
		foreach my $probe(sort keys %$hashMethyinfo){
			push(@probes,$probe) if($hashMethyinfo->{$probe}{'probe_chr'} eq $chr);
		}
		next if(@snvs<1 or @probes<1);
		foreach my $snv(@snvs){
			my $snv_pos = $hashSNVinfo->{$snv}{'snv_pos'};
			   $snv_pos = (split/\-/,$snv_pos)[0] if($snv_pos=~/\-/);
			foreach my $probe(@probes){
				my $probe_pos = $hashMethyinfo->{$probe}{'probe_pos'};
				   $probe_pos = (split/\-/,$probe_pos)[0] if($probe_pos=~/\-/);
				next if(abs($snv_pos-$probe_pos)>($snv_probe_dis/2)*10**6);
				my $dis =($snv_pos-$probe_pos)/(10**6);
				$hashDis{$snv}{$probe}=$dis;
				$count++;
			}
		}		
	}
	print "OK\n";
	my $endTime = localtime();
	print "###Note: SNP and Probe (distance between 1M) pair num : $count (finish time:$endTime)\n";
	return %hashDis;
}

sub generateFinalResult{
	my $config        = shift @_;
	my $dir           = shift @_;
	my $snvfiles      = shift @_;
	my $methyfiles    = shift @_;
	my $hashSNVinfo   = shift @_;
	my $hashMethyinfo = shift @_;
	my $hashDis       = shift @_;
	my $finalfile     = shift @_;
	my $DisOver       = shift @_;
	my $over_dis_sig_cutoff = (exists $config->{"over_dis_sig_cutoff"} and $config->{"over_dis_sig_cutoff"}=~/\w+/)?$config->{"over_dis_sig_cutoff"}:"FDR|0.05";
	my ($cutoff_type,$cutoff_value)=split/\|/,$over_dis_sig_cutoff;
	my $startTime     = localtime();
	print "###Start generating FinalResult at $startTime...\n";
	my $pm = Parallel::ForkManager->new(10);
	foreach my $i(0..@$snvfiles-1){		
		foreach my $j(0..@$methyfiles-1){
			$pm->start and next;
			my $input = "$dir/eQTL_result/$$snvfiles[$i]"."_".$$methyfiles[$j].".out";
			my $output1 = "$dir/eQTL_result/$$snvfiles[$i]"."_".$$methyfiles[$j].".out.filter";
			my $output2 = "$dir/eQTL_result/$$snvfiles[$i]"."_".$$methyfiles[$j].".out.Sig.filter";
			open OUT1,">$output1";
			open OUT2,">$output2";
			print "###Start reading $input\n";
			open IN,"$input";
			<IN>;
			while(<IN>){
				$_=~s/[\r\n]//;chomp;
				my @info = split/\t/,$_;
				my($snv,$probe,$pvalue,$FDR)=@info[0,1,4,5];
				my $snvid = $snv;
				   $snvid=~s/\(\S+\)//g;
				my $outstring = "$snvid\t$hashSNVinfo->{$snv}{'snv_chr'}\t$hashSNVinfo->{$snv}{'snv_pos'}\t$hashSNVinfo->{$snv}{'snv_ref'}\t$hashSNVinfo->{$snv}{'snv_alt'}\t";
				   $outstring.="$probe\t$hashMethyinfo->{$probe}{'probe_chr'}\t$hashMethyinfo->{$probe}{'probe_pos'}\t";
				   $outstring.=join"\t",@info[2..5];
				if(exists $hashDis->{$snv}{$probe}){					
					$outstring.="\t$hashDis->{$snv}{$probe}";
					print OUT1 "$outstring\n";
				}elsif($cutoff_type eq "FDR" and $FDR<$cutoff_value){
					print OUT2 "$outstring\n";
				}elsif($cutoff_type eq "pvalue" and $pvalue<$cutoff_value){
					print OUT2 "$outstring\n";
				}        			
			}
			close IN;
			close OUT1;
			close OUT2;
			$pm->finish;
		}				
	}
	$pm->wait_all_children;
	my $head1 = "snv\tsnv_chr\tsnv_pos\tsnv_ref\tsnv_alt\tprobe\tprobe_chr\tprobe_pos\tbeta\tt_stat\tp_value\tFDR\tDistance\n";
	my $head2 = "snv\tsnv_chr\tsnv_pos\tsnv_ref\tsnv_alt\tprobe\tprobe_chr\tprobe_pos\tbeta\tt_stat\tp_value\tFDR\n";
	system("cat $dir/eQTL_result/*.out.filter > $finalfile");
	system("cat $dir/eQTL_result/*.out.Sig.filter > $DisOver");
	system("sed -i '1i\\$head1' $finalfile");
	system("sed -i '1i\\$head2' $DisOver");
	print "OK\n";
}

sub outputExcel{
	my $config          = shift @_;
	my $FinalResultFile = shift @_;
	my $SigResultFile   = shift @_;
	my $excel_sig_cutoff= (exists $config->{"excel_sig_cutoff"} and $config->{"excel_sig_cutoff"}=~/\w+/)?$config->{"excel_sig_cutoff"}:"pvalue|0.05";
	my ($cutoff_type,$cutoff_value) = split/\|/,$excel_sig_cutoff;
	my $startTime       = localtime();
	print "###Start Generating Excel at $startTime...";
	my $outputExcel     = $FinalResultFile;
	$outputExcel =~s/result/Sig.xlsx/;
	my $workbook = Excel::Writer::XLSX->new($outputExcel);
    my %format   = format_run($workbook);
    my $sheet    = $workbook->add_worksheet("DisIn1M_Sig");
    my $row      = 0;
    open IN,$FinalResultFile;
    my $title = <IN>;
       $title =~s/[\r\n]//;chomp $title;
    my @titles = split/\t/,$title;
    foreach my $i(0..$#titles){
    	$sheet->write($row,$i,$titles[$i],$format{'title'});
    }
    $row++;
    while(<IN>){
    	$_=~s/[\r\n]//;chomp;
    	my @datas = split/\t/,$_;
    	my($pvalue,$FDR) = @datas[10,11];
    	next if($cutoff_type eq "pvalue" and $pvalue>=$cutoff_value);
    	next if($cutoff_type eq "FDR" and $FDR>=$cutoff_value);
    	foreach my $col(0..$#datas){
    		$sheet->write($row,$col,$datas[$col],$format{'normal'});
    	}
    	$row++;
    }
    close IN;
    
    # trans表格输出
    my $outputExcel_trans = $SigResultFile;
	$outputExcel_trans =~s/result/.xlsx/;
	my $workbook_trans = Excel::Writer::XLSX->new($outputExcel_trans);
    my %format_trans   = format_run($workbook_trans);
    my $sheet_trans    = $workbook_trans->add_worksheet("DisOver1M_Sig");
    my $row_trans      = 0;
    open TRANS,$SigResultFile;
    my $title_trans = <TRANS>;
       $title_trans =~s/[\r\n]//;chomp $title_trans;
    my @titles_trans = split/\t/,$title_trans;
    foreach my $i(0..$#titles_trans){
    	$sheet_trans->write($row_trans,$i,$titles_trans[$i],$format_trans{'title'});
    }
    $row_trans++;
    while(<TRANS>){
    	$_=~s/[\r\n]//;chomp;
    	my @datas = split/\t/,$_;
    	foreach my $col(0..$#datas){
    		$sheet_trans->write($row_trans,$col,$datas[$col],$format_trans{'normal'});
    	}
    	$row_trans++;
    }
    close TRANS;
    print "OK\n";
}

sub plotPoint{
	my $config          = shift @_;
	my $dir             = shift @_;
	my $FinalResultFile = shift @_;
	my $snv_probe_dis   = (exists $config->{"snv_probe_dis"} and $config->{"snv_probe_dis"}=~/\d+/)?$config->{"snv_probe_dis"}:1;
	my $xlim            = $snv_probe_dis/2;
	print "###Start plotting point figure\n";
	my $Rbin="/home/genesky/software/r/3.5.1/bin/R";
    my $RLib="/home/genesky/software/r/3.5.1/lib64/R/library";
    my $R = Statistics::R->new(bin => $Rbin);
    $R->startR;
    $R->send(qq` .libPaths("$RLib") \n`) ;
    $R->send(qq` library("ggplot2") \n`) ;
    $R->send(qq` point=read.table("$FinalResultFile",sep="\t",header=TRUE) \n`) ; 
    $R->send(qq` p=ggplot(point, aes(x=Distance, y=-log10(p_value))) + geom_point() + xlim(c(-$xlim,$xlim)) + xlab("SNP to Probe distance(Mb)") + ylab("-LOG10pvalue") \n`) ;
    $R->send(qq` p=p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) \n`) ; 
  
    # $R->send(qq` p = plot(x=point\$Distance,y=-log10(point\$p_value),pch=20,xlim=c(-0.5,0.5),xlab="SNP to Probe distance(Mb)",ylab="-LOG10pvalue") \n`) ;
    $R->send(qq` pdf("$dir/Dis_pvalue.pdf") \n`) ;
    $R->send(qq` p \n`) ;
    $R->send(qq` dev.off \n`) ;
    $R->send(qq` png("$dir/Dis_pvalue.png") \n`) ;
    $R->send(qq` p \n`) ;
    $R->send(qq` dev.off \n`) ;
    $R->stopR;
}

sub plotCircos{
	my $dir           = shift @_;
	my $SigResultFile = shift @_;
	my $SVCircosKT    = "/home/pub/database/Human/hg19/circos-karyotype/karyotype.human.hg19.rmMT.txt";
	my $Circos        = "/home/ganb/soft/circos-0.69-5/bin/circos";
	my $CircosDir     = "$dir/Circos";
	mkdir $CircosDir if(not -e $CircosDir);
	my $data          = "$CircosDir/Circos.data";
	my $config        = "$CircosDir/Circos.ini";
	my $startTime     = localtime();
	print "###Start Plotting Circos at $startTime...";
	generateCircosData($SigResultFile,$data);
	open CONFIG,">$config";
	print CONFIG "karyotype = $SVCircosKT\n";# 基因组信息
	# 基因组显示配置
	print CONFIG "<ideogram> \n";
	print CONFIG "\t<spacing> \n";
	print CONFIG "\t\tdefault = 0.005r \n";# 设置圈图中染色体之间的空隙大小，以下设置为每个空隙大小为周长的 0.5% 
	print CONFIG "\t</spacing> \n";
	print CONFIG "\tradius            = 0.85r \n";# 设定 ideograms 的位置，以下设定 ideograms 在图离圆心的 90% 处 
	print CONFIG "\tthickness         = 60p   \n";# 设定 ideograms 的厚度，可以使用 r（比例关系） 或 p（像素）作为单位 
	print CONFIG "\tstroke_color      = black \n";# 设定 ideograms 轮廓的颜色及其厚度。如果没有该参数或设定其厚度为0，则表示没有轮廓
	print CONFIG "\tstroke_thickness  = 2p    \n";# 设定 ideograms 轮廓的颜色及其厚度。如果没有该参数或设定其厚度为0，则表示没有轮廓
	print CONFIG "\tshow_label        = yes   \n";# 设定是否显示 label 。 label 对应着 karyotype 文件的第 4 列
	print CONFIG "\tlabel_font        = bold  \n";# label字体
	print CONFIG "\tlabel_radius      = 1r + 130p \n";# 设定 label 的位置 
	print CONFIG "\tlabel_size        = 60    \n";# 设定 label 的大小 
	print CONFIG "\tlabel_parallel    = yes   \n";# 设定 label 的字体方向，yes 是易于浏览的方向。
	print CONFIG "\tshow_bands        = yes   \n"; 
	print CONFIG "\tfill_bands        = yes   \n"; 
	print CONFIG "\tband_transparency = 0     \n"; 
	print CONFIG "</ideogram>   \n"; 
	# 显示基因组刻度尺
	print CONFIG "chromosomes_units = 1000000   \n"; # 1u代表的长度
	print CONFIG "show_ticks        = yes       \n"; # 是否显示ticks
	print CONFIG "show_tick_labels  = yes       \n"; # 是否显示ticks标签
	print CONFIG "<ticks>   \n";  
	print CONFIG "\tradius           = 1r    \n";  # 设定 ticks 的位置 
	print CONFIG "\tcolor            = black \n";  # 设定 ticks 的颜色 
	print CONFIG "\tthickness        = 2p    \n";  # 设定 ticks 的厚度
	print CONFIG "\tmultiplier       = 1e-6  \n";  # 设定 ticks' label 的值的计算。将该刻度对应位置的值 * multiplier 得到能展示到圈图上的 label 值。
	print CONFIG "\tformat           = %d    \n";  # label 值的格式化方法; %d 表示结果为整数
	print CONFIG "\t<tick> \n";   # 添加一个刻度
	print CONFIG "\t\tspacing        = 30u    \n"; # 没30u（chromosomes_units）添加一个刻度   
	print CONFIG "\t\tsize           = 15p    \n"; #    
	print CONFIG "\t\tshow_label     = yes    \n"; #    
	print CONFIG "\t\tlabel_font     = bold    \n"; #   
	print CONFIG "\t\tlabel_size     = 30p    \n"; #    
	print CONFIG "\t\tlabel_offset   = 10p    \n"; #    
	print CONFIG "\t\tformat         = %d    \n"; #    
	print CONFIG "\t</tick>    \n";   
	print CONFIG "</ticks>  \n";

	print CONFIG "<links>	  \n"; 
	print CONFIG "\t<link>	  \n"; 
	print CONFIG "\tfile          = $data 	 \n"; 
	print CONFIG "\tradius        = 0.94r 	 \n"; # 设置 link 曲线的半径
	print CONFIG "\tbezier_radius = 0r 	     \n"; # 设置贝塞尔曲线半径，该值设大后曲线扁平，使图像不太好看。
	print CONFIG "\tcolor         = black_a4 \n"; # 设置 link 曲线的颜色 
	print CONFIG "\tthickness     = 4 	     \n"; # 设置 link 曲线的厚度 
	print CONFIG "\t</link>	  \n"; 
	print CONFIG "</links>	  \n"; 

	print CONFIG "<image> \n"; 
	print CONFIG "\t<<include etc/image.conf>>  \n"; 
	print CONFIG "</image>  \n"; 
	print CONFIG "<<include etc/colors_fonts_patterns.conf>>  \n"; 
	print CONFIG "<<include etc/housekeeping.conf>>  \n"; 
	system("$Circos -conf $config -outputdir $dir -outputfile Circos")	;
	print "OK\n";
}

sub generateCircosData{
	my $SigResultFile = shift @_;
	my $data          = shift @_;
	open IN,$SigResultFile;
	open OUT,">$data";
	<IN>;
	while(<IN>){
		$_=~s/[\r\n]//;chomp;
		my @datas = split/\t/,$_;
		my($snv_chr,$snv_pos,$probe_chr,$probe_pos,$FDR)=@datas[1,2,6,7,11];
		$snv_pos=(split/\-/,$snv_pos)[0] if($snv_pos=~/\-/);
		$probe_pos=(split/\-/,$probe_pos)[0] if($probe_pos=~/\-/);
		print OUT "hs$snv_chr\t$snv_pos\t$snv_pos\ths$probe_chr\t$probe_pos\t$probe_pos\n";
	}
	close IN;
	close OUT;
}

sub plotBoxplot{
	my $dir           = shift @_;
	my $snv_probes    = shift @_;
	my $hashSNV       = shift @_;
	my $hashSNVinfo   = shift @_;
	my $hashMethy     = shift @_;
	my $hashMethyinfo = shift @_;
	my $samples       = shift @_;
	my $snvfiles      = shift @_;
	my $methyfiles    = shift @_;	
	my $plotDir       = $dir."/snv_probe_boxplot";
	mkdir $plotDir if(not -e $plotDir);	
	my $pdf  = "$dir/SNV_Probe.boxplot.pdf";
	my $Rbin="/home/genesky/software/r/3.5.1/bin/R";
    my $RLib="/home/genesky/software/r/3.5.1/lib64/R/library";
	my $R = Statistics::R->new(bin => $Rbin);
       $R->startR;
       $R->run( qq` .libPaths("$RLib") ` );
    $R->send(qq` pdf("$pdf") \n`) ;
    my @snv_probes = split/;/,$snv_probes;
	foreach my $snv_probe(@snv_probes){
		print "###Start plotting $snv_probe boxplot...";
		my($snv,$probe)=split/\|/,$snv_probe;
		if(not exists $hashSNV->{$snv} or not exists $hashMethyinfo->{$probe}){print "###Error: Cannot find $snv or $probe\n";next;}
		open OUT,">$plotDir/$snv"."_"."$probe.plot";
		print OUT "sample\tmethyl\tgeno\n";
		my $ref = $hashSNVinfo->{$snv}{'snv_ref'};
		my $alt = $hashSNVinfo->{$snv}{'snv_alt'};
		foreach my $sample(@$samples){
			my $genotype = $hashSNV->{$snv}{$sample};
			next if($genotype eq " ");
			$genotype = $ref.$ref if($genotype eq "0");
			$genotype = $ref.$alt if($genotype eq "1");
			$genotype = $alt.$alt if($genotype eq "2");		
			foreach my $gene(sort keys %$hashMethy){
				next if(not exists $hashMethy->{$gene}{$probe}{$sample} or $hashMethy->{$gene}{$probe}{$sample}!~/\w+/);
				my $methyValue = $hashMethy->{$gene}{$probe}{$sample};
				print OUT "$sample\t$methyValue\t$genotype\n";
			}
		}
		close OUT;
		my $pvalue=-1;
		foreach my $i(0..@$snvfiles-1){	
			last if($pvalue=~/\d+/ and $pvalue!=-1);	
			foreach my $j(0..@$methyfiles-1){
				last if($pvalue=~/\d+/ and $pvalue!=-1);
				my $input = "$dir/eQTL_result/$$snvfiles[$i]"."_".$$methyfiles[$j].".out";
				$pvalue = `awk '{if(\$1=="$snv" && \$2=="$probe")print \$5}' $input`;chomp $pvalue;				
			}
		}
		my $plotFile = "$plotDir/$snv"."_"."$probe.plot";
		$R->send(qq` data=read.table("$plotFile",sep="\t",header=TRUE) \n`) ;    	
    	$R->send(qq` boxplot(methyl~geno,data,xlab="$snv Genotype",ylab="Probe Methylation",main="$snv-$probe") \n`) ;
    	$R->send(qq` mtext("p.value = $pvalue",1,line=2,adj=1,col="blue") \n`) ;
    	print "OK\n";
	}
	$R->send(qq` dev.off \n`) ;	
}

sub eQTL{
	my $outputDir = shift @_;
	my $methyfile = shift @_;
	my $snvfile   = shift @_;
	my $covariatesfile = shift @_;	
	my @tempsnv = split/\//,$snvfile;
	my @tempmethy = split/\//,$methyfile;
	my $dir = $outputDir."/eQTL_result";
	mkdir $dir if(not -e $dir);
	my $output = $dir."/".$tempsnv[$#tempsnv]."_".$tempmethy[$#tempmethy].".out";
	if(-e $output and -s $output>0){print "###Note:$output had been generated before!!!\n";return;};
	print "###Star dealing with $snvfile and $methyfile > $output\n";
	my $Rbin="/home/genesky/software/r/3.5.1/bin/R";
    my $RLib="/home/genesky/software/r/3.5.1/lib64/R/library";
	my $R = Statistics::R->new(bin => $Rbin);
       $R->startR;
       $R->run( qq` .libPaths("$RLib") ` );
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
}

sub readsnv{
	my $config   = shift @_;
	my $samples  = shift @_;
	my $file     = shift @_;
	my $hash     = shift @_;
	my $hashinfo = shift @_;
	my $snv_maf  = (exists $config->{"snv_maf"} and $config->{"snv_maf"}=~/\d+/)?$config->{"snv_maf"}:0.05;
	my $startTime=localtime();
	print "###Start reading $file at $startTime...";
	open SNV,$file;
	my $title = <SNV>;
	$title=~s/[\r\n]//;chomp($title);
	my @titles = split/\t/,$title;
	my($snvindex)=grep{$titles[$_] eq 'SNV NO.'}0..$#titles;
	my($chrindex)=grep{$titles[$_] eq 'Chrs'}0..$#titles;
	my($posindex)=grep{$titles[$_] eq 'Position'}0..$#titles;
	my($refindex)=grep{$titles[$_] eq 'Ref Allele'}0..$#titles;
	my($altindex)=grep{$titles[$_] eq 'Alt Allele'}0..$#titles;
	my($rsindex)=grep{$titles[$_] eq 'SNP ID'}0..$#titles;
	my %index;
	foreach my $sample(@$samples){
		($index{$sample})=grep{$titles[$_] eq $sample}0..$#titles;
	}
	while(<SNV>){
		$_=~s/[\r\n]//;chomp;
		my($snv,$chr,$pos,$ref,$alt,$snpid)=(split/\t/,$_)[$snvindex,$chrindex,$posindex,$refindex,$altindex,$rsindex];
		$snpid=~s/^\s+//;
		$snv = $snpid if($snpid=~/^\w+/ and $snpid!~/\(/);
		$hashinfo->{$snv}{'snv_chr'}=$chr;
		$hashinfo->{$snv}{'snv_pos'}=$pos;
		$hashinfo->{$snv}{'snv_ref'}=$ref;
		$hashinfo->{$snv}{'snv_alt'}=$alt;
		my $refcount=0;my $altcount=0;
		foreach my $sample(@$samples){
			my $type = (split/\t/,$_)[$index{$sample}];
			next if($type eq " ");
			my($sref,$salt)=split/\//,$type;
			$refcount++ if($sref eq $ref);
			$refcount++ if($salt eq $ref);
			$altcount++ if($sref eq $alt);
			$altcount++ if($salt eq $alt);
		}
		my $sum = $refcount+$altcount;
		next if($refcount<$altcount and $refcount/$sum<$snv_maf);
		next if($refcount>$altcount and $altcount/$sum<$snv_maf);
		$hash->{$snv}{'SNV'}=$snv;
		foreach my $sample(@$samples){
			my $type = (split/\t/,$_)[$index{$sample}];
			if($type eq " "){$hash->{$snv}{$sample}=" ";next;}
			my($sref,$salt)=split/\//,$type;
			$hash->{$snv}{$sample}="0" if($sref eq $ref and $salt eq $ref);
			$hash->{$snv}{$sample}="1" if($sref eq $ref and $salt eq $alt);
			$hash->{$snv}{$sample}="2" if($sref eq $alt and $salt eq $alt);
		}		
	}
	close SNV;
	print "OK\n";
}

sub readMethy{
	my $config   = shift @_;
	my $samples  = shift @_;
	my $file     = shift @_;
	my $hash     = shift @_;
	my $hashinfo = shift @_;
	my $probe_cutoff = (exists $config->{"probe_cutoff"} and $config->{"probe_cutoff"}=~/\w+/)?$config->{"probe_cutoff"}:"pvalue|0.05";
	my ($cutoff_type,$cutoff_value) = split/\|/,$probe_cutoff;
	my $startTime=localtime();
	print "###Start reading $file at $startTime...";
	open METHY,"gzip -dc $file|";
	my $title = <METHY>;
	$title=~s/[\r\n]//;chomp($title);
	my @titles = split/\t/,$title;
	my($chrindex)=grep{$titles[$_] eq 'Chr'}0..$#titles;
	my($startindex)=grep{$titles[$_] eq 'Start'}0..$#titles;
	my($endindex)=grep{$titles[$_] eq 'End'}0..$#titles;
	my($geneindex)=grep{$titles[$_] eq 'Gene'}0..$#titles;
	my($pvalueindex)=grep{$titles[$_] eq 'pvalue'}0..$#titles;
	my($fdrindex)=grep{$titles[$_] eq 'BH.adjust'}0..$#titles;
	my %index;
	foreach my $sample(@$samples){
		($index{$sample})=grep{$titles[$_] eq $sample}0..$#titles;
	}
	while(<METHY>){
		$_=~s/[\r\n]//;chomp;
		my($probe,$gene,$pvalue,$chr,$start,$end,$FDR)=(split/\t/,$_)[0,$geneindex,$pvalueindex,$chrindex,$startindex,$endindex,$fdrindex];
		next if($cutoff_type eq "pvalue" and $pvalue>=$cutoff_value);
		next if($cutoff_type eq "FDR" and $FDR>=$cutoff_value);
		$hashinfo->{$probe}{'probe_chr'}=$chr;
		$hashinfo->{$probe}{'probe_pos'}=$start;
		$hash->{$gene}{$probe}{'probe'}=$probe;
		foreach my $sample(@$samples){
			my $value = (split/\t/,$_)[$index{$sample}];
			$hash->{$gene}{$probe}{$sample}=$value;
		}
	}
	close METHY;
	print "OK\n";
}

sub readCovar{
	my $samples  = shift @_;
	my $file     = shift @_;
	my $hash     = shift @_;
	print "###Start reading $file...";
	open IN,$file;
	my $title = <IN>;
	   $title =~s/[\r\n]//;chomp $title;
	my @titles = split/\t/,$title;
	my %index;
	foreach my $sample(@$samples){
		($index{$sample})=grep{$titles[$_] eq $sample}0..$#titles;
	}
	while(<IN>){
		$_=~s/[\r\n]//;chomp;
		my @datas = split/\t/,$_;
		my $covar = $datas[0];
		foreach my $sample(@$samples){
			if(not exists $index{$sample}){print "$sample has no covariates info!\n";next;}
			my $value = $datas[$index{$sample}];
			$hash->{$covar}{$sample}=$value;
		}
	}
	close IN;
	print "OK\n";
}

sub outputsnv{
	my $samples = shift @_;
	my $hash = shift @_;
	my $output = shift @_;
	my $Sigsnvs = shift @_;
	print "###Start generating $output\n";
	open OUT,">$output";
	my $sample = join"\t",@$samples;
	my $title = "snpid\t$sample\n";
	print OUT $title;
	my $snvnum=0;
	foreach my $snv(sort keys %$hash){
		next if(!($snv~~@$Sigsnvs));
		$snvnum++;
		my $outstring = "$hash->{$snv}{'SNV'}";
		foreach my $sample(@$samples){
			$outstring.="\t$hash->{$snv}{$sample}";
		}
		$outstring.="\n";
		print OUT $outstring;
	}
	close OUT;
	return $snvnum;
}

sub outputmethy{
	my $samples = shift @_;
	my $hash = shift @_;
	my $output = shift @_;
	print "###Start generating $output\n";
	open OUT,">$output";
	my $sample = join"\t",@$samples;
	my $title = "probe\t$sample\n";
	print OUT $title;
	foreach my $gene(sort keys %$hash){
		foreach my $probe(sort keys %{$hash->{$gene}}){
			my $outstring = "$hash->{$gene}{$probe}{'probe'}";
			foreach my $sample(@$samples){
				$outstring.="\t$hash->{$gene}{$probe}{$sample}";
			}
			$outstring.="\n";
			print OUT $outstring;
		}
	}
	close OUT;	
}

sub outputcovar{
	my $samples = shift @_;
	my $hash = shift @_;
	my $output = shift @_;
	print "###Start generating $output\n";
	open OUT,">$output";
	my $sample = join"\t",@$samples;
	my $title = "sample\t$sample\n";
	print OUT $title;
	foreach my $covar(sort keys %$hash){
		my $outstring = $covar;
		foreach my $sample(@$samples){
			$outstring.="\t$hash->{$covar}{$sample}";
		}
		$outstring.="\n";
		print OUT $outstring;
	}
	close OUT;
}

sub readConfig{
    my $configFile=shift @_;
    die "Lost $configFile !\n" if(not -e $configFile or -s $configFile < 0);
    my %config;
    open IN,$configFile;
    my $count=0;
    while(<IN>){
        $_=~s/[\r\n\s]//g;
        $count++;
        next if($_=~/^#/ or $_!~/\w/);
        my @datas=split /=/,$_;
        if($datas[0] eq "Group"){
        	$config{$datas[0]}{$datas[1]}=$datas[2];
        }elsif($datas[0] eq "Compare"){
        	$config{$datas[0]}{$datas[1]}=$count;
        }else{
           $config{$datas[0]}=$datas[1];
        }
    }
    close IN;
    checkConfig(\%config);# 检查配置文件样本是否有问题
    return %config;
}

sub checkConfig{
	my $config  = shift @_;
	my %hashSample;
	map{ $hashSample{$_}++; }(split /,/,$config->{'Sample'});
	foreach my $groupName(keys %{$config->{'Group'}}){
		my @groupSamples = split /,/,$config->{'Group'}{$groupName};
		my @losts;
		map{ push @losts,$_ if(!exists($hashSample{$_})); } @groupSamples;
		die "Sample Lost:$groupName=@losts;\n" if(@losts>0);# 分组里的样本不存在
	}
}

sub getCompair{
	my $config   = shift @_;
	my @compares = sort {$config->{'Compare'}{$a}<=>$config->{'Compare'}{$b}} keys %{$config->{'Compare'}};# 差异化分析分组
    return @compares;
}

sub showConfig{
	my $config      = shift @_;
	my $sample      = (exists($config->{"Sample"}))      ? $config->{"Sample"} : "[NA]";
	my $methylDir   = (exists($config->{"methyl"}))      ? $config->{"methyl"} : "[NA]";
	my $allsnv      = (exists($config->{"allsnv"}))      ? $config->{"allsnv"} : "[NA]";
	my $report      = (exists($config->{"report"}))      ? $config->{"report"} : "[NA]";

    my @samples = split/,/,$sample;
    my $sample_num = @samples;
	print "Sample:          $sample\n";
    print "SampleNum:       $sample_num\n";
	print "ReportDir:       $report\n";
	foreach my $groupName(keys %{$config->{'Group'}}){
		print "Group = $groupName = $config->{'Group'}{$groupName}\n";
	}
	foreach my $compare(sort {$config->{'Compare'}{$a}<=>$config->{'Compare'}{$b}} keys %{$config->{'Compare'}}){
		print "Compare = $compare\n";
	}
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
	
	return %format;
}