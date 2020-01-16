#!usr/bin/perl -w
use strict;
use warnings;
use Bio::SearchIO;
use Parallel::ForkManager;

###############################
# 识别PCR的R1/R2的merge序列的突变流程
#
# 软件功能详细说明
# （1）分析前提要求R1/R2已合并，即存在 output/sample/flash/sample.extendedFrags.fa
# （2）突变统计的结果见 report/status/sample.mutation.statistics.txt，每行四列数据依次为序列ID、序列、以序列自身作为描述标准的突变策略、以参考序列作为描述标准的突变策略
# （3）突变策略记录格式为 碱基位置（不计入比对的gap）:突变类型（错配为M、插入为I、缺失为D）:突变碱基情况，同一序列多个突变之间以逗号分隔
# （4）突变统计结果只记录存在突变的序列，且相同的序列只保留一条；同一条序列存在多个不同的比对情况的，只保留匹配长度最大情况下的突变策略
# （5）blast比对结果中，与参照片段方向相反的片段，实际是匹配上参考序列的反向互补序列，此类片段的突变策略仍以参考序列作为描述标准，而非参考序列的反向互补序列，且在突变策略起始位置以Trans作为标志


system("clear");
die "Usage perl $0 Config\n" if(@ARGV!=1);
my $hashConfig   = shift @ARGV;
my %hashConfig   = read_config($hashConfig);      ###读取配置文件
my $report_dir   = $hashConfig{'Report'}; 
my $output_dir   = $hashConfig{'Output'}; 
my $target_fasta = $hashConfig{"TargetFasta"};

die "Report Dir Error\n" if(not -e $report_dir or not -d $report_dir);
die "Output Dir Error\n" if(not -e $output_dir or not -d $output_dir);
die "Target Sequence lost" if(not -e $target_fasta or not -f $target_fasta or -s $target_fasta == 0);

show_config(\%hashConfig);

my $result_dir = "$report_dir/status";
mkdir $result_dir if(not -e $result_dir or not -d $result_dir);

# 软件
my $BlastPlus  = "/home/ganb/soft/ncbi-blast-2.7.1+/bin/blastn";

# 获取样本列表
my @samples = get_sample(\%hashConfig, 'case', 'control');
my %hashCondition; 
foreach my $sample(@samples){
    my $extendedFrags_fa = "$output_dir/$sample/flash/$sample.extendedFrags.fa";
    my $result_blast = "$output_dir/$sample/flash/$sample.extendedFrags.blast.out";
    my $result_statistics = "$result_dir/$sample.mutation.statistics.txt";
    # 已完成该步骤
    if(is_file_ok($result_blast) and is_file_ok($result_statistics)){
        $hashCondition{"Finish"}{$sample}++;
        next;            
    }
    # 原始数据没问题
    if(is_file_ok($extendedFrags_fa)){
        $hashCondition{"Good2Run"}{$sample}++;
        next;
    }
    # 原始数据丢失
    $hashCondition{"Error"}{$sample}++;
}

print "Mutation Statics Start\n";
my @sample_runs = sort keys %{$hashCondition{'Good2Run'}};
if(@sample_runs > 0){
    my $threshold = exists $hashConfig{"Mutation_Statics"} ? $hashConfig{"Mutation_Statics"} : 10;
    my $pm = Parallel::ForkManager->new($threshold);
   $pm->run_on_start(sub{my ($pid, $sample) = @_; process_bar_array($sample, \@sample_runs)});# 进度条
    foreach my $sample(@sample_runs){
       $pm->start($sample) and next;
        run_Mutation_Statics($output_dir, $result_dir, $target_fasta, $sample, $BlastPlus);
        $pm->finish;    
    }
    $pm->wait_all_children;
    print "Process Completed\n"
}else{
    print "[Note] Process None, for reProcess Delete Result\n";
}


sub run_Mutation_Statics{
    my $output_dir   = shift @_;
    my $result_dir   = shift @_;
    my $target_fasta = shift @_;
    my $sample       = shift @_;
    my $BlastPlus    = shift @_;
    my $extendedFrags_fa  = "$output_dir/$sample/flash/$sample.extendedFrags.fa";
    my $result_blast      = "$output_dir/$sample/flash/$sample.extendedFrags.blast.out";
    my $result_statistics = "$result_dir/$sample.mutation.statistics.txt";


    ### blast
    system ("$BlastPlus -task blastn -query $extendedFrags_fa -evalue 0.00001 -db $target_fasta -out $result_blast -num_threads 4 -num_alignments 2");

    my @query_name;
    my %hashcount;           ###计数
    my %hashdirection;       ###序列方向
    my %hashpos_query;       ###突变信息(query序列)
    my %hashpos_subject;      ###突变信息(subject序列)
    my $searchio = new Bio::SearchIO(-format => "blast",
                                    -file   => "$result_blast");
    while(my $result = $searchio -> next_result){
        while(my $hit = $result -> next_hit){
            while(my $hsp = $hit -> next_hsp){
                my $align            = $hsp->homology_string; ###比对情况
                my $query_string     = $hsp->query_string;    ###query序列
                my $subject_string   = $hsp->hit_string;      ###subject序列
                my $query_name       = $result->query_name;   ###query序列ID
                my $query_string_key = $query_string;
                $query_string_key    =~ s/-//g;
                next if($query_name ~~ @query_name);
                push @query_name,$query_name;
                next if(exists $hashcount{$query_string_key});
                next if($align!~/\s/);
                my $subject_start    = $hsp->start('hit');   ###query序列对应在subject序列上的起始位置
                my $strand           = $hsp->strand('hit');
                ###query序列与subject序列方向一致
                if($strand == 1){
                    my @pos_mutations;
                    while($align =~ /( )/g){
                        my $pos = pos($align);
                        push @pos_mutations,$pos;
                    }
                    my @query_strings   = split//,$query_string;               ###比对序列碱基
                    my @subject_strings = split//,$subject_string;             ###目序列碱基
                    my @pos_insert;                                            ###插入位点记录
                    my @pos_delete;                                            ###缺失位点记录
                    foreach my $pos(@pos_mutations){
                        my $base_query   = $query_strings[$pos-1];             ###当前qurey碱基
                        my $base_subject = $subject_strings[$pos-1];           ###当前subject碱基
                        my $insert_num   = @pos_insert;                        ###当前位置前插入碱基个数
                        my $delete_num   = @pos_delete;                        ###当前位置前缺失碱基个数

                        ###错配
                        if($base_query =~ /[ATCGatcg]/ and $base_subject =~ /[ATCGatcg]/){
                            $hashpos_query{$query_string_key}{$pos-$delete_num}{'mismatch'}                    = $base_query;
                            $hashpos_subject{$query_string_key}{$pos+$subject_start-1-$insert_num}{'mismatch'} = $base_query;
                        }

                        ###缺失
                        if($base_query =~ /\-/ and $base_subject =~ /[ATCGatcg]/){
                            next if($pos-1 ~~ @pos_delete);
                            my $pos_delete = $pos-1;
                            my $base_delete = $base_subject;
                            push @pos_delete,$pos_delete;
                            $hashpos_subject{$query_string_key}{$pos+$subject_start-1-$insert_num}{'delete'} = $base_subject;
                            while($query_strings[$pos] =~ /\-/ and $subject_strings[$pos] =~ /[ATCGatcg]/){
                                push @pos_delete,$pos;
                                $base_delete = $base_delete.$subject_strings[$pos];
                                $hashpos_subject{$query_string_key}{$pos+1+$subject_start-1-$insert_num}{'delete'} = $subject_strings[$pos];
                                $pos++;
                            }
                            $hashpos_query{$query_string_key}{$pos_delete-$delete_num}{'delete'} = $base_delete;
                        }

                        ###插入
                        if($base_query =~ /[ATCGatcg]/ and $base_subject =~ /\-/){
                            next if($pos-1 ~~ @pos_insert);
                            my $pos_insert = $pos-1;
                            my $base_insert = $base_query;
                            push @pos_insert,$pos_insert;
                        
                            ###记录插入多个碱基情况,合并作一个插入突变进行记录
                            while($query_strings[$pos] =~ /[ATCGatcg]/ and $subject_strings[$pos] =~ /\-/){
                                push @pos_insert,$pos;
                                $base_insert = $base_insert.$query_strings[$pos];
                                $pos++;
                            }
                            $hashpos_query{$query_string_key}{$pos_insert-$delete_num}{'insert'}                    = $base_insert;
                            $hashpos_subject{$query_string_key}{$pos_insert+$subject_start-1-$insert_num}{'insert'} = $base_insert;
                        }
                    }
                    $hashcount{$query_string_key} = $query_name;
                    $hashdirection{$query_string_key} = 1;
                    next;
                }
                
                ###query序列与subject序列方向相反，此时subject序列是参考序列的反向互补序列
                elsif($strand == -1){
                    my $query_string_trans   = reverse $query_string;
                    my $subject_string_trans = reverse $subject_string;
                    my $align_trans          = reverse $align;
                    $query_string_trans   =~ tr/ATCGatcg/TAGCtagc/;            ###query序列的反向互补序列
                    $subject_string_trans =~ tr/ATCGatcg/TAGCtagc/;            ###subject序列的反向互补序列

                    my @pos_mutations;
                    while($align =~ /( )/g){
                        my $pos = pos($align);
                        push @pos_mutations,$pos;
                    }
                    my @query_strings   = split//,$query_string;               ###比对序列碱基
                    my @subject_strings = split//,$subject_string;             ###目标序列碱基
                    my @pos_insert;                                            ###插入位点记录
                    my @pos_delete;                                            ###缺失位点记录

                    foreach my $pos(@pos_mutations){
                        my $base_query   = $query_strings[$pos-1];             ###当前qurey碱基
                        my $base_subject = $subject_strings[$pos-1];           ###当前subject碱基
                        # my $insert_num   = @pos_insert;                        ###当前位置前插入碱基个数
                        my $delete_num   = @pos_delete;                        ###当前位置前缺失碱基个数

                        ###错配
                        if($base_query =~ /[ATCGatcg]/ and $base_subject =~ /[ATCGatcg]/){
                            $hashpos_query{$query_string_key}{$pos-$delete_num}{'mismatch'} = $base_query;
                        }

                        ###缺失
                        if($base_query =~ /\-/ and $base_subject =~ /[ATCGatcg]/){
                            next if($pos-1 ~~ @pos_delete);
                            my $pos_delete  = $pos-1;
                            my $base_delete = $base_subject;
                            push @pos_delete,$pos_delete;
                            while($query_strings[$pos] =~ /\-/ and $subject_strings[$pos] =~ /[ATCGatcg]/){
                                push @pos_delete,$pos;
                                $base_delete = $base_delete.$subject_strings[$pos];
                                $pos++;
                            }
                            $hashpos_query{$query_string_key}{$pos_delete-$delete_num}{'delete'} = $base_delete;
                        }

                        ###插入
                        if($base_query =~ /[ATCGatcg]/ and $base_subject =~ /\-/){
                            next if($pos-1 ~~ @pos_insert);
                            my $pos_insert = $pos-1;
                            my $base_insert = $base_query;
                            push @pos_insert,$pos_insert;
                        
                            ###记录插入多个碱基情况,合并作一个插入突变进行记录
                            while($query_strings[$pos] =~ /[ATCGatcg]/ and $subject_strings[$pos] =~ /\-/){
                                push @pos_insert,$pos;
                                $base_insert = $base_insert.$query_strings[$pos];
                                $pos++;
                            }
                            $hashpos_query{$query_string_key}{$pos_insert-$delete_num}{'insert'} = $base_insert;
                        }
                    }

                    my @pos_mutations_trans;
                    while($align_trans =~ /( )/g){
                        my $pos = pos($align_trans);
                        push @pos_mutations_trans,$pos;
                    }
                    my @query_strings_trans   = split//,$query_string_trans;               ###比对序列碱基
                    my @subject_strings_trans = split//,$subject_string_trans;             ###目序列碱基
                    my @pos_insert_trans;                                                  ###插入位点记录
                    # my @pos_delete_trans;                                                  ###缺失位点记录
                    foreach my $pos(@pos_mutations_trans){
                        my $base_query   = $query_strings_trans[$pos-1];             ###当前qurey碱基
                        my $base_subject = $subject_strings_trans[$pos-1];           ###当前subject碱基
                        my $insert_num   = @pos_insert_trans;                        ###当前位置前插入碱基个数

                        ###错配
                        if($base_query =~ /[ATCGatcg]/ and $base_subject =~ /[ATCGatcg]/){
                            $hashpos_subject{$query_string_key}{$pos+$subject_start-1-$insert_num}{'mismatch'} = $base_query;
                        }

                        ###缺失
                        if($base_query =~ /\-/ and $base_subject =~ /[ATCGatcg]/){
                            $hashpos_subject{$query_string_key}{$pos+$subject_start-1-$insert_num}{'delete'} = $base_subject;
                        }

                        ###插入
                        if($base_query =~ /[ATCGatcg]/ and $base_subject =~ /\-/){
                            next if($pos-1 ~~ @pos_insert_trans);
                            my $pos_insert = $pos-1;
                            my $base_insert = $base_query;
                            push @pos_insert_trans,$pos_insert;
                        
                            ###记录插入多个碱基情况,合并作一个插入突变进行记录
                            while($query_strings_trans[$pos] =~ /[ATCGatcg]/ and $subject_strings_trans[$pos] =~ /\-/){
                                push @pos_insert_trans,$pos;
                                $base_insert = $base_insert.$query_strings_trans[$pos];
                                $pos++;
                            }
                            $hashpos_subject{$query_string_key}{$pos_insert+$subject_start-1-$insert_num}{'insert'} = $base_insert;
                        }
                    }
                    $hashcount{$query_string_key} = $query_name;
                    $hashdirection{$query_string_key} = -1;
                    next;
                }
            }
        }
    }


    open OUT,">$result_statistics";
    foreach my $query_string_key(sort keys %hashpos_query){
        my $query_name = $hashcount{$query_string_key};
        print OUT "$query_name\t$query_string_key";
        my $mutation_string = "\t";
        ###query序列的突变信息
        foreach my $pos(sort {$a<=>$b} keys %{$hashpos_query{$query_string_key}}){
            $mutation_string = $mutation_string."$pos:M:$hashpos_query{$query_string_key}{$pos}{'mismatch'}," if(exists $hashpos_query{$query_string_key}{$pos}{'mismatch'});
            $mutation_string = $mutation_string."$pos:I:$hashpos_query{$query_string_key}{$pos}{'insert'}," if(exists $hashpos_query{$query_string_key}{$pos}{'insert'});
            $mutation_string = $mutation_string."$pos:D:$hashpos_query{$query_string_key}{$pos}{'delete'}," if(exists $hashpos_query{$query_string_key}{$pos}{'delete'});
        }
		###为正向/反向互补序列添加标志，反向互补则为Trans
        my $strand = ($hashdirection{$query_string_key} == -1)? "Trans,":"";
        $mutation_string = $mutation_string."\t$strand";
        ###参考序列的突变信息
        foreach my $pos(sort {$a<=>$b} keys %{$hashpos_subject{$query_string_key}}){
            $mutation_string = $mutation_string."$pos:M:$hashpos_subject{$query_string_key}{$pos}{'mismatch'}," if(exists $hashpos_subject{$query_string_key}{$pos}{'mismatch'});
            $mutation_string = $mutation_string."$pos:I:$hashpos_subject{$query_string_key}{$pos}{'insert'}," if(exists $hashpos_subject{$query_string_key}{$pos}{'insert'});
            $mutation_string = $mutation_string."$pos:D:$hashpos_subject{$query_string_key}{$pos}{'delete'}," if(exists $hashpos_subject{$query_string_key}{$pos}{'delete'});
        }
        print OUT "$mutation_string\n";
    }
    close OUT;
}




sub direction{
    my $start = shift @_;
    my $end   = shift @_;
    my $strand = "+";
    $strand = "-" if($end >= $start)
}



# 展示配置内容
sub show_config{
    my $hashConfig = shift @_;
    my @case_samples    = get_sample($hashConfig, "case");
    my @control_samples = get_sample($hashConfig, "control");

    my $case_num    = @case_samples;
    my $control_num = @control_samples;
    my $output_dir  = (exists $hashConfig->{"Output"})     ? $hashConfig->{"Output"}     : "[NA]";
    my $report_dir  = (exists $hashConfig->{"Report"})     ? $hashConfig->{"Report"}     : "[NA]";
 
    print "case:       $case_num\n";
    print "control:    $control_num\n";
    print "Output:     $output_dir\n";
    print "Report:     $report_dir\n";
    print "\n";

}

# 获取当前项目分析的样本
sub get_sample{
    my $hashConfig = shift @_;
    my @need_types = @_;
    my @samples;

    foreach my $need_type(@need_types)
    {
        next if(!exists($hashConfig->{$need_type}));
        foreach my $sample(split /,/, $hashConfig->{$need_type})
        {
            push @samples, $sample if($sample=~/\w/);
        }
    }

    return @samples;
}

# 读取配置文件
sub read_config{
    my $config_file = shift @_;
    my %hashConfig;

    open CONFIG, $config_file or die "[ERR] Loss Config File,$config_file\n";
    while(my $_=<CONFIG>){
        $_=~s/[\s\t\r\n]//g;
        next if($_ =~/^#/ or $_!~/\w/);
        my ($name, $value) = split /=/, $_;
        next if(!defined $value);

        if($name eq 'case' or $name eq 'control')
        {
            $hashConfig{$name} .= "$value,";
        }
        else
        {
            $hashConfig{$name} = $value;
        }
    }
    close CONFIG;
    return %hashConfig;
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

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}
