use strict;
use warnings;
$|=1;
my $table2excel = "/home/pub/bin/NGS/chip/GATK4/tools/personal/table2excel.pl";
my $tassel = "/home/genesky/software/tassel/5/run_pipeline.pl";
my $input_dir = "/home/ganb/work/STR/20B0309B";
my $output_dir = "$input_dir/run_result";
my %hashData =  get_data();
make_dir($output_dir);

# 开始分析
foreach my $data(sort keys %hashData)
{
    my $dir = "$output_dir/$data";
    make_dir($dir);

    my $snp = $hashData{$data}{'snp'};
    my $pheno = $hashData{$data}{'pheno'};
    my $log = "$dir/run.log";

    # glm分析
    system("perl $tassel -fork1 -t $pheno -fork2 -h $snp -combine3 -input1 -input2 -intersect -FixedEffectLMPlugin -endPlugin -export $dir/glm_output 2>&1 |tee -a $log");

    # kingship
    my $kingship       = "$dir/king.txt";
    my $kingship_clean = "$dir/king.clean.txt";

    system("perl $tassel -fork1 -h $snp -KinshipPlugin  -method Normalized_IBS  -endPlugin -export $dir/king -exportType SqrMatrix 2>&1 |tee -a $log");

    # 表头特殊处理
    system("sed '1,2d' $kingship > $kingship_clean");
    my $sample_name = `cut -f 1 $kingship_clean|xargs`;
    $sample_name=~s/[\r\n]//g;
    my @sample_names = split /\s+/, $sample_name;
    my $head_new = join "\t", ('sample', @sample_names);
       $head_new.="\n";
    system("sed -i '1i\\$head_new' $kingship_clean");

    # mlm分析
    system("perl $tassel -fork1 -t $pheno -fork2 -h $snp -fork3 -k $dir/king.txt -combine4 -input1 -input2 -intersect -combine5 -input3 -input4 -mlm -export $dir/mlm_output 2>&1 |tee -a $log");

}

# 结果输出
my $readme = "./readme.txt";
foreach my $data(sort keys %hashData)
{   
    my $dir = "$output_dir/$data";
    
    my $glm1           = "$dir/glm_output1.txt";
    my $glm2           = "$dir/glm_output2.txt";
    my $kingship_clean = "$dir/king.clean.txt";
    my ($mlm1, $mlm2)  = get_mlm_result($dir, "mlm_output");
    
    print "output $output_dir/$data.xlsx\n";
    system("perl $table2excel -i $glm1,$glm2,$kingship_clean,$mlm1,$mlm2,$readme -s GLM1,GLM2,kingship,MLM1,MLM2,readme -o $output_dir/$data.xlsx");
}

sub get_mlm_result{
    my $dir    = shift @_;
    my $prefix = shift @_;
    my @files = glob "$dir/$prefix*";
    my %hashFile;
    foreach my $file(@files)
    {
        my ($id) = $file=~/$prefix(\d+)\.txt/;
        $hashFile{$id} = $file;
    }
    my @ids = sort {$b<=>$a} keys %hashFile;
    my $mlm1 = $hashFile{$ids[1]};
    my $mlm2 = $hashFile{$ids[0]};

    return ($mlm1, $mlm2);
}


sub get_data{
    my %hashData;
    # 耐盐数据
    $hashData{'data_set1_salt_snp_form'}{'snp'} = "$input_dir/data_set1_salt/output/snp_form.snp.hmp.txt";
    $hashData{'data_set1_salt_snp_form'}{'pheno'} = "$input_dir/data_set1_salt/output/snp_form.pheno.txt";

    $hashData{'data_set1_salt_snp_physiology'}{'snp'} = "$input_dir/data_set1_salt/output/snp_physiology.snp.hmp.txt";
    $hashData{'data_set1_salt_snp_physiology'}{'pheno'} = "$input_dir/data_set1_salt/output/snp_physiology.pheno.txt";

    $hashData{'data_set1_salt_ssr_form'}{'snp'} = "$input_dir/data_set1_salt/output/ssr_form.snp.hmp.txt";
    $hashData{'data_set1_salt_ssr_form'}{'pheno'} = "$input_dir/data_set1_salt/output/ssr_form.pheno.txt";

    $hashData{'data_set1_salt_ssr_physiology'}{'snp'} = "$input_dir/data_set1_salt/output/ssr_physiology.snp.hmp.txt";
    $hashData{'data_set1_salt_ssr_physiology'}{'pheno'} = "$input_dir/data_set1_salt/output/ssr_physiology.pheno.txt";

    # 耐冷数据
    $hashData{'data_set2_cold_snp_form'}{'snp'} = "$input_dir/data_set2_cold/output/snp_form.snp.hmp.txt";
    $hashData{'data_set2_cold_snp_form'}{'pheno'} = "$input_dir/data_set2_cold/output/snp_form.pheno.txt";

    $hashData{'data_set2_cold_snp_physiology'}{'snp'} = "$input_dir/data_set2_cold/output/snp_physiology.snp.hmp.txt";
    $hashData{'data_set2_cold_snp_physiology'}{'pheno'} = "$input_dir/data_set2_cold/output/snp_physiology.pheno.txt";

    $hashData{'data_set2_cold_ssr_form'}{'snp'} = "$input_dir/data_set2_cold/output/ssr_form.snp.hmp.txt";
    $hashData{'data_set2_cold_ssr_form'}{'pheno'} = "$input_dir/data_set2_cold/output/ssr_form.pheno.txt";

    $hashData{'data_set2_cold_ssr_physiology'}{'snp'} = "$input_dir/data_set2_cold/output/ssr_physiology.snp.hmp.txt";
    $hashData{'data_set2_cold_ssr_physiology'}{'pheno'} = "$input_dir/data_set2_cold/output/ssr_physiology.pheno.txt";

    foreach my $data(sort keys %hashData)
    {
        foreach my $type(keys %{$hashData{$data}})
        {
            my $file = $hashData{$data}{$type};
            die " lost file $file\n" if(is_file_ok($file) == 0);
        }
    }
    return %hashData;
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
# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}