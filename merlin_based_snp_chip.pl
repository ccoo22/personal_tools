$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
merlin 分析
Version: v1.0 2020-04-21
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
my $DEFAULT_CALLRATE     = 0.95;
my $DEFAULT_MAF          = 0.01;
my $DEFAULT_NSNP         = 2;
my $DEFAULT_DISEASE_ALLELE_FREQUENCY = 0.0001;
my $DEFAULT_HOMR_PENETRANCE          = 0.0001;
my $DEFAULT_HET_PENETRANCE           = 1;
my $DEFAULT_HOMA_PENETRANCE          = 1;
my $DEFAULT_STEP                     = 1;
my $DEFAULT_CM_SPECIES               = 'hg19';

my $DEFAULT_SOFT_PLINK   = "/home/genesky/software/plink/1.07/plink";
my $DEFAULT_SOFT_PLINK19 = "/home/genesky/software/plink/1.90_beta/plink";
my $DEFAULT_SOFT_MAPTHIN = "/home/genesky/software/mapthin/v1.11/mapthin";
my $DEFAULT_SOFT_MERLIN  = "/home/genesky/software/merlin/1.1.2/executables";
my %hashCMDB = ('hg19' => "/home/genesky/database/self_build_database/gwas/genetic_map/hg19/text/chr@.b37.gmap", 
                'hg38' => "/home/genesky/database/self_build_database/gwas/genetic_map/hg38/text/chr@.b38.gmap", 
    );


# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($input_plink, $sample_relation, $output_dir, $callrate, $maf, $nsnp, $update_pedigree_only, $update_cm, $cm_species, $disease_allele_frequency, $homr_penetrance, $het_penetrance, $homa_penetrance, $step, $SOFT_PLINK, $SOFT_MAPTHIN, $SOFT_MERLIN, $if_help);
GetOptions(
	"input_plink|i=s"     => \$input_plink,
	"sample_relation|s=s" => \$sample_relation,
	"output_dir|o=s"      => \$output_dir,

	"callrate=s"          => \$callrate,
	"maf=s"               => \$maf,
	"nsnp=s"              => \$nsnp,

	"update_pedigree_only!" => \$update_pedigree_only,




	"disease_allele_frequency=s" => \$disease_allele_frequency,
	"homr_penetrance=s"          => \$homr_penetrance,
	"het_penetrance=s"           => \$het_penetrance,
	"homa_penetrance=s"          => \$homa_penetrance,
	"step=s"                     => \$step,

	"update_cm!" => \$update_cm,
	"cm_species=s" => \$cm_species,

	"plink=s"             => \$SOFT_PLINK,
	"mapthin=s"           => \$SOFT_MAPTHIN,
	"merlin=s"            => \$SOFT_MERLIN,
	"help|h"              => \$if_help,
);
die "
Options: 必填

        --input_plink/-i                原始plink格式数据输入前缀， 例如sample.ped, sample.map文件， 则输入 -i sample
                                        注意： .map文件的第三列一定要有cm数据，否则位点都mapthin软件被过滤了。如果没有，也可以通过本软件的参数 --update_cm  --cm_species 做填充
        --sample_relation/-s            ped格式样本关系矩阵, 含有表头，前6列与plink的ped文件格式一致，第7列表示当前样本在input_plink中的样本编号（因为: 有可能input_plink文件中的样本名不合适，或者需要更换样本名.如果样本在input_plink中没有，第7列空白即可）
        --output_dir/-o                 结果输出路径

Options: 可选
        --update_pedigree_only         仅更新了sample_relation中的部分数值。 当流程运行过一次后，需要再次修改数值时，可以增加该参数，流程会跳过、简化一些分析过的步骤。
                                       如果是删除了部分样本，建议重新运行整个流程。
                                       注意：不要更新样本的ID
        --disease_allele_frequency     merlin分析 致病率        (default: $DEFAULT_DISEASE_ALLELE_FREQUENCY)
        --homr_penetrance              merlin分析 野生型致病率   (default: $DEFAULT_HOMR_PENETRANCE)
        --het_penetrance               merlin分析 杂合突变致病率 (default: $DEFAULT_HET_PENETRANCE)
        --homa_penetrance              merlin分析 纯合突变致病率 (default: $DEFAULT_HOMA_PENETRANCE)
        --step                         merlin分析 步长密度, 值越大，运行越慢   (default: $DEFAULT_STEP)

        --callrate                 位点callrate过滤, 保留 > callrate 的位点 (default: $DEFAULT_CALLRATE)
        --maf                      位点maf过滤, 保留 > maf 的位点 (default: $DEFAULT_MAF)
        --nsnp                     1cm取n个snp (default: $DEFAULT_NSNP)
        --update_cm                如果原始map文件的第三列 genetic distances 为0，则需要添加该参数进行填充。否则，--nsnp 参数过滤后就没有snp了。
        --cm_species               如果声明参数  --update_cm ，则需要指定物种版本，目前支持 hg19  hg38 两种  (default: '$DEFAULT_CM_SPECIES')
        --plink                    更改软件 plink 版本 (default: '$DEFAULT_SOFT_PLINK')
        --mapthin                  更改软件 mapthin 版本 (default: '$DEFAULT_SOFT_MAPTHIN')
        --merlin                   更改软件 merlin 版本 (default: '$DEFAULT_SOFT_MERLIN')
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $input_plink or not defined $sample_relation or not defined $output_dir);

$disease_allele_frequency   = $DEFAULT_DISEASE_ALLELE_FREQUENCY if (not defined $disease_allele_frequency);
$homr_penetrance            = $DEFAULT_HOMR_PENETRANCE if (not defined $homr_penetrance);
$het_penetrance             = $DEFAULT_HET_PENETRANCE if (not defined $het_penetrance);
$homa_penetrance            = $DEFAULT_HOMA_PENETRANCE if (not defined $homa_penetrance);
$step                       = $DEFAULT_STEP if (not defined $step);

$callrate   = $DEFAULT_CALLRATE if (not defined $callrate);
$maf        = $DEFAULT_MAF if (not defined $maf);
$nsnp       = $DEFAULT_NSNP if (not defined $nsnp);
$cm_species = $DEFAULT_CM_SPECIES if (not defined $cm_species);

$SOFT_PLINK   = $DEFAULT_SOFT_PLINK if (not defined $SOFT_PLINK);
$SOFT_MAPTHIN = $DEFAULT_SOFT_MAPTHIN if (not defined $SOFT_MAPTHIN);
$SOFT_MERLIN  = $DEFAULT_SOFT_MERLIN if (not defined $SOFT_MERLIN);
 
$output_dir = File::Spec->rel2abs($output_dir);  # 转换为绝对路径
mkdir $output_dir if(not -e $output_dir);

###################################################################### 初始化
my $DATA_TIME = `date +\"\%Y-\%m-\%d \%H:\%M.\%S\"`;
my $RUN_INFO = "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
[SET] 参数 disease_allele_frequency : $disease_allele_frequency
[SET] 参数 homr_penetrance : $homr_penetrance
[SET] 参数 het_penetrance : $het_penetrance
[SET] 参数 homa_penetrance : $homa_penetrance
[SET] 参数 step : $DEFAULT_STEP
[SET] 参数 callrate : $callrate
[SET] 参数 maf : $maf
[SET] 参数 nsnp : $nsnp
[SET] 软件 plink : $SOFT_PLINK
[SET] 软件 mapthin : $SOFT_MAPTHIN
[SET] 软件 merlin : $SOFT_MERLIN
";
open SAVE, ">>$output_dir/merlin_run.info"; print SAVE $SCRIPT_INFO.$RUN_INFO; close SAVE;
print $RUN_INFO;

###################################################################### 主程序

print "读入样本家系信息\n";
my %hashRelation = read_relation($sample_relation);
my %hashRelationID = read_relation_ID($sample_relation);
 
#########
# （1）筛选样本，并替换样本信息
#########
print "\n(1) 筛选需要的样本，并替换样本信息\n";
my $select_sample = "$output_dir/1.select_sample";
if(is_file_ok("$select_sample.ped", "$select_sample.map") == 1 and  defined $update_pedigree_only)
{
    print "[warning] 数据已经存在，且声明了 [update_pedigree_only] 参数，仅更新ped中家系的参数\n";
    update_pedigree($select_sample, \%hashRelationID, 'tab_split');

}elsif(is_file_ok("$select_sample.ped", "$select_sample.map") == 1)
{
    print "[warning] 数据已经存在，跳过\n";
}else
{
    select_sample(\%hashRelation, $input_plink, $select_sample);
}
 
#########
# （2）筛选常染色体
#########
print "\n(2) 筛选常染色体数据\n";
my $autosomes = "$output_dir/2.autosomes.txt";
my $select_autosome = "$output_dir/2.select_autosome";
if(is_file_ok("$select_autosome.ped", "$select_autosome.map") == 1 and  defined $update_pedigree_only)
{
    print "[warning] 数据已经存在，且声明了 [update_pedigree_only] 参数，仅更新ped中家系的参数\n";
    update_pedigree($select_autosome, \%hashRelationID, 'space_split');

}elsif(is_file_ok("$select_autosome.ped", "$select_autosome.map") == 1)
{
    print "[warning] 数据已经存在，跳过\n"
}else
{
    open AUTOSOME, ">$autosomes"; 
    map{ print AUTOSOME "$_\t1\t300000000\t$_\n"} (1..22); 
    close AUTOSOME;
    system("$SOFT_PLINK --file $select_sample --extract $autosomes --range  --recode --out $select_autosome --noweb");
}

#########
# （3）callrate 过滤
#########
print "\n(3) callrate $callrate 过滤\n";
my $select_callrate = "$output_dir/3.select_callrate";
if(is_file_ok("$select_callrate.ped", "$select_callrate.map") == 1 and  defined $update_pedigree_only)
{
    print "[warning] 数据已经存在，且声明了 [update_pedigree_only] 参数，仅更新ped中家系的参数\n";
    update_pedigree($select_callrate, \%hashRelationID, 'space_split');

}elsif(is_file_ok("$select_callrate.ped", "$select_callrate.map") == 1)
{
    print "[warning] 数据已经存在，跳过\n"
}else
{   
    my $geno = 1 - $callrate;
    system("$SOFT_PLINK --file $select_autosome --geno $geno --write-snplist --out $select_callrate --nonfounders --noweb");
    system("$SOFT_PLINK --file $select_autosome --extract $select_callrate.snplist  --recode --out $select_callrate --noweb");
}

#########
# （4）maf 过滤
#########
print "\n(4) maf $maf 过滤\n";
my $select_maf = "$output_dir/4.select_maf";
if(is_file_ok("$select_maf.ped", "$select_maf.map") == 1 and  defined $update_pedigree_only)
{
    print "[warning] 数据已经存在，且声明了 [update_pedigree_only] 参数，仅更新ped中家系的参数\n";
    update_pedigree($select_maf, \%hashRelationID, 'space_split');

}elsif(is_file_ok("$select_maf.ped", "$select_maf.map") == 1)
{
    print "[warning] 数据已经存在，跳过\n"
}else
{
    system("$SOFT_PLINK --file $select_callrate --maf $maf --write-snplist --out $select_maf   --nonfounders --noweb");
    system("$SOFT_PLINK --file $select_callrate --extract $select_maf.snplist  --recode --out $select_maf --noweb");
}

#########
# （5）插入缺失样本，构建完整家系
#########
print "\n(5) 插入缺失样本，构建完整家系\n";
my $insert_miss = "$output_dir/5.insert_miss";
if(is_file_ok("$insert_miss.ped", "$insert_miss.map") == 1 and  defined $update_pedigree_only)
{
    print "[warning] 数据已经存在，且声明了 [update_pedigree_only] 参数，仅更新ped中家系的参数\n";
    update_pedigree($insert_miss, \%hashRelationID, 'space_split');

}elsif(is_file_ok("$insert_miss.ped", "$insert_miss.map") == 1)
{
    print "[warning] 数据已经存在，跳过\n"
}else
{
    insert_sample(\%hashRelation, $select_maf, $insert_miss);
}

#########
# (6) 孟德尔遗传错误检测, 删除异常位点，否则merlin会报错
#########
print "\n(6) 孟德尔遗传错误检测, 删除异常位点，否则merlin会报错\n";
my $select_mendel = "$output_dir/6.select_mendel";
if(is_file_ok("$select_mendel.ped", "$select_mendel.map") == 1 and  defined $update_pedigree_only)
{
    print "[warning] 数据已经存在，且声明了 [update_pedigree_only] 参数，仅更新ped中家系的参数\n";
    update_pedigree($select_mendel, \%hashRelationID, 'space_split');

}elsif(is_file_ok("$select_mendel.ped", "$select_mendel.map") == 1)
{
    print "[warning] 数据已经存在，跳过\n"
}else
{
    system("$SOFT_PLINK --file $insert_miss --mendel --out $select_mendel.mendel --nonfounders  --noweb");
    system("less $select_mendel.mendel.lmendel|awk -F ' ' '{ if(\$3 > 0) print \$2  }'|sed '1d' > $select_mendel.mendel.error.snplist");
    system("$SOFT_PLINK --file $insert_miss --exclude $select_mendel.mendel.error.snplist  --recode --out $select_mendel --noweb");
}

#########
# (7)  每1cm取 2 个SNP
#########
print "\n(7) 每1cm取 $nsnp 个SNP\n";
my $select_1cm_nsnp = "$output_dir/7.select_1cm_nsnp";
if(is_file_ok("$select_1cm_nsnp.ped", "$select_1cm_nsnp.map") == 1 and  defined $update_pedigree_only)
{
    print "[warning] 数据已经存在，且声明了 [update_pedigree_only] 参数，仅更新ped中家系的参数\n";
    update_pedigree($select_1cm_nsnp, \%hashRelationID, 'space_split');

}elsif(is_file_ok("$select_1cm_nsnp.ped", "$select_1cm_nsnp.map") == 1)
{
    print "[warning] 数据已经存在，跳过\n"
}else
{   
    my $this_input = $select_mendel;
    chdir $output_dir;
    # 更新 cm 信息
    if(defined $update_cm)
    {   
        print "更新cm信息\n";
        system("$DEFAULT_SOFT_PLINK19 --file $this_input --cm-map $hashCMDB{$cm_species} --recode --out $select_1cm_nsnp.with_cms");
        $this_input = "$select_1cm_nsnp.with_cms";
    }
    system("$SOFT_MAPTHIN -t $nsnp $this_input.map $select_1cm_nsnp.map.mapthin");
    system("cut -f2 $select_1cm_nsnp.map.mapthin > $select_1cm_nsnp.map.mapthin.snplist");
    system("$SOFT_PLINK --file $this_input --extract $select_1cm_nsnp.map.mapthin.snplist  --recode --out $select_1cm_nsnp --noweb");
    chdir PWD;
}

#########
# (8) merlin 输入文件准备
#########
print "\n(8) merlin 输入文件准备\n";
my $final = "$output_dir/8.final";
if(is_file_ok("$final.dat", "$final.map", "$final.ped") == 1 and  defined $update_pedigree_only)
{
    print "[warning] 数据已经存在，且声明了 [update_pedigree_only] 参数，仅更新ped中家系的参数\n";
    update_pedigree($final, \%hashRelationID, 'space_split');

}elsif(is_file_ok("$final.dat", "$final.map", "$final.ped") == 1)
{
    print "[warning] 数据已经存在，跳过\n"
}else
{
    build_final_ped("$select_1cm_nsnp.ped", "$final.ped");  # 生成ped文件，主要是替换-9
    system("cut -f1,2,3 $select_1cm_nsnp.map > $final.map.tmp ");
    system("sed '1i\\CHROMOSOME\\tMARKER\\tPOSITION' $final.map.tmp > $final.map  && rm $final.map.tmp");
    system(" awk '{print \"M\\t\"\$2}' $select_1cm_nsnp.map > $final.dat.tmp ");
    system(" sed '1i\\A\\tX' $final.dat.tmp > $final.dat && rm  $final.dat.tmp");
}
 
# (9) merlin 分析
print "\n(9) merlin 分析\n";
my $merlin_result = "$output_dir/merlin";
my $model         = "$output_dir/merlin.model.txt";

system("echo X	$disease_allele_frequency	$homr_penetrance,$het_penetrance,$homa_penetrance	Rare_Dominant > $model");
print "[DATA/PED 文件检测]\n";
chdir $output_dir;
system("$SOFT_MERLIN/pedstats -d $final.dat -p $final.ped  ");  # 检测.ped,.dat文件是否有问题
chdir PWD;

print "[merlin 无参分析]\n";
system("$SOFT_MERLIN/merlin -d $final.dat -p $final.ped -m $final.map --steps $step --npl --markerNames --pdf --tabulate --prefix $merlin_result.noparametric");

print "[merlin 有参分析]\n";
system("$SOFT_MERLIN/merlin -d $final.dat -p $final.ped -m $final.map --model $model --steps $step  --markerNames --pdf  --tabulate --bits 50  --prefix $merlin_result.parametric");

###################################################################### 子程序
sub update_pedigree{
    my $plink_prefix   = shift @_;
    my $hashRelationID = shift @_;
    my $split_way      = shift @_;

    open PED, "$plink_prefix.ped";
    open PED_TMP, ">$plink_prefix.ped.tmp";

    while(<PED>)
    {
        my ($fid, $iid, $pid, $mid, $sex, $pheno, $tmp) = split /\s+/, $_, 7;
        print PED_TMP "$hashRelationID->{$iid}{$split_way}\t$tmp" if($split_way eq 'tab_split');
        print PED_TMP "$hashRelationID->{$iid}{$split_way} $tmp" if($split_way eq 'space_split');
    }
    close PED;
    close PED_TMP;
    # 替换
    system("mv $plink_prefix.ped.tmp $plink_prefix.ped");
}

sub build_final_ped{
    my $input_bed = shift @_;
    my $output_bed = shift @_;
    open INPUT, $input_bed;
    open OUTPUT, ">$output_bed";
    while(<INPUT>)
    {
        my ($fid, $iid, $pid, $mid, $sex, $pheno, $tmp) = split /\s+/, $_, 7;
        $pheno = 0 if($pheno eq '-9');  # merlin 只能接受0
        print OUTPUT "$fid $iid $pid $mid $sex $pheno $tmp";
    }
    close INPUT;
    close OUTPUT;
}


sub insert_sample{
    my $hashRelation  = shift @_;
    my $select_maf    = shift @_;
    my $select_insert = shift @_;

    system("cp $select_maf.map $select_insert.map");  # map不变，只需要拷贝即可

    # 所有样本都齐全，不需要补充缺失家系样本
    if(not exists $hashRelation{'FAKEDATA'})
    {
        system("cp $select_maf.ped $select_insert.ped");
        return;
    }

    # 制作缺失样本数据
    open INPUT, "$select_maf.ped";
    my $line1 = <INPUT>;
    my @heads = split /\s+/, $line1;
    close INPUT;

    open FAKEDATA, ">$select_insert.fakesample.ped";
    foreach my $sample(sort keys $hashRelation{'FAKEDATA'})
    {
        my @datas = ($hashRelation{'FAKEDATA'}{$sample});
        foreach(6..$#heads)
        {
            push @datas, 0;
        }
        print FAKEDATA (join " ", @datas) . "\n";
    }
    close FAKEDATA;

    # 合并
    system("cat $select_maf.ped $select_insert.fakesample.ped > $select_insert.ped");

}

sub select_sample{
    my $hashRelation = shift @_;
    my $input        = shift @_;
    my $output       = shift @_;

    my %hashSample = map{ ($_, 1) } keys %{$hashRelation{'WITHDATA'}};

    my $map = "$input.map";
    my $ped = "$input.ped";
    my $map_select = "$output.map";
    my $ped_select = "$output.ped";

    system("cp $map $map_select");
    open PED, $ped;
    open PED_SELECT, ">$ped_select";
    while(<PED>)
    {
        my ($fid, $iid, $pid, $mid, $sex, $pheno, $tmp) = split /\s+/, $_, 7;
        next if(not exists $hashRelation{'WITHDATA'}{$iid});

        print PED_SELECT "$hashRelation{'WITHDATA'}{$iid}\t$tmp";
        $hashSample{$iid}++;
    }
    close PED;
    close PED_SELECT;

    my @losts = grep{ $hashSample{$_} == 1 } keys %hashSample;
    die "[ERROR] 部分样本没有在输入文件找到，请仔细确认： @losts\n" if(scalar(@losts) > 0);
}

sub read_relation{
    my $relation = shift @_;
    my %hashRelation;
    open RELATION, $relation;
    <RELATION>;

    while(<RELATION>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        my ($fid, $iid, $pid, $mid, $sex, $pheno, $truename) = split /\s+/, $_;
        if(defined $truename and $truename =~ /\w/)
        {   # 样本有数据
            $hashRelation{'WITHDATA'}{$truename} = "$fid\t$iid\t$pid\t$mid\t$sex\t$pheno";  # 原始数据时\t拆分 
        }else
        {  # 样本没有数据，需要用缺失值填充
            $hashRelation{'FAKEDATA'}{$iid} = "$fid $iid $pid $mid $sex $pheno"; # plink转换的格式都是空格拆分，这里保持一致
        }
    }
    close RELATION;
    return %hashRelation;
}

sub read_relation_ID{
    my $relation = shift @_;
    my %hashRelationID;
    open RELATION, $relation;
    <RELATION>;

    while(<RELATION>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        my ($fid, $iid, $pid, $mid, $sex, $pheno, $truename) = split /\s+/, $_;
 
        $hashRelationID{$iid}{'tab_split'} = "$fid\t$iid\t$pid\t$mid\t$sex\t$pheno";
        $hashRelationID{$iid}{'space_split'} = "$fid $iid $pid $mid $sex $pheno";
    }
    close RELATION;
    return %hashRelationID;
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

# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}

