$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
代谢物分析
Version: v1.0 2020-4-15
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Parallel::ForkManager;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# Perl 自定义包
use lib SCRIPT_DIR. "../package";
use system_resource qw(check_system_resource set_thread);
use system_time qw(time_to_datetime);

# 变量
my $DEFAULT_THREAD        = 20;
my $DEFAULT_MEMORY        = 1;
my $THREAD_EXPECT         = 1;
my $DEFAULT_COMPOUND_DB   = "/home/genesky/database_new/self_build_database/metabolite_database/metaboanalyst.compound_db.rds";
my $DEFAULT_SYN_NMS       = "/home/genesky/database_new/self_build_database/metabolite_database/metaboanalyst.syn_nms.rds";
my $DEFAULT_KEGG_DB       = "/home/genesky/database_new/self_build_database/metabolite_database/kegg.hsa.rda";
my $DEFAULT_SOFT_R_SCRIPT = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $DEFAULT_SOFT_R_LIB    = "/home/genesky/software/r/3.5.1/lib64/R/library";
my $SOFT_REPLACE_MISSING  = SCRIPT_DIR."/metabolite/metabolite_replace_missing.r";  # 缺失值替换
my $SOFT_NORMALIZED       = SCRIPT_DIR."/metabolite/metabolite_normalized.r";  # 归一化
my $SOFT_NAME_MAP         = SCRIPT_DIR."/metabolite/metabolite_name_map.r";  # 代谢物名称与数据库id映射
my $SOFT_CLUSTERPROFILER  = SCRIPT_DIR."/metabolite/metabolite_clusterProfiler_metaboanalyst.r";  # kegg 富集分析
my $SOFT_PCA              = SCRIPT_DIR."/metabolite/metabolite_pca.r";  # PCA分析
my $SOFT_PLSDA            = SCRIPT_DIR."/metabolite/metabolite_plsda.r";  # plsda分析
my $SOFT_DENDROGRAM       = SCRIPT_DIR."/metabolite/metabolite_dendrogram.r";  # 进化树分析
my $SOFT_TTEST            = SCRIPT_DIR."/metabolite/metabolite_ttest.r";  # ttest分析
my $SOFT_ANOVA            = SCRIPT_DIR."/metabolite/metabolite_anova.r";  # anova分析
my $SOFT_HEATMAP          = SCRIPT_DIR."/metabolite/metabolite_heatmap.r";  # heatmap分析


# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($config_file, $metabolite_file, $sample_group, $output_dir, $already_normalized, $compound_db, $syn_nms, $kegg_db, $SOFT_R_SCRIPT, $SOFT_R_LIB, $thread, $force_cover, $no_check, $if_help);
GetOptions(
	"config|c=s"       => \$config_file,
	"metabolite|m=s"   => \$metabolite_file,
	"sample-group|s=s" => \$sample_group,
	"output-dir|o=s"   => \$output_dir,

	"already_normalized"   => \$already_normalized,
	
	"force-cover:s" => \$force_cover,
	"thread=i"      => \$thread,
	"compound_db=s" => \$compound_db,
	"syn_nms=s"     => \$syn_nms,
	"kegg_db=s"     => \$kegg_db,
	
	"Rscript=s" => \$SOFT_R_SCRIPT,
	"R_Lib=s"   => \$SOFT_R_LIB,
    "no-check"  =>\$no_check,
	"help|h"    => \$if_help,
);
die "
Options: 必填

        --config/-c                配置文件路径, 包含 Compare = Group1 : A, B 定义
        --metabolite/-m            代谢物原始数据矩阵（缺失数据填充、矫正之前的数据）。每一行是一个代谢物，每一列是一个样本，第一列是id编号，第二列是代谢物名称，第三列开始是样本实际代谢物数据。
                                   注意：请确保第一、二列的列名分别是 ID、Name
        --sample-group/-s          样本信息文件, 第1列为下机名称, 第2列为分析名称, 第3列开始为定义分组列。要有表头。
        --output-dir/-o            结果输出路径

Options: 可选
        --compound_db              代谢物注释数据库 (defaulg: $DEFAULT_COMPOUND_DB)
        --syn_nms                  代谢物同义词注释数据库，是对上一个数据库的补全 (defaulg: $DEFAULT_SYN_NMS)
        --kegg_db                  代谢物kegg数据库 (defaulg: $DEFAULT_KEGG_DB)
        --force-cover [sample,]    无参数, 删除所有原始结果文件, 完全覆盖;
                                   指定sample, 删除特定样本结果文件, 多个样本间','分隔
        --already_normalized       输入的代谢物矩阵已经做过缺失填充、矫正了。 当输入该参数时，流程跳过缺失填充、矫正步骤，直接使用输入数据做下游分析。
        --no-check                 不计算空闲CPU数量
        --thread                   使用线程数 (default: $DEFAULT_THREAD)
        --Rscript                  更改软件 Rscript 版本 (default: '$DEFAULT_SOFT_R_SCRIPT')
        --R_Lib                    更改软件 R Lib路径 (default: $DEFAULT_SOFT_R_LIB)
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $config_file or not defined $metabolite_file or not defined $sample_group or not defined $output_dir);
$thread          = $DEFAULT_THREAD if (not defined $thread);
$DEFAULT_MEMORY *= $thread/$DEFAULT_THREAD;
$compound_db     = $DEFAULT_COMPOUND_DB if (not defined $compound_db);
$syn_nms         = $DEFAULT_SYN_NMS if (not defined $syn_nms);
$kegg_db         = $DEFAULT_KEGG_DB if (not defined $kegg_db);
$SOFT_R_SCRIPT   = $DEFAULT_SOFT_R_SCRIPT if (not defined $SOFT_R_SCRIPT);
$SOFT_R_LIB      = $DEFAULT_SOFT_R_LIB if (not defined $SOFT_R_LIB);
$output_dir      = File::Spec->rel2abs($output_dir);

###################################################################### 初始化
foreach($config_file, $metabolite_file, $sample_group, $output_dir){
	die "\n[ERROR] 文件缺失, $_\n" if(not -e $_);
}
# 备份分组文件
system("cp $sample_group $output_dir/sample_group.txt");
$sample_group = "$output_dir/sample_group.txt";

my %config = read_config($config_file);  # 读入Compare信息
my %sample_group = read_sample_group($sample_group);  # 读入分组文件中每一列下每一组的样本信息（第3列开始是分组信息）
check_sample($metabolite_file, $sample_group); # 检查分组文件中的样本能否与丰度表对应上

# 分析分组
my %group;
my $id = 0;
foreach my $compare_info(sort{$config{$a} <=> $config{$b}} keys %config){
	my ($column_name, $group_list) = split /:/, $compare_info;
	die "\n[ERROR] --config 文件 Compare 格式错误, $compare_info\n" if(not defined $group_list);
	die "\n[ERROR] --sample-group 分组列名没有定义, $column_name\n" if(not exists $sample_group{$column_name});
	my @group_names = split /,/, $group_list;
	$group{$id}{"compare_name"} = "$column_name-".(join "_vs_", @group_names);
    my $sort = 0;
	foreach my $group_name(@group_names){
		die "\n[ERROR] --sample-group 分组列 $column_name 中 分组没有定义, $group_name\n" if(not exists $sample_group{$column_name}{$group_name});
		$group{$id}{"group"}{$group_name} = $sort;
		$group{$id}{"sample"}{$group_name} = $sample_group{$column_name}{$group_name};
        $sort++;
	}
	$id++;
}
foreach my $key(sort{$a <=> $b} keys %group){
	my $compare_name = $group{$key}{"compare_name"};
	die "\n[ERROR] 分组中发现样本名冲突, $compare_name\n" if(find_same_sample($group{$key}{"sample"}));
}

# 删除结果
if(defined $force_cover){
	if($force_cover =~ /\S/){
        # 查看哪个分组里的样本需要重新分析
		foreach my $key(sort{$a <=> $b} keys %group){
			my $compare_name = $group{$key}{"compare_name"};
			my @sample = map{split /,/, $group{$key}{"sample"}{$_}} keys %{$group{$key}{"sample"}};
			my $isdel = 0;
			foreach(split /,/, $force_cover){
				$isdel = 1 if($_ ~~ @sample);
			}
			if($isdel == 1){
				print "[WARNING] 覆盖结果, $compare_name\n";
				system "rm -rf $output_dir/$compare_name";
			}
		}
	}else{
		foreach my $key(sort{$a <=> $b} keys %group){
			my $compare_name = $group{$key}{"compare_name"};
			print "[WARNING] 覆盖结果, $compare_name\n";
			system "rm -rf $output_dir/$compare_name";
		}
	}
}

# 运行检测
my (@good2run, @good2run_name, @finish );
foreach my $key(sort{$a <=> $b} keys %group){
	my $compare_name = $group{$key}{"compare_name"};
	if(is_file_ok("$output_dir/$compare_name/finish.txt") == 1){  # 这里先不动他，等设计完成后，再确定任务完成的标识
		push @finish, $compare_name;
	}else{
		push @good2run, $key;
		push @good2run_name, $compare_name;
	}
}
die "\n[OK] 所有分析已运行过，finish.txt已存在，跳过\n" if(@good2run == 0);


# 运行参数
check_system_resource($thread, $DEFAULT_MEMORY) if(not defined $no_check);
my ($t_soft, $t_sample)= set_thread($thread, @good2run+0, $THREAD_EXPECT);
my $DATA_TIME = time_to_datetime(time);
my $RUN_INFO = "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
[SET] 软件 Rscript : $SOFT_R_SCRIPT
[SET] 软件 R 数据库 : $SOFT_R_LIB
[SET] 代谢物注释数据库 : $DEFAULT_COMPOUND_DB
[SET] 代谢物同义词注释数据库 : $DEFAULT_SYN_NMS
[SET] 代谢物kegg数据库 : $DEFAULT_KEGG_DB
[SET] CPU 倍数 $t_soft, 并行数量 $t_sample, CPU 总消耗 $thread, 内存总消耗 $DEFAULT_MEMORY
";
$RUN_INFO .= "\n[RUNTIME] 已完成分组, @finish\n\n" if(@finish > 0);
$RUN_INFO .= "\n[RUNTIME] 未完成分组, @good2run_name\n\n" if(@good2run_name > 0);
open SAVE, ">>$output_dir/metabolite_run.info"; print SAVE $SCRIPT_INFO.$RUN_INFO; close SAVE;
print $RUN_INFO;

###################################################################### 主程序
$ENV{"R_LIBS"} = $SOFT_R_LIB;
$DATA_TIME =~ s/\D+/\_/g;
$DATA_TIME =~ s/\D+/\_/g;
my $log_step = "$output_dir/log.step.$DATA_TIME.txt";  # 记录所有分析结果的情况，用于检查是否异常
my $log_run  = "$output_dir/log.run.$DATA_TIME.txt";   # 命令运行日志

# 根据输入数据是否经过缺失填充、矫正，选择分析
my $normalized_data = "$output_dir/metabolite_data.normalized.txt";
if(not defined $already_normalized)
{
    system "cp $metabolite_file $output_dir/metabolite_data.raw.txt";

    # 缺失值替换， 有许多方案，目前仅用一种
    my $replace_missing_data = "$output_dir/metabolite_data.replace.missing.txt";
    run_commond("[replace missing]", "$SOFT_R_SCRIPT $SOFT_REPLACE_MISSING -i $output_dir/metabolite_data.raw.txt -o $replace_missing_data", $log_run, 'new');
    check_result_file($log_step, 'new', 'replace missing', 'txt', $replace_missing_data);

    # 归一化， 有许多方案，目前仅用一种
    run_commond("[normalized]", "$SOFT_R_SCRIPT $SOFT_NORMALIZED -i  $output_dir/metabolite_data.replace.missing.txt -o $normalized_data", $log_run);
    check_result_file($log_step, 'add', 'normalized', 'txt', $normalized_data,);
}else{
    copy_normalized_data($metabolite_file, $normalized_data);  # 主要目的：替换第一、二列的名称为 ID、Name
    check_result_file($log_step, 'new', 'normalized', 'txt', $normalized_data,);
}

# 代谢物映射到数据库
my $metabolite_name            = "$output_dir/metabolite.name.txt";
my $metabolite_name_map        = "$output_dir/metabolite_name_map.txt";
my $metabolite_name_map_format = "$output_dir/metabolite_name_map.format.txt";
system "cut -f2 $normalized_data | sed '1d' > $metabolite_name";
run_commond("[metabolite name map]", "$SOFT_R_SCRIPT $SOFT_NAME_MAP -i $metabolite_name --compound_db $compound_db --syn_nms $syn_nms --output $metabolite_name_map", $log_run);
format_name_map($normalized_data, $metabolite_name_map, $metabolite_name_map_format);  # 结果整理一下
check_result_file($log_step, 'add', 'metabolite name map', 'txt', $metabolite_name_map_format,);
my $compound_count = `wc -l $metabolite_name_map_format`;
   ($compound_count) = split /\s+/, $compound_count;
   $compound_count -= 1;
 
# 分组分析
my $ForkManager = Parallel::ForkManager->new($t_sample);
foreach my $key(@good2run)
{
	my $pid = $ForkManager->start and next;
	
    my $compare_name    = $group{$key}{"compare_name"};
	my @group_names     = sort{$group{$key}{"group"}{$a} <=> $group{$key}{"group"}{$b}} keys %{$group{$key}{"group"}};
	my $compare_dir     = "$output_dir/$compare_name";  # 当前分组保存目录
	system "rm -rf $compare_dir && mkdir -p $compare_dir";
    my %hashGroup = parse_group(\@group_names, $group{$key}{"sample"}, $compare_dir);  # 分组文件拆分: (1) 所有样本 （2）两两分组的
    my $log_compare_step    = "$compare_dir/log.step.$DATA_TIME.txt";
    my $log_compare_run     = "$compare_dir/log.run.$DATA_TIME.txt";
    system("touch $log_compare_step");
    system("touch $log_compare_run");

    # 针对每一个分组文件做分析
    foreach my $sub_compare_name(sort keys %hashGroup)
    {   
        my $group_file   = $hashGroup{$sub_compare_name}{'File'};
        my $group_count  = $hashGroup{$sub_compare_name}{'GroupCount'};
        my $sample_count = $hashGroup{$sub_compare_name}{'SampleCount'};
        
        ###########
        # ttest
        ###########
        if($group_count == 2)
        {   
            my $ttest_dir = "$compare_dir/ttest/$sub_compare_name";
            system "mkdir -p $ttest_dir" if(not -e $ttest_dir);
            
            my $ttest_result = "$ttest_dir/ttest.txt";
            my $ttest_volcano = "$ttest_dir/ttest.volcano.pdf";
            run_commond("[$compare_name][$sub_compare_name] ttest", "$SOFT_R_SCRIPT  $SOFT_TTEST -i $normalized_data --sample_group_file $group_file -k 1,2 -o $ttest_result --volcano $ttest_volcano", $log_compare_run );
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] ttest", 'txt', $ttest_result, 'pdf', $ttest_volcano);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] ttest", '[Skip] group count is not equal 2'); }

        ###########
        # # anova
        ###########
        if($group_count > 2)
        {
            my $anova_dir = "$compare_dir/anova/$sub_compare_name";
            system "mkdir -p $anova_dir" if(not -e $anova_dir);

            my $anova_result = "$anova_dir/anova.txt";
            run_commond("[$compare_name][$sub_compare_name] anova", "$SOFT_R_SCRIPT  $SOFT_ANOVA -i $normalized_data --sample_group_file $group_file -k 1,2  -o $anova_result", $log_compare_run );
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] anova", 'txt', $anova_result);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] anova", '[Skip] group count is not over 2'); }

        ###########
        # # pca
        ###########
        if($sample_count > 2 and $compound_count > 2)
        {   
            my $pca_dir = "$compare_dir/pca/$sub_compare_name";
            system "mkdir -p $pca_dir" if(not -e $pca_dir);

            my $pca            = "$pca_dir/pca.pdf";
            my $pca_biplot     = "$pca_dir/pca_biplot.pdf";
            my $pca_coordinate = "$pca_dir/pca_coordinate.csv";
            run_commond("[$compare_name][$sub_compare_name] pca", "$SOFT_R_SCRIPT  $SOFT_PCA -i $normalized_data --sample_group_file $group_file --pcaplot $pca --biplot $pca_biplot --pca_coordinate $pca_coordinate", $log_compare_run );
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] pca", 'txt', $pca_coordinate, 'pdf', $pca, $pca_biplot);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] pca", '[Skip] sample count is not over 2 or compound_count is not over 2'); }

        ###########
        # # plsda
        ###########
        if($group_count > 1 and $sample_count > 2 and $compound_count > 2)
        {   
            my $plsda_dir = "$compare_dir/plsda/$sub_compare_name";
            system "mkdir -p $plsda_dir" if(not -e $plsda_dir);

            my $plsda            = "$plsda_dir/plsda.pdf";
            my $plsda_vip        = "$plsda_dir/plsda_vip.pdf";
            my $plsda_coordinate = "$plsda_dir/plsda_coordinate.txt";
            my $plsda_vip_txt    = "$plsda_dir/plsda_vip.txt";
            run_commond("[$compare_name][$sub_compare_name] plsda", "$SOFT_R_SCRIPT  $SOFT_PLSDA -i $normalized_data --sample_group_file $group_file --show_confidence --plsda $plsda --vip $plsda_vip --plsda_coordinate $plsda_coordinate --plsda_vip $plsda_vip_txt", $log_compare_run );
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] plsda", 'txt', $plsda_coordinate, $plsda_vip_txt, 'pdf', $plsda, $plsda_vip);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] plsda", '[Skip] group count is not over 1 and sample count is not over 2 or compound_count is not over 2'); }

        ###########
        # # dendrogram
        ###########
        if($sample_count > 1 and $compound_count > 1)
        {
            my $dendrogram_dir = "$compare_dir/dendrogram/$sub_compare_name";
            system "mkdir -p $dendrogram_dir" if(not -e $dendrogram_dir);

            my @dendrogram_results;
            # 距离计算方法
            foreach my $distance_method(qw(euclidean pearson spearman))  
            {   # 聚类计算方法
                foreach my $cluster_method(qw(ward.D single complete average))  
                {   
                    my $dendrogram = "$dendrogram_dir/dendrogram_$distance_method\_$cluster_method.pdf";
                    push @dendrogram_results, $dendrogram;
                    run_commond("[$compare_name][$sub_compare_name] dendrogram distance=$distance_method cluster=$cluster_method", "$SOFT_R_SCRIPT  $SOFT_DENDROGRAM -i $normalized_data --sample_group_file $group_file --clusterplot $dendrogram --smplDist $distance_method --clstDist $cluster_method", $log_compare_run );
                }
            }
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] dendrogram", 'pdf', @dendrogram_results);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] dendrogram", '[Skip] sample count is not over 1 or compound_count is not over 1'); }

        ###########
        # # heatmap
        ###########
        if($sample_count > 1 and $compound_count > 1)
        {   
            my $heatmap_dir = "$compare_dir/heatmap/$sub_compare_name";
            system "mkdir -p $heatmap_dir" if(not -e $heatmap_dir);
            my $heatmap = "$heatmap_dir/heatmap.pdf";
            run_commond("[$compare_name][$sub_compare_name] heatmap", "$SOFT_R_SCRIPT  $SOFT_HEATMAP -i $normalized_data --sample_group_file $group_file -o $heatmap --scale row --cluster_row --cluster_col --show_row --show_col", $log_compare_run );
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] heatmap", 'pdf', $heatmap);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] heatmap", '[Skip] sample count is not over 1 or compound_count is not over 1'); }

        ###########
        # # 差异代谢物ID获取：pvalue  < 0.05, vip > 1
        ###########
        my $diff_test_file = "";
           $diff_test_file = "$compare_dir/ttest/$sub_compare_name/ttest.txt" if(is_file_ok("$compare_dir/ttest/$sub_compare_name/ttest.txt") == 1);
           $diff_test_file = "$compare_dir/anova/$sub_compare_name/anova.txt" if(is_file_ok("$compare_dir/anova/$sub_compare_name/anova.txt") == 1);

        my $ttest_volcano = "$compare_dir/ttest/$sub_compare_name/ttest.volcano.pdf";

        my $plsda_vip      = "";
           $plsda_vip    = "$compare_dir/plsda/$sub_compare_name/plsda_vip.txt" if(is_file_ok("$compare_dir/plsda/$sub_compare_name/plsda_vip.txt") == 1);

        my $diff_dir = "$compare_dir/diff/$sub_compare_name";
        system "mkdir -p $diff_dir" if(not -e $diff_dir);
        my $diff_compound_id   = "$diff_dir/diff.compound.id.txt";  # 差异代谢物ID
        my $diff_compound_kegg = "$diff_dir/diff.compound.kegg.txt"; # 差异代谢物kegg ID
        my $diff_summary       = "$diff_dir/diff.txt"; # 差异信息汇总
        # 需要的参数： 代谢物ID映射文件、差异结果文件、差异类型、plsda_vip文件、输出代谢物ID、输出代谢物KEGG ID、差异结果汇总文件、pvalue阈值、vip阈值
        my $diff_compound_count = 0;
        if($diff_test_file ne "" and is_file_ok($diff_test_file) == 1 and $plsda_vip ne "" and is_file_ok($plsda_vip) == 1)
        {
            parse_diff_compound($metabolite_name_map_format, $diff_test_file, $plsda_vip, $diff_compound_id, $diff_compound_kegg, $diff_summary, 0.05, 1);
            system("cp $ttest_volcano $diff_dir/") if(is_pdf_ok($ttest_volcano) == 1);
            $diff_compound_count = `wc -l $diff_compound_id`;
            ($diff_compound_count) = split /\s+/, $diff_compound_count;
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] diff", 'txt', $diff_summary);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] diff", '[Skip] no ttest/anova file or no vip file'); }


        ###########
        # # 基于差异代谢物的绘图：pca
        ###########
        if($sample_count > 2 and $diff_compound_count > 2)
        {   
            my $pca_dir = "$compare_dir/pca/$sub_compare_name";
            system "mkdir -p $pca_dir" if(not -e $pca_dir);

            my $pca            = "$pca_dir/pca.diff_compound.pdf";
            my $pca_biplot     = "$pca_dir/pca_biplot.diff_compound.pdf";
            my $pca_coordinate = "$pca_dir/pca_coordinate.diff_compound.csv";
            run_commond("[$compare_name][$sub_compare_name] pca diff_compound", "$SOFT_R_SCRIPT  $SOFT_PCA -i $normalized_data --sample_group_file $group_file --id_file $diff_compound_id --pcaplot $pca --biplot $pca_biplot --pca_coordinate $pca_coordinate", $log_compare_run );
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] pca diff_compound", 'txt', $pca_coordinate, 'pdf', $pca, $pca_biplot);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] pca diff_compound", '[Skip] sample count is not over 2 or diff compound_count is not over 2'); }

        ###########
        # # 基于差异代谢物的绘图：plsda
        ###########
        if($group_count > 1 and $sample_count > 2 and $diff_compound_count > 2)
        {   
            my $plsda_dir = "$compare_dir/plsda/$sub_compare_name";
            system "mkdir -p $plsda_dir" if(not -e $plsda_dir);

            my $plsda            = "$plsda_dir/plsda.diff_compound.pdf";
            my $plsda_vip        = "$plsda_dir/plsda_vip.diff_compound.pdf";
            my $plsda_coordinate = "$plsda_dir/plsda_coordinate.diff_compound.txt";
            my $plsda_vip_txt    = "$plsda_dir/plsda_vip.diff_compound.txt";
            run_commond("[$compare_name][$sub_compare_name] plsda diff_compound", "$SOFT_R_SCRIPT  $SOFT_PLSDA -i $normalized_data --sample_group_file $group_file --id_file $diff_compound_id --show_confidence --plsda $plsda --vip $plsda_vip --plsda_coordinate $plsda_coordinate --plsda_vip $plsda_vip_txt", $log_compare_run );
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] plsda diff_compound", 'txt', $plsda_coordinate, $plsda_vip_txt, 'pdf', $plsda, $plsda_vip);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] plsda diff_compound", '[Skip] group count is not over 1 and sample count is not over 2 or diff compound_count is not over 2'); }

        ###########
        # # 基于差异代谢物的绘图: heatmap
        ###########
        if($sample_count > 1 and $diff_compound_count > 1)
        {   
            my $heatmap_dir = "$compare_dir/heatmap/$sub_compare_name";
            system "mkdir -p $heatmap_dir" if(not -e $heatmap_dir);

            my $heatmap = "$heatmap_dir/heatmap.diff_compound.pdf";
            run_commond("[$compare_name][$sub_compare_name] heatmap diff_compound", "$SOFT_R_SCRIPT  $SOFT_HEATMAP -i $normalized_data --sample_group_file $group_file --id_file $diff_compound_id  -o $heatmap --scale row --cluster_row --cluster_col --show_row --show_col", $log_compare_run );
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] heatmap diff_compound", 'pdf', $heatmap);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] heatmap diff_compound", '[Skip] sample count is not over 1 or diff compound_count is not over 1'); }

        ###########
        # # 基于差异代谢物的绘图: dendrogram
        ###########
        if($sample_count > 1 and $diff_compound_count > 1)
        {
            my $dendrogram_dir = "$compare_dir/dendrogram/$sub_compare_name";
            system "mkdir -p $dendrogram_dir" if(not -e $dendrogram_dir);

            my @dendrogram_results;
            # 距离计算方法
            foreach my $distance_method(qw(euclidean pearson spearman))  
            {   # 聚类计算方法
                foreach my $cluster_method(qw(ward.D single complete average))  
                {   
                    my $dendrogram = "$dendrogram_dir/dendrogram_$distance_method\_$cluster_method.diff_compound.pdf";
                    push @dendrogram_results, $dendrogram;
                    run_commond("[$compare_name][$sub_compare_name] dendrogram diff_compound distance=$distance_method cluster=$cluster_method", "$SOFT_R_SCRIPT  $SOFT_DENDROGRAM -i $normalized_data --sample_group_file $group_file --id_file $diff_compound_id --clusterplot $dendrogram --smplDist $distance_method --clstDist $cluster_method", $log_compare_run );
                }
            }
            check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] dendrogram diff_compound", 'pdf', @dendrogram_results);
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] dendrogram diff_compound", '[Skip] sample count is not over 1 or compound_count is not over 1'); }

        ###########
        # # kegg
        ###########
        if(is_file_ok($diff_compound_id) == 1)
        {   
            my $kegg_dir = "$compare_dir/kegg/$sub_compare_name";
            system "mkdir -p $kegg_dir" if(not -e $kegg_dir);

            # 根据ID，提取差异的通路代谢物ID
            if(is_file_ok($diff_compound_kegg) == 1)
            {
                run_commond("[$compare_name][$sub_compare_name] kegg", "$SOFT_R_SCRIPT  $SOFT_CLUSTERPROFILER -c $diff_compound_kegg --pathway_db $kegg_db --output_dir $kegg_dir/ --download_pathway ", $log_compare_run );
                check_result_file($log_compare_step, 'add', "[$compare_name][$sub_compare_name] kegg", 'txt', "$kegg_dir/kegg_enrichment.xls");
            }else{
                system("rm -r $kegg_dir");
                add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] kegg", '[Skip] no qualify significant compount');
            }
        }else{ add_process_log($log_compare_step, 'add', "[$compare_name][$sub_compare_name] kegg", '[Skip] no significant compount'); }
    }

    # 所有的步骤log追加
    system("cat $log_compare_step >> $log_step");
    system("cat $log_compare_run  >> $log_run");
    system("echo Finish > $compare_dir/finish.txt");

	$ForkManager->finish;
}
$ForkManager->wait_all_children;
my $error = `grep '\\[Error\\]' $log_step`;
if($error =~/\w/)
{
    print "\n[ERROR] 运行错误, $error\n";
    print "请核查文件：\n    $log_step \n    $log_run\n";
}else{
    print "\n[OK] 运行完成\n";
}


###################################################################### 子程序
sub run_commond{   
    my $tag          = shift @_;
    my $commond      = shift @_;
    my $run_log_file = shift @_;
    my $write_model  = "add";
       $write_model  = shift @_ if(exists $_[0]);

    system("echo start $DATA_TIME                           > $run_log_file") if($write_model eq "new");
    system("echo                                           >> $run_log_file");
    system("echo '##########' $tag $DATA_TIME '##########' >> $run_log_file");  # 把tag加入log文件
    system("echo $commond                                  >> $run_log_file");  # 把命令加入log文件
    system("$commond                                       >> $run_log_file 2>&1 ");  # 执行命令
}

sub check_result_file{
    my $process_log_file = shift @_;
    my $write_model      = shift @_;
    my $tag              = shift @_;
    my @results          = @_;

    my $file_type = "";
    my $is_ok     = 1;  # 默认全部正常
    my @error;
    foreach my $result(@results)
    {   
        # 文件类型确认
        if($result eq 'txt' or $result eq 'pdf')
        {
            $file_type = $result;
            next;
        }

        # 检查txt文件
        if($file_type eq 'txt' and is_file_ok($result) == 0)
        {
            $is_ok = 0;
            push @error, $result;
        }

        # 检查pdf文件
        if($file_type eq 'pdf' and is_pdf_ok($result) == 0)  # pdf文件现在先用file模式检测，后面增加pdf专用检测方案
        {
            $is_ok = 0;
            push @error, $result;
        }
    }

    my $condition = ($is_ok == 1) ? '[Good]' : "[Error]";
    
    system "echo $DATA_TIME $tag $condition @error >> $process_log_file" if($write_model eq 'add');  # 追加模式
    system "echo $DATA_TIME $tag $condition @error >  $process_log_file" if($write_model eq 'new');  # 新建模式
}

sub add_process_log{
    my $process_log_file = shift @_;
    my $write_model      = shift @_;
    my $tag              = shift @_;
    my $info             = shift @_;
    system "echo $DATA_TIME $tag $info  >> $process_log_file" if($write_model eq 'add');  # 追加模式
    system "echo $DATA_TIME $tag $info  >  $process_log_file" if($write_model eq 'new');  # 新建模式
}


sub parse_diff_compound{
    my $metabolite_name_map_format = shift @_;
    my $diff_test_file             = shift @_;
    my $plsda_vip                  = shift @_;
    my $diff_compound_id           = shift @_;
    my $diff_compound_kegg         = shift @_;
    my $diff_summary               = shift @_;
    my $pvalue_cutoff              = shift @_;
    my $vip_cutoff                 = shift @_; 

    my %hashMap  = read_matrix($metabolite_name_map_format, 0);
    my %hashTest = read_matrix($diff_test_file, 0);
    my %hashVip  = read_matrix($plsda_vip, 0);

    # map文件要提取的表头
    my @heads_map = @{$hashMap{'HEAD'}};
    my @heads_map_need = @heads_map[2..$#heads_map];

    # vip 文件要提取的表头
    my @heads_vip = @{$hashVip{'HEAD'}};
    my @heads_vip_need = @heads_vip[2..$#heads_vip];

    # test 表头
    my @heads_test = @{$hashTest{'HEAD'}};
    
    # 最终表头
    my @heads = (@heads_test[0..1], @heads_map_need, @heads_vip_need, @heads_test[2..$#heads_test]);

    open DIFF_ID, ">$diff_compound_id";
    open DIFF_KEGG, ">$diff_compound_kegg";
    open DIFF_SUMMARY, ">$diff_summary";
    print DIFF_SUMMARY (join "\t", @heads) . "\n";
    foreach my $id( sort {$hashMap{'DATA'}{$a}{'sort_value'}<=>$hashMap{'DATA'}{$b}{'sort_value'}}  keys %{$hashMap{'DATA'}})
    {
        my $pvalue = $hashTest{'DATA'}{$id}{'P_value'};
        my $vip    = $hashVip{'DATA'}{$id}{'vip scores'};
        my $kegg   = $hashMap{'DATA'}{$id}{'KEGG'};
        print DIFF_ID "$id\n"     if($pvalue < $pvalue_cutoff and $vip > $vip_cutoff);
        print DIFF_KEGG "$kegg\n" if($pvalue < $pvalue_cutoff and $vip > $vip_cutoff and $kegg =~ /^C\d+/);

        
        my @values_map = map{ $hashMap{'DATA'}{$id}{$_} } @heads_map_need;
        my @values_vip = map{ $hashVip{'DATA'}{$id}{$_} } @heads_vip_need;
        my @values_test_s = map{ $hashTest{'DATA'}{$id}{$_} } @heads_test[0..1];
        my @values_test_e = map{ $hashTest{'DATA'}{$id}{$_} } @heads_test[2..$#heads_test];

        my @values = (@values_test_s, @values_map, @values_vip, @values_test_e);
        print DIFF_SUMMARY (join "\t", @values) . "\n";
    }
    close DIFF_ID;
    close DIFF_KEGG;
    close DIFF_SUMMARY;

}

sub parse_group{
    my $group_names = shift @_;
    my $hashGSample = shift @_;
    my $thisdir     = shift @_;

    my %hashGroup;
    # (1) 所有样本的分组
    my $group_file = "$thisdir/sample_group.all_samples.txt";
    $hashGroup{'all_samples'}{'File'} = $group_file;  # 分组文件
	open GROUP_ALL, ">$group_file";
    my @combines;
	foreach my $index(0..$#$group_names)
    {   
        $hashGroup{'all_samples'}{'GroupCount'}++;  # 分组数量
        my $group_name = $$group_names[$index];
		foreach my $sample(split /,/, $hashGSample->{$group_name})
        {
			print GROUP_ALL "$sample\t$group_name\n";
            $hashGroup{'all_samples'}{'SampleCount'}++;
		}
        if($index < $#$group_names)
        {
            foreach my $index2(($index + 1)..$#$group_names)
            {
                push @combines, [$group_name, $$group_names[$index2]];
            }
        }
	}
	close GROUP_ALL;

    # (2) 两两对比的
    @combines = () if($hashGroup{'all_samples'}{'GroupCount'} < 3 );
    foreach my $combine(@combines)
    {   
        my $group_name_1 = $combine->[0];
        my $group_name_2 = $combine->[1];
        my $sub_name = "$group_name_1\_VS_$group_name_2";
        my $group_file = "$thisdir/sample_group.$sub_name.txt";

        $hashGroup{$sub_name}{'File'}       = $group_file;
        $hashGroup{$sub_name}{'GroupCount'} = 2;
        open GROUP, ">$group_file";
		foreach my $sample(split /,/, $hashGSample->{$group_name_1})
        {
			print GROUP "$sample\t$group_name_1\n";
            $hashGroup{$sub_name}{'SampleCount'}++;
		}
		foreach my $sample(split /,/, $hashGSample->{$group_name_2})
        {
			print GROUP "$sample\t$group_name_2\n";
            $hashGroup{$sub_name}{'SampleCount'}++;
		}
        close GROUP;

    }

    return %hashGroup;
}

sub format_name_map{
    my $metabolite_normalized = shift @_;
    my $name_map_raw          = shift @_;
    my $output                = shift @_;
    
    my %hashNormalized = read_matrix($metabolite_normalized, 0);
    my %hashMap        = read_matrix($name_map_raw, 0);

    my @heads = ('ID', @{$hashMap{'HEAD'}});
    open OUT, ">$output";
    print OUT (join "\t", @heads) . "\n";

    foreach my $id( sort {$hashNormalized{'DATA'}{$a}{'sort_value'}<=>$hashNormalized{'DATA'}{$b}{'sort_value'}} keys %{$hashNormalized{'DATA'}})
    {   
        my @values = ($id);
        my $query  = $hashNormalized{'DATA'}{$id}{'Name'};  # 代谢物原始名称
        foreach my $head(@{$hashMap{'HEAD'}})
        {
            my $value = exists $hashMap{'DATA'}{$query}{$head} ? $hashMap{'DATA'}{$query}{$head} : "";
            push @values, $value;
        }
        print OUT (join "\t", @values) . "\n";
    }
    close OUT;

}
sub copy_normalized_data{
    my $input  = shift @_;
    my $output = shift @_;

    my %hashData = read_matrix($input, 0);
    my @heads_input = @{$hashData{'HEAD'}};

    my @samples = @heads_input[2..$#heads_input];
    my @heads = ('ID', 'Name', @samples);

    open OUT, ">$output";
    print OUT (join "\t", @heads) . "\n";

    foreach my $id( sort {$hashData{'DATA'}{$a}{'sort_value'}<=>$hashData{'DATA'}{$b}{'sort_value'}} keys %{$hashData{'DATA'}})
    {   
        my $query  = $hashData{'DATA'}{$id}{$hashData{'HEAD'}->[1]};  # 代谢物原始名称
        my @sample_values = map{ $hashData{'DATA'}{$id}{$_} } @samples;
        my @values = ($id, $query, @sample_values);
        print OUT (join "\t", @values) . "\n";
    }
    close OUT;
}

sub read_config{
	my $file = shift;
	my %hash;
	my $sort = 0;
	open FILE, $file;
	while(my $line = <FILE>){
		$line =~ s/#.+//g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		my ($a, $b) = split /\s*=\s*/, $line, 2;
		next if(not defined $b);
		next if($a ne "Compare");
		$b =~ s/\s//g;
		$hash{$b} = $sort;
		$sort++;
	}
	close FILE;
	return %hash;
}

sub read_sample_group{
	my $file = shift;
	open FILE, $file;
	my $title = <FILE>;
	$title =~ s/[\r\n]//g;
	my @split_title = split /\t/, $title;
	my %hash;
	my $sort = 0;
	while(my $line = <FILE>){
		$line =~ s/[\r\n]//g;
		my @split_line = split /\t/, $line;
		for my $i(2..@split_line-1){
			$hash{$split_title[$i]}{$split_line[$i]} .= "$split_line[1],";
		}
	}
	close FILE;
	return %hash;
}

sub find_same_sample{
	my $hash = shift;
    my %hashSample;
    foreach my $group_name(keys %$hash)
    {   
        map{ $hashSample{$_}++ } split /,/, $hash->{$group_name};
    }
    my @dups = grep{ $hashSample{$_} > 1 } keys %hashSample;
    my $dup_count = scalar(@dups);
    return 1 if($dup_count > 0);
    return 0;
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
# 检验pdf文件是否为空
sub is_pdf_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file < 4400);
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

# 读取带有表头的矩阵
sub read_matrix{
    my $file     = shift @_;
    my @key_cols = @_;

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
        my @key_values = map{ $datas[$_] } @key_cols;
        my $key_value  = join "\t", @key_values;
        foreach my $col(0..$#heads)
        {
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $hashMatrix{'DATA'}{$key_value}{$heads[$col]} = $value;
        }
        $hashMatrix{'DATA'}{$key_value}{'sort_value'} = $.;
    }
    close FILE;
    return %hashMatrix;
}

sub check_sample{
    my $metabolite_file = shift @_;
    my $sample_group = shift @_;
    
    open META, $metabolite_file;
    my $mline1 = <META>;
       $mline1=~s/[\r\n]//g;
    my ($id, $name, @samples) = split /\t/, $mline1;
    my %hashSample = map{ ($_, 1) } @samples;
    close META;

    my @losts;
    open GROUP, $sample_group;
    <GROUP>;
    while(<GROUP>)
    {
        $_=~s/[\r\n]//g;
        my ($id, $sample) = split /\t/, $_;
        next if(exists $hashSample{$sample});
        push @losts, $sample;
    }
    die "[Error] 分组文件的部分样本在丰度表里缺失，请仔细核对！ 丢失样本： @losts\n" if(scalar(@losts) > 0);
}