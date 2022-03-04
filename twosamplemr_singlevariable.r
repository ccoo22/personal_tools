#!/home/genesky/software/r/3.6.1/bin/Rscript

library(docopt)

"Usage: twosamplemr_singlevariable.r -e <gwasid> -o <gwasid> -r <dir> -p <string> [--p1 <pvalue_cutoff> --r2 <r2_cutoff>  --sample_size_exposure <int>  --f_filter <numeric> --rlib <dir> --no_check_id --html ]

Options:
    -e, --exposure <gwasid>       暴露因素id， gwasid, 例如： ukb-b-5443 与免疫相关的gwas数据库 
    -o, --outcome <gwasid>        结局因素id， gwasid, 例如： ebi-a-GCST005647 与Amyotrophic lateral sclerosis肌萎缩性脊髓侧索硬化症 相关的gwas数据库 
                                  具体数据库请在 https://gwas.mrcieu.ac.uk/ 搜索
    -r, --result_dir <dir>        结果输出目录, 例如: ./
    -p, --prefix <string>         输出目录下，所有文件添加前缀, 例如： e_o

    --p1 <pvalue_cutoff>          exposure突变位点提取时的pvalue限制 [default: 5e-08]
                                  比较严格，但是有时候也会导致没有任何SNP位点满足条件，此时请适当降低阈值，但也不要太低，否则提取的SNP过多
    --r2 <r2_cutoff>              exposure突变位点提取时的r2限制 [default: 0.001]
    --sample_size_exposure <int>  暴露因素中的样本个数 [default: 0]
                                  计算F统计量时，要用到它。0 表示使用官方下载的表格中的 samplesize进行计算
                                  但，有时候，官方表格的 samplesize.exposure 是 NA，此时，可以自己声明这个数值（去https://gwas.mrcieu.ac.uk/搜索暴露因素的id，能找到Sample size信息）
    --f_filter <numeric>          根据F统计值过滤。 0 表示不过滤 [default: 0]
                                  通常认为，F统计值大于10的点比较可信。
    --no_check_id                 不检查exposure outcome id编号是否正确，因为要联网下载数据库，有点慢。可以取消掉。
    --html                        制作html版本。注意：如果选择了这个，所有的分析内容会被做两遍，第二遍是为了整合出html。
    --rlib <dir>                  R包路径 [default: /home/genesky/software/r/3.6.1/lib64/R/library]

                                  注意：需要联网获取数据，请保持网络畅通，偶尔官方网站也会宕机，请耐心等待，重试。
    
    " -> doc
# 以SNP为中介，检测暴露因素 与 结局因素 是否存在因果关系
opts        <- docopt(doc, version='甘斌，twosamplemr  单因素分析， 参考网站： https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html \n')
exposure    <- opts$exposure
outcome     <- opts$outcome
result_dir  <- opts$result_dir
prefix      <- opts$prefix
p1          <- as.numeric(opts$p1)
r2          <- as.numeric(opts$r2)
f_filter    <- as.numeric(opts$f_filter)
sample_size_exposure <- as.integer(opts$sample_size_exposure)
no_check_id <- opts$no_check_id
html        <- opts$html
rlib        <- opts$rlib
 
if(!is.null(rlib)) .libPaths(rlib)

# exposure = 'ukb-b-5443'
# outcome = 'ebi-a-GCST005647'

message("loading TwoSampleMR")
library(TwoSampleMR)
library(MRPRESSO)
library(ggplot2 )
set.seed(91)

### 获取gwas数据库列表
message('[process 0] check input id')
if(!no_check_id)
{
    ao <- available_outcomes()

    if(exposure %in% ao$id){
        message("exposure id 正常")
    }else{
        message("exposure is 错误. 在gwas数据库中库中不存在 ")
        q()
    }
    if(outcome %in% ao$id){
        message("outcome id 正常")
    }else{
        message("outcome is 错误. 在gwas数据库中库中不存在 ")
        q()
    }   
}else{
    message('no check input id')
}


### （1）读入数据
message('\n[process 1] extract data ')
# 暴露数据
exposure_dat <- extract_instruments(outcomes=exposure, p1 = p1, r2 = r2)  # 也有可能会出现重复SNP, 不过，这个函数已经考虑过这个问题了，不需要再额外处理
if(is.null(exposure_dat)){
    message("暴露因素没有提取到任何位点，请适当降低 p1/r2 过滤阈值， 当前使用的阈值为： p1 = ", p1, "; r2 = ", r2)
    q()
}else{
    message('暴露因素找到的snp位点数目为： ', nrow(exposure_dat))
    message('    snp id :', paste0(exposure_dat$SNP, collapse=','))
}

# 结局数据
outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcome )
outcome_dat <- outcome_dat[!duplicated(outcome_dat$SNP), ]  # 去掉重复，发现在部分数据库里会找到重复的数据
if(is.null(outcome_dat)){
    message("结局因素没有提取到任何位点，请适当降低 p1/r2 过滤阈值， 当前使用的阈值为： p1 = ", p1, "; r2 = ", r2)
    q()
}

### （2）数据整合清洗
message('\n[process 2] harmonise_data')
harmonise_dat <- harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat
)
last_snp_count = sum(harmonise_dat$mr_keep)  # 最终SNP数量

if(sample_size_exposure > 0)
{
    harmonise_dat$samplesize.exposure = sample_size_exposure
}

# 补充暴露因素F统计量
harmonise_dat$r2.exposure = harmonise_dat$eaf.exposure * (1 - harmonise_dat$eaf.exposure) * harmonise_dat$beta.exposure * harmonise_dat$beta.exposure
harmonise_dat$fstatistic.exposure = (harmonise_dat$samplesize.exposure - 2) * harmonise_dat$r2.exposure / (1 - harmonise_dat$r2.exposure)

if(f_filter > 0)
{
    message("根据F值过滤： ", f_filter)
    harmonise_dat = harmonise_dat[harmonise_dat$fstatistic.exposure >= f_filter, ]
}

write.table(harmonise_dat[harmonise_dat$mr_keep, ], file=paste0(result_dir, "/", prefix, ".harmonise_dat.txt"), row.names =F, quote=F, sep='\t')

message("SNP left: ", last_snp_count)
if(last_snp_count == 0){
    message("harmonise_data 后数据量为0， 请适当降低 pvalue过滤阈值，增加SNP数量。 当前使用的阈值为： ", p1)
    q()
}

## 稳健调整轮廓评分 Robust adjusted profile score
## 如果选择的SNPs是弱的工具变量，则mr函数会无法运行。此时，可以用 mr.raps 方法运行。它是MR的一个备选。
# TwoSampleMR函数mr_raps更新不及时，无法运行，自己根据需要重新写一遍
parameters = default_parameters()
mr_raps_result = mr.raps::mr.raps(b_exp = harmonise_dat$beta.exposure,
        b_out=harmonise_dat$beta.outcome,
        se_exp=harmonise_dat$se.exposure,
        se_out=harmonise_dat$se.outcome,
        over.dispersion=parameters$over.dispersion,
        loss.function=parameters$loss.function,
        diagnosis=FALSE,
        )
mr_raps_result = data.frame(mr_raps_result, nsnp=nrow(harmonise_dat))
mr_raps_result$or = exp(mr_raps_result$beta.hat)
mr_raps_result$L95 = exp(mr_raps_result$beta.hat - 1.96 * mr_raps_result$beta.se)
mr_raps_result$U95 = exp(mr_raps_result$beta.hat + 1.96 * mr_raps_result$beta.se)
write.table(mr_raps_result, file=paste0(result_dir, "/", prefix, ".mr_RAPS.txt"), row.names =F, quote=F, sep='\t', col.names=T)

### （3）MR分析
message('\n[process 3] Perform MR')
# 综合考虑所有SNP
message('    mr')
res <- mr(harmonise_dat)
res$or = exp(res$b)
res$L95 = exp(res$b - 1.96 * res$se)
res$U95 = exp(res$b + 1.96 * res$se)
write.table(res, file=paste0(result_dir, "/", prefix, ".mr.txt"), row.names =F, quote=F, sep='\t')

# 每个SNP单独考虑
message('    mr_singlesnp')
res_single <- mr_singlesnp(harmonise_dat)
write.table(res_single, file=paste0(result_dir, "/", prefix, ".mr_singlesnp.txt"), row.names =F, quote=F, sep='\t')

# 留一法，分别检测
message('    mr_leaveoneout')
res_loo <- mr_leaveoneout(harmonise_dat)
write.table(res_loo, file=paste0(result_dir, "/", prefix, ".mr_leaveoneout.txt"), row.names =F, quote=F, sep='\t')

# 因果方向是否正确， TRUE/FALSE
message('    directionality_test')
res_direction = tryCatch({
        directionality_test(harmonise_dat)
    }, warning = function(w){
        message('warning')
    }, error = function(e){
        message("directionality_test ERROR")
    })

if(!is.null(res_direction))
{
    write.table(res_direction, file=paste0(result_dir, "/", prefix, ".directionality_test.txt"), row.names =F, quote=F, sep='\t')
}

### （4）sensitivity analysis 分析
message('\n[process 4] sensitivity analysis')

# 工具变量异质性检测，如果p< 0.05, 说明这些SNP之间存在很强的异质性，建议剔除某些outcome的p值非常小的SNP.
message('    mr_heterogeneity')
res_heterogeneity <- mr_heterogeneity(harmonise_dat)
write.table(res_heterogeneity, file=paste0(result_dir, "/", prefix, ".mr_heterogeneity.txt"), row.names =F, quote=F, sep='\t')

# 水平多效性检测：p< 0.05, 说明存在水平多效性，表明因果关系不成立。
message('    mr_pleiotropy_test')
res_intercept <- mr_pleiotropy_test(harmonise_dat)
write.table(res_intercept, file=paste0(result_dir, "/", prefix, ".mr_pleiotropy_test.txt"), row.names =F, quote=F, sep='\t')

# MR-PRESSO 另一种检测水平多效性的方法。P值小于0.05，说明暴露因素和结局变量存在水平多效性
message('    mr_presso')
res_presso <- mr_presso(data = harmonise_dat[harmonise_dat$mr_keep, ], BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 1000,  SignifThreshold = 0.05)
presso_pvalue <- res_presso[['MR-PRESSO results']][['Global Test']]$Pvalue
presso_rssobs <- res_presso[['MR-PRESSO results']][['Global Test']]$RSSobs
presso_result = data.frame(pvalue=presso_pvalue, RSSobs=presso_rssobs)
write.table(presso_result, file=paste0(result_dir, "/", prefix, ".mr_presso.MR-PRESSON-results.txt"), row.names =F, quote=F, sep='\t', col.names=T)
write.table(res_presso[['Main MR results']], file=paste0(result_dir, "/", prefix, ".mr_presso.Main-MR-results.txt"), row.names =F, quote=F, sep='\t', col.names=T)


### （5） plot
message("\n[process 5] plot ")
# 散点图，横纵坐标是每一个SNP位点在暴露、结局因素中的effect值，并绘制每一个差异分析方法的拟合曲线，其斜率对应b值.是MR拟合结果
# 可以直观地看出来，随着exposure的增长，outcome的变化
p_scatter <- mr_scatter_plot(res, harmonise_dat)
ggplot2::ggsave(p_scatter[[1]], file=paste0(result_dir, "/", prefix, ".scatter.pdf"), width=7, height=7)

# 对每个SNP单独分析的结果effect size 绘图，方便查看每一个SNP的贡献.每一行代表一个SNP/差异方法， 横坐标是b值。 Wald ratio 算法。
# 检查每个位点对outcome的影响
p_forest <- mr_forest_plot(res_single)
ggplot2::ggsave(p_forest[[1]], file=paste0(result_dir, "/", prefix, ".single.forest.pdf"), width=7, height = 7 + last_snp_count * 0.01)

# funnel 图，每一个点代表一个snp, 横坐标是SNP的b值，纵坐标是1/se值， 垂直线是差异算法的b值
# 如果有点（SNP）偏离总体较远，说明它有异质性，建议去除这种位点
p_funnel <- mr_funnel_plot(res_single)
ggplot2::ggsave(p_funnel[[1]], file=paste0(result_dir, "/", prefix, ".single.funnel.pdf"), width=7, height=7)

# 留一法结果绘图，展示效果与但SNP一样.每一行代表去掉的SNP， 横坐标是b值
# 检查模型的稳健型
p_leaveoneout <- mr_leaveoneout_plot(res_loo)
ggplot2::ggsave(p_leaveoneout[[1]], file=paste0(result_dir, "/", prefix, ".leaveoneout.pdf"), width=7, height= 7 + last_snp_count * 0.01)

if(html)
{
    message('\n[process 6] make html. all analysis will be done again ')

    html_raw <- mr_report(dat=harmonise_dat, output_path = result_dir)
    file.rename(paste0(result_dir, "/", html_raw) , paste0(result_dir, "/", prefix, ".html"))  # 重命名
    file.remove(paste0(result_dir, "/mr_report.md"))  # 删除不要
    system(paste0(result_dir, "/figure"))  # 删除不要
}
                                                                                                                                                                                                                     