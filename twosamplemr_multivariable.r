#!/home/genesky/software/r/3.6.1/bin/Rscript

library(docopt)

"Usage: twosamplemr_singlevariable.r -e <gwasid_list> -o <gwasid> -r <dir> -p <string> [--p1 <pvalue_cutoff>  --rlib <dir> --no_check_id]

Options:
    -e, --exposure <gwasid_list>   暴露因素id， gwasid_list 提供2个及以上暴露因素，用于检测多个暴露因素的综合效应。 多个id用逗号分隔。
    -o, --outcome <gwasid>         结局因素id， gwasid, 例如： ebi-a-GCST005647 与Amyotrophic lateral sclerosis肌萎缩性脊髓侧索硬化症 相关的gwas数据库 
                                   具体数据库请在 https://gwas.mrcieu.ac.uk/ 搜索
    -r, --result_dir <dir>         结果输出目录, 例如: ./
    -p, --prefix <string>          输出目录下，所有文件添加前缀, 例如： e_o

    --p1 <pvalue_cutoff>           exposure突变位点提取时的pvalue限制 [default: 5e-08]
                                   比较严格，但是有时候也会导致没有任何SNP位点满足条件，此时请适当降低阈值，但也不要太低，否则提取的SNP过多
    --no_check_id                  不检查exposure outcome id编号是否正确，因为要联网下载数据库，有点慢。可以取消掉。
    --rlib <dir>                   R包路径 [default: /home/genesky/software/r/3.6.1/lib64/R/library]" -> doc
# 以SNP为中介，检测暴露因素 与 结局因素 是否存在因果关系
opts        <- docopt(doc, version='甘斌，twosamplemr  多因素因素分析， 参考网站： https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html \n')
exposure    <- opts$exposure
outcome     <- opts$outcome
result_dir  <- opts$result_dir
prefix      <- opts$prefix
p1          <- as.numeric(opts$p1)
no_check_id <- opts$no_check_id
rlib        <- opts$rlib

exposure = unlist(strsplit(exposure, ','))

if(!is.null(rlib)) .libPaths(rlib)

# exposure = 'ukb-b-5443,ukb-d-III_BLOOD_IMMUN'
# outcome = 'ebi-a-GCST005647'

message("loading TwoSampleMR")
library(TwoSampleMR )
set.seed(91)

### 获取gwas数据库列表
message('[process 0] check input id')
if(!no_check_id)
{
    ao <- available_outcomes()

    if(sum(!exposure %in% ao$id)){
        message("exposure is 错误. 在gwas数据库中库中不存在 ")
    }else{
        message("exposure id 正常")
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
exposure_dat <- mv_extract_exposures(exposure, pval_threshold = p1)  # 也有可能会出现重复SNP, 不过，这个函数已经考虑过这个问题了，不需要再额外处理
if(is.null(exposure_dat)){
    message("暴露因素没有提取到任何位点，请适当降低 p1 过滤阈值， 当前使用的阈值为： ", p1)
    q()
}else{
    message('SNP count in exposure data： ', nrow(exposure_dat))
}

# 结局数据
outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcome)
outcome_dat <- outcome_dat[!duplicated(outcome_dat$SNP), ]  # 去掉重复，发现在部分数据库里会找到重复的数据
 
### （2）数据整合清洗
message('\n[process 2] harmonise_data')
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

last_snp_count = nrow(mvdat$exposure_beta)  # 最终SNP数量
write.table(mvdat, file=paste0(result_dir, "/", prefix, ".harmonise_dat.txt"), row.names =F, quote=F, sep='\t')

message("SNP left: ", last_snp_count)
if(last_snp_count == 0){
    message("harmonise_data 后数据量为0， 请适当降低 pvalue过滤阈值，增加SNP数量。 当前使用的阈值为： ", p1)
    q()
}

### （3）MR分析
message('\n[process 3] Perform MR')
# 综合考虑所有SNP
res <- mv_multiple(mvdat)
write.table(res, file=paste0(result_dir, "/", prefix, ".mv_multiple.txt"), row.names =F, quote=F, sep='\t')

