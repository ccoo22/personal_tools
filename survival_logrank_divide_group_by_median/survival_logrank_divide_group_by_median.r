#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: survival_logrank_divide_group_by_median.r -s <file> -i <file> -o <dir> [ --input_format <string> --rlib <dir> ] 
Options:
    -s, --surv <file>               生存信息矩阵，第一列样本名，第二列time, 第三列 status（0 alive, 1 dead），有表头。不允许缺失。如果有缺失，会自动排除包含缺失的样本。
                                    示例：/home/zhangshuang/work/other_project/20B1217C-3/part2/data/suvival_data/Overall_Survival_Time.txt
    -i, --input <file>              样本的特征值矩阵，每一行对应一个特征，每一列对应一个样本，第一列是特征值ID, 第二列开始是样本的特征值。含有表头，tab分隔符。不允许缺失。如果有缺失，会自动排除包含缺失的样本。
                                    特征值必须是数值，可以是样本的表达量数据，也可以是样本的某个指标或打分
    -o, --output_dir <dir>          结果输出目录
    --input_format <string>         input文件格式 col/row。col表示input文件每一列对应一个样本；row表示exp文件每一行对应一个样本 [default: col]
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc


opts                   <- docopt(doc, version = 'Program : 通过特征值中位数得到样本分组的生存分析曲线绘图 \n          张爽 298\n')
surv                   <- opts$surv
input                  <- opts$input
output_dir             <- opts$output_dir
input_format           <- opts$input_format
rlib                   <- opts$rlib


# surv                   <- "/home/zhangshuang/work/research/survival_logrank_divide_group_by_median/example/input/Overall_Survival_Time.txt"
# input                  <- "/home/zhangshuang/work/research/survival_logrank_divide_group_by_median/example/input/input.txt"
# output_dir             <- "/home/zhangshuang/work/research/survival_logrank_divide_group_by_median/example/result"
# input_format           <- "col"
# rlib                   <- "/home/genesky/software/r/3.5.1/lib64/R/library"


.libPaths(rlib)
library(survival)
library(survminer)

##################
# (1) 读入生存数据
##################
message("[process] loading surv data")
data_surv = read.table(surv, header = T, row.names = 1, stringsAsFactors = FALSE, check.names = F, sep = '\t')  # 生存信息
data_surv = data_surv[, c(1,2)]
colnames(data_surv)[1:2] = c('time', 'status')
samples = rownames(data_surv)

##################
# (2) 读入表达量数据
##################
message("[process] loading input data")
data_input = read.table(input, head = T, row.names = 1, check.names = F, sep = "\t", quote = "")
if(input_format == 'col') data_input = t(data_input)

##################
# (3) 检查样本存在与否
##################
message("[process] check input sample ok")

if(sum(!samples %in% rownames(data_input)) > 0)
{
    lost_samples = samples[!samples %in% rownames(data_input)]
    message("[Error] surv 中的样本在 input 文件中缺失:", lost_samples)
    q()
}

##################
# (4) 提取分析的数据
##################
message("[process] build sample data matrix")
data_clean = data_input[samples, , drop = F]   # 仅保留需要的样本
data_clean <- as.matrix(data_clean)

# (5) 特征名称临时替换，防止存在空格、-等特殊字符，导致报错
feature_names_input         = colnames(data_clean)
feature_names_new           = paste0('feature', 1:ncol(data_clean))
colnames(data_clean)        = feature_names_new
names(feature_names_input)  = feature_names_new

##################
# (6) KM 分析
##################
message("[process] KM analysis")

count = 0
pic_count = 0
pvalue_data = matrix(0, nrow = ncol(data_clean), ncol = 3)
colnames(pvalue_data) = c('Feature', 'Sample Count', 'Pvalue')

pdf_file <- paste(output_dir, "/Survival.Logrank.pdf", sep = "")
pdf(pdf_file, width = 7)
for (i in 1:ncol(data_clean))
{
    feature = colnames(data_clean)[i]
    feature_input = feature_names_input[feature]
    message("[process Log-Rank] ", feature_input)
    count = count + 1
	data_analysis = cbind(data_surv, data_clean[, feature], stringsAsFactors = FALSE)
	colnames(data_analysis) = c("time", "status", feature)
	data_analysis <- data_analysis[complete.cases(data_analysis), ]
	data_analysis <- data_analysis[apply(data_analysis, 1, function(x){sum(x == '')}) == 0, ]

	data_analysis$group <- "Low"
	data_analysis$group[data_analysis[, feature] > median(data_analysis[, feature])] = "High"
	data_analysis_surv <- data_analysis[, c("time", "status", "group")]
    sample_count  = nrow(data_analysis_surv)

    # 只有一组数据，无法分析
    if(length(unique(data_analysis_surv[, 3])) == 1 )
    {
        pvalue_data[count, ] = c(feature_input, sample_count, 'NA')
        next
    }

	pic_count = pic_count + 1
	fit_survival_logrank       <- survfit(Surv(time, status) ~ group, data = data_analysis_surv)
	surv_diff_survival_logrank <- survdiff(Surv(time, status) ~ group, data = data_analysis_surv) # 差异分析
	pvalue_survival_logrank    <- pchisq(surv_diff_survival_logrank$chisq, length(surv_diff_survival_logrank$n)-1, lower.tail = FALSE) # 获取P值
	surv_info_survival_logrank <- surv_summary(fit_survival_logrank, data = data_analysis_surv) # 拟合详细结果 

    # 记录pvalue
    pvalue_data[count, ] = c(feature_input, sample_count, pvalue_survival_logrank)

	p <- ggsurvplot(fit_survival_logrank, data = data_analysis_surv,
        fun = NULL,
        pval = TRUE, 
        conf.int = FALSE,
        risk.table = TRUE, # Add risk table
        surv.median.line = "hv", # Specify median survival
        ggtheme = theme_bw(), # Change ggplot2 theme
        add.all = FALSE,
        palette = "npg", #杂志nature的配色
        legend.title = feature_input,
        legend = "top",
        xlab = "Time",
    )
    # 第一幅图要设置newpage = FALSE, 否则会多一个空白页
    result = tryCatch(print(p, newpage = pic_count !=1), error=function(e){cat("ERROR :",conditionMessage(e),"\n")}) # 注部分图可能无法绘制，报错
}
dev.off()

survival_logrank_result <- paste(output_dir, "/Survival.Logrank.pvalue.txt", sep = "")
write.table(pvalue_data, survival_logrank_result, row.names = F, col.names = T, quote = F, sep = '\t')

