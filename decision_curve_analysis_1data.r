# 系统 -> 初始化
R_VERSION <- "3.4.0"
R_SYS_VERSION <- R.version.string #返回R版本
if(!grepl(R_VERSION, R_SYS_VERSION)){
    cat(" -- [ERR] 必须使用 R", R_VERSION, "版本运行此脚本,", R_SYS_VERSION, "\n")
    q()
}
.libPaths("/home/ganb/soft/R-3.4.0/lib64/R/library")

# 检测 -> 脚本输入
init <- function(){
    tmp <- commandArgs();
    tmp <- tmp[grep("--file=", tmp)]
    tmp <- substring(tmp, 8, nchar(tmp))
    tmp_pos <- gregexpr("[\\/]", tmp)
    tmp_pos <- tmp_pos[[1]][length(tmp_pos[[1]])]
    return(list(substring(tmp, 1, tmp_pos), substring(tmp, tmp_pos+1, nchar(tmp))))
}
SCRIPT <- init();
SCRIPTDIR <- SCRIPT[[1]]
SCRIPT <- SCRIPT[[2]]

# 检测 -> 脚本输入
ARGS <- commandArgs(T)
if(length(ARGS) != 3){
    cat("
Program: Decision Curve Analysis，DCA
Version: v1.0
Contact: 129 甘斌

Usage:   Rscript", SCRIPT, "INPUTFile OUTDIR PREFIX

Options:

         INPUTFile       输入文件,包含表头，第一列样本名，第二列分组信息（只能有1、0）,第三列往后是用于逻辑回归的所有自变量
         OUTDIR          结果文件输出目录
         PREFIX          输出结果前缀
    \n");
    q()
}
INPUTFile <- normalizePath(ARGS[1]) #返回绝对路径
OUTDIR <- normalizePath(ARGS[2])
PREFIX <- ARGS[3] 

###################################################################### 主程序
library(DecisionCurve)

# (1) 数据读取
Data <- read.table(INPUTFile,
                         header = TRUE,
                         sep = "\t",
                         quote = "",
                         check.names = FALSE,
                         )
colnames(Data)[1] <- 'Sample' 
colnames(Data)[2] <- 'Group'
col_names <- colnames(Data)
independent_variables <- col_names[3:length(col_names)]

not_allow_character <- grep("-", independent_variables)
if(length(not_allow_character) > 0){
    cat('[ERROR] independent variables has unsupport character [-] : ', independent_variables[not_allow_character], "\n")
    q()
}
formula <- as.formula(paste('Group ~', paste(independent_variables , collapse = ' + ') ) )

# (2) 拟合
fit <- decision_curve(formula, 
    data = Data, 
    family = binomial(link ='logit'),
    thresholds= seq(0,1, by = 0.01),
    confidence.intervals = 0.95,
    study.design = 'case-control'
    )

#（3）绘图
pdf_file = paste(OUTDIR, '/', PREFIX, ".pdf", sep = "", collapse = NULL)
pdf(pdf_file, width = 12, height = 9)
plot_decision_curve(fit, 
           curve.names= c(PREFIX),
           cost.benefit.axis =FALSE,
           col = c('red'),
           confidence.intervals =FALSE,
           standardize = FALSE,
           )
dev.off()

cat('plot: ', pdf_file, "\n")

# (4) 数据输出
txt_file = paste(OUTDIR, '/', PREFIX, ".decision_curve.txt", sep = "", collapse = NULL)
output_data = list()

data_all <- subset(fit$derived.data, is.element(model, c("All", "None"))) # all 数据
data_fit <- fit$derived.data # 拟合数据
thresholds_count <- 101

thresholds =   100 * data_all[, 'thresholds'][1:thresholds_count]
all_sNB    =   data_all[, 'NB'][1:thresholds_count]  
fit_sNB    =   data_fit[, 'NB'][1:thresholds_count] 
reduction  = (fit_sNB - all_sNB) * 100 / (thresholds / (100 - thresholds)) # 公式参考地址 https://www.mskcc.org/sites/default/files/node/4511/documents/v3-worked-example-of-decision-curve-analysis-using-r.pdf

output_data[[paste(PREFIX, 'thresholds', sep = '.')]] = thresholds
output_data[[paste(PREFIX, 'ALL', sep = '.')]] = all_sNB
output_data[[paste(PREFIX, 'Fit', sep = '.')]] = fit_sNB
output_data[[paste(PREFIX, 'Reduction', sep = '.')]] = reduction


output_data = as.data.frame(output_data, check.names = FALSE)

write.table(output_data, 
    txt_file, 
    quote = FALSE, 
    sep = "\t", 
    row.names = FALSE,
    )
cat('write: ', txt_file, "\n")


# reduction绘图
#颜色
mycol <- c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
mycol <- colors()[mycol]
color_list = c('red', 'blue', 'green', mycol)

pdf_file = paste(OUTDIR, '/', PREFIX, ".reduction.pdf", sep = "", collapse = NULL)
pdf(pdf_file, width = 12, height = 9)
reduction <- pmax(0, reduction) 
plot(thresholds, reduction, col = color_list[1], type="l", lwd=2, xlab="Threshold probability (%)", ylab="Net reduction in interventions per 100 patients") 
lines(thresholds, reduction, col = color_list[1], type="l", lwd=2, xlab="Threshold probability (%)", ylab="Net reduction in interventions per 100 patients") 
dev.off()
