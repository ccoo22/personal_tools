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
if(length(ARGS) < 3){
    cat("
Program: Decision Curve Analysis，DCA
Version: v1.0
Contact: 129 甘斌

Usage:   Rscript", SCRIPT, "OUTDIR PREFIX INPUTFile1 INPUTFile2 ...

Options:

         OUTDIR          结果文件输出目录
         PREFIX          输出结果前缀
         INPUT           输入文件,包含表头，第一列样本名，第二列分组信息（只能有1、0）,第三列往后是用于逻辑回归的所有自变量         
    \n");
    q()
}
OUTDIR <- normalizePath(ARGS[1])
PREFIX <- ARGS[2] 
INPUTFile <- ARGS[3]  # 数据路径

###################################################################### 主程序

library(rms)
set.seed(91)

table = read.table(INPUTFile, head = T, row.names = 1, check.names = F)
colnames(table)[1] = 'y'
table_colnames = colnames(table)
genes = table_colnames[2:length(table_colnames)]

ddist <- datadist(table)
options(datadist='ddist')
# logistic回归
formula1 <- as.formula(paste('y ~',paste(genes , collapse = ' + ') ) )
f <- lrm(formula1, data = table, x=TRUE, y=TRUE)
lrm_result = f$coefficients
write.table(lrm_result, paste(OUTDIR, '/', PREFIX, '.nomogram.logistic.txt', sep = ''), quote = FALSE, row.names = TRUE, sep = '\t')
# nomogram
nom <- nomogram(f, 
    fun = plogis,
    fun.at = c(.001, .01, .05, seq(.15, .85, by=.25), .85, .95, .99, .999),
    lp = F,  
    funlabel = "Risk")
# fun.at = c(.001, .01, .05, seq(.1, .9, by=.1), .95, .99, .999) 默认
pdf(paste(OUTDIR, '/', PREFIX, '.nomogram.pdf', sep = ''), family="GB1", width = 10, height = 7)
plot(nom)
dev.off()

# 校正
pdf(paste(OUTDIR, '/', PREFIX, '.nomogram.calibrate.pdf', sep = ''), family="GB1")
cal <- calibrate(f,  B = 1000)
plot(cal)
dev.off()