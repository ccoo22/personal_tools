#!/usr/bin/env Rscript

library(docopt)
library(parallel)
"Usage: correlation_analysis.r [options] METH_FILE EXPR_FILE OUTPUT_FILE
Options:
   -c the cor.test method     [default: pearson]
   -p parallel                [default: 10]
Arguments:
  METH_FILE      the methylation file name
  EXPR_FILE      the expression file name
  OUTPUT_FILE    the output file" -> doc


opts              <- docopt(doc)
meth_file         <- opts$METH_FILE
expr_file         <- opts$EXPR_FILE
output_file       <- opts$OUTPUT_FILE
cor_method        <- opts$c
parallel          <- opts$p

# 相关计算函数
correlation <- function(meth_name){
    meth_data1 = meth_data[, meth_name] # 甲基化数据

    res <- matrix(0, nrow=length(expr_names), ncol=5) # 创建空矩阵，用于存放数据
    for(col in 1:length(expr_names))
    { 
      expr_name  = expr_names[col]
      expr_data1 = expr_data[, expr_name] # 表达量数据
      
      data <- cbind(meth_data1, expr_data1)
      data <- data.frame(data[complete.cases(data),]) # 去掉缺失值
      
      nmiss <- nrow(data)
      pvalue = NULL
      estimate = NULL

      # 样本量>=3才能计算
      if(nmiss > 2)
      {
        result <- cor.test(data[, 1], data[, 2], alternative = "two.sided", method = cor_method, conf.level = 0.95)
        pvalue <- result$p.value
        estimate <- result$estimate[[1]]
      }      

      res[col, ] = c(meth_name, expr_name, nmiss, pvalue, estimate)
    }

    res

}

# 读入数据
cat("read data\n");
meth_data <- read.table(meth_file, sep="\t", header=TRUE, row.name=1, check.names=FALSE, comment.char="")
expr_data <- read.table(expr_file, sep="\t", header=TRUE, row.name=1, check.names=FALSE, comment.char="")
meth_data = meth_data * 100 # 甲基化转为百分比形式

meth_names = colnames(meth_data)
expr_names = colnames(expr_data)

if(length(meth_names) < parallel ) parallel = length(meth_names)  # 线程数不要超过实际使用的数量,否则流程可能无法正常结束

# 并行运算
cat("start parallel\n")
cl <- makeForkCluster(nnodes = parallel) # 创建计算组，使用了所有的核心
clusterExport(cl, c('correlation', 'meth_data', 'expr_data', 'expr_names', 'cor_method')) # 通知计算组要用到的变量
res.p <- parLapply(cl, meth_names, function(x){ correlation(x) } ) # 并行运算 ，后台也是出现n个独立的R进程，都是用了100CPU
stopCluster(cl) # 解散计算组

# 数据unlist合并
res = do.call("rbind", res.p)
res = res[order(res[,4], decreasing=F), ] # 按照P值升序排列
colnames(res) = c('MethName', 'ExprName', 'NMISS', paste0(cor_method, '.Pvalue'), paste0(cor_method, '.Estimate'))

# 输出
write.table(res, output_file , sep = "\t", quote = F, row.names = F )
 
