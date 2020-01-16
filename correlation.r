# 检测 -> 脚本输入
ARGS <- commandArgs(T)
if(length(ARGS) < 2){
  cat("
Program: correlation 
Version: v1.0
Contact: 129 甘斌
      
Usage:   Rscript stringdb.r ARGS1 ARGS2 ARGS3
      
Options:
      
      INFILE            输入矩阵，3列，第一列注释，第二列样本1，第三列样本2
      OUTPDF            输出文件
      

      \n");
      q()
}

 

inputfile  = normalizePath(ARGS[1])
outputpdf  =  ARGS[2]

data = read.table(inputfile, head = T , row.names = 1, sep = '\t');
data = log10(data + 0.00000000001)
sample1 = colnames(data)[1]
sample2 = colnames(data)[2]

cor_count = round(cor(data[, 1], data[, 2], method = "pearson" ), 4)
z<- lm(data[, 2]~data[, 1]+1)


pdf(outputpdf);
# plot(x = data[, 1], y = data[, 2], xlab = paste(sample1, "(log10)"), ylab = paste(sample2, "(log10)") , pch = 16, cex = 0.4, xaxs = "i", yaxs = "i")
plot(x = data[, 1], y = data[, 2], xlab = paste(sample1, "(log10)"), ylab = paste(sample2, "(log10)") , pch = 16, xlim = c(0, 5), ylim = c(0, 5), cex = 0.4, xaxs = "i", yaxs = "i")
lines(data[, 1], fitted(z))#添加拟合值对x的散点图并连线
legend("topleft", legend = paste("r = ", cor_count))
dev.off()




