#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: plot_correlation_network.r -i <string> -o <string> -s <string> [-N <numeric> -c <numeric> -p <numeric> -m <string> --offer_data]
Options:
   -s , --sample_list <string>          分组内的样本名，以','分割
   -N <numeric>                         按丰度排序，取前n个otu进行分析，默认分析完整数据 [default: 100] 
   -c , --cor_cutoff <numeric>          cor阈值 [default: 0.8] 
   -p , --pvalue_cutoff <numeric>       pvalue阈值 [default: 0.05] 
   -m , --method <string>       	检验方法:spearman,pearson [default: pearson]
   --offer_data                         输出节点和边的信息表格，可以使用其他软件绘制网络图像
   -i , --input <string>                输入文件，第一列为id，第二列为图像中节点名称，后续每列为一个样本的表达量信息，包含表头
   -o , --output <string>               输出文件路径,结果包含cor.test表格和网络图像" -> doc


# 加载
opts                      <- docopt(doc, version = '网络图绘制\ncontact: 余晨阳 312\n')
sample_list               <- opts$sample_list
N                      	  <- opts$N
cor_cutoff                <- opts$cor
pvalue_cutoff             <- opts$pvalue
method                    <- opts$method
offer_data                <- opts$offer_data
input                     <- opts$input
output                    <- opts$output



if (method != 'spearman' & method != 'pearson'){
  stop("method只能选择spearman或pearson\n")
}

N <- as.numeric(N)
cor_cutoff <- as.numeric(cor_cutoff)
pvalue_cutoff <- as.numeric(pvalue_cutoff)

Data = read.table(input, header = T, sep = "\t", row.names = 1, comment.char = "", quote = "", fill = T)
setwd(output)

sample_names = unlist(strsplit(sample_list, ','))
OTUtype <- as.factor(Data[, 1])
names(OTUtype) = rownames(Data)

Data <- Data[, sample_names, drop = F]
OTUsum = rowSums(Data)
Data <- Data[OTUsum != 0, ]
OTUsum <- OTUsum[OTUsum != 0]
####top N
if (is.numeric(N) & N < nrow(Data)){
  Data <- Data[order(OTUsum, decreasing = TRUE),][1:N,]
}

SamplesOTU = Data[complete.cases(Data),]

# =================================== COR ===================================
# 相关性分析 计算两两OTU之间的相关性

g = combn(row.names(SamplesOTU), 2)
cor_test <- function(t){
  OTUpair = g[,t]
  result = cor.test(unlist(SamplesOTU[OTUpair[[1]],]), unlist(SamplesOTU[OTUpair[[2]],]), method = method, exact = F)
  cor = as.numeric(result$estimate)
  return(c(OTUpair[[1]], OTUpair[[2]], result$p.value, cor))
}

Result = as.data.frame(t(sapply(1:ncol(g), cor_test)), stringsAsFactors = F)
colnames(Result) = c('Source', 'Target', 'Pvalue', 'cor')
Result$Pvalue = as.numeric(Result$Pvalue)
Result$cor    = as.numeric(Result$cor)
write.table(Result,"cortest_result.xls",row.names=F, sep="\t", quote=F)

ResultForPlot = Result[(abs(Result$cor)>cor_cutoff & Result$Pvalue<pvalue_cutoff), ]

if (nrow(ResultForPlot) == 0){
  cat('[Warning]: 没有显著的结果，跳过network绘图！')
  options("show.error.messages" = F)
  stop()
}

# cor转换为logical类型,Ture表示cor>0,FALSE表示cor<0
ResultForPlot$cor = ResultForPlot$cor > 0

id = unique(c(ResultForPlot$Source, ResultForPlot$Target))
Size = OTUsum[id]
Type = as.character(OTUtype[id])

if (offer_data){
  write.table(ResultForPlot[,c(1,2,4)], "edge.csv", row.names=F, sep=",", quote=F)

  node = cbind(id,Size,Type)
  write.table(node, "node.csv", row.names=F, sep=",", quote=F)
}

# ========================== Prepare Data For igraph Plot ==========================
library(igraph)
nodes = data.frame(OTU = id, size = Size, taxon = Type)
links = ResultForPlot

nodes = nodes[order(nodes[,2], decreasing = TRUE),]

net   <- graph_from_data_frame(d = links,vertices = as.character(nodes[[1]]), directed = F)
V(net)$label <- as.character(nodes[[3]])
V(net)$size  <- nodes[[2]]

min_source_size <- min(nodes[[2]])
max_source_size <- max(nodes[[2]])

min_target_size <- 3
max_target_size <- 11
V(net)$size  <- (nodes[[2]] - min_source_size ) / (max_source_size - min_source_size) * (max_target_size - min_target_size) + min_target_size

#Edge_color
E(net)[links$cor]$color <- "#9ACD32"
E(net)[!links$cor]$color <- "#B22222"
E(net)$width = rep(0.3, nrow(links))

mycol = c("#CD1076", "#104E8B", "#DAA520", "#8B008B", "#6E8B3D", "#F08080", "#00B2EE", "#FFEFD5", "#E066FF", "#FF0000", 
              "#0000EE", "#C2601A", "#8B008B", "#8B0000", "#FFC1C1", "#98F5FF", "#FFFF00", "#8A2BE2", "#00FF7F", "#191970", 
              "#00EE00", "#8B4726", "#68228B", "#006400", "#FF4500", "#B4EEB4", "#4169E1", "#FFB90F", "#BF3EFF", "#B3EE3A")

names = levels(nodes[[3]])
cols = mycol[1:length(names)]
names(cols) = names

V(net)$color <- sapply(as.character(nodes[[3]]), function(t){cols[t]})

#plot
pdf("network_with_igraph.pdf")
plot(net, vertex.frame.color=NA,vertex.size = V(net)$size, vertex.label = NA, edge.label = NA, edge.color=E(net)$color, vertex.color=V(net)$color ,edge.curved=.1,layout = layout_nicely)
pos = xy.coords(-0.4,-1.2)
legend(pos,c('Positive correlation', 'Negative correlation'), lty=1,lwd = 5, ncol= 2,col=c("#9ACD32","#B22222"), cex=0.82,bty="n")
pos = xy.coords(-1.2,-1.2)
legend(pos, names, pch=21, ncol= 2, pt.cex = 1.3,pt.bg=cols, cex=0.82,bty="n")
dev.off()
