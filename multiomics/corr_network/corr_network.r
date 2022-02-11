#!/home/genesky/software/r/4.0.3/bin/Rscript

library(docopt)

"Usage: corr_network.r  --table_file <file> --sample_group <file> --case <string> --control <string> --output_dir <pdf file> [--method <string> --width <int>  --height <int> --filter_pvalue <numeric> --filter_cor <numeric> --N_node <integer> --N_edge <integer>  --show_vertex_label <degree>  --vertex_size_max <integer>  --vertex_size_min <integer>  --rlib <dir>]

Options:
    --table_file <file>            绘制多组学相关性热图的必要文件。含表头，每一行表示一个组学相关信息，第一列组学名称，第二列表达量文件。
                                   表达量文件：含表头，每一行是一个特征，每一列是一个样本，第一行是样本名，第一列是特征名。
    --sample_group <file>          样本分组文件。含表头，第一列样本名称，第二列分组名称，流程对--case --control指定组的样本进行后续相关性计算及绘图
    --case <string>                指定case组名称，即--sample_group第二列的某一分组名称
    --control <string>             指定control组名称，即--sample_group第二列的某一分组名称
    --output_dir <file>            输出路径
    --method <string>              相关性计算方法，可选 'pearson', 'spearman', 'kendall' [default: spearman]
    --filter_pvalue <numeric>      pvalue值过滤 [default: 0.05]
    --filter_cor <numeric>         相关性系数值过滤 [default: 0.8]
    --N_node <integer>             绘制网络图节点的数量 [default: 200]
    --N_edge <integer>             绘制网络图边的数量 [default: 3000]
    --show_vertex_label <degree>   显示节点的度>=某阈值的节点标签 [default: 3]
    --vertex_size_max <integer>    将相对丰度数值映射到绘制网络图节点的最大值 [default: 8]
    --vertex_size_min <integer>    将相对丰度数值映射到绘制网络图节点的最小值 [default: 4]
    --width <int>                  PDF宽度 [default: 10]
    --height <int>                 PDF高度 [default: 10]
    --rlib <dir>                   R包路径 [default: /home/genesky/software/r/4.0.3/lib64/R/library]" -> doc

opts           <- docopt(doc, version='igraph包 -- 绘制相关性网络图\n')
table_file     <-  opts$table_file
sample_group   <-  opts$sample_group
case           <-  opts$case
control        <-  opts$control
output_dir     <-  opts$output_dir
method         <-  opts$method
filter_pvalue  <- as.numeric(opts$filter_pvalue)
filter_cor     <- as.numeric(opts$filter_cor)
N_node         <- as.integer(opts$N_node)
N_edge         <- as.integer(opts$N_edge)
show_vertex_label <- as.integer(opts$show_vertex_label)
vertex_size_max   <- as.integer(opts$vertex_size_max)
vertex_size_min   <- as.integer(opts$vertex_size_min)
width          <-  as.integer(opts$width)
height         <-  as.integer(opts$height)  
rlib           <-  opts$rlib

# 加载R包
.libPaths(rlib)
library(psych) ##corr.test
library(igraph)
set.seed(123)

# 测试数据
# table_file = "/home/dongxj/work/research/plot/network/test/table_file.txt"
# sample_group = "/home/dongxj/work/research/plot/heatmap/test/sample_group.txt"
# case    = "Post"
# control = "Baseline"
# method  = 'spearman'
# output_dir = "/home/dongxj/work/research/plot/network/test/result"
# N_node  = 135
# N_edge  = 100


# 读入数据
table_data  <- read.table(table_file,  sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
sample_data <- read.table(sample_group, sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
colnames(table_data)[1]  <-  'exp_txt'
colnames(sample_data) <-  'group'

# 读取指定分组的样本
cases <- rownames(sample_data[sample_data[,'group']==case,,drop=F])
controls <- rownames(sample_data[sample_data[,'group']==control,,drop=F])
samples <- c(cases, controls)

# 读取表达量数据，特征存储到network_all_nodes:三列，c('enrich_group', 'class', 'relative_abundance')
data_name_exp <- list()
network_all_nodes  <- data.frame()
for (i in 1:nrow(table_data)) {
    data_name <- rownames(table_data)[i]
    data_exp  <- read.table(table_data[,'exp_txt'][i],  sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
    data_samples = colnames(data_exp)

    lost_samples = samples[!samples %in% data_samples]
    if (length(lost_samples) > 0)
    {
        message("[Error] case control分组文件中的样本在", data_name, "表达量文件中没有找到 : ", lost_samples)
        q()
    }
    ## 表达量数据，对--case --control指定组的样本进行后续相关性计算
    exp <- data_exp[,samples]
    data_name_exp[[i]] <- exp
    
    ## network.nodes  点文件
    data_log2FoldChange <- log(apply(exp[,cases], 1, mean)/apply(exp[,controls], 1, mean),2)
    data_relabundance   <- apply(exp, 1, mean)
    network_node <- data.frame(log2FoldChange = data_log2FoldChange, relative_abundance = data_relabundance)
    network_node <- network_node[rownames(exp),,drop=F] ##node文件与exp文件特征名称顺序保持一致
    network_node$enrich_group[network_node$log2FoldChange > 0] <- paste('enriched in', case, 'samples')
    network_node$enrich_group[network_node$log2FoldChange < 0] <- paste('enriched in', control, 'samples')
    network_node$class <- data_name
    network_node <- network_node[,c('enrich_group', 'class', 'relative_abundance')]
    network_all_nodes <- rbind(network_node, network_all_nodes)
}

# 两两计算相关性值、pvalue值，分别输出；相关性信息存储到network_all_links
network_all_links <- data.frame()
for (i in 1:(nrow(table_data)-1)) {
    for (j in (i+1):nrow(table_data) ){
        data1_name <- rownames(table_data)[i]
        data2_name <- rownames(table_data)[j]
        data1_exp    <- data_name_exp[[i]]
        data2_exp    <- data_name_exp[[j]]

        # 计算相关性值、pvalue值，分别输出
        result = matrix(nrow = nrow(data1_exp) * nrow(data2_exp), ncol = 7);       
        result_estimate = matrix(nrow = nrow(data1_exp), ncol = nrow(data2_exp))
        result_pvalue   = matrix(nrow = nrow(data1_exp), ncol = nrow(data2_exp))
        
        rownames(result_estimate) = rownames(data1_exp)
        colnames(result_estimate) = rownames(data2_exp)
        
        rownames(result_pvalue) = rownames(data1_exp)
        colnames(result_pvalue) = rownames(data2_exp)

        row = 0
        for(data1 in rownames(data1_exp))
        {
            for(data2 in rownames(data2_exp))
            {
                row  <- row + 1
                data <- cbind(name1 = as.numeric(data1_exp[data1, samples]), name2 = as.numeric(data2_exp[data2, samples]))
                data <- data.frame(data[complete.cases(data),]) #去掉缺失值
                data <- data[apply(data == "", 1, sum) == 0, ]  #去掉空值               
                nmiss <- nrow(data) #nmiss必须>1
                pmt <- NA
                cmt <- NA
                if(nmiss > 2){
                    cor <- corr.test(data[,1], data[,2], method = method, adjust = "none") ### corr.test()进行特征间相关性分析（按列数据进行相关性分析，文件内容需要转置）;如果不矫正，即adjust ="none"，其相关系数、P值和cor.test()结果一样。      
                    cmt <- cor$r
                    pmt <- cor$p
                }
                result[row, ] <- c(data1_name, data1, data2_name, data2, nmiss, pmt, cmt)
                result_estimate[data1, data2] <- cmt
                result_pvalue[data1, data2] <- pmt
            }
        }
        # 输出 相关性值、pvalue值
        result_estimate <- data.frame(id = rownames(result_estimate), result_estimate, check.names=F)
        result_pvalue   <- data.frame(id = rownames(result_pvalue), result_pvalue, check.names=F)
        corr_file     <- paste(output_dir, '/',data1_name, '_', data2_name, '_corr.txt', sep = '')
        write.table(result_estimate, corr_file,  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        pvalue_file   <- paste(output_dir, '/',data1_name, '_', data2_name, '_pvalue.txt', sep = '')
        write.table(result_pvalue, pvalue_file,  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

        ##所有的link
        network_all_links <- rbind(result, network_all_links)
    }
}
colnames(network_all_links) <- c('source_class', 'source', 'target_class', 'target',  'NMISS', 'Pvalue', 'Estimate')

##准备网络图数据，过滤条件：pvalue < filter_pvalue, cor > filter_cor，根据相关性值降序排序
network_all_links$Pvalue   <- as.numeric(network_all_links$Pvalue)
network_all_links$Estimate <- as.numeric(network_all_links$Estimate)
network_edge <-  network_all_links[network_all_links$Pvalue < filter_pvalue & abs(network_all_links$Estimate) > filter_cor, ]
network_edge <-  network_edge[order(abs(network_edge$Estimate),decreasing=T),] #相关性值降序
network_edge$direction[network_edge$Estimate > 0] <- 'plus correlation'
network_edge$direction[network_edge$Estimate < 0] <- 'minus correlation'
network_edge <- network_edge[,c('source','target', 'Pvalue', 'Estimate','direction')]
##过滤后的边数据输出
network_fiter_edge_file  <- paste(output_dir, '/network.links', sep = '')
write.table(network_edge, network_fiter_edge_file,  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
##过滤后的点数据输出
network_node <- network_all_nodes[unique(unlist(network_edge[, c("source", "target")])),]
network_fiter_node_file  <- paste(output_dir, '/network.nodes', sep = '')
network_output_node <- data.frame(name = rownames(network_node), network_node, check.names = F)
write.table(network_output_node, network_fiter_node_file,  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

##创建绘图数据
N_edge <- min(N_edge, nrow(network_edge)) ##确定边的数量
nodes  <- c()
i_row  <- 0
while(length(nodes) < N_node && i_row < N_edge){
	i_row <- i_row + 1
	nodes <- unique(c(nodes, unlist(network_edge[i_row, c("source", "target")])))
}
edge <- network_edge[1:i_row, , drop = F] ##提取边
node <- network_node[nodes, , drop = F]  ##提取点
net  <- graph_from_data_frame(d = edge, vertices = nodes, directed = F) ##创建graph

##设置点的颜色、形状
case_enrich_group <- paste('enriched in', case, 'samples')
control_enrich_group <- paste('enriched in', control, 'samples')
V(net)[node$enrich_group == case_enrich_group]$color    <- "red"
V(net)[node$enrich_group == control_enrich_group]$color <- "green"

#自定义三角形
mytriangle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }

    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v]
    }

    symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
            stars=cbind(vertex.size, vertex.size, vertex.size),
            add=TRUE, inches=FALSE)
}
add_shape("triangle", clip=shapes("circle")$clip,  plot=mytriangle)

#自定义菱形
mydiamond <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }

  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           symbols(x=x, y=y, bg=bg,
                   stars=matrix(size, nrow=1, ncol=nor),
                   add=TRUE, inches=FALSE)
         })
}
add_shape("diamond", clip=shape_noclip, plot=mydiamond, parameters=list(vertex.norays=4))

class  <- unique(node$class)
shapes <- c("circle", "triangle", "square", "diamond")
pch    <- c(1, 2, 0, 5)
for (i in 1:length(class))
{
    V(net)[node$class == class[i]]$shape <- shapes[i] ##形状
    V(net)[node$class == class[i]]$pch   <- pch[i]    ##pch点的样式
}

##设置边的颜色
E(net)[edge$direction == "plus correlation"]$color  <-  "#c53f2d"
E(net)[edge$direction == "minus correlation"]$color <-  "#26479b"

## 点的大小（相关性值）
degree <- degree(net)
relabundance <- node$relative_abundance
max_relabundance <- max(node$relative_abundance)
min_relabundance <- min(node$relative_abundance)

## 点的大小映射
new_vertex_size <- vertex_size_min +(vertex_size_max - vertex_size_min)*(relabundance - min_relabundance)/(max_relabundance - min_relabundance)
new_vertex_size <- as.numeric(format(new_vertex_size, digits = 2))

###绘图
locs <- layout_with_fr(net) # 统一节点布局
##第一张图：显示节点的度>=某阈值的节点标签
outputfile  <- paste(output_dir, '/corr_network.pdf', sep = '')
cairo_pdf(outputfile, width = width, height = height)
par(mai = c(2, 0.1, 0.1, 0.1))

plot(net, vertex.label.color='black', # 点标签字体颜色
            vertex.label.cex = 0.6,     # 点标签字体大小
            vertex.size = new_vertex_size,   # 点的大小
            vertex.label = ifelse(new_vertex_size>=show_vertex_label, rownames(node), ""), # 显示节点的度>=某阈值的节点标签
            vertex.shape = V(net)$shape, # 点的形状
            vertex.color=V(net)$color,   # 点的填充颜色
            edge.label = NA,             # 边的标签
            edge.color = E(net)$color, # 边的颜色
            edge.curved = 0,           # 边是否弯曲 
            layout = locs)   # 布局   

xy = par("usr") #绘图区域左下角和右上角的坐标
##plot绘图, pch点的样式，legend
pch <- as.numeric(unique(V(net)$pch))
legend(x=xy[1]+xinch(0.2,2), y=xy[3]-yinch(-0.15,4), legend = class, pch = pch, ncol = 1, pt.cex = 2, cex = 1.1, bty = "n")
## 分组legend
legend(x=xy[1]+xinch(0.2,2), y=xy[3]-yinch(0.5,4), legend = c(case_enrich_group, control_enrich_group), col = c("red", "green"), pch = c(21,21), ncol= 1, pt.cex = 2, pt.bg = c("red", "green"), cex = 1.2, bty = "n")
## 正负相关legend
legend(x=xy[1]+xinch(0,2), y=xy[3]-yinch(1,4), legend = c("plus correlation", "minus correlation" ), col = c("#c53f2d", "#26479b"), lty = c(1,1), ncol = 1, cex = 1.2,bty = "n")
dev.off()

##第二张图：不显示节点标签
outputfile2  <- paste(output_dir, '/corr_network_hide_name.pdf', sep = '')
cairo_pdf(outputfile2, width = width, height = height)
par(mai = c(2, 0.1, 0.1, 0.1))

plot(net, vertex.label.color='black', # 点标签字体颜色
            vertex.label.cex = 0.6,     # 点标签字体大小
            vertex.size = new_vertex_size,   # 点的大小
            vertex.label = "",           # 不显示点标签
            vertex.shape = V(net)$shape, # 点的形状
            vertex.color=V(net)$color,   # 点的填充颜色
            edge.label = NA,             # 边的标签
            edge.color = E(net)$color, # 边的颜色
            edge.curved = 0,           # 边是否弯曲 
            layout = locs)   # 布局   

xy = par("usr") #绘图区域左下角和右上角的坐标
##plot绘图, pch点的样式，legend
pch <- as.numeric(unique(V(net)$pch))
legend(x=xy[1]+xinch(0.2,2), y=xy[3]-yinch(-0.15,4), legend = class, pch = pch, ncol = 1, pt.cex = 2, cex = 1.1, bty = "n")
## 分组legend
legend(x=xy[1]+xinch(0.2,2), y=xy[3]-yinch(0.5,4), legend = c(case_enrich_group, control_enrich_group), col = c("red", "green"), pch = c(21,21), ncol= 1, pt.cex = 2, pt.bg = c("red", "green"), cex = 1.2, bty = "n")
## 正负相关legend
legend(x=xy[1]+xinch(0,2), y=xy[3]-yinch(1,4), legend = c("plus correlation", "minus correlation" ), col = c("#c53f2d", "#26479b"), lty = c(1,1), ncol = 1, cex = 1.2,bty = "n")
dev.off()
