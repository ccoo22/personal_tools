#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)

"Usage: metabolite_clusterProfiler_metaboanalyst.R  -c <file> --pathway_db <file>   --output_dir <dir> [--download_pathway]
Options:
   -c, --compound <file>          代谢物id列表文件, 没有表头，第一列必须是kegg代谢物id编号
   --pathway_db <file>            metaboanalyst下载的通路数据库，例如： ./metabolite_database/kegg.hsa.rda
   --output_dir <dir>             结果输出目录
   --download_pathway             下载通路图, 默认不下载，有点慢" -> doc


opts              <- docopt(doc, version='甘斌，metaboanalyst数据库的通路分析\n')
compound          <- opts$compound
pathway_db        <- opts$pathway_db 
output_dir        <- opts$output_dir
download_pathway  <- opts$download_pathway
kegg_dir          <- paste0(output_dir, "/KEGG_Pathway_Illustrations")

message("start metabolite_clusterProfiler_metaboanalyst.r")


library(plotly)
library(ggplot2)
library(pathview)


dir.create(output_dir, showWarnings = FALSE)

message("loading")
# 读入代谢物keggid列表
data_compound <- read.table(compound, header = F, sep = "\t" ,stringsAsFactors = F, quote = "", comment.char ="")
data_compound <- unique(toupper(data_compound[, 1])) 

# 加载metpa数据库
load(pathway_db)
current.mset <- metpa$mset.list;  # 所有通路列表，每个列表里是当前通路的代谢物向量
uniq.compound <- unique(unlist(current.mset, use.names=FALSE))  # 所有通路中包含的代谢物id
uniq.count <- length(uniq.compound);  # 代谢物id数量

# 去掉不存在的代谢物id
ora.vec <- data_compound[data_compound %in% uniq.compound]
q.size <- length(ora.vec);
if(q.size == 0)
{
	message("输入的代谢物ID没有匹配到kegg数据库，请仔细核查")
	q()
}

# 开始计数
message("calculating")
hits <- lapply(current.mset, function(x){x[x %in% ora.vec]});  # 每个通路匹配上的代谢物列表
hit.num <- unlist(lapply(hits, function(x){length(x)}), use.names=FALSE);  # 匹配上的数量
set.size <-length(current.mset);  # 通路数量
set.num <- unlist(lapply(current.mset, length), use.names=FALSE);  # 每个通路包含的代谢物数量

# 结果表格准备
res.mat<-matrix(0, nrow=set.size, ncol=12)
res.mat<- as.data.frame(res.mat, check.names = F, stringsAsFactors =F)
rownames(res.mat)<-names(current.mset)
colnames(res.mat)<-c("ID", "Description", "CompoundRatio", "BgRatio", "Pvalue", "-log(p)", "Holm adjust", "FDR", "Impact", "CompoundID", "Count", 'KEGGImage')
# colnames(res.mat)<-c("Total", "Expected", "Hits", "Raw p", "-log(p)", "Holm adjust", "FDR", "Impact");  # 官方原始的结果
	

res.mat[, 'ID']            <-names(current.mset)  # 通路id
res.mat[, 'Description']   <-names(metpa$path.ids)  # 通路描述
res.mat[, 'CompoundRatio'] <-paste(hit.num, q.size, sep = '/')   
res.mat[, 'BgRatio']       <-paste(set.num, uniq.count, sep = '/')   
res.mat[, 'Count']         <-hit.num   
# https://www.genome.jp/kegg-bin/show_pathway?hsa00260/C00631%09red
# 超几何分布
res.mat[, "Pvalue"]   <- phyper(hit.num-1, set.num, uniq.count-set.num, q.size, lower.tail=F)
res.mat[, "-log(p)"] <- -log(res.mat[, "Pvalue"])
  
# adjust for multiple testing problems
res.mat[, "Holm adjust"] <- p.adjust(res.mat[, "Pvalue"], "holm")
res.mat[, "FDR"]         <- p.adjust(res.mat[, "Pvalue"], "fdr")

# 拓扑打分 calculate the sum of importance
# 拓扑打分数据库（通路上，每一个代谢物都有一个分值，求和为1，有两种算法）
if("rbc.list" %in% names(metpa))
{
    # imp.list <- metpa$rbc;  # relative betweenness centrality
    # # imp.list <- metpa$dgr;  # out degree centrality
    res.mat[, "Impact"]     <- mapply(function(x, y){sum(x[y])}, metpa$rbc.list, hits);
}

# compoundid list
res.mat[, "CompoundID"] <- unlist( lapply(hits, function(x){paste(x, collapse = '/')}  ) )

# 通路图在线阅览地址
color_code <- unlist( lapply(hits, function(x){paste(x, collapse = '%09red/')     }  ) )
res.mat[, "KEGGImage"]  <- paste("https://www.genome.jp/kegg-bin/show_pathway?", res.mat[, 'ID'], "/", color_code, "%09red", sep = "")

# 去掉没有结果的通路
res.mat <- res.mat[hit.num > 0, , drop=FALSE];
# res.mat <- res.mat[!is.na(res.mat[, "Impact"]), , drop=FALSE]

# 按照p值、impact排序
if(nrow(res.mat) > 1){
    ord.inx <- order(res.mat[, "Pvalue"])
    res.mat <- res.mat[ord.inx,]
}

# 输出
kegg_enrichment <- paste(output_dir, "kegg_enrichment.xls", sep="/")
message("output: ", kegg_enrichment)
write.table(res.mat, kegg_enrichment, sep="\t", quote=F, row.names=F)

#############
# dotplot 图绘制
#############
kegg_pdf <- paste(output_dir, "kegg_dotplot.pdf", sep="/")
message("output: ", kegg_pdf)

# 最多绘制p值最显著的前10个
data_plot <- res.mat[order(res.mat[, "Pvalue"]), ]
if(nrow(data_plot) > 10) data_plot <- data_plot[1:10, ]

data_plot$CompoundRatio = unlist(lapply(as.character(data_plot[, 'CompoundRatio']), function(x){ res = as.integer( unlist(strsplit(x, split="/")) ); res[1]/res[2] })) # 计算比例
data_plot = data_plot[order(data_plot$CompoundRatio), ]  # 从小到大排序
data_plot$Description = factor(data_plot$Description, data_plot$Description)  # 转换成factor格式，从而确保绘图顺序

pdf(kegg_pdf, width = 10)
p <- ggplot(data_plot, aes(x=Description, y=CompoundRatio)) + geom_point(aes(col=Pvalue, size=Count)) + scale_color_continuous(low="red", high="blue", name = 'Pvalue', guide=guide_colorbar(reverse=TRUE)) + xlab(NULL) + scale_size(range=c(3, 8)) + coord_flip()
print(p)
dev.off()

###############
#  气泡图绘制
###############
# 当有impact统计时，才绘制气泡图
if( "rbc.list" %in% names(metpa))
{
    kegg_bubble_pdf <- paste(output_dir, "kegg_bubble.pdf", sep="/")
    kegg_bubble_html <- paste(output_dir, "kegg_bubble.html", sep="/")
    message("output: ", kegg_bubble_pdf)

    res.tmp <- res.mat[!is.na(res.mat[, "Impact"]), , drop=FALSE]
    x <- res.tmp[, 'Impact']
    y <- res.tmp[, '-log(p)']
    names(x) <- res.tmp[, 'Description']
    names(y) <- res.tmp[, 'Description']

    # first sort values based on p
    inx <- order(y, decreasing= T)
    x <- x[inx]
    y <- y[inx]

    # 圆圈的半径计算
    # set circle size according to impact
    # take sqrt to increase spread out
    sqx <- sqrt(x)
    min.x<- min(sqx, na.rm = TRUE)
    max.x <- max(sqx, na.rm = TRUE)
    if(min.x == max.x){ # only 1 value
        max.x = 1.5*max.x
        min.x = 0.5*min.x
    }

    maxR <- (max.x - min.x)/40
    minR <- (max.x - min.x)/160
    radi.vec <- minR+(maxR-minR)*(sqx-min.x)/(max.x-min.x)  
  
    # set background color according to y
    # 设定圆圈的背景色
    bg.vec <- heat.colors(length(y));

    # 绘制
    pdf(kegg_bubble_pdf)
    plot(x, y, type="n", axes=F, xlab="Pathway Impact", ylab="-log(p)")
    axis(1);
    axis(2);
    grid(col="blue");
    symbols(x, y, add = TRUE, inches = F, circles = radi.vec, bg = bg.vec, xpd=T)
    dev.off()

    # html 交互性图, 文件有点大 5M左右，暂时不绘制
    # data_plot = data.frame(pvalue = y, impact=x, txt=names(x) ) 
    # p <- plot_ly(data_plot, x = ~impact, y = ~pvalue, color =~pvalue,  size = ~impact,   text = ~txt) 
    # htmlwidgets::saveWidget(p, kegg_bubble_html, selfcontained = T)
}

###############
#  kegg图绘制
###############
if(download_pathway)
{
    dir.create(kegg_dir, showWarnings = FALSE)
    message("output: kegg pathway image ")
    setwd(kegg_dir)  # 跳转

    # 只下载p值显著的
    res.sig <- res.mat[res.mat[, 'Pvalue'] < 0.05, ]
    for(row in 1:nrow(res.sig))
    {   
        message("process: ", row, "/", nrow(res.sig))
        # （1）准备id输入数据
        compound <- unlist(strsplit(res.sig[row, 'CompoundID'], "/"))
        cpd.data <- rep(1, length(compound))
        names(cpd.data) <- compound

        # （2）准备通路id
        pathid_full <- res.sig[row, 'ID']
        species <- sub("\\d+","",pathid_full)  # 物种名称
        pathid  <- sub(species,"",pathid_full)  # 通路编号

        # （3）下载通路图、绘图
        pv.out <- pathview(cpd.data = cpd.data, pathway.id = pathid, species = species, out.suffix = "mark_compound", new.signature = F, plot.col.key = F, low='blue', mid='gray', high='red', map.null = F)
    
        # XML文件删除,客户用不到
        file.remove(paste0(pathid_full, ".xml"))
 
    }
}


message("finish metabolite_clusterProfiler_metaboanalyst.r")
