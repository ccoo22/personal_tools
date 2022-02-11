#!/home/genesky/software/r/4.0.3/bin/Rscript

library(docopt)
"Usage: enricher.r -g <file> -o <file> -d <file> [--rlib <dir>]
Options:
        -g, --gene <file>            差异基因列表,没有表头，一列基因名称
        -o, --output <dir>           输出路径

        -d, --database <file>        个人创建的GO数据库， Gene_ID  列是基因名称， GO列是基因所属的GO Term，多个GO Term用分号分隔，每个GO Term的格式为 'GO Term ID|GO Term描述|GO Term分类', 另外，务必保证Gene_ID没有重复.
        --rlib <dir>                 R包路径 [default: /home/genesky/software/r/4.0.3/lib64/R/library]" -> doc

opts     <- docopt(doc, version = 'Program : 基于自定义的通路数据库，通路富集分析 \n      甘斌 129\n')

gene_file       <- opts$gene
output_dir      <- opts$output
database_file   <- opts$database

library(clusterProfiler)
library(dplyr)
library(tidyr)
library(rlist)

data_gene = read.table(gene_file, header=F, sep = '\t', check.names=F, comment.char='', stringsAsFactors=F, quote = '')
go_data = read.table(database_file, header=T, sep = '\t', check.names=F, comment.char='', stringsAsFactors=F, quote = '')  
 
if(! 'Gene_ID' %in% colnames(go_data) || !'GO' %in% colnames(go_data))
{
    message('[Error] database文件 必须含有 Gene_ID、GO 列，且符合格式定义!')
    q()
}

# 数据库格式化
database = go_data %>% 
    select(Gene_ID, GO) %>% # 选择Pathway列
    separate_rows(GO, sep = ';', convert = F) %>%  # 按分号拆掉
    filter(GO != '')  %>%  # 却掉空值
    separate(col = GO, into = c('ID', 'Description', 'Ontology'), sep = '[|]', convert=FALSE) # 把GO列按照=拆分一下

# 清理输入的基因
gene = data_gene[, 1] 
gene = gene[!duplicated(gene)]  
gene = gene[gene %in% database$Gene_ID]

# ID Ontology映射关系
map_id_ontology = as.data.frame(database[, c('ID', 'Ontology')])
map_id_ontology = map_id_ontology[!duplicated(map_id_ontology$ID), ]
rownames(map_id_ontology) = map_id_ontology$ID


# 富集分析
enrich <- enricher(gene = gene,
                  TERM2GENE = database[c('ID','Gene_ID')],
                  TERM2NAME = database[c('ID','Description')],
                  pvalueCutoff = 1,
                  pAdjustMethod = 'BH',
                  qvalueCutoff = 1) 

result = as.data.frame(enrich)
result$Type = map_id_ontology[result$ID, 'Ontology']
result$Type[result$Type == 'biological_process'] = 'BP'
result$Type[result$Type == 'molecular_function'] = 'MF'
result$Type[result$Type == 'cellular_component'] = 'CC'

############
# BP/CC/MF top10 绘图
############
    ego <- result
    types       <- unique(ego$Type)
    DataForPlot <- list()
    for(j in 1:length(types)){
        DataForPlot[[j]] <-  ego[grep(types[j],ego$Type),]
        # 展示最显著的topN个
        TopN <- 10
        if(nrow(DataForPlot[[j]]) < TopN){
            DataForPlot[[j]] <- DataForPlot[[j]][c(1:nrow(DataForPlot[[j]])), c(2,5,9,10)]
            DataForPlot[[j]] <- DataForPlot[[j]][order(DataForPlot[[j]]$Count, decreasing = T), ]
        }else{
            DataForPlot[[j]]=DataForPlot[[j]][c(1:TopN),c(2,5,9,10)]
            DataForPlot[[j]]=DataForPlot[[j]][order(DataForPlot[[j]]$Count,decreasing = T),]
        }
    }

    GO <- list.rbind(DataForPlot)     # rbind all elements in a list by row

    go_bar <- paste(output_dir, "go_barplot.pdf", sep="/")
    color  <- colorRampPalette(c("Goldenrod1","Tomato1","MediumPurple4"))(length(unique(GO$Type)))
    color1 <- rep(color,as.numeric(table(GO$Type)))

    pdf(go_bar, width=15, height=8)
    if(max(nchar(GO$Description)) >= 100){
        layout(matrix(c(3,2,1), nrow = 1, byrow = T),widths = c(1.2,0.1,1))
    }else{
        layout(matrix(c(3,2,1), nrow = 1, byrow = T),widths = c(0.8,0.1,1))
    }

    par(mar = c(5,0.1, 1, 3))
    gobar <- barplot(GO$Count, plot=T,
            cex.lab=1, las=1, ps=0.5, border = F,
            xlim=c(0,max(GO$Count)*1.25),
            #xaxt= "n",
            cex.names=0.8, axis.lty=1,
            axes=TRUE,col=color1,horiz=T,
            mgp=c(0,-0.5,-1))
    text(cex = 1, y = gobar, x = GO$Count+max(GO$Count)/35, lab=c(GO$Count))
    #axis(side = 1, at =c(0,max(GO$Count)*1.25),labels = c("0","Gene Number"))
    title(xlab = "Gene Number",line = 1, cex.lab = 1.2)
    legend("topright", legend=unique(GO$Type),bty="n",fill=color,cex = 1.2)
    y1 = par("usr")[3]
    y2 = par("usr")[4]

    par(mar = c(5,0, 1, 0))
    p = log10(GO$pvalue)
    p_bar = barplot(p,horiz=T,col="PaleGreen3",border = NA, xaxt= "n")
    max = format(max(GO$pvalue),scientific=TRUE,digit=2)
    axis(side = 1, line = -1, at = c(0, unique(c(0,ceiling(min(p))))), labels = T)
    title(xlab = "log10(pvalue)",line = 0.91, cex.lab = 1)

    par(mar = c(5,3, 1, 1))
    a = plot(1:5,ylim = c(y1,y2),type = "n", xaxs = "i", yaxs = "i", axes=F, xlab="",ylab="")
    for(i in 1:length(GO$Description)){
        text(x = 5 ,y = gobar[i], labels = GO$Description[i], cex = 1.4, adj = 1)
    }
    dev.off()


#（2）表格输出
go_enrichment   <- paste(output_dir, "go_enrichment.xls", sep="/")
write.table(result, go_enrichment, sep="\t", quote=F, row.names=F)

