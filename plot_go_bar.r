#!/home/genesky/software/r/4.0.3/bin/Rscript
# 参数解析
parse_args <- function(){
    library("argparse")
    parser <- ArgumentParser(description='GO 富集分析绘图')

    parser$add_argument('-i', '--input', type='character', metavar='file', required=TRUE,
        help = "GO 富集分析结果文件，有表头，且必须含有：Description pvalue Count Type")

    parser$add_argument( '-o', '--output', type='character', metavar='file', required=TRUE,
        help = "pdf 输出文件")

    parser$add_argument('-t', '--top', type='integer', metavar='number', default=10,
        help = "对文件中pvalue最显著的前 n 个绘图 [默认: %(default)s]")

    args <- parser$parse_args()

    return(args)
}
args = parse_args()


# 读入
data_input <- read.table(args[['input']], header = T, sep = "\t" , check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
need_heads = c('Description', 'pvalue', 'Count', 'Type')

if(sum(!need_heads %in% colnames(data_input)) > 0){
    message("[Error] 输入文件缺失必要的表头。它必须含有表头：", paste(need_heads, collapse=','))
    q()
}


############
# BP/CC/MF top10 绘图
############
types       <- unique(data_input$Type)
DataForPlot <- list()
for(j in 1:length(types)){
    DataForPlot[[j]] <-  data_input[grep(types[j],data_input$Type),]
    # pvalue 排序
    DataForPlot[[j]] =  DataForPlot[[j]][order(DataForPlot[[j]]$pvalue), ]
    # 展示最显著的topN个
    TopN <- args[['top']]
    if(nrow(DataForPlot[[j]]) < TopN){
        DataForPlot[[j]] <- DataForPlot[[j]][c(1:nrow(DataForPlot[[j]])), need_heads]
        DataForPlot[[j]] <- DataForPlot[[j]][order(DataForPlot[[j]]$Count, decreasing = T), ]
    }else{
        DataForPlot[[j]]=DataForPlot[[j]][c(1:TopN),need_heads]
        DataForPlot[[j]]=DataForPlot[[j]][order(DataForPlot[[j]]$Count,decreasing = T),]
    }          
}
GO <- do.call(rbind, DataForPlot)
 
color  <- colorRampPalette(c("Goldenrod1","Tomato1","MediumPurple4"))(length(unique(GO$Type)))
color1 <- rep(color,as.numeric(table(GO$Type)))

pdf(args[['output']], width=15, height=8)
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

 


