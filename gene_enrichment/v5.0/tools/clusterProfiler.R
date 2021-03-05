#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: clusterProfiler.r  -i <file> -o <dir> --org <string> --orgdb <file> --orgdb_keggid <string> [--rlib <dir>]

Options:
    -i, --input <file>        输入基因列表文件，没有表头，第一列数据必须是geneSymbol, 第二列用于标记基因的颜色，可以不写，默认红色，值仅允许up/down，up红色、down蓝色。 
                              注意：基因大小写是严格区分的，一定要保证与NCBI一致，否则无法转换ENTREZID，进而在富集分析中丢失。
    -o, --output_dir <dir>    结果输出目录（脚本自动创建）
    --org <string>            kegg 物种简称，例如 人对应 hsa, 小鼠对应 mmu, 大鼠对应 rno。 可以在下述网址查找： https://www.genome.jp/kegg/catalog/org_list.html
    --orgdb <file>            orgdb 数据库，通过 AnnotationHub 下载得到，我整理了下载教程。 人数据库示例： /home/genesky/database_new/self_build_database/clusterprofiler/human.hsa.orgdb
    --orgdb_keggid <string>   orgdb 数据库中，哪一个数据列能对应上kegg gene id。 
                              目前支持： 
                                      ENTREZID    : 只要下面的链接能打开，则建议用这个 http://rest.kegg.jp/conv/org/ncbi-geneid
                                      其他         ： 当ENTREZID用不了的时候，可以使用orgdb其他列, 能够直接对应KEGG中基因编号的列名称， 例如 ALIAS。
    --rlib <dir>              R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，GO/KEGG 富集分析\n')
input           <- opts$input
output_dir      <- opts$output_dir
org             <- opts$org
orgdb           <- opts$orgdb
orgdb_keggid    <- opts$orgdb_keggid
rlib            <- opts$rlib
.libPaths(rlib)

dir.create(output_dir, showWarnings = FALSE)
 
library(AnnotationDbi)
library(clusterProfiler)
library(topGO)
library(rlist)
library(DOSE)  # enrichDO 需要用到这个包、数据库也存在DOSE安装目录中，仅适用于人
library(pathview)  # 通路图下载工具

############
# 加载基因、orgdb
############
message("load gene/orgdb")
data_input      <- read.table(input, stringsAsFactors=F, fill = T)
gene_symbol     <- data_input[,1]
gene_symbol_down <- c()
if(ncol(data_input) > 1) gene_symbol_down = data_input[grepl('down', data_input[,2], ignore.case = T), 1] # 下调基因symbol
org_db   <- loadDb(file = orgdb)
org_db
 
############
# 基因SYMBOL 转 ENTREZID
############
message("convert SYMBOL to ENTREZID")
eg   <- bitr(gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb=org_db)
rownames(eg) = eg[[2]]
gene <- eg[[2]]  # ENTREZID


############
# GO富集分析
############
message("start GO analysis")
message('CC')
ego_cc <- enrichGO(gene          = gene,
                   OrgDb         = org_db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   readable      = TRUE)
message('BP')
ego_bp <- enrichGO(gene          = gene,
                   OrgDb         = org_db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   readable      = TRUE)
message('MF')
ego_mf <- enrichGO(gene          = gene,
                   OrgDb         = org_db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   readable      = TRUE)
############
# BP/CC/MF 绘图
############
go_bp_pdf          <- paste(output_dir, "go.bp.pdf", sep="/")
go_cc_pdf          <- paste(output_dir, "go.cc.pdf", sep="/")
go_mf_pdf          <- paste(output_dir, "go.mf.pdf", sep="/")

pdf(go_bp_pdf)
plotGOgraph(ego_bp)
dev.off()

pdf(go_cc_pdf)
plotGOgraph(ego_cc)
dev.off()

pdf(go_mf_pdf)
plotGOgraph(ego_mf)
dev.off()

############
# GO 富集分析结果汇总输出
############
cc <- as.data.frame(ego_cc)
if(nrow(cc) >= 1) cc$Type <- "CC"

bp <- as.data.frame(ego_bp)
if(nrow(bp) >= 1) bp$Type <- "BP"

mf <- as.data.frame(ego_mf)
if(nrow(mf) >= 1) mf$Type <- "MF"
go_enrichment   <- paste(output_dir, "go_enrichment.xls", sep="/")
res <- rbind(bp, cc, mf)

if(nrow(res) >= 1) write.table(res, go_enrichment, sep="\t", quote=F, row.names=F)

############
# BP/CC/MF top10 绘图
############
ego <- res
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

go_bar <- paste(output_dir, "GO_barplot.pdf", sep="/")
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

############
# KEGG分析
############
message("start KEGG analysis")

# 基因搜索方式设置
keytype = 'ncbi-geneid'  # 默认ENTREZID
if(orgdb_keggid != 'ENTREZID')
{
    message("orgdb_keggid is not ENTREZID, convert symbol to kegg id as parameter 'orgdb_keggid'")
    eg   <- bitr(gene_symbol, fromType="SYMBOL", toType=orgdb_keggid, OrgDb=org_db)  
    rownames(eg) = eg[[2]]  # 注意：因为有的SYMBOL存在多个别名，导致eg的第一列存在重复的SYMBOL
    gene <- eg[[2]]   
    keytype = 'kegg'
}
# 下调的geneid，为通路图做准备
geneid_down = eg[eg[, 1] %in% gene_symbol_down, 2]  


# 富集分析
kk <- enrichKEGG(gene         = gene,
                 organism     = org,
                 keyType      = keytype,
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)
kk_enrich <- as.data.frame(kk)
# 找回geneID对应的geneSymbol
kk_enrich$geneSymbol <- unlist(lapply(1:length(kk_enrich$geneID), function(row){  
                                        paste(eg[unlist(strsplit(kk_enrich$geneID[row], "/")), 1], collapse = "/" )
                                        }))
# 添加在线查看地址，标记颜色
kk_enrich$kegg_link <- unlist(lapply(1:length(kk_enrich$geneID), function(row){
                                          geneid = unlist(strsplit(kk_enrich$geneID[row], "/"))
                                          down = geneid[geneid %in% geneid_down]  # 下调
                                          up   = geneid[!geneid %in% geneid_down]  # 上调

                                          gene_color = c()
                                          if(length(down) > 0) gene_color = c(gene_color, paste(down, '%09blue', sep = ''))
                                          if(length(up) > 0)   gene_color = c(gene_color, paste(up,    '%09red', sep = ''))
                                          gene_color = paste(gene_color, collapse = '/')

                                          # 最终地址
                                          paste0('https://www.kegg.jp/kegg-bin/show_pathway?', kk_enrich$ID[row], '/', gene_color)
                                        }))
# 结果输出
if(nrow(kk_enrich) > 0)
{   
    # （1）富集通路点图
    kegg_pdf <- paste(output_dir, "kegg_dotplot.pdf", sep="/")
    pdf(kegg_pdf,height=7,width=12)
    p <-dotplot(kk)
    print(p)
    dev.off()
    
    #（2）表格输出
    kegg_enrichment <- paste(output_dir, "kegg_enrichment.xls", sep="/")
    write.table(kk_enrich, kegg_enrichment, sep="\t", quote=F, row.names=F)

    # （3）下载通路图
    pathway_down_dir = paste(output_dir, "KEGG_Pathway_Illustrations", sep="/")
    dir.create(pathway_down_dir)
    setwd(pathway_down_dir)

    id_sigs   = kk_enrich$ID[kk_enrich$pvalue < 0.05]  # 显著通路
    count_sig = length(id_sigs)
    count     = 0
    for(id in id_sigs)
    {
        count = count + 1
        message("[dowload] ", count, "/", count_sig, "  ", id)

        # （3.1）需要标记颜色的geneid/geneid类型
        geneid           = unlist(strsplit(kk_enrich[id, 'geneID'], '/'))
        gene.data        = c(rep(1, length(geneid)))
        names(gene.data) = geneid
        gene.idtype = ifelse(orgdb_keggid != 'ENTREZID', 'KEGG', "entrez" ) 

        # （3.2）确定下调基因，设定颜色, -1 下调/ 1 上调
        if(length(geneid_down) > 0)
        {
            geneid_down_in_kegg = geneid_down[geneid_down %in% geneid]
            gene.data[geneid_down_in_kegg] = -1
        }

        # （3.3） 下载
        pathview(gene.data = gene.data, cpd.data = NULL, pathway.id =  sub(org,"",id), species = org, gene.idtype = gene.idtype, out.suffix = "mark_gene", map.null = F, new.signature = FALSE,  kegg.native = T, plot.col.key = FALSE, low = list(gene = "blue", cpd = "green"), high = list(gene = "red", cpd = "yellow"))
        
        # （3.4）XML文件删除,客户用不到
        file.remove(paste0(id, ".xml"))
    }
}

############
# DOSE Disease Ontology 分析
############
if(org == 'hsa')
{
    message("start DOSE Disease Ontology analysis")

    do = enrichDO(gene, 
             ont = "DO", 
             pvalueCutoff = 1, 
             qvalueCutoff = 1,
             readable = TRUE, # 自动把geneid 转成 symbol。因为只有人才能做，所以数据库简单
             )
    do_enrich <- as.data.frame(do)

    if(nrow(do_enrich) >= 1)
    {
        do_pdf <- paste(output_dir, "do_dotplot.pdf", sep="/")
        pdf(do_pdf, height=7, width=12)
        p <-dotplot(do)
        print(p)
        dev.off()

        do_enrichment <- paste(output_dir, "do_enrichment.xls", sep="/")
        write.table(do_enrich, do_enrichment, sep="\t", quote=F, row.names=F)
    }
}


