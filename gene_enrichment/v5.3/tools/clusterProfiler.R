#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: clusterProfiler.r  -i <file> -o <dir> --org <string> --orgdb <file> --orgdb_keggid <string> [--skip_go  --skip_kegg --skip_do  --rlib <dir>]

Options:
    -i, --input <file>        输入基因列表文件，没有表头，第一列数据必须是geneSymbol, 第二列用于标记基因的颜色，可以不写，默认红色，值仅允许up/down，up红色、down蓝色。 
                              注意：基因大小写是严格区分的，一定要保证与NCBI一致，否则无法转换ENTREZID，进而在富集分析中丢失。
    -o, --output_dir <dir>    结果输出目录（脚本自动创建）
    --org <string>            kegg 物种简称，例如 人对应 hsa, 小鼠对应 mmu, 大鼠对应 rno。 可以在下述网址查找： https://www.genome.jp/kegg/catalog/org_list.html
    --orgdb <file>            orgdb 数据库，通过 AnnotationHub 下载得到，我整理了下载教程。支持的物种查询： https://annotationhub.bioconductor.org/species 
                              人数据库示例： /home/genesky/database_new/self_build_database/clusterprofiler/human.hsa.orgdb
                           
    --orgdb_keggid <string>   orgdb 数据库中，哪一个数据列能对应上kegg gene id。 
                              目前支持： 
                                      ENTREZID    : 只要后面的链接能打开（注意org替换为物种缩写），则用这个。 http://rest.kegg.jp/conv/org/ncbi-geneid
                                      其他         ： 当ENTREZID用不了的时候，可以使用orgdb其他列, 这个列中的内容必须是KEGG中使用的基因编号， 例如 ALIAS。
    --skip_go                 不做GO分析
    --skip_kegg               不做kegg分析
    --skip_do                 不做do分析
    --rlib <dir>              R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，GO/KEGG 富集分析\n')
input           <- opts$input
output_dir      <- opts$output_dir
org             <- opts$org
orgdb           <- opts$orgdb
orgdb_keggid    <- opts$orgdb_keggid
skip_go         <- opts$skip_go
skip_kegg       <- opts$skip_kegg
skip_do         <- opts$skip_do
rlib            <- opts$rlib
.libPaths(rlib)
# input = '/home/pub/bin/NGS/chip/GATK4/tools/personal/gene_enrichment/v5.2/example_gene.txt'
# org = 'hsa'
# orgdb = '/home/genesky/database_new/self_build_database/clusterprofiler/human.hsa.orgdb'
# orgdb_keggid = 'ENTREZID'
dir.create(output_dir, showWarnings = FALSE)
 
library(AnnotationDbi)
library(clusterProfiler)
library(topGO)
library(rlist)
library(DOSE)  # enrichDO 需要用到这个包、数据库也存在DOSE安装目录中，仅适用于人
library(pheatmap)
library(dplyr)
library(KEGGREST)
library(jsonlite)
library(stringr)

############
# 加载基因、orgdb
############
message("load gene/orgdb")
data_input      <- read.table(input, stringsAsFactors=F, fill = T, header = FALSE)
gene_symbol     <- data_input[,1]
gene_symbol_down <- c()  # 下调基因symbol
if(ncol(data_input) > 1) gene_symbol_down = data_input[grepl('down', data_input[,2], ignore.case = T), 1] 
org_db   <- loadDb(file = orgdb)
org_db
 
############
# 基因SYMBOL 转 ENTREZID
############
message("convert SYMBOL to ENTREZID")
gene_map_symbol2entrezid   <- bitr(gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb=org_db)
gene_map_symbol2entrezid   <- gene_map_symbol2entrezid[!duplicated(gene_map_symbol2entrezid$SYMBOL), ]  # 可能存在同一个SYMBOL有多个entrezid
gene_entrezid <- gene_map_symbol2entrezid[, 'ENTREZID']  # ENTREZID


############
# GO富集分析
############
if(!skip_go)
{
    message("start GO analysis")
    message('CC')
    ego_cc <- enrichGO(gene          = gene_entrezid,
                       OrgDb         = org_db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       readable      = TRUE)
    message('BP')
    ego_bp <- enrichGO(gene          = gene_entrezid,
                       OrgDb         = org_db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       readable      = TRUE)
    message('MF')
    ego_mf <- enrichGO(gene          = gene_entrezid,
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

}else{
    message("skip GO analysis")
}


############
# KEGG分析
############
if(!skip_kegg)
{
    message("start KEGG analysis")
    
    # keggid 获取, 两列： symbol kegg
    if(orgdb_keggid == 'ENTREZID')  # ncbi-geneid 的情况下，我们把它转换为kegg-id
    {   
        message('convert entrezid to keggid')
        # entrezid 转 keggid
        gene_map_entrezid2keggid = bitr_kegg(gene_entrezid, 'ncbi-geneid', "kegg", org, drop = TRUE)
        colnames(gene_map_entrezid2keggid)[1] = 'ENTREZID'
        
        # symbol -> keggid 映射关系
        gene_map_symbol2keggid = dplyr::left_join(gene_map_symbol2entrezid, gene_map_entrezid2keggid, by = "ENTREZID")
        gene_map_symbol2keggid = gene_map_symbol2keggid[complete.cases(gene_map_symbol2keggid), ]  # 去缺失
        gene_map_symbol2keggid = gene_map_symbol2keggid[!duplicated(gene_map_symbol2keggid), ]  # 去重复
        rownames(gene_map_symbol2keggid) = gene_map_symbol2keggid[, 'kegg']
    
    }else{  # ALIAS情况下， 直接通过ORGDB获取kegg-id
        message('convert SYMBOL to keggid')
        gene_map_symbol2keggid   <- bitr(gene_symbol, fromType="SYMBOL", toType=orgdb_keggid, OrgDb=org_db)
        colnames(gene_map_symbol2keggid)[2] = 'kegg'
        rownames(gene_map_symbol2keggid) = gene_map_symbol2keggid[, 'kegg']
    }  
    
    # 下调的keggid，为通路图做准备
    keggid_down = gene_map_symbol2keggid[gene_map_symbol2keggid[, 'SYMBOL'] %in% gene_symbol_down, 'kegg']  
    
    # 富集分析
    kk <- enrichKEGG(gene         = gene_map_symbol2keggid[, 'kegg'],
                     organism     = org,
                     keyType      = 'kegg',
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)
    kk_enrich <- as.data.frame(kk)
    # 找回geneID对应的geneSymbol
    kk_enrich$geneSymbol <- unlist(lapply(1:length(kk_enrich$geneID), function(row){  
                                            paste(gene_map_symbol2keggid[unlist(strsplit(kk_enrich$geneID[row], "/")), 'SYMBOL'], collapse = "/" )
                                            }))
    # 添加在线查看地址，标记颜色,基于geneid的
    kk_enrich$keggLink <- unlist(lapply(1:length(kk_enrich$geneID), function(row){
                                              keggid      = unlist(strsplit(kk_enrich$geneID[row], "/"))
                                              keggid_down = keggid[keggid %in% keggid_down]  # 下调
                                              keggid_up   = keggid[!keggid %in% keggid_down]  # 上调
                                              
                                              # 颜色
                                              gene_color = c()
                                              if(length(keggid_down) > 0) gene_color = c(gene_color, paste(keggid_down, '%09blue', sep = ''))
                                              if(length(keggid_up) > 0)   gene_color = c(gene_color, paste(keggid_up,    '%09red', sep = ''))
                                              gene_color = paste(gene_color, collapse = '/')
    
                                              # 最终地址
                                              paste0('https://www.kegg.jp/kegg-bin/show_pathway?', kk_enrich$ID[row], '/', gene_color)
                                            }))
    # 添加KEGG通路的level_a/level_b 的分级注释
    # 实时下载关系数据库
    pathway_brite_txt = keggGet('br:br08901', 'json')
    pathway_brite_raw = parse_json(pathway_brite_txt)
    pathway_brite_df = matrix(nrow=0,ncol=4)
    colnames(pathway_brite_df) = c('level_a', 'level_b', 'level_c', 'pathway_id')
    # level A 循环
    for(level_a in 1:length(pathway_brite_raw$children))
    {
        level_a_name = pathway_brite_raw$children[[level_a]]$name
        level_a_children = pathway_brite_raw$children[[level_a]]$children
        # level B 循环
        for(level_b in 1:length(level_a_children))
        {
            level_b_name = level_a_children[[level_b]]$name
            level_b_children = level_a_children[[level_b]]$children
            # level C 循环
            for(level_c in 1:length(level_b_children))
            {
                level_c_name = level_b_children[[level_c]]$name
                level_c_name_split =  str_split_fixed(level_c_name, "\\s+", n=2)
                pathway_id = paste0(org, level_c_name_split[1])
                pathway_desc = level_c_name_split[2]
                pathway_brite_df = rbind(pathway_brite_df, c(level_a_name, level_b_name, pathway_desc, pathway_id))
            }
        }
    }
    pathway_brite_df  = as.data.frame(pathway_brite_df)
    kk_enrich = dplyr::left_join(kk_enrich, pathway_brite_df[, c('level_a', 'level_b', 'pathway_id')], by = c("ID" = "pathway_id"))
    rownames(kk_enrich) = kk_enrich$ID  # left_join 导致rownames丢失，这里补充上

    # 结果输出
    if(nrow(kk_enrich) > 0)
    {   
        # （1）富集通路点图
        kegg_pdf <- paste(output_dir, "kegg_dotplot.pdf", sep="/")
        pdf(kegg_pdf,height=7,width=12)
        p <-dotplot(kk, x = "GeneRatio", color = "pvalue", showCategory = 10)
        print(p)
        dev.off()
        
        #（2）表格输出
        kegg_enrichment <- paste(output_dir, "kegg_enrichment.xls", sep="/")
        write.table(kk_enrich, kegg_enrichment, sep="\t", quote=F, row.names=F)

        # (3) KEGG热图绘制
        sig_pathway   = kk_enrich$ID[kk_enrich$pvalue < 0.05]
        sig_kegg_gene = sort(unlist(lapply(kk_enrich$geneSymbol[kk_enrich$pvalue < 0.05], function(x){ strsplit(x, '/') })))
        sig_kegg_gene = sig_kegg_gene[!duplicated(sig_kegg_gene)]

        if(length(sig_pathway) > 1)
        {   
            # 构造绘图数据
            data_heatmap = matrix(0, nrow = length(sig_pathway), ncol = length(sig_kegg_gene))
            rownames(data_heatmap) = sig_pathway
            colnames(data_heatmap) = sig_kegg_gene
            for(pathway in sig_pathway)
            {   
                message(pathway)
                gene = unlist(strsplit(kk_enrich[pathway, 'geneSymbol'], '/'))
                data_heatmap[pathway, gene] = 1
            }

            # 绘图
            longest_colname <- max(nchar(colnames(data_heatmap)))
            longest_rowname <- max(nchar(rownames(data_heatmap)))
            pdf_width <- 0.15*ncol(data_heatmap)+0.09*longest_rowname
            pdf_height <- 0.15*nrow(data_heatmap)+0.1*longest_colname
            kegg_heatmap <- paste(output_dir, "kegg_heatmap.pdf", sep="/")
            
            pdf(kegg_heatmap, width = pdf_width, height = pdf_height, onefile=FALSE)
            pheatmap(data_heatmap, 
                border_color="white",
                color = c("lightgrey", "red"), 
                cluster_rows = FALSE, 
                cluster_cols = FALSE,
                cellwidth = 10, 
                cellheight = 10, 
                legend = FALSE)
            dev.off()
        }
    }else{
        message("no enrich result")
    }

}else{
    message("skip KEGG analysis")
}


############
# DOSE Disease Ontology 分析
############
if(org == 'hsa' & !skip_do)
{
    message("start DOSE Disease Ontology analysis")

    do = enrichDO(gene_entrezid, 
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
}else{
    message("skip DOSE Disease Ontology analysis")
}


