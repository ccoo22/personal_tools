#!/usr/bin/env Rscript

library(docopt)

"Usage: enrich_heatmap.R [options] INPUT OUTPUT 

Arguments:
  INPUT　　the input file name
  OUTPUT   the output file name "-> doc 

opts   <- docopt(doc)
input  <- opts$INPUT
output <- opts$OUTPUT


library(pheatmap)
data <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", quote="")       
ncol <- ncol(data)
nrow <- nrow(data)
longest_colname <- max(nchar(colnames(data)))
longest_rowname <- max(nchar(rownames(data)))

pdf_width <- 0.15*ncol+0.09*longest_rowname
pdf_height <- 0.15*nrow+0.1*longest_colname
# cat(ncol,nrow,longest_colname,longest_rowname,pdf_width,pdf_height,"\n")
# heatmap.2(data.matrix(data), 
# 	      Rowv = F, 
# 	      Colv = F,  
# 	      dendrogram = "none", 
# 	      trace      = "none", 
# 	      scale      = "none", 
# 	      colsep     = 1:ncol(data), 
# 	      rowsep     = 1:nrow(data), 
# 	      key        = F, 
# 	      col        = c("black", "red"),
# 	      lwid       = c(1,10),
# 	      lhei       = c(1,10),
# 	      margins    = c(5,17),
# 	      cexRow     = 1,
# 	      cexCol     = 1,
# 	      )

pdf(output, width = pdf_width, height = pdf_height, onefile=FALSE)
pheatmap(data, 
	    border_color="white",
	     color = c("lightgrey", "red"), 
	     cluster_rows = FALSE, 
	     cluster_cols = FALSE,
	     cellwidth = 10, 
	     cellheight = 10, 
	     legend = FALSE)
dev.off()


