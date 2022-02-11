#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: cluster_sample.r  -i <file> -o <dir> [--sample_file <file> --gene_file <file>  --mad_cutoff <numeric>  --maxk <int>  --reps <int> --no_normalize --do_log2 --cluster_method <string> --distance_method <string> ]
Options:
    -i, --input <file>                   一行对应一个特征，每一列对应一个样本，不要有缺失，有表头。第一列是特征名
    -o, --output_dir <dir>               输出路径
    --sample_file <file>                 指定分析的样本，一列数据，一行一个样本,没有表头。 如果不指定，默认用 input的所有样本
    --gene_file <file>                   指定分析的基因，一列数据，一行一个基因,没有表头。如果不指定，默认用 input的所有基因
    --mad_cutoff <numeric>               对输入的基因做MAD过滤，仅保留 > mad 值的基因，例如 0.5 。默认不过滤。 
    --maxk <int>                         最大聚类数量 [default: 9]
    --reps <int>                         resamplings, number of subsamples [default: 50]
    --no_normalize                       取消归一化。建议做归一化，对每一行数据做归一化处理，使不同特征之间也具有可比性，否则聚类效果不会很好看
    --do_log2                            对输入的数据做log2处理，然后再做归一化、cluster。 
    --cluster_method <string>            聚类算法： 'hc' heirarchical(hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, 'kmdist' for k-means upon distance matrices (former km option), or a function that returns a clustering [default: hc]
    --distance_method <string>           距离算法： 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski' or custom distance function. [default: pearson]
   " -> doc

opts                    <- docopt(doc, version = '基于特征，对样本聚类\n')
input                   <- opts$input
output_dir              <- opts$output_dir
no_normalize            <- opts$no_normalize
do_log2                 <- opts$do_log2
cluster_method          <- opts$cluster_method
distance_method         <- opts$distance_method
sample_file             <- opts$sample_file
gene_file               <- opts$gene_file
maxk                    <- as.integer(opts$maxk)
reps                    <- as.integer(opts$reps)
mad_cutoff              <- opts$mad_cutoff

library(ConsensusClusterPlus, quietly = TRUE)

 
# 读取数据，行为样本
message("read data")
data <- read.table(input, sep='\t', header = TRUE, row.names = 1, check.names=FALSE, stringsAsFactors = F) 
data <- as.matrix(data)
genes = rownames(data)
samples = colnames(data)


if(! is.null(sample_file))
{
    sample_file_data <- read.table(sample_file, sep='\t', header = FALSE, check.names=FALSE, stringsAsFactors = F) 
    samples = sample_file_data[,1]
    lost_samples = samples[!samples %in% colnames(data)]
    if(length(lost_samples) > 0)
    {
        message("指定的样本不存在：", paste(lost_samples, collapse=','))
        q()
    }
    message("指定分析 ", length(samples), " 个样本")
}else{
    message("使用input中的所有", ncol(data), " 个样本进行分析")
}
if(! is.null(gene_file))
{
    gene_file_data <- read.table(gene_file, sep='\t', header = FALSE, check.names=FALSE, stringsAsFactors = F) 
    genes = gene_file_data[,1]
    lost_genes = genes[!genes %in% rownames(data)]
    if(length(lost_genes) > 0)
    {
        message("指定的基因不存在：", paste(lost_genes, collapse=','))
        q()
    }
    message("指定分析 ", length(genes), " 个基因")
}else{
    message("使用input中的所有", nrow(data), " 个基因进行分析")
}

data = data[genes, samples]

if(! is.null(mad_cutoff))
{
    message('做 MAD 过滤基因， mad = ', mad_cutoff)
    mad_cutoff = as.numeric(mad_cutoff)
    mads = apply(data,1,mad)
    mad_genes = names(mads)[mads > mad_cutoff]
    if(length(mad_genes) <=2 )
    {
        message('经 MAD 过滤后，剩余基因不足2个，无法继续分析')
        q()
    }else{
        message('经 MAD 过滤后，剩余基因 数量： ', length(mad_genes))
        data = data[mad_genes, ]

        mad_genes = data.frame(gene=mad_genes, stringsAsFactors = F)
        file = paste0(output_dir, "/mad_filter_genes.txt")
        write.table(mad_genes, file, col.names = F, row.names=F, sep = "\t", quote=F)
    }
}else{
    message('基因不做 MAD 过滤')
}

if(do_log2)
{
    message('log2 转换')
    data = log2(data + 1)  # 要提前加 1, 防止存在负值
}

if(!no_normalize)
{   
    message("normalize data")
    data = sweep(data,1, apply(data,1,median,na.rm=T))
}


message("cluster")
results = ConsensusClusterPlus(d=data,maxK=maxk,reps=reps,pItem=0.8,pFeature=1, title=output_dir,clusterAlg=cluster_method,distance=distance_method,seed=1262118388.71279,plot="png",writeTable=FALSE)

message("output cluster file")
for(k in 2:maxk)
{
    sample_cluster = data.frame(sample = names(results[[k]][["consensusClass"]]), cluster = results[[k]][["consensusClass"]])
    file = paste0(output_dir, "/sample_cluster.k", k, ".txt")
    write.table(sample_cluster, file, col.names = T, row.names=F, sep = "\t", quote=F)
}

