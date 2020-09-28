#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: cluster_sample.r  -i <file> -o <dir> [--maxk <int>  --reps <int> --no_normalize --do_log2 --cluster_method <string> --distance_method <string> ]
Options:
    -i, --input <file>                   一行对应一个特征，每一列对应一个样本，不要有缺失，有表头。第一列是特征名
    -o, --output_dir <dir>               输出路径
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
maxk                    <- as.integer(opts$maxk)
reps                    <- as.integer(opts$reps)
 
library(ConsensusClusterPlus, quietly = TRUE)

 
# 读取数据，行为样本
message("read data")
data <- read.table(input, sep='\t', header = TRUE, row.names = 1, check.names=FALSE) 
data <- as.matrix(data)
if(do_log2)
{
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
