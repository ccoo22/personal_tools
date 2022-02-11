sheet	Title	Description
meQTL	snps	snp ID
meQTL	gene	gene ID。 注意：本结果中，连续型数据的表头前缀统一命名为gene，实际上可能是甲基化、ATAC peakID等数据，请已实际数据说明为准。
meQTL	statistic	pvalue的t值
meQTL	pvalue	线性回归pvalue
meQTL	FDR	FDR矫正
meQTL	beta	线性回归beta值
meQTL	sampleNum	当前snp-gene配对结果，非缺失样本数量
meQTL	snp_SampleNum	当前snp-gene配对结果，非缺失样本中，snp样本野生型、杂合、纯合样本的数量。注意：以当前分析群体的MAF对应的allele作为突变等位基因。
meQTL	snp_chr	snp的染色体
meQTL	snp_pos	snp的位置
meQTL	gene_chr	gene的染色体
meQTL	gene_left_pos	gene的左边界位置
meQTL	gene_right_pos	gene的右边界位置
meQTL	gene_mean	gene的数据均值
meQTL	gene_median	gene的数据中位数
meQTL	snp_gene_distance	snp_pos - gene_left
