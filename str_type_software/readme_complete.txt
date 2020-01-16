Sheet Name	标注与说明	
Reads	target	片段名称
Reads	str	测序STR,格式为：motif(motif数量)
Reads	样本	本片断测序得到的STR对应的reads数量
Freq	样本	本片断测序得到的STR对应的频率。红色表示频率大于分型阈值，橘黄色表示在噪音阈值和分型阈值之间。
SlipRatio	SlipRatio	片段前滑移比例
SlipRatio	Type	前滑移比例获取方式，SampleGenerate表示来自于样本，unitary_quadratic_equation表示来自于一元二次函数拟合
SlipRatio	SampleCount	计算滑移比例的样本数量
SlipRatio	SlipRatio_median	滑移比例中值
SlipRatio	SlipRatio_sd	滑移比例标准差
SlipRatio	SlipRatio_min	滑移比例最小值
SlipRatio	SlipRatio_max	滑移比例最大值
SlipRatio	Detail	滑移比例详情
Reads Correct	校正滑移后的reads数量	
Freq Correct	校正滑移后的频率	
Amplification Ratio	AMPR_mean	扩增比例均值
Amplification Ratio	AMPR_count	用于计算扩增比例的样本量
Amplification Ratio	AMPR_median	扩增比例均值
Amplification Ratio	AMPR_sd	扩增比例标准差
Amplification Ratio	AMPR_min	扩增比例最小值
Amplification Ratio	AMPR_max	扩增比例最大值
Amplification Ratio	Detail	扩增比例详情
Typing Relative	1..N	每个样本在每个片段上的STR分型结果（相对拷贝数）（1..N表示当前分型对应的motif数量）
Typing Absolute	1..N	每个样本在每个片段上的STR分型结果（绝对拷贝数）（1..N表示当前分型对应的motif数量）
Genotype	样本基因型	
Allele	样本等位基因	
