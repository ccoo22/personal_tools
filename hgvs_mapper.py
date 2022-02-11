#!/home/genesky/software/python/3.9.4/bin/python3
import os
import sys
import argparse
import re
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.variantmapper

def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(
        description="根据HGVS coding信息（例如：NM_175629.2:c.G2645A），找到基因组坐标 \n", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input', '-i', type=str, required=True,
                        help="输入文件，一行对应一个HGVS coding信息")

    parser.add_argument('--output', '-o', type=str,
                        help="输出文件，如果不写，则打印到屏幕上")

    parser.add_argument('--genome', '-g', type=str, default='human_hg38',
                        help="参考基因组版本， human_hg19/human_hg38")

    parser.add_argument('--refgene', '-r', type=str, required=True,
                        help="UCSC refgene.txt文件，用于查找输入的转录本的染色体编号。 例如：/home/genesky/database_new/ucsc/refgene/hg38-rmalt/annovar/hg38-rmalt_refGene.txt /home/genesky/database_new/ucsc/refgene/hg19/annovar/hg19_refGene.txt")

    parser.add_argument('--show_process', action="store_true",
                        help="显示进度")  # 布尔值
    args = parser.parse_args()

    return args


def get_chrom_ncbi_map(genome):
    hg19_map = {}
    hg19_map["chr1"] = "NC_000001.10"
    hg19_map["chr2"] = "NC_000002.11"
    hg19_map["chr3"] = "NC_000003.11"
    hg19_map["chr4"] = "NC_000004.11"
    hg19_map["chr5"] = "NC_000005.9"
    hg19_map["chr6"] = "NC_000006.11"
    hg19_map["chr7"] = "NC_000007.13"
    hg19_map["chr8"] = "NC_000008.10"
    hg19_map["chr9"] = "NC_000009.11"
    hg19_map["chr10"] = "NC_000010.10"
    hg19_map["chr11"] = "NC_000011.9"
    hg19_map["chr12"] = "NC_000012.11"
    hg19_map["chr13"] = "NC_000013.10"
    hg19_map["chr14"] = "NC_000014.8"
    hg19_map["chr15"] = "NC_000015.9"
    hg19_map["chr16"] = "NC_000016.9"
    hg19_map["chr17"] = "NC_000017.10"
    hg19_map["chr18"] = "NC_000018.9"
    hg19_map["chr19"] = "NC_000019.9"
    hg19_map["chr20"] = "NC_000020.10"
    hg19_map["chr21"] = "NC_000021.8"
    hg19_map["chr22"] = "NC_000022.10"
    hg19_map["chrX"] = "NC_000023.10"
    hg19_map["chrY"] = "NC_000024.9"

    hg38_map = {}
    hg38_map["chr1"] = "NC_000001.11"
    hg38_map["chr2"] = "NC_000002.12"
    hg38_map["chr3"] = "NC_000003.12"
    hg38_map["chr4"] = "NC_000004.12"
    hg38_map["chr5"] = "NC_000005.10"
    hg38_map["chr6"] = "NC_000006.12"
    hg38_map["chr7"] = "NC_000007.14"
    hg38_map["chr8"] = "NC_000008.11"
    hg38_map["chr9"] = "NC_000009.12"
    hg38_map["chr10"] = "NC_000010.11"
    hg38_map["chr11"] = "NC_000011.10"
    hg38_map["chr12"] = "NC_000012.12"
    hg38_map["chr13"] = "NC_000013.11"
    hg38_map["chr14"] = "NC_000014.9"
    hg38_map["chr15"] = "NC_000015.10"
    hg38_map["chr16"] = "NC_000016.10"
    hg38_map["chr17"] = "NC_000017.11"
    hg38_map["chr18"] = "NC_000018.10"
    hg38_map["chr19"] = "NC_000019.10"
    hg38_map["chr20"] = "NC_000020.11"
    hg38_map["chr21"] = "NC_000021.9"
    hg38_map["chr22"] = "NC_000022.11"
    hg38_map["chrX"] = "NC_000023.11"
    hg38_map["chrY"] = "NC_000024.10"

    if genome == 'human_hg19':
        return hg19_map
    else:
        return hg38_map


if __name__ == '__main__':
    args = set_and_parse_args()
    chrom_ncbi_map = get_chrom_ncbi_map(args.genome)


    # 读入refgene
    refgene = {}
    with open(args.refgene, 'r') as fh:
        for line in fh:
            line.strip()
            rowid, mrna, chrom, tmp = line.split('\t', 3)
            # 仅保留常规染色体
            if not re.search('_', chrom):
                refgene[mrna] = chrom

    # 读入hgvs文件
    hgvs_input = {}
    with open(args.input, 'r') as fh:
        for row, line in enumerate(fh):
            line = line.strip()
            line = re.sub('\s', '', line)
            if re.match('\w', line):
                hgvs_input[row] = {}
                hgvs_input[row]['hgvs'] = line
                hgvs_input[row]['mrna'] = line.split(':')[0].split('.')[0]  # 转录本名称，去除版本号
                if hgvs_input[row]['mrna'] in refgene:
                    hgvs_input[row]['chrom'] = refgene[hgvs_input[row]['mrna']]  # ucsc 染色体
                    hgvs_input[row]['chrom_ncbi'] = chrom_ncbi_map[hgvs_input[row]['chrom']] # NCBI染色体

    # 转换
    hgvsparser = hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    mapper = hgvs.variantmapper.VariantMapper(hdp)
    for row in sorted(hgvs_input.keys()):
        if args.show_process:
            print(hgvs_input[row]['hgvs'])
        if 'chrom_ncbi' in hgvs_input[row]:
            var = hgvsparser.parse_hgvs_variant(hgvs_input[row]['hgvs'])
            try:
                protein = mapper.c_to_p(var)
            except Exception as result:
                print("python错误 %s" % result)
                hgvs_input[row]['protein'] = "转换失败，请确认是否写错了转录组ID"
            else:
                hgvs_input[row]['protein'] = protein.format()
            
            try:
                genome = mapper.c_to_g(var, hgvs_input[row]['chrom_ncbi'])
            except Exception as result:
                print("python错误 %s" % result)
                hgvs_input[row]['genomic'] = "转换失败，请确认是否写错了转录组ID、或者选错了基因组版本"
            else:
                hgvs_input[row]['genomic'] = genome.format()
        else:
            hgvs_input[row]['chrom'] = ''
            hgvs_input[row]['chrom_ncbi'] = ''
            hgvs_input[row]['genomic'] = ''
            hgvs_input[row]['protein'] = ''
    
    # 输出
    headers = ['chrom', 'chrom_ncbi', 'hgvs', 'genomic', 'protein']

    if args.output != None:
        with open(args.output, 'w') as fh:
            fh.write('\t'.join(headers) + '\n')
            for row in sorted(hgvs_input.keys()):
                values = [ hgvs_input[row][header] for header in headers]
                fh.write('\t'.join(values) + '\n')
    else:
        print('\t'.join(headers))
        for row in sorted(hgvs_input.keys()):
            values = [ hgvs_input[row][header] for header in headers]
            print('\t'.join(values))



    

 
