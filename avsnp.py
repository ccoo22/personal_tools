#!/home/genesky/software/python/3.8.12/bin/python3
import argparse
import logging
import multiprocessing
import os
import re
import sys

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(
        description="pca plot with interactive output", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--dbsnp', type=str, required=True,
                        help="从dbsnp官网下载的vcf文件。例如 /home/genesky/database_new/dbsnp/155/GCF_000001405.25.gz。需要有.tbi索引")

    parser.add_argument('--rename_chrs', type=str, required=True,
                        help="需要重命名的染色体，以及需要提取的染色体。第一列：vcf染色体名称。第二列：重命名的染色体(也是最终要提取的染色体)。例如 /home/genesky/database_new/dbsnp/155/hg19.chrs.txt")

    parser.add_argument('--genome', '-g', type=str, required=True,
                        help="参考基因组 .fa 文件（注意：要包含 rename_chrs 第二列 对应的染色体编号。同时，也要有 .fai索引文件，由samtools faidx创建）。 例如 /home/genesky/database_new/self_build_database/ucsc/hg19_modify/hg19_modify.fa")
    
    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help="结果输出目录。自动创建。")

    parser.add_argument('--prefix', '-p', type=str, default="avsnp.txt",
                        help="输出文件的前缀，例如 hg19_avsnp155.txt  [默认：%(default)s]")
    
    parser.add_argument('--thread', '-t', type=int, default=25,
                        help="线程数量。 主要用于色体拆分处理加速  [默认：%(default)s]")
    
    parser.add_argument('--bcftools', type=str, default='/home/genesky/software/bcftools/1.16/bin/bcftools',
                        help="bcftools软件 [默认：%(default)s]")
    
    parser.add_argument('--convert2annovar', type=str, default='/home/genesky/software/annovar/2018Apr16/convert2annovar.pl',
                        help="annovar工具箱中的 convert2annovar.pl工具 [默认：%(default)s]")
    
    parser.add_argument('--createindex', type=str, default='/home/pub/bin/NGS/chip/GATK4/tools/personal/createAnnovarIDX/createindex.pl',
                        help="annovar索引创建工具 [默认：%(default)s]")  
    args = parser.parse_args()

    return args


def get_chrs(file):
    '''提取染色体编号'''
    chr_raws = []
    chr_renames = []
    with open(file, 'r') as f:
        for line in f:
            if not re.search('\w', line):
                continue
            chr_raw, chr_rename = line.strip().split("\t")
            chr_raws.append(chr_raw)
            chr_renames.append(chr_rename)
    return chr_raws, chr_renames


def create_dir(*dirs):
    '''批量创建目录'''
    for dir in dirs:
        if not os.path.exists(dir) :
            os.makedirs(dir)


def norm_annovar_convert(input, output):
    '''vcf文件转annovar格式'''
    fh_out = open(output, 'w')
    with os.popen(f'gzip -cd {input}') as fh:
        for line in fh:
            if not re.search('\w', line):
                continue
            if re.match('^#', line):
                continue
            values = line.split('\t', 5)
            chrom, pos, id, ref, alt = values[:5]
            start = int(pos)
            end = start + len(ref) - 1
            fh_out.write("\t".join([chrom, str(start), str(end), ref, alt, id]) + "\n")
    fh_out.close()


def indels_annovar_convert(input, output):
    '''convert2annovar生成的文件转annovar格式'''
    fh_out = open(output, 'w')
    with open(input, 'r') as fh:
        for line in fh:
            if not re.search('\w', line):
                continue
            if re.match('^#', line):
                continue
            values = line.split('\t', 8)
            chrom, start, end, ref, alt = values[:5]
            id = values[7]
            end = int(start) + len(ref) - 1
            fh_out.write("\t".join([chrom, str(start), str(end), ref, alt, id]) + "\n")
    fh_out.close()


def sort_file(input, output):
    '''按照坐标排序'''
    os.system(f"sort -k2n,2 -k3n,3 {input} > {output}")

def split_by_chr_and_order(output_dir, *inputs):
    '''文件合并、染色体拆分、排序'''
    log.info("数据库合并、染色体拆分")
    fh_chrom = {}
    with os.popen(f"""cat {' '.join(inputs)}""") as fh:
        for line in fh:
            chrom, tmp = line.split('\t', 1)
            if chrom not in fh_chrom:
                # 创建染色体输出句柄
                fh_chrom[chrom] = open(os.path.join(output_dir, chrom), 'w')
            
            fh_chrom[chrom].write(line)
    
    # 关闭句柄
    for chrom in fh_chrom.keys():
        fh_chrom[chrom].close()

    # 排序
    log.info("排序")
    result = []
    pool = multiprocessing.Pool(processes=10)
    for chrom in fh_chrom.keys():
        input = os.path.join(output_dir, chrom)
        output = os.path.join(output_dir, chrom + ".sort")
        res = pool.apply_async(sort_file, (input, output,))
        result.append(res)  # 先把结果对象保存
    pool.close()
    pool.join()
    for res in result:
        res.get()


def handle_each_chrom(chr_raw, chr_rename, input_vcf, output_dir, args):

    # bcftools index 主要目的是检查是否文件truncated
    chrom_dir = os.path.join(output_dir, chr_rename)
    create_dir(chrom_dir)
    
    log.info(f"[process 1] chr{chr_rename} 染色体名称替换、提取")
    rename_vcf = os.path.join(chrom_dir, f"chr{chr_rename}-1.rename.vcf.gz")
    os.system(f"{args.bcftools} annotate --rename-chrs {args.rename_chrs} -x INFO -r {chr_raw}  {input_vcf} -O z -o {rename_vcf}")
    os.system(f"{args.bcftools} index {rename_vcf} --tbi")

    log.info(f"[process 2] chr{chr_rename} 多态位点拆分")
    norm_vcf = os.path.join(chrom_dir, f"chr{chr_rename}-2.norm.vcf.gz")
    os.system(f"{args.bcftools} norm -m -both -f {args.genome} {rename_vcf} -O z -o {norm_vcf}")
    os.system(f"{args.bcftools} index {norm_vcf} --tbi")  
    
    log.info(f"[process 3] chr{chr_rename} 原始annovar文件制作")
    norm_to_annovar = os.path.join(chrom_dir, f"chr{chr_rename}-3.norm.to.annovar")
    norm_annovar_convert(norm_vcf, norm_to_annovar)
    
    log.info(f"[process 4] chr{chr_rename} indels 提取")
    indels_vcf = os.path.join(chrom_dir, f"chr{chr_rename}-4.indels.vcf.gz")
    os.system(f"{args.bcftools} view --types indels {norm_vcf} -O z -o {indels_vcf}")
    os.system(f"{args.bcftools} index {indels_vcf} --tbi")
    
    log.info(f"[process 5] chr{chr_rename} indels 转 annovar 原始格式")
    indels_to_annovar_raw = os.path.join(chrom_dir, f"chr{chr_rename}-5.indels.annovar.raw")
    os.system(f"perl {args.convert2annovar} --format vcf4old --includeinfo --comment --outfile {indels_to_annovar_raw} {indels_vcf}")

    log.info(f"[process 6] chr{chr_rename} annovar 特有的indels 文件制作")
    indels_to_annovar_clean = os.path.join(chrom_dir, f"chr{chr_rename}-6.indels.annovar.clean")
    indels_annovar_convert(indels_to_annovar_raw, indels_to_annovar_clean)
    
    log.info(f"[process 7] chr{chr_rename} snp/indels 合并、排序")
    snp_indel_final = os.path.join(chrom_dir, f"chr{chr_rename}-7.final")
    os.system(f" cat {norm_to_annovar} {indels_to_annovar_clean} |sort -k2n,2 -k3n,3 > {snp_indel_final}")

if __name__ == "__main__":
    args = set_and_parse_args()
    # 绝对路径
    args.output_dir = os.path.abspath(args.output_dir)
    # 临时目录
    tmp_dir = os.path.join(args.output_dir, 'tmp')
    create_dir(args.output_dir, tmp_dir)
    
    # 染色体
    chr_raws, chr_renames = get_chrs(args.rename_chrs)

    # 拆分染色体处理后续
    result = []
    pool = multiprocessing.Pool(processes=args.thread)
    for chr_raw, chr_rename in zip(chr_raws, chr_renames):
        res = pool.apply_async(handle_each_chrom, (chr_raw, chr_rename, args.dbsnp, tmp_dir, args))
        result.append(res)  # 先把结果对象保存
    pool.close()
    pool.join()
    for res in result:
        res.get()
        
    log.info(f"[process 8] 染色体数据合并")
    final_file = os.path.join(args.output_dir, args.prefix)
    os.system(f"cd {tmp_dir} && cat {' '.join([ './' + chrom +'/chr' + chrom + '-7.final'  for chrom in  chr_renames])} > {final_file}")
    
    log.info("[process 9] annovar索引创建")
    os.system(f"perl {args.createindex} {final_file} 100")
    
    log.info("完成")
    
