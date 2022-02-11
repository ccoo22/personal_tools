#!/home/genesky/software/python/3.6.7/bin/python3
import argparse
import os
import sys
import re
import pysam


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description='模拟二代测序数据，同时引入一定比例的突变\n注意：软件还不是很稳定，产生的序列PCR重复效果不好，使比对去重不平衡，最终使得到的突变频率与设定的相差甚远',formatter_class=argparse.RawTextHelpFormatter) 
    parser.add_argument('-g', '--genome', type=str, required=True, help="参考基因组fasta文件。\n注意：需要有samtools索引 .fai") 
    parser.add_argument('-s', '--snv', type=str, required=True, help='模拟数据的突变信息文件，特殊格式, 表头名称固定，由7列数据构成,tab分隔，必须是：chr	pos	ref	alt	depth	alt_freq	extend	mean_fragment_length\nchr  染色体\npos 位置\nref 参考等位基因\nalt 突变等位基因\ndepth 生成的reads数量（不能保证extend后的区域每个位点都完全覆盖）\nalt_freq 期望的突变等位基因频率\nextend 想要在突变位点两侧扩展多大的区域，建议 200，太长的话，目标位点的覆盖估计可能不足\nmean_fragment_length 每条reads的平均DNA片段长度，建议200') 
    parser.add_argument('-o', '--output_dir', type=str, required=True, help="结果输出目录, 脚本自动创建") 
    parser.add_argument('-p', '--prefix', type=str, required=True, help="fastq文件名前缀") 
    parser.add_argument('--art_illumina', type=str, default='/home/genesky/software/art/2016.06.05/bin/art_illumina', help="art_illumina软件路径（默认： /home/genesky/software/art/2016.06.05/bin/art_illumina）") 
    args = parser.parse_args()
    return args


def parse_snv_file(file):
    """   解析SNV文件    """
    heads_need = ['chr', 'pos', 'ref', 'alt', 'depth', 'alt_freq', 'extend', 'mean_fragment_length']
    snv_info = {}
    with open(file, 'r') as fh:
        heads = fh.readline().rstrip().split('\t')
        if  heads != heads_need:
            print('[Error] SNV文件的表头必须是', heads_need)
            sys.exit(1)
        for line in fh:
            if not re.search('\w', line):
                continue
            chrom, pos, ref, alt, depth, alt_freq, extend, mean_fragment_length = line.rstrip().split('\t')
            ref = ref.upper()
            alt = alt.upper()
            title = '_'.join([chrom, pos, ref, alt])

            pos = int(pos)
            depth = int(depth)
            alt_freq = float(alt_freq)
            extend = int(extend)
            mean_fragment_length = int(mean_fragment_length)
            
            snv_info[title] = {}
            snv_info[title]['title'] = title
            snv_info[title]['chrom'] = chrom
            snv_info[title]['pos'] = pos
            snv_info[title]['start'] = pos - extend  # 区域扩展
            snv_info[title]['end'] = pos + extend
            snv_info[title]['ref'] = ref
            snv_info[title]['alt'] = alt
            snv_info[title]['depth'] = depth
            snv_info[title]['depth_ref'] = int(depth * (1-alt_freq))
            snv_info[title]['depth_alt'] = int(depth * alt_freq)
            snv_info[title]['alt_freq'] = alt_freq
            snv_info[title]['extend'] = extend
            snv_info[title]['mean_fragment_length'] = mean_fragment_length
            snv_info[title]['region'] = chrom + ':' + str(snv_info[title]['start']) + '-' + str(snv_info[title]['end'])
            
    return snv_info


def create_fastq(dir, file_name, seq_name, seq, depth, mean_fragment_length, art_illumina):
    """  根据需要，创建fastq """
    os.chdir(dir)

    # 参考基因组生成
    fh = open(file_name, 'w')
    fh.write('>' + seq_name + '\n')
    fh.write(seq + '\n')
    fh.close()

    # 运行
    cmd = art_illumina + ' -ss HSXt -i ' + file_name + ' -M  -rs 1234567890 -l 150  -s 80 -o ' + file_name + ' -f ' + str(depth) + ' -m ' + str(mean_fragment_length) + ' > ' + file_name + '.run.log' + ' 2>&1'
    # print(cmd)
    os.system(cmd)


def count_coverage(aln_r1, aln_r2, extend, cover_txt):
    """ 统计目标位点覆盖深度 """
    count = {}
    for file in [aln_r1, aln_r2]:
        with open(file, 'r') as fh:
            for line in fh:
                if(not re.search('^>', line)):
                    continue
                chrom, name, pos, strand = re.split('\s+', line.rstrip())
                name = re.split('\/', name)[0]
                pos = int(pos)
                start = pos
                end = pos + 150 - 1
                # 因为目标位点处于正中间，故，不需要再考虑+-方向问题
                if(extend >= start and extend <= end):
                    count[name] = line
    
    # reads名称输出，方便后续查看
    fh = open(cover_txt, 'w')
    for name in count.keys():
        fh.write(name + '\t' + count[name])
    fh.close()

    return len(count.keys())


if __name__ == '__main__':
    args = set_and_parse_args()
    sample = args.prefix

    # 路径准备
    output_dir = os.path.abspath(args.output_dir)
    sample_dir = output_dir + '/' + sample
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    
    # SNV/基因组文件处理
    snv_info = parse_snv_file(args.snv)
    fa_genome = pysam.FastaFile(args.genome)

    # bed 覆盖信息
    bed = output_dir + "/" + sample + ".realign.bed"
    fh_bed = open(bed, 'w')

    # 拟合数据在SNV上的覆盖情况
    cover = output_dir + '/' + sample + ".coverage.txt"
    fh_cover = open(cover, 'w')
    heads = ['chrom', 'pos', 'ref', 'alt', 'depth', 'alt_freq', 'extend', 'mean_fragment_length', 'simulation_depth', 'simulation_depth_ref', 'simulation_depth_alt', 'simulation_alt_freq']
    fh_cover.write('\t'.join(heads) + '\n')

    # fastq制作
    for title in snv_info.keys():
        print(title)
        this_dir = sample_dir + '/' + title
        if not os.path.exists(this_dir):
            os.makedirs(this_dir)

        # 提取序列
        seq = fa_genome.fetch(region=snv_info[title]['region']).upper()
        seq_alt = seq[0:snv_info[title]['extend']] + snv_info[title]['alt'] + seq[snv_info[title]['extend']+1:]
        base_ref =  seq[snv_info[title]['extend']]
        if base_ref != snv_info[title]['ref']:
            print('[ERROR] 发现提取的序列中参考碱基与SNV文件中的ref不一致，请仔细确认', title)
        
        # 生成fastq
        # ref
        ref_reads_count = 0
        if(snv_info[title]['depth_ref'] > 0):
            create_fastq(this_dir, 'ref.fa', title + '_ref', seq, snv_info[title]['depth_ref'], snv_info[title]['mean_fragment_length'], args.art_illumina)
            ref_reads_count = count_coverage(this_dir + '/' + 'ref.fa1.aln', 
                                             this_dir + '/' + 'ref.fa2.aln',
                                             snv_info[title]['extend'],
                                             this_dir + '/' + 'ref.snv.cover.txt'
                                             )
        # alt
        alt_reads_count = 0
        if(snv_info[title]['depth_alt'] > 0):
            create_fastq(this_dir, 'alt.fa', title + '_alt', seq_alt, snv_info[title]['depth_alt'], snv_info[title]['mean_fragment_length'], args.art_illumina)
            alt_reads_count = count_coverage(this_dir + '/' + 'alt.fa1.aln', 
                                             this_dir + '/' + 'alt.fa2.aln',
                                             snv_info[title]['extend'],
                                             this_dir + '/' + 'alt.snv.cover.txt'
                                             )
        # 拟合数据SNV统计
        snv_info[title]['simulation_depth'] = ref_reads_count + alt_reads_count
        snv_info[title]['simulation_depth_ref'] = ref_reads_count
        snv_info[title]['simulation_depth_alt'] = alt_reads_count
        snv_info[title]['simulation_alt_freq'] = alt_reads_count / snv_info[title]['simulation_depth']

        # bed输出
        fh_bed.write('\t'.join([snv_info[title]['chrom'], str(snv_info[title]['start']), str(snv_info[title]['end']), '+', title]) + '\n')

        # 统计输出
        values = map(lambda x: str(snv_info[title][x]) ,heads)
        fh_cover.write('\t'.join(values) + '\n')

    fh_bed.close()
    fh_cover.close()

    # 合并、压缩 fastq
    fastq_r1 = output_dir + '/' + sample + '_R1.fastq.gz'
    fastq_r2 = output_dir + '/' + sample + '_R2.fastq.gz'
    os.system("cat " + sample_dir + "/*/*.fa1.fq | gzip  > " + fastq_r1)
    os.system("cat " + sample_dir + "/*/*.fa2.fq | gzip  > " + fastq_r2)

    # 结束
    print('finished')
    print(fastq_r1)
    print(fastq_r2)