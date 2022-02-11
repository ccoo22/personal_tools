#!/home/genesky/software/python/3.9.4/bin/python3
import os
import sys
import re
import time
import glob
import argparse
import logging
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

project_dir_db = {}
project_dir_db['wes'] = '/home/pub/output/WES'
project_dir_db['wgs'] = '/home/pub/output/WGS'
project_dir_db['wts'] = '/home/pub/output/WTS'
project_dir_db['mgs'] = '/home/pub/output/MGS'
project_dir_db['target'] = '/home/pub/output/Target'
project_dir_db['methyltarget'] = '/home/pub/output/Methylation'
project_dir_db_string = '\n'.join([ f"{project_type} {dir}" for project_type, dir in project_dir_db.items() ])

def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "genesky  项目垃圾数据清理。\n\n建议先运行check模式，检查待清理文件是否正确，然后再运行clean模式。\n注意事项：是否有清理权限\n"
            "运行示例：\n"
            "(1) 使用mysql, 检查2个月之前的WES垃圾\n./clean.py --month 2\n\n"
            "(2) 使用mysql, 清理2个月之前的WES垃圾\n./clean.py -m clean --month 2\n\n"
            "(3) 使用time,  检查3个月之前的WES垃圾\n./clean.py -s time --month 3 \n\n"
            "(4) 使用mysql, 检查指定目录的WES垃圾\n./clean.py -i /home/pub/output/Target \n\n"
            "(5) 使用mysql, 检查WTS垃圾\n./clean.py -t wts \n\n"
            "(6) 清理指定WES项目\n./clean.py -s project -t wes --project 21B0206A_WES,21AATI003_WES \n\n"
        ))

    parser.add_argument('--project_type', '-t', type=str, default='wes', choices=['wes', 'wgs', 'wts', 'mgs', 'target', 'methyltarget'],
                        help="数据类型. 每种类型要清理哪些文件已经定死了。 [wes]")

    parser.add_argument('--strategy', '-s', type=str, default='mysql', choices=['mysql', 'time', 'project'],
                        help="清理策略，共有三种。[mysql]\nmysql: 检测 month 个月之前的，且项目合同已完成的进行清理，只能在2号服务器运行； \ntime : 根据month指定的时间，清理超过该时间的项目； \nproject : 根据project指定的项目列表做清理。")

    parser.add_argument('--mode', '-m', type=str, default='check', choices=['check', 'clean'],
                        help="运行模式。 check 模式： 仅检测目标，并列出详细信息。 clean 模式： 直接删除 [check]")

    parser.add_argument('--input_dir', '-i', type=str, default=None,
                        help=f"输入路径. 默认不写，使用 project_type 对应的固定路径。也可以通过input_dir 自己指定路径 [None]\n几种数据类型对应的默认路径为:\n{project_dir_db_string}")

    parser.add_argument('--project', '-p', type=str, default=None,
                        help="需要清理的项目名称，多个项目用逗号分隔， 例如： 21B0206A_WES,21AATI003_WES")

    parser.add_argument('--month', type=float, default=1,
                        help="查找目录下时间超过month个月的目录进行清理 [1]")

    args = parser.parse_args()

    return args


def get_file_over_time(dir, month, type, show_info=True, time_return=False, sort_time=False):
    """[获取目录下，最后一次修改时间超过month个月的文件/目录]

    Args:
        dir ([str]): [目录]
        month ([number]): [时间长度，单位： 月。当前脚本按照一个月31天来计算。]
        type ([str]): [提取类型， file/dir/all  三种文件、目录、全部]
        show_info ([boole]): [是否屏幕打印提取的文件信息] boole
        time_return ([boole]): [返回的文件列表中，是否添加日期信息]
        sort_time ([boole]): [文件按照时间排序。默认 ASCII 排序]
    """    
    time_now = time.time()
    time_cutoff = month * 31 * 24 * 60 * 60
    file_list = []
    for file in os.listdir(dir):
        file_path = os.path.join(dir, file)
        if type == 'file' and os.path.isdir(file_path):
            continue
        if type == 'dir' and os.path.isfile(file_path):
            continue
        # 时间
        time_file = os.path.getmtime(file_path)
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time_file))
        if time_now - time_file < time_cutoff:
            continue

        file_list.append((file, time_str))

    if show_info:
        for file, time_str in sorted(file_list, key=lambda x:(x[1], x[0])):
            print(f"[find] {time_str}  {file}")
    if sort_time:
        file_list = sorted(file_list, key=lambda x:(x[1], x[0]))

    file_final = [x[0] for x in file_list]
    if time_return:
        file_final = file_list
    return file_final


def check_path_status(path_list):
    """[检查列表中的路径的日期、大小等信息，可以是文件、目录]
    Args:
        path_list ([list]): [每一个元素对应一个文件、目录的路径]]
    """    
    # 每个文件、目录向前
    result = []
    # 总文件大小
    size_total = 0
    # 路径没有丢失的文件、目录个数
    path_exists = 0
    for path in path_list:
        path_type = 'lost'
        if os.path.isdir(path):
            path_type = 'dir'
        elif os.path.isfile(path):
            path_type = 'file'
        
        time_info = ''
        size = 0
        if path_type == 'file':
            time_info = time.strftime('%Y-%m-%d %H:%M', time.localtime(os.path.getmtime(path)))
            size = round(os.path.getsize(path) / 1024 / 1024 / 1024, 4)
            path_exists += 1
        if path_type == 'dir':
            # 获取每一个子目录下的文件大小
            for path_sub, dirs, files in os.walk(path):
                for file in files:
                    fp = os.path.join(path_sub, file)
                    if os.path.exists(fp):
                        size += os.path.getsize(fp)
            size =  round(size / 1024 / 1024 / 1024, 4)
            path_exists += 1
        result.append((path_type, time_info, size, path))
        size_total += size
    return result, round(size_total, 4), path_exists


if __name__ == '__main__':
    args = set_and_parse_args()

    print(f"项目类型：{args.project_type}")
    print(f"项目选择模式：{args.strategy}")
    print(f"项目处理模式：{args.mode}")
    log_dir = './log'
    time_str = time.strftime('%Y_%m_%d_%H_%M_%S')
    log_file = f"{log_dir}/{args.project_type}.{time_str}.{args.mode}.log"
    
    # 路径
    if args.input_dir == None:
        args.input_dir = project_dir_db[args.project_type]

    print("(1) 确认项目 \n")
    projects_lost = []
    projects_ok = []
    # 指定项目
    if args.strategy == 'project':
        for project in args.project.split(','):
            project_dir = f"{args.input_dir}/{project}"
            if not os.path.isdir(project_dir):
                projects_lost.append(project_dir)
            else:
                projects_ok.append(project_dir)
    # 时间
    elif args.strategy == 'time':
        projects_name = get_file_over_time(args.input_dir, args.month, 'dir', show_info=True, time_return=False)
        projects_ok = [ f"{args.input_dir}/{project}" for project in projects_name]
    # mysql 已结束项目
    elif args.strategy == 'mysql':
        # 已结束
        finished = {}
        log_mysql_file = f"{log_dir}/{args.project_type}.{time_str}.mysql.{args.mode}.log"
        os.system(f'mysql -u bioinfo -pbioinfo -h 192.168.0.14 -e "select number,status from project where status > 1.5" newweb2018 > {log_mysql_file}')
        # with os.popen(f'mysql -u bioinfo -pbioinfo -h 192.168.0.14 -e "select number,status from project where status > 1.5" newweb2018', encoding='unicode_escape') as f:
        with open(log_mysql_file, 'r', encoding='unicode_escape') as f:
            f.readline()
            for line in f:
                project, status = line.rstrip().split('\t')
                finished[project] = 1
        # 检测 1个月之前的，而且项目已完成的
        for project, time_info in  get_file_over_time(args.input_dir, args.month, 'dir', show_info=False, time_return=True, sort_time=True):
            tmps = project.split('_')
            if tmps[0] in finished:
                print(f"[find] {time_info}  {project}")
                projects_ok.append(f"{args.input_dir}/{project}")

    if len(projects_lost) > 0:
        print(f"[Error] 部分指定项目没有找到, 请仔细核查：{projects_lost}")
        os.sys.exit(1)
    if len(projects_ok) == 0:
        print("没有发现任何满足条件的待删除项目")
        os.sys.exit(1)

    
    # (2) 确认删除列表
    print("\n(2) 确认项目删除列表（检查目录下所有以backup/output开头的目录）")
    clean_list = []
    count_project = 0
    count_not_lost = 0
    size_total = 0
    for project_dir in projects_ok:
        count_project += 1
        print(f"\r        {count_project}/{len(projects_ok)} {project_dir}           ", end='')
        path_tmp = []
        for file in os.listdir(project_dir):
            dir = f"{project_dir}/{file}"
            # WES 项目
            if args.project_type == 'wes':
                if re.match('backup', file) and os.path.isdir(dir):
                    path_tmp.append(f"{dir}/vcf/chr_split")
                    path_tmp.append(f"{dir}/library/homhit_tmp_dir")
                    path_tmp.append(f"{dir}/library/str_tmp_dir")
                    path_tmp.append(f"{dir}/library/vcf2annovar")
                    path_tmp.append(f"{dir}/sv/annotation/breakPointSeq.fa.blast")
                if re.match('output', file) and os.path.isdir(dir):
                    for sample in os.listdir(dir):
                        if os.path.isdir(f"{dir}/{sample}"):
                            path_tmp.append(f"{dir}/{sample}/{sample}_mt_sort.bam")
                            path_tmp.append(f"{dir}/{sample}/{sample}_mt_sort.bam.bai")
                            path_tmp.append(f"{dir}/{sample}/{sample}_mt_sort.rm_dup_flag.bam")
                            path_tmp.append(f"{dir}/{sample}/{sample}_sort_dup.bai")
                            path_tmp.append(f"{dir}/{sample}/{sample}_sort_dup.bam")
                            path_tmp.append(f"{dir}/{sample}/{sample}_sort_dup_realign.bai")
                            path_tmp.append(f"{dir}/{sample}/{sample}_sort_dup_realign.bam")
                            path_tmp.append(f"{dir}/{sample}/{sample}_sort_dup_realign.intervals")
                            path_tmp.append(f"{dir}/{sample}/{sample}_sort_dup_realign_recalibrator.table")
                            path_tmp.append(f"{dir}/{sample}/{sample}_sort.rm_dup_flag.bam")
            # WGS 项目
            elif args.project_type == 'wgs':
                if re.match('backup', file) and os.path.isdir(dir):
                    path_tmp.append(f"{dir}/vcf/chr_split")
                    path_tmp.append(f"{dir}/library/homhit_tmp_dir")
                    path_tmp.append(f"{dir}/library/str_tmp_dir")
                    path_tmp.append(f"{dir}/library/vcf2annovar")
                    path_tmp.append(f"{dir}/sv/annotation/breakPointSeq.fa.blast")
            # Target 项目(PCR/芯片捕获)
            elif args.project_type == 'target':
                if re.match('output', file) and os.path.isdir(dir):
                    path_tmp.append(dir)
            # methyltarget 项目
            elif args.project_type == 'methyltarget':
                if re.match('output', file) and os.path.isdir(dir):
                    path_tmp.append(dir)
            # WTS 项目
            elif args.project_type == 'wts':
                if re.match('backup', file) and os.path.isdir(dir):
                    # RNA
                    path_tmp.extend(glob.glob(f"{dir}/rna/qc/*_final_R*.fastq.gz"))
                    # M6A
                    path_tmp.extend(glob.glob(f"{dir}/qc/*_final_R*.fastq.gz"))
                if re.match('tmp', file) and os.path.isdir(dir):
                    path_tmp.append(dir)
            # MGS 项目
            elif args.project_type == 'mgs':
                if re.match('backup', file) and os.path.isdir(dir):
                    path_tmp.append(dir)
        # 检查每个文件的情况
        path_tmp_check, size, path_not_lost = check_path_status(path_tmp)
        count_not_lost += path_not_lost
        size_total += size
        if path_not_lost == 0:
            print(f"    目录下没有检测到任何可清理的数据")
        else:
            print(f'    目录下检测到 {path_not_lost} 个对象，共 {size}G 数据待删除')
        clean_list.extend(path_tmp_check)
    print(f'\n        共检测到 {count_not_lost} 个待删除的对象，预计可清理出 {round(size_total, 4)}G 空间')
    
    # (3) 清理
    if args.mode == 'check':
        print("(3) 保存待清理日志")
    elif args.mode == 'clean':
        print("(3) 开始清理")
    fh = open(log_file, 'w')
    fh.write('type\ttime\tsize\tfile_dir\n')
    path_count = 0
    for type, time_info, size, file_dir in clean_list:
        path_count += 1
        fh.write(f'{type}\t{time_info}\t{size}\t{file_dir}\n')
        if args.mode  == 'clean':
            print(f"\r清理进度 {path_count}/{len(clean_list)}", end='')
            if type != 'lost':
                os.system(f"rm -rf {file_dir}")
    fh.close()
    print(f"\n[complete] 文件日志 {log_file}")


