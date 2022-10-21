#!/home/genesky/software/python/3.8.12/bin/python3
import argparse
import xlrd
import os

def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(
        description="xlsx 表格转txt", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--input', '-i', type=str, required=True,
                        help="excel 输入文件路径")

    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help="结果输出目录")

    parser.add_argument('--sheet_name', '-s', type=str,
                        help="指定要输出的excel sheet名称，默认输出input中的所有sheet")
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = set_and_parse_args()
    eob = xlrd.open_workbook(args.input)
    sheet_names = eob.sheet_names()
    if args.sheet_name:
        if args.sheet_name in sheet_names:
            sheet_names = args.sheet_name
        else:
            print("[error] 你要求的sheet_name 不存在，请确认是否写错了！")
    if(not os.path.exists(args.output_dir)):
        os.makedirs(args.output_dir)
    for sheet_name in sheet_names:
        output_file = f"{args.output_dir}/{sheet_name}.txt"
        print(f'[process] {sheet_name} -> {output_file}')
        sob =eob.sheet_by_name(sheet_name)
        with open(output_file, 'w') as fh:
            for row in range(sob.nrows):
                row_value = sob.row_values(row)
                fh.write('\t'.join([str(x) for x in row_value]) + "\n")
