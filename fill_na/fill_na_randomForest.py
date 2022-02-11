#!/home/genesky/software/python/3.9.4/bin/python3
import os
import sys
import argparse
import pandas as pd
import logging
import fill_na_randomForest_test
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(
        description="随机森林回归 数据缺失值填充", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--input', '-i', type=str, required=True,
                        help="需要检测的矩阵，m*n 的矩阵，m个样本、n个特征；或者 m个特征、n个样本。第一行、第一列是特征名、样本名。")

    parser.add_argument('--output', '-o', type=str, required=True,
                        help="输出文件")

    parser.add_argument('--col_sample', action="store_true",
                        help="input m*n 矩阵是 m个特征、n个样本 [False]")

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = set_and_parse_args()

    data_raw = pd.read_table(args.input, sep = '\t', header = 0, index_col = 0)
    if args.col_sample:
        data_raw = data_raw.T

    # 填充
    data_fill = fill_na_randomForest_test.randomForest_fill_na(data_raw)
    # 输出
    if args.col_sample:
        data_fill = data_fill.T
    data_fill.to_csv(path_or_buf=args.output, sep='\t', columns=None, header=True, index=True, index_label=None,mode='w', encoding=None)
    

