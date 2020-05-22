#!/home/genesky/software/python/3.6.7/bin/python3

import plotly.express as px
import pandas as pd
import operator
import argparse 


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description="pca plot with interactive output") 
     
    parser.add_argument('--pca', '-p', type = str, required=True, help = "PCA坐标+其他信息列文件，第一列必须是样本名, 需要含有表头")
    parser.add_argument('--output', '-o', type = str, required=True, help = "pca绘图输出文件，例如： ./pca.html")

    parser.add_argument('--pc1', type = str, required=True, help = "pc1坐标表头名称")
    parser.add_argument('--pc2', type = str, required=True, help = "pc2坐标表头名称")
    parser.add_argument('--color', type = str, required=False, default=None, help = "定义颜色表头名称")
    parser.add_argument('--symbol', type = str, required=False, default=None, help = "定义形状表头名称")
    parser.add_argument('--size', type = str, required=False, default=None, help = "定义点的大小表头名称")
    parser.add_argument('--size_max', type = int, required=False, default=None, help = "定义点的最大大小")


    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = set_and_parse_args()

    # 读入pc数据，仅保留前两列
    print('[process] read data')
    pc_data = pd.read_table(args.pca, sep='\t', header = 0, index_col = 0)
    pc_data['sample'] = pc_data.index

    # 检查输入的列名是否有问题
    if args.pc1 not in  pc_data.columns:
        print("[Error] 输入的文件中不包含列： %s" % args.pc1)
        exit(0)

    if args.pc2 not in  pc_data.columns:
        print("[Error] 输入的文件中不包含列： %s" % args.pc2)
        exit(0)

    if args.color != None and args.color not in  pc_data.columns:
        print("[Error] 输入的文件中不包含列： %s" % args.color)
        exit(0)

    if args.symbol != None and args.symbol not in  pc_data.columns:
        print("[Error] 输入的文件中不包含列： %s" % args.symbol)
        exit(0)

    if args.size != None and args.size not in  pc_data.columns:
        print("[Error] 输入的文件中不包含列： %s" % args.size)
        exit(0)

    # 绘图
    print('[process] plot')
    fig = px.scatter(data_frame=pc_data, x = args.pc1, y = args.pc2, hover_name='sample', color=args.color, symbol=args.symbol, size=args.size, size_max=args.size_max)
    
    # 输出
    fig.write_html(args.output)
