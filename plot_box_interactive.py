#!/home/genesky/software/python/3.6.7/bin/python3

import plotly.express as px
import pandas as pd
import operator
import argparse 


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description="box plot with interactive output") 
     
    parser.add_argument('--input', '-i', type = str, required=True, help = "boxplot输入文件，第一作为x轴，第二列作为y轴数据。允许存在第三列，作为分组")
    parser.add_argument('--output', '-o', type = str, required=True, help = "boxplot绘图输出文件，例如： ./boxplot.html")

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = set_and_parse_args()

    # 读入pc数据，仅保留前两列
    print('[process] read data')
    boxplot_data = pd.read_table(args.input, sep='\t', header = 0, index_col = None)

    # 分组颜色设定
    color_name = None
    if len(boxplot_data.columns) >=3:
        color_name = boxplot_data.columns[2]

    # 绘图
    print('[process] plot')
    fig = px.box(boxplot_data, x=boxplot_data.columns[0], y=boxplot_data.columns[1], color=color_name)
    
    # 输出
    fig.write_html(args.output)
