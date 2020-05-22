#!/home/genesky/software/python/3.6.7/bin/python3

import plotly.graph_objects as go
import pandas as pd
import operator
import argparse 


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description="bar plot with interactive output") 
     
    parser.add_argument('--input', '-i', type = str, required=True, help = "barplot矩阵，每一列对应一个柱子，每一行对应一个分组。 第一列必须是分组名称，第一行必须是柱子横坐标名称")
    parser.add_argument('--output', '-o', type = str, required=True, help = "barplot绘图输出文件，例如： ./barplot.html")
    parser.add_argument('--barmode', '-m', type = str, default='stack', choices =['stack', 'group'], help = "barplot 绘图模式，支持 stack 和 group 两种模式，默认: stack")
    
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = set_and_parse_args()

    # 读入pc数据，仅保留前两列
    print('[process] read data')
    bar_data = pd.read_table(args.input, sep='\t', header = 0, index_col = 0)
    
    print("[process] plot")
    fig = go.Figure()

    for index in bar_data.index:
        print("          %s" % index)
        fig.add_trace(go.Bar(
            x=bar_data.columns.values.tolist(),
            y=bar_data.loc[[index],:].values.tolist()[0],
            name=index,
        ))
    fig.update_layout(barmode=args.barmode)
    # 输出
    fig.write_html(args.output)
