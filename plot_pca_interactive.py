#!/home/genesky/software/python/3.6.7/bin/python3

import plotly.express as px
import pandas as pd
import operator
import argparse 


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description="pca plot with interactive output") 
     
    parser.add_argument('--pca', '-p', type = str, required=True, help = "PCA坐标文件，三列数据：样本名、PC1、PC2, 需要含有表头")
    parser.add_argument('--group', '-g', type = str, required=True, help = "样本分组文件，2-3列数据：样本名、颜色分组、形状分组, 需要含有表头。默认第2列用于颜色，第3列用于形状，可以不输入第3列")
    parser.add_argument('--output', '-o', type = str, required=True, help = "pca绘图输出文件，例如： ./pca.html")
    
    args = parser.parse_args()

    return args


def check_index(data1, data2):
    """检查两个dataframe的index是否一致"""
    if len(data1.index) != len(data2.index):
        print("[Error] 输入的两个文件的样本数量不一致, 请仔细检查")
        exit(0)

    if sum(data1.index == data2.index) != len(data1.index):
        print("[Error] 输入的两个文件的样本名存在不一致的情况， 请仔细检查")
        exit(0)


if __name__ == '__main__':
    args = set_and_parse_args()

    # 读入pc数据，仅保留前两列
    print('[process] read data')
    pc_data = pd.read_table(args.pca, sep='\t', header = 0, index_col = 0)
    pc_data = pc_data.iloc[:, 0:2]
    pc_data.columns = ['PC1', 'PC2']
    pc_data = pc_data.sort_index()
    pc_data['sample'] = pc_data.index
    pc_data = pc_data.loc[:,['sample', 'PC1', 'PC2']]

    # 读入group数据 
    group_data = pd.read_table(args.group, sep='\t', header = 0, index_col = 0)
    group_data = group_data.sort_index()

    # 检查pc/group 的样本是否一致
    print('[process] check sample name')
    check_index(pc_data, group_data) 
    
    # 合并
    print('[process] merge data')
    merge_data = pd.merge(pc_data, group_data, left_index=True, right_index=True)
    
    # 确认颜色/形状分组
    color_name = merge_data.columns[3]
    symbol_name = None
    if len(merge_data.columns) >=5:
        symbol_name = merge_data.columns[4]

    # 绘图
    print('[process] plot')
    fig = px.scatter(data_frame=merge_data, x = 'PC1', y = 'PC2', hover_name='sample', color=color_name, symbol=symbol_name)
 
    # 输出
    fig.write_html(args.output)
