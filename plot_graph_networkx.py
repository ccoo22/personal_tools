#!/home/genesky/software/python/3.6.7/bin/python3

import networkx as nx 
import matplotlib.pyplot as plt
import pandas as pd
import argparse 


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description="graph plot ") 
     
    parser.add_argument('--edge', '-e', type = str, required=True, help = "边信息，第一、二列是两个顶点, 需要含有表头")
    parser.add_argument('--node', '-n', type = str, required=True, help = "顶点信息，第一列是顶点名称")
    parser.add_argument('--output', '-o', type = str, required=True, help = "绘图输出文件，例如： ./graph.pdf")
    parser.add_argument('--pdf_width', type = int, default=12.8, help = "pdf宽度， 默认： 12.8, 单位：inch")
    parser.add_argument('--pdf_height', type = int, default=9.6, help = "pdf高度， 默认： 9.6, 单位：inch, 可以通过增大该参数，减小顶点的大小")
    parser.add_argument('--layout', type = str, default='spring', help = "网络图布局设计， 默认： spring ", choices=["circular", "kamada_kawai", "planar", "random", "shell", "spring", "spiral"])
    parser.add_argument('--layout_sprint_k', type = float, default=None, help = "sprint layout 的参数k设置，默认 = 1/sqrt(n)， 值越大，点之间的距离越远")
    parser.add_argument('--show_label', action = 'store_true', help = "显示顶点名称")
    parser.add_argument('--font_size', type = int, default=5, help = "设置名称大小， 默认： 5")
    parser.add_argument('--node_color', type = str, default=None, help = "顶点颜色列名，node文件，该列数据只能是 red/blue等标准颜色名称，或者是16进制 #0088ff 类型的数据， 默认： #1f78b4")
    parser.add_argument('--node_size', type = str, default=None, help = "顶点大小列名，node文件， 该列数据只能是数字类型。 默认：300")
    parser.add_argument('--node_shape', type = str, default=None, help = "顶点形状列名，node文件， 该列数据只能是： o/v/^/</>/s/D/P/+ 等， 参考：https://matplotlib.org/3.2.1/api/markers_api.html。 默认： o ")
    parser.add_argument('--edge_color', type = str, default=None, help = "边的颜色列名，edge文件， 该列数据只能是：red/blue等标准颜色名称，或者是16进制 #0088ff 类型的数据, 默认： k ")
    parser.add_argument('--edge_width', type = str, default=1, help = "边的宽度, 默认： 1")
    
    args = parser.parse_args()

    return args


def check_graph(edge, node):
    """检查node是否包含了edge中的所有顶点"""

    # 取出edge中的顶点、去重
    edge_node = edge_data.iloc[:,0].to_list()
    edge_node.extend(edge_data.iloc[:,1].to_list())
    edge_node = set(edge_node)

    node_list = set(node.iloc[:,0].to_list())

    # 
    if node_list < edge_node:
        print("[Error] edge文件中包含未在node文件中声明的顶点")
        exit(0)        


if __name__ == '__main__':
    args = set_and_parse_args()

    # 读入pc数据，仅保留前两列
    print('[process] read data')
    edge_data = pd.read_table(args.edge, sep='\t', header = 0, index_col = None)
    node_data = pd.read_table(args.node, sep='\t', header = 0, index_col = None)

    edge_count = edge_data.shape[0]
    node_count = node_data.shape[0]

    # 检查是edge/node 异常
    print('[process] check node')
    check_graph(edge_data, node_data)


    # 创建
    print('[process] creat graph')
    G = nx.Graph()  # 无向网络图
    edge_list = [ (edge_data.iat[row, 0], edge_data.iat[row, 1]) for row in range(edge_count) ] 
    G.add_edges_from(edge_list)

    # layout
    pos=None
    fig = plt.figure(figsize=(args.pdf_width, args.pdf_height))
    # if args.layout == 'bipartite':
    #     top = nx.bipartite.sets(G)[0]
    #     pos = nx.bipartite_layout(G, top)

    if args.layout == 'circular':
        pos = nx.circular_layout(G)

    if args.layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)

    if args.layout == 'planar':
        pos = nx.planar_layout(G)

    if args.layout == 'random':
        pos = nx.random_layout(G)

    # if args.layout == 'rescale':
    #     pos = nx.rescale_layout(G)

    if args.layout == 'shell':
        pos = nx.shell_layout(G)

    if args.layout == 'spring':
        pos = nx.spring_layout(G, args.layout_sprint_k)

    if args.layout == 'spectral':
        pos = nx.spectral_layout(G)

    if args.layout == 'spiral':
        pos = nx.spiral_layout(G)

    # 绘图
    print('[process] plot')
    # 显示设置
    node_size = 300
    if args.node_size:
        size_map = { node_data.iat[row, 0] : node_data.at[row, args.node_size] for row in range(node_count) }
        node_size = [ size_map[node] for node in G.nodes]

    node_color = '#1f78b4'
    if args.node_color:
        color_map = { node_data.iat[row, 0] : node_data.at[row, args.node_color] for row in range(node_count) }
        node_color = [ color_map[node] for node in G.nodes]

    node_shape = 'o'
    if args.node_shape:
        shape_map = { node_data.iat[row, 0] : node_data.at[row, args.node_shape] for row in range(node_count) }
        node_shape = [ shape_map[node] for node in G.nodes]

    edge_color='k'
    if args.edge_color != None:
        edge_color_map = { (edge_data.iat[row, 0], edge_data.iat[row, 1]) : edge_data.at[row, args.edge_color] for row in range(edge_count) }
        edge_color = [ edge_color_map[edge] for edge in G.edges]


    nx.draw(G, with_labels=args.show_label, pos=pos, font_size=args.font_size, node_size=node_size, node_color=node_color, node_shape=node_shape, width=args.edge_width)
    fig.savefig(args.output)
 
