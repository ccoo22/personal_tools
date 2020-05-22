#!/home/genesky/software/python/3.6.7/bin/python3

import igraph 
import pandas as pd
import argparse 


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description="graph plot ") 
     
    parser.add_argument('--edge', '-e', type = str, required=True, help = "边信息，第一、二列是两个顶点, 需要含有表头")
    parser.add_argument('--node', '-n', type = str, required=True, help = "顶点信息，第一列是顶点名称")
    parser.add_argument('--output', '-o', type = str, required=True, help = "绘图输出文件，例如： ./graph.pdf")
    parser.add_argument('--pdf_width', type = int, default=600, help = "pdf宽度， 默认： 600， 单位：像素")
    parser.add_argument('--pdf_height', type = int, default=600, help = "pdf高度， 默认： 600， 单位：像素, 可以通过增大该参数，减小顶点的大小")
    parser.add_argument('--margin', type = int, default=20, help = "边框距离， 默认： 20")
    parser.add_argument('--layout', type = str, default='kamada_kawai', help = "网络图布局设计， 默认： kamada_kawai ", choices=["circle", "drl", "fruchterman_reingold", "fruchterman_reingold_3d", "grid_fruchterman_reingold", "kamada_kawai", "kamada_kawai_3d", "lgl", "random", "random_3d", "reingold_tilford", "reingold_tilford_circular", "sphere"])
    parser.add_argument('--show_label', action = 'store_true', help = "显示顶点名称")
    parser.add_argument('--label_size', type = int, default=5, help = "设置名称大小， 默认： 5")
    parser.add_argument('--vertex_color', type = str, default=None, help = "顶点颜色列名，node文件，该列数据只能是 red/blue等标准颜色名称，或者是16进制 #0088ff 类型的数据")
    parser.add_argument('--vertex_size', type = str, default=None, help = "顶点大小列名，node文件， 该列数据只能是数字类型。 数值 5-10基本够用")
    parser.add_argument('--vertex_shape', type = str, default=None, help = "顶点形状列名，node文件， 该列数据只能是： rectangle, circle, hidden, triangle-up, triangle-down")
    parser.add_argument('--edge_color', type = str, default=None, help = "边的颜色列名，edge文件， 该列数据只能是：red/blue等标准颜色名称，或者是16进制 #0088ff 类型的数据")
    parser.add_argument('--edge_width', type = str, default=None, help = "边的宽度列名，edge文件， 单位是像素, 数值宽度 1 基本够用")
    
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

    node_count = node_data.shape[0]

    # 检查是edge/node 异常
    print('[process] check node')
    check_graph(edge_data, node_data)

    # 基于igraph的特性，进行特殊处理
    # 顶点名称 -> ID 映射关系
    node_id_dict = { node_data.iat[id,0]: id for id in range(0, node_count)}  

    # 边映射关系生成
    edge_igraph = [ (node_id_dict[edge_data.iat[row, 0]], node_id_dict[edge_data.iat[row, 1]]) for row in range(0, edge_data.shape[0])]

    print('[process] plot')
    g = igraph.Graph()
    g.add_vertices(node_count)  # 添加顶点
    g.add_edges(edge_igraph)

    # 显示设置
    visual_style = {}
    visual_style["layout"] = g.layout(args.layout) 
    visual_style["vertex_label_size"] = args.label_size
    visual_style["bbox"] = (args.pdf_width, args.pdf_height)
    visual_style["margin"] = args.margin

    if args.show_label:
        visual_style["vertex_label"] = node_data.iloc[:,0].to_list()

    if args.vertex_color != None:
        visual_style["vertex_color"] = node_data.loc[:, args.vertex_color]

    if args.vertex_size != None:
        visual_style["vertex_size"] = node_data.loc[:, args.vertex_size]

    if args.vertex_shape != None:
        visual_style["vertex_shape"] = node_data.loc[:, args.vertex_shape]

    if args.edge_color != None:
        visual_style["edge_color"] = edge_data.loc[:, args.edge_color]

    if args.edge_width != None:
        visual_style["edge_width"] = edge_data.loc[:, args.edge_width]
    
    # 绘图
    igraph.plot(g, args.output, **visual_style)
 
