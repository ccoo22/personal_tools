#!/home/genesky/software/python/3.6.7/bin/python3
import plotly.graph_objects as go
import igraph 
import pandas as pd
import argparse 


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description="graph plot interactive(注：没有legend)") 
     
    parser.add_argument('--edge', '-e', type = str, required=True, help = "边信息，第一、二列是两个顶点, 需要含有表头")
    parser.add_argument('--node', '-n', type = str, required=True, help = "顶点信息，第一列是顶点名称")
    parser.add_argument('--output', '-o', type = str, required=True, help = "绘图输出文件，例如： ./graph.html")
    parser.add_argument('--layout', type = str, default='kamada_kawai', help = "网络图布局设计， 默认： kamada_kawai ", choices=["circle", "drl", "fruchterman_reingold", "fruchterman_reingold_3d", "grid_fruchterman_reingold", "kamada_kawai", "kamada_kawai_3d", "lgl", "random", "random_3d", "reingold_tilford", "reingold_tilford_circular", "sphere"])
    parser.add_argument('--vertex_color', type = str, default=None, help = "顶点颜色列名，node文件，该列数据只能是 red/blue/yellow等标准颜色名称，或者是16进制 #0088ff 类型的数据")
    parser.add_argument('--vertex_size', type = str, default=None, help = "顶点大小列名，node文件， 该列数据只能是数字类型。 数值 5-10基本够用")
    parser.add_argument('--vertex_shape', type = str, default=None, help = "顶点形状列名，node文件， 该列数据只能是： circle, square, diamond, cross, triangle-up, triangle-down, pentagon, hexagon, star, asterisk ......")
    parser.add_argument('--edge_color', type = str, default='black', help = "边的颜色，默认：black, 只能是：red/blue/yellow等标准颜色名称，或者是16进制 #0088ff 类型的数据")
    parser.add_argument('--edge_width', type = int, default=1, help = "边的宽度，默认：1 ")
    
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

    # 检查node数据是否重复
    if node_data.iloc[:,0].duplicated().sum() > 0:
        print("[warning] node文件存在重复node id, 脚本做去重处理，优先保留行id靠前的node id")
        dup = node_data.iloc[:,0].duplicated()
        node_data = node_data.loc[dup == False,:].copy()
        
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
    lay = g.layout(args.layout)   # 记录了每一个顶点布局后的x/y坐标轴，[(x0,y0), (x1,y1),...] 存储


    # 位置信息提取，便于html绘图

    position = {k: lay[k] for k in range(node_count)}  # 记录每一个顶点的x/y轴坐标， dict存储，key = index
    Y = [lay[k][1] for k in range(node_count)]  # 获取素有Y轴坐标
    M = max(Y)  # 最大的Y轴

    es = igraph.EdgeSeq(g) # sequence of edges
    E = [e.tuple for e in g.es] # list of edges, 建议一这段代码改为  g.get_edgelist(), 返回边信息

    # 取得每一个边两侧顶点对应的X坐标，Y坐标
    L = len(position)
    Xn = [position[k][0] for k in range(L)]
    Yn = [2*M-position[k][1] for k in range(L)] # 重新定义Y坐标轴
    Xe = []
    Ye = []
    for edge in E:
        Xe += [position[edge[0]][0], position[edge[1]][0], None]
        Ye += [2*M-position[edge[0]][1], 2*M-position[edge[1]][1], None]

    labels = node_data.iloc[:,0].to_list()


    # 绘图
    fig = go.Figure()

    # 添加线条
    # edge_color = None
    # if args.edge_color != None:
    #     edge_color = edge_data.loc[:, args.edge_color]

    # edge_width = None
    # if args.edge_width != None:
    #     edge_width = edge_data.loc[:, args.edge_width]

    fig.add_trace(go.Scatter(x=Xe,
                   y=Ye,
                   mode='lines',
                   name='lines',
                   line=dict(color=args.edge_color, width=args.edge_width),
                   hoverinfo='none',
                   showlegend=None
                   ))

    # 添加顶点/属性
    vertex_size = None
    if args.vertex_size != None:
        vertex_size = node_data.loc[:, args.vertex_size]

    vertex_shape = None
    if args.vertex_shape != None:
        vertex_shape = node_data.loc[:, args.vertex_shape]

    vertex_color = None
    if args.vertex_color != None:
        vertex_color = node_data.loc[:, args.vertex_color]

    fig.add_trace(go.Scatter(x=Xn,
                  y=Yn,
                  mode='markers',
                  name='scatter',
                  marker=dict(symbol=vertex_shape,
                                size=vertex_size,
                                color=vertex_color,    #'#DB4551',
                                line=dict(color='rgb(50,50,50)', width=1)
                                ),
                  text=labels,
                  hoverinfo='text',
                  opacity=0.8,
                  showlegend=True
                  ))
    # 输出
    fig.write_html(args.output)

