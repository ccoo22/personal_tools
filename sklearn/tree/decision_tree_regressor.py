#!/home/genesky/software/python/3.9.4/bin/python3
import os
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn import tree 
from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt
import graphviz
import joblib

def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(
        description="决策树（回归）\n参数都采用默认的, 遇到实际问题，需要手动调节优化参数\n可以不输入数据，默认调用测试数据集分析。\n", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help="分析结果输出目录")

    parser.add_argument('--model', '-m', type=str, default='train_test', choices=['all', 'train_test', 'train_test_cv'],
                        help="分析模式： \n all: 使用所有样本制作分类器\ntrain_test: 样本3 7 分，7成样本做训练集，3成样本做测试集\ntrain_test_cv:采用交叉验证法评估模型,默认10次")

    parser.add_argument('--feature', '-f', type=str,
                        help="特征文件，m*n的矩阵，m个样本，n个特征，第一列是样本名，第一行是特征名。不要有缺失值，否则流程自动删除带有缺失值的样本。必须是数字，脚本内部自动转换为 np.float32 格式。")

    parser.add_argument('--predict', '-p', type=str,
                        help="样本要预测的信息，包含两列信息，第一列样本名，第二列预测数据（连续型数字）。第一行是表头，名字随意。该文件中的样本必须都被feature文件包含")

    args = parser.parse_args()

    return args


def initialize_DecisionTreeRegressor():
    """[初始化 决策树回归模型]

    Returns:
        [sklearn.tree._classes.DecisionTreeRegressor]: [sklearn.tree._classes.DecisionTreeRegressor 对象]
    """    
    reg = tree.DecisionTreeRegressor(criterion='mse'  # "mse"使用均方误差mean squared error(MSE), “friedman_mse”使用费尔德曼均方误差, "mae"使用绝对平均误差MAE（mean absolute error）
            ,random_state = None  # 随机性， 可以取任意数值，把随机性固定下来，保证结果的可重复性， 例如： 0
            ,splitter='best'  # 分支时，节点特征的选择方式。 best / random。 best  选择效果最好的特征、随机选择效果最好的特征。也会增加随机性
            ,max_depth=None  # 限制树的最大深度，超过设定深度的树枝全部剪掉。剪枝策略。防止过拟合。例如： 3
            ,min_samples_split=2  # 一个节点必须要包含至少min_samples_split个训练样本，这个节点才允许被分枝
            ,min_samples_leaf=1  # 一个节点在分枝后的每个子节点都必须包含至少min_samples_leaf个训练样本否则分枝就不会发生
            ,max_features=None  # 限制分枝时考虑的特征个数，超过限制个数的特征都会被舍弃。和max_depth异曲同工,max_features是用来限制高维度数据的过拟合的剪枝参数，但其方法比较暴力,是直接限制可以使用的特征数量而强行使决策树停下的参数，在不知道决策树中的各个特征的重要性的情况下，强行设定这个参数可能会导致模型学习不足。
            ,min_impurity_decrease=0.0  # 限制信息增益的大小，信息增益小于设定数值的分枝不会发生
            )
    return reg


if __name__ == '__main__':
    args = set_and_parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # 数据准备
    # 特征矩阵m*n
    data = None
    # 样本连续型表型
    target = None
    # data列特征对应的别名向量
    feature_name = None

    if args.feature != None and args.predict != None:
        print('使用输入数据建立模型')
        data_raw = pd.read_table(args.feature, sep = '\t', header=0, index_col=0, low_memory=False)
        feature_name = list(data_raw.columns)
        target_raw = pd.read_table(args.predict, sep = '\t', header=0, index_col=0, low_memory=False)
        target_raw.columns = ['target']
        # 删除有缺失值的样本
        target_clean = target_raw[target_raw.isna().sum(axis=1) == 0]
        # 转换数据类型
        target_clean = target_clean.astype({'target': np.float64})
        # 排序(很重要)
        target_clean = target_clean.sort_values(by='target', axis=0)

        # data_raw 提取对齐
        data_clean = data_raw.loc[target_clean.index]

        # 再次删除缺失值
        notna_row = data_clean.isna().sum(axis=1) == 0
        data_clean = data_clean[notna_row]
        target_clean = target_clean[notna_row]

        # 最终版整理出来的数据
        data = data_clean.to_numpy()
        target = list(target_clean['target'])
    else:
        print("使用系统自带测试数据集建立模型")
        boston = load_boston()
        data, target = boston.data.copy(), boston.target.copy()
        feature_name = boston.feature_names
        # 排序， 方便绘图
        order = target.argsort()
        data = data[order,:]
        target = target[order]
    print(f"分析模式：{args.model}")

    if args.model == 'all':
        reg = initialize_DecisionTreeRegressor()
        # 训练
        reg = reg.fit(data, target)
        # 导入测试集，计算R2值，越大越好，最大值为1
        score = reg.score(data, target)
        # 预测
        predict_train= reg.predict(data)

        # 输出score
        with open(os.path.join(args.output_dir, 'score.txt'), 'w') as fh:
            fh.write(str(score))
        print(f"score={score}")

        # 每个特征的重要性
        feature_importance = pd.DataFrame(zip(feature_name, reg.feature_importances_), columns=['feature_name', 'feature_importances'])
        feature_importance.to_csv(os.path.join(args.output_dir, 'feature_importance.txt'), sep='\t',header=True, index=False)

        # 对树模型绘图
        dot_data = tree.export_graphviz(reg
            ,out_file = None
            ,feature_names= feature_name
            ,filled=True
            ,rounded=True
            )
        graph = graphviz.Source(dot_data)
        graph.render(filename="tree", directory=args.output_dir, cleanup=True, format='pdf')

        # 散点图
        png_train = os.path.join(args.output_dir, 'train.png')
        plt.figure()
        plt.scatter(x=np.arange(len(target)) + 1, y=target, c="darkorange", label="raw")
        plt.scatter(x=np.arange(len(predict_train)) + 1, y=predict_train, c="darkblue", label="predict")
        plt.legend()
        plt.savefig(png_train, format = 'png')

        # 模型保存，方便以后直接使用它来预测
        joblib_file = os.path.join(args.output_dir, 'model.pkl')
        joblib.dump(reg, joblib_file)
        # 模型加载
        # joblib_model = joblib.load(joblib_file)

    # train_test 模式
    if args.model == 'train_test':
        # 训练集、测试集拆分，同时把样本顺序打乱
        Xtrain, Xtest, Ytrain, Ytest = train_test_split(data, target, test_size=0.3, shuffle=True)
        # 排序，方便绘图
        order1 = Ytrain.argsort()
        Xtrain = Xtrain[order1, :]
        Ytrain = Ytrain[order1]

        order2 = Ytest.argsort()
        Xtest = Xtest[order2, :]
        Ytest = Ytest[order2]


        # 模型初始化
        reg = initialize_DecisionTreeRegressor()
        # 训练
        reg = reg.fit(Xtrain, Ytrain)
        # 导入测试集，计算准确度
        score = reg.score(Xtest, Ytest)
        # 预测
        predict_train = reg.predict(Xtrain)
        predict_test  = reg.predict(Xtest)

        # 输出score
        with open(os.path.join(args.output_dir, 'score.txt'), 'w') as fh:
            fh.write(str(score))
        print(f"score={score}")

        # 每个特征的重要性
        feature_importance = pd.DataFrame(zip(feature_name, reg.feature_importances_), columns=['feature_name', 'feature_importances'])
        feature_importance.to_csv(os.path.join(args.output_dir, 'feature_importance.txt'), sep='\t',header=True, index=False)

        # 对树模型绘图
        dot_data = tree.export_graphviz(reg
            ,out_file = None
            ,feature_names= feature_name
            ,filled=True
            ,rounded=True
            )
        graph = graphviz.Source(dot_data)
        graph.render(filename="tree", directory=args.output_dir, cleanup=True, format='pdf')
        
        # 散点图
        png_train = os.path.join(args.output_dir, 'train.png')
        plt.figure()
        plt.scatter(x=np.arange(len(Ytrain)) + 1, y=Ytrain, c="darkorange", label="raw")
        plt.scatter(x=np.arange(len(predict_train)) + 1, y=predict_train, c="darkblue", label="predict")
        plt.legend()
        plt.savefig(png_train, format = 'png')

        png_test = os.path.join(args.output_dir, 'test.png')
        plt.figure()
        plt.scatter(x=np.arange(len(Ytest)) + 1, y=Ytest, c="darkorange", label="raw")
        plt.scatter(x=np.arange(len(predict_test)) + 1, y=predict_test, c="darkblue", label="predict")
        plt.legend()
        plt.savefig(png_test, format = 'png')

        # 模型保存，方便以后直接使用它来预测
        joblib_file = os.path.join(args.output_dir, 'model.pkl')
        joblib.dump(reg, joblib_file)
        # 模型加载
        # joblib_model = joblib.load(joblib_file)

    # train_test cv 模式
    if args.model == 'train_test_cv':
        reg = initialize_DecisionTreeRegressor()
        cv = 10
        scores = cross_val_score(reg, data, target, cv=cv)
        score_str = '\n'.join([str(x) for x in scores])
        print(f"交叉验证 {cv} 次， 每次的score R2值为： \n{score_str}")
        score = scores.mean()
        # 输出score
        with open(os.path.join(args.output_dir, 'score.txt'), 'w') as fh:
            fh.write(str(score))
        print(f"score={score}")



# reg = tree.DecisionTreeRegressor(criterion='mse'  # "mse"使用均方误差mean squared error(MSE), “friedman_mse”使用费尔德曼均方误差, "mae"使用绝对平均误差MAE（mean absolute error）
#             ,random_state = 0  # 随机性， 可以取任意数值，把随机性固定下来，保证结果的可重复性， 例如： 0
#             ,splitter='best'  # 分支时，节点特征的选择方式。 best / random。 best  选择效果最好的特征、随机选择效果最好的特征。也会增加随机性
#             ,max_depth=3  # 限制树的最大深度，超过设定深度的树枝全部剪掉。剪枝策略。防止过拟合。例如： 3
#             ,min_samples_split=5  # 一个节点必须要包含至少min_samples_split个训练样本，这个节点才允许被分枝
#             ,min_samples_leaf=3  # 一个节点在分枝后的每个子节点都必须包含至少min_samples_leaf个训练样本否则分枝就不会发生
#             ,max_features=None  # 限制分枝时考虑的特征个数，超过限制个数的特征都会被舍弃。和max_depth异曲同工,max_features是用来限制高维度数据的过拟合的剪枝参数，但其方法比较暴力,是直接限制可以使用的特征数量而强行使决策树停下的参数，在不知道决策树中的各个特征的重要性的情况下，强行设定这个参数可能会导致模型学习不足。
#             ,min_impurity_decrease=0.0  # 限制信息增益的大小，信息增益小于设定数值的分枝不会发生
#             )
# Xtrain, Xtest, Ytrain, Ytest = train_test_split(data, target, test_size=0.3, shuffle=False)
# # 训练
# reg = reg.fit(Xtrain, Ytrain)
# # 导入测试集，计算准确度
# score = reg.score(Xtest, Ytest)
# # 预测
# predict_train = reg.predict(Xtrain)
# predict_test  = reg.predict(Xtest)

# plt.figure()
# plt.scatter(x=np.arange(len(Ytrain)) + 1, y=Ytrain, c="darkorange", label="raw")
# plt.scatter(x=np.arange(len(predict_train)) + 1, y=predict_train, c="darkblue", label="predict")
# plt.legend()


# plt.figure()
# plt.scatter(x=np.arange(len(Ytest)) + 1, y=Ytest, c="darkorange", label="raw")
# plt.scatter(x=np.arange(len(predict_test)) + 1, y=predict_test, c="darkblue", label="predict")
# plt.legend()



# rng = np.random.RandomState(1)
# X = np.sort(5 * rng.rand(80,1), axis=0)
# y = np.sin(X).ravel()
# order2 = y.argsort()
# X = X[order2,:]
# y = y[order2]
# y[::5] += 3 * (0.5 - rng.rand(16))
# regr_1 = tree.DecisionTreeRegressor(max_depth=2)
# regr_2 = tree.DecisionTreeRegressor(max_depth=5)
# regr_1.fit(X, y)
# regr_2.fit(X, y)
# predict_y1 = regr_1.predict(X)
# predict_y2 = regr_2.predict(X)
# plt.figure()
# plt.scatter(x=np.arange(len(y)) + 1, y=y, edgecolor="black",c="darkorange", label="data")
# plt.scatter(x=np.arange(len(predict_y1)) + 1, y=predict_y1, color="cornflowerblue",label="max_depth=2", linewidth=2)
# plt.scatter(x=np.arange(len(predict_y2)) + 1, y=predict_y2, color="yellowgreen", label="max_depth=5", linewidth=2)



# X_test = np.arange(0.0, 5.0, 0.01)[:, np.newaxis]
# y_test = np.sin(X_test).ravel()
# order3 = y_test.argsort()
# X_test = X_test[order3,:]
# y_test = y_test[order3]

# y_1 = regr_1.predict(X_test)
# y_2 = regr_2.predict(X_test)
# plt.figure()
# plt.scatter(x=np.arange(len(y_test)) + 1, y=y_test, edgecolor="black",c="darkorange", label="data")
# plt.scatter(x=np.arange(len(y_1)) + 1, y=y_1, color="cornflowerblue",label="max_depth=2", linewidth=2)
# plt.scatter(x=np.arange(len(y_2)) + 1, y=y_2, color="yellowgreen", label="max_depth=5", linewidth=2)


 


