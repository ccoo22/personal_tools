#!/home/genesky/software/python/3.9.4/bin/python3
import os
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn import ensemble 
from sklearn.datasets import load_wine
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
import joblib

def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(
        description="随机森林（分类）\n参数都采用默认的, 遇到实际问题，需要手动调节优化参数\n可以不输入数据，默认调用测试数据集分析。\n", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help="分析结果输出目录")

    parser.add_argument('--model', '-m', type=str, default='train_test', choices=['all', 'train_test', 'train_test_cv'],
                        help="分析模式： \n all: 使用所有样本制作分类器\ntrain_test: 样本3 7 分，7成样本做训练集，3成样本做测试集\ntrain_test_cv:采用交叉验证法评估模型,默认10次")

    parser.add_argument('--feature', '-f', type=str,
                        help="特征文件，m*n的矩阵，m个样本，n个特征，第一列是样本名，第一行是特征名。不要有缺失值，否则流程自动删除带有缺失值的样本。必须是数字，脚本内部自动转换为 np.float32 格式。")

    parser.add_argument('--classify', '-c', type=str,
                        help="样本分类信息，包含三列信息，第一列样本名，第二列样本分组id（数字 0  1 2 构成的分类id）, 第三列分类标签。第一行是表头，名字随意。该文件中的样本必须都被feature文件包含")

    parser.add_argument('--choose_class', type=str,
                        help="从classify中选择指定的几种分类样本用于分析（第二列数字编号），多个特征之间用逗号分隔。默认用所有的样本。")

    args = parser.parse_args()

    return args


def initialize_RandomForestClassifier():
    """[初始化 随机深林分类模型]

    Returns:
        [sklearn.ensemble._forest.RandomForestClassifier]: [sklearn.ensemble._forest.RandomForestClassifier 对象]
    """    
    rfc = ensemble.RandomForestClassifier(n_estimators = 100
            ,criterion='gini'  # 不纯度计算方法 gini / entropy  # 
            ,random_state = None  # 随机性， 可以取任意数值，把随机性固定下来，保证结果的可重复性， 例如： 0
            ,max_depth=None  # 限制树的最大深度，超过设定深度的树枝全部剪掉。剪枝策略。防止过拟合。例如： 3
            ,min_samples_split=2  # 一个节点必须要包含至少min_samples_split个训练样本，这个节点才允许被分枝
            ,min_samples_leaf=1  # 一个节点在分枝后的每个子节点都必须包含至少min_samples_leaf个训练样本否则分枝就不会发生
            ,max_features='auto'  # 限制分枝时考虑的特征个数，超过限制个数的特征都会被舍弃。和max_depth异曲同工,max_features是用来限制高维度数据的过拟合的剪枝参数，但其方法比较暴力,是直接限制可以使用的特征数量而强行使决策树停下的参数，在不知道决策树中的各个特征的重要性的情况下，强行设定这个参数可能会导致模型学习不足。
            ,min_impurity_decrease=0.0  # 限制信息增益的大小，信息增益小于设定数值的分枝不会发生
            ,class_weight=None  # 给不同的样本赋予不同的权重，完成样本标签平衡的参数。样本不平衡是指在一组数据集中，标签的一类天生占有很大的比例，比如说，在银行要判断“一个办了信用卡的人是否会违约”，就是是vs否（1%：99%）的比例
            ,min_weight_fraction_leaf=0.0
            )
    return rfc


if __name__ == '__main__':
    args = set_and_parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # 数据准备
    # 特征矩阵m*n
    data = None
    # 样本分类向量
    target = None
    # data列特征对应的别名向量
    feature_name = None
    # target分类对应的别名向量
    class_name = None
    if args.feature != None and args.classify != None:
        print('使用输入数据建立模型')
        data_raw = pd.read_table(args.feature, sep = '\t', header=0, index_col=0, low_memory=False)
        feature_name = list(data_raw.columns)
        target_raw = pd.read_table(args.classify, sep = '\t', header=0, index_col=0, low_memory=False)
        target_raw.columns = ['target', 'class_name']
        # 删除有缺失值的样本
        target_clean = target_raw[target_raw.isna().sum(axis=1) == 0]
        # 转换数据类型
        target_clean = target_clean.astype({'target': np.int32, 'class_name': 'str'})
        # 排序(很重要)
        target_clean = target_clean.sort_values(by='target', axis=0)
        # 选择指定分组的样本
        if args.choose_class != None:
            choose_class = [int(i) for i in args.choose_class.split(',')]
            target_clean = target_clean[target_clean['target'].isin(choose_class)]
        # data_raw 提取对齐
        data_clean = data_raw.loc[target_clean.index]

        # 再次删除缺失值
        notna_row = data_clean.isna().sum(axis=1) == 0
        data_clean = data_clean[notna_row]
        target_clean = target_clean[notna_row]

        # 最终版整理出来的数据
        data = data_clean.to_numpy()
        target = list(target_clean['target'])
        class_name = list(target_clean['class_name'].drop_duplicates())
    else:
        print("使用系统自带测试数据集建立模型")
        wine = load_wine()
        data, target = wine.data, wine.target
        feature_name = wine.feature_names
        class_name = list(wine.target_names)
    print(f"分析模式：{args.model}")

    if args.model == 'all':
        rfc = initialize_RandomForestClassifier()
        # 训练
        rfc = rfc.fit(data, target)
        # 导入测试集，计算准确度
        score = rfc.score(data, target)
        # 预测
        # predict_test = rfc.predict(Xtest)

        # 输出score
        with open(os.path.join(args.output_dir, 'score.txt'), 'w') as fh:
            fh.write(str(score))
        print(f"score={score}")

        # 每个特征的重要性
        feature_importance = pd.DataFrame(zip(feature_name, rfc.feature_importances_), columns=['feature_name', 'feature_importances'])
        feature_importance.to_csv(os.path.join(args.output_dir, 'feature_importance.txt'), sep='\t',header=True, index=False)

        # 模型保存，方便以后直接使用它来预测
        joblib_file = os.path.join(args.output_dir, 'model.pkl')
        joblib.dump(rfc, joblib_file)
        # 模型加载
        # joblib_model = joblib.load(joblib_file)

    # train_test 模式
    if args.model == 'train_test':
        # 训练集、测试集拆分，同时把样本顺序打乱
        Xtrain, Xtest, Ytrain, Ytest = train_test_split(data, target, test_size=0.3, shuffle =True)
        # 模型初始化
        rfc = initialize_RandomForestClassifier()
        # 训练
        rfc = rfc.fit(Xtrain, Ytrain)
        # 导入测试集，计算准确度
        score = rfc.score(Xtest, Ytest)
        # 预测
        # predict_test = rfc.predict(Xtest)

        # 输出score
        with open(os.path.join(args.output_dir, 'score.txt'), 'w') as fh:
            fh.write(str(score))
        print(f"score={score}")

        # 每个特征的重要性
        feature_importance = pd.DataFrame(zip(feature_name, rfc.feature_importances_), columns=['feature_name', 'feature_importances'])
        feature_importance.to_csv(os.path.join(args.output_dir, 'feature_importance.txt'), sep='\t',header=True, index=False)
        
        # 模型保存，方便以后直接使用它来预测
        joblib_file = os.path.join(args.output_dir, 'model.pkl')
        joblib.dump(rfc, joblib_file)
        # 模型加载
        # joblib_model = joblib.load(joblib_file)

    # train_test cv 模式
    if args.model == 'train_test_cv':
        rfc = initialize_RandomForestClassifier()
        cv = 10
        scores = cross_val_score(rfc, data, target, cv=cv)
        score_str = '\n'.join([str(x) for x in scores])
        print(f"交叉验证 {cv} 次， 每次的准确率为： \n{score_str}")
        score = scores.mean()
        # 输出score
        with open(os.path.join(args.output_dir, 'score.txt'), 'w') as fh:
            fh.write(str(score))
        print(f"score={score}")

