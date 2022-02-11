#!/home/genesky/software/python/3.9.4/bin/python3
import os
import sys
import argparse
import multiprocessing
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestRegressor
import logging
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(
        description="随机森林回归 数据缺失值填充测试。 \n注意：为了保证测试的正常进行。流程会自动把有缺失的样本、特征删除掉", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--input', '-i', type=str, required=True,
                        help="需要检测的矩阵，m*n 的矩阵，m个样本、n个特征；或者 m个特征、n个样本。第一行、第一列是特征名、样本名。")

    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help="输出路径，流程自动创建 ")

    parser.add_argument('--col_sample', action="store_true",
                        help="input m*n 矩阵是 m个特征、n个样本 [False]")

    parser.add_argument('--missing_rate_min', type=float, default=0.01,
                        help="最低缺失率,小数位不要超过2位 [0.01]")

    parser.add_argument('--missing_rate_max', type=float, default=0.7,
                        help="最高缺失率 [0.7]")

    parser.add_argument('--missing_rate_step', type=float, default=0.01,
                        help="缺失率步长 [0.01]")

    parser.add_argument('--thread', type=int, default=10,
                        help="并行线程 [10]")

    args = parser.parse_args()

    return args


def randomForest_fill_na(input, n_estimators=100):
    """[随机森林回归填充缺失值]

    Args:
        input ([pandas.core.frame.DataFrame]): [m*n矩阵，m个样本，n个特征]
        n_estimators ([int]): [随机森林中树的数量]
    """    
    # 随机森林回归填充
    X_missing_reg = input.copy()
    # 清楚行、列索引
    X_missing_reg = pd.DataFrame(X_missing_reg.values)
    sortindex = np.argsort(X_missing_reg.isnull().sum(axis=0)).values
    for i in sortindex:
        # 构建我们的新特征矩阵和新标签
        df = X_missing_reg
        fillc = df.iloc[:,i]
        # 没有缺失/都是缺失
        if fillc.isnull().sum() == 0 or fillc.isnull().sum() == fillc.shape[0]:
            continue 
        df = df.iloc[:,df.columns != i]
        # 在新特征矩阵中，对含有缺失值的列，进行0的填补
        df_0 =SimpleImputer(missing_values=np.nan, strategy='constant',fill_value=0).fit_transform(df)
        # 找出我们的训练集和测试集
        Ytrain = fillc[fillc.notnull()]
        Ytest = fillc[fillc.isnull()]
        Xtrain = df_0[Ytrain.index,:]
        Xtest = df_0[Ytest.index,:]
        # 用随机森林回归来填补缺失值, 100棵树
        rfc = RandomForestRegressor(n_estimators=n_estimators)
        rfc = rfc.fit(Xtrain, Ytrain)
        Ypredict = rfc.predict(Xtest)
        # 将填补好的特征返回到我们的原始的特征矩阵中
        X_missing_reg.loc[X_missing_reg.iloc[:,i].isnull(),i] = Ypredict
    X_missing_reg.index = input.index
    X_missing_reg.columns = input.columns
    return X_missing_reg


def fill_and_check_quality(input, missing_rate, output_dir):
    """[模拟缺失，并填充，检验]

    Args:
        input ([pandas.core.frame.DataFrame]): [不存在缺失的矩阵，m*n，m个样本，n个特征]
        missing_rate ([float]): [设定缺失率，例如 0.01]
        output_dir ([str]): [结果输出目录]
    """    
    n_samples = input.shape[0]
    n_features = input.shape[1]

    # 模拟数据缺失
    rng = np.random.RandomState(0)
    # 缺失数量
    n_missing = int(np.floor(n_samples * n_features * missing_rate))
    # 不放回随机抽样
    missing = rng.choice(n_samples * n_features, size=n_missing, replace=False)  
    # 得到列坐标
    missing_features = [int(i) for i in missing / n_samples] 
    # 得到行坐标
    missing_samples = missing % n_samples  
    # 设置缺失值
    X_missing = input.copy().values  # values 转换为numpy array 格式
    X_missing[missing_samples,missing_features] = np.nan  # 替换为np.nan  只有 numpy array 可以这么填充
    X_missing = pd.DataFrame(X_missing)

    # 随机森林填充
    X_fill_rf = randomForest_fill_na(X_missing)

    # 均值填充
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    X_fill_mean = imp_mean.fit_transform(X_missing)

    # 0填充
    imp_0 = SimpleImputer(missing_values=np.nan, strategy="constant",fill_value=0)
    X_fill_0 = imp_0.fit_transform(X_missing)

    # 实际值与填充值对比
    names = ['Regressor Imputation', 'Mean Imputation', '0 Imputation']

    fig = plt.figure()
    fig.set_size_inches(24, 8)
    count = 0
    value_true_raw = input.values[missing_samples, missing_features]
    mses = []
    for fill in [X_fill_rf.values, X_fill_mean, X_fill_0]:
        count += 1
        # 实际值与填充值线性绘图
        value_fill_raw = fill[missing_samples, missing_features]
        value_true = value_true_raw[value_fill_raw != np.nan]
        value_fill = value_fill_raw[value_fill_raw != np.nan]
        mse = sum((value_true - value_fill) ** 2) / len(value_true)
        mses.append(mse)
        plt.subplot(1, 3, count)
        plt.plot(value_true, value_fill, 'ro')
        plt.plot([0, 1], [0, 1], 'b-')
        plt.text(0, 1, f'MSE={mse}')
        plt.xlabel('value true') 
        plt.ylabel('value fill')
        plt.title(names[count - 1])

    file_prefix = f"missing_rate_{'%.2f' % missing_rate}"
    plt.savefig(os.path.join(output_dir, f"{file_prefix}.png"), format = 'png')

    fh = open(os.path.join(output_dir, f"{file_prefix}.mse.txt"), 'w')
    for name, mse in zip(names, mses):
        fh.write(name + '\t' + str(mse) + '\n')
    fh.close()



if __name__ == '__main__':
    args = set_and_parse_args()
    missing_rates = [ i / 100 for i in range(int(args.missing_rate_min * 100), int(args.missing_rate_max * 100), int(args.missing_rate_step * 100))]
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    log.info("[process 1] 数据读入、清理")

    data_raw = pd.read_table(args.input, sep = '\t', header = 0, index_col = 0)
    if args.col_sample:
        data_raw = data_raw.T

    # 选取缺失率较低的特征
    drop_feature = (data_raw.isnull().sum(axis=0) / data_raw.shape[0]) > 0.01
    data_drop_feature = data_raw[drop_feature[drop_feature == False].index]
    print(f"    删除缺失率超过0.01的特征 {drop_feature.sum()} 个")
    drop_feature[drop_feature].to_csv(path_or_buf=os.path.join(args.output_dir, 'drop_feature.csv'), sep=',', columns=None, header=True, index=True, index_label=None,mode='w', encoding=None)

    # 去掉有缺失的样本
    drop_sample = data_drop_feature.isnull().sum(axis=1) != 0
    data_clean = data_drop_feature.loc[drop_sample[drop_sample == False].index,]
    print(f"    删除存在缺失的样本 {drop_sample.sum()} 个")
    drop_sample[drop_sample].to_csv(path_or_buf=os.path.join(args.output_dir, 'drop_sample.csv'), sep=',', columns=None, header=True, index=True, index_label=None,mode='w', encoding=None)

    print(f"    剩余矩阵包含 {data_clean.shape[0]} 个样本， {data_clean.shape[1]} 个特征")

    log.info("[process 2] 对每种缺失率情况进行检验")
    result = []
    pool = multiprocessing.Pool(processes=args.thread)
    for missing_rate in missing_rates:
        res = pool.apply_async(fill_and_check_quality, (data_clean, missing_rate, args.output_dir))
        result.append(res)  # 先把结果对象保存
    pool.close()
    pool.join()
    for res in result:
        res.get()

    log.info("[process 3] 缺失率与 MSE 的线性关系")
    mse_dict = {}
    line_names = []
    count = 0
    for missing_rate in missing_rates:
        count += 1
        file_prefix = f"missing_rate_{'%.2f' % missing_rate}"
        with open(os.path.join(args.output_dir, f"{file_prefix}.mse.txt"), 'r') as f:
            for line in f:
                name, mse = line.rstrip().split('\t')
                if not name in mse_dict:
                    mse_dict[name] = []
                mse_dict[name].append(float(mse))
                if count == 1:
                    line_names.append(name)
    fig = plt.figure()
    fig.set_size_inches(8, 8)
    for name in line_names[:2]:
        plt.plot(missing_rates, mse_dict[name], label = name)
    plt.xlabel('missing rate') 
    plt.ylabel('MSE') 
    plt.legend()
    plt.savefig(os.path.join(args.output_dir, f"mse.png"), format = 'png')

