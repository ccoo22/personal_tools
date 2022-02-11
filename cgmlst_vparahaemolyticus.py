#!/home/genesky/software/python/3.8.12/bin/python3

import os
import glob
import argparse
from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.common.by import By


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description="做http://pubmlst.org/vparahaemolyticus物种的cgmlst分型", formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('--input', '-i', type = str, required=True, help = "物种的组装fasta序列")
    parser.add_argument('--output_dir', '-o', type = str, required=True, help = "结果输出目录, 例如 ./test")
    args = parser.parse_args()

    return args

args = set_and_parse_args()
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# 转绝对路径
args.input = os.path.abspath(args.input)
args.output_dir = os.path.abspath(args.output_dir)

print("[process] 初始化浏览器")
# 配置chrome浏览器的参数
options = webdriver.ChromeOptions()
options.add_argument('--headless')  # 确保无头
options.add_argument('--disable-gpu')  # 无需要gpu加速
options.add_argument('--no-sandbox')  # 无沙箱#
#禁止图片/css/javascript加载，加快html导入速度


# 初始化浏览器对象，依赖chromedriver软件，请提前安装
# 打开浏览器操作，速度比较慢
driver  = webdriver.Chrome(executable_path="/usr/bin/chromedriver", options=options)
 
# 设置浏览器窗口大小
driver.set_window_size(1920,1080)


# 在当前窗口打开一个网页
print("[process] 打开目标页面 https://pubmlst.org/bigsdb?db=pubmlst_vparahaemolyticus_seqdef&page=sequenceQuery")
driver.get('https://pubmlst.org/bigsdb?db=pubmlst_vparahaemolyticus_seqdef&page=sequenceQuery')
driver.save_screenshot(os.path.join(args.output_dir, '1.homepage.png'))


# 选择cgmlst分析
print("[process] 选择cgmlst分析，并上传fasta文件")
select_object = Select(driver.find_element(by='name', value='locus'))
select_object.select_by_value('SCHEME_3')
driver.save_screenshot(os.path.join(args.output_dir, '2.select_cgmlst.png'))

# 上传fasta文件
driver.find_element(by='id', value='fasta_upload').send_keys(args.input)  # 一定要给绝对路径
driver.save_screenshot(os.path.join(args.output_dir, '3.fasta_upload.png'))

# 单击，开始运行, 可能需要一点时间
print("[process] 开始运行，需要1-2分钟左右")
driver.find_element(by='name', value='submit').click()
# 每5秒检测一次 resultstable 元素是否出现，最多等待300秒，超时后中断运行
WebDriverWait(driver, 300, 5).until(EC.presence_of_element_located((By.ID, "resultstable")))  # 检查resultstable是否出现，如果出现了，则分析完成
driver.save_screenshot(os.path.join(args.output_dir, '4.finished.png'))

# 分析完成，开始下载
print("[process] 分析完成，开始下载")
add_all = driver.find_element(by='id', value='resultstable').find_elements_by_tag_name('a')
add_txt = add_all[-2].get_attribute('href')  # txt 文件下载地址
print(f"          下载地址 {add_txt}")
output_file = os.path.join(args.output_dir, 'cgmlst.txt')
os.system(f"wget {add_txt} -O {output_file} --no-check-certificate")

# 关闭浏览器, driver.close()  关闭当前页面，如果只有一个页面，则等同于quit()
driver.quit()

if os.path.exists(output_file) and os.path.getsize(output_file) > 1000:
    print("[OK] 分析完成")
else:
    print("[Error] 结果文件不存在或太小，可能分析失败")

 