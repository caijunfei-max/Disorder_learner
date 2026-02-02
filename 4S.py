# Disorder_learner

#导入必要的库
import numpy as np   #numpy数学计算
from pymatgen.core import Structure   #pymatgen处理晶体结构
from collections import Counter  #counter统计元素出现次数

file_path = r"D:\github\Disorder_learner\POSCAR"   #定义POSCAR文件路径，使用原始字符串避免转义问题

def calculate_sconfig(file_path):
    """
    file_path: POSCAR结构文件的路径
    Returns: 该POSCAR结构的构型熵Sconfig
    """
    # 读取POSCAR文件
    struct = Structure.from_file(file_path)
	#使用pymatgen读取POSCAR文件，创建structure对象
    
    # 阳离子位点的邻近阳离子元素分布 
    cation_elements = ['Li', 'Mn']  
    cation_species = []
    
    for site in struct:
        if site.species.element_composition.elements[0].symbol in cation_elements:
            cation_species.append(site.species.element_composition.elements[0].symbol)
    # 检查当前位点的元素是否在阳离子元素列表中
    # element_composition.elements[0].symbol 获取元素的化学符号
	# 如果是阳离子，将其元素符号添加到列表中
			
    element_counts = Counter(cation_species)   # 使用Counter统计每种阳离子元素出现的次数
    total_cations = len(cation_species)   # 计算阳离子总数
    
    # 计算构型熵
    R = 8.314  # 理想气体常数 J/mol·K
    entropy = 0.0   # 初始化理想气体常数和熵值
    
    for element, count in element_counts.items():
        x_i = count / total_cations     # 计算阳离子元素的摩尔分数
        if x_i > 0:
            entropy -= x_i * np.log(x_i)   # 只对正数进行计算，避免log(0)错误
    
<<<<<<< HEAD
    sconfig = R * entropy   # 将标准化熵值乘以理想气体常数得到实际熵值
=======
    sconfig_value = R * entropy   # 将标准化熵值乘以理想气体常数得到实际熵值
>>>>>>> 0114b2556407b6a141464e080a9e11aa98e19f5c
    
    # 输出计算结果
    print(f"阳离子总数: {total_cations}")
    print("阳离子元素分布:")
    for element, count in element_counts.items():
        print(f"  {element}: {count} 个 ({count/total_cations:.3f})")
    print(f"构型熵 Sconfig: {sconfig:.4f} J/mol·K")
    print(f"标准化值: {sconfig/R:.4f} R")
    
    return sconfig

# 调用函数
if __name__ == "__main__":
<<<<<<< HEAD
    calculate_sconfig(file_path)
=======
    calculate_sconfig(file_path)
>>>>>>> 0114b2556407b6a141464e080a9e11aa98e19f5c
