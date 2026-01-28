import numpy as np
from collections import Counter
from pymatgen.core import Structure

def calculate_sconfig(poscar_path, cation_elements=None):
    """
    计算晶体结构的构型熵
    
    Parameters:
    -----------
    poscar_path : str
        POSCAR文件路径
    cation_elements : list, optional
        阳离子元素符号列表，如['Li', 'Ni', 'Co', 'Mn']
        如果为None，则自动识别带正电元素
    
    Returns:
    --------
    float
        构型熵值 (J/mol·K)
    """
    # 1. 读取结构文件
    try:
        struct = Structure.from_file(poscar_path)
    except Exception as e:
        print(f"读取文件失败: {e}")
        return None
    
    # 2. 确定阳离子元素
    if cation_elements is None:
        # 暂时使用常见过渡金属作为示例
        cation_elements = ['Li', 'Na', 'K', 'Mg', 'Ca', 
                          'Ti', 'V', 'Cr', 'Mn', 'Fe', 
                          'Co', 'Ni', 'Cu', 'Zn', 'Pb']
    
    # 3. 收集阳离子元素
    cation_species = []
    for site in struct:
        # 获取位点上的元素
        elements = list(site.species.element_composition.elements)
        if elements and elements[0].symbol in cation_elements:
            cation_species.append(elements[0].symbol)
    
    if not cation_species:
        print("警告：未找到阳离子元素")
        return 0.0
    
    # 4. 统计元素分布
    element_counts = Counter(cation_species)
    total_cations = len(cation_species)
    
    # 5. 计算构型熵
    R = 8.314  # J/mol·K
    entropy = 0.0
    
    for element, count in element_counts.items():
        x_i = count / total_cations
        if x_i > 0:
            entropy -= x_i * np.log(x_i)
    
    sconfig_value = R * entropy
    
    # 6. 输出结果
    print(f"阳离子总数: {total_cations}")
    print("阳离子元素分布:")
    for element, count in element_counts.items():
        fraction = count / total_cations
        print(f"  {element}: {count}个 ({fraction:.3f})")
    print(f"构型熵 Sconfig: {sconfig_value:.4f} J/mol·K")
    print(f"标准化值: {sconfig_value/R:.4f} R")
    
    return sconfig_value

# 使用示例
if __name__ == "__main__":
    # 使用原始字符串处理Windows路径
    file_path = r"D:\github\Disorder_learner\POSCAR"
    result = calculate_sconfig(file_path, cation_elements=['Li', 'Ni', 'Co', 'Mn', 'Pb'])