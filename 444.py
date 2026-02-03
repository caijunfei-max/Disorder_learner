# Disorder_learner

# 导入必要的库
import os                                                                #os处理文件路径
import sys                                                               #环境交互
import glob                                                              #glob匹配目录下以POSCAR开头的文件
import numpy as np                                                       #numpy数学计算
from pymatgen.core import Structure                                      #pymatgen处理晶体结构
from collections import Counter, defaultdict                             #counter, defaultdict统计元素出现次数
from scipy.spatial import cKDTree                                        #cKDTree查找近邻原子
from scipy.stats import entropy                                          #熵计算

# 使用os.path.join进行灵活的路径管理，兼容交互式环境和脚本运行
try:
    # 获取当前脚本的目录（在脚本文件中运行）
    work_path = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # 在Jupyter Notebook交互式环境中，使用当前工作目录
    work_path = os.getcwd()
    print(f"注意：在交互式环境中运行，使用当前工作目录: {work_path}")

structure_dir = os.path.join(work_path, "structure_file")                 #结构文件目录

def find_poscar_files(structure_dir):
    """查找structure_file目录下所有以POSCAR开头的文件"""
    poscar_patterns = [
        os.path.join(structure_dir, "POSCAR*"),                           # POSCAR开头的文件作为输入结构
    ]
    
    poscar_files = []
    for pattern in poscar_patterns:
        poscar_files.extend(glob.glob(pattern))
    poscar_files = sorted(list(set(poscar_files)))
    return poscar_files

def calculate_sconfig(file_path, cutoff_radius=3.8, z_tolerance=0.2, target_coordination=6):
    """
    file_path: POSCAR结构文件的路径
	cutoff_radius: 截断半径设为3.8Å，用来区分近邻位点
	z_tolerance：c轴方向容差设为0.2，仅识别同一层ab平面内的近邻位点
	target_coordination：阳离子配位数设为6，取6个近邻阳离子
	Returns: 该POSCA结构的构型熵Sconfig
    """
    # 从文件路径提取结构名称
    file_name = os.path.basename(file_path)
    struct_name = file_name.replace("POSCAR", "").replace("_", " ").strip()
    if not struct_name:
        struct_name = "默认结构"
    
    # 检查文件是否存在
    if not os.path.exists(file_path):
        print(f"错误: 文件 {file_path} 不存在")
        return None
    
    # 使用pymatgen读取POSCAR文件，创建structure对象
    struct = Structure.from_file(file_path)
    
    # 定义阳离子元素
    all_cation_elements = ['Li', 'Na', 'K', 'Rb', 'Mg', 'Ca', 'Sr', 'Sc', 
                          'Y', 'Ti', 'Zr', 'V', 'Nb', 'Cr', 'Mo', 'Mn', 'Fe', 'Ru', 
                          'Co', 'Rh', 'Ni', 'Rd', 'Cu', 'Ag','Zn', 'Cd', 'Al', 'Ga', 
                          'In', 'Ge', 'Sn', 'Sb', 'Se', 'Te', 'Lu', 'Hf', 'Ta', 'W', 'Ir', 'Pt', 'Pb', 'Bi']
    
    
    # 1. 按照整个结构中各阳离子的数量计算——整体组成熵
    print("\n1. 组成构型熵")
    
    all_cation_species = []
    for site in struct:
        element = site.species.element_composition.elements[0].symbol
        if element in all_cation_elements:
            all_cation_species.append(element)
    # 检查当前位点的元素是否在阳离子元素列表中
    # element_composition.elements[0].symbol 获取元素的化学符号
	
    element_counts = Counter(all_cation_species)
    total_cations = len(all_cation_species)                                  # 计算阳离子总数
    
    print(f"阳离子总数: {total_cations}")
    for element, count in sorted(element_counts.items()):
        fraction = count / total_cations
        print(f"  {element}: {count}个 ({fraction:.3f})")
    
    R = 8.314                                                                # 理想气体常数 J/mol·K
    global_entropy = 0.0                                                     # 初始化理想气体常数和熵值
    for element, count in element_counts.items():
        x_i = count / total_cations                                          # 计算阳离子元素的摩尔分数
        if x_i > 0:
            global_entropy -= x_i * np.log(x_i)                              # 只对正数进行计算，避免log(0)错误
    sconfig_global = R * global_entropy                                      # 将标准化熵值乘以理想气体常数得到实际熵值
    
    print(f"  组成构型熵 (S_config): {sconfig_global:.4f} J/(mol·K)")
    
    # 2. 各阳离子位点的局域环境分析
    print("\n2. 阳离子位点的局域环境分析")
    print(f"截断半径: {cutoff_radius}Å, c方向容差: {z_tolerance}Å, 配位数: {target_coordination}")
    
    lattice = struct.lattice                                                 # 输出晶格参数
    a, b, c = lattice.a, lattice.b, lattice.c
    print(f"晶格参数: a={a:.3f}Å, b={b:.3f}Å, c={c:.3f}Å")
    
    layers_x = int(np.ceil(cutoff_radius / a)) + 1                           # 计算结构层数
    layers_y = int(np.ceil(cutoff_radius / b)) + 1
    
    cation_sites = []                                                        # 收集所有阳离子位点
    for i, site in enumerate(struct):
        element = site.species.element_composition.elements[0].symbol
        if element in all_cation_elements:
            cation_sites.append({
                'index': i,                                                  # 位点在结构中的索引
                'element': element,
                'frac_coords': site.frac_coords,                             # 分数坐标
                'coords': site.coords,                                       # 笛卡尔坐标,
            })
    
    print(f" {len(cation_sites)} 个阳离子位点")
    
    # 通过创建镜像位点以模拟晶体结构的周期性
    all_cation_coords = []
    all_cation_indices = []
    all_cation_elements_list = []
    all_cation_frac_z = []
    all_cation_original_coords = []
    
    for site in cation_sites:
        base_coords = site['coords']
        base_index = site['index']
        base_element = site['element']
        base_frac_z = site['frac_coords'][2]
        
        # 添加中心单元位点
        all_cation_coords.append(base_coords)
        all_cation_indices.append(base_index)
        all_cation_elements_list.append(base_element)
        all_cation_frac_z.append(base_frac_z)
        all_cation_original_coords.append(base_coords)
        
        # 在ab平面内创建周期镜像，跳过该中心位点
        for dx in range(-layers_x, layers_x + 1):                             # +1确保能覆盖截断半径内的所有可能近邻
            for dy in range(-layers_y, layers_y + 1):
                if dx == 0 and dy == 0:
                    continue
                
                mirror_coords = base_coords.copy()                            # 在ab平面内沿晶格向量方向平移
                mirror_coords[0] += dx * a
                mirror_coords[1] += dy * b
                
                all_cation_coords.append(mirror_coords)
                all_cation_indices.append(base_index)
                all_cation_elements_list.append(base_element)
                all_cation_frac_z.append(base_frac_z)
                all_cation_original_coords.append(base_coords)
    
    # 创建KDTree
    all_cation_coords = np.array(all_cation_coords)
    kdtree = cKDTree(all_cation_coords)
    
    all_environments = []
    site_details = []
    valid_sites = 0
    
    print(f"\n分析每个阳离子位点的局域环境...")
    
    for i, site in enumerate(cation_sites):
        site_index = site['index']
        site_element = site['element']
        site_frac_z = site['frac_coords'][2]
        site_coords = site['coords']
        
        try:
            # 收集近邻阳离子位点
            k = min(len(all_cation_coords), target_coordination * 3)
            distances, indices = kdtree.query(site_coords, k=k)
            
            # 筛选ab平面内的近邻位点，跳过该中心位点
            ab_plane_neighbors = []
            for j, idx in enumerate(indices):
                if j == 0:  
                    continue
                
                neighbor_coords = all_cation_coords[idx]
                neighbor_global_idx = all_cation_indices[idx]
                neighbor_element = all_cation_elements_list[idx]
                neighbor_frac_z = all_cation_frac_z[idx]
                distance = distances[j]
                
                # 根据截断半径进行距离筛选
                if distance > cutoff_radius:
                    continue
                
                # 根据c方向容差进行ab平面内筛选
                z_diff = abs(neighbor_frac_z - site_frac_z)
                z_diff_periodic = min(z_diff, 1 - z_diff) * c
                if z_diff_periodic > z_tolerance:
                    continue
                
                # 检查是否为阳离子
                if neighbor_element in all_cation_elements:
                    original_neighbor_coords = all_cation_original_coords[idx]
                    is_periodic_image = not np.allclose(neighbor_coords, original_neighbor_coords)
                    
                    ab_plane_neighbors.append({
                        'element': neighbor_element,
                        'distance': distance,
                        'site_index': neighbor_global_idx,
                        'z_diff': z_diff_periodic,
                        'is_periodic_image': is_periodic_image,
                    })
            
            # 按距离排序
            ab_plane_neighbors.sort(key=lambda x: x['distance'])
            physical_site_neighbors = defaultdict(list)
            for neighbor in ab_plane_neighbors:
                physical_site_neighbors[neighbor['site_index']].append(neighbor)
            
            # 为每个位点选择最近的镜像位点
            unique_neighbors = []
            for site_idx, mirrors in physical_site_neighbors.items():
                closest_mirror = min(mirrors, key=lambda x: x['distance'])
                unique_neighbors.append(closest_mirror)
            
            # 每个位点取6个近邻阳离子
            unique_neighbors.sort(key=lambda x: x['distance'])
            selected_neighbors = unique_neighbors[:target_coordination]
            if len(selected_neighbors) >= target_coordination:
                valid_sites += 1
                
                # 近邻阳离子位点分析
                neighbor_elements = [n['element'] for n in selected_neighbors]
                distances = [n['distance'] for n in selected_neighbors]
                z_diffs = [n['z_diff'] for n in selected_neighbors]
                periodic_count = sum(1 for n in selected_neighbors if n['is_periodic_image'])
                
                # 创建环境指纹
                element_counter = Counter(neighbor_elements)
                fingerprint = tuple(sorted(neighbor_elements))
                
                all_environments.append(fingerprint)
                site_details.append({
                    'site_index': site_index,
                    'element': site_element,
                    'fingerprint': fingerprint,
                    'neighbor_elements': neighbor_elements,
                    'neighbor_distances': distances,
                    'num_neighbors': len(selected_neighbors),
                    'periodic_neighbors': periodic_count,
                })
        
        except Exception as e:
            print(f"阳离子位点 {site_index} 分析失败: {e}")
            continue
    
    print(f"位点数 (有{target_coordination}个ab平面内阳离子近邻): {valid_sites}/{len(cation_sites)}")
    
    # 3. 局域无序构型熵计算
    print("\n3. 局域环境无序构型计算")
    
    if valid_sites == 0:
        print("警告: 没有找到足够的有效位点进行局域环境分析")
        return {
            'struct_name': struct_name,
            'file_name': file_name,
            'sconfig_global': sconfig_global,
            'local_disorder': 0.0,
            'valid_sites': 0,
            'total_sites': len(cation_sites)
        }
    
    environment_counts = Counter(all_environments)                                         # 统计不同局域环境类型的出现频率
    total_environments = len(all_environments)
    
    print(f"共 {len(environment_counts)} 种不同的局域环境")
    for i, (env, count) in enumerate(environment_counts.most_common(10)):                  # 显示前10种局域环境类型
        freq = count / total_environments
        env_str = ", ".join([f"{e}{env.count(e)}" for e in sorted(set(env))])
        print(f"  环境{i+1}: [{env_str}] - {count}个位点 ({freq:.3f})")
    if len(environment_counts) > 10:
        print(f"  ... 还有 {len(environment_counts) - 10} 种其他环境")
    
    # 计算每个位点的局域环境熵
    local_shannon_entropy = 0.0
    for env, count in environment_counts.items():
        p_env = count / total_environments
        if p_env > 0:
            local_shannon_entropy -= p_env * np.log(p_env)
    
    local_disorder = R * local_shannon_entropy
    
    print(f"  局域环境无序构型熵: {local_disorder:.4f} J/(mol·K)")
    
    # 4. 输出计算结果
    print(f"结构: {struct_name}")
    print(f"组成构型熵 (S_config_global): {sconfig_global:.4f} J/(mol·K)")
    print(f"局域环境无序构型熵 (S_config_local): {local_disorder:.4f} J/(mol·K)")
    
    results = {
        'struct_name': struct_name,
        'file_name': file_name,
        'sconfig_global': sconfig_global,
        'local_disorder': local_disorder,
        'valid_sites': valid_sites,
        'total_sites': len(cation_sites)
    }
    
    return results

def batch_process_structures():
    """批量处理所有POSCAR文件"""
    # 确保结构目录存在
    if not os.path.exists(structure_dir):
        print(f"创建结构文件目录: {structure_dir}")
        os.makedirs(structure_dir, exist_ok=True)
        print(f"请将POSCAR文件放置在以下目录: {structure_dir}")
        print(f"支持的命名格式:")
        print(f"  - POSCAR_结构名 (例如: POSCAR_perovskite)")
        print(f"  - 结构名_POSCAR (例如: perovskite_POSCAR)")
        print(f"  - POSCAR (单个文件)")
        return {}
    
    # 查找所有POSCAR文件
    poscar_files = find_poscar_files(structure_dir)
    
    if not poscar_files:
        print(f"在目录 {structure_dir} 中未找到POSCAR文件")
        print(f"当前工作目录: {work_path}")
        print(f"查找的目录: {structure_dir}")
        return {}
    
    print(f"找到 {len(poscar_files)} 个POSCAR文件:")
    for i, file_path in enumerate(poscar_files, 1):
        print(f"  {i}. {os.path.basename(file_path)}")
    
    all_results = {}
    for file_path in poscar_files:
        try:
            print(f"\n{'-'*60}")
            print(f" {os.path.basename(file_path)}")
            result = calculate_sconfig(file_path, cutoff_radius=3.8, z_tolerance=0.2, target_coordination=6)
            if result:
                all_results[result['struct_name']] = result
        except Exception as e:
            file_name = os.path.basename(file_path)
            print(f"处理文件 {file_name} 时出错: {e}")
            continue
    
    return all_results

def print_summary(all_results):
    """处理结果汇总"""
    if not all_results:
        print("\n没有成功处理任何结构文件")
        return
    
    print("="*100)
    print(f"{'结构名称':<20} {'文件名':<25} {'组成构型熵':<15} {'局域环境熵':<15} {'有效位点':<10}")
    print("-"*100)
    
    for struct_name, result in all_results.items():
        print(f"{struct_name:<20} {result['file_name']:<25} "
              f"{result['sconfig_global']:.4f} J/(mol·K)  "
              f"{result['local_disorder']:.4f} J/(mol·K)  "
              f"{result['valid_sites']}/{result['total_sites']}")
    
    print("-"*100)

# 调用函数
if __name__ == "__main__":
    print(f"工作目录: {work_path}")
    print(f"结构文件目录: {structure_dir}")
    
    all_results = batch_process_structures()
    print_summary(all_results)
    
    # 如果没有任何结果，提供更多帮助信息
    if not all_results:
        print(f"\n请确保:")
        print(f"1. 在 {structure_dir} 目录中放置了POSCAR文件")
        print(f"2. 文件命名符合支持的格式")
        print(f"3. 如果使用Jupyter Notebook，可能需要重启内核后运行")