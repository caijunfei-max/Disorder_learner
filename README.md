# Disorder\_learner

For coding of the classification of disorder in cathode intralyer metals.
该python代码用于计算输入POSCAR晶体结构构型熵（Sconfig）：

1.读取VASP格式的POSCAR结构文件
2.识别结构中的阳离子位点（仅能识别给出的元素种类）
3.统计阳离子的元素种类和数量
4.计算阳离子的构型熵
5.输出计算结果和统计信息

ZHY test

存疑：
1.当前代码仅可读入给定文件路径的POSCAR，是否不够灵活?
2.阳离子位点识别与配位构型的优化；
