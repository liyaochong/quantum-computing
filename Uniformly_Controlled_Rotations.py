# 本程序用于完成参考文献中的Uniformly Controlled Rotations
# 引用论文：
# Vartiainen J J, Bergholm V, Salomaa M M.
# Transformation of quantum states using uniformly controlled rotations[J].
# Quantum Information & Computation, 2012, 5(6):467-473.
from numpy import *
import numpy as np
import math
from qutip import *

# 定义一个格雷码的类
class GrayCode:
    def __init__(self):
        self.base = ["0", "1"]
    def getNext(self, prelist, z_or_one):
        output = []
        for code in prelist:
            new_code = "%s%s" % (z_or_one,code)
            output.append(new_code)
        if z_or_one == 1:
            output.reverse()
        return output
    def gray(self):
        haf1 = self.getNext(self.base, 0)
        haf2 = self.getNext(self.base, 1)
        ret = haf1 + haf2
        self.base = ret
    def getGray(self, n):
        for i in range(n-1):
            self.gray()
        return self.base

# 定义一个函数：将二进制字符串转为列表类型的向量形式
def str2vec(s):
    s_vec = []
    for i in range(len(s)):
        s_vec.append(int(s[i]))
    return s_vec

# 获得 k 位的二进制编码以及格雷码
def Code_Obtaintion(k):
    N = int(math.pow(2, k))
    B = []
    bin_str = '{:0' + str(k) + 'b}'
    for i in range(0, N):
        k_bin = bin_str.format(i)
        B.append(k_bin)
    # 获得 k 位的格雷码
    a = GrayCode()
    G = a.getGray(k)
    return  G, B

# 定义一个函数，获得变换之后的旋转门对应的参数
def Rotation_Gate_Parameters(k, kernel, G, B):
    N = int(math.pow(2, k))
    M = np.zeros((N, N))
    # 生成论文中的 M 矩阵
    for i in range(N):
        for j in range(N):
            g = G[i]
            b = B[j]
            g_vec = str2vec(g)
            b_vec = str2vec(b)
            t = np.dot(g_vec, b_vec)
            M[i][j] = (2 ** -k) * math.pow(-1, t)
    all_weight = np.tile(kernel, int(N/4))
    all_weight =  all_weight.reshape(N,1)
    theta = np.dot(M, all_weight)
    return theta

# 获得所有受控非门的控制位的location
# 将相邻的格雷码间比特不同的位置取出，存放至location变量之中
# 控制非门所对应的真正索引为 k - location[i]
def CNOT_index(G, k):
    N = int(math.pow(2, k))
    location = []
    for i in range(N - 1):
        str1 = G[i]
        # print(str1)
        str2 = G[i + 1]
        for j in range(len(str1)):
            if str1[j] == str2[j]:
                continue
            else:
                location.append(j)
    # 对于最后的一个码单独地与第一个码相比较，同样的存放至location
    str1 = G[0]
    str2 = G[N - 1]
    for j in range(len(str1)):
        if str1[j] == str2[j]:
            continue
        else:
            location.append(j)
    return location

# 主函数(程序运行开始的地方)
if __name__ == '__main__':
    # 共需要 k+1 个量子比特线路,第一条为目标位，其余的 k 条作为控制位
    # 定义卷积核的内部参数
    kernel = np.array([2 * np.pi / 3, 2 * np.pi / 4, 2 * np.pi / 5, 2 * np.pi / 6])
    # 设定 Uniformly Controlled Rotations 的参数 k
    k = 5
    # 获得 k 位的二进制编码以及格雷码
    G,B = Code_Obtaintion(k)
    # 得到变换之后的旋转门对应的参数
    theta = Rotation_Gate_Parameters(k, kernel, G, B)
    # 获得所有受控非门的控制位的location
    location = CNOT_index(G, k)
    # 创建量子电路
    N = int(math.pow(2, k))
    qc = QubitCircuit(k + 1)
    for i in range (N):
        # 变量angle为旋转门的参数，从转换结果对应位置取出
        angle = theta[i][0]
        # print(angle)
        qc.add_1q_gate("RY", qubits = [0], arg_value = angle)
        # 变量control为受控非门的控制点位置
        control = k - location[i]
        # print(control)
        qc.add_gate("CNOT", targets = 0, controls = control)
        # 生成pdf,查看电路图是否正确(此为测试使用,如在windows下安装的qutip,则会报错)
        # qc.png
    U_list = qc.propagators()
    U = gate_sequence_product(U_list, left_to_right=True)
    # 输出整个电路的酉矩阵
    # print(U)
    # 在电路图中所对应的从下至上是|t> |c3> |c2> |c1>
    # input_state = ket("t c3 c2 c1")
    # 制备输入量子态
    input_state = ket("000000")
    # 将输入量子态张至希尔伯特空间，输出
    print(input_state)
    # 应用电路至输入量子态，获得输出端量子态
    output_state = U * input_state
    # 打印输出量子态
    # 输出量子态由下至上为|t> |c3> |c2> |c1>
    print(output_state)
