# -*- coding: utf-8 -*-
"""
基于坡度、左右非对称条带范围、以及重叠率 η 的多波束测线迭代布设算法
思路：
1. 从西侧边界开始，放置第一条测线，使其左边缘正好在海域西界。
2. 根据上一条测线的右边缘位置 + 设定重叠率，迭代计算下一条测线中心位置。
3. 直到覆盖到东界为止。


"""

import math, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

# ------------------ 参数设置 ------------------
nm_to_m = 1852.0
L_ew = 4.0 * nm_to_m    # 东西向长度（米）
B_ns = 2.0 * nm_to_m    # 南北向长度（米）（此处没用到）
D0 = 110.0              # 海域中心（x=0）水深（米）
alpha_deg = 1.5         # 坡度（度）（东西向）
alpha_rad = math.radians(alpha_deg)
k = math.tan(alpha_rad) # 坡度斜率：水深 = D0 - k * x （x 向东为正，越东越浅）

theta_half_deg = 60.0   # 多波束半开角（度）
theta_half = math.radians(theta_half_deg)

eta = 0.15              # 目标重叠率（比例）

# ------------------ 条带几何函数 ------------------
def left_extent(D):
    """
    输入：中心水深 D（米）
    输出：测线中心到左边缘的水平距离（米）（负值）
    """
    t = math.tan(theta_half)
    return - D * t / (1 + k * t)

def right_extent(D):
    """
    输入：中心水深 D（米）
    输出：测线中心到右边缘的水平距离（米）（正值）
    """
    t = math.tan(theta_half)
    return D * t / (1 - k * t)

def swath_width(D):
    """条带总宽度 = 右边距 - 左边距"""
    return right_extent(D) - left_extent(D)

def depth_at(x):
    """某位置 x（米，中心为0）处的水深"""
    return D0 - k * x

# ------------------ 求第一条测线位置 ------------------
west_bound = -L_ew/2.0   # 西界
east_bound = L_ew/2.0    # 东界

def find_center_left_cover(target_bound):
    """
    找到第一条测线中心位置，使得它的左边缘刚好等于海域西界。
    用二分法求解。
    """
    a = west_bound
    b = east_bound
    for _ in range(60):
        m = 0.5*(a+b)
        val = m + left_extent(depth_at(m)) - target_bound
        if abs(val) < 1e-12:
            return m
        if (a + left_extent(depth_at(a)) - target_bound) * val <= 0:
            b = m
        else:
            a = m
    return 0.5*(a+b)

x0 = find_center_left_cover(west_bound)

# ------------------ 迭代布设测线 ------------------
centers = [x0]
max_iters = 500

for it in range(max_iters):
    xi = centers[-1]                    # 当前测线中心
    Di = depth_at(xi)                    # 当前测线中心水深
    prev_right = xi + right_extent(Di)   # 当前测线右边界位置

    # 如果已经覆盖到东界，停止
    if prev_right >= east_bound - 1e-6:
        break

    # 定义 f(x)：上一条测线右边缘 到 下一条测线左边缘 的距离 - 重叠控制量
    def f(x):
        D = depth_at(x)
        return prev_right - (x + left_extent(D)) - eta * swath_width(D)

    # 在区间 (xi, east_bound) 内找 f(x) = 0 的点
    a = xi + 1e-6
    b = east_bound
    fa = f(a)
    fb = f(b)

    if fa < 0:
        # 如果一开始就负，说明无法再放下一条
        break

    # 检查是否有符号变化，没有的话就用近似公式
    if fa * fb > 0:
        approx = xi + (1-eta) * swath_width(depth_at(xi))
        centers.append(approx)
        continue

    # 二分法求解
    for _ in range(60):
        m = 0.5*(a+b)
        fm = f(m)
        if fm == 0 or (b-a) < 1e-6:
            a = b = m
            break
        if fa * fm <= 0:
            b = m
            fb = fm
        else:
            a = m
            fa = fm

    x_next = 0.5*(a+b)
    centers.append(x_next)

# ------------------ 结果计算 ------------------
centers = np.array(centers)
Ns = len(centers)

lefts = centers + np.array([left_extent(depth_at(x)) for x in centers])
rights = centers + np.array([right_extent(depth_at(x)) for x in centers])
W_vals = np.array([swath_width(depth_at(x)) for x in centers])

covered_left = lefts.min()
covered_right = rights.max()

# 生成结果表
df = pd.DataFrame({
    '测线中心x坐标_米': centers,
    '中心水深_米': [round(depth_at(x),3) for x in centers],
    '左边界x坐标_米': lefts,
    '右边界x坐标_米': rights,
    '条带宽度_米': W_vals
})

# 输出结果
print("===== 测线布设结果 =====")
print(df.to_string(index=False))
print("\n总测线条数:", Ns)
print(f"覆盖范围: 从 {covered_left:.2f} m 到 {covered_right:.2f} m")
print(f"总覆盖长度: {(covered_right-covered_left)/1000:.3f} km")

# ------------------ 可视化图表 ------------------

# 1. 测线布设示意图
print("生成图表 1/3: 测线布设示意图")
fig1 = plt.figure(figsize=(12, 8))
colors = plt.cm.viridis(np.linspace(0, 1, Ns))

for i, (center, left, right, width) in enumerate(zip(centers, lefts, rights, W_vals)):
    # 绘制测线覆盖区域
    rect = Rectangle((left, i-0.3), width, 0.6, 
                     facecolor=colors[i], alpha=0.7, edgecolor='black', linewidth=1)
    plt.gca().add_patch(rect)
    # 绘制测线中心
    plt.plot(center, i, 'ro', markersize=8, markeredgecolor='black')

# 绘制海域边界
plt.axvline(west_bound, color='red', linestyle='--', linewidth=2, label='西界')
plt.axvline(east_bound, color='red', linestyle='--', linewidth=2, label='东界')

plt.xlabel('东西向距离 (m)')
plt.ylabel('测线编号')
plt.title('测线布设示意图', fontsize=14, fontweight='bold')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# 2. 水深剖面图
print("生成图表 2/3: 水深剖面图")
fig2 = plt.figure(figsize=(12, 8))
x_range = np.linspace(west_bound-1000, east_bound+1000, 1000)
depths = [depth_at(x) for x in x_range]

plt.plot(x_range, depths, 'b-', linewidth=3, label='海底地形')
plt.scatter(centers, [depth_at(x) for x in centers], 
           c=range(Ns), cmap='viridis', s=100, 
           edgecolors='black', linewidth=1, label='测线中心')

plt.xlabel('东西向距离 (m)')
plt.ylabel('水深 (m)')
plt.title('水深剖面与测线中心位置', fontsize=14, fontweight='bold')
plt.legend()
plt.grid(True, alpha=0.3)
plt.gca().invert_yaxis()  # 水深向下为正
plt.tight_layout()
plt.show()

# 3. 条带宽度变化图
print("生成图表 3/3: 条带宽度变化图")
fig3 = plt.figure(figsize=(12, 8))
plt.plot(centers, W_vals, 'o-', linewidth=2, markersize=8, 
         color='green', markeredgecolor='black')
plt.fill_between(centers, W_vals, alpha=0.3, color='green')

plt.xlabel('测线中心位置 (m)')
plt.ylabel('条带宽度 (m)')
plt.title('条带宽度随位置变化', fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3)

# 添加趋势线
z = np.polyfit(centers, W_vals, 1)
p = np.poly1d(z)
plt.plot(centers, p(centers), "r--", alpha=0.8, label=f'趋势线: y={z[0]:.3f}x+{z[1]:.1f}')
plt.legend()
plt.tight_layout()
plt.show()


print("所有图表生成完成！")
