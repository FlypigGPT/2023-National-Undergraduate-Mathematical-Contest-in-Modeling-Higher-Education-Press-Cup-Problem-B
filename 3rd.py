import numpy as np
import math

# 问题3：最短测量长度的测线设计，满足重叠率10%~20%
# --------------------------------------------------
# 海域：南北长2海里，东西宽4海里，中心水深110m，西深东浅，坡度1.5°，开角120°
# 目标：测线数量最少，总长度最短，重叠率10%~20%，完全覆盖

# 参数
beam_angle_deg = 120
slope_deg = 1.5
D0 = 110
theta = np.radians(beam_angle_deg / 2)
alpha = np.radians(slope_deg)

width_nm = 4  # 东西宽4海里
length_nm = 2 # 南北长2海里
width_m = width_nm * 1852
length_m = length_nm * 1852

# 东西向坐标，x=0为中心，x∈[-2, 2]海里
x_edges_nm = np.array([-2, 2])
x_edges_m = x_edges_nm * 1852

# 计算两端水深（西深东浅，符号为减）
D_west = D0 - x_edges_m[0] * np.tan(alpha)
D_east = D0 - x_edges_m[1] * np.tan(alpha)

# 计算两端覆盖宽度
W_west = 2 * D_west * np.tan(theta)
W_east = 2 * D_east * np.tan(theta)
W_min = min(W_west, W_east)

# 满足重叠率10%~20%的测线间距范围
S_max = 0.90 * W_min  # 重叠率10%
S_min = 0.80 * W_min  # 重叠率20%

# 取S_max保证测线数量最少且全覆盖
S = S_max
num_lines = math.ceil(width_m / S) + 1
actual_spacing = width_m / (num_lines - 1)

# 实际重叠率（用最小覆盖宽度计算，保证最差情况）
actual_overlap = (1 - actual_spacing / W_min) * 100

# 总测量长度
total_length = num_lines * length_m

print(f"东西两端水深: 西 {D_west:.2f} m, 东 {D_east:.2f} m")
print(f"东西两端覆盖宽度: 西 {W_west:.2f} m, 东 {W_east:.2f} m")
print(f"最小覆盖宽度: {W_min:.2f} m")
print(f"测线间距范围: {S_min:.2f} m ~ {S_max:.2f} m")
print(f"最优测线间距: {actual_spacing:.2f} m")
print(f"测线数量: {num_lines}")
print(f"实际重叠率: {actual_overlap:.2f}%")
print(f"总测量长度: {total_length/1000:.2f} km")
D_west = D0 - x_edges_m[0] * np.tan(alpha)
D_east = D0 - x_edges_m[1] * np.tan(alpha)
D_west = D0 - x_edges_m[0] * np.tan(alpha)
D_east = D0 - x_edges_m[1] * np.tan(alpha)