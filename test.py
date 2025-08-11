import numpy as np
import pandas as pd

# 问题2 多波束测深覆盖宽度建模思路
# --------------------------------------------------
# 已知：
#   多波束开角120°，坡度1.5°，中心点水深120m
#   测线方向与坡面法向夹角β
#   船距中心点距离x（单位：海里，1海里=1852米）
# 求：不同β、不同x下的覆盖宽度
#
# 主要变量：
#   α = 坡度（1.5°），θ = 60°，D0 = 120m
#   β = 测线与坡面法向夹角
#   x = 距中心点距离
#
# 1. 水深计算（x方向分量）：
#    D(x) = D0 + x * 1852 * tan(α) * cos(β)
# 2. 覆盖宽度：
#    W(x, β) = 2 * D(x) * tan(θ)

# 参数
beam_angle_deg = 120
slope_deg = 1.5
D0 = 120
theta = np.radians(beam_angle_deg / 2)
alpha = np.radians(slope_deg)

# x为海里，转为米
x_nm = np.array([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1])
x_m = x_nm * 1852


# β为度，转为弧度（0到360，步长45°）
# β为度，转为弧度
beta_deg = np.array([0, 45, 90, 135, 180, 225, 270, 315])

beta = np.radians(beta_deg)

# 结果矩阵
cover_width = np.zeros((len(beta), len(x_m)))

# 计算
for i, b in enumerate(beta):
    for j, x in enumerate(x_m):
        D_x = D0 + x * np.tan(alpha) * np.cos(b)
        W_x = 2 * D_x * np.tan(theta)
        cover_width[i, j] = W_x

# 构建DataFrame
df = pd.DataFrame(cover_width, index=beta_deg, columns=x_nm)
df.index.name = '测线方向夹角/°'
df.columns.name = '距离/海里'

# 打印结果
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
print('覆盖宽度/m:')
print(df.round(3))
