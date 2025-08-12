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
# 水深计算（x方向分量）：
#    D(x) = D0 - x * 1852 * tan(α) * cos(β)



import numpy as np
import pandas as pd

# 参数
beam_angle_deg = 120
slope_deg = 1.5
D0 = 120
theta_deg = beam_angle_deg
alpha_rad = np.radians(slope_deg)

# x为海里，转为米
x_nm = np.array([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1])
x_m = x_nm * 1852

# β为度，转为弧度
beta_deg = np.array([0, 45, 90, 135, 180, 225, 270, 315])
beta_rad = np.radians(beta_deg)

# 结果矩阵
cover_width = np.zeros((len(beta_rad), len(x_m)))

# 计算
for i, b in enumerate(beta_rad):
    for j, x in enumerate(x_m):
        # 水深公式
        D_x = D0 - x * np.tan(alpha_rad) * np.cos(b)
        # delta计算
        delta = np.arcsin(
            np.abs(np.sin(b) * np.sin(alpha_rad)) /
            np.sqrt(np.cos(b)**2 * np.cos(alpha_rad)**2 + np.sin(b)**2)
        )
        delta_deg = np.degrees(delta)
        # 覆盖宽度公式
        sin_theta2 = np.sin(np.radians(theta_deg / 2))
        sin1 = np.sin(np.radians((180 - theta_deg) / 2 + delta_deg))
        sin2 = np.sin(np.radians((180 - theta_deg) / 2 - delta_deg))
        W_x = D_x * sin_theta2 * (1 / sin1 + 1 / sin2)
        cover_width[i, j] = W_x

# 构建DataFrame
df = pd.DataFrame(cover_width, index=beta_deg, columns=x_nm)
df.index.name = '测线方向夹角/°'
df.columns.name = '距离/海里'

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
print('覆盖宽度/m:')
print(df.round(3))
