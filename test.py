import pandas as pd
import numpy as np


def get_width(B):
    # 初始化参数
    D_0 = 120  # 海底中心深度（单位：m）
    alpha = 1.5  # 坡度（单位：度）
    theta = 120  # 换能器的开角（单位：度）

    # 计算各距离点的深度
    D = D_0 - distances * np.tan(np.radians(alpha)) * np.cos(np.radians(180 - B))

    # 计算实际方位坡度（投影坡度）
    alpha_proj = np.arctan(abs(np.sin(np.radians(B))) * np.tan(np.radians(alpha))) * 180 / np.pi

    print(f"方向 {B}° 时水深：", D)

    # 覆盖宽度公式（考虑坡度和方向）
    W = D * np.sin(np.radians(theta / 2)) * (
            1 / np.sin(np.radians((180 - theta) / 2 + alpha_proj)) +
            1 / np.sin(np.radians((180 - theta) / 2 - alpha_proj))
    )
    print(f"方向 {B}° 时覆盖宽度：", W)
    return W


# 各测线距离中心点的距离（单位：海里，转米）
distances = np.array([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1])
distances = distances * 1852
print("测线距离（米）：", distances)

# 不同方向
angle = [0, 45, 90, 135, 180, 225, 270, 315]
W = []
for i in angle:
    W.append(get_width(i))