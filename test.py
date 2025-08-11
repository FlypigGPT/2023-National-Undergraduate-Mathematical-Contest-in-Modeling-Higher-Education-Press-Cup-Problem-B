import numpy as np
import pandas as pd


# 参数设置
beam_angle_deg = 120  # 多波束开角，度
slope_deg = 1.5       # 坡度，度
D0 = 70               # 中心点处海水深度，米
positions = np.array([-800, -600, -400, -200, 0, 200, 400, 600, 800])  # 测线中心点位置


# 角度转弧度
beam_angle = np.radians(beam_angle_deg)
slope = np.radians(slope_deg)
theta = beam_angle / 2  # 最大倾角


# 计算每个位置的水深、覆盖宽度和重叠率
depths = []
coverage_widths = []
overlap_rates = []
spacing = positions[1] - positions[0]  # 200m
for idx, x in enumerate(positions):
    D_x = D0 - x * np.tan(slope)
    W_x = 2 * D_x * np.tan(theta)
    depths.append(D_x)
    coverage_widths.append(W_x)
    if idx == 0:
        overlap_rates.append('—')
    else:
        overlap = (1 - spacing / W_x) * 100
        overlap_rates.append(f"{overlap:.2f}")


# 构建DataFrame
result = pd.DataFrame({
    '测线距中心点处的位置/m': positions,
    '海水深度/m': [f"{d:.2f}" for d in depths],
    '覆盖宽度/m': [f"{w:.2f}" for w in coverage_widths],
    '与前一条测线的重叠率/%': overlap_rates
})



# 打印结果
print(result)
