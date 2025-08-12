import pandas as pd
import numpy as np

# 初始化参数
D_0 = 70                         # 中心点水深（米）
theta = 120                      # 多波束开角（度）
alpha = 1.5                      # 坡度（度）
d = 200                          # 测线间距（米）

# 修正条带间距（考虑坡度与投影角度）
d_proj = d * np.sin(np.radians(90 + theta / 2)) / np.sin(np.radians(90 + alpha - theta / 2))

# 各测线距中心点的距离
distances = np.array([-800, -600, -400, -200, 0, 200, 400, 600, 800])

# 每个点的水深
D = D_0 - distances * np.tan(np.radians(alpha))

# 覆盖宽度（更精确公式，考虑坡度和条带角度）
W = D * np.sin(np.radians(theta / 2)) * (
    1 / np.sin(np.radians((180 - theta) / 2 + alpha)) +
    1 / np.sin(np.radians((180 - theta) / 2 - alpha))
)

# 重叠率
n = (1 - d_proj / W) * 100

# 创建 DataFrame 用于保存结果
df = pd.DataFrame({'测线距中心点处的距离/m': distances})
df['海水深度/m'] = D
df['覆盖宽度/m'] = W
df['与前一条测线的重叠率/%'] = n

print(df)