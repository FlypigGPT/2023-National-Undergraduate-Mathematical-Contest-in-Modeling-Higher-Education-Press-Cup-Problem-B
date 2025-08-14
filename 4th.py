import math, numpy as np, pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

# 设置图表样式
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

filename = "附件.xlsx"
L_ew_nm = 4.0
B_ns_nm = 5.0
beam_opening_deg = 120.0
eta_target = 0.20
overlap_threshold = 0.20  # >20%
# ==================================

nm_to_m = 1852.0
theta_half_deg = beam_opening_deg / 2.0

# ---------- 读取 Excel ----------
xls = pd.read_excel(filename, header=None)
xls = xls.dropna(axis=0, how='all').dropna(axis=1, how='all')
arr = xls.values.astype(object)

is_num = np.vectorize(lambda x: isinstance(x, (int, float, np.integer, np.floating)))(arr)
rows = np.where(is_num.any(axis=1))[0]
cols = np.where(is_num.any(axis=0))[0]
r0, r1 = rows[0], rows[-1]
c0, c1 = cols[0], cols[-1]
sub = arr[r0:r1 + 1, c0:c1 + 1]

subf = np.full(sub.shape, np.nan, dtype=float)
for i in range(sub.shape[0]):
    for j in range(sub.shape[1]):
        try:
            subf[i, j] = float(sub[i, j])
        except:
            subf[i, j] = np.nan

dfsub = pd.DataFrame(subf)
if dfsub.isna().values.any():
    dfsub = dfsub.interpolate(axis=0, limit_direction='both').interpolate(axis=1, limit_direction='both')
subf = dfsub.values

depth = subf.copy()
nrows, ncols = depth.shape
print(f"深度矩阵：{nrows} 行 x {ncols} 列 (单位 m)")

dx = (L_ew_nm * nm_to_m) / max(1, (ncols - 1))
dy = (B_ns_nm * nm_to_m) / max(1, (nrows - 1))
print(f"网格间距：dx={dx:.2f} m, dy={dy:.2f} m")

x_coords = np.linspace(-L_ew_nm * nm_to_m / 2.0, L_ew_nm * nm_to_m / 2.0, ncols)
y_coords = np.linspace(-B_ns_nm * nm_to_m / 2.0, B_ns_nm * nm_to_m / 2.0, nrows)


# ---------- 基本几何函数 ----------
def left_extent_geom(D, alpha_deg):
    denom = math.sin(math.radians(90.0 - theta_half_deg + alpha_deg))
    if abs(denom) < 1e-9:
        denom = 1e-9
    return - D * math.sin(math.radians(theta_half_deg)) / denom


def right_extent_geom(D, alpha_deg):
    denom = math.sin(math.radians(90.0 - theta_half_deg - alpha_deg))
    if abs(denom) < 1e-9:
        denom = 1e-9
    return D * math.sin(math.radians(theta_half_deg)) / denom


# ---------- 插值函数 ----------
def build_interp_funcs_for_axis(depth_mat, axis):
    funcs = []
    if axis == 'x':
        for i in range(nrows):
            funcs.append(interp1d(x_coords, depth_mat[i, :], bounds_error=False, fill_value="extrapolate"))
    else:
        for j in range(ncols):
            funcs.append(interp1d(y_coords, depth_mat[:, j], bounds_error=False, fill_value="extrapolate"))
    return funcs


# ---------- 计算每条带在沿线上每个样点的 left/right（返回二维数组） ----------
def per_center_left_right_arrays(center_xc, interp_funcs, axis):
    """
    返回两个数组 Ls（长度 = 沿线样点数），Rs，以及 Ds（深度），Alphas（坡度角度）
    每个沿样点上基于局部深度与差分估算 alpha 计算 left/right（全局坐标）
    """
    m = len(interp_funcs)
    Ls = np.zeros(m)
    Rs = np.zeros(m)
    Ds = np.zeros(m)
    Alphas = np.zeros(m)
    for ia in range(m):
        D = float(interp_funcs[ia](center_xc))
        eps = max(dx * 0.5, 1.0) if axis == 'x' else max(dy * 0.5, 1.0)
        Dp = float(interp_funcs[ia](center_xc + eps))
        Dm = float(interp_funcs[ia](center_xc - eps))
        k_local = (Dp - Dm) / (2.0 * eps)
        alpha_deg = math.degrees(math.atan(k_local))
        l = left_extent_geom(D, alpha_deg)
        r = right_extent_geom(D, alpha_deg)
        Ls[ia] = center_xc + l
        Rs[ia] = center_xc + r
        Ds[ia] = D
        Alphas[ia] = alpha_deg
    return Ls, Rs, Ds, Alphas


# ---------- 逐条迭代放置中心（与之前相同，但我们会存储每条带的 per-sample L/R 数组） ----------
def iterative_centers_with_per_sample(axis, eta=eta_target):
    if axis == 'x':
        cross_positions = x_coords
        interp_funcs = build_interp_funcs_for_axis(depth, 'x')
        west = x_coords[0];
        east = x_coords[-1]
        line_length = abs(y_coords[-1] - y_coords[0])
        along_count = len(interp_funcs)
    else:
        cross_positions = y_coords
        interp_funcs = build_interp_funcs_for_axis(depth, 'y')
        west = y_coords[0];
        east = y_coords[-1]
        line_length = abs(x_coords[-1] - x_coords[0])
        along_count = len(interp_funcs)

    def find_first_center():
        a = west;
        b = east
        for _ in range(80):
            m = 0.5 * (a + b)
            Lm, Rm, _, _ = global_left_right_for_center_simple(m, interp_funcs, axis)
            val = Lm - west
            if abs(val) < 1e-6:
                return m
            La, Ra, _, _ = global_left_right_for_center_simple(a, interp_funcs, axis)
            if (La - west) * val <= 0:
                b = m
            else:
                a = m
        return 0.5 * (a + b)

    # 简化函数：只取 min/max，用于定位第一条
    def global_left_right_for_center_simple(xc, interp_funcs_local, axis_local):
        Ls, Rs, Ds, Alphas = per_center_left_right_arrays(xc, interp_funcs_local, axis_local)
        return float(np.min(Ls)), float(np.max(Rs)), Ds, Alphas

    first = find_first_center()
    centers = [first]
    per_center_L = []
    per_center_R = []
    per_center_D = []
    per_center_alpha = []

    # 保存第一条的 per-sample L/R
    Ls0, Rs0, Ds0, Alphas0 = per_center_left_right_arrays(first, interp_funcs, axis)
    per_center_L.append(Ls0);
    per_center_R.append(Rs0);
    per_center_D.append(Ds0);
    per_center_alpha.append(Alphas0)

    for it in range(1000):
        xi = centers[-1]
        Lxi = per_center_L[-1];
        Rxi = per_center_R[-1]
        # 取当前带在沿线上的最右值作为 prev_right
        prev_right = float(np.max(Rxi))
        if prev_right >= east - 1e-6:
            break

        def f(x):
            Lx, Rx, _, _ = per_center_left_right_arrays(x, interp_funcs, axis)
            W = float(np.max(Rx) - np.min(Lx))
            return prev_right - float(np.min(Lx)) - eta * W

        a = xi + 1e-6;
        b = east
        fa = f(a);
        fb = f(b)
        if fa < 0:
            break
        if fa * fb > 0:
            approx = xi + (1 - eta) * (float(np.max(Rxi)) - float(np.min(Lxi)))
            centers.append(approx)
            Ls, Rs, Ds, Alphas = per_center_left_right_arrays(approx, interp_funcs, axis)
            per_center_L.append(Ls);
            per_center_R.append(Rs);
            per_center_D.append(Ds);
            per_center_alpha.append(Alphas)
            continue

        for _ in range(80):
            m = 0.5 * (a + b)
            fm = f(m)
            if abs(fm) < 1e-9 or (b - a) < 1e-6:
                a = b = m
                break
            if fa * fm <= 0:
                b = m;
                fb = fm
            else:
                a = m;
                fa = fm
        xnext = 0.5 * (a + b)
        centers.append(xnext)
        Ls, Rs, Ds, Alphas = per_center_left_right_arrays(xnext, interp_funcs, axis)
        per_center_L.append(Ls);
        per_center_R.append(Rs);
        per_center_D.append(Ds);
        per_center_alpha.append(Alphas)

    centers = np.array(centers)
    per_center_L = np.array(per_center_L)  # shape (n_centers, along_count)
    per_center_R = np.array(per_center_R)
    per_center_D = np.array(per_center_D)
    per_center_alpha = np.array(per_center_alpha)
    widths = per_center_R - per_center_L  # per-sample widths by center
    # compute aggregate width per center by max-min (conservative)
    agg_widths = np.max(per_center_R, axis=1) - np.min(per_center_L, axis=1)
    return centers, per_center_L, per_center_R, per_center_D, per_center_alpha, widths, agg_widths, line_length


# ---------- 逐点精确计算重叠率超限的沿线长度 ----------
def compute_exact_overlap_length(axis, centers, per_L, per_R, per_D, per_alpha, per_widths, agg_widths, line_length):
    """
    对于每一对相邻测线 (i, i+1)，对沿线每个样点 r：
      计算 overlap = max(0, min(R_i[r], R_j[r]) - max(L_i[r], L_j[r]))
      计算 overlap_rate = overlap / min(width_i[r], width_j[r])
      若 overlap_rate > overlap_threshold，则把该沿样点对应的沿线段长度累加（沿线采样间距 = line_length / (m-1)）
    最后返回所有对的累加和（单位：米）
    """
    n_centers = per_L.shape[0]
    m = per_L.shape[1]  # 沿线样点数
    # 沿线样点的间距（米）
    if m <= 1:
        along_spacing = line_length
    else:
        along_spacing = line_length / (m - 1)

    total_excess_length = 0.0
    # 对每一对相邻测线
    for i in range(n_centers - 1):
        L_i = per_L[i];
        R_i = per_R[i];
        W_i = per_widths[i]
        L_j = per_L[i + 1];
        R_j = per_R[i + 1];
        W_j = per_widths[i + 1]
        # 对沿线每个样点计算
        for r in range(m):
            w_i_r = W_i[r]
            w_j_r = W_j[r]
            # 若任一宽为零或非常小，跳过
            if (w_i_r <= 1e-6) or (w_j_r <= 1e-6):
                continue
            overlap = max(0.0, min(R_i[r], R_j[r]) - max(L_i[r], L_j[r]))
            if overlap <= 1e-9:
                continue
            overlap_rate = overlap / min(w_i_r, w_j_r)  # 以较窄带为基准（保守）
            if overlap_rate > overlap_threshold:
                total_excess_length += along_spacing
    return total_excess_length


# ---------- 主流程 ----------

# A: 测线走南北（中心沿 x）
centers_x, perLx, perRx, perDx, perAx, perWx, aggWx, line_len_x = iterative_centers_with_per_sample('x', eta_target)
metrics_x = {
    'n_centers': len(centers_x),
    'line_length': line_len_x,
    'total_survey_length': len(centers_x) * line_len_x
}
exact_overlap_x = compute_exact_overlap_length('x', centers_x, perLx, perRx, perDx, perAx, perWx, aggWx, line_len_x)


# 覆盖与漏测（用之前的粗方法）
def quick_coverage_stats(axis, perL, perR):
    cover_count = np.zeros_like(depth, dtype=int)
    if axis == 'x':
        for Ls, Rs in zip(perL, perR):
            col_l = np.searchsorted(x_coords, np.min(Ls), side='left')
            col_r = np.searchsorted(x_coords, np.max(Rs), side='right') - 1
            col_l = max(0, col_l);
            col_r = min(ncols - 1, col_r)
            if col_l <= col_r:
                cover_count[:, col_l:col_r + 1] += 1
    else:
        for Ls, Rs in zip(perL, perR):
            row_l = np.searchsorted(y_coords, np.min(Ls), side='left')
            row_r = np.searchsorted(y_coords, np.max(Rs), side='right') - 1
            row_l = max(0, row_l);
            row_r = min(nrows - 1, row_r)
            if row_l <= row_r:
                cover_count[row_l:row_r + 1, :] += 1
    covered_points = np.count_nonzero((cover_count > 0) & (~np.isnan(depth)))
    total_valid = np.count_nonzero(~np.isnan(depth))
    miss_fraction = 1.0 - covered_points / total_valid
    return miss_fraction


miss_x = quick_coverage_stats('x', perLx, perRx)

# B: 测线走东西（中心沿 y）
centers_y, perLy, perRy, perDy, perAy, perWy, aggWy, line_len_y = iterative_centers_with_per_sample('y', eta_target)
metrics_y = {
    'n_centers': len(centers_y),
    'line_length': line_len_y,
    'total_survey_length': len(centers_y) * line_len_y
}
exact_overlap_y = compute_exact_overlap_length('y', centers_y, perLy, perRy, perDy, perAy, perWy, aggWy, line_len_y)
miss_y = quick_coverage_stats('y', perLy, perRy)

# ---------- 输出 ----------
print("\n方案 A（测线走南北，中心沿 x）:")
print(
    f"  测线条数 = {metrics_x['n_centers']}, 每条长度 = {metrics_x['line_length']:.1f} m, 总长度 = {metrics_x['total_survey_length']:.1f} m")
print(f"  漏测率（近似） = {miss_x * 100:.3f} %")
print(f"  精确计算：重叠率 > {overlap_threshold * 100:.0f}% 的沿线累计长度 = {exact_overlap_x:.1f} m")

print("\n方案 B（测线走东西，中心沿 y）:")
print(
    f"  测线条数 = {metrics_y['n_centers']}, 每条长度 = {metrics_y['line_length']:.1f} m, 总长度 = {metrics_y['total_survey_length']:.1f} m")
print(f"  漏测率 = {miss_y * 100:.3f} %")
print(f"  精确计算：重叠率 > {overlap_threshold * 100:.0f}% 的沿线累计长度 = {exact_overlap_y:.1f} m")
print(f"\n按照总长度最短原则:")
if metrics_x['total_survey_length'] <= metrics_y['total_survey_length']:
    print("推荐方案：A（测线走南北）")
else:
    print("推荐方案：B（测线走东西）")


# === 图表生成函数 ===

def plot_depth_heatmap(depth_matrix, x_coords, y_coords):
    """生成深度热力图：显示海底地形"""
    plt.figure(figsize=(14, 10))

    # 创建掩码处理NaN值
    masked_depth = np.ma.masked_invalid(depth_matrix)

    # 创建热力图
    im = plt.imshow(masked_depth, cmap='viridis', aspect='auto',
                    extent=[x_coords[0], x_coords[-1], y_coords[0], y_coords[-1]],
                    origin='lower')

    # 添加颜色条
    cbar = plt.colorbar(im, ax=plt.gca(), shrink=0.8)
    cbar.set_label('Depth (m)', fontsize=12)

    # 设置标题和标签
    plt.title('Seabed Depth Heatmap', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('East-West Distance (m)', fontsize=12)
    plt.ylabel('North-South Distance (m)', fontsize=12)

    # 添加网格线
    plt.grid(True, alpha=0.3)

    # 添加统计信息
    valid_depths = depth_matrix[~np.isnan(depth_matrix)]
    stats_text = f'Depth Range: {valid_depths.min():.1f} - {valid_depths.max():.1f} m\n'
    stats_text += f'Mean Depth: {valid_depths.mean():.1f} m\n'
    stats_text += f'Terrain Relief: {valid_depths.max() - valid_depths.min():.1f} m\n'
    stats_text += f'Data Points: {len(valid_depths)}'

    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.show()


# 生成深度热力图
plot_depth_heatmap(depth, x_coords, y_coords)


def plot_survey_lines(depth_matrix, centers_x, centers_y, perLx, perRx, perLy, perRy, x_coords, y_coords):
    """生成测线布局图：显示设计的测线位置"""
    plt.figure(figsize=(16, 12))

    # 绘制深度背景（灰度）
    masked_depth = np.ma.masked_invalid(depth_matrix)
    plt.imshow(masked_depth, cmap='gray', aspect='auto',
               extent=[x_coords[0], x_coords[-1], y_coords[0], y_coords[-1]],
               origin='lower', alpha=0.7)

    # 绘制南北方向测线（中心沿x轴）
    for i, center_x in enumerate(centers_x):
        # 绘制测线中心线
        plt.axvline(x=center_x, color='red', linewidth=2, alpha=0.8,
                    label='North-South Lines' if i == 0 else "")

        # 绘制覆盖范围
        Ls = perLx[i]
        Rs = perRx[i]
        plt.fill_betweenx(y_coords, Ls, Rs, alpha=0.2, color='red')

        # 添加测线编号
        plt.text(center_x, y_coords[-1] * 0.95, f'Y{i + 1}',
                 rotation=90, fontsize=8, color='red', fontweight='bold',
                 ha='center', va='top')

    # 绘制东西方向测线（中心沿y轴）
    for i, center_y in enumerate(centers_y):
        # 绘制测线中心线
        plt.axhline(y=center_y, color='blue', linewidth=2, alpha=0.8,
                    label='East-West Lines' if i == 0 else "")

        # 绘制覆盖范围
        Ls = perLy[i]
        Rs = perRy[i]
        plt.fill_between(x_coords, Ls, Rs, alpha=0.2, color='blue')

        # 添加测线编号
        plt.text(x_coords[-1] * 0.95, center_y, f'X{i + 1}',
                 fontsize=8, color='blue', fontweight='bold',
                 ha='right', va='center')

    # 设置标题和标签
    plt.title('Multibeam Survey Lines Layout', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('East-West Distance (m)', fontsize=12)
    plt.ylabel('North-South Distance (m)', fontsize=12)

    # 添加图例
    plt.legend(loc='upper right', fontsize=10)

    # 添加统计信息
    stats_text = f'North-South Lines: {len(centers_x)}\n'
    stats_text += f'East-West Lines: {len(centers_y)}\n'
    stats_text += f'Total Lines: {len(centers_x) + len(centers_y)}\n'
    stats_text += f'NS Total Length: {metrics_x["total_survey_length"]:.1f} m\n'
    stats_text += f'EW Total Length: {metrics_y["total_survey_length"]:.1f} m'

    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


# 生成测线布局图
plot_survey_lines(depth, centers_x, centers_y, perLx, perRx, perLy, perRy, x_coords, y_coords)


def plot_3d_terrain(depth_matrix, x_coords, y_coords):
    """生成3D地形图：立体展示海底地形"""
    fig = plt.figure(figsize=(18, 14))
    ax = fig.add_subplot(111, projection='3d')

    # 创建坐标网格
    X, Y = np.meshgrid(x_coords, y_coords)

    # 处理NaN值
    Z = depth_matrix.copy()
    Z[np.isnan(Z)] = np.nanmin(Z)  # 将NaN替换为最小值

    # 创建3D表面图
    surf = ax.plot_surface(X, Y, Z, cmap='viridis',
                           linewidth=0, antialiased=True, alpha=0.8)

    # 添加等高线
    contour = ax.contour(X, Y, Z, levels=10, colors='black', alpha=0.5, linewidths=0.5)

    # 设置视角和标签
    ax.view_init(elev=30, azim=45)
    ax.set_xlabel('East-West Distance (m)', fontsize=12)
    ax.set_ylabel('North-South Distance (m)', fontsize=12)
    ax.set_zlabel('Depth (m)', fontsize=12)
    ax.set_title('3D Seabed Terrain', fontsize=16, fontweight='bold', pad=20)

    # 添加颜色条
    cbar = fig.colorbar(surf, ax=ax, shrink=0.8, aspect=20)
    cbar.set_label('Depth (m)', fontsize=12)

    # 添加统计信息
    valid_depths = depth_matrix[~np.isnan(depth_matrix)]
    stats_text = f'Depth Range: {valid_depths.min():.1f} - {valid_depths.max():.1f} m\n'
    stats_text += f'Mean Depth: {valid_depths.mean():.1f} m\n'
    stats_text += f'Terrain Relief: {valid_depths.max() - valid_depths.min():.1f} m\n'
    stats_text += f'Std Dev: {valid_depths.std():.1f} m'

    ax.text2D(0.02, 0.98, stats_text, transform=ax.transAxes,
              fontsize=10, verticalalignment='top',
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.show()


# 生成3D地形图
plot_3d_terrain(depth, x_coords, y_coords)



