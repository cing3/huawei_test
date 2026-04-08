import sys
import numpy as np

# ==========================================
# 计算几何与 NFP (No-Fit Polygon) 核心算法库
# ==========================================

def is_ccw(poly):
    """判断多边形是否为逆时针方向 (基于多边形符号面积)"""
    x = poly[:, 0]
    y = poly[:, 1]
    area = np.sum(x * np.roll(y, -1) - y * np.roll(x, -1))
    return area > 0

def clean_polygon(poly):
    """清理多边形的共线点和极近的重叠点，防止三角化失败"""
    if len(poly) < 3: return np.array(poly)
    cleaned = [poly[0]]
    for i in range(1, len(poly)):
        if np.hypot(poly[i][0] - cleaned[-1][0], poly[i][1] - cleaned[-1][1]) > 1e-7:
            cleaned.append(poly[i])
    if np.hypot(cleaned[-1][0] - cleaned[0][0], cleaned[-1][1] - cleaned[0][1]) < 1e-7:
        cleaned.pop()
    
    if len(cleaned) < 3: return np.array(cleaned)
    final_poly = []
    n = len(cleaned)
    for i in range(n):
        p_prev = cleaned[i - 1]
        p_curr = cleaned[i]
        p_next = cleaned[(i + 1) % n]
        cross = (p_curr[0] - p_prev[0]) * (p_next[1] - p_curr[1]) - (p_curr[1] - p_prev[1]) * (p_next[0] - p_curr[0])
        if abs(cross) > 1e-7:
            final_poly.append(p_curr)
            
    if len(final_poly) < 3: return np.array(cleaned)
    return np.array(final_poly)

def point_in_triangle(p, a, b, c):
    """判断点p是否在三角形abc内部或边界上"""
    v0, v1, v2 = c - a, b - a, p - a
    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)
    denom = dot00 * dot11 - dot01 * dot01
    if abs(denom) < 1e-14: return False
    invDenom = 1.0 / denom
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom
    return (u >= -1e-9) and (v >= -1e-9) and (u + v <= 1.0 + 1e-9)

def triangulate(poly):
    """耳切法 (Ear Clipping) 将简单多边形分解为多个凸三角形"""
    poly = clean_polygon(poly)
    if not is_ccw(poly):
        poly = poly[::-1]
        
    vertices = list(poly)
    triangles = []
    indices = list(range(len(vertices)))
    
    iterations = 0
    max_iters = len(vertices) * 10
    
    while len(indices) > 3 and iterations < max_iters:
        iterations += 1
        n_current = len(indices)
        ear_found = False
        for i in range(n_current):
            prev_i = indices[i - 1]
            curr_i = indices[i]
            next_i = indices[(i + 1) % n_current]
            
            p_prev, p_curr, p_next = vertices[prev_i], vertices[curr_i], vertices[next_i]
            
            # 判断顶点是否凸出
            cross = (p_curr[0] - p_prev[0]) * (p_next[1] - p_curr[1]) - (p_curr[1] - p_prev[1]) * (p_next[0] - p_curr[0])
            if cross <= 1e-9: continue
                
            # 判断其他顶点是否在此耳朵内部
            is_ear = True
            for j in range(n_current):
                test_i = indices[j]
                if test_i in (prev_i, curr_i, next_i): continue
                if point_in_triangle(vertices[test_i], p_prev, p_curr, p_next):
                    is_ear = False
                    break
                    
            if is_ear:
                triangles.append(np.array([p_prev, p_curr, p_next]))
                indices.pop(i)
                ear_found = True
                break
                
        if not ear_found:
            # 防死循环后备策略
            triangles.append(np.array([vertices[indices[-1]], vertices[indices[0]], vertices[indices[1]]]))
            indices.pop(0)

    if len(indices) == 3:
        triangles.append(np.array([vertices[indices[0]], vertices[indices[1]], vertices[indices[2]]]))
        
    return triangles

def convex_hull(points):
    """计算二维点集的凸包"""
    points = np.unique(points, axis=0)
    if len(points) <= 2: return points
    idx = np.lexsort((points[:, 1], points[:, 0]))
    points = points[idx]
    
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 1e-9: lower.pop()
        lower.append(p)
        
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 1e-9: upper.pop()
        upper.append(p)
        
    return np.array(lower[:-1] + upper[:-1])

def point_strictly_in_convex_poly(p, poly, eps=1e-7):
    """判断点是否严格在凸多边形内部（距离边界大于 eps）"""
    x, y = p
    n = len(poly)
    for i in range(n):
        p1, p2 = poly[i], poly[(i + 1) % n]
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        length = np.hypot(dx, dy)
        if length < 1e-12: continue
        cross = (dx * (y - p1[1]) - dy * (x - p1[0])) / length
        if cross <= eps:
            return False
    return True

def point_in_convex_poly(p, poly, eps=0.0):
    """判断点是否在凸多边形内部或边界上"""
    x, y = p
    n = len(poly)
    for i in range(n):
        p1, p2 = poly[i], poly[(i + 1) % n]
        cross = (p2[0] - p1[0]) * (y - p1[1]) - (p2[1] - p1[1]) * (x - p1[0])
        if cross < eps:
            return False
    return True

def extract_union_boundary_chunked(polygons):
    """基于向量化的线段相交测试，提取凸多边形并集的真实外边界"""
    edges = []
    for poly in polygons:
        n = len(poly)
        for i in range(n):
            edges.append((poly[i], poly[(i+1)%n]))
            
    E = np.array(edges)
    n_edges = len(E)
    if n_edges == 0: return np.array([])
    
    P1, P2 = E[:, 0, :], E[:, 1, :]
    D = P2 - P1
    
    t_intersections = [[] for _ in range(n_edges)]
    CHUNK_SIZE = 2000  # 分块处理避免内存爆炸
    
    for i in range(0, n_edges, CHUNK_SIZE):
        i_end = min(i + CHUNK_SIZE, n_edges)
        Di, P1i = D[i:i_end, np.newaxis, :], P1[i:i_end, np.newaxis, :]
        
        for j in range(0, n_edges, CHUNK_SIZE):
            j_end = min(j + CHUNK_SIZE, n_edges)
            Dj, P1j = D[np.newaxis, j:j_end, :], P1[np.newaxis, j:j_end, :]
            
            det = Di[..., 0] * Dj[..., 1] - Di[..., 1] * Dj[..., 0]
            DP = P1j - P1i
            
            num_t = DP[..., 0] * Dj[..., 1] - DP[..., 1] * Dj[..., 0]
            num_u = DP[..., 0] * Di[..., 1] - DP[..., 1] * Di[..., 0]
            
            valid = np.abs(det) > 1e-9
            t, u = np.zeros_like(det), np.zeros_like(det)
            t[valid], u[valid] = num_t[valid] / det[valid], num_u[valid] / det[valid]
            
            intersect = valid & (t >= -1e-9) & (t <= 1.0 + 1e-9) & (u >= -1e-9) & (u <= 1.0 + 1e-9)
            
            idx_i, idx_j = np.nonzero(intersect)
            for k in range(len(idx_i)):
                global_i = i + idx_i[k]
                t_intersections[global_i].append(t[idx_i[k], idx_j[k]])

    poly_boxes = np.array([
        [np.min(p[:,0]), np.max(p[:,0]), np.min(p[:,1]), np.max(p[:,1])] 
        for p in polygons
    ])
    
    exposed_segments = []
    
    for i in range(n_edges):
        t_vals = np.array(t_intersections[i])
        t_vals = np.clip(t_vals, 0.0, 1.0)
        t_vals = np.concatenate(([0.0, 1.0], t_vals))
        t_vals.sort()
        t_vals = t_vals[np.append([True], np.diff(t_vals) > 1e-7)]
        
        for k in range(len(t_vals) - 1):
            t_mid = (t_vals[k] + t_vals[k+1]) / 2.0
            midpoint = P1[i] + t_mid * D[i]
            
            x, y = midpoint
            possible_polys = np.nonzero(
                (x >= poly_boxes[:, 0]) & (x <= poly_boxes[:, 1]) &
                (y >= poly_boxes[:, 2]) & (y <= poly_boxes[:, 3])
            )[0]
            
            is_inside = False
            for idx in possible_polys:
                if point_strictly_in_convex_poly(midpoint, polygons[idx], eps=1e-7):
                    is_inside = True
                    break
            
            if not is_inside:
                p_start = P1[i] + t_vals[k] * D[i]
                p_end = P1[i] + t_vals[k+1] * D[i]
                if np.hypot(p_end[0]-p_start[0], p_end[1]-p_start[1]) > 1e-9:
                    exposed_segments.append((p_start, p_end))
                
    return np.array(exposed_segments)

# ==========================================
# 核心查询运算：向量化查询平移情况
# ==========================================

def calculate_mtv_query_vectorized(tx, ty, polygons, boundary_segments, poly_boxes):
    """计算单个平移向量的 MTV (向量化处理距离以防止超时)"""
    T = np.array([tx, ty], dtype=np.float64)
    
    # 1. 判断此时 A 和移动后的 B 是否重叠 (即 T 是否在 Minkowski Union 中)
    is_inside = False
    possible_polys = np.nonzero(
        (tx >= poly_boxes[:, 0]) & (tx <= poly_boxes[:, 1]) &
        (ty >= poly_boxes[:, 2]) & (ty <= poly_boxes[:, 3])
    )[0]
    
    for idx in possible_polys:
        if point_in_convex_poly(T, polygons[idx], eps=-1e-9):
            is_inside = True
            break
            
    # 无重叠直接返回 0,0
    if not is_inside:
        return 0.0, 0.0
        
    # 2. 存在重叠，寻找从 T 到并集外边界的最短逃逸向量 (高度向量化处理)
    if len(boundary_segments) == 0: return 0.0, 0.0
        
    P1 = boundary_segments[:, 0, :]
    P2 = boundary_segments[:, 1, :]
    
    AB = P2 - P1
    AP = T - P1
    AB_sq = np.sum(AB * AB, axis=1)
    
    valid = AB_sq > 1e-14
    t = np.zeros(len(boundary_segments))
    t[valid] = np.sum(AP[valid] * AB[valid], axis=1) / AB_sq[valid]
    t = np.clip(t, 0.0, 1.0)
    
    closest_pts = P1 + t[:, np.newaxis] * AB
    vecs = closest_pts - T
    dists_sq = np.sum(vecs * vecs, axis=1)
    
    best_idx = np.argmin(dists_sq)
    best_vec = vecs[best_idx]
    
    return best_vec[0], best_vec[1]

# ==========================================
# I/O 交互与主流程
# ==========================================

def main():
    # 读取第一行: 顶点数
    line1 = sys.stdin.readline().strip()
    if not line1: return
    na, nb = map(int, line1.split())
    
    # 读取多边形 A 和 B 的坐标
    coords_a_raw = list(map(float, sys.stdin.readline().split()))
    coords_b_raw = list(map(float, sys.stdin.readline().split()))
    
    poly_a = np.array(coords_a_raw, dtype=np.float64).reshape(-1, 2)
    poly_b = np.array(coords_b_raw, dtype=np.float64).reshape(-1, 2)
    
    # 读取第一阶段的 OK 信号
    sys.stdin.readline().strip()
    
    # === [关键预处理] 计算闵可夫斯基差并集边界 ===
    tri_a = triangulate(poly_a)
    tri_b = triangulate(poly_b)
    
    global_polys = []
    for ta in tri_a:
        for tb in tri_b:
            # 计算两三角形的 Minkowski Difference: ta (+) (-tb)
            pts = []
            for pa in ta:
                for pb in tb:
                    pts.append(pa - pb)
            hull = convex_hull(np.array(pts))
            global_polys.append(hull)
            
    # 提取完整的真实外边界
    global_boundary = extract_union_boundary_chunked(global_polys)
    
    # 预计算 Bounding Box 加速判断
    global_boxes = np.array([
        [np.min(poly[:,0]), np.max(poly[:,0]), np.min(poly[:,1]), np.max(poly[:,1])] 
        for poly in global_polys
    ])
    # ==========================================
    
    # 预处理完毕，向判题器发送 OK
    sys.stdout.write("OK\n")
    sys.stdout.flush()
    
    # 读取 K 个平移测试向量
    line_k = sys.stdin.readline().strip()
    if not line_k: return
    k = int(line_k)
    
    translations = []
    for _ in range(k):
        tx, ty = map(float, sys.stdin.readline().split())
        translations.append((tx, ty))
        
    # 读取数据完毕的 OK 信号
    sys.stdin.readline().strip()
    
    # 输出结果数量
    sys.stdout.write(f"{k}\n")
    
    # 进行高速判定
    for tx, ty in translations:
        mtv_x, mtv_y = calculate_mtv_query_vectorized(tx, ty, global_polys, global_boundary, global_boxes)
        sys.stdout.write(f"{mtv_x:.5f} {mtv_y:.5f}\n")
        
    # 结果输出完毕
    sys.stdout.write("OK\n")
    sys.stdout.flush()

    # 规范退出：防挂起设计，等待管道关闭(EOF)
    while True:
        if not sys.stdin.readline():
            break

if __name__ == "__main__":
    main()