import sys
import numpy as np

# ==========================================
# 健壮的 Token 提取器，无视多余换行、空格
# ==========================================
def get_tokens():
    for line in sys.stdin:
        for token in line.split():
            yield token

tokens = get_tokens()

def next_token():
    try:
        return next(tokens)
    except StopIteration:
        return None

def next_float():
    tok = next_token()
    if tok is None: return None
    return float(tok)

def next_int():
    tok = next_token()
    if tok is None: return None
    return int(tok)

def expect_ok():
    while True:
        tok = next_token()
        if tok is None: break
        if tok == "OK":
            break

# ==========================================
# 计算几何与 NFP (No-Fit Polygon) 核心算法库
# ==========================================

def is_ccw(poly):
    """判断多边形是否为逆时针方向"""
    x = poly[:, 0]
    y = poly[:, 1]
    area = np.sum(x * np.roll(y, -1) - y * np.roll(x, -1))
    return area > 0

def clean_polygon(poly):
    """清理多边形的共线点和极近的重叠点"""
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
            
            cross = (p_curr[0] - p_prev[0]) * (p_next[1] - p_curr[1]) - (p_curr[1] - p_prev[1]) * (p_next[0] - p_curr[0])
            if cross <= 1e-9: continue
                
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

def point_in_convex_poly(p, poly, margin=1e-9):
    """鲁棒的凸多边形内点判定"""
    x, y = p
    n = len(poly)
    for i in range(n):
        p1 = poly[i]
        p2 = poly[(i + 1) % n]
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        length = np.hypot(dx, dy)
        if length < 1e-12: continue
        cross = dx * (y - p1[1]) - dy * (x - p1[0])
        # cross / length 就是投影的法向距离
        if cross / length < -margin:
            return False
    return True

def is_point_in_NFP(P, global_polys, global_boxes, margin=1e-9):
    """判断点是否在 NFP 内或 NFP 的边界容差范围内"""
    possible = np.nonzero(
        (P[0] >= global_boxes[:, 0] - margin) & (P[0] <= global_boxes[:, 1] + margin) &
        (P[1] >= global_boxes[:, 2] - margin) & (P[1] <= global_boxes[:, 3] + margin)
    )[0]
    for idx in possible:
        if point_in_convex_poly(P, global_polys[idx], margin=margin):
            return True
    return False

def is_strictly_inside_NFP(T, global_polys, global_boxes):
    """采用 8 方向探测法判断 T 是否被严格包含在 NFP 内部"""
    angles = np.linspace(0, 2*np.pi, 8, endpoint=False)
    for a in angles:
        pt = T + 1e-8 * np.array([np.cos(a), np.sin(a)])
        if not is_point_in_NFP(pt, global_polys, global_boxes, margin=1e-10):
            return False
    return True

# ==========================================
# 核心查询：通过最近边界打分进行查询
# ==========================================
def calculate_mtv(tx, ty, segments, global_polys, global_boxes):
    T = np.array([tx, ty], dtype=np.float64)
    
    # 1. 严格过滤，如果不被严格包裹，MTV直接视为 (0, 0)
    if not is_strictly_inside_NFP(T, global_polys, global_boxes):
        return 0.0, 0.0

    if len(segments) == 0:
        return 0.0, 0.0

    P1 = segments[:, 0, :]
    P2 = segments[:, 1, :]
    
    AB = P2 - P1
    AP = T - P1
    AB_sq = np.sum(AB * AB, axis=1)
    
    valid = AB_sq > 1e-14
    t = np.zeros(len(segments))
    t[valid] = np.sum(AP[valid] * AB[valid], axis=1) / AB_sq[valid]
    t = np.clip(t, 0.0, 1.0)
    
    # 计算 T 到各候选线段的投影点距离
    closest_pts = P1 + t[:, np.newaxis] * AB
    vecs = closest_pts - T
    dists_sq = np.sum(vecs * vecs, axis=1)
    
    # 按照距离短的优先扫描，首次命中必然是最短逃逸向量
    sorted_indices = np.argsort(dists_sq)
    
    for idx in sorted_indices:
        dist_sq = dists_sq[idx]
        v = vecs[idx]
        dist = np.sqrt(dist_sq)
        if dist < 1e-11:
            # 自动跳过恰好踩在 NFP 内部隔板上的假目标
            continue
            
        d = v / dist
        # 沿着逃出方向迈出极其微小的一步，来确认外侧是否真的开阔
        P_test = closest_pts[idx] + 1e-7 * d
        
        if not is_point_in_NFP(P_test, global_polys, global_boxes, margin=1e-10):
            return v[0], v[1]
            
    return 0.0, 0.0

# ==========================================
# I/O 交互与主流程
# ==========================================

def main():
    na = next_int()
    if na is None: return
    nb = next_int()
    
    poly_a = np.zeros((na, 2), dtype=np.float64)
    for i in range(na):
        poly_a[i, 0] = next_float()
        poly_a[i, 1] = next_float()
        
    poly_b = np.zeros((nb, 2), dtype=np.float64)
    for i in range(nb):
        poly_b[i, 0] = next_float()
        poly_b[i, 1] = next_float()
        
    expect_ok()
    
    # === [关键预处理] 三角化与 Minkowski 差异分块 ===
    tri_a = triangulate(poly_a)
    tri_b = triangulate(poly_b)
    
    global_polys = []
    for ta in tri_a:
        for tb in tri_b:
            pts = []
            for pa in ta:
                for pb in tb:
                    pts.append(pa - pb)
            hull = convex_hull(np.array(pts))
            global_polys.append(hull)
            
    global_boxes = np.array([
        [np.min(poly[:,0]), np.max(poly[:,0]), np.min(poly[:,1]), np.max(poly[:,1])] 
        for poly in global_polys
    ])
    
    # === [关键预处理] 构建候选逃逸界限 ===
    segments = []
    # 策略 1: A 的边 与 B 的顶点 构成线段
    for i in range(na):
        a1, a2 = poly_a[i], poly_a[(i+1)%na]
        for j in range(nb):
            b1 = poly_b[j]
            segments.append((a1 - b1, a2 - b1))
            
    # 策略 2: A 的顶点 与 B 的边 构成线段
    for i in range(na):
        a1 = poly_a[i]
        for j in range(nb):
            b1, b2 = poly_b[j], poly_b[(j+1)%nb]
            segments.append((a1 - b1, a1 - b2))
            
    segments = np.array(segments)
    
    # 预处理完毕，向判题器发送 OK
    sys.stdout.write("OK\n")
    sys.stdout.flush()
    
    k = next_int()
    if k is None: return
    
    translations = np.zeros((k, 2), dtype=np.float64)
    for i in range(k):
        translations[i, 0] = next_float()
        translations[i, 1] = next_float()
        
    expect_ok()
    
    sys.stdout.write(f"{k}\n")
    
    # 高速向量化查表
    for i in range(k):
        tx, ty = translations[i, 0], translations[i, 1]
        mtv_x, mtv_y = calculate_mtv(tx, ty, segments, global_polys, global_boxes)
        sys.stdout.write(f"{mtv_x:.5f} {mtv_y:.5f}\n")
        
    sys.stdout.write("OK\n")
    sys.stdout.flush()

if __name__ == "__main__":
    main()