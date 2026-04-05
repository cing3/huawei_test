import sys
import numpy as np

EPS = 1e-9
EPS_STRICT = 1e-10
DIR_KEY_SCALE = 1e10
_CTX_CACHE = {}


def _cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _signed_area(poly):
    x = poly[:, 0]
    y = poly[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - y * np.roll(x, -1))


def _point_on_segment(p, a, b, eps=EPS):
    ab = b - a
    ap = p - a
    if abs(_cross(ab, ap)) > eps:
        return False
    dot = np.dot(ap, ab)
    if dot < -eps:
        return False
    if dot - np.dot(ab, ab) > eps:
        return False
    return True


def _point_in_triangle_inclusive(p, a, b, c, eps=EPS):
    v0 = b - a
    v1 = c - b
    v2 = a - c
    cp0 = _cross(v0, p - a)
    cp1 = _cross(v1, p - b)
    cp2 = _cross(v2, p - c)

    has_pos = (cp0 > eps) or (cp1 > eps) or (cp2 > eps)
    has_neg = (cp0 < -eps) or (cp1 < -eps) or (cp2 < -eps)
    return not (has_pos and has_neg)


def _sanitize_polygon(poly, eps=EPS):
    pts = [np.array(p, dtype=np.float64) for p in poly]
    if len(pts) >= 2 and np.linalg.norm(pts[0] - pts[-1]) <= eps:
        pts.pop()

    if not pts:
        return np.empty((0, 2), dtype=np.float64)

    cleaned = [pts[0]]
    for p in pts[1:]:
        if np.linalg.norm(p - cleaned[-1]) > eps:
            cleaned.append(p)

    if len(cleaned) >= 2 and np.linalg.norm(cleaned[0] - cleaned[-1]) <= eps:
        cleaned.pop()

    if len(cleaned) < 3:
        return np.array(cleaned, dtype=np.float64)

    changed = True
    while changed and len(cleaned) >= 3:
        changed = False
        n = len(cleaned)
        keep = []
        for i in range(n):
            prev = cleaned[(i - 1) % n]
            cur = cleaned[i]
            nxt = cleaned[(i + 1) % n]
            e1 = cur - prev
            e2 = nxt - cur
            if np.linalg.norm(e1) <= eps or np.linalg.norm(e2) <= eps:
                changed = True
                continue
            if abs(_cross(e1, e2)) <= eps and np.dot(e1, e2) >= 0:
                changed = True
                continue
            keep.append(cur)
        if len(keep) < 3:
            break
        cleaned = keep

    return np.array(cleaned, dtype=np.float64)


def _triangulate_ear_clipping(poly):
    poly = _sanitize_polygon(poly)
    n = len(poly)
    if n < 3:
        return []

    area = _signed_area(poly)
    if abs(area) <= EPS:
        return []

    ccw = area > 0
    idx = list(range(n))
    triangles = []

    def is_convex(i0, i1, i2):
        a = poly[i0]
        b = poly[i1]
        c = poly[i2]
        val = _cross(b - a, c - b)
        return val > EPS if ccw else val < -EPS

    max_iter = n * n
    it = 0

    while len(idx) > 3 and it < max_iter:
        it += 1
        ear_found = False
        m = len(idx)

        for k in range(m):
            i_prev = idx[(k - 1) % m]
            i_curr = idx[k]
            i_next = idx[(k + 1) % m]

            if not is_convex(i_prev, i_curr, i_next):
                continue

            a = poly[i_prev]
            b = poly[i_curr]
            c = poly[i_next]

            bad = False
            for j in idx:
                if j == i_prev or j == i_curr or j == i_next:
                    continue
                p = poly[j]
                if _point_in_triangle_inclusive(p, a, b, c, eps=EPS):
                    bad = True
                    break

            if bad:
                continue

            triangles.append(np.array([a, b, c], dtype=np.float64))
            del idx[k]
            ear_found = True
            break

        if not ear_found:
            break

    if len(idx) == 3:
        triangles.append(np.array([poly[idx[0]], poly[idx[1]], poly[idx[2]]], dtype=np.float64))

    return triangles


def _edge_normals(convex_poly):
    normals = []
    m = len(convex_poly)
    for i in range(m):
        e = convex_poly[(i + 1) % m] - convex_poly[i]
        l = np.hypot(e[0], e[1])
        if l <= EPS:
            continue
        n = np.array([-e[1], e[0]], dtype=np.float64) / l
        normals.append(n)
    return normals


def _project(poly, axis):
    vals = poly @ axis
    return np.min(vals), np.max(vals)


def _poly_bbox(poly):
    if len(poly) == 0:
        return np.array([0.0, 0.0]), np.array([0.0, 0.0])
    return np.min(poly, axis=0), np.max(poly, axis=0)


def _bbox_overlap(min_a, max_a, min_b, max_b, eps=EPS):
    if max_a[0] < min_b[0] - eps or max_b[0] < min_a[0] - eps:
        return False
    if max_a[1] < min_b[1] - eps or max_b[1] < min_a[1] - eps:
        return False
    return True


def _convex_overlap_strict(poly1, poly2):
    axes = _edge_normals(poly1) + _edge_normals(poly2)
    if not axes:
        return False

    for ax in axes:
        mn1, mx1 = _project(poly1, ax)
        mn2, mx2 = _project(poly2, ax)
        if not (mx1 > mn2 + EPS_STRICT and mx2 > mn1 + EPS_STRICT):
            return False
    return True


def _polygons_overlap_strict_by_tris(tris_a, tris_b):
    for ta in tris_a:
        for tb in tris_b:
            if _convex_overlap_strict(ta, tb):
                return True
    return False


def _interval_strict_overlap_along_dir(conv_a, conv_b, d):
    axes = _edge_normals(conv_a) + _edge_normals(conv_b)
    if not axes:
        return None

    low = -np.inf
    high = np.inf

    for ax in axes:
        pmin, pmax = _project(conv_a, ax)
        qmin, qmax = _project(conv_b, ax)

        a = float(np.dot(d, ax))
        L = (pmin - qmax) + EPS_STRICT
        U = (pmax - qmin) - EPS_STRICT

        if L >= U:
            return None

        if abs(a) <= EPS:
            if not (0.0 > L and 0.0 < U):
                return None
            continue

        t1 = L / a
        t2 = U / a
        lo = t1 if t1 < t2 else t2
        hi = t2 if t2 > t1 else t1

        if lo > low:
            low = lo
        if hi < high:
            high = hi
        if low >= high:
            return None

    if high <= 0.0:
        return None
    if low < 0.0:
        low = 0.0

    if low >= high:
        return None
    return (low, high)


def _min_non_overlap_t_along_dir(tris_a, tris_b_moved, d, current_best=np.inf):
    intervals = []

    for ta in tris_a:
        for tb in tris_b_moved:
            seg = _interval_strict_overlap_along_dir(ta, tb, d)
            if seg is None:
                continue
            s, e = seg

            if s >= current_best:
                continue
            if e > current_best:
                e = current_best

            if s < e:
                intervals.append((s, e))

    if not intervals:
        return np.inf

    intervals.sort(key=lambda x: (x[0], x[1]))

    cur = 0.0
    for s, e in intervals:
        if e <= cur + EPS:
            continue

        if cur <= s + EPS:
            return cur

        if cur < e - EPS:
            cur = e
            if cur >= current_best - EPS:
                return np.inf

    return cur


def _add_direction(vec, keys, dirs):
    n = np.hypot(vec[0], vec[1])
    if n <= EPS:
        return
    u = vec / n
    kx = int(round(u[0] * DIR_KEY_SCALE))
    ky = int(round(u[1] * DIR_KEY_SCALE))
    key = (kx, ky)
    if key not in keys:
        keys.add(key)
        dirs.append(u)


def _build_edge_normal_dirs(poly_a, moved_b):
    keys = set()
    dirs = []

    def add_edge_normals(poly):
        m = len(poly)
        for i in range(m):
            e = poly[(i + 1) % m] - poly[i]
            n = np.array([-e[1], e[0]], dtype=np.float64)
            _add_direction(n, keys, dirs)
            _add_direction(-n, keys, dirs)

    add_edge_normals(poly_a)
    add_edge_normals(moved_b)

    if not dirs:
        dirs = [np.array([1.0, 0.0], dtype=np.float64), np.array([0.0, 1.0], dtype=np.float64)]
    return dirs

def _build_candidate_dirs(poly_a, moved_b):
    keys = set()
    out = []

    # 1) 所有边法向（双向）
    for d in _build_edge_normal_dirs(poly_a, moved_b):
        _add_direction(d, keys, out)

    # 2) 所有顶点差方向（双向）
    for va in poly_a:
        for vb in moved_b:
            diff = va - vb
            _add_direction(diff, keys, out)
            _add_direction(-diff, keys, out)

    if not out:
        out = [np.array([1.0, 0.0], dtype=np.float64), np.array([0.0, 1.0], dtype=np.float64)]
    return out


def _is_convex_polygon(poly):
    n = len(poly)
    if n < 3:
        return False
    sign = 0
    for i in range(n):
        a = poly[i]
        b = poly[(i + 1) % n]
        c = poly[(i + 2) % n]
        cr = _cross(b - a, c - b)
        if abs(cr) <= EPS:
            continue
        s = 1 if cr > 0 else -1
        if sign == 0:
            sign = s
        elif sign != s:
            return False
    return True

def _get_context(poly_a, poly_b):
    key = (id(poly_a), id(poly_b))
    ctx = _CTX_CACHE.get(key)
    if ctx is not None:
        return ctx

    a_clean = _sanitize_polygon(poly_a)
    b_clean = _sanitize_polygon(poly_b)

    convex_a = _is_convex_polygon(a_clean)
    convex_b = _is_convex_polygon(b_clean)

    min_a, max_a = _poly_bbox(a_clean)
    min_b, max_b = _poly_bbox(b_clean)

    ctx = {
        "a_clean": a_clean,
        "b_clean": b_clean,
        "convex_a": convex_a,
        "convex_b": convex_b,
        "tris_a": _triangulate_ear_clipping(a_clean),
        "tris_b_base": _triangulate_ear_clipping(b_clean),
        "bbox_min_a": min_a,
        "bbox_max_a": max_a,
        "bbox_min_b": min_b,
        "bbox_max_b": max_b,
    }
    _CTX_CACHE[key] = ctx
    return ctx

def calculate_mtv(poly_a, poly_b, translation):
    shift = np.array(translation, dtype=np.float64)
    ctx = _get_context(poly_a, poly_b)

    a_clean = ctx["a_clean"]
    b_clean = ctx["b_clean"]
    moved_b = b_clean + shift

    if len(a_clean) < 3 or len(b_clean) < 3:
        return 0.0, 0.0

    min_b = ctx["bbox_min_b"] + shift
    max_b = ctx["bbox_max_b"] + shift
    if not _bbox_overlap(ctx["bbox_min_a"], ctx["bbox_max_a"], min_b, max_b):
        return 0.0, 0.0

    # 快路径：两者均为凸多边形（先判是否重叠，避免无重叠时构造大量候选方向）
    if ctx["convex_a"] and ctx["convex_b"]:
        if not _convex_overlap_strict(a_clean, moved_b):
            return 0.0, 0.0

        candidate_dirs = _build_candidate_dirs(a_clean, moved_b)
        best_t = np.inf
        best_v = np.array([0.0, 0.0], dtype=np.float64)
        for d in candidate_dirs:
            seg = _interval_strict_overlap_along_dir(a_clean, moved_b, d)
            if seg is None:
                continue
            _, e = seg
            t = e if e > 0.0 else 0.0
            if t < best_t:
                best_t = t
                best_v = d * t

        if not np.isfinite(best_t):
            return 0.0, 0.0
        return float(best_v[0]), float(best_v[1])

    # 通用路径：凹多边形（三角剖分并集）
    tris_a = ctx["tris_a"]
    tris_b_base = ctx["tris_b_base"]
    if not tris_a or not tris_b_base:
        return 0.0, 0.0

    tris_b = [tri + shift for tri in tris_b_base]

    if not _polygons_overlap_strict_by_tris(tris_a, tris_b):
        return 0.0, 0.0

    candidate_dirs = _build_candidate_dirs(a_clean, moved_b)
    best_t = np.inf
    best_v = np.array([0.0, 0.0], dtype=np.float64)

    for d in candidate_dirs:
        t = _min_non_overlap_t_along_dir(tris_a, tris_b, d, current_best=best_t)
        if not np.isfinite(t):
            continue
        if t < best_t:
            best_t = t
            best_v = d * t

    if not np.isfinite(best_t):
        return 0.0, 0.0

    return float(best_v[0]), float(best_v[1])


def main():
    line = sys.stdin.readline()
    if not line:
        return

    try:
        n1, n2 = map(int, line.split())
    except ValueError:
        return

    poly_a_list = []
    for _ in range(n1):
        x, y = map(float, sys.stdin.readline().split())
        poly_a_list.append((x, y))

    poly_b_list = []
    for _ in range(n2):
        x, y = map(float, sys.stdin.readline().split())
        poly_b_list.append((x, y))

    poly_a = np.array(poly_a_list, dtype=np.float64)
    poly_b = np.array(poly_b_list, dtype=np.float64)

    ok_response = sys.stdin.readline().strip()
    if ok_response != "OK":
        return

    print("OK")
    sys.stdout.flush()

    line_k = sys.stdin.readline().strip()
    if not line_k:
        return
    k = int(line_k)

    translations = []
    for _ in range(k):
        tx, ty = map(float, sys.stdin.readline().split())
        translations.append((tx, ty))

    ok_response = sys.stdin.readline().strip()
    if ok_response != "OK":
        return

    for t in translations:
        rx, ry = calculate_mtv(poly_a, poly_b, t)
        print(f"{rx:.5f} {ry:.5f}")
        sys.stdout.flush()

    print("OK")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
