#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

const double EPS = 1e-7;

struct Point {
    double x, y;
    Point operator-(const Point& o) const { return {x - o.x, y - o.y}; }
    Point operator+(const Point& o) const { return {x + o.x, y + o.y}; }
    Point operator*(double scalar) const { return {x * scalar, y * scalar}; }
};

typedef vector<Point> Polygon;

double dist(const Point& a, const Point& b) {
    return std::hypot(a.x - b.x, a.y - b.y);
}

double cross_product(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

double dot(const Point& v1, const Point& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

double length(const Point& v) {
    return std::hypot(v.x, v.y);
}

double area(const Polygon& p) {
    double a = 0;
    int n = p.size();
    for (int i = 0; i < n; ++i) {
        a += p[i].x * p[(i + 1) % n].y - p[(i + 1) % n].x * p[i].y;
    }
    return a;
}

void ensure_ccw(Polygon& p) {
    if (area(p) < 0) std::reverse(p.begin(), p.end());
}

Polygon remove_collinear(const Polygon& p) {
    Polygon res;
    int n = p.size();
    for (int i = 0; i < n; ++i) {
        Point prev = p[(i - 1 + n) % n];
        Point curr = p[i];
        Point next = p[(i + 1) % n];
        if (abs(cross_product(prev, curr, next)) > EPS) {
            res.push_back(curr);
        }
    }
    return res;
}

bool point_in_triangle(const Point& pt, const Point& a, const Point& b, const Point& c) {
    double cp1 = cross_product(a, b, pt);
    double cp2 = cross_product(b, c, pt);
    double cp3 = cross_product(c, a, pt);
    bool has_neg = (cp1 < -EPS) || (cp2 < -EPS) || (cp3 < -EPS);
    bool has_pos = (cp1 > EPS) || (cp2 > EPS) || (cp3 > EPS);
    return !(has_neg && has_pos);
}

vector<Polygon> triangulate(Polygon p) {
    ensure_ccw(p);
    p = remove_collinear(p);
    vector<Polygon> tris;
    while (p.size() > 3) {
        int n = p.size();
        bool ear_found = false;
        for (int i = 0; i < n; ++i) {
            int prev = (i - 1 + n) % n;
            int next = (i + 1) % n;
            if (cross_product(p[prev], p[i], p[next]) > EPS) { 
                bool empty = true;
                for (int j = 0; j < n; ++j) {
                    if (j == prev || j == i || j == next) continue;
                    if (point_in_triangle(p[j], p[prev], p[i], p[next])) {
                        empty = false; break;
                    }
                }
                if (empty) {
                    tris.push_back({p[prev], p[i], p[next]});
                    p.erase(p.begin() + i);
                    ear_found = true;
                    break;
                }
            }
        }
        if (!ear_found) {
            tris.push_back({p[0], p[1], p[2]});
            p.erase(p.begin() + 1);
        }
    }
    if (p.size() == 3) tris.push_back(p);
    return tris;
}

bool try_merge(const Polygon& p1, const Polygon& p2, Polygon& merged) {
    int e1 = -1, e2 = -1;
    int n1 = p1.size(), n2 = p2.size();
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            if (dist(p1[i], p2[(j + 1) % n2]) < EPS && dist(p1[(i + 1) % n1], p2[j]) < EPS) {
                e1 = i; e2 = j; break;
            }
        }
        if (e1 != -1) break;
    }
    if (e1 == -1) return false;

    Polygon temp;
    for (int i = 1; i < n1; ++i) temp.push_back(p1[(e1 + i) % n1]);
    for (int j = 1; j < n2; ++j) temp.push_back(p2[(e2 + j) % n2]);

    temp = remove_collinear(temp);
    if (temp.size() < 3) return false;

    int n = temp.size();
    for (int i = 0; i < n; ++i) {
        if (cross_product(temp[i], temp[(i + 1) % n], temp[(i + 2) % n]) < -EPS) return false;
    }
    merged = temp;
    return true;
}

vector<Polygon> convex_decomposition(Polygon p) {
    vector<Polygon> polys = triangulate(p);
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < polys.size() && !merged; ++i) {
            for (size_t j = i + 1; j < polys.size() && !merged; ++j) {
                Polygon m;
                if (try_merge(polys[i], polys[j], m)) {
                    polys.erase(polys.begin() + j);
                    polys.erase(polys.begin() + i);
                    polys.push_back(m);
                    merged = true;
                }
            }
        }
    }
    return polys;
}

bool compare_pts(const Point& a, const Point& b) {
    if (abs(a.x - b.x) > EPS) return a.x < b.x;
    return a.y < b.y;
}

Polygon convex_hull(vector<Point>& pts) {
    int n = pts.size(), k = 0;
    if (n <= 2) return pts;
    vector<Point> h(2 * n);
    sort(pts.begin(), pts.end(), compare_pts);
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cross_product(h[k - 2], h[k - 1], pts[i]) <= 0) k--;
        h[k++] = pts[i];
    }
    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross_product(h[k - 2], h[k - 1], pts[i]) <= 0) k--;
        h[k++] = pts[i];
    }
    h.resize(k - 1);
    return h;
}

Polygon minkowski_diff(const Polygon& A, const Polygon& B) {
    vector<Point> pts;
    for (const auto& a : A) {
        for (const auto& b : B) {
            pts.push_back({a.x - b.x, a.y - b.y});
        }
    }
    return convex_hull(pts);
}

bool in_convex(const Point& pt, const Polygon& poly) {
    int n = poly.size();
    for (int i = 0; i < n; ++i) {
        if (cross_product(poly[i], poly[(i + 1) % n], pt) < -1e-6) return false;
    }
    return true;
}

int n1 = 0, n2 = 0, m = 0;
Polygon polyA, polyB;
vector<Polygon> NFP_list;

struct Cand {
    Point v;
    double len;
    bool operator<(const Cand& o) const { return len < o.len; }
};

void PreProcess() {
    vector<Polygon> decompA = convex_decomposition(polyA);
    vector<Polygon> decompB = convex_decomposition(polyB);
    for (const auto& a : decompA) {
        for (const auto& b : decompB) {
            NFP_list.push_back(minkowski_diff(a, b));
        }
    }
}

Point SolveQuery(const Point& P) {
    bool inside = false;
    for (const auto& M : NFP_list) {
        if (in_convex(P, M)) {
            inside = true; break;
        }
    }
    if (!inside) return {0.0, 0.0}; // 已经不重合

    vector<Cand> C;
    for (const auto& M : NFP_list) {
        int n = M.size();
        for (int i = 0; i < n; ++i) {
            Point a = M[i];
            Point b = M[(i + 1) % n];
            
            // 顶点到P的向量
            Point v_vert = a - P;
            C.push_back({v_vert, length(v_vert)});

            // P到边的垂线投影向量
            Point ab = b - a;
            Point ap = P - a;
            double t = dot(ap, ab) / dot(ab, ab);
            if (t >= 0 && t <= 1) {
                Point proj = a + (ab * t);
                Point v_edge = proj - P;
                C.push_back({v_edge, length(v_edge)});
            }
        }
    }
    
    sort(C.begin(), C.end());

    for (const auto& cand : C) {
        if (cand.len < 1e-9) continue;
        
        Point dir = {cand.v.x / cand.len, cand.v.y / cand.len};
        // 往外微移一点，避免卡在精度边界上
        Point test_pt = P + cand.v + (dir * 1e-4);

        bool still_inside = false;
        for (const auto& M : NFP_list) {
            if (in_convex(test_pt, M)) {
                still_inside = true; break;
            }
        }
        if (!still_inside) return cand.v;
    }
    
    return {0.0, 0.0};
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (!(cin >> n1 >> n2)) return 0;

    polyA.resize(n1);
    for (int i = 0; i < n1; ++i) cin >> polyA[i].x >> polyA[i].y;

    polyB.resize(n2);
    for (int i = 0; i < n2; ++i) cin >> polyB[i].x >> polyB[i].y;

    string okResp;
    cin >> okResp;
    if (okResp != "OK") return 0;

    PreProcess();

    cout << "OK\n";
    cout.flush();

    cin >> m;
    vector<Point> queries(m);
    for (int i = 0; i < m; ++i) cin >> queries[i].x >> queries[i].y;

    cin >> okResp;
    if (okResp != "OK") return 0;

    cout << m << "\n";
    for (int i = 0; i < m; ++i) {
        Point res = SolveQuery(queries[i]);
        cout << fixed << setprecision(5) << res.x << " " << res.y << "\n";
    }

    cout << "OK\n";
    cout.flush();
    return 0;
}