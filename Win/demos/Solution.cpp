#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

const double EPS = 1e-7;

struct Point {
    double x, y;
    Point operator-(const Point& o) const { return {x - o.x, y - o.y}; }
    Point operator+(const Point& o) const { return {x + o.x, y + o.y}; }
    Point operator*(double s) const { return {x * s, y * s}; }
    double dot(const Point& o) const { return x * o.x + y * o.y; }
    double length() const { return std::hypot(x, y); }
    Point normalize() const {
        double l = length();
        return l < EPS ? *this : Point{x / l, y / l};
    }
    Point perp() const { return {-y, x}; }
};

double cross_product(const Point& o, const Point& a, const Point& b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// 求凸包，去重去共线，极大减少多边形边数
vector<Point> convex_hull(vector<Point>& pts) {
    int n = pts.size(), k = 0;
    if (n <= 2) return pts;
    vector<Point> h(2 * n);
    sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) {
        if (abs(a.x - b.x) > EPS) return a.x < b.x;
        return a.y < b.y;
    });
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cross_product(h[k - 2], h[k - 1], pts[i]) <= EPS) k--;
        h[k++] = pts[i];
    }
    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross_product(h[k - 2], h[k - 1], pts[i]) <= EPS) k--;
        h[k++] = pts[i];
    }
    h.resize(k - 1);
    return h;
}

// 预处理保存的投影轴信息，用于 O(1) 极速查询
struct AxisInfo {
    Point axis;
    double a_min, a_max;
    double b_base_min, b_base_max;
};

int n1 = 0, n2 = 0, m = 0;
vector<Point> polyA_raw, polyB_raw;
vector<AxisInfo> optimized_axes;
Point centerA, centerB_base;

void PreProcess() {
    vector<Point> hullA = convex_hull(polyA_raw);
    vector<Point> hullB = convex_hull(polyB_raw);

    auto get_center = [](const vector<Point>& p) {
        Point c{0, 0};
        for (const auto& v : p) { c.x += v.x; c.y += v.y; }
        return Point{c.x / p.size(), c.y / p.size()};
    };
    
    centerA = get_center(hullA);
    centerB_base = get_center(hullB);

    vector<Point> raw_axes;
    auto extract_axes = [&](const vector<Point>& hull) {
        int sz = hull.size();
        for (int i = 0; i < sz; ++i) {
            Point p1 = hull[i];
            Point p2 = hull[(i + 1) % sz];
            Point edge = p2 - p1;
            if (edge.length() > EPS) {
                raw_axes.push_back(edge.perp().normalize());
            }
        }
    };
    
    extract_axes(hullA);
    extract_axes(hullB);

    // 剔除重复/平行的轴，进一步压缩查询次数
    for (const auto& ax : raw_axes) {
        bool is_duplicate = false;
        for (const auto& ext : optimized_axes) {
            if (abs(ax.dot(ext.axis)) > 1.0 - EPS) { 
                is_duplicate = true; break; 
            }
        }
        if (!is_duplicate) {
            AxisInfo info;
            info.axis = ax;
            info.a_min = 1e15; info.a_max = -1e15;
            info.b_base_min = 1e15; info.b_base_max = -1e15;
            
            for (const auto& v : hullA) {
                double p = v.dot(ax);
                if (p < info.a_min) info.a_min = p;
                if (p > info.a_max) info.a_max = p;
            }
            for (const auto& v : hullB) {
                double p = v.dot(ax);
                if (p < info.b_base_min) info.b_base_min = p;
                if (p > info.b_base_max) info.b_base_max = p;
            }
            optimized_axes.push_back(info);
        }
    }
}

// 极速查询：没有任何 O(V) 的顶点遍历，全是 O(1) 的标量运算
Point SolveQuery(const Point& vec) {
    double min_overlap = 1e15;
    Point best_axis = {0, 0};

    for (const auto& ax : optimized_axes) {
        // 利用线性偏移常数特性，O(1) 更新多边形 B 的投影
        double shift = vec.dot(ax.axis);
        double b_min = ax.b_base_min + shift;
        double b_max = ax.b_base_max + shift;

        // 若在任意一个轴上分离，直接返回无重叠
        if (ax.a_max <= b_min + EPS || b_max <= ax.a_min + EPS) {
            return {0.0, 0.0};
        }

        double overlap = min(ax.a_max, b_max) - max(ax.a_min, b_min);
        if (overlap < min_overlap) {
            min_overlap = overlap;
            best_axis = ax.axis;
        }
    }

    // 确定推动方向：使 B 远离 A
    Point current_centerB = centerB_base + vec;
    Point dir = current_centerB - centerA;
    if (best_axis.dot(dir) < 0) {
        best_axis = best_axis * -1.0;
    }

    return best_axis * min_overlap;
}

int main() {
    // 开启极速 IO，防止 Docker 内产生不必要的耗时
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (!(cin >> n1 >> n2)) return 0;

    polyA_raw.resize(n1);
    for (int i = 0; i < n1; ++i) cin >> polyA_raw[i].x >> polyA_raw[i].y;

    polyB_raw.resize(n2);
    for (int i = 0; i < n2; ++i) cin >> polyB_raw[i].x >> polyB_raw[i].y;

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