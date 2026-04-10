// 强制开启 O3 优化、循环展开以及 AVX2 向量化指令集 (打榜标配)
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio> // 使用更快的 C 风格 IO

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

struct AxisInfo {
    double ax, ay; // 拆解 Point 以增加缓存命中率
    double a_min, a_max;
    double b_base_min, b_base_max;
    double initial_overlap; 
};

int n1 = 0, n2 = 0, m = 0;
vector<Point> polyA_raw, polyB_raw;
vector<AxisInfo> optimized_axes;
Point centerA, centerB_base;

void PreProcess() {
    vector<Point> hullA = convex_hull(polyA_raw);
    vector<Point> hullB = convex_hull(polyB_raw);

    auto get_center = [](const vector<Point>& p) {
        double cx = 0, cy = 0;
        for (const auto& v : p) { cx += v.x; cy += v.y; }
        return Point{cx / p.size(), cy / p.size()};
    };
    
    centerA = get_center(hullA);
    centerB_base = get_center(hullB);

    vector<Point> raw_axes;
    
    // 强制加入 X 轴和 Y 轴。它们等价于 AABB 包围盒检测，能够最快拦截大部分无碰撞情况
    raw_axes.push_back({1.0, 0.0});
    raw_axes.push_back({0.0, 1.0});

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

    for (const auto& ax : raw_axes) {
        bool is_duplicate = false;
        for (const auto& ext : optimized_axes) {
            if (abs(ax.x * ext.ax + ax.y * ext.ay) > 1.0 - EPS) { 
                is_duplicate = true; break; 
            }
        }
        if (!is_duplicate) {
            AxisInfo info;
            info.ax = ax.x;
            info.ay = ax.y;
            info.a_min = 1e15; info.a_max = -1e15;
            info.b_base_min = 1e15; info.b_base_max = -1e15;
            
            for (const auto& v : hullA) {
                double p = v.x * ax.x + v.y * ax.y;
                if (p < info.a_min) info.a_min = p;
                if (p > info.a_max) info.a_max = p;
            }
            for (const auto& v : hullB) {
                double p = v.x * ax.x + v.y * ax.y;
                if (p < info.b_base_min) info.b_base_min = p;
                if (p > info.b_base_max) info.b_base_max = p;
            }
            
            // 记录初始重叠度，用作排序优先级
            info.initial_overlap = (info.a_max < info.b_base_max ? info.a_max : info.b_base_max) - 
                                   (info.a_min > info.b_base_min ? info.a_min : info.b_base_min);
            optimized_axes.push_back(info);
        }
    }

    // 核心优化：将更容易分离的轴放在前面，触发提前 return
    sort(optimized_axes.begin(), optimized_axes.end(), [](const AxisInfo& a, const AxisInfo& b) {
        return a.initial_overlap < b.initial_overlap;
    });
}

// 内联极速查询函数
inline Point SolveQuery(double vx, double vy) {
    double min_overlap = 1e15;
    double best_ax = 0, best_ay = 0;

    int num_axes = optimized_axes.size();
    for (int i = 0; i < num_axes; ++i) {
        const auto& ax = optimized_axes[i];
        
        double shift = vx * ax.ax + vy * ax.ay;
        double b_min = ax.b_base_min + shift;
        double b_max = ax.b_base_max + shift;

        // 快逃！不碰撞直接中断
        if (ax.a_max <= b_min + EPS || b_max <= ax.a_min + EPS) {
            return {0.0, 0.0};
        }

        // 用三元运算替代 std::min/max，利于硬件生成 CMOV 指令
        double max_of_mins = (ax.a_min > b_min) ? ax.a_min : b_min;
        double min_of_maxs = (ax.a_max < b_max) ? ax.a_max : b_max;
        double overlap = min_of_maxs - max_of_mins;

        if (overlap < min_overlap) {
            min_overlap = overlap;
            best_ax = ax.ax;
            best_ay = ax.ay;
        }
    }

    double dir_x = (centerB_base.x + vx) - centerA.x;
    double dir_y = (centerB_base.y + vy) - centerA.y;
    
    if (best_ax * dir_x + best_ay * dir_y < 0) {
        return {-best_ax * min_overlap, -best_ay * min_overlap};
    }

    return {best_ax * min_overlap, best_ay * min_overlap};
}

int main() {
    // 极速 IO 宏拦截
    if (scanf("%d %d", &n1, &n2) != 2) return 0;

    polyA_raw.resize(n1);
    for (int i = 0; i < n1; ++i) scanf("%lf %lf", &polyA_raw[i].x, &polyA_raw[i].y);

    polyB_raw.resize(n2);
    for (int i = 0; i < n2; ++i) scanf("%lf %lf", &polyB_raw[i].x, &polyB_raw[i].y);

    char okResp[16];
    scanf("%s", okResp);

    PreProcess();

    printf("OK\n");
    fflush(stdout);

    scanf("%d", &m);
    vector<Point> queries(m);
    for (int i = 0; i < m; ++i) {
        scanf("%lf %lf", &queries[i].x, &queries[i].y);
    }

    scanf("%s", okResp);

    printf("%d\n", m);
    for (int i = 0; i < m; ++i) {
        Point res = SolveQuery(queries[i].x, queries[i].y);
        // 使用纯 C 的 printf 保证严苛情况下的最高输出效率
        printf("%.5f %.5f\n", res.x, res.y);
    }

    printf("OK\n");
    fflush(stdout);
    return 0;
}