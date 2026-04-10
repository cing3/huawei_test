#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring>

using namespace std;

const double EPS = 1e-6;

struct Point {
    double x, y;
    Point() : x(0), y(0) {}
    Point(double x_, double y_) : x(x_), y(y_) {}
    Point operator-(const Point& o) const { return Point(x - o.x, y - o.y); }
    Point operator+(const Point& o) const { return Point(x + o.x, y + o.y); }
    Point operator*(double s) const { return Point(x * s, y * s); }
    double dot(const Point& o) const { return x * o.x + y * o.y; }
    double cross(const Point& o) const { return x * o.y - y * o.x; }
    double length() const { return hypot(x, y); }
    Point normalize() const {
        double l = length();
        return l < EPS ? Point(0, 0) : Point(x / l, y / l);
    }
    Point perp() const { return Point(-y, x); }
    bool operator==(const Point& o) const {
        return fabs(x - o.x) < EPS && fabs(y - o.y) < EPS;
    }
};

// 叉积：o->a × o->b
double cross_product(const Point& o, const Point& a, const Point& b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// 凸包优化：Andrew算法，严格去重+去共线
vector<Point> convex_hull(vector<Point>& pts) {
    int n = pts.size();
    if (n <= 1) return pts;
    if (n == 2) {
        if (pts[0] == pts[1]) return {pts[0]};
        return pts;
    }

    sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) {
        return fabs(a.x - b.x) < EPS ? a.y < b.y - EPS : a.x < b.x - EPS;
    });

    vector<Point> hull;
    hull.reserve(n);
    // 下凸壳
    for (int i = 0; i < n; ++i) {
        while (hull.size() >= 2) {
            Point& a = hull[hull.size() - 2];
            Point& b = hull[hull.size() - 1];
            double cr = cross_product(a, b, pts[i]);
            if (cr < -EPS) break;
            if (fabs(cr) < EPS) hull.pop_back();
            else hull.pop_back();
        }
        hull.push_back(pts[i]);
    }

    // 上凸壳
    int lower_size = hull.size();
    for (int i = n - 2; i >= 0; --i) {
        while (hull.size() > lower_size) {
            Point& a = hull[hull.size() - 2];
            Point& b = hull[hull.size() - 1];
            double cr = cross_product(a, b, pts[i]);
            if (cr < -EPS) break;
            if (fabs(cr) < EPS) hull.pop_back();
            else hull.pop_back();
        }
        hull.push_back(pts[i]);
    }

    hull.pop_back(); // 移除重复的起点
    return hull;
}

// 多边形面积加权质心（比简单平均更精准）
Point polygon_centroid(const vector<Point>& poly) {
    int n = poly.size();
    if (n == 0) return Point(0, 0);
    if (n == 1) return poly[0];

    double area = 0.0;
    Point centroid(0, 0);
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        double cr = poly[i].x * poly[j].y - poly[j].x * poly[i].y;
        area += cr;
        centroid.x += (poly[i].x + poly[j].x) * cr;
        centroid.y += (poly[i].y + poly[j].y) * cr;
    }
    area *= 0.5;
    if (fabs(area) < EPS) { // 退化多边形（线/点）
        centroid = Point(0, 0);
        for (const Point& p : poly) {
            centroid.x += p.x;
            centroid.y += p.y;
        }
        centroid.x /= n;
        centroid.y /= n;
        return centroid;
    }
    centroid.x /= (6 * area);
    centroid.y /= (6 * area);
    return centroid;
}

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

    // 计算质心
    centerA = polygon_centroid(hullA);
    centerB_base = polygon_centroid(hullB);

    // 提取所有候选轴（凸包边的垂直法向量）
    vector<Point> raw_axes;
    auto extract_axes = [&](const vector<Point>& hull) {
        int sz = hull.size();
        for (int i = 0; i < sz; ++i) {
            Point p1 = hull[i];
            Point p2 = hull[(i + 1) % sz];
            Point edge = p2 - p1;
            if (edge.length() < EPS) continue; // 过滤零长度边
            Point axis = edge.perp().normalize();
            raw_axes.push_back(axis);
        }
    };
    extract_axes(hullA);
    extract_axes(hullB);

    // 轴去重：合并正反方向轴，减少遍历次数
    for (const Point& ax : raw_axes) {
        bool duplicate = false;
        for (const AxisInfo& info : optimized_axes) {
            double dot_val = ax.dot(info.axis);
            if (fabs(dot_val - 1.0) < EPS || fabs(dot_val + 1.0) < EPS) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) continue;

        AxisInfo info;
        info.axis = ax;
        info.a_min = 1e18;
        info.a_max = -1e18;
        info.b_base_min = 1e18;
        info.b_base_max = -1e18;

        // 预处理多边形A的投影范围
        for (const Point& p : hullA) {
            double proj = p.dot(ax);
            if (proj < info.a_min) info.a_min = proj;
            if (proj > info.a_max) info.a_max = proj;
        }
        // 预处理多边形B的基础投影范围
        for (const Point& p : hullB) {
            double proj = p.dot(ax);
            if (proj < info.b_base_min) info.b_base_min = proj;
            if (proj > info.b_base_max) info.b_base_max = proj;
        }
        optimized_axes.push_back(info);
    }
}

Point SolveQuery(const Point& vec) {
    double min_overlap = 1e18;
    Point best_axis(0, 0);

    for (const AxisInfo& ax : optimized_axes) {
        // O(1)计算平移后B的投影（利用线性偏移特性）
        double shift = vec.dot(ax.axis);
        double b_min = ax.b_base_min + shift;
        double b_max = ax.b_base_max + shift;

        // 轴分离：无重叠直接返回0向量
        if (ax.a_max + EPS <= b_min || b_max + EPS <= ax.a_min) {
            return Point(0, 0);
        }

        // 计算重叠长度
        double overlap = min(ax.a_max, b_max) - max(ax.a_min, b_min);
        if (overlap < min_overlap - EPS) {
            min_overlap = overlap;
            best_axis = ax.axis;
        }
    }

    // 无有效轴（理论上不会走到这）
    if (fabs(min_overlap - 1e18) < EPS) {
        return Point(0, 0);
    }

    // 确定推动方向：让B远离A
    Point current_centerB = centerB_base + vec;
    Point dir = current_centerB - centerA;
    if (best_axis.dot(dir) < -EPS) {
        best_axis = best_axis * -1.0;
    }

    return best_axis * min_overlap;
}

int main() {
    // 1. 读取多边形顶点
    if (scanf("%d %d", &n1, &n2) != 2) return 0;

    polyA_raw.resize(n1);
    for (int i = 0; i < n1; ++i) {
        scanf("%lf %lf", &polyA_raw[i].x, &polyA_raw[i].y);
    }
    polyB_raw.resize(n2);
    for (int i = 0; i < n2; ++i) {
        scanf("%lf %lf", &polyB_raw[i].x, &polyB_raw[i].y);
    }

    char okResp[10];
    if (scanf("%9s", okResp) != 1 || strcmp(okResp, "OK") != 0) {
        return 0;
    }

    // 2. 预处理
    PreProcess();
    printf("OK\n");
    fflush(stdout);

    // 3. 读取测试用例
    if (scanf("%d", &m) != 1) return 0;
    vector<Point> queries(m);
    for (int i = 0; i < m; ++i) {
        scanf("%lf %lf", &queries[i].x, &queries[i].y);
    }
    if (scanf("%9s", okResp) != 1 || strcmp(okResp, "OK") != 0) {
        return 0;
    }

    // 4. 计算并输出结果
    printf("%d\n", m);
    for (int i = 0; i < m; ++i) {
        Point res = SolveQuery(queries[i]);
        printf("%.5f %.5f\n", res.x, res.y);
        fflush(stdout);
    }
    printf("OK\n");
    fflush(stdout);

    return 0;
}