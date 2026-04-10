#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

struct Vector2D;
struct Polygon;

using namespace std;
const double EPS = 1e-6;

struct Vector2D {
    double x, y;

    Vector2D(double x = 0, double y = 0) : x(x), y(y) {}
    Vector2D operator-(const Vector2D& other) const { return Vector2D(x - other.x, y - other.y); }
    Vector2D operator+(const Vector2D& other) const { return Vector2D(x + other.x, y + other.y); }
    Vector2D operator*(double scalar) const { return Vector2D(x * scalar, y * scalar); }
    double Dot(const Vector2D& other) const { return x * other.x + y * other.y; }
    double Length() const { return std::sqrt(x * x + y * y); }
    Vector2D Normalize() const
    {
        double len = Length();
        if (len == 0) {
            return *this;
        }
        return Vector2D(x / len, y / len);
    }
    Vector2D Perp() const { return Vector2D(-y, x); }
};

struct Polygon {
    std::vector<Vector2D> vertices;

    Polygon() = default;
    Polygon(std::initializer_list<Vector2D> vts) : vertices(vts) {}
    Vector2D GetCenter() const
    {
        Vector2D center(0, 0);
        if (vertices.empty()) {
            return center;
        }
        for (const auto& v : vertices) {
            center = center + v;
        }
        return center * (1.0 / vertices.size());
    }
    void MoveByVec(const Vector2D& vec) {
        for (auto& v : vertices) {
            v = v + vec;
        }
    }
};

int n1 = 0, n2 = 0, m = 0;
Polygon polygon1;
Polygon polygon2;
vector<Vector2D> testCases;

struct Projection {
    double min, max;
};

Projection ProjectPolygon(const Polygon& poly, const Vector2D& axis)
{
    double minProj = poly.vertices[0].Dot(axis);
    double maxProj = minProj;

    for (size_t i = 1; i < poly.vertices.size(); ++i) {
        double proj = poly.vertices[i].Dot(axis);
        if (proj < minProj) {
            minProj = proj;
        }
        if (proj > maxProj) {
            maxProj = proj;
        }
    }
    return {minProj, maxProj};
}

Vector2D GenSolution(const Vector2D& vec)
{
    Polygon polyB = polygon2;
    polyB.MoveByVec(vec);

    double minOverlap = std::numeric_limits<double>::infinity();
    Vector2D smallestAxis;

    const Polygon* polygons[2] = {&polygon1, &polyB};

    for (int i = 0; i < 2; ++i) {
        const Polygon& currentPoly = *polygons[i];
        for (size_t j = 0; j < currentPoly.vertices.size(); ++j) {
            Vector2D p1 = currentPoly.vertices[j];
            Vector2D p2 = currentPoly.vertices[(j + 1) % currentPoly.vertices.size()];
            Vector2D edge = p2 - p1;

            Vector2D axis = edge.Perp().Normalize();

            Projection projA = ProjectPolygon(polygon1, axis);
            Projection projB = ProjectPolygon(polyB, axis);

            if (projA.max <= projB.min || projB.max <= projA.min) {
                return {0.0, 0.0};
            }

            double overlap = std::min(projA.max, projB.max) - std::max(projA.min, projB.min);

            if (overlap < minOverlap) {
                minOverlap = overlap;
                smallestAxis = axis;
            }
        }
    }
    Vector2D centerA = polygon1.GetCenter();
    Vector2D centerB = polyB.GetCenter();
    Vector2D dir = centerB - centerA;

    if (smallestAxis.Dot(dir) < 0.0) {
        smallestAxis = smallestAxis * -1.0;
    }

    return {smallestAxis * minOverlap};
}

// 选手在规定的时间内进行预处理，完成后返回OK
void PreProcess()
{
    // player can perform preprocessing here
}

int main()
{
    // =============== 1. read polygons ===================
    cin >> n1 >> n2;

    if (n1 <= 0 || n2 <= 0) {
        cerr << "Input data error: the number of vertices of both polygons should be greater than 2" << endl;
        return 1;
    }

    polygon1.vertices.resize(n1);
    for (int i = 0; i < n1; ++i) {
        cin >> polygon1.vertices[i].x >> polygon1.vertices[i].y;
    }

    polygon2.vertices.resize(n2);
    for (int i = 0; i < n2; ++i) {
        cin >> polygon2.vertices[i].x >> polygon2.vertices[i].y;
    }

    // wait for OK to ensure all polygon data is received
    string okResp;
    cin >> okResp;
    if (okResp != "OK") {
        cerr << "Input data error: waiting for OK after obtaining polygons but I get " << okResp << endl;
        return 0;
    }

    // ============== 2. preprocess ===================
    PreProcess();
    // send OK after finishing preprocessing
    cout << "OK" << endl;
    cout.flush();

    // ============== 3. read test data ===================
    cin >> m;
    testCases.resize(m);

    for (int i = 0; i < m; ++i) {
        cin >> testCases[i].x >> testCases[i].y;
    }

    // wait for OK to ensure all test cases are received
    cin >> okResp;
    if (okResp != "OK") {
        cerr << "Input data error: waiting for OK after that I have received all test points but I get " << okResp
             << endl;
        return 0;
    }

    // ================ 4. solve and output results ===================
    cout << m << endl;
    for (int i = 0; i < m; ++i) {
        const Vector2D& res = GenSolution(testCases[i]);
        cout << fixed << std::setprecision(5) << res.x << " " << res.y << endl;
        cout.flush();
    }

    // send OK after outputting all answer
    cout << "OK" << endl;
    cout.flush();

    return 0;
}