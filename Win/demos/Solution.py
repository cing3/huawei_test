import sys
import math
from typing import List, Tuple


# 常量定义
EPS = 1e-6


class Vector2D:
    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = x
        self.y = y

    def __sub__(self, other):
        return Vector2D(self.x - other.x, self.y - other.y)

    def __add__(self, other):
        return Vector2D(self.x + other.x, self.y + other.y)

    def __mul__(self, scalar: float):
        return Vector2D(self.x * scalar, self.y * scalar)

    def dot(self, other) -> float:
        return self.x * other.x + self.y * other.y

    def length(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y)

    def normalize(self):
        length = self.length()
        if length == 0:
            return self
        return Vector2D(self.x / length, self.y / length)

    def perp(self):
        return Vector2D(-self.y, self.x)


class Polygon:
    def __init__(self, vertices: List[Vector2D] = None):
        self.vertices = vertices if vertices is not None else []

    def get_center(self) -> Vector2D:
        center = Vector2D(0, 0)
        if not self.vertices:
            return center
        for v in self.vertices:
            center = center + v
        return center * (1.0 / len(self.vertices))

    def move_by_vec(self, vec: Vector2D):
        """将多边形平移指定向量"""
        for i in range(len(self.vertices)):
            self.vertices[i] = self.vertices[i] + vec


# 全局变量
n1 = 0
n2 = 0
m = 0
polygon1 = Polygon()
polygon2 = Polygon()
test_cases = []

class Projection:
    def __init__(self, min_val: float, max_val: float):
        self.min = min_val
        self.max = max_val


def project_polygon(poly: Polygon, axis: Vector2D) -> Projection:
    min_proj = poly.vertices[0].dot(axis)
    max_proj = min_proj

    for i in range(1, len(poly.vertices)):
        proj = poly.vertices[i].dot(axis)
        if proj < min_proj:
            min_proj = proj
        if proj > max_proj:
            max_proj = proj

    return Projection(min_proj, max_proj)


def gen_solution(vec: Vector2D) -> Vector2D:
    poly_b = Polygon(polygon2.vertices[:])
    poly_b.move_by_vec(vec)

    min_overlap = float('inf')
    smallest_axis = Vector2D(0, 0)

    polygons = [polygon1, poly_b]

    for i in range(2):
        current_poly = polygons[i]
        for j in range(len(current_poly.vertices)):
            p1 = current_poly.vertices[j]
            p2 = current_poly.vertices[(j + 1) % len(current_poly.vertices)]
            edge = p2 - p1

            axis = edge.perp().normalize()

            proj_a = project_polygon(polygon1, axis)
            proj_b = project_polygon(poly_b, axis)

            if proj_a.max <= proj_b.min or proj_b.max <= proj_a.min:
                return Vector2D(0.0, 0.0)

            overlap = min(proj_a.max, proj_b.max) - max(proj_a.min, proj_b.min)

            if overlap < min_overlap:
                min_overlap = overlap
                smallest_axis = axis

    center_a = polygon1.get_center()
    center_b = poly_b.get_center()
    dir_vec = center_b - center_a

    if smallest_axis.dot(dir_vec) < 0.0:
        smallest_axis = smallest_axis * -1.0

    return smallest_axis * min_overlap


def pre_process():
    # player can perform preprocessing here
    pass


def main():
    global n1, n2, m, polygon1, polygon2, test_cases

    # =============== 1. read polygons ===================
    line = sys.stdin.readline()
    try:
        n1, n2 = map(int, line.split())
    except ValueError:
        print('Input data error: wrong polygon vertex count.', file=sys.stderr)


    for _ in range(n1):
        x, y = map(float, sys.stdin.readline().split())
        polygon1.vertices.append(Vector2D(x, y))

    for _ in range(n2):
        x, y = map(float, sys.stdin.readline().split())
        polygon2.vertices.append(Vector2D(x, y))

    # wait for OK to ensure all polygon data is received
    ok_response = sys.stdin.readline().strip()
    if ok_response != "OK":
        print(f"Input data error: waiting for OK after obtaining polygons but I get {ok_response}", file=sys.stderr)
        return

    # ============== 2. preprocess ===================
    pre_process()
    # send OK after finishing preprocessing
    print("OK")
    sys.stdout.flush()

    # ============== 3. read test points ===================
    m = int(sys.stdin.readline())

    for _ in range(m):
        x, y = map(float, sys.stdin.readline().split())
        test_cases.append(Vector2D(x, y))

    # wait for OK to ensure all test cases are received
    ok_response = sys.stdin.readline().strip()
    if ok_response != "OK":
        print(f'Input data error: waiting for OK after that I have received all test points but I get {ok_response}', file=sys.stderr)
        return

    # ================ 4. solve and output results ===================
    print(m)
    sys.stdout.flush()
    for testCase in test_cases:
        res = gen_solution(testCase)
        print(f"{res.x:.5f} {res.y:.5f}")
        sys.stdout.flush()

    # send OK after outputting all answer
    print("OK")
    sys.stdout.flush()

if __name__ == "__main__":
    main()
