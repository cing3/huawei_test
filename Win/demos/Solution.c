#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EPS 1e-6

typedef struct {
    double x, y;
} Vector2D;

typedef struct {
    Vector2D* vertices;
    int size;
} Polygon;

typedef struct {
    double min, max;
} Projection;

int n1 = 0, n2 = 0, m = 0;
Polygon polygon1;
Polygon polygon2;
Vector2D* testCases;

Vector2D vector_sub(Vector2D a, Vector2D b) { return (Vector2D){a.x - b.x, a.y - b.y}; }
Vector2D vector_add(Vector2D a, Vector2D b) { return (Vector2D){a.x + b.x, a.y + b.y}; }
Vector2D vector_mul(Vector2D a, double scalar) { return (Vector2D){a.x * scalar, a.y * scalar}; }
double vector_dot(Vector2D a, Vector2D b) { return a.x * b.x + a.y * b.y; }
double vector_length(Vector2D a) { return sqrt(a.x * a.x + a.y * a.y); }
Vector2D vector_normalize(Vector2D a)
{
    double len = vector_length(a);
    if (len == 0) {
        return a;
    }
    return (Vector2D){a.x / len, a.y / len};
}
Vector2D vector_perp(Vector2D a) { return (Vector2D){-a.y, a.x}; }


Vector2D polygon_get_center(Polygon poly)
{
    Vector2D center = {0.0, 0.0};
    for (int i = 0; i < poly.size; i++) {
        center = vector_add(center, poly.vertices[i]);
    }
    return vector_mul(center, 1.0 / poly.size);
}

void polygon_move_by_vec(Polygon* poly, Vector2D vec)
{
    for (int i = 0; i < poly->size; i++) {
        poly->vertices[i] = vector_add(poly->vertices[i], vec);
    }
}

Projection project_polygon(Polygon poly, Vector2D axis)
{
    double minProj = vector_dot(poly.vertices[0], axis);
    double maxProj = minProj;

    for (int i = 1; i < poly.size; i++) {
        double proj = vector_dot(poly.vertices[i], axis);
        if (proj < minProj) {
            minProj = proj;
        }
        if (proj > maxProj) {
            maxProj = proj;
        }
    }
    return (Projection){minProj, maxProj};
}

double get_overlap(Projection projA, Projection projB)
{
    double overlap = fmin(projA.max, projB.max) - fmax(projA.min, projB.min);
    return overlap > 0 ? overlap : 0;
}

int is_overlap(Projection projA, Projection projB) { return !(projA.max <= projB.min || projB.max <= projA.min); }

Vector2D gen_solution(Vector2D vec)
{
    Polygon polyB;
    polyB.size = polygon2.size;
    polyB.vertices = (Vector2D*)malloc(sizeof(Vector2D) * polyB.size);
    if (polyB.vertices == NULL) {
        fprintf(stderr, "Memory allocation failed for polygon.vertices\n");
        return (Vector2D){0.0, 0.0};
    }

    for (int i = 0; i < polyB.size; i++) {
        polyB.vertices[i] = polygon2.vertices[i];
    }

    polygon_move_by_vec(&polyB, vec);

    double minOverlap = INFINITY;
    Vector2D smallestAxis = {0.0, 0.0};

    Polygon* polygons[2] = {&polygon1, &polyB};

    for (int i = 0; i < 2; i++) {
        Polygon* currentPoly = polygons[i];
        for (int j = 0; j < currentPoly->size; j++) {
            Vector2D p1 = currentPoly->vertices[j];
            Vector2D p2 = currentPoly->vertices[(j + 1) % currentPoly->size];
            Vector2D edge = vector_sub(p2, p1);

            Vector2D axis = vector_perp(edge);
            axis = vector_normalize(axis);

            Projection projA = project_polygon(polygon1, axis);
            Projection projB = project_polygon(polyB, axis);

            if (!is_overlap(projA, projB)) {
                free(polyB.vertices);
                return (Vector2D){0.0, 0.0};
            }

            double overlap = get_overlap(projA, projB);

            if (overlap < minOverlap) {
                minOverlap = overlap;
                smallestAxis = axis;
            }
        }
    }

    Vector2D centerA = polygon_get_center(polygon1);
    Vector2D centerB = polygon_get_center(polyB);
    Vector2D dir = vector_sub(centerB, centerA);

    if (vector_dot(smallestAxis, dir) < 0.0) {
        smallestAxis = vector_mul(smallestAxis, -1.0);
    }

    free(polyB.vertices);

    return vector_mul(smallestAxis, minOverlap);
}

// 预处理函数
void pre_process()
{
    // player can perform preprocessing here
}

int main()
{
    // =============== 1. read polygons ===================
    int ret = scanf("%d %d", &n1, &n2);
    if (ret != 2) {
        fprintf(stderr, "Input data error: can not get number of vertices of polygons.\n");
    }

    polygon1.size = n1;
    polygon2.size = n2;

    polygon1.vertices = (Vector2D*)malloc(sizeof(Vector2D) * n1);
    if (polygon1.vertices == NULL) {
        fprintf(stderr, "Memory allocation failed for polygon1.vertices\n");
        return 0;
    }
    for (int i = 0; i < n1; i++) {
        ret = scanf("%lf %lf", &polygon1.vertices[i].x, &polygon1.vertices[i].y);
        if (ret != 2) {
            fprintf(stderr, "Input data error: can not get polygon 1 data.\n");
        }
    }

    polygon2.vertices = (Vector2D*)malloc(sizeof(Vector2D) * n2);
    if (polygon2.vertices == NULL) {
        fprintf(stderr, "Memory allocation failed for polygon2.vertices\n");
        free(polygon1.vertices);
        polygon1.vertices = NULL;
        return 0;
    }
    for (int i = 0; i < n2; i++) {
        ret = scanf("%lf %lf", &polygon2.vertices[i].x, &polygon2.vertices[i].y);
        if (ret != 2) {
            fprintf(stderr, "Input data error: can not get polygon 2 data.\n");
        }
    }

    char okResp[10];
    ret = scanf("%9s", okResp);
    if (ret != 1) {
        fprintf(stderr, "Input data error: can not get OK after get polygons.\n");
    }
    if (strcmp(okResp, "OK") != 0) {
        fprintf(stderr, "Input data error: waiting for OK after obtaining polygons but I get %s\n", okResp);
        return 0;
    }

    // ============== 2. preprocess ===================
    pre_process();
    printf("OK\n");
    fflush(stdout);

    // ============== 3. read test data ===================
    ret = scanf("%d", &m);
    if (ret != 1) {
        fprintf(stderr, "Input data error: can not get number of test cases.\n");
    }

    testCases = (Vector2D*)malloc(sizeof(Vector2D) * m);
    if (testCases == NULL) {
        fprintf(stderr, "Memory allocation failed for testCases\n");
        free(polygon1.vertices);
        free(polygon2.vertices);
        polygon1.vertices = NULL;
        polygon2.vertices = NULL;
        return 0;
    }

    for (int i = 0; i < m; i++) {
        ret = scanf("%lf %lf", &testCases[i].x, &testCases[i].y);
        if (ret != 2) {
            fprintf(stderr, "Input data error: can not get test case %d\n", i);
        }
    }

    // wait for OK to ensure all test cases are received
    ret = scanf("%9s", okResp);
    if (ret != 1) {
        fprintf(stderr, "Input data error: can not get OK after get all test cases.\n");
    }
    if (strcmp(okResp, "OK") != 0) {
        fprintf(stderr, "Input data error: waiting for OK after that I have received all test points but I get %s\n",
                okResp);
        return 0;
    }

    // ================ 4. solve and output results ===================
    for (int i = 0; i < m; i++) {
        Vector2D res = gen_solution(testCases[i]);
        printf("%.5f %.5f\n", res.x, res.y);
        fflush(stdout);
    }

    // send OK after outputting all answer
    printf("OK\n");
    fflush(stdout);

    free(polygon1.vertices);
    free(polygon2.vertices);
    free(testCases);

    return 0;
}