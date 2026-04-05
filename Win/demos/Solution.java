import java.io.*;
import java.util.*;

class Vector2D {
    double x, y;

    public Vector2D() {
        this.x = 0;
        this.y = 0;
    }

    public Vector2D(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public Vector2D subtract(Vector2D other) {
        return new Vector2D(this.x - other.x, this.y - other.y);
    }

    public Vector2D add(Vector2D other) {
        return new Vector2D(this.x + other.x, this.y + other.y);
    }

    public Vector2D multiply(double scalar) {
        return new Vector2D(this.x * scalar, this.y * scalar);
    }

    public double dot(Vector2D other) {
        return this.x * other.x + this.y * other.y;
    }

    public double length() {
        return Math.sqrt(this.x * this.x + this.y * this.y);
    }

    public Vector2D normalize() {
        double len = length();
        if (len == 0) {
            return this;
        }
        return new Vector2D(this.x / len, this.y / len);
    }

    public Vector2D perpendicular() {
        return new Vector2D(-this.y, this.x);
    }

    public Vector2D deepCopy() {
        return new Vector2D(this.x, this.y);
    }
}

class Polygon {
    List<Vector2D> vertices;

    public Polygon() {
        this.vertices = new ArrayList<>();
    }

    public Polygon(List<Vector2D> vertices) {
        this.vertices = vertices;
    }

    public Vector2D getCenter() {
        Vector2D center = new Vector2D(0, 0);
        if (vertices.size() == 0) {
            return center;
        }
        for (Vector2D v : vertices) {
            center = center.add(v);
        }
        return center.multiply(1.0 / vertices.size());
    }

    public void moveByVector(Vector2D vec) {
        for (int i = 0; i < vertices.size(); i++) {
            vertices.set(i, vertices.get(i).add(vec));
        }
    }

    public Polygon deepCopy() {
        List<Vector2D> newVertices = new ArrayList<>();
        for (Vector2D vertex : vertices) {
            newVertices.add(vertex.deepCopy());
        }
        return new Polygon(newVertices);
    }
}

class Projection {
    double min, max;

    public Projection (double min, double max) {
        this.min = min;
        this.max = max;
    }
}

public class Solution {
    private static final double EPS = 1e-6;

    private static int n1 = 0, n2 = 0, m = 0;
    private static Polygon polygon1 = new Polygon();
    private static Polygon polygon2 = new Polygon();
    private static List<Vector2D> testCases = new ArrayList<>();

    private static Projection projectPolygon(Polygon poly, Vector2D axis) {
        if (poly.vertices.size() == 0) {
            System.err.println("Polygon is empty and return Projection(0,0)");
            return new Projection(0.0, 0.0);
        }
        double minProj = poly.vertices.get(0).dot(axis);
        double maxProj = minProj;

        for (int i = 1; i < poly.vertices.size(); i++) {
            double proj = poly.vertices.get(i).dot(axis);
            if (proj < minProj) {
                minProj = proj;
            }
            if (proj > maxProj) {
                maxProj = proj;
            }
        }
        return new Projection(minProj, maxProj);
    }

    private static Vector2D genSolution(Vector2D vec) {
        Polygon polyA = new Polygon(polygon1.vertices);
        Polygon polyB = polygon2.deepCopy();
        polyB.moveByVector(vec);

        double minOverlap = Double.POSITIVE_INFINITY;
        Vector2D smallestAxis = new Vector2D();

        Polygon[] polygons = {polyA, polyB};

        for (int i = 0; i < 2; i++) {
            Polygon currentPoly = polygons[i];
            for (int j = 0; j < currentPoly.vertices.size(); j++) {
                Vector2D p1 = currentPoly.vertices.get(j);
                Vector2D p2 = currentPoly.vertices.get((j + 1) % currentPoly.vertices.size());
                Vector2D edge = p2.subtract(p1);

                Vector2D axis = edge.perpendicular().normalize();

                Projection projA = projectPolygon(polyA, axis);
                Projection projB = projectPolygon(polyB, axis);

                if (projA.max <= projB.min || projB.max <= projA.min) {
                    return new Vector2D(0.0, 0.0);
                }

                double overlap = Math.min(projA.max, projB.max) - Math.max(projA.min, projB.min);

                if (overlap < minOverlap) {
                    minOverlap = overlap;
                    smallestAxis = axis;
                }
            }
        }

        Vector2D centerA = polyA.getCenter();
        Vector2D centerB = polyB.getCenter();
        Vector2D dir = centerB.subtract(centerA);

        if (smallestAxis.dot(dir) < 0.0) {
            smallestAxis = smallestAxis.multiply(-1.0);
        }

        return smallestAxis.multiply(minOverlap);
    }

    private static void preProcess() {
        // player can perform preprocessing here
    }

    public static void main(String[] args) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        PrintWriter writer = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out)));

        // =============== 1. read polygons ===================
        StringTokenizer st = new StringTokenizer(reader.readLine());
        n1 = Integer.parseInt(st.nextToken());
        n2 = Integer.parseInt(st.nextToken());

        for (int i = 0; i < n1; ++i) {
            st = new StringTokenizer(reader.readLine());
            double x = Double.parseDouble(st.nextToken());
            double y = Double.parseDouble(st.nextToken());

            polygon1.vertices.add(new Vector2D(x, y));
        }

        for (int i = 0; i < n2; ++i) {
            st = new StringTokenizer(reader.readLine());
            double x = Double.parseDouble(st.nextToken());
            double y = Double.parseDouble(st.nextToken());
            polygon2.vertices.add(new Vector2D(x, y));
        }

        // wait for OK to ensure all polygon data is received
        String okResp = reader.readLine().trim();
        if (okResp == null || !okResp.equalsIgnoreCase("OK")) {
            System.err.println("Input data error: waiting for OK after obtaining polygons but I get " + okResp);
            return;
        }

        // ============== 2. preprocess ===================
        preProcess();
        // send OK after finishing preprocessing
        writer.println("OK");
        writer.flush();

        // ============== 3. read test data ===================
        m = Integer.parseInt(reader.readLine().trim());

        for (int i = 0; i < m; ++i) {
            st = new StringTokenizer(reader.readLine());
            double x = Double.parseDouble(st.nextToken());
            double y = Double.parseDouble(st.nextToken());

            testCases.add(new Vector2D(x, y));
        }

        // wait for OK to ensure all test cases are received
        okResp = reader.readLine().trim();
        if (okResp == null || !okResp.equalsIgnoreCase("OK")) {
            System.err.println("Input data error: waiting for OK after that I have received all test points but I get  " + okResp);
            return;
        }

        // ================ 4. solve and output results ===================
        for (int i = 0; i < m; ++i) {
            Vector2D res = genSolution(testCases.get(i));
            writer.printf("%.5f %.5f%n", res.x, res.y);
            writer.flush();
        }

        // send OK after outputting all answer
        writer.println("OK");
        writer.flush();

        writer.close();
        reader.close();
    }
}