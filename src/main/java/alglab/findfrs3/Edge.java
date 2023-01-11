package alglab.findfrs3;

public class Edge implements Comparable<Edge> {

    int a, b, support;

    Edge(int a, int b, int support) {
        this.a = a;
        this.b = b;
        this.support = support;
    }

    public int compareTo(Edge x) {
        return Integer.compare(x.support, support);
    }
}
