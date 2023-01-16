
package alglab.findfrs3;

public class Variant implements Comparable<Variant> {

    int[][] subpaths;

    Variant(int[][] subpath) {
        this.subpaths = subpath;
    }

    public int compareTo(Variant x) {
        return Integer.compare(x.subpaths.length, subpaths.length);
    }
}