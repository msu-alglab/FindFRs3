

package alglab.findfrs3;

public class FR implements Comparable<FR> {

    int node, support;

    FR(int node, int support) {
        this.node = node;
        this.support = support;
    }

    public int compareTo(FR x) {
        return Integer.compare(x.support, support);
    }
    
}
