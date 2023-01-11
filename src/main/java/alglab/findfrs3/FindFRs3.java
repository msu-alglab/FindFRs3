/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Project/Maven2/JavaApp/src/main/java/${packagePath}/${mainClassName}.java to edit this template
 */
package alglab.findfrs3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import org.apache.commons.io.*;

public class FindFRs3 {

    static double alpha = 0.5;
    static int minSup = 20;

    static HashMap<Integer, Integer> nodeIDtoIndex;
    static int numNodes = 0;
    static int[] nodeLength;

    static int numPaths = 0;
    static String[] seqName;
    static int[][] seqPath;

    static TreeMap<Integer, TreeMap<Integer, Integer>> pairSup;

    // for union-find:
    static int[] parent;
    static int[] size;

    static int find(int i) {
        if (parent[i] != i) {
            parent[i] = find(parent[i]);
        }
        return parent[i];
    }

    static void union(int i, int j) {
        if (find(i) != find(j)) {
            size[find(i)] += size[find(j)];
            parent[find(j)] = find(i);
        }
    }

    static void readSeg() {
        HashMap<Integer, Integer> nodeIDtoLength = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> indextoNodeID = new HashMap<Integer, Integer>();
        try {
            LineIterator it = FileUtils.lineIterator(new File("yeast.k31.gfa3.cf_seg"), "UTF-8");
            try {
                while (it.hasNext()) {
                    String line = it.nextLine();
                    String[] lineSplit = line.split("\t");
                    int nodeID = Integer.parseInt(lineSplit[0]);
                    int length = lineSplit[1].length();
                    indextoNodeID.put(numNodes, nodeID);
                    nodeIDtoIndex.put(nodeID, numNodes);
                    numNodes++;
                    nodeIDtoLength.put(nodeID, length);
                }
            } finally {
                LineIterator.closeQuietly(it);
            }
        } catch (IOException ex) {
        }
        nodeLength = new int[numNodes];
        parent = new int[numNodes];
        size = new int[numNodes];

        for (int i = 0; i < numNodes; i++) {
            nodeLength[i] = nodeIDtoLength.get(indextoNodeID.get(i));
            parent[i] = i;
            size[i] = 1;
        }
    }

    static void readSeq() {
        HashMap<Integer, String> indexSeq = new HashMap<Integer, String>();
        HashMap<Integer, int[]> indexPath = new HashMap<Integer, int[]>();
        try {
            LineIterator it = FileUtils.lineIterator(new File("yeast.k31.gfa3.cf_seq"), "UTF-8");
            try {
                while (it.hasNext()) {
                    String line = it.nextLine();
                    String[] lineSplit = line.split("\t");
                    String seq = lineSplit[0].split(":")[2];
                    String[] pathStr = lineSplit[1].split(" ");
                    int[] path = new int[pathStr.length];
                    for (int i = 0; i < pathStr.length; i++) {
                        int node = Integer.parseInt(pathStr[i].substring(0, pathStr[i].length() - 1));
                        path[i] = nodeIDtoIndex.get(node);
                    }
                    indexSeq.put(numPaths, seq);
                    indexPath.put(numPaths, path);
                    numPaths++;
                }
            } finally {
                LineIterator.closeQuietly(it);
            }
        } catch (IOException ex) {
        }
        nodeIDtoIndex.clear();
        seqName = new String[numPaths];
        seqPath = new int[numPaths][];
        for (int i = 0; i < numPaths; i++) {
            seqName[i] = indexSeq.get(i);
            seqPath[i] = indexPath.get(i);
        }
    }

    static void processPath(int s) {
        int[] path = seqPath[s];
        int[] setPath = new int[path.length];
        for (int i = 0; i < path.length; i++) {
            setPath[i] = find(path[i]);
        }
        int a = -1, b = -1;
        int aStart = 0, bStart = 0;
        TreeSet<Integer> nodesSeen = new TreeSet<>();
        for (int i = 0; i < setPath.length; i++) {
            int nextSet = setPath[i];
            if (a == -1) {
                a = nextSet;
                aStart = i;
            } else if (b == -1) {
                b = nextSet;
                bStart = i;
            } else if (nextSet != a && nextSet != b) {
                nodesSeen.clear();
                for (int j = aStart; j < i; j++) {
                    nodesSeen.add(path[j]);
                }
                if (nodesSeen.size() > alpha * (size[a] + size[b])) {
                    //System.out.println(Math.min(a, b) + "-" + Math.max(a, b) + " supported by path " + s + "[" + aStart + "," + i + ")");
                    if (!pairSup.containsKey(Math.min(a, b))) {
                        pairSup.put(Math.min(a, b), new TreeMap<>());
                    }
                    if (!pairSup.get(Math.min(a, b)).containsKey(Math.max(a, b))) {
                        pairSup.get(Math.min(a, b)).put(Math.max(a, b), 0);
                    }
                    pairSup.get(Math.min(a, b)).put(Math.max(a, b), pairSup.get(Math.min(a, b)).get(Math.max(a, b)) + 1);
                }
                a = b;
                aStart = bStart;
                b = nextSet;
                bStart = i;
            }
        }
        if (b != -1) {
            nodesSeen.clear();
            for (int j = aStart; j < setPath.length; j++) {
                nodesSeen.add(path[j]);
            }
            if (nodesSeen.size() > alpha * (size[a] + size[b])) {
                //System.out.println(Math.min(a, b) + "-" + Math.max(a, b) + " supported by path " + s + "[" + aStart + "," + setPath.length + ")");
                if (!pairSup.containsKey(Math.min(a, b))) {
                    pairSup.put(Math.min(a, b), new TreeMap<>());
                }
                if (!pairSup.get(Math.min(a, b)).containsKey(Math.max(a, b))) {
                    pairSup.get(Math.min(a, b)).put(Math.max(a, b), 0);
                }
                pairSup.get(Math.min(a, b)).put(Math.max(a, b), pairSup.get(Math.min(a, b)).get(Math.max(a, b)) + 1);
            }
        }
    }

    public static void main(String[] args) {

        nodeIDtoIndex = new HashMap<>();
        readSeg();
        System.out.println("numNodes: " + numNodes);

        readSeq();
        System.out.println("numPaths: " + numPaths);

        pairSup = new TreeMap<>();

        int numMerges = 0;
        do {
            pairSup.clear();
            for (int s = 0; s < numPaths; s++) {
                if (s % 1000 == 0) {
                    System.out.println("processing path " + s);
                }
                processPath(s);
            }
            ArrayList<Edge> merges = new ArrayList<>();
            for (Integer A : pairSup.keySet()) {
                for (Integer B : pairSup.get(A).keySet()) {
                    if (pairSup.get(A).get(B) >= minSup) {
                        int a = A;
                        int b = B;
                        int sup = pairSup.get(A).get(B);
                        Edge e = new Edge(a, b, sup);
                        merges.add(e);
                    }
                }
            }
            Edge[] temp = new Edge[0];
            Edge[] m = merges.toArray(temp);
            merges.clear();
            Arrays.sort(m);
            numMerges = 0;
            boolean[] marked = new boolean[numNodes];
            for (int i = 0; i < m.length; i++) {
                if (!marked[m[i].a] && !marked[m[i].b]) {
                    marked[m[i].a] = marked[m[i].b] = true;
                    union(m[i].a, m[i].b);
                    numMerges++;
                    if (i < 100) {
                        System.out.println(m[i].a + "-" + m[i].b + " : " + m[i].support);
                    }
                }
            }
            System.out.println("number of merges: " + numMerges);

        } while (numMerges > 0);
    }
}
