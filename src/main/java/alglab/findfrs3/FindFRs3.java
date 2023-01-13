package alglab.findfrs3;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;
import org.apache.commons.io.*;
import org.apache.commons.cli.*;

public class FindFRs3 {

    static int k = 31;
    static double alpha = 0.5;
    static int minSup = 20;
    static String segFile = "";
    static String seqFile = "";
    static String bedFile = "";

    static HashMap<Integer, Integer> nodeIDtoIndex;
    static int numNodes = 0;
    static int[] nodeLength;

    static int numPaths = 0;
    static String[] seqName;
    static int[][] seqPath;
    static int[][] seqPathStarts;

    static ConcurrentHashMap<Integer, ConcurrentHashMap<Integer, AtomicInteger>> pairSup;
    static ConcurrentHashMap<Integer, AtomicInteger> sup;
    static TreeSet<Integer> frSet;
    static ConcurrentHashMap<Integer, ConcurrentLinkedQueue<int[]>> frSupPaths;
    static Variant[][] frVarSup;

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
        nodeIDtoIndex = new HashMap<>();
        HashMap<Integer, Integer> nodeIDtoLength = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> indextoNodeID = new HashMap<Integer, Integer>();
        try {
            LineIterator it = FileUtils.lineIterator(new File(segFile), "UTF-8");
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
            System.out.println(ex);
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
            LineIterator it = FileUtils.lineIterator(new File(seqFile), "UTF-8");
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

    static void processPathPairs(int s) {
        int[] path = seqPath[s];
        int a = -1, b = -1, aStart = 0, bStart = 0, nextSet = -1;
        TreeSet<Integer> nodesSeen = new TreeSet<>();
        for (int i = 0; i < path.length; i++) {
            if (i < path.length) {
                nextSet = find(path[i]);
            }
            if (a == -1) {
                a = nextSet;
                aStart = i;
            } else if (b == -1) {
                b = nextSet;
                bStart = i;
            }
            if ((nextSet != a && nextSet != b) || i == path.length) {
                nodesSeen.clear();
                for (int j = aStart; j < i; j++) {
                    nodesSeen.add(path[j]);
                }
                if (nodesSeen.size() > alpha * (size[a] + size[b])) {
                    pairSup.putIfAbsent(Math.min(a, b), new ConcurrentHashMap<>());
                    pairSup.get(Math.min(a, b)).putIfAbsent(Math.max(a, b), new AtomicInteger(0));
                    pairSup.get(Math.min(a, b)).get(Math.max(a, b)).incrementAndGet();
                }
                a = b;
                aStart = bStart;
                b = nextSet;
                bStart = i;
            }
        }
    }

    static void clusterNodes() {
        pairSup = new ConcurrentHashMap<>();
        int numMerges = 0;
        do {
            pairSup.clear();
            IntStream.range(0, numPaths).parallel().forEach(s -> {
                processPathPairs(s);
            });
            ArrayList<Edge> merges = new ArrayList<>();
            for (Integer A : pairSup.keySet()) {
                for (Integer B : pairSup.get(A).keySet()) {
                    if (pairSup.get(A).get(B).get() >= minSup) {
                        int a = A;
                        int b = B;
                        int sup = pairSup.get(A).get(B).get();
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
                }
            }
            System.out.println("number of merges: " + numMerges);
        } while (numMerges > 0);
    }

    static void processPathSupport(int s) {
        int[] path = seqPath[s];
        int a = -1, aStart = 0, nextSet = -1;
        TreeSet<Integer> nodesSeen = new TreeSet<>();
        for (int i = 0; i <= path.length; i++) {
            if (i < path.length) {
                nextSet = find(path[i]);
            }
            if (a == -1) {
                a = nextSet;
                aStart = i;
            }
            if (nextSet != a || i == path.length) {
                nodesSeen.clear();
                for (int j = aStart; j < i; j++) {
                    nodesSeen.add(path[j]);
                }
                if (nodesSeen.size() > alpha * size[a]) {
                    sup.putIfAbsent(a, new AtomicInteger(0));
                    sup.get(a).incrementAndGet();
                }
                a = nextSet;
                aStart = i;
            }
        }
    }

    static void processPathFindFRSupport(int s) {
        int[] path = seqPath[s];
        int a = -1, aStart = 0, nextSet = -1;
        TreeSet<Integer> nodesSeen = new TreeSet<>();
        for (int i = 0; i < path.length; i++) {
            if (i < path.length) {
                nextSet = find(path[i]);
            }
            if (a == -1) {
                a = nextSet;
                aStart = i;
            }
            if ((nextSet != a || i == path.length)) {
                if (frSet.contains(a)) {
                    nodesSeen.clear();
                    for (int j = aStart; j < i; j++) {
                        nodesSeen.add(path[j]);
                    }
                    if (nodesSeen.size() > alpha * size[a]) {
                        frSupPaths.putIfAbsent(a, new ConcurrentLinkedQueue<>());
                        int[] subpath = new int[3];
                        subpath[0] = s;
                        subpath[1] = aStart;
                        subpath[2] = i;
                        frSupPaths.get(a).add(subpath);
                    }
                }
                a = nextSet;
                aStart = i;
            }
        }
    }

    static void findFRs() {
        sup = new ConcurrentHashMap<>();
        IntStream.range(0, numPaths).parallel().forEach(s -> {
            processPathSupport(s);
        });
        frSet = new TreeSet<>();
        for (Integer I : sup.keySet()) {
            if (sup.get(I).get() >= minSup) {
                frSet.add(I);
            }
        }
        sup.clear();
        frSupPaths = new ConcurrentHashMap<>();
        IntStream.range(0, numPaths).parallel().forEach(s -> {
            processPathFindFRSupport(s);
        });
    }

    static void findFRVariants() {
        TreeMap<String, ArrayList<int[]>> variants = new TreeMap<>();
        int frNum = 0;
        int[] temp = new int[3];
        frVarSup = new Variant[frSet.size()][];
        for (Integer I : frSet) {
            variants.clear();
            for (int[] subpath : frSupPaths.get(I)) {
                String s = "";
                for (int j = subpath[1]; j < subpath[2]; j++) {
                    s += " " + seqPath[subpath[0]][j];
                }
                variants.putIfAbsent(s, new ArrayList<>());
                variants.get(s).add(subpath);
            }
            frVarSup[frNum] = new Variant[variants.keySet().size()];
            int v = 0;

            for (String s : variants.keySet()) {
                ArrayList<int[]> A = variants.get(s);
                frVarSup[frNum][v] = new Variant(new int[A.size()][]);
                int j = 0;
                for (int[] subpath : A) {
                    frVarSup[frNum][v].subpaths[j] = subpath;
                    j++;
                }
                v++;
            }
            Arrays.sort(frVarSup[frNum]);
            frNum++;
        }
        frSet.clear();
        frSupPaths.clear();
    }

    static void outputBED() {
        seqPathStarts = new int[numPaths][];
        for (int s = 0; s < numPaths; s++) {
            seqPathStarts[s] = new int[seqPath[s].length];
            int start = 0;
            for (int i = 0; i < seqPath[s].length; i++) {
                seqPathStarts[s][i] = start;
                start += nodeLength[seqPath[s][i]] - (k - 1);
            }
        }

        System.out.println("writing bed file");
        try {
            BufferedWriter bedOut = new BufferedWriter(new FileWriter(bedFile));
            for (int fr = 0; fr < frVarSup.length; fr++) {
                for (int v = 0; v < frVarSup[fr].length; v++) {
                    for (int i = 0; i < frVarSup[fr][v].subpaths.length; i++) {
                        int path = frVarSup[fr][v].subpaths[i][0];
                        int start = seqPathStarts[path][frVarSup[fr][v].subpaths[i][1]];
                        int stop = seqPathStarts[path][frVarSup[fr][v].subpaths[i][2] - 1] + nodeLength[seqPath[path][frVarSup[fr][v].subpaths[i][2] - 1]];
                        bedOut.write(seqName[path] // chrom
                                + "\t" + start // chromStart (starts with 0)
                                + "\t" + stop // chromEnd
                                + "\tfr" + fr + ":" + v // fr variant
                                + "\n");
                    }
                }
            }
            bedOut.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    public static void main(String[] args) {
        Options options = new Options();

        Option segO = new Option("n", "seg", true, "GFA3 seg file");
        segO.setRequired(true);
        options.addOption(segO);

        Option seqO = new Option("p", "seq", true, "GFA3 seq file");
        seqO.setRequired(true);
        options.addOption(seqO);

        Option bedO = new Option("b", "bed", true, "BED output file");
        bedO.setRequired(true);
        options.addOption(bedO);

        Option aO = new Option("a", "alpha", true, "alpha parameter");
        aO.setRequired(true);
        options.addOption(aO);

        Option kO = new Option("k", "kmer", true, "k-mer size");
        kO.setRequired(true);
        options.addOption(kO);

        Option minSO = new Option("m", "minsup", true, "minsup parameter");
        minSO.setRequired(true);
        options.addOption(minSO);

        HelpFormatter formatter = new HelpFormatter();
        try {
            CommandLineParser parser = new DefaultParser();

            CommandLine cmd = parser.parse(options, args);

            segFile = cmd.getOptionValue("seg");
            seqFile = cmd.getOptionValue("seq");
            bedFile = cmd.getOptionValue("bed");
            alpha = Double.parseDouble(cmd.getOptionValue("alpha"));
            k = Integer.parseInt(cmd.getOptionValue("kmer"));
            minSup = Integer.parseInt(cmd.getOptionValue("minsup"));

            readSeg();
            System.out.println("numNodes: " + numNodes);
            readSeq();
            System.out.println("numPaths: " + numPaths);

            clusterNodes();
            findFRs();
            findFRVariants();
            outputBED();

        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("FindFRs", options);
            System.exit(0);
        }

    }
}
