# FindFRs3
new FindFRs ideas from 2023
This version of the FindFRs improves memory efficiency versus the previous version.  It now expects
the compressed DeBruijn graph to be built using Cuttlefish (which generates a .seq file and .seg file,
the nodes and paths through the graph).  The other parameters, alpha, kappa and minSup are the same as
the previous version.  A new optional parameter, minLen, specifies the minimum average supporting path
length for an FR to be reported (so short FRs can be filtered out).
