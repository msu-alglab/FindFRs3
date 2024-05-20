# FindFRs3

This version of the FindFRs improves memory efficiency versus the previous version.  It expects
the compressed DeBruijn graph to be built using [Cuttlefish](https://github.com/COMBINE-lab/cuttlefish)
(with the -f3 option that generates a .seq file and .seg file, the nodes and paths through the graph).  
The other parameters, alpha (minimum penetrance), kappa (maximum insert size) and minSup (minimum support) are the same as
the previous version.  A new optional parameter, minLen, specifies the minimum average supporting path
length for an FR to be reported (so short FRs can be filtered out).
