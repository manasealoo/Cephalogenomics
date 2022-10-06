# Asymmetron genome assembly

## Pre-assembly analysis

```
../jellyfish/bin/jellyfish count -C -m 21 -s 1000M -t 10 ASY_*.fastq -o ASY_both.jf
../jellyfish/bin/jellyfish histo -t 10 ASY_both.jf > ASY_both.histo
```

The resulting distribution was inputted into GenomeScope with default options of 

```
p = 2
k = 21
```
and
```
Max k-mer coverage = -1
Average k-mer coverage for polyploid genome = -1
```
http://genomescope.org/genomescope2.0/analysis.php?code=Ah6wDuTUPCxsKDsFqEsF

##
