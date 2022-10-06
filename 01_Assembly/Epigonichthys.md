
# Epigonichthys genome assembly

## Pre-assembly analysis

```
../jellyfish/bin/jellyfish count -C -m 21 -s 1000M -t 10 EPI_*.fastq -o EPI_both.jf
../jellyfish/bin/jellyfish histo -t 10 EPI_both.jf > EPI_both.histo
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
http://qb.cshl.edu/genomescope/genomescope2.0/analysis.php?code=Y5BCeyMIUGNPBt1q07DK

## Genome assembly

```

```
