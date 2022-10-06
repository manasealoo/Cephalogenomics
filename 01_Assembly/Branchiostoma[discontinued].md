# _Branchiostoma lanceolatum_ (North Sea) genome assembly[annotation discontinued]

Note: this north sea specimen genome was discontinued for further analysis except to used to scaffold the transcriptome assembly.

## Pre-assembly analysis

### Nanopore reads

```
NanoPlot \
-t 8 \
-o ONT_stats \
-p Blnc_log \
--loglength \
--N50 \
--title "Blnc" \
--fastq Blnc_Feb20.fastq.gz
```

```
NanoPlot \
-t 8 \
-o ONT_stats \
-p Blnc \
-c Pastel1 \
--N50 \
--title "Blnc" \
--fastq Blnc_Feb20.fastq.gz
```

## Genome assembly

_B. lanceolatum_ (North Sea) Canu assembly with options to keep haplotypes
  
### Canu assembly (v2.1.1)
Using the option to keep haplotypes, `corOutCoverage=200 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50"`.
```
./canu-2.1.1/bin/canu \
-p Blnc_canu_poly1 \
-d Blnc_canu_poly1 \
genomeSize=500m \
-nanopore Blnc_Feb20.fastq.gz \
corOutCoverage=200 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50"
```
Where `Blnc_Feb20.fastq.gz` is the long reads.
