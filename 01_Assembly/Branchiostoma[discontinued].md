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

## Post-assembly processing

### Polishing using Racon v1.4.3

Mapping long reads onto the assembly using Minimap2

```
./minimap2 \
-x map-ont \
-t 12 \
asm.fasta \
../../Blnc_Feb20.fastq.gz | gzip -c - > asm.paf.gz
```

Where `asm.fasta` is the assembly and `Blnc_Feb20.fastq.gz` is the long reads.

Racon with default parameters and raw reads
```
./racon \
-t 15 \
../../Blnc_Feb20.fastq.gz \
asm.paf.gz \
asm.fasta > racon.fasta
```

Where `asm.fasta` is the assembly, `Blnc_Feb20.fastq.gz` is the long reads and `asm.paf.gz` is the alignment file.

### Haplotig purging

#### Mapping long reads using Minimap2

Map long reads
```
minimap2 -x map-ont -t 12 asm.fasta /hps/nobackup/research/marioni/sodai/Blnc_canu_poly1/Blnc_canu_poly1.trimmedReads.fasta.gz | gzip -c - > asm.paf.gz
```

Where asm.fasta is the polished assembly and `/hps/nobackup/research/marioni/sodai/Blnc_canu_poly1/Blnc_canu_poly1.trimmedReads.fasta.gz` is the full path to the Canu trimmed reads.

#### purge_dups v1.2.5

Calculate read depth histogram

```
./purge_dups/bin/pbcstat asm.paf.gz
```

Where `asm.paf.gz` is the output from the alignment step.

Calculate base-level read depth

```
./purge_dups/bin/calcuts -l 5 -m 28 -u 120 PB.stat > cutoffs 2>calcults.log
```

Here, I set the custom cutoff. The custom cutoffs are `5 27 27 28 28 120`

Split an assembly and do a self-self alignment

```
./purge_dups/bin/split_fa asm.fasta > asm.split
```

```
minimap2 -xasm20 -DP asm.split asm.split | gzip -c - > asm.split.self.paf.gz
```

Where `asm.split` is the split assembly.

Purge haplotigs and overlaps

```
./purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log
```

Where `cutoffs` is a file containing the manually calculated cutoffs.

Get purged primary and haplotig sequences

```
./purge_dups/bin/get_seqs dups.bed asm.fasta
```

Where `.bed` file `dups.bed` contains the coordinates for purging. Notice, `-e` was not included.
