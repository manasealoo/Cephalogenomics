# _Epigonichthys_ genome annotation

## BRAKER ET mode

### Align RNA-seq reads using Hisat2

Index the draft genome

```
hisat2-build \
EPI_2_scaf.fa \
EPI_2_scaf
```
Where `EPI_2_scaf.fa` is the scaffolded assembly (repeat-masked, haplotig-purged and polished).

Map the RNA-seq reads

```
hisat2 \
-x EPI_2_scaf \
-1 ../../RNA_preprocessing/EPI_RNA_R1_trimmed.fq.gz \
-2 ../../RNA_preprocessing/EPI_RNA_R2_trimmed.fq.gz \
-p 10 \
-S input.sam
```

Where `EPI_RNA_R1_trimmed.fq.gz` and `EPI_RNA_R2_trimmed.fq.gz` are the forward (R1) and reverse (R2) trimmed RNA-seq reads. For trimming protocol, see the [RNA-scaffolding step](https://github.com/LotharukpongJS/Cephalogenomics/blob/main/01_Assembly/Epigonichthys.md#rna-scaffolding).

Convert `.sam` to `.bam`

```
samtools view -bS input.sam > input.bam
```
Sort

```
samtools sort -n -@ 4 -m 2G input.bam -o input.sorted.bam
```

### BRAKER ET mode (BRAKER v2.1.6)

```
braker.pl \
--species=EPI \
--genome=EPI_2_scaf.fa \
--bam=input.sorted.bam \
--softmasking \
--cores=20 \
--useexisting
```

Where `input.sorted.bam` is the RNA-seq read alignment file and `EPI_2_scaf.fa` is the scaffolded assembly.  I used `--useexisting` because the AUGUSTUS training parameter file `EPI` already exists from a previously failed run.

## Protein sequence extraction

### Get CDS using getAnnoFasta.pl (BRAKER v2.1.6)

```
getAnnoFasta.pl ../braker/braker.gtf --seqfile ../EPI_2_scaf.fa
```

Where `EPI_2_scaf.fa` is the scaffolded genome and `braker.gtf` is the BRAKER annotation file.

### Use transeq from EMBOSS to translate CDS

```
transeq ../braker/braker.codingseq EPI.faa
```

Where `braker.codingseq` is the output of the previous step.
