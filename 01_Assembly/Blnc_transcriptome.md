# _Branchiostoma lanceolatum_ (from Banyuls-sur-Mer) genome-guided

Note: this north sea specimen transcriptome used the discontinued genome assembly to scaffold the transcriptome assembly.

## Transcriptome assembly

We used the reads from adults.

### Trinity v2.12.0 genome-guided mode

```
./Trinity \
--genome_guided_bam /hps/nobackup/research/marioni/sodai/braker_results/Blnc/igv_sorted.bam \
--genome_guided_max_intron 10000 \
--seqType fq \
--left ../RNA_preprocessing/Blnc_RNA_R1_trimmed.fq.gz \
--right ../RNA_preprocessing/Blnc_RNA_R2_trimmed.fq.gz \
--CPU 6 \
--max_memory 20G \
--output /hps/nobackup/research/marioni/sodai/trinity_results/trinity_guided_out_dir
```

Where `Blnc_RNA_R1_trimmed.fq.gz` and `Blnc_RNA_R2_trimmed.fq.gz` are the forward (R1) and reverse (R2) trimmed RNA-seq reads. These were trimmed using trimmomatic, i.e.

```
trimmomatic PE \
-threads 10 \
-trimlog fastqc2.log \
../Blnc_RNA_R1_raw.fastq \
../Blnc_RNA_R2_raw.fastq \
Blnc_RNA_R1_trimmed.fq.gz \
Blnc_RNA_R1_unpaired.fq.gz \
Blnc_RNA_R2_trimmed.fq.gz \
Blnc_RNA_R2_unpaired.fq.gz \
ILLUMINACLIP:TruSeq2-PE.fa:4:30:10 \
SLIDINGWINDOW:5:15
```

`igv_sorted.bam` is the sorted RNA-seq mapping (with the prefix `igv` because the sorting was initially done for IGV), i.e.
```
./hisat2 \
-x Blnc_2_scaf \
-1 ../../RNA_preprocessing/Blnc_RNA_R1_trimmed.fq.gz \
-2 ../../RNA_preprocessing/Blnc_RNA_R2_trimmed.fq.gz \
-p 10 \
-S input.sam
```
```
./samtools view -bS input.sam > input.bam
```
then
```
./samtools sort input.bam -o igv_sorted.bam
```
