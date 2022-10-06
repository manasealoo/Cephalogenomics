# _Branchiostoma lanceolatum_ (from Banyuls-sur-Mer) genome-guided

Note: this north sea specimen transcriptome used the discontinued genome assembly to scaffold the transcriptome assembly.

## Transcriptome assembly

I used the reads from adults.

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

Where `Blnc_RNA_R1_trimmed.fq.gz` and `Blnc_RNA_R2_trimmed.fq.gz` are the forward (R1) and reverse (R2) trimmed RNA-seq reads. `igv_sorted.bam` is the sorted RNA-seq mapping (which was done for IGV), i.e.
```
./samtools sort input.bam -o igv_sorted.bam
```
