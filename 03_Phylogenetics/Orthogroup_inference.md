# Orthogroup inference using OrthoFinder

Prior to the phylogenetic reconstruction and orthogroup analysis, we ran OrthoFinder on amino acid sequences extracted from the longest isoform of each gene. However, before that, we had to preprocess the protein sequences.

## Extracting the longest isoform

Aside from the Trinity assemblies (from which we obtained longest isoforms), we used KinFin to extract the longest isoforms of each gene. This allows us to obtain better orthogroups to reconstruct the species tree (= phylogeny). Scripts from KinFin were used for this process.

`filter_fastas_before_clustering.py` script was used to 'sanitise' the protein sequences, printing single line fasta sequences and removing proteins shorter than 30 aa. This script was performed on all proteome files, as shown in this example.

```
../../kinfin/scripts/filter_fastas_before_clustering.py \
-f ASY.asymmetron.faa > ASY.asymmetron.filtered.faa
```

From these 'sanitised' protein sequences, we extracted the longest isoforms. Most protein sequences (for species: _Branchiostoma belcheri_
_Branchiostoma floridae_, _Saccoglossus kowalevskii_, _Anneissia japonica_, _Lepisosteus oculatus_, _Latimeria chalumnae_, _Ciona_intestinalis_, _Callorhinchus milii_ and _Amblyraja radiata_) came from NCBI and these were processed in the following manner.

```
 # NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f bflor.branchiostoma_floridae.filtered.faa \
-g gffs/bflor.branchiostoma_floridae.gff \
-t NCBI
```
 
Protein sequences inferred from our genome annotation (i.e. using BRAKER; _Asymmetron_, _Epigonichthys_ and _B. lanceolatum 2018 version_) were processed as follows. The `-t AUGUSTUS` mode was also used for _Eptatretus burgeri_.

```
# AUGUSTUS
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f ASY.asymmetron.filtered.faa \
-g ASY.asymmetron.gff3 \
-t AUGUSTUS
```

Protein sequences inferred from our transcriptome assemblies (_B. lanceolatum_ and _Asymmetron_ from Yue 2014) were processed as follows.

```
# Transdecoder
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f ASY_Yue.asymmetron.filtered.faa \
-t custom \
--fdel _ \
--fs 1 \
--fe 3
```

Some sources of protein sequences had to be further processed before the `filter_isoforms_based_on_gff3.py` could be applied.

```
sed '/^>/s/p$/t/' botrylloides.botrylloides_leachii.filtered.pre.faa > botrylloides.botrylloides_leachii.filtered.faa
```
Then
```
# AUGUSTUS
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f botrylloides.botrylloides_leachii.filtered.faa \
-g gffs/botrylloides.botrylloides_leachii.gff \
-t AUGUSTUS
```

Using the `.txt` files containing the longest isoforms, the resulting fasta files were obtained. This is done according to the KinFin documentation.
```
grep -A1 -wFf bblch.longest_isoforms.txt \
bblch.branchiostoma_belcheri.filtered.faa \
> longest_isoforms/bblch.branchiostoma_belcheri.filtered.longest.faa
```

For the transcriptome assemblies, an additional script had to be fun.

```
sed '/^[^>]/s/*$//' < BlncHG_Trinity.branchiostoma_lanceolatum.filtered.longest_2.faa > BlncHG_Trinity.branchiostoma_lanceolatum.filtered.longest.faa
```

## Running OrthoFinder

Using the longest isoforms, we ran OrthoFinder.

```
orthofinder \
-f . \
-t 64
```

Where the current directory contains the fasta files with the longest isofroms. An automatically outputted species tree was created using STAG (which is one of the figures in the supplement). To obtain OrthoFinder's automatic concatenation analysis, `-M msa` was used. Note: these figures in the supplement do not form the main arguments for the revision of the amphioxus phylogeny. For this, please look at the ASTRAL and IQTREE analyses.
