# Phylogenetic analysis under the coalescence model

Here we used IQ-tree, MAFFT, phylopypruner and ASTRAL to infer the phylogenetic relationship between the three amphioxus genera.

## Multiple sequence alignment construction

MAFFT einsi was used on the hierarchical orthogroup (HOG) output from orthofinder. Node 0 was used.

```
find *.faa -type f -exec sed -i.bak "/^[^>]/s/U/X/g" {} \;
```

A job array was created to run einsi.

```
chmod u+x run_einsi.sh
./run_einsi.sh
```

Where `run_einsi.sh` consists of the following lines of code:
```
#!/usr/bin/env bash
files=(*.faa)
einsi --amino --thread 4 --quiet ${files[$LSB_JOBINDEX-1]} > ${files[$LSB_JOBINDEX-1]%.faa}.aln
```
We couldn't align one HOG (N0.HOG0000895) because it required more than 120GB of RAM.
 
## HOG phylogenetic reconstruction (gene tree)
 
IQ-tree was used to infer gene trees.
 
```
 #!/usr/bin/env bash
files=(*.aln)
iqtree -s ${files[$LSB_JOBINDEX-1]} -T 8 --seqtype AA -keep-ident --quiet

chmod u+x run_iqtree.sh
```

```
./run_iqtree.sh
```

We reran IQ-tree until gene trees were inferred.

## Quartet-based species tree reconstruction using paralogues and orthologues.

We ran ASTRAL-Pro without pruning any gene trees as follows. Note: for the main figures, the gene trees were pruned using phylopypruner before conducting the concatenation and the quartet-based analyses.

```
cp /hps/nobackup/research/marioni/sodai/phylogeny/run_7_phylogenomic_dataset_3/phylopypruner_input/*.aln.treefile .
cat /hps/nobackup/research/marioni/sodai/phylogeny/run_7_astral/pd3_input_all/*.aln.treefile > pd3_all_input_trees.tre

sed 's/|XP_[0-9]*\.[0-9]\:/\:/g' pd3_all_input_trees.tre | \
sed 's/|TRINITY_GG_[0-9]*_c[0-9]*_g[0-9]*_i[0-9]*\.p[0-9]*\:/\:/g' | \
sed 's/|TRINITY_DN[0-9]*_c[0-9]*_g[0-9]*_i[0-9]*\.p[0-9]*\:/\:/g' | \
sed 's/|g[0-9]*\.t[0-9]*\:/\:/g' | \
sed 's/|Boleac\.CG\.SB_v[0-9]*\.S[0-9]*\.g[0-9]*\.[0-9]*\.t\:/\:/g' | \
sed 's/|Eptbu[0-9]*\.t[0-9]*\:/\:/g' | \
sed 's/|NP_[0-9]*\.[0-9]\:/\:/g' > pd3_all_input_trees.cl.tr
```

Then

```
java \
-D'java.library.path=/hps/nobackup/research/marioni/sodai/A-pro/ASTRAL-MP/lib' \
-jar /hps/nobackup/research/marioni/sodai/A-pro/ASTRAL-MP/astral.1.1.5.jar \
-T 4 \
-i pd3_all_input_trees.cl.tre -o pd3_all_astral.tre 2> pd3_all_astral.lo
```

as well as 

```
java \
-D'java.library.path=/hps/nobackup/research/marioni/sodai/A-pro/ASTRAL-MP/lib' \
-jar /hps/nobackup/research/marioni/sodai/A-pro/ASTRAL-MP/astral.1.1.5.jar \
-T 4 \
-t 2 \
--exact \
-i pd3_all_input_trees.cl.tre -o pd3_all_astral_2.tre 2> pd3_all_astral_2.log
```

to get the full annotation of branch support.

## Gene-tree pruning

We used `phylopypruner` to obtain the single copy orthologes for the quartet-based ASTRAL and concatenation-based IQ-TREE analyses.

```
phylopypruner \
--dir phylopypruner_input \
--outgroup anneissia saccoglossus \
--threads 8 \
--min-len 50 \
--trim-lb 5 \
--min-support 0.75 \
--prune MI \
--min-taxa 8 \
--min-otu-occupancy 0.1 \
--min-gene-occupancy 0.1 \
--mask pdist \
--output phylopypruner_output_para
```

Where

These options were chosen to
* Remove sequences which are shorter than 50 bases
* Remove branches which are longer than 5 times the standard deviation of all branches
* Collapse nodes with less than 75% bootstrap support into polytomies
* Use the maximum inclusion (MI) algorithm for paralogy pruning; iteratively select the largest subtree which does not contain more than one sequence per taxa within that subtree
* Remove a taxon if that taxon is present within fewer than 10% of all output alignments
* Remove alignments where more than 90% of the positions consists of gaps

## Quartet-based species tree reconstruction using inferred single-copy orthologues.

ASTRAL III was used for this analysis and is the phylogenetic reconstruction (aside from the concatenation analysis) used in the main figures.

##

