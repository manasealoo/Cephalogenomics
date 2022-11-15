# Analysis based on Dollo's parsimony on gene family gain and loss.

KinFin scripts were used on the outputs of orthofinder.

```bash
echo '#IDX,TAXON' > config.txt
sed 's/: /,/g' SpeciesIDs.txt | \
    cut -f 1 -d"." \
    >> config.txt
```

```bash
xvfb-run ../../kinfin/kinfin \
--cluster_file Orthogroups.txt \
--config_file config.txt \
--sequence_ids_file SequenceIDs.txt \
--species_ids_file SpeciesIDs.txt \
--tree_file SpeciesTree_rooted.txt \
--fasta_dir longest_isoforms \
--plot_tree \
--outprefix run_7
```

Where `Orthogroups.txt`, `SequenceIDs.txt`, `SpeciesIDs.txt` and `SpeciesTree_rooted.txt` were obtained from the orthofinder run.

## Checking the RNA-seq support for the gene models

To test whether the new genes we see are directly supported by RNA-seq evidence we used a script in the BRAKER pipeline `selectsupportedsubsets.py` and used the `braker.gtf` and `hintsfile.gff`.

```bash
./selectSupportedSubsets.py --fullSupport fullSupport2.tsv --anySupport anySupport2.tsv --noSupport noSupport2.tsv braker.gtf hintsfile.gff
```

After processing the results, i.e. in _B. lanceolatum_ 2018 reannotation
```bash
for i in *Support2.tsv; do
grep "gene_id" $i | awk -F "\t" ' { print $9 } ' | awk -F ";" ' { print $1 } ' | uniq | sed 's/.*_//' | sed 's/"//' > ${i%.tsv}_geneID.tsv
done

for i in *.tsv; do
sort -u $i > ${i%.tsv}.sort.tsv
done

# Extract Blnc2018_re lines from Cephalochordate-restricted HOGs
grep "Blnc2018_re" ~/ceph_OG.txt | awk ' { print $1 "\t" $3 } ' > ~/braker_support/ceph_OG_BLNC2018_cl.tsv

cd ~/braker_support

grep "Blnc2018_re\." Orthogroups.tsv | awk -F "\t" ' { print $1 "\t" $4 } ' > OG_BLNC2018.tsv

grep -wFf BLNC2018/anySupport2_geneID.sort.tsv OG_BLNC2018.tsv | awk ' { print $1 } ' | sed 's/$/\t1/' > OG_BLNC2018_any.tsv
```
The results were plotted in R.

## Checking for hits elsewhere in the tree of life

Checking whether the cephalochordate-exclusive orthogroups have hits elsewhere in the tree of life.

```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
```

These files were concatenated as `uniprot_all.fasta`. Then, we made a diamond database.

```bash
diamond makedb --in uniprot_all.fasta -d uniprot_all
```

For complete orthogroups of _B. floridae_

```bash
diamond blastp --ultra-sensitive \
-d uniprot/uniprot_all \
-q new_OG_complete/new_OG_n6_complete.fa \
-o diamond_out_complete.tsv \
--outfmt 6 qseqid sseqid pident qcovhsp evalue \
--evalue 1e-5 \
-k0 \
--threads 16
```

Note: `-k0` was chosen to allow for unlimited hits.

The diamond output was processed to search for non-cephalochordate hits.

```bash
cat diamond_out_complete.tsv | grep -v '_BRAFL' | grep -v '_BRABE' | grep -v '_BRACL' | grep -v '_BRALA' | grep -v '_9BRAN' | \
grep 'bflor.XP' | awk '$3 >= 30' | awk '$4 >= 70' | \
sed 's/_.*//' | sort -u | wc -l
```

And for partial orthogroups of _B. floridae_

```bash
diamond blastp --ultra-sensitive \
-d uniprot/uniprot_all \
-q new_OG_partial/new_OG_n6_partial.fa \
-o diamond_out_partial.tsv \
--outfmt 6 qseqid sseqid pident qcovhsp evalue \
--evalue 1e-5 \
-k0 \
--threads 16
```

...similarly

```bash
cat diamond_out_partial.tsv | grep -v '_BRAFL' | grep -v '_BRABE' | grep -v '_BRACL' | grep -v '_BRALA' | grep -v '_9BRAN' | \
grep 'bflor.XP' | awk '$3 >= 30' | awk '$4 >= 70' | \
sed 's/_.*//' | sort -u | wc -l
```

The output numbers were compared to the number of complete and partial orthogroups that contained _B. floridae_.

Other thresholds were tested, i.e. by changing `awk '$3 >= {$pident_threshold}'` and `awk '$4 >= {$qcovhsp_threshold}'` as well as `awk '$5 >= {$evalue_threshold}'`. These results are in the supplementary.
