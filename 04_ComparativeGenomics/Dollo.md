# Analysis based on Dollo's parsimony on gene family gain and loss.

KinFin scripts were used on the outputs of orthofinder.

```
echo '#IDX,TAXON' > config.txt
sed 's/: /,/g' SpeciesIDs.txt | \
    cut -f 1 -d"." \
    >> config.txt
```

```
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
