# Gene age inference based on a phylostratigraphic framework

To assess further whether the gene-space novelty as detected in the emergence of new orthogroups presents a true biological signal, we ran an orthogonal analysis using [GenEra](https://www.biorxiv.org/content/10.1101/2022.07.07.498977v1.abstract).

We installed made a diamond database for the NCBI-NR according to the [documentation](https://github.com/josuebarrera/GenEra).
```bash
diamond makedb \
 --in nr \
 --db nr \
 --taxonmap prot.accession2taxid \
 --taxonnodes taxdump/nodes.dmp \
 --taxonnames taxdump/names.dmp \
 --memory-limit 100
```

We ran GenEra on the longest gene isoforms of B. floridae (taxid 7739), Epigonichthys (taxid 231028) and Asymmetron (taxid 223987). For example, for Epigonichthys
```bash
# EPI
genEra \
-q /nfs/research/marioni/sodai/FINAL_PROTEOME/EPI.epigonichthys.filtered.longest.faa \
-t 231028 \
-b nr \
-d taxdump \
-a CEPH/CEPH_proteins.tsv
```

Additional cephalochordate proteins from our assemblies/annotations were included in `-a CEPH/CEPH_proteins.tsv`.
