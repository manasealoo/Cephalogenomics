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

## Overview of gene function at the origin of cephalochordates

This was done for the GenEra results for _B. floridae_. We first extracted genes belonging to phylostrata of 'Branchiostomidae' (aka modern Cephalochordata).

```bash
grep "Branchiostomidae" 7739_phylostrata_assignation.tsv | awk ' { print $1 } ' > GO/7739_phylostrata_assignation_CEPH.tsv
```

Using [seqtk](https://github.com/lh3/seqtk), we extracted protein sequences belonging to genes (longest isoforms) at the Branchiostomidae phylostrata.

```bash
seqtk subseq bflor.branchiostoma_floridae.filtered.longest.faa \
GO/7739_phylostrata_assignation_CEPH.tsv > GO/7739_phylostrata_assignation_CEPH.faa

sed -i.bak "s/\-\-$//" 7739_phylostrata_assignation_CEPH.faa
```

We then ran interproscan.

```bash
interproscan
-i 7739_phylostrata_assignation_CEPH.faa \
-goterms
```
We then extracted all the GO terms that are associated with each genes, creating a tab-delimited table. Josué Barerra-Redondo kindly provided me with a script to do this, `get_GO.pl`.

```perl
#!/usr/bin/perl
use strict;
use warnings;

#USAGE: perl Get_GO.pl List1 List2 > output
# List1 = non-redundant list of genes in the genome/transcriptome/metagenome
# List2 = tsv output obtained from InterProScan

open(FILE1, $ARGV[0]) || die "\nThe List1 with non-redundant gene names does not exist.\n\n";

open(FILE2, $ARGV[1]) || die "\nThe List2 obtained from InterProScan does not exist.\n\n";

my @LISTA1=<FILE1>;
chomp @LISTA1;
close (FILE1);

my @GO=<FILE2>;
chomp @GO;
close (FILE2);


foreach my $result (@LISTA1){
	
	my @query = split (/\t/, $result);
	
	my @go=();

	foreach my $linea (@GO){
		
		my @homolog = split (/\t/, $linea);	

		if ($query[0] eq $homolog[0]){		
		
		push (@go, "$homolog[13]");
	
		}

	}
	print "$query[0]\t@go\n";
}
```

Finally, we looked at the top GO hits.

```bash
grep ">" 7739_phylostrata_assignation_CEPH.faa | sed 's/>//g' | awk '{ print $1 }' > lista_genes

perl get_GO.pl lista_genes 7739_phylostrata_assignation_CEPH.faa.tsv > topGO_temp.tab

while IFS=$'\t' read -r -a LINEA; do
        GO=$(echo ${LINEA[1]} | sed 's/|/\n/g' | sed 's/ /\n/g' | sort -u | sed '/^$/d' | sed ':a;N;$!ba;s/\n/,/g')
        GEN=$(echo "${LINEA[0]}")
        echo -e ${GEN}"\t"${GO} >> gene_family_universe.tsv
done < topGO_temp.tab

rm topGO_temp.tab
```
