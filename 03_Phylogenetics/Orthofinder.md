# Orthogroup inference using OrthoFinder

Prior to the phylogenetic reconstruction and orthogroup analysis, we ran OrthoFinder on amino acid sequences extracted from the longest isoform of each gene.

## Extracting the longest isoform

Aside from the Trinity assemblies (from which we obtained longest isoforms), we used KinFin to extract the longest isoforms of each gene. This allows us to obtain better orthogroups to reconstruct the species tree (= phylogeny).

```
/hps/nobackup/research/marioni/sodai/kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f ../BlncHG.branchiostoma_lanceolatum.faa \
-g ../../braker/braker.gtf
```

```
../../kinfin/scripts/filter_fastas_before_clustering.py \
-f ASY.asymmetron.faa > ASY.asymmetron.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f EPI.epigonichthys.faa > EPI.epigonichthys.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f ASY_Yue.asymmetron.faa > ASY_Yue.asymmetron.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f Blnc2018_re.branchiostoma_lanceolatum.faa > Blnc2018_re.branchiostoma_lanceolatum.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f BlncHG_Trinity.branchiostoma_lanceolatum.faa > BlncHG_Trinity.branchiostoma_lanceolatum.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f amblyraja.amblyraja_radiata.faa > amblyraja.amblyraja_radiata.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f callorhinchus.callorhinchus_milii.faa > callorhinchus.callorhinchus_milii.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f ciona.ciona_intestinalis.faa > ciona.ciona_intestinalis.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f eptatretus.eptatretus_burgeri.faa > eptatretus.eptatretus_burgeri.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f latimeria.latimeria_chalumnae.faa > latimeria.latimeria_chalumnae.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f lepisosteus.lepisosteus_oculatus.faa > lepisosteus.lepisosteus_oculatus.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f anneissia.anneissia_japonica.faa > anneissia.anneissia_japonica.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f saccoglossus.saccoglossus_kowalevskii.faa > saccoglossus.saccoglossus_kowalevskii.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f botrylloides.botrylloides_leachii.faa > botrylloides.botrylloides_leachii.filtered.faa


../../kinfin/scripts/filter_fastas_before_clustering.py \
-f bflor.branchiostoma_floridae.faa > bflor.branchiostoma_floridae.filtered.faa

../../kinfin/scripts/filter_fastas_before_clustering.py \
-f bblch.branchiostoma_belcheri.faa > bblch.branchiostoma_belcheri.filtered.faa


ls *filtered.faa


cp ../../braker_results/ASY/braker/braker.gtf ASY.asymmetron.gtf

sed -i -e 's/file_1_file_1_//g' ASY.asymmetron.gtf
sed -i -e 's/file_1_file_2_//g' ASY.asymmetron.gtf


bsub -M 8000 -n 4 -o ASY_gtf2gff.out -e ASY_gtf2gff.err \
"gtf2gff.pl < ASY.asymmetron.gtf --out=ASY.asymmetron.gff3 --gff3"

sed -i -e 's/_1//g' ASY.asymmetron.filtered.faa

sed -i -e 's/_1//g' ASY.asymmetron.faa

# AUGUSTUS
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f ASY.asymmetron.filtered.faa \
-g ASY.asymmetron.gff3 \
-t AUGUSTUS

cp ../../braker_results/EPI/braker/braker.gtf EPI.epigonichthys.gtf

sed -i -e 's/file_1_file_1_//g' EPI.epigonichthys.gtf
sed -i -e 's/file_1_file_2_//g' EPI.epigonichthys.gtf

bsub -M 8000 -n 4 -o EPI_gtf2gff.out -e EPI_gtf2gff.err \
"gtf2gff.pl < EPI.epigonichthys.gtf --out=EPI.epigonichthys.gff3 --gff3"





sed -i -e 's/_1//g' EPI.epigonichthys.filtered.faa

sed -i -e 's/_1//g' EPI.epigonichthys.faa

# AUGUSTUS
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f EPI.epigonichthys.filtered.faa \
-g EPI.epigonichthys.gff3 \
-t AUGUSTUS


cp ../../braker_results/Blnc2018/ET_adult_embryo/braker/braker.gtf Blnc2018_re.branchiostoma_lanceolatum.gtf

sed -i -e 's/file_1_file_1_//g' Blnc2018_re.branchiostoma_lanceolatum.gtf
sed -i -e 's/file_1_file_2_//g' Blnc2018_re.branchiostoma_lanceolatum.gtf

conda activate braker2-env

bsub -M 8000 -n 4 -o Blnc2018_re_gtf2gff.out -e Blnc2018_re_gtf2gff.err \
"gtf2gff.pl < Blnc2018_re.branchiostoma_lanceolatum.gtf --out=Blnc2018_re.branchiostoma_lanceolatum.gff3 --gff3"


(not done yet)

sed -i -e 's/_1//g' Blnc2018_re.branchiostoma_lanceolatum.faa

sed -i -e 's/_1//g' Blnc2018_re.branchiostoma_lanceolatum.filtered.faa


conda activate kinfin_env

# AUGUSTUS
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f Blnc2018_re.branchiostoma_lanceolatum.filtered.faa \
-g Blnc2018_re.branchiostoma_lanceolatum.gff3 \
-t AUGUSTUS



cp /hps/nobackup/research/marioni/sodai/trinity_results/trinity_Yue_adult_out_dir/Trinity.fasta.transdecoder_dir/longest_orfs.pep ASY_Yue.asymmetron.faa

cp /hps/nobackup/research/marioni/sodai/trinity_results/trinity_guided_out_dir/Trinity-GG.fasta.transdecoder_dir/longest_orfs.pep BlncHG_Trinity.branchiostoma_lanceolatum.faa


# Transdecoder
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f ASY_Yue.asymmetron.filtered.faa \
-t custom \
--fdel _ \
--fs 1 \
--fe 3

# Transdecoder
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f BlncHG_Trinity.branchiostoma_lanceolatum.filtered.faa \
-t custom \
--fdel _ \
--fs 1 \
--fe 3

Must download gff files for NCBI taxa

mkdir gffs

# NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f amblyraja.amblyraja_radiata.filtered.faa \
-g gffs/amblyraja.amblyraja_radiata.gff \
-t NCBI

# NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f callorhinchus.callorhinchus_milii.filtered.faa \
-g gffs/callorhinchus.callorhinchus_milii.gff \
-t NCBI

# NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f ciona.ciona_intestinalis.filtered.faa \
-g gffs/ciona.ciona_intestinalis.gff \
-t NCBI

# AUGUSTUS
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f eptatretus.eptatretus_burgeri.filtered.faa \
-g gffs/eptatretus.eptatretus_burgeri.gff \
-t AUGUSTUS

# NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f latimeria.latimeria_chalumnae.filtered.faa \
-g gffs/latimeria.latimeria_chalumnae.gff \
-t NCBI

# NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f lepisosteus.lepisosteus_oculatus.filtered.faa \
-g gffs/lepisosteus.lepisosteus_oculatus.gff \
-t NCBI

# NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f anneissia.anneissia_japonica.filtered.faa \
-g gffs/anneissia.anneissia_japonica.gff \
-t NCBI

# NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f saccoglossus.saccoglossus_kowalevskii.filtered.faa \
-g gffs/saccoglossus.saccoglossus_kowalevskii.gff \
-t NCBI

# AUGUSTUS
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f botrylloides.botrylloides_leachii.filtered.faa \
-g gffs/botrylloides.botrylloides_leachii.gff \
-t AUGUSTUS

[-] None of the GFF protein IDs seems to be in FASTA
Example of IDs in GFF:
Boleac.CG.SB_v3.S0.g00001.01.t,Boleac.CG.SB_v3.S0.g00002.01.t,Boleac.CG.SB_v3.S0.g00003.01.t,Boleac.CG.SB_v3.S0.g00004.01.t,Boleac.CG.SB_v3.S0.g00005.01.t,Boleac.CG.SB_v3.S0.g00006.01.t,Boleac.CG.SB_v3.S0.g00007.01.t,Boleac.CG.SB_v3.S0.g00008.01.t,Boleac.CG.SB_v3.S0.g00009.01.t,Boleac.CG.SB_v3.S0.g00010.01.t


nano test2.fasta
sed '/^>/s/p$/t/' test2.fasta > output2.fasta
works!

mv botrylloides.botrylloides_leachii.filtered.faa botrylloides.botrylloides_leachii.filtered.pre.faa
mv botrylloides.botrylloides_leachii.faa botrylloides.botrylloides_leachii.pre.faa


sed '/^>/s/p$/t/' botrylloides.botrylloides_leachii.filtered.pre.faa > botrylloides.botrylloides_leachii.filtered.faa

sed '/^>/s/p$/t/' botrylloides.botrylloides_leachii.pre.faa > botrylloides.botrylloides_leachii.faa


redoing 
# AUGUSTUS
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f botrylloides.botrylloides_leachii.filtered.faa \
-g gffs/botrylloides.botrylloides_leachii.gff \
-t AUGUSTUS

WORKS!

# NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f bflor.branchiostoma_floridae.filtered.faa \
-g gffs/bflor.branchiostoma_floridae.gff \
-t NCBI

# NCBI
../../kinfin/scripts/filter_isoforms_based_on_gff3.py \
-f bblch.branchiostoma_belcheri.filtered.faa \
-g gffs/bblch.branchiostoma_belcheri.gff \
-t NCBI


mkdir longest_isoforms


grep -A1 -wFf amblyraja.longest_isoforms.txt \
amblyraja.amblyraja_radiata.filtered.faa \
> longest_isoforms/amblyraja.amblyraja_radiata.filtered.longest.faa

grep -A1 -wFf anneissia.longest_isoforms.txt \
anneissia.anneissia_japonica.filtered.faa \
> longest_isoforms/anneissia.anneissia_japonica.filtered.longest.faa

grep -A1 -wFf ASY.longest_isoforms.txt \
ASY.asymmetron.filtered.faa \
> longest_isoforms/ASY.asymmetron.filtered.longest.faa

grep -A1 -wFf ASY_Yue.longest_isoforms.txt \
ASY_Yue.asymmetron.filtered.faa \
> longest_isoforms/ASY_Yue.asymmetron.filtered.longest.faa

grep -A1 -wFf bblch.longest_isoforms.txt \
bblch.branchiostoma_belcheri.filtered.faa \
> longest_isoforms/bblch.branchiostoma_belcheri.filtered.longest.faa

grep -A1 -wFf bflor.longest_isoforms.txt \
bflor.branchiostoma_floridae.filtered.faa \
> longest_isoforms/bflor.branchiostoma_floridae.filtered.longest.faa

grep -A1 -wFf Blnc2018_re.longest_isoforms.txt \
Blnc2018_re.branchiostoma_lanceolatum.filtered.faa \
> longest_isoforms/Blnc2018_re.branchiostoma_lanceolatum.filtered.longest.faa

grep -A1 -wFf BlncHG_Trinity.longest_isoforms.txt \
BlncHG_Trinity.branchiostoma_lanceolatum.filtered.faa \
> longest_isoforms/BlncHG_Trinity.branchiostoma_lanceolatum.filtered.longest.faa

grep -A1 -wFf botrylloides.longest_isoforms.txt \
botrylloides.botrylloides_leachii.filtered.faa \
> longest_isoforms/botrylloides.botrylloides_leachii.filtered.longest.faa

grep -A1 -wFf callorhinchus.longest_isoforms.txt \
callorhinchus.callorhinchus_milii.filtered.faa \
> longest_isoforms/callorhinchus.callorhinchus_milii.filtered.longest.faa

grep -A1 -wFf ciona.longest_isoforms.txt \
ciona.ciona_intestinalis.filtered.faa \
> longest_isoforms/ciona.ciona_intestinalis.filtered.longest.faa

grep -A1 -wFf EPI.longest_isoforms.txt \
EPI.epigonichthys.filtered.faa \
> longest_isoforms/EPI.epigonichthys.filtered.longest.faa

grep -A1 -wFf eptatretus.longest_isoforms.txt \
eptatretus.eptatretus_burgeri.filtered.faa \
> longest_isoforms/eptatretus.eptatretus_burgeri.filtered.longest.faa

grep -A1 -wFf latimeria.longest_isoforms.txt \
latimeria.latimeria_chalumnae.filtered.faa \
> longest_isoforms/latimeria.latimeria_chalumnae.filtered.longest.faa

grep -A1 -wFf lepisosteus.longest_isoforms.txt \
lepisosteus.lepisosteus_oculatus.filtered.faa \
> longest_isoforms/lepisosteus.lepisosteus_oculatus.filtered.longest.faa

grep -A1 -wFf saccoglossus.longest_isoforms.txt \
saccoglossus.saccoglossus_kowalevskii.filtered.faa \
> longest_isoforms/saccoglossus.saccoglossus_kowalevskii.filtered.longest.faa


Here, we have to run OrthoFinder on this file

/hps/nobackup/research/marioni/sodai/kinfin_results/run_6_input/longest_isoforms

mkdir run_6_analysis
```

## Running OrthoFinder
