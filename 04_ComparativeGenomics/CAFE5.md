# Analysis of gene family birth and death

CAFE uses a birth and death process to model gene gain and loss across a given phylogenetic tree. This analysis was done to further complement the findings of new orthogroups at the origin of cephalochordates (via Dollo's parsimony).

## Preprocessing

First, the orthogroup file was taken from `Orthogroups/Orthogroups.GeneCount.tsv` output of OrthoFinder (which was renamed to `run_7_OG.GeneCount.tsv`).

Orthogroup names were changed to `Family ID`. Then we ran

```bash
awk 'BEGIN { FS="\t"; OFS="\t" } { $1=$1 "\t" $1 } 1' run_7_OG.GeneCount.tsv > run_7_OG.GeneCount2.tsv
```

Next, we changed the name of the file to `run_7_OG.cafe_fam.tsv`. (Don't ask why, I just thought it was less confusing this way).

Using the APE package in R, we also dated and ultrametricised the concatenation tree, 

```R
y <- ape::read.tree(text = '((((bblch:0.0212143807,(bflor:0.0198708565,(Blnc2018_re:0.0080206352,BlncHG_Trinity:0.0050874146)100/100:0.0179084105)100/100:0.0067292076):0.0196086266,(EPI:0.0335587903,(ASY_Yue:0.0140559522,ASY:0.0162506137)100/100:0.0262314821)100/100:0.0094741904)100/100:0.389693067,((((latimeria:0.1033966438,lepisosteus:0.1340149553)100/100:0.0329687226,(callorhinchus:0.0760672984,amblyraja:0.0895193412)100/100:0.0634934855)100/100:0.1513918608,eptatretus:0.5187042999)100/100:0.219172228,(ciona:0.5148541349,botrylloides:0.5794070871)100/100:0.5284676027)100/100:0.1085351426)100/100:0.06412,(saccoglossus:0.3843243013,anneissia:0.5624068368)100/100:0.06412)Root;')

yt <- ape::chronos(y, model = "correlated")
```

This gives
```R
Setting initial dates...
Fitting in progress... get a first set of estimates
         (Penalised) log-lik = -9.068348 
Optimising rates... dates... -9.068348 
Optimising rates... dates... -9.068327 

log-Lik = -8.998056 
PHIIC = 106.05
```

The tree was dated with the following date at the origin of deuterostomes. (Discrete was chosen here because applying other models changed the branch length compared to the initial ultrametricisation.

```R
cal <- makeChronosCalib(yt, interactive = TRUE)
yc <- ape::chronos(yt, model = "discrete", calibration = cal)
```

where
```R
  node age.min age.max soft.bounds
1   17     530     636       FALSE
```

This gives
```R
Setting initial dates...
Fitting in progress... get a first set of estimates
         (Penalised) log-lik = -12.93759 
Optimising rates... frequencies... dates... -12.93759 
Optimising rates... frequencies... dates... -10.86958 
Optimising rates... frequencies... dates... -10.86957 
Optimising rates... frequencies... dates... -10.86956 

log-Lik = -10.86956 
PHIIC = 89.74 
```

and

```R
> is.rooted(yc)
[1] TRUE
> is.ultrametric(yc)
[1] TRUE
> is.binary(yc)
[1] TRUE
```
resulting in the following tree after removing the support values (for CAFE5 analysis). 

```R
((((bblch:26.80891516,(bflor:20.86078282,(Blnc2018_re:5.90689614,BlncHG_Trinity:5.90689614):14.95388667):5.948132345):21.63061117,(EPI:38.29041288,(ASY_Yue:14.21204082,ASY:14.21204082):24.07837206):10.14911345):473.6441633,((((latimeria:106.2193981,lepisosteus:106.2193981):31.10054534,(callorhinchus:75.14480955,amblyraja:75.14480955):62.17513387):164.9499229,eptatretus:302.2698663):167.4216333,(ciona:252.6353794,botrylloides:252.6353794):217.0561202):52.39219006):39.77890186,(saccoglossus:494.483958,anneissia:494.483958):67.37863353);
```

For others in the supplementary,

```R
ape::drop.tip()
```

was used to test the consistency of the results with removal of less reliable assemblies/annotations (e.g. transcriptome assemblies).

## Running CAFE5

To estimate lambda with no among family rate variation, we ran

```bash
/hps/software/users/marioni/sodai/CAFE5/bin/cafe5 \
-i run_7_OG.cafe_fam_filt.tsv \
-t deut_ultra_tree.txt \
-o res_2
```
The top 100 gene families with the largest difference between the max and min counts were removed as suggested in the [CAFE5 documentation](https://github.com/hahnlab/CAFE5#known-limitations).

To incorporate among family rate variation with both lambda and alpha estimated and three discrete gamma rate categories, we ran

```bash
/hps/software/users/marioni/sodai/CAFE5/bin/cafe5 \
-i run_7_OG.cafe_fam.tsv \
-t deut_ultra_tree.txt \
-k 3 \
-o res_k3_1
````

No gene families were filtered in this analysis.
