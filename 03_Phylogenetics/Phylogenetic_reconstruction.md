# Phylogenetic analysis under the coalescence model

Here we used IQ-tree, MAFFT, phylopypruner and ASTRAL to infer the phylogenetic relationship between the three amphioxus genera.

## Multiple sequence alignment construction

MAFFT einsi was used on the hierarchical orthogroup (HOG) output from orthofinder. Node 0 was used.

```bash
find *.faa -type f -exec sed -i.bak "/^[^>]/s/U/X/g" {} \;
```

A job array was created to run einsi.

```bash
chmod u+x run_einsi.sh
./run_einsi.sh
```

Where `run_einsi.sh` consists of the following lines of code:
```bash
#!/usr/bin/env bash
files=(*.faa)
einsi --amino --thread 4 --quiet ${files[$LSB_JOBINDEX-1]} > ${files[$LSB_JOBINDEX-1]%.faa}.aln
```
We couldn't align one HOG (N0.HOG0000895) because it required more than 120GB of RAM.
 
## HOG phylogenetic reconstruction (gene tree)
 
IQ-tree was used to infer gene trees.
 
```bash
 #!/usr/bin/env bash
files=(*.aln)
iqtree -s ${files[$LSB_JOBINDEX-1]} -T 8 --seqtype AA -keep-ident --quiet -B 1000

chmod u+x run_iqtree.sh
```

```bash
./run_iqtree.sh
```

We reran IQ-tree until gene trees were inferred.

## Quartet-based species tree reconstruction using paralogues and orthologues.

We ran ASTRAL-Pro without pruning any gene trees as follows. Note: for the main figures, the gene trees were pruned using phylopypruner before conducting the concatenation and the quartet-based analyses.

```bash
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

```bash
java \
-D'java.library.path=/hps/nobackup/research/marioni/sodai/A-pro/ASTRAL-MP/lib' \
-jar /hps/nobackup/research/marioni/sodai/A-pro/ASTRAL-MP/astral.1.1.5.jar \
-T 4 \
-i pd3_all_input_trees.cl.tre -o pd3_all_astral.tre 2> pd3_all_astral.lo
```

as well as 

```bash
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

```bash
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

ASTRAL III was used for this analysis (and is the phylogenetic reconstruction presented in the main figures along with the concatenation approach).

Gene trees were reconstructed from the pruned HOGs (phylopypruner outputs) using IQ-tree. Note: the phylopypruner output was renamed `PPP_aln_out`.

We created a file `run_iqtree_BS.sh`. (BS stands for bootstrap not what you might be thinking).
```bash
#!/usr/bin/env bash
files=(*.fasta)
iqtree -s ${files[$LSB_JOBINDEX-1]} -T 16 --seqtype AA -keep-ident --quiet -B 1000
```

```bash
chmod u+x run_iqtree_BS.sh
```

```bash
./run_iqtree.sh
```

Now, run ASTRAL

```bash
cat /hps/nobackup/research/marioni/sodai/phylogeny/run_7_phylogenomic_dataset_3/PPP_BS_aln_out/*treefile > PPP_BS_input_trees.tre
```

From the resulting trees, we removed the names of genes and only kept the names of the taxon identitfier.

```bash
sed 's/|XP_[0-9]*\.[0-9]\:/\:/g' PPP_BS_input_trees.tre | \
sed 's/|TRINITY_GG_[0-9]*_c[0-9]*_g[0-9]*_i[0-9]*\.p[0-9]*\:/\:/g' | \
sed 's/|TRINITY_DN[0-9]*_c[0-9]*_g[0-9]*_i[0-9]*\.p[0-9]*\:/\:/g' | \
sed 's/|g[0-9]*\.t[0-9]*\:/\:/g' | \
sed 's/|Boleac\.CG\.SB_v[0-9]*\.S[0-9]*\.g[0-9]*\.[0-9]*\.t\:/\:/g' | \
sed 's/|Eptbu[0-9]*\.t[0-9]*\:/\:/g' | \
sed 's/|NP_[0-9]*\.[0-9]\:/\:/g' > PPP_BS_input_trees.cl.tre
```

We then ran ASTRAL...

```bash
java -jar ../../Astral/astral.5.7.7.jar -i PPP_BS_input_trees.cl.tre -o PPP_BS_astral.tre 2> PPP_BS_astral.log
```

which resulted in the following tree

```bash
(anneissia,(saccoglossus,(((botrylloides,ciona)1:3.40147366291824,(eptatretus,((lepisosteus,latimeria)1:0.5110667093105757,(amblyraja,callorhinchus)1:2.044120473262542)1:1.769062324635397)1:1.7671768429741073)1:0.6250968033108345,((bblch,(bflor,(BlncHG_Trinity,Blnc2018_re)1:2.6570162975980915)1:0.4991383151300315)1:1.083883665366072,((ASY_Yue,ASY)1:1.949768020023296,EPI)1:0.17663702610154278)1:4.8248219636106215)1:0.7931137480077383):0.0);
```

Running ASTRAL with `-t 1` resulted in:

```bash
(anneissia,(saccoglossus,(((botrylloides,ciona)97.81:3.40147366291824,(eptatretus,((lepisosteus,latimeria)60.03:0.5110667093105757,(amblyraja,callorhinchus)91.39:2.044120473262542)88.66:1.769062324635397)88.64:1.7671768429741073)64.34:0.6250968033108345,((bblch,(bflor,(BlncHG_Trinity,Blnc2018_re)95.34:2.6570162975980915)59.54:0.4991383151300315)77.46:1.083883665366072,((ASY_Yue,ASY)90.53:1.949768020023296,EPI)44.14:0.17663702610154278)99.49:4.8248219636106215)69.86:0.7931137480077383):0.0);
```

and running ASTRAL with `-t 2` resulted in:

```bash
(anneissia,(saccoglossus,(((botrylloides,ciona)'[q1=0.9781441906165563;q2=0.011499017440432325;q3=0.010356791943011544;f1=2654.683333333334;f2=31.20833333333333;f3=28.10833333333333;pp1=1.0;pp2=0.0;pp3=0.0;QC=45;EN=2714.0]':3.40147366291824,(eptatretus,((lepisosteus,latimeria)'[q1=0.6002647053824846;q2=0.28080361948140986;q3=0.1189316751361058;f1=2140.54393939394;f2=1001.3457070707075;f3=424.1103535353533;pp1=1.0;pp2=3.68172050811551E-232;pp3=0.0;QC=24;EN=3566.0]':0.5110667093105757,(amblyraja,callorhinchus)'[q1=0.9139380341311863;q2=0.046229710035680174;q3=0.039832255833133795;f1=3122.926262626264;f2=157.96691919191915;f3=136.10681818181817;pp1=1.0;pp2=0.0;pp3=0.0;QC=24;EN=3417.0]':2.044120473262542)'[q1=0.8866259182736456;q2=0.05399498229043683;q3=0.05937909943591762;f1=2730.8078282828283;f2=166.30454545454543;f3=182.88762626262627;pp1=1.0;pp2=0.0;pp3=0.0;QC=44;EN=3080.0]':1.769062324635397)'[q1=0.8864311049425553;q2=0.04294134088981418;q3=0.07062755416763053;f1=2554.6944444444443;f2=123.75694444444447;f3=203.54861111111117;pp1=1.0;pp2=0.0;pp3=0.0;QC=72;EN=2882.0]':1.7671768429741073)'[q1=0.6433822091886608;q2=0.2012040916073174;q3=0.15541369920402176;f1=2193.9333333333334;f2=686.1059523809523;f3=529.9607142857142;pp1=1.0;pp2=0.0;pp3=0.0;QC=140;EN=3410.0]':0.6250968033108345,((bblch,(bflor,(BlncHG_Trinity,Blnc2018_re)'[q1=0.9534141776934929;q2=0.020189907266261;q3=0.026395915040246947;f1=4894.828388278393;f2=103.65498390498398;f3=135.51662781662782;pp1=1.0;pp2=0.0;pp3=0.0;QC=13;EN=5134.0]':2.6570162975980915)'[q1=0.5954129072497905;q2=0.1959017362115137;q3=0.20868535653869633;f1=3075.3076659451676;f2=1011.8324675324682;f3=1077.8598665223665;pp1=1.0;pp2=0.0;pp3=0.0;QC=24;EN=5165.0]':0.4991383151300315)'[q1=0.7746303597999609;q2=0.116844022992859;q3=0.10852561720718039;f1=4004.0643298059977;f2=603.9667548500881;f3=560.9689153439155;pp1=1.0;pp2=0.0;pp3=0.0;QC=81;EN=5169.0]':1.083883665366072,((ASY_Yue,ASY)'[q1=0.9053051892400363;q2=0.050761962384084675;q3=0.043932848375883723;f1=4640.594400044426;f2=260.205819180818;f3=225.19978077477998;pp1=1.0;pp2=0.0;pp3=0.0;QC=13;EN=5126.0]':1.949768020023296,EPI)'[q1=0.4413628148096929;q2=0.3171582711012769;q3=0.2414789140890303;f1=2279.638938492064;f2=1638.1224702380953;f3=1247.2385912698414;pp1=1.0;pp2=5.5195361057352476E-58;pp3=0.0;QC=72;EN=5165.0]':0.17663702610154278)'[q1=0.9948808409131442;q2=0.0030374977706438375;q3=0.0020816613162118777;f1=4250.130952380952;f2=12.976190476190474;f3=8.892857142857142;pp1=1.0;pp2=0.0;pp3=0.0;QC=168;EN=4272.0]':4.8248219636106215)'[q1=0.6985599780947401;q2=0.18662036691310088;q3=0.11481965499215893;f1=2672.6904761904757;f2=714.009523809524;f3=439.30000000000007;pp1=1.0;pp2=0.0;pp3=0.0;QC=49;EN=3826.0]':0.7931137480077383):0.0);
```

## Concatenation analysis

Pruned HOGs (phylopypruner outputs as detailed in the previous section) were used to construct the maximum likelihood species tree based on concatenation analysis.

The concatenated alignment (supermatrix), including the partition file for each gene, was reduced in size while retaining phylogenetic signal using MARE v0.1.2. This required the `supermatrix.fas` and `partition_data.txt` outputs of phylopypruner.

But to do this, we had to preprocess the input to MARE.

```bash
sed 's/AUTO\,/charset/' partition_data_MARE_1.txt | sed 's/$/ ;/' | sed 's/-/ - /' | sed 's/\.aln//' | sed ':a;N;$!ba;s/$/\n/g' > partition_data_MARE_2.cl.txt
```

where `partition_data_MARE_1.txt` was the renamed `partition_data.txt` file.

Wherefrom we ran MARE.

```bash
MARE \
../partition_data_MARE_2.cl.txt \
../supermatrix.fas \
-t 100 \
-d 30.0 > logfile.txt
```

The output of MARE was then processed to input into IQ-tree.

```bash
sed 's/;//' partition_PPP_BS_MARE_t100_d30.txt | sed 's/charset/charset,/' | sed 's/ - /-/' | sed 's/charset/AUTO/' > partition_PPP_BS_MARE_t100_d30.super_cl.txt
```

We then finally ran IQ-tree.

```bash
iqtree \
-s supermatrix_PPP_BS_MARE_t100_d30.fasta \
-p partition_PPP_BS_MARE_t100_d30.super_cl.txt \
-m MFP-MERGE \
-T 40 \
--seqtype AA \
-alrt 1000 \
-B 1000
```

Resulting in this tree:

```bash
(bblch:0.0212143807,(bflor:0.0198708565,(Blnc2018_re:0.0080206352,BlncHG_Trinity:0.0050874146)100/100:0.0179084105)100/100:0.0067292076,((EPI:0.0335587903,(ASY_Yue:0.0140559522,ASY:0.0162506137)100/100:0.0262314821)100/100:0.0094741904,((saccoglossus:0.3843243013,anneissia:0.5624068368)100/100:0.1282400014,((((latimeria:0.1033966438,lepisosteus:0.1340149553)100/100:0.0329687226,(callorhinchus:0.0760672984,amblyraja:0.0895193412)100/100:0.0634934855)100/100:0.1513918608,eptatretus:0.5187042999)100/100:0.2191722280,(ciona:0.5148541349,botrylloides:0.5794070871)100/100:0.5284676027)100/100:0.1085351426)100/100:0.3896930670)100/100:0.0196086266);
```
