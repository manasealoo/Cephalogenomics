# Analysis of phylogenetic discordance

We noticed that some gene trees disagreed with the _Branchiostoma_-sister topology derived from our study. We therefore tested the level of phylogenetic discordance in our data.

## Distribution of quartet scores

This was tested by running ASTRAL with `-t 2` as outlined in the `Phylogenetic_reconstruction.md`, i.e.

```bash
java -jar ../../Astral/astral.5.7.7.jar -i PPP_BS_input_trees.cl.tre -o PPP_BS_astral_t2.tre -t 2 2> PPP_BS_astral_t2.log
```

resulting in:

```bash
(anneissia,(saccoglossus,(((botrylloides,ciona)'[q1=0.9781441906165563;q2=0.011499017440432325;q3=0.010356791943011544;f1=2654.683333333334;f2=31.20833333333333;f3=28.10833333333333;pp1=1.0;pp2=0.0;pp3=0.0;QC=45;EN=2714.0]':3.40147366291824,(eptatretus,((lepisosteus,latimeria)'[q1=0.6002647053824846;q2=0.28080361948140986;q3=0.1189316751361058;f1=2140.54393939394;f2=1001.3457070707075;f3=424.1103535353533;pp1=1.0;pp2=3.68172050811551E-232;pp3=0.0;QC=24;EN=3566.0]':0.5110667093105757,(amblyraja,callorhinchus)'[q1=0.9139380341311863;q2=0.046229710035680174;q3=0.039832255833133795;f1=3122.926262626264;f2=157.96691919191915;f3=136.10681818181817;pp1=1.0;pp2=0.0;pp3=0.0;QC=24;EN=3417.0]':2.044120473262542)'[q1=0.8866259182736456;q2=0.05399498229043683;q3=0.05937909943591762;f1=2730.8078282828283;f2=166.30454545454543;f3=182.88762626262627;pp1=1.0;pp2=0.0;pp3=0.0;QC=44;EN=3080.0]':1.769062324635397)'[q1=0.8864311049425553;q2=0.04294134088981418;q3=0.07062755416763053;f1=2554.6944444444443;f2=123.75694444444447;f3=203.54861111111117;pp1=1.0;pp2=0.0;pp3=0.0;QC=72;EN=2882.0]':1.7671768429741073)'[q1=0.6433822091886608;q2=0.2012040916073174;q3=0.15541369920402176;f1=2193.9333333333334;f2=686.1059523809523;f3=529.9607142857142;pp1=1.0;pp2=0.0;pp3=0.0;QC=140;EN=3410.0]':0.6250968033108345,((bblch,(bflor,(BlncHG_Trinity,Blnc2018_re)'[q1=0.9534141776934929;q2=0.020189907266261;q3=0.026395915040246947;f1=4894.828388278393;f2=103.65498390498398;f3=135.51662781662782;pp1=1.0;pp2=0.0;pp3=0.0;QC=13;EN=5134.0]':2.6570162975980915)'[q1=0.5954129072497905;q2=0.1959017362115137;q3=0.20868535653869633;f1=3075.3076659451676;f2=1011.8324675324682;f3=1077.8598665223665;pp1=1.0;pp2=0.0;pp3=0.0;QC=24;EN=5165.0]':0.4991383151300315)'[q1=0.7746303597999609;q2=0.116844022992859;q3=0.10852561720718039;f1=4004.0643298059977;f2=603.9667548500881;f3=560.9689153439155;pp1=1.0;pp2=0.0;pp3=0.0;QC=81;EN=5169.0]':1.083883665366072,((ASY_Yue,ASY)'[q1=0.9053051892400363;q2=0.050761962384084675;q3=0.043932848375883723;f1=4640.594400044426;f2=260.205819180818;f3=225.19978077477998;pp1=1.0;pp2=0.0;pp3=0.0;QC=13;EN=5126.0]':1.949768020023296,EPI)'[q1=0.4413628148096929;q2=0.3171582711012769;q3=0.2414789140890303;f1=2279.638938492064;f2=1638.1224702380953;f3=1247.2385912698414;pp1=1.0;pp2=5.5195361057352476E-58;pp3=0.0;QC=72;EN=5165.0]':0.17663702610154278)'[q1=0.9948808409131442;q2=0.0030374977706438375;q3=0.0020816613162118777;f1=4250.130952380952;f2=12.976190476190474;f3=8.892857142857142;pp1=1.0;pp2=0.0;pp3=0.0;QC=168;EN=4272.0]':4.8248219636106215)'[q1=0.6985599780947401;q2=0.18662036691310088;q3=0.11481965499215893;f1=2672.6904761904757;f2=714.009523809524;f3=439.30000000000007;pp1=1.0;pp2=0.0;pp3=0.0;QC=49;EN=3826.0]':0.7931137480077383):0.0);
```

where 
* q1,q2,q3: these three values show quartet support (as defined in the description of -t 1) for the main topology, the first alternative, and the second alternative, respectively.
* f1, f2, f3: these three values show the total number of quartet trees in all the gene trees that support the main topology, the first alternative, and the second alternative, respectively.
* pp1, pp2, pp3: these three show the local posterior probabilities (as defined in the description of -t 4) for the main topology, the first alternative, and the second alternative, respectively.
* QC: this shows the total number of quartets defined around each branch (this is what our paper calls m).
* EN: this is the effective number of genes for the branch. If you don't have any missing data, this would be the number of branches in your tree. When there are missing data, some gene trees might have nothing to say about a branch. Thus, the effective number of genes might be smaller than the total number of genes.

The logic of `*1,*2,*3` stems from `RL|SO` `RS|LO` `RO|LS` configurations of quartet topologies. `R` and `L` are from the order in the newick file. For example, at the origin of cephalochordates, the annotations show `q1=0.4413628148096929;q2=0.3171582711012769;q3=0.2414789140890303;f1=2279.638938492064;f2=1638.1224702380953;f3=1247.2385912698414;pp1=1.0;pp2=5.5195361057352476E-58;pp3=0.0;QC=72;EN=5165.0`, which means

* 44.14% (2279.64) of quartets support ((ASY,EPI),BRA) with local PP 1.0
* 31.72% (1638.12) of quartets support ((BRA,EPI),ASY) with local PP 5.5*10^-58
* 24.15% (1247.24) of quartets support ((ASY,BRA),EPI) with local PP 0
* 5165 is the effective number of gene trees.

## Internode certainty

To further assess the discordance seen in our phylogeny, we examined the [Quartet-Based Computations of Internode Certainty](https://github.com/lutteropp/QuartetScores) where score values close to 1 (or –1) suggest that the given internal branch is supported (or contested) whereas score values close to 0 indicate high levels of incongruence or the lack of phylogenetic signal (according to the [paper](https://academic.oup.com/sysbio/article/69/2/308/5556115)).

After some preprocessing,
```bash
cat /hps/nobackup/research/marioni/sodai/phylogeny/run_7_phylogenomic_dataset_3/PPP_BS_aln_out/*.treefile > catML_BS_eval_new.tre

sed 's/|XP_[0-9]*\.[0-9]\:/\:/g' catML_BS_eval_new.tre | \
sed 's/|TRINITY_GG_[0-9]*_c[0-9]*_g[0-9]*_i[0-9]*\.p[0-9]*\:/\:/g' | \
sed 's/|TRINITY_DN[0-9]*_c[0-9]*_g[0-9]*_i[0-9]*\.p[0-9]*\:/\:/g' | \
sed 's/|g[0-9]*\.t[0-9]*\:/\:/g' | \
sed 's/|Boleac\.CG\.SB_v[0-9]*\.S[0-9]*\.g[0-9]*\.[0-9]*\.t\:/\:/g' | \
sed 's/|Eptbu[0-9]*\.t[0-9]*\:/\:/g' | \
sed 's/|NP_[0-9]*\.[0-9]\:/\:/g' > catML_BS_eval_new.cl.tre
```

The concatenation tree was used as the reference tree `catML_BS_reference.tre` (since the species tree topology is the same).

we ran QuartetScore

```bash
QuartetScores \
-t 8 \
-r catML_BS_reference.tre \
-e catML_BS_eval.cl.tre \
-o catML_BS_tree
```

This resulted in the following annotation.
```bash
(bflor:0.019831,(Blnc2018_re:0.008333,BlncHG_Trinity:0.004664):0.018063[qp-ic:0.809926;lq-ic:0.784249;eqp-ic:0.809926],(bblch:0.020873,(((ASY_Yue:0.014234,ASY:0.017148):0.026322[qp-ic:0.662782;lq-ic:0.576005;eqp-ic:0.619170],EPI:0.034693):0.009572[qp-ic:0.024226;lq-ic:0.022041;eqp-ic:0.024226],((saccoglossus:0.382325,anneissia:0.550996):0.126742[qp-ic:0.266062;lq-ic:0.221047;eqp-ic:0.266062],((((latimeria:0.104214,lepisosteus:0.133424):0.031916[qp-ic:0.183695;lq-ic:0.177893;eqp-ic:0.183695],(callorhinchus:0.076915,amblyraja:0.089525):0.063518[qp-ic:0.683322;lq-ic:0.651023;eqp-ic:0.683322]):0.151402[qp-ic:0.633298;lq-ic:0.605766;eqp-ic:0.633298],eptatretus:0.531021):0.2116[qp-ic:0.670829;lq-ic:0.648686;eqp-ic:0.670829],(ciona:0.523926,botrylloides:0.576242):0.515517[qp-ic:0.896867;lq-ic:0.871096;eqp-ic:0.896867]):0.105001[qp-ic:0.193464;lq-ic:0.165745;eqp-ic:0.193464]):0.388338[qp-ic:0.984377;lq-ic:0.963000;eqp-ic:0.976753]):0.019577[qp-ic:0.400349;lq-ic:0.383173;eqp-ic:0.400349]):0.006812[qp-ic:0.141911;lq-ic:0.126829;eqp-ic:0.141911]);
```

Note: it is also interesting to see the internode certainty scores at the base of chordates.

## Analysis of the bootstrap support of the individual gene trees

For this, we used [DiscoVista](https://github.com/esayyari/DiscoVista).

```bash
mkdir PPP_BS_para_discovista

cp ../PPP_BS_aln_out/*.treefile PPP_BS_para_discovista

for f in *aln_pruned_1.fasta.treefile; do
mv -- "$f" "${f%.aln_pruned_1.fasta.treefile}_1.aln_pruned.fasta.treefile"
done
for f in *aln_pruned_2.fasta.treefile; do
mv -- "$f" "${f%.aln_pruned_2.fasta.treefile}_2.aln_pruned.fasta.treefile"
done
```

We then created a file `run_clean_geneID.sh`

```bash
#!/usr/bin/env bash
files=(*.treefile)
sed 's/|XP_[0-9]*\.[0-9]\:/\:/g' ${files[$LSB_JOBINDEX-1]} | \
sed 's/|TRINITY_GG_[0-9]*_c[0-9]*_g[0-9]*_i[0-9]*\.p[0-9]*\:/\:/g' | \
sed 's/|TRINITY_DN[0-9]*_c[0-9]*_g[0-9]*_i[0-9]*\.p[0-9]*\:/\:/g' | \
sed 's/|g[0-9]*\.t[0-9]*\:/\:/g' | \
sed 's/|Boleac\.CG\.SB_v[0-9]*\.S[0-9]*\.g[0-9]*\.[0-9]*\.t\:/\:/g' | \
sed 's/|Eptbu[0-9]*\.t[0-9]*\:/\:/g' | \
sed 's/|NP_[0-9]*\.[0-9]\:/\:/g' > ${files[$LSB_JOBINDEX-1]%.treefile}.cl.tre

chmod u+x run_clean_geneID.sh
```

which we ran

```bash
./run_clean_geneID.sh
```

These gene trees were then moved
```bash
for f in *.cl.tre; do
mkdir ${f%.aln_pruned.fasta.cl.tre}
mkdir ${f%.aln_pruned.fasta.cl.tre}/${f%.aln_pruned.fasta.cl.tre}_PPP_test-FAA
mv $f ${f%.aln_pruned.fasta.cl.tre}/${f%.aln_pruned.fasta.cl.tre}_PPP_test-FAA/estimated_gene_trees.tree
done

rm *.treefile
```

Then we ran DiscoVista
```bash
singularity exec docker://esayyari/discovista discoVista.py -m 1 -c data/parameter/clades-def_4.txt -p PPP_BS_para_discovista -t 80 -o results_BS_trees_new
```

where the clades to be tested are defined as follows.
```bash
Clade Name	Clade Definition	Section Letter	Components	Show	Comments
All	"callorhinchus""+""lepisosteus""+""latimeria""+""eptatretus""+""amblyraja""+""bblch""+""BlncHG_Trinity""+""Blnc2018_re""+""bflor""+""ASY_Yue""+""ASY""+""EPI""+""ciona""+""botrylloides""+""anneissia""+""saccoglossus"	None		0	
CHOR	"callorhinchus""+""lepisosteus""+""latimeria""+""eptatretus""+""amblyraja""+""bblch""+""BlncHG_Trinity""+""Blnc2018_re""+""bflor""+""ASY_Yue""+""ASY""+""EPI""+""ciona""+""botrylloides"	chordate		1	
VERT	"callorhinchus""+""lepisosteus""+""latimeria""+""eptatretus""+""amblyraja"	chordate		1	
CEPH	"bblch""+""BlncHG_Trinity""+""Blnc2018_re""+""bflor""+""ASY_Yue""+""ASY""+""EPI"	amphioxus		1	
TUNI	"ciona""+""botrylloides"	chordate		1	
AMBU	"anneissia""+""saccoglossus"	chordate		1	
VERT/CEPH	"callorhinchus""+""lepisosteus""+""latimeria""+""eptatretus""+""amblyraja""+""bblch""+""BlncHG_Trinity""+""Blnc2018_re""+""bflor""+""ASY_Yue""+""ASY""+""EPI"	chordate		1	
VERT/TUNI	"callorhinchus""+""lepisosteus""+""latimeria""+""eptatretus""+""amblyraja""+""ciona""+""botrylloides"	chordate		1	
BRAN	"bblch""+""BlncHG_Trinity""+""Blnc2018_re""+""bflor"	amphioxus		1	
BFLO/BLNC	"BlncHG_Trinity""+""Blnc2018_re""+""bflor"	amphioxus		1	
BBLC/BLNC	"BlncHG_Trinity""+""Blnc2018_re""+""bblch"	amphioxus		1	
ASY/EPI	"ASY_Yue""+""ASY""+""EPI"	amphioxus		1	
EPI/BRAN	"bblch""+""BlncHG_Trinity""+""Blnc2018_re""+""bflor""+""EPI"	amphioxus		1	
ASY/BRAN	"bblch""+""BlncHG_Trinity""+""Blnc2018_re""+""bflor""+""ASY""+""ASY_Yue"	amphioxus		1	
BFLO/BBLC	"bflor""+""bblch"	amphioxus		1	
CEPH/TUNI	"bblch""+""BlncHG_Trinity""+""Blnc2018_re""+""bflor""+""ASY_Yue""+""ASY""+""EPI""+""ciona""+""botrylloides"	chordate		1	
```

## Likelihood-based phylogenetic signal in the supermatrix

## Polytomy test
