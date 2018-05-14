---
title: "Class12_DrugDevelopment"
author: "Cat"
date: "5/11/2018"
output: 
  html_document: 
    keep_md: yes
---



## Set up HIV-Pr for docking studies

Get the library packages and protein first (download)

```r
library(bio3d)
```


```r
file.name <- get.pdb("1hsg")
```

```
## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download
```

```r
file.name
```

```
## [1] "./1hsg.pdb"
```

```r
#. means that it lives here
```

Read this file in and trim out the protein and small mlcl ligand from everything else

```r
hiv <- read.pdb(file.name)
hiv
```

```
## 
##  Call:  read.pdb(file = file.name)
## 
##    Total Models#: 1
##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
## 
##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 172  (residues: 128)
##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
## 
##    Protein sequence:
##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
##       VNIIGRNLLTQIGCTLNF
## 
## + attr: atom, xyz, seqres, helix, sheet,
##         calpha, remark, call
```

```r
#reading the file you can see there's lots of non-protein/nucleic residues...
```


Get the ligand first...

```r
ligand <- trim.pdb(hiv, "ligand")
ligand
```

```
## 
##  Call:  trim.pdb(pdb = hiv, "ligand")
## 
##    Total Models#: 1
##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
## 
##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 45  (residues: 1)
##      Non-protein/nucleic resid values: [ MK1 (1) ]
## 
## + attr: atom, helix, sheet, seqres, xyz,
##         calpha, call
```

```r
#Now we have 45 non-protein atoms...
```

Repeat for protein...

```r
protein <- trim.pdb(hiv, "protein")
protein
```

```
## 
##  Call:  trim.pdb(pdb = hiv, "protein")
## 
##    Total Models#: 1
##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
## 
##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 0  (residues: 0)
##      Non-protein/nucleic resid values: [ none ]
## 
##    Protein sequence:
##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
##       VNIIGRNLLTQIGCTLNF
## 
## + attr: atom, helix, sheet, seqres, xyz,
##         calpha, call
```

```r
#Have 1514 atoms, same as the original info sheet, but now non-protein atoms are gone. It's imprtant to check your logic liek this as you go along
```

Now we will use write function to make new files


```r
write.pdb(ligand, "1hsg_ligand.pdb")
write.pdb(protein, "1hsg_protein.pdb")
#view these and see that everything check out
```

##Now time to use AutoDock.
Installation can be a little tricky, so we do it together...that was a trip. Now to analysis...

##Analysis of AutoDock output
In order to visualize the docks and compare to the crystal conformation of the ligand we will
process the all.pdbqt to a PDB format file that can be loaded into VMD. To do this we will use R
and the Bio3D package.

Open bio3d packages


```r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
```


```r
#create an output file of the results
write.pdb(res, "results.pdb")
```

The results.pdb file is now the collected results from AutoDock. When this is visualized by the PLAY mode, it cycles through all the docking orientations of the ligand (the MK compound) in the active site/pocket of the HIV-protease.

We displayed the protein in cartoon representation, MK1 ligand in licorice (also called stick
representation) and the docked conformations (results) as licorice as well. PLayed around and it was awesome!! Can see how well the MK1 cmpd fits; it overlays almost directly the highest scored output of the AutoDock results.

>**Q6** RMSD based on non hydrogen atoms

To assess the results quantitatively we will calculate the RMSD (root mean square distance)
between each of the docking results and the known crystal structure using the bio3d package. Read the original ligand with added hydrogens that we produced earlier, and use the rmsd() function to compare to your docking results:


```r
res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("ligand.pdbqt")
rmsd(ori,res)
```

```
##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
## [11]  4.318  6.249 11.084  8.929
```

```r
#This is telling us that the first output has the best score
```


#5.2 Revist Search and retrieve Adenylate kinase structures

Below we perform a blast search of the PDB database to identify related structures to our query
Adenylate kinase UniProt sequence (with identifier P69441). In this particular example we use
the function get.seq() to fetch the query sequence directly from the UniProt sequence
database.


```r
aa <- get.seq("P69441")
```

```
## Warning in get.seq("P69441"): Removing existing file: seqs.fasta
```

```r
# Blast or hmmer search
b <- blast.pdb(aa)
```

```
##  Searching ... please wait (updates every 5 seconds) RID = FKMHCV92015 
##  
##  Reporting 90 hits
```

Function plot.blast() facilitates the visualization and filtering of the Blast results. It will
attempt to set a seed position to the point of largest drop-off in normalized scores (i.e. the
biggest jump in E-values). In this particular case we specify a cutoff (after initial plotting) of 225
to include only the relevant E.coli structures:


```r
# Plot a summary of search results to parse best hits out visually
hits <- plot(b)
```

```
##   * Possible cutoff values:    198 -3 
##             Yielding Nhits:    12 90 
## 
##   * Chosen cutoff value of:    198 
##             Yielding Nhits:    12
```

![](Class12_DockingExercises_files/figure-html/unnamed-chunk-11-1.png)<!-- -->



```r
# Print out accessions of best hits
hits$pdb.id
```

```
##  [1] "1AKE_A" "4X8M_A" "4X8H_A" "3HPR_A" "1E4V_A" "5EJE_A" "1E4Y_A"
##  [8] "3X2S_A" "4K46_A" "4NP6_A" "3GMT_A" "4PZL_A"
```


```r
# Print out info of blast results
b$hit.tbl
```

```
##         queryid subjectids identity alignmentlength mismatches gapopens
## 1  Query_239639     1AKE_A  100.000             214          0        0
## 2  Query_239639     4X8M_A   99.533             214          1        0
## 3  Query_239639     4X8H_A   99.533             214          1        0
## 4  Query_239639     3HPR_A   99.533             214          1        0
## 5  Query_239639     1E4V_A   99.533             214          1        0
## 6  Query_239639     5EJE_A   99.065             214          2        0
## 7  Query_239639     1E4Y_A   99.065             214          2        0
## 8  Query_239639     3X2S_A   98.598             214          3        0
## 9  Query_239639     4K46_A   73.239             213         57        0
## 10 Query_239639     4NP6_A   72.642             212         58        0
## 11 Query_239639     3GMT_A   62.500             216         75        1
## 12 Query_239639     4PZL_A   57.346             211         86        2
## 13 Query_239639     5G3Y_A   55.505             218         88        2
## 14 Query_239639     5G3Z_A   50.459             218         99        2
## 15 Query_239639     5G40_A   49.541             218        101        2
## 16 Query_239639     5X6J_A   50.000             218         98        3
## 17 Query_239639     2C9Y_A   53.723             188         83        1
## 18 Query_239639     1S3G_A   49.541             218         99        3
## 19 Query_239639     1AK2_A   52.660             188         85        1
## 20 Query_239639     3BE4_A   48.611             216        102        3
## 21 Query_239639     1AKY_A   46.119             219        108        3
## 22 Query_239639     3AKY_A   46.119             219        108        3
## 23 Query_239639     3FB4_A   48.165             218        104        2
## 24 Query_239639     4QBI_A   47.248             218        106        2
## 25 Query_239639     1DVR_A   45.205             219        110        3
## 26 Query_239639     3DKV_A   49.772             219         99        3
## 27 Query_239639     3DL0_A   48.165             218        104        2
## 28 Query_239639     1ZIN_A   45.413             218        110        2
## 29 Query_239639     2P3S_A   47.248             218        106        2
## 30 Query_239639     2EU8_A   47.248             218        106        2
## 31 Query_239639     1P3J_A   47.248             218        106        2
## 32 Query_239639     4QBF_A   49.772             219         99        3
## 33 Query_239639     2ORI_A   47.248             218        106        2
## 34 Query_239639     5X6I_A   46.789             218        107        2
## 35 Query_239639     2QAJ_A   47.005             217        106        2
## 36 Query_239639     2OO7_A   46.789             218        107        2
## 37 Query_239639     2OSB_A   46.789             218        107        2
## 38 Query_239639     4MKF_A   46.789             218        107        2
## 39 Query_239639     3TLX_A   44.393             214        106        3
## 40 Query_239639     4MKH_A   48.624             218        101        3
## 41 Query_239639     4QBH_A   45.872             218        109        2
## 42 Query_239639     4TYQ_A   48.165             218        102        3
## 43 Query_239639     4QBG_B   47.248             218        104        3
## 44 Query_239639     4TYP_C   47.248             218        104        3
## 45 Query_239639     4JKY_A   44.037             218        103        5
## 46 Query_239639     2RGX_A   43.578             218        104        4
## 47 Query_239639     4JLO_A   43.578             218        104        5
## 48 Query_239639     1ZAK_A   42.326             215        112        3
## 49 Query_239639     1ZD8_A   43.915             189         96        3
## 50 Query_239639     2AK3_A   44.324             185        101        2
## 51 Query_239639     4NTZ_A   38.532             218        119        4
## 52 Query_239639     2AR7_A   41.304             184        102        3
## 53 Query_239639     3NDP_A   40.761             184        103        3
## 54 Query_239639     1P4S_A   39.785             186         77        2
## 55 Query_239639     2CDN_A   39.785             186         77        2
## 56 Query_239639     3L0P_A   32.735             223        131        7
## 57 Query_239639     5X6L_A   35.784             204         98        3
## 58 Query_239639     2XB4_A   32.735             223        131        7
## 59 Query_239639     5X6K_A   35.294             204         99        3
## 60 Query_239639     3ADK_A   36.066             183         89        3
## 61 Query_239639     3UMF_A   33.333             186         92        3
## 62 Query_239639     1Z83_A   34.973             183         91        3
## 63 Query_239639     3CM0_A   34.434             212        106        5
## 64 Query_239639     1UKY_A   27.962             211        117        5
## 65 Query_239639     1TEV_A   31.963             219        109        7
## 66 Query_239639     2BWJ_A   30.851             188         98        4
## 67 Query_239639     1QF9_A   27.907             215        117        5
## 68 Query_239639     1RKB_A   40.000              35         21        0
## 69 Query_239639     3IIJ_A   40.000              35         21        0
## 70 Query_239639     5YBV_B   34.545              55         33        2
## 71 Query_239639     5YBV_A   34.545              55         33        2
## 72 Query_239639     4HBD_A   34.545              55         33        2
## 73 Query_239639     5JZV_A   40.000              35         21        0
## 74 Query_239639     5NGH_A   40.000              40         24        0
## 75 Query_239639     1Q3T_A   36.842              38         24        0
## 76 Query_239639     2GA8_A   53.846              26          8        1
## 77 Query_239639     4CVN_A   38.235              34         21        0
## 78 Query_239639     4XRU_A   52.632              19          9        0
## 79 Query_239639     1XE4_A   31.667              60         35        1
## 80 Query_239639     4II9_A   31.667              60         35        1
## 81 Query_239639     1XF8_A   43.333              30         17        0
## 82 Query_239639     1NE9_A   43.333              30         17        0
## 83 Query_239639     3GKR_A   43.333              30         17        0
## 84 Query_239639     4XRP_A   52.632              19          9        0
## 85 Query_239639     3DC0_A   38.636              44         27        0
## 86 Query_239639     5TUE_A   36.364              55         30        1
## 87 Query_239639     1J6P_A   40.541              37         21        1
## 88 Query_239639     1UA7_A   38.636              44         27        0
## 89 Query_239639     1BAG_A   38.636              44         27        0
## 90 Query_239639     1P1M_A   40.541              37         21        1
##    q.start q.end s.start s.end    evalue bitscore positives mlog.evalue
## 1        1   214       1   214 8.46e-157    432.0    100.00 359.3705104
## 2        1   214       1   214 1.51e-156    432.0    100.00 358.7911649
## 3        1   214       1   214 8.38e-156    430.0     99.53 357.0774266
## 4        1   214       1   214 1.19e-155    430.0     99.53 356.7267361
## 5        1   214       1   214 1.26e-155    430.0     99.53 356.6695777
## 6        1   214       1   214 3.77e-155    429.0     99.07 355.5736144
## 7        1   214       1   214 2.18e-154    427.0     99.07 353.8187794
## 8        1   214       1   214 3.65e-154    426.0     98.60 353.3033772
## 9        1   213       1   213 9.80e-116    329.0     84.98 264.8174884
## 10       2   213       5   216 5.44e-114    325.0     84.43 260.8009215
## 11       2   211      10   225  4.25e-90    265.0     71.30 205.7857394
## 12       2   209      26   235  1.01e-86    256.0     74.41 198.0123677
## 13       1   214       1   213  1.50e-76    230.0     68.81 174.5910020
## 14       1   214       1   213  2.88e-73    221.0     69.27 167.0309215
## 15       1   214       1   213  5.89e-71    216.0     68.35 161.7102856
## 16       1   213       1   212  1.08e-68    210.0     65.60 156.4988253
## 17       1   184      17   204  5.43e-68    209.0     69.68 154.8838472
## 18       1   213       1   212  5.58e-68    208.0     65.14 154.8565975
## 19       1   184      17   204  1.40e-67    207.0     70.21 153.9367290
## 20       2   213       7   217  3.76e-67    206.0     68.06 152.9487823
## 21       1   214       5   218  5.32e-67    206.0     65.75 152.6017279
## 22       1   214       5   218  6.47e-67    205.0     65.30 152.4060251
## 23       1   214       1   213  2.38e-66    204.0     65.14 151.1035156
## 24       1   214       1   213  1.15e-64    199.0     65.60 147.2256840
## 25       1   214       5   218  2.49e-64    199.0     64.84 146.4531632
## 26       1   214       1   213  1.34e-63    197.0     66.67 144.7701912
## 27       1   214       1   213  5.88e-63    195.0     66.97 143.2913041
## 28       1   214       1   213  6.55e-63    195.0     65.60 143.1833958
## 29       1   214       1   213  9.91e-63    194.0     66.97 142.7693165
## 30       1   214       1   213  1.04e-62    194.0     66.97 142.7210551
## 31       1   214       1   213  1.40e-62    194.0     66.97 142.4238035
## 32       1   214       1   213  2.24e-62    194.0     66.21 141.9537999
## 33       1   214       1   213  3.51e-62    193.0     66.97 141.5046597
## 34       1   214       1   213  6.67e-62    192.0     66.51 140.8626559
## 35       1   213       1   212  7.20e-62    192.0     66.82 140.7861947
## 36       1   214       1   213  8.38e-62    192.0     66.51 140.6344279
## 37       1   214       1   213  1.43e-61    192.0     66.51 140.1000162
## 38       1   214       1   213  1.59e-61    191.0     66.51 139.9939567
## 39       2   211      31   235  1.01e-60    190.0     64.95 138.1451552
## 40       1   213       3   214  1.03e-60    189.0     65.60 138.1255468
## 41       1   214       1   213  2.06e-60    189.0     64.68 137.4323996
## 42       1   213       1   212  2.53e-60    188.0     65.60 137.2268863
## 43       1   213       1   212  4.43e-59    185.0     65.60 134.3641209
## 44       1   213       1   212  8.33e-59    184.0     65.14 133.7326570
## 45       1   214       1   203  3.82e-56    177.0     66.51 127.6045148
## 46       1   214       1   203  4.38e-56    177.0     65.60 127.4677165
## 47       1   214       1   203  1.31e-55    176.0     66.51 126.3721530
## 48       1   214       6   209  2.63e-54    173.0     63.72 123.3726112
## 49       1   185       8   190  2.05e-50    164.0     64.55 114.4114149
## 50       1   185       7   189  2.81e-50    163.0     65.41 114.0960702
## 51       1   213       6   213  1.00e-46    154.0     62.39 105.9189143
## 52       1   182      28   207  3.54e-45    150.0     64.13 102.3522025
## 53       1   182       6   185  2.96e-44    148.0     63.59 100.2285548
## 54       1   182       1   155  4.92e-39    133.0     56.99  88.2075101
## 55       1   182      21   175  8.85e-39    133.0     56.99  87.6204012
## 56       1   209       1   218  1.89e-31    115.0     54.26  70.7435611
## 57       3   205      13   184  2.31e-31    114.0     52.45  70.5428904
## 58       1   209       1   218  2.50e-31    114.0     54.26  70.4638472
## 59       3   205      13   184  3.86e-31    113.0     52.45  70.0294707
## 60       3   184      12   167  1.59e-30    111.0     54.10  68.6138188
## 61       3   185      32   188  8.40e-30    110.0     56.45  66.9493211
## 62       3   184      12   167  9.85e-30    109.0     54.64  66.7900813
## 63       3   214       7   185  4.75e-29    107.0     50.94  65.2168231
## 64       3   210      18   196  5.14e-25     97.8     53.55  55.9275742
## 65       3   213       6   192  3.64e-24     95.5     49.77  53.9700586
## 66       3   187      15   173  3.76e-22     90.1     49.47  49.3324531
## 67       3   213       9   189  1.46e-20     85.9     47.44  45.6732654
## 68       2    36       6    40  3.90e-01     31.6     57.14   0.9416085
## 69       2    36      13    47  4.80e-01     31.6     57.14   0.7339692
## 70     153   204      59   113  6.20e-01     31.6     49.09   0.4780358
## 71     153   204      56   110  6.60e-01     31.2     49.09   0.4155154
## 72     153   204      74   128  7.20e-01     31.2     49.09   0.3285041
## 73       2    36      39    73  7.30e-01     31.2     57.14   0.3147107
## 74     134   173     102   141  7.80e-01     30.8     55.00   0.2484614
## 75       1    38      17    54  2.60e+00     29.6     57.89  -0.9555114
## 76       3    24      27    52  4.30e+00     28.9     69.23  -1.4586150
## 77       1    34      12    45  5.60e+00     28.5     52.94  -1.7227666
## 78       3    21      14    32  6.00e+00     28.5     84.21  -1.7917595
## 79     110   169      80   133  6.30e+00     28.5     45.00  -1.8405496
## 80     110   169      81   134  6.40e+00     28.5     45.00  -1.8562980
## 81     110   139      80   109  6.40e+00     28.5     63.33  -1.8562980
## 82     110   139      80   109  6.50e+00     28.5     63.33  -1.8718022
## 83     110   139      81   110  6.50e+00     28.5     63.33  -1.8718022
## 84       3    21      14    32  6.70e+00     28.5     84.21  -1.9021075
## 85      27    70       1    44  6.80e+00     28.5     45.45  -1.9169226
## 86     116   170      82   131  6.80e+00     28.5     43.64  -1.9169226
## 87      98   133     345   381  7.00e+00     28.5     54.05  -1.9459101
## 88      27    70       1    44  7.10e+00     28.5     45.45  -1.9600948
## 89      27    70       4    47  7.50e+00     28.5     45.45  -2.0149030
## 90      98   133     333   369  8.70e+00     28.1     54.05  -2.1633230
##    pdb.id    acc
## 1  1AKE_A 1AKE_A
## 2  4X8M_A 4X8M_A
## 3  4X8H_A 4X8H_A
## 4  3HPR_A 3HPR_A
## 5  1E4V_A 1E4V_A
## 6  5EJE_A 5EJE_A
## 7  1E4Y_A 1E4Y_A
## 8  3X2S_A 3X2S_A
## 9  4K46_A 4K46_A
## 10 4NP6_A 4NP6_A
## 11 3GMT_A 3GMT_A
## 12 4PZL_A 4PZL_A
## 13 5G3Y_A 5G3Y_A
## 14 5G3Z_A 5G3Z_A
## 15 5G40_A 5G40_A
## 16 5X6J_A 5X6J_A
## 17 2C9Y_A 2C9Y_A
## 18 1S3G_A 1S3G_A
## 19 1AK2_A 1AK2_A
## 20 3BE4_A 3BE4_A
## 21 1AKY_A 1AKY_A
## 22 3AKY_A 3AKY_A
## 23 3FB4_A 3FB4_A
## 24 4QBI_A 4QBI_A
## 25 1DVR_A 1DVR_A
## 26 3DKV_A 3DKV_A
## 27 3DL0_A 3DL0_A
## 28 1ZIN_A 1ZIN_A
## 29 2P3S_A 2P3S_A
## 30 2EU8_A 2EU8_A
## 31 1P3J_A 1P3J_A
## 32 4QBF_A 4QBF_A
## 33 2ORI_A 2ORI_A
## 34 5X6I_A 5X6I_A
## 35 2QAJ_A 2QAJ_A
## 36 2OO7_A 2OO7_A
## 37 2OSB_A 2OSB_A
## 38 4MKF_A 4MKF_A
## 39 3TLX_A 3TLX_A
## 40 4MKH_A 4MKH_A
## 41 4QBH_A 4QBH_A
## 42 4TYQ_A 4TYQ_A
## 43 4QBG_B 4QBG_B
## 44 4TYP_C 4TYP_C
## 45 4JKY_A 4JKY_A
## 46 2RGX_A 2RGX_A
## 47 4JLO_A 4JLO_A
## 48 1ZAK_A 1ZAK_A
## 49 1ZD8_A 1ZD8_A
## 50 2AK3_A 2AK3_A
## 51 4NTZ_A 4NTZ_A
## 52 2AR7_A 2AR7_A
## 53 3NDP_A 3NDP_A
## 54 1P4S_A 1P4S_A
## 55 2CDN_A 2CDN_A
## 56 3L0P_A 3L0P_A
## 57 5X6L_A 5X6L_A
## 58 2XB4_A 2XB4_A
## 59 5X6K_A 5X6K_A
## 60 3ADK_A 3ADK_A
## 61 3UMF_A 3UMF_A
## 62 1Z83_A 1Z83_A
## 63 3CM0_A 3CM0_A
## 64 1UKY_A 1UKY_A
## 65 1TEV_A 1TEV_A
## 66 2BWJ_A 2BWJ_A
## 67 1QF9_A 1QF9_A
## 68 1RKB_A 1RKB_A
## 69 3IIJ_A 3IIJ_A
## 70 5YBV_B 5YBV_B
## 71 5YBV_A 5YBV_A
## 72 4HBD_A 4HBD_A
## 73 5JZV_A 5JZV_A
## 74 5NGH_A 5NGH_A
## 75 1Q3T_A 1Q3T_A
## 76 2GA8_A 2GA8_A
## 77 4CVN_A 4CVN_A
## 78 4XRU_A 4XRU_A
## 79 1XE4_A 1XE4_A
## 80 4II9_A 4II9_A
## 81 1XF8_A 1XF8_A
## 82 1NE9_A 1NE9_A
## 83 3GKR_A 3GKR_A
## 84 4XRP_A 4XRP_A
## 85 3DC0_A 3DC0_A
## 86 5TUE_A 5TUE_A
## 87 1J6P_A 1J6P_A
## 88 1UA7_A 1UA7_A
## 89 1BAG_A 1BAG_A
## 90 1P1M_A 1P1M_A
```


```r
# Reading and writing FASTA files ourselves...
write.fasta(aa, file="ruby.fa")
aa2 <- read.fasta("ruby.fa")
aa2
```

```
##                                1        .         .         .         .         50 
## [Truncated_Name:1]sp|P69441.   MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAK
##                                1        .         .         .         .         50 
## 
##                               51        .         .         .         .         100 
## [Truncated_Name:1]sp|P69441.   DIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAG
##                               51        .         .         .         .         100 
## 
##                              101        .         .         .         .         150 
## [Truncated_Name:1]sp|P69441.   INVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTG
##                              101        .         .         .         .         150 
## 
##                              151        .         .         .         .         200 
## [Truncated_Name:1]sp|P69441.   EELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTK
##                              151        .         .         .         .         200 
## 
##                              201        .   214 
## [Truncated_Name:1]sp|P69441.   PVAEVRADLEKILG
##                              201        .   214 
## 
## Call:
##   read.fasta(file = "ruby.fa")
## 
## Class:
##   fasta
## 
## Alignment dimensions:
##   1 sequence rows; 214 position columns (214 non-gap, 0 gap) 
## 
## + attr: id, ali, call
```


```r
# Fetch PDBs: SPlit by chain, compress to make faster
files <- get.pdb(hits$pdb.id, path="pdbs", split=TRUE, gzip=TRUE)
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/1AKE.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/4X8M.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/4X8H.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/3HPR.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/1E4V.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/5EJE.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/1E4Y.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/3X2S.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/4K46.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/4NP6.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/3GMT.pdb.gz exists. Skipping download
```

```
## Warning in get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE):
## pdbs/4PZL.pdb.gz exists. Skipping download
```

```
##   |                                                                         |                                                                 |   0%  |                                                                         |=====                                                            |   8%  |                                                                         |===========                                                      |  17%  |                                                                         |================                                                 |  25%  |                                                                         |======================                                           |  33%  |                                                                         |===========================                                      |  42%  |                                                                         |================================                                 |  50%  |                                                                         |======================================                           |  58%  |                                                                         |===========================================                      |  67%  |                                                                         |=================================================                |  75%  |                                                                         |======================================================           |  83%  |                                                                         |============================================================     |  92%  |                                                                         |=================================================================| 100%
```


Next we use the pdbaln() function to align fit all the PDB structures.


```r
# Align & superimpose structures
pdbs <- pdbaln(files, fit=TRUE)
```

```
## Reading PDB files:
## pdbs/split_chain/1AKE_A.pdb
## pdbs/split_chain/4X8M_A.pdb
## pdbs/split_chain/4X8H_A.pdb
## pdbs/split_chain/3HPR_A.pdb
## pdbs/split_chain/1E4V_A.pdb
## pdbs/split_chain/5EJE_A.pdb
## pdbs/split_chain/1E4Y_A.pdb
## pdbs/split_chain/3X2S_A.pdb
## pdbs/split_chain/4K46_A.pdb
## pdbs/split_chain/4NP6_A.pdb
## pdbs/split_chain/3GMT_A.pdb
## pdbs/split_chain/4PZL_A.pdb
##    PDB has ALT records, taking A only, rm.alt=TRUE
## ...   PDB has ALT records, taking A only, rm.alt=TRUE
## ..   PDB has ALT records, taking A only, rm.alt=TRUE
## ...   PDB has ALT records, taking A only, rm.alt=TRUE
## ....
## 
## Extracting sequences
## 
## pdb/seq: 1   name: pdbs/split_chain/1AKE_A.pdb 
##    PDB has ALT records, taking A only, rm.alt=TRUE
## pdb/seq: 2   name: pdbs/split_chain/4X8M_A.pdb 
## pdb/seq: 3   name: pdbs/split_chain/4X8H_A.pdb 
## pdb/seq: 4   name: pdbs/split_chain/3HPR_A.pdb 
##    PDB has ALT records, taking A only, rm.alt=TRUE
## pdb/seq: 5   name: pdbs/split_chain/1E4V_A.pdb 
## pdb/seq: 6   name: pdbs/split_chain/5EJE_A.pdb 
##    PDB has ALT records, taking A only, rm.alt=TRUE
## pdb/seq: 7   name: pdbs/split_chain/1E4Y_A.pdb 
## pdb/seq: 8   name: pdbs/split_chain/3X2S_A.pdb 
## pdb/seq: 9   name: pdbs/split_chain/4K46_A.pdb 
##    PDB has ALT records, taking A only, rm.alt=TRUE
## pdb/seq: 10   name: pdbs/split_chain/4NP6_A.pdb 
## pdb/seq: 11   name: pdbs/split_chain/3GMT_A.pdb 
## pdb/seq: 12   name: pdbs/split_chain/4PZL_A.pdb
```


```r
pdbs
```

```
##                                 1        .         .         .         40 
## [Truncated_Name:1]1AKE_A.pdb    ----------MRIILLGAPGAGKGTQAQFIMEKYGIPQIS
## [Truncated_Name:2]4X8M_A.pdb    ----------MRIILLGAPGAGKGTQAQFIMEKYGIPQIS
## [Truncated_Name:3]4X8H_A.pdb    ----------MRIILLGAPGAGKGTQAQFIMEKYGIPQIS
## [Truncated_Name:4]3HPR_A.pdb    ----------MRIILLGAPGAGKGTQAQFIMEKYGIPQIS
## [Truncated_Name:5]1E4V_A.pdb    ----------MRIILLGAPVAGKGTQAQFIMEKYGIPQIS
## [Truncated_Name:6]5EJE_A.pdb    ----------MRIILLGAPGAGKGTQAQFIMEKYGIPQIS
## [Truncated_Name:7]1E4Y_A.pdb    ----------MRIILLGALVAGKGTQAQFIMEKYGIPQIS
## [Truncated_Name:8]3X2S_A.pdb    ----------MRIILLGAPGAGKGTQAQFIMEKYGIPQIS
## [Truncated_Name:9]4K46_A.pdb    ----------MRIILLGAPGAGKGTQAQFIMAKFGIPQIS
## [Truncated_Name:10]4NP6_A.pdb   --------NAMRIILLGAPGAGKGTQAQFIMEKFGIPQIS
## [Truncated_Name:11]3GMT_A.pdb   ----------MRLILLGAPGAGKGTQANFIKEKFGIPQIS
## [Truncated_Name:12]4PZL_A.pdb   TENLYFQSNAMRIILLGAPGAGKGTQAKIIEQKYNIAHIS
##                                           **^*****  *******  *  *^ *  ** 
##                                 1        .         .         .         40 
## 
##                                41        .         .         .         80 
## [Truncated_Name:1]1AKE_A.pdb    TGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKE
## [Truncated_Name:2]4X8M_A.pdb    TGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKE
## [Truncated_Name:3]4X8H_A.pdb    TGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKE
## [Truncated_Name:4]3HPR_A.pdb    TGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKE
## [Truncated_Name:5]1E4V_A.pdb    TGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKE
## [Truncated_Name:6]5EJE_A.pdb    TGDMLRAAVKSGSELGKQAKDIMDACKLVTDELVIALVKE
## [Truncated_Name:7]1E4Y_A.pdb    TGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKE
## [Truncated_Name:8]3X2S_A.pdb    TGDMLRAAVKSGSELGKQAKDIMDCGKLVTDELVIALVKE
## [Truncated_Name:9]4K46_A.pdb    TGDMLRAAIKAGTELGKQAKSVIDAGQLVSDDIILGLVKE
## [Truncated_Name:10]4NP6_A.pdb   TGDMLRAAIKAGTELGKQAKAVIDAGQLVSDDIILGLIKE
## [Truncated_Name:11]3GMT_A.pdb   TGDMLRAAVKAGTPLGVEAKTYMDEGKLVPDSLIIGLVKE
## [Truncated_Name:12]4PZL_A.pdb   TGDMIRETIKSGSALGQELKKVLDAGELVSDEFIIKIVKD
##                                 ****^*  ^* *^ **   *  ^*   ** *  ^^ ^^*^ 
##                                41        .         .         .         80 
## 
##                                81        .         .         .         120 
## [Truncated_Name:1]1AKE_A.pdb    RIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFD
## [Truncated_Name:2]4X8M_A.pdb    RIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFD
## [Truncated_Name:3]4X8H_A.pdb    RIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFD
## [Truncated_Name:4]3HPR_A.pdb    RIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFD
## [Truncated_Name:5]1E4V_A.pdb    RIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFD
## [Truncated_Name:6]5EJE_A.pdb    RIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFD
## [Truncated_Name:7]1E4Y_A.pdb    RIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFD
## [Truncated_Name:8]3X2S_A.pdb    RIAQEDSRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFD
## [Truncated_Name:9]4K46_A.pdb    RIAQDDCAKGFLLDGFPRTIPQADGLKEVGVVVDYVIEFD
## [Truncated_Name:10]4NP6_A.pdb   RIAQADCEKGFLLDGFPRTIPQADGLKEMGINVDYVIEFD
## [Truncated_Name:11]3GMT_A.pdb   RLKEADCANGYLFDGFPRTIAQADAMKEAGVAIDYVLEID
## [Truncated_Name:12]4PZL_A.pdb   RISKNDCNNGFLLDGVPRTIPQAQELDKLGVNIDYIVEVD
##                                 *^   *   *^* ** **** **  ^   *^ ^**^^* * 
##                                81        .         .         .         120 
## 
##                               121        .         .         .         160 
## [Truncated_Name:1]1AKE_A.pdb    VPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTG
## [Truncated_Name:2]4X8M_A.pdb    VPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTG
## [Truncated_Name:3]4X8H_A.pdb    VPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTG
## [Truncated_Name:4]3HPR_A.pdb    VPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDGTG
## [Truncated_Name:5]1E4V_A.pdb    VPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTG
## [Truncated_Name:6]5EJE_A.pdb    VPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTG
## [Truncated_Name:7]1E4Y_A.pdb    VPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTG
## [Truncated_Name:8]3X2S_A.pdb    VPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTG
## [Truncated_Name:9]4K46_A.pdb    VADSVIVERMAGRRAHLASGRTYHNVYNPPKVEGKDDVTG
## [Truncated_Name:10]4NP6_A.pdb   VADDVIVERMAGRRAHLPSGRTYHVVYNPPKVEGKDDVTG
## [Truncated_Name:11]3GMT_A.pdb   VPFSEIIERMSGRRTHPASGRTYHVKFNPPKVEGKDDVTG
## [Truncated_Name:12]4PZL_A.pdb   VADNLLIERITGRRIHPASGRTYHTKFNPPKVADKDDVTG
##                                 *    ^^^*^ *** *  *** **  ^*****  *** ** 
##                               121        .         .         .         160 
## 
##                               161        .         .         .         200 
## [Truncated_Name:1]1AKE_A.pdb    EELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGN
## [Truncated_Name:2]4X8M_A.pdb    EELTTRKDDQEETVRKRLVEWHQMTAPLIGYYSKEAEAGN
## [Truncated_Name:3]4X8H_A.pdb    EELTTRKDDQEETVRKRLVEYHQMTAALIGYYSKEAEAGN
## [Truncated_Name:4]3HPR_A.pdb    EELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGN
## [Truncated_Name:5]1E4V_A.pdb    EELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGN
## [Truncated_Name:6]5EJE_A.pdb    EELTTRKDDQEECVRKRLVEYHQMTAPLIGYYSKEAEAGN
## [Truncated_Name:7]1E4Y_A.pdb    EELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGN
## [Truncated_Name:8]3X2S_A.pdb    EELTTRKDDQEETVRKRLCEYHQMTAPLIGYYSKEAEAGN
## [Truncated_Name:9]4K46_A.pdb    EDLVIREDDKEETVLARLGVYHNQTAPLIAYYGKEAEAGN
## [Truncated_Name:10]4NP6_A.pdb   EDLVIREDDKEETVRARLNVYHTQTAPLIEYYGKEAAAGK
## [Truncated_Name:11]3GMT_A.pdb   EPLVQRDDDKEETVKKRLDVYEAQTKPLITYYGDWARRGA
## [Truncated_Name:12]4PZL_A.pdb   EPLITRTDDNEDTVKQRLSVYHAQTAKLIDFYRNFSSTNT
##                                 * *  * ** *^ *  **  ^   *  ** ^*         
##                               161        .         .         .         200 
## 
##                               201        .         .      227 
## [Truncated_Name:1]1AKE_A.pdb    T--KYAKVDGTKPVAEVRADLEKILG-
## [Truncated_Name:2]4X8M_A.pdb    T--KYAKVDGTKPVAEVRADLEKILG-
## [Truncated_Name:3]4X8H_A.pdb    T--KYAKVDGTKPVAEVRADLEKILG-
## [Truncated_Name:4]3HPR_A.pdb    T--KYAKVDGTKPVAEVRADLEKILG-
## [Truncated_Name:5]1E4V_A.pdb    T--KYAKVDGTKPVAEVRADLEKILG-
## [Truncated_Name:6]5EJE_A.pdb    T--KYAKVDGTKPVAEVRADLEKILG-
## [Truncated_Name:7]1E4Y_A.pdb    T--KYAKVDGTKPVAEVRADLEKILG-
## [Truncated_Name:8]3X2S_A.pdb    T--KYAKVDGTKPVAEVRADLEKILG-
## [Truncated_Name:9]4K46_A.pdb    T--QYLKFDGTKAVAEVSAELEKALA-
## [Truncated_Name:10]4NP6_A.pdb   T--QYLKFDGTKQVSEVSADIAKALA-
## [Truncated_Name:11]3GMT_A.pdb   E-------NGLKAPA-----YRKISG-
## [Truncated_Name:12]4PZL_A.pdb   KIPKYIKINGDQAVEKVSQDIFDQLNK
##                                          *                  
##                               201        .         .      227 
## 
## Call:
##   pdbaln(files = files, fit = TRUE)
## 
## Class:
##   pdbs, fasta
## 
## Alignment dimensions:
##   12 sequence rows; 227 position columns (204 non-gap, 23 gap) 
## 
## + attr: xyz, resno, b, chain, id, ali, resid, sse, call
```

```r
#pritn out sequence alignment independt of BLAST website
```


```r
#view(pdbs)
#in developer mode, you can see the video in the markdown
```

PCS to see patterns

```r
# Perform PCA & plot the results
pc.xray <- pca(pdbs)
plot(pc.xray)
```

![](Class12_DockingExercises_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

```r
#tons of variance captures. Scree plot looks great.
```


```r
# Visualize first principal component
pc1 <- mktrj(pc.xray, pc=1, file="pc_1.pdb")
```

In VMD, chose the “Drawing Method” Tube and “Coloring Method” Index. Then clicked the play button shown below to animate the structure and **visualize the major structural variations along PC1**.


```r
#view(pc1)
#in developer mode, you can see the video in the markdown
```


```r
#To get the interactive video into the Markdown html output, use developer version of bio3d
#plot.pc <- view(pcl)
#rgl::rglwidget(elementId = "plot.pc")
```

