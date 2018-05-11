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














