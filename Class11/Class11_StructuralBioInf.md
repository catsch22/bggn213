---
title: "Class11_StructuralBioInf"
author: "Cat"
date: "5/9/2018"
output: 
  html_document: 
    keep_md: yes
---



## Exploring PBD database
Downloaded CSV from “Analyze” -> “PDB Statistics” > “by Experimental Method and Molecular Type
Now as Q of what oercent structures by experimental method. (Remember to keep the intermediate file to show the output once in github online)


```r
p <- read.csv("DataExportSummary.csv", row.names=1)
#allows experimental method
p
```

```
##                     Proteins Nucleic.Acids Protein.NA.Complex Other  Total
## X-Ray                 117481          1919               6011    10 125421
## NMR                    10708          1243                249     8  12208
## Electron Microscopy     1546            31                542     0   2119
## Other                    215             4                  6    13    238
## Multi Method             116             4                  2     1    123
```


```r
sum(p$Total)
```

```
## [1] 140109
```

```r
percent <- (p$Total / sum(p$Total)) * 100
percent
```

```
## [1] 89.51673340  8.71321614  1.51239392  0.16986775  0.08778879
```


```r
names(percent) <-  row.names(p)
percent
```

```
##               X-Ray                 NMR Electron Microscopy 
##         89.51673340          8.71321614          1.51239392 
##               Other        Multi Method 
##          0.16986775          0.08778879
```

Here we can see that majority are X-ray in the database




