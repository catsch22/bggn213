---
title: "BGGN Class13_Genome"
author: "Cat"
date: "5/16/2018"
output: 
  html_document: 
    keep_md: yes
---



## Analyzing Mexican ancesttry in asthma SNPs

Downloaded dataset of asthma SNPs present in Mexican/LA patients

Read 1000 genome project data for MXL dataset:


```r
genotype <- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
```

Looking at this, second column is what we are interested in...


```r
#acess the second column
genotype[,2]
```

```
##  [1] A|A G|G A|A G|G G|G A|G A|G A|A A|G A|A G|A A|A A|A G|G A|A A|G A|G
## [18] A|G A|G G|A A|G G|G G|G G|A G|G A|G A|A A|A A|G A|A A|G G|A G|G A|A
## [35] A|A A|A G|A A|G A|G A|G A|A G|A A|G G|A G|A A|A A|A A|G A|A A|A A|G
## [52] A|G A|A G|A A|A G|A A|G A|A G|A A|G G|G A|A G|A A|G
## Levels: A|A A|G G|A G|G
```


```r
#a better way to get information
table(genotype[,2])
```

```
## 
## A|A A|G G|A G|G 
##  22  21  12   9
```


```r
#to get percentages...divided by 64 (use nrow so if data changes, can still have correct info)
table(genotype[,2]) / nrow(genotype) *100
```

```
## 
##     A|A     A|G     G|A     G|G 
## 34.3750 32.8125 18.7500 14.0625
```

RESULT: About 14% of this patient dataset have the homozygous GG SNP

##FASTQ ANALYSIS


```r
#install.packages("seqinr")
#install.packages("gtools")
```


```r
#identify the numerical values of ASCII-based score
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred
```

```
##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31
```

```r
## D D D D C D E D C D D D D B B D D D C C @
## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31
```


```r
#What is the probability 
prob <- 10**(-phred/10)
prob
```

```
##            D            D            D            D            C 
## 0.0003162278 0.0003162278 0.0003162278 0.0003162278 0.0003981072 
##            D            E            D            C            D 
## 0.0003162278 0.0002511886 0.0003162278 0.0003981072 0.0003162278 
##            D            D            D            B            B 
## 0.0003162278 0.0003162278 0.0003162278 0.0005011872 0.0005011872 
##            D            D            D            C            C 
## 0.0003162278 0.0003162278 0.0003162278 0.0003981072 0.0003981072 
##            @ 
## 0.0007943282
```

##Section 4

Is there a pattern in gene expression patterns based on genotype...?


```r
geno <- read.table("rs8067378_ENSG00000172057.6.webarchive.txt")
```


```r
#get info on this dataset... can see that 121 of the 200+ samples have the GG SNP
summary(geno)
```

```
##      sample     geno          exp        
##  HG00096:  1   A/A:108   Min.   : 6.675  
##  HG00097:  1   A/G:233   1st Qu.:20.004  
##  HG00099:  1   G/G:121   Median :25.116  
##  HG00100:  1             Mean   :25.640  
##  HG00101:  1             3rd Qu.:30.779  
##  HG00102:  1             Max.   :51.518  
##  (Other):456
```

Is the GG value (RPKM) the same or different as the AA? Statistically? 


```r
#tell me T/F value of which ones are GG
geno$geno == "G/G"
```

```
##   [1] FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE
##  [12] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE
##  [23]  TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE FALSE
##  [34] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
##  [45] FALSE  TRUE  TRUE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE
##  [56]  TRUE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
##  [67] FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE  TRUE
##  [78] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
##  [89]  TRUE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
## [100] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE
## [111]  TRUE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE FALSE
## [122] FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE
## [133] FALSE FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE
## [144] FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE
## [155] FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE
## [166]  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE FALSE
## [177] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [188] FALSE FALSE  TRUE FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE
## [199]  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
## [210] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
## [221] FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
## [232]  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE
## [243] FALSE FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE  TRUE
## [254]  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE
## [265] FALSE FALSE  TRUE  TRUE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE
## [276] FALSE FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE  TRUE FALSE
## [287]  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE
## [298] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE
## [309] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE  TRUE
## [320] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE
## [331] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
## [342] FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
## [353] FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE
## [364]  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
## [375]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE  TRUE
## [386] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE
## [397] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [408] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [419] FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
## [430] FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE
## [441] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
## [452] FALSE FALSE  TRUE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE
```



```r
#tell you the values of GG
geno$exp[geno$geno == "G/G"]
```

```
##   [1] 18.25141 17.67473 18.55622 23.10383 30.94554 21.14387 18.39547
##   [8] 12.02809 17.44761 29.82254 23.01983 13.42470 22.65437 11.07445
##  [15] 28.35841 28.79371 27.08956 16.11138 26.61928 30.18323 19.40790
##  [22] 19.52301 26.56808 17.34076 10.74263 16.66051 29.01720 20.69333
##  [29] 21.15677 18.58691 19.04962 22.81974 32.01142 21.12823 18.61268
##  [36] 19.37093 31.42162 16.67764 19.08659 21.55001  8.29591 12.58869
##  [43] 17.34109 28.23642 19.99979 25.55413 24.45672 23.53572 22.48273
##  [50] 14.66862 33.95602 18.26466 16.06661 17.32504 19.14766 12.57599
##  [57] 22.28749 17.29261 24.18141 16.07627 14.80495 23.46573 28.97074
##  [64] 27.78837 23.92355  9.55902 12.35836 22.53910 21.98118 16.40569
##  [71] 25.21931 24.32857 19.42882 26.56993 13.34557 16.60507 24.85165
##  [78] 21.56943 23.95528 16.18962 22.53720 26.04123  6.67482 20.07363
##  [85] 19.76527 18.50772 20.14146 18.07151  6.94390 22.14277 14.23742
##  [92] 19.85388 27.73467 19.02064 14.49816 26.78940 20.84709 10.77316
##  [99] 12.82128 16.90256 29.60045 14.81945 17.46326 23.26922 21.39806
## [106] 18.06320 15.91528 24.80823 26.04514 18.28089 23.24907 17.91118
## [113] 21.09502 24.74366 27.40521 24.85772 23.08482 16.56929 16.69044
## [120] 25.08880 32.78519
```


```r
#get the numeric values of SNP expression values
summary( geno$exp[geno$geno == "A/A"])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   11.40   27.02   31.25   31.82   35.92   51.52
```

```r
summary( geno$exp[geno$geno == "G/A"])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 
```

```r
summary( geno$exp[geno$geno == "G/G"])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   6.675  16.903  20.074  20.594  24.457  33.956
```


But now we want to PLOT this, therefore use the boxplot function, but need to tweak it...Formula interface vs data input

##Let's make a boxplot


```r
#add the notch to visualize that the SNP levels are indeed very different from eachother
boxplot(exp ~ geno , data =geno, notch = TRUE)
```

![](BGGN_Class13_GenomeAnalysis_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


```r
#install.packages("ggplot2")
```


```r
library(ggplot2)
```


```r
#we are telling ggplot what we have and what type of geometry we want it to look like to visualize the data
ggplot(geno, aes(geno, exp)) + geom_boxplot()
```

![](BGGN_Class13_GenomeAnalysis_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

#Switchign to Example.rmd

```r
## Histogram of the exp column with ggplot2. Note that we defined the read table as geno, not expr, which is used in the class website rmd document.
ggplot(geno, aes(exp, fill = geno)) + geom_density(alpha = 0.2)
```

![](BGGN_Class13_GenomeAnalysis_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


```r
# Boxplot with the data shown. Note that we defined the read table as geno, not expr, which is used in the class website rmd document.
ggplot(geno, aes(geno, exp, fill=geno)) + 
  geom_boxplot(notch=TRUE, outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.4)
```

![](BGGN_Class13_GenomeAnalysis_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

To do analysis exactly in example rmd, call it as follows:

```r
#url <- "https://bioboot.github.io/bggn213_S18/class-material/rs8067378_ENSG00000172057.6.txt"
#expr <- read.table(url)
```

