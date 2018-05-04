---
title: "BioinfomaticsClass7"
output: 
  html_document: 
    keep_md: yes
---



# Functions again!

We can source any file of R code with the 'source'() function


```r
source("http://tinyurl.com/rescale-R")
```

Let's make sure things are here...can see them in environment, or use ls() to list what is in envrionment:


```r
ls()
```

```
##  [1] "both_na"         "both_na2"        "both_na3"       
##  [4] "df1"             "df2"             "df3"            
##  [7] "gene_intersect"  "gene_intersect2" "gene_intersect3"
## [10] "gene_intersect4" "rescale"         "rescale2"
```

Check our rescale function 'rescale()' is working


```r
rescale(1:10)
```

```
##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
##  [8] 0.7777778 0.8888889 1.0000000
```

Tried with character, yields an error

```r
rescale( c(1:10, "string"))
```

In order to see what's wrong:
#  the !is will provide more information of why the function when wrong
# if( !is.numeric(x) ) {
    stop("Input x should be numeric", call.=FALSE)
  }

Let's check if 'resacle2()' does any better. Indeed, now the error message is more helpful


```r
rescale2( c(1:10, "string"))
```

## EXAMPLE 2: Function for finding missing valuesin two datasets
Want to make 'both_na()' function to call out vectors whth missing info
First make simple definitions


```r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA,3,NA,3, 4)
is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```

```r
is.na(y)
```

```
## [1]  TRUE FALSE  TRUE FALSE FALSE
```

```r
!is.na(x)
```

```
## [1]  TRUE  TRUE FALSE  TRUE FALSE
```


```r
which(is.na(x))
```

```
## [1] 3 5
```

```r
#tells you which indeces of the vector have NA
which(!is.na(x))
```

```
## [1] 1 2 4
```

```r
#tells you the indeces that DO NOT have a value of NA
```


```r
z <- 3
if(z > 5) {
  print("more")
} else {print("less")
  }
```

```
## [1] "less"
```


```r
sum( is.na(x) )
```

```
## [1] 2
```

```r
#now know how MANY indeces in the vector qualify as NA
## so recap...
is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```

```r
is.na(y)
```

```
## [1]  TRUE FALSE  TRUE FALSE FALSE
```

```r
#but want to know when it occurs in BOTH.
```


```r
is.na(x) & is.na(y)
```

```
## [1] FALSE FALSE  TRUE FALSE FALSE
```

```r
#The & asks if its true in BOTH vectors, and returns it as a vector of TRUE/FALSE
#The answer TRUE corresponds to the third value in both x and y vectors
#wrap in sum to tell me how many total are TRUE for NA in both vecotrs (by basically turnign the T/F into 0/1 and coutning up the 1s)
sum( is.na(x) & is.na(y) )
```

```
## [1] 1
```

OK now we have snippet of code that we know works!

```r
sum( is.na(x) & is.na(y) )
```

```
## [1] 1
```

The function can start fmor this snippet

```r
both_na <- function(x, y) {
  sum( is.na(x) & is.na(y) )
}
```

Test it

```r
both_na(x, y)
```

```
## [1] 1
```

## "Eejit-proofing"
Now we want to test on other situations because better to catch now...


```r
x <-  c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)

both_na(x, y1)
```

```
## [1] 2
```

```r
#this works BUT if used with y2...length of y2 doesn't match lenght of x, but the error message might not be clear. Therefore want to add code that will return a more understandable error
#both_na(x, y2)
#"[1] 2 longer object length is not a multiple of shorter object length[1] 3"
```

Adding better error message

```r
function(x, y) {
  ## Check for NA elements in both input vectors and don't allow re-cycling 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  sum( is.na(x) & is.na(y) )
}
```

```
## function(x, y) {
##   ## Check for NA elements in both input vectors and don't allow re-cycling 
##   if(length(x) != length(y)) {
##     stop("Input x and y should be vectors of the same length", call.=FALSE)
##   }
##   sum( is.na(x) & is.na(y) )
## }
```


```r
#both_na2(x, y2)
#yields error bc different length values of vectors x and y2
```


```r
function(x, y) {
  ## Print some info on where NA's are as well as the number of them 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  na.in.both <- ( is.na(x) & is.na(y) )
  #asking which ones in both vectors have NA
  na.number  <- sum(na.in.both)
  #sums the amount of shared NAs
  na.which   <- which(na.in.both)
  #identify which ones

  message("Found ", na.number, " NA's at position(s):", 
          paste(na.which, collapse=", ") ) 
  
  return( list(number=na.number, which=na.which) )
  #return a list with numbers that define how many NA instanaces in each vector
}
```

```
## function(x, y) {
##   ## Print some info on where NA's are as well as the number of them 
##   if(length(x) != length(y)) {
##     stop("Input x and y should be vectors of the same length", call.=FALSE)
##   }
##   na.in.both <- ( is.na(x) & is.na(y) )
##   #asking which ones in both vectors have NA
##   na.number  <- sum(na.in.both)
##   #sums the amount of shared NAs
##   na.which   <- which(na.in.both)
##   #identify which ones
## 
##   message("Found ", na.number, " NA's at position(s):", 
##           paste(na.which, collapse=", ") ) 
##   
##   return( list(number=na.number, which=na.which) )
##   #return a list with numbers that define how many NA instanaces in each vector
## }
```


```r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

ans <- both_na3(x, y)
```

```
## Found 1 NA's at position(s):3
```


```r
ans$which
```

```
## [1] 3
```

```r
#after this, can access the ans in the environment and get the entire table of info
```

##EXAMPLE 3: Find genes that overlap btwn datasets

First check that df1 and df2 exist but looking for them in global environment (we already downloaded them, so yes)

```r
#want ID columns only and save as x
x <- df1$IDs
y <- df2$IDs
#result = vector of IDs

x
```

```
## [1] "gene1" "gene2" "gene3"
```

```r
y
```

```
## [1] "gene2" "gene4" "gene3" "gene5"
```

NOW we want to search for existing functions by ?intersect. If we 'intersect(x,y)', gives gene2 and gene3

```r
intersect(x,y)
```

```
## [1] "gene2" "gene3"
```

But I want mroe info, like the numbers for each intersecting genes etc. Intersect not enough. Use '%in%' after searching for intersection functions.

```r
x %in% y
```

```
## [1] FALSE  TRUE  TRUE
```

```r
#checking if the elements of vector x are also in vector y. Spits out a vector of logics TRUE/FALSE
y %in% x
```

```
## [1]  TRUE FALSE  TRUE FALSE
```

```r
#this will be in relation to y instead of x, therefore
#we want this position information to pull out all columns in larger datasets
```

We can use the logical output in '%in%' to get at our data. Now also use 'cbind()' function to take the inputs of two or more vector columns

```r
x %in% y
```

```
## [1] FALSE  TRUE  TRUE
```

```r
x[x %in% y]
```

```
## [1] "gene2" "gene3"
```

```r
#tells you the geneIDs within x that have values interscting with y
y[y %in% x]
```

```
## [1] "gene2" "gene3"
```

```r
cbind(x[x %in% y], y[y %in% x] )
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```

Test in other setting


```r
cbind( c("Hello", "Help"), c("Please", "Me"))
```

```
##      [,1]    [,2]    
## [1,] "Hello" "Please"
## [2,] "Help"  "Me"
```

```r
#Here we can see that the vector with hello and help in it was added with the vector of PLease and Help; in the output, vector1 is now in column1, whereas vector2 in now in column2
rbind( c("Hello", "Help"), c("Please", "Me"))
```

```
##      [,1]     [,2]  
## [1,] "Hello"  "Help"
## [2,] "Please" "Me"
```

```r
#this bound them but kept the vectors as rows
```

Now we have our working snippet, called it gene_intersect

```r
gene_intersect <- function(x, y) {
   cbind( x[ x %in% y ], y[ y %in% x ] )
}
```

Test it on x and y

```r
gene_intersect(x, y)
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```

But now we want to work with real data...change the input to take from our dataframes

```r
gene_intersect2(df1, df2)
```

```
##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
## 2 gene2   1                               -2
## 3 gene3   1                                1
```

It's working but the output is disgusting. And the code is not very clear. Hard-wired to ONLY look for columns named IDs. What if the column names were GeneID or different names, same idea?

```r
#geneintersect3 has better output, but the code is not as clear, esp if you come back to it later and forgot what you were trying to do.
gene_intersect3(df1, df2)
```

```
##     IDs exp exp2
## 2 gene2   1   -2
## 3 gene3   1    1
```

```r
function(df1, df2, gene.colname="IDs") { 

  df1.name <- df1[,gene.colname]
  df2.name <- df2[,gene.colname]

  df1.inds <- df1.name %in% df2.name
  df2.inds <- df2.name %in% df1.name

   cbind( df1[ df1.inds, ], 
          exp2=df2[ df2.inds, "exp"] )
}
```

```
## function(df1, df2, gene.colname="IDs") { 
## 
##   df1.name <- df1[,gene.colname]
##   df2.name <- df2[,gene.colname]
## 
##   df1.inds <- df1.name %in% df2.name
##   df2.inds <- df2.name %in% df1.name
## 
##    cbind( df1[ df1.inds, ], 
##           exp2=df2[ df2.inds, "exp"] )
## }
```

```r
gene_intersect4(df1, df2)
```

```
##     IDs exp exp2
## 2 gene2   1   -2
## 3 gene3   1    1
```

```r
gene_intersect4(df1, df3)
```

```
## Warning in data.frame(..., check.names = FALSE): row names were found from
## a short variable and have been discarded
```

```
##     IDs exp exp2
## 1 gene2   1   -2
## 2 gene2   1   NA
```

##MERGE FUNCTION. 

The previous results are correct according to how we wrote the function, but error message could be more explicit

```r
# Additional features we could add
# - Catch and stop when user inputs weird things
# - Use different colnames for matching in df1 and df2,
# - Match based on the content of multiple columns,
# - Optionally return rows not in df1 or not in df2 with NAs
# - Optionally sort results by matching column
# - etc...
merge(df1, df2, by="IDs")
```

```
##     IDs exp.x exp.y
## 1 gene2     1    -2
## 2 gene3     1     1
```

