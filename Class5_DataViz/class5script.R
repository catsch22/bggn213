---
  title: "BGGN Class 5"
author: Catherine Schrankel
date: April 18, 2018
output: html_document
---
  


---
  title: "BGGN Class 5"
output:
  html_document:
  code_folding: hide
---

# Bioinformatics Class 5
# Plots

x <- rnorm(1000,0)

summary(x)

# lets see this data as a graph
boxplot(x)

# Good old historgrams
hist(x)

# Section 1 from lab sheet
baby <- read.table("bggn213_05_rstats/weight_chart.txt", header = TRUE)

plot(baby, type="l", pch=15, cex=1, ylim=c(0,12), lwd=3, lty=6, xlab="age (months)")
# Now we are adding the parameters in section 1A to change the layout or parameters of the graph

# Moving to SectionB
feat <- read.table("bggn213_05_rstats/feature_counts.txt", sep="\t", header = TRUE)
# "/t" didn't work, tried \t
# then labeled this whole thing as feature by "feat <-", and changed header to TRUE
# overall we are arguing with the read function to open the data and look at it and change thigns where needed

par(mar=c(5, 8, 4, 2))
barplot(feat$Count, names.arg=feat$Feature, horiz = TRUE, las=2, ylab = "A title", main = "Some title")
# par()$mar to see the default margin values, since we want to change one
# par(mar=c(xyz)) to override the default settings, keep it above the barplot

# moving to Section2
read.table("bggn213_05_rstats/male_female_counts.txt", sep = "\t", header = TRUE)
#duplicate row.names not allowed ERROR, bc seperator fucked up, therefore change all spaces to tab
# also getting tired of typign in same commands such as sep, header, etc
# look at code of function by typing function without the ()
# use read.delim instead

mfcount <- read.delim("bggn213_05_rstats/male_female_counts.txt")
# read.csv uses seperator commas etc
# now want to plot it

barplot(mfcount$Count, col=rainbow(4))
# margins still the same bc overrode the parameters; used broom to clear it so back to default settings
# col=red to add red everywhere in all columns

rainbow(10)
# rainbow(10) is a vector, outputs colors you can use, but lists as the color ID (hexadeciaml strings)
# col =rainbow(10) makes all boxes cycle throgu 10 colors of rainbow

#mcols <- rainbow(10) will hardwire 1o colors to use, BUT to link it the size of dataset, use the nrow function
mycols <- heat.colors(nrow(mfcount))

barplot(mfcount$Count, col=mycols)
# now this will add color up to 10, the number of counts in the data set
# Changing the color function from rainbow to heat.colors or others shows different color schemes
# if you use col=c("blue2", red2"), you can vector the color on the M/F

#Section2B
expr <- read.delim("bggn213_05_rstats/up_down_expression.txt")
View(expr)
plot(expr$Condition1, expr$Condition2, col=expr$State)
nrow(expr)
#plotting RNA-Seq data! Use nrow to determine how many datapoints... is over 5000
unique(expr$State)
#tells you unique set of states. Used the $ to call State, but could have asked for other things

table(expr$State)
# very useful trick to show you how many genes are in each category/state: see that the output is 72 down, 4997 unchanging, 127 up
# by defualt, the plot output will color each stage based on the color pallete set in the parameters (par)
#pallete() will show what colors are there: jsut fills in the first three colors to the alphabetically -listed states

# Need to change the default palette such that the first three colors are blue (to correspond to down), grey (unchanging), red (up)
# double check the palette is now the scheme you want
palette(c("blue", "grey", "red"))
palette()

#replotting this indeed changed the color of plots, yay!plo
plot(expr$Condition1, expr$Condition2, col=expr$State)

#Section2C: leaving until next time (will be an example of when someone else has written an R function/script to plot data, and you want to FIX it to your liking)

