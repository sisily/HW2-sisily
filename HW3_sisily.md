# HW2-sisily
Shirley Chang  
3/1/2015  

## Assignment

Reproduce the Figure 2 of the following paper: 
http://www.ncbi.nlm.nih.gov/pubmed/23220997

1. Get the data from GEO
2. Normalize the data (if necessary)
3. Use limma to test for differential expression
4. Display the results using a heatmap [Hint: Use the pheatmap package]

## Set up


```r
# Set some global knitr options
library("knitr")
opts_chunk$set(tidy=FALSE, cache=TRUE, messages=FALSE, fig.width=5, fig.height=7)
```

### Install the needed packages.


```r
packages <- c("GEOquery", "reshape", "limma", "pheatmap", "gplots")
source("http://bioconductor.org/biocLite.R")
```

```
## Bioconductor version 3.0 (BiocInstaller 1.16.1), ?biocLite for help
```

```r
for (pkg in packages)
{
    require(pkg, character.only = TRUE) || biocLite(pkg) 
}
```

```
## Loading required package: reshape
```

### Load all packages that we need:


```r
library(GEOquery)
library(reshape)
library(limma)
library(pheatmap)
```

## Process the data that we want

### Step 1: get the data


```r
# Create the data folder if it does not already exist
datadir <- "./Data/GEO/"
dir.create(file.path(datadir), showWarnings = FALSE, recursive = TRUE)

# Construct the data file path from the accession code and folder path
accession <- "GSE40812"
datafile <- paste(c(datadir, accession, "_series_matrix.txt.gz"), collapse = "")

# Download the datafile if it does not already exist and load into eset
if (file.exists(datafile)) {
    gds <- getGEO(filename = datafile) # getGEO returns an "S4" object
} else {
    gds <- getGEO(accession, destdir = datadir)[[1]]  # getGEO returns a "list"
}
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40812/matrix/
## Found 1 file(s)
## GSE40812_series_matrix.txt.gz
## File stored at: 
## ./Data/GEO//GPL10558.soft
```

```r
pd<-pData(gds)
```

### Step 2: only use "Monocyte-derived Macrophage" data


```r
mmpd <- pd[pd$source_name_ch1=="Monocyte-derived Macrophage",]
rownames(mmpd)
```

```
##  [1] "GSM1002366" "GSM1002367" "GSM1002368" "GSM1002369" "GSM1002370"
##  [6] "GSM1002371" "GSM1002372" "GSM1002373" "GSM1002374" "GSM1002375"
## [11] "GSM1002376" "GSM1002377" "GSM1002378" "GSM1002379" "GSM1002380"
## [16] "GSM1002381" "GSM1002382" "GSM1002383" "GSM1002384" "GSM1002385"
## [21] "GSM1002386" "GSM1002387" "GSM1002388" "GSM1002389" "GSM1002390"
## [26] "GSM1002391" "GSM1002392" "GSM1002393" "GSM1002394" "GSM1002395"
## [31] "GSM1002396" "GSM1002397" "GSM1002398" "GSM1002399" "GSM1002400"
## [36] "GSM1002401" "GSM1002402" "GSM1002403" "GSM1002404" "GSM1002405"
```

```r
mmeset <- gds[,rownames(mmpd)]
```

### Step 3: Sanitize the data


```r
mmpd$subj <- substring(gsub("^[A-Z][A-Z][0-9]+_","",mmpd[,'title']),1,4)
mmpd$HCV <- gsub(".*: ", "", mmpd$characteristics_ch1)
mmpd$HCV <- factor(ifelse(mmpd$HCV =="Neg", 1, 2))
mmpd$Poly_IC <- tolower(gsub(".*: ","", mmpd$characteristics_ch1.2))
mmpd$Poly_IC<- factor(ifelse(mmpd$Poly_IC == "mock",1,2))
```

### Step 4: use limma to design matrix for treatment

Get the probes that are differentially expressed between the "treatment" groups (p<0.05 and a fold change of > 1.5)


```r
mm_tx <- model.matrix(~Poly_IC,mmpd)
fit_tx <- lmFit(exprs(mmeset), mm_tx)
ebay_tx <- eBayes(fit_tx)
top_tx <- topTable(ebay_tx, coef="Poly_IC2", number=Inf)
probes_tx <- top_tx[top_tx$adj.P.Val < 0.05 & abs(top_tx$logFC)>log2(1.5), ]
nrow(probes_tx)# getting 1146 probes
```

```
## [1] 1146
```

### Step 5: subset the pData by Poly_IC value



```r
info <- mmpd[,c('geo_accession','subj','Poly_IC','HCV')]
# Poly_IC = 1 is mock
info$Poly_IC <- ifelse(info$Poly_IC == 1, -1, 1) 
exprs.probes_tx <- exprs(mmeset)[rownames(probes_tx),]
# Each subject has 2 GEO numbers: one corresponds to a mock and the other 
# a poly value (-1 is mock and 1 is poly)
info.cast <- cast(info, subj~geo_accession, value = "Poly_IC")
rownames(info.cast) <- info.cast[,'subj']
info.cast[is.na(info.cast)] <- 0
exprs.probes_tx.diff <- exprs.probes_tx%*%t(info.cast)
nrow(exprs.probes_tx.diff); ncol(exprs.probes_tx.diff)
```

```
## [1] 1146
```

```
## [1] 20
```

### Step 6: design matrix again for HCV 

Find the probes that are differentially expressed between HCV+ anf HCV- (p<0.1)

The expression data is ordered by subject ID.


```r
info.ordered <- unique(info[,c('subj','HCV')])
info.ordered <- info.ordered[order(info.ordered$subj),] 
```


```r
mm_HCV <- model.matrix(~HCV,info.ordered)
fit_HCV <- lmFit(exprs.probes_tx.diff, mm_HCV)
ebay_HCV <- eBayes(fit_HCV)
Top_HCV <- topTable(ebay_HCV, coef="HCV2", number=Inf)
```


```r
figure_probes <- Top_HCV[Top_HCV$P.Value < 0.1,]
figure_probes <- rownames(figure_probes)
length(figure_probes) # probe number =43
```

```
## [1] 43
```

## Reproduce the Figure 2 of the paper using heatmap: 

### Step 1: calculate the z value of each probes


```r
figure_probes.exprs <- exprs.probes_tx[figure_probes,]

z_scores <- function(row){
  row <- t(as.matrix(row))
  xbar <- mean(row)
  sdev <- sd(row)
  z_scores <- rep(0,ncol(row))
  for(i in 1:ncol(row)){
    z_scores[i] <- (row[i]-xbar)/sdev
  }
  return(z_scores)
}

z <- rep(0,40) 
    for(i in 1:nrow(figure_probes.exprs)){
  z <- rbind(z,z_scores(figure_probes.exprs[i,]))
    }

z <- z[-1,]
```

### Step 2: use pheatmap to create the figure 

Rename the columns, sort the columns, and create the heatmap.


```r
# Use descriptive column names instead of creating the labeled dendrogram
colnames(z) <- paste(mmpd$Poly_IC, mmpd$HCV, mmpd$subj, sep="_")
# Sort by column name (which is by treatment, infection status, and subject)
z <- z[,order(colnames(z))]
# Produce the heatmap using pheatmap, preserving column and row order
pheatmap(z, cluster_rows=F, cluster_cols=F, legend=TRUE)
```

![](HW3_sisily_files/figure-html/unnamed-chunk-13-1.png) 

Produce some alternate heatmaps for practice.


```r
# Use alternate color palatte
hmcols<-colorRampPalette(c("red", "orange", "lightyellow"))(20)
pheatmap(z, cluster_rows=F, cluster_cols=F, legend=TRUE, color=hmcols)
```

![](HW3_sisily_files/figure-html/unnamed-chunk-14-1.png) 

```r
# Remove grey border
pheatmap(z, cluster_rows=F, cluster_cols=F, legend=TRUE, color=hmcols, border_color=NA)
```

![](HW3_sisily_files/figure-html/unnamed-chunk-14-2.png) 

```r
# Use "scale" instead of custom calculated z-score scaling
colnames(figure_probes.exprs) <- paste(mmpd$Poly_IC, mmpd$HCV, mmpd$subj, sep="_")
figure_probes.exprs <- figure_probes.exprs[,order(colnames(figure_probes.exprs))]
pheatmap(figure_probes.exprs, cluster_rows=F, cluster_cols=F, scale="row", legend=TRUE, color=hmcols, show_rownames=FALSE, border_color=NA)
```

![](HW3_sisily_files/figure-html/unnamed-chunk-14-3.png) 

```r
# Set clustering_distance_rows to euclidean
pheatmap(figure_probes.exprs, cluster_rows=T, cluster_cols=F, scale="row", legend=TRUE, color=hmcols, clustering_method="ward", show_rownames=FALSE, treeheight_row=0, clustering_distance_rows="euclidean", border_color=NA)
```

```
## The "ward" method has been renamed to "ward.D"; note new "ward.D2"
```

![](HW3_sisily_files/figure-html/unnamed-chunk-14-4.png) 

```r
# Set clustering_distance_rows to euclidean manually
drows <- dist(figure_probes.exprs, method = "euclidean")
pheatmap(figure_probes.exprs, cluster_rows=T, cluster_cols=F, scale="row", legend=TRUE, color=hmcols, clustering_method="ward", show_rownames=FALSE, treeheight_row=0, clustering_distance_rows=drows, border_color=NA)
```

```
## The "ward" method has been renamed to "ward.D"; note new "ward.D2"
```

![](HW3_sisily_files/figure-html/unnamed-chunk-14-5.png) 

```r
# Using heatmap.2
library(gplots)
hclust.ward <- function(x) hclust(x, method="ward.D2")
dist.eucl <- function(x) dist(x, method="euclidean")
heatmap.2(figure_probes.exprs, scale="row", dendrogram = "none",  
          symkey=FALSE, density.info="none", trace="none", 
          cexRow=0.7, cexCol=1, margins=c(10,0), labRow=FALSE, 
          lmat=rbind( c(0, 3, 4), c(2,1,1) ), lwid=c(1, 3, 2), 
          hclustfun=hclust.ward, distfun=dist.eucl, Colv=NA, 
          key=TRUE, keysize=1.0
          )
```

![](HW3_sisily_files/figure-html/unnamed-chunk-14-6.png) 
