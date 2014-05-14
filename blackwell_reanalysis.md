
GEO data re-analysis
========================================================




Array analysis for Colin Ewald (collin.ewald@joslin.harvard.edu) at Keith Blackwell Lab: HMS Genetics and HSCI. Contact Meeta Mistry (mmistry@hsph.harvard.edu) for additional details. Request from client was:

> I am trying to show that in their published data set (that is on geo omnibus) the most overrepresented GO term are collagens. Hence, all I would like is a second opinion whether I did the re-analysis of their published counts right.

The previous study found:
> Two gene lists (897 up and 1111 down) for C.elegans. There are even more differentially expressed genes for the daf-2 vs N2 comparison, but since we compared the C.elegans genes with mouse only orthologous genes were consideredin this list.

#### The data
- 2 Illumina platforms: M. musculus and C. Elegans
- 14 samples
- 3 different count matrices
  1. C.elegans N2 var Bristol: WT vs Mutant (6 samples)
  2. M. musculus MEF: WT vs irs -/- knockout cells
  3. M. musculus MEF: insr+/- -lox vs insr +/- knockout cells

DE was evaluated using DESeq and edgeR; both use a model based on the negative binomial distribution. BH used for controlling the FDR. Transcripts with padj < 0.01 by both packages were assigned as DE. *Comparison #1 is all we are interested in.*


## Bioconductor and R libraries used

```r
library(edgeR)
library(DESeq)
library(limma)
library(gplots)
library(lattice)
library(plyr)
library(ggplot2)
library(RColorBrewer)

source("~/R/scripts/useful functions/convenience.R")
```


### Get variables
- get base directory for analyses
- specify data and results directories
- specify column headers used in metadata file


```r
baseDir = getwd()
dataDir = paste(baseDir, "/data", sep = "")
resultsDir = paste(baseDir, "/results", sep = "")
metaDir = paste(dataDir, "/meta", sep = "")
```



## Load data

```r
counts <- file.path(dataDir, "GSE36041_GSM879792-GSM879797_ce_Sample1-6_gene_counts.csv.gz")
read.counts <- data.frame(read.table(gzfile(counts), header = T, sep = ","))

# get unique gene count matrix
collapse.counts <- ddply(read.counts, .(gene), colwise(mean))
row.names(collapse.counts) <- collapse.counts[, 1]
collapse.counts <- collapse.counts[, -1]

samples <- ncol(collapse.counts)
genes <- nrow(collapse.counts)

# Get data and df for DE
counts <- na.omit(collapse.counts)
group <- factor(c(rep("WT", 5), rep("MUT", 5)))
```



## classic edgeR
-----

```r
yv <- DGEList(counts = floor(counts), group = group)
yv <- calcNormFactors(yv)

# calculate dispersion
disp <- estimateGLMCommonDisp(yv)
disp <- estimateGLMTrendedDisp(disp)
et <- exactTest(disp)

s <- summary(d <- decideTestsDGE(et))
colnames(s) <- "No. Genes"
rownames(s) <- c("down-regulated", "unchanged", "up-regulated")
kable(s, format = "html")
```

<table>
 <thead>
  <tr>
   <th align="left">   </th>
   <th align="right"> No. Genes </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> down-regulated </td>
   <td align="right">  4044 </td>
  </tr>
  <tr>
   <td align="left"> unchanged </td>
   <td align="right"> 14928 </td>
  </tr>
  <tr>
   <td align="left"> up-regulated </td>
   <td align="right">   939 </td>
  </tr>
</tbody>
</table>





## Collagen genes
A list of 174 collagen genes for C.elegans were identified using the search term 'col-' in [Wormbase](http://www.wormbase.org/), all of which are present in our count matrix. From the set of 174 genes, 74 appear to be significantly differentially expressed between wildtype and mutant. The heatmap below displays expression data from 174 collagen genes. Genes are ordered by FDR with more significant genes at the top. The large majority of significant genes appear to be down-regulated in the mutant cells.


```r
# Get collagen genes
cgene.file <- read.delim(file.path(dataDir, "col-_gene_WS241.csv"), header = T, 
    sep = ",", row.names = 1)
cgenes <- as.character(cgene.file$label[which(cgene.file$label %in% rownames(counts))])

# Set threshold
p.cutoff <- 0.001
logfc.cutoff <- 2

# heatmap of only collagen genes that are significantly DE
toptable <- topTags(et, n = nrow(counts), sort.by = "logFC")$table
toptable$threshold <- as.factor(abs(toptable$logFC) > logfc.cutoff & toptable$FDR < 
    p.cutoff)

cgene.table <- toptable[cgenes, ]
cgene.table <- cgene.table[order(cgene.table$FDR, decreasing = T), ]

# plot heatmap of all collagen genes
colors <- rev(brewer.pal(6, "YlOrRd"))
heatmap(cpm(yv$counts)[rownames(cgene.table), ], Rowv = NA, Colv = NA, scale = "row", 
    col = colors)
```

<img src="figure/collagen genes.png" title="plot of chunk collagen genes" alt="plot of chunk collagen genes" width="800px" />


## Checking housekeeping genes
Next step is to explore a couple of genes of interest. I picked thirteen housekeeping genes from a recent publication (["Selection of Reliable Reference Genes in C. elegans..."](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0031849). The pattern of expression change is not similar to that of the collagen genes, also nly one of these housekeeping genes is significant in the edgeR results (csq-1).


```r
# Subset expression data to genes of interest
hkgenes <- c("act-1", "cdc-42", "pmp-3", "eif-3.C", "actin", "act-2", "csq-1", 
    "Y45F10D.4", "tba-1", "mdh-1", "ama-1", "F35G12.2", "rbd-1")

hk <- cpm(yv$counts)[which(row.names(yv$counts) %in% hkgenes), ]
heatmap(hk, Rowv = NA, Colv = NA, scale = "row", col = colors)
```

<img src="figure/hk.png" title="plot of chunk hk" alt="plot of chunk hk" width="800px" />




