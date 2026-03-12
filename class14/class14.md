# Class 14 Mini Project
Karolina Navarro (PID:A19106745)

## Background

Our data for today comes from a HOX gene knock-out study

## Data Import

We have 2 key input files: counts and metadata

``` r
library(DESeq2)
```

``` r
colData <- read.csv("GSE37704_metadata.csv", row.names=1)
countData <- read.csv("GSE37704_featurecounts.csv", row.names=1)
```

## Check and tidy

``` r
head(countData)
```

                    length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ENSG00000186092    918         0         0         0         0         0
    ENSG00000279928    718         0         0         0         0         0
    ENSG00000279457   1982        23        28        29        29        28
    ENSG00000278566    939         0         0         0         0         0
    ENSG00000273547    939         0         0         0         0         0
    ENSG00000187634   3214       124       123       205       207       212
                    SRR493371
    ENSG00000186092         0
    ENSG00000279928         0
    ENSG00000279457        46
    ENSG00000278566         0
    ENSG00000273547         0
    ENSG00000187634       258

``` r
colData
```

                  condition
    SRR493366 control_sirna
    SRR493367 control_sirna
    SRR493368 control_sirna
    SRR493369      hoxa1_kd
    SRR493370      hoxa1_kd
    SRR493371      hoxa1_kd

We need to remove the first “length” column from `countData` to have a
1:1 correspondance with `colData` rows.

``` r
rownames(colData) == colnames(countData)
```

    Warning in rownames(colData) == colnames(countData): longer object length is
    not a multiple of shorter object length

    [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE

## Remove zero count genes

Some genes (rows) have no count data (i.e. zero values). We should
remove these before any further analysis

``` r
countData <- countData[, -1]   # remove length column

head(countData)
```

                    SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
    ENSG00000186092         0         0         0         0         0         0
    ENSG00000279928         0         0         0         0         0         0
    ENSG00000279457        23        28        29        29        28        46
    ENSG00000278566         0         0         0         0         0         0
    ENSG00000273547         0         0         0         0         0         0
    ENSG00000187634       124       123       205       207       212       258

``` r
dim(countData)
```

    [1] 19808     6

``` r
to.keep <- rowSums(countData) > 0
countData <- countData[to.keep, ]
```

## Setup for DESeq

``` r
dds <- DESeqDataSetFromMatrix(countData = countData, 
                       colData = colData, 
                       design = ~condition)
```

    Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    design formula are characters, converting to factors

``` r
dds <- DESeq(dds)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

## Get Results

``` r
res <- results(dds)
```

``` r
head(res)
```

    log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    Wald test p-value: condition hoxa1 kd vs control sirna 
    DataFrame with 6 rows and 6 columns
                     baseMean log2FoldChange     lfcSE       stat      pvalue
                    <numeric>      <numeric> <numeric>  <numeric>   <numeric>
    ENSG00000279457   29.9136      0.1792571 0.3248216   0.551863 5.81042e-01
    ENSG00000187634  183.2296      0.4264571 0.1402658   3.040350 2.36304e-03
    ENSG00000188976 1651.1881     -0.6927205 0.0548465 -12.630158 1.43990e-36
    ENSG00000187961  209.6379      0.7297556 0.1318599   5.534326 3.12428e-08
    ENSG00000187583   47.2551      0.0405765 0.2718928   0.149237 8.81366e-01
    ENSG00000187642   11.9798      0.5428105 0.5215598   1.040744 2.97994e-01
                           padj
                      <numeric>
    ENSG00000279457 6.86555e-01
    ENSG00000187634 5.15718e-03
    ENSG00000188976 1.76549e-35
    ENSG00000187961 1.13413e-07
    ENSG00000187583 9.19031e-01
    ENSG00000187642 4.03379e-01

## Volcano plot

``` r
library(ggplot2)

ggplot(res) +
aes(log2FoldChange,
    -log(padj)) +
  geom_point()
```

    Warning: Removed 1237 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](class14_files/figure-commonmark/unnamed-chunk-12-1.png)

Let’s add some color to this plot alomg with cutoff lines for
fold-change and P-value

``` r
mycols <- rep("gray", nrow(res))
mycols[ abs(res$log2FoldChange) > 2] <- "blue"
mycols[ res$padj > 0.01] <- "red"
head(mycols)
```

    [1] "red"  "gray" "gray" "gray" "red"  "red" 

``` r
ggplot(res) +
aes(log2FoldChange,
    -log(padj)) +
  geom_point(col = mycols) +
geom_vline(xintercept = c(-2,2)) +
  geom_hline(yintercept = -log(0.01))
```

    Warning: Removed 1237 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](class14_files/figure-commonmark/unnamed-chunk-14-1.png)

## Add Annotation

``` r
library(org.Hs.eg.db)
```

    Loading required package: AnnotationDbi

``` r
library(AnnotationDbi)

columns(org.Hs.eg.db)
```

     [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
     [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
    [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    [26] "UNIPROT"     

## MapIDs

``` r
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=rownames(res),
                      keytype = "ENSEMBL",
                     column="SYMBOL")
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=rownames(res),
                      keytype = "ENSEMBL",
                     column="ENTREZID")
```

    'select()' returned 1:many mapping between keys and columns

## Save Results

``` r
write.csv(res,file = "results_annotated.csv")
```

## PA

``` r
library(gage)
```

``` r
library(gageData)
library(pathview)
```

    ##############################################################################
    Pathview is an open source software package distributed under GNU General
    Public License version 3 (GPLv3). Details of GPLv3 is available at
    http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    formally cite the original Pathview paper (not just mention it) in publications
    or products. For details, do citation("pathview") within R.

    The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ##############################################################################

``` r
data(kegg.sets.hs)
```

``` r
foldchanges <- res$log2FoldChange
names(foldchanges) <- res$entrez
head(res$entrez)
```

    ENSG00000279457 ENSG00000187634 ENSG00000188976 ENSG00000187961 ENSG00000187583 
                 NA        "148398"         "26155"        "339451"         "84069" 
    ENSG00000187642 
            "84808" 

``` r
keggres <- gage(foldchanges, gsets=kegg.sets.hs)
```

``` r
head(keggres$less)
```

                                                      p.geomean stat.mean
    hsa04110 Cell cycle                            8.995727e-06 -4.378644
    hsa03030 DNA replication                       9.424076e-05 -3.951803
    hsa05130 Pathogenic Escherichia coli infection 1.405864e-04 -3.765330
    hsa03013 RNA transport                         1.246882e-03 -3.059466
    hsa03440 Homologous recombination              3.066756e-03 -2.852899
    hsa04114 Oocyte meiosis                        3.784520e-03 -2.698128
                                                          p.val       q.val
    hsa04110 Cell cycle                            8.995727e-06 0.001889103
    hsa03030 DNA replication                       9.424076e-05 0.009841047
    hsa05130 Pathogenic Escherichia coli infection 1.405864e-04 0.009841047
    hsa03013 RNA transport                         1.246882e-03 0.065461279
    hsa03440 Homologous recombination              3.066756e-03 0.128803765
    hsa04114 Oocyte meiosis                        3.784520e-03 0.132458191
                                                   set.size         exp1
    hsa04110 Cell cycle                                 121 8.995727e-06
    hsa03030 DNA replication                             36 9.424076e-05
    hsa05130 Pathogenic Escherichia coli infection       53 1.405864e-04
    hsa03013 RNA transport                              144 1.246882e-03
    hsa03440 Homologous recombination                    28 3.066756e-03
    hsa04114 Oocyte meiosis                             102 3.784520e-03

``` r
pathview(gene.data=foldchanges, pathway.id="hsa04110")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory C:/Users/anaka/OneDrive/Desktop/BIMM 143 R WORKSPACE/bimm143_github_redo/class14

    Info: Writing image file hsa04110.pathview.png

![](hsa04110.pathview.png)

``` r
keggrespathways <- rownames(keggres$greater)[1:5]

# Extract the 8 character long IDs part of each string
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
```

    [1] "hsa04060" "hsa05323" "hsa05146" "hsa05332" "hsa04640"

``` r
pathview(gene.data=foldchanges, pathway.id="hsa04060")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory C:/Users/anaka/OneDrive/Desktop/BIMM 143 R WORKSPACE/bimm143_github_redo/class14

    Info: Writing image file hsa04060.pathview.png

![](hsa04060.pathview.png)

``` r
pathview(gene.data=foldchanges, pathway.id="hsa05323")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory C:/Users/anaka/OneDrive/Desktop/BIMM 143 R WORKSPACE/bimm143_github_redo/class14

    Info: Writing image file hsa05323.pathview.png

![](hsa05323.pathview.png)

``` r
pathview(gene.data=foldchanges, pathway.id="hsa05146")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory C:/Users/anaka/OneDrive/Desktop/BIMM 143 R WORKSPACE/bimm143_github_redo/class14

    Info: Writing image file hsa05146.pathview.png

![](hsa05146.pathview.png)

``` r
pathview(gene.data=foldchanges, pathway.id="hsa05332")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory C:/Users/anaka/OneDrive/Desktop/BIMM 143 R WORKSPACE/bimm143_github_redo/class14

    Info: Writing image file hsa05332.pathview.png

![](hsa05332.pathview.png)

``` r
pathview(gene.data=foldchanges, pathway.id="hsa04640")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory C:/Users/anaka/OneDrive/Desktop/BIMM 143 R WORKSPACE/bimm143_github_redo/class14

    Info: Writing image file hsa04640.pathview.png

![](hsa04640.pathview.png)

``` r
data(go.sets.hs)
data(go.subs.hs)

# Focus on Biological Process subset of GO
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres = gage(foldchanges, gsets=gobpsets)

lapply(gobpres, head)
```

    $greater
                                                 p.geomean stat.mean        p.val
    GO:0007156 homophilic cell adhesion       8.519724e-05  3.824205 8.519724e-05
    GO:0002009 morphogenesis of an epithelium 1.396681e-04  3.653886 1.396681e-04
    GO:0048729 tissue morphogenesis           1.432451e-04  3.643242 1.432451e-04
    GO:0007610 behavior                       1.925222e-04  3.565432 1.925222e-04
    GO:0060562 epithelial tube morphogenesis  5.932837e-04  3.261376 5.932837e-04
    GO:0035295 tube development               5.953254e-04  3.253665 5.953254e-04
                                                  q.val set.size         exp1
    GO:0007156 homophilic cell adhesion       0.1951953      113 8.519724e-05
    GO:0002009 morphogenesis of an epithelium 0.1951953      339 1.396681e-04
    GO:0048729 tissue morphogenesis           0.1951953      424 1.432451e-04
    GO:0007610 behavior                       0.1967577      426 1.925222e-04
    GO:0060562 epithelial tube morphogenesis  0.3565320      257 5.932837e-04
    GO:0035295 tube development               0.3565320      391 5.953254e-04

    $less
                                                p.geomean stat.mean        p.val
    GO:0048285 organelle fission             1.536227e-15 -8.063910 1.536227e-15
    GO:0000280 nuclear division              4.286961e-15 -7.939217 4.286961e-15
    GO:0007067 mitosis                       4.286961e-15 -7.939217 4.286961e-15
    GO:0000087 M phase of mitotic cell cycle 1.169934e-14 -7.797496 1.169934e-14
    GO:0007059 chromosome segregation        2.028624e-11 -6.878340 2.028624e-11
    GO:0000236 mitotic prometaphase          1.729553e-10 -6.695966 1.729553e-10
                                                    q.val set.size         exp1
    GO:0048285 organelle fission             5.841698e-12      376 1.536227e-15
    GO:0000280 nuclear division              5.841698e-12      352 4.286961e-15
    GO:0007067 mitosis                       5.841698e-12      352 4.286961e-15
    GO:0000087 M phase of mitotic cell cycle 1.195672e-11      362 1.169934e-14
    GO:0007059 chromosome segregation        1.658603e-08      142 2.028624e-11
    GO:0000236 mitotic prometaphase          1.178402e-07       84 1.729553e-10

    $stats
                                              stat.mean     exp1
    GO:0007156 homophilic cell adhesion        3.824205 3.824205
    GO:0002009 morphogenesis of an epithelium  3.653886 3.653886
    GO:0048729 tissue morphogenesis            3.643242 3.643242
    GO:0007610 behavior                        3.565432 3.565432
    GO:0060562 epithelial tube morphogenesis   3.261376 3.261376
    GO:0035295 tube development                3.253665 3.253665

``` r
sig_genes <- res[res$padj <= 0.05 & !is.na(res$padj), "symbol"]
print(paste("Total number of significant genes:", length(sig_genes)))
```

    [1] "Total number of significant genes: 8147"

``` r
res_filtered <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

sig_up <- res_filtered[res_filtered$padj <= 0.05 & 
                       res_filtered$log2FoldChange > 1, "symbol"]

sig_down <- res_filtered[res_filtered$padj <= 0.05 & 
                         res_filtered$log2FoldChange < -1, "symbol"]
```

``` r
write.table(sig_genes, file="significant_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
```

## GO Online

(Optional)
