# Class 12 - Population Scale Analysis
Karolina Navarro (PID: A19106745)

\##Section 4: Population Scale Analysis

One sample is obviously not enough to know what is happening in a
population. You are interested in assessing genetic differences on a
population scale. So, you processed about ~230 samples and did the
normalization on a genome level. Now, you want to find whether there is
any association of the 4 asthma-associated SNPs (rs8067378…) on ORMDL3
expression.

How many samples do we have?

``` r
expr <- read.table("rs8067378_ENSG00000172057.6.txt",
                   header = TRUE,
                   stringsAsFactors = FALSE)
head(expr)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
nrow(expr)
```

    [1] 462

``` r
table(expr$geno)
```


    A/A A/G G/G 
    108 233 121 

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
expr %>%
  group_by(geno) %>%
  summarise(
    n = n(),
    median_expression = median(exp, na.rm = TRUE)
  )
```

    # A tibble: 3 × 3
      geno      n median_expression
      <chr> <int>             <dbl>
    1 A/A     108              31.2
    2 A/G     233              25.1
    3 G/G     121              20.1

``` r
library(ggplot2)
```

Let’s make a boxplot

``` r
ggplot(expr) + aes(geno, exp, fill=geno) +
  geom_boxplot(notch=TRUE) +
  theme_bw()
```

![](class12_files/figure-commonmark/unnamed-chunk-6-1.png)

Given the relative expression values provided, I can infer that A/A has
the highest median expression and G/G has the lowest median expression
value. This hints at the idea that there is a clear stepwise decrease in
gene expression as the G allele increases in dosage. Based on this plot,
the SNP appears to influence expression levels ORMDL3.
