# speeding up science metatranscriptomics binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zeyaxue/speeding-up-science-binder/master?urlpath=rstudio)

This repo contains data and codes to visualize taxa and KEGG functions as bar plots or heatmaps (with clustering). 

* Thumbnail of expected Heatmap  
[![Binder](https://raw.githubusercontent.com/zeyaxue/speeding-up-science-binder/master/figs/taxa_heat_thumb.png)]


* Thumbnail of expected bar plot   

```{r echo=FALSE, out.width='50%'}
knitr::include_graphics("https://raw.githubusercontent.com/zeyaxue/speeding-up-science-binder/master/figs/unnamed-chunk-8-1.png")
```
 
* P.S. The demonstration shown here is using data from [this paper](https://aem.asm.org/content/84/1/e02026-17.short)

## Introduction
The starting point of the workflow is 
  + A count table, normalized or not. Looks like this:
```{r echo=FALSE}
TabTPM <- read.table(file.path("example_data/sample_TPM.tsv"),
                     header = TRUE, sep = "\t")
head(TabTPM)
```
  + A annotation or taxonomy table 
```{r echo=FALSE}
Tabanno <- read.table(file.path("example_data/sample_annotation_classifications.tsv"),
                      header = TRUE, sep = "\t", na.strings = "<NA>")
head(Tabanno)
```
  + A sample metainfo table 
```{r echo=FALSE, warning=FALSE}
samdf <- read.csv(file.path("example_data/Samdf.csv"))
head(samdf)
```
