# speeding up science metatranscriptomics binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zeyaxue/speeding-up-science-binder/master?urlpath=rstudio)

This repo contains data and codes to visualize taxa and KEGG functions as bar plots or heatmaps (with clustering). 

-   Thumbnail of expected Heatmap

<img src="https://raw.githubusercontent.com/zeyaxue/speeding-up-science-binder/master/figs/taxa_heat_thumb.png" width="50%" />

-   Thumbnail of expected bar plot

<img src="https://raw.githubusercontent.com/zeyaxue/speeding-up-science-binder/master/figs/unnamed-chunk-8-1.png" width="50%" />

-   P.S. The demonstration shown here is using data from [this paper](https://aem.asm.org/content/84/1/e02026-17.short)

Introduction
------------

The starting point of the workflow is + A count table, normalized or not. Looks like this:

    ##      CDS_ID     Low_A    Low_B   High_A   High_B
    ## 1   CDS_100 10.980509 8.364332 28.15585 0.000000
    ## 2 CDS_10011  2.284764 1.522854 10.25240 0.000000
    ## 3 CDS_10013  9.833839 6.554506  0.00000 0.000000
    ## 4 CDS_10015  7.194288 2.501831  0.00000 0.000000
    ## 5 CDS_10016  0.000000 0.000000 44.98916 0.000000
    ## 6 CDS_10017 18.738610 6.516398  0.00000 7.725487

-   A annotation or taxonomy table

<!-- -->

    ##      CDS_ID KO_ID  Domain        Phylum           Class              Order
    ## 1  CDS_2383  <NA> Archaea Euryarchaeota    Methanococci    Methanococcales
    ## 2  CDS_4184  <NA> Archaea Euryarchaeota Methanomicrobia  Methanosarcinales
    ## 3  CDS_6601  <NA> Archaea Euryarchaeota Methanomicrobia  Methanosarcinales
    ## 4  CDS_6904  <NA> Archaea Euryarchaeota    Methanococci    Methanococcales
    ## 5  CDS_7927  <NA> Archaea Euryarchaeota    Halobacteria       Natrialbales
    ## 6 CDS_15355  <NA> Archaea Euryarchaeota Methanomicrobia Methanomicrobiales
    ##                  Family              Genus                     Species
    ## 1 Methanocaldococcaceae Methanocaldococcus Methanocaldococcus_infernus
    ## 2    Methanosarcinaceae     Methanosarcina      Methanosarcina_barkeri
    ## 3    Methanosarcinaceae     Methanosarcina      Methanosarcina_barkeri
    ## 4      Methanococcaceae      Methanococcus     Methanococcus_vannielii
    ## 5         Natrialbaceae      Haloterrigena      Haloterrigena_jeotgali
    ## 6   Methanomicrobiaceae     Methanolacinia  Methanolacinia_petrolearia

-   A sample metainfo table

<!-- -->

    ##   SampleID Place Group
    ## 1   High_A  High     A
    ## 2   High_B  High     B
    ## 3    Low_A   Low     A
    ## 4    Low_B   Low     B

