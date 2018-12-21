## Diversity Outcross data from Karen Svenson and colleagues

### Source

Primary genotype files available at figshare, <https://doi.org/10.6084/m9.figshare.7359542.v1>


### Files

- [`svenson.zip`](do.zip) - Zip file containing
  - `svenson.json` - JSON control file
  - `svenson_covar.csv` - covariate data (individuals x covariates)
  - `svenson_geno*.csv` - genotype data (markers x individuals) for
    `*` = 1, ..., 19, X
  - `MM/MM_gmap*.csv` - Genetic map of markers (positions in cM)
  - `MM/MM_pmap*.csv` - Physical map of markers (positions in NCBI38/mm10 Mbp)
  - `MM/MM_foundergeno?.csv` - Founder genotype data (markers x individuals)

- [`svenson_chrXint.csv`](svenson_chrXint.csv) - Average intensities
  for SNPs on the X chromosome (markers x individuals)

- [`svenson_chrYint.csv`](svenson_chrYint.csv) - Average intensities
  for SNPs on the Y chromosome (markers x individuals)


### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).


### Citation

Gatti DM, Simecek P, Somes L, Jeffery CT, Vincent MJ, Choi K, Chen X,
Churchill GA, Svenson KL (2017) The effects of sex and diet on
physiology and liver gene expression in Diversity Outbred mice.
bioRxiv [doi:10.1101/098657](https://doi.org/10.1101/098657)



### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Svenson291/svenson.zip")
do <- read_cross2(file, quiet=FALSE)
```
