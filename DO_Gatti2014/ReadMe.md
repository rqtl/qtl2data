## Diversity Outcross data from Gatti et al (2014)

### Source

<http://phenome.jax.org/projects/Gatti2>

Founder genotypes from <ftp://ftp.jax.org/MUGA/>


### Files

- [`do.json`](do.json) - Control file in JSON format
- [`do_covar.csv`](do_covar.csv) - covariate data (individuals x
  covariates)
- [`do_pheno.csv`](do_pheno.csv) - phenotype data (individuals x
  phenotypes)
- [`do_geno.csv`](do_geno.csv) - genotype data (markers x individuals)
- [`do_foundergeno.csv`](do_foundergeno.csv) - founder genotype data
  (markers x founders)
- [`do_gmap.csv`](do_gmap.csv) - Genetic map of markers (positions in
  cM)
- [`do_pmap.csv`](do_pmap.csv) - Physical map of markers (positions in
  NCBI38/mm10 Mbp)

The data are also available as a zip file, [`do.zip`](do.zip).


### X and Y chromosome SNP intensities

Detailed SNP array data are available at the
[Mouse Phenome Database](https://phenome.jax.org/projects/Gatti2).

We've extracted the overall intensities for SNPs on the X and Y
chromosomes, see [`XYint/`](XYint/).


### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).


### Citation

Gatti DM, Svenson KL, Shabalin A, Wu L-Y, Valdar W, Simecek P, Goodwin
N, Cheng R, Pomp D, Palmer A, Chesler EJ, Broman KW, Churchill GA
(2014)
[Quantitative trait locus mapping methods for Diversity Outbred mice](https://doi.org/10.1534/g3.114.013748).
G3 (Bethesda) 4:1623-1633


### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Gatti2014/do.zip")
do <- read_cross2(file)
```
