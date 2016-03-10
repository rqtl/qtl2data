## Diversity Outcross data from Recla et al (2014)

### Source

<http://phenome.jax.org/db/q?rtn=projects/projdet&reqprojid=285>

Founder genotypes from <ftp://ftp.jax.org/MUGA/>

This is the same data as in [`DO_Recla`](../DO_Recla), but reduced to
three chromosomes, one phenotype, and with reduced markers.

### Files

- [`DOex.json`](DOex.json) - Control file in JSON format
- [`DOex_covar.csv`](DOex_covar.csv) - covariate data (individuals x
  covariates)
- [`DOex_pheno.csv`](DOex_pheno.csv) - phenotype data (individuals x
  phenotypes)
- [`DOex_geno.csv`](DOex_geno.csv) - genotype data (markers x individuals)
- [`do_foundergeno.csv`](DOex_foundergeno.csv) - founder genotype data
  (markers x founders)
- [`DOex_gmap.csv`](DOex_gmap.csv) - Genetic map of markers (positions in
  cM)
- [`DOex_pmap.csv`](DOex_pmap.csv) - Physical map of markers (positions in
  NCBI38/mm10 Mbp)

The data are also available as a zip file, [`DOex.zip`](DOex.zip).

Also included are some derived calculations:

- [`DOex_genoprobs.rds`](DOex_genoprobs) - Genotype probabilities
  calculated with `qtl2geno::calc_genoprob()`
- [`DOex_alleleprobs.rds`](DOex_alleleprobs) - Allele probabilities
  calculated from [`DOex_genoprobs.rds`](DOex_genoprobs) and collapsed
  to alleles with `qtl2geno::genoprob_to_alleleprob()`

### File format

See the [R/qtl2 input file format](http://kbroman.org/qtl2/assets/vignettes/input_files.html).


### Citations

Recla JM, Robledo RF, Gatti DM, Bult CJ, Churchill GA, Chesler EJ (2014)
[Precise genetic mapping and integrative bioinformatics in Diversity Outbred mice reveals Hydin as a novel pain gene](http://www.ncbi.nlm.nih.gov/pubmed/24700285).
Mamm Genome 25:211-222

### Use with [R/qtl2](http://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2geno)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DOex/DOex.zip")
DOex <- read_cross2(file)
```

You can load pre-calculated genotype probabilities (~19 MB) as follows:

```r
tmpfile <- tempfile()
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DOex/DOex_genoprobs.rds")
download.file(file, tmpfile)
pr <- readRDS(tmpfile)
unlink(tmpfile)
```

You can load pre-calculated allele probabilities (~5 MB) as follows:

```r
tmpfile <- tempfile()
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DOex/DOex_alleleprobs.rds")
download.file(file, tmpfile)
apr <- readRDS(tmpfile)
unlink(tmpfile)
```
