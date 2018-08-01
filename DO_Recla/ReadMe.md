## Diversity Outcross data from Recla et al (2014) and Logan et al. (2013)

### Source

<https://phenome.jax.org/projects/Recla1>

Founder genotypes from <ftp://ftp.jax.org/MUGA/>

Additional phenotypes at
<https://phenome.jax.org/projects/Chesler4>
(in particular, <http://phenomedoc.jax.org/MPD_projdatasets/Chesler4.csv>)


### Files

- [`recla.json`](recla.json) - Control file in JSON format
- [`recla_covar.csv`](recla_covar.csv) - covariate data (individuals x
  covariates)
- [`recla_pheno.csv`](recla_pheno.csv) - phenotype data (individuals x
  phenotypes)
- [`recla_geno.csv`](recla_geno.csv) - genotype data (markers x individuals)
- [`do_foundergeno.csv`](recla_foundergeno.csv) - founder genotype data
  (markers x founders)
- [`recla_gmap.csv`](recla_gmap.csv) - Genetic map of markers (positions in
  cM)
- [`recla_pmap.csv`](recla_pmap.csv) - Physical map of markers (positions in
  NCBI38/mm10 Mbp)

The data are also available as a zip file, [`recla.zip`](recla.zip).

### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).


### Citations

Recla JM, Robledo RF, Gatti DM, Bult CJ, Churchill GA, Chesler EJ (2014)
[Precise genetic mapping and integrative bioinformatics in Diversity Outbred mice reveals Hydin as a novel pain gene](https://www.ncbi.nlm.nih.gov/pubmed/24700285).
Mamm Genome 25:211-222

Logan RW, Robledo RF, Recla JM, Philip VM, Bubier JA, Jay JJ, Harwood
C, Wilcox T, Gatti DM, Bult CJ, Churchill GA, Chesler EJ (2013)
[High-precision genetic mapping of behavioral traits in the diversity outbred mouse population](https://www.ncbi.nlm.nih.gov/pubmed/23433259).
Genes Brain Behav 12:424-437


### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Recla/recla.zip")
recla <- read_cross2(file)
```
