## Diversity Outcross data from Recla et al (2019)

### Source

- Phenotype data from <https://phenome.jax.org/projects/Recla2>
- Genotype data from <https://dodb.jax.org> (Carol Bult, Genetic
  mapping in Diversity Outbred mice identifies a Trpa1 variant
influencing late...)

The genotype data were subset to the 275 mice for which phenotype data was
available, and to the markers with mapped locations in the revised
MegaMUGA annotations at <https://github.com/kbroman/MUGAarrays>.
There were 83 markers missing from the genotype data, so the available
data is the subset of 74,906 MegaMUGA markers that map to chr 1-19 or X and
are present in the data.


### Files

- [`recla2.json`](recla2.json) - Control file in JSON format
- [`recla2_covar.csv`](recla2_covar.csv) - covariate data (individuals x
  covariates)
- [`recla2_pheno.csv`](recla2_pheno.csv) - phenotype data (individuals x
  phenotypes)
- [`recla2_geno.csv`](recla2_geno.csv) - genotype data (markers x individuals)
- [`recla2_foundergeno.csv`](recla2_foundergeno.csv) - founder genotype data
  (markers x founders)
- [`recla2_gmap.csv`](recla2_gmap.csv) - Genetic map of markers (positions in
  cM)
- [`recla2_pmap.csv`](recla2_pmap.csv) - Physical map of markers (positions in
  NCBI38/mm10 Mbp)

The data are also available as a zip file, [`recla2.zip`](recla2.zip).

### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).


### Citations

Recla JM, Bubier JA, Gatti DM, Ryan JL, Long KH, Robledo RF, Glidden
NC, Hou G, Churchill GA, Maser RS, Zhang ZW, Young EE, Chesler EJ,
Bult CJ (2019) Genetic mapping in Diversity Outbred mice identifies a
_Trpa1_ variant influencing late-phase formalin response.
Pain 160:1740-1753
[doi:10.1097/j.pain.0000000000001571](https://doi.org/10.1097/j.pain.0000000000001571)

### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Recla2/recla2.zip")
recla2 <- read_cross2(file)
```
