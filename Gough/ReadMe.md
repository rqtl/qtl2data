## Data from Gray et al. (2015) on a mouse intercross, Gough Island &times; WSB

### Source

This is an F2 intercross between mice from Gough Island and WSB.
Two sibling pairs from two partially-inbred lines derived from Gough
Island Mice were crossed to WSB and pairs of F1s intermated.
The phenotype is body weight measured weekly from 1-16 weeks.
The phenotype data includes smoothed body weights and estimated first
derivatives of body weight.

The phenotypes include raw body weights weekly for 16 weeks, smoothed
body weights, and the first derivative of body weight (in g/week).

Also available at the
[Mouse Phenome Database](https://phenome.jax.org/projects/Payseur1).

### Files

- [`gough.json`](gough.json), the control file ([JSON format](http://www.json.org/))
- [`gough_geno.csv`](gough_geno.csv), genotype data
- [`gough_gmap.csv`](gough_gmap.csv), genetic map
- [`gough_pmap.csv`](gough_pmap.csv), physical map
- [`gough_covar.csv`](gough_covar.csv), covariate data
- [`gough_pheno.csv`](gough_pheno.csv), phenotype data
- [`gough_phenocovar.csv`](gough_phenocovar.csv), phenotype covariates
  (e.g., the week at which the phenotypes were measured)

The data are also available as a zip file, [`gough.zip`](gough.zip).

### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).

### Citation

Gray MM, Parmenter MD, Hogan CA, Ford I, Cuthbert RJ, Ryan PG, Broman
KW, Payseur BA (2015)
[Genetics of rapid and extreme size evolution in island mice](https://www.ncbi.nlm.nih.gov/pubmed/26199233).
Genetics 201:213-228
[doi: 10.1534/genetics.115.177790](https://doi.org/10.1534/genetics.115.177790)

### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/Gough/gough.zip")
gough <- read_cross2(file)
```
