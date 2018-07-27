## Data from Alan Attie B6xBTBR intercross

### Source

Data from a B6 &times; BTBR intercross to investigate obesity and
diabetes. Taken from
[Mouse Phenome Database](https://phenome.jax.org/projects/Attie1)
(plus three clinical phenotypes considered in
[Broman et al. (2015)](https://www.ncbi.nlm.nih.gov/pubmed/26290572):
10 week insulin and binary indicators for each of agouti and tufted coats).

### Files

- [`b6btbr.json`](b6btbr.json), the control file ([JSON format](http://www.json.org/))
- [`b6btbr_geno.csv`](b6btbr_geno.csv), genotype data
- [`b6btbr_gmap.csv`](b6btbr_gmap.csv), genetic map
- [`b6btbr_pmap.csv`](b6btbr_pmap.csv), physical map
- [`b6btbr_covar.csv`](b6btbr_covar.csv), covariate data
- [`b6btbr_pheno.csv`](b6btbr_pheno.csv), phenotype data

The data are also available as a zip file, [`b6btbr.zip`](b6btbr.zip).

In addition, there are gene expression microarray data on 6 tissues;
each is a zip file containing a single csv file:

- [`b6btbr_microarray_annotation.csv.zip`](b6btbr_microarray_annotation.csv.zip),
  microarray annotation information
- [`b6btbr_adipose.csv.zip`](b6btbr_adipose.csv),
  gene expression data for adipose
- [`b6btbr_gastroc.csv.zip`](b6btbr_gastroc.csv),
  gene expression data for gastroc
- [`b6btbr_hypo.csv.zip`](b6btbr_hypo.csv),
  gene expression data for hypo
- [`b6btbr_islet.csv.zip`](b6btbr_islet.csv),
  gene expression data for islet
- [`b6btbr_kidney.csv.zip`](b6btbr_kidney.csv),
  gene expression data for kidney
- [`b6btbr_liver.csv.zip`](b6btbr_liver.csv),
  gene expression data for liver

Each of the array data files is about 50 GB; the annotion file is
about 4 GB.

### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).

### Citations

- Tian J, Keller MP, Oler AT, Rabaglia ME, Schueler KL, Stapleton DS,
  Broman AT, Zhao W, Kendziorski C, Yandell BS, Hagenbuch B, Broman
  KW, Attie AD (2015)
  [Identification of the bile acid transporter _Slco1a6_ as a candidate gene that broadly affects gene expression in mouse pancreatic islets](https://www.ncbi.nlm.nih.gov/pubmed/26385979).
  Genetics 201:1253-1262
  [doi:10.1534/genetics.115.179432](https://doi.org/10.1534/genetics.115.179432)

- Broman KW, Keller MP, Broman AT, Kendziorski C, Yandell BS, Sen &#346;,
  Attie AD (2015)
  [Identification and correction of sample mix-ups in expression genetic data: A case study](https://www.ncbi.nlm.nih.gov/pubmed/26290572).
  G3 5:2177-2186
  [doi:10.1534/g3.115.019778](https://doi.org/10.1534/g3.115.019778)

### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/B6BTBR/b6btbr.zip")
b6btbr <- read_cross2(file)
```

To read the microarray data, use `read_pheno()`. You can read the data
for a single tissue, as a matrix:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/B6BTBR/b6btbr_islet.csv.zip")
islet <- read_pheno(file)
```

Alternatively, read the data for a single tissue plus the annotation
information, creating a list with `pheno` and `phenocovar`:

```r
library(qtl2)
url <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/B6BTBR/"
islet_list <- read_pheno(paste0(url, "b6btbr_islet.csv.zip"),
                         paste0(url, "b6btbr_microarray_annotation.csv.zip"))
```
