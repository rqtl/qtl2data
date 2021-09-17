## Data on Maize MAGIC lines

This is data on a set of maize MAGIC lines. Largely 8 founders, but a
ninth founder was introduced for some lines due to incompatibility
between a pair of founders.

### Source

Dell'Acqua et al. (2015) Genetic properties of the MAGIC maize
population: a new platform for high definition QTL mapping in _Zea
mays_. Genome Biol. 16:167
[doi:10.1186/s13059-015-0716-z](https://doi.org/10.1186/s13059-015-0716-z)

The raw data are available as supplemental tables and at Figshare:

- Founders genotypes (with replicates)
  <https://doi.org/10.6084/m9.figshare.1437453>;
  direct download
  <https://ndownloader.figshare.com/files/2096751>

- MAGIC lines' genotypes
  <https://doi.org/10.6084/m9.figshare.1437449>;
  direct download
  <https://ndownloader.figshare.com/files/2096747>

- Founders sequence-based genotypes (imputed)
  <https://doi.org/10.6084/m9.figshare.1425350>;
  direct download
  <https://ndownloader.figshare.com/files/2081430>

- [Table S1](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0716-z/MediaObjects/13059_2015_716_MOESM1_ESM.xlsx):
  details of maize founder lines

- [Table S2](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0716-z/MediaObjects/13059_2015_716_MOESM2_ESM.xlsx): MM population breeding design

- [Table S3](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0716-z/MediaObjects/13059_2015_716_MOESM8_ESM.xlsx): genetic map

- [Table S6](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0716-z/MediaObjects/13059_2015_716_MOESM13_ESM.xlsx): phenotypes (pollen shed, PS; plant height PH, ear height EH, transformed grain yield GYrad)


### Files

The following files comprise the data reformatted for
[R/qtl2](https://kbroman.org/rqtl2). We're treating the data as a
recombinant inbred lines with nine founders. For the founder
genotypes, we combined the array-based data (which had a lot of
missing or heterozygous calls for the 9th founder, CML91) with the
sequence-based imputed genotypes.

- [`maize_magic.json`](maize_magic.json), the control file ([JSON format](https://json.org))
- [`maize_magic_geno.csv`](maize_magic_geno.csv), genotype data for
  the MAGIC lines
- [`maize_magic_foundergeno.csv`](maize_magic_geno.csv), genotype data
  for the 9 founder lines
- [`maize_magic_gmap.csv`](maize_magic_gmap.csv), genetic map
- [`maize_magic_pmap.csv`](maize_magic_pmap.csv), physical map
- [`maize_magic_crossinfo.csv`](maize_magic_covar.csv), cross
  information (number of generations followed by values indicating the
  relative contributions of the nine founders)
- [`maize_magic_pheno.csv`](maize_magic_pheno.csv), phenotype data
- [`maize_magic_phenocovar.csv`](maize_magic_pheno.csv), phenotype
  covariates, with just descriptions of the phenotype

The data are also available as a zip file, [`maize_magic.zip`](maize_magic.zip).


### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).


### Citation

Dell’Acqua M, et al. (2015) Genetic properties of the MAGIC maize
population: a new platform for high definition QTL mapping in _Zea
mays_. Genome Biol. 16:167
[doi:10.1186/s13059-015-0716-z](https://doi.org/10.1186/s13059-015-0716-z)


### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/main/MaizeMAGIC/maize_magic.zip")
mm <- read_cross2(file)
```
