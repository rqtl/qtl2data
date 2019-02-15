## Data from Gnan et al. (2014) 19-way Arabidopsis MAGIC

### Source

Data from [Gnan et al.
(2014)](https://doi.org/10.1534/genetics.114.170746), on seed size and
number in the 19-way Arabidopsis MAGIC lines of [Kover et al.
(2009)](https://doi.org/10.1371/journal.pgen.1000551).

The raw data were included as supplements: the genotypes in
[File S1](http://www.genetics.org/lookup/suppl/doi:10.1534/genetics.114.170746/-/DC1/genetics.114.170746-6.xls)
and the phenotypes in
[File S2](http://www.genetics.org/lookup/suppl/doi:10.1534/genetics.114.170746/-/DC1/genetics.114.170746-3.xls).

There is genotype data (at a total of 1260 markers on 5 chromosomes) on
703 MAGIC lines. 677 of the lines have phenotype data, for 8 phenotypes.

We don't have a genetic map of the markers, but we do have a physical
map, in two forms: positions in the "TAIR8 build" were present in the
Gnan et al. (2014) data (though just for 1024 markers), and positions
in the "TAIR9 build" are from <> (for 1251 markers).

For related data, see
<http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/magic.html> and
the supplementary material from [Kover et al.
(2009)](https://doi.org/10.1371/journal.pgen.1000551).


### File conversion

- [`R/convert_data.R`](R/convert_data.R) for converting the data to
  [R/qtl2](https://kbroman.org/qtl2) format.

- [`Makefile`](Makefile) to automate and document the data conversions.

- I also make use of a Python script,
  [xlsx2csv.py](https://github.com/dilshod/xlsx2csv); the `Makefile`
  includes instructions for downloading and using it.




### Files

- [`arabmagic_tair8.json`](arabmagic_tair8.json) - the JSON control
  file, pointing to the TAIR8 physical map
- [`arabmagic_tair9.json`](arabmagic_tair9.json) - the JSON control
  file, pointing to the TAIR9 physical map
- [`arabmagic_foundergeno.csv`](arabmagic_foundergeno.csv) - the founder genotype data
- [`arabmagic_geno.csv`](arabmagic_geno.csv) - the MAGIC genotype data
- [`arabmagic_pheno.csv`](arabmagic_pheno.csv) - the MAGIC phenotype data
- [`arabmagic_pmap.csv`](arabmagic_pmap.csv) - physical map of the markers
- [`arabmagic_tair8.zip`](arabmagic_tair8.zip) - all of the data
  zipped into one file, with the TAIR8 physical map
- [`arabmagic_tair9.zip`](arabmagic_tair9.zip) - all of the data
  zipped into one file, with the TAIR9 physical map




### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).




### Citations

- Gnan S, Priest A, Kover PX (2014) The genetic basis of natural
  variation in seed size and seed number and their trade-off using
  _Arabidopsis thaliana_ MAGIC lines. Genetics 198:1751-1758.
  [doi:10.1534/genetics.114.170746](https://doi.org/10.1534/genetics.114.170746)

- Kover PX, Valdar W, Trakalo J, Scarcelli N, Ehrenreich IM,
  Purugganan MD, Durrant C, Mott R (2009) A multiparent advanced
  generation inter-cross to fine-map quantitative traits in
  _Arabidopsis thaliana_. PLoS Genetics 5:e1000551.
  [doi:10.1371/journal.pgen.1000551](https://doi.org/10.1371/journal.pgen.1000551)




### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows (for the TAIR9
version of the physical map; replace 9 with 8 for the TAIR8 version):

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/ArabMAGIC/arabmagic_tair9.zip")
arab <- read_cross2(file)
```
