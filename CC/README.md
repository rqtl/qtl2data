## Collaborative Cross mice

These are genotypes and related information for 69 mouse
Collaborative Cross lines, derived from the supplemental data for
Srivastava et al. (2017) Genomes of the Mouse Collaborative Cross. Genetics 206:537-556,
[doi:10.1534/genetics.116.198838](https://doi.org/10.1534/genetics.116.198838).

The supplementary data is at Zenodo,
[doi:10.5281/zenodo.377036](https://doi.org/10.5281/zenodo.377036).

**Warning**: There are some problems here, for example the published
data suggest that CC037 has Y chromosome and a significant portion of
the X chromosome from the parental strain NOD (D), but the way R/qtl2
currently handles these 8-way RIL data, this is not possible.
There are two other strains with similar problems on the X chromosome:
CC031 and CC056.

We've used the following:

- `genotypes.zip`

  - contains `MRCAgenotypes.csv`, `SEQgenotypes.csv`, and
    `CC018genotypes.csv`; the latter maybe should have included the
    GigaMUGA genotypes? We'll use just `SEQgenotypes.csv`, and will
    extract the markers in the GigaMUGA and/or MegaMUGA data.

- `Prob36.zip`

  - contains a `.csv.gz` file for each strain which has marker,
    chromosome, build38 bp position, and then the 36 genotypes

- `SupplementalData.zip`

  - contains `CCStrains.csv` with columns `Strain`, `N_Founders`, `ChrY`, `Mitochondria`


The founder genotypes are from FigShare
([doi:10.6084/m9.figshare.5404762.v2](https://doi.org/10.6084/m9.figshare.5404762.v2))
and combine markers from the MegaMUGA and the GigaMUGA arrays.

The script [`R/convert_cc_data.R`](R/convert_cc_data.R) both downloads
and converts the data into a form useful for R/qtl2.



### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).



### Citation

Srivastava et al. (2017) Genomes of the Mouse Collaborative Cross. Genetics 206:537-556,
[doi:10.1534/genetics.116.198838](https://doi.org/10.1534/genetics.116.198838).



### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/CC/cc.zip")
cc <- read_cross2(file)
```
