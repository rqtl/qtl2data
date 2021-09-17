## Collaborative Cross mice

These are genotypes and related information for 69 mouse
Collaborative Cross lines, derived from the supplemental data for
Srivastava et al. (2017) Genomes of the Mouse Collaborative Cross. Genetics 206:537-556,
[doi:10.1534/genetics.116.198838](https://doi.org/10.1534/genetics.116.198838).

The supplementary data is at Zenodo,
[doi:10.5281/zenodo.377036](https://doi.org/10.5281/zenodo.377036).

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

In addition, [Table S2 at the journal
website](http://www.genetics.org/highwire/filestream/438137/field_highwire_adjunct_files/10/TableS2.xlsx)
contains important information, including the funnel codes for 55 of
the 69 strains. This also contains mtDNA and Y chromosome columns,
which are not entirely consistent with the `CCStrains.csv` file at
zenodo.

We'll use the information in Table S2 and will insert arbitrary cross
direction information for the lines where the funnel codes are not
available.

**Note**: Three strains (CC031, CC037, and CC056) have some problems
where the founder that contributed their Y chromosome is also present
on their X chromosome. This is not supposed to happen, and is
inconsistent with the way R/qtl2 handles 8-way RIL data. So we've
fudged the cross information for these strains so that the X
chromosome genotype probabilities won't be messed up. Also, a number
of strains are indicated to have been formed from only a subset of the
eight founders. Our cross information will still include all eight
founders, because the `"risib8"` cross type requires it.

The founder genotypes are from FigShare
([doi:10.6084/m9.figshare.5404762.v2](https://doi.org/10.6084/m9.figshare.5404762.v2))
and combine markers from the MegaMUGA and the GigaMUGA arrays. We use
the array annotation information from <https://github.com/kbroman/MUGAarrays>.

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
               "qtl2data/main/CC/cc.zip")
cc <- read_cross2(file)
```
