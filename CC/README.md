## Collaborative Cross mice

These are genotypes and related information for the mouse
Collaborative Cross, derived from the supplemental data for
Srivastava et al. (2017) Genomes of the Mouse Collaborative Cross. Genetics 206:537-556,
[doi:10.1534/genetics.116.198838](https://doi.org/10.1534/genetics.116.198838).

The supplementary data is at Zenodo,
[doi:10.5281/zenodo.377036](https://doi.org/10.5281/zenodo.377036).

We'll use the following:

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


Also of interest:

- `GigaHaps.zip`, `MRCAsHaps.zip`, `SeqHaps.zip`

   - A `.hap` file for each strain which is a CSV file with some
     header information and then two rows for each chromosome with
     alleles and locations of breakpoints

The script [`R/convert_cc_data.R`](R/convert_cc_data.R) both downloads
and converts the data into a form useful for R/qtl2.