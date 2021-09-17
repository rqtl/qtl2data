## BXD mouse recombinant inbred lines

Genotype and phenotype data for the mouse BXD recombinant inbred
lines, from [GeneNetwork](http://gn2.genenetwork.org).

### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).



### Citation

Need to add a citation.



### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/main/BXD/bxd.zip")
bxd <- read_cross2(file)
```
