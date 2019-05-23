## BXD data

Downloaded from genenetwork.org

- genotypes:
  - select "Type: Genotypes" and "Data Set: BXD Genotypes"
  - click on "Info" button next to "Data Set: BXD Genotypes"
  - There'll be a link to the data on the right under "Download
    datasets and supplementary data files"
  - Using the xlsx file,
    [`BXD_Geno-19Jan2017_forGN.xlsx`](http://datafiles.genenetwork.org/download/GN600/BXD_Geno-19Jan2017_forGN.xlsx)

- `BXD.geno`:
  - is what is *actually* used in GeneNetwork for mapping
  - grabbed using the Rest API, <http://gn2-zach.genenetwork.org/api/v_pre1/genotypes/BXD>
  - has same marker order problem as the xlsx file

    ```
    x <- read.table("BXD.geno", stringsAsFactors=FALSE, skip=21, sep="\t", header=TRUE)
    # chr 16: 251 ... swap 251/252
    # chr 18: 12 ... swap 12/13
    # chr 3: 177, 290  ... swap 177/178, swap 290/291
    # chr 5: 195 ... actually 195 belongs after 197
    ```

- phenotypes:
  - Rob Williams emailed them to me
  - He had selected "Type: phenotypes" and then "*" under "Get Any:"
  - Then click the "Search" button
  - Then click "Select" to select all
  - Then "Download Table"
  - Want to strip off the phenotypes that don't have PubMed IDs
