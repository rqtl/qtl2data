# download the raw data files for Dell'acqua et al (2015)

library(here)

# raw data dir
rawdata_dir <- here("Rawdata")

if(!dir.exists(rawdata_dir)) {
    dir.create(rawdata_dir)
}

# URLs of files
urls <- c("https://ndownloader.figshare.com/files/2096751",
          "https://ndownloader.figshare.com/files/2096747",
#         "https://ndownloader.figshare.com/files/2081430",
          "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0716-z/MediaObjects/13059_2015_716_MOESM1_ESM.xlsx",
          "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0716-z/MediaObjects/13059_2015_716_MOESM2_ESM.xlsx",
          "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0716-z/MediaObjects/13059_2015_716_MOESM8_ESM.xlsx",
          "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0716-z/MediaObjects/13059_2015_716_MOESM13_ESM.xlsx")

# file names
files <- c("MMfounders.geno.gz",      # founder genotypes
           "MMlines.geno.gz",         # MAGIC lines' genotypes
#          "MM.imputed.SNP.vcf.bz2",  # imputed founder genotypes
           "TableS1.xlsx",            # details of founder lines
           "TableS2.xlsx",            # MM population breeding design
           "TableS3.xlsx",            # genetic map
           "TableS6.xlsx")            # phenotypes

# download the files, if not already done
stopifnot(length(urls) == length(files))
for(i in seq_along(urls)) {
    file <- file.path(rawdata_dir, files[i])
    if(!file.exists(file)) {
        download.file(urls[i], file)
    }
}
