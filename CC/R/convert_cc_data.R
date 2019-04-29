# download and convert the Collaborative Cross data
#
# supplemental data for Srivastava et al. (2017) Genomes of the Mouse
# Collaborative Cross. Genetics 206:537-556, doi:10.1534/genetics.116.198838
#
# available at Zenodo, doi:10.5281/zenodo.377036

# required libraries
library(data.table)
library(qtl2)
library(qtl2convert)
library(broman)
library(readxl)
set.seed(83763628)


# create RawData/ if not available
rawdata_dir <- "../RawData"
if(!dir.exists(rawdata_dir)) {
    dir.create(rawdata_dir)
}



# download Prob36.zip if not available
zenodo_url <- "https://zenodo.org/record/377036/files/"
zenodo_postpend <- "?download=1"
prob_file <- "Prob36.zip"
distant_file <- paste0(zenodo_url, prob_file, zenodo_postpend)
local_file <- file.path(rawdata_dir, prob_file)
if(!file.exists(local_file)) {
    message("Downloading Prob36.zip")
    download.file(distant_file, local_file)
}

# create directory for data if it doesn't exist
prob_dir <- file.path(rawdata_dir, "Prob36")
if(!dir.exists(prob_dir)) {
    dir.create(prob_dir)
}

# unzip if it hasn't been unzipped
gzfile <- file.path(prob_dir, "CC001-Uncb38V01.csv.gz")
csvfile <- file.path(prob_dir, "CC001-Uncb38V01.csv")
if(!file.exists(gzfile) && !file.exists(csvfile)) {
    message("unzipping Prob36.zip")
    unzip(local_file, exdir=prob_dir)
}

# gunzip all of the Prob36 files
if(!file.exists(csvfile)) {
    files <- list.files(prob_dir, pattern=".csv.gz$")

    message("unzipping the probability files")

    for(file in files) {
        system(paste("gunzip", file.path(prob_dir, file)))
    }
}

# load the Prob36 files and determine X, Y, and M genotypes
message("reading the probability files")
files <- list.files(prob_dir, pattern=".csv$")
strains <- sub("\\.csv$", "", files)
probs <- setNames(vector("list", length(strains)), strains)
for(i in seq_along(files)) {
    probs[[i]] <- data.table::fread(file.path(prob_dir, files[i]), data.table=FALSE)
}

##############################
# guess the rest of the cross order
##############################
message("inferring cross order")
mprob <- t(sapply(probs, function(a) colMeans(a[a[,2]=="M", paste0(LETTERS, LETTERS)[1:8]])))
yprob <- t(sapply(probs, function(a) colMeans(a[a[,2]=="Y", paste0(LETTERS, LETTERS)[1:8]])))
xprob <- t(sapply(probs, function(a) colMeans(a[a[,2]=="X", paste0(LETTERS, LETTERS)[1:8]])))
# a bunch where we can't tell Y or M

# also need the supplementary data file
supp_file <- "SupplementalData.zip"
distant_file <- paste0(zenodo_url, supp_file, zenodo_postpend)
local_file <- file.path(rawdata_dir, supp_file)
if(!file.exists(local_file)) {
    message("downloading SupplementalData.zip")
    download.file(distant_file, local_file)
}

# extract just the CCStrains.csv file
csv_file <- "SupplmentalData/CCStrains.csv"
if(!file.exists(csv_file)) {
    message("unzipping SupplementalData.zip")
    unzip(local_file, csv_file, exdir=rawdata_dir)
}
csv_file <- file.path(rawdata_dir, csv_file)
ccstrains <- data.table::fread(csv_file, data.table=FALSE)

# check that the strain names are the same
stopifnot( all( paste0(sub("/", "-", ccstrains$Strain, fixed=TRUE), "b38V01") ==
                strains))

# download table S2
url <- "http://www.genetics.org/highwire/filestream/438137/field_highwire_adjunct_files/10/TableS2.xlsx"
file <- basename(url)
local_file <- file.path("..", "RawData", file)
if(!file.exists(local_file)) {
    download.file(url, local_file)
}
tab <- as.data.frame(readxl::read_xlsx(local_file))

# check that the strain names are the same
stopifnot( all( paste0(sub("/", "-", tab$Strain, fixed=TRUE), "b38V01") ==
                strains))

funnel <- tab[,"Funnel Code"]
mtdna <- tab[,"Mitochondria"]
ychr <- tab[,"Chromosome Y"]

# determine cross orders; use orders in Table S2 when available
cross_info <- matrix(ncol=8, nrow=length(strains))
dimnames(cross_info) <- list(ccstrains$Strain, LETTERS[1:8])
for(i in which(!is.na(funnel))) {
    cross_info[i,] <- match(unlist(strsplit(funnel[i], "")), LETTERS)
}

# when not available, use the values in table S2
use_mtdna <- mtdna; use_mtdna[mtdna=="A/D" | mtdna=="D/A"] <- "A" # use A when A/D or D/A
cross_info[is.na(cross_info[,1]),1] <- match(use_mtdna, LETTERS[1:8])[is.na(cross_info[,1])]
cross_info[is.na(cross_info[,8]),8] <- match(ychr, LETTERS[1:8])[is.na(cross_info[,8])]

# problems: CC013, CC023, CC027
# CC013 M = Y = E [genotypes say Y could be B or C, too] --- change Y to B

cross_info["CC013/GeniUnc", 8] <- 2

# check cases where mitochondrial genotype is clear
max_mprob <- apply(mprob, 1, max)
wh_mprob <- apply(mprob, 1, which.max)
stopifnot( all(cross_info[max_mprob>0.9,1] == wh_mprob[max_mprob > 0.9]) )

stopifnot( all( cross_info[,1] != cross_info[,8] ))

# are there any obvious differences?

# mitochondria
stopifnot( all( sort(max_mprob[wh_mprob != cross_info[,1]]) < 0.5))

# now the Y chr
max_yprob <- apply(yprob, 1, max)
wh_yprob <- apply(yprob, 1, which.max)
stopifnot( all( sort(max_yprob[wh_yprob != cross_info[,8]]) < 0.5))

# founders with largest X chr probabilities
wh_xprob <- t(apply(xprob, 1, order, decreasing=TRUE))

# most common X chr genotype that's not the mtDNA one in the 3rd slot
cross_info[is.na(cross_info[,3]),3] <- wh_xprob[is.na(cross_info[,3]),1]
stopifnot( all(cross_info[,1] != cross_info[,3]) )

# sort other X chr probs...put three most probable in the 2nd, 5th, and 6th slots
for(i in which(is.na(cross_info[,2]))) {
    z <- (wh_xprob[i, ] %wnin% cross_info[i, c(1,3,8)])
    cross_info[i, c(2,5,6)] <- sample(z[1:3])
    cross_info[i, c(4,7)] <- sample(z[4:5])
}

# further problems:
# CC031 Y chr is B but this is clearly on the X chromosome
# CC037 Y chr is D but this is clearly on the X chromosome
# CC056 Y chr is E but this is clearly on the X chromosome
# ... swap these with one of 2,5,6

cross_info["CC031/GeniUnc", ] <- c(1,6,5,3,7,2,4,8)
cross_info["CC037/TauUnc", ] <- c(3,4,8,7,6,1,5,2)
cross_info["CC056/GeniUnc", ] <- c(3,8,1,7,5,4,6,2)

stopifnot( all(apply(cross_info, 1, sort) == 1:8) )

message("writing cross and covariate data")

# write the cross information to a file
write2csv(cbind(id=rownames(cross_info), cross_info),
          "../cc_crossinfo.csv",
          comment=paste("Cross information for Collaborative Cross (CC) lines inferred from",
                        "genotypes from Srivastava et al. (2017)",
                        "doi:10.1534/genetics.116.198838,",
                        "data at doi:10.5281/zenodo.377036"),
          overwrite=TRUE)


# write covariate info with M and Y as inferred
covar <- data.frame(id=rownames(cross_info),
                    mitochondria=ccstrains$Mitochondria,
                    Ychr=ccstrains$ChrY,
                    n_founders=ccstrains$N_Founders,
                    origin=ccstrains$Origin,
                    stringsAsFactors=FALSE)

write2csv(covar, "../cc_covar.csv",
          comment=paste("Covariate information for Collaborative Cross (CC) lines inferred from",
                        "genotypes from Srivastava et al. (2017)",
                        "doi:10.1534/genetics.116.198838,",
                        "data at doi:10.5281/zenodo.377036"),
          overwrite=TRUE)

# download genotypes.zip if not available
genotype_file <- "genotypes.zip"
distant_file <- paste0(zenodo_url, genotype_file, zenodo_postpend)
local_file <- file.path(rawdata_dir, genotype_file)
if(!file.exists(local_file)) {
    message("downloading genotypes.zip")
    download.file(distant_file, local_file)
}

# create directory for data if it doesn't exist
genotype_dir <- file.path(rawdata_dir, "genotypes")
if(!dir.exists(genotype_dir)) {
    dir.create(genotype_dir)
}

# unzip if it hasn't been unzipped
file <- file.path(genotype_dir, "SEQgenotypes.csv")
if(!file.exists(file)) {
    message("unzipping genotypes")
    unzip(local_file, exdir=genotype_dir)
}

# read genotypes
message("reading genotypes")
g <- data.table::fread(file, data.table=FALSE)
g[g=="N" | g=="H"] <- NA


# founder genotypes from figshare
# https://doi.org/10.6084/m9.figshare.5404762.v2
fg_url <- "https://ndownloader.figshare.com/files/13623080"
local_file <- file.path(rawdata_dir, "MMnGM_processed_files.zip")
if(!file.exists(local_file)) {
    message("downloading founder genotypes")
    download.file(fg_url, local_file)
}

# unzip the founder genotype data
if(!file.exists(file.path("..", "MMnGM", "MMnGM_info.csv"))) {
    message("unzipping founder genotypes")
    unzip(local_file, exdir="..")
}

# load the allele codes
fga <- read_csv(file.path("..", "MMnGM", "MMnGM_allelecodes.csv"))
fga <- fga[fga$chr %in% c(1:19,"X"),]

# omit markers not in founder genotypes
g <- g[g$marker %in% rownames(fga),]
fga <- fga[rownames(fga) %in% g$marker,]

# reorder rows of genotype data to match fga
g <- g[match(g$marker, rownames(fga)),]

# cut down to just the genotypes
rownames(g) <- g[,1]
g <- g[,-(1:3)]

# strip off individual IDs from column names
colnames(g) <- paste0(sapply(strsplit(colnames(g), "Unc"), "[", 1), "Unc")
stopifnot( all( colnames(g) == ccstrains$Strain ))

# encode genotypes
message("encoding genotypes")
g <- encode_geno(g, fga[,c("A","B")], cores=0)

# omit markers with no data
g <- g[rowSums(g!="-") > 0, , drop=FALSE]

# write genotypes to files, one chromosome at a time
message("writing genotypes")
for(chr in c(1:19,"X")) {
    file <- paste0("../cc_geno", chr, ".csv")
    gsub <- g[rownames(g) %in% rownames(fga)[fga$chr==chr], , drop=FALSE]
    write2csv(cbind(marker=rownames(gsub), gsub), file,
              overwrite=TRUE,
              comment=paste("Chromosome", chr, "genotypes",
                        "for Collaborative Cross (CC) lines inferred from",
                        "genotypes from Srivastava et al. (2017)",
                        "doi:10.1534/genetics.116.198838,",
                        "data at doi:10.5281/zenodo.377036"))
}

# create JSON file
chr <- c(1:19,"X")
message("Write control file")
write_control_file("../cc.json", crosstype="risib8",
                   geno_file=paste0("cc_geno", chr, ".csv"),
                   founder_geno_file=paste0("MMnGM/MMnGM_foundergeno", chr, ".csv"),
                   gmap_file=paste0("MMnGM/MMnGM_gmap", chr, ".csv"),
                   pmap_file=paste0("MMnGM/MMnGM_pmap", chr, ".csv"),
                   covar_file="cc_covar.csv",
                   crossinfo_file="cc_crossinfo.csv",
                   geno_codes=c("A"=1, "B"=3),
                   alleles=LETTERS[1:8],
                   xchr="X",
                   geno_transposed=TRUE,
                   founder_geno_transposed=TRUE,
                   description=paste("Data for Collaborative Cross (CC) lines",
                                     "from Srivastava et al. (2017)",
                                     "doi:10.1534/genetics.116.198838,",
                                     "data at doi:10.5281/zenodo.377036"),
                   overwrite=TRUE)

# create zip file
message("Create zip file")
zip_datafiles("../cc.json", overwrite=TRUE)
