# convert maize MAGIC data for R/qtl2

library(here)
library(qtl2)
library(qtl2convert)
library(data.table)
library(readxl)
library(parallel)
n_cores <- parallel::detectCores()

rawdata_dir <- here("Rawdata")

##############################
# founder genotypes
##############################
founder_gfile <- file.path(rawdata_dir, "foundergeno_imputed.csv")
fg <- read_csv(founder_gfile)

# pull out founder names
founders <- colnames(fg)

# find unique genotypes
ug <- find_unique_geno(fg)

##############################
# MAGIC genotypes
##############################
magic_gfile <- file.path(rawdata_dir, "MMlines.geno")
cleanup <- FALSE
if(!file.exists(magic_gfile)) {
    system(paste("gunzip -k", magic_gfile))
    cleanup <- TRUE
}

# load data
g <- data.table::fread(magic_gfile, data.table=FALSE)
if(cleanup) unlink(magic_gfile)

rownames(g) <- g[,1]
g <- t(g[,-1])

stopifnot(all(rownames(g) == rownames(fg)))

# encode genotypes
fg <- encode_geno(fg[!is.na(ug[,1]),], ug[!is.na(ug[,1]),], cores=n_cores)
g <- encode_geno(g[!is.na(ug[,1]),], ug[!is.na(ug[,1]),], cores=n_cores)

##############################
# maps
##############################
map <- as.data.frame( readxl::read_excel(file.path(rawdata_dir, "TableS3.xlsx")) )
# drop the uninformative markers
map <- map[map$Marker %in% rownames(g), ]
stopifnot(all(rownames(g) == map$Marker))

pmap <- map[,1:3]
colnames(pmap) <- c("marker", "chr", "pos_Mbp_RefGenV3")

gmap <- map[,c(1:2,4)]
colnames(pmap) <- c("marker", "chr", "pos_cM")


##############################
# phenotype data
##############################
pheno <- as.data.frame( readxl::read_excel(file.path(rawdata_dir, "TableS6.xlsx")) )

# drop parents phenotypes
pheno <- pheno[!(pheno$Line %in% founders),]


##############################
# cross info
##############################
tabS2 <- as.data.frame( readxl::read_excel(file.path(rawdata_dir, "TableS2.xlsx")) )
family <- tabS2[,"8-ways (#)"]
cross <- tabS2[,"8-ways"]
stopifnot(all(is.na(family) == is.na(cross)))
cross <- cross[!is.na(family)]
family <- family[!is.na(family)]

# get rid of cross direction information
for(char in c("[", "(", "x", ")", "]", "/", "+")) {
    cross <- gsub(char, " ", cross, fixed=TRUE)
}

# count number of times each parent appears in each cross
parents <- strsplit(cross, "\\s+")
parents <- lapply(parents, function(a) table(factor(a[a != ""], levels=founders)))

# turn into a matrix
crossinfo <- matrix(unlist(parents), byrow=TRUE, ncol=length(parents[[1]]))
dimnames(crossinfo) <- list(family, names(parents[[1]]))

# expand to a matrix for all of the MAGIC lines
full_crossinfo <- matrix(ncol=ncol(crossinfo), nrow=ncol(g))
dimnames(full_crossinfo) <- list(colnames(g), rownames(full_crossinfo))
family <- sapply(strsplit(rownames(full_crossinfo), "_"), "[", 1)
ufamily <- rownames(crossinfo)
for(i in 1:nrow(crossinfo)) {
    for(j in which(family==ufamily[i])) {
        full_crossinfo[j,] <- crossinfo[i,]
    }
}

# multiply counts of parents by 100, and allow small chance for the missing parents
full_crossinfo <- full_crossinfo*100
full_crossinfo[full_crossinfo==0] <- 1
full_crossinfo <- cbind(ngen=4, full_crossinfo)



##############################
# write data
##############################

description <- "Maize MAGIC from Dell'Acqua et al (2015) doi:10.1186/s13059-015-0716-z"

## write maps
write2csv(pmap, here("maize_magic_pmap.csv"),
          comment=paste("Physical map for", description),
          overwrite=TRUE)

write2csv(gmap, here("maize_magic_gmap.csv"),
          comment=paste("Genetic map for", description),
          overwrite=TRUE)

## write founder genotypes
write2csv(fg, here("maize_magic_foundergeno.csv"),
          comment=paste("Founder genotypes for", description),
          overwrite=TRUE, row.names="id")

## write MAGIC genotypes
write2csv(g, here("maize_magic_geno.csv"),
          comment=paste("Genotypes for", description),
          overwrite=TRUE, row.names="id")

## write cross information
write2csv(full_crossinfo, here("maize_magic_crossinfo.csv"),
          comment=paste("Cross info for", description),
          overwrite=TRUE, row.names="id")

## write phenotypes
write2csv(pheno, here("maize_magic_pheno.csv"),
          comment=paste("Phenotypes for", description),
          overwrite=TRUE, row.names=NULL)

## write phenotype covariates
phenocovar <- data.frame(pheno=colnames(pheno)[-1],
                         description=c("pollen shed",
                                       "plant height",
                                       "ear height",
                                       "transformed grain yield"))
write2csv(phenocovar, here("maize_magic_phenocovar.csv"),
          comment=paste("Phenotype covariates for", description),
          overwrite=TRUE, row.names=NULL)

## write control file
write_control_file(here("maize_magic.json"),
                   crosstype="genril9",
                   geno_file="maize_magic_geno.csv",
                   founder_geno_file="maize_magic_foundergeno.csv",
                   gmap_file="maize_magic_gmap.csv",
                   pmap_file="maize_magic_pmap.csv",
                   pheno_file="maize_magic_pheno.csv",
                   crossinfo_file="maize_magic_crossinfo.csv",
                   phenocovar_file="maize_magic_phenocovar.csv",
                   geno_codes=c(A=1, B=3),
                   geno_transposed=TRUE,
                   founder_geno_transposed=TRUE,
                   description=description,
                   overwrite=TRUE)

zip_datafiles(here("maize_magic.json"), overwrite=TRUE)
