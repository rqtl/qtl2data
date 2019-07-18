# prepare the BXD data for R/qtl2

library(qtl2convert)
library(readxl)
library(qtl2)

# read the genotypes (and maps)
g <- readxl::read_excel("../RawData/BXD_Geno-19Jan2017_forGN.xlsx", skip=20, col_types="text")
g <- as.data.frame(g[g$Chr %in% c(1:19,"X"),])

# grab genetic and physical maps
gmap <- data.frame(g[,c("Locus","Chr","cM_BXD")])
colnames(gmap) <- c("marker", "chr", "pos")
rownames(gmap) <- gmap[,"marker"]

pmap <- data.frame(g[,c("Locus","Chr","Mb_mm10")])
colnames(pmap) <- c("marker", "chr", "pos")
rownames(pmap) <- pmap[,"marker"]

# fix differences in order; take physical map as correct and swap genetic positions
# chr 3: swap rs6301139, UNC5348732
gmap[c("rs6301139", "UNC5348732"), "pos"] <- rev(gmap[c("rs6301139", "UNC5348732"), "pos"])
# chr 3: swap rs30452147, rs31526621
gmap[c("rs30452147", "rs31526621"), "pos"] <- rev(gmap[c("rs30452147", "rs31526621"), "pos"])
# chr 5: move rs52382600 two down (after rs3659098, rs31502646)
gmap[c("rs52382600", "rs3659098", "rs31502646"), "pos"] <-
    gmap[c("rs31502646", "rs52382600", "rs3659098"), "pos"]
# chr 16: swap rs4219893, rs4219613
gmap[c("rs4219893", "rs4219613"), "pos"] <- rev(gmap[c("rs4219893", "rs4219613"), "pos"])
# chr 18: move rs6377403 two earlier (before UNCHS045463, UNCHS045465)
gmap[c("UNCHS045463", "UNCHS045465", "rs6377403"), "pos"] <-
    gmap[c("UNCHS045465", "rs6377403", "UNCHS045463"), "pos"]

# reorder genetic and physical map by physical order
pmap <- pmap[order(factor(pmap$chr, c(1:19,"X")), as.numeric(pmap$pos)),]
gmap <- gmap[rownames(pmap),]

# write genetic and physical maps
write2csv(gmap, "../bxd_gmap.csv", comment="Genetic map for BXD data", overwrite=TRUE)
write2csv(pmap, "../bxd_pmap.csv", comment="Physical map (mm10 Mbp) for BXD data", overwrite=TRUE)

# grab genotypes
g <- g[, !(colnames(g) %in% c("Chr", "Mb_mm9", "Mb_mm10", "cM_BXD",
                              "B6D2F1", "D2B5F1", "C57BL/6J", "DBA/2J"))]
colnames(g)[colnames(g)=="Locus"] <- "marker"
rownames(g) <- g[,"marker"]

# reorder markers as in the maps
g <- g[rownames(pmap),]

# write the genotypes
write2csv(g, "../bxd_geno.csv", comment="Genotypes for BXD data", overwrite=TRUE)

# write cross info (BxD for all strains)
crossinfo <- data.frame(id=names(g)[-1], cross_direction=rep("BxD", ncol(g)-1))
write2csv(crossinfo, "../bxd_crossinfo.csv", overwrite=TRUE,
          comment=paste0("Cross info for BXD data\n#",
                         "(all lines formed from cross between female B and male D)"))

# read the phenotypes
ph <- scan("../RawData/phenotypes_export_2018-02-27.txt", what=character(), sep="\n", skip=10)
ph <- strsplit(ph, "\t")
# look for "^BXD.+_Value$" for columns with the mean phenotypes
# (7) Description: seems blank
# (9) PubMed_ID: omit those with blank or "None"
# (10) Phenotype: really a description
#
# too many rows: (all with 633 rows)
# 1130 2425 2426 2427 2428 2429 2437 2438 2439 2440 2441 2442 2443 4981 5085 5175 5176 5182 5191 5211 5234
# seems like delete 10th field in 5182 and 26th in the rest
ph[[5182]] <- ph[[5182]][-10]
for(i in which(sapply(ph, length) != 632))
    ph[[i]] <- ph[[i]][-26]
stopifnot( all(sapply(ph, length) == 632) )

# columns with all blanks:
blank_cols <- which(sapply(seq_along(ph[[2]]), function(i) all(sapply(ph[-1], "[", i) == "")))
# 6-8, 12-25
ph <- lapply(ph, function(a) a[-blank_cols])

# columns with no variation
n_el <- sapply(seq_along(ph[[2]]), function(i) length(unique(sapply(ph[-1], "[", i))))
# no variation in values: Species, Cross, Database

# column 5: RecordID
# column 7: Phenotype (description)
# columns with values: grep("^BXD.+_Value$", ph[[1]])

phenocovar <- data.frame(id=vapply(ph, "[", "", 5)[-1],
                         description=vapply(ph, "[", "", 7)[-1])

phecol <- grep("^BXD.+_Value$", ph[[1]])
pheno <- matrix(nrow=length(phecol), ncol=length(ph)-1)
dimnames(pheno) <- list(sub("_Value$", "", ph[[1]][phecol]),
                        phenocovar$id)
for(i in seq_along(phecol)) {
    tmp <- vapply(ph, "[", "", phecol[i])[-1]
    tmp[tmp=="None"] <- NA
    pheno[i,] <- as.numeric(tmp)
}

# keep only the phenotypes with >= 10 lines
keep <- (colSums(!is.na(pheno)) >= 10)
phenocovar <- phenocovar[keep,]
pheno <- pheno[,keep]

# commas to semicolons, in descriptions
phenocovar$description <- gsub(",", ";", phenocovar$description)

# write to files
write2csv(phenocovar, "../bxd_phenocovar.csv", overwrite=TRUE,
          "Phenotype covariates (metadata) for BXD phenotype data")
write2csv(cbind(id=rownames(pheno), pheno), "../bxd_pheno.csv", overwrite=TRUE,
          "BXD phenotype data")

# write json file
write_control_file("../bxd.json",
                   description="BXD mouse data from GeneNetwork",
                   crosstype="risib",
                   geno_file="bxd_geno.csv",
                   geno_transposed=TRUE,
                   geno_codes=list(B=1, D=2),
                   xchr="X",
                   gmap_file="bxd_gmap.csv",
                   pmap_file="bxd_pmap.csv",
                   pheno_file="bxd_pheno.csv",
                   phenocovar_file="bxd_phenocovar.csv",
                   crossinfo_file = "bxd_crossinfo.csv",
                   crossinfo_codes = c("BxD"=0),
                   alleles=c("B", "D"),
                   overwrite=TRUE)

# test reading the data
x <- read_cross2("../bxd.json")
