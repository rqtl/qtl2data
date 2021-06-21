# prepare Bult/Recla data for R/qtl2
# downloaded genotypes in R/qtl2 format from dodb.jax.org
# downloaded phenotypes from http://phenomedoc.jax.org/MPD_projdatasets/Recla2.csv
# publication doi: 10.1097/j.pain.0000000000001571

# required packages
library(data.table)
library(qtl2convert)
library(qtl2)
library(here)

# extract genotypes
zipfile <- here("Orig", "bult_recla_megamuga.zip")
rawdir <- here("RawData")
if(!dir.exists(rawdir)) dir.create(rawdir)
unzip(zipfile, exdir=rawdir)

# download MM annotations (drop markers that didn't map to chr 1-19,X
url <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/master/UWisc/mm_uwisc_v1.csv"
file <- here("Orig", basename(url))
if(!file.exists(file)) download.file(url, file)
mm_uwisc <- data.table::fread(file, data.table=FALSE)
mm_uwisc <- mm_uwisc[!is.na(mm_uwisc$chr) & !is.na(mm_uwisc$bp_mm10) & mm_uwisc$chr %in% c(1:19,"X"),]

# download phenotype data
url <- "http://phenomedoc.jax.org/MPD_projdatasets/Recla2.csv"
file <- here("Orig", basename(url))
if(!file.exists(file)) download.file(url, file)
phe <- data.table::fread(file, data.table=FALSE)
rownames(phe) <- paste0("CJB-DO-", phe$id)

# read genotype data (note: only 77,725 markers)
g <- data.table::fread(here("RawData/Bult-Bult_Recla-MegaMUGA_geno.csv"), data.table=FALSE)
rownames(g) <- g[,1]
g <- g[,colnames(g) %in% rownames(phe)]
g <- g[rownames(g) %in% mm_uwisc$marker,]

# read founder genotype data (also only 77,725 markers)
fg <- data.table::fread(here("RawData/Bult-Bult_Recla-MegaMUGA_foundergeno.csv"), data.table=FALSE)
rownames(fg) <- fg$marker
fg <- fg[rownames(fg) %in% mm_uwisc$marker,-1]

# tag for data
datatag <- "for Recla et al. (2019) DO data (https://doi.org/10.1097/j.pain.0000000000001571)"

# covariate information
covar <- cbind(phe[,c("sex", "am_or_pm")], dogen=8)
covar$sex <- toupper(covar$sex)
qtl2convert::write2csv(covar, here("recla2_covar.csv"),
                       row.names="id", overwrite=TRUE,
                       comment=paste("Covariate information", datatag))

# phenotype data
phe <- phe[,c("formalin_score", "formalin_score_sqrt")]
qtl2convert::write2csv(phe, here("recla2_pheno.csv"),
                       row.names="id", overwrite=TRUE,
                       comment=paste("Phenotypes", datatag))

# revised genetic and physical maps (reduce to markers in genotype data)
mm_uwisc <- mm_uwisc[mm_uwisc$marker %in% rownames(g),]
mm_uwisc$Mbp_mm10 <- mm_uwisc$bp_mm10/1e6
qtl2convert::write2csv(mm_uwisc[,c("marker", "chr", "cM_cox")],
                       here("recla2_gmap.csv"),
                       comment=paste("Genetic map", datatag),
                       overwrite=TRUE)
qtl2convert::write2csv(mm_uwisc[,c("marker", "chr", "Mbp_mm10")],
                       here("recla2_pmap.csv"),
                       comment=paste("Physical map (mm10)", datatag),
                       overwrite=TRUE)

# write genotype data
qtl2convert::write2csv(g, here("recla2_geno.csv"),
                       comment=paste("Genotype data", datatag),
                       row.names="marker", overwrite=TRUE)

# write founder genotype data
qtl2convert::write2csv(fg, here("recla2_foundergeno.csv"),
                       comment=paste("Founder genotype data", datatag),
                       row.names="marker", overwrite=TRUE)

# fix json file (covar and pheno -> v2)
qtl2::write_control_file(here("recla2.json"),
                         crosstype="do",
                         geno_file = "recla2_geno.csv",
                         founder_geno_file = "recla2_foundergeno.csv",
                         gmap_file = "recla2_gmap.csv",
                         pmap_file = "recla2_pmap.csv",
                         pheno_file = "recla2_pheno.csv",
                         covar_file = "recla2_covar.csv",
                         sex_covar = "sex",
                         sex_codes = list(F="female", M="male"),
                         geno_codes = list("1"="1", "2"="2", "3"="3"),
                         xchr = "X",
                         crossinfo_covar = "dogen",
                         geno_transposed=TRUE,
                         founder_geno_transposed=TRUE,
                         alleles = LETTERS[1:8],
                         description=paste("Control file", datatag),
                         overwrite=TRUE)




# zip files
qtl2::zip_datafiles(here("recla2.json"), overwrite=TRUE)
