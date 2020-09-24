# for the Dell'Acqua et al (2015) maize MAGIC population
# combine the array-based on sequence based founder genotypes

# in MMfounders.geno, CML91 is highly heterozygous; that's not the case in the imputed data
# the imputed data has the problem that the alleles may be from the opposite strand
# also, markers may appear twice

library(here)
library(data.table)
library(qtl2convert)
library(parallel)
n_cores <- parallel::detectCores()

# load array-based founder genotypes
founder_gfile <- here("Rawdata", "MMfounders.geno")
cleanup <- FALSE
if(!file.exists(founder_gfile)) {
    system(paste("gunzip -k", founder_gfile))
    cleanup <- TRUE
}

# load data
fg <- data.table::fread(founder_gfile, data.table=FALSE)
rownames(fg) <- fg$LINE
fg <- fg[,colnames(fg) != "LINE"]
if(cleanup) unlink(founder_gfile)

# pull out founder names
founders <- sapply(strsplit(rownames(fg), "\\."), "[", 1)
ufounders <- unique(founders)

# get consensus genotypes
ufg <- matrix(nrow=ncol(fg), ncol=length(ufounders))
dimnames(ufg) <- list(colnames(fg), ufounders)
for(u in ufounders) {
    ufg[,u] <- find_consensus_geno(t(fg[founders==u,]), cores=n_cores, na.strings=c("H", "N"))
}

# load "VCF file" (though not really a VCF)
vcf_file <- here("Rawdata", "MM.imputed.SNP.vcf")
cleanup <- FALSE
if(!file.exists(vcf_file)) {
    system(paste0("bzip2 -kd ", vcf_file, ".bz2"))
    cleanup <- TRUE
}
fgimp <- data.table::fread(vcf_file, data.table=FALSE)
unlink(vcf_file)
# chromosome:position (in bp)
vcf_pos <- paste(fgimp$"#CHROM", fgimp$POS, sep=":")

# load map file
map <- as.data.frame( readxl::read_excel( here("Rawdata", "TableS3.xlsx") ) )
map$"Position (RefGen V3)" <- round(map$"Position (RefGen V3)" * 1e6)
map_pos <- paste(map$Chr, map$"Position (RefGen V3)", sep=":")

# subset imputed markers to those on array
fgimp <- fgimp[vcf_pos %in% map_pos,]
vcf_pos <- vcf_pos[vcf_pos %in% map_pos]

# combine common rows
tab <- table(vcf_pos)
dups <- names(tab)[tab>1]
to_omit <- NULL
for(dup in dups) {
    wh <- which(vcf_pos==dup)
    # pull out duplicate rows
    dup_rows <- fgimp[wh,]

    # save just highest quality
    dup_rows <- dup_rows[dup_rows$QUAL == max(dup_rows$QUAL),]

    # find consensus genotype calls
    if(nrow(dup_rows)>1) {
        cg <- qtl2convert::find_consensus_geno(t(dup_rows[,10:18]), na.strings=c("./.", "0/1", "1/0"))
        fgimp[wh[1],10:18] <- cg
    }

    # keep track of rows to omit
    to_omit <- c(to_omit, wh[-1])
}

# drop duplicate rows
fgimp <- fgimp[-to_omit,]

# replace 0/0 and 1/1 with REF and ALT
for(i in 10:18) {
    fgimp[!is.na(fgimp[,i]) & fgimp[,i]=="0/0",i] <- fgimp[!is.na(fgimp[,i]) & fgimp[,i]=="0/0","REF"]
    fgimp[!is.na(fgimp[,i]) & fgimp[,i]=="1/1",i] <- fgimp[!is.na(fgimp[,i]) & fgimp[,i]=="1/1","ALT"]
}


# focus on the common markers
fgimp_pos <- paste(fgimp$"#CHROM", fgimp$POS, sep=":")
markers <- map$Marker[match(fgimp_pos, map_pos)]
stopifnot(all(map$Marker == rownames(ufg)))
fgsub <- ufg[markers,]

# reorder fgimp columns
fgimp[,10:18] <- fgimp[,colnames(fgsub)]
colnames(fgimp)[10:18] <- colnames(fgsub)

# try to match alleles
ug <- find_unique_geno(fgsub)
mismatch <- !is.na(ug[,1]) & !((ug[,1] == fgimp$REF & ug[,2]==fgimp$ALT) |
                               (ug[,1] == fgimp$ALT & ug[,2]==fgimp$REF))
# ignore C/G or A/T cases because we can't tell
omit <- NULL
for(i in which(mismatch)) {
    if(paste0(fgimp[i,"REF"], fgimp[i,"ALT"]) %in% c("CG", "GC", "AT", "TA",
                                                     "A.", "T.", "C.", "G.")) {
        omit <- c(omit, i)
    }
}
# drop the data for this one
fgimp[omit,10:18] <- NA
mismatch[omit] <- FALSE

swapped <- c(A="T",C="G",G="C",T="A")

# swap the others
for(i in which(mismatch)) {
    ref <- fgimp[i,"REF"]
    alt <- fgimp[i,"ALT"]
    g <- fgimp[i,10:18]

    g[!is.na(g) & g==ref] <-  fgimp[i,"REF"] <- swapped[ref]
    g[!is.na(g) & g==alt] <- fgimp[i,"ALT"] <- swapped[alt]
    fgimp[i,10:18] <- g
}


# for data missing from array data, plug in data from imputed
dont_sub <- is.na(ug[,1])
for(i in 1:9) {
    imp <- fgimp[,i+9]
    imp[imp=="./." | imp=="0/1"] <- NA
    tosub <- is.na(fgsub[,i]) & !dont_sub
    fgsub[tosub,i] <- imp[tosub]
}

# plug back into ufg
ufg[markers,] <- fgsub

# write to CSV file
write2csv(ufg, here("Rawdata", "foundergeno_imputed.csv"),
          comment=paste("combined founder genotypes for Maize MAGIC from ",
                        "Dell'Acqua et al (2015) doi:10.1186/s13059-015-0716-z"),
          row.names="id", overwrite=TRUE)
