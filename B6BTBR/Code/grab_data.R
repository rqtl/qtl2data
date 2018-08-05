# grab the B6BTBR data and organize for R/qtl2

library(data.table)

# full data at
#     https://phenomedoc.jax.org/QTL_Archive/attie_2015/Attie_2015_eqtl_clean.zip
# unzip into Clean/

# genotypes
g <- fread("Clean/genotypes_clean.csv", data.table=FALSE)
g <- g[-(1:2),] # drop chr and position rows

covar <- g[,1:3] # MouseNum, Sex, and pgm (Sex Male/Female; pgm all 1's
rownames(covar) <- covar[,1]
g <- g[,-(2:3)] # drop Sex and pgm


# batches
tissues <- c("adipose", "gastroc", "hypo", "islet", "kidney", "liver")
for(tissue in tissues) {
    batch <- fread(paste0("Clean/", tissue, "_batch_clean.csv"), data.table=FALSE)
    dates <- colnames(batch)[-1]
    batches <- rep("", nrow(batch))
    names(batches) <- batch[,1]
    for(i in 2:ncol(batch))
        batches[batch[,i]==1] <- dates[i-1]
    batches[batches==""] <- "other"
    batches <- batches[names(batches) %in% covar$MouseNum]
    covar <- cbind(covar, rep(NA, nrow(covar)))
    covar[names(batches),ncol(covar)] <- batches
    colnames(covar)[ncol(covar)] <- paste0(tissue, "_batch")
}

# write covariates
write.table(covar, file="b6btbr_covar.csv", row.names=FALSE, col.names=TRUE,
            quote=FALSE, sep=",")

# genetic and physical maps
pmap <- fread("Clean/markers_physical_map.csv", data.table=FALSE, colClasses=rep("character", 3))
gmap <- fread("Clean/markers_genetic_map.csv", data.table=FALSE, colClasses=rep("character", 2))
gmap <- cbind(pmap[,1,drop=FALSE], gmap)
gmap[is.na(match(gmap$marker, colnames(g))),]
write.table(gmap, file="b6btbr_gmap.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(pmap, file="b6btbr_pmap.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

g <- g[,!is.na(match(colnames(g), c(gmap$marker, "MouseNum")))]
write.table(g, file="b6btbr_geno.csv", row.names=FALSE, col.names=TRUE,
            quote=FALSE, sep=",") # codes BB, BR, RR, -

# array data
for(tissue in tissues)
    system(paste0("mv Clean/", tissue, "_mlratio_clean.csv ", tissue, ".csv"))

# microarray anootation
annot <- fread("Clean/microarray_annot.csv", data.table=FALSE)
annot <- annot[,c("a_gene_id", "chr", "pos.Mb", "pos.cM", "pStart", "pEnd", "pStrand",
                  "officialgenesymbol", "SNPids", "NumofSNPs", "gobp.id", "gobp.term",
                  "gocc.id", "gocc.term", "gomf.id", "gomf.term", "path_id", "kegg.term",
                  "probe_use", "probesequence")]
write.table(annot, "b6btbr_microarray_annotation.csv", row.names=FALSE, col.names=TRUE, sep=",",
            quote=TRUE)

# clinical phenotypes
load("~/Projects/Attie/GoldStandard/FinalData/lipomics_final_rev2.RData")
phe <- lipomics[,c("MouseNum", "INSULIN (ng/ml) 10 wk", "AGOUTI COAT (Y/N)", "TUFT COAT (Y/N)")]
phe[,2] <- log10(phe[,2])
colnames(phe)[2] <- "log10_insulin_10wk"
phe[,3] <- match(phe[,3], c("B","T"))-1
colnames(phe)[3] <- "agouti_tan"
phe[,4] <- match(phe[,4], c("N","Y"))-1
colnames(phe)[4] <- "tufted"

phe <- phe[!is.na(match(phe$MouseNum, covar$MouseNum)),]
rownames(phe) <- phe[,1]
phe <- phe[rownames(covar),]

write.table(phe, "b6btbr_pheno.csv", col.names=TRUE, row.names=FALSE,
            quote=FALSE, sep=",")
