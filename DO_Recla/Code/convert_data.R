## convert DO data from Recla et al (2014)
## for use with R/qtl2
#
# DO data from https://phenome.jax.org/projects/Recla1
# founder data from ftp://ftp.jax.org/MUGA/
# Additional phenotypes at http://phenomedoc.jax.org/MPD_projdatasets/Chesler4.csv

library(data.table)
library(qtl2)
library(qtl2convert)

phe <- data.table::fread("../Orig/Recla2013_pheno.csv",
                         data.table=FALSE, colClasses="character")
covar <- phe[,c(1,3:6)]
phe <- phe[,c(1,7:8)]
colnames(covar)[1] <- colnames(phe)[1] <- "id"
colnames(phe)[2:3] <- c("latency_to_respond", "log_latency")

# Chesler4: seems to be the same mice
phe2 <- data.table::fread("../Orig/Chesler4_table.csv", skip=5,
                          data.table=FALSE, colClasses="character")
phe2 <- phe2[phe2$strain=="J:DO",]
phe2 <- phe2[match(phe$id, phe2$id),]
covar$ngen <- substr(phe2$generation, 2, 2) # no. DO generations
covar$coat_color <- phe2$coat_color

stopifnot(all(covar$Subgroup == phe2$group))
stopifnot(all((covar$Sex=="female" & phe2$sex=="f") | (covar$Sex=="male" & phe2$sex=="m")))

# latency in the Recla data is the same as HC_latency in the Chesler data
phe <- phe2[,-c(1:2, 4:6)]

# genotypes
g <- data.table::fread("../Orig/Recla2013_geno.csv", data.table=FALSE, colClasses="character")
rownames(g) <- g[,1]
g <- g[,-(1:3)]
g[g=="--"] <- NA
g <- as.matrix(g)
ng <- apply(g, 1, function(a) length(unique(a[!is.na(a)]))) # number of unique genotypes
g <- g[ng == 3,]

# keep just the individuals that have both phenotypes and genotypes
gid <- colnames(g)
pid <- rownames(phe) <- rownames(covar) <- phe[,1]
pid <- pid[!is.na(match(pid, gid))]
g <- g[,pid]
phe <- phe[pid,]
covar <- covar[pid,]

# write phenotype and covariate data
write.table(phe, file="../recla_pheno.csv", quote=FALSE, sep=",",
            row.names=FALSE, col.names=TRUE)
write.table(covar, file="../recla_covar.csv", quote=FALSE, sep=",",
            row.names=FALSE, col.names=TRUE)

# founder info
load("../Orig/muga_code.Rdata")
load("../Orig/muga_snps.Rdata")
load("../Orig/muga_geno.Rdata")

# fix weirdness in muga_code names
names(muga_code) <- gsub("_X", "", names(muga_code))
stopifnot(all(names(muga_code) == colnames(muga_geno)))
stopifnot(all(rownames(muga_geno) == muga_snps[,1]))

founders2keep <- names(muga_code)[muga_code %in% paste0(LETTERS, LETTERS)[1:8]]
gmap <- muga_snps[,c("SNP_ID", "Chr", "cM")]
pmap <- muga_snps[,c("SNP_ID", "Chr", "Mb_NCBI38")]
colnames(gmap) <- colnames(pmap) <- c("marker", "chr", "pos")
muga_geno <- muga_geno[,founders2keep]
muga_code <- muga_code[founders2keep]

# marker names
marg <- rownames(muga_geno)
marm <- as.character(gmap[,1])
stopifnot(all(marg==marm))
rownames(gmap) <- rownames(pmap) <- marm

# I've dropped some markers from the DO data (with <3 unique genotypes); otherwise they're all there
mar_do <- rownames(g)
sum(is.na(match(mar_do, marm))) # 0
sum(is.na(match(marm, mar_do))) # had dropped some

# drop markers not in DO  data
muga_geno <- muga_geno[mar_do,]
gmap <- gmap[mar_do,]
pmap <- pmap[mar_do,]
g <- g[mar_do,]

# drop markers on M, Y, or P, or missing cM or Mbp
keep <- gmap$chr %in% c(1:19, "X") & !is.na(gmap$pos)
gmap <- gmap[keep,]
pmap <- pmap[rownames(gmap),]
pmap <- pmap[!is.na(pmap$pos),]
gmap <- gmap[rownames(pmap),]
muga_geno <- muga_geno[rownames(gmap),]
g <- g[rownames(gmap),]

# omit H's and N's in founder data
muga_geno[muga_geno=="N" | muga_geno=="H"] <- NA

founder_geno <- matrix(ncol=8, nrow=nrow(muga_geno))
dimnames(founder_geno) <- list(rownames(muga_geno), LETTERS[1:8])

pat <- paste0(LETTERS, LETTERS)[1:8]
for(i in 1:8) {
    founder_geno[,i] <- find_consensus_geno(muga_geno[,muga_code==pat[i]])
}

# count unique values in each row
nug <- apply(founder_geno, 1,  function(a) length(unique(a[!is.na(a)])))

# drop markers that don't have 2 observed genotypes in the founders
founder_geno <- founder_geno[nug==2,]
g <- g[rownames(founder_geno),]

# the pair of alleles at each marker
alleles <- find_unique_geno(founder_geno)

# convert genotypes using those alleles
founder_geno <- encode_geno(founder_geno, alleles)

# convert DO genotypes, but in parallel
g <- encode_geno(g, alleles)

# drop markers from gmap and pmap
gmap <- gmap[rownames(g),]
pmap <- pmap[rownames(g),]

# reorder markers
chr <- as.character(pmap$chr)
chr[chr=="X"] <- "20"
o <- order(as.numeric(as.character(chr)), pmap$pos)
gmap <- gmap[o,]
pmap <- pmap[o,]
g <- g[o,]
founder_geno <- founder_geno[o,]

# save to files
write.table(gmap, "../recla_gmap.csv", sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
write.table(pmap, "../recla_pmap.csv", sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE)

tmp <- cbind(marker=rownames(founder_geno), founder_geno)
write.table(tmp, "../recla_foundergeno.csv", sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE)

tmp <- cbind(marker=rownames(g), g)
write.table(tmp, "../recla_geno.csv", sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE)


# control file
write_control_file("../recla.json",
                   description="DO data from Recla et al. (2014) Mamm Genome 25:211-222",
                   crosstype="do",
                   geno_file="recla_geno.csv",
                   geno_transposed=TRUE,
                   founder_geno_file="recla_foundergeno.csv",
                   founder_geno_transposed=TRUE,
                   geno_codes=list("A"=1,"H"=2,"B"=3),
                   pheno_file="recla_pheno.csv",
                   pheno_transposed=FALSE,
                   covar_file="recla_covar.csv",
                   sex_covar="Sex",
                   sex_codes=list(female="female", male="male"),
                   xchr="X",
                   crossinfo_covar="ngen",
                   gmap_file="recla_gmap.csv",
                   pmap_file="recla_pmap.csv",
                   alleles=LETTERS[1:8],
                   overwrite=TRUE)
