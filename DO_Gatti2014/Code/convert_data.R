## convert DO data from Gatti et al (2014)
## for use with R/qtl2
#
# DO data from https://phenome.jax.org/projects/Gatti2
# founder data from ftp://ftp.jax.org/MUGA/

library(data.table)
library(qtl2)
library(qtl2convert)


# load phenotype file
phe <- data.table::fread("../Orig/Gatti_2014_pheno.csv",
                         data.table=FALSE, colClasses="character")
covar <- phe[,c(1:3)]
phe <- phe[,c(1,4:5)]
colnames(covar)[1] <- colnames(phe)[1] <- "id"

# load sample info (including generation number)
gen <- data.table::fread("../Orig/jaxwest_samples.csv",
                         colClasses="character", data.table=FALSE)
gen[,1] <- sub("\\-", ".", gen[,1])
m <- match(phe[,1], gen[,1])
stopifnot(all(!is.na(m)))
covar$ngen <- gen[m, "DO Generation"]

# birth dates
bd <- gen[m, "Birth Date"]
bd <- vapply(strsplit(bd, "-"), function(a) paste0("20", a[3], "-", a[1], "-", a[2]), "")
covar$birth_date <- bd

# write phenotypes
cat(paste0("# phenotype data from Gatti et al. (2014) G3 4:1623-1633\n",
           "# nrow ", nrow(phe), "\n",
           "# ncol ", ncol(phe), "\n"), file="../do_pheno.csv")
write.table(phe, file="../do_pheno.csv", quote=FALSE, sep=",",
            row.names=FALSE, col.names=TRUE, append=TRUE)

# write covariates
cat(paste0("# covariate data from Gatti et al. (2014) G3 4:1623-1633\n",
           "# nrow ", nrow(covar), "\n",
           "# ncol ", ncol(covar), "\n"), file="../do_covar.csv")
write.table(covar, file="../do_covar.csv", quote=FALSE, sep=",",
            row.names=FALSE, col.names=TRUE, append=TRUE)

# genotypes
g <- data.table::fread("../Orig/Gatti_2014_geno.csv", data.table=FALSE, colClasses="character")
rownames(g) <- g[,1]
g <- g[,-1]
g[g=="--"] <- NA
ng <- apply(g, 1, function(a) length(unique(a[!is.na(a)]))) # number of unique genotypes
g <- g[ng == 3,]

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

# I've dropped some markers from the DO data; otherwise they're all there
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

# function to infer founder genotype
pat <- paste0(LETTERS, LETTERS)[1:8]
for(i in 1:8) {
    founder_geno[,i] <- find_consensus_geno(muga_geno[,muga_code==pat[i]])
}

# count unique values in each row
nug <- apply(founder_geno, 1, function(a) length(unique(a[!is.na(a)])))

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
chr[chr=="X"] <- 20
o <- order(as.numeric(as.character(chr)), pmap$pos)
gmap <- gmap[o,]
pmap <- pmap[o,]
g <- g[o,]
founder_geno <- founder_geno[o,]

# save to files
cat(paste0("# genetic marker map from Gatti et al. (2014) G3 4:1623-1633\n",
           "# nrow ", nrow(gmap), "\n",
           "# ncol ", ncol(gmap), "\n"), file="../do_gmap.csv")
write.table(gmap, "../do_gmap.csv", sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE, append=TRUE)
cat(paste0("# physical marker map from Gatti et al. (2014) G3 4:1623-1633\n",
           "# nrow ", nrow(pmap), "\n",
           "# ncol ", ncol(pmap), "\n"), file="../do_pmap.csv")
write.table(pmap, "../do_pmap.csv", sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE, append=TRUE)

tmp <- cbind(marker=rownames(founder_geno), founder_geno)
cat(paste0("# founder genotype data from Gatti et al. (2014) G3 4:1623-1633\n",
           "# nrow ", nrow(tmp), "\n",
           "# ncol ", ncol(tmp), "\n"), file="../do_foundergeno.csv")
write.table(tmp, "../do_foundergeno.csv", sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE, append=TRUE)

tmp <- cbind(marker=rownames(g), g)
cat(paste0("# genotype data from Gatti et al. (2014) G3 4:1623-1633\n",
           "# nrow ", nrow(tmp), "\n",
           "# ncol ", ncol(tmp), "\n"), file="../do_geno.csv")
write.table(tmp, "../do_geno.csv", sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE, append=TRUE)


# control file
write_control_file("../do.json",
                   description="DO data from Gatti et al. (2014) G3 4:1623-1633",
                   crosstype="do",
                   geno_file="do_geno.csv",
                   geno_transposed=TRUE,
                   founder_geno="do_foundergeno.csv",
                   founder_geno_transposed=TRUE,
                   geno_codes=list(A="1",H="2",B="3"),
                   pheno_file="do_pheno.csv",
                   pheno_transposed=FALSE,
                   covar_file="do_covar.csv",
                   sex_covar="Sex",
                   sex_codes=c(F="female", M="male"),
                   xchr="X",
                   crossinfo_covar="ngen",
                   gmap_file="do_gmap.csv",
                   pmap_file="do_pmap.csv",
                   overwrite=TRUE)
