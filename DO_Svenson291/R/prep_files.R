# create other data files, including the control file

library(qtl2)
chr <- c(1:19, "X")
write_control_file("svenson.json",
                   crosstype="do",
                   description="MegaMUGA genotypes for 291 DO mice from Karen Svenson and colleages",
                   founder_geno_file=paste0("MM/MM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("MM/MM_gmap", chr, ".csv"),
                   pmap_file=paste0("MM/MM_pmap", chr, ".csv"),
                   geno_file=paste0("svenson_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   covar_file="svenson_covar.csv",
                   sex_covar="sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="ngen",
                   overwrite=TRUE)

# IDs 276-325 are from 1st file (99 samples); generation 8?
# IDs 326-425 are from 2nd file (192 samples); generation 11?

ids <- colnames(read.table("svenson_geno19.csv", comment.char="#",
                           header=TRUE, nrows=1, sep=","))[-1]
sex <- substr(ids, 1, 1)
num <- as.numeric(substr(ids, 2, nchar(ids)))
ngen <- rep(8, length(num))
ngen[num >= 326] <- 11
covar <- cbind(id=ids, sex=sex, ngen=ngen)
write.table(covar, "svenson_covar.csv", sep=",",
            row.names=FALSE, col.names=TRUE, quote=FALSE)
