# Prepare a smallish set of DO example data

# Full Recla data
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Recla/recla.zip")
DOex <- read_cross2(file)

# reduce to 3 chromosomes
DOex <- DOex[, c("2","3","X")]

# just one phenotype
DOex$pheno <- DOex$pheno[,5,drop=FALSE]

# subset map to 0.5 cM
sub <- reduce_markers(DOex$gmap, 0.5)
DOex <- pull_markers(DOex, unlist(lapply(sub, names)))

# write files
### genotypes
g <- DOex$geno[[1]]
for(i in 2:length(DOex$geno))
    g <- cbind(g, DOex$geno[[i]])
g <- as.data.frame(g)
for(i in 1:ncol(g))
    g[,i] <- sub("1","A",sub("2","H",sub("3","B",sub("0","-", as.character(g[,i])))))
write.table(cbind(ind=rownames(g),g), file="DOex_geno.csv",
            sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

### founder genotypes
fg <- DOex$founder_geno[[1]]
for(i in 2:length(DOex$founder_geno))
    fg <- cbind(fg, DOex$founder_geno[[i]])
fg <- as.data.frame(fg)
for(i in 1:ncol(g))
    fg[,i] <- sub("1","A",sub("2","H",sub("3","B",sub("0","-", as.character(fg[,i])))))
fg <- t(fg)
write.table(cbind(marker=rownames(fg),fg), file="DOex_foundergeno.csv",
            sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

### phenotypes
write.table(cbind(ind=rownames(DOex$pheno), DOex$pheno), file="DOex_pheno.csv",
            sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

### covariates
write.table(cbind(ind=rownames(DOex$covar), DOex$covar), file="DOex_covar.csv",
            sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

### genetic map
gmaptab <- qtl::map2table(DOex$gmap)
write.table(cbind(marker=rownames(gmaptab), gmaptab), file="DOex_gmap.csv",
            sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

### physical map
pmaptab <- qtl::map2table(DOex$pmap)
write.table(cbind(marker=rownames(pmaptab), pmaptab), file="DOex_pmap.csv",
            sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

write_control_file(output_file="DOex.json",
                   crosstype="do",
                   geno_file="DOex_geno.csv",
                   geno_codes=c(A=1, H=2, B=3),
                   founder_geno_file="DOex_foundergeno.csv",
                   founder_geno_transposed=TRUE,
                   gmap_file="DOex_gmap.csv",
                   pmap_file="DOex_pmap.csv",
                   pheno_file="DOex_pheno.csv",
                   covar_file="DOex_covar.csv",
                   sex_covar="Sex",
                   sex_codes=list(male="male", female="female"),
                   xchr="X",
                   crossinfo_covar="ngen",
                   alleles=LETTERS[1:8],
                   description=paste0("Example DO data, from Recla et al. (2014) Mamm ",
                                      "Genome 25:211-222; subsetted to 3 chr, 1 phenotype, ",
                                      "and reduced markers"),
                   overwrite=TRUE)

zip_datafiles("DOex.json", overwrite=TRUE)
