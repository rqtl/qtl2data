# grab X and Y chromosome SNP intensities for Gatti et al (2014) data

# download Gatti_2014.zip from https://phenome.jax.org/projects/Gatti2
# extract and grab FinalReport.txt and SNP_Map.txt

# load detailed genotypes and SNP map
int <- data.table::fread("FinalReport.txt", skip=9, data.table=FALSE)
map <- data.table::fread("SNP_Map.txt", , data.table=FALSE)

# snps on the X and Y
Xsnp <- map$Name[map$Chromosome=="X"]
Ysnp <- map$Name[map$Chromosome=="Y"]

# individual IDs
id <- unique(int$Sample.ID)

# grab just those rows
Xdat <- int[int$SNP.Name %in% Xsnp,]
Ydat <- int[int$SNP.Name %in% Ysnp,]

# grab just those rows and reorganize as matrices
#  (taking the sum of two allele intensities)
# [ this is rather slow :( ]
Xint <- matrix(ncol=length(Xsnp), nrow=length(id))
Yint <- matrix(ncol=length(Ysnp), nrow=length(id))
dimnames(Xint) <- list(id, Xsnp)
dimnames(Yint) <- list(id, Ysnp)
for(i in id) { # some arrays done twice; some three times
     tmp <- Xdat[Xdat$SNP.Name %in% Xsnp & Xdat$Sample.ID==i,] # grab rows
    tmp <- tapply(tmp$X, tmp$SNP.Name, mean) + tapply(tmp$Y, tmp$SNP.Name, mean) # take averages then sum
    Xint[i,Xsnp] <- tmp[Xsnp] # fill in

    tmp <- Ydat[Ydat$SNP.Name %in% Ysnp & Ydat$Sample.ID==i,] # grab rows
    tmp <- tapply(tmp$X, tmp$SNP.Name, mean) + tapply(tmp$Y, tmp$SNP.Name, mean) # take averages then sum
    Yint[i,Ysnp] <- tmp[Ysnp] # fill in
}

write.table(data.frame(id=rownames(Xint), Xint, stringsAsFactors=FALSE),
            "Gatti2014_Xint.csv", sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(data.frame(id=rownames(Yint), Yint, stringsAsFactors=FALSE),
            "Gatti2014_Yint.csv", sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
