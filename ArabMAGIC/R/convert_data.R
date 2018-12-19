# convert Gnan et al (2014) data for R/qtl2

# needed packages
library(data.table)
library(qtl2convert)
library(qtl2)


##############################
# Genotype Data
##############################

# read genotype data for magic lines
g <- data.table::fread("../RawData/magic_genotypes.csv", data.table=FALSE,
                       na.strings=c("ND", "--"))
# line names as row names
rownames(g) <- g[,1]
# transpose and drop line names columns
g <- t(g[,-1])



# read genotype data for founders
y <- data.table::fread("../RawData/founder_genotypes.csv", skip=1, data.table=FALSE)
# marker names as rownames
rownames(y) <- y[,6]



# determine alleles present at each marker
alleles <- apply(g, 1, function(a) { u <- unique(a[!is.na(a)]);
    onechar <- u[nchar(u)==1]
    twochar <- unlist(strsplit(u[nchar(u)==2], ""))
    unique(c(onechar, twochar)) })
alleles <- t(alleles) # matrix with 1260 rows and two columns



# grab just the founder genotypes
fg <- as.matrix(y[,18 + 0:18])


# swap founder alleles when column 10 is not missing
to_switch <- rownames(y)[y[,10] != ""]

# switch the founder alleles where necessary
fg_2sw <- fg[to_switch,]
fg_sw <- fg_2sw
fg_sw[fg_2sw=="A"] <- "T"
fg_sw[fg_2sw=="C"] <- "G"
fg_sw[fg_2sw=="G"] <- "C"
fg_sw[fg_2sw=="T"] <- "A"
fg[to_switch,] <- fg_sw

# omit values that aren't among the codes
for(i in 1:nrow(fg)) {
    fg[i, !(fg[i,] %in% alleles[i,])] <- NA
}


##############################
# Physical map
##############################

# Download this from Richard Mott's website, as it provides with TAIR10 genome coordinates
map <- lapply(1:5, function(i){
  # Map file address
  map_file <- paste0("http://mtweb.cs.ucl.ac.uk/mus/www/POOLING/ARABIDOPSIS/FOUNDER/GENOTYPES/chr",
                     i, ".MAGIC.map")
  # Reading data
  data.table::fread(map_file, data.table = FALSE)
})
map <- do.call("rbind", map)

# Retain only relevant columns
map <- map[, c("marker", "chromosome", "bp")]
names(map) <- c("marker", "chr", "pos")
rownames(map) <- map$marker
map$pos <- map$pos/1e6

# retain only common markers
map <- map[rownames(map) %in% rownames(g), ]
g <- g[rownames(g) %in% rownames(map), ]
fg <- fg[rownames(fg) %in% rownames(map), ]

# make sure markers are in the same order
stopifnot( all(map$marker == rownames(g)) )
stopifnot( all(map$marker == rownames(fg)) )

# convert the genotypes to A,H,B where A = most frequent allele in the founders
minor_allele <- apply(fg, 1, function(a) names(sort(table(a)))[1])

# subset the alleles
alleles <- alleles[rownames(fg),]

# swap columns in alleles so that the second column is the minor allele
sw <- which(minor_allele == alleles[,1])
alleles[sw, ] <- alleles[sw, 2:1]

# encode the genotypes
gg <- encode_geno(g, alleles)
fgg <- encode_geno(fg, alleles)

# add MAGIC. to the ids in the genotype data
colnames(gg) <- paste0("MAGIC.", colnames(gg))

# write the files
qtl2convert::write2csv(map, "../arabmagic_pmap.csv",
                       comment=paste("Arabidopsis MAGIC physical map with positions in Mbp,",
                                     "Gnan et al (2014) 10.1534/genetics.114.170746"),
                       overwrite=TRUE)

qtl2convert::write2csv(cbind(marker=rownames(gg), gg), "../arabmagic_geno.csv",
                       comment=paste("Arabidopsis MAGIC genotype data,",
                                     "Gnan et al (2014) 10.1534/genetics.114.170746"),
                       overwrite=TRUE)

qtl2convert::write2csv(cbind(marker=rownames(fgg), fgg), "../arabmagic_foundergeno.csv",
                       comment=paste("Arabidopsis MAGIC founder genotype data,",
                                     "Gnan et al (2014) 10.1534/genetics.114.170746"),
                       overwrite=TRUE)


##############################
# Read and Write the Phenotypes
##############################

# read and write the phenotype data
phe <- data.table::fread("../RawData/phenotypes.csv", data.table=FALSE)

# the line IDs in the phenotype data are crap, because the following lines were missing
# so we need to re-number the lines
missing_lines <- c( 13,  76, 281, 329, 411,    445, 536, 540, 550, 560,
                   564, 568, 569, 572, 587,    589, 596, 603, 640, 648,
                   649, 650, 663, 665, 667,    678)
lines <- paste0("MAGIC.", 1:703)[-missing_lines]
phe <- phe[seq_along(lines),]
rownames(phe) <- phe[,1] <- lines

# clean up phenotype names
cn <- colnames(phe)
cn[1] <- "id"
cn <- tolower(cn)
cn <- gsub(" ", "_", cn, fixed=TRUE)
cn <- gsub(".", "_", cn, fixed=TRUE)
cn <- gsub("__", "_", cn, fixed=TRUE)
colnames(phe) <- cn

# write the phenotypes
qtl2convert::write2csv(phe, "../arabmagic_pheno.csv",
                       comment=paste("Arabidopsis MAGIC phenotype data,",
                                     "Gnan et al (2014) 10.1534/genetics.114.170746"),
                       overwrite=TRUE)

##############################
# Write the Control File
##############################

# control file
write_control_file("../arabmagic.json",
                   crosstype="magic19",
                   geno_file="arabmagic_geno.csv",
                   founder_geno_file="arabmagic_foundergeno.csv",
                   pheno_file="arabmagic_pheno.csv",
                   pmap_file="arabmagic_pmap.csv",
                   gmap_file="arabmagic_pmap.csv", # take physical map also as genetic map
                   geno_codes=c(A=1, H=2, B=3),
                   geno_transposed=TRUE,
                   founder_geno_transposed=TRUE,
                   description="Arabidopsis MAGIC data, Gnan et al (2014) 10.1534/genetics.114.170746",
                   overwrite=TRUE)
