# convert genotype probabilities from Srivastava et al (2017) as R/qtl2 probs object + map
#
# supplemental data for Srivastava et al. (2017) Genomes of the Mouse
# Collaborative Cross. Genetics 206:537-556, doi:10.1534/genetics.116.198838
# available at Zenodo, doi:10.5281/zenodo.377036

library(data.table)
library(qtl2)
library(qtl2convert)

prob_dir <- "../RawData/Prob36"
files <- list.files(prob_dir, pattern=".csv$")
strains <- sub("b38V01.csv$", "", files)

message("Reading probabilities")
v <- lapply(files, function(file) data.table::fread(file.path(prob_dir, file),
                                                    data.table=FALSE))

# grab map
map <- v[[1]][,1:3]
map[,3] <- map[,3]/1e6
map <- map[map$chromosome %in% c(1:19,"X"),]
pmap <- map_df_to_list(map, chr_column="chromosome", pos_column="position(B38)")

probs <- vector("list", 20)
names(probs) <- c(1:19,"X")
message("Reorganizing probabilities")
for(chr in names(probs)) {
    probs[[chr]] <- array(dim=c(length(v), ncol(v[[1]])-3, length(pmap[[chr]])))
    dimnames(probs[[chr]]) <- list(strains, colnames(v[[1]])[-(1:3)], names(pmap[[chr]]))

    for(i in seq_along(v)) {
        probs[[chr]][i,,] <- t(v[[i]][v[[i]][,2]==chr, -(1:3)])
    }
}

# drop all but the first 8 genotypes
# force to sum to 1 with no missing values
message("Reducing to 8 states")
probs8 <- probs
for(chr in names(probs)) {
    probs8[[chr]] <- probs[[chr]][, 1:8, ] + 1e-8
    for(i in 1:nrow(probs[[chr]]))
        probs8[[chr]][i,,] <- t(t(probs8[[chr]][i,,]) / colSums(probs8[[chr]][i,,]))
}

attr(probs8, "crosstype") <- "risib8"
attr(probs8, "is_x_chr") <- setNames(rep(c(FALSE,TRUE), c(19,1)), c(1:19,"X"))
attr(probs8, "alleles") <- LETTERS[1:8]
attr(probs8, "alleleprobs") <- FALSE
class(probs8) <- c("calc_genoprob", "list")

message("Saving to files")
saveRDS(probs8, "cc_rawprobs.rds")
saveRDS(pmap, "cc_rawprobs_pmap.rds")
