# calculate genotype probabilities for the CC

library(qtl2)

cc_rds <- "cc.rds"
if(file.exists(cc_rds)) {
    cc <- readRDS(cc_rds)
} else {
    cc <- read_cross2("../cc.zip", quiet=FALSE)
    saveRDS(cc, cc_rds)
}

gmap <- insert_pseudomarkers(cc$gmap, step=0.2, stepwidth="max")
pmap <- interp_map(gmap, cc$gmap, cc$pmap)
saveRDS(pmap, "cc_probs_pmap.rds")

pr <- calc_genoprob(cc, gmap, err=0.002, map_function="c-f", cores=0)
saveRDS(pr, "cc_probs.rds")
