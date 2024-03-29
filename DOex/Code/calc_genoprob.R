# calc genoprobs for DOex data

library(qtl2)
DOex <- read_cross2("DOex.zip")

pr <- calc_genoprob(DOex, error_prob=0.002, map_function="c-f", cores=0)
saveRDS(pr, "DOex_genoprobs.rds")

pr2 <- pr[,"2"]
saveRDS(pr2, "DOex_genoprobs_2.rds")

apr <- genoprob_to_alleleprob(pr, cores=0)
saveRDS(apr, "DOex_alleleprobs.rds")
