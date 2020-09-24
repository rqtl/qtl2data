R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

maize_magic.json: R/prep_data.R Rawdata/foundergeno_imputed.csv
	R CMD BATCH $(R_OPTS) $<

Rawdata/foundergeno_imputed.csv: R/combine_foundergeno.R Rawdata/MMfounders.geno.gz
	R CMD BATCH $(R_OPTS) $<

Rawdata/MMfounders.geno.gz: R/download_data.R
	R CMD BATCH $(R_OPTS) $<
