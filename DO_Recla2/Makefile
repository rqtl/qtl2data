R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

recla2.zip: Code/prepData.R Orig/bult_recla_megamuga.zip
	R CMD BATCH $(R_OPTS) $<
