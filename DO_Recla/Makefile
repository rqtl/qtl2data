all: recla.json recla.zip

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

Orig/Recla2013_pheno.csv:
	[ -d Orig ] || mkdir Orig
	curl https://phenomedoc.jax.org/QTL_Archive/recla_2013/$(@F) > $@

Orig/Recla2013_geno.csv:
	[ -d Orig ] || mkdir Orig
	curl https://phenomedoc.jax.org/QTL_Archive/recla_2013/$(@F) > $@

Orig/Chesler4_table.csv:
	[ -d Orig ] || mkdir Orig
	curl http://phenomedoc.jax.org/MPD_projdatasets/Chesler4.csv > $@

Orig/muga_code.Rdata:
	[ -d Orig ] || mkdir Orig
	curl ftp://ftp.jax.org/MUGA/$(@F) > $@

Orig/muga_geno.Rdata:
	[ -d Orig ] || mkdir Orig
	curl ftp://ftp.jax.org/MUGA/$(@F) > $@

Orig/muga_snps.Rdata:
	[ -d Orig ] || mkdir Orig
	curl ftp://ftp.jax.org/MUGA/$(@F) > $@

recla.json: Code/convert_data.R Orig/Recla2013_geno.csv Orig/Recla2013_pheno.csv Orig/muga_code.Rdata Orig/muga_geno.Rdata Orig/muga_snps.Rdata Orig/Chesler4_table.csv
	cd Code;R $(R_OPTS) -e 'source("$(<F)", echo=TRUE)'

recla.zip: recla.json
	zip recla.zip recla.json recla_*.csv
