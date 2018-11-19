R_OPTS=--no-save --no-restore --no-init-file --no-site-file

svenson.zip: R/zip_data.R svenson_geno19.csv svenson.json
	R $(R_OPTS) -e "source('$<')"

svenson.json: R/prep_files.R svenson_geno19.csv
	R $(R_OPTS) -e "source('$<')"
# also creates svenson_covar.csv

svenson_geno19.csv: R/geneseek2qtl2.R \
					RawData/Jackson_Lab_Mouse_30jan2013_FinalReport.txt.gz \
					RawData/Jackson_Lab-Ciciotte_Mouse_05jun2013_FinalReport.txt.gz \
					MM/MM_allelecodes.csv
	R $(R_OPTS) -e "source('$<')"
# also creates svenson_geno*.csv for other chr
# plus svenson_chrXint.csv and Data/svenson_chrYint.csv
