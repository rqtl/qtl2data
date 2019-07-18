bxd.zip: bxd.json
	R -e "qtl2::zip_datafiles('bxd.json', overwrite=TRUE)"
bxd.json: R/convertData.R RawData/BXD_Geno-19Jan2017_forGN.xlsx RawData/phenotypes_export_2018-02-27.txt
	cd R;R -e "source('$(<F)')"
