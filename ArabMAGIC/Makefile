R_OPTS=--no-save --no-restore --no-init-file --no-site-file

arabmagic.zip: arabmagic.json
	R -e "qtl2::zip_datafiles('$<', overwrite=TRUE)"

arabmagic.json: R/convert_data.R \
				RawData/phenotypes.csv \
				RawData/magic_genotypes.csv \
				RawData/founder_genotypes.csv
	cd R;R $(R_OPTS) -e "source('$(<F)')"

RawData/phenotypes.csv: RawData/genetics.114.170746-3_phenotypes.xlsx Python/xlsx2csv.py
	Python/xlsx2csv.py $< --sheet 1 > $@

RawData/magic_genotypes.csv: RawData/genetics.114.170746-6_genotypes.xlsx Python/xlsx2csv.py
	Python/xlsx2csv.py $< --sheet 1 > $@

RawData/founder_genotypes.csv: RawData/genetics.114.170746-6_genotypes.xlsx Python/xlsx2csv.py
	Python/xlsx2csv.py $< --sheet 2 > tmp
	head -n1262 tmp > $@
	rm tmp

RawData/genetics.114.170746-3_phenotypes.xlsx:
	[ -d RawData ] || mkdir RawData
	wget http://www.genetics.org/lookup/suppl/doi:10.1534/genetics.114.170746/-/DC1/genetics.114.170746-3.xls
	mv genetics.114.170746-3.xls $@

RawData/genetics.114.170746-6_genotypes.xlsx:
	[ -d RawData ] || mkdir RawData
	wget http://www.genetics.org/lookup/suppl/doi:10.1534/genetics.114.170746/-/DC1/genetics.114.170746-6.xls
	mv genetics.114.170746-6.xls $@

Python/xlsx2csv.py:
	[ -d Python ] || mkdir Python
	wget https://raw.githubusercontent.com/dilshod/xlsx2csv/master/xlsx2csv.py
	mv xlsx2csv.py Python/
	chmod +x $@
