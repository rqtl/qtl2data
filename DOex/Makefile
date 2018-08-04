R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

all: DOex.json DOex_genoprobs.rds

# this also creates a bunch of other data files, including DOex.zip
DOex.json: Code/construct_example.R
	R $(R_OPTS) -e "source('$<')"

# this also creates DOex_alleleprobs.rds
DOex_genoprobs.rds: Code/calc_genoprob.R DOex.json
	R $(R_OPTS) -e "source('$<')"
