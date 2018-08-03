R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

DOex.json: Code/construct_example.R
	R $(R_OPTS) -e "source('$<')"
