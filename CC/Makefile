R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

cc.json: R/convert_cc_data.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"
