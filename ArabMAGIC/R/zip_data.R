# zip the Arabidopsis MAGIC data, in two version (TAIR8 and TAIR9)
# - first subset the _geno and _foundergeno data to the markers with positions in the corresponding map build

library(qtl2)

set.seed(20190215)
tmpdir <- file.path(tempdir(), paste(sample(letters, 20, replace=TRUE), collapse=""))
tmpdir <- sub("//", "/", tmpdir)
if(!dir.exists(tmpdir)) dir.create(tmpdir)

for(build in c(8, 9)) {
    json_file <- paste0("arabmagic_tair", build, ".json")
    file.copy(file.path("..", json_file), tmpdir)
    file.copy("../arabmagic_pheno.csv", tmpdir)
    pmap_file <- paste0("../arabmagic_pmap_tair", build, ".csv")
    file.copy(pmap_file, tmpdir)

    markers <- rownames(read_csv(pmap_file))

    for(g in c("geno", "foundergeno")) {
        file <- paste0("../arabmagic_", g, ".csv")
        desc <- sub("^# ", "", readLines(file, n=1))
        dat <- read_csv(file)
        dat <- dat[markers,]
        qtl2convert::write2csv(cbind(marker=rownames(dat), dat),
                               file.path(tmpdir, basename(file)),
                               comment=desc, overwrite=TRUE)

    }

    zip_datafiles(file.path(tmpdir, json_file))
    file.copy(file.path(tmpdir, sub(".json", ".zip", json_file)), "..")
}

unlink(file.path(tmpdir, list.files(tmpdir)))
