# download and convert the Collaborative Cross data
#
# supplemental data for Srivastava et al. (2017) Genomes of the Mouse
# Collaborative Cross. Genetics 206:537-556, doi:10.1534/genetics.116.198838
#
# available at Zenodo, doi:10.5281/zenodo.377036


# create RawData/ and Data/ if they're not available
rawdata_dir = "../RawData"
data_dir = "../Data"
if(!dir.exists(rawdata_dir)) {
    dir.create(rawdata_dir)
}
if(!dir.exists(data_dir)) {
    dir.create(data_dir)
}

# download genotypes.zip if not available
zenodo_url = "https://zenodo.org/record/377036/files/"
zenodo_postpend = "?download=1"
genotype_file = "genotypes.zip"
distant_file = paste0(zenodo_url, genotype_file, zenodo_postpend)
local_file = file.path(rawdata_dir, genotype_file)
if(!file.exists(local_file)) {
    download.file(distant_file, local_file)
}

# create directory for data if it doesn't exist
genotype_dir = file.path(rawdata_dir, "genotypes")
if(!dir.exists(genotype_dir)) {
    dir.create(genotype_dir)
}

# unzip if it hasn't been unzipped
file = file.path(genotype_dir, "SEQgenotypes.csv")
if(!file.exists(file)) {
    unzip(local_file, exdir=genotype_dir)
}


# download Prob36.zip if not available
prob_file = "Prob36.zip"
distant_file = paste0(zenodo_url, prob_file, zenodo_postpend)
local_file = file.path(rawdata_dir, prob_file)
if(!file.exists(local_file)) {
    download.file(distant_file, local_file)
}

# create directory for data if it doesn't exist
prob_dir = file.path(rawdata_dir, "Prob36")
if(!dir.exists(prob_dir)) {
    dir.create(prob_dir)
}

# unzip if it hasn't been unzipped
gzfile = file.path(prob_dir, "CC001-Uncb38V01.csv.gz")
csvfile = file.path(prob_dir, "CC001-Uncb38V01.csv")
if(!file.exists(gzfile) && !file.exists(csvfile)) {
    unzip(local_file, exdir=prob_dir)
}

# gunzip all of the Prob36 files
if(!file.exists(csvfile)) {
    files <- list.files(prob_dir, pattern=".csv.gz$")

    for(file in files) {
        system(paste("gunzip", file.path(prob_dir, file)))
    }
}

# load the Prob36 files and determine X, Y, and M genotypes
files <- list.files(prob_dir, pattern=".csv$")
strains <- sub("\\.csv$", "", files)
probs <- setNames(vector("list", length(strains)), strains)
for(i in seq_along(files)) {
    probs[[i]] <- data.table::fread(file.path(prob_dir, files[i]), data.table=FALSE)
}

# guess the rest of the cross order

# get the founder genotypes for MMnGM markers
