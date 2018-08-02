# split vector 1:n into pieces for parallel processing
vector4parallel <-
    function(n, cores=0)
{
    if(cores==0) cores <- parallel::detectCores()

    m <- as.list(data.frame(matrix(1:(ceiling(n/cores)*cores), ncol=cores)))
    lapply(m, function(a) a[a<=n])

}
