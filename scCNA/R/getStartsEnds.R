getStartsEnds <-
function(window,
                          chr,
                          path=paste0("/srv/data/vanloo/mtarabichi/sc-all/code/variable_",
                                      window,
                                      "_150_bowtie"))
{
    t <- read.table(path,header=T)
    subt <- t[t$CHR==chr,]
    starts <- c(1,as.numeric(as.character(subt$END[-c(length(subt$END))]))+1)
    ends <- as.numeric(as.character(subt$END))
    list(starts=starts,ends=ends)
}
