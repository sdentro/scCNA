getStartsEnds.badbins <-
function(window,
                                   chr,
                                   pathWindows=paste0("/srv/data/vanloo/mtarabichi/sc-all/code/variable_",
                                               window,
                                               "_150_bowtie"),
                                   pathBadBins=paste0("/srv/data/vanloo/mtarabichi/sc-all/code/badbins_variable_",
                                                      window,
                                                      "_150_bowtie"))
{
    badbins <- read.table(pathBadBins)[,1]
    t <- read.table(pathWindows,header=T)
    t <- t[-c(badbins),]
    subt <- t[t$CHR==chr,]
    starts <- c(1,as.numeric(as.character(subt$END[-c(length(subt$END))]))+1)
    ends <- as.numeric(as.character(subt$END))
    list(starts=starts,ends=ends)
}
