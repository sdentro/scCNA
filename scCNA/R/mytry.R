mytry <-
function(x,retVal=NA,...)
{
    out <- try(x,silent=T,...)
    if(inherits(out,"try-error")) return(retVal)
    out
}
