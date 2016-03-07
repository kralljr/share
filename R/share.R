#' share
#'
#' \code{share} estimates major sources of pollution using the share method
#'
#' This function estimates major sources of pollution across
#' multiple ambient monitors.  Other outputs of share include
#' a matrix to guide pooling short-term health effects of pollution
#' to conduct regional and national studies.
#' Works on list where each element of the list is a 
#' dataframe corresponding to one ambient monitor.  
#' The first column of each dataframe is date and all subsequent 
#' columns are concentrations of chemical constituents.
#' 
#' @param data data frame of daily constituent concentrations with date as first column
#' @param cut cutoff for eigenvalues (see nmsource), default is 1. 
#' @param nmsources number of major sources.  If null, uses number of eigenvalues of the correlation matrix greater than cut.
#' @param thres cutoff for angle between local and major sources
#' @param ndays number of days of data necessary to apply PCA for each monitor
#' @export
#' @examples
#' data(consConc)
#' share(consConc)
share <- function(data, cut = 1, nmsources = NULL, thres = pi/4, ndays = 50) {
    
    #Apply vPCA to each monitor and across all monitors
    vPCA <- outervPCA(data, ndays, cut, nmsources)	
    major.sig <- vPCA$major.sig
    source.sig <- vPCA$source.sig
    
    #Match local and regional source signatures
    match1 <- matchfun(source.sig, major.sig, thres = thres)
    ang <- match1$angle
    match <- match1$match
    
    
    reorder1 <- sapply(source.sig, function(x) {
        matchfun(list(major.sig), x, thres = thres)[[1]]
    })
    
    
    #find sources at each site
    share1 <- lapply(match, whichCS)
    shareREORG <- lapply(reorder1, whichCS)
    
    #find sources across all sites
    Sources <- unique(unlist(share1))
    
    share <- list()
    
    share$major.sig <- major.sig
    share$source.sig <- source.sig
    share$Sources <- Sources
    share$share <- share1
    share$angle <- angle
    share$match <- match
    share$shareREORG <- shareREORG
    
    class(share) <- "share"
    share
}
