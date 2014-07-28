#' Basic functions for share
#' 
#' \code{ppart} Gets positive part of a vector
#'
#' @param vec numeric data vector
ppart <- function(vec) {
    wh0 <- which(vec < 0)
    if(length(wh0) > 0) {
        vec[wh0] <- 0
    }
    vec
}

#' \code{whichCS} Find matches from matchfun
#'
#' @param locreg matrix of nsources by nmsources of matches
#' @export
whichCS <- function(locreg) {
    lr <- suppressWarnings(apply(locreg, 1, function(x) min(which(x > 0))))
    lr
    
}
