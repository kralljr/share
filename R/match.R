
#' \code{matchfun} Matches major and local source signatures
#'
#' @param major.sig major source signatures 
#' @param source.sig local source signatures from all monitors
#' @param thres cutoff for angle between local and major sources
#' @export
matchfun <- function(source.sig, major.sig, thres = pi/4) {
    
    
    match1 <- list()
    ang <- list()
    
    #For each monitor
    for(i in 1 : length(source.sig)) {
        
        #Get angle
        ang1 <- angle(major.sig, source.sig[[i]])
        
        #Add columns if wrong shape
        nc <- ncol(source.sig[[i]])
        ns <- nc - ncol(major.sig) 
        if(ns > 0) {
            adds <- matrix(rep(100, nc * ns), nrow = nc, ncol = ns)
            ang1 <- cbind(ang1, adds)
        }
        
        
        #Make large for exceed threshold
        whT <- which(ang1 > thres, arr.ind = T)
        if(length(whT) > 0) {
            for(k in 1 : nrow(whT)) {
                ang1[whT[k, 1], whT[k, 2]] <- 100
            }
        }
        
        #Match
        mins <- as.vector(solve_LSAP(ang1))
        
        
        xs <- seq(1, nc)
        mins <- cbind(xs, mins)
        sq <- max(c(nc, ncol(major.sig)))
        match <- matrix(0, nrow = sq, ncol = sq)
        for(j in 1 : nrow(mins)) {
            #if not NA and angle less than 50 assign match
            if(!is.na(mins[j, 2])) {
                if(ang1[j, mins[j, 2]] <  50) {
                    match[mins[j, 1], mins[j, 2]] <- 1
                }
            }
        }
        match <- match[1 : nc, 1 : ncol(major.sig)]
        
        
        ang[[i]] <- ang1
        match1[[i]] <- match
    }
    
    out <- list()
    out$match <- match1
    out$angle <- ang
    out
}
#' \code{angle} Computes angle between major and local source signatures
#'
#' @param major.sig major source signatures 
#' @param source.sig1 local source signature from one monitor
#' @export
#' @examples
#' data(simdat)
#' share1 <- share(simdat)
#' angle(share1$major.sig, share1$source.sig)
angle <- function(major.sig, source.sig1) {
    
    cp <- abs(t(source.sig1) %*% major.sig)
    #compute norms
    l2 <- sqrt(rowSums(t(source.sig1)^2)) %*% 
        t(sqrt(colSums(major.sig^2)))
    
    #get angles    
    temp <- cp / l2
    ang1 <- acos(temp)
    
    ang1
}

