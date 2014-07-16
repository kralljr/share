#' Find regional sources
#'
#' \code{outervPCA} Applies vPCA locally and across all monitors
#'
#' Applies vPCA to data from each monitor, and then
#' a second vPCA to the PCs from the first stage to find regional sources
#'
#' @param data data frame of daily constituent concentrations with date as first column
#' @param ndays number of days of data necessary to apply PCA
#' @param nmsource number of major sources.  If null, uses number of eigenvalues of the correlation matrix greater than cut.
#' @param cut cutoff for eigenvalues (see nmsource), default is 1.
#' #param ... other arguments
outervPCA <- function(data, ndays = 50, cut = 1, nmsource = NULL) {
    
    N <- length(data)
    
    #Setup
    source.sig <- vector(mode = "list", length = N)
    dates <- vector(mode = "list", length = N)
    nums <- 0
    scores <- list()
    
    #For each monitor
    for ( i in 1 : N) {    
        
        #Set data
        dat1 <- data[[i]]
        if(tolower(colnames(dat1)[1]) != "date") {
            stop("First column is not date")
        }
        dates <- dat1[, 1]
        dat1 <- dat1[, -1]
        
        #Get complete cases
        dates[[i]] <- dates[complete.cases(dat1)]
        dat1 <- dat1[complete.cases(dat1), ]
        nc <- ncol(dat1)
        
        #Create concatenated matrix
        if(i == 1) {
            datall <- vector(, length = nc)
            V <- matrix(nrow = nc, ncol = 1)
        }
        datall <- rbind(datall, dat1)
        
        
        #Perform if enough data
        if(nrow(dat1) >= ndays) {
            
            #get varimax and unrotated loadings
            vpc <- vPCAf(data = dat1, cut = cut, ...)
            
            #Set up outcome
            source.sig[[i]] <- vpc$vpc
            V <- cbind(V, vpc$pc)
            scores[[i]] <- vpc$scores
            
        }else{
            nums <- c(nums, i)
        }
    }#end loop over sites
    
    V <- V[, -1]
    
    #Get major sources
    datall <- datall[-1, ]
    mvpc <- vPCAf(t(V), nsources = nmsource, cut = cut)
    mvpc <- mvpc$vpc
    
    #eliminate sites with not enough data
    if(length(nums) > 1) {
        nums <- nums[-1]
        source.sig <- source.sig[-nums]
    }
    
    reg$source.sig <- source.sig
    reg$mvpc <- mvpc
    reg$V <- V
    reg$datall <- datall
    reg$dates <- dates
    reg$scores <- scores
    
    reg
}

#' \code{vPCA} applies PCA with a varimax rotation
#'
#' Applies PCA with a varimax rotation to data from a single monitor
#'
#' @param data data frame of daily constituent concentrations 
#' @param nsource number of sources.  If null, uses number of eigenvalues of the correlation matrix greater than cut.
#' @param cut cutoff for eigenvalues (see nsource), default is 1.
#' #param ... other arguments
vPCA <- function(data, nsource = NULL, cut = 1, ...) {
    
    #Perform PCA
    pc1 <- prcomp(data, retx = T, ...)
    
    #Scale eigenvectors by sdev
    rots <- sweep(pc1$rot, 2, pc1$sdev, "*")
    
    #Specify number of sources
    if(is.null(nsource)) {
        nsource <- length(which(pc1$sdev > cut))
    }
    
    #Apply varimax rotation
    vpc <- varimax(rots[, 1 : nsource], normalize = T)
    
    vpc$vpc <- vpc
    vpc$nsource <- nsource
    vpc$pc <- pc1
    
    class(vpc) <- "vpca"
    vpc

}



#' \code{vPCAf} Applies PCA and vPCA and formats outcome
#'
#' @param data data frame of daily constituent concentrations 
#' @param nsource number of sources.  If null, uses number of eigenvalues of the correlation matrix greater than cut.
#' @param cut cutoff for eigenvalues (see nsource), default is 1.
vPCAf <- function(data, nsources = NULL, cut = 1) {
    
    #Number of constituents
    nc <- ncol(data)
    
    #Standardize data
    temp <- stdize1(data)
    dat1 <- temp[[1]]
    wh0 <- temp[[2]]
    
    #Apply vPCA, PCA
    vpc <- vPCA(data = dat1, cut = cut, nsources = nsources)
    pc <- vpc$pc
    nsource <- vpc$nsource
    scores <- pc$x[, 1 : nsource]
    pc <- pc$rot[1: nc, 1 : nsource]
    
    vpc <- pr1$vpc
    vpc <- vpc$load[1: nc, 1 : nsource]
    
    vpc1 <- matrix(rep(NA, nc * nsource), nrow = nc, ncol = nsource)
    pc1 <- vpc1
    if(length(wh0 ) > 0) {
        vpc1[-wh0, ] <- vpc
        pc1[-wh0, ] <- pc
    }else{
        vpc1 <- vpc
        pc1 <- pc
    }
    
    vpc$vpc <- vpc1
    vpc$pc <- pc1
    vpc$nsource <- nsource
    vpc$scores <- scores
    vpc$data <- dat1
    
    vpc
}








