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
outervPCA <- function(data, ndays = 50, cut = 1, nmsource = NULL) {
    
    N <- length(data)
    
    #Setup
    source.sig <- vector(mode = "list", length = N)
    dates <- vector(mode = "list", length = N)
    nums <- 0
    scores <- list(length = N)
    
    #For each monitor
    for ( i in 1 : N) {    
        #Set data
        dat1 <- data[[i]]
        if(class(dat1[, 1]) != "Date") {
            stop("First column is not date")
        }
        
        dates1 <- dat1[, 1]
        dat1 <- dat1[, -1]
        
        #Get complete cases
        dates[[i]] <- dates1[complete.cases(dat1)]
        dat1 <- dat1[complete.cases(dat1), ]
        nc <- ncol(dat1)
        
        #Create concatenated matrix
        if(i == 1) {
            V <- matrix(nrow = nc, ncol = 1)
        }
        
        
        #Perform if enough data
        if(nrow(dat1) >= ndays) {
            
            #get varimax and unrotated loadings
            vpc <- vPCAf(data = dat1, i = i, cut = cut)
            
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
    major.sig <- vPCAf(t(V), "major", nsources = nmsource, cut = cut)
    major.sig <- major.sig$vpc
    
    #eliminate sites with not enough data
    if(length(nums) > 1) {
        nums <- nums[-1]
        source.sig <- source.sig[-nums]
    }
    
    reg <- list()
    
    reg$source.sig <- source.sig
    reg$major.sig <- major.sig
    reg$V <- V
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
vPCA <- function(data, nsources = NULL, cut = 1) {
    
    #Perform PCA
    pc1 <- prcomp(data, retx = T)
    
    #Scale eigenvectors by sdev
    rots <- sweep(pc1$rot, 2, pc1$sdev, "*")
    
    #Specify number of sources
    if(is.null(nsources)) {
        nsources <- length(which(pc1$sdev > cut))
    }
    
    #Apply varimax rotation
    vpc1 <- varimax(rots[, 1 : nsources], normalize = T)
    
    
    vpc <- list()
    vpc$vpc <- vpc1
    vpc$nsources <- nsources
    vpc$pc <- pc1
    
    class(vpc) <- "vpca"
    vpc

}



#' \code{vPCAf} Applies PCA and vPCA and formats outcome
#'
#' @param data data frame of daily constituent concentrations 
#' @param nsource number of sources.  If null, uses number of eigenvalues of the correlation matrix greater than cut.
#' @param cut cutoff for eigenvalues (see nsource), default is 1.
vPCAf <- function(data, i, nsources = NULL, cut = 1) {
    
    #Standardize data
    temp <- stdize1(data, i)
    dat1 <- temp[[1]]
    wh0 <- temp[[2]]
    
    #Number of constituents
    nc <- ncol(dat1)
    nc1 <- ncol(data)
    
    #Apply vPCA, PCA
    vpc <- vPCA(data = dat1, cut = cut, nsources = nsources)
    pc <- vpc$pc
    nsources <- vpc$nsources
    scores <- pc$x[, 1 : nsources]
    pc <- pc$rot[1: nc, 1 : nsources]
    
    vpc <- vpc$vpc
    vpc <- vpc$load[1: nc, 1 : nsources]
    
    vpc1 <- matrix(rep(0, nc1 * nsources), nrow = nc1, ncol = nsources)
    pc1 <- vpc1
    if(length(wh0 ) > 0) {
        vpc1[-wh0, ] <- vpc
        pc1[-wh0, ] <- pc
    }else{
        vpc1 <- vpc
        pc1 <- pc
    }
    
    vpc <- list()
    vpc$vpc <- vpc1
    vpc$pc <- pc1
    vpc$nsources <- nsources
    vpc$scores <- scores
    vpc$data <- dat1
    
    vpc
}








