#' Functions for mAPCA
#'
#' \code{mAPCA} Apply mAPCA to pollution data
#'
#' @param data list of data frames of daily constituent concentrations with date as first column
#' @param tots list of totals (PM2.5) corresponding to each monitor
#' @param nsource number of sources.  If null, uses number of eigenvalues of the correlation matrix greater than cut.
#' @param adjust method of handles adjustment
#' @param mdl list of MDLs for each monitor if applicable
#' @param cut cutoff for eigenvalues (see nsource), default is 1.
#' @param lim Number of days to include data
#' @export
#' @examples
#' data(nycdat)
#' mAPCA(nycdat)
mAPCA <- function(data, tots = NULL, nsources = NULL, adjust = NULL, mdl = NULL, 
    cut = 1, lim = 50) {
    
    N <- length(data)
    
    #setup output
    dates <- list() 
    list.none <- 0
    k <- 1
    #for each monitor
    for (i in 1 : N) {
        
        #Get data and date for monitor i
        dat1 <- data[[i]]
        
        if(class(dat1[, 1]) != "Date") {
            stop("First column is not date")
        }
        dates[[i]] <- dat1[, 1]
        dat1 <- dat1[, -1]
        
        #Get PM2.5 for monitor i
        if(!is.null(tots)){
            pm25 <- tots[[i]]
        }else{
            pm25 <- rowSums(dat1)
        }
        
        
        #concat total data
        dattot <- data.frame(dat1, pm25)
        
        #complete.cases
        dates[[i]] <- dates[[i]][complete.cases(dattot)]
        pm25 <- pm25[complete.cases(dattot)]
        dat1 <- dat1[complete.cases(dat1), ]
        
        #if enough data
        if( nrow(dat1) > lim ) {
            
            
            #create concatenated data
            if (k == 1) {
                stackd <- dat1
                datestot <- dates[[1]]
                tots1 <- pm25
                if(!is.null(mdl)) {mdls <- mdl[[1]]}
            }else{
                stackd <- rbind(stackd, dat1)
                datestot <- c(datestot, dates[[i]])
                tots1 <- c(tots1, pm25)
                if(!is.null(mdl)) {mdls <- rbind(mdls, mdl[[i]])}
            }
            k <- k + 1
            
        #if not enough data
        } else{

            list.none <- c(list.none, i)
            cat("Not enough data, monitor = ", monitors[i], "\n")
            cat("Number of days = ", nrow(dattemp), "\n")
            
        }#end else for > 50 data
        
    }#end loop over monitors
    

    list.none <- list.none[-1]
    
    #get apcs from stacked data
    stackd <- data.frame(datestot, stackd)
    colnames(stackd)[1] <- "date"
    
    
    if(is.null(names(data))) {
        mons <- seq(1, length(data))
    }else{
        mons <- names(data)
    }
    nrows <- sapply(dates, length)
    mons <- unlist(mapply(function(x, y) rep(x, each = y), 
                   mons, nrows, SIMPLIFY = F))
    
    apca <- apca(stackd, tots1, nsources, adjust, mdls, 
        cut, type = "mapca", mons)

    
    #get output		
    list(apca = apca, list.none = list.none,
         dates = dates, datestot = datestot)
}



