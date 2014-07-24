#' Find sources SHared Across a REgion (SHARE)
#'
#' \code{sharehealth} Estimates regional health effects using SHARE
#'
#' This function 1) Applies SHARE, 2) Estimates sources using
#' APCA at each monitor, 3) Estimates health effects for each community, 
#' and 4) Estimates regional health effects for each source
#' 
#' @param consdata list of data frames of daily constituent concentrations with date as first column for each monitor
#' @param healthdata list of data frames corresponding to health counts and covariates with date as first column for each community
#' @param match vector of length equal to consdata indicating which monitors correspond to which elements in healthdata
#' @param tots list of vectors corresponding to total concentrations (total PM2.5) for each day and monitor.  If null, uses \code{rowSums(data)}
#' @param method regional source apportionment method (SHARE or mAPCA)
#
#' @export
sharehealth <- function(consdata, healthdata, match, formula, 
    lag, groupvar, print = F, 
    cut = 1, thres = pi/4,
    tots = NULL, method = "SHARE") {
    
    out <- list()
    
    
    #get SHARE and source concentrations, averaged over community
    source1 <- getsources(consdata, cut, thres, method)
    sourcec <- source1$sources
    out$iqr <- source1$iqr
    
    regcoef <- gethealth(sourcec, healthdata, formula, lag)
    regcoefREG <- combhealth(regcoef, print = print)
    
    out$regcoef <- regcoefREG
    out$iqrinc <- apply(regcoefREG, 2, percinc, iqrs = iqrs)
    
    out
}



combhealth <- function(regcoef, print = F) {
    
    regcoef <- ldply(regcoef, data.frame)
    sources <- unique(regcoef$source)
    
    regcoefREG <- matrix(nrow = length(sources), ncol = 2)
    rownames(regcoefREG) <- sources
    colnames(regcoefREG) <- c("est", "se")
        
    for(i in 1 : length(sources)) {
        regcoef1 <- regcoef[which(regcoef$source == sources[i]), c("est", "se")]
        
        #combine results
        if(nrow(regcoef1) > 1) {
            tln1 <- tlnise(Y = regcoef1$est, V = regcoef1$se^2, 
                prnt = print)
            if(tln1$converge == "no") {
                tln1 <- tlnise(Y = regcoef1$est, V = regcoef1$se^2, 
                    prnt = print, maxiter = 5000)
                if(tln1$converge == "no") {
                    cat("\nDid not converge: source", sources[i], "\n")
                }
            }
            
            regcoefREG[i, ] <- tln1$gamma[1:2]
            
        }else{
            
            regcoefREG[i, ] <- regcoef1
        }
    }
    
    regcoefREG
    
    
}

#groupvar is list is factors to group lag by
#formula includes source, can include ns
gethealth <- function(sourcec, healthdata, formula, lag, groupvar = NULL) {
    
    regcoef <- list()
    for(i in 1 : length(sourcec)) {
        
        #merge health and source data
        colnames(sourcec[[i]])[1] <- tolower(colnames(sourcec[[i]])[1])
        colnames(healthdata[[i]])[1] <- tolower(colnames(healthdata[[i]])[1])
        merged <- merge(sourcec[[i]], healthdata[[i]], by = "date", all = T)
        
        cnALL <- colnames(sourcec[[i]])[-1]
        
        #remove local sources
        whINF <- which(substr(cnALL, 7) == "I")
        if(length(whINF) > 0 ) {
            cnREG <- cnALL[-whINF]
        }
        
        regcoef[[i]] <- matrix(nrow = length(cnREG), ncol = 3)
        rownames(regcoef[[i]]) <- cnREG
        colnames(regcoef[[i]]) <- c("source", "est", "se")
        regcoef[[i]]$source <- cnREG
        
        for(j in 1 : length(cnREG)) {
            #get rid of other sources
            cnREM <- cnALL[-which(cnALL == cnREG[j])]
            merged1 <- merged[, -which(colnames(merged) %in% cnREM)]
            colnames(merged1)[which(colnames(merged1) == cnREG[j])] <- "source"
            merged1[, "source"] <- Lag(merged1[, "source"], k = lag, group = groupvar)
            
            glm1 <- glm(eval(formula), data = merged1, family = "quasipoisson")
            regcoef[[i]][j, c("est", "se")] <- summary(glm1)$coef["source", c(1, 2)]
        }
    }
    
    regcoef
    
}




getsources <- function(consdata, match = NULL, cut = 1, thres = pi/4, tots = NULL, 
    method = "SHARE") {
   
    if(tolower(method) == "share") {
        
        #perform SHARE
        share1 <- share(data = consdata, cut = cut, thres = thres)
        share <- share1$share
        sources <- share1$Sources
        reg <- share1$major.sig
        
        #apply APCA to each monitor
        sourcec <- list()
        for(i in 1 : length(consdata)) {
            nf1 <- length(share[[i]])
            temp <- apca(data = consdata[[i]], nsources = nf1, tots = tots)$conc
            colnames(temp) <- c("date", paste0("source", share[[i]]))
            sourcec[[i]] <- temp
        }
        
        
        out <- list(share = share, sources = sources, reg = reg,
            sources = sourcec)
        
    }else if(tolower(method) == "mAPCA") {
        
        #perform mAPCA
        mapca <- mAPCA(data = consdata, lim = 50, tots = tots)$apca
        mapcasource1 <- mapca$conc
        mapca <- as.matrix(mapca[["vmax"]]$load)
        
        #get list of mapca results by monitor
        cn <- which(colnames(mapcasource1) == "mons")
        mons <- unique(mapcasource1$mons)
        sourcec <- list()
        cn <- c("date", paste0("source", seq(1, ncol(mapcasource1) - 1)))
        
        for(i in 1 : length(mons)) {
            
             temp <- mapcasource1[which(mapcasource1$mons == mons[i]), -cn]
             colnames(temp) <- cn
             sourcec[[i]] <- temp
        }
        
        
        out <- list(mapca = mapca, sources = sourcec)
        
    }else{
        stop("Method is not recognized.  Must be either mAPCA or SHARE")
    }
    
    #reorder source concentration results
    if(!is.null(match)) {
        out$sourcec <- combsource(out$sourcec, match)
    }
    out$summary <- getsummary(out$sourcec)

    out
}



sumfun <- function(vec) {
    out <- c(summary(vec, na.rm = T), IQR(vec, na.rm  =T))
    names(out)[length(out)] <- "IQR"
    out
}


getsummary <- function(sourcec) {
    #summary stats for each monitor
    sumstats <- lapply(sourcec, function(x) t(apply(x, 2, sumfun)))
    iqrs <- median()
    dat <- ldply(sourcec, data.frame)
    aggregate(dat, aggregate(dat, by = list(dat$date), 
        FUN= "mean", na.rm = T)[, -1])
    
}

#columns of sources must be named by source according to SHARE or whatever
#need first column to be date for each
combsource <- function(sources, match) {
    
    #which duplicated
    dups <- match[which(duplicated(match))]
    comms <- sort(unique(match))
    ind <- 1 * (comms %in% dups)
    
    sources1 <- list()
    for(i in 1 : length(comms)) {
        
        #if duplicated monitors in community,  average
        if(ind[i] == 1) {
            #find which monitors belong in city dups[i]
            dat <- sources[[which(match == dups[i])]]
            
            #average by date
            dat <- ldply(dat, data.frame)
            colnames(dat)[1] <- tolower(colnames(dat)[1])
            sources1[[i]] <- aggregate(dat, by = list(dat$date), 
                 FUN= "mean", na.rm = T)[, -1]
            
            
        #otherwise, take that one monitor
        }else{
            
            sources1[[i]] <- sources[[which(match == comms[i])]]
        }
    }
    names(sources1) <- comms
    sources1
    
}






#create list of scores by city
#conc is estimated source conc for each monitor
#share is what is shared between monitors
getallsource <- function(conc, share, mons1, allmons) {
    
    namesDUP <- names(which(sapply(mons1, 
                                   function(x) length(x)) != 1))
    
    #average sources in same city
    scoresMerge <- list()	
    for(i in 1: length(namesDUP)) {
        # if(i == 8) {browser()}
        # print(i)
        scoresMerge[[i]] <- dupfun(namesDUP[i], conc, mons1, allmons, share)
    }
    
    
    #for each city
    k <- 1
    scoresAll <- list()
    shareN <- list()
    for(i in 1 : length(mons1)) {
        
        city <- names(mons1)[i]
        # print(city)
        if(city %in% namesDUP) {
            scoresAll[[i]] <- scoresMerge[[k]]
            
            cn <- colnames(scoresMerge[[k]])[-1]
            # cnt <- sapply(strsplit(cn, "\\."), function(x) x[2])
            shareN[[i]] <- as.numeric(substring(cn, 7))
            
            k <- k + 1
        }else{
            mon1 <- mons1[[i]]
            whMon <- which(allmons %in% mon1)
            temp <- conc[[whMon]]
            cnt <- colnames(temp)[-1]
            cnt <- sapply(strsplit(cnt, "\\."), function(x) x[2])
            colnames(temp) <- c("Date", paste0("source", cnt))
            scoresAll[[i]] <- temp
            shareN[[i]] <- share[[whMon]]
        }	
    }
    names(scoresAll) <- names(mons1)
    list(scoresAll, shareN)
}




getsources <- function(dat, nf1, type1, bstar1 = NULL, stdrow = T,
                       shares = NULL, tots = NULL) {
    
    if(type1 == "true") {
        sourceconc <- dat
        
    }else if(type1 == "apca") {
        temp <- abspca(dat, tots, nfactors = nf1, bstar1 = bstar1)
        sourceconc <- temp[[1]]
        
    }
    
    if(!is.null(shares)) {
        sourceconc1 <- try(sourceconc[, shares])
        if(class(sourceconc1) == "try-error") {
            browser()
        }else{
            sourceconc <- sourceconc1
        }
    }
    sourceconc
    
    
}






tlncomb <- function(ests, Sources, share, sinkf = NULL, names = cities) {
    tln <- list()
    
    if(!is.null(sinkf)) {
        sink(file = sinkf)
    }
    
    names1 <- list()
    Sources <- sort(Sources)
    for(i in Sources) {
        tln[[i]] <- list()
        tln[[i]][[1]] <- matrix(nrow = 3, ncol = 2)	
        tln[[i]][[2]] <- list()
        tln[[i]][[3]] <- vector()
        
        whin <- sapply(share, function(x) {
            whin <- which(x == as.character(i))
            ifelse(length(whin) > 0, whin, 0)
        })
        whn0 <- which(whin != 0)
        names1[[i]] <- names(whn0)
        
        if(length(whn0) > 0) {
            for(lag in 1 : 3) {
                order <- matrix(nrow = length(ests[[lag]]), ncol = 2)
                
                for(j in whn0) {
                    order[j, ] <- ests[[lag]][[j]][whin[j], ]
                    
                }
                rownames(order) <- names
                order <- order[complete.cases(order), ]
                
                if(lag == 1) {
                    tln[[i]][[4]] <- array(dim = c(nrow(order), 2, 3))
                    dimnames(tln[[i]][[4]]) <- list(rownames(order), 
                                                    c("est", "se"), 
                                                    paste("lag", seq(0, 2)))
                }
                
                print(i)
                print(c("TLN:", i))
                temp <- tlnise(Y = order[, 1], V = order[, 2]^2, 
                               maxiter = 5000)
                
                tln[[i]][[1]][lag, ] <- temp$gamma[1:2]
                tln[[i]][[3]][lag] <- temp$A
                tln[[i]][[2]][[lag]] <- order
                tln[[i]][[4]][, , lag] <- cbind(temp$theta, temp$SDtheta)
                
                lagUSE <- lag - 1
            }#end loop over lag
            
            
        }else if(length(whn0) == 1){#end whn0
            temp <- ests[[whn0]][whin[whn0], ]
            temp <- t(as.matrix(temp))
            rownames(temp) <- names1[[i]]
            
            tln[[i]][[1]][lag, ] <- temp
            tln[[i]][[2]][[lag]] <- temp
            
        }else{
            
            sink()
            stop("no source error")
        }
        
        rownames(tln[[i]][[1]]) <- paste0("lag", seq(0, 2))
        names(tln[[i]][[2]]) <- paste0("lag", seq(0, 2))
        names(tln[[i]][[3]]) <- paste0("lag", seq(0, 2))
    }#end loop over source
    
    sink()
    list(tln, names1)
}


#nameC is name of city
#conc is matrix of dates and estimated source conc
#share is share matrix from domatchsim
dupfun <- function(nameC, conc, mons1, allmons, share = share) { 
    
    #find monitors corresponding to city
    mons <- mons1[[nameC]]
    whMon <- which(allmons %in% mons)
    conc <- conc[whMon]
    shareMon <- share[whMon]
    
    
    #find unique sources in this city
    unsource <- sort(unique(unlist(shareMon)))
    whInf <- which(is.infinite(unsource))
    if(length(whInf) > 0) {
        unsource <- unsource[-whInf]
    }
    
    #merge all sources
    scoresMERGE <- list()
    for(j in 1 : length(unsource)) {
        # print(j)
        whT <- sapply(shareMon, function(x) (unsource[j] %in% x))
        scoresJ <- conc[whT]
        wS <- sapply(shareMon[whT], function(x) which(x == unsource[j]))
        whT <- which(whT == T)
        
        
        #get first monitor info
        #add one for date
        scoresJall <- scoresJ[[1]][, c(1, wS[1] + 1)]
        colnames(scoresJall) <- c("Date", paste0("source", unsource[j]))
        
        #if more than 1 monitor, merge data
        if(length(whT) > 1) {
            for(k in 2 : length(whT)) {
                
                scoresJk <- scoresJ[[k]][, c(1, wS[k] + 1)]
                scoresJall <- merge(scoresJall, scoresJk, 
                                    by = "Date", all.x = T, all.y = T)
            }
            
            #average source concentrations
            temp <- scoresJall[, -1]
            temp <- apply(temp, 1, mean, na.rm = T)
            scoresMERGE[[j]] <- data.frame(scoresJall[, 1], temp)
            colnames(scoresMERGE[[j]]) <- c("Date", paste0("source", unsource[j]))
            
        }else{
            
            scoresMERGE[[j]] <- scoresJall
        }
        
        
        if(j == 1) {
            tempJ <- scoresMERGE[[j]]
        }else{
            tempJ <- merge(tempJ, scoresMERGE[[j]], 
                           by = "Date", all.x = T, all.y = T)
            
        }
        
    }#end loop over j
    tempJ
    
}


#function to read TLNise convergence from file
readTLNconv <- function(filename) {
    text <- scan(filename, what = "string")
    
    
    #find location of all convergence items
    whTLN <- which(text == "TLN:")
    whTLN <- c(whTLN, whTLN + 1, whTLN + 2)
    whC <- which(text == "Converged")
    whC <- c(whC, whC + 2)
    whNC <- which(text == "converge")
    text2 <- text[sort(c(whTLN, whC, whNC))]
    whTLN2 <- which(text2 == "TLN:")
    
    
    #for each iteration (whTLN2), concat convergence info
    outs <- vector()
    for(i in 1 : (length(whTLN2) - 1)) {
        
        start <- whTLN2[i]
        stops <- whTLN2[i + 1] - 1
        outs[i] <- paste(text2[(start : stops)], collapse = "/")
    }
    
    #find non-convergence
    outs2 <- outs[which(substring(outs, 10, 10) == "c")]
    
    #make output nice
    outs2 <- strsplit(substr(outs2, 6, 8), "/")
    nonconv <- sapply(outs2, function(x) paste0("source=", 
                                                paste(x, collapse = ", method=")))
    
    list(nonconv, paste(text2, collapse = ""))	
    
}
