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
#' @param list list of length equal to healthdata indicating which monitors correspond to which elements in healthdata
#' @param tots list of vectors corresponding to total concentrations (total PM2.5) for each day and monitor.  If null, uses \code{rowSums(data)}
#' @param method regional source apportionment method (SHARE or mAPCA)
#
#' @export
sharehealth <- function(consdata, healthdata = NULL, list = NULL, 
    formula = NULL, iqrs = NULL,
    lag = NULL, groupvar = NULL, print = F, 
    cut = 1, thres = pi/4,
    tots = NULL, method = "SHARE") {
    
    out <- getsources(consdata, type = tolower(method), tots = tots, list = list)
    
    
    if(!is.null(healthdata)) {
        
        
        #get source concentrations, averaged over community
        if(!is.null(list)) {
           
            out$sourceorder <- combsource(out$sources, list)
        }

        
        #get health effects
#         regcoef <- gethealth(sourcec, healthdata, formula, lag)
#         regcoefREG <- combhealth(regcoef, print = print)
#         
#         out$regcoef <- regcoefREG
#         
#         #iqr increase
#         if(is.null(iqrs)) {
#             iqrs <- out$summary["IQR"]
#         }
#         
#         out$iqrinc <- apply(regcoefREG, 2, percinc, iqrs = iqrs)
#         
    }
    
    

    
    out
}


combsource <- function(outsource, list) {
    
    lens <- sapply(list, length)
    wh <- unlist(list[which(lens > 1)])
    namesU <- unique(names(wh))
    
    sources <- list()
    for(i in 1 : length(list)) {
        if(lens[i] == 1) {
            sources[[i]] <- outsource[list[[i]]]
        }else{
            source1 <- outsource[list[[i]]]
            source1 <- ldply(source1, data.frame)[, -1]
            sources[[i]] <- aggregate(source1, by = list(source1$date), 
                                 FUN= "mean", na.rm = T)[, -1]
        }
    }
    names(sources) <- names(list)
    sources
}





getsources <- function(data, type = "share", tots = tots, list = NULL) {
    
    
    if(type == "share") {        
        share1 <- share(data)
        share <- share1$share
        
        sources <- list()
        for(i in 1 : length(data.rr)) {
            if(is.null(tots)) {
                tots1 <- tots
            }else{
                tots1 <- tots[[i]]
            }
            temp <- apca(data[[i]], tots = tots1)$conc
            temp <- data.frame(data[[i]][, 1], temp)
            colnames(temp) <- c("date", paste0("source", share[[i]]))
            sources[[i]] <- temp
        }
        
        major.sig <- share1$major.sig
        
    }else if(type == "mapca") {
        mapca1 <- mAPCA(data, tots = tots)
        
        apca <- mapca1[["apca"]]
        sourcesM <- apca$conc
        dates <- as.Date(substr(rownames(sourcesM), 1, 10))
        #need to get list
        
        sources <- list()
        unmon <- names(data)
        whM <- which(colnames(sourcesM) == "mons")
        for(i in 1 : length(unmon)) {
            whR <- which(sourcesM$mons == unmon[i])
            temp1 <- sourcesM[whR, -whM]
            temp <- data.frame(dates[whR], temp1)
            colnames(temp) <- c("date", colnames(temp1))
            sources[[i]] <- temp
        }
        
        
        major.sig <- apca[["vmax"]][["loadings"]][1 : (ncol(data[[1]]) - 1), ]
        
        share <- matrix(rep(seq(1, ncol(sources[[i]]) - 1), length(data)), 
            nrow = length(data), byrow = T)
        share <- sapply(apply(share, 1, list), 
            function(x) x[[1]], simplify = F)
    }
    
    #names(sources) <- names(data)
    
    out <- list()
    if(is.null(list)) {
        list <- as.list(names(data))
        names(list) <- names(data)
    }
    
    names(sources) <- names(data)
    
    out$summary <- getsummary(sources, list, major.sig)
    out$sources <- sources
    out$share <- share
    out$major.sig <- major.sig
    
    out
}


sumfun <- function(vec) {
    out <- c(summary(vec, na.rm = T), IQR(vec, na.rm  =T))
    names(out)[length(out)] <- "IQR"
    out
}



#get summary for each monitor
getsummary <- function(sourcec, list, vmax) {
    #get rid of dates
    sourcec <- lapply(sourcec, function(x) {
        x2 <- x[, -1]
        colnames(x2) <- colnames(x)[-1]
        x2
        })
    #summary stats for each monitor
    sumstats <- lapply(sourcec, function(x) t(apply(x, 2, sumfun)))

    #median over unique sources
    unsource <- unique(unlist(sapply(sourcec, colnames, simplify = F)))
    unsource <- sort(unsource[unsource != "sourceInf"])
    
    sumall <- matrix(nrow = length(unsource), ncol = 9)
    colnames(sumall) <- c(colnames(sumstats[[1]]), "monitor", "counties")
    rownames(sumall) <- unsource
    
    #correspond counties and monitors
    list <- unlist(list)
    source1 <- sort(as.numeric(substr(unsource, 7, 7)))
    cons <- vector()
    for(j in source1) {
        
        #find which column
        wh1 <- sapply(sourcec, function(x) {
            s1 <- paste0("source", j)
            if(s1 %in% colnames(x)) {
                which(colnames(x) == s1)
            }else{
                0
            }})
        #get rid of missing
        sumstats1 <- sumstats[wh1 != 0]
        cn <- names(sourcec)[wh1 != 0]
        wh1 <- wh1[wh1 != 0]
        
        #find median summary
        sumj <- mapply(function(x, y) x[y, ], sumstats1, wh1)
        sumall[j, 1 : 7] <- apply(sumj, 1, median, na.rm = T)
        sumall[j, 8] <- ncol(sumj)
        
        
        counties <- unique(names(list)[which(list %in% cn)])
        sumall[j, 9] <- length(counties)
        
        vmax1 <- abs(vmax[, j])
        vmax1 <- sort(vmax1, decreasing = T)
        cons[j] <- paste(names(vmax1)[which(vmax1 > 0.4)], collapse = ", ")
    }

    sumall <- data.frame(sumall, cons)
    sumall    
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
            tln1 <- tlniseC(Y = regcoef1$est, V = regcoef1$se^2, 
                prnt = print)
            if(tln1$converge == "no") {
		tln1 <- tlniseC(Y = regcoef1$est, V = regcoef1$se^2, 
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













#' \code{getallsource} Averages source concentrations for monitors in the same community
#' 
#' @param conc list of source concentrations for each monitor
#' @param share list of sources present at each monitor
#' @param mons1 list of communities where each element is vector of monitors in that community
#' @param allmons 
#
#' @export
getallsource <- function(conc, share, mons1) {
    
    namesDUP <- names(which(sapply(mons1, 
                                   function(x) length(x)) != 1))
    
    allmons <- names(conc)
    #average sources in same city
    scoresMerge <- list()	
    for(i in 1: length(namesDUP)) {
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
            
            #rename sources
            cn <- colnames(scoresMerge[[k]])[-1]
            shareN[[i]] <- as.numeric(substring(cn, 7))
            
            k <- k + 1
        }else{
        	#get source concentrations
            mon1 <- mons1[[i]]
            whMon <- which(allmons %in% mon1)
            temp <- conc[[whMon]]
            
            #rename columns
            cnt <- colnames(temp)[-1]
            cnt <- sapply(strsplit(cnt, "\\."), function(x) x[2])
            colnames(temp) <- c("Date", paste0("source", cnt))
            
            #save output
            scoresAll[[i]] <- temp
            shareN[[i]] <- share[[whMon]]
        }	
    }
    names(scoresAll) <- names(mons1)
    list(scores = scoresAll, share = shareN)
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
                temp <- tlniseC(Y = order[, 1], V = order[, 2]^2, 
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
        # find which monitors have source j
        whT <- sapply(shareMon, function(x) (unsource[j] %in% x))
        scoresJ <- conc[whT]
        wS <- sapply(shareMon[whT], function(x) which(x == unsource[j]))
        whT <- which(whT == T)
        
        
        #get first monitor info
        #add one for date
        scoresJall <- scoresJ[[1]][, c(1, wS[1] + 1)]
        colnames(scoresJall) <- c("Date", paste0("source", unsource[j]))
        
        #if more than 1 monitor, merge data by date
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
            
        #if only one monitor measures, use that monitor
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

