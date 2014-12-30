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
    formula = NULL, 
    lag = NULL, groupvar = NULL, print = F, 
    cut = 1, thres = pi/4,
    tots = NULL, method = "SHARE") {
    
    out <- getsources(consdata, type = tolower(method), tots = tots, list = list)
    
    
    if(!is.null(healthdata)) {
        
        
        #get source concentrations, averaged over community
        if(!is.null(list)) {
           
            sourcec <- combsource(out$sources, list)
            out$sourceorder <- sourcec
        }else{
            
            sourcec <- out$sources
        }

        
        #get health effects
        regcoef <- list()
        for(i in 1 : length(healthdata)) {
            print(i)
            regcoef[[i]] <- get.hosp(dat = healthdata[[i]], sources = sourcec[[i]], 
                lag = lag, formu = form1, gv = gv)
        }
        
        regcoefREG <- combhealth(regcoef, print = print)
        sources <- as.numeric(substr(rownames(regcoefREG), 7, 7))
        regcoefREG <- regcoefREG[order(sources), ]
        
        out$hcounty <- regcoef
        out$regcoef <- regcoefREG
        
        #iqr increase (IQRS from data)
        iqrs <- out$summary["IQR"]
        
        lb <- regcoefREG[, 1] - 1.96 * regcoefREG[, 2]
        lb <- percinc(lb, scale = iqrs)
        ub <- regcoefREG[, 1] + 1.96 * regcoefREG[, 2]
        ub <- percinc(ub, scale = iqrs)
        
        est <- percinc(regcoefREG[, 1], scale = iqrs)
        
        percincci <- data.frame(est, lb, ub)
        colnames(percincci) <- c("est", "lb", "ub")
        out$iqrinc <- percincci
        
    }
    
    

    
    out
}


# dat is hosp dat for one county
# sources is sources for one county
# lag is lag
get.hosp <- function(dat, sources, lag, formu, outcome = "cardio", gv = "agecat") {
    
    #merge with sources
    colnames(dat) <- tolower(colnames(dat))
    colnames(sources) <- tolower(colnames(sources))
    merged <- merge(dat, sources, all.x = T, 
                    by = "date")
    
    #get lag info	
    #which columns are sources
    whSource <- which(substr(colnames(merged), 1, 4) == "sour")
    lagSource <- vector(, length = nrow(merged))
    
    temp2 <- rep(NA, nrow(merged))
    #for each source, #lag by agecat
    for(i in 1 : length(whSource)) {
        
        sour <- merged[, whSource[i]]	
        #for each source with at least 1 day of data
        if(length(which(!is.na(sour))) > 0) {
            
            #if group variable
            if(!is.null(gv)) {
                gv1 <- merged[, gv]
            }else{
                gv1 <- gv
            }
            
            temp <- Lag(sour, k = lag, group = gv1)
            lagSource <- cbind(lagSource, temp)
        
        #else all NA    
        }else{
            lagSource <- cbind(lagSource, temp2)
        }
    }
    lagSource <- lagSource[, -1]
    
    #add in lagged sources
    merged[, whSource] <- lagSource
    
    #number of years of data (should be 12)
    years <- length(unique(substr(merged$date, 1, 4)))
    
    #get formula
    formUSE1 <- paste0(formu, years)
    
    #set up estimates
    ests <- matrix(nrow = length(whSource), ncol = 2)
    #for each source
    for(l in 1 : length(whSource)) {
        # print(c("l", l))
        covar1 <- paste(colnames(merged)[whSource[l]], collapse = "+")
        formUSE <- paste0(outcome, " ~", covar1, "+", formUSE1, ")")
        
        #run model
        options(warn = 2)
        glm1 <- try(glm(formula = eval(formUSE), 
                        data = merged, family = "quasipoisson",
                        offset = log(denom)), silent = T)
        options(warn = 1)
        tf <- suppressWarnings(class(glm1) == "try-error")
        if(tf[1]) {
            glm1 <- try(glm(formula = eval(formUSE), 
                data = merged, family = "quasipoisson",
                offset = log(denom), control = list(epsilon = .00001)), 
                silent = T)
        }
        
        #save results
        tf <- suppressWarnings(class(glm1)  != "try-error")
        if(tf[1]) {
            out1 <- summary(glm1)$coef
            whS <- which(substr(rownames(out1), 1, 4) == "sour")
            out <- out1[whS, c(1, 2)]
        }else{
            out <- c(NA, NA)
        }
        ests[l, ] <- out
        
    }
    
    rownames(ests) <- colnames(merged)[whSource]
    colnames(ests) <- c("est", "se")
    ests
}




combsource <- function(outsource, list) {
    
    lens <- sapply(list, length)
    wh <- unlist(list[which(lens > 1)])
    namesU <- unique(names(wh))
    
    sources <- list()
    for(i in 1 : length(list)) {
        if(lens[i] == 1) {
            sources[[i]] <- outsource[[list[[i]]]]
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
        for(i in 1 : length(data)) {
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
    
    regcoef <- sapply(regcoef, function(x) {
        x <- data.frame(rownames(x), x)
        colnames(x) <- c("source", "est", "se")
        x
        }, simplify = F)
    regcoef <- ldply(regcoef, data.frame)
    sources <- unique(regcoef$source)
    sources <- sources[substr(sources, 7, 7) != "i"]
    regcoefREG <- matrix(nrow = length(sources), ncol = 2)
    rownames(regcoefREG) <- sources
    colnames(regcoefREG) <- c("est", "se")
        
    for(i in 1 : length(sources)) {
        regcoef1 <- regcoef[which(regcoef$source == sources[i]), c("est", "se")]
        
        #combine results
        if(nrow(regcoef1) > 1) {
            tln1 <- tlniseC(Y = regcoef1$est, V = regcoef1$se^2, 
                prnt = print, brief = 2)
            if(tln1$converge == "no") {
		tln1 <- tlniseC(Y = regcoef1$est, V = regcoef1$se^2, 
                    prnt = print, maxiter = 5000, brief = 2)
                if(tln1$converge == "no") {
                    cat("\nDid not converge: source", sources[i], "\n")
                }
            }
            
            regcoefREG[i, ] <- tln1$gamma[1:2]
            
        }else{
            
            regcoefREG[i, ] <- as.matrix(regcoef1)
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



#' \code{percinc} Percent increase
#' 
#' @param beta regression coefficients from poisson model
#' @param scale amount for percent increase
#' @export
percinc <- function(beta, scale = 10) {
    
    if(!is.na(beta[1])) {
        100 * (exp(beta * scale ) - 1)
    }else{
        beta
    }
}