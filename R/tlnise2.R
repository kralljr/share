#' share
#'
#' \code{tlnise} tlnise function from tlnise package
#'
#' Reproduces tlnise function from tlnise package, but gives convergence output.
#' See ?tlnise for more info.
#' 
#' @export
tlniseC <- function (Y, V, w = NA, V0 = NA, prior = NA, 
	N = 1000, seed = NULL, 
          Tol = 1e-06, maxiter = 1000, intercept = TRUE, labelY = NA, 
          labelYj = NA, labelw = NA, digits = 4, brief = 1, prnt = TRUE) 
{
    out.chk <- tlnise:::checkcon(Y, V, w, intercept, prior, prnt)
    Y <- out.chk$Y
    V <- out.chk$V
    w <- out.chk$w
    J <- out.chk$J
    p <- out.chk$p
    q <- out.chk$q
    r <- p * q
    prior <- out.chk$prior
    if (prnt) 
        print(paste("******** Prior Parameter =", prior), quote = FALSE)
    if (missing(V0)) 
        V0 <- rowMeans(V, dims = 2)
    if (max(abs(V0 - diag(p))) < Tol) {
        Ys <- Y
        Vs <- V
        rtV0 <- diag(p)
    }
    else {
        newvars <- tlnise:::standard.f(Y, V, V0)
        Ys <- newvars$Y
        Vs <- newvars$V
        rtV0 <- newvars$rtVo
    }
    if (prnt) {
        cat("\n")
        print("******** Locating Posterior Mode ************ ", 
              quote = FALSE)
        cat("\n")
    }
    Astart <- diag(p)
    out.mode <- tlnise:::postmode.f(Ys, Vs, w, rtV0, prior, Astart, Tol, 
                           maxiter, J, p, q, r)
    modeB0 <- solve(diag(p) + out.mode$newA)
    lfmode <- out.mode$lf
    modeA <- rtV0 %*% out.mode$newA %*% rtV0
    df <- J - q + prior
    Sigma <- modeB0/(df - p - 1)
    eigS <- eigen(Sigma, symmetric = TRUE)
    d <- 1/eigS[[1]]
    if (p == 1) {
        Siginv <- matrix(1/Sigma)
        rtSig <- matrix(sqrt(Sigma))
    }
    else {
        Siginv <- eigS[[2]] %*% diag(d) %*% t(eigS[[2]])
        rtSig <- eigS[[2]] %*% diag(sqrt(1/d))
    }
    ld1 <-  tlnise:::ldet(modeB0)
    tr1 <- tlnise:::tr(modeB0 %*%  Siginv)
    lf0mode <- (df - p - 1) * ld1/2 - tr1/2
    adj <- lf0mode - lfmode
    iter <- out.mode[[16]]
    if (prnt) {
        if (iter < maxiter) {
            print(paste("Converged in", iter, "EM iterations."), 
                  quote = FALSE)
        }
        else {
            print(paste("Did not converge in maxiter =", maxiter, 
                        "EM steps."), quote = FALSE)
        }
        print("Posterior mode of B0:", quote = FALSE)
        cat("\n")
        print(modeB0, digits)
        cat("\n")
        print(paste("lf(modeB0) =", signif(lfmode, digits), "; lf0(modeB0) =", 
                    signif(lf0mode, digits), "; adj =", signif(adj, digits)), 
              quote = FALSE)
    }
    if (is.null(seed)) 
        seed <- ceiling(runif(1) * 1e+08)
    if (seed > 0) 
        seed <- -seed
    if (prnt) {
        cat("\n")
        print("******** Drawing Constrained Wisharts ******** ", 
              quote = FALSE)
        cat("\n")
    }
    pd <- pchisq(d, df)
    if (min(pd) > 1 - Tol) {
        pd <- rep(1, p)
        d <- rep(0, p)
    }
    outU <- tlnise:::rscwish(N, p, df, d, pd, seed)
    Uvals <- outU$Ua
    nvec <- outU$nvec
    nrej <- outU$nrej
    if (prnt) {
        cat("\n")
        print(paste("CWish acceptance rate =", N, "/", nvec[1], 
                    "=", signif(N/nvec[1], digits)), quote = FALSE)
    }
    B0vals <- tlnise:::mammult(rtSig, Uvals, t(rtSig))
    if (prnt) {
        cat("\n")
        print("******** Processing Draws ************** ", quote = FALSE)
        cat("\n")
    }
    out.lf <- tlnise:::lfB0.f(B0vals, Ys, Vs, w, rtV0, df, Siginv, N, 
                     J, p, q, r, prior, adj)
    lr <- out.lf$lrv
    lf <- out.lf$lfv
    lf0 <- out.lf$lf0v
    avewt <- mean(exp(lr - max(lr)))
    if (prnt) {
        print(paste("Average scaled importance weight =", signif(avewt, 
                                                                 digits)), quote = FALSE)
        cat("\n")
    }
    meanA <- rtV0 %*% (out.lf$meanA) %*% rtV0
    if (p == 1) {
        rtA <- sqrt(meanA)
    }
    else {
        rtA <- sqrt(diag(meanA))
    }
    if (prnt) {
        print("Posterior mean estimate of A:", quote = FALSE)
        print(meanA, digits)
        cat("\n")
        print("Between-group SD estimate:", quote = FALSE)
        print(rtA, digits)
        cat("\n")
    }
    if (missing(labelY)) 
        labelY <- 1:J
    if (missing(labelYj)) 
        labelYj <- 1:p
    theta <- rtV0 %*% out.lf$thetahat
    Vtheta <- tlnise:::mammult(rtV0, out.lf$Vthetahat, rtV0)
    if (p == 1) {
        SDtheta <- sqrt(c(Vtheta))
        theta <- c(theta)
        names(SDtheta) <- names(theta) <- labelY
    }
    else {
        SDtheta <- 0 * Y
        for (j in 1:J) SDtheta[, j] <- sqrt(diag(Vtheta[, , j]))
        dimnames(theta) <- dimnames(SDtheta) <- list(labelYj, 
                                                     labelY)
        dimnames(Vtheta) <- list(labelYj, labelYj, labelY)
    }
    gamma <- out.lf$gamhat
    Dgamma <- out.lf$Dgamhat
    GammaMat <- cbind(gamma, sqrt(diag(Dgamma)))
    GammaMat <- cbind(GammaMat, GammaMat[, 1]/GammaMat[, 2])
    if (missing(labelw)) 
        labelw <- seq(0, r - 1, 1)
    dimnames(GammaMat) <- list(labelw, c("est", "se", "est/se"))
    names(gamma) <- labelw
    dimnames(Dgamma) <- list(labelw, labelw)
    
    converge <- ifelse(iter < maxiter, "yes", "no")
    if (brief == 0) 
        out <- list(gamma = GammaMat, theta = theta, SDtheta = SDtheta, 
                    A = meanA, rtA = rtA)
    if (brief == 1) 
        out <- list(gamma = GammaMat, theta = theta, SDtheta = SDtheta, 
                    A = meanA, rtA = rtA, Dgamma = Dgamma, Vtheta = Vtheta, 
                    B0 = B0vals, lr = lr)
    if (brief == 2) 
        out <- list(gamma = GammaMat, theta = theta, SDtheta = SDtheta, 
                    A = meanA, rtA = rtA, Dgamma = Dgamma, Vtheta = Vtheta, 
                    B0 = B0vals, lr = lr, lf = lf, lf0 = lf0, df = df, 
                    Sigma = Sigma, nvec = nvec, nrej = nrej, converge = converge)
    if (brief > 2) 
        out <- list(out.mode, outU = outU, out.lf = out.lf)
    out
}
