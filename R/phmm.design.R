### ======================================
### phmm.design
### ======================================
## Adapted from bayesurvreg.design in 
## bayesSurv package (Arnost Komarek, 2005)
##
##  - extract the design information
##
## 01/25/2008
##
phmm.design <- function (m, formula, random, data)
{
    tempF <- c("", "formula", "data", "subset", "na.action")
    mF <- m[match(tempF, names(m), nomatch = 0)]
    mF[[1]] <- as.name("model.frame")
    special <- c("cluster")
    TermsF <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
    mF$formula <- TermsF
    mF <- eval(mF, parent.frame())
    Y <- model.extract(mF, "response")
    Yinit <- Y
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object. ")
    type <- attr(Y, "type")
    if (type == "counting") 
        stop("Invalid survival type ('counting' is not implemented). ")
    n <- nrow(Y)
    nY <- ncol(Y)
    cluster <- attr(TermsF, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
        tempc <- untangle.specials(TermsF, "cluster", 1:10)
        ord <- attr(TermsF, "order")[tempc$terms]
        if (any(ord > 1)) 
            stop("Cluster can not be used in an interaction")
        cluster <- strata(mF[, tempc$vars], shortlabel = TRUE)
        dropx <- tempc$terms
    }
    if (length(dropx)) {
        newTermsF <- TermsF[-dropx]
        attr(newTermsF, "intercept") <- attr(TermsF, "intercept")
    }
    else {
        newTermsF <- TermsF
    }
    attr(newTermsF, "intercept") <- 1
    Zinit <- model.matrix(newTermsF, mF)
    rnamesZ <- row.names(Zinit)
    cnamesZ <- colnames(Zinit)
    nZ <- ncol(Zinit) - 1
    cnamesZ <- cnamesZ[-1]
    if (nZ) {
        Z <- Zinit
        Z <- Z[, -1]
        attr(Z, "assign") <- attr(Zinit, "assign")[-1]
        attr(Z, "contrasts") <- attr(Zinit, "contrasts")
    }
    else {
        Z <- NULL
    }
    indb <- if (nZ) 
        rep(-1, nZ)
    else 0
    randomInt <- FALSE
    nrandom <- 0
    if (!missing(random)) {
        if (!length(cluster)) 
            stop("You have to indicate clusters when you want to include some random effects. ")
        tempR <- c("", "random", "data", "subset", "na.action")
        mR <- m[match(tempR, names(m), nomatch = 0)]
        mR[[1]] <- as.name("model.frame")
        names(mR)[2] <- "formula"
        TermsR <- if (missing(data)) 
            terms(random)
        else terms(random, data = data)
        lTR <- length(attr(TermsR, "variables"))
        if (lTR == 1 & !attr(TermsR, "intercept")) {
            names.random <- character(0)
        }
        else {
            if (lTR == 1 & attr(TermsR, "intercept")) {
                randomInt <- TRUE
                names.random <- character(0)
            }
            else {
                mR$formula <- TermsR
                mR <- eval(mR, parent.frame())
                if (attr(TermsR, "intercept")) {
                  randomInt <- TRUE
                  names.random <- colnames(model.matrix(TermsR, mR))[-1]
				  W <- model.matrix(TermsR, mR)[-1,]
                }
                else {
                  names.random <- colnames(model.matrix(TermsR, mR))
				  W <- model.matrix(TermsR, mR)
                }
            }
            nrandom <- 1 * randomInt + length(names.random)
            #if (sum(names.random %in% cnamesZ) != nrandom - 1 * 
            #    randomInt) 
            #    stop("Each random effect has to have also its fixed counterpart.")
            find.indeces <- function(all.eff) {
                where <- names.random %in% all.eff
                if (!sum(where)) 
                  return(-1)
                if (sum(where) > 1) 
                  stop("Error, contact the author.")
                index <- (1:length(names.random))[where]
                if (!randomInt) 
                  index <- index - 1
                return(index)
            }
            if (nZ) 
                indb <- as.numeric(apply(matrix(cnamesZ, ncol = 1), 
                  1, find.indeces))
        }
    }
    else {
        names.random <- character(0)
    }
    if (randomInt) 
        names.random <- c("(Intercept)", names.random)
    nfixed <- nZ # - (nrandom - 1 * randomInt) # want to be able to include
    n.factors <- 0                             # both fixed and random and random
    n.in.factors <- NULL                       # effects for same covariate
    factors <- NULL
    if (nZ) {
        temp <- attr(Z, "assign")
        if (length(temp) == 1) 
            factors <- 0
        else {
            factors <- numeric(length(temp))
            n.in.factors <- numeric(0)
            temp <- temp - c(0, temp[1:(length(temp) - 1)])
            i <- length(temp)
            while (i >= 1) {
                if (temp[i] == 0) {
                  n.factors <- n.factors + 1
                  factors[i] <- n.factors
                  n.in.factor <- 1
                  while (temp[i - 1] == 0) {
                    i <- i - 1
                    factors[i] <- n.factors
                    n.in.factor <- n.in.factor + 1
                  }
                  i <- i - 1
                  factors[i] <- n.factors
                  n.in.factors <- c(n.in.factor + 1, n.in.factors)
                }
                else n.in.factors <- c(1, n.in.factors)
                i <- i - 1
            }
        }
        if (length(temp) != nZ) 
            stop("Something is wrong, contact the author.")
    }
    if (length(cluster)) {
        ordering <- order(cluster)
        Y <- Y[ordering, ]
        cluster <- cluster[ordering]
        rnamesZ <- rnamesZ[ordering]
        if (nZ) {
            namesZ <- cnamesZ
            if (nZ == 1) 
                Z <- matrix(Z[ordering], ncol = 1)
            else Z <- as.matrix(Z[ordering, ])
            colnames(Z) <- namesZ
        }
        if (!missing(random)) {
            namesw <- colnames(W)
            if (lTR == 1) 
                W <- matrix(W[ordering], ncol = 1)
            else W <- as.matrix(W[ordering, ])
            colnames(W) <- namesw
        }
        ncluster <- length(attr(cluster, "levels"))
        helpf <- function(cl) {
            return(sum(cluster %in% attr(cluster, "levels")[cl]))
        }
        nwithin <- apply(matrix(1:ncluster, ncol = 1), 1, "helpf")
    }
    else {
        if (nZ) {
            namesZ <- cnamesZ
            if (nZ == 1) 
                Z <- matrix(Z, ncol = 1)
            else Z <- as.matrix(Z[, ])
            colnames(Z) <- namesZ
            cluster <- 1:n
            ncluster <- n
            nwithin <- rep(1, n)
        }
    }
    if (type == "interval") {
        stop("Interval censoring is not supported.")
    }
    if (!all(is.finite(Y))) 
        stop("Invalid survival times for this distribution (infinity on log-scale not allowed). ")
    names.cluster = cluster
	cluster = rep(1:ncluster,table(cluster))
    design <- list(n = n,
				   ncluster = ncluster,
				   nwithin = nwithin,
				   nY = nY,
				   nZ = nZ,
				   nfixed = nfixed,
				   nrandom = nrandom,
				   randomInt = randomInt,
				   Y = Y,
				   Z = Z,
				   W = W,
				   Yinit = Yinit,
				   Zinit = Zinit,
				   cluster = cluster,
				   names.cluster = names.cluster,
				   indb = indb,
				   rnamesZ = rnamesZ,
				   names.random = names.random,
				   factors = factors,
				   n.factors = n.factors,
				   n.in.factors = n.in.factors)
    return(design)
}
