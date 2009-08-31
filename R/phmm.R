phmm <- function (formula, random, data = parent.frame(), subset, 
	na.action = na.fail, Sigma="identity", varcov = "diagonal",
	NINIT = 10, VARSTART=1, MAXSTEP=100, CONVERG=90, Gbs=100, Gbsvar=1000, verbose=FALSE) 
{
	call <- match.call(expand.dots = TRUE)
    m <- match.call(expand.dots = FALSE)
    des <- phmm.design(m, formula, random, data)
	
	if(Sigma=="identity") Sigma0 = diag(1,des$nrandom)
	invSigma = rbind(0,cbind(0,solve(Sigma0)))
	detSigma = det(Sigma0)
	Sigma    = rbind(0,cbind(0,solve(Sigma0)))
	
	if(varcov=="diagonal"){ varcov = 2
	}else{ stop(paste("\nUnknown structure specified for var-covariance matrix:",
		varcov))}

    if(verbose){
		cat("\nProportional Hazards Mixed-Effects Model fit with MCMC-EM\n")
		date0 <- date()
	}

	fit <- .C("phmm", 
		X = as.numeric(c(0,des$Y[,1])), 
		Z = rbind(0, as.matrix(des$Z)), 
		W = rbind(0, as.matrix(des$W)),
		delta = as.integer(c(0,des$Y[,2])), 
		cluster = as.integer(c(0,des$cluster)),
		varcov = as.integer(varcov),
		Sigma = Sigma,
		invSigma = invSigma,
		detSigma = as.double(detSigma),
		NINIT=as.integer(10),
		MAXSTEP=as.integer(MAXSTEP),
		CONVERG=as.integer(CONVERG),
		emstep = as.integer(0),
		Gbs = as.integer(Gbs),
		Gbsvar = as.integer(Gbsvar),
		n = as.integer(des$n),
		nfixed = as.integer(des$nfixed),
		nrandom = as.integer(des$nrandom),
		bhat = double(des$nclust*des$nrandom),
		sdbhat = double(des$nclust*des$nrandom),
		steps = double((MAXSTEP+1)*(des$nfixed+des$nrandom)),
		var = double((des$nfixed+des$nrandom+des$n+1)^2),
		verbose = as.integer(verbose),
		bridgeC = double(1),
		laplacetot = double(1),
		llaplace = double(1), 
		limport = double(1), 
		lbridge = double(1),
		lambda = double(des$n+1),
		Lambda = double(des$n+1),
		PACKAGE="phmm" )

	if(varcov==2) fit$varcov = "diagonal"
	
	class(fit) <- "phmm"

	fit$call <- call
	fit$formula <- formula
	fit$random <- random
	fit$verbose <- as.logical(fit$verbose)
	fit$steps <- matrix(fit$steps, nrow=MAXSTEP+1, byrow=TRUE)
	colnames(fit$steps) <- c(colnames(des$Z), paste("var(", des$names.random, ")", sep=""))
	rownames(fit$steps) <- 0:MAXSTEP
	
	fit$var <- matrix(fit$var, nrow=des$nfixed+des$nrandom+des$n+1, byrow=FALSE)
	fit$var <- fit$var[1 + 1:(des$nfixed+des$nrandom), 1 + 1:(des$nfixed+des$nrandom)]
	rownames(fit$var) <- colnames(fit$var) <- c(colnames(des$Z), des$names.random)
	fit$varFix <- fit$var[1:des$nfixed, 1:des$nfixed]
	
	fit$Sigma0 <- Sigma0
	fit$Sigma <- matrix(fit$Sigma, nrow=des$nrandom+1, byrow=FALSE)
	fit$Sigma <- fit$Sigma[1 + 1:des$nrandom, 1 + 1:des$nrandom]
	if(des$nrandom==1){ names(fit$Sigma)=paste("var(", des$names.random, ")", sep="")
	}else{
		colnames(fit$Sigma) = des$names.random
		rownames(fit$Sigma) = des$names.random
	}
	fit$invSigma <- matrix(fit$invSigma, nrow=des$nrandom+1, byrow=FALSE)
	fit$invSigma <- fit$invSigma[1 + 1:des$nrandom, 1 + 1:des$nrandom]

	fit$coefficients = fit$steps[MAXSTEP+1,1:fit$nfixed]
	names(fit$coefficients) <- colnames(des$Z)
	fit$bhat <- matrix(fit$bhat, nrow=des$nclust, byrow=TRUE)
	fit$sdbhat <- matrix(fit$sdbhat, nrow=des$nclust, byrow=TRUE)
	
	rownames(fit$bhat) <- rownames(fit$sdbhat) <- unique(des$cluster) 
	colnames(fit$bhat) <- colnames(fit$sdbhat) <- colnames(des$W)	
	
	fit$X = des$Y[,1]
	fit$Z = as.matrix(des$Z)
	fit$W = as.matrix(des$W)
	fit$delta = des$Y[,2] 
	fit$cluster = des$cluster
	
	fit$Lambda = fit$Lambda[-1]
	fit$lambda = fit$lambda[-1]
	
	fit$bhat.long <- as.matrix(merge(
	   cbind(cluster=as.numeric(rownames(fit$bhat)), as.matrix(fit$bhat)),
	   cbind(cluster=as.numeric(as.character(fit$cluster)))))[,-1]
	
	if(fit$nrandom > 2) fit$lbridge <- NULL
	fit$loglik <- c(Conditional=loglik.cond(fit),
		Laplace=fit$llaplace,RIS=fit$limport,BS=fit$lbridge)
		
	fit$linear.predictors <- linear.predictors(fit)
	
	fit <- fit[c('X', 'Z', 'W', 'delta', 'cluster', 'varcov', 'Sigma', 
	   'invSigma', 'detSigma', 'NINIT', 'MAXSTEP', 'CONVERG', 'emstep', 
	   'Gbs', 'Gbsvar', 'n', 'nfixed', 'nrandom', 'bhat', 'bhat.long', 
	   'sdbhat', 'steps', 'var', 'verbose', 'loglik', 'lambda', 'Lambda', 
	   'call', 'formula', 'random', 'varFix', 'coefficients', 
	   'linear.predictors')]
	fit <- c(fit, des["Y"])
	
	class(fit) <- "phmm"
    return(fit)
}

loglik.cond <- function (x) UseMethod("loglik.cond")
loglik.cond.phmm <- function(x){
	#Function to compute conditional log-likelihood
	phmm.cond.loglik(time=x$X, delta=x$delta, z=x$Z, beta=x$coef, w=x$W, b=as.matrix(x$bhat.long))
}

phmm.cond.loglik <- function(time, delta, z, beta, w, b){
	#Function to compute conditional log-likelihood
    z <- as.matrix(z)
    wb <- matrix(0,nrow=nrow(w),ncol=ncol(w))
    for(i in 1:ncol(w)){
	  if(length(w[,i])!=length(b[,i])) stop("length(w[,i])!=length(b[,i])")
	  wb[,i] <- w[,i]*b[,i]
	  }
    if(!is.null(dim(wb))) wb <- apply(wb,1,sum)
    numerator <- exp(z%*%beta+wb)
    denominator <- unlist(lapply(time,
      FUN=function(x){
          sum(exp((z%*%beta+wb)[time>=x]))
          }))
    sum(ifelse(delta,1,0)*log(numerator/denominator))
}

linear.predictors <- function (x) UseMethod("linear.predictors")
linear.predictors.phmm <- function(x){
    #Function to compute linear predictors
    #cluster is assumed to be sorted
	
	z=x$Z; beta=x$coef; w=x$W; b=as.matrix(x$bhat.long)
	
    z <- as.matrix(z)
    wb <- matrix(0,nrow=nrow(w),ncol=ncol(w))
    for(i in 1:ncol(w)){
	  if(length(w[,i])!=length(b[,i])) stop("length(w[,i])!=length(b[,i])")
	  wb[,i] <- w[,i]*b[,i]
	  }
    if(!is.null(dim(wb))) wb <- apply(wb,1,sum)
    return(z%*%beta+wb)
    }

# creates pseudo data for poisson PHMM model
# from object of class "phmm"
pseudoPoisPHMM <- function (x) UseMethod("pseudoPoisPHMM")
pseudoPoisPHMM.phmm <- function(x){	
	dd <- cbind(x$cluster, x$Z, x$W)
	group <- apply(dd,1,paste,collapse="XX")
	groups <- unique(group)
	dd <- cbind(x$X, x$delta, x$linear.predictors, dd)
	colnames(dd) <- c("time", "delta", "linear.predictors",
		"cluster",
		paste("z", 1:x$nfixed, sep=''),
		paste("w", 1:x$nrandom, sep=''))
	ddext <- c()
	for(t in sort(unique(dd[dd[,"delta"]==1,"time"]))){
		tdd <- c()
		for(i in 1:length(groups)){
			tdd <- rbind(tdd,
			unlist(c(t,sum(dd[group==groups[i],"time"]>=t),
			sum(dd[group==groups[i]&dd[,"delta"]==1,"time"]==t),
			dd[group==groups[i]&!duplicated(group==groups[i]),
				c("cluster",
				paste("z", 1:x$nfixed, sep=''),
				paste("w", 1:x$nrandom, sep=''),
				"linear.predictors")])))
			}
		ddext <- rbind(ddext,tdd)
		}
	colnames(ddext) <- c('time','N','m',"cluster",
				paste("z", 1:x$nfixed, sep=''),
				paste("w", 1:x$nrandom, sep=''),
				"linear.predictors")
	times <- sort(unique(dd[dd[,"delta"]==1,"time"]))
	timematrix <- matrix(0,nrow(ddext),ncol=length(times))
	colnames(timematrix) <- paste("t",1:length(times),sep='')
	for(i in 1:length(times) ){		
		timematrix[ddext[,"time"]==times[i],paste("t",i,sep='')] <- 1
	}
	ddext <- cbind(ddext,timematrix)
	ddext <- data.frame(ddext)
	ddext <- ddext[ddext$N!=0,]
	ddext <- ddext[order(ddext$cluster),]
	return(ddext)
}

AIC.phmm <- function(object, ..., k=2){
	if(object$varcov=="diagonal"){ 
		return(-2*object$loglik+k*(object$nrandom+object$nfixed))
	}else{ stop(paste("\nUnknown structure specified for var-covariance matrix:",
		object$varcov))}
}

cAIC <- function(object, ..., k=2){
	-2*object$loglik["Conditional"]+k*traceHat(object)
}
	
print.phmm <-
 function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nProportional Hazards Mixed-Effects Model fit by MCMC-EM\n")
  cat("  Model:", deparse(x$call$formula),"\n")
  cat("  Data:", deparse( x$call$data ), "\n")
  cat("  Log-likelihood:\n")
	print(x$loglik[!is.null(x$loglik)], digits = digits, ...)
  Ffix <- paste("~", strsplit(deparse(x$call$formula), "~")[[1]][2], sep="")
  cat("\nFixed effects:", Ffix, "\n")
  print(coef(x), digits = digits, ...)
  cat("\n")
  cat("Random effects:", deparse(x$call$random), "\n")
  cat("Estimated variance-covariance matrix:\n")
  print(x$Sigma, digits = digits, ...)
#  cat("Variance-Covariance:\n")
#	print(x$var, ...)
  cat("\nNumber of Observations:", x$n)
  cat("\nNumber of Groups: ", nrow(x$bhat))
  cat("\n\n")
}

print.summary.phmm <- 
 function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nProportional Hazards Mixed-Effects Model fit by MCMC-EM\n")
  cat("  Model:", deparse(x$call$formula),"\n")
  cat("  Data:", deparse( x$call$data ), "\n")
  cat("  Log-likelihood:\n")
	print(x$loglik[!is.null(x$loglik)], digits = digits, ...)
  Ffix <- paste("~", strsplit(deparse(x$call$formula), "~")[[1]][2], sep="")
  cat("\nFixed effects:", Ffix, "\n")
  print(x$coef, digits = digits, ...)
  cat("\n")
  cat("Random effects:", deparse(x$call$random), "\n")
  cat("Estimated variance-covariance matrix:\n")
  print(x$Sigma, digits = digits, ...)
#  cat("Variance-Covariance:\n")
#	print(x$var, ...)
	cat("\nNumber of Observations:", x$n)
  cat("\nNumber of Groups: ", nrow(x$bhat))
  cat("\n\n")
}

summary.phmm <-
 function(object, ...)
{
	object$coefficients = cbind(Estimate = object$coef,
	   Std.Error = sqrt(ifelse(diag(object$var)<0, NA, diag(object$var)))[1:object$nfixed])
	class(object)<-"summary.phmm"
	return(object)
}

plot.phmm <-
 function(x, ...)
{
	x = as.data.frame(x$steps)
	colnames(x) = make.names(colnames(x))
	fm = paste(paste(colnames(x), collapse=' + '), "EM.Step", sep=" ~ ")
	x$EM.Step = as.numeric(rownames(x))
	xyplot(formula(fm), data=x, type="l", allow.multiple=TRUE, outer=TRUE, scales="free", ...)
}