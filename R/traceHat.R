traceHat.default <- function(z, w, cluster, Sigma, fitted){ 
	# z -- covs for fixed effects
	# w -- covs for random effects
	# Sigma -- varcov matrix for random effects
	# fitted -- linear part
	Sigma <- as.matrix(Sigma)
	nreff <- nrow(Sigma)
	
	z <- as.matrix(z)
	w <- as.matrix(w)	
	z <- as.matrix(z[order(cluster),])
	w <- as.matrix(w[order(cluster),])
	fitted <- fitted[order(cluster)]	
	cluster <- sort(cluster)
	nclust <- length(table(cluster))
	cluster <- rep(1:nclust,table(cluster))
	X <- z
	
	Z <- matrix(0,nrow(w),ncol(w)*nclust)
	for(i in 1:nclust){
		Z[cluster==i, (i-1)*ncol(w)+1:ncol(w)] <- w[cluster==i,]
		}
	D <- matrix(0,nclust*nreff,nclust*nreff)
	for(i in 1:nclust){
		D[(i-1)*nreff+1:nreff,(i-1)*nreff+1:nreff] <- Sigma
		}
	W <- diag(exp(fitted))
	U <- solve(D)

	if(nrow(Sigma)==1&Sigma[1,1]==0){
		return(sum(diag(solve(t(X)%*%W%*%X%*%t(X)%*%W%*%X))))
	}else{
		return(sum(diag(solve(
		rbind(cbind(t(X)%*%W%*%X,t(X)%*%W%*%Z),
			  cbind(t(Z)%*%W%*%X,t(Z)%*%W%*%Z+U)))%*%
		rbind(cbind(t(X)%*%W%*%X,t(X)%*%W%*%Z),
			  cbind(t(Z)%*%W%*%X,t(Z)%*%W%*%Z)))))
	}
}

traceHat.direct <- function(x, method="direct"){ 
	# use Ha, Lee, MacKenzie 2007 framework
	# to compute traceHat directly
	# z -- covs for fixed effects
	# w -- covs for random effects
	# time -- min(T,C)
	# Sigma -- varcov matrix for random effects
	# fitted -- linear part
	time = x$X; delta = x$delta; z = x$Z
	w = x$W; b = x$bhat.long; Sigma=as.matrix(x$Sigma)
	fitted = x$linear.predictors
	eventtimes = unique(sort(x$X[x$delta==1]))
	
	cluster = as.numeric(as.character(x$cluster))
	nclust = length(unique(cluster))
	
	xx <- cbind(ID = 1:length(time), time = time, delta=delta)
	Lambda=cumsum(x$lambda)[!duplicated(sort(x$X),fromLast=TRUE)]
	lambda=Lambda-c(0,Lambda[1:(length(Lambda)-1)])
	Lambda <- cbind(Lambda=Lambda,
					lambda=lambda,
	                time=unique(sort(x$X)))
	xx <- merge(xx, Lambda, by="time", all.x=TRUE)
	Lambda <- xx[order(xx$ID),"Lambda"]
	lambda <- xx[order(xx$ID),"lambda"]

	xxx <- xx[xx$delta==1, c("time", "lambda")]
	ulambda <- xxx[!duplicated(xxx$time), "lambda"]

	time <- time[order(cluster)]
	delta <- delta[order(cluster)]
	if(ncol(z)==1){ z <- as.matrix(z[order(cluster)])
		}else{ z <- as.matrix(z[order(cluster),]) }
	if(ncol(w)==1){ 
		w <- as.matrix(w[order(cluster)])
		b <- as.matrix(b[order(cluster)])
		bhat <- b[!duplicated(sort(cluster))]
		}else{ 
			w <- as.matrix(w[order(cluster),]) 
			b <- as.matrix(b[order(cluster),])
			bhat <- b[!duplicated(sort(cluster)),]
			}
	
	cluster <- sort(cluster)
	nclust <- length(table(cluster))
	cluster <- rep(1:nclust,table(cluster))

	X <- z
	Z <- matrix(0,x$n,x$nrandom*nclust)
	if(x$nrandom==1){
	for(i in unique(cluster)){
		Z[cluster==i, (i-1)*x$nrandom+1] <- w[cluster==i]
		}
	}else{
	for(i in unique(cluster)){
		Z[cluster==i, (i-1)*x$nrandom+1:x$nrandom] <- w[cluster==i,]
		}
	}
	D <- matrix(0,nclust*x$nrandom,nclust*x$nrandom)
	for(i in 1:nclust){
		D[(i-1)*x$nrandom+1:x$nrandom,(i-1)*x$nrandom+1:x$nrandom] <- Sigma
		}

	W3 <- diag(as.vector(exp(fitted)))
	B <- diag(Lambda)
	W1 <- W3%*%B
	
	M <- outer(time,eventtimes,FUN=function(x,y) ifelse(x>=y,1,0))

	dk <- rep(0,length(eventtimes))
	for(k in 1:length(eventtimes))
		dk[k] <- sum(delta[time==eventtimes[k]])
	C <- diag(dk/(ulambda^2))
	W2 <- (W3%*%M)%*%solve(C)%*%t(W3%*%M)
	W <- W1-W2
	
	if(nrow(Sigma)==1&Sigma[1,1]==0){ 
		return(sum(diag(solve(t(X)%*%W%*%X%*%t(X)%*%W%*%X))))
	}else{
		U <- solve(D)
		J <- U%*%as.vector(t(bhat))%*%t(as.vector(t(bhat)))%*%U
		if(method=="direct")
			return(sum(diag(solve(
			rbind(cbind(t(X)%*%W%*%X,t(X)%*%W%*%Z),
			  cbind(t(Z)%*%W%*%X,t(Z)%*%W%*%Z+U)))%*%
			rbind(cbind(t(X)%*%W%*%X,t(X)%*%W%*%Z),
			  cbind(t(Z)%*%W%*%X,t(Z)%*%W%*%Z+J)))))
		else if(method=="HaLee")
			return(sum(diag(solve(
			rbind(cbind(t(X)%*%W%*%X,t(X)%*%W%*%Z),
			  cbind(t(Z)%*%W%*%X,t(Z)%*%W%*%Z+U)))%*%
			rbind(cbind(t(X)%*%W%*%X,t(X)%*%W%*%Z),
			  cbind(t(Z)%*%W%*%X,t(Z)%*%W%*%Z)))))
	}
}
	
traceHat <- function(x, method="direct"){ 
	if(! method %in% c("direct","pseudoPois","HaLee")) 
		stop("Undefined traceHat method.")
	if(method=="pseudoPois"){
	xx <- pseudoPoisPHMM(x)
	
	xx$ID <- 1:nrow(xx)
	Lambda=cumsum(x$lambda)[!duplicated(sort(x$X),fromLast=TRUE)]
	Lambda <- cbind(Lambda=Lambda,
	                time=unique(sort(x$X)))
	xx <- merge(xx, Lambda, by="time", all.x=TRUE)
	xx <- xx[order(xx$ID),]
	
	neventtimes <- length(unique(x$X[x$delta==1]))
	
	return(traceHat.default(xx[,c(paste("t", 1:neventtimes, sep=''),
						    paste("z", 1:x$nfixed, sep=''))],
							xx[,paste("w", 1:x$nrandom, sep='')], 
							xx$cluster, x$Sigma, 
							xx$linear.predictors + xx$Lambda + log(xx$N)))
	}else if(method=="direct") return(traceHat.direct(x, method="direct")) else
	if(method=="HaLee") return(traceHat.direct(x, method="HaLee"))
}
