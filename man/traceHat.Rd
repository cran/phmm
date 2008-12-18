\name{ traceHat}
\alias{ traceHat}
\alias{ traceHat.default}
\title{
Trace of the "hat" matrix from PHMM-MCEM fit
}
\description{
Compute trace of the "hat" matrix from PHMM-MCEM fit using pseudo poisson generalized linear mixed-effects model (GLMM).
}
\usage{
 traceHat(x, method="direct")
 traceHat.default(z, w, cluster, Sigma, fitted)
}
\arguments{
\item{x}{an object of class \code{phmm},}
\item{method}{acceptable values are \code{"direct"}, \code{"pseudoPois"}, \code{"HaLee"} or \code{coxph} (no random effects),}
\item{z}{matrix of fixed effects covariates,}
\item{w}{matrix of random effects covariates,}
\item{cluster}{integer valued cluster indices,}
\item{Sigma}{variance-covariance matrix of fit,}
\item{fitted}{linear predictors (fixed and random).}
}
\value{
The trace of the "hat" matrix which can be used as a measure of complexity of the model.
}

\references{
Breslow, NE, Clayton, DG. (1993). Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association, Vol. 88, No. 421, pp. 9-25.

Ha, ID, Lee, Y, MacKenzie, G. (2007). Model Selection for multi-component frailty models. Statistics in Medicine, Vol. 26, pp. 4790-4807.

Whitehead, J. (1980). Fitting Cox's Regression Model to Survival Data using GLIM. Journal of the Royal Statistical Society. Series C, Applied statistics, 29(3), 268-. }

\seealso{\code{\link{phmm}}, \code{\link{AIC.phmm}}}

\examples{\dontrun{
N <- 100
B <- 100
n <- 50
nclust <- 5
clusters <- rep(1:nclust,each=n/nclust)
beta0 <- c(1,2)
set.seed(13)
#generate phmm data set
Z <- cbind(Z1=sample(0:1,n,replace=TRUE),
           Z2=sample(0:1,n,replace=TRUE),
           Z3=sample(0:1,n,replace=TRUE))
b <- cbind(rep(rnorm(nclust),each=n/nclust),rep(rnorm(nclust),each=n/nclust))
Wb <- matrix(0,n,2)
for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
Wb <- apply(Wb,1,sum)
T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]\%*\%beta0-Wb)
C <- runif(n,0,1)
time <- ifelse(T<C,T,C)
event <- ifelse(T<=C,1,0)
mean(event)
phmmdata <- data.frame(Z)
phmmdata$cluster <- clusters
phmmdata$time <- time
phmmdata$event <- event

fit.phmm <- phmm(Surv(time, event)~Z1+Z2+cluster(cluster),
   ~-1+Z1+Z2, phmmdata, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
   NINIT = 10, MAXSTEP = 100, CONVERG=90)

# Same data can be fit with lmer,
# though the correlation structures are different.
poisphmmdata <- pseudoPoisPHMM(fit.phmm)

library(lme4)
fit.lmer <- lmer(m~-1+as.factor(time)+z1+z2+
  (-1+w1+w2|cluster)+offset(log(N)), 
  poisphmmdata, family=poisson)

fixef(fit.lmer)[c("z1","z2")]
fit.phmm$coef

VarCorr(fit.lmer)$cluster
fit.phmm$Sigma

logLik(fit.lmer)
fit.phmm$loglik

traceHat(fit.phmm)
}}
\keyword{survival}