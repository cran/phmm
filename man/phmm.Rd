\name{phmm}
\alias{phmm}
\title{Proportional Hazards Model with Mixed Effects}
\description{
Fits a proportional hazards regression model incorporating random effects.
The function implements an EM algorithm using Markov Chain Monte Carlo (MCMC) at the E-step as described in Vaida and Xu (2000).
}
\usage{
phmm(formula, random, data, subset, na.action, Sigma, varcov, NINIT,
	VARSTART, MAXSTEP, CONVERG, Gbs, Gbsvar, verbose)
}
\arguments{
  \item{formula}{model formula for the fixed part of the model, i.e. the
    part that specifies \eqn{\beta'x_{ij}}{beta'x[ij]}. See
    \code{\link[survival]{survreg}} for further details. Intercept is implicitely
    included in the model by estimation of the error distribution. As a
    consequence \code{-1} in the model formula does not have any effect
    on the model. The left-hand side of the \code{formula} must be an object
    created using \code{\link[survival]{Surv}}.
	If \code{random} is used then the formula must contain
    an identification of clusters in the form \code{cluster(id)}, where
    \code{id} is a name of the variable that determines clusters.
  }
  \item{random}{formula for the random part of the model, i.e. the
    part that specifies \eqn{w_{ij}'b_{i}}{w[ij]'b[i]}. If omitted,
    no random part is included in the model. E.g. to specify the model with a
    random intercept, say \code{random=~1}.
  }
  \item{data}{optional data frame in which to interpret the variables
    occuring in the formulas.}
  \item{subset}{subset of the observations to be used in the fit.}
  \item{na.action}{function to be used to handle any \code{NA}s in the
    data. The user is discouraged to change a default value
    \code{na.fail}.}
  \item{Sigma}{initial covariance matrix for the random effects. Defaults to
  "identity".}
  \item{varcov}{constraint on \code{Sigma}. Currently only \code{"diagonal"} is supported.}
  \item{NINIT}{number of starting values supplied to Adaptive Rejection Metropolis Sampling (ARMS) algorithm.}
  \item{VARSTART}{starting value of the variances of the random effects.}
  \item{MAXSTEP}{number of EM iterations.}
  \item{CONVERG}{iteration after which Gibbs sampling size changes from Gbs to Gbsvar.}
  \item{Gbs}{initial Gibbs sampling size (until CONVERG iterations).}
  \item{Gbsvar}{Gibbs sampling size after CONVERG iterations.}
  \item{verbose}{Set to \code{TRUE} to print EM steps.}
}
\value{
The function produces an object of class "phmm" consisting of: \item{steps}{ a matrix of estimates at each EM step;}
\item{bhat}{ empirical Bayes estimates of expectation of random effects;}
\item{sdbhat}{ empirical Bayes estimates of standard deviation of random effects;}
\item{coef}{ the final parameter estimates for the fixed effects;}
\item{var}{ the estimated variance-covariance matrix;}
\item{loglik}{ a vector of length four with the conditional log-likelihood and marginal log-likelihood estimated by Laplace approximation, reciprocal importance sampling, and bridge sampling (only implemented for \code{nreff} < 3);}
\item{lambda}{ the estimated baseline hazard;}
\item{Lambda}{ the estimated cumulative baseline hazard.}
}
\details{
The proportional hazards model with mixed effects is equipped to handle clustered survival data. The model generalizes the usual frailty model by allowing log-linearl multivariate random effects. The software can only handle random effects from a multivariate normal distribution. Maximum likelihood estimates of the regression parameters and variance components is gotten by EM algorithm, with Markov chain Monte Carlo (MCMC) used in the E-step.
}

\references{Gilks, W. R. and Wild, P. (1992) Adaptive rejection sampling for Gibbs sampling. Applied Statistics 41, pp 337-348.

Vaida, F. and Xu, R. 2000. "Proportional hazards model with random effects", \emph{Statistics in Medicine,} 19:3309-3324.

Xu, R, Gamst, A, Donohue, M, Vaida, F, \& Harrington, DP. 2006. Using Profile Likelihood for Semiparametric Model Selection with Application to Proportional Hazards Mixed Models. \emph{Harvard University Biostatistics Working Paper Series,} Working Paper 43.
}
\seealso{
\code{\link[survival]{survfit}}, \code{\link[survival]{Surv}}.}
\examples{
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
summary(fit.phmm)
}
\keyword{survival}
