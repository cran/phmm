\name{phmm.design}
\alias{phmm.design}
\title{
PHMM Design
}
\description{
Internal function for extracting design matrix from call to \code{phmm}. Code adapted from \code{bayessurvreg1.help} function in the \code{bayesSurv} package.
}
\usage{
phmm.design(m, formula, random, data)
}
\arguments{
\item{m}{\code{match.call} from call to \code{phmm}.}
\item{formula}{formula component from call to \code{phmm}.}
\item{random}{random component from call to \code{phmm}.}
\item{data}{data component from call to \code{phmm}.}
}
\value{
  \item{n}{number of observations (in the case of bivariate data, this
    is a number of single observations, i.e. \eqn{2\times \mbox{sample
	size}}{2*sample size}) included in the dataset}
  \item{ncluster}{number of clusters included in the dataset. In the
    case of bivariate data this is equal to the number of bivariate
    observations. If there are no random effects included in the model
    and if the observations are not bivariate then \code{ncluster = n}}
  \item{nwithin}{a~vector of length equal to \code{ncluster} with
    numbers of observations within each cluster. In the case of
    bivariate observations this is a~vector filled with 2's, if there are
    no random effects and if the observations are not bivariate then
    this is a~vector filled with 1's}
  \item{nY}{number of columns in the response matrix \eqn{Y}{Y}. This is
    equal to 2 if there are no interval-censored observations and equal
    to 3 if there is at least one interval censored observation in the
    dataset}
  \item{nZ}{number of columns in the design matrix \eqn{Z}{Z}. Note that
    the matrix \eqn{Z}{Z} contains covariates for both fixed and random
    effects}
  \item{nfixed}{number of fixed effects involved in the model. Note that
    possible intercept is always removed from the model}
  \item{nrandom}{number of random effects in the model, possible random
    intercept included}
  \item{randomInt}{\code{TRUE}/\code{FALSE} indicating whether the
    random intercept is included in the model}
  \item{Y}{response matrix. Its last column is always equal to the
    status indicator (1 for exactly observed event times, 0 for
    right-censored observations, 2 for left-censored observations, 3 for
    interval-censored observations).}
  \item{Z}{design matrix containing covariates for fixed effects.}
  \item{W}{design matrix containing covariates for random effects.}
  \item{Yinit}{response matrix extracted from \code{formula} using \code{\link[stats]{model.extract}}}
  \item{Zinit}{design matrix extracted from \code{formula} using
    \code{\link[stats]{model.matrix}} function}
  \item{cluster}{a~vector of length \code{n} with identifications of
    clusters (as given by \code{cluster} in \code{formula})}
  \item{indb}{a~vector of length \code{nZ} identifying fixed and random
    effects. \code{indb[j] = -1} if the \eqn{j}{j}th column of matrix
    \eqn{Z}{Z} is a fixed effects. it is equal to \eqn{l}{l} if the
    \eqn{j}{j}th column of matrix \eqn{Z}{Z} corresponds to the
    \eqn{l}{l}th random effect (in C++ indexing)}
  \item{rnames.Z}{row names of \code{Zinit}}
  \item{names.random}{column names of the \eqn{Z}{Z} matrix corespning
    to the random effects. If there is the random intercept in the
    model, the first component of this vector is equal to "(Intercept)"}
  }  
\seealso{
\code{\link{phmm}}, \code{\link[bayesSurv]{bayessurvreg1.help}}.}
\references{
Arnost Komarek (2007). bayesSurv: Bayesian Survival Regression with Flexible Error and Random Effects Distributions. R package version 0.5-9. http://www.karlin.mff.cuni.cz/~komarek
}
\keyword{survival}