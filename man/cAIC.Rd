\name{ cAIC}
\alias{ cAIC}
\title{
Conditional Akaike Information Criterion for PHMM
}
\description{
Function calculating the conditional Akaike information criterion (Vaida & Blanchard 2005) for PHMM fitted model objects, according to the formula \eqn{-2*log-likelihood + k*rho}{-2*log-likelihood + k*rho}, where \eqn{rho}{rho} represents the "effective degrees of freedom" in the sense of Hodges and Sargent (2001). The function uses the log-likelihood conditional on the estimated random effects; and trace of the "hat matrix", using the generalized linear mixed model formulation of PHMM, to estimate \eqn{rho}{rho}. The default k = 2, conforms with the usual AIC.}
\usage{
 cAIC(object, ..., k = 2)
}
\arguments{
\item{object}{a fitted PHMM model object of class \code{phmm},}
\item{...}{optionally more fitted model objects.}
\item{k}{numeric, the penalty per parameter to be used; the default k = 2 conforms with the classical AIC.}
}
\value{
Returns a numeric value of the cAIC corresonding to the PHMM fit.
}

\references{
Vaida, F, & Blanchard, S. 2005. Conditional Akaike information for mixed-effects models. Biometrika, 92(2), 351-.

Breslow, NE, Clayton, DG. (1993). Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association, Vol. 88, No. 421, pp. 9-25.

Whitehead, J. (1980). Fitting Cox's Regression Model to Survival Data using GLIM. Journal of the Royal Statistical Society. Series C, Applied statistics, 29(3), 268-.

Hodges, JS, & Sargent, DJ. 2001. Counting degrees of freedom in hierarchical and other richly-parameterised models. Biometrika, 88(2), 367-.}

\seealso{\code{\link{phmm}}, \code{\link[stats]{AIC}}}

\examples{
\dontrun{
data(e1582)
e1582.fit <- phmm(Surv(time, event)~z1+z2+z3+z4+z5+cluster(cluster),
   ~-1+z1, e1582, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
   NINIT = 10, MAXSTEP = 50, CONVERG=40)
summary(e1582.fit)
cAIC(e1582)
}
}
\keyword{survival}