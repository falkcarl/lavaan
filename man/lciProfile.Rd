\name{lciProfile}
\alias{lciProfile}
\title{Profile of difference test for likelihood-based confidence intervals}
\description{Computes values of the chi-square difference test for various null hypothesis values for a parameter or function of model parameters.}
\usage{
lciProfile(object, label, diff.method="default",grid=NULL)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{label}{A string corresponding to a label for the parameter estimate or quantity computed with the \code{":="} operator (e.g., see \code{\link{model.syntax}}). Valid labels for the current model can be inspected using \code{\link{parTable}}.}
\item{diff.method}{Passed as the \code{"method"} argument to \code{\link{lavTestLRT}} to compute difference test.}
\item{grid}{An optional vector indicating at what values the parameter (or function of model parameters) should be constrained for the purposes of computing the difference test. If NULL, plus and minus 2 standard deviations from the estimate will be used.}
}
\author{Carl F. Falk}
\details{This function is intended to be useful for diagnostic purposes when computing likelihood-based confidence intervals (e.g., using \code{\link{lci}}). For specified values along a grid, chi-square difference tests for the parameter (or function of model parameters) of interest (specified by \code{label}) will be computed. This function uses \code{\link{lavTestLRT}} to compute the difference test, and is intended for model comparisons where the difference test follows a chi-square reference distribution under the null hypothesis of no difference.
}
\value{
A list containing the following:

\item{grid}{Values of the grid that were used.}
\item{D}{Values of the difference test at each point along the grid}
}

\examples{
\dontrun{
# fit the Holzinger and Swineford (1939) example
# label for second loading (l2)
# label for factor variance (vVis)
# label for error variances (psi1 and psi2)
# and computed standardized loadings (lam1 and lam2)
HS.model <- ' visual  =~ x1 + l2*x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 
              
              visual~~ vVis*visual
              x1 ~~ psi1*x1
              x2 ~~ psi2*x2
              lam1 := 1/sqrt(1^2*vVis + psi1)*sqrt(vVis)
              lam2 := l2/sqrt(l2^2*vVis + psi2)*sqrt(vVis)'

# Normal theory
fit <- cfa(HS.model, data=HolzingerSwineford1939, se="standard")

test<-lciProfile(fit,"l2")
plot(test$grid,test$D,type="l",ylab="Difference Test",xlab="Parameter Value")

# "Robust"
fit.robust <- cfa(HS.model, data=HolzingerSwineford1939, se="standard",
  estimator="MLM", meanstructure=TRUE)

grid<-seq(0,.812)
testS2000<-lciProfile(fit,"l2",diff.method="satorra.2000",grid=grid)
plot(testS2000$grid,testS2000$D,type="l",ylab="Difference Test",xlab="Parameter Value")

testSB2001<-lciProfile(fit,"l2",diff.method="satorra.bentler.2001",grid=grid)
lines(testSB2001$grid,testSB2001$D,lty=2)

testSB2010<-lciProfile(fit,"l2",diff.method="satorra.bentler.2010",grid=grid)
lines(testSB2010$grid,testSB2010$D,lty=3)

}
}
