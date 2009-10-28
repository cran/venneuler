\name{venneuler}
\alias{venneuler}
\title{
Calculates Venn and Euler Diagram
}
\description{
\code{venneuler} calculates a Venn diagram from a set specification.
}
\usage{
venneuler(combinations, weights, ...)
}
\arguments{
  \item{combinations}{
This can be one of: 
\itemize{
\item a character vector (specifies disjoint class combinations as
class names separated by the ampersand \code{&} character --
e.g. \code{c("A","B","A&B")})
\item a named numeric vector (names specify class combinations and
values specify weights -- e.g. \code{c(A=1, B=2, `A&B`=0.5)})
\item a character matrix of two columns (specifies mapping of
elements to sets -- elements in the first column and set names in the
second column, \code{weights} argument is ignored)
\item a logical or numeric matrix whose columns represent sets and
co-occurrence is defined by non-zero (rep. \code{TRUE}) values in rows
(weight for a row being 1 for logical matrices or the row sum for
numeric matrices).}
For convenience data frames can be passed instead of matrices and they
will be coerced using \code{as.matrix()}.
}
  \item{weights}{
If \code{combinations} is a character vector then this argument
specifies the associated weights. It is ignored in all other cases.
}
  \item{\dots}{
Additional arguments (currently unused).
}
}
%%\details{
%%}
\value{
An object of the class \code{VennDiagram} with following components:
\item{centers}{centers of the circles (columns are \code{x} and \code{y} coordinates)}
\item{diameters}{diameters of the circles}
\item{colors}{colors of the circles as values between 0 and 1}
\item{labels}{labels of the circles}
\item{residuals}{residuals (percentage difference between input intersection area and fitted intersection area)}
\item{stress}{stress value for solution}
\item{stress01}{.01 critical value for stress based on random data}
\item{stress05}{.05 critical value for stress based on random data}
}
%%\references{
%%}
\author{Lee Wilkinson <leland.wilkinson@gmail.com>, R package: Simon Urbanek <simon.urbanek@r-project.org>}
%%\note{
%%}
%%\seealso{
%%}
\examples{
vd <- venneuler(c(A=0.3, B=0.3, C=1.1, "A&B"=0.1, "A&C"=0.2, "B&C"=0.1 ,"A&B&C"=0.1))
plot(vd)
# same as c(A=1, `A&B&C`=1, C=1)
m <- data.frame(elements=c("1","2","2","2","3"), sets=c("A","A","B","C","C"))
v <- venneuler(m)
plot(v)
m <- as.matrix(data.frame(A=c(1.5, 0.2, 0.4, 0, 0),
                          B=c(0  , 0.2, 0  , 1, 0),
                          C=c(0  , 0  , 0.3, 0, 1)))
# without weights
v <- venneuler(m > 0)
plot(v)
# with weights
v <- venneuler(m)
plot(v)
}
\keyword{multivariate}
\keyword{hplot}
