\name{plot.VennDiagram}
\alias{plot.VennDiagram}
\title{
\code{plot} method for Venn diagrams
}
\description{
  Plots the Venn diagram returned by \code{\link{venneuler}}.
}
\usage{
\method{plot}{VennDiagram}(x, col, col.fn = function(col) hcl(col * 360, 130, 60),
   alpha = 0.3, main = NULL, edges = 200, border = NA, col.txt = 1,
   cex = 1, lwd = 1, lty = 1, font = NULL, family = "", \dots)
}
\arguments{
  \item{x}{object of the class \code{VennDiagram} as returned from the
    \code{\link{venneuler}()} function.}
  \item{col}{optional, vector of colors (as accepted by the graphics
    system) to use. The colors are recycled if necessary (so
    passing a scalar will result in all circles having the same color).
    If not specified, colors are obtained by calling \code{col.fn} on the
    \code{colors} component of the \code{x} object.}
  \item{col.fn}{function taking one argument (numeric vector of values
    between 0 and 1), returning a vector of colors of the same
    length. It is not used if the \code{col} argument is specified.}
  \item{alpha}{numeric, value of the alpha channel(s) for the colors
    (hence their opacity). It will override any alpha channel
    information in the color specification, recycling as neeed. If set
    to \code{NA} then no aplha adjustment to the colors is performed.}
  \item{main}{passed to \code{title()}}
  \item{edges}{integer scalar, specifies the number of edges to use when
    drawing circles}
  \item{border}{color of the border for each circle (recycled) or
    \code{NULL} if no border is to be drawn}
  \item{lwd}{line width used to draw borders of the circles}
  \item{lty}{line type used to draw borders of the circles}
  \item{col.txt}{passed as \code{col} to \code{text()} for text labels in the circle centers}
  \item{cex}{passed to \code{text()} for text labels in the circle centers}
  \item{font}{passed to \code{text()} for text labels in the circle centers}
  \item{family}{passed to \code{text()} for text labels in the circle centers}
  \item{\dots}{any further arguments passed to \code{title()}}
}
%\details{
%}
\value{
  Returns \code{NULL} invisibly.
}
\author{
  Simon Urbanek
}
\seealso{
\code{\link{venneuler}}
}
\examples{
vd <- venneuler(c(A=0.3, B=0.3, C=1.1, "A&B"=0.1, "A&C"=0.2, "B&C"=0.1 ,"A&B&C"=0.1))
plot(vd, border=1, lwd = c(1,1,3), cex=2)
}
\keyword{multivariate}
\keyword{hplot}
