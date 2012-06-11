\name{pip2d}
\alias{pip2d}
\title{
Test for Point Containment in 2D Polygon
}
\description{
Tests whether points are contained within a two-dimensional polygon.
}
\usage{pip2d(Vertices,Queries)}
\arguments{
    \item{Vertices}{N by 2 matrix containing the XY coordinates of N vertices of the polygon}
    \item{Queries}{P by 2 matrix containing the XY coordinates of P points
	to be tested for containment in the polygon defined by 'Vertices'}
}
\details{
The XY coordinates of the vertices are stored \emph{in order} in the
matrix \code{Vertices}. It is assumed that the last vertex listed in the matrix
is connected to the first vertex, so that the polygon does not have a ''hole''.
}
\value{
Returns a vector containing P values, one for each of the P points listed in
the \code{Queries} matrix.

'1' indicates that the point is contained in the polygon.

'0' indicates that the point lies exactly on the surface of the polygon.

'-1' indicates that the point lies outside the polygon.

'-3' (error) indicates that the \code{Vertices} matrix didn't have two columns

'-6' (error) indicates that the \code{Queries} matrix didn't have two columns

'-8' (error) indicates computational error not otherwise specified
}
\note{
The polygon defined by \code{Vertices} \emph{must} be "non-leaky";
i.e., it must define an "inside" versus "outside" and must not contain any holes.
}
\references{
W.P. Horn and D.L. Taylor, \emph{A theorem to determine the spatial containment of a point in a planar polygon}, Computer Vision, Graphics and Image Processing, vol. 45, pp. 106-116,1989.

S. Nordbeck and B. Rysedt, \emph{Computer cartography point-in-polygon programs}, BIT, vol. 7, pp. 39-64, 1967.

J.A. Baerentzen and H. Aanaes, \emph{Signed distance computation using the angle weighted pseudo-normal}, IEEE Trans. Visualization and Computer Graphics, vol. 11, no. 3, pp. 243-253, May/June 2005.

J. Liu, Y.Q. Chen, J.M. Maisog, G. Luta, \emph{A new point containment test algorithm for polygon composed of huge number of triangles}, Computer-Aided Design, Volume 42, Issue 12, December 2010, Pages 1143-1150.

\url{http://ptinpoly.pbworks.com/}
}
\examples{
#-------------------------------------------
# Load sample data defining a comb, spiral, and fractal. 
data(comb)
data(spiral)
data(fractal)

# Plot the comb, spiral, and fractal.
plot(rbind(comb,comb[1,]),type="l")
plot(rbind(spiral,spiral[1,]),type="l")
plot(rbind(fractal,fractal[1,]),type="l")

# Generate 3333 random test points.
set.seed(1902)
n       <- 3333
x1      <- rnorm(n) ; x2 <- rnorm(n)
X       <- cbind(x1,x2)
queries <- as.matrix(X)

# Check whether test points are contained in the comb, spiral, and factal.
# Most of these points will lie outside the polygons.
containment1 = pip2d(comb,queries);
containment2 = pip2d(spiral,queries);
containment3 = pip2d(fractal,queries);
}
\keyword{methods}
