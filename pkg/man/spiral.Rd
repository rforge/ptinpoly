\name{spiral}
\alias{spiral}
\docType{data}
\title{
Sample Data: Spiral Polygon
}
\description{
Sample data defining an enclosed spiral 2D polygon.
}
\usage{data(spiral)}
\examples{
# Load polygon.
data(spiral)

# Plot the polygon.
plot(rbind(spiral,spiral[1,]),type="l")

# Generate 3333 random test points.
set.seed(1902)
n       <- 3333
x1      <- rnorm(n) ; x2 <- rnorm(n)
X       <- cbind(x1,x2)
queries <- as.matrix(X)

# Check whether test points are contained in the polygon.
# Most of these points will lie outside the polygon.
containment = pip2d(spiral,queries);
}
\keyword{datasets}
