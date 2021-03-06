\name{comb}
\alias{comb}
\docType{data}
\title{
Sample Data: Comb Polygon
}
\description{
Sample data defining an enclosed comb 2D polygon.
}
\usage{data(comb)}
\examples{
# Load polygon.
data(comb)

# Plot the polygon.
plot(rbind(comb,comb[1,]),type="l")

# Generate 3333 random test points.
set.seed(1902)
n       <- 3333
x1      <- rnorm(n) ; x2 <- rnorm(n)
X       <- cbind(x1,x2)
queries <- as.matrix(X)

# Check whether test points are contained in the polygon.
# Most of these points will lie outside the polygon.
containment = pip2d(comb,queries);
}
\keyword{datasets}
