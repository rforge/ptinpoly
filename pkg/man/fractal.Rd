\name{fractal}
\alias{fractal}
\docType{data}
\title{
Sample Data: Fractal Polygon
}
\description{
Sample data defining an enclosed fractal 2D polygon.
}
\usage{data(fractal)}
\examples{
# Load polygon.
data(fractal)

# Plot the polygon.
plot(rbind(fractal,fractal[1,]),type="l")

# Generate 3333 random test points.
set.seed(1902)
n       <- 3333
x1      <- rnorm(n) ; x2 <- rnorm(n)
X       <- cbind(x1,x2)
queries <- as.matrix(X)

# Check whether test points are contained in the polygon.
# Most of these points will lie outside the polygon.
containment = pip2d(fractal,queries);
}
\keyword{datasets}
