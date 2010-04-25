\name{triMesh}
\alias{triMesh}
\docType{data}
\title{
Sample Data: Complex Polyhedron
}
\description{
Sample data defining a much more complex polyhedron than the simple cube example.
}
\usage{data(triMesh)}
\format{
\code{Triangles3D} object as defined by package \bold{misc3d}.
}
\note{
In order to use the \code{pip3d} function on
this data, it must first be converted to the \dfn{vertices-faces} format, as demonstrated
in the \bold{Examples} section below.
}
\examples{
# Load polyhedron.
data(triMesh)

# If desired, this complex polyhedron can be rendered for visualization, e.g.:
# library(geometry)
# library(rgl)
# library(misc3d)
drawScene.rgl(triMesh)

# Extract three blocks from triMesh.
block1 = triMesh$v1
block2 = triMesh$v2
block3 = triMesh$v3

# Convert from 3-Block representation to vertices-faces representation.
# Note double square brackets.
vertsFaces = blocks2vf(block1,block2,block3)
Vertices   = vertsFaces[[1]]
Faces      = vertsFaces[[2]]

# For a shorter way to convert this object to vertices-faces representation,
# see the Examples section in the documentation for the 'pip3d' function.

# Generate 3333 random test points.
set.seed(1902)
n       <- 3333
x1      <- rnorm(n) ; x2 <- rnorm(n) ; x3 <- rnorm(n)
X       <- cbind(x1,x2,x3)
Queries <- as.matrix(X)

# Check whether test points are contained in the polyhedron.
# Most of these points will lie outside the polyhedron.
containment = pip3d(Vertices,Faces,Queries);
}
\keyword{datasets}
