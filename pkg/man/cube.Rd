\name{cube}
\alias{cube}
\alias{faces}
\alias{verts}
\alias{queries}
\docType{data}
\title{
Sample Data: Simple Cube
}
\description{
This is sample data that defines a simple cube: eight vertices,
and twelve triangles that make up the six faces.

Also included is an example matrix \code{queries} of five test points.
The first test point is contained within the cube, the second
through fourth test points lie exactly on the surface of the
cube, and the fifth test point lies outside the cube.
}
\usage{
data(verts)
data(faces)
data(queries)
}
\format{
\code{verts} is an 8 by 3 matrix containing the XYZ coordinates of the vertices of a simple cube.

\code{faces} is a 12 by 3 matrix containing the indices of the vertices defining the twelve
triangular faces making up the surface of the cube.

\code{queries} is a 5 by 3 matrix containing the XYZ coordinates of five test points to be
tested for containment within the cube.
}
\examples{
# Load sample data defining a simple cube. 
data(verts)
data(faces)

# Also load sample data for five test points.
data(queries)

# Test whether each point in 'queries' is contained in
# the simple cube defined by 'verts' and 'faces'.
# This should return "1  0  0  0 -1".
containment = pip3d(verts,faces,queries);
}
\keyword{datasets}
