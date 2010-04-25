\name{vf2blocks}
\alias{vf2blocks}
\title{
Convertor from Vertices-Faces Representation to 3-Block Representation
}
\description{
Converts a polyhedron from the \dfn{vertices-faces} representation to the \dfn{three-block} representation.
}
\usage{vf2blocks(Vertices,Faces)}
\arguments{
    \item{Vertices}{N by 3 matrix containing the XYZ coordinates of N vertices}
    \item{Faces}{M by 3 matrix containing indices of vertices defining M faces}
}
\details{
The values in the \code{Face} matrix must be integers with values running from 1 to N,
where N is the number of vertices.  A value of '1' in this matrix, for example,
represents the 1st vertex, i.e., the vertex defined by the first row in the
matrix \code{Vertices}.

This function is the inverse of the \code{blocks2vf} function.
}
\value{
Returns a list of three matrices.
The first is an M by 3 matrix containing the XYZ coordinates of vertex 1 of the M faces of the polyhedron.
The second and third are similarly M by 3 matrices, but contain the XYZ coordinates of
vertices 2 and 3 of the faces.  See the example below.
}
\examples{
# Load example data.
data(verts)
data(faces)

# Use vf2blocks to convert from vertices-faces representation to 3-block representation.
# Note double square brackets.
blocks = vf2blocks(verts,faces)
block1 = blocks[[1]]
block2 = blocks[[2]]
block3 = blocks[[3]]
}
\keyword{methods}
