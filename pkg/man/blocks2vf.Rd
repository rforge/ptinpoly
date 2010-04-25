\name{blocks2vf}
\alias{blocks2vf}
\title{
Convertor from 3-Block Representation Representation to Vertices-Faces
}
\description{
Converts a polyhedron from the \dfn{three-block} representation to the \dfn{vertices-faces} representation.
}
\usage{blocks2vf(Block1,Block2,Block3)}
\arguments{
    \item{Block1}{M by 3 matrix containing the XYZ coordinates of vertex 1 of the M faces of the polyhedron}
    \item{Block2}{M by 3 matrix containing the XYZ coordinates of vertex 2 of the M faces of the polyhedron}
    \item{Block3}{M by 3 matrix containing the XYZ coordinates of vertex 3 of the M faces of the polyhedron}
}
\value{
Returns a list of two matrices.
The first is a N by 3 matrix containing the XYZ coordinates of the N vertices of the polyhedron.
The second M by 3 matrix containing indices of the vertices defining the M faces.  See the example below.
}
\details{
The values in the first output matrix can be floating point integers, representing the
XYZ coordinates of the vertices of the polyhedron.

The values in the second output matrix will be integers with values running from 1 to N,
where N is the number of vertices.  A value of '1' in this matrix, for example,
represents the 1st vertex, i.e., the vertex defined by the first row in the
matrix \code{Vertices}.  Each row in this matrix defines a triangular face in the polyhedron.

This function is the inverse of the \code{vf2blocks} function.
}
\note{
This function requires the \bold{misc3d} library.
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

# Now use blocks2vf to convert back to vertices-faces representation.
# 'verts2' and 'faces2' should encode the same polyhedron as the
# original 'verts' and 'faces', although perhaps in a different order.
# Note double square brackets.
vertsFaces = blocks2vf(block1,block2,block3)
verts2 = vertsFaces[[1]]
faces2 = vertsFaces[[2]]
}
\keyword{methods}
