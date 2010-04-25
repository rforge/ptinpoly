\name{pip3d}
\alias{pip3d}
\title{
Test for Point Containment in 3D Polyhedron
}
\description{
Tests whether points are contained within a three-dimensional polyhedron.
}
\usage{pip3d(Vertices,Faces,Queries)}
\arguments{
    \item{Vertices}{N by 3 matrix containing the XYZ coordinates of N vertices of the polyhedron}
    \item{Faces}{M by 3 matrix containing the indices of the three vertices defining the M triangular faces of the polyhedron}
    \item{Queries}{P by 3 matrix containing the XYZ coordinates of P points
	to be tested for containment in the polyhedron defined by 'Vertices' and 'Faces'}
}
\details{
The values in the \code{Faces} matrix must be integers with values running from 1 to N,
where N is the number of vertices.  A value of '1' in this matrix, for example,
represents the 1st vertex, i.e., the vertex defined by the first row in the
matrix \code{Vertices}.
}
\value{
Returns a vector containing P values, one for each of the P points listed in
the \code{Queries} matrix.

'1' indicates that the point is contained in the polyhedron.

'0' indicates that the point lies exactly on the surface of the polyhedron.

'-1' indicates that the point lies outside the polyhedron.
}
\note{
The polyhedron defined by \code{Vertices} and \code{Faces} \emph{must} be "non-leaky";
i.e., it must define an "inside" versus "outside" and must not contain any holes.
Running this function on a defective polyhedron \emph{will} cause R to crash.

For an example of external software that could potentially be used to fix defective polyhedra,
see, e.g., PolyMender (\url{http://www1.cse.wustl.edu/~taoju/code/polymender.htm}).
}
\references{
W.P. Horn and D.L. Taylor, \emph{A theorem to determine the spatial containment of a point in a planar polygon}, Computer Vision, Graphics and Image Processing, vol. 45, pp. 106-116,1989.

S. Nordbeck and B. Rysedt, \emph{Computer cartography point-in-polygon programs}, BIT, vol. 7, pp. 39-64, 1967.

J.A. Baerentzen and H. Aanaes, \emph{Signed distance computation using the angle weighted pseudo-normal}, IEEE Trans. Visualization and Computer Graphics, vol. 11, no. 3, pp. 243-253, May/June 2005.

J. Liu, Y.Q. Chen, J.M. Maisog, G. Luta, \emph{A new point containment test algorithm for polyhedron composed of huge number of triangles}, paper in preparation.

\url{http://ptinpoly.pbworks.com/}
}
\examples{
#-------------------------------------------
# Simple Cube example.

# Load sample data defining a simple cube. 
data(verts)
data(faces)

# Also load sample data for five test points.
data(queries)

# Test whether each point in 'queries' is contained in
# the simple cube defined by 'verts' and 'faces'.
# This should return "1  0  0  0 -1".
containment = pip3d(verts,faces,queries);

#-------------------------------------------
# Torus example.

# Make a donut-shaped polyhedron.
# library(misc3d)
torus <- parametric3d(fx = function(u,v) (1+0.25*cos(v))*cos(u),
                      fy = function(u,v) (1+0.25*cos(v))*sin(u),
                      fz = function(u,v) 0.25*sin(v),
                      u = seq(0,2*pi,length.out=10),
                      v = seq(0,2*pi,length.out=10),
                      engine = "none", color="orange", alpha=0.25)

# If desired, this torus can be rendered for visualization, e.g.:
# library(geometry)
# library(rgl)
# library(misc3d)
drawScene.rgl(torus)

# Convert the torus to vertices-faces representation.
ve       <- misc3d:::t2ve(torus)
Vertices <- t(ve$vb)
Faces    <- t(ve$ib)

# Generate 3333 random test points.
set.seed(1902)
n       <- 3333
x1      <- rnorm(n) ; x2 <- rnorm(n) ; x3 <- rnorm(n)
X       <- cbind(x1,x2,x3)
Queries <- as.matrix(X)

# Check whether test points are contained in the torus.
# Most of these points will lie outside the torus.
containment = pip3d(Vertices,Faces,Queries);

#-------------------------------------------
# Complex polyhedron example.

# Load polyhedron.
data(triMesh)

# If desired, this complex polyhedron can be rendered for visualization, e.g.:
# drawScene.rgl(triMesh)

# Convert to vertices-faces representation.
ve       <- misc3d:::t2ve(triMesh)
Vertices <- t(ve$vb)
Faces    <- t(ve$ib)

# Generate 3333 random test points.
set.seed(1902)
n       <- 3333
x1      <- rnorm(n) ; x2 <- rnorm(n) ; x3 <- rnorm(n)
X       <- cbind(x1,x2,x3)
Queries <- as.matrix(X)

# Check whether test points are contained in the polyhedron.
# Most of these points will lie outside the polyhedron.
containment = pip3d(Vertices,Faces,Queries);

#-------------------------------------------
# A warning.
#
# If you remove one of the faces of the cube, the resulting cube is no 
# longer "non-leaky".  Running 'pip3d' on the resulting defective
# polyhedron will cause R to crash.
#
#    badcube     = faces[1:11,]
#    containment = pip3d(verts,badcube,queries); # crash
#
}
\keyword{methods}
