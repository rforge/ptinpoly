blocks2vf = function(Block1,Block2,Block3) {
    # Conversion from "block" representation to Vertices-Faces representation:
    # the vertices matrix contains XYZ coordinates of all vertices,
    # while the faces matrix contains integers indicating the three vertices
    # defining each face.

    # Basic checks of Block1 input argument.
    blkDims1 = dim(Block1)
    numFaces = blkDims1[1] ;
    numCols1 = blkDims1[2] ;
	if ( numCols1 != 3 ) {
	    print("blocks2vf(): Number of columns in Block1 must be 3!")
	    return(0);
	}

    # Basic checks of Block2 input argument.
    blkDims2 = dim(Block2)
	numRows2 = blkDims2[1] ;
    numCols2 = blkDims2[2] ;
	if ( numCols2 != 3 ) {
	    print("blocks2vf(): Number of columns in Block2 must be 3!")
	    return(0);
	}
	if ( numRows2 != numFaces ) {
	    print("blocks2vf(): Number of rows in Block2 must be the same as the number of rows in Block1!")
	    return(0);
	}

    # Basic checks of Block3 input argument.
    blkDims3 = dim(Block3)
    numRows3 = blkDims3[1] ;
    numCols3 = blkDims3[2] ;
	if ( numCols3 != 3 ) {
	    print("blocks2vf(): Number of columns in Block3 must be 3!")
	    return(0);
	}
	if ( numRows3 != numRows2 ) {
	    print("blocks2vf(): Number of rows in Block3 must be the same as the number of rows in Block2!")
	    return(0);
	}

	# Make sure misc3d is loaded
	if ( ! require("misc3d") ) {
	    print("blocks2vf(): package misc3d is required but missing!")
	    return(0);
	}

	# Use misc3d function to perform the conversion.
    ta    = makeTriangles(Block1,Block2,Block3)
    ve    = misc3d:::t2ve(ta)
    verts = t(ve$vb)
    faces = t(ve$ib)
    return (list(verts,faces))
}
