pip3d = function(Vertices,Faces,Queries) {
    # Output is a NUMERIC vector.

    # Basic checks of Vertices input argument.
    vertDims = dim(Vertices)
    numColsV = vertDims[2] ;
    if ( numColsV != 3 ) {
        print("pip3d(): Number of columns in Vertices must be 3!")
        return(-3);
    }

    # Basic checks of Faces input argument.
    faceDims = dim(Faces)
    numColsF = faceDims[2] ;
    if ( numColsF != 3 ) {
        print("pip3d(): Number of columns in Faces must be 3!")
        return(-4);
    }

    if ( min(Faces) < 1 ) {
        print("pip3d(): Values in Faces must be greater than 0!")
        return(-5);
    }

    # Basic checks of Queries input argument.
    querDims = dim(Queries)
    numColsQ = querDims[2] ;
    if ( numColsQ != 3 ) {
        print("pip3d(): Number of columns in Queries must be 3!")
        return(-6);
    }

    # Invoke theptinpoly C++ code (for convenience, it has
    # the same name as this R function, but this wasn't necessary)
    Output =.C("pip3d",
        vts    = as.double(Vertices),
        nVts   = as.integer(nrow(Vertices)),
        fcs    = as.integer(Faces),
        nFcs   = as.integer(nrow(Faces)),
        qrs    = as.double(Queries),
        nQrs   = as.integer(nrow(Queries)),
        Result = as.vector(rep(0,nrow(Queries)),mode="integer")
    )

    # Return the Results data.
    return(Output$Result)
}
