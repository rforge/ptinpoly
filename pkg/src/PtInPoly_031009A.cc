#include "R.h"     // R functions
//#include "Rmath.h" // R math

void PIP3D_jianfei_cpp(double *vertices, int *numV,
                       int    *faces,    int *numF,
                       double *query,    int *numQ,
                       int    *result);

void PIP2D_jianfei_cpp(double *vertices, int *numV,
                       double *query,    int *numQ,
                       int    *result);

//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Functions Passed to C++ from R must be passed in C extern format
// All variables are passed to C by reference (pointers);
// All output of functions is "void" (adjustments made via reference change)
//
extern "C" {
    void pip3d(double *vertices, int *numV,
               int    *faces,    int *numF,
               double *query,    int *numQ,
               int    *result)              {

        // Invoke wrapper function to C++ code.
        PIP3D_jianfei_cpp(vertices, numV,
                        faces,    numF,
                        query,    numQ,
                        result);

        return;
    }

    void pip2d(double *vertices, int *numV,
               double *query,    int *numQ,
               int    *result)              {

        // Invoke wrapper function to C++ code.
        PIP2D_jianfei_cpp(vertices, numV,
                          query,    numQ,
                          result);

        return;
    }
}
