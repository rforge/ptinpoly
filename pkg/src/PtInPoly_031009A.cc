///////////////////////////////
// XDemo_main.cc
// C interaction file
//#include "XDemo.h"
#include "R.h" // R functions
#include "Rmath.h" // R math

//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Below are my own includes and preprocessor macros.

void PIP_jianfei_cpp(double *vertices, int *numV,
                     int    *faces,    int *numF,
                     double *query,    int *numQ,
                     int    *result);
                     
// Set USE_VISUAL_CPP if compiling with some system that
// doesn't have the rint() function, e.g. MS Visual Studio
#undef USE_VISUAL_CPP
// Otherwise, unset this macro variable.

#include <stdlib.h>
#define VERBOSE 0

#if !defined M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

typedef double xyzCoordinate;

#include <stddef.h>    // For size_t

//-------------------------------------------------------------------
//-------------------------------------------------------------------

// Functions Passed to C++ from R must be passed in C extern format
// All variables are passed to C by reference (pointers);
// All output of functions is "void" (adjustments made via reference change)
extern "C" {
    // Memory allocation for doubles.
    double *double_Calloc(int num_items);
    double **double_Calloc2(int num_items1, int num_items2);

    // For ran2() function from Numerical Recipes in C
    double ran2(long *idum);
    long idum; /* Random seed for ran2() */

    int numVertices;
    int numFaces;
    int numQueries;

#ifdef USE_VISUAL_CPP
    int IsEven(double inputVal);
    double newRint(double inputVal);
#endif

#if 0    
    void DemoAutoCor(double *RetV, int *pLwant, double *InputVec, int *pLengthInput) {
        double *AutoOutput = XDemo(InputVec, pLengthInput[0]);
        int ii;
        int MaxTake = pLwant[0];
        if ( (int) floor(pLengthInput[0] / 4) < MaxTake) {
            MaxTake = (int) floor(pLengthInput[0] / 4);
        }
        for (ii = 0; ii <MaxTake; ii++) {
            RetV[ii] =AutoOutput[ii];
        }
        Free(AutoOutput); // Free Memory created by function
        double UseLessNormal = .4 + rnorm(0.0, 1.0) * 2;
        // Completely Useless Generation of Normal
        Rprintf("DemoAutoCor:: Completely Useless Normal = %.3f\n", UseLessNormal);
        R_FlushConsole();
        R_ProcessEvents();
        return; // Return Nothing.
    }
#endif

    void JoeTest(double *RetV, int *pLwant, double *InputVec, int *pLengthInput) {
        int numElements = pLengthInput[0];
        int i;
        for (i=0; i<numElements; i++)
            Rprintf("%g\n",InputVec[i]);
        return;
    }

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

    double ran2(long *idum)
    {
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        double temp;

        if (*idum <= 0) {
            if (-(*idum) < 1) *idum=1;
            else *idum = -(*idum);
            idum2=(*idum);
            for (j=NTAB+7;j>=0;j--) {
                k=(*idum)/IQ1;
                *idum=IA1*(*idum-k*IQ1)-k*IR1;
                if (*idum < 0) *idum += IM1;
                if (j < NTAB) iv[j] = *idum;
            }
            iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
    }
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

    /* (C) Copr. 1986-92 Numerical Recipes Software ,12$*1(.}12&V. */

    #include <stdio.h>
    #include <math.h>
    #define EXIT_FAILURE 1
    #define X 0
    #define Y 1
    #define Z 2
    #define MAX_INT   2147483647 
    typedef enum { MY_FALSE, MY_TRUE } ;
    
    #define DIM 3                  /* Dimension of points */
    typedef int    tPointi[DIM];   /* Type integer point */
    typedef double tPointd[DIM];   /* Type double point */
    #define PMAX 10000             /* Max # of pts */
    int check = 0;
#if 0
    tPointd Vertices[PMAX];        /* All the points */
    tPointi Faces[PMAX];           /* Each triangle face is 3 indices */
    tPointd Box[PMAX][2];          /* Box around each face */
#else
    double *Vertices, *Box;
    int *Faces;
#endif

    /*---------------------------------------------------------------------
    Function prototypes.
    ---------------------------------------------------------------------*/
    char 	InPolyhedron( int F, tPointd q, tPointd bmin, tPointd bmax, int radius );
    char    SegPlaneInt( int *Triangle, tPointd q, tPointd r, tPointd p, int *m );
    int     PlaneCoeff( int *T, tPointd N, double *D );
    void    Assigndi( tPointd p, tPointi a );
    int     ReadVertices( void );
    int     ReadFaces( void );
//    void    NormalVec( tPointd q, tPointd b, tPointd c, tPointd N );
    void    NormalVec( double *q, double *b, double *c, tPointd N );
    double  DotVertices( double *q, tPointd d );
    double  DotTpointd( tPointd q, tPointd d );
    void    SubVec( tPointd q, tPointd b, tPointd c );
    char    InTri3D( int *T, int m, tPointd p );
    char    InTri2D( tPointd Tp[3], tPointd pp );
    int     AreaSign( tPointd q, tPointd b, tPointd c );
    char    SegTriInt( int *Triangle, tPointd q, tPointd r, tPointd p );
    char    InPlane( int *Triangle, int m, tPointd q, tPointd r, tPointd p);
    int     VolumeSign( tPointd a, tPointd b, tPointd c, tPointd d );
    char    SegTriCross( int *Triangle, tPointd q, tPointd r );
    int  	ComputeBox( int F, tPointd bmin, tPointd bmax );
    void 	RandomRay( tPointd ray, int radius );
    void 	AddVec( tPointd q, tPointd ray );
    int  	InBox( tPointd q, tPointd bmin, tPointd bmax );
    char 	BoxTest ( int n, tPointd a, tPointd b );
    void 	PrintPoint( tPointd q );
    int	irint( double x);

    /*
    This code is described in "Computational Geometry in C" (Second Edition),
    Chapter 7.  It is not written to be comprehensible without the
    explanation in that book.

    Compile:    gcc -o inhedron inhedron.c -lm (or simply: make)
    Run (e.g.): inhedron < i.8

    Written by Joseph O'Rourke, with contributions by Min Xu.
    Last modified: April 1998
    Questions to orourke@cs.smith.edu.
    --------------------------------------------------------------------
    This code is Copyright 1998 by Joseph O'Rourke.  It may be freely
    redistributed in its entirety provided that this copyright notice is
    not removed.
    --------------------------------------------------------------------
    */


    /*
    This function returns a char:
    'V': the query point a coincides with a Vertex of polyhedron P.
    'E': the query point a is in the relative interior of an Edge of polyhedron P.
    'F': the query point a is in the relative interior of a Face of polyhedron P.
    'i': the query point a is strictly interior to polyhedron P.
    'o': the query point a is strictly exterior to( or outside of) polyhedron P.
    */
    char InPolyhedron( int F, tPointd q, tPointd bmin, tPointd bmax, int radius )
    {
        tPointd r;  /* Ray endpoint. */
        tPointd p;  /* Intersection point; not used. */
        int f, k = 0, crossings = 0;
        char code = '?';

#if VERBOSE
        Rprintf("Entered InPolyhedron()...\n");
        R_FlushConsole();
        R_ProcessEvents();
#endif
        /* If query point is outside bounding box, finished. */
        if ( !InBox( q, bmin, bmax ) )
            return 'o';
#if VERBOSE
        Rprintf("In InPolyhedron(), returned from InBox()...\n");
        R_FlushConsole();
        R_ProcessEvents();
#endif

LOOP:
        while( k++ < F ) {
            crossings = 0;

            RandomRay( r, radius ); 
            AddVec( q, r );
#if VERBOSE
            Rprintf("Ray endpoint: (%d,%d,%d)\n", r[0],r[1],r[2] );
            R_FlushConsole();
            R_ProcessEvents();
#endif

            for ( f = 0; f < F; f++ ) {  /* Begin check each face */
                if ( BoxTest( f, q, r ) == '0' ) {
                    code = '0';
#if VERBOSE
                    Rprintf("BoxTest = 0!\n");
                    R_FlushConsole();
                    R_ProcessEvents();
#endif
                }
                else code = SegTriInt( Faces+f, q, r, p );
#if VERBOSE
                Rprintf( "Face = %d: BoxTest/SegTriInt returns %c\n\n", f, code );
                R_FlushConsole();
                R_ProcessEvents();
#endif

                /* If ray is degenerate, then goto outer while to generate another. */
                if ( code == 'p' || code == 'v' || code == 'e' ) {		 
#if VERBOSE
                    Rprintf("Degenerate ray\n");
                    R_FlushConsole();
                    R_ProcessEvents();
#endif
                    goto LOOP;
                }

                /* If ray hits face at interior point, increment crossings. */
                else if ( code == 'f' ) {
                    crossings++;
#if VERBOSE
                    Rprintf( "crossings = %d\n", crossings );
                    R_FlushConsole();
                    R_ProcessEvents();
#endif
                }

                /* If query endpoint q sits on a V/E/F, return that code. */
                else if ( code == 'V' || code == 'E' || code == 'F' )
                    return( code );

                /* If ray misses triangle, do nothing. */
                else if ( code == '0' )
                    ;

                else 
                    error("Error in InPolyhedron()!\n");

            } /* End check each face */

            /* No degeneracies encountered: ray is generic, so finished. */
            break;

        } /* End while loop */

#if VERBOSE
        Rprintf( "Crossings = %d\n", crossings );
        R_FlushConsole();
        R_ProcessEvents();
#endif
        /* q strictly interior to polyhedron iff an odd number of crossings. */
        if( ( crossings % 2 ) == 1 )
            return   'i';
        else return 'o';
    }

    int ComputeBox( int F, tPointd bmin, tPointd bmax )
    {
        int i, j; //   int i, j, k;
        double radius;

        for( i = 0; i < F; i++ )
            for( j = 0; j < DIM; j++ ) {
                if( Vertices[i+j*numVertices] < bmin[j] )
                    bmin[j] = Vertices[i+j*numVertices];
                if( Vertices[i+j*numVertices] > bmax[j] ) 
                    bmax[j] = Vertices[i+j*numVertices];
            }

            radius = sqrt( pow( (double)(bmax[X] - bmin[X]), 2.0 ) +
                pow( (double)(bmax[Y] - bmin[Y]), 2.0 ) +
                pow( (double)(bmax[Z] - bmin[Z]), 2.0 ) );
#if VERBOSE
            Rprintf("radius = %lf\n", radius);
            R_FlushConsole();
            R_ProcessEvents();
#endif
            return irint( radius +1 ) + 1;
    }

    /* Return a random ray endpoint */
    void RandomRay( tPointd ray, int radius )
    {
        double x, y, z, w, t;

#if VERBOSE
        Rprintf("Entered RandomRay()...\n");
        R_FlushConsole();
        R_ProcessEvents();
#endif

        /* Generate a random point on a sphere of radius 1. */
        /* the sphere is sliced at z, and a random point at angle t
        generated on the circle of intersection. */
        /*  z = 2.0 * (double) random() / MAX_INT - 1.0;
        t = 2.0 * M_PI * (double) random() / MAX_INT; */
        z = 2.0 * ran2(&idum) - 1.0;
        t = 2.0 * M_PI * ran2(&idum);
        w = sqrt( 1 - z*z );
        x = w * cos( t );
        y = w * sin( t );

        ray[X] = irint ( radius * x );
        ray[Y] = irint ( radius * y );
        ray[Z] = irint ( radius * z );

        /*printf( "RandomRay returns %6d %6d %6d\n", ray[X], ray[Y], ray[Z] );*/
    }

    void AddVec( tPointd q, tPointd ray )
    {
        int i;

#if VERBOSE
        Rprintf("Entered AddVec()...\n");
        R_FlushConsole();
        R_ProcessEvents();
#endif

        for( i = 0; i < DIM; i++ )
            ray[i] = q[i] + ray[i];
    }

    int InBox( tPointd q, tPointd bmin, tPointd bmax )
    {
        //  int i;
#if VERBOSE
        Rprintf("Entered InBox()...\n");
        R_FlushConsole();
        R_ProcessEvents();
#endif

        if( ( bmin[X] <= q[X] ) && ( q[X] <= bmax[X] ) &&
            ( bmin[Y] <= q[Y] ) && ( q[Y] <= bmax[Y] ) &&
            ( bmin[Z] <= q[Z] ) && ( q[Z] <= bmax[Z] ) )
            return MY_TRUE;
        return MY_FALSE;
    }

    /*---------------------------------------------------------------------
    'p': The segment lies wholly within the plane.
    'q': The q endpoint is on the plane (but not 'p').
    'r': The r endpoint is on the plane (but not 'p').
    '0': The segment lies strictly to one side or the other of the plane.
    '1': The segement intersects the plane, and 'p' does not hold.
    ---------------------------------------------------------------------*/
    char	SegPlaneInt( int *T, tPointd q, tPointd r, tPointd p, int *m)
    {
        tPointd N; double D;
        tPointd rq;
        double num, denom, t;
        int i;

#if VERBOSE
        Rprintf("Entered SegPlaneInt()...\n");
        R_FlushConsole();
        R_ProcessEvents();
#endif

        *m = PlaneCoeff( T, N, &D );
        /*printf("m=%d; plane=(%lf,%lf,%lf,%lf)\n", m, N[X],N[Y],N[Z],D);*/
        num = D - DotTpointd( q, N );
        SubVec( r, q, rq );
        denom = DotTpointd( rq, N );
        /*printf("SegPlaneInt: num=%lf, denom=%lf\n", num, denom );*/

        if ( denom == 0.0 ) {  /* Segment is parallel to plane. */
            if ( num == 0.0 )   /* q is on plane. */
                return 'p';
            else
                return '0';
        }
        else
            t = num / denom;
        /*printf("SegPlaneInt: t=%lf \n", t );*/

        for( i = 0; i < DIM; i++ )
            p[i] = q[i] + t * ( r[i] - q[i] );

        if ( (0.0 < t) && (t < 1.0) )
            return '1';
        else if ( num == 0.0 )   /* t == 0 */
            return 'q';
        else if ( num == denom ) /* t == 1 */
            return 'r';
        else return '0';
    }

    /*---------------------------------------------------------------------
    Computes N & D and returns index m of largest component.
    ---------------------------------------------------------------------*/
    int	PlaneCoeff( int *T, tPointd N, double *D )
    {
        int i;
        double t;              /* Temp storage */
        double biggest = 0.0;  /* Largest component of normal vector. */
        int m = 0;             /* Index of largest component. */

        NormalVec( Vertices+T[0*numFaces], Vertices+T[1*numFaces], Vertices+T[2*numFaces], N );
        /*printf("PlaneCoeff: N=(%lf,%lf,%lf)\n", N[X],N[Y],N[Z]);*/
        *D = DotVertices( Vertices+T[0*numFaces], N );

        /* Find the largest component of N. */
        for ( i = 0; i < DIM; i++ ) {
            t = fabs( N[i] );
            if ( t > biggest ) {
                biggest = t;
                m = i;
            }
        }
        return m;
    }

    /*---------------------------------------------------------------------
    a - b ==> c.
    ---------------------------------------------------------------------*/
    void    SubVec( tPointd a, tPointd b, tPointd c )
    {
        int i;

        for( i = 0; i < DIM; i++ )
            c[i] = a[i] - b[i];
    }

    /*---------------------------------------------------------------------
    Returns the dot product of the two input vectors.
    ---------------------------------------------------------------------*/
    double	DotVertices( double *a, tPointd b )
    {
        int i;
        double sum = 0.0;

        for( i = 0; i < DIM; i++ )
            sum += a[i*numVertices] * b[i];

        return  sum;
    }

    double	DotTpointd( tPointd a, tPointd b )
    {
        int i;
        double sum = 0.0;

        for( i = 0; i < DIM; i++ )
            sum += a[i] * b[i];

        return  sum;
    }

    /*---------------------------------------------------------------------
    Compute the cross product of (b-a)x(c-a) and place into N.
    ---------------------------------------------------------------------*/
    void	NormalVec( double *a, double *b, double *c, tPointd N )
    {
        N[X] = ( c[Z*numVertices] - a[Z*numVertices] ) * ( b[Y*numVertices] - a[Y*numVertices] ) -
            ( b[Z*numVertices] - a[Z*numVertices] ) * ( c[Y*numVertices] - a[Y*numVertices] );
        N[Y] = ( b[Z*numVertices] - a[Z*numVertices] ) * ( c[X*numVertices] - a[X*numVertices] ) -
            ( b[X*numVertices] - a[X*numVertices] ) * ( c[Z*numVertices] - a[Z*numVertices] );
        N[Z] = ( b[X*numVertices] - a[X*numVertices] ) * ( c[Y*numVertices] - a[Y*numVertices] ) -
            ( b[Y*numVertices] - a[Y*numVertices] ) * ( c[X*numVertices] - a[X*numVertices] );
    }

    /* Assumption: p lies in the plane containing T.
    Returns a char:
    'V': the query point p coincides with a Vertex of triangle T.
    'E': the query point p is in the relative interior of an Edge of triangle T.
    'F': the query point p is in the relative interior of a Face of triangle T.
    '0': the query point p does not intersect (misses) triangle T.
    */

    char 	InTri3D( int *T, int m, tPointd p )
    {
        int i;           /* Index for X,Y,Z           */
        int j;           /* Index for X,Y             */
        int k;           /* Index for triangle vertex */
        tPointd pp;      /* projected p */
        tPointd Tp[3];   /* projected T: three new vertices */

        /* Project out coordinate m in both p and the triangular face */
        j = 0;
        for ( i = 0; i < DIM; i++ ) {
            if ( i != m ) {    /* skip largest coordinate */
                pp[j] = p[i];
                for ( k = 0; k < 3; k++ )
                    Tp[k][j] = Vertices[T[k*numFaces]+i*numVertices];
                j++;
            }
        }
        return( InTri2D( Tp, pp ) );
    }

    char 	InTri2D( tPointd Tp[3], tPointd pp )
    {
        int area0, area1, area2;

        /* compute three AreaSign() values for pp w.r.t. each edge of the face in 2D */
        area0 = AreaSign( pp, Tp[0], Tp[1] );
        area1 = AreaSign( pp, Tp[1], Tp[2] );
        area2 = AreaSign( pp, Tp[2], Tp[0] );

#if VERBOSE
        Rprintf("area0=%d  area1=%d  area2=%d\n",area0,area1,area2);
        R_FlushConsole();
        R_ProcessEvents();
#endif

        if ( ( area0 == 0 ) && ( area1 > 0 ) && ( area2 > 0 ) ||
             ( area1 == 0 ) && ( area0 > 0 ) && ( area2 > 0 ) ||
             ( area2 == 0 ) && ( area0 > 0 ) && ( area1 > 0 ) ) 
            return 'E';

        if ( ( area0 == 0 ) && ( area1 < 0 ) && ( area2 < 0 ) ||
             ( area1 == 0 ) && ( area0 < 0 ) && ( area2 < 0 ) ||
             ( area2 == 0 ) && ( area0 < 0 ) && ( area1 < 0 ) )
            return 'E';                 

        if ( ( area0 >  0 ) && ( area1 > 0 ) && ( area2 > 0 ) ||
             ( area0 <  0 ) && ( area1 < 0 ) && ( area2 < 0 ) )
            return 'F';

        if ( ( area0 == 0 ) && ( area1 == 0 ) && ( area2 == 0 ) )
            error("Error in InTriD\n");

        if ( ( area0 == 0 ) && ( area1 == 0 ) ||
             ( area0 == 0 ) && ( area2 == 0 ) ||
             ( area1 == 0 ) && ( area2 == 0 ) )
            return 'V';

        else  
            return '0';  
    }

    int     AreaSign( tPointd a, tPointd b, tPointd c )  
    {
        double area2;

        area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
                ( c[0] - a[0] ) * (double)( b[1] - a[1] );

        /* The area should be an integer. */
        if      ( area2 >  0.5 ) return  1;
        else if ( area2 < -0.5 ) return -1;
        else                     return  0;
    }                            

    char    SegTriInt( int *T, tPointd q, tPointd r, tPointd p )
    {
        int code = '?';
        int m = -1;
#if VERBOSE
        Rprintf("Entered SegTriInt()...\n");
        R_FlushConsole();
        R_ProcessEvents();
#endif

        code = SegPlaneInt( T, q, r, p, &m );

#if VERBOSE
        Rprintf("SegPlaneInt code=%c, m=%d; p=(%lf,%lf,%lf)\n", code,m,p[X],p[Y],p[Z]);
        R_FlushConsole();
        R_ProcessEvents();
#endif

        if      ( code == '0')
            return '0';
        else if ( code == 'q')
            return InTri3D( T, m, q );
        else if ( code == 'r')
            return InTri3D( T, m, r );
        else if ( code == 'p' )
            return InPlane( T, m, q, r, p );
        else if ( code == '1' )
            return SegTriCross( T, q, r );
        else /* Error */
            return code;
    }

    char	InPlane( int *T, int m, tPointd q, tPointd r, tPointd p)
    {
        /* NOT IMPLEMENTED */
        return 'p';
    }

    /*---------------------------------------------------------------------
    The signed volumes of three tetrahedra are computed, determined
    by the segment qr, and each edge of the triangle.  
    Returns a char:
    'v': the open segment includes a vertex of T.
    'e': the open segment includes a point in the relative interior of an edge
    of T.
    'f': the open segment includes a point in the relative interior of a face
    of T.
    '0': the open segment does not intersect triangle T.
    ---------------------------------------------------------------------*/

    char SegTriCross( int *T, tPointd q, tPointd r )
    {
        int vol0, vol1, vol2;

#if VERBOSE
        Rprintf("Entered SegTriCross()...\n");
        R_FlushConsole();
        R_ProcessEvents();
#endif

        vol0 = VolumeSign( q, Vertices+T[0*numFaces], Vertices+T[1*numFaces], r ); 
        vol1 = VolumeSign( q, Vertices+T[1*numFaces], Vertices+T[2*numFaces], r ); 
        vol2 = VolumeSign( q, Vertices+T[2*numFaces], Vertices+T[0*numFaces], r );

#if VERBOSE
        Rprintf( "SegTriCross:  vol0 = %d; vol1 = %d; vol2 = %d\n",vol0, vol1, vol2 ); 
        R_FlushConsole();
        R_ProcessEvents();
#endif

        /* Same sign: segment intersects interior of triangle. */
        if ( ( ( vol0 > 0 ) && ( vol1 > 0 ) && ( vol2 > 0 ) ) || 
             ( ( vol0 < 0 ) && ( vol1 < 0 ) && ( vol2 < 0 ) ) )
            return 'f';

        /* Opposite sign: no intersection between segment and triangle */
        if ( ( ( vol0 > 0 ) || ( vol1 > 0 ) || ( vol2 > 0 ) ) &&
             ( ( vol0 < 0 ) || ( vol1 < 0 ) || ( vol2 < 0 ) ) )
            return '0';

        else if ( ( vol0 == 0 ) && ( vol1 == 0 ) && ( vol2 == 0 ) )
            error("Error 1 in SegTriCross\n");

        /* Two zeros: segment intersects vertex. */
        else if ( ( ( vol0 == 0 ) && ( vol1 == 0 ) ) || 
                  ( ( vol0 == 0 ) && ( vol2 == 0 ) ) || 
                  ( ( vol1 == 0 ) && ( vol2 == 0 ) ) )
            return 'v';

        /* One zero: segment intersects edge. */
        else if ( ( vol0 == 0 ) || ( vol1 == 0 ) || ( vol2 == 0 ) )
            return 'e';

        else
            error("Error 2 in SegTriCross\n");

        // Should NEVER get this far!  Return 'X' indicating a problem!
        return 'X';
    }

    int 	VolumeSign( tPointd a, double *b, double *c, tPointd d )
    { 
        double vol;
        double ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
        double bxdx, bydy, bzdz, cxdx, cydy, czdz;

        ax = a[X];
        ay = a[Y];
        az = a[Z];
        bx = b[X*numVertices];
        by = b[Y*numVertices];
        bz = b[Z*numVertices];
        cx = c[X*numVertices]; 
        cy = c[Y*numVertices];
        cz = c[Z*numVertices];
        dx = d[X];
        dy = d[Y];
        dz = d[Z];

        bxdx=bx-dx;
        bydy=by-dy;
        bzdz=bz-dz;
        cxdx=cx-dx;
        cydy=cy-dy;
        czdz=cz-dz;
        vol =   (az-dz) * (bxdx*cydy - bydy*cxdx)
            + (ay-dy) * (bzdz*cxdx - bxdx*czdz)
            + (ax-dx) * (bydy*czdz - bzdz*cydy);

        /* The volume should be an integer. */
        if      ( vol > 0.5 )   return  1;
        else if ( vol < -0.5 )  return -1;
        else                    return  0;
    }

    /*
    This function returns a char:
    '0': the segment [ab] does not intersect (completely misses) the 
    bounding box surrounding the n-th triangle T.  It lies
    strictly to one side of one of the six supporting planes.
    '?': status unknown: the segment may or may not intersect T.
    */
    char BoxTest ( int n, tPointd a, tPointd b )
    {
        int i; /* Coordinate index */
        double w;

#if VERBOSE
        Rprintf("Entered BoxTest()...\n");
        R_FlushConsole();
        R_ProcessEvents();
#endif

        for ( i=0; i < DIM; i++ ) {
            w = Box[ n*6 + 0*3 + i]; /* min: lower left */
//            w = Box[ n ][0][i]; /* min: lower left */
            if ( (a[i] < w) && (b[i] < w) ) return '0';
            w = Box[ n*6 + 1*3 + i]; /* max: upper right */
//            w = Box[ n ][1][i]; /* max: upper right */
            if ( (a[i] > w) && (b[i] > w) ) return '0';
        }
        return '?';
    }


    /* irint not available in some libraries, so... */
#ifdef USE_VISUAL_CPP

    // Determine whether a value is even or odd.  Needed for newRint().
    int IsEven(double inputVal) {
        if ( 2.0*(double)(((int)inputVal)/2) == inputVal )
            return 1;
        else
            return 0;
    }

    // Implementation of rint().
    // Algorithm from http://en.wikipedia.org/wiki/Rounding#Round-to-even_method
    double newRint(double inputVal) {
        /* Decide which is the last digit to keep.  */
        double a,b;
        b = modf(inputVal,&a);

        /* Input is exactly zero. */
        if ( inputVal == 0.0 )
            return 0.0;

        /* Positive numbers */
        else if ( inputVal > 0.0 ) {
            /* Increase it by 1 if the next digit is 6 or more, or a 5 followed by one or more non-zero digits.  */
            if ( b > 0.5 )
                return a + 1.0;

            /* Leave it the same if the next digit is 4 or less  */
            if ( b < 0.5 )
                return a;

            /* Otherwise, if all that follows the last digit is a 5 and possibly trailing zeroes;
            * then increase the rounded digit if it is currently odd; else, if it is already even, leave it alone.  */
            if ( IsEven(a) )
                return a;
            else
                return a + 1.0;

            /* Negative numbers */
        } else {
            /* DEcrease it by 1 if the next digit is 6 or more, or a 5 followed by one or more non-zero digits.  */
            if ( b < -0.5 )
                return a - 1.0;

            /* Leave it the same if the next digit is 4 or less  */
            if ( b > -0.5 )
                return a;

            /* Otherwise, if all that follows the last digit is a 5 and possibly trailing zeroes;
            * then DEcrease the rounded digit if it is currently odd; else, if it is already even, leave it alone.  */
            if ( IsEven(a) )
                return a;
            else
                return a - 1.0;
        }
    }

    int	irint( double x )
    {
        return (int) newRint( x );
    }

#else
    int	irint( double x )
    {
        return (int) rint( x );
    }
#endif

    void PIPorourke(double *vertices, int *numV,
                     int    *faces,    int *numF,
                     double *query,    int *numQ,
                     char   **result)              {
        numVertices = numV[0];
        numFaces    = numF[0];
        numQueries  = numQ[0];

        int i, j, k;
        double w;
        tPointd q, bmin, bmax;
        int radius;

        // Allocate memory for Box[].
        Box = double_Calloc(6*numFaces); // Store two XYZ coordinatess per Face

        /*
        srandom( (int) time( (long *) 0 ) ); 
        */
        idum = -1357; /* Random seed */
        ran2(&idum);   /* Initialize ran2() */

        // n = ReadVertices();
        Vertices = vertices;

        // F = ReadFaces();
        Faces = faces;
        for ( i = 0; i < numFaces; i++ ) {
            // Decrement vertex indices by 1, since O'Rourke's code starts indexing at 0.
            Faces[i+0*numFaces]--;
            Faces[i+1*numFaces]--;
            Faces[i+2*numFaces]--;

            /* Compute bounding box. */
            /* Initialize to first vertex. */
            for ( j=0; j < 3; j++ ) {
                Box[i*6 + 0*3 + j] = Vertices[ Faces[i+0*numFaces] + j*numVertices];
                Box[i*6 + 1*3 + j] = Vertices[ Faces[i+1*numFaces] + j*numVertices];
//                Box[i][0][j] = Vertices[ Faces[i][0] ][j];
//                Box[i][1][j] = Vertices[ Faces[i][1] ][j];
            }
            /* Check k=1,2 vertices of face. */
            for ( k=1; k < 3; k++ )
                for ( j=0; j < 3; j++ ) {
                    w = Vertices[ Faces[i+k*numFaces] + j*numVertices];
                    if ( w < Box[i*6 + 0*3 + j] ) Box[i*6 + 0*3 + j] = w;
                    if ( w > Box[i*6 + 1*3 + j] ) Box[i*6 + 1*3 + j] = w;
//                    w = Vertices[ Faces[i][k] ][j];
//                    if ( w < Box[i][0][j] ) Box[i][0][j] = w;
//                    if ( w > Box[i][1][j] ) Box[i][1][j] = w;
                }
        }

        /* Initialize the bounding box */
        for ( i = 0; i < DIM; i++ )
            bmin[i] = bmax[i] = Vertices[0+i*numVertices];
//            bmin[i] = bmax[i] = Vertices[0][i];
        radius = ComputeBox( numVertices, bmin, bmax );

        // Loop over queries, feed each to InPolyhedron()
        for (i=0; i<numQueries; i++) {
            q[0] = query[i+0*numQueries];
            q[1] = query[i+1*numQueries];
            q[2] = query[i+2*numQueries];
            result[i][0] = InPolyhedron( numFaces, q, bmin, bmax, radius );
            result[i][1] = '\0';
        }

        // Add back the 1 so that indexing in Faces[] starts at 1 again.
        // Release memory and return.
        for ( i = 0; i < numFaces; i++ ) {
            Faces[i+0*numFaces]++;
            Faces[i+1*numFaces]++;
            Faces[i+2*numFaces]++;
        }
        Free(Box);
        return;
    }

    void pip3d(double *vertices, int *numV,
               int    *faces,    int *numF,
               double *query,    int *numQ,
               int    *result)              {
//Rprintf("numV = %d\n",*numV);
//Rprintf("numF = %d\n",*numF);
//Rprintf("numQ = %d\n",*numQ);
        // Invoke wrapper function to C++ code.
        PIP_jianfei_cpp(vertices, numV,
                        faces,    numF,
                        query,    numQ,
                        result);

        return;
    }

    /* ************************************************************************** */
    double *double_Calloc(int num_items)
    {
        double *array=NULL;
        array=(double *)Calloc(num_items,double);
        if (array==NULL) {
            error("double_Calloc: unable to allocate %d bytes!\n",num_items*sizeof(double));

        }
        return array;
    }
}
