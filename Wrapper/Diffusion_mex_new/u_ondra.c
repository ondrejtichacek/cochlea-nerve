/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include "mex.h"
#include "math.h"

void uu(double *J, double *E, double *V, mwSize n, mwSize nn, mwSize N, mwSize nt, mwSize nr)
{
    mwSize mm = fmin(n, N);
    // mwSize m = nt - mm; // v1
    mwSize m = N - mm; // v2

    mwSize i;
    mwSize j;
    mwSize k;
    mwSize ind;

    k = nn - mm;

    for (i=0; i<nr; i++) {
        V[i] = 0.0;
        //ind = i*nt + m;
        for (j=0; j<mm; j++) {
            //ind++;
            //ind = i*nt + (m+j);
            ind = i + (m+j)*nr;
            V[i] += J[k+j] * E[ind];
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *J;               /* 1xN input matrix */
    double *E;               /* 1xN input matrix */
    double n;                /* input scalar */
    double nn;               /* input scalar */
    double N;                /* input scalar */
    double nt;               /* input scalar */
    size_t nr;               /* size of matrix */
    double *V;               /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Seven inputs required.");
    }
    if(nlhs!=0) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Zero output required.");
    }
    /* make sure the first input argument is scalar */
//     if( !mxIsDouble(prhs[0]) || 
//          mxIsComplex(prhs[0]) ||
//          mxGetNumberOfElements(prhs[0])!=1 ) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
//     }
    
    /* make sure the second input argument is type double */
//     if( !mxIsDouble(prhs[1]) || 
//          mxIsComplex(prhs[1])) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
//     }
    
    /* check that number of rows in second input argument is 1 */
//     if(mxGetM(prhs[1])!=1) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
//     }

    /* create a pointer to the real data in the input matrix  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    J = mxGetDoubles(prhs[0]);
    E = mxGetDoubles(prhs[1]);
    V = mxGetDoubles(prhs[6]);
    #else
    J = mxGetPr(prhs[0]);
    E = mxGetPr(prhs[1]);
    V = mxGetPr(prhs[6]);
    #endif

    /* get the value of the scalar inputs  */
    n = mxGetScalar(prhs[2]);
    nn = mxGetScalar(prhs[3]);
    N = mxGetScalar(prhs[4]);
    nt = mxGetScalar(prhs[5]);

    /* get nrows of the input matrix */
    nr = mxGetM(prhs[1]);
    //nr = mxGetN(prhs[1]);

    /* create the output matrix */
    // plhs[0] = mxCreateDoubleMatrix((mwSize)nr,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    // V = mxGetDoubles(plhs[0]);
    #else
    // V = mxGetPr(plhs[0]);
    #endif

    /* call the computational routine */
    uu(J, E, V, (mwSize)n, (mwSize)nn, (mwSize)N, (mwSize)nt, (mwSize)nr);
}
