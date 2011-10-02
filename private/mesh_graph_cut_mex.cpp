// $Id: mesh_graph_cut_mex.cpp 171 2009-10-22 13:23:06Z gramfort $
// $LastChangedBy: gramfort $
// $LastChangedDate: 2009-10-22 15:23:06 +0200 (Thu, 22 Oct 2009) $
// $Revision: 171 $

#include "mex.h"
#include "maxflow/maxflow.h"

#include <limits>

void mexFunction( int  nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check number of inout arguments - must be 6 */
    if  ( (nrhs  != 6) & (nrhs  != 7)) {
        mexErrMsgTxt("MinCut-MaxFlow mex: Wrong number of input arguments");
    }

    int npoints,nedges,ntimes;
    float* space_edges;
    float* lambda_time;
    float* data_cost;

    /* check input width: must be 1 int element */
    npoints = mxGetScalar(prhs[0]);

    /* check input #edges: must 1 int element */
    nedges = mxGetScalar(prhs[1]);

    /* check input #times: must 1 int element */
    ntimes = mxGetScalar(prhs[2]);

    data_cost = (float*)mxGetData(prhs[3]);
    space_edges = (float*)mxGetData(prhs[4]);
    lambda_time = (float*)mxGetData(prhs[5]);

    float verbose = 1.0f; // true
    if(nrhs == 7) {
        verbose = ((float*)mxGetData(prhs[6]))[0];
    }

    if(verbose) {
        mexPrintf("v : %d :: e : %d :: t : %d\n",npoints,nedges,ntimes);
    }

    // ============================
    // = start graph construction =
    // ============================
    typedef Graph<float,float,float> GraphType;
    GraphType *g = new GraphType(npoints*ntimes,nedges*ntimes+npoints*(ntimes-1));

    for(int i = 0; i < npoints*ntimes; ++i) {
        g -> add_node();
    }

    // Handle T-Links
    for(int i = 0; i < npoints*ntimes; ++i) {
        // mexPrintf("v : %f :: t : %f\n",*(data_cost+2*i),*(data_cost+2*i+1));
        g -> add_tweights( i, *(data_cost+2*i), *(data_cost+2*i+1) );
    }

    // Handle temporal edges
    for(int t = 0; t < ntimes-1; ++t) {
        for(int i = 0; i < npoints; ++i) {
            g -> add_edge( t*npoints+i, (t+1)*npoints+i, *(lambda_time+i), *(lambda_time+i));
        }
    }

    // Handle spacial edges
    for(int t = 0; t < ntimes; ++t) {
        for(int i = 0; i < nedges; ++i) {
            g -> add_edge( t*npoints+(int)space_edges[3*i+0], t*npoints+(int)space_edges[3*i+1],
                           space_edges[3*i+2],space_edges[3*i+2]);
        }
    }

    plhs[0] = mxCreateDoubleMatrix(npoints,ntimes,mxREAL);
    double* labels = mxGetPr(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* flow = mxGetPr(plhs[1]);

    if(verbose) {
        mexPrintf("-- Starting MaxFlow\n");
    }
    *flow = (double) g -> maxflow();
    if(verbose) {
        mexPrintf("-- MaxFlow done !\n");
    }

    for(int i = 0; i < npoints*ntimes; ++i) {
        if (g->what_segment(i) == GraphType::SOURCE) {
            labels[i] = 1.0;
        } else {
            labels[i] = 0.0;
        }
    }

    delete g;

    return;
}
