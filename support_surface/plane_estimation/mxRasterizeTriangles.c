/*
 mxRasterizeTriangles.c
 Mex routine that can be used to rasterize the triangles of a DelaunayTriangulation into an image.
 This is useful in contexts where you have computed the Delaunay Triangulation of a set of image features.
 This rasterization avoids trimesh which can be very slow when you have thousands of edges. Note that you are producing
 a rasterized variant which will have it's own discretization issues.
*/

#include "mex.h"
#include <math.h>
#include <stdlib.h>

#define MSG_ID "mxRasterizeTriangles:arg"

/**************************************************************************************************/

/* Some relevant globals */

mwSize nrows, ncols, npts, ntri;
double *uv_ptr, *out_ptr, bg_color;
int *min_col_buffer, *max_col_buffer;

/**************************************************************************************************/

typedef struct {
    int u, v;
} point;

typedef struct {
    point P1, P2;
} edge;

/**************************************************************************************************/

point makePoint (int idx)
{
    point thePoint;
    
    /* Note the -1's to account for 0 based indexing */
    /* thePoint.u = round(uv_ptr[idx - 1]) - 1;
    thePoint.v = round(uv_ptr[idx+npts - 1]) - 1; */
    thePoint.u = (int)(uv_ptr[idx - 1] + 0.5) - 1;
    thePoint.v = (int)(uv_ptr[idx+npts - 1] + 0.5) - 1;
    
    return thePoint;
}

edge makeEdge (point P, point Q)
{
    edge theEdge;
    
    /* Make sure that the point with the smallest v coordinate is P1
     This simplifies subsequent rasterization */
    if (P.v < Q.v) {
        theEdge.P1 = P;
        theEdge.P2 = Q;
    } else {
        theEdge.P1 = Q;
        theEdge.P2 = P;
    }
    
    return theEdge;
}

int edgeHeight (edge *theEdge) {
    return (theEdge->P2.v - theEdge->P1.v);
}

void updateLineExtents (int row, int col)
{
    if ( (row >= 0) && (row < nrows) ) {
        if (col < min_col_buffer[row]) min_col_buffer[row] = col;
        if (col > max_col_buffer[row]) max_col_buffer[row] = col;
    }
        
}

void DrawSpans (int row1, int row2, int label)
{
    int row, col, min_col, max_col;
    
    for (row = row1; row <= row2; ++row) {
        min_col = min_col_buffer[row];
        max_col = max_col_buffer[row];
        
        if (min_col < 0) min_col = 0;
        if (max_col >= ncols) max_col = ncols-1;
        
        for (col = min_col; col <= max_col; ++col)
            out_ptr[row + (col*nrows)] = label;
    }
}

/* rasterizeEdge - this function fills in the min_col_buffer and max_col_buffer arrays which determine the
 limits of the triangle on each row */

void rasterizeEdge(edge *theEdge) {
    int x1 = theEdge->P1.u;
    int y1 = theEdge->P1.v;
    
    int x2 = theEdge->P2.u;
    int y2 = theEdge->P2.v;
    
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int length, i, x, y, idx;
    
    length = (dx > dy) ? dx : dy;
    
    if (length == 0) {
        
        updateLineExtents (y1, x1);
        
    } else {
        
        for (i = 0; i <= length; ++i) {
            
            /* Here we are trying to compute raster coordinates using integer arithmetic in such a way that
             we get the rounding correct: round (a/b) = floor ((a/b) + (1/2)) = floor ((2*a + b)/(2*b))
             NOTE that this only works if a and b are >= 0 */
            
            x = ( 2*(x1*i + x2*(length-i)) + length ) / (2*length);
            y = ( 2*(y1*i + y2*(length-i)) + length ) / (2*length);
            
            updateLineExtents (y, x);
            
        }
        
    }
}


void rasterizeTriangle (edge *edge1, edge *edge2, edge *edge3, int label)
{
    int height1 = edgeHeight(edge1);
    int height2 = edgeHeight(edge2);
    int height3 = edgeHeight(edge3);
    edge *tall_edge, *short_edge1, *short_edge2;
    int i, row1, row2;
    
    /* Find the tallest edge */
    if ( (height1 >= height2) && (height1 >= height3) ) {
        tall_edge   = edge1;
        short_edge1 = edge2;
        short_edge2 = edge3;
    } else if ( (height2 >= height1) && (height2 >= height3) ) {
        tall_edge   = edge2;
        short_edge1 = edge1;
        short_edge2 = edge3;
    } else {
        tall_edge   = edge3;
        short_edge1 = edge2;
        short_edge2 = edge1;
    }
    
    row1 = tall_edge->P1.v;
    row2 = tall_edge->P2.v;
    
    if (row1 < 0) row1 = 0;
    if (row2 >= nrows) row2 = nrows-1;
    
    /* Clear the min and max idx buffers */
    for (i = row1; i <= row2; ++i) {
        min_col_buffer[i] = ncols + 1;
        max_col_buffer[i] = -1;
    }
    
    /* render the edges to get the triangle extents on each row */
    rasterizeEdge (edge1);
    rasterizeEdge (edge2);
    rasterizeEdge (edge3);
    
    /* Go through and fill in the spans */
    DrawSpans (row1, row2, label);
}

/**************************************************************************************************/

/*
 * This function tests whether the specified argument is a 2D double array or scalar. You can
 * use the two arguments, nrows and ncols, to test the size of the array or set them to 0 if you don't care.
 */

double *checkArgDoubleArray(const mxArray *prhs[], int arg_index, char *error_msg, int nrows, int ncols) {
    const mxArray *in = prhs[arg_index];
    
    if (!mxIsDouble(in) || mxIsComplex(in) || (mxGetNumberOfDimensions(in) != 2)) {
        mexErrMsgIdAndTxt(MSG_ID, error_msg);
    }
    
    if ((nrows >= 1) && (mxGetM(in) != nrows)) {
        mexErrMsgIdAndTxt(MSG_ID, error_msg);
    }
    
    if ((ncols >= 1) && (mxGetN(in) != ncols)) {
        mexErrMsgIdAndTxt(MSG_ID, error_msg);
    }
    
    return (mxGetPr(in));
}


                
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) 
{
    const mxArray *tri;
    int i, *tri_ptr;
    point P1, P2, P3;
    edge edge1, edge2, edge3;
    
    
    /* function signature:
    // out = mxRasterizeTriangles (nrows, ncols, uv, tri, bg_color)
    // nrows, ncols - the size of the output image
    // uv - an n x 2 array of the point coordinates as doubles
    // tri - an m x 3 array of indices into the uv array
    // bg_color - the background color that the output is cleared to.
    
    // Note that we are assuming the coordinate system associated with an image starting at
    // the upper left corner of the array */
    
    /* check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("mxRasterizeTriangles:nrhs", "Five inputs required: nrows, ncols, uv, tri, bg_color");
    }
    
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("mxRasterizeTriangles:nlhs", "One output required.");
    }
    
    /* nrows = (mwSize) round(*checkArgDoubleArray(prhs, 0, "First argument nrows must be a scalar", 1, 1));
    ncols = (mwSize) round(*checkArgDoubleArray(prhs, 1, "Second argument ncols must be a scalar", 1, 1)); */
    nrows = (mwSize) (int)(*checkArgDoubleArray(prhs, 0, "First argument nrows must be a scalar", 1, 1) + 0.5 );
    ncols = (mwSize) (int)(*checkArgDoubleArray(prhs, 1, "Second argument ncols must be a scalar", 1, 1) + 0.5 );
    
    uv_ptr = checkArgDoubleArray(prhs, 2, "Third argument must be an npts x 2 array", 0, 2);
    npts = mxGetM (prhs[2]);
    
    /* Check the tri argument */
    tri = prhs[3];
    if (!mxIsClass(tri, "int32")) {
        mexErrMsgIdAndTxt("mxRasterizeTriangles:tri", "tri argument needs to be of class int32");
    }
    
    if (mxGetNumberOfDimensions(tri) != 2 || mxGetN(tri) != 3) {
        mexErrMsgIdAndTxt("mxRasterizeTriangles:tri", "tri should be an ntri x 3 array");
    }
    
    ntri = mxGetM(tri);
    tri_ptr = mxGetData(tri);
    
    bg_color = *checkArgDoubleArray(prhs, 4, "Fifth argument bg_color must be a scalar", 1, 1);
    
    
    if ( (nrows <= 0) || (ncols <= 0) ) mexErrMsgIdAndTxt(MSG_ID, "nrows <= 0 or ncols <= 0");
    
    /* We use these buffers to store the maximum and minimum triangle extents on each row */
    min_col_buffer = mxMalloc (nrows * sizeof(int));
    max_col_buffer = mxMalloc (nrows * sizeof(int));
    
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    out_ptr  = mxGetPr(plhs[0]);
    
    /* Clear the output to the specified background color */
    for (i=0; i < (nrows*ncols); ++i) out_ptr[i] = bg_color;
    
    /* Main computational Loop - Visit each triangle */
    for (i=0; i < ntri; ++i) {
        
        /* Extract the vertices */
        P1 = makePoint (tri_ptr[i]);
        P2 = makePoint (tri_ptr[i+ntri]);
        P3 = makePoint (tri_ptr[i+2*ntri]);
        
        /* Make the edges */
        edge1 = makeEdge (P1, P2);
        edge2 = makeEdge (P2, P3);
        edge3 = makeEdge (P3, P1);
        
        /* Draw the triangle */
        rasterizeTriangle (&edge1, &edge2, &edge3, i+1);
    
    }
    
    mxFree (min_col_buffer);
    mxFree (max_col_buffer);

}
