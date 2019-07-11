/*
// mxCircleStats.c
// Mex routine computes the mean color of pixel regions defined by circles in the input image
*/

#include "mex.h"
#include <math.h>

#define MSG_ID "mxCircleStats:arg"

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
    mwSize nrows, ncols, ntri, ndims;
    const mwSize *dims;
    double *CC_ptr, *RCC_ptr;
    double *mean_pixel_ptr, *pixel_count_ptr;
    double *image_ptr;
    double cx, cy, radius, radius2, dx, dy;
    int min_row, max_row, min_col, max_col, i, j, row, col, count, index;
    const mxArray *image;
    
    /* function signature:
    // [mean_pixel, pixel_count] = mxCircleStats (CC, RCC, image)
    // CC is an ntri x 2 array of doubles
    // RCC is an ntri x 1 array of doubles
    // image is an nrows x ncols x ndims double array - ndims image channels
    // mean_pixel will be ntri x ndims array of doubles
    // pixel_count will be an ntri x 1 array of doubles */
    
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("mxCircleStats:nrhs", "Three inputs required CC, RCC and image");
    }
    
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("mxCircleStats:nlhs", "Two Outputs required.");
    }
    
    CC_ptr = checkArgDoubleArray(prhs, 0, "First argument CC must be an ntri x 2 array", 0, 2);
    ntri = mxGetM(prhs[0]);
    
    RCC_ptr = checkArgDoubleArray(prhs, 1, "Second argument RCC must be an ntri x 1 array", ntri, 1);
    
    
    /* Check the third argument */
    image = prhs[2];
    
    if (!mxIsClass(image, "double")) {
        mexErrMsgIdAndTxt("mxCircleStats:image", "image argument needs to be of class double");
    }
    
    dims = mxGetDimensions(image);
    
    if (mxGetNumberOfDimensions(image) > 3) {
        mexErrMsgIdAndTxt("mxLabelConnectedComponents:image", "image should be an nrows x ncols x ndims array");
    }
    
    nrows = dims[0];
    ncols = dims[1];
    
    if (mxGetNumberOfDimensions(image) == 3)
        ndims = dims[2];
    else
        ndims = 1;
    
    image_ptr = mxGetData(image);
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(ntri, ndims, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(ntri, 1, mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    mean_pixel_ptr  = mxGetPr(plhs[0]);
    pixel_count_ptr = mxGetPr(plhs[1]);
    
    /* Initialize the means to zero */
    for (i=0; i < ntri; ++i)
        for (j=0; j < ndims; ++j)
            mean_pixel_ptr[i + (j*ntri)] = 0.0;
    
    /* Main computational Loop */
    for (i=0; i < ntri; ++i) {

        /* Get the center and radius
        // Note the -1 to account for 1 based indexing */
        cx = CC_ptr[i] - 1;
        cy = CC_ptr[i+ntri] - 1;
        radius = RCC_ptr[i];
        
        
        /* Find the bounds in the image to check */
        min_col = (int)(cx - radius);
        if (min_col < 0) min_col = 0;
        
        max_col = (int)(cx + radius);
        if (max_col >= ncols) max_col = ncols-1;
        
        min_row = (int)(cy - radius);
        if (min_row < 0) min_row = 0;
        
        max_row = (int)(cy + radius);
        if (max_row >= nrows) max_row = nrows-1;
        
        
        /* Note the 0.25 fudge factor which is there to ensure that the 3 corners of the
        // triangle are included */
        radius2 = radius*radius + 0.25;
        
        /* Initialize */
        count = 0;
        
        /* Loop over the area summing colors */
        for (row = min_row; row <= max_row; ++row) {
            for (col = min_col; col <= max_col; ++col) {
                dx = col - cx;
                dy = row - cy;
                
                if ( (dx*dx + dy*dy) < radius2 ) {
                    index = row + (col*nrows);
                    
                    for (j=0; j < ndims; ++j)
                        mean_pixel_ptr[i+(j*ntri)] += image_ptr[index + j*(nrows*ncols)];
                    
                    ++count;
                }
            }
        }
        
        pixel_count_ptr[i] = count;
        
        /* We should have at least 3 pixels in each triangle - one for each vertex. */
        if (count < 3) {
            mexErrMsgIdAndTxt("mxCircleStats:count", "Found a triangle with fewer than 3 pixels");
        }

        /* Divide by count. */
        for (j=0; j < ndims; ++j) mean_pixel_ptr[i+(j*ntri)] /= count;
        
    }
}
