/*
// mxFindPlanes.c
// Mex routine computes the label of the pixel regions defined by triangles in the input image.
*/

#include "mex.h"
#include <math.h>
#include <stdlib.h>

#define MSG_ID "mxFindPlanes:arg"


/**************************************************************************************************/

typedef struct {
    int label, row, col;
    double u, v, w;
} pixel;

typedef struct {
    double nx, ny, nz;
} plane_normal;

typedef struct {
    plane_normal n;
    int count;  /* number of pixels in this plane */
    pixel *start_pixel;
    int merged_label;
} plane;


/** Globals ***************************************************************************************/

int nrows, ncols, min_points, ransac_trials, max_merged_planes, max_segments;
double cu, cv, focal_length, inlier_threshold, outlier_ratio, inlier_ratio, dotpthreshold;

struct {
    int nplanes;
    plane *thePlanes;
} Planes;

/*** Sorting functions *******************************************************/

static int compare_pixel_labels (const void *p1, const void *p2)
{
    int l1 = ((const pixel *) p1)->label;
    int l2 = ((const pixel *) p2)->label;
    
    if (l1 < l2) return +1;
    if (l1 > l2) return -1;
    return 0;
}


static int compare_plane_counts (const void *p1, const void *p2)
{
    int n1 = ((const plane *) p1)->count;
    int n2 = ((const plane *) p2)->count;
    
    if (n1 < n2) return +1;
    if (n1 > n2) return -1;
    return 0;
}

/**************************************************************************************************/
/* Routines to support least squares fitting */

struct {
    int count;
    double A[3][3];
    double b[3];
} LSQ;

void clearLSQ()
{
    int i, j;
    LSQ.count = 0;
    
    for (i=0; i < 3; ++i) {
        for (j=0; j < 3; ++j)
            LSQ.A[i][j] = 0.0;
        LSQ.b[i] = 0.0;
    }
}

void updateLSQ (pixel *p) 
{
    double u = p->u;
    double v = p->v;
    double w = p->w;
    
    LSQ.A[0][0] += u*u; LSQ.A[0][1] += u*v; LSQ.A[0][2] += u;
    LSQ.A[1][0] += v*u; LSQ.A[1][1] += v*v; LSQ.A[1][2] += v;
    LSQ.A[2][0] +=   u; LSQ.A[2][1] +=   v; LSQ.A[2][2] += 1;
    
    LSQ.b[0] += w*u;
    LSQ.b[1] += w*v;
    LSQ.b[2] += w;
    
    ++(LSQ.count);
}

void updateLSQplane (plane *thePlane)
{
    int i;
    
    for (i=0; i < thePlane->count; ++i)
        updateLSQ (thePlane->start_pixel + i);
}

plane_normal solveLSQ ()
{
    double det, C[3][3];
    plane_normal n;
    
    /* Compute matrix of cofactors transposed */
    C[0][0] =  (LSQ.A[1][1]*LSQ.A[2][2] - LSQ.A[1][2]*LSQ.A[2][1]);
    C[1][0] = -(LSQ.A[1][0]*LSQ.A[2][2] - LSQ.A[1][2]*LSQ.A[2][0]);
    C[2][0] =  (LSQ.A[1][0]*LSQ.A[2][1] - LSQ.A[1][1]*LSQ.A[2][0]);
    
    C[0][1] = -(LSQ.A[0][1]*LSQ.A[2][2] - LSQ.A[0][2]*LSQ.A[2][1]);
    C[1][1] =  (LSQ.A[0][0]*LSQ.A[2][2] - LSQ.A[0][2]*LSQ.A[2][0]);
    C[2][1] = -(LSQ.A[0][0]*LSQ.A[2][1] - LSQ.A[0][1]*LSQ.A[2][0]);
    
    C[0][2] =  (LSQ.A[0][1]*LSQ.A[1][2] - LSQ.A[0][2]*LSQ.A[1][1]);
    C[1][2] = -(LSQ.A[0][0]*LSQ.A[1][2] - LSQ.A[0][2]*LSQ.A[1][0]);
    C[2][2] =  (LSQ.A[0][0]*LSQ.A[1][1] - LSQ.A[0][1]*LSQ.A[1][0]);
    
    /* Compute the determinant using some entries from the cofactor matrix */
    det = LSQ.A[0][0]*C[0][0] + LSQ.A[0][1]*C[1][0] + LSQ.A[0][2]*C[2][0];
    
    n.nx = (C[0][0]*LSQ.b[0] + C[0][1]*LSQ.b[1] + C[0][2]*LSQ.b[2]) / det;
    n.ny = (C[1][0]*LSQ.b[0] + C[1][1]*LSQ.b[1] + C[1][2]*LSQ.b[2]) / det;
    n.nz = (C[2][0]*LSQ.b[0] + C[2][1]*LSQ.b[1] + C[2][2]*LSQ.b[2]) / det;
    
    return n;
}

/**************************************************************************************************/

plane_normal select3 (pixel *pixels, int npts)
{
    int samples[3];
    
    samples[0] = rand() % npts;
    
    do {
        samples[1] = rand() % npts;
    } while (samples[1] == samples[0]);
    
    do {
        samples[2] = rand() % npts;
    } while ((samples[2] == samples[0]) || (samples[2] == samples[1]));
    
    clearLSQ();
    
    updateLSQ(pixels + samples[0]);
    updateLSQ(pixels + samples[1]);
    updateLSQ(pixels + samples[2]);
    
    return solveLSQ();
}

double compute_residual (pixel *p, plane_normal n)
{
    return ( fabs( ((p->u)*n.nx + (p->v)*n.ny + n.nz) - (p->w) ) );
}

int count_inliers (pixel *pixels, int npts, plane_normal n)
{
    int i, count;
    pixel *ptr;
    
    for (count=0, i=0, ptr=pixels; i < npts; ++i, ++ptr) {
        if (compute_residual(ptr, n) < inlier_threshold) ++count;
    }
    
    return count;
}

/* Recursively fit planes to image segments */
void find_planes (pixel *pixels, int npts)
{
    int i, score, best_score, inliers;
    plane_normal normal, best_normal;
    pixel *ptr1, *ptr2, temp;
    
    if (npts < min_points) return;
    
    best_score = 0;
    for (i=0; i < ransac_trials; ++i) {
        
        /* Select 3 points and fit a plane */
        normal = select3(pixels, npts);
        
        score = count_inliers (pixels, npts, normal);
        
        if (score > best_score) {
            best_score = score;
            best_normal = normal;
        }        
    }
    
    /* Find inliers and refit */
    clearLSQ();
    for (i=0, ptr1=pixels; i < npts; ++i, ++ptr1) {
        if (compute_residual(ptr1, best_normal) < inlier_threshold)
            updateLSQ(ptr1);
    }
    normal = solveLSQ();
    
    /* Find final inliers and compact them together at the top of the array */
    for (i=0, ptr1=pixels, ptr2=pixels, inliers=0; i < npts; ++i, ++ptr1) {
        if (compute_residual(ptr1, normal) < inlier_threshold) {
         
            if (ptr1 != ptr2) {
                temp = *ptr2; *ptr2 = *ptr1; *ptr1 = temp;
            }
            
            ++ptr2;
            
            ++inliers;
        }
    }
    
    if (inliers >= min_points) {
        
        /* Add a plane */
        Planes.thePlanes[Planes.nplanes].n              = normal;
        Planes.thePlanes[Planes.nplanes].count          = inliers;
        Planes.thePlanes[Planes.nplanes].start_pixel    = pixels;
        Planes.thePlanes[Planes.nplanes].merged_label   = 0;

        ++(Planes.nplanes);
        
        /* Decide whether to fit the outliers recursively */
        if ((npts-inliers) > outlier_ratio*npts) {
            find_planes(pixels+inliers, (npts-inliers));
        }
    }
    
}

double dotproduct (plane_normal n1, plane_normal n2)
{
    double nx1 = n1.nx, ny1 = n1.ny, nz1 = n1.nz;
    double nx2 = n2.nx, ny2 = n2.ny, nz2 = n2.nz;
    
    return fabs ( (nx1*nx2 + ny1*ny2 + nz1*nz2) / sqrt ((nx1*nx1 + ny1*ny1 + nz1*nz1) * (nx2*nx2 + ny2*ny2 + nz2*nz2)) );
}

void merge_planes ()
{
    int i, j, pass, merged_label, inlier_planes;
    plane *thePlane;
    plane_normal n;
    
    merged_label = 1;
    
    for (i=0; i < Planes.nplanes; ++i) {
        /* Find the next unmerged plane */
        if (Planes.thePlanes[i].merged_label == 0) {
            
            n = Planes.thePlanes[i].n;
            
            Planes.thePlanes[i].merged_label = merged_label;
            
            for (pass=0; pass < 3; ++pass) {
                
                clearLSQ();
                
                updateLSQplane (Planes.thePlanes + i);
                
                for (j=0, thePlane=Planes.thePlanes, inlier_planes=0; j < Planes.nplanes; ++j, ++thePlane) {
                    if (Planes.thePlanes[j].merged_label <= 0) {
                        if ( (dotproduct (n, thePlane->n) > dotpthreshold) &&
                             (count_inliers (thePlane->start_pixel, thePlane->count, n) > inlier_ratio*thePlane->count) ) {
                            /* Note that we use a label of -1 to annotate inliers */
                            thePlane->merged_label = -1;
                            updateLSQplane (Planes.thePlanes + j);
                            ++inlier_planes;
                        } else
                            thePlane->merged_label = 0;
                    }
                }
                
                /* Refit the plane */
                n = solveLSQ();
                
                if (inlier_planes == 0) break;
            }
            
            /* Mark all the current inlier planes - those with merged_labels < 0 */
            for (j=0, thePlane=Planes.thePlanes; j < Planes.nplanes; ++j, ++thePlane)
                if (thePlane->merged_label < 0)
                    thePlane->merged_label = merged_label;
            
            ++merged_label;
            
            if (merged_label > max_merged_planes) break;
        }
    }
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

/**************************************************************************************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    const mxArray* labels;
    mwSize dims[2];
    double *W_ptr;
    /*plane_normal *normals_ptr; */
    int i, j, row, col, npixels;
    unsigned int *labels_ptr, *plane_labels_ptr, *merged_plane_labels_ptr;
    pixel *thePixels, *pixel_ptr;
    

    /* function signature:
    // [merged_plane_labels, plane_labels] = mxFindPlanes (W, labels, focal_length, min_points, ransac_trials, inlier_threshold,
    //                                                     outlier_ratio, inlier_ratio, dotpthreshold, max_merged_planes)
    //  Inputs:
    //   W - nrows x ncols array of doubles 1/depth - NaNs are used to mark missing measurements
    //   labels - nrows x ncols array of uint32 the initial segmentation labels
    //   focal_length - focal_length of the sensor
    //   min_points - minimum number of points in a planar segment
    //   ransac_trials - number of ransac trials to use
    //   inlier_threshold - threshold used to decide if a point is an inlier to a hypothesis
    //   outlier_ratio - fraction of leftover points that triggers a recursive refit
    //   inlier_ratio - ratio used in plane merging to decide whether to merge planes
    //   dotpthreshold - dot product threshold - used in merging to decide if planes are sufficiently similair
    //   max_merged_planes - maximum number of merged planes
    //
    // Outputs:
    //  merged_plane_labels - nrows x ncols array of uint32 - labels of merged planes
    //  plane_labels - nrows x ncols array of uint32 labels of extracted planar segments
    //
    // Note that from the labels you can easily extract the uvw coordinates and the plane normals */
    
    /* check for proper number of arguments */
    if(nrhs!=10) {
        mexErrMsgIdAndTxt("mxFindPlanes:nrhs", "Ten inputs required: W, labels, focal_length, min_points, ransac_trials, inlier_threshold, outlier_ratio, inlier_ratio, dotpthreshold, max_merged_planes");
    }
    
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("mxFindPlanes:nlhs", "Two outputs required: [merged_plane_labels, plane_labels]");
    }
    
    W_ptr = checkArgDoubleArray(prhs, 0, "First argument W must be a 2D array", 0, 0);
    
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    
    /* We assume that the center of the image is the focal point */
    cu = 0.5 * ncols;
    cv = 0.5 * nrows;
    
    /* Check the second argument */
    labels = prhs[1];
    
    if (!mxIsClass(labels, "uint32")) {
        mexErrMsgIdAndTxt("mxFindPlanes:labels", "labels argument needs to be of class uint32");
    }
    
    if (mxGetNumberOfDimensions(labels) != 2) {
        mexErrMsgIdAndTxt("mxFindPlanes:labels", "labels should be an nrows x ncols array");
    }
    
    if ( (mxGetM(labels) != nrows) || (mxGetN(labels) != ncols) ) {
        mexErrMsgIdAndTxt("mxFindPlanes:labels", "labels should be the same size as the W array");
    }
    
    labels_ptr = mxGetData(labels);
    
    
    focal_length        = *checkArgDoubleArray (prhs, 2, "Third argument, focal_length, must be a scalar", 1, 1);
    
    min_points          = (int)(*checkArgDoubleArray (prhs, 3, "Fourth argument, min_points, must be a scalar", 1, 1));
    
    ransac_trials       = (int)(*checkArgDoubleArray (prhs, 4, "Fifth argument, ransac_trials, must be a scalar", 1, 1));
    
    inlier_threshold    = *checkArgDoubleArray (prhs, 5, "Sixth argument, inlier_threshold, must be a scalar", 1, 1);
    
    outlier_ratio       = *checkArgDoubleArray (prhs, 6, "Seventh argument, outlier_ratio, must be a scalar", 1, 1);
    
    inlier_ratio        = *checkArgDoubleArray (prhs, 7, "Eighth argument, inlier_ratio, must be a scalar", 1, 1);
    
    dotpthreshold       = *checkArgDoubleArray (prhs, 8, "Ninth argument, dotpthreshold, must be a scalar", 1, 1);
    
    max_merged_planes   = (int)(*checkArgDoubleArray (prhs, 9, "Tenth argument, max_merged_planes, must be a scalar", 1, 1));
    

    /* create the output matrices */
    /* Note that the documentation says that the entries in these output matrices are initialized to zero. */
    
    
    /* create the output matrices */
    dims[0] = nrows;
    dims[1] = ncols;

    
    plhs[0] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
  /*  plhs[2] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[4] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    normal_x                = mxGetData(plhs[2]);
    normal_y                = mxGetData(plhs[3]);
    normal_z                = mxGetData(plhs[4]); */
    
    /* get a pointer to the real data in the output matrix */
    merged_plane_labels_ptr = mxGetData(plhs[0]);
    plane_labels_ptr        = mxGetData(plhs[1]);


    /* Compute the maximum number of planar segments */
    max_segments = (nrows*ncols) / min_points;
    
    Planes.nplanes = 0;
    Planes.thePlanes = mxMalloc (max_segments * sizeof(plane));

            
    /* First pass create a list of pixels   */  
    thePixels = mxMalloc ( nrows * ncols * sizeof(pixel));
    
    i = 0;
    j = 0;
    
    for (col=0; col < ncols; ++col) {
        for (row=0; row < nrows; ++row) {
            
            if (!mxIsNaN(W_ptr[i])) {
                thePixels[j].label = labels_ptr[i];
                thePixels[j].row = row;
                thePixels[j].col = col;
                thePixels[j].u = (col - cu)/focal_length;
                thePixels[j].v = (row - cv)/focal_length;
                thePixels[j].w = W_ptr[i];
                
                ++j;
            }
            
            ++i;
        }
    }
    
    npixels = j;
    
    /* Sort the pixels by label */
    qsort (thePixels, npixels, sizeof(pixel), compare_pixel_labels);
    
    /* Second pass - fit planes to segments recursively */
 
    i = 0;
    while (i < npixels) {
        /* iterate down the list until you hit the end or encounter a pixel with a different label */
        for (j=i+1; (j < npixels) && (thePixels[j].label == thePixels[i].label); ++j);

        find_planes(thePixels+i, (j-i));

        i = j;
    }
    
    /* Sort the planes in order of decreasing population */
    qsort (Planes.thePlanes, Planes.nplanes, sizeof(plane), compare_plane_counts);
    
    /* Third pass - merge the ssgments */
    merge_planes();
    
    /* Create the output */
    for (i=0; i < Planes.nplanes; ++i) {
        for (j=0, pixel_ptr = Planes.thePlanes[i].start_pixel; j < Planes.thePlanes[i].count; ++j, ++pixel_ptr)
            plane_labels_ptr[(pixel_ptr->col)*nrows + (pixel_ptr->row)] = (i+1);
        
        /* Note the (i+1) to account for zero offset - pixels that are not part of a plane should have a 0 label */
    }
    
    for (i=0; i < Planes.nplanes; ++i) {
        if (Planes.thePlanes[i].merged_label) {
            for (j=0, pixel_ptr = Planes.thePlanes[i].start_pixel; j < Planes.thePlanes[i].count; ++j, ++pixel_ptr)
                merged_plane_labels_ptr[(pixel_ptr->col)*nrows + (pixel_ptr->row)] = Planes.thePlanes[i].merged_label;
                /*Create a new array to store the normals, re-compile the code*/
               /* plane_normal n = Planes.thePlanes[i].n; 
                normal_x[(pixel_ptr->col)*nrows + (pixel_ptr->row)] = Planes.thePlanes[i].n.nx;
                normal_y[(pixel_ptr->col)*nrows + (pixel_ptr->row)] = Planes.thePlanes[i].n.ny;
                normal_z[(pixel_ptr->col)*nrows + (pixel_ptr->row)] = Planes.thePlanes[i].n.nz; */
                /*
                double normal_x = n.nx; 
                normals_ptr[(pixel_ptr->col)*nrows + (pixel_ptr->row)] = Planes.thePlanes[i].n; */
                

        }
    }
    
   /* 
    pl_dims[0] = counter + 1;
    pl_dims[1] = 1;
    plhs[2] = mxCreateNumericArray(2, pl_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericArray(2, pl_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[4] = mxCreateNumericArray(2, pl_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[5] = mxCreateNumericArray(2, pl_dims, mxDOUBLE_CLASS, mxREAL);
    normal_x                = mxGetData(plhs[2]);
    normal_y                = mxGetData(plhs[3]);
    normal_z                = mxGetData(plhs[4]);
    lab                     = mxGetData(plhs[5]);
    
    for (i=0; i < Planes.nplanes; ++i) {
        if (Planes.thePlanes[i].merged_label){
            normal_x[i+1] = Planes.thePlanes[i].n.nx;
            normal_y[i+1] = Planes.thePlanes[i].n.ny;
            normal_z[i+1] = Planes.thePlanes[i].n.nz;
           lab[i+1] = Planes.thePlanes[i].merged_label;
            
        }
    }
    */
    

    mxFree(thePixels);
    mxFree(Planes.thePlanes);
}
