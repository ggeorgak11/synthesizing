/*
// mxMergeTriangles.c
// Mex routine that starts with a weighted graph derived from the Delaunay Triangulation and repeatedly
// merges regions based on a computed merge_cost inspired by normalized cuts
*/

#include "mex.h"
#include <math.h>

#define MSG_ID "mxMergeTriangles:arg"

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

/******************************************************************************************/

#define MIN(a, b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a, b) ( ((a) > (b)) ? (a) : (b) )

/******************************************************************************************/

/*
// Type declarations
// */

typedef struct edge_tag {
    /* The indices of the two nodes joined by the edge */
    int idx1, idx2;
    
    double merge_cost;    /* Average cost of this edge, cost per unit length */
    double total_length;  /* Total length of this edge in pixels in the image */
    
    /* This edge will appear in two adjacency lists - one for node1 and another for node2 */
    struct edge_tag *next1, *prev1;
    struct edge_tag *next2, *prev2;
    
    /* We need to store the current index of the edge in the heap to make updating and deleting more efficient. */
    int heap_index;
} edge;

typedef struct {
    /* Could maintain a count of how many edges emanate from this node. It could be used to slightly optimize
    /* the node merging operation */
    
    /* Pointer to a list of all of the edges impinging on this node */
    struct edge_tag *edge_list;
} node;

/******************************************************************************************/

/*
// Data structures that we use to allocate nodes, edges etc.
// */

struct {
    int count;
    node *buf;
} node_buffer;

struct {
    int count;
    edge *buf;
} edge_buffer;

struct {
    int count;  /* total number of edges currently in the heap */
    edge **elts;
} edge_heap;

/* Used as a table during merging to indicate which nodes we have edges to already */
edge **adjacency_table;

/******************************************************************************************/

/* Globals */

/* ntri - Number of nodes in the graph of Delaunay triangles */
mwSize ntri;

/* nedges - Maximum number of edges - the real number of edges may be less than this */
mwSize nedges;

/******************************************************************************************/

/* Some basic operations on edges and nodes */

int other_idx(edge *theEdge, int idx) {
    return ( (theEdge->idx1 == idx) ? theEdge->idx2 : theEdge->idx1 );
}

void update_index (edge *theEdge, int old_idx, int new_idx) {
    if (theEdge->idx1 == old_idx)
        theEdge->idx1 = new_idx;
    else
        theEdge->idx2 = new_idx;
}

struct edge_tag *get_next_edge (edge *theEdge, int idx) {
    return (theEdge->idx1 == idx) ? theEdge->next1 : theEdge->next2;
}

struct edge_tag *get_prev_edge (edge *theEdge, int idx) {
    return (theEdge->idx1 == idx) ? theEdge->prev1 : theEdge->prev2;
}


void set_next_edge (edge *theEdge, int idx, struct edge_tag *ptr) {
    if (theEdge->idx1 == idx)
        theEdge->next1 = ptr;
    else
        theEdge->next2 = ptr;
}

void set_prev_edge (edge *theEdge, int idx, struct edge_tag *ptr) {
    if (theEdge->idx1 == idx)
        theEdge->prev1 = ptr;
    else
        theEdge->prev2 = ptr;
}

/******************************************************************************************/

/* Routines to manipulate adjacency lists */

void add_to_adjacency_list(int node_index, edge *theEdge) {
    node *theNode = node_buffer.buf + node_index;
    edge *head_elt = theNode->edge_list;
                
    set_next_edge (theEdge, node_index, head_elt);
    set_prev_edge (theEdge, node_index, NULL);
    
    if (head_elt) {
        set_prev_edge (head_elt, node_index, theEdge);
    }
    
    theNode->edge_list = theEdge;
}

void remove_from_adjacency_list(int node_index, edge *theEdge) {
    node *theNode = node_buffer.buf + node_index;
    edge *prev = get_prev_edge (theEdge, node_index);
    edge *next = get_next_edge (theEdge, node_index);
    
    if (prev) {
        set_next_edge(prev, node_index, next);
    } else {
        /* If the prev field is NULL we are looking at the start of the list */
        theNode->edge_list = next;
    }
    
    if (next) {
        set_prev_edge (next, node_index, prev);
    }
}

/******************************************************************************************/

void heap_swap(int idx1, int idx2) {
    edge* temp = edge_heap.elts[idx1];
    edge_heap.elts[idx1] = edge_heap.elts[idx2];
    edge_heap.elts[idx2] = temp;
    
    edge_heap.elts[idx1]->heap_index = idx1;
    edge_heap.elts[idx2]->heap_index = idx2;
}


/* These routines swap elements as necessary to restore the min_heap property */
void heapify_up(int k) {
    double merge_cost = edge_heap.elts[k]->merge_cost;
    int parent;
    
    while (k > 0) {
        parent = (k-1) / 2;
        if (edge_heap.elts[parent]->merge_cost > merge_cost) {
            heap_swap(parent, k);
            k = parent;
        } else
            break;
    }
}

void heapify_down(int k) {
    double merge_cost = edge_heap.elts[k]->merge_cost;
    int child;
    
    while (1) {
        
        /* Find the child with the minimum merge cost if any. Note that it is important that you
        // swap with the minimum cost child to ensure that you preserve the min_heap property or
        // something with a larger value could end up the parent of something with a smaller value.
        // If there are no children then break */
        child = 2*k + 1;
        
        if (child >= edge_heap.count)
            break;
        
        if ( ((child+1) < edge_heap.count) && (edge_heap.elts[child+1]->merge_cost < edge_heap.elts[child]->merge_cost) )
            child = child + 1;
        
        if (edge_heap.elts[child]->merge_cost < merge_cost) {
            heap_swap(child, k);
            k = child;
        } else
            break; /* If you don't swap you are all done. */
    }
}

void add_to_edge_heap(edge *theEdge) {
    edge_heap.elts[edge_heap.count] = theEdge;
    theEdge->heap_index = edge_heap.count;
    heapify_up(edge_heap.count);
    ++edge_heap.count;
}

void remove_from_edge_heap(edge *theEdge) { 
    /* resist the temptation to eliminate this idx since theEdge->heap_index is changed
    // by  the subsequent heap_swap call. */
    int idx = theEdge->heap_index;
    double new_merge_cost, old_merge_cost;
    
    --(edge_heap.count);
    
    old_merge_cost = edge_heap.elts[idx]->merge_cost;
    new_merge_cost = edge_heap.elts[edge_heap.count]->merge_cost;
    
    heap_swap(idx, edge_heap.count); /* Swap with the last element in the heap */
    
    if (new_merge_cost < old_merge_cost)
        heapify_up(idx);    /* percolate up */
    else
        heapify_down(idx);  /* percolate down */
}


/******************************************************************************************/


void InitNodes(int *SN_ptr) {
    int i;
    node *node_ptr;
    
    for (i=0, node_ptr = node_buffer.buf; i < ntri; ++i, ++node_ptr) {
        node_ptr->edge_list = NULL;
    }
}


void InitEdges(int *SN_ptr, double *edge_weights_ptr, double *edge_lengths_ptr) {
    int i, j, idx;
    edge *theEdge;
    
    for (i=0; i < ntri; ++i) {
        for (j=0; j < 3; ++j) {
            idx = SN_ptr[i + j*ntri];
            
            /* Check for case where idx <= 0 indicating a missing edge */
            if (idx > 0) {
                
                if (idx > ntri) {
                    mexErrMsgIdAndTxt("mxMergeTriangles:SN", "index entry in SN too large");
                }
                
                idx = idx - 1; /* account for 1 based indexing */
                
                /* Note that we only add an edge the first time it is encounter since we check whether
                // (i < idx). We assume that the next version of this edge in the table has the same info but we don't check it. */
                if (i < idx) {
                    /* Create the edge */
                    theEdge = edge_buffer.buf + edge_buffer.count;
                    ++edge_buffer.count;
                    
                    theEdge->idx1 = i;
                    theEdge->idx2 = idx;
                    
                    theEdge->merge_cost   = edge_weights_ptr[i + j*ntri];
                    theEdge->total_length = edge_lengths_ptr[i + j*ntri];
                    
                    if (theEdge->merge_cost < 0) {
                        mexErrMsgIdAndTxt("mxMergeTriangles:edge_weights", "edge_weight entry negative");
                    }
                    
                    if (theEdge->total_length < 0) {
                        mexErrMsgIdAndTxt("mxMergeTriangles:edge_lengths", "edge_length entry is < 0");
                    }

                    theEdge->prev1 = NULL;
                    theEdge->next1 = NULL;
                    
                    theEdge->prev2 = NULL;
                    theEdge->next2 = NULL;
                    
                    /* Add it to i's edge_list */
                    add_to_adjacency_list(i, theEdge);
                    
                    /* add it to idx's edge_list */
                    add_to_adjacency_list(idx, theEdge);
                    
                    /* add it to the edge_heap */
                    add_to_edge_heap(theEdge);
                }
                
            }
        }
    }
}


/******************************************************************************************/

void merge_edges (edge *edge1, edge *edge2) {
    double merge_cost1, merge_cost2;
    double length1, length2;
    double new_merge_cost, new_length;
    
    /* merge edge2 into edge1 */
    
    /* Update merge costs 
    // Update lengths
    // Update heap up or down depending on change */
    
    merge_cost1 = edge1->merge_cost;
    merge_cost2 = edge2->merge_cost;
    
    length1 = edge1->total_length;
    length2 = edge2->total_length;
    
    /* New total edge length */
    new_length = length1 + length2;
    
    /* Convex combination of the two edge costs so the results should be somewhere in between */
    new_merge_cost = (length1*merge_cost1 + length2*merge_cost2) / new_length;
    
    /* This is put in just in case numerical error causes us to have a merge cost that is smaller than either component
    // which is mathematically impossible.
    // new_merge_cost = MAX(new_merge_cost, MIN(merge_cost1, merge_cost2)); */
        
    edge1->merge_cost = new_merge_cost;
    edge1->total_length = new_length;
    
    /* Update the heap */
    if (new_merge_cost < merge_cost1)
        heapify_up (edge1->heap_index);   /* merge_cost decreased */
    else
        heapify_down (edge1->heap_index); /* merge_cost increased */
}

/* Code for merging nodes - tricky stuff */

void merge_nodes(int idx1, int idx2) {
    int head_idx;
    node *node1, *node2;
    edge *list1, *list2, *head_elt, *clear_list;
    
    if (idx1 == idx2) return;
    
    node1 = node_buffer.buf + idx1;
    node2 = node_buffer.buf + idx2;
    
    list1 = node1->edge_list;
    list2 = node2->edge_list;
    
    /* Clear node1's adjacency list */
    node1->edge_list = NULL;
    
    /* Merge the lists of edges */

    /* We tackle list2 first    */         
    while (list2) {
        head_elt = list2;
        list2 = get_next_edge(list2, idx2); /* Pop list2 */
        
        head_idx = other_idx(head_elt, idx2);
        
        /* Update indices to replace references to idx2 with idx1 */
        update_index (head_elt, idx2, idx1);
        
        /* store in adjacency_table for later reference */
        adjacency_table[head_idx] = head_elt;
        
        /* Add to new adjacency list */
        add_to_adjacency_list (idx1, head_elt);
        
        /* Note that since we are not changing the merge cost yet we do not update the heap. */
    }
    
    /* Retain start of list of things that need to be cleared from the adjacency_table */
    clear_list = node1->edge_list;
    
    while (list1) {
        head_elt = list1;
        list1 = get_next_edge(list1, idx1); /* Pop list1 */
        
        head_idx = other_idx(head_elt, idx1);
        
        if (adjacency_table[head_idx]) {
            
            /* We have a node that is linked to both idx1 and idx2 so we merge the edges */
            merge_edges (adjacency_table[head_idx], head_elt);
            
            /* Remove the edge from the other nodes adjacency list */
            remove_from_adjacency_list(head_idx, head_elt);
            
            /* Remove the edge from the heap */
            remove_from_edge_heap(head_elt);
            
        } else {
            /* Add to new adjacency list */
            add_to_adjacency_list (idx1, head_elt);
            
            /* Note that since we are not changing the merge cost we do not update the heap. */
        }
        
    }
    
    /* Clear the adjacency_table - note that we only clear things that were used */
    while (clear_list) {
        adjacency_table[other_idx(clear_list, idx1)] = NULL;
        clear_list = get_next_edge(clear_list, idx1);
    }
   
}

/******************************************************************************************/


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {
    int *SN_ptr;
    double *edge_weights_ptr, *edge_lengths_ptr, *merge_cost_ptr, *merge_length_ptr;
    int *old_labels_ptr, *new_labels_ptr;
    mwSize dims[2];
    const mxArray *SN;
    edge *best_edge;
    int i;
    
    /* function signature: 
    /* [old_labels, new_labels, merge_costs, merge_lengths] = mxMergeTriangles (SN, edge_weights, edge_lengths) */
    
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("mxMergeTriangles:nrhs", "Three inputs required: SN, edge_weights and edge_lengths");
    }
    
    if(nlhs!=4) {
        mexErrMsgIdAndTxt("mxMergeTriangles:nlhs", "Four outputs required.");
    }
    
    /* Check the second argument */
    edge_weights_ptr = checkArgDoubleArray(prhs, 1, "Second argument must be an ntri x 3 array", 0, 3);
    
    ntri = mxGetM(prhs[1]);

    /* Check the third argument */
    edge_lengths_ptr = checkArgDoubleArray(prhs, 2, "Third argument must be an ntri x 3 array", ntri, 3);
    
    /* Check the first argument - note that NaN entries will be converted to zeros and that the entries
    // in the SN array are indexed from 1 rather than 0. */
    SN = prhs[0];
    if (!mxIsClass(SN, "int32")) {
        mexErrMsgIdAndTxt("mxMergeTriangles:SN", "SN argument needs to be of class int32");
    }
    
    if (mxGetNumberOfDimensions(SN) != 2 || mxGetM(SN) != ntri || mxGetN(SN) != 3) {
        mexErrMsgIdAndTxt("mxMergeTriangles:SN", "SN should be an ntri x 3 array");
    }
    
    SN_ptr = mxGetData(SN);

    /* nedges this is (ntri * 3) / 2 since each edge will only be recorded once we are assuming a symmetric graph where the
    // w_ij = w_ji for all i, j. */
    nedges = ((ntri * 3) / 2) + 1;
    
    /* create the output matrices */
    dims[0] = ntri-1;
    dims[1] = 1;
    
    plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    
    plhs[2] = mxCreateDoubleMatrix(ntri-1, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(ntri-1, 1, mxREAL);

    
    old_labels_ptr = mxGetData(plhs[0]);
    new_labels_ptr = mxGetData(plhs[1]);
    
    merge_cost_ptr   = mxGetData(plhs[2]);
    merge_length_ptr = mxGetData(plhs[3]);

    /*
    // Preallocate storage for the graph
    // */
    
    node_buffer.count = 0;
    node_buffer.buf = mxMalloc( ntri * sizeof(node) );
    
    edge_buffer.count = 0;
    /* Note that we only nee (ntri * 3)/2 edges since each edge is recorded once not twice since we are assuming an undirected graph */
    edge_buffer.buf = mxMalloc( nedges * sizeof(edge) );
    
    edge_heap.count = 0;
    /* In the heap we store pointers to the edges rather than the edges themselves */
    edge_heap.elts = mxMalloc( nedges * sizeof(edge*) );
    
    /* In the adjacency_table we store pointers to the edges rather than the edges themselves */
    adjacency_table = mxCalloc( ntri, sizeof(edge*) );
    
    /*
    // Initialize all of the data structures
    // */
    
    /* Initialize the nodes */
    InitNodes(SN_ptr);
    
    /* Initialize the edges and add them to the heap */
    InitEdges(SN_ptr, edge_weights_ptr, edge_lengths_ptr);
    
    /*
    // Main computational loop
    // */
    
    for (i=0; i < (ntri-1); ++i) {
        
        /* The edge at the beginning of the min_heap is the edge with smallest merge_cost */
        best_edge = edge_heap.elts[0];
        
        /* record the merge - translate the results to 1 based indexing */
        new_labels_ptr[i]   = best_edge->idx1 + 1;
        old_labels_ptr[i]   = best_edge->idx2 + 1;
        merge_cost_ptr[i]   = best_edge->merge_cost;
        merge_length_ptr[i] = best_edge->total_length;
        
        /* Debug
        // mexPrintf("iteration %d : merging %d and %d with %8.4f cost\n", i, old_labels_ptr[i],
        //        new_labels_ptr[i], merge_cost_ptr[i]); */
        
        /* Remove the edge from both adjacency lists and the heap */
        remove_from_adjacency_list (best_edge->idx1, best_edge);
        remove_from_adjacency_list (best_edge->idx2, best_edge);
        
        remove_from_edge_heap(best_edge);
        
        merge_nodes(best_edge->idx1, best_edge->idx2);
    }

    
    /*
    // Free storage for the graph
    // */
    
    mxFree(node_buffer.buf);
    mxFree(edge_buffer.buf);
    mxFree(edge_heap.elts);
    mxFree(adjacency_table);
}
