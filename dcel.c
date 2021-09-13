#include "dcel.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

/* Need to know about watchTowers to store them in the face. */
#include "watchtowerStruct.h"
/* Need to talk to watchTowers to store face in them. */
#include "watchtowerStruct.c"

#define INITIALVERTICES 4
#define INITIALEDGES 4
#define INITIALFACES 1
#define NOVERTEX (-1)
#define NOEDGE (-1)

#define DIR_UNDECIDED (0)
#define INSIDE (1)
#define OUTSIDE (-1)
#define NODIAMETER (-1)

#define MAX_BUFFER_SIZE 256

struct halfEdge;
struct vertex;

struct vertex {
    double x;
    double y;
};

struct face {
    // Any half-edge in the face.
    struct halfEdge *he;
    // The watchtower stored for the face.
    struct watchtowerStruct *wt;
};

struct split {
    int startEdge;
    int endEdge;
    struct vertex startSplitPoint;
    struct vertex endSplitPoint;
    int verticesSpecified;
    int enforceFace;
    int expectedFace;
};

struct halfEdge {
    int startVertex;
    int endVertex;
    struct halfEdge *next;
    struct halfEdge *prev;
    struct halfEdge *pair;
    int face;
    int edge;
};

struct edge {
    // One of the half-edges corresponding to this edge, used for splits.
    struct halfEdge *halfEdge;
};

struct DCEL {
    struct edge *edges;
    int edgesUsed;
    int edgesAllocated;
    struct face *faces;
    int facesUsed;
    int facesAllocated;
    struct vertex *vertices;
    int verticesUsed;
    int verticesAllocated;
};


struct bisector {
    double grad;
    double mp_x;
    double mp_y;
    int finiteGrad;
};

/* Gets a point at least distance away from the midpoint of the bisector given. */
void getBisectorPoint(double distance, struct bisector *b, double *x, double *y);

void getBisectorPoint(double distance, struct bisector *b, double *x, double *y){

    // for infinite gradient case, just subtract distance from mp_y to get point
    if (! b->finiteGrad) {
        *x = b->mp_x;
        *y = b->mp_y - distance;
    }
    /* for finite gradient case, use equations and solve simultaneously:
     *
     *      alpha^2 + beta^2 = distance^2    (i.e. Pythagoras theorem, b was taken)
     *      beta / alpha = grad
     *
     *      // solved for beta in terms of alpha
     *      beta = grad * alpha
     *
     *      // solve for alpha (we know grad and distance)
     *      alpha^2 + (grad * alpha)^2 = distance^2
     *      alpha^2 + grad^2 * alpha^2 = distance^2
     *      alpha^2 * (1 + grad^2) = distance^2
     *      alpha^2 = distance^2 / (1 + grad^2)
     *      alpha = sqrt( distance^2 / (1 + grad^2) )
     *
     *      // solve for beta
     *      beta = grad * alpha
     *
     */
    else {

        // sign variable ensures the sign of distance determines if
        // alpha and beta are subtracted or added to get the new point
        int sign = 1;
        if (distance < 0) {
            sign = -1;
        }
        double alpha = sign * ( sqrt(pow(distance, 2) / (1 + pow(b->grad, 2))) );
        double beta = alpha * b->grad;
        *x = b->mp_x - alpha;
        *y = b->mp_y - beta;
    }
}

struct bisector *getBisector(double x1, double y1, double x2, double y2);

struct bisector *getBisector(double x1, double y1, double x2, double y2){

    struct bisector *p = (struct bisector*)malloc(sizeof(struct bisector));
    assert(p);

    // store midpoint of the pair of points in bisector struct
    p->mp_x = (x1 + x2) / 2.0;
    p->mp_y = (y1 + y2) / 2.0;

    // store whether the gradient is finite or infinite, store gradient if finite
    if (y2 - y1 == 0) {
        p->finiteGrad = 0;
    }
    else {
        p->finiteGrad = 1;
        p->grad = -((x2 - x1) / (y2 - y1));
    }

    return p;
}

struct bisector *readNextBisector(FILE *bisectorfile){

    double x1, y1, x2, y2;
    
    if(fscanf(bisectorfile, "%lf %lf %lf %lf", &x1, &y1, &x2, &y2) != 4){
        return NULL;
    }
    
    return getBisector(x1, y1, x2, y2);
}

char *getBisectorEquation(struct bisector *b){

    if(! b){
        return NULL;
    }

    char *returnString = (char*)malloc(MAX_BUFFER_SIZE);
    assert(returnString);

    if(! b->finiteGrad){

        int stringLength = snprintf(returnString, MAX_BUFFER_SIZE, "x = %lf",
                                    b->mp_x);
        free(returnString);
        returnString = (char *) malloc(sizeof(char) * (stringLength + 1));
        assert(returnString);
        sprintf(returnString, "x = %lf", b->mp_x);

    } else {

        int stringLength = snprintf(returnString, MAX_BUFFER_SIZE,
                                    "y = %lf * (x - %lf) + %lf", b->grad, b->mp_x,
                                    b->mp_y);
        free(returnString);
        returnString = (char *) malloc(sizeof(char) * (stringLength + 1));
        assert(returnString);
        sprintf(returnString,
                "y = %lf * (x - %lf) + %lf", b->grad, b->mp_x,
                b->mp_y);
    }

    return returnString;
}

void freeBisector(struct bisector *bisector){
    if(bisector){
        free(bisector);
    }
}

enum intersectType;

enum intersectType {
    DOESNT_INTERSECT = 0,  // Doesn't intersect
    INTERSECT = 1,         // Intersects
    SAME_LINE_OVERLAP = 2, // Lines are the same
    ENDS_OVERLAP = 3       // Intersects at exactly one point (endpoint)
};

struct intersection {
    int edge_i;
    struct vertex vertex_i;
    int edge_f;
    struct vertex vertex_f;
    int multipleIntersects;
};

/* 
This intersection is based on code by Joseph O'Rourke and is provided for use in 
COMP20003 Assignment 2.

The approach for intersections is:
- Use the bisector to construct a finite segment and test it against the half-edge.
- Use O'Rourke's segseg intersection (https://hydra.smith.edu/~jorourke/books/ftp.html)
    to check if the values overlap./intersectio
*/
/* 
    Generates a segment with each end at least minLength away in each direction 
    from the bisector midpoint. Returns 1 if b intersects the given half-edge
    on this segment, 0 otherwise. Sets the intersection point to the given x, y
    positions.
*/

/* Returns -1, 0 or 1, based on the area enclosed by the three points. 0 corresponds
    to no area enclosed.
*/
int areaSign(double sx, double sy, double ex, double ey, double x, double y);

/* Returns 1 if the point (x, y) is in the line from s(x, y) to e(x, y), 0 otherwise. */
int collinear(double sx, double sy, double ex, double ey, double x, double y);
    
int collinear(double sx, double sy, double ex, double ey, double x, double y){
    /* Work out area of parallelogram - if it's 0, points are in the same line. */
    if (areaSign(sx, sy, ex, ey, x, y) == 0){
        return 1;
    } else {
        return 0;
    }
}

int areaSign(double sx, double sy, double ex, double ey, double x, double y){
    double areaSq;
    /* |AB x AC|^2, squared area */
    /* See https://mathworld.wolfram.com/CrossProduct.html */
    areaSq = (ex - sx) * (y  - sy) -
            (x  - sx) * (ey - sy);
    
    if(areaSq > 0.0){
        return 1;
    } else if(areaSq == 0.0){
        return 0;
    } else {
        return -1;
    }
}

/* Returns 1 if point (x, y) is between (sx, sy) and (se, se) */
int between(double sx, double sy, double ex, double ey, double x, double y);

int between(double sx, double sy, double ex, double ey, double x, double y){
    if(sx != ex){
        /* If not vertical, check whether between x. */
        if((sx <= x && x <= ex) || (sx >= x && x >= ex)){
            return 1;
        } else {
            return 0;
        }
    } else {
        /* Vertical, so can't check _between_ x-values. Check y-axis. */
        if((sy <= y && y <= ey) || (sy >= y && y >= ey)){
            return 1;
        } else {
            return 0;
        }
    }
}

enum intersectType parallelIntersects(double heSx, double heSy, double heEx, double heEy,
        double bSx, double bSy, double bEx, double bEy, double *x, double *y);

enum intersectType parallelIntersects(double heSx, double heSy, double heEx, double heEy,
        double bSx, double bSy, double bEx, double bEy, double *x, double *y){
    if(!collinear(heSx, heSy, heEx, heEy, bSx, bSy)){
        /* Parallel, no intersection so don't set (x, y) */
        return DOESNT_INTERSECT;
    }
    /* bS between heS and heE */
    if(between(heSx, heSy, heEx, heEy, bSx, bSy)){
        *x = bSx; 
        *y = bSy;
        return SAME_LINE_OVERLAP;
    }
    /* bE between heS and heE */
    if(between(heSx, heSy, heEx, heEy, bEx, bEy)){
        *x = bEx;
        *y = bEy;
        return SAME_LINE_OVERLAP;
    }
    /* heS between bS and bE */
    if(between(bSx, bSy, bEx, bEy, heSx, heSy)){
        *x = heSx;
        *y = heSy;
        return SAME_LINE_OVERLAP;
    }
    /* heE between bS and bE */
    if(between(bSx, bSy, bEx, bEy, heEx, heEy)){
        *x = heEx; 
        *y = heEy;
        return SAME_LINE_OVERLAP;
    }
    
    return DOESNT_INTERSECT;
}

enum intersectType intersects(struct halfEdge *he, struct bisector *b, 
        struct DCEL *dcel, double minLength, double *x, double *y);

enum intersectType intersects(struct halfEdge *he, struct bisector *b, 
        struct DCEL *dcel, double minLength, double *x, double *y){

    /* Half-edge x, y pair */
    double heSx = dcel->vertices[he->startVertex].x;
    double heSy = dcel->vertices[he->startVertex].y;
    double heEx = dcel->vertices[he->endVertex].x;
    double heEy = dcel->vertices[he->endVertex].y;
    
    /* Bisector x, y pair */
    double bSx;
    double bSy;
    double bEx;
    double bEy;

    // get start and end points for line segment inside bisector line equation
    getBisectorPoint(minLength, b, &bSx, &bSy);
    getBisectorPoint(-1 * minLength, b, &bEx, &bEy);

    /* Parametric equation parameters */
    double t1;
    double t2;
    /* Numerators for X and Y coordinate of intersection. */
    double numeratorX;
    double numeratorY;
    /* Denominators of intersection coordinates. */
    double denominator;
    
    /*
    See http://www.cs.jhu.edu/~misha/Spring20/15.pdf
    for explanation and intuition of the algorithm here.
    x_1 = heSx, y_1 = heSy    |    p_1 = heS
    x_2 = heEx, y_2 = heEy    |    q_1 = heE
    x_3 = bSx , y_3 = bSy     |    p_2 =  bS
    x_4 = bEx , y_4 = bEy     |    q_2 =  bE
    ----------------------------------------
    So the parameters t1 and t2 are given by:
    | t1 |   | heEx - heSx  bSx - bEx | -1  | bSx - heSx |
    |    | = |                        |     |            |
    | t2 |   | heEy - heSy  bSy - bEy |     | bSy - heSy |
    
    Hence:
    | t1 |       1     | bSy - bEy        bEx - bSx |  | bSx - heSx |
    |    | = --------- |                            |  |            |
    | t2 |    ad - bc  | heSy - heEy    heEx - heSx |  | bSy - heSy |
    
        where 
        a = heEx - heSx
        b = bSx  -  bEx
        c = heEy - heSy
        d = bSy  -  bEy
    */
    
    /* Here we calculate ad - bc */
    denominator = heSx * (bEy  -  bSy) +
            heEx * (bSy  -  bEy) +
            bEx  * (heEy - heSy) +
            bSx  * (heSy - heEy);
    
    if(denominator == 0){
        /* In this case the two are parallel */
        return parallelIntersects(heSx, heSy, heEx, heEy, bSx, bSy, bEx, bEy, x, y);
    }
    
    /*
    Here we calculate the top row.
    | bSy - bEy        bEx - bSx |  | bSx - heSx |
    |                            |  |            |
    |                            |  | bSy - heSy |
    */
    numeratorX = heSx * (bEy  -  bSy) +
            bSx  * (heSy -  bEy) +
            bEx  * (bSy  - heSy);
    
    /*
    Here we calculate the bottom row.
    |                            |  | bSx - heSx |
    |                            |  |            |
    | heSy - heEy    heEx - heSx |  | bSy - heSy |
    */
    numeratorY = -(heSx * (bSy  -  heEy) +
            heEx * (heSy -  bSy) +
            bSx  * (heEy  - heSy));
    
    /* Use parameters to convert to the intersection point */
    t1 = numeratorX/denominator;
    t2 = numeratorY/denominator;
    *x = heSx + t1 * (heEx - heSx);
    *y = heSy + t1 * (heEy - heSy);
    
    /* Make final decision - if point is on segments, parameter values will be
    between 0, the start of the line segment, and 1, the end of the line segment.
    */
    if (0.0 < t1 && t1 < 1.0 && 0.0 < t2 && t2 < 1.0){
        return INTERSECT;
    } else if(t1 < 0.0 || 1.0 < t1 || t2 < 0.0 || 1.0 < t2){
        /* s or t outside of line segment. */
        return DOESNT_INTERSECT;
    } else {
        /* 
        ((numeratorX == 0) || (numeratorY == 0) || 
         (numeratorX == denominator) || (numeratorY == denominator))
        */
        return ENDS_OVERLAP;
    }
}

char *getIntersectionString(struct intersection *intersection){

    if(! intersection){
        return NULL;
    }

    char *returnString = (char*)malloc(MAX_BUFFER_SIZE);
    assert(returnString);

    if(! intersection->multipleIntersects){

        int stringLength = snprintf(returnString, MAX_BUFFER_SIZE,
                                    "From Edge %d (%lf, %lf) to Edge %d (%lf, %lf)",
                                    intersection->edge_i, intersection->vertex_i.x, intersection->vertex_i.y,
                                    intersection->edge_f, intersection->vertex_f.x, intersection->vertex_f.y);
        free(returnString);
        returnString = (char *) malloc(sizeof(char) * (stringLength + 1));
        assert(returnString);
        sprintf(returnString, "From Edge %d (%lf, %lf) to Edge %d (%lf, %lf)",
                intersection->edge_i, intersection->vertex_i.x, intersection->vertex_i.y,
                intersection->edge_f, intersection->vertex_f.x, intersection->vertex_f.y);
    } else {

        int stringLength = snprintf(returnString, MAX_BUFFER_SIZE,
                                    "From Edge %d (%lf, %lf) to Edge %d (%lf, %lf)",
                                    intersection->edge_i, intersection->vertex_i.x, intersection->vertex_i.y,
                                    intersection->edge_f, intersection->vertex_f.x, intersection->vertex_f.y);
        free(returnString);
        returnString = (char *) malloc(sizeof(char) * (stringLength + 1));
        assert(returnString);
        sprintf(returnString, "From Edge %d (%lf, %lf) to Edge %d (%lf, %lf)",
                intersection->edge_i, intersection->vertex_i.x, intersection->vertex_i.y,
                intersection->edge_f, intersection->vertex_f.x, intersection->vertex_f.y);
    }
    return returnString;
}

void freeIntersection(struct intersection *intersection){
    if(! intersection){
        return;
    }
    free(intersection);
}

struct DCEL *newDCEL();

struct DCEL *newDCEL(){
    /* Setup DCEL. */
    struct DCEL *dcel = (struct DCEL *) malloc(sizeof(struct DCEL));
    assert(dcel);
    dcel->edges = NULL;
    dcel->edgesUsed = 0;
    dcel->edgesAllocated = 0;
    dcel->faces = NULL;
    dcel->facesUsed = 0;
    dcel->facesAllocated = 0;
    dcel->vertices = NULL;
    dcel->verticesUsed = 0;
    dcel->verticesAllocated = 0;
    
    return dcel;
}

/* Allocate a new halfEdge and return it. */
struct halfEdge *newHalfEdge();

struct halfEdge *newHalfEdge(){
    struct halfEdge *he = (struct halfEdge *) malloc(sizeof(struct halfEdge));
    assert(he);
    he->startVertex = NOVERTEX;
    he->endVertex = NOVERTEX;
    he->next = NULL;
    he->prev = NULL;
    he->pair = NULL;
    he->face = NOFACE;
    he->edge = NOEDGE;
    return he;
}

/* Returns INSIDE if the points is on the INSIDE of the vector pair by the CW winding 
    order, OUTSIDE if it is OUTSIDE by the CW winding order and DIR_UNDECIDED if the 
    point lies on the vector between the points v1 and v2. */
int getRelativeDir(double x, double y, struct vertex *v1, struct vertex *v2);

/* Check there's space for another vertex in the DCEL, or increase the allocated space. */
void ensureSpaceForVertex(struct DCEL *dcel);

void ensureSpaceForVertex(struct DCEL *dcel){
    if(! (dcel->vertices)){
        dcel->vertices = (struct vertex *) 
                malloc(sizeof(struct vertex) * INITIALVERTICES);
        assert(dcel->vertices);
        dcel->verticesAllocated = INITIALVERTICES;
    } else if((dcel->verticesUsed + 1) > dcel->verticesAllocated){
        dcel->vertices = (struct vertex *) realloc(dcel->vertices, 
                sizeof(struct vertex) * dcel->verticesAllocated * 2);
        assert(dcel->vertices);
        dcel->verticesAllocated = dcel->verticesAllocated * 2;
    }
}

/* Check there's space for another edge in the DCEL, or increase the allocated space. */
void ensureSpaceForEdge(struct DCEL *dcel);

void ensureSpaceForEdge(struct DCEL *dcel){
    if(! (dcel->edges)){
        dcel->edges = (struct edge *) 
                malloc(sizeof(struct edge) * INITIALEDGES);
        assert(dcel->edges);
        dcel->edgesAllocated = INITIALEDGES;
    } else if((dcel->edgesUsed + 1) > dcel->edgesAllocated){
        dcel->edges = (struct edge *) realloc(dcel->edges, 
                sizeof(struct edge) * dcel->edgesAllocated * 2);
        assert(dcel->edges);
        dcel->edgesAllocated = dcel->edgesAllocated * 2;
    }
}

/* Check there's space for another face in the DCEL, or increase the allocated space. */
void ensureSpaceForFace(struct DCEL *dcel);

void ensureSpaceForFace(struct DCEL *dcel){
    if(! (dcel->faces)){
        dcel->faces = (struct face *) 
                malloc(sizeof(struct face) * INITIALFACES);
        assert(dcel->faces);
        dcel->facesAllocated = INITIALFACES;
    } else if((dcel->facesUsed + 1) > dcel->facesAllocated){
        dcel->faces = (struct face *) realloc(dcel->faces, 
                sizeof(struct face) * dcel->facesAllocated * 2);
        assert(dcel->faces);
        dcel->facesAllocated = dcel->facesAllocated * 2;
    }
}

/* 
    Add an edge from the startVertex index vertex to the endVertex index. Only fills 
    one half-edge as other half-edges will always be added through geometry construction.
*/
void addEdge(struct DCEL *dcel, int startVertex, int endVertex);

/* Add a face to the DCEL given using the given halfEdge and sets the face. */
void addFace(struct DCEL *dcel, struct halfEdge *he);

void addEdge(struct DCEL *dcel, int startVertex, int endVertex){
    ensureSpaceForEdge(dcel);
    
    int newEdge = dcel->edgesUsed;
    
    struct halfEdge *newHE = newHalfEdge();
    newHE->startVertex = startVertex;
    newHE->endVertex = endVertex;
    // newHE->next = NULL;
    // newHE->prev = NULL;
    // newHE->pair = NULL;
    // newHE->face = NOFACE;
    newHE->edge = newEdge;
    
    (dcel->edges)[newEdge].halfEdge = newHE;
    
    dcel->edgesUsed = dcel->edgesUsed + 1;
}

void addFace(struct DCEL *dcel, struct halfEdge *he){
    ensureSpaceForFace(dcel);
    (dcel->faces)[dcel->facesUsed].he = he;
    /* Set the face in the half-edges. */
    he->face = dcel->facesUsed;
    
    (dcel->faces)[dcel->facesUsed].wt = NULL;
    
    struct halfEdge *current = he->next;
    while(current != he){
        current->face = dcel->facesUsed;
        current = current->next;
    }
    
    dcel->facesUsed = dcel->facesUsed + 1;
}

struct DCEL *readPolygonFile(char *polygonfileName){
    struct DCEL *dcel = newDCEL();
    
    FILE *polygonFile = fopen(polygonfileName, "r");
    assert(polygonFile);
    double x;
    double y;
    
    int startVertex = NOVERTEX;
    int endVertex = NOVERTEX;
    
    /* Used to finish off the polygon in the first face. */
    int firstVertex = NOVERTEX;
    int firstEdge = NOEDGE;
    
    while(fscanf(polygonFile, "%lf %lf", &x, &y) == 2){
        ensureSpaceForVertex(dcel);
        (dcel->vertices)[dcel->verticesUsed].x = x;
        (dcel->vertices)[dcel->verticesUsed].y = y;
        dcel->verticesUsed = dcel->verticesUsed + 1;
        if(startVertex == NOVERTEX){
            startVertex = dcel->verticesUsed - 1;
            firstVertex = startVertex;
        } else if(endVertex == NOVERTEX) {
            /* First edge */
            endVertex = dcel->verticesUsed - 1;
            firstEdge = dcel->edgesUsed;
            addEdge(dcel, startVertex, endVertex);
        } else {
            /* Start from last vertex. */
            startVertex = endVertex;
            endVertex = dcel->verticesUsed - 1;
            addEdge(dcel, startVertex, endVertex);
            /* Connect last edge added to newest edge */
            ((dcel->edges)[dcel->edgesUsed - 2].halfEdge)->next = 
                    (dcel->edges)[dcel->edgesUsed - 1].halfEdge;
            /* Connect newest edge to last edge added */
            ((dcel->edges)[dcel->edgesUsed - 1].halfEdge)->prev = 
                    (dcel->edges)[dcel->edgesUsed - 2].halfEdge;
        }
    }
    
    assert(firstEdge != NOEDGE);
    /* Finalise polygon by adding edge back to first vertex. */
    int finalEdge = dcel->edgesUsed;
    addEdge(dcel, endVertex, firstVertex);
    /* Connect previous edge to this edge. */
    ((dcel->edges)[dcel->edgesUsed - 2].halfEdge)->next = 
            (dcel->edges)[dcel->edgesUsed - 1].halfEdge;
    /* Connect newest edge to last edge added */
    ((dcel->edges)[dcel->edgesUsed - 1].halfEdge)->prev = 
            (dcel->edges)[dcel->edgesUsed - 2].halfEdge;
    
    /* Connect final edge back to start edge. */
    ((dcel->edges)[finalEdge].halfEdge)->next = (dcel->edges)[firstEdge].halfEdge;
    /* Connect start edge back to final edge. */
    ((dcel->edges)[firstEdge].halfEdge)->prev = (dcel->edges)[finalEdge].halfEdge;
    
    /* Add face to DCEL - could be any edge we constructed, so may as well be the first. */
    addFace(dcel, (dcel->edges)[firstEdge].halfEdge);
    if(polygonFile){
        fclose(polygonFile);
    }
        
    return dcel;
}

struct split *readNextSplit(FILE *splitfile){
    int firstEdge;
    int secondEdge;
    if(fscanf(splitfile, "%d %d", &firstEdge, &secondEdge) != 2){
        return NULL;
    }
    struct split *split = (struct split *) malloc(sizeof(struct split));
    assert(split);
    split->startEdge = firstEdge;
    split->endEdge = secondEdge;
    split->verticesSpecified = 0;
    split->enforceFace = 0;
    // split->expectedFace = 0;
    return split;
}

void freeSplit(struct split *split){
    if(split){
        free(split);
    }
}

/* Returns 1 if vertices are sufficiently close, 0 otherwise. */
int vertexMatch(struct vertex *v1, struct vertex *v2);

int vertexMatch(struct vertex *v1, struct vertex *v2){
    if(v1->x != v2->x){
        return 0;
    }
    if(v1->y != v2->y){
        return 0;
    }
    return 1;
}

struct halfEdge* applySplit(struct split *split, struct DCEL *dcel){
    int isAdjacent;
    double midpointX;
    double midpointY;
    struct halfEdge *startHE;
    struct halfEdge *endHE;
    struct halfEdge *newJoinHE = NULL;
    struct halfEdge *newJoinHEPair;
    struct halfEdge *newStartHEToMid;
    struct halfEdge *newStartHEToMidPair;
    struct halfEdge *newMidHEToEnd;
    struct halfEdge *newMidHEToEndPair;
    /* Temporary holders for old pair edges */
    struct halfEdge *oldStartPairPrev;
    struct halfEdge *oldEndPairNext;
    /* Temporary holder for old pairs */
    struct halfEdge *oldStartPair;
    struct halfEdge *oldEndPair;
    
    int newVertexMidStart;
    int newVertexMidEnd;
    /* The vertex representing the end of the original starting edge */
    int oldVertexStart;
    /* The vertex representing the start of the original ending edge */
    int oldVertexEnd;
    
    /* Each split creates exactly 3 edges, so we can set up space for these now. */
    int joinEdge;
    int newStartEdge;
    int newEndEdge;
    
    /* Work out what half-edges we need to use. */
    startHE = (dcel->edges)[split->startEdge].halfEdge;
    endHE = (dcel->edges)[split->endEdge].halfEdge;

    /* See if we have the required property before modifying DCEL. */
    /* Uncomment to require face is present in both edges before doing split. */
    if(split->enforceFace){
        if(! (startHE->face == split->expectedFace ||
        (startHE->pair && startHE->pair->face == split->expectedFace)) ||
        ! (endHE->face == split->expectedFace ||
        (endHE->pair && endHE->pair->face == split->expectedFace))
        ){
            int sf = startHE->face;
            int ef = endHE->face;
            int spf = NOFACE;
            int epf = NOFACE;
            if(startHE->pair){
                spf = startHE->pair->face;
            }
            if(endHE->pair){
                epf = endHE->pair->face;
            }
            fprintf(stderr, "Warning: Split was requested with enforceFace ON, "
                            "start edge faces (%d, %d) or end edge faces (%d, %d)"
                            "don't have expected face %d, no split performed\n",
                            sf, spf, ef, epf, split->expectedFace);
            return newJoinHE;
        }
    }

    ensureSpaceForEdge(dcel);
    joinEdge = dcel->edgesUsed;
    dcel->edgesUsed = dcel->edgesUsed + 1;

    ensureSpaceForEdge(dcel);
    newStartEdge = dcel->edgesUsed;
    dcel->edgesUsed = dcel->edgesUsed + 1;

    ensureSpaceForEdge(dcel);
    newEndEdge = dcel->edgesUsed;
    dcel->edgesUsed = dcel->edgesUsed + 1;

    /* Get vertices for MidStart and MidEnd */
    ensureSpaceForVertex(dcel);
    newVertexMidStart = dcel->verticesUsed;
    dcel->verticesUsed = dcel->verticesUsed + 1;

    ensureSpaceForVertex(dcel);
    newVertexMidEnd = dcel->verticesUsed;
    dcel->verticesUsed = dcel->verticesUsed + 1;

    /* Set midpoint of start */
    double startX = (dcel->vertices)[startHE->startVertex].x;
    double startY = (dcel->vertices)[startHE->startVertex].y;
    double endX = (dcel->vertices)[startHE->endVertex].x;
    double endY = (dcel->vertices)[startHE->endVertex].y;
    double startMidX = (startX + endX) / 2.0;
    double startMidY = (startY + endY) / 2.0;

    /* Reusing the start vertex means we don't have to create new half-edges for the
    start side of the split. */
    int startVertexOverlap = MIDEDGE;
    if(split->verticesSpecified){
        /* See if vertex needs to be reused */
        if(vertexMatch(&(dcel->vertices)[startHE->endVertex],
                       &split->startSplitPoint)){
            /* Don't need to create a new vertex for endVertex, so use
                newVertexMidStart as newVertexMidEnd */
            dcel->verticesUsed = dcel->verticesUsed - 1;
            newVertexMidEnd = newVertexMidStart;

            newVertexMidStart = startHE->endVertex;
            startVertexOverlap = ENDOFEDGE;
        } else if(vertexMatch(&(dcel->vertices)[startHE->startVertex],
                              &split->startSplitPoint)) {
            /* Don't need to create a new vertex for endVertex, so use
                newVertexMidStart as newVertexMidEnd */
            dcel->verticesUsed = dcel->verticesUsed - 1;
            newVertexMidEnd = newVertexMidStart;

            newVertexMidStart = startHE->startVertex;
            startVertexOverlap = STARTOFEDGE;
        } else {
            (dcel->vertices)[newVertexMidStart].x = split->startSplitPoint.x;
            (dcel->vertices)[newVertexMidStart].y = split->startSplitPoint.y;
        }
    } else {
        (dcel->vertices)[newVertexMidStart].x = (startX + endX) / 2.0;
        (dcel->vertices)[newVertexMidStart].y = (startY + endY) / 2.0;
    }


    /* Set midpoint of end */
    startX = (dcel->vertices)[endHE->startVertex].x;
    startY = (dcel->vertices)[endHE->startVertex].y;
    endX = (dcel->vertices)[endHE->endVertex].x;
    endY = (dcel->vertices)[endHE->endVertex].y;
    double endMidX = (startX + endX) / 2.0;
    double endMidY = (startY + endY) / 2.0;

    /* Reusing the start vertex means we don't have to create new half-edges for the
    start side of the split. */
    int endVertexOverlap = MIDEDGE;
    if(split->verticesSpecified){
        /* See if vertex needs to be reused */
        if(vertexMatch(&(dcel->vertices)[endHE->startVertex],
                       &split->endSplitPoint)){
            /* Don't need to create a new vertex. */
            dcel->verticesUsed = dcel->verticesUsed - 1;

            newVertexMidEnd = endHE->startVertex;
            endVertexOverlap = STARTOFEDGE;
        } else if(vertexMatch(&(dcel->vertices)[endHE->endVertex],
                              &split->endSplitPoint)){
            /* Don't need to create a new vertex. */
            dcel->verticesUsed = dcel->verticesUsed - 1;

            newVertexMidEnd = endHE->endVertex;
            endVertexOverlap = ENDOFEDGE;
        } else {
            (dcel->vertices)[newVertexMidEnd].x = split->endSplitPoint.x;
            (dcel->vertices)[newVertexMidEnd].y = split->endSplitPoint.y;
        }
    } else {
        (dcel->vertices)[newVertexMidEnd].x = (startX + endX) / 2.0;
        (dcel->vertices)[newVertexMidEnd].y = (startY + endY) / 2.0;
    }


    /* Get point halfway between both midpoints */
    double x1 = startMidX;
    double x2 = endMidX;
    double y1 = startMidY;
    double y2 = endMidY;
    midpointX = (x1 + x2) / 2.0;
    midpointY = (y1 + y2) / 2.0;

    /* Work out whether on correct side. */
    struct vertex *v1 = &((dcel->vertices)[startHE->startVertex]);
    struct vertex *v2 = &((dcel->vertices)[startHE->endVertex]);
    /* Switch this if statement to get face enforcement. */
    //if((getRelativeDir(midpointX, midpointY, v1, v2) == OUTSIDE &&
    //    (! split->enforceFace)) ||
    //   (split->enforceFace && startHE->face != split->expectedFace)){
    if(getRelativeDir(midpointX, midpointY, v1, v2) == OUTSIDE){
        startHE = startHE->pair;
        /* Apply reflections to side of half-edge if not on middle of edge */
        if(startVertexOverlap == ENDOFEDGE){
            startVertexOverlap = STARTOFEDGE;
        } else if(startVertexOverlap == STARTOFEDGE){
            startVertexOverlap = ENDOFEDGE;
        }
    }
    v1 = &((dcel->vertices)[endHE->startVertex]);
    v2 = &((dcel->vertices)[endHE->endVertex]);
    /* Switch this if statement to get face enforcement. */
    //if((getRelativeDir(midpointX, midpointY, v1, v2) == OUTSIDE &&
    //    (! split->enforceFace)) ||
    //    (split->enforceFace && endHE->face != split->expectedFace)){
    if(getRelativeDir(midpointX, midpointY, v1, v2) == OUTSIDE){
        endHE = endHE->pair;
        /* Apply reflections to side of half-edge if not on middle of edge */
        if(endVertexOverlap == ENDOFEDGE){
            endVertexOverlap = STARTOFEDGE;
        } else if(endVertexOverlap == STARTOFEDGE){
            endVertexOverlap = ENDOFEDGE;
        }
    }

    /* Work out whether edges are adjacent. */
    if(startHE->next == endHE){
        isAdjacent = 1;
    } else {
        isAdjacent = 0;
    }

    /* Store old prev and next from start and end edges for convenience */
    struct halfEdge *oldEndPrev = endHE->prev;
    struct halfEdge *oldStartNext = startHE->next;
    oldVertexEnd = endHE->startVertex;
    oldVertexStart = startHE->endVertex;

    /* Update vertices of endHE and startHE */
    endHE->startVertex = newVertexMidEnd;
    startHE->endVertex = newVertexMidStart;

    /* Add bridging edges */
    newJoinHE = newHalfEdge();

    newJoinHE->startVertex = newVertexMidStart;
    newJoinHE->endVertex = newVertexMidEnd;
    newJoinHE->next = endHE;
    endHE->prev = newJoinHE;
    newJoinHE->prev = startHE;
    startHE->next = newJoinHE;
    newJoinHE->pair = NULL; // Will be set later
    /* joinHE is same face as startHE and endHE */
    newJoinHE->face = startHE->face;
    newJoinHE->edge = joinEdge;

    /* Set joinEdge to relevant halfEdge */
    (dcel->edges)[joinEdge].halfEdge = newJoinHE;

    newJoinHEPair = newHalfEdge();
    /* Pair is in opposite direction. */
    newJoinHEPair->startVertex = newVertexMidEnd;
    newJoinHEPair->endVertex = newVertexMidStart;
    newJoinHEPair->next = NULL; // Will join to new HEs
    newJoinHEPair->prev = NULL; // Will join to new HEs
    newJoinHEPair->pair = newJoinHE;
    newJoinHE->pair = newJoinHEPair;
    newJoinHEPair->face = NOFACE; // Will be new face set later
    newJoinHEPair->edge = joinEdge;

    /* Set up what we can of new edges */
    newStartHEToMid = newHalfEdge();
    newStartHEToMid->startVertex = newVertexMidStart;
    newStartHEToMid->endVertex = oldVertexStart;
    newStartHEToMid->next = NULL; // Different setting based on adjacency, set below.
    newStartHEToMid->prev = newJoinHEPair;
    newJoinHEPair->next = newStartHEToMid;
    newStartHEToMid->pair = NULL; // Will be set up later if needed.
    newStartHEToMid->face = NOFACE; // Will be new face set later
    newStartHEToMid->edge = newStartEdge;

    /* Set newStartEdge to relevant halfEdge */
    (dcel->edges)[newStartEdge].halfEdge = newStartHEToMid;

    newMidHEToEnd = newHalfEdge();
    newMidHEToEnd->startVertex = oldVertexEnd;
    newMidHEToEnd->endVertex = newVertexMidEnd;
    newMidHEToEnd->next = newJoinHEPair;
    newJoinHEPair->prev = newMidHEToEnd;
    newMidHEToEnd->prev = NULL; // Different setting based on adjacency, set below.
    newMidHEToEnd->pair = NULL; // Will be set up later if needed.
    newMidHEToEnd->face = NOFACE;
    newMidHEToEnd->edge = newEndEdge;

    /* Set newEndEdge to relevant halfEdge */
    (dcel->edges)[newEndEdge].halfEdge = newMidHEToEnd;

    /* If either start or end HEs have paired Half-Edges, we also need to split those. */
    if(startHE->pair){
        oldStartPairPrev = startHE->pair->prev;
        oldStartPair = startHE->pair;

        newStartHEToMidPair = newHalfEdge();
        /* Reverse of pair */
        newStartHEToMidPair->startVertex = oldVertexStart;
        newStartHEToMidPair->endVertex = newVertexMidStart;
        newStartHEToMidPair->next = oldStartPair;
        newStartHEToMidPair->prev = oldStartPairPrev;
        startHE->pair->prev = newStartHEToMidPair;
        oldStartPair->prev = newStartHEToMidPair;
        oldStartPair->startVertex = newVertexMidStart;
        oldStartPairPrev->next = newStartHEToMidPair;
        newStartHEToMid->pair = newStartHEToMidPair;
        newStartHEToMidPair->pair = newStartHEToMid;
        newStartHEToMidPair->face = startHE->pair->face;
        newStartHEToMidPair->edge = newStartEdge;
    } else {
        newStartHEToMidPair = NULL;
    }
    if(endHE->pair){
        oldEndPairNext = endHE->pair->next;
        oldEndPair = endHE->pair;

        newMidHEToEndPair = newHalfEdge();
        newMidHEToEndPair->startVertex = newVertexMidEnd;
        newMidHEToEndPair->endVertex = oldVertexEnd;
        newMidHEToEndPair->next = oldEndPairNext; // endHE->pair ?
        oldEndPair->next = newMidHEToEndPair;
        oldEndPairNext->prev = newMidHEToEndPair; // Next?
        oldEndPair->endVertex = newVertexMidEnd;
        newMidHEToEndPair->prev = oldEndPair;
        newMidHEToEnd->pair = newMidHEToEndPair;
        newMidHEToEndPair->pair = newMidHEToEnd;
        newMidHEToEndPair->face = endHE->pair->face;
        newMidHEToEndPair->edge = newEndEdge;
    } else {
        newMidHEToEndPair = NULL;
    }

    /* Set up remaining edges. */
    if(isAdjacent){
        newStartHEToMid->next = newMidHEToEnd;
        newMidHEToEnd->prev = newStartHEToMid;
    } else {
        /* Edges are old start and end edges (maybe the same edge). */
        newStartHEToMid->next = oldStartNext;
        oldStartNext->prev = newStartHEToMid;
        newMidHEToEnd->prev = oldEndPrev;
        oldEndPrev->next = newMidHEToEnd;
    }

    /* Setup new face. */
    addFace(dcel, newJoinHEPair);

    /* Check if face has overwritten other face */
    int joinFace = startHE->face;
    if((dcel->faces)[joinFace].he->face != joinFace){
        (dcel->faces)[joinFace].he = startHE;
    }

    /* 20210907 - Added complexity:
        Due to intersection points being possible anywhere along each edge, we have a
        number of added complexities.

        1. The point in the starting edge split may be at the start of the starting edge.
        2. The point in the starting edge split may be at the end of the starting edge.
        3. The point in the end edge split may be at the start of the ending edge.
        4. The point in the end edge split may be at the end of the starting edge.

        This can be handled two ways:
        1. The edges can simply not be created.
        2. The edges can be deleted.

        Of these, (2) presents the cleanest approach
            We have two independent pairs of cases:
            (1a) The point in the starting edge split may be at the end of the starting
                edge:

            In this case we connect the joiningHEPair to startHEtoMid->next and back,
                the edge is then marked NOFACE and deleted, we can verify the face
                pointer. We likewise connect through startHEtoMid->pair in a similar way
                and mark it as NOFACE.

            (1b) The point in the starting edge split may be at the start of the starting
                edge:

            In this case we connect the startHE->prev to joiningHE and back, the edge is
                then marked NOFACE and deleted, we can verify the face pointer. We
                likewise connect through startHE->pair in a similar way and mark it as
                NOFACE, potentially verifying the face pointer.

            (2a) The point in the ending edge split may be at the start of the ending
                edge:

            In this case we connect the endHEtoMid->prev to joiningHEPair and back,
                the edge is then marked NOFACE and deleted, we can verify the face
                pointer. We likewise connect through endHetoMid->pair in a similar way
                and mark it as NOFACE.

            (2b) The point in the ending edge split may be at the end of the ending edge:

            In this case we connect the joiningHE to endHE->next and back, the edge is
                then marked NOFACE and deleted, we can verify the face pointer. We
                likewise connect through endHE->pair in a similar way and mark it as
                NOFACE.

            If both (1b) and (2b) co-occur, we have one additional case:
                if endHE->next is startHE, we instead connect the joiningHE next and prev
                to itself. This creates a zero area face.

            If both (1a) and (2a) co-occur, we have one additional case:
                if endHEtoMid->next is startHEtoMid, we instead connect the
                joiningHE->pair next and prev to itself. This creates a zero area face.
    */
    if(startVertexOverlap == STARTOFEDGE && endVertexOverlap == ENDOFEDGE){
        /* 1b and 2b co-occur */
        struct halfEdge *nextHE = endHE->next;
        if(nextHE == startHE){
            nextHE = newJoinHE;
        }
        struct halfEdge *prevHE = startHE->prev;
        if(prevHE == endHE){
            prevHE = newJoinHE;
        }

        newJoinHE->next = nextHE;
        newJoinHE->prev = prevHE;
        startHE->face = NOFACE;
        endHE->face = NOFACE;
        nextHE->prev = newJoinHE;
        prevHE->next = newJoinHE;
        if((dcel->faces)[joinFace].he->face == NOFACE){
            (dcel->faces)[joinFace].he = newJoinHE;
        }

        int deletedFace;
        /* Handle pairs */
        if(startHE->pair){
            nextHE = startHE->pair->next;
            prevHE = startHE->pair->prev;
            nextHE->prev = prevHE;
            prevHE->next = nextHE;

            deletedFace = startHE->pair->face;
            startHE->pair->face = NOFACE;
            if((dcel->faces)[deletedFace].he->face == NOFACE){
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
        }
        if(endHE->pair){
            nextHE = endHE->pair->next;
            prevHE = endHE->pair->prev;
            nextHE->prev = prevHE;
            prevHE->next = nextHE;

            deletedFace = endHE->pair->face;
            endHE->pair->face = NOFACE;
            if((dcel->faces)[deletedFace].he->face == NOFACE){
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
        }
    } else if(startVertexOverlap == ENDOFEDGE && endVertexOverlap == STARTOFEDGE){
        /* 1a and 2a co-occur */
        struct halfEdge *nextHE = newStartHEToMid->next;
        if(nextHE == newMidHEToEnd){
            nextHE = newJoinHEPair;
        }
        struct halfEdge *prevHE = newMidHEToEnd->prev;
        if(prevHE == newStartHEToMid){
            prevHE = newJoinHEPair;
        }

        newJoinHEPair->next = nextHE;
        newJoinHEPair->prev = prevHE;
        newStartHEToMid->face = NOFACE;
        newMidHEToEnd->face = NOFACE;
        nextHE->prev = newJoinHEPair;
        prevHE->next = newJoinHEPair;
        if((dcel->faces)[newJoinHEPair->face].he->face == NOFACE){
            (dcel->faces)[newJoinHEPair->face].he = newJoinHEPair;
        }
        int deletedFace;
        /* Handle pairs */
        if(newStartHEToMid->pair){
            nextHE = newStartHEToMid->pair->next;
            prevHE = newStartHEToMid->pair->prev;
            nextHE->prev = prevHE;
            prevHE->next = nextHE;

            deletedFace = newStartHEToMid->pair->face;
            newStartHEToMid->pair->face = NOFACE;
            if((dcel->faces)[deletedFace].he->face == NOFACE){
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
        }
        if(newMidHEToEnd->pair){
            nextHE = newMidHEToEnd->pair->next;
            prevHE = newMidHEToEnd->pair->prev;
            nextHE->prev = prevHE;
            prevHE->next = nextHE;

            deletedFace = newMidHEToEnd->pair->face;
            newMidHEToEnd->pair->face = NOFACE;
            if((dcel->faces)[deletedFace].he->face == NOFACE){
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
        }
    } else if(startVertexOverlap != MIDEDGE || endVertexOverlap != MIDEDGE){
        if(startVertexOverlap == ENDOFEDGE){
            /* 1a */
            struct halfEdge *nextHE = newStartHEToMid->next;
            struct halfEdge *prevHE = newStartHEToMid->prev;
            int deletedFace = newStartHEToMid->face;
            prevHE->next = nextHE;
            nextHE->prev = prevHE;
            newStartHEToMid->face = NOFACE;
            if((dcel->faces)[deletedFace].he->face == NOFACE){
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
            if(newStartHEToMid->pair){
                deletedFace = newStartHEToMid->pair->face;
                nextHE = newStartHEToMid->pair->next;
                prevHE = newStartHEToMid->pair->prev;
                newStartHEToMid->pair->face = NOFACE;

                prevHE->next = nextHE;
                nextHE->prev = prevHE;
                if((dcel->faces)[deletedFace].he->face == NOFACE){
                    /* Valid to use either prevHE or nextHE */
                    (dcel->faces)[deletedFace].he = prevHE;
                }
            }
        } else if(startVertexOverlap == STARTOFEDGE){
            /* 1b */
            struct halfEdge *nextHE = startHE->next;
            struct halfEdge *prevHE = startHE->prev;
            int deletedFace = startHE->face;
            prevHE->next = nextHE;
            nextHE->prev = prevHE;
            startHE->face = NOFACE;
            if((dcel->faces)[deletedFace].he->face == NOFACE){
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
            if(startHE->pair){
                deletedFace = startHE->pair->face;
                nextHE = startHE->pair->next;
                prevHE = startHE->pair->prev;
                startHE->pair->face = NOFACE;

                prevHE->next = nextHE;
                nextHE->prev = prevHE;
                if((dcel->faces)[deletedFace].he->face == NOFACE){
                    /* Valid to use either prevHE or nextHE */
                    (dcel->faces)[deletedFace].he = prevHE;
                }
            }
        }
        if(endVertexOverlap == STARTOFEDGE){
            /* 2a */
            struct halfEdge *nextHE = newMidHEToEnd->next;
            struct halfEdge *prevHE = newMidHEToEnd->prev;
            int deletedFace = newMidHEToEnd->face;
            prevHE->next = nextHE;
            nextHE->prev = prevHE;
            newMidHEToEnd->face = NOFACE;
            if((dcel->faces)[deletedFace].he->face == NOFACE){
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
            if(newMidHEToEnd->pair){
                deletedFace = newMidHEToEnd->pair->face;
                nextHE = newMidHEToEnd->pair->next;
                prevHE = newMidHEToEnd->pair->prev;
                newMidHEToEnd->pair->face = NOFACE;

                prevHE->next = nextHE;
                nextHE->prev = prevHE;
                if((dcel->faces)[deletedFace].he->face == NOFACE){
                    /* Valid to use either prevHE or nextHE */
                    (dcel->faces)[deletedFace].he = prevHE;
                }
            }
        } else if(endVertexOverlap == ENDOFEDGE){
            /* 2b */
            struct halfEdge *nextHE = endHE->next;
            struct halfEdge *prevHE = endHE->prev;
            int deletedFace = endHE->face;
            prevHE->next = nextHE;
            nextHE->prev = prevHE;
            endHE->face = NOFACE;
            if((dcel->faces)[deletedFace].he->face == NOFACE){
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
            if(endHE->pair){
                deletedFace = endHE->pair->face;
                nextHE = endHE->pair->next;
                prevHE = endHE->pair->prev;
                endHE->pair->face = NOFACE;

                prevHE->next = nextHE;
                nextHE->prev = prevHE;
                if((dcel->faces)[deletedFace].he->face == NOFACE){
                    /* Valid to use either prevHE or nextHE */
                    (dcel->faces)[deletedFace].he = prevHE;
                }
            }
        }
    }

    return newJoinHE;

}

void freeDCEL(struct DCEL *dcel){
    if(! dcel){
        return;
    }
    int i;
    if(dcel->edges){
        for(i = 0; i < dcel->edgesUsed; i++){
            if((dcel->edges)[i].halfEdge){
                if(((dcel->edges)[i]).halfEdge->pair){
                    /* Free if edge has two halves. */
                    free(((dcel->edges)[i]).halfEdge->pair);
                }
                free(((dcel->edges)[i]).halfEdge);
            }
        }
        free(dcel->edges);
    }
    if(dcel->faces){
        /* All edges are freed above, so no need to free each edge here. */
        free(dcel->faces);
    }
    if(dcel->vertices){
        free(dcel->vertices);
    }
    free(dcel);
}

int getFaceCount(struct DCEL *dcel){
    if(!dcel){
        return 0;
    } else {
        return dcel->facesUsed;
    }
}

int getRelativeDir(double x, double y, struct vertex *v1, struct vertex *v2){
    /* Here we're doing a simple half-plane check against the vector v1->v2. */
    double x1 = v1->x;
    double x2 = v2->x;
    double y1 = v1->y;
    double y2 = v2->y;
    if(x1 == x2 && y1 == y2){
        /* Same point. */
        return DIR_UNDECIDED;
    } else if(x1 == x2){
        /* y = c line */
        /* Work out whether line is going up or down. */
        if(y2 > y1){
            if(x > x1){
                return INSIDE;
            } else if(x < x1){
                return OUTSIDE;
            } else {
                return DIR_UNDECIDED;
            }
        } else {
            if(x < x1){
                return INSIDE;
            } else if(x > x1){
                return OUTSIDE;
            } else {
                return DIR_UNDECIDED;
            }
        }
    } else if(y1 == y2){
        /* x = c line */
        /* Work out whether line is going left or right. */
        if(x2 > x1){
            if(y < y1){
                return INSIDE;
            } else if(y > y1){
                return OUTSIDE;
            } else {
                return DIR_UNDECIDED;
            }
        } else {
            if(y > y1){
                return INSIDE;
            } else if(y < y1){
                return OUTSIDE;
            } else {
                return DIR_UNDECIDED;
            }
        }
    }
    
    /* 
        x1, x2, y1, y2 distinct, so see whether point being tested is
        above or below gradient line.
    */
    double m = (y2 - y1)/(x2 - x1);
    double c = y1 - m*x1;
    
    double predictedY = x * m + c;
    double residual = y - predictedY;
    
    /*
        Being inside or outside the polygon depends on the direction
        the half-edge is going.
    */
    if(x2 > x1){
        if(residual < 0){
            return INSIDE;
        } else if(residual > 0){
            return OUTSIDE;
        } else {
            return DIR_UNDECIDED;
        }
    } else {
        if(residual > 0){
            return INSIDE;
        } else if(residual < 0){
            return OUTSIDE;
        } else {
            return DIR_UNDECIDED;
        }
    }
};

/* 
Takes an established direction and a new direction, returns 1 if
the direction matches the decidedDirection or if the direction is
undecided.
 */
int directionOrUndecided(int decidedDirection, int direction);

int directionOrUndecided(int decidedDirection, int direction){
    if(direction == decidedDirection || direction == DIR_UNDECIDED){
        return 1;
    } else {
        return 0;
    }
}

int inFace(struct DCEL *dcel, double x, double y, int faceIndex){
    if(dcel->facesUsed < faceIndex || ! (dcel->faces)[faceIndex].he){
        return OUTSIDE;
    }
    struct halfEdge *start = (dcel->faces)[faceIndex].he;
    int first = 1;
    int direction = DIR_UNDECIDED;
    
    struct halfEdge *current = start;
    while(start != current || first){
        if(direction == DIR_UNDECIDED){
            /* Doesn't matter where the point is until we find it on one side or the 
            other. */
            direction = getRelativeDir(x, y, &(dcel->vertices)[current->startVertex],
                                       &(dcel->vertices)[current->endVertex]);
        } else {
            if(! directionOrUndecided(direction, 
                                      getRelativeDir(x, y, &(dcel->vertices)[current->startVertex],
                                                     &(dcel->vertices)[current->endVertex]))){
                /* If the point is on the different side of any edge, it be inside 
                    the face, because the face is convex. */
                return 0;
            }
        }
        current = current->next;
        first = 0;
    }
    
    return 1;
}

int getDCELPointCount(struct DCEL *dcel){
    if(!dcel){
        return 0;
    }
    return dcel->verticesUsed;
}

double getDCELVertexX(struct DCEL *dcel, int vertex){
    return (dcel->vertices)[vertex].x;
}

double getDCELVertexY(struct DCEL *dcel, int vertex){
    return (dcel->vertices)[vertex].y;
}

int getDCELEdgeCount(struct DCEL *dcel){
    if(!dcel){
        return 0;
    }
    return dcel->edgesUsed;
}

int getDCELEdgeVertexStart(struct DCEL *dcel, int edge){
    if(!dcel){
        return 0;
    }
    return (dcel->edges)[edge].halfEdge->startVertex;
}

int getDCELEdgeVertexEnd(struct DCEL *dcel, int edge){
    if(!dcel){
        return 0;
    }
    return (dcel->edges)[edge].halfEdge->endVertex;
}

int getDCELEdgeVertexPairStart(struct DCEL *dcel, int edge){
    if(!dcel){
        return 0;
    }
    return (dcel->edges)[edge].halfEdge->pair->startVertex;
}

int getDCELEdgeVertexPairEnd(struct DCEL *dcel, int edge){
    if(!dcel){
        return 0;
    }
    return (dcel->edges)[edge].halfEdge->pair->endVertex;
}

int DCELhasEdge(struct DCEL *dcel, int edge){
    if((dcel->edges)[edge].halfEdge->face != NOFACE){
        return 1;
    } else {
        return 0;
    }
}

int DCELhasEdgePair(struct DCEL *dcel, int edge){
    if((dcel->edges)[edge].halfEdge->face == NOFACE){
        return 0;
    }
    if((dcel->edges)[edge].halfEdge->pair){
        return 1;
    } else {
        return 0;
    }
}

struct intersection *getIntersection(struct bisector *b, struct DCEL *dcel, int face,
        double minLength){

    struct halfEdge *curr = dcel->faces[face].he;
    struct intersection *intersection = (struct intersection*)malloc(sizeof(struct intersection));
    assert(intersection);
    int *tempEdge = &intersection->edge_i, startVert = curr->startVertex, intersectsFound = 0;
    struct vertex intersectPoint, *tempVert = &intersection->vertex_i;

    intersection->multipleIntersects = 0;


    do {

        if (intersects(curr, b, dcel, minLength, &intersectPoint.x, &intersectPoint.y) != DOESNT_INTERSECT) {

            intersectsFound++;

            // store int bool to determine if there were 2 intersections points
            // otherwise there is a single point (not used but might be useful)
            if (intersectsFound == 2) {
                intersection->multipleIntersects = 1;
            }
            // as polygons are strictly convex, there should be max 2 intersections
            else if (intersectsFound >= 2) {
                fprintf(stderr, "too many intersects for face ind %d\n", face);
                break;
            }

            // store intersection point info is either initial or final edge & vertex
            *tempEdge = curr->edge;
            tempVert->x = intersectPoint.x;
            tempVert->y = intersectPoint.y;

            // set temp pointers to final edge & vertex, as there should be max
            // 2 intersections this is sufficient to fill the intersection struct
            tempEdge = &intersection->edge_f;
            tempVert = &intersection->vertex_f;
        }

        curr = curr->next;

    } while(curr->startVertex != startVert);

    return intersection;
}

double getDiameter(struct DCEL *dcel, int faceIndex){

    struct halfEdge *outer = dcel->faces[faceIndex].he;

    if (outer->face != faceIndex) {
        outer = outer->pair;
    }

    struct halfEdge *inner = outer->next;
    struct vertex v1, v2;
    int initStartVert = outer->startVertex;
    double max = 0, distance;

    // pairwise calculate distance between every vertex in the face
    do {

        do {

            v1 = dcel->vertices[outer->startVertex];
            v2 = dcel->vertices[inner->startVertex];

            // store the greatest distance
            if (( distance = sqrt( pow(v2.x - v1.x, 2) + pow(v2.y - v1.y, 2) ) ) > max) {
                max = distance;
            }

            inner = inner->next;

        } while (inner->startVertex != initStartVert);

        outer = outer->next;
        inner = outer->next;

    } while (outer->startVertex != initStartVert);

    return max;
}

void incrementalVoronoi(struct DCEL *dcel, struct watchtowerStruct *wt){

    // case for the first watchtower
    if (dcel->facesUsed == 1 && ! dcel->faces[0].wt) {

        dcel->faces[0].wt = wt;
        wt->face = 0;

    }
    // general case for adding a 2nd or above watchtower, starts by executing
    // a split along the bisector in the face the new watchtower resides
    else {

        // find the face the new watchtower is inside
        int nextFace;
        for(nextFace=0; nextFace<dcel->facesUsed; nextFace++){
            if(inFace(dcel, wt->x, wt->y, nextFace)){
                break;
            }
        }

        // the index of the new face once complete
        int newFace = dcel->facesUsed;

        printf("init face: %d\n", nextFace);

        struct halfEdge *startJoinEdge = NULL;
        struct halfEdge *currJoinEdge = NULL;
        struct halfEdge *nextFaceEdge = NULL;
        struct halfEdge *clockwiseBorderEdge = NULL;
        struct halfEdge *counterClockwiseBorderEdge = NULL;
        int connectsToBorder = 0;
        int i = 0;

        // starting from the face in which the new watchtower is located, execute a split, if the starting edge of
        // the split has a twin, repeat with the face to which that twin edge belongs, this iteration goes clockwise
        // about all faces with this twin edge; if returning to the starting face, break, if it reaches a face where
        // the starting edge has no twin, it has reached the border, continue counterclockwise from starting face
        // until it reaches the border again, than break
        do {

            // if nextFaceEdge is null, we have reached the border, continue splits from starting face but now going
            // counterclockwise, or if we were already doing that just break from the loop
            if (! nextFaceEdge && startJoinEdge != NULL) {
                if (! connectsToBorder) {
                    connectsToBorder = 1;
                    clockwiseBorderEdge = currJoinEdge->next;
                    nextFaceEdge = startJoinEdge->prev->pair;
                    printf("connectsToBorder switched to true\n");
                }
                if (! nextFaceEdge && connectsToBorder) {
                    counterClockwiseBorderEdge = currJoinEdge->prev;
                    printf("exited from border case\n");
                    break;
                }
            }

            // set face index for next split, if statement prevents this from happening on the first loo[
            // where nextFaceEdge is NULL and nextFace has been set as the initial face
            if (nextFaceEdge != NULL) {
                nextFace = nextFaceEdge->face;
            }

            printf("loop iteration: %d, face: %d, connectsToBorder: %d\n",
                   i++, nextFace, connectsToBorder);

            // split along bisector between intersects
            currJoinEdge = executeBisectorIntersectsSplit(dcel, wt, nextFace)->pair;

            printf("split completed succesfully\n");

            // if on first loop, set startJoinEdge (used for border edge case and cleanup)
            if (! startJoinEdge) {
                startJoinEdge = currJoinEdge;
            }

            // set half edge pointer for twin edge of end edge of last split
            if (! connectsToBorder) {
                nextFaceEdge = currJoinEdge->next->pair;
            }
            else {
                nextFaceEdge = currJoinEdge->prev->pair;
            }

        } while(nextFaceEdge != startJoinEdge->prev);

        // linkup new faces and cleanup garbage edges
        if (startJoinEdge == currJoinEdge) {
            // watchtower only interferes with the watchtower of the starting face
            // so no cleanup or linkage required, doesnt account for special case
            // where watchtower does interfere with more than just the watchtower
            // of the starting face, but the bisector only intersects with null
            // (border) half edges
        }
        // general case for 3rd or greater watchtower
        else {

            // for border edge case, clean up border edges in new face
            if (connectsToBorder) {

                int counterClockwiseSkipped = 0;

                do {

                    // vertices on the border created by the splits have straight angles, and after linkage
                    // become redundant, for each of these cases link border half edge to edge after next,
                    // skipped half edge has face set to NOFACE
                    if (clockwiseBorderEdge->next->pair != NULL) {

                        // set redundant internal edges for cleanup
                        clockwiseBorderEdge->next->pair->face = NOFACE;
                        clockwiseBorderEdge->next->face = NOFACE;

                        // link current border edge and border edge after next, set skipped edge for cleanup
                        clockwiseBorderEdge->next = clockwiseBorderEdge->next->pair->next->next;
                        clockwiseBorderEdge->next->prev->face = NOFACE;
                        if (clockwiseBorderEdge->next->prev == counterClockwiseBorderEdge) {
                            counterClockwiseSkipped = 1;
                        }
                        clockwiseBorderEdge->next->prev = clockwiseBorderEdge;
                        clockwiseBorderEdge->endVertex = clockwiseBorderEdge->next->startVertex;
                    }
                    // only set next pointer on loops where no redundant vertex is found, so multiple
                    // collinear redundant vertices can have the edges between them collapsed into one
                    else {
                        clockwiseBorderEdge = clockwiseBorderEdge->next;
                    }

                } while (clockwiseBorderEdge != counterClockwiseBorderEdge && ! counterClockwiseSkipped);
            }

            struct halfEdge *curr = startJoinEdge;
            struct halfEdge *temp = NULL;

            do {

                // if next half edge has a pair which is part of the newly created faces, then link current half edge,
                // which must be the join edge of a split that occurred in this incrementVoronoi call, to the join edge
                // of the next split in a clockwise direction; set redundant internal edges face index to NOFACE
                if (curr->next->pair != NULL && (curr->next->pair->face == NOFACE || curr->next->pair->face - newFace >= 0)) {

                    temp = curr->next;
                    i = 0;

                    // traverse through current face, cleaning up any non newly created internal edges, that is,
                    // edges that were present before any splits that were effected in this incrementVoronoi call
                    do {
                        if (i >= 2) {
                            temp->prev->face = NOFACE;
                        }
                        i++;
                        temp = temp->prev;
                    } while (temp != curr->prev);

                    // set face indices for cleanup
                    curr->next->pair->face = NOFACE;
                    curr->next->face = NOFACE;

                    // linking join edges
                    curr->next = curr->next->pair->next;
                    curr->next->prev = curr;
                }

                // set face index of current edge to newFace
                curr->face = newFace;
                curr = curr->next;

            } while(curr != startJoinEdge);
        }

        // store watchtower in face and vice versa
        dcel->faces[newFace].wt = wt;
        wt->face = newFace;

        // set redundant face half edge pointers to null and cleanup faces
        for (i=dcel->facesUsed-1; i>newFace; i--) {
            dcel->faces[i].he = NULL;
        }
        dcel->facesUsed = newFace + 1;

    }

}

struct halfEdge* executeBisectorIntersectsSplit(struct DCEL *dcel, struct watchtowerStruct *wt, int face){

    // calculate bisector, get intersects, create split struct from intersects
    struct bisector *currBisector = getBisector(dcel->faces[face].wt->x, dcel->faces[face].wt->y, wt->x, wt->y);
    struct intersection *currIntersection = getIntersection(currBisector, dcel, face, DEFAULTMINLENGTH);
    struct split *currSplit = getSplitFromIntersect(currIntersection);
    struct halfEdge* joinEdge;

    // make sure split order results in new watchtower being in new face
    if(getRelativeDir(wt->x, wt->y, &currSplit->startSplitPoint, &currSplit->endSplitPoint) == INSIDE) {
        reverseSplit(currSplit);
    }

    // retrieve new face's join edge for use in incrementVoronoi multi-face linkage
    joinEdge = applySplit(currSplit, dcel);

    freeBisector(currBisector);
    freeIntersection(currIntersection);
    freeSplit(currSplit);

    return joinEdge;
}

struct split *getSplitFromIntersect(struct intersection *intersection){

    struct split *split = (struct split*)malloc(sizeof(struct split));
    assert(split);

    // create split from intersect data
    split->startSplitPoint = intersection->vertex_i;
    split->endSplitPoint = intersection->vertex_f;
    split->startEdge = intersection->edge_i;
    split->endEdge = intersection->edge_f;
    split->verticesSpecified = 1;
    split->enforceFace = 0;
    // split->expectedFace =

    return split;
}

void reverseSplit(struct split* split){

    // reverse start and end edges and points so split happens in reverse
    int tempIndex = split->startEdge;
    struct vertex tempVertex = split->startSplitPoint;

    split->startEdge = split->endEdge;
    split->startSplitPoint = split->endSplitPoint;
    split->endEdge = tempIndex;
    split->endSplitPoint = tempVertex;

}

/////////////////////
// DEBUG FUNCTIONS //
/////////////////////

void printBisector(struct bisector *b){

    // print out bisector data for finite and infinite gradient cases
    if (! b->finiteGrad) {
        printf("bisector with infinite grad, midpoint at %lf, %lf\n", b->mp_x, b->mp_y);
    }
    else {
        printf("bisector with finite grad, grad of %lf, midpoint at %lf, %lf\n", b->grad, b->mp_x, b->mp_y);
    }
}

void printDcel(struct DCEL *dcel) {

    int i;

    printf("\ndcel data:\n\n");

    // print array size for array of vertices, edges and faces
    printf("array size of veritces array: %d\n", (int)dcel->verticesAllocated);
    printf("array size of edges array: %d\n", (int)dcel->edgesAllocated);
    printf("array size of faces array: %d\n", (int)dcel->facesAllocated);

    // print number of faces, edges and vertices in the DCEL
    printf("\n# of faces: %d\n", dcel->facesUsed);
    printf("# of edges: %d\n", dcel->edgesUsed);
    printf("# of vertices: %d\n", dcel->verticesUsed);

    // print info for every face, edge and vertex in the DCEL
    for (i=0; i<dcel->facesUsed; i++) {
        printFace(dcel, i);
    }

    for (i=0; i<dcel->edgesUsed; i++) {
        printEdge(dcel, i);
    }

    for (i=0; i<dcel->verticesUsed; i++) {
        printVertex(dcel, i);
    }
}

// print the indices of every half edge in a face
void printFace(struct DCEL *dcel, int faceIndex) {

    struct halfEdge *currHalfEdge = dcel->faces[faceIndex].he;
    int initVertIndex = currHalfEdge->startVertex;
    printf("\nface index: %d\n\n", faceIndex);

    // print every half edge in a face
    do {
        printHalfEdge(dcel, currHalfEdge);
        currHalfEdge = currHalfEdge->next;
        printf("\n");
    } while (currHalfEdge->startVertex != initVertIndex);
}

// print the indices in a half edge
void printHalfEdge(struct DCEL *dcel, struct halfEdge *halfEdge) {

    printf("face index %d\n", halfEdge->face);
    printf("edge index: %d\n", halfEdge->edge);
    printf("starting "); printVertex(dcel, halfEdge->startVertex);
    printf("ending "); printVertex(dcel, halfEdge->endVertex);
}

// print all half edges at a given edge index
void printEdge(struct DCEL *dcel, int edgeIndex) {

    printf("edge at index %d has half edges:\n\n", edgeIndex);
    printHalfEdge(dcel, dcel->edges[edgeIndex].halfEdge);
    printf("\n");

    // print twin half edge if it exists
    if (dcel->edges[edgeIndex].halfEdge->pair != NULL) {
        printHalfEdge(dcel, dcel->edges[edgeIndex].halfEdge->pair);
        printf("\n");
    }
}

// print coordinates of a vertex
void printVertex(struct DCEL *dcel, int vertexIndex) {

    printf("vertex x,y = %lf,%lf at index: %d\n",
           dcel->vertices[vertexIndex].x,
           dcel->vertices[vertexIndex].y,
           vertexIndex);
}

void printSplit(struct split *split) {
    printf("from edge %d to %d\n", split->startEdge, split->endEdge);
    printf("starting vertex %lf, %lf\n", split->startSplitPoint.x, split->startSplitPoint.y);
    printf("ending vertex %lf, %lf\n", split->endSplitPoint.x, split->endSplitPoint.y);
}