/* 
This is a C implementation of a KDtree, used in neighbour searches.
It is translated from c++, the original code can be found in the
Numerical Recipees, third edition
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include "argsort.h"
#define BIG 1e99


double inline square(double a) {return a*a;}

typedef struct Point 
{
    int DIM;
    double *x;
}Point;

typedef struct Box
// Base element of the tree
{
    Point hi,lo;
    int mom, dau1, dau2, ptlo, pthi;
}Box;

typedef struct Tree
{
    int N;         // Number of points in the tree
    int DIM;       // Dimension of the system
    Point *pts;    // Array of Points with coordinates
    int *ptindx;   // Keep track of the position of each point
    int *rptindx;  // reverse index: rptindx[ptindx[i]]=i
    int nboxes;    // Total number of boxes in the tree
    Box *boxes;    // Array of boxes
}Tree;

Point create_point(int DIM, double *coords);
Point *build_points(double *coords, int N, int DIM);
Box init_Box(Point lo, Point hi, int mom, int dau1, int dau2, int ptlo, int pthi);
void swap(int *array, int i, int j);
double dist(Point a, Point b);
double distBox(Box b, Point p);
double disti(Tree tree, int j, int k);
int selecti(const int k, int *indx, int n, double *arr);
void sift_down(double *heap, int *ndx, int nn);
Tree KDtree(double *coords, int N, int DIM);
int locatePoint(Tree tree, Point pt);
int locateMember(Tree tree, int jpt);
int nearest(Tree tree, Point pt);
void nnearest(Tree tree, int jpt, int *nn, double *dn, int n);
int locatenear(Tree tree, Point pt, double r, int *list, int nmax);


//------------------------------------------------------------------------------------------------
//--------------------   Structure management: creation, initialization       --------------------
//------------------------------------------------------------------------------------------------


Point create_point(int DIM, double *coords)
// Creates a struct point from a double array of its coordinates
// ----> If a point that does not end up in the tree is created, its
//       coordinates should be freed afterwards: free(pt.x)
{
    Point pt;
    int i;
    pt.DIM=DIM;
    pt.x=malloc(DIM*sizeof(double));
    for (i=0;i<DIM;i++)   pt.x[i]=coords[i];
    return pt;
}

Point *build_points(double *coords, int N, int DIM)
// Creates an array of struct points from a 1D array of coordinates
{
    int i,j;
    double *tmp_coords=malloc(DIM*sizeof(double));
    Point *pts=malloc(N*sizeof(Point));
    for(i=0;i<N;i++){
        for (j=0;j<DIM;j++) tmp_coords[j]=coords[N*j+i];   
        pts[i]=create_point(DIM, tmp_coords);
    }
    return pts;
}


Box init_Box(Point lo, Point hi, int mom, int dau1, int dau2, int ptlo, int pthi)
// Creates a box.
{
    int DIM=lo.DIM;
    Box b;
    b.lo=create_point(DIM,lo.x);
    b.hi=create_point(DIM,hi.x);
    b.mom=mom; b.dau1=dau1; b.dau2=dau2; b.ptlo=ptlo; b.pthi=pthi;
    return b;
}



//------------------------------------------------------------------------------------------------
//-----------------------------   Distances functions       ----------------------------------
//------------------------------------------------------------------------------------------------



double dist(Point a, Point b)
// Gives distance between 2 points
{
    int i, DIM=a.DIM;
    double dist=0;
    for(i=0;i<DIM;i++)  dist+=square(a.x[i]-b.x[i]);
    return sqrt(dist);
}

double distBox(Box b, Point p)
// Gives distance of a point from a box
{
    int DIM=p.DIM;
    double dd = 0;
    int i;
    for (i=0; i<DIM; i++){
        if (p.x[i]<b.lo.x[i]) dd += square(p.x[i]-b.lo.x[i]);
        if (p.x[i]>b.hi.x[i]) dd += square(p.x[i]-b.hi.x[i]);
    }
    return sqrt(dd);
}

double disti(Tree tree, int j, int k) 
//Return distance between point j and k in the tree
{
    if (j == k) return BIG;
    else return dist(tree.pts[j], tree.pts[k]);
}



//------------------------------------------------------------------------------------------------
//-----------------------------   Array manipulation: swap, sort      ----------------------------
//------------------------------------------------------------------------------------------------


void swap(int *array, int i, int j)
// swap elements i and j in array.
{
    int tmp;
    tmp=array[i];
    array[i]=array[j];
    array[j]=tmp;
}

int selecti(const int k, int *indx, int n, double *arr)
// All the work of the tree is done here: the array arr is provided, then
// partially sorted through an index array indx (arr itself is untouched).
// At the end of the routine:  arr[indx[0..k-1]] < arr[indx[k]] < arr[indx[k+1..n-1]]  
{
    int i,ia,ir,j,l,mid;
    double a;

    l=0;
    ir=n-1;
    for (;;) {
        if (ir <= l+1) {
            if (ir == l+1 && arr[indx[ir]] < arr[indx[l]])
                swap(indx,l,ir);
            return indx[k];
        } else {
            mid=(l+ir) >> 1;
            swap(indx,mid,l+1);
            if (arr[indx[l]] > arr[indx[ir]]) swap(indx,l,ir);
            if (arr[indx[l+1]] > arr[indx[ir]]) swap(indx,l+1,ir);
            if (arr[indx[l]] > arr[indx[l+1]]) swap(indx,l,l+1);
            i=l+1;
            j=ir;
            ia = indx[l+1];
            a=arr[ia];
            for (;;) {
                do i++; while (arr[indx[i]] < a);
                do j--; while (arr[indx[j]] > a);
                if (j < i) break;
                swap(indx,i,j);
            }
            indx[l+1]=indx[j];
            indx[j]=ia;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
        }
    }
}

void sift_down(double *heap, int *ndx, int nn) 
// This is a sorting routine, related to the heap-sort algorithm.
// It is used by the nnearest routine
{
    int n = nn - 1;
    int j,jold,ia;
    double a;
    a = heap[0];
    ia = ndx[0];
    jold = 0;
    j = 1;
    while (j <= n) {
        if (j < n && heap[j] < heap[j+1]) j++;
        if (a >= heap[j]) break;
        heap[jold] = heap[j];
        ndx[jold] = ndx[j];
        jold = j;
        j = 2*j + 1;
    }
    heap[jold] = a;
    ndx[jold] = ia;
}



//------------------------------------------------------------------------------------------------
//-----------------         MAIN FUNCTION: Kdtree building          ------------------------------
//------------------------------------------------------------------------------------------------
 

Tree KDtree(double *coords, int N, int DIM)
// Build a KD tree from a 1D double array containing all the coordinates: [x[...],y[...],..etc]
// The dimension and total number of points have to be provided.
{
    int jtmp,i,nboxes,ntmp,m,k,kk,j,nowtask,jbox,np,tmom,tdim,ptlo,pthi;
    int *hp;
    double *cp;
    int taskmom[50], taskdim[50];
    int *ptindx=malloc(N*sizeof(int));
    int *rptindx=malloc(N*sizeof(int));
    Point *pts = build_points(coords,N,DIM);

    for (k=0; k<N; k++) ptindx[k] = k;
    m = 1;
    for (ntmp = N; ntmp; ntmp >>= 1) {
        m <<= 1;
    }
    nboxes = 2*N - (m >> 1);
    if (m < nboxes) nboxes = m;
    nboxes--;
    Box *boxes = malloc(nboxes*sizeof(Box));
    for (j=0, kk=0; j<DIM; j++, kk += N) {
        for (k=0; k<N; k++) coords[kk+k] = pts[k].x[j];
    }
    
    Point lo,hi;
    lo.DIM=hi.DIM=DIM;
    lo.x=malloc(DIM*sizeof(double));
    hi.x=malloc(DIM*sizeof(double));
    for(i=0;i<DIM;i++) { lo.x[i]=-BIG; hi.x[i]=BIG; }
    
    boxes[0]=init_Box(lo,hi,0,0,0,0,N-1);
    jbox = 0;
    taskmom[1] = 0;
    taskdim[1] = 0;
    nowtask = 1;
    while (nowtask) {
        tmom = taskmom[nowtask];
        tdim = taskdim[nowtask--];
        ptlo = boxes[tmom].ptlo;
        pthi = boxes[tmom].pthi;
        hp = &ptindx[ptlo];
        cp = &coords[tdim*N];
        np = pthi - ptlo + 1;
        kk = (np-1)/2;
        (void) selecti(kk,hp,np,cp);
        hi = create_point(DIM,boxes[tmom].hi.x);
        lo = create_point(DIM,boxes[tmom].lo.x);
        hi.x[tdim] =  coords[tdim*N + hp[kk]];
        lo.x[tdim] =  coords[tdim*N + hp[kk]];
        jbox++;
        boxes[jbox] = init_Box(boxes[tmom].lo,hi,tmom,0,0,ptlo,ptlo+kk);
        jbox++;
        boxes[jbox] = init_Box(lo,boxes[tmom].hi,tmom,0,0,ptlo+kk+1,pthi);
        boxes[tmom].dau1 = jbox-1;
        boxes[tmom].dau2 = jbox;
        if (kk > 1) {
            taskmom[++nowtask] = jbox-1;
            taskdim[nowtask] = (tdim+1) % DIM;
        }
        if (np - kk > 3) {
            taskmom[++nowtask] = jbox;
            taskdim[nowtask] = (tdim+1) % DIM;
        }
    }
    for (j=0; j<N; j++) rptindx[ptindx[j]] = j;
    Tree tree = {N, DIM, pts, ptindx, rptindx, nboxes, boxes};
    
    return tree;
}


//------------------------------------------------------------------------------------------------
//-----------------------------   Locate functions      ----------------------------------
//------------------------------------------------------------------------------------------------


int locatePoint(Tree tree, Point pt) 
// " In which box of the tree is the arbitrary pt lying ?"
{
    int nb,d1,jdim;
    nb = jdim = 0;
    while (tree.boxes[nb].dau1) {
        d1 = tree.boxes[nb].dau1;
        if (pt.x[jdim] <= tree.boxes[d1].hi.x[jdim]) nb=d1;
        else nb=tree.boxes[nb].dau2;
        jdim = ++jdim % tree.DIM;
    }
    return nb;
}

int locateMember(Tree tree, int jpt) 
// " In which box of the tree is the member point lying ?"
{
    int nb,d1,jh;
    jh = tree.rptindx[jpt];
    nb = 0;
    while (tree.boxes[nb].dau1) {
        d1 = tree.boxes[nb].dau1;
        if (jh <= tree.boxes[d1].pthi) nb=d1;
        else nb = tree.boxes[nb].dau2;
    }
    return nb;
}



//------------------------------------------------------------------------------------------------
//------------------------        KDtree applications         ------------------------------------
//------------------------------------------------------------------------------------------------



int nearest(Tree tree, Point pt) 
// Returns the index of the closest neighbour of the arbitrary point x
{
    int i,k,nrst,ntask;
    int task[50];
    double dnrst = BIG, d;
    k = locatePoint(tree,pt);
    for (i=tree.boxes[k].ptlo; i<=tree.boxes[k].pthi; i++) {
        d = dist(tree.pts[tree.ptindx[i]],pt);
        if (d < dnrst) {
            nrst = tree.ptindx[i];
            dnrst = d;
        }
    }
    task[1] = 0;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (distBox(tree.boxes[k],pt) < dnrst) {
            if (tree.boxes[k].dau1) {
                task[++ntask] = tree.boxes[k].dau1;
                task[++ntask] = tree.boxes[k].dau2;
            } else {
                for (i=tree.boxes[k].ptlo; i<=tree.boxes[k].pthi; i++) {
                    d = dist(tree.pts[tree.ptindx[i]],pt);
                    if (d < dnrst) {
                        nrst = tree.ptindx[i];
                        dnrst = d;
                    }
                }
            }
        }
    }
    return nrst;
}



void nnearest(Tree tree, int jpt, int *nn, double *dn, int n)
// finds the n neighbours of the point of index jpt, returns their
// indexes in nn and their distances to jpt in dn
{
    int j,i,k,ntask,kp;
    int task[50];
    double d;
    if (jpt > tree.N -1) {printf("Error, index of point out of bounds.\n"); return;}
    if (n > tree.N-1) { printf("Error: too many neighbors requested\n"); return;}
    for (i=0; i<n; i++) dn[i] = BIG;
    kp = tree.boxes[locateMember(tree,jpt)].mom;
    while (tree.boxes[kp].pthi - tree.boxes[kp].ptlo < n) kp = tree.boxes[kp].mom;
    for (i=tree.boxes[kp].ptlo; i<=tree.boxes[kp].pthi; i++) {
        if (jpt == tree.ptindx[i]) continue;
        d = disti(tree,tree.ptindx[i],jpt);
        if (d < dn[0]) {
            dn[0] = d;
            nn[0] = tree.ptindx[i];
            if (n>1) sift_down(dn,nn,n);
        }
    }
    task[1] = 0;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (k == kp) continue;
        if (distBox(tree.boxes[k],tree.pts[jpt]) < dn[0]) {
            if (tree.boxes[k].dau1) {
                task[++ntask] = tree.boxes[k].dau1;
                task[++ntask] = tree.boxes[k].dau2;
            } else {
                for (i=tree.boxes[k].ptlo; i<=tree.boxes[k].pthi; i++) {
                    d = disti(tree,tree.ptindx[i],jpt);
                    if (d < dn[0]) {
                        dn[0] = d;
                        nn[0] = tree.ptindx[i];
                        if (n>1) sift_down(dn,nn,n);
                    }
                }
            }
        }
    }
    //    We now sort the neighbours from closest to furthest:
    int *idx=malloc(n*sizeof(int));
    sort(n,dn,idx);
    permutation_int(n,idx,1,nn);
    return;
}


int locatenear(Tree tree, Point pt, double r, int *list, int nmax) 
// Given an arbitrary point of coordinates x, returns in list the indexes
// of all members of the tree lying within r of it. nmax is the maximum number of neighbours
// which will be returned.
{
    int k,i,nb,nbold,nret,ntask,jdim,d1,d2;
    int task[50];
    nb = jdim = nret = 0;
    if (r < 0.0) { printf("Error: radius must be nonnegative"); return 1;}
    while (tree.boxes[nb].dau1) {
        nbold = nb;
        d1 = tree.boxes[nb].dau1;
        d2 = tree.boxes[nb].dau2;
        if (pt.x[jdim] + r <= tree.boxes[d1].hi.x[jdim]) nb = d1;
        else if (pt.x[jdim] - r >= tree.boxes[d2].lo.x[jdim]) nb = d2;
        jdim = ++jdim % tree.DIM;
        if (nb == nbold) break;
    }
    task[1] = nb;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (distBox(tree.boxes[k],pt) > r) continue;
        if (tree.boxes[k].dau1) {
            task[++ntask] = tree.boxes[k].dau1;
            task[++ntask] = tree.boxes[k].dau2;
        } else {
			for (i=tree.boxes[k].ptlo; i<=tree.boxes[k].pthi; i++) {
                            if (dist(tree.pts[tree.ptindx[i]],pt) <= r && nret < nmax)
                                list[nret++] = tree.ptindx[i];
                            if (nret == nmax) return nmax;
			}
        }
    }
    return nret;
}



/* int main() */
/* { */
/*     int N=100; */
/*     int DIM=2; */
/*     double *x=malloc(DIM*N*sizeof(double)); */
/*     int i; */
/*     Tree tree; */
/*     srand(time(NULL)); */

/*     for(i=0;i<DIM*N;i++) */
/*     { */
/*         x[i]=(float)rand()/(float)RAND_MAX; */
/*     } */

/*     tree=KDtree(x,N,DIM); */

/*     int j=0,nb=5; */
/*     int *nn=malloc(nb*sizeof(int)); */
/*     double *dist=malloc(nb*sizeof(double)); */
/*     (void)nnearest(tree,j,nn,dist,nb); */
/*     for(i=0;i<nb;i++) printf(" %d ",nn[i]); */
/*     printf("\n"); */

/*     FILE *out=fopen("outcoords","w"); */
/*     for(i=0;i<nb;i++) fprintf(out, "%d ",nn[i]); */
/*     fprintf(out,"\n"); */

/*     for(i=0;i<N;i++){ */
/*         fprintf(out,"%lf  %lf\n",tree.pts[i].x[0],tree.pts[i].x[1]); */
/*     } */

/*     fclose(out); */

/*     printf("box 0 dau1: %d\n",tree.boxes[0].dau1); */


/* } */
