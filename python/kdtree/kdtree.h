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
