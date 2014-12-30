#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#define plop printf("plop\n");
#define plop1 printf("plop1\n");
#define plop2 printf("plop2\n");
#define plop3 printf("plop3\n");
#define plop4 printf("plop4\n");
#define plop5 printf("plop5\n");

 
static int compar (const void *a, const void *b);
void sort(int n, double *arr, int *idx);
void permutation_double(int array_size, int *idx, int count, ...);
void permutation_int(int array_size, int *idx, int count, ...);


double *base_arr;
static int compar (const void *a, const void *b)
{
  int aa = *((int *) a), bb = *((int *) b);
  if (base_arr[aa] < base_arr[bb])
    return -1;
  if (base_arr[aa] == base_arr[bb])
    return 0;
  if (base_arr[aa] > base_arr[bb])
    return 1;
}
 
void sort(int n, double *arr, int *idx)
{
  int i,j;
  double *copy = malloc (sizeof (double) * n);
  for (i = 0; i < n; i++)
    {
      idx[i] = i;
    }
  base_arr = arr;
  qsort (idx, n, sizeof (int), compar);
  for(j=0;j<n;j++)
  {
      copy[j]=arr[j];
  }
  for(j=0;j<n;j++)
  {
      arr[j]=copy[idx[j]];
  }
}


void permutation_double(int array_size, int *idx, int count, ...)
{
    va_list ap;
    int i,j;
    va_start (ap, count);         /* Initialize the argument list. */
    for (i = 0; i < count; i++)
    {
        double *array = va_arg (ap, double*);    /* Get the next argument value. */
        double *copy = malloc(array_size * sizeof(double));
        for(j=0;j<array_size;j++)
            copy[j]=array[j];
        for(j=0;j<array_size;j++)
            array[j]=copy[idx[j]];
    }   
    va_end (ap);                  /* Clean up. */
}

void permutation_int(int array_size, int *idx, int count, ...)
{
    va_list ap;
    int i,j;
    va_start (ap, count);         /* Initialize the argument list. */
    for (i = 0; i < count; i++)
    {
        int *array = va_arg (ap, int*);    /* Get the next argument value. */
        int *copy = malloc(array_size * sizeof(int));
        for(j=0;j<array_size;j++)
            copy[j]=array[j];
        for(j=0;j<array_size;j++)
            array[j]=copy[idx[j]];
    }   
    va_end (ap);                  /* Clean up. */
}



/* int main() */
/* { */

/*     int i,N=10; */
/*     double *x = malloc(N * sizeof(double)); */
/*     double *y = malloc(N * sizeof(double)); */
/*     double *z = malloc(N * sizeof(double)); */
/*     int *idx = malloc(N * sizeof(int)); */
    
/*     srand(time(NULL)); */
/*     for (i=0;i<N;i++) */
/*     {   x[i] = rand(); */
/*         y[i] = i; */
/*         z[i]=i*i; */
/*     } */
    
/*     for (i=0;i<N;i++) */
/*     { */
/*         printf("%lf %lf  %lf\n",x[i],y[i],z[i]); */
/*     } */

/*     sort(N,x,idx); */
/*     permutation(N,idx,2,y,z); */
/*     printf("Sorted.\n"); */

/*     for (i=0;i<N;i++) */
/*     { */
/*         printf("%lf %lf  %lf\n",x[i],y[i],z[i]); */
/*     }     */


/* } */



