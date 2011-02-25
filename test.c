#include <stdio.h>
#include <stdlib.h>
#include "hilbert.h"

#define NP 32768


void getkeycube(float x,float y,float z, float dx, int *keymin, int *keymax){

  int ix,iy,iz;
  

}



int main(){

  unsigned nd=3;
  unsigned nb=5,nB=4;
  bitmask_t c1[8*sizeof(bitmask_t)],c2[8*sizeof(bitmask_t)];

  int i,j;
  int ci,cj;
  int kmin,kmax;
  bitmask_t imin,imax;
  int level;

  FILE *fp;
  fp=fopen("out.dat","w");

/*   printf("kmin,kmax?\n"); */
/*   scanf("%d %d",&kmin,&kmax); */
/*   imin=kmin; */
/*   imax=kmax; */


  for(i=0;i<NP;i++){
    hilbert_i2c(nd,nb,(bitmask_t)i,c1);
    //printf("%d %d\n",(unsigned)c1[0],(unsigned)c1[1]);
    fprintf(fp,"%d %d %d\n",(unsigned)c1[0],(unsigned)c1[1],(unsigned)c1[2],(i>=kmin)&&(i<=kmax));
  }
  fclose(fp);

/*   int dx; */
/*   printf("lower corner ?\n"); */
/*   scanf("%d%d",&ci,&cj); */
/*   printf("dx ?\n"); */
/*   scanf("%d",&dx); */
  
/*   int min,max; */

/*   c1[0]=ci; */
/*   c1[1]=cj; */
/*   min=(unsigned)(hilbert_c2i(nd,nb,c1)); */
/*   max=min; */

/*   c1[0]=ci; */
/*   c1[1]=cj+dx; */
/*   j=((unsigned)(hilbert_c2i(nd,nb,c1))); */
/*   min=(j<min?j:min); */
/*   max=(j>max?j:max); */

/*   c1[0]=ci+dx; */
/*   c1[1]=cj+dx; */
/*   j=((unsigned)(hilbert_c2i(nd,nb,c1))); */
/*   min=(j<min?j:min); */
/*   max=(j>max?j:max); */

/*   c1[0]=ci+dx; */
/*   c1[1]=cj; */
/*   j=((unsigned)(hilbert_c2i(nd,nb,c1))); */
/*   min=(j<min?j:min); */
/*   max=(j>max?j:max); */

/*   if(((min<=kmax)&&(min>=kmin))||((max<=kmax)&&(max>=kmin))){ */
/*     printf("the square x=%d y=%d dx=%d is flaged\n",ci,cj,dx);} */
/*   else{ */
/*     printf("the square x=%d y=%d dx=%d is NOT flaged\n",ci,cj,dx);} */

  return 0;
}


