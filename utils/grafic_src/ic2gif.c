/*
Write gif file.
*/

#include <stdio.h>
#define BITS_PER_PIXEL	8
#include <math.h>
#include <malloc.h>

main(int argc, char *argv[])


{
	int width,height,depth,i,codesize,fputshort(),compress(),idat;
	int np1,np2,np3,n1c,n2c,n3c,nref,ibc;
	int m1t,m2t,m3t,m1off,m2off,m3off;
	float *delta,*tdelta,dat,dmin,dmax;
	char *image,*timage,tchar;
	float dx,astart,omegam,omegav,h0,xoff,yoff,zoff;
	int hanning;
	double drat,mrat;
	FILE *fp;
	int fd;
#include "cmap.h"

/* Check command line arguments. */
/*
	if (!argv[1]) {
		fprintf(stderr,"usage: raw2gif file1.raw file2.gif\n");
		return(-1);
	}
	if (!argv[2]) {
		fprintf(stderr,"usage: raw2gif file1.raw file2.gif\n");
		return(-1);
	}
*/

/* Read raw file written as fortran binary. */
	fp = fopen(argv[1],"r");
	if (!fp) {
		fprintf(stderr,"Cannot open file %s\n",argv[1]);
		return(-1);
	}
    /* Skip 4 bytes start of record (written by fortran). */
	fseek(fp,4,SEEK_CUR);
	fread(&np1,sizeof(int),1,fp);
	fread(&np2,sizeof(int),1,fp);
	fread(&np3,sizeof(int),1,fp);
	fread(&dx,sizeof(float),1,fp);
	fread(&xoff,sizeof(float),1,fp);
	fread(&yoff,sizeof(float),1,fp);
	fread(&zoff,sizeof(float),1,fp);
//	fread(&m1t,sizeof(int),1,fp);
//	fread(&m2t,sizeof(int),1,fp);
//	fread(&m3t,sizeof(int),1,fp);
//	fread(&m1off,sizeof(int),1,fp);
//	fread(&m2off,sizeof(int),1,fp);
//	fread(&m3off,sizeof(int),1,fp);
//	fread(&hanning,sizeof(int),1,fp);
	fread(&astart,sizeof(float),1,fp);
	fread(&omegam,sizeof(float),1,fp);
	fread(&omegav,sizeof(float),1,fp);
	fread(&h0,sizeof(float),1,fp);
	printf("np1, np2, np3 = %d %d %d\n",np1,np2,np3);
	printf("dx, xoff, yoff, zoff = %g %g %g %g\n",dx,xoff,yoff,zoff);
//	printf("m1t, m2t, m3t = %d %d %d\n",m1t,m2t,m3t);
//	printf("m1off, m2off, m3off = %d %d %d\n",m1off,m2off,m3off);
//	printf("hanning = %d\n",hanning);
	printf("astart, omegam, omegav, h0 = %g %g %g %g\n",astart,omegam,omegav,h0);
	width=np1;
	height=np2;
	printf("width, height= %d %d\n",width,height);
    /* Skip 8 bytes end+start of record (written by fortran). */
	fseek(fp,8,SEEK_CUR);
	printf("Enter plane (0-255): ");
	scanf("%d",&i);

	fseek(fp,i*sizeof(float)*(width*height+2),SEEK_CUR);
	delta = (float *)malloc(sizeof(float)*width*height);
	tdelta=delta;
	image = (char *)malloc((unsigned)(sizeof(char)*width*height));
	timage=image;
	fread(delta,sizeof(float),width*height,fp);
	fclose(fp);
	if (ferror(fp)) {
		fprintf(stderr,"Problem reading file %s\n",argv[1]);
		return(-1);
	}

printf("Finding min, max\n");
/* Find min, max of delta. */
	dmin=1.e20;
	dmax=-1.e20;
	for (i=0;i<width*height;i++) {
		dat=*delta++;
		if (dmin > dat) dmin=dat;
		if (dmax < dat) dmax=dat;
	}
	printf("Min, max delta values = %g %g\n",dmin,dmax);
	delta=tdelta;
	printf("Enter dmax: \n");
	scanf("%g",&dmax);
	dmin=-dmax;
//	dmax=1.0e-4;
//	dmin=1.e-8*dmax;
//	dmin=-11.96;
//	dmax=-dmin;
	mrat=dmax/dmin;
	if (mrat < 0.0) mrat=-mrat;
/* Logarithmic rescaling of the data. */
	for (i=0;i<width*height;i++) {
		idat=0;
		dat=(*delta++);
		idat=(int)255.999*(dat-dmin)/(dmax-dmin);
//		drat=dat/dmin;
//		if (drat < 0.0) drat=-drat;
//		if (drat < dmin) drat=dmin;
//		idat=(int)255.999*log(drat)/log(mrat);
		if (idat > 255) idat=255;
		if (idat < 0) idat=0;
		*image++=(unsigned char)idat;
	}
	image=timage;

/* Open gif file and begin filling it. */

	fp = fopen(argv[2],"w");
	if (!fp) {
		fprintf(stderr,"Cannot open file %s\n",argv[2]);
		return(-1);
	}
	fwrite("GIF87a",1,6,fp);/* GIF signature. */

/* Screen Descriptor section. */

	fputshort(width,fp);	/* Screen width. */
	fputshort(height,fp);	/* Screen height. */
	i = ( 0x80 |  (BITS_PER_PIXEL-1) | ((BITS_PER_PIXEL-1)*16));
	fputc(i,fp);		/* Color descriptor byte. */
	fputc(0,fp);		/* Background color byte. */
	fputc(0,fp);		/* Zero byte to end screen descriptor. */

	for (i=0;i<768;i++)
		fputc(cmap[i],fp);	/* Color map. */

/* Image Descriptor section. */
	fputc(0x2c,fp);		/* Image separator character. */
	i = 0;
	fputshort(i,fp);	/* Image left. */
	fputshort(i,fp);	/* Image top. */
	fputshort(width,fp);	/* Image width. */
	fputshort(height,fp);	/* Image height. */
	fputc(0,fp);		/* Color descriptor byte: 0x00 for global color map, sequential ordering. */

/* Raster data section. */
	codesize = BITS_PER_PIXEL;
	fputc(codesize,fp);	/* Initial code size */
	compress(codesize+1,fp,image,width*height);
	fputc(0,fp);		/* End of data stream */

	fputc(0x3b,fp);		/* GIF terminator byte. */

	if (ferror(fp)) {
		fprintf(stderr,"There was an error writing the GIF file\n");
		return(-1);
	}
	free(image);
	fclose(fp);
	return(0);

}

int fputshort(int word, FILE *fp)

{
/* writes a 16-bit integer in GIF order (LSB first) */
	fputc((word&0xff),fp);
	fputc(((word>>8)&0xff),fp);
}
