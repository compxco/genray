/******************************************************************************
**	NAME		XDATA.C
**	AUTHOR		Sheryl M. Glasser
**
**	DESCRIPTION
**		
**
**	Copyright (c) Toptools SCF 1990.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

void main(void);
void writeheader(long,long,long,float*);
void bufwrite(int);
void dcon(void);
void polar(void);

long *buf;
FILE *file;
char outfile[20]="";
long acode='A',xcode = 'X',icode='I',vcode='V';
long npsi=11, nm=50, nn=4, nv=3, nri=2;
long npsi0=2, nm0=1;
long nchar = 0;			/* chars in label */
long version = 0, flags = 0;
int loops_use_headers=0, loops_inited=0;;

#define psize (p-buf)
/*-----------------------------------------------------------------------------
|	main
-----------------------------------------------------------------------------*/
void main()
{
  char c;
  printf("Enter test data type D=Dcon, P=Polar: ");
  scanf("%c",&c);
  c = toupper(c);
  if (c=='D') dcon();
  else if (c=='P') polar();
  if (*outfile) printf("%s written\n",outfile);
  else printf("Invalid data specification\n");
}

/*-----------------------------------------------------------------------------
|	dcon
-----------------------------------------------------------------------------*/
void dcon()
{
  long ipsi, m, n, v;
  long *p,size, npsi0, nm0;
  float off,height, *q;
  double theta;
  
  strcpy(outfile,"dcon.bin");
  file = fopen(outfile,"wb");
  
  buf = (long *)malloc(5000);
/*----- version=0, nloop=5, no labels, only 3 longs ==> no headers -----*/
  loops_use_headers = 1;
  writeheader(5L, nchar, flags, NULL);
  
  writeheader(icode, npsi, npsi0, NULL);
  writeheader(icode, nm, nm0, NULL);
  writeheader(icode, nn, 0L, NULL);
  writeheader(vcode, nv, 0L, NULL);
  writeheader(icode, nri, 0L, NULL);
  
  size = (long)(nm * (nm0 + nn*nv*nri) + npsi0) * sizeof(long);
  for(ipsi=0L; ipsi<npsi; ipsi++)
    { 
      printf("ipsi=%d\n",ipsi);
      off = (float)ipsi;			/* Offset 0..9 */
      p=buf;
      *p++ = (long)size;			/* Record header */
      q = (float *)p;
      *q++ = (float)ipsi;			/* header variables ipsi */
      *q++ = (float)(1.-exp((double)ipsi - 4.));
      for(m=0; m<nm; m++)
        {
	  *q++ = (float)m;			/* header variable m */
          for(n=0; n<nn; n++)
	    {
	      height = (float)(n+3);		/* Height = 3..6 */
	      height += n*off*.2;
	      for(v=1L; v<=nv; v++)
	        {
		  theta = 2*3.14159 * (double)(v*m) / (double)(nm-1); 
		  *q++ = height * cos(theta) + off;
		  *q++ = height * sin(theta) + off;
		}
	    }
        }
      p = (long *)q;
      *p++ = size;
      bufwrite(psize);
    }
  free(buf);
  fclose(file);
}

/*-----------------------------------------------------------------------------
|	polar
|	* Distance source at (ds,phi) to (r,the) =
|	  dist: (r*ct, r*st) to (ds*cp, ds*sp)
-----------------------------------------------------------------------------*/
void polar()
{
long i,ir,iz,nr,ns,nz,itheta,ntheta,ndata,*p;
float rmax,tmax,f, radius[100], angle[3];
double theta, dtheta, r,dr,ds,ct,st;
double ds1,dds;
double f1,f2,f3;
float fmin,fmax;

strcpy(outfile,"polar.bin");
file = fopen(outfile,"wb");

nr = 15; rmax = (float)3;
ns = 3;
ds1 = 1.;
dds = .1;
nz = 6;

dr = rmax/(float)(nr-1);
for(i=0; i<nr; i++) radius[i] = (float)i * dr;

tmax = (float)(2. * 3.14159);
ntheta = 40;
dtheta = (double)tmax/ntheta;
angle[0] = (float)0;
angle[1] = (float)360;
angle[2] = angle[1] / ntheta;

ndata = nz * nr * ntheta;
buf = (long *)malloc((size_t)((ndata+2) * sizeof(long)));

loops_use_headers = 0;
writeheader(4L, nchar, flags, NULL);

writeheader(icode, nz, 0L, NULL);
writeheader(acode, nr, 0L, radius);
writeheader(xcode, ntheta, 0L, angle);
writeheader(vcode, 1L, 0L, NULL);

p=buf;					/* Data */
*p++ = ndata * sizeof(long);
for(iz=0,ds=ds1; iz<nz; iz++, ds+=dds)
for(ir=0; ir<nr; ir++)
  for(itheta=0; itheta<ntheta; itheta++)
    {
      theta = itheta*dtheta;
      r = ir * dr;
      ct = cos(theta);
      st = sin(theta);
      if (ns == 1)
	f = (float)(100. * exp(-r*r));
      else
	{
	  f1 = -r*r - 2.*r*ds*ct - ds*ds;
	  f2 = -r*r + 2.*r*ds*ct - ds*ds;
	  f3 = -r*r + 2.*r*ds*st - ds*ds;
	  f = (float)(100. * (exp(f1) + exp(f2) + exp(f3)));
	}
      if (ir==0 && itheta==0) fmin=fmax=f;
      else if (f<fmin) fmin=f;
      else if (f>fmax) fmax=f;
      *p++ = *(long *)&f;
    }
*p++ = ndata * sizeof(long);
bufwrite(psize);
free(buf);
fclose(file);
printf("Min, Max = %f %f\n",fmin,fmax);
}

/*-----------------------------------------------------------------------------
|	writeheader
|	n1 = count of loop; nh = # at head of loop
|	1st rec: 4 longs ==> yes headers, 3 ==> no
-----------------------------------------------------------------------------*/
void writeheader(long code, long n1, long nh, float *qa)
{
  long *p, *ql, nbytes;
  float *qf, *q2, *pf;
  
  nbytes = loops_use_headers ? 16L : 12L;	/* count for 1st record */
  if (loops_inited)				/* count for subsequent recs */
    {
      nbytes = loops_use_headers ? 12L : 8L;
      if (code==acode) nbytes += n1 * sizeof(long);
      if (code==xcode) nbytes +=  3 * sizeof(long);
    }
  
  p=buf;
  *p++ = nbytes;
  if (!loops_inited) *p++ = version;
  *p++ = code;					/* 1st rec: # loops */
  *p++ = n1;					/* 1st rec: # label chars */
  if (loops_use_headers) *p++ = nh;		/* 1st rec: flags */

  if (code==xcode) { n1 = 3; code = acode; }
  if (code==acode)
    {
      ql = (long *)qa;
      pf = (float *)p;
      for(qf=qa, q2=qf+n1; qf<q2; pf++, qf++)
	{
	  *p++ = *ql++;
	}
    }
  *p++ = nbytes;
  bufwrite(psize);
  loops_inited = 1;
}

/*-----------------------------------------------------------------------------
|	bufwrite
-----------------------------------------------------------------------------*/
void bufwrite(int n)
{
int i;
char *p,temp;

p = (char *)buf;
#ifdef HOME
for(i=0; i<n; i++,p+=4)
  {
    temp = *p; *p = *(p+3); *(p+3) = temp;
    temp = *(p+1); *(p+1) = *(p+2); *(p+2) = temp;
  }
#endif
fwrite(buf,sizeof(long),n,file);
}
