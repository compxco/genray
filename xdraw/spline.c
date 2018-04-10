/***********************************************************************
**  NAME:	spline.c
**  AUTHOR:	Sheryl M. Glasser and Alan H. Glasser
**
**  DESCRIPTION:
**	fits functions to cubic splines.
***********************************************************************/
#include <stdio.h>
#include <malloc.h>
#include "gendefs.h"
#include "spline.h"

static void abort_spline(void);
static void fit(float *, float *, float *);
static float fourfit(float *, float *, int);
static void getcoeff(float *x, float *f,
		     float *a, float *b, float *c, float *d);
static void getcoeff2(float *x, float *f, float *by, float *cy);

int splinetype = 0;	     
static int crudespline=0;

#define FZ (float)0
#define F1 (float)1
#define F2 (float)2
#define F3 (float)3
#define F6 (float)6

#define evalf(f0,f1, g0,g1, z,z1, dx) \
   f0 * z1*z1 * (F3-F2*z1) + f1 * z*z * (F3-F2*z) \
   + (g0*z1 - g1*z) * z * z1 * dx
#define evalg(f0,f1, g0,g1, z,z1, dx) \
   (F6*z*z1*(f1-f0)/dx + g0*z1*(F3*z1-F2) + g1*z*(F3*z-F2))
#define evals(f0,f1, g0,g1, z,z1, dx) \
   ((F6*(f1-f0)*(z1-z)/dx -g0*(F6*z1-F2) + g1*(F6*z-F2)))/dx

#define value4(a,b,c,d,dx) a*dx*dx*dx + b*dx*dx + c*dx + d
#define deriv4(a,b,c,d,dx) F3*a*dx*dx + F2*b*dx + c

/*============================================================================
|		Begin Interface Functions
============================================================================*/
static int ns;
static float *x, *u, *f, *g, *t;
static int parameterize;

/*----------------------------------------------------------------------------
|	init_spline
----------------------------------------------------------------------------*/
int init_spline(int n, int param)
{
  int ni;
  parameterize = param;
  if (n <= 2) return 0;	
  if (splinetype<0 || splinetype>1) splinetype=1;
  ns = n;
  ni = parameterize ? 5 : 3;
  x = (float *)calloc(ni*ns, sizeof(float));
  if (x) { f = x+ns; g = f+ns; u = g + ns; t = u + ns; }
  else abort_spline();
  return ( x != NULL );
}

/*----------------------------------------------------------------------------
|	zap_spline
----------------------------------------------------------------------------*/
void zap_spline()
{
  if (x) free(x);
  x = f = g = NULL;
}

/*-----------------------------------------------------------------------------
|	abort_spline
-----------------------------------------------------------------------------*/
static void abort_spline()
{
  static int inited=0;	
  if (!inited)
    xprintf("Abort spline: Insufficient space for allocating arrays\n");
  inited = 1;
  zap_spline();
}

/*----------------------------------------------------------------------------
|	add_to_spline
----------------------------------------------------------------------------*/
void add_to_spline(int i, float xi, float yi)
{
  float dx;
  static float dx0;
  static int inited=0;

  if (i >= ns) return;
  *(x + i)   = xi;
  *(f + i)   = yi;
  if (parameterize) *(t + i) = (float)i;
  if (i>0 && !parameterize)
    {
      dx = *(x+i) - *(x+i-1);
      if (i==1) dx0 = dx;
      else if (dx*dx0 <= (float)0)
        {
	  if (!inited)
	    {
	      xprintf(\
"Function not single valued.  Solving for \
spline on %d out of %d points\n", i, ns);
	      inited=1;
            }
	  ns = i;
	}
    }
  crudespline = (ns<4) ? 1 : splinetype;
}

/*----------------------------------------------------------------------------
|	fit_spline
----------------------------------------------------------------------------*/
void fit_spline()
{
  if (crudespline) return;
  if (ns==3) { abort_spline(); return; }
  if (!parameterize)
    fit(x, f, g);
  else
    {
      fit(t, x, u);
      fit(t, f, g);
    }
}

/*-----------------------------------------------------------------------------
|	fit
-----------------------------------------------------------------------------*/
static void fit( float *x, float *f, float *g)
{
  int i, ni;
  float *z, *u11, *delta, l, d;
  
  ni = ns - 1;					/* nb intervals */
  
  g[0]  = fourfit(x,      f,      0);		/* g on ends from 4-fit */
  g[ni] = fourfit(x+ni-3, f+ni-3, 1);

  u11 = (float *)calloc(2*ni, sizeof(float));	/* allocate work array */
  if (u11==NULL) { abort_spline(); return; }
  delta = u11+ni;				/* work array assignments */
  z = g;
  
  for(i=0; i<ni; i++)				/* (data is 0...ni) */
    {
      d = x[i+1] - x[i];
      delta[i] = F1 / d;
    }
  
  for(i=1; i<ni; i++)
    {
      z[i] = (f[i+1] - f[i]) * delta[i] * delta[i] +
	     (f[i] - f[i-1]) * delta[i-1] * delta[i-1];
      z[i] *= F3;
      if (i==1)    z[i] -= g[0] * delta[i-1];
      if (i==ni-1) z[i] -= g[ni] * delta[i];
      u11[i] = F2 * (delta[i] + delta[i-1]);
      if (i > 1)
        {
	  l = delta[i-1] / u11[i-1];
	  z[i] -= l * z[i-1];
	  u11[i] -= l * delta[i-1];
        }
    }

  g[ni-1] = z[ni-1] / u11[ni-1];
  for (i=ni-2; i>0; i--)
    g[i] = (z[i] - g[i+1]*delta[i]) / u11[i];

  free(u11);
}

/*----------------------------------------------------------------------------
|	fourfit
|	* IN:  x, f = pointers to 1st of 4 points; dx = offset to calculate p
|	* OUT: the derivitive at the desired point
----------------------------------------------------------------------------*/
static float fourfit(float *x, float *f, int right)
{
  float a, b, c, d, dx, fi, gi;
  getcoeff(x+1, f+1, &a, &b, &c, &d);
  dx = right ? *(x+3) : *x;
  dx -= *(x+1);
  fi/* = value4(a,b,c,d,dx)*/ ;		/* DEBUG! */
  gi = deriv4(a,b,c,d,dx);
  return gi;
}

/*----------------------------------------------------------------------------
|	getcoeff
|	*   f = a x**3 + bx**2 + cx + d
|	*   IN: x, f = pointer to the 2nd of 4 points
----------------------------------------------------------------------------*/
static void getcoeff(float *x, float *f,
		     float *a, float *b, float *c, float *d)
{
  float M[9], B[9], det, x0, f0, *p, dx;
  static int which=0;
  int i, i0, j;
#define Det(B)  (B[0]*B[4]*B[8] + B[1]*B[5]*B[6] + B[2]*B[3]*B[7] - \
	         B[2]*B[4]*B[6] - B[5]*B[7]*B[0] - B[8]*B[3]*B[1])

  x0 = *x;
  f0 = *d = *f;
  for(i=0; i<9; i+=3)
    {
      j = (i==0) ? -1 : i/3;
      dx = *(x+j) - x0;
      M[i+0] = dx*dx*dx; M[i+1] = dx*dx; M[i+2] = dx;
    }
  det = Det(M);
  for(i0=0; i0<3; i0++)
    {
      for(i=0; i<9; i++) B[i] = M[i];
      B[i0+0] = *(f-1) - f0;
      B[i0+3] = *(f+1) - f0;
      B[i0+6] = *(f+2) - f0;
      p=a; if (i0==1) p=b; else if (i0==2) p=c;
      *p = Det(B) / det;
    }
  /*printf("%s: %g %g %g %g --> %g %g %g %g\n", (which==0)?"X":"Y",
	 *(x-1), *x, *(x+1), *(x+2), *a, *b, *c, *d);*/
  which = which ? 0 : 1;
}

static void getcoeff2(float *x, float *f, float *by, float *cy)
{
  float dx0, dx1, det;	
  dx0 = x[0] - x[1];
  dx1 = x[2] - x[1];
  det = dx0*dx0 * dx1 - dx1*dx1 * dx0;
  *by = ( (f[0]-f[1])*dx1 - (f[2]-f[1])*dx0 ) / det;
  *cy = ( dx0*dx0*(f[2]-f[1]) - dx1*dx1*(f[0]-f[1]) ) / det;
}

/*----------------------------------------------------------------------------
|	eval_spline
|	* evaluates splines in interval i (i=0..ns-2)
|	* z = j/ninterval, j in 1...ninterval-1
|	* crudespline: 0=use all, 1=4-fit, 2=p's from 4-fit, 3=p's from f's
----------------------------------------------------------------------------*/
int eval_spline(int i, float z, float *xout, float *yout, int first)
{
  int j;
  static float x0, ax,bx,cx,dx, ay,by,cy,dy;
  float dx0, dt, z1;

  if (g==NULL || i>=ns-1) return 0;

  z1 = F1 - z;

  if (!crudespline && !parameterize)
    {
      dx0 = x[i+1] - x[i];
      *xout = x[i] + z * dx0;
      *yout = evalf(f[i], f[i+1], g[i],g[i+1], z,z1, dx0);
    }
  else if (!crudespline && parameterize)
    {
      dx0 = t[i+1] - t[i];
      *xout = evalf(x[i],x[i+1], u[i], u[i+1], z,z1, dx0);
      *yout = evalf(f[i],f[i+1], g[i], g[i+1], z,z1, dx0);
    }

  else				/* x,y from 4-fit */
    {
      if (first && ns==3 && i==0)
        {
          ax = ay = FZ;
	  dx = x[1]; dy = f[1];
	  x0 = x[1];
	  if (!parameterize) getcoeff2(x, f, &cy,&dy);
	  else { getcoeff2(t, x, &cx, &dx); getcoeff2(t, f, &cy, &dy); }
	}
      else if (first && ns>3)
        {
          j = (i==0) ? i+1 : i;			/* j=index to 2nd of 4 */
	  if (i==ns-2) j=i-1;
	  x0 = x[j];
	  if (!parameterize) getcoeff(x+j, f+j, &ay, &by, &cy, &dy);
	  else
	    {
	      getcoeff(t+j, x+j, &ax, &bx, &cx, &dx);
	      getcoeff(t+j, f+j, &ay, &by, &cy, &dy);
	    }
        }
      if (!parameterize)
        {
          *xout = x[i] + z * (x[i+1] - x[i]);
          dt = *xout - x0;
          *yout = value4(ay,by,cy,dy, dt);
        }
      else
        {
	  dt = t[i] + z * (t[i+1] - t[i]);
	  *xout = value4(ax,bx,cx,dx, dt);
	  *yout = value4(ay,by,cy,cy, dt);
	}
    }

  return 1;
}

