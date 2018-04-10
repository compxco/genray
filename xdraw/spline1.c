/***********************************************************************
**  NAME:	spline1.c
**		spline1.c = spline.c saved shortly after beginning to
**		use the fact that nl = nu = 1.  Before removing
**		the top loop in blus()
**  AUTHOR:	Sheryl M. Glasser and Alan H. Glasser
**
**  DESCRIPTION:
**	fits functions to cubic splines.
**	Reference: H. Spaeth, "Spline Algorithms for Curves and Surfaces,"
**	Translated from the German by W. D. Hoskins and H. W. Sager.
**	Utilitas Mathematica Publishing Inc., Winnepeg, 1974.
**	
**  USAGE FOR XDRAW
**	-S option.  Load arrays xs and fs (if 2 parameterized variables,
**	will have to do xs, fs, gs).  call splfit, loads fs1.
**	Then for each pair of coordinates, determine if need
**	intermediate points, if so call spleval for those points.
***********************************************************************/
#include <stdio.h>
#include <malloc.h>

#define Float double
#define ind(i,j) nd*(j)+(i-nl)

int nl=1, nu=1;
/*-----------------------------------------------------------------------------
|	splfit
|	* IN:  arrays xs, fs - size m+1 (0 to m)
|	* OUT: derivitives in fs1
-----------------------------------------------------------------------------*/
void splfit(Float *xs, Float *fs, Float *fs1, int m)
{
  int j,k,endmode;
  Float *wk;
  endmode = 0;

/*-------- allocate & load work array -------*/
  wk = (Float *)calloc(4*m+6, sizeof(Float));	/* 1..4m+5 used */
  splfac(xs, wk, m);

/*-------- compute boundary derivatives (polynomial fit, nearest 4 pts) */
  if(endmode == 0)
    {
      fs1[0] = fs1[m] = 0.;
      k=4*m-2;
      for(j=0; j<=3; j++)
	{
	  fs1[0] += wk[j+k]*fs[j];
	  fs1[m] += wk[j+k+4]*fs[m-j];
	}
    }

/*-------- compute first derivatives, interior -------*/
  for(j=1; j<m; j++)
    fs1[j] = 3 * ((fs[j+1]-fs[j]) * wk[j+1] + (fs[j]-fs[j-1]) * wk[j]);

/*-------- extrapolation boundary conditions --------*/
  if(endmode == 0)
    {
      fs1[1]   = fs1[1]   - fs1[0]/(xs[1]-xs[0]);
      fs1[m-1] = fs1[m-1] - fs1[m]/(xs[m]-xs[m-1]);
    }
/*-------- not-a-knot boundary conditions */
  else
    {
      fs1[1]   = fs1[1] - (2*fs[1]-fs[0]-fs[2])*2*wk[1];
      fs1[m-1] = fs1[m-1]+(2*fs[m-1]-fs[m]-fs[m-2])*2*wk[m];
    }

/*-------- solution --------*/
  blus(wk[m+1], fs1[1], m-1, 1);

/*-------- compute not-a-knot boundary nodes --------*/
  if(endmode == 1)
    {
      fs1[0] = (2*(2*fs[1]-fs[0]-fs[2]) + (fs1[1]+fs1[2])*(xs[2]-xs[1])
		-fs1[1]*(xs[1]-xs[0]))/(xs[1]-xs[0]);
      fs1[m) = (2*(fs[m-2]+fs[m]-2*fs[m-1]) +(fs1[m-1]+fs1[m-2])*(xs[m-1]-xs[m-2])
		-fs1[m-1]*(xs[m]-xs[m-1]))/(xs[m]-xs[m-1]);
    }

  free(wk);
}

/*-----------------------------------------------------------------------------
|	splfac
|	* IN:  xs = 0..m, wk = 1..4*m+5
|	* loads matrix for cubic spline fitting
-----------------------------------------------------------------------------*/
void splfac(Float *xs, Float *wk, int m)
{      
  int j,k,endmode;
  endmode = 0;

/*-------- compute inverse differences --------*/
  for(j=1; j<=m; j++)			/* 1..m */
    {
      wk[j] = 1./(xs[j]-xs[j-1]);
      wk[j] *= wk[j];
    }

/*-------- compute interior matrix --------*/
  for(k=m+2,j=1; j<=m-1; j++, k+=3)	/* m+1..(m+1)+3m-4 */
    {
      wk[k-1] = wk[j];
      wk[k]   = 2*(wk[j]+wk[j+1]);
      wk[k+1] = wk[j+1];
    }

/*-------- coefficients for extrapolation left boundary conditions ---*/
  if(endmode == 0)
    {
      k = 4*m-2;
      wk[k] = (xs[0]*( 3*xs[0]-2*(xs[1]+xs[2]+xs[3]) ) +
	       xs[1]*xs[2]+xs[1]*xs[3]+xs[2]*xs[3]) /
		 ((xs[0]-xs[1])*(xs[0]-xs[2])*(xs[0]-xs[3]));

      wk[k+1] = ((xs[2]-xs[0])*(xs[3]-xs[0])) /
	((xs[1]-xs[0])*(xs[1]-xs[2])*(xs[1]-xs[3]));

      wk[k+2] = ((xs[0]-xs[1])*(xs[3]-xs[0])) /
	((xs[0]-xs[2])*(xs[1]-xs[2])*(xs[3]-xs[2]));

      wk[k+3]=((xs[1]-xs[0])*(xs[2]-xs[0])) /
	((xs[3]-xs[0])*(xs[3]-xs[1])*(xs[3]-xs[2]));

/*-------- coefficients for extrapolation right boundary conditions ---*/
      wk[k+4]=(xs[m]*(3*xs[m]-2*(xs[m-1]+xs[m-2]+xs[m-3])) +
	       xs[m-1]*xs[m-2]+xs[m-1]*xs[m-3]+xs[m-2]*xs[m-3]) /
		 ((xs[m]-xs[m-1])*(xs[m]-xs[m-2])*(xs[m]-xs[m-3]));

      wk[k+5]=((xs[m-2]-xs[m])*(xs[m-3]-xs[m])) /
	((xs[m-1]-xs[m])*(xs[m-1)-xs[m-2])*(xs[m-1)-xs[m-3]));

      wk[k+6]=((xs[m]-xs[m-1))*(xs[m-3]-xs[m])) /
	((xs[m]-xs[m-2])*(xs[m-1)-xs[m-2])*(xs[m-3]-xs[m-2]));

      wk[k+7]=((xs[m-1)-xs[m])*(xs[m-2]-xs[m])) /
	((xs[m-3]-xs[m])*(xs[m-3]-xs[m-1))*(xs[m-3]-xs[m-2]));
    }

/*-------- not-a-knot boundary conditions --------*/
  else
    {
      k = m+2;
      wk[k]   += (xs[2]+xs[0]-2*xs[1])*wk[1);
      wk[k+1] += (xs[2]-xs[1])*wk[1];
      k = k+3*(m-2);
      wk[k]   += (2*xs[m-1]-xs[m-2]-xs[m])*wk[m];
      wk[k-1] += (xs[m-1]-xs[m-2])*wk[m];
    }

/*-------- factor matrix --------*/
  bluf(&wk[m+1], m-1);
}


/*-----------------------------------------------------------------------------
|	spleval
|	* evaluates real cubic spline function
-----------------------------------------------------------------------------*/
void spleval(Float *xs, Float *fs, Float *fs1,
	     Float *f, Float *f1, Float *f2, Float *f3,
	     Float x, int m, int n, int nqty, int mode)

{
  int i,j;
  Float d,z,z1;
  /*save j*/

/*-------- find cubic spline interval --------*/
      j = max(j,0);
      j = min(j,m-1);
      do { j--; } while(x<xs[j] && j>0);
      do { j++; } while(x>=xs[j+1] && j<m-1);

/*-------- evaluate offset and related quantities ------*/
      d = xs[j+1]-xs[j];
      z = (x-xs[j])/d;
      z1 = 1-z;

/*-------- evaluate output values --------*/

  for (i=1; i<=nqty; i++)
    f[i] = fs[j,i]*z1*z1*(3-2*z1) + fs[j+1,i]*z*z*(3-2*z) +
      d*z*z1*(fs1[j,i]*z1-fs1[j+1,i]*z);

  if(mode == 0) return;

/*-------- evaluate first derivatives --------*/

  for(i=1; i<=nqty; i++)
    f1[i] = 6*(fs[j+1,i]-fs[j,i])*z*z1/d +
      fs1[j,i]*z1*(3*z1-2)+fs1[j+1,i]*z*(3*z-2);

  if(mode == 1) return;

/*------ evaluate second derivatives. ------*/

  for(i=1; i<=nqty; i++)
    f2[i] = (6*(fs[j+1,i]-fs[j,i])*(z1-z)/d -
	     fs1[j,i]*(6*z1-2)+fs1[j+1,i]*(6*z-2))/d;

  if(mode == 2) return;

/*-------- evaluate third derivatives ------*/

  for(i=1; i<=nqty; i++)
    f3[i]=(12*(fs[j,i]-fs[j+1,i])/d+6*(fs1[j,i]+fs1[j+1,i]))/(d*d);

}


/*-----------------------------------------------------------------------------
|	bluf
|	* IN: the tri-diagonal array elements in wk
|	* performs band storage LU factorization
-----------------------------------------------------------------------------*/
void bluf(Float *a, int n)
{
  int i,j,k,jmin,jmax;
  int nd, kmin;

/*      real*8 a(-nl:nu,n) ==> a(-1..1, m-1)*/

  nl = nu = 1;
  nd = nu+nl+1;
  for(i=1; i<=n; i++)			/* loop over rows 1..m-1 */
    {
      jmin = max(1-i, -nl);		/* jmin = 0 or -1 */
      jmax = min(n-i, nu);		/* jmax = 0 or 1 */
      for(j=jmin; j<=-1; j++)		/* ........compute lower elements */
	{
          kmin = max(jmin, j-nu);	/* kmin = -1 */
	  for(k=kmin; k<=j-1; k++)	/* (never enter this loop) */
	    a[ind(j,i)] -= a[ind(k,i)]*a[ind(j-k,i+k)];
	  a[ind(j,i)] *= a[ind(0,i+j)];
	}
      kmin = max(jmin, -nu);		/* kmin = 0 */
      for(k=kmin; k<=-1; k++)		/* ........compute diagonal element */
	a[ind(0,i)] -= a[ind(k,i)]*a(ind(-k,i+k));
      a[ind(0,i)] = 1./a[ind(0,i)];
    }
}

/*-----------------------------------------------------------------------------
|	blus
|	* performs band storage LU solution
-----------------------------------------------------------------------------*/
void blus(Float *a, Float *x, int n, int ns)
{
  int i,j,is;
      /*real*8 a(-nl:nu,n),x(n,*)*/

  for(is=1; is<=ns; is++)			/* down sweep */
    {
      for(i=1; i<=n; i++)
	{
	  for(j=max(1-i,-nl); j<=-1; j++)
	    x[i,is] -= a[j,i]*x[i+j,is];
	}
      for(i=n; i>=1; i--)			/* up sweep */
        {
          for(j=1; j<=min(n-i,nu); j++)
            x[i,is] -= a[j,i]*x[i+j,is];
          x[i,is] *= a[0,i];
        }         
    }
}

