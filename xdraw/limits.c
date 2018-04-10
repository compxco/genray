/******************************************************************************
**  NAME      LIMITS.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Read in data from .in and .bin, load structures
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#define O_BINARY 0
#define PATHSEP '/'
#define STATENAME "state.xdraw"

#else
#include <sys\types.h>
#include <sys\stat.h>
#ifdef DOS
#include <graph.h>
#endif
#include <io.h>
#include <conio.h>
#include <dos.h>
#include "xlib.h"
#include "xutil.h"
#define PATHSEP '\\'
#define STATENAME "state.xdr"
#endif

#ifndef O_RDONLY
#define O_RDONLY 0			/* for penning */
#endif

#include "gendefs.h"
#include "curves.h"
#include "xinit.h"
#include "xtools.h"

static int var_exists(int index);

extern LOOP loop[];		/* Loop structure of the data */
extern int nloop;
extern VARIABLE_ID varid[];	/* Variables used for plotting */
extern int nvar;
extern CURVE_SET curveset[];	/* Describes the plots */
extern int ncset;
extern NODE *nodelist;
extern int nnode, inode, ivar, ncount_equal;
extern float *buf;
extern int param[];

int ldebug = 0;
static int getavg=0;

#ifdef USEMAIN
/*-----------------------------------------------------------------------------
|	main (for debugging)
-----------------------------------------------------------------------------*/
void main(int argc, char *argv[])
{
  init(argc,argv);
}
#endif

/*=============================================================================
**                  INITIALIZATION
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	d_abort, s_abort
-----------------------------------------------------------------------------*/
void d_abort(char *s, int i, int j)
{
  printf("Abort: "); printf(s,i,j); printf(".\n"); exit(0);
}

void s_abort(char *s, char *m, int i, int j)
{
  printf("Abort: "); printf(s,m,i,j); printf(".\n"); exit(0);
}
void i_abort(char *s, int i, int j)
{
  printf("%s (variable ", s);
  if (i & LOOPBIT) printf("I");
  printf("%d in curve %d)\n", s, i&INDBITS, j);
  printf("Warning (may be relevant):\n");
  printf("    Graph titles that begin with a number (0-9) or with the\n");
  printf("    letter I (or i) must be prefixed with a quote (\").\n");
  exit(0);
}

/*-----------------------------------------------------------------------------
|	get_ncount, step_over_ncount
|	* e.g. if loop struc = (t, psi, sol'n step, var)
|	  and nt=5, npsi=4, t current = 2
|	  then ncurve --> 2 --> 2*4+0 == the number of curves to step over
-----------------------------------------------------------------------------*/
int get_ncount(void)
{
   int i, ncurve;
   LOOP *lp;
   for(i=0,ncurve=0,lp=loop; i<inode; i++,lp++)
      ncurve = ncurve * lp->count + lp->current;
   return ncurve;
}

Long step_over_ncount(int i0, int n)
{
   int i;
   Long off0;
   for(i=i0,off0=0; i<n; i++)
      off0 += nodelist[i].nsep;
   return off0;
}

/*-----------------------------------------------------------------------------
|	get_limits
-----------------------------------------------------------------------------*/
void get_limits(CURVE_SET *cp)
{
  CVAR *vp;
  limits(buf, cp, &cp->ix);
  limits(buf, cp, &cp->iy);
  if (cp->gtype != 'G')
    {
      vp = &cp->iz;
      if (vp->index < 0) vp->index = 0;
      vp->hloop = v_is_header(vp->index);
      limits(buf, cp, &cp->iz);
      if (!(cp->flags & FORCED) ||
	   (vp->min < cp->force && vp->max < cp->force) ||
	   (vp->min > cp->force && vp->max > cp->force))
         cp->force = (vp->min + vp->max) / 2;
    }
}

/*-----------------------------------------------------------------------------
|	limits. Finds the limits to be used for the graph.
|	getavg: used to get weighted average and rms for contour
|               currently always 0
-----------------------------------------------------------------------------*/
void limits(float *buf, CURVE_SET *cp, CVAR *xp)
{
#ifdef UNIX
  float *f;
#else
  float huge *f;
#endif
  float dp, xmin, xmax, v;
  int ifam,istep;
  long off0;
  int ncurve, ncurve0;
  float xsum, xsumsq;
  Int n;
  int i,j,ns,nf;
  CURVE_SET *cpp;
  CVAR *xpp, *ypp, *zpp, *wpp;
  LOOP *ql,*qf,*qs,*qv,*lp;
  static int dcall=-1, call=0;
  static FILE *file=NULL;
  int debugging;

/*----- limits for loop type I,X,A,H -----*/
  i = var_exists(xp->index);
  if (i == nvar) i_abort("Indexed variable not found",
	  xp->index, cp - curveset);
  if (xp->label == NULL) xp->label = (varid+i)->name;

  debugging = (dcall>=0) && ldebug;
  if (debugging)
    {
      if (!file) file = fopen("debug.dat", "wt");
      fprintf(file,"******* Min, max for %s, index ", xp->label);
      if (xp->index & LOOPBIT) fprintf(file,"i");
      fprintf(file,"%d ********\n", xp->index & INDBITS);
      fflush(stdout);
    }

  if (!(xp->index & LOOPBIT)) goto VARLIM;

  i = xp->index & INDBITS;
  if (i < 0 || i >=nloop || i==ivar)
    d_abort("Index %d in curve %d out of bounds",i, cp - curveset);

  ql = loop+i;
  f = ql->extra;
  if (ql->ltype=='I')
    {
      xmin = (float)0;
      xmax = (float)(ql->count - 1);
    }
  else if (ql->ltype=='X')
    {
      xmin = *f;
      xmax = *(f+1);
    }
  else if (ql->ltype=='A')
    {
      xmin = *f;
      xmax = *(f + ql->count -1);
    }
  else if (ql->ltype=='H')
    {
      xp->off0 = 0;			/* 'H' always outer loop, for now */
      xmin = *(buf + xp->off0);
      xmax = *(buf + xp->off0 + ql->sep * (ql->count - 1));
    }
  goto LDONE;

/*----- limits for V -----*/
VARLIM:
  qf = loop + cp->lfaml; nf = qf->count;
  qs = loop + cp->lstep; ns = qs->count;
  qv = loop + ivar;
  qv->current = xp->index - (int)qf->ih0;
  qf->current = qs->current = 0;
  xp->dstep = qs->use_sep ? qs->sep : 0;
  xp->dfaml = qf->use_sep ? qf->sep : 0;

  if (xp->hloop >= 0) nf = xp->hloop + 1;

/*------- get off0, load 'current' ------*/
  off0 = 0;
  ncurve = 1;
  for(i=0,lp=loop; i<nloop; i++,lp++)			/* load 'current' */
      if (cp->param[i] >= 0) lp->current = cp->param[i];

  if (inode == -1)					/* get off0 - no 'N' */
    for(i=0,lp=loop; i<nloop; i++,lp++)
    {
      if ((int)lp->ih0 + (int)lp->hcount > xp->index)  /* if subhead variable..*/
	{                                              /* ..at the current level..*/
	  off0 += (xp->index - (int)lp->ih0); break;   /* ..then done! */
  	}
      off0 += lp->sep * lp->current;
      if (i > 0) off0 += (lp-1)->hcount;
    }
  else							/* get off0 - yes 'N' */
  {
    ncurve = get_ncount();
    off0 = step_over_ncount(0, ncurve);
    off0 += xp->index;
  }
  xp->off0 = off0;

  for(i=0; i<ncset; i++)		/* See if this var already there */
    {					/* ncset= # already entered */
      cpp = curveset+i;
      if (cpp==cp) continue;
      xpp = &cpp->ix;
      ypp = &cpp->iy;
      zpp = &cpp->iz;
      wpp = NULL;
      if      (xpp->index==xp->index) wpp = xpp;
      else if (ypp->index==xp->index) wpp = ypp;
      else if (zpp && zpp->index==xp->index) wpp = zpp;
      else continue;

      if      (cpp->lstep==cp->lstep && cpp->lfaml==cp->lfaml) ;
      else if (cpp->lstep==cp->lfaml && cpp->lfaml==cp->lstep) ;
      else continue;

      for(j=0; j<nloop; j++)
	if (cp->param[j]>=0 && cp->param[j]!=cpp->param[j]) break;
      if (j==nloop)
        {
	  xp->min = xmin = wpp->min;
	  xp->max = xmax = wpp->max;
	  if (cp->gtype != 'G' && cpp->gtype == 'G') goto LDONE;
	  return;
        }
    }
  for(ifam=0; ifam<nf; ifam++)			/* Obtain min, max */
    {
      if (ifam==0) off0=xp->off0;
      else if (inode == -1) off0 += xp->dfaml;
      else
        {
	  ncurve0 = ncurve;
	  loop[cp->lfaml].current++;
	  ncurve = get_ncount();
	  off0 += step_over_ncount(ncurve0, ncurve);
	  ns = nodelist[ncurve].ncount;
	}
      f = addfloat(buf, off0);
      if (debugging) {
	fprintf(file,"Family %d, at %lx, step by %ld\n", ifam,f,xp->dstep);
	fflush(stdout);
        }

      for(istep=0; istep<ns; istep++)
	{
	  v = *f;
	  if (xp->dstep > 0) f += xp->dstep;
	  else f += nodelist[istep].nsep;	/* Kludge */
	  if (debugging && (ifam < 2 || ifam > nf-2)){
	    fprintf(file,"%9.5f ",v);
	    if ((istep%10)==9 || istep==ns-1) fprintf(file,"\n");
	    fflush(stdout);
	    }
	  if (ifam==0 && istep==0)
	    {
	      xmin = xmax = xsum = v;
	      xsumsq = v*v;
	      n=1;
	    }
	  else
	    {
	      if (v < xmin) xmin=v;
	      if (v > xmax) xmax=v;
	      if (getavg)
		{
		  xsum += v;
		  xsumsq += v * v;
		  n++;
		}
	    }
	}
    }
	
/*----- set margins -----*/
LDONE:
  if (debugging)fprintf(file,"Min, max (before margins) = %8.3f, %8.3f\n",
			xmin, xmax);
  if (dcall==call) { dcall=-1; fclose(file); }
  else if (dcall>=0) call++;

  if (xmax != xmin)			/* set margins */
    dp = (float) .05 *(xmax - xmin);
  else if (xmin)
    dp = xmin * (float) .5;	/* min not 0 --> dp = min/2 */
  else
    dp = (float) .5;		/* min yes 0 --> dp = 1/2 */

  if (!getavg && cp->gtype != 'G' && xp != &cp->iz)
    {
      xmin -= dp;
      xmax += dp;
    }

  xp->min = xmin;
  xp->max = xmax;

  if (getavg)
    {
      xp->avg = xsum / n;
      xp->rms = (float) sqrt(xsumsq / n - xp->avg * xp->avg);
    }
}

/*-----------------------------------------------------------------------------
|	var_exists
-----------------------------------------------------------------------------*/
static int var_exists(int index)
{
  int i;
  VARIABLE_ID *vi;
  for(i=0, vi=varid;			/* xp->index must exist in varid[] */
      i<nvar && index != vi->index; i++,vi++) ;
  return i;
}

/*-----------------------------------------------------------------------------
|	getparam
|	param[] initially all 0, explicitly set via -i or edit 'O'
-----------------------------------------------------------------------------*/
void getparam(CURVE_SET *cp)
{
  int i, ixl, iyl, lstep;
  LOOP *ql;

/*------ find appropriate lstep ------*/
  ixl = cp->ix.hloop;
  iyl = cp->iy.hloop;
  if (cp->ix.index & LOOPBIT) ixl = cp->ix.index & INDBITS;
  if (cp->iy.index & LOOPBIT) iyl = cp->iy.index & INDBITS;

  if (ixl<0 && iyl<0)				/* ... Both are variables */
    {
      for(i=0,ql=loop; i<nloop; i++,ql++)
	{
	  if (i==cp->lfaml) continue;
	  if (ql->ltype=='N')			/* For now, REQUIRE that...*/
	    {lstep = i; break; }		/* have an N if V1 vs. V2 */
	}
      if (i==nloop) lstep=nloop-2;
      if (cp->lfaml < 0)
	cp->lfaml = (lstep==0) ? 1 : 0;
    }
  else if (iyl < 0)	lstep = ixl;
  else if (ixl < 0)	lstep = iyl;
  else if (ixl <= iyl)	lstep = ixl;
  else			lstep = iyl;

  if (cp->gtype=='G' && cp->lfaml < 0)		/* KLUDGE! */
    cp->lfaml = (lstep==0) ? 1 : 0;

  if (cp->lstep < 0) cp->lstep = lstep;
  else if (cp->lstep > lstep)
    d_abort("Error in lstep for i%d, curve %d",cp->lstep, cp-curveset);
  if (cp->gtype != 'G' && cp->lfaml < 0)
    cp->lfaml = (lstep==ixl) ? iyl : ixl;

  for(i=0; i<nloop; i++)
    cp->param[i] = param[i];

  cp->param[cp->lstep] = -1;
  if (cp->lfaml >= 0) cp->param[cp->lfaml] = -2;
  cp->param[ivar]      = -3;
}

/*-----------------------------------------------------------------------------
|	initparam: set param[] = 0
|		   cp=NULL if reading .in, non-NULL if editing
-----------------------------------------------------------------------------*/
void initparam(CURVE_SET *cp)
{
  int i;
  for(i=0; i<nloop; i++)
    {
      param[i] = 0;
      if (cp != NULL && cp->param[i] >= 0)
	param[i] = cp->param[i];
    }
}

