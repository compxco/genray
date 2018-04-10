/******************************************************************************
**  NAME      XCONTOUR.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Handles contour plots for xdraw
**      Notation: "ir" = x-axis, "iz" = y-axis
**
**  Copyright (c) GlassWare 1992.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifdef UNIX
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#define Char short

#else
#include "xlib.h"
#include "xutil.h"
#define Char char
#endif

#include "gendefs.h"

#define XCONTOUR		/* used in xcontour.h */

typedef struct
  {				/* (locate here because hxdraw.h uses) */
    byte ix, iy;
    byte where, used;
    int x, y;
  } POINTS;

typedef struct
  {
    byte ix, iy;
    byte where;
  } FLAG;

typedef struct
  {
    int level;
    Char rsign, asign;
    float ax, bx, cx;
    float ay, by, cy;
  } BIN;

typedef struct
  {
    float a, b, c, g;
    int flag;
  } QCOEFF;

#ifndef DOS
#define POINTS_ POINTS
#define BIN_ BIN
#define NCOLOR ncurve
#define Malloc(i,j) malloc(i*j)

#else
#define POINTS_ POINTS huge
#define BIN_ BIN huge
#define NCOLOR 13
#define Malloc(i,j) halloc(i,j)
#endif

#include "curves.h"
#include "xcontour.h"
#include "xdraw.h"		/* for freadf() */
#include "xtools.h"
#include "setcolor.h"
#include "ps.h"

#define FZ (float)0

static int gtype, polar;	/* WARNING! formerly external, now not set */
static int conlabels;
extern int ps_modeon;
extern int redrawflag;
extern float *buf;
extern LOOP loop[];
extern int nloop;
extern CURVE_SET curveset[];

static float *p00;
static float *p10,*p20,*pm0;
static float *p01,*p02,*p0m;

static CURVE_SET *lastcp, *current_cp;
static long off0, xstep, ystep;
static int mr, mz, show_values=0;

int debug_index = -1;
static int debug1 = 0;		/* set if debug_index >= 0; calls savelevel */
static int debug2 = 0;		/* calls savegrid */
static double pi, twopi, dtheta;
static LOOP *rlp;

/*=============================================================================
**                    SUPPORTING FUNCTIONS
**===========================================================================*/
extern Display *mydisplay;
extern Window redraw_win;
extern GC redraw_gc;
extern float xscale, xoffset, yscale, yoffset;
static float rscale, roffset, zscale, zoffset;

static int isep = -1;
static int ncurve=16, last_ncurve, hzcurve=0;
extern int default_ncurve;
static int nrz, nrza;
static float psi0, dpsi, dpsi0;
static float force;
static int forced = 0;
static int quad_always = 1;
static int i0 = 0, kz0, kzn, kr0, krn;

#define PS_AT           0x1	/* bits for ps->where */
#define PS_RIGHT        0x2
#define PS_ABOVE        0x4
#define PS_FROM         0x1	/* bits for ps->used */
#define PS_TO           0x2
#define PS_USED        (PS_FROM | PS_TO)

static POINTS_ *list;
static BIN_ *bin;

/*-----------------------------------------------------------------------------
|	get_dpsi
-----------------------------------------------------------------------------*/
static void get_dpsi(CURVE_SET *cp, float *zmin, float *zmax)
{
  float min, max, f, delta, df;
  int i, icrit;
  float dist, distmin, psi;
  char text[80];
  
  min = cp->iz.min;
  max = cp->iz.max;
  delta = max - min;
  if (ncurve == 1)
    {
      dpsi = delta;
      psi0 = forced ? cp->force : (min + max) / 2.;
      goto DONE;
    }
  if (!forced)
    {
      df = (ncurve & 1) ? delta / (float)(ncurve+1) : delta / (float)ncurve;
      min += df / 2.;
      max -= df / 2.;
      delta = max - min;
      dpsi = delta / (float)(ncurve - 1); psi0 = min;
      goto DONE;
    }

  df = cp->forcemin * delta;
  if (min < force && min+df <= force) min += df;
  df = cp->forcemax * delta;
  if (max > force && max-df >= force) max -= df;
  
  f = cp->force;			/* ....Forced! */
  if (f < min) f = min;
  if (f > max) f = max;
  
  dpsi = delta / (float)(ncurve-1);		/* trial dpsi */
  icrit=-1; distmin=delta;			/* find which i is closest...*/
  for(i=0; i<ncurve; i++)			/* ... to f */
    {
      psi = min + i * dpsi;
      dist = fabs(psi-f);
      if (dist < distmin) { distmin=dist; icrit=i; }
    }
  if ((icrit == 0 || f > icrit*dpsi) && icrit != ncurve-1)
    dpsi = (max-f)/(float)(ncurve-1-icrit);
  else
    dpsi = (f-min)/(float)icrit;

  psi0 = f - (float)icrit * dpsi;
  
  if (hzcurve)
    {
      dpsi0 = dpsi / (hzcurve + 1);
      hzcurve = 2 * hzcurve + 1;
    }
DONE:
  *zmin = min;
  *zmax = max;
  if (show_values)
    {
       sprintf(text,"Min, max of data: %g,%g\n", cp->iz.min, cp->iz.max);
       xprintf(text);
       sprintf(text,"Min, max used:    %g %g\n", min, max);
       xprintf(text);
       sprintf(text,"%d levels, delta = %g\n", ncurve, dpsi);
       xprintf(text);
    }
}

/*-----------------------------------------------------------------------------
|	redraw1.  draw contours
|	xmin..ymax are limits according to current zoom
-----------------------------------------------------------------------------*/
void redraw1(CURVE_SET *cp, float xmin, float xmax, float ymin, float ymax)
{
  float *p, *pr, *pa, *p0;
  float min, max, psi, slope;
  float dr, dz, fr, fz, rdelt, zdelt, r1a, z1a;
  int count[50];
  CVAR *xp, *yp;
  int i, ii, j, ir, iz, r1, z1, r2, z2;
  static int inited = 0;
  double theta, rad;
  float x,y,dx;
  BIN_ *b;
  BIN_ *b2;
  char text[80];

  current_cp = cp;
  off0 = cp->iz.off0;
  xstep = cp->iz.dfaml;
  ystep = cp->iz.dstep;

  xp = &cp->ix;
  yp = &cp->iy;
  i = xp->index & INDBITS;  mr = loop[i].count;
  i = yp->index & INDBITS;  mz = loop[i].count;
  nrz = mr * mz;

  rlp = &loop[i];
  pi = atan2(0., -1.);			/* Some values in case polar */
  twopi = 2. * pi;
  dx = (float).1;			/* for solid; should reflect mz! */
  dtheta = twopi / (double)mr;

  i = (cp->gtype == 'S') ? 0 : 1;
  fr = (float)(mr-i);
  fz = (float)(mz-i);
  dr = (xp->max - xp->min) / fr;
  dz = (yp->max - yp->min) / fz;

  rdelt = rscale = dr * xscale;		/* Use these ==> ir, iz = type X */
  zdelt = zscale = dz * yscale;
  roffset = xp->min * xscale + xoffset;
  zoffset = yp->min * yscale + yoffset;
  
  gtype = cp->gtype;
  polar = cp->flags & POLAR;
  conlabels = cp->flags & LABELF;
  ncurve = cp->ncurve;
  if (ncurve==0) ncurve = cp->ncurve = default_ncurve;
  forced = cp->flags & FORCED;
  force = cp->force;		/* if not forced, this was set to avg */
  hzcurve = cp->hzcurve;

  /*-------- find index limits for tracing contour ------*/
  kr0 = (int)((xmin - xp->min) * fr / (xp->max - xp->min));  kr0 -= 2;
  krn = (int)((xmax - xp->min) * fr / (xp->max - xp->min));  krn += 2;
  kz0 = (int)((ymax - yp->min) * fz / (yp->max - yp->min));  kz0 -= 2;
  kzn = (int)((ymin - yp->min) * fz / (yp->max - yp->min));  kzn += 2;

  if (kr0 < 0) kr0=0;
  if (kz0 < 0) kz0=0;
  if (krn >mr) krn=mr;
  if (kzn >mz) kzn=mz;

  if (polar) { kzn = mz; krn = mr; }

/*--------- Get contour values ------------*/
  get_dpsi(cp, &min, &max);
						     
/*--------- Allocate, load list[], bin[] -------*/
  if (!inited)				/* Allocate list[] and bin[] */
    {
      nrza = nrz;
      list = (POINTS_ *)Malloc((long)nrz,sizeof(POINTS));
      bin = (BIN_ *)Malloc((long)nrz,sizeof(BIN));
      if (list == NULL || bin == NULL)
        {
	  printf("Unable to allocate auxiliary array%s\n",
		  list ? "." : "s.");
	  exit(1);
        }
    }
#ifdef UNIX
  else if (nrza < nrz)
    {
      nrza = nrz;
      list = (POINTS_ *)realloc(list, (long)nrz * sizeof(POINTS));
      bin = (BIN_ *)realloc(bin, (long)nrz * sizeof(BIN));
    }
#endif

  if (!inited || ncurve!=last_ncurve ||		/* Load bin[] for points */
      cp != lastcp)
    {
      lastcp = cp;
      p0 = buf+off0;
      for (iz=0,b=bin; iz<mz; iz++,p0+=ystep)
      for (ir=0,p=p0; ir<mr; ir++,p+=xstep,b++)
	{
	  if	  (*p < min)	b->level = 0;	/* .....level */
	  else if (*p > max)	b->level = (int) ((max - min) / dpsi);
	  else			b->level = (int) ((*p - min) / dpsi);

	  b->rsign = b->asign = 0;		/* .....rsign and asign */
	  pr = p+xstep;
	  pa = p+ystep;
	  
	  if (ir == mr-1) slope = FZ;
	  else
	    {
	      slope = *pr - *p;
	      if      (slope < FZ) b->rsign = -1;
	      else if (slope > FZ) b->rsign = 1;
	    }
	  if (iz == mz-1) slope = FZ;
	  else
	    {
	      slope = *pa - *p;
	      if      (slope < FZ) b->asign = -1;
	      else if (slope > FZ) b->asign = 1;
	    }
	}
    }
  b2 = bin + nrz;
  inited = 1;
  last_ncurve = ncurve;

/*------ Draw solid area */

  if (gtype == 'S')
    {
      if (ps_modeon)
	return;
      for (i = 0; i <= ncurve; i++)
	count[i] = 0;

      ir = mr - 1;
      iz = -1;
      r1a = xmin * xscale + xoffset;
      z1a = ymax * yscale + yoffset;
      for (b = bin; b < b2; b++)		/* (i=0; i<nrz; i++) */
	{
	  ir++;
	  if (ir >= mr)
	    {
	      ir = 0;
	      iz++;
	      r1 = (int)r1a;
	      fz = z1a + (float)iz * zdelt;
	      z1 = (int)fz;
	      z2 = (int)(fz + zdelt);
	    }
	  else
	    r1 = r2;
	  r2 = (int)(r1a + (float)(ir+1) * rdelt);
	  
	  if (polar)
	    {
	      /*if (iz==0 && ir>0) continue;*/
	      theta = (double)(ir * 2. * 3.14159 / mr);
	      rad = (double)(iz * cp->iz.max / (mz-1));
	      x = (float)(rad * cos(theta));
	      y = (float)(rad * sin(theta));
	      r1 = (int)((x-dx) * xscale + xoffset);
	      r2 = (int)((x+dx) * xscale + xoffset);
	      z1 = (int)((y-dx) * yscale + yoffset);
	      z2 = (int)((y+dx) * yscale + yoffset);
	    }

	  i = b->level;
	  count[i]++;
	  setcolor(i, ncurve, isep);		/* color for solid */

	  XFillRectangle(mydisplay, redraw_win, redraw_gc,
			 r1, z2, (r2 - r1), (z1 - z2));
	}
    }

/*------ Draw curves */

  else
    {
      rand_label(0);
      for (i = 0; i < ncurve; i++)
	{
	  psi = psi0 + i * dpsi;
	  ii = (int) ((psi - min) / dpsi);
	  setcolor(i, ncurve, isep);	/* color for contour */
	  if (debug_index == ii)
	    { debug1 = 1;
	      printf("See debug.dat: savelevel() for psi=%f\n", psi); }
	  if (i == isep && hzcurve)
	    {
	      for (psi -= dpsi - dpsi0, j = 0; j < hzcurve; j++, psi += dpsi0)
		drawcontour(psi, ii);
	      psi -= dpsi;
	      continue;
	    }
	  if (show_values)
	    {
	      sprintf(text, "Contour %d: %g\n", ii, psi);
	      xprintf(text);
	    }
	  drawcontour(psi, ii);
	  debug1 = 0;
	}
    }
  show_values = 0;
}

/*-----------------------------------------------------------------------------
|	drawcontour
-----------------------------------------------------------------------------*/
void drawcontour(float psi, int ii)
{
#ifdef UNIX
  POINTS *p, *p1, *p2, *pnow, *pfrom1, *pfrom2, *q, *qln, *qlf, *qlft;
  POINTS *qlf0, *qlft0, *psave, *qsave1, *qsave2;
#else
  POINTS huge *p,huge *p1, huge *p2, huge *pnow, huge *pfrom1, huge *pfrom2;
  POINTS huge *q,huge *qln,huge *qlf,huge *qlft;
  POINTS huge *qlf0,huge *qlft0,huge *psave,huge *qsave1,huge *qsave2;
#endif
  int ir, iz, jr, jz, first, pass;
  int n, nvec, ivec, labeli, lx, ly;
  int ind, ind2, di, dj, dii;
  byte where, used;
  char text[80];
  FLAG ilist[8], *il, *il2;
  BIN_ *b;
  BIN_ *b0;

  /*printf("Contour %d: %g\n", ii, psi);*/
  b0 = bin + kz0 * mr;
  for (iz=kz0,p=list; iz<kzn; iz++,b0+= mr)	/* set flag for all cells */
    {						/* ir,iz containing psi */
      for (ir = kr0; ir < krn; ir++)
	{
	  b = b0 + ir;
	  if (debug2)
	    savegrid(ir, iz);
	  dii = b->level - ii;
	  if (abs(dii) <= 2) goto MAYBE;
	  if (ir < mr-1 && ((b+1)->level-ii)*dii <= 0)
	    goto MAYBE;
	  if (iz < mz-1 && ((b+mr)->level-ii)*dii <= 0)
	    goto MAYBE;
	  continue;
	MAYBE:
	  p = has_psi(psi, ir, iz, p1=p);	/* load into list */
	}
    }
  p2 = p;
  nvec = p2 - list;
  ivec = -1;
  labeli = rand_label(nvec);
  
  if (debug1)
    savelevel(1, p2 - list, NULL);	/* Open file, print entire list */

  if (gtype == 'D')			/* Dots --> draw points in list, */
    {					/* don't try to connect */
      for (p = list; p < p2; p++)
	{
	  jr = p->x;
	  jz = p->y;
	  if (!ps_modeon)
	    XDrawLine(mydisplay, redraw_win, redraw_gc, jr, jz, jr, jz);
	  else
	    {
	      ps_moveto(jr, jz);
	      ps_lineto(jr, jz);
	    }
	}
      if (ps_modeon)
	ps_stroke();
      return;
    }

/*----------- list[] = all gridpoints bordering psi ----*/
/*            next line starts with q: a totally unused point */
  quad_always = (gtype != 'L');
  quad_always = 0;
  p1 = list - 1;
START:					/* Get 1st point on a contour */
  for (q=p1+1; q<p2 && q->used; q++);	/* = totally unused point */
  if (q >= p2)				/* If none found, done! */
    {
      if (debug1)
	savelevel(0, p2 - list, NULL);	/* print new list, close */
      if (conlabels)
        {
	  sprintf(text, "%g", psi);
	  XDrawString(mydisplay, redraw_win, redraw_gc, lx, ly,
		      text, strlen(text));
	}
      return;
    }
  if (debug1)
    savelevel(2, 0, "New line\n");

  qlf0 = qlft0 = 0;
  pass = 0;
  first = -1;
/*--------------- Seek next point q ------*/
  for (p1=q,pnow=pfrom1=0;;)
    {
      qln = qlf = qlft = 0;
      pfrom2 = pfrom1;
      pfrom1 = pnow;
      pnow = q;						/* pnow=last pt drawn */
      nextpoint((pass==3)?0:first, pnow, pfrom1);	/* init angle test */

      il2 = surround(ilist, pnow);	/* ilist[]=rel indices of...*/
      if (debug1)			/* ...points surrounding pnow */
	savelevel(2, pnow - list, "Line from %d ");
      n = il2 - ilist;
      if (ps_modeon && first)
	ps_moveto(pnow->x, pnow->y);

      q = (first || pass==3) ? list : p1;
      for (; q < p2; q++)		/* Search for all q's that are... */
	{				/* ..reasonably nearby (di,dj ok) */
	  if (q==pnow || q==pfrom1 || q==pfrom2)
	    continue;
	  dj = q->iy - pnow->iy;	/* get q = likely candidate */
	  if (dj < -1) continue;	/* try only i,j in -1..1 */
	  if (dj > 1)  break;
	  di = q->ix - pnow->ix;
	  if (polar && (q->ix==0 || pnow->ix==0) &&
	      (di==mr-1 || di==1-mr));
	  else if (di < -1 || di > 1)
	    continue;
	  if (debug1)
	    savelevel(2, q - list, " Try %d");

	  ind = q->ix + 256 * q->iy;
	  where = q->where;
	  used = q->used;
	  for (il = ilist; il < il2; il++)	/* find it in ilist[]... */
	    {
	      ind2 = il->ix + 256 * il->iy;
	      if (ind != ind2)
		continue;
	      if (!(where & il->where))
		continue;
	      if (!used)			/* ...and set qln etc */
		qln = nextpoint(1, q, qln);	/* ..based on "empty"..*/
	      else if (used != PS_USED)		/* ..or min angle */
		qlf = nextpoint(2, q, qlf);
	      else
		qlft = nextpoint(3, q, qlft);
	    }
	}

      if (first)
	{
	  qlf0 = qlf;
	  qlft0 = qlft;
	}
      if (qln)
	q = qln;		/* connect to unused pt */
      else if (qlf)
	q = qlf;		/* connect- close contour */
      else if (qlft)
	q = qlft;		/* connect to prior contour */
      else			/* no connection: check if bkwrd from... */
	{			/* beginning, if not break out */
	  if (debug1)
	    savelevel(2, 0, "\nBreak\n");
	  goto ONE_MORE;
	}
      if (debug1)
	savelevel(2, q-list, " Using %d\n");
      pnow->used |= PS_FROM;
      q->used |= PS_TO;
      if (!ps_modeon)
        {
	  XDrawLine(mydisplay, redraw_win, redraw_gc,
		    pnow->x, pnow->y, q->x, q->y);
	  ivec++;
	  if (conlabels && (ivec==labeli || ivec==0)) { lx = q->x; ly = q->y; }
        }
      else
	ps_lineto(q->x, q->y);

      if (pass == 0)			/* Save info re:1st pt, reverse dir */
	{
	  psave = pnow;
	  qsave1 = q;
	  qsave2 = 0;
	  pass = 1;
	}
      else if (pass == 1)
	{
	  pass = 2;
	  qsave2 = q;
	}

      first = 0;
      if (!qln)			/* done if connect was already used */
	{
	ONE_MORE:		/* but first see if start of line connects too */
	  if (pass == 3) break;
	  if (first)
	    {
	      XDrawPoint(mydisplay, redraw_win, redraw_gc, pnow->x, pnow->y);
	      break;
	    }
	  pass = 3;		/* pass=3 ==> reverse of original direction */
	  q = psave;
	  pnow = qsave1;
	  pfrom1 = qsave2;
	  first = -1;
	}
    }
  if (ps_modeon)
    ps_stroke();
  goto START;
}

/*-----------------------------------------------------------------------------
|	nextpoint
|	* return q of possible next point with minimum angle
|	* k=-1 or 0 ==> initialize ("first" true or false)
|	* k=1,2,3 ==> e.g. closing contour (qln,qlf,qlft)  
-----------------------------------------------------------------------------*/
POINTS_ *nextpoint(int k, POINTS_ * q, POINTS_ * qln)
{
  static double diff[3], angle1[3], angle2[3];
  double angle, newdiff;
  static int inited[3];
  static POINTS_ *pnow;
  static POINTS_ *pfrom;
  static int first;
  int i;
#define xangle(f,p) atan2((double)(f->y - p->y),(double)(f->x - p->x))
#define tangle(f,p) f->y==p->y && f->x==p->x
#define tprint(f,p,sf,sp) printf("%s=%d, %s=%d\n", sf, f-list, sp, p-list)

  if (k <= 0)		/* first or next */
    {
      first = k;
      pnow = q;
      pfrom = qln;
      for (i = 0; i < 3; i++)
	inited[i] = 0;
      return (pnow);
    }

  i = k - 1;
  if (qln == NULL)
    qln = q;			/* If don't have one yet, take it */
  else if (!first)		/* Else if already have one, and.. */
    {				/* point is not first on line.. */
      if (!inited[i])
	{
	  /*if (tangle(pnow, pfrom))
	    tprint(pnow,pfrom,"pnow","pfrom");
	  if (tangle(qln, pnow))
	    tprint(qln, pnow, "qln","pnow");*/
	  angle1[i] = xangle(pnow, pfrom);	/* atan returns -pi to pi */
	  angle2[i] = xangle(qln, pnow);
	  diff[i] = anglecmp(angle1[i], angle2[i]);
	  inited[i]++;
	}
      angle = xangle(q, pnow);
      /*if (tangle(q,pnow))
	tprint(q,pnow,"q","pnow");*/
      newdiff = anglecmp(angle1[i], angle);
      if (newdiff < diff[i])
	{
	  diff[i] = newdiff;
	  angle2[i] = angle;
	  qln = q;
	}
    }
  return (qln);
}

/*-----------------------------------------------------------------------------
|	anglecmp
-----------------------------------------------------------------------------*/
double anglecmp(double a1, double a2)
{
  double d;
  d = fabs(a2 - a1);		/* value 0 to 2*pi */
  if (d > pi)
    d = twopi - d;
  return (d);
}

/*-----------------------------------------------------------------------------
|    surround -- load array 'ilist' with indices of points surround 'ind'
-----------------------------------------------------------------------------*/
FLAG *surround(FLAG * ilist, POINTS * p)
{
  FLAG *ls;
  int iz,ir, i,j, i1,i2, j1,j2, ii,jj, ni;
  int rpolar, apolar;
  byte *q;
  static int inited = 0;
  static byte at[9];
  static byte right[6];
  static byte above[6];
#define X    PS_AT
#define R    PS_RIGHT
#define A    PS_ABOVE

  if (!inited)				/* load these arrays once only */
    {
      q = at;
      *q++ = X | R | A;
      *q++ = X | R;
      *q++ = X | A;
      *q++ = X | A;
      *q++ = 0;
      *q++ = X | A;
      *q++ = X | R;
      *q++ = X | R;
      *q++ = X;

      q = right;
      *q++ = X | R | A;		/* 0,-1 */
      *q++ = X | A;		/* 1,-1 */
      *q++ = A;			/* 0,0 */
      *q++ = A;			/* 1,0 */
      *q++ = X | R;		/* 0,1 */
      *q++ = X;			/* 1,1 */

      q = above;
      *q++ = X | R | A;		/* -1,0 */
      *q++ = R;			/*  0,0 */
      *q++ = X | A;		/*  1,0 */
      *q++ = X | R;		/* -1,1 */
      *q++ = R;			/*  0,1 */
      *q++ = X;			/*  1,1 */
      inited++;
    }

  rpolar = apolar = 0;
  if (p==NULL) printf("NULL pointer\n");
  if (p->where & PS_AT)
    {
      q = at;
      i1 = -1;
      j1 = -1;
    }
  else if (p->where & PS_RIGHT)
    {
      q = right;
      i1 = 0;
      j1 = -1;
      rpolar=polar;
    }
  else if (p->where & PS_ABOVE)
    {
      q = above;
      i1 = -1;
      j1 = 0;
      apolar=polar;
    }

  iz = p->iy;
  ir = p->ix;
  i2 = j2 = 1;
  ni = i2 - i1 + 1;
  for (j=j1, ls=ilist; j<=j2; j++)
    {
      for (i=i1; i<=i2; i++, q++)
	{
	  ii = i;
	  jj = j;
	  if (rpolar)
	    {
	      if (iz==0) continue;		/* no cases of this */
	      if (ir==mr-1 && i==1) ii = -ir + 0;
	      if (ii!=i) goto ADDIT;
	    }
	  else if (apolar)
	    {
	      if (ir==0 && i<0)	    ii = -ir + mr-1;
	      if (ir==mr-1 && i==1) ii = -ir + 0;
	      if (ii!=i) goto ADDIT;
	    }

	  if (ir == 0 && i < 0)
	    continue;
	  if (iz == 0 && j < 0)
	    continue;
	  if (ir == (mr-1) && (i==1 || (i==0 && *q==R)))
	    continue;
	  if (iz == (mz-1) && (j==1 || (j==0 && *q==A)))
	    continue;

	  if (!*q) continue;		/* skip at[0,0] */

	ADDIT:
	  ls->ix = (byte) (ir + ii);
	  ls->iy = (byte) (iz + jj);
	  ls->where = *q;
	  ls++;
	}
    }
  return (ls);
}

/*-----------------------------------------------------------------------------
|	has_psi
|	* see if psi occurs in ir,iz...ir+1,iz or ir,iz...ir,iz+1
|	* if yes, load into q, increment q
|	* save_linear, save_extremum generate quad coefficients and
|	  increment q as a counter only!
-----------------------------------------------------------------------------*/
POINTS_ *has_psi(float psi, int ir, int iz, POINTS_ * q)
{
  float delt, delta, psi0, xa,ya,ra,ta;
  int irz, verysmall;
  QCOEFF kx, ky;
  POINTS_ *q0, *q2, *p, *psame;
  BIN_ *b;

  irz = getp00(ir,iz);		/* load p00,p01,p10 etc */
  b = bin + irz;
  kx.flag = ky.flag = 0;
  q0 = q;
  psi0 = *p00;
  delta = psi - psi0;
  if (*p00 == psi)		/* ......... *p00 == psi */
    {
      q->ix = (byte) ir;
      q->iy = (byte) iz;
      q->where = PS_AT;
      q->used = 0;
      q->x = (int) ((float) ir * rscale + roffset);
      q->y = (int) ((float) iz * zscale + zoffset);
      q++;
      return (q);
    }

  if (*p00 < psi)				/* ......... *p00 < psi */
    {
      if (polar && iz==0) ;
      else if (polar || ir<(mr-1))		/* r: look for *p10 > psi */
	{
	  if (*p10 == psi) ;
	  else if (*p10 > psi)
	    q = save_linear(ir, mr, 1, &kx, q);
	  else if (ir == 0 || ir == mr - 2);
	  else if ((b - 1)->rsign <= 0 || (b + 1)->rsign >= 0);
	  else
	    q = save_extremum(ir, mr, 1, psi, 1, &kx, q);
	}
      if (iz < (mz - 1))			/* z: look for *p01 > psi */
	{
	  if (*p01 == psi );
	  else if (*p01 > psi)
	    q = save_linear(iz, mz, mr, &ky, q);
	  else if (iz == 0 || iz == mz - 2);
	  else if ((b - mr)->asign <= 0 || (b + mr)->asign >= 0);
	  else
	    q = save_extremum(iz, mz, mr, psi, 1, &ky, q);
	}
    }
  else						/* ......... *p00 > psi */
    {
      if (polar && iz==0) ;
      else if (polar || ir<(mr-1))		/* r: look for *p10 < psi */
	{
	  if (*p10 == psi) ;
	  else if (*p10 < psi)
	    q = save_linear(ir, mr, 1, &kx, q);
	  else if (ir == 0 || ir == mr - 2);
	  else if ((b - 1)->rsign >= 0 || (b + 1)->rsign <= 0);
	  else
	    q = save_extremum(ir, mr, 1, psi, -1, &kx, q);
	}
      if (iz < (mz - 1))			/* z: look for *p01 < psi */
	{
	  if (*p01 == psi) ;
	  else if (*p01 < psi)
	    q = save_linear(iz, mz, mr, &ky, q);
	  else if (iz == 0 || iz == mz - 2);
	  else if ((b - mr)->asign >= 0 || (b + mr)->asign <= 0);
	  else
	    q = save_extremum(ir, mr, mr, psi, -1, &ky, q);
	}
    }

  q2 = q;			/* q is a count of how many to add */
  for (q=p=q0; q<q2; q++)	/* We have quad coeffs, now add to list */
    {				/* calculate delt, use ir+delt or iz+delt */
      p->ix = (byte) ir;
      p->iy = (byte) iz;
      p->used = 0;
      if (q->where == PS_RIGHT)
	{
	  delt = (*p10 == psi0) ? FZ : delta / (*p10 - psi0);
	  if (kx.flag) delt = quadratic(&kx, psi, delt);
	}
      else
	{
	  delt = (*p01 == psi0) ? FZ : delta / (*p01 - psi0);
	  if (ky.flag) delt = quadratic(&ky, psi, delt);
	}

      /*if (fabs(delt)<.000001 || fabs(delt-(float)1.)<.000001)
	printf("%d,%d:%d has delt = %f\n", ir,iz,q->where,delt);*/
      verysmall = (delt <= .00001) || ((1.-delt) < .00001);
      p->where = q->where;
      if (!polar)
        {
	  xa = (float)ir;
	  ya = (float)iz;
	  if (p->where == PS_RIGHT) xa += delt;
	  else if (p->where != PS_AT) ya += delt;
	  p->x = (int)(xa * rscale + roffset);
	  p->y = (int)(ya * zscale + zoffset);
	  if (verysmall && (psame = same_as(p)) )
	    { psame->where = PS_AT; p--; }
	}
      else
	{
	  ta = (float)ir;
	  ra = get_loopval(rlp, iz);
	  if (p->where == PS_RIGHT) ta += delt;
	  else if (p->where != PS_AT)
	    ra += delt * (get_loopval(rlp, iz+1) - ra);
	  ta *= dtheta;
	  xa = ra * cos(ta);
	  ya = ra * sin(ta);
	  p->x = (int)(xa * xscale + xoffset);
	  p->y = (int)(ya * yscale + yoffset);
        }
      p++;
    }

  return (p);
}

/*-----------------------------------------------------------------------------
|	same_as
-----------------------------------------------------------------------------*/
POINTS_ *same_as(POINTS_ *p)
{
  POINTS_ *q;
  int ix, iy;
  ix = p->ix;
  iy = p->iy;
  for(q=p-1; q>=list && q->iy>=iy-1; q--)
    {
      if (q->ix >= ix-1 && q->ix <= ix+1 && (q->x==p->x && q->y==p->y))
	return q;
    }
  return 0;
}

/*-----------------------------------------------------------------------------
|	getp00
-----------------------------------------------------------------------------*/
int getp00(int ir,int iz)
{
  int irz;
  long sign;
  irz = mr * iz + ir;
  
  p00 = buf + off0 + ir * xstep + iz * ystep;
  pm0 = p00 - xstep;
  p10 = p00 + xstep;
  p20 = p10 + xstep;

  p0m = p00 - ystep;
  p01 = p00 + ystep;
  p02 = p01 + ystep;
  if (!polar) return(irz);
			  /* "iz"=r, "ir"=theta */
  sign = (ir < mr/2) ? xstep : -xstep;
  if (iz==0)
    p0m = p01 + sign * mr/2;

  if (ir==0) pm0 = p00 + xstep * (mr -1);
  if (ir==mr-1) p10 = p00 - xstep * (mr-1);
  p20 = (ir!=mr-2) ? p10 + xstep : p00 - xstep * (mr-2);
  
  return(irz);
}

/*-----------------------------------------------------------------------------
|	save_linear
-----------------------------------------------------------------------------*/
POINTS_ *save_linear(int ir, int mr, int n, QCOEFF *kp, POINTS_ * q)
{
  int right;
  right = (n==1);
  if (quad_always && (ir > 0 && ir < mr - 2))
    {
      getquadcoeff(n, ir, mr, kp);
      kp->flag++;
    }
  else
    {
      kp->a = (float) 0;
      kp->b = right ? *p10 - *p00 : *p01 - *p00;
      kp->c = *p00;
    }
  q->where = right ? (byte) PS_RIGHT : (byte) PS_ABOVE;
  return (++q);
}

/*-----------------------------------------------------------------------------
|	save_extremum
-----------------------------------------------------------------------------*/
POINTS_ *save_extremum(int ir, int mr, int n, float psi,
		       int how, QCOEFF *k, POINTS_ * q)
{
  int ok;
/*if (ir==0 || ir==mr-2) return(q);
 *if (b1sign <= 0 || b2sign >= 0) return(q);
 */

  getquadcoeff(n, ir, mr, k);
  k->g = k->c - k->b * k->b / (4. * k->a);
  k->flag++;
  if (how == 1)
    ok = (k->b > FZ && (float) 2 * k->a < -k->b && k->g >= psi);
  else
    ok = (k->b < FZ && (float) 2 * k->a > -k->b && k->g <= psi);
  if (ok)
    {
      q->where = (n == 1) ? (byte) PS_RIGHT : (byte) PS_ABOVE;
      q++;
    }
  return (q);
}

/*-----------------------------------------------------------------------------
|   getquadcoeff
-----------------------------------------------------------------------------*/
void getquadcoeff(int n, int ir, int mr, QCOEFF * k)
{
  float a1, a2, b1, b2, c1;
  static int cubic = 0;
  float af0,af1,af2,af3;
  
  af1 = *p00;
  if (n==1)
    {
      af2 = *p10;
      af3 = *p20;
      af0 = *pm0;
    }
  else
    {
      af2 = *p01;
      af3 = *p02;
      af0 = *p0m;
    }

  if (!cubic)
    goto QUAD;
  k->a = (af3 - 3 * af2 + 3 * af1 - af0) / 6.;
  k->b = (3 * af2 - 6 * af1 + 3 * af0) / 6.;
  k->c = (-af3 + 6 * af2 - 3 * af1 - 2 * af0) / 6.;
  return;

QUAD:
  c1 = (af1 + af2) / 2.;
  b1 = (af2 - af0) / 4.;
  b2 = (af3 - af1) / 4.;
  a1 = (af2 - 2 * af1 + af0) / 4.;
  a2 = (af3 - 2 * af2 + af1) / 4.;

  k->c = c1 - b2 + a2;
  k->b = b1 + b2 - 2. * a2;
  k->a = a1 + a2;
}

/*-----------------------------------------------------------------------------
|	quadratic
|	* find point in (0..1) where value = psi
|	* x0 = linear likely guess (psi-psi1)/(psi2-psi1)
-----------------------------------------------------------------------------*/
float quadratic(QCOEFF * k, float psi, float x0)
{
  float f, df;
  float x, xx, dxx, err;

  k->c -= psi;
  for (x = x0;;)
    {
      f = k->a * x * x + k->b * x + k->c;
      if (f == FZ) return(x);
      df = 2 * k->a * x + k->b;
      xx = x;
      dxx = f / df;
      x = xx - dxx;
      err = (x - xx) / (x + xx);
      if (err < FZ)
	err = -err;
      if (err < (float) .0001)
	return (x);
    }
}

#ifdef DEAD_CODE
/*-----------------------------------------------------------------------------
|	get_cont_minmax
-----------------------------------------------------------------------------*/
void get_cont_minmax(int how, int i, float *xn, float *xx,
				     float *yn, float *yx)
{
  float rfac, zfac;
  struct CONT *p;
  p = cdata + i;
  
  if (how == 1 && !polar)
    {
      *xn = *yx = (float) 0;
      *xx = (float) (p->mr - 1);
      *yn = (float) (p->mz - 1);
    }
  else if (how==1 && polar)
    {
      *xn = *yx = -(float)(p->mz-1);
      *xx = *yn =  (float)(p->mz-1);
    }
  else
    {
      rfac = (p->rmax - p->rmin) / (float) (p->mr - 1);
      zfac = (p->zmax - p->zmin) / (float) (p->mz - 1);
      *xn = *xn * rfac + p->rmin;
      *xx = *xx * rfac + p->rmin;
      *yn = *yn * zfac + p->zmin;
      *yx = *yx * zfac + p->zmin;
    }
}
#endif

/*-----------------------------------------------------------------------------
|	new_ncurve
-----------------------------------------------------------------------------*/
int new_ncurve(CURVE_SET *cp, char how)
{
  if (cp->gtype == 'G') return(0);
  if (how >= 'a' && how <= 'z') how += ('A'-'a');
  if (how == 'D') cp->ncurve *= 2;
  else		  cp->ncurve /= 2;
  if (cp->ncurve < 1) cp->ncurve = 1;
  redrawflag = 1;
  return (1);
}
/*-----------------------------------------------------------------------------
|	contour_values
-----------------------------------------------------------------------------*/
void contour_values(CURVE_SET *cp)
{
   if (cp->gtype != 'G')
      show_values = redrawflag = 1;
}

/*=============================================================================
**                  DEBUG FUNCTIONS
**===========================================================================*/
/*-----------------------------------------------------------------------------
|   savegrid
|       write grid.dat = contains values on the grid.  if debug2 set.
-----------------------------------------------------------------------------*/
void savegrid(int ir, int iz)
{
  static FILE *f2;
  float xpsi;
  int ipsi,i,jz;
  static char text[] = "     ";
  static char format[] = "%5d";	/* # chars in max, below */
  static float mul = (float) 1000.;
  static float max = (float) 10000;
  static int irmin=0, irmax = 21;
  float pmax;

  if (ir == mr - 1 && iz == mz - 1)
    {
      fclose(f2);
      debug2 = 0;
      return;
    }

  if (ir == 0 && iz == 0)
    {
      f2 = fopen("grid.dat", "wt");
      pmax = current_cp->iy.max;			/* (assumes values are all positive) */
      if (pmax > (float)999) mul=(float)1;
      else if (pmax > (float)99) mul=(float)10;
      else if (pmax > (float)9)  mul=(float)100;
      else mul=(float)1000;
    }
  if (ir == 0)
    fprintf(f2, "\n");
  if (ir >= irmin && ir < irmax)	/* (Can't fit ALL on the page) */
    {
      jz = mz-1-iz;
      if (!polar) i = jz * mr + ir;
      else if (jz==0) i=0;
      else i = 1 + (jz-1)*mr + ir;
      xpsi = *(buf + i) * mul;		/* mul ==> # decimal places keep */
      if (xpsi > FZ)
	xpsi += (float) .5;
      else
	xpsi -= (float) .5;
      if (xpsi >= max || xpsi <= -max)
	fprintf(f2, text);
      else
	{
	  ipsi = (int) (xpsi);
	  fprintf(f2, format, ipsi);
	}
    }
}

/*-----------------------------------------------------------------------------
|	savelevel -- from debug1, writes to disk all points at current level
-----------------------------------------------------------------------------*/
void savelevel(int how, int index, char *format)
{
#ifdef UNIX
  POINTS *p, *p2;
#else
  POINTS huge *p, huge * p2;
#endif
  static FILE *file;
  char text[5];

  if (how == 1)
    {
      p2 = list + index;
      file = fopen("debug.dat", "wt");
      fprintf(file, "Summary of relevant vertices A,R edges\n");
      fprintf(file, "Format is index: ix,iy,<where> ==> x,y\n");
      for (p = list; p < p2; p++)
	{
	  strcpy(text, "    ");
	  if (p->where & 1)  *(text) = 'V';
	  if (p->where & 2)  *(text + 1) = 'R';
	  if (p->where & 4)  *(text + 2) = 'A';
	  if (p->where & 0xf8) *(text + 3) = 'O';
	  fprintf(file, "%3d: %3d,%3d<%s> ==> %4d,%4d\n",
		  p - list, p->ix, p->iy, text, (int) p->x, (int) p->y);
	}
    }
  else if (how == 0)
    {
      p2 = list + index;
      fprintf(file, "Done\n");
      for (p = list; p < p2; p++)
	{
	  strcpy(text, "  ");
	  if (p->used & PS_FROM) *text = 'F';
	  if (p->used & PS_TO)   *(text + 1) = 'T';
	  fprintf(file, "%d:<%s> \n", p - list, text);
	}
      fclose(file);
    }

  else
    {
      fprintf(file, format, index);
      fclose(file);
      file = fopen("debug.dat", "at");
    }
}

#ifdef DEAD_CODE
/*-----------------------------------------------------------------------------
|   psicount
-----------------------------------------------------------------------------*/
void psicount(float *buf, long bufsize, float psi0, float dpsi,
	      int ncurve, float *count)
{
  int i;
  float *p, *buf2;
  long nvalues;

  for (i = 0; i <= ncurve; i++)
    count[i] = FZ;

  nvalues = bufsize / 4L;
  buf2 = buf + nvalues;

  for (p = buf; p < buf2; p++)
    {
      i = (int) ((*p - psi0) / dpsi);
      count[i]++;
    }
}
#endif
