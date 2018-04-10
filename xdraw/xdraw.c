/******************************************************************************
**  NAME      XDRAW.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Usage: draw [options] [suffix]
**      options: -d = draw (draw.in), -c = cntour (cont.in)
**
**  Copyright (c) GlassWare 1992.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <time.h>

#ifdef UNIX
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

#else
#include <dos.h>
#include "xlib.h"
#include "xutil.h"
#include "xatom.h"
#endif

#define XDRAW_SOURCE
#include "gendefs.h"
#include "curves.h"
#include "xdraw.h"
#include "xinit.h"
#include "xcontour.h"
#include "xtools.h"
#include "setcolor.h"
#include "menu.h"
#include "ps.h"

static void draw_ten(int, int, int, int, float);
static void drawticks(int is_y,
		float xlim[], float xmax, float ylim[], float ymin,
		float xscale, float xoffset, float yscale, float yoffset);
static char *ticks(float *, float *, float *, int *, float *);

float *buf;
int xterm;
int default_ncol = 3, figures = 4;
int default_ncurve = 16;
int exitflag = 0;
int redrawflag;

Window dialog_win;
GC dialog_gc;
#ifdef MSWINDOWS
#include <windows.h>
#include "wdraw.h"
#define expose_event() winDummy()
#define process_event() winEvent()
int dialogwindow = 1;		/* Windows MUST have its own window */
#else
#define expose_event() event()
#define process_event() event()
int dialogwindow = 0;		/* XTerm can use normal dialog window */
#endif

char title[1000];
int titledrawn=0, showcomment=0;

LOOP loop[20];
int nloop=0;

NODE *nodelist;
int nnode=0;			/* TOTAL count of nodelist */
int inode = -1, ivar = -1;
int ncount_equal = 0;

CURVE_SET curveset[NCS];	/* Describes the plots */
CURVE_SET *cp, *cp2;
int ncset=0;

LOOP loop[20];
int nloop, iloop;

int ftype;			/* 0=graph, 1=contour, 2=self-determined */
extern int ps_modeon;
extern unsigned int xscreen, yscreen;

#ifdef UNIX
#define incre(p) p++
#define FP_inc(p)    p++;
#define FP_incn(p,n) p+=n;
#define millisec(c) c /= (CLOCKS_PER_SEC / 1000)

#else
#define incre(p) if (FP_OFF(p)>=0xfffc) \
		{ FP_SEG(p)+=FP_OFF(p)/16;FP_OFF(p)&=0xf; } p++
#define FP_adjust(p) {if (FP_OFF(p) & 0xff00) \
		  { FP_SEG(p) += FP_OFF(p)/16; FP_OFF(p) &= 0xf; }}
#define FP_inc(p)    {FP_adjust(p); p++;}
#define FP_incn(p,n) {FP_adjust(p); p+=n;}
#define millisec(c) c = (c * 1000) / CLK_TCK
#endif

/*-----------------------------------------------------------------------------
|    main
-----------------------------------------------------------------------------*/
void main(argc, argv)
int argc;
char *argv[];
{
  CURVE_SET *cp;
  int nwin, ncol, nwinplus;
  char text[200];

  init(argc, argv);				/* read input data */
  for(cp=curveset,nwin=0; cp<cp2; cp++)
    if (cp->which=='1')
      {
	if (nwin<MAXWINDOW) nwin++;
	else { cp2 = cp; break; }
      }

  ncol = default_ncol;
  if (nwin >= 9 && dialogwindow) { dialogwindow = 2; nwinplus = nwin; }
  else nwinplus = nwin + (dialogwindow ? 1 : 0);
  if ((ncol==2 && nwinplus>4) || ncol>2) ncol = 3;
  else if (ncol==2 && nwinplus==1) ncol = 1;
  menu_setwindows(nwin);
  xterm = opendisplay(title, ncol);

  for (cp = curveset; cp < cp2; cp++)		/* create, draw each window */
    {
      if (cp->which=='1')
        {
	  parse_title(cp, cp->title, text);
	  cp->window = makewindow(cp->title, 0);
	}
      else cp->window = (cp-1)->window;
      if (cp<cp2-1 && (cp+1)->which=='2') continue;

      if (xterm) expose_event();		/* Draw (Expose event) */
      else postscript(cp->title);		/* no xterm: write ps */
    }

  if (dialogwindow)				/* UNIX or Windows */
    dialog_win =  makewindow("XDraw Dialogue", 1);
  else init_menus();				/* Motif or Sheryl's menu */
  if (!ncset) ;					/* ncset=#curves=set in xinit.c */
  else while (!exitflag) process_event();	/* process events, Q= exit */
  if (xterm) closedisplay();
}

/*=============================================================================
**		ARRAY INDEXING UTILITIES
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	get_loopval
-----------------------------------------------------------------------------*/
float get_loopval(LOOP *qs, int i)
{
  float *p;
  p = qs->extra;
  if (qs->ltype=='I') return((float)i);
  else if (qs->ltype=='X') return(*p + i * *(p+2));
  else if (qs->ltype=='A') return(*(p+i));
  else return((float)0);
}

/*=============================================================================
**                    SUPPORTING FUNCTIONS
**===========================================================================*/
Display *mydisplay;		/* set via redraw()-->get_expose_info() */
Window redraw_win;
int myscreen;
GC redraw_gc;
XRectangle clipbox, bigbox;
float clipx1, clipx2, clipy1, clipy2;	/* for xdraw's own clip of BIG numbers*/
unsigned int width, height;
int font_height, monochrome=0;
int wx0, wy0, w, h, wx2, wy2;
int npair, pairs[20], preserve_aspect;
float xscale, xoffset, yscale, yoffset;
static char ps_msg[] = "Generating PostScript file...";

/*-----------------------------------------------------------------------------
|	get_box_info
-----------------------------------------------------------------------------*/
int get_box_info(Window win, XRectangle *bbx, XRectangle *cbx,
	         float *xminp, float *yminp, float *xmaxp, float *ymaxp)
{
  int icurve, i, i0;
  int wxb, wya, wyb, border;
  unsigned int depth;
  float xmin,xmax,ymin,ymax, dp;
  float xmin0, xmax0, ymin0, ymax0;
  CURVE_SET *cp, *cq;
  int zoomed;
  static int pass=0;

  if (!ps_modeon) get_properties(win, &width, &height, &depth);
  else { width = 1264; height=984; depth=4; }		/* condor full screen */
  monochrome = (depth==1);
  /*printf("bits per pixel = %d\n", depth);*/

  icurve = i0 = xterm ? get_graph(win, 0) : get_graph((Window)0, pass++);
  if (icurve==-1) return(-1);
  cp = curveset + icurve;

  border = (int) (.02 * ((width >= height) ? width : height));
  wx0 = 4 * font_height;
  wxb = width - border;
  wya = border + 2 * font_height;
  wyb = height - border;
  if (cp->subtitle && !ps_modeon) wyb -= font_height;
#ifndef UNIX
  wyb -= font_height;
#endif
  wy0 = height - wyb;
  w = wxb - wx0;
  h = wyb - wya;
  wx2 = wxb;
  wy2 = wyb;

  bigbox.x = bigbox.y = 0;
  bigbox.width = width;
  bigbox.height = height;

  clipbox.x = wx0;	clipx1 = (float)clipbox.x;
  clipbox.y = wy0;	clipy1 = (float)clipbox.y;
  clipbox.width = w;	clipx2 = clipx1 + (float)clipbox.width;
  clipbox.height = h;	clipy2 = clipy1 + (float)clipbox.height;
  givebox(wx0, wy0, w, h, width, height);

/*------ Get min, max for all curves in redraw_win */
  if (cp->gtype == 'B') return(-1);

  for(cq=cp,i=icurve; cq<cp2 && (i==i0 || cq->which=='2'); i++,cq++)
    {
      pairs[i-i0] = i;
      if (i==i0 || cq->ix.min < xmin) xmin = cq->ix.min;
      if (i==i0 || cq->ix.max > xmax) xmax = cq->ix.max;
      if (i==i0 || cq->iy.min < ymax) ymax = cq->iy.min;
      if (i==i0 || cq->iy.max > ymin) ymin = cq->iy.max;
    }
  npair = i-i0;
  if (cp->flags & POLAR)
    {
      xmax = ymin =  cp->iy.max;
      xmin = ymax = -xmax;
    }

  if (cp->gtype == 'G')					/* set margins */
    {
      if ((cp->forcemin != FZ || cp->forcemax != FZ) &&
	   cp->forcemin != cp->forcemax && !(cp->flags & E_IGNORE))
        {
	  ymax = cp->forcemin;
	  if (cp->forcemax > cp->forcemin)
	    ymin = cp->forcemax;
        }
      if ((cp->flags & E_VIS) && (cp->ncurve != -1 || cp->hzcurve != -1))
	{
          for(i=0; i<npair; i++)
	    redraw00(icurve, i, npair, &xmin, &xmax, &ymin, &ymax);
	}
      if (xmax!=xmin) dp = (float) .05 *(xmax - xmin);
      else if (xmin)  dp = xmin * FH;
      else dp = FH;
      xmin -= dp;
      xmax += dp;

      if (ymax!=ymin) dp = (float) .05 *(ymin - ymax);
      else if (ymin)  dp = ymin * FH;			/* ymin > ymax */
      else dp = FH;
      ymin += dp;
      ymax -= dp;
    }

/*------ Get new xmin..ymax if zoomed or preserve aspect */
  preserve_aspect = cp->flags & ASPECT;
  xmin0 = xmin; xmax0 = xmax;
  ymin0 = ymin; ymax0 = ymax;
  zoomed = get_zoomed_limits(cp, &xmin, &ymin, &xmax, &ymax) ;
  if (preserve_aspect) get_aspected_limits(win, w, h,
      &xmin, &ymin, &xmax, &ymax, xmin0, ymin0, xmax0, ymax0);
  *xminp = xmin; *xmaxp = xmax;
  *yminp = ymin; *ymaxp = ymax;
  if (bbx != &bigbox)  *bbx = bigbox;
  if (cbx != &clipbox) *cbx = clipbox;
  return(icurve);
}

/*-----------------------------------------------------------------------------
|	get_aspected_limits
|	* Only called when zoom flag NOT set
|	* Calculates aspected limits, sets new zoom borders
-----------------------------------------------------------------------------*/
void get_aspected_limits(Window win, int w, int h,
			 float *xmin, float *ymin, float *xmax, float *ymax,
			 float xmin1, float ymin1, float xmax1, float ymax1)
{
  float xc, yc, dx, dy;
  float aspect_d, aspect_w;

  dx = *xmax - *xmin;
  dy = *ymin - *ymax;
  aspect_w = (float)h / (float)w;		/* aspect of box */
  aspect_d = dy / dx;				/* aspect of data */
  if (aspect_d > aspect_w)			/* data is tall */
    {
      xc = (*xmax + *xmin) / 2;
      dx = dy / aspect_w;
      *xmin = xc - dx/2;
      *xmax = xc + dx/2;
    }
  if (aspect_d < aspect_w)			/* data is wide */
    {
      yc = (*ymax + *ymin) / 2;
      dy = dx * aspect_w;
      *ymax = yc - dy/2;
      *ymin = yc + dy/2;
    }
  set_aspected_limits(win, 0, xmin1, xmax1, *xmin, *xmax);
  set_aspected_limits(win, 1, ymin1, ymax1, *ymin, *ymax);
}

/*-----------------------------------------------------------------------------
|	get_zoomed_limits
|	* IN:  data's autoscaled unzoomed pointers to xmin..ymax
|	* OUT: adjusted by fractional zoom window fx1..fy2
-----------------------------------------------------------------------------*/
int get_zoomed_limits(CURVE_SET *cp, float *xmin, float *ymin, float *xmax, float *ymax)
{
  float fx1, fy1, fx2, fy2, dx, dy;
  int retval=0;
  if (get_zoomedclip(cp->window, &fx1,&fy1,&fx2,&fy2))
    {
      dx = *xmax - *xmin;
      dy = *ymin - *ymax;
      *xmax = *xmin + fx2 * dx;      *xmin = *xmin + fx1 * dx;
      *ymax = *ymin - fy1 * dy;      *ymin = *ymin - fy2 * dy;
      retval = 1;
    }
  return retval;
}

/*-----------------------------------------------------------------------------
|	draw_frame
-----------------------------------------------------------------------------*/
void draw_frame()
{
  int ixoff, iyoff;
  float ixmin, ixmax, iymin, iymax;
  
  if (!ps_modeon)
    {
      clear_if_needed( redraw_win, redraw_gc, bigbox.width, bigbox.height);
      XSetClipRectangles(mydisplay, redraw_gc, 0, 0, &bigbox, 1, Unsorted);
      XDrawRectangle(mydisplay, redraw_win, redraw_gc, wx0, wy0, w, h);
      XSetClipRectangles(mydisplay, redraw_gc, 0, 0, &clipbox, 1, Unsorted);
    }
  else if (ps_modeon)
    {
      xmessage(wx0 + font_height, wy0 + font_height, ps_msg);
      ps_adjustpage(height, wx0, wy0, w, h);
      ps_color(FZ, FZ, FZ);
      ps_clip(0, 0, 0, 0);
      ps_rectangle(wx0, wy0, w, h);
      ps_clip(clipbox.x, clipbox.y, clipbox.width, clipbox.height);
    }

/*------ draw vertical line */
  if (cp->gtype != 'G') goto ENDLINE;
  ixmin = (clipx1 - xoffset) / xscale;
  ixmax = (clipx2 - xoffset) / xscale;
  iymin = (clipy2 - yoffset) / yscale;
  iymax = (clipy1 - yoffset) / yscale;

  ixoff = (int) xoffset;
  iyoff = (int) yoffset;
  if ((float)0 >= ixmin && (float)0 <= ixmax)
    {
      if (!ps_modeon)
	XDrawLine(mydisplay, redraw_win, redraw_gc, ixoff, wy0,
		  ixoff, wy0 + h);
      else
	ps_line(ixoff, wy0, ixoff, wy0 + h);
    }

/*------ draw horizontal line */
  if ((float)0 >= iymin && (float)0 <= iymax)
    {
      if (!ps_modeon)
	XDrawLine(mydisplay, redraw_win, redraw_gc, wx0, iyoff,
		  wx0 + w, iyoff);
      else
	ps_line(wx0, iyoff, wx0 + w, iyoff);
    }

ENDLINE:
  if (ps_modeon)
    ps_thinline();
}

/*-----------------------------------------------------------------------------
|	redraw. redraws all contents of window on expose events.
-----------------------------------------------------------------------------*/
void redraw()
{
  int i;
  float xmin, xmax, ymin, ymax, xlim2[2], ylim2[2];
  char *p;

/*------ get properties of window */
  get_expose_info(&mydisplay, &redraw_win, &redraw_gc, &font_height);
  i = get_box_info(redraw_win, &bigbox, &clipbox, &xmin, &ymin, &xmax, &ymax);
  if (i==-1) return;
#ifndef UNIX
  if (bigbox.width < xscreen*3/4 && !titledrawn && showcomment)
    {
      for(p = title; *p; p += strlen(p) + 1)
	  xPrintf1(p, (p==title)?1:0);
      titledrawn = 1;
    }
#endif

  cp = curveset + i;

  /*------ Get scale factors and offsets */
  getscale(wx0, w, xmin, xmax, &xscale, &xoffset);
  getscale(wy0, h, ymin, ymax, &yscale, &yoffset);

  xoffset += FH;	/* roundoff */
  yoffset += FH;

  draw_frame();		/* draw frame, horizontal & vertical lines */

  if (cp->gtype=='G')				/*.....Draw Graph */
    {
      for(i=0; i<npair; i++)
	redraw0(pairs[i], i, npair, xmin, xmax, ymin, ymax);
    }

  else if (cp->gtype == 'P')			/*.....Test Palette (PC) */
    testpalette(wx0, wy0, w, h);

  else						/*.....Contour */
    redraw1(cp, xmin, xmax, ymin, ymax);

  if (!ps_modeon)
    XSetClipRectangles(mydisplay, redraw_gc, 0, 0, &bigbox, 1, Unsorted);
  else
    {
      ps_color(FZ, FZ, FZ);
      ps_clip(0, 0, 0, 0);
    }

/*------ draw tick marks on axes */
  xlim2[0] = xmin;
  xlim2[1] = xmax;
  ylim2[0] = ymax;
  ylim2[1] = ymin;

  label(xlim2, xscale, xoffset,	ylim2, yscale, yoffset, cp);

  if (ps_modeon)
    xmessage(wx0 + font_height, wy0 + font_height, ps_msg);
}

/*-----------------------------------------------------------------------------
|    getscale
-----------------------------------------------------------------------------*/
void getscale(int wx0, int dwx, float xmin, float xmax,
	      float *xscale, float *xoffset)
{
  *xscale = dwx / (xmax - xmin);
  *xoffset = wx0 - *xscale * xmin;
}

#ifdef DEAD_CODE
void getminmax(int i, int ns, float x, float y)
{
  static float xn,yn,xx,yx;
  static int inited=0;
  if (inited) return;
  if (i==0) { xn = xx = x; yn = yx = y; }
  else
    {
      if (x < xn) xn = x;  if (x > xx) xx = x;
      if (y < yn) yn = y;  if (y > yx) yx = y;
    }
  if (i == ns-1) { printf("xmin..ymax = %f,%f ... %f,%f\n",
			xn, yn, xx, yx); inited++; }
}
#endif

/*-----------------------------------------------------------------------------
|    get_graph
-----------------------------------------------------------------------------*/
int get_graph(Window mywindow, int i0)
{
  CURVE_SET *cp;
  int i, n;
  if (mywindow)
    for (i=i0, cp=curveset+i; cp < cp2 && cp->window != mywindow; i++, cp++);
  else
    {
      for (i=0,n=-1,cp=curveset; cp < cp2; i++, cp++)
	{
	  if (cp->which=='1') n++;
	  if (n==i0) break;
	}
      i = n;
    }
  return ((cp == cp2) ? -1 : i);
}

/*-----------------------------------------------------------------------------
|	toggle_markers
-----------------------------------------------------------------------------*/
void toggle_markers(CURVE_SET *cp1, int *mark)
{
  CURVE_SET *cp;
  if (cp1->gtype != 'G') return;

  for(cp=cp1; cp<cp2; cp++)
    {
      if (cp > cp1 && cp->which=='1') break;
      if (cp->flags & MARKERS)
	{ cp->flags &= (~MARKERS); *mark = 0; }
      else
	{ cp->flags |= MARKERS; *mark = 1; }
    }
  redrawflag = 1;
}

/*-----------------------------------------------------------------------------
|	toggle_aspect
-----------------------------------------------------------------------------*/
void toggle_aspect(CURVE_SET *cp1, int *aspect)
{
  CURVE_SET *cp;
  /*if (cp1->gtype != 'G') return;*/

  redrawflag=0;
  for(cp=cp1; cp<cp2; cp++)
    {
      if (cp > cp1 && cp->which=='1') break;
  /*----- only toggle unzoomed window ---------*/
      if (get_zoomedclip(cp->window, NULL,NULL,NULL,NULL) & 1)
	continue;
      if (cp->flags & ASPECT)
  /*----- aspect off: need set v->clipped to 0 ------*/
	{ cp->flags &= (~ASPECT); *aspect = 0; newclip(-1,0,0); redrawflag=1; }
      else
	{ cp->flags |= ASPECT; *aspect = 1; redrawflag=1;}
    }
}

/*-----------------------------------------------------------------------------
|	toggle_extrema_use
-----------------------------------------------------------------------------*/
void toggle_extrema_use(CURVE_SET *cp1)
{
  CURVE_SET *cp;
  int single, force, eflag;
  if (cp1->gtype != 'G') return;

  redrawflag=0;
  for(cp=cp1; cp<cp2; cp++)
    {
      if (cp > cp1 && cp->which=='1') break;
  /*----- only toggle if unzoomed ---------*/
      if (get_zoomedclip(cp->window, NULL,NULL,NULL,NULL) & 1)
	continue;
      single = cp->flags & SINGLE;
      force = (cp->forcemin != cp->forcemax);
      if (single || force)
	{
	  redrawflag = 1;
	  eflag = cp->flags & (E_VIS | E_IGNORE);
	  if (single && force)
	    {
	      eflag += E_VIS;
	      if (eflag > E_IGNORE) eflag = 0;
            }
	  else if (single) eflag = E_VIS - eflag;
          else if (force)  eflag = E_IGNORE - eflag;
	  cp->flags = (cp->flags & ~(E_VIS | E_IGNORE)) | eflag;
	}
    }
}

/*-----------------------------------------------------------------------------
|	toggle_fill
-----------------------------------------------------------------------------*/
void toggle_fill(CURVE_SET *cp1)
{
  CURVE_SET *cp;
  for (cp=cp1; cp<cp2; cp++)
    {
      if (cp > cp1 && cp->which=='1') break;
      if (!(cp->flags & DOTS))
	{
	  if (cp->flags & FILL) cp->flags &= ~FILL;
	  else cp->flags |= FILL;
	  redrawflag=1;
	}
    }
}

/*-----------------------------------------------------------------------------
|	toggle_single
-----------------------------------------------------------------------------*/
void toggle_single(CURVE_SET *cp1, char c)
{        
  int i, n, nt, *which;
  CURVE_SET *cp;
  char text[100];
  for(cp=cp1; cp<cp2; cp++)
    {
      if (cp > cp1 && cp->which=='1') break;
      if (c=='0') cp->flags &= ~SINGLE;
      else if (c=='1' || c=='2')
	{
	  if (!(cp->flags & SINGLE)) cp->ncurve=cp->hzcurve=-1;
	  cp->flags |= SINGLE;
	  which = (c=='1') ? &cp->ncurve : &cp->hzcurve;
	  (*which)++;
	  n = (npair>1) ? npair : (loop+cp->lfaml)->count;
	  if (*which == n) *which=0;
	  if (cp->ncurve >= 0 && cp->hzcurve >= 0 &&
	      cp->ncurve==cp->hzcurve) (*which)++;
	  if (*which >= n) *which=0;
	}
      else if ((c=='<' || c=='>') && nloop > 3)
	{
	  for(i=0; i<nloop && (n=cp->param[i])<0 ; i++) ;
	  if (i < nloop)
	    {
	      nt = loop[i].count;
	      if (c=='>') n++; else n--;
	      if (n < 0 ) n = nt - 1;
	      else if (n >= nt) n = 0;
	      cp->param[i] = n;
	      get_limits(cp);
	      if (parse_title(cp, cp->title, text))
		XChangeProperty(mydisplay, redraw_win, XA_WM_NAME, XA_WM_NAME,
			    8, PropModeReplace, (unsigned char *)text, 1);
	    }
	}
    }
  redrawflag=1;
}

/*-----------------------------------------------------------------------------
|	toggle_flabel
-----------------------------------------------------------------------------*/
void toggle_flabel(CURVE_SET *cp)
{
   if (cp->flags & LABELF) cp->flags &= ~LABELF;
   else cp->flags |= LABELF;
   redrawflag = 1;
}

/*-----------------------------------------------------------------------------
|	draw_graphtype
-----------------------------------------------------------------------------*/
void draw_graphtype(char *type)
{
  char *p;
  for (cp = curveset; cp < cp2; cp++)
    {
      for (p = type; *p && cp->gtype != *p; p++);
      if (!*p)
	continue;
      XClearArea(mydisplay, cp->window, 0, 0, 0, 0, True);
    }
}

/*=============================================================================
**                LABELS
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	ticks. computes positions of tick marks
|	* IN:  lim[] = min, max of data.
|	* OUT: new[] = integer start, end of ticks (may be outside)
|	*      majorz = separation, iexp=exponent
-----------------------------------------------------------------------------*/
static char *ticks(float *lim, float *new, float *dx, int *iexp, float *offset)
{
  int i, n, ii, ni, n1, n2;
  float ratio, min, max, amin, amax;
  float xmin, xmax, x;
  double t, tl, ti, diff, cutoff, expt, expfac;
  int pass;
#define NI 19		/* WARN! keep as sizeof(interval) */
  static double interval[] = {
    .001, .002, .005, .01,  .02,  .05, .1, .2, .5,
    1.,  2.,  5.,  10., 20., 50., 100., 200., 500., 1000.
    };
  static char *format[] = { "%.3lf", "%.2lf", "%.1lf", "%.0lf" } ;
#define neartick(min, i, n1p) modf((double)min/interval[i], &ti); *n1p = (int)ti

  *offset = FZ;
  if (figures<2) figures=2;
  if (figures>5) figures=5;
  cutoff = 2.5 * pow(10., -(double)figures);
  xmin = lim[0];
  xmax = lim[1];
  pass = 0;
START:  
  pass++;
  if (pass > 10) { printf("Crash in 'ticks'\n"); exit(0); }

  amin = (float) fabs(xmin);
  amax = (float) fabs(xmax);

  if (amin != FZ && amax != FZ)		/* replace small endpoints with 0 */
    {
      ratio = amax / amin;		/* e.g. 1..22, 2..44 etc. --> */
      if (ratio > (float) 22)		/*      0..22, 0..44 */
	xmin = FZ;
      if (ratio < (float) (1. / 22.))
	xmax = FZ;
    }

  t = (double)((amax > amin) ? amax: amin);	/* compute exponent */
  tl = log10(t);				/* e.g. 80-->1.9 ,.008-->-2.1*/
  expt = floor(tl);				/* e.g. 1.9-->1,  -2.1-->-3*/
  i = (int)expt;				/* exponent! */
  if (i > 0)
    {
      if (i <= 2) i = 0;	/* reduce exponent */
      /*else i -= 2;*/		/* (can give 6000 = 600 x 10^1) */
      expt = (double)i;
    }
  *iexp = i;
  expfac = pow(10., -expt);			/* e.g. 1-->.1,   -3-->1000*/

  diff = (float)expfac * (xmax - xmin);		/* see if need offset */
  if (diff < cutoff)
    {
   NEWOFF:	    
      tl = log10(diff);
      expt = floor(tl);
      expfac = pow(10., -expt);
      t = (double)(xmax+xmin) * expfac / 2.;
      for(;;)
	{
	  if (t>0. && t<1.) { t *= 10.; expfac *= 10.; }
	  else if (t<0. && t>-1.) { t *= 10.; expfac *= 10.; }
	  else if (t>1000.)  { t /= 10.; expfac /= 10.; }
	  else if (t<-1000.) { t /= 10.; expfac /= 10.; }
	  else break;
	}
      /* *offset = (float)(floor(t) / expfac);*/
      *offset = (float)( t / expfac);
      xmin -= *offset;
      xmax -= *offset;
      goto START;
    }
  min = xmin * (float)expfac;			/* get min > -9.99 and... */
  max = xmax * (float)expfac;			/* ... max <  9.99 */
  i = (figures > 3) ? 0 : 3;
  for(ni=0; i<NI; i++)
    {
      neartick(min, i, &n1);
      neartick(max, i, &n2);
      if (n1==0 && n2==0) break;
      n = n2 - n1;
      if (ni==0) ;			/* 1st one: allow anything */
      else if (ni>=8 && n<=8) ;		/* 10 is too many always */
      else if (i>=9 && n>=4 && n<=8) ;	/* 1 significant figure, 8 ok */
      else if (i>=6 && i<9 &&		/* 2 or 1 significant figure... */
	       n<=6 && n>=4) ;		/* ... 4 to 6 is perfect */
      else continue;
      ni = n; ii = i;
    }

  neartick(min, ii, &n1);
  neartick(max, ii, &n2);
  new[0] = (float)(n1-1) * (float)(interval[ii] / expfac) ;
  new[1] = (float)(n2+1) * (float)(interval[ii] / expfac) ;
  *dx = (float)(interval[ii] / expfac);
  for(x=new[0]+*offset; x<=new[1]+*offset; x+=*dx)
    {
      if (x >= lim[0] && x <= lim[1]) n++;
    }
  if (n <= 2) goto NEWOFF;
  ii /= 3;
  if (ii > 3) ii = 3;
  return format[ii];
}

static int lsp, wsp;

/*-----------------------------------------------------------------------------
|   label. places tick marks and numbers on graphs
-----------------------------------------------------------------------------*/
void label(float *xlim, float xscale, float xoffset,
	   float *ylim, float yscale, float yoffset, CURVE_SET *cp)
{
  float xleft, ytop;
  char string[100], *sp, *xname, *yname;
  int x0, y0, xc, yc, i;

  xname = cp->ix.label;
  yname = cp->iy.label;
  if (ps_modeon)
    {
      ps_save(1);
      ps_normalline();
      ps_setlabelfont(0);
    }

/*----------- draw and label tick marks */
  if (!ps_modeon)
    XSetForeground(mydisplay, redraw_gc, white());
  xleft = xlim[0];
  ytop  = ylim[1];

  drawticks(0, xlim, xlim[1], ylim, ylim[0], xscale, xoffset, yscale, yoffset);
  drawticks(1, ylim, ytop,    xlim, xleft,   yscale, yoffset, xscale, xoffset);

  if (ps_modeon) ps_save(0);

/*----------------- draw title on x-axis */
  if (ps_modeon) ps_setlabelfont(1);
  sp = xname;
  lsp = strlen(sp);
  xc = (int) (xscale * (xlim[0] + xlim[1]) / 2. + xoffset);
  y0 = (int) (yscale * ylim[0] + 2.0 * font_height + yoffset);
  if (!ps_modeon)
    {
      wsp = textwidth(sp, lsp);
      x0 = xc - wsp / 2;
      XDrawString(mydisplay, redraw_win, redraw_gc, x0, y0, sp, lsp);
    }
  else
    ps_centerstring(xc, y0, 0, sp, lsp);

/*--------------- draw title on y-axis */
  sp = yname;
  lsp = strlen(sp);
  x0 = 10;
  yc = (int) (yscale * (ylim[0] + ylim[1]) / 2.0 + yoffset);
  y0 = yc - lsp * font_height / 2;
  if (y0 < font_height)
    y0 = font_height;
  if (!ps_modeon)
    for (i = 0; i < lsp; i++)
      {
	wsp = textwidth(sp, 1);		/* width of 1st char */
#ifndef UNIX
	wsp = 0;
#endif
	XDrawString(mydisplay, redraw_win, redraw_gc, x0 - wsp / 2, y0, sp++, 1);
	y0 += font_height;
      }
  else
    ps_centervertically(0, yc, 1, sp, lsp);

/*------ draw subtitle */
  if (ps_modeon) return;		/* subtitle for ps: see event(), 'p' */
  if (!parse_subtitle(cp, string)) return;

  sp = string;
  lsp = strlen(sp);
  /* xc = ... already set */
#ifdef UNIX
  y0 = font_height;			/* Looks good on condor */
#else
  y0 = 2 * font_height;			/* looks good on PC */
#endif
  if (!ps_modeon)
    {
      wsp = textwidth(sp, lsp);
      x0 = xc - wsp / 2;
      XDrawString(mydisplay, redraw_win, redraw_gc, x0, y0, sp, lsp);
    }
  else
    ps_centerstring(xc, y0, 0, sp, lsp);
}

/*-----------------------------------------------------------------------------
|	drawticks
-----------------------------------------------------------------------------*/
static void drawticks(int is_y,
		float xlim[], float xmax, float ylim[], float ymin,
	       	float xscale, float xoffset, float yscale, float yoffset)
{
  float xnew[2], offlx, x, x1, dx, dy;
  int iexpx, xv, yv, xv1, yv1, xv2, yv2, del, xc, yc;
  char string[40], *format;
  
  format = ticks(xlim, xnew, &dx, &iexpx, &offlx);
  dy = (ylim[1] - ylim[0]) * .02;		/* height of ticks in y... */
  del = (int) (yscale * dy + FH);		/* ...and in screen coords */
  /*yv =(int) (yscale * ymin + yoffset);*//* bottom of ticks, screen coord*/
  yv = is_y ? wx0 : wy0 + h;

  for (x1 = xnew[0]; x1 <= xnew[1]; x1 += dx)
    {
      if (fabs((double)(x1/dx)) < .001) x1 = (float)0;
      x = x1 + offlx;
      if (x >= xlim[0] && x <= xmax)
        {
	  xv = (int)(xscale * x + xoffset);
	  if (!is_y)				/* tick is xv1,yv1..xv2,yv2*/
	    {
	      xv1 = xv2 = xv;
	      yv1 = yv; yv2 = yv+del;
            }
	  else
	    {
	      yv1 = yv2 = xv;
	      xv1 = yv; xv2 = yv + del;
            }
	  sprintf(string, format, (x-offlx) * pow(10., (double) (-iexpx)));
	  lsp = strlen(string);
	  if (!ps_modeon)
	    {
	      XDrawLine(mydisplay, redraw_win, redraw_gc, xv1, yv1, xv2, yv2);
	      wsp = textwidth(string, lsp);
	      if (!is_y)
	        {
	          xc = xv - wsp / 2;		/* x1 of string */
	          if (xc + wsp > wx2) xc = wx2 - wsp;
		  yc = yv + font_height;
	        }
	      else
	        {
		  xc = xv1 - (wsp + 5);
		  yc = yv1 + (int)(.3 * font_height);
	        }
	      XDrawString(mydisplay, redraw_win, redraw_gc, xc,yc, string, lsp);
	    }
	  else
	    {
	      ps_line(xv1, yv1, xv2, yv2);
	      if (!is_y) ps_centerstring(xv1, yv1, -1, string, lsp);
	      else   ps_centervertically(xv1, yv1, -1, string, lsp);
	    }
        }
    }
/*--------- write exponents on x- and y- axis */
  if (iexpx != 0 || offlx != FZ)
    {
      xv = is_y ? 5 : (int) (xscale * xlim[1] + xoffset);
      if (is_y)
#ifdef UNIX      
        yv = 2 * font_height;
#else
        yv = 3 * font_height / 2;
#endif
      else	      
        yv = (int) (yscale * ylim[0] + yoffset + 2.25*font_height);/* was 2.0*/
      draw_ten(is_y+1, xv, yv, iexpx, offlx);
    }
}

/*-----------------------------------------------------------------------------
|   draw_ten -- draw the power of ten on the axis
|	WARN! for ps, should let IT calculate line width wsp
|	and perhaps want text vertical
-----------------------------------------------------------------------------*/
static void draw_ten(int x_or_y, int xv, int yv, int iexp, float off)
{
  int yvh, n;
  char string[80], *p, *p2;
  
  p = string;
  if (iexp)
    {
      sprintf(p, "x10%d",iexp); p += strlen(p); p2 = p;
    }
  if (off)
    {
      sprintf(p,"+ %g", off);
    }
  n = textwidth(string, strlen(string));
  if (x_or_y == 1) xv -= n;

  if (!ps_modeon)
    {
      p = string;
      yvh = yv - font_height / 2;
      if (!strncmp(p,"x10",3))
        {
	  XDrawString(mydisplay, redraw_win, redraw_gc, xv, yv, p, 3);
          n = textwidth(p, 3);    xv += n; p += 3;
	  XDrawString(mydisplay, redraw_win, redraw_gc, xv, yvh, p, p2-p);
	  n = textwidth(p, p2-p); xv += n; p = p2;
        }
      if (off)
	XDrawString(mydisplay, redraw_win, redraw_gc, xv, yv, p,strlen(p));
    }
  else if (iexp || off)
    ps_power(x_or_y, string);
}

