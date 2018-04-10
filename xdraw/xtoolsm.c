/******************************************************************************
**  NAME        XTOOLS
**  AUTHOR        Sheryl M. Glasser
**
**  DESCRIPTION
**     Supporting functions for xwindows
**
**  Copyright (c) CounterPoint Graphics 1993.  All rights reserved.
******************************************************************************/
/*
Bugs, other
1.  test_window in makewindow after XMapRaised: last x,y still 0,0; ok in
    execute_motif.  Find a way to do an expose event or something.

*/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#ifdef DOS
#include <graph.h>
#endif

#ifdef UNIX
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>

#else
#include <dos.h>
#include "xlib.h"
#include "xutil.h"
extern void Drawcursor(void);
extern void Drawcrosshair(int, int, int, int, int, int);
#endif

#ifdef MOTIF
#include <Xm/Xm.h>
#include <Xm/DrawingA.h>
#endif

#include "gendefs.h"
#include "curves.h"		/* (xcontour.h uses CURVE_SET) */
#include "xtools.h"
#include "xdraw.h"
#include "xcontour.h"
#include "xinit.h"
#include "xedit.h"
#include "setcolor.h"
#include "ps.h"

#ifdef MOTIF
Widget topLevel;
XtAppContext app;
#endif
/*=============================================================================
**            DEFINITIONS AND VARIABLES
**===========================================================================*/

VIEW view[MAXWINDOW];
extern Display *mydisplay;
extern int myscreen;

static Window parent;
static XSizeHints myhint;
unsigned long myforeground, mybackground;

Window root;
Font font;
XFontStruct *font_struct;
extern int font_height;
Colormap cmap;

int nwindow = 0;
static char *maintitle;
static int ncol, icol=0, base_dx, base_dy, maxwindow;
#ifdef UNIX
static int xsep = 20, ysep = 45;
#else
static int xsep = 20, ysep = 20;
#endif
unsigned int xscreen=1280, yscreen=1024;	/* set by XGetGeometry() */
static int argc;
static char **argv;

extern char datafile[];
extern CURVE_SET curveset[];

/*=============================================================================
**                    DISPLAY, WINDOW, GC INITIALIZATION
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	give_command_args -- save argc, argv for Motif
-----------------------------------------------------------------------------*/
void give_command_args(int argcc, char **argvv)
{
  argc = argcc;
  argv = argvv;
}

/*-----------------------------------------------------------------------------
|   opendisplay. initializes display_in, etc.
-----------------------------------------------------------------------------*/
#ifdef UNIX
#define FPT1 12
#define FPT2 18
#else
#define FPT1 18
#define FPT2 22
#endif

opendisplay(char *title_in, int n)
{
/*------ declarations */
  int x0, y0;
  unsigned int bwidth, depth;
  char finalname[100];

/*------- copy input data */
  maintitle = title_in;

/*------- initialization */
#if defined(XCALLS)
  mydisplay = XOpenDisplay("");
#elif defined(MSWINDOWS)
  mydisplay = XOpenDisplay(maintitle);
#elif defined(MOTIF)
  XtSetLanguageProc(NULL, NULL, NULL);
  topLevel = XtVaAppInitialize(&app, "XMenu", NULL,0,
			       &argc, argv, NULL, NULL);
  mydisplay = XtDisplay(topLevel);
#endif

  if (NOX)
    {
      printf("WARNING -- XTerm display not found.\n");
      goto NO_XTERM;
    }

  myscreen = DefaultScreen(mydisplay);
  cmap = DefaultColormap(mydisplay, myscreen);

/*------- set up font */
  if (!getgoodfont(finalname, FPT1, FPT2,	/* MSWINDOWS: dummy */
         "times-medium-r-normal")) return 0;
  font_struct = XLoadQueryFont(mydisplay, finalname);
  font = (*font_struct).fid;

/*------- default pixel values */
  mybackground = BlackPixel(mydisplay, myscreen);
  myforeground = WhitePixel(mydisplay, myscreen);

  root = parent = DefaultRootWindow(mydisplay);
  XGetGeometry(mydisplay, root, &root, &x0, &y0,
	       &xscreen, &yscreen, &bwidth, &depth);

NO_XTERM:
  ncol = n;
  maxwindow = ncol * ncol;
  if (ncol == 1)
    {
      base_dx = xscreen * 2 / 3;
      base_dy = yscreen * 2 / 3;
    }
  else
    {
      base_dx = (xscreen + 1 - (ncol + 1) * xsep) / ncol;
      base_dy = (yscreen + 1 - (ncol + 1) * ysep) / ncol;
    }
  return (!NOX);
}

/*-----------------------------------------------------------------------------
|	get_screensize
-----------------------------------------------------------------------------*/
void get_screensize(unsigned int *dx, unsigned int *dy)
{
  *dx = xscreen; *dy = yscreen;
}

/*-----------------------------------------------------------------------------
|	closedisplay
-----------------------------------------------------------------------------*/
void closedisplay()
{
#ifndef MOTIF
  int i;
  for (i = nwindow - 1; i >= 0; i--)
    {
      XDestroyWindow(mydisplay, view[i].window);
      XFreeGC(mydisplay, view[i].gc);
    }
  XCloseDisplay(mydisplay);
#endif
}

/*-----------------------------------------------------------------------------
|   makewindow. creates window and gc
-----------------------------------------------------------------------------*/
Window makewindow(char *title, int has_menu)
{
  int argc = 0;		/* for XSetStandardProperties */
  int xwindow, ywindow;
  char **argv, *p;
  unsigned long fore, back;
  Window topwin, win;
  GC gc;
  VIEW *v;
#ifdef MOTIF
  Widget shell, area;
  char text[8];
#endif
  extern GC dialog_gc;

/*------- default program-specified window position and size */
  nextwindow(&xwindow, &ywindow);
  icol++;
  myhint.width =  base_dx;
  myhint.height = base_dy;
  myhint.x = xwindow;	/* windows are positioned... */
  myhint.y = ywindow;	/* ..incrementally */
  topwin = parent;
  back = has_menu? myforeground : mybackground;
  fore = has_menu? mybackground : myforeground;
  myhint.flags = PPosition | PSize;

/*-------- window creation */

  v = view + nwindow;
  p = title;
  if (!NOX)
    {
#if defined(XCALLS)
      argv = NULL;		/* (unused) */
      XWarpPointer(mydisplay, None, root, 0,0,0,0, myhint.x, myhint.y);
      win = v->window =
	XCreateSimpleWindow(mydisplay, topwin, myhint.x, myhint.y,
			    myhint.width, myhint.height, 5, fore, back);
      XSetStandardProperties(mydisplay, win,
			   p, p, None,	/* None,icon titles; icon pixmap */
			   argv, argc, &myhint);
#elif defined(MSWINDOWS)
      win = v->window = xCreateMainWindow(p, has_menu, myhint.x, myhint.y,
					  myhint.width, myhint.height);
#elif defined(MOTIF)
      shell =
	XtVaAppCreateShell(NULL, "XDraw", topLevelShellWidgetClass,
			   mydisplay, XmNx, myhint.x, XmNy, myhint.y,
			   XmNwidth, myhint.width, XmNheight, myhint.height,
			   XmNtitle, title, NULL);
      sprintf(text, "Win%d", nwindow);
      area =
	XtVaCreateWidget(text, xmDrawingAreaWidgetClass, shell, NULL);
      enableWidget(nwindow, shell, area, text);
      win = v->window = XtWindow(area);
#endif
    }

  nwindow++;
  if (NOX)
    return ((Window) 0);

/*---------- gc creation and initialization */
  gc = v->gc = XCreateGC(mydisplay, win, (long) 0, 0);
  if (has_menu) dialog_gc = gc;
  XSetBackground(mydisplay, gc, back);	/* (MSWINDOWS these unused) */
  XSetForeground(mydisplay, gc, fore);
  XSetFont(mydisplay, gc, font);

/*---------- input event selection */
  set_winmask(win);

/*---------- generate expose event */
  XMapRaised(mydisplay, win);
#ifdef MOTIF
  if (nwindow==1) test_window("makewindow", win, shell);
#endif

  return (win);
}
/*-----------------------------------------------------------------------------
|	nextwindow *** ERROR here
-----------------------------------------------------------------------------*/
void nextwindow(int *xw, int *yw)
{
  int i,n;
  n = ncol * ncol;
  i = (icol >= n) ? n-1 : icol;
  *xw = (i % ncol) * (base_dx + xsep);
  *yw = (i / ncol) * (base_dy + ysep);
}

/*-----------------------------------------------------------------------------
|   getview, getwindow
-----------------------------------------------------------------------------*/
int getview(Window w)
{
  VIEW *v;
  int i;
  for (i = 0, v = view; i < nwindow && v->window != w; i++, v++);
  return (i);
}

Window getwindow(int i)
{
  return((view+i)->window);
}

unsigned long white()
{
  return (myforeground);
}

/*-----------------------------------------------------------------------------
|   get_properties -- called from redraw
-----------------------------------------------------------------------------*/
void get_properties(Window w, unsigned int *pw, unsigned int *ph,
		    unsigned int *pd)
{
  Window root;
  int x0, y0;
  unsigned int bwidth, status;

  if (NOX)
    {
      *pw = myhint.width;
      *ph = myhint.height;
     }
  else
    status = XGetGeometry(mydisplay, w, &root,&x0, &y0, pw, ph, &bwidth, pd);
}

#ifndef UNIX
/*-----------------------------------------------------------------------------
|	addfloat
-----------------------------------------------------------------------------*/
float huge *addfloat(float *buf, long ix0)
{
  float huge *p;
  long ix;
  p = buf;
  if (ix0 < 0) ;
  else if (ix0 < 0x4000) p += ix0;
  else {
    ix = ix0 & 0xfffff000L;
	p += (ix0 & 0xfff);
	FP_SEG(p) += (int)(ix >> 2);
	}
  return(p);
}
#endif		/* of ndef UNIX */
/*-----------------------------------------------------------------------------
|	textwidth
-----------------------------------------------------------------------------*/
int textwidth(char *string, int n)
{
  if (string==0) return 0;
  if (font_struct==0) return strlen(string);
  return (XTextWidth(font_struct, string, n));
}

