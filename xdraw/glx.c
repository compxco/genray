/******************************************************************************
**  NAME            GLX.C
**  AUTHOR          Sheryl M. Glasser
**
**  DESCRIPTION
**
**
**  Copyright (c) GlassWare 1994.  All rights reserved.
******************************************************************************/
#include <string.h>
#include <stdio.h>		/* for NULL? */
#ifdef UNIX
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>

#else
#include "xlib.h"
#include "xutil.h"
#endif

#include "glx.h"
#include "device.h"
#include "setcolor.h"		/* for getcolor() */
#include "xtools.h"		/* for get_eventmask() */
#include "gendefs.h"

static Display *mydisplay;
static XSizeHints menuhint;
static Window root, menuwin;
static XEvent *evp;
static GC gc;
static Font font;
static XFontStruct *font_struct;
static unsigned long fore, back;
static int xtext, ytext;
static long lastcolor;
static int eventmask;

#define NONE 0
#define LINE 1
#define POLY 2
#define SIZE 100

static int mode = 0;
static int npoints = 0;
static XPoint points[SIZE];
static int descender=0, fontheight=0;
static unsigned int mdx, mdy;
static int xscreen, yscreen;

/*-----------------------------------------------------------------------------
|	init_glx
-----------------------------------------------------------------------------*/
void init_glx(Display *d, Window r, unsigned long f, unsigned long b,
	      Font fon, XFontStruct *fonst, XEvent *p, int dx, int dy)
{
  mydisplay = d;
  root = r;
  fore = f;
  back = b;
  font = fon;
  font_struct = fonst;
  evp = p;
  eventmask = get_eventmask();
  xscreen = dx;
  yscreen = dy;
}

/*-----------------------------------------------------------------------------
|	prefposition
-----------------------------------------------------------------------------*/
void prefposition(int mx1, int mx2, int my1, int my2)
{
  menuhint.x = mx1; menuhint.width  = mx2 - mx1;
  menuhint.y = my2; menuhint.height = my2 - my1;
  menuhint.flags = PPosition | PSize;
}

/*-----------------------------------------------------------------------------
|	winopen
-----------------------------------------------------------------------------*/
long winopen(char *title)
{
  int argc = 0;
  char **argv;
  int x,y;
  unsigned int bw, depth;
  Window root2;

  if (NOX) return(0L);
  XWarpPointer(mydisplay, None, root, 0, 0, 0, 0, menuhint.x, menuhint.y);
  menuwin =
    XCreateSimpleWindow(mydisplay, root, menuhint.x, menuhint.y,
			menuhint.width, menuhint.height, 5, fore, back);
  XSetStandardProperties(mydisplay, menuwin, title, title, None,
			 argv, argc, &menuhint);
  gc = XCreateGC(mydisplay, menuwin, (long)0, 0);
  XSetBackground(mydisplay, gc, back);
  XSetForeground(mydisplay, gc, fore);
  XSetFont(mydisplay, gc, font);
  XSelectInput(mydisplay, menuwin, (long)eventmask);
  XMapRaised(mydisplay, menuwin);
  XGetGeometry(mydisplay, menuwin, &root2, &x, &y, &mdx, &mdy, &bw, &depth);
  return(menuwin);
}

/*-----------------------------------------------------------------------------
|	winclose
-----------------------------------------------------------------------------*/
void winclose(long win)
{
  win;
}

/*-----------------------------------------------------------------------------
|	getorigin
|	* puts out 0,0 for ox,oy after 1st couple of times: useless
-----------------------------------------------------------------------------*/
void getorigin(long *x, long *y)
{
  Window root2;
  int ox,oy;
  unsigned int dx,dy,bw,depth;
  Status ok;
  ok = XGetGeometry(mydisplay, menuwin, &root2,
		    &ox, &oy, &dx, &dy, &bw, &depth);
  *x = (long)ox;
  *y = (long)(yscreen - oy - mdy);
  *x = *y = 0L;
  /*printf("Origin: win=%lx from root=%lx, at %d %d, ok=%d\n",
	 menuwin,root2,ox,oy, ok);*/
}

/*-----------------------------------------------------------------------------
|	maxsize, minsize, winconstraints
-----------------------------------------------------------------------------*/
void maxsize(long dx, long dy) { dx; dy; }
void minsize(long dx, long dy) { dx; dy; }
void winconstraints(void) {}
/*-----------------------------------------------------------------------------
|	windepth, winpop, winpush
-----------------------------------------------------------------------------*/
int windepth (long win) { win; return(0); }
void winpop(void) { XRaiseWindow(mydisplay, menuwin); }
void winpush(void) { XLowerWindow(mydisplay, menuwin); }
/*-----------------------------------------------------------------------------
|	winset, winget
-----------------------------------------------------------------------------*/
#ifdef UNIX
long winget(void) { return 0; }
void winset(long win) { win; }

#else
long winget(void)
{
  int xw,yw, xr,yr, mask;
  Window curswin, root1;
  XQueryPointer(mydisplay, root,
	  &root1, &curswin, &xr, &yr, &xw, &yw, &mask);
  return(curswin);
}

void winset(long win)
{
  int xw,yw, dx,dy;
  int bwidth, depth;
  XGetGeometry(mydisplay, win, &root, &xw, &yw,
	       &dx, &dy, &bwidth, &depth);
#ifndef UNIX
  xw = -2;
#endif
  XMoveResizeWindow(mydisplay, win, xw,yw, dx,dy);
}
#endif
/*-----------------------------------------------------------------------------
|	cpack
-----------------------------------------------------------------------------*/
void cpack(long color)
{
  int r,g,b;
  r = (int)((color & 0xff) << 8);
  g = (int)(((color>>8) & 0xff) << 8);
  b = (int)(((color>>16) & 0xff) << 8);
  lastcolor = getcolor(r,g,b);
}

long glx_index(long color)
{
  cpack(color);
  return(lastcolor);
}

static void set_gc_color(void)
{
  XSetForeground(mydisplay, gc, lastcolor);
  XSetBackground(mydisplay, gc, lastcolor);
}

/*-----------------------------------------------------------------------------
|	clear
-----------------------------------------------------------------------------*/
void clear(void)
{
  XSetBackground(mydisplay, gc, lastcolor);
  XClearWindow(mydisplay, menuwin);
}

/*-----------------------------------------------------------------------------
|	rectfi, recti
-----------------------------------------------------------------------------*/
void rectfi(int x1, int y1, int x2, int y2)
{
  int y;
  set_gc_color();
  if (y1 > y2) { y=y1; y1=y2; y2=y; }
  XFillRectangle(mydisplay, menuwin, gc, x1, mdy-y2, x2-x1, y2-y1);
}

void recti(int x1, int y1, int x2, int y2)
{
  int y;
  set_gc_color();
  if (y1 > y2) { y=y1; y1=y2; y2=y; }
  XDrawRectangle(mydisplay, menuwin, gc, x1, mdy-y2, x2-x1, y2-y1);
}

/*-----------------------------------------------------------------------------
|	drawstr, cmov2i, getstrwidth
-----------------------------------------------------------------------------*/
void drawstr(char *s)
{
  int dx, n;
  set_gc_color();
  n = strlen(s);
  XDrawString(mydisplay, menuwin, gc, xtext, ytext, s, n);
  dx = XTextWidth(font_struct, s, n);
  xtext += dx;
}

void cmov2i(int x, int y)
{
  xtext = x;
  ytext = mdy-y;
}

long getstrwidth(char *text)
{
  int width;
  width = XTextWidth(font_struct, text, strlen(text));
  return((long)width);
}

/*-----------------------------------------------------------------------------
|	getdescender, getheight
-----------------------------------------------------------------------------*/
static void fontextent(void)
{
  char alpha[128];
  int i, ascent, dir;
  XCharStruct overall;

  if (descender && fontheight) return;
  for(i=32; i<128; i++) alpha[i] = (char)i;
  XTextExtents(font_struct, alpha+32, 128 - 32, &dir,
	       &ascent, &descender, &overall);
  fontheight = ascent + descender;

}
int getdescender(void)
{
  fontextent();
  return(descender);
}

int getheight(void)
{
  fontextent();
  return(fontheight);
}

/*-----------------------------------------------------------------------------
|	bgnpolygon, endpolygon
-----------------------------------------------------------------------------*/
void bgnpolygon(void)
{
  mode = POLY;
  npoints = 0;
}

void endpolygon(void)
{
  set_gc_color();
  XFillPolygon(mydisplay, menuwin, gc, points, npoints,
	       Convex, CoordModeOrigin);
  mode = NONE;
}
/*-----------------------------------------------------------------------------
|	bgnline, endline
-----------------------------------------------------------------------------*/
void bgnline(void)
{
  mode = LINE;
  npoints = 0;
}

void endline(void)
{
  set_gc_color();
  XDrawLines(mydisplay, menuwin, gc, points, npoints, CoordModeOrigin);
  mode = NONE;
}

/*-----------------------------------------------------------------------------
|	v2i
|	(if want short int, should use v2s!)
-----------------------------------------------------------------------------*/
void v2i(long xy[])
{
  if (npoints==SIZE) return;
  points[npoints].x = (int)xy[0];
  points[npoints].y = mdy - (int)xy[1];
  npoints++;
}

/*-----------------------------------------------------------------------------
|	popmatrix
-----------------------------------------------------------------------------*/
void popmatrix(void) {}
void qdevice(int dev) { dev; }

/*-----------------------------------------------------------------------------
|	getvaluator
|	* values returned are relative to upper left
-----------------------------------------------------------------------------*/
int getvaluator(int dev)
{
  static int rx, ry, wx, wy;
  if (dev==MOUSEX)
    {
      rx = evp->xbutton.x_root;
      ry = evp->xbutton.y_root;
      wx = evp->xbutton.x;
      wy = evp->xbutton.y;
      return(wx);
    }
  if (dev==MOUSEY) return(mdy - wy);
  return(0);
}

/*-----------------------------------------------------------------------------
|	scale_font
-----------------------------------------------------------------------------*/
void scale_font(int scale) { scale; }

#ifdef NOT_YET
zap_events
add_event
#endif
