/******************************************************************************
**	NAME		XWINLIB.C
**	AUTHOR		Sheryl M. Glasser
**
**	DESCRIPTION
**		Input from Xutils.h, X.h, Xatom.h
**
**	Copyright (c) Toptools SCF 1990.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bios.h>
#include <dos.h>
#include <conio.h>
#include <malloc.h>

#include "xlib.h"
#include "xutil.h"
#include "xatom.h"
#include <windows.h>
#include "wdraw.h"

#define MAX 20

typedef unsigned char byte;

HPEN hwhitePen, hblackPen, hcurrentPen;
HBRUSH hblackBrush, hwhiteBrush;
HWND dialogWnd;
unsigned long rgbColor;

#define FDX	40			/* dx for font0 */
#define FDY	60			/* dy for font0 */
#define FS	16			/* sep for font0 */
#define ROOTWIN 10
#define XSCREEN 1279
#define YSCREEN 1023

static int tdx=FDX, tdy=FDY, tsep=FS;
static int current_textsize=-1, textwidth=16;
static int current_index=-1;
static int White;

static int gin_enabled=0,is_mouse=-1;
static int enternotify_win = 0;
static char ginstring[20];
static int xcurs,ycurs,ncurs=0;
static int xscreen,yscreen;
static Cursor cursid[20],current_curs=0,new_curs=0;

static int ngc=0;
static struct _XGC xgc[MAX], *current_gc=0;

struct WIN
{
Window id;
HWND hWnd;
HDC hdc;
unsigned char used,expose,big;
};

static Window current_id;
static int nwin=0, current_view=-1;
static int wx1,wy1;
static struct WIN win[MAX];		// root in win[0] is never drawn
static char fontname[100],*fontnamep[5];
static XFontStruct fontstruct;
static long eventmask;
static XColor bkgd;

static struct _XDisplay display;
Screen screenlist[2];

#define F 0xff
#define gc2dc(gc) (HDC)gc->gid
/*=============================================================================
**			SERVICE FUNCTIONS
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	getwin
-----------------------------------------------------------------------------*/
static int getwin(Window w)
{
struct WIN *p,*p2;
int i;
for(p=win,p2=p+nwin; p<p2; p++)
	if (p->used && p->id==w) break;
i = (int)((p<p2) ? p - win : 1);
return i;
}

/*=============================================================================
**			INITIALIZE AND TERMINATE DISPLAY
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	XOpenDisplay, XCloseDisplay
-----------------------------------------------------------------------------*/
Display *XOpenDisplay(const char *aName)
{
   register struct WIN *p;
   Screen *s;
   winRegister(aName);		// (exits if can't register)
   
   xscreen = GetSystemMetrics(SM_CXFULLSCREEN);
   yscreen = GetSystemMetrics(SM_CYFULLSCREEN);
   xscreen--;
   yscreen--;

   current_gc=0;
   p=&win[0];				// initialize default (root) window
   p->id = (Window)ROOTWIN;		// window id's will be 10,11,etc.
   p->used=1;
   p->big=0;
   p->expose=0;				// Don't redraw root window
   nwin++;

// Set screen to grey

   s = display.screens = screenlist;
   s->root = p->id;

   hwhitePen   = GetStockObject(WHITE_PEN);
   hblackPen   = GetStockObject(BLACK_PEN);
   hwhiteBrush = GetStockObject(WHITE_BRUSH);
   hblackBrush = GetStockObject(BLACK_BRUSH);

   // set cursid[0] as value 1, define cursor for ROOTWIN as cursid[0]

   return(&display);
}

int XCloseDisplay(Display *disp)
{
	disp++;
	return(0);
}

/*=============================================================================
**			WINDOWS
**===========================================================================*/
//-----------------------------------------------------------------
//	xCreateMainWindow
//	* replaces XCreateSimpleWindow()
//	* calls CreateMainWindow(), which calls CreateWindow()
//-----------------------------------------------------------------
Window xCreateMainWindow(char *name, int has_menu,
			int maxWin, int x, int y, int dx, int dy)
{
   HWND hWnd;
   int w;
   register struct WIN *p;

   x++; y++; dx++; dy++;
   hWnd = CreateMainWindow(name, has_menu, maxWin);
   for(w=0,p=win; w<nwin && p->used; w++,p++) ;		// find empty win[]
   if (w==nwin) nwin++;
   p->used=1;
   p->expose=0;
   p->big=0;
   p->id = (p-1)->id +1;
   p->hWnd = hWnd;
   if (has_menu) dialogWnd = hWnd;
   return(p->id);
}

/*-----------------------------------------------------------------------------
|	XChangeProperty -- change window name
-----------------------------------------------------------------------------*/
XChangeProperty(Display *d, Window w, Atom property, Atom type,
	int format, int mode, const unsigned char *data, int nelements)
{
	register struct WIN *p;
	d++; type++; format++; mode++; nelements++;
	if (property != XA_WM_NAME) return 0;
	p = win + getwin(w);
	SetWindowText(p->hWnd, (LPSTR)data);
	return 0;
}

/*-----------------------------------------------------------------------------
|	XDestroyWindow
-----------------------------------------------------------------------------*/
XDestroyWindow(Display *d, Window w)
{
struct WIN *p;
HWND hWnd;
HDC hdc;
d++;
p = win + getwin(w);
hWnd = p->hWnd;
hdc = p->hdc;
ReleaseDC(hWnd, hdc);
DestroyWindow(hWnd);
return(1);
}

XMoveResizeWindow(Display *d, Window w, int x,int y,
				  unsigned int dx, unsigned int dy)
{ d++; w++; x++; y++; dx++; dy++; return 1; }

/*-----------------------------------------------------------------------------
|	  XGetGeometry
-----------------------------------------------------------------------------*/
Status XGetGeometry(Display *d, Drawable w, Window *root,
	int *x, int *y, unsigned int *dx, unsigned int *dy,
	unsigned int *wbord, unsigned int *depth)
{
register struct WIN *p;
RECT rect;
int i;
d++; root++;
i=getwin((Window)w);
p = win + i;
GetWindowRect(p->hWnd, &rect);
*x = rect.left;
*y = rect.top;
*dx = rect.right - *x;
*dy = rect.bottom - *y;
*wbord = 1;
*depth = 0;
return(0);
}

/*-----------------------------------------------------------------------------
|	XMapRaised -- generate an expose event (post WM_PAINT)
-----------------------------------------------------------------------------*/
XMapRaised(Display *d, Window w)
{
d++;
UpdateWindow((win+getwin(w))->hWnd);
return(1);
}

/*-----------------------------------------------------------------------------
|	XClearArea -- clear a rectangle (set rect, post WM_PAINT)
-----------------------------------------------------------------------------*/
int XClearArea(Display *d, Window w, int x, int y,
	       unsigned int width, unsigned int height, Bool expose)
{
d++;
if (x==0 && y==0 && width==0 && height==0)	
   InvalidateRect((win+getwin(w))->hWnd, NULL, expose?TRUE:FALSE);
return 1;
}


/*-----------------------------------------------------------------------------
|	XRaiseWindow, XLowerWindow
-----------------------------------------------------------------------------*/
XRaiseWindow(Display *d, Window w)
	{ d++; BringWindowToTop((win+getwin(w))->hWnd); return 1; }
XLowerWindow(Display *d, Window w) { d++; w++; return 1; }

/*-----------------------------------------------------------------------------
|	XSetClipRectangles
-----------------------------------------------------------------------------*/
XSetClipRectangles(Display *d, GC g, int x, int y,
		   XRectangle *r, int n, int o)
{
int x1,y1,x2,y2;
HRGN hRgn;
d++; g++; o++; x++; y++;
if (n==0) return 1;
x1 = r->x;				// this one if logorg active
y1 = r->y;
x2 = x1 + r->width;
y2 = y1 + r->height;
hRgn = CreateRectRgn(x1, y1, x2, y2);
SelectClipRgn((HDC)g->gid, hRgn);
return 1;
}

/*=============================================================================
**				GRAPHICS CONTEXT
**===========================================================================*/
GC XCreateGC(Display *d, Drawable w, unsigned long mask, XGCValues *v)
{
   register int i;
   struct WIN *p;
   GC g;
   HDC hdc;
   HWND hwnd;
   d++; w++;

   p = win + getwin(w);
   for(i=0,g=xgc; i<ngc && g->dirty; i++,g++) ;
   if (i==ngc) ngc++;
   g->dirty=1;
   g->values.foreground = (mask || GCForeground) ? v->foreground : 0;
   hcurrentPen = hwhitePen;
   hwnd = p->hWnd;
   hdc = GetDC(hwnd);
   p->hdc = hdc;
   SelectObject(hdc, hwhitePen);
   SelectObject(hdc, hblackBrush);
   SelectObject(hdc, GetStockObject(SYSTEM_FONT));
   if (hwnd==dialogWnd)
   {
      SetMapMode(hdc, MM_TEXT);
   }
   else
   {
      SetMapMode(hdc, MM_ANISOTROPIC);
      SetWindowExt(hdc, 1280, 1024);
      SetViewportExt(hdc, GetDeviceCaps(hdc,HORZRES),
		          GetDeviceCaps(hdc,VERTRES));
      SetViewportOrg(hdc, 0, 0);
   }
   g->gid = (GContext)hdc;
   current_gc = g;
   return(g);
}

#ifdef DEAD_CODE
XFreeGC(Display *d, GC gc)
{
d;gc;
if (gc-xgc < ngc)
	gc->dirty=0;
return(1);
}
#endif

XSetFont(Display *d, GC gc, Font f)
{
d++;
gc->values.font = f;
//ready_textsize(gc);
return(1);
}

XChangeGC(Display *d, GC gc, unsigned long mask, XGCValues *v)
{
//xor_ends(1);
d++; gc++; mask++;
xormode((v->function == GXxor) ? (byte)1 : (byte)0);
return 1;
}

/*=============================================================================
**				FONT INFO
**===========================================================================*/
static char *doubledash(const char *pat)
{
register char *p;
for(p=(char *)pat; p;)					// find "--" in pattern
	{
	p=strchr(p,'-');
	if (p && *(++p)=='-') break;
	}
if (p) p++;
return(p);
}

XFontStruct *XLoadQueryFont(Display *d, const char *name)
{
/*register char *p;
Font i;
if (p=doubledash(name))	sscanf(p,"%ld",&i);
else	i=12;
fontstruct.fid=i;*/
d++; name++;
return(&fontstruct);
}

char **XListFonts(Display *d, const char *pat, int maxnames, int *count)
{
#ifdef DEAD_CODE	
register char *p, *q;
strcpy(fontname,"Standard--");		// only fonts here are Standard--xx
*count=0;
if (p=doubledash(pat))
	{
	if (q=strchr(p,'*')) *q=0;
	strcat(fontname,p);
	*count=1;
	}
#endif	
d++; pat++; maxnames++; count++;
return(fontnamep);
}

XFreeFontNames(char **list)
{
//ok -- list was not allocated
list++;
return 1;
}

XTextWidth(XFontStruct *f, const char *s, int count)
{
long w;
f++;
w = GetTextExtent(win->hdc, s, count);
return (int)w & 0xffff;
}

/*=============================================================================
**			DRAW -- note 0,0 is upper left
**===========================================================================*/
XDrawPoint(Display *d, Drawable w, GC gc, int x1,int y1)
{
HDC hdc;
d++; w++;
hdc = gc2dc(gc);
SetPixel(hdc, x1, y1, rgbColor);
return(1);
}

XDrawLine(Display *d, Drawable w, GC gc, int x1,int y1, int x2, int y2)
{
HDC hdc;
d++; w++;
hdc = gc2dc(gc);
MoveTo(hdc, x1, y1);
LineTo(hdc, x2, y2);
return(1);
}

XDrawLines(Display *d, Drawable w, GC gc, XPoint points[],
	   int npoints, int mode)
{
  int i,x,y;
  XPoint *p;
  HDC hdc;
  d++; w++; mode++;
  hdc = gc2dc(gc);
  for(i=0, p=&points[0]; i<npoints; i++,p++)
    {
      x = p->x;
      y = p->y;
      if (i==0)
	MoveTo(hdc, x, y);
      if (i>0)
	LineTo(hdc, x, y);
    }
    return 1;
}

XDrawRectangle(Display *d, Drawable w, GC gc, int x, int y,
	unsigned int dx, unsigned int dy)
{
   int x2,y2;
   HDC hdc;
   d++; w++;
   hdc = gc2dc(gc);
   x2=x+dx;
   y2=y+dy;
   Rectangle(hdc, x, y, x2, y2);
   return 1;
}

XFillRectangle(Display *d, Drawable w, GC gc, int x, int y,
	unsigned int dx, unsigned int dy)
{
  int x2,y2;
  HBRUSH hBrush;
  HDC hdc;
  d++; w++;
  hdc = gc2dc(gc);
  hBrush = GetStockObject(GRAY_BRUSH);	// Temporary!
  SelectObject(hdc, hBrush);
  x2=x+dx;
  y2=y+dy;
  Rectangle(hdc, x, y, x2, y2);
  SelectObject(hdc, GetStockObject(NULL_BRUSH));
  return 1;
}

#ifdef DEAD_CODE
XFillArc(Display *d, Drawable w, GC gc, int x, int y,
	unsigned int dx, unsigned int dy, int a1, int a2)
{
  struct WIN *p;
  int x2,y2;
  if (a2 != a1 + 360*64) return;
  p = win + ready_draw((Window)w, gc);
  x2 = x+dx;
  y2 = y+dy;
  pc_coord(x,y);
  pc_coord(x2,y2);
  _ellipse(_GFILLINTERIOR, x,y,x2,y2);
}

XFillPolygon(Display *d, Drawable w, GC gc, XPoint points[], int npoints,
	     int shape, int mode)
{
  int x,y,i;
  int xmin,ymin,xmax,ymax;
  if (shape==Complex) return;
  XDrawLines(d,w,gc,points,npoints,0);
  xmin = xscreen; xmax = 0;
  ymin = yscreen; ymax = 0;
  for(i=x=y=0; i<npoints-1; i++)	// centered x,y
    {					// ignore last point == first
      x = points[i].x;
      y = points[i].y;
		if (x<xmin) xmin=x;
      if (y<ymin) ymin=y;
      if (x>xmax) xmax=x;
      if (y>ymax) ymax=y;
    }
  pc_coord(xmin, ymin);
  pc_coord(xmax, ymax);
  x = (xmin + xmax) / 2;
  y = (ymin + ymax) / 2;
  x = xmin + 1;				// VERY specific
  _floodfill(x,y,current_index);
}
#endif

XDrawString(Display *d,Drawable w,GC gc,
	int x, int y, const char *s, int length)
{
char text[100];
HDC hdc;
d++; w++;
hdc = gc2dc(gc);
strncpy(text,s,length);
*(text+length)=0;
TextOut(hdc, x, y, text, strlen(text));
return 1;
}


/*=============================================================================
**				COLORS
**===========================================================================*/
Colormap XDefaultColormap(Display *d, int s)
{
d++; s++;
return(0);
}

XAllocColor(Display *d, Colormap m, XColor *x)
{
   d++; m++;
   x->pixel = rgbColor = RGB(x->red, x->green, x->blue);
   return 1;
}

XSetForeground(Display *d, GC gc, unsigned long fore)
{
	d++;
	gc->values.foreground = fore;
	if (hcurrentPen != hwhitePen && hcurrentPen != hblackPen)
		DeleteObject(hcurrentPen);
	if (fore==0) hcurrentPen = hblackPen;
	else if (fore==0xffffffL) hcurrentPen = hwhitePen;
   else hcurrentPen = CreatePen(PS_SOLID, 3,
	RGB(fore>>16, (fore>>8)&F, fore&F));
   SelectObject(gc2dc(gc), hcurrentPen);
   return(1);
}

XSetBackground(Display *d, GC gc, unsigned long bkgd)
{
d++;
gc->values.background = bkgd;
return(1);
}

/*=============================================================================
**				EVENTS, I/O
**===========================================================================*/
static int status=0, laststatus;
static HWND EvHWnd;
static int EvType;
static char EvChar;
//------------------------------------------------------
//	createXEvent
//------------------------------------------------------
void createXEvent(HWND hWnd, int type, char c)
{
   EvHWnd = hWnd;
   EvType = type;
   EvChar = c;
}

//------------------------------------------------------
//	XNextEvent
//------------------------------------------------------
XNextEvent(Display *d, XEvent *x)
{
x->xany.display = &display;
x->xany.window = (win + getwin(EvHWnd))->id;
x->xany.type = EvType;
if (EvType == KeyPress) x->xkey.keycode = (int)EvChar;
d++; return 1;
}

XRefreshKeyboardMapping(XMappingEvent *x)
{
x++; return 1;
}

XSelectInput(Display *d, Window w, long event_mask)
{
d++; w++;
eventmask=event_mask;
return 1;
}

XLookupString(XKeyEvent *x, char *buf, int nbytes,		// in Xutil.h
	KeySym *k, XComposeStatus *cs)
{
nbytes++; k++; cs++;
*buf = x->keycode & 0xff;
if (*buf) return(1);
*(buf+1) = x->keycode>>8;
return(2);
}


/*=============================================================================
**					CURSOR CONTROL
**===========================================================================*/
static RECT getWindowRect(HWND hWnd)
{
   WINDOWPLACEMENT place;
   RECT rect;
   GetWindowPlacement(hWnd, &place);
   if (place.showCmd == SW_SHOWMAXIMIZED)
   {
      rect.left  = place.ptMaxPosition.x;
      rect.top   = place.ptMaxPosition.y;
      rect.right = rect.left;
      rect.bottom= rect.top;
   }
   else if (place.showCmd == SW_SHOWNORMAL)
   {
      return place.rcNormalPosition;
   }
   else rect.top = rect.bottom = rect.left = rect.right = 0;
   return rect;
}

XWarpPointer(Display *d, Window src, Window dst, int srcx, int srcy,
			 unsigned int srcw, unsigned int srch, int x, int y)
{
   struct WIN *p;
   RECT rect;
   d++; src++; srcx++; srcy++; srcw++; srch++;
   p = win + getwin(dst);
   rect = getWindowRect(p->hWnd);
   Drawcursor();
   xcurs = rect.left + x;
   ycurs = rect.top + y;
   SetCursorPos(xcurs, ycurs);
   Drawcursor();
   return 1;
}

Bool XQueryPointer(Display *d, Window w, Window *root, Window *child,
int *xr,int *yr, int *xw,int *yw, unsigned int *mask)
{
   struct WIN *p;
   int i;
   POINT pos;
   RECT rect;
   GetCursorPos(&pos);
   *xr = pos.x;
   *yr = pos.y;
   *root  = NULL;
   for(i=0,p=win; i<nwin; i++,p++)
   {
      rect = getWindowRect(p->hWnd);
      if (rect.top==rect.bottom && rect.left==rect.right) break;
      if (pos.x >= rect.left && pos.x <= rect.right &&
	  pos.y >= rect.top  && pos.y <= rect.bottom) break;
   }
   *xw = pos.x - rect.left;
   *yw = pos.y - rect.top;
   *child = p->id;
   d++; w++; mask++; return 1;
}

int XDefineCursor(Display *d, Window w, Cursor curs)
{
int dx,dy;
struct WIN *p;
#ifdef DEAD_CODE
if (gin_enabled)
	{
	new_curs=curs; return(0);
	}
else is_mouse =mouse(&dx,&dy,&status);
if (!is_mouse) return(0);
p=win+getwin(w);
p=win+getwin((Window)ROOTWIN);
xcurs = p->x + p->dx/2;
ycurs = p->y + p->dy/2;
drawcurs(curs);
current_curs = new_curs = curs;
gin_enabled=1;
#endif
d++; w++; curs++; return(0);
}

int XUndefineCursor(Display *d, Window w)
{
new_curs=cursid[0];
//gin_enabled=0;
//drawcurs();
}

Cursor XCreateFontCursor(Display *d, unsigned int shape)
{
int i;
d++;
for(i=0; i<ncurs; i++)
	if (cursid[i]==shape) return(cursid[i]);
cursid[ncurs++]=(Cursor)shape;
return(cursid[ncurs-1]);
}

