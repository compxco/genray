/******************************************************************************
**	NAME		XMSLIB.C
**	AUTHOR		Sheryl M. Glasser
**
**	DESCRIPTION
**		Input from Xutils.h, X.h, Xatom.h
**
**	Copyright (c) Toptools SCF 1990.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <bios.h>
#include <graph.h>
#include <dos.h>
#include <conio.h>
#include <malloc.h>

#include "xlib.h"
#include "xutil.h"
#include "xatom.h"

#define MAX 20

typedef unsigned char byte;

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
int x,y,dx,dy;					/* window based on 1279,1023 */
int back,fore,wbord;
unsigned char used,expose,big;
char title[40];
};

static Window current_id;
static int nwin=0, current_view=-1;
static int wx1,wy1;
static struct WIN win[MAX];					// root in win[0] is never drawn
static char fontname[100],*fontnamep[5];
static XFontStruct fontstruct;
static long eventmask;
static XColor bkgd;

static struct _XDisplay display;
Screen screenlist[2];

#define F 0xff
#define H 0x80
#define S 0x40

static int offct = 1*64;
static int colortable[] =
{
 0, 0,0,0,  1, S,S,S,  2, F,F,F,  3, F,F,0,	// black,grey,white,yellow
 4, H,F,0,  5, S,F,0,  6, 0,F,0,  7, 0,F,S,	// ..green
 8, 0,F,H,  9, 0,F,F, 10, 0,H,F, 11, 0,S,F,	// ..cyan
12, 0,0,F, 13, S,0,F, 14, H,0,F, 15, F,0,F,	// ..blue..magenta

 0, 0,0,0,  1, S,S,S,  2, F,F,F,  3, F,F,H,	// black,dkgrey,white,ltyel
 4, F,F,0,  5, H,F,S,  6, S,F,0,  7, 0,F,H,	// yellow,grn1,grn2,ltaqua
 8, 0,F,F,  9, H,F,F, 10, 0,H,F, 11, 0,S,F,	// cyan,ltblue,blue,dkblue
12, H,0,F, 13, H,0,H, 14, F,0,F, 15, F,H,F	// purple,dkmag,magenta,ltpink
};

/*=============================================================================
**			SERVICE FUNCTIONS ("PCTEK")
**===========================================================================*/
static long xsmul,ysmul,xsdiv,ysdiv,xsadd,ysadd;

#define pc_coord1(x,y) { x=(int)(((long)x*xsmul+xsadd)/xsdiv); \
			 y=(int)(((long)y*ysmul+ysadd)/ysdiv); }
#define pc_coord(x,y) { x=(int)(((long)x*xsmul)/xsdiv); \
			y=(int)(((long)y*ysmul)/ysdiv); }

#define select_view(i) vview(i,0)
#define renew_view(i)  vview(i,1)
#define delete_view(i) vview(i,2)
#define move(x,y) {pc_coord(x,y);_moveto(x,y);}
#define draw(x,y) {pc_coord(x,y);_lineto(x,y);}

/*-----------------------------------------------------------------------------
|	ready_textsize
-----------------------------------------------------------------------------*/
static void ready_textsize(GC gc)
{
  current_textsize = (int)gc->values.font;
  tdy = 52 * current_textsize / 24;	// get dx,dy,sep for the current
  tdx = 39 * current_textsize / 24;	// point for standard viewport
  tsep = tdy -tdx;			// 0,0..4095,3276
  textwidth = tdy;			// program has all chars same size
  tdy = tdy / 3;
  tdx = tdx / 3;
  tsep = tdy - tdx;
}

/*-----------------------------------------------------------------------------
|	vview
|	* how: 0=select_view, 1=renew_view, 2=delete_view
-----------------------------------------------------------------------------*/
static void vview(int i, int how)
{
  struct WIN *q,*p;
  int x1,y1,x2,y2,dx,dy;
  int cw,n,x,y;

  if (how==0 && current_view==i) return;
  current_view=i;
  q = p = win+i;
  if (q->big) p=win;
  current_id=q->id;
  x1=p->x;
  y1=p->y;
  x2=p->x + p->dx;
  y2=p->y + p->dy;
  pc_coord(x1,y1);
  pc_coord(x2,y2);
  if (how==0)						// select_view()
    {
      _setlogorg(0,0);
      _setcliprgn(x1,y1,x2,y2);
      _setlogorg(x1,y1);
      wx1=p->x; wy1=p->y;
      //_setviewport(x1,y1,x2,y2);
    }
  else if (how==2)				// delete_view()
    {
      q->fore=q->back=0;
      renew_view(i);
      q->used=0;
    }
  else							// renew_view()
    {
      dx = x2-x1;
      dy = y2-y1;
      _setcolor(q->back);
      _rectangle(_GFILLINTERIOR,0,0,dx,dy);
      _setcolor(q->fore);
      _rectangle(_GBORDER,0,0,dx,dy);
      ready_textsize(current_gc);
      cw = ((long)textwidth * 1280L *2L) / 4096L;
      n=strlen(q->title);
      x = p->dx - (n+1)*cw/2;
      if (x<0) x=cw/2/2;
      y = 3*cw/2/2;
      current_gc->values.foreground = q->fore;
      XDrawString(&display,q->id,current_gc,x,y,
		  q->title,n);		// WARN! should have a p->gc
    }
}

/*-----------------------------------------------------------------------------
|	colormap
-----------------------------------------------------------------------------*/
static void colormap(int n, int *p)
{
static long far list[16];
int i,j;
char *q;
for(i=0; i<n; i+=4)
	{
	j=*p++;
	q=(char *)&list[j];
	*q++ = (*p++) >> 2;
	*q++ = (*p++) >> 2;
	*q++ = (*p++) >> 2;
	}
_remapallpalette(list);
}

/*-----------------------------------------------------------------------------
|	ggetcolor
|		r,g,b are in 0..65535
-----------------------------------------------------------------------------*/
unsigned long ggetcolor(unsigned int r, unsigned int g, unsigned int b)
{
XColor xcolor;
Colormap cmap;
xcolor.red   = r;
xcolor.blue  = b;
xcolor.green = g;
XAllocColor(&display, cmap, &xcolor);
return(xcolor.pixel);
}

/*-----------------------------------------------------------------------------
|	getwin
-----------------------------------------------------------------------------*/
static int getwin(Window w)
{
struct WIN *p,*p2;
for(p=win,p2=p+nwin; p<p2; p++)
	if (p->used && p->id==w) break;
return((p<p2) ? p - win : 1);
}

/*-----------------------------------------------------------------------------
|	ready_draw
-----------------------------------------------------------------------------*/
static int ready_draw(Window w, GC gc)
{
  if (current_id != w)
    {
      select_view(getwin(w));
#ifdef DEAD_CODE
      current_view = getwin(w);
      current_id = w;
      select_view(current_view);
#endif
    }
  if (current_index != (int)gc->values.foreground)
    {
      current_index = (int)gc->values.foreground;
      _setcolor(current_index);
    }
  return(current_view);
}

/*-----------------------------------------------------------------------------
|	gtext
-----------------------------------------------------------------------------*/
void gtext(int x1, int y1, char *text)
{
   unsigned char *p, *q;
   int dx, dy, y0,xa,ya;
   char flag;
   extern int font0adr[];
   int far *pi;

   pi=font0adr;
   FP_SEG(q) = FP_SEG(pi);
   for (p = text; *p; p++)
   {
      FP_OFF(q) = font0adr[*p - 32];
      y0 = *q++;
      for (; ; )
      {
	 dx = (int)(*q++);
	 dy = (int)(*q++);
	 flag = (byte)(dy & 0xc0);
	 dy = (dy & 0x3f)-y0;
	 xa = x1 + dx*tdx/FDX;
	 ya = y1 - dy*tdy/FDY;
	 if (flag & 0x40)
		 draw(xa, ya)		// NOTE: absence of semicolon is OK
	 else
		 move(xa, ya)
	 if (flag & 0x80) break;
      }
      x1 += (FDX + FS) * tdx/FDX;
   }
}

/*-----------------------------------------------------------------------------
|	
-----------------------------------------------------------------------------*/
void monomode()
{
union REGS reg;
reg.h.ah=0;					// set mode
reg.h.al=3;
int86(0x10,&reg,&reg);

reg.h.ah=2;					// set cursor position
reg.h.bh=0;
reg.x.dx=0;
int86(0x10,&reg,&reg);

reg.h.ah=1;					// set cursor type
reg.x.cx=0x25b7;
int86(0x10,&reg,&reg);
}

void colorcard()
{
unsigned char far *p;
unsigned char byte;
union REGS reg;
FP_OFF(p) = 0x10;			// equipf
FP_SEG(p) = 0x40;
byte = *p & 0x30;
if (byte==0x30)
	{
	*p = (*p & 0xcf) | 0x20;	// change monitors
	monomode();
	}
}

/*-----------------------------------------------------------------------------
|	xormode
|	B$ResetEGA + (FA4-F85) = desired address to change    (mov ax,3)
-----------------------------------------------------------------------------*/
void xormode(int enable)
{
unsigned char byte;
int far *p;
static int old=0x3b8,new=0x7eb;
extern int _QCEnsScrSwap();
extern int far *ms_adr;

byte = enable ? (unsigned char)0x18 : 0;
outpw(0x3ce,3);
outp(0x3cf,byte);
p = ms_adr;
FP_OFF(p) += (0xfa4-0xf85);
if (enable) { old=*p; *p=new; }
else *p=old;
}

/*-----------------------------------------------------------------------------
|	drawcurs
-----------------------------------------------------------------------------*/
void drawcurs(Cursor curs)
{
int x1,x2,y1,y2,xc,yc,n;
#define N1 3
#define N2 6
n = (curs==cursid[0]) ? N1 : N2;
xc=xcurs;
yc=ycurs;
pc_coord(xc,yc);
x1=xc-n; if (x1<0) x1=0;
y1=yc-n; if (y1<0) y1=0;
x2=xc+n; if (x2>xscreen) x2=xscreen;
y2=yc+n; if (y2>yscreen) y2=yscreen;
_setcolor(White);
current_index=White;
if (curs==cursid[1])
	{
	x1=xc; y1=yc;
	}
else if (curs==cursid[2])
	{
	x2=xc; y2=yc;
	}
xormode(1);
_moveto(x1,yc);
_lineto(x2,yc);
_moveto(xc,y1);
_lineto(xc,y2);
xormode(0);
}

/*-----------------------------------------------------------------------------
|	Drawcurs
-----------------------------------------------------------------------------*/
void Drawcursor()
{
int i;
i = current_view;
select_view(0);
drawcurs(current_curs);
select_view(i);
}

/*-----------------------------------------------------------------------------
|	Drawcrosshair
-----------------------------------------------------------------------------*/
void Drawcrosshair(int x,int y,int x1,int y1,int x2,int y2)
{
static char far hline[324],far vline[1924];
char far *p, far *q;
int far *pi, far *qi;
static int nx=0,ny=0;
int nh,nv,ch,cv,i,j,i2;

pc_coord(x,y);
pc_coord(x1,y1);
pc_coord(x2,y2);
if ((x2-x1+1)!=nx || (y2-y1+1)!=ny)
	{
	nx=(x2-x1+1);
	ny=(y2-y1+1);
	nh=_imagesize(0,0,nx-1,0);
	nv=_imagesize(0,0,0,ny-1);
	p=hline; pi = (int *)p;
	q=vline; qi = (int *)q;
	*pi++ = nx; *pi = 1;
	*qi++ = 1;  *qi = ny;
	i2=(nx-1)/8;
	for(j=0,p+=4; j<4; j++)				// horiz line: 4 bit planes
		{
		ch= (j==1) ? 0xff: 0;
		for(i=0; i<=i2; i++)			// 0..i2 bytes
			{
			if (i==i2 && j==1) ch=0x80;
			*p++ = ch;
			}
		}
	
	for(j=0,q+=4; j<ny; j++)
		{
		for(i=0; i<4; i++)
			{
			cv= (i==1) ? 0x80: 0;
			*q++=cv;
			}
		}
	}
_putimage(x1,y,hline,_GXOR);
_putimage(x,y1,vline,_GXOR);
}

/*-----------------------------------------------------------------------------
|	mouse
-----------------------------------------------------------------------------*/
int mouse(int *dx, int *dy, int *status)
{
union REGS reg;
int ok,i;
if (is_mouse==-1)
	{
	for(i=ok=0; i<10 && !ok; i++)
		{
		reg.x.ax=0;
		ok=int86(0x33,&reg,&reg);
		}
	is_mouse = ok? 1 : 0;
	if (!ok) return(0);
	}
if (!is_mouse) return;
reg.x.ax=3; int86(0x33,&reg,&reg);
*status = reg.x.bx;
reg.x.ax = 11; int86(0x33,&reg,&reg);
*dx = reg.x.cx;
*dy = reg.x.dx;
return(1);
}
/*=============================================================================
**			INITIALIZE AND TERMINATE DISPLAY
**===========================================================================*/
Display *XOpenDisplay(const char *name)
{
register struct WIN *p;
XColor c;
Colormap m;
Screen *s;
struct videoconfig config;
name;

colorcard();
_setvideomode(_VRES16COLOR);
_getvideoconfig(&config);
_setlogorg(0,0);

xscreen = config.numxpixels;
yscreen = config.numypixels;
xsmul = (long)xscreen/32L;		// x' = x * 640/1280 = x * 20/40
xsdiv = (long)1280/32L;
xsadd =  xsdiv/2L;
ysmul = (long)yscreen/32L;		// y' = y * 480/1024 = y * 15/32
ysdiv = (long)1024/32L;
ysadd = ysdiv/2L;
xscreen--;
yscreen--;

current_gc=0;
p=&win[0];						// initialize default (root) window
p->id = (Window)ROOTWIN;		// window id's will be 10,11,etc.
p->x=p->y=0;
p->dx=1279;
p->dy=1023;
p->used=1;
p->big=0;
p->expose=0;					// Don't redraw root window
nwin++;

colormap(64,colortable+offct);			// initialize color map
bkgd.red = bkgd.green = bkgd.blue = 0x4000;		// clear screen to grey
XAllocColor(&display, m, &bkgd);
_setcolor(bkgd.pixel);
_rectangle(_GFILLINTERIOR,0,0,xscreen,yscreen);

display.default_screen = 0;
s = display.screens = screenlist;
s->root = p->id;

c.red = c.green = c.blue = 0xffff;			// init display.screens.white_pixel
XAllocColor(&display, m, &c);
White = s->white_pixel = c.pixel;

fontnamep[0]=fontname;
XCreateFontCursor(&display,1);
XDefineCursor(&display,(Window)ROOTWIN,cursid[0]);
return(&display);

}

int XCloseDisplay(Display *disp)
{
disp;

_setvideomode(_ERESNOCOLOR);
monomode();
return(0);
}

/*=============================================================================
**			WINDOWS
**===========================================================================*/
Window XDefaultRootWindow(Display *d)
{
d;
return(win[0].id);
}

Window XCreateSimpleWindow(Display *d, Window parent, int x, int y,
	unsigned int dx, unsigned int dy, unsigned int border_width,
	unsigned long border_color, unsigned long bkgd)
{
int w,i,x1,y1,x2,y2;
long parentx0,parenty0,parentdx,parentdy;
register struct WIN *p;
	
d;
p = win + getwin(parent);
parentx0 = (long)p->x;
parenty0 = (long)p->y;
parentdx = (long)p->dx;
parentdy = (long)p->dy;
if (parentx0 + parentdx > 1279) parentx0 = 1279 - parentdx;
if (parenty0 + parentdy > 1023) parenty0 = 1023 - parentdy;

for(w=0,p=win; w<nwin && p->used; w++,p++) ;		// find empty win[]
if (w==nwin) nwin++;
p->used=1;
p->expose=0;
p->big=0;
p->id = (p-1)->id +1;
p->fore = (int)border_color;
p->back = (int)bkgd;
p->wbord=border_width;

//dx = (int)((long)dx * parentdx/1279L);
//dy = (int)((long)dy * parentdy/1023L);
//x1 = parentx0 + x*(long)parentdx/1279L;
//y1 = parenty0 + y*(long)parentdy/1023L;
x1=x; y1=y;
p->x=x1; p->dx=dx;							// x,y,dx,dy are in 1279,1023
p->y=y1; p->dy=dy;
select_view(p-win);
return(p->id);
}

XSetStandardProperties(Display *d, Window w, const char *application,
	const char *icon, Pixmap pix, char **argv, int argc,
	XSizeHints *h)
{
register struct WIN *p;
p = win + getwin(w);
strcpy(p->title,application);
}

XChangeProperty(Display *d, Window w, Atom property, Atom type,
	int format, int mode, const unsigned char *data, int nelements)
{
register struct WIN *p;
if (property != XA_WM_NAME) return;
p = win + getwin(w);
strcpy(p->title, data);
}

XDestroyWindow(Display *d, Window w)
{
d;
delete_view(getwin(w));
return(1);
}

XMoveResizeWindow(Display *d, Window w, int x,int y,
				  unsigned int dx, unsigned int dy)
{
int i;
struct WIN *p,*q,*qbig;
q=qbig=0;
for(i=0,p=win; i<nwin; i++,p++)
	{
	if (!p->used) continue;
	if (p->id==w) q=p;
	if (p->big) qbig=p;
	}
i=q-win;
if (x==0 && q>win)				/* x=0 ==> F10 ==> make full screen */
	{
	if (qbig) qbig->big=0;
	q->big=1;
	Drawcursor();
	current_view = -1;
	select_view(i);
	renew_view(i);
	Drawcursor();
	XMapRaised(d,w);
	}
else if (x==-1 && qbig)				/* x=-1 ==> default size */
	{
	i=qbig-win;
	select_view(i);
	_setcolor(bkgd.pixel);
	_rectangle(_GFILLINTERIOR,0,0,xscreen,yscreen);
	qbig->big=0;
	current_view = -1;
	for(i=nwin-1; i>=0; i--)
		{
		p=win+i;
		if (p->id==(Window)ROOTWIN || !p->used) continue;
		select_view(i);
		renew_view(i);
		p->expose=1;
		}
	Drawcursor();
	}
else if (x==-2) select_view(i);		/* from winset() */
}


Status XGetGeometry(Display *d, Drawable w, Window *root,
	int *x, int *y, unsigned int *dx, unsigned int *dy,
	unsigned int *wbord, unsigned int *depth)
{
register struct WIN *p;
int i;
i=getwin((Window)w);
p = win + i;
if (p->big)
	p = win + getwin((Window)ROOTWIN);
select_view(i);
*x = p->x;
*y = p->y;
*dx = p->dx;
*dy = p->dy;
*wbord = 1;
*depth = 0;
return(0);
}

XMapRaised(Display *d, Window w)
{
// set priority of window w to top
register struct WIN *p;
int x,y,n,cw;
d;
current_view = getwin(w);
current_id = w;
select_view(current_view);
Drawcursor();
renew_view(current_view);
Drawcursor();
p = win + current_view;
p->expose=1;
return(1);
}

XClearWindow(Display *d, Window w)
{
// set priority of window w to top
register struct WIN *p;
int x,y,n,cw;
d;
current_view = getwin(w);
current_id = w;
select_view(current_view);
Drawcursor();
renew_view(current_view);
Drawcursor();
p = win + current_view;
return(1);
}

int XClearArea(Display *d, Window win, int x, int y,
               unsigned int w, unsigned int h, Bool expose)
{
XMapRaised(d,win);
}

XRaiseWindow(Display *d, Window w) { d; w; }
XLowerWindow(Display *d, Window w) { d; w; }

XSetClipRectangles(Display *d, GC g, int x, int y,	XRectangle *r,
				   int n, int o)
{
int x1,y1,x2,y2;
d;g;o;
if (n==0) return;
x1 = x + r->x;
y1 = y + r->y;
x1 = r->x;				// this one if logorg active
y1 = r->y;
x1 += wx1;
y1 += wy1;
x2 = x1 + r->width;
y2 = y1 + r->height;
pc_coord(x1,y1);
pc_coord(x2,y2);
_setcliprgn(x1,y1,x2,y2);
}

/*=============================================================================
**				GRAPHICS CONTEXT
**===========================================================================*/
GC XCreateGC(Display *d, Drawable w, unsigned long mask, XGCValues *v)
{
register int i;
GC g;
d;w;
for(i=0,g=xgc; i<ngc && g->dirty; i++,g++) ;
if (i==ngc) ngc++;
g->dirty=1;
current_gc = g;
g->values.foreground = (mask || GCForeground) ? v->foreground : 0;
return(g);
}

XFreeGC(Display *d, GC gc)
{
d;gc;
if (gc-xgc < ngc)
	gc->dirty=0;
return(1);
}

XSetForeground(Display *d, GC gc, unsigned long fore)
{
d;
gc->values.foreground = fore;
return(1);
}

XSetBackground(Display *d, GC gc, unsigned long bkgd)
{
d;
gc->values.background = bkgd;
return(1);
}

XSetFont(Display *d, GC gc, Font f)
{
d;
gc->values.font = f;
//ready_textsize(gc);
return(1);
}

XChangeGC(Display *d, GC gc, unsigned long mask, XGCValues *v)
{
//xor_ends(1);
xormode((v->function == GXxor) ? (byte)1 : (byte)0);
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
register char *p;
Font i;
if (p=doubledash(name))
	sscanf(p,"%ld",&i);
else
	i=12;
fontstruct.fid=i;
return(&fontstruct);
}

char **XListFonts(Display *d, const char *pat, int maxnames, int *count)
{
register char *p, *q;
strcpy(fontname,"Standard--");		// only fonts here are Standard--xx
*count=0;
if (p=doubledash(pat))
	{
	if (q=strchr(p,'*')) *q=0;
	strcat(fontname,p);
	*count=1;
	}
return(fontnamep);
}

XFreeFontNames(char **list)
{
//ok -- list was not allocated
}

XTextWidth(XFontStruct *f, const char *s, int count)
{
long w;
f;s;
/*------eg if textwidth is 18 point = 39 in 0..4096,
		then it's 39 * 1280 / 4096 in 0..1280 */
w = ((long)textwidth * 1280L) / 4096L;
return(count *(int)w);
}

XTextExtents(XFontStruct *f, const char *s, int count, int *dir,
	     int *ascent, int *descent, XCharStruct *overall)
{
int w;
f;s;
/*------eg if textwidth is 18 point = 39 in 0..4096,
	then it's 39 * 1280 / 4096 in 0..1280 */
w = (int)((long)textwidth * 1280L) / 4096L;
*ascent = w;			// (so far, use ascent + descent)
*descent = w/2;
return(1);
}

/*=============================================================================
**			DRAW -- note 0,0 is upper left
**===========================================================================*/
XDrawPoint(Display *d, Drawable w, GC gc, int x1,int y1)
{
struct WIN *p;
int x2,y2;
d;
p = win + ready_draw((Window)w, gc);
x2=x1; y2=y1;
move(x1,y1);
draw(x2,y2);		/* because pc_coord changes x,y */
return(1);
}

XDrawLine(Display *d, Drawable w, GC gc, int x1,int y1, int x2, int y2)
{
struct WIN *p;
d;
p = win + ready_draw((Window)w, gc);
move(x1,y1);
draw(x2,y2);
return(1);
}

XDrawLines(Display *d, Drawable w, GC gc, XPoint points[],
	   int npoints, int mode)
{
  int i,x,y;
  XPoint *p;
  d;
  ready_draw((Window)w, gc);
  for(i=0, p=&points[0]; i<npoints; i++,p++)
    {
      x = p->x;
      y = p->y;
      if (i==0)
        move(x, y);		/* (remember move() is a macro {} */
      if (i>0)
	draw(x, y);
    }
}

XDrawRectangle(Display *d, Drawable w, GC gc, int x, int y,
	unsigned int dx, unsigned int dy)
{
struct WIN *p;
int x2,y2;
p = win + ready_draw((Window)w,gc);
x2=x+dx;
y2=y+dy;
pc_coord(x,y);
pc_coord(x2,y2);
_moveto(x,y);
_lineto(x2,y);
_lineto(x2,y2);
_lineto(x,y2);
_lineto(x,y);
}

XFillRectangle(Display *d, Drawable w, GC gc, int x, int y,
	unsigned int dx, unsigned int dy)
{
struct WIN *p;
int x2,y2;
p = win + ready_draw((Window)w,gc);
x2=x+dx;
y2=y+dy;
pc_coord(x,y);
pc_coord(x2,y2);
_rectangle(_GFILLINTERIOR,x,y,x2,y2);
}

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

XDrawString(Display *d,Drawable w,GC gc,
	int x, int y, const char *s, int length)
{
struct WIN *p;
char text[100];

p = win + ready_draw((Window)w,gc);
if (current_textsize != (int)gc->values.font)
	ready_textsize(gc);
strncpy(text,s,length);
*(text+length)=0;
gtext(x,y,text);
}


/*=============================================================================
**				IMAGE MAP
**===========================================================================*/
XImage *XGetImage(Display *d, Drawable w, int x, int y,
		  unsigned int width, unsigned int height,
		  unsigned long plane_mask, int format)
{
  int x2,y2;
  XImage *q;
  char *p;
  long nbytes;
  
  if (format!=XYPixmap) return;
  x2 = x + width;
  y2 = y + height;
  select_view(getwin(w));
  pc_coord(x,y);
  pc_coord(x2,y2);
  nbytes = _imagesize(x,y,x2,y2);
  p = (char *)malloc(nbytes);
  q = (XImage *)malloc(sizeof(XImage));
  q->data = p;
  _getimage(x,y,x2,y2,p);
  return(q);
}

XPutImage(Display *d, Drawable w, GC gc, XImage *buf, int x0, int y0,
	  int x, int y, unsigned int width, unsigned int height)
{
  select_view(getwin(w));
  pc_coord(x,y);
  _putimage(x,y,buf->data,_GXOR);
}

getscreencoords(int x, int y, int *sx, int *sy)
{
  pc_coord(x,y);
  *sx = x;
  *sy = y;
}
getworldcoords(int x, int y, int *wx, int *wy)
{
  x = (int)((long)(x+xsmul/2)*xsdiv/xsmul);
  y = (int)((long)(y+ysmul/2)*ysdiv/ysmul);
  *wx = x;
  *wy = y;
}

/*=============================================================================
**				COLORS
**===========================================================================*/
Colormap XDefaultColormap(Display *d, int s)
{
return(0);
}

XAllocColor(Display *d, Colormap m, XColor *x)
{
  register int *p, *p2;
  int i, err, minerr, r,g,b;
  minerr=1024;
  for(p=colortable+offct,p2=p+64; p<p2; p+=4)
    {
      r = (x->red >> 8) & 0xff;
      g = (x->green >> 8 ) & 0xff;
      b = (x->blue >> 8) & 0xff;
      r -= *(p+1);	if (r<0) r=-r;
      g -= *(p+2);	if (g<0) g=-g;
      b -= *(p+3);	if (b<0) b=-b;
#ifdef DEAD_CODE
      err = (r>g) ? r: g;
      if (b>err) err=b;
#endif
      err = r + g + b;
      if (err<minerr)
	{
	  minerr=err;
	  i=*p;
	}
    }
  x->pixel = (long)i;
  return(1);
}

XStoreColor(Display *d, Colormap m, XColor *c)
{
int i,i0,j;
unsigned short *p;
static int newval[]={0,S,H,F};
d; m;
i0 = c->pixel*4 + offct;
colortable[i0]=c->pixel;
for(i=1,p=&c->red; i<4; i++,p++)
  {
    j=*p;
    colortable[i0 + i] = newval[j];
  }
colormap(4,colortable+i0);
}

/*=============================================================================
**				EVENTS, I/O
**===========================================================================*/
static int status=0, laststatus;

int curswindow()
{
int i,x1,y1,x2,y2;
struct WIN *p;
for(i=1; i<nwin; i++)
	{
	p=win+i;
	if (p->big) return(i);
	}
for(i=nwin-1; i>0; i--)
	{
	p = win+i;
	x1 = p->x;
	y1 = p->y;
	x2 = x1+p->dx;
	y2 = y1+p->dy;
	if (xcurs>=x1 && xcurs<=x2 && ycurs>=y1 && ycurs<=y2) break;
	}
return(i);
}

XPending(Display *d)
{
int i,n;
struct WIN *p;
for(i=n=0, p=win; i<nwin; i++,p++)
  if (p->expose) n++;
if (_bios_keybrd(_KEYBRD_READY)) n++;
return(n);
}

XNextEvent(Display *d, XEvent *x)
{
int i,x1,y1,x2,y2,xnew,ynew,dx,dy;
struct WIN *p;
char *q;

x->xany.display = &display;

for(;;)										// Wait for an event!
  {
    for(i=0; i<nwin; i++)		// seek first time exposure event
      {
	p=win+i;
	if (!p->expose) continue;
	x->xexpose.type = Expose;
	x->xexpose.window  = p->id;
	x->xexpose.count = 0;
	p->expose=0;
	return;
      }
    
    if (enternotify_win)
      {
	x->xcrossing.type = EnterNotify;
	x->xcrossing.window = enternotify_win;
	enternotify_win = 0;
	return;
      }

    if (_bios_keybrd(_KEYBRD_READY))	// Key hit
      {
	x->xkey.keycode = _bios_keybrd(_KEYBRD_READ);
	x->xkey.state = 0;
	x->type = KeyPress;
	i=curswindow();
	p=win+i;
	x->xkey.display=&display;
	x->xkey.window = p->id;
	x->xkey.x_root = xcurs;
	x->xkey.y_root = ycurs;
	x->xkey.x = xcurs - p->x;
	x->xkey.y = ycurs - p->y;
	return;
      }

    if (gin_enabled)				// enabled in XDefineCursor
      {
	laststatus=status;
	mouse(&dx,&dy,&status);
	status &= 1;				// ONLY button 1, for now
	if (dx || dy)
	  {
	    xnew=xcurs+dx;
	    ynew=ycurs+dy;
	    if (xnew<0) xnew=0;
	    if (ynew<0) ynew=0;
	    if (xnew>XSCREEN) xnew=XSCREEN;
	    if (ynew>YSCREEN) ynew=YSCREEN;
	    if (xnew!=xcurs || ynew!=ycurs)
	      {
		/*i = current_view*/;
		i=curswindow();
		p = win + i;
		select_view(0);
		drawcurs(current_curs);
		xcurs=xnew;
		ycurs=ynew;
		drawcurs(new_curs);
		current_curs = new_curs;
		select_view(i);
		x->xbutton.x_root = xcurs;
		x->xbutton.y_root = ycurs;
		x->xbutton.x = xcurs - p->x;
		x->xbutton.y = ycurs - p->y;
			
		if (eventmask & PointerMotionMask)
		  {
		    x->xmotion.type = MotionNotify;
		    x->xmotion.window = p->id;
		    return;
		  }
	      }
	  }
	else if (status != laststatus)
	  {
	    i=curswindow();
	    select_view(i);
	    p=win+i;
	    x->xbutton.type = status ? ButtonPress : ButtonRelease;
	    x->xbutton.display=&display;
	    x->xbutton.window = p->id;
	    x->xbutton.x_root = xcurs;
	    x->xbutton.y_root = ycurs;
	    x->xbutton.x = xcurs - p->x;
	    x->xbutton.y = ycurs - p->y;
	    //x->xbutton.state = 
	    x->xbutton.button = Button1;
	    return;
	  }
      }
  }
}

XRefreshKeyboardMapping(XMappingEvent *x)
{
}

XSelectInput(Display *d, Window w, long event_mask)
{
eventmask=event_mask;
}

XLookupString(XKeyEvent *x, char *buf, int nbytes,		// in Xutil.h
	KeySym *k, XComposeStatus *cs)
{
*buf = x->keycode & 0xff;
if (*buf) return(1);
*(buf+1) = x->keycode>>8;
return(2);
}


/*=============================================================================
**					CURSOR CONTROL
**===========================================================================*/
XWarpPointer(Display *d, Window src, Window dst, int srcx, int srcy,
			 unsigned int srcw, unsigned int srch, int x, int y)
{
struct WIN *p;
Drawcursor();
p = win + getwin(dst);
xcurs = p->x + x;
ycurs = p->y + y;
Drawcursor();
enternotify_win = dst;
}

Bool XQueryPointer(Display *d, Window w, Window *root, Window *child,
int *xr,int *yr, int *xw,int *yw, unsigned int *mask)
{
struct WIN *p;
int i;
i = curswindow();
p = win + i;
/*p=win+getwin(w);*/
*xw = *xr = xcurs;
*yw = *yr = ycurs;
*root  = win[0].id;
*child = p->id;
if (!p->big) { *xw -= p->x; *yw -= p->y; }
select_view(i);
}

int XDefineCursor(Display *d, Window w, Cursor curs)
{
int dx,dy;
struct WIN *p;
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
return(0);
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
for(i=0; i<ncurs; i++)
	if (cursid[i]==shape) return(cursid[i]);
cursid[ncurs++]=(Cursor)shape;
return(cursid[ncurs-1]);
}

/*-----------------------------------------------------------------------------
|	xPrintf, xPrintf1
-----------------------------------------------------------------------------*/
#define ROW1 21
#define COL1 dialog_col1
#define TCOLOR 2
#define BCOLOR 7
#define TCOLORX 3	/* 1 xor 2 */
#define BCOLORX 6	/* 1 xor 7 */

int dialog_col1 = 53;
void xPrintf(const char *s, int a, int b, int c, int d,
			    int e, int f, int g, int h)
{
static int row=1, col=1;
char text[256];
struct rccoord rc;

_settextwindow(ROW1,COL1, 34,79);	/* r1,c1, r2,c2 */
_wrapon(_GWRAPON);
_settextcolor(TCOLOR);

if (*s=='\b') { strcpy(text, " "); col--; }
else sprintf(text, s, a,b,c,d,e,f,g,h);
_settextposition(row, col);
_outtext(text);
if (*s=='\b') _settextposition(row, col);
rc = _gettextposition();
row = rc.row;
col = rc.col;
}

void xPrintf1(char *text, int start)
{
static int row=1, col=1;
int nchar;
struct rccoord rc;
char *p1, *p2, *p, s[80];

strcpy(s, text);
if ((p1 = strchr(s, '\n')) != NULL) *p1 = 0;

_settextwindow(ROW1,2, 34,COL1-2);	/* r1,c1, r2,c2 */
_wrapon(_GWRAPON);
xormode(1);
if (start) row = col = 1;
_settextposition(row, col);
for(p=s;;)
{
   p1 = strchr(p, '~');
   if (p1) p2 = strchr(p1+1, '~');
   if (p1 && p2) *p1 = *p2 = 0;
   else p1 = p2 = NULL;
   if (!p1 || p1>p)
   {
      _settextcolor(TCOLORX);
      _outtext(p);
   }
   if (p1)
   {
      _settextcolor(BCOLORX);
      _outtext(p1+1);
      *p1 = *p2 = '~';
      p = p2 + 1;
   }
   else break;
}
xormode(0);
rc = _gettextposition();
nchar = COL1 - 2 - rc.col;
for(p=s; p<s+nchar; *p++ = ' ') ;
*p = 0;
_outtext(s);
rc = _gettextposition();
row = rc.row;
col = rc.col;
}
