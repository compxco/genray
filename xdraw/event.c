/******************************************************************************
**  NAME	EVENT.C
**  AUTHOR	Sheryl M. Glasser
**
**  DESCRIPTION
**     Supporting functions for xwindows
**
**  Copyright (c) CounterPoint Graphics 1993.  All rights reserved.
******************************************************************************/
/*
Bugs, zoom + aspect
1.  givebox(): parameters set not necessarily for window currently being zoomed

Current problems, menu version:
1.  Zoom enabled: only receives MotionNotify after some other event,
    e.g. EnterNotify.  So crosshair doesn't move.  Yet the mask bit IS set.
2.  Doesn't respond to menu for ZZ or ZO
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include <fcntl.h>
#include <math.h>*/
#ifdef DOS
#include <graph.h>
#endif

#ifdef UNIX
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>
#define zprintf if (zdebug) printf

#else
#include <dos.h>
#include "xlib.h"
#include "xutil.h"
extern void Drawcursor(void);
extern void Drawcrosshair(int, int, int, int, int, int);
#define zprintf if (zdebug) xPrintf
#endif
extern int dialogwindow;
extern Window dialog_win;
int input_enabled = 0;

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
extern void postscript(char *);	/* near heap problem, try elim .h's*/

#ifdef USE_MENU
#include "device.h"
#include "menuwin.h"
#include "menu.h"		/* for set_selected_iwin */
#include "glx.h"
static int zdebug=1;
#else
static int zdebug=0;
#endif				/* ...of USE_MENU */

/*=============================================================================
**            DEFINITIONS AND VARIABLES
**===========================================================================*/
extern Display *mydisplay;
extern int myscreen;
extern Colormap cmap;
extern Window root;
extern Font font;
extern XFontStruct *font_struct;

static XEvent myevent;
static KeySym mykey;
static Window redraw_win, event_win=-1;
static GC redraw_gc, event_gc;

extern int nwindow;
extern VIEW view[];
extern int font_height;
extern unsigned long myforeground, mybackground;
extern CURVE_SET curveset[];
extern int exitflag, redrawflag, titledrawn;
extern unsigned int xscreen, yscreen;
int xhair_type = XHAIR_XOR;

/*=============================================================================
**                      VARIABLES FOR ZOOM
**===========================================================================*/
int zoom_on = 0, coord_on = 0;
float zoomfac = (float).05, xfromcoord = 0, yfromcoord=0;
static int zoom_count = 0, coord_count = 0, zoomi, coordi;
static int rubber = 1;
static Window zoomwin;
static GC zoomgc;
static float fxbox1, fybox1, fdxbox, fdybox, zoomdx, zoomdy;
static int xcurs, ycurs;
static int cursor_defined = 0;
static Cursor curs1, curs2;

#ifdef UNIX
#define CURS1 XC_draft_large
#define CURS2 XC_top_left_arrow
/*#define CURS1 XC_top_left_corner*/
/*#define CURS2 XC_bottom_right_corner*/
#else
#define CURS1 8
#define CURS2 9
#endif

#ifdef USE_MENU
#define PRESSMASK (ButtonPressMask | ButtonReleaseMask | KeyPressMask)
#define WHEREMASK (LeaveWindowMask | EnterWindowMask)

#else			/* ...of USE_MENU */
#define PRESSMASK (ButtonPressMask | KeyPressMask)
#define WHEREMASK LeaveWindowMask
#endif			/* ...of USE_MENU */

#define MOTIONMASK (PointerMotionMask | PointerMotionHintMask)
#define EVENTMASK3 (PRESSMASK | WHEREMASK | ExposureMask)
#define EVENTMASK4 (EVENTMASK3 | MOTIONMASK)

/*=============================================================================
**		FUNCTIONS FOR MENU VERSION
**===========================================================================*/
#ifndef USE_MENU
void set_menu_visibility(int i) { i; }
void menu_setwindows(int i) { i; }

#else			/* .....of #ifndef USE_MENU, so USE_MENU is defined */
/*-----------------------------------------------------------------------------
|	init_menuXlibrary
-----------------------------------------------------------------------------*/
void init_menuXlibrary()
{
  unsigned long back;
  back = glx_index(get_menubkgd());
  init_glx(mydisplay, root, myforeground, back,
	   font, font_struct, &myevent, xscreen, yscreen);
}
#endif			/* ......of USE_MENU undefined, defined */


/*-----------------------------------------------------------------------------
|	init_menus, clear_if_needed: non-Motif
-----------------------------------------------------------------------------*/
#ifndef MOTIF
#ifndef USE_MENU
  void init_menus() {}
#endif

void clear_if_needed(Window w, GC gc, unsigned int dx, unsigned int dy)
		    {w;gc;dx;dy;}

/*-----------------------------------------------------------------------------
|	init_menus, clear_if_needed: Motif
-----------------------------------------------------------------------------*/
#else			/* .....MOTIF */
void init_menus()
{
  int x, y;
  nextwindow(&x, &y);
  execute_motif(x, y, nwindow, mydisplay);
}

void clear_if_needed(Window w, GC gc, unsigned int dx, unsigned int dy)
{
  xrectangle(dx,dy,BlackPixel(mydisplay, myscreen),
	     WhitePixel(mydisplay, myscreen));
}

void test_window(char *s, Window win, Widget wid)
{
  unsigned int dx, dy, b, d, nchildren;
  int x0, y0;
  Window root, parent, *children;
  Position x,y;
  Dimension w, h;

  XGetGeometry(mydisplay, win, &root, &x0, &y0, &dx, &dy, &b, &d);
  printf ("\n%s: Window %lx, root %lx has x,y = %d %d, size %d %d\n",
	  s, win, root, x0, y0, dx, dy);
  if (wid != 0)
    {
      XtVaGetValues(wid, XmNx, &x, XmNy, &y, XmNwidth, &w, XmNheight, &h, NULL);
      printf ("%s: Widget %lx, Window %lx has x,y = %d %d, size %d %d\n",
	      s, wid, win, x, y, w, h);
    }
  XQueryTree(mydisplay, win, &root, &parent, &children, &nchildren);
  /*printf("%s: window=%lx, root=%lx, parent=%lx, nchildren=%d\n",
	 s, win, root, parent, nchildren);*/
  XGetGeometry(mydisplay, parent, &root, &x0, &y0, &dx, &dy, &b, &d);
  printf ("%s: Window %lx has x,y = %d %d, size %d %d\n",
	  s, parent, x0, y0, dx, dy);
  XFree(children);
}

/*-----------------------------------------------------------------------------
|	set_expose: for calls to get_expose_info
-----------------------------------------------------------------------------*/
Window set_expose(int i)
{
  myevent.xexpose.display = mydisplay;
  myevent.xexpose.window = view[i].window;
  return view[i].window;
}
#endif			/* ......of MOTIF undefined, defined */


/*=============================================================================
**                        EVENTS
**===========================================================================*/

/*-----------------------------------------------------------------------------
|	get_expose_info
-----------------------------------------------------------------------------*/
void get_expose_info(Display **d, Window *w, GC *gc, int *fonthp)
{
  int iwindow;

  if (NOX) return;
  mydisplay = *d = myevent.xexpose.display;
  redraw_win = *w = myevent.xexpose.window;
  iwindow = getview(redraw_win);	/* get which window is redrawing */
  if (iwindow < nwindow)
    redraw_gc = *gc = view[iwindow].gc;
  *fonthp = font_height;
}

/*-----------------------------------------------------------------------------
|	
-----------------------------------------------------------------------------*/
void set_winmask(Window win)  
{
#ifndef MOTIF
  XSelectInput(mydisplay, win, EVENTMASK3);
#endif
}
/*-----------------------------------------------------------------------------
|	readkey, parsekey
-----------------------------------------------------------------------------*/
int parsekey()
{
  int i, ic;
  char text[20];
  text[0] = text[1] = 0;
  i = XLookupString(&myevent.xkey, text, 10, &mykey, 0);
  if (i < 1) ic=0;
  /*else if (i==0) ic = *(int *)&mykey; (extra line from xtools.c--ok? */
  else ic = *(int *)text;
  return(ic);
}

int readkey()
{
  zprintf("In readkey()\n");
  for(;;)
    {
      XNextEvent(mydisplay, &myevent);	/* read the next event */
      if (myevent.type == KeyPress) break;
    }
  return(parsekey());
}

/*-----------------------------------------------------------------------------
|   xmessage -- display or erase a message in the current window during redraw
-----------------------------------------------------------------------------*/
void xmessage(int x, int y, char *s)
{
  if (NOX)
    return;
  set_xor(redraw_gc, 1);
  XSetForeground(mydisplay, redraw_gc, myforeground);
  XDrawString(mydisplay, redraw_win, redraw_gc, x, y, s, strlen(s));
  set_xor(redraw_gc, 0);
}

void xrectangle(unsigned int dx, unsigned int dy, int bkgd, int fore)
{
  XSetForeground(mydisplay, redraw_gc, BlackPixel(mydisplay, myscreen));
  XFillRectangle(mydisplay, redraw_win, redraw_gc, 0, 0, dx, dy);
  XSetForeground(mydisplay, redraw_gc, WhitePixel(mydisplay, myscreen));
}

/*-----------------------------------------------------------------------------
|	print_event
-----------------------------------------------------------------------------*/
void print_event(int etype, int iwin, long event_win)
  {
#ifdef UNIX
    static int inited=0;
    XWindowAttributes att;
    if (!inited) { printf("Display=%lx\n",mydisplay); inited=1; }
    printf("------- ");
    if      (etype==MotionNotify)  printf("Motion ");
    else if (etype==Expose)        printf("Expose ");
    else if (etype==ButtonPress)   printf("Press  ");
    else if (etype==ButtonRelease) printf("Release");
    else if (etype==EnterNotify)   printf("Enter  ");
    else if (etype==LeaveNotify)   printf("Leave  ");
    else if (etype==KeyPress)      printf("Key");
    else                           printf("Unknown");

    XGetWindowAttributes(mydisplay, event_win, &att);
    printf(" in window %d, ID %lx, masks %lx,%lx\n",
	   iwin, event_win, att.your_event_mask, att.all_event_masks);
    fflush(stdout);
#else
    etype; iwin; event_win;
#endif
  }

/*-----------------------------------------------------------------------------
|	keystrokes
-----------------------------------------------------------------------------*/
static void keystrokes(void);
static void keystrokes()
{
  xprintf("***** Keystroke Summary *********************************\n");
  xprintf("|  q ......quit\n");
  xprintf("|  p ......print window (create .ps file)\n");
  xprintf("|  t ......display descriptive title of window\n");
  xprintf("|  k ......this message\n");
  xprintf("|  \n");
  xprintf("|  m ......toggle markers: on/off\n");
  xprintf("|  a ......toggle aspect ratio: preserved/auto-scale\n");
  xprintf("|  e ......toggle Extrema, y from:  all curves/visible curves/-E values\n");
  xprintf("|  f ......toggle fill: on/off\n");
  xprintf("|  l ......toggle automatic labels on curves (random positions)\n");
  xprintf("|  1 ......draw one curve of family (first,next,...)\n");
  xprintf("|  2 ......draw a second curve of family (first,next,...)\n");
  xprintf("|  0 ......draw all curves of family\n");
  xprintf("|  # ......prompt for one curve of family\n");
  xprintf("|  >,<.....draw next, previous step of outer loop\n");
  xprintf("|  @ ......edit mode (partially implemented)\n");
  xprintf("|  \n");
  xprintf("|  z ......zoom (use crosshairs)\n");
  xprintf("|  o ......original, unzoomed size\n");
  xprintf("|  +,- ....zoom up or down slightly\n");
  xprintf("|  \n");
  xprintf("|  c ......get-coordinates mode on/off (use crosshairs)\n");
  xprintf("|  s ......get-slope mode on/off (use crosshairs)\n");
  xprintf("|  r ......get-ratio mode on/off (use crosshairs)\n");
  xprintf("|  \n");
  xprintf("|  d,h ....double, halve number of lines in contour plot\n");
  xprintf("|  v ......print value of contour lines\n");
  xprintf("********************************************************\n");
}

#ifndef USE_MENU

#ifndef MOTIF
/*-----------------------------------------------------------------------------
|	abort_zoom
-----------------------------------------------------------------------------*/
static int abort_zoom(void);
static int abort_zoom()
{
  int retval;
  retval = 1;
  if      (zoom_on  && event_win != zoomwin) zoom(zoomwin, ZABORT);
  else if (coord_on && event_win != zoomwin) coord(coord_on);
  else retval = 0;
  return retval;
}

/*-----------------------------------------------------------------------------
|	event -- for keystroke version
|	* wait for an event, process it, return
-----------------------------------------------------------------------------*/
void event()
{
  char c;
  int iwin, ic, dummy, n, result;
  int xw, yw, xr, yr, x;
  unsigned int keys_buttons;
  Window rw, cw, oldw;
  static char editing=0;
  static char inited=0;
  CURVE_SET *cp;
  char text[200];

  if (NOX) exit(0);
  XNextEvent(mydisplay, &myevent);	/* read the next event */
  oldw = event_win;
  event_win = myevent.xany.window;
  /*print_event(myevent.type, getview(event_win), event_win);*/
  
  switch (myevent.type)
    {
    case Expose:		/* repaint window on expose events */
      if (myevent.xexpose.count == 0)	/* draw the curve */
	{
#ifdef UNIX
	  if (event_win == dialog_win) redraw_dialog();
	  else redraw();
#else
	  Drawcursor();
	  redraw();
	  Drawcursor();
#endif
	  if (!inited && (dialogwindow==0 || event_win==dialog_win))
	      { xprintf("Press <k> for Keystroke summary.\n"); inited=1; }

	}
      break;

    case MappingNotify:			/* process keyboard mapping changes */
      XRefreshKeyboardMapping(&myevent.xmapping);
      break;

    case ButtonPress:			/* process mouse-button presses */
      if (abort_zoom()) ;		/* click in another win: disable */
      else if (zoom_on) zoom(zoomwin, ZGET);
      else if (coord_on) coord(8);
      break;

    case LeaveNotify:
      if (zoom_on) zoom(zoomwin, ZIGNORE);	/* Gives message for debug */
      break;

    case MotionNotify:
      if ((zoom_on || coord_on) && event_win == zoomwin)
      {
        XQueryPointer(mydisplay, event_win,
		    &rw, &cw, &xr, &yr, &xw, &yw, &keys_buttons);
        crosshair(xcurs, ycurs);
        xcurs = xw;
        ycurs = yw;
        crosshair(xw, yw);
      }
      break;

    case KeyPress:		/* process keyboard input */
      ic = parsekey();
      c = *(char *)&ic;
      n = get_graph(event_win, 0);
      cp = curveset + n;
      iwin = getview(event_win);
      redrawflag = 0;
      
      abort_zoom();		/* key in another win: disable */

/*-------- Output commands */
      if (input_enabled)
      {
	 *text = c; *(text+1) = 0;
	 if (xinput(text, &ic) == 1)
	 {
	    cp->flags |= SINGLE;
	    cp->ncurve = ic - 1;
	    cp->hzcurve = -1;
	    toggle_single(cp, '1');
	 }
      }
      else if (c == '#') xinput("Enter desired index: ", NULL);
      else if (c == 'q') exitflag = 1;
      else if (c == 'p' && iwin < nwindow)
	{
	  parse_title(cp, cp->title, text);
	  if (cp->subtitle)
	    {
	      strcat(text, "; ");
	      parse_subtitle(cp, text+strlen(text));
	    }
	  postscript(text);
	}
      else if (c=='k') keystrokes();
      else if (c=='t') print_title(iwin, cp);

/*--------- Zoom commands */
      else if (zoom_on  && event_win==zoomwin)
      {
        if (c == 'z') zoom(zoomwin, ZABORT);
      }
      else if (coord_on && event_win==zoomwin &&
	       (c!='c' && c!='s' && c!='r'));

      else if (c == 'z') zoom(event_win, ZENABLE);
      else if (c == 'o') zoom(event_win, ZORIGINAL);
      else if (c == 'c') coord(1);
      else if (c == 's') coord(2);
      else if (c == 'r') coord(3);
      else if (c=='+' || c=='-') addzoom(event_win, c=='+');

/*----------- Display commands */
      else if (c=='l') toggle_flabel(cp);
      else if (c=='d' || c=='h') new_ncurve(cp, c);
      else if (c=='v') contour_values(cp);
      else if (c == 'm') toggle_markers(cp, &dummy);
      else if (c == 'a') toggle_aspect(cp, &dummy);
      else if (c == 'f') toggle_fill(cp);
      else if (c == 'e') toggle_extrema_use(cp);

/*--------- Editor commands */
      else if (c=='1' || c=='2' || c=='0') toggle_single(cp,c);
      else if (c=='>' || c=='<') toggle_single(cp,c);
      else if (c == '@')
	{
	  if (editcurve(iwin,cp,NULL))
	    XClearArea(mydisplay, event_win, 0, 0, 0, 0, True);
	}


#ifndef UNIX
      else if (c==0 && ic!=0 && iwin < nwindow)
	{
	  c = *((char *)&ic + 1);
	  if	  (c == 0x44)    x = 0;		/* F10 */
	  else if (c == 0x43)    x =-1;		/* F9 */
	  else c = 0;
	  if (c)
	    XMoveResizeWindow(mydisplay, event_win, x, x, 0, 0);
          if (c && x==-1) titledrawn=0;
	}
#endif
      if (redrawflag)
	XClearArea(mydisplay, event_win, 0, 0, 0, 0, True);
      break;
    }				/* switch (myevent.type) */
}				/* while (done==0) */

#else			/* ...of not MOTIF */
/*-----------------------------------------------------------------------------
|	event -- for MOTIF version
|	* All action in execute_motif(), called from init_menus()
-----------------------------------------------------------------------------*/
void event() {  exitflag = 1; }

#endif			/* ...of not MOTIF */

#else			/* ...USE_MENU is defined */

/*-----------------------------------------------------------------------------
|	
-----------------------------------------------------------------------------*/
int get_eventmask()
{
return EVENTMASK3;
}

Window get_op_window(Window win)
{
  return(win_is_menu(event_win) ? win : event_win);
}

/*-----------------------------------------------------------------------------
|	enable_window
|	* Preparation for Zoom, Coord, Slope
|	* If cursor is in menu, need to warp to selected window.
|	  This will first generate a LeaveNotify, EnterNotify event
|	* THEN enable appropriate function
-----------------------------------------------------------------------------*/
Window enable_window(Window win)
{
  Window root;
  int x0, y0, xw, yw;
  unsigned int dx, dy, bwidth, depth;
  if (win != event_win)
  {
    XGetGeometry(mydisplay, win, &root,
	  &x0, &y0, &dx, &dy, &bwidth, &depth);
    xw = dx / 2;
    yw = dy / 2;
    XWarpPointer(mydisplay, None, win, 0,0,0,0, xw,yw);
    do {
      XNextEvent(mydisplay, &myevent);
    }
    while (myevent.type != EnterNotify);
    event_win = myevent.xany.window;
  }
  return(event_win);
}

/*-----------------------------------------------------------------------------
|	event - for menu version
|		call when enable expose for each window, then call repeatedly
|	Verify: QueryPointer returns same as .xbutton.x etc.
|		value of .button
|		only gives events if in a valid window
|		ButtonPress only when it's a change
-----------------------------------------------------------------------------*/
#ifdef UNIX
#define redraw_exposed() redraw()
#else
#define redraw_exposed() { Drawcursor(); redraw(); Drawcursor(); }
#endif

void event()
{
  char c;
  int iwin, i, ic, etype, enter;
  Window oldw;
  short int state;
  int button;

  if (NOX) exit(0);
  zprintf("waiting.."); fflush(stdout);
  XNextEvent(mydisplay, &myevent);	/* read the next event */
  oldw = event_win;
  event_win = myevent.xany.window;
  etype = myevent.type;
  iwin = getview(event_win);		/* index to window (NOT curveset) */
					/* WARN! could be menu window! */
  print_event(etype,iwin,event_win);

#ifndef UNIX
  if (etype != Expose && oldw != event_win && win_is_menu(event_win))
    menu_ievent(event_win);
#endif

  if (etype==Expose && myevent.xexpose.count==0)
    {
      if (win_is_menu(event_win)) redrawmenu(NULL, 0);
      else redraw_exposed();	/* getproperties() sets redraw_win */
    }

  else if (etype==MappingNotify)	/* process keyboard mapping changes */
    XRefreshKeyboardMapping(&myevent.xmapping);

  else if (etype==ButtonPress || etype==ButtonRelease)
    {
      button = myevent.xbutton.button;
      if (button==Button1 || button!=Button2)
	button = LEFTMOUSE;
      else
	button = RIGHTMOUSE;
      state = (etype==ButtonPress) ? DOWN : UP;
      if (menu_event(button, state)) ;
      else if (state != DOWN) ;
      else if (zoom_on)
        {
	  if (button==LEFTMOUSE && event_win==zoomwin)	/* WARN! can be Abort */
	    zoom(zoomwin, ZGET);
	  else zoom(zoomwin, ZDISABLE);
	}
      else if (button==LEFTMOUSE) set_selected_iwin(iwin);
    }

  else if (etype==EnterNotify || etype==LeaveNotify)
    {
      /*menu_event(INPUTCHANGE, iwin);*/
      enter = (etype==EnterNotify);
      menu_ievent(enter ? event_win : 0L);
    }

  else if (etype==MotionNotify && zoom_on && event_win==zoomwin)
    {
      crosshair(xcurs, ycurs);
      xcurs = myevent.xbutton.x;
      ycurs = myevent.xbutton.y;
      crosshair(xcurs, ycurs);
    }

  else if (etype==KeyPress)		/* process keyboard input */
    {
      ic = parsekey();
      c = *(char *)&ic;
#ifndef UNIX
      if (c==0 && ic!=0 && iwin < nwindow)
	{
	  c = *((char *)&ic + 1);
	  if	  (c == 0x44) i = 0;		/* F10 */
	  else if (c == 0x43) i = -1;		/* F9 */
	  else	c = 0;
	  if (c) XMoveResizeWindow(mydisplay, event_win, i, i, 0, 0);
	}
      else
#endif
      menu_event(KEYBD, (short int)c);
    }
}

/*-----------------------------------------------------------------------------
|	clear_redraw
-----------------------------------------------------------------------------*/
void clear_redraw(int iwin)
{
  XClearArea(mydisplay, view[iwin].window, 0, 0, 0, 0, True);
}

#endif			/* ...of USE_MENU */

/*=============================================================================
**				ZOOM AND COORD
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	zoom
|	* 1=ZENABLE, 2=ZGET=process button hit while enabled, 3=ZABORT
|	* 0=ZDISABLE, -1=ZORIGINAL, -2=ZIGNORE=Leave event during zoom
-----------------------------------------------------------------------------*/
void zoom(Window nextwin, int enable)
{
  /*Window nextwin;*/
  int i, xw, yw, xr, yr;
  unsigned keys_buttons;
  Window rw, cw;
  VIEW *v;

  zprintf("Zoom(%d): zoomwin=%lx, myevent.xmotion.win=%lx\n",
	  enable, zoomwin, nextwin /*myevent.xmotion.window*/);
  if (enable==ZIGNORE) return;		/* ignore LeaveEvent */

  /*nextwin = myevent.xkey.window;*/
  i = getview(nextwin);
  if (i == nwindow) return;
  v = view + i;

  XQueryPointer(mydisplay, myevent.xmotion.window,
		&rw, &cw, &xr, &yr, &xw, &yw, &keys_buttons);

  if (enable==ZDISABLE || enable==ZORIGINAL || enable==ZABORT)
    {
      if (zoom_on && zoomwin != 0)
	{
	  XUndefineCursor(mydisplay, zoomwin);
	  XSelectInput(mydisplay, zoomwin, EVENTMASK3);
	}
      else if (zoom_on && zoomwin == 0)
	zprintf("Wierd! zoom_on TRUE, zoomwin==0!\n");
      if (enable == ZORIGINAL)
      {
        v->clipped = 0;
        XClearArea(mydisplay, nextwin, 0, 0, 0, 0, True);
      }
      else if (enable==ZABORT)
      {
	/*newclip(-1, xw, yw);*/
        if (zoom_count > 0)
	  XClearArea(mydisplay, nextwin, 0, 0, 0, 0, True);
	else crosshair(xcurs, ycurs);
      }
      zoom_on = zoom_count = 0;
    }

  else if (enable == ZGET)	/* Process button hit while enabled */
    {
      if (nextwin != zoomwin) return;
      ztest(&xw, &yw);	/* WARN! */
      xcurs = xw;
      ycurs = yw;
      if (zoom_count++ == 0)
	{
	  XDefineCursor(mydisplay, zoomwin, curs2);
	  if (xhair_type==XHAIR_CURSOR) xhair_type = XHAIR_1ST_LINE;
	  newclip(0, xcurs, ycurs);
	  crosshair(-1,0);
	  crosshair(xw, yw);	/* "erase" leaves crosshair on move */
	  if (xhair_type==XHAIR_1ST_LINE) xhair_type = XHAIR_CURSOR;
	}
      else
	{
	  zoom(zoomwin,ZDISABLE);	/* (no need to erase crosshair) */
	  newclip(1, xcurs, ycurs);
	  XClearArea(mydisplay, zoomwin, 0, 0, 0, 0, True);
	}
    }

  else				/* Enable crosshair */
    {
      if (!cursor_defined) define_cursor(nextwin);
      zoom(zoomwin,ZDISABLE);
      zoomwin = nextwin;
      zoomi = getview(zoomwin);
      zoomgc = view[zoomi].gc;

      XDefineCursor(mydisplay, zoomwin, curs1);
      XSelectInput(mydisplay, zoomwin, EVENTMASK4);
      crosshair(-1, 0);	/* Initialize */
      crosshair(xcurs = xw, ycurs = yw);
      zoom_on = 1;
    }
}

/*-----------------------------------------------------------------------------
|	define_cursor
-----------------------------------------------------------------------------*/
void define_cursor(Window nextwin)
{
   static Pixmap crosspix=0;
#define XX 0x80
   static unsigned char cross[] = {
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0xff,0xff,0xff,0xff,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0
    };
#ifdef UNIX
   XColor fore, bkgd;
   unsigned bestdx, bestdy;

   if (xhair_type==XHAIR_CURSOR)
   {
      XQueryBestCursor(mydisplay, nextwin, xscreen, yscreen,
		       &bestdx, &bestdy);
      /*printf("Best cursor size %d %d\n", bestdx, bestdy);*/
      if (bestdx > 32) bestdx = 32;
      crosspix = XCreatePixmapFromBitmapData(mydisplay, nextwin, (char *)cross,
					     31, 31, 1, 0, 1);
      bkgd.pixel = mybackground;
      fore.pixel = myforeground;
      XQueryColor(mydisplay, cmap, &bkgd);
      XQueryColor(mydisplay, cmap, &fore);
      curs1 = curs2 = XCreatePixmapCursor(mydisplay, crosspix, crosspix,
					  &fore, &bkgd, 15, 15);
    }
   else
#endif
   {
      curs1 = XCreateFontCursor(mydisplay, CURS1);
      curs2 = XCreateFontCursor(mydisplay, CURS2);
   }
   cursor_defined = 1;
   nextwin;
}

/*-----------------------------------------------------------------------------
|	tell_zoom -- prepare to call functions from motif
-----------------------------------------------------------------------------*/
void tell_zoom(Window zwin, int which, int flag)
{
   zoomi = which;
   zoomwin = view[which].window;
   zoomgc  = view[which].gc;
   coord_on = flag;
   /*myevent.xany.window = zwin;*/	/* should be same as zoomwin */
   zwin;
   myevent.xany.window = zoomwin;	/* WARN! zwin=0 from xmotif */
}

void tell_key(XEvent *p)
{
  myevent.xany.window = p->xany.window;
  myevent.xkey = p->xkey;
}

/*-----------------------------------------------------------------------------
|   newclip - called by zoom() to get current clip
|   * input xcurs,ycurs is relative to window, i.e. 0..width
|   * i: 0=1st point, 1=2nd point, -1=original
|   * gx, gy = fraction (0.0 ... 1.0) within inner box
-----------------------------------------------------------------------------*/
void newclip(int i, int xcurs, int ycurs)
{
  Window root;
  int x0, y0;
  VIEW *v;
  unsigned int bwidth, depth, zdx, zdy;
  float gx, gy, t;
  static float gx0, gy0;

#define fswap(i,j) { t = v->f[i]; v->f[i]=v->f[j]; v->f[j]=t; }

  v = view + zoomi;
  if (i ==-1) { v->clipped = 0; return; }

  XGetGeometry(mydisplay, zoomwin, &root,	/* Geometry can change! */
	       &x0, &y0, &zdx, &zdy, &bwidth, &depth);
  zoomdx = (float)zdx;
  zoomdy = (float)zdy;

  /*----- get fractional position of cursor relative to unclipped -----*/
  gx = (float) (xcurs - fxbox1 * zoomdx) / (fdxbox * zoomdx);
  gy = (float) (ycurs - fybox1 * zoomdy) / (fdybox * zoomdy);
  /*if (i==0) printf("\n");
    printf("Cursor fractionally at %f,%f in box", gx,gy);*/
  
  if (v->clipped)
    {
      gx = v->f[0] + gx * (v->f[2] - v->f[0]);
      gy = v->f[3] + gy * (v->f[1] - v->f[3]);
      /*printf("-->%f,%f", gx,gy);*/
    }
  /*printf("\n");*/

  if (i == 0)
    {
      gx0 = gx;
      gy0 = gy;
    }
  else
    {
      v->f[0] = gx0;
      v->f[1] = gy0;
      v->f[2] = gx;
      v->f[3] = gy;
      if (v->f[0] > v->f[2]) fswap(0, 2);
      if (v->f[1] < v->f[3]) fswap(1, 3);
      v->clipped = 1;
    }
}

/*-----------------------------------------------------------------------------
|   givebox -- called by redraw() in preparation for possible zoom
-----------------------------------------------------------------------------*/
void givebox(int x0, int y0, int w, int h, int ww, int wh)
{
  zoomdx = (float)ww;
  zoomdy = (float)wh;
  fxbox1 = (float) x0 / zoomdx;
  fybox1 = (float) y0 / zoomdy;
  fdxbox = (float) w / zoomdx;
  fdybox = (float) h / zoomdy;
}

/*-----------------------------------------------------------------------------
|	set_aspected_limits
|	* IN:  xmin1, xmax1 = autoscale data limits, xmin, xmax = aspected
|	* which: 0=x, 1=y
-----------------------------------------------------------------------------*/
void set_aspected_limits(Window win, int which, float xmin1, float xmax1,
			 float xmin, float xmax)
{
  VIEW *v;
  float min, max, temp;
  v = view + getview(win);
  if (! (v->clipped & 2))
    {
      min = (xmin - xmin1) / (xmax1 - xmin1);
      max = (xmax - xmin1) / (xmax1 - xmin1);
      if (which==1) { temp = min; min = max; max = temp; }
      v->f[0+which] = min;
      v->f[2+which] = max;
      if (which==1) v->clipped |= 2;
#ifdef DEAD_CODE
      printf("Aspect%s:\
 %f ... %f --> %f ... %f ==> Fraction = %f ... %f\n",
	     (which==0)?"x":"y", xmin1,xmax1,xmin,xmax,
	     v->f[0+which], v->f[2+which]);
      if (which==0) printf("f[] before: ");
      else printf("f[] after: ");
      printf("%f ... %f\n", v->f[1], v->f[3]);
#endif
    }
}

/*-----------------------------------------------------------------------------
|   get_zoomedclip - called by redraw() to get current clip
|   e.g. xmin --> xmin + fx1 * (xmax-xmin)
-----------------------------------------------------------------------------*/
int get_zoomedclip(Window w, float *fx1, float *fy1, float *fx2, float *fy2)
{
  VIEW *v;
  int i;

  i = getview(w);
  v = view + i;
  if (i == nwindow || !v->clipped)
    return (0);
  if (fx1 != NULL)
    {
      *fx1 = v->f[0];
      *fx2 = v->f[2];
      *fy1 = v->f[1];
      *fy2 = v->f[3];
    }
  return (v->clipped);
}

/*-----------------------------------------------------------------------------
|	addzoom
|	* inc = TRUE for bigger image ==> smaller surrounding box
-----------------------------------------------------------------------------*/
void addzoom(Window w, int inc)
{
  float df,diff,avg;
  VIEW *v;
  df = inc ? -zoomfac : zoomfac;
  df = (float) 1 + df / (float)2;
  v = view + getview(w);
  if (!v->clipped)
    {
       v->f[0] = v->f[3] = (float)0;
       v->f[2] = v->f[1] = (float)1;
       v->clipped = 1;
    }
  avg =  (v->f[0] + v->f[2]) / (float)2;
  diff = (v->f[2] - v->f[0]) / (float)2;
  v->f[2] = avg + df * diff;
  v->f[0] = avg - df * diff;

  avg =	 (v->f[1] + v->f[3]) / (float)2;
  diff = (v->f[1] - v->f[3]) / (float)2;
  v->f[1] = avg + df * diff;
  v->f[3] = avg - df * diff;
  redrawflag=1;
}

/*-----------------------------------------------------------------------------
|	show_coord
|	* coord_on: 1=coordinate 2=slope, 3=ratio
-----------------------------------------------------------------------------*/
void show_coord(int xcurs, int ycurs, int count)
{
  float xmin, ymin, xmax, ymax, dx,dy, x,y;
  static float x1,y1;
  XRectangle bigbox, clipbox;
  static char *point[] = { "point on ", "y-value for " };
  static char *slope[] = { "slope", "ratio" };
  static char *which[] = { "1st ", "2nd "};
  char text[80], *pt;
  int i, j, n;
  
  get_box_info(zoomwin, &bigbox, &clipbox, &xmin, &ymin, &xmax, &ymax);
  dx = (float)(xcurs - clipbox.x) / (float)clipbox.width;
  dy = (float)(ycurs - clipbox.y) / (float)clipbox.height;
  x = xmin + dx * (xmax - xmin);
  y = ymin + dy * (ymax - ymin);

  i = coord_on - 2;		/* 0: slope, 1: ratio */
  j = (count==1) ? 0 : 1;
  *text = 0;
  if (coord_on==1)
    sprintf(text, "Coordinate = ");
  else if (coord_on==2 || coord_on==3)
    sprintf(text, "%s%s%s = ", which[j], point[i], slope[i]);
  pt = text + strlen(text);
  if (coord_on != 3) sprintf(pt, "%.4g, %.4g", x, y);
  else sprintf(pt, "%.4g", y);
  xfromcoord = x;
  yfromcoord = y;
  strcat(text, "\n");
  xprintf(text);

  if (coord_on == 1) ;			/* 'c' mode */
  else if (count==2)			/* 's' or 'r' mode */
  {
    n = strlen(which[j]) + strlen(point[i]);
    pt = text;
    for (j=0,*pt++='\n'; j<n; j++,*pt++=' ');
    sprintf(pt, "%s = ", slope[i]);
    pt += strlen(pt);
    if (coord_on==2)
    {
      dx = x - x1;
      dy = y - y1;
      if (dx == FZ) sprintf(pt, "infinity");
      else sprintf(pt, "%g", dy / dx);
    }
    else if (coord_on==3)
    {
      if (y==FZ) sprintf(pt, "infinity");
      else sprintf(pt, "%g", (float)y1 / (float)y);
    }
    strcat(text, "\n");
    xprintf(text);
#ifndef MOTIF
    if (coord_on==2 || coord_on==3)	/* disable 's' or 'r', redraw win */
      coord(coord_on);
#endif
  }
  x1 = x;
  y1 = y;
}

/*-----------------------------------------------------------------------------
|	coord
|	* 1=enable coord, 2=enable slope, 3=enable ratio, 8=process hit
-----------------------------------------------------------------------------*/
void coord(int enable)
{
  Window nextwin;
  static int cursor_defined = 0;
  int i, old_coord_on;
  int xw, yw, xr, yr;
  unsigned keys_buttons;
  Window rw, cw;
  VIEW *v;

  nextwin = myevent.xkey.window;
  i = getview(nextwin);
  if (i == nwindow) return;
  v = view + i;

  if (enable != 8 && coord_on)		/* it's a toggle, on-->off */
    {
      XUndefineCursor(mydisplay, zoomwin);
      XSelectInput(mydisplay, zoomwin, EVENTMASK3);
      old_coord_on = coord_on;
      coord_on = coord_count = 0;
      /*XClearArea(mydisplay, nextwin, 0, 0, 0, 0, True);*/
      crosshair(xcurs, ycurs);
      if (old_coord_on == enable) return;
    }

  XQueryPointer(mydisplay, myevent.xmotion.window,
		&rw, &cw, &xr, &yr, &xw, &yw, &keys_buttons);
  xcurs = xw;
  ycurs = yw;

  if (enable == 8)		/* Process button hit while enabled */
    {
      if (nextwin != zoomwin) return;
      if (coord_count >= 2) coord_count = 0;
      coord_count++;
      show_coord(xcurs,ycurs, coord_count);
    }

  else				/* Enable! (off-->on) */
  {
    zoomwin = nextwin;
    coordi = getview(zoomwin);
    zoomgc = view[coordi].gc;
    coord_on = enable;
    if (!cursor_defined) define_cursor(nextwin);
    XDefineCursor(mydisplay, zoomwin, curs1);
    XSelectInput(mydisplay, zoomwin, EVENTMASK4);
    crosshair(-1,0);	/* Initialize */
    crosshair(xcurs, ycurs);
  }
}

/*-----------------------------------------------------------------------------
|	crosshair
|	* normally called in pairs: erase the old, draw the new
|	* crosshair types: 1=normal, all others for cases where no xor
|	*   0=save under, draw line; 3=user cursor
-----------------------------------------------------------------------------*/
void crosshair(int x, int y)
{
  int x1, y1, x2, y2;
  unsigned int dx, dy;
  static int state = -1;
  static XImage *xih=0, *xiv=0;
  static XImage *xiH=0, *xiV=0;

  if (!rubber) return;
  if (x==-1) { state = 0; return; }

  x1 = (int)(fxbox1 * zoomdx);
  y1 = (int)(fybox1 * zoomdy);
  x2 = (int)((fxbox1 + fdxbox) * zoomdx);
  y2 = (int)((fybox1 + fdybox) * zoomdy);
  dx = x2-x1+1;
  dy = y2-y1+1;
  if (x < x1 || x > x2 || y < y1 || y > y2) return;

#ifdef UNIX
  if (xhair_type==XHAIR_CURSOR) ;
  else if (xhair_type==XHAIR_XOR || xhair_type==XHAIR_1ST_LINE)
    {
      if (xhair_type==XHAIR_XOR) set_xor(zoomgc, 1);
      XDrawLine(mydisplay, zoomwin, zoomgc, x1, y, x2, y);
      XDrawLine(mydisplay, zoomwin, zoomgc, x, y1, x, y2);
    }

  else
    {
  /* This method involved saving H & V lines as images, then
   * using XPutImage.  Problem if resize window: didn't resize lines
   */
      if (xiH == 0)
	{
	  xih = XGetImage(mydisplay, zoomwin, x1,y, dx,1, AllPlanes,XYPixmap);
	  xiv = XGetImage(mydisplay, zoomwin, x, y1,1, dy,AllPlanes,XYPixmap);
	  XDrawLine(mydisplay, zoomwin, zoomgc, x1, y, x2, y);
	  XDrawLine(mydisplay, zoomwin, zoomgc, x, y1, x, y2);
	  xiH = XGetImage(mydisplay, zoomwin, x1,y, dx,1, AllPlanes,XYPixmap);
	  xiV = XGetImage(mydisplay, zoomwin, x, y1,1, dy,AllPlanes,XYPixmap);
	  XPutImage(mydisplay, zoomwin, zoomgc, xih, 0,0,x1,y, dx,1);
	  XPutImage(mydisplay, zoomwin, zoomgc, xiv, 0,0,x,y1, 1,dy);
	}

#ifdef TRY_XORP
  /* This image was TRY_SETP with the hope that XPutImage() responds
   zoom* to set_xor().  It didn't!
   */
      set_xor(zoomgc, 1);
      XPutImage(mydisplay, zoomwin, zoomgc, xiH, 0,0,x1,y, dx,1);
      XPutImage(mydisplay, zoomwin, zoomgc, xiV, 0,0,x,y1, 1,dy);
#endif
/*
Note about bad response time on quetzal.  Can have "state=0" here and
nothing on screen for a while.  Then quick flash on an intermediate crosshair,
then the final one steady.  Means it takes a while for the X-function
events to be implemented.
*/
      if (state==1)
	{
	  XPutImage(mydisplay, zoomwin, zoomgc, xih, 0,0,x1,y, dx,1);
	  XPutImage(mydisplay, zoomwin, zoomgc, xiv, 0,0,x,y1, 1,dy);
	  state=0;
	}
      else
	{
	  xih = XGetImage(mydisplay, zoomwin, x1,y, dx,1, AllPlanes,XYPixmap);
	  xiv = XGetImage(mydisplay, zoomwin, x, y1,1, dy,AllPlanes,XYPixmap);
	  if (xhair_type==XHAIR_LINE)
	    {
	      XDrawLine(mydisplay, zoomwin, zoomgc, x1, y, x2, y);
	      XDrawLine(mydisplay, zoomwin, zoomgc, x, y1, x, y2);
	    }
	  else
	    {
	      XPutImage(mydisplay, zoomwin, zoomgc, xiH, 0,0,x1,y, dx,1);
	      XPutImage(mydisplay, zoomwin, zoomgc, xiV, 0,0,x,y1, 1,dy);
	    }
	  state=1;
	}
    }

  set_xor(zoomgc, 0);

#else			/* UNIX is not defined */
  Drawcrosshair(x, y, x1, y1, x2, y2);
#endif			/* of UNIX defined */
}
