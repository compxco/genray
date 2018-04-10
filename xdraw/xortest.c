/******************************************************************************
**  NAME      xortest.c
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      draw 2 lines - one "set", one "xor
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

extern void redraw(void);

Display *mydisplay;
int myscreen;
GC mygc;
Colormap cmap;
unsigned long myforeground, mybackground;
Window root, mywin;
XSizeHints hint;

/*-----------------------------------------------------------------------------
|    main
-----------------------------------------------------------------------------*/
void main(argc, argv)
int argc;
char *argv[];
{
  unsigned int xscreen, yscreen, bwidth, depth;
  int x0, y0;
  XEvent myevent;
  static char *title = "Test xor";

  mydisplay = XOpenDisplay("");
  if (mydisplay == NULL) { printf("No xterm\n"); exit(0); }
  myscreen = DefaultScreen(mydisplay);
  cmap = DefaultColormap(mydisplay, myscreen);
  mybackground = BlackPixel(mydisplay, myscreen);
  myforeground = WhitePixel(mydisplay, myscreen);
  root = DefaultRootWindow(mydisplay);
  XGetGeometry(mydisplay, root, &root, &x0, &y0,
	       &xscreen, &yscreen, &bwidth, &depth);

  hint.x = hint.y = 0;
  hint.width = hint.height = yscreen / 3;
  hint.flags = PPosition | PSize;
  mywin = XCreateSimpleWindow(mydisplay, root, hint.x, hint.y,
			      hint.width, hint.height, 5,
			      myforeground, mybackground);
  XSetStandardProperties(mydisplay, mywin, title, title, None,
			 argv, argc, &hint);

  mygc = XCreateGC(mydisplay, mywin, (long)0, 0);
  XSetBackground(mydisplay, mygc, mybackground);
  XSetForeground(mydisplay, mygc, myforeground);

  XSelectInput(mydisplay, mywin, KeyPressMask | ExposureMask);
  XMapRaised(mydisplay, mywin);
  printf(
"Description: Window contains two diagonal white lines.\n\
Upper-left-to-lower-right is drawn using GXcopy.\n\
Upper-right-to-lower-left is drawn using GXxor.\n\n\
Press any key with graphics window as active window to end.\n");

  for(;;)
    {
      XNextEvent(mydisplay, &myevent);
      if (myevent.type == Expose) redraw();
      else if (myevent.type == KeyPress) break;
    }

  XDestroyWindow(mydisplay, mywin);
  XFreeGC(mydisplay, mygc);
  XCloseDisplay(mydisplay);
}

/*-----------------------------------------------------------------------------
|	redraw
-----------------------------------------------------------------------------*/
void redraw()
{
  int x1, y1;
  XGCValues values;

  x1 = hint.width;
  y1 = hint.height;
  XDrawLine(mydisplay, mywin, mygc, 0, 0, x1, y1);

  values.function = GXxor;
  XChangeGC(mydisplay, mygc, GCFunction, &values);
  XDrawLine(mydisplay, mywin, mygc, 0, y1, x1, 0);
  values.function = GXcopy;
  XChangeGC(mydisplay, mygc, GCFunction, &values);
}
