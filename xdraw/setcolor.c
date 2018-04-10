/******************************************************************************
**  NAME        SETCOLOR.C
**  AUTHOR        Sheryl M. Glasser
**
**  DESCRIPTION
**     Supporting functions for xwindows
**
**  Copyright (c) CounterPoint Graphics 1993.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef UNIX
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

#else
#include <dos.h>
#include <graph.h>
#include "xlib.h"
#include "xutil.h"
#endif

#include "setcolor.h"
#include "ps.h"

#define XPRINTF
#include "gendefs.h"

extern Display *mydisplay;
extern GC redraw_gc;
extern int ps_modeon;
extern Font font;
extern XFontStruct *font_struct;
extern int font_height;
extern Colormap cmap;
extern int dialogwindow;
extern Window dialog_win;
extern GC dialog_gc;

static XColor xcolor;

/*-----------------------------------------------------------------------------
|   set_xor
-----------------------------------------------------------------------------*/
void set_xor(GC gc, int enable)
{
  XGCValues values;
  unsigned long mask;
  values.function = enable ? GXxor : GXcopy;
  mask = GCFunction;
  XChangeGC(mydisplay, gc, mask, &values);
}

#define C1 (float)60
#define C2 (float)300
#define DC (float)(C2-C1)
#define HDC (float)(DC/2)

/*-----------------------------------------------------------------------------
|   setcolor -- for k=0 to ncurve-1
|               k: 0 to ncurve-1
|               how: -1=linear, 0..n = isep ==> asinh or exp
-----------------------------------------------------------------------------*/
void setcolor(int k, int ncurve, int isep)
{
  float hue, dhue;

  hue = C1;
  if (ncurve == 1)
    dhue = (float) 0;
  else if (isep == -1)
    dhue = DC * (float) k / (float) (ncurve - 1);
  else
    dhue = fancycolor(k, ncurve, isep);

  hue += dhue;
  if (!ps_modeon) hlsrgb(hue, (float) .5, (float) 1.);
  else ps_color(hue, (float).5, (float)1.);
}

/*-----------------------------------------------------------------------------
|   fancycolor
|       dhue =  A(1-exp(-B(i-isep)) + C         i>=isep
|            = -A(1-exp( B(i-isep)) + C         i< isep
|       i    =  0, n-1, infinity, isep+di
|       dhue =  0,  DC,    A,       fA
-----------------------------------------------------------------------------*/
float fancycolor(int i, int n, int isep)
{
  static int inited = 0;
  static float a, c;
  static double b;
  float a1, a2;
  double f;
  int less, di;

  if (n < 4)
    {
      a2 = (n <= 1) ? (float) 0 : (float) i / (float) (n - 1);
      if (i == isep)
	a2 = HDC;
      return (a2);
    }
  if (inited != n)
    {
      di = n / 4;
      f = .5;
      f = 1. / (1. - f);
      b = log(f) / (double) di;
      a1 = (float) exp(-(double) (b * (n - 1 - isep)));
      a2 = (float) exp(-(double) (b * isep));
      a = DC / ((float) 2 - a1 - a2);
      c = a * ((float) 1 - a2);
      inited = n;
    }
  less = (i <= isep);
  f = b * (double) (less ? i - isep : isep - i);
  a1 = a * (float) (1. - exp(f));
  a2 = less ? c - a1 : c + a1;
  return (a2);

#ifdef DEAD_CODE
  if (isep < 0)
    dhue = DC * (float) k / (float) (ncurve - 1);
  else if (k <= isep)
    dhue = HDC * (float) k / (float) isep;
  else
    dhue = HDC + (DC - HDC) * (float) (k - isep) / (float) (ncurve - 1 - isep);
#endif
}
/*-----------------------------------------------------------------------------
|    getcolor.  return a color index corresponding to r,g,b
-----------------------------------------------------------------------------*/
unsigned long getcolor(int r, int g, int b)
{
  xcolor.red = r;
  xcolor.blue = b;
  xcolor.green = g;
  if (!NOX)
    XAllocColor(mydisplay, cmap, &xcolor);
  return (xcolor.pixel);
}

/*--------------------------------------------------------------------
|	get_xcolor - need to limit to 256 colors (e.g.)
|	* some XAllocColor returned status 0
--------------------------------------------------------------------*/
static int get_xcolor(float f)
{
  double r;
  r = floor(255. * (double)f + .5);
  r = 65535. * r / 255.;
  return (int)r;
}

/*-----------------------------------------------------------------------------
|    hlsrgb. Converts hls color coordinates to rgb.
-----------------------------------------------------------------------------*/
void hlsrgb(float hue, float lightness, float saturation)
{
  float m1, m2;
  Status status;

  /* specify color */
  if (lightness <= (float).5)
    m2 = lightness * (saturation + 1.);
  else
    m2 = lightness + saturation - lightness * saturation;

  m1 = 2 * lightness - m2;

  xcolor.red   = get_xcolor(value(m1, m2, hue + 120.));
  xcolor.green = get_xcolor(value(m1, m2, hue));
  xcolor.blue  = get_xcolor(value(m1, m2, hue - 120.));
  xcolor.flags = DoRed | DoGreen | DoBlue;
  if (xcolor.red==0 && xcolor.green==0 && xcolor.blue==0)
    xcolor.red = xcolor.green = xcolor.blue = 65535;

  /* set color */
  status = XAllocColor(mydisplay, cmap, &xcolor);
  if (status)
    XSetForeground(mydisplay, redraw_gc, xcolor.pixel);
}

/*-----------------------------------------------------------------------------
|   value. auxiliary function used by hlsrgb
-----------------------------------------------------------------------------*/
static float value(float n1, float n2, float hue)
{
  while (hue >= (float)360.)
    hue -= (float)360.;
  while (hue < (float)0.)
    hue += (float)360.;
  if (hue < (float)60.)
    return (n1 + (n2 - n1) * hue / (float)60.);
  else if (hue < (float)180.)
    return (n2);
  else if (hue < (float)240.)
    return (n1 + (n2 - n1) * ((float)240. - hue) / (float)60.);
  else
    return (n1);
}

/*-----------------------------------------------------------------------------
|   testpalette
-----------------------------------------------------------------------------*/
void testpalette(int x1, int y1, int dx, int dy)
{
#ifdef DOS
  int i, j, k, x, y, new[4];
  char text[10], *p;
  Colormap cmap;
  XColor c;

  dy /= 8;
  dx /= 8;
  x1 /= 2;
  y1 /= 2;
  for (j = k = 0; j < 4; j++)
    {
      y = j * dy + y1;
      for (i = 0; i < 4; i++, k++)
	{
	  x = i * dx + x1;
	  _setcolor(k);
	  _rectangle(_GFILLINTERIOR, x, y, x + dx, y + dy);
	}
    }

  _settextwindow(16, 0, 21, 60);
  _settextposition(0, 0);
  for (;;)
    {				/* tweak the color table */
      printf("Enter iRGB: ");
      scanf("%s", text);
      for (i = 0, p = text; i < 4; i++, p++)
	{
	  for (; *p == ' '; p++);	/* skip spaces, if any */
	  j = *p;
	  if (j == 0)
	    return;		/* abort if not all there */
	  if (i == 0)
	    {
	      if (j >= 'a')		j += (10 - 'a');
	      else if (j >= 'A')	j += (10 - 'A');
	      else			j -= '0';
	      if (j < 0 || j > 15)	j = 1;
	      k = new[i] = j;
	      c.pixel = j;
	    }
	  else
	    {
	      j -= '0';
	      if (j < 0 || j > 3)
		j = 1;
	      new[i] = j;
	    }
	}
      c.red = new[1];
      c.green = new[2];
      c.blue = new[3];
      XStoreColor(mydisplay, cmap, &c);
    }
#endif
}

/*-----------------------------------------------------------------------------
|   getgoodfont ().  Finds an appropiate font name for use in opendisplay.
|   fonttitle: Name of font returned.
|   lowpoint..highpoint: font height
|   Returns 0 if not found, 1=found
-----------------------------------------------------------------------------*/
getgoodfont(char *fonttitle, int lowpoint, int highpoint,
	    char *searchfontname)
{
#ifdef MSWINDOWS
  return 1;
#else
  /* declarations */
  int numfonts,			/*1- there is a font 0- there is no font */
    pointctr,			/*point size being searched for */
    ERRORlevel;			/* 1 GREAT 2 WRONG NAME 3 WRONG SIZE & NAME
                               0 BUST: completely failed (as if not called)*/

  char fsearchstr[80],		/*what are we looking for?*/
  **myfontlist;			/*what we found*/

  /* find what fonts are availible, that are between lowpoint
     and highpoint points */

#define maxfonts  1		/* # of fonts to put on list (don't change)*/

  numfonts = 0;
  pointctr = highpoint;

  /*start search with fsearchstr */
  ERRORlevel = 1;
  myfontlist = XListFonts(mydisplay, "*", maxfonts, &numfonts);	/* get ANY font */
  numfonts = 0;

/*
Hunt until you findsomething, or until you have exhausted all points in range.
You are starting with any font in myfontlist.
Create pattern names like "*times-medium-r-normal*--12*"
(The (*) are probably wild cards.)  Test if it's there.
*/
  while (!numfonts && pointctr >= lowpoint)
    {
      sprintf(fsearchstr, "*%s*--%d*", searchfontname, pointctr);	/* create desired name */
      XFreeFontNames(myfontlist);
      myfontlist = XListFonts(mydisplay, fsearchstr, maxfonts, &numfonts);
      pointctr--;
    }

/*
Try again for the right size, if the right name didn't exist
Creating pattern name like "*--12*"
*/

  if (!numfonts)		/* didn't find desired name? */
    {
      ERRORlevel = 2;
      pointctr = highpoint;
    }

  while (!numfonts && pointctr >= lowpoint)	/* (BUG! SG) try ANY pattern */
    {
      sprintf(fsearchstr, "*--%d*", pointctr);
      XFreeFontNames(myfontlist);
      myfontlist = XListFonts(mydisplay, fsearchstr, maxfonts, &numfonts);
      pointctr--;
    }

  /*correct for final pointctr-- (BUG! SG) */
  pointctr++;

  /*can't find something?  No problem.
    Select any font that's available.*/

  if (!numfonts)
    {
      ERRORlevel = 3;
      XFreeFontNames(myfontlist);
      myfontlist = XListFonts(mydisplay, fsearchstr, maxfonts, &numfonts);
      pointctr = 12;
    }

  /*what is the string name of the font?*/
  strcpy(fonttitle, "ERROR- no fonts exist");
  if (numfonts > 0)
    strcpy(fonttitle, *myfontlist);
  else
    ERRORlevel = 0;

  XFreeFontNames(myfontlist);

  font_height = pointctr;
  return ERRORlevel;
#endif
}

/*=============================================================================
|			Dialogue Area
=============================================================================*/
#define xprintf1(s) xprintf(s,0,0,0,0,0,0,0,0)
/*-----------------------------------------------------------------------------
|	xinput
|	* s = prompt or next char
|	v* return 0=processed char, 1=<cr>, -1=<esc>
-----------------------------------------------------------------------------*/
int xinput(char *s, int *value)
{
#define NT 10
   static char text[NT+1];
   char c;
   int n;
   extern int input_enabled;
   
   if (*(s+1))
   {
      xprintf1(s);
      input_enabled = 1;
      return 0;
   }
   
   c = *s;
   if (c == 27 || c == '\n' || c == '\r')
   {
      input_enabled = 0;
      xprintf1("\n");
      if (c == 27 || !*s) return -1;
      sscanf(text, "%d", value);
      *text = 0;
      return 1;
   }
   n = strlen(text);
   if (c == '\b')
   {
      if (n == 0) return 0;
      *(text + n - 1) = 0;
      xprintf1(s);
      return 0;
   }
   if (c < '0' || c > '9' || n == NT) return 0;
   strcat(text, s);
   xprintf1(s);
   return 0;
}

/*-----------------------------------------------------------------------------
|	xprintf etc.
-----------------------------------------------------------------------------*/
#ifndef DOS
#define MAXROWS 25
#define NPERROW 80
static char buf[MAXROWS*NPERROW];
static int head=0, tail=-1;
static int xline=0, sep;
static int nrows_win, nrows_buf=0, nrows_vis=0;
static unsigned int width, height;

#define get_dlgy(row) (font_height * (row+1) + sep * row)
/*-----------------------------------------------------------------------------
|	get_dialogrows
-----------------------------------------------------------------------------*/
int get_dialogrows()
{
  int row;
  unsigned int depth;
  get_properties(dialog_win, &width, &height, &depth);
  sep = font_height / 4;
  nrows_win = (height - sep) / (font_height + sep);
  if (nrows_win > MAXROWS) nrows_win = MAXROWS;
  else if (nrows_win < 1) nrows_win = 1;

  if (tail==-1) return 0;
  if (head==0) nrows_buf = tail - head + 1;
  else nrows_buf = MAXROWS;
  nrows_vis = (nrows_win < nrows_buf) ? nrows_win : nrows_buf;
  row = nrows_vis - 1;
  if (xline==0) row++;
  return row;
}

/*-----------------------------------------------------------------------------
|	redraw_dialog
-----------------------------------------------------------------------------*/
void redraw_dialog()
{
  int i, row, nrows, y;
  char *p;
  if (tail==-1) return;

  get_dialogrows();
  i = tail - nrows_vis + 1;
  if (i < 0) i += MAXROWS;
  for(p=buf+i*NPERROW, row=0;; p+=NPERROW, row++)
    {
      y = get_dlgy(row);
      XDrawString(mydisplay, dialog_win, dialog_gc, 0, y, p, strlen(p));
      if (i==tail) break;
      if (++i==MAXROWS) i=0;
    }
}

/*-----------------------------------------------------------------------------
|	add_dialog
|	* input from xPrintf has no '\n', or one at end
|	* updates head, tail, newline, xline
-----------------------------------------------------------------------------*/
static void add_dialog(char *text);
static void add_dialog(char *text)
{
  char *p;
  static int newline=1;
  if (newline)
    {
      if (tail==-1) tail=0;
      else if (head==0 && tail < MAXROWS-1) tail++;
      else
	{
	  if (++head==MAXROWS) head=0;
	  if (++tail==MAXROWS) tail=0;
	}
      p = buf + tail * NPERROW;
      *p = 0; xline = 0;
    }
  else
    {
      if (*text == '\b') *(text + strlen(text) - 1) = 0;
      p = buf + tail * NPERROW;
      xline = XTextWidth(font_struct, p, strlen(p));
    }
  if (*text != '\b') strcat(p, text);
  newline = (strchr(text, '\n') != 0);
  if (newline) *(strchr(p, '\n')) = 0;
}

/*-----------------------------------------------------------------------------
|	xprintf  -- (note x,y = lower left corner of text)
|	* Bugs: 1. Ghosts of crosshair!
-----------------------------------------------------------------------------*/
void xprintf(const char *format, int a, int b, int c, int d,
	     int e, int f, int g, int h)
{
  int y, y1, y2, n, row;
  char text[100], *p, *q;
  extern unsigned long myforeground,mybackground;

  sprintf(text, format, a,b,c,d,e,f,g,h);
  if (dialogwindow==0)
    {
      printf(text); fflush(stdout);
    }
  else
    {
      row = get_dialogrows();
      for(p=text;;)
	{
	  q = strchr(p, '\n');
	  if (q) { c = *(q+1); *(q+1) = 0; }
	  add_dialog(p);
	  if (row >= nrows_win)		/* Scroll */
	    {
	      y1 = get_dlgy(1) - font_height;
	      y2 = get_dlgy(row - 1);
	      XCopyArea(mydisplay, dialog_win, dialog_win, dialog_gc,
		0,y1, width, y2-y1, 0, get_dlgy(0)-font_height);
	      XClearArea(mydisplay, dialog_win, 0, y2-font_height,
		width, font_height, False);
	      row--;
	    }
	  y  = get_dlgy(row);
	  n = strlen(p);
	  if (*p == '\b') *p = ' ';
	  XDrawString(mydisplay, dialog_win, dialog_gc, xline, y, p, n);
	  if (q==NULL || c==0) break;
	  p = q + 1; *p = c;
	  xline = 0; row++;
	}
    }
}
#endif		/* of not DOS */
