/******************************************************************************
**  NAME	MENU.C
**  AUTHOR	Sheryl M. Glasser
**
**  DESCRIPTION
**		
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
**  Changes, to return to goshawk:
**  1.  retval removed
**  2.  enablemenus now takes mtitle[], x1, y1; no more popup, set_menutype
**  3.  need label formats include %s
**  4.  remove convert to uppercase: done in menuwin
******************************************************************************/
#include <stdio.h>

#ifdef UNIX
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#else
#include "xlib.h"
#include "xutil.h"
#endif

#include "menu.h"
#include "menuwin.h"
#include "gendefs.h"
#include "curves.h"
#include "xdraw.h"
#include "xinit.h"		/* for parse_title */
#include "xtools.h"		/* for ZENABLE etc., zoom() */
#include "xcontour.h"
#include "xedit.h"
#include "ps.h"

extern CURVE_SET curveset[];
extern int ncset;
extern int exitflag, redrawflag;

static int selected_iwin = 0;
static CURVE_SET *selected_cp;
static Window selected_win;

/*----- menus --------------------*/
char mtitle[] = "xdraw menu";

struct MENU {
  char *text; int *arg; long *subp;
  } ;

#define ci(x) (int *)&x[0]
#define fi(fx) (int *)&fx

static long file, select, display, zoomm, values, edit, mymenu, menuopt;
static struct MENU mmain[] = {				/* mymenu */
  "^File",		NULL,		&file,
  "^Select Window (%d)",&selected_iwin,	&select,
  "^Display",		NULL,		&display,
  "^Zoom",		NULL,		&zoomm,
  "^Values",		NULL,		&values,
  "^Edit",		NULL,		&edit,
  "^Menu",		NULL,		&menuopt,
  "^Quit",		NULL,		NULL
  } ;

static struct MENU mmenu[] = {
  "^Raise", NULL, NULL,
  "^Lower", NULL, NULL
  } ;

static long state;
extern char prefix[];
static struct MENU mfile[] = {				/* file */
  "^Read %s?.in",		ci(prefix),  &state,
  "^Write %s?.in",		ci(prefix),  &state,
  "Write ^Postscript File",	NULL,    NULL
  } ;

static int selsize = 9;
static struct MENU msel[] = {
  "^0", NULL, NULL, "^1", NULL, NULL, "^2", NULL, NULL,
  "^3", NULL, NULL, "^4", NULL, NULL, "^5", NULL, NULL,
  "^6", NULL, NULL, "^7", NULL, NULL, "^8", NULL, NULL
  } ;

static long contour;
static int marker=0, aspect=0;
static struct MENU mdispl[] = {				/* display */
  "^Markers}On|Off",			&marker, NULL,
  "Preserve ^Aspect Ratio|^Auto Scale", &aspect, NULL,
  "^Contour Settings",			NULL,	 &contour,
  } ;

extern int zoom_on;
extern float zoomfac;
static struct MENU mzoom[] = {
  "^Zoom", &zoom_on, NULL,
  "^Original Size",   NULL, NULL,
  "^Bigger",  NULL, NULL,
  "^Smaller", NULL, NULL,
  "^Increment (%4.2f)", fi(zoomfac), NULL
  } ;
static struct MENU mvalu[] = {
  "^Coordinate", NULL, NULL,
  "^Slope", NULL, NULL
  };

static int ncurve;
static struct MENU mcont[] = {				/* contour */
  "^Halve Lines (%d)",	&ncurve, NULL,
  "^Double Lines (%d)",	&ncurve, NULL
  } ;

static struct MENU mstate[] = {				/* state */
  "^0. %s0.in", ci(prefix),NULL,
  "^1. %s1.in", ci(prefix),NULL,
  "^2. %s2.in", ci(prefix),NULL,
  "^3. %s3.in", ci(prefix),NULL,
  "^x. %s.in",  ci(prefix),NULL
  } ;

static char flabel[30]="Ftemp",xlabel[30]="Xtemp",ylabel[30]="Ytemp";
static long next;
static struct MENU medit[] = {
  "^X (%s)",	  ci(xlabel), &next,
  "^Y (%s)",	  ci(ylabel), &next,
  "^Family (%s)", ci(flabel), &next,
  "^Other",	  NULL,	NULL,
  "^Go",	  NULL, NULL,
  "^Abort",	  NULL, NULL
  } ;

static char nextvar[30]="Next", prevvar[30]="Prev", oldvar[30]="Old";
static struct MENU mnext[] = {
  "^Next (%s)",     ci(nextvar), NULL,
  "^Previous (%s)", ci(prevvar), NULL,
  "^Original (%s)", ci(oldvar), NULL
  } ;

static long create_submenu(struct MENU mu[], int n);
#define mcount(mmain) (sizeof(mmain) / sizeof(struct MENU))

void do_menus(char *);
/*-----------------------------------------------------------------------------
|	init_menus
-----------------------------------------------------------------------------*/
void init_menus()
{
  struct MENU *p;
  int i,n, x1, y1;

  nextwindow(&x1, &y1);
#ifndef UNIX
  set_menubkgd(0, 0x80,0,0xff);		/* window bkgd: purple */
  set_menubkgd(2, 0x80,0,0xff);		/* highlight bkgd: purple */
  set_menufore(2, 0xff,0xff,0xff);	/* highlight fore: white */
  set_menufore(4, 0xff,0x80,0xff);
#endif

  init_menuXlibrary();

  state    = create_submenu(mstate,mcount(mstate));	/* from file read */
  contour  = create_submenu(mcont, mcount(mcont));	/* from display */
  zoomm    = create_submenu(mzoom, mcount(mzoom));	/* ditto */
  values   = create_submenu(mvalu, mcount(mvalu));
  next     = create_submenu(mnext, mcount(mnext));

  display  = create_submenu(mdispl,mcount(mdispl));	/* from mymenu */
  file     = create_submenu(mfile, mcount(mfile));
  select   = create_submenu(msel, selsize);
  menuopt  = create_submenu(mmenu, mcount(mmenu));
  edit     = create_submenu(medit, mcount(medit));

  mymenu = newmenu();
  n = mcount(mmain);
  for(i=0,p=mmain; i<n; i++,p++)
    addtomenu(mymenu, p->text, p->arg, p->subp);

  enable_menus(mymenu, mtitle, x1, y1, 0, do_menus);
  selected_iwin = 0;
  selected_cp = curveset;
  selected_win = getwindow(selected_iwin);
}

/*-----------------------------------------------------------------------------
|	create_submenu
-----------------------------------------------------------------------------*/
static long create_submenu(struct MENU mu[], int m)
{
  long submenu;
  int i;
  struct MENU *p;

  submenu = newmenu();
  for(i=0,p=mu; i<m; i++,p++)
    addtomenu(submenu, p->text, p->arg, p->subp);
  return(submenu);
}

/*-----------------------------------------------------------------------------
|	menu_ncol
-----------------------------------------------------------------------------*/
void menu_setwindows(int n)
{
  selsize = n;
}
/*-----------------------------------------------------------------------------
|	toggle
-----------------------------------------------------------------------------*/
static void toggle(int *var, int n)
{
  (*var)++;
  if (*var >= n) *var=0;

}

/*-----------------------------------------------------------------------------
|	set_selected_iwin
-----------------------------------------------------------------------------*/
void set_selected_iwin(int iwin)
{
  int i,n;
  CURVE_SET *cp;

  for(i=n=0,cp=curveset; i<ncset; i++,cp++)
    {
      if (cp->which != '1') continue;
      if (n==iwin) break;
      n++;
    }
  selected_iwin = n;
  selected_cp = cp;
  selected_win = getwindow(selected_iwin);
  marker = (cp->flags & MARKERS) ? 1 : 0;
  aspect = (cp->flags & ASPECT)  ? 1 : 0;
  menu_root();
}

/*-----------------------------------------------------------------------------
|	do_menus
-----------------------------------------------------------------------------*/
void do_menus(char *mname)
{
  int found;
  char *p,*p2,*p3,*p4;
  Window op_window;
  char text[200];

  if (!*mname) return;
  /*printf("Menu keyword %s\n",mname);*/
  p = mname;
  p2 = p+1;
  p3 = p2+1;
  p4 = p3+1;

  found = 1;
  redrawflag = 0;
  op_window = get_op_window(selected_win);	/* selected_win or event_win */

  if	  (*p == 'Q') exitflag=1;
  else if (*p == 'S') set_selected_iwin((int)(*p2 - '0'));

  else if (*p == 'F')
    {
      if      (*p2=='P')
	{
	  parse_title(selected_cp, selected_cp->title, text);
	  postscript(text);
	}
      else if (*p2=='W') ;
    }

  else if (*p == 'E') editcurve(selected_iwin, selected_cp, p2);

  else if (*p == 'D')
    {
      if      (*p2 == 'M') toggle_markers(selected_cp, &marker);
      else if (*p2 == 'A') toggle_aspect(selected_cp, &aspect);
    }

  else if (*p == 'Z')
    {
      if      (*p2=='Z')   zoom(enable_window(op_window), ZENABLE);
      else if (*p2=='O')   zoom(op_window, ZORIGINAL);
      else if (*p2=='B' || *p2=='S') addzoom(op_window, *p2=='B');
    }

  else if (*p == 'C')
    {
      if (*p2=='D' || *p2=='H') new_ncurve(selected_cp, *p2);
    }

  else if (*p == 'M')
    {
      if      (*p2=='R') raise_menu();
      else if (*p2=='L') lower_menu();
    }

  else found=0;

  if (redrawflag) clear_redraw(selected_iwin);
}
