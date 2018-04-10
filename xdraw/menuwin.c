/******************************************************************************
**  NAME	MENUWIN.C
**  AUTHOR	Sheryl M. Glasser
**
**  DESCRIPTION: Interface for menus
**
**	void set_menubkgd(bkgd)
**	long newmenu()
**	void addtomenu(id,text,&arg,submenu)
**	void enablemenus(main_id, title, x, y, func)
**	void exec_menu()
**	void init_menupop(graphics_win,key)
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
**
**  Changes, to return to condor (taken from glmenu.c 2-11-94):
**  1.  retval removed: from struc, addtomenu
**  2.  enable_menus now takes mtitle[], x,y; NOTE - must alter enablemenus
**  3.  need label formats include %s
**  4.  get rid of popup and manywin
**  5.  zap set_menutype(), get_menutype(), append_pupcode(),
**      reset_popup(), add_retvals(), get_mname
**  6.  no popup ==> don't need menuexec() as intermediary
**  7.  .value becomes .submenu; addtomenu arg=sub, not &sub
**  8.  zap getmenustring(), give mname to eg. do_menu()
**  9.  menu_event takes short int! But x version needs long(?)
** 10.  on SGI must compile with #define SGI
** 11.  set_menubkgd args now int
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>		/* for exit() */
#include <fcntl.h>
#include <malloc.h>
#include <string.h>
#include <malloc.h>

#include "menuwin.h"
#include "device.h"

#ifndef UNIX
extern void Drawcursor(void);
#endif

#ifdef SGI
#include "event.h"
#include <gl.h>
#define USE_EVENT 1

#else
#ifdef UNIX
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#define WBORD 5
#define ONE 1

#else
#include "xlib.h"
#include "xutil.h"
#define WBORD 20
#define ONE 2
#endif

#include "glx.h"
#endif
#define mprintf(s) { printf(s);fflush(stdout); }	/* debug menu printf */
#define mprintf2(s,a,b) { printf(s,a,b);fflush(stdout); }
#define mprintf1(s,a) { printf(s,a);fflush(stdout); }
/*==============================================================================
|		DEFINITIONS AND VARIABLES
==============================================================================*/
extern int zdebug;
static int zapcurrent=0;
static int shadow=1;
static int menu_visible = 1;

static void (*menuexec) (char *) = NULL;

#define cBLACK	  0
#define cWHITE   -1
#define cGREY    0xff404040
#define iGREY    0xff808080
#define iiGREY   0xffc0c0c0
#define cCYAN    0xff404000
#define iCYAN    0xffffff00
#define iMAGENTA 0xffff00ff
#define cBLUE    0xffffc000

#define FRAC .5

static long wbkgd=0xff504040;
static long wfore=0xffffffff;		/* (unused) */

static long mbkgd=cWHITE, mfore=cBLACK;		/* bkgd, fore of normal menu */
static long hbkgd=iiGREY, hfore=cBLACK;		/* of highlit item */

static long bold=iMAGENTA, hbold=iMAGENTA;
static long mshadow=cBLACK;

#define DYH1 4
#define DYH2 0

#define ytext(i)    (yindent + i*sep + (i+1)*fontheight)
#define xhi1(x1)    (x1+xframe+ONE)
#define xhi2(x2)    (x2-xframe-ONE)
#define yhi1(y)     (y-DYH1)
#define yhi2(y)     (y+fontheight+DYH2)

#define NUL	    '\0'

struct LABEL {
  char *text;		/* text string, including '\' separating states */
  char key[2];		/* Key press character for this label */
  char type;		/* 'x', 'm' */
  int *state;		/* 0,1 etc */
  int laststate;	/* state used the last time drawn */
  int nstates;		/* number of states */
  long submenu;		/* submenu id */
  struct LABEL *next;	/* next entry in menu, NULL if done */
} ;

#define NL sizeof(struct LABEL)

struct MENU {
  long id;		/* menu identifier */
  long win;		/* window identifier */
  struct LABEL *first;	/* first menu label */
  int n;		/* label count */
  int width;		/* width of menu, from strwidth() */
  int height;		/* height of menu */
  int x0,y0;		/* upper left corner */
  int level;		/* nesting level: 1,2 etc */
  int visible;		/* 0 if not visible, 1 if visible */
  int lit;		/* which label is highlighted, -1 if none */
  int oldlit;		/* which label was lit before */
  struct MENU *from;	/* if submenu, pointer to calling menu */
} ;

long cursorwin;		/* Global! */

static int nmenu = 0;
static struct MENU menulist[20], *menu2;
static struct MENU *mptr, *mainp, *topp;
static char mname[20];

static int winx1=0,winy2=0;	/* coordinates for main menu UL */
static int xbord=0,ybord=0;	/* sep main menu from window */
static int xindent,yindent;	/* sep text from menu boundary */
static int xframe=0,yframe=0;	/* inner line of menu */
static int xshad=0,yshad=0;
static int fontheight, fontdesc, sep;
static int ndialog=0;
static char *dtext;

static int dxmax,dymax;
static int upevent_ok=0;
static long otherwin=0;		/* load via inform_menu() */
static int popkey=0;

#ifdef SGI
static Matrix idm = {
  1.0, 0.0, 0.0, 0.0,
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0,
  0.0, 0.0, 0.0, 1.0
  } ;
#endif

static struct MENU *getmenuptr(long id);
static void getlabel(char *s, int *n, int max, char type, char *label,
		     int count);
static int formatted(char *);
static int statecount(char *);
void exec_menu(void);
static void openmenu(long menu_id, char *title);
static void drawmenu(struct MENU *p, int all);
static int selected_index(void);
static void enter_menu(long win);
static void zapmenu(struct MENU *p);
static int  zapdaughters(void);

static void menu_lm(void *dummy, int state);
static void menukey(void *, short int);
static void pushmenu(void *a, int b);

static struct MENU *topwin(void);
static void go_menu(int ind);
static void getbox(int);

#define dprintf printf
#define strCat(lab,p,p2) {c=*p2; *p2=NUL; strcat(lab,p); *p2=c; p=p2+1;}

/*==============================================================================
|		FUNCTIONS FOR DEFINING MENUS
==============================================================================*/
/*-----------------------------------------------------------------------------
|	set_menu_visibility
-----------------------------------------------------------------------------*/
void set_menu_visibility(int v)
{
  menu_visible = v;
}

/*-----------------------------------------------------------------------------
|	set_menubkgd, set_menufore
|	R,G,B, values 0 to 255
|	0=total menu, 1=submenu, 2=highlight, (3=bold, 4=hbold)
-----------------------------------------------------------------------------*/
void set_menubkgd(int which, int r, int g, int b)
{
  long bkgd;
  bkgd = 0xff000000 | (r & 0xff) | (g & 0xff) * 0x100 |
          (b & 0xff) * 0x10000;
  if (which==0) wbkgd = bkgd;
  if (which==1) mbkgd = bkgd;
  if (which==2) hbkgd = bkgd;
}

void set_menufore(int which, int r, int g, int b)
{
  long fore;
  fore = 0xff000000 | (r & 0xff) | (g & 0xff) * 0x100 |
          (b & 0xff) * 0x10000;
  if (which==0) wfore = fore;
  if (which==1) mfore = fore;
  if (which==2) hfore = fore;
  if (which==3)  bold = fore;
  if (which==4) hbold = fore;
}

unsigned long get_menubkgd() { return(wbkgd); }

/*-----------------------------------------------------------------------------
|	newmenu -- initialize a new menu
-----------------------------------------------------------------------------*/
long newmenu(void)
{
  struct MENU *p;
  long id;
  static int inited=0;

  if (!inited)
    {
      fontdesc = (int)getdescender();
      fontheight = (int)getheight() - fontdesc;
      sep = fontheight/2;
      xbord = ybord = WBORD;
      xframe = yframe = sep;
      xindent = xframe+10;
      yindent = yframe+5;
      /*dprintf("Font height=%ld, descender=%ld\n",fontheight,fontdesc);*/
      inited++;
    }

  p = &menulist[nmenu];
  id = nmenu+100;
  p->id = id;
  p->win = -1;
  p->from = NULL;
  p->first = NULL;
  p->n = 0;
  p->width = 0;
  p->height = -sep;		/* ultimately, n*fontheight+(n+1)*sep */
  p->x0 = p->y0 = -1;
  p->visible = p->level = 0;
  p->lit = p->oldlit = -1;
  nmenu++;
  menu2 =p+1;
  return(p->id);
}
/*-----------------------------------------------------------------------------
|	addtomenu -- add a line to menu 'id'
-----------------------------------------------------------------------------*/
void addtomenu(long id, char *text, int *arg, long *subp)
{
  struct LABEL *q, *qn;
  struct MENU *p, *pm;
  int i;
  long w, extra;
  char label[80], type, *k;
  static char btext[] = "^";

  if ((p=getmenuptr(id)) == NULL)		/* 'id' must be defined */
    return;
  scale_font(12);
  extra = getstrwidth(btext);
  extra=0;

  if (subp==NULL) type = 'x';
  else if (*subp!=0) type = 'm';
  else return;					/* submenu not inited: zap */

  q = (struct LABEL *)malloc((size_t)NL);	/* allocate label */
  if (p->first==NULL) p->first = q;		/* menu points to it */
  else
    {
      for(qn=p->first; qn->next; qn=qn->next) ;
      qn->next = q;
    }
  p->n++;
						/* initialize label */
  q->text = text;
  q->key[0] = *text;
  q->key[1] = 0;
  k = strchr(text,'^');
  if (k)
    {
      q->key[0] = *(k+1);
      k = strchr(k+1,'^');
      if (k) q->key[1] = *(k+1);
    }
  if (q->key[0]>='a' && q->key[0]<='z') q->key[0] += ('A'-'a');
  if (q->key[1]>='a' && q->key[1]<='z') q->key[1] += ('A'-'a');

  q->type = type;
  q->state = arg;
  q->nstates = statecount(text);
  q->laststate = arg ? *arg : 0;
  q->submenu = subp ? *subp : 0;
  q->next = NULL;

  if (type=='m')			/* if calls submenu, set 'from' */
    {
      for(pm=menulist; pm<menu2 && pm->id!=*subp; pm++) ;
      if (pm<menu2) pm->from = p;
    }

  p->height += (fontheight+sep);

  for(i=0; i<q->nstates; i++)		/* get max size of label */
    {					/* try all states of label */
      getlabel(text,q->state,q->nstates,type,label, i);
      w = getstrwidth(label);
      if (k!=NULL) w -= extra;
      w += 2*extra;
      if (w > p->width) p->width=(int)w;
    }

}

/*-----------------------------------------------------------------------------
|	statecount
-----------------------------------------------------------------------------*/
static int statecount(char *s)
{
  char *p1,*p2;
  int n;
#ifdef DEAD
  printf("In statecount() for ");
  for (p1=s; *p1; p1++)
    { if (*p1>=255) printf("(hex)%x",*p1); else printf("%c",*p1); }
  printf("\n");
#endif
  if (formatted(s)) return(1);
  if ((p2=strchr(s,'|'))==NULL) return(1);
  for(n=2,p1=p2+1; p2=strchr(p1,'|'); n++,p1=p2+1) ;
  return(n);
}

/*-----------------------------------------------------------------------------
|	getmenuptr -- get pointer in menulist[] having ->id
-----------------------------------------------------------------------------*/
static struct MENU *getmenuptr(long id)
{
  struct MENU *p;
  for(p=menulist; p<menu2 && p->id != id; p++) ;
  return((p==menu2) ? NULL : p);
}

/*-----------------------------------------------------------------------------
|	win_is_menu
-----------------------------------------------------------------------------*/
int win_is_menu(long w)
{
  return(menu_visible ? w==mainp->win : 0);
}

/*=============================================================================
|		OPEN, DRAW MENUS
=============================================================================*/
/*-----------------------------------------------------------------------------
|	openmenu -- open, display, enable a menu
|	from enablemenus() and LEFTMOUSE event in 'm' menu
-----------------------------------------------------------------------------*/
static void openmenu(long menu_id, char *title)
{
  int n;
  struct MENU *p;
  int mx1,mx2,my1,my2;
  long dx,dy;

  if ((p=getmenuptr(menu_id)) == NULL) return;
  topp = p;

  if (!menu_visible || p!=mainp)
    p->win = p->id;

  else if (p==mainp)			/* Open a new window if needed */
    {
      n = p->n;
      mx1 = p->x0 + winx1;
      my2 = p->y0 + winy2;
      mx2 = mx1 + dxmax;
      my1 = my2 - dymax;

      dx = mx2-mx1;
      dy = my2-my1;
      prefposition(mx1,mx2,my1,my2);
      maxsize(dx,dy);		/* doesn't work!! */
      minsize(dx,dy);
      p->win = winopen(title);				/* open it */
      winconstraints();
      if (zdebug) printf("Menu window = %lx\n", p->win);
    }


  /*dprintf("Opening window %ld %ld\n",p->id,p->win);*/

#ifndef SGI				/* XWindows: if main window.. */
  if (p!=mainp)				/* .. wait for expose event */
#endif
  drawmenu(p, 1);					/* draw it */

#ifdef USE_EVENT
  if (p==mainp)						/* enable events */
    {
      add_event(p->win, LEFTMOUSE, DOWN, menu_lm, NULL);
      add_event(p->win, LEFTMOUSE, UP,   menu_lm, NULL);
      if (popkey)
	{
	  add_event(ANY, popkey, UP, pushmenu, NULL);
	  qdevice(popkey);
	}
      add_window_event(p->win,REDRAW, redrawmenu);

#define MISC_EVENT
#ifdef MISC_EVENT
      add_window_event(p->win,WINTHAW,redrawmenu); 	/* no events came! */
      add_window_event(p->win,WINFREEZE,redrawmenu);
      add_window_event(p->win,REDRAWICONIC,redrawmenu);
      add_window_event(p->win,PIECECHANGE,redrawmenu);
      add_window_event(p->win,DEPTHCHANGE,redrawmenu); /* yes event came */
						       /* was a REDRAW too */
#endif
    }
#endif
  p->visible = 1;
}

/*-----------------------------------------------------------------------------
|	init_menupop
-----------------------------------------------------------------------------*/
void init_menupop(long w, int k)
{
  otherwin = w;
  popkey = k;
}

/*-----------------------------------------------------------------------------
|	pushmenu
-----------------------------------------------------------------------------*/
static void pushmenu(void *a, int b)
{
  int gd,md;
  long cw;
  a; b;
  gd = windepth(otherwin);
  md = windepth(mainp->win);
  cw = winget();
  /*printf("Current window=%ld, Depth graphics=%d, menu=%d\n",cw,gd,md);*/
  if (md > gd) raise_menu();
  else lower_menu();
}

void raise_menu()
{
  winset(mainp->win);
  winpop();
}

void lower_menu()
{
  winset(mainp->win);
  winpush();
}

/*-----------------------------------------------------------------------------
|	redraw_dialog
-----------------------------------------------------------------------------*/
static void redraw_dialog(void)
{
}

/*-----------------------------------------------------------------------------
|	dialog
-----------------------------------------------------------------------------*/
void dialog(int row, char *text)
{
  row; text;
}

/*-----------------------------------------------------------------------------
|	drawmenu -- draw menu 'p'
|	'all': 1 draws entire menu, 0 draws oldlit and lit
-----------------------------------------------------------------------------*/
static void drawmenu(struct MENU *p, int all)
{
  int i,j, x,y, x1,y1,x2,y2, dflag;
  struct LABEL *q;
  char label[80];
  long oldwin, color, bcolor;
  static int needed=1;
  char *k;
  static char key[] = " ";
  float sx,sy;
  long arrow[2], arrow2[8];
  static int Arrow[] = { 0,0, 10,5, 0,10, 0,0 };
  static int dArrow[]= { -2,-2, 2,0, -2,2,  -2,-2};
  static int inited=0,solida=1,xright=0;
  static float ascale = (float).75;

  sx; sy;
  if (p==NULL) return;

  if (!menu_visible)
    {
      for(q=p->first; q; q=q->next)
	if (q->state) q->laststate = *q->state;
      return;
    }

  if (!inited)
    {
      for(j=0; j<8; j++) Arrow[j] *= (fontheight * ascale/10.);
      inited=1;
    }

  oldwin=winget();			/* save old window */
  winset(mainp->win);			/* set params for menu win */

#ifdef SGI
  if (needed)
    {
      RGBmode();
      singlebuffer();
      pushmatrix();
      loadmatrix(idm);
      sx = dxmax;
      sy = dymax;
      ortho2(-(float).5, sx, -(float).5, sy);
      cpack(wbkgd);
      gconfig();		/* WARN! window clears to black here! */
    }
#endif

  x1 = p->x0;
  y2 = p->y0 + dymax;
  x2 = x1 + p->width;
  y1 = y2 - p->height;

#ifndef UNIX
  Drawcursor();
#endif
  if (all)				/* clear window, draw 1st */
    {
      if (p==mainp)
	{
	  cpack(wbkgd);
	  clear();
	}
      cpack(mbkgd);
      rectfi(x1,y1,x2,y2);
      cpack(mfore);
      recti(x1+ONE,y1+ONE,x2-ONE,y2-ONE);
      recti(x1+xframe,y1+yframe,x2-xframe,y2-yframe);
      xright = x2 - xframe -2 - fontheight;
      if (shadow)
	{
	  cpack(mshadow);
	  rectfi(x1+yindent, y1-ONE, x2+yindent, y1-yindent);
	  rectfi(x2+ONE, y1-ONE, x2+yindent, y2-yindent);
	}
    }
  scale_font(12);

  for(i=0,q=p->first; q; i++,q=q->next)		/* draw labels */
    {
      x = x1 + xindent;
      y = y2 - ytext(i);
      getlabel(q->text, q->state, q->nstates, q->type, label, -1);
      if (q->state)
	q->laststate = *q->state;

      dflag = (all || i==p->lit || i==p->oldlit);
      if (!dflag) continue;

      color = mfore;
      bcolor = bold;
      if (i==p->lit)				/* highlight, if needed */
        {
	  cpack(hbkgd);
	  rectfi(xhi1(x1), yhi1(y), xhi2(x2), yhi2(y));
	  color = hfore;
	  bcolor = hbold;
	}
      else if (!all && i==p->oldlit)
	{
	  cpack(mbkgd);
	  rectfi(xhi1(x1), yhi1(y), xhi2(x2), yhi2(y));
	}
      cpack(color);

      cmov2i(x,y);
      k = strchr(label,'^');
      if (k==NULL)
	drawstr(label);
      else
	{
	  *k = 0;
	  if (k>label) drawstr(label);
	  *key = *(k+1);
	  cpack(bcolor);
	  drawstr(key);
	  cpack(color);
	  drawstr(k+2);
	}
      if (q->type=='m')
	{
	  cpack(color);			/* mfore or hfore (was wbkgd) */
	  if (solida) bgnpolygon();
	  else bgnline();
	  for(j=0; j<8;)
	    {
	      arrow[0] = Arrow[j] + xright;
	      arrow2[j] = arrow[0] + dArrow[j];
	      j++;
	      arrow[1] = Arrow[j] + y;
	      arrow2[j] = arrow[1] + dArrow[j];
	      j++;
	      v2i(arrow);
	    }
	  if (solida) endpolygon();
	  else endline();
#ifdef DEAD_CODE
	  cpack(bold);
	  bgnline();
	  for(j=0; j<8; j+=2) v2i(arrow2+j);
	  endline();
#endif
	}
      if (i==p->lit) cpack(mbkgd);
    }

#ifndef UNIX
  Drawcursor();
#endif
  if (needed)
    popmatrix();				/* restore old window */
  winset(oldwin);
}

/*-----------------------------------------------------------------------------
|	redrawmenu -- for REDRAW events
-----------------------------------------------------------------------------*/
void redrawmenu(void *arg1, int state)
{
  struct MENU *p,*p0;
  int i,n;
  long x1,y1,oldwin;

  state;
  if (arg1)
    {
      n = *(int *)arg1;
      if (n!=REDRAW && n!=DEPTHCHANGE) /* 528 and 543 ok */
	printf("Redraw from event %d\n",n);
    }

  oldwin = winget();
  winset(mainp->win);

  if (menu_visible) getorigin(&x1,&y1);
#ifdef DEAD_CODE
  getsize(&dx,&dy);		/* did NOT work! */
  dxmax = dx;
  dymax = dy;
#endif
  winx1 = x1;
  winy2 = y1+dymax;

  winset(oldwin);
  p0 = topwin();
  for(n=p0->level; n ;n--)
    {
      for(i=1,p=p0; i<n; i++) p=p->from;
      drawmenu(p, 1);
    }
  redraw_dialog();
}

/*-----------------------------------------------------------------------------
|	getstate -- get the label *pp..p2 corresponding to state n
|	input string pointer *pp must point past the prefix, if any
|	Return p2=NULL means last option
-----------------------------------------------------------------------------*/
static char *getstate(char **pp, int n)
{
  int i;
  char *p2, *p0, *p, *p20;
  p0 = p = *pp;
  for(i=0;; i++,p=p2+1)			/* advance to label of state n */
    {
      p2 = strchr(p,'|');		/* delimiter of i-th label p..p2 */
      if (i==0) p20=p2;
      if (p2==NULL) break;		/* n too big or last item */
      if (i==n) break;
    }
  if (i<n) { p=p0; p2=p20; }
  *pp = p;
  return(p2);
}

/*-----------------------------------------------------------------------------
|	getlabel -- load label[] with correct label for current state
-----------------------------------------------------------------------------*/
static void getlabel(char *s, int *n, int max, char type, char *label,
		     int count)
{
  char *p,*p2,*pc,*pl,*p0, c;
  int nn, arrow, options, prefix, suffix;
  static char arrowtext[] = "\256";

#define add_option() p2=getstate(&p,nn); \
if (p2==NULL) { if (suffix) p2=pc; else p2=p+strlen(p); } \
strCat(label,p,p2);

#define arrowbug(n) if (arrow && p2) printf("  %d:%s\n",n,p2); \
else if (arrow && !p2) printf("  %d gives NULL\n",n);

  if ((nn=formatted(s)) && n != NULL)
    {
      if (nn==1) sprintf(label, s, *n);
      else if (nn==2) sprintf(label, s, (char *)n);
      else if (nn==3) sprintf(label, s, *(float *)n);
      goto LDONE;
    }

  p = s;				/* p at beginning of string */
  pc = strchr(p,'}');			/* pc after prefix */
  if (!pc) pc = strchr(p,'>');
  pl = strchr(p,'|');			/* pl after 1st option */
  p2 = p+strlen(p);			/* p2: end of line */
  arrow = (pc && *pc=='>');		/* flag if label uses an arrow */
  options = (n!=NULL && pl!=NULL);
  prefix = (pc && pl && pc<pl);
  suffix = (pc && pl && pc>pl);

  if (!options)
    {
      strcpy(label,s);
      goto LDONE;
    }

  *label = '\0';
  if (prefix)				/* get prefix of options, if any */
    {
      strCat(label,p,pc);		/* write prefix p..pc */
      if (!arrow) strcat(label," ");
      else strcat(label,":");
      p = p0 = pc+1;			/* now p = 1st option */
      pc = strchr(p,'}');		/* and pc = suffix pointer */
      suffix = (pc!=NULL);
    }

  if (options)				/* There's an integer and options */
    {					/* e.g. "Spin}On|Off", &spin */
      nn = (count==-1) ? *n : count;
      if (nn<0 || nn>=max) nn=0;
      add_option();
    }

  if (arrow)				/* append arrow and 2nd option */
    {
      strcat(label,arrowtext);
      p = p0;
      nn++; if (nn>=max) nn=0;
      add_option();
    }

  if (suffix)				/* Suffix */
    {
      strcat(label," ");
      strcat(label,pc+1);
    }

LDONE:
  if (type=='m') strcat(label," ");	/* space preceds > to submenu */
}

/*-----------------------------------------------------------------------------
|	formatted
-----------------------------------------------------------------------------*/
static int formatted(char *s)
{
  char *p,*p2;
  for(p=s; *p;p=p2+1)			/* Test for formatted text, eg. %d */
    {
      p2=strchr(p,'%');
      if (p2++ == NULL) return(0);
      for(p=p2;*p=='.' || (*p>='0' && *p<='9'); p++) ;
      if (*p=='d') return(1);	/* (currently only support %d) */
      if (*p=='s') return(2);
      if (*p=='f') return(3);
    }
  return(0);
}

/*=============================================================================
|		EXECUTION, EVENT PROCESSING
=============================================================================*/
/*-----------------------------------------------------------------------------
|	enablemenus
|	all menus and labels should be defined first
|	  mmain: the id of the main menu
|	  func:  the function to be called for an executable selection
-----------------------------------------------------------------------------*/
void enablemenus(long mmain, char *title, int x, int y,
		 void (*func) (char *))
{
  long oldwin;
  static int inited=0;
  void (*tmpfn)();

  tmpfn = func;
  if (inited) return;
  inited=1;

  menuexec = func;
  mainp = getmenuptr(mmain);
  *mname=0;
  winx1 = x; winy2 = y;

  getbox(0);					/* set all ->x0..height */
#ifdef USE_EVENT
  set_inputchange(enter_menu);			/* enable enter_menu() */
  add_event(ANY, KEYBD, ANY, menukey, NULL);
  qdevice(KEYBD);
#endif

  oldwin=winget();
  openmenu(mainp->id, title);			/* open,display,enable menu */
  winset(oldwin);

}

/*-----------------------------------------------------------------------------
|	menu_event, menu_ievent
|	* menu_ievent() is for xwindows, INPUTCHANGE
-----------------------------------------------------------------------------*/
int menu_event(int dev, short int state)
{
  int in;
  if (dev==KEYBD) { menukey(NULL, state); return(1); }
  if (!menu_visible) return(0);

  in = (cursorwin==mainp->win);

  if (dev==INPUTCHANGE)
    cursorwin = state;
  else if (dev==LEFTMOUSE && in)
    menu_lm(NULL, state);
  else if (dev==REDRAW && in)
    redrawmenu(NULL,state);
  else return(0);
  return(1);
}

void menu_ievent(long newwin)
{
  if (zdebug) printf("menu_ievent() changes cursorwin to %lx\n", newwin);
  cursorwin = newwin;
}

long get_cursorwin()
{
  return(cursorwin);
}

/*-----------------------------------------------------------------------------
|	enable_menus
|	this is companion to menu_event, doesn't need event.c
|	all menus and labels should be defined first
|	* mmain: the id of the main menu
|	* x,y: upper left corner of menu, rel. to 0,0 = lower left
|	  x=-1 ==> y=0..8 for UL, UC, UR, ... LR
-----------------------------------------------------------------------------*/
void enable_menus(long mmain, char *title, int x, int y, int ndialog,
		  void (*func)(char *))
{
  long oldwin;
  static int inited=0;

  if (inited) return;
  inited=1;

  menuexec = func;
  mainp = getmenuptr(mmain);

  getbox(ndialog);				/* set all ->x0..height */
  winx1 = x; winy2 = y;

  qdevice(KEYBD);
  qdevice(INPUTCHANGE);
  qdevice(REDRAW);
  qdevice(LEFTMOUSE);
  if (popkey) qdevice(popkey);

  oldwin=winget();
  openmenu(mainp->id, title);			/* open,display,enable menu */
  winset(oldwin);

  *mname=0;
}

/*-----------------------------------------------------------------------------
|	getbox
|	For each menu, set level, width,height,x0,y0
|	* width initially strwidth()+extra
|	* height initially n*(fontheight+sep) - n
|	* x0,y0 = UL of box within win
-----------------------------------------------------------------------------*/
static void getbox(int mdialog)
{
  struct MENU *p,*pf,*p2;
  struct LABEL *q;
  int i,n, nmax,x2,y1;
  int xmin,xmax,ymin,ymax;
  for(p=menulist, nmax=0; p<menu2; p++)		/* For each menu, get.. */
    {						/* ..nesting level */
      for(p->level=1,p2=p; p2=p2->from; p->level++) ;
      if (p->level > nmax) nmax = p->level;
    }

  xmax = ymin = 0;
  for(n=1; n<=nmax; n++)			/* For each menu at.. */
    for(p=menulist; p<menu2; p++)		/* ..a given level */
      {
	if (p->level != n) continue;
	p->width  += 2*xindent;			/* get width, height */
	p->height += 2*yindent;
	if (n==1)				/* get x0,y0 */
	  {
	    p->x0 = xmin = xbord;		/* ..for main menu */
	    p->y0 = ymax =-ybord;
	  }
	else					/* ..for daughters */
	  {
	    pf = p->from;
	    p->x0 = pf->x0 + FRAC * pf->width;
	    for(i=0,q=pf->first;		/* get line in parent */
		q && q->submenu!=p->id; i++,q=q->next) ;
	    p->y0 = pf->y0 - ytext(i) + fontheight + yindent;
	  }
	x2 = p->x0 + p->width;
	y1 = p->y0 - p->height;
	if (x2 > xmax) xmax = x2;
	if (y1 < ymin) ymin = y1;
      }
  dxmax = xmax - xmin + 2*xbord;
  dymax = ymax - ymin + 2*ybord;
  if (ndialog=mdialog)
    {
      dymax += ndialog * (fontheight + sep);
      dtext = malloc(80 * ndialog);
    }
}

/*-----------------------------------------------------------------------------
|	exec_menu
|	executable selection foundfrom LEFTMOUSE UP event in a menu
|	(option can be called from simulation mode; else unneeded,
|	 should call (*menuexec) directly in menu_lm())
-----------------------------------------------------------------------------*/
void exec_menu()
{
  if (!menu_visible) printf("Executing <%s>\n",mname);
  (*menuexec)(mname);
  if (!menu_visible) menu_root();
}					/* else see menu_lm() */

/*-----------------------------------------------------------------------------
|	menu_lm -- LEFTMOUSE process
-----------------------------------------------------------------------------*/
static void menu_lm(void *dummy, int state)
{
  static int keydown;
  int ind;
  struct MENU *p;

  dummy;

  /*---------- DOWN event ------------*/
  if (state==DOWN)
    {
       keydown = upevent_ok = 1;
       return;
     }

  /*---------- UP event -------------*/
  if (!upevent_ok) return;		/* set on DOWN, clrd in enter_menu() */
  if (!keydown) return;
  keydown = 0;
  ind = selected_index();		/* sets mptr, returns offset in menu */
  if (ind >= 0)				/* Process the selection */
    go_menu(ind);
  else if (ind!=-1)			/* In a menu, but missed */
    {					/* back to main */
      for(p=menulist; p<menu2; p++)
	{
	  p->visible= (p==mainp) ? 1 : 0;
	  p->lit = p->oldlit = -1;
	}
      *mname = 0;
      redrawmenu(NULL,0);
    }
}

/*-----------------------------------------------------------------------------
|	go_menu -- bookkeeping, execution for a menu event
-----------------------------------------------------------------------------*/
static void go_menu(int ind)
{
  struct LABEL *q;
  struct MENU *p;
  char *m;
  int i, n, backup, flag;

  for(q=mptr->first, i=0; q; q=q->next)	/* count to the item to get q */
    {					/* (mptr==menu cursor is in) */
      if (i==ind) break;
      i++;
    }
  if (!q || mptr->lit==i) return;
  mptr->oldlit = mptr->lit;

  flag=0;
  if (mptr->oldlit != -1 && mptr->oldlit != i)	/* If another is lit... */
    {
      flag = zapdaughters();			/* .. zap its submenus */
      if (flag) redrawmenu(NULL,0);
    }

  mptr->oldlit = mptr->lit;
  mptr->lit = i;
  topp = mptr;
  drawmenu(mptr, 0);				/* highlight parent win */

  for(n=0,p=mptr; p->from!=NULL; p=p->from)	/* ready mname[] */
    n++;
  *(mname+n) = 0;
  m = mname+strlen(mname);			/* append to mname */
  *m = q->key[0];
  if (!*m) *m='~';
  *(m+1) = 0;

  if (q->type == 'm')				/* open submenu, enable.. */
    openmenu(q->submenu, NULL);			/* .. and reset topp */

  else						/* or execute */
    {
      exec_menu();

      backup = (mptr != mainp) && zapcurrent;
      if (backup)
	{
	  p = mptr->from;
	  zapmenu(mptr);
	  mptr = p;
	}
      mptr->oldlit = mptr->lit;
      mptr->lit = -1;
      topp = mptr;

      if (mptr!=mainp)
      for(p=mptr->from; p; p=p->from)		/* See if any parents changed */
	{
	  for(q=p->first; q; q=q->next)
	    {
	      if (q->state!=NULL && *(q->state)!=q->laststate)
		goto DRAW;
	    }
	}

      if (!backup)		/* Nothing zapped, just draw changed labels */
	drawmenu(mptr, 0);
      else			/* Redraw entire structure, bottom up */
      DRAW:
	redrawmenu(NULL,0);
      mptr->oldlit = -1;

      if (backup) mptr = NULL;			/* Hope for INPUTCHANGE? */
    }
}

/*=============================================================================
|			MISC FOR PROCESS EVENTS
=============================================================================*/
/*-----------------------------------------------------------------------------
|	selected_index
-----------------------------------------------------------------------------*/
static int selected_index()
{
  int i, y, mx,my;
  struct MENU *p;
  struct LABEL *q;
  int x1,x2,y1,y2,n;

  mx = (int)getvaluator(MOUSEX);
  my = (int)getvaluator(MOUSEY);	/* value is over entire screen */

  p = topwin();
  n = p->level;

  for(; n; n--,p=p->from)		/* test which menu, if in at all */
    {
      x1 = p->x0 + winx1;
      y2 = p->y0 + winy2;
      x2 = x1 + p->width;
      y1 = y2 - p->height;
      if (mx>=x1 && mx<=x2 &&		/* Test if in menu at all */
	  my>=y1 && my<=y2) break;
    }
  if (!n) return(-2);
  mptr=p;

  for(q=mptr->first, i=0; q; q=q->next)
    {
      y = y2 - ytext(i);
      if (my >= yhi1(y) && my <= yhi2(y))
	return(i);
      i++;
    }
  return(-1);
}

/*-----------------------------------------------------------------------------
|	topwin -- return pointer to top menu window
-----------------------------------------------------------------------------*/
static struct MENU *topwin()
{
  struct MENU *p,*p1;
  int n;
  for(p1=menulist, n=0; p1<menu2; p1++)
    {
      if (!p1->visible || p1->level<=n) continue;
      p = p1;
      n = p->level;
    }
  return(p);
}

/*-----------------------------------------------------------------------------
|	menu_backup, menu_root
-----------------------------------------------------------------------------*/
void menu_backup()
{
  menukey(NULL, '\b');
}

void menu_root()
{
  struct MENU *p;
  for(p=topwin(); p!=mainp; p=mptr)
    {
      zapmenu(p);
      mptr = topp = p->from;
      mptr->lit = -1;
    }
  redrawmenu(NULL,0);
}

/*-----------------------------------------------------------------------------
|	menukey
-----------------------------------------------------------------------------*/
static void menukey(void *arg, short int state)
{
  char c;
  struct MENU *p;
  struct LABEL *q;
  int ind;
  arg;
  
  c = (char)state;
  if (state>='a' && state<='z')			/* make uppercase */
    c += ('A'-'a');
  if (c=='!') exit(0);			/* Panic button! */

#ifdef DEAD_CODE
  printf("Menukey detects ");
  if (state>=' ' && state<=127)
    printf("<%c>\n", c);
  else
    printf("<%x>\n",state);
#endif

  p = topwin();
  if (c=='\b'/* && p!=mainp*/)			/* backspace */
    {
      zapmenu(p);
      mptr = topp = p->from;
      mptr->lit = -1;
      redrawmenu(NULL,0);
    }
  else for(; p; p=p->from)			/* find menu, index w/ char */
    {
      for(q=p->first,ind=0; q; q=q->next)
        {
	  if (q->key[0]==c || q->key[1]==c)
	    {
	      mptr = p;
	      go_menu(ind);
	      return;
	    }
	  ind++;
        }
    }
}

/*-----------------------------------------------------------------------------
|	enter_menu
|	* called when INPUTCHANGE event (see call to set_inputchange())
|	* does bookkeeping when enter or leave a menu
|	* win = current cursor window
-----------------------------------------------------------------------------*/
static void enter_menu(long win)
{
  struct MENU *p;
  long id;

  upevent_ok=0;			/* Key was down, disable next up */
				/* Find the menu assoc with win */
  for(p=menulist; p<menu2 && p->win!=win; p++) ;

  id = (p>=menu2) ? -1 : p->id;
  /*dprintf("Entering window ID=%ld, WIN=%ld\n", id, win);*/

  if (p==menu2)			/* New window is NOT a menu */
    mptr = NULL;
  else
    mptr = p;			/* set mptr as current window */

}

/*-----------------------------------------------------------------------------
|	print_menulist
-----------------------------------------------------------------------------*/
static void print_menulist(void)
{
  struct MENU *p;
  for(p=menulist; p<menu2; p++)
    {
      printf("ID=%ld, WIN=%ld, FROM=",p->id,p->win);
      if (p->from) printf("%ld",(p->from)->id);
      else printf("[null]");
      printf(", LIT=%d\n", p->lit);
    }
}

/*-----------------------------------------------------------------------------
|	zapdaughters
-----------------------------------------------------------------------------*/
static int zapdaughters()
{
  struct MENU *p, *daughter[100], *mom;
  int i,n,nl;

  /*printf("Zapping daughters of Menu %ld, window %ld\n", mptr->id, mptr->win);*/

  for(n=0, mom=mptr; n<100;)
    {
      nl = n;
      for(p=menulist; p<menu2; p++)	/* find a menu from current 'mom'*/
	{
	  if (p->from != mom || p->win==-1) continue;
	  /*dprintf("Stacking daughter menu %ld, window %ld\n", p->id, p->win);*/
	  daughter[n++] = mom = p;
	  break;
	}
      if (n==nl) break;
    }

  for(i=0; i<n; i++)
    zapmenu(daughter[i]);
  return(n);
}

/*-----------------------------------------------------------------------------
|	zapmenu -- close a menu window, disable associated events
-----------------------------------------------------------------------------*/
static int waiting=0;

static void zapmenu(struct MENU *p)
{
  if (p->win <= 0) return;
  /*printf("Zapping window ID=%ld, WIN=%ld ",p->id,p->win);*/
  if (p != mainp)
    {
      /*printf("(not main)\n");*/
      /*pause1(p->win,REDRAW);*/		/* Wait for WINCLOSE event */
#ifdef USE_EVENT
      zap_events(p->win);
#endif
      p->win = -1;
      p->visible = 0;
      p->lit = p->oldlit = -1;
    }
  else
    {
      /*printf("(yes main, so just redraw, lit--> -1!)\n");*/
      p->lit = -1;
      drawmenu(p, 0);
    }
}


#ifdef DEAD_CODE
/*=============================================================================
**		DIALOG WINDOW
**===========================================================================*/
static int dialog_on = 0, nlines=0;
static long dialog_win = 0;
static char dline[80*10];
static int diadx=0, diady=0;

/*-----------------------------------------------------------------------------
|	redraw_dialog
-----------------------------------------------------------------------------*/
static void redraw_dialog()
{
  int i;
  char *p;

  winset(dialog_win);
  cpack(mbkgd);
  rectfi(0,0,diadx,diady);
  cpack(mfore);
  for(i=0,p=dline; i<nlines; i++,p+=80)
    {
      cmov2i(xindent, diady - ytext(i));
      drawstr(p);
    }
}

/*-----------------------------------------------------------------------------
|	open_dialog
-----------------------------------------------------------------------------*/
static void open_dialog(int nrow, char *title)
{
  int x1,x2,y1,y2;

  if (dialog_on)
    {
      if (nrow) return;
      dialog_on = 0;
      return;
    }

  winset(mainp->win);
  getorigin(&x1,&y1);
  diadx = (int)dxmax;
  diady = nrow * (fontheight + sep);
  x2 = x1 + diadx;
  y2 = y1 + diady;
  prefposition(x1,x2,y1,y2);
  maxsize(diadx,diady);		/* doesn't work!! */
  minsize(diadx,diady);
  dialog_win = winopen(title);	/* open it */
  winconstraints();

#ifdef USE_EVENT
  add_window_event(dialog_win,REDRAW, redraw_dialog);
#endif
  dialog_on = 1;
}

/*-----------------------------------------------------------------------------
|	dialog
-----------------------------------------------------------------------------*/
void dialog(int row, char *text)
{
  int n;
  char *p;

  n = strlen(text);
  p = dline + 80 * row;
  if (n < 80) strcpy(p, text);
  else { strncpy(p,text, 79); *(p+79)=0; }
  if (row > nlines-1) nlines = row+1;
}
#endif
