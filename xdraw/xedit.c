/******************************************************************************
**  NAME      XEDIT.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      (Separated from xinit.c because ms compiler heap space error)
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#define O_BINARY 0

#else
#include <sys\types.h>
#include <sys\stat.h>
#ifdef DOS
#include <graph.h>
#endif
#include <io.h>
#include <conio.h>
#include <dos.h>
#include "xlib.h"
#include "xutil.h"
#endif

#include "gendefs.h"
#include "curves.h"
#include "xinit.h"
#include "xtools.h"
static void ptitle(int, CURVE_SET *);

extern int param[];
extern LOOP loop[];		/* Loop structure of the data */
extern int nloop;
extern VARIABLE_ID varid[];	/* Variables used for plotting */
extern int nvar;

static char newtitle[200];
static int from_menus=1;
static int editing=0;

extern CURVE_SET curveset[];	/* Describes the plots */
extern CURVE_SET *cp, *cp2;
extern int ncset;

#define is_white(c) (c==' ' || c=='\t')
#define is_eol(c)   (c=='\0' || c=='\n')
#define is_crlf(c)  (c=='\n' || c=='\r' || c==27)
#define rdlen(s) strlen(fgets(s,LNAME,frd))

/*=============================================================================
**			Edit
**===========================================================================*/
#ifndef DOS
#define outtext(s) printf(s); fflush(stdout);
/*void settextcolor(int i) { i; }*/	/* defined in xinit.c*/

#else
#define outtext _outtext
#define settextcolor _settextcolor
#endif

#define crlf() sprintf(text,"%c\n",c); outtext(text)
#define TEXT0 2
#define TEXT1 8
#define TEXT2 14

/*-----------------------------------------------------------------------------
|	get_descriptive_title
-----------------------------------------------------------------------------*/
char *get_descriptive_title(int iwin, CURVE_SET *p)
{
  int flag;
  flag = from_menus;
  from_menus = 1;
  ptitle(iwin, p);
  from_menus = flag;
  return newtitle;
}

/*-----------------------------------------------------------------------------
|	ptitle -- print descriptive title of curve
-----------------------------------------------------------------------------*/
static void ptitle(int iwin, CURVE_SET *p)
{
  char text[80];
  int j, jx,jy,jf;
  int indexf, undef;
  CVAR *xp, *yp;

  if (!from_menus)				/* print window number */
    {
      settextcolor(TEXT2);
      sprintf(text,"Window %d: ",iwin);
      outtext(text);
    }

  printvar(0, NULL, newtitle);		/* Initialize newtitle[] */
  yp = &p->iy;
  xp = &p->ix;
  indexf = p->lfaml | LOOPBIT;
					/* Add text to newtitle[] */
  jy = printvar(yp->index, "%s vs. ", newtitle);
  jx = printvar(xp->index, "%s, ", newtitle);
  jf = printvar(indexf, "all %s", newtitle);

  for(j=undef=0; j<nloop; j++)
    {
      if (p->param[j]<0) continue;
      strcpy(text, ", ");		/* Create format for 'other' */
      strcat(text, "%s = ");
      printvar(j | LOOPBIT, text, newtitle);
      if (p->param[j] < 0)
	{ printvar(-1, "UNDEFINED", text); undef++; }
      else sprintf(text,"%d",p->param[j]);
      strcat(newtitle,text);
    }
  if (p->flags & SINGLE)
    {
      sprintf(text,", curve index=%d\n",p->ncurve);
      strcat(newtitle,text);
    }

  if (!from_menus)
    {
      outtext(newtitle);
      outtext("\n");
      settextcolor(TEXT0);
    }
}

/*-----------------------------------------------------------------------------
|	print_title
-----------------------------------------------------------------------------*/
void print_title(int iwin, CURVE_SET *p)
{
  ptitle(iwin, p);
  printf(newtitle);
  printf("\n");
}

/*-----------------------------------------------------------------------------
|	pmenu -- print main menu (all=1), or prompt and return a char
-----------------------------------------------------------------------------*/
char pmenu(int all)
{
  VARIABLE_ID *p;
  int i,ic;
  char text[30], c;
  static int showvar=0;

  if (all)
    {
      settextcolor(TEXT1);
      if (showvar)
        {
          outtext("Variables: ");
          for(i=0,p=varid; i<nvar; i++,p++)
            {
	      outtext(p->name);
	      if (i<nvar-1) {outtext(", ");}
	      else {outtext("\n");}
	      }
	  }
      outtext("*********************************************************************\n");
      outtext("| Menu: Y=change Y, X=change X, F=change Family, O=set Other variable\n");
      outtext("|       1=one family member, 2=all family members\n");
      outtext("|       M=list this Menu, G=Go (use new settings), Q=Quit (with old)\n");
      outtext("*********************************************************************\n");
      settextcolor(TEXT0);
      return(0);
    }

  outtext("  Y,X,F,O, 1,2, M,G,Q ==> ");
  for(;;)
    {
      ic = readkey();
      c = *(char *)&ic;
      c = toupper(c);
      if (c=='X' || c=='Y' || c=='F' || c=='O') break;
      if (c=='1' || c=='2') break;
      if (c=='M' || c=='G' || c=='Q') break;
      if (is_crlf(c)) { c=' ';  break; }
    }
  if (c=='F') strcpy(text,"Family\n");
  else if (c=='O') strcpy(text,"Other\n");
  else sprintf(text,"%c\n", c);
  outtext(text);
  return(c);
}

/*-----------------------------------------------------------------------------
|	newcode
|	  'which': bit 0==>variables, bit 1==>loop variables
|	  3=x or y, 2=f or o
-----------------------------------------------------------------------------*/
int newcode(int which)
{
  int i,j, n,m, isloop;
  char text[200], s[20], c;
  int allowed[20];
  VARIABLE_ID *q;

AGAIN:
  strcpy(text,"  Enter ");
  for(i=n=0,q=varid; i<nvar; i++,q++)
    {
      isloop = (q->index & LOOPBIT);
      if (!isloop && !(which & 1)) continue;
      if ( isloop && !(which & 2)) continue;
      if (n) strcat(text,", ");
      sprintf(s,"%d=%s", i, q->name);		/* eg. "6=psifac" */
      allowed[n++] = i;
      if (strlen(text) + strlen(s) > 60)
	{
	  strcat(text,"\n"); outtext(text);
	  strcpy(text,"  ");
	}
      strcat(text, s);
    }
  strcat(text,": ");
  outtext(text);
  for(;;)
    {
      c = (char)readkey();
      if (is_crlf(c)) {	outtext("\n"); return(-1); }
      if (c<'0' || c>'9') continue;
      i = (int)(c-'0');
      if (i<nvar) break;
    }
  crlf();

  for(j=0; j<n && allowed[j]!=i; j++);
  if (j==n) goto AGAIN;

  m = (varid+i)->index;
  if (which==2) m &= INDBITS;
  return(m);
}

/*-----------------------------------------------------------------------------
|	readtext
-----------------------------------------------------------------------------*/
int readtext(char *prompt, char *format, char *s)
{
  int n, n0, is_int, col,ic;
  char *p, c;

  outtext(prompt);
  col = strlen(prompt) - 1;
  n0 = col - 1;
  is_int = (*(format+1)=='d');
  p = s; *p = 0;			/* p points to final NUL */
  for(n=0;;)
    {
      ic = readkey();
      c = *(char *)&ic;
      fflush(stdout);
      if (is_crlf(c)) break;
      if (c=='\b')
        {
	  if (n==0) continue;
	  col = n0 + n;
	  n--; p--; *p=0;
	  outtext("\b");
	  continue;
	}
      else if (n>=78) continue;
      else if (!is_int || c=='-');
      else if (c<'0' || c>'9') continue;
      n++; *p++ = c; *p = 0;
      col = n0 + n;
      outtext(p-1);
    }
  fflush(stdin);
  outtext("\n");
  return(p>s);				/* 0 if null string */
}

/*-----------------------------------------------------------------------------
|	readint -- returns #fields read, value in *val
-----------------------------------------------------------------------------*/
int readint(char *prompt, int *valp, int i1, int i2)
{
  char text[80];
  int n, val;
  for(;;)
    {
#ifdef DEAD_CODE
      n = scanf("%d", &val);
      if (n) sprintf(text,"%d",val);
#endif
      n = readtext(prompt, "%d", text);
      if (n==0) return(0);
      sscanf(text, "%d", &val);
      if (val<i1 || val>i2) continue;
      *valp = val;
      return(1);
    }
}

/*-----------------------------------------------------------------------------
|	newval
-----------------------------------------------------------------------------*/
int newval(int ind, CURVE_SET *cp)
{
  LOOP *lp;
  VARIABLE_ID *q;
  int j, val, *p;
  char text[80];

  p = &cp->param[ind];			/* cp gives old value */
  lp = loop + ind;
  val = *p;
  ind |= LOOPBIT;

/*------ Get variable name assoc. with 'ind' */
  for(j=0,q=varid; j<nvar && q->index != ind; j++,q++) ;
  sprintf(text,"  Enter %s (0..%d, currently %d): ",
	  q->name, lp->count - 1, val);
  readint(text, &val,0,lp->count-1);
  return(val);
}

/*-----------------------------------------------------------------------------
|	newrange
-----------------------------------------------------------------------------*/
int newrange(CURVE_SET *cp)
{
  LOOP *lp;
  char text[80];
  int val;
  lp = loop + cp->lfaml;
  sprintf(text,"Enter index of desired curve in family (0..%d): ",
	  lp->count-1);
  readint(text, &val, 0, lp->count-1);
  return(val);
}

/*-----------------------------------------------------------------------------
|	change_curve_parameter
|	* code: X,Y,F,O,1,2
-----------------------------------------------------------------------------*/
void change_curve_parameter(char code, CURVE_SET *newcs,
			    int *newx, int *newy)
{
  int ind, j;
  if (code=='X' && (ind=newcode(3)) != -1) {
    *newx=1; newcs->ix.index = ind;
  }
  else if (code=='Y' && (ind=newcode(3)) != -1) {
    *newy=1; newcs->iy.index = ind;
  }
  else if (code=='F' && (ind=newcode(2)) != -1) {
    newcs->lfaml = ind;
    }
  else if (code=='O' && (ind=newcode(2)) != -1 &&
	   (j=newval(ind, newcs)) != -1) {
    param[ind] = j;
  }
  else if (code=='1' && newcs->gtype=='G' &&
	   (ind=newrange(newcs)) != -1) {
    newcs->ncurve=ind; newcs->flags |= SINGLE;
  }
  else if (code=='2') newcs->flags &= ~SINGLE;

  getparam(newcs);
}

/*-----------------------------------------------------------------------------
|	editcurve
|	* mstring is used with "menu" version of xdraw, set from_menus=1;
-----------------------------------------------------------------------------*/
int editcurve(int iwin, CURVE_SET *cp, char *mstring)
{
  static CURVE_SET newcs;
  char code;
  int newx, newy, ok;
  static int inited=0;
  static int x=10, y=20;

#ifdef DOS
  if (!inited)
    {
      _settextwindow(20, 1, 34, 79);
      _settextposition(1,1);
      inited = 1;
    }
  settextcolor(TEXT0);
#endif

  from_menus = (mstring==0) ? 0 : 1;
  newcs = *cp;
  initparam(cp);				/* Load param[] */

  if (!from_menus) pmenu(1);			/* print menu */
  ptitle(iwin, &newcs);				/* print description */
  if (!editing)
    {
      xmessage(x,y,newtitle);
      editing = 1;
    }

  for(newx=newy=0;;)
    {
      ok = curvetest(&newcs);
      code = from_menus ? mstring[1] : pmenu(0);	/* returns ' ' for nul */
      if (code=='M')
	{ if (!from_menus) pmenu(1); continue; }
      else if (code=='G') break;
      else if (code=='Q')
	{
	  editing=0; xmessage(x,y,newtitle);
	  ok=0; break;
	}
      else change_curve_parameter(code, &newcs, &newx, &newy);

      xmessage(x,y,newtitle);
      ptitle(iwin, &newcs);			/* print description */
      xmessage(x,y,newtitle);
      if (from_menus) return(0);
    }

  if (!from_menus)
    { outtext("Done with edit mode");
      if (ok) {outtext(", enabling new parameters"); }
      outtext(".\n");
    }
  if (!ok) return(0);
  *cp = newcs;
  get_limits(cp);
  return(1);
}
/*=============================================================================
|	Functions for Label File I/O
=============================================================================*/
static lx1, ly1, nline, boxwidth;
static char *boxedtextp;
/*-----------------------------------------------------------------------------
|	init_label -- read and store label info from .in
-----------------------------------------------------------------------------*/
void init_label(CURVE_SET *cp, char *info, int boxed)
{
  struct LABEL *q, *q1;
  char *p, *p1, *p2, *pstart[4], align;
  int ncomma, nquote, n;

  q = (struct LABEL *)malloc(sizeof(struct LABEL));
  q1 = cp->label;
  if (q1==NULL) cp->label = q;
  else
    {
      for (; q1->next; q1 = q1->next) ;
      q1->next = q;
    }
  q->next = NULL;

  pstart[0] = info;
  for(p=info,ncomma=0,nquote=0; nquote < 2; p++)
    {
      if (*p==' ' && nquote<1)
	d_abort("Space found while parsing -l or -b option",0,0);
      if (*p=='"' && ncomma<3) break;
      if (*p==' ' && (ncomma<3 || nquote<1)) break;
      if (*p==',' && nquote<1)
	{
	  ncomma++; if (ncomma==3 && *(p+1)!='"') break;
	  pstart[ncomma] = p+1;
	}
      if (*p=='"') { nquote++; if (nquote==1) p1 = p+1;  else p2 = p; }
    }
  if (ncomma < 3)
    s_abort("Insufficient arguments found for -%s option", boxed?"b":"l",0,0);
  if (nquote<2)
    s_abort("Text must be in quotes for -%s option",boxed?"b":"l",0,0);
  if (!boxed)
    {
      sscanf(info, "%d", &q->color);
      if (q->color < LBL_WHITE)	q->color = LBL_WHITE;
    }
  else
    {
      sscanf(info, "%c", &align);
      q->color = LBL_ALIGNL;
      if      (align=='C') q->color = LBL_ALIGNC;
      else if (align=='R') q->color = LBL_ALIGNR;
    }
  sscanf(pstart[1], "%f", &q->x);
  sscanf(pstart[2], "%f", &q->y);
  q->boxed = boxed;
  n = p2 - p1;
  q->text = (char *)malloc(n+1);
  strncpy(q->text, p1, n);
  *(q->text + n) = 0;
}

/*-----------------------------------------------------------------------------
|	get_newline
-----------------------------------------------------------------------------*/
char *get_newline(char *s)
{
  char *p;
  for(p=s; p=strchr(p,'\\'); p++)
    {
      if (*(p+1)=='n') return p;
    }
  return s+strlen(s);
}

/*-----------------------------------------------------------------------------
|	get_label
|	* k = index number (family ind/color) to search for, -1 ==> white
|	* starts search at q, returns q->next
-----------------------------------------------------------------------------*/
struct LABEL *get_label(int k, struct LABEL *q, int *lx, int *ly,
	       char **text, int *nline, int *align,
	       float xscale, float xoffset, float yscale, float yoffset)
{
  float x,y;
  int width, nl;
  char *p, *p2, delim;
  
  *text = NULL;
  *nline = 1;
  for(; q; q = q->next)			/* find next label w/ k */
    {
      if (q->color==k || (k==LBL_WHITE && q->color<=k)) break;
    }
  
  if (q == NULL) return q;
  *align = q->color;
  
  x = xscale * q->x + xoffset; *lx = float_to_int(x);
  y = yscale * q->y + yoffset; *ly = float_to_int(y);
  *text = q->text;
  nl = 1;
  if (q->boxed)
    {
      for(p=q->text,nl=0,delim=1; delim ; p=p2+2, nl++)
	{
	  p2 = get_newline(p);
	  delim = *p2; *p2 = 0;
	  width = textwidth(p, strlen(p));
	  if (p==q->text || width > boxwidth)
	    {
	      boxwidth = width; boxedtextp = p;
	    }
	  *p2 = delim;
	}
    }
  *nline = nl;
  return q->next;
}

/*-----------------------------------------------------------------------------
|	get_labelbox
-----------------------------------------------------------------------------*/
void get_labelbox(int *widthp, char **s)
{
  *s = boxedtextp;
  *widthp = boxwidth;
}



