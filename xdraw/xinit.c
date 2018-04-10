/******************************************************************************
**  NAME      XINIT.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Read in data from .in and .bin, load structures
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <X11/Xlib.h>
/*#include <X11/Xutil.h>
  #include <X11/Xatom.h>*/
#define O_BINARY 0
#define PATHSEP '/'
#define STATENAME "state.xdraw"

#else
#include <sys\types.h>
#include <sys\stat.h>
#include <io.h>
#include <conio.h>
#include <dos.h>
#include "xlib.h"
/*#include "xutil.h"*/
#define PATHSEP '\\'
#define STATENAME "state.xdr"
#endif

#ifndef O_RDONLY
#define O_RDONLY 0			/* for penning */
#endif

#include "gendefs.h"
#ifdef DOS
#include <graph.h>
#endif
#include "curves.h"
#include "xinit.h"
#include "ps.h"
#include "xcontour.h"
#include "xtools.h"
#include "menuwin.h"
#include "xdraw.h"
#include "xedit.h"		/* for init_label */

static void freadf(FILE *, char *, char *);
static void readstate(void);
static FILE *get_inname(int argc, char *argv[]);
static int parsevar(char *);
static void testvar(CURVE_SET *);
static void read_asc_header(FILE *frd);
static void use_asc_header(void);
static int var_exists(int index);
static void set_indexed_value(char *p, char *fname);
static char *set_subtitle(char *, CURVE_SET *, char, char *);

#define QUOTE '\"'

extern float *buf;
static unsigned Long bufsize;
static char namebuf[NCHAR];
static char newtitle[100];
static char newlabel[100];
static char basepath[100];
static char xpolar[20]="X",ypolar[20]="Y";

extern int ftype;	/* 0=graph, 1=contour, 2=self-determined, 3,4 */
extern int gtype;	/* graph in window: gGCQLBD etc */
static int version;	/* 0=Alan's draw.in, 1=Sheryl's */
extern int default_ncol, default_ncurve;
extern int pscolor, pstitle, splinetype, showcomment, figures;
extern int labeldy, labelsep;
extern float psfontfac;
extern int xhair_type;
extern int dialogwindow;

#ifndef UNIX
extern int dialog_col1;
#endif

static char datapath[100];
char prefix[20];
int zdebug = 0;			/* Only while zoom error - see xtools */
int ztestoption=0;
extern int debug_index, ldebug, verifying;

extern char title[];

extern LOOP loop[];		/* Loop structure of the data */
extern int nloop;
static int iloop, ih0=0;
static int counting_outer=0, i_outer=0;
int param[MAXLOOP];
int outerloop_added = 0;
int outerloop_count[] = { 0,0,0,0,0,0,0,0,0,0 };

static int use_hcount=0, stringbytes=0;
static char *dummytitle = " ";
extern NODE *nodelist;
extern int nnode, inode, ivar, ncount_equal;
#define NBLOCK 20

#define MAXVAR 128
VARIABLE_ID varid[MAXVAR];	/* Variables used for plotting */
int nvar=0;

extern CURVE_SET curveset[];	/* Describes the plots */
extern CURVE_SET *cp, *cp2;
extern int ncset;

#define is_white(c) (c==' ' || c=='\t')
#define is_eol(c)   (c=='\0' || c=='\n')
#define is_crlf(c)  (c=='\n' || c=='\r' || c==27)

#ifdef USEMAIN
/*-----------------------------------------------------------------------------
|	main (for debugging)
-----------------------------------------------------------------------------*/
void main(int argc, char *argv[])
{
  init(argc,argv);
}
#endif

/*-----------------------------------------------------------------------------
|	ztest
|	* testing = 0 via codeview --> use xtest[], ytest[] as zoom
|	* ztestoption: 1=print zoom coords, 0=no option,
|	  2 ==> array of xtest[], ytest[] (xw = char *)
-----------------------------------------------------------------------------*/
void ztest(int *xw, int *yw)
{
   static int testing=-1;	/* normal value -1 */
   static int nprint=0;
   char *p,*q;
   int i;

   static int xtest[12], ytest[12];
   if (ztestoption==1)
   {
      xprintf("%d %d   ", *xw, *yw);
      if (nprint == 1) xprintf("\n");
      nprint = nprint ? 0 : 1;
   }
   else if (ztestoption==2)
   {
      for(i=0, p=(char *)xw; q=skip_to_arg(p,2*i+2,0); i++)
        sscanf(q, "%d %d", xtest+i, ytest+i);
      *(xtest+i) = -1;
      testing = ztestoption = 0;
   }
   else if (testing>=0)
   {
      if (xtest[testing]==-1) testing=0;
      if (testing==0) xprintf("Warning! autozoom in effect\n");
      *xw = xtest[testing];
      *yw = ytest[testing];
      testing++;
   }
}

/*=============================================================================
**                  INITIALIZATION
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	parsevar
-----------------------------------------------------------------------------*/
static int parsevar(char *p)
{
  int mask, i;
  mask = 0;
  if (*p=='i' || *p=='I') { mask = LOOPBIT; p++; }
  sscanf(p,"%d", &i);
  i |= mask;
  return(i);
}

static int variable_ok(int index)
{
  int i, result;
  result = 1;
  i = index & ~LOOPBIT;
  if (index & LOOPBIT)
    {
      if (i < 0 || i >= nloop) result = 0;
    }
  else if (i < 0 || i > loop[ivar].count + loop[ivar].ih0) result = 0;
  return result;
}

static void testvar(CURVE_SET *cp)
{
  int ok1, ok2;
  ok1 = variable_ok(cp->ix.index);
  ok2 = variable_ok(cp->iy.index);
  if (ok1 && ok2) return;
  i_abort("Inappropriate variable found",
	  ok2 ? cp->ix.index : cp->iy.index, cp - curveset);
}

/*-----------------------------------------------------------------------------
|	parse_title, parse_subtitle
|	* returns parsed version of title or subtitle in text
|	* <i0> variables expanded, SINGLE expanded
-----------------------------------------------------------------------------*/
int parse_title(CURVE_SET *cp, char *title, char *text)
{
  char *p, *p2;
  int var, i;
  var = -1;
  p = strchr(title,'<');
  if (p) p2 = strchr(p,'>');
  if (p && p2) var = 0;
  if (var != -1)
  {
    *p2 = '\0';
    i = parsevar(p+1);
    *p2 = '>';
    if (!variable_ok(i) || ! (i & LOOPBIT)) var = -1;
    else var = i & ~LOOPBIT;
  }
  if (var == -1) strcpy(text, title);
  else
  {
    strncpy(text, title, p-title);
    sprintf(text+(p-title), "%d", cp->param[var]);
    strcat(text, p2+1);
  }
  return( (var==-1) ? 0 : 1);
}

int parse_subtitle(CURVE_SET *cp, char *string)
{
  if (!cp->subtitle) return 0;
  if ((cp->flags & SINGLE) && load_single_label(cp, string)) ;
  else parse_title(cp, cp->subtitle, string);
  return 1;
}

/*-----------------------------------------------------------------------------
|	axis_label
-----------------------------------------------------------------------------*/
static void axis_label(char *p, CVAR *xp)
{
  char text[20], *q;
  static int inited=0;
  if (!inited) { *newlabel='\0'; inited++; }
  for(q=text; !is_white(*p) && !is_eol(*p); ) *q++ = *p++;
  *q = '\0';
  for(q=newlabel; *q && strcmp(q, text); q += strlen(q) +1) ;
  if (!*q) { strcpy(q, text); *(q+strlen(q)+1) = '\0'; }
  xp->label = q;
}

/*-----------------------------------------------------------------------------
|    freadf
|	version=0: read a header line, then a value
|	version=1: read a label on current line to (:), then a value
|	both cases skip over blank lines
-----------------------------------------------------------------------------*/
static void freadf(FILE * frd, char *format, char *s)
{
  char text[80], *p;
  int i, nline;
  nline = version ? 1 : 2;	/* version 0 has to read a header */
  for (i = 0; i < nline;)	/* Read 1 or 2 lines, plus blank lines */
    {
      fgets(text, 80, frd);
      *strchr(text, '\n') = 0;
      if (*text != 0)
	i++;
    }
  if (version == 1)
    {
      if (!strcmp(format,"%l")) { strcpy(s,text); return; }
      p = strchr(text, ':');
      for (p++; isspace(*p); p++);
      if (!strcmp(format, "%s"))
	strcpy(s, p);
      else if (!strcmp(format, "%c"))
	*s = *p;
      else if (!strcmp(format, "%d"))
	sscanf(p, format, (int *) s);
      else if (!strcmp(format, "%f"))
	sscanf(p, format, (float *) s);
    }
  else
    {
      if (!strcmp(format, "%s"))
	strcpy(s, text);
      else
	sscanf(text, format, s);
    }
}

static int blankline(char *s)
{
  char *p;
  for(p=s; *p==' ' || *p=='\t'; p++);
  return (*p < ' ');
}

/*-----------------------------------------------------------------------------
|	readstate
-----------------------------------------------------------------------------*/
static void readstate()
{
  FILE *frd;
  char ddelim,*p,fname[LNAME];
  char *inipath, ininame[LNAME+100];
  int i;
#define code(s) !strcmp(fname, s)

  ddelim = PATHSEP;
  *datapath=0;				/* init variables */
  strcpy(prefix, "draw");

  p = ininame; *p = 0;
  inipath = getenv("XDRAWPATH");
  if (inipath)
    {
      strcat(ininame, inipath);
      p = ininame + strlen(ininame);
      if (*(p-1) != PATHSEP) { *p++ = PATHSEP; *p = 0; }
    }

  strcat(ininame, "xdraw.ini");
  frd = fopen(ininame, "rt");
  if (frd==NULL)
    {
      strcpy(p, STATENAME);
      if ((frd = fopen(STATENAME, "rt"))==NULL) return;
    }

/*-----xdraw.ini or state.xdraw exists */
  for (;;)
    {
      p = fgets(fname, LNAME, frd);		/* read a line */
      if (!p || !strncmp(fname,"END",3))
	break;
      if (*fname==';') continue;		/* comment line */
      if (!(p=strchr(fname,':'))) continue;	/* look for (:) */

      for (*p++=0; is_white(*p); p++);		/* skip white */
      if (is_eol(*p)) continue;

      if (code("path"))				/* "path" */
	{
	  strcpy(datapath, p);
	  p = datapath + strlen(datapath) - 1;
	  if (*p != ddelim)
	    {
	      *p++ = ddelim; *p = 0;
	    }
        }
      else if (code("prefix"))			/* "prefix" */
	sscanf(p, "%s", prefix);
      else if (code("crosshair"))		/* "crosshair" */
	sscanf(p, "%d", &xhair_type);
      else if (code("splinetype"))		/* "splinetype" */
	sscanf(p, "%d", &splinetype);
      else if (code("figures"))			/* # "figures" for labels */
	sscanf(p, "%d", &figures);
      else if (code("markerstyle"))		/* "markerstyle" */
	setmarkerstyle(p);
      else if (code("showcomment"))		/* "showcomment" */
	sscanf(p, "%d", &showcomment);
      else if (code("menu"))			/* "menu" */
        {
	  sscanf(p, "%d", &i);
	  set_menu_visibility(i);
        }
      else if (code("pscolor"))			/* "pscolor" */
	sscanf(p, "%d", &pscolor);
      else if (code("pstitle"))			/* "pstitle" */
	sscanf(p, "%d", &pstitle);	      
      else if (code("psfontsize"))		/* "psfontsize" */
        {
	  sscanf(p, "%f", &psfontfac);
	  if (psfontfac < (float).5) psfontfac = (float).5;
	  else if (psfontfac > (float)1.5) psfontfac = (float)1.5;
        }
      else if (code("pstest"))
        sscanf(p, "%d %d", &labeldy, &labelsep);

      else if (code("ncontour"))		/* "ncontour" */
	sscanf(p, "%d", &default_ncurve);
      else if (code("big"))			/* "big, was fillscreen" */
        { sscanf(p, "%d", &default_ncol);
	  default_ncol = default_ncol ? 2 : 3; }
      else if (code("xpolar"))			/* "xpolar" */
	sscanf(p, "%s", xpolar);
      else if (code("ypolar"))			/* "ypolar" */
	sscanf(p, "%s", ypolar);

      else if (code("debug"))			/* "debug" for limits */
	sscanf(p, "%d", &ldebug);
      else if (code("zdebug"))			/* "zdebug" for zoom */
	sscanf(p, "%d", &zdebug);
      else if (code("debugindex"))		/* "debugindex" */
	sscanf(p, "%d", &debug_index);		/* >=0: index to savelevel*/
      else if (code("verify"))			/* "verify" */
	sscanf(p, "%d", &verifying);
      else if (code("autozoom"))			/* "autozoom" */
      {
	sscanf(p, "%d", &ztestoption);
	if (ztestoption > 1) ztest((int *)p, NULL);
      }
#ifndef UNIX
      else if (code("dialogsize"))
      { sscanf(p, "%d", &dialog_col1); dialog_col1 = 80 - dialog_col1; }
#endif
      else if (code("dialog"))
	{ sscanf(p, "%d", &dialogwindow); if (dialogwindow) dialogwindow=1; }
    }
  fclose(frd);
}

/*-----------------------------------------------------------------------------
|	get_inname
|	* parse input arg, get filename
|	* read 1st line, get ftype
-----------------------------------------------------------------------------*/
static FILE *get_inname(int argc, char *argv[])
{
  int i,suffix;
  FILE *frd;
  char *p, filename[100], text[150];

  *basepath = 0;
  for (i = 1, ftype = suffix = 0; i < argc; i++)
    {
      p = argv[i];
      if (*p++ == '-')				/* ------ Dash options */
	{
	  if      (*p == 'c')  ftype = 1;
	  else if (*p == 'p')  ps_suffix(p + 1);
	  else if (*p == 'd')  dialogwindow = 1;
  	  else if (*p == 'b')  { strcpy(basepath, p + 1); strcat(basepath,"/"); }
	}
      else suffix = i;
    }
  strcpy(filename, datapath);
  if (*prefix=='*' || *prefix=='?')
    {  if (!suffix) strcat(filename, "draw"); }
  else strcat(filename, prefix);
  if (suffix)
    strcat(filename, argv[suffix]);
  strcat(filename, ".in");			/* ------ Reading drawx.in */
  if (!(frd = fopen(filename, "r")))
    s_abort("Cannot open %s", filename,0,0);

/*------ Read 1st line of draw.in */
  fgets(text, LNAME, frd);			/* "filename" or "Type" */
  version = (strncmp(text, "Type", 4)) ? 0 : 1;
  if (version == 1)
    {
      ftype = 0;
      p = strchr(text, ':');			/* get type G,C,I etc */
      if (p != NULL)
	for (p++; *p == ' ' || *p == '\t'; p++);
      for(; *p && *p!='\n'; p++)
	{
	  if (*p>='a' && *p<='z') *p -= ('a'-'A');
	  if      (*p == 'C') ftype=1;		/* Contour data format */
	  else if (*p == 'G') ftype=0;		/* original Graphics format */
	  else if (*p == 'I') ftype=2;		/* multi-dim (Indexed) */
	  else if (*p == 'H') ftype=3;		/* 'G' with Header in .in */
	  else if (*p == 'L') ftype=4;		/* 'G' with 1 header rec */
	  else s_abort("Undefined char in 1st line %s\n",filename,0,0);
	}
    }
  return(frd);
}

/*-----------------------------------------------------------------------------
|   skip_to_arg -- skip to the n-th argument in line text (1st arg: n=1)
|	           if endprompt non-zero (eg ':') skip it first
|		   returns NULL if no more
-----------------------------------------------------------------------------*/
char *skip_to_arg(char *text, int n, char endprompt)
{
  register char *p;
  int i, open_quote;

  if (endprompt) p = strchr(text, endprompt);	/* skip past endprompt */
  else p = NULL;
  if (!p) p = text;
  else p++;

  for (i = 0; i < n; i++)
    {
      open_quote = 0;
      for (; is_white(*p); p++);		/* skip leading spaces */
      if (is_eol(*p)) return (NULL);
      if (i == n - 1) return (p);
      for (; !is_white(*p) && !is_eol(*p); p++)		/* skip over arg */
	{
	  if (*p=='"') for(p++; *p!='"' && !is_eol(*p); p++) ;
	}
    }
  return (NULL);
}

/*-----------------------------------------------------------------------------
|	readline
-----------------------------------------------------------------------------*/
static int readline(FILE *frd, char *fname, int size);
static int readline(FILE *frd, char *fname, int size)
{
  char *p, *q;
  int leading_done;
  for(*fname = 0; ;)
    {
      p = fname + strlen(fname);
      if (fgets(p, p-fname+size, frd)==0) return 0;
      q = strchr(fname, '\n');
      if (q) *q=0;
      for(q=p, leading_done=0; *q; q++)
	{
	  if (!leading_done && (*q==' ' || *q=='\t')) continue;
	  leading_done = 1;
	  *p++ = *q;
	}
      *p = 0;
      if (p==fname) return 0;
      if (*(p-1) != '\\') return 1;
      *(p-1) = 0;
    }
}

/*-----------------------------------------------------------------------------
|	init
-----------------------------------------------------------------------------*/
void init(int argc, char *argv[])
{
#define LLNAME (LNAME*10)
  char fname[LLNAME], *nameptr, *p, *q, *qs;
  VARIABLE_ID *qv;
  CVAR *vp;
  int i,j;
  int more, next, already;
  int usecolor;
  char code, *p1, *p2, c;
  FILE *frd;

/*------ Read state.xdr or state.xdraw, open drawx.in */

  give_command_args(argc, argv);	/* save argc,argv for Motif */
  usecolor = 1;
  readstate();				/* Read state.xdr or state.xdraw */
  frd=get_inname(argc,argv);		/* get eg."draw1.in", read ftype */

  if (ftype==0)				/* Read ascii header if any */
    read_asc_header(frd);

  for(;;)				/* skip to "filename(s)" */
    {
      fgets(fname,LNAME,frd);
      if (strlen(fname) > 1) break;
    }

/*------ read binary file names and open input binary files */
  for (i=0,buf=NULL,bufsize=0;;i++)		/* read data file name(s) */
    {
      if (*basepath) strcpy(fname, basepath);	/* from -b */
      else strcpy(fname, datapath);
      p = fname + strlen(fname);
      fgets(p,LNAME,frd);
      if (*p == ';') continue;
      if (strlen(p) <= 1) break;
      *strchr(p, '\n') = 0;
      if (*basepath && (q=strchr(p,'/')))
      {
	 for(q++; strchr(q,'/'); q=strchr(q,'/')+1);
	 strcpy(p, q);
      }

      if (i==0) ps_dataname(p);
      
      buf = binread(fname, &bufsize);		/* read data (.bin file) */
      if (*basepath && i>0) continue;		/* -b: only read 1 file */
      if (buf == NULL)
      {
        printf("Abort!\n");
	exit(1);
      }
    }
  if (outerloop_added) use_asc_header();
  loop_structure();

/*------ Read rest of drawx.in */
  nameptr = namebuf;
  if (version == 0)				/* Read usecolor */
    freadf(frd, "%d", (char *) &usecolor);
  freadf(frd, "%s", title);			/* Read global title */
  for(p=title+strlen(title)+1,more=1;;more=0)
  {
    fgets(nameptr,LNAME,frd);
    if (blankline(nameptr)) break;
    if (more) { strcat(title, "\n"); p++; }
    strcpy(p, nameptr);
    p += strlen(p) +1;
  }
  *p = 0;

/*------ Read variable names ------*/
  fgets(nameptr,LNAME,frd);			/* Skip "variable names" */
  for (nvar=already=0,qv=varid; ;nvar++,qv++)
    {
      if (!already) fgets(nameptr, LNAME, frd);
      p1 = skip_to_arg(nameptr, 1, '\0');	/* 1st arg=0,1,2.. in seq */
      if (p1 == NULL) break;
      p2 = skip_to_arg(nameptr, 2, '\0');	/* 2nd arg=name, eg. x,y,theta*/
      qv->index = parsevar(p1);
      qv->name = p2;
      already = (p2=strchr(p2, '\\')) != NULL;
      if (!already)
        {
	  *(p = strchr(nameptr, '\n')) = 0;
	  nameptr = p + 1;
        }
      else
        {
	  *p2 = 0;
          nameptr = read_flabels(frd, p2+1, qv->index);
        }
    }

/*------ read info about each plot ------*/
  fgets(fname, LNAME, frd);			/* "ix iy" or "plot type" */
  for(ncset=more=0; ncset<NCS; ncset++)
    {
      cp=curveset+ncset;
      if (more==0)
	{
	  if (!readline(frd, fname, LLNAME)) break;
	  next = 2;
	}
      if (*fname==';')			/* skip comment line */
	{ ncset--; continue; }
      p = fname;

      if (*p>='0' && *p<='9') ;			/* Zap old type C file */
      else if (*p=='i' || *p=='I') ;
      else
        s_abort("Error specifying x-variable for graph spec:\n%s\n",
		fname,ncset,0);

      cp->title = cp->subtitle = NULL;
      cp->label = NULL;
      cp->window = 0;
      cp->flags = 0;
      cp->mcount = -1;
      cp->gtype = 'G';				/* if contour, need -c */
      cp->which = (char)(more? '2' : '1');
      cp->lstep = cp->lfaml = -1;
      cp->fskip = 1;
      for(i=0; i<nloop; cp->param[i++]=0) ;
      cp->iz.index = -1;
      cp->forcemin = cp->forcemax = (float)0;

      vp = &cp->ix;				/* Read ix */
      p=fname;
      vp->index = parsevar(p);
      vp->hloop = v_is_header(vp->index);

      p = skip_to_arg(fname,next,0);		/* Read iy */
      for(q=p+1; !is_white(*q) && !is_eol(*q); q++) ;
      if (*(q-1)=='.')      cp->flags |= DOTS;
      else if (*(q-1)=='&') cp->flags |= SPLINE;
      if (cp->flags & (DOTS|SPLINE)) *(q-1)=' ';
      vp = &cp->iy;
      vp->index = parsevar(p);
      vp->hloop = v_is_header(vp->index);
      testvar(cp);
      for(j=next+1,more=0;;j++)			/* If more iy's,get... */
        {					/* ...more,next for next pass */
	  p = skip_to_arg(fname,j,0);
	  if (p==NULL || *p=='-' || *p==QUOTE) break;
	  if (*p=='i' || *p=='I');
	  else if (*p>='0' && *p<='9');
	  else break;
	  if (!more) { more=1; next=j; }
	}

      initparam(NULL);				/* param[i]=0 */
      if (p && *p=='-')				/*........ any dash options? */
	for(;;j++)				/* get next dash option */
        {
	  p = skip_to_arg(fname,j,0);
	  if (p==NULL || *p != '-') break;
	  p++;					/* p on 1st char after dash*/
	  code = *p++;				/* p on 1st char after code*/
	  if	  (code=='s') cp->lstep = parsevar(p) & INDBITS;
	  else if (code=='f') cp->lfaml = parsevar(p) & INDBITS;
	  else if (code=='j') sscanf(p,"%d", &cp->fskip);
	  else if (code=='a') cp->flags |= ASPECT;
	  else if (code=='C') cp->flags |= COLORS;
	  else if (code=='p') cp->flags |= POLAR;
	  else if (code=='M')
	    { if (isdigit(*p)) sscanf(p,"%d", &cp->mcount); else cp->mcount=0;}
	  else if (code=='x') axis_label(p, &cp->ix);
	  else if (code=='y') axis_label(p, &cp->iy);
	  else if (code=='z') cp->iz.index = parsevar(p) & INDBITS;
	  else if (code=='t' || code=='T')
	    nameptr = set_subtitle(p, cp, code, nameptr);
	  else if (code=='X') { sscanf(p,"%f",&cp->force); cp->flags|=FORCED; }
	  else if (code=='E')		/* -E<min>,<max>, can have either... */
	    {				/* ...absent, or max < min */
	      q=strchr(p, ',');
	      for(qs=p; *qs && *qs>' '; qs++) ;		/* qs=1st white */
	      if (q && qs>q)
	        {
		  if (q > p  ) sscanf(p,  "%f", &cp->forcemin);
		  if (qs> q+1) sscanf(q+1,"%f", &cp->forcemax);
	        }
	    }
	  else if (code=='e')   sscanf(p,"%d",&cp->hzcurve);
	  else if (code=='L')   sscanf(p,"%d",&cp->ncurve);
	  else if (code=='F') cp->flags |= FILL;
	  else if (code=='A') { sscanf(p,"%d",&cp->ncurve);/* timestep here */
				cp->flags |= ANIMATE; }
	  else if (code=='#') { sscanf(p,"%d",&cp->ncurve);
				cp->flags |= SINGLE; }
	  else if (code=='r') cp->flags |= REVERSE;
	  else if (code=='c')
	    {
	      cp->gtype = 'C';
	      c = toupper(*p);
	      if (c=='S' || c=='D' || c=='L') cp->gtype = *p++;
	      if (cp->iz.index < 0) cp->iz.index = 0;
	    }

	  else if (code=='i') set_indexed_value(p, fname);
	  else if (code=='l')			/* labels on curves */
	    {
	      cp->flags |= LABELF;
	      if (*p > ' ') init_label(cp, p, 0);
	    }
	  else if (code=='b' && *p>' ')		/* boxed labels */
	    {
	      cp->flags |= LABELF;
	      init_label(cp, p, 1);
	    }
	}

      getparam(cp);			/* choose lstep, set cp->param[] */
      if (!curvetest(cp))		/* test if all OK */
	d_abort("Illegal combination found in curve %d",cp-curveset,0);

      if (p /*&& cp->which=='1'*/)		/* Read title, if any */
	{
	  if (*p == QUOTE) p++;			/* skip quote char */
	  strcpy(nameptr, p);
	  p = nameptr;
	  nameptr = strchr(nameptr, 0) + 1;
	}
      else if (cp->title != NULL)
	p = cp->title;
      else if (title != NULL)			/* no explicit title, use.. */
	p = title;				/* .. general label */
      else p = dummytitle;
	       
      /*if (cp->which=='1')*/ cp->title = p;

      if (cp->flags & POLAR)
	{
	  if (cp->ix.label==NULL) cp->ix.label = xpolar;
	  if (cp->iy.label==NULL) cp->iy.label = ypolar;
	}
      get_limits(cp);
    }
  cp2 = cp;
  ncset = cp2 - curveset;
}

/*-----------------------------------------------------------------------------
|	set_subtitle
-----------------------------------------------------------------------------*/
static char *set_subtitle(char *p, CURVE_SET *cp, char code, char *nameptr)
{
  char *q, delim;
  int n;
  if (*p==' ' || *p=='\t') return p;
  if (*p=='"') delim = *p++;
  else delim = '\0';
  for(q=p++;; p++)
    {
      if (*p==delim || *p=='\n') break;
      else if (!delim && (*p==' ' || *p=='\t')) break;
    }
  n = p - q;
  strncpy(nameptr, q, n);
  *(nameptr + n) = '\0';
  if (code=='t') cp->subtitle = nameptr;
  else cp->title = nameptr;
  p = nameptr + n + 1;
  return p;
}

/*-----------------------------------------------------------------------------
|	set_indexed_value
-----------------------------------------------------------------------------*/
static void set_indexed_value(char *p, char *fname)
{
  char *q, delim, text[20];
  int index, maxind;

  q = strchr(p,':');
  if (q==NULL) q=strchr(p,'=');
  if (q==NULL) d_abort("-i option in curve %d must include a (:) or (=)",
		                    cp - curveset,0);
  delim = *q; *q = 0;
  sscanf(p, "%d", &index);
  sprintf(text, "%s%c", p-2, delim);
  *q = delim;
  maxind = (inode==-1) ? nloop-1 : inode-1;
  if (index<0 || index>maxind)
    {
      sprintf(fname,"%s option in curve %d uses bad index %d (max value %d)",
	      text,cp - curveset,index,maxind);
      s_abort("%s",fname,0,0);
    }
  sscanf(q+1, "%d", &param[index]);
  if (param[index] >= loop[index].count)
    {
      sprintf(fname,"Illegal value %s%d option in curve %d (max value %s%d)",
	      text, param[index],cp-curveset,text,loop[index].count-1);
      s_abort("%s",fname,0,0);
    }
}

/*-----------------------------------------------------------------------------
|	v_is_header
-----------------------------------------------------------------------------*/
int v_is_header(int index)
{
  int i,j;
  LOOP *lp;
  if (index & LOOPBIT) return(-1);
  for(i=0,lp=loop; i<nloop; i++,lp++)
    {
      j = index - (int)lp->ih0;
      if (j < 0 ||  j >= (int)lp->hcount) continue;
      return(i);
    }
  return(-1);
}

/*-----------------------------------------------------------------------------
|	use_asc_header
-----------------------------------------------------------------------------*/
static void use_asc_header()
{
  LOOP *lp;
  int i, n, total;
  for(i=total=1; i<=outerloop_added; i++)
  {
    n = loop[i].count = outerloop_count[i];
    total *= n;
  }
  loop[0].count /= n;
}

/*-----------------------------------------------------------------------------
|	read_asc_header
-----------------------------------------------------------------------------*/
static void read_asc_header(FILE *frd)
{
  char text[150], *p, *q;
  int index, value, i;
  for(;;)
    {
      fgets(text, 150, frd);
      p = strchr(text, ':');
      if (p==NULL) return;
      if (!strncmp(text, "Outer", 5))
        {
          sscanf(p+1, "%d", &outerloop_added);
	  for(i=0; i<=outerloop_added; i++)
	    outerloop_count[i] = i ? 1 : 0;
        }
      else if (!strncmp(text, "Loop", 4) && outerloop_added)
        {
	  for(p++; p && (q=strchr(p,'='));)
	    {
	      *q = ' ';
	      sscanf(p, "%d %d", &index, &value);
	      p = strchr(q+1, ' ');
	      if (index <= outerloop_added && value > 0)
		outerloop_count[index] = value;
	    }
	}
    }
}

/*=============================================================================
**			Edit
**===========================================================================*/
#ifndef DOS
#define outtext(s) printf(s); fflush(stdout);
void settextcolor(int i) { i; }

#else
#define outtext _outtext
#define settextcolor _settextcolor
#endif

#define TEXT0 2
#define TEXT1 8
#define TEXT2 14

/*-----------------------------------------------------------------------------
|	curvetest
-----------------------------------------------------------------------------*/
int curvetest(CURVE_SET *cp)
{
  int indexf, error, i;
  CVAR *xp, *yp, *zp;
  char text[100];

  settextcolor(TEXT2);

/*------ Error test */
  error = 0;
  xp = &cp->ix;
  yp = &cp->iy;
  zp = &cp->iz;
  if (cp->gtype != 'G')
    {
      if (!(xp->index & LOOPBIT) || !(yp->index & LOOPBIT))
	error |= 0x10;
      if (zp->index & LOOPBIT) error |= 20;
    }
  else
    {
      indexf = cp->lfaml | LOOPBIT;
      if (indexf==xp->index) error |= 1;
      if (indexf==yp->index) error |= 2;
      if ((xp->index & LOOPBIT) && (yp->index & LOOPBIT)) error |= 4;
    }

  for(i=1; i<0x40; i*=2)
    {
      printvar(0, NULL, text);			/* Init text */
      if (!(error & i)) continue;
      printvar(-1, "WARNING! ", text);		/* Append to text */
      if (i & 3)
	{
	  printvar(indexf, "%s is used for both family and ", text);
	  if (i==1) printvar(-1, "X.\n", text);
	  else if (i==2) printvar(-1, "Y.\n", text);
	}
      else if (i==4)
	{
	  printvar(yp->index, "Y and X (%s and ", text);
	  printvar(xp->index, "%s) are both index variables.\n", text);
	}
      else if (i==0x10)
	strcat(text,"Contour x and y must both be index variables.\n");
      else if (i==0x20)
	strcat(text,"Contour must have a non-index z specified.\n");
      else continue;
      outtext(text);
    }

  settextcolor(TEXT0);
  return(error==0);
}

/*-----------------------------------------------------------------------------
|	printvar -- create a section of newtitle[]
|	using label assoc. w/ 'index'.  format==NULL ==> initialize
-----------------------------------------------------------------------------*/
int printvar(int index, char *format, char *newtitle)
{
  int i;
  VARIABLE_ID *q;
  char s[10],name[30];
  char *p;

  if (format==NULL) { *newtitle=0; return(0); }
  p = newtitle + strlen(newtitle);
  if (index==-1) { strcpy(p,format); return(0); }

  if (index & LOOPBIT) sprintf(s,"i%d",index & INDBITS);
  else sprintf(s,"%d",index);

  for(i=0,q=varid; i<nvar && q->index!=index; i++,q++) ;
  if (i==nvar) sprintf(name,"<%s>",s);
  else strcpy(name, (varid+i)->name);

  sprintf(p, format, name);
  return(i);
}

