/******************************************************************************
**  NAME      CURVES.H
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      structures LOOP and CURVE_SET xdraw.c, xinit.c
**
**  Copyright (c) GlassWare 1994.  All rights reserved.
******************************************************************************/

#ifdef UNIX
#define in1 int
#else
#define in1 unsigned char	/* (want LOOP = 8 ints long) */
#endif

typedef struct {		/*----- Loop structure */
    char ltype;			/* I,N,X,A,V,H */
    char use_sep;		/* if not '1', use nnode[] instead of count */
    in1 hcount;			/* # variables as header to loop */
    in1 ih0;			/* index of first header variable (.in uses) */
    int count;			/* # non-header elements in loop */
    long sep;			/* count to next index at this level */    
    int current;		/* temp: current (fixed) value */
    float *extra;
    char *labels;		/* individual family member names from .in */
    } LOOP;

typedef struct {		/*----- struct in CURVE_SET for var */
  int index;			/* index from draw.in; if 'i', add 0x1000 */
  int hloop;			/* it's a var in a header; index of loop */
  long off0, dstep, dfaml;
  float min,max;
  float avg,rms;
  char *label;
} CVAR;

typedef struct {
  int ncount;
  long nsep;
} NODE;

struct LABEL {
  int color, boxed;
  float x,y;
  char *text;
  struct LABEL *next;	/* We REALLY want (struct LABEL *), but no compile */
  } ;

#define LBL_WHITE -1
#define LBL_ALIGNL -2
#define LBL_ALIGNC -3
#define LBL_ALIGNR -4

#define MAXLOOP 8
typedef struct {		/* Info about what's drawn in each window */
    char gtype;			/* N,G,g,C,S,D,B,K */
    char which;			/* '1' = first in the window, else 2nd, etc. */
    Window window;
    int flags;			/* MARKERS,COLORS,ASPECT,FILL, etc */
    int lstep;			/* index, offset in loop[], for step */
    int lfaml;			/* index, offset in loop[], for family */
    int fskip;
    int param[MAXLOOP];		/* >=0 ==>fixed value, -1=step, -2=family */
    CVAR ix,iy;
    CVAR iz;			/* for contours only, corresponds to lfaml */
    char *title;
    char *subtitle;
    struct LABEL *label;
    int mcount;			/* -1:markers all +; 0:cycle; n>0:n no line*/
    float force;		/* force one specific value for contour */
    float forcemin;		/* G: forcemin, forcemax used for ymin, ymax */
    float forcemax;
    int ncurve, hzcurve;	/* c:#contour lines; G:1st,2nd family members */
  } CURVE_SET;

typedef struct {	/*----- struct for each variable */
  int index;		/* index assigned in draw.in; if 'i', add 0x1000 */
  char *name;
} VARIABLE_ID;

/* CURVE_SET
 * type: N=unused, G=Graph, C=Contour,
 *       S=Solid, D=Dots on contour, B=Big
 */
#define MARKERS 0x1		/* flag bits */
#define COLORS	0x2		/* Combined graphs in window change color */
#define ASPECT	0x4		/* Graphs in window use data aspect ratio */
#define DOTS	0x8
#define POLAR   0x10
#define FORCED	0x20
#define FILL    0x40
#define ANIMATE 0x80
#define REVERSE 0x100
#define SINGLE  0x200
#define LABELF	0x400
#define SPLINE	0x800
#define E_VIS	0x1000		/* Extrema from visible family members only */
#define E_IGNORE 0x2000		/* Extrema: use / ignore -E option if any */

#define NNAME 50
#define LNAME 80
#define MAXWINDOW 9
#define NCHAR (LNAME*(NNAME+MAXWINDOW+1))	/* size of namebuf[] */
#define NCS 1024					/* size of curveset[] */

#define LOOPBIT 0x1000
#define INDBITS 0x0FFF

