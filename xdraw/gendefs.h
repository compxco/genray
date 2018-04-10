/******************************************************************************
**  NAME      GENDEFS.H
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      structures LOOP and CURVE_SET xdraw.c, xinit.c
**
**  Copyright (c) GlassWare 1994.  All rights reserved.
******************************************************************************/
#ifdef UNIX
#define Int int
#define Long int
typedef int byte;
#else
#define Int long
#define Long long
typedef unsigned char byte;
#endif

#define FZ (float)0
#define FH (float).5

#ifdef UNIX
#define addfloat(bx,ix0) bx+ix0;
#endif

#define float_to_int(fx1) (int) ( (fx1>=FZ) ? (fx1+FH) : (fx1-FH) )

#define dprintf printf
#define NOX (mydisplay==NULL)

#if defined(MOTIF)
#elif defined(MSWINDOWS)
#else
#define XCALLS
#ifndef UNIX
#define DOS
#endif
#endif

#if defined(DOS)
#define xprintf xPrintf		/* in xmslib.c */
#else
#ifndef XPRINTF			/* dialog fns in setcolor.c */
extern void xprintf(const char *, ...);
#else
extern void xprintf(const char *, int, int, int, int, int, int, int, int);
#endif				/* end XPRINTF */
#endif				/* end not DOS */

extern int xinput(char *, int *);

