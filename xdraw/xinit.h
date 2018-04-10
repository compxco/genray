/******************************************************************************
**	NAME		XINIT.H
**	AUTHOR		Sheryl M. Glasser
**
**	DESCRIPTION
**		
**
**	Copyright (c) Toptools SCF 1990.  All rights reserved.
******************************************************************************/
extern void get_limits(CURVE_SET *cp);
extern int parse_title(CURVE_SET *cp, char *title, char *text);
extern int parse_subtitle(CURVE_SET *, char *);

#ifndef XDRAW_SOURCE
extern void limits(float *, CURVE_SET *, CVAR *);
extern void initparam(CURVE_SET *cp);
extern void getparam(CURVE_SET *cp);
extern int printvar(int index, char *format, char *newtitle);
extern int curvetest(CURVE_SET *cp);
extern char *read_flabels(FILE *frd, char *text, int index);
extern float *binread(char *fname, unsigned Long *sizep);
extern void loop_structure(void);
extern void d_abort(char *, int, int);
extern void s_abort(char *, char *, int,int);
extern void i_abort(char *, int, int);
extern int get_ncount(void);
extern Long step_over_ncount(int, int);
extern int asc_header_found(FILE *, int *);
extern int v_is_header(int);
extern char *skip_to_arg(char *, int, char);
extern void ztest(int *, int *);
#endif
