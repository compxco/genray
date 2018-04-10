/******************************************************************************
**	NAME		XEDIT.H
**	AUTHOR		Sheryl M. Glasser
**
**	DESCRIPTION
**		
**
**	Copyright (c) Toptools SCF 1990.  All rights reserved.
******************************************************************************/
extern int editcurve(int iwin, CURVE_SET *cp, char *mname2);
extern void print_title(int iwin, CURVE_SET *p);
extern char *get_descriptive_title(int iwin, CURVE_SET *p);
extern void init_label(CURVE_SET *, char *, int);
extern struct LABEL *get_label(int, struct LABEL *, int *, int *,
		      char **, int *, int *,  float,float,float,float);
extern void get_labelbox(int *, char **);
extern char *get_newline(char *s);
