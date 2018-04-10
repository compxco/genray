/******************************************************************************
**  NAME        PSPACK
**  AUTHOR      Sheryl M. Glasser
**
**  DESCRIPTION
**      Usage: pspack [-r<nr>] [-c<nc>] [output file] [input file....input file]
**      Packs input postscript files into a combined output file
**	r,c options specify desired # rows, columns to lay out the page in
**	If unspecified, defaults are given in table layout[] in ps.c
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>

#else
#include <sys\types.h>
#include <sys\stat.h>
#include <io.h>
#endif

#include "pspack.h"

#define FZ (float)0

FILE *file, *frd;
int wdx, wdy, wbord;
extern int ps_aspect;
char outfile[40] = "";

/*-----------------------------------------------------------------------------
|   main
-----------------------------------------------------------------------------*/
void main(int argc, char **argv)
{
  int n, i, i0, nrows, ncols, lbord;
  float xborder;
  char fname[40], names[200], text[100];
  char *namesp[25];		/* max pictures, actually 9 */
  char *p;

  nrows = ncols = 0;
  xborder = (float)1;
  for(i=1; i<argc; i++)
    {
      p = *(argv+i);
      if (*p++ != '-') break;
      if (*p == 'h')
	{
	  printf("Usage:   pspack [options] [output_file [input_file(s)]]\n");
	  printf("Options:\n");
	  printf("   -r<n> specifies a layout having n rows.  e.g. -r4\n");
	  printf("   -c<n> specifies a layout having n columns.  e.g. -r3\n");
	  printf("   -b<x> specifies border width of x inches.  e.g. -b.5\n");
	  printf("   -w    make graphs wide as appropriate for the layout.\n");
	  printf("   -a    use A4 paper\n");
	  printf("   -h    help\n");
	  exit(0);
	}
      else if (*p == 'r') sscanf(p+1, "%d", &nrows);
      else if (*p == 'c') sscanf(p+1, "%d", &ncols);
      else if (*p == 'b') sscanf(p+1, "%f", &xborder);
      else if (*p == 'a') ps_country(0);
      else if (*p == 'w') ps_aspect = 0;
    }

  i0 = i;
  if (i0 == argc)
    input("Enter output file name [pxdraw.ps]: ", outfile);
  else
    strcpy(outfile, argv[i0]);

  if (!strlen(outfile))
    strcpy(outfile, "pxdraw.ps");
  if (!strchr(outfile, '.'))
    strcat(outfile, ".ps");

  i0++;
  if (argc > i0) n = argc - i0;
  else if (nrows > 0 && ncols > 0) n = nrows * ncols;
  else if (nrows > 0) n = nrows * nrows;
  else if (ncols > 0) n = ncols * ncols;
  else n = get_maxps();
  for (i = 0, p = names; i < n; i++)
    {
      if (argc > i0) namesp[i] = argv[i + i0];
      else
	{
	  input("Enter input file name: ", p);
	  if (!strlen(p)) break;
	  namesp[i] = p;
	  p += strlen(p) + 1;
	}
    }
  if ((n = i) == 0) exit(0);
  if (nrows == 0 && ncols > 0) nrows = n;
  else if (ncols == 0 && nrows > 0) ncols = n;
  if (nrows > 0 && ncols > 0) ps_layout(n, nrows, ncols);
  lbord = (int)(xborder * (float)72);
  
  file = fopen(outfile, "wt");
  fprintf(file, "%c! Postscript packed file, %d subfiles\n\n", '%', n);
  for (i = 0; i < n; i++)
    {
      strcpy(fname, namesp[i]);
      if (*fname == '-') { printf("Blank entry\n"); continue; }
      if (strcmp(fname + strlen(fname) - 3, ".ps"))
	strcat(fname, ".ps");
      printf(fname);
      if (!(frd = fopen(fname, "rt")))
	{
	  printf(" not found.\n");
	  continue;
	}
      else
	printf("\n");
      readsize();				/* get wdx, wdy, wbord */
      writetrf(frd, file, i, n, wdx, wdy, wbord, lbord);
      for (; fgets(text, 100, frd);)
	{
	  if (i < n - 1 && !strncmp(text, "showpage", 8))
	    continue;
	  fputs(text, file);
	}
      fprintf(file, "grestore\n\n");
      fclose(frd);
    }
  fclose(file);
  printf("Output file %s\n", outfile);
}

/*-----------------------------------------------------------------------------
|   readsize -- open file, read header lines, scan for wdx, wdy, border
|	As set in xdraw, wbord = 72, wdx = 468 = 8.5*72-2*wbord
-----------------------------------------------------------------------------*/
void readsize()
{
  char text[100];
  char *p;
#define readsize(xp) p=strchr(p,'=')+1; sscanf(p,"%d",xp)

  for (;;)
    {
      fgets(text, 100, frd);
      if (!strncmp(text, "%!", 2)) continue;		/* 1st line has %! */
      if (strncmp(text, "%%", 2)) return;		/* next 3 lines have %% */
      for (p = text + 2; *p == ' '; p++);		/* skip leading spaces */
      if (strncmp(p, "Size", 2)) continue;		/* get size */
      readsize(&wdx);
      readsize(&wdy);
      readsize(&wbord);
      return;
    }
}

/*-----------------------------------------------------------------------------
|	input
-----------------------------------------------------------------------------*/
void input(char *prompt, char *p)
{
  char *q;
  printf(prompt);
  gets(p);
  if (q = strchr(p, '\n'))
    *q = 0;
}

/*-----------------------------------------------------------------------------
|	redraw (dummy - ps.c uses)
-----------------------------------------------------------------------------*/
void redraw() {}
void xPrintf() {}
