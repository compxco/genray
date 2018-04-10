/******************************************************************************
**  NAME            PS.C
**  AUTHOR          Sheryl M. Glasser
**
**  DESCRIPTION
**
**
**  Copyright (c) GlassWare 1992.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ps.h"
#include "gendefs.h"

int ps_modeon = 0;
int pscolor = 0, pstitle = 1, ps_aspect = 1;
static float psfontsize = (float)1.75;
float psfontfac = (float)1;
int verifying = 0;
static int verify_count = 0;
static int ps_initialized = 0;
static char showstring[100];

static char *ps_fontname, ps_filename[20], *title;
static FILE *ps_file;

static int pagenumber = 0;
static int wx0, wy0, wdx, wdy;
static int oldwin_dy;
static float iborder = (float)1;	/* desired border width, in inches */
static int bsep, tsep, lsep;		/* see description for ps_scale() */

static int ticklabel = 32, expheight = 24, axistitle = 48;
static int titleheight = 64;
static float paperwidth = (float)8.5, paperheight = (float)11;
static int font_height;
static int afh, lfh, tfh;
static char suf = 'a';
static float Scale;		/* here for debug output only */
static int dataname[80];

struct LAYOUT
  {
    int npic;
    char orient;
    int rows, cols;
  };
struct LAYOUT layout[] =
{
  1, 'L', 1, 1,   2, 'P', 2, 1,
  3, 'L', 2, 2,   4, 'L', 2, 2,
  5, 'P', 3, 2,   6, 'P', 3, 2,
  7, 'L', 3, 3,   8, 'L', 3, 3,
  9, 'L', 3, 3
};

#define NLAY (sizeof(layout)/sizeof(struct LAYOUT))

#define LABELDY 28	/* Empirical! MSU was 22, 4. stringwidth gave 0! */
#define LABELSEP 16
int labeldy = LABELDY, labelsep = LABELSEP;
static int labeldeltay;
#define scalelabel(x) (int)((float)(x) * psfontfac * (float)2)

/*-----------------------------------------------------------------------------
|   postscript -- initiate a redraw() in postscript mode
-----------------------------------------------------------------------------*/
extern void redraw(void);
void postscript(char *gtitle)
{
  char *p;
  p = (char *)dataname;		/* penning had problems with ps_init */
  title = gtitle;
  if (!pstitle) titleheight = 0;
  ps_init("default", "Times-Roman", p);
  redraw();
  ps_showpage();
  ps_close();
}

/*-----------------------------------------------------------------------------
|	ps_dataname
-----------------------------------------------------------------------------*/
void ps_dataname(char *s)
{
  char *p = (char *)&dataname[0];
  strcpy(p,s);
}

/*-----------------------------------------------------------------------------
|   ps_suffix
-----------------------------------------------------------------------------*/
void ps_suffix(char *p)
{
  char c;
  c = *p;
  if (c >= 'A' && c <= 'Z')
    c += ('a' - 'A');
  if ((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9'))
    suf = c;
}


/*-----------------------------------------------------------------------------
|   ps_init
-----------------------------------------------------------------------------*/
void ps_init(char *filename, char *fontname, char *data)
{
  struct tm *ps_time;
  time_t ps_clock;
  char text[80], *p;

  ps_modeon = 1;
  if (ps_initialized) return;

  labeldeltay = labeldy + labelsep;
  labeldy = scalelabel(labeldy);
  labelsep = scalelabel(labelsep);
  labeldeltay = scalelabel(labeldeltay);

  ps_fontname = fontname;
  if (strcmp(ps_fontname, "default") == 0)
    ps_fontname = "Times-Roman";

  p = (strcmp(filename, "default") == 0) ? "xdraw" : filename;
  sprintf(ps_filename, "%s%c.ps", p, suf + (char) pagenumber);

  ps_file = fopen(ps_filename, "wt");
  if (ps_file == NULL)
    {
      ps_modeon = 0;
      exit(0);
    };

  tzset();
  time(&ps_clock);
  ps_time = localtime(&ps_clock);
  strcpy(text, asctime(ps_time));
  if (p = strchr(text, '\n'))
    *p = 0;

  fprintf(ps_file, "%s%s\n", "%!  PostScript/Source:...", data);
  fprintf(ps_file, "%s%s\n", "%%  Title:...............", title);
  fprintf(ps_file, "%s%s\n", "%%  CreationDate:........", text);

  pagenumber++;
  ps_initialized = 1;
}

/*-----------------------------------------------------------------------------
|   ps_adjustpage
|	new page: height = height of outermost window
|		  bx1...bheight = interior box
-----------------------------------------------------------------------------*/
void ps_adjustpage(int height, int bx1, int by1, int bwidth, int bheight)
{
  int xt, yt;
  static int outerbox=0, inited=0;

  if (!inited) expheight = ticklabel*3/4;
  ps_scale(height, bx1, by1, bwidth, bheight);

  ps_font(axistitle);
  ps_normalline();
  fprintf(ps_file, "\n%s%d\n\n", "%----Page: ", pagenumber);

/*------Draw Box Around Window, Center Title Above It */
  if (outerbox)
    {
      fprintf(ps_file, "\n%s\n", "%----Draw rectangle around window");
      fprintf(ps_file, "newpath\n");
      fprintf(ps_file, "%d %d moveto %d 0 rlineto\n0 %d rlineto %d 0 rlineto\n",
	      wx0, wy0, wdx, wdy, -wdx);
      fprintf(ps_file, "closepath stroke\n");
      xt = wx0 + wdx / 2;		/* xcenter if outer box drawn */
    }
  else
    xt = bx1 + bwidth / 2;		/* xcenter if outer box not drawn */
  yt = oldwin_dy - (wy0 + wdy);		/* use same ycenter --> better spacing */

  if (pstitle)
    {
      fprintf(ps_file, "gsave\n");
      ps_font(titleheight);
      ps_centerstring(xt, yt, 1, title, strlen(title));
      fprintf(ps_file, "grestore\n");
    }

  fprintf(ps_file, "\n%s\n", "/print1       % edit this in to see top of stack");
  fprintf(ps_file, "{ 7 string cvs show ( ) show } def\n\n");
}

/*-----------------------------------------------------------------------------
|	ps_scale -- properly scale and rotate the picture
|	height = dy of screen-window, eg. 1024/2. Incoming y in 0..height
|	bx1,by1,bwidth,bheight: define interior box
|
|	Transformations: x' = M x
|	For landscape, M = T(540,72) R(90) S(1.5) T(-wx0,-wy0)
|	Successive transformations are post-multiplied to create M.
|	So want the transformations to appear in the ps file as last one first.
|
|	bsep = separation between window and box, eg top or right
|	lsep = separation between eg. box and tick-labels
|	tsep = separation between tick-labels and axis-title
|
|	wdx = new width of entire window in screen units
|	    = bsep+axistitle+tsep+ticklabel+lsep+bwidth+bsep
|	scale*wdx = scale*(bwidth+afh+lfh) +tsep+lsep+bsep+bsep
|	scale = 6.5*72/wdx, for portrait with 1" borders,
|	2 equations, 2 unknowns (scale and wdx)
|
|	font size 96 points = 2.2 cm at scale factor 1.0
-----------------------------------------------------------------------------*/
void ps_scale(int height, int bx1, int by1, int bwidth, int bheight)
{
  int lborder, rborder, sep;
  int paper_g_and_t, screen_g_and_t;
  float scale, pgunits, a;
  int zdx,zdy,zbord;
  static int inited=0;
/*#define sizemodify(x) (x * (32 + 4 * psbigfont)) / 32*/
#define sizemodify(x) (int)((float)x * psfontsize * psfontfac)  

  oldwin_dy = height;
  lborder = (int) (iborder * 72);			/* 1" ==> 72 */
  rborder = (int) ((paperwidth - iborder) * 72);	/* 1" ==> 540 */
  pgunits = (float) (rborder - lborder);		/* 1" ==> 6.5*72 */

  lsep = bsep = 72 / 12;			/* some sizes, page units */
  tsep = 72 / 8;
  if (!inited)
    {
      ticklabel = sizemodify(ticklabel);
      axistitle = sizemodify(axistitle);
      titleheight = sizemodify(titleheight);
      expheight = sizemodify(expheight);
      inited++;
    }

  a = (float)(2.2*72./(96.*2.54));	/* empirical # gives text sizes */
  afh = (int)(a * (float) axistitle);	/* size of axis text, screen units */
  lfh = (int)(a * (float) ticklabel);
  tfh = (int)(a * (float) titleheight);

  paper_g_and_t  = (int)pgunits - (bsep+tsep+lsep+bsep);
  screen_g_and_t = bwidth + afh + lfh;
  scale = (float)paper_g_and_t / (float)screen_g_and_t;

  tsep /= scale;				/* sizes in screen units */
  lsep /= scale;
  bsep /= scale;

  sep = bsep + tsep + lsep + afh + lfh;		/* window to axis - screen u's */

  wdx = (int) (pgunits / scale);		/* scrn size page-margins */
  wdy = sep + bsep + bheight;
  wx0 = bx1 - sep;				/* LL of page - screen units */
  wy0 = height - by1 - bheight - sep;

  Scale = scale;
  
  zdx = (int)pgunits;
  zdy = (int)((wdy + tsep + tfh) * scale);
  zbord = lborder;
  fprintf(ps_file, "%s  Size: x=%d y=%d border=%d\n\n", "%%", zdx, zdy, zbord);

/*-------- create landscape centered transformation */
  writetrf(NULL, ps_file, 0,1, zdx,zdy,zbord, zbord); /* landscape, centered */

/*-------- create portrait, bottom of page transformation */  
  fprintf(ps_file, "%d %d translate\n", lborder, lborder);
  fprintf(ps_file, "%f %f scale\n", scale, scale);
  fprintf(ps_file, "%d %d translate\n", -wx0, -wy0);
}

void cprintf(char *s, int i)
{
  float x;
  x = i * Scale;		/* inches * 72 */
  x = x *2.54 / 72.;
  printf("%s=%f\n",s,x);
}

/*-----------------------------------------------------------------------------
|	ps_color
-----------------------------------------------------------------------------*/
void ps_color(float hue, float l, float s)
{
  if (!pscolor) return;
  hue += (float)30;
  if (hue > (float)360) hue -= (float)360;
  hue /= (float)360;
  if (l > (float)0) l = (float)1;
  fprintf(ps_file, "%f %f %f sethsbcolor\n", hue, s, l);
}

/*-----------------------------------------------------------------------------
|	ps_comment
-----------------------------------------------------------------------------*/
void ps_comment(char *s)
{
  fprintf(ps_file, "%s%s\n", "%----", s);
}

/*=============================================================================
**                  LINES
**===========================================================================*/
/*-----------------------------------------------------------------------------
|   ps_line
-----------------------------------------------------------------------------*/
void ps_line(int x1, int y1, int x2, int y2)
{
  y1 = oldwin_dy - y1;
  y2 = oldwin_dy - y2;
  fprintf(ps_file, "%d %d moveto  %d %d lineto stroke\n", x1, y1, x2, y2);
}

/*-----------------------------------------------------------------------------
|   ps_stroke
-----------------------------------------------------------------------------*/
void ps_stroke()
{
  fprintf(ps_file, "stroke\n");
  verify_count = 0;
}

/*-----------------------------------------------------------------------------
|   ps_moveto, ps_lineto
-----------------------------------------------------------------------------*/
void ps_moveto(int x1, int y1)
{
  y1 = oldwin_dy - y1;
  fprintf(ps_file, "%s\n", "%----Polyline");
  fprintf(ps_file, "%d %d moveto\n", x1, y1);
}

void ps_lineto(int x1, int y1)
{
  y1 = oldwin_dy - y1;
  fprintf(ps_file, "%d %d lineto\n", x1, y1);
}

/*-----------------------------------------------------------------------------
|   ps_fmoveto, ps_flineto
-----------------------------------------------------------------------------*/
void ps_fmoveto(float x1, float y1)
{
  y1 = (float) oldwin_dy - y1;
  fprintf(ps_file, "%s\n", "%----Polyline");
  fprintf(ps_file, "%g %g moveto\n", x1, y1);
}

void ps_flineto(float x1, float y1)
{
  if (!verifying)
    {
      y1 = (float) oldwin_dy - y1;
      fprintf(ps_file, "%g %g lineto\n", x1, y1);
    }
  else fprintf(ps_file, "%.3f %.3f data %d\n", x1, y1, verify_count++);
}

/*-----------------------------------------------------------------------------
|   ps_rectangle
-----------------------------------------------------------------------------*/
void ps_rectangle(int x1, int y1, int w, int h)
{
  int y2;
  y2 = oldwin_dy - y1;
  y1 = y2 - h;
  fprintf(ps_file, "%s", "%----Rectangle\n");
  fprintf(ps_file, "\
newpath \n\
%d %d moveto \n\
%d %d lineto \n\
%d %d lineto \n\
%d %d lineto \n\
closepath \n\
stroke \n", x1, y1, x1 + w, y1, x1 + w, y1 + h, x1, y1 + h);
}

/*-----------------------------------------------------------------------------
|   ps_clip
-----------------------------------------------------------------------------*/
void ps_clip(int x1, int y1, int w, int h)
{
  int y2;
  if (w == 0 && h == 0)
    {
      x1 = wx0;
      w = wdx;
      y1 = wy0;
      h = wdy;
      y2 = y1 + h;
    }
  else
    {
      y2 = oldwin_dy - y1;
      y1 = y2 - h;
    }

  fprintf(ps_file, "%s", "%----Clip begin\n");
  fprintf(ps_file, "\
initclip newpath\n\
%d %d moveto \n\
%d %d lineto \n\
%d %d lineto \n\
%d %d lineto \n\
closepath\n\
clip\n\
newpath\n", x1, y1, x1 + w, y1, x1 + w, y2, x1, y2);
}

/*-----------------------------------------------------------------------------
|   unclip
-----------------------------------------------------------------------------*/
void ps_unclip()
{
  fprintf(ps_file, "%s", "%----Clip done \ngrestore \n");
}

/*-----------------------------------------------------------------------------
|   ps_thinline, ps_normalline (0 = REALLY thin)
-----------------------------------------------------------------------------*/
void ps_thinline()
{
  int w;
  w = 2; /* was 3*/
  if (!pscolor && psfontfac<(float)1) w = 1;
  fprintf(ps_file, "%d setlinewidth\n", w);
}

void ps_normalline()
{
  fprintf(ps_file, "2 setlinewidth\n");
}

/*=============================================================================
**                  TEXT
**===========================================================================*/
/*-----------------------------------------------------------------------------
|   ps_font -- initialize font (findfont, scalefont, setfont)
-----------------------------------------------------------------------------*/
void ps_font(int newfont_height)
{
  font_height = newfont_height;
  fprintf(ps_file,
	  "/%s findfont %d scalefont setfont\n", ps_fontname, font_height);
}

/*-----------------------------------------------------------------------------
|   ps_setlabelfont
-----------------------------------------------------------------------------*/
void ps_setlabelfont(int which)
{
  ps_font(which ? axistitle : ticklabel);
}


/*-----------------------------------------------------------------------------
|   goodstring -- replace ( and ) with \( and \)
-----------------------------------------------------------------------------*/
static void goodstring(char *string, int stringlen)
{
  char *p, *q;
  int i;
  p = string;
  q = showstring;

  for (i = 0; i < stringlen; i++)
    {
      if (*p == '(' || *p == ')')
	*q++ = '\\';
      *q++ = *p++;
    }

  *q = 0;
}

/*-----------------------------------------------------------------------------
|   ps_draw_label
-----------------------------------------------------------------------------*/
void ps_draw_label(int lx1, int ly1, int lineno, char *s)
{
  ly1 = oldwin_dy - ly1;
  goodstring(s, strlen(s));
  fprintf(ps_file, "%d %d ", lx1, ly1);
  if (lineno==0) fprintf(ps_file, "%s\n","%----Begin boxed text");
  if (lineno > 0)
#ifndef LABELDY
    fprintf(ps_file, "(%s) stringwidth exch pop %d mul add ", showstring, lineno);
#else
    fprintf(ps_file, "%d %d mul sub ", labeldeltay, lineno);
#endif
  fprintf(ps_file, "moveto\n(%s) show\n", showstring);
}

/*-----------------------------------------------------------------------------
|   ps_draw_labelbox
-----------------------------------------------------------------------------*/
void ps_draw_labelbox(int lx1, int ly1, int nline, char *longtext)
{
  int boxheight, boxy0;
  ly1 = oldwin_dy - ly1;
  
  fprintf(ps_file, "%s", "%----Box around labels\n");
  /*boxheight = nline * labeldeltay + 3 * labelsep;
  boxy0 = labeldeltay + labelsep;*/
  boxheight = nline * labeldeltay + labelsep;
  boxy0 = labeldeltay;
  goodstring(longtext, strlen(longtext));
  fprintf(ps_file, "%d %d moveto ", lx1, ly1);
  fprintf(ps_file, "%d %d rmoveto\n", -labelsep, boxy0);
  fprintf(ps_file, "(%s) stringwidth\n", showstring);
#ifdef LABELDY
  fprintf(ps_file, "pop %d add %d ", 2*labelsep, boxheight);
#endif
							/* stack=h,w */
#ifdef DEAD_CODE
  fprintf(ps_file, "dup -%d exch rmoveto\n", labelsep);	/* to UL */
  fprintf(ps_file, "%d mul %d add\n", nline, labelsep);
  fprintf(ps_file, "exch %d add exch\n", 2*labelsep);
#endif  
							/* stack=total h,w */
  fprintf(ps_file, "dup neg 0 exch rlineto\n");		/* to LL */
  fprintf(ps_file, "exch dup 0 rlineto exch\n");	/* to LR */
  fprintf(ps_file, "0 exch rlineto\n");			/* to UR */
  fprintf(ps_file, "neg 0 rlineto stroke\n\n");
}

/*-----------------------------------------------------------------------------
|   ps_rshowstring -- draw string relative
-----------------------------------------------------------------------------*/
void ps_rshowstring(int xpos, int ypos, char *string, int stringlen)
{
  ypos = -ypos;
  goodstring(string, stringlen);
  if (xpos && ypos) fprintf(ps_file, "%d %d rmoveto\n", xpos, ypos);
  fprintf(ps_file, "(%s) show\n", showstring);
}

/*-----------------------------------------------------------------------------
|   ps_centerstring -- draw string centered
|       'above': 1 = main graph title, above;
|		-1 = tick-label,below; x,y = start of tick
|		 0 = axis label: above wy0
-----------------------------------------------------------------------------*/
void ps_centerstring(int xpos, int ypos, int above, char *string, int stringlen)
{
  ypos = oldwin_dy - ypos;
  if	  (above == 0)  ypos  = wy0 + bsep;		/* axis title */
  else if (above ==  1) ypos += tsep;			/* title */
  else if (above == -1) ypos -= (lsep + lfh);		/* tick label */
  goodstring(string, stringlen);
  fprintf(ps_file, "%s", "%----Center String\n");
  fprintf(ps_file, "%d (%s) stringwidth pop 2 div sub %d moveto\n",
	  xpos, showstring, ypos);
  fprintf(ps_file, "(%s) show\n", showstring);
}

/*-----------------------------------------------------------------------------
|   ps_centervertically -- draw vertical string, centered
|		right: 1 = axis title, -1 = tick labels
-----------------------------------------------------------------------------*/
void ps_centervertically(int xpos, int ypos, int right,
			 char *string, int stringlen)
{
  ypos = oldwin_dy - ypos;
  if (right == 1)
    xpos = wx0 + bsep + afh;
  else if (right == -1)
    xpos -= lsep;
  goodstring(string, stringlen);
  fprintf(ps_file, "%s", "%----Center String Vertically\n");
  fprintf(ps_file, "gsave\n%d %d translate\n90 rotate\n", xpos, ypos);
  fprintf(ps_file, "(%s) stringwidth pop 2 div neg 0 moveto\n",
	  showstring);
  fprintf(ps_file, "(%s) show\n", showstring);
  fprintf(ps_file, "grestore\n");
}

/*-----------------------------------------------------------------------------
|	ps_power
|	* x_or_y: 1=x, 2=y;  string = e.g. "x10-1+.46"
|	* Note: don't need goodstring(), since power has no ()
-----------------------------------------------------------------------------*/
void ps_power(int x_or_y, char *string)
{
  char *p, *pp;
  int lftx, ritx, topy, boty, up;
  lftx = wx0 + bsep;
  ritx = wx0 + wdx - bsep;
  boty = wy0 + bsep;
  topy = wy0 + wdy - bsep;

  fprintf(ps_file, "%s", "%----Power - Right or Top Justify\n");
  fprintf(ps_file, "gsave\n");
  if (x_or_y == 1)
    fprintf(ps_file, "%d %d translate\n", ritx, boty);
  else
    fprintf(ps_file, "%d %d translate 90 rotate\n", lftx + afh, topy);

  fprintf(ps_file, "(%s) stringwidth pop neg 0 moveto\n", string);
  p = string;
  pp = strchr(p, '+');
  if (!strncmp(p,"x10", 3))
    {
      fprintf(ps_file, "(x10) show\n");
      ps_font(expheight); p += 3;
      if (pp) *pp = 0;
      up = afh * 4 / 5;
      fprintf(ps_file, "0 %d rmoveto (%s) show ", up, p);
      if (pp)
        {
	  *pp = '+'; ps_font(ticklabel);
	  fprintf(ps_file, "0 %d rmoveto ", -up);
        }
    }
  if (pp) fprintf(ps_file, "(%s) show", pp);
  fprintf(ps_file, "\ngrestore\n");
}

/*=============================================================================
**                  OTHER
**===========================================================================*/
/*-----------------------------------------------------------------------------
|   ps_save
-----------------------------------------------------------------------------*/
void ps_save(int save)
{
  static int old_height;
  if (save)
    old_height = font_height;
  else
    font_height = old_height;
  fprintf(ps_file, "%s\n", save ? "gsave" : "grestore");
}

/*-----------------------------------------------------------------------------
|   ps_translate
-----------------------------------------------------------------------------*/
void ps_translate(int xpos, int ypos)
{
  fprintf(ps_file, "grestore \n %d %d translate\ngsave\n", xpos, ypos);
}

/*-----------------------------------------------------------------------------
|   ps_showpage
-----------------------------------------------------------------------------*/
void ps_showpage()
{
  static int pagenumber = 2;

  fprintf(ps_file, "showpage\n");
  fflush(ps_file);
  ps_modeon = 0;
/*    ps_0 (); */
}

/*-----------------------------------------------------------------------------
|   ps_close
-----------------------------------------------------------------------------*/
void ps_close()
{
  fclose(ps_file);
  ps_initialized = 0;
  ps_modeon = 0;
/*  xprintf("PostScript file %s written.\n", ps_filename);  */
}
/*=============================================================================
|	Functions called by pspack
=============================================================================*/
/*-----------------------------------------------------------------------------
|   writetrf -- called by xdraw AND pspack
|   * IN: e.g. wdx=(8.5-2)*72, wdy=equivalent page dy of pic, wbord=72
|	  i-th of n pictures
|   * pspack: wdx...wbord are values used in the original xdrawa.ps etc.
|   * layout: 2*border + (ncol-1)*sep + ncol * picdx = 8.5 * 72 (portrait)
|        or.. 2*border + (nrow-1)*sep + nrow * picdy = 11 * 72
-----------------------------------------------------------------------------*/
void writetrf(FILE *frd, FILE *fwr, int i, int n,
	      int wdx, int wdy, int wbord, int border)
{
  struct LAYOUT *p, *p2;
  float xinch, yinch, picdx, picdy, picdx2, picdy2;
  float xempty, yempty;
  int portrait, ok1, ok2, sep;
  int xsep, ysep, tx, ty, xborder, yborder;
  int row, col, nrow, ncol, dxpic, dypic;
  float xscale, yscale;
  char text[100];

  sep = (int)(.25 * 72);
  wdx = (int)( (float)wdx * paperwidth / 8.5);		/* USA-->European */
  wdy = (int)( (float)wdy * paperheight / 11.);

  for (p = layout, p2 = p + NLAY; p < p2 && p->npic != n; p++);
  portrait = (p->orient == 'P');
  nrow = p->rows;
  ncol = p->cols;
  xborder = yborder = border;

  xinch = (float) (portrait ? paperwidth : paperheight);
  yinch = (float) (portrait ? paperheight : paperwidth);
  picdx = (xinch * 72 - 2 * border - (ncol - 1) * sep) / ncol;
  picdy = (yinch * 72 - 2 * border - (nrow - 1) * sep) / nrow;
  picdx2 = picdy * (float) wdx / (float) wdy;
  picdy2 = picdx * (float) wdy / (float) wdx;
  ok1 = (picdy2 <= picdy);
  ok2 = (picdx2 <= picdx);
  if (ps_aspect)
    {
      if (ok1 && ok2)
        {
          if (picdx >= picdx2)	picdy = picdy2;
	  else	picdx = picdx2;
        }
      else if (ok1) picdy = picdy2;
      else picdx = picdx2;
    }
  xscale = yscale = picdx / wdx;
  if (!ps_aspect) yscale = picdy / wdy;
  dxpic = (int) picdx;
  dypic = (int) picdy;

  xempty = xinch * 72 - ncol * picdx;
  yempty = yinch * 72 - nrow * picdy;
  xsep = ysep = 0;
  if (ncol > 1)
    xsep = (int) ((xempty - 2 * border) / (ncol - 1));
  if (nrow > 1)
    ysep = (int) ((yempty - 2 * border) / (nrow - 1));
  if (ncol == 1)
    xborder = (int)(xempty / (float)2);
  else if (xsep > border)
    xsep = xborder = (int)(xempty / (float)(ncol + 1));
  if (nrow == 1)
    yborder = (int)(yempty / (float)2);
  else if (ysep > border)
    ysep = yborder = (int)(yempty / (float)(nrow + 1));

  row = i / ncol;
  col = i - row * ncol;
  fprintf(fwr, "%c----------- Position picture %d at %d,%d\n", '%', i, row, col);
  if (frd)
    {
      fprintf(fwr, "gsave\n");
      for(i=0; i<8; i++) fgets(text, 100, frd);
    }

  if (portrait)
    {
      tx = xborder + col * (dxpic + xsep);
      ty = yborder + (nrow - 1 - row) * (dypic + ysep);
    }
  else
    {
      tx = yborder + row * (dypic + ysep) + dypic;
      ty = xborder + col * (dxpic + xsep);
    }

  fprintf(fwr, "%d %d translate\n", tx, ty);
  if (!portrait)
    fprintf(fwr, "90 rotate\n");
  fprintf(fwr, "%f %f scale\n", xscale, yscale);
  fprintf(fwr, "%d %d translate\n", -wbord, -wbord);
  fprintf(fwr, "%c------------------------------------------\n\n", '%');
}

/*-----------------------------------------------------------------------------
|	get_maxps, ps_layout, ps_country
-----------------------------------------------------------------------------*/
int get_maxps() { return NLAY; }

void ps_layout(int n, int rows, int cols)
{
  layout[0].npic = n;
  layout[0].orient = 'P';
  layout[0].rows = rows;
  layout[0].cols = cols;
}

void ps_country(int usa)
{
  if (!usa) { paperwidth = (float)8.2; paperheight = (float)11.6; }
}


