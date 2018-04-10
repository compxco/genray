/**************************************************************
**	NAME:	xmotif.c
**************************************************************/
/*
  BUGS:
  Expose redraws window twice sometimes, e.g.
    resize a window bigger
  xdraw.c call to clear_if_needed() shouldn't be needed: 
    find out why won't clear the normal way, e.g. zoom then original
  newevent.xkey.state is hard-coded as 4
  Coord mode on with full-screen window, down to normal window -->
    doesn't redraw.

  IMPLEMENT:
  Click in a window ==> select it
  Get keystrokes working when cursor is in a drawing window
  DISPLAY menu: make different for G, C?
    Nope, if a G can be displayed with c option
  Use a translation table somewhere
  Dialog box for selecting a particular single curve
  Text: make bigger for ctribm1, find those fonts!
  Keep a submenu open after implementation, not as tear-off
  NEXT button: figure out how to add / remove from menu
  CONTOUR LINES menu: give a pull-right
  SINGLE CURVE menu: give pull-right: next, prev, select,
     one psifac / of family, one time / outermost, all
  Coord mode with full screen window: would like to see output!
  Implement keystrokes in windows (not just ctrl-key)
*/

#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <Xm/Xm.h>
#include <X11/Intrinsic.h>
#include <Xm/PushB.h>
#include <Xm/CascadeB.h>
#include <Xm/ToggleB.h>
#include <Xm/RowColumn.h>
#include <Xm/DrawingA.h>
#include <Xm/MessageB.h>
#include <Xm/SelectioB.h>
#include <Xm/XmStrDefs.h>

#include "curves.h"
#include "xtools.h"
#include "xdraw.h"
#include "xdialog.h"

extern void redraw(void);

void fileCall(Widget, XtPointer, XtPointer);
void displayCall(Widget, XtPointer, XtPointer);
void selectCall(Widget, XtPointer, XtPointer);
void exposeCall(Widget, XtPointer, XtPointer);
void biggerCall(Widget, XtPointer, XtPointer);
void zoomCall(Widget, XtPointer, XtPointer);
void keyCall(Widget, XtPointer, XtPointer);
void editCall(Widget, XtPointer, XtPointer);

static void init_motif(int xa, int ya, int nwin);
static void disableZoom(void);
static void motifZoom(Widget, XtPointer, XEvent *);

extern Widget topLevel;
extern XtAppContext app;
extern CURVE_SET curveset[];
extern int redrawflag;

static int selected=0;
static Display *mydisplay;
static Widget menubar, menushell;
static Widget fileSub, selcSub, dispSub, zoomSub, editSub;

typedef struct {
  Widget shell, area, button;
  char text[8];
} GW;

static GW g[10];
static int nwindow;
static Dimension bordw, bordh;

typedef struct {
  Widget button;
  Widget *parent;
  char *id, *label0, *label1;
  char mne;
  void *(func);
} ITEM;

static ITEM mark = { NULL,&dispSub,"Mark","Markers on", "Markers off",
					   'M', displayCall };
static ITEM aspc = { NULL,&dispSub,"Aspc","Preserve Aspect Ratio","Auto Scale",
					   'A', displayCall };
static ITEM alll = { NULL,&dispSub,"Alll","Single Curve", "All Curves",
					   'C', displayCall };
static ITEM cont = { NULL,&dispSub,"Cont","Contour Lines", "",
		       			   'L', displayCall };
static ITEM next = { NULL,&dispSub,"Next","Next Curve", "",
					   'N', displayCall };

static ITEM zoom = { NULL,&zoomSub,"Zoom","Begin Zoom",
		                          "Abort Zoom",       'Z', zoomCall };
static ITEM valu = { NULL,&zoomSub,"Valu","Get Coordinates...",
					  "Done Coordinates...", 'C', zoomCall };
static ITEM slop = { NULL,&zoomSub,"Slop","Get Slope",
		                          "Abort Slope",      'S', zoomCall };
static ITEM rato = { NULL,&zoomSub,"Rato","Get Ratio",
		       	                  "Abort Ratio",      'R', zoomCall };
				  
static ITEM orig = { NULL,&zoomSub,"Orig","Original size", "", 'O', biggerCall };
static ITEM bigg = { NULL,&zoomSub,"Bigg","Bigger",        "", 'B', biggerCall };
static ITEM smal = { NULL,&zoomSub,"Smal","Smaller",       "", 'm', biggerCall };

static ITEM prnt = { NULL,&fileSub,"Prnt","Print",   "", 'P', fileCall };
static ITEM quit = { NULL,&fileSub,"Quit","Quit...",    "", 'Q', fileCall };

static ITEM edit = { NULL,&editSub,"Edit","Configure","",'C', editCall };

static ITEM *buttonList[] = {
  &mark, &aspc, &cont, &alll, &next,
  &zoom, &orig, &valu, &slop, &rato, &bigg, &smal,
  &prnt, &quit, &edit
  } ;

#define NITEM (sizeof(buttonList) / sizeof(ITEM *) )

Widget selectButton[9];

#define DefineCall(b,f,i) XtAddCallback(b,XmNactivateCallback,f,(XtPointer)i);
#define RadioCall(b,f,i) XtAddCallback(b,XmNvalueChangedCallback,f,(XtPointer)i);
#define ExposeCall(b,f,i) XtAddCallback(b,XmNexposeCallback,f,(XtPointer)i);
#define KeyCall(b,f,i) XtAddCallback(b,XmNinputCallback, f, (XtPointer)i);

static Widget zoom_caller = NULL, zoom_area = NULL;
static int zoom_count = 0;
static int xcurs, ycurs, xfirst, yfirst;
/*#define ZOOM_MASK (ButtonPressMask | ButtonReleaseMask | \
                   PointerMotionMask | LeaveNotifyMask)*/
#define ZOOM_MASK (ButtonReleaseMask | PointerMotionMask)

/*-----------------------------------------------------------------------------
|	trace
-----------------------------------------------------------------------------*/
/*#define TRACE*/
#ifdef TRACE
int nev=0, ev[5000];
static void trace(int type, int show)
{
  int etype, ok, i;

  if (event->type == MotionNotify) etype=1;
  else if (event->type == ButtonPress) etype=2;
  else if (event->type == ButtonRelease) etype=3;
  else etype=4;
  ok = 1;
  if (etype==1 && nev && ev[nev-1]==0x11) ok=0;
  else ev[nev++] = etype;
}
#endif
/*=============================================================================
**		Functions Called From xtools.c
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	enableWidget
|	* called from makewindow() BEFORE buttons "Select" are created
-----------------------------------------------------------------------------*/
void enableWidget(int i, Widget shell, Widget area, char *text)
{
  g[i].shell = shell;
  g[i].area = area;
  strcpy(g[i].text, text);
  XtManageChild(area);
  XtRealizeWidget(shell);
  XtVaSetValues(area, XtVaTypedArg, XmNbackground,
                      XmRString, "black", strlen("black"), NULL);
}

/*-----------------------------------------------------------------------------
|	execute_motif -- from init_menus()
-----------------------------------------------------------------------------*/
void execute_motif(int xa, int ya, int nwin, Display *display)
{
  
  mydisplay = display;
  test_window("execute_motif", XtWindow(g[0].shell), g[0].shell);
  init_motif(xa, ya, nwin);
  XtAppMainLoop(app);
}

void get_bord()
{
  XtVaGetValues(g[0].shell, XmNx, &bordw, XmNy, &bordh, NULL);
  printf("bord = %d %d\n", bordw, bordh);
}

/*=============================================================================
**		init_motif() and its Service Functions
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	CreateSubmenu
-----------------------------------------------------------------------------*/
Widget CreateSubmenu(char *submenuId, char *buttonId, char *label, char c)
{
  Widget submenu;
  XmString text;

  submenu = XmCreatePulldownMenu(menubar, submenuId, NULL, 0);
  if (!strcmp(submenuId, "selcSub"))
    XtVaSetValues(submenu,
		  XmNradioBehavior, True, XmNradioAlwaysOne, True, NULL);
  else if (!strcmp(submenuId, "dispSub") ||
	   !strcmp(submenuId, "zoomSub"))
      XtVaSetValues(submenu, XmNtearOffModel, XmTEAR_OFF_ENABLED, NULL);

  text = XmStringCreateLocalized(label);
  XtVaCreateManagedWidget(buttonId, xmCascadeButtonWidgetClass,
			  menubar, XmNlabelString, text,
			  XmNmnemonic, c, XmNsubMenuId, submenu,
			  NULL);
  XmStringFree(text);
  return submenu;
}

/*-----------------------------------------------------------------------------
|	CreateButton
-----------------------------------------------------------------------------*/
Widget CreateButton(Widget menu, char *ButtonId, char *label, char c)
{
  Widget button;
  XmString text;
  int i;

  text = XmStringCreateLocalized(label);

  if (menu != selcSub)
    {
      button = XtVaCreateManagedWidget(ButtonId, xmPushButtonWidgetClass,
	      menu, XmNlabelString, text, XmNmnemonic, c, NULL);
    }
  else
    {
      button = XtVaCreateManagedWidget(ButtonId, xmToggleButtonWidgetClass,
	       menu, XmNlabelString, text, XmNmnemonic, c,
	       XmNradioBehavior, True, XmNset, (c=='0' ? True : False),
	       NULL);
      i = (int)('c'-'0');
      if (i >= nwindow) i = nwindow - 1;
      selectButton[i] = button;
    }

  XmStringFree(text);
  return button;
}

/*-----------------------------------------------------------------------------
|	TitleMenu -- title of menu window changes when select graphics window
-----------------------------------------------------------------------------*/
static void TitleMenu(void)
{
    char text[80];
    sprintf(text, "xdraw, current window = %d", selected);
    XtVaSetValues(menushell, XmNtitle, text, NULL);
}

/*-----------------------------------------------------------------------------
|	AddAccelerator
-----------------------------------------------------------------------------*/
static void AddAccelerator(Widget button, char mne)
{
  char text[80], key[80];
  XmString accel_text;

  sprintf(key, "Ctrl<Key>%c", mne);
  sprintf(text, "Ctrl+%c", mne);
  accel_text = XmStringCreateLocalized(text);
  XtVaSetValues(button, XmNaccelerator, key,
		XmNacceleratorText, accel_text, NULL);
  XmStringFree(accel_text);
}

/*-----------------------------------------------------------------------------
|	init_motif
-----------------------------------------------------------------------------*/
static void init_motif(int xa, int ya, int nwin)
{
  Widget shell, area;
  ITEM *p;
  Position x,y;
  char text[30];
  int i;

  nwindow = nwin;
  x = xa;
  y = ya;
  shell = XtVaAppCreateShell(NULL, "XDraw", topLevelShellWidgetClass,
	     mydisplay, XmNx, x, XmNy, y, NULL);

/*----- Main menu and first submenu -----*/
  menubar = XmVaCreateSimpleMenuBar(shell, "menubar", NULL);
  menushell = shell;
  TitleMenu();

  fileSub = CreateSubmenu("fileSub", "File", "File", 'F');
  selcSub = CreateSubmenu("selcSub", "Select", "Select", 'S');
  dispSub = CreateSubmenu("dispSub", "Display", "Display", 'D');
  zoomSub = CreateSubmenu("zoomSub", "Zoom", "Zoom", 'Z');
  editSub = CreateSubmenu("editSub", "Edit", "Edit", 'E');

  for(i=0; i<NITEM; i++)
    {
      p = buttonList[i];
      p->button = CreateButton(*(p->parent), p->id, p->label0, p->mne);
      DefineCall(p->button, p->func, i);
    }
  AddAccelerator(quit.button, quit.mne);
  AddAccelerator(valu.button, valu.mne);

  for(i=0; i<nwin; i++)		/* Make labels for selection menu */
    {
      sprintf(text, "Window %d", i);
      g[i].button = CreateButton(selcSub, g[i].text, text, *(text+7));
      RadioCall(g[i].button, selectCall, i);
      ExposeCall(g[i].area,  exposeCall, i);
      KeyCall(g[i].area, keyCall, i);
      AddAccelerator(g[i].button, *(text+7));
    }

  XtManageChild(menubar);
  XtRealizeWidget(shell);
  XtUnmapWidget(next.button);
}

/*=============================================================================
**			Callback Functions
=============================================================================*/
/*-----------------------------------------------------------------------------
|	newLabel
-----------------------------------------------------------------------------*/
static void newLabel(Widget w, int flag, char *zero, char *one)
{
  char *text;
  if (flag) text = XmStringCreateLocalized(one);
  else text = XmStringCreateLocalized(zero);
  XtVaSetValues(w, XmNlabelString, text, NULL);
  XmStringFree(text);
}

#define get_cp() (curveset + get_graph(0, selected))
/*------------------------------------------------------------
|	displayCall
------------------------------------------------------------*/
void displayCall(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int ident, flag, singleToAll;
  CURVE_SET *cp;
  
  disableZoom();
  cp = get_cp();

  if (widget == mark.button)		 	/*----- Markers */
    {
      toggle_markers(cp, &flag);
      newLabel(mark.button, flag, mark.label0, mark.label1);
    }
  else if (widget == aspc.button)		/*----- Aspect */
    {
      toggle_aspect(cp, &flag);
      newLabel(aspc.button, flag, aspc.label0, aspc.label1);
    }
  else if (widget == alll.button)
    {
      flag = (cp->flags & SINGLE);
      toggle_single(cp, flag ?  '2' : '1');
      newLabel(alll.button, flag ? 0 : 1, alll.label0, alll.label1);
      if (flag) XtUnmapWidget(next.button);
      else
	  XtMapWidget(next.button);
    }
  else if (widget == next.button)
    toggle_single(cp, '1');

  if (redrawflag)
    {
      set_expose(selected);
      redraw();
    }
}

/*-----------------------------------------------------------------------------
|	setAllLabels
-----------------------------------------------------------------------------*/
void setAllLabels()
{
  CURVE_SET *cp;
  int flag;
  cp = get_cp();
  flag = (cp->flags & MARKERS) ? 1 : 0;
  newLabel(mark.button, flag, mark.label0, mark.label1);
  flag = (cp->flags & ASPECT) ? 1 : 0;
  newLabel(mark.button, flag, mark.label0, mark.label1);
  TitleMenu();
}

/*------------------------------------------------------------
|	selectCall
|	* behavior for radio and select new one ==> 
|	  calls twice: unset the old, set the new
------------------------------------------------------------*/
void selectCall(Widget w, XtPointer client_data, XtPointer call_data)
{
  int old_selected, set;
  char text[30];
  Widget button;

  /*printf("In selectCall()\n");*/
  disableZoom();
  old_selected = selected;
  selected = (int)client_data;
  if (call_data)
    set = ((XmToggleButtonCallbackStruct *)call_data)->set;

  if (selected != old_selected)
  {
     if (call_data == NULL)
       {
	 XtVaSetValues(selectButton[selected],     XmNset, True, NULL);
	 /*XtVaSetValues(selectButton[old_selected], XmNset, False, NULL);*/
       }
     setAllLabels();
     XtVaSetValues(selcSub, XmNinitialFocus, g[selected].button);
  }
}

/*-----------------------------------------------------------------------------
|	keyCall
|	* inputCallback: events in graphics windows send you here
|	* possible events: INPUT (keys, mouse), EXPOSE, RESIZE
-----------------------------------------------------------------------------*/
void keyCall(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int which, ic;
  char c, text[80], *q;
  int k,s;
  KeySym ks;
  KeyCode kc;
  XmDrawingAreaCallbackStruct *p;
  XEvent *event, newevent;

  p = (XmDrawingAreaCallbackStruct *)call_data;
  event = p->event;
  which = (int)client_data;

  if (p->reason == XmCR_INPUT)
    {
      if (event->xany.type == KeyPress)
	{
	  tell_key(event); ic = parsekey(); c = *(char *)&ic;
	  k = event->xkey.keycode;
	  s = event->xkey.state;
	  ks = XLookupKeysym((XKeyEvent *)event, 0);
	  q = XKeysymToString(ks);
	  kc = XKeysymToKeycode(mydisplay, ks);
	  /*printf("Keycode = %lx (%s) --> %lx, state=%lx\n",k,q,kc,s);*/
	}
      if (event->xany.type == KeyPress)
	{
	  event->xany.window = XtWindow(menushell);
	  XPutBackEvent(mydisplay, event);
	}

      else if (event->xany.type == ButtonPress &&
	       (!zoom_caller || event->xany.window != XtWindow(zoom_area)))
	{
	  newevent = *event;
	  newevent.xkey.window = XtWindow(menushell);
	  newevent.xkey.type = KeyPress;
	  *text = which + '0';
	  *(text+1) = '\0';
	  ks = XStringToKeysym(text);
	  kc = XKeysymToKeycode(mydisplay, ks);
	  newevent.xkey.keycode =  kc;
	  newevent.xkey.state = 4;	/* WARN */
	  XSendEvent(mydisplay, XtWindow(menushell), True,
		     KeyPressMask, &newevent);
	}
    }
}

/*-----------------------------------------------------------------------------
|	biggerCall
-----------------------------------------------------------------------------*/
void biggerCall(Widget widget, XtPointer client_data, XtPointer call_data)
{
  disableZoom();
  if (widget == orig.button)
  {
    tell_zoom(0, selected, 0);
    newclip(-1,0,0);
  }
  else addzoom(XtWindow(g[selected].area), widget == bigg.button);
  set_expose(selected);
  redraw();
}

/*------------------------------------------------------------
|	fileCall
------------------------------------------------------------*/
void fileCall(Widget widget, XtPointer client_data, XtPointer call_data)
{
  CURVE_SET *cp;
  char text[100];

  disableZoom();
  cp = get_cp();
  if (widget == prnt.button)
    {
      parse_title(cp, cp->title, text);
      postscript(text);
    }
  else if (widget == quit.button)
    { printf("It's been nice knowing you!\n"); exit(0); }
}

/*-----------------------------------------------------------------------------
|	editcall
-----------------------------------------------------------------------------*/
void editCall(Widget widget, XtPointer client_data, XtPointer call_data)
{
  Widget edit_dialog, shell;
  int x, y;
  unsigned int wdx, wdy;
  Dimension ew, eh, dx, dy;
  Position x0, y0;

  shell = g[selected].shell;
  edit_dialog = CreateEditDialog(widget, selected, mydisplay);
  get_screensize(&wdx, &wdy);

  XtVaGetValues(edit_dialog, XmNwidth, &ew, XmNheight, &eh, NULL);
  XtVaGetValues(shell, XmNx, &x0, XmNy, &y0,
		XmNwidth, &dx, XmNheight, &dy, NULL);
  printf("shell %lx, x,y = %d %d, size = %d %d\n", shell, x0, y0, dx, dy);
  printf("border dx,dy = %d %d\n", bordw, bordh);
  printf("dialog dx,dy = %d %d\n", ew, eh);
  x = x0 + bordw;
  if (x0 + ew >= wdx)
    x = wdx - ew - bordw;
  y = y0 + dy + bordw + bordh;
  if (y + eh + bordw >= wdy) y = wdy - eh - bordw;
  XtMoveWidget(edit_dialog, x, y);
  XtPopup(edit_dialog, XtGrabExclusive);
}

/*-----------------------------------------------------------------------------
|	exposeCall
-----------------------------------------------------------------------------*/
void exposeCall(Widget w, XtPointer client_data, XtPointer call_data)
{
  int ident;
  ident = (int)client_data;
  set_expose(ident);
  redraw();
}

/*=============================================================================
**			Zoom, Value etc. Functions
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	disableZoom
-----------------------------------------------------------------------------*/
#define zenable(fl,z) flag=fl; newLabel(z.button,1,z.label0,z.label1)
#define zdisable(z) if (zoom_caller==z.button) \
                    newLabel(z.button,0,z.label0,z.label1)

static void disableZoom()
{
   if (zoom_caller != NULL)
   {
      zdisable(zoom);		/* newLabel() as appropriate */
      zdisable(valu);
      zdisable(slop);
      zdisable(rato);
      zoom_caller = NULL;
      crosshair(xcurs, ycurs);
      if (zoom_count >= 1) crosshair(xfirst, yfirst);
      XtRemoveEventHandler(zoom_area, ZOOM_MASK, False, 
			   (XtEventHandler)motifZoom, NULL);
   }
}

/*-----------------------------------------------------------------------------
|	zoomCall -- when click on Zoom, Coordinate, etc.
-----------------------------------------------------------------------------*/
void zoomCall(Widget widget, XtPointer client_data, XtPointer call_data)
{
  int xr, yr;
  Window rw, cw, zwin;
  unsigned int keys_buttons;
  int flag;
  int x0, y0, xw, yw;
  unsigned int dx, dy, bwidth, depth;

  if (zoom_caller != NULL)		/* If already on, disable mode */
    disableZoom();
  else
    {
      zoom_caller = widget;		/* Save the button that called it */
      zoom_count = 0;
      zoom_area = g[selected].area;
      if (widget == zoom.button) { zenable(0, zoom); }
      else if (widget == slop.button) { zenable(2, slop); }
      else if (widget == rato.button) { zenable(3, rato); }
      else { zenable(1, valu); }

      zwin = XtWindow(zoom_area);
      XGetGeometry(mydisplay, zwin, &rw, &x0, &y0, &dx, &dy, &bwidth, &depth);
      xw = dx / 2;
      yw = dy / 2;
      XWarpPointer(mydisplay, None, zwin, 0,0,0,0, xw,yw);
      XQueryPointer(mydisplay, zwin, &rw, &cw, &xr, &yr,
		    &xcurs, &ycurs, &keys_buttons);
      tell_zoom(zwin, selected, flag);
      crosshair(xcurs, ycurs);
      XtInsertEventHandler(zoom_area, ZOOM_MASK, False,
			   (XtEventHandler)motifZoom, NULL, XtListHead);
    }
}

/*-----------------------------------------------------------------------------
|	motifZoom -- handles motion & button events for zoom etc.
-----------------------------------------------------------------------------*/
static void motifZoom(Widget widget, XtPointer client_data, XEvent *event)
{
  Window rw, cw;
  int xr,yr, xw,yw;
  unsigned int keys_buttons;
  static int armed=0;

   XQueryPointer(mydisplay, event->xany.window,
		&rw, &cw, &xr, &yr, &xw, &yw, &keys_buttons);

   if (event->type == MotionNotify)
   {
      crosshair(xcurs, ycurs);
      crosshair(xw, yw);
      xcurs = xw; ycurs = yw;	/* This WAS outside this if */
   }

   if (event->type == ButtonRelease)
   {
      zoom_count++;
      if (zoom_count == 1)
      {
	 xfirst = xcurs; yfirst = ycurs;
	 if (zoom_caller != valu.button) crosshair(xcurs, ycurs);
      }
      if (zoom_caller == zoom.button)		/* ......Zoom */
      {
	 newclip((zoom_count==1) ? 0 : 1, xcurs, ycurs);
         if (zoom_count == 2)
	 {
	    XClearArea(mydisplay, XtWindow(zoom_area), 0,0,0,0,True);
	    disableZoom();
	 }
      }
      else					/* ......Valu,Slop,etc */
      {
	 if (zoom_count == 2) disableZoom();
	 show_coord(xcurs, ycurs, zoom_count);
	 fflush(stdout);
	 if (zoom_caller == valu.button || zoom_count == 2)
	    zoom_count = 0;
      }
   }
}


