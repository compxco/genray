/**************************************************************
**	NAME:	xdialog.c
**************************************************************/
/*
   Add text entry: which family member if single
   "OK" ==> redraw edited curve
   NOTE for lstep: (a) if x & y = variables, lstep = prev node?
      not necessarily, if all of node loop=same (b) else lstep = x
   ctribm1 dotted lines in lists are wierd

   Position of menu box: x = +1 bordw; y never is above; x on right
      is at right of screen, not of current shell
   You can select from main menu while dialog menu is up!!
*/
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/PanedW.h>
#include <Xm/PushBG.h>
#include <Xm/ToggleBG.h>
#include <Xm/LabelG.h>
#include <Xm/List.h>
#include <Xm/RowColumn.h>
#include <Xm/TextF.h>
#include <Xm/SelectioB.h>
#include <Xm/MessageB.h>
#include <Xm/Separator.h>
#include <Xm/Frame.h>
#include <Xm/XmStrDefs.h>

#include "xdialog.h"
#include "curves.h"
#include "xedit.h"

typedef unsigned char byte;
static Widget AddLabel(Widget, String);
static Widget AddText(Widget, String);
static Widget AddText(Widget, String);
static Widget AddList(Widget, int);
static Widget AddRadio(Widget, int, String *, int, int);
static Widget AddSeparator(Widget form2, int vert, int n);
static Widget CreateSelectionBox(Widget, int, char *);
static void ActionButton(Widget, String, int, int, Boolean, int);
static void PromptBox(Widget, char *, int);
static void ErrorBox(Widget, char *, char *);
static Widget GetWiderWidget(Widget label1, Widget list1);
static void AddVRLists(Widget form2, Widget *x, Widget *y);
static void AddVarMenu(Widget form2, Widget *x, Widget *y);

static void buttonCall(Widget, XtPointer, XtPointer);
static void listCall(Widget, XtPointer, XtPointer);
static void roleCall(Widget, XtPointer, XtPointer);
static void xyCall(Widget, XtPointer, XtPointer);
static void otherCall(Widget, XtPointer, XtPointer);
static void toggleCall(Widget, XtPointer, XtPointer);

extern int ncset, nvar, nloop;
extern CURVE_SET curveset[];
extern VARIABLE_ID varid[];
extern LOOP loop[];
extern int redrawflag;

static Widget dialog, descriptive, desctext, rolelist;
static Widget varlist=NULL, indlist=NULL;
static selected_role = -1;
static int nui = 0;		/* nb undefined indexed */
static int selected_win;
static CURVE_SET new_curve;
extern int param[];

#define Awid  XmATTACH_WIDGET
#define AOwid XmATTACH_OPPOSITE_WIDGET
#define Aform XmATTACH_FORM

#define Atop XmNtopAttachment
#define Abot XmNbottomAttachment
#define Alft XmNleftAttachment
#define Arit XmNrightAttachment
#define Atopw Atop,Awid,XmNtopWidget
#define Abotw Abot,Awid,XmNbottomWidget
#define Alftw Alft,Awid,XmNleftWidget
#define AlftwO Alft,AOwid,XmNleftWidget

static char *role_list[] = {
  "X-axis", "Y-axis", "Step", "Family", "Parameter"
  };

static char *single[] = { "Single curve", "All" };
static char *contour[] = { "Y-vs-X", "Contour" };
#define NSING (sizeof(single)/sizeof(char *))
#define NCONT (sizeof(contour)/sizeof(char *))

#define MSHIFT 0x100
#define SING 0
#define CONT 1

/*=============================================================================
**		Functions Begin
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	GetTopShell
-----------------------------------------------------------------------------*/
Widget GetTopShell(Widget w)
{
  while (w && !XtIsWMShell(w)) w = XtParent(w);
  return w;
}

/*-----------------------------------------------------------------------------
|	CreateEditDialog
-----------------------------------------------------------------------------*/
Widget CreateEditDialog(Widget parent, int selected, Display *disp)
{
  Widget pane, form1, form2, form3, action, widget;
  Widget label1, label2, label3, label4, xsep, ysep;
  Widget radio1, radio2;
  char text[280], *p;
  int i;

  widget = GetTopShell(parent);
  dialog = XtVaCreatePopupShell("EditShell", xmDialogShellWidgetClass,
		widget, XmNdeleteResponse, XmDESTROY,
		NULL);
  pane = XtVaCreateWidget("pane", xmPanedWindowWidgetClass, dialog,
		XmNsashWidth,1, XmNsashHeight,1, NULL);

/*--------- Titles */
  form1   = XtVaCreateWidget("titl", xmFormWidgetClass, pane, NULL);

  selected_win = selected;
  i = get_graph(0, selected);
  sprintf(text, "Window %d: %s", selected, curveset[i].title);
  label1 = AddLabel(form1, text);

  new_curve = *(curveset+i);
  for(i=0; i<nloop; i++) param[i] = new_curve.param[i];

/*-------- Work Region */
  form2   = XtVaCreateWidget("work", xmFormWidgetClass, pane, NULL);

  sprintf(text, "\"%s\" changes to...",
	  p = get_descriptive_title(selected, &new_curve));
  label2 = AddLabel(form2, text);
  desctext = AddText(form2, p);
  descriptive = XtParent(desctext);
  XtVaSetValues(label2, Atop, Aform, Alft, Aform, NULL);
  XtVaSetValues(descriptive, Atopw, label2, Alft, Aform, NULL);

  /*AddVRLists(form2, &xsep, &ysep);*/
  AddVarMenu(form2, &xsep, &ysep);

  /*label3 = AddLabel(form2, "Show Family:");
    radio1 = AddRadio(form2, SING, single, NSING,
		    (new_curve.flags & SINGLE)?0:1);
  */
  label4 = AddLabel(form2, "Graph Type:");
  radio2 = AddRadio(form2, CONT, contour,NCONT,
		    (new_curve.gtype=='G')?0:1);
  XtVaSetValues(label4, Atopw, ysep, Alftw, xsep, NULL);
  XtVaSetValues(radio2, Atopw, label4, AlftwO, label4, NULL);

  form3   = XtVaCreateWidget("ctrl", xmFormWidgetClass, pane, NULL);

  action = XtVaCreateWidget("actn", xmFormWidgetClass, pane,
		XmNfractionBase, 5, NULL);
  ActionButton(action, "OK",     1, 2, True,  0);
  ActionButton(action, "Cancel", 3, 4, False, 1);

  XtManageChild(form1);
  XtManageChild(form2);
  XtManageChild(action);
  XtManageChild(pane);
  return dialog;
}

/*-----------------------------------------------------------------------------
|	AddVRLists: lists at left of work area = Variables, Roles
-----------------------------------------------------------------------------*/
void AddVRLists(Widget form2, Widget *x, Widget *y)
{
  Widget label1, label2, list1, list2, xsep1, xsep2, ysep, wider;
  label1 = AddLabel(form2, "Variables:");
  label2 = AddLabel(form2, "Role:");
  list1 = AddList(form2, 0);
  list2 = AddList(form2, 1);
  xsep1 = AddSeparator(form2, XmVERTICAL, 20);
  xsep2 = AddSeparator(form2, XmVERTICAL, 20);
  ysep  = AddSeparator(form2, XmHORIZONTAL, 20);
  XtVaSetValues(ysep, Atopw, descriptive, Alft, Aform, NULL);

  wider = GetWiderWidget(label1, list1);
  XtVaSetValues(label1, Atopw, ysep, Alft, Aform, NULL);
  XtVaSetValues(xsep1, Atopw, ysep, Alftw, wider, NULL);
  XtVaSetValues(label2, Atopw, ysep, Alftw, xsep1, NULL);
  wider = GetWiderWidget(label2, list2);
  XtVaSetValues(xsep2, Atopw, ysep, Alftw, wider, NULL);

  XtVaSetValues(XtParent(list1), Atopw, label1, Alft, Aform, NULL);
  XtVaSetValues(XtParent(list2), Atopw, label2,	AlftwO, label2, NULL);
  *x = xsep2;
  *y = ysep;
}

/*-----------------------------------------------------------------------------
|	AddVarMenu: list at left of work area = "X-variable = %s" etc
-----------------------------------------------------------------------------*/
void AddVarMenu(Widget form2, Widget *x, Widget *y)
{
  Widget label1, label2, xsep, ysep, wider;
  label1 = AddLabel(form2, "Select item to change");
  label2 = AddLabel(form2, "(starred items mean inconsistent)");
  rolelist = AddList(form2, 2);
  xsep = AddSeparator(form2, XmVERTICAL, 20);
  ysep  = AddSeparator(form2, XmHORIZONTAL, 20);
  XtVaSetValues(ysep, Atopw, descriptive, Alft, Aform, NULL);

  wider = GetWiderWidget(label1, label2);
  XtVaSetValues(label1, Atopw, ysep, Alft, Aform, NULL);
  XtVaSetValues(label2, Atopw, label1, Alft, Aform, NULL);
  XtVaSetValues(XtParent(rolelist), Atopw, label2, Alft, Aform, NULL);
  wider = GetWiderWidget(wider, rolelist);
  XtVaSetValues(xsep, Atopw, ysep, Alftw, wider, NULL);

  *x = xsep;
  *y = ysep;
}

/*-----------------------------------------------------------------------------
|	GetWiderWidget
-----------------------------------------------------------------------------*/
Widget GetWiderWidget(Widget label1, Widget list1)
{
  Dimension wl1, wl2;
  Widget wider;

  XtVaGetValues(label1, XmNwidth, &wl1, NULL);
  XtVaGetValues(list1,  XmNwidth, &wl2, NULL);
  wider = (wl1 >= wl2) ? label1 : list1;
  return wider;
}

/*=============================================================================
**			Basic Widgets for Dialog
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	AddSeparator
-----------------------------------------------------------------------------*/
Widget AddSeparator(Widget form2, int vert, int n)
{
  Widget sep;
  sep = XmCreateSeparator(form2, "sep", NULL, 0);
  XtVaSetValues(sep, XmNorientation, vert, XmNseparatorType,
		XmNO_LINE, (vert==XmVERTICAL)?XmNwidth:XmNheight, n, NULL);
  XtManageChild(sep);
  return sep;
}

/*-----------------------------------------------------------------------------
|	getlbl -- called by get_variables_menu
-----------------------------------------------------------------------------*/
void getlbl(int index, char *s)
{
   int i;
   char text[80];
   *s = '\0';
   for(i=0; i<nvar && varid[i].index != index; i++) ;
   if (i<nvar)
   {
      strcpy(s, varid[i].name);
      if (index & LOOPBIT) strcat(s,"=");
   }
   if (index & LOOPBIT) sprintf(s+strlen(s),"<i%d>", index & INDBITS);
}

/*-----------------------------------------------------------------------------
|	get_variables_menu
-----------------------------------------------------------------------------*/
char *get_variables_menu(int *np)
{
  char *inames, *p, s[100];
  int i,n,nother,count,*q;
  static char *lbl[] = {
    " X-variable = %s", " Y-variable = %s", " Index/step-along-curves = %s",
    " Index/family-of-curves = %s", " Index/specified-value: %s"
    } ;

  nother = nloop - 3;
  inames = (char *)malloc( (4 + nother) * 2 * strlen(lbl[4]));
  p = inames;
#define addlbl(i,j) getlbl(j,s); sprintf(p,lbl[i],s); p+=strlen(p)+1
#define addlbl2(i,j,k) addlbl(i,j); sprintf(p,"=%d",k); p+=strlen(p)+1
  addlbl(0, new_curve.ix.index);
  addlbl(1, new_curve.iy.index);
  addlbl(2, new_curve.lstep | LOOPBIT);
  addlbl(3, new_curve.lfaml | LOOPBIT);
  for(i=n=0,q=new_curve.param; i<nloop; i++,q++)
  {
    if (*q >= 0) { addlbl2(4, i | LOOPBIT, *q); n++; }
  }

  if (n != nother) printf("Error: inconsistent nb 'other'\n");
  *np = n + 4;
  return inames;
}

/*-----------------------------------------------------------------------------
|	AddList
|	* which: 0==>all variables, 1==> all roles, 2==>role+variable+compatible
-----------------------------------------------------------------------------*/
Widget AddList(Widget form, int which)
{
  XmStringTable str_list;
  Widget list;
  int i, j, n, m;
  char *p, *inames;
  static XtCallbackProc calls[] = { listCall, roleCall, xyCall, otherCall } ;
  static int initialized[] = { 0, 0, 0, 0 };

  if (which==0)
    {
      n = nvar;
      m = strlen("<i00>") + 1;
      inames = (char *)malloc(nloop * m);
      for(i=nui=0; i<nloop; i++)
	{
	  if (loop[i].ltype == 'V') continue;
	  for(j=0; j<nvar; j++)
	    if (varid[j].index == (i | LOOPBIT)) break;
	  if (j==nvar)
	    sprintf(inames + m * (nui++), "<i%d>", i);
	}
      n += nui;
    }
  else if (which==1)
    n = sizeof(role_list)/sizeof(char *);
  else inames = get_variables_menu(&n);
  str_list = (XmStringTable) XtMalloc(n * sizeof(XmString *));

  for(i=0; i<n; i++)
    {
      if (which==1) p = role_list[i];
      else if (which==2) p = (i==0) ? inames : p + strlen(p) +1;
      else if (which==0 && i < nvar) p = varid[i].name;
      else p = inames + (i-nvar)*m;
      str_list[i] = XmStringCreateLocalized(p);
    }
  if (which==0 || which==2) free(inames);

  list = XmCreateScrolledList(form, "List", NULL, 0);
  XtVaSetValues(list, XmNvisibleItemCount, (n <= 5) ? n : 5,
		XmNitemCount, n, XmNitems, str_list,
		XmNselectionPolicy, XmSINGLE_SELECT,
		/*XmNborderWidth, 10,*/
		NULL);
  for(i=0; i<n; i++) XmStringFree(str_list[i]);

  XtAddCallback(list, XmNsingleSelectionCallback, calls[which], NULL);
  XtManageChild(list);
  return list;
}

/*-----------------------------------------------------------------------------
|	AddLabel
-----------------------------------------------------------------------------*/
Widget AddLabel(Widget form, String text)
{
  XmString str;
  Widget label;

  str = XmStringCreateLocalized(text);
  label = XtVaCreateManagedWidget("label", xmLabelGadgetClass, form,
		XmNalignment, XmALIGNMENT_BEGINNING,
		XmNlabelString, str, NULL);
  XmStringFree(str);
  return label;
}

/*-----------------------------------------------------------------------------
|	AddText
-----------------------------------------------------------------------------*/
Widget AddText(Widget form, String text)
{
  Widget label, frame;
  XmString str;

  frame = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, form,
				  XmNshadowType, XmSHADOW_ETCHED_IN, NULL);
  str = XmStringCreateLocalized(text);
  label = XtVaCreateManagedWidget("label", xmLabelGadgetClass, frame,
		XmNalignment, XmALIGNMENT_BEGINNING,
		XmNlabelString, str, NULL);
  XmStringFree(str);
#ifdef DEAD_CODE
  label = XtVaCreateManagedWidget("text", xmTextFieldWidgetClass, form,
		XmNeditable, False, XmNcursorPositionVisible, False,
		XmNvalue, text, XmNcolumns, strlen(text)+10, NULL);
#endif
  return label;
}

/*-----------------------------------------------------------------------------
|	AddRadio
-----------------------------------------------------------------------------*/
static Widget AddRadio(Widget form, int which, String *s, int ns,
		       int whichset)
{
  Widget radio, button, frame;
  int i, j;
  XmString str;
  Boolean set;

  frame = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, form,
				  XmNshadowType, XmSHADOW_ETCHED_IN, NULL);
  radio = XmCreateRadioBox(frame, "radio", NULL, 0);
  XtVaSetValues(radio, XmNmarginHeight, 0, NULL);
  /*XtVaSetValues(frame, XmNnumChildren, 1, XmNchildren, &radio);*/
  j = MSHIFT * which;
  for(i=0; i<ns; i++)
    {
      button = XtVaCreateManagedWidget("But", xmToggleButtonGadgetClass,
				       radio, NULL);
      XtAddCallback(button, XmNvalueChangedCallback, toggleCall,
		    (XtPointer)(j | i));
      str = XmStringCreateLocalized(s[i]);
      set = (i==whichset) ? True : False;
      XtVaSetValues(button, XmNlabelString, str, XmNset, set,
		    XmNfillOnSelect, True, XmNselectColor,0, NULL);
      XmStringFree(str);
    }
  XtManageChild(radio);
  /*XtManageChild(frame);*/
  /*return radio;*/
  return frame;
}
/*-----------------------------------------------------------------------------
|	ActionButton
-----------------------------------------------------------------------------*/
void ActionButton(Widget action, String s, int i1, int i2, 
		  Boolean state, int which)
{
  Widget button;
  XmString str;
  str = XmStringCreateLocalized(s);
  button = XtVaCreateManagedWidget(s, xmPushButtonGadgetClass, action,
		XmNtopAttachment,   XmATTACH_FORM,
		XmNbottomAttachment,XmATTACH_FORM,
		XmNleftAttachment,  XmATTACH_POSITION,
		XmNrightAttachment, XmATTACH_POSITION,
		XmNleftPosition, i1, XmNrightPosition, i2,
		XmNshowAsDefault, state,
		XmNdefaultButtonShadowThickness, 1,
		XmNlabelString, str,
		NULL);
  XmStringFree(str);
  XtAddCallback(button, XmNactivateCallback, buttonCall, (XtPointer)which);
}

void ShowNewTitle()
{
  char text[100];
  sprintf(text, "%s",
	  get_descriptive_title(selected_win, &new_curve));
  XtVaSetValues(desctext, XmNvalue, text, NULL);
}

/*=============================================================================
**			Callbacks
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	listCall -- selected = one of the variables
-----------------------------------------------------------------------------*/
void listCall(Widget w, XtPointer client_data, XtPointer call_data)
{
  XmListCallbackStruct *p;
  char *choice;
  int which, ind;

  p = (XmListCallbackStruct *)call_data;
  if (p->selected_item_count == 0) return;
  XmStringGetLtoR(p->item, XmFONTLIST_DEFAULT_TAG, &choice);
  which = p->item_position - 1;
  printf("Count %d\n", p->selected_item_count);
  if (selected_role==-1) printf("Role not defined for %s\n", choice);
  else
    {
      printf("Selected item %d: %s\n", which, choice);
      if (which >= nvar)
	{
	  *(strchr(choice, '>')) = '\0';
	  sscanf(strchr(choice, '<'), "%d", &ind);
	  ind |= LOOPBIT;
	}
      else ind = varid[which].index;

      if (!(ind & LOOPBIT) && selected_role > 1)
	ErrorBox(XtParent(w), "Variable for %s must be a loop variable",
		 role_list[selected_role]);
      else
	{
	  switch(selected_role)
	    {
	    case 0: new_curve.ix.index = ind; break;
	    case 1: new_curve.iy.index = ind; break;
	    case 2: new_curve.lstep = ind; break;
	    case 3: new_curve.lfaml = ind; break;
	    case 4: PromptBox(w, choice, ind); break;
	    }
	  if (selected_role != 4) ShowNewTitle();
	}
    }
  XtFree(choice);
}

/*-----------------------------------------------------------------------------
|	toggleCall
-----------------------------------------------------------------------------*/
void toggleCall(Widget w, XtPointer client_data, XtPointer call_data)
{
  int whichMenu, whichItem;
  whichItem = (int)client_data;
  whichMenu = whichItem / MSHIFT;
  whichItem -= whichMenu * MSHIFT;
  if (whichMenu==SING)
    {
      if (whichItem==0) new_curve.flags |= SINGLE;
      else new_curve.flags &= ~SINGLE;
    }
  else if (whichMenu==CONT)
    {
      new_curve.gtype = (whichItem==0) ? 'G' : 'C';
    }
}

static void xyCall(Widget w, XtPointer client_data, XtPointer call_data)
{
  printf("xyCall\n");
}

static void otherCall(Widget w, XtPointer client_data, XtPointer call_data)
{
  printf("otherCall\n");
}

/*-----------------------------------------------------------------------------
|	roleCall -- selected = one of the roles
-----------------------------------------------------------------------------*/
void roleCall(Widget w, XtPointer client_data, XtPointer call_data)
{
  XmListCallbackStruct *p;
  char *choice;
  p = (XmListCallbackStruct *)call_data;
  if (p->selected_item_count == 0)
    { selected_role = -1; return; }
  XmStringGetLtoR(p->item, XmFONTLIST_DEFAULT_TAG, &choice);
  selected_role = p->item_position - 1;
  printf("Selected role %d: %s\n", selected_role, choice);
  XtFree(choice);
}

/*-----------------------------------------------------------------------------
|	buttonCall -- selected = e.g. OK, Cancel
-----------------------------------------------------------------------------*/
void buttonCall(Widget w, XtPointer client_data, XtPointer call_data)
{
  int which, i;
  which = (int)client_data;
  printf("Button %d\n", which);
  if (which==0)
    {
      getparam(&new_curve);
      i = get_graph(0, selected_win);
      curveset[i] = new_curve;
      redrawflag = 1;
    }
  XtDestroyWidget(dialog);
}

/*=============================================================================
**		Dialog for read a parameter
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	CreateSelectionBox
-----------------------------------------------------------------------------*/
Widget CreateSelectionBox(Widget parent, int error, char *text)
{
  Widget sb, help;
  XmString str;

  if (error)
    {
      sb = XmCreateErrorDialog(parent, "error", NULL, 0);
      help = XmMessageBoxGetChild(sb, XmDIALOG_HELP_BUTTON);
      XtUnmanageChild(
	  XmMessageBoxGetChild(sb, XmDIALOG_CANCEL_BUTTON));
    }
  else
    {
      sb = XmCreatePromptDialog(parent, "prompt", NULL, 0);
      help = XmSelectionBoxGetChild(sb, XmDIALOG_HELP_BUTTON);
    }

  str = XmStringCreateLocalized(text);
  XtVaSetValues(sb, error? XmNmessageString : XmNselectionLabelString, str,
		XmNautoUnmanage, False, XmNtextColumns, 1, NULL);
  XmStringFree(str);  
  XtAddCallback(sb, error? XmNokCallback : XmNcancelCallback,
		    (XtCallbackProc)XtDestroyWidget, NULL);

  XtUnmanageChild(help);
  XtManageChild(sb);
  return sb;
}

/*-----------------------------------------------------------------------------
|	ErrorBox
-----------------------------------------------------------------------------*/
void ErrorBox(Widget parent, char *format, char *s)
{
  char text[100];
  Widget error;

  sprintf(text, format, s);
  error = CreateSelectionBox(parent, 1, text);
  XtPopup(XtParent(error), XtGrabNone);
}

/*-----------------------------------------------------------------------------
|	valueCall -- callback from prompt for OK
-----------------------------------------------------------------------------*/
void valueCall(Widget w, XtPointer client_data, XtPointer call_data)
{
  int i, ind;
  char *text;

  XmSelectionBoxCallbackStruct *p;
  p = (XmSelectionBoxCallbackStruct *)call_data;
  XmStringGetLtoR(p->value, XmFONTLIST_DEFAULT_TAG, &text);
  ind = (int)client_data;
  printf("Data: %s\n", text);
  sscanf(text, "%d", &i);
  XtFree(text);
  if (i >= 0 && i < loop[ind].count)
    {
      param[ind] = i;
      XtDestroyWidget(w);
      ShowNewTitle();
    }
  else ErrorBox(XtParent(w), "Value out of range", NULL);
}

/*-----------------------------------------------------------------------------
|	PromptBox
-----------------------------------------------------------------------------*/
void PromptBox(Widget parent, char *choice, int ind)
{
  Widget prompt;
  char text[100];
  int i;
  
  i = ind & INDBITS;
  sprintf(text, "Value for %s [0..%d]", choice, loop[i].count-1);
  prompt = CreateSelectionBox(parent, 0, text);
  XtAddCallback(prompt, XmNokCallback, valueCall, (XtPointer)i);
  XtPopup(XtParent(prompt), XtGrabNone);
}

