/*************************************************************************
**  NAME	wdraw.c
**  AUTHOR	Sheryl M. Glasser
**                                                   
**  DESCRIPTION
**	WinMain for Windows version of XDraw
**	Link also with xwinlib.c
**	Calls main(), which calls XOpenDisplay(), which inits Windows
*************************************************************************/

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xlib.h"
#include "xutil.h"
#include "xatom.h"
#include "wdraw.h"

int PASCAL WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance,
		    LPSTR lpszCmdLine, int nCmdShow );
LRESULT CALLBACK WndEvent(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam);
BOOL FAR PASCAL _export MsgEvent(HWND, UINT, UINT, LONG);
extern void createXEvent(HWND, int, int);
extern void event(void);

static HWND dialogWnd;
static int sx, sy;
static HINSTANCE hInstance, hPrevInstance;
static int nCmdShow;

WPARAM lastWParam = 0;
char appName[80];

extern int exitflag, xscreen, yscreen;
//-------------------------------------------
//	WinMain
//-------------------------------------------
int PASCAL WinMain( HINSTANCE hInst, HINSTANCE hPrev,
		    LPSTR lpszCmdLine, int cmdShow )
{
   int argc;
   char *argv[20], text[100], *p, *q;
   extern int main(int argc, char *argv[]);

   strcpy(text, lpszCmdLine);
   argc = 0; argv[argc++] = lpszCmdLine;
   for(p=text; *p; p+=strlen(p))
   {
      q = strchr(++p, ' ');
      if (q) *q = 0;
      argv[argc++] = p;
   }
   hInstance = hInst;
   hPrevInstance = hPrev;
   nCmdShow = cmdShow;
   main(argc, argv);
   return lastWParam;
}

//----------------------------------------------------
//	wRegister -- called from XOpenDisplay
//----------------------------------------------------
void winRegister(const char *aName)
{
   WNDCLASS wndclass;

   strcpy(appName, aName);
   if (!hPrevInstance)
   {
      wndclass.style         = CS_HREDRAW | CS_VREDRAW;
      wndclass.lpfnWndProc   = WndEvent;
      wndclass.cbClsExtra    = 0;
      wndclass.cbWndExtra    = 0;
      wndclass.hInstance     = hInstance;
      wndclass.hIcon         = LoadIcon( NULL, IDI_APPLICATION );
      wndclass.hCursor       = LoadCursor( NULL, IDC_ARROW );
      wndclass.hbrBackground = (HBRUSH)GetStockObject( WHITE_BRUSH );
      wndclass.lpszMenuName  = NULL;
      wndclass.lpszClassName = appName;

      if ( ! RegisterClass( &wndclass ) )  exit(0);
   }
}

//-----------------------------------------------------------------
//	CreateMainWindow -- called from xCreateMainWindow
//-----------------------------------------------------------------
HWND CreateMainWindow(char *title, int has_menu, int maxWin)
{
   HWND hWnd;
   HMENU hMenu;
   static int nwin=0;
   static int nrow=3, ncol=3, dx, dy;
   int x, y, row, col;

   hMenu = has_menu ? LoadMenu(hInstance, "wdrawMenu") : 0;
   if (nwin==0)
   {
      if (maxWin < 4) nrow = ncol = 2;
      else nrow = ncol = 3;
      dx = xscreen / ncol;
      dy = yscreen / nrow;
   }

   row = nwin / ncol;
   col = nwin - row * ncol;
   x = col * dx;
   y = row * dy;
   nwin++;
   hWnd = CreateWindow( appName, title, WS_OVERLAPPEDWINDOW,
			x, y, dx, dy, NULL, hMenu, hInstance, appName);
   if (!hWnd ) exit(0);
   ShowWindow(hWnd, nCmdShow);		// Display Window on screen
   					// (XMapRaised will paint)
   if (has_menu) dialogWnd = hWnd;
   return hWnd;
}

//----------------------------------------------------------------------
//	winEvent, winDummy == process_event(), expose_event()
//----------------------------------------------------------------------
void winEvent()
{
   MSG msg;
   while (GetMessage(&msg, NULL, 0, 0))
   {
      TranslateMessage(&msg);
      DispatchMessage(&msg);
   }
}

void winDummy() {}

//-------------------------------------------------------
//	WndEvent
//-------------------------------------------------------
LRESULT CALLBACK WndEvent(HWND hWnd, UINT msgId, WPARAM w, LPARAM l )
{
   int gotcha = 1;
   switch (msgId)
   {
      case WM_PAINT:
	 createXEvent(hWnd, Expose, 0); break;

      case WM_COMMAND:
	 switch(w)
	 {
	    case WDR_QUIT:
		SendMessage(hWnd,WM_CLOSE,0,0L); return 0;
	    case WDR_MARKER:
		createXEvent(hWnd, KeyPress, 'm'); break;

	    case WM_DESTROY: PostQuitMessage( 0 ); return 0;
	    default: gotcha = 0; break;
	 }
      default: gotcha = 0;
   }

   if (!gotcha) return DefWindowProc( hWnd, msgId, w, l );
   event();
   return 0;
}

//-------------------------------------------------------
//	MsgEvent -- messages from a dialog window
//-------------------------------------------------------
BOOL FAR PASCAL _export MsgEvent(HWND hWnd, UINT msgId, WPARAM w, LPARAM l )
{
   if (l==1) ;
   switch (msgId)
   {
		case WM_INITDIALOG:	return TRUE;
		case WM_COMMAND:
			switch (w)
			{
				case IDOK:
				case IDCANCEL:
					EndDialog(hWnd, 0); return TRUE;
			}
		break;
	}
	lastWParam = w;
	return FALSE;
}

//-------------------------------------------
//	gPaint
//-------------------------------------------
static void gPaint( HWND hWnd )
{
	PAINTSTRUCT ps;
	HDC hdc;
	HPEN hPen;
	HBRUSH hBrush;
	HRGN hRgn;
	char hello[100];
	static int x1=20, y1=20, x2=150, y2=140;

	hdc = BeginPaint( hWnd, &ps );
	sprintf(hello, "GRAPHICS WINDOW!  Size %d,%d", sx,sy);
	TextOut(hdc, 0, 0, hello, strlen(hello));

	hBrush = GetStockObject(BLACK_BRUSH);
	SelectObject(hdc, hBrush);
	Rectangle(hdc, x1,y1,x2,y2);
	DeleteObject(hBrush);

	hRgn = CreateRectRgn(x1,y1,x2,y2);
	SelectObject(hdc, hRgn);			// Enable clipping

	hPen = GetStockObject(WHITE_PEN);
	SelectObject(hdc, hPen);
	MoveTo(hdc, 30,30);
	LineTo(hdc, 100,30);
	DeleteObject(hPen);
	hPen = CreatePen(PS_SOLID, 3, RGB(0,255,255));
	SelectObject(hdc, hPen);
	LineTo(hdc, 100,100);
	DeleteObject(hPen);
	hPen = CreatePen(PS_SOLID, 1, RGB(0, 255, 0));
	SelectObject(hdc, hPen);
	LineTo(hdc, 200,200);

	DeleteObject(hPen);
	DeleteObject(hRgn);
	EndPaint( hWnd, &ps );
}

//-------------------------------------------
//	mPaint
//-------------------------------------------
static void mPaint( HWND hWnd )
{
	HDC hDC;
	PAINTSTRUCT ps;
	char hello[100];

	sprintf(hello, "Hello, World! ");
	hDC = BeginPaint( hWnd, &ps );
	TextOut(hDC, 0, 0, hello, strlen(hello));
	EndPaint( hWnd, &ps );
}



