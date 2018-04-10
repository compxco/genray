/*************************************************************************
**  NAME	wtest.c
**  AUTHOR	Sheryl M. Glasser
**                                                   
**  DESCRIPTION
**	Test program for Windows
*************************************************************************/

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtest.h"

int PASCAL WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance,
		    LPSTR lpszCmdLine, int nCmdShow );
LRESULT CALLBACK WndEvent(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam);
BOOL FAR PASCAL _export MsgEvent(HWND, UINT, UINT, LONG);
void winRegister(const char *, int);
HWND CreateMainWindow(char *, char *, int, int);
void gPaint(void);
void dPaint(void);
void DrawCrosshair(HWND, int);


HWND grafWnd, dialogWnd;
HDC currentDC;

int xscreen, yscreen, caption, xbord, ybord, xframe, yframe;

int xcurs=0, ycurs=0, crosshair_enabled = 0;

BYTE xhairHbits = 0xff, xhairVbits = 0x80;
HBITMAP xhairH, xhairV;

char *dtext;
char not_gwin[] = "Cursor not in graphics window." ;
char text[80];
HINSTANCE hInstance, hPrevInstance;
int nCmdShow;

//-------------------------------------------
//	WinMain
//-------------------------------------------
int PASCAL WinMain( HINSTANCE hInst, HINSTANCE hPrev,
		    LPSTR lpszCmdLine, int cmdShow )
{
   MSG msg;
   hInstance = hInst;
   hPrevInstance = hPrev;
   nCmdShow = cmdShow;
   dtext = not_gwin;
   xhairH = CreateBitmap(8, 1, 1, 1, &xhairHbits);
   xhairV = CreateBitmap(1, 1, 1, 1, &xhairVbits);

   winRegister("wtest_g", BLACK_BRUSH);
   winRegister("wtest_d", WHITE_BRUSH);
   grafWnd =   CreateMainWindow("Graphics Window",  "wtest_g", 0, 2);
   dialogWnd = CreateMainWindow("Dialog and Menus", "wtest_d", 1, 2);
   UpdateWindow(grafWnd);
   UpdateWindow(dialogWnd);

   while (GetMessage(&msg, NULL, 0, 0))
   {
      TranslateMessage(&msg);
      DispatchMessage(&msg);
   }
   return msg.wParam;
}

//----------------------------------------------------
//	winRegister -- called from XOpenDisplay
//----------------------------------------------------
void winRegister(const char *name, int brush)
{
   WNDCLASS wndclass;

   xscreen = GetSystemMetrics(SM_CXFULLSCREEN);
   yscreen = GetSystemMetrics(SM_CYFULLSCREEN);
   caption = GetSystemMetrics(SM_CYCAPTION);
   xbord   = GetSystemMetrics(SM_CXBORDER);
   ybord   = GetSystemMetrics(SM_CYBORDER);
   xframe  = GetSystemMetrics(SM_CXFRAME);
   yframe  = GetSystemMetrics(SM_CYFRAME);
   if (!hPrevInstance)
   {
      wndclass.style         = CS_HREDRAW | CS_VREDRAW;
      wndclass.lpfnWndProc   = WndEvent;
      wndclass.cbClsExtra    = 0;
      wndclass.cbWndExtra    = 0;
      wndclass.hInstance     = hInstance;
      wndclass.hIcon         = LoadIcon( NULL, IDI_APPLICATION );
      wndclass.hCursor       = LoadCursor(NULL, IDC_ARROW);
      wndclass.hbrBackground = (HBRUSH)GetStockObject( brush );
      wndclass.lpszMenuName  = NULL;
      wndclass.lpszClassName = name;

      if ( ! RegisterClass( &wndclass ) )  exit(0);
   }
}

//-----------------------------------------------------------------
//	CreateMainWindow 
//-----------------------------------------------------------------
HWND CreateMainWindow(char *title, char *wclass, int has_menu, int maxWin)
{
   HWND hWnd;
   HMENU hMenu;
   static int nwin=0;
   static int nrow=3, ncol=3, dx, dy;
   int x, y, row, col;

   hMenu = has_menu ? LoadMenu(hInstance, "wtestMenu") : 0;
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
   hWnd = CreateWindow( wclass, title, WS_OVERLAPPEDWINDOW,
			x, y, dx, dy, NULL, hMenu, hInstance, NULL);
   if (!hWnd ) exit(0);
   ShowWindow(hWnd, nCmdShow);		// Display Window on screen
   return hWnd;
}

//-------------------------------------------------------
//	CreateMainDC
//-------------------------------------------------------
HDC CreateMainDC(HWND hwnd, HPEN hpen, HBRUSH hbrush)
{
   HDC hdc;
   hdc = GetDC(hwnd);
   SelectObject(hdc, GetStockObject(hpen));
   SelectObject(hdc, GetStockObject(hbrush));
   return hdc;
}

//-------------------------------------------------------
//	WndEvent
//-------------------------------------------------------
LRESULT CALLBACK WndEvent(HWND hWnd, UINT msgId, WPARAM w, LPARAM l )
{
   PAINTSTRUCT ps;
   POINT pt;
   switch (msgId)
   {
      case WM_PAINT:
	 if (hWnd != grafWnd && hWnd != dialogWnd) break;
	 currentDC = BeginPaint(hWnd, &ps);
	 if (hWnd == grafWnd) gPaint();
	 else if (hWnd == dialogWnd) dPaint();
	 EndPaint(hWnd, &ps);
	 return 0;

      case WM_MOUSEMOVE:
	 if (!crosshair_enabled) break;
	 //if (hWnd != grafWnd) break;
	 GetCursorPos(&pt);
	 xcurs = pt.x;
	 ycurs = pt.y;
	 if (hWnd != grafWnd) dtext = not_gwin;
	 DrawCrosshair(hWnd, 2);
	 return 0;

      case WM_DESTROY: PostQuitMessage( 0 ); return 0;

      case WM_COMMAND:
	 switch(w)
	 {
	    case WT_QUIT:
		SendMessage(hWnd,WM_CLOSE,0,0L); return 0;
	    case WT_ARROW:
		if (crosshair_enabled)
		{
		   DrawCrosshair(grafWnd, 0);
		   crosshair_enabled = 0;
		}
		return 0;
	    case WT_CROSS:
		if (!crosshair_enabled)
		{
		   crosshair_enabled = 1;
		   DrawCrosshair(hWnd, 1);
		}
		return 0;

	    default: break;
	 }
      default: break;
   }

   return DefWindowProc( hWnd, msgId, w, l );
}

//-------------------------------------------------------
//	MsgEvent -- messages from a dialog window
//-------------------------------------------------------
BOOL FAR PASCAL _export MsgEvent(HWND hWnd, UINT msgId, WPARAM w, LPARAM l )
{
   l++;
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
   return FALSE;
}

//-------------------------------------------
//	gPaint
//-------------------------------------------
static void gPaint()
{
   HDC hdc, hdcMem;
   HPEN hPen;
   HBRUSH hBrush;
   HRGN hRgn;
   char hello[100];
   static int x1=20, y1=20, x2=150, y2=140;

   hdc = currentDC;
   sprintf(hello, "GRAPHICS WINDOW!  Size %d,%d", xscreen, yscreen);
   TextOut(hdc, 0, 0, hello, strlen(hello));

   hBrush = GetStockObject(GRAY_BRUSH);
   SelectObject(hdc, hBrush);
   Rectangle(hdc, x1,y1,x2,y2);
   DeleteObject(hBrush);

   hRgn = CreateRectRgn(x1,y1,x2,y2);
   //SelectObject(hdc, hRgn);			// Enable clipping

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
   hdcMem = CreateCompatibleDC(currentDC);
   SelectObject(hdcMem, xhairH);
   StretchBlt(currentDC, 0, 100, 200, 1,
	      hdcMem,    0,   0, 8, 1, SRCINVERT);
   SelectObject(hdcMem, xhairV);
   StretchBlt(currentDC, 150, 0, 1, 200,
	      hdcMem,    0,   0, 1, 1, SRCINVERT);
}

//-------------------------------------------
//	mPaint
//-------------------------------------------
static void dPaint()
{
   HDC hDC;
   char hello[100];

   sprintf(hello, "Hello, World!");
   hDC = currentDC;
   SelectObject(hDC, GetStockObject(ANSI_VAR_FONT));
   TextOut(hDC, 0, 0, hello, strlen(hello));
   TextOut(hDC, 0,(dtext==not_gwin)?20:40, dtext, strlen(dtext));
}


//------------------------------------------------------
//	DrawCrosshair -- 0=disable,1=enable,2=move
//------------------------------------------------------
void DrawCrosshair(HWND hWnd, int draw)
{
   static int drawn;
   static int oldxcurs=0, oldycurs=0;
   int x, y;
   PAINTSTRUCT ps;
   HDC hdcMem;
   RECT rect, rectw;

   GetClientRect(grafWnd, &rect);
   GetWindowRect(grafWnd, &rectw);
   x = xcurs - rectw.left - xframe;
   y = ycurs - rectw.top  - caption - yframe;
   if (hWnd == grafWnd)
   {
       sprintf(text,
       "Cursor at %d %d! Window at %d %d, size %d %d",
	  x, y, rectw.left, rectw.top, rect.right, rect.bottom);
       dtext = text;
   }
   InvalidateRect(dialogWnd, NULL, TRUE);	// message

   if (draw==0 && !drawn) return;
   else if (draw==1)
   {
      drawn = (hWnd == grafWnd);
      if (!drawn) return;
   }
   else if (draw==2 && !drawn && hWnd == grafWnd)
   {
      DrawCrosshair(hWnd, 1); return;
   }
   else if (draw==2 && hWnd != grafWnd) return;
   
   currentDC = BeginPaint(grafWnd, &ps);
   hdcMem = CreateCompatibleDC(currentDC);

   SelectObject(hdcMem, xhairH);
   if (draw != 1)
      StretchBlt(currentDC, 0, oldycurs, rect.right, 1,
	 hdcMem, 0, 0, 8, 1, SRCINVERT);
   if (draw != 0)
      StretchBlt(currentDC, 0, y, rect.right, 1,
	 hdcMem, 0, 0, 8, 1, SRCINVERT);

   SelectObject(hdcMem, xhairV);
   if (draw != 1)
      StretchBlt(currentDC, oldxcurs, 0, 1, rect.bottom,
	 hdcMem, 0, 0, 1, 1, SRCINVERT);
   if (draw != 0)	 
      StretchBlt(currentDC, x, 0, 1, rect.bottom,
	 hdcMem, 0, 0, 1, 1, SRCINVERT);
   if (draw) { oldxcurs = x; oldycurs = y; }	 
   EndPaint(grafWnd, &ps);
   DeleteDC(hdcMem);
}

