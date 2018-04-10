/***************************************
*	wdraw.h
***************************************/
#define WDR_QUIT	1
#define WDR_PS		2
#define WDR_MARKER	3

extern HWND CreateMainWindow(char *, int, int);
extern LRESULT CALLBACK WndEvent(HWND hWnd, UINT msgId, WPARAM w, LPARAM l );
extern void winRegister(const char *);
extern void winDummy(void);
extern void winEvent(void);
