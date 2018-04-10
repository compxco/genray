/*-----------------------------------------------------------------------------
|	glx.h
-----------------------------------------------------------------------------*/
extern void init_glx(Display *, Window, unsigned long, unsigned long,
		     Font, XFontStruct *, XEvent *, int, int);
extern void prefposition(int mx1, int mx2, int my1, int my2);
extern long winopen(char *p);
extern void minsize(long dx, long dy);
extern void winconstraints(void);
extern void maxsize(long,long);
extern void minsize(long,long);
extern int windepth (long win);
extern void winpop(void);
extern void winpush(void);
extern long winget(void);
extern void winset(long win);

extern void cpack(long color);
extern long glx_index(long color);
extern void clear(void);
extern void rectfi(int x1, int y1, int x2, int y2);
extern void recti(int x1, int y1, int x2, int y2);

extern void cmov2i(int x, int y);
extern void v2i(long xy[]);
extern void drawstr(char *s);
extern void bgnpolygon(void);
extern void endpolygon(void);
extern void bgnline(void);
extern void endline(void);

extern void popmatrix(void);
extern void getorigin(long *x, long *y);
extern void qdevice(int dev);
extern int getvaluator(int dev);

extern int getdescender(void);
extern int getheight(void);
extern void scale_font(int scale);
extern long getstrwidth(char *text);

#ifdef NOT_YET
zap_events
add_event
#endif
