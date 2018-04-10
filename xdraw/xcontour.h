/*-----------------------------------------------------------------------------
|	xcontour.h
-----------------------------------------------------------------------------*/
extern  void redraw1(CURVE_SET *cp,
		     float xmin, float xmax, float ymin, float ymax );
extern	int new_ncurve(CURVE_SET *cp, char how);
extern  void contour_values(CURVE_SET *cp);
extern	void get_contlim(float *xlim, float *ylim);

#ifdef XCONTOUR
extern  void drawcontour(float psi,int ii);
extern  FLAG *surround(FLAG *ilist, POINTS *p);
extern	POINTS_ *has_psi(float psi, int ir, int iz, POINTS_ *q);
extern	POINTS_ *same_as(POINTS_ *p);
extern	POINTS_ *nextpoint(int k,POINTS_ *q,POINTS_ *qln);
extern	double anglecmp(double a1, double a2);
extern	POINTS_ *save_linear(int ir,int mr,int n,QCOEFF *kp,POINTS_ *q);
extern	POINTS_ *save_extremum(int ir,int mr,int n,float psi,int how,
	QCOEFF *k,POINTS_ *q);
extern	void getquadcoeff(int n,int ir,int mr, QCOEFF *k);
extern  int getp00(int,int);
extern	float quadratic(QCOEFF *k, float psi, float x0);
extern	void savelevel(int how, int index, char *format);
extern  void savegrid(int ir,int iz);
extern void psicount(float *buf, long bufsize,
		     float psi0,float dpsi,int ncurve,float *count);
#endif


