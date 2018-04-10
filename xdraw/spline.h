/*-----------------------------------------------------------------------------
|	spline.h
-----------------------------------------------------------------------------*/
extern int init_spline(int, int);
extern void add_to_spline(int, float, float);
extern int eval_spline(int, float, float *, float *, int);
extern void fit_spline(void);
extern void zap_spline(void);
