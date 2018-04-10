/*-----------------------------------------------------------------------------
|	xdraw.h
|	If your "vanilla" c-compiler doesn't like prototypes,
|	you must compile with acc or gcc.
-----------------------------------------------------------------------------*/
extern  void main(int argc,char * *argv);
extern  int get_graph(Window w, int i0);
extern  void getscale(int wx0, int dwx, float xmin, float xmax,
                      float *xscale, float *xoffset);
extern	float get_loopval(LOOP *qs, int i);
extern	int  get_box_info(Window, XRectangle *, XRectangle *,
			  float *, float *, float *, float *);
extern	int  get_zoomed_limits(CURVE_SET *, float *, float *, float *, float *);
extern	void get_aspected_limits(Window, int, int,
				 float *, float *, float *, float *,
				 float,float,float,float);
extern  void draw_frame(void);
extern  void redraw(void );
extern  void redraw0(int, int, int, float, float, float, float);
extern  void redraw00(int, int, int, float *, float *, float *, float *);
extern	void drawmarker(int x1,int y1,float fx1,float fy1);
extern	void setmarkerstyle(char *);
extern  void drawline(float, float, float, float, int, int, int);
extern	void draw_label(int, int, int, char *);
extern	void draw_misc_labels(CURVE_SET *);
extern  void draw_graphtype(char *type);
extern	void toggle_markers(CURVE_SET *cp, int *m);
extern	void toggle_aspect(CURVE_SET *cp, int *a);
extern	void toggle_extrema_use(CURVE_SET *cp1);
extern	void toggle_fill(CURVE_SET *cp);
extern	void toggle_single(CURVE_SET *cp, char c);
extern	void toggle_flabel(CURVE_SET *cp);
extern  char *skip_to_arg(char *text,int n,char delim);
extern	void init(int argc,char *argv[]);
extern void label(float *xlim, float xscale, float xoffset,
	  float *ylim, float yscale, float yoffset, CURVE_SET *cp);
extern int load_single_label(CURVE_SET *cp, char *string);
extern int rand_label(int);
