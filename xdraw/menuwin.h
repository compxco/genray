/******************************************************************************
**  NAME	MENUWIN.H
**  AUTHOR	Sheryl M. Glasser
**
**  DESCRIPTION
**	Interface for menus
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
******************************************************************************/
extern long newmenu(void);
extern void addtomenu(long id, char *text, int *arg, long *submenu);
extern void enablemenus(long mmain, char *title, int x, int y,
			void (*func) (char *));
extern void enable_menus(long mmain, char *title, int x, int y,
			 int ndialog, void (*func) (char *));
extern int win_is_menu(long w);
extern void set_menu_visibility(int i);
extern int menu_event(int dev, short int state);
extern void menu_ievent(long newwin);
extern long get_cursorwin(void);
extern void menu_backup(void);
extern void menu_root(void);
extern void redrawmenu(void *, int);
extern void set_menubkgd(int, int, int, int);
extern void set_menufore(int, int, int, int);
extern unsigned long get_menubkgd(void);
extern void raise_menu(void);
extern void lower_menu(void);

#define ANY	-1
#define UP	0
#define DOWN	1


