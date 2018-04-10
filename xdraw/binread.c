/******************************************************************************
**  NAME      BINREAD.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Read in data from .in and .bin, load structures
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#define O_BINARY 0

#else
#include <sys\types.h>
#include <sys\stat.h>
#ifdef DOS
#include <graph.h>
#endif
#include <io.h>
#include <conio.h>
#include <dos.h>
#include "xlib.h"
#include "xutil.h"
#endif

#ifndef O_RDONLY
#define O_RDONLY 0			/* for penning */
#endif

#include "gendefs.h"
#include "curves.h"
#include "xinit.h"

static int  update_outerloop(void);
static void update_nodeloop(int);
static int start_new_rec(unsigned Long count, unsigned Long count0,
		 unsigned Long rec, unsigned Long rec0, int nhead);

extern int ftype;	/* 0=graph, 1=contour, 2=self-determined, 3,4 */

extern LOOP loop[];		/* Loop structure of the data */
extern int nloop;
extern int outerloop_added;
static int iloop;
static int counting_outer=0, i_outer=0;
extern int param[MAXLOOP];
static int use_hcount=0, stringbytes=0;
static int invert_bytes;

extern float *buf;
static unsigned long bufsize;

extern NODE *nodelist;
extern int nnode, inode, ivar, ncount_equal;
#define NBLOCK 20

/*=============================================================================
**		READ BIN FILES
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	init_type_g
|	n = # "outer" loops (parameterization levels ABOVE node level)
-----------------------------------------------------------------------------*/
void init_type_g(int, int);
void init_type_g( int n , int count)
{
  LOOP *q;
  int i;

  nloop = n + 2;
  nloop += outerloop_added;
  for (q=loop, i=0; i<nloop-2; i++,q++)
  {
    q->ltype = 'I'; q->count = 0; q->use_sep=0;   q->extra = NULL;
  }
  q->ltype = 'N'; q->count = 0;   q->use_sep='1'; q->extra=NULL; q++;
  q->ltype = 'V'; q->count=count; q->use_sep='1'; q->extra=NULL;
  nodelist = (NODE *)malloc( NBLOCK * sizeof(NODE) * n);
  nodelist->ncount = (ftype==0) ? 1 : 0;
  stringbytes = use_hcount = i_outer = 0;
  counting_outer = 1;
  inode = nloop - 2;
  ivar  = nloop - 1;
}

/*-----------------------------------------------------------------------------
|	initloops -- after read 1st rec
-----------------------------------------------------------------------------*/
static int initloops(Long *p, Long count)
{
  LOOP *q, *qr, *qz;
  int nhead;

  q = loop;
  nloop = 3;
  iloop = 0;
  ivar=-1;
  nnode = 0;

  if (ftype==0 || ftype==4)
    {
      init_type_g(ftype==0 ? 1 : (int)*p, (int)count);
      nhead = (ftype==0) ? 0 : 1;
    }

  else if (ftype==1)
    {
      if (count==3 && *(p+2)==0)	/* flag for time steps */
        {
	  q->ltype = 'I';
	  q->use_sep = '1';
	  q->count = 0;
	  q->extra = NULL;
	  q++;
	}
      qz = q; qr = q+1; q += 2;
      qr->ltype = qz->ltype = 'X';
      qr->use_sep  = qz->use_sep  = '1';
      qr->count = (int)(*p++);
      qz->count = (int)(*p);
      qr->extra = (float *)malloc(3 * sizeof(float));
      qz->extra = (float *)malloc(3 * sizeof(float));
      q->ltype = 'V'; q->count = 1; q->use_sep='1'; q->extra = NULL;
      nhead = 1; ivar = q - loop; nloop = ivar + 1;
      if (nloop == 4) counting_outer = qr->count * qz->count * q->count;
    }
  else					/* Type I: version, nloop, stringbytes*/
    {
      nloop = nhead = (int)*(p+1);
      if (*p==0 && count==3)
	  stringbytes = use_hcount = 0;
      else
	{
	  stringbytes = (int)*(p+2); use_hcount = 1;
	}
    }

  return(nhead);
}

/*-----------------------------------------------------------------------------
|	addloop -- end of header record, enter loop into structure
-----------------------------------------------------------------------------*/
static void addloop(Long *p)
{
  LOOP *q;
  float *f, *f1;
  int i, n;

  if (ftype==1)
    {
      f = (float *)p;
      q = loop + nloop - 2;
      f1 = q->extra;			/* rmin, rmax = inner loop */
      *f1++ = *f++;
      *f1++ = *f++;
      *f1   = (*(f1-1) - *(f1-2)) / (float)(q->count-1);

      q = loop + nloop - 3;
      f1 = q->extra;			/* zmin, zmax = outer loop */
      *f1++ = *f++;
      *f1++ = *f++;
      *f1   = (*(f1-1) - *(f1-2)) / (float)(q->count-1);
    }
  else if (ftype==4)
    {
      for(i=0,q=loop; i<inode; i++,q++)
      {
	q->count = (int)(*p++);
	if (i>0) counting_outer *= q->count;
	else q->count = 0;			/* ALWAYS count outer loop */
      }
    }

  if (ftype != 2) return;
  n=0;
  q = loop + iloop;
  q->ltype = invert_bytes ?			/* 1 char, but must take up 4 bytes */
             *(char *)p : *((char *)p +3);
  p++;

  q->use_sep = '1';
  q->count = (int)(*p++);		/* (type 'N' must include count=0) */
  q->hcount = q->ih0 = 0;

  i_outer = 0;
  if (iloop > 0) q->ih0 = (q-1)->ih0 + (q-1)->hcount;
  else if (iloop==0 && q->count==0) counting_outer=1;
  if (use_hcount) q->hcount = (in1)(*p++);

  q->extra = NULL;
  f1 = (float *)p;

  if (q->ltype=='X')
    {
      f = q->extra = (float *)malloc(3*sizeof(float));
      n = 2;
      *f++ = *f1++;		/* x1 */
      *f++ = *f1++;		/* x2 */
      *f   = *f1;		/* dx */
      if (*f==FZ && q->count > 1)
	*f = (*(f-1) - *(f-2)) / (q->count-1);
    }
  else if (q->ltype=='A')
    {
      f = q->extra = (float *)malloc(q->count*sizeof(float));
      n = q->count;
      for(i=0; i<n; i++) *f++ = *f1++;
    }
  else if (q->ltype=='N' && iloop>0) (q-1)->use_sep = 0;
  else if (q->ltype=='N') d_abort("Inappropriate level for Node loop",0,0);
  else if (q->ltype=='V' && ivar<0) ivar = iloop;
  else if (q->ltype=='V') d_abort("Loop contains too many type V",0,0);
  iloop++;
}

/*-----------------------------------------------------------------------------
|	update_outerloop -- for ftype==0 or 4, null record
-----------------------------------------------------------------------------*/
static int update_outerloop()
{
  LOOP *nq;
  int retval;
  i_outer++;
  if (i_outer!=counting_outer) retval=0;
  else { i_outer = 0; retval = 1; }
  if (inode==-1) return(retval);	/* see if exists node loop */
  nq = loop + inode;

  nnode++;
  if (nnode==1) nq->count = nodelist->ncount;
  if ((nnode % NBLOCK) == 0)
    nodelist = (NODE *)realloc(nodelist, (nnode+NBLOCK)*sizeof(NODE));
  (nodelist + nnode)->ncount = 0;
  return(retval);
}

/*-----------------------------------------------------------------------------
|	update_nodeloop - for ftype=0.  Update # nodes, at end of record
-----------------------------------------------------------------------------*/
static void update_nodeloop( int increase)
{
  if (increase) (nodelist + nnode)->ncount += 1;
  else		(nodelist + nnode)->ncount -= 1;
}

#ifndef UNIX
/*------------------------------------------------------------------------------
|	bigread -- for PC, read data in blocks
------------------------------------------------------------------------------*/
static void bigread(int fd, unsigned char *p, unsigned Long newsize)
{
  unsigned Long count;
  unsigned Long nread;
  Long off;

  for(count=newsize;;)
    {
      nread = (count < 0x10000L) ? count : 0xff00;
      off = FP_OFF(p);
      if (off+(Long)nread > 0x10000L)
      	{
      	FP_SEG(p) += (int)(off>>4);
      	FP_OFF(p) &= 0xf;
      	}
      nread = (unsigned Long)read(fd, p, (unsigned int)nread);
      p += nread;
      count -= nread;
      if (count==0 || nread==0) break;
    }
}
#endif

/*-----------------------------------------------------------------------------
|	go_writenew
|	* q non-NULL ==> writes 4 bytes at q
|	* n >= 0 ==> writes a header rec to outer loop, value n**2
|	* n = -1 ==> close file
|	* n = -2 ==> write L header records
-----------------------------------------------------------------------------*/
void go_writenew(int n, unsigned char *q)
{
#ifndef UNIX
   static int fw=0, count=0, nhead=0;
   float value;
   char *p;
   unsigned char buf[28];
   int i;
   
   if (q)
     {
       write(fw, q, 4L); count += 4;
       return;
     }
   if (n==-1)
     {
       close(fw); printf("%ld bytes written to temp.bin\n",count);
       return;
     }
   if (fw==0)
     fw = open("temp.bin", O_RDWR|O_BINARY|O_CREAT, S_IREAD|S_IWRITE);

   if (n==-2)
     {
       for(i=0; i<28; i++) buf[i]=0;
#ifdef UNIX			/* note this can't happen */
       buf[3] = buf[11] = 4;
       buf[7] = 2;
       buf[15] = buf[27] = 8;
       buf[23] = 5;
#endif
       write(fw, buf, 28L);
       count += 28;
       return;
     }

   nhead++;
   if (nhead > 10)
     {
       printf("Abort: Header count exceeds maximum\n");
       return;
     }
   printf("New header records, count = %d\n", count);
   for(i=0; i<20; i++) buf[i]=0;
   value = (float)n * (float)n;
   p = (char *)&value;
   buf[3] = buf[11] = 4;
   if (invert_bytes)
   {
      buf[4] = *(p+3); buf[5] = *(p+2);
      buf[6] = *(p+1); buf[7] = *p;
   }
   else
   {
      buf[7] = *(p+3); buf[6] = *(p+2);
      buf[5] = *(p+1); buf[4] = *p;
   }
   write(fw, buf, 20L);
   count += 20;
#endif
}

/*-----------------------------------------------------------------------------
|	binread. Read from a .bin file, handle fd.
|	*sizep = size of buf already allocated, in bytes
|	nh = number of records to be considered as header (G=0, C=2)
|
|	Format of the input data:
|	*  "words" are 4 bytes long: float or long int
|	*  records are: <nbytes>, ...<nbytes of data>..., <nbytes>
|	*  draw.c had eg. 3 parameters per record, many records per curve,
|	   a curve was terminated by a record with nbytes=0, use size
|	   of file to terminate
|	*  contour.c had 1 parameter, all data in the data record.
|	   other curves followed, but the program did not use.
|	   12, m,n,o, 12, 16, rmin..zmax, 16, m*n*4, data... , m*n*4, ..
|
|	Format of the output data:
|	* header records: <nbytes>, ...<nbytes of data>...
|	* data records: <nparam>,<nbytes>, ...<nbytes*nparam of data>...
|	* if multiple curves, end when nparam=0.
|
|	* contour: 0C, <mr,mz,mpsi>, 10, <rmin..zmax>,
|		   4, <nparam>, data, mparam, data,...
-----------------------------------------------------------------------------*/
float *binread(char *fname, unsigned Long *sizep)
{
  unsigned Long oldsize, newsize;
  float *buf0, temp;
  unsigned char *q;
  Int i;
  int fd, nhead, oflag, result;
  static int nhead0;
  unsigned Long rec, rec0, count, request;
  Long *p0, *longp, count0;
  int ntry;
#ifdef UNIX
  unsigned char *p;
  float *f;
  Int *data;
#else
  unsigned char huge *p;
  float huge *f, huge *fsrc;
  unsigned long adr1,adr2;
  static Int huge *data = NULL;
#endif
#define skip_endcount i+=4; p+=4; rec++; if (writenew) go_writenew(0,p); 
  static int debugg = 500;
  static int writenew=0;	/* 1 for bal5, 2 for bal4 */
  oflag = O_RDONLY;
#ifndef UNIX
  oflag |= O_BINARY;
#endif

  fd = open(fname, oflag, 0);				/* eg. bal3.bin */
  						/* O_BIN required(!) on PC */
  if (fd == -1)
    {							/* (can have many */
      printf("Abort: Cannot open %s.\n", fname);	/* input files?) */
      exit(0);
    }

/*------ determine file size, allocate buffer, and read binary data */
  oldsize = *sizep;
  newsize = lseek(fd, 0L, SEEK_END);
#ifndef UNIX
  /*newsize &= 0xffffL;*/		/* WARNING! because compiler does cwd! */
#endif
  newsize += 4;			/* one more word, in case need for <nparams> */
  ntry = 0;
  
  if (oldsize == 0)
#ifdef UNIX
    buf0 = (float *) malloc(newsize + 8);
#else
AGAIN:
    buf0 = (float *) halloc(newsize - ntry*1000 + 8, 1);
#endif
  else
    {
      request = oldsize + newsize + 8;
#ifdef UNIX
      buf0 = (float *) realloc(buf, (size_t)request);
#else
      buf0 = (float *) halloc(request,1);
      for(i=0,f=buf0,fsrc=buf; i<oldsize; i+=sizeof(float)) *f++ = *fsrc++;
      hfree(buf);
#endif
    }

  if (buf0 == 0)
    {
#ifndef UNIX
      if (ntry < 20) { ntry++; goto AGAIN; }
#endif
      printf("Insufficient memory for data.\n");
      exit(1);
    }
#ifndef UNIX
  if (ntry > 0) { printf("buf allocation reduced by %d,000\n", ntry);
    printf("Press any key..."); getchar(); fflush(stdin);
    }
#endif
  buf = buf0 + oldsize / sizeof(float);

/*------ Read ENTIRE file into buf+1 */
  lseek(fd, 0L, SEEK_SET);
  p = (unsigned char *)(buf+1);
#ifdef UNIX
  read(fd, p, newsize);
#else
  bigread(fd,p,newsize);
#endif
  data = (Int *)buf;			/* calc data 2 steps so huge ok */
  data += newsize / sizeof(float) - 2;
  if (*data != 0 || *(data + 1) != 0)	/* check for termination null rec */
    {
      *(data + 2) = 0;
      *(data + 3) = 0;
      newsize += 8;
    }

/*------- Parse / pack the data */
  request = *(Long *)(buf+1);
  invert_bytes = (request & 0xffff0000) != 0;
  f = buf;
  q = (unsigned char *) &temp;
  rec = rec0 = count = count0 = nhead = 0;

  if (writenew==1) go_writenew(0, NULL);	/* write value of psi0 */
  else if (writenew==2) go_writenew(-2, NULL);	/* write 2 headers */

  for (i=4, p=(unsigned char *)buf + i; i < newsize; i+=4, p+=4)
    {				/* (compiler ok on p+=4, bad on p=buf+i) */
#ifdef STRINGS
      if (rec<nhead)
	{
	  for(pc=(char *)f; *p; i++) *pc++ = *p++;
	  *pc++ = *p++;
	  datahead = f = (float *)pc;
	}
#endif

      if (invert_bytes)
        {
	  request = *(Long *)p;
          *(q + 0) = *(p + 3);	/* bytes are written reversed!! */
	  *(q + 1) = *(p + 2);
	  *(q + 2) = *(p + 1);
	  *(q + 3) = *(p + 0);
	  request = *(Long *)q;
        }
    
      else temp = *(float *) p;

      if (writenew) go_writenew(0, p);

      longp = (Long *)p;		/* debug; current "temp" */

      if (count == 0)			/* ----- if starting a new record */
	{
	  count = *(Int *)&temp / 4L;	/* get new count for record */
	  if (!(result = start_new_rec(count,count0,rec,rec0,nhead)))
	    break;
	  count0 = count;
	  p0 = (Long *)f;
	  if (count==0)
	    {
	      rec0 = rec;
	      skip_endcount;
	      if (result == 2 && writenew==1) go_writenew(nnode+1,NULL);
	    }
	  continue;
	}

      *f++ = temp;
      if (--count == 0)			/* ----- if ending record */
	{
	  if (rec==0)
	    {
	    if (oldsize==0) nhead0 = initloops(p0,count0);
	    else if (ftype==0) nodelist[nnode].ncount = 1;
	    nhead = nhead0;
            }
	  else if (rec <= nhead)
	    {
	      if (oldsize==0) addloop(p0);
	      if (rec == nhead)
	        {
		  f = buf;
		  if (ivar==-1) d_abort("No type V found in data file",0,0);
		}
	    }
	  else if (ftype==0 || ftype==4)	/* (generalize type N later) */
	    update_nodeloop(1);
	  else if (ftype==1 && !counting_outer) break;
	  skip_endcount;
	}
    }

  if (writenew) go_writenew(-1, NULL);

#ifdef UNIX
  *sizep = (f - buf0) * 4;
#else
  adr1 = FP_SEG(f);
  adr1 = (adr1<<4) + FP_OFF(f);
  adr2 = FP_SEG(buf0);
  adr2 = (adr2<<4) + FP_OFF(buf0);
  *sizep = adr1 - adr2;
#endif
  close(fd);
  return(buf0);
}

/*-----------------------------------------------------------------------------
|	start_new_rec
-----------------------------------------------------------------------------*/
static int start_new_rec(unsigned Long count, unsigned Long count0,
		unsigned Long rec, unsigned Long rec0, int nhead)
{
  static int need_varcount=0;
  int uses_nullrec;
  int j, retval;

  uses_nullrec = (ftype==0 || ftype==4);

  if (count==0 && (count0==0 || !uses_nullrec))	/* Double null or single */
    return 0;

  retval = 1;
  if (count==0)					/* single null on G,L */
    {
      if (rec>nhead && rec-rec0 <= 2)		/* if loop header... */
	{					/* ...add to loop[] */
	  for(j=0; j<inode && loop[j].hcount!=0; j++) ;
	  if (j<inode)
	    {
	      loop[j].hcount = (in1)count0;
	      if (j==0) loop[j].ih0 = 0;
	      else loop[j].ih0 = loop[j-1].ih0 + loop[j-1].hcount;
	      need_varcount=1;
	    }
	  update_nodeloop(0);
	}
      else					/* else count outer loop */
	{
	  if (update_outerloop()) loop->count++;
	  retval = 2;
	}
    }

  else if (uses_nullrec && need_varcount)
    {
      loop[ivar].count = (int)count;
      need_varcount = 0;
    }

  else if ( !uses_nullrec && rec>nhead && counting_outer &&
	   (ftype==1 || ftype==2) )
    loop->count++;
  
  if (rec==nhead+1 && ftype==4)		/* L give size of extra loop */
    (loop+nloop-1)->count = (int)count;
  
  return retval;
}

/*-----------------------------------------------------------------------------
|	loop_structure
-----------------------------------------------------------------------------*/
void loop_structure()
{
  int i, ncount;
  LOOP *lp;
  NODE *np;

  for(i=0,lp=loop+nloop-1; i<nloop; i++,lp--)
    {
      if (i==0) lp->sep = 1;
      else lp->sep = (lp+1)->sep * (lp+1)->count + lp->hcount;
      /*----- Note: any lp->sep above a type Node is invalid! */
      if (lp->ltype=='H')
        {
	  if ((lp+1)->ltype != 'N') lp->sep++;
	  else for(i=0; i<nnode; i++) (nodelist + i)->ncount += 1;
	}
    }

  if (inode != -1)
  {
    for(i=0,lp=loop+inode,np=nodelist,ncount_equal=1; i<nnode; i++,np++)
    {
      if (i==0) ncount = np->ncount;		/* See if equal count */
      else if (np->ncount != ncount) ncount_equal = 0;
      np->nsep = np->ncount * lp->sep;		/* Obtain nsep */
      np->nsep += (lp-1)->hcount;		/* Kludge */
    }
    if (ncount_equal)
    {
      for(i=0,lp=loop; i<inode; i++,lp++) lp->use_sep = '1';
      loop[inode].ltype = 'I';
      inode = -1;
    }
  }
}


