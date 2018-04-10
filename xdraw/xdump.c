/******************************************************************************
**  NAME            XDUMP.C
**  AUTHOR          Sheryl M. Glasser
**
**  DESCRIPTION
**
**
**  Copyright (c) GlassWare 1994.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <fcntl.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#define O_BINARY 0

#else
#include <io.h>
#endif

extern void showblock(int i0, int ni);

char formcode;
int fd;
long ntotal;
long buf[1024];
int nullrec = 0;

/*-----------------------------------------------------------------------------
|	
-----------------------------------------------------------------------------*/
main(int argc, char *argv[])
{
  char fname[100];
  int n, ni, i0;
  long m;
  char action;

  if (argc < 2)
    {
      printf("Enter input file name: ");
      scanf("%s", fname);
    }
  else strcpy(fname, argv[1]);
  fd = open(fname, O_RDWR|O_BINARY, 0);
  if (fd==-1) { printf("Cannot open %s\n",fname); exit(0); }
  
  ntotal = lseek(fd, 0L, SEEK_END);
  printf("------- %s has %ld bytes (%lx hex)\n", fname, ntotal, ntotal);
  lseek(fd, 0L, SEEK_SET);
  ntotal /= sizeof(long);

  if (argc >= 4)
    {
      sscanf(argv[2], "%d", &ni);
      sscanf(argv[3], "%c", &formcode);
    }
  else
    {
      printf("Enter # lines per output block (decimal): ");
      scanf("%d", &ni); fflush(stdin);
      printf("Enter x=hex, f=float: ");
      scanf("%c", &formcode); fflush(stdin);
    }
  ni *= 8;

  for(i0=0; i0<(int)ntotal; i0+=ni)
    {
      showblock(i0, ni);
      for(;;)
	{
	  printf("n=next block, c=change a data value, q=quit: ");
	  scanf("%c", &action);
	  fflush(stdin); fflush(stdout);
	  if (action=='n' || action=='c' || action=='q') break;
	}
      if (action=='q') break;
      if (action=='n') continue;
      printf("Enter offset to change (hex bytes), or -1: ");
      scanf("%x",&n);
      if (n >= 0)
	{
	  printf("Enter new value (hex), or -1 [%lx]: ", buf[n/4]);
	  scanf("%x", &m);
	  if (m != -1)
	    {
	      lseek(fd, (long)n, SEEK_SET);
	      write(fd, (char *)&m, 4);
	    }
	}
    }
  close(fd);
  printf("%d null records encountered\n", nullrec);
  return 0;
}

/*-----------------------------------------------------------------------------
|	showblock
|	* ni = # longs per block (multiple of 8)
|	* i0 = offset of first long
-----------------------------------------------------------------------------*/
void showblock(int i0, int ni)
{
  int i,j,n, newrow, count;
  long m;
  float v;
  long val;
  unsigned char *p, *q;

  n =  (i0+ni-1 < (int)ntotal ) ? ni : (int)ntotal - i0;
  for(i=count=0; i<n; i+=512)
    {
      m = (i+512 <= n) ? 512 : n-i;
      read(fd, (char *)buf, m * (long)sizeof(long));
      newrow = 1;
      for(j=0; j<(int)m; j++)
	{
#ifdef UNIX
	  if (newrow) printf("%06x: ", (i+j+i0)*sizeof(long));
	  val = buf[j];
#else
	  if (newrow) printf("%05x: ", (i+j+i0)*sizeof(long));
	  p = (unsigned char *)&val;
	  q = (unsigned char *)&buf[j];
	  *p=*(q+3); *(p+1)=*(q+2); *(p+2)=*(q+1); *(p+3)=*q;
#endif
	  if (count==0 || count==1)
	    {
	      if (count==0) count = (int)val/sizeof(long) + 1;
	      else count--;
	      if (val != 0L) printf("**  %04x ", val);
	      else         { printf("........ "); nullrec++; }
	    }
	  else
	    {
	      v = *(float *)&buf[j];
	      if (formcode!='f') printf("%08x ", val);
	      else printf("%8.3f ", v);
	      count--;
	    }
	  newrow=0;
	  if (j==(int)m-1 || (j%8==7))
	    {
	      printf("\n"); newrow=1;
	    }
	}
    }
}
