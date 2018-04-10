/*************************************************************************
**  NAME	smaller.c
**  USAGE	'smaller infile writef'
**		  writef: 0=report, 1=write
**
**  DESCRIPTION
**	1. Read data/esoln.bin (128 curves, length 64, 9 var per rec)
**	2. Output esoln2.bin (24 curves, 3 var per rec,
**	   length 8 @ 16, 8 @ 32, 8 @ 64
*************************************************************************/
#include <stdio.h>
#include <ctype.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int init(int argc, char *argv[]);

int writef = 0, in = 0, out = 0;
int nvar_new = 2, var_index[10], nfamily;
/*---------------------------------------------
|	main
---------------------------------------------*/
void main(int argc, char *argv[])
{
  int i, j, k, m, n, nvar, ncurve, length, length0, nread, writeme;
  int total_curves, ncurves_written;
  long size, count, buf[20], nullrec[2], newcount;
  static int curve_length[] = { 8, 16, 32 };
  static int nskip = 1, var_per_rec = 3, curves_per_length = 3;

  if (!init(argc, argv))
    {
      printf("Init error\n"); return;
    }
  printf("File handle = %d\n", in);
  if (writef)
    total_curves = nfamily;

  size = (long)lseek(in, 0L, SEEK_END);
  printf("File size %ld\n", size);
  lseek(in, 0L, SEEK_SET);

  nvar = ncurve = length = length0 = 0;
  nullrec[0] = nullrec[1] = 0;
  ncurves_written = 0;
  for(i=0; i < size; i+= count+8)
    {
      nread = read(in, &count, 4L);
      if (nread==0) { printf("Error: read 0 bytes, i=%d\n", i); break;}
      m = count / 4;
      if (i==0) { nvar = m; printf("%d variables\n", nvar); }
      else if (m && m != nvar) printf("Error! nvar=%d, m=%d\n", nvar, m);
      if (m) nread = read(in, buf+1, count);
      if (nread==0) { printf("Read error\n"); break; }
      read(in, &count, 4L);

      if (writef && ncurve < nfamily)
	{
	  if (count==0)
	    {
	      write(out, nullrec, 8);
	      ncurves_written++;
	      if (ncurves_written == nfamily) break;
	    }
	  else
	    {
	      for(i=1; i<=nvar_new; i++) buf[i] = buf[var_index[i]+1];
	      buf[0] = buf[nvar_new+1] = nvar_new * 4;
	      write(out, buf, buf[0]+8);
	    }
	}

      if (count == 0)
	{
	  printf("Curve %d: length %d\n", ncurve, length);
	  ncurve++;
	  length0 = length; length = 0;
	}
      else length++;
    }
  close(in);  
  if (writef) close(out);
}

/*------------------------------------------------------------------
|	init
------------------------------------------------------------------*/
int init(int argc, char *argv[])
{
  char text[80], *p;
  int i;
  strcpy(text, "data/");
  p = text + strlen(text);
  if (argc < 2)
    {
      printf("Enter file name: ");
      scanf("%s", p);
      printf("Enter 0=no write, 1=write: ");
      scanf("%d", &writef);
    }
  else
    {
      sscanf(argv[2], "%d", &writef);
      strcpy(p, argv[1]);
    }
  in = open(text, O_RDONLY, 0);
  if (in == -1)
    { printf("%s not found\n", text); return 0; }
  if (writef)
    {
      strcat(text, ".sml");
      printf("Output file is %s\n", text);
      out = open(text, O_RDWR|O_CREAT|O_TRUNC);
      printf("Enter desired number of variables: ");
      scanf("%d", &nvar_new);
      for(i=0; i<nvar_new; i++)
	{
	  printf("Index for variable %d: ", i);
	  scanf("%d", var_index+i);
        }
      printf("Enter desired number of family members: ");
      scanf("%d", &nfamily);
    }
  return 1;
}
