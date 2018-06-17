      subroutine contour2d(nxmax,nymax,f,x,y,plot_name,x_name,y_name,
     .n_contour,contour,name_param,param,n_param)
c-----print contours of 2D array f(1:nxmax,1:nymax)
c
c     input
c      fmap(1:nxmax,1:nymax) 2D array of the function f(x,y)
c      nxmax,nymax are the dimentions of f
c      x(1:nxmax)  array of x variable
c      y(1:nymax)  array of y variable
c      plot_name  the name of the top of plot
c      x_name     the name of x axis  
c      y_name     the name of y axis
c      n_contour  the number of contours
c      contour(n_contour)      work array for the contour's values
c      n_param is the number of parameters
c      name_param(n_param)     are the names of parameters
c      param(n_param)          are the values of parameters

      implicit none
c-----input
      integer nxmax,nymax,n_contour,n_param
      real f(nxmax,nymax),x(nxmax),y(nymax),contour(n_contour)
      character*(*)plot_name,x_name,y_name
      character*(*)name_param(n_param)
      real param(*)
      real y_e,x_e,cnpar,T
c-----external
      real rbound !convert double to real bounded 1.e-33 and 1.e+33
c-----local
      real tr(6),RILIN
      real xmin,xmax
      real ymin,ymax
      real fmin,fmax
      integer i,j
      
c      CHARACTER*7 LABEL
      CHARACTER*2 LABEL

      integer INTVAL,MININT
      character*72 text    
               
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
         WRITE(*,*)'contour2d: TRYING TO PGPAGE at rank=',myrank
         STOP
      endif     

c-----Advance plotter to a new page or panel, clearing the screen if
c     necessary.
      CALL PGPAGE  

c-----Change the size and position of the viewport, specifying
c     the viewport in normalized device coordinates.  Normalized
c     device coordinates run from 0 to 1 in each dimension.
c     SUBROUTINE PGSVP (XLEFT, XRIGHT, YBOT, YTOP)
c     REAL XLEFT, XRIGHT, YBOT, YTOP

c      CALL PGSVP(.2,.8,.3,.95) 
      CALL PGSVP(.2,.8,.5,.95) 

c-----Change the window in world coordinate space that is to be mapped on
c     to the viewport.  Usually PGSWIN is called automatically by PGENV,
c     but it may be called directly by the user.
c     SUBROUTINE PGSWIN (X1, X2, Y1, Y2)
c     REAL X1, X2, Y1, Y2
      call rminmx(x,1,nxmax,1,xmin,xmax)
      if (xmin .eq. xmax) xmin=.9*xmax-1.e-20
      call rminmx(y,1,nymax,1,ymin,ymax)
      if (ymin .eq. ymax) ymin=.9*xmax-1.e-20
      CALL PGSWIN (xmin,xmax,ymin,ymax)

c-----PGBOX -- draw labeled frame around viewport
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)

c-----PGLAB write labels  for x,y axis and top of plot
      CALL PGLAB(x_name,y_name,plot_name)
 
c-----draw a contour map of an array.  The map is truncated if
c     necessary at the boundaries of the viewport.  Each contour line
c     is drawn with the current line attributes (color index, style, and
c     width); except that if argument NC is positive (see below), the line
c     style is set by PGCONT to 1 (solid) for positive contours or 2
c     (dashed) for negative contours.
c
c     Arguments:
c     A      (input) : data array.
c     IDIM   (input) : first dimension of A.
c     JDIM   (input) : second dimension of A.
c     I1, I2 (input) : range of first index to be contoured (inclusive).
c     J1, J2 (input) : range of second index to be contoured (inclusive).
c     C      (input) : array of NC contour levels; dimension at least NC.
c     NC     (input) : +/- number of contour levels (less than or equal
c                      to dimension of C). If NC is positive, it is the
c                      number of contour levels, and the line-style is
c                      chosen automatically as described above. If NC is
c                      negative, it is minus the number of contour
c                      levels, and the current setting of line-style is
c                      used for all the contours.
c     TR     (input) : array defining a transformation between the I,J
c                      grid of the array and the world coordinates.
c                      The world coordinates of the array point A(I,J)
c                      are given by:
c                      X = TR(1) + TR(2)*I + TR(3)*J
c                      Y = TR(4) + TR(5)*I + TR(6)*J
c                      Usually TR(3) and TR(5) are zero - unless the
c                      coordinate transformation involves a rotation or
c                      shear.
c
c     SUBROUTINE PGCONT (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, TR)
c     INTEGER IDIM, JDIM, I1, J1, I2, J2, NC
c     REAL A(IDIM,JDIM), C(*), TR(6)
cSm030515 these lines work for g77 fortran under linux on PC 
      fmin=1.e+30
      fmax=-1.e+30
cSm030515 these lines work for digital fortan  under windows on PC
c      fmin=1.e+38
c      fmax=-1.e+38
      do i=1,nxmax
         do j=1,nymax
            fmin=MIN(f(i,j),fmin)
            fmax=MAX(f(i,j),fmax)
         enddo
      enddo
          
c      write(*,*)'contour.f contour2d fmin,fmax',fmin,fmax      
      do i=1,n_contour
         contour(i)=fmin+(i-1)*(fmax-fmin)/real(n_contour-1)
c         contour(i)=1.e-2*fmin+
c     .   (i-1)*1.e-2*(fmax-fmin)/real(n_contour-1)
      enddo

      do i=1,n_contour-1
         if ((contour(i).le.0.).and.(contour(i+1).gt.0.)) then
            contour(i)=0.
         endif
         write(*,*)'contour2d i,contour(i)', i,contour(i)
      enddo
      
 10   continue

      tr(1)=xmin-(xmax-xmin)/real(nxmax-1)
      tr(2)=(xmax-xmin)/real(nxmax-1)
      tr(3)=0.
      tr(4)=ymin-(ymax-ymin)/real(nymax-1)
      tr(5)=0.
      tr(6)=(ymax-ymin)/real(nymax-1)
      
      call PGCONT(f,nxmax,nymax,1,nxmax,1,nymax,contour,n_contour,tr)

c---------------------------------------------------------------
c     SUBROUTINE PGCONL (A, IDIM, JDIM, I1, I2, J1, J2, C, TR,
c    1                   LABEL, INTVAL, MININT)
c     INTEGER IDIM, JDIM, I1, J1, I2, J2, INTVAL, MININT
c     REAL A(IDIM,JDIM), C, TR(6)
c     CHARACTER*(*) LABEL
c     Label a contour map drawn with routine PGCONT. Routine PGCONT should
c     be called first to draw the contour lines, then this routine should be
c     called to add the labels. Labels are written at intervals along the
c     contour lines, centered on the contour lines with lettering aligned
c     in the up-hill direction. Labels are opaque, so a part of the under-
c     lying contour line is obscured by the label. Labels use the current
c     attributes (character height, line width, color index, character
c     font).

c     The first 9 arguments are the same as those supplied to PGCONT, and
c     should normally be identical to those used with PGCONT. Note that
c     only one contour level can be specified; tolabel more contours, call
c     PGCONL for each level.

c     The Label is supplied as a character string in argument LABEL.

c     The spacing of labels along the contour is specified by parameters
c     INTVAL and MININT. The routine follows the contour through the
c     array, counting the number of cells that the contour crosses. The
c     first label will be written in the MININT'th cell, and additional
c     labels will be written every INTVAL cells thereafter. A contour
c     that crosses less than MININT cells will not be labelled. Some
c     experimentation may be needed to get satisfactory results; a good
c     place to start is INTVAL=20, MININT=10.

c     Arguments:
c     A      (input) : data array.
c     IDIM   (input) : first dimension of A.
c     JDIM   (input) : second dimension of A.
c     I1, I2 (input) : range of first index to be contoured (inclusive).
c     J1, J2 (input) : range of second index to be contoyred 
c     C      (input) : the level of the contour to be labelled (one of the
c                      values given to PGCONT).
c      TR    (input) : array defining a transformation between the I,J
c                      grid of the array and the world coordinates.
c                      The world coordinates of the array point A(I,J)
c                      are given by:
c                      X = TR(1) + TR(2)*I + TR(3)*J
c                      Y = TR(4) + TR(5)*I + TR(6)*J
c                      Usually TR(3) and TR(5) are zero - unless the
c                      coordinate transformation involves a rotation or
c                      shear.
c   LABEL  (input) : character strings to be used to label the specified
c                    contour. Leading and trailing blank spaces are
c                    ignored.
c   INTVAL (input) : spacing along the contour between labels, in
c                    grid cells.
c   MININT (input) : contours that cross less than MININT cells
c                    will not be labelled.
      INTVAL=100
      MININT=40

      INTVAL=20
      MININT=10
c      do i=1,n_contour,2
       do i=1,n_contour,1

c         write(LABEL,'(1pe7.1)')contour(i)
         write(LABEL,'(i2)')i
c         WRITE(*,*)'i,contour(i),LABEL',i,contour(i),LABEL
         call  PGCONL(f,nxmax,nymax,1,nxmax,1,nymax,contour(i),
     .                tr,LABEL,INTVAL,MININT)
      enddo

      write(text,560)plot_name
 560  format("Contour values:",A)
      RILIN=8.
      CALL PGMTXT('B',RILIN,-0.2,0.0,text)

 570  format(1X,A,1pe10.3)
      do i=1,n_param
        rilin=rilin+1
        write(text,570)name_param(i),param(i)
        CALL PGMTXT('B',RILIN,-0.2,0.0,text)    
      enddo

      j=0
c      do i=1,n_contour,2
       do i=1,n_contour,1
         j=j+1
         write(text,150)i,contour(i)
         CALL PGMTXT('B',rilin+j,0.,0.,text)
      enddo

 150  format('i=',i3,'contour(I)=',1pe8.1)
      return
      end 

      real function rbound(r8)
c
c     Converts a real*8 argument to a real number,
c     equal to 0. (if r8=0.) or,
c     bounded in absolute value by 1.e-33 and 1.e+33.
c     This can be used to convert real*8 numbers to
c     real numbers, and to keep the resulting numbers
c     within the specified bounds.  This is necessary
c     for the PGPLOT library running on a 32-bit machine.
c     (1.e-35 was found to be too small in some cases,
c      on the DEC Alpha).
c     For a 64-bit machine, one might consider appropriate
c     adjustment of em33/ep33.
c
      real*8 r8,r8sign,r8abs
      real*8 em33,ep33,zero,one
      data em33/1.d-33/, ep33/1.d+33/, zero/0.d0/, one/1.d0/

      r8abs=abs(r8)
      if (r8abs.ne.zero) then
         r8sign=sign(one,r8)
         r8abs=min(r8abs,ep33)
         rbound=r8sign*max(r8abs,em33)
      else
         rbound=0.
      endif

      return
      end

      subroutine aminmx(x,n1,n2,n,xmin,xmax)
      double precision x,xmin,xmax
      dimension x(n2)
      xmin=+1.d100
      xmax=-1.d100
      do 1  i=n1,n2,n
      xmin=dmin1(xmin,x(i))
 1    xmax=dmax1(xmax,x(i))
      return
      end


      subroutine rminmx(x,n1,n2,n,xmin,xmax)
      real x,xmin,xmax
      dimension x(n2)
cSm030515 these lines work for g77 under linux at PC
      xmin=+1.e30
      xmax=-1.e30
cSm030515 these lines work for digital fortran under windows at PC
c      xmin=+1.e38
c      xmax=-1.e38
      do 1  i=n1,n2,n
      xmin=amin1(xmin,x(i))
 1    xmax=amax1(xmax,x(i))
      return
      end


c====================================================================
c====================================================================
      subroutine plotinit
c-----initialise the otput file plot.ps for pgplot
      character*100 line,line_
      integer PGOPEN     

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
         WRITE(*,*)'plotinit: TRYING TO PGOPEN at rank=',myrank
         STOP
      endif     
c-----Open a graphics device for PGPLOT output
c     INTEGER FUNCTION PGOPEN (DEVICE)
c     CHARACTER*(*) DEVICE
c     Here 
c     the name of otput file is plot.ps,
c     the device has /VPS tipe

c      write(*,*) 'PLOTINIT before PGOPEN'
      ier=PGOPEN('plot.ps/VCPS')   
c      ier=PGOPEN('plot.ps/VPS')   
c      ier=PGOPEN('plot.ps/PS') 
c      ier=PGOPEN('?')
      write(*,*)'in plotinit after PGOPEN: ier=',ier
      if (ier .LE. 0 ) then
         WRITE(*,*)'ier should be positive but ier.le.0 ier=',ier
         STOP
      endif
 
c     PGSCI -- set color index
      CALL PGSCI(1)
      write(*,*)'in plotinit/PGOPEN()  ier=1 is OK: ier=',ier
           
 1000  format(a80)
 1001  format(a80,'$')
      
      return
      end

c=====================================================================
c=====================================================================
      subroutine plotend
c-----close the work with pgplot    

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
         WRITE(*,*)'plotend/PGCLOS: TRYING TO close at rank=',myrank
         STOP
      endif     
   
c     Close the currently selected graphics device.
      CALL PGCLOS
      write(*,*)'plotend after PGCLOS'

      close(unit=100)
      return
      end

c=====================================================================
c=====================================================================
      subroutine plot1dt(x,y,nxa,nya,nx,ny,ll,ymin,ymax,ilabx,ilaby)
ctest
cmnt   This routine gives one plot of j=1,ny curves specified by
cmnt   (y(i[,j]),i=1,nx[(j)]), versus (x(i[,j]),i=1,nx[(j)]).
cmnt   x and/or y may be singly or doubly dimensioned.
cmnt   If x is doubly dimensioned, then set nxa = first
cmnt   dimension of x;  else set nxa = 0.
cmnt   Similarly for y and nya.
cmnt   If either y or x are singly dimensioned then the
cmnt   first dimension of all the vectors is nx(1).
cmnt   ll specifies data transformation:
cmnt      'linlin$' linear in x, linear in y
cmnt      'linlog$' linear in x, log in y  ,......
cmnt      'loglin$',  and $loglog$'
cmnt   ymin,ymax if .ne.0., specify min and max of yaxis.
c
      implicit double precision (a-h,o-z)
      dimension x(*),y(*),nx(*)
      character*(*) ilabx,ilaby
      character *4 text
      character(len=*) :: ll

c     Conversion to real function for PGPLOT
      REAL RBOUND
c     PGPLOT REAL Variables:
c      parameter (maxref=1000)
      parameter (maxref=10000)
      REAL RILIN
      REAL RPG1,RPG2
      REAL RX(maxref),RY(maxref)

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
         WRITE(*,*)'plot1d: TRYING TO PGPAGE at rank=',myrank
         STOP
      endif     

c      write(*,*)'begin of plot1dt'
c      write(*,*)'ilabx',ilabx
c      write(*,*)'ilaby',ilaby
c      write(*,*)'nx(1)',nx(1)
c      write(*,*)'x',(x(i),i=1,nx(1))
c      write(*,*)'nya',nya

      if(nya.eq.0) then
          call aminmx(y,1,nx(1),1,fminy,fmaxy)
          if (fminy .eq. fmaxy) fminy=.9*fmaxy-1.e-20
      else
          fminy=1.d100
          fmaxy=-1.d100
          do 100  n=1,ny
             in=1+(n-1)*nya
             nnx=nx(1)
             if(nxa*nya.ne.0) nnx=nx(n)
             call aminmx(y(in),1,nnx,1,fmin1,fmax1)
             fminy=dmin1(fminy,fmin1)
             fmaxy=dmax1(fmaxy,fmax1)
  100     continue
          if (fminy .eq. fmaxy) fminy=.9*fmaxy-1.e-20
      endif

c      write(*,*)'fminy,fmaxy',fminy,fmaxy

      if(nxa.eq.0) then
          call aminmx(x,1,nx(1),1,fminx,fmaxx)
          if (fminx .eq. fmaxx) fminx=.9*fmaxx-1.e-20
      else
          fminx=1.d100
          fmaxx=-1.d100
          do 101  n=1,ny
          in=1+(n-1)*nxa
          nnx=nx(1)
          if(nxa*nya.ne.0)  nnx=nx(n)
          call aminmx(x(in),1,nnx,1,fmin1,fmax1)
          fminx=dmin1(fminx,fmin1)
          fmaxx=dmax1(fmaxx,fmax1)
  101     continue
          if (fminx .eq. fmaxx) fminx=.9*fmaxx-1.e-20
      endif
      fmaxy_total=fmaxy
      fminy_total=fminy
      if(ymax.ne.0.d0)  fmaxy=ymax
      if(ymin.ne.0.d0)  fminy=ymin

c      write(*,*)' fmaxy,fmin', fmaxy,fmin
c      write(*,*)'ny',ny 
   
      do 200 n=1,ny
         inx=1+(n-1)*nxa
         iny=1+(n-1)*nya
         nnx=nx(1)
         if(nxa*nya.ne.0)  nnx=nx(n)
         nl=n
         if(n.gt.10)  nl=10
        
c         write(*,*)'nnx',nnx,'inx',inx
         DO I=1,NNX
            RX(I)=RBOUND(X(INX-1+I))
c            write(*,*)'i,rx(i),inx-1+i,x(inx-1+i)',
c     .      i,rx(i),inx-1+i,x(inx-1+i)
            RY(I)=RBOUND(Y(INY-1+I))
c            write(*,*)'i,ry(i),y(i)',i,ry(i),y(i)
         ENDDO
c         write(*,*)'rx',(RX(I),i=1,nnx)
         if (n.eq.1) then
           
             CALL PGPAGE
             CALL PGSLS(1) !
             CALL PGSVP(.2,.8,.3,.95)
c              write(*,*) 'PLOT1D: fminx,fmaxx,fminy,fmaxy',
c     +                    fminx,fmaxx,fminy,fmaxy
             CALL PGSWIN(rbound(fminx),rbound(fmaxx),
     +            rbound(fminy),rbound(fmaxy))
             CALL PGBOX('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
             CALL PGLAB(ilabx,ilaby, ' ' )          
             CALL PGLINE(nnx,RX(inx),RY(inx))
         endif
         if (n.eq.3) then
c------------vper1
             CALL PGPAGE
             CALL PGSVP(.2,.8,.3,.95)
c             write(*,*) 'PLOT1D: fminx,fmaxx,fminy,fmaxy',
c     +                    fminx,fmaxx,fminy,fmaxy
             CALL PGSWIN(rbound(fminx),rbound(fmaxx),
     +            rbound(fminy),rbound(fmaxy))
             CALL PGBOX('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
             CALL PGLAB(ilabx,'vper1', ' ' )
c             write(*,*)'plot1dt vper1 RY',RY
c             write(*,*)'plot1dt vper1 RX',RX 
             CALL PGLINE(nnx,RX(inx),RY(inx))
         endif !3

         if (n.eq.4) then
c------------zeta
             CALL PGPAGE
             CALL PGSVP(.2,.8,.3,.95)
c             write(*,*) 'PLOT1D: fminx,fmaxx,fminy,fmaxy',
c     +                    fminx,fmaxx,fminy,fmaxy
             CALL PGSWIN(rbound(fminx),rbound(fmaxx),
     +            rbound(fminy_total),rbound(fmaxy_total))
             CALL PGBOX('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
             CALL PGLAB(ilabx, 'zeta', ' ' )
             CALL PGLINE(nnx,RX(inx),RY(inx))
         endif !3
  200 continue
c      write(*,*) 'PLOT1Dt: HERE)'

      return
      end



      subroutine plot1dt_param(x,y,nxa,nya,nx,ny,ll,ymin,ymax,
     &plot_name,ilabx,ilaby,n_param,name_param,param)
ctest
cmnt   This routine gives one plot of j=1,ny curves specified by
cmnt   (y(i[,j]),i=1,nx[(j)]), versus (x(i[,j]),i=1,nx[(j)]).
cmnt   x and/or y may be singly or doubly dimensioned.
cmnt   If x is doubly dimensioned, then set nxa = first
cmnt   dimension of x;  else set nxa = 0.
cmnt   Similarly for y and nya.
cmnt   If either y or x are singly dimensioned then the
cmnt   first dimension of all the vectors is nx(1).
cmnt   ll specifies data transformation:
cmnt      'linlin$' linear in x, linear in y
cmnt      'linlog$' linear in x, log in y  ,......
cmnt      'loglin$',  and $loglog$'
cmnt   ymin,ymax if .ne.0., specify min and max of yaxis.
c      plot_name               is the name of the top of plot
c      n_param                 is the number of parameters
c      name_param(n_param)     are the names of parameters
c      param(n_param)          are the values of parameters
c--------------------------------------------------------------
      implicit double precision (a-h,o-z)
cSAP80709
c      dimension x(1),y(1),nx(1)
      dimension x(*),y(*),nx(*)
      character*(*)plot_name,ilabx,ilaby
      integer  n_param                 !is the number of parameters
      character*(*)name_param(n_param) !are the names of parameters
      real param(*)                    !are the values of parameters

      character*72 text_param             
      character *4 text
      character(len=*) :: ll

c     Conversion to real function for PGPLOT
      REAL RBOUND
c     PGPLOT REAL Variables:
c      parameter (maxref=1000)
      parameter (maxref=10000)
      REAL RILIN
      REAL RPG1,RPG2
      REAL RX(maxref),RY(maxref)
               
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
         WRITE(*,*)'plot1dt_param: TRYING TO PGPAGE at rank=',myrank
         STOP
      endif     

c      write(*,*)'begin of plot1dt'
c      write(*,*)'ilabx ',ilabx
c      write(*,*)'ilaby ',ilaby
c      write(*,*)'nx(1)',nx(1)
c      write(*,*)'x',(x(i),i=1,nx(1))

      if(nya.eq.0) then
          call aminmx(y,1,nx(1),1,fminy,fmaxy)
          if (fminy .eq. fmaxy) fminy=.9*fmaxy-1.e-20
      else
          fminy=1.d100
          fmaxy=-1.d100
          do 100  n=1,ny
             in=1+(n-1)*nya
             nnx=nx(1)
             if(nxa*nya.ne.0) nnx=nx(n)
             call aminmx(y(in),1,nnx,1,fmin1,fmax1)
             fminy=dmin1(fminy,fmin1)
             fmaxy=dmax1(fmaxy,fmax1)
  100     continue
          if (fminy .eq. fmaxy) fminy=.9*fmaxy-1.e-20
      endif
      if(nxa.eq.0) then
          call aminmx(x,1,nx(1),1,fminx,fmaxx)
          if (fminx .eq. fmaxx) fminx=.9*fmaxx-1.e-20
      else
          fminx=1.d100
          fmaxx=-1.d100
          do 101  n=1,ny
          in=1+(n-1)*nxa
          nnx=nx(1)
          if(nxa*nya.ne.0)  nnx=nx(n)
          call aminmx(x(in),1,nnx,1,fmin1,fmax1)
          fminx=dmin1(fminx,fmin1)
          fmaxx=dmax1(fmaxx,fmax1)
  101     continue
          if (fminx .eq. fmaxx) fminx=.9*fmaxx-1.e-20
      endif
      fmaxy_total=fmaxy
      fminy_total=fminy
      if(ymax.ne.0.d0)  fmaxy=ymax
      if(ymin.ne.0.d0)  fminy=ymin

c      write(*,*)'ny',ny 
      CALL PGSLS(1) ! set usial full line 
      CALL PGSCI(1) ! set black line
      do 200 n=1,ny
         inx=1+(n-1)*nxa
         iny=1+(n-1)*nya
         nnx=nx(1)
         if(nxa*nya.ne.0)  nnx=nx(n)
         nl=n
         if(n.gt.10)  nl=10
        
c         write(*,*)'nnx',nnx,'inx',inx
         DO I=1,NNX
            RX(I)=RBOUND(X(INX-1+I))
c            write(*,*)'i,rx(i),inx-1+i,x(inx-1+i)',
c     .      i,rx(i),inx-1+i,x(inx-1+i)
            RY(I)=RBOUND(Y(INY-1+I))
c            write(*,*)'i,ry(i),y(i)',i,ry(i),y(i)
         ENDDO
c         write(*,*)'rx',(RX(I),i=1,nnx)
         if (n.eq.1) then
             CALL PGPAGE
             CALL PGSVP(.2,.8,.3,0.80) !070710 new
c             CALL PGSVP(.2,.8,.3,.95) !070710 old
c             CALL PGSVP(.2,.8,.5,.95)

c              write(*,*) 'PLOT1D: fminx,fmaxx,fminy,fmaxy',
c     +                    fminx,fmaxx,fminy,fmaxy
             CALL PGSWIN(rbound(fminx),rbound(fmaxx),
     +            rbound(fminy),rbound(fmaxy))
             CALL PGBOX('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
c             CALL PGLAB(ilabx,ilaby, ' ' ) 
             CALL PGLAB(ilabx,ilaby,plot_name)                    
             CALL PGLINE(nnx,RX(inx),RY(inx))
         endif
         if (n.eq.3) then
c------------vper1
             CALL PGPAGE
             CALL PGSVP(.2,.8,.3,.95)
c             CALL PGSVP(.2,.8,.5,.95)
c             write(*,*) 'PLOT1D: fminx,fmaxx,fminy,fmaxy',
c     +                    fminx,fmaxx,fminy,fmaxy
             CALL PGSWIN(rbound(fminx),rbound(fmaxx),
     +            rbound(fminy),rbound(fmaxy))
             CALL PGBOX('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
             CALL PGLAB(ilabx,'vper1', ' ' )
c             write(*,*)'plot1dt vper1 RY',RY
c             write(*,*)'plot1dt vper1 RX',RX 
             CALL PGLINE(nnx,RX(inx),RY(inx))
         endif !3

         if (n.eq.4) then
c------------zeta
             CALL PGPAGE
c             CALL PGSVP(.2,.8,.3,.95)
             CALL PGSVP(.2,.8,.5,.95)
c             write(*,*) 'PLOT1D: fminx,fmaxx,fminy,fmaxy',
c     +                    fminx,fmaxx,fminy,fmaxy
             CALL PGSWIN(rbound(fminx),rbound(fmaxx),
     +            rbound(fminy_total),rbound(fmaxy_total))
             CALL PGBOX('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
             CALL PGLAB(ilabx, 'zeta', ' ' )
             CALL PGLINE(nnx,RX(inx),RY(inx))
         endif !3
  200 continue

c      write(*,*)'plot1dt_param n_param=',n_param
c      write(*,*)'plotqdt_param name_param',(name_param(i),i=1,n_param)
c      write(*,*)'plotqdt_param param',(param(i),i=1,n_param)
c      write(*,*)'rilin',rilin

 570  format(1X,A,1pe10.3)
      rilin=8.
      do i=1,n_param
        rilin=rilin+1
c        write(*,*)'plot1dt_param i=',i,'rilin=',rilin
        write(text_param,570)name_param(i),param(i)
cc       CALL PGMTXT('B',RILIN,-0.2,0.0,text_param)         
        CALL PGMTXT('B',RILIN,-0.1,0.0,text_param)    
      enddo

c      write(*,*) 'PLOT1Dt: HERE)'

      return
      end

      subroutine contour2d_S(nxmax_a,nymax_a,n_contour_a,
     .nxmax,nymax,f,x,y,plot_name,x_name,y_name,
     .n_contour,contour,name_param,param,n_param)
c-----print contours of 2D array f(1:nxmax,1:nymax)
c
c     input
c      fmap(1:nxmax_a,1:nymax_a) 2D array of the function f(x,y)
c      nxmax_a,nymax_a are maximal dimensions of f,x,y
c      x(1:nxmax_a)  array of x variable
c      y(1:nymax_a)  array of y variable
c      nxmax,  range of the first and second index to be contourd (inclusive)
c      nymax,   range of are the boindaitsarray of x variable
c      plot_name  the name of the top of plot
c      x_name     the name of x axis  
c      y_name     the name of y axis
c      n_contour  the number of contours
c      n_contour_a  the maximal  number of contours
c      contour(n_contour_a)      work array for the contour's values
c      n_param is the number of parameters
c      name_param(n_param)     are the names of parameters
c      param(n_param)          are the values of parameters

      implicit none
c-----input
      integer nxmax,nymax,n_contour,n_param,
     &nxmax_a,nymax_a,n_contour_a
      real f(nxmax_a,nymax_a),x(nxmax_a),y(nymax_a),contour(n_contour_a)
      character*(*)plot_name,x_name,y_name
      character*(*)name_param(n_param)
      real param(*)
      real y_e,x_e,cnpar,T
c-----external
      real rbound !convert double to real bounded 1.e-33 and 1.e+33
c-----local
      real tr(6),RILIN
      real xmin,xmax
      real ymin,ymax
      real fmin,fmax
      integer i,j
      
c      CHARACTER*7 LABEL
      CHARACTER*2 LABEL

      integer INTVAL,MININT
      character*72 text          

               
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
         WRITE(*,*)'contour2d_S: TRYING TO PGPAGE at rank=',myrank
         STOP
      endif     

      if(nxmax.gt.nxmax_a) then
        WRITE(*,*)'in contour2d_S nxmax >nxmax_a'
        WRITE(*,*)'It should be nxmax.le.nxmax_a'
        WRITE(*,*)' nxmax,nxmax_a',nxmax,nxmax_a
        WRITE(*,*)'Please change arguments'
        STOP 'contour2d_S'
      endif

      if(nymax.gt.nymax_a) then
        WRITE(*,*)'in contour2d_S nymax > nymax_a'
        WRITE(*,*)'It should be nymax.le.nymax_a'
        WRITE(*,*)'nymax,nymax_a',nymax,nymax_a
        WRITE(*,*)'Please change arguments'
        WRITE(*,*)'nymax,nymax_a',nymax,nymax_a
        STOP 'contour2d_S'
      endif

      if(n_contour.gt.n_contour_a) then
        WRITE(*,*)'in contour2d_S n_contour > n_contour_a'
        WRITE(*,*)'It should be n_contour.le.n_contour_a'
        WRITE(*,*)'n_contour,n_contour_a',n_contour,n_contour_a
        WRITE(*,*)'Please change arguments'
        STOP 'contour2d_S'
      endif

c-----Advance plotter to a new page or panel, clearing the screen if
c     necessary.
      CALL PGPAGE  

      CALL PGSLW(4)
c-----Change the size and position of the viewport, specifying
c     the viewport in normalized device coordinates.  Normalized
c     device coordinates run from 0 to 1 in each dimension.
c     SUBROUTINE PGSVP (XLEFT, XRIGHT, YBOT, YTOP)
c     REAL XLEFT, XRIGHT, YBOT, YTOP

c      CALL PGSVP(.2,.8,.3,.95) 
c      CALL PGSVP(.2,.8,.5,.95) !070711 old 
      CALL PGSVP(.2,.8,.5,.80) !070711 new

c-----Change the window in world coordinate space that is to be mapped on
c     to the viewport.  Usually PGSWIN is called automatically by PGENV,
c     but it may be called directly by the user.
c     SUBROUTINE PGSWIN (X1, X2, Y1, Y2)
c     REAL X1, X2, Y1, Y2
      call rminmx(x,1,nxmax,1,xmin,xmax)
      if (xmin .eq. xmax) xmin=.9*xmax-1.e-20
      call rminmx(y,1,nymax,1,ymin,ymax)
      if (ymin .eq. ymax) ymin=.9*xmax-1.e-20
      CALL PGSWIN (xmin,xmax,ymin,ymax)

c-----PGBOX -- draw labeled frame around viewport
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)

c-----PGLAB write labels  for x,y axis and top of plot
      CALL PGLAB(x_name,y_name,plot_name)
 
c-----draw a contour map of an array.  The map is truncated if
c     necessary at the boundaries of the viewport.  Each contour line
c     is drawn with the current line attributes (color index, style, and
c     width); except that if argument NC is positive (see below), the line
c     style is set by PGCONT to 1 (solid) for positive contours or 2
c     (dashed) for negative contours.
c
c     Arguments:
c     A      (input) : data array.
c     IDIM   (input) : first dimension of A.
c     JDIM   (input) : second dimension of A.
c     I1, I2 (input) : range of first index to be contoured (inclusive).
c     J1, J2 (input) : range of second index to be contoured (inclusive).
c     C      (input) : array of NC contour levels; dimension at least NC.
c     NC     (input) : +/- number of contour levels (less than or equal
c                      to dimension of C). If NC is positive, it is the
c                      number of contour levels, and the line-style is
c                      chosen automatically as described above. If NC is
c                      negative, it is minus the number of contour
c                      levels, and the current setting of line-style is
c                      used for all the contours.
c     TR     (input) : array defining a transformation between the I,J
c                      grid of the array and the world coordinates.
c                      The world coordinates of the array point A(I,J)
c                      are given by:
c                      X = TR(1) + TR(2)*I + TR(3)*J
c                      Y = TR(4) + TR(5)*I + TR(6)*J
c                      Usually TR(3) and TR(5) are zero - unless the
c                      coordinate transformation involves a rotation or
c                      shear.
c
c     SUBROUTINE PGCONT (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, TR)
c     INTEGER IDIM, JDIM, I1, J1, I2, J2, NC
c     REAL A(IDIM,JDIM), C(*), TR(6)
cSm030515 these lines work for g77 fortran under linux on PC 
      fmin=1.e+30
      fmax=-1.e+30
cSm030515 these lines work for digital fortan  under windows on PC
c      fmin=1.e+38
c      fmax=-1.e+38
      do i=1,nxmax
         do j=1,nymax
            fmin=MIN(f(i,j),fmin)
            fmax=MAX(f(i,j),fmax)
         enddo
      enddo
          
c      write(*,*)'contour.f conrour2d fmin,fmax',fmin,fmax      
      do i=1,n_contour
         contour(i)=fmin+(i-1)*(fmax-fmin)/real(n_contour-1)
c         contour(i)=1.e-2*fmin+
c     .   (i-1)*1.e-2*(fmax-fmin)/real(n_contour-1)
      enddo

      do i=1,n_contour-1
         if ((contour(i).le.0.).and.(contour(i+1).gt.0.)) then
            contour(i)=0.
         endif
      enddo
      
 10   continue

      tr(1)=xmin-(xmax-xmin)/real(nxmax-1)
      tr(2)=(xmax-xmin)/real(nxmax-1)
      tr(3)=0.
      tr(4)=ymin-(ymax-ymin)/real(nymax-1)
      tr(5)=0.
      tr(6)=(ymax-ymin)/real(nymax-1)
       
      call PGCONT(f,nxmax_a,nymax_a,1,nxmax,1,nymax,contour,n_contour,
     &tr)

c---------------------------------------------------------------
c     SUBROUTINE PGCONL (A, IDIM, JDIM, I1, I2, J1, J2, C, TR,
c    1                   LABEL, INTVAL, MININT)
c     INTEGER IDIM, JDIM, I1, J1, I2, J2, INTVAL, MININT
c     REAL A(IDIM,JDIM), C, TR(6)
c     CHARACTER*(*) LABEL
c     Label a contour map drawn with routine PGCONT. Routine PGCONT should
c     be called first to draw the contour lines, then this routine should be
c     called to add the labels. Labels are written at intervals along the
c     contour lines, centered on the contour lines with lettering aligned
c     in the up-hill direction. Labels are opaque, so a part of the under-
c     lying contour line is obscured by the label. Labels use the current
c     attributes (character height, line width, color index, character
c     font).

c     The first 9 arguments are the same as those supplied to PGCONT, and
c     should normally be identical to those used with PGCONT. Note that
c     only one contour level can be specified; tolabel more contours, call
c     PGCONL for each level.

c     The Label is supplied as a character string in argument LABEL.

c     The spacing of labels along the contour is specified by parameters
c     INTVAL and MININT. The routine follows the contour through the
c     array, counting the number of cells that the contour crosses. The
c     first label will be written in the MININT'th cell, and additional
c     labels will be written every INTVAL cells thereafter. A contour
c     that crosses less than MININT cells will not be labelled. Some
c     experimentation may be needed to get satisfactory results; a good
c     place to start is INTVAL=20, MININT=10.

c     Arguments:
c     A      (input) : data array.
c     IDIM   (input) : first dimension of A.
c     JDIM   (input) : second dimension of A.
c     I1, I2 (input) : range of first index to be contoured (inclusive).
c     J1, J2 (input) : range of second index to be contoyred 
c     C      (input) : the level of the contour to be labelled (one of the
c                      values given to PGCONT).
c      TR    (input) : array defining a transformation between the I,J
c                      grid of the array and the world coordinates.
c                      The world coordinates of the array point A(I,J)
c                      are given by:
c                      X = TR(1) + TR(2)*I + TR(3)*J
c                      Y = TR(4) + TR(5)*I + TR(6)*J
c                      Usually TR(3) and TR(5) are zero - unless the
c                      coordinate transformation involves a rotation or
c                      shear.
c   LABEL  (input) : character strings to be used to label the specified
c                    contour. Leading and trailing blank spaces are
c                    ignored.
c   INTVAL (input) : spacing along the contour between labels, in
c                    grid cells.
c   MININT (input) : contours that cross less than MININT cells
c                    will not be labelled.
      INTVAL=100
      MININT=40

      INTVAL=20
      MININT=10
c      do i=1,n_contour,2
       do i=1,n_contour,1

c         write(LABEL,'(1pe7.1)')contour(i)
         write(LABEL,'(i2)')i
c         WRITE(*,*)'i,contour(i),LABEL',i,contour(i),LABEL
         call  PGCONL(f,nxmax_a,nymax_a,1,nxmax,1,nymax,contour(i),
     .                tr,LABEL,INTVAL,MININT)
      enddo

      write(text,560)plot_name
 560  format("Contour values:",A)
c      RILIN=8.
      RILIN=5.
      CALL PGMTXT('B',RILIN,-0.2,0.0,text)

 570  format(1X,A,1pe10.3)
      do i=1,n_param
        rilin=rilin+1
        write(text,570)name_param(i),param(i)
        CALL PGMTXT('B',RILIN,-0.2,0.0,text)    
      enddo

      j=0
c      do i=1,n_contour,2
c       do i=1,n_contour,1
c         j=j+1
c         write(text,150)i,contour(i)
c         CALL PGMTXT('B',rilin+j,0.,0.,text)
c      enddo

c 150  format('i=',i3,'contour(I)=',1pe8.1)

      do i=1,n_contour/2
         j=j+1
         write(text,160)2*i-1,contour(2*i-1),2*i,contour(2*i)
c         CALL PGMTXT('B',rilin+j,0.,0.,text)   !070711 old
         CALL PGMTXT('B',rilin+j,-0.2,0.,text) !070711 new
      enddo

 160  format('i=',i2,'contour(I)=',1pe8.1,'i=',i2,'contour(I)=',1pe8.1)

      return
      end 


c!!!!!!!!!!!
      subroutine plot1dt_marker_param(x,y,nxa,nya,nx,ny,
     *n_marker,x_marker,y_marker,n_symbol_marker,    
     &ll,ymin,ymax,
     &plot_name,ilabx,ilaby,n_param,name_param,param)
ctest
cmnt   This routine gives one plot of j=1,ny curves specified by
cmnt   (y(i[,j]),i=1,nx[(j)]), versus (x(i[,j]),i=1,nx[(j)]).
cmnt   x and/or y may be singly or doubly dimensioned.
cmnt   If x is doubly dimensioned, then set nxa = first
cmnt   dimension of x;  else set nxa = 0.
cmnt   Similarly for y and nya.
cmnt   If either y or x are singly dimensioned then the
cmnt   first dimension of all the vectors is nx(1).
cmnt   ll specifies data transformation:
cmnt      'linlin$' linear in x, linear in y
cmnt      'linlog$' linear in x, log in y  ,......
cmnt      'loglin$',  and $loglog$'
cmnt   ymin,ymax if .ne.0., specify min and max of yaxis.
c      plot_name               is the name of the top of plot
c      n_param                 is the number of parameters
c      name_param(n_param)     are the names of parameters
c      param(n_param)          are the values of parameters
c
c      Drawing Markers:
c     *n_marker         number of points to be marked,
c      x_marker         x and y coordinates of the points,
c      y_marker
c      n_symbol_marker  The number of the symbol to be used to 
c                       mark the points. 
c                        =9 circle with a central dot 
c--------------------------------------------------------------
      implicit double precision (a-h,o-z)
cSAP080709
c      dimension x(1),y(1),nx(1)
      dimension x(*),y(*),nx(*)
      character*(*)plot_name,ilabx,ilaby
      integer  n_param                 !is the number of parameters
      character*(*)name_param(n_param) !are the names of parameters
      real param(*)                    !are the values of parameters

      character*72 text_param             
      character *4 text
      character(len=*) :: ll

c-----for Drawing Markers:
      integer n_marker        ! number of points to be marked,
      real  x_marker(*),y_marker(*)
      integer n_symbol_marker

      

c     Conversion to real function for PGPLOT
      REAL RBOUND
c     PGPLOT REAL Variables:
c      parameter (maxref=1000)
      parameter (maxref=10000)
      REAL RILIN
      REAL RPG1,RPG2
      REAL RX(maxref),RY(maxref)
               
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
      WRITE(*,*)'plot1dt_marker_param:TRYING TO PGPAGE at rank',myrank
         STOP
      endif     

c      write(*,*)'begin of plot1dt plot1dt_marker_param'
c      write(*,*)'ilabx',ilabx
c      write(*,*)'ilaby',ilaby
c      write(*,*)'nx(1)',nx(1)
c      write(*,*)'x',(x(i),i=1,nx(1))
      if(nya.eq.0) then
          call aminmx(y,1,nx(1),1,fminy,fmaxy)
          if (fminy .eq. fmaxy) fminy=.9*fmaxy-1.e-20
      else
          fminy=1.d100
          fmaxy=-1.d100
          do 100  n=1,ny
             in=1+(n-1)*nya
             nnx=nx(1)
             if(nxa*nya.ne.0) nnx=nx(n)
             call aminmx(y(in),1,nnx,1,fmin1,fmax1)
             fminy=dmin1(fminy,fmin1)
             fmaxy=dmax1(fmaxy,fmax1)
  100     continue
          if (fminy .eq. fmaxy) fminy=.9*fmaxy-1.e-20
      endif
      if(nxa.eq.0) then
          call aminmx(x,1,nx(1),1,fminx,fmaxx)
          if (fminx .eq. fmaxx) fminx=.9*fmaxx-1.e-20
      else
          fminx=1.d100
          fmaxx=-1.d100
          do 101  n=1,ny
          in=1+(n-1)*nxa
          nnx=nx(1)
          if(nxa*nya.ne.0)  nnx=nx(n)
          call aminmx(x(in),1,nnx,1,fmin1,fmax1)
          fminx=dmin1(fminx,fmin1)
          fmaxx=dmax1(fmaxx,fmax1)
  101     continue
          if (fminx .eq. fmaxx) fminx=.9*fmaxx-1.e-20
      endif
      fmaxy_total=fmaxy
      fminy_total=fminy
      if(ymax.ne.0.d0)  fmaxy=ymax
      if(ymin.ne.0.d0)  fminy=ymin

c      write(*,*)'ny',ny 
      do 200 n=1,ny
         inx=1+(n-1)*nxa
         iny=1+(n-1)*nya
         nnx=nx(1)
         if(nxa*nya.ne.0)  nnx=nx(n)
         nl=n
         if(n.gt.10)  nl=10
        
c         write(*,*)'nnx',nnx,'inx',inx
         DO I=1,NNX
c           write(*,*)'i,inx-1+i',i,inx-1+i
            RX(I)=RBOUND(X(INX-1+I))
c            write(*,*)'i,rx(i),inx-1+i,x(inx-1+i)',
c     .      i,rx(i),inx-1+i,x(inx-1+i)
            RY(I)=RBOUND(Y(INY-1+I))
c            write(*,*)'i,ry(i),y(i)',i,ry(i),y(i)
         ENDDO
c         write(*,*)'rx',(RX(I),i=1,nnx)
         if (n.eq.1) then
             CALL PGPAGE
c             CALL PGSVP(.2,.8,.3,.95)
c             CALL PGSVP(.2,.8,.5,.95)  !070710 old
             CALL PGSVP(.2,.8,.4,0.80) !070710 new

c              write(*,*) 'PLOT1D: fminx,fmaxx,fminy,fmaxy',
c     +                    fminx,fmaxx,fminy,fmaxy
             CALL PGSWIN(rbound(fminx),rbound(fmaxx),
     +            rbound(fminy),rbound(fmaxy))
             CALL PGBOX('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
c             CALL PGLAB(ilabx,ilaby, ' ' ) 
             CALL PGLAB(ilabx,ilaby,plot_name)                    
             CALL PGLINE(nnx,RX(inx),RY(inx))

c------------Drawing Markers
ctest
             write(*,*)'contour.f in  plot1dt_marker_param' 
             write(*,*)'n_marker',n_marker
             write(*,*)'n_symbol_marker',n_symbol_marker
             do i=1,n_marker
               write(*,*)'i,x_marker(i),y_marker(i) ',
     &                    i,x_marker(i),y_marker(i) 
             enddo
cendtest
    
             CALL PGPT(n_marker,x_marker,y_marker,n_symbol_marker)

         endif
         if (n.eq.3) then
c------------vper1
             CALL PGPAGE
c             CALL PGSVP(.2,.8,.3,.95)
             CALL PGSVP(.2,.8,.5,.95)
c             write(*,*) 'PLOT1D: fminx,fmaxx,fminy,fmaxy',
c     +                    fminx,fmaxx,fminy,fmaxy
             CALL PGSWIN(rbound(fminx),rbound(fmaxx),
     +            rbound(fminy),rbound(fmaxy))
             CALL PGBOX('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
             CALL PGLAB(ilabx,'vper1', ' ' )
c             write(*,*)'plot1dt vper1 RY',RY
c             write(*,*)'plot1dt vper1 RX',RX 
             CALL PGLINE(nnx,RX(inx),RY(inx))
         endif !3

         if (n.eq.4) then
c------------zeta
             CALL PGPAGE
c             CALL PGSVP(.2,.8,.3,.95)
             CALL PGSVP(.2,.8,.5,.95)
c             write(*,*) 'PLOT1D: fminx,fmaxx,fminy,fmaxy',
c     +                    fminx,fmaxx,fminy,fmaxy
             CALL PGSWIN(rbound(fminx),rbound(fmaxx),
     +            rbound(fminy_total),rbound(fmaxy_total))
             CALL PGBOX('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
             CALL PGLAB(ilabx, 'zeta', ' ' )
             CALL PGLINE(nnx,RX(inx),RY(inx))
         endif !3
  200 continue

      write(*,*)'plot1dt_param n_param=',n_param
      write(*,*)'plot1dt_param name_param',(name_param(i),i=1,n_param)
      write(*,*)'plot1dt_param param',(param(i),i=1,n_param)
c      write(*,*)'riline',riline

 570  format(1X,A,1pe11.3)
      rilin=8.
      do i=1,n_param
        rilin=rilin+1
c        write(*,*)'plot1dt_param i=',i,'rilin=',rilin
        write(text_param,570)name_param(i),param(i)
c        write(*,*)'i,text_param',i,text_param
c       CALL PGMTXT('B',RILIN,-0.2,0.0,text_param)         
        CALL PGMTXT('B',RILIN,-0.1,0.0,text_param)    
      enddo

c      write(*,*) 'PLOT1Dt: HERE)'

      return
      end


c=================
      subroutine contour2d_test(nxmax,nymax,nxmaxa,nymaxa,
     &f,x,y,plot_name,
     &x_name,y_name,
     &n_contour,contour,
     &name_param,param,n_param)
c-----print contours of 2D array f(1:nxmax,1:nymax)
c
c     input
c      fmap(1:nxmax,1:nymax) 2D array of the function f(x,y)
c      nxmax,nymax are the dimentions of f
c      nxmaxa,nymaxa are max of    nxmax,nyma
c      x(1:nxmax)  array of x variable
c      y(1:nymax)  array of y variable
c      plot_name  the name of the top of plot
c      x_name     the name of x axis  
c      y_name     the name of y axis
c      n_contour  the number of contours
c      contour(n_contour)      work array for the contour's values
c      n_param is the number of parameters
c      name_param(n_param)     are the names of parameters
c      param(n_param)          are the values of parameters

      implicit none
c-----input
      integer nxmax,nymax,nxmaxa,nymaxa,
     &n_contour,n_param
      real f(nxmaxa,nymaxa),x(nxmaxa),y(nymaxa),
     &contour(n_contour)
      character*(*)plot_name,x_name,y_name
      character*(*)name_param(n_param)
      real param(*)
c-----external
      real rbound !convert double to real bounded 1.e-33 and 1.e+33
c-----local
      real tr(6),RILIN
      real xmin,xmax
      real ymin,ymax
      real fmin,fmax
      integer i,j
      
c      CHARACTER*7 LABEL
      CHARACTER*2 LABEL

      integer INTVAL,MININT
      character*72 text       

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
         WRITE(*,*)'contour2d_test: TRYING TO PGPAGE at rank=',myrank
         STOP
      endif     
      
      write(*,*)' contour2d_test nxmax,nymax',nxmax,nymax
      do i=1,nxmax
         do j=1,nymax
            write(*,*)'i,j,f(i,j)',i,j,f(i,j)
         enddo
      enddo


c-----Advance plotter to a new page or panel, clearing the screen if
c     necessary.
      CALL PGPAGE  

c-----Change the size and position of the viewport, specifying
c     the viewport in normalized device coordinates.  Normalized
c     device coordinates run from 0 to 1 in each dimension.
c     SUBROUTINE PGSVP (XLEFT, XRIGHT, YBOT, YTOP)
c     REAL XLEFT, XRIGHT, YBOT, YTOP

c      CALL PGSVP(.2,.8,.3,.95) 
      CALL PGSVP(.2,.8,.5,.95) 

c-----Change the window in world coordinate space that is to be mapped on
c     to the viewport.  Usually PGSWIN is called automatically by PGENV,
c     but it may be called directly by the user.
c     SUBROUTINE PGSWIN (X1, X2, Y1, Y2)
c     REAL X1, X2, Y1, Y2
      call rminmx(x,1,nxmax,1,xmin,xmax)
      if (xmin .eq. xmax) xmin=.9*xmax-1.e-20
      call rminmx(y,1,nymax,1,ymin,ymax)
      if (ymin .eq. ymax) ymin=.9*xmax-1.e-20
      CALL PGSWIN (xmin,xmax,ymin,ymax)

c-----PGBOX -- draw labeled frame around viewport
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)

c-----PGLAB write labels  for x,y axis and top of plot
      CALL PGLAB(x_name,y_name,plot_name)
 
c-----draw a contour map of an array.  The map is truncated if
c     necessary at the boundaries of the viewport.  Each contour line
c     is drawn with the current line attributes (color index, style, and
c     width); except that if argument NC is positive (see below), the line
c     style is set by PGCONT to 1 (solid) for positive contours or 2
c     (dashed) for negative contours.
c
c     Arguments:
c     A      (input) : data array.
c     IDIM   (input) : first dimension of A.
c     JDIM   (input) : second dimension of A.
c     I1, I2 (input) : range of first index to be contoured (inclusive).
c     J1, J2 (input) : range of second index to be contoured (inclusive).
c     C      (input) : array of NC contour levels; dimension at least NC.
c     NC     (input) : +/- number of contour levels (less than or equal
c                      to dimension of C). If NC is positive, it is the
c                      number of contour levels, and the line-style is
c                      chosen automatically as described above. If NC is
c                      negative, it is minus the number of contour
c                      levels, and the current setting of line-style is
c                      used for all the contours.
c     TR     (input) : array defining a transformation between the I,J
c                      grid of the array and the world coordinates.
c                      The world coordinates of the array point A(I,J)
c                      are given by:
c                      X = TR(1) + TR(2)*I + TR(3)*J
c                      Y = TR(4) + TR(5)*I + TR(6)*J
c                      Usually TR(3) and TR(5) are zero - unless the
c                      coordinate transformation involves a rotation or
c                      shear.
c
c     SUBROUTINE PGCONT (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, TR)
c     INTEGER IDIM, JDIM, I1, J1, I2, J2, NC
c     REAL A(IDIM,JDIM), C(*), TR(6)
    
cSm030515 these lines work for g77 fortran under linux on PC 
      fmin=1.e+30
      fmax=-1.e+30
cSm030515 these lines work for digital fortan  under windows on PC
c      fmin=1.e+38
c      fmax=-1.e+38
      do i=1,nxmax
         do j=1,nymax
            fmin=MIN(f(i,j),fmin)
            fmax=MAX(f(i,j),fmax)
         enddo
      enddo
          
c      write(*,*)'contour.f contour2d fmin,fmax',fmin,fmax      
      do i=1,n_contour
         contour(i)=fmin+(i-1)*(fmax-fmin)/real(n_contour-1)
c         contour(i)=1.e-2*fmin+
c     .   (i-1)*1.e-2*(fmax-fmin)/real(n_contour-1)
      enddo

      do i=1,n_contour-1
         if ((contour(i).le.0.).and.(contour(i+1).gt.0.)) then
            contour(i)=0.
         endif
      enddo
      
 10   continue

      tr(1)=xmin-(xmax-xmin)/real(nxmax-1)
      tr(2)=(xmax-xmin)/real(nxmax-1)
      tr(3)=0.
      tr(4)=ymin-(ymax-ymin)/real(nymax-1)
      tr(5)=0.
      tr(6)=(ymax-ymin)/real(nymax-1)

      write(*,*)'n_contour',n_contour
      do i=1,n_contour
         write(*,*)'i,contour(i)',i,contour(i)
      enddo
      write(*,*)'tr',tr
  
      call PGCONT(f,nxmax,nymax,1,nxmax,1,nymax,contour,n_contour,tr)

      return
c---------------------------------------------------------------
c     SUBROUTINE PGCONL (A, IDIM, JDIM, I1, I2, J1, J2, C, TR,
c    1                   LABEL, INTVAL, MININT)
c     INTEGER IDIM, JDIM, I1, J1, I2, J2, INTVAL, MININT
c     REAL A(IDIM,JDIM), C, TR(6)
c     CHARACTER*(*) LABEL
c     Label a contour map drawn with routine PGCONT. Routine PGCONT should
c     be called first to draw the contour lines, then this routine should be
c     called to add the labels. Labels are written at intervals along the
c     contour lines, centered on the contour lines with lettering aligned
c     in the up-hill direction. Labels are opaque, so a part of the under-
c     lying contour line is obscured by the label. Labels use the current
c     attributes (character height, line width, color index, character
c     font).

c     The first 9 arguments are the same as those supplied to PGCONT, and
c     should normally be identical to those used with PGCONT. Note that
c     only one contour level can be specified; tolabel more contours, call
c     PGCONL for each level.

c     The Label is supplied as a character string in argument LABEL.

c     The spacing of labels along the contour is specified by parameters
c     INTVAL and MININT. The routine follows the contour through the
c     array, counting the number of cells that the contour crosses. The
c     first label will be written in the MININT'th cell, and additional
c     labels will be written every INTVAL cells thereafter. A contour
c     that crosses less than MININT cells will not be labelled. Some
c     experimentation may be needed to get satisfactory results; a good
c     place to start is INTVAL=20, MININT=10.

c     Arguments:
c     A      (input) : data array.
c     IDIM   (input) : first dimension of A.
c     JDIM   (input) : second dimension of A.
c     I1, I2 (input) : range of first index to be contoured (inclusive).
c     J1, J2 (input) : range of second index to be contoyred 
c     C      (input) : the level of the contour to be labelled (one of the
c                      values given to PGCONT).
c      TR    (input) : array defining a transformation between the I,J
c                      grid of the array and the world coordinates.
c                      The world coordinates of the array point A(I,J)
c                      are given by:
c                      X = TR(1) + TR(2)*I + TR(3)*J
c                      Y = TR(4) + TR(5)*I + TR(6)*J
c                      Usually TR(3) and TR(5) are zero - unless the
c                      coordinate transformation involves a rotation or
c                      shear.
c   LABEL  (input) : character strings to be used to label the specified
c                    contour. Leading and trailing blank spaces are
c                    ignored.
c   INTVAL (input) : spacing along the contour between labels, in
c                    grid cells.
c   MININT (input) : contours that cross less than MININT cells
c                    will not be labelled.
      INTVAL=100
      MININT=40

      INTVAL=20
      MININT=10
c      do i=1,n_contour,2
       do i=1,n_contour,1

c         write(LABEL,'(1pe7.1)')contour(i)
         write(LABEL,'(i2)')i
c         WRITE(*,*)'i,contour(i),LABEL',i,contour(i),LABEL
         call  PGCONL(f,nxmax,nymax,1,nxmax,1,nymax,contour(i),
     .                tr,LABEL,INTVAL,MININT)
      enddo

      write(text,560)plot_name
 560  format("Contour values:",A)
      RILIN=8.
      CALL PGMTXT('B',RILIN,-0.2,0.0,text)

 570  format(1X,A,1pe10.3)
      do i=1,n_param
        rilin=rilin+1
        write(text,570)name_param(i),param(i)
        CALL PGMTXT('B',RILIN,-0.2,0.0,text)    
      enddo

      j=0
c      do i=1,n_contour,2
       do i=1,n_contour,1
         j=j+1
         write(text,150)i,contour(i)
         CALL PGMTXT('B',rilin+j,0.,0.,text)
      enddo

 150  format('i=',i3,'contour(I)=',1pe8.1)
      return
      end 

      subroutine contour2d_test_1(nxmax,nymax,nxmaxa,nymaxa,f)
      implicit none
      include 'param.i'
        
c-----input
      integer nxmax,nymax,nxmaxa,nymaxa

      real f(nxmaxa,nymaxa)
      
      integer i,j

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
         WRITE(*,*)'contour2d_test_1:TRYING TO PGPAGE at rank=',myrank
         STOP
      endif     
      
      
      write(*,*)' contour2d_test_1 nxmax,nymax',nxmax,nymax
      do i=1,nxmax
         do j=1,nymax
            write(*,*)'i,j,f(i,j)',i,j,f(i,j)
         enddo
      enddo

      return
      end

      subroutine contour2d_1(nxmax,nymax,f,x,y,plot_name,x_name,y_name,
     .n_contour,contour,name_param,param,n_param)
c-----print contours of 2D array f(1:nxmax,1:nymax)
c     small difference with contour2d
c     In contour2d_1 following lines are commented
c      do i=1,n_contour-1
c         if ((contour(i).le.0.).and.(contour(i+1).gt.0.)) then
c            contour(i)=0.
c         endif
c         write(*,*)'contour2d i,contour(i)', i,contour(i)
c      enddo

c     input
c      fmap(1:nxmax,1:nymax) 2D array of the function f(x,y)
c      nxmax,nymax are the dimentions of f
c      x(1:nxmax)  array of x variable
c      y(1:nymax)  array of y variable
c      plot_name  the name of the top of plot
c      x_name     the name of x axis  
c      y_name     the name of y axis
c      n_contour  the number of contours
c      contour(n_contour)      work array for the contour's values
c      n_param is the number of parameters
c      name_param(n_param)     are the names of parameters
c      param(n_param)          are the values of parameters

      implicit none
c-----input
      integer nxmax,nymax,n_contour,n_param
      real f(nxmax,nymax),x(nxmax),y(nymax),contour(n_contour)
      character*(*)plot_name,x_name,y_name
      character*(*)name_param(n_param)
      real param(*)
      real y_e,x_e,cnpar,T
c-----external
      real rbound !convert double to real bounded 1.e-33 and 1.e+33
c-----local
      real tr(6),RILIN
      real xmin,xmax
      real ymin,ymax
      real fmin,fmax
      integer i,j
      
c      CHARACTER*7 LABEL
      CHARACTER*2 LABEL

      integer INTVAL,MININT
      character*72 text             

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0)then
         WRITE(*,*)'contour2d_1: TRYING TO PGPAGE at rank=',myrank
         STOP
      endif     
      
c-----Advance plotter to a new page or panel, clearing the screen if
c     necessary.
      CALL PGPAGE  

c-----Change the size and position of the viewport, specifying
c     the viewport in normalized device coordinates.  Normalized
c     device coordinates run from 0 to 1 in each dimension.
c     SUBROUTINE PGSVP (XLEFT, XRIGHT, YBOT, YTOP)
c     REAL XLEFT, XRIGHT, YBOT, YTOP

c      CALL PGSVP(.2,.8,.3,.95) 
      CALL PGSVP(.2,.8,.5,.95) 

c-----Change the window in world coordinate space that is to be mapped on
c     to the viewport.  Usually PGSWIN is called automatically by PGENV,
c     but it may be called directly by the user.
c     SUBROUTINE PGSWIN (X1, X2, Y1, Y2)
c     REAL X1, X2, Y1, Y2
      call rminmx(x,1,nxmax,1,xmin,xmax)
      if (xmin .eq. xmax) xmin=.9*xmax-1.e-20
      call rminmx(y,1,nymax,1,ymin,ymax)
      if (ymin .eq. ymax) ymin=.9*xmax-1.e-20
      CALL PGSWIN (xmin,xmax,ymin,ymax)

c-----PGBOX -- draw labeled frame around viewport
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)

c-----PGLAB write labels  for x,y axis and top of plot
      CALL PGLAB(x_name,y_name,plot_name)
 
c-----draw a contour map of an array.  The map is truncated if
c     necessary at the boundaries of the viewport.  Each contour line
c     is drawn with the current line attributes (color index, style, and
c     width); except that if argument NC is positive (see below), the line
c     style is set by PGCONT to 1 (solid) for positive contours or 2
c     (dashed) for negative contours.
c
c     Arguments:
c     A      (input) : data array.
c     IDIM   (input) : first dimension of A.
c     JDIM   (input) : second dimension of A.
c     I1, I2 (input) : range of first index to be contoured (inclusive).
c     J1, J2 (input) : range of second index to be contoured (inclusive).
c     C      (input) : array of NC contour levels; dimension at least NC.
c     NC     (input) : +/- number of contour levels (less than or equal
c                      to dimension of C). If NC is positive, it is the
c                      number of contour levels, and the line-style is
c                      chosen automatically as described above. If NC is
c                      negative, it is minus the number of contour
c                      levels, and the current setting of line-style is
c                      used for all the contours.
c     TR     (input) : array defining a transformation between the I,J
c                      grid of the array and the world coordinates.
c                      The world coordinates of the array point A(I,J)
c                      are given by:
c                      X = TR(1) + TR(2)*I + TR(3)*J
c                      Y = TR(4) + TR(5)*I + TR(6)*J
c                      Usually TR(3) and TR(5) are zero - unless the
c                      coordinate transformation involves a rotation or
c                      shear.
c
c     SUBROUTINE PGCONT (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, TR)
c     INTEGER IDIM, JDIM, I1, J1, I2, J2, NC
c     REAL A(IDIM,JDIM), C(*), TR(6)
cSm030515 these lines work for g77 fortran under linux on PC 
      fmin=1.e+30
      fmax=-1.e+30
cSm030515 these lines work for digital fortan  under windows on PC
c      fmin=1.e+38
c      fmax=-1.e+38
      do i=1,nxmax
         do j=1,nymax
            fmin=MIN(f(i,j),fmin)
            fmax=MAX(f(i,j),fmax)
         enddo
      enddo
          
c      write(*,*)'contour.f contour2d fmin,fmax',fmin,fmax      
      do i=1,n_contour
         contour(i)=fmin+(i-1)*(fmax-fmin)/real(n_contour-1)
c         contour(i)=1.e-2*fmin+
c     .   (i-1)*1.e-2*(fmax-fmin)/real(n_contour-1)
      enddo

      do i=1,n_contour-1
c         if ((contour(i).le.0.).and.(contour(i+1).gt.0.)) then
c            contour(i)=0.
c         endif
c         write(*,*)'contour2d i,contour(i)', i,contour(i)
      enddo
      
 10   continue

      tr(1)=xmin-(xmax-xmin)/real(nxmax-1)
      tr(2)=(xmax-xmin)/real(nxmax-1)
      tr(3)=0.
      tr(4)=ymin-(ymax-ymin)/real(nymax-1)
      tr(5)=0.
      tr(6)=(ymax-ymin)/real(nymax-1)
      
      call PGCONT(f,nxmax,nymax,1,nxmax,1,nymax,contour,n_contour,tr)

c---------------------------------------------------------------
c     SUBROUTINE PGCONL (A, IDIM, JDIM, I1, I2, J1, J2, C, TR,
c    1                   LABEL, INTVAL, MININT)
c     INTEGER IDIM, JDIM, I1, J1, I2, J2, INTVAL, MININT
c     REAL A(IDIM,JDIM), C, TR(6)
c     CHARACTER*(*) LABEL
c     Label a contour map drawn with routine PGCONT. Routine PGCONT should
c     be called first to draw the contour lines, then this routine should be
c     called to add the labels. Labels are written at intervals along the
c     contour lines, centered on the contour lines with lettering aligned
c     in the up-hill direction. Labels are opaque, so a part of the under-
c     lying contour line is obscured by the label. Labels use the current
c     attributes (character height, line width, color index, character
c     font).

c     The first 9 arguments are the same as those supplied to PGCONT, and
c     should normally be identical to those used with PGCONT. Note that
c     only one contour level can be specified; tolabel more contours, call
c     PGCONL for each level.

c     The Label is supplied as a character string in argument LABEL.

c     The spacing of labels along the contour is specified by parameters
c     INTVAL and MININT. The routine follows the contour through the
c     array, counting the number of cells that the contour crosses. The
c     first label will be written in the MININT'th cell, and additional
c     labels will be written every INTVAL cells thereafter. A contour
c     that crosses less than MININT cells will not be labelled. Some
c     experimentation may be needed to get satisfactory results; a good
c     place to start is INTVAL=20, MININT=10.

c     Arguments:
c     A      (input) : data array.
c     IDIM   (input) : first dimension of A.
c     JDIM   (input) : second dimension of A.
c     I1, I2 (input) : range of first index to be contoured (inclusive).
c     J1, J2 (input) : range of second index to be contoyred 
c     C      (input) : the level of the contour to be labelled (one of the
c                      values given to PGCONT).
c      TR    (input) : array defining a transformation between the I,J
c                      grid of the array and the world coordinates.
c                      The world coordinates of the array point A(I,J)
c                      are given by:
c                      X = TR(1) + TR(2)*I + TR(3)*J
c                      Y = TR(4) + TR(5)*I + TR(6)*J
c                      Usually TR(3) and TR(5) are zero - unless the
c                      coordinate transformation involves a rotation or
c                      shear.
c   LABEL  (input) : character strings to be used to label the specified
c                    contour. Leading and trailing blank spaces are
c                    ignored.
c   INTVAL (input) : spacing along the contour between labels, in
c                    grid cells.
c   MININT (input) : contours that cross less than MININT cells
c                    will not be labelled.
      INTVAL=100
      MININT=40

      INTVAL=20
      MININT=10
c      do i=1,n_contour,2
       do i=1,n_contour,1

c         write(LABEL,'(1pe7.1)')contour(i)
         write(LABEL,'(i2)')i
c         WRITE(*,*)'i,contour(i),LABEL',i,contour(i),LABEL
         call  PGCONL(f,nxmax,nymax,1,nxmax,1,nymax,contour(i),
     .                tr,LABEL,INTVAL,MININT)
      enddo

      write(text,560)plot_name
 560  format("Contour values:",A)
      RILIN=8.
      CALL PGMTXT('B',RILIN,-0.2,0.0,text)

 570  format(1X,A,1pe10.3)
      do i=1,n_param
        rilin=rilin+1
        write(text,570)name_param(i),param(i)
        CALL PGMTXT('B',RILIN,-0.2,0.0,text)    
      enddo

      j=0
c      do i=1,n_contour,2
       do i=1,n_contour,1
         j=j+1
         write(text,150)i,contour(i)
         CALL PGMTXT('B',rilin+j,0.,0.,text)
      enddo

 150  format('i=',i3,'contour(I)=',1pe8.1)
      return
      end 
