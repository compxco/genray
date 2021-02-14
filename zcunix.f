c     
c
      real*8 function taper(x,x0,x1,x2)
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none !YuP[2020-01-14]
      real*8 x,x0,x1,x2, pi, xx  !YuP[2020-01-14]

c     A multiplication factor between 0. and 1., giving a tapered
c     (i.e., smoothed) box function.
c
c     For abs(x-x0) .lt. x1/2., gives 1.
c     For (x1+x2)/2 .gt. abs(x-x0) .gt. x1/2., gives 
c                                  0.5*(1.+cos(abs(x-x1)/x2*pi))
c     For abs(x-x0) .ge. (x1+x2)/2., gives 0.
c     
c     x1 and x2 are presumed .ge.0.

      data pi/3.141592653589793/

      xx=abs(x-x0)
      
      if (xx.lt. x1/2.) then
        taper=1.
      elseif (xx.gt.x1/2. .and. xx.lt.(x1+x2)/2.) then
        taper=0.5*(1.+cos(xx/x2*pi))
      else
        taper=0.0
      endif

      return
      end
c 
c     subroutine aminmx(array,ifirst,ilast,istride,amin,amax,
c     *indmin,indmax)
c     real array(1)
c     length = 1+(ilast+1-ifirst)/istride
c     if(ilast.lt.ifirst+(length-1)*istride) length=length-1
cc    k=istride
c     m=ismin(length,array(ifirst),istride)
c     m=ifirst+(m-1)*istride
cc    m=ifirst
cc    do 1 i=ifirst+k,ilast,k
cc1   if(array(i).le.array(m)) m=i
c     amin=array(m)
c     m=ismax(length,array(ifirst),istride)
c     m=ifirst+(m-1)*istride
cc    m=ifirst
cc    do 2 i=ifirst+k,ilast,k
cc2   if(array(i).ge.array(m)) m=i
c     amax=array(m)
c     return
c     end

c***********************************************************************
c     SPLINE PACKAGE
c***********************************************************************
c     package cubspl         note--documentation for individual routines
c     follows the general package information
c
c     latest revision        january 1985
c
c     purpose:         to perform one and two-dimensional cubic spline
c     interpolation with choice of boundary
c     conditions.  the function and selected
c     derivatives may be evaluated at any point where
c     interpolation is required.  in the
c     two-dimensional case, the given data points
c     must be on a rectangular grid, which need not
c     be equally spaced.  the package cubspl contains
c     seven routines.
c
c     subroutine coeff1
c     computes the coefficients for one-dimensional
c     cubic spline interpolation using one of
c     the following boundary conditions at
c     each end of the range.
c     . second derivative given at boundary.
c     . first derivative given at boundary.
c     . periodic boundary condition.
c     . first derivative determined by fitting a
c     cubic to the four points nearest to the
c     boundary.
c
c     subroutine terp1
c     using the coefficients computed by coeff1,
c     this routine evaluates the function and/or
c     first and second derivatives at any point
c     where interpolation is required.
c
c     subroutine coeff2
c     computes the coefficients for two-dimensional
c     bicubic spline interpolation with the same
c     choice of boundary conditions as for coeff1.
c
c     function terp2
c     using the coefficients produced by coeff2,
c     this subroutine evaluates the function or a
c     selected derivative at any point where
c     two-dimensional interpolation is required.
c
c     subroutine trip
c     a simple periodic, tridiagonal linear
c     equation solver used by coeff1.
c
c     subroutine search
c     performs a binary search in a one-dimensional
c     floating point table arranged in ascending
c     order.  this routine is called by terp1 and
c     terp2.
c
c     subroutine intrp
c     given coefficients provided by coeff1 and the
c     position of the interpolation point in the
c     independent variable table, this routine
c     performs one-dimensional interpolation for
c     the function value, first and second
c     derivative, as desired.  this routine is
c     called by terp1 and terp2.
c
c     usage:      for one-dimensional cubic spline interpolation,
c     the user first calls coeff1 by
c
c     call coeff1 (n,x,f,w,iop,int,wk)
c
c     this subroutine returns the coefficients
c     needed for the subsequent interpolation
c     in the array w.  the user then calls
c     subroutine terp1 by
c
c     call terp1 (n,x,f,w,y,int,tab,itab)
c
c     in order to compute the value of the
c     function and/or its derivatives.  the user
c     specifies y, the value of the independent
c     variable where the interpolation is to be
c     performed.  the interpolated values are
c     returned in tab.  the parameters
c     n,x,f,w, and int must not be changed
c     between the subroutine calls.
c
c     for two-dimensional cubic spline interpolation
c     the user first calls coeff2 by
c
c     call coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,
c     idm,ibd,wk)
c
c     this subroutine returns the coefficients
c     needed for the subsequent interpolation in
c     the arrays fxx, fyy, fxxyy.  the user then
c     calls the routine terp2 by
c
c     r = terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,
c     idm,ixd,iyd)
c
c     in order to perform the interpolation at the
c     point specified by the values of xb and yb.
c     depending on the values input in ixd and iyd,
c     the routine returns an interpolated value
c     for the function or one of its partial
c     derivatives.  the parameters nx,x,ny,y,f,fxx,
c     fyy,fxxyy, and idm must not be changed
c     between the calls to coeff2 and terp2.
c
c     special conditions     tables of independent variables must be
c     arranged in ascending order.  for
c     two-dimensional interpolation, the data points
c     must lie on a rectangular mesh.
c
c     i/o                    none
c
c     precision              single
c
c     files                  cray machines.
c
c     language               fortran
c
c     history                this package is based on the routines
c     la e 102a, spl1d1
c     la e 103a, spl1d2
c     la e 104a, spl2d1
c     la e 105a, spl2d2
c     of the los alamos cubic spline package written
c     by thomas j. jordan and bertha fagan, 1968.
c     the routines have been streamlined and
c     standardized.  the algorithm for
c     two-dimensional interpolation is considerably
c     modified.
c
c     algorithm:      for one-dimensional interpolation, the cubic
c     spline is expressed in terms of the function
c     values and second derivatives at the data
c     points.  the second derivatives are
c     determined from a tridiagonal linear system
c     which describes the continuity of the first
c     derivative and incorporates the given
c     boundary conditions.  coeff1  sets up this
c     system and calls subroutine trip to solve it.
c
c     the cubic segment between two adjacent
c     tabular points is constructed from the
c     function values and second derivatives at
c     these points.  these provide the four
c     constants needed to define the cubic
c     uniquely.  from this cubic, values of the
c     function and its first and second
c     derivatives are readily determined at any
c     intermediate point.  one-dimensional
c     interpolation is performed by the routine
c     terp1.  for two-dimensional interpolation,
c     the bicubic spline is described in terms of
c     values of f,fxx,fyy, and fxxyy  at each
c     point on the given two-dimensional
c     rectangular grid of data points.  here f
c     is the function value,
c
c     fxx = (d/dx)**2*f
c
c     and so on.  the coefficients are determined
c     by coeff2, which uses successive
c     applications of coeff1.
c
c     1.  the array fxx is determined by applying
c     coeff1 to f along each line in the
c     x-direction.
c
c     2.  the array fyy is determined by applying
c     coeff1 to f along each line in the
c     y-direction.
c
c     3.  fxxyy is determined on the constant y
c     boundaries by applying coeff1 to fyy
c     along these boundaries.
c
c     4.  the remainder of the array fxxyy is
c     determined by applying coeff1 to fxx
c     along each line in the y-direction.
c
c     the bicubic within any rectangular element
c     of the grid is constructed from the values
c     of f,fxx,fyy,fxxyy at the four corners.
c     these provide the 16 constants necessary
c     to define the bicubic uniquely.  to find
c     the value of f corresponding to a point
c     (xb,yb) within the elementary rectangle,
c     (x(i),y(j)),(x(i+1),y(j)),(x(i),y(j+1)),
c     (x(i+1),y(j+1)), five one dimensional
c     interpolations are performed.
c
c     1.  f at (xb,y(j)) and (xb,y(j+1)) are
c     found by interpolating f in the
c     x-direction using fxx. (two interpolations)
c
c     2.  fyy at (xb,y(j)) and (xb,y(j+1)) are
c     found by interpolating fyy in the
c     x-direction using fxxyy. (two
c     interpolations.)
c
c     3.  finally f at (xb,yb) is found by
c     interpolating between f(xb,y(j)) and
c     f(xb,y(j+1)) in the y-direction using
c     values of fyy(xb,y(j)) and fyy(xb,y(j+1))
c     obtained above. (one interpolation).
c
c     two-dimensional interpolation is performed
c     in terp2.
c
c     references             for greater detail, see j.l.walsh,
c     j.h.ahlberg, e.n.nilsen, best approximation
c     properties of the spline fit, journal of
c     mathematics and mechanics, vol.ii(1962),
c     225-234.
c
c     t.l. jordan, smoothing and multivariable
c     interpolation with splines, los alamos
c     report, la-3137, 1965.
c
c     accuracy               near machine accuracy was obtained when
c     interpolating a cubic in one dimension
c     or a bicubic in two dimensions.
c
c     portability            fully portable with respect to fortran 66.
c***********************************************************************
c
c     subroutine coeff1 (n,x,f,w,iop,int,wk)
c
c     dimension of           x(n),f(int*(n-1)+1),w(int*(n-1)+1),iop(2),
c     arguments              wk(3*n+1)
c
c     purpose:        subroutine coeff1 computes the coefficients
c     for one-dimensional cubic spline
c     interpolation using one of the following
c     boundary conditions at each end of the
c     range
c
c     .  second derivatives given at boundary
c     .  first derivative given at boundary
c     .  periodic boundary conditions
c     .  first derivative calculated by fitting a
c     cubic to the four points nearest to the
c     boundary
c
c     note that terp1 must be called to perform
c     the interpolation.
c
c     usage                  call coeff1 (n,x,f,w,iop,int,wk)
c
c     arguments
c
c     on input               n
c     the number of data points.  n must be at
c     least 4.
c
c     x
c     table of n independent variable values in
c     ascending order.  dimension of x in calling
c     program must be at least n.
c
c     f
c     table of n corresponding dependent variable
c     values.  the values are separated by interval
c     int.  this is usually unity for
c     one-dimensional interpolation.  other values
c     are useful when coeff1 is incorporated in a
c     two-dimensional interpolation scheme (see
c     coeff2).  dimension of f in the calling
c     program must be at least (int*(n-1)+1).
c
c     iop
c     two element integer array defining boundary
c     conditions at x(1) and x(n) according to the
c     following code.
c
c     for iop(1)
c     = 1  second derivative given at x(1).  place
c     value of second derivative in w(1)
c     before call to coeff1.
c     = 2  first derivative given at x(1).  place
c     value of first derivative in w(1) before
c     call.
c     = 3  periodic boundary condition.  x(1) and
c     x(n) are equivalent points.  f(1) and
c     f(int*(n-1)+1) are equal.
c     = 4  the first derivative at x(1) is
c     calculated by fitting a cubic to points
c     x(1) through x(4).
c     similarly, iop(2) defines the boundary
c     condition at x(n).  when iop(2) = 1 (or 2),
c     the value of the second (or first) derivative
c     must be placed in w(int*(n-1)+1).  note that
c     if iop(1) = 3, consistency demands that
c     iop(2) = 3 also.
c
c     int
c     spacing in f and w tables.  for
c     one-dimensional interpolation this will
c     usually be unity.
c
c     wk
c     work area of dimension at least (3*n+1).
c
c     on output              w
c     table of second derivatives corresponding to
c     given x and f values.  the separation of
c     tabular entries is int (see above).
c     dimension of w in the calling program must be
c     at least (int*(n-1)+1).
c
c     the arrays x, f, w are used as input for the
c     routine terp1, which performs interpolation
c     at a given value of the independent variable.
c
c     timing:       the timing is linearly proportional to n, the
c     number of data points.
c***********************************************************************
c
c     subroutine terp1 (n,x,f,w,y,int,tab,itab)
c
c
c     dimension of           x(n),f(int*(n-1)+1),w(int*(n-1)+1),tab(3),
c     arguments              itab(3)
c
c     purpose                using the coefficients computed by coeff1,
c     this routine evaluates the function and/or
c     first and second derivatives at any point.
c
c     usage                  call terp1 (n,x,f,w,y,int,tab,itab)
c
c     arguments
c
c     on input               n
c     the number of data points.  n must be at
c     least 4.
c
c     x
c     table of n independent variable values in
c     ascending order.  dimension of x in the
c     calling program must be at least n.
c
c     f
c     table of n corresponding dependent variable
c     values separated by interval int, usually
c     unity for one-dimensional interpolation.
c     dimension of f in the calling program must be
c     at least (int*(n-1)+1).
c
c     w
c     table of second derivatives computed by
c     coeff1.  the separation of tabular entries is
c     int.  dimension of w in the calling program
c     must be at least (int*(n-1)+1).
c
c     y
c     value of the independent variable at which
c     interpolation is required.  if y lies outside
c     the range of the table, extrapolation takes
c     place.
c
c     int
c     spacing of tabular entries in f and w arrays.
c     this is usually unity for one-dimensional
c     interpolation.
c
c     itab
c     three element integer array defining
c     interpolation to be performed at y.
c     if itab(1) = 1, the function value is
c     returned in tab(1).
c     if itab(2) = 1, the first derivative is
c     returned in tab(2).
c     if itab(3) = 1, the second derivative is
c     returned in tab(3).
c     if itab(i) = 0 for i = 1, 2 or 3, the
c     corresponding function value or derivative is
c     not computed and tab(i) is not referenced.
c
c     on output              tab
c     three element array in which interpolated
c     function value, first and second derivatives
c     are returned as dictated by itab (see above).
c
c     timing:       this procedure is fast.  the maximum time for
c     the binary search is proportional to alog(n).
c     the time for function and derivative evaluation
c     is independent of n.
c***********************************************************************
c
c     subroutine coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ibd,wk)
c
c
c     dimensions:        x(nx),y(ny),f(idm,ny),fxx(idm,ny),fyy(idm,ny),
c     arguments              fxxyy(idm,ny),ibd(4),wk(3*max0(nx,ny)+1)
c     (idm must be .ge. nx)
c
c     purpose:        subroutine coeff2 computes the coefficients
c     for two-dimensional bicubic spline
c     interpolation with the same choice of
c     boundary conditions as for coeff1.  terp2
c     is called to perform interpolation.
c
c     usage:       call coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ibd,wk)
c
c     arguments
c
c     on input               nx
c     number of grid points in the x-direction.  nx
c     must be at least 4.
c
c     x
c     table of nx values of the first independent
c     variable arranged in ascending order.
c     dimension of x in the calling program must be
c     at least nx.
c
c     ny
c     number of grid points in the y-direction.  ny
c     must be at least 4.
c
c     y
c     table of ny values of the second independent
c     variable arranged in ascending order.
c     dimension of y in the calling program must be
c     at least ny.
c
c     f
c     two dimensional array of function values at
c     the grid points defined by the arrays x and
c     y.  dimension of f in the calling program is
c     (idm, nyy) where
c     idm .ge. nx
c     nyy .ge. ny
c
c     idm
c     first dimension in the calling program of
c     arrays f, fxx, fyy, fxxyy.  idm must be at
c     least nx.
c
c     ibd
c     four element integer array defining boundary
c     conditions according to the following code.
c     for ibd(1)
c     = 1  the second derivative of f with respect
c     to x is given at (x(1),y(j)) for
c     j = 1,ny,1.  values of this second
c     derivative must be placed in fxx(1,j)
c     for j = 1,ny,1, before calling coeff2.
c     = 2  the first derivative of f with respect
c     to x is given at (x(1),y(j)) for
c     j = 1,ny,1.  values of the derivative
c     must be placed in fxx(1,j) for
c     j = 1,ny,1 before calling coeff2.
c     = 3  periodic boundary condition in the
c     x-direction.  (x(1),y(j)) and
c     and (x(nx),y(j)) are equivalent points
c     for j = 1,ny,1.  f(1,j) and f(nx,j) are
c     equal.
c     = 4  the first derivative of f with respect
c     to x at (x(1),y(j)) is computed by
c     fitting a cubic to f(1,j) through f(4,j)
c     for j = 1,ny,1.
c
c     similarly, ibd(2) defines the boundary
c     condition at (x(nx),y(j)) for j = 1,ny,1.
c     when ibd(2) = 1 (or 2) the values of the
c     second (or first) derivative of f with
c     respect to x are placed in fxx(nx,j) for
c     j = 1,ny,1.
c     note that if ibd(1) = 3, consistency
c     requires that ibd(2) = 3 also.
c     for ibd(3)
c     = 1  the second derivative of f with respect
c     to y is given at (x(i),y(1)).  place
c     values of the derivative in fyy(i,1) for
c     i = 1,nx,1 before calling coeff2.
c     = 2  the first derivative of f with respect
c     to y is given at (x(i),y(1)).  values of
c     this derivative must be placed in
c     fyy(i,1) for i = 1,nx,1 before calling
c     coeff2.
c     = 3  periodic boundary condition in the
c     y-direction.  (x(i),y(1)) and
c     (x(i),y(ny)) are equivalent points.
c     f(i,1) and f(i,ny) are equal.
c     = 4  the first derivative of f with respect
c     to y at (x(i),y(1)) is computed by
c     fitting a cubic to f(i,1) through f(i,4)
c     for i = 1,nx,1.
c
c     similary, ibd(4) defines the boundary
c     condition at (x(i),y(ny)) for i = 1,nx,1 and
c     given derivative values are placed in
c     fyy(i,ny).
c     note that consistency demands that if
c     ibd(3) = 3, then ibd(4) = 3 also.
c
c     wk
c     work area of dimension at least
c     (3*max0(nx,ny)+1)
c
c     on output              fxx
c     array of second derivatives of f with respect
c     to x computed by coeff2.  fxx(i,j) is
c     derivative at (x(i),y(j)).  as for f,
c     dimension of fxx in the calling program is
c     (idm,nyy).
c
c     fyy
c     array of second derivatives of f with respect
c     to y computed by coeff2.  dimension of fyy in
c     the calling program is (idm,nyy).
c
c     fxxyy
c     array of fourth derivatives
c     (d/dx)**2*(d/dy)**2*f, computed by coeff2.
c     dimension of fxxyy in the calling program is
c     (idm,nyy).
c
c     the arrays x, y, f, fxx, fyy, fxxyy are used as
c     input for the routine terp2 which performs
c     interpolation at required values of the two
c     independent variables.
c
c     timing                 the timing is proportional to nx*ny.
c***********************************************************************
c
c     function terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd)
c
c
c     dimensions:        x(nx),y(ny),f(idm,ny),fxx(idm,ny),fyy(idm,ny),
c     arguments              fxxyy(idm,ny))
c     (idm must be .ge. nx)
c
c     purpose                using the coefficients produced by coeff2,
c     this routine evaluates the function on a
c     selected derivative of any point where
c     two-dimensional interpolation is required.
c
c     usage:          r = terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,
c     ixd,iyd)
c
c     arguments
c
c     on input               xb, yb
c     values of the independent variables, x and y,
c     at which interpolation is required.
c
c     nx
c     number of grid points in the x-direction.  nx
c     must be at least 4.
c
c     x
c     table of nx values of the independent
c     variable, x, arranged in ascending order.
c     dimension of x in the calling program must be
c     at least nx.
c
c     ny
c     number of grid points in the y-direction.  ny
c     must be at least 4.
c
c     y
c     table of ny values of the independent
c     variable, y, arranged in ascending order.
c     dimension of y in the calling program must be
c     at least ny.
c
c     f
c     two-dimensional array of function values at
c     grid points defined by the arrays x and y.
c     dimension of f in the calling program is
c     (idm,nyy), where
c     idm .ge. nx
c     nyy .ge. ny
c
c     fxx
c     array of second derivatives of f with respect
c     to x computed by coeff2.  dimension of fxx in
c     the calling program is (idm,nyy).  see under
c     f above.
c
c     fyy
c     array of second derivatives of f with respect
c     to y computed by coeff2.  dimension of fyy in
c     the calling program is (idm,nyy).
c
c     fxxyy
c     array of fourth derivatives,
c     (d/dx)**2*(d/dy)**2*f, computed by coeff2.
c     dimension of fxxyy in the calling program is
c     (idm,nyy).
c
c     idm
c     first dimension in the calling program of
c     arrays f, fxx, fyy and fxxyy,
c     idm .ge. nx
c
c     ixd, iyd
c     define derivative to be returned by the
c     function terp2.  ixd, iyd may each take the
c     the values 0, 1, 2.  the derivative returned
c     is (d/dx)**ixd*(d/dy)**iyd*f.
c     note that if ixd = iyd = 0, the function
c     value itself is returned.
c
c     timing                 this procedure is fast.  the maximum
c     time for the binary search is proportional to
c     alog(nx*ny).  the time for function evaluation
c     is independent of n.
c***********************************************************************


      subroutine coeff1 (n,x,f,w,iop,int,wk)
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      integer n,iop,int
      integer i,i1,i2,ii,j0,j1,j2,j3,j4,jm,ml,mk,nn,in,iw1,index
      real*8 x,f,w,wk,y2,b2,a12,a13,a14,a23,a24,a34
      
cSm030221 
c     dimension       x(n)       ,f(n*int)       ,w(n*int)      ,iop(2),
c     1  wk(n,*)
      dimension       x(*)       ,f(*)       ,w(*)      ,iop(2),
     1  wk(n,*)
!YuP cdir$ nobounds
      logical q8q4
      save q8q4
      data q8q4 /.true./
c
c     arithmetic statenent function used to locate entries in 
c     f and w arrays:
c
      ii(index)=(index-1)*int+1
c
c
c
c
c
c
c     the following call is for gathering statistics at ncar
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     start to set up tridiagonal system
c
      j0 = 1
      do 101 i=2,n
        jm = j0
        j0 = j0+int
        wk(i,1) = x(i)-x(i-1)
        wk(i,2) = (f(j0)-f(jm))/wk(i,1)
        wk(i,3) = wk(i,1)/6.
        wk(i,1) = wk(i,1)/3.
 101  continue

c      write(*,*)'coeff1 wk(i,1)',(wk(i,1),i=1,n)
c      write(*,*)'coeff1 wk(i,2)',(wk(i,2),i=1,n) 
c      write(*,*)'coeff1 wk(i,3)',(wk(i,3),i=1,n)

      nn = n
      mk = iop(1)
      ml = iop(2)
c
c     apply boundary conditions at boundary 1
c

c      write(*,*)'coeff1 mk,ml',mk,ml

      go to (102,103,104,105),mk
c
c     second derivative given at boundary 1
c
 102  continue
      wk(2,2) = wk(3,2)-wk(2,2)-wk(2,3)*w(1)
      wk(2,3) = 0.
      wk(2,1) = wk(2,1)+wk(3,1)
      i1 = 2
      nn = nn-1
      go to 106
c
c     first derivative given at boundary 1
c
 103  continue
      wk(1,2) = wk(2,2)-w(1)
      wk(2,2) = wk(3,2)-wk(2,2)
      wk(1,3) = 0.
      wk(1,1) = wk(2,1)
      wk(2,1) = wk(2,1)+wk(3,1)
      i1 = 1
      go to 106
c
c     periodic boundary condition
c
 104  continue
      y2 = wk(2,2)
      b2 = wk(2,1)
      wk(2,2) = wk(3,2)-wk(2,2)
      wk(2,1) = wk(3,1)+wk(2,1)
      i1 = 2
      nn = nn-1
      go to 106
c
c     first derivative at boundary 1 from 4 point interpolation.
c
 105  continue
      a12 = x(1)-x(2)
      a13 = x(1)-x(3)
      a14 = x(1)-x(4)
      a23 = x(2)-x(3)
      a24 = x(2)-x(4)
      a34 = x(3)-x(4)
      j1 = 1
      j2 = j1+int
      j3 = j2+int
      j4 = j3+int
      w(1)    = (1./a12+1./a13+1./a14)*f(j1)-
     1  a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)-
     2  a12*a13/(a14*a24*a34)*f(j4)
      go to 103
c     compute tridiagonal arrays
 106  continue
      i2 = n-2
      do 107 i=3,i2
        wk(i,2) = wk(i+1,2)-wk(i,2)
        wk(i,1) = wk(i+1,1)+wk(i,1)
 107  continue
c
c     apply boundary conditions at boundary 2.
c
      in = ii(n)
      go to (108,109,110,111),ml
c
c     second derivative given at boundary 2.
c
 108  continue
      wk(n-1,2) = wk(n,2)-wk(n-1,2)-wk(n,3)*w(in)
      wk(n,3) = 0.
      wk(n-1,1) = wk(n-1,1)+wk(n,1)
      nn = nn-1
      go to 112
c
c     first derivative given at boundary 2.
c
 109  continue
      wk(n-1,2) = wk(n,2)-wk(n-1,2)
      wk(n,2) = -wk(n,2)+w(in)
      wk(n-1,1) = wk(n-1,1)+wk(n,1)
      wk(1,4) = 0.
      go to 112
c
c     periodic boundary condition
c
 110  continue
      wk(n-1,2) = wk(n,2)-wk(n-1,2)
      wk(n,2) = y2-wk(n,2)
      wk(n-1,1) = wk(n-1,1)+wk(n,1)
      wk(n,1) = wk(n,1)+b2
      wk(1,4) = wk(2,3)
      go to 112
c
c     first derivative at boundary 2 from 4 point interpolation.
c
 111  continue
      a12 = x(n)-x(n-1)
      a13 = x(n)-x(n-2)
      a14 = x(n)-x(n-3)
      a23 = x(n-1)-x(n-2)
      a24 = x(n-1)-x(n-3)
      a34 = x(n-2)-x(n-3)
      j1 = in
      j2 = j1-int
      j3 = j2-int
      j4 = j3-int
      w(in)   = (1./a12+1./a13+1./a14)*f(j1)-
     1  a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)-
     2  a12*a13/(a14*a24*a34)*f(j4)
      go to 109
 112  continue
      iw1 = ii(i1)
      call trip (nn,wk(i1,3),wk(i1,1),wk(i1+1,3),wk(i1,2),w(iw1),int)
      go to (114,114,113,114),mk
 113  continue
      w(1) = w(in)
 114  continue
!YuP cdir$ bounds
      return
      end


      subroutine coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ibd,wk)
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      integer i,j,ny,nx,idm,ibd,iloc,jloc
      real*8 x,y,f,fxx,fyy,fxxyy,wk
c
cSm030221
c      dimension       x(nx)       ,y(ny)       ,f(idm,ny)  ,fxx(idm,ny),
c     1  fyy(idm,ny) ,fxxyy(idm,ny)           ,ibd(4)     ,
c     2  iloc(2)    ,jloc(2)
      dimension       x(*)       ,y(*)       ,f(idm,*)  ,fxx(idm,*),
     1  fyy(idm,*) ,fxxyy(idm,*)           ,ibd(4)     ,
     2  iloc(2)    ,jloc(2),
cSAP091201
     &  wk(*)
      logical q8q4
      save q8q4
      data q8q4 /.true./
      data iloc(1),iloc(2),jloc(1),jloc(2)/1,1,4,4/
c     the following call is for gathering statistics at ncar
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     compute fxx
c
      do 101 j=1,ny

c        write(*,*)'coeff 2 101 j',j
c        write(*,*)'coeff2 nx',nx
c        write(*,*)'coeff2 x',x
c        write(*,*)'coeff2 f(i,j)',(f(i,j),i=1,nx)
c        write(*,*)'coeff2 ibd',ibd
c        write(*,*)'coeff2 wk',wk
cBH020822
cBH020822c        call coeff1 (nx,x,f(1,j),fxx(1,j),ibd(1),1,wk)
         call coeff1 (nx,x,f(1,j),fxx(1,j),ibd(1),1,wk)
cBH020822
c        write(*,*)'coeff2 fxx(1,1)',fxx(1,1)

c
 101  continue
    
c
c     compute fyy
c
      do 102 i=1,nx
        call coeff1 (ny,y,f(i,1),fyy(i,1),ibd(3),idm,wk)
 102  continue
c
c     check for periodic boundary condition in both directions
c
      if (ibd(1) .eq. 3) go to 103
      if (ibd(3) .eq. 3) go to 105
c
c     calculate fxxyy along left and right boundaries
c
      call coeff1 (ny,y,fxx(1,1),fxxyy(1,1),jloc,idm,wk)
      call coeff1 (ny,y,fxx(nx,1),fxxyy(nx,1),jloc,idm,wk)
      go to 106
 103  continue
c
c     periodic in x direction . calculate fxxyy along lower and upper
c     boundaries.
c
      call coeff1 (nx,x,fyy(1,1),fxxyy(1,1),ibd(1),1,wk)
      call coeff1 (nx,x,fyy(1,ny),fxxyy(1,ny),ibd(1),1,wk)
c
c     calculate remaining fxxyy
c
      do 104 i=1,nx
        call coeff1 (ny,y,fxx(i,1),fxxyy(i,1),iloc,idm,wk)
 104  continue
      go to 108
 105  continue
c
c     periodic in y direction. calculate fxxyy along left and right
c     boundaries.
c
      call coeff1 (ny,y,fxx(1,1),fxxyy(1,1),ibd(3),idm,wk)
      call coeff1 (ny,y,fxx(nx,1),fxxyy(nx,1),ibd(3),idm,wk)
 106  continue
c
c     calculate remaining fxxyy
c
      do 107 j=1,ny
        call coeff1 (nx,x,fyy(1,j),fxxyy(1,j),iloc,1,wk)
 107  continue
 108  continue
      return
      end


      subroutine intrp (n,x,f,w,y,i,int,tab,itab)
      implicit integer (i-n), real*8 (a-h,o-z)
cSm030221
c      dimension       x(i+1)    ,f(i*int+1)    ,w(i*int+1)  ,tab(3)
c     -  ,itab(3)
      dimension       x(*)    ,f(*)    ,w(*)  ,tab(3)
     -  ,itab(3)

c
c     arithmetic statement function used to locate entries in 
c     f and w arrays:
c
      ii(index)=(index-1)*int+1
c
c     perform interpolation or extrapolation
c
      flk = x(i+1)-x(i)
      flp = x(i+1)-y
      fl0 = y-x(i)

c      write(*,*)'intrp flk,flp,fl0', flk,flp,fl0

      i0 = ii(i)
      ip = i0+int

c      write(*,*)'i0,ip,itab(1)',i0,ip,itab(1)

      if (itab(1) .ne. 0) go to 101
      go to 102
 101  continue
c
c     calculate f(y)
c
      a = (w(i0)*flp**3+w(ip)*fl0**3)/(6.*flk)

c      write(*,*)'i0,ip,w(i0),w(ip),flp,fl0,flk,a',
c     +i0,ip,w(i0),w(ip),flp,fl0,flk,a

      b = (f(ip)/flk-w(ip)*flk/6.)*fl0
      c = (f(i0)/flk-w(i0)*flk/6.)*flp
      tab(1) = a+b+c

c      write(*,*)'101 a,b,c tab(1)',a,b,c,tab(1)

 102  continue
      if (itab(2) .ne. 0) go to 103
      go to 104
 103  continue
c
c     calculate first derivative at y
c
      a = (w(ip)*fl0**2-w(i0)*flp**2)/(2.*flk)
      b = (f(ip)-f(i0))/flk
      c = (w(i0)-w(ip))*flk/6.
      tab(2) = a+b+c
c      write(*,*)'103 a,b,c,tab(2)',a,b,c,tab(2)
 104  continue
      if (itab(3) .ne. 0) go to 105
      go to 106
 105  continue
c
c     calculate second derivative at y
c
      tab(3) = (w(i0)*flp+w(ip)*fl0)/flk
c      write(*,*)'105 tab(3),flk',tab(3),flk
 106  continue
      return
      end


      subroutine search (xbar,x,n,i)
      implicit integer (i-n), real*8 (a-h,o-z)
      dimension       x(n)
      data b/.69314718/
c
c     if xbar is outside range of x table extrapolate
c
      if (xbar .gt. x(2)) go to 101
      i = 1
      return
 101  continue
      if (xbar .lt. x(n-1)) go to 102
      i = n-1
      return
 102  continue
c
c     find maximum power of two less than n
c
      m = int((alog(float(n)))/b)
      i = 2**m
      if (i .ge. n) i = i/2
      k = i
      nm1 = n-1
c
c     conduct binary search.
c
 103  continue
      k = k/2
      if (xbar .ge. x(i)) go to 104
      i = i-k
      go to 103
 104  continue
      if (xbar .le. x(i+1)) return
      i = min0(i+k,nm1)
      go to 103
      end


      subroutine terp1 (n,x,f,w,y,int,tab,itab)
      implicit integer (i-n), real*8 (a-h,o-z)
cSm030221
c      dimension       x(n)       ,f(n*int)       ,w(n*int)    ,tab(3),
c     1  itab(3)
      dimension       x(*)       ,f(*)       ,w(*)    ,tab(3),
     1  itab(3)

c     the following call is for gathering statistics at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     perform search
c
      call search (y,x,n,i)
c
c     carry out interpolation (or extrapolation)
c
      call intrp (n,x,f,w,y,i,int,tab,itab)
      
      return
      end


      real*8 function terp2 
     +     (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd)
      implicit integer (i-n), real*8 (a-h,o-z)
cSm030221
c      dimension       x(nx)      ,y(ny)      ,f(idm,ny)  ,fxx(idm,ny),
c     1  fyy(idm,ny) ,fxxyy(idm,ny)           ,ff(2)      ,
c     2  ww(2)      ,tab(3)     ,itab(3)
      dimension       x(*)      ,y(*)      ,f(idm,*)  ,fxx(idm,*),
     1  fyy(idm,*) ,fxxyy(idm,*)           ,ff(2)      ,
     2  ww(2)      ,tab(3)     ,itab(3)

c     the following call is for gathering statistics 
c       on library use at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     search in x and y arrays.
c
      call search (xb,x,nx,i)
      call search (yb,y,ny,j)
c
c     interpolate in x direction
c
      do 101 i1=1,3
        itab(i1) = 0
 101  continue
      i1 = ixd+1
      itab(i1) = 1
      do 102 j1=1,2
        jj = j+j1-1
        call intrp (n,x,f(1,jj),fxx(1,jj),xb,i,1,tab,itab)
        ff(j1) = tab(i1)
        call intrp (n,x,fyy(1,jj),fxxyy(1,jj),xb,i,1,tab,itab)
        ww(j1) = tab(i1)
 102  continue
c
c     interpolate in y direction
c
      do 103 j1=1,3
        itab(j1) = 0
 103  continue
      j1 = iyd+1
      itab(j1) = 1
      call intrp (2,y(j),ff,ww,yb,1,1,tab,itab)
      terp2 = tab(j1)
      return
c
c     revision history---
c
c     june 1977        replaced non-standard statement functions and
c     subscripts to enhance portability.
c
c     january 1978     deleted references to the  *cosy  cards, moved
c     the revision histories to appear before the
c     final end card, and moved the initial comment
c     cards to appear after the first subroutine card
c     and changed  itab  from logical to integer in
c     subroutine intrp and corrected problem with
c     version numbers in one statistics call
c-----------------------------------------------------------------------
      end


      subroutine searche (xbar,x,n,i,dx)
      implicit integer (i-n), real*8 (a-h,o-z)
cSm030221
c      dimension       x(n)
      dimension       x(*)

      data b/.69314718/
c
c     if xbar is outside range of x table extrapolate
c
      if (xbar .gt. x(2)) go to 101
      i = 1
      return
 101  continue
      if (xbar .lt. x(n-1)) go to 102
      i = n-1
      return
 102  continue
c..................................................................
c     This version knows data is evenly spaced with spacing dx
c..................................................................

      i=(xbar-x(1))/dx+1

      return
      end


      subroutine terp1e (n,x,f,w,y,int,tab,itab,dx)
      implicit integer (i-n), real*8 (a-h,o-z)
cSm030221
c      dimension       x(n)       ,f(n*int)   ,w(n*int)    ,tab(3)     ,
c     1  itab(3)
      dimension       x(*)       ,f(*)   ,w(*)    ,tab(3)     ,
     1  itab(3)

c     the following call is for gathering statistics at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     perform search
c
      call searche (y,x,n,i,dx)
c
c     carry out interpolation (or extrapolation)
c
      call intrp (n,x,f,w,y,i,int,tab,itab)
      return
      end


      real*8 function t2
     +     (dx,dy,xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd)
      implicit integer (i-n), real*8 (a-h,o-z)
cSm030221
c      dimension       x(nx)     ,y(ny)     ,f(idm,ny)  ,fxx(idm,ny) ,
c     1  fyy(idm,ny)  ,fxxyy(idm,ny)        ,ff(2)      ,
c     2  ww(2)      ,tab(3)     ,itab(3)
      dimension       x(*)     ,y(*)     ,f(idm,*)  ,fxx(idm,*) ,
     1  fyy(idm,*)  ,fxxyy(idm,*)        ,ff(2)      ,
     2  ww(2)      ,tab(3)     ,itab(3)
c     the following call is for gathering statistics 
c       on library use at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     search in x and y arrays.
c
      call searche (xb,x,nx,i,dx)
      call searche (yb,y,ny,j,dy)
c
c     interpolate in x direction
c
      do 101 i1=1,3
        itab(i1) = 0
 101  continue
      i1 = ixd+1
      itab(i1) = 1
      do 102 j1=1,2
        jj = j+j1-1
        call intrp (n,x,f(1,jj),fxx(1,jj),xb,i,1,tab,itab)
        ff(j1) = tab(i1)
        call intrp (n,x,fyy(1,jj),fxxyy(1,jj),xb,i,1,tab,itab)
        ww(j1) = tab(i1)
 102  continue
c
c     interpolate in y direction
c
      do 103 j1=1,3
        itab(j1) = 0
 103  continue
      j1 = iyd+1
      itab(j1) = 1
      call intrp (2,y(j),ff,ww,yb,1,1,tab,itab)
      t2 = tab(j1)
      return
c
c     revision history---
c
c     june 1977        replaced non-standard statement functions and
c     subscripts to enhance portability.
c
c     january 1978     deleted references to the  *cosy  cards, moved
c     the revision histories to appear before the
c     final end card, and moved the initial comment
c     cards to appear after the first subroutine card
c     and changed  itab  from logical to integer in
c     subroutine intrp and corrected problem with
c     version numbers in one statistics call
c-----------------------------------------------------------------------
      end


      subroutine trip (n,a,b,c,y,z,int)
      implicit integer (i-n), real*8 (a-h,o-z)
cSm030221
c      dimension       a(n)       ,b(n)       ,c(n)       ,y(n)       ,
c     1  z(n*int)
      dimension       a(*)       ,b(*)       ,c(*)       ,y(*)       ,
     1  z(*)
c
c     arithmetic statement function used to locate entries in array z.
c
      ii(index)=(index-1)*int+1
c
c     gaussian elimination
c
      bn = b(n)
      yn = y(n)
      v = c(n)
      y(1) = y(1)/b(1)
      a(1) = a(1)/b(1)
      b(1) = c(1)/b(1)
      nm2 = n-2
      do 101 j=2,nm2
        den = b(j)-a(j)*b(j-1)
        b(j) = c(j)/den
        y(j) = (y(j)-a(j)*y(j-1))/den
        a(j) = -a(j)*a(j-1)/den
        bn = bn-v*a(j-1)
        yn = yn-v*y(j-1)
        v = -v*b(j-1)
 101  continue
      den = b(n-1)-a(n-1)*b(n-2)
      b(n-1) = (c(n-1)-a(n-1)*a(n-2))/den
      y(n-1) = (y(n-1)-a(n-1)*y(n-2))/den
      bn = bn-v*a(n-2)
      yn = yn-v*y(n-2)
      v = a(n)-v*b(n-2)
c     back substitution
      iin = ii(n)
      z(iin) = (yn-v*y(n-1))/(bn-v*b(n-1))
      iin2 = ii(n-1)
      z(iin2) = y(n-1)-b(n-1)*z(iin)
      nm1 = n-1
      in = ii(n)
      do 102 j=2,nm1
        k = n-j
        ik = ii(k)
        ikt = ik+int
 102  z(ik) = y(k)-b(k)*z(ikt)-a(k)*z(in)
      return
      end
c***********************************************************************
c     END OF SPLINES
c***********************************************************************
c
c
      real*8 function erf(xxx)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c------------------------------------------------------
c     This routine computes the ERROR FUNCTION.
c------------------------------------------------------
      common /tusq/ tusqpi
      dimension a(5)
      sign=1.
      if (xxx .lt. 0.) sign=-1.
      xcg=sign*xxx
      x2=xcg*xcg
      if (xcg .ge. .6) go to 20
      sum=xcg
      term=xcg
      kmax=6
      do 10 k=1,kmax
        t1=float(k)
        t2=float(2*k+1)/float(2*k-1)
        term=-term*x2/(t1*t2)
        sum=sum+term
 10   continue
      erf=tusqpi*sum
      erf=sign*erf
      return
 20   continue
      p=.3275911
      a(1)=.225836846
      a(2)=-.252128668
      a(3)=1.25969513
      a(4)=-1.287822453
      a(5)=.94064607
      eta=1./(1.+p*xcg)
      phip=tusqpi*exp(-x2)
      term=(((a(5)*eta+a(4))*eta+a(3))*eta+a(2))*eta+a(1)
      erf=1.-term*eta*phip
      erf=sign*erf
      return
      end
c######date01jan1984     copyright ukaea, harwell.
c######aliasim01ad


      real*8 function im01ad(la,a,inc)
      implicit integer (i-n), real*8 (a-h,o-z)
      real*8 a(la*inc)
      kount = 0
      do 100 k=1,la
        ipos = (k-1)*inc + 1
        if (a(ipos).eq.0.0d0) kount = kount + 1
 100  continue
      im01ad = kount
      return
      end


      double precision function terp2p
     + (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd
     1                 ,isrch)
c
c Modified terp2 by adding isrch:  Binary search of grid for
c                    interpolation point is carried out only if
c                    isrch=1.  Otherwise, it is assumed that values
c                    have been generated by a previous call.
c
c
c dimension of           x(nx),y(ny),f(idm,ny),fxx(idm,ny),fyy(idm,ny),
c arguments              fxxyy(idm,ny))
c                        (idm must be .ge. nx)
c
c latest revision        february 1974
c
c usage                  r = terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,
c                                   ixd,iyd)
c
c arguments
c
c on input               xb, yb
c                          values of the independent variables, x and y,
c                          at which interpolation is required.
c
c                        nx
c                          number of grid points in the x-direction.  nx
c                          must be at least 4.
c
c                        x
c                          table of nx values of the independent
c                          variable, x, arranged in ascending order.
c                          dimension of x in the calling program must be
c                          at least nx.
c
c                        ny
c                          number of grid points in the y-direction.  ny
c                          must be at least 4.
c
c                        y
c                          table of ny values of the independent
c                          variable, y, arranged in ascending order.
c                          dimension of y in the calling program must be
c                          at least ny.
c
c                        f
c                          two-dimensional array of function values at
c                          grid points defined by the arrays x and y.
c                          dimension of f in the calling program is
c                          (idm,nyy), where
c                              idm .ge. nx
c                              nyy .ge. ny
c
c                        fxx
c                          array of second derivatives of f with respect
c                          to x computed by coeff2.  dimension of fxx in
c                          the calling program is (idm,nyy).  see under
c                          f above.
c
c                        fyy
c                          array of second derivatives of f with respect
c                          to y computed by coeff2.  dimension of fyy in
c                          the calling program is (idm,nyy).
c
c                        fxxyy
c                          array of fourth derivatives,
c                          (d/dx)**2*(d/dy)**2*f, computed by coeff2.
c                          dimension of fxxyy in the calling program is
c                          (idm,nyy).
c
c                        idm
c                          first dimension in the calling program of
c                          arrays f, fxx, fyy and fxxyy,
c                              idm .ge. nx
c
c                        ixd, iyd
c                          define derivative to be returned by the
c                          function terp2.  ixd, iyd may each take the
c                          the values 0, 1, 2.  the derivative returned
c                          is (d/dx)**ixd*(d/dy)**iyd*f.
c                            note that if ixd = iyd = 0, the function
c                            value itself is returned.
c
c space required         220 (octal) = 144 (decimal)
c
c timing                 this procedure is fast.  the maximum time for
c                        the binary search is proportional to
c                        alog(nx*ny).  the time for function evaluation
c                        is independent of n.  for a 21 x 21 grid of
c                        data points, an average time for an
c                        interpolation on the ncar cdc 7600 is about
c                        .29 milliseconds.
c
c
c
c
      !implicit double precision (a-h,o-z)
      implicit none
      integer i,j,i1,j1,jj,n,ny,nx,isrch,ixd,iyd,idm,isav,jsav,itab
      real*8 x,xb,y,yb,f,ff,fxx,fyy,fxxyy,ww,tab
cSm030221
c      dimension       x(1)       ,y(1)       ,f(idm,1)   ,fxx(idm,1) ,
c     1                fyy(idm,1) ,fxxyy(idm,1)           ,ff(2)      ,
c     2                ww(2)      ,tab(3)     ,itab(3)
      dimension       x(*)       ,y(*)       ,f(idm,*)   ,fxx(idm,*) ,
     1                fyy(idm,*) ,fxxyy(idm,*)           ,ff(2)      ,
     2                ww(2)      ,tab(3)     ,itab(3)
c the following call is for gathering statistics at ncar
c
c search in x and y arrays.
c
      SAVE

      if(isrch.eq.1)  then
        call search (xb,x,nx,i)
        call search (yb,y,ny,j)
        isav=i
        jsav=j
      endif

c
c interpolate in x direction
c
      do 101 i1=1,3
         itab(i1) = 0
  101 continue
      i1 = ixd+1
      itab(i1) = 1
      do 102 j1=1,2
         jj = jsav+j1-1
c         write(*,*)'jj,j1,jsave,fxx(1,jj)',jj,j1,jsave,fxx(1,jj)
         call intrp (n,x,f(1,jj),fxx(1,jj),xb,isav,1,tab,itab)
         ff(j1) = tab(i1)
c         write(*,*)'102 i1,j1 ff(j1)',i1,j1, ff(j1)
         call intrp (n,x,fyy(1,jj),fxxyy(1,jj),xb,isav,1,tab,itab)
         ww(j1) = tab(i1)
  102 continue
c
c interpolate in y direction
c
      do 103 j1=1,3
         itab(j1) = 0
  103 continue
      j1 = iyd+1
      itab(j1) = 1
      call intrp (2,y(jsav),ff,ww,yb,1,1,tab,itab)
      terp2p = tab(j1)

c      write(*,*)'xb,yb,terp2p',xb,yb,terp2p

      return
      end



c---------------------------------------------------------------
c     modified spline functions to use maximal dimensions

      subroutine coeff2_Sm(nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ibd,wk,nxa)
cSm030221
c     nxa is the maximal value of nx.
c     nx is dimension of arrays f,fxx,fyy,fxxyy(nxa,nya) and wk
c--------------------------------------------------------
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      integer i,j,ny,nx,nxa, iloc,jloc, ibd,idm
      real*8 x,y,f,fxx,fyy,fxxyy,wk
c
cSm030220
c      dimension       x(nx)       ,y(ny)       ,f(idm,ny)  ,fxx(idm,ny),
c     1  fyy(idm,ny) ,fxxyy(idm,ny)           ,ibd(4)     ,
c     2  iloc(2)    ,jloc(2)
      dimension       x(*)       ,y(*)       ,f(nxa,*)  ,fxx(nxa,*),
     1  fyy(nxa,*) ,fxxyy(nxa,*)           ,ibd(4)     ,
     2  iloc(2)    ,jloc(2),
cSAP091201
     &  wk(*)

      logical q8q4
      save q8q4
      data q8q4 /.true./
      data iloc(1),iloc(2),jloc(1),jloc(2)/1,1,4,4/
c     the following call is for gathering statistics at ncar
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     compute fxx
c
      do 101 j=1,ny

cBH020822
cBH020822c        call coeff1 (nx,x,f(1,j),fxx(1,j),ibd(1),1,wk)
         call coeff1 (nx,x,f(1,j),fxx(1,j),ibd(1),1,wk)
cBH020822
c        write(*,*)'coeff2 fxx(1,1)',fxx(1,1)

c
 101  continue
    
c
c     compute fyy
c
      do 102 i=1,nx
        call coeff1 (ny,y,f(i,1),fyy(i,1),ibd(3),idm,wk)
 102  continue
c
c     check for periodic boundary condition in both directions
c
      if (ibd(1) .eq. 3) go to 103
      if (ibd(3) .eq. 3) go to 105
c
c     calculate fxxyy along left and right boundaries
c
      call coeff1 (ny,y,fxx(1,1),fxxyy(1,1),jloc,idm,wk)
      call coeff1 (ny,y,fxx(nx,1),fxxyy(nx,1),jloc,idm,wk)
      go to 106
 103  continue
c
c     periodic in x direction . calculate fxxyy along lower and upper
c     boundaries.
c
      call coeff1 (nx,x,fyy(1,1),fxxyy(1,1),ibd(1),1,wk)
      call coeff1 (nx,x,fyy(1,ny),fxxyy(1,ny),ibd(1),1,wk)
c
c     calculate remaining fxxyy
c
      do 104 i=1,nx
        call coeff1 (ny,y,fxx(i,1),fxxyy(i,1),iloc,idm,wk)
 104  continue
      go to 108
 105  continue
c
c     periodic in y direction. calculate fxxyy along left and right
c     boundaries.
c
      call coeff1 (ny,y,fxx(1,1),fxxyy(1,1),ibd(3),idm,wk)
      call coeff1 (ny,y,fxx(nx,1),fxxyy(nx,1),ibd(3),idm,wk)
 106  continue
c
c     calculate remaining fxxyy
c
      do 107 j=1,ny
        call coeff1 (nx,x,fyy(1,j),fxxyy(1,j),iloc,1,wk)
 107  continue
 108  continue
      return
      end


      real*8 function terp2_Sm 
     +     (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd,nxa)
cSm030221
c     nxa is the maximal value of nx.
c     nx is dimension of arrays f,fxx,fyy,fxxyy(nxa,nya)
c----------------------------------------------------------
      implicit integer (i-n), real*8 (a-h,o-z)
c
cSm030220
c      dimension       x(nx)      ,y(ny)      ,f(idm,ny)  ,fxx(idm,ny),
c     1  fyy(idm,ny) ,fxxyy(idm,ny)           ,ff(2)      ,
c     2  ww(2)      ,tab(3)     ,itab(3)
      dimension       x(*)      ,y(*)      ,f(nxa,*)  ,fxx(nxa,*),
     1  fyy(nxa,*) ,fxxyy(nxa,*)           ,ff(2)      ,
     2  ww(2)      ,tab(3)     ,itab(3)
c     the following call is for gathering statistics 
c       on library use at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     search in x and y arrays.
c
      call search (xb,x,nx,i)
      call search (yb,y,ny,j)
c
c     interpolate in x direction
c
      do 101 i1=1,3
        itab(i1) = 0
 101  continue
      i1 = ixd+1
      itab(i1) = 1
      do 102 j1=1,2
        jj = j+j1-1
        call intrp (n,x,f(1,jj),fxx(1,jj),xb,i,1,tab,itab)
        ff(j1) = tab(i1)
        call intrp (n,x,fyy(1,jj),fxxyy(1,jj),xb,i,1,tab,itab)
        ww(j1) = tab(i1)
 102  continue
c
c     interpolate in y direction
c
      do 103 j1=1,3
        itab(j1) = 0
 103  continue
      j1 = iyd+1
      itab(j1) = 1
      call intrp (2,y(j),ff,ww,yb,1,1,tab,itab)
      terp2_Sm = tab(j1)
      return
c
c     revision history---
c
c     june 1977        replaced non-standard statement functions and
c     subscripts to enhance portability.
c
c     january 1978     deleted references to the  *cosy  cards, moved
c     the revision histories to appear before the
c     final end card, and moved the initial comment
c     cards to appear after the first subroutine card
c     and changed  itab  from logical to integer in
c     subroutine intrp and corrected problem with
c     version numbers in one statistics call
c-----------------------------------------------------------------------
      end

      real*8 function t2_Sm
     +     (dx,dy,xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd,nxa)
cSm030221
c     nxa is the maximal value of nx.
c     nx is dimension of arrays f,fxx,fyy,fxxyy(nxa,nya)
c-------------------------------------------------------
      implicit integer (i-n), real*8 (a-h,o-z)
c
cSmo30220
c      dimension       x(nx)     ,y(ny)     ,f(idm,ny)  ,fxx(idm,ny) ,
c     1  fyy(idm,ny)  ,fxxyy(idm,ny)        ,ff(2)      ,
c     2  ww(2)      ,tab(3)     ,itab(3)
      dimension       x(*)     ,y(*)     ,f(nxa,*)  ,fxx(nxa,*) ,
     1  fyy(nxa,*)  ,fxxyy(nxa,*)        ,ff(2)      ,
     2  ww(2)      ,tab(3)     ,itab(3)
c     the following call is for gathering statistics 
c       on library use at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     search in x and y arrays.
c
      call searche (xb,x,nx,i,dx)
      call searche (yb,y,ny,j,dy)
c
c     interpolate in x direction
c
      do 101 i1=1,3
        itab(i1) = 0
 101  continue
      i1 = ixd+1
      itab(i1) = 1
      do 102 j1=1,2
        jj = j+j1-1
        call intrp (n,x,f(1,jj),fxx(1,jj),xb,i,1,tab,itab)
        ff(j1) = tab(i1)
        call intrp (n,x,fyy(1,jj),fxxyy(1,jj),xb,i,1,tab,itab)
        ww(j1) = tab(i1)
 102  continue
c
c     interpolate in y direction
c
      do 103 j1=1,3
        itab(j1) = 0
 103  continue
      j1 = iyd+1
      itab(j1) = 1
      call intrp (2,y(j),ff,ww,yb,1,1,tab,itab)
      t2_Sm = tab(j1)
      return
c
c     revision history---
c
c     june 1977        replaced non-standard statement functions and
c     subscripts to enhance portability.
c
c     january 1978     deleted references to the  *cosy  cards, moved
c     the revision histories to appear before the
c     final end card, and moved the initial comment
c     cards to appear after the first subroutine card
c     and changed  itab  from logical to integer in
c     subroutine intrp and corrected problem with
c     version numbers in one statistics call
c-----------------------------------------------------------------------
      end




      double precision function terp2p_Sm
     + (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd
     1                 ,isrch,nxa)
cSm030221
c     nxa is the maximal value of nx.
c     nx is dimension of arrays f,fxx,fyy,fxxyy(nxa,nya)
c----------------------------------------------------------
c Modified terp2 by adding isrch:  Binary search of grid for
c                    interpolation point is carried out only if
c                    isrch=1.  Otherwise, it is assumed that values
c                    have been generated by a previous call.
c
c
c dimension of           x(nx),y(ny),f(idm,ny),fxx(idm,ny),fyy(idm,ny),
c arguments              fxxyy(idm,ny))
c                        (idm must be .ge. nx)
c
c latest revision        february 1974
c
c usage                  r = terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,
c                                   ixd,iyd)
c
c arguments
c
c on input               xb, yb
c                          values of the independent variables, x and y,
c                          at which interpolation is required.
c
c                        nx
c                          number of grid points in the x-direction.  nx
c                          must be at least 4.
c
c                        x
c                          table of nx values of the independent
c                          variable, x, arranged in ascending order.
c                          dimension of x in the calling program must be
c                          at least nx.
c
c                        ny
c                          number of grid points in the y-direction.  ny
c                          must be at least 4.
c
c                        y
c                          table of ny values of the independent
c                          variable, y, arranged in ascending order.
c                          dimension of y in the calling program must be
c                          at least ny.
c
c                        f
c                          two-dimensional array of function values at
c                          grid points defined by the arrays x and y.
c                          dimension of f in the calling program is
c                          (idm,nyy), where
c                              idm .ge. nx
c                              nyy .ge. ny
c
c                        fxx
c                          array of second derivatives of f with respect
c                          to x computed by coeff2.  dimension of fxx in
c                          the calling program is (idm,nyy).  see under
c                          f above.
c
c                        fyy
c                          array of second derivatives of f with respect
c                          to y computed by coeff2.  dimension of fyy in
c                          the calling program is (idm,nyy).
c
c                        fxxyy
c                          array of fourth derivatives,
c                          (d/dx)**2*(d/dy)**2*f, computed by coeff2.
c                          dimension of fxxyy in the calling program is
c                          (idm,nyy).
c
c                        idm
c                          first dimension in the calling program of
c                          arrays f, fxx, fyy and fxxyy,
c                              idm .ge. nx
c
c                        ixd, iyd
c                          define derivative to be returned by the
c                          function terp2.  ixd, iyd may each take the
c                          the values 0, 1, 2.  the derivative returned
c                          is (d/dx)**ixd*(d/dy)**iyd*f.
c                            note that if ixd = iyd = 0, the function
c                            value itself is returned.
c
c space required         220 (octal) = 144 (decimal)
c
c timing                 this procedure is fast.  the maximum time for
c                        the binary search is proportional to
c                        alog(nx*ny).  the time for function evaluation
c                        is independent of n.  for a 21 x 21 grid of
c                        data points, an average time for an
c                        interpolation on the ncar cdc 7600 is about
c                        .29 milliseconds.
c
c
c
c
      !implicit double precision (a-h,o-z)
      implicit none
      integer i,j,i1,j1,jj,n,ny,nx,nxa,isrch,ixd,iyd,idm,isav,jsav,itab
      real*8 x,xb,y,yb,f,ff,fxx,fyy,fxxyy,ww,tab
cSm030220
c      dimension       x(1)       ,y(1)       ,f(idm,1)   ,fxx(idm,1) ,
c     1                fyy(idm,1) ,fxxyy(idm,1)           ,ff(2)      ,
c     2                ww(2)      ,tab(3)     ,itab(3)
      dimension       x(*)       ,y(*)       ,f(nxa,*)   ,fxx(nxa,*) ,
     1                fyy(nxa,*) ,fxxyy(nxa,*)           ,ff(2)      ,
     2                ww(2)      ,tab(3)     ,itab(3)
c the following call is for gathering statistics at ncar
c
c search in x and y arrays.
c
      SAVE

      if(isrch.eq.1)  then
        call search (xb,x,nx,i)
        call search (yb,y,ny,j)
        isav=i
        jsav=j
      endif

c
c interpolate in x direction
c
      do 101 i1=1,3
         itab(i1) = 0
  101 continue
      i1 = ixd+1
      itab(i1) = 1
      do 102 j1=1,2
         jj = jsav+j1-1
c         write(*,*)'jj,j1,jsave,fxx(1,jj)',jj,j1,jsave,fxx(1,jj)
              
         call intrp (n,x,f(1,jj),fxx(1,jj),xb,isav,1,tab,itab)
         ff(j1) = tab(i1)
c         write(*,*)'102 i1,j1 ff(j1)',i1,j1, ff(j1)
         call intrp (n,x,fyy(1,jj),fxxyy(1,jj),xb,isav,1,tab,itab)
         ww(j1) = tab(i1)
  102 continue
c
c interpolate in y direction
c
      do 103 j1=1,3
         itab(j1) = 0
  103 continue
      j1 = iyd+1
      itab(j1) = 1
      call intrp (2,y(jsav),ff,ww,yb,1,1,tab,itab)
      terp2p_Sm = tab(j1)

c      write(*,*)'xb,yb,terp2p',xb,yb,terp2p

      return
      end
c
      subroutine lookup(x,xarray,length,weightu,weightl,lement)
      !implicit integer (i-n), double precision (a-h,o-z)
      implicit none !YuP[2020-01-14]
      integer length,lement,luf  !YuP[2020-01-14]
      real*8 x,xarray,weightu,weightl  !YuP[2020-01-14]
c..................................................................
c     This routine uses luf to do a table look up. Then it interpolates
c     to gain a bit of accuracy.
c     x is the argument; xarray is the monotonic array; length
c     is the length of the array. lement is the first index such
c     that xarray(lement).gt.x.
c     If x falls outside bounds of the array, weightl/weightr set
c     to give constant value at relevant bounds of the array.
c..................................................................

      save
      dimension xarray(*) 

c      write(*,*)'in lookup before luf x,length',x,length

      lement=luf(x,xarray,length)

cSAP080713BH080714
c      write(*,*)'in lookup after luf lement',lement
      if(lement.lt.2) then 
         write(*,*)'WARNING: in lookup lement.lt.2'
         write(*,*)'the code will set lement=2'       
         lement=2
         weightl=1.d0
         weightu=0.d0
         goto 10
      endif

      if(lement.gt.length) then 
         write(*,*)'WARNING: in lookup lement.gt.length'
         write(*,*)'the code will set lement=lenthg'       
         lement=length
         weightl=0.d0
         weightu=1.d0
         goto 10
      endif

      weightl=(xarray(lement)-x)/(xarray(lement)-xarray(lement-1))
      weightu=1.d0-weightl

 10   continue

      return
      end


      integer function luf(x,table,n)
      !implicit integer (i-n), double precision (a-h,o-z)
      implicit none !YuP[2020-01-14]
      integer n,i_c,i_l,i_r  !YuP[2020-01-14]
      real*8 x,table  !YuP[2020-01-14]
c
c     THIS ROUTINE SHOULD BE A BINARY SEARCH.  IT NEEDS
C        WORK!
c     luf(x,table,n) (MATHLIB) which is a function returning the index
c        of the first element in the table that is greater than x.
c        Elements must be strictly increasing. x.gt.table(n)==>n+1.
c
c     find first index such that table(luf).gt.x
c
c
      dimension table(n)
c
c      do i=1,n
c         if (table(i) .gt. x) go to 10
c      end do

c 10   continue
c      luf = 1 if x.lt.table(1) and luf=n+1 if x>ge.table(n)
c      luf = i
      

cSm011227 a binary search------------------------
c       luf_line=luf
       if(x.ge.table(n))then
         i_c=n+1
         goto 30
       endif

       if(x.lt.table(1))then
         i_c=0
         goto 30
       endif

      i_l=1 !initialization
      i_r=n !initialization
 20   continue
      if ((i_r-i_l).eq.1) then
         i_c=i_r
         goto 30
      else
         i_c=0.5*(i_l+i_r) !the central point
      endif

      if (x.lt.table(i_c)) then
         i_r=i_c
      else
         i_l=i_c
      endif

      goto 20      
 30   continue
      luf = i_c
       
c      if (luf.ne.luf_line) then
c      write(*,*)'binary luf,luf_line',luf,luf_line
c      write(*,*)'n,x,table(1),table(n)',n,x,table(1),table(n)
c      endif
c------------------------------------------

c
      return
      end


      subroutine lin_interp(x,f,ix,xx,ff,ixx)
      !implicit integer (i-n), double precision (a-h,o-z)
      implicit none !YuP[2020-01-14]
      integer ix,ixx,ii,lx  !YuP[2020-01-14]
      real*8 x,f,xx,ff  !YuP[2020-01-14]
      real*8 weightu,weightl

c this routine linearly interpolates the values of array f on x
c to the values ff on xx. (BobH, 050522).

c Input:  x(1:ix),f(1:ix),ix, where f is f(x)
c         xx(1:ixx),ixx

c If xx(i).le.x(1), then ff(i)=f(1)   [for any i]
c If xx(i).ge.x(ix), then ff(i)=f(ix) [for any i]

c Uses subroutines lookup and luf

      dimension x(ix),f(ix),xx(ixx),ff(ixx)


c      write(*,*)'zcunix.f lin_interp ix,ixx',ix,ixx

      do ii=1,ixx

cSAP080713
c         write(*,*)'lin_interp bef lookup ii,ixx',ii,ixx
c         if(ii.lt.ixx) then
c             write(*,*)'xx(ii),x(ii+1)',xx(ii),x(ii+1)
c         endif

c         write(*,*)'before lookup xx(ii)',xx(ii)

         call lookup(xx(ii),x,ix,weightu,weightl,lx)

c         write(*,*)'after lookup:ix,weightu,weightl,lx',
c     &                           ix,weightu,weightl,lx

c$$$         if (lx.gt.ix) then
c$$$            lx=ix
c$$$            weightl=0.d0
c$$$            weightu=1.d0
c$$$         elseif (lx.le.1) then
c$$$            lx=2
c$$$            weightl=1.d0
c$$$            weightu=0.d0
c$$$         endif
         
         ff(ii)=weightl*f(lx-1)+weightu*f(lx)
         
      enddo

c      write(*,*)'lin_interp:ix,x',ix,x
c      write(*,*)'lin_interp:,f',f
c      write(*,*)'lin_interp:ixx,xx',ixx,xx
c      write(*,*)'lin_interp:,ff',ff
      
      return

      end


c!!!!!!!!!!!!!!!test
      subroutine terp1_Sm (n,x,f,w,y,int,tab,itab)
      implicit integer (i-n), real*8 (a-h,o-z)
cSm030221
      dimension       x(n)       ,f(n*int)       ,w(n*int)    ,tab(3),
     1  itab(3)
c      dimension       x(*)       ,f(*)       ,w(*)    ,tab(3),
c     1  itab(3)

c     the following call is for gathering statistics at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
c
c     perform search
c
      call search (y,x,n,i)
c
c     carry out interpolation (or extrapolation)
c

c      call intrp (n,x,f,w,y,i,int,tab,itab)
      call intrp_Sm (n,x,f,w,y,i,int,tab,itab)

      return
      end

      subroutine intrp_Sm (n,x,f,w,y,i,int,tab,itab)
      implicit integer (i-n), real*8 (a-h,o-z)
cSm030221
c      dimension       x(i+1)    ,f(i*int+1)    ,w(i*int+1)  ,tab(3)
c     -  ,itab(3)
      dimension       x(*)    ,f(*)    ,w(*)  ,tab(3)
     -  ,itab(3)

c
c     arithmetic statement function used to locate entries in 
c     f and w arrays:
c
      ii(index)=(index-1)*int+1
c
c     perform interpolation or extrapolation
c
      flk = x(i+1)-x(i)
      flp = x(i+1)-y
      fl0 = y-x(i)

c      write(*,*)'intrp flk,flp,fl0', flk,flp,fl0

      i0 = ii(i)
      ip = i0+int

c      write(*,*)'i0,ip,itab(1)',i0,ip,itab(1)

      if (itab(1) .ne. 0) go to 101
      go to 102
 101  continue
c
c     calculate f(y)
c
      a = (w(i0)*flp**3+w(ip)*fl0**3)/(6.*flk)

c      write(*,*)'i0,ip,w(i0),w(ip),flp,fl0,flk,a',
c     +i0,ip,w(i0),w(ip),flp,fl0,flk,a

      b = (f(ip)/flk-w(ip)*flk/6.)*fl0
      c = (f(i0)/flk-w(i0)*flk/6.)*flp
      tab(1) = a+b+c

c      write(*,*)'101 a,b,c tab(1)',a,b,c,tab(1)

 102  continue
      if (itab(2) .ne. 0) go to 103
      go to 104
 103  continue
c
c     calculate first derivative at y
c
      a = (w(ip)*fl0**2-w(i0)*flp**2)/(2.*flk)
      b = (f(ip)-f(i0))/flk
      c = (w(i0)-w(ip))*flk/6.
      tab(2) = a+b+c
c      write(*,*)'103 a,b,c,tab(2)',a,b,c,tab(2)
 104  continue
      if (itab(3) .ne. 0) go to 105
      go to 106
 105  continue
c
c     calculate second derivative at y
c
      tab(3) = (w(i0)*flp+w(ip)*fl0)/flk
c      write(*,*)'105 tab(3),flk',tab(3),flk
 106  continue
      return
      end
