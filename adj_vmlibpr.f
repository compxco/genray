

      subroutine mcsevu(h, x1, y, ny, yc, ic, ns, s, xs, ier)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  11:00:56 6/14/01  
C...Switches: -p4 -yb
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ny, ic, ns, ier
      real h, x1
      real, dimension(ny) :: y
      real, dimension(ic,3) :: yc
      real, dimension(ns) :: s, xs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(1) :: is
      integer :: i
      real, dimension(ns) :: d
      real :: zero
C-----------------------------------------------
c--------------------
c  evaluation of cubic spline
c   given y(ny),yc(ny,3) at x(j)=x1+(j-1)*h   j=1...ny
c   evaluate s(1:ns) at xs(1:ns)
c    xs(1:ns) should be less than x(ny)
c----------------------------------
      data zero/0.0/
      is(:ns) = (xs - x1)/h
      is(:ns) = min0(ny - 2,is(:ns))
      d = xs - is(:ns)*h - x1
      is(:ns) = is(:ns) + 1
 
c     do i=1,ns
c       wk(i)=yc(is(i),3)
c     enddo
c     do i=1,ns
c       s(i)=wk(i)*d(i)
c     enddo
c     do i=1,ns
c       wk(i)=yc(is(i),2)
c     enddo
c     do i=1,ns
c       s(i)=(s(i)+wk(i))*d(i)
c     enddo
c     do i=1,ns
c       wk(i)=yc(is(i),1)
c     enddo
c     do i=1,ns
c       s(i)=(s(i)+wk(i))*d(i)
c     enddo
c     do i=1,ns
c       wk(i)=y(is(i))
c     enddo
c     do i=1,ns
c       s(i)=s(i)+wk(i)
c     enddo
      s = ((yc(is(:ns),3)*d+yc(is(:ns),2))*d+yc(is(:ns),1))*d+
     1   y(is(:ns))
      return 
      end subroutine mcsevu

      subroutine mcsevu1(h, x1, y, ny, yc, ic, ns, s, xs, ier)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  11:00:56 6/14/01  
C...Switches: -p4 -yb
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ny, ic, ns, ier
      real h, x1
      real, dimension(ny) :: y
      real, dimension(ic,3) :: yc
      real, dimension(ns) :: s, xs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(ns) :: is
      integer :: i
      real, dimension(ns) :: d, wk
C-----------------------------------------------
c---------------------------
c   given y(1:ny),yc(1:ny,1:3) at x(i)=x1+(i-1)*h
c   find s(1:ns)=first derivative of function at xs(1:ns)
c-----------------------------
      is = (xs - x1)/h
      is = min0(ny - 2,is)
      d = xs - (is*h + x1)
      is = is + 1
      wk = yc(is,3)
      s = 3.*wk*d
      wk = yc(is,2)
      s = (s + 2.*wk)*d
      wk = yc(is,1)
      s = s + wk
      return 
      end subroutine mcsevu1

