!
!
      subroutine bcast(a,val,n)
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none !YuP[2020-01]
      real*8 val,a  !YuP[2020-01]
      integer i,n   !YuP[2020-01]
!..................................................................
!     Temporary bcast routine until I can find UNICOS equivalent
!..................................................................
      dimension a(n)
      do 100 i=1,n
        a(i)=val
 100  continue
      return
      end
!
!
      subroutine ibcast(ia,ival,n)
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none !YuP[2020-01]
      integer i,n,ival,ia   !YuP[2020-01]
!..................................................................
!     Temporary bcast routine until I can find UNICOS equivalent
!..................................................................
      dimension ia(n)
      do 100 i=1,n
        ia(i)=ival
 100  continue
      return
      end

! NME bcast routine for complex arrays
      subroutine ccast(c,cval,n)
      !implicit integer (i-n), complex*16 (c)
      implicit none !YuP[2020-01]
      complex*16 cval,c  !YuP[2020-01]
      integer i,n   !YuP[2020-01]
      dimension c(n)
      do 100 i=1,n
         c(i)=cval
 100  continue
      return
      end

      subroutine r4bcast(a,val,n)
      !implicit integer (i-n), real*4 (a-h,o-z)
      implicit none !YuP[2020-01]
      real*4 val,a  !YuP[2020-01] ![2020-09-05] Corrected to real*4
      integer i,n   !YuP[2020-01]
!..................................................................
!     Temporary bcast routine until I can find UNICOS equivalent
!..................................................................
      dimension a(n)
      do 100 i=1,n
        a(i)=val
 100  continue
      return
      end
!
