c
c
      subroutine bcast(a,val,n)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     Temporary bcast routine until I can find UNICOS equivalent
c..................................................................

      dimension a(n)
      do 100 i=1,n
        a(i)=val
 100  continue
      return
      end
c
c
      subroutine ibcast(ia,ival,n)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     Temporary bcast routine until I can find UNICOS equivalent
c..................................................................

      dimension ia(n)
      do 100 i=1,n
        ia(i)=ival
 100  continue
      return
      end

c NME bcast routine for complex arrays
      subroutine ccast(c,cval,n)
      implicit integer (i-n), complex*16 (c)
      dimension c(n)
      do 100 i=1,n
         c(i)=cval
 100  continue
      return
      end

      subroutine r4bcast(a,val,n)
      implicit integer (i-n), real*4 (a-h,o-z)

c..................................................................
c     Temporary bcast routine until I can find UNICOS equivalent
c..................................................................

      dimension a(n)
      do 100 i=1,n
        a(i)=val
 100  continue
      return
      end
c
