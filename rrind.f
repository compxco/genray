c
c
c
c      double precision function rrind(npar,nper,omega,omegpe,nsigma) 
cSAP100921
      double precision function rrind(npar,nper,omega,omegpe)
      implicit double precision (a-h,o-z)
      double precision npar,npar2,nper,nper2,n2
       
      write(*,*)' function rrind: npar,nper,omega,omegpe',
     &                            npar,nper,omega,omegpe
c
c     calculate ray refractive index
c
      eps1=1.d0-omegpe*omegpe/(omega*omega-1.d0)
      eps2=-omegpe*omegpe/(omega*(omega*omega-1.d0))
      eps3=1.d0-omegpe*omegpe/(omega*omega)
c
      npar2=npar*npar
c
      nper2=nper*nper
      n2=nper2+npar2

      write(*,*)' function rrind 1'
c
c     calculate ray refractive index squared
c
      rn=dsqrt(n2)
      st=nper/rn
      ct=npar/rn
      st2=st*st
      ct2=ct*ct
      var1=eps1-n2
      var2=eps3-nper2
      var3=eps1-npar2
      eps22=eps2*eps2
      dnumer=-st*ct*var1*var2+st*ct*var3*var1-eps22*st*ct
     1       +n2*var1*st*ct*(ct2-st2)
      ddenom=eps22*st2+n2*n2*st2*ct2-var1*var2*ct2-var3*var2
     1       -st2*var3*var1-2.0*n2*st2*ct2*var1
      dndt=rn*dnumer/ddenom
      droot=dsqrt(1.d0+(dndt/rn)*(dndt/rn))

      write(*,*)' function rrind 3'

      deriv=( var1*var2*(st2-ct2)+2.d0*rn*dndt*var2*st*ct
     1       +2.d0*rn*dndt*st2*st*ct*var1+2.0*st2*ct2*var1*n2
     2       +(ct2-st2)*var3*var1-
     1       2.d0*rn*dndt*st*ct*(ct2*var1+var3)
     3       +2.d0*st2*ct2*n2*var1
     4       +eps22*(st2-ct2)
     5       +2.d0*rn*dndt*st*ct*(ct2-st2)*(var1-n2)
     6       +n2*var1*(ct2*ct2-3.d0*ct2*st2+st2*st2-
     1       3.d0*st2*ct2) )/ddenom

      write(*,*)' function rrind 4'

      deriv=deriv - (dnumer/(ddenom*ddenom))*(2.d0*eps22*st*ct
     1      +n2*n2*2.d0*st*ct*(ct2-st2)+4.d0*n2*rn*dndt*ct2*st2
     2      +2.d0*ct*st*var1*var2+
     1      2.d0*rn*dndt*ct2*(var2+var1*st2)
     3      +n2*2.d0*st*ct*ct2*var1+2.d0*rn*dndt*
     1      (ct2*var2+st2*var3)
     4      +2.d0*st*ct*(var3-var2)*n2
     5      -2.d0*ct*st*(var3*var1+n2*st2*var1)
     6      +2.d0*rn*dndt*st2*(ct2*var1+var3)
     7      -4.d0*rn*dndt*st2*ct2*(var1-n2)
     8      -4.d0*st*ct*n2*var1*(ct2-st2) )

      write(*,*)' function rrind 5'
c
      ddenom=-st/droot+ct*(dndt/rn)/droot+st*deriv/droot
     1 - (ct+(dndt/rn)*st)*(dndt/rn)*deriv/(droot*droot*droot)
      rri2=n2*st*droot/ddenom
      rri2=dabs(rri2)
c
      rrind=dsqrt(rri2)

      write(*,*)' rrind', rrind
      write(*,*)'dsqrt(n2),st,droot,ddenom',dsqrt(n2),st,droot,ddenom
      write(*,*)'st,ct,dndt/rn,deriv',st,ct,dndt/rn,deriv
      write(*,*)'st/droot+ct*(dndt/rn)/droot+st*deriv/droot',
     &           st/droot+ct*(dndt/rn)/droot+st*deriv/droot
      write(*,*)'(ct+(dndt/rn)*st)*(dndt/rn)*deriv/(droot**3)',
     &           (ct+(dndt/rn)*st)*(dndt/rn)*deriv/(droot**3)
      return
      end

