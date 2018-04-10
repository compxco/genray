
c
c  This module contains subroutines from adj
c
c----------------------------------------------------------

      function current (du, ua, c, nmax, dth0, th0a, imax, imax1,
     1   fm, chi, cur)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(kind=dp) current
      integer nmax, imax, imax1
      real(kind=dp) du, c, dth0
      real(kind=dp), dimension(0:nmax - 1) :: ua
      real(kind=dp), dimension(0:imax - 1) :: th0a
      real(kind=dp), dimension(0:nmax - 1) :: fm
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: chi
      real(kind=dp), dimension(0:imax - 1) :: cur
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, i
      real(kind=dp) :: pi
ctest
      integer :: j
      real(kind=dp) :: chi_norm
C-----------------------------------------------

      pi = 4*datan(1.0d0)
      cur(:imax-1) = 0

      do n = 0, nmax - 1
         cur(:imax-1) = cur(:imax-1) + chi(:imax-1,n)*fm(n)*ua(n)**3/
     1      dsqrt(1 + (ua(n)/c)**2)           
      end do

      current=dot_product(dsin(th0a(:imax-1))*dcos(th0a(:imax-1)),cur(:
     1   imax-1))

      current = 4*pi*du*dth0*current

      return
      end function current


      subroutine decomp(chi, nmax, nmax1, imax, imax1, lmaxa, lmaxa1,
     1   th0a, dth0, legp, chil)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, nmax1, imax, imax1, lmaxa, lmaxa1
      real(kind=dp) dth0
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: chi
      real(kind=dp), dimension(0:imax - 1) :: th0a
      real(kind=dp), dimension(0:imax1 - 1,lmaxa1) :: legp
      real(kind=dp), dimension(0:nmax1 - 1,lmaxa1) :: chil
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: la, l, n, i
      real(kind=dp) :: sn
C-----------------------------------------------

      do la = 1, lmaxa
         l = max(2*la - 1,0)

         chil(:nmax-1,la) = 0

         do i = 0, imax - 1
            sn = dsin(th0a(i))*legp(i,la)
            chil(:nmax-1,la) = chil(:nmax-1,la) + chi(i,:nmax-1)*sn
         end do

         chil(:nmax-1,la) = (2*l + 1)*dth0*chil(:nmax-1,la)

      end do

      return
      end subroutine decomp


      function deltab (phi, eps)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(kind=dp) phi, eps
      real(kind=dp) deltab
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------

      deltab = 2*eps*dsin(phi/2)**2/(1 + eps*dcos(phi))

      return
      end function deltab

      subroutine error(du, ua, umax, nmax, dth0, th0a, th0max, imax,
     1   imax1, resid, chi, abserr, relerr, aerr, rerr)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, imax, imax1
      real(kind=dp) du, umax, dth0, th0max, abserr, relerr
      real(kind=dp), dimension(0:nmax - 1) :: ua
      real(kind=dp), dimension(0:imax - 1) :: th0a
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: resid, chi
      real(kind=dp), dimension(0:imax - 1) :: aerr, rerr
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, n
C-----------------------------------------------

      aerr(:imax-1) = 0
      rerr(:imax-1) = 0

      do n = 0, nmax - 1
         aerr(:imax-1) = aerr(:imax-1) + resid(:imax-1,n)**2*(ua(n)**2)
         rerr(:imax-1) =rerr(:imax-1)+(resid(:imax-1,n)/max(1.0E-16_dp,
     1      abs(chi(:imax-1,n))))**2*(ua(n)**2)
      end do

      abserr = dot_product(dsin(th0a(:imax-1)),aerr(:imax-1))
      relerr = dot_product(dsin(th0a(:imax-1)),rerr(:imax-1))

      abserr = dsqrt(abserr*du*dth0/(umax**3/3*(1 -dcos(th0max))))
      relerr = dsqrt(relerr*du*dth0/(umax**3/3*(1 -dcos(th0max))))

      return
      end subroutine error


      subroutine hinit(eps, kmax, lmaxa1, harr, garr, dthp, dl, bb)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kmax, lmaxa1
      real(kind=dp) eps, dthp
      real(kind=dp), dimension(((lmaxa1 - 1)
     1 *lmaxa1*(2*lmaxa1-1))/6+(lmaxa1- 1)*lmaxa1+lmaxa1):: harr
      real(kind=dp),dimension((((lmaxa1+3)-1)*
     1   (lmaxa1+3))/2+(lmaxa1+1))::garr
      real(kind=dp), dimension(0:kmax - 1) :: dl, bb
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k, la, la1, la2, l, l1
      real(kind=dp) :: b, pi
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(kind=dp) , EXTERNAL :: deltab
C-----------------------------------------------

      pi = 4*datan(1.0d0)

C
C  Initialize |harr| and set the edge values of |garr| to zero.
C
  
      do la = 1, lmaxa1
         do la1 = 1, la
            harr((la-1)*la*(2*la-1)/6+(la1-1)*la+1:la*(la*(2*(la-1)-1)+1
     1         )/6+la1*la) = 0
         end do
      end do

      do la = -2, lmaxa1 - 1
         garr((((la+3)-1)*(la+3))/2+(0+1)) = 0
      end do
      do la = -1, lmaxa1 - 1
         garr((((la+3)-1)*(la+3))/2+(la+1+1)) = 0
         garr((((la+3)-1)*(la+3))/2+(la+2+1)) = 0
      end do
      garr((((1+3)-1)*(1+3))/2+(1+1)) = 1

c      write(*,*)'hinit 0 garr',garr

      do k = 0, kmax - 1
c     phi=2*pi*(k+0.25e0)/kmax
         b = bb(k)
 
C  Use the recursion relation for 'g' to fill garr
C
         do la = 2, lmaxa1
            l = max(2*la - 1,0)
            do la1 = 1, la
               l1 = max(2*la1 - 1,0)
               garr((((la+3)-1)*(la+3))/2+(la1+1)) = (2*(2*l - 3)*(l**2
     1             - 3*l + 1)*garr((((la-1+3)-1)*(la-1+3))/2+(la1+1))-(l
     2            -3)*(l-2)*(2*l-1)*garr((((la-2+3)-1)*(la-2+3))/2+(la1+
     3            1)))/(l*(l - 1)*(2*l - 5)) + (2*l - 3)*(2*l - 1)*b*((
     4            l1 + 1)*(l1 + 2)*garr((((la-1+3)-1)*(la-1+3))/2+(la1+1
     5            +1))/((2*l1+3)*(2*l1+5))-2*(l1**2+l1-1)*garr((((la-1+3
     6            )-1)*(la-1+3))/2+(la1+1))/((2*l1-1)*(2*l1+3))+l1*(l1-1
     7            )*garr((((la-1+3)-1)*(la-1+3))/2+(la1-1+1))/((2*l1-3)*
     8            (2*l1-1)))/(l*(l - 1))          
            end do
         end do

c         write(*,*)'hinit 1 garr',garr

         do la = 1, lmaxa1
            do la1 = 1, la
               harr((la-1)*la*(2*la-1)/6+(la1-1)*la+1:la*(la*(2*(la-1)-1
     1            )+1)/6+la1*la) = harr((la-1)*la*(2*la-1)/6+(la1-1)*la+
     2            1:la*(la*(2*(la-1)-1)+1)/6+la1*la) + b*garr(((la+3)*(
     3            la+3)-3-la)/2+la1+1)*garr((la*(la+5)+6)/2+2:la+1+(la*(
     4            la+5)+6)/2)*dl(k)                
            end do
         end do
      end do

      do la = 1, lmaxa1
         do la1 = 1, la
            do la2 = 1, la
               harr(((la-1)*la*(2*la-1))/6+(la1-1)*la+la2) = (2*max(2*la
     1             - 1,0) + 1)/real(2*max(2*la2 - 1,0) + 1)*harr(((la-1)
     2            *la*(2*la-1))/6+(la1-1)*la+la2)
            end do
         end do
      end do

c      write(*,*)'hinit 3 harr',harr

      return
      end subroutine hinit


      subroutine initjy(ua, c, nmax, nmax1, lmaxa, lmaxa1, legjy, jarray
     1   , djarray, lmax1)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, nmax1, lmaxa, lmaxa1, lmax1
      real(kind=dp) c
      real(kind=dp), dimension(0:nmax - 1) :: ua
      real(kind=dp), dimension(0:nmax1 - 1,0:6,0:3,0:lmaxa1) :: legjy
      real(kind=dp), dimension((-lmax1) - 1:lmax1,
     1   0:6) :: jarray, djarray
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, mult, la, l, a
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL legjn
C-----------------------------------------------

      do n = 0, nmax - 1
         call legjn(ua(n),c,max(2*lmaxa-1,0),jarray,djarray,lmax1)
         do la = 0, lmaxa
            l = max(2*la - 1,0,0)
            mult = (-1)**(l + 1)
            legjy(n,:6,0,la) = jarray(l,:6)
            legjy(n,:6,1,la) = mult*jarray((-l)-1,:6)
            legjy(n,:6,2,la) = djarray(l,:6)
            legjy(n,:6,3,la) = mult*djarray((-l)-1,:6)
         end do
      end do

      return
      end subroutine initjy


      subroutine initp(x, imax, imax1, lmaxa, lmaxa1, legp, leg0)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer imax, imax1, lmaxa, lmaxa1
      real(kind=dp), dimension(0:imax - 1) :: x
      real(kind=dp), dimension(0:imax1 - 1,lmaxa1) :: legp
      real(kind=dp), dimension(0:imax - 1) :: leg0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, l, la
C-----------------------------------------------

      if (lmaxa < 1) return

      la = 1
      leg0(:imax-1) = 1
      legp(:imax-1,la) = x(:imax-1)

      do la = 2, lmaxa
         l = max(2*la - 1,0)
         leg0(:imax-1) = (((2*l) - 3)*x(:imax-1)*legp(:imax-1,la-1)-(l-2
     1      )*leg0(:imax-1))/(l - 1)
         legp(:imax-1,la) = (((2*l) - 1)*x(:imax-1)*leg0(:imax-1)-(l-1)*
     1      legp(:imax-1,la-1))/l
      end do

      return
      end subroutine initp


      subroutine isotrop(ua, du, c, fm, chil, nmax, nmax1, legjy, psid,
     1   ze, zi, duug, fu, dtt, work)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, nmax1
      real(kind=dp) du, c, ze, zi
      real(kind=dp), dimension(0:nmax - 1) :: ua, fm, chil
      real(kind=dp), dimension(0:nmax1 - 1,0:6,0:3) :: legjy
      real(kind=dp), dimension(0:nmax1 - 1,0:4,0:1) :: psid
      real(kind=dp), dimension(0:nmax) :: duug
      real(kind=dp), dimension(0:nmax - 1) :: fu, dtt
      real(kind=dp), dimension(0:nmax + (nmax1 + 1)*7*2 - 1) :: work
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: p0 = 0
      integer, parameter :: p02 = 1
      integer, parameter :: p022 = 2
      integer, parameter :: p1 = 3
      integer, parameter :: p11 = 4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n
      real(kind=dp) :: pi, g, u
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL pot
C-----------------------------------------------
      call pot (ua, du, c, fm, chil, nmax, nmax1, legjy, psid, work(0),
     1   work(nmax))

      pi = 4*datan(1.0d0)
      do n = 0, nmax - 1
         u = ua(n)
         g = dsqrt(1 + (u/c)**2)
         duug(n) = ze*4*pi*g/u*(2*g**2*psid(n,p02,1)-u*psid(n,p0,0)-8*g
     1      **2/c**2*psid(n,p022,1)+8*u/c**4*psid(n,p022,0))
         dtt(n) = ze*4*pi/(g*u)*((-g**2*psid(n,p02,1))-u/c**2*psid(n,p02
     1      ,0)+4*g**2/c**2*psid(n,p022,1)-4*u/c**4*psid(n,p022,0)) + zi
     2      *g/(2*u)
         fu(n) = ze*4*pi*g*((-psid(n,p1,1))+2/c**2*psid(n,p11,1))

      end do



      duug(nmax) = (3*duug(nmax-1)-duug(nmax-2))/2
      duug(nmax-1:1:(-1)) = (duug(nmax-1:1:(-1))+duug(nmax-2:0:(-1)))/2
      duug(0) = duug(0)

      return
      end subroutine isotrop


      function lam1 (eps, th0, kmax, dthp, dl, b)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(kind=dp) lam1
      integer kmax
      real(kind=dp) eps, th0, dthp
      real(kind=dp), dimension(0:kmax - 1) :: dl, b
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k
      real(kind=dp) :: deltabv, sum, costh
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(kind=dp) , EXTERNAL :: deltab
C-----------------------------------------------

c     pi=4*datan(1.0d0)
      sum = 0
      do k = 0, kmax - 1
c     phi=2*pi*(k+0.25e0)/kmax
         deltabv = b(k) - 1.d0

c----------------------------------------------------------------
c       costh=cos(theta)/cos(theta_0)
c---------------------------------------------------------------
        if (dabs(deltabv).lt.1.d-12) then
            costh=1.d0
         else  
           costh = dsqrt(1 - min(1.0d0,deltabv*dtan(th0)**2))
         endif 

c         costh = dsqrt(1 - min(1.0_dp,deltabv*dtan(th0)**2))
         sum = sum + 1/costh*dl(k)
      end do
      lam1 = sum

      return
      end function lam1


      function lam2 (eps, th0, kmax, dthp, dl, b)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(kind=dp) lam2
      integer kmax
      real(kind=dp) eps, th0, dthp
      real(kind=dp), dimension(0:kmax - 1) :: dl, b
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k
      real(kind=dp) :: deltabv, sum, costh
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(kind=dp) , EXTERNAL :: deltab
C-----------------------------------------------
c      write(*,*)'in lam2 eps,th0,kmax,dthp',eps,th0,kmax,dthp
c     pi=4*datan(1.0d0)
      sum = 0
      do k = 0, kmax - 1
c     phi=2*pi*(k+0.25e0)/kmax
         deltabv = b(k) - 1. 
c----------------------------------------------------------------
c       costh=cos(theta)/cos(theta_0)
c---------------------------------------------------------------
        if (dabs(deltabv).lt.1.d-12) then
            costh=1.d0
         else  
           costh = dsqrt(1 - min(1.0d0,deltabv*dtan(th0)**2))
         endif 

c         costh = dsqrt(1 - min(1.0E0_dp,deltabv*tan(th0)**2)) 
c         write(*,*)'b(k),deltabv,th0,dtan(th0),tan(th0)',
c     &              b(k),deltabv,th0,dtan(th0),tan(th0)
c         write(*,*)'deltabv*dtan(th0)**2',deltabv*tan(th0)**2
c         write(*,*)'min(1.0E0_dp,deltabv*dtan(th0)**2)',
c     &              min(1.0E0_dp,deltabv*dtan(th0)**2) 

cSAP070812       
c         sum = sum + costh/(1 + deltabv)*dl(k)/b(k)
         sum = sum + costh/(1 + deltabv)*dl(k)

c         write(*,*)'lam2 k,costh,deltabv,dl(k),b(k),dl(k)/b(k)',
c     &   k,costh,deltabv,dl(k),b(k),dl(k)/b(k)
c         write(*,*)'lam2 sum',sum
      end do
      lam2 = sum

      return
      end function lam2


      subroutine legjn(u, c, lmax, jarray, djarray, lmax1)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer lmax, lmax1
      real(kind=dp) u, c
      real(kind=dp), dimension((-lmax1) - 1:lmax1,
     1   0:6) :: jarray, djarray
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: j0 = 0
      integer, parameter :: j1 = 1
      integer, parameter :: j2 = 2
      integer, parameter :: j02 = 3
      integer, parameter :: j11 = 4
      integer, parameter :: j22 = 5
      integer, parameter :: j022 = 6
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, a, amax
      real(kind=dp) :: scale, g
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL legjs
C-----------------------------------------------

      amax = 2
      call legjs (u/c, lmax + 2, amax, jarray, lmax1)

      scale = 1
      do l = 1, lmax + 2
         scale = scale*u/(2*l + 1)
         jarray(l,:amax) = scale*jarray(l,:amax)
      end do

      scale = 1
      do l = -1, (-lmax) - 2, -1
         scale = scale/u*(2*l + 3)
         jarray(l,:amax) = scale*jarray(l,:amax)
      end do

      do l = (-lmax) - 2, lmax
         jarray(l,j02) = u/2*jarray(l+1,1)
         jarray(l,j11) = u/2*jarray(l+1,0)
         jarray(l,j22) = u/2*jarray(l+1,1) + (u/c)**2/2*jarray(l+2,0)
         jarray(l,j022) = u**2/8*jarray(l+2,0)
      end do

      g = dsqrt(1 + (u/c)**2)
      do a = 0, 6
         do l = (-lmax) - 1, lmax
            djarray(l,a) = jarray(l-1,a)/g - (l + 1)/u*jarray(l,a)
         end do
      end do

      return
      end subroutine legjn


      subroutine legjs(z, lmax, amax, jarray, lmax1)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer lmax, amax, lmax1
      real(kind=dp) z
      real(kind=dp), dimension((-lmax1) - 1:lmax1,0:amax) :: jarray
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, a, lstart
      real(kind=dp) :: z2, g, ja0, ja1, ja2, ratio, sigma
C-----------------------------------------------

      if (amax<0 .or. lmax<0) return
      if (lmax > lmax1) return

      z2 = z**2
      g = dsqrt(1 + z2)

      do a = 1, amax
C**************a=0 is excluded
         ja2 = 1
C********ja(a-1,a)
         ja1 = g
C**********ja(a-2,a)
         do l = a - 3, lmax - 1, -1
C
C   Get things started by iterating down to $|lmax|-1$ if need be.
C   This is from equation (A29)
C
            ja0=g*ja1-z2*(l-a+2)*(l+a+2)*ja2/((2*l+3)*(2*l+5))
            ja2 = ja1
            ja1 = ja0
         end do
         l = min(lmax,a - 1)
         jarray(l,a) = ja2
         l = l - 1
         jarray(l,a) = ja1
         do l = min(lmax,a - 1) - 2, max((-a),(-lmax) - 1), -1
            jarray(l,a) = g*jarray(l+1,a) - z2*(l - a + 2)*(l + a + 2)*
     1         jarray(l+2,a)/((2*l + 3)*(2*l + 5))
         end do
      end do

C
C   Calculate the scaled Legendre functions for the other easy case, namely
C   $l<-a$.
C
      do a = 0, min(amax,lmax)
         l = (-a) - 1
         jarray(l,a) = 1
         l = l - 1
         if (l >= (-lmax) - 1) jarray(l,a) = g
         do l = (-a) - 3, (-lmax) - 1, -1
            jarray(l,a) = g*jarray(l+1,a) - z2*(l - a + 2)*(l + a + 2)*
     1         jarray(l+2,a)/((2*l + 3)*(2*l + 5))
         end do
      end do

C
C  Calculate scaled Legendre functions for -a.<= l <=a
C
C
C   Calculate $\ja{l\ge0,0}$.  Depending on the magnitude of the argument we
C  either use backwards recursion (for $z \le |lmax|+1$) or forwards
C  recursion using $\ja{-1,0}=1$ and $\ja{0,0} = (\sinh^{-1} z)/z$.
C
      a = 0
      if (abs(z) <= lmax + 1) then
C
C  Estimate where to start the recursion.  The 10 is a extra guard.
C
         lstart=lmax+log(1.0E-16)/(2*log((1.0E-16+abs(z))/(g+1)))+10
         ratio = 2/(g + 1)
         do l = lstart, lmax - 1, -1
            ratio=1/(g-z2*(l-a+2)*(l+a+2)*ratio/((2*l+3)*(2*l+5)))
C***At this point, |ratio| is approximation for $\ja{l+1,a}/\ja{l,a}$. **
         end do
         jarray(lmax,a) = ratio
         ja0 = 1
         do l = lmax - 2, -1, -1
            jarray(l+1,a) = ja0
            ja0 = g*jarray(l+1,a) - z2*(l - a + 2)*(l + a + 2)*jarray(l+
     1         2,a)/((2*l + 3)*(2*l + 5))
         end do
         ja0 = 1/ja0
         jarray(0:lmax,a) = ja0*jarray(0:lmax,a)
      else
C********** Use forward recurrence for $z>|lmax|+1$. **
         sigma = log(g + abs(z))
C******$\sinh^{-1}(\abs z)$ **
         l = 0
         jarray(l,a) = sigma/abs(z)
C********* Starting value.  We assume
C            that $\ja{-1,0}$ has already been set to $1$. */


         do l = 1, lmax
C***********Equation (A29)***************
            jarray(l,a) = (g*jarray(l-1,a)-jarray(l-2,a))*(2*l - 1)*(2*l
     1          + 1)/(z2*(l - a)*(l + a))
         end do
      endif
C
C   Do the other cases, i.e., $a>0$, $l\ge a$.
C
C  < Do the other cases, i.e., $a>0$, $l\ge a$ @>=
C
      do a = 1, min(lmax,amax)
         l = lmax
         jarray(l,a) = ((2*l + 1)*jarray(l-1,a-1)-(l+1-a)*g*jarray(l,a-1
     1      ))/(l + a)
         l = lmax - 1
         if (l == 0) then
C
C    implies that a=1. Use ja(-1,0)=1
C
            ja0 = ((2*l + 1) - (l + 1 - a)*g*jarray(l,a-1))/(l + a)
         else
            ja0=((2*l+1)*jarray(l-1,a-1)-(l+1-a)*g*jarray(l,a-1))/(l+a)
         endif
         do l = lmax - 2, a - 1, -1
            jarray(l+1,a) = ja0
            ja0 = g*jarray(l+1,a) - z2*(l - a + 2)*(l + a + 2)*jarray(l+
     1         2,a)/((2*l + 3)*(2*l + 5))
         end do
         ja0 = 1/ja0
         jarray(a:lmax,a) = ja0*jarray(a:lmax,a)
      end do

      return
      end subroutine legjs


      subroutine matinit(du, ua, ug, nmax, nmax1, dth0, th0a, th0g, imax
     1   , imax1, duug, fu, dtt, lam1a, lam2g, uinv, th0inv, rho)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, nmax1, imax, imax1
      real(kind=dp) du, dth0, rho
      real(kind=dp), dimension(0:nmax - 1) :: ua
      real(kind=dp), dimension(0:nmax) :: ug
      real(kind=dp), dimension(0:imax - 1) :: th0a
      real(kind=dp), dimension(0:imax) :: th0g
      real(kind=dp), dimension(0:nmax) :: duug
      real(kind=dp), dimension(0:nmax - 1) :: fu, dtt
      real(kind=dp), dimension(0:imax - 1) :: lam1a
      real(kind=dp), dimension(0:imax) :: lam2g
      real(kind=dp), dimension(0:nmax1 - 1,-1:1) :: uinv
      real(kind=dp), dimension(0:imax1 - 1,-1:1) :: th0inv
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, i
C-----------------------------------------------
C
C  Define matrices for inversion of diffusion and friction operator.  These
C  are one-dimensional because the operator is independent of $\theta_0$.
C
      uinv(:nmax-1,-1) = ug(:nmax-1)**2*duug(:nmax-1)/((ua(:nmax-1)**2)*
     1   (du**2))
      uinv(:nmax-1,1) = ug(1:nmax)**2*duug(1:nmax)/((ua(:nmax-1)**2)*(du
     1   **2))
      uinv(:nmax-1,0) = (-uinv(:nmax-1,-1)) - uinv(:nmax-1,1)
      uinv(:nmax-1,-1) = uinv(:nmax-1,-1) - fu(:nmax-1)/(2*du)
      uinv(:nmax-1,1) = uinv(:nmax-1,1) + fu(:nmax-1)/(2*du)
      n = 0
      uinv(n,0) = uinv(n,0) - uinv(n,-1)
      uinv(n,-1) = 0
      n = nmax - 1
      uinv(n,-1) = uinv(n,-1) - rho*uinv(n,1)
      uinv(n,0) = uinv(n,0) + (1 + rho)*uinv(n,1)
      uinv(n,1) = 0

      uinv(:nmax-1,-1) = uinv(:nmax-1,-1)*(ua(:nmax-1)**2)/dtt(:nmax-1)
      uinv(:nmax-1,0) = uinv(:nmax-1,0)*(ua(:nmax-1)**2)/dtt(:nmax-1)
      uinv(:nmax-1,1) = uinv(:nmax-1,1)*(ua(:nmax-1)**2)/dtt(:nmax-1)
C
C   Define matrices for inversion of pitch-angle scattering operator.  These
C  are one-dimensional because, after divided out $D_{\theta\theta,0}/u^2$,
C  the operator is independent of $u$.
C
      th0inv(:imax-1,-1) =dsin(th0g(:imax-1))*lam2g(:imax-1)/((lam1a(:
     1   imax-1)*(dsin(th0a(:imax-1))))*(dth0**2))

cSAP080320 in lam1a it was bellow:  lam1a(:imax-
c     - was at 73 position
      th0inv(:imax-1,1) = dsin(th0g(1:imax))*lam2g(1:imax)/((lam1a(:imax
     1   -1)*(dsin(th0a(:imax-1))))*(dth0**2))
      th0inv(:imax-1,0) = (-th0inv(:imax-1,-1)) - th0inv(:imax-1,1)
      i = 0
      i = imax - 1  

c      write(*,*)'th0inv(imax-1,-1),th0inv(imax-1,1),th0inv(imax-1,0)',
c     &th0inv(imax-1,-1),th0inv(imax-1,1),th0inv(imax-1,0)
c      write(*,*)'th0inv(imax-1,-1),th0inv(imax-1,1),th0inv(imax-1,0)',
c     &th0inv(imax-1,-1),th0inv(imax-1,1),th0inv(imax-1,0)
c      write(*,*)'dsin(th0g(imax)),lam2g(imax)',
c     &dsin(th0g(imax)),lam2g(imax)
c      write(*,*)'dsin(th0a(imax-1)),lam1a(imax-1)',
c     &dsin(th0a(imax-1)),lam1a(imax-1)
     
      th0inv(i,0) = th0inv(i,0) - th0inv(i,1)

c      write(*,*)'th0inv(imax-1,0)',th0inv(imax-1,0)

      th0inv(i,1) = 0

c      write(*,*)'in matinit  th0inv'

c      do i=1,imax-1
c         write(*,*)'i,th0g(i),th0a(i),lam2g(i),lam1a(i)',
c     &              i,th0g(i),th0a(i),lam2g(i),lam1a(i)
c         write(*,*)'th0inv(i,-1),th0inv(i,0),th0inv(i,1)',
c     &              th0inv(i,-1),th0inv(i,0),th0inv(i,1)
c      enddo

C* :86 *
*line 1805 "adj.web"


      return
      end subroutine matinit



      subroutine maxwel(ua, du, c, nmax, fm, t)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax
      real(kind=dp) du, c, t
      real(kind=dp), dimension(0:nmax - 1) :: ua, fm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n
      real(kind=dp) :: sum, pi
C-----------------------------------------------

      fm(:nmax-1)=dexp((-ua(:nmax-1)**2/(dsqrt(1+(ua(:nmax-1)/c)**2)+1)/
     1   t))
      sum = dot_product(ua(:nmax-1)**2,fm(:nmax-1))
      pi = 4*datan(1.0d0)
      sum = 4*pi*du*sum
      fm(:nmax-1) = fm(:nmax-1)/sum

      return
      end subroutine maxwel

      subroutine operinit
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------

      open(unit=20, file='adjout', form='unformatted', status='unknown')
      open(unit=29, file='adjinp', form='unformatted', status='old')

      return
      end subroutine operinit


      subroutine pot(ua,du,c,fm,chil,nmax,nmax1,legjy,psid,fla,ijy)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, nmax1
      real(kind=dp) du, c
      real(kind=dp), dimension(0:nmax - 1) :: ua, fm, chil
      real(kind=dp), dimension(0:nmax1 - 1,0:6,0:3) :: legjy
      real(kind=dp), dimension(0:nmax1 - 1,0:4,0:1) :: psid
      real(kind=dp), dimension(0:nmax - 1) :: fla
      real(kind=dp), dimension(0:nmax1,0:6,0:1) :: ijy
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: p0 = 0
      integer, parameter :: p02 = 1
      integer, parameter :: p022 = 2
      integer, parameter :: p1 = 3
      integer, parameter :: p11 = 4
      integer, parameter :: j0 = 0
      integer, parameter :: j1 = 1
      integer, parameter :: j2 = 2
      integer, parameter :: j02 = 3
      integer, parameter :: j11 = 4
      integer, parameter :: j22 = 5
      integer, parameter :: j022 = 6
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, a
C-----------------------------------------------


C   Define various indefinite integrals needed in the definition of the
C  potentials.  First do the integration to cell edges and then interpolate
C  to the cell centers.

      fla(:nmax-1) = fm(:nmax-1)*chil(:nmax-1)*ua(:nmax-1)**2/dsqrt(1 +(
     1   ua(:nmax-1)/c)**2)*du/2

      do a = 0, 6
         ijy(1:nmax,a,0) = fla(:nmax-1)*legjy(:nmax-1,a,0)
         ijy(:nmax-1,a,1) = fla(:nmax-1)*legjy(:nmax-1,a,1)
         ijy(0,a,0) = 0
         ijy(nmax,a,1) = 0
      end do

      do n = 1, nmax
         ijy(n,:6,0) = ijy(n-1,:6,0) + ijy(n,:6,0)
      end do
      do n = nmax - 1, 0, -1
         ijy(n,:6,1) = ijy(n+1,:6,1) + ijy(n,:6,1)
      end do

      ijy(:nmax-1,:6,0) = ijy(:nmax-1,:6,0) + ijy(1:nmax,:6,0)
      ijy(:nmax-1,:6,1) = ijy(:nmax-1,:6,1) + ijy(1:nmax,:6,1)
C
C  Define various indefinite integrals needed in the definition of
C   the potentials. First do the integration to cell edges and then
C   interpolate to the cell centers.
C
      psid(:nmax-1,p0,0) = legjy(:nmax-1,j0,1)*ijy(:nmax-1,j0,0) + ijy(:
     1   nmax-1,j0,1)*legjy(:nmax-1,j0,0)
      psid(:nmax-1,p02,0) = legjy(:nmax-1,j0,1)*ijy(:nmax-1,j02,0) +
     1   legjy(:nmax-1,j02,1)*ijy(:nmax-1,j2,0) + ijy(:nmax-1,j0,1)*
     2   legjy(:nmax-1,j02,0) + ijy(:nmax-1,j02,1)*legjy(:nmax-1,j2,0)
      psid(:nmax-1,p022,0) = legjy(:nmax-1,j0,1)*ijy(:nmax-1,j022,0) +
     1   legjy(:nmax-1,j02,1)*ijy(:nmax-1,j22,0) + legjy(:nmax-1,j022,1)
     2   *ijy(:nmax-1,j2,0) + ijy(:nmax-1,j0,1)*legjy(:nmax-1,j022,0) +
     3   ijy(:nmax-1,j02,1)*legjy(:nmax-1,j22,0) + ijy(:nmax-1,j022,1)*
     4   legjy(:nmax-1,j2,0)
      psid(:nmax-1,p1,0) = legjy(:nmax-1,j1,1)*ijy(:nmax-1,j1,0) + ijy(:
     1   nmax-1,j1,1)*legjy(:nmax-1,j1,0)
      psid(:nmax-1,p11,0) = legjy(:nmax-1,j1,1)*ijy(:nmax-1,j11,0) +
     1   legjy(:nmax-1,j11,1)*ijy(:nmax-1,j1,0) + ijy(:nmax-1,j1,1)*
     2   legjy(:nmax-1,j11,0) + ijy(:nmax-1,j11,1)*legjy(:nmax-1,j1,0)

      psid(:nmax-1,p0,1) = legjy(:nmax-1,j0,3)*ijy(:nmax-1,j0,0) + ijy(:
     1   nmax-1,j0,1)*legjy(:nmax-1,j0,2)
      psid(:nmax-1,p02,1) = legjy(:nmax-1,j0,3)*ijy(:nmax-1,j02,0) +
     1   legjy(:nmax-1,j02,3)*ijy(:nmax-1,j2,0) + ijy(:nmax-1,j0,1)*
     2   legjy(:nmax-1,j02,2) + ijy(:nmax-1,j02,1)*legjy(:nmax-1,j2,2)
      psid(:nmax-1,p022,1) = legjy(:nmax-1,j0,3)*ijy(:nmax-1,j022,0) +
     1   legjy(:nmax-1,j02,3)*ijy(:nmax-1,j22,0) + legjy(:nmax-1,j022,3)
     2   *ijy(:nmax-1,j2,0) + ijy(:nmax-1,j0,1)*legjy(:nmax-1,j022,2) +
     3   ijy(:nmax-1,j02,1)*legjy(:nmax-1,j22,2) + ijy(:nmax-1,j022,1)*
     4   legjy(:nmax-1,j2,2)
      psid(:nmax-1,p1,1) = legjy(:nmax-1,j1,3)*ijy(:nmax-1,j1,0) + ijy(:
     1   nmax-1,j1,1)*legjy(:nmax-1,j1,2)
      psid(:nmax-1,p11,1) = legjy(:nmax-1,j1,3)*ijy(:nmax-1,j11,0) +
     1   legjy(:nmax-1,j11,3)*ijy(:nmax-1,j1,0) + ijy(:nmax-1,j1,1)*
     2   legjy(:nmax-1,j11,2) + ijy(:nmax-1,j11,1)*legjy(:nmax-1,j1,2)

      return
      end subroutine pot


      subroutine reactall(ua, du, c, fm, t, chil, nmax, nmax1, imax,
     1   imax1, lam1a, legjy, psid, legp, harr, ze, lmaxa, lmaxa1, react
     2   , reactl, chila, work)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, nmax1, imax, imax1, lmaxa, lmaxa1
      real(kind=dp) du, c, t, ze
      real(kind=dp), dimension(0:nmax - 1) :: ua, fm
      real(kind=dp), dimension(0:nmax1 - 1,lmaxa1) :: chil
      real(kind=dp), dimension(0:imax - 1) :: lam1a
      real(kind=dp), dimension(0:nmax1 - 1,0:6,0:3,lmaxa1) :: legjy
      real(kind=dp), dimension(0:nmax1 - 1,0:4,0:1) :: psid
      real(kind=dp), dimension(0:imax1 - 1,lmaxa1) :: legp
      real(kind=dp), dimension(((lmaxa1 - 1)*
     1  lmaxa1*(2*lmaxa1-1))/6+(lmaxa1-1)*lmaxa1+lmaxa1):: harr
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: react
      real(kind=dp), dimension(0:nmax - 1) :: reactl, chila
      real(kind=dp), dimension(0:nmax + (nmax1 + 1)*7*2 - 1) :: work
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: la, la1, la2, n, i
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL reaction
C-----------------------------------------------

      react(:imax-1,:nmax-1) = 0

      do la1 = 1, lmaxa
         reactl(:nmax-1) = 0
         do la = la1, lmaxa
            chila(:nmax-1) = 0
            do la2 = 1, la
               chila(:nmax-1) = chila(:nmax-1) + harr(((la-1)*la*(2*la-1
     1            ))/6+(la1-1)*la+la2)*chil(:nmax-1,la2)
            end do
            call reaction (ua, du, c, fm, t, chila, nmax, nmax1, legjy(0
     1         ,0,0,la), psid, ze, reactl, max(2*la - 1,0), work)
         end do

         do n = 0, nmax - 1
            react(:imax-1,n) = react(:imax-1,n) + reactl(n)*legp(:imax-1
     1         ,la1)
         end do
      end do

      do n = 0, nmax - 1
         react(:imax-1,n) = react(:imax-1,n)/lam1a(:imax-1)
      end do

      return
      end subroutine reactall


      subroutine reaction(ua, du, c, fm, t, chil, nmax, nmax1, legjy,
     1   psid, ze, reactl, l, work)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, nmax1, l
      real(kind=dp) du, c, t, ze
      real(kind=dp), dimension(0:nmax - 1) :: ua, fm, chil
      real(kind=dp), dimension(0:nmax1 - 1,0:6,0:3) :: legjy
      real(kind=dp), dimension(0:nmax1 - 1,0:4,0:1) :: psid
      real(kind=dp), dimension(0:nmax - 1) :: reactl
      real(kind=dp), dimension(0:nmax + (nmax1 + 1)*7*2 - 1) :: work
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: p0 = 0
      integer, parameter :: p02 = 1
      integer, parameter :: p022 = 2
      integer, parameter :: p1 = 3
      integer, parameter :: p11 = 4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n
      real(kind=dp) :: pi, u, g
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL pot
C-----------------------------------------------

      call pot (ua, du, c, fm, chil, nmax, nmax1, legjy, psid, work(0),
     1   work(nmax))

      pi = 4*datan(1.0d0)
      do n = 0, nmax - 1
         u = ua(n)
         g = dsqrt(1 + (u/c)**2)
         reactl(n) = reactl(n) + ze*4*pi*(fm(n)*chil(n)/g-u/t*psid(n,p1,
     1      1)-2/(c**2*g)*psid(n,p1,0)+2*u/(c**2*t)*psid(n,p11,1)+u/t*
     2      psid(n,p0,1)-(u**2/(g*t**2)-1/t)*psid(n,p0,0)+(2*g*u/t**2-2*
     3      u/(c**2*t))*psid(n,p02,1)-(l*(l+1)/(g*t**2)-2/(c**2*t))*psid
     4      (n,p02,0)-8*g*u/(c**2*t**2)*psid(n,p022,1)+4*((l+2)*(l-1)/g+
     5      2*g)/(c**2*t**2)*psid(n,p022,0))
      end do

      return
      end subroutine reaction


      subroutine residue(ua, nmax, nmax1, csth0a, imax, imax1, chi, dtt
     1   , c, lam1a, uinv, th0inv, resid)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, nmax1, imax, imax1
      real(kind=dp) c
      real(kind=dp), dimension(0:nmax - 1) :: ua
      real(kind=dp), dimension(0:imax - 1) :: csth0a
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: chi
      real(kind=dp), dimension(0:nmax - 1) :: dtt
      real(kind=dp), dimension(0:imax - 1) :: lam1a
      real(kind=dp), dimension(0:nmax1 - 1,-1:1) :: uinv
      real(kind=dp), dimension(0:imax1 - 1,-1:1) :: th0inv
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: resid
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, i
      real(kind=dp) chi_nm_i,chi_n_i,chi_np_i,chi_n_im,chi_n_ip
C-----------------------------------------------
C   Add in electric field term
C  $$\frac{v \cos\theta_0}{\lambda_1}. $$
C
      do n = 0, nmax - 1
         resid(:imax-1,n) = resid(:imax-1,n) + csth0a(:imax-1)/lam1a(:
     1      imax-1)*ua(n)/dsqrt(1 + (ua(n)/c)**2)
      end do
C
C   Normalize to pitch-angle diffusion coefficient and add in
C   differential terms.
C
CSAP080222
c      do n = 0, nmax - 1
c         resid(:imax-1,n) = resid(:imax-1,n)*ua(n)**2/dtt(n) + chi(:imax
c    1      -1,n-1)*uinv(n,-1) + chi(:imax-1,n)*uinv(n,0) + chi(:imax-1,
c    2      n+1)*uinv(n,1) + chi(-1:imax-2,n)*th0inv(:imax-1,-1) + chi(:
c     3      imax-1,n)*th0inv(:imax-1,0) + chi(1:imax,n)*th0inv(:imax-1,1
c     4      )
c      end do

      do n = 0, nmax - 1
         do i =0, imax - 1

            if (n.eq.0) then 
               chi_nm_i=0.d0 
            else
               chi_nm_i=chi(i,n-1)
            endif

            chi_n_i=chi(i,n)

            if (n.eq.(nmax-1)) then 
               chi_np_i=0.d0 
            else
               chi_np_i=chi(i,n+1)
            endif

            if (i.eq.0) then 
               chi_n_im=0.d0 
            else
               chi_n_im=chi(i-1,n)
            endif

            if (i.eq.imax) then 
               chi_n_ip=0.d0 
            else
               chi_n_ip=chi(i+1,n)
            endif
            
            resid(i,n) = resid(i,n)*ua(n)**2/dtt(n) +
     1      chi_nm_i*uinv(n,-1) + chi_n_i*uinv(n,0) +
     2      chi_np_i*uinv(n,1)  + chi_n_im*th0inv(i,-1) +
     3      chi_n_i*th0inv(i,0) + chi_n_ip*th0inv(i,1)

 
         end do
      end do

      return
      end subroutine residue


      subroutine runa(t, ze, nmax, umax, imax, kmax, dt,
     1   tmax, alpha, rho, lmax, tskip, npsi0, dla, ba, alenb, thetps,
     2   thetpe, psimx, psimn, iout3, iout4, iout5, iout7, iout8, iswchi
     3   , aerrmx, rerrmx,
     &   dens_averaged_ar)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
C
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-------------------------------------------   
       integer
     & nmax,       !number of momentum mesh points : ua(0:nmax-1)  
     & imax,       !number of pitch angle mesh points: th0a(0:imax1-1)
     & kmax,       !number of integration points for field line integration.
                   !kmax ~ 1000 
     & tmax,       !is the maximal number of time steps
     & lmax,       !is the maximal number of Legandre harmonic
                   !lmax=1 corresponds to collision with a Maxwellian.
                   !          Usually lmax=21
     & tskip,
     & npsi0,      !number of radial points
     & iout3,      !number of input file='adjinp',
     & iout4,      !number of output file='adjout',
     & iout5,      !write (iout5, 3000) psis(np), dene, teme
     & iout7,
     & iout8,      !if(np=1)write(iout8,3030)((chi(i,n),i=0,imax-1),n=0,nmax-1)
     & iswchi      !=0  chi(:imax-1,:nmax-1) = 0
                   !else
                   !read(iout7,3030) ((chi(i,n),i=0,imax-1),n=0,nmax-1)

      real(kind=dp)
     & t,       ! is the temperature in units of m*u_n**2
                ! If you want velocities normalized to sqrt(T/m), set t=1.   
     & ze,      ! is the electron charge state. It is usually =1, 
                ! but for Lorentz limit, ion charge Zi ==> infinity,
                ! set  ze=0 and zi=1
     & umax,    ! is the maximal value of u_0. 
                ! The momentum grid is spaced between 0 and umax.
     & dt,      ! is the time step
     & alpha,   ! is a weighting factor in the time stepping algorithm.
                !       =0 for explicit algorithm (usually =0.55)
     & rho,     ! governs the treatment  of the boundary at u=umax
                ! for rho=1 the second derivative from the adjoint function chi
                ! d^2(chi)/dx^**2(at u=umax)=0 
                ! for rho=0 the first derivative =0 
     & psimx,
     & psimn,
     & aerrmx,
     & rerrmx,
     & dens_averaged_ar(1:npsi0) ! bounce averaged density from
                                 ! the relativistic maxwelian distribution
 
c     integer nmax, imax, kmax, tmax, lmax, tskip, npsi0, iout3, iout4,
c    1   iout5, iout7, iout8, iswchi
c     real(kind=dp) t, c2, ze, zi, eps, umax, dt, alpha, rho, psimx,
c    1    psimn,
c    1   aerrmx, rerrmx
      real(kind=dp), dimension(0:kmax - 1) :: dla, ba
      real(kind=dp), dimension(npsi0) :: alenb
      real(kind=dp), dimension(0:kmax) :: thetps
      real(kind=dp), dimension(0:kmax - 1) :: thetpe
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, np, nmax1, imax1, lmaxa, lmaxa1, worksize, pwork,
     1   pua, pug, pth0a, pth0g, pcsth0a, pchi, plegjy, pchil, plegp,
     2   pduug, pfu, pdtt, pfm, ppsid, plam1a, plam2g, pharr, puinv,
     3   pth0inv, presid, presida, ierc, ierd, memsize
      real(kind=dp), dimension(npsi0) :: bmax, bmin, psis
      real(kind=dp) :: dene, teme
      real(kind=dp), allocatable, save, dimension(:) :: memory
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL runb
      real*8 rhopsi
c-----local
      real(kind=dp) :: rho_loc,unorm,zi,c2,eps
C-----------------------------------------------
c
c   Common block for Comdynray
c

C
C
C
C  indices for ua,ug,th0a,th0g,csth0a,chi,legjy,chil,legp,duug,
C               fu,dtt,fm,psid,lam1a,lam2g,harr,uinv,th0inv,resid,
C               resida

C
C  allocate memory
C

c      write(*,*)'adj_sub.f in runa nmax,imax,npsi0',nmax,imax,npsi0

      nmax1 = nmax
      imax1 = imax
      if (mod(nmax1,2) == 1) nmax1 = nmax1 + 1
      if (mod(imax1,2) == 0) imax1 = imax1 + 1

c      write(*,*)'adj_sub.f in runa nmax1,imax1',nmax1,imax1

C  |nmax1| and |imax1| are the allocated dimensions for two-dimensional
C  arrays.These might differ from |nmax| and |imax| on the Crays, where
C   strides are inefficiently handled.
C   To avoid memory conflicts on the Crays
C   Make even for accessing
C              |ijy(0:nmaxa1,0:6,0:1)| by row in |potential|
C    Make odd for accessing
C        |chi(0:imax1-1,0:nmax-1)| by row in |decomp|


C |lmaxa1| is the allocated dimension for the |la| index.  This is
C  usually |lmaxa| but must be at least $1$ to ensure that dimension
C  specifications such as |1:lmaxa1| are legal.  Also $H_{1,1,1}$ is needed
C        even when |lmaxa=0|.


C  Round |lmax| down to an odd number
C  |lmaxa1| needs to be at least 1
C
      if (mod(lmax,2) == 0) lmax = lmax - 1
      lmaxa = (lmax + 1)/2
      lmaxa1 = lmaxa
      if (lmaxa1 == 0) lmaxa1 = 1

C  count memory size
C
      memsize = 0
      worksize = 0
C
      pua = memsize
      memsize = memsize + nmax
C
      pug = memsize
      memsize = memsize + nmax + 1
C
      pth0a = memsize
      memsize = memsize + imax
C
      pth0g = memsize
      memsize = memsize + imax + 1
C
      pcsth0a = memsize
      memsize = memsize + imax
C
      pchi = memsize
      memsize = memsize + imax1*nmax
C
      plegjy = memsize
      memsize = memsize + nmax1*7*4*(lmaxa1 + 1)
C
C Figure the workspace for |initjy| (|jarray| and |djarray|).
C
      worksize = max(worksize,2*(max(2*lmaxa - 1,0) + 2 + 1)*7*2)
C
C   Workspace for |potential| (|fla| and |ijy|)
C
      worksize = max(worksize,nmax + (nmax1 + 1)*7*2)
C
C  indices for chil, legp
C
      pchil = memsize
      memsize = memsize + nmax1*lmaxa1

      plegp = memsize
      memsize = memsize + imax1*lmaxa1
C
C  Workspace for |initp|  (|leg0|).
C
      worksize = max(worksize,imax)
C
C  Indices for duug,fu,dtt,fm,psid

      pduug = memsize
      memsize = memsize + nmax + 1

      pfu = memsize
      memsize = memsize + nmax

      pdtt = memsize
      memsize = memsize + nmax

      pfm = memsize

      memsize = memsize + nmax

      ppsid = memsize
      memsize = memsize + nmax1*5*2
C
C  Workspace for |reactall| (|reactl| and |chila|) and for |potential|
C  (|fla| and |ijy|).
C
      worksize = max(worksize,nmax*2 + nmax + (nmax1 + 1)*7*2)
C
C  Indices of lam1a,lam2g,harr
C
      plam1a = memsize
      memsize = memsize + imax

      plam2g = memsize
      memsize = memsize + imax + 1

      pharr = memsize
      memsize = memsize + (((lmaxa1 - 1)*lmaxa1*(2*lmaxa1 - 1))/6 + (
     1   lmaxa1 - 1)*lmaxa1 + lmaxa1)
C
C  Workspace for |hinit| (|garr|).
C
      worksize = max(worksize,(((lmaxa1 + 3) - 1)*(lmaxa1 + 3))/2 + (
     1   lmaxa1 + 1))
C
C  Workspace for |error| (|aerr| and |rerr|).
C
      worksize = max(worksize,2*imax)
C
C  Indices for uinv,th0inv,resid,resida
C
      puinv = memsize
      memsize = memsize + nmax1*3

      pth0inv = memsize
      memsize = memsize + imax1*3

      presid = memsize
      memsize = memsize + imax1*nmax

      presida = memsize
      memsize = memsize + imax1*nmax
C
C  Workspace for |current| (|cur|).
C
      worksize = max(worksize,imax)
C
C  Sum memsize and worksize
C
      pwork = memsize
      memsize = memsize + worksize

c      write(*,*)'in runa memsize',memsize

C
C  allocate memory
C
c      call hpalloc(memptr,memsize,ierc,1)
      if (.not.allocated(memory)) then
c         ALLOCATE(memory(memsize),stat=ierc)
          ALLOCATE(memory(0:memsize),stat=ierc)
         if (ierc /= 0) then
            print *, 'memory could not be allocated in runa'
            return
         endif
      endif
c      call bcast (memory, memsize)
      memory = 0.

c      write(*,*)'in runa before read thetps iout3',iout3

      read (iout3, 2020) thetps

c      write(*,*)'in runa thetps',thetps

      read (iout3, 2020) thetpe

c      write(*,*)'in runa thetpe',thetpe
     
c      write(*,*)'in runa npsi0',npsi0

      do np = 1, npsi0   
         write(*,*)'runa np',np
         
         read (iout3, 2010) psis(np)

c         write(*,*)'runa psis(np)',psis(np) 
         rho_loc=rhopsi(psis(np))
         write(*,*)'rhopsi',rho_loc

         read (iout3, 2020) dla

c         write(*,*)'runa kmax',kmax
c         write(*,*)'runa dla',dla
c         do n=0,kmax-1
c           write(*,*)'n,dla(n)',n,dla(n)
c         enddo

c         write(*,*)'ba(0:nthpp0-1) !b/b0along b line kmax=nthp0'

         read (iout3, 2020) ba

c         do n=0,kmax-1
c           write(*,*)'n,ba(n)',n,ba(n)
c         enddo

         read (iout3, 2010) alenb(np), bmax(np), bmin(np)

c         write(*,*)'runa np,alenb(np), bmax(np), bmin(np)',
c     &                   np,alenb(np), bmax(np), bmin(np)
         
         call subpar (psis(np), dene, teme, zi)

c         write(*,*)'runa   psis(np), dene, teme, zi',
c     &                     psis(np), dene, teme, zi

         write (iout5, 3000) psis(np), dene, teme
c
c   teme in Kev
c
         c2 = 511.925/teme
c         write(*,*)'c2=511.925/teme ',c2

cSAP070913
         unorm=dsqrt(teme)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
         c2=(2.99792458d10/unorm)**2 !(clight/unorm)**2
c         write(*,*)'unorm,c2',unorm,c2

c         stop 'runa before runb'

c       if(np.ne.2) goto10 !for test

        write(*,*)'in runa before runb'
c        write(*,*)'pua,pug,pth0a,pth0g',pua,pug,pth0a,pth0g
c        write(*,*)'pcsth0a,pchi,plegjy',pcsth0a,pchi,plegjy
c        write(*,*)'pchil,plegp,pduug,pfu',pchil,plegp,pduug,pfu
c        write(*,*)'pdtt,pfm,ppsid,plam1a',pdtt,pfm,ppsid,plam1a
c        write(*,*)'plam2g,pharr,puinv,pth0inv',
c     &             plam2g,pharr,puinv,pth0inv
c        write(*,*)'presid,presida,pwork',presid,presida,pwork

        call runb (t, c2, ze, zi, eps, nmax, umax, imax, kmax, dt,tmax
     1      , alpha, rho, lmax, tskip, nmax1, imax1, lmaxa, lmaxa1,
     2      worksize, memory(pua), memory(pug), memory(pth0a), memory(
     3      pth0g), memory(pcsth0a), memory(pchi), memory(plegjy),
     4      memory(pchil), memory(plegp), memory(pduug), memory(pfu),
     5      memory(pdtt), memory(pfm), memory(ppsid), memory(plam1a),
     6      memory(plam2g), memory(pharr), memory(puinv), memory(pth0inv
     7      ), memory(presid), memory(presida), memory(pwork), npsi0,
     8      dla, ba, alenb, thetps, thetpe, np, bmax, bmin, iout3, iout4
     9      , iout5, iout7, iout8, iswchi, aerrmx, rerrmx,
     &      dens_averaged_ar(np))

c         write(*,*)'runa after runb'
c         write(*,*)'no,dens_averaged_ar(np)',np,dens_averaged_ar(np)
c         stop 'runa after runb'

C
C  Deallocate memory
C

c10        continue !for test
      end do
      close(iout5)

c      call hpdeallc(memptr,ierc,0)
c     if (allocated(memory)) then
c        print *, 'memory still allocated'
c        DEALLOCATE(memory, stat=ierd)
c        print *, 'ierd=', ierd
c     endif
c     print *, 'memory deallocated'
c-----------------------------
 2010 format(4(2x,e20.12))
 2020 format(5(2x,e20.12))
 3000 format(4(2x,e20.12))
c 2010 format(4(2x,e24.16))
c 2020 format(5(2x,e24.16))
c 3000 format(4(2x,e24.16))

      return
      end subroutine runa

      subroutine runb(t, c2, ze, zi, eps, nmax, umax, imax, kmax, dt,
     1   tmax, alpha, rho, lmax, tskip, nmax1, imax1, lmaxa, lmaxa1,
     2   worksize, ua, ug, th0a, th0g, csth0a, chi, legjy, chil, legp,
     3   duug, fu, dtt, fm, psid, lam1a, lam2g, harr, uinv, th0inv,
     4   resid, resida, work, npsi0, dla, ba, alenb, thetps, thetpe, np
     5   , bmax, bmin, iout3, iout4, iout5, iout7, iout8, iswchi, aerrmx
     6   , rerrmx,
     &   dens_averaged)

C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
C
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nmax, imax, kmax, tmax, lmax, tskip, nmax1, imax1, lmaxa,
     1   lmaxa1, worksize, npsi0, np, iout3, iout4, iout5, iout7, iout8
     2   , iswchi
      real(kind=dp) t, c2, ze, zi, eps, umax, dt, alpha, rho,
     1    aerrmx, rerrmx
      real(kind=dp), dimension(0:nmax - 1) :: ua
      real(kind=dp), dimension(0:nmax) :: ug
      real(kind=dp), dimension(0:imax - 1) :: th0a
      real(kind=dp), dimension(0:imax) :: th0g
      real(kind=dp), dimension(0:imax - 1) :: csth0a
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: chi
      real(kind=dp), dimension(0:nmax1 - 1,0:6,0:3,0:lmaxa1) :: legjy
      real(kind=dp), dimension(0:nmax1 - 1,lmaxa1) :: chil
      real(kind=dp), dimension(0:imax1 - 1,lmaxa1) :: legp
      real(kind=dp), dimension(0:nmax) :: duug
      real(kind=dp), dimension(0:nmax - 1) :: fu, dtt, fm
      real(kind=dp), dimension(0:nmax1 - 1,0:4,0:1) :: psid
      real(kind=dp), dimension(0:imax - 1) :: lam1a
      real(kind=dp), dimension(0:imax) :: lam2g
      real(kind=dp), dimension(((lmaxa1 - 1)*lmaxa1
     1  *(2*lmaxa1 - 1))/6+(lmaxa1-1)*lmaxa1+lmaxa1):: harr
      real(kind=dp), dimension(0:nmax1 - 1,-1:1) :: uinv
      real(kind=dp), dimension(0:imax1 - 1,-1:1) :: th0inv
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: resid, resida
      real(kind=dp), dimension(0:worksize - 1) :: work
      real(kind=dp), dimension(0:kmax - 1) :: dla, ba
      real(kind=dp), dimension(npsi0) :: alenb
      real(kind=dp), dimension(0:kmax) :: thetps
      real(kind=dp), dimension(0:kmax - 1) :: thetpe
      real(kind=dp), dimension(npsi0) :: bmax, bmin
      real(kind=dp)  ::  dens_averaged
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: itim, n, i, time,j
      real(kind=dp)::dthp,du,dth0,th0max,pi,c,deltb,abserr,relerr,curv,
     1   cond,e
cSm070914
      integer la 
cSAP080321 for using md03uaf
      integer  ifail,ierc
c      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: work_2d
      real(kind=dp), allocatable, save, dimension(:,:) :: work_2d
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(kind=dp) , EXTERNAL :: deltab, lam1, lam2, current,
     &averaged_density,current1
      EXTERNAL initjy, initp, maxwel, isotrop, hinit, matinit, decomp,
     1   reactall, residue, error, solve
C-----------------------------------------------
C
C  Velocity space meshes
C

C
C   Array holding Legendre functions
C
C
C  Legendre decomposition of chi
C
C
C

C
C
C   Arrays to hold various flux averaged quantities.
C
C
C  Matrices for inversion
C
C
C  Time advancement;  residues
C
C
C   Current and conductivity
C  The conductivity is define as $J_{0\parallel}/E_{0\parallel}$ where
C    $$E_{0\parallel} = B_0\frac{\int E_{\parallel} dl}{\int B\,dl}.$$
C    Substituting $E_\parallel = (B_\zeta/B) \Phi/(2\pi R)$ and
C    $\Phi=TQ/L$ gives $E_\parallel = T/\int b \,dl/L$.
C

c      write(*,*)'runb t,c2,ze,zi,eps,nmax',t,c2,ze,zi,eps,nmax
c      write(*,*)'runb umax,imax,kmax,dt,tmax,alpha,rho',
c     &umax,imax,kmax,dt,tmax,alpha,rho
c      write(*,*)'runb lmax,tskip,nmax1,imax1,lmaxa,lmaxa1',
c     &lmax,tskip,nmax1,imax1,lmaxa,lmaxa1

C
C  Initialize velocity meshes, chi etc.
C
      du = umax/nmax
      do n = 0, nmax - 1
         ua(n) = (n + 0.5d0)*du
         ug(n) = n*du
      end do
      ug(nmax) = umax 

c      write(*,*)'runb ua=',(ua(n),n=0,nmax-1)
c      write(*,*)'runb ug=',(ug(n),n=0,nmax)


      pi = 4*datan(1.0d0)
      deltb = (bmax(np)-bmin(np))/bmin(np)

c      write(*,*)'runb np,bmax(np),bmin(np),deltb',
c     & np,bmax(np),bmin(np),deltb

      th0max = datan2(1.0_dp,dsqrt(deltb))

c      write(*,*)'runb th0max =',th0max

      eps = deltb/(2.d0 + deltb)

c      write(*,*)'runb eps',eps
      
      dth0 = th0max/imax
      do i = 0, imax - 1
         th0a(i) = (i + 0.5d0)*dth0
         csth0a(i) =dcos(th0a(i))
         th0g(i) = i*dth0
      end do
      th0g(imax) = th0max

c      write(*,*)'runb th0a(i)',(th0a(n),n=0,imax-1)
c      write(*,*)'th0g(i)',(th0g(n),n=0,imax-1)
c      write(*,*)'csth0a(i)',(csth0a(n),n=0,imax-1)

C
C  Set up trial chi-function
C
      if (np == 1) then
         if (iswchi == 0) then
            chi(:imax-1,:nmax-1) = 0
         else
            read (iout7, 3030) ((chi(i,n),i=0,imax - 1),n=0,nmax - 1)
         endif
      endif
C
C   Fill the arrays
C
      c = dsqrt(c2)
c      write(*,*)'runb c,c2',c,c2        
      lmaxa = (lmax + 1)/2

c      stop 'runb  before  call initjy'

c      write(*,*)'runb before  call initjy lmaxa',lmaxa

      call initjy (ua, c, nmax, nmax1, lmaxa, lmaxa1, legjy, work(0),
     1   work(2*(max(2*lmaxa-1,0)+2+1)*7), max(2*lmaxa - 1,0) + 2)

c      write(*,*)'legjy',legjy

C
C  Initialize legp
C
      call initp (csth0a, imax, imax1, lmaxa, lmaxa1, legp, work)

c      write(*,*)'legp',legp

C
C  Initialize isotropic terms.
C
 
c     write(*,*)'runb before  call maxwel c,t',c,t

      call maxwel (ua, du, c, nmax, fm, t)

c      write(*,*)'fm',fm

      chil(:nmax-1,1) = 1
     
      call isotrop (ua, du, c, fm, chil(0,1), nmax, nmax1, legjy(0,0,0,0
     1   ), psid, ze, zi, duug, fu, dtt, work)

c      write(*,*)'duug ',duug
c      write(*,*)'fu ',fu
c      write(*,*)'dtt ',dtt
    
     
C  Initialize lam1a,lam2g,harr
C
      dthp = thetpe(1) - thetpe(0)
 
c      write(*,*)'dla',dla
c      write(*,*)'ba',ba

      do i = 0, imax - 1
c         write(*,*)'i,th0a(i),kmax,dthp',i,th0a(i),kmax,dthp
         lam1a(i) = lam1(eps,th0a(i),kmax,dthp,dla(0),ba(0))
c         write(*,*)'runb i,lam1a(i)', i,lam1a(i)
      end do

c      write(*,*)'eps,kmax,dthp',eps,kmax,dthp
c      write(*,*)'dla',dla
c      write(*,*)'ba',ba

      do i = 0, imax
c         write(*,*)'i,th0g(i)',i,th0g(i)
         lam2g(i) = lam2(eps,th0g(i),kmax,dthp,dla(0),ba(0))
c         write(*,*)'runb i,lam2g(i)',i,lam2g(i)
      end do
      
c      do i=0,imax-1
c        write(*,*)'i,lam1a(i),lam2g(i)',i,lam1a(i),lam2g(i)
c      end do

c      write(*,*)'runb before hinit'
c      write(*,*)'eps, kmax, lmaxa1',eps, kmax, lmaxa1
c      write(*,*)'dla',dla
c      write(*,*)'ba',ba
       
c      write(*,*)'runb before dens fm',fm
c      write(*,*)'lam1a',lam1a
c      write(*,*)'ua',ua
c      write(*,*)'th0a',th0a
c      write(*,*)'runb before averaged_density'
      dens_averaged=averaged_density(nmax,imax,ua,th0a,fm,lam1a)

      write(*,*)'in runb dens_averaged',dens_averaged
      call hinit (eps, kmax, lmaxa1, harr, work, dthp, dla(0), ba(0))

c      write(*,*)'harr',harr    

C
C  Initialize matrices for inversion.
C

c      write(*,*)'runb before matinit'

      call matinit (du, ua, ug, nmax, nmax1, dth0, th0a, th0g, imax,
     1   imax1, duug, fu, dtt, lam1a, lam2g, uinv, th0inv, rho)

c      write(*,*)'uinv ',uinv
c      write(*,*)'th0inv',th0inv

C
C  time evolution of solution
C
      itim = 0
      do time = 0, tmax

         write(*,*)'time,itim',time,itim

         if (itim /= 1) then
C
C  Calculate reaction term.
C
 
c      write(*,*)'runb before decomp'

            call decomp (chi, nmax, nmax1, imax, imax1, lmaxa, lmaxa1,
     1         th0a, dth0, legp, chil)

c            write(*,*)'chil',chil 

c           write(*,*)'runb before reactall'

            call reactall (ua, du, c, fm, t, chil, nmax, nmax1, imax,
     1         imax1, lam1a, legjy(0,0,0,1), psid, legp, harr, ze, lmaxa
     2         , lmaxa1, resid, work(0), work(nmax), work(2*nmax))
             
c           write(*,*)'reactl'
c            do i=0,nmax-1
c               write(*,*)'i,work(i)',i,work(i)
c            enddo 

c            write(*,*)'chila'
c            do i=nmax,2*nmax-1
c               write(*,*)'i,work(i)',i,work(i)
c            enddo

C

C
C  Calculate residue and its norm
C
            call residue (ua, nmax, nmax1, csth0a, imax, imax1, chi, dtt
     1         , c, lam1a, uinv, th0inv, resid)
            if (mod(time,tskip) == 0) call error (du, ua, umax, nmax,
     1         dth0, th0a, th0max, imax, imax1, resid, chi, abserr,
     2         relerr, work(0), work(imax))
C
C  Advance chi
C

cSAP080321
             ifail=0 
             if (.not.allocated(work_2d)) then
c                ALLOCATE(memory(memsize),stat=ierc)
                 ALLOCATE(work_2d(0:imax1 - 1,0:nmax - 1),stat=ierc)
                 if (ierc /= 0) then
                   print *, 'memory could not be allocated in runb'
                   return
                 endif
             endif
     
c             write(*,*)'before  md03uaf'
      call md03uaf(imax,nmax,imax1,th0inv(0,-1),th0inv(0,0),th0inv(0,1),
     &uinv(0,-1),uinv(0,0),uinv(0,1),-dt,alpha,time,resid,resida,work_2d
     &,ifail)            
c             write(*,*)'after  md03uaf'

             if(ifail.NE.0)then
               write(6,*)'md03uaf failed: ',ifail,time
      
      
               write(*,*)'runb before solve 1'

            call solve (resida, imax1, nmax, 1, imax, uinv(0,-1), uinv(0
     1         ,0), uinv(0,1), resid, dt, alpha)


            write(*,*)'runb before solve 2'

            call solve (resid, 1, imax, imax1, nmax, th0inv(0,-1),
     1         th0inv(0,0), th0inv(0,1), resida, dt, alpha)
        end if
c            do i=1,imax-1
c              write(*,*)'i,th0inv(i,-1),th0inv(i,0),th0inv(i,1)',
c     &        i,th0inv(i,-1),th0inv(i,0),th0inv(i,1)
c            enddo

c            write(*,*)'resida'
c            do n=1,nmax-1
c               do i=1,imax-1
c                  write(*,*)'n,i,resida(i,n)',n,i,resida(i,n) 
c               enddo
c            enddo

c            write(*,*)'resid'
c             do n=1,nmax-1
c               do i=1,imax-1
c                  write(*,*)'n,i,resid(i,n)',n,i,resid(i,n) 
c               enddo
c            enddo

c            write(*,*)'1 chi'
c            do n=1,nmax-1
c               do i=1,imax-1
c                  write(*,*)'n,i,chi(i,n)',n,i,chi(i,n) 
c               enddo
c            enddo

            chi(:imax-1,:nmax-1) = chi(:imax-1,:nmax-1) + dt*resid(:imax
     1         -1,:nmax-1)

c            write(*,*)'2 chi'
c            do n=1,nmax-1
c               do i=1,imax-1
c                  write(*,*)'n,i,chi(i,n)',n,i,chi(i,n) 
c               enddo
c            enddo

c            stop 'runb after solve'

C
C  compute current and conductivity.
C

c            write(*,*)'mod(time,tskip)',mod(time,tskip)
           
            if (mod(time,tskip) == 0) then
               curv = current(du,ua,c,nmax,dth0,th0a,imax,imax1,fm,chi,
     1            work)
               e = t/harr(((1-1)*1*(2*1-1))/6+(1-1)*1+1)

               write(*,*)'t,harr(((1-1)*1*(2*1-1))/6+(1-1)*1+1)',
     &         t,harr(((1-1)*1*(2*1-1))/6+(1-1)*1+1)

               cond = curv*zi/(t*dsqrt(t)*e)
               write (iout4, 4020) time, cond, abserr, relerr
 4020          format(2x,i5,3x,e20.12,2x,e16.8,2x,e16.8)
 

            endif
c            write(*,*)'curv,zi,t,dsqrt(t),e',curv,zi,t,dsqrt(t),e
            write(*,*)'curv,cond',curv,cond
            write(*,*)'abserr,aerrmx',abserr,aerrmx
            write(*,*)'relerr,rerrmx',relerr,rerrmx

            if (abserr<aerrmx .and. relerr<rerrmx) itim = 1
         endif
      end do
C
C  Write out results.
C
      write (iout5, 3010) nmax, imax, kmax, lmax
      write (iout5, 3020) t, c2, ze, zi, eps
      write(iout5,3020)umax,th0max,du,dth0,rho,deltb,bmin(np),bmax(np)
      write (iout5, 3040) dt, tmax, alpha
      write (iout5, 3020) cond, abserr, relerr
      write (iout5, 3030) ua
      write (iout5, 3030) ug
      write (iout5, 3030) th0a
      write (iout5, 3030) th0g
      write (iout5, 3030) lam1a
      write (iout5, 3030) fm
      write (iout5, 3030) ((chi(i,n),i=0,imax - 1),n=0,nmax - 1)
      if(np==1)write(iout8,3030)((chi(i,n),i=0,imax-1),n=0,nmax-1)

c      write(*,*)'runb end chi'
c       do n=0,nmax-1
c         do i=0,imax-1
c            write(*,*)'n,i,chi(i,n)',n,i,chi(i,n)
c         enddo
c       enddo
        curv = current1(du,ua,c,nmax,dth0,th0a,imax,imax1,fm,chi,
     1            work)
               e = t/harr(((1-1)*1*(2*1-1))/6+(1-1)*1+1)

               write(*,*)'t,harr(((1-1)*1*(2*1-1))/6+(1-1)*1+1)',
     &         t,harr(((1-1)*1*(2*1-1))/6+(1-1)*1+1)

               cond = curv*zi/(t*dsqrt(t)*e)
        write(*,*)'end runb'
        write(*,*)'np',np
c        write(*,*)'chi(0,130)',chi(0,130)
c        write(*,*)'du,c,nmax,dth0,imax,imax1',du,c,nmax,dth0,imax,imax1
c        write(*,*)'ua',ua
c        write(*,*)'th0a',th0a
c        write(*,*)'fm',fm 
c        write(*,*)'chi',chi 
        write(*,*)'end curv,zi,t,dsqrt(t),e',curv,zi,t,dsqrt(t),e
        write(*,*)'curv,cond',curv,cond
c       stop 'end runb'
c-----------------------------
c
 3010 format(1x,4i8)

 3020 format(4(1x,e20.12))
 3030 format(5(1x,e20.12))
 3040 format(1x,e20.12,2x,i8,2x,e20.12)
c 3020 format(4(1x,e24.16))
c 3030 format(5(1x,e24.16))
c 3040 format(1x,e24.16,2x,i8,2x,e20.16)

C
      return
      end subroutine runb


      subroutine solve(x, ns, n, ms, m, a, b, c, z, dt, alpha)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
C
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ns, n, ms, m
      real(kind=dp) dt, alpha
      real(kind=dp), dimension(0:ms - 1,0:*) :: x
      real(kind=dp), dimension(0:n - 1) :: a, b, c
      real(kind=dp), dimension(0:ms - 1,0:*) :: z
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, k
      real(kind=dp) :: dt2, den
C-----------------------------------------------

      dt2 = -alpha*dt

      i = 0
      den = 1 + dt2*b(i)

cSAP080320
      write(*,*)'in solve i=0 ns,n,ms,m',ns,n,ms,m
      do k=1,m
c         write(*,*)'k,z(ns*i,k-1)',k,z(ns*i,k-1)
         x(ns*i,k-1) = z(ns*i,k-1)/den
         z(ns*i,k-1) = -dt2*c(i)/den
c         write(*,*)'x(ns*i,k-1)',x(ns*i,k-1)
      enddo

c      x(ns*i,:m-1) = z(ns*i,:m-1)/den
c      z(ns*i,:m-1) = -dt2*c(i)/den

      do i = 1, n - 2
         do k = 0, m - 1
            WRITE(*,*)'i,k,ns*i',i,k,ns*i
            den = 1 + dt2*(b(i)+a(i)*z(ns*(i-1),k))
            x(ns*i,k) = (z(ns*i,k)-dt2*a(i)*x(ns*(i-1),k))/den
            z(ns*i,k) = -dt2*c(i)/den
         end do
      end do

      i = n - 1
      do k = 0, m - 1
         den = 1 + dt2*(b(i)+a(i)*z(ns*(i-1),k))
         x(ns*i,k) = (z(ns*i,k)-dt2*a(i)*x(ns*(i-1),k))/den
         z(ns*i,k) = 0
      end do





      do i = n - 2, 0, -1
         x(ns*i,:m-1) = z(ns*i,:m-1)*x(ns*(i+1),:m-1) + x(ns*i,:m-1)
      end do

      return
      end subroutine solve

      subroutine subadj(nmax,umax,imax,nthp0,dt,tmax,alpha,ze,t,
     1   rho, lmax, iout3, iout5, iout4, iout7, iout8, iswchi, aerrmx,
     2   rerrmx,
     &   npsi0,dens_averaged_ar)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
C
C  calculates the Spitzer function chi(theta,ua) and output to adjout
C
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer
     & nmax,    !number of momentum mesh points :ua
     & imax,    !number of pitch angle mesh points: th0a
     & nthp0,   !number of integration points for field line integration.
                ! ~ 1000 
     & tmax,    !is the maximal number of time steps
     & lmax,    !is the maximal number of Legandre harmonic
                !lmax=1 corresponds to collision with a Maxwellian.
                !          Usually lmax=21
     & iout3,   !number of input file='adjinp',
     & iout5,   !number of output file='adjout',
     & iout4,   !previously was used to write cpu time file='scrach',
     & iout7,
     & iout8,   !if(np=1)write(iout8,3030)((chi(i,n),i=0,imax-1),n=0,nmax-1)
     & iswchi,   !=0  chi(:imax-1,:nmax-1) = 0
                !else
                !read(iout7,3030) ((chi(i,n),i=0,imax-1),n=0,nmax-1)
     & npsi0    !number of radial points for adj function calculations 
      real(kind=dp) 
     & umax,    ! is the maximal value of u_0. 
                ! The momentum grid is spaced between 0 and umax.
     & dt,      ! is the time step
     & alpha,   ! alpha is a weighting factor in the time stepping algorithm.
                !       =0 for explicit algorithm (usually =0.55)
     & rho,     ! governs the treatment  of the boundary at u=umax
                ! for rho=1 the second derivative from the adjoint function chi
                ! d^2(chi)/dx^**2(at u=umax)=0 
                ! for rho=0 the first derivative =0 
     & aerrmx,
     & rerrmx,
     & ze,      !ze is the electron charge state. It is usually =1, 
                !but for Lorentz limit, ion charge Zi ==> infinity,
                ! set  ze=0 and zi=1
     & t,        ! the temperature in units of m*u_n**2
                ! If you want velocities normalized to sqrt(T/m), set t=1.
c-----output
     &dens_averaged_ar (1:npsi0) !bounce averaged density from the maxwellian distribution

c      integer nmax, imax, nthp0, tmax, lmax, iout3, iout5, iout4, iout7
c     1   , iout8, iswchi
c      real(kind=dp) umax, dt, alpha, rho, aerrmx, rerrmx
C----------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer::tskip, dlaptr, baptr, alnptr, ndum
c      integer::npsi0
      integer::thsptr, theptr, ii, iercd

      real(kind=dp) :: psimx, psimn, c2, cpu, cpu0
      real(kind=dp),allocatable,save, dimension(:) :: dum


      real(kind=dp) rho_in,temp_kev,unorm
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL runa

      real*8 temperho,psi_rho
C-----------------------------------------------
c
c   Common block for Comdynray
c
c     pointer(dumptr,dum(1))
C
C
 
      write(*,*)'adj_sub.f in subadj input arguments'
      write(*,*)'nmax, umax, imax, nthp0, dt,ze',
     &nmax, umax, imax, nthp0, dt,ze
      write(*,*)'tmax,alpha,rho,lmax,aerrmx,rerrmx,iswchi',
     &tmax,alpha,rho,lmax,aerrmx,rerrmx,iswchi
      write(*,*)'iout3, iout5, iout4, iout7, iout8',
     &iout3, iout5, iout4, iout7, iout8
     
      write(*,*)'adj_sub.f  subadj  before read 2000 iout3=',iout3
      write(*,*)'adj_sub npsi0',npsi0
cSAP070727
      open(iout3, file='adjinp', status='old')

c      write(*,*)'adj_sub.f  subadj after open  iout3'

      read (iout3, 2000) nthp0, npsi0

c      write(*,*)'subadj  after read nthp0, npsi0',nthp0, npsi0

c      write(*,*)'subadj  before read 2010 iout3=',iout3

      read (iout3, 2010) psimx, psimn

c      write(*,*)'subadj  after read psimx, psimn', psimx, psimn

      ii = 1
      dlaptr = ii
      ii = ii + nthp0
      baptr = ii
      ii = ii + nthp0
      alnptr = ii
      ii = ii + npsi0
      thsptr = ii
      ii = ii + nthp0 + 1
      theptr = ii
      ii = ii + nthp0
      ndum = ii

c      write(*,*)'adj_sub.f in subadj'
c      write(*,*)'dlaptr,baptr,alnptr,thsptr,theptr,ndum',
c     &           dlaptr,baptr,alnptr,thsptr,theptr,ndum

c      call hpalloc(dumptr,ndum,iercd,1)
      if (.not.allocated(dum)) then
         ALLOCATE(dum(ndum),stat=iercd)
         if (iercd /= 0) then
            print *, 'dum could not be allocated in subadj'
            return
         endif
      endif
      write (iout5, 2000) nthp0, npsi0
      write (iout5, 2010) psimx, psimn
   
      tskip = 10    !time to print out current and errors
!      call second (cpu0)
      
      write(*,*)'adj_sub.f in subadj before runa nmax,umax,imax,npsi0',
     &nmax, umax, imax,npsi0
      write(*,*)'subadj iout3,iout4,iout5,iout7,iout8',
     &iout3,iout4,iout5,iout7,iout8

      call runa (t, ze, nmax, umax, imax, nthp0, dt, tmax,
     1   alpha, rho, lmax, tskip, npsi0, dum(dlaptr), dum(baptr), dum(
     2   alnptr), dum(thsptr), dum(theptr), psimx, psimn, iout3, iout4,
     3   iout5, iout7, iout8, iswchi, aerrmx, rerrmx,
     &   dens_averaged_ar)

       write(*,*)'after runa dens_averaged_ar', dens_averaged_ar

!      call second (cpu)
  
c      write(*,*)'adj_sub.f in subadj after runa'

!      write (iout4, 4030) cpu - cpu0
c      call hpdeallc(dumptr,iercd,0)
      DEALLOCATE(dum)
 4030 format(1x,'adj CPU time: ',f8.4,' secs.')
c----------------------
 2000 format(1x,2i8)
 2010 format(4(2x,e20.12))

c 2010 format(4(2x,e24.16))
      return
      end subroutine subadj

      subroutine xtrema(imx, xmx, ymx, dymx, nwk, wkt, wk, wk1, ier)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
C
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer imx, nwk, ier
      real(kind=dp) xmx, ymx, dymx
      real(kind=dp), dimension(nwk) :: wkt, wk
      real(kind=dp), dimension(nwk,3) :: wk1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: it
      real(kind=dp), dimension(1) :: x0, y0, dy0, x1, y1, dy1
      real(kind=dp) :: h, dx
C-----------------------------------------------
      x0(1) = wkt(imx)
      y0(1) = wk(imx)
      h = wkt(2) - wkt(1)
      call mcsevu1 (h, wkt(1), wk, nwk, wk1, nwk, 1, dy0, x0, ier)
      dx = 0.5*h
      it = 0
      ier = 0
   50 continue
      x1(1) = x0(1) + dx
      call mcsevu (h, wkt(1), wk, nwk, wk1, nwk, 1, y1, x1, ier)
      call mcsevu1 (h, wkt(1), wk, nwk, wk1, nwk, 1, dy1, x1, ier)
      dx = (x1(1)-x0(1))*dy1(1)/(dy0(1)-dy1(1))
      x0(1) = x1(1)
      y0(1) = y1(1)
      dy0(1) = dy1(1)
      it = it + 1
      if (it > 50) then
         ier = 19
         return
      endif
      if (abs(dx)<1.E-8*h .or. abs(dy0(1))<1.E-8) go to 100
      go to 50
  100 continue
      xmx = x0(1)
      ymx = y0(1)
      dymx = dy0(1)
      return
      end subroutine xtrema


      subroutine subpar(psival, dene, teme, zef)
c-----given psi-value, evaluate densities and temperatures,and zeffective
c     inside plasma: small radius =< 1

      implicit none
        
c-----input
      real*8 
     & psival
   
c-----output
      real*8
     & dene, teme, zef


c-----local
      real*8 rholoc

c-----external
      real*8  rhopsi,densrho,temperho,zeffrho
       
     
      rholoc=rhopsi(psival)

      dene=densrho(rholoc,1) !10**13_cm**-3

      teme=temperho(rholoc,1)!KeV

      zef=zeffrho(rholoc)

      return
      end



      SUBROUTINE MD03UAF(N1,N2,N1M,ax,bx,cx,ay,by,cy,mult,APARAM,IT,R,
     $    WRKSP1,WRKSP2,IFAIL)
      implicit none 
cSm020909 chage real to real*8 
      real*8  APARAM,mult
      INTEGER           IFAIL, IT, N1, N1M, N2
 
      real*8  ax(n1),bx(n1),cx(n1),ay(n2),by(n2),cy(n2),
     *                  R(N1M,N2), WRKSP1(N1M,N2),
     *                  WRKSP2(N1M,N2)
 
      real*8  ALM, ALPHA, CS, RIBJL, RIJB, SB, SBSEIJ, SC,
     *                  SCSFIJ, SD, SEIBJL, SEIJB, SFIBJL, SFIJB, XKS1,
     *                  XKS2
      INTEGER           I, IB, IERROR, IL, IS, J, JB, JL, KS, KS1, KS2,
     *                  N1M1, N1P1, N2M1, N2P1
 
      real*8  ALP(9)
 
      INTRINSIC         MOD
 
      DATA              ALP(1), ALP(2), ALP(3), ALP(4), ALP(5), ALP(6),
     *      ALP(7), ALP(8), ALP(9)/1.0d0 , 0.625d0 , 0.25d0 ,
     *      0.875d0 , 0.5d0 , 0.125d0 , 0.75d0 ,
     *      0.375d0 , 0.0d0 /
 

cSAP080321
c      write(*,*)'in   md03uaf'
c      write(*,*)'n1,n2.n1m',n1,n2,n1m

      IERROR = 0

      IF (N1.LE.1) GO TO 220
      IF (N2.LE.1) GO TO 220

      IF (N1M.LT.N1) GO TO 240 

      N1M1 = N1 - 1
      N2M1 = N2 - 1
      N1P1 = N1 + 1
      N2P1 = N2 + 1

      KS = MOD(IT-1,2) + 1

      IF (KS.EQ.0) KS = 2

c      write(*,*)'n1m1,n2m1,n1p1,n2p1,ks',n1m1,n2m1,n1p1,n2p1,ks

      KS1 = 2 - KS
      XKS1 = KS1
      KS2 = KS - 1
      XKS2 = KS2

      IS = MOD(IT-1,18)
      IF (IS.LT.0) IS = IS + 18
      IS = IS/2 + 1

      ALM = 2.d0 *APARAM/((N1M1*N1M1+N2M1*N2M1))

c      write(*,*)'ks1,xks1,ks2,xks2,is,alm',ks1,xks1,ks2,xks2,is,alm       

      IF (APARAM.LE.0.0d0 ) GO TO 260

      IF (ALM.GT.1.0d0 ) GO TO 280

      ALPHA = 1.d0  - ALM**ALP(IS)

c      write(*,*)'alp',alp
c      write(*,*)'alp(is)',alp(is) 
c      write(*,*)'alpha',alpha

      DO 160 JL = 1, N2

c        write(*,*)'n2,jl',n2,jl

         J = KS1*JL + KS2*(N2P1-JL)
         
         JB = J - KS1 + KS2

c        write(*,*)'J,lb',j,jb

         DO 140 I = 1, N1

c           write(*,*)'i',i

            IB = I - 1
            CS = (1.0d0 +mult*(bx(I)+by(J))) 

c            write(*,*)'IB,bx(i),by(j),CS', IB,bx(i),by(j),CS           
c            write(*,*)'jb,n2p1',jb,n2p1

            IF ((JB.EQ.0) .OR. (JB.EQ.N2P1)) GO TO 20
c            write(*,*)'i,jb',i,jb
          
c            write(*,*)'i,jb,WRKSP1(I,JB)',i,jb,WRKSP1(I,JB)
            SEIJB = WRKSP1(I,JB)
c            write(*,*)'SEIJB',SEIJB
c            write(*,*)'i.jb,WRKSP2(I,JB)',i,jb,WRKSP2(I,JB)
            SFIJB = WRKSP2(I,JB)
c            write(*,*)'SFIJB',SFIJB

c            write(*,*)'i.jb,R(I,JB)',i,jb,R(I,JB)
            RIJB = R(I,JB)

c            write(*,*)'SEIJB,SFIJB,RIJB',SEIJB,SFIJB,RIJB

            GO TO 40

   20       CONTINUE

            SEIJB = 0.0d0 
            SFIJB = 0.0d0 
            RIJB = 0.0d0 
c            write(*,*)'after 20 SEIJB SFIJB RIJB,I',SEIJB,SFIJB,RIJB,I
   40       CONTINUE

            IF (I.EQ.1) GO TO 60

c            write(*,*)'after 20 IB,J',IB,J
            SEIBJL = WRKSP1(IB,J)
c            write(*,*)'after 20 IB,J,WRKSP1(IB,J)',IB,J,WRKSP1(IB,J)
            SFIBJL = WRKSP2(IB,J)
c            write(*,*)'after 20 WRKSP2(IB,J)',WRKSP2(IB,J)
            RIBJL = R(IB,J)
c            write(*,*)'after 20 R(IB,J)',R(IB,J)
            GO TO 80

   60       CONTINUE

            SEIBJL = 0.0d0 
            SFIBJL = 0.0d0 
            RIBJL = 0.0d0 
   80       CONTINUE

c            write(*,*)'after 80 XKS1,mult,XKS2,ALPHA,SEIJB',
c     &                          XKS1,mult,XKS2,ALPHA,SEIJB

c            write(*,*)'after 80 J ay(J),cy(J)', J,ay(J),cy(J)

            SB = (XKS1*(mult*ay(J)) +XKS2*(mult*cy(J)) )/
     *              (1.d0 +ALPHA*SEIJB)
          
            SC = (mult*ax(I)) /(1.d0 +ALPHA*SFIBJL)

c            write(*,*)'SB,SC',SB,SC

            SBSEIJ = SB*SEIJB
            SCSFIJ = SC*SFIBJL

            SD = 1.d0 /(-SB*SFIJB-SC*SEIBJL+CS+ALPHA*(SBSEIJ+SCSFIJ))

c            write(*,*)'i,j,cx(I)',i,j,cx(I)
c            write(*,*)'WRKSP1(I,J)', WRKSP1(I,J)
            WRKSP1(I,J) = ((mult*cx(I)) -ALPHA*SBSEIJ)*SD
c            write(*,*)'WRKSP1(I,J)', WRKSP1(I,J)

c            write(*,*)'WRKSP2(I,J)', WRKSP2(I,J)
            WRKSP2(I,J) = (XKS1*(mult*cy(J)) +XKS2*(mult*ay(J)) -
     *           ALPHA*SCSFIJ)*SD

            R(I,J) = (R(I,J)-SB*RIJB-SC*RIBJL)*SD

            GO TO 120

  100       CONTINUE

            WRKSP1(I,J) = 0.0d0 
            WRKSP2(I,J) = 0.0d0 

  120       CONTINUE

  140    CONTINUE
  160 CONTINUE

      DO 200 JL = 1, N2

         J = KS1*(N2P1-JL) + KS2*JL
         JB = J + KS1 - KS2

         DO 180 IL = 1, N1

            I = N1P1 - IL
            IF ((JB.NE.0) .AND. (JB.NE.N2P1)) R(I,J) = R(I,J) -
     *          WRKSP2(I,J)*R(I,JB)

            IF (I.NE.N1) R(I,J) = R(I,J) - WRKSP1(I,J)*R(I+1,J)

  180    CONTINUE
  200 CONTINUE

      IFAIL = 0
      RETURN 

  220 CONTINUE
      IERROR = 1
      GO TO 300

  240 CONTINUE
      IERROR = 2
      GO TO 300

  260 CONTINUE
      IERROR = 3
      GO TO 300

  280 CONTINUE
      IERROR = 4

  300 CONTINUE

      RETURN

      END
