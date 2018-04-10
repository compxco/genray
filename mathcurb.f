      subroutine bsplcof (m, x, y, t, c, mp, wk1, wk2, wkr, ipvt, ifail)
c
      implicit none
c
      character what_id*45
      save      what_id
      data      what_id/'@(#)math.f         11 Apr 95  General Atomics'/
c
      integer m,mp,ifail
      real x(mp),y(mp),t(mp+4),c(mp)
      real wk1(mp,*),wk2(mp,*),wkr(mp)
      integer ipvt(mp)
c
c     uses subroutines 'bsplvb', 'seigm'
c     returns b-spline coefficients for (x(i),y(i)),i=1,m
c     wk1,wk2,wkr, and ipvt are working arrays
c     on a successful exit,ifail will be 0.
c     see M.G. Cox, J. Inst. Math Appl.[15] 95, (1975), also
c     C. de Boor 'A Practical guide to Splines' Chapter 13.
c ------------------------------------------------------------------------------
c
      integer n,nm2,np1
      parameter (n=4,nm2=n-2,np1=n+1)    ! n=4 for cubic spline
      integer    stderr
      parameter (stderr=6)
      real tol
      parameter (tol=1.e-10)
      integer i,j,info,imsg
      real bi(n),relerr
      character errmsg(2)*80
      save errmsg
      data errmsg(1)
     .   /'bsplcof:  singular matrix in determining c'/
      data errmsg(2)
     .   /'bsplcof:  iterative improvement fails in determining c'/
c
      do i=1,n     ! set the knots at ends
         t(i)=x(1)
         t(i+m)=x(m)
      end do
      do i=np1,m
         t(i)=x(i-n+2)
      end do
c
      do j=1,m
         do i=1,m
            wk1(i,j)=0.
         end do
      end do
      wk1(1,1)=1.
      call bsplvb(t,n,1,x(2),n,bi)
      do j=1,n
         wk1(2,j)=bi(j)
      end do
      do i=3,m-2
         call bsplvb(t,n,1,x(i),i+2,bi)
         do j=1,n
            wk1(i,i+j-nm2)=bi(j)
         end do
      end do
      call bsplvb(t,n,1,x(m-1),m,bi)
      do j=1,n
         wk1(m-1,m-n+j)=bi(j)
      end do
      wk1(m,m)=1.
c
      call sgeim1 (wk1,m,mp,y,c,relerr,tol,wk2,wkr,ipvt,info)
c
      if (info .gt. 0) then
         imsg=1
      else if (info .lt. 0) then
         imsg=2
      else
         ifail=info
         return
      end if
      write (stderr,'(a)')errmsg(imsg)
      write (stderr,'(a,i5,e15.6)')'info, relerr = ',info,relerr
      write (stderr,'(a)')'input knots (i,x,y):'
      write (stderr,'(i5,2e15.6)')(i,x(i),y(i),i=1,m)
      if (ifail .eq. 0)
     .  call STOP ('subroutine BSPLCOF: unfinished coding', 201)
      ifail=info
      return
c
      end

      subroutine bsplint(t,bcoef,n,k,x,bvalue,jcal)
c
      implicit none
c
      real t(*),bcoef(*),x,bvalue(*)
      integer n,k,jcal
c     uses subroutine 'interv'
c     returns value and the derivatives of b-spline rep. at x
c     jcal=1, value alone
c     jcal>1, value plus 1st, (jcal-1) derivatives.
c ------------------------------------------------------------------------------
      integer kmax
      parameter (kmax=20)
      integer i,mflag,km1,jcmin,imk,j,jcmax,nmi,kmj,ilo,jj,jc
      integer jderiv
      real dl(kmax),dr(kmax),aj(kmax),ajs(kmax),fkmj
c
      if (jcal .lt. 1)jcal=1
c
      call interv(t,n+k,x,i,mflag)
      if (mflag .lt. 0 .or. mflag .gt. 1)  return
      km1=k-1
      if (km1 .eq. 0) then
         bvalue(1)=bcoef(i)
         do j=2,jcal-1
            bvalue(j)=0.
         end do
         return
      end if
c
      jcmin=1
      imk=i-k
      if (imk .lt. 0) then
         jcmin=1-imk
         do j=1,i
            dl(j)=x-t(i+1-j)
         end do
         do j=i,km1
            aj(k-j)=0.
            dl(j)=dl(i)
         end do
      else
         do j=1,km1
            dl(j)=x-t(i+1-j)
         end do
      end if
c
      jcmax=k
      nmi=n-i
      if (nmi .lt. 0) then
         jcmax=k+nmi
         do j=1,jcmax
            dr(j)=t(i+j)-x
         end do
         do j=jcmax,km1
            aj(j+1)=0.
            dr(j)=dr(jcmax)
         end do
      else
         do j=1,km1
            dr(j)=t(i+j)-x
         end do
      end if
      do jc=jcmin,jcmax
         aj(jc)=bcoef(imk+jc)
      end do
c
      do j=1,k
         ajs(j)=aj(j)
      end do
c
      jderiv=0
 101  continue
      if (jderiv .ge. min(jcal,k))  return
      if (jderiv .ne. 0) then
         do j=1,jderiv
            kmj=k-j
            fkmj=float(kmj)
            ilo=kmj
            do jj=1,kmj
               aj(jj)=((aj(jj+1)-aj(jj))/(dl(ilo)+dr(jj)))*fkmj
               ilo=ilo-1
            end do
         end do
      end if
      if (jderiv .ne. km1) then
         do j=jderiv+1,km1
            kmj=k-j
            ilo=kmj
            do jj=1,kmj
               aj(jj)=(aj(jj+1)*dl(ilo)+aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
               ilo=ilo-1
            end do
         end do
      end if
      bvalue(jderiv+1)=aj(1)
      do j=1,k
         aj(j)=ajs(j)
      end do
      jderiv=jderiv+1
      go to 101
c
      end

      subroutine bsplvb(t,jhigh,index,x,left,biatx)
c
      implicit none
c
      integer jhigh,index,left
      real t(*),x,biatx(*)
c     returns nonvanishing b-splines at x of order jhigh when index=1.
c     t(*) is a nondecreasing knot sequence.
c     t(left) .le. x .le. t(left+1)
c     For the case of index=2 see the reference;
c     C. de Boor's 'A practical Guide to Splines', p134.
c ------------------------------------------------------------------------------
      integer jmax
      parameter (jmax=20)
      integer stderr
      parameter (stderr=6)
      integer i,j,jp1
      real deltar(jmax),deltal(jmax),s,term
      save j,deltar,deltal
      data j/1/
c
      if (x .lt. t(left) .or. x .gt. t(left+1)) then
         write (stderr,'(a)')'bsplvs:  x < t(left) .or. t > t(left+1)'
         write (stderr,'(a)')'x,t(left),t(left+1),left:'
         write (stderr,*)x,t(left),t(left+1),left
         call STOP ('subroutine BSPLVB: unfinished coding', 202)
      end if
c
      if (index .eq. 1) then
         j=1
         biatx(1)=1.
         if (j .eq. jhigh)  return
      end if
c
 20   continue
      jp1=j+1
      deltar(j)=t(left+j)-x
      deltal(j)=x-t(left+1-j)
      s=0.
      do i=1,j
         term=biatx(i)/(deltar(i)+deltal(jp1-i))
         biatx(i)=s+deltar(i)*term
         s=deltal(jp1-i)*term
      end do
      biatx(jp1)=s
      j=jp1
      if (j .lt. jhigh)  go to 20
c
      return
c
      end

      function carlsnc(x,y,ifail)
c
      implicit none
c
      real carlsnc,x,y
      integer ifail
c     returns the carlson elliptic integral: Rc(x,y)
c     x >=0 and y .not.= 0
c ------------------------------------------------------------------------------
c
      include 'rrange.i'    ! define rsmall, rbig, srsmall
c
      real tiny,big
      parameter (tiny=5.*rsmall, big=0.2*rbig)
      real srtiny,tnbg,comp1,comp2
      parameter (srtiny=2.236*srsmall,tnbg=tiny*big,comp1=2.236/srtiny,
     .     comp2=tnbg**2/25.)
      real third
      parameter (third=1./3.)
      real errtol
      parameter (errtol=0.0012)
      real c1,c2,c3,c4
      parameter (c1=0.3 ,c2=1./7., c3=3./8., c4=9./22.)
      real xt,yt,mu,lamda
      real w,s
      integer stderr,imsg
      parameter (stderr=6)
      character errmsg(3)*80
      save errmsg
      data errmsg(1)/'carlsnc:  x < 0. .or. y = 0.'/
      data errmsg(2)/'carlsnc:  x+abs(y) too small .or. too big'/
      data errmsg(3)/'carlsnc:  y is too negative .or. x is too small'/
c
      if (x .lt. 0. .or. y .eq. 0.) then
         imsg=1
         go to 9999
      end if
      if ((x+abs(y)) .lt. tiny .or. (x+abs(y)) .gt. big) then
         imsg=2
         go to 9999
      end if
      if (y. lt. -comp1 .and. x .gt. 0. .and. x .lt. comp2) then
         imsg=3
         go to 9999
      end if
c
      if (y .gt. 0.) then
         xt=x
         yt=y
         w=1.
      else
         xt=x-y
         yt=-y
         w = SQRT (x)/SQRT (xt)
      end if
c
 101  continue
      lamda=yt+2.*SQRT (xt*yt)
      xt=0.25*(xt+lamda)
      yt=0.25*(yt+lamda)
      mu=third*(xt+2.*yt)
      s=(yt-mu)/mu
      if (abs(s).gt. errtol)  go to 101
c
      carlsnc=w*(1.+s**2*(c1+s*(c2+s*(c3+s*c4))))/SQRT (mu)
      ifail=0
      return
c
 9999 continue
      write (stderr,'(a)')errmsg(imsg)
      write (stderr,'(a)')'x, y:'
      write (stderr,*)x,y
      if (ifail .eq. 0)
     .  call STOP ('function CARLSNC: unfinished coding', 203)
      ifail=imsg
      return
c
      end

      function carlsnd(x,y,z,ifail)
c
      implicit none
c
      real carlsnd,x,y,z
      integer ifail
c     returns the carlson's elliptical integeral Rd(x,y,z)
c     x,y >=0 and at most one is zero, z>0
c ------------------------------------------------------------------------------
c
      include 'rrange.i'    ! define rsmall, rbig, ttrbig, ttrsmall
c
      real errtol
      parameter (errtol=0.0015)
      real tiny,big
      parameter (tiny=2./ttrbig, big=0.1*errtol/ttrsmall)
      real c1,c2,c3,c4,c5,c6
      parameter (c1=3./7. ,c2=1./3., c3=3./22., c4=3./11., c5=3./13.,
     .     c6=3./13.)
      integer stderr
      parameter (stderr=6)
      real xt,yt,zt,mu,lamda
      real srz,dx,dy,dz
      real sum,fac
      real s2,s3,s4,s5
      integer imsg
      character errmsg(3)*80
      save errmsg
      data errmsg(1)/'carlsnd:  x .or. y < 0 .or. z <= 0'/
      data errmsg(2)/'carlsnd:  input arguments underflow'/
      data errmsg(3)/'carlsnd:  input arguments overflow'/
c
      if (MIN (x,y) .lt. 0. .or. z .le. 0.) then
         imsg=1
         go to 9999
      end if
      if (MIN (x+y,z) .lt. tiny) then
         imsg=2
         go to 9999
      end if
      if (MAX (x,y,z) .gt. big) then
         imsg=3
         go to 9999
      end if
c
      xt=x
      yt=y
      zt=z
      sum=0.
      fac=1.
 101  continue
      srz = SQRT (zt)
      lamda = SQRT (xt*yt)+(SQRT (yt)+SQRT (xt))*srz
      sum=sum+fac/((zt+lamda)*srz)
      fac=fac*0.25
      xt=0.25*(xt+lamda)
      yt=0.25*(yt+lamda)
      zt=0.25*(zt+lamda)
      mu=0.2*(xt+yt+3.*zt)
      dx=(mu-xt)/mu
      dy=(mu-yt)/mu
      dz=(mu-zt)/mu
      if (MAX (abs(dx),abs(dy),abs(dz)) .gt. errtol)  go to 101
      s2=(dx**2+dy**2+3.*dz**2)/4.
      s3=(dx**3+dy**3+3.*dz**3)/6.
      s4=(dx**4+dy**4+3.*dz**4)/8.
      s5=(dx**5+dy**5+3.*dz**5)/10.
      carlsnd=3.*sum+fac*(1.+s2*(c1+c3*s2+c5*s3)+c2*s3+c4*s4+c6*s5)/
     .     (SQRT (mu)*mu)
      ifail=0
      return
c
 9999 continue
      write (stderr,'(a)')errmsg(imsg)
      write (stderr,'(a)')'x, y, z:'
      write (stderr,*)x,y,z
      if (ifail .eq. 0)
     .  call STOP ('function CARLSND: unfinished coding', 204)
      ifail=imsg
      return
c
      end

      function carlsnf(x,y,z,ifail)
c
      implicit none
c
      real carlsnf,x,y,z
      integer ifail
c     returns the carlson elliptical integeral: Rf(x,y,z)
c     x,y,z >=0 and at most one is zero.
c ------------------------------------------------------------------------------
c
      include 'rrange.i'    ! define rsmall, rbig
c
      real tiny,big,errtol
      parameter (tiny=5.*rsmall, big=0.2*rbig, errtol=0.0025)
      real c1,c2,c3,c4,third
      parameter (c1=1./24. ,c2=0.1, c3=3./44., c4=1./14., third=1./3.)
      integer stderr
      parameter (stderr=6)
      real xt,yt,zt,mu,lamda,xtmp
      real dx,dy,dz,e2,e3
      integer imsg
      character errmsg(3)*80
      save errmsg
      data errmsg(1)/'carlsnf:  smallest input argument < 0.'/
      data errmsg(2)/'carlsnf:  input arguments underflow'/
      data errmsg(3)/'carlsnf:  input arguments overflow'/
c
c     arrange the array s.t.  xt < yt < zt
      xt = MIN (x,y)
      yt = MAX (x,y)
      zt = MAX (yt,z)
      yt = MIN (yt,z)
      xtmp=xt
      xt = MIN (xtmp,yt)
      yt = MAX (xtmp,yt)
c
      if (xt .lt. 0.) then
         imsg=1
         go to 9999
      end if
      if ((xt+yt) .lt. tiny) then
         imsg=2
         go to 9999
      end if
      if (zt .gt. big) then
         imsg=3
         go to 9999
      end if
c
 101  continue
      lamda = SQRT (xt*yt)+SQRT (yt*zt)+SQRT (zt*xt)
      xt=0.25*(xt+lamda)
      yt=0.25*(yt+lamda)
      zt=0.25*(zt+lamda)
      mu=third*(xt+yt+zt)
      dx=(mu-xt)/mu
      dy=(mu-yt)/mu
      dz=(mu-zt)/mu
      if (MAX (abs(dx),abs(dy),abs(dz)) .gt. errtol)  go to 101
      e2=dx*dy-dz**2
      e3=dx*dy*dz
      carlsnf=(1.+(c1*e2-c2-c3*e3)*e2+c4*e3)/SQRT (mu)
      ifail=0
      return
c
 9999 continue
      write (stderr,'(a)')errmsg(imsg)
      write (stderr,'(a)')'x, y, z:'
      write (stderr,*)x,y,z
      if (ifail .eq. 0)
     .  call STOP ('function CARLSNF: unfinished coding', 205)
      ifail=imsg
      return
c
      end

      function carlsnj(x,y,z,p,ifail)
c
      implicit none
c
      real carlsnj,x,y,z,p
      integer ifail
c     returns the carlson elliptical integeral: Rj(x,y,z,p)
c     x,y,z >=0 and at most one is zero and p .not.=0
c     uses functions 'carlsnc' 'carlsnf'
      real carlsnc,carlsnf
c ------------------------------------------------------------------------------
c
      include 'rrange.i'    ! define trsmall, trbig
c
      real third
      parameter (third=1.0/3.0)
      real tiny,big
      parameter (tiny=2.*trsmall, big=0.2*trbig)
      real errtol
      parameter (errtol=0.0015)
      real c1,c2,c3,c4,c5,c6
      parameter (c1=3./7. ,c2=1./3., c3=3./22., c4=3./11.,
     .     c5=3./13., c6=3./13.)
      integer stderr
      parameter (stderr=6)
      real xt,yt,zt,pt,a,b,rho,tau,rcx
      real sum,fac
      real alpha,beta,mu,lamda
      real srx,sry,srz,dx,dy,dz,dp
      real s2,s3,s4,s5
      real xtmp
      integer imsg,ig
      character errmsg(3)*80
      save errmsg
      data errmsg(1)/'carlsnj:  x, y, or z < 0 .or. p = 0'/
      data errmsg(2)/'carlsnj:  input arguments underflow'/
      data errmsg(3)/'carlsnj:  input arguments overflow'/
c
c     arrange the array s.t. xt < yt < zt
      xt = MIN (x,y)
      yt = MAX (x,y)
      zt = MAX (yt,z)
      yt = MIN (yt,z)
      xtmp=xt
      xt = MIN (xtmp,yt)
      yt = MAX (xtmp,yt)
c
      if (xt .lt. 0. .or. p .eq. 0.) then
         imsg=1
         go to 9999
      end if
      if (MIN (xt+yt,abs(p)) .lt. tiny) then
         imsg=2
         go to 9999
      end if
      if (MAX (zt,abs(p)) .gt. big) then
         imsg=3
         go to 9999
      end if
c
      ig=0
      if (p .gt. 0.) then
         pt=p
      else
         a=1./(yt-p)
         b=a*(zt-yt)*(yt-xt)
         pt=yt+b
         rho=xt*zt/yt
         tau=p*pt/yt
         rcx=carlsnc(rho,tau,ig)
      end if
c
      sum=0.
      fac=1.
 101  continue
      srx = SQRT (xt)
      sry = SQRT (yt)
      srz = SQRT (zt)
      lamda=srx*(sry+srz)+sry*srz
      alpha=(pt*(srx+sry+srz)+srx*sry*srz)**2
      beta=pt*(pt+lamda)**2
      sum=sum+fac*carlsnc(alpha,beta,ig)
      fac=fac*0.25
      xt=0.25*(xt+lamda)
      yt=0.25*(yt+lamda)
      zt=0.25*(zt+lamda)
      pt=0.25*(pt+lamda)
      mu=0.2*(xt+yt+zt+2.*pt)
      dx=(mu-xt)/mu
      dy=(mu-yt)/mu
      dz=(mu-zt)/mu
      dp=(mu-pt)/mu
      if (MAX (abs(dx),abs(dy),abs(dz),abs(dp)) .gt. errtol)  go to 101
      s2=(dx**2+dy**2+dz**2+2.*dp**2)/4.
      s3=(dx**3+dy**3+dz**3+2.*dp**3)/6.
      s4=(dx**4+dy**4+dz**4+2.*dp**4)/8.
      s5=(dx**5+dy**5+dz**5+2.*dp**5)/10.
      carlsnj=3.*sum+fac*(1.+s2*(c1+c3*s2+c5*s3)+s3*c2+s4*c4+s5*c6)/
     .     (mu*SQRT (mu))
      if (p .lt. 0.)carlsnj=a*(b*carlsnj+3.*(rcx-carlsnf(xt,yt,zt,ig)))
      ifail=0
      return
c
 9999 continue
      write (stderr,'(a)')errmsg(imsg)
      write (stderr,'(a)')'x, y, z, p:'
      write (stderr,*)x,y,z,p
      if (ifail .eq. 0)
     .  call STOP ('function CARLSNJ: unfinished coding', 206)
      ifail=imsg
      return
c
      end

      subroutine dvcpr (n,fcni,fcnj,fcnb,xa,xb,ngmax,ngrid,ip,ir,tol,x,
     .                   y,iy,abt,par,work,iwork,ier)
c
c     imsl routines for two point boundray value problems using
c     deferred correction techniques of V. Pereyra
c ------------------------------------------------------------------------------
c                                  specifications for arguments
c
      integer            n,ngmax,ngrid,ip,ir,iy,ier,iwork(1)
      real               xa,xb,tol,x(1),y(iy,1),abt(1),par(5),work(1)
c                                  specifications for local variables
      integer            im17,im2,im3,im4,i,jerror,j,m10,m11,m12,m14,
     .                   m15,m16,m17,m2,m3,m4,m5,m6,m7,m8,m9,mmax2,mmax,
     .                   mtnmax,mt
cSm990901
      character*6        name

      external           fcni,fcnj,fcnb
c                                  first executable statement
      mmax = n
      mtnmax = n*ngmax
      mmax2 = 2*mmax
      m2 = 1+mmax
      m3 = m2+mmax**2
      m4 = m3+mmax**2
      m5 = m4+ngmax
      mt = mtnmax*mmax
      m6 = m5+mt
      m7 = m6+mt
      m8 = m7+mt
      m9 = m8+mtnmax
      m10 = m9+mtnmax
      m11 = m10+mtnmax
      m12 = m11+ngmax
      m14 = m12+mtnmax
      m15 = m14+mtnmax
      m16 = m15+mmax2*mmax
      m17 = m16+mtnmax
      im2 = 1+mtnmax
      im3 = im2+mtnmax
      im4 = im3+ngmax
      if (par(1).eq.0.0) par(4) = 0.0
      if (par(4).eq.0.0)  go to 15
      im17 = m17
      do 10 i=1,ngrid
         do 5 j=1,n
            work(im17) = y(j,i)
            im17 = im17+1
    5    continue
   10 continue
   15 continue
      call dvcps(mmax,n,ngmax,ngrid,ip,ir,mtnmax,mmax2,xa,xb,par,tol,x,
     .work(m17),abt,work(1),work(m2),work(m3),work(m4),work(m5),work(m6)
     .,work(m7),work(m8),work(m9),work(m10),work(m11),work(m12),
     .work(m14),work(m15),work(m16),iwork(1),iwork(im2),iwork(im3),
     .iwork(im4),fcni,fcnj,fcnb,jerror)
      im17 = m17
      do 25 i=1,ngrid
         do 20 j=1,n
            y(j,i) = work(im17)
            im17 = im17+1
   20    continue
   25 continue
      ier = 0
      if (jerror.gt.0) ier = jerror+128
 
cSm990901
      name='dvcpr ' 
      if (ier.gt.0) call uertst(ier,name)
c     if (ier.gt.0) call uertst(ier,6hdvcpr )
      return
c
      end

      subroutine dvcps (mmax,m,nmax,n,p,r,mtnmax,mmax2,a,b,par,tol,x,y,
     .                   abt,alpha,a1,b1,ej,a2,c2,del,uu,res,f,hx,sk,
     .                   gradf,aux,xau,ic,ir,iqj,ica,fcni,fcnj,fcnb,
     .                   jerror)
c
c                                  specifications for arguments
c
      integer            mmax,m,nmax,n,p,r,mtnmax,mmax2,jerror,ic(nmax,
     .                   mmax),ir(nmax,mmax),iqj(nmax),ica(mmax)
      real               a,b,par(5),tol,x(nmax),y(mtnmax),abt(mmax),
     .                   alpha(mmax),a1(mmax,mmax),b1(mmax,mmax),
     .                   ej(nmax),a2(mtnmax,mmax),c2(mtnmax,mmax),
     .                   del(mmax,mtnmax),uu(mtnmax),res(mtnmax),
     .                   f(mtnmax),hx(nmax),sk(mtnmax),gradf(mtnmax),
     .                   aux(mmax2,mmax),xau(mtnmax)
c
c                                  specifications for common /newt /
      common /newt/      inwt,nu,casi
      integer            inwt,nu
      logical            casi
c                                  specifications for common /c1 /
      common /c1/        epsnu,cont
      real               epsnu
      logical            cont
c                                  specifications for local variables
c
      integer            i1,i7,i8,icon,ierror,ifin,ifi,iim,iip,ii,ipm,
     .                   iprint,ip,iq,itemp,i,ji,j,k11,kii1,kii,kij,ki,
     .                   kk,kmax,k,l,mpnm,mpn,n1,n2,nma,nold,npu,np,
     .                   ntem,ntop,nu2,nin,nout
      real               aa(50),alg,auxi,bb(50),bma,c3,c(50),deleps,
     .                   epbar,eps1,epsmac,eps,errnew,errold,e,hcua,hi,
     .                   h,rabs,reps,sig02,tem,te,u1,uun,vain,xn,xte,
     .                   xxn,yki,z1,z2,am16,u2
      logical            lin,sing
      external           fcni,fcnj,fcnb
      data               reps/0.710543e-14/
c                                  first executable statement
      call ugetio(1,nin,nout)
      casi = .false.
      sing = .false.
      if (m.gt.1 .and. m.le.mmax .and. n.gt.3 .and. n.le.nmax .and.
     .p.ge.0 .and. p.lt.m .and. r.ge.0 .and. p+r.le.m)  go to 5
      jerror = 1
      return
c                                  initialization
    5 mpn = m*n
c
      epsmac = 10.*reps
c
      epsmac = MAX (epsmac,1.e-4*tol)
      if (par(1).gt.0.)  go to 10
      deleps = 0.
      iprint = 0
      vain = 0.
      lin = .false.
      go to 15
   10 deleps = par(2)
      iprint = par(3)+0.5
      vain = par(4)
      lin = par(5).ne.0.
c
   15 n1 = n-1
      epsnu = 0.
      if (deleps.eq.0.) epsnu = 1.
      rabs = 1.
      jerror = 0
      epbar = 1.e10
      cont = .false.
      icon = 3
      nold = 1
      nma = nmax
      am16 = nmax-16
      ntop = MAX (.9*nmax,am16)
      bma = 0.07
      bb(1) = 1.
      do 20 i=2,50
   20 bb(i) = 0.
      do 25 nu=1,20
         nu2 = 2*nu+1
         aa(nu2) = -nu/(2.**(2*nu-1)*nu2)
         aa(nu2+1) = 0.
   25 continue
      if (vain.gt.0.)  go to 40
c                                  first approximation for y and x
      x(1) = a
      x(n) = b
      h = (b-a)/n1
      do 30 i=2,n1
   30 x(i) = a+(i-1)*h
      do 35 i=1,mpn
   35 y(i) = 0.
c
   40 h = 0.
      do 45 i=1,n1
         hx(i) = x(i+1)-x(i)
         if (hx(i).gt.h) h = hx(i)
   45 continue
      hcua = h**2
      if (iprint.ne.0) write (nout,50) tol, h
   50 format (11h tolerance=, e12.2, 8h max. h=, e12.2)
c                                  main body of pasva3
      nu = 0
      eps = MAX (epsmac,h**2*.1)
c                                  enter after mesh change
   55 errold = 1.0e20
      kmax = max0(15,(n-2)/2)
      mpn = m*n
      n1 = n-1
      mpnm = m*n1
      ipm = mpnm+p+1
      c3 = 0.8
      do 60 i=1,mpn
   60 sk(i) = 0.
      if (nu.eq.0)  go to 80
c                                  after mesh change we have to
c                                    initialize sk if nu .gt. 0
      do 65 i7=1,n
         i8 = (i7-1)*m+1
         call fcni(m,x(i7),y(i8),f(i8))
   65 continue
      call fcnb(m,y(1),y(i8),alpha)
      call dvcpt(nu,2,2,n1,m,aa,x,f,res,ierror)
      do 75 i=1,n1
         ki = (i-1)*m
         do 70 j=1,m
            kij = ki+j
            sk(kij) = hx(i)*res(kij)
   70    continue
   75 continue
      if (nu.lt.kmax)  go to 80
      nu = nu-1
      go to 230
c                                  nu is too large go to refine the
c                                    mesh. newton iteration
   80 if (epsnu.ge.1.)  go to 90
      eps1 = eps
   85 eps = MAX (rabs,eps)
   90 call dvcpv(m,n,p,r,alpha,a1,b1,x,y,lin,a2,c2,del,jerror,iprint,
     .eps,ir,ic,uu,res,mmax,mtnmax,nmax,mmax2,f,hx,sk,gradf,aux,ica,xau,
     .fcni,fcnj,fcnb)
      if (jerror.ne.3)  go to 95
      return
   95 if (epsnu.ge.1.)  go to 150
c                                  continuation-- choose step and new
c                                    initial profile
      cont = .true.
      do 100 i7=1,n
         i8 = (i7-1)*m+1
         call fcni(m,x(i7),y(i8),f(i8))
  100 continue
      call fcnb(m,y(1),y(i8),alpha)
      cont = .false.
      do 105 i=1,p
  105 uu(i) = -alpha(i)
      do 115 i=1,n1
         ii = (i-1)*m
         iip = ii+p
         iim = ii+m
         do 110 j=1,m
  110    uu(iip+j) = .5*hx(i)*(f(iim+j)+f(ii+j))
  115 continue
      ip = p+1
      do 120 i=ip,m
  120 uu(n1*m+i) = -alpha(i)
      call dvcpy(a2,c2,del,uu,m,n,p,r,ir,ic,uu,mtnmax,mmax,nmax,xau)
      if (inwt-1) 125, 125, 130
  125 deleps = 2.*deleps
      go to 135
  130 deleps = MAX (.01,deleps/(inwt-1))
  135 do 140 i=1,mpn
  140 y(i) = y(i)+deleps*uu(i)
      epsnu = MIN (epsnu+deleps,1.)
      if (iprint.ne.0) write (nout,145) deleps, epsnu, eps
  145 format (11h del-epsnu=, e15.7, 7h epsnu=, e15.7, 5h eps=, e15.7)
      if (epsnu.lt.1.)  go to 85
      eps = eps1
      go to 90
c                                  correction and error control starts
  150 call dvcpt(nu+1,2,2,n1,m,aa,x,f,res,ierror)
      if (ierror.eq.1)  go to 230
      do 160 i=1,n1
         ii = (i-1)*m
         do 155 j=1,m
            ki = ii+j
            auxi = res(ki)*hx(i)
            res(ki) = sk(ki)-auxi
            sk(ki) = auxi
  155    continue
  160 continue
      if (icon.le.12)  go to 275
  165 if (p.eq.0)  go to 175
      do 170 i=1,p
  170 uu(i) = 0.
  175 do 180 i=ipm,mpn
  180 uu(i) = 0.
      do 185 i=1,mpnm
  185 uu(i+p) = res(i)
      call dvcpw(m,n,p,r,x,y,a1,b1,a2,c2,del,.true.,sing,ir,ic,uu,uu,
     .lin,mmax,mtnmax,nmax,mmax2,hx,gradf,aux,ica,xau,fcnj,fcnb)
c
c                                  estimate for max. absolute error (by
c                                    components)
      icon = 15
      errnew = 0.
      do 190 j=1,m
  190 abt(j) = 0.
      do 200 i=1,n
         kk = (i-1)*m
         do 195 j=1,m
            u1 = abs(uu(kk+j))
            u2 = u1/MAX (1.0,abs(y(kk+j)))
            if (u1.gt.abt(j)) abt(j) = u1
            if (u2.gt.errnew) errnew = u2
  195    continue
  200 continue
      k = nu+1
      if (iprint.eq.0)  go to 215
      write (nout,205) (abt(j),j=1,m)
      write (nout,210) errnew, nu
  205 format (30h estimated error by components/(10e12.3))
  210 format (17h estimated error=, e12.3, 14h in correction, i4)
  215 if (errnew.gt.tol)  go to 220
      jerror = 0
      return
c                                  precision achieved if not enough
c                                    points , refine the mesh
  220 if (nu+1.ge.kmax)  go to 275
      if (errnew.le.0.1*errold)  go to 225
      if (errnew.gt.c3*errold)  go to 230
      c3 = 0.5*c3
c                                  either keep correcting .
  225 errold = errnew
      eps = MAX (epsmac,1.e-3*errold)
      nu = nu+1
      go to 90
c                                  or refine the mesh, unless jerror = 4
  230 if (jerror.ne.4)  go to 235
      return
  235 if (n.le.ntop)  go to 240
c                                  too many grid points
      return
  240 eps = MAX (epsmac,1.e-3*errold)
      if (errold.le.1.)  go to 245
      epbar = MIN (1.0,1.e-2*errold)
      go to 250
  245 epbar = .01*errold
  250 icon = 0
      nold = n
      bma = 1.
      if (nu.lt.1)  go to 275
      do 260 i=1,n1
         ii = (i-1)*m
         do 255 j=1,m
            ki = ii+j
            sk(ki) = res(ki)+sk(ki)
  255    continue
  260 continue
      nu = nu-1
      if (nu.eq.0)  go to 275
      call dvcpt(nu,2,2,n1,m,aa,x,f,res,ierror)
      do 270 i=1,n1
         ii = (i-1)*m
         do 265 j=1,m
            ki = ii+j
            res(ki) = res(ki)*hx(i)-sk(ki)
  265    continue
  270 continue
c                                  mesh variation equidistribution of
c                                    the l2 norm of the error for the
c                                    o(h**(2*nu+2)) method.
  275 icon = icon+1
      if (iprint.ne.0) write (nout,280)
  280 format (18h ---step change---)
      alg = 1.5
      sig02 = 1./(2*nu+2)
      tem = 0.
      uun = 0.
      do 290 i=1,n1
         te = 0.
         ki = (i-1)*m
         do 285 j=1,m
            k11 = ki+j
            z1 = abs(y(k11))
            z2 = abs(res(k11))
            if (z1.gt.tem) tem = z1
            if (z2.gt.te) te = z2
  285    continue
         ej(i) = (te/hx(i))**sig02
         uun = uun+ej(i)
  290 continue
c
      if (icon.gt.1 .and. nold.gt.1)  go to 295
      epbar = MAX (MIN (epbar,tem*bma),tol)
      e = epbar**sig02
  295 iq = 0
      n2 = n-2
      i = 0
  300 i = i+1
      iqj(i) = ej(i)/e-0.33
      iq = iq+iqj(i)
      if (i.le.n2)  go to 300
      ifin = .04*n
      nma = min0(nmax-n,70)
      if (iprint.ne.0) write (nout,305) iq, nma
  305 format (12h new points=, i4, 6h nmax=, i4)
      if (iq.le.ifin .or. nma.le.0)  go to 385
      if (iq.le.nma)  go to 315
c                                  we attempt to diminish the number of
c                                    points to be introduced
      if (alg.lt.1.09)  go to 385
      alg = alg-.1
      xn = n+iq
      xte = n+nma
      e = e*xn/MIN (alg*n,xte)
      if (iprint.eq.0)  go to 295
      write (nout,310) e, alg
  310 format (7h level=, e15.5, 5h alg=, e15.5)
      go to 295
c                                  construct new mesh
  315 j = 2
      do 320 i=1,mpn
  320 uu(i) = y(i)
      sk(1) = a
      npu = 2*nu+3
      do 370 i=1,n1
         kii1 = i*m
         ifi = j+iqj(i)
         hi = hx(i)/(iqj(i)+1)
         do 365 l=j,ifi
            sk(l) = x(i)+hi*(l-j+1)
            kii = (l-1)*m
            if (iqj(i).gt.0)  go to 330
            do 325 k=1,m
  325       y(kii+k) = uu(kii1+k)
            go to 365
  330       if (i.le.nu+2)  go to 340
            if (i.gt.n-(nu+3))  go to 345
            np = nu+2
  335       itemp = i
            call dvcpu(itemp,npu,np,c,bb,x,sk(l))
            go to 350
  340       np = i
            go to 335
  345       np = i-n+npu
            go to 335
  350       do 360 k=1,m
               yki = 0.
               do 355 ji=1,npu
  355          yki = yki+c(ji)*uu((i-np+ji-1)*m+k)
               ki = kii+k
               y(ki) = yki
  360       continue
  365    continue
         j = ifi+1
  370 continue
c
      casi = .false.
      n = n+iq
      n1 = n-1
      ki = m*n1
      do 375 i=1,m
  375 y(ki+i) = uu(mpnm+i)
      h = 0.
      do 380 i=2,n1
         i1 = i-1
         hx(i1) = sk(i)-sk(i1)
         if (hx(i1).gt.h) h = hx(i1)
         ki = i1*m+1
         x(i) = sk(i)
  380 continue
c
      x(n) = b
      hx(n1) = x(n)-x(n1)
      if (hx(n1).gt.h) h = hx(n1)
      hcua = h**2
      if (abs(epbar-tem*bma).lt.1.e-5) epbar = 1.e10
      if (icon.ge.5) icon = 15
      go to 55
c                                  end of mesh selection
c
  385 if (nold.eq.1 .and. alg.lt.1.5 .and. 2*n-1.le.nmax)  go to 405
      if (n.gt.nold)  go to 165
      if (alg.lt.1.45)  go to 390
      alg = 1.4
      xn = n+iq
      xxn = nmax
      e = e*xn/MIN (alg*n,xxn)
      if (iprint.eq.0)  go to 295
      write (nout,310) e, alg
      go to 295
  390 if (nu.gt.0)  go to 415
c                                  bisection
      if (iprint.eq.0)  go to 400
      write (nout,395)
  395 format (23h bisection is performed)
  400 ntem = 2*n-1
      if (ntem.le.nmax)  go to 405
c                                  too many grid points
      return
  405 do 410 i=1,n1
  410 iqj(i) = 1
      iq = n1
      icon = 0
      go to 315
  415 nu = nu-1
      icon = 0
c                                  since we are going back we must
c                                    recompute the estimate of the
c                                    error.
      mpnm = m*n1
      do 420 i=1,mpnm
  420 sk(i) = res(i)+sk(i)
      if (nu.gt.0)  go to 435
      do 430 i=1,n1
         ki = (i-1)*m
         do 425 j=1,m
            kij = ki+j
            res(kij) = -sk(kij)
  425    continue
  430 continue
      go to 275
  435 call dvcpt(nu,2,2,n1,m,aa,x,f,res,ierror)
      do 445 i=1,n1
         hi = hx(i)
         ki = (i-1)*m
         do 440 j=1,m
            kij = ki+j
            res(kij) = hi*res(kij)-sk(kij)
  440    continue
  445 continue
      go to 275
c                                  end of mesh variation
c
      end

      subroutine dvcpt (k,p,q,n,m,a,x,y,s,ierror)
c
c                                  specifications for arguments
      integer            k,p,q,n,m,ierror
      real               a(1),x(1),y(1),s(1)
c                                  specifications for local variables
c
      integer            i1,iad,ii1,ii,im,ite,it,i,j,kk1,kk,kmid1,
     .                   kmidp1,kmid,l,nf
      real               acum,c(50)
c                                  error exit
c                                  first executable statement
      if (k.gt.(n+1-q)/p .or. p.lt.1 .or. k.lt.1)  go to 65
      if (q.eq.0)  go to 10
      do 5 i=1,q
    5 a(i) = 0.
   10 kk1 = q+p*k
      kk = kk1-1
      kmid = kk1/2
      ierror = 0
      kmid1 = kmid-1
c                                  unsymmetric approximation left
c                                    boundary
      if (kmid1.lt.1)  go to 30
      do 25 i=1,kmid1
         ite = i
         call dvcpu(ite,kk1,ite,c,a,x,.5*(x(i+1)+x(i)))
         im = (i-1)*m
         do 20 l=1,m
            acum = 0.
            do 15 j=1,kk1
   15       acum = acum+c(j)*y((j-1)*m+l)
            s(im+l) = acum
   20    continue
   25 continue
c                                  center range
   30 nf = n+1-kk1+kmid
      do 45 i=kmid,nf
         ite = i
         call dvcpu(ite,kk1,kmid,c,a,x,.5*(x(i+1)+x(i)))
         i1 = i-1
         ii = i1-kmid
         it = i1*m
         do 40 l=1,m
            acum = 0.
            do 35 j=1,kk1
   35       acum = acum+c(j)*y((ii+j)*m+l)
            s(it+l) = acum
   40    continue
   45 continue
c                                  right boundary
      kmidp1 = kmid+1
      ii = n-kk
      ii1 = ii-1
      do 60 i=kmidp1,kk
         iad = ii+i
         ite = i
         call dvcpu(iad,kk1,ite,c,a,x,.5*(x(iad+1)+x(iad)))
         it = (ii1+i)*m
         do 55 l=1,m
            acum = 0.
            do 50 j=1,kk1
   50       acum = acum+c(j)*y((ii1+j)*m+l)
            s(it+l) = acum
   55    continue
   60 continue
c                                  regular exit
      return
   65 ierror = 1
      return
c
      end

      subroutine dvcpu (i0,n,np,c,bb,x,xbar)
c
c                                  specifications for arguments
      integer            i0,n,np
      real               c(1),bb(1),x(1),xbar
c                                  specifications for local variables
      integer            i,jm1,j,km1,k,ll,n1,nn
      real               alf(50),hi
c                                  first executable statement
      do 5 i=1,n
    5 c(i) = bb(i)
      do 10 i=1,n
         hi = x(i0+1)-x(i0)
         alf(i) = (x(i0-np+i)-xbar)/hi
   10 continue
      nn = n-1
      n1 = n+1
      do 20 i=1,nn
         ll = n-i
         do 15 j=1,ll
            k = n1-j
            c(k) = c(k)-alf(i)*c(k-1)
   15    continue
   20 continue
      do 30 i=1,nn
         k = n-i
         km1 = k+1
         do 25 j=km1,n
            c(j) = c(j)/(alf(j)-alf(j-k))
            jm1 = j-1
            c(jm1) = c(jm1)-c(j)
   25    continue
   30 continue
c
      return
c
      end

      subroutine dvcpv (m,n,p,r,alpha,a1,b1,x,y,lin,a2,c2,del,jerror,
     .                   iprint,eps,ir,ic,uu,res,mmax,mtnmax,nmax,mmax2,
     .                   f,hx,sk,gradf,aux,ica,xau,fcni,fcnj,fcnb)
c
c                                  specifications for arguments
c
      integer            m,n,p,r,jerror,iprint,mmax,mtnmax,nmax,mmax2,
     .                   ir(nmax,mmax),ic(nmax,mmax),ica(mmax),nin,nout
      real               alpha(mmax),a1(mmax,mmax),b1(mmax,mmax),x(nmax)
     .                   ,y(mtnmax),a2(mtnmax,mmax),c2(mtnmax,mmax),
     .                   del(mmax,mtnmax),eps,uu(mtnmax),res(mtnmax),
     .                   f(mtnmax),hx(nmax),sk(mtnmax),gradf(mtnmax),
     .                   aux(mmax2,mmax),xau(mtnmax)
      logical            lin
c                                  specifications for common /c1 /
      common /c1/        epsnu,cont
      real               epsnu
      logical            cont
c                                  specifications for common /newt /
      common /newt/      inwt,nu,casi
      integer            inwt,nu
      logical            casi
c                                  specifications for local variables
c
      integer            i1,i7,i8,i,j,k1jm,k1j,k1,mpnm,mpn,mp,p1
      real               dt,gnor,pnor,rabs1,rabs,reold,rol1,rol,scpr,
     .                   te,tmin,tn,t,tpm8
      logical            sing,step
      external           fcnj,fcnb,fcni
c                                  error exit jerror = 3 non
c                                    convergence. newton iteration
c                                    starts
c                                  first executable statement
      call ugetio(1,nin,nout)
      dt = 1.
      mp = m-p
      p1 = p+1
      t = 1.
      te = 1.
      tpm8 = 2.**(-8)
      tmin = tpm8
      step = .false.
      scpr = 0.
      mpn = m*n
      do 5 i=1,mpn
    5 uu(i) = 0.
      mpnm = m*(n-1)
      inwt = 0
      rol = 1.e20
      reold = 1.0e20
      if (iprint.eq.0)  go to 20
      write (nout,10) n
   10 format (7h ngrid=, i5)
      write (nout,15) nu, eps
   15 format (12h correction=, i4, 33h residual for newton should be .l,
     .2he., e14.7)
c                                  residual computation
   20 rabs = 0.
      do 25 i7=1,n
         i8 = (i7-1)*m+1
         call fcni(m,x(i7),y(i8),f(i8))
   25 continue
      call fcnb(m,y(1),y(i8),alpha)
      if (p.eq.0)  go to 35
      do 30 i=1,p
         res(i) = -alpha(i)
         rabs = rabs+res(i)**2
   30 continue
c
   35 do 45 i=2,n
         i1 = i-1
         k1 = i1*m
         do 40 j=1,m
            k1j = k1+j
            k1jm = k1j-m
            res(k1j-mp) = -y(k1j)+y(k1jm)+.5*hx(i1)*(f(k1j)+f(k1jm))
     .      +sk(k1jm)
            rabs = rabs+res(k1j-mp)**2
   40    continue
   45 continue
c
      do 50 j=p1,m
         res(mpnm+j) = -alpha(j)
         rabs = rabs+alpha(j)**2
   50 continue
      rabs1 = SQRT (rabs)
      if (iprint.eq.0)  go to 65
      write (nout,55) inwt, rabs1
   55 format (18h newton iteration=, i3, 20h max. abs. residual=, e14.7)
      write (nout,60) (res(j),j=1,mpn)
   60 format (11h residuals=/(1x, 6e12.3))
   65 if (inwt.eq.0 .and. (epsnu.eq.1. .or. epsnu.eq.0.))  go to 140
c
c                                  check for convergence
      if (rabs1.gt.eps)  go to 70
      if (epsnu.lt.1.)  return
c                                  jacobian is kept constant until next
c                                    mesh change. step and angle
c                                    control are shortcircuited.
      casi = .true.
      return
c                                  newton exit.
c                                  convergence or too many iterations.
c                                  newton test in order to avoid cycling
c
   70 if ((reold-rabs.ge..5*t*scpr .and. inwt.lt.15) .or. (inwt.eq.1))
     .go to 140
      if (inwt.eq.15)  go to 115
c                                  step control starts
      if (.not.(casi .or. lin))  go to 80
      if (lin)  go to 75
      casi = .false.
      go to 145
   75 if (rabs1.le.100.*eps)  go to 120
      go to 135
   80 if (step)  go to 95
      step = .true.
      if (te.le.1.)  go to 85
      tn = te
      go to 90
   85 tn = 100.*te
      if (tn.gt.1.)  go to 95
   90 tmin = tn*tpm8
      go to 100
   95 tn = .5*t
      if (tn.lt.tmin)  go to 115
      rol = rabs
  100 dt = tn-t
      t = tn
      do 105 i=1,mpn
  105 y(i) = y(i)+dt*uu(i)
      if (iprint.eq.0)  go to 20
      write (nout,110) t, tmin
  110 format (20h step control=t - dt, 2e15.7)
      go to 20
c
  115 rol1 = SQRT (rol)
      if (rol1.ge.100.*eps)  go to 135
      rabs1 = rol1
  120 do 125 i=1,mpn
  125 y(i) = y(i)-dt*uu(i)
      if (iprint.ne.0) write (nout,130) nu, rabs1, eps
  130 format (20h warning--correction, i4/26h residual in newton could ,
     .18honly be reduced to, e12.2, 17h instead of below,
     .e12.2/56h if upon termination tol has not been reached the proble,
     .1hm, 32h requires a longer computer word)
c                                  newton did not quite reach the
c                                    tolerance. further mesh
c                                    refinements are not allowed.
      jerror = 4
      if (epsnu.lt.1.)  return
c                                  jacobian is kept constant until next
c                                    mesh change. step and angle
c                                    control are shortcircuited.
      casi = .true.
      return
c                                  we assume divergence and return
c                                    error exit 3
  135 jerror = 3
      return
c
  140 if (inwt.eq.4) casi = .false.
  145 call dvcpw(m,n,p,r,x,y,a1,b1,a2,c2,del,casi,sing,ir,ic,uu,res,lin,
     .mmax,mtnmax,nmax,mmax2,hx,gradf,aux,ica,xau,fcnj,fcnb)
      if (sing .and. iprint.ne.0) write (nout,150)
  150 format (21h jacobian is singular)
      reold = rabs
      scpr = 1.e-10
      if (casi .or. lin)  go to 175
      tmin = tpm8
      gnor = 0.
      t = 1.
      step = .false.
      scpr = 0.
      pnor = 0.
      do 155 i=1,mpn
         pnor = pnor+uu(i)**2
         gnor = gnor+gradf(i)**2
         scpr = scpr+gradf(i)*uu(i)
  155 continue
      if (pnor.eq.0.)  go to 115
c                                  we check if the direction uu is of
c                                    descent (as it should) and also if
c                                    the identity (gradf,uu)=rabs is
c                                    approximately verified. if either
c                                    one of these checks fail, we take
c                                    gradf as our next search
c                                    direction, since the above
c                                    indicates that uu is an unreliable
c                                    direction due to ill-conditioning
c                                    in the jacobian.
      if (scpr) 165, 165, 160
  160 if (abs(scpr-rabs).le..1*rabs)  go to 170
  165 scpr = -1.
  170 te = 1.
c                                  approximate solution is corrected
  175 do 180 i=1,mpn
         if (scpr.le.0.) uu(i) = gradf(i)
         y(i) = y(i)+uu(i)
  180 continue
      if (scpr.le.0.) scpr = gnor
      t = 1.
      inwt = inwt+1
      go to 20
c
      end

      subroutine dvcpw (m,n,p,r,x,y,a1,b1,a,c,del,casi,sing,ir,ic,uu,
     .    res,lin,mmax,mtnmax,nmax,mmax2,hx,gradf,aux,ica,xau,fcnj,fcnb)
c
c                                  specifications for arguments
c
      integer            m,n,p,r,mmax,mtnmax,nmax,mmax2,ir(nmax,mmax),
     .                   ic(nmax,mmax),ica(1)
      real               x(nmax),y(mtnmax),a1(mmax,mmax),b1(mmax,mmax),
     .                   a(mtnmax,mmax),c(mtnmax,mmax),del(mmax,mtnmax),
     .                   uu(mtnmax),res(mtnmax),hx(nmax),gradf(mtnmax),
     .                   aux(mmax2,mmax),xau(mtnmax)
      logical            casi,sing,lin
c                                  specifications for common /newt /
      common             /newt/ inwt,nu,cas1
      integer            inwt,nu
      logical            cas1
c                                  specifications for local variables
c
      integer            i1j,i1,i2j,i2,i3p,i3,i4,i7,i8,i9,ii,im,ip,it,i,
     .                   j7,jj,j,k,m1,mmp,n1,p1,pm,r1
      real               econd,h1,hak,hbk,reps,te
      external           fcnj, fcnb
      data               reps/0.710543e-14/
c                                  casi = .true. no decomposition
c                                  first executable statement
      if (casi)  go to 165
      p1 = p-1
      n1 = n-1
      m1 = m-p
      mmp = m1-1
      pm = p+1
      r1 = r-1
c                                  construction of a,c,del
      do 15 i=1,n
         jj = (i-1)*m
         i8 = (i-1)*m+1
         call fcnj(m,x(i),y(i8),del)
         do 10 k=1,m
            do 5 j=1,m
               c(jj+j,k) = del(j,k)
    5       continue
   10    continue
   15 continue
      econd = SQRT (reps)
      do 30 i7=1,m
         hak = MAX (y(i7),0.1)*econd
         y(i7) = y(i7)+hak
         call fcnb(m,y(1),y(i8),del(1,2))
         y(i7) = y(i7)-2.0*hak
         call fcnb(m,y(1),y(i8),del(1,1))
         y(i7) = y(i7)+hak
         do 20 j7=1,m
            a1(j7,i7) = (del(j7,2)-del(j7,1))/(2.0*hak)
   20    continue
         i9 = i8+i7-1
         hbk = MAX (y(i9),0.1)*econd
         y(i9) = y(i9)+hbk
         call fcnb(m,y(1),y(i8),del(1,2))
         y(i9) = y(i9)-2.0*hbk
         call fcnb(m,y(1),y(i8),del(1,1))
         y(i9) = y(i9)+hbk
         do 25 j7=1,m
            b1(j7,i7) = (del(j7,2)-del(j7,1))/(2.0*hbk)
   25    continue
   30 continue
      call fcnj(m,x(n),y(i8),del)
      if (p.eq.0)  go to 45
      do 40 i=1,p
         do 35 j=1,m
   35    a(i,j) = a1(i,j)
   40 continue
   45 do 85 ii=1,n1
         h1 = .5*hx(ii)
         i1 = (ii-1)*m+1
         i2 = i1+mmp
         it = 0
         do 55 i=i1,i2
            ip = i+p
            it = it+1
            do 50 j=1,m
               a(ip,j) = -h1*c(i,j)
               if (it.eq.j) a(ip,j) = a(ip,j)-1.
   50       continue
   55    continue
         i3 = i2+m
         i4 = i1+p1
         it = 0
         if (p.eq.0)  go to 70
         do 65 i=i1,i4
            it = it+1
            do 60 j=1,m
               im = i+m
               a(im,j) = -h1*c(i3+it,j)
               c(i,j) = -h1*c(i2+it,j)
               if (it+m1.ne.j)  go to 60
               a(im,j) = a(im,j)+1.
               c(i,j) = c(i,j)-1.
   60       continue
   65    continue
   70    it = 0
         do 80 i=i1,i2
            ip = i+p
            it = it+1
            im = i+m
            do 75 j=1,m
               c(ip,j) = -h1*c(im,j)
               if (it.eq.j) c(ip,j) = c(ip,j)+1.
   75       continue
   80    continue
   85 continue
c
      i1 = n1*m+pm
      i2 = n*m
      it = p
      do 95 i=i1,i2
         it = it+1
         do 90 j=1,m
   90    a(i,j) = b1(it,j)
   95 continue
      if (r.eq.0)  go to 110
      do 105 i=1,r
         do 100 j=1,m
  100    del(i,j) = a1(i+p,j)
  105 continue
  110 if (lin)  go to 160
c                                  computation of gradient
      it = 0
      do 145 i=1,n
         i1 = (i-1)*m
         i2 = i1-m1
         i3 = i*m
         do 140 j=1,m
            it = it+1
            te = 0.
            do 115 jj=1,m
               i1j = i1+jj
               te = te+a(i1j,j)*res(i1j)
  115       continue
            if (i.eq.n)  go to 125
            do 120 jj=1,p
  120       te = te+c(i1+jj,j)*res(i3+jj)
  125       if (i.eq.1)  go to 135
            do 130 jj=1,m1
               i2j = i2+jj
               te = te+c(i2j,j)*res(i2j)
  130       continue
  135       gradf(i1+j) = te
  140    continue
  145 continue
c
      if (r.eq.0)  go to 160
      i3p = n1*m+p
      do 155 j=1,m
         te = 0.
         do 150 i=1,r
  150    te = te+del(i,j)*res(i3p+i)
         gradf(j) = gradf(j)+te
  155 continue
  160 call dvcpx(a,c,del,m,n,p,r,ir,ic,sing,ica,aux,mtnmax,mmax,nmax,
     .mmax2)
      if (sing)  go to 170
  165 call dvcpy(a,c,del,res,m,n,p,r,ir,ic,uu,mtnmax,mmax,nmax,xau)
  170 return
c
      end

      subroutine dvcpx (a,c,del,m,n,p,r,ir,ic,sing,ica,aux,mtnmax,mmax,
     .                  nmax,mmax2)
c
c                                  specifications for arguments
c
      integer            m,n,p,r,mtnmax,mmax,nmax,mmax2,ir(nmax,mmax),
     .                   ic(nmax,mmax),ica(mmax)
      real               a(mtnmax,mmax),c(mtnmax,mmax),del(mmax,mtnmax),
     .                   aux(mmax2,mmax)
      logical            sing
c                                  specifications for local variables
c
      integer            i01,i0,i1s,i1,i2,i3,i4,i7,i8,ii1,ii2,ii,ini,
     .                   ipiv,ipp,ip,is1,ism,isp,is,ite,ixi,ixm,ixp,ix,
     .                   i,j1,jp,jsx,js,jx,j,ks,m1,mp,n1,p1
      real               ak,te,xmu,y
c                                  first executable statement
      sing = .false.
      p1 = p+1
      mp = m+p
      n1 = n-1
c                                  main loop
      do 275 ii=1,n1
         ix = (ii-1)*m
         ixm = ix+m
         ixi = ix-m
         do 10 i=1,m
            ic(ii,i) = i
            ir(ii,i) = i
            do 5 j=1,m
    5       aux(i,j) = a(ix+i,j)
   10    continue
         if (p) 65, 65, 15
   15    do 25 i=1,p
            do 20 j=1,m
   20       aux(i+m,j) = c(ix+i,j)
   25    continue
c                                  reduction of first p rows of a1(ii)
         do 60 i=1,p
            te = 0.
            do 30 j=i,m
               y = abs(aux(i,j))
               if (y.le.te)  go to 30
               te = y
               ipiv = j
   30       continue
            if (te.eq.0)  go to  335
c                                  column interchanges
            if (ipiv.eq.i)  go to 45
            ite = ic(ii,i)
            ic(ii,i) = ic(ii,ipiv)
            ic(ii,ipiv) = ite
            do 35 is=1,mp
               te = aux(is,i)
               aux(is,i) = aux(is,ipiv)
               aux(is,ipiv) = te
   35       continue
            if (ii.eq.1)  go to 45
            i1 = ixi+p1
            i2 = ixi+m
            do 40 is=i1,i2
               te = c(is,i)
               c(is,i) = c(is,ipiv)
               c(is,ipiv) = te
   40       continue
c                                  factorization
   45       ak = 1./aux(i,i)
            i1 = i+1
            do 55 is=i1,mp
               xmu = aux(is,i)*ak
               aux(is,i) = xmu
               do 50 js=i1,m
   50          aux(is,js) = aux(is,js)-xmu*aux(i,js)
   55       continue
   60    continue
c                                  reduction of remaining part of
c                                    a1(ii)
   65    do 110 j=p1,m
            te = 0.
            do 70 i=j,mp
               y = abs(aux(i,j))
               if (y.le.te)  go to 70
               te = y
               ipp = i
   70       continue
            if (te.eq.0.)  go to  335
c                                  row interchanges
            if (ipp.eq.j)  go to 95
            ipiv = ipp-p
            jp = j-p
            ite = ir(ii,jp)
            ir(ii,jp) = ir(ii,ipiv)
            ir(ii,ipiv) = ite
            do 75 js=1,m
               te = aux(j,js)
               aux(j,js) = aux(ipp,js)
               aux(ipp,js) = te
   75       continue
            ii1 = ix+j
            ii2 = ix+ipp
            if (ipp.gt.m)  go to 85
            do 80 js=1,m
               te = c(ii1,js)
               c(ii1,js) = c(ii2,js)
               c(ii2,js) = te
   80       continue
            go to 95
   85       do 90 js=1,m
               te = c(ii1,js)
               c(ii1,js) = a(ii2,js)
               a(ii2,js) = te
   90       continue
c                                  factorization
   95       if (j.eq.m)  go to 110
            ak = 1./aux(j,j)
            j1 = j+1
            do 105 is=j1,mp
               xmu = aux(is,j)*ak
               aux(is,j) = xmu
               do 100 js=j1,m
  100          aux(is,js) = aux(is,js)-xmu*aux(j,js)
  105       continue
  110    continue
c                                  restoring b1(ii+1) and al(ii)
         if (p.eq.0)  go to 125
         i1 = m-p+1
         do 120 is=i1,m
            ip = ir(ii,is)
            if (ip.eq.is)  go to 120
            isp = is+p+ixi
            ipp = ip+p+ix
            do 115 js=1,m
  115       c(isp,js) = a(ipp,js)
  120    continue
  125    do 135 i=1,m
            do 130 j=1,m
  130       a(ix+i,j) = aux(i,j)
  135    continue
c
         is = ix+1
         isp = ix+p
         do 140 js=1,m
  140    ica(js) = ic(ii,js)
         js = 1
  145    if (js.ge.m)  go to 180
         ip = ica(js)
         if (ip.ne.js)  go to 150
         js = js+1
         go to 145
c
  150    ini = js
         ica(ini) = ini
c                                  column interchanges for b1(ii+1)
         if (p.eq.0)  go to 165
  155    do 160 i1=is,isp
            te = c(i1,ini)
            c(i1,ini) = c(i1,ip)
            c(i1,ip) = te
  160    continue
  165    if (r.eq.0)  go to 175
c                                  column interchanges of d(ii)
         ipp = ix+ip
         jx = ix+ini
         do 170 is1=1,r
            te = del(is1,ipp)
            del(is1,ipp) = del(is1,jx)
            del(is1,jx) = te
  170    continue
  175    ini = ip
         ip = ica(ip)
         ica(ini) = ini
         if (ip.ne.js)  go to 155
         js = js+1
         go to 145
c
  180    if (r.eq.0)  go to 230
c                                  solution of de(ii)*al(ii)=-
c                                    delta(ii)
         i1 = ix+1
         do 210 is=1,r
            do 195 js=i1,ixm
               te = del(is,js)
               jsx = js-ix
               if (js.eq.i1)  go to 190
               j1 = js-1
               do 185 ks=i1,j1
  185          te = te-del(is,ks)*a(ks,jsx)
  190          del(is,js) = te/a(js,jsx)
  195       continue
            m1 = m-1
            do 205 js=1,m1
               j = ixm-js
               te = del(is,j)
               j1 = j+1
               do 200 ks=j1,ixm
  200          te = te-del(is,ks)*a(ks,j-ix)
               del(is,j) = te
  205       continue
  210    continue
c                                  de(ii+1)=-de(ii)*ga(ii)
         ixp = ix+p1
         do 225 is=1,r
            do 220 js=1,m
               te = 0.
               do 215 ks=ixp,ixm
  215          te = te-del(is,ks)*c(ks,js)
               del(is,ixm+js) = te
  220       continue
  225    continue
c                                  solution of be(ii+1)*al(ii)=b(ii+1)
  230    i1 = ix+1
         if (p.eq.0)  go to 275
         i2 = ix+p
         do 270 is=i1,i2
            do 245 js=1,m
               te = c(is,js)
               if (js.eq.1)  go to 240
               j1 = js-1
               do 235 ks=1,j1
  235          te = te-c(is,ks)*a(ix+ks,js)
  240          c(is,js) = te/a(ix+js,js)
  245       continue
            m1 = m-1
            do 255 js=1,m1
               j = m-js
               te = c(is,j)
               j1 = j+1
               do 250 ks=j1,m
  250          te = te-c(is,ks)*a(ix+ks,j)
               c(is,j) = te
  255       continue
c                                  al(ii+1)=a(ii+1)-be(ii+1)*ga(ii)
            ism = is+m
            do 265 js=1,m
               te = a(ism,js)
               do 260 ks=p1,m
  260          te = te-c(is,ks)*c(ix+ks,js)
               a(ism,js) = te
  265       continue
  270    continue
  275 continue
c                                  complete computation of alfa(n)
      i2 = n1*m
      i1 = i2+p
      if (r.eq.0)  go to 290
      do 285 is=1,r
         i1s = i1+is
         do 280 js=1,m
  280    a(i1s,js) = a(i1s,js)+del(is,i2+js)
  285 continue
c                                  l u decomposition of alfa(n) (column
c                                    pivoting)
  290 i4 = i1-m+1
      i1 = i2+1
      i3 = n*m-1
      i8 = i3+1
      do 295 i=1,m
  295 ic(n,i) = i
      do 330 i=i1,i3
         i0 = i-i2
         te = 0.
         do 300 j=i0,m
            y = abs(a(i,j))
            if (y.le.te)  go to 300
            te = y
            ipiv = j
  300    continue
         if (te.eq.0.)  go to  335
c                                  column interchanges
         if (ipiv.eq.i0)  go to 315
         ite = ic(n,i0)
         ic(n,i0) = ic(n,ipiv)
         ic(n,ipiv) = ite
         do 305 is=i1,i8
            te = a(is,i0)
            a(is,i0) = a(is,ipiv)
            a(is,ipiv) = te
  305    continue
         do 310 is=i4,i2
            te = c(is,i0)
            c(is,i0) = c(is,ipiv)
            c(is,ipiv) = te
  310    continue
c                                  factorization
  315    ak = 1./a(i,i0)
         i01 = i0+1
         i7 = i+1
         do 325 is=i7,i8
            xmu = a(is,i0)*ak
            a(is,i0) = xmu
            do 320 js=i01,m
  320       a(is,js) = a(is,js)-xmu*a(i,js)
  325    continue
  330 continue
c
      return
  335 sing = .true.
      return
c
      end

      subroutine dvcpy (a,c,del,y,m,n,p,r,ir,ic,u,mtnmax,mmax,nmax,x)
c
c                                  specifications for arguments
c
      integer            m,n,p,r,mtnmax,mmax,nmax,ir(nmax,mmax),ic(nmax,
     .                   mmax)
      real               a(mtnmax,mmax),c(mtnmax,mmax),del(mmax,mtnmax),
     .                   y(mtnmax),u(mtnmax),x(mtnmax)
c
c                                  specifications for local variables
c
      integer            i0j,i0,i11,i1,i2,ii,ixn,ixp,ix,i,j1,jm,j,k1,k,
     .                   m1,mn,n1m,n1,n2,p1,p2
      real               te
c                                  first executable statement
      m1 = m-1
      n2 = n+1
      n1 = n-1
      p1 = p-1
      p2 = p+1
c                                  row interchanges on the right hand
c                                    side
      if (p.eq.0)  go to 10
      do 5 i=1,p
    5 x(i) = y(i)
   10 mn = m*n1
      do 15 i=1,m
   15 x(mn+i) = y(mn+i)
      do 25 i=1,n1
         ix = (i-1)*m+p
         do 20 j=1,m
            ixp = ix+ir(i,j)
            x(ix+j) = y(ixp)
   20    continue
   25 continue
c                                  solve l * y = x
      do 40 i=2,n
         i1 = m*(i-2)+1
         i2 = i1+p1
         i11 = i1-1
         if (p.eq.0)  go to 45
         do 35 j=i1,i2
            jm = j+m
            te = x(jm)
            do 30 k=1,m
   30       te = te-c(j,k)*x(i11+k)
            x(jm) = te
   35    continue
   40 continue
c
   45 n1m = n1*m
      ixn = n1m+p
      if (r.eq.0)  go to 60
      do 55 ii=1,r
         i1 = ixn+ii
         te = x(i1)
         do 50 j=1,n1m
   50    te = te-del(ii,j)*x(j)
         x(i1) = te
   55 continue
c                                  solve u * z = y
   60 do 100 i=1,n
         ii = n2-i
         i0 = (ii-1)*m
         i1 = i0+2
         i2 = ii*m
         if (i.eq.1)  go to 75
         do 70 j=p2,m
            i0j = i0+j
            te = x(i0j)
            do 65 k=1,m
   65       te = te-c(i0j,k)*x(i2+k)
            x(i0j) = te
   70    continue
   75    do 85 j=i1,i2
            k1 = j-i1+1
            te = x(j)
            do 80 k=1,k1
   80       te = te-a(j,k)*x(i0+k)
            x(j) = te
   85    continue
         x(i2) = x(i2)/a(i2,m)
         do 95 j=1,m1
            j1 = i2-j
            te = x(j1)
            k1 = m-j+1
            do 90 k=k1,m
   90       te = te-a(j1,k)*x(i0+k)
            x(j1) = te/a(j1,j1-i0)
   95    continue
  100 continue
c                                  interchanges in x
      do 110 i=1,n
         ix = (i-1)*m
         do 105 j=1,m
            ixp = ix+ic(i,j)
            u(ixp) = x(ix+j)
  105    continue
  110 continue
c
      return
c
      end

      subroutine gauleg (x, w, n)
c
      implicit none
c
      integer n
      real x(n),w(n)
c     returns the abscissas and weights for n-point
c     gauss-legendre integration
c ------------------------------------------------------------------------------
      real pi,eps
      parameter (pi=3.141592654,eps=3.e-14)
      integer i,j,m
      real z,p1,p2,p3,pp,z1
c
      m=(n+1)/2
      do i=1,m
         z=cos(pi*(i-0.25)/(n+0.5))
 1       continue
         p1=1.
         p2=0.
         do j=1,n
            p3=p2
            p2=p1
            p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j
         end do
         pp=n*(z*p1-p2)/(z*z-1.)
         z1=z
         z=z1-p1/pp
         if (abs(z-z1) .gt. eps)  go to 1
         x(i)=-z
         x(n+1-i)=z
         w(i)=2./((1.-z*z)*pp*pp)
         w(n+1-i)=w(i)
      end do
c
      return
c
      end

      subroutine interv (t, nt, x, left, mflag)
c
      implicit none
c
      real t(*),x
      integer nt,left,mflag
c
c     returns mflag=0 and the index left for t(left) <= x < t(left+1)
c     mflag=-1 left=1 if x < t(1)
c     mflag= 1 left=j where t(j) < t(nt) if t(nt) = x
c     mflag= 2 left=nt+1if t(nt) < x
c ------------------------------------------------------------------------------
c
      integer ilo,ihi,middle,j
c
      if (x .lt. t(1)) then
         left=1
         mflag=-1
         return
      else if (x .eq. t(nt)) then
         do j=nt-1,1,-1
            if (x .gt. t(j))  go to 5
         end do
 5       continue
         left=j
         mflag=1
         return
      else if (x .gt. t(nt)) then
         left=nt+1
         mflag=2
      else
         mflag=0
         ilo=0
         ihi=nt+1
 10      if (ihi-ilo .gt. 1) then
            middle=(ilo+ihi)/2
            if (x .ge. t(middle)) then
               ilo=middle
            else
               ihi=middle
            end if
            go to 10
         end if
      end if
c
      left=ilo
      return
c
      end

      function isamax1 (n, sx, incx)
c
      implicit none
c
      integer isamax1,n,incx
      real sx(1)
c     returns the index of element having the max absolute value
c ------------------------------------------------------------------------------
      integer ix,i
      real smax
c
      isamax1=0
      if (n .lt. 1)  return
      isamax1=1
      if (n .eq. 1)  return
c
      if (incx .ne. 1) then
         ix=1
         smax=abs(sx(1))
         ix=ix+incx
         do i=2,n
            if (abs(sx(ix)) .gt. smax) then
               isamax1=i
               smax=abs(sx(ix))
            end if
            ix=ix+incx
         end do
         return
      end if
c
      smax=abs(sx(1))
      do i=2,n
         if (abs(sx(i)).gt.smax) then
            isamax1=i
            smax=abs(sx(i))
         end if
      end do
      return
c
      end

      function qgauleg (func, a, b, n, ifail)
c
      implicit none
c
      real qgauleg,func,a,b
      integer n,ifail
      external func
c
c     returns n-point gauss-legendre quadrature of 'func'
c     from a to b.  n should be less than nmax =64
c ------------------------------------------------------------------------------
c
      integer nmax
      parameter (nmax=64)
      integer stderr
      parameter (stderr=6)
      integer nn,ns,j
      real x(nmax),w(nmax),xx(nmax),y(nmax),ss,xm,xr
      character errmsg*80
      save ns,x,w,errmsg
      data ns/-1/
      data errmsg/'warning from qgauleg:  n is set to nmax =64'/
c
      nn=n
      if (nn .gt. nmax) then
         nn=nmax
         ifail=1
         write (stderr,'(a)')errmsg
      end if
      if (ns .ne. nn) then
         ns=nn
         call gauleg(x,w,nn)
      end if
c
      xm=0.5*(b+a)
      xr=0.5*(b-a)
      do j=1,nn
         xx(j)=xm+xr*x(j)
         y(j)=func(xx(j))
      end do
      ss=0.
      do j=1,nn
         ss=ss+w(j)*y(j)
      end do
      qgauleg=ss*xr
      return
c
      end

      function sasum1 (n, sx, incx)
c
      implicit none
c
      real sasum1,sx(1)
      integer n,incx
c     returns the sum of absolute values
c     see linpack users' guide pD.2
c ------------------------------------------------------------------------------
      real stemp
      integer nincx,m,mp1
      integer i
c
      sasum1=0.
      stemp=0.
      if (n .le. 0)  return
      if (incx .ne. 1) then
         nincx=n*incx
         do i=1,nincx,incx
            stemp=stemp+abs(sx(i))
         end do
         sasum1=stemp
         return
      end if
c
      m=mod(n,6)
      if (m .ne. 0) then
         do i=1,m
            stemp=stemp+abs(sx(i))
         end do
         if (n .lt. 6) then
            sasum1=stemp
            return
         end if
      end if
      mp1=m+1
      do i=mp1,n,6
         stemp=stemp+abs(sx(i))+abs(sx(i+1))+abs(sx(i+2))
     .        +abs(sx(i+3))+abs(sx(i+4))+abs(sx(i+5))
      end do
      sasum1=stemp
      return
c
      end

      subroutine saxpy1 (n, sa, sx, incx, sy, incy)
c
      implicit none
c
      integer n,incx,incy
      real sa,sx(1),sy(1)
c     returns a constant times a vector plus a vector
c     see linpack users' guide pD.3.
c ------------------------------------------------------------------------------
      integer ix,iy,m,i,mp1
c
      if (n .le. 0)  return
      if (sa .eq. 0.)  return
c
      if (incx .eq. 1 .and. incy .eq. 1)  go to 20
      ix=1
      iy=1
      if (incx .lt. 0)ix=(-n+1)*incx+1
      if (incy .lt. 0)iy=(-n+1)*incy+1
      do i=1,n
         sy(iy)=sy(iy)+sa*sx(ix)
         ix=ix+incx
         iy=iy+incy
      end do
      return
 20   continue
c
      m=mod(n,4)
      if (m .ne. 0) then
         do i=1,m
            sy(i)=sy(i)+sa*sx(i)
         end do
         if (n .lt. 4)  return
      end if
      mp1=m+1
      do i=mp1,n,4
         sy(i)=sy(i)+sa*sx(i)
         sy(i+1)=sy(i+1)+sa*sx(i+1)
         sy(i+2)=sy(i+2)+sa*sx(i+2)
         sy(i+3)=sy(i+3)+sa*sx(i+3)
      end do
      return
c
      end

      function sdot1 (n, sx, incx, sy, incy)
c
      implicit none
c
      real sdot1,sx(1),sy(1)
      integer n,incx,incy
c     returns the dot product od two vectors
c     see linpack users' guide pD.5.
c ------------------------------------------------------------------------------
      integer ix,iy,i,m,mp1
c
      sdot1=0.
      if (n .le. 0.)  return
      if (incx .ne. 1 .or. incy .ne. 1) then
         ix=1
         iy=1
         if (incx .lt. 0)incx=(-n+1)*incx+1
         if (incy .lt. 0)incy=(-n+1)*incy+1
         do i=1,n
            sdot1=sdot1+sx(ix)*sy(iy)
            ix=ix+incx
            iy=iy+incy
         end do
      else
         m=mod(n,5)
         if (m .ne. 0) then
            do i=1,m
               sdot1=sdot1+sx(i)*sy(i)
            end do
            if (n .lt. 5)  return
         end if
         mp1=m+1
         do i=mp1,n,5
            sdot1=sdot1+sx(i)*sy(i)+sx(i+1)*sy(i+1)+sx(i+2)*sy(i+2)
     .           +sx(i+3)*sy(i+3)+sx(i+4)*sy(i+4)
         end do
      end if
      return
c
      end

      real function sdsdot1 (n, x, incx, y, incy, c)
c
      implicit none
c
      integer  n, incx, incy, i
      real     x(incx,*), y(incy,*), c, sum
c
      sum = 0.0
      if (n .gt. 0) then
        do i=1,n
          sum = sum + x(1,i) * y(1,i)
        end do
      end if
      sum     = sum + c
      sdsdot1 = sum
      return
c
      end

      subroutine sgefa1(a,lda,n,ipvt,info)
c
      implicit none
c
      integer lda,n,ipvt(1),info
      real a(lda,1)
c
c     uses function 'isamax1'
c
      integer isamax1
c
c     uses subroutines 'sscal1','saxpy1'
c     performs LU decomposition using gaussian elimination with
c     partial pivoting
c     see linpack users' guide pC.3
c ------------------------------------------------------------------------------
c
      integer nm1,k,l,kp1,j
      real t
c
      if (n .le. 0) then
         call STOP ('subroutine SGAFE1: invalid input n <= 0', 215)
      end if
      if (n .eq. 1) then
         info=0
         ipvt(1)=1
         return
      end if
c
      info=0
      nm1=n-1
      do k=1,nm1
         l = isamax1 (n-k+1,a(k,k),1)+k-1 ! find pivot element in column
         ipvt(k)=l
         if (a(l,k) .eq. 0.) then
            info = k                     ! column already triangularized
         else
            if (l .ne. k) then
               t=a(l,k)
               a(l,k)=a(k,k)
               a(k,k)=t
            end if
            t=-1./a(k,k)
            call sscal1 (n-k, t, a(k+1,k), 1) ! store column of L-matrix
            kp1=k+1
            do j=kp1,n         ! do row elimination with column indexing
               t=a(l,j)
               if (l .ne. k) then
                  a(l,j)=a(k,j)
                  a(k,j)=t
               end if
               call saxpy1(n-k,t,a(k+1,k),1,a(k+1,j),1)
            end do
         end if
      end do
c
      ipvt(n)=n
      if (a(n,n) .eq. 0.)info=n
      return
c
      end

      subroutine sgeim1 (a,lda,n,b,x,relerr,tol,aa,r,ipvt,info)
c
      implicit none
c
      integer lda,n,info
      real a(lda,1),b(n),x(n),relerr,tol,aa(lda,1),r(n)
      integer ipvt(n)
c     uses subroutines 'sgefa1','sgesl1'
c     uses functions 'sdsdot1','sasum1'
      real sdsdot1,sasum1
c     solves a.x = b using iterative improvement
c     aa and r are working arrays.
c     relerr returns an estimate of relative error in the solution
c     obtained by 'sgesl1'.
c     tol specifys the tolerance; above it the iterative improvement
c     procedure will be preformed.
c     ipvt contains pivot information, Specifically,
c     ipvt(k) is the index of the k-th pivot row.
c     on a successful run, info returns 0.  A negative value of info
c     indicates iterative improvement process does not converge.  A positive
c     value indicates that 'sgefa' detects exact singularity.
c     see linpack users' guide p1.9.
c ------------------------------------------------------------------------------
c
      integer itermax
      parameter (itermax=20)
      real xnorm,rnorm,t,xnormp
      integer i,j,iter
c
      do j=1,n
         do i=1,n
            aa(i,j)=a(i,j)
         end do
         x(j)=b(j)
      end do
c
      call sgefa1(aa,lda,n,ipvt,info)
      if (info .ne. 0) then    ! failure in sgefa1
         return
      end if
      call sgesl1(aa,lda,n,ipvt,x,0)    ! 0 for solving the aa.x = b
      xnorm = sasum1 (n, x, 1)
      if (xnorm .eq. 0.0)  return       ! is not required to improved
c
      relerr=0.
      xnormp=xnorm+tol
      do iter=1,itermax
         do i=1,n
            r(i)=sdsdot1(n,a(i,1),lda,x(1),1,-b(i))
         end do
         call sgesl1(aa,lda,n,ipvt,r,0)
         do i=1,n
            x(i)=x(i)-r(i)
         end do
         rnorm=sasum1(n,r,1)
         if (iter .eq. 1)relerr=rnorm/xnorm
         t=xnorm+rnorm
         if (t .le. xnormp)  return
      end do
c
      info=-1
      return
c
      end


      subroutine sgesl1(a,lda,n,ipvt,b,job)
c
      implicit none
c
      integer lda,n,ipvt(1),job
      real a(lda,1),b(1)
c     uses function 'sdot1'
      real sdot1
c     uses subroutine 'saxpy1'
c     solves the real system a.x = b (job=0),
c     or a(transpose).x = b (job .not.=0).
c     see linpack users' guide p.C.7.
c ------------------------------------------------------------------------------
      integer nm1,k,l,kb
      real t
c
      nm1=n-1
      if (job .eq. 0) then
         if (nm1 .ge. 1) then
            do k=1,nm1    ! solve L.y = b
               l=ipvt(k)
               t=b(l)
               if (l .ne. k) then
                  b(l)=b(k)
                  b(k)=t
               end if
               call saxpy1(n-k,t,a(k+1,k),1,b(k+1),1)
            end do
         end if
         do kb=1,n    ! solve U.x = y
            k=n+1-kb
            b(k)=b(k)/a(k,k)
            t=-b(k)
            call saxpy1(k-1,t,a(1,k),1,b(1),1)
         end do
      else                        ! solve a(transpose).x = b
         do k=1,n                 ! solve U(transpose).t = b
            t=sdot1(k-1,a(1,k),1,b(1),1)
            b(k)=(b(k)-t)/a(k,k)
         end do
         if (nm1 .lt. 1)  return
         do kb=1,nm1
            k=n-kb
            b(k)=b(k)+sdot1(n-k,a(k+1,k),1,b(k+1),1)
            l=ipvt(k)
            if (l .ne. k) then
               t=b(l)
               b(l)=b(k)
               b(k)=t
            end if
         end do
      end if
      return
c
      end

      subroutine sscal1(n,sa,sx,incx)
c
      implicit none
c
      integer n,incx
      real sa,sx(1)
c     scales a vector by a constant
c     see linpack users' guide pD.10.
c ------------------------------------------------------------------------------
      integer i,nincx,m,mp1
c
      if (n .le. 0)  return
      if (incx .ne. 1) then
         nincx=n*incx
         do i=1,nincx,incx
            sx(i)=sa*sx(i)
         end do
         return
      end if
c
      m=mod(n,5)
      if (m .ne. 0) then
         do i=1,m
            sx(i)=sa*sx(i)
         end do
         if (n .lt. 5)  return
      end if
      mp1=m+1
      do i=mp1,n,5
         sx(i)=sa*sx(i)
         sx(i+1)=sa*sx(i+1)
         sx(i+2)=sa*sx(i+2)
         sx(i+3)=sa*sx(i+3)
         sx(i+4)=sa*sx(i+4)
      end do
      return
c
      end

