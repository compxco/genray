c
c***************************************************************
c     This set of subroutines was coded up by Cary B. Forest   *
c     using the hot plasma dispersion tensor in Stix's book,   *
c     (Thomas Howard Stix, "Waves in Plasmas",                 *
c      Publ. by Amer. Inst. of Physics, 1992)                  *
c     and improved by A.P. Smirnov, 1999.                      *
c***************************************************************
c

c      double precision FUNCTION RTBIS(FUNC,X1,X2,XACC)
      !implicit none
c      double precision func,x1,x2,dx,xmid,fmid,xacc,f
c      integer*4 j,jmax
c      external func
c      PARAMETER (JMAX=40)
c      FMID=FUNC(X2)
c      F=FUNC(X1)
c      IF(F*FMID.GE.0.d0) PAUSE 'Root must be bracketed for bisection.'
c      IF(F.LT.0.d0)THEN
c        RTBIS=X1
c        DX=X2-X1
c      ELSE
c        RTBIS=X2
c        DX=X1-X2
c      ENDIF
c      DO 11 J=1,JMAX
c        DX=DX*.5d0
c        XMID=RTBIS+DX
c        FMID=FUNC(XMID)
c        IF(FMID.LE.0.d0)RTBIS=XMID
c        IF(dabs(DX).LT.XACC .OR. FMID.EQ.0.d0) RETURN
c11    CONTINUE
c      PAUSE 'too many bisections'
c      END

      SUBROUTINE CZETA(CX,CZ0,CZ1,IERR)
C
C ==========================================================
C
C  Plasma dispersion function of complex argument
C  and its first derivative (up to the 3d possible)
C
C  Input:      CX       - argument (complex)
C
C  Output:     CZ0      - Z(x)     (complex)
C              CZ1      - dZ/dx    (complex)
C              IERR     - failure test (integ.)
C
C              IERR = 1  -  convergence obtained
C              IERR = 0  -  convergence not obtained
C
C ==========================================================
C

      implicit none

      double precision H2N(30),A(7),ZSQPI,y,x,x2,xmodsq,xx,yy,h,q

      integer*4 N,jmax,kmax,nu,I,IERR
      DOUBLE COMPLEX CRE,CIM,CZERO,cw,cw2,cw1,cw3
      DOUBLE COMPLEX CX,CZ1,CZ0,CX2,CXQ,CQ
C
      DATA    CRE    /(1.0D0, 0.0D0)/,
     +        CIM    /(0.0D0, 1.0D0)/,
     +        CZERO  /(0.0D0, 0.0D0)/
C
      DATA ZSQPI /1.7724538509D0/
C
      DATA A  /6.0D0,              30.0D0,            157.5D0,
     +      9.45D+02,        6.496875D+03,    5.0675625D+04,
     +           4.4341171875D+05 /
C
C ==========================================================
C  Preliminaries
C
      IERR = 1
      X = DREAL(CX)
      Y = DIMAG(CX)
      X2 = X*X
      XMODSQ = X2 + Y*Y
      CW = CZERO
      CW2 = -CIM*CX
C
      IF(XMODSQ.LT.35.7D0)  THEN
C
C  Gautschi Algorithm
C
         CW1 = CZERO
         X = DABS(X)
         Y = DABS(Y)
C
         IF(X.LT.5.33 .AND. Y.LT.4.29)  THEN
C
            XX = X/5.33
            YY = Y/4.29
            H = (1.D0-YY)*DSQRT(1.D0-XX*XX)
            Q = 6.D0 + 23.D0*H
            N = IDINT(Q + 0.5D0)
            Q = 9.D0 + 21.D0*H
            NU = IDINT(Q + 0.5D0)
            JMAX = NU - N
            KMAX = N + 1
            N = NU + 2
            H = 1.6D0*H
C
            IF(Y.LT.0.D0)  THEN
               H = -H
               CW = ZSQPI*CDEXP(-CX*CX)
            END IF
C
            CW2 = CW2 + H*CRE
            CW3 = CZERO
C
            H = 2.D0*H
            H2N(1) = 1.
            DO 10  I=2,KMAX
   10          H2N(I) = H2N(I-1)*H
C
            DO 20  I=1,JMAX
               N = N-1
   20          CW3 = 0.5D0/(CW2+DFLOAT(N)*CW3)
C
            DO 30  I=1,KMAX
               N = N-1
               CW3 = 0.5D0/(CW2+DFLOAT(N)*CW3)
   30          CW1 = CW3*(H2N(N)+CW1)
C
         ELSE
C
            N = 10
            DO 40  I=1,9
               N = N-1
   40          CW1 = 0.5D0/(CW2+DFLOAT(N)*CW1)
C
         END IF
C
         CZ0 = 2.D0*(CIM*CW1+CW)
C
         IF(Y.EQ.0.D0)
     >      CZ0 = DCMPLX(DBLE(CZ0),ZSQPI*DEXP(-X2))
C
         CZ1 = -2.D0*(CRE+CX*CZ0)
C
      ELSE
C
C  Asymptotic expansion
C
         CX2 = -CX*CX
         CQ = CX2 + 0.5D0*CRE
         CXQ = 2.D0*CX*(CQ+CRE)
C
         IF(Y.LT.0.)  THEN
C
            XX = DBLE(CX2)
            YY = DIMAG(CX2)
            IF(XX.GT.100.D0  .OR. YY.GT.1.D14)  THEN
               IERR = 0
               RETURN
            END IF
            IF(X.GT.-30.D0)
     >         CW = 8.D0*ZSQPI*CDEXP(CX2)*CXQ*CIM
C
         ELSE IF(Y.EQ.0.)  THEN
            IF(X.GT.-30.D0)
     >         CW = 4.D0*ZSQPI*CDEXP(CX2)*CXQ*CIM
         END IF
C
         CW1 = -1./CX2
         NU = INT(35.7D0/DSQRT(XMODSQ) + 1.94D0)
         N = NU + 1
         CW3 = CZERO
C
         DO 110  I=1,NU
            N = N - 1
  110       CW3 = CW1*(A(N)+CW3)
C
         CW3 = CW3*CW1 + CW
         CW2 = (4.D0*CRE-CQ*CW3)/CXQ
         CZ1 = (0.5D0*CX*CW2-CRE)/CQ
         CZ0 = -(CRE+0.5D0*CZ1)/CX
C
      END IF
C
      RETURN
      END SUBROUTINE CZETA

      subroutine CZETA0(n_parallel,CX,CZ0,CZ1,IERR)
c     Z_0=CZ0 and its derivative d(Z_0)/d(CX)=CZ1 calculations
c     using the plasma dispersion function Z(cx) (Stix book p. 202 (82))
      implicit none
c     input
      double precision n_parallel      
      double complex CX
c     output
      integer*4 IERR
      double complex CZ0,CZ1
c     local
      double complex cp,cpz0,cpz1
      integer j, jmax

      double complex CZ2,cpz2
      double precision X

c     used:  CZETA(CX,CZ0,CZ1,IERR)

c-----CZETA is from forest code
c      if (n_parallel.ge.0.d0) then
c        call CZETA(CX,CZ0,CZ1,IERR)
c        if (IERR.eq.0) then
c           write(*,*)'CZETA IERR=0 npar >0 CX=',CX
c        endif
c      else
c        cp=-CX
c        call CZETA(cp,cpz0,cpz1,IERR) 
c        if (IERR.eq.0) then
c           write(*,*)'CZETA IERR=0 npar <0 cp=(-CX)=',cp
c        endif
c        CZ0=-cpz0
c        CZ1=cpz1
c      endif
c      write(*,*)'forest CZ0,CZ1',CZ0,CZ1

c-----zfun  is from torray code
c      if (n_parallel.ge.0.d0) then
c        call zfun (CX,CZ0)
c        CZ1=-2.d0*(1.d0+CX*CZ0)
c      else
c        cp=-CX
c        call zfun (cp,cpz0)
c        CZ0=-cpz0
c        CZ1=-2.d0*(1.d0+CX*CZ0)
c      endif

c-----zfunc and pfunc from Brambilla LH/FW code
c      complex w,z,zp
c      w=cmplx(CX)
c      if (n_parallel.ge.0.d0) then
c        call zfunc(w,z,zp)        
c        CZ0=z
c        CZ1=zp
c      else
c        call zfunc(-w,z,zp)     
c        CZ0=-z
cc test       CZ1=-2.d0*(1.d0+CX*CZ0)
cc        write(*,*)'w,CZ0,CZ1,zp',w,CZ0,CZ1,zp
c        CZ1=zp
c      endif

c-----zfunc and pfunc from curray code
c      jmax=200
c      do j=1,jmax
c         x=-5.+ 10.*(j-1.)/(jmax-1.)
c         call ZFUN_cur(x, CZ0, CZ1,CZ2) 
c         write(*,'(5e12.3)')x,CZ0,CZ1
c      enddo
c      pause
      
       ! x=dble(CX)  ! YuP 120430 added; skip n_parallel part below
       ! call ZFUN_cur(X, CZ0, CZ1, CZ2)   ! CZ2 is not used here 
       ! return

      if (n_parallel.ge.0.d0) then  
        x=dble(CX)
        call ZFUN_cur(X, CZ0, CZ1, CZ2)   ! CZ2 is not used here 
      else
        cp=-CX
        x=dble(cp)
        call ZFUN_cur(x, cpz0, cpz1, cpz2)  ! cpz2 is not used here 
c        write(*,*)'forest x,cpz0,cpz1',x,cpz0,cpz1
        CZ0=-cpz0
        CZ1=cpz1
        !CZ2=cpz2 !?
      endif
cc      write(*,*)'curray CZ0,CZ1',CZ0,CZ1

      return
      end subroutine CZETA0



      subroutine DHOT_s(mass,X,Yc,T_av,tpop,vflow,
     .nll_in,np_in,iherm,K)
c Revised by YuP [2020-08-24]
c-----------------------------------------------------------
c     calculates sucseptibilities K(3,3) for the given hot (non-relativistic)
c     plasma specie, T.Stix, Waves in plasmas,(1992), p.258 (57)  make

c------------------------------------------------------------
c     INPUTS:
c      mass - the mass of the given specie (in electron mass) 
c      X = (fpe/f)**2
c      Y = fce/f ! for electron Ye should has the sign opposite to Yi
c      T_av=average temperature=Te(1+2*tpop)/3 in ev.
c           Here Te is the parallel temperature in eV  
c      tpop = T_perp/T_parallel  
c      vflow    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c      iherm =1 hermitian dielectric tensor, 2-full     
c     OUTPUT:
c      K(3,3):  the nine components of sucseptibilities
c               (for the given specie)
c          evaluated at (mass,X,Y,Te,tpop,vflow,Nll,np)
c-----------------------------------------------------------

      implicit none
c     input
      double precision mass
      double precision  X, Y, Yc, T_av, tpop, vflow
      double precision np_in,nll_in
      integer iherm      
      integer n        !n.lt.nharmmax is the number of ECR harmonics
      integer n_new
      integer nmax     ! max number of ECR harmonics in arrays
      parameter(nmax=300)
 
c     output
      double complex K(3,3)
c     locals
      double precision np,nll, nll_adj
      double precision nps,nlls,c
      double precision te, beta
      double precision pi,k4,vte,lambda,ri01, omega,rho_larm, one_lambda
      real*8 nll_min

      double precision dr,di,xne(nmax),drp, ox,resn, resn_min, vkw_min,
     + vkw_in, vkw_adj, err_vkw
      double precision Ibess(-(2*nmax+1):(2*nmax+1)),
     &Ibessp(-(2*nmax+1):(2*nmax+1)),
     &Ibesspp(-(2*nmax+1):(2*nmax+1))

      double complex   dp
      double complex A(-(2*nmax+1):(2*nmax+1)),B(-(2*nmax+1):(2*nmax+1))
      double complex i
      double complex  Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz

      double complex cx,cz0, cz0_vkw, cz0_resn_vkw, cz1,cd,cxp
      integer q,qmin,qmax,ierr,j

      integer n_harm0,n_harm1,n_harm2 ! the boundaries for the harmonic numbers
      integer istop   
      real*8 qIbess_lambda !local: q*Ibess/lambda
      
c      istop=1, it is possible to calculate the susceptibility tensor
c      istop=0, it is impossible to calculate the susceptibilty tensor
      pi = 4*datan(1.d0)
      c= 2.99792458d10          !speed of light
      k4 = 4.19d7               !constant for elec therm vel
      i = ( 0.0d0,1.0d0)        !imaginary number

c-----tpop =  ratio of Te_perp / Te_parallel 
c-----vflow = flow velocity of electrons
      te=   3.d0*T_av/(1.d0+2.d0*tpop) !longitudinal           
      vte=  k4*dsqrt(2.0d0*te/mass)    !vth= sqrt(2.0*kT/m) longitudinal
      beta= vte / c 
      
!YuP      call npnllmin(nll_in,np_in,nll,np) !in DHOT_s: Not called anymore.
      nll=  nll_in ! original Npar
      nll_min=1.d-7 !Must be a positive value
      ! YuP[11-2016]: Adjust the value of nll:
      if (dabs(nll_in).lt.nll_min) then
         if (nll_in.ge.0.d0) then
            nll= nll_min
         else 
            nll=-nll_min
         endif 
         !Adjusted, with lower-limit value of |nll|, 
         !to avoid a jump in plasma dispersion function at resonance  
         !(jump in Real part of CZ0) when Npar~0.
         !However, for calculation of damping (Im part of CZ0),  
         !we will use the original nll_in, i.e. using vkw_in, see below. 
      endif
      
      np=   max(0.d0,np_in)
      !YuP[2020-09-02] Note: this subr.DHOT_s is called by func.dhot_sum,
      ! which may use negative np (perpendicular refr.index) values.
      ! See subr.hotnp which searches the root of hot dispersion determinant,
      ! and it may use a range of np from negative values.
      ! In original version, np would be adjusted by npnllmin() to a small
      ! positive value (1.d-7). Now, we adjusted it to be a non-negative value.
      
      nps=  np**2   ! nperp^2
      nlls= nll**2

      Y=Yc  ! Y=omega_c/omega == n_harm0
      ! Note: B=0 points exist in FRC, so that Y can approach 0.
      ! But if ray was not absorbed at multiple resonances
      ! when approaching to B=0 point, 
      ! it would not be absorbed anyway, 
      ! so just pass through such B=0 point, using nmax harmonics.
      if (abs(Y).lt. 1.d0/dfloat(nmax-4)) then 
         !YuP120430 omega/omega_ce > nmax   Need more harmonics
         !omega=2*pi*60.d9
         !rho_larm= vte/(Yc*omega) ! [cm] Larmor radius of ion species
         lambda= 0.5d0*tpop*nps*(beta/Yc)**2 ! 0.5 k_perp^2 * rho_gyro^2
         write(*,'(a,3e12.3)')
     +    'DHOT_s: wce/w, lambda, 1/(nmax-4) =', 
     +     Yc, lambda, 1.d0/dfloat(nmax-4)
         write(*,*)'DHOT_s: omega/omega_ce >nmax-4  Need more harmonics'
         write(*,*)'DHOT_s: POSSIBLE B~0 POINT. ADJUSTING Y to 1/nmax'
         !stop 'DHOT_s:  Increase nmax'
         !write(*,*)'in DHOT_s before Y-adjusted:   Y,nmax=',Y,nmax
         Y= dsign(1.d0/dfloat(nmax-4), Y)
         !write(*,*)'in DHOT_s after Y-adjusted:    Y,nmax=',Y,nmax
         !pause !!!
      endif
  
     
c-----lambda= 0.5 k_perp^2 * rho_gyro^2
      if(abs(Yc).ne.0.d0)then
        lambda= 0.5d0*tpop*nps*(beta/Yc)**2 ! 0.5 k_perp^2 * rho_gyro^2
      else ! Yc=0 which means B=0 null point.  Use "adjusted" Y:
        lambda= 0.5d0*tpop*nps*(beta/Y)**2 
      endif
      
      if(abs(tpop*nps*beta**2).ne.0.d0)then ! Basically, when lambda>0
        one_lambda= Yc**2/(0.5d0*tpop*nps*beta**2) !=1/lambda, used below
      else ! could happen at far edge, where Vte=0 and so beta=0
        !When lambda=0, use this:
        one_lambda=0.d0 !Be careful how it is used, see below.
      endif
      
      vkw_in= nll_in*beta  ! same as Vthermal*Kpar/omega, using original nll
      vkw_adj= nll*beta ! Adjusted, with lower-limit value of |nll|, 
      !to avoid a jump in plasma dispersion function (Real part of CZ0) 
      !when Npar~0.
      !However, for calculation of damping (Im part of CZ0), we will use 
      !the original nll_in, i.e. using vkw_in, see ZFUN_vkw below. 
      
      ! Initialize:
      Kxx = (0.0D0,0.0D0)
      Kyy = (0.0D0,0.0D0)
      Kzz = (0.0D0,0.0D0)
      !if(abs(nll*beta).ge. 1.d-2)then
      !Kzz = dcmplx(vflow/(c*tpop), 0.0D0) *2*X/(nll*beta**2)
      ! YuP: this term (vflow/c) cancels with 
      ! (-vflow/c) in B0 (q=0) term
      !endif
      Kxy = (0.0D0,0.0D0)
      Kxz = (0.0D0,0.0D0)
      Kyz = (0.0D0,0.0D0)
     
c   Now, calculate the modified bessel functions 
c   Remember that these are actually exp(-lambda) * In(lambda)
      do j = 1,nmax
         xne(j) = 0.0D0
      enddo
      do j = -(2*nmax+1),(2*nmax+1)
         ibess(j) = 0.0D0
         ibessp(j) = 0.0D0
         ibesspp(j) = 0.0D0
         A(j) = dcmplx(0.0D0,0.0D0)
         B(j) = dcmplx(0.0D0,0.0D0)
      enddo

c-----calculations of the max and minimal numbers of the harmonics
!      istop=1 ! initialize
!      call harmon_z(X,Y,nll,vflow,c,beta,lambda,nmax,
!     .n_harm1,n_harm2,n,istop)
!
!      if (istop.eq.0) then
!         write(*,*)'DHOT_s istop=0 it is impossible to calculate' 
!         write(*,*)'the harmonics of susceptibily tensor  1/Y=',1/Y,nll
!         write(*,*)'n_harm1,n_harm2,n,nmax =',n_harm1,n_harm2,n,nmax
!         pause !!!  
!      endif      

      ! YuP[Oct-2014] 
      ! Alternative definition of summation range over harmonics
      n_harm0= ceiling(1.d0/Y)  ! can be negative
      n_harm0= min(n_harm0,nmax-2) 
      if(lambda.le. 0.4d0)then
         ! small lambda:  0 < lambda < 0.4 
         n_harm2=  iabs(n_harm0) + 3
         if(n_harm2.gt.nmax)then
            write(*,'(a,e12.3,i7)')
     +       'DHOT_s: n_harm2>nmax;  lambda,n_harm2=',lambda,n_harm2
            !pause !!!
         endif
         n_harm2=  min(n_harm2,nmax-1) ! not to exceed nmax-1
         n_harm1= -n_harm2
      elseif(lambda.le.10.d0)then
         ! 0.4 < lambda < 10.
         n_harm2=  iabs(n_harm0) + 8
         if(n_harm2.gt.nmax)then
            write(*,'(a,e12.3,i7)')
     +       'DHOT_s: n_harm2>nmax;  lambda,n_harm2=',lambda,n_harm2
            !pause !!!
         endif
         n_harm2=  min(n_harm2,nmax-1) ! not to exceed nmax-1
         n_harm1= -n_harm2
      else
         ! large lambda: lambda > 10.
         n_harm2=  iabs(n_harm0) + int(3*sqrt(lambda))
         if(n_harm2.gt.nmax)then
            write(*,'(a,2e12.3,i15)')
     +       'DHOT_s: n_harm2>nmax;  omega/omega_c, lambda, n_harm2=',
     +                                 1.0/Y,       lambda, n_harm2
            !pause !!!
         endif
         n_harm2= min(n_harm2,nmax-1) ! not to exceed nmax-1
         !n_harm2= min(n_harm2,20) ! not to exceed 20
         n_harm2= max(n_harm2, iabs(n_harm0)+3) ! but not less than w/wc
         n_harm1= -n_harm2
      endif
      !Tests O-X-EBW, with cases of large lambda (up to lambda=160.) 
      ! and w/wc ~ 2 (edge) to 1(core):  
      ! results are same as with original harmon_z method,
      ! but the range is much smaller (not more than n_harm2=40).
          
c yup      if (istop.eq.0) goto 10
!      n_new=n
!      do j=1,n 
!        if (dabs(xne(j)).lt.1.d-14) then
!          !skip terms that are too small because of small exp(-lambda)*In
!           istop=1 ! initialize
!           call harmon_z(X,Y,nll,vflow,c,beta,lambda,j+2,
!     .     n_harm1,n_harm2,n_new,istop)
!           goto 9          
!        endif
!      enddo
! 9    continue
!      n=n_new


 10   continue

c      write(*,*)'dhot_s j, xne(j)',j,xne(j)
c      write(*,*)'dhot_s n,n_harm1,n_harm_2',n,n_harm1,n_harm2
c      write(*,*)'dhot_s before goto 20 istop=',istop     
cyup      if (istop.eq.0) goto 20
      ! Be sure to include q=0 harmonic (plus 2 extra):
      qmin=min(-2,n_harm1)
      qmax=max(n_harm2,2)
      n= max(iabs(qmin),iabs(qmax))

c      write(*,'(a, 5e13.4,3i6)') 
c     + 'DHOT_s nper2,beta*nll,X,Y,lambda, qmin,qmax,n',
c     +           nps,beta*nll,X,Y,lambda, qmin,qmax,n
           
      if (n+1.gt.nmax) then
         write(*,*)'DHOT_s:  n+1>nmax' 
         write(*,*)'qmin,qmax,n+1,nmax =',qmin,qmax,n+1,nmax
         stop !!!  
      endif      
      
      call ibess0(lambda,ri01)
      call ibessn(lambda,n+1,ri01,xne) ! xne== In(lambda)/exp(lambda) in DHOT_s

      err_vkw=1.d-6 !for checking:  sqrt(resn*resn+vkw*vkw) .LE. err_vkw
      !
      do q  = qmin,qmax
         resn= 1.d0 -q*Yc - nll*vflow/c   ! resonance term  
         !--------------
         if (q.eq.0) then
            Ibess(q+n+1)  = ri01   !== I0/exp(lambda)
            Ibessp(q+n+1) = xne(1) !== I0'/exp() = I1/exp()
c            Ibesspp(q+n+1)= xne(2)+xne(1)/lambda ! Not used?
            !YuP[2020-09-02] Below, we will use a combination q*I_q/lambda.
            !For q=0, it is simply 0 at all lambda:
            qIbess_lambda=0.d0
         elseif(iabs(q).eq.1)then
            Ibess(q+n+1)  = xne(iabs(q))
            !-- Special care for the derivative I1' when lambda-->0.
            !Note: one_lambda== 1/lambda= Yc**2/(0.5d0*tpop*nps*beta**2)
            !      and it was set to be 0 when lambda=0.
            !YuP[2020-09-02] Note that for nps=0 (Nperp=0), 
            ! we get lambda=0. In this case, for I1'
            ! we should use asymptotic expansion of I1(lambda):
            !  I1'(x) = I2(x) + I1(x)/x = [x-->0] = 
            !         = 0     + (x/2)/x == 1/2
            if(lambda.lt.1.d-7)then
              Ibessp(q+n+1)=0.5d0
              qIbess_lambda=q*0.5d0 ! q*I_q(lambda)/lambda  at lambda-->0
            else ! lambda>0, which means one_lambda is properly set
              ! Use the usual I_n' = I_{n+1} + I_n*(n/x)
              Ibessp(q+n+1) = xne(1+iabs(q))
     .                      + iabs(q)* xne(iabs(q)) *one_lambda
              qIbess_lambda= q*Ibess(q+n+1)*one_lambda
            endif
         else ! q=+/-2 or higher
            Ibess(q+n+1)  = xne(iabs(q))
            Ibessp(q+n+1) = xne(1+iabs(q))
     .                    + iabs(q)* xne(iabs(q)) *one_lambda
            !Note that for higher q, e.g. q=2, we have at small x:
            ! I2' = I3(x) +I2(x)*(2/x) --> 
            !     =   0   +(x/2)^2/2! *(2/x) --> 0
            !So, using one_lambda->0 is justified.
            qIbess_lambda= q*Ibess(q+n+1)*one_lambda ! q*I_q(lambda)/lambda
         endif
         !vkw_adj== nll*beta, where nll was adjusted
         !vkw_in=nll_in*beta with original nll_in
         call ZFUN_vkw(resn, vkw_in,vkw_adj, CZ0_vkw, CZ0_resn_vkw) 
         ! no 1/npar divergence in CZ0_vkw or CZ0_resn_vkw
         !In case  |resn| > 10*|vkw|  including npar=0 case,
         ! cz0_vkw ~~ -1/resn
         ! cz0_resn_vkw ~ -1. (No divergence)
         A(q+n+1)= (tpop-1.d0) + ( CZ0_resn_vkw*tpop + CZ0_vkw*q*Yc )
         ! In this shape, A has no 1/npar divergence 
         ! (but may have 1/resn divergence from cz0_vkw)
         Kxx=  Kxx + X*q*qIbess_lambda * A(q+n+1) 
         Kxy=  Kxy -i*q*X *(Ibess(q+n+1)-Ibessp(q+n+1))*A(q+n+1)
         Kyy=  Kyy + X*( q*qIbess_lambda +
     .         2.d0*lambda*(Ibess(q+n+1)-Ibessp(q+n+1)) )*A(q+n+1)
!         write(*,*)'DHOT_s: np,q,n=',
!     &    np,q,n, Ibess(q+n+1), Ibessp(q+n+1), A(q+n+1)
         !if(abs(nll*beta).lt.abs(resn)*1.d-2)then ! 1d-2 or 1d-3 ~same result
         if(dsqrt(resn*resn+vkw_adj*vkw_adj) .LE. err_vkw)then
            ! Both npar~0 and resn~0 
            ! do nothing: assume no contribution to Kzz
            !write(*,*)'DHOT_s: sqrt(resn*resn+vkw*vkw)<err_vkw',resn,vkw_adj
            !pause
         elseif( abs(vkw_adj).lt.0.001*abs(resn) )then ! |vkw/resn|<<1, no damping
            !YuP case of nll->0, including nll=0.    
            !(in general, case of |vkw/resn|<<1 )
            ! The following is valid when also vflow=0:
            ox= vkw_adj/resn 
            B(q+n+1)=-( (resn+q*Yc)*vflow/c + 
     +                  (1.d0-q*Yc)*(tpop+q*Yc/resn)*0.5*ox*beta  )/resn 
            Kzz= Kzz+ Ibess(q+n+1)*(tpop+q*Yc/resn)*(-X)
            !if(abs(nll).lt.0.02 .and. abs(resn).lt.1.d-3)then
            !  write(*,*)'DHOT_s: abs(vkw)<0.001*abs(resn)',resn,vkw_adj,Kzz
            !endif
         else ! General case:  nll.ne.0  (resn can be any value here)
            ! Absorption becomes large when |resn/vkw|~1 or less
            ! which can be written as abs(vkw)>abs(resn).
            B(q+n+1)= (1.0 + CZ0_resn_vkw*(resn*tpop+q*Yc) 
     +                  + resn*(tpop-1.d0)   )/nll  
     +              + ( A(q+n+1)-1.d0 )*vflow/c
            Kzz= Kzz+(1.d0-q*Yc)*Ibess(q+n+1)*B(q+n+1)*2*X/(nll*beta**2)
            !if(abs(nll).lt.0.02 .and. abs(resn).lt.1.d-3)then
            !  write(*,*)'  DHOT_s: abs(vkw)>0.001*abs(resn)',resn,vkw
            !  write(*,*)'    DHOT_s:Kzz         ', Kzz
            !  write(*,*)'    DHOT_s:CZ0_vkw     ', CZ0_vkw
            !  write(*,*)'    DHOT_s:CZ0_resn_vkw', CZ0_resn_vkw
            !  !if(abs(vkw)*0.1>abs(resn)) pause
            !endif
         endif
         Kxz= Kxz +   X*np*  qIbess_lambda * B(q+n+1)/Yc
         Kyz= Kyz + i*X*np*(Ibess(q+n+1)-Ibessp(q+n+1)) * B(q+n+1)/Yc
      enddo ! q  = qmin,qmax
          
 20   continue ! handle for istop=0

   
 30   continue
 
      Kyx = -1.0D0* Kxy
      Kzx =  Kxz
      Kzy = -1.0D0* Kyz

      K(1,1) =  Kxx
      K(2,2) =  Kyy
      K(3,3) =  Kzz
      K(1,2) =  Kxy
      K(1,3) =  Kxz
      K(2,1) =  Kyx
      K(2,3) =  Kyz
      K(3,1) =  Kzx
      K(3,2) =  Kzy

c     Hermitian part
      if (iherm.eq.1) then
       Kxx = 0.5D0 * ( K(1,1) + dconjg(K(1,1) )) !kxx
       Kyy = 0.5D0 * ( K(2,2) + dconjg(K(2,2) )) !kyy
       Kzz = 0.5D0 * ( K(3,3) + dconjg(K(3,3) )) !kzz
       Kxy = 0.5D0 * ( K(1,2) + dconjg(K(2,1) )) !kxy
       Kxz = 0.5D0 * ( K(1,3) + dconjg(K(3,1) )) !kxz
       Kyz = 0.5D0 * ( K(2,3) + dconjg(K(3,2) )) !kyz
       Kyx = 0.5D0 * ( K(2,1) + dconjg(K(1,2) )) !kyx
       Kzx = 0.5D0 * ( K(3,1) + dconjg(K(1,3) )) !kzx
       Kzy = 0.5D0 * ( K(3,2) + dconjg(K(2,3) )) !kzy
      endif

      K(1,1) =  Kxx
      K(2,2) =  Kyy
      K(3,3) =  Kzz
      K(1,2) =  Kxy
      K(1,3) =  Kxz
      K(2,1) =  Kyx
      K(2,3) =  Kyz
      K(3,1) =  Kzx
      K(3,2) =  Kzy

      return
      end subroutine DHOT_s

      



      subroutine DKtens(X,Y,T_av,tpop,vflow,nll_in,np_in,  !Not called
     & K,dK,dd,ddn)
c
c     This function returns the derivatives of the Hermitian
c     part of the hot, nonrelativistic dielectric tensor
c     with respect to  X, Y, Te_av,tpop,vflow, nll_in, and np_in
c     and the derivative of the dispersion function wrt
c     X,Y,Te_av,nll_in, and np_in
c
c     INPUTS:
c
c      X = (fpe/f)**2
c      Y =  fce/f  !Y_e should have the sign opposite to Y_i  
c     T_av=average temperature=Te(1+2*tpop)/3 in ev
c     Te = parallel electron temperature in eV.   
c     tpop = T_perp/T_parallel  (future)
c     vflow (future)    cm/sec
c     K(3,3)
c     nll - parallel index of refraction n.
c     np - perpendicular index of refraction n.

c     OUTPUT:  
c
c     dK(*,*,1)  (the partial derivatives K wrt np
c     dK(*,*,2)  (the partial derivatives K wrt nll
c     dK(*,*,3)  (the partial derivatives K wrt  X)
c     dK(*,*,4)  (the partial derivatives K wrt  Y)
c     dK(*,*,5)  (the partial derivatives K wrt  Te_av) [1/eV]
c     dK(*,*,6)  (the partial derivatives K wrt  tpop)
c     dK(*,*,7)  (the partial derivatives K wrt  vflow) [sec/cm]

c
c     derivatives of the determinant of the hermitian part of the
c     dielectric tensor:
c
c     DD(1)  d Dh(w,k) / d np
c     DD(2)  d Dh(w,k) / d nll
c     DD(3)  d Dh(w,k) / d X
c     DD(4)  d Dh(w,k) / d Y
c     DD(5)  d Dh(w,k) / d Te_av
c     DD(6)  d Dh(w,k) / d tpop
c     DD(7)  d Dh(w,k) / d vflow
c
c     derivatives of the determinant of the dielectric tensor
c
c     DDN    d D(w,k) / d np   

      implicit none

      integer nmax     ! max number of ECR harmonics in arrays
      parameter(nmax=300)

      double precision np_in,nll_in,np,nll,npi
      double precision nps,nlls,c
      double precision  X, Y,T_av,te, Xp, Yp,beta,betap,nllp,npp
      double precision pi,k4,vte,lambda,lambdap,tpop,tpopp,ri01,
     .vflow,vflowp,betat_av,betatpop
cSAP090204
c      double precision dr,di,xne(nmax),drp,Ibpp(nmax)
c      double precision Ibess(nmax),Ibessp(nmax),Ibesspp(nmax),Ibp(nmax)
c      double complex A(nmax),B(nmax),Ap(nmax),Bp(nmax)
      double precision dr,di,xne(nmax),drp,Ibpp(-(2*nmax+1):(2*nmax+1))
      double precision Ibess(-(2*nmax+1):(2*nmax+1)),
     &Ibessp(-(2*nmax+1):(2*nmax+1)),
     &Ibesspp(-(2*nmax+1):(2*nmax+1)),Ibp(-(2*nmax+1):(2*nmax+1))
      double complex A(-(2*nmax+1):(2*nmax+1)),
     &B(-(2*nmax+1):(2*nmax+1)),Ap(-(2*nmax+1):(2*nmax+1)),
     &Bp(-(2*nmax+1):(2*nmax+1))


      double complex i,dK(3,3,7),K(3,3),dd(7),dp,ddn
      double complex dKt_av,dKtpop,ddt_av,ddtpop
      double complex kxx,kxy,kxz,kyx,kyy,kyz,kzx,kzy,kzz
      double complex kxxp,kxyp,kxzp,kyxp
      double complex kyyp,kyzp,kzxp,kzyp,kzzp
      double precision vflowdc
      double complex cx,cz0,cz1,cd,cxp,cz1_p
      integer*4 q,n,ierr,j,jj
cSm060816
      integer*4 n_new
        
      integer n_harm0,n_harm1,n_harm2 !the limits for the harmonic numbers
      integer istop

      pi = 4.d0*datan(1.d0)
      c= 2.99792458d10                !speed of light
      k4 = 4.19d7               !constant for elec therm vel
      i = ( 0.0d0,1.0d0)        !imaginary number
c      n = 20                   !max number of the cyclotron harmonics
      te=3.d0*T_av/(1.d0+2.d0*tpop) !longitudinal temperature
    
      vte =k4*dsqrt(2.0d0*te)         !vte = sqrt(2.0*kTe/me)
      beta = vte / c
      vflowdc=vflow/c

      call npnllmin(nll_in,np_in,nll,np) !in subroutine DKtens
     
      nps = np**2
      nlls = nll**2

c     lambda = 0.5 k_perp^2 * rho_e^2

      lambda = nps * beta**2 * tpop / 2.0d0 /  Y**2
cSAP090204
      do j = 1,nmax
         xne(j) = 0.0D0
c         ibess(j) = 0.0D0
c         ibessp(j) = 0.0D0
c         ibp(j) = 0.0D0
c         ibpp(j) = 0.0D0
c         A(j) = dcmplx(0.0D0,0.0D0)
c         Ap(j) = dcmplx(0.0D0,0.0D0)
c         B(j) = dcmplx(0.0D0,0.0D0)
c         Bp(j) = dcmplx(0.0D0,0.0D0)
      enddo
      do j = -(2*nmax+1),(2*nmax+1)
          ibess(j) = 0.0D0
          ibessp(j) = 0.0D0
          ibp(j) = 0.0D0
          ibpp(j) = 0.0D0
          A(j) = dcmplx(0.0D0,0.0D0)
          Ap(j) = dcmplx(0.0D0,0.0D0)
          B(j) = dcmplx(0.0D0,0.0D0)
          Bp(j) = dcmplx(0.0D0,0.0D0)
      enddo

c-----calulations of the max and minimal numbers of the harmonics
c      write(*,*)'DKtens before harmon_s nmax',nmax
      call harmon_z(X,Y,nll,vflow,c,beta,lambda,nmax,
     .n_harm1,n_harm2,n,istop)

c      write(*,*)'DKtens after harmon_s nmax',nmax
      if (istop.eq.0) then
         write(*,*)'DKtens istop=0 it is impossible to calulate' 
         write(*,*)'the harmonics of susceptibily tensor n_harm0 nmax='
     .   ,nmax
      endif
      
c   Now, calculate the modified bessel functions 
c   Remember that these are actually exp(-lambda) * In(lambda)

cyup      if (istop.eq.0) goto 10
      if(istop.eq.0)then
      if(n_harm1.lt.0 .and. n_harm2.lt.0)then
         n_harm1=-nmax
         n_harm2=n_harm1+2
         n=min(n,nmax)
         else
         n_harm1=nmax-2
         n_harm2=n_harm1+2
         n=min(n,nmax)
      endif
      endif
      
      call ibess0(lambda,ri01)
      call ibessn(lambda,n+2,ri01,xne) !in DKtens

cSAP090306
cSAP090204
c      n=n_new
      n_new=n
      do j=1,n
        if (dabs(xne(j)).lt.1.d-14) then
cSm060816           
          call harmon_z(X,Y,nll,vflow,c,beta,lambda,j+2,
     .    n_harm1,n_harm2,n_new,istop)

c          write(*,*)'DKtens after harmon_z nmax,n,n_new',nmax,n,n_new
c          goto 10
           goto 9
        endif
      enddo

cSm060816
 9    continue
      n=n_new

 10   continue
c      write(*,*)'DKtens j, xne(j)',j,xne(j)
c      write(*,*)'DKtens n_harm1,n_harm2', n_harm1,n_harm2

c  We need to take derivatives of D(w,k) with respect to
c  w,np,nll,  X, and  Y

c  Derivative wrt np

cSm990823
      betat_av= 0.5D0 * beta /  T_av     !d(beta_||)/d(T_av||)
      betatpop= -beta /(1.d0+2.d0*tpop)  !d(beta_||)/d(tpop||)

      do jj = 0,7 
         j = jj
      if (jj.eq.0) j = 1

      npp =  0.0D0
      nllp = 0.0D0
      Xp = 0.0D0
      Yp = 0.0D0
cSm990823    
c      betap= 0.50D0 * beta /te !  d(beta_||)/d(T_||)
      betap=0.d0
      tpopp=0.d0
      vflowp=0.d0  

      if (j.eq.1) npp    = 1.0D0
      if (j.eq.2) nllp   = 1.0D0
      if (j.eq.3)  Xp    = 1.0D0
      if (j.eq.4)  Yp    = 1.0D0
      if (j.eq.5) betap  = 1.0D0
      if (j.eq.6) tpopp  = 1.0D0
      if (j.eq.7) vflowp = 1.0D0

cSm990728
      lambdap = 2.0D0*lambda*(npp/np+betap/beta-Yp/Y+0.5d0*tpopp/tpop)
 
cFill the arrays for both positive and negative harmonics 

cwith A(1)     -> -n harmonic
c     A(n + 1) ->  0 harmonic
c     A(2n+ 1) -> +n harmonic
c     etc for B,Ibess...

c      do q = -n,n
cyup       if (istop.eq.0) goto 20 
       do q = n_harm1,n_harm2 

cSm990728
         cx = (1.d0 - nll*vflow/c - q *  Y) / nll / beta
cSm990728
         cxp = -q*Yp/nll/beta  
     .   -nllp * (1.d0-q*Y-nll*vflow/c) / nll**2 / beta 
     .   -nllp*vflow/c/nll/beta
     .   -betap* (1.d0-q*Y-nll*vflow/c) / nll    / beta**2
     .   -vflowp/beta

c       call CZETA(CX,CZ0,CZ1,IERR)           
        call CZETA0(nll,CX,CZ0,CZ1,IERR)      
cSm990728
        cz1_p=cz1
        cz1 = cz1*cxp

c   Use the Stix notation of A(n) and B(n)
c
c     My functions A(n) = omega * Astix(n)
c           and    B(n) = omega * Bstix(n) / c


         A(q+n+1)=(tpop-1.d0)+
     .   ((1.0D0-nll*vflow/c-q*Y)*tpop+q*Y)*cz0/nll/beta 

         B(q+n+1)  =  ((1.D0-nll*vflow/c)+ (1.D0-q*Y)*A(q+n+1)) /nll

         Ap(q+n+1)=((-1.0D0*nllp*vflow/c-q*Yp)*tpop+q*yp)*cz0/nll/beta 
     .+  ((1.D0-nll*vflow/c-q*Y)*tpop + q*Y)* cz1 /nll / beta  
     .-  ((1.D0-nll*vflow/c-q*Y)*tpop + q*Y)* cz0 *nllp /nll**2 /beta 
     .-  ((1.D0-nll*vflow/c-q*Y)*tpop + q*Y)* cz0 *betap /nll / beta**2
     .+  tpopp*(1.d0+(1.D0-nll*vflow/c-q*Y)*cz0/nll/beta)
     .+  vflowp*(-nll*tpop*cz0/nll/beta)

 

      Bp(q+n+1)=-nllp*
     . ((1.D0-q*Y)*A(q+n+1)+1.D0)/nll**2 
     . -Yp*q*A(q+n+1)/nll
     . -vflowp
     . +(1.D0-q*Y)*Ap(q+n+1)/nll 

c     Remember that Ibess = exp(-lambda)*I(n,lambda)

         if (q.eq.0)  Ibess(q+n+1) = ri01
         if (q.ne.0)  Ibess(q+n+1) = xne(Iabs(q))

         if (q.eq.0)  Ibessp(q+n+1) = xne(1) 
         if (q.ne.0)  Ibessp(q+n+1)=xne(1+iabs(q))
     .                 +Iabs(q)* xne(iabs(q)) / lambda

         if (q.eq.0)  Ibesspp(q+n+1) = xne(2)+xne(1)/lambda 
cSm990818
c         if (q.ne.0)  Ibesspp(q+n+1)=xne(2+iabs(q)) 
c     .     +(2.0D0*Iabs(q)+ 1)*xne(iabs(q)+1 )/lambda
c     .     +(Iabs(q)*Iabs(q)/lambda - Iabs(q)/lambda)*xne(iabs(q))
         if (q.ne.0)  Ibesspp(q+n+1)=xne(2+iabs(q)) 
     .     +(2.0D0*Iabs(q)+ 1)*xne(iabs(q)+1 )/lambda
     .     +(Iabs(q)*Iabs(q) - Iabs(q))/lambda**2*xne(iabs(q))

c     Ibp  =  l' * d/dl ( exp(-l) *  I(n,lambda) )
c     Ibpp  = l' * d/dl ( exp(-l) * I'(n,lambda) )

         Ibp(q+n+1) =  lambdap* (Ibessp(q+n+1)  - Ibess(q+n+1))
         Ibpp(q+n+1) = lambdap* (Ibesspp(q+n+1) - Ibessp(q+n+1))

      enddo
 20   continue
      if (jj.ne.0) then

c  use the hermitian part
      Kxx = 0.5D0 * ( K(1,1) + dconjg(K(1,1) )) !kxx
      Kyy = 0.5D0 * ( K(2,2) + dconjg(K(2,2) )) !kyy
      Kzz = 0.5D0 * ( K(3,3) + dconjg(K(3,3) )) !kxx
      Kxy = 0.5D0 * ( K(1,2) + dconjg(K(2,1) )) !kxy
      Kxz = 0.5D0 * ( K(1,3) + dconjg(K(3,1) )) !kxz
      Kyz = 0.5D0 * ( K(2,3) + dconjg(K(3,2) )) !kxz
      Kyx = 0.5D0 * ( K(2,1) + dconjg(K(1,2) )) !kyx
      Kzx = 0.5D0 * ( K(3,1) + dconjg(K(1,3) )) !kzx
      Kzy = 0.5D0 * ( K(3,2) + dconjg(K(2,3) )) !kxz
      endif 

      if (jj.eq.0) then
c     complete tensor for damping calculation
      Kxx = K(1,1)  !kxx
      Kyy = K(2,2)  !kyy
      Kzz = K(3,3)  !kxx
      Kxy = K(1,2)  !kxy
      Kxz = K(1,3)  !kxz
      Kyz = K(2,3)  !kxz
      Kyx = K(2,1)  !kyx
      Kzx = K(3,1)  !kzx
      Kzy = K(3,2)  !kxz
      end if

      Kxxp = (0.0D0,0.0D0)
      Kyyp = (0.0D0,0.0D0)
cSm990818
c      Kzzp = (0.0D0,0.0D0)
      Kzzp= dcmplx(-nllp*2.d0*X*vflow/c/(nll**2*beta**2*tpop)
     .             +Xp*2.d0*vflow/c/(nll*beta**2*tpop)
     .             -betap*4.d0*X*vflow/c/(nll*beta**3*tpop)
     .             +vflowp*2.d0*X/(nll*beta**2*tpop)
     .             -tpopp*X*2.d0*vflow/c/(nll*beta**2*tpop**2),0.d0)
      Kxyp = (0.0D0,0.0D0)
      Kxzp = (0.0D0,0.0D0)
      Kyzp = (0.0D0,0.0D0)

c  There is no shortcut: just do it!!!!
c  Derivative of all the other terms

c      do q = -n,n
cyup      if (istop.eq.0) goto 30
      do q  = n_harm1,n_harm2

      Kxxp = Kxxp +  Xp* q**2 * Ibess(q+n+1) * A(q+n+1) / lambda 
     . +  X* q**2 * Ibp(q+n+1) * A(q+n+1) / lambda 
     . +  X* q**2 * Ibess(q+n+1) * Ap(q+n+1) / lambda 
     . -  X* q**2 * Ibess(q+n+1) * A(q+n+1) * lambdap / lambda**2

      Kxyp= Kxyp-i*q* Xp*(Ibess(q+n+1)-Ibessp(q+n+1))*A(q+n+1)
     . - i*q* X*(Ibp(q+n+1)-Ibpp(q+n+1))*A(q+n+1)
     . - i*q* X*(Ibess(q+n+1)-Ibessp(q+n+1))*Ap(q+n+1)

      Kyyp=Kyyp+ Xp*
     . ( q**2*Ibess(q+n+1)/lambda+2.0D0*lambda*Ibess(q+n+1)
     .  -2.d0*lambda*Ibessp(q+n+1) )* A(q+n+1) +
     .  X*(q**2*Ibp(q+n+1)/lambda-lambdap*q**2 * Ibess(q+n+1)/lambda**2   
     .  + 2.0D0*lambdap*Ibess(q+n+1)  +  2.0D0*lambda*Ibp(q+n+1)
     .  -2.d0*lambdap*Ibessp(q+n+1) -2.d0*lambda*Ibpp(q+n+1))* A(q+n+1)
     .  + X *(q**2*Ibess(q+n+1)/lambda+2.0D0*lambda*Ibess(q+n+1)
     . -2.d0*lambda*Ibessp(q+n+1) )* Ap(q+n+1)

      Kxzp= Kxzp+ Xp*np*q*Ibess(q+n+1) * B(q+n+1) / lambda /  Y
     . + X*npp*q*Ibess(q+n+1) * B(q+n+1) / lambda /  Y
     . + X*np*q*Ibp(q+n+1) * B(q+n+1) / lambda /  Y
     . + X*np*q*Ibess(q+n+1) * Bp(q+n+1) / lambda /  Y
     . -lambdap*  X*np*q*Ibess(q+n+1) * B(q+n+1) / lambda**2 /  Y
     . - Yp *  X*np*q*Ibess(q+n+1) * B(q+n+1) / lambda /  Y**2

      Kyzp=Kyzp+i* Xp*np*(Ibess(q+n+1)-Ibessp(q+n+1))*B(q+n+1)/ Y
     . +i* X*npp*(Ibess(q+n+1)-Ibessp(q+n+1))*B(q+n+1)/ Y
     . +i* X*np*(Ibp(q+n+1)-Ibpp(q+n+1))*B(q+n+1)/ Y
     . +i* X*np*(Ibess(q+n+1)-Ibessp(q+n+1))*Bp(q+n+1)/ Y
     . - Yp*i* X*np*(Ibess(q+n+1)-Ibessp(q+n+1))*B(q+n+1)/ Y**2

      Kzzp= Kzzp
     .+2.0D0* Xp*(1.D0-q* Y)*Ibess(q+n+1)*B(q+n+1)/nll/beta**2
     .-2.0D0* X*q* Yp*Ibess(q+n+1)*B(q+n+1)/nll/beta**2
     .+2.0D0* X*(1.D0-q* Y)*Ibp(q+n+1)*B(q+n+1)/nll/beta**2
     .+2.0D0* X*(1.D0-q* Y)*Ibess(q+n+1)*Bp(q+n+1)/nll/beta**2
     .-nllp*2.0D0* X*(1.D0-q* Y)*Ibess(q+n+1)*
     . B(q+n+1)/nlls/beta**2
     .-2.0D0*betap*2.0* X*(1.D0-q* Y)*Ibess(q+n+1)*B(q+n+1)/nll
     . /beta**2

      enddo
 30   continue

      Kyxp = -1.0*Kxyp
      Kzxp = Kxzp
      Kzyp = -1.0*Kyzp

c     The conjugation and derivative operators should commute 
c     ? check CBF
c

      if (jj.ne.0) then
cSm990820
        dK(1,1,j) = kxxp 
        dK(2,2,j) = kyyp 
        dk(3,3,j) = kzzp 
        dk(1,2,j) = kxyp
        dk(1,3,j) = kxzp
        dk(2,1,j) = kyxp
        dk(2,3,j) = kyzp
        dk(3,1,j) = kzxp
        dk(3,2,j) = kzyp
        kxxp = 0.5D0 * ( dK(1,1,j) + dconjg(dK(1,1,j))) !kxx
        kyyp = 0.5D0 * ( dK(2,2,j) + dconjg(dK(2,2,j))) !kyy
        kzzp = 0.5D0 * ( dK(3,3,j) + dconjg(dK(3,3,j))) !kzz
        kxyp = 0.5D0 * ( dK(1,2,j) + dconjg(dK(2,1,j))) !kxy
        kxzp = 0.5D0 * ( dK(1,3,j) + dconjg(dK(3,1,j))) !kxz
        kyxp = 0.5D0 * ( dK(2,1,j) + dconjg(dK(1,2,j))) !kyx
        kyzp = 0.5D0 * ( dK(2,3,j) + dconjg(dK(3,2,j))) !kyz
        kzxp = 0.5D0 * ( dK(3,1,j) + dconjg(dK(1,3,j))) !kzx
        kzyp = 0.5D0 * ( dK(3,2,j) + dconjg(dK(2,3,j))) !kzy

c       kxxp = 0.5D0 * ( kxxp + dconjg(kxxp )) !kxx
c       kyyp = 0.5D0 * ( kyyp + dconjg(kyyp )) !kyy
c       kzzp = 0.5D0 * ( kzzp + dconjg(kzzp )) !kzz
c       kxyp = 0.5D0 * ( kxyp + dconjg(kyxp )) !kxy
c       kxzp = 0.5D0 * ( kxzp + dconjg(kzxp )) !kxz
c       kyxp = 0.5D0 * ( kyxp + dconjg(kxyp )) !kyx
c       kyzp = 0.5D0 * ( kyzp + dconjg(kzyp )) !kyz
c       kzxp = 0.5D0 * ( kzxp + dconjg(kxzp )) !kzx
c       kzyp = 0.5D0 * ( kzyp + dconjg(kyzp )) !kzy
      endif

      dK(1,1,j) = kxxp 
      dK(2,2,j) = kyyp 
      dK(3,3,j) = kzzp 
      dK(1,2,j) = kxyp
      dK(1,3,j) = kxzp
      dk(2,1,j) = kyxp
      dK(2,3,j) = kyzp
      dK(3,1,j) = kzxp
      dK(3,2,j) = kzyp
     
c   Now, compute the derivative of the D(w,k) wrt all 7

      dp = 
     . (Kxxp-2.0D0*nll*nllp) *   (Kyy-nlls-nps)   *  (Kzz-nps)
     .+(Kxx-nlls)*(Kyyp-2.0D0*nllp*nll-2.0D0*np*npp)*     (Kzz-nps)
     .+(Kxx-nlls)*   ( Kyy-nlls-nps)  *        ( Kzzp-2.0D0*np*npp)

      dp = dp 
     .+ Kxyp* Kyz *( Kzx  +np*nll) 
     .+ Kxy * Kyzp*( Kzx  +np*nll) 
     .+ Kxy * Kyz *( Kzxp +nll*npp + np*nllp) 

      dp = dp 
     .+ (Kxzp+npp*nll+np *nllp )* Kyx * Kzy
     .+ (Kxz +np*nll)             * Kyxp* Kzy
     .+ (Kxz +np*nll)             * Kyx * Kzyp

      dp = dp
     .-(Kzxp+ npp*nll+np*nllp)*( Kyy-nlls-nps)      *( Kxz+np*nll)
     .-(Kzx+np*nll) * (Kyyp-2.0D0*(np*npp+nll*nllp))*( Kxz+np*nll)
     .-(Kzx+np*nll)*( Kyy-nlls-nps)     *( Kxzp+nll*npp + np*nllp)

      dp = dp 
     .- Kzyp * Kyz  * ( Kxx - nlls )
     .- Kzy  * Kyzp * ( Kxx - nlls )
     .- Kzy  * Kyz  * ( Kxxp - 2.0D0*nll*nllp)

      dp = dp 
     .-(Kzzp - 2.0D0*np)* Kyx * Kxy
     .-(Kzz  - nps)     * Kyxp* Kxy
     .-(Kzz  - nps)     * Kyx * Kxyp

      if (jj.eq.0)  ddn = dp
   
      DD(j) = dp

c      write(*,*)'j,DD(j)',j,DD(j)

      end do
      

      do jj=1,3
        do j=1,3
          dKt_av=dK(jj,j,5)*betat_av               ! dK/dT_av
          dKtpop=dK(jj,j,5)*betatpop+dK(jj,j,6)    ! dK/dtpop
          dK(jj,j,5)=dKt_av
          dK(jj,j,6)=dKtpop
          dK(jj,j,7)= dK(jj,j,7)/c
        enddo
      enddo
      ddt_av=DD(5)*betat_av               ! dD/dT_av
      ddtpop=DD(5)*betatpop+DD(6)         ! dD/dtpop
      DD(5)=ddt_av                        ! [1/ev]
      DD(6)=ddtpop
      DD(7)=DD(7)/c                       ! [sec/cm]


      return
      end subroutine DKtens

      subroutine DKtens_s(mass,X,Y,T_av,tpop,vflow,nll_in,np_in,
     .dK)
c     Calculates the derivatives dK(3,3,8) of the hot non-relativietic
c     sucseptibilities for the given plasma specie.
c      
c     This function returns the derivatives dK(3,3,jj) of the Hermitian
c     part sucseptibilities with respect to 
c     np_in,nll_in,X,Y,Te_av,tpop,vflow (jj=1,7)
c
c     and the derivative dK(3,3,8) 
c     of the complete (Hermitian+anti-Hermitian )
c     sucseptibilies with respect np_in=N_perp (jj=8)
c
c     INPUTS:
c     mass the mass of the given specie (in the electron mass)
c     X = (fpe/f)**2    (for the given specie)
c     Y =  fce/f       the algebraic value  (for the given specie)
c     T_av=average temperature=Te(1+2*tpop)/3 in ev (for the given specie)
c          Te is the parallel temperature of the given specie in eV . 
c     tpop = T_perp/T_parallel (for the given specie)
c     vflow (future)    cm/sec (for the given specie)
c
c     nll_in - parallel index of refraction n.
c     np_in -  perpendicular index of refraction n.

c     OUTPUT:  
c
c     dK(*,*,1)  (the partial derivatives hermitian K_s wrt np
c     dK(*,*,2)  (the partial derivatives hermitian K_s wrt nll
c     dK(*,*,3)  (the partial derivatives hermitian K_s wrt  X)
c     dK(*,*,4)  (the partial derivatives hermitian K_s wrt  Y)
c     dK(*,*,5)  (the partial derivatives hermitian K_s wrt  Te_av) [1/eV]
c     dK(*,*,6)  (the partial derivatives hermitian K_s wrt  tpop)
c     dK(*,*,7)  (the partial derivatives hermitian K_s wrt  vflow)[sec/cm]
c
c     dK(*,*,8)  (the partial derivatives complete  K_s wrt np)   

      implicit none
 
      integer nmax
      parameter(nmax=300)

      double precision mass
      double precision np_in,nll_in,np,nll,npi
      double precision nps,nlls,c
      double precision X,Y,T_av,te,Xp,Yp,beta,betap,nllp,npp
      double precision pi,k4,vte,lambda,lambdap,tpop,tpopp,ri01,
     .vflow,vflowp,betat_av,betatpop
cSAP090204
c      double precision dr,di,xne(nmax),drp,Ibpp(nmax)
c      double precision Ibess(nmax),Ibessp(nmax),Ibesspp(nmax),Ibp(nmax)
c      double complex A(nmax),B(nmax),Ap(nmax),Bp(nmax)
      double precision dr,di,xne(nmax),drp,Ibpp(-(2*nmax+1):(2*nmax+1))
      double precision Ibess(-(2*nmax+1):(2*nmax+1)),
     &Ibessp(-(2*nmax+1):(2*nmax+1)),
     &Ibesspp(-(2*nmax+1):(2*nmax+1)),Ibp(-(2*nmax+1):(2*nmax+1))
      double complex A(-(2*nmax+1):(2*nmax+1)),
     &B(-(2*nmax+1):(2*nmax+1)),
     &Ap(-(2*nmax+1):(2*nmax+1)),Bp(-(2*nmax+1):(2*nmax+1))

      double complex i,dK(3,3,8)
      double complex dKt_av,dKtpop
      double complex kxxp,kxyp,kxzp,kyxp
      double complex kyyp,kyzp,kzxp,kzyp,kzzp
      double precision vflowdc
      double complex cx,cz0,cz1,cd,cxp,cz1_p
      integer*4 q,n,ierr,j,jj, qmin,qmax
cSm060816
      integer*4 n_new
      integer n_harm0,n_harm1,n_harm2 ! the limits for the harmonic numbers
      integer istop
        
      pi = 4.d0*datan(1.d0)
      c= 2.99792458d10                    !speed of light
      k4 = 4.19d7                   !constant for elec therm vel
      i = ( 0.0d0,1.0d0)            !imaginary number
      n = 20                        !max number of the cyclotron harmonics
      te=3.d0*T_av/(1.d0+2.d0*tpop) !longitudinal temperature
    
      vte =k4*dsqrt(2.0d0*te/mass)  !vte = sqrt(2.0*kTe/me/mass)
      beta = vte / c
      vflowdc=vflow/c

      call npnllmin(nll_in,np_in,nll,np) !in subroutine DKtens_s
     
      nps = np**2
      nlls= nll**2

c     lambda = 0.5 k_perp^2 * rho_larm^2

      lambda = nps * beta**2 * tpop / 2.0d0 /  Y**2
cSAP090204
      do j = 1,nmax
         xne(j) = 0.0D0
c         ibess(j) = 0.0D0
c         ibessp(j) = 0.0D0
c         ibp(j) = 0.0D0
c         ibpp(j) = 0.0D0
c         A(j) = dcmplx(0.0D0,0.0D0)
c         Ap(j) = dcmplx(0.0D0,0.0D0)
c         B(j) = dcmplx(0.0D0,0.0D0)
c         Bp(j) = dcmplx(0.0D0,0.0D0)
      enddo
      do j = -(2*nmax+1),(2*nmax+1)
c         xne(j) = 0.0D0
         ibess(j) = 0.0D0
         ibessp(j) = 0.0D0
         ibp(j) = 0.0D0
         ibpp(j) = 0.0D0
         A(j) = dcmplx(0.0D0,0.0D0)
         Ap(j) = dcmplx(0.0D0,0.0D0)
         B(j) = dcmplx(0.0D0,0.0D0)
         Bp(j) = dcmplx(0.0D0,0.0D0)
      enddo
c-----calulations of the max and minimal numbers of the harmonics
c      write(*,*)'dk_tens_s bef harmon_z nmax',nmax        
      call harmon_z(X,Y,nll,vflow,c,beta,lambda,nmax,
     .n_harm1,n_harm2,n,istop)
        
c      write(*,*)'dktens_s aft harmon_z n_harm1,n_harm2,n', 
c     &            n_harm1,n_harm2,n

      if (istop.eq.0) then
         write(*,*)'DKtens_s istop=0 it is impossible to calculate' 
         write(*,*)'the harmonics for susceptibily tensor n_harm0,nmax'
     .   ,nmax
         pause
      endif
    
c   Now, calculate the modified bessel functions 
c   Remember that these are actually exp(-lambda) * In(lambda)
      if (istop.eq.0) goto 10
      call ibess0(lambda,ri01)
      call ibessn(lambda,n+2,ri01,xne) !in DKtens_s

cSAP090306
cSAP090204
c      n=n_new 
       n_new=n
      do j=1,n

c         write(*,*)'j,xne(j)',j,xne(j)

        if (dabs(xne(j)).lt.1.d-14) then
c          write(*,*)'dktens_s nmax,j',j
cSm060816
          call harmon_z(X,Y,nll,vflow,c,beta,lambda,j+2,
     .    n_harm1,n_harm2,n_new,istop)


c          write(*,*)'dktens_s after harmon_z nmax,j',nmax,j
c          goto 10
          goto 9
        endif
      enddo

cSm060816
 9    continue
      n=n_new

 10   continue
c      write(*,*)'DKtens_s j, xne(j)',j,xne(j)
c      write(*,*)'DKtens_s n_harm1,n_harm2',n_harm1,n_harm2

      betat_av= 0.5D0 * beta /  T_av     !d(beta_||)/d(T_av||)
      betatpop= -beta /(1.d0+2.d0*tpop)  !d(beta_||)/d(tpop||)

      do jj = 0,7 
         j = jj
         if (jj.eq.0) j = 1

         npp =  0.0D0
         nllp = 0.0D0
         Xp = 0.0D0
         Yp = 0.0D0
         betap=0.d0
         tpopp=0.d0
         vflowp=0.d0  

         if (j.eq.1) npp    = 1.0D0
         if (j.eq.2) nllp   = 1.0D0
         if (j.eq.3)  Xp    = 1.0D0
         if (j.eq.4)  Yp    = 1.0D0
         if (j.eq.5) betap  = 1.0D0
         if (j.eq.6) tpopp  = 1.0D0
         if (j.eq.7) vflowp = 1.0D0

         lambdap = 2.0D0*lambda*
     .   (npp/np+betap/beta-Yp/Y+0.5d0*tpopp/tpop)
 
cFill the arrays for both positive and negative harmonics 

cwith A(1)     -> -n harmonic
c     A(n + 1) ->  0 harmonic
c     A(2n+ 1) -> +n harmonic
c     etc for B,Ibess...

c      do q = -n,n
      if (istop.eq.0) goto 20
      
      !!YuP[11-2016] Be sure to include q=0 harmonic (plus 2 extra):
      qmin=min(-2,n_harm1)
      qmax=max(n_harm2,2)
      n= max(iabs(qmin),iabs(qmax))

c      write(*,*)'DKtens_s n.n_harm1,n_harm2',n,n_harm1,n_harm2

      do q = qmin,qmax !YuP[11-2016] n_harm1,n_harm2

         cx = (1.d0 - nll*vflow/c - q *  Y) / nll / beta

         cxp = -q*Yp/nll/beta  
     .   -nllp * (1.d0-q*Y-nll*vflow/c) / nll**2 / beta 
     .   -nllp*vflow/c/nll/beta
     .   -betap* (1.d0-q*Y-nll*vflow/c) / nll    / beta**2
     .   -vflowp/beta
           
         call CZETA0(nll,CX,CZ0,CZ1,IERR)      

         cz1_p=cz1
         cz1 = cz1*cxp

c   Use the Stix notation of A(n) and B(n)
c
c     My functions A(n) = omega * Astix(n)
c           and    B(n) = omega * Bstix(n) / c

c         write(*,*)'DKtens_s,q,n,q+n+1',q,n,(q+n+1)

         A(q+n+1)=(tpop-1.d0)+
     .   ((1.0D0-nll*vflow/c-q*Y)*tpop+q*Y)*cz0/nll/beta 

         B(q+n+1)  =  ((1.D0-nll*vflow/c)+ (1.D0-q*Y)*A(q+n+1)) /nll

         Ap(q+n+1)=((-1.0D0*nllp*vflow/c-q*Yp)*tpop+q*yp)*cz0/nll/beta 
     .+  ((1.D0-nll*vflow/c-q*Y)*tpop + q*Y)* cz1 /nll / beta  
     .-  ((1.D0-nll*vflow/c-q*Y)*tpop + q*Y)* cz0 *nllp /nll**2 /beta 
     .-  ((1.D0-nll*vflow/c-q*Y)*tpop + q*Y)* cz0 *betap /nll / beta**2
     .+  tpopp*(1.d0+(1.D0-nll*vflow/c-q*Y)*cz0/nll/beta)
     .+  vflowp*(-nll*tpop*cz0/nll/beta)

 
      Bp(q+n+1)=-nllp*
     . ((1.D0-q*Y)*A(q+n+1)+1.D0)/nll**2 
     . -Yp*q*A(q+n+1)/nll
     . -vflowp
     . +(1.D0-q*Y)*Ap(q+n+1)/nll 

c     Remember that Ibess = exp(-lambda)*I(n,lambda)

         if (q.eq.0)  Ibess(q+n+1) = ri01
         if (q.ne.0)  Ibess(q+n+1) = xne(Iabs(q))

         if (q.eq.0)  Ibessp(q+n+1) = xne(1) 
         if (q.ne.0)  Ibessp(q+n+1)=xne(1+iabs(q))
     .                 +Iabs(q)* xne(iabs(q)) / lambda

         if (q.eq.0)  Ibesspp(q+n+1) = xne(2)+xne(1)/lambda 
         if (q.ne.0)  Ibesspp(q+n+1)=xne(2+iabs(q)) 
     .     +(2.0D0*Iabs(q)+ 1)*xne(iabs(q)+1 )/lambda
     .     +(Iabs(q)*Iabs(q) - Iabs(q))/lambda**2*xne(iabs(q))

c     Ibp  =  l' * d/dl ( exp(-l) *  I(n,lambda) )
c     Ibpp  = l' * d/dl ( exp(-l) * I'(n,lambda) )

         Ibp(q+n+1) =  lambdap* (Ibessp(q+n+1)  - Ibess(q+n+1))
         Ibpp(q+n+1) = lambdap* (Ibesspp(q+n+1) - Ibessp(q+n+1))

      enddo ! loop over q
 20   continue

      Kxxp = (0.0D0,0.0D0)
      Kyyp = (0.0D0,0.0D0)
      Kzzp= dcmplx(-nllp*2.d0*X*vflow/c/(nll**2*beta**2*tpop)
     .             +Xp*2.d0*vflow/c/(nll*beta**2*tpop)
     .             -betap*4.d0*X*vflow/c/(nll*beta**3*tpop)
     .             +vflowp*2.d0*X/(nll*beta**2*tpop)
     .             -tpopp*X*2.d0*vflow/c/(nll*beta**2*tpop**2),0.d0)
      Kxyp = (0.0D0,0.0D0)
      Kxzp = (0.0D0,0.0D0)
      Kyzp = (0.0D0,0.0D0)

c  There is no shortcut: just do it!!!!
c  Derivative of all the other terms

c      do q = -n,n
      if (istop.eq.0) goto 30
      do q = qmin,qmax !YuP[11-2016] n_harm1,n_harm2       

      Kxxp = Kxxp +  Xp* q**2 * Ibess(q+n+1) * A(q+n+1) / lambda 
     . +  X* q**2 * Ibp(q+n+1) * A(q+n+1) / lambda 
     . +  X* q**2 * Ibess(q+n+1) * Ap(q+n+1) / lambda 
     . -  X* q**2 * Ibess(q+n+1) * A(q+n+1) * lambdap / lambda**2

      Kxyp= Kxyp-i*q* Xp*(Ibess(q+n+1)-Ibessp(q+n+1))*A(q+n+1)
     . - i*q* X*(Ibp(q+n+1)-Ibpp(q+n+1))*A(q+n+1)
     . - i*q* X*(Ibess(q+n+1)-Ibessp(q+n+1))*Ap(q+n+1)

      Kyyp=Kyyp+ Xp*
     . ( q**2*Ibess(q+n+1)/lambda+2.0D0*lambda*Ibess(q+n+1)
     .  -2.d0*lambda*Ibessp(q+n+1) )* A(q+n+1) +
     .  X*(q**2*Ibp(q+n+1)/lambda-lambdap*q**2 * Ibess(q+n+1)/lambda**2   
     .  + 2.0D0*lambdap*Ibess(q+n+1)  +  2.0D0*lambda*Ibp(q+n+1)
     .  -2.d0*lambdap*Ibessp(q+n+1) -2.d0*lambda*Ibpp(q+n+1))* A(q+n+1)
     .  + X *(q**2*Ibess(q+n+1)/lambda+2.0D0*lambda*Ibess(q+n+1)
     . -2.d0*lambda*Ibessp(q+n+1) )* Ap(q+n+1)

      Kxzp= Kxzp+ Xp*np*q*Ibess(q+n+1) * B(q+n+1) / lambda /  Y
     . + X*npp*q*Ibess(q+n+1) * B(q+n+1) / lambda /  Y
     . + X*np*q*Ibp(q+n+1) * B(q+n+1) / lambda /  Y
     . + X*np*q*Ibess(q+n+1) * Bp(q+n+1) / lambda /  Y
     . -lambdap*  X*np*q*Ibess(q+n+1) * B(q+n+1) / lambda**2 /  Y
     . - Yp *  X*np*q*Ibess(q+n+1) * B(q+n+1) / lambda /  Y**2

      Kyzp=Kyzp+i* Xp*np*(Ibess(q+n+1)-Ibessp(q+n+1))*B(q+n+1)/ Y
     . +i* X*npp*(Ibess(q+n+1)-Ibessp(q+n+1))*B(q+n+1)/ Y
     . +i* X*np*(Ibp(q+n+1)-Ibpp(q+n+1))*B(q+n+1)/ Y
     . +i* X*np*(Ibess(q+n+1)-Ibessp(q+n+1))*Bp(q+n+1)/ Y
     . - Yp*i* X*np*(Ibess(q+n+1)-Ibessp(q+n+1))*B(q+n+1)/ Y**2

      Kzzp= Kzzp
     .+2.0D0* Xp*(1.D0-q* Y)*Ibess(q+n+1)*B(q+n+1)/nll/beta**2
     .-2.0D0* X*q* Yp*Ibess(q+n+1)*B(q+n+1)/nll/beta**2
     .+2.0D0* X*(1.D0-q* Y)*Ibp(q+n+1)*B(q+n+1)/nll/beta**2
     .+2.0D0* X*(1.D0-q* Y)*Ibess(q+n+1)*Bp(q+n+1)/nll/beta**2
     .-nllp*2.0D0* X*(1.D0-q* Y)*Ibess(q+n+1)*
     . B(q+n+1)/nlls/beta**2
     .-2.0D0*betap*2.0* X*(1.D0-q* Y)*Ibess(q+n+1)*B(q+n+1)/nll
     . /beta**2

      enddo !loop over q
 30   continue

      Kyxp = -1.0*Kxyp
      Kzxp = Kxzp
      Kzyp = -1.0*Kyzp

c     The conjugation and derivative operators should commute 
c     ? check CBF
c

      if (jj.ne.0) then
c-------Hermitian part
        dK(1,1,j) = kxxp 
        dK(2,2,j) = kyyp 
        dK(3,3,j) = kzzp 
        dK(1,2,j) = kxyp
        dK(1,3,j) = kxzp
        dK(2,1,j) = kyxp
        dK(2,3,j) = kyzp
        dK(3,1,j) = kzxp
        dK(3,2,j) = kzyp
        kxxp = 0.5D0 * ( dK(1,1,j) + dconjg(dK(1,1,j))) !kxx
        kyyp = 0.5D0 * ( dK(2,2,j) + dconjg(dK(2,2,j))) !kyy
        kzzp = 0.5D0 * ( dK(3,3,j) + dconjg(dK(3,3,j))) !kzz
        kxyp = 0.5D0 * ( dK(1,2,j) + dconjg(dK(2,1,j))) !kxy
        kxzp = 0.5D0 * ( dK(1,3,j) + dconjg(dK(3,1,j))) !kxz
        kyxp = 0.5D0 * ( dK(2,1,j) + dconjg(dK(1,2,j))) !kyx
        kyzp = 0.5D0 * ( dK(2,3,j) + dconjg(dK(3,2,j))) !kyz
        kzxp = 0.5D0 * ( dK(3,1,j) + dconjg(dK(1,3,j))) !kzx
        kzyp = 0.5D0 * ( dK(3,2,j) + dconjg(dK(2,3,j))) !kzy

        dK(1,1,j) = kxxp 
        dK(2,2,j) = kyyp 
        dK(3,3,j) = kzzp 
        dK(1,2,j) = kxyp
        dK(1,3,j) = kxzp
        dk(2,1,j) = kyxp
        dK(2,3,j) = kyzp
        dK(3,1,j) = kzxp
        dK(3,2,j) = kzyp
     
      else
c-------complete tensor for damping calculation
        dK(1,1,8) = kxxp 
        dK(2,2,8) = kyyp 
        dK(3,3,8) = kzzp 
        dK(1,2,8) = kxyp
        dK(1,3,8) = kxzp
        dk(2,1,8) = kyxp
        dK(2,3,8) = kyzp
        dK(3,1,8) = kzxp
        dK(3,2,8) = kzyp
      endif
      
      enddo ! loop over jj      

      do jj=1,3
        do j=1,3
          dKt_av=dK(jj,j,5)*betat_av               ! dK/dT_av
          dKtpop=dK(jj,j,5)*betatpop+dK(jj,j,6)    ! dK/dtpop
          dK(jj,j,5)=dKt_av
          dK(jj,j,6)=dKtpop
          dK(jj,j,7)= dK(jj,j,7)/c
        enddo
      enddo

      return
      end subroutine DKtens_s
      
      

      subroutine Ddhot(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
     .,nll_in,np_in,k_sum,dd,ddnp_h,ddnll_h,ddnp)
c
c     This function returns the derivatives of the Hermitian
c     part of the hot, nonrelativistic dispersion function
c     with respect to X,Y,Te_av,tpop,vflow,nll_in, and np_in
c
c     INPUTS:
c      nbulk the total number of plasma species
c      mass_ar - the masses  of the plasma species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=average temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar = T_perp/T_parallel  
c      vflow_ar    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c      k_sum(3,3) -  the nine components of the dielectric tensor
c                    evaluated at (mass,X,Y,Te,tpop,vflow,nll,np)
c     OUTPUT:
c      DD(1-5,nbulk) derivatives from the dispersion function (for hot
c                  Hermitian dielectric tensor ) with respect
c                  X_js,Y_js,T_av_js,tpop_js,vflow_js (js=1,nbulk)
c      DD(1,js)=dD/dX_js
c      DD(2,js)=dD/dY_js
c      DD(3,js)=dD/dT_av_js
c      DD(4,js)=dD/dtpop_js
c      DD(5,js)=dD/dvflow_js
c
c      derivatives from the dispersion function (for hot
c      Hermitian dielectric tensor ) with respect
c      ddnp_h  = dD/DN_perpendicular
c      ddnll_h = dD/dN_parallel 
c
c      derivative from the dispersion function (for hot
c      complete dielectric tensor ) with respect
c      ddnp = dD/DN_perpendicular
c-----------------------------------------------------------
      implicit none
c-----input
      integer nbulk
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double complex k_sum(3,3)
      double precision np_in,nll_in
c-----external DKtens_s
c-----output
c     derivatives from the hermitian dispersion function D
c     dd(1,js)=dD/dX_js,dd(2,js)=dD/dY_js,dd(3,js)=dD/dT_av_js
c     dd(4,js)=dD/dtpop_js,dd(5,js)=dD/dVflow_js
c     ddnp_h=dD/N_perp,ddnll_h=dD/dN_parallel
c     derivartive from the complete dispersion function 
c     ddnp=dD/N_perp
      double complex dd(5,*)  
      double complex ddnll_h,ddnp,ddnp_h 
c-----locals
      double precision np,nll,nlls,nps
      double complex dk(3,3,8),k_herm(3,3)
      double complex ddeps(3,3),ddeps_h(3,3)
      double complex ddnll
      double complex dp,dp1,dp2,dp8
      double complex ddp(8)
      integer js,j,i1,i2

!YuP      call npnllmin(nll_in,np_in,nll,np) ! YuP[2020-08-24]] Having nll=0 is ok here
      nll=nll_in !YuP[2020-08-24] Having nll=0 is ok here 
      np=max(0.d0,np_in) !YuP[2020-08-24] Having np=0 is ok here
      nlls=nll*nll
      nps=np*np

c-----Now, compute the derivatives of the D(w,k)  
c     for Hermitian tensor

c     compute k_herm the hermitian part of the dielectric tensor k_sum

      call herm(k_sum,k_herm)

c     compute the derivatives from the dispersion function D
c     for the hermitian hot tensor k_herm
c     with respect dielectric tensor components: ddeps_h
c     with respect N_perpendicular ,N_parallel : ddnp_h,ddnll_h)

      call dddeps(k_herm,nll,np,ddeps_h,ddnll_h,ddnp_h) !ok to have nll=0
     
c     compute the derivatives from the dispersion function D
c     for the complete hot tensor k_sum
c     with respect dielectric tensor components: ddeps
c     with respect N_perpendicular ,N_parallel : ddnp,ddnll)

      call dddeps(k_sum,nll,np,ddeps,ddnll,ddnp) !ok to have nll=0

      dp1=dcmplx(0.d0,0.d0) ! initialization
      dp2=dcmplx(0.d0,0.d0) ! initialization  
      dp8=dcmplx(0.d0,0.d0) ! initialization

      do js=1,nbulk ! the loop over the plasma species
c--------calculations of the derivatives dK(3,3,8) from the
c        sucseptibilities K_js for js specie (i1=1,2,3, i2=1,2,3):
c        dK(i1,i2,1)_js=dK(i1,i2)_js/dN_perp
c        dK(i3,i2,2)_js=dK(i1,i2)_js/dN_parallel
c        dK(i1,i2,3)_js=dK(i1,i2)_js/dX_js
c        dK(i1,i2,4)_js=dK(i1,i2)_js/dY_js
c        dK(i3,i2,5)_js=dK(i1,i2)_js/dT_av_js
c        dK(i1,i2,6)_js=dK(i1,i2)_js/dtpop_js
c        dK(i1,i2,7)_js=dK(i1,i2)_js/dV_js
c        dK(i1,i2,8)_js=dK(i1,i2)_js/dN_per
c        dK(i1,i2,3,jj=1-7) for Hermitian part of the tensor
c        dK(i1,i2.jj=8)=dD/dN_perp for the complete tensor

         call DKtens_s(mass_ar(js),x_ar(js),y_ar(js),t_av_ar(js),
     .   tpop_ar(js),vflow_ar(js),nll_in,np_in,dk)

         do j = 1,8
c----------Now, compute the derivative of the D(w,k) with respect
c          all 8 variables  
c          for Hermitian and complete (j=8) tensor
c          Sum_j{i1=1,3, i2=1,3}[dD/deps(i1,i2)*deps_js(i1,i2)/d(variable_j)]
c          variable_j=(j=1 n_perp),(j=2 n_par),(j=3 X_js),(j=4 Y_js),
c          (j=5 T_av_js), (j=6 t_pop_js) ,(j=7 v_js), (j=8 n_perp)
c          the first 7 derivatives (j=1,..,7) are calculated from Hermitian D,
c          the last 8th derivative is calculated from the complete D

           dp=dcmplx(0.d0,0.d0) ! initialization
                    
	   do i1=1,3
	     do i2=1,3
               
	       if (j.ne.8) then
	         dp=dp+ddeps_h(i1,i2)*dk(i1,i2,j)
               else
	         dp=dp+ddeps(i1,i2)*dk(i1,i2,8)
	       endif
               
	     enddo
	   enddo

           if ((j.ge.3).and.(j.le.7)) dd(j-2,js) = dp

c          summation by all js=1,...,nbulk species 
           if  (j.eq.1) dp1=dp1+dp 
           if  (j.eq.2) dp2=dp2+dp 
           if  (j.eq.8) dp8=dp8+dp 

         enddo ! j
        
      enddo  ! js

      ddnp_h  = ddnp_h  + dp1    ! DD_hermitian/DN_perp
      ddnll_h = ddnll_h + dp2    ! DD_hermitian/DN_parallel
      ddnp    = ddnp    + dp8    ! DD_complete/DN_perp
                    
cyup 200  continue ! not used


      return
      end subroutine Ddhot




      subroutine dddeps(eps,nll,np,ddeps,ddnll,ddnp)
c     calculates the derivatives from the determinant 
c     D(eps,N_parrallel,N_perpendicular)) 
c     with respect eps(i,j) - ddeps(i,j)
c     with respect N_perp - ddnll and 
c     with respect N_perp - ddnp and 
c
c     input
c       eps - double complex dielectric tensor 
c       nll - N_parallel
c       np  - N_perpendicular
c     output
c       ddeps(i,j) double complex dD/deps(i,j)
c       ddnll      double complex dD/dN_parallel
c       ddnp       double complex dD/dN_parallel
c-----------------------------------------------

      implicit none
c-----input
      double precision nll,np
      double complex eps(3,3)
c-----output
      double complex ddeps(3,3),ddnll,ddnp
c-----locals
      double complex Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz
      double precision nps,nlls

      nps=np*np
      nlls=nll*nll

      Kxx = eps(1,1)
      Kyy = eps(2,2)
      Kzz = eps(3,3)
      Kxy = eps(1,2)
      Kxz = eps(1,3)
      Kyx = eps(2,1)
      Kyz = eps(2,3)
      Kzx = eps(3,1)
      Kzy = eps(3,2)

     
      ddeps(1,1)=(Kyy-nlls-nps) * (Kzz-nps) 
     .-  Kzy * Kyz

      ddeps(1,2)= Kyz * (Kzx+np*nll) 
     .- (Kzz - nps) * Kyx 

      ddeps(1,3)= Kyx * Kzy
     .- (Kzx+np*nll) * (Kyy-nlls-nps)

      ddeps(2,1)= (Kxz+np*nll) * Kzy
     .- (Kzz - nps) * Kxy

      ddeps(2,2)=(Kxx-nlls) * (Kzz-nps) 
     .- (Kzx+np*nll) * (Kxz+np*nll)

      ddeps(2,3)=Kxy * (Kzx+np*nll) 
     .-  Kzy * (Kxx-nlls)

      ddeps(3,1)=Kxy * Kyz  
     .- (Kyy-nlls-nps) * (Kxz+np*nll)

      ddeps(3,2)=(Kxz+np*nll) * Kyx 
     .-  Kyz * (Kxx-nlls)

      ddeps(3,3)=(Kxx-nlls) * (Kyy-nlls-nps)
     .-  Kyx * Kxy

      ddnll    =(-2.d0*nll) * (Kyy-nlls-nps) * (Kzz-nps)
     .+         (Kxx-nlls) * (-2.d0*nll) * (Kzz-nps)

     .+  Kxy * Kyz * np 

     .+ np * Kyx * Kzy

     .- np * (Kyy-nlls-nps) * (Kxz+np*nll)
     .- (Kzx+np*nll) * (-2.d0*nll) * (Kxz+np*nll)
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * np

     .-  Kzy * Kyz * (-2.d0*nll)

      ddnp     =(Kxx-nlls) * (-2.d0*np) * (Kzz-nps)
     .+         (Kxx-nlls) * (Kyy-nlls-nps) * (-2.d0*np )

     .+  Kxy * Kyz * nll 

     .+ nll * Kyx * Kzy

     .- (nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
     .- (Kzx+np*nll) * (-2.d0*np) * (Kxz+np*nll)
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (nll)

     .- (- 2.d0*np) * Kyx * Kxy

      return
      end subroutine dddeps


      double precision function ddwrap(log_nper)
c-----calculates Real(hot dispersion function)
c     INPUT:
c       log_nper=ln(Re(N_perpendicular))
c       the input data from common /nperpcom/
      implicit none
    
c-----input
      double precision log_nper
c-----local
      double precision nperp
      double complex  K(3,3),d
      include 'param.i'
      include 'nperpcom.i'
c-----external: dhot_sum
      double complex dhot_sum

      nperp = dexp(log_nper)
      d=dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     .vflowc,nllc,nperp,1,K)
      ddwrap = dreal(d) 

      return
      end function ddwrap


      double complex function  
     . dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     . vflow_ar,nll_in,np_in,iherm,K_sum)
c-----------------------------------------------------------
c     calculates the complex hot plasma dispersion function and
c     the hot dielectric tensor  K_sum(3,3)
c     for the hot electron and ion plasma
c     T.Stix, Waves in plasmas,(1992), p.252 (40), p.258 (57)
c     K_sum=1+sum{s}K_s  
c-----------------------------------------------------------
c     INPUTS:
c      nbulk - the total number of the plasma species 
c      mass_ar - the masses of the plasma  species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=average temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar = T_perp/T_parallel  
c      vflow_ar    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c      iherm =1 hermitian dielectric tensor, 2-full
c     OUTPUT:
c      K_sum(3,3) -  the nine components of the dielectric tensor
c                    evaluated at (mass,X,Y,Te,tpop,vflow,Nll,np)
c      dhotsum    -  hot double complex dispersion function
c-----------------------------------------------------------
      implicit none
c     input
      integer nbulk,iherm   
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double precision np_in,nll_in, np_min
c-----output
      double complex K_sum(3,3)
c-----local
      double precision np,nll,nps,nlls
      double complex K(3,3) ! sucseptibilities 
      double complex Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz

      integer j,j1,j2
c-----external DHOT_s,zeroK

      np_min=1.d-7

c-----initialization of the dielectric tensor K_sum(3,3)
      call zeroK(K_sum)
  
      do j=1,nbulk            
        call DHOT_s(mass_ar(j),x_ar(j),y_ar(j),t_av_ar(j),tpop_ar(j),
     .  vflow_ar(j),nll_in,np_in,iherm,K)    
         
        do j1=1,3
          do j2=1,3
            K_sum(j1,j2)=K_sum(j1,j2)+K(j1,j2)
          enddo
         enddo
       
      enddo

      K_sum(1,1) = K_sum(1,1) + dcmplx(1.0D0,0.0D0)
      K_sum(2,2) = K_sum(2,2) + dcmplx(1.0D0,0.0D0)
      K_sum(3,3) = K_sum(3,3) + dcmplx(1.0D0,0.0D0)
      
      Kxx = K_sum(1,1)
      Kyy = K_sum(2,2)
      Kzz = K_sum(3,3)
      Kxy = K_sum(1,2)
      Kxz = K_sum(1,3)
      Kyx = K_sum(2,1)
      Kyz = K_sum(2,3)
      Kzx = K_sum(3,1)
      Kzy = K_sum(3,2)

c-----hot dispersion function
!YuP      call npnllmin(nll_in,np_in,nll,np) !in func.dhot_sum: no need to adjust here
      nll=nll_in ! YuP 120430 ! Having nll=0 is ok here
      np=max(0.d0,np_in)   ! YuP 120430
      !np=max(np,np_min) !YuP[2020-08-31] Optionally: Set a lower limit for np
      !                 ! (as done in subr.npnllmin)
      
      nlls=nll**2
      nps=np**2

      dhot_sum= (Kxx-nlls) * (Kyy-nlls-nps) * (Kzz-nps)
     .+  Kxy * Kyz * (Kzx+np*nll) 
     .+ (Kxz+np*nll) * Kyx * Kzy
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
     .-  Kzy * Kyz * (Kxx-nlls)
     .- (Kzz - nps) * Kyx * Kxy

      !write(*,*)'in dhot_sum nlls,nps',nlls,nps
      !write(*,*)'Kxx,Kxy,Kxz',Kxx,Kxy,Kxz 
      !write(*,*)'Kyx,Kyy,Kyz',Kyx,Kyy,Kyz
      !write(*,*)'Kzx,Kzy,Kzz',Kzx,Kzy,Kzz 
      !write(*,*)'dhot_sum=',dhot_sum

      return
      end function dhot_sum

      subroutine dnd(z,r,phi,cnz,cnr,cm,
     .dnpdz,dnpdr,dnpdphi,dnpdcnz,dnpdcnr,dnpdcm,
     .dnlldz,dnlldr,dnlldphi,dnlldcnz,dnlldcnr,dnlldcm,
     .dnpdw,dnlldw)
c----------------------------------------------------------
c     calculates the derivatives from N_perp=np and N_parallel=nll
c     on z,r,phi,cnz=N_z,cnr=N_r,cm=N_phi/r
c
c     it uses bz,br,bphi from common one.i
c
c     b(z,r,phi) should be run before this subroutine.
c     b calculates magnetic field bz,br,bphi and derivatives from b^ by z,r,phi 
c     b put these data to common one.i
c--------------------------------------------------------- 
      implicit none
      include 'param.i'
      include 'one.i'
c     input
      double precision z,r,phi,cnz,cnr,cm
c     output
      double precision cnp,nll,
     .dnpdr,dnpdz,dnpdphi,dnpdcnz,dnpdcnr,dnpdcm,       ! derivatives from N_perpendicular
     .dnlldr,dnlldz,dnlldphi,dnlldcnz,dnlldcnr,dnlldcm, !derivatives from N_parallel
     &dnpdw,dnlldw
c       
c     locals
    
      nll=(cnz*bz+cnr*br+cm*bphi/r)/bmod
      cnp=dsqrt(cnz**2+cnr**2+(cm/r)**2-nll**2)
      nll=(cnz*bz+cnr*br+cm*bphi/r)/bmod
 
c-----derivatives from N_parallel    
      dnlldz=(cnz*dbzdz+cnr*dbrdz+cm*dbphdz/r)/bmod-nll/bmod*dbmdz     
      dnlldr=(cnz*dbzdr+cnr*dbrdr+cm*(dbphdr/r-bphi/r**2))/bmod-
     .       nll/bmod*dbmdr     
      dnlldphi=(cnz*dbzdph+cnr*dbrdph+cm*dbpdph/r)/bmod-nll/bmod*dbmdph
      !YuP: careful! There is {dbphdz,dbphdr,dbpdph}. There is no dbphdph !
      !Tested: set dnlldphi to 0 (tor.symmetry) --> same results.
      
      dnlldcnz=bz/bmod
      dnlldcnr=br/bmod
      dnlldcm=bphi/(r*bmod)
      dnlldw=-nll/frqncy  !dN_parallel/d_omega

c-----derivatives from N_perpendicular  
      dnpdz=-nll*dnlldz/cnp
      dnpdr=-(cm**2/r**3+nll*dnlldr)/cnp
      dnpdphi=-nll*dnlldphi/cnp
      dnpdcnz=(cnz-nll*dnlldcnz)/cnp
      dnpdcnr=(cnr-nll*dnlldcnr)/cnp
      dnpdcm=(cm/r**2-nll*dnlldcm)/cnp
      dnpdw=-cnp/frqncy   !dN_perpendicular/d_omega
c	write(*,*)'dnd dnlldcnz,dnpdcnz',dnlldcnz,dnpdcnz
      return
      end subroutine dnd




      subroutine herm(k_total,k_herm)
c     calculates the hermitian part of the complex matrix
c     k_total(3,3) 
c     input
c        k_total(3,3) double complex matrix
c      output
c         k_herm(3,3) hermitian part of k_total
c---------------------------------------------	
      implicit none
c     input
	double complex k_total(3,3)
c     output
      double complex k_herm(3,3)
c     locals
      integer i1,i2

	do i1=1,3
	  do i2=1,3
	    k_herm(i1,i2)=0.5d0*(k_total(i1,i2)+dconjg(k_total(i2,i1)))
	  enddo
      enddo

      return
      end subroutine herm



      subroutine aherm(k_total,k_aherm)
c     calculates anti-hermitian part of the complex matrix
c     k_total(3,3) 
c     input
c        k_total(3,3) double complex matrix
c      output
c         k_aherm(3,3) anti hermitian part of k_total
c---------------------------------------------	
      implicit none
c     input
      double complex k_total(3,3)
c     output
      double complex k_aherm(3,3)
c     locals
      integer i1,i2
      
      do i1=1,3
	 do i2=1,3
	   k_aherm(i1,i2)=0.5d0*(k_total(i1,i2)-dconjg(k_total(i2,i1)))
	 enddo
      enddo

      return
      end subroutine aherm





c    991223
      subroutine hotdervs(nbulk,u,wf,dddcnz,dddcnr,dddcm,
     .dddz,dddr,dddph,dddw)
c---------------------------------------------------------
c     analytical calculation of the derivatives for the ray-tracing
c     equations from the non-relativistic hot plasma
c     with ions and electrons
c--------------------------------------------------------- 
      implicit none
cSm030226
      include 'param.i'
c-----input
      double precision u(6)  ! z,r,phi,Nz,Nr,M
      double precision wf    ! wave frequency
      integer nbulk          ! the total number of the plasma components
c-----output
      double precision dddcnz,dddcnr,dddcm,
     .dddz,dddr,dddph,dddw   !derivatives from D
c-----uses
      external dhot_sum,Ddhot,b,x,y,rhof,tempe,
     .dxdz,dxdr,dxdphi,dydz,dydr,dydphi,
     .temperho,tpoprho,vflowrho,
     .dtempdz,dtempdr,dtpopdz,dtpopdr,dvflowdz,dvflowdr
      double complex dhot_sum
      double precision b,x,y,rhof,tempe,
     .temperho,tpoprho,vflowrho,
     .dxdz,dxdr,dxdphi,dydz,dydr,dydphi,
     .dtempdz,dtempdr,dtpopdz,dtpopdr,dvflowdz,dvflowdr

c-----locals

c      integer nbulka
cSm030226
c      parameter (nbulka=5)
c      parameter (nbulka=nbulka)
      double complex dd(5,nbulka) !dD/d(X_s,Y_s,Tav_s,tpop_as,V_s)
cSm030226
c      double precision t_kev,nll,np
      double precision t_kev,nll,nperp

      double complex K_sum(3,3),dK(3,3,7),ddnp_h,ddnll_h,ddnp,d
      double precision z,r,phi,cnz,cnr,cm,bmod,rholoc
      integer i
      double precision dxdze,dxdre,dxdphie,dydze,dydre,dydphie, 
     .dtempdze,dtempdre,dtpopdze,dtpopdre,dvflowze,dvflowre,
     .dtempdpe,dtpopdpe,dvflowpe,
     .dxdwe,dydwe,
     .dnpdz,dnpdr,dnpdphi,dnpdcnz,dnpdcnr,dnpdcm,
     .dnlldz,dnlldr,dnlldphi,dnlldcnz,dnlldcnr,dnlldcm,
     .dnpdw,dnlldw

cSm030225
c      integer ncomp1
c      parameter (ncomp1=5) 
c      double precision mass_ar(ncomp1),x_ar(ncomp1),y_ar(ncomp1),
c     .t_av_ar(ncomp1),tpop_ar(ncomp1),vflow_ar(ncomp1)

      double precision mass_ar(nbulka),x_ar(nbulka),y_ar(nbulka),
     .t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka)


cfor test only
      double precision step,zp,zm,rholocp,xep,yep,te_keVp,tep,
     .tpopep,vflowep,rholocm,xem,yem,te_keVm,tem,
     .tpopem,vflowem,rp,rm,phip,phim,cnzp,cnzm,cnrp,cnrm,cmp,cmm,
     .nllp,nllm,cnpp,cnpm
      double complex dp,dm,ddnum(7),deld
      double precision rdp,rdm,rdeld,rddnum(7)
      integer j
      double precision ddtest(7),ktest(3,3),ddntest

ctest 
      integer k      

      step=1.d-5
      
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)

      bmod=b(z,r,phi)
      rholoc=rhof(z,r,phi)

c      write(*,*)'hotdervs z,r,phi,bmod',z,r,phi,bmod
c      write(*,*)'cnz,cnr,cm',cnz,cnr,cm

cSm030226
c     if(nbulk.gt.ncomp1) then
c          write(*,*)'in forest.f in hotdervs nbulk.gt.ncomp1'
c          write(*,*)'put the value of parameter ncomp1.ge.nbulk'
c          write(*,*)'parameter ncomp is given hotdervs'
c          write(*,*)'ncomp1,nbulk',ncomp1,nbulk
c          stop
c     endif
      if(nbulk.gt.nbulka) then
          write(*,*)'in forest.f in hotdervs nbulk.gt.nbulka'
          write(*,*)'put the value of parameter nbulka.ge.nbulk'
          write(*,*)'parameter nbulka is given in param.i in hotdervs'
          write(*,*)'nbulka,nbulk',nbulka,nbulk
          stop
      endif

c-----The initilization mass_ar
      call put_mass(mass_ar,nbulk)
      
c-----nll=(cnz*bz+cnr*br+cm*bphi/r)/bmod
c     nperp=dsqrt(cnz**2+cnr**2+(cm/r)**2-nll**2)
c     calculation nll=N_parallel and nperp=N_perp
c     for given z,r,phi,cnz,cnr,cm,
      call nllcnp_s(z,r,phi,cnz,cnr,cm,nll,nperp)

ctest numer deriv from dhot
c      do i=1,nbulk                  

c          x_ar(i)=x(z,r,phi,i)
c          y_ar(i)=y(z,r,phi,i)
c          if(i.eq.1) y_ar(1)=-y_ar(1)
c          t_keV=tempe(z,r,phi,i)   !keV averaged temperature
c          t_av_ar(i)=t_keV*1.d3    !eV
c          tpop_ar(i)=tpoprho(rholoc,i)
c	  vflow_ar(i)=vflowrho(rholoc,i) !cm/sec
	
c     enddo
      
c      cnpp=nperp+step
c      cnpm=nperp-step

c     write(*,*)'nperp,step,cnpp,cnpm',nperp,step,cnpp,cnpm
     
c      dp= dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
c     .vflow_ar,nll,cnpp,1,K_sum)

c     write(*,*)'hotdervs before dhot_sum dm'

c      dm= dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
c     .vflow_ar,nll,cnpm,1,K_sum)

c     write(*,*)'hotdervs after dhot_sum dm'

c      ddnum(1)=(dp-dm)/(2.d0*step)
c      deld=dp-dm
c      rdp=dreal(dp)
c      rdm=dreal(dm)
c      rdeld=rdp-rdm
c      rddnum(1)=(rdp-rdm)/(2.d0*step) 

c      write(*,*)'num dhotdnp '
c      write(*,*)'cnpp,cnpm',cnpp,cnpm
c      write(*,*)'dp,dm',dp,dm
c     write(*,*)'deld',deld
c     write(*,*)'rdp,rdm,rdeld',rdp,rdm,rdeld
c      write(*,*)'ddnum(1)',ddnum(1)
c      write(*,*)'rddnum(1)',rddnum(1)
 
c      nllp=nll+step
c      nllm=nll-step

c      write(*,*)'before nllp dhot '

c      dp= dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
c     .vflow_ar,nllp,nperp,1,K_sum)
c      dm= dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
c     .vflow_ar,nllm,nperp,1,K_sum)

c      ddnum(2)=(dp-dm)/(2.d0*step)
c      deld=dp-dm
c      rdp=dreal(dp)
c      rdm=dreal(dm)
c      rdeld=rdp-rdm
c      rddnum(2)=(rdp-rdm)/(2.d0*step) 

c      write(*,*)'num dhotdnll '
c	write(*,*)'nllp,nllm',nllp,nllm
c 	write(*,*)'dp,dm',dp,dm
c 	write(*,*)'deld',deld
c 	write(*,*)'rdp,rdm,rdeld',rdp,rdm,rdeld
c      write(*,*)'ddnum(2)',ddnum(2)
c      write(*,*)'rddnum(2)',rddnum(2)
 
ctestend

c        write(*,*)'hotdervs nbulk',nbulk
	do i=1,nbulk                  

          x_ar(i)=x(z,r,phi,i)
          y_ar(i)=y(z,r,phi,i)
cSm000324
          if(i.eq.1) y_ar(1)=-y_ar(1)
          t_keV=tempe(z,r,phi,i)   !keV averaged temperature
          t_av_ar(i)=t_keV*1.d3    !eV
          tpop_ar(i)=tpoprho(rholoc,i)
	  vflow_ar(i)=vflowrho(rholoc,i) !cm/sec
	
	enddo

c	do i=1,nbulk                  
c         write(*,*)'i,x_ar(i),y_ar(i),t_av_ar(i)',
c     &   i,x_ar(i),y_ar(i),t_av_ar(i)
c         write(*,*)'mass_ar(i),tpop_ar(i),vflow_ar(i)',
c     &   mass_ar(i),tpop_ar(i),vflow_ar(i)	
c	enddo
c        write(*,*)'nll,nperp',nll,nperp          
        
	d= dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .vflow_ar,nll,nperp,1,K_sum)
c        write(*,*)'d',d
c        write(*,*)'K_sum',K_sum

      call Ddhot(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
     .,nll,nperp,K_sum,dd,ddnp_h,ddnll_h,ddnp)

c       write(*,*)'ddnp_h,ddnll_h,ddnp',ddnp_h,ddnll_h,ddnp
  
c       do i=1,nbulk
c        write(*,*)'i,nbulk,nbulka',i,nbulk,nbulka
c        do k=1,5
c         write(*,*)'k,dd(k,i)',k,dd(k,i)
c        enddo
c       enddo 

      call dnd(z,r,phi,cnz,cnr,cm,
     .dnpdz,dnpdr,dnpdphi,dnpdcnz,dnpdcnr,dnpdcm,
     .dnlldz,dnlldr,dnlldphi,dnlldcnz,dnlldcnr,dnlldcm,
     .dnpdw,dnlldw)

      dddz = ddnp_h*dnpdz+ddnll_h*dnlldz
      dddr = ddnp_h*dnpdr+ddnll_h*dnlldr
      dddph= ddnp_h*dnpdphi+ddnll_h*dnlldphi
      dddcnz=ddnp_h*dnpdcnz+ddnll_h*dnlldcnz
      dddcnr=ddnp_h*dnpdcnr+ddnll_h*dnlldcnr
      dddcm =ddnp_h*dnpdcm +ddnll_h*dnlldcm
      
      dddw=ddnp_h*dnpdw+ddnll_h*dnlldw

c      write(*,*)'dddz,dddr,dddph',dddz,dddr,dddph
c      write(*,*)'dddcnz,dddcnr,dddcm',dddcnz,dddcnr,dddcm
c      write(*,*)'dddw',dddw
c      write(*,*)'ddnp_h,ddnll_h',ddnp_h,ddnll_h
c      write(*,*)'dnpdz,dnpdr,dnpdphi',dnpdz,dnpdr,dnpdphi
c      write(*,*)'dnlldz,dnlldr,dnlldphi',dnlldz,dnlldr,dnlldphi
c      write(*,*)'dnpdcnz,dnpdcnr,dnpdcm',dnpdcnz,dnpdcnr,dnpdcm
c      write(*,*)'dnlldcnz,dnlldcnr,dnlldcm',dnlldcnz,dnlldcnr,dnlldcm




      do i=1,nbulk

        dxdze=dxdz(z,r,phi,i)
        dxdre=dxdr(z,r,phi,i)
        dxdphie=dxdphi(z,r,phi,i)

        dydze=dydz(z,r,phi,i)
        dydre=dydr(z,r,phi,i)
        dydphie=dydphi(z,r,phi,i)
        if(i.eq.1) then ! negative Y for electrons
           dydze=-dydze
           dydre=-dydre
           dydphie=-dydphie
        endif 

        dtempdze=dtempdz(z,r,phi,i)
        dtempdre=dtempdr(z,r,phi,i)
        dtempdpe=0.d0			      !d(T_average)/d(phi)

        dtpopdze=dtpopdz(z,r,phi,i)
        dtpopdre=dtpopdr(z,r,phi,i)
        dtpopdpe=0.d0				  !d(tpop)/d(phi)

        dvflowze=dvflowdz(z,r,phi,i)
        dvflowre=dvflowdr(z,r,phi,i)
        dvflowpe=0.d0				  !d(vflow)\d(phi)

        dddz = dddz+
cSm021114
c     .       dd(1,i)*dxdze+dd(2,i)*dydze+dd(3,i)*dtempdze+
     .       dd(1,i)*dxdze+dd(2,i)*dydze+dd(3,i)*dtempdze*1.d3+
     .       dd(4,i)*dtpopdze+dd(5,i)*dvflowze
        
c        write(*,*)'i,dd(1,i),dd(2,i),dd(3,i),dd(4,i),dd(5,i)',
c     &             i,dd(1,i),dd(2,i),dd(3,i),dd(4,i),dd(5,i)
c        write(*,*)'dxdze,dydze,dtempdze',dxdze,dydze,dtempdze
c        write(*,*)'dtpopdze,dvflowze',dtpopdze,dvflowze
c        write(*,*)'dddz',dddz  

        dddr = dddr+
cSm021114
c     .       dd(1,i)*dxdre+dd(2,i)*dydre+dd(3,i)*dtempdre+
     .       dd(1,i)*dxdre+dd(2,i)*dydre+dd(3,i)*dtempdre*1.d3+
     .       dd(4,i)*dtpopdre+dd(5,i)*dvflowre
        dddph= dddph+
cSm021114
c     .       dd(1,i)*dxdphie+dd(2,i)*dydphie+dd(3,i)*dtempdpe+
     .       dd(1,i)*dxdphie+dd(2,i)*dydphie+dd(3,i)*dtempdpe*1.d3+
     .       dd(4,i)*dtpopdpe+dd(5,i)*dvflowpe
      
	dxdwe=-2.d0*x_ar(i)/wf
        dydwe=-y_ar(i)/wf
        dddw=dddw+dd(1,i)*dxdwe+dd(2,i)*dydwe

      enddo

c      write(*,*)'dddz,dddr,dddph',dddz,dddr,dddph
c      write(*,*)'dddcnz,dddcnr,dddcm',dddcnz,dddcnr,dddcm
c      write(*,*)'dddw',dddwc

c      write(*,*)'hotdervs end'

      return
      end subroutine hotdervs

      double complex function hotnp(nbulk,ibw,cnper_0,
     .cnper2p,cnper2m,K,iraystop)
c     calculates n_perp( and the hot dielectric tensor K)
c     from the hot non-relativistic electron+ions
c     dispersion function.
c-------------------------------------------------------------------
c     INPUTS:
c     nbulk              the total number of plasma species
c     ibw  =0 the calculation the roots for O or X modes
c          =1 the calculation the root for BW
c     cnper_0 the approximation of N_perp, O or X mode from the cold plasma
!             according to ioxm parameter
c     cnper2p,cnper2m are N_perp**2 roots of the cold plasma dispersion 
c     the data from common block /nperpcomp/
c     massc(nbulka)       the masses of the plasma species
c     nllc=N_parallel 
c     xc(nbulka) = (fpe/f)**2
c     yc(nbulka) = fce/f
c     tc(nbulka) = averaged temperature in eV
c     nllc - parallel index of refraction n.
c     tpopc(nbulka)=T_perp/T_parallel
c     vflowc(nbulka)=flow velocity
c     thus function calls dhot and uses bisection to find the hot root
c     i.e. dhot_sum(nbulk,mass_ar,x_ar,y_ar,tav_ar,tpop_ar,
c     . vflow_ar,nll_in,np_in,iherm,K)=0
c
c     OUTPUTS:
c       function hotnp=nperp (double complex)
c       K double complex dielectric tensor
c       iraystop=0 the root have been found,
c       iraystop=1 the root was not found
c-------------------------------------------------------------
      implicit none
c-----input
      integer nbulk,ibw 
      double precision cnper_0,cnper2p,cnper2m
      include 'param.i'
      include 'nperpcom.i'
c-----external 
      double  precision rtbis,ddwrap
      double complex dhot_sum
      external rtbis,ddwrap
c-----output
      integer iraystop
      double complex K(3,3)
c-----local      
      double complex i,dK(3,3,7),dd(7),d,ddn
      integer ntry,j
      double precision x1,x2,xacc,nperp_x,nperp,npi,delta_n
      double precision cnpermax,p1,p2 
      double precision x,h_nperp,d_left,d_right,n_perp_l,n_perp_r

cfor test
      double precision step,x_test,d_test
    
      integer imode,jdelta_n,j_max,jmax,imoden,irootm,irootp
      integer icheck_d_left_good, icheck_d_right_good !local
c-----the number to extend the interval for the determination of the 
c     n_perpendicular root from the hot dispersion function
c     The hot plasma (X or O mode) root will be searched at the interval
c     [N_perpendicular_left, N_perpendicular_right]
c     N_perp_left =N_perp_left_cold_plasma  - jdelta_n*delta_n
c     N_perp_right=N_perp_right_cold_plasma + jdelta_n*delta_n
c     delta_n=N_perp_right_cold_plasma-N_perp_left_cold_plasma
c
c     j_max is the number of steps in the N_perp direction
c     in the interval delta_n, to make the table D(N_perp_j) 
c     for the hot plasma root calculations. 
      jdelta_n=10
      j_max=100        

      write(*,*)'in hotnp nllc=',nllc
      i = (0.0d0,1.0d0)
      xacc=1.d-14  ! accuracy of the root calculations
      ntry=200
      
      
      
c-----calculations of O or X mode , using the cold plasma root
c     cnper_0 as the initial approximation 
      if (ibw.eq.0) then !-------------------------------------ibw=0
         irootm=1
         if(cnper2m.lt.0.d0) then
           cnper2m=0.d0
           irootm=-1
         endif

         irootp=1
         if(cnper2p.lt.0.d0) then
           cnper2p=0.d0
           irootp=-1
         endif

         !Initial range for the root search with bisections:
         n_perp_l=dmin1(dsqrt(dabs(cnper2m)),dsqrt(dabs(cnper2p)))
         n_perp_r=dmax1(dsqrt(dabs(cnper2m)),dsqrt(dabs(cnper2p)))
         if(n_perp_l.le.0.d0 .and. n_perp_r.le.0.d0)then
            !write(*,*)' hotnp-1: n_perp_l, n_perp_r=', n_perp_l,n_perp_r
            ! YuP 120501
            n_perp_l=abs(nllc) ! just a guess; to be apart from 0.
            n_perp_r=1. ! just a wild guess
            !write(*,*)' hotnp-2: n_perp_l, n_perp_r=', n_perp_l,n_perp_r
         endif

ctest 070730
c         write(*,*)'in hotnp ibw=',ibw
!         write(*,*)'in hotnp cnper_m,cnper_p',
!     +              dsqrt(dabs(cnper2m)),dsqrt(dabs(cnper2p))
         write(*,*)'in hotnp n_perp_l,n_perp_r',n_perp_l,n_perp_r

         d_left= dreal(dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     &                 vflowc,nllc,n_perp_l,1,K))
         d_right=dreal(dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     &                 vflowc,nllc,n_perp_r,1,K))
         write(*,*)'hotnp d_left,d_right',d_left,d_right
         
         !YuP[2020-11-25] Adjusted logic for searching the proper root,
         !which is problematic in near vacuum (when two roots are very close).
         icheck_d_left_good=  0 ! to be determined/adjusted
         icheck_d_right_good= 0 ! to be determined/adjusted
         if(dabs(d_left).lt.1.d-7)then ! Good, within accuracy
            icheck_d_left_good=  1
         endif
         !But also check the other root - maybe it gives D()<1e-7, 
         !and maybe it is even better
         if(dabs(d_right).lt.1.d-7)then
            icheck_d_right_good= 1 
         endif
         !Now the logic itself:
         if(icheck_d_left_good+icheck_d_right_good.eq.2)then
           !both roots are good. Select one that is closer 
           !to target cold value cnper_0 
           !(from input: O or X mode according to ioxm parameter)
           if(dabs(n_perp_l-cnper_0).lt.dabs(n_perp_r-cnper_0))then
             nperp=n_perp_l
             goto 30 !->root is found (got lucky)
           else
             nperp=n_perp_r
             goto 30 !->root is found (got lucky)
           endif
         elseif(icheck_d_left_good.eq.1)then !only this one is good
           nperp=n_perp_l
           goto 30 !->root is found (got lucky)
         elseif(icheck_d_right_good.eq.1)then !only this one is good
           nperp=n_perp_r
           goto 30 !->root is found (got lucky)
         else ! none is good
           continue
         endif
         
c_end_test 070730

         delta_n=dabs(n_perp_r - n_perp_l)
         if (delta_n.lt.1.d-5) then ! two roots are almost same
            !YuP120430 Try small interval around 0.5*(n_perp_l+n_perp_r)
            jmax=1000
            nperp= 0.5*(n_perp_l+n_perp_r)
            x1=nperp*0.95
            x2=nperp*1.05
            h_nperp=(x2-x1)/(jmax-1)
            d_left=dreal(dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     .                   vflowc,nllc,x1,1,K))
            do j=2,jmax
               x2=x1+h_nperp
               d_right=dreal(dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     .                   vflowc,nllc,x2,1,K))
               !write(*,*)'x  D(x)=',x2,d_right 
               
               icheck_d_left_good=  0 ! to be determined/adjusted
               icheck_d_right_good= 0 ! to be determined/adjusted
               if(dabs(d_left).lt.1.d-7)then ! Good, within accuracy
                icheck_d_left_good=  1
               endif
               !But also check the other root - maybe it gives D()<1e-7, 
               !and maybe it is even better
               if(dabs(d_right).lt.1.d-7)then
                icheck_d_right_good= 1 
               endif
               !Now the logic itself:
               if(icheck_d_left_good+icheck_d_right_good.eq.2)then
                 !both roots are good. Select one that is closer 
                 !to target cold value cnper_0 
                 !(from input: O or X mode according to ioxm parameter)
                 if(dabs(x1-cnper_0).lt.dabs(x2-cnper_0))then
                   nperp=x1
                   goto 30 !->root is found (got lucky)
                 else
                   nperp=x2
                   goto 30 !->root is found (got lucky)
                 endif
               elseif(icheck_d_left_good.eq.1)then !only this one is ok
                 nperp=x1
                 goto 30 !->root is found (got lucky)
               elseif(icheck_d_right_good.eq.1)then !only this one is ok
                 nperp=x2
                 goto 30 !->root is found (got lucky)
               else ! none is good
                 continue
               endif

               if (d_left*d_right.le.0.d0) then
                  goto 10 !-> root is between x1 and x2; find exact root
               endif
               x1=x2 ! for the next step
               d_left=d_right ! for the next step
            enddo
         endif

c--------determination of the boundaries for the hot plasma dispersion
c        for the n_perp solver
         imode=0  !initialization
         imoden=0 !the number of hot plasma root
         
         x=dmin1(n_perp_l,n_perp_r)
        
         icheck_d_left_good=  0 ! to be determined/adjusted
         icheck_d_right_good= 0 ! to be determined/adjusted

         if (dabs(n_perp_l-cnper_0).lt.1.d-7) then
c--------imode=1, cold plasma mode (target value cnper_0) coincides with the left
c                 boundary of cold of plasma interval
            imode=1
            d_left=dreal(dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     .                   vflowc,nllc,cnper_0,1,K))
            write(*,*)'hotnp: n_perp_l=cnper_0=',cnper_0,
     &                '  D(cnper_0)=',d_left 
            if (dabs(d_left).lt.1.d-7) then
               icheck_d_left_good=  1
               nperp=cnper_0 ! or n_perp_l ?
               goto 30 !->root is found (got lucky)
            endif
         endif
         
         if (dabs(n_perp_r-cnper_0).lt.1.d-7) then
c--------imode=2, cold plasma mode (target value cnper_0) coincides with the right
c                 boundary of cold of plasma interval
            imode=2
            d_right=dreal(dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     .                   vflowc,nllc,cnper_0,1,K))
            write(*,*)'hotnp: n_perp_r=cnper_0=',cnper_0,
     &                '  D(cnper_0)=',d_right 
            if (dabs(d_right).lt.1.d-7) then
               icheck_d_right_good= 1
               nperp=cnper_0 ! or n_perp_r ?
               goto 30 !->root is found (got lucky)
            endif
         endif
        
        !pause
        
         if((irootm.eq.-1).or.(irootp.eq.-1)) then
c          we have only one positive cold plasma root N_perp
           imode=1
         endif

         write(*,*)'hotnp imode',imode
         write(*,*)'hotnp n_perp_l,cnper_0,n_perp_r',
     .   n_perp_l,cnper_0,n_perp_r
         
         h_nperp=delta_n/dfloat(j_max) 
         !note: delta_n= abs(n_perp_r - n_perp_l), j_max=100, jdelta_n=10
         
         nperp=x+(-jdelta_n*j_max)*h_nperp !can start with -10*delta_n
         !note: x=dmin1(n_perp_l,n_perp_r)  Can be 0.
         write(*,*)'forest hotnp x,jdelta_n,j_max,nperp',
     +   x,jdelta_n,j_max,nperp
c         write(*,*)'hotnp nbulkc,nllc',nbulkc,nllc
c         write(*,*)'hotnp massc',massc
c         write(*,*)'hotnp xc',xc
c         write(*,*)'hotnp yc',yc
c         write(*,*)'hotnp tec',tec
c         write(*,*)'hotnp tpopc',tpopc
c         write(*,*)'hotnp wflowc',vflowc
       
         if (nperp.lt.0.d0) then
           write(*,*)'nperp<0'
           nperp=1.d-5
           write(*,*)'nperp',nperp
         endif

         d_left=dreal(dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     .      vflowc,nllc,nperp,1,K))
c         write(*,*)'hotnp d_left',d_left 
         if (dabs(d_left).lt.1.d-7) goto 30 !->root is found (got lucky)
  
         do j=-jdelta_n*j_max+1,(jdelta_n+1)*j_max

            nperp=x+j*h_nperp

            if (nperp.lt.0.d0) then
               !-YuP Added:  re-define d_left using nperp=0.
               d_left=dreal(dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     .                      vflowc,nllc,0.d0,1,K)) 
               goto 20 !-> next j
            endif
 
            d_right=dreal(dhot_sum(nbulkc,massc,xc,yc,tec,tpopc,
     .                    vflowc,nllc,nperp,1,K))
            if (dabs(d_right).lt.1.d-7) goto 30 !->root is found (got lucky)
            !write(*,*)'j,nperp,d_left,d_right', j,nperp,d_left,d_right

            if(d_left*d_right.lt.0.d0) then
c-------------the hot dispersion function has the different signs 
c             at the interval  nperp={x+(j-1)*h_nperp, x+j*h_nperp}
              imoden=imoden+1 ! the number of the hot plasma root
              !write(*,*)'d_left*d_right.lt.0.d0 imoden=',imoden
              d_left=d_right
              if (imoden.eq.imode) then
                write(*,*)'hotnp imoden=imode'
                x1=nperp-h_nperp
                x2=nperp
                goto 10 !-> finished with j-scan
                !Found such x1 and x2 that d_left*d_right<0
              endif
            else
              d_left=d_right
            endif
            
 20      continue    
         enddo ! j-scan

 10      continue ! exit handle from j-loop, when 
          
         
         if (x1.lt.0.d0) x1=1.d-5
         write(*,*)'hotnp before rtbis x1,x2,xacc',x1,x2,xacc
         x1=dlog(x1)
         x2=dlog(x2)
         
         nperp=dexp(rtbis(ddwrap,x1,x2,xacc)) ! O or X mode root
         write(*,*)'hotnp after rtbis; O or X modes nperp=',nperp   
         !pause

 30      continue  ! exit handle for "root is found"

         d=dhot_sum(nbulk,massc,xc,yc,tec,tpopc,
     .    vflowc,nllc,nperp,2,K)

         write(*,*)'in hotnp O or X mode 2 complete d=',d   
         d=dhot_sum(nbulk,massc,xc,yc,tec,tpopc,
     .    vflowc,nllc,nperp,1,K)
         write(*,*)'in hotnp O or X mode 1 hermitian d=',d   
ccc         iraystop=0
      endif ! ibw=0 -------------------------------------------ibw=0
      
      
c-----calculation of the BW
      if(ibw.eq.1)then !---------------------------------------ibw=1
c       x1,x2 initialization
        cnpermax=dmax1(cnper2m,cnper2p)
        if (cnpermax.gt.0.d0) then 
           nperp_x=dsqrt(cnpermax)
        else
           nperp_x=0.d0
        endif

c       if ((cnper2m.lt.0.d0).and.(cnper2p.lt.0.d0)) nperp_x=0.d0 

        write(*,*)'forest.f nperp_x',nperp_x

        x1 = dlog(nperp_x+1.d-2)
        write(*,*)'forest.f hotnp x1',x1
c        x2 = dlog(1.d0)         !nperp
        x2 = x1+dlog(1.d0)         !nperp

c       expand the range of x2 until the first root is found
        iraystop=0

        do j = 1, ntry
          if (j.eq.ntry) then
            iraystop=1
            write(*,*)'hotnperp could not find ebw root'
            return
          endif
          if (ddwrap(x2)*ddwrap(x1).lt.0) go to 12
          x1 = x2
          if (ddwrap(x2)*ddwrap(x1).ge.0) x2 = x2 +dlog(2.0D0)
        end do

 12     continue

        write(*,*)'forest.f hotnp before dexp(rtbis)'
        nperp=dexp(rtbis(ddwrap,x1,x2,xacc))
        write(*,*)'in hotnp nperp=',nperp   
        d=dhot_sum(nbulk,massc,xc,yc,tec,tpopc,
     .    vflowc,nllc,nperp,2,K)

        write(*,*)'in hotnp 2 complete d=',d   
        d=dhot_sum(nbulk,massc,xc,yc,tec,tpopc,
     .    vflowc,nllc,nperp,1,K)
        write(*,*)'in hotnp 1 hermitian d=',d
   
      endif   !ibw=1 BW ---------------------------------------ibw=1
       
c-----------
      npi=0.d0
      hotnp = dcmplx(nperp,npi)

c      write(*,*)'in hotnp=',hotnp
      return
      end function hotnp




      double complex function hotnperp(z,r,phi,nll,cnteta,cnphi,
     .K,iraystop)
c     calculates n_perp( and the hot dielectric tensor K)
c     from the hot non-relativistic electron+ions
c     dispresion function.
c     if ibw=0 it calculates the O and X modes , using the cold plasma roots
c              as the initial approximation
c     if ibw=1 it calculates BW root N_perp_bw> N_perp_O and N_perp_X
c--------------------------------------------------------------------
c     INPUTS:    
c     z,r,phi - space cordinates
c     nll     - parallel index of refraction n
c     cnteta  _ N component parallel to the poloidal direction
c     cnphi    _ N component parallel to the toroidal direction
c    
c     OUTPUTS:
c       function hotnperp=N_perpendicular (double complex)
c       K double complex dielectric tensor
c       iraystop=0 the root have been found,
c       iraystop=1 the root was not found
c-----------------------------------------------------    
      !implicit integer (i-n), real*8 (a-h,o-z) 
      implicit none
      real*8 cnr,cnz,cm,cnphi,cnteta,cnpar2,cntang2
      real*8 cnper2p,cnper2m,cnper_0
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'nperpcom.i'

c-----input
      double precision z,r,phi
      double precision nll ! N_parallel

c-----external 
      double complex hotnp
      double precision b,x,y,tempe,tpoprho,vflowrho
      external hotnp,b,x,y,tempe,tpoprho,vflowrho
c-----local  
     
      double precision mass_ar(nbulka),x_ar(nbulka),y_ar(nbulka),
     .t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka)
      integer j,id_old

c-----output      
      double complex K(3,3)
      integer iraystop

      cnpar2=nll*nll
      cntang2=cnteta*cnteta+cnphi*cnphi

c-----calculation of two roots N_perp**2(N_parallel)
c     from the cold electron+ions plasma dispersion function
      id_old=id
      id=1            ! to use the cold plasma dispersion
      cnper2p=-1.d0
      cnper2m=-1.d0 

      !write(*,*)'in hotnperp before npernpar'
      call npernpar(z,r,phi,nll,cnper2p,cnper2m)     
      !write(*,*)'in hotnperp after npernpar'


      id=id_old      
      write(*,*)'forest.f in hotnperp nll,cnper2p,cnper2m',
     .nll,cnper2p,cnper2m
      write(*,*)'htnperp ibw=',ibw
      
      iraystop=0
      cnper_0=0.d0 ! YuP[08-2017] Initialize

c-----choosing of the cold mode N_perp=cnper_0 according to ioxm parameter      
      if(ibw.eq.0) then !it is not the direct launch of Bernstein wave
        id_old=id
        id=1            ! to use the cold plasma dispersion

        !write(*,*)'hotnperp:before cninit12'
        call cninit12(z,r,phi,nll,cnteta,cnphi,
     .               cnz,cnr,cm,iraystop)
        !write(*,*)'hotnperp:After cninit12',iraystop

        id=id_old  

        if(iraystop.eq.0) then
          cnper_0=dsqrt(cnz*cnz+cnr*cnr+(cm/r)**2-nll*nll)
          write(*,*)'hotnperp after cninit12 cnper_0,cnper_02',
     *    cnper_0,cnper_0**2
        else
          write(*,*)'in hotnperp iraystop=1, hotnperp could not find 
     .    the root for the cold plasma mode'
        endif
        !YuP: what is cnper_0 in case of ibw=1?
      endif ! ibw=0

c-----initialization common/nperpcom/, it is in nperpcom.i
      nllc=nll
      nbulkc=nbulk
cSm030226
      if(nbulka.lt.nbulk) then
        write(*,*)'in forest.f in hotnperp nbulka.lt.nbulk'          
        write(*,*)'nbulka=',nbulka,'nbulk=',nbulk
	write(*,*)'change parameter nbulka in file param.i'  
        stop
      endif

c      write(*,*)'hotnperp before loop j nbulk=',nbuklk

      do j=1,nbulk
        massc(j)=dmas(j)
        xc(j)=x(z,r,phi,j)
        yc(j)=y(z,r,phi,j)
cSm000324
        if(j.eq.1) yc(1)=-yc(1) ! negative Y=(omega_ce/omega) for electrons
        tec(j)=tempe(z,r,phi,j)*1.d+3 !(eV) averaged temperature
        tpopc(j)=tpoprho(rho,j)
        vflowc(j)=vflowrho(rho,j)         
      enddo

c      write(*,*)'hotnperp after  loop j=1,nbulk'1
      write(*,*)'in hotnperp before hotnp cnper_0**2,cnper2p,cnper2m',
     .cnper_0**2,cnper2p,cnper2m
      
      hotnperp=hotnp(nbulk,ibw,cnper_0,cnper2p,cnper2m,K,iraystop)
      
      write(*,*)'after hotnp hotnperp=',hotnperp,'iraystop=',iraystop
      return
      end function hotnperp


c        subroutine ibess0 
c        purpose                                                     
c            compute the modified bessel function Io of order zero    
c                                                                    
c        usage                                                       
c            call ibess0(x,ri0)                                          
c                                                                    
c        description of parameters                                   
c            x    -given argument of the bessel function I of order 0
c            ri0  -resultant value of the bessel function I of order 0

CBF   This function returns exp(-x) * Io

c                                                                     
c        remarks                                                      
c            large values of the argument may cause overflow in the   
c            builtin exp-function                                     
c                                                                     
c        subroutines and function subprograms required                
c           none                                                      
c                                                                     
c        method                                                       
c           polynomial approximations given by e.e. allen are used for
c           calculation.                                              
c           for reference see                                         
c           m. abramowitz and i.a. stegun,'handbook of mathematical   
c           functions', u.s. department of commerce, national bureau of
c           standards applied mathematics series, 1966, p.378.         
c                                                                      
c     .................................................................
c                                                                      
      subroutine ibess0(x,ri0)
      implicit none
      double precision x,ri0,expxx,z,pa
      ri0=dabs(x)
      if(ri0-3.75d0)1,1,2
    1 z=x*x*7.111111d-2
      expxx = dexp(-ri0) ! Here ri0 is not large (no underflow)
      ri0=((((( 4.5813d-3*z+3.60768d-2)*z+2.659732d-1)*z+1.206749d0)*z
     &   +3.089942d0)*z+3.515623d0)*z+1.d0
      ri0=ri0*expxx
      return
    2 z=3.75d0/ri0 !Here ri0 can be very large, but we don't have exp(-x)
      pa=1.0d0/dsqrt(ri0)*((((((((3.92377d-3*z-1.647633d-2)*z
     & +2.635537d-2)*z-2.057706d-2)*z+9.16281d-3)*z-1.57565d-3)*z
     & +2.25319d-3)*z+1.328592d-2)*z+3.989423d-1)
      ri0 = pa
      return
      end subroutine ibess0
c
c     .................................................................
c
c        subroutine ibessn
c
c        purpose
c           compute the modified bessel functions i for orders 1 to n  
c
c        usage
c           call ibessn(x,n,zi,ri) !->  I_n(x)*exp(-x)
c
c        description of parameters
c           x     -given argument of the bessel functions i    
c           n     -given maximum order of bessel functions i       
c           zi    -given value of bessel function i of order zero  
c
c           ri    -resultant vector of dimension n, containing the 
!                  YuP:  I_n(x)*exp(-x)
c
c
c        remarks
c           the value of zi may be calculated using subroutine i0.     
c           using a different value has the effect that all values of  
c           bessel functions i are multiplied by the  factor zi/i(0,x) 
c           where i(0,x) is the value of i for order 0 and argument x. 
c           this may be used disadvantageously if only the ratios of i 
c           for different orders are required.
c
c        subroutines and function subprograms required
c           none
c
c        method
c           the values are obtained using backward recurrence relation 
c           technique. the ratio i(n+1,x)/i(n,x) is obtained from a    
c           continued fraction.
c           for reference see
c           g. blanch,'numerical evaluation of continued fractions',   
c           siam review, vol.6,no.4,1964,pp.383-421.
c
c     .................................................................
c
      subroutine ibessn(x,n,zi,ri)
      implicit none 
      double precision x,zi,ri,a0,a1,b0,fn,fi,a,b,b1,q0,q1,an
      integer*4 i,j,n,k
      dimension ri(*)

      if(n)10,10,1
    1 fn=dfloat(n+n)
      q1=x/fn
      if(dabs(x)-5.d-4)6,6,2
    2 a0=1.d0
      a1=0.d0
      b0=0.d0
      b1=1.d0
      fi=fn
    3 fi=fi+2.d0
      an=fi/dabs(x)
      a=an*a1+a0
      b=an*b1+b0
      a0=a1
      b0=b1
      a1=a
      b1=b
      q0=q1
      q1=a/b
      if(dabs((q1-q0)/q1)-1.d-6)4,4,3
    4 if(x)5,6,6
    5 q1=-q1
    6 k=n
    7 q1=x/(fn+x*q1)
      ri(k)=q1
      fn=fn-2.d0
      k=k-1
      if(k)8,8,7
    8 fi=zi
      do 9 i=1,n
      fi=fi*ri(i)
    9 ri(i)=fi
   10 return
      end subroutine ibessn


      subroutine nllcnp_s(z,r,phi,cnz,cnr,cm,nll,cnp)
c     for testing calcultes N_perp=cnp and N_parallel=nll
c     use bz,br,bphi,bmod from common one.i
      implicit none
      include 'param.i'
      include 'one.i'
c     input
      double precision z,r,phi,cnz,cnr,cm
c     output
      double precision nll,cnp
c-----externals
      double precision b 
      bmod=b(z,r,phi)
      nll=(cnz*bz+cnr*br+cm*bphi/r)/bmod
      cnp=dsqrt(dabs(cnz**2+cnr**2+(cm/r)**2-nll**2))
      return
      end subroutine nllcnp_s
c=====================================================================
 
      subroutine npnllmin(nll_in,np_in,nll_cor,np_cor)
c     gives the minimal values for n_paral and n_perp near the zero value
      implicit none
c     input      
      double precision nll_in,np_in !parallel and perpendicular N components 
c     output
      double precision nll_cor,np_cor !corrected value of N_par and N_perp
c     local parameters
      double precision nll_min,np_min

c      nll_min=1.d-4   ! should be positive
c      np_min=1.d-3    ! should be positive
      nll_min=1.d-7   ! should be positive
      np_min=1.d-7    ! should be positive
cSAP090725
c      nll_min=1.d-11   ! should be positive

      np_cor=np_in
      nll_cor=nll_in
      !--------------------------------------
      if (np_in.le.dabs(np_min)) then
        np_cor = dabs(np_min)
        write(*,*)'npnllmin:Nperp: np_in, np_cor',np_in,np_cor
        !pause
      endif
      !--------------------------------------
      if (dabs(nll_in).le.dabs(nll_min)) then
         if (nll_in.ge.0.d0) then
            nll_cor=dabs(nll_min)
         else 
            nll_cor=-dabs(nll_min)
         endif 
         write(*,*)'npnllmin:Npar: nll_in, nll_cor',nll_in,nll_cor
         !pause
      endif
      !--------------------------------------
      return
      end subroutine npnllmin
c=====================================================================
      
      subroutine put_mass(mass_ar,nbulk)
c     puts mass_ar(nbulk)=dmass(nbulka) from commom/ioons.i/
      
	implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
	include 'ions.i'
c-----input
      integer nbulk   
c-----output
      double precision mass_ar(*)
c-----local 
      integer i	
	   
cSm030226
      if(nbulk.gt.nbulka) then
	   write(*,*)'in forest.f in hotdervs nbulk.gt.nbulka'
	   write(*,*)'put the value of parameter nbulka.ge.nbulk'
	   write(*,*)'parameter nbulka is given  in /param.i/'
	   stop
      endif

      do i=1,nbulk
	 mass_ar(i)=dmas(i)
      enddo
	
      return
      end subroutine put_mass




      



      subroutine xmode(nll,x,y,nperp_x,icut)
c     calculates N_perp(N_parallel) for x mode in cold electron plasma
      implicit none

c     input
      double precision nll,x,y 
c     nll=N_parallel
c     x=(omega_pl_e/omega)**2
c     y=(omega_c_e/omega)   ! the sign of Ye is the same as for ions Y-i            

c     output
      double precision nperp_x
      integer icut
c     nperp_x=N_perp for X mode
c     icut =0 if x_mode exists, =1 x_mode can not exist
      
c     local
      double precision ncld(2),nz,az2,an2(2)
      integer icuto(2)
      nz=nll  
      az2=nz**2 
      call coldm(nz,az2,x,y,an2,ncld,icuto)

      
      if(icuto(1).eq.1) write(*,*)'in xmode cutoff omode'  
      if(icuto(2).eq.2) then
         write(*,*)'in xmode cutoff xmode'
         icut=1
      else
         icut=0
         nperp_x=ncld(2)
      endif

      return
      end subroutine xmode


      subroutine harmon_z(X,Y,nll,vflow,c,beta,lambda,nmax,
     .n_harm1,n_harm2,n,istop)
c-----calculates:
c     The number n_harm0 of the main resonance gyro-harmonic
c     from the condition 1-N_parallel(V/C)=n_harm0*Y  
c     The boundary of gyro-harmonic interval h_harm1,n_harm_2
c     as the roots of the following equation  
c     abs{(1-N_parallel(V/C)+n_harm*Y)/(N_parallel*beta_parallel)}=1/(2*eps_z)
c     Index: istop=1, it is possible to calculate the susceptibility tensor
c            istop=0, it is impossible to calculate the susceptibilty tensor
c     n      the total number of the used Bessell functions 
      implicit none
c-----input
      double precision Y      ! Y_e should have the sign opposite to Y_i
      double precision X
      double precision nll    ! N_parallel
      double precision beta   ! vte_longitudinal/c 
      double precision lambda !0.5 k_perp**2 *rho_gyro**2
      double precision vflow  ! the drift velocity (cm/sec)
      double precision c      ! the light velocity (cm/sec)
      integer nmax            ! the max number of possible Bessel functions
c-----output  
      integer n_harm0,n_harm1,n_harm2 
      integer istop,n                   
c-----locals
      double precision ph1,ph2,ph11,ph22
      double precision eps_z    ! the abs accuracy for Z_0 function
      double precision eps_I    ! the amplitude for  modified bessel functions
      double precision ln_I     !=ln(I) the estimation of bessel function 
c     n_delta is the maximal  difference between the main harmonic n_harm0 and h_harm1 and h_harm2
      integer n_delta
      parameter (n_delta=5) 

      eps_z=1.d-14
      eps_I=1.d-14
      
      eps_z=1.d-7
      eps_I=1.d-7

      eps_z=eps_z/X    
      istop=1
      
c      write(*,*)'harmon_z nmax,X,Y',nmax,X,Y
      
      if (Y.eq.0.d0) then
        n_harm1=0
        n_harm2=nmax-2
        n_harm0=n_harm2/2
      else
        n_harm0=idint((1.d0-nll*(vflow/c))/Y)
c        write(*,*)'harmon_z nll,Y,vflow,n_harm0',nll,Y,vflow,n_harm0

c-------The estimation of the modified bessel function I_(n_harm0)
c        write(*,*)'harm_z lamda',lambda
        ln_I=abs(n_harm0)*dlog(lambda/abs(n_harm0))+
     .       lambda+0.5d0*dlog(0.5d0*abs(n_harm0))
 
        if(ln_I.lt.dlog(eps_I))then
c--------- modified Bessel function I_n_harm0<eps_I
c           write(*,*)'in harmon_z: modified Bessel I_n_harm0<eps_I'
c           write(*,*)'n_harm0=',n_harm0,'nmax=',nmax
c           write(*,*)'log(I_n_harm0)=',ln_I,'eps_I',eps_I
c          istop=0
c          goto 10
        endif
           
c        write(*,*)'nll,beta', nll,beta          
        ph1=(1.d0-nll*(vflow/c)-nll*beta/(2.d0*eps_z))/Y
        ph2=(1.d0-nll*(vflow/c)+nll*beta/(2.d0*eps_z))/Y

c-------the corrections to get the integer number I smaller
c       than the maximal integer for g77
        if(ph1.lt.-2147483640.d0) ph1=-2147483640.d0
        if(ph1.gt. 2147483640.d0) ph1= 2147483640.d0
        if(ph2.lt.-2147483640.d0) ph2=-2147483640.d0
        if(ph2.gt. 2147483640.d0) ph2= 2147483640.d0

        n_harm1=idint(dmin1(ph1,ph2))-1
        n_harm2=idint(dmax1(ph1,ph2))+1
c        write(*,*)'harmon_z n_harm1,n_harm2',n_harm1,n_harm2
c        write(*,*)'harmon_z Y,beta,nll,nll*beta/(2.d0*eps_z)',
c     .  Y,beta,nll,nll*beta/(2.d0*eps_z)
c        write(*,*)'harmon_z n_harm1,n_harm2',n_harm1,n_harm2
c        n=max(abs(n_harm1),abs(n_harm2))
  
      endif

      if(nmax.lt.3) then
         write(*,*)'in harmon_z the total number of the possible Bessel'
         write(*,*)' functions nmax <3. Increase it in the routine,'
         write(*,*)' which calls harmon_z.'
         stop
      endif
     
c      write(*,*)'nharmon_z(0) n_harm1,n_harm2',n_harm1,n_harm2
c-----(1) 0 =< n_harm1 =< n_harm2
      if ((0.le.n_harm1).and.(0.le.n_harm2)) then
c      write(*,*)'nharmon_z(1)n_harm1,n_harm2,nmax',n_harm1,n_harm2,nmax
         if((n_harm2).le.(nmax-2)) then
c-----------n_harm2 =<  nmax-2          
            n=n_harm2
c            write(*,*)'nharm 1.a n,n_harm2,n_harm1',n,n_harm2,n_harm1
            goto 10
         endif

         if((n_harm1.le.(nmax-2)).and.((nmax-2).lt.n_harm2)) then
c----------n_harm1 =< nmax-2 < n_harm2
           n_harm2=nmax-2
           n=n_harm2
c           write(*,*)'nharm 1.b n,n_harm2,n_harm1',n,n_harm2,n_harm1
           goto 10
         endif
          
         if((nmax-2).lt.n_harm1) then
c----------nmax-2 < n_harm1
           write(*,*)'nharm 1c istop',istop
           istop=0      
           goto 10
         endif

      endif !(1)

c-----(2) n_harm1 < 0 =< n_harm2
      if ((n_harm1.lt.0).and.(n_harm2.ge.0)) then
          
c--------(2.1) 0 < | n_harm1 | =< n_harm2
         if(abs(n_harm1).le.n_harm2) then
   
c            write(*,*)'(2.1)'

           if(n_harm2.le.(nmax-2)) then
c------------n_harm2 =< nmax-2
             n=n_harm2

c            write(*,*)' 2.1 1 harmon_z:n,n_harm1,n_harm2', 
c     &                   n,n_harm1,n_harm2   
             goto 10
           endif

           if((abs(n_harm1).le.(nmax-2)).and.
     .        ((nmax-2).lt.n_harm2)) then
c------------| n_harm1 | =< nmax-2 < n_harm2
             n_harm2=nmax-2
             n=n_harm2

c             write(*,*)' 2.1 2 harmon_z:n,n_harm1,n_harm2', 
c     &                   n,n_harm1,n_harm2   
             goto 10
           endif

           if((nmax-2).lt.abs(n_harm1)) then
c------------nmax-2 < | n_harm1 |
             n_harm1=-(nmax-2)
             n_harm2=nmax-2
             n=n_harm2   

c             write(*,*)' 2.1 3 harmon_z:n,n_harm1,n_harm2', 
c     &                   n,n_harm1,n_harm2   

             goto 10
           endif

c           write(*,*)'harmon_z:n,n_harm1,n_harm2', n,n_harm1,n_harm2

         endif !(2.1)

c--------(2.2) 0 < n_harm2 < | n_harm1 | 
         if(n_harm2.lt.abs(n_harm1)) then
c           write(*,*)'harmon_z 2.2'
           if(abs(n_harm1).le.(nmax-2)) then
c------------| n_harm1 | =< nmax-2
             n=abs(n_harm1)
             goto 10
           endif

           if((n_harm2.le.(nmax-2)).and.
     .        ((nmax-2).lt.abs(n_harm1))) then
c------------n_harm2  =< nmax-2 < | n_harm1 !
             n_harm1=-(nmax-2)
             n=abs(n_harm1)
             goto 10
           endif

           if((nmax-2).lt.n_harm2) then
c------------nmax-2 <  n_harm2
             n_harm1=-(nmax-2)
             n_harm2=nmax-2
             n=n_harm2    
             goto 10
           endif

         endif !(2.2)

      endif    !(2)

c-----(3)  n_harm1 =< n_harm2 < 0
      if ((n_harm1.lt.0).and.(n_harm2.lt.0)) then

         if(abs(n_harm1).le.(nmax-2)) then
c-----------| n_harm1 | =<  nmax-2          
            n=abs(n_harm1)
            goto 10
         endif

         if((abs(n_harm2).le.(nmax-2)).and.
     .      ((nmax-2).lt.abs(n_harm1))) then
c----------| n_harm2 | =< nmax-2 < | n_harm1 |
           n_harm1=-(nmax-2)
           n=abs(n_harm1)
           goto 10
         endif

         if((nmax-2).lt.abs(n_harm2)) then
c----------nmax-2 < | n_harm2|
c      write(*,*)'nmax-2<|n_harm2| n_harm1,2,nmax=',n_harm1,n_harm2,nmax
c      write(*,*) nll*beta/(2.d0*eps_z), Y
           istop=0      
           goto 10
         endif

      endif !(3)



 10   continue

c      write(*,*)'harmon_z Y,n_harm0,n_harm1,n_harm2,n',
c     .Y,n_harm0,n_harm1,n_harm2,n
     
c      if (abs(n_harm2-n_harm0).gt.n_delta) n_harm2=n_harm0+n_delta
c      if (abs(n_harm1-n_harm0).gt.n_delta) n_harm1=n_harm0-n_delta

c      n=7 

c      if(istop.eq.0) then,n
c         write(*,*)'in harmon_z istop=0, n_harm0',n_harm0
c         n_harm2=n_harm0+n_delta        
c         n_harm1=n_harm0-n_delta
c        write(*,*)'in harmon_z istop=0, n_harm1,n_harm0,n_harm2',
c     &  n_harm1,n_harm0,n_harm2
c        istop=1
c      endif
   
      return
      end subroutine harmon_z


      subroutine set_nperpcom(nll,nbulk,z,r,phi,dmas)
c-----set the common block /nperpcom.i for the function hotnp
      implicit none
c-----input 
      integer nbulk !is the number of the plasma species
      double precision nll !the parallel refructive index N_par
      double precision z,r,phi
      double precision dmas(*) 
c-----externals
      double precision x,y,tempe,tpoprho,vflowrho,fpsi,rhopsi
c-----output the data for common block 
      include 'param.i' 
      include 'nperpcom.i'
c-----local 
      integer j
      double precision rho,psi

      psi=fpsi(r,z)    ! or psif(z,r) ?
      rho=rhopsi(psi) !small radius

      nbulkc=nbulk
      nllc=nll

      if(nbulka.lt.nbulk) then
        write(*,*)'in forest.f in set_nperpcom nbulka.lt.nbulk'          
        write(*,*)'nbulka=',nbulka,'nbulk=',nbulk
	write(*,*)'change parameter nbulka in file param.i'
        stop  
      endif

      do j=1,nbulk
         massc(j)=dmas(j)
         xc(j)=x(z,r,phi,j)
         yc(j)=y(z,r,phi,j)
         if(j.eq.1) yc(1)=-yc(1) ! negative Y=(omega_ce/omega) for electrons
         tec(j)=tempe(z,r,phi,j)*1.d+3 !(eV) averaged temperature
         tpopc(j)=tpoprho(rho,j)
         vflowc(j)=vflowrho(rho,j)         
       enddo
      
       return
       end subroutine set_nperpcom
!

       subroutine dhot_sum_ei(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     . vflow_ar,nll_in,np_in,K_sum,dhot_sum_t,dhot_sum_e,dhot_sum_i)
c-----------------------------------------------------------
c     calculates the hot dielectric tensor  K_sum(3,3)
c     for the hot electron and ion plasma
c     T.Stix, Waves in plasmas,(1992), p.252 (40), p.258 (57)
c     K_sum=1+sum{s}K_s  
c-----------------------------------------------------------
c     INPUTS:
c      nbulk - the total number of the plasma species 
c      mass_ar - the masses of the plasma  species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=average temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar = T_perp/T_parallel  
c      vflow_ar    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c     OUTPUT:
c      K_sum(3,3) -  the nine components of the dielectric tensor
c                    evaluated at (mass,X,Y,Te,tpop,vflow,Nll,np)
c      dhot_sum    -  hot double complex dispersion function
c      dhot_sum_e
c      dhot_sum_e
c-----------------------------------------------------------
      implicit none
c     input
      integer nbulk,iherm   
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double precision np_in,nll_in
c-----output
      double complex K_sum(3,3)
      double complex dhot_sum_t,dhot_sum_e,dhot_sum_i
c-----local
      double precision np,nll,nps,nlls
      double complex K(3,3),K_h(3,3)! sucseptibilities
      double complex Ka_e(3,3), Ka_i(3,3) ! anty Hermitian electron and ions sucseptibilities
      double complex Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz
      double complex Kxx_e,Kxy_e,Kxz_e,Kyx_e,Kyy_e,Kyz_e,
     +Kzx_e,Kzy_e,Kzz_e
      double complex Kxx_i,Kxy_i,Kxz_i,Kyx_i,Kyy_i,Kyz_i,
     +Kzx_i,Kzy_i,Kzz_i

      integer j,j1,j2
c-----external DHOT_s,zeroK


c-----initialization of the dielectric tensor K_sum(3,3)
      call zeroK(K_sum)
      call zeroK(Ka_e)
      call zeroK(Ka_i)
      
      do j=1,nbulk      
       
        call DHOT_s(mass_ar(j),x_ar(j),y_ar(j),t_av_ar(j),tpop_ar(j),
     .  vflow_ar(j),nll_in,np_in,1,K_h) !hermitian tensor
       
        call DHOT_s(mass_ar(j),x_ar(j),y_ar(j),t_av_ar(j),tpop_ar(j),
     .  vflow_ar(j),nll_in,np_in,2,K)   !total tensor

        do j1=1,3
          do j2=1,3
            K_sum(j1,j2)=K_sum(j1,j2)+K_h(j1,j2) !hermitian 
            if (j.eq.1) then ! electron anty-Hermitian tensor
               Ka_e(j1,j2)=0.5*(K(j1,j2)-dconjg(K(j2,j1)))
            else             ! ions anty-Hermitian twensor
               Ka_i(j1,j2)=Ka_i(j1,j2)+0.5*(K(j1,j2)-dconjg(K(j2,j1)))
            endif
          enddo
        enddo
       
      enddo

      K_sum(1,1) = K_sum(1,1) + dcmplx(1.0D0,0.0D0)
      K_sum(2,2) = K_sum(2,2) + dcmplx(1.0D0,0.0D0)
      K_sum(3,3) = K_sum(3,3) + dcmplx(1.0D0,0.0D0)
      
      Kxx = K_sum(1,1)
      Kyy = K_sum(2,2)
      Kzz = K_sum(3,3)
      Kxy = K_sum(1,2)
      Kxz = K_sum(1,3)
      Kyx = K_sum(2,1)
      Kyz = K_sum(2,3)
      Kzx = K_sum(3,1)
      Kzy = K_sum(3,2)

      Kxx_e =Kxx+ Ka_e(1,1)
      Kyy_e =Kyy+ Ka_e(2,2)
      Kzz_e =Kzz+ Ka_e(3,3)
      Kxy_e =Kxy+ Ka_e(1,2)
      Kxz_e =Kxz+ Ka_e(1,3)
      Kyx_e =Kyx+ Ka_e(2,1)
      Kyz_e =Kyz+ Ka_e(2,3)
      Kzx_e =Kzx+ Ka_e(3,1)
      Kzy_e =Kzy+ Ka_e(3,2)

      Kxx_i =Kxx+ Ka_i(1,1)
      Kyy_i =Kyy+ Ka_i(2,2)
      Kzz_i =Kzz+ Ka_i(3,3)
      Kxy_i =Kxy+ Ka_i(1,2)
      Kxz_i =Kxz+ Ka_i(1,3)
      Kyx_i =Kyx+ Ka_i(2,1)
      Kyz_i =Kyz+ Ka_i(2,3)
      Kzx_i =Kzx+ Ka_i(3,1)
      Kzy_i =Kzy+ Ka_i(3,2)

c-----hot dispersion function
!YuP      call npnllmin(nll_in,np_in,nll,np) ! YuP[2020-08-24] in dhot_sum_ei: Having nll=0 is ok here
      nll=nll_in !YuP[2020-08-24] Having nll=0 is ok here 
      np=max(0.d0,np_in) !YuP[2020-08-24] Having np=0 is ok here
      nlls=nll**2
      nps=np**2

      dhot_sum_t =(Kxx-nlls) * (Kyy-nlls-nps) * (Kzz-nps)
     .+  Kxy * Kyz * (Kzx+np*nll) 
     .+ (Kxz+np*nll) * Kyx * Kzy
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
     .-  Kzy * Kyz * (Kxx-nlls)
     .- (Kzz - nps) * Kyx * Kxy

      dhot_sum_e =(Kxx_e-nlls) * (Kyy_e-nlls-nps) * (Kzz_e-nps)
     .+  Kxy_e * Kyz_e * (Kzx_e+np*nll) 
     .+ (Kxz_e+np*nll) * Kyx_e * Kzy_e
     .- (Kzx_e+np*nll) * (Kyy_e-nlls-nps) * (Kxz_e+np*nll)
     .-  Kzy_e * Kyz_e * (Kxx_e-nlls)
     .- (Kzz_e - nps) * Kyx_e * Kxy_e

      dhot_sum_i =(Kxx_i-nlls) * (Kyy_i-nlls-nps) * (Kzz_i-nps)
     .+  Kxy_i * Kyz_i * (Kzx_i+np*nll) 
     .+ (Kxz_i+np*nll) * Kyx_i * Kzy_i
     .- (Kzx_i+np*nll) * (Kyy_i-nlls-nps) * (Kxz_i+np*nll)
     .-  Kzy_i * Kyz_i * (Kxx_i-nlls)
     .- (Kzz_i - nps) * Kyx_i * Kxy_i


c      write(*,*)'dhot_sum_e=',dhot_sum_e

      return
      end subroutine dhot_sum_ei


      subroutine solv_nperp_hot(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     . vflow_ar,nll_in,np_in,iter_max,nperp_real,nperp_imag)
c     solve the complex hot plasma dispersion relation
c     calculate ReN_perp and ImN_perp using Newton method
c-----------------------------------------------------------
c     INPUTS:
c      nbulk - the total number of the plasma species 
c      mass_ar - the masses of the plasma  species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=average temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar = T_perp/T_parallel  
c      vflow_ar    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c      iherm =1 hermitian dielectric tensor, 2-full
c      iter_max the maximal number of iterations
c     OUTPUT:
c      nperp_real    real and imaginary parts of the refructve index
c      nperp_imag    root of the dispersion relation
c-----------------------------------------------------------
      implicit none
c     input
      integer nbulk,iter_max
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double precision np_in,nll_in
c-----output
      double precision nperp_real,nperp_imag
c-----local
      integer iherm,i
      double precision np,nll,nps,nlls
      double complex K_sum(3,3)
      double precision cnper_bot,cnper_top,cnprim_bot,cnprim_top,step
      double precision cnper_p,cnper_m,cnprim_p,cnprim_m,dp_r,dm_r,
     .  dp_i,dm_i,d_rd_d_rn,d_id_d_rn,d_rd_d_in,d_id_d_in,d_r_bot,
     .  d_i_bot,det,det_r,det_i,eps,deltn,eps1,alpha,
     .  d_r_top,d_i_top,norm_bot,norm_top 
      double complex d
           
c-----external DHOT_s,zeroK
      double complex dhot_sum,dhot_sum_c

c      write(*,*)'solv_nperp_hot'

      eps=1.d-5
      eps1=1.d-12
      eps1=1.d-6
      step=1.d-6

      cnper_bot=np_in
c      cnprim_bot=0.d0
      cnprim_bot=nperp_imag
      iherm=2
      i=1
c      write(*,*)'iter_max',iter_max
      d=dhot_sum_c(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     . vflow_ar,nll_in,cnper_bot,cnprim_bot,iherm,K_sum)
      d_r_bot=dreal(d)
      d_i_bot=dimag(d)
      norm_bot=dsqrt(d_r_bot**2+d_i_bot**2)
 10   continue


c      write(*,*)'i,cnper_bot,cnprim_bot,d',i,cnper_bot,cnprim_bot,d
         
      if(norm_bot.lt.eps) then
c        the iterations were finished
         goto 20
      endif         

      cnper_p=cnper_bot+step
      cnper_m=cnper_bot-step
      d= dhot_sum_c(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     . vflow_ar,nll_in,cnper_p,cnprim_bot,iherm,K_sum)
      dp_r=dreal(d)
      dp_i=dimag(d)
      d= dhot_sum_c(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     . vflow_ar,nll_in,cnper_m,cnprim_bot,iherm,K_sum)
      dm_r=dreal(d)
      dm_i=dimag(d)

c      write(*,*)'i,cnper_p,cnper_m,d',i,cnper_p,cnper_m,d

      d_rd_d_rn=(dp_r-dm_r)/(2.d0*step)
      d_id_d_rn=(dp_i-dm_i)/(2.d0*step)
     
      cnprim_p=cnprim_bot+step
      cnprim_m=cnprim_bot-step
      
      d= dhot_sum_c(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     . vflow_ar,nll_in,cnper_bot,cnprim_p,iherm,K_sum)

c      write(*,*)'cnprim_p,cnprim_m,d',cnprim_p,cnprim_m,d

      dp_r=dreal(d)
      dp_i=dimag(d)
      d= dhot_sum_c(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     . vflow_ar,nll_in,cnper_bot,cnprim_m,iherm,K_sum)
      dm_r=dreal(d)
      dm_i=dimag(d)
      d_rd_d_in=(dp_r-dm_r)/(2.d0*step)
      d_id_d_in=(dp_i-dm_i)/(2.d0*step)

c     Newton method
c
c     d_rd_d_rn*(cnper_top-cnper_bot)+d_rd_d_in*(cnprim_top-cnprim_bot)=-d_r_bot
c     d_id_d_rn*(cnper_top-cnper_bot)+d_id_d_in*(cnprim_top-cnprim_bot)=-d_i_bot
c
      det=d_rd_d_rn*d_id_d_in-d_rd_d_in*d_id_d_rn
      det_r=-(d_r_bot*d_id_d_in-d_rd_d_in*d_i_bot)
      det_i=-(d_rd_d_rn*d_i_bot-d_id_d_rn*d_r_bot)

      alpha=1.d0
c      alpha=0.1d0
 30   continue
      cnper_top=cnper_bot+alpha*det_r/det
      cnprim_top=cnprim_bot+alpha*det_i/det

c      if (cnprim_top.lt.0.d0) then
c         alpha=-0.5*alpha
c         write(*,*)'cnprim_bot,alpha,cnprim_top',
c     .              cnprim_bot,alpha,cnprim_top
c         goto 30
c      endif

      d=dhot_sum_c(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     . vflow_ar,nll_in,cnper_top,cnprim_top,iherm,K_sum)
      d_r_top=dreal(d)
      d_i_top=dimag(d)
      norm_top=dsqrt(d_r_top**2+d_i_top**2)

c      write(*,*)'i,cnper_top,cnprim_top,d',i,cnper_top,cnprim_top,d
c      write(*,*)'norm_bot,norm_top,eps',norm_bot,norm_top,eps
     
      if (norm_top.gt.norm_bot) then
         if (alpha.gt.(eps*1.d2))then
            alpha=0.5*alpha
c           write(*,*)'alpha',alpha
            goto 30
         endif
      endif

      deltn=dsqrt((cnper_top-cnper_bot)**2+(cnprim_top-cnprim_bot)**2)
      cnper_bot=cnper_top
      cnprim_bot=cnprim_top
      norm_bot=norm_top
c      write(*,*)'forest.f i,deltn,norm_top',i,deltn,norm_top
      if (deltn.lt.eps1) then
c        iterations were finished
c         write(*,*)'deltn',deltn
         goto 20         
      endif

      i=i+1
  
      if (i.gt.iter_max)then
        write(*,*)'Newton hot N_perp i>iter_max',iter_max
        goto 20
      endif

      goto 10

 20   continue
      nperp_real=cnper_top
      nperp_imag=cnprim_top
      return
      end subroutine solv_nperp_hot


      double complex function  
     . dhot_sum_c(nbulk, mass_ar, x_ar, y_ar, t_av_ar, tpop_ar,
     . vflow_ar, nll_in, np_in, cnprim, iherm, K_sum)
c-----------------------------------------------------------
c     calculates the hot dielectric tensor  K_sum(3,3)
c     and the hot non-relativistic dispersion function dhotsum_c
c     for the hot electron and ion plasma
c     T.Stix, Waves in plasmas,(1992), p.252 (40), p.258 (57)
c     K_sum=1+sum{s}K_s  
c     It uses Hermitian(iherm=1) or full (iherm=2) dielectric tensor
c     and complex N_perp={np_in,cnprim}
c-----------------------------------------------------------
c     INPUTS:
c      nbulk - the total number of the plasma species 
c      mass_ar - the masses of the plasma  species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=av erage temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar = T_perp/T_parallel  
c      vflow_ar    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - Real part of the perpendicular index of refraction n.
c      cnprim  Imaginary part of N_perp
c      iherm =1 hermitian dielectric tensor, 2-full
c     OUTPUT:
c      K_sum(3,3) -  the nine components of the dielectric tensor
c                    evaluated at (mass,X,Y,Te,tpop,vflow,Nll,np)
c      dhotsum_c    -  hot double complex dispersion function
c-----------------------------------------------------------
      implicit none
c     input
      integer nbulk,iherm   
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double precision np_in,nll_in,cnprim
c-----output
      double complex K_sum(3,3)
c-----local
      double precision nll,nlls,np_r
      double complex np,nps
      double complex K(3,3) ! sucseptibilities 
      double complex Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz

      integer j,j1,j2
c-----external DHOT_s,zeroK


c-----initialization of the dielectric tensor K_sum(3,3)
      call zeroK(K_sum)

      do j=1,nbulk      
        call DHOT_s_c(mass_ar(j),x_ar(j),y_ar(j),t_av_ar(j),tpop_ar(j),
     .  vflow_ar(j),nll_in,np_in,cnprim, iherm,K)

        do j1=1,3
          do j2=1,3
            K_sum(j1,j2)=K_sum(j1,j2)+K(j1,j2)
          enddo
         enddo
       
      enddo

      K_sum(1,1) = K_sum(1,1) + dcmplx(1.0D0,0.0D0)
      K_sum(2,2) = K_sum(2,2) + dcmplx(1.0D0,0.0D0)
      K_sum(3,3) = K_sum(3,3) + dcmplx(1.0D0,0.0D0)
      
      Kxx = K_sum(1,1)
      Kyy = K_sum(2,2)
      Kzz = K_sum(3,3)
      Kxy = K_sum(1,2)
      Kxz = K_sum(1,3)
      Kyx = K_sum(2,1)
      Kyz = K_sum(2,3)
      Kzx = K_sum(3,1)
      Kzy = K_sum(3,2)

c-----hot dispersion function
      
!YuP      call npnllmin(nll_in,np_in,nll,np_r) ! YuP[2020-08-24] in dhot_sum_c: Having nll=0 is ok here
      !Note: np_r is not used.
      nll=nll_in !YuP[2020-08-24] Having nll=0 is ok here 
      np=dcmplx(np_in,cnprim)
      !np=dcmplx(max(0.d0,np_in),cnprim) !YuP[2020-08-24] Set lower limit on real part
      nlls=nll**2
      nps=np**2

      dhot_sum_c =(Kxx-nlls) * (Kyy-nlls-nps) * (Kzz-nps)
     .+  Kxy * Kyz * (Kzx+np*nll) 
     .+ (Kxz+np*nll) * Kyx * Kzy
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
     .-  Kzy * Kyz * (Kxx-nlls)
     .- (Kzz - nps) * Kyx * Kxy


      return
      end function dhot_sum_c






      subroutine hot_disp(comp_nperp,comp_d)
c     calculate the complex hot dispersion function comp_d
c     from the complex N_perp comp_nperp
c     the input plasma parameters are in commom /hot_plasma/
c     INPUTS:
c      nbulk - the total number of the plasma species 
c      mass_ar - the masses of the plasma  species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=av erage temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar = T_perp/T_parallel  
c      vflow_ar    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - Real part of the perpendicular index of refraction n.
c      cnprim  Imaginary part of N_perp
c      iherm =1 hermitian dielectric tensor, 2-full
c     OUTPUT:
c      K_sum(3,3) -  the nine components of the dielectric tensor
c                    evaluated at (mass,X,Y,Te,tpop,vflow,Nll,np)
c      comp_d     -  hot double complex dispersion function
c------------------------------------------
c-----input
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
cSm030514     
      double complex comp_nperp !complex N_perpendicular     
      integer nbulk_l,iherm_l
      double precision mass_ar_l(nbulka),x_ar_l(nbulka),y_ar_l(nbulka),
     .t_av_ar_l(nbulka),tpop_ar_l(nbulka),
     .vflow_ar_l(nbulka),nll_in_l
      double complex K_sum_l(3,3)
cSm030514
c      common /hot_plasma/nbulk_l,mass_ar_l,x_ar_l,y_ar_l,
c     .t_av_ar_l,tpop_ar_l,
c     .vflow_ar_l,nll_in_l,K_sum_l,iherm_l
      common /hot_plasma/mass_ar_l,x_ar_l,y_ar_l,
     .t_av_ar_l,tpop_ar_l,
     .vflow_ar_l,nll_in_l,K_sum_l,iherm_l,nbulk_l

c-----output
      double complex comp_d !complex dispersion function

c-----local
      double precision np_in,cnprim
     
c-----external 
      double complex dhot_sum_c
      cnprim=dimag(comp_nperp)   
      np_in=dreal(comp_nperp)
      comp_d=dhot_sum_c(nbulk_l,mass_ar_l,x_ar_l,y_ar_l,t_av_ar_l,
     .tpop_ar_l,vflow_ar_l,nll_in_l,np_in,cnprim,iherm_l,K_sum_l)
       
      return
      end subroutine hot_disp





      subroutine solv_nperp_hot_grad(nbulk,mass_ar,x_ar,y_ar,t_av_ar,
     . tpop_ar,vflow_ar,nll_in,cnprim_in,np_in,cnper_new,cnprim_new)
c-----solve hot plasma complex dispersion relation
c     D(Re(N_perp),Im(N_perp))=0 using the minimizasion of the function
c     J(Re(N_perp),Im(N_perp))=(ReD)**2+(ImD)**2 
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
         
c-----input
      double precision np_in,cnprim_in !the initial iteration for complex N_perp
      integer nbulk,iherm
      double precision  mass_ar(nbulka),x_ar(nbulka),y_ar(nbulka),
     .t_av_ar(nbulka),tpop_ar(nbulka),
     .vflow_ar(nbulka),nll_in
      double complex K_sum(3,3)
    
c-----output
      double precision cnper_new,cnprim_new !complex N_perp
    
c-----local
      integer  nbulk_l,iherm_l
      double precision mass_ar_l(nbulka),x_ar_l(nbulka),y_ar_l(nbulka),
     .t_av_ar_l(nbulka),tpop_ar_l(nbulka),
     .vflow_ar_l(nbulka),nll_in_l
      double complex K_sum_l(3,3)
cSm030514
c      common /hot_plasma/nbulk_l,mass_ar_l,x_ar_l,y_ar_l,
c     .t_av_ar_l,tpop_ar_l,
c     .vflow_ar_l,nll_in_l,K_sum_l,iherm_l
      common /hot_plasma/mass_ar_l,x_ar_l,y_ar_l,
     .t_av_ar_l,tpop_ar_l,
     .vflow_ar_l,nll_in_l,K_sum_l,iherm_l,nbulk_l

      double complex cz0 !initial iteration for complex N_perp 
      double complex comp_nperp !complex N_perpendicular
      integer i

c-----set the data in common block hot_plasma
      nbulk_l=nbulk
      do i=1,nbulk
         mass_ar_l(i)=mass_ar(i)
         x_ar_l(i)=x_ar(i)
         y_ar_l(i)=y_ar(i)
         t_av_ar_l(i)=t_av_ar(i)
         tpop_ar_l(i)=tpop_ar(i)
         vflow_ar_l(i)= vflow_ar(i)        
      enddo
      nll_in_l=nll_in !Nparallel

c-----solve the dispersion relation
      cz0=dcmplx(np_in,cnprim_in)  
c      call croot(cz0,comp_nperp)
      call grad_min(cz0,comp_nperp)
      cnper_new=dreal(comp_nperp)
      cnprim_new=dimag(comp_nperp)
       
      return
      end subroutine solv_nperp_hot_grad



      subroutine DHOT_s_c(mass,X,Y,T_av,tpop,vflow,
     .nll_in,np_in,cnprim,iherm,K)
c-----------------------------------------------------------
c     calculates sucseptibilities K(3,3) for the given hot (non-relativistic)
c     plasma specie, T.Stix, Waves in plasmas,(1992), p.258 (57)  make

c------------------------------------------------------------
c     INPUTS:
c      mass - the mass of the given specie (in electron mass) 
c      X = (fpe/f)**2
c      Y = fce/f ! for electron Ye should has the sign opposite to Yi
c      T_av=average temperature=Te(1+2*tpop)/3 in ev.
c           Here Te is the parallel temperature in eV  
c      tpop = T_perp/T_parallel  
c      vflow    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - Real part of the perpendicular refractition index =ReN_perp
c      cnprim - Image part of the perpendicular refractition index =ImN_perp
c      iherm =1 hermitian dielectric tensor, 2-full     
c     OUTPUT:
c      K(3,3):  the nine components of sucseptibilities
c               (for the given specie)
c          evaluated at (mass,X,Y,Te,tpop,vflow,Nll,np)
c-----------------------------------------------------------

      implicit none
c     input
      double precision mass
      double precision  X, Y,T_av,tpop,vflow
      double precision np_in,nll_in,cnprim
      integer iherm      
      integer n        !n.lt.nharmmax is the number of ECR harmonics
cSm060815
      integer n_new       
        
      integer nmax     ! max number of ECR harmonics in arrays
      parameter(nmax=300)


c     output
      double complex K(3,3)
c     locals
      double precision np,nll
      double precision nps,nlls,c
      double precision te, beta
      double precision pi,k4,vte,dlambda,dri01
cSm020410
      double complex lambda,cnps,cnp
      integer ize,ncalc

      double precision dr,di,dxne(nmax),drp
c      double precision Ibess(nmax),Ibessp(nmax),Ibesspp(nmax)
cSm020410
cSAP090204
c      double complex Ibess(nmax),Ibessp(nmax),Ibesspp(nmax)
      double complex Ibess(-(2*nmax+1):(2*nmax+1)),
     &Ibessp(-(2*nmax+1):(2*nmax+1)),Ibesspp(-(2*nmax+1):(2*nmax+1))
      double precision br(nmax),bi(nmax) !real and imaginary modified
                                         ! bessel func. I_n      
      double complex ri01,xne(nmax)

      double complex   dp
cSAP090204
c      double complex A(nmax),B(nmax)
      double complex A(-(2*nmax+1):(2*nmax+1)),
     &B(-(2*nmax+1):(2*nmax+1))
      double complex i
      double complex  Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz

      double complex cx,cz0,cz1,cd,cxp
      integer q,ierr,j

      integer n_harm0,n_harm1,n_harm2 ! the boundaries for the harmonic numbers
      integer istop   
c-----external
c     beslci1(x,y,nb,ize,br,bi,ncalc)
c     ctk--calculation of Bessel function with a complex argument (x,y)
c     br(nb),bi(nb)  
c     y=0,ize=0 >J_0(x)=br(1),J_1(x)=br(2)
c     y=0,ize=1 >I_0(x)=br(1),I_1(x)=br(2),... ,I_n(x)=br(n+1),
c     n=0,1,2,...; nb.ge.1
C     x,y,br,bi - 


c      istop=1, it is possible to calculate the susceptibility tensor
c      istop=0, it is impossible to calculate the susceptibilty tensor

  
      pi = 4*datan(1.d0)
      c= 2.99792458d10                !speed of light
      k4 = 4.19d7               !constant for elec therm vel
      i = ( 0.0d0,1.0d0)        !imaginary number

c-----tpop =  ratio of Te_perp / Te_parallel 
c-----vflow = flow velocity of electrons
      te=3.d0*T_av/(1.d0+2.d0*tpop) !longitudinal

c       write(*,*)'dhot_s_c longitudinal te,mass',te,mass      
ctest
 27   format(2D25.17)
c      write(*,*)'dhot_s_c tpop,vflow'
c      write(*,27)tpop,vflow
           
      vte =k4*dsqrt(2.0d0*te/mass)   !vte = sqrt(2.0*kTe/me) longitudinal
      beta = vte / c                 !longitudinal
     
c      write(*,*)'dhot vte, beta'
c      write(*,27)vte,beta
c      write(*,*)'dhot_s_c befor npnllmin nll_in,np_in',
c     .nll_in,np_in
      call npnllmin(nll_in,np_in,nll,np) !in subroutine DHOT_s_c

      nps = np**2
cSm020410
      cnp=dcmplx(np,cnprim)
      cnps=cnp*cnp

      nlls = nll**2
 
c-----lambda = 0.5 k_perp^2 * rho_e^2
cSm020410
      dlambda = nps*beta**2*tpop/2.0d0/Y**2
      lambda =cnps*beta**2*tpop/2.0d0/Y**2
c      write(*,*)'dhot_s_c cnps,beta,tpop,Y,lambda',cnps,beta,tpop,Y,lambda

 26   format(d25.17) 
c      write(*,*)'dhot_s_c nps,beta,tpop,Y,lambda',nps,beta,tpop,Y,lambda


c   Now, calculate the modified bessel functions 
c   Remember that these are actually exp(-lambda) * In(lambda)

cSAP09204
      do j = 1,nmax
         dxne(j)   = 0.0D0
c         ibess(j) = 0.0D0
c         ibessp(j)= 0.0D0
         Ibess(j)  = dcmplx(0.d0,0.d0)
         Ibessp(j) =dcmplx(0.d0,0.d0)
         Ibesspp(j)=dcmplx(0.d0,0.d0)
         br(j)=0.d0
         bi(j)=0.d0
         xne(j)=dcmplx(0.d0,0.d0)
         A(j) = dcmplx(0.0D0,0.0D0)
         B(j) = dcmplx(0.0D0,0.0D0)
      enddo
      do j = -(2*nmax+1),(2*nmax+1)
c         dxne(j)   = 0.0D0
c         ibess(j) = 0.0D0
c         ibessp(j)= 0.0D0
         Ibess(j)  = dcmplx(0.d0,0.d0)
         Ibessp(j) =dcmplx(0.d0,0.d0)
         Ibesspp(j)=dcmplx(0.d0,0.d0)
c         br(j)=0.d0
c         bi(j)=0.d0
c         xne(j)=dcmplx(0.d0,0.d0)
         A(j) = dcmplx(0.0D0,0.0D0)
         B(j) = dcmplx(0.0D0,0.0D0)
      enddo

c-----calulations of the max and minimal numbers of the harmonics
c      write(*,*)'in DHOT_s_c before harmon_z nmax',nmax     

cSm030515
      call harmon_z(X,Y,nll,vflow,c,beta,dlambda,nmax, ! dlambda=Real(lambda) 
     .n_harm1,n_harm2,n,istop)

c      write(*,*)'in DHOT_s_c 3843 after harmon_z mass,n_harm1,n_harm2,n
c     &,isto1p', mass,n_harm1,n_harm2,n,istop

      if (istop.eq.0) then
         write(*,*)'DHOT_s_c istop=0 it is impossible to calculate' 
         write(*,*)'the harmonics of susceptibily tensor n_harm0 nmax=',
     .   nmax
      endif      
       
      if (istop.eq.0) goto 10
      call ibess0(dlambda,dri01)        !dlambda=Re(lambda)
      call ibessn(dlambda,n+1,dri01,dxne)  !-> I_n(x)*exp(-x)  in DHOT_s_c
      !YuP[2020-01-23] The dri01 from the above two calls are not used anymore]
      !Why switching to beslci1() below? Because of complex lambda?

      ize=1 !for the modified bessel functions
      call beslci1(dreal(lambda),dimag(lambda),n+2,ize,br,bi,ncalc)
      ri01=dcmplx(br(1),bi(1))*cdexp(-lambda) !I_0*exp(-lambda)
      !YuP[2020-01-23] Potential underflow problem from exp(-lambda)

c      write(*,*)'n,lambda',n,lambda
      !YuP[2020-01-23] Better calculate exp(-lambda) once, outside of j-loop
      do j=1,n+1
        xne(j)=dcmplx(br(1+j),bi(1+j))*cdexp(-lambda) !I_j*exp(-lambda) j>0
c        write(*,*)'j,br(1+j),bi(1+j),cdexp(-lambda)',
c     &             j,br(1+j),bi(1+j),cdexp(-lambda) 
c        write(*,*)'n,j,xne(j)',n,j,xne(j)
      enddo

cSAP090306
cSAP090204 
c      n=n_new
      n_new=n
      do j=1,n

cSm020410 
c        write(*,*)'dhot_s_c j,xne(j)',j,xne(j)

        if (cdabs(xne(j)).lt.1.d-14) then
     
            
          call harmon_z(X,Y,nll,vflow,c,beta,dlambda,j+2,
     .    n_harm1,n_harm2,n_new,istop)

c          write(*,*)'dhot_s_c after harmon_z n,n_new',n_new
cSm060816
c          goto 10
          goto 9
        endif
      enddo

cSm060816
 9    continue
      n=n_new
      
 10   continue
c      write(*,*)'dhot_s_c j, xne(j)',j,xne(j)
c       write(*,*)'dhot_s_c n,n_harm1,n_harm_2',n,n_harm1,n_harm2
c      do q  = -n,n
c      write(*,*)'dhot_s_c before goto 20 istop=',istop     
      if (istop.eq.0) goto 20

      do q  = n_harm1,n_harm2 

         cx = (1.d0 - nll*vflow/c - q *  Y) / nll / beta

c         write(*,*)'DHOT_s_c cx q',q 
c         write(*,27)dreal(cx),dimag(cx)
           
         call CZETA0(nll,CX,CZ0,CZ1,IERR)
c         write(*,*)'DHOT_s_c q,IERR,cz0',q,IERR,cz0      

c---- Use the Stix notation of A(n) and B(n)
c
c---- My functions A(n) = omega * Astix(n)
c           and    B(n) = omega * Bstix(n) / c

       if(dabs(dfloat(q+n+1)).gt.(2*nmax+1)) then
         write(*,*)'Dhot_s_c before A: q,n,q+n+1',q,n,q+n+1
       endif

        A(q+n+1)  =  (tpop - 1.d0) +  
     .   ((1.0D0-nll*vflow/c-q* Y) * tpop + q*  Y ) * cz0 /nll / beta

         B(q+n+1)  =  (
     .                 (1.D0 - q*Y) *A(q+n+1)
     .                 +(1.D0 - nll*vflow/c)   
     .                )/nll
                  
         if (q.eq.0)  Ibess(q+n+1) = ri01
         if (q.ne.0)  Ibess(q+n+1) = xne(Iabs(q))
     

         if (q.eq.0)  Ibessp(q+n+1) = xne(1) 
         if (q.ne.0)  Ibessp(q+n+1)=xne(1+iabs(q))
     .                 +Iabs(q)* xne(iabs(q)) / lambda
     .                          
         if (q.eq.0)  Ibesspp(q+n+1) = xne(2)+xne(1)/lambda 
         if (q.ne.0)  Ibesspp(q+n+1)=xne(2+iabs(q)) 
     .     +(2.0D0*Iabs(q)+ 1)*xne(iabs(q)+1 )/lambda
     .     +(Iabs(q)*Iabs(q)/lambda**2 - Iabs(q)/lambda**2)*xne(iabs(q))

      enddo
 20   continue

       Kxx = (0.0D0,0.0D0)
       Kyy = (0.0D0,0.0D0)
       Kzz = (0.0D0,0.0D0)
       Kzz = dcmplx(2.d0 * X * vflow/c / nll / (beta**2 * tpop),0.0D0)
       Kxy = (0.0D0,0.0D0)
       Kxz = (0.0D0,0.0D0)
       Kyz = (0.0D0,0.0D0)


c      do q = -n,n

      if (istop.eq.0) goto 30 
      do q  = n_harm1,n_harm2
        
        Kxx=  Kxx+ X * q**2.0D0 * Ibess(q+n+1) * A(q+n+1) / lambda

        Kxy=  Kxy-i* q *  X *(Ibess(q+n+1)-Ibessp(q+n+1))*A(q+n+1)


        Kyy=  Kyy+ X*(q**2*Ibess(q+n+1)/lambda +
     .  2.0D0*lambda*Ibess(q+n+1)
     *     - 2.d0*lambda*Ibessp(q+n+1))* A(q+n+1)
        Kxz=  Kxz + X* np*q*Ibess(q+n+1) * B(q+n+1) / lambda /  Y
        Kyz= Kyz+i* X*np*(Ibess(q+n+1)-Ibessp(q+n+1))*B(q+n+1)/ Y
        Kzz=  Kzz+ 2.0D0* X*(1.D0 - q* Y)*Ibess(q+n+1)*B(q+n+1)
     *      /nll/beta**2

      enddo
 30   continue

cSm991110
c      Kxx = Kxx + dcmplx(1.0D0,0.0D0)
c      Kyy = Kyy + dcmplx(1.0D0,0.0D0)
c      Kzz = Kzz + dcmplx(1.0D0,0.0D0)

      Kyx = -1.0D0* Kxy
      Kzx =  Kxz
      Kzy = -1.0D0* Kyz

      K(1,1) =  Kxx
      K(2,2) =  Kyy
      K(3,3) =  Kzz
      K(1,2) =  Kxy
      K(1,3) =  Kxz
      K(2,1) =  Kyx
      K(2,3) =  Kyz
      K(3,1) =  Kzx
      K(3,2) =  Kzy

c     Hermitian part
      if (iherm.eq.1) then
       Kxx = 0.5D0 * ( K(1,1) + dconjg(K(1,1) )) !kxx
       Kyy = 0.5D0 * ( K(2,2) + dconjg(K(2,2) )) !kyy
       Kzz = 0.5D0 * ( K(3,3) + dconjg(K(3,3) )) !kxx
       Kxy = 0.5D0 * ( K(1,2) + dconjg(K(2,1) )) !kxy
       Kxz = 0.5D0 * ( K(1,3) + dconjg(K(3,1) )) !kxz
       Kyz = 0.5D0 * ( K(2,3) + dconjg(K(3,2) )) !kxz
       Kyx = 0.5D0 * ( K(2,1) + dconjg(K(1,2) )) !kyx
       Kzx = 0.5D0 * ( K(3,1) + dconjg(K(1,3) )) !kzx
       Kzy = 0.5D0 * ( K(3,2) + dconjg(K(2,3) )) !kxz
      endif

      K(1,1) =  Kxx
      K(2,2) =  Kyy
      K(3,3) =  Kzz
      K(1,2) =  Kxy
      K(1,3) =  Kxz
      K(2,1) =  Kyx
      K(2,3) =  Kyz
      K(3,1) =  Kzx
      K(3,2) =  Kzy

      return
      end subroutine DHOT_s_c

      



cZ_function from Torray
       subroutine zfun (z,fu)
cProlog

      implicit none
c
c explicit type declaration 7/12/01 (RAJ)
      !doubleprecision to real*8 !YuP[2020-01-27] 
      real*8 delta, yi, osqpi, tpi, w, conoi, x, y, e, f, c, d,
     & tpiod, g, yy, h, za, zr, zi, a, b, den, oden
      integer n, nm3, i, itest, no2, no2p1, npi, nmi
c
c [1] (4550)zfun,     1-aug-83,  edit by a.kluge(user 4550 @dma)
c [1] release plasma dispersion function for use with toray
c
c *********************
c       routine which evaluates the plasma dispersion function.  uses
c       numerical integration (absolute value of the complex argument, z,
c       less than 5) or asymptotic expansion.
c *********************
c
      !YuP[2020-01-27] doublecomplex z, fu, temp1, temp2, z2, tpiiod
      complex*16 z, fu, temp1, temp2, z2, tpiiod !YuP[2020-01-27]
      
c
      dimension c(21), d(21), e(21), f(21), w(21)
c
      data delta/0.5d0/, yi/-1.0d0/, n/21/, itest/0/, osqpi/0.
     &56418958355d0/, tpi/6.28318530718d0/
c
      if (itest .eq. 1) go to 1
c
      itest = 1
c
c *********************
c       define weights and store constants used for integration.
c       weights are derived from a 3-point integration scheme
c *********************
c
      nm3 = n - 3
      conoi = delta * osqpi
      w(1) = conoi * 0.375d0
      w(2) = conoi * 1.1666666667d0
      w(3) = conoi * 0.95833333333d0
c
      do 20 i = 1,3
c
   20 w(n-i+1) = w(i)
c
      do 30 i = 4,nm3
c
   30 w(i) = conoi
c
      no2 = n/2
c
      no2p1 = no2 + 1
c
      x = 0.0d0
c
      y = yi
c
      e(no2p1) = x
c
      f(no2p1) = y
c
      temp1 = dcmplx(x,y)
c
      temp2 = EXP (-temp1 * temp1)
c
      c(no2p1) = dble(temp2) * w(no2p1)
c
      d(no2p1) = dimag(temp2) * w(no2p1)
c
      do 200 i = 1,no2
      x = delta * FLOAT (i)
      npi = no2p1 + i
      nmi = no2p1 - i
      temp1 = dcmplx(x,y)
      temp2 = EXP (-temp1 * temp1)
      c(npi) = dble(temp2) * w(npi)
      c(nmi) = dble(temp2) * w(nmi)
      d(npi) = dimag(temp2) * w(npi)
      d(nmi) = -dimag(temp2) * w(nmi)
      e(nmi) = -x
      e(npi) = x
      f(npi) = y
      f(nmi) = y
  200 continue
c
      tpiod = tpi/delta
      tpiiod = dcmplx(0.0d0,tpiod)
c
c *********************
c       begin calculations
c *********************
c
    1 g = dble(z)
      yy = dimag(z)
      h = ABS (yy)
      za = z * dconjg(z)
      if (za .ge. 25.0d0) go to 5
      z2 = z * z
c
c *********************
c       numerical integration.
c       f = 1/sqrt (pi)*sum of...w(i)*EXP (-x(i)**2)/(x(i)-z)...i = 1,n.
c       integration is along a line x(i) in the complex plane, where
c       the imaginary part of x(i) = yi and the difference between
c       successive real parts of x(i) = delta.  limits of integration
c       are from -delta*n/2 to delta*n/2.
c       compute the integral by taking the sum from 1 to n of the
c       constants divided by x(i)-z.  uses real arithmetic.
c *********************
c
      zr = 0.0d0
      zi = 0.0d0
c
      do 7 i = 1,n
      a = e(i) - g
      b = f(i) - h
      den = a * a + b * b
      oden = 1.0d0/den
      zr = zr + (a * c(i) + b * d(i)) * oden
      zi = zi + (a * d(i) - b * c(i)) * oden
    7 continue
c
c *********************
c       add the correction term
c *********************
c
      fu = dcmplx(zr,zi) + (0.0d0,-3.5449077018d0)* EXP (-z2 - tpiod*(h 
     &- yi) + tpiiod * g)
c
      if (yy .ge. 0.0d0) go to 6
c
c *********************
c       imaginary part of argument is negative.
c *********************
c
      fu = dconjg(fu) + (0.0d0,3.5449077018d0) * EXP (-z2)
c
      go to 6
c
c *********************
c     magnitude of argument is greater than 5, use asymptotic expansion
c *********************
c
    5 call aexpan (z,g,h,yy,fu)
    6 return
c
      end subroutine zfun


      subroutine aexpan (z, g, h, yy, fu)
cProlog

      implicit none
c
c **********************
c routine which computes the plasma dispersion function using
c asymptotic expansion.  if the imaginary part of the argument, yy,
c is equal to zero real arithmetic is used
c **********************
c
      !YuP[2020-01-27] doublecomplex z, fu, a, z2, oz2
      complex*16 z, fu, a, z2, oz2 !YuP[2020-01-27]
c
c  explicit type declaration added 7/12/01  RAJ
      !YuP[2020-01-27] doubleprecision yy, en, h, g, x2, ox2, f, b, c
      real*8 yy, en, h, g, x2, ox2, f, b, c !YuP[2020-01-27] 
      integer n,i
      data n/8/
c
      if (yy .eq. 0.0d0) go to 10
c
c **********************
c     complex arithmetic
c **********************
c
      z2 = z * z
      oz2 = 1.0d0/z2
      fu = -1.0d0/z
      a = fu
      en = 0.5d0
c
      do i=1,n
      a = en * a * oz2
      fu = fu + a
      en = en + 1.0d0
      enddo
c
      if (yy .gt. 0.0d0) go to 30
      if (h .gt. sqrt (g * g + 172.0d0)) go to 20
      fu = fu + (0.0d0,3.5449077018d0) * EXP (-z2)
      go to 30
c
c **********************
c       real arithmetic. error stop to avoid overflow
c **********************
c
   20 call STOP ('subroutine AEXPAN: impending floating overflow', 1)
c
c **********************
c       real arithmetic
c **********************
c
   10 x2 = g * g
      ox2 = 1.0d0/x2
      f = -1.0d0/g
      b = f
      en = 0.5d0
c
      do i=1,n
      b = en * b * ox2
      f = f + b
      en = en + 1.0d0
      enddo
c
      c = 1.7724538509d0 * EXP (-x2)
      fu = dcmplx(f, c)
   30 return
c
      end subroutine aexpan
      

      double precision function fomega_hot(cnz,cnr,cm,z,r,phi)
c-----It calculates the dispersion function and dielectric tensor reps 
c     for non-relativistic hot plasma 
c     with full (Hermitian + anti-Hermitian) dielctric tensor using
c     Berstein-Friedland formula.
c     The hot plasma complex dielectric tensor reps(3,3) will be in eps.i
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      integer i
      real*8 r,z,phi,hf,hw,df,step,t_kev,fplus,hfrqnc
      real*8 frqncpl,frqncmn,fminus,cnpar,cnper,cn2
      real*8 hamiltmp,hamiltmq
      real*8 cnr,cnrplus,cnz,cnzplus,cm,cmplus
      real*8 cnparplus,cnplus2,cnperplus,cnrminus,cnzminus
      real*8 cnparminus,cnminus2,cmminus,cnperminus
      include 'param.i'
      include 'one.i'
      include 'ions.i'	
      include 'eps.i'   
c-----external b
      double precision b,x,y,tempe,tpoprho,vflowrho
c-----------------------
     
      double precision wp(nbulka),vp(nbulka)

      double precision x_ar_plus(nbulka),y_ar_plus(nbulka),
     &x_ar_minus(nbulka),y_ar_minus(nbulka),
     &t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka)

      double complex K_sum(3,3),dhot_sum,dplus,dminus  

cSm030426
c----for the dielectric tensor calculations      
      double precision x_ar(nbulka),y_ar(nbulka)
      double complex d

      step=1.d-5
      hw=step*1.0d0
      hfrqnc=hw
      hfrqnc=hfrqnc*frqncy
      bmod=b(z,r,phi)
      do 11 i=1,nbulk
        vp(i)=v(i)
        wp(i)=w(i)
 11   continue
      frqncpl=frqncy+hfrqnc
      df=frqncy/frqncpl
      do 12 i=1,nbulk
        v(i)=vp(i)*df* df
        w(i)=wp(i)*df
 12   continue

      do i=1,nbulk                  
         x_ar_plus(i)=x(z,r,phi,i)
         y_ar_plus(i)=y(z,r,phi,i)
         if(i.eq.1) y_ar_plus(1)=-y_ar_plus(1)
         t_keV=tempe(z,r,phi,i)   !keV averaged temperature
         t_av_ar(i)=t_keV*1.d3    !eV
         tpop_ar(i)=tpoprho(rho,i)
	 vflow_ar(i)=vflowrho(rho,i) !cm/sec
      enddo
c************************************************
      cnrplus=cnr*df
      cnzplus=cnz*df
      cmplus=cm*df

      cnparplus=(cnrplus*br+cnzplus*bz+(cmplus/r)*bphi)/bmod
      cnplus2=cnrplus**2+cnzplus**2+(cmplus/r)**2
      cnperplus=dsqrt(cnplus2-cnparplus**2)
c      write(*,*)'forest fomega before dplus'
c      write(*,*)'nbulk, cnparplus,cnperplus',nbulk,cnparplus,cnperplus
c      do i=1,nbulk
c        write(*,*)'i,dmas(i),x_ar_plus(i),y_ar_plus(i)',
c     &  i,dmas(i),x_ar_plus(i),y_ar_plus(i)
c        write(*,*)'t_av_ar(i),tpop_ar(i),vflow_ar(i)',
c     &  t_av_ar(i),tpop_ar(i),vflow_ar(i)
c      enddo

      dplus= dhot_sum(nbulk,dmas,x_ar_plus,y_ar_plus,t_av_ar,tpop_ar,
     &vflow_ar,cnparplus,cnperplus,2,K_sum)
c      write(*,*)'dplus',dplus 
      hamiltmp=dimag(dplus)
      hamiltmq=dreal(dplus)
      fplus=hamiltmp*hamiltmp+hamiltmq*hamiltmq
c*************************************************
c----------------------------------------------------------
      frqncmn=frqncy-hfrqnc
      df=frqncy/frqncmn
      do 15 i=1,nbulk
        v(i)=vp(i)*df*df
        w(i)=wp(i)*df
 15   continue

      do i=1,nbulk                  
         x_ar_minus(i)=x(z,r,phi,i)
         y_ar_minus(i)=y(z,r,phi,i)
         if(i.eq.1) y_ar_minus(1)=-y_ar_minus(1)
      enddo
c************************************************
      cnrminus=cnr*df
      cnzminus=cnz*df
      cmminus=cm*df

      cnparminus=(cnrminus*br+cnzminus*bz+(cmminus/r)*bphi)/bmod
      cnminus2=cnrminus**2+cnzminus**2+(cmminus/r)**2
      cnperminus=dsqrt(cnminus2-cnparminus**2)
c      write(*,*)'forest fomega before dpminus'
c      write(*,*)'cnparminus,cnperminus',cnparminus,cnperminus
c      do i=1,nbulk
c        write(*,*)'i,dmas(i),x_ar_minus(i),y_ar_minus(i)',
c     &  i,dmas(i),x_ar_minus(i),y_ar_minus(i)
c        write(*,*)'t_av_ar(i),tpop_ar(i),vflow_ar(i)',
c     &  t_av_ar(i),tpop_ar(i),vflow_ar(i)
c      enddo
      dminus= dhot_sum(nbulk,dmas,x_ar_minus,y_ar_minus,t_av_ar,
     &tpop_ar,
     &vflow_ar,cnparminus,cnperminus,2,K_sum)
c      write(*,*)'dminus',dminus 
      hamiltmp=dimag(dminus)
      hamiltmq=dreal(dminus)
      fminus=hamiltmp*hamiltmp+hamiltmq*hamiltmq
c*************************************************
      fomega_hot=(fplus-fminus)/(2.0d0*hw)
c      write(*,*)'forest.f fomega_hot',fomega_hot 
c-----------------------------------------------------------
      do 14 i=1,nbulk
        v(i)=vp(i)
        w(i)=wp(i)
 14   continue


cSm030426 
c-----calculations of the hot plasma complex dielectric tensor reps
c-----It will be in 'eps.i'
      do i=1,nbulk                  
         x_ar(i)=x(z,r,phi,i)
         y_ar(i)=y(z,r,phi,i)
         if(i.eq.1) y_ar(1)=-y_ar(1)
         t_keV=tempe(z,r,phi,i)   !keV averaged temperature
         t_av_ar(i)=t_keV*1.d3    !eV
         tpop_ar(i)=tpoprho(rho,i)
	 vflow_ar(i)=vflowrho(rho,i) !cm/sec
      enddo
 
      cnpar=(cnr*br+cnz*bz+(cm/r)*bphi)/bmod
      cn2=cnr**2+cnz**2+(cm/r)**2
      cnper=dsqrt(cn2-cnpar**2)

      d= dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &vflow_ar,cnpar,cnper,2,reps)
     
      return
      end function fomega_hot




 
      subroutine hot_roots_solver(nbulk,T_ar,tpop_ar,vflow_ar,
     &X_ar,Y_ar,cnpar,
     &N_perp_root_max,n_points,  
     &N_perp_root_ar,n_hot_roots,
     &e_xyz_nperp_root_ar,e_mode_root_ar)
c-------------------------------------------------------------------------------
c     calculates roots N_perp_ar(k=1,..,n_hot_roots),
c     of the hermitian hot plasma dispersion function
c
c     D_hot(T_ar,X_ar,Y_ar,N_par,N_perp_root)=0
c      
c     in the interval 0 < N_perp_root< N_perp_root_max.
c
c     To find roots the interval (0,N_perp_root_max) will be divided by 
c     n_points equal bins 
c     N_perp_bin(i)=[N_perp_root_max/(n_points-1)]*(i-1), i=1,..., n_points.
c     
c     The root will be calculated in each bin, where the dispersion function
c     has different signs at bins boundaries
c
c     Electric filed polarisation e_nperp_root_ar(1-3,k) will be calculated 
c     for each root N_perp_root_ar(k)
c     e_xyz_nperp_root_ar(1,k) E_x compomnent
c     e_xyz_nperp_root_ar(2,k) E_y compomnent
c     e_xyz_nperp_root_ar(3,k) E_z compomnent
c     e_mode_root_ar(1,k) E_plus=cdabs(cex+i*cey)/sqrt(2)
c     e_mode_root_ar(2,k) E_minus=cdabs(cex-i*cey)/sqrt(2)
c     e_mode_root_ar(3,k) E_parallel=cdabs(cez*cnpar+cex*dreal(cnx))/
c                                    dsqrt(cnpar**2+dreal(cnx)**2)
c------------------------------------------------------------------------------
      implicit none 
      
      include 'param.i'
      include 'ions.i'
      include 'nperpcom.i'
      include 'eps.i'
      include 'cefield.i'

c-----input
      integer
     &nbulk,             ! number of plasma components
     &n_points           !the number of equal bins at the interval
                         !(0,N_perp_root_max)   

      double precision
     & cnpar,            ! parallel refractive index
     & T_ar(nbulka),     ! temperature KEV for each plasma specie
     & tpop_ar(nbulka),  ! T_perp/T_parallel for each plasma specie
     & vflow_ar(nbulka), ! drift velocity (cm/sec) for each plasma specie
     & X_ar(nbulka),     ! (omega_p/omega)**2  for each plasma specie
     & Y_ar(nbulka),     ! omega_c/omega  the algebraic variable for each plasma specie
     & N_perp_root_max   ! maximal value of the perpendicular refractive index
                         ! used to find roots       
c-----output
      integer 
     & n_hot_roots       !The total number of calculated roots

      double precision  
     & N_perp_root_ar(n_hot_roots_a) ! roots of D_hot(N_perpendicular)=0

      double complex e_xyz_nperp_root_ar(3,n_hot_roots_a) !electric field
                                                          !polarization
                                                          !Ex,Ey,Ez
                                                          !at each N_perp_root(i)
      double precision  e_mode_root_ar(3,n_hot_roots_a)  !electric field
                                                         !polarization at
                                                         ! E+,E-,E_parallel
                                                         !each N_perp_root(i)

c--------------------------------------------------------
c                                                 e_x(i)=e_nperp_root(1,i)
c                                                 e_y(i)=e_nperp_root(2,i)
c                                                 e_z(i)=e_nperp_root(3,i)
c-----------------------------------------------------------

c-----external 
      double precision rtbis
      double complex dhot_sum,d_hot_nper
      external d_hot_nper
c-----locals
      integer i,k1,k2
      double precision cnper,step,dleft,dright,xacc,x1,x2,
     *ex,ey,ez,eplus,eminus,epar
      double complex K(3,3),cnx



c-----set data to  common/nperpcom/
      nllc=cnpar
      do i=1,nbulk
         tec(i)= T_ar(i)
         tpopc(i)= tpop_ar(i)
         vflowc(i)=vflow_ar(i)
         xc(i)=X_ar(i)
         yc(i)=Y_ar(i)
      enddo
c--------------------------------------

      xacc=1.d-14  ! accuracy of the root calculations
      
      n_hot_roots=0
      step=N_perp_root_max/dfloat(n_points-1)
 
      write(*,*)'hot_roots_solver: n_points,step',n_points,step
     
      cnper=0.d0
      dleft= dreal(dhot_sum(nbulk,dmas,x_ar,y_ar,t_ar,tpop_ar,
     .      vflow_ar,cnpar,cnper,1,K))

      write(*,*)'hot_roots_solver: nbulk,cnpar',nbulk,cnpar
      do i=1,nbulk
        write(*,*)'=i,dmas(i),x_ar(i),y_ar(i)',i,dmas(i),x_ar(i),y_ar(i)
        write(*,*)'t_ar(i),tpop_ar(i),vflow_ar(i)',
     &             t_ar(i),tpop_ar(i),vflow_ar(i)
      enddo

      do i=1,n_points
        cnper=step*(i-1)
        dright=dreal(dhot_sum(nbulk,dmas,x_ar,y_ar,t_ar,tpop_ar,
     .      vflow_ar,cnpar,cnper,1,K))
c         write(*,*)'forest.f i,cnper,dright',i,cnper,dright
c        write(*,*)'forest.f i,cnper,dleft,dright',i,cnper,dleft,dright

        if ((dleft*dright).lt.0.d0) then
           n_hot_roots=n_hot_roots+1
           if(n_hot_roots.gt.n_hot_roots_a) then
              write(*,*)'n_hot_roots.gt.n_hot_rots_a'
              write(*,*)'n_hot_roots,n_hot_roots_a',
     &                   n_hot_roots,n_hot_roots_a
              write(*,*)'Please increase parameter n_hot_rots_a'
              write(*,*)'in param.i'
              write(*,*)'and recompile the code'
              stop 'forest,f insubroutine hot_roots_solver'
           endif  
           x1= cnper-step
           x2= cnper

           N_perp_root_ar(n_hot_roots)=rtbis(d_hot_nper,x1,x2,xacc)

c------------------------------------------------------------------
c          electric field for hot plasma tensor
c------------------------------------------------------------------
c          cnx=dcmplx(cnper,cnprim)
           cnx=dcmplx( N_perp_root_ar(n_hot_roots),0.d0)

           do k1=1,3
             do k2=1,3
              reps(k1,k2)=K(k1,k2) ! put hot tensor to common /eps.i/
             enddo
           enddo
           call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
           e_xyz_nperp_root_ar(1,n_hot_roots)=cex
           e_xyz_nperp_root_ar(2,n_hot_roots)=cey
           e_xyz_nperp_root_ar(3,n_hot_roots)=cez
           e_mode_root_ar(1,n_hot_roots)=eplus
           e_mode_root_ar(2,n_hot_roots)=eminus
           e_mode_root_ar(3,n_hot_roots)=epar

c           write(*,*)'cex',cex
c           write(*,*)'cey',cey
c           write(*,*)'cez',cez
c           write(*,*)'eplus',eplus
c           write(*,*)'eminus',eminus
c           write(*,*)'epar',epar

        endif
        dleft=dright
     

      enddo

      return
      end subroutine hot_roots_solver

      
      double precision function d_hot_nper(cnper)
c-----calculates Real(hot dispersion function)
c     INPUT:
c       cnper=Re(N_perpendicular))
c       the input data from common /nperpcom.i/
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'nperpcom.i'
c-----input
      double precision cnper
c-----local
      double complex  K(3,3),d
     
c-----external: dhot_sum
      double complex dhot_sum
   

      d=dhot_sum(nbulk,dmas,xc,yc,tec,tpopc,
     .vflowc,nllc,cnper,1,K)
      d_hot_nper = dreal(d) 

      return
      end function d_hot_nper


      subroutine calculate_hot_nperp_roots(z,r,phi,cnpar,
     &n_hot_roots,N_perp_root_ar)
c---------------------------------------------------------------------
c     calculates hot plasma dispersion function roots
c     N_perp and polarization at given point z,r,phi
c     and cnpar=N_parallel
c----------------------------------------------------------
      implicit none
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'      
c-----input
      double precision
     &z,r,phi,           !space coordinates
     &cnpar              !N_parallel
         
c-----output
      integer n_hot_roots              !number of hot roots
      double precision 
     &N_perp_root_ar(n_hot_roots_a)    !hot plasma roots
      double complex e_xyz_nperp_root_ar(3,n_hot_roots_a) ! electric filed 
                                                      ! for ecach root
                               ! e_nperp_root_ar(1,k) E_x compomnent
                               ! e_nperp_root_ar(2,k) E_y compomnent
                               ! e_nperp_root_ar(3,k) E_z compomnent

      double precision  e_mode_root_ar(3,n_hot_roots_a)  !electric field
                                                         !polarization 
                                                         ! E+,E-,E_parallel
                                                         !at each N_perp_root(i)]
                               !e_mode_root_ar(1,k)      ! E+
                               !e_mode_root_ar(2,k)      ! E-
                               !e_mode_root_ar(3,k)      ! E||


c-----externals
      double precision b,x,y,tempe,tpoprho,vflowrho
c-----locals
      integer j
      double precision x_ar(nbulka),y_ar(nbulka),T_ar_ev(nbulka),     
     &tpop_ar(nbulka),vflow_ar(nbulka)![sm/sec]

      character*60  text
      character*6 format
      character*2 format1

      bmod=b(z,r,phi)
      do j=1,nbulk
        x_ar(j)=x(z,r,phi,j)
        y_ar(j)=y(z,r,phi,j)
        if(j.eq.1) y_ar(1)=-y_ar(1) ! negative Y=(omega_ce/omega)
                                    ! for electrons
        T_ar_ev(j)=tempe(z,r,phi,j)*1.d+3 !(eV) averaged temperature
        tpop_ar(j)=tpoprho(rho,j)
        vflow_ar(j)=vflowrho(rho,j)         
      enddo

      call hot_roots_solver(nbulk,T_ar_ev,tpop_ar,vflow_ar,
     &x_ar,y_ar,cnpar,
     &cN_perp_root_max,n_points_root,  
     &N_perp_root_ar,n_hot_roots,
     &e_xyz_nperp_root_ar,e_mode_root_ar)

      if(myrank.ne.0) return
      
      write(*,*)'Number of roots:n_hot_roots',n_hot_roots 
      open(10,file='hot_roots_at_given_point.dat')
      format='d21.15'
      format1='i2'
      text='(1X,"r=",'//format//')'
      write(10,fmt=text)r
      text='(1X,"z=",'//format//')'
      write(10,fmt=text)z
      text='(1X,"phi=",'//format//')'
      write(10,fmt=text)phi
      text='(1X,"cnpar=",'//format//')'
      write(10,fmt=text)cnpar
      text='(1X,"The code searched roots at the interval")'
      write(10,fmt=text)
      text='(1X,"0 < N_perp < cN_perp_root_max ")'
      write(10,fmt=text)
      text='(1X,"using n_points_root uniform N_perp mesh")'
      write(10,fmt=text)
      text='(1X,"cN_perp_root_max=",'//format//')'
      write(10,fmt=text)cN_perp_root_max
      text='(1X," n_points_root=",'//format1//')'
      write(10,fmt=text)n_points_root

      do j=1,n_hot_roots
        write(*,*)'===============hot root================='
        write(*,*)'j,N_perp_root_ar(j)',j,N_perp_root_ar(j)
        write(*,*)'EX e_xyz_nperp_root_ar(1,j)',e_xyz_nperp_root_ar(1,j)
        write(*,*)'EY e_xyz_nperp_root_ar(2,j)',e_xyz_nperp_root_ar(2,j)
        write(*,*)'EZ e_xyz_nperp_root_ar(3,j)',e_xyz_nperp_root_ar(3,j)
        write(*,*)'E+ e_mode_root_ar(1,j)',e_mode_root_ar(1,j)
        write(*,*)'E- e_mode_root_ar(2,j)',e_mode_root_ar(2,j)
        write(*,*)'E|| e_mode_root_ar(3,j)',e_mode_root_ar(3,j)
        write(*,*)'========================================'
c-----------------------------------------------------------------------
c       write roots and polarization to the file:
c       hot_roots_at_given_point.dat
c------------------------------------------------------------------
        text='(1X,"========================================")'
        write(10,fmt=text)
        text='(1X,"root number=",'//format1//')'
        write(10,fmt=text)j
        text='(1X,"root N_perpendicular=",'//format//')'
        write(10,fmt=text)N_perp_root_ar(j)
        text='(1X,"ReEx=",'//format//')'
        write(10,fmt=text) dreal(e_xyz_nperp_root_ar(1,j))
        text='(1X,"ImEx=",'//format//')'
        write(10,fmt=text) dimag(e_xyz_nperp_root_ar(1,j))
        text='(1X,"ReEy=",'//format//')'
        write(10,fmt=text) dreal(e_xyz_nperp_root_ar(2,j))
        text='(1X,"ImEy=",'//format//')'
        write(10,fmt=text) dimag(e_xyz_nperp_root_ar(2,j))
        text='(1X,"ReEz=",'//format//')'
        write(10,fmt=text) dreal(e_xyz_nperp_root_ar(3,j))
        text='(1X,"ImEz=",'//format//')'
        write(10,fmt=text) dimag(e_xyz_nperp_root_ar(3,j))
        text='(1X,"|E+|=",'//format//')'
        write(10,fmt=text)e_mode_root_ar(1,j)
        text='(1X,"|E-|=",'//format//')'
        write(10,fmt=text)e_mode_root_ar(2,j)
        text='(1X,"|E|||=",'//format//')'
        write(10,fmt=text)e_mode_root_ar(3,j)

      enddo

      close(10)

      return
      end subroutine calculate_hot_nperp_roots


      subroutine rho_ini_hot_nperp_roots(r_edge,z_edge,phi_edge,cnpar)     
c     &rho_ini,z_ini,r_ini)
c-----finds the small radius rho_ini > rho_min_find_hot_nperp_roots
c     at the vector rho^ where
c     hot plasma dispersdion function D_hot(npar) has three roots.
c     The vector rho^ is starting at the edge point (r_edge,z_edge,phi_edge),
c     and directed to the magnetic axis O(rma,zma,phi_edge)
c      
      implicit none
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'three.i'
c-----input
      double precision
     &r_edge,z_edge,phi_edge,          ! edge point coordinates
     &cnpar                            ! N_parallel
     
c-----output
      double precision  rho_ini(n_hot_roots_a), ! normalized small radius of M point 
                                                ! where D_hot(Nperp)=0 has 'i' roots
     &z_ini(n_hot_roots_a),r_ini(n_hot_roots_a),
      !space coordinates of M point
 
     &N_perp_roots_diff_points(n_hot_roots_a,n_hot_roots_a) 
      ! hot plasma roots in n_hot_roots_a points
      
      double complex e_xyz_diff_points(3,n_hot_roots_a,n_hot_roots_a) 
                     ! electric filed Ex,Ey,Ez
                     ! for each root at i=1,... n_hot_roots_a points
                     ! e_xyz_diff_points(1,k,i) = E_x compomnent
                     ! e_xyz_diff_points(2,k,i) = E_y compomnent
                     ! e_xyz_diff_points(3,k,i) = E_z compomnent

      double precision e_mode_diff_points(3,n_hot_roots_a,n_hot_roots_a)  
                      !electric field polarization E+,E-,E_parallel
                      !at each N_perp_root(k)  an i=1,... n_hot_roots_a points
                      !e_mode_diff_points(1,k,i) = E+
                      !e_mode_diff_points(2,k,i) = E-
                      !e_mode_diff_points(3,k,i) = E||

c-----external
      double precision psi_rho,b,x,y,tempe,tpoprho,vflowrho,psif,rhopsi

c-----local
      double precision x_ar(nbulka),y_ar(nbulka),T_ar_ev(nbulka),     
     &tpop_ar(nbulka),vflow_ar(nbulka)![sm/sec]
      double precision psi,rho_loc,zaxis,raxis,theta,
     &costheta,sintheta,step,z,r,rho_edge

      integer i,j,k,i_found_n_hot_roots(n_hot_roots_a),n_hot_roots
 
      double precision N_perp_root_ar(n_hot_roots_a)    !hot plasma roots
      
      double complex e_xyz_nperp_root_ar(3,n_hot_roots_a) ! electric filed 
                                                          ! for each root
                               ! e_nperp_root_ar(1,k) E_x compomnent
                               ! e_nperp_root_ar(2,k) E_y compomnent
                               ! e_nperp_root_ar(3,k) E_z compomnent

      double precision  e_mode_root_ar(3,n_hot_roots_a)  !electric field
                                                         !polarization 
                                                         ! E+,E-,E_parallel
                                                         !at each N_perp_root(i)]
                               !e_mode_root_ar(1,k)      ! E+
                               !e_mode_root_ar(2,k)      ! E-
                               !e_mode_root_ar(3,k)      ! E||

      character*60  text
      character*6 format

      do i=1,n_hot_roots_a 
        i_found_n_hot_roots(i)=0
        do k=1,n_hot_roots_a 
          N_perp_roots_diff_points(k,i)=0.d0
        enddo
      enddo 

c-----find magnetic axis (psi=0, theta=0)
      call zr_psith(psimag,0.d0,zaxis,raxis)
      
      if(myrank.eq.0) then
        write(*,*)'rho_ini_hot_nperp_roots'
        write(*,*)'zaxis,raxis',zaxis,raxis
        write(*,*)'yma,xma',yma,xma
      endif ! myrank=0

      pi=4.d0*datan(1.d0)
      costheta=(r_edge-raxis)/dsqrt((z_edge-zaxis)**2+(r_edge-raxis)**2)
      sintheta=(z_edge-zaxis)/dsqrt((z_edge-zaxis)**2+(r_edge-raxis)**2)
      if(sintheta.gt.0d0) then
         theta=dacos(costheta)
      else
         theta=2*pi-dacos(costheta)
      endif  

      step=rho_step_find_hot_nperp_roots 

      pi=4.d0*datan(1.d0)

      psi=psif(z_edge,r_edge)
      rho_edge=rhopsi(psi)           !initialization

      rho_loc= rho_edge
 10   continue

      psi=psi_rho(rho_loc)

c      write(*,*)'forest.f  rho_loc= ',rho_loc

      call zr_psith(psi,theta,z,r)

      bmod=b(z,r,phi_edge)
      do j=1,nbulk
        x_ar(j)=x(z,r,phi_edge,j)
        y_ar(j)=y(z,r,phi_edge,j)
        if(j.eq.1) y_ar(1)=-y_ar(1) ! negative Y=(omega_ce/omega)
                                    ! for electrons
        T_ar_ev(j)=tempe(z,r,phi_edge,j)*1.d+3 !(eV) averaged temperature
        tpop_ar(j)=tpoprho(rho,j)
        vflow_ar(j)=vflowrho(rho,j)         
      enddo

      call hot_roots_solver(nbulk,T_ar_ev,tpop_ar,vflow_ar,
     &x_ar,y_ar,cnpar,
     &cN_perp_root_max,n_points_root,  
     &N_perp_root_ar,n_hot_roots,
     &e_xyz_nperp_root_ar,e_mode_root_ar)

      if(myrank.eq.0) then
        write(*,*)' rho_ini_hot_nperp_roots after hot_roots_solver'
        write(*,*)'n_hot_roots',n_hot_roots
        write(*,*)'n_hot_roots_a',n_hot_roots_a
      endif ! myrank=0

      do i=1,n_hot_roots_a

         if(myrank.eq.0) then
         write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)
         endif ! myrank=0

         if ((n_hot_roots.eq.i).and.(i_found_n_hot_roots(i).eq.0))then
c---------------------------------------------------------
c          coordinates of the first point along the line
c          where the solution D_hot(N_perp)=0 had 'i' roots
c-------------------------------------------------------           
           z_ini(i)=z
           r_ini(i)=r
           psi=psif(z,r)
           rho_ini(i)=rhopsi(psi)
           i_found_n_hot_roots(i)=1
           do k=1,n_hot_roots_a 
             N_perp_roots_diff_points(k,i)=N_perp_root_ar(k)

             do j=1,3
                e_xyz_diff_points(j,k,i)= e_xyz_nperp_root_ar(j,k)
                e_mode_diff_points(j,k,i)=e_mode_root_ar(j,k)
             enddo
 
           enddo
         endif
      enddo

      if  (i_found_n_hot_roots(3).eq.1) goto 20 ! point with three roots
                                              ! was found 
      rho_loc=rho_loc-step
        
      if (rho_loc.lt.rho_min_find_hot_nperp_roots) then 
         if(myrank.eq.0) then
           write(*,*)'forest.f rho_ini_hot_nperp_roots'
           write(*,*)'did not find 3 roots at interval'
           write(*,*)'(rho_min_find_hot_nperp_roots,rho_edge)'
         endif ! myrank=0
         goto 20 
c         stop
      endif

      go to 10                   
    
 20   continue

      if(myrank.ne.0) return
      
      write(*,*)'open.f before open10' 
      open(10,file='find_hot_roots.dat')
      format='d21.15'
      text='(1X,"r_edge=",'//format//')'
      write(10,fmt=text)r_edge
      text='(1X,"z_edge=",'//format//')'
      write(10,fmt=text)z_edge
      text='(1X,"cnpar=",'//format//')'
      write(10,fmt=text)cnpar

      do i=1,n_hot_roots_a !loop over space points along rho

         write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)

         if ((i_found_n_hot_roots(i).eq.1).and.(i.eq.1))then
            write(*,*)'i_found_n_hot_roots(i).eq.1)'
          write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)

       text='(1X,"one root was found along small radius in point with")'
            write(10,fmt=text)
            text='(1X,"rho=",'//format//')'
            write(10,fmt=text)rho_ini(i)
            text='(1X,"r=",'//format//')'
            write(10,fmt=text)r_ini(i)
            text='(1X,"z=",'//format//')'
            write(10,fmt=text)z_ini(i)
c-------------------------------------------------------------------------------            
            text='(1X,"root N_perp=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(1,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,1,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,1,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,1,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,1,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,1,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,1,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,1,i)
c-----------------------------------------------------------------------------

         endif

         if ((i_found_n_hot_roots(i).eq.1).and.(i.eq.2))then
           write(*,*)'i_found_n_hot_roots(i) i.eq.2)'
         write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)

          text='(1X,"two roots was found along small radius in point")'
            write(10,fmt=text)
            text='(1X,"rho=",'//format//')'
            write(10,fmt=text)rho_ini(i)
            text='(1X,"r=",'//format//')'
            write(10,fmt=text)r_ini(i)
            text='(1X,"z=",'//format//')'
            write(10,fmt=text)z_ini(i)
c------------------------------------------------------------------------
            text='(1X,"root N_perp(1)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(1,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,1,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,1,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,1,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,1,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,1,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,1,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,1,i)
c------------------------------------------------------------------------
            text='(1X,"root N_perp(2)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(2,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,2,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,2,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,2,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,2,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,2,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,2,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,2,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,2,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,2,i)
c--------------------------------------------------------------------------
         endif

         if ((i_found_n_hot_roots(i).eq.1).and.(i.eq.3))then
          write(*,*)'i_found_n_hot_roots(i) i.eq.3)'
        write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)

        text='(1X,"three roots was found along small radius in point")'
            write(10,fmt=text)
            text='(1X,"rho=",'//format//')'
            write(10,fmt=text)rho_ini(i)
            text='(1X,"r=",'//format//')'
            write(10,fmt=text)r_ini(i)
            text='(1X,"z=",'//format//')'
            write(10,fmt=text)z_ini(i)
c-------------------------------------------------------------------------
            text='(1X,"root N_perp(1)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(1,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,1,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,1,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,1,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,1,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,1,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,1,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,1,i)
c---------------------------------------------------------------------------
            text='(1X,"root N_perp(2)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(2,i)
             text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,2,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,2,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,2,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,2,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,2,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,2,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,2,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,2,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,2,i)
c---------------------------------------------------------------------------
            text='(1X,"root N_perp(3)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(3,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,3,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,3,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,3,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,3,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,3,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,3,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,3,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,3,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,3,i)
         endif
      enddo !i=1,n_hot_roots_a !
      close(10)
  
      return
      end subroutine rho_ini_hot_nperp_roots



      subroutine hot_nperp_muller(nbulk,mass_ar,x_ar,y_ar,t_av_ar,
     &tpop_ar,vflow_ar,cnpar,cnper,cnprim)
c     INPUTS:
c      nbulk - the total number of the plasma species 
c      mass_ar - the masses of the plasma  species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=average temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar = T_perp/T_parallel  
c      vflow_ar    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.

      implicit none
      include 'param.i' 
c-----input
      integer nbulk,iherm   
      double precision mass_ar(nbulka),x_ar(nbulka),y_ar(nbulka),
     &t_av_ar(nbulka),
     &tpop_ar(nbulka),vflow_ar(nbulka)
      double precision cnpar

c-----output
      real*8 cnper,cnprim

c-----data for common /com_hot_nperp_muller/
      
      real*8 mass_ar_l,x_ar_l,y_ar_l,t_av_ar_l,
     &tpop_ar_l,vflow_ar_l,cnpar_l

      integer nbulk_l,iherm_l

      common /com_hot_nperp_muller/ 
     &mass_ar_l(nbulka),x_ar_l(nbulka),y_ar_l(nbulka),t_av_ar_l(nbulka),
     &tpop_ar_l(nbulka),vflow_ar_l(nbulka),
     &cnpar_l,nbulk_l

c-----locals
      integer i

      real*8 errabs
      integer nsig,nknown,nrts,nguess,itmax,nnew,info(2),ier
      complex*16 roots(2)
c-----external
      complex*16 hot_disper_for_muller
      external hot_disper_for_muller

c-----set data to common /com_hot_nperp_muller/

      do i=1,nbulk
         mass_ar_l(i)=mass_ar(i)
         x_ar_l(i)=x_ar(i)
	 y_ar_l(i)=y_ar(i)
	 t_av_ar_l(i)=t_av_ar(i)     ! ev 
         tpop_ar_l(i)=tpop_ar(i)
         vflow_ar_l(i)=vflow_ar(i)
      enddo
      cnpar_l=cnpar
      nbulk_l=nbulk

      errabs=1.d-6
      nsig=6
      nknown=0
      nrts=1                 ! just look for one root
      nguess=nrts 
      nnew=nrts
      itmax=50

      roots(1)=dcmplx(cnper,0.0d0)
      !write(*,*)' hot_nperp_muller before call muller'
      call muller (hot_disper_for_muller,
     & errabs,nsig,nknown,nguess,nnew,roots,
     & itmax,info,ier)

      !write(*,*)' hot_nperp_muller after call muller'

      cnprim=abs(imag(roots(1)))
      cnper=abs(dble(roots(1)))

      return
      end subroutine hot_nperp_muller
      
      
      
      

      complex*16  function hot_disper_for_muller(cnperp,dum)
      !For iabsorp=4
      implicit none
      include 'param.i'

c-----data from common /com_hot_nperp_muller/
      real*8 mass_ar_l,x_ar_l,y_ar_l,t_av_ar_l,
     &tpop_ar_l,vflow_ar_l,cnpar_l
      integer nbulk_l,iherm_l
      common /com_hot_nperp_muller/
     &mass_ar_l(nbulka),x_ar_l(nbulka),y_ar_l(nbulka),t_av_ar_l(nbulka),
     &tpop_ar_l(nbulka),vflow_ar_l(nbulka),
     &cnpar_l,nbulk_l
   
      real*8 
     &mass_ar_ll(nbulka),x_ar_ll(nbulka),y_ar_ll(nbulka),
     &t_av_ar_ll(nbulka),
     &tpop_ar_ll(nbulka),vflow_ar_ll(nbulka),
     &cnpar_ll,cnperp_ll
      integer nbulk_ll
c-----input
      complex*16 cnperp,dum(1:3,1:3)

c-----extrnals
      complex*16 dhot_sum_c

c-----locals
      integer i
      complex*16 reps_l(3,3)
      !write(*,*)'in hot_disper_for_muller'
c      write(*,*)'nbulk_l',nbulk_l
c      write(*,*)'mass_ar_l',mass_ar_l
c      write(*,*)'x_ar_l',x_ar_l
      !write(*,*)'y_ar_l',y_ar_l
c      write(*,*)'tpop_ar_l',tpop_ar_l
c      write(*,*)'t_av_ar_l',t_av_ar_l
c      write(*,*)'vflow_ar_l',vflow_ar_l
      !write(*,*)'cnpar_l',cnpar_l
      !write(*,*)'cnperp',cnperp

      nbulk_ll=nbulk_l
      do i=1,nbulk_ll
         mass_ar_ll(i)=mass_ar_l(i)
         x_ar_ll(i)=x_ar_l(i)
         y_ar_ll(i)=y_ar_l(i)
         t_av_ar_ll(i)=t_av_ar_l(i)
         tpop_ar_ll(i)=tpop_ar_l(i)
         vflow_ar_ll(i)=vflow_ar_l(i) !YuP[2020-01-21] small bug 
      enddo
      cnpar_ll=cnpar_l

c      hot_disper_for_muller=dhot_sum_c(nbulk_l,mass_ar_l,
c     &x_ar_l,y_ar_l,t_av_ar_l,
c     &tpop_ar_l,vflow_ar_l,cnpar_l,cnperp,2,reps_l)
      cnperp_ll=real(cnperp)
      hot_disper_for_muller=dhot_sum_c(nbulk_ll, mass_ar_ll,
     &x_ar_ll, y_ar_ll, t_av_ar_ll,
     &tpop_ar_ll, vflow_ar_ll, cnpar_ll, cnperp_ll, 0.d0, 2, reps_l)

c      write(*,*)'hot_disper_for_muller',hot_disper_for_muller
    

      return
      end function hot_disper_for_muller
