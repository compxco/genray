
c Following is from the brambilla LH and FW code
ctk--calculation of Bessel function with a complex argument (x,y)
c     br(nb),bi(nb)  
c     y=0,ize=0 >J_0(x)=br(1),J_1(x)=br(2)
c     y=0,ize=1 >I_0(x)=br(1),I_1(x)=br(2),... ,I_n(x)=br(n+1),
c     n=0,1,2,...; nb.ge.1
C     x,y,br,bi - real

      subroutine beslci(x,y,nb,ize,br,bi,ncalc)
c
c*****     *****     *****     *****     *****     *****     *****
c
c     journal of research  national bureau of standards
c     b-mathematical sciences   vol.77-b  nos.3-4  july-dec. 1973  p118
c     routine to calculate bessel functions j and i
c     of complex argument and integer order
c                  input variables
c     x real part of complex argument
c     y imaginary part of complex argument
c     nb a positive integer designating highest order to be calculated
c     ize     zero for j*s    one for i*s
c     br  for normal exit,br contains real part of solution vector
c     bi  for normal exit,bi contains imaginary part of solution vector
c     normal exit if  ncalc=nb
c     br-bi-ncalc need not be initialized
c                  machine dependent constants
c     nsig  decimal significance desired. set to ifix(alog10(2)*nbit+1)
c     nbit  number of bits in the mantissa
c     the relative truncation error is limited to t=.5*10**-nsig for
c     order greater than abs(z).
c     for order less than abs(z) (general test),the relative error is
c     limited to t for function values of magnitude at least 1.
c     the absolute error is limited to t for smaller values.
c     nten  largest integer k such that 10**k is machine representable.
c     largez  upper limit on the magnitude of z. if abs(z)=n, at least
c     n iterations of the backward recursion will be executed.
c     exparg  largest argument that the library exp routine can handle.
c                  error returns
c     let g denote either i or j.
c     in case of an error, ncalc.ne.nb and not all g*s are calculated
c     to the desired accuracy.
c     if ncalc.lt.0, an argument is out of range. nb.le.0 or ize is
c     neither 0 nor 1 or ize=0 and abs(y).gt.exparg, or ize=1 and
c     abs(x).gt.exparg. in this case,the vectors br and bi are not
c     calculated, and ncalc is set to min0(nb,0)-1 so ncalc.ne.nb.
c     nb.gt.ncalc.gt.0 will occur if nb.gt.magz and abs(g-sub-nb-of-z/g
c     -sub-magx-of-z).lt.10.**(nten/2), i.e. nb is much greater than
c     magz. in this case, br(n) and bi(n) are calculated to the desired
c     accuracy for n.le.ncalc, but for ncalc.lt.n.le.nb, precision  is
c     lost. if n.gt.ncalc and abs(g(ncalc-1)/g(n-1)).eq.10**-k, then
c     the last k significant figures of g(n-1) (=br(n)+i*bi(n)) are
c     erroneous. if the user wishes to calculate g(n-1) to higher
c     accuracy, he should use an asymptotic formula for large order.
c
c*****     *****     *****     *****     *****     *****     *****
c
      dimension br(*),bi(*)
c      data nsig,nten,largez,exparg/15,307,10000,700./
      data nsig,nten,largez,exparg/15,38,10000,81./
      tempar=sqrt(x*x+y*y)
      magz=ifix(tempar)
      if(nb.gt.0.and.magz.le.largez.and.((ize.eq.0.and.abs(y).le.exparg)
     c.or.(ize.eq.1.and.abs(x).le.exparg))) go to 1
c     error return    z,nb,or ize is out of range
      ncalc=min0(nb,0)-1
      return
    1 sign1=1-2*ize
      ncalc=nb
c     use 2-term ascending series for small z
      if(tempar**4.lt.0.1**nsig) go to 50
c     initialize the calculation of the p*s
      nbmz=nb-magz
      n=magz+1
      if(abs(x).lt.abs(y)) go to 2
      zinvr=1.0/(x+y*y/x)
      zinvi=-y*zinvr/x
      go to 3
    2 zinvi=-1.0/(y+x*x/y)
      zinvr=-x*zinvi/y
    3 plastr=1.0
      plasti=0.0
      pr=sign1*(n+n)*zinvr
      pi=sign1*(n+n)*zinvi
      test=2.0*(10.0)**nsig
      m=0
      if(nbmz.lt.3) go to 6
c     calculate p*s until n=nb-1.  check for possible overflow.
      tover=10.0**(nten-nsig)
      nstart=magz+2
      nend=nb-1
      do 5 n=nstart,nend
      poldr=plastr
      poldi=plasti
      plastr=pr
      plasti=pi
      pr=sign1*((n+n)*(plastr*zinvr-plasti*zinvi)-poldr)
      pi=sign1*((n+n)*(plasti*zinvr+plastr*zinvi)-poldi)
      if((pr/tover)**2+(pi/tover)**2-1.0) 5,5,7
    5 continue
      n=nend
c     calculate special significance test for nbmz.gt.2.
      tempbi=amax1(abs(pr),abs(pi))
      tempbi=tempbi*sqrt(2.0*(10.0)**nsig*sqrt(((pr/tempbi)**2+(pi
     c/tempbi)**2)*((plastr/tempbi)**2+(plasti/tempbi)**2)))
      test=amax1(test,tempbi)
c     calculate p*s until significance test is passed
    6 n=n+1
      poldr=plastr
      poldi=plasti
      plastr=pr
      plasti=pi
      pr=sign1*((n+n)*(plastr*zinvr-plasti*zinvi)-poldr)
      pi=sign1*((n+n)*(plasti*zinvr+plastr*zinvi)-poldi)
      if((pr/test)**2+(pi/test)**2.lt.1.0) go to 6
      if(m.eq.1) go to 12
c     calculate strict variant of significance test,and
c     calculate p*s until this test is passed.
      m=1
      tempbi=amax1(abs(pr),abs(pi))
      tempbr=sqrt(((pr/tempbi)**2+(pi/tempbi)**2)/
     c((plastr/tempbi)**2+(plasti/tempbi)**2))
      tempbi=(n+1)/tempar
      if(tempbr+1.0/tempbr.gt.2.0*tempbi) tempbr=tempbi+sqrt(tempbi**2
     c-1.0)
      test=test/sqrt(tempbr-1.0/tempbr)
      if((pr/test)**2+(pi/test)**2-1.0) 6,12,12
    7 nstart=n+1
c     to avoid overflow, normalize p*s by dividing by tover.
c     calculate p*s until unnormalized p would overflow.
      pr=pr/tover
      pi=pi/tover
      plastr=plastr/tover
      plasti=plasti/tover
      psaver=pr
      psavei=pi
      tempcr=plastr
      tempci=plasti
      test=10.0**(2*nsig)
    8 n=n+1
      poldr=plastr
      poldi=plasti
      plastr=pr
      plasti=pi
      pr=sign1*((n+n)*(plastr*zinvr-plasti*zinvi)-poldr)
      pi=sign1*((n+n)*(plasti*zinvr+plastr*zinvi)-poldi)
      if(pr**2+pi**2.le.test) go to 8
c     calculate backward test,and find ncalc,the highest n
c     such that the test is passed.
      tempbr=sqrt((plastr**2+plasti**2)/(poldr**2+poldi**2))
      tempbi=n/tempar
      if(tempbr+1.0/tempbr.gt.2.0*tempbi) tempbr=tempbi+sqrt(tempbi**2
     c-1.0)
      test=0.5*(1.0-1.0/tempbr**2)/10.0**nsig
      test=((plastr**2+plasti**2)*test)*((poldr**2+poldi**2)*test)
      pr=plastr*tover
      pi=plasti*tover
      n=n-1
      nend=min0(nb,n)
      do 9 ncalc=nstart,nend
      poldr=tempcr
      poldi=tempci
      tempcr=psaver
      tempci=psavei
      psaver=sign1*((n+n)*(tempcr*zinvr-tempci*zinvi)-poldr)
      psavei=sign1*((n+n)*(tempci*zinvr+tempcr*zinvi)-poldi)
      if((psaver**2+psavei**2)*(tempcr**2+tempci**2)-test) 9,9,10
    9 continue
      ncalc=nend+1
   10 ncalc=ncalc-1
c     the coefficient of b(n) in the normalized sum is
c     m*sqrt(-1)**imag,where m=-2,0,or2,and imag is 0 or 1.
c     calculate recursion rules for m and imag,and initialize them.
   12 n=n+1
      tempbr=ize*x+(1-ize)*y
      ipos=0
      if(tempbr) 13,14,13
   13 ipos=ifix(1.1*tempbr/abs(tempbr))
   14 mrecur=4*((2+ize+ipos)/2)-3-2*(ize+ipos)
      k=2+ipos+2*ize*ipos**2-ize
      l=n-4*(n/4)
      mlast=2+8*((k*l)/4)-4*((k*l)/2)
      if(ipos.eq.0.and.(l.eq.1.or.l.eq.3)) mlast=0
      l=l+3-4*((l+3)/4)
      m=2+8*((k*l)/4)-4*((k*l)/2)
      if(ipos.eq.0.and.(l.eq.1.or.l.eq.3)) m=0
      imrecr=(1-ize)*ipos**2
      imag=imrecr*(l-2*(l/2))
c     initialize the backward recursion and the normalization sum
      tempbr=0.0
      tempbi=0.0
      if(abs(pi).gt.abs(pr)) go to 15
      tempar=1.0/(pr+pi*(pi/pr))
      tempai=-(pi*tempar)/pr
      go to 16
   15 tempai=-1.0/(pi+pr*(pr/pi))
      tempar=-(pr*tempai)/pi
   16 if(imag.ne.0) go to 17
      sumr=m*tempar
      sumi=m*tempai
      go to 18
   17 sumr=-m*tempai
      sumi=m*tempar
   18 nend=n-nb
      if(nend) 26,22,19
c     recur backward via difference equation calculating (but not
c     storing) br(n) and bi(n) until n=nb
   19 do 21 l=1,nend
      n=n-1
      tempcr=tempbr
      tempci=tempbi
      tempbr=tempar
      tempbi=tempai
      pr=(n+n)*zinvr
      pi=(n+n)*zinvi
      tempar=pr*tempbr-pi*tempbi-sign1*tempcr
      tempai=pr*tempbi+pi*tempbr-sign1*tempci
      imag=(1-imag)*imrecr
      k=mlast
      mlast=m
      m=k*mrecur
      if(imag.ne.0) go to 20
      sumr=sumr+m*tempar
      sumi=sumi+m*tempai
      go to 21
   20 sumr=sumr-m*tempai
      sumi=sumi+m*tempar
   21 continue
c     store  br(nb),bi(nb)
   22 br(n)=tempar
      bi(n)=tempai
      if(n.gt.1) go to 23
c     nb=1.  since 2*tempar and 2*tempai were added to sumr and sumi
c     respectively,we must subaxis(i)ract tempar and tempai
      sumr=sumr-tempar
      sumi=sumi-tempai
      go to 35
c     calculate and store br(nb-1),bi(nb-1)
   23 n=n-1
      pr=(n+n)*zinvr
      pi=(n+n)*zinvi
      br(n)=pr*tempar-pi*tempai-sign1*tempbr
      bi(n)=pr*tempai+pi*tempar-sign1*tempbi
      if(n.eq.1) go to 34
      imag=(1-imag)*imrecr
      k=mlast
      mlast=m
      m=k*mrecur
      if(imag.ne.0) go to 24
      sumr=sumr+m*br(n)
      sumi=sumi+m*bi(n)
      go to 30
   24 sumr=sumr-m*bi(n)
      sumi=sumi+m*br(n)
      go to 30
c     n.lt.nb,so store br(n),bi(n) and set higher orders zero
   26 br(n)=tempar
      bi(n)=tempai
      nend=-nend
      do 27 l=1,nend
      br(n+l)=0.0
   27 bi(n+l)=0.0
   30 nend=n-2
      if(nend.eq.0) go to 33
c     calculate via difference equation and store br(n),bi(n)
c     until n=2
      do 32 l=1,nend
      n=n-1
      pr=(n+n)*zinvr
      pi=(n+n)*zinvi
      br(n)=pr*br(n+1)-pi*bi(n+1)-sign1*br(n+2)
      bi(n)=pr*bi(n+1)+pi*br(n+1)-sign1*bi(n+2)
      imag=(1-imag)*imrecr
      k=mlast
      mlast=m
      m=k*mrecur
      if(imag.ne.0) go to 31
      sumr=sumr+m*br(n)
      sumi=sumi+m*bi(n)
      go to 32
   31 sumr=sumr-m*bi(n)
      sumi=sumi+m*br(n)
   32 continue
c     calculate and store br(1),bi(1)
   33 br(1)=2.0*(br(2)*zinvr-bi(2)*zinvi)-sign1*br(3)
      bi(1)=2.0*(br(2)*zinvi+bi(2)*zinvr)-sign1*bi(3)
   34 sumr=sumr+br(1)
      sumi=sumi+bi(1)
c     calculate normalization factor. tempar+i*tempai
   35 if(ize.eq.1) go to 36
      tempcr=ipos*y
      tempci=-ipos*x
      go to 37
   36 tempcr=ipos*x
      tempci=ipos*y
   37 tempcr=exp(tempcr)
      tempbr=cos(tempci)
      tempbi=sin(tempci)
      if(abs(sumr).lt.abs(sumi)) go to 38
      tempci=sumi/sumr
      tempcr=(tempcr/sumr)/(1.0+tempci*tempci)
      tempar=tempcr*(tempbr+tempbi*tempci)
      tempai=tempcr*(tempbi-tempbr*tempci)
      go to 39
   38 tempci=sumr/sumi
      tempcr=(tempcr/sumi)/(1.0+tempci*tempci)
      tempar=tempcr*(tempbr*tempci+tempbi)
      tempai=tempcr*(tempbi*tempci-tempbr)
c     normalize
   39 do 40 n=1,nb
      tempbr=br(n)*tempar-bi(n)*tempai
      bi(n)=br(n)*tempai+bi(n)*tempar
   40 br(n)=tempbr
      return
c     two-term ascending series for small z
   50 tempar=1.0
      tempai=0.0
      tempcr=0.25*(x*x-y*y)
      tempci=0.5*x*y
      br(1)=1.0-sign1*tempcr
      bi(1)=-sign1*tempci
      if(nb.eq.1) go to 52
      do 51 n=2,nb
      tempbr=(tempar*x-tempai*y)/(n+n-2)
      tempai=(tempar*y+tempai*x)/(n+n-2)
      tempar=tempbr
      tempbr=n
      br(n)=tempar*(1.0-sign1*tempcr/tempbr)+tempai*tempci/tempbr
   51 bi(n)=tempai*(1.0-sign1*tempcr/tempbr)-tempar*tempci/tempbr
   52 return
      end




      subroutine besj(x,n,bj,d,ier)
c-----------------------------------------------------------------------
c     i.b.m. bessel function j
c     x  - argument
c     n  - order (integer .ge. 0)
c     bj - j(x)
c     d  - relative allowable error (input)
c     ier- error switch, if ier=0 all is ok
c-----------------------------------------------------------------------
      bj=0.
      if(n) 10,20,20
10     ier=1
      return
20    if(x) 30,30,31
30    ier=2
      bj=0.
      if(n.eq.0) bj=1.
      return
31    if(x-15.) 32,32,34
32    ntest=20.+10.*x-x*x/3
      go to 36
34    ntest=90.+x/2.
36    if(n-ntest) 40,38,38
38    ier=4
      return
40    ier=0
      n1=n+1
      bprev=0.
      if(x-5.) 50,60,60
50    ma=x+6.
      go to 70
60    ma=1.4*x+60./x
70    mb=n+ifix(x)/4+2
      mzero=max0(ma,mb)
      mmax=ntest
100    do 190 m=mzero,mmax,3
      fm1=1.e-28
      fm=0.
      alpha=0.
      if(m-(m/2)*2)120,110,120
110   jt=-1
      go to 130
120   jt=1
130   m2=m-2
      do 160 k=1,m2
      mk=m-k
      bmk=2.*float(mk)*fm1/x-fm
      fm=fm1
      fm1=bmk
      if(mk-n-1) 150,140,150
140   bj=bmk
150   jt=-jt
      s=1+jt
160   alpha=alpha+bmk*s
      bmk=2.*fm1/x-fm
      if(n) 180,170,180
170   bj=bmk
180   alpha=alpha+bmk
      bj=bj/alpha
      if(abs(bj-bprev)-abs(d*bj))200,200,190
190   bprev=bj
      ier=3
200   return
      end
c
c


c88888888888888888888888888888888888888888888888888888888888888888888888888


      SUBROUTINE DBESJ(X,ALPHA,N,Y,NZ)
C***BEGIN PROLOGUE  DBESJ
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  851111   (YYMMDD)
C***CATEGORY NO.  C10A3
C***KEYWORDS  BESSEL FUNCTION,DOUBLE PRECISION,J BESSEL FUNCTION,
C             SPECIAL FUNCTION
C***AUTHOR  AMOS, D. E., (Sandia National Laboratories, Albuquerque)
C           DANIEL, S. L., (Sandia National Laboratories, Albuquerque)
C           WESTON, M. K., (Sandia National Laboratories, Albuquerque)
C***PURPOSE  Compute an N member sequence of J Bessel functions
C            J/SUB(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA
C            and X. (At most 14 digits.)
C***DESCRIPTION
C
C     Written by D. E. Amos, S. L. Daniel and M. K. Weston, January 1975
C
C     References
C         SAND-75-0147
C
C         CDC 6600 Subroutines IBESS and JBESS for Bessel Functions
C         I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0  by D. E. Amos, S. L.
C         Daniel, M. K. Weston. ACM Trans Math Software,3,pp 76-92
C         (1977)
C
C         Tables of Bessel Functions of Moderate or Large Orders,
C         NPL Mathematical Tables, Vol. 6, by F. W. J. Olver, Her
C         Majesty's Stationery Office, London, 1962.
C
C     Abstract  **** a double precision routine ****
C         DBESJ computes an N member sequence of J Bessel functions
C         J/sub(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA and X.
C         A combination of the power series, the asymptotic expansion
C         for X to infinity and the uniform asymptotic expansion for
C         NU to infinity are applied over subdivisions of the (NU,X)
C         plane.  For values of (NU,X) not covered by one of these
C         formulae, the order is incremented or decremented by integer
C         values into a region where one of the formulae apply. Backward
C         recursion is applied to reduce orders by integer values except
C         where the entire sequence lies in the oscillatory region.  In
C         this case forward recursion is stable and values from the
C         asymptotic expansion for X to infinity start the recursion
C         when it is efficient to do so. Leading terms of the series and
C         uniform expansion are tested for underflow.  If a sequence is
C         requested and the last member would underflow, the result is
C         set to zero and the next lower order tried, etc., until a
C         member comes on scale or all members are set to zero.
C         Overflow cannot occur.
C
C         The maximum number of significant digits obtainable
C         is the smaller of 14 and the number of digits carried in
C         double precision arithmetic.
C
C         DBESJ calls DASYJY, DJAIRY, DLNGAM, D1MACH, I1MACH, XERROR
C
C     Description of Arguments
C
C         Input      X,ALPHA are double precision
C           X      - X .GE. 0.0D0
C           ALPHA  - order of first member of the sequence,
C                    ALPHA .GE. 0.0D0
C           N      - number of members in the sequence, N .GE. 1
C
C         Output     Y is double precision
C           Y      - a vector whose first N components contain
C                    values for J/sub(ALPHA+K-1)/(X), K=1,...,N
C           NZ     - number of components of Y set to zero due to
C                    underflow,
C                    NZ=0   , normal return, computation completed
C                    NZ .NE. 0, last NZ components of Y set to zero,
C                             Y(K)=0.0D0, K=N-NZ+1,...,N.
C
C     Error Conditions
C         Improper input arguments - a fatal error
C         Underflow  - a non-fatal error (NZ .NE. 0)
C***REFERENCES  CDC 6600 SUBROUTINES IBESS AND JBESS FOR BESSEL
C                 FUNCTIONS I(NU,X) AND J(NU,X), X .GE. 0, NU .GE. 0,
C                 BY D. E. AMOS, S. L.DANIEL, M. K. WESTON,  ACM
C                 TRANSACTIONS ON MATHEMATICALSOFTWARE, VOL. 3,
C                 PP. 76-92 (1977).
C***ROUTINES CALLED  D1MACH,DASYJY,DJAIRY,DLNGAM,I1MACH,XERROR
C***END PROLOGUE  DBESJ
      EXTERNAL DJAIRY
      INTEGER I,IALP,IDALP,IFLW,IN,INLIM,IS,I1,I2,K,KK,KM,KT,N,NN,
     1        NS,NZ
      INTEGER I1MACH
      DOUBLE PRECISION AK,AKM,ALPHA,ANS,AP,ARG,COEF,DALPHA,DFN,DTM,
     1           EARG,ELIM1,ETX,FIDAL,FLGJY,FN,FNF,FNI,FNP1,FNU,
     2           FNULIM,GLN,PDF,PIDT,PP,RDEN,RELB,RTTP,RTWO,RTX,RZDEN,
     3           S,SA,SB,SXO2,S1,S2,T,TA,TAU,TB,TEMP,TFN,TM,TOL,
     4           TOLLN,TRX,TX,T1,T2,WK,X,XO2,XO2L,Y,SLIM,RTOL
      DOUBLE PRECISION D1MACH, DLNGAM
      DIMENSION Y(*), TEMP(3), FNULIM(2), PP(4), WK(7)
      DATA RTWO,PDF,RTTP,PIDT                    / 1.34839972492648D+00,
     1 7.85398163397448D-01, 7.97884560802865D-01, 1.57079632679490D+00/
      DATA  PP(1),  PP(2),  PP(3),  PP(4)        / 8.72909153935547D+00,
     1 2.65693932265030D-01, 1.24578576865586D-01, 7.70133747430388D-04/
      DATA INLIM           /      150            /
      DATA FNULIM(1), FNULIM(2) /      100.0D0,     60.0D0     /
C***FIRST EXECUTABLE STATEMENT  DBESJ
      
cBH      write(*,*)'DBESJ:x,alpha,n,y(1),nzero',x,alpha,n,y(1),nz
cBH
cBH      write(*,*)'DBESJ:D1MACH(1:5)',D1MACH(1),D1MACH(2),D1MACH(3),
cBH     +          D1MACH(4),D1MACH(5)


      NZ = 0
      KT = 1
      NS=0
C     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
C     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
      TA = D1MACH(3)
      TOL = DMAX1(TA,1.0D-15)
      I1 = I1MACH(14) + 1
      I2 = I1MACH(15)
      TB = D1MACH(5)

cBH      write(*,*)'DBESJ:I1,I2,TB',I1,I2,TB
 
      ELIM1 = 2.303D0*(DBLE(FLOAT(-I2))*TB-3.0D0)
      RTOL=1.0D0/TOL
      SLIM=D1MACH(1)*RTOL*1.0D+3
C     TOLLN = -LN(TOL)
      TOLLN = 2.303D0*TB*DBLE(FLOAT(I1))
      TOLLN = DMIN1(TOLLN,34.5388D0)
      IF (N-1) 720, 10, 20
   10 KT = 2
   20 NN = N
      IF (X) 730, 30, 80
   30 IF (ALPHA) 710, 40, 50
   40 Y(1) = 1.0D0
      IF (N.EQ.1) RETURN
      I1 = 2
      GO TO 60
   50 I1 = 1
   60 DO 70 I=I1,N
        Y(I) = 0.0D0
   70 CONTINUE
      RETURN
   80 CONTINUE
      IF (ALPHA.LT.0.0D0) GO TO 710
C
      IALP = INT(SNGL(ALPHA))
      FNI = DBLE(FLOAT(IALP+N-1))
      FNF = ALPHA - DBLE(FLOAT(IALP))
      DFN = FNI + FNF
      FNU = DFN
      XO2 = X*0.5D0
      SXO2 = XO2*XO2
C
C     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
C     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
C     APPLIED.
C
      IF (SXO2.LE.(FNU+1.0D0)) GO TO 90
      TA = DMAX1(20.0D0,FNU)
      IF (X.GT.TA) GO TO 120
      IF (X.GT.12.0D0) GO TO 110
      XO2L = DLOG(XO2)
      NS = INT(SNGL(SXO2-FNU)) + 1
      GO TO 100
   90 FN = FNU
      FNP1 = FN + 1.0D0
      XO2L = DLOG(XO2)
      IS = KT
      IF (X.LE.0.50D0) GO TO 330
      NS = 0
  100 FNI = FNI + DBLE(FLOAT(NS))
      DFN = FNI + FNF
      FN = DFN
      FNP1 = FN + 1.0D0
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 330
  110 ANS = DMAX1(36.0D0-FNU,0.0D0)
      NS = INT(SNGL(ANS))
      FNI = FNI + DBLE(FLOAT(NS))
      DFN = FNI + FNF
      FN = DFN
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 130
  120 CONTINUE
      RTX = DSQRT(X)
      TAU = RTWO*RTX
      TA = TAU + FNULIM(KT)
      IF (FNU.LE.TA) GO TO 480
      FN = FNU
      IS = KT
C
C     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
C
  130 CONTINUE
      I1 = IABS(3-IS)
      I1 = MAX0(I1,1)
      FLGJY = 1.0D0
      CALL DASYJY(DJAIRY,X,FN,FLGJY,I1,TEMP(IS),WK,IFLW)
      IF(IFLW.NE.0) GO TO 380
      GO TO (320, 450, 620), IS
  310 TEMP(1) = TEMP(3)
      KT = 1
  320 IS = 2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF(I1.EQ.2) GO TO 450
      GO TO 130
C
C     SERIES FOR (X/2)**2.LE.NU+1
C
  330 CONTINUE
      GLN = DLNGAM(FNP1)
      ARG = FN*XO2L - GLN
      IF (ARG.LT.(-ELIM1)) GO TO 400
      EARG = DEXP(ARG)
  340 CONTINUE
      S = 1.0D0
      IF (X.LT.TOL) GO TO 360
      AK = 3.0D0
      T2 = 1.0D0
      T = 1.0D0
      S1 = FN
      DO 350 K=1,17
        S2 = T2 + S1
        T = -T*SXO2/S2
        S = S + T
        IF (DABS(T).LT.TOL) GO TO 360
        T2 = T2 + AK
        AK = AK + 2.0D0
        S1 = S1 + FN
  350 CONTINUE
  360 CONTINUE
      TEMP(IS) = S*EARG
      GO TO (370, 450, 610), IS
  370 EARG = EARG*FN/XO2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IS = 2
      GO TO 340
C
C     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
C     UNDERFLOW CAN ONLY OCCRFOR NS=0 SINCE THE ORDER MUST BE LARGER
C     THAN 36. THEREFORE, NS NEE NOT BE TESTED.
C
  380 Y(NN) = 0.0D0
      NN = NN - 1
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 440, 390, 130
  390 KT = 2
      IS = 2
      GO TO 130
  400 Y(NN) = 0.0D0
      NN = NN - 1
      FNP1 = FN
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 440, 410, 420
  410 KT = 2
      IS = 2
  420 IF (SXO2.LE.FNP1) GO TO 430
      GO TO 130
  430 ARG = ARG - XO2L + DLOG(FNP1)
      IF (ARG.LT.(-ELIM1)) GO TO 400
      GO TO 330
  440 NZ = N - NN
cBH      write(*,*)'**0**DBESJ:Y(1)',Y(1)
      RETURN
C
C     BACKWARD RECURSION SECTION
C
  450 CONTINUE
      IF(NS.NE.0) GO TO 451
      NZ = N - NN
      IF (KT.EQ.2) GO TO 470
C     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(NN) = TEMP(1)
      Y(NN-1) = TEMP(2)
cBH      write(*,*)'**0.1**DBESJ:Y(1)',Y(1)
      IF (NN.EQ.2) RETURN
  451 CONTINUE
      TRX = 2.0D0/X
      DTM = FNI
      TM = (DTM+FNF)*TRX
      AK=1.0D0
      TA=TEMP(1)
      TB=TEMP(2)
      IF(DABS(TA).GT.SLIM) GO TO 455
      TA=TA*RTOL
      TB=TB*RTOL
      AK=TOL
  455 CONTINUE
      KK=2
      IN=NS-1
      IF(IN.EQ.0) GO TO 690
      IF(NS.NE.0) GO TO 670
      K=NN-2
      DO 460 I=3,NN
        S=TB
        TB = TM*TB - TA
        TA=S
        Y(K)=TB*AK
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
        K = K - 1
  460 CONTINUE
      RETURN
  470 Y(1) = TEMP(2)
cBH      write(*,*)'**0.2**DBESJ:Y(1)',Y(1)
      RETURN
C
C     ASYMPTOTIC EXPANSION FOR X TO INFINITY WITH FORWARD RECURSION IN
C     OSCILLATORY REGION X.GT.MAX(20, NU), PROVIDED THE LAST MEMBER
C     OF THE SEQUENCE IS ALSO IN THE REGION.
C
  480 CONTINUE
      IN = INT(SNGL(ALPHA-TAU+2.0D0))
      IF (IN.LE.0) GO TO 490
      IDALP = IALP - IN - 1
      KT = 1
      GO TO 500
  490 CONTINUE
      IDALP = IALP
      IN = 0
  500 IS = KT
      FIDAL = DBLE(FLOAT(IDALP))
      DALPHA = FIDAL + FNF
      ARG = X - PIDT*DALPHA - PDF
      SA = DSIN(ARG)
      SB = DCOS(ARG)
      COEF = RTTP/RTX
      ETX = 8.0D0*X
  510 CONTINUE
      DTM = FIDAL + FIDAL
      DTM = DTM*DTM
      TM = 0.0D0
      IF (FIDAL.EQ.0.0D0 .AND. DABS(FNF).LT.TOL) GO TO 520
      TM = 4.0D0*FNF*(FIDAL+FIDAL+FNF)
  520 CONTINUE
      TRX = DTM - 1.0D0
      T2 = (TRX+TM)/ETX
      S2 = T2
      RELB = TOL*DABS(T2)
      T1 = ETX
      S1 = 1.0D0
      FN = 1.0D0
      AK = 8.0D0
      DO 530 K=1,13
        T1 = T1 + ETX
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = -T2*AP/T1
        S1 = S1 + T2
        T1 = T1 + ETX
        AK = AK + 8.0D0
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = T2*AP/T1
        S2 = S2 + T2
        IF (DABS(T2).LE.RELB) GO TO 540
        AK = AK + 8.0D0
  530 CONTINUE
  540 TEMP(IS) = COEF*(S1*SB-S2*SA)
      IF(IS.EQ.2) GO TO 560
      FIDAL = FIDAL + 1.0D0
      DALPHA = FIDAL + FNF
      IS = 2
      TB = SA
      SA = -SB
      SB = TB
      GO TO 510
C
C     FORWARD RECURSION SECTION
C
  560 IF (KT.EQ.2) GO TO 470
      S1 = TEMP(1)
      S2 = TEMP(2)
      TX = 2.0D0/X
      TM = DALPHA*TX
      IF (IN.EQ.0) GO TO 580
C
C     FORWARD RECUR TO INDEX ALPHA
C
      DO 570 I=1,IN
        S = S2
        S2 = TM*S2 - S1
        TM = TM + TX
        S1 = S
  570 CONTINUE
      IF (NN.EQ.1) GO TO 600
      S = S2
      S2 = TM*S2 - S1
      TM = TM + TX
      S1 = S
  580 CONTINUE
C
C     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1
C
      Y(1) = S1
      Y(2) = S2
      IF (NN.EQ.2) RETURN
      DO 590 I=3,NN
        Y(I) = TM*Y(I-1) - Y(I-2)
        TM = TM + TX
  590 CONTINUE
cBH      write(*,*)'**0.3**DBESJ:Y(1)',Y(1)
      RETURN
  600 Y(1) = S2

cBH      write(*,*)'**1**DBESJ:Y(1)',Y(1)
      RETURN
C
C     BACKWARD RECURSION WITH NORMALIZATION BY
C     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
C
  610 CONTINUE
C     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      AKM = DMAX1(3.0D0-FN,0.0D0)
      KM = INT(SNGL(AKM))
      TFN = FN + DBLE(FLOAT(KM))
      TA = (GLN+TFN-0.9189385332D0-0.0833333333D0/TFN)/(TFN+0.5D0)
      TA = XO2L - TA
      TB = -(1.0D0-1.5D0/TFN)/TFN
      AKM = TOLLN/(-TA+DSQRT(TA*TA-TOLLN*TB)) + 1.5D0
      IN = KM + INT(SNGL(AKM))
      GO TO 660
  620 CONTINUE
C     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      GLN = WK(3) + WK(2)
      IF (WK(6).GT.30.0D0) GO TO 640
      RDEN = (PP(4)*WK(6)+PP(3))*WK(6) + 1.0D0
      RZDEN = PP(1) + PP(2)*WK(6)
      TA = RZDEN/RDEN
      IF (WK(1).LT.0.10D0) GO TO 630
      TB = GLN/WK(5)
      GO TO 650
  630 TB=(1.259921049D0+(0.1679894730D0+0.0887944358D0*WK(1))*WK(1))
     1 /WK(7)
      GO TO 650
  640 CONTINUE
      TA = 0.5D0*TOLLN/WK(4)
      TA=((0.0493827160D0*TA-0.1111111111D0)*TA+0.6666666667D0)*TA*WK(6)
      IF (WK(1).LT.0.10D0) GO TO 630
      TB = GLN/WK(5)
  650 IN = INT(SNGL(TA/TB+1.5D0))
      IF (IN.GT.INLIM) GO TO 310
  660 CONTINUE
      DTM = FNI + DBLE(FLOAT(IN))
      TRX = 2.0D0/X
      TM = (DTM+FNF)*TRX
      TA = 0.0D0
      TB = TOL
      KK = 1
      AK=1.0D0
  670 CONTINUE
C
C     BACKWARD RECUR UNINDEXED
C
      DO 680 I=1,IN
        S = TB
        TB = TM*TB - TA
        TA = S
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
  680 CONTINUE
C     NORMALIZATION
      IF (KK.NE.1) GO TO 690
      S=TEMP(3)
      SA=TA/TB
      TA=S
      TB=S
      IF(DABS(S).GT.SLIM) GO TO 685
      TA=TA*RTOL
      TB=TB*RTOL
      AK=TOL
  685 CONTINUE
      TA=TA*SA
      KK = 2
      IN = NS
      IF (NS.NE.0) GO TO 670
  690 Y(NN) = TB*AK
      NZ = N - NN
cBH      write(*,*)'**1.1**DBESJ:Y(1)',Y(1)
      IF (NN.EQ.1) RETURN
      K = NN - 1
      S=TB
      TB = TM*TB - TA
      TA=S
      Y(K)=TB*AK
cBH      write(*,*)'**1.2**DBESJ:Y(1)',Y(1)
      IF (NN.EQ.2) RETURN
      DTM = DTM - 1.0D0
      TM = (DTM+FNF)*TRX
      K=NN-2
C
C     BACKWARD RECUR INDEXED
C
      DO 700 I=3,NN
        S=TB
        TB = TM*TB - TA
        TA=S
        Y(K)=TB*AK
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
        K = K - 1
  700 CONTINUE
cBH      write(*,*)'**2**DBESJ:Y(1)',Y(1)
      RETURN
C
C
C
  710 CONTINUE
      CALL XERROR('DBESJ - ORDER, ALPHA, LESS THAN ZERO.', 37, 2, 1)
      RETURN
  720 CONTINUE
      CALL XERROR('DBESJ - N LESS THAN ONE.', 24, 2, 1)
      RETURN
  730 CONTINUE
      CALL XERROR('DBESJ - X LESS THAN ZERO.', 25, 2, 1)
      RETURN
      END
      SUBROUTINE DJAIRY(X,RX,C,AI,DAI)
C***BEGIN PROLOGUE  DJAIRY
C***REFER TO  DBESJ,DBESY
C
C     1-2-74
C                  DJAIRY computes the Airy function AI(X)
C                   and its derivative DAI(X) for DASYJY
C
C                                   INPUT
C
C         X - Argument, computed by DASYJY, X unrestricted
C        RX - RX=DSQRT(DABS(X)), computed by DASYJY
C         C - C=2.*(DABS(X)**1.5)/3., computed by DASYJY
C
C                                  OUTPUT
C
C        AI - Value of function AI(X)
C       DAI - Value of the derivative DAI(X)
C
C                                Written by
C
C                                D. E. Amos
C                               S. L. Daniel
C                               M. K. Weston
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DJAIRY
C
      INTEGER I, J, M1, M1D, M2, M2D, M3, M3D, M4, M4D, N1, N1D, N2,
     1 N2D, N3, N3D, N4, N4D
      DOUBLE PRECISION A,AI,AJN,AJP,AK1,AK2,AK3,B,C,CCV,CON1,CON2,
     1 CON3, CON4, CON5, CV, DA, DAI, DAJN, DAJP, DAK1, DAK2, DAK3,
     2 DB, EC, E1, E2, FPI12, F1, F2, RTRX, RX, SCV, T, TEMP1, TEMP2,
     3 TT, X
      DIMENSION AJP(19), AJN(19), A(15), B(15)
      DIMENSION AK1(14), AK2(23), AK3(14)
      DIMENSION DAJP(19), DAJN(19), DA(15), DB(15)
      DIMENSION DAK1(14), DAK2(24), DAK3(14)
      SAVE N1, N2, N3, N4, M1, M2, M3, M4, FPI12, CON1, CON2, CON3,
     1 CON4, CON5, AK1, AK2, AK3, AJP, AJN, A, B,
     2 N1D, N2D, N3D, N4D, M1D, M2D, M3D, M4D, DAK1, DAK2, DAK3,
     3 DAJP, DAJN, DA, DB
      DATA N1,N2,N3,N4/14,23,19,15/
      DATA M1,M2,M3,M4/12,21,17,13/
      DATA FPI12,CON1,CON2,CON3,CON4,CON5/
     1 1.30899693899575D+00, 6.66666666666667D-01, 5.03154716196777D+00,
     2 3.80004589867293D-01, 8.33333333333333D-01, 8.66025403784439D-01/
      DATA AK1(1), AK1(2), AK1(3), AK1(4), AK1(5), AK1(6), AK1(7),
     1     AK1(8), AK1(9), AK1(10),AK1(11),AK1(12),AK1(13),
     2     AK1(14)         / 2.20423090987793D-01,-1.25290242787700D-01,
     3 1.03881163359194D-02, 8.22844152006343D-04,-2.34614345891226D-04,
     4 1.63824280172116D-05, 3.06902589573189D-07,-1.29621999359332D-07,
     5 8.22908158823668D-09, 1.53963968623298D-11,-3.39165465615682D-11,
     6 2.03253257423626D-12,-1.10679546097884D-14,-5.16169497785080D-15/
      DATA AK2(1), AK2(2), AK2(3), AK2(4), AK2(5), AK2(6), AK2(7),
     1     AK2(8), AK2(9), AK2(10),AK2(11),AK2(12),AK2(13),AK2(14),
     2     AK2(15),AK2(16),AK2(17),AK2(18),AK2(19),AK2(20),AK2(21),
     3     AK2(22),AK2(23) / 2.74366150869598D-01, 5.39790969736903D-03,
     4-1.57339220621190D-03, 4.27427528248750D-04,-1.12124917399925D-04,
     5 2.88763171318904D-05,-7.36804225370554D-06, 1.87290209741024D-06,
     6-4.75892793962291D-07, 1.21130416955909D-07,-3.09245374270614D-08,
     7 7.92454705282654D-09,-2.03902447167914D-09, 5.26863056595742D-10,
     8-1.36704767639569D-10, 3.56141039013708D-11,-9.31388296548430D-12,
     9 2.44464450473635D-12,-6.43840261990955D-13, 1.70106030559349D-13,
     1-4.50760104503281D-14, 1.19774799164811D-14,-3.19077040865066D-15/
      DATA AK3(1), AK3(2), AK3(3), AK3(4), AK3(5), AK3(6), AK3(7),
     1     AK3(8), AK3(9), AK3(10),AK3(11),AK3(12),AK3(13),
     2     AK3(14)         / 2.80271447340791D-01,-1.78127042844379D-03,
     3 4.03422579628999D-05,-1.63249965269003D-06, 9.21181482476768D-08,
     4-6.52294330229155D-09, 5.47138404576546D-10,-5.24408251800260D-11,
     5 5.60477904117209D-12,-6.56375244639313D-13, 8.31285761966247D-14,
     6-1.12705134691063D-14, 1.62267976598129D-15,-2.46480324312426D-16/
      DATA AJP(1), AJP(2), AJP(3), AJP(4), AJP(5), AJP(6), AJP(7),
     1     AJP(8), AJP(9), AJP(10),AJP(11),AJP(12),AJP(13),AJP(14),
     2     AJP(15),AJP(16),AJP(17),AJP(18),
     3     AJP(19)         / 7.78952966437581D-02,-1.84356363456801D-01,
     4 3.01412605216174D-02, 3.05342724277608D-02,-4.95424702513079D-03,
     5-1.72749552563952D-03, 2.43137637839190D-04, 5.04564777517082D-05,
     6-6.16316582695208D-06,-9.03986745510768D-07, 9.70243778355884D-08,
     7 1.09639453305205D-08,-1.04716330588766D-09,-9.60359441344646D-11,
     8 8.25358789454134D-12, 6.36123439018768D-13,-4.96629614116015D-14,
     9-3.29810288929615D-15, 2.35798252031104D-16/
      DATA AJN(1), AJN(2), AJN(3), AJN(4), AJN(5), AJN(6), AJN(7),
     1     AJN(8), AJN(9), AJN(10),AJN(11),AJN(12),AJN(13),AJN(14),
     2     AJN(15),AJN(16),AJN(17),AJN(18),
     3     AJN(19)         / 3.80497887617242D-02,-2.45319541845546D-01,
     4 1.65820623702696D-01, 7.49330045818789D-02,-2.63476288106641D-02,
     5-5.92535597304981D-03, 1.44744409589804D-03, 2.18311831322215D-04,
     6-4.10662077680304D-05,-4.66874994171766D-06, 7.15218807277160D-07,
     7 6.52964770854633D-08,-8.44284027565946D-09,-6.44186158976978D-10,
     8 7.20802286505285D-11, 4.72465431717846D-12,-4.66022632547045D-13,
     9-2.67762710389189D-14, 2.36161316570019D-15/
      DATA A(1),   A(2),   A(3),   A(4),   A(5),   A(6),   A(7),
     1     A(8),   A(9),   A(10),  A(11),  A(12),  A(13),  A(14),
     2     A(15)           / 4.90275424742791D-01, 1.57647277946204D-03,
     3-9.66195963140306D-05, 1.35916080268815D-07, 2.98157342654859D-07,
     4-1.86824767559979D-08,-1.03685737667141D-09, 3.28660818434328D-10,
     5-2.57091410632780D-11,-2.32357655300677D-12, 9.57523279048255D-13,
     6-1.20340828049719D-13,-2.90907716770715D-15, 4.55656454580149D-15,
     7-9.99003874810259D-16/
      DATA B(1),   B(2),   B(3),   B(4),   B(5),   B(6),   B(7),
     1     B(8),   B(9),   B(10),  B(11),  B(12),  B(13),  B(14),
     2     B(15)           / 2.78593552803079D-01,-3.52915691882584D-03,
     3-2.31149677384994D-05, 4.71317842263560D-06,-1.12415907931333D-07,
     4-2.00100301184339D-08, 2.60948075302193D-09,-3.55098136101216D-11,
     5-3.50849978423875D-11, 5.83007187954202D-12,-2.04644828753326D-13,
     6-1.10529179476742D-13, 2.87724778038775D-14,-2.88205111009939D-15,
     7-3.32656311696166D-16/
      DATA N1D,N2D,N3D,N4D/14,24,19,15/
      DATA M1D,M2D,M3D,M4D/12,22,17,13/
      DATA DAK1(1), DAK1(2), DAK1(3), DAK1(4), DAK1(5), DAK1(6),
     1     DAK1(7), DAK1(8), DAK1(9), DAK1(10),DAK1(11),DAK1(12),
     2    DAK1(13),DAK1(14)/ 2.04567842307887D-01,-6.61322739905664D-02,
     3-8.49845800989287D-03, 3.12183491556289D-03,-2.70016489829432D-04,
     4-6.35636298679387D-06, 3.02397712409509D-06,-2.18311195330088D-07,
     5-5.36194289332826D-10, 1.13098035622310D-09,-7.43023834629073D-11,
     6 4.28804170826891D-13, 2.23810925754539D-13,-1.39140135641182D-14/
      DATA DAK2(1), DAK2(2), DAK2(3), DAK2(4), DAK2(5), DAK2(6),
     1     DAK2(7), DAK2(8), DAK2(9), DAK2(10),DAK2(11),DAK2(12),
     2     DAK2(13),DAK2(14),DAK2(15),DAK2(16),DAK2(17),DAK2(18),
     3     DAK2(19),DAK2(20),DAK2(21),DAK2(22),DAK2(23),
     4     DAK2(24)        / 2.93332343883230D-01,-8.06196784743112D-03,
     5 2.42540172333140D-03,-6.82297548850235D-04, 1.85786427751181D-04,
     6-4.97457447684059D-05, 1.32090681239497D-05,-3.49528240444943D-06,
     7 9.24362451078835D-07,-2.44732671521867D-07, 6.49307837648910D-08,
     8-1.72717621501538D-08, 4.60725763604656D-09,-1.23249055291550D-09,
     9 3.30620409488102D-10,-8.89252099772401D-11, 2.39773319878298D-11,
     1-6.48013921153450D-12, 1.75510132023731D-12,-4.76303829833637D-13,
     2 1.29498241100810D-13,-3.52679622210430D-14, 9.62005151585923D-15,
     3-2.62786914342292D-15/
      DATA DAK3(1), DAK3(2), DAK3(3), DAK3(4), DAK3(5), DAK3(6),
     1     DAK3(7), DAK3(8), DAK3(9), DAK3(10),DAK3(11),DAK3(12),
     2    DAK3(13),DAK3(14)/ 2.84675828811349D-01, 2.53073072619080D-03,
     3-4.83481130337976D-05, 1.84907283946343D-06,-1.01418491178576D-07,
     4 7.05925634457153D-09,-5.85325291400382D-10, 5.56357688831339D-11,
     5-5.90889094779500D-12, 6.88574353784436D-13,-8.68588256452194D-14,
     6 1.17374762617213D-14,-1.68523146510923D-15, 2.55374773097056D-16/
      DATA DAJP(1), DAJP(2), DAJP(3), DAJP(4), DAJP(5), DAJP(6),
     1     DAJP(7), DAJP(8), DAJP(9), DAJP(10),DAJP(11),DAJP(12),
     2     DAJP(13),DAJP(14),DAJP(15),DAJP(16),DAJP(17),DAJP(18),
     3     DAJP(19)        / 6.53219131311457D-02,-1.20262933688823D-01,
     4 9.78010236263823D-03, 1.67948429230505D-02,-1.97146140182132D-03,
     5-8.45560295098867D-04, 9.42889620701976D-05, 2.25827860945475D-05,
     6-2.29067870915987D-06,-3.76343991136919D-07, 3.45663933559565D-08,
     7 4.29611332003007D-09,-3.58673691214989D-10,-3.57245881361895D-11,
     8 2.72696091066336D-12, 2.26120653095771D-13,-1.58763205238303D-14,
     9-1.12604374485125D-15, 7.31327529515367D-17/
      DATA DAJN(1), DAJN(2), DAJN(3), DAJN(4), DAJN(5), DAJN(6),
     1     DAJN(7), DAJN(8), DAJN(9), DAJN(10),DAJN(11),DAJN(12),
     2     DAJN(13),DAJN(14),DAJN(15),DAJN(16),DAJN(17),DAJN(18),
     3     DAJN(19)        / 1.08594539632967D-02, 8.53313194857091D-02,
     4-3.15277068113058D-01,-8.78420725294257D-02, 5.53251906976048D-02,
     5 9.41674060503241D-03,-3.32187026018996D-03,-4.11157343156826D-04,
     6 1.01297326891346D-04, 9.87633682208396D-06,-1.87312969812393D-06,
     7-1.50798500131468D-07, 2.32687669525394D-08, 1.59599917419225D-09,
     8-2.07665922668385D-10,-1.24103350500302D-11, 1.39631765331043D-12,
     9 7.39400971155740D-14,-7.32887475627500D-15/
      DATA DA(1),  DA(2),  DA(3),  DA(4),  DA(5),  DA(6),  DA(7),
     1     DA(8),  DA(9),  DA(10), DA(11), DA(12), DA(13), DA(14),
     2     DA(15)          / 4.91627321104601D-01, 3.11164930427489D-03,
     3 8.23140762854081D-05,-4.61769776172142D-06,-6.13158880534626D-08,
     4 2.87295804656520D-08,-1.81959715372117D-09,-1.44752826642035D-10,
     5 4.53724043420422D-11,-3.99655065847223D-12,-3.24089119830323D-13,
     6 1.62098952568741D-13,-2.40765247974057D-14, 1.69384811284491D-16,
     7 8.17900786477396D-16/
      DATA DB(1),  DB(2),  DB(3),  DB(4),  DB(5),  DB(6),  DB(7),
     1     DB(8),  DB(9),  DB(10), DB(11), DB(12), DB(13), DB(14),
     2     DB(15)          /-2.77571356944231D-01, 4.44212833419920D-03,
     3-8.42328522190089D-05,-2.58040318418710D-06, 3.42389720217621D-07,
     4-6.24286894709776D-09,-2.36377836844577D-09, 3.16991042656673D-10,
     5-4.40995691658191D-12,-5.18674221093575D-12, 9.64874015137022D-13,
     6-4.90190576608710D-14,-1.77253430678112D-14, 5.55950610442662D-15,
     7-7.11793337579530D-16/
C***FIRST EXECUTABLE STATEMENT  DJAIRY
      IF (X.LT.0.0D0) GO TO 90
      IF (C.GT.5.0D0) GO TO 60
      IF (X.GT.1.20D0) GO TO 30
      T = (X+X-1.2D0)*CON4
      TT = T + T
      J = N1
      F1 = AK1(J)
      F2 = 0.0D0
      DO 10 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK1(J)
        F2 = TEMP1
   10 CONTINUE
      AI = T*F1 - F2 + AK1(1)
C
      J = N1D
      F1 = DAK1(J)
      F2 = 0.0D0
      DO 20 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK1(J)
        F2 = TEMP1
   20 CONTINUE
      DAI = -(T*F1-F2+DAK1(1))
      RETURN
C
   30 CONTINUE
      T = (X+X-CON2)*CON3
      TT = T + T
      J = N2
      F1 = AK2(J)
      F2 = 0.0D0
      DO 40 I=1,M2
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK2(J)
        F2 = TEMP1
   40 CONTINUE
      RTRX = DSQRT(RX)
      EC = DEXP(-C)
      AI = EC*(T*F1-F2+AK2(1))/RTRX
      J = N2D
      F1 = DAK2(J)
      F2 = 0.0D0
      DO 50 I=1,M2D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK2(J)
        F2 = TEMP1
   50 CONTINUE
      DAI = -EC*(T*F1-F2+DAK2(1))*RTRX
      RETURN
C
   60 CONTINUE
      T = 10.0D0/C - 1.0D0
      TT = T + T
      J = N1
      F1 = AK3(J)
      F2 = 0.0D0
      DO 70 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK3(J)
        F2 = TEMP1
   70 CONTINUE
      RTRX = DSQRT(RX)
      EC = DEXP(-C)
      AI = EC*(T*F1-F2+AK3(1))/RTRX
      J = N1D
      F1 = DAK3(J)
      F2 = 0.0D0
      DO 80 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK3(J)
        F2 = TEMP1
   80 CONTINUE
      DAI = -RTRX*EC*(T*F1-F2+DAK3(1))
      RETURN
C
   90 CONTINUE
      IF (C.GT.5.0D0) GO TO 120
      T = 0.4D0*C - 1.0D0
      TT = T + T
      J = N3
      F1 = AJP(J)
      E1 = AJN(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 100 I=1,M3
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + AJP(J)
        E1 = TT*E1 - E2 + AJN(J)
        F2 = TEMP1
        E2 = TEMP2
  100 CONTINUE
      AI = (T*E1-E2+AJN(1)) - X*(T*F1-F2+AJP(1))
      J = N3D
      F1 = DAJP(J)
      E1 = DAJN(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 110 I=1,M3D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DAJP(J)
        E1 = TT*E1 - E2 + DAJN(J)
        F2 = TEMP1
        E2 = TEMP2
  110 CONTINUE
      DAI = X*X*(T*F1-F2+DAJP(1)) + (T*E1-E2+DAJN(1))
      RETURN
C
  120 CONTINUE
      T = 10.0D0/C - 1.0D0
      TT = T + T
      J = N4
      F1 = A(J)
      E1 = B(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 130 I=1,M4
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + A(J)
        E1 = TT*E1 - E2 + B(J)
        F2 = TEMP1
        E2 = TEMP2
  130 CONTINUE
      TEMP1 = T*F1 - F2 + A(1)
      TEMP2 = T*E1 - E2 + B(1)
      RTRX = DSQRT(RX)
      CV = C - FPI12
      CCV = DCOS(CV)
      SCV = DSIN(CV)
      AI = (TEMP1*CCV-TEMP2*SCV)/RTRX
      J = N4D
      F1 = DA(J)
      E1 = DB(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 140 I=1,M4D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DA(J)
        E1 = TT*E1 - E2 + DB(J)
        F2 = TEMP1
        E2 = TEMP2
  140 CONTINUE
      TEMP1 = T*F1 - F2 + DA(1)
      TEMP2 = T*E1 - E2 + DB(1)
      E1 = CCV*CON5 + 0.5D0*SCV
      E2 = SCV*CON5 - 0.5D0*CCV
      DAI = (TEMP1*E1-TEMP2*E2)*RTRX
      RETURN
      END
      DOUBLE PRECISION FUNCTION DLNGAM(X)
C***BEGIN PROLOGUE  DLNGAM
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A
C***KEYWORDS  ABSOLUTE VALUE,DOUBLE PRECISION,GAMMA FUNCTION,LOGARITHM,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. logarithm of the absolute value of the
C            Gamma function
C***DESCRIPTION
C
C DLNGAM(X) calculates the double precision logarithm of the
C absolute value of the gamma function for double precision
C argument X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9LGMC,DGAMMA,DINT,XERROR
C***END PROLOGUE  DLNGAM
      EXTERNAL DGAMMA
      DOUBLE PRECISION X, DXREL, PI, SINPIY, SQPI2L, SQ2PIL, XMAX,
     1  Y, DINT, DGAMMA, D9LGMC, D1MACH,TEMP
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA SQPI2L / +.2257913526 4472743236 3097614947 441 D+0    /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA XMAX, DXREL / 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DLNGAM
      IF (XMAX.NE.0.D0) GO TO 10
      TEMP = 1.0D0/DLOG(D1MACH(2))
      XMAX = TEMP * D1MACH(2)
      DXREL = DSQRT (D1MACH(4))
C
 10   Y = DABS (X)
      IF (Y.GT.10.D0) GO TO 20
C
C DLOG (DABS (DGAMMA(X)) ) FOR DABS(X) .LE. 10.0
C
      DLNGAM = DLOG (DABS (DGAMMA(X)) )
      RETURN
C
C DLOG ( DABS (DGAMMA(X)) ) FOR DABS(X) .GT. 10.0
C
 20   IF (Y.GT.XMAX) CALL XERROR ( 'DLNGAM  DABS(X) SO BIG DLNGAM OVERFL
     1OWS', 39, 2, 2)
C
      IF (X.GT.0.D0) DLNGAM = SQ2PIL + (X-0.5D0)*DLOG(X) - X + D9LGMC(Y)
      IF (X.GT.0.D0) RETURN
C
      SINPIY = DABS (DSIN(PI*Y))
      IF (SINPIY.EQ.0.D0) CALL XERROR ( 'DLNGAM  X IS A NEGATIVE INTEGER
     1', 31, 3, 2)
C
      IF (DABS ((X-DINT(X-0.5D0))/X).LT.DXREL) CALL XERROR ( 'DLNGAM  AN
     1SWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', 68, 1
     2, 1)
C
      DLNGAM = SQPI2L + (X-0.5D0)*DLOG(Y) - X - DLOG(SINPIY) - D9LGMC(Y)
      RETURN
C
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes an error (diagnostic) message.
C***DESCRIPTION
C     Abstract
C        XERROR processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed, containing
C                no more than 72 characters.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C
C     Examples
C        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
C        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
C                    43,2,1)
C        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
C    1ULLY COLLAPSED.',65,3,0)
C        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes error message allowing 2 integer and two real
C            values to be included in the message.
C***DESCRIPTION
C     Abstract
C        XERRWV processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C        In addition, up to two integer values and two real
C        values may be printed along with the message.
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C        NI    - number of integer values to be printed. (0 to 2)
C        I1    - first integer value.
C        I2    - second integer value.
C        NR    - number of real values to be printed. (0 to 2)
C        R1    - first real value.
C        R2    - second real value.
C
C     Examples
C        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
C    1   1,NUM,0,0,0.,0.)
C        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
C    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
C                    XGETUA
C***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
C     GET FLAGS
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Records that an error occurred.
C***DESCRIPTION
C     Abstract
C        Record that this error occurred.
C
C     Description of Parameters
C     --Input--
C       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
C       except that when NMESSG=0 the tables will be
C       dumped and cleared, and when NMESSG is less than zero the
C       tables will be dumped and not cleared.
C     --Output--
C       ICOUNT will be the number of times this message has
C       been seen, or zero if the table has overflowed and
C       does not contain this message specifically.
C       When NMESSG=0, ICOUNT will not be altered.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 Mar 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
C     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
C     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE XGETUA(IUNITA,N)
C***BEGIN PROLOGUE  XGETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Returns unit number(s) to which error messages are being
C            sent.
C***DESCRIPTION
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END



      SUBROUTINE DASYJY(FUNJY,X,FNU,FLGJY,IN,Y,WK,IFLW)
C***BEGIN PROLOGUE  DASYJY
C***REFER TO  DBESJ,DBESY
C***ROUTINES CALLED  D1MACH,I1MACH
C***DESCRIPTION
C
C                 DASYJY computes Bessel functions J and Y
C               for arguments X.GT.0.0 and orders FNU .GE. 35.0
C               on FLGJY = 1 and FLGJY = -1 respectively
C
C                                  INPUT
C
C      FUNJY - External function JAIRY or YAIRY
C          X - Argument, X.GT.0.0D0
C        FNU - Order of the first Bessel function
C      FLGJY - Selection flag
C              FLGJY =  1.0D0 gives the J function
C              FLGJY = -1.0D0 gives the Y function
C         IN - Number of functions desired, IN = 1 or 2
C
C                                  OUTPUT
C
C         Y  - A vector whose first IN components contain the sequence
C       IFLW - A flag indicating underflow or overflow
C                    return variables for BESJ only
C      WK(1) = 1 - (X/FNU)**2 = W**2
C      WK(2) = DSQRT(DABS(WK(1)))
C      WK(3) = DABS(WK(2) - DATAN(WK(2)))  or
C              DABS(LN((1 + WK(2))/(X/FNU)) - WK(2))
C            = DABS((2/3)*ZETA**(3/2))
C      WK(4) = FNU*WK(3)
C      WK(5) = (1.5*WK(3)*FNU)**(1/3) = DSQRT(ZETA)*FNU**(1/3)
C      WK(6) = DSIGN(1.,W**2)*WK(5)**2 = DSIGN(1.,W**2)*ZETA*FNU**(2/3)
C      WK(7) = FNU**(1/3)
C
C                                  Written by
C                                  D. E. Amos
C
C     Abstract   **** A Double Precision Routine ****
C         DASYJY implements the uniform asymptotic expansion of
C         the J and Y Bessel functions for FNU.GE.35 and real
C         X.GT.0.0D0. The forms are identical except for a change
C         in sign of some of the terms. This change in sign is
C         accomplished by means of the flag FLGJY = 1 or -1. On
C         FLGJY = 1 the Airy functions AI(X) and DAI(X) are
C         supplied by the external function JAIRY, and on
C         FLGJY = -1 the Airy functions BI(X) and DBI(X) are
C         supplied by the external funtion YAIRY.
C***END PROLOGUE  DASYJY
      INTEGER I, IFLW, IN, J, JN,JR,JU,K, KB,KLAST,KMAX,KP1, KS, KSP1,
     * KSTEMP, L, LR, LRP1, ISETA, ISETB
      INTEGER I1MACH
      DOUBLE PRECISION ABW2, AKM, ALFA, ALFA1, ALFA2, AP, AR, ASUM, AZ,
     * BETA, BETA1, BETA2, BETA3, BR, BSUM, C, CON1, CON2,
     * CON3,CON548,CR,CRZ32, DFI,ELIM, DR,FI, FLGJY, FN, FNU,
     * FN2, GAMA, PHI,  RCZ, RDEN, RELB, RFN2,  RTZ, RZDEN,
     * SA, SB, SUMA, SUMB, S1, TA, TAU, TB, TFN, TOL, TOLS, T2, UPOL,
     *  WK, X, XX, Y, Z, Z32
      DOUBLE PRECISION D1MACH
      DIMENSION Y(*), WK(*), C(65)
      DIMENSION ALFA(26,4), BETA(26,5)
      DIMENSION ALFA1(26,2), ALFA2(26,2)
      DIMENSION BETA1(26,2), BETA2(26,2), BETA3(26,1)
      DIMENSION GAMA(26), KMAX(5), AR(8), BR(10), UPOL(10)
      DIMENSION CR(10), DR(10)
      EQUIVALENCE (ALFA(1,1),ALFA1(1,1))
      EQUIVALENCE (ALFA(1,3),ALFA2(1,1))
      EQUIVALENCE (BETA(1,1),BETA1(1,1))
      EQUIVALENCE (BETA(1,3),BETA2(1,1))
      EQUIVALENCE (BETA(1,5),BETA3(1,1))
      SAVE TOLS,CON1, CON2, CON3, CON548, AR, BR, C,
     1 ALFA1, ALFA2, BETA1, BETA2, BETA3, GAMA
      DATA TOLS            /-6.90775527898214D+00/
      DATA CON1,CON2,CON3,CON548/
     1 6.66666666666667D-01, 3.33333333333333D-01, 1.41421356237310D+00,
     2 1.04166666666667D-01/
      DATA  AR(1),  AR(2),  AR(3),  AR(4),  AR(5),  AR(6),  AR(7),
     A      AR(8)          / 8.35503472222222D-02, 1.28226574556327D-01,
     1 2.91849026464140D-01, 8.81627267443758D-01, 3.32140828186277D+00,
     2 1.49957629868626D+01, 7.89230130115865D+01, 4.74451538868264D+02/
      DATA  BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     A      BR(9), BR(10)  /-1.45833333333333D-01,-9.87413194444444D-02,
     1-1.43312053915895D-01,-3.17227202678414D-01,-9.42429147957120D-01,
     2-3.51120304082635D+00,-1.57272636203680D+01,-8.22814390971859D+01,
     3-4.92355370523671D+02,-3.31621856854797D+03/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3       -2.08333333333333D-01,        1.25000000000000D-01,
     4        3.34201388888889D-01,       -4.01041666666667D-01,
     5        7.03125000000000D-02,       -1.02581259645062D+00,
     6        1.84646267361111D+00,       -8.91210937500000D-01,
     7        7.32421875000000D-02,        4.66958442342625D+00,
     8       -1.12070026162230D+01,        8.78912353515625D+00,
     9       -2.36408691406250D+00,        1.12152099609375D-01,
     A       -2.82120725582002D+01,        8.46362176746007D+01,
     B       -9.18182415432400D+01,        4.25349987453885D+01,
     C       -7.36879435947963D+00,        2.27108001708984D-01,
     D        2.12570130039217D+02,       -7.65252468141182D+02,
     E        1.05999045252800D+03,       -6.99579627376133D+02/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3        2.18190511744212D+02,       -2.64914304869516D+01,
     4        5.72501420974731D-01,       -1.91945766231841D+03,
     5        8.06172218173731D+03,       -1.35865500064341D+04,
     6        1.16553933368645D+04,       -5.30564697861340D+03,
     7        1.20090291321635D+03,       -1.08090919788395D+02,
     8        1.72772750258446D+00,        2.02042913309661D+04,
     9       -9.69805983886375D+04,        1.92547001232532D+05,
     A       -2.03400177280416D+05,        1.22200464983017D+05,
     B       -4.11926549688976D+04,        7.10951430248936D+03,
     C       -4.93915304773088D+02,        6.07404200127348D+00,
     D       -2.42919187900551D+05,        1.31176361466298D+06,
     E       -2.99801591853811D+06,        3.76327129765640D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65)/
     3       -2.81356322658653D+06,        1.26836527332162D+06,
     4       -3.31645172484564D+05,        4.52187689813627D+04,
     5       -2.49983048181121D+03,        2.43805296995561D+01,
     6        3.28446985307204D+06,       -1.97068191184322D+07,
     7        5.09526024926646D+07,       -7.41051482115327D+07,
     8        6.63445122747290D+07,       -3.75671766607634D+07,
     9        1.32887671664218D+07,       -2.78561812808645D+06,
     A        3.08186404612662D+05,       -1.38860897537170D+04,
     B        1.10017140269247D+02/
      DATA ALFA1(1,1), ALFA1(2,1), ALFA1(3,1), ALFA1(4,1), ALFA1(5,1),
     1     ALFA1(6,1), ALFA1(7,1), ALFA1(8,1), ALFA1(9,1), ALFA1(10,1),
     2     ALFA1(11,1),ALFA1(12,1),ALFA1(13,1),ALFA1(14,1),ALFA1(15,1),
     3     ALFA1(16,1),ALFA1(17,1),ALFA1(18,1),ALFA1(19,1),ALFA1(20,1),
     4     ALFA1(21,1),ALFA1(22,1),ALFA1(23,1),ALFA1(24,1),ALFA1(25,1),
     5     ALFA1(26,1)     /-4.44444444444444D-03,-9.22077922077922D-04,
     6-8.84892884892885D-05, 1.65927687832450D-04, 2.46691372741793D-04,
     7 2.65995589346255D-04, 2.61824297061501D-04, 2.48730437344656D-04,
     8 2.32721040083232D-04, 2.16362485712365D-04, 2.00738858762752D-04,
     9 1.86267636637545D-04, 1.73060775917876D-04, 1.61091705929016D-04,
     1 1.50274774160908D-04, 1.40503497391270D-04, 1.31668816545923D-04,
     2 1.23667445598253D-04, 1.16405271474738D-04, 1.09798298372713D-04,
     3 1.03772410422993D-04, 9.82626078369363D-05, 9.32120517249503D-05,
     4 8.85710852478712D-05, 8.42963105715700D-05, 8.03497548407791D-05/
      DATA ALFA1(1,2), ALFA1(2,2), ALFA1(3,2), ALFA1(4,2), ALFA1(5,2),
     1     ALFA1(6,2), ALFA1(7,2), ALFA1(8,2), ALFA1(9,2), ALFA1(10,2),
     2     ALFA1(11,2),ALFA1(12,2),ALFA1(13,2),ALFA1(14,2),ALFA1(15,2),
     3     ALFA1(16,2),ALFA1(17,2),ALFA1(18,2),ALFA1(19,2),ALFA1(20,2),
     4     ALFA1(21,2),ALFA1(22,2),ALFA1(23,2),ALFA1(24,2),ALFA1(25,2),
     5     ALFA1(26,2)     / 6.93735541354589D-04, 2.32241745182922D-04,
     6-1.41986273556691D-05,-1.16444931672049D-04,-1.50803558053049D-04,
     7-1.55121924918096D-04,-1.46809756646466D-04,-1.33815503867491D-04,
     8-1.19744975684254D-04,-1.06184319207974D-04,-9.37699549891194D-05,
     9-8.26923045588193D-05,-7.29374348155221D-05,-6.44042357721016D-05,
     1-5.69611566009369D-05,-5.04731044303562D-05,-4.48134868008883D-05,
     2-3.98688727717599D-05,-3.55400532972042D-05,-3.17414256609022D-05,
     3-2.83996793904175D-05,-2.54522720634871D-05,-2.28459297164725D-05,
     4-2.05352753106481D-05,-1.84816217627666D-05,-1.66519330021394D-05/
      DATA ALFA2(1,1), ALFA2(2,1), ALFA2(3,1), ALFA2(4,1), ALFA2(5,1),
     1     ALFA2(6,1), ALFA2(7,1), ALFA2(8,1), ALFA2(9,1), ALFA2(10,1),
     2     ALFA2(11,1),ALFA2(12,1),ALFA2(13,1),ALFA2(14,1),ALFA2(15,1),
     3     ALFA2(16,1),ALFA2(17,1),ALFA2(18,1),ALFA2(19,1),ALFA2(20,1),
     4     ALFA2(21,1),ALFA2(22,1),ALFA2(23,1),ALFA2(24,1),ALFA2(25,1),
     5     ALFA2(26,1)     /-3.54211971457744D-04,-1.56161263945159D-04,
     6 3.04465503594936D-05, 1.30198655773243D-04, 1.67471106699712D-04,
     7 1.70222587683593D-04, 1.56501427608595D-04, 1.36339170977445D-04,
     8 1.14886692029825D-04, 9.45869093034688D-05, 7.64498419250898D-05,
     9 6.07570334965197D-05, 4.74394299290509D-05, 3.62757512005344D-05,
     1 2.69939714979225D-05, 1.93210938247939D-05, 1.30056674793963D-05,
     2 7.82620866744497D-06, 3.59257485819352D-06, 1.44040049814252D-07,
     3-2.65396769697939D-06,-4.91346867098486D-06,-6.72739296091248D-06,
     4-8.17269379678658D-06,-9.31304715093561D-06,-1.02011418798016D-05/
      DATA ALFA2(1,2), ALFA2(2,2), ALFA2(3,2), ALFA2(4,2), ALFA2(5,2),
     1     ALFA2(6,2), ALFA2(7,2), ALFA2(8,2), ALFA2(9,2), ALFA2(10,2),
     2     ALFA2(11,2),ALFA2(12,2),ALFA2(13,2),ALFA2(14,2),ALFA2(15,2),
     3     ALFA2(16,2),ALFA2(17,2),ALFA2(18,2),ALFA2(19,2),ALFA2(20,2),
     4     ALFA2(21,2),ALFA2(22,2),ALFA2(23,2),ALFA2(24,2),ALFA2(25,2),
     5     ALFA2(26,2)     / 3.78194199201773D-04, 2.02471952761816D-04,
     6-6.37938506318862D-05,-2.38598230603006D-04,-3.10916256027362D-04,
     7-3.13680115247576D-04,-2.78950273791323D-04,-2.28564082619141D-04,
     8-1.75245280340847D-04,-1.25544063060690D-04,-8.22982872820208D-05,
     9-4.62860730588116D-05,-1.72334302366962D-05, 5.60690482304602D-06,
     1 2.31395443148287D-05, 3.62642745856794D-05, 4.58006124490189D-05,
     2 5.24595294959114D-05, 5.68396208545815D-05, 5.94349820393104D-05,
     3 6.06478527578422D-05, 6.08023907788436D-05, 6.01577894539460D-05,
     4 5.89199657344698D-05, 5.72515823777593D-05, 5.52804375585853D-05/
      DATA BETA1(1,1), BETA1(2,1), BETA1(3,1), BETA1(4,1), BETA1(5,1),
     1     BETA1(6,1), BETA1(7,1), BETA1(8,1), BETA1(9,1), BETA1(10,1),
     2     BETA1(11,1),BETA1(12,1),BETA1(13,1),BETA1(14,1),BETA1(15,1),
     3     BETA1(16,1),BETA1(17,1),BETA1(18,1),BETA1(19,1),BETA1(20,1),
     4     BETA1(21,1),BETA1(22,1),BETA1(23,1),BETA1(24,1),BETA1(25,1),
     5     BETA1(26,1)     / 1.79988721413553D-02, 5.59964911064388D-03,
     6 2.88501402231133D-03, 1.80096606761054D-03, 1.24753110589199D-03,
     7 9.22878876572938D-04, 7.14430421727287D-04, 5.71787281789705D-04,
     8 4.69431007606482D-04, 3.93232835462917D-04, 3.34818889318298D-04,
     9 2.88952148495752D-04, 2.52211615549573D-04, 2.22280580798883D-04,
     1 1.97541838033063D-04, 1.76836855019718D-04, 1.59316899661821D-04,
     2 1.44347930197334D-04, 1.31448068119965D-04, 1.20245444949303D-04,
     3 1.10449144504599D-04, 1.01828770740567D-04, 9.41998224204238D-05,
     4 8.74130545753834D-05, 8.13466262162801D-05, 7.59002269646219D-05/
      DATA BETA1(1,2), BETA1(2,2), BETA1(3,2), BETA1(4,2), BETA1(5,2),
     1     BETA1(6,2), BETA1(7,2), BETA1(8,2), BETA1(9,2), BETA1(10,2),
     2     BETA1(11,2),BETA1(12,2),BETA1(13,2),BETA1(14,2),BETA1(15,2),
     3     BETA1(16,2),BETA1(17,2),BETA1(18,2),BETA1(19,2),BETA1(20,2),
     4     BETA1(21,2),BETA1(22,2),BETA1(23,2),BETA1(24,2),BETA1(25,2),
     5     BETA1(26,2)     /-1.49282953213429D-03,-8.78204709546389D-04,
     6-5.02916549572035D-04,-2.94822138512746D-04,-1.75463996970783D-04,
     7-1.04008550460816D-04,-5.96141953046458D-05,-3.12038929076098D-05,
     8-1.26089735980230D-05,-2.42892608575730D-07, 8.05996165414274D-06,
     9 1.36507009262147D-05, 1.73964125472926D-05, 1.98672978842134D-05,
     1 2.14463263790823D-05, 2.23954659232457D-05, 2.28967783814713D-05,
     2 2.30785389811178D-05, 2.30321976080909D-05, 2.28236073720349D-05,
     3 2.25005881105292D-05, 2.20981015361991D-05, 2.16418427448104D-05,
     4 2.11507649256221D-05, 2.06388749782171D-05, 2.01165241997082D-05/
      DATA BETA2(1,1), BETA2(2,1), BETA2(3,1), BETA2(4,1), BETA2(5,1),
     1     BETA2(6,1), BETA2(7,1), BETA2(8,1), BETA2(9,1), BETA2(10,1),
     2     BETA2(11,1),BETA2(12,1),BETA2(13,1),BETA2(14,1),BETA2(15,1),
     3     BETA2(16,1),BETA2(17,1),BETA2(18,1),BETA2(19,1),BETA2(20,1),
     4     BETA2(21,1),BETA2(22,1),BETA2(23,1),BETA2(24,1),BETA2(25,1),
     5     BETA2(26,1)     / 5.52213076721293D-04, 4.47932581552385D-04,
     6 2.79520653992021D-04, 1.52468156198447D-04, 6.93271105657044D-05,
     7 1.76258683069991D-05,-1.35744996343269D-05,-3.17972413350427D-05,
     8-4.18861861696693D-05,-4.69004889379141D-05,-4.87665447413787D-05,
     9-4.87010031186735D-05,-4.74755620890087D-05,-4.55813058138628D-05,
     1-4.33309644511266D-05,-4.09230193157750D-05,-3.84822638603221D-05,
     2-3.60857167535411D-05,-3.37793306123367D-05,-3.15888560772110D-05,
     3-2.95269561750807D-05,-2.75978914828336D-05,-2.58006174666884D-05,
     4-2.41308356761280D-05,-2.25823509518346D-05,-2.11479656768913D-05/
      DATA BETA2(1,2), BETA2(2,2), BETA2(3,2), BETA2(4,2), BETA2(5,2),
     1     BETA2(6,2), BETA2(7,2), BETA2(8,2), BETA2(9,2), BETA2(10,2),
     2     BETA2(11,2),BETA2(12,2),BETA2(13,2),BETA2(14,2),BETA2(15,2),
     3     BETA2(16,2),BETA2(17,2),BETA2(18,2),BETA2(19,2),BETA2(20,2),
     4     BETA2(21,2),BETA2(22,2),BETA2(23,2),BETA2(24,2),BETA2(25,2),
     5     BETA2(26,2)     /-4.74617796559960D-04,-4.77864567147321D-04,
     6-3.20390228067038D-04,-1.61105016119962D-04,-4.25778101285435D-05,
     7 3.44571294294968D-05, 7.97092684075675D-05, 1.03138236708272D-04,
     8 1.12466775262204D-04, 1.13103642108481D-04, 1.08651634848774D-04,
     9 1.01437951597662D-04, 9.29298396593364D-05, 8.40293133016090D-05,
     1 7.52727991349134D-05, 6.69632521975731D-05, 5.92564547323195D-05,
     2 5.22169308826976D-05, 4.58539485165361D-05, 4.01445513891487D-05,
     3 3.50481730031328D-05, 3.05157995034347D-05, 2.64956119950516D-05,
     4 2.29363633690998D-05, 1.97893056664022D-05, 1.70091984636413D-05/
      DATA BETA3(1,1), BETA3(2,1), BETA3(3,1), BETA3(4,1), BETA3(5,1),
     1     BETA3(6,1), BETA3(7,1), BETA3(8,1), BETA3(9,1), BETA3(10,1),
     2     BETA3(11,1),BETA3(12,1),BETA3(13,1),BETA3(14,1),BETA3(15,1),
     3     BETA3(16,1),BETA3(17,1),BETA3(18,1),BETA3(19,1),BETA3(20,1),
     4     BETA3(21,1),BETA3(22,1),BETA3(23,1),BETA3(24,1),BETA3(25,1),
     5     BETA3(26,1)     / 7.36465810572578D-04, 8.72790805146194D-04,
     6 6.22614862573135D-04, 2.85998154194304D-04, 3.84737672879366D-06,
     7-1.87906003636972D-04,-2.97603646594555D-04,-3.45998126832656D-04,
     8-3.53382470916038D-04,-3.35715635775049D-04,-3.04321124789040D-04,
     9-2.66722723047613D-04,-2.27654214122820D-04,-1.89922611854562D-04,
     1-1.55058918599094D-04,-1.23778240761874D-04,-9.62926147717644D-05,
     2-7.25178327714425D-05,-5.22070028895634D-05,-3.50347750511901D-05,
     3-2.06489761035552D-05,-8.70106096849767D-06, 1.13698686675100D-06,
     4 9.16426474122779D-06, 1.56477785428873D-05, 2.08223629482467D-05/
      DATA GAMA(1),   GAMA(2),   GAMA(3),   GAMA(4),   GAMA(5),
     1     GAMA(6),   GAMA(7),   GAMA(8),   GAMA(9),   GAMA(10),
     2     GAMA(11),  GAMA(12),  GAMA(13),  GAMA(14),  GAMA(15),
     3     GAMA(16),  GAMA(17),  GAMA(18),  GAMA(19),  GAMA(20),
     4     GAMA(21),  GAMA(22),  GAMA(23),  GAMA(24),  GAMA(25),
     5     GAMA(26)        / 6.29960524947437D-01, 2.51984209978975D-01,
     6 1.54790300415656D-01, 1.10713062416159D-01, 8.57309395527395D-02,
     7 6.97161316958684D-02, 5.86085671893714D-02, 5.04698873536311D-02,
     8 4.42600580689155D-02, 3.93720661543510D-02, 3.54283195924455D-02,
     9 3.21818857502098D-02, 2.94646240791158D-02, 2.71581677112934D-02,
     1 2.51768272973862D-02, 2.34570755306079D-02, 2.19508390134907D-02,
     2 2.06210828235646D-02, 1.94388240897881D-02, 1.83810633800683D-02,
     3 1.74293213231963D-02, 1.65685837786612D-02, 1.57865285987918D-02,
     4 1.50729501494096D-02, 1.44193250839955D-02, 1.38184805735342D-02/
C***FIRST EXECUTABLE STATEMENT  DASYJY
      TA = D1MACH(3)
      TOL = DMAX1(TA,1.0D-15)
      TB = D1MACH(5)
      JU = I1MACH(15)
      IF(FLGJY.EQ.1.0D0) GO TO 6
      JR = I1MACH(14)
      ELIM = 2.303D0*TB*(DBLE(FLOAT(-JU))-DBLE(FLOAT(JR)))
      GO TO 7
    6 CONTINUE
      ELIM = 2.303D0*(TB*DBLE(FLOAT(-JU))-3.0D0)
    7 CONTINUE
      FN = FNU
      IFLW = 0
      DO 170 JN=1,IN
        XX = X/FN
        WK(1) = 1.0D0 - XX*XX
        ABW2 = DABS(WK(1))
        WK(2) = DSQRT(ABW2)
        WK(7) = FN**CON2
        IF (ABW2.GT.0.27750D0) GO TO 80
C
C     ASYMPTOTIC EXPANSION
C     CASES NEAR X=FN, DABS(1.-(X/FN)**2).LE.0.2775
C     COEFFICIENTS OF ASYMPTOTIC EXPANSION BY SERIES
C
C     ZETA AND TRUNCATION FOR A(ZETA) AND B(ZETA) SERIES
C
C     KMAX IS TRUNCATION INDEX FOR A(ZETA) AND B(ZETA) SERIES=MAX(2,SA)
C
        SA = 0.0D0
        IF (ABW2.EQ.0.0D0) GO TO 10
        SA = TOLS/DLOG(ABW2)
   10   SB = SA
        DO 20 I=1,5
          AKM = DMAX1(SA,2.0D0)
          KMAX(I) = INT(SNGL(AKM))
          SA = SA + SB
   20   CONTINUE
        KB = KMAX(5)
        KLAST = KB - 1
        SA = GAMA(KB)
        DO 30 K=1,KLAST
          KB = KB - 1
          SA = SA*WK(1) + GAMA(KB)
   30   CONTINUE
        Z = WK(1)*SA
        AZ = DABS(Z)
        RTZ = DSQRT(AZ)
        WK(3) = CON1*AZ*RTZ
        WK(4) = WK(3)*FN
        WK(5) = RTZ*WK(7)
        WK(6) = -WK(5)*WK(5)
        IF(Z.LE.0.0D0) GO TO 35
        IF(WK(4).GT.ELIM) GO TO 75
        WK(6) = -WK(6)
   35   CONTINUE
        PHI = DSQRT(DSQRT(SA+SA+SA+SA))
C
C     B(ZETA) FOR S=0
C
        KB = KMAX(5)
        KLAST = KB - 1
        SB = BETA(KB,1)
        DO 40 K=1,KLAST
          KB = KB - 1
          SB = SB*WK(1) + BETA(KB,1)
   40   CONTINUE
        KSP1 = 1
        FN2 = FN*FN
        RFN2 = 1.0D0/FN2
        RDEN = 1.0D0
        ASUM = 1.0D0
        RELB = TOL*DABS(SB)
        BSUM = SB
        DO 60 KS=1,4
          KSP1 = KSP1 + 1
          RDEN = RDEN*RFN2
C
C     A(ZETA) AND B(ZETA) FOR S=1,2,3,4
C
          KSTEMP = 5 - KS
          KB = KMAX(KSTEMP)
          KLAST = KB - 1
          SA = ALFA(KB,KS)
          SB = BETA(KB,KSP1)
          DO 50 K=1,KLAST
            KB = KB - 1
            SA = SA*WK(1) + ALFA(KB,KS)
            SB = SB*WK(1) + BETA(KB,KSP1)
   50     CONTINUE
          TA = SA*RDEN
          TB = SB*RDEN
          ASUM = ASUM + TA
          BSUM = BSUM + TB
          IF (DABS(TA).LE.TOL .AND. DABS(TB).LE.RELB) GO TO 70
   60   CONTINUE
   70   CONTINUE
        BSUM = BSUM/(FN*WK(7))
        GO TO 160
C
   75   CONTINUE
        IFLW = 1
        RETURN
C
   80   CONTINUE
        UPOL(1) = 1.0D0
        TAU = 1.0D0/WK(2)
        T2 = 1.0D0/WK(1)
        IF (WK(1).GE.0.0D0) GO TO 90
C
C     CASES FOR (X/FN).GT.DSQRT(1.2775)
C
        WK(3) = DABS(WK(2)-DATAN(WK(2)))
        WK(4) = WK(3)*FN
        RCZ = -CON1/WK(4)
        Z32 = 1.5D0*WK(3)
        RTZ = Z32**CON2
        WK(5) = RTZ*WK(7)
        WK(6) = -WK(5)*WK(5)
        GO TO 100
   90   CONTINUE
C
C     CASES FOR (X/FN).LT.DSQRT(0.7225)
C
        WK(3) = DABS(DLOG((1.0D0+WK(2))/XX)-WK(2))
        WK(4) = WK(3)*FN
        RCZ = CON1/WK(4)
        IF(WK(4).GT.ELIM) GO TO 75
        Z32 = 1.5D0*WK(3)
        RTZ = Z32**CON2
        WK(7) = FN**CON2
        WK(5) = RTZ*WK(7)
        WK(6) = WK(5)*WK(5)
  100   CONTINUE
        PHI = DSQRT((RTZ+RTZ)*TAU)
        TB = 1.0D0
        ASUM = 1.0D0
        TFN = TAU/FN
        RDEN=1.0D0/FN
        RFN2=RDEN*RDEN
        RDEN=1.0D0
        UPOL(2) = (C(1)*T2+C(2))*TFN
        CRZ32 = CON548*RCZ
        BSUM = UPOL(2) + CRZ32
        RELB = TOL*DABS(BSUM)
        AP = TFN
        KS = 0
        KP1 = 2
        RZDEN = RCZ
        L = 2
        ISETA=0
        ISETB=0
        DO 140 LR=2,8,2
C
C     COMPUTE TWO U POLYNOMIALS FOR NEXT A(ZETA) AND B(ZETA)
C
          LRP1 = LR + 1
          DO 120 K=LR,LRP1
            KS = KS + 1
            KP1 = KP1 + 1
            L = L + 1
            S1 = C(L)
            DO 110 J=2,KP1
              L = L + 1
              S1 = S1*T2 + C(L)
  110       CONTINUE
            AP = AP*TFN
            UPOL(KP1) = AP*S1
            CR(KS) = BR(KS)*RZDEN
            RZDEN = RZDEN*RCZ
            DR(KS) = AR(KS)*RZDEN
  120     CONTINUE
          SUMA = UPOL(LRP1)
          SUMB = UPOL(LR+2) + UPOL(LRP1)*CRZ32
          JU = LRP1
          DO 130 JR=1,LR
            JU = JU - 1
            SUMA = SUMA + CR(JR)*UPOL(JU)
            SUMB = SUMB + DR(JR)*UPOL(JU)
  130     CONTINUE
          RDEN=RDEN*RFN2
          TB = -TB
          IF (WK(1).GT.0.0D0) TB = DABS(TB)
          IF(RDEN.LT.TOL) GO TO 131
          ASUM = ASUM + SUMA*TB
          BSUM = BSUM + SUMB*TB
          GO TO 140
  131     IF(ISETA.EQ.1) GO TO 132
          IF(DABS(SUMA).LT.TOL) ISETA=1
          ASUM=ASUM+SUMA*TB
  132     IF(ISETB.EQ.1) GO TO 133
          IF(DABS(SUMB).LT.RELB) ISETB=1
          BSUM=BSUM+SUMB*TB
  133     IF(ISETA.EQ.1 .AND. ISETB.EQ.1) GO TO 150
  140   CONTINUE
  150   TB = WK(5)
        IF (WK(1).GT.0.0D0) TB = -TB
        BSUM = BSUM/TB
C
  160   CONTINUE
        CALL FUNJY(WK(6), WK(5), WK(4), FI, DFI)
        TA=1.0D0/TOL
        TB=D1MACH(1)*TA*1.0D+3
        IF(DABS(FI).GT.TB) GO TO 165
        FI=FI*TA
        DFI=DFI*TA
        PHI=PHI*TOL
  165   CONTINUE
        Y(JN) = FLGJY*PHI*(FI*ASUM+DFI*BSUM)/WK(7)
        FN = FN - FLGJY
  170 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION D9LGMC(X)
C***BEGIN PROLOGUE  D9LGMC
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  COMPLETE GAMMA FUNCTION,CORRECTION FACTOR,
C             DOUBLE PRECISION,GAMMA FUNCTION,LOGARITHM,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the  d.p. log Gamma correction factor for
C            X .GE. 10. so that DLOG(DGAMMA(X)) = DLOG(DSQRT(2*PI)) +
C            (X-5.)*DLOG(X) - X + D9LGMC(X)
C***DESCRIPTION
C
C Compute the log gamma correction factor for X .GE. 10. so that
C DLOG (DGAMMA(X)) = DLOG(DSQRT(2*PI)) + (X-.5)*DLOG(X) - X + D9lGMC(X)
C
C Series for ALGM       on the interval  0.          to  1.00000E-02
C                                        with weighted error   1.28E-31
C                                         log weighted error  30.89
C                               significant figures required  29.81
C                                    decimal places required  31.48
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DCSEVL,INITDS,XERROR
C***END PROLOGUE  D9LGMC
      DOUBLE PRECISION X, ALGMCS(15), XBIG, XMAX, DCSEVL, D1MACH
      DATA ALGMCS(  1) / +.1666389480 4518632472 0572965082 2 D+0      /
      DATA ALGMCS(  2) / -.1384948176 0675638407 3298605913 5 D-4      /
      DATA ALGMCS(  3) / +.9810825646 9247294261 5717154748 7 D-8      /
      DATA ALGMCS(  4) / -.1809129475 5724941942 6330626671 9 D-10     /
      DATA ALGMCS(  5) / +.6221098041 8926052271 2601554341 6 D-13     /
      DATA ALGMCS(  6) / -.3399615005 4177219443 0333059966 6 D-15     /
      DATA ALGMCS(  7) / +.2683181998 4826987489 5753884666 6 D-17     /
      DATA ALGMCS(  8) / -.2868042435 3346432841 4462239999 9 D-19     /
      DATA ALGMCS(  9) / +.3962837061 0464348036 7930666666 6 D-21     /
      DATA ALGMCS( 10) / -.6831888753 9857668701 1199999999 9 D-23     /
      DATA ALGMCS( 11) / +.1429227355 9424981475 7333333333 3 D-24     /
      DATA ALGMCS( 12) / -.3547598158 1010705471 9999999999 9 D-26     /
      DATA ALGMCS( 13) / +.1025680058 0104709120 0000000000 0 D-27     /
      DATA ALGMCS( 14) / -.3401102254 3167487999 9999999999 9 D-29     /
      DATA ALGMCS( 15) / +.1276642195 6300629333 3333333333 3 D-30     /
      DATA NALGM, XBIG, XMAX / 0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  D9LGMC
      IF (NALGM.NE.0) GO TO 10
      NALGM = INITDS (ALGMCS, 15, SNGL(D1MACH(3)) )
      XBIG = 1.0D0/DSQRT(D1MACH(3))
      XMAX = DEXP (DMIN1(DLOG(D1MACH(2)/12.D0), -DLOG(12.D0*D1MACH(1))))
C
 10   IF (X.LT.10.D0) CALL XERROR ( 'D9LGMC  X MUST BE GE 10', 23, 1, 2)
      IF (X.GE.XMAX) GO TO 20
C
      D9LGMC = 1.D0/(12.D0*X)
      IF (X.LT.XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS,
     1  NALGM) / X
      RETURN
C
 20   D9LGMC = 0.D0
      CALL XERROR ( 'D9LGMC  X SO BIG D9LGMC UNDERFLOWS', 34, 2, 1)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION DCSEVL(X,A,N)
C***BEGIN PROLOGUE  DCSEVL
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C3A2
C***KEYWORDS  CHEBYSHEV,FNLIB,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluate the double precision N-term Chebyshev series A
C            at X.
C***DESCRIPTION
C
C Evaluate the N-term Chebyshev series A at X.  Adapted from
C R. Broucke, Algorithm 446, C.A.C.M., 16, 254 (1973).
C W. Fullerton, C-3, Los Alamos Scientific Laboratory.
C
C       Input Arguments --
C X    double precision value at which the series is to be evaluated.
C A    double precision array of N terms of a Chebyshev series.  In
C      evaluating A, only half of the first coefficient is summed.
C N    number of terms in array A.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  DCSEVL
C
       DOUBLE PRECISION A(N),X,TWOX,B0,B1,B2
C***FIRST EXECUTABLE STATEMENT  DCSEVL
       IF(N.LT.1)CALL XERROR( 'DCSEVL  NUMBER OF TERMS LE 0', 28, 2,2)
       IF(N.GT.1000) CALL XERROR ( 'DCSEVL  NUMBER OF TERMS GT 1000',
     1   31, 3, 2)
       IF ((X.LT.-1.D0) .OR. (X.GT.1.D0)) CALL XERROR ( 'DCSEVL  X OUTSI
     1DE (-1,+1)', 25, 1, 1)
C
       TWOX = 2.0D0*X
       B1 = 0.D0
       B0=0.D0
       DO 10 I=1,N
         B2=B1
         B1=B0
         NI = N - I + 1
         B0 = TWOX*B1 - B2 + A(NI)
 10    CONTINUE
C
       DCSEVL = 0.5D0 * (B0-B2)
C
       RETURN
      END
      DOUBLE PRECISION FUNCTION DGAMMA(X)
C***BEGIN PROLOGUE  DGAMMA
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A
C***KEYWORDS  COMPLETE GAMMA FUNCTION,DOUBLE PRECISION,GAMMA FUNCTION,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. complete Gamma function.
C***DESCRIPTION
C
C DGAMMA(X) calculates the double precision complete gamma function
C for double precision argument X.
C
C Series for GAM        on the interval  0.          to  1.00000E+00
C                                        with weighted error   5.79E-32
C                                         log weighted error  31.24
C                               significant figures required  30.00
C                                    decimal places required  32.05
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9LGMC,DCSEVL,DGAMLM,DINT,INITDS,XERROR
C***END PROLOGUE  DGAMMA
      DOUBLE PRECISION X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX,
     1  XMIN, Y,  DINT, D9LGMC, DCSEVL, D1MACH
C
      DATA GAM CS(  1) / +.8571195590 9893314219 2006239994 2 D-2      /
      DATA GAM CS(  2) / +.4415381324 8410067571 9131577165 2 D-2      /
      DATA GAM CS(  3) / +.5685043681 5993633786 3266458878 9 D-1      /
      DATA GAM CS(  4) / -.4219835396 4185605010 1250018662 4 D-2      /
      DATA GAM CS(  5) / +.1326808181 2124602205 8400679635 2 D-2      /
      DATA GAM CS(  6) / -.1893024529 7988804325 2394702388 6 D-3      /
      DATA GAM CS(  7) / +.3606925327 4412452565 7808221722 5 D-4      /
      DATA GAM CS(  8) / -.6056761904 4608642184 8554829036 5 D-5      /
      DATA GAM CS(  9) / +.1055829546 3022833447 3182350909 3 D-5      /
      DATA GAM CS( 10) / -.1811967365 5423840482 9185589116 6 D-6      /
      DATA GAM CS( 11) / +.3117724964 7153222777 9025459316 9 D-7      /
      DATA GAM CS( 12) / -.5354219639 0196871408 7408102434 7 D-8      /
      DATA GAM CS( 13) / +.9193275519 8595889468 8778682594 0 D-9      /
      DATA GAM CS( 14) / -.1577941280 2883397617 6742327395 3 D-9      /
      DATA GAM CS( 15) / +.2707980622 9349545432 6654043308 9 D-10     /
      DATA GAM CS( 16) / -.4646818653 8257301440 8166105893 3 D-11     /
      DATA GAM CS( 17) / +.7973350192 0074196564 6076717535 9 D-12     /
      DATA GAM CS( 18) / -.1368078209 8309160257 9949917230 9 D-12     /
      DATA GAM CS( 19) / +.2347319486 5638006572 3347177168 8 D-13     /
      DATA GAM CS( 20) / -.4027432614 9490669327 6657053469 9 D-14     /
      DATA GAM CS( 21) / +.6910051747 3721009121 3833697525 7 D-15     /
      DATA GAM CS( 22) / -.1185584500 2219929070 5238712619 2 D-15     /
      DATA GAM CS( 23) / +.2034148542 4963739552 0102605193 2 D-16     /
      DATA GAM CS( 24) / -.3490054341 7174058492 7401294910 8 D-17     /
      DATA GAM CS( 25) / +.5987993856 4853055671 3505106602 6 D-18     /
      DATA GAM CS( 26) / -.1027378057 8722280744 9006977843 1 D-18     /
      DATA GAM CS( 27) / +.1762702816 0605298249 4275966074 8 D-19     /
      DATA GAM CS( 28) / -.3024320653 7353062609 5877211204 2 D-20     /
      DATA GAM CS( 29) / +.5188914660 2183978397 1783355050 6 D-21     /
      DATA GAM CS( 30) / -.8902770842 4565766924 4925160106 6 D-22     /
      DATA GAM CS( 31) / +.1527474068 4933426022 7459689130 6 D-22     /
      DATA GAM CS( 32) / -.2620731256 1873629002 5732833279 9 D-23     /
      DATA GAM CS( 33) / +.4496464047 8305386703 3104657066 6 D-24     /
      DATA GAM CS( 34) / -.7714712731 3368779117 0390152533 3 D-25     /
      DATA GAM CS( 35) / +.1323635453 1260440364 8657271466 6 D-25     /
      DATA GAM CS( 36) / -.2270999412 9429288167 0231381333 3 D-26     /
      DATA GAM CS( 37) / +.3896418998 0039914493 2081663999 9 D-27     /
      DATA GAM CS( 38) / -.6685198115 1259533277 9212799999 9 D-28     /
      DATA GAM CS( 39) / +.1146998663 1400243843 4761386666 6 D-28     /
      DATA GAM CS( 40) / -.1967938586 3451346772 9510399999 9 D-29     /
      DATA GAM CS( 41) / +.3376448816 5853380903 3489066666 6 D-30     /
      DATA GAM CS( 42) / -.5793070335 7821357846 2549333333 3 D-31     /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA NGAM, XMIN, XMAX, DXREL / 0, 3*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DGAMMA
      IF (NGAM.NE.0) GO TO 10
      NGAM = INITDS (GAMCS, 42, 0.1*SNGL(D1MACH(3)) )
C
      CALL DGAMLM (XMIN, XMAX)
      DXREL = DSQRT (D1MACH(4))
C
 10   Y = DABS(X)
      IF (Y.GT.10.D0) GO TO 50
C
C COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
C GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
C
      N = X
      IF (X.LT.0.D0) N = N - 1
      Y = X - DBLE(FLOAT(N))
      N = N - 1
      DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
C
      IF (N.GT.0) GO TO 30
C
C COMPUTE GAMMA(X) FOR X .LT. 1.0
C
      N = -N
      IF (X.EQ.0.D0) CALL XERROR ( 'DGAMMA  X IS 0', 14, 4, 2)
      IF (X.LT.0.0 .AND. X+DBLE(FLOAT(N-2)).EQ.0.D0) CALL XERROR ( 'DGAM
     1MA  X IS A NEGATIVE INTEGER', 31, 4, 2)
      IF (X.LT.(-0.5D0) .AND. DABS((X-DINT(X-0.5D0))/X).LT.DXREL) CALL
     1  XERROR ( 'DGAMMA  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NE
     2GATIVE INTEGER', 68, 1, 1)
C
      DO 20 I=1,N
        DGAMMA = DGAMMA/(X+DBLE(FLOAT(I-1)) )
 20   CONTINUE
      RETURN
C
C GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
C
 30   DO 40 I=1,N
        DGAMMA = (Y+DBLE(FLOAT(I))) * DGAMMA
 40   CONTINUE
      RETURN
C
C GAMMA(X) FOR DABS(X) .GT. 10.0.  RECALL Y = DABS(X).
C
 50   IF (X.GT.XMAX) CALL XERROR ( 'DGAMMA  X SO BIG GAMMA OVERFLOWS',
     1  32, 3, 2)
C
      DGAMMA = 0.D0
      IF (X.LT.XMIN) CALL XERROR ( 'DGAMMA  X SO SMALL GAMMA UNDERFLOWS'
     1  , 35, 2, 1)
      IF (X.LT.XMIN) RETURN
C
      DGAMMA = DEXP ((Y-0.5D0)*DLOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
      IF (X.GT.0.D0) RETURN
C
      IF (DABS((X-DINT(X-0.5D0))/X).LT.DXREL) CALL XERROR ( 'DGAMMA  ANS
     1WER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER'  , 61, 1, 1)
C
      SINPIY = DSIN (PI*Y)
      IF (SINPIY.EQ.0.D0) CALL XERROR ( 'DGAMMA  X IS A NEGATIVE INTEGER
     1', 31, 4, 2)
C
      DGAMMA = -PI/(Y*SINPIY*DGAMMA)
C
      RETURN
      END
      FUNCTION INITDS(DOS,NOS,ETA)
C***BEGIN PROLOGUE  INITDS
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C3A2
C***KEYWORDS  CHEBYSHEV,DOUBLE PRECISION,INITIALIZE,
C             ORTHOGONAL POLYNOMIAL,SERIES,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Initializes the d.p. properly normalized orthogonal
C            polynomial series to determine the number of terms needed
C            for specific accuracy.
C***DESCRIPTION
C
C Initialize the double precision orthogonal series DOS so that INITDS
C is the number of terms needed to insure the error is no larger than
C ETA.  Ordinarily ETA will be chosen to be one-tenth machine precision
C
C             Input Arguments --
C DOS    dble prec array of NOS coefficients in an orthogonal series.
C NOS    number of coefficients in DOS.
C ETA    requested accuracy of series.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  INITDS
C
      DOUBLE PRECISION DOS(NOS)
C***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS.LT.1) CALL XERROR ( 'INITDS  NUMBER OF COEFFICIENTS LT 1',
     1 35, 2, 2)
C
      ERR = 0.
      DO 10 II=1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(SNGL(DOS(I)))

cBH        write(*,*)'initds:err,dos(i),ABS(SNGL(DOS(I))):',err,dos(i),
cBH     +            ABS(SNGL(DOS(I)))

        IF (ERR.GT.ETA) GO TO 20
 10   CONTINUE
C
 20   IF (I.EQ.NOS) CALL XERROR ( 'INITDS  ETA MAY BE TOO SMALL', 28,
     1  1, 2)
      INITDS = I
C
      RETURN
      END
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Symbolic dump (should be locally written).
C***DESCRIPTION
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  23 May 1979
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C***BEGIN PROLOGUE  J4SAVE
C***REFER TO  XERROR
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                 = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                 = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                 = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                 = 6 Refers to the 2nd unit for error messages
C                 = 7 Refers to the 3rd unit for error messages
C                 = 8 Refers to the 4th unit for error messages
C                 = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C    Adapted from Bell Laboratories PORT Library Error Handler
C     Latest revision ---  23 MAY 1979
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERABT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Aborts program execution and prints error message.
C***DESCRIPTION
C     Abstract
C        ***Note*** machine dependent routine
C        XERABT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG and NMESSG are as in XERROR, except that NMESSG may
C        be zero, in which case no message is being supplied.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
C***BEGIN PROLOGUE  XERCTL
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Allows user control over handling of individual errors.
C***DESCRIPTION
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCTL.
C        If the user has provided his own version of XERCTL, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        MESSG1 - the first word (only) of the error message.
C        NMESSG - same as in the call to XERROR or XERRWV.
C        NERR   - same as in the call to XERROR or XERRWV.
C        LEVEL  - same as in the call to XERROR or XERRWV.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Prints error messages.
C***DESCRIPTION
C     Abstract
C        Print the Hollerith message in MESSG, of length NMESSG,
C        on each file indicated by XGETUA.
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
C     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
C***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE DGAMLM(XMIN,XMAX)
C***BEGIN PROLOGUE  DGAMLM
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A,R2
C***KEYWORDS  COMPLETE GAMMA FUNCTION,DOUBLE PRECISION,GAMMA FUNCTION,
C             LIMITS,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. minimum and maximum bounds for X in
C            GAMMA(X).
C***DESCRIPTION
C
C Calculate the minimum and maximum legal bounds for X in gamma(X).
C XMIN and XMAX are not the only bounds, but they are the only non-
C trivial ones to calculate.
C
C             Output Arguments --
C XMIN   double precision minimum legal value of X in gamma(X).  Any
C        smaller value of X might result in underflow.
C XMAX   double precision maximum legal value of X in gamma(X).  Any
C        larger value of X might cause overflow.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,XERROR
C***END PROLOGUE  DGAMLM
      DOUBLE PRECISION XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD, D1MACH
C***FIRST EXECUTABLE STATEMENT  DGAMLM
      ALNSML = DLOG(D1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = DLOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML)
     1    / (XMIN*XLN+0.5D0)
        IF (DABS(XMIN-XOLD).LT.0.005D0) GO TO 20
 10   CONTINUE
      CALL XERROR ( 'DGAMLM  UNABLE TO FIND XMIN', 27, 1, 2)
C
 20   XMIN = -XMIN + 0.01D0
C
      ALNBIG = DLOG (D1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = DLOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG)
     1    / (XMAX*XLN-0.5D0)
        IF (DABS(XMAX-XOLD).LT.0.005D0) GO TO 40
 30   CONTINUE
      CALL XERROR ( 'DGAMLM  UNABLE TO FIND XMAX', 27, 2, 2)
C
 40   XMAX = XMAX - 0.01D0
      XMIN = DMAX1 (XMIN, -XMAX+1.D0)
C
      RETURN
      END
c
c
c it is from brambilla code. Smirnov970118 changed it to double precision
ctk--calculation of Bessel function with a complex argument (x,y)
c     br(nb),bi(nb)  
c     y=0,ize=0 >J_0(x)=br(1),J_1(x)=br(2)
c     y=0,ize=1 >I_0(x)=br(1),I_1(x)=br(2),... ,I_n(x)=br(n+1),
c     n=0,1,2,...; nb.ge.1
C     x,y,br,bi - real

      subroutine beslci1(x,y,nb,ize,br,bi,ncalc)
c
c*****     *****     *****     *****     *****     *****     *****
c
c     journal of research  national bureau of standards
c     b-mathematical sciences   vol.77-b  nos.3-4  july-dec. 1973  p118
c     routine to calculate bessel functions j and i
c     of complex argument and integer order
c                  input variables
c     x real part of complex argument
c     y imaginary part of complex argument
c     nb a positive integer designating highest order to be calculated
c     ize     zero for j*s    one for i*s
c     br  for normal exit,br contains real part of solution vector
c     bi  for normal exit,bi contains imaginary part of solution vector
c     normal exit if  ncalc=nb
c     br-bi-ncalc need not be initialized
c                  machine dependent constants
c     nsig  decimal significance desired. set to ifix(alog10(2)*nbit+1)
c     nbit  number of bits in the mantissa
c     the relative truncation error is limited to t=.5*10**-nsig for
c     order greater than abs(z).
c     for order less than abs(z) (general test),the relative error is
c     limited to t for function values of magnitude at least 1.
c     the absolute error is limited to t for smaller values.
c     nten  largest integer k such that 10**k is machine representable.
c     largez  upper limit on the magnitude of z. if abs(z)=n, at least
c     n iterations of the backward recursion will be executed.
c     exparg  largest argument that the library exp routine can handle.
c                  error returns
c     let g denote either i or j.
c     in case of an error, ncalc.ne.nb and not all g*s are calculated
c     to the desired accuracy.
c     if ncalc.lt.0, an argument is out of range. nb.le.0 or ize is
c     neither 0 nor 1 or ize=0 and abs(y).gt.exparg, or ize=1 and
c     abs(x).gt.exparg. in this case,the vectors br and bi are not
c     calculated, and ncalc is set to min0(nb,0)-1 so ncalc.ne.nb.
c     nb.gt.ncalc.gt.0 will occur if nb.gt.magz and abs(g-sub-nb-of-z/g
c     -sub-magx-of-z).lt.10.**(nten/2), i.e. nb is much greater than
c     magz. in this case, br(n) and bi(n) are calculated to the desired
c     accuracy for n.le.ncalc, but for ncalc.lt.n.le.nb, precision  is
c     lost. if n.gt.ncalc and abs(g(ncalc-1)/g(n-1)).eq.10**-k, then
c     the last k significant figures of g(n-1) (=br(n)+i*bi(n)) are
c     erroneous. if the user wishes to calculate g(n-1) to higher
c     accuracy, he should use an asymptotic formula for large order.
c
c*****     *****     *****     *****     *****     *****     *****
      implicit double precision (a-h,o-z)
      dimension br(*),bi(*)
c      data nsig,nten,largez,exparg/15,307,10000,700./
c      data nsig,nten,largez,exparg/15,38,10000,81./
      data nsig,nten,largez,exparg/15,38,10000,308./
      tempar=dsqrt(x*x+y*y)
      magz=int(tempar)
      if(nb.gt.0.and.magz.le.largez.and.
     +((ize.eq.0.and.dabs(y).le.exparg)
     c.or.(ize.eq.1.and.dabs(x).le.exparg))) go to 1
c     error return    z,nb,or ize is out of range
      ncalc=min0(nb,0)-1
      return
    1 sign1=1-2*ize
      ncalc=nb
c     use 2-term ascending series for small z
      if(tempar**4.lt.0.1d0**nsig) go to 50
c     initialize the calculation of the p*s
      nbmz=nb-magz
      n=magz+1
      if(dabs(x).lt.dabs(y)) go to 2
      zinvr=1.0d0/(x+y*y/x)
      zinvi=-y*zinvr/x
      go to 3
    2 zinvi=-1.0d0/(y+x*x/y)
      zinvr=-x*zinvi/y
    3 plastr=1.0d0
      plasti=0.0d0
      pr=sign1*(n+n)*zinvr
      pi=sign1*(n+n)*zinvi
      test=2.0d0*(10.0d0)**nsig
      m=0
      if(nbmz.lt.3) go to 6
c     calculate p*s until n=nb-1.  check for possible overflow.
      tover=10.0d0**(nten-nsig)
      nstart=magz+2
      nend=nb-1
      do 5 n=nstart,nend
      poldr=plastr
      poldi=plasti
      plastr=pr
      plasti=pi
      pr=sign1*((n+n)*(plastr*zinvr-plasti*zinvi)-poldr)
      pi=sign1*((n+n)*(plasti*zinvr+plastr*zinvi)-poldi)
      if((pr/tover)**2+(pi/tover)**2-1.0d0) 5,5,7
    5 continue
      n=nend
c     calculate special significance test for nbmz.gt.2.
      tempbi=dmax1(dabs(pr),dabs(pi))
      tempbi=tempbi*dsqrt(2.0d0*(10.0d0)**nsig*dsqrt(((pr/tempbi)**2+(pi
     c/tempbi)**2)*((plastr/tempbi)**2+(plasti/tempbi)**2)))
      test=dmax1(test,tempbi)
c     calculate p*s until significance test is passed
    6 n=n+1
      poldr=plastr
      poldi=plasti
      plastr=pr
      plasti=pi
      pr=sign1*((n+n)*(plastr*zinvr-plasti*zinvi)-poldr)
      pi=sign1*((n+n)*(plasti*zinvr+plastr*zinvi)-poldi)
      if((pr/test)**2+(pi/test)**2.lt.1.0d0) go to 6
      if(m.eq.1) go to 12
c     calculate strict variant of significance test,and
c     calculate p*s until this test is passed.
      m=1
      tempbi=dmax1(dabs(pr),dabs(pi))
      tempbr=dsqrt(((pr/tempbi)**2+(pi/tempbi)**2)/
     c((plastr/tempbi)**2+(plasti/tempbi)**2))
      tempbi=(n+1)/tempar
      if(tempbr+1.0d0/tempbr.gt.2.0d0*tempbi) then
        tempbr=tempbi+dsqrt(tempbi**2
     c-1.0d0)
      endif

      test=test/dsqrt(tempbr-1.0d0/tempbr)
      if((pr/test)**2+(pi/test)**2-1.0d0) 6,12,12
    7 nstart=n+1
c     to avoid overflow, normalize p*s by dividing by tover.
c     calculate p*s until unnormalized p would overflow.
      pr=pr/tover
      pi=pi/tover
      plastr=plastr/tover
      plasti=plasti/tover
      psaver=pr
      psavei=pi
      tempcr=plastr
      tempci=plasti
      test=10.0d0**(2*nsig)
    8 n=n+1
      poldr=plastr
      poldi=plasti
      plastr=pr
      plasti=pi
      pr=sign1*((n+n)*(plastr*zinvr-plasti*zinvi)-poldr)
      pi=sign1*((n+n)*(plasti*zinvr+plastr*zinvi)-poldi)
      if(pr**2+pi**2.le.test) go to 8
c     calculate backward test,and find ncalc,the highest n
c     such that the test is passed.
      tempbr=dsqrt((plastr**2+plasti**2)/(poldr**2+poldi**2))
      tempbi=n/tempar
      if(tempbr+1.0d0/tempbr.gt.2.0d0*tempbi) then
        tempbr=tempbi+dsqrt(tempbi**2
     c-1.0d0)
      endif

      test=0.5d0*(1.0d0-1.0d0/tempbr**2)/10.0d0**nsig
      test=((plastr**2+plasti**2)*test)*((poldr**2+poldi**2)*test)
      pr=plastr*tover
      pi=plasti*tover
      n=n-1
      nend=min0(nb,n)
      do 9 ncalc=nstart,nend
      poldr=tempcr
      poldi=tempci
      tempcr=psaver
      tempci=psavei
      psaver=sign1*((n+n)*(tempcr*zinvr-tempci*zinvi)-poldr)
      psavei=sign1*((n+n)*(tempci*zinvr+tempcr*zinvi)-poldi)
      if((psaver**2+psavei**2)*(tempcr**2+tempci**2)-test) 9,9,10
    9 continue
      ncalc=nend+1
   10 ncalc=ncalc-1
c     the coefficient of b(n) in the normalized sum is
c     m*sqrt(-1)**imag,where m=-2,0,or2,and imag is 0 or 1.
c     calculate recursion rules for m and imag,and initialize them.
   12 n=n+1
      tempbr=ize*x+(1-ize)*y
      ipos=0
      if(tempbr) 13,14,13
   13 ipos=int(1.1d0*tempbr/dabs(tempbr))
   14 mrecur=4*((2+ize+ipos)/2)-3-2*(ize+ipos)
      k=2+ipos+2*ize*ipos**2-ize
      l=n-4*(n/4)
      mlast=2+8*((k*l)/4)-4*((k*l)/2)
      if(ipos.eq.0.and.(l.eq.1.or.l.eq.3)) mlast=0
      l=l+3-4*((l+3)/4)
      m=2+8*((k*l)/4)-4*((k*l)/2)
      if(ipos.eq.0.and.(l.eq.1.or.l.eq.3)) m=0
      imrecr=(1-ize)*ipos**2
      imag=imrecr*(l-2*(l/2))
c     initialize the backward recursion and the normalization sum
      tempbr=0.0d0
      tempbi=0.0d0
      if(dabs(pi).gt.dabs(pr)) go to 15
      tempar=1.0d0/(pr+pi*(pi/pr))
      tempai=-(pi*tempar)/pr
      go to 16
   15 tempai=-1.0d0/(pi+pr*(pr/pi))
      tempar=-(pr*tempai)/pi
   16 if(imag.ne.0) go to 17
      sumr=m*tempar
      sumi=m*tempai
      go to 18
   17 sumr=-m*tempai
      sumi=m*tempar
   18 nend=n-nb
      if(nend) 26,22,19
c     recur backward via difference equation calculating (but not
c     storing) br(n) and bi(n) until n=nb
   19 do 21 l=1,nend
      n=n-1
      tempcr=tempbr
      tempci=tempbi
      tempbr=tempar
      tempbi=tempai
      pr=(n+n)*zinvr
      pi=(n+n)*zinvi
      tempar=pr*tempbr-pi*tempbi-sign1*tempcr
      tempai=pr*tempbi+pi*tempbr-sign1*tempci
      imag=(1-imag)*imrecr
      k=mlast
      mlast=m
      m=k*mrecur
      if(imag.ne.0) go to 20
      sumr=sumr+m*tempar
      sumi=sumi+m*tempai
      go to 21
   20 sumr=sumr-m*tempai
      sumi=sumi+m*tempar
   21 continue
c     store  br(nb),bi(nb)
   22 br(n)=tempar
      bi(n)=tempai
      if(n.gt.1) go to 23
c     nb=1.  since 2*tempar and 2*tempai were added to sumr and sumi
c     respectively,we must subaxis(i)ract tempar and tempai
      sumr=sumr-tempar
      sumi=sumi-tempai
      go to 35
c     calculate and store br(nb-1),bi(nb-1)
   23 n=n-1
      pr=(n+n)*zinvr
      pi=(n+n)*zinvi
      br(n)=pr*tempar-pi*tempai-sign1*tempbr
      bi(n)=pr*tempai+pi*tempar-sign1*tempbi
      if(n.eq.1) go to 34
      imag=(1-imag)*imrecr
      k=mlast
      mlast=m
      m=k*mrecur
      if(imag.ne.0) go to 24
      sumr=sumr+m*br(n)
      sumi=sumi+m*bi(n)
      go to 30
   24 sumr=sumr-m*bi(n)
      sumi=sumi+m*br(n)
      go to 30
c     n.lt.nb,so store br(n),bi(n) and set higher orders zero
   26 br(n)=tempar
      bi(n)=tempai
      nend=-nend
      do 27 l=1,nend
      br(n+l)=0.0d0
   27 bi(n+l)=0.0d0
   30 nend=n-2
      if(nend.eq.0) go to 33
c     calculate via difference equation and store br(n),bi(n)
c     until n=2
      do 32 l=1,nend
      n=n-1
      pr=(n+n)*zinvr
      pi=(n+n)*zinvi
      br(n)=pr*br(n+1)-pi*bi(n+1)-sign1*br(n+2)
      bi(n)=pr*bi(n+1)+pi*br(n+1)-sign1*bi(n+2)
      imag=(1-imag)*imrecr
      k=mlast
      mlast=m
      m=k*mrecur
      if(imag.ne.0) go to 31
      sumr=sumr+m*br(n)
      sumi=sumi+m*bi(n)
      go to 32
   31 sumr=sumr-m*bi(n)
      sumi=sumi+m*br(n)
   32 continue
c     calculate and store br(1),bi(1)
   33 br(1)=2.0d0*(br(2)*zinvr-bi(2)*zinvi)-sign1*br(3)
      bi(1)=2.0d0*(br(2)*zinvi+bi(2)*zinvr)-sign1*bi(3)
   34 sumr=sumr+br(1)
      sumi=sumi+bi(1)
c     calculate normalization factor. tempar+i*tempai
   35 if(ize.eq.1) go to 36
      tempcr=ipos*y
      tempci=-ipos*x
      go to 37
   36 tempcr=ipos*x
      tempci=ipos*y
   37 tempcr=dexp(tempcr)
      tempbr=dcos(tempci)
      tempbi=dsin(tempci)
      if(dabs(sumr).lt.dabs(sumi)) go to 38
      tempci=sumi/sumr
      tempcr=(tempcr/sumr)/(1.0d0+tempci*tempci)
      tempar=tempcr*(tempbr+tempbi*tempci)
      tempai=tempcr*(tempbi-tempbr*tempci)
      go to 39
   38 tempci=sumr/sumi
      tempcr=(tempcr/sumi)/(1.0d0+tempci*tempci)
      tempar=tempcr*(tempbr*tempci+tempbi)
      tempai=tempcr*(tempbi*tempci-tempbr)
c     normalize
   39 do 40 n=1,nb
      tempbr=br(n)*tempar-bi(n)*tempai
      bi(n)=br(n)*tempai+bi(n)*tempar
   40 br(n)=tempbr
      return
c     two-term ascending series for small z
   50 tempar=1.0d0
      tempai=0.0d0
      tempcr=0.25d0*(x*x-y*y)
      tempci=0.5d0*x*y
      br(1)=1.0d0-sign1*tempcr
      bi(1)=-sign1*tempci
      if(nb.eq.1) go to 52
      do 51 n=2,nb
      tempbr=(tempar*x-tempai*y)/(n+n-2)
      tempai=(tempar*y+tempai*x)/(n+n-2)
      tempar=tempbr
      tempbr=n
      br(n)=tempar*(1.0d0-sign1*tempcr/tempbr)+tempai*tempci/tempbr
   51 bi(n)=tempai*(1.0d0-sign1*tempcr/tempbr)-tempar*tempci/tempbr
   52 return
      end
