C     ******************besk************************
c     *      subroutine for Bessel function K_N(x) *
c     *      input data:                           *
c     *      x - argument                          *
c     *      N -index                              *
c     *      output data:                          *
c     *      bk - Bessel function                  *
c     *      This code was in Mazzucato code       *
c***************************************************
      SUBROUTINE BESK(X,N,BK,IER)
c	:::::::larisa::::::::::::::
	implicit double precision(a-h,o-z)
c	::::::::::::::::::::::::::
      DIMENSION T(12)
c	""""""""""""""""""""""""""
c	write(*,*)'beg. of besk'
c	""""""""""""""""""""""""""

      BK=.0d0
      IF(X.GT.150.d0)RETURN
      IF(N)10,11,11
   10 IER=1
      RETURN
   11 IF(X)12,12,20
   12 IER=2
      RETURN
   20 IF(X-170.0d0)22,22,21
   21 IER=3
      RETURN
   22 IER=0
      IF(X-1.d0)36,36,25
   25 A=dEXP(-X)
      B=1.d0/X
      C=dSQRT(B)
      T(1)=B
      DO 26 L=2,12
   26 T(L)=T(L-1)*B
      IF(N-1)27,29,27
   27 G0=A*(1.25331414d0-.15666418d0*T(1)+.088111278d0*T(2)
     *-.091390954d0*T(3)
     2+.13445962d0*T(4)-.22998503d0*T(5)+.37924097d0*T(6)
     *-.52472773d0*T(7)
     3+.55753684d0*T(8)-.42626329d0*T(9)+.21845181d0*T(10)
     *-.066809767d0*T(11)
     4+.009189383d0*T(12))*C
      IF(N)20,28,29
   28 BK=G0
      RETURN
   29 G1=A*(1.2533141d0+.46999270d0*T(1)-.14685830d0*T(2)
     *+.12804266d0*T(3)
     2-.17364316d0*T(4)+.28476181d0*T(5)-.45943421d0*T(6)+
     *.62833807d0*T(7)
     3-.66322954d0*T(8)+.50502386d0*T(9)-.25813038d0*T(10)
     *+.078800012d0*T(11)
     4-.010824177d0*T(12))*C
      IF(N-1)20,30,31
   30 BK=G1
      RETURN
   31 DO 35 J=2,N
      GJ=2.d0*(dFLOAT(J)-1.d0)*G1/X+G0
      IF(GJ-1.0d35)33,33,32
   32 IER=4
      GO TO 34
   33 G0=G1
   35 G1=GJ
   34 BK=GJ
      RETURN
   36 B=X/2.d0
      A=.57721566d0+DLOG(B)
      C=B*B
      IF(N-1)37,43,37
   37 G0=-A
      X2J=1.d0
      FACT=1.d0
      HJ=.0d0
      DO 40 J=1,6
      RJ=1.d0/dFLOAT(J)
      X2J=X2J*C
      FACT=FACT*RJ*RJ
      HJ=HJ+RJ
   40 G0=G0+X2J*FACT*(HJ-A)
      IF(N)43,42,43
   42 BK=G0
      RETURN
   43 X2J=B
      FACT=1.d0
      HJ=1.d0
      G1=1.d0/X+X2J*(.5d0+A-HJ)
      DO 50 J=2,8
      X2J=X2J*C
      RJ=1.d0/dFLOAT(J)
      FACT=FACT*RJ*RJ
      HJ=HJ+RJ
   50 G1=G1+X2J*FACT*(.5d0+(A-HJ)*dFLOAT(J))
      IF(N-1)31,52,31
   52 BK=G1

c	""""""""""""""""""""""""""
c	write(*,*)'end of besk'
c	""""""""""""""""""""""""""


      RETURN
      END



