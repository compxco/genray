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

      KS1 = 2 - KS
      XKS1 = KS1
      KS2 = KS - 1
      XKS2 = KS2

      IS = MOD(IT-1,18)
      IF (IS.LT.0) IS = IS + 18
      IS = IS/2 + 1

      ALM = 2.d0 *APARAM/((N1M1*N1M1+N2M1*N2M1))

      IF (APARAM.LE.0.0d0 ) GO TO 260

      IF (ALM.GT.1.0d0 ) GO TO 280

      ALPHA = 1.d0  - ALM**ALP(IS)

      DO 160 JL = 1, N2

         J = KS1*JL + KS2*(N2P1-JL)
         JB = J - KS1 + KS2

         DO 140 I = 1, N1
            IB = I - 1
            CS = (1.0d0 +mult*(bx(I)+by(J))) 
            IF ((JB.EQ.0) .OR. (JB.EQ.N2P1)) GO TO 20 
            SEIJB = WRKSP1(I,JB)
            SFIJB = WRKSP2(I,JB)
            RIJB = R(I,JB)
            GO TO 40

   20       CONTINUE

            SEIJB = 0.0d0 
            SFIJB = 0.0d0 
            RIJB = 0.0d0 
   40       CONTINUE

            IF (I.EQ.1) GO TO 60

            SEIBJL = WRKSP1(IB,J)
            SFIBJL = WRKSP2(IB,J)
            RIBJL = R(IB,J)
            GO TO 80

   60       CONTINUE

            SEIBJL = 0.0d0 
            SFIBJL = 0.0d0 
            RIBJL = 0.0d0 
   80       CONTINUE

            SB = (XKS1*(mult*ay(J)) +XKS2*(mult*cy(J)) )/
     *              (1.d0 +ALPHA*SEIJB)
            SC = (mult*ax(I)) /(1.d0 +ALPHA*SFIBJL)

            SBSEIJ = SB*SEIJB
            SCSFIJ = SC*SFIBJL

            SD = 1.d0 /(-SB*SFIJB-SC*SEIBJL+CS+ALPHA*(SBSEIJ+SCSFIJ))

            WRKSP1(I,J) = ((mult*cx(I)) -ALPHA*SBSEIJ)*SD
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
