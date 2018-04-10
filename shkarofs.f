cSubroutines from the Shkarofsky (mildly relativistic dispersion) 
	SUBROUTINE DF(X,Y,LR,G,U,WLD,WK,RP,BMAG,DY,TEMP,BM,GRB,GRD,
     1	 GRT,GRBIJ,E,DR,DI,DP,DH,DK,DW,DX,DKI)
C
C Calls up RQ and PRE for each N = 0, +-1, +-2 and adds up the N contribution
C to give E(6,9),the six dielectric tensor elements and their eight derivative
C To obtain the last three derivatives, one uses a model for the magnetic field
C vector,  magnitude and derivative, density and temperature and derivatives. 
C Calls EDET to provide the Real DR and Imaginary DI parts of the dispersion 
C determinant D  and its derivatives.
	COMPLEX V(7),DV(7),D2V(8),D3V(8),D4V(9),R(14),E(6,9),EI,
     1	 ZERO,DP,DH,DK,DW,DX(3),DKI(3),CMPLX,HA(6)
	DIMENSION DNBDX(3),DNPDX(3),WK(3),BM(3),GRB(3),
     1	 GRD(3),GRT(3),GRBIJ(3,3),DNBDK(3),DNPDK(3)
	COMMON/VVV/V,DV,D2V,D3V,D4V
	WKB=G/WLD
	SU=SQRT(U)
	Y2=Y*Y
	YD=2.*Y
	B=SU*G
	CL=RP*RP/(Y2*U)
	WLDS=WLD*WLD
	EI=CMPLX(0.0,1.0)
	ZERO=CMPLX(0.0,0.0)
	DO 5 I=1,3
	DNPDX(I)=0.0
	DNBDX(I)=0.0
    5	DKI(I)=ZERO
	DO 6 I=1,6
	DO 6 J=1,9
    6	E(I,J)=ZERO
	MF=6
	IF(LR.EQ.0) MF=3
	DO 11 L=1,9
	N=L-5
	YN=Y*N
	DELTA=1.-YN
	A=U*DELTA
	CALL FQ(A,B)
	DO 10 M=1,MF
	ML=M
    	CALL RQ(Y,ML,G,U,N,R)
	CALL PRE(X,CL,ML,R,G,U,N,HA)
	DO 9 I=1,6
	SG=1.
	IF((N.LT.0).AND.((I.EQ.2) .OR. (I.EQ.3))) SG=-1.
    9	E(I,ML)=E(I,ML)+SG*HA(I)
   10 	CONTINUE
   11	CONTINUE
	E(1,1)=1.+E(1,1)
	E(4,1)=1.+E(4,1)
	E(6,1)=1.+E(6,1)
        DO 24 I=1,6
   24   E(I,2)=E(I,2)/(Y*SU)
	DO 15 ML=1,6
	IF((ML.GT.3).AND.(LR.EQ.0)) GO TO 30
	E(2,ML)=EI*E(2,ML)
	E(5,ML)=EI*E(5,ML)
   15	CONTINUE
	DO 22 II=1,3
	I=II
	DNBDK(I)=WLD*BM(I)/BMAG
	DNPDK(I)=WLDS*(WK(I)-WKB*BM(I)/BMAG)/RP
	DNBDR=0.
	DO 16 JJ=1,3
	J=JJ
   16 	DNBDR=DNBDR+WLD*WK(J)*(BMAG*GRBIJ(J,I)-BM(J)*GRB(I))/(BMAG**2)
	DNPDX(I)=-DNBDR*G/RP
	DNBDX(I)=DNBDR
   22	CONTINUE
	E(1,1)=E(1,1)-1.
	E(4,1)=E(4,1)-1.
	E(6,1)=E(6,1)-1.
	DO 23 I=1,6
	DO 23 J=1,3
	E(I,6+J)=E(I,2)*DNPDX(J)+E(I,3)*DNBDX(J)+E(I,1)*GRD(J)/DY
     1	 +E(I,5)*GRT(J)/TEMP+E(I,6)*GRB(J)/BMAG
   23	CONTINUE
	E(1,1)=E(1,1)+1.
	E(4,1)=E(4,1)+1.
	E(6,1)=E(6,1)+1.
   30	CALL EDET(E,RP,G,LR,DNBDX,DNPDX,DR,DI,DP,DH,DK,DW,DX)
	IF(LR.EQ.0) RETURN
	DO 35 II=1,3
	DKI(II)=DH*DNBDK(II)+DP*DNPDK(II)
   35	CONTINUE
	RETURN
	END

	SUBROUTINE EDET(E,U1,U3,LR,DNBDX,DNPDX,DR,DI,DP,DH,DK,DW,DX)
C
c Obtains from DF the values of E (dielectric tensor elements and derivatives),
C U1=Nperp, U3=Npar, DNBDX = dNpar/dr and DNPDX = dNperp/dr; then it calculate
C the determinant and its derivatives as mentioned above. Only kperp, kpar and
C k derivatives are provided when LR=0
C
	COMPLEX E(6,9),A,B,C,D,DA,DB,DC(6),DP,DH,DK,DW,UE,E2,DX(3),
     1    CMPLX,ZERO
	COMMON/ALL/IBUG,IQ,IIO,CFE0,PFE0,ST0,SSGN,PI,WLD,IIQ,ND125
	DIMENSION DNBDX(3),DNPDX(3)
c	COMMON/HHH/
	U12=U1*U1
	U32=U3*U3
	U2=U12+U32
	U13=2*U1*U3
	UE=U3*E(5,1)-U1*E(2,1)
	E2=E(1,1)*E(6,1)-E(3,1)**2
	ZERO=CMPLX(0.,0.)
	A=U12*E(1,1)+U13*E(3,1)+U32*E(6,1)
	B=U2*E2+UE**2
	C=E2*E(4,1)+E(6,1)*(E(2,1)**2)+(E(1,1)*E(5,1)
     +	  +2.*E(2,1)*E(3,1))*E(5,1)
	D=(U2-E(4,1))*A-B+C
	D=D/E(1,1)
	DR=REAL(D)
	DI=AIMAG(D)
	MA=8
	IF(LR.EQ.0) MA=2
	DO 10 I=1,MA
	IF((I.EQ.4).OR.(I.EQ.5)) GO TO 10
	J=I+1
	L=I
	IF(I.GT.3) L=I-2
	DC(L)=(E(1,J)*E(6,1)+E(6,J)*E(1,1)-2.*E(3,J)*E(3,1))*E(4,1)
     +	 +E(5,1)*(E(1,J)*E(5,1)+2.*E(5,J)*E(1,1)+E(3,J)*E(2,1)
     +	 +2.*E(2,J)*E(3,1)) +E2*E(4,J)+E(2,1)*(E(6,J)*E(2,1)
     +	 +2.*E(2,J)*E(6,1)+E(3,J)*E(5,1)+2.*E(5,J)*E(3,1))
   10	CONTINUE
	DA=U12*E(1,2)+2.*U1*E(1,1)+U13*E(3,2)+2.*U3*E(3,1)+U32*E(6,2)
	DB=2.*UE*(U3*E(5,2)-U1*E(2,2)-E(2,1))+2.*U1*E2
     +	  +U2*(E(1,2)*E(6,1)+E(1,1)*E(6,2)-2.*E(3,1)*E(3,2))
	DP=(U2-E(4,1))*DA+(2.*U1-E(4,2))*A-DB+DC(1)
	DP=(DP-D*E(1,2))/E(1,1)
	DA=U12*E(1,3)+U13*E(3,3)+2.*U1*E(3,1)+U32*E(6,3)+2.*U3*E(6,1)
	DB=2.*UE*(U3*E(5,3)+E(5,1)-U1*E(2,3))+2.*U3*E2
     +	  +U2*(E(1,3)*E(6,1)+E(1,1)*E(6,3)-2.*E(3,1)*E(3,3))
	DH= (U2-E(4,1))*DA+(2.*U3-E(4,3))*A-DB+DC(2)
	DH=(DH-D*E(1,3))/E(1,1)
	DK=DH*U3+DP*U1
	IF(LR.EQ.0) GO TO 30
	IF((IQ.LT.1).OR.(ND125.EQ.0)) NNN=1
	DA=(E(1,4)-2.*E(1,1))*U12+U13*(E(3,4)-2.*E(3,1))
     +	  +(E(6,4)-2.*E(6,1))*U32
	DB=2.*UE*(U3*(E(5,4)-E(5,1))-U1*(E(2,4)-E(2,1)))
     +	  +U2*((E(1,4)-E(1,1))*E(6,1)-2.*(E(3,4)-E(3,1))*E(3,1)
     +	  +(E(6,4)-E(6,1))*E(1,1))
	DW= (U2-E(4,1))*DA-(2.*U2+E(4,4))*A-DB+DC(3)
	DW=(DW-D*E(1,4))/E(1,1)
  	DO 20 I=1,3
	J=I+6
	DA=U1*2.*DNPDX(I)*E(1,1)+U12*E(1,J)+U13*E(3,J)+U32*E(6,J)
     +	+2.*(U3*DNPDX(I)+U1*DNBDX(I))*E(3,1)+U3*2.*DNBDX(I)*E(6,1)
	DB=2.*UE*(DNBDX(I)*E(5,1)+U3*E(5,J)
     +	  -DNPDX(I)*E(2,1)-U1*E(2,J))+U2*(E(1,J)*E(6,1)
     +	  +E(1,1)*E(6,J)-2.*E(3,1)*E(3,J))
	DX(I)=(U2-E(4,1))*DA-E(4,J)*A-DB+DC(I+3)
	DX(I)=(DX(I)-D*E(1,J))/E(1,1)
   20	CONTINUE

	IF(IBUG.GT.2)THEN
	  IF(NNN.EQ.1) WRITE(11,25) DW,DX(1),DX(2),DX(3) 
   25	     FORMAT(' COMPLEX PARTS',
     +	       T20,'DW  =',2(E10.3,2X),T50,'DX1 =',2(E10.3,2X),/,
     +	       T20,'DX2 =',2(E10.3,2X),T50,'DX3 =',2(E10.3,2X))
        ENDIF
	RETURN
   30	DW=ZERO

	DO 40 I=1,3
   40	DX(I)=ZERO
	RETURN
	END
C

	SUBROUTINE RQ(Y,ML,G,U,N,R)
C
C calculates the R terms.
C
	COMPLEX R(14),V(7),DV(7),D2V(8),D3V(8),D4V(9),ZERO,CMPLX
	COMMON/VVV/V,DV,D2V,D3V,D4V
	YN=Y*N
	DELTA=1.-YN
	G2=G*G
	GG=G2*U
	ZERO=CMPLX(0.0,0.0)
	GO TO (40,50,60,70,80,90) ML
   40	DO 41 I=1,5
   41	R(I)=V(I+2)
	DO 42 I=6,9
   42	R(I)=D2V(I-1)
	DO 43 I=10,13
   43	R(I)=DV(I-6)
	R(14)=D2V(4)
	RETURN
   50	R(1)=ZERO
	DO 51 I=2,5
   51	R(I)=2.*(I-1)*V(I+2)
	DO 52 I=6,9
   52	R(I)=(2*I-10)*D2V(I-1)
	DO 53 I=10,13
   53	R(I)=(2*I-19)*DV(I-6)
	R(14)=ZERO
	RETURN
   60	DO 61 I=1,5
   61	R(I)=G*U*D2V(I+3)
C        G*U IS USED HERE FOR ML=3 INSTEAD OF GG, BECAUSE WE EVALUATE
C        d(EPSILON)/d(G), NOT Gd(EPSILON)/d(G). FOR THE OTHER R(I) IN
C        #62-63, A SIMILAR G FACTOR IS DIVIDED OUT IN SUBROUTINE PRE.
	DO 62 I=6,9
   62	R(I)=GG*D4V(I)+2.*D2V(I-1)
	DO 63 I=10,13
   63	R(I)=GG*D3V(I-5)+DV(I-6)
	R(14)=GG*D4V(5)+2.*D2V(4)
	RETURN
   70	DO 71 I=1,5
        IF(N.EQ.0) R(I)=-GG*D2V(I+3)-2.*V(I+2)
   71	IF(N.NE.0) R(I)=-U*(YN*DV(I+2)+G2*D2V(I+3))-2.*V(I+2)
	DO 72 I=6,9
        IF(N.EQ.0) R(I)=-GG*D4V(I)-4.*D2V(I-1)
   72	IF(N.NE.0) R(I)=-U*(YN*D3V(I-1)+G2*D4V(I))-4.*D2V(I-1)
	DO 73 I=10,13
        IF(N.EQ.0) R(I)=-GG*D3V(I-5)-3.*DV(I-6)
   73	IF(N.NE.0) R(I)=-U*(YN*D2V(I-6)+G2*D3V(I-5))-3.*DV(I-6)
	R(14)=-GG*D4V(5)-4.*D2V(4)
	RETURN
   80	R(1)=GG*0.5*D2V(5)-2.5*DV(4)
        R(2)=U*(DELTA*DV(4)-0.5*G2*D2V(5))
        R(6)=U*(DELTA*D3V(5)-0.5*G2*D4V(6))-D2V(5)
        R(7)=U*(DELTA*D3V(6)-0.5*G2*D4V(7))
        R(10)=U*(DELTA*D2V(4)-0.5*G2*D3V(5))-DV(4)
        R(11)=U*(DELTA*D2V(5)-0.5*G2*D3V(6))
 	DO 81 I=3,5
   81	R(I)=U*(DELTA*DV(I+2)-0.5*G2*D2V(I+3))+(I-2)*V(I+2)
	DO 82 I=8,9
   82	R(I)=U*(DELTA*D3V(I-1)-0.5*G2*D4V(I))+(I-7)*D2V(I-1)
	DO 83 I=12,13
   83   R(I)=U*(DELTA*D2V(I-6)-0.5*G2*D3V(I-5))+(I-11)*DV(I-6)
	R(14)=U*(D3V(4)-0.5*G2*D4V(5))-2.*D2V(4)
        RETURN
   90   IF(N.EQ.0) GO TO 100
        R(1)=U*YN*DV(3)
        DO 91 I=2,5
   91   R(I)=U*YN*DV(I+2)-2.*(I-1)*V(I+2)
        DO 92 I=6,9
   92   R(I)=U*YN*D3V(I-1)-2.*(I-5)*D2V(I-1)
        DO 93 I=10,13
   93   R(I)=U*YN*D2V(I-6)-(2*I-19)*DV(I-6)
	R(14)=ZERO
	RETURN
  100   R(1)=ZERO
        DO 101 I=2,5
  101   R(I)=-2.*(I-1)*V(I+2)
        DO 102 I=6,9
  102   R(I)=-2.*(I-5)*D2V(I-1)
        DO 103 I=10,13
  103   R(I)=-(2*I-19)*DV(I-6)
        R(14)=ZERO
        RETURN
	END
C

	SUBROUTINE PRE(X,CL,ML,R,G,U,N,HA)
C
C It combines the various R terms in the dielectric tensor elements or their 
C derivatives for a given N in W - N*Wc. 
C
c	COMPLEX R(14),ZERO,HA(6),H4R1,H4R2
	COMPLEX R(14),ZERO,HA(6),H4R1,H4R2,CMPLX
	CT=0.625
	X5=X/2.0
	H4=X5*U
	SCL=SQRT(CL)
	SU=SQRT(U)
	CL2=CL*CL
	CL3=CL2*CL
	G2=G*G
        IF(ML.EQ.3) G2=G
C       FOR ML=3, WE EVALUATE d(..)/dG, NOT Gd(..)/dG. 
C       A SIMILAR G FACTOR WAS DIVIDED OUT IN SUBROUTINE
C       RQ FOR R(I),I=1,5.
	ZERO=CMPLX(0.,0.)
	L=N+5
	CTT=7./24.
        CTT2=CTT*2.
        H4R1=ZERO
        H4R2=ZERO
        H2=H4*SCL
        H5=H4*SU
        IF(ML.EQ.2) GO TO 5
C       FOR ML=2, WE EVALUATE d(..)/d(SCL), NOT (SCL)d(..)/d(SCL). 
        H4R1=H4*R(1)
        H4R2=H4*(R(1)+U*G2*R(14))
        H2=H4*CL
        H5=H4*SU*SCL
    5   CONTINUE
        IF(ML.NE.3) H5=H5*G
	H6=H5*CL/2.
        HC3=H2*CL*3./8.
        HC4=H2*CL2/12.
	GO TO (10,20,30,40,50,40,30,20,10) L
   40	CONTINUE
        HA(1)=-H4R1+H2*(   R(2)   -CT*CL*R(3)+CL2*R(4)*CTT)
	HA(4)=-H4R1+H2*(3.*R(2)-4.625*CL*R(3)+CL2*R(4)*97./24.)
	HA(2)= H4R1-H2*(2.*R(2)-3.*CT*CL*R(3)+CL2*R(4)*7./6.)
	HA(6)=     -H2*(   R(2)-CL*R(3)+CL2*R(4)*CT-CL3*R(5)*CTT
     +              +U*G2*(R(6)-CL*R(7)+CL2*R(8)*CT-CL3*R(9)*CTT))
	HA(3)=-H5*(R(10)   -CL*R(11)   +CT*CL2*R(12)-CL3*R(13)*CTT)
	HA(5)=-H5*(R(10)-2.*CL*R(11)+3.*CT*CL2*R(12)-CL3*R(13)*7./6.)
  	RETURN
   30	HA(1)=-H2*(R(2)    -CL*R(3)+CL2*R(4)*CTT2)
	HA(4)=-H2*(R(2) -2.*CL*R(3)+CL2*R(4)*25./12.)
	HA(2)= H2*(R(2)-1.5*CL*R(3)+CL2*R(4)*7./6.)
	HA(6)=    -0.25*H2*CL*(R(3) -CL*R(4)+CL2*R(5)*CTT2
     +                  +U*G2*(R(7) -CL*R(8)+CL2*R(9)*CTT2))
	HA(3)=-H6*(R(11)    -CL*R(12)+CL2*R(13)*CTT2)
	HA(5)=-H6*(R(11)-1.5*CL*R(12)+CL2*R(13)*7./6.)
	RETURN
   50	DO 60 I=1,3
   60	HA(I)=ZERO
	HA(4)=-4.*H2*(R(2)-1.5*CL*R(3)+1.25*CL2*R(4))
	HA(5)=H5*(2.*R(10)-3.*CL*R(11)+2.5*CL2*R(12)-CL3*R(13)*35./24.)
	HA(6)=-2.*H4R2
     1    +2.*H2*(R(2)-CL*R(3)*0.75+CL2*R(4)*5./12.-CL3*R(5)*35./192.  
     2     +U*G2*(R(6)-CL*R(7)*0.75+CL2*R(8)*5./12.-CL3*R(9)*35./192.))
	RETURN
   20   HA(1)=-HC3*(R(3)-CL*R(4))
        HA(4)=-HC3*(R(3)-CL*R(4)*5./3.)
        HA(2)= HC3*(R(3)-CL*R(4)*4./3.)
        HA(6)=-(HC4/2.)*(R(4)-CL*R(5)+U*G2*(R(8)-CL*R(9)))
        HA(3)=-H5*CL2*(R(12)   -CL*R(13))/8.
        HA(5)=-H5*CL2*(R(12)/8.-CL*R(13)/6.)
        RETURN
   10   HA(1)=-HC4*R(4)
        HA(4)= HA(1)
        HA(2)=-HA(1)
        HA(6)=-(HC4*CL/16.)*(R(5)+U*G2*R(9))
        HA(3)=-(H5*CL3/48.)*R(13)
        HA(5)=HA(3)
       RETURN
       END
        FUNCTION FQM(Q,M,A,B)
C
C Provides the Real contribution to Fq in FQ1 as a function of A=mu*delta 
C and B=kpar*c**2/W*Vth
C M order of derivative and subscript Q=q.
C
c        IMPLICIT DOUBLE PRECISION(A-E,G-H,O-Z)
c	REAL*8 FACM
        EPS=1.E-5
        TEO=1.
        P2=B*B/2.
        IF(P2.EQ.A) GO TO 10
        TNO=1.
        TN=1.
        SUMT=1.
        SUMA=0.
        R=P2/(A*A)
        DO 4 N=1,50
        XN=N
        SUM=TN
        DO 2 NS=1,50
        XS=NS
        TEN=-TEO*(XS-2.+2.*XN+M)*(Q+XS+XN-2.)/(XS*A)
        IF(ABS(TEN).GT.ABS(TEO)) GO TO 3
        IF(ABS(TEN/SUM).LT.EPS) GO TO 3
        SUM=SUM+TEN
        TEO=TEN
    2   CONTINUE
    3   TN=TNO*R*(2.*XN+M-1.)*(2.*XN+M)/XN
        SUMT=SUMT+TN
        SUMA=SUMA+SUM
        IF(ABS(TN).GT.ABS(TNO)) GO TO 5
        IF(ABS(TN/SUMT).LT.EPS) GO TO 5
        TNO=TN
        TEO=TN
    4   CONTINUE
    5   CONTINUE
	FACM=FAC(M)
        FQM=SUMA*FACM/(A**(M+1))
        RETURN
 10     CONTINUE
        SUM=1.
        DO 20 NS=1,50
        XS=NS
        TEN=TEO*(XS+M)*(-Q+XS+M+1.)/(XS*A)
        IF(ABS(TEN).GT.ABS(TEO)) GO TO 30
        IF(ABS(TEN/SUM).LT.EPS) GO TO 30
        SUM=SUM+TEN
        TEO=TEN
   20   CONTINUE
   30   CONTINUE
	FACM=FAC(M)
        FQM=SUM*FACM/(A**(M+1))
        RETURN
        END
C
        SUBROUTINE FQ1(A,B)
C
Calculates Fq and derivatives when mu*delta > 100 and |kpar*Vth/(W-n*Wc0|<0.038
C using a double expansion in these arguments. 
C The Imaginary part is also deduced.
C
        COMPLEX V(7),DV(7),D2V(8),D3V(8),D4V(9),CMPLX,
     1	 PHI,XX(2),ZE,ZSQ,CC(9,5),CST(10),CEXPP
        COMMON/VVV/V,DV,D2V,D3V,D4V
        COMMON/CCCC/CFACT
C SMIRNOV 970103 BEG
	CFACT=1.E35
C SMIRNOV 970103 END
        AFACT=ALOG(CFACT)
        BFACT=1.
        PSI=B/SQRT(2.)
        P2=B*B/2.
        PH=P2-A
        ZERO=0.0
        DO 2 K=2,7
        Q=K-0.5
        M=0
        RV=FQM(Q,M,A,B)
    2   V(K)=CMPLX(RV,ZERO)
        DO 3 K=3,7
        Q=K-0.5
        M=1
        RV=FQM(Q,M,A,B)
    3   DV(K)=CMPLX(RV,ZERO)
        DO 4 K=4,8
        Q=K-0.5
        M=2
        RV=FQM(Q,M,A,B)
    4   D2V(K)=CMPLX(RV,ZERO)
        DO 5 K=4,8
        Q=K-0.5
        M=3
        RV=FQM(Q,M,A,B)
    5   D3V(K)=CMPLX(RV,ZERO)
        DO 6 K=5,9
        Q=K-0.5
        M=4
        RV=FQM(Q,M,A,B)
    6   D4V(K)=CMPLX(RV,ZERO)
        IF(PH.LT.ZERO) GO TO 7
        IF(PH.EQ.ZERO) RETURN
        PHI=CMPLX(SQRT(PH),ZERO)
        GO TO 9
    7   RI=-SQRT(-PH)
        PHI=CMPLX(ZERO,RI)
    9   XX(1)=PSI-PHI
        XX(2)=-PSI-PHI	
        IF((AIMAG(XX(1)).GT.ZERO).AND.(AIMAG(XX(2)).GT.ZERO)) RETURN
        P4=P2*P2
        DO 46 K=1,2
        ZE=XX(K)
        ZSQ=ZE*ZE
        IF(K.EQ.1) KK=1
        IF(K.EQ.2) KK=-1
        PSI=KK*PSI
        P3=P2*PSI
        P5=P2*P3
        X=REAL(ZE)
        Y=AIMAG(ZE)
        XA=ABS(X)
        YA=ABS(Y)
        IF(Y.GT.ZERO) GO TO 46
        CS=COS(2.*X*Y)
        SS=SIN(2.*X*Y)
        SR=SIGN(BFACT,CS)
        SI=SIGN(BFACT,SS)
        IF(CS.EQ.ZERO) SR=ZERO
        IF(SS.EQ.ZERO) SI=ZERO
        CST(1)=CMPLX(ZERO,3.544907702)
        IF(Y.EQ.ZERO) CST(1)=0.5*CST(1)
        CST(2)=-2.*ZE*CST(1)
        DO 10  L=3,5
   10   CST(L)=-2.*(L-2)*CST(L-2)-2.*ZE*CST(L-1)
        XXX=EXP(-AFACT/8.0)
        IF(P2.LT.XXX) GO TO 140
        CC(1,1)=-CST(1)/(2.*PHI)
        CC(2,1)=-CST(1)/(2.*PSI)
        DO 11 L=3,7
   11   CC(L,1)=(PH*CC(L-2,1)-(L-2.5)*CC(L-1,1))/P2
        CC(2,2)=CST(2)/(4.*PSI*PHI)
        CC(3,2)=(CST(2)-CST(1)/PSI)/(4.*P2)
        DO 12 L=4,7
   12   CC(L,2)=(-(L-2.5)*CC(L-1,2)+CC(L-2,1)+PH*CC(L-2,2))/P2
        CC(3,3)=(-CST(3)+CST(2)/PSI)/(8.*PHI*P2)
        CC(4,3)=(-CST(3)+3.*CST(2)/PSI-3.*CST(1)/P2)/(8.*P3)
        DO 13 L=5,8
   13   CC(L,3)=(-(L-2.5)*CC(L-1,3)+2.*CC(L-2,2)+PH*CC(L-2,3))/P2
        CC(4,4)=(CST(4)-3.*CST(3)/PSI+3.*CST(2)/P2)/(16.*PHI*P3)
        CC(5,4)=(CST(4)-6.*CST(3)/PSI+15.*CST(2)/P2-15.*CST(1)/P3)
     1   /(16.*P4)
        DO 14 L=6,8
   14   CC(L,4)=(-(L-2.5)*CC(L-1,4)+3.*CC(L-2,3)+PH*CC(L-2,4))/P2
        CC(5,5)=(-CST(5)+6.*CST(4)/PSI-15.*CST(3)/P2+15.*CST(2)/P3)
     1	   /(32.*PHI*P4)
        CC(6,5)=(-CST(5)+10.*CST(4)/PSI-45.*CST(3)/P2+105.*CST(2)/P3
     1	    -105.*CST(1)/P4)/(32.*P5)
        DO 15 L=7,9
   15   CC(L,5)=(-(L-2.5)*CC(L-1,5)+4.*CC(L-2,4)+PH*CC(L-2,5))/P2
        GO TO 150
  140   CC(1,1)=-CST(1)/PHI
        CC(2,1)=-CST(2)
        DO 16 L=3,7
   16   CC(L,1)=-A*CC(L-1,1)/(FLOAT(L)-1.5)
        CC(2,2)=CST(3)/(2.*PHI)
        CC(3,2)=CST(4)/6.
        DO 17 L=4,7
   17   CC(L,2)=(CC(L-1,1)-A*CC(L-1,2))/(FLOAT(L)-1.5)
        DO 18 L=6,10
   18   CST(L)=-2.*(L-2)*CST(L-2)-2.*ZE*CST(L-1)
        CC(3,3)=-CST(5)/(12.*PHI)
        CC(4,3)=-CST(6)/60.
        DO 19 L=5,8
   19   CC(L,3)=(2.*CC(L-1,2)-A*CC(L-1,3))/(FLOAT(L)-1.5)
        CC(4,4)=CST(7)/(120.*PHI)
        CC(5,4)=CST(8)/840.
        DO 20 L=6,8
   20   CC(L,4)=(3.*CC(L-1,3)-A*CC(L-1,4))/(FLOAT(L)-1.5)
        CC(5,5)=-CST(9)/(1680.*PHI)
        CC(6,5)=-CST(10)/15120.
        DO 21 L=7,9
   21   CC(L,5)=(4.*CC(L-1,4)-A*CC(L-1,5))/(FLOAT(L)-1.5)
  150   DO 41 II=2,6
        DO 41 J=1,5
        I=II+1
        IF((J.EQ.3).OR.(J.EQ.4)) I=II+2
        IF(J.EQ.5) I=II+3
        CAB=CABS(CC(I,J))
        CAB=ALOG(CAB)
        TOT=-XA*XA+YA*YA+CAB
        IF(TOT.GT.AFACT) GO TO 32
        CC(I,J)=CC(I,J)*CEXPP(-ZSQ)
        GO TO 33
   32   CC(I,J)=CC(I,J)*CMPLX(SR,-SI)
        CCA=CABS(CC(I,J))
        CC(I,J)=CC(I,J)/CCA
        CC(I,J)=CC(I,J)*CFACT
   33   CONTINUE
   41   CONTINUE
        DO 45 I=2,6
         V(I+1)= V(I+1)+CC(I+1,1)
        DV(I+1)=DV(I+1)+CC(I+1,2)
        D2V(I+2)=D2V(I+2)+CC(I+2,3)
        D3V(I+2)=D3V(I+2)+CC(I+2,4)
        D4V(I+3)=D4V(I+3)+CC(I+3,5)
   45   CONTINUE
        IF(P2.LT.XXX) GO TO 50
   46   CONTINUE
   50   CONTINUE
        RETURN
        END
C
	SUBROUTINE FQ(A,B)
C
C Calculates the slightly relativistic plasma dispersion Fq and up to its 
C fourth derivative as a function of A= mu*delta and B=kpar*c**2/W*Vth. 
C Both real and Imaginary parts are obtained. 
C For small arguments DZ and DD are used.
C
	COMMON/VVV/V,DV,D2V,D3V,D4V
	COMMON/CCCC/CFACT
	COMPLEX V(7),DV(7),D2V(8),D3V(8),D4V(9),CMPLX,DZ,PHI,X(2),Z(2)
     +  ,ZP,ZM,D1(2),D1P,D1M,D2(2),D2P,D2M,D3(2),D3P,D3M,D4(2),D4P,D4M
     +  ,ZO
	ZERO=0.0
	IF(ABS(A).LT.100.) GO TO 1
	IF(ABS(B/A).LT.0.038) GO TO 500
   1	PSI=B/SQRT(2.)
	P2=B*B/2.
	P3=P2*PSI
	P4=P2*P2
	P5=P4*PSI
	PH=P2-A
	BB=B
	IF(P2.LT.(5.8E-3/ALOG(CFACT))) BB=ZERO
	ZO=CMPLX(ZERO,ZERO)
	IF(PH.EQ.ZERO) GO TO 19
	IF(PH.LT.ZERO) GO TO 3
	PHI=CMPLX(SQRT(PH),ZERO)
	GO TO 4
    3	RI=-SQRT(-PH)
	PHI=CMPLX(ZERO,RI)
    4	IF(BB.EQ.ZERO) GO TO 29
	X(1)=PSI-PHI
	X(2)=-PSI-PHI
	CALL DD(X,Z,D1,D2,D3,D4)
	ZP=Z(1)+Z(2)
	ZM=Z(1)-Z(2)
	IF(PH.GT.ZERO) GO TO 50
	ZP=CMPLX(ZERO,AIMAG(ZP))
	ZM=CMPLX(REAL(ZM),ZERO)
   50	V(1)=-ZP/(2.*PHI)
	V(2)=-ZM/(2.*PSI)
	DO 10 I=3,7
   10	V(I)=(1.+PH*V(I-2)-(I-2.5)*V(I-1))/P2
	D1P=D1(1)+D1(2)
	D1M=D1(1)-D1(2)
	IF(PH.GT.ZERO) GO TO 51
	D1P=CMPLX(REAL(D1P),ZERO)
	D1M=CMPLX(ZERO,AIMAG(D1M))
   51	DV(2)=D1M/(4.*PSI*PHI)
	DV(3)=(D1P-ZM/PSI)/(4.*P2)
	DO 11 I=4,7
   11	DV(I)=(-(I-2.5)*DV(I-1)+V(I-2)+PH*DV(I-2))/P2
	D2P=D2(1)+D2(2)
	D2M=D2(1)-D2(2)
	IF(PH.GT.ZERO) GO TO 52
	D2P=CMPLX(ZERO,AIMAG(D2P))
	D2M=CMPLX(REAL(D2M),ZERO)
   52	D2V(3)=(-D2P+D1M/PSI)/(8.*PHI*P2)
	D2V(4)=(-D2M+3.*D1P/PSI-3.*ZM/P2)/(8.*P3)
	DO 12 I=5,8
   12	D2V(I)=(-(I-2.5)*D2V(I-1)+2.*DV(I-2)+PH*D2V(I-2))/P2
	D3P=D3(1)+D3(2)
	D3M=D3(1)-D3(2)
	IF(PH.GT.ZERO) GO TO 53
	D3P=CMPLX(REAL(D3P),ZERO)
	D3M=CMPLX(ZERO,AIMAG(D3M))
   53	D3V(4)=(D3M-3.*D2P/PSI+3.*D1M/P2)/(16.*PHI*P3)
	D3V(5)=(D3P-6.*D2M/PSI+15.*D1P/P2-15.*ZM/P3)/(16.*P4)
	DO 13 I=6,8
   13	D3V(I)=(-(I-2.5)*D3V(I-1)+3.*D2V(I-2)+PH*D3V(I-2))/P2
	D4P=D4(1)+D4(2)
	D4M=D4(1)-D4(2)
	IF(PH.GT.ZERO) GO TO 54
	D4P=CMPLX(ZERO,AIMAG(D4P))
	D4M=CMPLX(REAL(D4M),ZERO)
   54	D4V(5)=(-D4P+6.*D3M/PSI-15.*D2P/P2+15.*D1M/P3)/(32.*PHI*P4)
	D4V(6)=(-D4M+10.*D3P/PSI-45.*D2M/P2+105.*D1P/P3-105.*ZM/P4)/
     1  (32.*P5)
	DO 14 I=7,9
   14	D4V(I)=(-(I-2.5)*D4V(I-1)+4.*D3V(I-2)+PH*D4V(I-2))/P2
	RETURN
   19	IF(BB.EQ.ZERO) GO TO 39
        IF(ABS(A).GT.20.0) GO TO 500
	X(1)=PSI
	X(2)=-PSI
	CALL DD(X,Z,D1,D2,D3,D4)
	ZM=Z(1)-Z(2)
	ZM=CMPLX(REAL(ZM),ZERO)
	V(2)=-ZM/(2.*PSI)
	DO 20 I=3,7
   20	V(I)=(1.-(I-2.5)*V(I-1))/P2
	D1P=D1(1)+D1(2)
	D1P=CMPLX(REAL(D1P),ZERO)
	DV(3)=(D1P-ZM/PSI)/(4.*P2)
	DO 21 I=4,7
   21	DV(I)=(-(I-2.5)*DV(I-1)+V(I-2))/P2
	D2M=D2(1)-D2(2)
	D2M=CMPLX(REAL(D2M),ZERO)
	D2V(4)=(-D2M+3.*D1P/PSI-3.*ZM/P2)/(8.*P3)
	DO 22 I=5,8
   22	D2V(I)=(-(I-2.5)*D2V(I-1)+2.*DV(I-2))/P2
	D3P=D3(1)+D3(2)
	D3P=CMPLX(REAL(D3P),ZERO)
	D3V(5)=(D3P-6.*D2M/PSI+15.*D1P/P2-15.*ZM/P3)/(16.*P4)
	DO 23 I=6,8
   23	D3V(I)=(-(I-2.5)*D3V(I-1)+3.*D2V(I-2))/P2
	D4M=D4(1)-D4(2)
	D4M=CMPLX(REAL(D4M),ZERO)
	D4V(6)=(-D4M+10.*D3P/PSI-45.*D2M/P2+105.*D1P/P3-105.*ZM/P4)/
     1 	 (32.*P5)
	DO 24 I=7,9
   24	D4V(I)=(-(I-2.5)*D4V(I-1)+4.*D3V(I-2))/P2
	GO TO 45
   29	V(1)=-DZ(-PHI,0)/PHI
	V(2)=-DZ(-PHI,1)
	DO 30 I=3,7
   30	V(I)=(1.-A*V(I-1))/(FLOAT(I)-1.5)
	DV(2)=DZ(-PHI,2)/(2.*PHI)
	DV(3)=DZ(-PHI,3)/6.
	DO 31 I=4,7
   31	DV(I)=(V(I-1)-A*DV(I-1))/(FLOAT(I)-1.5)
	D2V(3)=-DZ(-PHI,4)/(12.*PHI)
	D2V(4)=-DZ(-PHI,5)/60.
	DO 32 I=5,8
   32	D2V(I)=(2.*DV(I-1)-A*D2V(I-1))/(FLOAT(I)-1.5)
	D3V(4)=DZ(-PHI,6)/(120.*PHI)
	D3V(5)=DZ(-PHI,7)/840.
	DO 33 I=6,8
   33	D3V(I)=(3.*D2V(I-1)-A*D3V(I-1))/(FLOAT(I)-1.5)
	D4V(5)=-DZ(-PHI,8)/(1680.*PHI)
	D4V(6)=-DZ(-PHI,9)/15120.
	DO 34 I=7,9
   34	D4V(I)=(4.*D3V(I-1)-A*D4V(I-1))/(FLOAT(I)-1.5)
	RETURN
   39	DO 40 I=2,7
   40	V(I)=1./(FLOAT(I)-1.5)
	DO 41 I=3,7
   41	DV(I)=V(I-1)/(FLOAT(I)-1.5)
	DO 42 I=4,8
   42	D2V(I)=2.*DV(I-1)/(FLOAT(I)-1.5)
	DO 43 I=5,8
   43	D3V(I)=3.*D2V(I-1)/(FLOAT(I)-1.5)
	DO 44 I=6,9
   44	D4V(I)=4.*D3V(I-1)/(FLOAT(I)-1.5)
   45	V(1)=ZO
	DV(2)=ZO
	D2V(3)=ZO
	D3V(4)=ZO
	D4V(5)=ZO
	RETURN
  500 	CALL FQ1(A,B)
	RETURN
        END 
C
	COMPLEX FUNCTION DZ(ZE,M)
       	COMPLEX ZE,ZSQ,SUM,TERM,FMULT,TERME,ZM,CC(10),CEXPP,CMPLX
	COMMON/BBBB/CC,SUM
	COMMON/CCCC/CFACT
	AFACT=ALOG(CFACT)
	BFACT=1.0
	EPS=1.0E-6
	X=REAL(ZE)
	Y=AIMAG(ZE)
	XA=ABS(X)
	YA=ABS(Y)
	AZS=XA*XA+YA*YA
	IF(AZS.LT.26.8) GO TO 300
   20	ZSQ=ZE*ZE
	FN=5.0
	SUM=CMPLX(1.,0.)
	TERM=SUM
	FMULT=0.5/ZSQ
   21	TERME=TERM
	TERM=TERM*FMULT*(FN+M-3.)*(FN+M-4.)/(FN-3.)
	ZM=TERM/TERME
	IF(CABS(ZM).GT.1.0) GO TO 22
	SUM=SUM+TERM
	FN=FN+2.
	RET=ABS(REAL(TERM))
	RIT=ABS(AIMAG(TERM))
	IF((RET.GT.(EPS*ABS(REAL(SUM)))).AND.(RET.NE.0.)) GO TO 21
	IF((RIT.GT.(EPS*ABS(AIMAG(SUM)))).AND.(RIT.NE.0.)) GO TO 21
   22	SUM=FAC(M)*SUM/((-ZE)**(M+ 1))
	IF(Y.GT.0.) GO TO 31
	CS=COS(2.*X*Y)
	SS=SIN(2.*X*Y)
	SR=SIGN(BFACT,CS)
	SI=SIGN(BFACT,SS)
	IF(CS.EQ.0.0) SR=0.0
	IF(SS.EQ.0.0) SI=0.0
	CC(1)=CMPLX(0.0,3.544907702)
	IF(Y.EQ.0.) CC(1)=0.5*CC(1)
        CC(2)=-2.*ZE*CC(1)
	DO 23 L=3,M+1
   23   CC(L)=-2.*(L-2)*CC(L-2)-2.*ZE*CC(L-1)
	DO 32 L=1,M+1
	CAB=CABS(CC(L))
	CAB=ALOG(CAB)
	TOT=-XA*XA+YA*YA+CAB
	IF(TOT.GT.AFACT) GO TO 33
        CC(L)=CC(L)*CEXPP(-ZSQ)
	GO TO 32
   33	CC(L)=CC(L)*CMPLX(SR,-SI)
	CCA=CABS(CC(L))
	CC(L)=CC(L)/CCA
	CC(L)=CC(L)*CFACT
   32	CONTINUE
   	DZ=SUM+CC(M+1)
	RETURN
   31	DZ=SUM
	RETURN
  300   IF((AZS.GE.(16.+M*2.7)).AND.(M.GT.4)) GO TO 20	     
    	CALL ZZ(X,Y,RE,RI)
	CC(1)=CMPLX(RE,RI)
	CC(2)=-2.*(1.+ZE*CC(1))
	IF(M.LT.2) GO TO 320
	DO 310 L=3,M+1
  310	CC(L)=-2.*(L-2)*CC(L-2)-2.*ZE*CC(L-1)
  320	DZ=CC(M+1)
	RETURN
	END
C
	COMPLEX FUNCTION CEXPP(Z)   
Calculates the exponential of a complex argument
	COMPLEX CMPLX,ZA,CEXP,Z
	COMMON/CCCC/CFACT
	AFACT=ALOG(CFACT)
	RZ=REAL(Z)
	AZ=AIMAG(Z)
	IF(RZ.GT.AFACT) RZ=AFACT
	IF(RZ.LT.-AFACT) RZ=-AFACT
	ZA=CMPLX(RZ,AZ)
	CEXPP=CEXP(ZA)
	RETURN
	END
C
	FUNCTION EXPP(Z)           ! calculates exponential for a real argument
	COMMON/CCCC/CFACT
	AFACT=ALOG(CFACT)
	IF(Z.GT.AFACT) Z=AFACT
	IF(Z.LT.-AFACT) Z=-AFACT
	EXPP=EXP(Z)
	RETURN
	END
C
	FUNCTION FAC(N)             ! calculates the factorial function.
	F=1.
	IF(N.LE.1) GO TO 20
	DO 10 I=2,N
   10	F=I*F
   20	FAC=F
	RETURN
	END
C
	SUBROUTINE ZZ(X,Y,ZZR,ZZI)
C
C It calculates the Real ZZR and Imaginary ZZI parts of the 
C plasma dispersion function Z(X+iY) by calling WZ.
C
	X1=ABS(X)
	Y1=ABS(Y)
        CT=1.7724538509055
	CALL WZ(X1,Y1,WZR1,WZI1)
	IF((X.GE.0.).AND.(Y.GE.0.)) GO TO 1
	IF((X.LE.0.).AND.(Y.GE.0.)) GO TO 2
	A=2.*X1*Y1
	B=-(X1*X1-Y1*Y1)
	ABR=2.*EXPP(B)*COS(A)
	ABI=-2.*EXPP(B)*SIN(A)
	WZR1=ABR-WZR1
	WZI1=ABI-WZI1
	IF((X.LE.0.).AND.(Y.LE.0.)) GO TO 1
    2	WZI1=-WZI1
    1	ZZR=-CT*WZI1
	ZZI=CT*WZR1
	RETURN
	END
C
	SUBROUTINE WZ(X,Y,RE,IM)
C
C It calculates the Real RE and Imaginary parts of Z/(i*sqrt(pi))
C
	EQUIVALENCE(EPSH,EPSL,EPSY)
	REAL IM,LAMBDA
	INTEGER CAPN
	LOGICAL B
	EPSH=1.E-12
	IF((Y.LT.4.29).AND.(X.LT.5.33)) GO TO 10
	H=0.
	CAPN=0
	NU=8
	GO TO 20
   10	S=(1.-Y/4.29)*SQRT(1.-X*X/28.41)
	H=1.6*S
	H2=2.*H
	CAPN=6.+23.*S+.5
	NU=9.+21.*S+.5
	LAMBDA=H2**CAPN
   20	B=((H.EQ.0.).OR.(LAMBDA.LT.EPSL))
   	RR=0.
	RI=0.
	SR=0.
	SI=0.
	NUP=NU+1
	DO 100 I=1,NUP
	N=NUP-I
	NP1=N+1
	TR=Y+H+NP1*RR
	TI=X-NP1*RI
	C=.5/(TR*TR+TI*TI)
	RR=C*TR
	RI=C*TI
	IF(.NOT.((H.GT.0.).AND.(N.LE.CAPN))) GO TO 100
	TR=LAMBDA+SR
	SR=RR*TR-RI*SI
	SI=RI*TR+RR*SI
	LAMBDA=LAMBDA/H2
  100	CONTINUE
	CC=1.12837916709551
	IF(Y.LT.EPSY) GO TO 120
	IF(B) GO TO 110
	RE=SR*CC
	GO TO 130
  110	RE=RR*CC
	GO TO 130
  120	RE=EXPP(-X*X)
  130	IF(B) GO TO 140
	IM=SI*CC
	GO TO 150
  140	IM=RI*CC
  150	RETURN
	END		 
C
	SUBROUTINE DD(X,Z,D,D2,D3,D4)
C
C It provides for two X arguments, the two values of Z and their first four 
C D,D2,D3,D4 derivatives needed to give Fq when the argument is not too large 
C and obtained by calling subroutine DZ. 
C For large arguments, the order of derivation is Z'''',Z''',Z'',Z' and Z. 
C For small arguments one first obtain Z and then calculates the derivatives.
C
	COMPLEX X(2),XS,Z(2),D(2),D2(2),D3(2),D4(2),DZ,CC(10),SUM
	COMMON/BBBB/CC,SUM
	DO 8 K=1,2
	L=1
	RI=AIMAG(X(K))
	XS=X(K)*X(K)
	IF((CABS(X(K))**2).LT.26.8) GO TO 7
	D4(K)=DZ(X(K),4)
	IF(RI.GT.0.0) GO TO 5
	D4(K)=SUM
	L=0
    5	D3(K)=(24.-X(K)*(2.*XS-3.)*D4(K))/(3.-12.*XS+4.*XS*XS)
	D2(K)=-(4.+(XS-0.5)*D3(K))/(X(K)*(2.*XS-3.))
	D(K)=(2.-X(K)*D2(K))/(2.*XS-1.)
	Z(K)=-(1.+0.5*D(K))/X(K)
	IF(L.EQ.0) GO TO 6
	GO TO 8
    6	D4(K)=D4(K)+CC(5)
	D3(K)=D3(K)+CC(4)
	D2(K)=D2(K)+CC(3)
	D(K)=D(K)+CC(2)
	Z(K)=Z(K)+CC(1)
	GO TO 8
    7	Z(K)=DZ(X(K),0)
	D(K)=-2.*(1.+X(K)*Z(K))
	D2(K)=-2.*Z(K)-2.*X(K)*D(K)
	D3(K)=-4.*D(K)-2.*X(K)*D2(K)
	D4(K)=-6.*D2(K)-2.*X(K)*D3(K)
    8   CONTINUE	 
    	RETURN
	END
C
	SUBROUTINE EPOLAR(G,RP,E,BM,BMAG,WK,WKDTB,ER,ERAT,E3)
C
C Given G = Npar, RP = Nperp,E(6,9) from subroutine DF, BM = B/Bo,
C BMAG = |B|/Bo, Wk = k and WKDTB = Kpar|B|/Bo, one obtains ER, the complex 
C components of the electric polarization vector in a coordinate system 
C with the third component along B, ERAT, the complex component in the 
C (x,y,z) coordinate system and E3 = |Ex|,|Ey|,|Ez|. 
C
	COMPLEX E(6,9),ETERM(4),
     1	 ERATIO(6),ER(3),ERAT(3),AII,CMPLX
	DIMENSION EXH(3),EYH(3),EZH(3),E3(3),WK(3),
     1	  BM(3)
C             ,GRB(3),GRD(3),GRT(3),GRBIJ(3,3)
	COMMON/ALL/IBUG,IQ,IIO,CFE0,PFE0,T0,SGN,PI,WLD,IIQ,ND125
	G2=G*G
	RP2=RP*RP
	GRP=G*RP
	RPN2=G2+RP2
	ETERM(1)=E(1,1)-G2
	ETERM(2)=E(3,1)+GRP
	ETERM(3)=E(4,1)-RPN2
	ETERM(4)=E(6,1)-RP2
	ERATIO(1)=E(5,1)**2+ETERM(3)*ETERM(4)
	ERATIO(2)=E(5,1)*ETERM(2)+E(2,1)*ETERM(4)
	ERATIO(3)=E(2,1)*E(5,1)-ETERM(2)*ETERM(3)
	ERATIO(4)=-(ETERM(1)*E(5,1)+E(2,1)*ETERM(2))
	ERATIO(5)=E(2,1)**2+ETERM(1)*ETERM(3)
	ERATIO(6)=ETERM(2)**2-ETERM(1)*ETERM(4)
	DO 20 I=1,3
   20	EZH(I)=BM(I)/BMAG
	DO 22 I=1,3
	J1=I+1
	J2=I+2
	IF(J1.GT.3) J1=J1-3
	IF(J2.GT.3) J2=J2-3
   22	EYH(I)=BM(J1)*WK(J2)-BM(J2)*WK(J1)
	ARG=EYH(1)**2+EYH(2)**2+EYH(3)**2
	IF(ARG.EQ.0.) GO TO 27
	EYMAG=SQRT(ARG)
	GO TO 28
   27	EYH(1)=0.0
	ARG1=SQRT(EZH(1)**2+EZH(3)**2)
	ARG2=SQRT(EZH(2)**2+EZH(3)**2)
	EYH(2)=EZH(3)/ARG2
	EYH(3)=-EZH(2)/ARG2
	EXH(1)=EZH(3)/ARG1
	EXH(2)=0.0
	EXH(3)=-EZH(1)/ARG1
	GO TO 33
   28	EXMAG=0.0
	DO 30 I=1,3
	EXH(I)=WK(I)*(BMAG**2)-BM(I)*WKDTB
   30	EXMAG=EXMAG+EXH(I)**2
	EXMAG=SQRT(EXMAG)
	DO 31 I=1,3
	EYH(I)=EYH(I)/EYMAG
   31	EXH(I)=EXH(I)/EXMAG
   33	CONTINUE
	AII=CMPLX(0.0,-1.0)
	IF((ABS(G).LT.1.E-3).OR.(ABS(1.0-RPN2).LT.0.1)) GO TO 40
	DO 34 I=1,3
   34	ER(I)=ERATIO(I)
	GO TO 60
   40	CONTINUE
   	IF(SGN.LT.0.) GO TO 51
	DO 41 I=1,3
   41	ER(I)=-ERATIO(I+2)
	GO TO 60
   51	ER(1)=ERATIO(2)*AII
	ER(2)=ERATIO(6)*AII
	ER(3)=ERATIO(4)*AII
   60	CONTINUE
	DO 70 I=1,3
   70	ERAT(I)=ER(1)*EXH(I)+ER(2)*EYH(I)+ER(3)*EZH(I)
	DO 80 I=1,3
   80	E3(I)=(REAL(ERAT(I)))**2+(AIMAG(ERAT(I)))**2
	PMAG=SQRT(E3(1)+E3(2)+E3(3))
	DO 90 I=1,3
	ER(I)=ER(I)/PMAG
	ERAT(I)=ERAT(I)/PMAG
   90	E3(I)=(SQRT(E3(I)))/PMAG
       	RETURN
	END
C

      double precision function
     * dshkarof(xe,te,ye,cnpar,cnper,friq,hami)
c-------------------------------------------------------------------
c     calculates real and imaginary parts of Shkarofsky dispersion function
c     using the given xe,te,ye,cnpar,cnper
c--------------------------------------------------------------------
      double precision xe,te,ye,cnpar,cnper,friq,hamr,hami

cBobH990117      REAL X,Y,LR,G,U,WLD,RP,BMAG,DY,TEMP,DR,DI
      REAL X,Y,G,U,WLD,RP,BMAG,DY,TEMP,DR,DI
      REAL WFMU,FMU,T0
      REAL WK(3),BV(3),GRB(3),GRD(3),GRT(3),GRBIJ(3,3)
      COMPLEX DP,DH,DK,DW,DX(3),DKI(3),E(6,9)
      REAL BM(3)
      REAL PI
cSm030426
c     to put the dielectric tensor
      include 'eps.i'
      integer i,j

      PI=4.0*atan(1.0)
      LR=0
      X=xe
      Y=ye
      G=cnpar
      FMU=friq   !GHZ
      WFMU = 2*PI*FMU
      WLD= 29.979246/WFMU   ! c/W   2.997 10 e10 cm/s  // fmu 10 e9
      TEMP=te*1000.	    ! TEMP     kev
      T0=1.		    ! T0  (ev) ?
      U=1./(TEMP*T0*1.95693415E-6)
      RP=cnper
      DY=1.		    ! density ?
      BMAG=1.
      BM(1)=0.
      BM(2)=0.
      BM(3)=1.
      write(*,*)'shkarof X,Y,TEMP',X,Y,TEMP
      write(*,*)'shkarof G,RP',G,RP
      call DF(X,Y,LR,G,U,WLD,WK,RP,BMAG,DY,TEMP,BM,GRB,GRD,
     1	 GRT,GRBIJ,E,DR,DI,DP,DH,DK,DW,DX,DKI)
cSm030514	
c      hamr=DR2
	hamr=DR
      hami=DI
      dshkarof=hamr
      WRITE(*,*)'dshkarof DR,DI',DR,DI

cSm030426
c-----put the dielectric tensor to reps(3,3),It will be in eps.i
 
      reps(1,1)=dcmplx(E(1,1))
      reps(1,2)=dcmplx(E(2,1))
      reps(1,3)=dcmplx(E(3,1))

      reps(2,1)=-reps(1,2)
      reps(2,2)=dcmplx(E(4,1))
      reps(2,3)=dcmplx(E(5,1))

      reps(3,1)=reps(1,3)
      reps(3,2)=-reps(2,3)
      reps(3,3)=dcmplx(E(6,1))
c---------------------------------------------

      return
      end

      subroutine absorpsh(cnpar,cnper,xe,ye,te_kev,friq,cnprim)
c-------------------------------------------------------
c     absorbtion for Sharofsky tensor
c     calculation of the numerical derivative from the
c     real part of the hermitian dispersion function
c     d(dispersion_function)/d(N_perp_real)
c-------------------------------------------------------
      implicit none
c-----input
      double precision cnpar,cnper,xe,ye,te_kev,friq
      external dshkarof 
      double precision dshkarof 
c-----output
      double precision cnprim
c-----local
      double precision step,cnperp,cnperm,hamr,hamrp,hamrm,dhamrdnr,
     *hamimag,hami,hamip,hamim
      double complex chamilt

      hamr=dshkarof(xe,te_kev,ye,cnpar,cnper,friq,hamimag)    
        
      step=1.0d-7
      if(dabs(cnper).gt.1.d0) step=cnper*step
            
      cnperp=cnper+step
      hamrp=dshkarof(xe,te_kev,ye,cnpar,cnperp,friq,hamip)    
               
      cnperm=cnper-step
      hamrm=dshkarof(xe,te_kev,ye,cnpar,cnperm,friq,hamim)    
                
      dhamrdnr=(hamrp-hamrm)/(2.d0*step)
      cnprim=-hamimag/dhamrdnr ! imaginary part of N_perp

      write(*,*)'shkarofsk: hamimag,dhamrdnr,cnprim',
     *hamimag,dhamrdnr,cnprim
      
      return
      end

      subroutine efiedshk (cnpar,cnper,xe,ye,te_kev,friq,cex,cey,cez,
     *ex,ey,ez)
c-------------------------------------------------------
c     calculation of the electric field for Sharofsky tensor
c     electric field polarisations: ex,ey,ex
c     using the given xe,te,ye,real(cnpar),real(cnper)
      implicit double precision (a-h,o-z)
c-----input
      double precision cnpar,cnper,xe,ye,te_kev,friq
c-----output      
      double precision ex,ey,ez
      double complex cex,cey,cez

cBobH990117      REAL X,Y,LR,G,U,WLD,RP,BMAG,DY,TEMP,DR,DI
      REAL X,Y,G,U,WLD,RP,BMAG,DY,TEMP,DR,DI
      REAL WFMU,FMU,T0
      REAL WK(3),BV(3),GRB(3),GRD(3),GRT(3),GRBIJ(3,3)
      COMPLEX DP,DH,DK,DW,DX(3),DKI(3),E(6,9)
      REAL BM(3)
      REAL PI
            
cSmirnov970321
      REAL E3(3),WKDTB
      COMPLEX ER(3),ERAT(3)
      

      PI=4.0*atan(1.0)
      LR=0
      X=xe
      Y=ye
      G=cnpar
      FMU=friq   !GHZ
      WFMU = 2*PI*FMU
      WLD= 29.979246/WFMU         ! c/W   2.997 10 e10 cm/s  // fmu 10 e9
      TEMP=te_kev*1000.d0 	!| TEMP     ev
      T0=1.		!| T0  (ev) ?
      U=1./(TEMP*T0*1.95693415E-6)
      RP=cnper
      DY=1.		!| density ?
cSmirnov970321
      BMAG=1.
      BM(1)=0.
      BM(2)=0.
      BM(3)=1.
      call DF(X,Y,LR,G,U,WLD,WK,RP,BMAG,DY,TEMP,BM,GRB,GRD,
     1 GRT,GRBIJ,E,DR,DI,DP,DH,DK,DW,DX,DKI)
      CALL EPOLAR(G,RP,E,BM,BMAG,WK,WKDTB,ER,ERAT,E3)
      
      ex=E3(1)
      ey=E3(2)
      ez=E3(3)
      cex=ER(1)
      cey=ER(2)
      cez=ER(3)

      write(*,*)'shkarof ERAT(1)',ERAT(1)
      write(*,*)'shkarof ERAT(2)',ERAT(2)
      write(*,*)'shkarof ERAT(3)',ERAT(3)

      write(*,*)'shkarof cex',cex
      write(*,*)'shkarof cey',cey
      write(*,*)'shkarof cez',cez

      WRITE(*,*)'shkarof ex,ey,ez',ex,ey,ez

      return
      end











 
