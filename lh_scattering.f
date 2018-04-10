!     Nicola Bertelli, PPPL, 6 April 2012
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
      MODULE parameters_lhscat
      
      IMPLICIT NONE  

!      INCLUDE 'param.i'
!      INCLUDE 'one.i'

      DOUBLE PRECISION, PARAMETER:: savg = 0.1D0  !!, zcorr = 66.0  !247.5D0  [dimensionless being multiplied by "a" (minor radius)]
      DOUBLE PRECISION, PARAMETER:: zcorr = 3.14D0 ! [cm^{-1}] given by zcorr = 2 * pi / lambda, where lambda ~ 2 cm see Wallace PhD thesis p. 185
      DOUBLE PRECISION, PARAMETER:: valid = 1.D0, dcx = 272.25D0
      DOUBLE PRECISION, PARAMETER:: deltan0 = 0.3333333
      DOUBLE PRECISION, PARAMETER:: clight = 2.99792458D10  ! [cm/s]
!      DOUBLE PRECISION, PARAMETER:: omlh = 2.D0 * pi * frqncy

      END MODULE parameters_lhscat
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------





!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

      SUBROUTINE lhscattering(u, s_total_new, s_total_old)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    THIS SUBROUTINE EVALUATES THE LOWER HYBRID SCATTERING IN THE EDGE PLASMA              !!
!!    ACCORDING TO THE MODEL DESCRIBED IN P. BONOLI AND E. OTT Phys. Fluids 25              !!
!!    (1982) 359. SEE ALSO P. BONOLI'S PhD thesis and PRL (PRL, 46 (1981) 424)              !!
!!    SEE ALSO: OTT, PHYS. FLUIDS 22 (1979) 1732                                            !!
!!              OTT ET AL., PHYS. FLUIDS 23 (1980) 1031                                     !!
!!              HUI ET AL., NF 21 (1981) 339                                                !!
!!              VAHALA ET AL., PHYS. FLUIDS B 4 (1992) 4033.                                !! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!---------------------------------------------------------------------------------------------
!     INPUT & OUTPUT OF THIS SUBROUTINE
!     INPUT  = array "u" and array "prmt"
!     OUTPUT = array "u" and array "prmt"
!
!     In particular, using the GENRAY's labels,
!     u(1) = z, u(2) = r, u(3) = phi, u(4) = cnz , u(5) = cnr, u(6) = cm 
!     where z = vertical axis, r = major radius, phi = toroidal angle
!     cnz = z-component of the refractive index
!     cnr = r-component of the refractive index
!     cm  = r * N_phi = major radius X the phi-component of the refractive index.
!  
!     (r, phi, z) are cylindrical space coordinates 
!     (cnr, cm, cnz) are conjugate coordinates for the refractive index
!      See GENRAY manual (www.compxco.com/Genray_manual.pdf) page 6.
!
!  
!     prmt is a array of dimension 9 but in this subroutine we need only the component prmt(6)
!     prmt(6) = [meter] poloidal or total distance for results output (see genray.in in the
!     "numercal" namelist).
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
!                            !!!!!!  IMPORTANT  !!!!!!!
!     in this subroutine we need the TOTAL distance in order to advance the rays
!     correctly. In other words, when we run GENRAY including this scattering model we need
!     to set "i_output = 2"  in the namelist "numercal" and so, first, "irkmeth = 2".   
!---------------------------------------------------------------------------------------------



      USE parameters_lhscat  
      IMPLICIT NONE
cSAP120511
c      INCLUDE 'F90DIR/param.i'
c      INCLUDE 'F90DIR/one.i'
c      INCLUDE 'F90DIR/rho.i'
      INCLUDE 'param.i'
      INCLUDE 'one.i'
      INCLUDE 'rho.i'

      DOUBLE PRECISION, DIMENSION(6), INTENT(INOUT):: u
      DOUBLE PRECISION, INTENT(INOUT):: s_total_new
c SAP120517
c      DOUBLE PRECISION, INTENT(IN):: s_total_old   
      real*8, INTENT(INOUT):: s_total_old
!!      DOUBLE PRECISION, INTENT(INOUT):: us
!!      DOUBLE PRECISION, DIMENSION(9), INTENT(INOUT):: prmt
!!      DOUBLE PRECISION, INTENT(OUT):: ds
!      DOUBLE PRECISION:: ds
      DOUBLE PRECISION:: deltanx, dsyrms 
      DOUBLE PRECISION:: skp, skpm, rk
      DOUBLE PRECISION:: slength, func1, func2, func3
      DOUBLE PRECISION:: angle, rpaxy, rpaxx, sumfs, sumff, sumsf,sumss
      DOUBLE PRECISION:: probxx, probxy, probtot, xnum, xdum
      DOUBLE PRECISION:: rpa, deltat, tout, told
      DOUBLE PRECISION:: z, r, phi, wkz, wkr, wkph, omlh, al, amin, kvac
      DOUBLE PRECISION:: qr, qz, qph, wkn
      DOUBLE PRECISION:: xbeta, beta
      DOUBLE PRECISION:: nper2m, nper2p
cSAP120603
      DOUBLE PRECISION:: t, dtau, vgmod
      INTEGER:: rootmode

      DOUBLE PRECISION:: s_out, ds, delta_s_total, s_total_new_init
      DOUBLE PRECISION, DIMENSION(6) :: deru
cSAP120603
c      DOUBLE PRECISION, DIMENSION(1):: iden
      INTEGER, DIMENSION(1):: iden

      DOUBLE PRECISION, SAVE:: s_old = 0.d0

cSAP120604
      IF(s_total_old.lt.1.d-14) s_old=0.d0



!     OMEGA_LH = 2 * PI * FREQUENCY  ("frqncy" IS GIVEN IN "one.i" IN GHz)
      omlh = 2.D0 * pi * frqncy * 1.e+9


!     IN ORDER TO FIND THE RADIAL COORDINATE: ONE WAY IS SQRT( TOROIDAL FLUX / (PI * B0) )
!     THE EASIEST WAY IS SQRT(AREA) WHERE AREA = AREATOT. I USE THIS LAST ONE FOR THE FIRST MOMENT.
!     "torftot" (tesla*m^2) = TOROIDAL FLUX GIVEN IN "rho.i"
!     "areatot" (cm^2)= TOTAL CROSS-SECTIONAL AREA OF LCFS GIVEN IN "rho.i"
!     SEE BOB'S EMAIL ON APRIL 4 2012.

      amin = sqrt( areatot ) * 100.d0 

      kvac = omlh / clight  ! [cm^{-1}]
      al = clight / (omlh * amin)


!     IDENTIFY THE WAVE ROOTS OF THE DISPERSION RELATION (n-perp**2): 
!     IF WE ARE IN THE SLOW ("nper2m") OR FAST MODE ("nper2p")
!     "rootmode" IS A FLAG TO DISTINGUISH THE TWO MODES
!     FAST MODE ("rootmode = 1") & SLOW MODE ("rootmode = 2")  
 

      call wsort (u, nper2m, nper2p, rootmode)





!     CALCULATE THE SPATIAL PROFILE OF DENSITY FLUCTUATIONS.
!
!     "deltanx" IS "\delta n/N" PROFILE AS A FUNCTION OF "rho" ("deltan0" AND "dcx" ARE CONSTANT
!     AND GIVEN IN THE MODULE PARAMETERS AND "rho" IS GIVEN IN "one.i")

      deltanx = deltan0 * EXP( - dcx * ( rho - 1.0 )**2 )
!!!      deltanx = deltan0 * EXP( - dcx * ( rho - 1.20 )**2 )
      deltanx = MAX( deltanx, 0.001d0 )
      dsyrms = deltanx**2

!!     if (rho.ge.0.8d0.and.rho.le.1.d0) then
!!         deltanx = 1.35*rho-1.05
!!      else 
!!         return
!!      endif   






!     "rootmode = 1" THE WAVE ROOT CORRESPONDS TO THE FAST MODE THEREFORE
!     THE POSSIBLE SCATTERING EVENTS ARE FAST-SLOW OR FAST-FAST

!     IF STATEMENT FOR THE THREE CASES
      IF( rootmode.EQ.1 ) THEN




!        FAST-SLOW WAVE SCATTERING EVENT



!        "skp" = PERPENDICULAR REFRACTIVE INDEX (FAST WAVE)
         skp = SQRT( nper2m )





!        "skpm" = PERPENDICULAR REFRACTIVE INDEX (slow mode) AFTER THE SCATTERING
         skpm = SQRT( nper2p )


!!         rk = ( skp / ( zcorr * al) )**2
         rk = ( skp * kvac / zcorr )**2 ![dimensionless]

!        IN "rk" THERE IS THE FACTOR "al" = c/(omega * a) IN ORDER TO OBTAIN
!        (i) k_perp FROM N_perp AND 
!        (ii) TO GIVE DIMENSION TO "zcorr" BEING DIMENSIONLESS BEING EQUAL  z_0 * a 


!        "func2" = THE SCATTERING LENGTH
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         slength = func2
         IF( rk.LT.1.0 ) rpaxy = slength * zcorr
         IF( rk.GE.1.0 ) rpaxy = slength * zcorr / rk
!!!!         IF( rpaxy.LT.valid ) RETURN

         IF( rpaxy.LT.valid ) THEN
            ds = 0.d0
         ENDIF
        IF( rpaxy.LT.valid ) RETURN 

!        "func1" = INTEGRAL OF THE PROBABILITY (FAST-SLOW) FOR THE TIME STEP
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         sumfs = func1


!        FAST-FAST WAVE SCATTERING EVENT


!        "skpm" = PERPENDICULAR REFRACTIVE INDEX (FAST MODE) AFTER THE SCATTERING 

         skpm = SQRT( nper2m )

!        "func2" = THE SCATTERING LENGTH
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         slength = func2
         IF( rk.LT.1.0 ) rpaxx = slength * zcorr
         IF( rk.GE.1.0 ) rpaxx = slength * zcorr / rk
!!!!         IF( rpaxx.LT.valid ) RETURN

         IF( rpaxx.LT.valid ) THEN
            ds = 0.d0
         ENDIF
         IF( rpaxx.LT.valid ) RETURN

         rpa = MIN( rpaxx, rpaxy )

!        "func1" = INTEGRAL OF THE PROBABILITY (FAST-FAST) FOR THE TIME STEP
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         sumff = func1

!        TIME STEP 
         deltat = savg / ( sumff + sumfs )
!         tout = told + deltat
!         dtau = MIN( dtau, deltat )
!         IF( t.LT.tout ) RETURN
!         told = t
!         dtau = deltat


!        CALL RSIDE1 IN ORDER TO GET THE GROUP VELOCITY (NORMALIZED TO THE SPEED OF LIGHT)
         CALL rside1(0, u, deru)
cSAP120511
c         vgmod = clight * sqrt( deru(1)**2 + deru(2)**2 + ( u(2) * deru(3) )**2 )
         vgmod = clight * sqrt( deru(1)**2 + deru(2)**2 +
     &                         ( u(2) * deru(3) )**2 )
!!!!         vgmod = sqrt( deru(1)**2 + deru(2)**2 + ( u(2) * deru(3) )**2 )

!        "ds" = THE SPACE STEP FOR ADVANCING THE RAY 
!        (FACTOR 100 BECAUSE "vgmod" AND "dtau" ARE IN CGS BUT "prmt(6)"
!        IS EXPRESSED IN METERS)

         ds = vgmod * deltat / 100.0  ![m]

!         s_total_new_init = s_total_new

!! CONDITIONS FOR SCATTERING 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         s_out = s_old + ds
         delta_s_total = s_total_new - s_total_old
         delta_s_total = MIN(delta_s_total, ds)
         IF (s_total_new .lt. s_out ) THEN
         RETURN
         ELSEIF (s_total_new .gt. s_out ) THEN
         s_old = s_total_new
         s_total_new = s_total_old + delta_s_total 
cSAP120517
         s_total_old= s_old
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!         us = us + ds

!!!!         ds = vgmod * dtau    
              
!!!!         prmt(6) = MIN( prmt(6), ds )
!!!!         us = us + ds
         print*,'ROOTMODE =1: dtau and ds', deltat, ds


!        PROBABILITY TO HAVE FAST-FAST SCATTERING
         probxx = deltat * sumff
!        PROBABILITY TO HAVE FAST-SLOW SCATTERING
         probxy = deltat * sumfs
!        TOTAL PROBABILITY
         probtot = probxx + probxy


!        "xnum" = RANDOM NUMBER
         CALL RANDOM_NUMBER( xnum )

!        IF "xnum" IS LARGER THAN THE TOTAL PROBABILITY ("probtot") THEN NO SCATTERING
!        IF "xnum" IS SMALLER THAN THE PROBABILITY TO HAVE FAST-FAST SCATTERING ("probxx")
!        THEN THE PERPENDICULAR REFRACTIVE INDEX AFTER THE SCATTERING IS EQUAL TO 
!        THE PERPENDICULAR REFRACTIVE INDEX OF THE FAST MODE
!        IF "xnum" IS LARGER THAN THE PROBABILITY TO HAVE FAST-FAST SCATTERING ("probxx")
!        AND SMALLER THAN THE TOTAL PROBABILITY THEN THE PERPENDICULAR REFRACTIVE INDEX
!        AFTER THE SCATTERING IS EQUAL TO THE PERPENDICULAR REFRACTIVE INDEX OF THE SLOW MODE


!         IF ( xnum.GT.probtot ) THEN
!            ds = 0 
!         ENDIF

         IF ( xnum.GT.probtot ) RETURN   ! NO SCATTERING

         skpm = skp 

         IF ( xnum.GT.probxx.AND.xnum.LE.probtot ) THEN
            skpm = SQRT( nper2p )
         ENDIF
   



!     "rootmode = 2" THE WAVE ROOT CORRESPONDS TO THE SLOW MODE AND
!     THE FAST ROOT IS EVANESCENT ("nper2m" < 0). THEREFORE THE
!     ONLY SCATTERING EVENT IS SLOW-SLOW 


      ELSEIF( rootmode.EQ.2.AND.nper2m.LE.0.0 ) THEN



!        SLOW-SLOW WAVE SCATTERING EVENT



!        "skp" = PERPENDICULAR REFRACTIVE INDEX (SLOW WAVE)
         skp = SQRT( nper2p )

!        "skpm" = PERPENDICULAR REFRACTIVE INDEX (SLOW MODE) AFTER THE SCATTERING
         skpm = skp


!!         rk = ( skp / ( zcorr * al) )**2
         rk = ( skp * kvac / zcorr )**2 ![dimensionless]

!        IN "rk" THERE IS THE FACTOR "al" = c/(omega * a) IN ORDER TO OBTAIN
!        (i) k_perp FROM N_perp AND 
!        (ii) TO GIVE DIMENSION TO "zcorr" BEING DIMENSIONLESS BEING EQUAL  z_0 * a 


!        "func2" = THE SCATTERING LENGTH
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         slength = func2
         IF( rk.LT.1.0 ) rpa = slength * zcorr
         IF( rk.GE.1.0 ) rpa = slength * zcorr / rk
!!!!         IF( rpa.LT.valid ) RETURN

!         IF( rpa.LT.valid ) THEN
!            ds = 0.d0
!         ENDIF
         IF( rpa.LT.valid ) RETURN   

!        "func1" = INTEGRAL OF THE PROBABILITY (SLOW-SLOW)FOR THE TIME STEP
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         sumss = func1

!        TIME STEP
 
         deltat = savg / sumss
!         tout = told + deltat

!         dtau = MIN( dtau, deltat )
!         IF( t.LT.tout ) RETURN
!         told = t
!         dtau = deltat


!        CALL RSIDE1 IN ORDER TO GET THE GROUP VELOCITY (NORMALIZED TO THE SPEED OF LIGHT)
         CALL rside1(0, u, deru)
cSAP120511
c         vgmod = clight * sqrt( deru(1)**2 + deru(2)**2 + ( u(2) * deru(3) )**2 )
         vgmod = clight * sqrt( deru(1)**2 + deru(2)**2 + 
     &                         ( u(2) * deru(3) )**2 )
!!!!         vgmod = sqrt( deru(1)**2 + deru(2)**2 + ( u(2) * deru(3) )**2 )

!        "ds" = THE SPACE STEP FOR ADVANCING THE RAY 
!        (FACTOR 100 BECAUSE "vgmod" AND "dtau" ARE IN CGS BUT "prmt(6)"
!        IS EXPRESSED IN METERS)
         ds = vgmod * deltat / 100.0  ![m]


!! CONDITIONS FOR SCATTERING 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         s_out = s_old + ds
         delta_s_total = s_total_new - s_total_old
         delta_s_total = MIN(delta_s_total, ds)
         IF (s_total_new .lt. s_out ) THEN
         RETURN
         ELSEIF (s_total_new .gt. s_out ) THEN
         s_old = s_total_new
         s_total_new = s_total_old + delta_s_total 
cSAP120517
         s_total_old= s_old
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!         us = us + ds

!!!!         ds = vgmod * dtau 

!!!!         prmt(6) = MIN( prmt(6), ds )
!!!!         us = us + ds
         print*,'ROOTMODE =2 & FAST < 0: dtau and ds', deltat, ds


!        PROBABILITY TO HAVE SLOW-SLOW SCATTERING
         probxx = deltat * sumss

!        "xnum" = RANDOM NUMBER
         CALL RANDOM_NUMBER( xnum )

!        IF "xnum" IS LARGER THAN THE PROBABILITY TO HAVE SLOW-SLOW SCATTERING ("probxx")
!        THEN NO SCATTERING OTHERWISE SLOW-SLOW SCATTERING

!         IF( xnum.GT.probxx ) THEN
!            ds = 0.d0
!         ENDIF   

         IF( xnum.GT.probxx ) RETURN




!     "rootmode = 2" THE WAVE ROOT CORRESPONDS TO THE SLOW MODE AND
!     THE FAST ROOT IS PROPAGATIVE ("nper2m" > 0). THEREFORE THE
!     SCATTERING EVENTS ARE SLOW-SLOW AND SLOW-FAST


      ELSEIF( rootmode.EQ.2.AND.nper2m.GT.0.0 ) THEN


!        SLOW-FAST WAVE SCATTERING EVENT



!        "skp" = PERPENDICULAR REFRACTIVE INDEX (SLOW WAVE)
         skp = sqrt( nper2p )

!        "skpm" = PERPENDICULAR REFRACTIVE INDEX (FAST MODE) AFTER THE SCATTERING
         skpm = sqrt( nper2m )


!!         rk = ( skp / ( zcorr * al) )**2
         rk = ( skp * kvac / zcorr )**2 ![dimensionless]

!        IN "rk" THERE IS THE FACTOR "al" = c/(omega * a) IN ORDER TO OBTAIN
!        (i) k_perp FROM N_perp AND 
!        (ii) TO GIVE DIMENSION TO "zcorr" BEING DIMENSIONLESS BEING EQUAL  z_0 * a 


!        "func2" = THE SCATTERING LENGTH
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         slength = func2
         IF( rk.LT.1.0 ) rpaxy = slength * zcorr
         IF( rk.GE.1.0 ) rpaxy = slength * zcorr / rk
!!!!         IF( rpaxy.lt.valid ) RETURN

!         IF( rpaxy.LT.valid ) THEN
!            ds = 0.d0
!         ENDIF
         IF( rpaxy.lt.valid ) RETURN

!        "func1" = INTEGRAL OF THE PROBABILITY (SLOW-FAST)FOR THE TIME STEP
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         sumsf = func1


!        SLOW-SLOW WAVE SCATTERING EVENT


!        "skpm" = PERPENDICULAR REFRACTIVE INDEX (SLOW MODE) AFTER THE SCATTERING
         skpm = sqrt( nper2p )

!        "func2" = THE SCATTERING LENGTH
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         slength = func2
         IF( rk.LT.1.0 ) rpaxx = slength * zcorr
         IF( rk.GE.1.0 ) rpaxx = slength * zcorr / rk
         IF( rpaxx.LT.valid ) RETURN

!         IF( rpaxx.LT.valid ) THEN
!            ds = 0.d0
!         ENDIF
         IF( rpaxx.lt.valid ) RETURN

         rpa = MIN( rpaxx, rpaxy )

!        "func1" = INTEGRAL OF THE PROBABILITY (SLOW-SLOW) FOR THE TIME STEP
         CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, 0.D0)
         sumss = func1

!        TIME STEP 

         deltat = savg / ( sumss + sumsf )
!         tout = told + deltat

!         dtau = MIN( dtau, deltat )
!         IF( t.LT.tout ) RETURN
!         told = t

!         dtau = deltat



!        CALL RSIDE1 IN ORDER TO GET THE GROUP VELOCITY (NORMALIZED TO THE SPEED OF LIGHT)
         CALL rside1(0, u, deru)
cSAP120511
c         vgmod = clight * sqrt( deru(1)**2 + deru(2)**2 + ( u(2) * deru(3) )**2 )
         vgmod = clight * sqrt( deru(1)**2 + deru(2)**2 + 
     &                          ( u(2) * deru(3) )**2 )
!!!!!         vgmod = sqrt( deru(1)**2 + deru(2)**2 + ( u(2) * deru(3) )**2 )

!        "ds" = THE SPACE STEP FOR ADVANCING THE RAY 
!        (FACTOR 100 BECAUSE "vgmod" AND "dtau" ARE IN CGS BUT "prmt(6)"
!        IS EXPRESSED IN METERS)

         ds = vgmod * deltat / 100.0 ![m]

!! CONDITIONS FOR SCATTERING 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         s_out = s_old + ds
         delta_s_total = s_total_new - s_total_old
         delta_s_total = MIN(delta_s_total, ds)
         IF (s_total_new .lt. s_out ) THEN
         RETURN
         ELSEIF (s_total_new .gt. s_out ) THEN
         s_old = s_total_new
         s_total_new = s_total_old + delta_s_total 
cSAP120517
         s_total_old= s_old
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!         us = us + ds

!!         ds = vgmod * dtau 


!!!         prmt(6) = MIN( prmt(6), ds )
!!!         us = us + ds
         print*,'ROOTMODE =2 & FAST > 0: dtau and ds', deltat, ds


!        PROBABILITY TO HAVE SLOW-SLOW SCATTERING
         probxx = deltat * sumss
!        PROBABILITY TO HAVE SLOW-FAST SCATTERING
         probxy = deltat * sumsf
!        TOTAL PROBABILITY
         probtot = probxx + probxy

!        "xnum" = RANDOM NUMBER
         CALL RANDOM_NUMBER( xnum )

!        IF "xnum" IS LARGER THAN THE TOTAL PROBABILITY ("probtot") THEN NO SCATTERING
!        IF "xnum" IS SMALLER THAN THE PROBABILITY TO HAVE SLOW-SLOW SCATTERING ("probxx")
!        THEN THE PERPENDICULAR REFRACTIVE INDEX AFTER THE SCATTERING IS EQUAL TO 
!        THE PERPENDICULAR REFRACTIVE INDEX OF THE SLOW MODE
!        IF "xnum" IS LARGER THAN THE PROBABILITY TO HAVE SLOW-SLOW SCATTERING ("probxx")
!        AND SMALLER THAN THE TOTAL PROBABILITY THEN THE PERPENDICULAR REFRACTIVE INDEX
!        AFTER THE SCATTERING IS EQUAL TO THE PERPENDICULAR REFRACTIVE INDEX OF THE FAST MODE

!         IF( xnum.GT.probtot ) THEN
!            ds = 0.d0
!         ENDIF   

         IF( xnum.GT.probtot ) RETURN
            
         skpm = skp

         IF( xnum.GT.probxx.AND.xnum.LE.probtot ) THEN
            skpm = sqrt( nper2m )
         ENDIF   



      
      ENDIF 
!     IF STATEMENT FOR THE THREE CASES  





!     GBETA SUBROUTINE EVALUATES ANGLE OF ROTATION, "beta", OF THE 
!     PERPENDICULAR COMPONENT OF THE WAVE VECTOR
      call gbeta ( u, skp, skpm, dsyrms, beta )


!     POSITION AND K-VECTOR COMPONENTS
      z = u(1)
      r = u(2)
      phi = u(3)
      wkz = u(4) * ( omlh / clight )
      wkr = u(5) * ( omlh / clight )
      wkph = u(6) / r * ( omlh / clight ) 



!     CALCULATE THE VALUES OF "qr,qth,qph" FOR A RANDOM SCATTERING ANGLE.
!     "qr,qth,qph" ARE THE COMPONENTS OF THE K-VECTOR AFTER THE SCATTERING 


      call rotate (u, qr, qz, qph, beta, skp, skpm, wkr, wkz, wkph)


!     NEW COMPONENT OF THE REFRACTIVE INDEX AFTER THE SCATTERING
      u(4) = ( clight / omlh ) * qz
      u(5) = ( clight / omlh ) * qr
      u(6) = r * ( clight / omlh ) * qph



!     "beta" = rotation angle
      xbeta = beta * 180.D0 / pi


      END SUBROUTINE lhscattering
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------







!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
      SUBROUTINE wsort (u, nper2m, nper2p, rootmode)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    THE wsort SUBROUTINE                                                                 !!
!!    IDENTIFIES N_PERP FOR SLOW AND FAST WAVE AND COMPARED WITH THE                       !! 
!!    ACTUAL VALUE OF N_PERP GIVEN BY RAY TRACING EQS.                                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      USE parameters_lhscat
      IMPLICIT NONE
CSAP120511
c      INCLUDE 'F90DIR/param.i'
c      INCLUDE 'F90DIR/one.i'      
      INCLUDE 'param.i'
      INCLUDE 'one.i'      

      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: u
cSAP120603
      DOUBLE PRECISION, INTENT(OUT):: nper2m, nper2p
      INTEGER, INTENT(OUT):: rootmode

      DOUBLE PRECISION :: z, r, phi
!      DOUBLE PRECISION :: bz, br, bphi, bmod
!      DOUBLE COMPLEX :: cnper2p, cnper2m
      DOUBLE PRECISION :: npar, nper, nper2, nref, bmag
      DOUBLE PRECISION, DIMENSION(2) :: df
      INTEGER, DIMENSION(1):: iden

!!      REAL, EXTERNAL:: b  

!     POSITION      
      z = u(1)  
      r = u(2)
      phi = u(3)

!      bmod=b(z,r,phi)
!     PARALLEL REFRACTIVE INDEX (FROM RAY TRACING SOLUTION) (MAGNETIC FIELD COMPONENTS COME FROM "one.i")     
      npar = ( u(4) * bz + u(5) * br + u(6) * bphi / u(2) ) / bmod 
!     REFRACTIVE INDEX (FROM RAY TRACING SOLUTION)
      nref = sqrt( u(4)**2 + u(5)**2 + u(6)**2 / u(2)**2 )
!     PERPENDICULAR REFRACTIVE INDEX (FROM RAY TRACING SOLUTION)
      nper = sqrt( nref**2 - npar**2 )
      nper2 = nper * nper
      

!     "npernpar" EVALUATES THE TWO ROOTS FROM THE DISPERSION RELATION
!     "nperp2p" = n_perp^2 with "+" sign: SLOW MODE
!     "nperp2m" = n_perp^2 with "-" sign: FAST MODE
      CALL npernpar(z,r,phi,npar,nper2p,nper2m)


      df(1) = ABS( nper2m - nper2 )
      df(2) = ABS( nper2p - nper2 )

!     "MINLOC" PROVIDES THE INDEX OF THE ARRAY "df" WHERE THE DIFFERENCE IS SMALLER
      iden = MINLOC(df)

!     "rootmode = 1" (FAST MODE) or "rootmode = 2" (SLOW MODE)
      rootmode = iden(1)

      END SUBROUTINE wsort
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------







!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
      SUBROUTINE gbeta (u, skp, skpm, dsyrms, beta)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                         !!
!!    THE SUBROUTINE GBETA EVALUATES THE BETA ANGLE BETWEEN k AND k' (AFTER SCATTERING)      !!                                                                                               
!!                                                                                         !!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      USE parameters_lhscat  
      IMPLICIT NONE
cSAP120511
c      INCLUDE 'F90DIR/param.i'
c      INCLUDE 'F90DIR/one.i'
      INCLUDE 'param.i'
      INCLUDE 'one.i'

      DOUBLE PRECISION, INTENT(IN):: skp, skpm, dsyrms
      DOUBLE PRECISION, INTENT(OUT):: beta
      DOUBLE PRECISION:: func1, func2, func3
      DOUBLE PRECISION:: angle, pbo, pbx, xnum, yran 
      DOUBLE PRECISION, DIMENSION(6):: u 

!     TO INCLUDE COMMENTS WITH LINK TO THE EQUATIONS OF PAUL'S THESIS AND PAPER

      angle = 0.0

      CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, angle)

      pbo = func3

!     "xdum" = RANDOM NUMBER,  0 =< xdum =< 1
      CALL RANDOM_NUMBER(xnum)

!     RANDOM ANGLE OF ROTATION
      angle = pi * ( 2. * xnum - 1.0 )

      CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, angle)

      pbx = func3

!     Y_r = RANDOM NUMBER. SEE PAUL'S THESIS P. 76
      CALL RANDOM_NUMBER(xnum)
      yran = xnum

!!      IF (yran <= pbx/pbo) THEN
!!         beta = angle
!!      ELSE
      DO WHILE ( yran > pbx/pbo )
      CALL RANDOM_NUMBER(xnum) 
      angle = pi * ( 2. * xnum - 1.0 )
!
      CALL func123(u, skp, skpm, dsyrms, func1, func2, func3, angle)
      pbx = func3
      CALL RANDOM_NUMBER(xnum)
      yran= xnum
      END DO


      beta = angle 
!!      END IF 


      END SUBROUTINE gbeta
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------









!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
cSAP120511
c      SUBROUTINE rotate (u, qr, qz, qph, beta, skp, &
c              skpm, wkr, wkz, wkph) 
      SUBROUTINE rotate (u, qr, qz, qph, beta, skp,
     &         skpm, wkr, wkz, wkph) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    THE SUBROUTINE ROTATE EVALUATES THE COMPONENTS OF THE K' VECTOR (k_phi, k_theta AND k_r)    !!
!!    GIVEN K VECTOR, r, AND BETA ANGLE. IN PARTICULAR, SEE EQS. 8(a-c) OF BONOLI AND OTT,        !!
!!    Phys. Fluids 25 (1982) 359.                                                                 !!
!!    qph = k'_phi, qth = k'_theta and qr = k'_r                                                  !!
!!    wkpr = k_parallel, wkp = k_perp and wkpm = k'_perp.                                         !! 
!!    Note that prime indicates the k-components                                                  !!
!!    after the scattering event                                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      USE parameters_lhscat
      IMPLICIT NONE
cSAP120511
c      INCLUDE 'F90DIR/param.i'
c      INCLUDE 'F90DIR/one.i'
      INCLUDE 'param.i'
      INCLUDE 'one.i'

      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: u         
      DOUBLE PRECISION, INTENT(IN):: beta, skp, skpm, wkr, wkz, wkph
      DOUBLE PRECISION, INTENT(OUT):: qph, qz, qr
      DOUBLE PRECISION:: cosbta, sinbta, wkpr, wkp, wkpm, omlh
      DOUBLE PRECISION:: co1, co2, co3, co4
      DOUBLE PRECISION:: skpr  
!        DOUBLE PRECISION:: br, bz, bphi, bmod

!!        REAL, EXTERNAL:: b
       


!     OMEGA_LH = 2 * PI * FREQUENCY  ("frqncy" IS GIVEN IN "one.i" IN GHz)
      omlh = 2.D0 * pi * frqncy * 1.e+9


!!          bmod = b(u(1),u(2),u(3))
!     PARALLEL REFRACTIVE INDEX (FROM RAY TRACING SOLUTION) (MAGNETIC FIELD COMPONENTS COME FROM "one.i")     
      skpr = ( u(4) * bz + u(5) * br + u(6) * bphi / u(2) ) / bmod

!     COS and SIN OF THE ANGLE OF ROTATION "beta"
      cosbta = cos( beta )
      sinbta = sin( beta )


!     PARALLEL COMPONENT OF THE K-VECTOR 
      wkpr = skpr * ( omlh / clight)

!     PERPENDICULAR COMPONENT OF THE K-VECTOR 
      wkp = skp * ( omlh / clight)

!     PERPENDICULAR COMPONENT OF THE K-VECTOR AFTER THE SCATTERING
      wkpm = skpm * ( omlh / clight)

!     IT IS IMPORTANT TO NOTE THAT IN GENRAY WE USE THE CYLINDIRCAL COORDINATE R,phi,z INSTEAD OF
!     THE TOROIDAL COORDINATE r,theta,phi. AS A CONSEQUENCE WE HAVE TO SUBSITUTE CONSISTENTLY
!     THE COMPONENT OF K AND B. So "br" is the "b_R" of GENRAY, "bth" is "b_z" in GENRAY and
!     bph is b_phi in GENRAY. THE SAME THING HAPPENS FOR "qr, qth, qph".

!     co1 = k_phi - k_par * B_phi / |B|
      co1 = wkph - wkpr * bphi / bmod

!     co2 = k_z * B_r / |B| - k_r * B_z / |B| 
      co2 = ( wkz * br - wkr * bz ) / bmod

!     co3 = k_phi * B_r / |B| - k_r * B_phi / |B| 
      co3 = ( wkph * br - wkr * bphi ) / bmod

!     co4 = k_perp * k'_perp * cos(beta) + k_par^2
      co4 = wkp * wkpm * cosbta + wkpr**2


!     EQS. (8a), (8b) and (8c) OF PHYS. FLUIDS, 25 (1982) 366, RESPECTIVELY
!     "qph = k'_phi, qth = k'_theta and qr = k'_r"
!     NOTE THAT NOW "qth" IS SUBSTITUTED WITH "qz". THE SAME THING FOR "bth"
!     WITH "bz" FOR THE REASON EXPLAINED ABOVE
cSAP051211
c      qph = wkpr * bphi / bmod + co1 * cosbta * &
c            wkpm / wkp + co2 * sinbta * wkpm / wkp
      qph = wkpr * bphi / bmod + co1 * cosbta * 
     &       wkpm / wkp + co2 * sinbta * wkpm / wkp

      qz = ( br *co4 / bmod - wkr * wkpr - qph * co3 ) / co2
cSAP120511
c      qr = bmod * ( wkpr - qz * bz / bmod - qph * &
c           bphi / bmod ) / br
      qr = bmod * ( wkpr - qz * bz / bmod - qph * 
     &     bphi / bmod ) / br

      END SUBROUTINE rotate
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------








!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
cSAP120511
c        SUBROUTINE func123(u, skp, skpm, dsyrms, func1, func2, func3, &
c         angle)
        SUBROUTINE func123(u, skp, skpm, dsyrms, func1, func2, func3,
     &   angle)
!     &      angle_arg)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE FUNC123 EVALUATES THE SCATTERING EVENT PROBABILITY FOLLOWING THE BONOLI'S           !!
!!   MODEL                                                                                          !!
!!      BLA BLA BLA                                                                                 !!
!!                                                                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!        USE dielectric_tensor
        USE parameters_lhscat
        IMPLICIT NONE
cSAP120511
c        INCLUDE 'F90DIR/param.i'
c        INCLUDE 'F90DIR/one.i'
c        INCLUDE 'F90DIR/eps.i'
c        INCLUDE 'F90DIR/rho.i'
        INCLUDE 'param.i'
        INCLUDE 'one.i'
        INCLUDE 'eps.i'
        INCLUDE 'rho.i'

        DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: u        
        DOUBLE PRECISION, INTENT(IN):: skp, skpm, dsyrms
        DOUBLE PRECISION, INTENT(IN):: angle
!        DOUBLE PRECISION, INTENT(IN), OPTIONAL:: angle_arg
        DOUBLE PRECISION, INTENT(out) :: func1, func2, func3
        DOUBLE PRECISION :: skp2, skpm2, skpr2
        DOUBLE PRECISION :: wkp, wkpm, arg, cosb, sinb, bio, bi1
        DOUBLE PRECISION :: uxx, uxy, uzz, hxx, hxy, hzz, wpc
        DOUBLE PRECISION :: dtor, dtorm
        DOUBLE PRECISION :: aa, aam, bb, bbm
        DOUBLE PRECISION :: u1, u2, urr, uii
        DOUBLE PRECISION :: vj1, vj2, vj3, xj, yj, hdw, hdwm
        DOUBLE PRECISION :: theta, wkr, wkth, asp, wkph, beta, wkz
        DOUBLE PRECISION :: vgperp, vgperps, vgdotb
        DOUBLE PRECISION :: exx, exy, ezz
!        DOUBLE PRECISION :: bz, br, bphi, bmod, 
        DOUBLE PRECISION :: omlh, cosl, al, amin, kvac
         DOUBLE PRECISION :: qr, qz, qph
        DOUBLE PRECISION :: gauss12, gauss3, const1, const2, const3
        DOUBLE PRECISION, DIMENSION(6) :: uu, deruu
        DOUBLE PRECISION, DIMENSION(3) :: vgroup, bf
        DOUBLE COMPLEX, DIMENSION(3,3):: resp
        DOUBLE COMPLEX, DIMENSION(3,3):: dwepsdw
        DOUBLE COMPLEX:: cflown
        DOUBLE PRECISION:: wdepsxxdw, wdepsxydw, wdepszzdw
        DOUBLE PRECISION:: z, r, phi, cnz, cnr, cm
        DOUBLE PRECISION:: skpr  
        INTEGER:: i

!!        REAL, EXTERNAL:: b




      omlh = 2.D0 * pi * frqncy * 1.e+9

! in order to find the radial coordinate: one way is sqrt( Toroidal flux / (pi * B0) )
! the easiest way is sqrt(area) where area = areatot. I use this last one for the first moment.
! torftot (Tesla*m^2) = Toroidal flux given in rho.i
! areatot (cm^2)= Total cross-sectional area of LCFS given in rho.i
! see Bob's email on April 4 2012.
      amin = sqrt( areatot ) * 100.d0

      kvac = omlh / clight  ![cm^{-1}]
      al = clight / (omlh * amin) 

      wkp=skp / al
      wkpm = skpm / al 
!!      arg = 2. * wkp * wkpm / zcorr**2
! there is the factor "al" = c/(omega * a) in order to obtain
! (i) k_perp from N_perp and (ii) to give dimension to "zcorr" being dimensionless
! being equal  z_0 * a

      arg = 2. * skp * skpm * kvac**2 / zcorr**2



!!          bmod = b(u(1),u(2),u(3))
! n_parallel from GENRAY
      skpr = ( u(4) * bz + u(5) * br + u(6) * bphi / u(2) ) / bmod

!      IF(PRESENT(angle_arg))THEN
!         angle = angle_arg
!      ELSE
!         angle = 0.D0
!      END IF

      cosb = cos( angle )
      sinb = sin( angle )

      call bessel_nicola(bio,bi1,arg)

      skpr2 = skpr**2
      skp2 = skp**2
      skpm2 = skpm**2

!
!...  compute the elements of the mobility tensor "u"
!

!
!  THE COMPONENTS OF THE DIELECTRIC TENSOR HAVE TO COME FROM GENRAY: 
!  
!

      exx = real(reps(1,1))
      exy = abs(aimag(reps(1,2)))
      ezz = real(reps(3,3))

      uxx = exx - 1.0
      uxy = exy
      uzz = ezz - 1.0
!
!...  compute elements of the "w" derivativew tensor "h"
!
      z = u(1) 
      r   = u(2) 
      phi = u(3) 
      cnz = u(4) 
      cnr = u(5) 
      cm  = u(6) 

! FLOWN is a GENRAY subroutine which evauated the wave energy density. I added an
! argument "dwepsdw" which corresponds to the derivative of the product between 
! omega and the dielectric tensor with respect to the omega

      call flown(z,r,phi,cnz,cnr,cm,cflown,dwepsdw)

      wdepsxxdw = ( real(dwepsdw(1,1)) - exx )  
      wdepsxydw = ( abs(aimag(dwepsdw(1,2))) - exy )  
      wdepszzdw = ( real(dwepsdw(3,3)) - ezz )  
 
      hxx = exx + wdepsxxdw / 2.d0
      hxy = exy + wdepsxydw / 2.d0
      hzz = ezz + wdepszzdw / 2.d0
     

!
!...  compute elements of the polarization vectors
!
      dtor = ( skpr2 - exx ) * ( skpr2 + skp2 - exx ) - exy**2
      dtorm = ( skpr2 - exx ) * ( skpr2 + skpm2 - exx ) - exy**2

      aa = skpr * skp * ( skpr2 + skp2 - exx ) / dtor
      aam = skpr * skpm * ( skpr2 + skpm2 - exx ) / dtorm
      bb = skpr * skp * exy / dtor
      bbm = skpr * skpm * exy / dtorm
!
!...  compute elements of the interaction kernel "v"
!
      u1 = aa * aam + bb * bbm
      u2 = aa * bbm + aam * bb
      urr = uxx * u1 + uxy * u2
      uii = uxy * u1 + uxx * u2
cS130511AP
c      vj1 = ( uzz**2 + urr**2 ) * bio + 2. * uzz * urr * bi1 + &
c             ( uii**2 - urr**2 ) * bi1 / arg

c      vj2 = ( uzz**2 + uii**2 ) * ( bio - bi1 ) + 2. * uzz * urr * &
c           ( bi1 * ( 1. + 1. / arg ) - bio ) + ( urr**2 - uii**2 ) * &
c           ( bio * ( 1. + 1. / arg ) - bi1 * ( 1. + 1. / arg + 2. / &
c             arg**2 ) )
      vj1 = ( uzz**2 + urr**2 ) * bio + 2. * uzz * urr * bi1 + 
     &        ( uii**2 - urr**2 ) * bi1 / arg

      vj2 = ( uzz**2 + uii**2 ) * ( bio - bi1 ) + 2. * uzz * urr * 
     &      ( bi1 * ( 1. + 1. / arg ) - bio ) + ( urr**2 - uii**2 ) * 
     &      ( bio * ( 1. + 1. / arg ) - bi1 * ( 1. + 1. / arg + 2. / 
     &        arg**2 ) )

      xj = uzz + cosb * urr
      yj = sinb * uii
      vj3 = xj**2 + yj**2


      hdw = hxx * ( aa**2 + bb**2 ) + 2. * aa * bb * hxy + hzz
      hdwm = hxx * ( aam**2 + bbm**2 ) + 2. * aam * bbm * hxy + hzz

!
!...  calculate the perpendicular group velocity.
!


!      rho = y(1)
!      theta = y(2)
!      wkr = y(4)
!      wkth = y(5) / rho
!      asp = 1. /eps + rho * cos( theta )
!      wkph = wkn/asp



      wkz = u(4) * ( omlh / clight )
      wkr = u(5) * ( omlh / clight )
      wkph = u(6) / r * ( omlh / clight ) 


      beta = 0.0
      CALL rotate (u,qr,qz,qph,beta,skp,skpm,wkr,wkz,wkph)

!
!  FROM GENRAY FOR THE GROUP VELOCITY: TO CHECK NORMALIZATION AND COORDINATE SYSTEM!!!!
!
!





       uu(1) = z    
       uu(2) = r    
       uu(3) = phi  
       uu(4) = ( clight / omlh ) * qz
       uu(5) = ( clight / omlh ) * qr
       uu(6) = r * ( clight / omlh ) * qph


! RSIDE1 is a GENRAY subroutine which evaluates the right hand side terms
! of the ray tracing equations. So "deru" corresponds to the derivative of
! the dispersion function with respect to the refractive index, i.e, "deru"
! is the group velocity. In GENRAY the group velocity is normalized to c!

       CALL rside1(0, uu, deruu)

! group velocity components (NORMALIZED TO SPEED OF LIGHT in GENRAY)      
      vgroup(1) = deruu(1)
      vgroup(2) = deruu(2)
      vgroup(3) = r * deruu(3)

! CHECK IF I CAN USE THE MODULE one_from_genray OR I NEED TO RECALCULATE THE B FIELD USING THE 
! ROUTINE "b.f" IN GENRAY
!!      bmod  = b(u(1),u(2),u(3))
      bf(1) = bz 
      bf(2) = br
      bf(3) = bphi


! projection of the group velocity on the magnetic field multiplited by bmod (= module of the magnetic field)
      vgdotb = 0.0d0
      DO i = 1,3
        vgdotb = vgdotb + vgroup(i) * bf(i) 
      ENDDO

! square of the perpendicular (to the magnetic field) group velocity
      vgperps = 0.0d0
      DO i = 1,3
        vgperps = vgperps + ( vgroup(i) - vgdotb *bf(i) / bmod**2 )**2
      ENDDO

! perpendicular (to the magnetic field) group velocity
      vgperp = clight * SQRT(vgperps)

    


!      gauss12 = exp( - ( wkp**2 + wkpm**2 ) / zcorr**2 )
cSAP120511
c      gauss12 = exp( - ( skp**2 * kvac**2 + skpm**2 * kvac**2 ) / zcorr**2 )
      gauss12 = exp( - ( skp**2 * kvac**2 + skpm**2 * kvac**2 ) / 
     &                 zcorr**2 )

!      cosl = wkp**2 + wkpm**2 - 2. * wkp * wkpm * cosb
cSAP120511
c      cosl = skp**2 * kvac**2 + skpm**2 * kvac**2 - 2. * skp * skpm * kvac**2 * cosb
      cosl = skp**2 * kvac**2 + skpm**2 * kvac**2 - 2. * skp * skpm *
     &                                                  kvac**2 * cosb
      gauss3 = exp( - cosl / zcorr**2 )  
    
!      const1 = 0.5 * wkpm * dsyrms * omlh**2 * gauss12 / &
!               ( hdw * hdwm * vgperp * zcorr**2 )
cSAP120511
c      const1 = 0.5 * skpm * kvac * dsyrms * omlh**2 * gauss12 / &
c               ( hdw * hdwm * vgperp * zcorr**2 )
      const1 = 0.5 * skpm * kvac * dsyrms * omlh**2 * gauss12 / 
     7          ( hdw * hdwm * vgperp * zcorr**2 )

!      const2 = wkpm * dsyrms * omlh**2 * gauss12 / &
!               ( hdw * hdwm * vgperp * zcorr**2 )
cSAP120511
c      const2 = skpm * kvac * dsyrms * omlh**2 * gauss12 / &
c               ( hdw * hdwm * vgperp * zcorr**2 )
      const2 = skpm * kvac * dsyrms * omlh**2 * gauss12 / 
     &          ( hdw * hdwm * vgperp * zcorr**2 )

!      const3 = 0.5 * wkpm * dsyrms * omlh**2 * gauss3 / &
!               ( hdw * hdwm * vgperp * zcorr**2 )
cSAP120511
c      const3 = 0.5 * skpm * kvac * dsyrms * omlh**2 * gauss3 / &
c               ( hdw * hdwm * vgperp * zcorr**2 )
      const3 = 0.5 * skpm * kvac * dsyrms * omlh**2 * gauss3 / 
     &          ( hdw * hdwm * vgperp * zcorr**2 )

      func1 = 2. * pi * const1 * vj1
      func2 = vgperp / ( pi * const2 * vj2 )
      func3 = vj3 * const3


      END SUBROUTINE func123  
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------








!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
c      SUBROUTINE bessel_nicola1(bio, bi1, arg)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THE SUBROUTINE BESSEL EVALUATES THE MODIFIED BESSEL FUNCTIONS      !!
!! WHICH ARE NECESSARY FOR THE SCATTERING MODEL (SEE PAUL'S NOTES)    !!
!! IN PARTICULAR IN THE SUBROUTINE FUNC123 FOR THE EVALUATION OF THE  !!
!! PROBABILITY                                                        !!
!! bio = I_0(x)                                                       !!
!! bi1 = I_1(x)                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
c       IMPLICIT NONE

c       DOUBLE PRECISION, INTENT(IN):: arg
c       DOUBLE PRECISION, INTENT(OUT):: bio, bi1
c       DOUBLE PRECISION:: S18AEF, S18AFF
c       INTEGER:: IFAIL


c       EXTERNAL:: S18AEF, S18AFF


c       IFAIL = -1
c       bio = S18AEF( arg, IFAIL )
c       bi1 = S18AFF( arg, IFAIL )
  
c      END SUBROUTINE bessel_nicola1 
!------------------------------------------------------------------------------------------------------     
!------------------------------------------------------------------------------------------------------
      SUBROUTINE bessel_nicola(bio, bi1, arg)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THE SUBROUTINE BESSEL EVALUATES THE MODIFIED BESSEL FUNCTIONS      !!
!! WHICH ARE NECESSARY FOR THE SCATTERING MODEL (SEE PAUL'S NOTES)    !!
!! IN PARTICULAR IN THE SUBROUTINE FUNC123 FOR THE EVALUATION OF THE  !!
!! PROBABILITY                                                        !!
!! bio = I_0(x)                                                       !!
!! bi1 = I_1(x)                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
       IMPLICIT NONE

       DOUBLE PRECISION, INTENT(IN):: arg
       DOUBLE PRECISION, INTENT(OUT):: bio, bi1
c       DOUBLE PRECISION:: S18AEF, S18AFF
c       INTEGER:: IFAIL

       integer ize,ncalc
       real*8, dimension(1:2)  :: br,bi
c       EXTERNAL::  beslci1
c       IFAIL = -1
c       bio = S18AEF( arg, IFAIL )
c       bi1 = S18AFF( arg, IFAIL )

       ize=1 !for the modified bessel functions
       call beslci1(arg,0.d0,2,ize,br,bi,ncalc)
       bio=br(1)
       bi1=br(2)
      END SUBROUTINE bessel_nicola 


 
