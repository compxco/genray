
      subroutine scatperp(bz,br,bphi,r,cnz,cnr,cm,scatd)
c
c  Scatters perp wavenumber by an angle given by a normal random
c  number generator with standard deviation dsqrt(2.*scatd)
c   Input: br,bz,bphi  magnetic field  components
c          cnr,cnz,cm  refractive index  components
c          r           the big radius
c          scatd       diffusion coefficient multiplied by the time step
c ************ ************ ************ ************ ************
c   Output:cnr,cnz,cm  refractive index components after
c                      n_perp scattering
c ************ ************ ************ ************ ************
      implicit none
      double precision bz,br,bphi,r,cnz,cnr,cm,scatd,
     +cnphi,scalbn,bmod2,bmod,cnpar_z,cnpar_r,cnpar_p,
     +cnper_z,cnper_r,cnper_p,cnper, 
     +elx_z,elx_r,elx_p,
     +ely_z,ely_r,ely_p,
     +elz_z,elz_r,elz_p,
     +deltheta,cs,sn,
     +cnper1lx,cnper1ly,cnper1lz,
     +cnper1_z,cnper1_r,cnper1_p
      real ranorm,fseed
c      write(*,*)'bz,br,bphi'
c      write(*,*)bz,br,bphi
c      write(*,*)'cnz,cnr,cm'
c      write(*,*)cnz,cnr,cm
      write(*,*)'scatperp ,scatd=',scatd
      cnphi=cm/r
c      write(*,*)'cnphi',cnphi
      scalbn=bz*cnz+br*cnr+bphi*cnphi
      bmod2=bz*bz+br*br+bphi*bphi
      bmod=dsqrt(bmod2)
c      write(*,*)'scalbn,bmod2,bmod'
c      write(*,*)scalbn,bmod2,bmod
ctest
c      write(*,*)'1 cnpar',scalbn/bmod
c-----------------------------------------------------------------------
c     Calculation of the components of N_parallel vector in (r,z,phi) frame
      cnpar_z=bz*scalbn/bmod2
      cnpar_r=br*scalbn/bmod2
      cnpar_p=bphi*scalbn/bmod2

c      write(*,*)'cnpar_z,cnpar_r,cnpar_p'
c      write(*,*)cnpar_z,cnpar_r,cnpar_p
ctest
c      write(*,*)'2 cnpar',dsqrt(cnpar_z**2+cnpar_r**2+cnpar_p**2)
c-----------------------------------------------------------------------
c     Calculation the components of N_perp vector (in r,z,phi frame)
      cnper_z=cnz-cnpar_z
      cnper_r=cnr-cnpar_r
      cnper_p=cnphi-cnpar_p
c      write(*,*)'cnper_z,cnper_r,cnper_p'
c      write(*,*)cnper_z,cnper_r,cnper_p
      cnper=dsqrt(cnper_z*cnper_z+cnper_r*cnper_r+cnper_p*cnper_p)
c      write(*,*)'cnper',cnper

      if (cnper.eq.0.d0) goto 10
c-----------------------------------------------------------------------
c     The local coordinate system (elx,ely,elz)
c     unit vector elz is parallel to the magnetic field
c     unit vector elx is parallel to the n_perp vector
c     unit vector ely=[elz*elx] is a vector multiplication
c-----------------------------------------------------------------------
c     The coordinates  of the local coordinate system in the 
c     frame (z,r,phi): el_x={ elx_z,elx_r,elx_p}
c     ely={ ely_z,ely_r,ely_p} ,el_z={ elz_z,elz_r,elz_p}   

      elx_z=cnper_z/cnper     
      elx_r=cnper_r/cnper     
      elx_p=cnper_p/cnper
     
      elz_z=bz/bmod     
      elz_r=br/bmod     
      elz_p=bphi/bmod

      ely_z=elz_r*elx_p-elz_p*elx_r     
      ely_r=-(elz_z*elx_p-elz_p*elx_z)     
      ely_p=elz_z*elx_r-elz_r*elx_z
     
c-----------------------------------------------------------------------
c     theta    is the polar angle around the axis parallel
c              to the magnetic field
c     deltheta is the deviation of the polar angle after the scattering

      deltheta=dsqrt(2.d0*scatd)*ranorm(fseed)
      write(*,*)'scatperp deltheta=',deltheta 
      cs=dcos(deltheta)
      sn=dsin(deltheta)

c     The scattered N_perp1 vector in {elx,ely,elz} frame

      cnper1lx=cnper*cs
      cnper1ly=cnper*sn
      cnper1lz=0.d0

c     The scattered N_perp1= vector in old {z,r,phi} frame

      cnper1_z=cnper1lx*elx_z+cnper1ly*ely_z
      cnper1_r=cnper1lx*elx_r+cnper1ly*ely_r
      cnper1_p=cnper1lx*elx_p+cnper1ly*ely_p
      
c     Recalculation of the refractive index components 
c     cnz,cnr,cm after the scattering
      cnz=cnpar_z+cnper1_z
      cnr=cnpar_r+cnper1_r
      cm=(cnpar_p+cnper1_p)*r

 10   continue
      return   
      end
 
      function ranorm(fseed)
c     generate one normal (0,1) pseudo random number using a method
c     based on central limit theorem and power residue method
c     output becomes truly normal as k (below) goes to infinity.
c     general formula for the deviate is
c        y=(sum of x(i),i=1,k) -k/2. )/ sqrt(k/12.)
c        where x(i) are are uniformally distributed on 0,1.
c     method is borrowed thru ibm. they use k=12.
c
      parameter (k=12,rtkd12=1.0000000)
      
      
c     rtkd12=sqrt(k/12.)
c
      a=0.
c *** note:  this loop is vector hazard.
c     option assert (hazard)
      do 1  i=1,k
c     y=ranf(fseed)

cSm 
      ix=17    
      y=rands_1(ix)

    1 a=a+y
      ranorm=(a-0.5*k)/rtkd12

      return
      end


      REAL FUNCTION rands_1(IXINPUT)

c     Purpose:
c         generation of pseudorandom number with the 
c         uniform distribution at the interval (0,1)
C     The integer input variable IX should be in : 0< IXINPUT < 32000
c     At the cyclic work rands function generates the sequence of
c     N<10**13  random numbers
c     Appl.Statistic, V31, p.188-190, 1982(?3)
c-----------------------------------------------------------------
c     Input:
c           The input integer variable 0< IXINPUT <32000 should be odd.
c           IX initializes rands1 function at the first its call. 
c-----------------------------------------------------------------

	INTEGER IX
cSm060725
c	INTEGER*2 IY,IZ
c	INTEGER*2 C2,C35,C63,C170,C171,C172,C176,C177,C178,
c    1		C269,C307,C323,CZERO
	INTEGER IY,IZ
	INTEGER C2,C35,C63,C170,C171,C172,C176,C177,C178,
     1		C269,C307,C323,CZERO

	REAL F30269,F30307,F30323,ONE
cBH060726        integer iysave,
        integer iysave
cIXINPUT
        save iysave,IX,IY,IZ
c	COMMON /RAND_/IY,IZ
	DATA C2/2/,C35/35/,C63/63/,C170/170/,C171/171/,C172/172/,
     1		C176/176/,C177/177/,C178/178/,C269/30269/,C307/30307/,
     2		C323/30323/,CZERO/0/,F30269/30269.0E0/,F30307/30307.0/,
     3		F30323/30323.0E0/,ONE/1.0E0/
         
c-----the first call of rands function
c     the variable IX will be set equal to the input variable IXINPUT
      rands_1 = 0
      if (iysave.ne.2) then 
        IX=IXINPUT
        IY=31
        IZ=89
      endif
c------------------------------------
      iysave=2

c        write(*,*)'1_IX,IY,IZ',IX,IY,IZ

	IX = C171*MOD(IX,C177) - C2*(IX/C177)
	IY = C172*MOD(IY,C176) - C35*(IY/C176)
	IZ = C170*MOD(IZ,C178) - C63*(IZ/C178)
	IF (IX .LT. CZERO) IX = IX + C269
	IF (IY .LT. CZERO) IY = IY + C307
	IF (IZ .LT. CZERO) IZ = IZ + C323
	rands = AMOD(FLOAT(IX)/F30269 + FLOAT(IY)/F30307 +
     1		FLOAT(IZ)/F30323,ONE)
cSm040708
        rands_1=rands

        IXINPUT=IX
c        write(*,*)'2_IX,IY,IZ',IX,IY,IZ
	RETURN
	END
c	BLOCK DATA ISEED2_
c	INTEGER*2 IY,IZ
c	COMMON /RAND_/IY,IZ
c	DATA IY,IZ/1,1/
c	END





