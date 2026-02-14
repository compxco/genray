

c     ********************** absorpfw*********************
c     *                      -----                       *
c     *  this subroutine calculates the imaginary part   *
c     *  of refractive index using a formula from:   	 *
c     *  S.C.Chiu,V.S.Chan,R.W.Harvey,M.Porkolab	 *
c     *  Theory of fast wave current drive for tokamak   *
c     *  plasma,					 *
c     *  Nuclear fusion,1989,vol.29,No.12, p.2175-2186   *
c     *                                                  *
c     *  Smirnov021018: Two corrections are introduced:  *
c     *    Eq. (6) for eps23 has "-" sign                *
c     *    Substitute Z_0(zeta)|Stix for all Z(zeta_e)   *
c     ****************************************************
c it uses the real*4 accuracy in the modified Bessel function
c which gives the zero ions damping (not good. See absorpf1, below).

c         input parameters: z,r,phi  (m,m,radian)
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)    
c                           tempe-electron temperature  (keV)     
c                           dense-electron density (10*13/cm**3	  
c                           tempiar(2:nbulka) ions temperature(kev)
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)
c                           
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c                           cnprim_i -imaginary part of N_perp    
c                                     (ion absorption)

c         The input parameter nb (the max order of the modified Bessel
c         functions) is set inside absorpfw using operator 'parameter'. 
c-----------------------------------------------------------------------

      subroutine absorpfw(z,r,phi,cnpar,cnper,
     1tempe,dense,tempiar,
     1nbulk,bmod,cnprim_e,cnprim_i)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      dimension tempiar(nbulka)
      double precision mu,mu2
      parameter (nb=100) ! the max order of the modifed Bessel function
      real br,bi
      dimension br(nb),bi(nb)
      complex argz,zf,zfp
      real rpr,rpim

c      write(*,*)'begin absorpfw '
c      write(*,*)'cnpar=',cnpar,'cnper=',cnper
c      write(*,*)'tempe=',tempe,'dense=',dense
c      write(*,*)'bmod=',bmod

c********************************************************************
c electron absorption:
c     Imag(N_perp)=abs(N_perp)*sqrt(pi)/4*beta_e*ksi_e*exp(-ksi_e**2)*G
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(2Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(tempe)	    !tempe keV
c----------------------------------------
c     cksi_e=omega/(k_par*ve)=cvac/(N_par*ve)
      cksi_e=cvac/(cnpar*ve)
c      cksi_e=dabs(cksi_e)
  
c----------------------------------------
c     beta_e=p_e/(btot**2/(8pi)),p_e electron pressure
c     p_e=dense*tempe
c     bmod                 ! tl
      btot=bmod*1.d4       ! gauss
      p_e=tempe*dense*1.6d+4	!(dyne/cm)
      p_b=btot**2/(8.d0*pi)	!magnetic pressure  (dyne/cm)
      beta_e=p_e/p_b
      mu=510.94d0/tempe
      mu2=mu*mu
      
c----------------------------------------
c     G=(N_par**2-S)/A2+
c       {1+omeg_pe**2*Z_r/(omega_e**2*abs(Z)**2*(N_par**2-S)}*
c       (m_e*c**2/T_e)**2*(omega_e*D/omega*abs(eps33))**2/(A1*A2)
c----------------------------------------
c     p_=P,d_=D,s_=S,  P,D,S are like in the article
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      s_=1.d0+xe/(ye*ye)
      d_=0.d0
      do i=2,nbulk
        xi=x(z,r,phi,i)
        yi=y(z,r,phi,i)
        p=1.d0/(1.d0-yi*yi)
        s_=s_-xi*p
c       d_=d_+xi*yi*p
        d_=d_+xi*p/yi
      enddo
      p_=-xe
      write(*,*)'xe,xi',xe,xi
      write(*,*)'ye,yi',ye,yi
      write(*,*)'s_,d_,p_',s_,d_,p_
c---------------------------------------
      rpr=real(cksi_e)
      argz=cmplx(rpr,0.)
      call zfunc(argz,zf,zfp)
      rpr=real(zf)
      zf_r=dble(rpr)      !real(Z(cksi_e))
      rpim=aimag(zf)      
      zf_im=dble(rpim)	  !imag(Z(cksi_e))
      cksi_e2=cksi_e*cksi_e
      z_r=2.d0*cksi_e2*(1.d0+cksi_e*zf_r) ! article formula (18)
      z_im=2.d0*cksi_e2*cksi_e*zf_im	  ! article formula (18)
      
      cmodez_2=z_r*z_r+z_im*z_im
      eps33_r=xe*z_r			  ! article formula (5)
      eps33_im=xe*z_im			  ! article formula (5)
      eps33_2=eps33_r*eps33_r+eps33_im*eps33_im
      cnpar2=cnpar*cnpar
      p1=cnpar2-s_			  ! for (15),(16),(17)
      dp1=1.d0/p1			  ! for	(15),(16),(17)
      z_rdz2=z_r/(cmodez_2)
      p2=xe*z_rdz2*dp1/(ye*ye)
      a1_=p1+p2*cnpar2			  ! (16)
      a2_=p1+p2*(cnpar2+s_)		  ! (17)
      g_=p1/a2_+(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
      cnprim_e=0.25d0*spi*cnper*beta_e*cksi_e*dexp(-cksi_e2)*g_	!(13)
       
      cn1=p1/a2_ !(15)
      cn2=(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
     
      cnpert2=(p1*p1-d_*d_)/(-p1+d_*d_*eps33_r*cnpar2/
     1(-eps33_2*p1))			  !(13)	
      cnpert=dsqrt(dabs(cnpert2))	  !(13)
	   
c end of the electron absorption calculations
c***********************************************************************
c ion absorption:
c-----------------------------------------------------------------
c     cksi_il=(1-l*yi)/(N_par*Vi/cvac)

c     eps11_i, eps22_i, eps12_r, - only with ion contributions
      eps11_i=0.d0
      eps22_i=0.d0
      eps12_r=0.d0

      do i=2,nbulk
c	 vi=(sqrt(2Te/me)),    vi in (cm/sec),Ti(keV)
	 vi=1.87d9*dsqrt(tempiar(i)/dmas(i))
	 xi=x(z,r,phi,i)
	 yi=y(z,r,phi,i)
         cmui=0.5d0*(cnper*vi/(yi*cvac))**2
	 expmui=dexp(-cmui)
         cksi_0=cvac/(cnpar*vi)
         cksi_0=dabs(cksi_0)

	 p22=xi*cksi_0*expmui	       ! for (2),(3),(4)
	 p11=p22/cmui		       ! for (2)
	 rpr=cmui
c	 nb=100
c------------------------------------------
c        izi=1 to calculate modified Bessel functions
c        BESSEL FUNCTIONS AND ARGUMENTS ARE REAL
	 ize=1
         call beslci(rpr,0.,nb,ize,br,bi,ncalc)
	 p=10.d0*dabs(cnpar)*vi/cvac
	 l0=1/yi
	 l1=(1-p)/yi
	 l2=(1+p)/yi+1
	 if(l2.gt.(nb-2)) then
	   WRITE(*,*)'in absorpwf l2>(nb-2)  Increase
     1	   parameter nb in the subroutine absorpwf' 
	   STOP
	 endif !l2=(nb-2)

	 if(l1.gt.(nb-2)) then
     	   WRITE(*,*)'in absorpwf l1>(nb-2)  Increase
     1	   parameter nb in the subroutine absorpwf' 
	   STOP
	 endif ! l1=(nb-2)

	 if(l1.lt.-(nb-2)) then
     	   WRITE(*,*)'in absorpwf l1<-(nb-2)  Increase
     1	   parameter nb in the subroutine absorpwf' 
	   STOP
	 endif !l1=-(nb-2)

	 do l=l1,l2
c          cksi_l=(1-l*yi)/(N_par*Vi/cvac)
	   if(l.ge.0) then
	     lmode=l
	   else
	    lmode=-l
	   endif
	   besm_l=br(lmode+1)	! |l| order modified Bessel function
	   besm_lp=br(lmode+2)
	   if (l.eq.0)then
	     cksi_l=cksi_0
             cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c	     besmp_l is a derivative from modified Bessel function I_0
c            by argument cmui :I_0 prime=I_1
	     besmp_l=besm_lp
	     eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
	   else
	     cksi_l=(1.d0-l*yi)*cvac/(cnpar*vi)
	     cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
	     eps11_i=eps11_i+p11*l*l*besm_l*zfi_l		   !(2)
c	     besmp_l is a derivative from modified Bessel function I_l
c            by argument cmui :I_l prime=0.5*(I_(l-1)+I_(l+1))
	     besm_lm=br(lmode)
	     besmp_l=0.5d0*(besm_lp+besm_lm)
	     eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
	     eps12_r=eps12_r+p22*l*(besm_l-besmp_l)*zfi_l	   !(4)
	   endif
         enddo !l
      enddo !i
      eps22_i=eps22_i+eps11_i					   !(3)
      cnprim_i=0.5d0*cnper*(
     1 (eps11_i-2.d0*d_*eps33_r*eps12_r/eps33_2)/a2_+
     1 (-p1*(eps22_i+eps11_i)-2.d0*d_*eps12_r)*a1_/((p1*p1-d_*d_)*a2_))
c      write(*,*)'after cnprim_i',cnprim_i

c end of the ion absorption calculations

c***********************************************************************
      return
      end

c     ********************** absorpf1**********************
c     *                      -----                       *
c     *  this subroutine calculates imaginary part       *
c     *  of refractive index using formula from:   	 *
c     *  S.C.Chiu,V.S.Chan,R.W.Harvey,M.Porkolab	 *
c     *  Theory of fast wave current drive for tokamak   *
c     *  plasma,					 *
c     *  Nuclear fusion,1989,vol.29,No.12, p.2175-2186   *
c     !!the differens with absorpwf in the modified bessel function 
c     !!the real*8 bessel function it gives the good ions absorption
c     *                                                  *
c     *  Smirnov021018: Two corrections are introduced:  *
c     *    Eq. (6) for eps23 has "-" sign                *
c     *    Substitute Z_0(zeta)|Stix for all Z(zeta_e)   *
c     ****************************************************
c         input parameters: z,r,phi  (m,m,radian) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c                           tempe- electron temperature  (keV)
c                           dense- electron density (10*13/cm**3
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)
c                           tempiar(2:nbulka) ions temperature(kev)
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)

c         The input parameter nb (the max number of the modified Bessel
c         functions) is set inside absorpfw using operator 'parameter'. 
      
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c                           cnprim_i -imaginary part of N_perp    
c                                     (ion absorption)  	  
c-----------------------------------------------------------------*
c    it uses beslci1 modified double precision Bessel functions   *
c------------------------------------------------------------------
      subroutine absorpf1(z,r,phi,cnpar,cnper,
     1tempe,dense,tempiar,
     1nbulk,bmod,cnprim_e,cnprim_i)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      dimension tempiar(nbulka)
      double precision mu,mu2
      parameter (nb=100) ! the max order of the modifed Bessel function
      dimension br(nb),bi(nb)
      complex argz,zf,zfp
      real rpr,rpim

c********************************************************************
c electron absorption:
c     Imag(N_perp)=abs(N_perp)*sqrt(pi)/4*beta_e*ksi_e*exp(-ksi_e**2)*G
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(2Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(tempe)	    !tempe keV
c----------------------------------------
c     cksi_e=omega/(k_par*ve)=cvac/(N_par*ve)
      cksi_e=cvac/(cnpar*ve)
      cksi_e=dabs(cksi_e)
c----------------------------------------
c     beta_e=p_e/(btot**2/(8pi)),p_e electron pressure
c     p_e=dense*tempe
c     bmod      ! tl
      btot=bmod*1.d4       ! gauss
      p_e=tempe*dense*1.6d+4	!(dyne/cm)
      p_b=btot**2/(8.d0*pi)	!magnetic pressure  (dyne/cm)
      beta_e=p_e/p_b

      mu=510.94d0/tempe
      mu2=mu*mu
c----------------------------------------
c     G=(N_par**2-S)/A2+
c       {1+omeg_pe**2*Z_r/(omega_e**2*abs(Z)**2*(N_par**2-S)}*
c       (m_e*c**2/T_e)**2*(omega_e*D/omega*abs(eps33))**2/(A1*A2)
c----------------------------------------
c     p_=P,d_=D,s_=S,  P,D,S are like in the article
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      s_=1.d0+xe/(ye*ye)
      d_=0.d0
      do i=2,nbulk
        xi=x(z,r,phi,i)
        yi=y(z,r,phi,i)
        p=1.d0/(1.d0-yi*yi)
        s_=s_-xi*p
        d_=d_+xi*p/yi
      enddo
      p_=-xe
c---------------------------------------
      rpr=real(cksi_e)
      argz=cmplx(rpr,0.)

      call zfunc(argz,zf,zfp)

      rpr=real(zf)
      zf_r=dble(rpr)      !real(Z(cksi_e))
      rpim=aimag(zf)      
      zf_im=dble(rpim)	  !imag(Z(cksi_e))
      cksi_e2=cksi_e*cksi_e
      z_r=2.d0*cksi_e2*(1.d0+cksi_e*zf_r) ! article formula (18)
      z_im=2.d0*cksi_e2*cksi_e*zf_im	  ! article formula (18)
      cmodez_2=z_r*z_r+z_im*z_im
      eps33_r=xe*z_r			  ! article formula (5)
      eps33_im=xe*z_im			  ! article formula (5)
      eps33_2=eps33_r*eps33_r+eps33_im*eps33_im
      cnpar2=cnpar*cnpar
      p1=cnpar2-s_			  ! for (15),(16),(17)
      dp1=1.d0/p1			  ! for	(15),(16),(17)
      z_rdz2=z_r/(cmodez_2)
      p2=xe*z_rdz2*dp1/(ye*ye)
      a1_=p1+p2*cnpar2			  ! (16)
      a2_=p1+p2*(cnpar2+s_)		  ! (17)
      g_=p1/a2_+(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
      cnprim_e=0.25d0*spi*cnper*beta_e*cksi_e*dexp(-cksi_e2)*g_	!(13)

      cn1=p1/a2_ !(15)
      cn2=(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)

c end of the electron absorption calculations
c***********************************************************************
c ion absorption:
c-----------------------------------------------------------------------
c     cksi_il=(1-l*yi)/(N_par*Vi/cvac)

c     eps11_i, eps22_i, eps12_r, - only with ion contributions
      eps11_i=0.d0
      eps22_i=0.d0
      eps12_r=0.d0
      do i=2,nbulk
c	 vi=(sqrt(2Te/me)),    vi in (cm/sec),Ti(keV)
	 vi=1.87d9*dsqrt(tempiar(i)/dmas(i))
	 xi=x(z,r,phi,i)
	 yi=y(z,r,phi,i)
         cmui=0.5d0*(cnper*vi/(yi*cvac))**2
	 expmui=dexp(-cmui)
         cksi_0=cvac/(cnpar*vi)
         cksi_0=dabs(cksi_0)

	 p22=xi*cksi_0*expmui	       ! for (2),(3),(4)
	 p11=p22/cmui		       ! for (2)
c	 nb=100
c------------------------------------------
c        izi=1 to calculate modified Bessel functions
c        BESSEL FUNCTIONS AND ARGUMENTS ARE REAL
	 ize=1
         call beslci1(cmui,0.d0,nb,ize,br,bi,ncalc)
	 p=10.d0*dabs(cnpar)*vi/cvac
	 l0=1/yi
	 l1=(1-p)/yi
	 l2=(1+p)/yi+1
	 if(l2.gt.(nb-2)) then
	   WRITE(*,*)'in absorpf1 l2>(nb-2)  Increase
     1	   parameter nb in the subroutine absorpf1' 
	   STOP
	 endif 

	 if(l1.gt.(nb-2)) then
     	   WRITE(*,*)'in absorpf1 l1>(nb-2)  Increase
     1	   parameter nb in the subroutine absorpf1' 
	   STOP
	 endif 

	 if(l1.lt.-(nb-2)) then
     	   WRITE(*,*)'in absorpf1 l1<-(nb-2)  Increase
     1	   parameter nb in the subroutine absorpf1' 
	   STOP
	 endif 

	 do l=l1,l2
c          cksi_l=(1-l*yi)/(N_par*Vi/cvac)
	   if(l.ge.0) then
	     lmode=l
	   else
	    lmode=-l
	   endif
	   besm_l=br(lmode+1)	! |l| order modified Bessel function
	   besm_lp=br(lmode+2)
	   if (l.eq.0)then
	     cksi_l=cksi_0
             cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c	     besmp_l is a derivative from modified Bessel function I_0
c            by argument cmui :I_0 prime=I_1
	     besmp_l=besm_lp
	     eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
	   else
	     cksi_l=(1.d0-l*yi)*cvac/(cnpar*vi)
	     cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
	     eps11_i=eps11_i+p11*l*l*besm_l*zfi_l		   !(2)
c	     besmp_l is a derivative from modified Bessel function I_l
c            by argument cmui :I_l prime=0.5*(I_(l-1)+I_(l+1))
	     besm_lm=br(lmode)
	     besmp_l=0.5d0*(besm_lp+besm_lm)
	     eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
	     eps12_r=eps12_r+p22*l*(besm_l-besmp_l)*zfi_l	   !(4)

	   endif
         enddo !l
      enddo !i

      eps22_i=eps22_i+eps11_i					   !(3)

      cnprim_i=0.5d0*cnper*(
     1 (eps11_i-2.d0*d_*eps33_r*eps12_r/eps33_2)/a2_+
     1 (-p1*(eps22_i+eps11_i)-2.d0*d_*eps12_r)*a1_/((p1*p1-d_*d_)*a2_))
c      write(*,*)'after cnprim_i',cnprim_i

c end of the ion absorption calculations

c***********************************************************************
      return
      end



c     ********************** absorpfd**********************
c     *                      -----                       *
c     *  this subroutine calculates imaginary part       *
c     *  of refractive index using formula from:   	 *
c     *  S.C.Chiu,V.S.Chan,R.W.Harvey,M.Porkolab	 *
c     *  Theory of fast wave current drive for tokamak   *
c     *  plasma,					 *
c     *  Nuclear fusion,1989,vol.29,No.12, p.2175-2186   *
c     !!the differens with absorpwf in the modified bessel function 
c     !!the real*8 bessel function it gives the good ions absorption
c     *                                                  *
c     *  Smirnov021018: Two corrections are introduced:  *
c     *    Eq. (6) for eps23 has "-" sign                *
c     *    Substitute Z_0(zeta)|Stix for all Z(zeta_e)   *
c     ****************************************************
c         input parameters: z,r,phi  (m,m,radian) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c                           tempe- electron temperature  (keV)
c                           dense- electron density (10*13/cm**3
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)
c                           tempiar(2:nbulka) ions temperature(kev)
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)

c         The input parameter nb (the max number of the modified Bessel
c         functions) is set inside absorpfw using operator 'parameter'. 
      
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c                           cnprim_i -imaginary part of N_perp    
c                                     (ion absorption) 
c
cBH041009:  Added breakdown of ion damping into components, cnprim_s(2:nbulk)
c 	  
c-----------------------------------------------------------------*
c    it uses beslci1 modified double precision Bessel functions   *
c------------------------------------------------------------------
      subroutine absorpfd(z,r,phi,cnpar,cnper,
     1tempe,dense,tempiar,
     1nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      include 'eps.i'

      dimension cnprim_s(*)   !BH041009, for separate ion contributions
                              !cnprim_s(1:nbulk), cnprim_s(1)=cnprim_e
      dimension tempiar(nbulka)
      double precision mu,mu2
cSAP091015
c      parameter (nb=100) ! the max order of the modifed Bessel function
cyup      parameter (nb=200) ! the max order of the modifed Bessel function
      parameter (nb=500) ! the max order of the modifed Bessel function
      dimension br(nb),bi(nb)

c      complex argz,zf,zfp
c      real rpr,rpim
      double complex argz,zf,zfp
      double precision rpr,rpim

      double complex ci,ckappae2,czero
         
      czero=dcmplx(0.d0,0.d0)
c********************************************************************
c electron absorption:
c     Imag(N_perp)=abs(N_perp)*sqrt(pi)/4*beta_e*ksi_e*exp(-ksi_e**2)*G
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(2Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(tempe)	    !tempe keV
c----------------------------------------
c     cksi_e=omega/(k_par*ve)=cvac/(N_par*ve)
      cksi_e=cvac/(cnpar*ve)
cSm021022
c      cksi_e=dabs(cksi_e)
c----------------------------------------
c     beta_e=p_e/(btot**2/(8pi)),p_e electron pressure
c     p_e=dense*tempe
c     bmod      ! tl
      btot=bmod*1.d4       ! gauss
      p_e=tempe*dense*1.6d+4	!(dyne/cm)
      p_b=btot**2/(8.d0*pi)	!magnetic pressure  (dyne/cm)
      beta_e=p_e/p_b

      mu=510.94d0/tempe
      mu2=mu*mu
c----------------------------------------
c     G=(N_par**2-S)/A2+
c       {1+omeg_pe**2*Z_r/(omega_e**2*abs(Z)**2*(N_par**2-S)}*
c       (m_e*c**2/T_e)**2*(omega_e*D/omega*abs(eps33))**2/(A1*A2)
c----------------------------------------
c     p_=P,d_=D,s_=S,  P,D,S are like in the article
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      s_=1.d0+xe/(ye*ye)
      d_=0.d0
      do i=2,nbulk
        xi=x(z,r,phi,i)
        write(*,*)'absorpfd i,xi',i,xi
        yi=y(z,r,phi,i)
        p=1.d0/(1.d0-yi*yi)
        s_=s_-xi*p
        d_=d_+xi*p/yi
      enddo
      p_=-xe
c---------------------------------------
c      rpr=real(cksi_e)
c      argz=cmplx(rpr,0.)
      argz=dcmplx(cksi_e,0.d0)


c      call zfunc(argz,zf,zfp)
cSm021022
c      call czeta(argz,zf,zfp,ierror)
      call CZETA0(cnpar,argz,zf,zfp,ierror)
c      rpr=real(zf)
c      rpr=dble(zf)
      rpr=dreal(zf)
      
      zf_r=dble(rpr)      !real(Z(cksi_e))
c      rpim=aimag(zf)      
      rpim=dimag(zf)      
      zf_im=dble(rpim)	  !imag(Z(cksi_e))

      cksi_e2=cksi_e*cksi_e
      z_r=2.d0*cksi_e2*(1.d0+cksi_e*zf_r) ! article formula (18)
      z_im=2.d0*cksi_e2*cksi_e*zf_im	  ! article formula (18)
      cmodez_2=z_r*z_r+z_im*z_im
      eps33_r=xe*z_r			  ! article formula (5)
      eps33_im=xe*z_im			  ! article formula (5)
      eps33_2=eps33_r*eps33_r+eps33_im*eps33_im
      cnpar2=cnpar*cnpar
      p1=cnpar2-s_			  ! for (15),(16),(17)
      dp1=1.d0/p1			  ! for	(15),(16),(17)
      z_rdz2=z_r/(cmodez_2)
      p2=xe*z_rdz2*dp1/(ye*ye)
      a1_=p1+p2*cnpar2			  ! (16)
      a2_=p1+p2*(cnpar2+s_)		  ! (17)
      g_=p1/a2_+(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
      cnprim_e=0.25d0*spi*cnper*beta_e*cksi_e*dexp(-cksi_e2)*g_	!(13)
cSm021021
      cnprim_e=dabs(cnprim_e)     

      cn1=p1/a2_ !(15)
      cn2=(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)

c end of the electron absorption calculations
      ci=dcmplx(0.d0,1.d0)
      cmue=0.5d0*(cnper*ve/(ye*cvac))**2
      ckappae2=xe*2.d0*cmue*cksi_e*(zf_r+ci*zf_im)  ! formula (7)

c       call tensrcld(z,r,phi)
   
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)

cSm021018
c-----tensor for the electrons with the hot correction
c     The hot correction is in eps(2,2),eps(3,3),eps(2,3)
c-----------------------------------------------------
      do i=1,3
         do j=1,3
         reps(i,j)=czero
         enddo
      enddo

      reps(1,1)=1.d0+xe/ye**2 !part of formula(2)

cSm021018 
c      reps(1,2)=-ci*xe/ye      !part of formula(4) 
cSm050207
c      reps(1,2)=xe/ye      !part of formula(4)
cSm080916
      reps(1,2)=-ci*xe/ye      !part of formula(4)


      reps(2,2)=reps(1,1)+ckappae2                  !part of formula (3)
      reps(3,3)=xe*(z_r+ci*z_im)                            !formula (5)
c      reps(2,3)=-ci*0.5*dsqrt(2.d0*cmue)*reps(3,3)/cksi_e  !formula (6)
cSm021018
      reps(2,3)=-ci*0.5*dsqrt(2.d0*cmue)*reps(3,3)/cksi_e
cSm021022
c      if(cnpar.lt.0.d0) reps(2,3)=-reps(2,3)
         
      cnprim_s(1)=cnprim_e

c      do i=2,nbulk
c        xi=x(z,r,phi,i)
c        yi=y(z,r,phi,i)
c        reps(1,1)=reps(1,1)+xi/(yi**2-1.d0)
c        reps(2,2)=reps(2,2)+xi/(yi**2-1.d0)
c        reps(1,2)=reps(1,2)-ci*yi*xi/(yi**2-1.d0)
c      enddo
c***********************************************************************
c ion absorption:
c-----------------------------------------------------------------
c     cksi_il=(1-l*yi)/(N_par*Vi/cvac)

c     eps11_i, eps22_i, eps12_r, - only with ion contributions
      cnprim_i=0.d0
      do i=2,nbulk

         eps11_i=0.d0
         eps22_i=0.d0
         eps12_r=0.d0

cSm050209 The elements eps11_r, eps22_r, eps12_i are used only
c        in the dielectric tensor reps()
c        They are not used in absorption (cnprim) calculations
         eps11_r=0.d0
         eps22_r=0.d0
         eps12_i=0.d0

c	 vi=(sqrt(2Te/me)),    vi in (cm/sec),Ti(keV)
	 vi=1.87d9*dsqrt(tempiar(i)/dmas(i))
	 xi=x(z,r,phi,i)
	 yi=y(z,r,phi,i)
         cmui=0.5d0*(cnper*vi/(yi*cvac))**2
cSAP091015
         write(*,*)'in absorpfd i,cnper,yi,vi/cvac,cmui',
     &                          i,cnper,yi,vi/cvac,cmui

	 expmui=dexp(-cmui)
         cksi_0=cvac/(cnpar*vi)
cSm050209
c         cksi_0=dabs(cksi_0)

	 p22=xi*cksi_0*expmui	       ! for (2),(3),(4)
	 p11=p22/cmui		       ! for (2)
c	 nb=100
c------------------------------------------
c        izi=1 to calculate modified Bessel functions
c        BESSEL FUNCTIONS AND ARGUMENTS ARE REAL
	 ize=1
         call beslci1(cmui,0.d0,nb,ize,br,bi,ncalc)

         write(*,*)'in absorpfd cmui,nb,ncalc',cmui,nb,ncalc
c         write(*,*)'br',br
c         write(*,*)'bi',bi
         do k=1,nb
            write(*,*)'k,expmui,br(k)*expmui,br(k)',
     &                 k,expmui,br(k)*expmui,br(k)
         enddo

	 p=10.d0*dabs(cnpar)*vi/cvac
	 l0=1/yi
	 l1=(1-p)/yi
	 l2=(1+p)/yi+1

         write(*,*)'l1,l2',l1,l2
cSAP091015
c	 if(l2.gt.98) then
	 if(l2.gt.(nb-2)) then
cSAP091015
           write(*,*)'l2,l1,cnpar,vi/cvac,yi,p',
     &                l2,l1,cnpar,vi/cvac,yi,p

	   WRITE(*,*)'in absorpfd l2>nb-2  Increase
     1	   parameter nb in the subroutine absorpfd' 
	   STOP
	 endif 

	 if(l1.gt.nb-2) then
     	   WRITE(*,*)'in absorpfd l1>nb-2  Increase
     1	   parameter nb in the subroutine absorpfd' 
	   STOP
	 endif 

	 if(l1.lt.-(nb-2)) then
     	   WRITE(*,*)'in absorpfd l1<-(nb-2)  Increase
     1	   parameter nb in the subroutine absorpfd' 
	   STOP
	 endif 

	 do l=l1,l2
c          cksi_l=(1-l*yi)/(N_par*Vi/cvac)
	   if(l.ge.0) then
	     lmode=l
	   else
	    lmode=-l
	   endif
	   besm_l=br(lmode+1)	! |l| order modified Bessel function
	   besm_lp=br(lmode+2)
	   if (l.eq.0)then
	     cksi_l=cksi_0
             cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c-------------------------------------------------------------------------
cSm050209    The calculatios for eps11_r,eps22_r,eps12_i
c            They are not used in absorption calculations
             argz=dcmplx(cksi_l,0.d0)             
             call CZETA0(cnpar,argz,zf,zfp,ierror)

             write(*,*)'absorpfd i,l,argz,zfi_l,imag(zf)',
     &                           i,l,argz,zfi_l,dimag(zf)

             zfr_l=dreal(zf)
c-------------------------------------------------------------------------
c	     besmp_l is a derivative from modified Bessel function I_0
c            by argument cmui :I_0 prime=I_1
	     besmp_l=besm_lp
	     eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)

cSm050209    The calculatios for eps22_r
c            They are not used in absorption calculations
             eps22_r=eps22_r+2.d0*p22*cmui*(besm_l-besmp_l)*zfr_l   !(3)
	   else
	     cksi_l=(1.d0-l*yi)*cvac/(cnpar*vi)
	     cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c-------------------------------------------------------------------------
cSm050209    The calculatios for eps11_r,eps22_r,eps12_i
c            They are not used in absorption calculations
             argz=dcmplx(cksi_l,0.d0)  
                        
             call CZETA0(cnpar,argz,zf,zfp,ierror)

             write(*,*)'absorpfd i,l,argz,zfi_l,zf',
     &                           i,l,argz,zfi_l,zf


             zfr_l=dreal(zf)
             zfi_l=dimag(zf)
c-------------------------------------------------------------------------
	     eps11_i=eps11_i+p11*l*l*besm_l*zfi_l		   !(2)
c	     besmp_l is a derivative from modified Bessel function I_l
c            by argument cmui :I_l prime=0.5*(I_(l-1)+I_(l+1))
	     besm_lm=br(lmode)
	     besmp_l=0.5d0*(besm_lp+besm_lm)
	     eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
	     eps12_r=eps12_r+p22*l*(besm_l-besmp_l)*zfi_l	   !(4)

cSm050209    The calculatios for epss11_r,eps22_r,eps12_i
c            They are not used in absorption calculations
             eps11_r=eps11_r+p11*l*l*besm_l*zfr_l                  !(2)
	     eps22_r=eps22_r+2.d0*p22*cmui*(besm_l-besmp_l)*zfr_l  !(3)
	     eps12_i=eps12_i-p22*l*(besm_l-besmp_l)*zfr_l	   !(4)

	   endif
         enddo !l  (bessel function order)

      eps22_i=eps22_i+eps11_i					   !(3)


      cnprim_s(i)=0.5d0*cnper*(
     1 (eps11_i-2.d0*d_*eps33_r*eps12_r/eps33_2)/a2_+
     1 (-p1*(eps22_i+eps11_i)-2.d0*d_*eps12_r)*a1_/((p1*p1-d_*d_)*a2_))

      cnprim_i=cnprim_i+cnprim_s(i)

c end of the ion absorption calculations

c-----add the hot ions terms to the tensor
cSm050209
c      reps(1,1)=reps(1,1)+eps11_i
c      reps(2,2)=reps(2,2)+eps22_i
c      reps(1,2)=reps(1,2)+ci*eps12_r

      reps(1,1)=reps(1,1) + ci*eps11_i + eps11_r
      reps(2,2)=reps(2,2) + ci*eps22_i + eps22_r
      reps(1,2)=reps(1,2) + eps12_r    + ci*eps12_i


      enddo !i  (ion species)


      write(*,*)'subroutine absorpfd after cnprim_i',cnprim_i

      reps(2,1)=-reps(1,2)
      reps(3,1)=reps(1,3)
      reps(3,2)=-reps(2,3)

c***********************************************************************
      return
      end




c     ********************** absorpfw_pinsker  ***********
c     *                      -----                       *
c     *  this subroutine calculates imaginary part       *
c     *  of refractive index using formula from:   	 *
c     *  R.i.Pinsker,M.Porkolab                 	 *
c     *  e-mail: pinsker@fusion.gat.com                  *
c     *   05/02/10
c     ****************************************************
c         input parameters: z,r,phi  (m,m,radian) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c                           tempe- electron temperature  (keV)
c                           dense- electron density (10*13/cm**3
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)
c                           tempiar(2:nbulka) ions temperature(kev)
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)

c         The input parameter nb (the max number of the modified Bessel
c         functions) is set inside absorpfw using operator 'parameter'. 
      
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c                           cnprim_i -imaginary part of N_perp    
c                                     (ion absorption) 
c
cBH041009:  Added breakdown of ion damping into components, cnprim_s(2:nbulk)
c 	  
c-----------------------------------------------------------------*
c    it uses beslci1 modified double precision Bessel functions   *
c------------------------------------------------------------------
      subroutine absorpfw_pinsker(z,r,phi,cnpar,cnper,
     1tempe,dense,tempiar,
     1nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      include 'eps.i'

      dimension cnprim_s(1)   !BH041009, for separate ion contributions
                              !cnprim_s(1:nbulk), cnprim_s(1)=cnprim_e
      dimension tempiar(nbulka)
      double precision mu,mu2
      parameter (nb=100) ! the max order of the modifed Bessel function
      dimension br(nb),bi(nb)

c      complex argz,zf,zfp
c      real rpr,rpim
       double complex argz,zf,zfp
       double precision rpr,rpim

      double complex ci,ckappae2,czero
         
      czero=dcmplx(0.d0,0.d0)
c********************************************************************
c electron absorption:
c     Imag(N_perp)=abs(N_perp)*sqrt(pi)/4*beta_e*ksi_e*exp(-ksi_e**2)*G
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(2Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(tempe)	    !tempe keV
c----------------------------------------
c     cksi_e=omega/(k_par*ve)=cvac/(N_par*ve)
      cksi_e=cvac/(cnpar*ve)
cSm021022
c      cksi_e=dabs(cksi_e)
c----------------------------------------
c     beta_e=p_e/(btot**2/(8pi)),p_e electron pressure
c     p_e=dense*tempe
c     bmod      ! tl
      btot=bmod*1.d4       ! gauss
      p_e=tempe*dense*1.6d+4	!(dyne/cm)
      p_b=btot**2/(8.d0*pi)	!magnetic pressure  (dyne/cm)
      beta_e=p_e/p_b

      mu=510.94d0/tempe
      mu2=mu*mu
c----------------------------------------
c     G=(N_par**2-S)/A2+
c       {1+omeg_pe**2*Z_r/(omega_e**2*abs(Z)**2*(N_par**2-S)}*
c       (m_e*c**2/T_e)**2*(omega_e*D/omega*abs(eps33))**2/(A1*A2)
c----------------------------------------
c     p_=P,d_=D,s_=S,  P,D,S are like in the article
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      s_=1.d0+xe/(ye*ye)
      d_=0.d0
      do i=2,nbulk
        xi=x(z,r,phi,i)
        yi=y(z,r,phi,i)
        p=1.d0/(1.d0-yi*yi)
        s_=s_-xi*p
        d_=d_+xi*p/yi
      enddo
      p_=-xe
c---------------------------------------
      argz=dcmplx(cksi_e,0.d0)

      call CZETA0(cnpar,argz,zf,zfp,ierror)

      rpr=dreal(zf)
      zf_r=dble(rpr)      !real(Z(cksi_e))     
      rpim=dimag(zf)      
      zf_im=dble(rpim)	  !imag(Z(cksi_e))

      cksi_e2=cksi_e*cksi_e
      z_r=2.d0*cksi_e2*(1.d0+cksi_e*zf_r) ! article formula (18)
      z_im=2.d0*cksi_e2*cksi_e*zf_im	  ! article formula (18)
      cmodez_2=z_r*z_r+z_im*z_im
      eps33_r=xe*z_r			  ! article formula (5)
      eps33_im=xe*z_im			  ! article formula (5)
      eps33_2=eps33_r*eps33_r+eps33_im*eps33_im
      cnpar2=cnpar*cnpar
      p1=cnpar2-s_			  ! for (15),(16),(17)
      dp1=1.d0/p1			  ! for	(15),(16),(17)
      z_rdz2=z_r/(cmodez_2)
      p2=xe*z_rdz2*dp1/(ye*ye)
      a1_=p1+p2*cnpar2			  ! (16)
      a2_=p1+p2*(cnpar2+s_)		  ! (17)
      g_=p1/a2_+(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
      cnprim_e=0.25d0*spi*cnper*beta_e*cksi_e*dexp(-cksi_e2)*g_	!(13)

      cnprim_e=dabs(cnprim_e)     

      cn1=p1/a2_ !(15)
      cn2=(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)

c end of the electron absorption calculations

      ci=dcmplx(0.d0,1.d0)
      cmue=0.5d0*(cnper*ve/(ye*cvac))**2
      ckappae2=xe*2.d0*cmue*cksi_e*(zf_r+ci*zf_im)  ! formula (7)

c       call tensrcld(z,r,phi)
   
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)

cSm021018
c-----tensor for the electrons with the hot correction
c     The hot correction is in eps(2,2),eps(3,3),eps(2,3)
c-----------------------------------------------------
      do i=1,3
         do j=1,3
         reps(i,j)=czero
         enddo
      enddo

      reps(1,1)=1.d0+xe/ye**2                       !part of formula(2)
      reps(1,2)=xe/ye                               !part of formula(4)
      reps(2,2)=reps(1,1)+ckappae2                  !part of formula (3)
      reps(3,3)=xe*(z_r+ci*z_im)                            !formula (5)
      reps(2,3)=-ci*0.5*dsqrt(2.d0*cmue)*reps(3,3)/cksi_e   !formula (6)
         
      cnprim_s(1)=cnprim_e

c      do i=2,nbulk
c        xi=x(z,r,phi,i)
c        yi=y(z,r,phi,i)
c        reps(1,1)=reps(1,1)+xi/(yi**2-1.d0)
c        reps(2,2)=reps(2,2)+xi/(yi**2-1.d0)
c        reps(1,2)=reps(1,2)-ci*yi*xi/(yi**2-1.d0)
c      enddo
c***********************************************************************
c ion absorption:
c-----------------------------------------------------------------
c     cksi_il=(1-l*yi)/(N_par*Vi/cvac)

c     eps11_i, eps22_i, eps12_r, - only with ion contributions

      cnprim_i=0.d0

      dL=1.d0-xe/(1.d0+dabs(ye))
      dR=1.d0-xe/(1.d0-dabs(ye))
      do i=2,nbulk
        xi=x(z,r,phi,i)
        yi=y(z,r,phi,i)
        dL=dL-xi/(1.d0-yi)
        dR=dR-xi/(1.d0+yi)
      enddo
      dS=0.5d0*(dL+dR)

      do i=2,nbulk
         
c	 vi=(sqrt(2Te/me)),    vi in (cm/sec),Ti(keV)
	 vi=1.87d9*dsqrt(tempiar(i)/dmas(i))
	 xi=x(z,r,phi,i)
	 yi=y(z,r,phi,i)
         cmui=0.5d0*(cnper*vi/(yi*cvac))**2
	 expmui=dexp(-cmui)
         cksi_0=cvac/(cnpar*vi)

	 p22=xi*cksi_0*expmui	        
	 p11=p22/cmui		      
c	 nb=100
c------------------------------------------
c        izi=1 to calculate modified Bessel functions
c        BESSEL FUNCTIONS AND ARGUMENTS ARE REAL
	 ize=1
         call beslci1(cmui,0.d0,nb,ize,br,bi,ncalc)
	 p=10.d0*dabs(cnpar)*vi/cvac
	 l0=1/yi
	 l1=(1-p)/yi
	 l2=(1+p)/yi+1

c         write(*,*)'absorpfw_pinsker l1,l2,1/yi',l1,l2,1./yi

	 if(l2.gt.(nb-2)) then
	   WRITE(*,*)'in absorpfw_pinsker l2>(nb-2)  Increase
     1	   parameter nb in the subroutine absorpfd' 
	   STOP
	 endif 

	 if(l1.gt.(nb-2)) then
     	   WRITE(*,*)'in absorpfw_pinsker l1>(nb-2)  Increase
     1	   parameter nb in the subroutine absorpfd' 
	   STOP
	 endif 

	 if(l1.lt.-(nb-2)) then
     	   WRITE(*,*)'in absorpfw_pinsker l1<-(nb-2)  Increase
     1	   parameter nb in the subroutine absorpf1' 
	   STOP
	 endif 

         cnprim_s(i)=0.d0
	 do l=l1,l2
c          cksi_l=(1-l*yi)/(N_par*Vi/cvac)
	   if(l.ge.0) then
	     lmode=l
	   else
	    lmode=-l
	   endif
	   besm_l=br(lmode+1)	! |l| order modified Bessel function
	   besm_lp=br(lmode+2)  ! |l|+1 order  modified Bessel function
	   besm_lm=br(lmode)    ! |l|-1 order  modified Bessel function

	   cksi_l=(1.d0-l*yi)*cvac/(cnpar*vi)
	   cksi_l2=cksi_l*cksi_l
           zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c-------------------------------------------------------------------------
cSm050209  The calculatios for eps11_r,eps22_r,eps12_i
c          They are not used in absorption calculations
           argz=dcmplx(cksi_l,0.d0)             
           call CZETA0(cnpar,argz,zf,zfp,ierror)
c           write(*,*)'absorpfw_pinskerc i,l,zfi_l,dimag(zf),dimag(zf)',
c     &                                  i,l,zfi_l,dimag(zf),dimag(zf)
           zfr_l=dreal(zf)
           zfi_l=dimag(zf)
c-------------------------------------------------------------------------
	   eps11_i=eps11_i+p11*l*l*besm_l*zfi_l		   !(2)
c	   besmp_l is a derivative from modified Bessel function I_l
c          by argument cmui :I_l prime=0.5*(I_(l-1)+I_(l+1))
	   besm_lm=br(lmode)
	   besmp_l=0.5d0*(besm_lp+besm_lm)
	   eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
	   eps12_r=eps12_r+p22*l*(besm_l-besmp_l)*zfi_l	   !(4)

cSm050209  The calculatios for epss11_r,eps22_r,eps12_i
c          They are not used in absorption calculations
           eps11_r=eps11_r+p11*l*l*besm_l*zfr_l                  !(2)
	   eps22_r=eps22_r+2.d0*p22*cmui*(besm_l-besmp_l)*zfr_l  !(3)
	   eps12_i=eps12_i-p22*l*(besm_l-besmp_l)*zfr_l	   !(4)


           delta_l_plus=besm_l-besm_lp
           delta_l_minus=besm_lm-besm_l

c       write(*,*)'pinsker i,l',i,l,cnpar
c       write(*,*)'delta_l_minus,(l-1)**2,((dR-cnpar2)/(dS-cnpar2))**2'
c   &             ,delta_l_minus,(l-1)**2,((dR-cnpar2)/(dS-cnpar2))**2
c       write(*,*)'delta_l_plus,(l+1)**2,((dL-cnpar2)/(dS-cnpar2))**2'
c   &             ,delta_l_plus,(l+1)**2,((dL-cnpar2)/(dS-cnpar2))**2
c       write(*,*)'delta_l_plus-delta_l_minus,cmui',
c   &              delta_l_plus-delta_l_minus,cmui
c      write(*,*)'0.25d0*spi*xi*dabs(cksi_0)*expmui/
c   &                 cnper)*dexp(-cksi_l2)'

           cnprim_s(i)=cnprim_s(i)+(0.25d0*spi*xi*dabs(cksi_0)*expmui/
     &                 dabs(cnper))*dexp(-cksi_l2)*(
     &                 l*((dR-cnpar2)/(dS-cnpar2))**2*delta_l_minus+
     &                 l*((dL-cnpar2)/(dS-cnpar2))**2*delta_l_plus+
     &                 2.d0*cmui*(delta_l_plus-delta_l_minus))
         enddo !l  (bessel function order)

      eps22_i=eps22_i+eps11_i					   !(3)
     
      cnprim_i=cnprim_i+cnprim_s(i)
c end of the ion absorption calculations

c-----add the hot ions terms to the tensor

      reps(1,1)=reps(1,1) + ci*eps11_i + eps11_r
      reps(2,2)=reps(2,2) + ci*eps22_i + eps22_r
      reps(1,2)=reps(1,2) + eps12_r    + ci*eps12_i


      enddo !i  (ion species)


c      write(*,*)'after cnprim_i',cnprim_i

      reps(2,1)=-reps(1,2)
      reps(3,1)=reps(1,3)
      reps(3,2)=-reps(2,3)

c***********************************************************************
      return
      end









c     ********************** absorpfw_pinsker_1 **********
c     *                      -----                       *
c     *  this subroutine calculates imaginary part       *
c     *  of refractive index for fast wave               *
c     *  using formula:   	                         *   
c     *  1) for ion absorption from                      *
c     *  R.i.Pinsker,M.Porkolab                          *
c     *  e-mail: pinsker@fusion.gat.com                  *
c     *   05/02/11                                       *  
c     *  2) for electron absorption from                 *
c     *  S.C.Chiu,V.S.Chan,R.W.Harvey,M.Porkolab	 *
c     *  Theory of fast wave current drive for tokamak   *
c     *  plasma,					 *
c     *  Nuclear fusion,1989,vol.29,No.12, p.2175-2186   *
c     ****************************************************
c         input parameters: z,r,phi  (m,m,radian) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c                           tempe- electron temperature  (keV)
c                           dense- electron density (10*13/cm**3
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)
c                           frqncy is a wave frequyncy GHZ
c                           tempiar(2:nbulka) ions temperature(kev)
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)

c         The input parameter nb (the max number of the modified Bessel
c         functions) is set inside absorpfw using operator 'parameter'. 
      
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c                           cnprim_i -imaginary part of N_perp    
c                                     (ion absorption) 
c
cBH041009:  Added breakdown of ion damping into components, cnprim_s(2:nbulk)
c 	  
c-----------------------------------------------------------------*
c    it uses beslci1 modified double precision Bessel functions   *
c------------------------------------------------------------------
      subroutine absorpfw_pinsker_1(z,r,phi,cnpar,cnper,
     1tempe,dense,tempiar,frqncy,
     1nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      include 'eps.i'
          
      dimension cnprim_s(1)   !BH041009, for separate ion contributions
                              !cnprim_s(1:nbulk), cnprim_s(1)=cnprim_e
      dimension tempiar(nbulka)
      double precision mu,mu2
      parameter (nb=100) ! the max order of the modifed Bessel function
      dimension br(nb),bi(nb)

c      complex argz,zf,zfp
c      real rpr,rpim
       double complex argz,zf,zfp
       double precision rpr,rpim

      double complex ci,ckappae2,czero
         
      czero=dcmplx(0.d0,0.d0)
c********************************************************************
c electron absorption:
c     Imag(N_perp)=abs(N_perp)*sqrt(pi)/4*beta_e*ksi_e*exp(-ksi_e**2)*G
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(2Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(tempe)	    !tempe keV
c----------------------------------------
c     cksi_e=omega/(k_par*ve)=cvac/(N_par*ve)
      cksi_e=cvac/(cnpar*ve)
cSm021022
c      cksi_e=dabs(cksi_e)
c----------------------------------------
c     beta_e=p_e/(btot**2/(8pi)),p_e electron pressure
c     p_e=dense*tempe
c     bmod      ! tl
      btot=bmod*1.d4       ! gauss
      p_e=tempe*dense*1.6d+4	!(dyne/cm)
      p_b=btot**2/(8.d0*pi)	!magnetic pressure  (dyne/cm)
      beta_e=p_e/p_b

      mu=510.94d0/tempe
      mu2=mu*mu
c----------------------------------------
c     G=(N_par**2-S)/A2+
c       {1+omeg_pe**2*Z_r/(omega_e**2*abs(Z)**2*(N_par**2-S)}*
c       (m_e*c**2/T_e)**2*(omega_e*D/omega*abs(eps33))**2/(A1*A2)
c----------------------------------------
c     p_=P,d_=D,s_=S,  P,D,S are like in the article
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      s_=1.d0+xe/(ye*ye)
      d_=0.d0
      do i=2,nbulk
        xi=x(z,r,phi,i)
        yi=y(z,r,phi,i)
        p=1.d0/(1.d0-yi*yi)
        s_=s_-xi*p
        d_=d_+xi*p/yi
      enddo
      p_=-xe
c---------------------------------------
      argz=dcmplx(cksi_e,0.d0)

      call CZETA0(cnpar,argz,zf,zfp,ierror)

      rpr=dreal(zf)
      zf_r=dble(rpr)      !real(Z(cksi_e))     
      rpim=dimag(zf)      
      zf_im=dble(rpim)	  !imag(Z(cksi_e))

      cksi_e2=cksi_e*cksi_e
      z_r=2.d0*cksi_e2*(1.d0+cksi_e*zf_r) ! article formula (18)
      z_im=2.d0*cksi_e2*cksi_e*zf_im	  ! article formula (18)
      cmodez_2=z_r*z_r+z_im*z_im
      eps33_r=xe*z_r			  ! article formula (5)
      eps33_im=xe*z_im			  ! article formula (5)
      eps33_2=eps33_r*eps33_r+eps33_im*eps33_im
      cnpar2=cnpar*cnpar
      p1=cnpar2-s_			  ! for (15),(16),(17)
      dp1=1.d0/p1			  ! for	(15),(16),(17)
      z_rdz2=z_r/(cmodez_2)
      p2=xe*z_rdz2*dp1/(ye*ye)
      a1_=p1+p2*cnpar2			  ! (16)
      a2_=p1+p2*(cnpar2+s_)		  ! (17)
      g_=p1/a2_+(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
      cnprim_e=0.25d0*spi*cnper*beta_e*cksi_e*dexp(-cksi_e2)*g_	!(13)

      cnprim_e=dabs(cnprim_e)     

      cn1=p1/a2_ !(15)
      cn2=(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)

c end of the electron absorption calculations

      ci=dcmplx(0.d0,1.d0)
      cmue=0.5d0*(cnper*ve/(ye*cvac))**2
      ckappae2=xe*2.d0*cmue*cksi_e*(zf_r+ci*zf_im)  ! formula (7)

c       call tensrcld(z,r,phi)
   
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)

cSm021018
c-----tensor for the electrons with the hot correction
c     The hot correction is in eps(2,2),eps(3,3),eps(2,3)
c-----------------------------------------------------
      do i=1,3
         do j=1,3
         reps(i,j)=czero
         enddo
      enddo

      reps(1,1)=1.d0+xe/ye**2                       !part of formula(2)
      reps(1,2)=xe/ye                               !part of formula(4)
      reps(2,2)=reps(1,1)+ckappae2                  !part of formula (3)
      reps(3,3)=xe*(z_r+ci*z_im)                            !formula (5)
      reps(2,3)=-ci*0.5*dsqrt(2.d0*cmue)*reps(3,3)/cksi_e   !formula (6)
         
      cnprim_s(1)=cnprim_e

c      do i=2,nbulk
c        xi=x(z,r,phi,i)
c        yi=y(z,r,phi,i)
c        reps(1,1)=reps(1,1)+xi/(yi**2-1.d0)
c        reps(2,2)=reps(2,2)+xi/(yi**2-1.d0)
c        reps(1,2)=reps(1,2)-ci*yi*xi/(yi**2-1.d0)
c      enddo
c***********************************************************************
c ion absorption:
c-----------------------------------------------------------------
c     cksi_il=(1-l*yi)/(N_par*Vi/cvac)

c     eps11_i, eps22_i, eps12_r, - only with ion contributions

      cnprim_i=0.d0

      dL=1.d0-xe/(1.d0+dabs(ye))
      dR=1.d0-xe/(1.d0-dabs(ye))
      do i=2,nbulk
        xi=x(z,r,phi,i)
        yi=y(z,r,phi,i)
        dL=dL-xi/(1.d0-yi)
        dR=dR-xi/(1.d0+yi)
      enddo
      dS=0.5d0*(dL+dR)

      do i=2,nbulk
         
c	 vi=(sqrt(2Te/me)),    vi in (cm/sec),Ti(keV)
	 vi=1.87d9*dsqrt(tempiar(i)/dmas(i))
	 xi=x(z,r,phi,i)
	 yi=y(z,r,phi,i)
         cmui=0.5d0*(cnper*vi/(yi*cvac))**2

         rho_larm=dsqrt(2.d0)*1.02d2*dsqrt(dmas(i)/1837.d0)/charge(i)*
     &            dsqrt(tempiar(i)*1.d+3)/btot   !cm
         c_kperp=cnper*2.d0*pi*frqncy*1.d9/cvac         !1/cm
         !write(*,*)'rho_larm,c_kperp',rho_larm,c_kperp
         clambda=0.5d0*(c_kperp*rho_larm)**2

         clambda_old=0.5d0*(cnper*vi/(yi*cvac))**2

         !write(*,*)'pinsker_1 clambda_old,clambda', clambda_old,clambda
     
         exp_lambda=dexp(-clambda)

         cksi_0=cvac/(cnpar*vi)

	 p22=xi*cksi_0*exp_lambda	        
	 p11=p22/clambda		      
c	 nb=100
c------------------------------------------
c        izi=1 to calculate modified Bessel functions
c        BESSEL FUNCTIONS AND ARGUMENTS ARE REAL
	 ize=1

         call beslci1(clambda,0.d0,nb,ize,br,bi,ncalc)

	 p=10.d0*dabs(cnpar)*vi/cvac
	 l0=1/yi
	 l1=(1-p)/yi
	 l2=(1+p)/yi+1

         !write(*,*)'absorpfw_pinsker l1,l2,1/yi',l1,l2,1./yi

	 if(l2.gt.(nb-2)) then
	   WRITE(*,*)'in absorpfw_pinsker_1 l2>(nb-2)  Increase
     1	   parameter nb' 
	   STOP
	 endif 

	 if(l1.gt.(nb-2)) then
     	   WRITE(*,*)'in absorpfw_pinsker_1 l1>(nb-2)  Increase
     1	   parameter nb' 
	   STOP
	 endif 

	 if(l1.lt.-(nb-2)) then
     	   WRITE(*,*)'in absorpfw_pinsker_1 l1<-(nb-2)  Increase
     1	   parameter nb' 
	   STOP
	 endif 

         cnprim_s(i)=0.d0
       eps11_i=0.d0 !YuP[2022-02-21] Was not initialized
       eps22_i=0.d0 !YuP[2022-02-21] Was not initialized
       eps12_i=0.d0 !YuP[2022-02-21] Was not initialized 
       eps11_r=0.d0 !YuP[2022-02-21] Was not initialized
       eps22_r=0.d0 !YuP[2022-02-21] Was not initialized
       eps12_r=0.d0 !YuP[2022-02-21] Was not initialized 
	 do l=l1,l2
c          cksi_l=(1-l*yi)/(N_par*Vi/cvac)
	   if(l.ge.0) then
	     lmode=l
	   else
	    lmode=-l
	   endif
	   besm_l=br(lmode+1)	! |l| order modified Bessel function
	   besm_lp=br(lmode+2)  ! |l|+1 order  modified Bessel function
	   besm_lm=br(lmode)    ! |l|-1 order  modified Bessel function

	   cksi_l=(1.d0-l*yi)*cvac/(cnpar*vi)
	   cksi_l2=cksi_l*cksi_l
           zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c-------------------------------------------------------------------------
cSm050209  The calculatios for eps11_r,eps22_r,eps12_i
c          They are not used in absorption calculations
           argz=dcmplx(cksi_l,0.d0)             
           call CZETA0(cnpar,argz,zf,zfp,ierror)
c           write(*,*)'absorpfw_pinskerc i,l,zfi_l,dimag(zf),dimag(zf)',
c     &                                  i,l,zfi_l,dimag(zf),dimag(zf)
           zfr_l=dreal(zf)
           zfi_l=dimag(zf)
c-------------------------------------------------------------------------
	   eps11_i=eps11_i+p11*l*l*besm_l*zfi_l		   !(2)
c	   besmp_l is a derivative from modified Bessel function I_l
c          by argument cmui :I_l prime=0.5*(I_(l-1)+I_(l+1))
	   besm_lm=br(lmode)
	   besmp_l=0.5d0*(besm_lp+besm_lm)
	   eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
	   eps12_r=eps12_r+p22*l*(besm_l-besmp_l)*zfi_l	   !(4)

cSm050209  The calculatios for epss11_r,eps22_r,eps12_i
c          They are not used in absorption calculations
           eps11_r=eps11_r+p11*l*l*besm_l*zfr_l                  !(2)
	   eps22_r=eps22_r+2.d0*p22*cmui*(besm_l-besmp_l)*zfr_l  !(3)
	   eps12_i=eps12_i-p22*l*(besm_l-besmp_l)*zfr_l	   !(4)


           delta_l_plus=besm_l-besm_lp
           delta_l_minus=besm_lm-besm_l

c       write(*,*)'pinsker i,l',i,l,cnpar
c       write(*,*)'delta_l_minus,(l-1)**2,((dR-cnpar2)/(dS-cnpar2))**2'
c   &             ,delta_l_minus,(l-1)**2,((dR-cnpar2)/(dS-cnpar2))**2
c       write(*,*)'delta_l_plus,(l+1)**2,((dL-cnpar2)/(dS-cnpar2))**2'
c   &             ,delta_l_plus,(l+1)**2,((dL-cnpar2)/(dS-cnpar2))**2
c       write(*,*)'delta_l_plus-delta_l_minus,cmui',
c   &              delta_l_plus-delta_l_minus,cmui
c      write(*,*)'0.25d0*spi*xi*dabs(cksi_0)*expmui/
c   &                 cnper)*dexp(-cksi_l2)'

           cnprim_s(i)=cnprim_s(i)+(0.25d0*spi*xi*cvac*exp_lambda/
     &                 (vi*cnper*dabs(cnpar)))*dexp(-cksi_l2)*(
     &                 l*((dR-cnpar2)/(dS-cnpar2))**2*delta_l_minus+
     &                 l*((dL-cnpar2)/(dS-cnpar2))**2*delta_l_plus+
     &                 2.d0*cmui*(delta_l_plus-delta_l_minus))
         enddo !l  (bessel function order)

      eps22_i=eps22_i+eps11_i					   !(3)
     
      cnprim_i=cnprim_i+cnprim_s(i)
c end of the ion absorption calculations

c-----add the hot ions terms to the tensor

      reps(1,1)=reps(1,1) + ci*eps11_i + eps11_r
      reps(2,2)=reps(2,2) + ci*eps22_i + eps22_r
      reps(1,2)=reps(1,2) + eps12_r    + ci*eps12_i


      enddo !i  (ion species)


c      write(*,*)'after cnprim_i',cnprim_i

      reps(2,1)=-reps(1,2)
      reps(3,1)=reps(1,3)
      reps(3,2)=-reps(2,3)

c***********************************************************************
      return
      end



c     ********************** absorp_hot ***** **********
c     *  this subroutine calculates imaginary part       *
c     *  of refractive index using formula from:   	 *
c     *  2k^_i_s*(P^+T^)=P_abs_s {Stix p.74 (17)}        *
c     *  and hot plasma polarization and flux            *
c     *   05/02/16                                       *
c     ****************************************************
c         input parameters: z,r,phi  (m,m,radian) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c                           nbulk is the number of the plasma species
c
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c                           cnprim_i -imaginary part of N_perp    
c                                     (ion absorption) 
c                           cnprim_s(nbulk) imaginary parts of N_perp_s
c                                           for all species
c
c                           It calculates the electric field polarization
c                           using the hot plasma tensor.The polarization
c                           (cex,cey,cez) will be in common /efield.i/ 
c------------------------------------------------------------------
      subroutine absorp_hot(z,r,phi,cnpar,cnper,nbulk,
     &cnprim_e,cnprim_i,cnprim_s)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      include 'eps.i'
      include 'cefield.i'
      
c-----input 
      integer nbulk
      real*8  z,r,phi,cnper,cnpar
c-----output
      real*8 cnprim_e,cnprim_i,cnprim_s(nbulka)
    
c-----local
      real*8 x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     & tpop_ar(nbulka),vflow_ar(nbulka),
     & cnprim,te,power_abs_s(nbulka),omega,bmod,pi,
     & T_perp_r,d_reps_herm_d_n_perp,
     & ex,ey,ez,eplus,eminus,epar,pp,step,cnper_p,cnper_m
     & Poynting_vector_perp_r
     
      complex*16 dhot, ! hot plasma dispersion function
     & K_s(3,3),     ! susceptibilties of the given specie  
     & cnx,          ! complex n_perp
     & ci,
     & kappa_s(3,3), !anti-Hermitian part of eps for one specie
     & T_perp_c,power_abs_s_c,
     & reps_p_c(3,3),reps_m_c(3,3),d_reps_herm_d_n_perp_c,
     & cbx,cby,cbz,ce(3),cec(3),cb(3),cbc(3),
     & Poynting_vector_x_c,
     & T_perp_r_old

      integer iherm_loc,i,j,j1,j2
 
      logical first
     
c-----external
      real*8 b,x,y,tempe,tpoprho,vflowrho
      complex*16 dhot_sum

      data first /.true./
      save first

      if (first) T_perp_r_old=0.d0
      first=.false.      

      pi=4.d0*datan(1.d0) !YuP[2018-10-13][2020-01-27] BUG was: pi=4.d0*dtan(1.d0)
c----------------------------------------------
c     calculate hot plasma tensor reps(). 
c     reps() will be in common /eps/ in eps.i file.
c---------------------------------------------
      bmod=b(z,r,phi) ! small radius rho calculations.
                      ! rho will be in common /one.i/

      do i=1,nbulk
         x_ar(i)=x(z,r,phi,i)
	 y_ar(i)=y(z,r,phi,i)
         if(i.eq.1) y_ar(1)=-y_ar(1)
	 te=tempe(z,r,phi,i)        ! kev
	 t_av_ar(i)=te*1000.d0      ! ev 
         tpop_ar(i)=tpop_zrp(z,r,phi,i) !YuP[2024-08-14] was tpoprho(rho,i)
         vflow_ar(i)=vflowrho(rho,i)
      enddo

      iherm_loc=2
      dhot=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &vflow_ar,cnpar,cnper,iherm_loc,reps)
c--------------------------------------------------------------
c     Calculate electric field polarization cex,cey,cez
c     using hot plasma dielectric tensor reps()
c     and real N_perp=cnper.
c     The polarization will be in common /cefield/ in cefield.i file.
c-------------------------------------------------------------------   
      cnx=dcmplx(cnper,0.d0)
      call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
      ce(1)=cex
      ce(2)=cey   
      ce(3)=cez
c-------------------------------------------------------------------
c     complex components of the wave magnetic field
c
c     ce - wave electric field
c     cec- comlex conjugate wave electric field
c     cb - wave magnetic field
c     cbc- comlex conjugate wave magnetic field
c-------------------------------------------------------------------
      cbx=-cnpar*cey
      cby=-cnper*cez+cnpar*cex
      cbz=cnper*cey
 
      cb(1)=cbx
      cb(2)=cby
      cb(3)=cbz
      
      do j=1,3
         cec(j)=dconjg(ce(j))
 	 cbc(j)=dconjg(cb(j))
      enddo      
c------------------------------------------------------------------
c     Calculate apsorped power for all species
c    
c     Power=omega/(8Pi)[E~(i) . eps_a_herm(i,j) 'E~(j)] 
c------------------------------------------------------------------
c     The code calculates absorped powes for all species
c     using anti-hermitian part of the susceptibiltis for
c     seperate species. 
c     The code uses normalized Power=Power/omega
c------------------------------------------------------------------
c     1. Calculate susceptibilities K_s() for all species.      
c     2. Calculate absorped power  power_abs_s(i) for
c        all "i=1,...,nbulk" species.      
c-----------------------------------------------------------------
      ci=dcmplx(0.d0,1.d0)

      cnprim=0.d0
c      omega=2.d0*pi*frqncy*1.d9       !frqncy in GHZ

      do i=1,nbulk
       
         call DHOT_s_c(dmas(i),x_ar(i),y_ar(i),t_av_ar(i),tpop_ar(i),
     &   vflow_ar(i),cnpar,cnper,cnprim,iherm_loc,K_s)

         power_abs_s_c=0.d0

         do j1=1,3
           do j2=1,3
              kappa_s(j1,j2)=0.5d0*(K_s(j1,j2)-dconjg(K_s(j2,j1)))/ci
              power_abs_s_c=power_abs_s_c+cec(j1)*kappa_s(j1,j2)*ce(j2)
           enddo
         enddo

        power_abs_s_c=power_abs_s_c/(8.d0*pi)
 
        power_abs_s(i)=dreal(power_abs_s_c)
c         write(*,*)'absorpfw_hot i,power_abs_s_c,power_abs_s(i)',
c     &                           i,power_abs_s_c,power_abs_s(i)
      enddo !i
      write(*,*)'power_abs_s',power_abs_s
c------------------------------------------------------------------
c     calculations of the Poynting flux 
c
c                      c_light
c     Poynting_flux^=  ------- [ E^~ (vectr*) B^+ E^ (vectr*) B^~ ]
c                       16 Pi 
c-------------------------------------------------------------------
c     x direction of the Stix sytem is parallel to N_perp^.
c     So, x component of the Poynting flux is parallel to N_perp^.
c
c     Poynting_vector_x_c is complex x component of 
c     Poynting flux in the Stix coordinates. 
c------------------------------------------------------------
c     In the code it was used the normalized variable 
c     Poynting_vector_perp_r= Poynting /c_light
c-------------------------------------------------------------
      pp=1.d0/(16.d0*pi)
      Poynting_vector_x_c = pp*((cec(2)*cb(3)-cec(3)*cb(2))+
     &                          (ce(2)*cbc(3)-ce(3)*cbc(2)))     

      Poynting_vector_perp_r=dreal(Poynting_vector_x_c)
c-------------------------------------------------------------------
c     calculations T=T_vector the flux of nonelectromagnetic energy
c     Stix book p.74 (19)
c           omega        d                c_light       d
c     T^= - ----- E^~.----(eps_h^^).E^= - -------- E^~. ---(eps_h^^) .E^
c           16Pi         dk^              16Pi          dN^
c
c     For calculations of Im (K_perp) we need T^_perp.
c     T^perp it is parallel to Re( k^_perp)
c
c                c_light          d
c     T^_perp= - -------- E^~(i). -----------(eps_h(i,j) . E(j)
c                 16Pi            dRe(N_perp)              
c   
c-------------------------------------------------------------------
c     In the code it was used the normalized variable
c     T_perp_r = (T/c_light) 
c------------------------------------------------------------------
      step=1.d-7
c      step=1.d-3     

c      cnper_p=cnper+step
c      cnper_m=cnper-step

      cnper_p=cnper*(1.d0+step)
      cnper_m=cnper*(1.d0-step)


      dhot=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &vflow_ar,cnpar,cnper_p,1,reps_p_c)

      dhot=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &vflow_ar,cnpar,cnper_m,1,reps_m_c)

      T_perp_c=dcmplx(0.d0,0.d0)

      do j1=1,3
         do j2=1,3
c           d_reps_herm_d_n_perp=(dreal(reps_p_c(j1,j2))-
c     &                           dreal(reps_m_c(j1,j2)))/(2.d0*step)

c              d_reps_herm_d_n_perp=(dreal(reps_p_c(j1,j2))-
c     &                      dreal(reps_m_c(j1,j2)))/(2.d0*step*cnper)

            d_reps_herm_d_n_perp_c=(reps_p_c(j1,j2)-
     &                            reps_m_c(j1,j2))/(2.d0*step*cnper)
   
c           T_perp_c=T_perp_c+cec(j1)*d_reps_herm_d_n_perp*ce(j2)
            T_perp_c=T_perp_c+cec(j1)*d_reps_herm_d_n_perp_c*ce(j2)
 
         enddo
      enddo

      T_perp_r=-dreal(T_perp_c)/(16.d0*pi)
c      write(*,*)'T_perp_c,T_perp_r',T_perp_c,T_perp_r 

      if(T_perp_r.lt.0.d0) then
        write(*,*)'*********************************************'
        write(*,*)' WARNING: T_perp_r.lt.0.d0 in absorp_hot'
        write(*,*)' Using T_perp_r from the previous point'
        write(*,*)'*********************************************'
        T_perp_r=T_perp_old
      endif 
      T_perp_old=T_perp_r

c     test 
c      T_perp_r=0.d0
c     end test

      do i=1,nbulk
         cnprim_s(i)=0.5d0*power_abs_s(i)/
     &               (T_perp_r+Poynting_vector_perp_r)

         if(cnprim_s(i).lt.0.d0) then
            write(*,*)'!!!cnprim_s<0 i,cnprim_s(i)', i,cnprim_s(i) 
         endif

      enddo

     
      write(*,*)'T_perp_r,Poynting_vector_perp_r,T/P',
     &T_perp_r,Poynting_vector_perp_r,T_perp_r/Poynting_vector_perp_r


      write(*,*)'cnprim_s',cnprim_s
     

      cnprim_e=cnprim_s(1)

      cnprim_i=0.d0
      do i=2,nbulk
         cnprim_i=cnprim_i+cnprim_s(i)
      enddo

c--------------------------------------------------------
c     calculate comples hot plasma dielectric tensor reps
c     reps will be in common /eps.i/
c--------------------------------------------------------
      dhot=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &vflow_ar,cnpar,cnper,2,reps)


      return
      end


     







c     ********************** absorpfw_pinsker_2  **********
c     *                      -----                       *
c     *  this subroutine calculates imaginary part       *
c     *  of refractive index using formula from:   	 *
c     *  R.i.Pinsker,M.Porkolab                 	 *
c     *  e-mail: pinsker@fusion.gat.com                  *
c     *   05/02/10
c     ****************************************************
c         input parameters: z,r,phi  (m,m,radian) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c                           tempe- electron temperature  (keV)
c                           dense- electron density (10*13/cm**3
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)
c                           tempiar(2:nbulka) ions temperature(kev)
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)

c         The input parameter nb (the max number of the modified Bessel
c         functions) is set inside absorpfw using operator 'parameter'. 
      
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c                           cnprim_i -imaginary part of N_perp    
c                                     (ion absorption) 
c
cBH041009:  Added breakdown of ion damping into components, cnprim_s(2:nbulk)
c 	  
c     electron and ion terms are calculated by the same set of fomula
c-----------------------------------------------------------------*
c    it uses beslci1 modified double precision Bessel functions   *
c------------------------------------------------------------------
      subroutine absorpfw_pinsker_2(z,r,phi,cnpar,cnper,
     1tempe,dense,tempiar,
     1nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      include 'eps.i'

      dimension cnprim_s(1)   !BH041009, for separate ion contributions
                              !cnprim_s(1:nbulk), cnprim_s(1)=cnprim_e
      dimension tempiar(nbulka)
      double precision mu,mu2
      parameter (nb=100) ! the max order of the modifed Bessel function
      dimension br(nb),bi(nb)

      double complex argz,zf,zfp
      double precision rpr,rpim

      double complex ci,ckappae2,czero
      
      write(*,*)'pinsker2 z,r,phi,cnpar,cnper',z,r,phi,cnpar,cnper
      write(*,*)'pinsker2 tempe,dense,tempiar',tempe,dense,tempiar
      write(*,*)'pinsker2 nbulk,bmod', nbulk,bmod

      czero=dcmplx(0.d0,0.d0)
c********************************************************************
c electron absorption:
c     Imag(N_perp)=abs(N_perp)*sqrt(pi)/4*beta_e*ksi_e*exp(-ksi_e**2)*G
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(2Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(tempe)	    !tempe keV
c----------------------------------------
c     cksi_e=omega/(k_par*ve)=cvac/(N_par*ve)
      cksi_e=cvac/(cnpar*ve)
cSm021022
c      cksi_e=dabs(cksi_e)
c----------------------------------------
c     beta_e=p_e/(btot**2/(8pi)),p_e electron pressure
c     p_e=dense*tempe
c     bmod      ! tl
      btot=bmod*1.d4       ! gauss
      p_e=tempe*dense*1.6d+4	!(dyne/cm)
      p_b=btot**2/(8.d0*pi)	!magnetic pressure  (dyne/cm)
      beta_e=p_e/p_b

      mu=510.94d0/tempe
      mu2=mu*mu
c----------------------------------------
c     G=(N_par**2-S)/A2+
c       {1+omeg_pe**2*Z_r/(omega_e**2*abs(Z)**2*(N_par**2-S)}*
c       (m_e*c**2/T_e)**2*(omega_e*D/omega*abs(eps33))**2/(A1*A2)
c----------------------------------------
c     p_=P,d_=D,s_=S,  P,D,S are like in the article
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      s_=1.d0+xe/(ye*ye)
      d_=0.d0
      do i=2,nbulk
        xi=x(z,r,phi,i)
        yi=y(z,r,phi,i)
        p=1.d0/(1.d0-yi*yi)
        s_=s_-xi*p
        d_=d_+xi*p/yi
      enddo
      p_=-xe
c---------------------------------------
      argz=dcmplx(cksi_e,0.d0)

      call CZETA0(cnpar,argz,zf,zfp,ierror)

      rpr=dreal(zf)
      zf_r=dble(rpr)      !real(Z(cksi_e))     
      rpim=dimag(zf)      
      zf_im=dble(rpim)	  !imag(Z(cksi_e))

      cksi_e2=cksi_e*cksi_e
      z_r=2.d0*cksi_e2*(1.d0+cksi_e*zf_r) ! article formula (18)
      z_im=2.d0*cksi_e2*cksi_e*zf_im	  ! article formula (18)
      cmodez_2=z_r*z_r+z_im*z_im
      eps33_r=xe*z_r			  ! article formula (5)
      eps33_im=xe*z_im			  ! article formula (5)
      eps33_2=eps33_r*eps33_r+eps33_im*eps33_im
      cnpar2=cnpar*cnpar
      p1=cnpar2-s_			  ! for (15),(16),(17)
      dp1=1.d0/p1			  ! for	(15),(16),(17)
      z_rdz2=z_r/(cmodez_2)
      p2=xe*z_rdz2*dp1/(ye*ye)
      a1_=p1+p2*cnpar2			  ! (16)
      a2_=p1+p2*(cnpar2+s_)		  ! (17)
      g_=p1/a2_+(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
c      cnprim_e=0.25d0*spi*cnper*beta_e*cksi_e*dexp(-cksi_e2)*g_	!(13)

c      cnprim_e=dabs(cnprim_e)     

      cn1=p1/a2_ !(15)
      cn2=(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)

c end of the electron absorption calculations

      ci=dcmplx(0.d0,1.d0)
      cmue=0.5d0*(cnper*ve/(ye*cvac))**2
      ckappae2=xe*2.d0*cmue*cksi_e*(zf_r+ci*zf_im)  ! formula (7)

c       call tensrcld(z,r,phi)
   
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)

cSm021018
c-----tensor for the electrons with the hot correction
c     The hot correction is in eps(2,2),eps(3,3),eps(2,3)
c-----------------------------------------------------
      do i=1,3
         do j=1,3
         reps(i,j)=czero
         enddo
      enddo

      reps(1,1)=1.d0+xe/ye**2                       !part of formula(2)
      reps(1,2)=xe/ye                               !part of formula(4)
      reps(2,2)=reps(1,1)+ckappae2                  !part of formula (3)
      reps(3,3)=xe*(z_r+ci*z_im)                            !formula (5)
      reps(2,3)=-ci*0.5*dsqrt(2.d0*cmue)*reps(3,3)/cksi_e   !formula (6)
         
c      cnprim_s(1)=cnprim_e

c      do i=2,nbulk
c        xi=x(z,r,phi,i)
c        yi=y(z,r,phi,i)
c        reps(1,1)=reps(1,1)+xi/(yi**2-1.d0)
c        reps(2,2)=reps(2,2)+xi/(yi**2-1.d0)
c        reps(1,2)=reps(1,2)-ci*yi*xi/(yi**2-1.d0)
c      enddo
c***********************************************************************
c ion absorption:
c-----------------------------------------------------------------
c     cksi_il=(1-l*yi)/(N_par*Vi/cvac)

c     eps11_i, eps22_i, eps12_r, - only with ion contributions

      cnprim_i=0.d0

      dL=1.d0-xe/(1.d0+dabs(ye))
      dR=1.d0-xe/(1.d0-dabs(ye))
      do i=2,nbulk
        xi=x(z,r,phi,i)
        yi=y(z,r,phi,i)
        dL=dL-xi/(1.d0-yi)
        dR=dR-xi/(1.d0+yi)
      enddo
      dS=0.5d0*(dL+dR)

      write(*,*)'pinsker_2 nbulk',nbulk  
      do i=1,nbulk
         
c--------vi=(sqrt(2Te/me)),    vi in (cm/sec),Ti(keV)
         if (i.eq.1) then
            vi=1.87d9*dsqrt(tempe)	    !tempe keV
         else  
	    vi=1.87d9*dsqrt(tempiar(i)/dmas(i))
         endif
 
	 xi=x(z,r,phi,i)
	 yi=y(z,r,phi,i)
         if(i.eq.1)  yi=-yi !electron term
         
         cmui=0.5d0*(cnper*vi/(yi*cvac))**2
         write(*,*)'pinsker2 i,cmui', i,cmui
	 expmui=dexp(-cmui)
         cksi_0=cvac/(cnpar*vi)

	 p22=xi*cksi_0*expmui	        
	 p11=p22/cmui		      
c	 nb=100
c------------------------------------------
c        izi=1 to calculate modified Bessel functions
c        BESSEL FUNCTIONS AND ARGUMENTS ARE REAL
	 ize=1
         call beslci1(cmui,0.d0,nb,ize,br,bi,ncalc)
	 p=10.d0*dabs(cnpar)*vi/cvac
	 l0=1/yi
	 l1=(1-p)/yi
	 l2=(1+p)/yi+1

         write(*,*)'absorpfw_pinsker2 l1,l2,1/yi',l1,l2,1./yi

	 if(l2.gt.(nb-2)) then
	   WRITE(*,*)'in absorpfw_pinsker2 l2>(nb-2)  Increase
     1	   parameter nb in the subroutine absorpfw_pinsker2' 
	   STOP
	 endif 

	 if(l1.gt.(nb-2)) then
     	   WRITE(*,*)'in absorpfw_pinsker2 l1>(nb-2)  Increase
     1	   parameter nb in the subroutine absorpfw_pinsker2' 
	   STOP
	 endif 

	 if(l1.lt.-(nb-2)) then
     	   WRITE(*,*)'in absorpfw_pinsker2 l1<-(nb-2)  Increase
     1	   parameter nb in the subroutine absorpfw_pinsker2' 
	   STOP
	 endif 

         cnprim_s(i)=0.d0
       eps11_i=0.d0 !YuP[2022-02-21] Was not initialized
	 do l=l1,l2
c          cksi_l=(1-l*yi)/(N_par*Vi/cvac)
	   if(l.ge.0) then
	     lmode=l
	   else
	    lmode=-l
	   endif
	   besm_l=br(lmode+1)	! |l| order modified Bessel function
	   besm_lp=br(lmode+2)  ! |l|+1 order  modified Bessel function
	   besm_lm=br(lmode)    ! |l|-1 order  modified Bessel function

	   cksi_l=(1.d0-l*yi)*cvac/(cnpar*vi)
	   cksi_l2=cksi_l*cksi_l
           zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c-------------------------------------------------------------------------
cSm050209  The calculatios for eps11_r,eps22_r,eps12_i
c          They are not used in absorption calculations
           argz=dcmplx(cksi_l,0.d0)             
           call CZETA0(cnpar,argz,zf,zfp,ierror)
c           write(*,*)'absorpfw_pinsker2 i,l,zfi_l,dimag(zf),dimag(zf)',
c     &                                  i,l,zfi_l,dimag(zf),dimag(zf)
           zfr_l=dreal(zf)
           zfi_l=dimag(zf)
c-------------------------------------------------------------------------
	   eps11_i=eps11_i+p11*l*l*besm_l*zfi_l		   !(2)
c	   besmp_l is a derivative from modified Bessel function I_l
c          by argument cmui :I_l prime=0.5*(I_(l-1)+I_(l+1))
	   besm_lm=br(lmode)
	   besmp_l=0.5d0*(besm_lp+besm_lm)
	   eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
	   eps12_r=eps12_r+p22*l*(besm_l-besmp_l)*zfi_l	   !(4)

cSm050209  The calculatios for epss11_r,eps22_r,eps12_i
c          They are not used in absorption calculations
           eps11_r=eps11_r+p11*l*l*besm_l*zfr_l                  !(2)
	   eps22_r=eps22_r+2.d0*p22*cmui*(besm_l-besmp_l)*zfr_l  !(3)
	   eps12_i=eps12_i-p22*l*(besm_l-besmp_l)*zfr_l	   !(4)


           delta_l_plus=besm_l-besm_lp
           delta_l_minus=besm_lm-besm_l

c       write(*,*)'pinsker2 i,l',i,l,cnpar
c       write(*,*)'delta_l_minus,(l-1)**2,((dR-cnpar2)/(dS-cnpar2))**2'
c   &             ,delta_l_minus,(l-1)**2,((dR-cnpar2)/(dS-cnpar2))**2
c       write(*,*)'delta_l_plus,(l+1)**2,((dL-cnpar2)/(dS-cnpar2))**2'
c   &             ,delta_l_plus,(l+1)**2,((dL-cnpar2)/(dS-cnpar2))**2
c       write(*,*)'delta_l_plus-delta_l_minus,cmui',
c   &              delta_l_plus-delta_l_minus,cmui
c      write(*,*)'0.25d0*spi*xi*dabs(cksi_0)*expmui/
c   &                 cnper)*dexp(-cksi_l2)'

           cnprim_s(i)=cnprim_s(i)+(0.25d0*spi*xi*dabs(cksi_0)*expmui/
     &                 dabs(cnper))*dexp(-cksi_l2)*(
     &                 l*((dR-cnpar2)/(dS-cnpar2))**2*delta_l_minus+
     &                 l*((dL-cnpar2)/(dS-cnpar2))**2*delta_l_plus+
     &                 2.d0*cmui*(delta_l_plus-delta_l_minus))
           write(*,*)'pinsk_2 l,cnprim_s(i)', l,cnprim_s(i)
         enddo !l  (bessel function order)
         
      eps22_i=eps22_i+eps11_i					   !(3)
     
      if(i.gt.1) cnprim_i=cnprim_i+cnprim_s(i)
      write(*,*)'pinsk_2 i,cnprim_i,cnprim_s(i)',i,cnprim_i,cnprim_s(i)
c     end of the ion absorption calculations

c-----add the hot ions terms to the tensor

      reps(1,1)=reps(1,1) + ci*eps11_i + eps11_r
      reps(2,2)=reps(2,2) + ci*eps22_i + eps22_r
      reps(1,2)=reps(1,2) + eps12_r    + ci*eps12_i


      enddo !i  (ion species)

      cnprim_e=cnprim_s(1)


      write(*,*)'pinsker_2 cnprim_i',cnprim_i

      reps(2,1)=-reps(1,2)
      reps(3,1)=reps(1,3)
      reps(3,2)=-reps(2,3)

c***********************************************************************
      return
      end


!!!!!

c     ********************** absorpfd**********************
c     *                      -----                       *
c     *  this subroutine calculates imaginary part       *
c     *  of refractive index using formula from:   	 *
c     *  S.C.Chiu,V.S.Chan,R.W.Harvey,M.Porkolab	 *
c     *  Theory of fast wave current drive for tokamak   *
c     *  plasma,					 *
c     *  Nuclear fusion,1989,vol.29,No.12, p.2175-2186   *
c     !!the differens with absorpwf in the modified bessel function 
c     !!the real*8 bessel function it gives the good ions absorption
c     *                                                  *
c     *  Smirnov021018: Two corrections are introduced:  *
c     *    Eq. (6) for eps23 has "-" sign                *
c     *    Substitute Z_0(zeta)|Stix for all Z(zeta_e)   *
c     ****************************************************
c         input parameters: z,r,phi  (m,m,radian) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c                           tempe- electron temperature  (keV)
c                           dense- electron density (10*13/cm**3
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)
c                           tempiar(2:nbulka) ions temperature(kev)
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)

c         The input parameter nb (the max number of the modified Bessel
c         functions) is set inside absorpfw using operator 'parameter'. 
      
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c                           cnprim_i -imaginary part of N_perp    
c                                     (ion absorption) 
c
cBH041009:  Added breakdown of ion damping into components, cnprim_s(2:nbulk)
c 	  
c-----------------------------------------------------------------*
c    it uses beslci1 modified double precision Bessel functions   *
c
c     It uses    call ibess0 and ibessn instead of beslci1
c     to calculate I_n(x)*exp(-x)  
c     In is a modified double precision Bessel function
c------------------------------------------------------------------
      subroutine absorpfd_091016(z,r,phi,cnpar,cnper,
     1tempe,dense,tempiar,
     1nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      include 'eps.i'
      include 'writencdf.i' !to access freqcy, for printout

      dimension cnprim_s(*)   !BH041009, for separate ion contributions
                              !cnprim_s(1:nbulk), cnprim_s(1)=cnprim_e
      dimension tempiar(nbulka)
      double precision mu,mu2
cSAP091015
c      parameter (nb=100) ! the max order of the modified Bessel function
cyup      parameter (nb=200) ! the max order of the modified Bessel function
      parameter (nb=500) ! the max order of the modified Bessel function
      dimension br(nb),bi(nb)
      real*8 ri01        !modified Bessel zero oder I_0(x)*exp(-x)     
      real*8 xne(nb+1)   !array of modified Bessel oder n=(1,nb+1)  I_n(x)*exp(-x)  
c      complex argz,zf,zfp
c      real rpr,rpim
       double complex argz,zf,zfp
       double precision rpr,rpim

      double complex ci,ckappae2,czero
         
      czero=dcmplx(0.d0,0.d0)
c********************************************************************
c electron absorption:
c     Imag(N_perp)=abs(N_perp)*sqrt(pi)/4*beta_e*ksi_e*exp(-ksi_e**2)*G
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(2Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(tempe)	    !tempe keV
c----------------------------------------
c     cksi_e=omega/(k_par*ve)=cvac/(N_par*ve)
      cksi_e=cvac/(cnpar*ve)
cSm021022
c      cksi_e=dabs(cksi_e)
c----------------------------------------
c     beta_e=p_e/(btot**2/(8pi)),p_e electron pressure
c     p_e=dense*tempe
c     bmod      ! tl
      btot=bmod*1.d4       ! gauss
      p_e=tempe*dense*1.6d+4	!(dyne/cm)
      p_b=btot**2/(8.d0*pi)	!magnetic pressure  (dyne/cm)
      beta_e=p_e/p_b

      mu=510.94d0/tempe
      mu2=mu*mu
c----------------------------------------
c     G=(N_par**2-S)/A2+
c       {1+omeg_pe**2*Z_r/(omega_e**2*abs(Z)**2*(N_par**2-S)}*
c       (m_e*c**2/T_e)**2*(omega_e*D/omega*abs(eps33))**2/(A1*A2)
c----------------------------------------
c     p_=P,d_=D,s_=S,  P,D,S are like in the article
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      s_=1.d0+xe/(ye*ye)
      d_=0.d0
      do i=2,nbulk
        xi=x(z,r,phi,i)
c        write(*,*)'absorpfd_091016 i,xi',i,xi
        yi=y(z,r,phi,i)
        p=1.d0/(1.d0-yi*yi)
        s_=s_-xi*p
        d_=d_+xi*p/yi
      enddo
      p_=-xe
      cnpar2=cnpar*cnpar
      p1=cnpar2-s_			  ! for (15),(16),(17)
      dp1=1.d0/p1			  ! for	(15),(16),(17)
c---------------------------------------

      if(abs(cksi_e).gt. 1.d5)then !YuP[2022-08-17]correction
          !(Omega/|Kpar| is far above Vte)==> Set cnprim_e to zero.
          zf_r= -1.d0/cksi_e
          zf_im= 0.d0
          z_r= 0.d0
          z_im=0.d0
          eps33_r= 0.d0
          eps33_im=0.d0
          eps33_2= 0.d0
          p2= 0.d0
          a1_=p1 !+p2*cnpar2			  ! (16)
          a2_=p1 !+p2*(cnpar2+s_)		  ! (17)
          cnprim_e=0.d0
          !In fact, we could reduce this limit to abs(cksi_e).gt.10. -
          !There is no imaginary part in cpz0 output
          !of ZFUN_cur(X,cpz0,cpz1,cpz2), where X==cksi_e
          !See notes below.
      else ! Imaginary part can be present
          argz=dcmplx(cksi_e,0.d0)
          call CZETA0(cnpar,argz,zf,zfp,ierror)
          rpr=dreal(zf)
          zf_r=dble(rpr)      !real(Z(cksi_e))
          rpim=dimag(zf)      
          zf_im=dble(rpim)	  !imag(Z(cksi_e))
          cksi_e2=cksi_e*cksi_e
          z_r= 2.d0*cksi_e2*(1.d0+cksi_e*zf_r) ! article formula (18)
          !YuP[2022-08-17]: If |cksi_e|>>1 (Omega/|Kpar| is far above Vte)
          ! then  CZETA0() gives zf_r= -1/cksi_e and zf_im=0.
          ! In such a case,  (1.d0+cksi_e*zf_r) becomes zero,
          ! hence z_r=0, then cmodez_2=0 and eps33_2=0,
          ! which causes 1/0 below.
          ! In fact, we don't need to worry about case |cksi_e|>>1
          ! because there can be no damping (there is no imaginary part in zf
          ! output from CZETA0 [cpz0 output of ZFUN_cur(X,cpz0,cpz1,cpz2)]
          ! when |X|>10. where X==cksi_e). 
          ! Simply set cnprim_e to zero in this case.
          z_im=2.d0*cksi_e2*cksi_e*zf_im	  ! article formula (18)
          cmodez_2= z_r*z_r + z_im*z_im !YuP: Both z_r and z_im can be 0
          !at large |cksi_e|>>1.  Then, eps33_2 becomes zero.
          eps33_r=xe*z_r			  ! article formula (5)
          eps33_im=xe*z_im			  ! article formula (5)
          eps33_2=eps33_r*eps33_r+eps33_im*eps33_im
          z_rdz2=z_r/(cmodez_2)
          p2=xe*z_rdz2*dp1/(ye*ye)
          a1_=p1+p2*cnpar2			  ! (16)
          a2_=p1+p2*(cnpar2+s_)		  ! (17)
          g_=p1/a2_+(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
          cnprim_e=0.25d0*spi*cnper*beta_e*cksi_e*dexp(-cksi_e2)*g_	!(13)
          cnprim_e=dabs(cnprim_e)     
          cn1=p1/a2_ !(15) Not used?
          cn2=(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
      endif ! abs(cksi_e).gt. 1.d5 !YuP[2022-08-17]correction
c end of the electron absorption calculations

      !YuP: Next few lines - for tests only?
      ci=dcmplx(0.d0,1.d0)
      cmue=0.5d0*(cnper*ve/(ye*cvac))**2
      ckappae2=xe*2.d0*cmue*cksi_e*(zf_r+ci*zf_im)  ! formula (7)
c       call tensrcld(z,r,phi)
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
c-----tensor for the electrons with the hot correction
c     The hot correction is in eps(2,2),eps(3,3),eps(2,3)
c-----------------------------------------------------
      do i=1,3
         do j=1,3
         reps(i,j)=czero
         enddo
      enddo
      reps(1,1)=1.d0+xe/ye**2 !part of formula(2)
      reps(1,2)=-ci*xe/ye      !part of formula(4)
      reps(2,2)=reps(1,1)+ckappae2                  !part of formula (3)
      reps(3,3)=xe*(z_r+ci*z_im)                            !formula (5)
c      reps(2,3)=-ci*0.5*dsqrt(2.d0*cmue)*reps(3,3)/cksi_e  !formula (6)
      reps(2,3)=-ci*0.5*dsqrt(2.d0*cmue)*reps(3,3)/cksi_e
c      if(cnpar.lt.0.d0) reps(2,3)=-reps(2,3)
         
      cnprim_s(1)=cnprim_e

c      do i=2,nbulk
c        xi=x(z,r,phi,i)
c        yi=y(z,r,phi,i)
c        reps(1,1)=reps(1,1)+xi/(yi**2-1.d0)
c        reps(2,2)=reps(2,2)+xi/(yi**2-1.d0)
c        reps(1,2)=reps(1,2)-ci*yi*xi/(yi**2-1.d0)
c      enddo


c***********************************************************************
c ion absorption:
c-----------------------------------------------------------------
c     cksi_il=(1-l*yi)/(N_par*Vi/cvac)

c     eps11_i, eps22_i, eps12_r, - only with ion contributions
 
c     write(*,*)'absorpfd_091016 before ion absorption'

      cnprim_i=0.d0
      do i=2,nbulk

         eps11_i=0.d0
         eps22_i=0.d0
         eps12_r=0.d0

cSm050209 The elements eps11_r, eps22_r, eps12_i are used only
c        in the dielectric tensor reps()
c        They are not used in absorption (cnprim) calculations
         eps11_r=0.d0
         eps22_r=0.d0
         eps12_i=0.d0

c	 vi=(sqrt(2Te/me)),    vi in (cm/sec),Ti(keV)
	 vi=1.87d9*dsqrt(tempiar(i)/dmas(i))
	 xi=x(z,r,phi,i)
	 yi=y(z,r,phi,i)
         cmui=0.5d0*(cnper*vi/(yi*cvac))**2

cSAP091015
c         write(*,*)'in absorpfd i,cnper,yi,vi/cvac,cmui',
c     &                          i,cnper,yi,vi/cvac,cmui
     
         !--- Alternatively, for printout:
         if(i.eq.2)then !print for first ion species
         rho_larm=dsqrt(2.d0)*1.02d2*dsqrt(dmas(i)/1837.d0)/charge(i)*
     &            dsqrt(tempiar(i)*1.d+3)/btot   !cm
         c_kperp=cnper*6.28*freqcy/cvac         !1/cm
         clambda=0.5d0*(c_kperp*rho_larm)**2
!         if(clambda.gt.0.2)then
!          write(*,'(a,1p5e10.3)')
!     &    'R,Z[m],rho_L[cm],Kper[1/cm],0.5(Kper*rho_L)^2',
!     &      r, z, rho_larm, c_kperp, clambda
!         endif
         endif
         !--- done printout; Verified: clambda==cmui

	 expmui=dexp(-cmui)
         cksi_0=cvac/(cnpar*vi)
cSm050209
c         cksi_0=dabs(cksi_0)

cSAP091016
c	 p22=xi*cksi_0*expmui	       ! for (2),(3),(4)
	 p22=xi*cksi_0   	       ! for (2),(3),(4)

	 p11=p22/cmui		       ! for (2)

c	 nb=100
c------------------------------------------
c        izi=1 to calculate modified Bessel functions
c        BESSEL FUNCTIONS AND ARGUMENTS ARE REAL
	 ize=1

cSAP091016
c         call beslci1(cmui,0.d0,nb,ize,br,bi,ncalc)         
c         write(*,*)'in absorpfd cmui,nb,ncalc',cmui,nb,ncalc
c         write(*,*)'br',br
c         write(*,*)'bi',bi

         call ibess0(cmui,ri01)
         call ibessn(cmui,nb,ri01,xne)

c         write(*,*)'ri01',ri01
c         write(*,*)'xne',xne

         do k=1,nb
            if (k.eq.1) then
               br(1)=ri01
c               write(*,*)'k,expmui,br(1),br(1)/expmui',
c     &                    k,expmui,br(1),br(1)/expmui
            else 
               br(k)=xne(k-1)
c               write(*,*)'k,expmui,br(k),br(k)/expmui',
c     &               k,expmui,br(k),br(k)/expmui
            endif
         enddo
c----------------------------------------------------

	 p=10.d0*dabs(cnpar)*vi/cvac
	 l0=1/yi
	 l1=(1-p)/yi
	 l2=(1+p)/yi+1

c         write(*,*)'l1,l2',l1,l2

cSAP091015
c	 if(l2.gt.98) then
	 if(l2.gt.(nb-2)) then
cSAP091015
           write(*,*)'absorpfd_091016: l2,l1,cnpar,vi/cvac,yi,p',
     &                l2,l1,cnpar,vi/cvac,yi,p

	   WRITE(*,*)'absorpfd_091016: l2>nb-2   increase
     1	   the parameter nb>200 in the subroutine absorpfd_091016' 
	   STOP
	 endif 
cSAP091016
c	 if(l1.gt.98) then
	 if(l1.gt.nb-2) then
     	   WRITE(*,*)'absorpfd_091016: l1>nb-2   increase
     1	   the parameter nb>200 in the subroutine absorpfd_091016' 
	   STOP
	 endif 
cSAP091016
c	 if(l1.lt.-98) then
	 if(l1.lt.-(nb-2)) then
     	   WRITE(*,*)'absorpfd_091016: l1<-(nb02)   increase
     1	   the parameter nb>200 in the subroutine absorpfd_091016' 
	   STOP
	 endif 

	 do l=l1,l2
c          cksi_l=(1-l*yi)/(N_par*Vi/cvac)
	   if(l.ge.0) then
	     lmode=l
	   else
	    lmode=-l
	   endif
	   besm_l=br(lmode+1)	! |l| order modified Bessel function
	   besm_lp=br(lmode+2)

c           write(*,*)'besm_l,besm_lp',besm_l,besm_lp

	   if (l.eq.0)then
	     cksi_l=cksi_0
             cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c-------------------------------------------------------------------------
cSm050209    The calculatios for eps11_r,eps22_r,eps12_i
c            They are not used in absorption calculations
             argz=dcmplx(cksi_l,0.d0)             
             call CZETA0(cnpar,argz,zf,zfp,ierror)

c             write(*,*)'absorpfd i,l,argz,zfi_l,imag(zf)',
c     &                           i,l,argz,zfi_l,dimag(zf)

             zfr_l=dreal(zf)
c-------------------------------------------------------------------------
c	     besmp_l is a derivative from modified Bessel function I_0
c            by argument cmui :I_0 prime=I_1
	     besmp_l=besm_lp
	     eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)

cSm050209    The calculatios for eps22_r
c            They are not used in absorption calculations
             eps22_r=eps22_r+2.d0*p22*cmui*(besm_l-besmp_l)*zfr_l   !(3)
	   else ! |l|>0
	     cksi_l=(1.d0-l*yi)*cvac/(cnpar*vi)
	     cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c-------------------------------------------------------------------------
cSm050209    The calculatios for eps11_r,eps22_r,eps12_i
c            They are not used in absorption calculations
             argz=dcmplx(cksi_l,0.d0)  
                        
             call CZETA0(cnpar,argz,zf,zfp,ierror)

c             write(*,*)'absorpfd i,l,argz,zfi_l,zf',
c     &                           i,l,argz,zfi_l,zf


             zfr_l=dreal(zf)
             zfi_l=dimag(zf)
c-------------------------------------------------------------------------
	     eps11_i=eps11_i+p11*l*l*besm_l*zfi_l		   !(2)
c	     besmp_l is a derivative from modified Bessel function I_l
c            by argument cmui :I_l prime=0.5*(I_(l-1)+I_(l+1))
	     besm_lm=br(lmode)
	     besmp_l=0.5d0*(besm_lp+besm_lm)
	     eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
	     eps12_r=eps12_r+p22*l*(besm_l-besmp_l)*zfi_l	   !(4)

cSm050209    The calculatios for epss11_r,eps22_r,eps12_i
c            They are not used in absorption calculations
             eps11_r=eps11_r+p11*l*l*besm_l*zfr_l                  !(2)
	     eps22_r=eps22_r+2.d0*p22*cmui*(besm_l-besmp_l)*zfr_l  !(3)
	     eps12_i=eps12_i-p22*l*(besm_l-besmp_l)*zfr_l	   !(4)

	   endif
         enddo !l  (bessel function order)

      eps22_i=eps22_i+eps11_i					   !(3)
      
      eps_33_r2=0.d0 !and eps33_r=0 in this case.
      if(eps33_2.ne.0.d0)then !YuP[2022-08-17] added, to avoid 0/0
        eps_33_r2= eps33_r/eps33_2
      endif

!      write(*,*)'i, eps11_i-2.d0*d_*eps12_r*eps_33_r2, a2_=',
!     &  i, eps11_i-2.d0*d_*eps12_r*eps_33_r2, a2_

      cnprim_s(i)=0.5d0*cnper*(
     1 (eps11_i - 2.d0*d_*eps12_r*eps_33_r2)/a2_+
     1 (-p1*(eps22_i+eps11_i)-2.d0*d_*eps12_r)*a1_/((p1*p1-d_*d_)*a2_))

      cnprim_i=cnprim_i+cnprim_s(i)

c end of the ion absorption calculations

c-----add the hot ions terms to the tensor
cSm050209
c      reps(1,1)=reps(1,1)+eps11_i
c      reps(2,2)=reps(2,2)+eps22_i
c      reps(1,2)=reps(1,2)+ci*eps12_r

      reps(1,1)=reps(1,1) + ci*eps11_i + eps11_r
      reps(2,2)=reps(2,2) + ci*eps22_i + eps22_r
      reps(1,2)=reps(1,2) + eps12_r    + ci*eps12_i


      enddo !i  (ion species)


c      write(*,*)'subroutine absorpfd_091016 after cnprim_i',cnprim_i

      reps(2,1)=-reps(1,2)
      reps(3,1)=reps(1,3)
      reps(3,2)=-reps(2,3)

c***********************************************************************
      return
      end

