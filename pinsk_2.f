
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
      implicit double precision (a-h,o-z)
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

	 if(l2.gt.98) then
	   write(*,*)'in absorpfw_pinsker2 l2>98 it should increase
     1	   the parameter nb>100 in the subroutine absorpfd' 
	   stop
	 endif 

	 if(l1.gt.98) then
     	   write(*,*)'in absorpfw_pinsker2 l1>98 it should increase
     1	   the parameter nb>100 in the subroutine absorpfd' 
	   stop
	 endif 

	 if(l1.lt.-98) then
     	   write(*,*)'in absorpfw_pinsker2 l1<-98 it should increase
     1	   the parameter nb>100 in the subroutine absorpf1' 
	   stop
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
      if(i.gt.1) then
       reps(1,1)=reps(1,1) + ci*eps11_i + eps11_r
       reps(2,2)=reps(2,2) + ci*eps22_i + eps22_r
       reps(1,2)=reps(1,2) + eps12_r    + ci*eps12_i
      endif

      enddo !i  (ion species)

      cnprim_e=cnprim_s(1)


      write(*,*)'pinsker_2 cnprim_i',cnprim_i

      reps(2,1)=-reps(1,2)
      reps(3,1)=reps(1,3)
      reps(3,2)=-reps(2,3)

c***********************************************************************
      return
      end
