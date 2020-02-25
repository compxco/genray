c        ********************** absorplh**********************
c        *                      -----                       *
c        *  this subroutine calculates the imaginary part   *
c        *  of refractive index using formula from:   	    *
c        *  P.Bonoli, Linear theory of lower hybrid heating *
c        *  IEEE Transaction on Plasma Science,Vol. PS-12,  *
c        *  No.2,(1984), p 95-107                           *
c        *  See also Paul T. Bonoli and Ronald C. Englade,  *
c        *  Phys. Fluids 29, 2937 (1986).                   *
c        ****************************************************
c         input parameters: u(6)  -z,r,phi,n_z,n_r,m  
c                           cnpar -N_par			  
c                           cnper -N_pre (real part)   
c                           tempe -electron temperature  (keV)  
c                           dense -electron density   (10*13/cm**3
c                           tempiar -ions temperature  (keV)      *
c                           b_z,b_r,b_phi,bmod-magnetic field	  *
c                           nbulk - number of plasma components	  *
c                           frqncy- wave friquency [GHz]	  *
c                           zeff  - plasma charge      		  *
c         output parameter: cnprim_e -imaginary part of N_perp    *
c                             (uncollisional electron absorption) *
c                           cnprim_i -imaginary part of N_perp    *
c                                (uncollisional ion absorption)   *
c                           cnprim_cl-imagionary part of N_perp	  *
c                                collisional(ions and electron)	  *
c-----------------------------------------------------------------*
      subroutine absorplh(u,cnpar,cnper,tempe,dense,tempiar
     1 ,b_z,b_r,b_phi,nbulk,bmod,frqncy,zeff,
     1 cnprim_e,cnprim_i,cnprim_cl)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      dimension u(6),deruu(6)
      dimension tempiar(nbulka)
********************************************************************
c electron absorption (Landau damping, formula (21a) in the Bonoli article)
c Di_e=2*sqrt(pi)*(omegape/omega)**2*(N_par*N_perp)**2*abs(x_oe**3)
c      *exp(-x_oe**2)
c--------------------------------------------------------------------
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
      cnpar2=cnpar*cnpar
      cnpar4=cnpar2*cnpar2
      cnper2=cnper*cnper
      cnper4=cnper2*cnper2
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c----------------------------------------
c     ve=(sqrt(Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(0.5d0*tempe)
      sqrt2=dsqrt(2.d0)
c----------------------------------------
c     x_oe=omega/((sqrt(2))*k_par*ve)=cvac/(sqrt(2)*N_par*ve)
      x_oe=cvac/(sqrt2*cnpar*ve)
      x_oe2=x_oe*x_oe
      x_oe3=x_oe2*dabs(x_oe)
      ye=y(z,r,phi,1)
      xe=x(z,r,phi,1)
c--------------------------------------
c     di_e is imaginary part of the dispersion function(uncollisional)
      di_e=2.d0*spi*xe*cnper2*cnpar2*x_oe3*dexp(-x_oe2)
c*********************************************************************
c     di_i is imaginary part of the dispersion function(uncollisional),
c     formula (21b) in the Bonoli article.
c----------------------------------------------------------
      di_i=0.d0

c      write(*,*)' absorplh.f nbulk',nbulk
      do i=2,nbulk
	tempi=tempiar(i)
c        write(*,*)'i,tempi,dmas(i)',i,tempi,dmas(i)
c-------vi=(sqrt(Ti/mi)),    vi in (cm/sec),Ti(keV)
        vi=1.87d9*dsqrt(0.5d0*tempi/dmas(i))
       
c-------x_oi=omega/((sqrt(2))*k_par*vi)=cvac/(sqrt(2)*N_par*vi)
        x_oi=cvac/(sqrt2*cnpar*vi)
        x_oi2=x_oi*x_oi
        x_oi3=x_oi2*dabs(x_oi)
        x_oi4=x_oi2*x_oi2
        yi=y(z,r,phi,i)
        xi=x(z,r,phi,i)

        sum=0.d0
        
c        write(*,*)'absorplh.f i,spi,xi,cnper4,vi',i,spi,xi,cnper4,vi
c        write(*,*)'absorplh.f x_oi,x_oi2,x_oi3,x_oi4',
c     &x_oi,x_oi2,x_oi3,x_oi4
         


        pdi_i=2.d0*spi*xi*cnper4*x_oi3*dexp(-x_oi2)
c        write(*,*)'pdi_i',pdi_i
        psum=yi*dabs(x_oi)/spi
	if (pdi_i.lt.1.d-10) goto 20
        n0=1/yi
        pn=10.d0*sqrt2*vi*dabs(cnpar)/cvac
        x_ni0=(1.d0-n0*yi)*cvac/(sqrt2*cnpar*vi)
        n1=(1-pn)/yi-1
        x_ni1=(1.d0-n1*yi)*cvac/(sqrt2*cnpar*vi)
        n2=(1+pn)/yi+1
        x_ni2=(1.d0-n2*yi)*cvac/(sqrt2*cnpar*vi)
        
c        write(*,*)'absorplh.f n1,n2',n1,n2 
        do n=n1,n2
c--------------------------------------
c         x_ni=(omega-n*omegas_ci)/sqrt(2)*k_par*vi)
          x_ni=(1.d0-n*yi)*cvac/(sqrt2*cnpar*vi)
          x_ni2=x_ni*x_ni
          if (x_ni2.gt.100)then
            goto 10
          endif
          sum=sum+dexp(-x_ni2)
 10       continue
        enddo !n
c        write(*,*)'absorplh.f,pdi_i,sum',pdi_i,sum
 20     sum=sum*psum

        di_i=di_i+pdi_i*sum

      enddo !i
c*********************************************************************
c     di_ic is the imaginary part of the dispersion function(collisional)
c     Collisional, ions and electron damping coefficients: (cgs units)
c     formula (26) from the Bonoli article
c     1.47e-9 = (4/3)*sqrt(2*pi)*16*q**4/(sqrt(me)*(1.6e-9 erg/keV)**1.5)
      rnuei = 1.47d-9*1.d+13*dense*zeff/tempe**1.5 ! [1/sec]
      omega=2.d0*pi*frqncy*1.d9	  ! [1/sec], frqncy in GHz
      di_ic=rnuei*(xe*cnper2/(ye*ye)+xe*cnpar2)*cnper2/omega
c--------------------------------------------------------------
c     calculation of (dD_real/dN_par)
c     D is from the P.Bonoli article  , formula (7),(9) and (17) 
      xi=x(z,r,phi,2)
      epsper=1.d0+xe/(ye*ye)-xi
      epspar=1.d0-xe-xi
      epsxy=xe/ye
c------------------------------------------
c p6 calculations
      ye2=ye*ye
      ye4=ye2*ye2
      pe=ve/cvac
      pe2=pe*pe
c     vi=(sqrt(Ti/mi)),    vi in (cm/sec),Ti(keV)
      tempi=tempiar(2)
      vi=1.87d9*dsqrt(0.5d0*tempi/dmas(2))
      xi=x(z,r,phi,2)
      pi=vi/cvac
      pi2=pi*pi
      p6=-(3.d0*xi*pi2+0.75d0*xe/ye4*pe2)
c------------------------------------------
      p4=epsper
      p2=(epsper+epspar)*(cnpar2-epsper)+epsxy*epsxy
      p0=epspar*((cnpar2-epsper)*(cnpar2-epsper)-epsxy*epsxy)
c------------------------------------------

      dddnper=(6.d0*p6*cnper4+4.d0*p4*cnper2+2.d0*p2)*cnper
      cnprim_e=-(di_e/dddnper)
      cnprim_i=-(di_i/dddnper)
      cnprim_cl=-(di_ic/dddnper) !   *coll_mult !BH: added 191017, Needs passing
      cnprim_e=dabs(cnprim_e)
      cnprim_i=dabs(cnprim_i)
      cnprim_cl=dabs(cnprim_cl)

c      write(*,*)'lh cnprim_e,cnprim_i,cnprim_cl',
c     1 cnprim_e,cnprim_i,cnprim_cl

      return
      end



      subroutine absorp_collisional(tempe,dense,frqncy,zeff,
     & v_gr_perp,coll_mult,
     & cnprim_cl)
c      calculates the imaginary part of refractive index
c      due to collisional electron-ion damping 
c      using the formula 
c      Im(N)=(nu_ei/(v_gr_perp))*clight/omega)
c
c      input parameters: 
c                        tempe -electron temperature  (keV)  
c                        dense -electron density   (10*13/cm**3	  
c                        frqncy- wave friquency [GHz]	  
c                        zeff  - plasma charge
c        		 v_gr_perp -perpendicular to the magnetic field
c                                   group velocity nirmalized by clight
c      output parameter: 
c                        cnprim_cl-imaginary part of N_perp	  
c                                collisional(ions and electron)	  
      implicit none

c-----input
      real*8 tempe,dense,frqncy,zeff,v_gr_perp,coll_mult

c-----output
      real*8 cnprim_cl

c-----local
      real*8 arg,cln,r_nu_ei,pi,omega
      
c-----Coulomb Logarithm = cln      
      arg=(1.d3/tempe)*dsqrt(10.d0*dense) !tempe KeV, 
                                           !dense 10**13/sm**3   
      cln=24.d0-dlog(arg)

c-----electron-ion collision rate =r_nu_ei [1/sec] 
      r_nu_ei=dense/(tempe*dsqrt(tempe))*zeff*cln*0.919d3 ![1/sec]

c      write(*,*)'dense,tempe,zeff,cln,r_nu_ei',
c     &dense,tempe,zeff,cln,r_nu_ei

c-----wave angle frequency 
      pi=4.d0*datan(1.d0)
      omega=2.d0*pi*frqncy*1.d9	  ! [1/sec], frqncy in GHz

      cnprim_cl=dabs(r_nu_ei/v_gr_perp/omega)*coll_mult
c      write(*,*)'cln,r_nu_ei,v_gr_perp,omega,cnprim_cl',
c     & cln,r_nu_ei,v_gr_perp,omega,cnprim_cl
      
      return
      end
