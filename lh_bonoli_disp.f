cMADE
c dddrz1.f 
c        call  ono_dervs(nbulk,u,wf,dddcnz,dddcnr,dddcm,
cMADE
c hamilt1.f
c        call ono_tens(dmas,X_ar,Y_ar,T_av_ar,nbulk,
cMADE
c rside1.f:
c           call  ono_dervs(nbulk,u,wf,dddcnz,dddcnr,dddcm,
c
c11/11/14
c     The dielectric tensor, its derivatives for
c     LH wave conditions using Paul T. Bonoli and Ronald C. Englade
c     approximaton: Simulation model for lower hybrid current drive,
c     Phys. Fluids 29(9), September 1986, pp. 2937-2950
c
c     The signs of the gyro frePaul T. Bonoli and Ronald C. Englade     
c     So, Y_e and Y_i has the same sign.    
      subroutine bonoli_Dispersion(
     & X_ar,Y_ar,T_ar, mass_ar,nbulk,
     & npar_in,nperp_in,D_0)
c     for LH wave approximation using Paul T. Bonoli and Ronald C. Englade
c     approximation: Simulation model for lower hybrid current drive,
c     Phys. Fluids 29(9), September 1986, pp. 2937-2950
c
c     gives D_0 Bonoli disperstion functiona
c------------------------------------------------------------
c     INPUTS:
c      nbulk is a number of the plasma cpecies
c      mass_ar(nbulk) - the mass of the given specie (in electron mass) 
c      X_ar(nbulk) = (fp/f)**2
c      Y_ar(nbulk) = fc/f ! for electron Ye should has the sign opposite to Yi
c      T_ar(nbulk)=temperature in eV  
c      npar_in - parallel index of refraction n.
c      nperp_in - perpendicular index of refraction n.
c
c     OUTPUT:
c      D_0      dispersion function evaluated at (mass,X,Y,T,Nll,np)
c               here T in [eV]d
c-----------------------------------------------------------
      implicit none
      include 'param.i'

c-----input
      integer nbulk
      real*8 X_ar(*),Y_ar(*),T_ar(*), mass_ar(*),
     & npar_in,nperp_in,r,z,phi
      integer i_deriv

c-----output
      real*8 D_0
    
c-----locals
      real*8 c,k4,npar,nperp,npars,nperps,nperp4,nperp6,
     &eps_perp,eps_par,eps_xy,
     &vt,vt_dc,
     &P_0,P_2,P_4,P_6,     
     &beta_p6,beta_p4,beta_p2,beta_p0
      
      integer i

    
c      write(*,*)'in  subroutine bonoli_Dispersion'
c      write(*,*)'npar_in,nperp_in',npar_in,nperp_in

c      write(*,*)'input arguments nbulk',nbulk
c      do i=1,nbulk
c        write(*,*)'i,mass_ar(i),X_ar(i),Y_ar(i),T_ar(i)',
c     &             i,mass_ar(i),X_ar(i),Y_ar(i),T_ar(i)
c      enddo

c      write(*,*)'npar_in,nperp_in',npar_in,nperp_in

ctest111120
      beta_p6=1.d0
      beta_p4=1.d0
      beta_p2=1.d0
      beta_p0=1.d0
c      beta_p6=0.d0
c      beta_p4=0.d0
c      beta_p2=0.d0

      c = 3.0d10                    !speed of light
      k4 = 4.19d7                   !constant for elec therm vel
 
      npar=npar_in
      nperp=nperp_in

      npars = npar**2
      nperps = nperp**2
      nperp4 = nperps**2
      nperp6 = nperp4*nperps
c-------------------------------------------------------------------
c     dielectric tesor calculations: eps_perp, eps_par, eps_xy
c     article formula 16(a-c)
c-------------------------------------------------------------------
c     electron terms
      eps_perp=1.d0+X_ar(1)/y_ar(1)**2  
      eps_par=1.d0-x_ar(1)
      eps_xy=x_ar(1)/y_ar(1) 
c     ion terms
      do i=2,nbulk
         eps_perp=eps_perp-x_ar(i)
         eps_par=eps_par-x_ar(i)
      enddo
c      write(*,*)'eps_perp,eps_perp,eps_xy',eps_perp,eps_perp,eps_xy
c--------------------------------------------------------------------
c     coefficients of dispersion function D_0 calculations:
c     P_0,P_2,P_4,P_6
c     article formula 15
c-------------------------------------------------------------------- 
      P_0=eps_par*((npars-eps_perp)**2 - eps_xy**2)
      P_2=(eps_perp+eps_par)*(npars-eps_perp)+eps_xy**2
      P_4=eps_perp
c--------------------------------------------------------------
      P_6=0.d0   
c--------------------------------------------------------------

      c = 3.0d10                    !speed of light
      k4 = 4.19d7                   !constant for elec therm v

c-----electron thermal term
      vt =k4*dsqrt(2.0d0*t_ar(1)/mass_ar(1)) !vt= sqrt(2.0*kT[eV]/mass) 
                                             !cm/sec
      vt_dc = vt / c
     
      P_6 = -(3.d0/8.d0)*(X_ar(1)/Y_ar(1)**4)*vt_dc**2
c-----plus ion's thermal term
      do i=2,nbulk 
         vt =k4*dsqrt(2.0d0*t_ar(i)/mass_ar(i)) !vt= sqrt(2.0*kT[eV]/mass) 
                                                !cm/sec 
         vt_dc = vt / c
         P_6 = P_6 - (3.d0/2.d0)*X_ar(i)*vt_dc**2
      enddo

c      P_6=0.d0   
c-------------------------------------------------------------
c     real dispersion function from article: formula 15
c--------------------------------------------------------------
      D_0 = P_6 * nperp6 + P_4 * nperp4 + P_2 * nperps + P_0
  
c      write(*,*)'bonoli_Disp P_6,P_4,P_2,P_0', P_6,P_4,P_2,P_0
c      write(*,*)'bonoli_Disp nperp6,nperp4,nperps,D_0',
c     &                       nperp6,nperp4,nperps,D_0
      return
      end

 
      subroutine bonoli_dervs(u,wf,
     &     dddz,dddr,dddphi,
     &     dddcnz,dddcnr,dddcm,dddw)   
c---------------------------------------------------------
c     analytical derivatives from Bonoli dispersion function
c---------------------------------------------------------
c     INPUT:
c     nbulk number of plasma species
c     u(6) ray coordinates
c     wf wave frequency
c     OUTPUT
c     Derivatives from dispersion function
c     dddz,dddr,dddphi
c     dddcnz,dddcnr,dddcm,dddw
c-----------------------------------------------  
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions_nml.i'
c-----input
      real*8 
     &u(6),
     &wf
c-----output
      real*8
     &dddz,dddr,dddphi,
     &dddcnz,dddcnr,dddcm,dddw
c-----externals
      real*8 b,cn,gamma1,x,y,tempe,
     &dxdr,dxdz,dxdphi,dydr,dydz,dydphi

c-----locals
      real*8  
     &z,r,phi,cnz,cnr,cm,
     &gam1,ds,dc,dcn,
     &npar,nperp,npars,nperps,nperp4,nperp6,
     &eps_perp,eps_par,eps_xy,
     &k4,c,vt,vt_dc,P_0,P_2,P_4,P_6,D_0,
     &d_nperp_dr,d_nperp_dz,d_nperp_dphi,
     &d_npar_dr,d_npar_dz,d_npar_dphi, 
     &d_nperp_dcnr,d_nperp_dcnz,d_nperp_dcm,  
     &d_npar_dcnr,d_npar_dcnz,d_npar_dcm,  
     &d_nperp_dw,d_npar_dw,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw,
     & d_P_6_dr,d_P_6_dz,d_P_6_dphi,
     & d_P_4_dr,d_P_4_dz,d_P_4_dphi,
     & d_P_2_dr,d_P_2_dz,d_P_2_dphi,
     & d_P_0_dr,d_P_0_dz,d_P_0_dphi,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw

      real*8, dimension(1:nbulk) :: x_ar,y_ar,t_ar
      real*8, dimension(1:nbulk) :: d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar 
      real*8, dimension(1:nbulk) :: d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar 
      integer i

c      write(*,*)'in bonoli_dervs'

      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)

      bmod=b(z,r,phi)  
      gam1=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam1)
      dc=dcos(gam1)
      dcn=cn(r,cnz,cnr,cm)
      npar=dcn*dc
      nperp=dcn*ds

      npars = npar**2
      nperps = nperp**2
      nperp4 = nperps**2
      nperp6 = nperp4*nperps

c      write(*,*)'npar,nperp',npar,nperp
c      write(*,*)'nperp6,nperp4,nperps',nperp6,nperp4,nperps
c------------------------------------
      do i=1,nbulk
         x_ar(i) = x(z,r,phi,i)
         y_ar(i) = y(z,r,phi,i)
         t_ar(i) = tempe(z,r,phi,i)*1.d3 !eV
         d_x_dr_ar(i)= dxdr(z,r,phi,i)
         d_x_dz_ar(i)= dxdz(z,r,phi,i)
         d_x_dphi_ar(i)=dxdphi(z,r,phi,i)
         d_y_dr_ar(i)= dydr(z,r,phi,i)
         d_y_dz_ar(i)= dydz(z,r,phi,i)
         d_y_dphi_ar(i)=dydphi(z,r,phi,i)
c         write(*,*)'i,x_ar(i),y_ar(i),t_ar(i)',
c     &              i,x_ar(i),y_ar(i),t_ar(i)
c         write(*,*)'d_x_dr_ar(i),d_x_dz_ar(i),d_x_dphi_ar(i)',
c     &              d_x_dr_ar(i),d_x_dz_ar(i),d_x_dphi_ar(i)
c         write(*,*)'d_y_dr_ar(i),d_y_dz_ar(i),d_y_dphi_ar(i)',
c     &              d_y_dr_ar(i),d_y_dz_ar(i),d_y_dphi_ar(i)
      enddo

c-------------------------------------------------------------------
c     dielectric tesor calculations: eps_perp, eps_par, eps_xy
c     article formula 16(a-c)
c-------------------------------------------------------------------
c     electron terms
      eps_perp=1.d0+X_ar(1)/y_ar(1)**2  
      eps_par=1.d0-x_ar(1)
      eps_xy=x_ar(1)/y_ar(1) 
c     ion terms
      do i=2,nbulk
         eps_perp=eps_perp-x_ar(i)
         eps_par=eps_par-x_ar(i)
      enddo

c--------------------------------------------------------------------
c     coefficients of dispersion function D_0 calculations:
c     P_0,P_2,P_4,P_6
c     article formula 15
c-------------------------------------------------------------------- 
      P_0=eps_par*((npars-eps_perp)**2 - eps_xy**2)
      P_2=(eps_perp+eps_par)*(npars-eps_perp)+eps_xy**2
      P_4=eps_perp
c--------------------------------------------------------------
      P_6=0.d0   
c--------------------------------------------------------------
      c = 3.0d10                    !speed of light
      k4 = 4.19d7                   !constant for elec therm v
c-----electron thermal term
      vt =k4*dsqrt(2.0d0*t_ar(1)/dmas(1)) !vt= sqrt(2.0*kT[eV]/mass) 
                                             !cm/sec
      vt_dc = vt / c
     
      P_6 = -(3.d0/8.d0)*(X_ar(1)/Y_ar(1)**4)*vt_dc**2
c-----plus ion's thermal term
      do i=2,nbulk 
         vt =k4*dsqrt(2.0d0*t_ar(i)/dmas(i)) !vt= sqrt(2.0*kT[eV]/mass) 
                                                !cm/sec 
         vt_dc = vt / c
         P_6 = P_6 - (3.d0/2.d0)*X_ar(i)*vt_dc**2
      enddo

c      P_6=0.d0   
c-------------------------------------------------------------
c     real dispersion function from article: formula 15
c--------------------------------------------------------------
      D_0 = P_6 * nperp6 + P_4 * nperp4 + P_2 * nperps + P_0

c      write(*,*)'bonoli_dervs P_6,P_4,P_2,P_0', P_6,P_4,P_2,P_0
c      write(*,*)'bonoli_drvs nperp6,nperp4,nperps,D_0',
c     &                       nperp6,nperp4,nperps,D_0

ctest-------------------------------------------------------------
c       call Numerical_der_dn_dw(r,z,phi,cnr,cnz,cm,
c     & d_npar_dw,d_nperp_dw)
c       write(*,*)'numerical d_npar_dw,d_nperp_dw',d_npar_dw,d_nperp_dw
ctest end---------------------------------------------------
      call dnd(z,r,phi,cnz,cnr,cm,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw)  
c       write(*,*)'analytical deriv from N'
c       write(*,*)'z,r,phi,cnz,cnr,cm',z,r,phi,cnz,cnr,cm
c       write(*,*)'d_nperp_dz,d_nperp_dr,d_nperp_dphi',
c     &            d_nperp_dz,d_nperp_dr,d_nperp_dphi
c       write(*,*)'d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm',
c     &            d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm     
c       write(*,*)'d_npar_dz,d_npar_dr,d_npar_dphi',
c     &            d_npar_dz,d_npar_dr,d_npar_dphi
c       write(*,*)'d_npar_dcnz,d_npar_dcnr,d_npar_dcm',
c     &            d_npar_dcnz,d_npar_dcnr,d_npar_dcm
c       write(*,*)'d_nperp_dw,d_npar_dw',d_nperp_dw,d_npar_dw
  
       call d_eps_d_rzphi_nrnznm(wf,x_ar,y_ar,
     &d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar, 
     &d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw)



c      write(*,*)'d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi',
c     &           d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi
c      write(*,*)'d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi',
c     &           d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi
c      write(*,*)'d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi',
c     &           d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi
c      write(*,*)'d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm',
c     &           d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm
c      write(*,*)'d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm',
c     &           d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm
c      write(*,*)'d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm',
c     &           d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm
c      write(*,*)'d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw',
c     &           d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c-----test numerical derivatives from eps
c      call Numerical_deriv_of_eps(r,z,phi,cnr,cnz,cm,
c     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
c     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
c     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,    
c     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw)

c      write(*,*)'Numerical derivatives of eps'

c      write(*,*)'d_eps_par_dr,d_eps_par_dz',
c     &           d_eps_par_dr,d_eps_par_dz
c      write(*,*)'d_eps_perp_dr,d_eps_perp_dz',
c     &           d_eps_perp_dr,d_eps_perp_dz
c      write(*,*)'d_eps_xy_dr,d_eps_xy_dz',
c     &           d_eps_xy_dr,d_eps_xy_dz
c      write(*,*)'d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw',
c     &           d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c-------------------------------------------------------------------
      call d_P0246_d_rzphi_nrnzcm(
     &npar,nperp,npars,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw, 
     &eps_perp,eps_par,eps_xy,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw,
     & d_P_6_dr,d_P_6_dz,d_P_6_dphi,
     & d_P_4_dr,d_P_4_dz,d_P_4_dphi,
     & d_P_2_dr,d_P_2_dz,d_P_2_dphi,
     & d_P_0_dr,d_P_0_dz,d_P_0_dphi,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm, 
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw)

c      write(*,*)'d_P_6_dr,d_P_6_dz,d_P_6_dphi',
c     &           d_P_6_dr,d_P_6_dz,d_P_6_dphi
c      write(*,*)'d_P_4_dr,d_P_4_dz,d_P_4_dphi',
c     &           d_P_4_dr,d_P_4_dz,d_P_4_dphi
c      write(*,*)'d_P_2_dr,d_P_2_dz,d_P_2_dphi',
c     &           d_P_2_dr,d_P_2_dz,d_P_2_dphi
c      write(*,*)'d_P_0_dr,d_P_0_dz,d_P_0_dphi',
c     &           d_P_0_dr,d_P_0_dz,d_P_0_dphi
c      write(*,*)'d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm',
c     &           d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm
c      write(*,*)'d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm',
c     &           d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm
c      write(*,*)'d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm',
c     &           d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm
c      write(*,*)'d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm',
c     &           d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm
c      write(*,*)'d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw',
c     &           d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw

      call d_P_6_d_rzphi_w(r,z,phi,wf,
     &d_P_6_dr,d_P_6_dz,d_P_6_dphi,d_P_6_dw)

c-----test of derivatives of P0246 by numerical differention
c      call  Numerical_deriv_of_P0246(r,z,phi,cnr,cnz,cm,
c     &d_P_0_dr,d_P_2_dr,d_P_4_dr,d_P_6_dr,
c     &d_P_0_dz,d_P_2_dz,d_P_4_dz,d_P_6_dz,
c     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
c     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
c     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
c     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
c     &d_P_0_dw,d_P_2_dw,d_P_4_dw,d_P_6_dw)

c      write(*,*)'numrerical derivatives of P0246'
c      write(*,*)'d_P_6_dr,d_P_6_dz',
c     &           d_P_6_dr,d_P_6_dz

c      write(*,*)'d_P_4_dr,d_P_4_dz',
c     &           d_P_4_dr,d_P_4_dz

c      write(*,*)'d_P_2_dr,d_P_2_dz',
c     &           d_P_2_dr,d_P_2_dz
c      write(*,*)'d_P_0_dr,d_P_0_dz',
c     &           d_P_0_dr,d_P_0_dz

c      write(*,*)'d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm',
c     &           d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm
c      write(*,*)'d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm',
c     &           d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm
c      write(*,*)'d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm',
c     &           d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm
c      write(*,*)'d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm',
c     &           d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm
c      write(*,*)'d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw',
c     &           d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw

c      stop 'lh_bonoli_disp.f'
c-------------------------------------------------------------------
      dddr= d_P_6_dr*nperp6 +6.d0*P_6 * nperp4*nperp*d_nperp_dr+
     &      d_P_4_dr*nperp4 +4.d0*P_4 * nperps*nperp*d_nperp_dr+
     &      d_P_2_dr*nperps +2.d0*P_2 * nperp*d_nperp_dr+
     &      d_P_0_dr

      dddz= d_P_6_dz*nperp6 +6.d0*P_6 * nperp4*nperp*d_nperp_dz+
     &      d_P_4_dz*nperp4 +4.d0*P_4 * nperps*nperp*d_nperp_dz+
     &      d_P_2_dz*nperps +2.d0*P_2 * nperp*d_nperp_dz+
     &      d_P_0_dz 
    
      dddphi= d_P_6_dphi*nperp6  +6.d0*P_6 * nperp4*nperp*d_nperp_dphi+
     &        d_P_4_dphi*nperp4 +4.d0*P_4 * nperps*nperp*d_nperp_dphi+
     &        d_P_2_dphi*nperps +2.d0*P_2 * nperp*d_nperp_dphi+
     &        d_P_0_dphi

      dddcnr= d_P_6_dcnr*nperp6 +6.d0*P_6 * nperp4*nperp*d_nperp_dcnr+
     &      d_P_4_dcnr*nperp4   +4.d0*P_4 * nperps*nperp*d_nperp_dcnr+
     &      d_P_2_dcnr*nperps   +2.d0*P_2 * nperp*d_nperp_dcnr+
     &      d_P_0_dcnr

      dddcnz= d_P_6_dcnz*nperp6 +6.d0*P_6 * nperp4*nperp*d_nperp_dcnz+
     &      d_P_4_dcnz*nperp4   +4.d0*P_4 * nperps*nperp*d_nperp_dcnz+
     &      d_P_2_dcnz*nperps   +2.d0*P_2 * nperp*d_nperp_dcnz+
     &      d_P_0_dcnz

      dddcm= d_P_6_dcm*nperp6  +6.d0*P_6 * nperp4*nperp*d_nperp_dcm+
     &      d_P_4_dcm*nperp4   +4.d0*P_4 * nperps*nperp*d_nperp_dcm+
     &      d_P_2_dcm*nperps   +2.d0*P_2 * nperp*d_nperp_dcm+
     &      d_P_0_dcm

      dddw= d_P_6_dw* nperp6   +6.d0*P_6 * nperp4*nperp*d_nperp_dw+
     &      d_P_4_dw*nperp4    +4.d0*P_4 * nperps*nperp*d_nperp_dw+
     &      d_P_2_dw*nperps    +2.d0*P_2 * nperp*d_nperp_dw+
     &      d_P_0_dw

c      write(*,*)'***d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw',
c     &              d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw  
c      write(*,*)'P_6,P_4,P_2,P_0',P_6,P_4,P_2,P_0
c      write(*,*)'nperp,nperps,nperp4,nperp6',nperp,nperps,nperp4,nperp6
c      write(*,*)'d_nperp_dw',d_nperp_dw
c      write(*,*)'d_P_6_dw* nperp6', d_P_6_dw* nperp6
c      write(*,*)'6.d0*P_6 * nperp4*nperp*d_nperp_dw',
c     &           6.d0*P_6 * nperp4*nperp*d_nperp_dw

c      write(*,*)'d_P_4_dw*nperp4',d_P_4_dw*nperp4
c      write(*,*)'4.d0*P_4 * nperps*nperp*d_nperp_dw',
c     &           4.d0*P_4 * nperps*nperp*d_nperp_dw

c      write(*,*)'d_P_2_dw*nperps',d_P_2_dw*nperps
c      write(*,*)'2.d0*P_2 * nperp*d_nperp_dw',
c     &           2.d0*P_2 * nperp*d_nperp_dw

c      write(*,*)' d_P_0_dw',d_P_0_dw

c      write(*,*)'bonoli analitical derivatives'
c      write(*,*)'dddcnz,dddcnr.dddcm',dddcnz,dddcnr,dddcm
c      write(*,*)'dddz,dddr,dddphi',dddz,dddr,dddphi
c      write(*,*)'dddw', dddw

c------test
c      write(*,*)'before dddw_analitic'
c      call dddw_analitic(u,wf,dddw,P_0,P_2,P_4,P_6)
c      write(*,*)'dddw', dddw
c      write(*,*)'after dddw_analitic'

c      stop 'lh_bonoli_disp.f'
      return
      end

      subroutine d_P0246_d_rzphi_nrnzcm(
     &npar,nperp,npars,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw, 
     &eps_perp,eps_par,eps_xy,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw,
     & d_P_6_dr,d_P_6_dz,d_P_6_dphi,
     & d_P_4_dr,d_P_4_dz,d_P_4_dphi,
     & d_P_2_dr,d_P_2_dz,d_P_2_dphi,
     & d_P_0_dr,d_P_0_dz,d_P_0_dphi,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw)

c-----calculates derivatives of coefficients P0,P2,P4,P6
c     with respect to r,z,phi,nr,nz,cm
c------------------------------------------------------------
c     INPUT:
c     &npar,nperp,npars,
c     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
c     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
c     &d_npar_dz,d_npar_dr,d_npar_dphi,
c     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
c     &d_nperp_dw,d_npar_dw, 
c     &eps_perp,eps_par,eps_xy,
c     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
c     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
c     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
c     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
c     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
c     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
c     d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c     OUTPUT
c     d_P_6_dr,d_P_6_dz,d_P_6_dphi,
c     d_P_4_dr,d_P_4_dz,d_P_4_dphi,
c     d_P_2_dr,d_P_2_dz,d_P_2_dphi,
c     d_P_0_dr,d_P_0_dz,d_P_0_dphi,
c     d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
c     d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
c     d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
c     d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
c     d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw
c-----------------------------------------------------
       implicit none
c-----input
      real*8 z,r,phi,cnz,cnr,cm ,
     &npar,nperp,npars,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw,
     &eps_perp,eps_par,eps_xy,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c-----output
      real*8
     & d_P_6_dr,d_P_6_dz,d_P_6_dphi,
     & d_P_4_dr,d_P_4_dz,d_P_4_dphi,
     & d_P_2_dr,d_P_2_dz,d_P_2_dphi,
     & d_P_0_dr,d_P_0_dz,d_P_0_dphi,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm , 
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw

c-----locals

      d_P_0_dr=d_eps_par_dr*((npars-eps_perp)**2 - eps_xy**2)+
     &eps_par*(2*(npars-eps_perp)*(2.d0*npar*d_npar_dr-d_eps_perp_dr)-
     &         2.d0*eps_xy*d_eps_xy_dr)

      d_P_0_dz=d_eps_par_dz*((npars-eps_perp)**2 - eps_xy**2)+
     &eps_par*(2*(npars-eps_perp)*(2.d0*npar*d_npar_dz-d_eps_perp_dz)-
     &         2.d0*eps_xy*d_eps_xy_dz)

      d_P_0_dphi=d_eps_par_dphi*((npars-eps_perp)**2 - eps_xy**2)+
     &eps_par*(2*(npars-eps_perp)*(2.d0*npar*d_npar_dphi-
     &                             d_eps_perp_dphi)-
     &         2.d0*eps_xy*d_eps_xy_dphi)

      d_P_2_dr=(d_eps_perp_dr+d_eps_par_dr)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2*npar*d_npar_dr-d_eps_perp_dr)+
     &                                    2.d0*eps_xy*d_eps_xy_dr

      d_P_2_dz=(d_eps_perp_dz+d_eps_par_dz)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2*npar*d_npar_dz-d_eps_perp_dz)+
     &                                    2.d0*eps_xy*d_eps_xy_dz

      d_P_2_dphi=(d_eps_perp_dphi+d_eps_par_dphi)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2*npar*d_npar_dphi-d_eps_perp_dphi)+
     &                                    2.d0*eps_xy*d_eps_xy_dphi

      d_P_4_dr=d_eps_perp_dr
      d_P_4_dz=d_eps_perp_dz
      d_P_4_dphi=d_eps_perp_dphi
c-----------------------------
      d_P_6_dr=0.d0
      d_P_6_dz=0.d0
      d_P_6_dphi=0.d0

c----------------------------------
      
      d_P_0_dcnr=d_eps_par_dcnr*((npars-eps_perp)**2 - eps_xy**2)+
     & eps_par*(2.d0*(npars-eps_perp)*
     &          (2.d0*npar*d_npar_dcnr-d_eps_perp_dcnr)-
     &          2.d0*eps_xy*d_eps_xy_dcnr)
     
      d_P_0_dcnz=d_eps_par_dcnz*((npars-eps_perp)**2 - eps_xy**2)+
     & eps_par*(2.d0*(npars-eps_perp)*
     &          (2.d0*npar*d_npar_dcnz-d_eps_perp_dcnz)-
     &          2.d0*eps_xy*d_eps_xy_dcnz)
     
      d_P_0_dcm=d_eps_par_dcm*((npars-eps_perp)**2 - eps_xy**2)+
     & eps_par*(2.d0*(npars-eps_perp)*
     &          (2.d0*npar*d_npar_dcm-d_eps_perp_dcm)-
     &          2.d0*eps_xy*d_eps_xy_dcm)
     
      d_P_2_dcnr=(d_eps_perp_dcnr+d_eps_par_dcnr)*(npars-eps_perp)+
     &(eps_perp+eps_par)*(2.d0*npar*d_npar_dcnr-d_eps_perp_dcnr)+
     &2.d0*eps_xy*d_eps_xy_dcnr

      d_P_2_dcnz=(d_eps_perp_dcnz+d_eps_par_dcnz)*(npars-eps_perp)+
     &(eps_perp+eps_par)*(2.d0*npar*d_npar_dcnz-d_eps_perp_dcnz)+
     &2.d0*eps_xy*d_eps_xy_dcnz

      d_P_2_dcm=(d_eps_perp_dcm+d_eps_par_dcm)*(npars-eps_perp)+
     &(eps_perp+eps_par)*(2.d0*npar*d_npar_dcm-d_eps_perp_dcm)+
     &2.d0*eps_xy*d_eps_xy_dcm
 
      d_P_4_dcnr=d_eps_perp_dcnr   
      d_P_4_dcnz=d_eps_perp_dcnz
      d_P_4_dcm=d_eps_perp_dcm


c-----------------------------
      d_P_6_dcnr=0.d0
      d_P_6_dcnz=0.d0
      d_P_6_dcm=0.d0 
c----------------------------------

      d_P_0_dw=d_eps_par_dw*((npars-eps_perp)**2 - eps_xy**2)+
     &          eps_par*(2.d0*(npars-eps_perp)*
     &         (2.d0*npar*d_npar_dw-d_eps_perp_dw)-
     &          2.d0*eps_xy*d_eps_xy_dw)

      d_P_2_dw=(d_eps_perp_dw+d_eps_par_dw)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2.d0*npar*d_npar_dw-d_eps_perp_dw)+
     &         2.d0*eps_xy*d_eps_xy_dw


      d_P_4_dw=d_eps_perp_dw


      d_P_6_dw=0.d0

c----------------------------------

c      write(*,*)'d_P0246_d_rzphi_nrnzcm'
c      write(*,*)'npar,nperp,npars',npar,nperp,npars
c      write(*,*)'&d_nperp_dw,d_npar_dw',d_nperp_dw,d_npar_dw
c      write(*,*)'eps_perp,eps_par,eps_xy',eps_perp,eps_par,eps_xy
c      write(*,*)'d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw',
c     &           d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c      write(*,*)'d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw',
c     &           d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw
      return
      end

      subroutine d_eps_d_rzphi_nrnznm(wf,x_ar,y_ar,
     &d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar, 
     &d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar, 
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw)
c---------------------------------------------------------
c     calculates derivatives of dielectric tensor
c     with respect to r,z,phi,cnr,cnz,cm
c------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      
c-----input
      real*8 wf,  
     &x_ar(*),y_ar(*), 
     &d_x_dr_ar(*),d_x_dz_ar(*),d_x_dphi_ar(*), 
     &d_y_dr_ar(*),d_y_dz_ar(*),d_y_dphi_ar(*) 
c-----output
      real*8
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c-----locals
      integer i
      real*8, dimension(1:nbulk) :: d_eps_par_d_x_ar,d_eps_par_d_y_ar
      real*8, dimension(1:nbulk) :: d_eps_perp_d_x_ar,d_eps_perp_d_y_ar
      real*8, dimension(1:nbulk) :: d_eps_xy_d_x_ar,d_eps_xy_d_y_ar

      real*8 dxdw,dydw

      do i=1,nbulk
         if (i.eq.1) then
c-----------electron term        
            d_eps_perp_d_x_ar(i)    =  1.d0/y_ar(i)**2  
            d_eps_par_d_x_ar(i)     = -1.d0
            d_eps_xy_d_x_ar(i)      =  1.d0/y_ar(i)

            d_eps_perp_d_y_ar(i)    =  -2.d0*x_ar(i)/y_ar(i)**3   
            d_eps_par_d_y_ar(i)     =   0.d0
            d_eps_xy_d_y_ar(i)      =   -x_ar(i)/y_ar(i)**2
         else
c-----------ion terms
            d_eps_perp_d_x_ar(i)    =   -1.d0
            d_eps_par_d_x_ar(i)     =   -1.d0
            d_eps_xy_d_x_ar(i)      =   0.d0

            d_eps_perp_d_y_ar(i)    =   0.d0 
            d_eps_par_d_y_ar(i)     =   0.d0
            d_eps_xy_d_y_ar(i)      =   0.d0
           
         endif
      enddo !nbulk

      d_eps_par_dr=0.d0
      d_eps_par_dz=0.d0
      d_eps_par_dphi=0.d0
      d_eps_perp_dr=0.d0
      d_eps_perp_dz=0.d0
      d_eps_perp_dphi=0.d0
      d_eps_xy_dr=0.d0
      d_eps_xy_dz=0.d0
      d_eps_xy_dphi=0.d0
      d_eps_par_dcnr=0.d0
      d_eps_par_dcnz=0.d0
      d_eps_par_dcm=0.d0
      d_eps_perp_dcnr=0.d0
      d_eps_perp_dcnz=0.d0
      d_eps_perp_dcm=0.d0
      d_eps_xy_dcnr=0.d0
      d_eps_xy_dcnz=0.d0
      d_eps_xy_dcm=0.d0
      d_eps_par_dw=0.d0
      d_eps_perp_dw=0.d0
      d_eps_xy_dw=0.d0

      do i=1,nbulk

         d_eps_par_dr=d_eps_par_dr + d_eps_par_d_x_ar(i)*d_x_dr_ar(i)+
     &                               d_eps_par_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_par_dz=d_eps_par_dz + d_eps_par_d_x_ar(i)*d_x_dz_ar(i)+
     &                               d_eps_par_d_y_ar(i)*d_y_dz_ar(i)

         d_eps_par_dphi=d_eps_par_dphi+
     &                  d_eps_par_d_x_ar(i)*d_x_dphi_ar(i)+
     &                  d_eps_par_d_y_ar(i)*d_y_dphi_ar(i)


         d_eps_perp_dr=d_eps_perp_dr+d_eps_perp_d_x_ar(i)*d_x_dr_ar(i)+
     &                               d_eps_perp_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_perp_dz=d_eps_perp_dz+d_eps_perp_d_x_ar(i)*d_x_dz_ar(i)+
     &                               d_eps_perp_d_y_ar(i)*d_y_dz_ar(i)

         d_eps_perp_dphi=d_eps_perp_dphi+
     &                   d_eps_perp_d_x_ar(i)*d_x_dphi_ar(i)+
     &                   d_eps_perp_d_y_ar(i)*d_y_dphi_ar(i)


         d_eps_xy_dr=d_eps_xy_dr+d_eps_xy_d_x_ar(i)*d_x_dr_ar(i)+
     &                           d_eps_xy_d_y_ar(i)*d_y_dr_ar(i)


         d_eps_xy_dz=d_eps_xy_dz+d_eps_xy_d_x_ar(i)*d_x_dz_ar(i)+
     &                           d_eps_xy_d_y_ar(i)*d_y_dz_ar(i)

         d_eps_xy_dphi=d_eps_xy_dphi+
     &                 d_eps_xy_d_x_ar(i)*d_x_dphi_ar(i)+
     &                 d_eps_xy_d_y_ar(i)*d_y_dphi_ar(i)

         d_eps_par_dcnr=0.d0
         d_eps_par_dcnz=0.d0
         d_eps_par_dcm=0.d0

         d_eps_perp_dcnr=0.d0        
         d_eps_perp_dcnz=0.d0
         d_eps_perp_dcm=0.d0
        
         d_eps_xy_dcnr=0.d0
         d_eps_xy_dcnz=0.d0
         d_eps_xy_dcm=0.d0

         dxdw=-2.d0*x_ar(i)/wf 
         dydw=-y_ar(i)/wf 
            
         d_eps_par_dw=d_eps_par_dw + d_eps_par_d_x_ar(i)*dxdw+
     &                               d_eps_par_d_y_ar(i)*dydw

         d_eps_perp_dw=d_eps_perp_dw+d_eps_perp_d_x_ar(i)*dxdw+
     &                               d_eps_perp_d_y_ar(i)*dydw

         d_eps_xy_dw=d_eps_xy_dw +   d_eps_xy_d_x_ar(i)*dxdw+
     &                               d_eps_xy_d_y_ar(i)*dydw

      enddo

      return
      end


c++++++++++++++++++++++++++++++++++++++++++++++++++++++
c      Following functions are tor test of derivatives
c------------------------------------------------------
      subroutine calc_P0246(r,z,phi,cnr,cnz,cm,P_0,P_2,P_4,P_6)
c-------------------------------------------------
c     calculate coefficients P1246 of the Bonoli dispersion function
c---------------------------------------------------------------
c      INPUT:
c     r,z,phi,cnr,cnr,cm
c      OUTPUT:       
c     P_0,P_2,P_4,P_6
c----------------------------------------------------------------- 
      implicit none
      include 'param.i'
      include 'one.i'
      
c-----input
      real*8 r,z,phi,cnr,cnz,cm
c-----output
      real*8 P_0,P_2,P_4,P_6
c-----externals
      real*8 b,gamma1,x,y,tempe,cn
c-----locals
      real*8 
     &gam1,ds,dc,npar,nperp,dcn,npars,
     &eps_perp,eps_par,eps_xy

      integer i
       
      real*8, dimension(1:nbulk) :: x_ar,y_ar,t_ar

      bmod=b(z,r,phi)
      gam1=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam1)
      dc=dcos(gam1)
      dcn=cn(r,cnz,cnr,cm)
      npar=dcn*dc
      nperp=dcn*ds

      npars=npar**2

      do i=1,nbulk
         x_ar(i) = x(z,r,phi,i)
         y_ar(i) = y(z,r,phi,i)
         t_ar(i) = tempe(z,r,phi,i)*1.d3 !eV
      enddo

      call eps_bonoli(x_ar,y_ar,eps_perp,eps_par,eps_xy) 

      P_0=eps_par*((npars-eps_perp)**2 - eps_xy**2)
      P_2=(eps_perp+eps_par)*(npars-eps_perp)+eps_xy**2
      P_4=eps_perp

      P_6=0.d0   
      
      return
      end


       subroutine eps_bonoli(x_ar,y_ar,eps_perp,eps_par,eps_xy) 
c------calculte Bonoli tensor
c      INPUT:
c      OUTPUT:
c      eps_perp,eps_par,eps_xy 
       implicit none
       include 'param.i'
       include 'one.i'
c------input
       real*8 x_ar(*),y_ar(*)
c------output
       real*8  eps_perp,eps_par,eps_xy 
c------locals
       integer i
c-------------------------------------------------------------------
c      dielectric tesor calculations: eps_perp, eps_par, eps_xy
c      article formula 16(a-c)
c-------------------------------------------------------------------
c      electron terms
       eps_perp=1.d0+X_ar(1)/y_ar(1)**2  
       eps_par=1.d0-x_ar(1)
       eps_xy=x_ar(1)/y_ar(1) 
c      ion terms
       do i=2,nbulk
         eps_perp=eps_perp-x_ar(i)
         eps_par=eps_par-x_ar(i)
       enddo

       return
       end
 

      subroutine Numerical_deriv_of_P0246(r,z,phi,cnr,cnz,cm,
     &d_P_0_d_r,d_P_2_d_r,d_P_4_d_r,d_P_6_d_r,
     &d_P_0_d_z,d_P_2_d_z,d_P_4_d_z,d_P_6_d_z,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
     &d_P_0_d_w,d_P_2_d_w,d_P_4_d_w,d_P_6_d_w)
c----------------------------------------------------------
c     calculates numerical derivatives of P0.P2,P4,P6
c     with respectr to r,z,phi
c-----------------------------------------------
c     INPUT:
c     r,z,phi
c     OUTPUT:
c     d_P_0_d_r,d_P_2_d_r,d_P_4_d_r,d_P_6_d_r,
c     d_P_0_d_z,d_P_2_d_z,d_P_4_d_z,d_P_6_d_z
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 r,z,phi,cnr,cnz,cm
c-----output
      real*8 
     &d_P_0_d_r,d_P_2_d_r,d_P_4_d_r,d_P_6_d_r,
     &d_P_0_d_z,d_P_2_d_z,d_P_4_d_z,d_P_6_d_z,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
     &d_P_0_d_w,d_P_2_d_w,d_P_4_d_w,d_P_6_d_w
c-----locals 
      integer i
      real*8 step,
     &rp,rm,zp,zm,phip,phim,cnrp,cnrm,cnzp,cnzm,cmp,cmm,
     &P_0_p,P_2_p,P_4_p,P_6_p,
     &P_0_m,P_2_m,P_4_m,P_6_m,
     &hw,hfrqnc,frqncpl,df,frqncmn,
     &cnzplus,cnzminus,
     &cnrplus,cnrminus,
     &cmplus,cmminus
      real*8, dimension(1:nbulk) :: vp,wp
c-----externals
      real*8 b
      step=1.d-5

      hw=step*1.0d0
      hfrqnc=hw
      hfrqnc=hfrqnc*frqncy

      rp=r+step
      rm=r-step

      call calc_P0246(rp,z,phi,cnr,cnz,cm,P_0_p,P_2_p,P_4_p,P_6_p)
      call calc_P0246(rm,z,phi,cnr,cnz,cm,P_0_m,P_2_m,P_4_m,P_6_m)

      d_P_0_d_r=(P_0_p-P_0_m)/(2.d0*step)
      d_P_2_d_r=(P_2_p-P_2_m)/(2.d0*step)
      d_P_4_d_r=(P_4_p-P_4_m)/(2.d0*step)   
      d_P_6_d_r=(P_6_p-P_6_m)/(2.d0*step)
 
      zp=z+step
      zm=z-step
     
      call calc_P0246(r,zp,phi,cnr,cnz,cm,P_0_p,P_2_p,P_4_p,P_6_p)
      call calc_P0246(r,zm,phi,cnr,cnz,cm,P_0_m,P_2_m,P_4_m,P_6_m)

      d_P_0_d_z=(P_0_p-P_0_m)/(2.d0*step)
      d_P_2_d_z=(P_2_p-P_2_m)/(2.d0*step)
      d_P_4_d_z=(P_4_p-P_4_m)/(2.d0*step)   
      d_P_6_d_z=(P_6_p-P_6_m)/(2.d0*step)

      cnrp=cnr+step
      cnrm=cnr-step

      call calc_P0246(r,z,phi,cnrp,cnz,cm,P_0_p,P_2_p,P_4_p,P_6_p)
      call calc_P0246(r,z,phi,cnrm,cnz,cm,P_0_m,P_2_m,P_4_m,P_6_m)

      d_P_0_dcnr=(P_0_p-P_0_m)/(2.d0*step)
      d_P_2_dcnr=(P_2_p-P_2_m)/(2.d0*step)
      d_P_4_dcnr=(P_4_p-P_4_m)/(2.d0*step)   
      d_P_6_dcnr=(P_6_p-P_6_m)/(2.d0*step)


      cnzp=cnz+step
      cnzm=cnz-step

      call calc_P0246(r,z,phi,cnr,cnzp,cm,P_0_p,P_2_p,P_4_p,P_6_p)
      call calc_P0246(r,z,phi,cnr,cnzm,cm,P_0_m,P_2_m,P_4_m,P_6_m)

      d_P_0_dcnz=(P_0_p-P_0_m)/(2.d0*step)
      d_P_2_dcnz=(P_2_p-P_2_m)/(2.d0*step)
      d_P_4_dcnz=(P_4_p-P_4_m)/(2.d0*step)   
      d_P_6_dcnz=(P_6_p-P_6_m)/(2.d0*step)


      cmp=cm+step
      cmm=cm-step

      call calc_P0246(r,z,phi,cnr,cnz,cmp,P_0_p,P_2_p,P_4_p,P_6_p)
      call calc_P0246(r,z,phi,cnr,cnz,cmm,P_0_m,P_2_m,P_4_m,P_6_m)

      d_P_0_dcm=(P_0_p-P_0_m)/(2.d0*step)
      d_P_2_dcm=(P_2_p-P_2_m)/(2.d0*step)
      d_P_4_dcm=(P_4_p-P_4_m)/(2.d0*step)   
      d_P_6_dcm=(P_6_p-P_6_m)/(2.d0*step)

      bmod=b(z,r,phi)
      do i=1,nbulk
	  vp(i)=v(i)
	  wp(i)=w(i)
      enddo 

      frqncpl=frqncy+hfrqnc
      df=frqncy/frqncpl
      do i=1,nbulk
	  v(i)=vp(i)*df* df
	  w(i)=wp(i)*df
      enddo 
      cnrplus=cnr*df
      cnzplus=cnz*df
      cmplus=cm*df
      call calc_P0246(r,z,phi,cnrplus,cnzplus,cmplus,
     &                P_0_p,P_2_p,P_4_p,P_6_p)


      frqncmn=frqncy-hfrqnc
      df=frqncy/frqncmn
      do i=1,nbulk
	  v(i)=vp(i)*df*df
	  w(i)=wp(i)*df
      enddo
      cnrminus=cnr*df
      cnzminus=cnz*df
      cmminus=cm*df
      call calc_P0246(r,z,phi, cnrminus,cnzminus,cmminus,
     &P_0_m,P_2_m,P_4_m,P_6_m)

      d_P_0_d_w=(P_0_p-P_0_m)/(2.d0*hfrqnc)
      d_P_2_d_w=(P_2_p-P_2_m)/(2.d0*hfrqnc)
      d_P_4_d_w=(P_4_p-P_4_m)/(2.d0*hfrqnc)
      d_P_6_d_w=(P_6_p-P_6_m)/(2.d0*hfrqnc)

      do i=1,nbulk
	  v(i)=vp(i)
	  w(i)=wp(i)
      enddo

      return
      end

      subroutine Numerical_deriv_of_eps(r,z,phi,cnr,cnz,cm,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,    
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw)
c----------------------------------------------------------
c     calculates numerical derivatives of Bonoli eps
c     with respectr to r,z,phi,cnr,cnz,cm,
c-----------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'

c-----input
      real*8 r,z,phi,cnr,cnz,cm
c-----output
      real*8  
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c-----externals
      real*8 b,x,y
c-----locals
      real*8 step,
     &rp,rm,zp,zm,phip,phim,cnrp,cnrm,cnzp,cnzm,cmp,cmm,
     &eps_perp_p,eps_par_p,eps_xy_p,
     &eps_perp_m,eps_par_m,eps_xy_m,
     &hw,hfrqnc,frqncpl,df,frqncmn

      integer i
      real*8, dimension(1:nbulk) :: x_ar,y_ar,vp,wp


      step=1.d-5

      hw=step*1.0d0
      hfrqnc=hw
      hfrqnc=hfrqnc*frqncy

      rp=r+step
      bmod=b(z,rp,phi)
      do i=1,nbulk
         x_ar(i) = x(z,rp,phi,i)
         y_ar(i) = y(z,rp,phi,i)
      enddo

      call eps_bonoli(x_ar,y_ar,eps_perp_p,eps_par_p,eps_xy_p)

      rm=r-step
      bmod=b(z,rm,phi)
      do i=1,nbulk
         x_ar(i) = x(z,rm,phi,i)
         y_ar(i) = y(z,rm,phi,i)
      enddo
      call eps_bonoli(x_ar,y_ar,eps_perp_m,eps_par_m,eps_xy_m)

      d_eps_par_dr=(eps_par_p-eps_par_m)/(2.d0*step)    
      d_eps_perp_dr=(eps_perp_p-eps_perp_m)/(2.d0*step)
      d_eps_xy_dr=(eps_xy_p-eps_xy_m)/(2.d0*step)
c--------------------------------------------------------------
      zp=z+step
      bmod=b(zp,r,phi)
      do i=1,nbulk
         x_ar(i) = x(zp,r,phi,i)
         y_ar(i) = y(zp,r,phi,i)
      enddo
      call eps_bonoli(x_ar,y_ar,eps_perp_p,eps_par_p,eps_xy_p)

      zm=z-step
      bmod=b(zm,r,phi)
      do i=1,nbulk
         x_ar(i) = x(zm,r,phi,i)
         y_ar(i) = y(zm,r,phi,i)
      enddo
      call eps_bonoli(x_ar,y_ar,eps_perp_m,eps_par_m,eps_xy_m)
    
      d_eps_par_dz=(eps_par_p-eps_par_m)/(2.d0*step)    
      d_eps_perp_dz=(eps_perp_p-eps_perp_m)/(2.d0*step)
      d_eps_xy_dz=(eps_xy_p-eps_xy_m)/(2.d0*step)

      bmod=b(z,r,phi)
      do i=1,nbulk
	  vp(i)=v(i)
	  wp(i)=w(i)
      enddo

      frqncpl=frqncy+hfrqnc
      df=frqncy/frqncpl
      do i=1,nbulk
	  v(i)=vp(i)*df* df
	  w(i)=wp(i)*df
      enddo
      do i=1,nbulk
         x_ar(i) = x(zm,r,phi,i)
         y_ar(i) = y(zm,r,phi,i)
      enddo
      call eps_bonoli(x_ar,y_ar,eps_perp_p,eps_par_p,eps_xy_p)
  
      frqncmn=frqncy-hfrqnc
      df=frqncy/frqncmn
      do i=1,nbulk
	  v(i)=vp(i)*df*df
	  w(i)=wp(i)*df
      enddo
      do i=1,nbulk
         x_ar(i) = x(zm,r,phi,i)
         y_ar(i) = y(zm,r,phi,i)
      enddo
      call eps_bonoli(x_ar,y_ar,eps_perp_m,eps_par_m,eps_xy_m)

      d_eps_par_dw=(eps_par_p-eps_par_m)/(2.d0*hfrqnc)  
      d_eps_perp_dw=(eps_perp_p-eps_perp_m)/(2.d0*hfrqnc)
      d_eps_xy_dw=(eps_xy_p-eps_xy_m)/(2.d0*hfrqnc)

      do i=1,nbulk
	  v(i)=vp(i)
	  w(i)=wp(i)
      enddo

      return
      end

      subroutine Numerical_der_dn_dw(r,z,phi,cnr,cnz,cm,
     & d_npar_dw,d_nperp_dw)
c-----calculates numerical derivatives
c     d_npar_dw and d_nperp_dw
c--------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 r,z,phi,cnr,cnz,cm
c-----output
      real*8 d_npar_dw,d_nperp_dw
c-----locals
      real*8 step,
     &hw,hfrqnc,frqncpl,df,frqncmn,
     &gam1,ds,dc,dcn,
     &cnzp,cnrp,cmp,cnzm,cnrm,cmm,
     &npar_p,nperp_p,npar_m,nperp_m
c-----external 
      real*8 b,gamma1,cn
      integer i

      real*8, dimension(1:nbulk) :: vp,wp 

      bmod=b(z,r,phi)  
      do i=1,nbulk
	  vp(i)=v(i)
	  wp(i)=w(i)
      enddo 

      step=1.d-5
      hw=step*1.0d0
      hfrqnc=hw
      hfrqnc=hfrqnc*frqncy
      
   
      frqncpl=frqncy+hfrqnc
      df=frqncy/frqncpl
      do i=1,nbulk
	  v(i)=vp(i)*df* df
	  w(i)=wp(i)*df
      enddo 
      cnzp=cnz*df
      cnrp=cnr*df
      cmp= cm*df

      gam1=gamma1(z,r,phi,cnzp,cnrp,cmp)
      ds=dsin(gam1)
      dc=dcos(gam1) 

      dcn=cn(r,cnzp,cnrp,cmp)
      npar_p=dcn*dc
      nperp_p=dcn*ds

c-------------------------------------
      frqncmn=frqncy-hfrqnc
      df=frqncy/frqncmn
      do i=1,nbulk
	  v(i)=vp(i)*df* df
	  w(i)=wp(i)*df
      enddo 
      cnzm=cnz*df
      cnrm=cnr*df
      cmm= cm*df

      gam1=gamma1(z,r,phi,cnzm,cnrm,cmm)
      ds=dsin(gam1)
      dc=dcos(gam1) 

      dcn=cn(r,cnzm,cnrm,cmm)
      npar_m=dcn*dc
      nperp_m=dcn*ds
c--------------------------------------------------
      d_npar_dw=(npar_p- npar_m)/(2.d0*hfrqnc)
      d_nperp_dw=(nperp_p- nperp_m)/(2.d0*hfrqnc)
c------------------------------------------------------
      do i=1,nbulk
	  v(i)=vp(i)
	  w(i)=wp(i)
      enddo

      return
      end


      subroutine d_P_6_d_rzphi_w(r,z,phi,wf,
     &d_P_6_dr,d_P_6_dz,d_P_6_dphi,d_P_6_dw)
c-----------------------------------------------------------------
c     Calculate derivatives: d_P_0/d_r,d_P_0/d_z,d_P_0/d_phi,d_P_0/d_w
c-------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'

c-----input
      real*8 r,z,phi,wf
c-----output
      real*8 d_P_6_dr,d_P_6_dz,d_P_6_dphi,d_P_6_dw
c-----externals
      real*8 b,x,y,tempe,dxdr,dxdz,dxdphi,dydr,dydz,dydphi,
     &dtempdr,dtempdz
c-----locals
      integer i

      real*8 c,k4,vt,vt_dc,p,vt_dc_s,dv_dr,dv_dz,dv_dphi
      real*8, dimension(1:nbulk) :: x_ar,y_ar,t_ar,
     &        d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar,
     &        d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar,
     &        d_t_dr_ar,d_t_dz_ar,d_t_dphi_ar

      bmod=b(z,r,phi)
  
      do i=1,nbulk
         x_ar(i) = x(z,r,phi,i)
         y_ar(i) = y(z,r,phi,i)
         t_ar(i) = tempe(z,r,phi,i)*1.d3 !eV
         d_x_dr_ar(i)= dxdr(z,r,phi,i)
         d_x_dz_ar(i)= dxdz(z,r,phi,i)
         d_x_dphi_ar(i)=dxdphi(z,r,phi,i)
         d_y_dr_ar(i)= dydr(z,r,phi,i)
         d_y_dz_ar(i)= dydz(z,r,phi,i)
         d_y_dphi_ar(i)=dydphi(z,r,phi,i) 
         d_t_dr_ar(i)= dtempdr(z,r,phi,i)*1.d3 
         d_t_dz_ar(i)= dtempdz(z,r,phi,i)*1.d3 
         d_t_dphi_ar(i)=0.d0
      enddo

      d_P_6_dr=0.d0
      d_P_6_dz=0.d0
      d_P_6_dphi=0.d0
      d_P_6_dw=0.d0

      c = 3.0d10                    !speed of light
      k4 = 4.19d7                   !constant for elec therm v
c-----electron terms
      vt =k4*dsqrt(2.0d0*t_ar(1)/dmas(1)) !vt= sqrt(2.0*kT[eV]/mass) 
                                             !cm/sec 
      vt_dc = vt / c
      vt_dc_s= vt_dc*vt_dc

      p=k4/dsqrt(2.d0*t_ar(1)*dmas(1))

      dv_dr=p* d_t_dr_ar(1)
      dv_dz=p* d_t_dz_ar(1)
      dv_dphi=0.d0

      d_P_6_dr=- (3.d0/8.d0)*(
     & (d_x_dr_ar(1)/y_ar(1)**4)*vt_dc_s -
     & (4.d0*x_ar(1)/y_ar(1)**5)*d_y_dr_ar(1)*vt_dc_s +
     &  (x_ar(1)/y_ar(1)**4)*(2.d0*vt_dc*(dv_dr/c)))

      d_P_6_dz=- (3.d0/8.d0)*(
     & (d_x_dz_ar(1)/y_ar(1)**4)*vt_dc_s -
     & (4.d0*x_ar(1)/y_ar(1)**5)*d_y_dz_ar(1)*vt_dc_s +
     &  (x_ar(1)/y_ar(1)**4)*(2.d0*vt_dc*(dv_dz/c)))

      d_P_6_dw=- (3.d0/8.d0)*
     &     (2.d0*x_ar(1)/y_ar(1)**4)*vt_dc_s/wf

c-----ion terms
      do i=2,nbulk
         vt =k4*dsqrt(2.0d0*t_ar(i)/dmas(i)) !vt= sqrt(2.0*kT[eV]/mass) 
                                                !cm/sec 
         vt_dc = vt / c
         vt_dc_s= vt_dc*vt_dc

         p=k4/dsqrt(2.d0*t_ar(i)*dmas(i))

         dv_dr=p* d_t_dr_ar(i)
         dv_dz=p* d_t_dz_ar(i)
         dv_dphi=0.d0

         d_P_6_dr = d_P_6_dr - (3.d0/2.d0)*(
     &                         d_x_dr_ar(i)*vt_dc_s +
     &                         x_ar(i)*(2.d0*vt_dc*(dv_dr/c)))

         
         d_P_6_dz = d_P_6_dz - (3.d0/2.d0)*(
     &                         d_x_dz_ar(i)*vt_dc_s +
     &                         x_ar(i)*(2.d0*vt_dc*(dv_dz/c)))

         d_P_6_dw = d_P_6_dw -(3.d0/2.d0)*
     &                         (-2.d0*x_ar(i)*vt_dc_s)/wf
      enddo

      return
      end

c==========================================================
c     test dddw
c----------------------------------------------------------
      subroutine dddw_analitic(u,wf,dddw,P_0,P_2,P_4,P_6)
c-----analitical dddw
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions_nml.i'
c-----input
      real*8 
     &u(6),
     &wf,
     &P_0,P_2,P_4,P_6
c-----output
      real*8 dddw
c-----locals
      real*8  
     &z,r,phi,cnz,cnr,cm, 
     &gam1,ds,dc,dcn,
     &npar,nperp,npars,                                  
     &nperp6,nperp4,nperps,    
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw,
cBH180620     &d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar,              !
cBH180620     &d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar,              !
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &eps_perp,eps_par,eps_xy, 
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw,
     &d_P_0_dw,d_P_2_dw,d_P_4_dw,d_P_6_dw            !
       
      integer i
c-----externals
      real*8 b,x,y,gamma1,cn
      real*8, dimension(1:nbulk) :: x_ar,y_ar!
cBH180620
      real*8, dimension(1:nbulk) :: d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar 
      real*8, dimension(1:nbulk) :: d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar 

      write(*,*)'in dddw_analitic'

      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)
      
      bmod=b(z,r,phi)  
      gam1=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam1)
      dc=dcos(gam1)
      dcn=cn(r,cnz,cnr,cm)
      npar=dcn*dc
      nperp=dcn*ds

      npars = npar**2
      nperps = nperp**2
      nperp4 = nperps**2
      nperp6 = nperp4*nperps

      call  Numerical_der_dn_dw(r,z,phi,cnr,cnz,cm,
     & d_npar_dw,d_nperp_dw)

      write(*,*)'numerical derivatives d_nperp_dw,d_npar_dw'
      write(*,*)'d_nperp_dw,d_npar_dw',d_nperp_dw,d_npar_dw

      call dnd(z,r,phi,cnz,cnr,cm,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw)

      write(*,*)'analitical derivatives d_nperp_dw,d_npar_dw'
      write(*,*)'d_nperp_dw,d_npar_dw',d_nperp_dw,d_npar_dw
c-------------------------------------------------------------      
      call Numerical_deriv_of_eps(r,z,phi,cnr,cnz,cm,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,    
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw)
      write(*,*)'numerical derivatives d_eps_dw'
      write(*,*)'d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw',
     &           d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw

      do i=1,nbulk
         x_ar(i) = x(z,r,phi,i)
         y_ar(i) = y(z,r,phi,i)
      enddo

      call eps_bonoli(x_ar,y_ar,eps_perp,eps_par,eps_xy) 

      call d_eps_d_rzphi_nrnznm(wf,x_ar,y_ar,
     &d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar, 
     &d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar, 
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw)

      write(*,*)'analitical derivatives d_eps_dw'

      write(*,*)'d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw',
     &           d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw

c-------------------------------------------------------------
      call  Numerical_deriv_of_P0246_dw(r,z,phi,cnr,cnz,cm,
     &d_P_0_dw,d_P_2_dw,d_P_4_dw,d_P_6_dw)

      write(*,*)'numerical derivatives of P0246_dw'  
      write(*,*)'d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw',
     &           d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw

      call d_P0246_dw_analitic(
     &npar,nperp,npars,   
     &d_nperp_dw,d_npar_dw, 
     &eps_perp,eps_par,eps_xy,     
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw,
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw)

      write(*,*)'analitical derivatives of P0246_dw'  
      write(*,*)'d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw',
     &           d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw
c-------------------------------------------------------------
c       d_P_6_dw=d_P_6_dw*frqncy
c       d_P_4_dw=d_P_4_dw*frqncy
c       d_P_2_dw=d_P_2_dw*frqncy
c       d_P_0_dw=d_P_0_dw*frqncy
      dddw= d_P_6_dw* nperp6   +6.d0*P_6 * nperp4*nperp*d_nperp_dw+
     &      d_P_4_dw*nperp4    +4.d0*P_4 * nperps*nperp*d_nperp_dw+
     &      d_P_2_dw*nperps    +2.d0*P_2 * nperp*d_nperp_dw+
     &      d_P_0_dw
      dddw=dddw*frqncy
      write(*,*)'P_6,P_4,P_2,P_0',P_6,P_4,P_2,P_0
      write(*,*)'nperp,nperps,nperp4,nperp6',nperp,nperps,nperp4,nperp6
      write(*,*)'d_nperp_dw',d_nperp_dw
      write(*,*)'d_P_6_dw* nperp6', d_P_6_dw* nperp6
      write(*,*)'6.d0*P_6 * nperp4*nperp*d_nperp_dw',
     &           6.d0*P_6 * nperp4*nperp*d_nperp_dw

      write(*,*)'d_P_4_dw*nperp4',d_P_4_dw*nperp4
      write(*,*)'4.d0*P_4 * nperps*nperp*d_nperp_dw',
     &           4.d0*P_4 * nperps*nperp*d_nperp_dw

      write(*,*)'d_P_2_dw*nperps',d_P_2_dw*nperps
      write(*,*)'2.d0*P_2 * nperp*d_nperp_dw',
     &           2.d0*P_2 * nperp*d_nperp_dw

      write(*,*)' d_P_0_dw',d_P_0_dw
   
      write(*,*)'danalitic ddw', dddw

      return
      end

      subroutine Numerical_deriv_of_P0246_dw(r,z,phi,cnr,cnz,cm,   
     &d_P_0_d_w,d_P_2_d_w,d_P_4_d_w,d_P_6_d_w)
c----------------------------------------------------------
c     calculates numerical derivatives of P0.P2,P4,P6
c     with respect to dw
c-----------------------------------------------
c     INPUT:c     INPUT:
c     r,z,phi
c     OUTPUT: 
c     d_P_0_d_w,d_P_2_d_w,d_P_4_d_w,d_P_6_d_w
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 r,z,phi,cnr,cnz,cm
c-----output
      real*8 
     &d_P_0_d_w,d_P_2_d_w,d_P_4_d_w,d_P_6_d_w

c-----locals 
      integer i
      real*8 step,
     &hw,hfrqnc,frqncpl,df,frqncmn,
     &cnzplus,cnzminus,
     &cnrplus,cnrminus,
     &cmplus,cmminus,   
     &P_0_p,P_2_p,P_4_p,P_6_p,
     &P_0_m,P_2_m,P_4_m,P_6_m
      real*8, dimension(1:nbulk) :: vp,wp
c-----externals
      real*8 b
      step=1.d-5

      hw=step*1.0d0
      hfrqnc=hw
      hfrqnc=hfrqnc*frqncy
      bmod=b(z,r,phi)
      do i=1,nbulk
	  vp(i)=v(i)
	  wp(i)=w(i)
      enddo 

      frqncpl=frqncy+hfrqnc
      df=frqncy/frqncpl
      do i=1,nbulk
	  v(i)=vp(i)*df* df
	  w(i)=wp(i)*df
      enddo 
      cnrplus=cnr*df
      cnzplus=cnz*df
      cmplus=cm*df

      call calc_P0246(r,z,phi,cnrplus,cnzplus,cmplus,
     &                P_0_p,P_2_p,P_4_p,P_6_p)
      frqncmn=frqncy-hfrqnc
      df=frqncy/frqncmn
      do i=1,nbulk
	  v(i)=vp(i)*df*df
	  w(i)=wp(i)*df
      enddo
      cnrminus=cnr*df
      cnzminus=cnz*df
      cmminus=cm*df
      call calc_P0246(r,z,phi, cnrminus,cnzminus,cmminus,
     &P_0_m,P_2_m,P_4_m,P_6_m)

      d_P_0_d_w=(P_0_p-P_0_m)/(2.d0*hfrqnc)
      d_P_2_d_w=(P_2_p-P_2_m)/(2.d0*hfrqnc)
      d_P_4_d_w=(P_4_p-P_4_m)/(2.d0*hfrqnc)
      d_P_6_d_w=(P_6_p-P_6_m)/(2.d0*hfrqnc)

      do i=1,nbulk
	  v(i)=vp(i)
	  w(i)=wp(i)
      enddo

      return
      end


      subroutine d_P0246_dw_analitic(
     &npar,nperp,npars,   
     &d_nperp_dw,d_npar_dw, 
     &eps_perp,eps_par,eps_xy,     
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw,
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw)
c -----test analitical derivatives  d_P0246_dw 
c     INPUT:
c     &npar,nperp,npars,
c     &d_nperp_dw,d_npar_dw,
c     &eps_perp,eps_par,eps_xy, 
c     d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c     OUTPUT
c     d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw
      implicit none
c-----input
      real*8 z,r,phi,cnz,cnr,cm ,
     &npar,nperp,npars,
     &d_nperp_dw,d_npar_dw,
     &eps_perp,eps_par,eps_xy,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw

c-----output
      real*8 
     &d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw

      write(*,*)'d_P0246_dw_analitic'
      write(*,*)'npar,nperp,npars',npar,nperp,npars
      write(*,*)'&d_nperp_dw,d_npar_dw',d_nperp_dw,d_npar_dw
      write(*,*)'eps_perp,eps_par,eps_xy',eps_perp,eps_par,eps_xy
      write(*,*)'d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw',
     &           d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw

      d_P_0_dw=d_eps_par_dw*((npars-eps_perp)**2 - eps_xy**2)+
     &          eps_par*(2.d0*(npars-eps_perp)*
     &         (2.d0*npar*d_npar_dw-d_eps_perp_dw)-
     &          2.d0*eps_xy*d_eps_xy_dw)


      d_P_2_dw=(d_eps_perp_dw+d_eps_par_dw)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2.d0*npar*d_npar_dw-d_eps_perp_dw)+
     &         2.d0*eps_xy*d_eps_xy_dw


      d_P_4_dw=d_eps_perp_dw

      d_P_6_dw=0.d0

      write(*,*)'d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw',
     &           d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw

      return
      end
