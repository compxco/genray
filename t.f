      subroutine bonoli_dervs(wf,
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
 
      call dnd(z,r,phi,cnz,cnr,cm,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw)  

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

      dddw= d_P_6_dw* nperp6     +6.d0*P_6 * nperp4*nperp*d_nperp_dw+
     &      d_P_4_dw*nperp4    +4.d0*P_4 * nperps*nperp*d_nperp_dw+
     &      d_P_2_dw*nperps    +2.d0*P_2 * nperp*d_nperp_dw+
     &      d_P_0_dw

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

      d_P_2_dr=(d_eps_perp_dz+d_eps_par_dz)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2*npar*d_npar_dz-d_eps_perp_dz)+
     &                                    2.d0*eps_xy*d_eps_xy_dz

      d_P_2_dr=(d_eps_perp_dphi+d_eps_par_dphi)*(npars-eps_perp)+
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
            d_eps_perp_d_x_ar(i)    =  1.d0/y_ar(i)  
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

         d_eps_par_dr=d_eps_par_dr + d_eps_perp_d_x_ar(i)*d_x_dr_ar(i)+
     &                               d_eps_perp_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_par_dz=d_eps_par_dz + d_eps_par_d_x_ar(i)*d_x_dr_ar(i)+
     &                               d_eps_par_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_par_dphi=d_eps_par_dphi+
     &                  d_eps_par_d_x_ar(i)*d_x_dr_ar(i)+
     &                  d_eps_par_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_perp_dr=d_eps_perp_dr+d_eps_perp_d_x_ar(i)*d_x_dr_ar(i)+
     &                               d_eps_perp_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_perp_dz=d_eps_perp_dz+d_eps_perp_d_x_ar(i)*d_x_dr_ar(i)+
     &                               d_eps_perp_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_perp_dphi=d_eps_perp_dphi+
     &                   d_eps_perp_d_x_ar(i)*d_x_dr_ar(i)+
     &                   d_eps_perp_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_xy_dr=d_eps_xy_dr+d_eps_xy_d_x_ar(i)*d_x_dr_ar(i)+
     &                           d_eps_xy_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_xy_dz=d_eps_xy_dz+d_eps_xy_d_x_ar(i)*d_x_dr_ar(i)+
     &                           d_eps_xy_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_xy_dphi=d_eps_xy_dphi+
     &                 d_eps_xy_d_x_ar(i)*d_x_dr_ar(i)+
     &                 d_eps_xy_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_par_dcnr=d_eps_par_dcnr+
     &                  d_eps_par_d_x_ar(i)*d_x_dr_ar(i)+
     &                  d_eps_par_d_x_ar(i)*d_x_dr_ar(i)

         d_eps_par_dcnz=d_eps_par_dcnz+
     &                  d_eps_par_d_x_ar(i)*d_x_dr_ar(i)+
     &                  d_eps_par_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_par_dcm=d_eps_par_dcm+d_eps_par_d_x_ar(i)*d_x_dr_ar(i)+
     &                 d_eps_par_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_perp_dcnr=d_eps_perp_dcnr+
     &                   d_eps_perp_d_x_ar(i)*d_x_dr_ar(i)+
     &                   d_eps_perp_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_perp_dcnz=d_eps_perp_dcnz+
     &                   d_eps_perp_d_x_ar(i)*d_x_dr_ar(i)+
     &                   d_eps_perp_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_perp_dcm=d_eps_perp_dcm+
     &                  d_eps_perp_d_x_ar(i)*d_x_dr_ar(i)+
     &                  d_eps_perp_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_xy_dcnr=d_eps_xy_dcnr+d_eps_xy_d_x_ar(i)*d_x_dr_ar(i)+
     &                 d_eps_xy_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_xy_dcnz=d_eps_xy_dcnz+d_eps_xy_d_x_ar(i)*d_x_dr_ar(i)+
     &                 d_eps_xy_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_xy_dcm=d_eps_xy_dcm+d_eps_xy_d_x_ar(i)*d_x_dr_ar(i)+
     &                d_eps_xy_d_y_ar(i)*d_y_dr_ar(i)

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
 











