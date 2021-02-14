
      subroutine lh_ql_flux(n_parallel,n_perp,u0,umin,t,temp_kev,
     &b_ratio,e_y,e_z, 
cSAP080411
     &y_loc,
     &sin_theta0_res,cos_theta0_res,dfm_du0,
     &Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0,gamma,unorm)
c--------------------------------------------------------------------------------
c     calculate LH quasi linear flux and 
c     under integral functions for power and CD calculations
c-------------------------------------------------------------------------------
      implicit none
c-----input
      real*8
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicula refractive index

     &u0,                  !normalized momentum u/unorm 
c     &umax,                !maximal momentum normalized by unorm
c                           !unorm=dsqrt(temp_kev/t)
     &umin,                !minimal momentum at the resonace curve
                           !normalized by unorm
     &temp_kev,            !electron temperature [kev]
     &t,                   !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1.
     &b_ratio,             !B(l)/B(l=0)=B(l)/B_min 
cSAP080411
     &y_loc                !|omega_c_e/omega|
      complex*16
     &e_y,e_z              !Stix frame electric filed polarization

c-----output
      real*8 sin_theta0_res,cos_theta0_res,dfm_du0,gamma,
     &Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0, 
     &unorm  !sqrt(T/m) cm/sec
   
c-----locals
      real*8 cos_theta_res,sin_theta_res,
     &clight,u0_dev_c,
     &ksi_res,j_bessel_1,j_bessel_0,
     &D_n0_res,
     &fm,mass_e,k_1_kev,theta_relativistic,pi,
     &gamma_min

      complex*16 theta_0_big_res,i_image

      integer nz
c-----external
c     DBESJ
      

      i_image=dcmplx(0.d0,1.d0)
      clight=2.99792458d10              !light speed [cm/sec]
      pi=4.d0*datan(1.d0)

      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
      u0_dev_c=(u0*unorm)/clight
      gamma=dsqrt(1.d0+u0_dev_c**2)
      gamma_min=dsqrt(1.d0+(umin*unorm/clight)**2)
c      write(*,*)'lh_ql_flux t,temp_kev,unorm',t,temp_kev,unorm
c-------------------------------------------------------------------------------
c     pitch angle at the resonance curve at the local poloidal point
c-------------------------------------------------------------------------------
      cos_theta_res=gamma/(n_parallel*u0_dev_c)
      sin_theta_res=dsqrt(1.d0- cos_theta_res**2)
c      write(*,*)'lh_ql_flux u0,n_parallel,cos_theta_res,sin_theta_res',
c     & u0,n_parallel,cos_theta_res,sin_theta_res
c-------------------------------------------------------------------------------
c     pitch angle at the resonance curve at the midplane b=bmin
c-------------------------------------------------------------------------------
c      write(*,*)'lh_ql_flux b_ratio',b_ratio

      sin_theta0_res=sin_theta_res/dsqrt(b_ratio)
      if (cos_theta_res.ge.0d0) then
        cos_theta0_res=dsqrt(1.d0-sin_theta0_res**2)
      else
        cos_theta0_res=-dsqrt(1.d0-sin_theta0_res**2)
      endif
c      write(*,*)'lh_ql_flux ucos_theta0_res,sin_theta0_res',
c     & cos_theta0_res,sin_theta0_res
c--------------------------------------------------------------------------------
c     Bessel function agrgument
c-------------------------------------------------------------------------------
cSAP080411      ksi_res=n_perp*u0_dev_c*sin_theta_res/gamma
      ksi_res=n_perp*u0_dev_c*sin_theta_res/y_loc
c------------------------------------------------------------------------------
c      nz controles bessel function calculations
c      nz should be =0. If nz=1 the bessel function(dBESJ) will be J=0
      call DBESJ( ksi_res,0.d0,1,j_bessel_0,nz)
c      write(*,*)'lh_ql_flux after dbesj 0 nz,j_bessel_0',nz,j_bessel_0
      call DBESJ( ksi_res,1.d0,1,j_bessel_1,nz)
c      write(*,*)'lh_ql_flux after dbesj 1 nz,j_bessel',nz,j_bessel_1

cSAP080411     
c      theta_0_big_res=-i_image*e_y*j_bessel_1+
      theta_0_big_res=+i_image*e_y*j_bessel_1+
     &                j_bessel_0*cos_theta_res*e_z/sin_theta_res

cSAP101223
c      write(*,*)'lh_ql_flux e_y,e_z',e_y,e_z 
c      write(*,*)'theta_0_big_res',theta_0_big_res
c      write(*,*)'theta_0_big_res*dconjg(theta_0_big_res)',
c     &theta_0_big_res*dconjg(theta_0_big_res)

      D_n0_res=0.5d0*pi*theta_0_big_res*dconjg(theta_0_big_res)
c      write(*,*)' D_n0_res', D_n0_res
cSAP080415
c      Gamma_a_LH_u0= (gamma/u0)*D_n0_res*dabs(cos_theta_res)*
      Gamma_a_LH_u0=- (gamma/u0)*D_n0_res*dabs(cos_theta_res)*
     &                (n_parallel*u0_dev_c*sin_theta_res/gamma)**2
 
c      write(*,*)'Gamma_a_LH_u0',Gamma_a_LH_u0

      Gamma_a_LH_theta_0=-Gamma_a_LH_u0*
     &                     (sin_theta0_res/cos_theta0_res)

      dg_d_theta0=dabs((n_parallel*u0_dev_c/gamma)*
     &                  b_ratio*sin_theta0_res*cos_theta0_res/
     &                  cos_theta_res) 

      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)

      theta_relativistic=mass_e*clight**2/(k_1_kev*temp_kev)
  
cSAP070831
c      fm=dexp(-theta_relativistic*gamma)
      fm=dexp(-theta_relativistic*(gamma-gamma_min)) 
      dfm_du0=-fm*theta_relativistic*u0_dev_c/gamma
c      write(*,*)'(gamma-gamma_min),fm ',(gamma-gamma_min),fm 
      return
      end

     

      subroutine CD_adj_LH_efficiency_1(n_parallel,n_perp,
     &n_radial0,theta_pol,e_y,e_z,
     &unorm,
     &lh_cd_efficiency)
c------------------------------------------------------
c     calculate LH CD efficiency using adj chi function
c     at one radial point psis(n_radial0) from chi radial mesh 
c------------------------------------------------------
      implicit none
      include 'param.i'
      include 'adj.i'
c-----input
      integer
     & n_radial0           ! is a number of radial chi mesh

      real*8      
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicular refractive index
     &theta_pol,           !poloidal angle [radian]
     &length               !field line length normalized to  
c     &t,                   !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1. 
      complex*16 
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1
c-----output
      real*8
     &pow_dens,            !Power density normalized to m*u_norm**2*n/tau_n 
     &lh_cd_dens,          !LH current drive density normalized to q*u_norm*n
                           !Here n is the density  
     &lh_cd_efficiency,
     &unorm                !sqrt(T/m) cm/sec
c-----external
      real*8 temperho,densrho,
     &bmax_psi,bmin_psi,psi_rho,rhopsi


c-----local
      integer i,n,n_radial,idm,nmin
      real*8 du,u0,     
     &rho_small,           !small radius
     &temp_kev,            !electron temperatura [KeV]
     &dense,               !electron dencity 10**13/cm**3     
     &b_ratio,             !B/B_min
cSAP080411
     &z,r,phi,             !bellow phi is set arbitratry phi=0
     &y_loc,               !electron y at point (z,r,phi) 
     &sin_theta0_res,cos_theta0_res,dfm_du0,theta0_res,theta0_odd,
     &Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0,
     &gamma,pi,bmin,bmax,deltb,psi,th0max,sin_trap,
     &d_chi_du0,d_chi_dtheta0,u0max,chi,
     &weight_radial,
     &umin,               !minimal u0 at the resonance LH hyperbola
     &clight,              !light speed [cm/sec]
     &arg1,cln,
cfor test
     &del_power,del_cd, w1,w2,w3,w4
     
c-----external
      real*8 terp2p_Sm,b_d_bmin,y
       
      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]

      !write(*,*)'in CD_adj_LH_efficiency_1 ,e_y,e_z,',e_y,e_z

c      do n_radial=1,npsi0
c         write(*,*)'psis',psis(n_radial)
c         rho_small=rhopsi(psis(n_radial))
c         write(*,*)'rho_small',rho_small
c      enddo

      rho_small=rhopsi(psis(n_radial0))
      psi=psis(n_radial0)

      write(*,*)'n_radial0=',n_radial0

      write(*,*)'psis(n_radial0),rho_small',psis(n_radial0),rho_small

      temp_kev=temperho(rho_small,1)
      dense=densrho(rho_small,1) 
      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
      bmax= bmax_psi(psi)
      bmin= bmin_psi(psi)
	     
      deltb = (bmax-bmin)/bmin
c      write(*,*)'bmax,mbin,deltb',bmax,bmin,deltb   
      th0max = atan2(1.0d0,sqrt(deltb))
c      write(*,*)'th0max',th0max
      sin_trap=dsqrt(bmin/bmax)

c      write(*,*)'sin_trap=',sin_trap,'dsin(th0max)=',dsin(th0max)

      du = umax/nmax_chi

      pow_dens=0.d0
      lh_cd_dens=0.d0
      write(*,*)'clight,n_parallel,unorm',clight,n_parallel,unorm

cSAP111014
      if((n_parallel**2-1.d0).lt.0.d0) then   
        write(*,*)'In subroutine CD_adj_LH_efficiency_1'
        write(*,*)'n_parallel**2-1.d0.lt.0.d0'
        write(*,*)'This subroutine works for LH case' 
        write(*,*)'n_parallel**2 >1 at ieffic=5'
        write(*,*)'Please set ieffic=6 in the input file'
        write(*,*)'genray.in or genday.dat'
        stop 'in subroutine CD_adj_LH_efficiency_1'
      endif

      umin=clight/dsqrt(n_parallel**2-1.d0)/unorm ! normalized umin      
      write(*,*)'umin',umin,'umin/dsqrt(2.d0)',umin/dsqrt(2.d0)
      nmin=(umin/du+0.5d0)
      nmin=(dabs(umin)/du-0.5d0)
      nmin=nmin+1
      write(*,*)'in CD_adj_LH_efficiency_1 nmin',nmin

      b_ratio=b_d_bmin(n_radial0,theta_pol)
c      write(*,*)'n_radial0,theta_pol,b_ratio',
c     &n_radial0,theta_pol,b_ratio

      !YuP [June 2015] In some cases, b_ratio=0 (print-out: at rho_small=0)
      !Skip calc. of current drive in such case:
      if(b_ratio.le.0.d0)then ! (probably b_ratio.lt.1.d0 is better)
         !WRITE(*,*)'CD_adj_LH_efficiency_1: b_ratio=',b_ratio
         !WRITE(*,*)'CD_adj_LH_efficiency_1: rho_small=',rho_small
         !pause
           lh_cd_efficiency=0.d0
           return
      endif

cSAP080411
      phi=0.d0
      call zr_psith(psi,theta_pol,z,r)
      y_loc=y(z,r,phi,1)

      write(*,*)'in CD_adj_LH_efficiency_1  nmin,nmax_chi-1',
     &nmin,nmax_chi-1

      do n=nmin,nmax_chi-1
         u0=(n + 0.5E0)*du
        
c        write(*,*)'in CD_adj_LH_efficiency_1,e_y,e_z,',e_y,e_z

        call lh_ql_flux(n_parallel,n_perp,u0,umin,t,temp_kev,
     &                   b_ratio,e_y,e_z,
cSAP080411
     &                   y_loc,
     &                   sin_theta0_res,cos_theta0_res, dfm_du0,
     &                   Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0,
     &                   gamma,unorm)

c         write(*,*)'n,Gamma_a_LH_u0,sin_theta0_res,cos_theta0_res',
c     &              n,Gamma_a_LH_u0,sin_theta0_res,cos_theta0_res
c         write(*,*)'dfm_du0,gamma,dg_d_theta0,dfm_du0,gamma,dg_d_theta0
c     &            ',dfm_du0,gamma,dg_d_theta0,dfm_du0,gamma,dg_d_theta0

c         pow_dens=pow_dens+du*(u0**4/gamma**2)*
c     &            sin_theta0_res*cos_theta0_res*Gamma_a_LH_u0*dfm_du0/
c     &            dg_d_theta0
  
c        del_power=du*(u0**4/gamma**2)*
c     &            sin_theta0_res*cos_theta0_res*Gamma_a_LH_u0*dfm_du0/
c     &            dg_d_theta0

         pow_dens=pow_dens+du*u0**2*sin_theta0_res*
     &            (u0/gamma)*
cSAP080415 in lambda it should be used |cos_theta0_res|
c     &            (u0/gamma*cos_theta0_res)*Gamma_a_LH_u0*dfm_du0/
     &            (u0/gamma*dabs(cos_theta0_res))*Gamma_a_LH_u0*dfm_du0/
     &            dg_d_theta0
 
        del_power=du*u0**2*sin_theta0_res*
     &            (u0/gamma)*
cSAP080415 in lambda it should be used |cos_theta0_res|
c     &            (u0/gamma*cos_theta0_res)*Gamma_a_LH_u0*dfm_du0/
     &            (u0/gamma*dabs(cos_theta0_res))*Gamma_a_LH_u0*dfm_du0/
     &            dg_d_theta0

c        write(*,*)'del_power',del_power

c         write(*,*)'cos_theta0_res,Gamma_a_LH_u0,dfm_du0,dg_d_theta0',
c     &              cos_theta0_res,Gamma_a_LH_u0,dfm_du0,dg_d_theta0
c         write(*,*)'del_power=',del_power,'pow_dens=',pow_dens

         if (dabs(cos_theta0_res).le.1.d0) then
            theta0_res=dacos(cos_theta0_res)
         else
            write(*,*)'CD_adj_LH_efficiency_1 cos_theta0_res',
     &      cos_theta0_res
            if (cos_theta0_res.gt.1.d0) then
              theta0_res=0.d0
            else
              theta0_res=pi
            endif
         endif
cSAP101223
c         write(*,*)'theta0_res',theta0_res 
 
cSAP090827               
c         idm=imax_chi_a
         idm=imax_chi

c         write(*,*)'sin_theta0_res,sin_trap',sin_theta0_res,sin_trap
c         write(*,*)'dsin(th0a_chi(imax_chi,n_radial0))',
c     &              dsin(th0a_chi(imax_chi,n_radial0))

         if((theta0_res.le.th0a_chi(imax_chi,n_radial0)).or.
     &     (theta0_res.ge.pi-th0a_chi(imax_chi,n_radial0))) then
c         if (sin_theta0_res.lt.sin_trap) then
c-----------passing particles
c            write(*,*)'theta0_res,u0',theta0_res,u0
            if (theta0_res.le.th0a_chi(imax_chi,n_radial0)) then
cSAP090129
               if (i_chi_interpolation.eq.2) then
                 call  interpolation_chi(2,n_radial0,u0,theta0_res,
     &           chi,d_chi_du0,d_chi_dtheta0)
               endif
cSAP090129
               if (i_chi_interpolation.eq.1) then
c-----------------spline interpolation
                 d_chi_du0=terp2p_Sm(theta0_res,u0,
     &           imax_chi,th0a_chi(1,n_radial0),
     &           nmax_chi,ua,
     &           chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &           chi_uu(1,1,n_radial0),
     &           chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)
            
                 d_chi_dtheta0=terp2p_Sm(theta0_res,u0,
     &           imax_chi,th0a_chi(1,n_radial0),
     &           nmax_chi,ua,
     &           chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &           chi_uu(1,1,n_radial0),
     &           chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)
               endif
c               write(*,*)'d_chi_du0,d_chi_dtheta0',
c     &                    d_chi_du0,d_chi_dtheta0
            else ! theta0_res.ge.pi-th0a_chi(imax_chi,n_radial0)
            
               theta0_odd=pi-theta0_res 
cSAP090129
               if (i_chi_interpolation.eq.2) then
                 call  interpolation_chi(2,n_radial0,u0,theta0_odd,
     &           chi,d_chi_du0,d_chi_dtheta0)
c              write(*,*)'interpolation d_chi_du0,d_chi_dtheta0',
c     &                    d_chi_du0,d_chi_dtheta0
                 d_chi_du0=-d_chi_du0
               endif
cSAP090129
               if (i_chi_interpolation.eq.1) then
c-----------------spline interpolation
                  d_chi_du0=-terp2p_Sm(theta0_odd,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)
            
                  d_chi_dtheta0=terp2p_Sm(theta0_odd,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)
               endif
               
            endif ! theta0_res.ge.pi-th0a_chi(imax_chi,n_radial0)
            
c                          write(*,*)'d_chi_du0,d_chi_dtheta0',
c     &                    d_chi_du0,d_chi_dtheta0

c            lh_cd_dens=lh_cd_dens+du*(u0**3/gamma)*
c     &                 cos_theta0_res*Gamma_a_LH_u0*
c     &                 (d_chi_du0-d_chi_dtheta0/(u0*cos_theta0_res))*
c     &                 dfm_du0/dg_d_theta0

c            del_cd=du*(u0**3/gamma)*
c     &                 cos_theta0_res*Gamma_a_LH_u0*
c     &                 (d_chi_du0-d_chi_dtheta0*sin_theta0_res/
c     &                 (u0*cos_theta0_res))*
c     &                 dfm_du0/dg_d_theta0

            lh_cd_dens=lh_cd_dens+du*u0**2*sin_theta0_res*
cSAP080415 in lambda it should be used |cos_theta0_res|
c     &           (u0/gamma*cos_theta0_res)*Gamma_a_LH_u0*
     &           (u0/gamma*dabs(cos_theta0_res))*Gamma_a_LH_u0*
     &           (d_chi_du0-
     &           (d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))*
     &                 dfm_du0/dg_d_theta0
c      if(abs(rho_small-0.444).le.0.001)then
       !w1=theta0_res
       !w2=th0a_chi(imax_chi,n_radial0)
       !w3=pi-th0a_chi(imax_chi,n_radial0)
       !w4=d_chi_du0
       !!w=(d_chi_du0-(d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))
c      WRITE(*,'(a,i5,5e11.3)')'in CD_adj_LH_efficiency_1',
c     &  n, rho_small, temp_kev, du, n_parallel, umin
c      !pause
c      endif


            del_cd=du*u0**2*sin_theta0_res*
cSAP080415 in lambda it should be used |cos_theta0_res|
c     &           (u0/gamma*cos_theta0_res)*Gamma_a_LH_u0*
     &           (u0/gamma*dabs(cos_theta0_res))*Gamma_a_LH_u0*
     &           (d_chi_du0-
     &           (d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))*
     &                 dfm_du0/dg_d_theta0

c            write(*,*)'del_cd=',del_cd
          
            
c            write(*,*)'(d_chi_du0-'//
c     & '(d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))',
c     & (d_chi_du0-(d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))   
             
c            write(*,*)'pow_dens,lh_cd_dens=',pow_dens,lh_cd_dens        

         endif 
      enddo ! n=nmin,nmax_chi-1

      pow_dens=2.d0*pi*pow_dens
      lh_cd_dens=-2.d0*pi*LH_CD_dens
c      write(*,*)'pow_dens',pow_dens
c      write(*,*)'lh_cd_dens',lh_cd_dens
     
      if(pow_dens.ne.0.d0) then
        lh_cd_efficiency=lh_cd_dens/pow_dens
      else
        lh_cd_efficiency=0.d0
      endif
               
c-----cln -Coulomb Ln
      arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense)	!temp KeV, dense 10**13/cm**3
      cln=24.d0-dlog(arg1)
      lh_cd_efficiency=lh_cd_efficiency*
     &3.84d0*temp_kev/(cln*dense)	       !  (A/m**2)/(W/m**3)
c      write(*,*)'temp_kev,dense,cln,lh_cd_efficiency',
c     &temp_kev,dense,cln,lh_cd_efficiency
      lh_cd_efficiency=lh_cd_efficiency*1.d-5  !  A/cm**2/(erg/(sec*cm**3)
      write(*,*)'temp_kev,dense,cln,lh_cd_efficiency',
     &temp_kev,dense,cln,lh_cd_efficiency
      return
      end

      subroutine splcoef_chi_1
c-------------------------------------------------------------
c     calculate spline coefficients for chi 2D function
c     for all radial points
c--------------------------------------------------------------

      implicit none
      include 'param.i'
      include 'adj.i'
c-----input
      integer
c     & imax_chi,   !number of points in momentum mesh
c     & nmax_chi,   !number of points in pitch angle mesh
     & n_radial    !number of radial point  

c-----local  
cSAP090827
      integer chi_nwka,istat  
c     in general,need chi_nwka=1+3*max(nmax_chi_a,imax_chi_a)
c      parameter (chi_nwka=1+3*nmax_chi_a)   
c     parameter (chi_nwka=1+3*imax_chi_a)   

c      real*8 
c     &chi_3d(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a), !3D chi function
c     &psis(1:npsi0_a),dene_chi(1:npsi0_a),teme_chi(1:npsi0_a), 
c     &bmin_chi(npsi0_a),bmax_chi(npsi0_a),
c     &umax_chi(1:npsi0_a),
c     &th0max_chi(npsi0_a),th0a_chi(1:imax_chi_a,1:npsi0_a),
c     &wk_chi(chi_nwka)
c-----second derivativers
c     &chi_tt(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a), !d^2chi/dpitch
c     &chi_uu(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a), ! d^2f/du
c-----fourth derivatives                         
c      chi_ttuu(1:imax_chi_a,1:nmax+chi_a,1:npsi0_a), !d^4chi/(d/dpitch)**2*(d/dx)**2

      real*8, pointer :: wk_chi(:) !1:(chi_nwka=1+3*max(nmax_chi_a,imax_chi_a))

      integer ibd(4),idm
      integer i,n


cfor test
       real*8 u0,theta0,chi_spline,norm,pi,theta0_odd
c-----externals
       real*8 terp2p_Sm
cendtest
c------------------------------------------------------------------
c     check parameter chi_nwka value 
c------------------------------------------------------------------
c      if (chi_nwka.lt.(1+3*nmax_chi_a)) then
c         write(*,*)'**************************************' 
c         write(*,*)'lh_ql_flux.f in  splcoef_chi'
c         write(*,*)'chi_nwka.lt.(1+3*nmax_chi_a)'
c         write(*,*)'it shoud be chi_nwka=1+3*max(nmax_chi_a,imax_chi_a)'
c         write(*,*)'chi_nwka, nmax_chi_a,imax_chi_a',
c     &              chi_nwka, nmax_chi_a,imax_chi_a
c         write(*,*)'Please change parameter in csplcoef_chi'
c         write(*,*)'for parameter (chi_nwka=1+3*nmax_chi_a)'
c         write(*,*)'and recompile code'  
c         write(*,*)'**************************************' 
c         stop
c      endif

c      if (chi_nwka.lt.(1+3*imax_chi_a)) then
c         write(*,*)'**************************************' 
c         write(*,*)'lh_ql_flux.f in  splcoef_chi'
c         write(*,*)'chi_nwka.lt.(1+3*imax_chi_a)'
c         write(*,*)'it shoud be chi_nwka=1+3*max(nmax_chi_a,imax_chi_a)'
c         write(*,*)'chi_nwka, nmax_chi_a,imax_chi_a',
c     &              chi_nwka, nmax_chi_a,imax_chi_a
c         write(*,*)'Please change parameter in csplcoef_chi'
c         write(*,*)'for parameter (chi_nwka=1+3*imax_chi_a)'
c         write(*,*)'and recompile code'  
c         write(*,*)'**************************************' 
c         stop
c      endif
 
 
c      chi_3d(1:imax_chi,1:nmax_chi,np)=chi(0:imax_chi-1,0:nmax_chi-1)
c      psis(n_radial)=psis
c      dene_chi(1:npsi0a)=dene
c      teme_chi(1:npsi0a)=teme
c      bmin_chi(n_radial)=bmin
c      bmax_chi(n_radial)=bmax
c      th0max_chi(n_radial)=th0max
c      th0a_chi(1:imax_chi_a,n_radial)=th0a(0:imax_chi_a-1)

c---------------------------------------------------
c     allocate pointer wk_chi
c--------------------------------------------------
      chi_nwka=1+3*max0(nmax_chi,imax_chi)
      allocate( wk_chi(1:chi_nwka),STAT=istat)
      call bcast(wk_chi,0.d0,SIZE(wk_chi))
 
c-----creates 2D spline coefficient for chi function
      ibd(1)=2
      ibd(2)=2
      ibd(3)=4
      ibd(4)=4

cSAP090827
c      idm=imax_chi_a
      idm=imax_chi

c      write(*,*)'splcoef_chi umax,nmax_chi',umax,nmax_chi
c      write(*,*)' splcoef_chi ua',ua
      do n=1,nmax_chi
        ua(n)=(n-0.5d0)*umax/nmax_chi
c        write(*,*)'n,ua(n)',n,ua(n)
      enddo
c      write(*,*)' splcoef_chi 1 ua',ua

      do n_radial=1,npsi0
        do n=1,nmax_chi
          chi_tt(1,n,n_radial)=0.d0
          chi_tt(imax_chi,n,n_radial)=0.d0
        enddo
c        write(*,*)' splcoef_ch n_radial',n_radial
c        write(*,*)' th0a_chi(1,n_radial)'
c        write(*,*) (th0a_chi(i,n_radial),i=1,imax_chi)
  
        call coeff2_Sm(imax_chi,th0a_chi(1,n_radial),nmax_chi,ua,
     &   chi_3d(1,1,n_radial),
     &   chi_tt(1,1,n_radial),chi_uu(1,1,n_radial),
     &   chi_ttuu(1,1,n_radial),
cSAP090827
c     &   idm,ibd,wk_chi,imax_chi_a)
     &   idm,ibd,wk_chi,imax_chi)
      enddo

c---------------------------------------------------------
      deallocate (wk_chi,STAT=istat)
c-----------------------------------------------------------

ctest spline
      idm=imax_chi
      norm=0.d0 
      do n_radial=1,npsi0
        do n=1,nmax_chi
           u0=ua(n)  
           do i=1,imax_chi
              theta0=th0a_chi(i,n_radial)
              chi_spline=terp2p_Sm(theta0,u0,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP090827
c     &        chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi)
cSAP101223
c              write(*,*)'n_radial,n,i,u0,theta0',n_radial,n,i,u0,theta0
c              write(*,*)'chi_3d(i,n,n_radial),chi_spline',
c     &                   chi_3d(i,n,n_radial),chi_spline
              norm=norm+dabs(chi_3d(i,n,n_radial)-chi_spline)
cc-------------dchi_d_theta
              chi_spline=terp2p_Sm(theta0,u0,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP090827
c     &        chi_ttuu(1,1,n_radial),idm,1,0,1,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,1,0,1,imax_chi)
cSAP101223
c              write(*,*)'spline dchi_dtheta=',  chi_spline
              if ((i.ne.1).and.(i.ne.imax_chi)) then
              chi_spline=(chi_3d(i+1,n,n_radial)-chi_3d(i-1,n,n_radial))
     &                 /(th0a_chi(i+1,n_radial)-th0a_chi(i-1,n_radial))
              else
                if (i.eq.1) then
                   chi_spline=(-1.5d0*chi_3d(1,n,n_radial)
     &                       + 2.d0*chi_3d(2,n,n_radial)
     &                       -0.5d0*chi_3d(3,n,n_radial))/
     &                 (th0a_chi(2,n_radial)-th0a_chi(1,n_radial))
                else
                   chi_spline=(1.5d0*chi_3d(imax_chi,n,n_radial)
     &                       - 2.d0*chi_3d(imax_chi-1,n,n_radial)
     &                       +0.5d0*chi_3d(imax_chi-2,n,n_radial))/
     &                 (th0a_chi(2,n_radial)-th0a_chi(1,n_radial))
                endif
              endif
              write(*,*)'difference dchi_dtheta=',  chi_spline
cc-------------dchi_d_u
              chi_spline=terp2p_Sm(theta0,u0,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP090827 
c    &        chi_ttuu(1,1,n_radial),idm,0,1,1,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,0,1,1,imax_chi)

              write(*,*)'spline dchi_du=',  chi_spline

              if ((n.ne.1).and.(n.ne.nmax_chi)) then
              chi_spline=(chi_3d(i,n+1,n_radial)-chi_3d(i,n-1,n_radial))
     &                 /(ua(n+1)-ua(n-1))
              else
                if (n.eq.1) then
                   chi_spline=(-1.5d0*chi_3d(i,1,n_radial)
     &                       + 2.d0*chi_3d(i,2,n_radial)
     &                       -0.5d0*chi_3d(i,3,n_radial))/
     &                        (ua(2)-ua(1))
                else
                   chi_spline=(1.5d0*chi_3d(i,nmax_chi,n_radial)
     &                       - 2.d0*chi_3d(i,nmax_chi-1,n_radial)
     &                       +0.5d0*chi_3d(i,nmax_chi-2,n_radial))/
     &                       (ua(2)-ua(1))
                endif
              endif
              write(*,*)'difference dchi_du=',  chi_spline
cc-------------dchi_d_u
           enddo
        enddo
      enddo 

      write(*,*)'norm=',norm

      pi=4.d0*datan(1.d0)
      n_radial= 8
      n_radial= 1
c      theta0=3.079727452065930 
      theta0=3.075248747886528d0
      theta0=6.675884388880000d-02
c      u0=9.670000000000000d0
      u0=4.235000000000000d0

      if((theta0.gt.th0a_chi(imax_chi,n_radial)).and.
     &   (theta0.lt.pi-th0a_chi(imax_chi,n_radial))) then
         chi_spline=0.d0
         write(*,*)'theta in trapped area'
         write(*,*)'chi and its derivatives  =0'
      else      
         if(theta0.le.th0a_chi(imax_chi,n_radial)) then
            chi_spline=terp2p_Sm(theta0,u0,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP090827
c     &        chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi)
cSAP101223
c            write(*,*)'n_radial,n,i,u0,theta0',n_radial,n,i,u0,theta0
c            write(*,*)'chi_spline=',chi_spline

            chi_spline=terp2p_Sm(theta0,u0,
     &          imax_chi,th0a_chi(1,n_radial),
     &          nmax_chi,ua,
     &          chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &          chi_uu(1,1,n_radial),
cSAP090827
c     &          chi_ttuu(1,1,n_radial),idm,0,1,1,imax_chi_a)
     &          chi_ttuu(1,1,n_radial),idm,0,1,1,imax_chi)

cSAP101223
c            write(*,*)'check n_radial,theta0,u0=',n_radial,theta0,u0
c            write(*,*)'spline dchi_du=',  chi_spline
            chi_spline=terp2p_Sm(theta0,u0,
     &          imax_chi,th0a_chi(1,n_radial),
     &          nmax_chi,ua,
     &          chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &          chi_uu(1,1,n_radial),
cSAP090827
c     &          chi_ttuu(1,1,n_radial),idm,1,0,1,imax_chi_a)
     &          chi_ttuu(1,1,n_radial),idm,1,0,1,imax_chi)


            write(*,*)'spline dchi_dtheta=',  chi_spline
            chi_spline=terp2p_Sm(theta0,u0,
     &          imax_chi,th0a_chi(1,n_radial),
     &          nmax_chi,ua,
     &          chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &          chi_uu(1,1,n_radial),
cSAP090827
c     &          chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi_a)
     &          chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi)

            write(*,*)'spline chi=',  chi_spline
         else
c-----------theta0.ge.pi-th0a_chi(imax_chi,n_radial)
            theta0_odd=pi-theta0
            write(*,*)'theta0_odd',theta0_odd
            chi_spline=-terp2p_Sm(theta0_odd,u0,
     &          imax_chi,th0a_chi(1,n_radial),
     &          nmax_chi,ua,
     &          chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &          chi_uu(1,1,n_radial),
cSAP090827
c     &          chi_ttuu(1,1,n_radial),idm,0,1,1,imax_chi_a)
     &          chi_ttuu(1,1,n_radial),idm,0,1,1,imax_chi)

            write(*,*)'check n_radial,theta0,u0=',n_radial,theta0,u0
            write(*,*)'spline dchi_du=',  chi_spline

            chi_spline=terp2p_Sm(theta0_odd,u0,
     &          imax_chi,th0a_chi(1,n_radial),
     &          nmax_chi,ua,
     &          chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &          chi_uu(1,1,n_radial),
cSAP090827
c     &          chi_ttuu(1,1,n_radial),idm,1,0,1,imax_chi_a)
     &          chi_ttuu(1,1,n_radial),idm,1,0,1,imax_chi)

            write(*,*)'spline dchi_dtheta=',  chi_spline
            chi_spline=-terp2p_Sm(theta0_odd,u0,
     &          imax_chi,th0a_chi(1,n_radial),
     &          nmax_chi,ua,
     &          chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &          chi_uu(1,1,n_radial),
cSAP090827
c     &          chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi_a)
     &          chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi)


            write(*,*)'spline chi=',  chi_spline
         endif
      endif

      return
      end



      subroutine read_iout5_chi
c-------------------------------------------------------------------
c     reads output file iout5 with calculated chi function
c     and some additional arrays
c-------------------------------------------------------------------
      use kind_spec

      implicit none

c-----input      
      include 'param.i'
      include 'adj.i'

c-----local
      integer n_radial,nmax,imax,kmax,lmax 

      integer i,n
        
      real(kind=dp), dimension(0:nmax_chi - 1) :: ua1  
      real(kind=dp), dimension(0:nmax_chi) :: ug
      real(kind=dp), dimension(0:imax_chi - 1) :: th0a1
      real(kind=dp), dimension(0:imax_chi) :: th0g
      real(kind=dp), dimension(0:nmax_chi - 1) :: fm
      real(kind=dp), dimension(0:imax_chi-1,0:nmax_chi-1) :: chi      
C-----------------------------------------------
   
      real(kind=dp)  psimx,psimn,cond,abserr,relerr,c2,zi,eps
      real(kind=dp)  th0max,du,dth0,deltb

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
      open(iout5, file='adjout', status='old') ! in read_iout5_chi
c      write(*,*)'read_iout5_chi iout5=',iout5
      !endif  !On myrank=0   ! myrank=0
 
 2000 format(1x,2i8)
 2010 format(4(2x,e20.12))
 2020 format(5(2x,e20.12))
 3000 format(4(2x,e20.12))
 3010 format(1x,4i8)
 3020 format(4(1x,e20.12)) 
 3030 format(5(1x,e20.12))
 3040 format(1x,e20.12,2x,i8,2x,e20.12)

      !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
      read (iout5, 2000) nthp0, npsi0 ! in read_iout5_chi
c      write(*,*)' nthp0, npsi0', nthp0, npsi0
      read (iout5, 2010) psimx, psimn
c      write(*,*)' psimx, psimn', psimx, psimn
      !endif  !On myrank=0   ! myrank=0

      do n_radial=1,npsi0 

         !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         write(*,*)'lh_ql_flux.f read_iout5_chi n_radial=', n_radial
         read (iout5, 3000) psis(n_radial),dene_chi(n_radial),
     &                      teme_chi(n_radial)
         write(*,*)'psis(n_radial),dene_chi(n_radial),teme_chi(n_radial)
     &   ',psis(n_radial),dene_chi(n_radial),teme_chi(n_radial)
         read (iout5, 3010) nmax, imax, kmax, lmax
c         write(*,*)' nmax, imax, kmax, lmax', nmax, imax, kmax, lmax
         read (iout5,3020) t, c2, ze, zi, eps
c         write(*,*)'t, c2, ze, zi, eps',t, c2, ze, zi, eps
         !endif  !On myrank=0   ! myrank=0

         zi_chi(n_radial)=zi ! MPI_BCAST -ed

         !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         read (iout5,3020)umax,th0max,du,dth0,rho_,deltb,
     &   bmin_chi(n_radial),bmax_chi(n_radial)
         read (iout5, 3040) dt, tmax, alpha
c         write(*,*)' dt, tmax, alpha', dt, tmax, alpha
         read (iout5, 3020) cond, abserr, relerr
c         write(*,*)'cond, abserr, relerr',cond, abserr, relerr
c         write(*,*)'before read ua iout5=',iout5
         read (iout5, 3030) ua1           !ua
c         write(*,*)'after read ua iout5=',iout5
c         write (*, 3030) ua1 
c         write(*,*)'before read ug iout5=',iout5
         read (iout5, 3030) ug
c         write(*,*)'ug'
c         write(*,3030) ug
c         write(*,*)'before read th0a iout5=',iout5
         read (iout5, 3030) th0a1             !th0a
c         write(*,*)'th0a'
c         write(*, 3030) th0a1   
         th0a_chi(1:imax_chi,n_radial)= th0a1(0:imax_chi-1) ! MPI_BCAST -ed
         !endif  !On myrank=0   ! myrank=0
	 
         !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
c         write(*,*)'before read th0g iout5=',iout5
         read (iout5, 3030) th0g                                    !th0g
c         write(*,*)'th0g'
c         write (*, 3030) th0g
c         write(*,*)'before read lam1a iout5=',iout5
         read (iout5, 3030) th0a1             !lam1a
c         write(*,*)'lam1a'
c         write(*, 3030) th0a1 
c         write(*,*)'before read fm iout5=',iout5
         read (iout5, 3030) fm             !fm         
c         write(*,*)'fm'
c         write(*,3030) fm
c         write(*,*)'before read chi iout5=',iout5
         read (iout5, 3030) ((chi(i,n),i=0,imax_chi-1),n=0,nmax_chi-1)
c         write(*,*)'iout5 chi'
c         write(*,3030) ((chi(i,n),i=0,imax_chi-1),n=0,nmax_chi-1)
         !endif  !On myrank=0   ! myrank=0

         do i=1,imax_chi
            do n=1,nmax_chi
               chi_3d(i,n,n_radial)=chi(i-1,n-1) ! in read_iout5_chi ! MPI_BCAST -ed
c               write(*,*)'i,n,n_radial,chi_3d(i,n,n_radial)',
c     &                    i,n,n_radial,chi_3d(i,n,n_radial)
            enddo
         enddo

c         chi_3d(1:imax_chi,1:nmax_chi,n_radial)=
c     &   chi(0:imax_chi-1,0:nmax_chi-1)

      enddo

      do n=1,nmax_chi
         ua(n)=ua1(n-1) ! MPI_BCAST -ed
      enddo
     
      !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
      close (iout5) ! in read_iout5_chi
      !endif  !On myrank=0   ! myrank=0
      return
      end


      subroutine find_psis_points(psi_loc,
     & n_radial_minus,n_radial_plus,weight_radial)
c------------------------------------------------------------
c     find the points of chi radial mesh
c      
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input
      real*8 psi_loc           !poloidal flux 
c-----output
      integer n_radial_minus,n_radial_plus
      real*8 weight_radial
c-----local
      integer n_radial

c      write(*,*)'find_psis_points psi_loc, psi_loc',psi_loc
c      write(*,*)'psis',psis

      if (psi_loc.le.psis(1)) then
         n_radial_plus=1
         n_radial_minus=1
         goto 10
      else
         if (psi_loc.ge.psis(npsi0)) then
            n_radial_plus=npsi0
            n_radial_minus=npsi0
            goto 10
         else
            do n_radial=2,npsi0
              if (psi_loc.lt.psis(n_radial)) then
                  n_radial_plus=n_radial
                  n_radial_minus=n_radial-1
                  goto 10             
              endif
            enddo
         endif
      endif  
     
 10   continue

c      write(*,*)'psi_loc',psi_loc
c      write(*,*)'npsi0',npsi0
c      write(*,*)'psis(n_radial)',(psis(n_radial),n_radial=1,npsi0)
c      write(*,*)' n_radial_plus n_radial_minus',
c     & n_radial_plus,n_radial_minus

      if((n_radial_plus.gt.1).and.(n_radial_minus.lt.npsi0)) then
         weight_radial=(psi_loc-psis(n_radial_minus))/
     &              (psis(n_radial_plus)-psis(n_radial_minus))      
      else
         if (n_radial_plus.eq.1) then
             weight_radial=0.d0
         else
             weight_radial=1.d0
         endif
      endif 
  
c      write(*,*)'weight_radial',weight_radial

      return
      end

      subroutine CD_adj_LH_efficiency(n_parallel,n_perp,
     &psi_loc,theta_pol,e_y,e_z,
     &lh_cd_efficiency)
c-----------------------------------------------------------------
c     calculate LH CD efficiency using adj chi function
c-----------------------------------------------------------------       
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input
      real*8      
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicular refractive index
     &psi_loc,             !polidal flux 
     &theta_pol,           !poloidal angle
     &length               !field line length normalized to  
c     &t,                  !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1. 
      complex*16  
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1
c-----output
      real*8 
     &lh_cd_efficiency        
c-----local
      integer  n_radial_minus,n_radial_plus
      real*8 weight_radial,lh_cd_efficiency_minus,
     &lh_cd_efficiency_plus,unorm_minus,unorm_plus,
     &v_ph
 
c      write(*,*)' CD_adj_LH_efficiency e_y,e_z', e_y,e_z
c      write(*,*)' CD_adj_LH_efficiency before find_psis_points'
c      write(*,*)' n_radial_minus,n_radial_plus,weight_radial',
c     &             n_radial_minus,n_radial_plus,weight_radial
      v_ph=1.d0 !to test analytical LH flux CD_adj_LH_efficiency_anal_1
 
      call  find_psis_points(psi_loc,
     & n_radial_minus,n_radial_plus,weight_radial)

c      WRITE(*,'(a,2e11.3,2i5)')
c     + ' CD_adj_LH_efficiency after find_psis_points',
c     + psi_loc, n_parallel, n_radial_minus,n_radial_plus
      !pause

      if((n_radial_minus.lt.npsi0).and.(n_radial_plus.gt.1)) then

         call CD_adj_LH_efficiency_1(n_parallel,n_perp,
     &   n_radial_minus,theta_pol,e_y,e_z,
     &   unorm_minus,
     &   lh_cd_efficiency_minus)

         write(*,*)'eff lh_cd_efficiency_minus',lh_cd_efficiency_minus

c         call CD_adj_LH_efficiency_anal_1(v_ph,n_parallel,
c     &   n_radial_minus,theta_pol,
c     &   unorm_minus,
c     &   lh_cd_efficiency_minus)
c         write(*,*)'anal lh_cd_efficiency_minus',lh_cd_efficiency_minus

         call CD_adj_LH_efficiency_1(n_parallel,n_perp,
     &   n_radial_plus,theta_pol,e_y,e_z,
     &   unorm_plus,
     &   lh_cd_efficiency_plus)
         write(*,*)'eff lh_cd_efficiency_plus',lh_cd_efficiency_plus

c         call CD_adj_LH_efficiency_anal_1(v_ph,n_parallel,
c     &   n_radial_plus,theta_pol,
c     &   unorm_plus,
c     &   lh_cd_efficiency_plus)
c          write(*,*)'anal lh_cd_efficiency_plus',lh_cd_efficiency_plus

         lh_cd_efficiency=lh_cd_efficiency_plus*weight_radial+
     &   lh_cd_efficiency_minus*(1.d0-weight_radial)

      else
        if (n_radial_plus.eq.1) then  
           call CD_adj_LH_efficiency_1(n_parallel,n_perp,
     &     1,theta_pol,e_y,e_z,
     &     unorm_plus,
     &     lh_cd_efficiency)
            write(*,*)'eff lh_cd_efficiency 1',lh_cd_efficiency

c           call CD_adj_LH_efficiency_anal_1(v_ph,n_parallel,
c     &     1,theta_pol,
c     &      unorm_plus,
c     &     lh_cd_efficiency)
c           write(*,*)'anal lh_cd_efficiency 1',lh_cd_efficiency
        else
           if (n_radial_minus.eq.npsi0) then
              call CD_adj_LH_efficiency_1(n_parallel,n_perp,
     &        npsi0,theta_pol,e_y,e_z,
     &        unorm_minus,
     &        lh_cd_efficiency)
           write(*,*)'eff lh_cd_efficiency npsi0',lh_cd_efficiency
c           call CD_adj_LH_efficiency_anal_1(v_ph,n_parallel,
c     &     npsi0,theta_pol,
c     &     unorm_plus,
c     &     lh_cd_efficiency)
c           write(*,*)'anal lh_cd_efficiency npsi0',lh_cd_efficiency

           endif
        endif
      endif

      write(*,*)'in CD_adj_LH_efficiency lh_cd_efficiency=',
     &lh_cd_efficiency

      return
      end

      real*8 function b_d_bmin(n_radial,theta_poloidal)
c------------------------------------------------------------
c     calculate B/B_min at the given n_radial point and 
c     poloidal angle
c-------------------------------------------------------------      
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input
      integer n_radial        !is the number of the flux point at chi flux mesh
 
      real*8 theta_poloidal    !poloidal angle radian  
                               !0 =< theta_poloidal =< 2pi  
  
c-----local
      integer nth
      real*8 pi

      pi=4.d0*datan(1.d0)
    
c      write(*,*)'lh_ql_flux.f in b_d_bmin n_radial,theta_poloidal',
c     &n_radial,theta_poloidal
    

      if (theta_poloidal.lt.0.d0) then
         write(*,*)'lh_ql_flux.f in function b_d_bmin'
         write(*,*)'theta_poloidal.lt.0'
         write(*,*)'theta_poloidal=',theta_poloidal
         write(*,*)'it should be 0 =< theta_poloidal =< 2pi'
         stop 'b_d_bmin'
      endif

      if (theta_poloidal.gt.2.d0*pi) then
         write(*,*)'lh_ql_flux.f in function b_d_bmin'
         write(*,*)'theta_poloidal.gt.2*pi'
         write(*,*)'theta_poloidal=',theta_poloidal
         write(*,*)'it should be 0 =< theta_poloidal =< 2pi'
         stop 'b_d_bmin'
      endif

      if (theta_poloidal.le.thetae_1(0)) then 
         b_d_bmin=ba_2d(0,n_radial)     
      else
cSm070921
c         if (theta_poloidal.ge.thetae_1(0)) then 

c         write(*,*)'in b_d_bmin theta_poloidal,thetae_1(nthp0-1)',
c     &              theta_poloidal,thetae_1(nthp0-1)

         if (theta_poloidal.ge.thetae_1(nthp0-1)) then 
            b_d_bmin=ba_2d(nthp0-1,n_radial)   
         else
            do nth=1,nthp0-1

c               write(*,*)'nth,theta_poloidal,thetae_1(nth)',
c     &                   nth,theta_poloidal,thetae_1(nth)

               if (theta_poloidal.le.thetae_1(nth)) then
                  b_d_bmin=ba_2d(nth-1,n_radial)+
     &            (theta_poloidal-thetae_1(nth-1))*
     &            (ba_2d(nth,n_radial)-ba_2d(nth-1,n_radial))/
     &            (thetae_1(nth)-thetae_1(nth-1))

c       write(*,*)'ba_2d(nth-1,n_radial)',ba_2d(nth-1,n_radial)
c       write(*,*)'ba_2d(nth,n_radial)',ba_2d(nth,n_radial)
c       write(*,*)'(theta_poloidal-thetae_1(nth-1))',
c     &            (theta_poloidal-thetae_1(nth-1))
c       write(*,*)'(ba_2d(nth,n_radial)-ba_2d(nth-1,n_radial))',
c     &            (ba_2d(nth,n_radial)-ba_2d(nth-1,n_radial))
c       write(*,*)' (thetae_1(nth)-thetae_1(nth-1))',
c     &             (thetae_1(nth)-thetae_1(nth-1))
c       write(*,*)'b_d_bmin',b_d_bmin
cSm070921
                  goto 10
               endif
            enddo
         endif
      endif 
cSm070921
 10   continue
 
c     write(*,*)'lh_ql_flux.f in b_d_bmin=',b_d_bmin

      return
      end


      subroutine splcoef_chi
c-------------------------------------------------------------
c     calculate spline coefficients for chi 2D function
c     for all radial points
c--------------------------------------------------------------

      implicit none
      include 'param.i'
      include 'adj.i'
c-----input
      integer
c     & imax_chi,   !number of points in momentum mesh
c     & nmax_chi,   !number of points in pitch angle mesh
     & n_radial    !number of radial point  

c-----local  
cSAP090827
      integer chi_nwka,istat  
c     in general,need chi_nwka=1+3*max(nmax_chi_a,imax_chi_a)
c      parameter (chi_nwka=1+3*nmax_chi_a)   
c     parameter (chi_nwka=1+3*imax_chi_a)   

c      real*8 
c     &chi_3d(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a), !3D chi function
c     &psis(1:npsi0_a),dene_chi(1:npsi0_a),teme_chi(1:npsi0_a), 
c     &bmin_chi(npsi0_a),bmax_chi(npsi0_a),
c     &umax_chi(1:npsi0_a),
c     &th0max_chi(npsi0_a),th0a_chi(1:imax_chi_a,1:npsi0_a),
c     &wk_chi(chi_nwka)
c-----second derivativers
c     &chi_tt(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a), !d^2chi/dpitch
c     &chi_uu(1:imax_chi_a,1:nmax_chi_a,1:npsi0_a), ! d^2f/du
c-----fourth derivatives                         
c      chi_ttuu(1:imax_chi_a,1:nmax+chi_a,1:npsi0_a), !d^4chi/(d/dpitch)**2*(d/dx)**2


      real*8, pointer :: wk_chi(:) !1:(chi_nwka=1+3*max(nmax_chi_a,imax_chi_a))

      integer ibd(4),idm
      integer i,n

cfor test 
      real*8 u0,theta0,chi_spline,norm,
     & norm_dchi_du,norm_dchi_dtheta,dchi_dtheta_spline,dchi_du_spline,
     & chi,d_chi_d_u_out,d_chi_d_theta_out,hu,htheta
c-----externals
       real*8 terp2p,terp2p_Sm
cendtest
c      write(*,*)'splcoef_chi begin'
c------------------------------------------------------------------
c     check parameter chi_nwka value 
c------------------------------------------------------------------
c      if (chi_nwka.lt.nmax_chi_a) then
c         write(*,*)'**************************************' 
c         write(*,*)'lh_ql_flux.f in  splcoef_chi'
c         write(*,*)'chi_nwka.lt.nmax_chi_a'
c         write(*,*)'it shoud be chi_nwka=1+3*max(nmax_chi_a,imax_chi_a)'
c         write(*,*)'chi_nwka, nmax_chi_a,imax_chi_a',
c     &              chi_nwka, nmax_chi_a,imax_chi_a
c         write(*,*)'Please change parameter in csplcoef_chi'
c         write(*,*)'for parameter (chi_nwka=1+3*nmax_chi_a)'
c         write(*,*)'and recompile code'  
c         write(*,*)'**************************************' 
c         stop
c      endif

c      if (chi_nwka.lt.imax_chi_a) then
c         write(*,*)'**************************************' 
c         write(*,*)'lh_ql_flux.f in  splcoef_chi'
c         write(*,*)'chi_nwka.lt.imax_chi_a'
c         write(*,*)'it shoud be chi_nwka=1+3*max(nmax_chi_a,imax_chi_a)'
c         write(*,*)'chi_nwka, nmax_chi_a,imax_chi_a',
c     &              chi_nwka, nmax_chi_a,imax_chi_a
c         write(*,*)'Please change parameter in csplcoef_chi'
c         write(*,*)'for parameter (chi_nwka=1+3*imax_chi_a)'
c         write(*,*)'and recompile code'  
c         write(*,*)'**************************************' 
c         stop
c      endif
 
 
c      chi_3d(1:imax_chi,1:nmax_chi,np)=chi(0:imax_chi-1,0:nmax_chi-1)
c      psis(n_radial)=psis
c      dene_chi(1:npsi0a)=dene
c      teme_chi(1:npsi0a)=teme
c      bmin_chi(n_radial)=bmin
c      bmax_chi(n_radial)=bmax
c      th0max_chi(n_radial)=th0max
c      th0a_chi(1:imax_chi_a,n_radial)=th0a(0:imax_chi_a-1)

c---------------------------------------------------
c     allocate pointer wk_chi
c--------------------------------------------------
      chi_nwka=1+3*max0(nmax_chi,imax_chi)
      allocate( wk_chi(1:chi_nwka),STAT=istat)
      call bcast(wk_chi,0.d0,SIZE(wk_chi))
 
c-----creates 2D spline coefficient for chi function
      ibd(1)=2
      ibd(2)=2
      ibd(3)=4
      ibd(4)=4

cSAP090827
c      idm=imax_chi_a
      idm=imax_chi

c      write(*,*)'splcoef_chi umax,nmax_chi',umax,nmax_chi
c      write(*,*)' splcoef_chi ua',ua
      do n=1,nmax_chi
        ua(n)=(n-0.5d0)*umax/nmax_chi
c        write(*,*)'n,ua(n)',n,ua(n)
      enddo
c      write(*,*)' splcoef_chi 1 ua',ua

      do n_radial=1,npsi0
        do n=1,nmax_chi
          chi_tt(1,n,n_radial)=0.d0
          chi_tt(imax_chi,n,n_radial)=0.d0
        enddo

c        write(*,*)' splcoef_ch n_radial',n_radial
c        write(*,*)' th0a_chi(i,n_radial)'
c        write(*,*) (th0a_chi(i,n_radial),i=1,imax_chi)
  
        call coeff2_Sm(imax_chi,th0a_chi(1,n_radial),nmax_chi,ua,
     &   chi_3d(1,1,n_radial),
     &   chi_tt(1,1,n_radial),chi_uu(1,1,n_radial),
     &   chi_ttuu(1,1,n_radial),
cSAP090827
c     &   idm,ibd,wk_chi,imax_chi_a)
     &   idm,ibd,wk_chi,imax_chi)
      enddo
c--------------------------------------------
      deallocate(wk_chi,STAT=istat)
c------------------------------------------

c      write(*,*)' th0a_chi(1,n_radial)'
c      write(*,*) (th0a_chi(i,n_radial),i=1,imax_chi)

c      do n_radial=1,npsi0
c          write(*,*)'in splcoef_chi 1 n_radial',n_radial
c          do i=1,imax_chi
c            write(*,*)'i,th0a_chi(i,n_radial)',
c     &                 i,th0a_chi(i,n_radial)
c          enddo
c      enddo

      goto100

ctest spline
      idm=imax_chi
      norm=0.d0
      norm_dchi_du=0.d0
      norm_dchi_dtheta=0.d0
  
      hu=0.5d0*(ua(2)-ua(1))
      do n_radial=1,npsi0
        htheta=0.5d0*(th0a_chi(2,n_radial)-th0a_chi(1,n_radial))
        do n=1,nmax_chi
           u0=ua(n)  
           do i=1,imax_chi
              theta0=th0a_chi(i,n_radial)
              chi_spline=terp2p_Sm(theta0,u0,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP090827
c     &        chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi)
              write(*,*)'n_radial,n,i',n_radial,n,i
              write(*,*)'chi_3d(i,n,n_radial),chi_spline',
     &                   chi_3d(i,n,n_radial),chi_spline

              norm=norm+dabs(chi_3d(i,n,n_radial)-chi_spline)

  

cc-------------dchi_d_theta
cSAP081227
c              chi_spline=terp2p_Sm(theta0,u0,
              dchi_dtheta_spline=terp2p_Sm(theta0,u0,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP090827
c     &        chi_ttuu(1,1,n_radial),idm,1,0,0,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,1,0,0,imax_chi)

              write(*,*)'spline dchi_dtheta=',dchi_dtheta_spline

              if ((i.ne.1).and.(i.ne.imax_chi)) then
              chi_spline=(chi_3d(i+1,n,n_radial)-chi_3d(i-1,n,n_radial))
     &                 /(th0a_chi(i+1,n_radial)-th0a_chi(i-1,n_radial))
              else
                if (i.eq.1) then
                   chi_spline=(-1.5d0*chi_3d(1,n,n_radial)
     &                       + 2.d0*chi_3d(2,n,n_radial)
     &                       -0.5d0*chi_3d(3,n,n_radial))/
     &                 (th0a_chi(2,n_radial)-th0a_chi(1,n_radial))
                else
                   chi_spline=(1.5d0*chi_3d(imax_chi,n,n_radial)
     &                       - 2.d0*chi_3d(imax_chi-1,n,n_radial)
     &                       +0.5d0*chi_3d(imax_chi-2,n,n_radial))/
     &                 (th0a_chi(2,n_radial)-th0a_chi(1,n_radial))
                endif
              endif

              norm_dchi_dtheta=norm_dchi_dtheta+
     &                         dabs(dchi_dtheta_spline-chi_spline)
              write(*,*)'difference dchi_dtheta=',  chi_spline
c-------------dchi_d_u
cSAP081227
c              chi_spline=terp2p_Sm(theta0,u0,
              dchi_du_spline=terp2p_Sm(theta0,u0,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP080927
c     &        chi_ttuu(1,1,n_radial),idm,0,1,0,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,0,1,0,imax_chi)
              write(*,*)'spline dchi_du=', dchi_du_spline

              if ((n.ne.1).and.(n.ne.nmax_chi)) then
              chi_spline=(chi_3d(i,n+1,n_radial)-chi_3d(i,n-1,n_radial))
     &                 /(ua(n+1)-ua(n-1))
              else
                if (n.eq.1) then
                   chi_spline=(-1.5d0*chi_3d(i,1,n_radial)
     &                       + 2.d0*chi_3d(i,2,n_radial)
     &                       -0.5d0*chi_3d(i,3,n_radial))/
     &                        (ua(2)-ua(1))
                else
                   chi_spline=(1.5d0*chi_3d(i,nmax_chi,n_radial)
     &                       - 2.d0*chi_3d(i,nmax_chi-1,n_radial)
     &                       +0.5d0*chi_3d(i,nmax_chi-2,n_radial))/
     &                       (ua(2)-ua(1))
                endif
              endif
              write(*,*)'difference dchi_du=',  chi_spline
ctest   
              norm_dchi_du=norm_dchi_du+
     &                   dabs(dchi_du_spline-chi_spline)

              call  interpolation_chi(2,n_radial,u0,theta0,
     &        chi,d_chi_d_u_out,d_chi_d_theta_out)

             write(*,*)'chi,d_chi_d_u_out,d_chi_d_theta_out',
     &                  chi,d_chi_d_u_out,d_chi_d_theta_out

c-------------------in intermediate points                        
               chi_spline=terp2p_Sm(theta0+htheta,u0+hu,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP090827
c     &        chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,0,0,1,imax_chi)

              write(*,*)'1chi_spline',chi_spline
c-------------dchi_d_theta
              chi_spline=terp2p_Sm(theta0+htheta,u0+hu,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP090827
c     &        chi_ttuu(1,1,n_radial),idm,1,0,0,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,1,0,0,imax_chi)
              write(*,*)'1spline dchi_dtheta=',  chi_spline
c-------------dchi_d_u
              chi_spline=terp2p_Sm(theta0+htheta,u0+hu,
     &        imax_chi,th0a_chi(1,n_radial),
     &        nmax_chi,ua,
     &        chi_3d(1,1,n_radial),chi_tt(1,1,n_radial),
     &        chi_uu(1,1,n_radial),
cSAP090827
c     &        chi_ttuu(1,1,n_radial),idm,0,1,0,imax_chi_a)
     &        chi_ttuu(1,1,n_radial),idm,0,1,0,imax_chi)

              write(*,*)'1spline dchi_du=',  chi_spline

           
              call  interpolation_chi(2,n_radial,u0+hu,theta0+htheta,
     &        chi,d_chi_d_u_out,d_chi_d_theta_out)
             write(*,*)'1 chi,d_chi_d_u_out,d_chi_d_theta_out',
     &                  chi,d_chi_d_u_out,d_chi_d_theta_out
           enddo
        enddo
      enddo 

      write(*,*)'norm=',norm
      write(*,*)'norm_dchi_du=',norm_dchi_du
      write(*,*)'norm_dchi_dtheta=',norm_dchi_dtheta

 100  continue
c      stop 'lh_ql_flux.f splcoef_chi 100'
      return
      end




      subroutine interpolation_chi(ientry,n_radial0,u0,theta0,
     &chi,d_chi_d_u_out,d_chi_d_theta_out) 

      implicit none
      include 'param.i'
      include 'adj.i'

c-----input
      integer
     & ientry,  !=1 to create arrays of derivatives
                !=2 to calculate chi function and derivatives 
     & n_radial0 !number of radial point at radial chi mesh
 
      real*8
     &u0,       !momentum p/(m*unorm)
     &theta0    !pitch angle
c-----output
      real*8   
     &chi,d_chi_d_u_out,d_chi_d_theta_out 
c-----locals
      integer i,n,n_radial,
     & i_minus,i_plus,n_minus,n_plus

      real*8
     &theta0_loc,u0_loc,h_u,h_theta

c      write(*,*)'interpolation_chi ientry,n_radial0,u0,theta0',
c     & ientry,n_radial0,u0,theta0
c      write(*,*)'ua',ua    
c     do n_radial=1,npsi0 
c         do i=1,imax_chi
c             write(*,*)'n_radial,i,th0a_chi(i,n_radial)',
c     &                  n_radial,i,th0a_chi(i,n_radial)
c          enddo
c      enddo
  
      if (ientry.eq.1) then
c-----------------------------------------------------------------
c       calculate arrays of derivatives: d_chi_d_u,d_chi_d_theta
c-----------------------------------------------------------------
        do n_radial=1,npsi0
          do n=1,nmax_chi
            do i=1,imax_chi
c              write(*,*)'n_radial,n,i',n_radial,n,i
c-------------d_chi_d_theta--------------------------------------------
              if((i.ne.1).and.(i.ne.imax_chi)) then
               d_chi_d_theta(i,n,n_radial)=
     &          (chi_3d(i+1,n,n_radial)-chi_3d(i-1,n,n_radial))
     &          /(th0a_chi(i+1,n_radial)-th0a_chi(i-1,n_radial))
              else
                if(i.eq.1) then
                  d_chi_d_theta(1,n,n_radial)=
     &            (-1.5d0*chi_3d(1,n,n_radial)
     &            + 2.d0*chi_3d(2,n,n_radial)
     &            -0.5d0*chi_3d(3,n,n_radial))/
     &            (th0a_chi(2,n_radial)-th0a_chi(1,n_radial))
                else
                  d_chi_d_theta(imax_chi,n,n_radial)=
     &            (1.5d0*chi_3d(imax_chi,n,n_radial)
     &            -2.d0*chi_3d(imax_chi-1,n,n_radial)
     &            +0.5d0*chi_3d(imax_chi-2,n,n_radial))/
     &            (th0a_chi(2,n_radial)-th0a_chi(1,n_radial))
                endif
              endif

c-------------d_chi_d_u-------------------------------------------

              if((n.ne.1).and.(n.ne.nmax_chi)) then
                d_chi_d_u(i,n,n_radial)=
     &          (chi_3d(i,n+1,n_radial)-chi_3d(i,n-1,n_radial))
     &          /(ua(n+1)-ua(n-1))
c               write(*,*)'n_radial,n,i,ua(n+1),ua(n-1),ua(n+1)-ua(n-1)'
c    &                  ,n_radial,n,i,ua(n+1),ua(n-1),ua(n+1)-ua(n-1)
c             write(*,*)'chi_3d(i,n+1,n_radial),chi_3d(i,n-1,n_radial)'
c     &                  ,chi_3d(i,n+1,n_radial),chi_3d(i,n-1,n_radial)
c              write(*,*)'d_chi_d_u(i,n,n_radial)',
c     &                   d_chi_d_u(i,n,n_radial)
              else
                if(n.eq.1) then
                  d_chi_d_u(i,1,n_radial)=
     &            (-1.5d0*chi_3d(i,1,n_radial)
     &            + 2.d0*chi_3d(i,2,n_radial)
     &            -0.5d0*chi_3d(i,3,n_radial))/
     &            (ua(2)-ua(1))
c              write(*,*)'n_radial,n,i,ua(2),ua(1),ua(2)-ua(1)'
c     &                  ,n_radial,n,i,ua(2),ua(1),ua(2)-ua(1)
c              write(*,*)'chi_3d(i,1,n_radial),chi_3d(i,2,n_radial),
c     &chi_3d(i,3,n_radial)'
c     &                  ,chi_3d(i,1,n_radial),chi_3d(i,2,n_radial)
c     &,chi_3d(i,3,n_radial)
c              write(*,*)'d_chi_d_u(i,1,n_radial)',
c     &                   d_chi_d_u(i,1,n_radial)
                else
                  d_chi_d_u(i,nmax_chi,n_radial)=
     &            (1.5d0*chi_3d(i,nmax_chi,n_radial)
     &            - 2.d0*chi_3d(i,nmax_chi-1,n_radial)
     &            +0.5d0*chi_3d(i,nmax_chi-2,n_radial))/
     &            (ua(2)-ua(1))
                endif
              endif

            enddo ! i
          enddo   ! n
        enddo     ! n_radial

c        do n_radial=1,npsi0
c          do n=1,nmax_chi
c            do i=1,imax_chi
c               write(*,*)'n_radial,n,i',n_radial,n,i
c               write(*,*)'chi,d_chi_d_u,d_chi_d_theta',
c     &         chi_3d(i,n,n_radial),
c     &         d_chi_d_u(i,n,n_radial),d_chi_d_theta(i,n,n_radial)
c            enddo ! i
c          enddo   ! n
c        enddo     ! n_radial

c      stop 'interpolation_chi ientry=1'
  
      endif     !ientry=1

      if (ientry.eq.2) then
c-----------------------------------------------------------------------
c        calculate chi function and its derivatives in the given point
c        (u0,theta0) at n_radial0 with small polidal flux psis(n_radial0) 
c-----------------------------------------------------------------------
c        find number of point at u0 mesh
c-----------------------------------------------------------------------
         u0_loc=u0

         if(u0.le.ua(1)) then 
            n_minus=1
            n_plus =2
            u0_loc=ua(1)
         else
           if(u0.ge.ua(nmax_chi)) then
             n_minus=nmax_chi-1     
             n_plus =nmax_chi
             u0_loc=ua(nmax_chi)
           else
             do n=2,nmax_chi
                if (u0.le.ua(n)) then
                  n_minus=n-1
                  n_plus =n
                  go to 10
               endif
             enddo
 10          continue
           endif
         endif
c----------------------------------------------------------------------
c        find number of point at theta mesh
c-----------------------------------------------------------------------
         theta0_loc=theta0
c         write(*,*)'theta0,th0a_chi(1,n_radial0)',
c     &              theta0,th0a_chi(1,n_radial0)
         if(theta0.le.th0a_chi(1,n_radial0)) then 
            i_minus=1
            i_plus =2
            theta0_loc=th0a_chi(1,n_radial0)
c            write(*,*)'i=1 i_minus,i_plus',i_minus,i_plus
         else
c            write(*,*)'theta0,th0a_chi(imax_chi,n_radial0)',
c     &      theta0,th0a_chi(imax_chi,n_radial0)

           if(theta0.ge.th0a_chi(imax_chi,n_radial0)) then
             i_minus=imax_chi-1     
             i_plus =imax_chi
             theta0_loc=th0a_chi(imax_chi,n_radial0)
           else
cSAP081224
c             do i=2,nmax_chi
             do i=2,imax_chi
c                write(*,*)'i,theta0,th0a_chi(i,n_radial0)',
c     &                     i,theta0,th0a_chi(i,n_radial0)
                if (theta0.le.th0a_chi(i,n_radial0)) then
                  i_minus=i-1
                  i_plus =i
                  goto 20
               endif
             enddo
 20          continue
           endif
         endif
c--------------------------------------------------------------
         h_u=ua(2)-ua(1)
         h_theta=th0a_chi(2,n_radial0)-th0a_chi(1,n_radial0)
c         write(*,*)'th0a_chi(2,n_radial0),th0a_chi(1,n_radial0)',
c     &   th0a_chi(2,n_radial0),th0a_chi(1,n_radial0)

c         write(*,*)'n_radial0,i_minus,i_plus,n_minus,n_plus',
c     &   n_radial0,i_minus,i_plus,n_minus,n_plus
c         write(*,*)'chi_3d(i_plus,n_plus,n_radial0)', 
c     &   chi_3d(i_plus,n_plus,n_radial0)   
c         write(*,*)'d_chi_d_u(i_plus,n_plus,n_radial0)', 
c     &   d_chi_d_u(i_plus,n_plus,n_radial0)   
  
c        write(*,*)'d_chi_d_theta(i_plus,n_plus,n_radial0)', 
c     &   d_chi_d_theta(i_plus,n_plus,n_radial0)   
  

c         write(*,*)'u0_loc,theta_loc',u0_loc,theta0_loc
c         write(*,*)'h_u,h_theta',h_u,h_theta
c         write(*,*)'th0a_chi(i_plus,n_radial0),theta0_loc',
c     &   th0a_chi(i_plus,n_radial0),theta0_loc

c         write(*,*)'th0a_chi(i_minus,n_radial0)',
c     &    th0a_chi(i_minus,n_radial0)

c         write(*,*)'ua(n_plus),u0_loc,ua(n_minus)',
c     &              ua(n_plus),u0_loc,ua(n_minus)

c         write(*,*)'th0a_chi(i_plus,n_radial0)-theta0_loc',
c     &              th0a_chi(i_plus,n_radial0)-theta0_loc
c         write(*,*)'theta0_loc-th0a_chi(i_minus,n_radial0)',
c     &   theta0_loc-th0a_chi(i_minus,n_radial0)
c         write(*,*)' (ua(n_plus)-u0_loc)', (ua(n_plus)-u0_loc)
c         write(*,*)' (u0_loc-ua(n_minus))', (u0_loc-ua(n_minus))
cSAP101223
c         write(*,*)'u0,theta0',u0,theta0
c         write(*,*)'th0a_chi(imax_chi,n_radial0)',
c     &              th0a_chi(imax_chi,n_radial0)

c         write(*,*)'i_minus,n_minus,n_radial0',i_minus,n_minus,n_radial0
c         write(*,*)'i_plus',i_plus
         chi=(chi_3d(i_minus,n_minus,n_radial0)*
     &        (th0a_chi(i_plus,n_radial0)-theta0_loc)
     &   +chi_3d(i_plus,n_minus,n_radial0)*
     &        (theta0_loc-th0a_chi(i_minus,n_radial0)))*
     &   (ua(n_plus)-u0_loc)
     &   +(chi_3d(i_minus,n_plus,n_radial0)*
     &        (th0a_chi(i_plus,n_radial0)-theta0_loc)
     &   +chi_3d(i_plus,n_plus,n_radial0)*
     &       (theta0_loc-th0a_chi(i_minus,n_radial0)))*
     &   (u0_loc-ua(n_minus))

         chi=chi/(h_u*h_theta)
c         write(*,*)'chi',chi 

c         write(*,*)'d_chi_d_u(i_minus,n_minus,n_radial0)',
c     &              d_chi_d_u(i_minus,n_minus,n_radial0)

c         write(*,*)'d_chi_d_u(i_plus,n_minus,n_radial0)',
c     &              d_chi_d_u(i_plus,n_minus,n_radial0)
 
c        write(*,*)'d_chi_d_u(i_minus,n_plus,n_radial0)',
c     &             d_chi_d_u(i_minus,n_plus,n_radial0)

c        write(*,*)'d_chi_d_u(i_plus,n_plus,n_radial0)',
c     &             d_chi_d_u(i_plus,n_plus,n_radial0)

c        write(*,*)'(th0a_chi(i_plus,n_radial0)-theta0_loc)',
c     & (th0a_chi(i_plus,n_radial0)-theta0_loc)
c        write(*,*)'(theta0_loc-th0a_chi(i_minus,n_radial0))',
c     &(theta0_loc-th0a_chi(i_minus,n_radial0))

c      write(*,*)'(ua(n_plus)-u0_loc)',(ua(n_plus)-u0_loc)
c      write(*,*)'(u0_loc-ua(n_minus))',(u0_loc-ua(n_minus))

         d_chi_d_u_out=(d_chi_d_u(i_minus,n_minus,n_radial0)*
     &        (th0a_chi(i_plus,n_radial0)-theta0_loc)
     &   +d_chi_d_u(i_plus,n_minus,n_radial0)*
     &        (theta0_loc-th0a_chi(i_minus,n_radial0)))*
     &   (ua(n_plus)-u0_loc)
     &   +(d_chi_d_u(i_minus,n_plus,n_radial0)*
     &        (th0a_chi(i_plus,n_radial0)-theta0_loc)
     &   +d_chi_d_u(i_plus,n_plus,n_radial0)*
     &       (theta0_loc-th0a_chi(i_minus,n_radial0)))*
     &   (u0_loc-ua(n_minus))

         d_chi_d_u_out=d_chi_d_u_out/(h_u*h_theta)
c         write(*,*)'d_chi_d_u_out',d_chi_d_u_out

         d_chi_d_theta_out=(d_chi_d_theta(i_minus,n_minus,n_radial0)*
     &        (th0a_chi(i_plus,n_radial0)-theta0_loc)
     &   +d_chi_d_theta(i_plus,n_minus,n_radial0)*
     &        (theta0_loc-th0a_chi(i_minus,n_radial0)))*
     &   (ua(n_plus)-u0_loc)
     &   +(d_chi_d_theta(i_minus,n_plus,n_radial0)*
     &        (th0a_chi(i_plus,n_radial0)-theta0_loc)
     &   +d_chi_d_theta(i_plus,n_plus,n_radial0)*
     &       (theta0_loc-th0a_chi(i_minus,n_radial0)))*
     &   (u0_loc-ua(n_minus))
         d_chi_d_theta_out= d_chi_d_theta_out/(h_u*h_theta)
      endif     !ientry=2

      return
      end

c-----Maxwellian function f_m(u)~exp(-mc**2*gamma/T)
      subroutine maxwel_adj(nmax_chi_a,ua,c,nmax,fm,t)
c-------------------------------------------------------------------
c     input 
c     ua is the normalized momentum mesh p/(m*unorm) 
c     du is the mesh step
c     nmax is the number of mesh points
c     nmax_chi_a maximal value of nmax
c     t    is the temperature in units of m*u_n**2
c          If you want velocities normalized to sqrt(T/m), set t=1.
c     c    cligth/unorm
c--------------------------------------------------------------------
c     output
c     fm   relativistic Maxwellian function array
c          normalized to 4*pi*integral[fm*du*u**2]=1
c--------------------------------------------------------
      implicit none

      integer nmax_chi_a,nmax,n
      real*8 du,c,t,sum,pi,ua(0:nmax_chi_a-1),fm(0:nmax_chi_a-1)
      sum=0
      du=ua(2)-ua(1)
c      write(*,*)'max_adj nmax'
      do n=0,nmax-1
        fm(n)=exp(-ua(n)**2/(sqrt(1+(ua(n)/c)**2)+1)/t)
c        write(*,*)'n,fm(n),ua(n)',n,fm(n),ua(n)
        sum=sum+ua(n)**2*fm(n)
      end do
      pi=4*atan(1.0d0)
      sum=4*pi*du*sum
      do n=0,nmax-1
      fm(n)=fm(n)/sum
      end do
      return
      end

      subroutine d_chi(du,nmax,dth0,imax,imax1,
     &chi,rho_)
c-----calculate derivatives from adjoint function chi
c     output: dchi_du and dchi_dth0 wil be in common flux_w.i
c
c     df/dv_parallel=(1/v**2)*d[v**2*cos(theta)*f]/dv 
c                    -1/(v**2*sin(theta))*d[v*sin**2(theta)*f)]/d(theta)=
c                   =cos(theta)*df/dv-sin(theta)/v*df/d(theta)     
      implicit none

c      include 'param_add.i' 
c      include 'flux_w.i'

      include 'param.i'

      integer nmax,imax,imax1,n,i
      real*8 du,dth0,chi,rho_, dchi_du, dchi_dth0
      dimension  chi(0:imax1-1,0:nmax-1)
      dimension  dchi_du(0:imax1-1,0:nmax-1)
      dimension  dchi_dth0(0:imax1-1,0:nmax-1)

      write(*,*)'d_chi nmax,imax,imax1',nmax,imax,imax1
      
      do i=0,imax-1
        n=0
        dchi_du(i,n)=(chi(i,n+1)+chi(i,n))/(2.d0*du)
        do n=1,nmax-2    
           dchi_du(i,n)=(chi(i,n+1)-chi(i,n-1))/(2.d0*du)
        enddo
        n=nmax-1
        dchi_du(i,n)=(1+rho_)*(chi(i,n)-chi(i,n-1))/(2.d0*du)                  
      end do

      do n=0,nmax-1
        i=0
cSm020923
        dchi_dth0(i,n)=(chi(i+1,n)-chi(i,n))/(2.d0*dth0)
        do i=1,imax-2    
           dchi_dth0(i,n)=(chi(i+1,n)-chi(i-1,n))/(2.d0*dth0)
        enddo
        i=imax-1
        dchi_dth0(i,n)=(-chi(i-1,n))/(2.d0*dth0)                  
      end do

      return
      end           

      subroutine CD_adj_LH_efficiency_anal_1(v_ph,n_parallel,
     &n_radial0,theta_pol,
     &unorm,
     &lh_cd_efficiency)
c------------------------------------------------------
c     calculate LH CD efficiency for anallytical 
c     flux using adj chi function
c     at one radial point psis(n_radial0) from chi radial mesh 
c------------------------------------------------------
      implicit none
      include 'param.i'
      include 'adj.i'
c-----input
      integer
     & n_radial0           ! is a number of radial ch mesh

      real*8      
     &v_ph,                !parallel V_ph/unorm
     &n_parallel,   
     &theta_pol,           !poloidal angle [radian]
     &length               !field line length normalized to  
c     &t,                   !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1. 
    
c-----output
      real*8
     &pow_dens,            !Power density normalized to m*u_norm**2*n/tau_n 
     &lh_cd_dens,          !LH current drive density normalized to q*u_norm*n
                           !Here n is the density  
     &lh_cd_efficiency,
     &unorm                !sqrt(T/m) cm/sec
c-----external
      real*8 temperho,densrho,
     &bmax_psi,bmin_psi,psi_rho,rhopsi


c-----local
      integer i,n,n_radial,idm,nmin
      real*8 du,u0,     
     &rho_small,           !small radius
     &temp_kev,            !electron temperatura [KeV]
     &dense,               !electron dencity 10**13/cm**3     
     &b_ratio,             !B/B_min
     &sin_theta0_res,cos_theta0_res,dfm_du0,theta0_res,theta0_odd,
     &Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0,
     &gamma,pi,bmin,bmax,deltb,psi,th0max,sin_trap,
     &d_chi_du0,d_chi_dtheta0,u0max,chi,
     &weight_radial,
     &umin,               !minimal u0 at the resonance LH hyperbola
     &clight,              !light speed [cm/sec]
     &arg1,cln,
cSAP090827
cSAP090503
c     &fm(0:nmax_chi_a-1),c,dchi_dupar,
     &c,dchi_dupar,
cfor test
     &del_power,del_cd,cos_theta_res,sin_theta_res,Gamma_lh_anal,
     &Gamma_lh_anal_u0,Gamma_lh_anal_theta0
         
cSAP090827
      real*8, pointer :: fm(:)   !0:nmax_chi-1
      integer istat
c-----external
      real*8 terp2p_Sm,b_d_bmin
 
c---------------------------------------------------
c     allocate pointer fm 
c--------------------------------------------------
      allocate( fm(0:nmax_chi-1),STAT=istat)
      call bcast(fm,0.d0,SIZE(fm))
c----------------------------------------------
      
      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]
 
    
      
      rho_small=rhopsi(psis(n_radial0))
      psi=psis(n_radial0)
      write(*,*)'##CD_adj_LH_efficiency_anal_1'
      write(*,*)'n_radial0=',n_radial0

      write(*,*)'psis(n_radial0),rho_small',psis(n_radial0),rho_small

      temp_kev=temperho(rho_small,1)
      dense=densrho(rho_small,1) 
      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m) 
      c=clight/unorm
      
c      write(*,*)'ua',ua 
cSAP090827     
c      call  maxwel_adj(nmax_chi_a,ua,c,nmax_chi,fm,t)
      call  maxwel_adj(nmax_chi,ua,c,nmax_chi,fm,t)

      bmax= bmax_psi(psi)
      bmin= bmin_psi(psi)	 
      
      deltb = (bmax-bmin)/bmin
      th0max = atan2(1.0d0,sqrt(deltb))
      
      sin_trap=dsqrt(bmin/bmax)
      write(*,*)'sin_trap=',sin_trap,'dsin(th0max)=',dsin(th0max)

      du = umax/nmax_chi

      pow_dens=0.d0
      lh_cd_dens=0.d0

      umin=v_ph
      write(*,*)'anal temp_kev,t',temp_kev,t
      write(*,*)'clight,temp_kev/t,unorm,n_parallel',
     &clight,temp_kev/t,unorm,n_parallel
      write(*,*)'anal umin=c/n_parallel=',c/n_parallel
cSAP070908
      umin=c/dsqrt(n_parallel**2-1.d0)

      v_ph=umin    

      write(*,*)'anal umin',umin,'umin/dsqrt(2.d0)',umin/dsqrt(2.d0)

      nmin=(dabs(umin)/du+0.5d0)
      nmin=(dabs(umin)/du-0.5d0)

      nmin=nmin+1
      u0=(nmin + 0.5E0)*du

c      write(*,*)'umin,du,(dabs(umin)/du-0.5d0),nmin,u0',
c     &umin,du,(dabs(umin)/du-0.5d0),nmin,u0

      b_ratio=b_d_bmin(n_radial0,theta_pol)

c     write(*,*)'in CD_adj_LH_efficiency_anal_1  nmin,nmax_chi-1',
c     &nmin,nmax_chi-1

      do n=nmin,nmax_chi-1
         u0=(n + 0.5E0)*du

c         write(*,*)'n,u0',n,u0

         gamma=dsqrt(1.d0+(u0/c)**2) 

c         write(*,*)'gamma,c/n_parallel,u0',
c     &             gamma,c/n_parallel,u0

         cos_theta_res=gamma*c/(u0*n_parallel)       
         sin_theta_res=dsqrt(1.d0-cos_theta_res**2)
         sin_theta0_res= sin_theta_res/b_ratio
         cos_theta0_res=dsqrt(1.d0-sin_theta0_res**2)
     
c         write(*,*)'b_ratio',b_ratio
c         write(*,*)'cos_theta_res,sin_theta_res', 
c     &   cos_theta_res,sin_theta_res
c         write(*,*)'cos_theta0_res,sin_theta0_res', 
c     &   cos_theta0_res,sin_theta0_res

c        write(*,*)'fm(n)',fm(n)

        if (cos_theta_res.lt.0d0) cos_theta0_res=-cos_theta0_res
 
        Gamma_lh_anal=-fm(n)*(u0/gamma)*cos_theta_res

c        write(*,*)' Gamma_lh_anal',Gamma_lh_anal

        Gamma_lh_anal_u0=Gamma_lh_anal*cos_theta_res

c        write(*,*)' Gamma_lh_anal_u0', Gamma_lh_anal_u0

c        Gamma_lh_anal_theta0=-Gamma_lh_anal*sin_theta_res*
c     &  (sin_theta0_res/cos_theta0_res)/(sin_theta_res/cos_theta_res)
        Gamma_lh_anal_theta0=-Gamma_lh_anal_u0*
     &                        (sin_theta0_res/cos_theta0_res)
   

        dg_d_theta0=(u0/gamma)*
     &  b_ratio*sin_theta0_res*cos_theta0_res/cos_theta_res

c        write(*,*)'dg_d_theta0',dg_d_theta0
              
c        call lh_ql_flux(n_parallel,n_perp,u0,umin,t,temp_kev,
c     &                   b_ratio,e_y,e_z,
c     &                   sin_theta0_res,cos_theta0_res, dfm_du0,
c     &                   Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0,
c     &                   gamma,unorm)

c         write(*,*)'n,Gamma_a_LH_u0,sin_theta0_res,cos_theta0_res',
c     &              n,Gamma_a_LH_u0,sin_theta0_res,cos_theta0_res
c         write(*,*)'dfm_du0,gamma,dg_d_theta0,dfm_du0,gamma,dg_d_theta0
c     &            ',dfm_du0,gamma,dg_d_theta0,dfm_du0,gamma,dg_d_theta0

c         pow_dens=pow_dens+du*(u0**4/gamma**2)*
c     &            sin_theta0_res*cos_theta0_res*Gamma_a_LH_u0*dfm_du0/
c     &            dg_d_theta0


c         del_power=du*(u0**4/gamma**2)*
c     &            sin_theta0_res*cos_theta0_res*Gamma_a_LH_u0*dfm_du0/
c     &            dg_d_theta0

c       pow_dens=pow_dens+du*u0**2*(u0/gamma*cos_theta_res)*
c     &                 Gamma_lh_anal/dg_d_theta0

       pow_dens=pow_dens+du*u0**2*sin_theta0_res*
     &                   (u0/gamma*cos_theta_res)*
     &                   Gamma_lh_anal/dg_d_theta0


       del_power=du*u0**2*sin_theta0_res*
     &                   (u0/gamma*cos_theta_res)*
     &                   Gamma_lh_anal/dg_d_theta0

c         write(*,*)'del_power=',del_power

         theta0_res=dacos(cos_theta0_res)                                                                                                                                  
cSAP090827                                    
c         idm=imax_chi_a
         idm=imax_chi

c         write(*,*)'sin_theta0_res,sin_trap',sin_theta0_res,sin_trap
c         write(*,*)'dsin(th0a_chi(imax_chi,n_radial0))',
c     &              dsin(th0a_chi(imax_chi,n_radial0))
c         write(*,*)'theta0_res,th0a_chi(imax_chi,n_radial0)',
c     &              theta0_res,th0a_chi(imax_chi,n_radial0)
c         write(*,*)'pi-th0a_chi(imax_chi,n_radial0)',
c     &              pi-th0a_chi(imax_chi,n_radial0)
     
         if((theta0_res.le.th0a_chi(imax_chi,n_radial0)).or.
     &     (theta0_res.ge.pi-th0a_chi(imax_chi,n_radial0))) then          
c         if (sin_theta0_res.lt.sin_trap) then
c-----------passing particles
c            write(*,*)'theta0_res,u0',theta0_res,u0
            if (theta0_res.le.th0a_chi(imax_chi,n_radial0)) then

cSAP090129
              if (i_chi_interpolation.eq.1) then
c-----------------spline interpolation
                 d_chi_du0=terp2p_Sm(theta0_res,u0,
     &           imax_chi,th0a_chi(1,n_radial0),
     &           nmax_chi,ua,
     &           chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &           chi_uu(1,1,n_radial0),
cSAP090827
c     &           chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi_a)
     &           chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)

                 d_chi_dtheta0=terp2p_Sm(theta0_res,u0,
     &           imax_chi,th0a_chi(1,n_radial0),
     &           nmax_chi,ua,
     &           chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &           chi_uu(1,1,n_radial0),
cSAP090827
c     &           chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi_a)      
     &           chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)      
              endif
cSAP090129
              if (i_chi_interpolation.eq.2) then
                 call  interpolation_chi(2,n_radial0,u0,theta0_res,
     &           chi,d_chi_du0,d_chi_dtheta0)
               endif
c              write(*,*)'interpolation d_chi_du0,d_chi_dtheta0',
c     &                    d_chi_du0,d_chi_dtheta0
c               write(*,*)'d_chi_du0,d_chi_dtheta0',
c     &                    d_chi_du0,d_chi_dtheta0
            else
               theta0_odd=pi-theta0_res 
cSAP090129
            if (i_chi_interpolation.eq.1) then
c--------------spline interpolation

               d_chi_du0=-terp2p_Sm(theta0_odd,u0,
     &         imax_chi,th0a_chi(1,n_radial0),
     &         nmax_chi,ua,
     &         chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &         chi_uu(1,1,n_radial0),
cSAP090827
c     &         chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi_a)
     &         chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)      
 
               d_chi_dtheta0=terp2p_Sm(theta0_odd,u0,
     &         imax_chi,th0a_chi(1,n_radial0),
     &         nmax_chi,ua,
     &         chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &         chi_uu(1,1,n_radial0),
cSAP090827
c     &         chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi_a)
     &         chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)
              endif
cSAP090129
              if (i_chi_interpolation.eq.2) then
                  call  interpolation_chi(2,n_radial0,u0,theta0_odd,
     &            chi,d_chi_du0,d_chi_dtheta0)
c              write(*,*)'interpolation d_chi_du0,d_chi_dtheta0',
c     &                    d_chi_du0,d_chi_dtheta0
                  d_chi_du0=-d_chi_du0
              endif

            endif
c                          write(*,*)'d_chi_du0,d_chi_dtheta0',
c     &                    d_chi_du0,d_chi_dtheta0

               lh_cd_dens=lh_cd_dens+du*u0**2*sin_theta0_res*
     &                         Gamma_lh_anal_u0*(d_chi_du0-
     &         (d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))
     &                  /dg_d_theta0
               

               del_cd=du*u0**2*sin_theta0_res*
     &                         Gamma_lh_anal_u0*(d_chi_du0-
     &         (d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))
     &                  /dg_d_theta0

c               write(*,*)'u0,sin_theta0_res,Gamma_lh_anal_u0',
c     &                    u0,sin_theta0_res,Gamma_lh_anal_u0
c               write(*,*)'d_chi_du0,d_chi_dtheta0',
c     &                    d_chi_du0,d_chi_dtheta0

c            write(*,*)'del_cd=',del_cd
c            write(*,*)'del_cd/del_power',del_cd/del_power   

c            write(*,*)'du,u0,gamma, cos_theta0_res,Gamma_a_LH_u0',
c     &                 du,u0,gamma, cos_theta0_res,Gamma_a_LH_u0
          
c             write(*,*)'lh_cd_dens=',lh_cd_dens
         endif 
      enddo

      pow_dens=2.d0*pi*pow_dens
      lh_cd_dens=2.d0*pi*LH_CD_dens

c      write(*,*)'pow_dens',pow_dens
c      write(*,*)'lh_cd_dens',lh_cd_dens
     
      if(pow_dens.ne.0.d0) then
        lh_cd_efficiency=lh_cd_dens/pow_dens
      else
        lh_cd_efficiency=0.d0
      endif

      write(*,*)'in CD_adj_LH_efficiency_anal_1 lh_cd_efficiency=', 
     &lh_cd_efficiency
         
c-----cln -Coulomb Ln
      arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense)	!temp KeV, dense 10**13/cm**3
      cln=24.d0-dlog(arg1)
      lh_cd_efficiency=lh_cd_efficiency*
     &3.84d0*temp_kev/(cln*dense)	       !  (A/m**2)/(W/m**3)
      lh_cd_efficiency=lh_cd_efficiency*1.d-5  !  A/cm**2/(erg/(sec*cm**3)

c----------------------------------------------------
      deallocate(fm,STAT=istat)
c-----------------------------------------------------
      return
      end


c=============================================================================
c     for u_perp integration
c=============================================================================    
    
      subroutine CD_adj_efficiency_1(n_radial0,
     &theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &unorm, 
     &cd_efficiency)
c------------------------------------------------------
c     calculate CD efficiency using adj chi function
c     at one radial point psis(n_radial0) from chi radial mesh 
c------------------------------------------------------
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input
      integer
     & n_radial0           ! is a number of radial ch mesh

      real*8  
     &theta_pol,           !poloidal angle [radian]    
     &n_parallel,          !parallel refractive index
     &n_perp               !perpendicular refractive index
    
      complex*16 
     &e_x,
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1
c-----output
      real*8
     &unorm,               !sqrt(T/m) cm/sec
     &cd_efficiency

c-----external
      real*8 temperho,densrho,
     &bmax_psi,bmin_psi,psi_rho,rhopsi,y,b_d_bmin

c-----local
      integer n_radial,n_harm_adj,
c     &i_resonance_curve_integration_method,
     &kmax,
     &i_calculate_CD   ! =1 calculate power and CD
                       ! at the given adj radial mesh point
                       ! with number n_radial0 
      real*8 u0,     
     &rho_small,           !small radius
     &temp_kev,            !electron temperatura [KeV]
     &dense,               !electron dencity 10**13/cm**3     
     &b_ratio,             !B/B_min
     &pi,bmin,bmax,deltb,psi,th0max,sin_trap,
     &d_chi_du0,d_chi_dtheta0,u0max,chi,
     &clight,             !light speed [cm/sec]
     &arg1,cln,
     &z,r,phi,            !bellow phi is set arbitratry phi=0
     &y_loc,              !electron y at point (z,r,phi)     
     &power_nharm_adj,CD_nharm_adj, !under integral functions
     &pow_dens,            !Power density normalized to m*u_norm**2*n/tau_n 
     &cd_dens,              !current drive density normalized to q*u_norm*n
                           !Here n is the density  
     &p_par_rl(2),         !for check max and min p_par
     &k_1_kev,             !egrs in 1 KeV      (erg)
     &charge_electron      !electron charge (statcoulomb)

      k_1_kev=1.6022d-9          !egrs in 1 KeV      (erg)
      charge_electron=4.8032d-10 !electron charge (statcoulomb)
      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]

      write(*,*)'CD_adj_efficiency_1 n_radial0',n_radial0

c      do n_radial=1,npsi0
c         write(*,*)'psis',psis(n_radial)
c         rho_small=rhopsi(psis(n_radial))
c         write(*,*)'rho_small',rho_small
c      enddo

      rho_small=rhopsi(psis(n_radial0))
      psi=psis(n_radial0)

      write(*,*)'n_radial0=',n_radial0

      write(*,*)'psis(n_radial0),rho_small',psis(n_radial0),rho_small

      temp_kev=temperho(rho_small,1)
      dense=densrho(rho_small,1) 
      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
      bmax= bmax_psi(psi)
      bmin= bmin_psi(psi)
	     
      deltb = (bmax-bmin)/bmin
c      write(*,*)'bmax,mbin,deltb',bmax,bmin,deltb   
      th0max = atan2(1.0d0,sqrt(deltb))
c      write(*,*)'th0max',th0max  mass_e=9.1094d-28    !electron rest mass (g)
c      k_1_kev=1.6022d-9                                !egrs in 1 KeV      (erg)
      sin_trap=dsqrt(bmin/bmax)

c      write(*,*)'sin_trap=',sin_trap,'dsin(th0max)=',dsin(th0max)

      pow_dens=0.d0
      cd_dens=0.d0

c      umin=clight/dsqrt(n_parallel**2-1.d0)/unorm ! normalized umin      
c      write(*,*)'umin',umin,'umin/dsqrt(2.d0)',umin/dsqrt(2.d0)
      
      b_ratio=b_d_bmin(n_radial0,theta_pol)

      write(*,*)'n_radial0,theta_pol,b_ratio',
     &n_radial0,theta_pol,b_ratio

      phi=0.d0
      call zr_psith(psi,theta_pol,z,r)
      y_loc=y(z,r,phi,1)

c      i_resonance_curve_integration_method_adj=3
    
      do n_harm_adj=n_harm_adj_min,n_harm_adj_max
c-------------------------------------------------------------------
c       the loop over used gyro-harmonics 
c----------------------------------------------------------------- 
       write(*,*)'CD_adj_efficiency_1 n_harm_adj=',n_harm_adj

c     & write(*,*)'i_resonance_curve_integration_method_adj=',
c     & i_resonance_curve_integration_method_adj


c---------------------------------------------------------------------
c       for check min,max p_parallel
 
        call root_res(0.d0,n_parallel,n_harm_adj,y_loc,kmax,p_par_rl)

c        write(*,*)'p_perp=0,n_harm_adj,n_parallel,y_loc,kmax,p_par_rl',
c     &   n_harm_adj,n_parallel,y_loc,kmax,p_par_rl

c------------------------------------------------------------------------
c        write(*,*)'before call intgr_relt_power_cd_adj'
c        write(*,*)'n_radial0,b_ratio',n_radial0,b_ratio
c        write(*,*)'y_loc,sin_trap,temp_kev',y_loc,sin_trap,temp_kev
c        write(*,*)'e_x,e_y,e_z',e_x,e_y,e_z
c        write(*,*)'n_harm_adj,n_parallel,n_perp',
c     &             n_harm_adj,n_parallel,n_perp
c        write(*,*)'unorm,umax',unorm,umax
c        write(*,*)'n_relt_intgr_adj',n_relt_intgr_adj
c        write(*,*)'i_resonance_curve_integration_method_adj',
c     &             i_resonance_curve_integration_method_adj
c        write(*,*)'epsi_adj',epsi_adj

        i_calculate_CD=1 ! calculate power and only C) 
                         ! at given adj radial mes point
                         ! with number n_radial0 
        call intgr_relt_power_cd_adj(i_calculate_CD,
     & n_radial0,b_ratio,
     &  y_loc,sin_trap,temp_kev, 
     &  e_x,e_y,e_z,  
     &  n_harm_adj,n_parallel,n_perp, 
     &  unorm,umax,
     &  n_relt_intgr_adj,i_resonance_curve_integration_method_adj,
     &  epsi_adj,
     &  power_nharm_adj,CD_nharm_adj)

        pow_dens=pow_dens +power_nharm_adj
        cd_dens=cd_dens+CD_nharm_adj

        write(*,*)'power_nharm_adj',power_nharm_adj
        write(*,*)'Cd_nharm_adj',CD_nharm_adj

      enddo !n_harm

      pow_dens=2.d0*pi*pow_dens
      cd_dens=-2.d0*pi*cd_dens
c      write(*,*)'pow_dens',pow_dens
c      write(*,*)'cd_dens',cd_dens
     
      if(pow_dens.ne.0.d0) then
        cd_efficiency=cd_dens/pow_dens
      else
        cd_efficiency=0.d0
      endif

      write(*,*)'in CD_adj_efficiency_1 cd_efficiency=',
     &cd_efficiency
         
c-----cln -Coulomb Ln
      arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense)	!temp KeV, dense 10**13/cm**3
      cln=24.d0-dlog(arg1)

c-----adj efficiency is normalized to eta_0_adj=e*tau_n/(m*u_t) 
c     tau_n=T**(3/2)*sqrt(m)/(4*pi*e**4*density*cln)
c     u_t=cqrt(T/m)
c     eta_0_adj=T/[4*pi*e**3*n*cln]|cgs=T_kev/(n_13*cln)*
c           [k_1_kev/(4*pi*charge_electron**3*10**13)]           !cgs
c      coef=1.6022d+8/(4*pi*4.8032**3)          !cgs  (statampere/cm**2)/(erg/(sec*cm**3))
c      coef=1.6022d-1/(4*pi*4.8032**3*3)=3.8352d-5!     (A/cm**2)/(erg/(sec*cm**3))
c      eta_0_adj=(temp_kev/(cln*dense))*coef

      cd_efficiency=cd_efficiency*
     &3.84d0*temp_kev/(cln*dense)	       !  (A/m**2)/(W/m**3)
c      write(*,*)'temp_kev,dense,cln,cd_efficiency',
c     &temp_kev,dense,cln,cd_efficiency
      cd_efficiency=cd_efficiency*1.d-5  !  (A/cm**2)/(erg/(sec*cm**3))
      write(*,*)'temp_kev,dense,cln,cd_efficiency',
     &temp_kev,dense,cln,cd_efficiency


      return
      end

     

      subroutine intgr_relt_power_cd_adj(
     &i_calculate_CD,
     &n_radial0,b_ratio,y_loc,
     &sin_trap,temp_kev, 
     &e_x,e_y,e_z,  
     &n_harm_adj,n_parallel,n_perp,
     &unorm,umax,
     &n_relt_intgr_adj,i_resonance_curve_integration_method_adj,epsi,
     &power_nharm_adj,CD_nharm_adj)
c-----------------------------------------------------------------
c     calculates power density  power_nharm and 
c     CD density  CD_nharm_adj for one harmonic using ADJ function
c     as  ther integral for the relativistic electrons 
c     integral(0<=p0_perp<=p0_perp_max)f_n_harm(p0_perp)
c     f_n_harm={sum_k=1,2,f_nk(p0_perp)
c     for the EC harmonic with number 'n_harm_adj'
c----------------------------------------------------------------- 
c     input
c        i_calculate_CD   ! =1 calculate power and CD 
c                         ! at given adj radial mesh point
c                         ! with number n_radial0
c                         !
c                         ! =0, calculate power only (no CD) 
c                         ! using input arguments computed at
c                         ! given space point
c                         ! In this case n_radial0 can be arbitary
c 
c       n_radial0    is a number of radial chi mesh
c       b_ratio      =B/B_min
c       y_loc        =|omega_ce|/omega positive value for the electron charge
c       sin_trap     =sqrt(bmin/bmax) at n_radial0 flux surface
c       temp_kev     electron temperature [KeV]
c       e_x,e_y,e_z  Stix frame complex electric filed polarization
c       n_harm_adj        - EC harmonic number
c       n_parallel          N_parallel
c       n_perp      N_perpendicular
c       unorm    = dsqrt(temp_kev/(m_e*t)) normalization velocity [sm/sec]  
c       umax       maximal momentum normalized by unorm
c       n_relt_intgr_adj the number of points for the integration over p0_perp
c----------------------------------------------------------
c     i_resonance_curve_interation_method_adj
c
c     i_resonance_curve_integration_method_adj=1 angle integration
c     now it works for ellipse case only                             c
c
c     i_resonance_curve_integration_method_adj=2 rectangle formula
c     for p0_perp integration
c     
c     i_resonance_curve_integration_method_adj=3 trapezoidal formula
c     for p0_perp integration
c
c     i_resonance_curve_integration_method_adj=4 !adaptive Simpson 
c     for p0_perp integration
c
c       epsi is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------ 
c     
c     output
c       power_nharm_adj,CD_nharm_adj =integral 
c       over 0<=p0_perp=<p0_perp_max      
c----------------------------------------------------------------    
c      The integration method is specified by the variable
c      i_resonance_curve_integration_method_adj.
c      Now this variable is set inside this subroutine.
c----------------------------------------------------------------

      implicit none
c-----input
      real*8 b_ratio,n_parallel,n_perp,y_loc,sin_trap,temp_kev,
     &epsi,
     &unorm,umax 

      integer n_harm_adj,n_relt_intgr_adj,
     &i_resonance_curve_integration_method_adj,n_radial0,
     &i_calculate_CD     

      complex*16
     &e_x,e_y,e_z              !Stix frame electric filed polarization
c-----output
      real*8 power_nharm_adj,CD_nharm_adj

c-----local
      real*8 h,p0_perp,p_perp,p_t,clight,p_int,
     &p_perp_max,  ! max value of the perpendicular momentum divided by mc
                    ! on the resonance ellipse
     &p0_perp_max, ! max value of the perpendicular momentum p0_perp divided by mc                   
     &f_power_n_harm,f_CD_n_harm, !under integral functions
                                   !for ADJ calculations
     &theta_temperature,           !mc**2/T
     &mass_e,             !electron rest mass (g)
     &k_1_kev,            !egrs in 1 KeV      (erg)
     &bk2,                 !Mackdonalds function bk2=K_2(theta)*EXP(theta)	
     &p_loc,pi,
c----to test the normalization
     &p1,dens_test,fm,gamma
c------------------------------------------------------------------
c     for integration along ellipse by the angle
c------------------------------------------------------------------
      real*8 rme,rtem0,xint,cs,sn,thet1,thet2,p_par,
     &p_par_pl,p_perp_min,p_perp_pl,v0dc,vmax1dc,vpar0dc,
     &vpar1dc,vpar2dc,dth,vt0,cdvt,tem0,vmaxdt
c-----------------------------------------------------------------

      integer j,jmax,
     &ires   ! ires     =0 no resonance,=1 ellipse,=2 parabola,=3 hyperbola  
     
     
c     external
c     ec_cond,power_CD_nharm_under_integral_functions
      real*8 temperho
 
c      write(*,*)'intgr_relt_power_cd_adj'
c      write(*,*)'i_resonance_curve_integration_method_dj',
c     & i_resonance_curve_integration_method_adj


      jmax=n_relt_intgr_adj

      mass_e=9.1094d-28         !electron rest mass (g)
      clight=2.99792458d10      !light speed [sm/sec]
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)
      theta_temperature=mass_e*clight**2/(k_1_kev*temp_kev) ! m*c**2/T_e


      call ec_cond(n_harm_adj,y_loc,n_parallel,ires,p_perp_max)

c      write(*,*)'after ec_cond n_parallel,ires,p_perp_max',
c     &          n_parallel,ires,p_perp_max

cSAP090126
      if (ires.eq.0) then !no resonance
         power_nharm_adj=0.d0
         CD_nharm_adj=0.d0
        return
      endif

      if (ires.eq.1) then
c-------ellipse case
cSAP090126
        p_perp_max=dsqrt((n_parallel**2+(n_harm_adj*y_loc)**2-1.d0))/
     &             (1.d0-n_parallel**2) !p/mc

cSAP081115
c        p_perp_max=p_perp_max*clight/unorm     !maximal u_perp normalized to sqrt(T/m)
      else
c-------ires=2,3 hyperbola or parabola cases
c        p_perp_max=umax   ! maximal u_perp normalized to sqrt(T/m)
cSAP081115
        p_perp_max=umax*unorm/clight !p/mc
      endif

c      write(*,*)'ires=',ires,'p_perp_max=', p_perp_max
 
      p0_perp_max=p_perp_max/(b_ratio) !p/mc                

c      write(*,*)'b_ratio,p0_perp_max=',b_ratio,p0_perp_max
cSAP081201
c      if (p0_perp_max.gt.umax) p0_perp_max=umax  
      if (p0_perp_max.gt.(unorm*umax/clight)) 
     &    p0_perp_max=(unorm*umax/clight)    !p/mc
  
c      write(*,*)'umax,p0_perp_max=',umax,p0_perp_max

      h=p0_perp_max/(n_relt_intgr_adj)   ! step of integration over p_perp
                                         ! p/mc
c-----calculations of the integrals over p_perp

      power_nharm_adj=0.d0
      CD_nharm_adj=0.d0

c----------------------------------------------------------
c     i_resonance_curve_interation_method_adj
c
c     i_resonance_curve_integration_method_adj=1 angle integration
c     now it works for ellipse case only                     
c     i_resonance_curve_integration_method_adj=2 rectangle formula
c     for p0_perp integration
c     
c     i_resonance_curve_integration_method_adj=3 trapezoidal formula
c     for p0_perp integration
c     
c------------------------------------------------------------ 
c      i_resonance_curve_integration_method_adj=1 does not work
c      i_resonance_curve_integration_method_adj=2 !rectangle
c      i_resonance_curve_integration_method_adj=3 !trapezoidal 
c      i_resonance_curve_integration_method_adj=4 !adaptive Simpson 
     
      goto (1,2,3,4) i_resonance_curve_integration_method_adj

 1    continue
cSm060327       
      if(dabs(n_parallel).ge.1.d0) goto 3 !to trapezoidal formula 
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     begin
c     This case now works for ellipce case only
c     |N_parallel| <1
c--------------------------------------------------------------------
      rme=9.1094d-28 
      call get_rtem0_from_one(rtem0)
      tem0=temperho(0.d0,1)*rtem0
      vt0=dsqrt(2.d0*tem0*1.6022d-9/rme) !(cm/sec)the  central
                                         ! thermal velocity
      clight=2.99792458d10
      cdvt=clight/vt0
      vmaxdt=dsqrt(rtem0)
      call ec_condh(n_harm_adj,y_loc,n_parallel,vmaxdt/cdvt,ires,v0dc,
     &vmax1dc,vpar0dc,
     &vpar1dc,vpar2dc,thet1,thet2)
     
      power_nharm_adj=0.d0
      CD_nharm_adj=0.d0
   
c-----Subdivide theta range of integration
      dth=(thet2-thet1)/(jmax-1)
      do j=1,jmax
         xint=thet1+(j-1)*dth
         cs=dcos(xint)
         sn=dsin(xint)
         p_par=(vpar0dc-v0dc*cs)       !vper/c
         p_perp=vmax1dc*sn             !vpar/c
         if(p_perp.lt.1.d-12) p_perp=1.d-3
        
         p0_perp=p_perp/dsqrt(b_ratio)
         
         call power_CD_nharm_under_integral_functions(i_calculate_CD,
     &   n_radial0,
     &   p0_perp,n_parallel,n_perp,theta_temperature,b_ratio,
     &   e_x,e_y,e_z,y_loc,temp_kev, 
     &   unorm, sin_trap,
     &   n_harm_adj,f_power_n_harm,f_CD_n_harm)

         f_power_n_harm = f_power_n_harm*p0_perp
         f_CD_n_harm =    f_CD_n_harm*p0_perp

         p_int=1.d0
        
         if((j.ne.jmax).and.(j.ne.1)) then         
           sn=dsin(xint+dth)
           p_perp_pl=vmax1dc*sn             !vper/c
           sn=dsin(xint-dth)
           p_perp_min=vmax1dc*sn            !vper/c  
           h=0.5d0*(p_perp_pl-p_perp_min)
         else
           if(j.eq.1) then
             sn=dsin(xint+dth)
             p_perp_pl=vmax1dc*sn             !vper/c
             h=0.5d0*(p_perp_pl-p_perp)
           else
             !j=jmax
             sn=dsin(xint-dth)
             p_perp_min=vmax1dc*sn             !vper/c
             h=0.5d0*(p_perp-p_perp_min)
           endif
         endif
                 
         power_nharm_adj=power_nharm_adj+h*f_power_n_harm
         CD_nharm_adj=CD_nharm_adj+h*f_CD_n_harm
      enddo 
      goto 20
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     end
c--------------------------------------------------------------------

 2    continue
c-------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     begin
c------------------------------------------------------------------- 
      do j=1,jmax
        p0_perp=h*(j-0.5d0) !p/mc

c        write(*,*)'j,p0_perp',j,p0_perp 
 
       call power_CD_nharm_under_integral_functions(i_calculate_CD,
     &  n_radial0,
     &  p0_perp,n_parallel,n_perp,theta_temperature,b_ratio,
     &  e_x,e_y,e_z,y_loc,temp_kev, 
     &  unorm, sin_trap,
     &  n_harm_adj,f_power_n_harm,f_CD_n_harm)

c        write(*,*)'f_power_n_harm,f_CD_n_harm',
c     &             f_power_n_harm,f_CD_n_harm

        power_nharm_adj=power_nharm_adj+f_power_n_harm*p0_perp*h
        CD_nharm_adj=CD_nharm_adj+f_CD_n_harm*p0_perp*h
               
      enddo 

     
      goto 20
c ------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     end
c-------------------------------------------------------------------

 3    continue
cSm060306
c -------------------------------------------------------------------
c     integration along the resonance curve using trapezoidal
c     formula
c     begin
c--------------------------------------------------------------------

c      write(*,*)'3 jmax',jmax

      do j=0,jmax

c      write(*,*)'j',j

c         p_int=1.d0
c         if((j.eq.1).or.(j.eq.jmax)) p_int=0.5d0

         p0_perp=h*j
        
c         write(*,*)'h,p0_perp=',h,p0_perp
         if(p0_perp.lt.1.d-3) p0_perp=1.d-3
         if(j.eq.jmax) p0_perp=p0_perp-h*1.d-3
c         write(*,*)'p0_perp=',p0_perp

         if (j.eq.1) then
             p_int=0.25d0*h**2
         else
            if(j.eq.jmax) then
               p_int=0.5d0*h*(p0_perp-0.5d0*h)
            else
               p_int=h*p0_perp
            endif
         endif  

c         write(*,*)'bef power_CD_nharm_under_integral_functions p0_perp
c     &   ',p0_perp,'p_int',p_int

         call power_CD_nharm_under_integral_functions(i_calculate_CD,
     &   n_radial0,
     &   p0_perp,n_parallel,n_perp,theta_temperature,b_ratio,
     &   e_x,e_y,e_z,y_loc,temp_kev,
     &   unorm, sin_trap,
     &   n_harm_adj,f_power_n_harm,f_CD_n_harm)

c         write(*,*)'after power_CD_nharm f_power_n_harm,f_CD_n_harm',
c     &              f_power_n_harm,f_CD_n_harm

         power_nharm_adj=power_nharm_adj+p_int*f_power_n_harm
         CD_nharm_adj=CD_nharm_adj+p_int*f_CD_n_harm
    
      enddo
 
      goto20
c-------------------------------------------------------------------
c     end integration along the resonance curve using trapezoidal
c     formula
c     end
c-------------------------------------------------------------------

 4    continue
cSm070417
c -------------------------------------------------------------------
c     integration along the resonance curve using 
c     the adaptive Simpson function
c     begin
c--------------------------------------------------------------------
c       write(*,*)'before calc_integral_array_by_adaptive_simpson'

      call  calc_power_CD_integral_array_by_adaptive_simpson
     &(i_calculate_CD,n_radial0, 
     & b_ratio,y_loc,sin_trap,temp_kev,   
     &e_x,e_y,e_z,
     &n_harm_adj,
     &n_parallel,n_perp,theta_temperature,unorm,umax,
     &n_relt_intgr_adj,
     &p0_perp_max,epsi,power_nharm_adj,CD_nharm_adj)
 
c      write(*,*)'after calc_integral_array_by_adaptive_simpson'
c----------------------------------------------------------------------                  
 20   continue

c-----normalization of the relativisatic Maxwellian distribution  
c     calculation the Mackdonalds function bk2=K_2(theta)*EXP(theta)	
      call besk2as(theta_temperature,bk2)
c-----The Maxwellian function was calculated in subroutine ql_flux_nharm
c     using formula fm=dexp(-theta_relativistic*(gamma-1.d0))
c
c     The derivative d(f_maxw)/d_u was calculated there  using  
c     dfm_du0=-fm*theta_relativistic*u0_dev_c/gamma/(clight/unorm)
c
c     It was used theta_relativistic=theta_temperature in ql_flux_nharm.
c
c     This distribution should be normalized to get the unit density
c     4*pi*integral{0,infinity}[fm*(u/mc)**2*d(u/mc)]=1
c     using formula:
c     fm=theta_temperature*dexp(theta_temperature*(1.d0-gamma))/ 
c            (4.d0*pi*bk2) 
      pi=4.d0*datan(1.d0)
      p_loc = theta_temperature/(4.d0*pi*bk2)
c      write(*,*)'temp_kev,theta_temperature,bk2,p_loc',
c     &           temp_kev,theta_temperature,bk2,p_loc
      power_nharm_adj = power_nharm_adj*p_loc
      CD_nharm_adj = CD_nharm_adj*p_loc

c-----check this normalization
c      dens_test=0.d0
c      h=(umax*unorm/clight)/(jmax-1)
c      write(*,*)'jmax,umax,unorm,clight',jmax,umax,unorm,clight
c      write(*,*)'(umax*unorm/clight),h',(umax*unorm/clight),h
c      do j=1,jmax 
c         p1=h*j          !p/mc
c         gamma=dsqrt(1.d0+p1**2)
c         fm=dexp(-theta_temperature*(gamma-1.d0)) 
c         dens_test=dens_test+fm*p1**2*h 
c      enddo
c      dens_test=4.d0*pi*dens_test
c      write(*,*)'dens_terst',dens_test
c      dens_test=dens_test*p_loc
c      write(*,*)'1 dens_test',dens_test
c      This test gave unit density
c--------------------------------------------
      
 10   continue
      
      return
      end
      
      subroutine power_CD_nharm_under_integral_functions(i_calculate_CD,
     &n_radial0,
     &p0_perp,nll,n_perp,theta_temperature,b_ratio,e_x,e_y,e_z,y_loc,
     &temp_kev,unorm,sin_trap,
     &n_harm_adj,f_power_n_harm,f_CD_n_harm)
c----------------------------------------------------------------
c     Calculates under integral functions
c     f_power_n_harnm,f_CD_n_harm for ADJ calculations
c     at one radial point psis(n_radial0) from chi radial mesh 
c     input
c       i_calculate_CD =1 !calculate f_power_n_harnm,f_CD_n_harm
c                         !for ADJ calculations
c                         !at given adj radial mesh point with number
c                         !n_radial0
c       i_calculate_CD =0 !calculate only f_power_n_harm
c                         !for QL absorption
c                         !at given adj radial point psi or rho
c                         !all input arguments are calculated at this rho
c
c       n_radial0    is a number of the point at the radial chi mesh
c       p0_perp  - momentum divided by (mc)
c       nll     - N_parallel
c       n_perp  - N_perpendicular 
c       theta_temperature   - mc**2/T
c       b_ratio - B(l)/b_min
c       e_x,e_y,e_z  Stix frame complex electric filed polarization
c       y_loc     = |omega_ce|/omega for the electron rest mass. 
c                   It is positive   
c       temp_kev     electron temperature [kev]
c       unorm        normalization velocity sqrt(temp_kev/t) [sm/sec]
c       sin_trap     dsqrt(bmin/bmax)
c       n_harm_adj  - number of the given resonance harmonic
c
c     output
c      f_power_n_harm,f_CD_n_harm  under integral functions
c                   f_power_n_harm/
c          {m*unorm**2*[q^2*E^2/(m^2*unorm)]*[[fm]/unorm]*[omega/unorm]}
c          Here fm=dexp(-theta_relativistic*(gamma-1.d0)) 
c      f_CD_n_harm                 for power and CD ADSJ calculations
c-------------------------------------------------------
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input        
      real*8 p0_perp,nll,n_perp,theta_temperature,b_ratio,
     &y_loc,temp_kev,unorm,sin_trap    

      integer n_harm_adj,n_radial0,i_calculate_CD

      complex*16
     &e_x,e_y,e_z              !Stix frame electric filed polarization

c-----external root_res
      real*8 terp2p_Sm

c-----output
      real*8 f_power_n_harm,f_CD_n_harm
  
c-----local
      real*8 gamma, p_par_rl(2),ps,pi,
     &dgdp_par,p_perp,p0_par,clight,u0,
     &Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &d_g_n_harm_d_u0_parallel,dfm_du0,
     &sin_theta0_res,cos_theta0_res,coef_trap,
     &theta0_res,theta0_odd,
     &d_chi_d_u0,d_chi_d_theta0,chi,
     &p

      integer k,kmax,idm
c     kmax - total number of the resonace condition root p_perp_k(p_perp)
      
    

      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]
      p_perp=p0_perp*dsqrt(b_ratio) 

c      write(*,*)'power_CD_nharm_under_integral_functions p_perp',p_perp
c      write(*,*)'p0_perp,b_ratio',p0_perp,b_ratio
c-------------------------------------------------------------------
c     calculations roots p_par_rl of the resonance condition
c     gamma=N_par*p_par_rl+n_harm*Y
c-------------------------------------------------------------------
      call root_res(p_perp,nll,n_harm_adj,y_loc,kmax,p_par_rl)

c      write(*,*)'p_perp,n_harm_adj,nll,y_loc,kmax,p_par_rl',
c     &p_perp, n_harm_adj,nll,y_loc,kmax,p_par_rl

c-----initialize 
      f_power_n_harm=0.d0
      f_CD_n_harm=0.d0  

cSAP090827
c      idm=imax_chi_a
      idm=imax_chi

      if (kmax.eq.0) goto 20 !no roots
   
      do k=1,kmax ! the sum over all P_parallel_k(perp) roots 

c         write(*,*)'power_CD_nharm_under_integral_functions k', k
c         write(*,*)'p_perp,p_par_rl(k)*p',p_perp,p_par_rl(k)
c         write(*,*)'b_ratio',b_ratio

c         ps=p_perp*p_perp+p_par_rl(k)*p_par_rl(k)       
c         write(*,*)'ps',ps

         p=p_par_rl(k)**2*b_ratio+p_perp**2*(b_ratio-1.d0)

c         p0_par=dsqrt((p_par_rl(k)**2*b_ratio+p_perp**2*(b_ratio-1.d0))
c     &          /b_ratio)

cSAP9090812
         p=p_par_rl(k)**2*b_ratio+p_perp**2*(b_ratio-1.d0)

         if (p.gt.0.d0) then
             p0_par=dsqrt(p/b_ratio)
         else
             p0_par=0.d0
         endif

         ps=p0_perp*p0_perp+p0_par*p0_par

c         write(*,*)'p0_perp,p0_par,ps',p0_perp,p0_par,ps
        
         u0=dsqrt(ps)*clight/unorm  ! normalized electron momentum per mass p/(m*unorm)

c        write(*,*)'power_CD_nharm_under_integral_func u0,umax',
c     &              u0,umax

         if(u0.gt.umax) goto 10

c--------------------------------------------------------------------------
c       calculates: 
c       dfm_du0,      derivative from Maxwellian distribution 
c                     fm=dexp(-theta_relativistic*gamma)
c                     d(f_maxw)/d(u/unorm)=dfm_du0=
c                     =-fm*theta_relativistic*u0_dev_c/gamma
c                              /(clight/unorm)    
c                      u0 [u/unorm]
c
c       Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,!Gamma/(q^2*E^2/(m^2*unorm)
c       d_g_n_harm_d_u0_parallel,                !(dg/du)/(omega/unorm)
c       gamma

c         write(*,*)'before ql_flux_nharm'

         call  ql_flux_nharm(n_harm_adj,
     &   nll,n_perp,u0,t,temp_kev,unorm,
     &   b_ratio,e_x,e_y,e_z,y_loc,
     &   sin_theta0_res,cos_theta0_res,dfm_du0,
     &   Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &   d_g_n_harm_d_u0_parallel,
     &   gamma)


c         write(*,*)'Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0',
c     &              Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0

         if (sin_theta0_res.gt.sin_trap) then
           coef_trap=0.d0
         else     
           coef_trap=1.d0
         endif

c         write(*,*)'lh_ql_flux.f sin_theta0_res,sin_trap,coef_trap',
c     &              sin_theta0_res,sin_trap,coef_trap

         if (dabs(cos_theta0_res).le.1.d0) then
            theta0_res=dacos(cos_theta0_res)
         else

c           write(*,*)'lh_ql_flux.f power_CD_nharm_under_integral_func'
c           write(*,*)'cos_theta0_res',cos_theta0_res

            if (cos_theta0_res.gt.1) then
               theta0_res=0.d0
            else
               theta0_res=pi
            endif
         endif

c         write(*,*)'power_CD_nharm_under_integral_functions theta0_res'
c     &             ,theta0_res
c         write(*,*)'th0a_chi(imax_chi,n_radial0)',
c     &              th0a_chi(imax_chi,n_radial0)
c         write(*,*)'sin_theta0_res,sin_trap,coef_trap',
c     &   sin_theta0_res,sin_trap,coef_trap

cSAP090812
         d_chi_d_u0=0.d0 ! initialize
         d_chi_d_theta0=0.d0 ! initialize
         
         if(i_calculate_CD.eq.0) goto 30 ! for power calculations only
                                         ! it can be problem with pointers for 
                                         ! chi and its derivatives at i_adj=0
                                         ! th0a_chi, ua
                                         ! May be the lines with chi,
                                         ! its derivatives and
                                         ! th0a_chi, ua to put 
                                         ! in the separate subroutine?

         if((theta0_res.le.th0a_chi(imax_chi,n_radial0)).or.
     &     (theta0_res.ge.pi-th0a_chi(imax_chi,n_radial0))) then          
c         if (sin_theta0_res.lt.sin_trap) then
c-----------passing particles
c            write(*,*)'theta0_res,u0',theta0_res,u0
           if (theta0_res.le.th0a_chi(imax_chi,n_radial0)) then

cSAP090130
               if (i_chi_interpolation.eq.1) then
c-----------------spline interpolation

                  d_chi_d_u0=terp2p_Sm(theta0_res,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)
                  d_chi_d_theta0=terp2p_Sm(theta0_res,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)
               endif

cSAP081224
cSAP090130
               if (i_chi_interpolation.eq.2) then
                 call  interpolation_chi(2,n_radial0,u0,theta0_res,
     &           chi,d_chi_d_u0,d_chi_d_theta0)
               endif

           else
               theta0_odd=pi-theta0_res 
cSAP0801208 
c--------------chi function is odd according theta equal to pi/2 
cSAP090130 
               if (i_chi_interpolation.eq.1) then
c-----------------spline interpolation
             
                  d_chi_d_u0=-terp2p_Sm(theta0_odd,u0,   
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)
             
                  d_chi_d_theta0=terp2p_Sm(theta0_odd,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)
               endif
cSAP081224 
cSAP090130
               if (i_chi_interpolation.eq.2) then
                  call  interpolation_chi(2,n_radial0,u0,theta0_odd,
     &            chi,d_chi_d_u0,d_chi_d_theta0)
                  d_chi_d_u0=-d_chi_d_u0
               endif
           endif

cSAP081208
c           f_CD_n_harm=f_CD_n_harm+coef_trap*(u0/gamma)*cos_theta0_res*
c     &               dfm_du0*(
c     &               Gamma_a_n_harm_u0*d_chi_d_u0+
c     &               Gamma_a_n_harm_theta_0*d_chi_d_theta0/u0)/
c     &               dabs(d_g_n_harm_d_u0_parallel)
         endif 
c         write(*,*)'i_chi_interpolation,d_chi_d_u0,d_chi_d_theta0==',
c     +              i_chi_interpolation,d_chi_d_u0,d_chi_d_theta0
cSAP081208
cSAP090727 due to Gamma=-Gamma_a*...
c         f_CD_n_harm=f_CD_n_harm+coef_trap*(u0/gamma)*
         f_CD_n_harm=f_CD_n_harm-coef_trap*(u0/gamma)*
cSAP090727 due to Lambdad=(u_0/gam_0)dabs(cos(theta_0))*tau_ B/L     
c     &              cos_theta0_res*
     &               dabs(cos_theta0_res)*
     &               dfm_du0*(Gamma_a_n_harm_u0*d_chi_d_u0+
     &               Gamma_a_n_harm_theta_0*d_chi_d_theta0/u0)/
     &               dabs(d_g_n_harm_d_u0_parallel)


c         write(*,*)'coef_trap,Gamma_a_n_harm_u0,d_chi_d_u0',
c     &              coef_trap,Gamma_a_n_harm_u0,d_chi_d_u0

c         write(*,*)'Gamma_a_n_harm_theta_0,d_chi_d_theta0,f_CD_n_harm',
c     &              Gamma_a_n_harm_theta_0,d_chi_d_theta0,f_CD_n_harm

 30      continue ! power computing

cSAP090727 due to Gamma=-Gamma_a*...
c         f_power_n_harm=f_power_n_harm+(u0/gamma)**2*
         if (dabs(cos_theta0_res).gt.1.d-12)  then
            f_power_n_harm=f_power_n_harm-(u0/gamma)**2*
cSAP090727 due to Lambdad=(u_0/gam_0)dabs(cos(theta_0))*tau_ B/L     
c     &     cos_theta0_res*
     &                   dabs(cos_theta0_res)*
     &                   Gamma_a_n_harm_u0*dfm_du0/
     &                   dabs(d_g_n_harm_d_u0_parallel)
         endif

c         write(*,*)'u0,gamma,cos_theta0_res, Gamma_a_n_harm_u0',
c     &              u0,gamma,cos_theta0_res, Gamma_a_n_harm_u0

c          write(*,*)'dfm_du0,d_g_n_harm_d_u0_parallel,f_power_n_harm',
c     &              dfm_du0,d_g_n_harm_d_u0_parallel,f_power_n_harm
          
 10      continue       
      enddo !kmax

 20   continue

      return
      end


c     subroutine power_CD_n_calc_theta(n_radial0,
c     &b_ratio,y_loc,temp_kev,e_x,e_y,e_z,unorm,umax,
c     &p0_perp,p0_par,y,nll,n_perp,theta_temperature,n_harm_adj,
c     &f_power_n_harm,f_CD_n_harm)
cc     calculates under integral functions: f_power_n_harm,f_CD_n_harm
cc     input
cc       n_radial0    is a number of radial chi mesh
cc       b_ratio      =B/B_min
cc       y_loc        =|omega_ce|/omega positive value
cc       sin_trap     =sqrt(bmin/bmax) at n_radial0 flux surface
cc       temp_kev,    electron temperatue [KeV]
cc       e_x,e_y,e_z  Stix frame complex electric filed polarization
cc       unorm    = dsqrt(temp_kev/t) normalization velocity [sm/sec]  
cc       umax       maximal momentum normalized by unorm
cc       p_perp  - perpendicular momentum divided by (mc)
cc       p_par     parallel momentum divided by (mc) along resonance curve
cc       y       = omega_ce/omega for the electron rest mass    
cc       nll     - N_parallel
cc       np      - N_perpendicular 
cc       theta_temperature   - mc**2/T
cc       n_harm_adj   number of the given resonance harmonic
cc     output
cc       f_power_n_harm,f_CD_n_harm  under integral functionns
cc-------------------------------------------------------
      !implicit none
cc-----input        
c      real*8 p0_perp,y,nll,n_perp,theta_temperature,p0_par,b_ratio,
c     &y_loc,sin_trap,temp_kev,unorm,umax

c      complex*16 e_x,e_y,e_z
  
c      integer n_harm_adj,n_radial0

cc-----external root_res,     

cc-----output
c      real*8 f_power_n_harm,f_CD_n_harm

cc-----local
c      real*8 gamma, p_par_rl(2),eps,ps,ps_max,pi,
c     +dgdp_par,p_perpc    
        
c      integer k,kmax,i,j
cc     kmax - total number of the resonace condition root p_par_k(p_perp)
      
c      pi=4.d0*datan(1.d0)
cc-----calculations ot the roots p_par_rl of the resonance condition
cc-----gamma=N_par*p_par_rl+nY

c      call root_res(p_perp,nll,n_harm_adj,y_loc,kmax,p_par_rl)

cc-----initialize 
c      f_power_n_harm=0.d0
c      f_CD_n_harm=0.d0

c      if (kmax.eq.0) goto 20
  
c      ps=p_perp*p_perp+p_par_rl(k)*p_par_rk(k)
c      gamma=dsqrt(1.d0+ps)

c      call power_CD_nharm_under_integral_functions(n_radial0,
c     &p0_perp,nll,n_perp,theta_temperature,b_ratio,e_x,e_y,e_z,
c     &y_loc,temp_kev,unorm,sin_trap,
c     &n_harm_adj,f_power_n_harm,f_CD_n_harm)
  
cc-----resonance condition uses the delta function with argument
cc     g_delta=gamma-nll*p_par-n*y
cc     the derivative from this argument d(g_delta)/dp_par=dgdp_par

c      dgdp_par=dabs((p_par-nll*gamma)/gamma)     

c 10   continue
       
c 20   continue

c      return
c      end

      subroutine ql_flux_nharm(n_harm,
     &n_parallel,n_perp,u0,t,temp_kev,unorm,
     &b_ratio,e_x,e_y,e_z,y_loc,
     &sin_theta0_res,cos_theta0_res,dfm_du0,
     &Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &d_g_n_harm_d_u0_parallel,
     &gamma)
c--------------------------------------------------------------------------------
c     calculate quasi linear flux ,
c     derivative d_g_n_harm/d_u0_parallel, gamma
c     for power and CD calculations
c     for one given hyro-harmonic:  n_harm
c-------------------------------------------------------------------------------
      implicit none
c-----input
      integer
     &n_harm               !the number of hyro-harmonic

      real*8
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicula refractive index
     &u0,                  !normalized momentum p/(m*unorm) , 
     &t,                   !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1.
     &temp_kev,            !electron temperature [kev]
     &unorm,               !normalization velocity sqrt(temp_kev/(m*t)) [sm/sec]
     &b_ratio,             !B(l)/B(l=0)=B(l)/B_min 
     &y_loc                !|omega_ce|/omega it is positive
     
      complex*16
     &e_x,e_y,e_z              !Stix frame electric filed polarization

c-----output
      real*8
     & sin_theta0_res,cos_theta0_res,
     &dfm_du0, ! derivative from ther relativistic Maxwellian distribution
               !fm=exp[(-mc^2/T)(gamma-1)]
               !dfm_du0 = d_fm/d(u/c)/(c/unorm)
     &gamma,
     &Gamma_a_n_harm_u0,  ! (Gamma_a*tau_B/L_B)*unorm
                          ! (q^2/m^2)/|E|^2
                          ! /
     &Gamma_a_n_harm_theta_0,
     &d_g_n_harm_d_u0_parallel  !(d_g/d_u_par)*(c/omega)
    
c-----locals
      real*8 cos_theta_res,sin_theta_res,
     &clight,u0_dev_c,u0_parallel_dev_c,
     &u0_perp_dev_c,
     &ksi_res,
     &j_bessel_n_harm,j_bessel_n_harm_minus,j_bessel_n_harm_plus,
     &j_bessel(1:3),d_n_harm,
     &D_n_harm_res,
     &fm,mass_e,k_1_kev,theta_relativistic,pi,
     &sigma,
     & bk2    ! Mackdonalds function bk2=K_2(theta)*EXP(theta)	

      complex*16 theta_n_harm_big_res,i_image

      integer nz,n_harm_abs
c-----external
c     DBESJ

c      write(*,*) 'ql_flux_nharm'

      i_image=dcmplx(0.d0,1.d0)
      clight=2.99792458d10              !light speed [cm/sec]
      pi=4.d0*datan(1.d0)

      u0_dev_c=(u0*unorm)/clight

c      write(*,*)'u0,unorm,clight.u0_dev_c',u0,unorm,clight,u0_dev_c

      gamma=dsqrt(1.d0+u0_dev_c**2)

c      gamma_min=dsqrt(1.d0+(umin*unorm/clight)**2)
c      write(*,*)'ql_flux t,temp_kev,unorm',t,temp_kev,unorm
c-------------------------------------------------------------------------------
c     pitch angle at the resonance curve at the local poloidal point
c-------------------------------------------------------------------------------
c      write(*,*)'lh_ql_flux.f gamma,n_harm,y_loc,n_parallel,u0_dev_c',
c     &gamma,n_harm,y_loc,n_parallel,u0_dev_c

      cos_theta_res=(gamma-n_harm*y_loc)/(n_parallel*u0_dev_c) !

      if((1.d0- cos_theta_res**2).lt.0.d0) then 
         sin_theta_res=0.d0
      else
         sin_theta_res=dsqrt(1.d0- cos_theta_res**2)
      endif

c      write(*,*)'ql_flux u0,n_parallel,cos_theta_res,sin_theta_res',
c     & u0,n_parallel,cos_theta_res,sin_theta_res
c-------------------------------------------------------------------------------
c     pitch angle at the resonance curve at the midplane b=bmin
c-------------------------------------------------------------------------------
c      write(*,*)'ql_flux b_ratio',b_ratio

      sin_theta0_res=sin_theta_res/dsqrt(b_ratio)
      if (cos_theta_res.ge.0d0) then
        if ((1.d0-sin_theta0_res**2).lt.0.d0) then
            cos_theta0_res=0.d0
        else
            cos_theta0_res=dsqrt(1.d0-sin_theta0_res**2)
        endif
      else 
         if ((1.d0-sin_theta0_res**2).lt.0.d0) then
            cos_theta0_res=0.d0
         else
            cos_theta0_res=-dsqrt(1.d0-sin_theta0_res**2)
         endif
      endif

      u0_parallel_dev_c=u0_dev_c*cos_theta0_res
      u0_perp_dev_c=u0_dev_c*sin_theta0_res

c      write(*,*)'ql_flux cos_theta0_res,sin_theta0_res',
c     & cos_theta0_res,sin_theta0_res
c--------------------------------------------------------------------------------
c     Bessel function agrgument
c-------------------------------------------------------------------------------
      ksi_res=n_perp*u0_dev_c*sin_theta_res/y_loc
c      write(*,*)'ql_flux_nharm n_perp,u0_dev_c,sin_theta_res,y_loc,
c     &ksi_res',n_perp,u0_dev_c,sin_theta_res,y_loc,ksi_res
c      write(*,*)'u0_perp_dev_c*clight/unorm',
c     &u0_perp_dev_c*clight/unorm
c------------------------------------------------------------------------------
c      nz controles bessel function calculations
c      nz should be =0. If nz=1 sum bessel function(dBESJ) will be J=0

c      write(*,*)'n_harm',n_harm
      d_n_harm=n_harm
c        write(*,*)'n_harm,d_n_harm',n_harm,d_n_harm
      if (n_harm.ge.1) then

c         write(*,*)'befor j_bessel: ksi_res,d_n_harm', 
c     &                              ksi_res,d_n_harm

         call DBESJ( ksi_res,d_n_harm,3,j_bessel,nz)

c        write(*,*)'j_bessel,nz',j_bessel,nz

         j_bessel_n_harm_minus=j_bessel(1)
         j_bessel_n_harm=j_bessel(2)
         j_bessel_n_harm_plus=j_bessel(3)
      else
         if(n_harm.eq.0)then
           call DBESJ( ksi_res,1.d0,1,j_bessel_n_harm_plus,nz)
c          write(*,*)'ql_flux after dbesj 1 nz,j_bessel_n_harm_plus',
c          nz,j_bessel_n_harm_plus',
           call DBESJ( ksi_res,0.d0,1,j_bessel_n_harm,nz)
           j_bessel_n_harm_minus=-j_bessel_n_harm_plus
         else
c----------n_harm.lt.0
           n_harm_abs=-n_harm
           call DBESJ( ksi_res,dfloat(n_harm_abs),3,j_bessel,nz)  
           j_bessel_n_harm_minus=(-1)**(n_harm_abs+1)*j_bessel(3)
           j_bessel_n_harm=(-1)**n_harm_abs*j_bessel(2)
           j_bessel_n_harm_plus=(-1)**(n_harm_abs-1)*j_bessel(1)
         endif
      endif

cSAP080411     
c      theta_n_harm_big_res=0.5d0*(e_x-i_image*e_y)*j_bessel_n_harm_plus
c     &                   +0.5d0*(e_x+i_image*e_y)*j_bessel_n_harm_minus
      theta_n_harm_big_res=0.5d0*(e_x-i_image*e_y)*j_bessel_n_harm_minus
     &                   +0.5d0*(e_x+i_image*e_y)*j_bessel_n_harm_plus
     &                 +e_z*j_bessel_n_harm*cos_theta_res/sin_theta_res


cSAP081216
c       theta_n_harm_big_res=
c     &                 +e_z*j_bessel_n_harm*cos_theta_res/sin_theta_res
c      test for n_harm=0 case
c      theta_n_harm_big_res=i_image*e_y*j_bessel_n_harm_plus
c     &                 +e_z*j_bessel_n_harm*cos_theta_res/sin_theta_res

c      write(*,*)'j_bessel_n_harm_minus,j_bessel_n_harm,
c     &j_bessel_n_harm_plus',
c     &j_bessel_n_harm_minus,j_bessel_n_harm,
c     &j_bessel_n_harm_plus

c      write(*,*)'cos_theta_res,cos_theta_res/sin_theta_res',
c     &           cos_theta_res,cos_theta_res/sin_theta_res

c      write(*,*)'ql_flux e_x,e_y,e_z',e_x,e_y,e_z 
c      write(*,*)'theta_n_harm_big_res',theta_n_harm_big_res

      D_n_harm_res=0.5d0*pi*
     &         theta_n_harm_big_res*dconjg(theta_n_harm_big_res)
 
c     write(*,*)' D_n_harm_res', D_n_harm_res
  
c-----Gamma_a are normalized to (q^2*E^2)/(m^2*unorm)
      Gamma_a_n_harm_u0= gamma/(u0*dabs(cos_theta_res))*D_n_harm_res*
     &     ((n_parallel*u0_dev_c*sin_theta_res*cos_theta_res+
     &       n_harm*y_loc*sin_theta_res)/gamma)**2

c      write(*,*)'Gamma_a_u0',Gamma_a_u0

      Gamma_a_n_harm_theta_0= gamma/(u0*dabs(cos_theta_res))*
     & D_n_harm_res*
     & (sin_theta0_res/cos_theta0_res)/(sin_theta_res/cos_theta_res)*
     & (-((n_parallel*u0_dev_c*sin_theta_res/gamma)**2-  
     &    (n_harm*y_loc/gamma)**2)*sin_theta_res*cos_theta_res+
     &   (n_parallel*u0_dev_c*sin_theta_res/gamma)*(n_harm*y_loc/gamma)
     &    *(cos_theta_res**2-sin_theta_res**2))

      sigma=1.d0
      if (u0_parallel_dev_c.lt.0.d0) sigma=-1.d0
   
c-----d(g)/d(u0_parallel) is normalized to (omega/unorm)
      d_g_n_harm_d_u0_parallel=(
     &( sigma*n_parallel*
     &  dsqrt(u0_parallel_dev_c**2-(b_ratio-1.d0)*u0_perp_dev_c**2)/
     &  gamma**2+n_harm*y_loc/gamma**2)*(u0_parallel_dev_c/gamma)-
     &sigma*n_parallel*u0_parallel_dev_c/(gamma*
     & dsqrt(u0_parallel_dev_c**2-u0_perp_dev_c**2*(b_ratio-1.d0)))
     & )/(clight/unorm)
    

      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)

      theta_relativistic=mass_e*clight**2/(k_1_kev*temp_kev)
  
cSAP070831
c      fm=dexp(-theta_relativistic*gamma)
c      write(*,*)'temp_kev,theta_relativistic,fm',
c     &           temp_kev,theta_relativistic,fm
     
c      fm=dexp(-theta_relativistic*(gamma-gamma_min))
c  
c------relativistic Maxwellian distribution 
       fm=dexp(-theta_relativistic*(gamma-1.d0)) 
    
cSAP081201
c      dfm_du0=-fm*theta_relativistic*u0_dev_c/gamma
c------d(f_maxw)/d(u/unorm)
c      the derivative d(fm)/d(u0) is normailized to [fm]/unorm
       dfm_du0=-fm*theta_relativistic*u0_dev_c/gamma
     &  /(clight/unorm)
c      write(*,*)'(gamma-gamma_min),fm ',(gamma-gamma_min),fm 
       
      return
      end


c     Modificated method for vector function.
c     
c     Original description and functions are from:
c
C PAGE 187-190: NUMERICAL MATHEMATICS AND COMPUTING, CHENEY/KINCAID, 1985
C
C FILE: SIMP.FOR
C
C ADAPTIVE SCHEME FOR SIMPSON'S RULE (SIMP,ASMP,PUSH,POP,FCN)
c============================================================================

      subroutine fcn_adj_vector(p0_perp,fcn)
c---------------------------------------------------------------------------
c     culculates under integral function for 
c     power and CD calulations using 
c---------------------------------------------------------------------------

      implicit none
      integer k_max !max number of elements
      parameter (k_max=2) 
c----------------------------------------------------------------
c     input
c----------------------------------------------------------------- 
      real*8 p0_perp !under integral function argument
c-----------------------------------------------------------------
c     output
c-----------------------------------------------------------------
      real*8 fcn(k_max) !array of under integral functions values
c-----------------------------------------------------------------
      integer n_harm_adj_l,n_radial0_l,i_calculate_CD_l

      real*8 y_l,n_parallel_l,n_perp_l,theta_temperature_l,b_ratio_l,
     &temp_kev_l,unorm_l,sin_trap_l

      complex*16 e_x_l,e_y_l,e_z_l

      common /fcn_adj_input/
     &y_l,n_parallel_l,n_perp_l,theta_temperature_l,
     &temp_kev_l,unorm_l,sin_trap_l,b_ratio_l,
     &e_x_l,e_y_l,e_z_l,
     &n_harm_adj_l,n_radial0_l,i_calculate_CD_l

     
c-----------------------------------------------------------------
c     locals
c----------------------------------------------------------------
      real*8 
     &f_power_n_harm,f_CD_n_harm ! under integral functions


      call power_CD_nharm_under_integral_functions(
     &i_calculate_CD_l,
     &n_radial0_l,
     &p0_perp,n_parallel_l,n_perp_l,theta_temperature_l,b_ratio_l,
     &e_x_l,e_y_l,e_z_l,y_l,
     &temp_kev_l,unorm_l,sin_trap_l,
     &n_harm_adj_l,f_power_n_harm,f_CD_n_harm)

      FCN(1) =f_power_n_harm*p0_perp
      if (i_calculate_CD_l.eq.1) then
         FCN(2) =f_CD_n_harm*p0_perp 
      else
         FCN(2) =0.d0
      endif

      return
      end


 
   
      subroutine  calc_power_CD_integral_array_by_adaptive_simpson
     &(i_calculate_CD,n_radial0, 
     & b_ratio,y_loc,sin_trap,temp_kev,   
     &e_x,e_y,e_z,
     &n_harm_adj,
     &n_parallel,n_perp,theta_temperature,unorm,umax,
     &n_relt_intgr_adj,
     &p0_perp_max,epsi,power_nharm_adj,CD_nharm_adj)      
c------------------------------------------------------------------------
c     calculate intergals power and CD (power_nharm_adj,CD_nharm_adj)    
c     using ADJ and Sypson adaptive method
c     for the EC harmonic with number 'n_harm_adj'
c-------------------------------------------------------------------------
       implicit none
c------------------------------------------------------------------------
c      input
c--------------------------------------------------------------------------
c        i_calculate_CD   ! =1 calculate power and CD 
c                         ! at given adj radial mesh point
c                         ! with number n_radial0
c                         !
c                         ! =0, calculate power only (no CD) 
c                         ! using input arguments computed at
c                         ! given space point
c                         ! In this case n_radial0 can be arbitary
c
c       n_radial0    is a number of radial chi mesh
c       b_ratio      =B/B_min
c       y_loc        =|omega_ce|/omega for the electron rest mass 
c                     It has positive value
c       sin_trap     =sqrt(bmin/bmax) at n_radial0 flux surface
c       temp_kev,    electron temperature [KeV] 
c       e_x,e_y,e_z  Stix frame complex electric filed polarization
c       n_harm_adj        - EC harmonic number
c       n_parallel          N_parallel
c       n_perp      N_perpendicular      
c       theta_temperature    = mc**2/T
c       unorm    = dsqrt(temp_kev/t) normalization velocity [sm/sec]  
c       umax       maximal momentum normalized by unorm
c       n_relt_intgr_adj the number of points for the integration over p0_perp
c       p0_perp_max is  the maximal boundary of integration
c       epsi is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------------------
      integer n_radial0,n_harm_adj, n_relt_intgr_adj,
     & i_calculate_CD

      real*8 b_ratio,y_loc,sin_trap,temp_kev,n_parallel,n_perp,    
     &theta_temperature,unorm,umax,epsi,p0_perp_max
      
      complex*16 e_x,e_y,e_z

c------------------------------------------------------------------------
c     output
c------------------------------------------------------------------------
      real*8 power_nharm_adj,CD_nharm_adj
c-------------------------------------------------------------------------
c     locals
c--------------------------------------------------------------------------
      integer k_max !max number of elements
      parameter (k_max=2)
      integer ist,ifs,k   
c      PARAMETER (IST=5,IFS=16)
c      PARAMETER (IST=9,IFS=256)  
c       PARAMETER (IST=13,IFS=10048)
c       PARAMETER (IST=16,IFS=100048)
c       PARAMETER (IST=17,IFS=200048)!070723-old
      PARAMETER (IST=20,IFS=200048)

      REAL*8 STACK(IST,3),FSTACK(IFS,3)
      real*8 a,b,sum(k_max)
      integer lvmax,iflag,i
      real*8 fcn(k_max)
    
      integer n_harm_adj_l,n_radial0_l,i_calculate_CD_l

      real*8 y_l,n_parallel_l,n_perp_l,theta_temperature_l,b_ratio_l,
     &temp_kev_l,unorm_l,sin_trap_l

      complex*16 e_x_l,e_y_l,e_z_l
     
      common /fcn_adj_input/
     &y_l,n_parallel_l,n_perp_l,theta_temperature_l,
     &temp_kev_l,unorm_l,sin_trap_l,b_ratio_l,
     &e_x_l,e_y_l,e_z_l,
     &n_harm_adj_l,n_radial0_l,i_calculate_CD_l
 
      EXTERNAL FCN_adj_vector
c      DATA EPSI/5.0d-5/, LVMAX/4/
c      DATA EPSI/5.0d-5/, LVMAX/7/
c      DATA EPSI/5.0d-5/, LVMAX/15/
c      DATA EPSI/5.0d-5/, LVMAX/12/
c      DATA EPSI/5.0d-5/, LVMAX/16/
c       DATA EPSI/5.0d-6/, LVMAX/16/  !070722-new
c       DATA EPSI/1.0d-6/, LVMAX/16/  !070722-new

c       DATA EPSI/1.0d-3/, LVMAX/16/ !070722-old
c      DATA EPSI/1.0d-3/, LVMAX/19/
      DATA LVMAX/16/  !070722-new



c-----------------------------------------------------------------
c     set common /fcn_adj_input/
c-----------------------------------------------------------------
      y_l=y_loc
      n_parallel_l=n_parallel
      n_perp_l=n_perp
      theta_temperature_l=theta_temperature
      n_harm_adj_l=n_harm_adj
      n_radial0_l=n_radial0     
      temp_kev_l=temp_kev
      unorm_l=unorm
      sin_trap_l=sin_trap
      b_ratio_l=b_ratio
      e_x_l=e_x
      e_y_l=e_y
      e_z_l=e_z
      n_harm_adj_l=n_harm_adj
      n_radial0_l=n_radial0
      i_calculate_CD_l=i_calculate_CD
c-----------------------------------------------------------------
c     integration boundaries
c------------------------------------------------------------------
c      A = 0.0d0
      A = 1.d-5
c      A = 1.d-3

      B = p0_perp_max-1.d-5
     
      CALL ASMP(FCN_adj_vector,k_max,fcn,
     &A,B,EPSI,LVMAX,STACK,IST,FSTACK,IFS,
     &SUM,IFLAG) 



      IF(IFLAG .EQ. 0) THEN 
c$$$        PRINT 4   
      ELSE
c        PRINT 5   
        write(*,*)  '          WITH BAD SUBINTERVALS:' 
        DO 2 I=1,IFLAG      
c          PRINT 6,FSTACK(I,1),FSTACK(I,2),INT(FSTACK(I,3))
          write(*,*) FSTACK(I,1),FSTACK(I,2),INT(FSTACK(I,3))
    2   CONTINUE  
      END IF      
   3  FORMAT(//5X,'APPROXIMATE INTEGRAL =',E22.14/)   
   4  FORMAT(10X,'WITH NO BAD SUBINTERVALS')    
   5  FORMAT(10X,'WITH BAD SUBINTERVALS:')      
   6  FORMAT(10X,'[',F10.5,',',F10.5,']',2X,'LEVEL =',I5) 

c-------------------------------------------------------

      power_nharm_adj=Sum(1) 
      if (i_calculate_CD.eq.1) then
         CD_nharm_adj=Sum(2)
      else 
         CD_nharm_adj=0.d0
      endif
     
      return
      END 
     
      subroutine CD_adj_efficiency(n_parallel,n_perp,
     &psi_loc,theta_pol,e_x,e_y,e_z,
     &cd_efficiency)
c-----------------------------------------------------------------
c     calculate CD efficiency using adj chi function
c-----------------------------------------------------------------       
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input
      real*8      
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicular refractive index
     &psi_loc,             !polidal flux 
     &theta_pol,           !poloidal angle
     &length               !field line length normalized to  
c     &t,                  !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1. 
      complex*16 
     &e_x,
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1
c-----output
      real*8 
     &cd_efficiency        
c-----local
      integer  n_radial_minus,n_radial_plus
      real*8 weight_radial,cd_efficiency_minus,
     &cd_efficiency_plus,unorm_minus,unorm_plus
     
 
c      write(*,*)' CD_adj_efficiency e_x,e_y,e_z', e_x,e_y,e_z
c      write(*,*)' CD_adj_efficiency before find_psis_points'
c      write(*,*)' n_radial_minus,n_radial_plus,weight_radial',
c     &             n_radial_minus,n_radial_plus,weight_radial
  
      call  find_psis_points(psi_loc,
     & n_radial_minus,n_radial_plus,weight_radial)

c      write(*,*)' CD_adj_efficiency after find_psis_points'
cSAP110124 for test only of the first point at rho=0 
c      if (n_radial_minus.eq.1)  n_radial_minus=2
c      if (n_radial_plus.eq.1) n_radial_plus=2
cend test

      if((n_radial_minus.lt.npsi0).and.(n_radial_plus.gt.1)) then

         
         call CD_adj_efficiency_1(n_radial_minus,
     &   theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &   unorm_minus, 
     &   cd_efficiency_minus)
         write(*,*)'eff cd_efficiency_minus',cd_efficiency_minus

          call CD_adj_efficiency_1(n_radial_plus,
     &   theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &   unorm_plus, 
     &   cd_efficiency_plus)
         write(*,*)'eff cd_efficiency_plus',cd_efficiency_plus


         cd_efficiency=cd_efficiency_plus*weight_radial+
     &   cd_efficiency_minus*(1.d0-weight_radial)

      else
        if (n_radial_plus.eq.1) then  
         
           call CD_adj_efficiency_1(1,
     &     theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &     unorm_plus, 
     &     cd_efficiency)

            write(*,*)'eff cd_efficiency 1',cd_efficiency
        else
           if (n_radial_minus.eq.npsi0) then
        
              call CD_adj_efficiency_1(npsi0,
     &        theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &        unorm_minus, 
     &        cd_efficiency)

              write(*,*)'eff cd_efficiency_ npsi0',cd_efficiency

           endif
        endif
      endif

      write(*,*)'in CD_adj_efficiency cd_efficiency=',
     &cd_efficiency

      return
      end

c######################################################1
c for test cases
c######################################################
c=============================================================================
c     for u_perp integration test
c=============================================================================    
    
      subroutine CD_adj_efficiency_1_test(n_radial0,
     &theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &unorm, 
     &cd_efficiency)
c------------------------------------------------------
c     calculate CD efficiency using adj chi function
c     at one radial point psis(n_radial0) from chi radial mesh 
c------------------------------------------------------
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input
      integer
     & n_radial0           ! is a number of radial ch mesh

      real*8  
     &theta_pol,           !poloidal angle [radian]    
     &n_parallel,          !parallel refractive index
     &n_perp               !perpendicular refractive index
    
      complex*16 
     &e_x,
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1
c-----output
      real*8
     &unorm,               !sqrt(T/m) cm/sec
     
     &cd_efficiency

c-----external
      real*8 temperho,densrho,
     &bmax_psi,bmin_psi,psi_rho,rhopsi,y,b_d_bmin

c-----local
      integer n_radial,n_harm_adj,
c     &i_resonance_curve_integration_method
     &i_calculate_CD     ! =1 calculate power and CD 
                         ! at given adj radial mesh point
                         ! with number n_radial0 
      real*8 u0,     
     &rho_small,           !small radius
     &temp_kev,            !electron temperatura [KeV]
     &dense,               !electron dencity 10**13/cm**3     
     &b_ratio,             !B/B_min
     &pi,bmin,bmax,deltb,psi,th0max,sin_trap,
     &d_chi_du0,d_chi_dtheta0,u0max,chi,
     &clight,             !light speed [cm/sec]
     &arg1,cln,
     &z,r,phi,            !bellow phi is set arbitratry phi=0
     &y_loc,              !electron y at point (z,r,phi)     
     &power_nharm_adj,CD_nharm_adj, !under integral functions
     &pow_dens,            !Power density normalized to m*u_norm**2*n/tau_n 
     &cd_dens              !current drive density normalized to q*u_norm*n
                           !Here n is the density  
      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]

      write(*,*)'CD_adj_efficiency_1_test n_radial0',n_radial0

c      do n_radial=1,npsi0
c         write(*,*)'psis',psis(n_radial)
c         rho_small=rhopsi(psis(n_radial))
c         write(*,*)'rho_small',rho_small
c      enddo

      rho_small=rhopsi(psis(n_radial0))
      psi=psis(n_radial0)

      write(*,*)'n_radial0=',n_radial0

      write(*,*)'psis(n_radial0),rho_small',psis(n_radial0),rho_small

      temp_kev=temperho(rho_small,1)
      dense=densrho(rho_small,1) 
      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
      bmax= bmax_psi(psi)
      bmin= bmin_psi(psi)
	     
      deltb = (bmax-bmin)/bmin
c      write(*,*)'bmax,mbin,deltb',bmax,bmin,deltb   
      th0max = atan2(1.0d0,sqrt(deltb))
c      write(*,*)'th0max',th0max  mass_e=9.1094d-28    !electron rest mass (g)
c      k_1_kev=1.6022d-9                                !egrs in 1 KeV      (erg)
      sin_trap=dsqrt(bmin/bmax)

      write(*,*)'sin_trap=',sin_trap,'dsin(th0max)=',dsin(th0max)

      pow_dens=0.d0
      cd_dens=0.d0

c      umin=clight/dsqrt(n_parallel**2-1.d0)/unorm ! normalized umin      
c      write(*,*)'umin',umin,'umin/dsqrt(2.d0)',umin/dsqrt(2.d0)
      
      b_ratio=b_d_bmin(n_radial0,theta_pol)
      write(*,*)'n_radial0,theta_pol,b_ratio',
     &n_radial0,theta_pol,b_ratio

      phi=0.d0
      call zr_psith(psi,theta_pol,z,r)
      y_loc=y(z,r,phi,1)

c      i_resonance_curve_integration_method_adj=3
    
      do n_harm_adj=n_harm_adj_min,n_harm_adj_max
c-------------------------------------------------------------------
c       the loop over used gyro-harmonics 
c----------------------------------------------------------------- 
       write(*,*)'CD_adj_efficiency_1 n_harm_adj=',n_harm_adj,
     & 'i_resonance_curve_integration_method_adj=',
     & i_resonance_curve_integration_method_adj
  
c        i_calculate_CD   ! =1 calculate power and CD 
c                         ! at given adj radial mesh point
c                         ! with number n_radial0 
c        It is not used in test version                          
        call intgr_relt_power_cd_adj_test(
     &  n_radial0,b_ratio,
     &  y_loc,sin_trap,temp_kev, 
     &  e_x,e_y,e_z,  
     &  n_harm_adj,n_parallel,n_perp, 
     &  unorm,umax,
     &  n_relt_intgr_adj,i_resonance_curve_integration_method_adj,
     &  epsi_adj,
     &  power_nharm_adj,CD_nharm_adj)
 
        pow_dens=pow_dens +power_nharm_adj
        cd_dens=cd_dens+CD_nharm_adj
 
      enddo !n_harm

      pow_dens=2.d0*pi*pow_dens
      cd_dens=-2.d0*pi*cd_dens
      write(*,*)'pow_dens',pow_dens
      write(*,*)'cd_dens',cd_dens
     
      if(pow_dens.ne.0.d0) then
        cd_efficiency=cd_dens/pow_dens
      else
        cd_efficiency=0.d0
      endif

      write(*,*)'in CD_adj_efficiency_1 cd_efficiency=',
     &cd_efficiency
         
c-----cln -Coulomb Ln
      arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense)	!temp KeV, dense 10**13/cm**3
      cln=24.d0-dlog(arg1)
      cd_efficiency=cd_efficiency*
     &3.84d0*temp_kev/(cln*dense)	       !  (A/m**2)/(W/m**3)
c      write(*,*)'temp_kev,dense,cln,cd_efficiency',
c     &temp_kev,dense,cln,cd_efficiency
      cd_efficiency=cd_efficiency*1.d-5  !  A/cm**2/(erg/(sec*cm**3)
      write(*,*)'temp_kev,dense,cln,cd_efficiency',
     &temp_kev,dense,cln,cd_efficiency
      return
      end

     

      subroutine intgr_relt_power_cd_adj_test(n_radial0,b_ratio,y_loc,
     &sin_trap,temp_kev, 
     &e_x,e_y,e_z,  
     &n_harm_adj,n_parallel,n_perp,
     &unorm,umax,
     &n_relt_intgr_adj,i_resonance_curve_integration_method_adj,epsi,
     &power_nharm_adj,CD_nharm_adj)
c-----------------------------------------------------------------
c     calculates power density  power_nharm and 
c     CD density  CD_nharm_adj for one harmonic using ADjJ function
c     as  ther integral for the relativistic electrons 
c     integral(0<=p0_perp<=p0_perp_max)f_n_harm(p0_perp)
c     f_n_harm={sum_k=1,2,f_nk(p0_perp)
c     for the EC harmonic with number 'n_harm_adj'
c----------------------------------------------------------------- 
c     input
c       n_radial0    is a number of radial chi mesh
c       b_ratio      =B/B_min
c       y_loc        =|omega_ce|/omega positive value
c       sin_trap     =sqrt(bmin/bmax) at n_radial0 flux surface
c       temp_kev     electron temperature [KeV]
c       e_x,e_y,e_z  Stix frame complex electric filed polarization
c       n_harm_adj        - EC harmonic number
c       n_parallel          N_parallel
c       n_perp      N_perpendicular
c       y_loc    = |omega_ce|/omega for the electron rest mass
c       unorm    = dsqrt(temp_kev/t) normalization velocity [sm/sec]  
c       umax       maximal momentum normalized by unorm
c       n_relt_intgr_adj the number of points for the integration over p0_perp
c----------------------------------------------------------
c     i_resonance_curve_interation_method_adj
c
c     i_resonance_curve_integration_method_adj=1 angle integration
c     now it works for ellipse case only                             c
c
c     i_resonance_curve_integration_method_adj=2 rectangle formula
c     for p0_perp integration
c     
c     i_resonance_curve_integration_method_adj=3 trapezoidal formula
c     for p0_perp integration
c
c     i_resonance_curve_integration_method_adj=4 !adaptive Simpson 
c     for p0_perp integration
c
c       epsi is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------ 
c     
c     output
c       power_nharm_adj,CD_nharm_adj =integral 
c       over 0<=p0_perp=<p0_perp_max      
c----------------------------------------------------------------    
c      The integration method is specified by the variable
c      i_resonance_curve_integration_method_adj.
c      Now this variable is set inside this subroutine.
c----------------------------------------------------------------

      implicit none
c-----input
      real*8 b_ratio,n_parallel,n_perp,y_loc,sin_trap,temp_kev,
     &epsi,
     &unorm,umax 

      integer n_harm_adj,n_relt_intgr_adj,
     &i_resonance_curve_integration_method_adj,n_radial0

      complex*16
     &e_x,e_y,e_z              !Stix frame electric filed polarization
c-----output
      real*8 power_nharm_adj,CD_nharm_adj

c-----local
      real*8 h,p0_perp,p_perp,p_t,clight,p_int,
     &p_perp_max,  ! max value of the perpendicular momentum divided by mc
                    ! on the resonance ellipse
     &p0_perp_max, ! max value of the perpendicular momentum p0_perp divided by mc                   
     &f_power_n_harm,f_CD_n_harm, !under integral functions
                                   !for ADJ calculations
     &theta_temperature,           !mc**2/T
     &mass_e,             !electron rest mass (g)
     &k_1_kev            !egrs in 1 KeV      (erg)

c------------------------------------------------------------------
c     for integration along ellipse by the angle
c------------------------------------------------------------------
      real*8 rme,rtem0,xint,cs,sn,thet1,thet2,p_par,
     &p_par_pl,p_perp_min,p_perp_pl,v0dc,vmax1dc,vpar0dc,
     &vpar1dc,vpar2dc,dth,vt0,cdvt,tem0,vmaxdt
c-----------------------------------------------------------------

      integer j,jmax,
     &ires   ! ires     =0 no resonance,=1 ellipse,=2 parabola,=3 hyperbola  

c-----to plot QL flux versus p_perp
      real*8 Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &D_n_harm_res,d_chi_d_u0,d_chi_d_theta0,d_g_n_harm_d_u0_parallel
      integer kmax
      
      real*8, dimension(n_relt_intgr_adj) :: Gamma_a_n_harm_u0_ar
      real*8, dimension(n_relt_intgr_adj) :: Gamma_a_n_harm_theta_0_ar
      real*8, dimension(n_relt_intgr_adj) :: D_n_harm_res_ar
      real*8, dimension(n_relt_intgr_adj) :: p_perp_d_vt_ar 
      real*8, dimension(n_relt_intgr_adj) :: p_perp_d_vt_ar_s
      real*8, dimension(n_relt_intgr_adj) :: p_perp_d_vt_ar_4
      real*8, dimension(n_relt_intgr_adj) :: power_spectrum
      real*8, dimension(n_relt_intgr_adj) :: cur_spectrum
      real*8, dimension(n_relt_intgr_adj) :: d_chi_d_u0_ar
      real*8, dimension(n_relt_intgr_adj) :: d_chi_d_theta0_ar
      real*8, dimension(n_relt_intgr_adj) :: d_g_n_harm_d_u0_parallel_ar
      real*8, dimension(n_relt_intgr_adj) :: kmax_ar
c-----end to plot QL flux versus p_perp

c     external
c     ec_cond,power_CD_nharm_under_integral_functions
      real*8 temperho
 
c      write(*,*)'intgr_relt_power_cd_adj'
c      write(*,*)'i_resonance_curve_integration_method_adj',
c     & i_resonance_curve_integration_method_adj


      jmax=n_relt_intgr_adj

      mass_e=9.1094d-28         !electron rest mass (g)
      clight=2.99792458d10      !light speed [sm/sec]
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)
      theta_temperature=mass_e*clight**2/(k_1_kev*temp_kev) ! m*c**2/T_e

      call ec_cond(n_harm_adj,y_loc,n_parallel,ires,p_perp_max)

cSAP090126
      if (ires.eq.0) then !no resonance
         power_nharm_adj=0.d0
         CD_nharm_adj=0.d0
        return
      endif

      if (ires.eq.1) then
c-------ellipse case
cSAP090126
c        p_perp_max=(n_parallel**2+(n_harm_adj*y_loc)**2-1.d0)/
        p_perp_max=dsqrt((n_parallel**2+(n_harm_adj*y_loc)**2-1.d0))/
     &             (1.d0-n_parallel**2) !p/mc
cSAP081115
c        p_perp_max=p_perp_max*clight/unorm     !maximal u_perp normilized to sqrt(T/m)
      else
c-------ires=2,3 hyperbola or parabola cases
c        p_perp_max=umax 
cSAP081115
        p_perp_max=umax*unorm/clight !p/mc        !
      endif

c      write(*,*)'ires=',ires,'p_perp_max=', p_perp_max
 
      p0_perp_max=p_perp_max/(b_ratio)                

c      write(*,*)'b_ratio,p0_perp_max=',b_ratio,p0_perp_max

      if (p0_perp_max.gt.umax) p0_perp_max=umax  

c      write(*,*)'umax,p0_perp_max=',umax,p0_perp_max

      h=p0_perp_max/(n_relt_intgr_adj)   ! step of integration over p_perp

c-----calculations of the integrals over p_perp

      power_nharm_adj=0.d0
      CD_nharm_adj=0.d0

c----------------------------------------------------------
c     i_resonance_curve_interation_method_adj
c
c     i_resonance_curve_integration_method_adj=1 angle integration
c     now it works for ellipse case only                             c
c     i_resonance_curve_integration_method_adj=2 rectangle formula
c     for p0_perp integration
c     
c     i_resonance_curve_integration_method_adj=3 trapezoidal formula
c     for p0_perp integration
c     
c------------------------------------------------------------ 
c      i_resonance_curve_integration_method_adj=1
c      i_resonance_curve_integration_method_adj=2
c      i_resonance_curve_integration_method_adj=3 !trapezoidal 
c      i_resonance_curve_integration_method_adj=4 !adaptive Simpson 
     
      goto (1,2,3,4) i_resonance_curve_integration_method_adj

 1    continue
cSm060327       
      if(dabs(n_parallel).ge.1.d0) goto 3 !to trapezoidal formula 
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     begin
c     This case now works for ellipce case only
c     |N_parallel| <1
c--------------------------------------------------------------------
      rme=9.1094d-28 
      call get_rtem0_from_one(rtem0)
      tem0=temperho(0.d0,1)*rtem0
      vt0=dsqrt(2.d0*tem0*1.6022d-9/rme) !(cm/sec)the  central
                                         ! thermal velocity
      clight=2.99792458d10
      cdvt=clight/vt0
      vmaxdt=dsqrt(rtem0)
      call ec_condh(n_harm_adj,y_loc,n_parallel,vmaxdt/cdvt,ires,v0dc,
     &vmax1dc,vpar0dc,
     &vpar1dc,vpar2dc,thet1,thet2)
     
      power_nharm_adj=0.d0
      CD_nharm_adj=0.d0
   
c-----Subdivide theta range of integration
      dth=(thet2-thet1)/(jmax-1)
      do j=1,jmax
         xint=thet1+(j-1)*dth
         cs=dcos(xint)
         sn=dsin(xint)
         p_par=(vpar0dc-v0dc*cs)       !vper/c
         p_perp=vmax1dc*sn             !vpar/c
         if(p_perp.lt.1.d-12) p_perp=1.d-3
        
c         call power_CD_n_calc_theta(p_perp,p_par,y_loc,n_parallel,
c     &   n_perp,theta_temperature,n_harm_adj,
c     &   f_power_n_harm,f_CD_n_harm)

         p_int=1.d0
        
         if((j.ne.jmax).and.(j.ne.1)) then         
           sn=dsin(xint+dth)
           p_perp_pl=vmax1dc*sn             !vper/c
           sn=dsin(xint-dth)
           p_perp_min=vmax1dc*sn            !vper/c  
           h=0.5d0*(p_perp_pl-p_perp_min)
         else
           if(j.eq.1) then
             sn=dsin(xint+dth)
             p_perp_pl=vmax1dc*sn             !vper/c
             h=0.5d0*(p_perp_pl-p_perp)
           else
             !j=jmax
             sn=dsin(xint-dth)
             p_perp_min=vmax1dc*sn             !vper/c
             h=0.5d0*(p_perp-p_perp_min)
           endif
         endif
                 
         power_nharm_adj=power_nharm_adj+h*f_power_n_harm
         CD_nharm_adj=CD_nharm_adj+h*f_CD_n_harm
      enddo 
      goto 20
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     end
c--------------------------------------------------------------------

 2    continue
c-------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     begin
c------------------------------------------------------------------- 
      do j=1,jmax
        p0_perp=h*(j-0.5)

c        write(*,*)'j,p0_perp',j,p0_perp 

        call power_CD_nharm_under_integral_functions_test(n_radial0,
     &  p0_perp,n_parallel,n_perp,theta_temperature,b_ratio,
     &  e_x,e_y,e_z,y_loc,temp_kev, 
     &  unorm, sin_trap, 
     &  kmax,Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &  D_n_harm_res,d_chi_d_u0,d_chi_d_theta0,
     &  d_g_n_harm_d_u0_parallel,
     &  n_harm_adj,f_power_n_harm,f_CD_n_harm)

c-------to plot QL diffision flux
        kmax_ar(j)= kmax
        Gamma_a_n_harm_u0_ar(j)=Gamma_a_n_harm_u0
        Gamma_a_n_harm_theta_0_ar(j)=Gamma_a_n_harm_theta_0
        D_n_harm_res_ar(j)=D_n_harm_res
        p_perp_d_vt_ar(j)=p0_perp*clight/unorm   !p_perp/(m*unorm)
        p_perp_d_vt_ar_s(j)=p_perp_d_vt_ar(j)**2 ![p_perp/(m*unorm)]**2
        p_perp_d_vt_ar_4(j)=p_perp_d_vt_ar(j)**4 ![p_perp/(m*unorm)]**4
        power_spectrum(j)=f_power_n_harm
        cur_spectrum(j)=f_CD_n_harm
        d_chi_d_u0_ar(j)=d_chi_d_u0
        d_chi_d_theta0_ar(j)=d_chi_d_theta0
        d_g_n_harm_d_u0_parallel_ar(j)= d_g_n_harm_d_u0_parallel

c        write(*,*)'j,p0_perp/unorm',j,p_perp_d_vt_ar(j)

c        write(*,*)'f_power_n_harm,f_CD_n_harm',
c     &             f_power_n_harm,f_CD_n_harm

        power_nharm_adj=power_nharm_adj+f_power_n_harm
        CD_nharm_adj=CD_nharm_adj+f_CD_n_harm
 
c       write(*,*)'f_CD_n_harm/f_power_n_harm',
c     &             f_CD_n_harm/f_power_n_harm

        if (f_power_n_harm.ne.0.d0) then
c          write(*,*)'f_CD_n_harm/f_power_n_harm',
c     &                f_CD_n_harm/f_power_n_harm
        endif
               
      enddo 

      power_nharm_adj=power_nharm_adj*h
      CD_nharm_adj=CD_nharm_adj*h

c      write(*,*)'before plot1dt n_radial0',n_radial0
      call plot1dt(p_perp_d_vt_ar,D_n_harm_res_ar,0,0,jmax,1,'linlin',
     &0.d0,0.d0,'(p_perp/vt)','D_n_harm')

      call plot1dt(p_perp_d_vt_ar_s,D_n_harm_res_ar,0,0,jmax,1,'linlin',
     &0.d0,0.d0,'(p_perp/vt)**2','D_n_harm')

      call plot1dt(p_perp_d_vt_ar_s,Gamma_a_n_harm_theta_0_ar,0,0,jmax,1
     &,'linlin',
     &0.d0,0.d0,'(p_perp/vt)**2','Gamma_a_n_harm_theta_0') 

      call plot1dt(p_perp_d_vt_ar_s,Gamma_a_n_harm_u0_ar,0,0,jmax,1,
     &'linlin',
     &0.d0,0.d0,'(p_perp/vt)**2','Gamma_a_n_harm_u0')

c      call plot1dt(p_perp_d_vt_ar_4,Gamma_a_n_harm_theta_0_ar,0,0,jmax,1
c     &,'linlin',
c     &0.d0,0.d0,'(p_perp/vt)**4','Gamma_a_n_harm_theta_0') 

c      call plot1dt(p_perp_d_vt_ar_4,Gamma_a_n_harm_u0_ar,0,0,jmax,1,
c     &'linlin',
c     &0.d0,0.d0,'(p_perp/vt)**4','Gamma_a_n_harm_u0')





      call plot1dt(p_perp_d_vt_ar,power_spectrum,0,0,jmax,1,
     &'linlin',
     &0.d0,0.d0,'(p_perp/vt)','power_spectrum')

      call plot1dt(p_perp_d_vt_ar,cur_spectrum,0,0,jmax,1,
     &'linlin',
     &0.d0,0.d0,'(p_perp/vt)','cur_spectrum') 

c      call plot1dt(p_perp_d_vt_ar,kmax_ar,0,0,jmax,1,
c     &'linlin',
c     &0.d0,0.d0,'(p_perp/vt)','kmax')  

      call plot1dt(p_perp_d_vt_ar,d_chi_d_u0_ar,0,0,jmax,1,
     &'linlin',
     &0.d0,0.d0,'(p_perp/vt)','d_chi_d_u0_ar')


      call plot1dt(p_perp_d_vt_ar,d_chi_d_theta0_ar,0,0,jmax,1,
     &'linlin',
     &0.d0,0.d0,'(p_perp/vt)','d_chi_d_theta0_ar')

      call plot1dt(p_perp_d_vt_ar,d_g_n_harm_d_u0_parallel_ar,0,0,jmax,1
     &,'linlin',
     &0.d0,0.d0,'(p_perp/vt)','d_g_n_harm_d_u0_parallel_ar')

      goto 20
c ------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     end
c-------------------------------------------------------------------

 3    continue
cSm060306
c -------------------------------------------------------------------
c     integration along the resonance curve using trapezoidal
c     formula
c     begin
c--------------------------------------------------------------------
      do j=0,jmax
         p_int=1.d0
         if((j.eq.1).or.(j.eq.jmax)) p_int=0.5d0
         p0_perp=h*j
         if(p0_perp.lt.1.d-3) p0_perp=1.d-3
         if(j.eq.jmax) p0_perp=p0_perp-h*1.d-3

          call power_CD_nharm_under_integral_functions_test(n_radial0,
     &    p0_perp,n_parallel,n_perp,theta_temperature,b_ratio,
     &    e_x,e_y,e_z,y_loc,temp_kev,
     &    unorm, sin_trap,
     &    kmax,Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &    D_n_harm_res,d_chi_d_u0,d_chi_d_theta0,
     &    d_g_n_harm_d_u0_parallel,
     &    n_harm_adj,f_power_n_harm,f_CD_n_harm)

        power_nharm_adj=power_nharm_adj+p_int*f_power_n_harm
        CD_nharm_adj=CD_nharm_adj+p_int*f_CD_n_harm
    
      enddo
 
      power_nharm_adj=power_nharm_adj*h
      CD_nharm_adj=CD_nharm_adj*h
     
      goto20
c-------------------------------------------------------------------
c     end integration along the resonance curve using trapezoidal
c     formula
c     end
c-------------------------------------------------------------------

 4    continue
cSm070417
c -------------------------------------------------------------------
c     integration along the resonance curve using 
c     the adaptive Simpson function
c     begin
c--------------------------------------------------------------------
c      write(*,*)'before calc_integral_array_by_adaptive_simpson'

c      call  calc_power_CD_integral_array_by_adaptive_simpson
c     &(y_loc,n_parallel,n_perp,theta_temperature,n_harm_adj,
c     &p0_perp_max,epsi,power_nharm_adj,CD_nharm_adj) 

cSAP110316
c      call  calc_power_CD_integral_array_by_adaptive_simpson
c     &(n_radial0, 
c     & b_ratio,y_loc,sin_trap,temp_kev,   
c     &e_x,e_y,e_z,
c     &n_harm_adj,
c     &n_parallel,n_perp,theta_temperature,unorm,umax,
c     &n_relt_intgr_adj,
c     &p0_perp_max,epsi,power_nharm_adj,CD_nharm_adj) 
      call  calc_power_CD_integral_array_by_adaptive_simpson
     &(1,n_radial0, 
     & b_ratio,y_loc,sin_trap,temp_kev,   
     &e_x,e_y,e_z,
     &n_harm_adj,
     &n_parallel,n_perp,theta_temperature,unorm,umax,
     &n_relt_intgr_adj,
     &p0_perp_max,epsi,power_nharm_adj,CD_nharm_adj) 
c----------------------------------------------------------------------                  
 20   continue
 10   continue

      return
      end
      

      subroutine power_CD_nharm_under_integral_functions_test
     &(n_radial0,
     &p0_perp,nll,n_perp,theta_temperature,b_ratio,e_x,e_y,e_z,y_loc,
     &temp_kev,unorm,sin_trap,
     &kmax,Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &D_n_harm_res,d_chi_d_u0,d_chi_d_theta0,
     &d_g_n_harm_d_u0_parallel,
     &n_harm_adj,f_power_n_harm,f_CD_n_harm)
c----------------------------------------------------------------
c     Calculates under integral functions
c     f_power_n_harnm,f_CD_n_harm for ADJ calculations
c     at one radial point psis(n_radial0) from chi radial mesh 
c     input
c       n_radial0    is a number of the point at the radial chi mesh
c       p0_perp  - momentum divided by (mc)
c       nll     - N_parallel
c       n_perp  - N_perpendicular 
c       theta_temperature   - mc**2/T
c       b_ratio - B(l)/b_min
c       e_x,e_y,e_z  Stix frame complex electric filed polarization
c       y_loc     = |omega_ce|/omega for the electron rest mass. 
c                   It is positive   
c       temp_kev     electron temperature [kev]
c       unorm        normalization velocity sqrt(temp_kev/t) [sm/sec]
c       sin_trap     dsqrt(bmin/bmax)
c       n_harm_adj  - number of the given resonance harmonic
c
c     output
c      f_power_n_harm,f_CD_n_harm  under integral functions
c      f_CD_n_harm                 for power and CD ADSJ calculations
c-------------------------------------------------------
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input        
      double precision p0_perp,nll,n_perp,theta_temperature,b_ratio,
     &y_loc,temp_kev,unorm,sin_trap    

      integer n_harm_adj,n_radial0

      complex*16
     &e_x,e_y,e_z              !Stix frame electric filed polarization

c-----external root_res
      double precision terp2p_Sm

c-----output
      double precision f_power_n_harm,f_CD_n_harm

c-----local
      double precision gamma, p_par_rl(2),ps,pi,
     &dgdp_par,p_perp,p0_par,clight,u0,
     &Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &D_n_harm_res,
     &d_g_n_harm_d_u0_parallel,dfm_du0,
     &sin_theta0_res,cos_theta0_res,coef_trap,
     &theta0_res,theta0_odd,
     &d_chi_d_u0,d_chi_d_theta0,chi

c-----test chi slpine
      real*8, dimension(1:imax_chi,1:nmax_chi) :: chi_2d    
      real*8, dimension(1:imax_chi,1:nmax_chi) :: chi_tt_2d  
      real*8, dimension(1:imax_chi,1:nmax_chi) :: chi_uu_2d   
      real*8, dimension(1:imax_chi,1:nmax_chi) :: chi_ttuu_2d 
      real*8, dimension(1:(1+3*nmax_chi)) :: wk_chi_2d
      real*8, dimension(1:imax_chi) :: th0a_chi_1d
      integer ibd(4),i,n 
c-----end_test chi spline

      integer k,kmax,idm
c     kmax - total number of the resonace condition root p_perp_k(p_perp)
      
    

      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]
      p_perp=p0_perp*dsqrt(b_ratio) 

c      write(*,*)'power_CD_nharm_under_integral_functions p_perp',p_perp
      
c-------------------------------------------------------------------
c     calculations roots p_par_rl of the resonance condition
c     gamma=N_par*p_par_rl+n_harm*Y
c-------------------------------------------------------------------
      call root_res(p_perp,nll,n_harm_adj,y_loc,kmax,p_par_rl)

c      write(*,*)'p_perp,n_harm_adj,nll,y_loc,kmax,p_par_rl',
c     &p_perp, n_harm_adj,nll,y_loc,kmax,p_par_rl

c-----initialize 
      f_power_n_harm=0.d0
      f_CD_n_harm=0.d0  

cSAP090827
c      idm=imax_chi_a
      idm=imax_chi

c------test 2d pchi spline
        idm=imax_chi
        ibd(1)=2
        ibd(2)=2
        ibd(3)=4
        ibd(4)=4   
        do i=1,imax_chi
         do n=1,nmax_chi
           chi_2d(i,n)=chi_3d(i,n,n_radial0)
         enddo
         th0a_chi_1d(i)=th0a_chi(i,n_radial0)
        enddo

        call coeff2_Sm(imax_chi,th0a_chi_1d,nmax_chi,ua,
     &   chi_2d,
     &   chi_tt_2d,chi_uu_2d,
     &   chi_ttuu_2d,
     &   idm,ibd,wk_chi_2d,imax_chi)
       
c------end test 2d chi spline
      if (kmax.eq.0) goto 20 !no roots
   
      do k=1,kmax ! the sum over all P_parallel_k(perp) roots 
   
         ps=p_perp*p_perp+p_par_rl(k)*p_par_rl(k)

         write(*,*)'ps',ps

         p0_par=dsqrt((p_par_rl(k)**2*b_ratio+p_perp**2*(b_ratio-1.d0))
     &          /b_ratio)
ctest
         ps=p0_perp*p0_perp+p0_par*p0_par
         write(*,*)'p0s',ps
cendtest
        
         u0=dsqrt(ps)*clight/unorm  ! normalized electron momentum per mass p/(m*unorm)

         write(*,*)'k,u0,umax',k,u0,umax

         if(u0.gt.umax) goto 10

         call  ql_flux_nharm_test(n_harm_adj,
     &   nll,n_perp,u0,t,temp_kev,unorm,
     &   b_ratio,e_x,e_y,e_z,y_loc,
     &   sin_theta0_res,cos_theta0_res,dfm_du0,
     &   Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &   D_n_harm_res,
     &   d_g_n_harm_d_u0_parallel,
     &   gamma)

         if (sin_theta0_res.gt.sin_trap) then
           coef_trap=0.d0
         else     
           coef_trap=1.d0
         endif

c------end test 2d chi spline
         theta0_res=dacos(cos_theta0_res)
c         write(*,*)'sin_theta0_res,sin_trap,coef_trap',
c     &   sin_theta0_res,sin_trap,coef_trap
c--------initialization: for trapped particles if will give zero dertivatives
         d_chi_d_theta0=0.d0
         d_chi_d_u0=0.d0

         if((theta0_res.le.th0a_chi(imax_chi,n_radial0)).or.
     &     (theta0_res.ge.pi-th0a_chi(imax_chi,n_radial0))) then          
c         if (sin_theta0_res.lt.sin_trap) then
c-----------passing particles
c            write(*,*)'theta0_res,u0',theta0_res,u0
           if (theta0_res.le.th0a_chi(imax_chi,n_radial0)) then
cSAP090130
               if (i_chi_interpolation.eq.2) then
                  call  interpolation_chi(2,n_radial0,u0,theta0_res,
     &            chi,d_chi_d_u0,d_chi_d_theta0)
c                  write(*,*)'interpolation: d_chi_d_u0,d_chi_d_theta0',
c     &                    d_chi_d_u0,d_chi_d_theta0
               endif
cSAP090130
               if (i_chi_interpolation.eq.1) then
c-----------------spline interpolation
                  d_chi_d_u0=terp2p_Sm(theta0_res,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)

                  d_chi_d_theta0=terp2p_Sm(theta0_res,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)

              endif
           else
               theta0_odd=pi-theta0_res
c               write(*,*)'theta0_res,theta0_odd,u0',
c     &                     theta0_res,theta0_odd,u0

cSAP090130
               if (i_chi_interpolation.eq.2) then
                  call  interpolation_chi(2,n_radial0,u0,theta0_odd,
     &            chi,d_chi_d_u0,d_chi_d_theta0)
                  d_chi_d_u0=-d_chi_d_u0
               endif

c               write(*,*)'interpolation d_chi_d_u0,d_chi_d_theta0',
c     &                    d_chi_d_u0,d_chi_d_theta0

cSAP090130
               if (i_chi_interpolation.eq.1) then
c-----------------spline interpolation      

                  d_chi_d_u0=-terp2p_Sm(theta0_odd,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)

                  d_chi_d_theta0=terp2p_Sm(theta0_odd,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)

c                 write(*,*)'spline d_chi_d_u0,d_chi_d_theta0',
c     &                    d_chi_d_u0,d_chi_d_theta0

c-----------------test 2d chi spline
c                 d_chi_d_u0=-terp2p_Sm(theta0_odd,u0,   
c     &           imax_chi,th0a_chi_1d,
c     &           nmax_chi,ua,
c     &           chi_2d,chi_tt_2d,
c     &           chi_uu_2d,
c     &           chi_ttuu_2d,idm,0,1,1,imax_chi)
            
c                 d_chi_d_theta0=terp2p_Sm(theta0_odd,u0,
c     &           imax_chi,th0a_chi_1d,
c     &           nmax_chi,ua,
c     &           chi_2d,chi_tt_2d,
c     &           chi_uu_2d,
c     &           chi_ttuu_2d,idm,1,0,1,imax_chi)
c                 write(*,*)'spline d_chi_d_u0,d_chi_d_theta0',
c     &                    d_chi_d_u0,d_chi_d_theta0
              endif  ! if(i_chi_interpolation.eq.1)

           endif
cSAP081208
c           f_CD_n_harm=f_CD_n_harm+coef_trap*(u0/gamma)*cos_theta0_res*
c     &               dfm_du0*(Gamma_a_n_harm_u0*d_chi_d_u0+
c     &                Gamma_a_n_harm_theta_0*d_chi_d_theta0/u0)/
c     &               dabs(d_g_n_harm_d_u0_parallel)
         endif 
cSAP090727 due to Gamma=-Gamma_a*...
c         f_power_n_harm=f_power_n_harm+(u0/gamma)**2*
         f_power_n_harm=f_power_n_harm-(u0/gamma)**2*
cSAP090727 due to Lambdad=(u_0/gam_0)dabs(cos(theta_0))*tau_ B/L     
c     &              cos_theta0_res*
     &               dabs(cos_theta0_res)*
     &      Gamma_a_n_harm_u0*dfm_du0/dabs(d_g_n_harm_d_u0_parallel)
         
c         write(*,*)'u0,gamma,cos_theta0_res, Gamma_a_n_harm_u0',
c     &              u0,gamma,cos_theta0_res, Gamma_a_n_harm_u0

c         write(*,*)'dfm_du0,d_g_n_harm_d_u0_parallel,f_power_n_harm',
c     &              dfm_du0,d_g_n_harm_d_u0_parallel,f_power_n_harm
cSAP090727 due to Gamma=-Gamma_a*...
c         f_CD_n_harm=f_CD_n_harm+coef_trap*(u0/gamma)*
         f_CD_n_harm=f_CD_n_harm-coef_trap*(u0/gamma)*
cSAP090727 due to Lambdad=(u_0/gam_0)dabs(cos(theta_0))*tau_ B/L     
c     &              cos_theta0_res*
     &               dabs(cos_theta0_res)*
     &               dfm_du0*(Gamma_a_n_harm_u0*d_chi_d_u0+
     &               Gamma_a_n_harm_theta_0*d_chi_d_theta0/u0)/
     &               dabs(d_g_n_harm_d_u0_parallel)


c         write(*,*)'coef_trap,Gamma_a_n_harm_u0,d_chi_d_u0',
c     &              coef_trap,Gamma_a_n_harm_u0,d_chi_d_u0

c         write(*,*)'Gamma_a_n_harm_theta_0,d_chi_d_theta0,f_CD_n_harm',
c     &              Gamma_a_n_harm_theta_0,d_chi_d_theta0,f_CD_n_harm



 10      continue       
      enddo !kmax

 20   continue

      return
      end


c     subroutine power_CD_n_calc_theta(n_radial0,
c     &b_ratio,y_loc,temp_kev,e_x,e_y,e_z,unorm,umax,
c     &p0_perp,p0_par,y,nll,n_perp,theta_temperature,n_harm_adj,
c     &f_power_n_harm,f_CD_n_harm)
cc     calculates under integral functions: f_power_n_harm,f_CD_n_harm
cc     input
cc       n_radial0    is a number of radial chi mesh
cc       b_ratio      =B/B_min
cc       y_loc        =|omega_ce|/omega positive value
cc       sin_trap     =sqrt(bmin/bmax) at n_radial0 flux surface
cc       temp_kev,    electron temperatue [KeV]
cc       e_x,e_y,e_z  Stix frame complex electric filed polarization
cc       unorm    = dsqrt(temp_kev/t) normalization velocity [sm/sec]  
cc       umax       maximal momentum normalized by unorm
cc       p_perp  - perpendicular momentum divided by (mc)
cc       p_par     parallel momentum divided by (mc) along resonance curve
cc       y       = omega_ce/omega for the electron rest mass    
cc       nll     - N_parallel
cc       np      - N_perpendicular 
cc       theta_temperature   - mc**2/T
cc       n_harm_adj   number of the given resonance harmonic
cc     output
cc       f_power_n_harm,f_CD_n_harm  under integral functionns
cc-------------------------------------------------------
      !implicit none
cc-----input        
c      real*8 p0_perp,y,nll,n_perp,theta_temperature,p0_par,b_ratio,
c     &y_loc,sin_trap,temp_kev,unorm,umax

c      complex*16 e_x,e_y,e_z
  
c      integer n_harm_adj,n_radial0

cc-----external root_res,     

cc-----output
c      real*8 f_power_n_harm,f_CD_n_harm

cc-----local
c      real*8 gamma, p_par_rl(2),eps,ps,ps_max,pi,
c     +dgdp_par,p_perpc    
        
c      integer k,kmax,i,j
cc     kmax - total number of the resonace condition root p_par_k(p_perp)
      
c      pi=4.d0*datan(1.d0)
cc-----calculations ot the roots p_par_rl of the resonance condition
cc-----gamma=N_par*p_par_rl+nY

c      call root_res(p_perp,nll,n_harm_adj,y_loc,kmax,p_par_rl)

cc-----initialize 
c      f_power_n_harm=0.d0
c      f_CD_n_harm=0.d0

c      if (kmax.eq.0) goto 20
  
c      ps=p_perp*p_perp+p_par_rl(k)*p_par_rk(k)
c      gamma=dsqrt(1.d0+ps)

c      call power_CD_nharm_under_integral_functions(n_radial0,
c     &p0_perp,nll,n_perp,theta_temperature,b_ratio,e_x,e_y,e_z,
c     &y_loc,temp_kev,unorm,sin_trap,
c     &n_harm_adj,f_power_n_harm,f_CD_n_harm)
  
cc-----resonance condition uses the delta function with argument
cc     g_delta=gamma-nll*p_par-n*y
cc     the derivative from this argument d(g_delta)/dp_par=dgdp_par

c      dgdp_par=dabs((p_par-nll*gamma)/gamma)     

c 10   continue
       
c 20   continue

c      return
c      end

      subroutine ql_flux_nharm_test(n_harm,
     &n_parallel,n_perp,u0,t,temp_kev,unorm,
     &b_ratio,e_x,e_y,e_z,y_loc,
     &sin_theta0_res,cos_theta0_res,dfm_du0,
     &Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &D_n_harm_res,
     &d_g_n_harm_d_u0_parallel,
     &gamma)
c--------------------------------------------------------------------------------
c     calculate quasi linear flux ,
c     derivative d_g_n_harm/d_u0_parallel, gamma
c     for power and CD calculations
c     for one given hyro-harmonic:  n_harm
c-------------------------------------------------------------------------------
      implicit none
c-----input
      integer
     &n_harm               !the number of hyro-harmonic

      real*8
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicula refractive index
     &u0,                  !normalized momentum p/(m*unorm) , 
     &t,                   !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1.
     &temp_kev,            !electron temperature [kev]
     &unorm,               !normalization velocity sqrt(temp_kev/(m*t)) [sm/sec]
     &b_ratio,             !B(l)/B(l=0)=B(l)/B_min 
     &y_loc                !|omega_ce|/omega it is positive
     
      complex*16
     &e_x,e_y,e_z              !Stix frame electric filed polarization

c-----output
      real*8 sin_theta0_res,cos_theta0_res,dfm_du0,gamma,
     &Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &d_g_n_harm_d_u0_parallel
   
c-----locals
      real*8 cos_theta_res,sin_theta_res,
     &clight,u0_dev_c,u0_parallel_dev_c,
     &u0_perp_dev_c,
     &ksi_res,
     &j_bessel_n_harm,j_bessel_n_harm_minus,j_bessel_n_harm_plus,
     &j_bessel(1:3),d_n_harm,
     &D_n_harm_res,
     &fm,mass_e,k_1_kev,theta_relativistic,pi,
     &sigma 

      complex*16 theta_n_harm_big_res,i_image

      integer nz,n_harm_abs
c-----external
c     DBESJ

      i_image=dcmplx(0.d0,1.d0)
      clight=2.99792458d10              !light speed [cm/sec]
      pi=4.d0*datan(1.d0)

      u0_dev_c=(u0*unorm)/clight
  
      gamma=dsqrt(1.d0+u0_dev_c**2)
c      gamma_min=dsqrt(1.d0+(umin*unorm/clight)**2)
c      write(*,*)'ql_flux t,temp_kev,unorm',t,temp_kev,unorm
c-------------------------------------------------------------------------------
c     pitch angle at the resonance curve at the local poloidal point
c-------------------------------------------------------------------------------
      cos_theta_res=(gamma-n_harm*y_loc)/(n_parallel*u0_dev_c) !
      sin_theta_res=dsqrt(1.d0- cos_theta_res**2)
c      write(*,*)'ql_flux u0,n_parallel,cos_theta_res,sin_theta_res',
c     & u0,n_parallel,cos_theta_res,sin_theta_res
c-------------------------------------------------------------------------------
c     pitch angle at the resonance curve at the midplane b=bmin
c-------------------------------------------------------------------------------
c      write(*,*)'ql_flux b_ratio',b_ratio

      sin_theta0_res=sin_theta_res/dsqrt(b_ratio)
      if (cos_theta_res.ge.0d0) then
        cos_theta0_res=dsqrt(1.d0-sin_theta0_res**2)
      else
        cos_theta0_res=-dsqrt(1.d0-sin_theta0_res**2)
      endif

      u0_parallel_dev_c=u0_dev_c*cos_theta0_res
      u0_perp_dev_c=u0_dev_c*sin_theta0_res
c      write(*,*)'ql_flux ucos_theta0_res,sin_theta0_res',
c     & cos_theta0_res,sin_theta0_res
c--------------------------------------------------------------------------------
c     Bessel function agrgument
c-------------------------------------------------------------------------------
      ksi_res=n_perp*u0_dev_c*sin_theta_res/y_loc
      write(*,*)'ksi_res',ksi_res
c------------------------------------------------------------------------------
c      nz controles bessel function calculations
c      nz should be =0. If nz=1 sum bessel function(dBESJ) will be J=0

      d_n_harm=n_harm-1.d0
      if (n_harm.ge.1) then
         call DBESJ( ksi_res,d_n_harm,3,j_bessel,nz)
         j_bessel_n_harm_minus=j_bessel(1)
         j_bessel_n_harm=j_bessel(2)
         j_bessel_n_harm_plus=j_bessel(3)
      else
         if(n_harm.eq.0)then
           call DBESJ( ksi_res,1.d0,1,j_bessel_n_harm_plus,nz)
c          write(*,*)'ql_flux after dbesj 1 nz,j_bessel_n_harm_plus',
c          nz,j_bessel_n_harm_plus',
           call DBESJ( ksi_res,0.d0,1,j_bessel_n_harm,nz)
           j_bessel_n_harm_minus=-j_bessel_n_harm_plus
         else
c----------n_harm.lt.0
           n_harm_abs=-n_harm
           call DBESJ( ksi_res,dfloat(n_harm_abs)-1.d0,3,j_bessel,nz)  
           j_bessel_n_harm_minus=(-1)**(n_harm_abs+1)*j_bessel(3)
           j_bessel_n_harm=(-1)**n_harm_abs*j_bessel(2)
           j_bessel_n_harm_plus=(-1)**(n_harm_abs-1)*j_bessel(1)
         endif
      endif
cSAP080411      
c      theta_n_harm_big_res=0.5d0*(e_x-i_image*e_y)*j_bessel_n_harm_plus
c     &                   +0.5d0*(e_x+i_image*e_y)*j_bessel_n_harm_minus
      theta_n_harm_big_res=0.5d0*(e_x-i_image*e_y)*j_bessel_n_harm_minus
     &                   +0.5d0*(e_x+i_image*e_y)*j_bessel_n_harm_plus
     &                 +e_z*j_bessel_n_harm*cos_theta_res/sin_theta_res

      write(*,*)'ql_flux e_x,e_y,e_z',e_x,e_y,e_z 
      write(*,*)'j_bessel_n_harm_plus,j_bessel_n_harm_minus',
     &           j_bessel_n_harm_plus,j_bessel_n_harm_minus
      write(*,*)'j_bessel_n_harm,cos_theta_res,sin_theta_res',
     &           j_bessel_n_harm,cos_theta_res,sin_theta_res
      write(*,*)'theta_n_harm_big_res',theta_n_harm_big_res

      D_n_harm_res=0.5d0*pi*
     &         theta_n_harm_big_res*dconjg(theta_n_harm_big_res)
      write(*,*)' D_n_harm_res', D_n_harm_res
  
      Gamma_a_n_harm_u0= gamma/(u0*dabs(cos_theta_res))*D_n_harm_res*
     &     ((n_parallel*u0_dev_c*sin_theta_res*cos_theta_res+
     &       n_harm*y_loc*sin_theta_res)/gamma)**2

      write(*,*)'Gamma_a_n_harm_u0=',Gamma_a_n_harm_u0

      Gamma_a_n_harm_theta_0= gamma/(u0*dabs(cos_theta_res))*
     & D_n_harm_res*
     & (sin_theta0_res/cos_theta0_res)/(sin_theta_res/cos_theta_res)*
     & (-((n_parallel*u0_dev_c*sin_theta_res/gamma)**2-
     &    (n_harm*y_loc/gamma)**2)*sin_theta_res*cos_theta_res+
     &   (n_parallel*u0_dev_c*sin_theta_res/gamma)*(n_harm*y_loc/gamma)
     &    *(cos_theta_res**2-sin_theta_res**2))
 
      write(*,*)'Gamma_a_n_harm_theta_0=', Gamma_a_n_harm_theta_0

      sigma=1.d0
      if (u0_parallel_dev_c.lt.0.d0) sigma=-1.d0
   
      d_g_n_harm_d_u0_parallel=(
     &( sigma*n_parallel*
     &  dsqrt(u0_parallel_dev_c**2-(b_ratio-1.d0)*u0_perp_dev_c**2)/
     &  gamma**2+n_harm*y_loc/gamma**2)*(u0_parallel_dev_c/gamma)-
     &sigma*n_parallel*u0_parallel_dev_c/(gamma*
     & dsqrt(u0_parallel_dev_c**2-u0_perp_dev_c**2*(b_ratio-1.d0)))
     & )/(clight/unorm)

      write(*,*)'d_g_n_harm_d_u0_parallel',d_g_n_harm_d_u0_parallel  

      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)

      theta_relativistic=mass_e*clight**2/(k_1_kev*temp_kev)
  
cSAP070831
c     fm=dexp(-theta_relativistic*gamma)
      fm=dexp(-theta_relativistic*(gamma-1.d0))
      write(*,*)'temp_kev,theta_relativistic,fm',
     &           temp_kev,theta_relativistic,fm
     
c      fm=dexp(-theta_relativistic*(gamma-gamma_min)) 
      dfm_du0=-fm*theta_relativistic*u0_dev_c/gamma
c      write(*,*)'(gamma-gamma_min),fm ',(gamma-gamma_min),fm 
      return
      end


c     Modificated method for vector function.
c     
c     Original description and functions are from:
c
C PAGE 187-190: NUMERICAL MATHEMATICS AND COMPUTING, CHENEY/KINCAID, 1985
C
C FILE: SIMP.FOR
C
C ADAPTIVE SCHEME FOR SIMPSON'S RULE (SIMP,ASMP,PUSH,POP,FCN)
c============================================================================

      subroutine fcn_adj_vector_test(p0_perp,fcn)
c---------------------------------------------------------------------------
c     culculates under integral function for 
c     power and CD calulations using 
c---------------------------------------------------------------------------

      implicit none
      integer k_max !max number of elements
      parameter (k_max=2) 
c----------------------------------------------------------------
c     input
c----------------------------------------------------------------- 
      real*8 p0_perp !under integral function argument
c-----------------------------------------------------------------
c     output
c-----------------------------------------------------------------
      real*8 fcn(k_max) !array of under integral fuctions values
c-----------------------------------------------------------------
      integer n_harm_adj_l,n_radial0_l

      real*8 y_l,n_parallel_l,n_perp_l,theta_temperature_l,b_ratio_l,
     &temp_kev_l,unorm_l,sin_trap_l

      complex*16 e_x_l,e_y_l,e_z_l

      common /fcn_adj_input/
     &y_l,n_parallel_l,n_perp_l,theta_temperature_l,
     &temp_kev_l,unorm_l,sin_trap_l,b_ratio_l,
     &e_x_l,e_y_l,e_z_l,
     &n_harm_adj_l,n_radial0_l

     
c-----------------------------------------------------------------
c     locals
c----------------------------------------------------------------
      real*8 
     &f_power_n_harm,f_CD_n_harm, ! under integral functions
     &Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &D_n_harm_res,d_chi_d_u0,d_chi_d_theta0,
     &d_g_n_harm_d_u0_parallel

      integer kmax

      call power_CD_nharm_under_integral_functions_test(n_radial0_l,
     &p0_perp,n_parallel_l,n_perp_l,theta_temperature_l,b_ratio_l,
     &e_x_l,e_y_l,e_z_l,y_l,
     &temp_kev_l,unorm_l,sin_trap_l,
     &kmax,Gamma_a_n_harm_u0,Gamma_a_n_harm_theta_0,
     &D_n_harm_res,d_chi_d_u0,d_chi_d_theta0,
     &d_g_n_harm_d_u0_parallel,
     &n_harm_adj_l,f_power_n_harm,f_CD_n_harm)

      FCN(1) =f_power_n_harm 
      FCN(2) =f_CD_n_harm
     
      return
      end


 
   
      subroutine  calc_power_CD_integral_array_by_adaptive_simpson_test
     &(n_radial0, 
     & b_ratio,y_loc,sin_trap,temp_kev,   
     &e_x,e_y,e_z,
     &n_harm_adj,
     &n_parallel,n_perp,theta_temperature,unorm,umax,
     &n_relt_intgr_adj,
     &p0_perp_max,epsi,power_nharm_adj,CD_nharm_adj)      
c------------------------------------------------------------------------
c     calculate intergals power and CD (power_nharm_adj,CD_nharm_adj)    
c     using ADJ and Sypson adaptive method
c     for the EC harmonic with number 'n_harm_adj'
c-------------------------------------------------------------------------
       implicit none
c------------------------------------------------------------------------
c      input
c--------------------------------------------------------------------------
c       n_radial0    is a number of radial chi mesh
c       b_ratio      =B/B_min
c       y_loc        =|omega_ce|/omega for the electron rest mass 
c                     It has positive value
c       sin_trap     =sqrt(bmin/bmax) at n_radial0 flux surface
c       temp_kev,    electron temperature [KeV] 
c       e_x,e_y,e_z  Stix frame complex electric filed polarization
c       n_harm_adj        - EC harmonic number
c       n_parallel          N_parallel
c       n_perp      N_perpendicular      
c       theta_temperature    = mc**2/T
c       unorm    = dsqrt(temp_kev/t) normalization velocity [sm/sec]  
c       umax       maximal momentum normalized by unorm
c       n_relt_intgr_adj the number of points for the integration over p0_perp
c       p0_perp_max is  the maximal boundary of integration
c       epsi is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------------------
      integer n_radial0,n_harm_adj, n_relt_intgr_adj

      real*8 b_ratio,y_loc,sin_trap,temp_kev,n_parallel,n_perp,    
     &theta_temperature,unorm,umax,epsi,p0_perp_max
      
      complex*16 e_x,e_y,e_z

c------------------------------------------------------------------------
c     output
c------------------------------------------------------------------------
      real*8 power_nharm_adj,CD_nharm_adj 
c-------------------------------------------------------------------------
c     locals
c--------------------------------------------------------------------------
      integer k_max !max number of elements
      parameter (k_max=2)
      integer ist,ifs,k   
c      PARAMETER (IST=5,IFS=16)
c      PARAMETER (IST=9,IFS=256)  
c       PARAMETER (IST=13,IFS=10048)
c       PARAMETER (IST=16,IFS=100048)
c       PARAMETER (IST=17,IFS=200048)!070723-old
      PARAMETER (IST=20,IFS=200048)

      REAL*8 STACK(IST,3),FSTACK(IFS,3)
      real*8 a,b,sum(k_max)
      integer lvmax,iflag,i
      real*8 fcn(k_max)
    
      integer n_harm_adj_l,n_radial0_l

      real*8 y_l,n_parallel_l,n_perp_l,theta_temperature_l,b_ratio_l,
     &temp_kev_l,unorm_l,sin_trap_l

      complex*16 e_x_l,e_y_l,e_z_l
     
      common /fcn_adj_input/
     &y_l,n_parallel_l,n_perp_l,theta_temperature_l,
     &temp_kev_l,unorm_l,sin_trap_l,b_ratio_l,
     &e_x_l,e_y_l,e_z_l,
     &n_harm_adj_l,n_radial0_l
 
      EXTERNAL FCN_adj_vector
c      DATA EPSI/5.0d-5/, LVMAX/4/
c      DATA EPSI/5.0d-5/, LVMAX/7/
c      DATA EPSI/5.0d-5/, LVMAX/15/
c      DATA EPSI/5.0d-5/, LVMAX/12/
c      DATA EPSI/5.0d-5/, LVMAX/16/
c       DATA EPSI/5.0d-6/, LVMAX/16/  !070722-new
c       DATA EPSI/1.0d-6/, LVMAX/16/  !070722-new

c       DATA EPSI/1.0d-3/, LVMAX/16/ !070722-old
c      DATA EPSI/1.0d-3/, LVMAX/19/
      DATA LVMAX/16/  !070722-new



c-----------------------------------------------------------------
c     set common /fcn_adj_input/
c-----------------------------------------------------------------
      y_l=y_loc
      n_parallel_l=n_parallel
      n_perp_l=n_perp
      theta_temperature_l=theta_temperature
      n_harm_adj_l=n_harm_adj
      n_radial0_l=n_radial0     
      temp_kev_l=temp_kev
      unorm_l=unorm
      sin_trap_l=sin_trap
      b_ratio_l=b_ratio
      e_x_l=e_x
      e_y_l=e_y
      e_z_l=e_z
      n_harm_adj_l=n_harm_adj
      n_radial0_l=n_radial0

c-----------------------------------------------------------------
c     integration boundaries
c------------------------------------------------------------------
c      A = 0.0d0
      A = 1.d-5
c      A = 1.d-3

      B = p0_perp_max-1.d-5
     
      CALL ASMP(FCN_adj_vector,k_max,fcn,
     &A,B,EPSI,LVMAX,STACK,IST,FSTACK,IFS,
     &SUM,IFLAG) 



      IF(IFLAG .EQ. 0) THEN 
c$$$        PRINT 4   
      ELSE
c        PRINT 5   
        write(*,*) '          WITH BAD SUBINTERVALS:'  
        DO 2 I=1,IFLAG      
c          PRINT 6,FSTACK(I,1),FSTACK(I,2),INT(FSTACK(I,3))
          write(*,*) FSTACK(I,1),FSTACK(I,2),INT(FSTACK(I,3))
    2   CONTINUE  
      END IF      
   3  FORMAT(//5X,'APPROXIMATE INTEGRAL =',E22.14/)   
   4  FORMAT(10X,'WITH NO BAD SUBINTERVALS')    
   5  FORMAT(10X,'WITH BAD SUBINTERVALS:')      
   6  FORMAT(10X,'[',F10.5,',',F10.5,']',2X,'LEVEL =',I5) 

c-------------------------------------------------------

      power_nharm_adj=Sum(1) 
      CD_nharm_adj=Sum(2)

     
      return
      END 
     
      subroutine CD_adj_efficiency_test(n_parallel,n_perp,
     &psi_loc,theta_pol,e_x,e_y,e_z,
     &cd_efficiency)
c-----------------------------------------------------------------
c     calculate CD efficiency using adj chi function
c-----------------------------------------------------------------       
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input
      real*8      
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicular refractive index
     &psi_loc,             !polidal flux 
     &theta_pol,           !poloidal angle
     &length               !field line length normalized to  
c     &t,                  !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1. 
      complex*16 
     &e_x,
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1
c-----output
      real*8 
     &cd_efficiency        
c-----local
      integer  n_radial_minus,n_radial_plus
      real*8 weight_radial,cd_efficiency_minus,
     &cd_efficiency_plus,unorm_minus,unorm_plus
     
 
c      write(*,*)' CD_adj_efficiency e_x,e_y,e_z', e_x,e_y,e_z
c      write(*,*)' CD_adj_efficiency before find_psis_points'
c      write(*,*)' n_radial_minus,n_radial_plus,weight_radial',
c     &             n_radial_minus,n_radial_plus,weight_radial
  
      call  find_psis_points(psi_loc,
     & n_radial_minus,n_radial_plus,weight_radial)

      write(*,*)' CD_adj_efficiency_test after find_psis_points'
      write(*,*)'n_radial_minus,n_radial_plus',
     &n_radial_minus,n_radial_plus

      if((n_radial_minus.lt.npsi0).and.(n_radial_plus.gt.1)) then

         
         call CD_adj_efficiency_1_test(n_radial_minus,
     &   theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &   unorm_minus, 
     &   cd_efficiency_minus)
         write(*,*)'eff cd_efficiency_minus',cd_efficiency_minus

          call CD_adj_efficiency_1_test(n_radial_plus,
     &   theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &   unorm_plus, 
     &   cd_efficiency_plus)
         write(*,*)'eff cd_efficiency_plus',cd_efficiency_plus


         cd_efficiency=cd_efficiency_plus*weight_radial+
     &   cd_efficiency_minus*(1.d0-weight_radial)

      else
        if (n_radial_plus.eq.1) then  
         
           call CD_adj_efficiency_1_test(1,
     &     theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &     unorm_plus, 
     &     cd_efficiency)

            write(*,*)'eff cd_efficiency 1',cd_efficiency
        else
           if (n_radial_minus.eq.npsi0) then
        
              call CD_adj_efficiency_1_test(npsi0,
     &        theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &        unorm_minus, 
     &        cd_efficiency)

              write(*,*)'eff cd_efficiency_ npsi0',cd_efficiency

           endif
        endif
      endif

      write(*,*)'in CD_adj_efficiency cd_efficiency=',
     &cd_efficiency

      return
      end


c===lh_test

      subroutine lh_ql_flux_test(n_parallel,n_perp,u0,umin,t,temp_kev,
     &b_ratio,e_y,e_z,
cSAP080411
     &                   y_loc,
     &sin_theta0_res,cos_theta0_res,dfm_du0,
     &Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0,gamma,unorm)
c--------------------------------------------------------------------------------
c     calculate LH quasi linear flux and 
c     under integral functions for power and CD calculations
c-------------------------------------------------------------------------------
      implicit none
c-----input
      real*8
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicula refractive index
     &u0,                  !normalized momentum u/unorm 
c     &umax,                !maximal momentum normalized by unorm
c                           !unorm=dsqrt(temp_kev/t)
     &umin,                !minimal momentum at the resonace curve
                           !normalized by unorm
     &temp_kev,            !electron temperature [kev]
     &t,                   !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1.
     &b_ratio,             !B(l)/B(l=0)=B(l)/B_min 
cSAP080411
     &y_loc                !|omega_c_e/omega|
      complex*16
     &e_y,e_z              !Stix frame electric filed polarization

c-----output
      real*8 sin_theta0_res,cos_theta0_res,dfm_du0,gamma,
     &Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0, 
     &unorm  !sqrt(T/m) cm/sec
   
c-----locals
      real*8 cos_theta_res,sin_theta_res,
     &clight,u0_dev_c,
     &ksi_res,j_bessel_1,j_bessel_0,
     &D_n0_res,
     &fm,mass_e,k_1_kev,theta_relativistic,pi,
     &gamma_min

      complex*16 theta_0_big_res,i_image

      integer nz
c-----xternal
c     DBESJ

      i_image=dcmplx(0.d0,1.d0)
      clight=2.99792458d10              !light speed [cm/sec]
      pi=4.d0*datan(1.d0)

      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
      u0_dev_c=(u0*unorm)/clight
      gamma=dsqrt(1.d0+u0_dev_c**2)
      gamma_min=dsqrt(1.d0+(umin*unorm/clight)**2)
c      write(*,*)'lh_ql_flux t,temp_kev,unorm',t,temp_kev,unorm
c-------------------------------------------------------------------------------
c     pitch angle at the resonance curve at the local poloidal point
c-------------------------------------------------------------------------------
      cos_theta_res=gamma/(n_parallel*u0_dev_c)
      sin_theta_res=dsqrt(1.d0- cos_theta_res**2)
c      write(*,*)'lh_ql_flux u0,n_parallel,cos_theta_res,sin_theta_res',
c     & u0,n_parallel,cos_theta_res,sin_theta_res
c-------------------------------------------------------------------------------
c     pitch angle at the resonance curve at the midplane b=bmin
c-------------------------------------------------------------------------------
c      write(*,*)'lh_ql_flux b_ratio',b_ratio

      sin_theta0_res=sin_theta_res/dsqrt(b_ratio)
      if (cos_theta_res.ge.0d0) then
        cos_theta0_res=dsqrt(1.d0-sin_theta0_res**2)
      else
        cos_theta0_res=-dsqrt(1.d0-sin_theta0_res**2)
      endif
c      write(*,*)'lh_ql_flux ucos_theta0_res,sin_theta0_res',
c     & cos_theta0_res,sin_theta0_res
c--------------------------------------------------------------------------------
c     Bessel function agrgument
c-------------------------------------------------------------------------------
cSAP080411      ksi_res=n_perp*u0_dev_c*sin_theta_res/gamma
      ksi_res=n_perp*u0_dev_c*sin_theta_res/y_loc
c------------------------------------------------------------------------------
c      nz controles bessel function calculations
c      nz should be =0. If nz=1 the bessel function(dBESJ) will be J=0
      call DBESJ( ksi_res,0.d0,1,j_bessel_0,nz)
c      write(*,*)'lh_ql_flux after dbesj 0 nz,j_bessel_0',nz,j_bessel_0
      call DBESJ( ksi_res,1.d0,1,j_bessel_1,nz)
c      write(*,*)'lh_ql_flux after dbesj 1 nz,j_bessel',nz,j_bessel_1

cSAP080411      
c      theta_0_big_res=-i_image*e_y*j_bessel_1+
      theta_0_big_res=i_image*e_y*j_bessel_1+
     &                j_bessel_0*cos_theta_res*e_z/sin_theta_res

c      write(*,*)'lh_ql_flux e_y,e_z',e_y,e_z 
c      write(*,*)'theta_0_big_res',theta_0_big_res

      D_n0_res=0.5d0*pi*theta_0_big_res*dconjg(theta_0_big_res)
c      write(*,*)' D_n0_res', D_n0_res
  
cSAP080415
c      Gamma_a_LH_u0= (gamma/u0)*D_n0_res*dabs(cos_theta_res)*
      Gamma_a_LH_u0=- (gamma/u0)*D_n0_res*dabs(cos_theta_res)*
     &                (n_parallel*u0_dev_c*sin_theta_res/gamma)**2

c      write(*,*)'Gamma_a_LH_u0',Gamma_a_LH_u0

      Gamma_a_LH_theta_0=-Gamma_a_LH_u0*
     &                     (sin_theta0_res/cos_theta0_res)

      dg_d_theta0=dabs((n_parallel*u0_dev_c/gamma)*
     &                  b_ratio*sin_theta0_res*cos_theta0_res/
     &                  cos_theta_res) 

      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)

      theta_relativistic=mass_e*clight**2/(k_1_kev*temp_kev)
  
cSAP070831
c      fm=dexp(-theta_relativistic*gamma)
      fm=dexp(-theta_relativistic*(gamma-gamma_min)) 
      dfm_du0=-fm*theta_relativistic*u0_dev_c/gamma
c      write(*,*)'(gamma-gamma_min),fm ',(gamma-gamma_min),fm 
      return
      end

      subroutine CD_adj_LH_efficiency_1_test(n_parallel,n_perp,
     &n_radial0,theta_pol,e_y,e_z,
     &unorm,
     &lh_cd_efficiency)
c------------------------------------------------------
c     calculate LH CD efficiency using adj chi function
c     at one radial point psis(n_radial0) from chi radial mesh 
c------------------------------------------------------
      implicit none
      include 'param.i'
      include 'adj.i'
c-----input
      integer
     & n_radial0           ! is a number of radial chi mesh

      real*8      
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicular refractive index
     &theta_pol,           !poloidal angle [radian]
     &length               !field line length normalized to  
c     &t,                   !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1. 
      complex*16 
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1
c-----output
      real*8
     &pow_dens,            !Power density normalized to m*u_norm**2*n/tau_n 
     &lh_cd_dens,          !LH current drive density normalized to q*u_norm*n
                           !Here n is the density  
     &lh_cd_efficiency,
     &unorm                !sqrt(T/m) cm/sec
c-----external
      real*8 temperho,densrho,
     &bmax_psi,bmin_psi,psi_rho,rhopsi


c-----local
      integer i,n,n_radial,idm,nmin
      real*8 du,u0,     
     &rho_small,           !small radius
     &temp_kev,            !electron temperatura [KeV]
     &dense,               !electron dencity 10**13/cm**3     
     &b_ratio,             !B/B_min
cSAP080411
     &z,r,phi,             !bellow phi is set arbitratry phi=0
     &y_loc,               !electron y at point (z,r,phi)
     &sin_theta0_res,cos_theta0_res,dfm_du0,theta0_res,theta0_odd,
     &Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0,
     &gamma,pi,bmin,bmax,deltb,psi,th0max,sin_trap,
     &d_chi_du0,d_chi_dtheta0,u0max,chi,
     &weight_radial,
     &umin,               !minimal u0 at the resonance LH hyperbola
     &clight,              !light speed [cm/sec]
     &arg1,cln,
cfor test
     &del_power,del_cd
c-----external
      real*8 terp2p_Sm,b_d_bmin,y
       
      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]

c      write(*,*)'in CD_adj_LH_efficiency_1,e_y,e_z,',e_y,e_z

c      do n_radial=1,npsi0
c         write(*,*)'psis',psis(n_radial)
c         rho_small=rhopsi(psis(n_radial))
c         write(*,*)'rho_small',rho_small
c      enddo

      rho_small=rhopsi(psis(n_radial0))
      psi=psis(n_radial0)

      write(*,*)'n_radial0=',n_radial0

      write(*,*)'psis(n_radial0),rho_small',psis(n_radial0),rho_small

      temp_kev=temperho(rho_small,1)
      dense=densrho(rho_small,1) 
      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
      bmax= bmax_psi(psi)
      bmin= bmin_psi(psi)
	     
      deltb = (bmax-bmin)/bmin
c      write(*,*)'bmax,mbin,deltb',bmax,bmin,deltb   
      th0max = atan2(1.0d0,sqrt(deltb))
c      write(*,*)'th0max',th0max
      sin_trap=dsqrt(bmin/bmax)

c      write(*,*)'sin_trap=',sin_trap,'dsin(th0max)=',dsin(th0max)

      du = umax/nmax_chi

      pow_dens=0.d0
      lh_cd_dens=0.d0

      umin=clight/dsqrt(n_parallel**2-1.d0)/unorm ! normalized umin      
c      write(*,*)'umin',umin,'umin/dsqrt(2.d0)',umin/dsqrt(2.d0)
      nmin=(umin/du+0.5d0)
      nmin=(dabs(umin)/du-0.5d0)
      nmin=nmin+1


      b_ratio=b_d_bmin(n_radial0,theta_pol)
c      write(*,*)'n_radial0,theta_pol,b_ratio',
c     &n_radial0,theta_pol,b_ratio
cSAP080411
      phi=0.d0
      call zr_psith(psi,theta_pol,z,r)
      y_loc=y(z,r,phi,1)

      write(*,*)'in CD_adj_LH_efficiency_1  nmin,nmax_chi-1',
     &nmin,nmax_chi-1

      do n=nmin,nmax_chi-1
         u0=(n + 0.5E0)*du
        
c        write(*,*)'in CD_adj_LH_efficiency_1,e_y,e_z,',e_y,e_z

        call lh_ql_flux_test(n_parallel,n_perp,u0,umin,t,temp_kev,
     &                   b_ratio,e_y,e_z,
cSAP080411
     &                   y_loc,
     &                   sin_theta0_res,cos_theta0_res, dfm_du0,
     &                   Gamma_a_LH_u0,Gamma_a_LH_theta_0,dg_d_theta0,
     &                   gamma,unorm)

         write(*,*)'n,Gamma_a_LH_u0,sin_theta0_res,cos_theta0_res',
     &              n,Gamma_a_LH_u0,sin_theta0_res,cos_theta0_res
         write(*,*)'dfm_du0,gamma,dg_d_theta0,dfm_du0,gamma,dg_d_theta0
     &            ',dfm_du0,gamma,dg_d_theta0,dfm_du0,gamma,dg_d_theta0

c         pow_dens=pow_dens+du*(u0**4/gamma**2)*
c     &            sin_theta0_res*cos_theta0_res*Gamma_a_LH_u0*dfm_du0/
c     &            dg_d_theta0
  
c        del_power=du*(u0**4/gamma**2)*
c     &            sin_theta0_res*cos_theta0_res*Gamma_a_LH_u0*dfm_du0/
c     &            dg_d_theta0

         pow_dens=pow_dens+du*u0**2*sin_theta0_res*
     &            (u0/gamma)*
cSAP080415 in lambda it should be used |cos_theta0_res|
c     &            (u0/gamma*cos_theta0_res)*Gamma_a_LH_u0*dfm_du0/
     &            (u0/gamma*dabs(cos_theta0_res))*Gamma_a_LH_u0*dfm_du0/
     &            dg_d_theta0
 
        del_power=du*u0**2*sin_theta0_res*
     &            (u0/gamma)*
cSAP080415 in lambda it should be used |cos_theta0_res|
c     &            (u0/gamma*cos_theta0_res)*Gamma_a_LH_u0*dfm_du0/
     &            (u0/gamma*dabs(cos_theta0_res))*Gamma_a_LH_u0*dfm_du0/
     &            dg_d_theta0

        write(*,*)'del_power',del_power

         write(*,*)'cos_theta0_res,Gamma_a_LH_u0,dfm_du0,dg_d_theta0',
     &              cos_theta0_res,Gamma_a_LH_u0,dfm_du0,dg_d_theta0
         write(*,*)'del_power=',del_power,'pow_dens=',pow_dens

         theta0_res=dacos(cos_theta0_res)

cSAP090827
c         idm=imax_chi_a
         idm=imax_chi

         write(*,*)'sin_theta0_res,sin_trap',sin_theta0_res,sin_trap
         write(*,*)'dsin(th0a_chi(imax_chi,n_radial0))',
     &              dsin(th0a_chi(imax_chi,n_radial0))

         if((theta0_res.le.th0a_chi(imax_chi,n_radial0)).or.
     &     (theta0_res.ge.pi-th0a_chi(imax_chi,n_radial0))) then          
c         if (sin_theta0_res.lt.sin_trap) then
c-----------passing particles
c            write(*,*)'theta0_res,u0',theta0_res,u0
            if (theta0_res.le.th0a_chi(imax_chi,n_radial0)) then  
cSAP090130
               if (i_chi_interpolation.eq.2) then
                  call  interpolation_chi(2,n_radial0,u0,theta0_res,
     &            chi,d_chi_du0,d_chi_dtheta0)
                 write(*,*)'interpolation: d_chi_du0,d_chi_dtheta0',
     &                    d_chi_du0,d_chi_dtheta0
               endif

cSAP090130
               if (i_chi_interpolation.eq.1) then
c-----------------spline interpolation   

                  d_chi_du0=terp2p_Sm(theta0_res,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)

                  d_chi_dtheta0=terp2p_Sm(theta0_res,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)

                  write(*,*)'d_chi_du0,d_chi_dtheta0',
     &                      d_chi_du0,d_chi_dtheta0
               endif     !if (i_chi_interpolation.eq.1)

            else
               theta0_odd=pi-theta0_res 
cSAP090130
               if (i_chi_interpolation.eq.2) then
                  call  interpolation_chi(2,n_radial0,u0,theta0_odd,
     &            chi,d_chi_du0,d_chi_dtheta0)
                  d_chi_du0=-d_chi_du0
                  write(*,*)'interpolation d_chi_du0,d_chi_dtheta0',
     &                    d_chi_du0,d_chi_dtheta0
               endif

               if (i_chi_interpolation.eq.1) then
c-----------------spline interpolation  
                  d_chi_du0=-terp2p_Sm(theta0_odd,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,0,1,1,imax_chi)

            
                  d_chi_dtheta0=terp2p_Sm(theta0_odd,u0,
     &            imax_chi,th0a_chi(1,n_radial0),
     &            nmax_chi,ua,
     &            chi_3d(1,1,n_radial0),chi_tt(1,1,n_radial0),
     &            chi_uu(1,1,n_radial0),
cSAP090827
c     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi_a)
     &            chi_ttuu(1,1,n_radial0),idm,1,0,1,imax_chi)
               endif
            endif
                          write(*,*)'d_chi_du0,d_chi_dtheta0',
     &                    d_chi_du0,d_chi_dtheta0

c            lh_cd_dens=lh_cd_dens+du*(u0**3/gamma)*
c     &                 cos_theta0_res*Gamma_a_LH_u0*
c     &                 (d_chi_du0-d_chi_dtheta0/(u0*cos_theta0_res))*
c     &                 dfm_du0/dg_d_theta0

c            del_cd=du*(u0**3/gamma)*
c     &                 cos_theta0_res*Gamma_a_LH_u0*
c     &                 (d_chi_du0-d_chi_dtheta0*sin_theta0_res/
c     &                 (u0*cos_theta0_res))*
c     &                 dfm_du0/dg_d_theta0

            lh_cd_dens=lh_cd_dens+du*u0**2*sin_theta0_res*
cSAP080415 in lambda it should be used |cos_theta0_res|
c     &           (u0/gamma*cos_theta0_res)*Gamma_a_LH_u0*
     &           (u0/gamma*dabs(cos_theta0_res))*Gamma_a_LH_u0*
     &           (d_chi_du0-
     &           (d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))*
     &                 dfm_du0/dg_d_theta0

            del_cd=du*u0**2*sin_theta0_res*
cSAP080415 in lambda it should be used |cos_theta0_res|
c     &           (u0/gamma*cos_theta0_res)*Gamma_a_LH_u0*
     &           (u0/gamma*dabs(cos_theta0_res))*Gamma_a_LH_u0*
     &           (d_chi_du0-
     &           (d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))*
     &                 dfm_du0/dg_d_theta0

           write(*,*)'del_cd=',del_cd
          
            
            write(*,*)'(d_chi_du0-'//
     & '(d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))',
     & (d_chi_du0-(d_chi_dtheta0/u0)*(sin_theta0_res/cos_theta0_res))   
             
            write(*,*)'pow_dens,lh_cd_dens=',pow_dens,lh_cd_dens        

         endif 
      enddo

      pow_dens=2.d0*pi*pow_dens
      lh_cd_dens=-2.d0*pi*LH_CD_dens
      write(*,*)'pow_dens',pow_dens
      write(*,*)'lh_cd_dens',lh_cd_dens
     
      if(pow_dens.ne.0.d0) then
        lh_cd_efficiency=lh_cd_dens/pow_dens
      else
        lh_cd_efficiency=0.d0
      endif

      write(*,*)'in CD_adj_LH_efficiency_1 lh_cd_efficiency=',
     &lh_cd_efficiency
         
c-----cln -Coulomb Ln
      arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense)	!temp KeV, dense 10**13/cm**3
      cln=24.d0-dlog(arg1)
      lh_cd_efficiency=lh_cd_efficiency*
     &3.84d0*temp_kev/(cln*dense)	       !  (A/m**2)/(W/m**3)
c      write(*,*)'temp_kev,dense,cln,lh_cd_efficiency',
c     &temp_kev,dense,cln,lh_cd_efficiency
      lh_cd_efficiency=lh_cd_efficiency*1.d-5  !  A/cm**2/(erg/(sec*cm**3)
      write(*,*)'temp_kev,dense,cln,lh_cd_efficiency',
     &temp_kev,dense,cln,lh_cd_efficiency
      return
      end

      subroutine CD_adj_LH_efficiency_test(n_parallel,n_perp,
     &psi_loc,theta_pol,e_y,e_z,
     &lh_cd_efficiency)
c-----------------------------------------------------------------
c     calculate LH CD efficiency using adj chi function
c-----------------------------------------------------------------       
      implicit none
      include 'param.i'
      include 'adj.i'

c-----input
      real*8      
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicular refractive index
     &psi_loc,             !polidal flux 
     &theta_pol,           !poloidal angle
     &length               !field line length normalized to  
c     &t,                  !is the temperature in units of m*unorm**2
                           !If velocities is normalized sqrt(T/m), t=1. 
      complex*16  
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1
c-----output
      real*8 
     &lh_cd_efficiency        
c-----local
      integer  n_radial_minus,n_radial_plus
      real*8 weight_radial,lh_cd_efficiency_minus,
     &lh_cd_efficiency_plus,unorm_minus,unorm_plus,
     &v_ph
 
c      write(*,*)' CD_adj_LH_efficiency e_y,e_z', e_y,e_z
c      write(*,*)' CD_adj_LH_efficiency before find_psis_points'
c      write(*,*)' n_radial_minus,n_radial_plus,weight_radial',
c     &             n_radial_minus,n_radial_plus,weight_radial
      v_ph=1.d0 !to test analytical LH flux CD_adj_LH_efficiency_anal_1
 
      call  find_psis_points(psi_loc,
     & n_radial_minus,n_radial_plus,weight_radial)

c      write(*,*)' CD_adj_LH_efficiency after find_psis_points'

      if((n_radial_minus.lt.npsi0).and.(n_radial_plus.gt.1)) then

         call CD_adj_LH_efficiency_1_test(n_parallel,n_perp,
     &   n_radial_minus,theta_pol,e_y,e_z,
     &   unorm_minus,
     &   lh_cd_efficiency_minus)

         write(*,*)'eff lh_cd_efficiency_minus',lh_cd_efficiency_minus

c         call CD_adj_LH_efficiency_anal_1(v_ph,n_parallel,
c     &   n_radial_minus,theta_pol,
c     &   unorm_minus,
c     &   lh_cd_efficiency_minus)
c         write(*,*)'anal lh_cd_efficiency_minus',lh_cd_efficiency_minus

         call CD_adj_LH_efficiency_1_test(n_parallel,n_perp,
     &   n_radial_plus,theta_pol,e_y,e_z,
     &   unorm_plus,
     &   lh_cd_efficiency_plus)
         write(*,*)'eff lh_cd_efficiency_plus',lh_cd_efficiency_plus

c         call CD_adj_LH_efficiency_anal_1(v_ph,n_parallel,
c     &   n_radial_plus,theta_pol,
c     &   unorm_plus,
c     &   lh_cd_efficiency_plus)
c          write(*,*)'anal lh_cd_efficiency_plus',lh_cd_efficiency_plus

         lh_cd_efficiency=lh_cd_efficiency_plus*weight_radial+
     &   lh_cd_efficiency_minus*(1.d0-weight_radial)

      else
        if (n_radial_plus.eq.1) then  
           call CD_adj_LH_efficiency_1_test(n_parallel,n_perp,
     &     1,theta_pol,e_y,e_z,
     &     unorm_plus,
     &     lh_cd_efficiency)
            write(*,*)'eff lh_cd_efficiency 1',lh_cd_efficiency

c           call CD_adj_LH_efficiency_anal_1(v_ph,n_parallel,
c     &     1,theta_pol,
c     &      unorm_plus,
c     &     lh_cd_efficiency)
c           write(*,*)'anal lh_cd_efficiency 1',lh_cd_efficiency
        else
           if (n_radial_minus.eq.npsi0) then
              call CD_adj_LH_efficiency_1_test(n_parallel,n_perp,
     &        npsi0,theta_pol,e_y,e_z,
     &        unorm_minus,
     &        lh_cd_efficiency)
           write(*,*)'eff lh_cd_efficiency npsi0',lh_cd_efficiency
c           call CD_adj_LH_efficiency_anal_1(v_ph,n_parallel,
c     &     npsi0,theta_pol,
c     &     unorm_plus,
c     &     lh_cd_efficiency)
c           write(*,*)'anal lh_cd_efficiency npsi0',lh_cd_efficiency

           endif
        endif
      endif

      write(*,*)'in CD_adj_LH_efficiency lh_cd_efficiency=',
     &lh_cd_efficiency

      return
      end

      subroutine DC_electric_field_adj_conductivity
c------------------------------------------------------------------------
c     calculate DC electric field conductivity sigma_E_adj(npsi0)
c     at all npsi0 flux surfaces psis(npsi0) using adj function chi
c     The result sigma_E_adj(npsi0) will be in adj_no_nml.i
c------------------------------------------------------------------------
c
c     Current density:
c     j=4*pi*integral{0<u_0<umax}integra{0<theta_0<theta_trap}
c            [chi*f_maxwell*u_0**3/sqrt(1+u**2/c**2)*
c             cos(theta_0)*sin(theta_0)]*d_u_0*d_theta_0/dens_maxwell
c
c     Electron density from Maxwellian distribytion function f_maxwell:
c     dens_maxwell= 4*pi*integral{0<u_0<umax}
c            [f_maxwell*u_0**2]*d_u_0
c
c     Conductivity:
c     sigma_E_adj=j/E
c
c     ELectric field:
c     E=T/Q
c  
c     T is the temperature
c
c     Q=integral{0<l<L}[(B_toroidal/B)/(2piR)]dl  is the safety factor
c
c     L is the field line length
c
c------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'adj.i'
c-----locals
cSAP090827
c      real*8 current_E_adj,cur(0:imax_chi_a-1),du,dth0,unorm,clight,
      real*8 current_E_adj,du,dth0,unorm,clight,
     & currrent_E_adj,
c     &chi_2d(0:imax_chi_a-1,0:nmax_chi_a-1),
c     &chi_2d(0:imax_chi_a,0:nmax_chi_a-1),
c     & fm(0:nmax_chi_a-1),
c     & th0a_loc(0:imax_chi_a-1),ua_loc(0:nmax_chi-1),
     & arg1,cln,mass_e,k_1_kev,charge_electron,pi,tau_ei,sigma_Spitzer,
     & tau_n,tau_e

      integer imax1,n_psi,i,n,istat

      real*8, pointer :: cur(:)       ! (0:imax_chi-1) 
      real*8, pointer :: chi_2d(:,:)  ! (0:imax_chi,0:nmax_chi-1),
      real*8, pointer :: fm(:)        ! (0:nmax_chi-1)
      real*8, pointer :: th0a_loc(:)  ! (0:imax_chi-1)      
      real*8, pointer :: ua_loc(:)    ! (0:nmax_chi-1)
      


c-----externals
      real*8 current,dens_maxwell,DC_E_field,rhopsi,current1,current2

      write(*,*)'DC_electric_field_adj_conductivity'
c      write(*,*)'nmax_chi,imax_chi',nmax_chi,imax_chi
c      write(*,*)'ua',ua

c---------------------------------------------------
c     allocate pointers 
c--------------------------------------------------
      allocate(cur(0:imax_chi-1),STAT=istat)
      call bcast(cur,0.d0,SIZE(cur))

      allocate(chi_2d(0:imax_chi,0:nmax_chi-1),STAT=istat)
      call bcast(chi_2d,0.d0,SIZE(chi_2d))

      allocate(fm(0:nmax_chi-1),STAT=istat)
      call bcast(fm,0.d0,SIZE(fm))

      allocate(th0a_loc(0:imax_chi-1),STAT=istat)
      call bcast(th0a_loc,0.d0,SIZE(th0a_loc))

      allocate(ua_loc(0:nmax_chi-1),STAT=istat)
      call bcast(ua_loc,0.d0,SIZE(ua_loc))

c------------------------------------------------      

      pi=4.d0*datan(1.d0)
      du=ua(2)-ua(1)
    
      imax1 = imax_chi
      if (mod(imax1,2) == 0) imax1 = imax1 + 1

      write(*,*)'imax1',imax1

      mass_e=9.1094d-28          !electron rest mass (g)
      k_1_kev=1.6022d-9          !egrs in 1 KeV      (erg)
      charge_electron=4.8032d-10 !electron charge (statcoulomb)

      do n_psi=1,npsi0
         write(*,*)'n_psi',n_psi   
         unorm=dsqrt(teme_chi(n_psi)/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
         clight=2.99792458d10/unorm

c-----------------------------------------------------
c        cln -Coulomb Ln
c-----------------------------------------------------
         arg1 = 1.d3/teme_chi(n_psi)*dsqrt(10.d0*dene_chi(n_psi)) !temp KeV,
                                                                  !dense 10**13/cm**3
         cln=24.d0-dlog(arg1)
         write(*,*)'cln',cln
c---------------------------------------------------------------------
c        tau_ei                   electron-ion momentum trancfer time
c---------------------------------------------------------------------
c        tau_ei=3/(4*sqrt(2*pi))*T**(3/2)*m**1/2/(Zeff*e**4*dens*Cln) [cgs]
c--------------------------------------------------------------------
         write(*,*)'zi_chi(n_psi)',zi_chi(n_psi)
       

         tau_ei=3.d0/(4.d0*dsqrt(2.d0*pi))*
     &   (1.6022d0*teme_chi(n_psi))*dsqrt(1.6022d0*teme_chi(n_psi))*
     &   dsqrt(9.1094d0)*
cSAP100112
c     &   (1.d-1)/
     &   dsqrt(1.d-1)/
     &   (zi_chi(n_psi)*4.8032d0**4*dene_chi(n_psi)*cln) ! sec
         write(*,*)'tau_ei',tau_ei

c         tau_ei=3.d0/(4.d0*dsqrt(2.d0*pi))*
c     &   (k_1_kev*teme_chi(n_psi))**1.5d0*dsqrt(mass_e)/
c     &   (zi_chi(n_psi)*charge_electron**4*dene_chi(npsi)*1.d13*cln) ! sec
c         write(*,*)'tau_ei',tau_ei

c-----------------------------------------------------------------
c        tau_n character time used in Karney code
c------------------------------------------------------------------
c        tau_n=4*pi*epsolin_0**2*T**3/2*m**1/2/(e**4*dens*Cln) [SI]=
c              =T**3/2*m**1/2/(4*pi*e**4*dens*Cln)  [cgs]
c-------------------------------------------------------------------
         tau_n=1.d0/(4.d0*pi)*
     &   (1.6022d0*teme_chi(n_psi))*dsqrt(1.6022d0*teme_chi(n_psi))*
     &   dsqrt(9.1094d0)*
cSAP100112
c     &   (1.d-1)/
     &   dsqrt(1.d-1)/
     &   (4.8032d0**4*dene_chi(n_psi)*cln) ! sec
         write(*,*)'tau_n',tau_n
c-------------------------------------------------------------------
c        from Killeen book 1.6022d-9   
         tau_e=3.44d5*
     &   (1.6022d0*teme_chi(n_psi))*dsqrt(1.6022d0*teme_chi(n_psi))*
cSAP100112
c     &   (1.d-1)/
     &   dsqrt(1.d-1)/
     &   (dene_chi(n_psi)*cln)            !sec
         write(*,*)'tau_e',tau_e
c---------------------------------------------
         sigma_Spitzer=2.d0*
     &   dene_chi(n_psi)*1.d13*charge_electron**2*tau_ei/mass_e ! sec*-1
         write(*,*)'sigma_Spitzer',sigma_Spitzer

         write(*,*)'sigma_Spitzer',sigma_Spitzer
c----------------------------------------------------------------------
c         write(*,*)'unorm,clight',unorm,clight
         do i=0,imax_chi-1
           th0a_loc(i)=th0a_chi(i+1,n_psi) 
c           write(*,*)'i,th0a_loc(i)',i,th0a_loc(i)
         enddo
         dth0= th0a_loc(1)-th0a_loc(0)

         do n=0,nmax_chi-1
            ua_loc(n)=ua(n+1)
c            write(*,*)'n,ua_loc(n)',n,ua_loc(n)
         enddo

c         write(*,*)'chi_3d(1,1,n_psi)',chi_3d(1,1,n_psi)

         do n=1,nmax_chi
           do i=1,imax_chi
             chi_2d(i-1,n-1)=chi_3d(i,n,n_psi)
           enddo
         enddo

c         write(*,*)'chi_2d(0,0)',chi_2d(0,0)
cSAP090827
cSAP090503 nmax_chi_a
c         call maxwel_adj(nmax_chi_a,ua_loc,clight,nmax_chi,fm,t)
         call maxwel_adj(nmax_chi,ua_loc,clight,nmax_chi,fm,t)

         write(*,*)'dens_averaged_ar(n_psi)',dens_averaged_ar(n_psi)
         write(*,*)'npsi0,n_psi',npsi0,n_psi
         write(*,*)'dene_chi(n_psi)',dene_chi(n_psi)
c         fm(:nmax_chi-1)=dene_chi(n_psi-1)*fm(:nmax_chi-1)/
c     &                   dens_averaged_ar(n_psi)

c         write(*,*)'fm',fm
c        write(*,*)'before chi_2d',chi_2d 
         !write(*,*)'before current chi_2d(0,nmax_chi-1)',chi_2d(0,nmax_chi-1) 
c         currrent_E_adj=current(du,ua_loc,clight,nmax_chi,dth0,
c     &                     th0a_loc,imax_chi,imax1,fm,chi_2d,cur) !q*dens*u_n
         currrent_E_adj=current1(du,ua_loc,clight,nmax_chi,dth0,
     &                     th0a_loc,imax_chi,imax1,fm,chi_2d,cur) !q*dens*u_n
c         write(*,*)'du,clight,nmax_chi,dth0,imax_chi,imax1',
c     &              du,clight,nmax_chi,dth0,imax_chi,imax1

c        write(*,*)'ua_loc',ua_loc
c        write(*,*)'th0a_loc',th0a_loc
c        write(*,*)'fm',fm 
c        write(*,*)'chi_2d',chi_2d 
         write(*,*)'currrent_E_adj',currrent_E_adj
c         write(*,*)'Q_safety_adj(n_psi)',Q_safety_adj(n_psi)

c         sigma_E_adj(n_psi)= currrent_E_adj/
c     &        Q_safety_adj(n_psi)


c         write(*,*)'Q_safety_adj(n_psi)',Q_safety_adj(n_psi)
c         write(*,*)'sigma_E_adj(n_psi)',sigma_E_adj(n_psi)

c----------------------------------------------------------------------
c
c  Spitzer conductivity
c----------------------------------------------------------------------
         sigma_E_adj(n_psi)= currrent_E_adj*bbar(n_psi)*zi_chi(n_psi)

         write(*,*)'bbar(n_psi),teme_chi(n_psi)',
     &              bbar(n_psi),teme_chi(n_psi)
         write(*,*)'currrent_E_adj*bbar(n_psi)*zi_chi(n_psi)',
     & sigma_E_adj(n_psi)

          sigma_E_adj(n_psi)=sigma_E_adj(n_psi)/
     &                       (3.d0*dsqrt(2.d0*pi))

         write(*,*)'sigma_E_adj(n_psi)',sigma_E_adj(n_psi)
c         write(*,*)'sigma_Spitzer',sigma_Spitzer

         rho_adj(n_psi)=rhopsi(psis(n_psi))
         write(*,*)'n_psi,sigma_E_adj(n_psi),rho_adj(n_psi)',
     &              n_psi,sigma_E_adj(n_psi),rho_adj(n_psi)
      enddo   !n_psi

c-------------------------------------------------------------
      deallocate(cur,STAT=istat)
      deallocate(chi_2d,STAT=istat)
      deallocate(fm,STAT=istat)
      deallocate(th0a_loc,STAT=istat)
      deallocate(ua_loc,STAT=istat)
c---------------------------------------------------------


      return
      end

      real*8 function averaged_density(nmax,imax,ua,th0a,fm,lam1a)
c-------------------------------------------------------------
c     calculates flux averaged density from maxwellian distriburtion fm
c     dens0=2pi*integral{0,pi}d_theta_0*sin(theta_0)
c               integral{0,infinity}u_0**2d_u_0
c               lambda(rho,u_0,theta_0)*f_maxwell(T,u_0)
c-------------------------------------------------------------
      implicit none
c-----input
      integer
     & nmax,                       !number of points in momentum mesh
     & imax                        !number of points in pitsch angle mesh
      real*8
     &ua(0:nmax - 1),              !momentum mesh 
     &th0a(0:imax - 1),            !pitch angle mesh         
     &fm(0:nmax-1),                ! relativistic maxwellian distribution
     &lam1a(0:imax)                !cos(theta_0)*integral[dl/cos(theta)]
                                   !cos(theta)=dsqrt(1-b(l)sin(theta_0)**2)
c-----locals
      integer i,n
      real*8 pi,du,dth0,sum_theta,sum_u

      pi=4.d0*datan(1.d0)
      
      du=ua(2)-ua(1)
      dth0=th0a(2)-th0a(1)
 
      write(*,*)'in fuction averaged_density nmax,imax',nmax,imax
c      write(*,*)'ua',ua
c      write(*,*)'th0a',th0a
!      write(*,*)'du,dth0,pi,imax',du,dth0,pi,imax  
      sum_theta=0.d0
      do i= 1, imax-1
!         write(*,*)'i,lam1a(i),th0a(i),dsin(th0a(i))',
!     &               i,lam1a(i),th0a(i),dsin(th0a(i))
         sum_theta=sum_theta+lam1a(i)*dsin(th0a(i))
      enddo
      sum_theta=sum_theta*2.d0 
      !write(*,*)'sum_theta',sum_theta
      !write(*,*)'nmax',nmax
      !write(*,*)'ua',ua
      !write(*,*)'fm',fm 

      sum_u = dot_product(ua(:nmax-1)**2,fm(:nmax-1))

      !write(*,*)'sum_u',sum_u     

      averaged_density=2.d0*pi*sum_theta*sum_u*du*dth0
      write(*,*)'averaged_density', averaged_density

      return
      end

      function current1 (du, ua, c, nmax, dth0, th0a, imax, imax1,
     1   fm, chi, cur)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(kind=dp) current
      integer nmax, imax, imax1
      real(kind=dp) du, c, dth0
      real(kind=dp), dimension(0:nmax - 1) :: ua
      real(kind=dp), dimension(0:imax - 1) :: th0a
      real(kind=dp), dimension(0:nmax - 1) :: fm
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: chi
      real(kind=dp), dimension(0:imax - 1) :: cur
      real(kind=dp)  current1 !SAP080101
      real(kind=dp), dimension(0:imax - 1) :: cur1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, i
      real(kind=dp) :: pi
ctest
      integer :: j
      real(kind=dp) :: chi_norm
C-----------------------------------------------

      pi = 4*datan(1.0d0)
      cur(:imax-1) = 0.d0
      cur1(:imax-1) = 0.d0

c      write(*,*)'current1 imax,imax1,nmax',imax,imax1,nmax

c      write(*,*)'in current1 chi(0,130)',chi(0,130) 

      do n = 0, nmax - 1
         cur(:imax-1) = cur(:imax-1) + chi(:imax-1,n)*fm(n)*ua(n)**3/
     1      dsqrt(1 + (ua(n)/c)**2)           
      end do

      do i=0,imax-1
        do n=0,nmax-1
          cur1(i) = cur1(i) + chi(i,n)*fm(n)*ua(n)**3/
     1      dsqrt(1 + (ua(n)/c)**2)
c          write(*,*)'i,n,ua(n),fm(n),chi(i,n)',i,n,ua(n),fm(n),chi(i,n)
c          write(*,*)'cur1(i)',cur1(i)  
        enddo
      enddo

c      do i=0,imax
c        write(*,*)'i,cur(i),cur1(i)',i,cur(i),cur1(i)
c      enddo   

      current1=0.d0  
      do i=0,imax-1
         current1=current1+cur1(i)*dsin(th0a(i))*dcos(th0a(i))
      enddo

      current=dot_product(dsin(th0a(:imax-1))*dcos(th0a(:imax-1)),cur(:
     1   imax-1)) 
c      write(*,*)'current,current1',current,current1
       
      current = 4*pi*du*dth0*current

      current1 = 4*pi*du*dth0*current1
      write(*,*)'current,current1',current,current1

      return
      end function current1

     
   



      subroutine read_iout3_chi
c-------------------------------------------------------------------
c     reads output file iout3 with calculated 
c     additional arrays for adj problem
c-------------------------------------------------------------------
      use kind_spec

      implicit none

c-----input      
      include 'param.i'
      include 'adj.i'

c-----local
      integer n_psi,n,nth
        
      real*8, dimension(nthp0  + 1) :: thetas
      real*8, dimension(0:nthp0 - 1) :: thetae,dla,ba

      real*8   psimx,psimn,psi,alenb, bmax, bmin
 
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      
      !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
       open(iout3, file='adjinp', status='old') !in read_out3_chi
      !endif  !On myrank=0   ! myrank=0

c      write(*,*)'read_iout3_chi iout3=',iout3
 
 2000 format(1x,2i8)
 2010 format(4(2x,e20.12))
 2020 format(5(2x,e20.12))

      !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
      read (iout3, 2000) nthp0, npsi0 !in read_out3_chi
      write(*,*)'read_iout3_chi nthp0, npsi0', nthp0, npsi0
      read (iout3, 2010) psimx, psimn
      write(*,*)'read_iout3_chi psimx, psimn', psimx, psimn
      read (iout3, 2020) (thetas(n),n=1,nthp0 + 1)
      do nth = 1, nthp0 + 1
         write(*,*)'nth,thetas(nth)',nth,thetas(nth)
      enddo
      read (iout3, 2020) (thetae(n),n=0,nthp0 - 1) !in read_out3_chi
      !endif  !On myrank=0   ! myrank=0

      do n=0,nthp0 - 1
cSAP110114
         write(*,*)'n,thetae(n)',n,thetae(n)
         thetae_1(n)=thetae(n)
         write(*,*)'thetae_1(n)',thetae_1(n)
      enddo

      do n_psi=1,npsi0 

         !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         read (iout3, 2010) psi !in read_out3_chi
         !endif  !On myrank=0   ! myrank=0
        
         write(*,*)' n_psi,psi', n_psi,psi
         write(*,*)'nthp0',nthp0

         !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         read (iout3, 2020) (dla(n),n=0,nthp0 - 1) !in read_out3_chi
         !endif  !On myrank=0   ! myrank=0

         do n=0,nthp0-1
           write(*,*)'n,dla(n)',n,dla(n)
         enddo

         !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         read (iout3, 2020) (ba(n),n=0,nthp0 - 1) !in read_out3_chi
         !endif  !On myrank=0   ! myrank=0

         do n=0,nthp0-1
           write(*,*)'n,ba(n)',n,ba(n)
         enddo
        
         do nth=0,nthp0-1
            ba_2d(nth,n_psi)=ba(nth)
         enddo    

         !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         read(iout3, 2010) alenb, bmax, bmin !in read_out3_chi
         !endif  !On myrank=0   ! myrank=0

c         write (*,*)'alenb, bmax, bmin',alenb, bmax, bmin

      enddo !n_psi

      !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
      close(iout3) !in read_out3_chi
      !endif  !On myrank=0   ! myrank=0

      return
      end


      subroutine absorbed_power_using_ql_flux(cnpar,cnper,z_m,r_m,phi,
     &fluxn,power,del_s_poloidal,
     &absorbed_power)

      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'cefield.i'    

c-----input
      real*8
     &cnpar,cnper,       !parallel and perpendicular refractive index
     &z_m,r_m,           !space ray coordinates [meter]
     &phi,               !toroidal angle or the ray point [radian]
     &power,             !power in ray channel [erg/sec]
     &del_s_poloidal     !poloidal length of the ray element []

      real*8 fluxn       !power flux at unit electric field |E|=1 
                         !flux=B~(i)*B(i)+E~(i)*d(omega*eps(i,j))/domega*E(j)
                         !fluxn=0.5*dreal(flux*dconjg(flux))*vgrpol 
                         !vgrpool is a poloidal group velocity
                         !normalized to clight   
c-----output          
      real*8
     &absorbed_power     !erg/sec the power abdsorbed at the ray element
   
c-----external
      real*8
     &psi_rho,
     &bmin_psi

c-----locals
      real*8
     &rho_loc,cos_theta_pol,sin_theta_pol,theta_pol,
     &psi_loc,unorm,
     &cd_efficiency,pow_dens_0_tilda,cd_dens,
     &b_pol,bmin,
     &clight,s_poloidal_tilda,
     &absorbed_power_dev_s_pol               !erg/[sec*cm] the power abdsorbed at the ray element 
                                             !of poloidal length delta_s_pol [cm]
                                             
c     &vgr_poloidal
c     &length_b

      integer n_theta_pol,
     &n_radial !can be arbitrary < n_radial_a

      pi=4.d0*atan(1.d0)
c      write(*,*)'pi',pi
      rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2) 

      write(*,*)'in absorbed_power_using_ql_flux rho_loc',rho_loc

      if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc = Inf
      cos_theta_pol=1.d0
      sin_theta_pol=0.d0
      else
      cos_theta_pol=(r_m-xma)/rho_loc
      sin_theta_pol=(z_m-yma)/rho_loc
      endif

c      write(*,*)'cos_theta_pol,sin_theta_pol',
c     &           cos_theta_pol,sin_theta_pol


      if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
      if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
      if (sin_theta_pol.ge.0.d0) then
         theta_pol=+dacos(cos_theta_pol)
      else  
cSAP081203
c         theta_pol=-dacos(cos_theta_pol)

c---------it should be 0 =< theta_pol< 2pi
         theta_pol=2.d0*pi-dacos(cos_theta_pol)
      endif 

c      write(*,*)'theta_pol',theta_pol

      if (theta_pol.lt.0.d0) then
         theta_pol=theta_pol+2.d0*pi
      endif

      if (theta_pol.gt.2.d0*pi) then
         n_theta_pol=theta_pol/(2.d0*pi)
c         write(*,*)'theta_pol,pi,n_theta_pol',
c     &              theta_pol,pi,n_theta_pol
         theta_pol=theta_pol-2.d0*pi*n_theta_pol
c         write(*,*)' theta_pol',theta_pol
      endif

      n_radial=1
      psi_loc=psi_rho(rho)

c      write(*,*)'in in absorbed_power_using_ql_flux rho',rho
c     write(*,*)'before QL_power_1'

      call QL_power_1(
     &psi_loc,
     &theta_pol,cnpar,cnper,cex,cey,cez,
     &unorm,
     &pow_dens_0_tilda)

c      write(*,*)'after QL_power_1'

      bmin= bmin_psi(psi_loc)
c-----lengt along B field: length_b [meters]
c      call calc_length_b(psi_loc,length_b)
c-----total group velocity
      clight=2.99792458d10              !light speed [cm/sec]

c      vgr_poloidal=dsqrt(vgr(1)**2+vgr(2)**2)! poloidal group velocity /clight
c      vgr_poloidal=vgr_poloidal*clight !sm/sec]

      s_poloidal_tilda= fluxn*clight/(8.d0*pi) 

c      write(*,*)'after CD_adj_efficiency_power_1'
c      write(*,*)'fluxn,s_poloidal_tilda',fluxn,s_poloidal_tilda,
c     &'[cm/sec]'
c      write(*,*)'bmod.bmin,bmod/bmin',bmod,bmin,bmod/bmin
c      write(*,*)'pow_dens_0_tilda',pow_dens_0_tilda,'[1/sec]'

      absorbed_power_dev_s_pol=power*(bmod/bmin)*
     &         (pow_dens_0_tilda/s_poloidal_tilda)
           
c      write(*,*)'power,absorbed_power_dev_s_pol',
c     &           power,absorbed_power_dev_s_pol

      absorbed_power=absorbed_power_dev_s_pol*del_s_poloidal ![erg/sec]
c      write(*,*)'del_s_poloidal,absorbed_power',
c     &           del_s_poloidal,absorbed_power
      return
      end





      subroutine QL_power_1(
     &psi_in,
     &theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &unorm, 
     &pow_dens_0_tilda)
c------------------------------------------------------
c     calculates bounce averaged power:  pow_dens 
c     at unit wave elecric field vector |E|=1
c------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one_nml.i'
      include 'adj.i'
      include 'write.i'
      
c-----input
 
      real*8  
     &theta_pol,           !poloidal angle [radian]    
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicular refractive index 
     &psi_in               !poloidal flux
      complex*16 
     &e_x,
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1

c-----output
      real*8
     &unorm,               !sqrt(T/m) cm/sec
     &cd_efficiency,
     &pow_dens_0_tilda,     ![1/sec], see manual: absorbed power calculation
c                          !Power density normalized to m*u_norm**2*n/tau_n 
     &cd_dens              !current drive density normalized to q*u_norm*n
                           !Here n is the density  
c-----external
      real*8 temperho,densrho,
     &bmax_psi,bmin_psi,psi_rho,rhopsi,y,b

c-----local
      integer n_radial0,n_harm_adj,
c     &i_resonance_curve_integration_method,
     &kmax,
     &i_calculate_CD     ! =0, calculate power only (no CD) 
                         ! using input arguments computed at
                         ! given space point
                         ! In this case n_radial0 can be arbitary
      real*8 u0,     
     &rho_small,           !small radius
     &temp_kev,            !electron temperatura [KeV]
     &dense,               !electron dencity 10**13/cm**3     
     &b_ratio,             !B/B_min
     &pi,bmin,bmax,deltb,psi,th0max,sin_trap,bmod,
     &d_chi_du0,d_chi_dtheta0,u0max,chi,
     &clight,             !light speed [cm/sec]
     &arg1,cln,
     &z,r,phi,            !bellow phi is set arbitratry phi=0
     &y_loc,              !electron y at point (z,r,phi)     
     &power_nharm_adj,CD_nharm_adj, !under integral functions
     

     &p_par_rl(2),         !for check max and min p_par
     &k_1_kev,             !egrs in 1 KeV      (erg)
     &charge_electron,     !electron charge (statcoulomb)
     &mass_e,              !electron mass [g]
     &theta_temperature,   !m_e*c^2/T_e
     &pow_dens,            
                           !using QL flux
     &omega                !2*pi*frqncy*1.d9 [1/sec]

      k_1_kev=1.6022d-9          !egrs in 1 KeV      (erg)
      charge_electron=4.8032d-10 !electron charge (statcoulomb)
      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]    
      mass_e=9.1094d-28         !electron rest mass (g)

c      write(*,*)'QL_power_1'


      rho_small=rhopsi(psi_in)
c      write(*,*)'rho_small',rho_small

      psi=psi_in

      temp_kev=temperho(rho_small,1)
      dense=densrho(rho_small,1) 
      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
      bmax= bmax_psi(psi)
      bmin= bmin_psi(psi)
	     
      deltb = (bmax-bmin)/bmin
c      write(*,*)'bmax,bmin,deltb',bmax,bmin,deltb   
      th0max = atan2(1.0d0,sqrt(deltb))
c      write(*,*)'th0max',th0max  mass_e=9.1094d-28    !electron rest mass (g)
c      k_1_kev=1.6022d-9                                !egrs in 1 KeV      (erg)
      sin_trap=dsqrt(bmin/bmax)

c      write(*,*)'sin_trap=',sin_trap,'dsin(th0max)=',dsin(th0max)

      pow_dens=0.d0
      cd_dens=0.d0

c      umin=clight/dsqrt(n_parallel**2-1.d0)/unorm ! normalized umin      
c      write(*,*)'umin',umin,'umin/dsqrt(2.d0)',umin/dsqrt(2.d0)
      
c      b_ratio=b_d_bmin(n_radial0,theta_pol)

      phi=0.d0
      call zr_psith(psi_in,theta_pol,z,r)

      bmod=b(z,r,phi)  
c      b_ratio=bmod/bmin
      b_ratio=dmax1(bmod/bmin,1.d0)
c      write(*,*)'theta_pol,b_ratio',
c     &theta_pol,b_ratio

      y_loc=y(z,r,phi,1)

c      i_resonance_curve_integration_method_adj=3
    
      do n_harm_adj=n_harm_adj_min,n_harm_adj_max
c-------------------------------------------------------------------
c       the loop over used gyro-harmonics 
c----------------------------------------------------------------- 
       write(*,*)'QL_power_1 n_harm_adj=',n_harm_adj

c     & write(*,*)'i_resonance_curve_integration_method_adj=',
c     & i_resonance_curve_integration_method_adj


c---------------------------------------------------------------------
c       for check min,max p_parallel
c       roots p_par_rl(2) at resonance curve: parallel momentum divided by (mc)
        call root_res(0.d0,n_parallel,n_harm_adj,y_loc,kmax,p_par_rl)

c        write(*,*)'p_perp=0,n_harm_adj,n_parallel,y_loc,kmax,p_par_rl',
c     &   n_harm_adj,n_parallel,y_loc,kmax,p_par_rl

c------------------------------------------------------------------------
c        write(*,*)'before call intgr_relt_power_cd_adj'
c        write(*,*)'n_radial0,b_ratio',n_radial0,b_ratio
c        write(*,*)'y_loc,sin_trap,temp_kev',y_loc,sin_trap,temp_kev
c        write(*,*)'e_x,e_y,e_z',e_x,e_y,e_z
c        write(*,*)'n_harm_adj,n_parallel,n_perp',
c     &             n_harm_adj,n_parallel,n_perp
c        write(*,*)'unorm,umax',unorm,umax
c        write(*,*)'n_relt_intgr_adj',n_relt_intgr_adj
c        write(*,*)'i_resonance_curve_integration_method_adj',
c     &             i_resonance_curve_integration_method_adj
c        write(*,*)'epsi_adj',epsi_adj
c        write(*,*)'before intgr_relt_power_cd_adj'

        i_calculate_CD=0 ! calculate power only (no CD) 
                         ! using input arguments computed at
                         ! given space point
                         ! In this case n_radial0 can be arbitary
  
        call intgr_relt_power_cd_adj(i_calculate_CD,
     &  n_radial0,b_ratio,
     &  y_loc,sin_trap,temp_kev, 
     &  e_x,e_y,e_z,  
     &  n_harm_adj,n_parallel,n_perp, 
     &  unorm,umax,
     &  n_relt_intgr_adj,i_resonance_curve_integration_method_adj,
     &  epsi_adj,
     &  power_nharm_adj,CD_nharm_adj)

c        write(*,*)'after intgr_relt_power_cd_adj power_nharm_adj',
c     &             power_nharm_adj

        pow_dens=pow_dens +power_nharm_adj
        cd_dens=cd_dens+CD_nharm_adj

        write(*,*)'power_nharm_adj',power_nharm_adj
        write(*,*)'Cd_nharm_adj',CD_nharm_adj

      enddo !n_harm

      pow_dens=2.d0*pi*pow_dens
      cd_dens=-2.d0*pi*cd_dens

c      write(*,*)'pow_dens',pow_dens
c      write(*,*)'cd_dens',cd_dens
     
      if(pow_dens.ne.0.d0) then
        cd_efficiency=cd_dens/pow_dens
      else
        cd_efficiency=0.d0
      endif

      write(*,*)'in CD_adj_efficiency_1 cd_efficiency=',
     &cd_efficiency
         
c-----cln -Coulomb Ln
      arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense)	!temp KeV, dense 10**13/cm**3
      cln=24.d0-dlog(arg1)

c-----adj efficiency is normalized to eta_0_adj=e*tau_n/(m*u_t) 
c     tau_n=T**(3/2)*sqrt(m)/(4*pi*e**4*density*cln)
c     u_t=cqrt(T/m)
c     eta_0_adj=T/[4*pi*e**3*n*cln]|cgs=T_kev/(n_13*cln)*
c           [k_1_kev/(4*pi*charge_electron**3*10**13)]           !cgs
c      coef=1.6022d+8/(4*pi*4.8032**3)          !cgs  (statampere/cm**2)/(erg/(sec*cm**3))
c      coef=1.6022d-1/(4*pi*4.8032**3*3)=3.8352d-5!     (A/cm**2)/(erg/(sec*cm**3))
c      eta_0_adj=(temp_kev/(cln*dense))*coef

      cd_efficiency=cd_efficiency*
     &3.84d0*temp_kev/(cln*dense)	       !  (A/m**2)/(W/m**3)
c      write(*,*)'temp_kev,dense,cln,cd_efficiency',
c     &temp_kev,dense,cln,cd_efficiency
      cd_efficiency=cd_efficiency*1.d-5  !  (A/cm**2)/(erg/(sec*cm**3))
      write(*,*)'temp_kev,dense,cln,cd_efficiency',
     &temp_kev,dense,cln,cd_efficiency

c---------------------------------------------------------------------
cSAP081201
c     
c     The used relativistic Maxwellian distribution was normalized
c     to  unit density electron density:
c     4*pi*integral{f_m*(p/mc)^2*d(p/mc)}=1
c---------------------------------------------------------------------
  
      omega=2.d0*pi*frqncy*1.d9 ![1/sec]
c      write(*,*)'frqncy,omega',frqncy,omega
c---- power_dens_0_tilda has dimension [1/sec]
c      write(*,*)'pow_dens',pow_dens
c      write(*,*)'dense*1.d13',dense*1.d13
c      write(*,*)'unorm,clight,unorm/clight',unorm,clight,unorm/clight
      pow_dens_0_tilda=pow_dens*
     &dense*1.d13*      !the power was calculated for unit density
                        !This term transforms power to local density
     &(unorm/clight)*   !it was the integration by (p_perp/mc)*d(p_perp/mc)
     &charge_electron**2/mass_e/omega !in subroutine intgr_relt_power_cd_adj
                        !This transforms this integration to
                        !(p_perp/(m*unorm)*)d(p_perp/(m*unorm))
                        !and uses all other normalizations
                        !in under the integral term
c      write(*,*)'charge_electron**2/mass_e/omega [cm^3/sec] ',
c     &           charge_electron**2/mass_e/omega
c       write(*,*)'pow_dens_0_tilda',pow_dens_0_tilda,'[1/sec]'
      return
      end

     

      subroutine calc_length_b(psi,length_b)
c-----calculate length at flux surface along B field line
c     when the poloidal angle is changed at 2*PI
      implicit none 
      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'three.i'
      include 'rho.i'
c      real*8 length_b(npsi) !length along B
c-----input
      real*8 psi !polidal flyx
c-----output
      real*8  length_b ! length of B field line
c-----externals
      real*8 b
c-----locals
      integer j,i
      real*8 theta,z,r,phi,b_pol,dl_pol,dl_b,psi_loc
      logical first

c----------------------------------------------------------
c     to use spline functions from zcunix.f:
c     coeff1 and terp1
      integer i1p(2)
      real*8, dimension(1:3*npsi+1) :: work_l
      integer itabl(3)
      real*8  tabl(3)

      data first /.true./
      save first

          
      if (first) then
     
         length_b_ar(1)=0.d0  !length along B field
    
         phi=0.d0   
         do j=2,npsi           
            psi_loc=arpsi(j)            
            length_b_ar(j)=0.d0  !initialization of the length 
                           !along B field
	    do i=1,nteta
	       theta=arteta(i)
               call zr_psith(psi_loc,theta,z,r)
               bmod=b(z,r,phi)
               b_pol=dsqrt(bz**2+br**2)

c               write(*,*)'j,i,bmod,b_pol,bmod/b_pol',
c     &                    j,i,bmod,b_pol,bmod/b_pol

               dl_pol=dsqrt((rpsi(j,i+1)-rpsi(j,i))**2+
     &             (zpsi(j,i+1)-zpsi(j,i))**2)
               dl_b=dl_pol*bmod/b_pol        
               length_b_ar(j)=length_b_ar(j)+dl_b  ![meter]
            enddo !i 
           
         enddo !j
 
c------------------------------------------------------------------
c        spline coefficient calculation for length_b(psi) using spline
c        from zcunix.f
c---------------------------------------------------------------
         i1p(1)=4 !first derivatives at the left boundary will be fitted

         i1p(2)=4 !first derivatives at the right boundary will be fitted

         call coeff1(npsi,arpsi,length_b_ar,d2_length_b_psi,i1p,1,
     &               work_l)
 
         first = .false.
      else
         if (psi.lt.psimag) then
            length_b=0.d0
            goto 10
         endif
c---------------------------------------------------
c        b field line calculation using spline from zcunix.f
c---------------------------------------------------
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0
          
         call terp1(npsi,arpsi,length_b_ar,d2_length_b_psi,
     &             psi,1,tabl,itabl)
         
         length_b=tabl(1)
      endif

 10   continue
      return
      end
     



