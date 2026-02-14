
      subroutine emiss_coef(aK,fluctcur,ce,flux,v_group,omega,nr,
     +iabsorp_collisional,coll_mult,tempe,dense,zeff,
     +al_emis,j_emis)  
c     al_emis (1/cm)
c     j_emis  (erg*sec/cm**3)
c      
c-----calculates the emission coefficient j_emis(erg*sec/cm**3) 
c     and absorption al_emis (1/cm) 
c
c     see the article
c     R.W.Harvey at al., Phys.Fluids B 5(2),Febrary 1993, page 446
c     j_emis is the power radiated by the plasma 
c     per unit volume-radian frequency-steradian,
c
c     formula (3)
c     j_emis=pi*Nr**2(omega/clight)**2*vectorE^*tensorG*vectorE/S
c
c     vector S(omega,k) is the flux energy density per frequency and
c     per unit volume of k space,
c
c     formula(4)
c     S=vectorVgroup*U   
c
c     vectorVgroup is the wave group velocity
c
c     U is the spectral energy density
c
c     formula (5)
c     U=1/(8*pi)(vectorB^*vectorB+
c     +vectorE^*d(omega*eps_hermitian)/domega*vectorE) 
c     
c     al_emis is absorption coefficient, the inverse damping length
c     along the ray of rf wave energy
c 
c     formula (2)
c     al_emis=(omega/(4*pi))vectorE^*eps_anti-herm*vectorE/abs(S) 
c 
c     G is the correlation tensor    formula(7)
c
c     if iabsorp_collisional=1 it will add the collisional absorption
c                              al_emis_collisional=nu_ei/(Vgroup) (1/cm)
c                              in the absorption coeficient
c                              nu_ei is the electron-ion collision frequency
c     coll_mult  = multiplier for coll absorption expression (default=1.0)
c------------------------------------------------------------------
      implicit none

c-----input
      double complex aK(3,3)      ! anti-hermitian relativistic tensor
      double complex fluctcur(3,3)! correltion tensor G for the fluctuating current(erg)
      double complex ce(3)        ! wave (normalized) electric field in Stix frame
c     wave energy flux for the normalized field |E|=1
c     flux=(vectorB*vector_B^*+vector_E*d(omega*tensor_eps_herm)/dw*vector_E)
      double complex flux         ! wave energy flux for the normalized field |E|=1
      double precision omega      ! wave angle frequency (1/sec) =2*pi*f
      double precision v_group    ! wave group velocity [cm/sec]
      double precision nr         ! ray refructive index
      integer iabsorp_collisional ! =0 do not calculate collisional absorption
                                  ! =1 to calculte collisonal absorption
      double precision coll_mult  ! multiplier for coll absorption
      double precision tempe,     !electron temperature  (keV)  
     & dense,                     !electron density   (10*13/cm**3	 
     & zeff                       !plasma charge
c-----output
      double precision al_emis !(1/sm)
      double precision j_emis  !(erg*sec/sm**3)
      
c-----local
      double complex cec(3),cal_emis,cj_emis
      integer i,j
      double precision s,pi,clight,absflux,p
      double precision arg,cln,r_nu_ei,al_emis_collisional

      pi=4*datan(1.d0)
      clight=2.99792458d10          !cm/sec
c--------------------------------------------------------------
c     ce - wave electric field
c     cec- comlex conjugate wave electric field
c----------------------------------------------------------
      do i=1,3
         cec(i)=dconjg(ce(i))
      enddo
    
      p=flux*dconjg(flux)    
      absflux=dsqrt(p)
c      write(*,*)'emiss_coeff flux,absflux',flux,absflux

c     s is  the energy flux density per frequency and per unit volume in k space
      

      s=dabs(v_group*absflux/(8.d0*pi))
      s=s*4.d0*pi 
         
c      write(*,*)'emiss_coeff v_group,absflux,s',v_group,absflux,s
            
      cal_emis=dcmplx(0.0d00,0.0d00)
      cj_emis=dcmplx(0.0d00,0.0d00)

c      write(*,*)'emiss_coef aK(i,j)',aK(1,1),aK(1,2),aK(1,3)
c      write(*,*)aK(2,1),aK(2,2),aK(2,3)
c      write(*,*)aK(3,1),aK(3,2),aK(3,3)

c      write(*,*)'emiss_coef fluctcur',fluctcur(1,1)*omega,
c     +fluctcur(1,2)*omega,fluctcur(1,3)*omega
c      write(*,*)fluctcur(2,1)*omega,
c     +fluctcur(2,2)*omega,fluctcur(2,3)*omega
c      write(*,*)fluctcur(3,1)*omega,
c     +fluctcur(3,2)*omega,fluctcur(3,3)*omega

c      write(*,*)'emiss_coef ce',ce
 
      do i=1,3          
         do j=1,3
 	    cal_emis=cal_emis+cec(i)*aK(i,j)*ce(j)
            cj_emis =cj_emis +cec(i)*fluctcur(i,j)*ce(j)
c            write(*,*)'emiss_coef i,j',i,j
c            write(*,*)'cec(i),ce(j),aK(i,j)',cec(i),ce(j),aK(i,j)
c            write(*,*)'cec(i)*aK(i,j)*ce(j)',cec(i)*aK(i,j)*ce(j)
c            write(*,*)'fluctcur(i,j)',fluctcur(i,j)
c            write(*,*)'cj_emis',cj_emis
         enddo
      enddo

c      write(*,*)'emiss_coef 0 cal_emis',cal_emis
c      write(*,*)'emiss_coef 0 cj_emis',cj_emis

      cal_emis=cal_emis*omega/(4.d0*pi)/s
      cj_emis=cj_emis*pi*(nr*omega/clight)**2/s

c      write(*,*)'omega,s,omega/(4.d0*pi)/s',omega,s,omega/(4.d0*pi)/s
      
c      write(*,*)'nr,clight,pi*(nr*omega/clight)**2/s',
c     +nr,clight,pi*(nr*omega/clight)**2/s

c      write(*,*)'emiss_coef  cal_emis',cal_emis
c      write(*,*)'emiss_coef  cj_emis',cj_emis

      p=cal_emis*dconjg(cal_emis)  
      al_emis=dsqrt(p)
      al_emis=al_emis*4.d0*pi
c      write(*,*)'emiss_coeff cal_emis,al_emis',cal_emis,al_emis    
      p=cj_emis*dconjg(cj_emis)
      j_emis=dsqrt(p)*omega
      j_emis=j_emis*4.d0*pi
c      write(*,*)'emiss_coeff nr,cj_emis,j_emis',nr,cj_emis,j_emis  
      if (iabsorp_collisional.eq.0) then
         al_emis_collisional=0.d0
      endif
      if (iabsorp_collisional.eq.1) then
c-------------------------------------------------------------------------------------
c        collisional absorption calculations
c-------------------------------------------------------------------------------------
c--------Coulomb Logarithm = cln      
         arg=(1.d3/tempe)*dsqrt(10.d0*dense) !tempe KeV, 
                                             !dense 10**13/sm**3   
         cln=24.d0-dlog(arg)
c--------electron-ion collision rate =r_nu_ei [1/sec] 
         r_nu_ei=dense/(tempe*dsqrt(tempe))*zeff*cln*0.919d3 ![1/sec]
         al_emis_collisional=2.d0*r_nu_ei/(v_group)*coll_mult
      endif
c      write(*,*)'al_emis,al_emis_collisional',
c     &al_emis,al_emis_collisional
c-----total absorption
      al_emis= al_emis+al_emis_collisional
 
      return
      end


      subroutine emis_tens(X,Y,T_kev,nll_in,np_in,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +i_resonance_curve_integration_method,epsi,
     +i_fkin,r,z,phi,
     +aK,fluctcur)
c-----calulates the anti-hermitian relativistic dielectric tensor aK
c     for electron plasma
c     and the correlation tensor fluctcur=G  for the fluctuating current density
c
c     see the article
c     R.W.Harvey at all Phys.Fluids B 5(2),Febrary 1993, page 446
c     G is the correlation tensor    formula(7)
c     aK  formula (6)
c
c     INPUTS:
c
c      X = (fpe/f)**2
c      Y = fce/f
c      T_kev  - electron temperature
c      nll_in - parallel index of refraction N.
c      np_in  - perpendicular refractive index N
c      n_relt_harm1  - min number of used cyclotron jn harmonics,
c      n_relt_harm2  - max number of used cyclotron jn harmonics,
c                      n_relt_harm1 <= jn <= n_relt_harm2 
c      n_relt_intgr - the number of points for the integration over p_perp
c      i_fkin =0 the usage of the analytic relativistic distribution
c             =1 the usage of the numerical 3D distribution from diskf or netcdfnm.nc 
c                written be CQL3D code or created analytically at mesh points
c      r      - the major radius (it is used for i_fkin=1)
c      z      - the vertical coordinate  (it is used for i_fkin=1)
c      phi    - toroidal angle
c----------------------------------------------------------------------
c     i_resonance_curve_interation_method
c
c     i_resonance_curve_integration_method=1 angle integration
c     now it works for ellipse case only                             c
c
c     i_resonance_curve_integration_method=2 rectangle formula
c     for p_perp integration
c     
c     i_resonance_curve_integration_method=3 trapezoidal formula
c     for p_perp integration
c
c     i_resonance_curve_integration_method=4 !adaptive Simpson 
c     for p_perp integration
c
c     epsi is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------ 
c     OUTPUT:
c      aK(3,3):  the nine components of the anti-hermition
c                part dielectric relativistic tensor
c                evaluated at (X,Y,Te,nll,np,n)
c      
c      fluctcur(3,3) the correlation tensor fluctcur=G (egr)  
c                    for the fluctuating current density
c      

      implicit none
c     input 
      double precision X,Y,T_kev,nll_in,np_in,epsi
      integer n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,
     &i_resonance_curve_integration_method

      double precision r,z,phi

c     output
      double complex aK(3,3),fluctcur(3,3)
      
c     local
      double precision c,mass_e,k_1_kev,nll,np,nlls,p_perp0,theta,pi,
     .dens,xpidens,xpidensmcs,p,psi,rho
      double complex integral(3,3)
      double complex fluctcur_n(3,3)  

      integer jn,ires,jn_negative

c     external zeroK, npnllmin, ec_cond, intgr_rl,fdens_fdist
      double precision fdens_fdist,fpsi,rhopsi,densrho
 
      c =2.99792458d10          !speed of light     (cm/sec)
      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)      

      theta=mass_e*c**2/(k_1_kev*T_kev)
      pi=4.d0*datan(1.d0)

      call npnllmin(nll_in,np_in,nll,np)
      nlls = nll**2

c-----initialization aK=0
      call zeroK(aK)
      call zeroK(fluctcur)      
c-----the loop over the cyclotron harmonics

      write(*,*)'emission.f in emis_tens'

      do jn=n_relt_harm1,n_relt_harm2               
        jn_negative=-jn
c-------control the EC resonance condition
cSm060315         
c        call ec_cond(jn,Y,nll,ires,p_perp0)
        call ec_cond(jn_negative,Y,nll,ires,p_perp0)

c        write(*,*)'emission.f emis_tens jn,ires,p_perp0',jn,ires,p_perp0
        
        if(ires.eq.0) then
c---------no resonace
          goto 10
        endif

cSm060315
c       call intgr_emis(jn,nll,np,Y,theta,ires,p_perp0,n_relt_intgr,     
        call intgr_emis(jn_negative,nll,np,Y,theta,ires,p_perp0,
     +  n_relt_intgr,i_resonance_curve_integration_method,epsi,
     +  i_fkin,r,z,phi,
     +  integral,fluctcur_n)
c        write(*,*)'emis_tens after integr_emis integral',integral
c        write(*,*)'emis_tens after integr_emis fluctcur_n',fluctcur_n

        aK(1,1)=aK(1,1)+integral(1,1)
        aK(1,2)=aK(1,2)+integral(1,2)
        aK(1,3)=aK(1,3)+integral(1,3)
        aK(2,2)=aK(2,2)+integral(2,2)
        aK(2,3)=aK(2,3)+integral(2,3)
        aK(3,3)=aK(3,3)+integral(3,3)

        fluctcur(1,1)=fluctcur(1,1)+fluctcur_n(1,1)
        fluctcur(1,2)=fluctcur(1,2)+fluctcur_n(1,2)
        fluctcur(1,3)=fluctcur(1,3)+fluctcur_n(1,3)
        fluctcur(2,2)=fluctcur(2,2)+fluctcur_n(2,2)
        fluctcur(2,3)=fluctcur(2,3)+fluctcur_n(2,3)
        fluctcur(3,3)=fluctcur(3,3)+fluctcur_n(3,3)


10      continue      
      enddo   !jn

      if (i_fkin.eq.0) then       
        dens=1.d0
      else      
c        dens=fdens_fdist(r,z,phi)
cSm011228 call the density of the bulk plasma to reduce the CPU time
        psi=fpsi(r,z)
        rho=rhopsi(psi)
        dens=densrho(rho,1)
c        write(*,*)'emision_tens dens',dens
c        dens=1.d0
      endif
      
c-----normalization for the unit electron density
      xpidens=X*pi/dens 
      
      aK(1,1)=-xpidens*aK(1,1)
      aK(1,2)=-xpidens*aK(1,2)
      aK(1,3)=-xpidens*aK(1,3)
      aK(2,2)=-xpidens*aK(2,2)
      aK(2,3)=-xpidens*aK(2,3)
      aK(3,3)=-xpidens*aK(3,3)
     
c      aK(2,1)=-aK(1,2)
c      aK(3,1)= aK(1,3)
c      aK(3,2)=-aK(2,3)
cSm060705  

      aK(2,1)=conjg(aK(1,2))
      aK(3,1)=conjg(aK(1,3))
      aK(3,2)=conjg(aK(2,3))

      xpidensmcs=X*pi/(2.d0*pi)**5*(mass_e*c**2)/dens

      fluctcur(1,1)=fluctcur(1,1)*xpidensmcs
      fluctcur(1,2)=fluctcur(1,2)*xpidensmcs
      fluctcur(1,3)=fluctcur(1,3)*xpidensmcs
      fluctcur(2,2)=fluctcur(2,2)*xpidensmcs
      fluctcur(2,3)=fluctcur(2,3)*xpidensmcs
      fluctcur(3,3)=fluctcur(3,3)*xpidensmcs
      
      fluctcur(2,1)=-fluctcur(1,2)
      fluctcur(3,1)= fluctcur(1,3)
      fluctcur(3,2)=-fluctcur(2,3)
      

      return
      end

      subroutine intgr_emis(n,nll,np,Y,theta,ires,p_perp0,
     +n_relt_intgr,i_resonance_curve_integration_method,epsi,
     +i_fkin,r,z,phi,
     +integral,fluctcur_n)
c-----------------------------------------------------------------
c     calculates the matrix: double complex integrals
c     for the relativistic electron plasma 
c     I(n)^=integral(0<=p_perp<=p_perp0)g_n, g_n={sum_k=1,2,G_nk(p_perp)
c     fluctcar(n)^=integral(0<=p_perp<=p_perp0)g_fluctcur_n.
c               g_fluctcur_n={sum_k=1,2,g_fluctcur_nk(p_perp)
c     for the EC harmonic with number 'n'
c----------------------------------------------------------------- 
c     input
c       n        - EC harmonic number
c       nll      - N_parallel
c       np       - N_perpendicular
c       Y        = omega_ce/omega for the electron rest mass
c       theta    = mc**2/T
c       ires     =0 no resonance,=1 ellipse,=2 parabola,=3 hyperbole 
c       n_relt_intgr - the number of points for the integration over p_perp
c       p_perp0  - max value of the perpendicular momentum divided by mc
c                  on the resonanse ellipse
c       i_fkin   =0 the usage of the analytic relativistic distribution
c                =1 the usage of the numerical 3D distribution from diskf 
c                   written be CQL3D code
c       r        - the major radius (it is used for i_fkin=1)
c       z        - the vertical coordinate  (it is used for i_fkin=1)
c       phi        - toroidal angle (it is used for i_fkin=1)
c     i_resonance_curve_interation_method
c
c     i_resonance_curve_integration_method=1 angle integration
c     now it works for ellipse case only                             c
c
c     i_resonance_curve_integration_method=2 rectangle formula
c     for p_perp integration
c     
c     i_resonance_curve_integration_method=3 trapezoidal formula
c     for p_perp integration
c
c     i_resonance_curve_integration_method=4 !adaptive Simpson 
c     for p_perp integration
c
c     epsi is the accuracy used in integration by adaptive Simpson
c    
c------------------------------------------------------------ 
c     output
c       integral(3,3) double complex integral from G_nk over 0<=p_perp=<p_perpmax
c       fluctcur(3,3) double complec integral from fluct_cur_nk over 0<=p_perp=<p_perpmax
c----------------------------------------------------------------    
c      The integration method is specified by the variable
c      i_resonance_curve_integration_method.
c      Now this variable is set inside this subroutine.
c----------------------------------------------------------------
      implicit none
c     input
      double precision nll,np,Y,theta,p_perp0,epsi
      integer n,ires,n_relt_intgr,i_fkin,
     &i_resonance_curve_integration_method

      double precision r,z,phi
      double precision vnormloc,massloc ! cm/sec, g
      COMMON /dskin1/vnormloc,massloc
c     output
      double complex integral(3,3)
      double complex fluctcur_n(3,3) 
c     local
cSm060725
      double precision eps,p_permax,h,p_perp,p,p_t,clight,p_int,
c     double precision eps,p_permax,h,p_perp,p,p_t,jmax,clight,p_int,
     &vmax_d_vt
      double complex i,g_n(3,3),g_fluctcur_n(3,3)
cSm060725
      integer j,jmax

c-----for integration along ellipse by angle
      double precision rme,rtem0,vper,xint,cs,sn,thet1,thet2,p_par,
     & p_par_min,p_par_pl,p_perp_min,p_perp_pl,v0dc,vmax1dc,vpar0dc,
     +vpar1dc,vpar2dc,dth,vt0,cdvt,tem0,dvper,vmaxdt

c     external
c     ec_cond, zeroK, g_n 
      double precision temperho

      i = ( 0.0d0,1.0d0)        !imaginary number     
      jmax=n_relt_intgr

cSm060310
      if(ires.eq.0)goto10 ! no resonance
    
      vmax_d_vt=10.d0
      call p_perp_max_calc(i_fkin,theta,n,Y,nll,vmax_d_vt,
     &vnormloc,p_permax,ires)

      if(ires.eq.4) then
c--------the resonance curve is outside the grid
         call zeroK(integral)
         call zeroK(fluctcur_n)
         goto 10
      else
         h=p_permax/(n_relt_intgr)   ! step of integration over p_perp
      endif

c     calculations of the integrals over p_perp
      call zeroK(integral)
      call zeroK(fluctcur_n)


c----------------------------------------------------------
c     i_resonance_curve_interation_method
c
c     i_resonance_curve_interation_method=1 angle integration
c     now it works for ellipse case only                             c
c     i_resonance_curve_interation_method=2 rectangle formula
c     for p_perp integraion
c     
c     i_resonance_curve_interation_method=3 trapezoidal formula
c     for p_perp integraion
c
c     i_resonance_curve_integration_method=4 !adaptive Simpson 
c     for p_perp integration
c------------------------------------------------------------ 
     
      goto (1,2,3,4) i_resonance_curve_integration_method

 1    continue      
cSm060327 
      if(dabs(nll).ge.1.d0) goto 3 !to trapezoidal formula 
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     begin
c-------------------------------------------------------------------- 
      clight=2.99792458d10
      rme=9.1094d-28
      call get_rtem0_from_one(rtem0)
      tem0=temperho(0.d0,1)*rtem0
      vt0=dsqrt(2.d0*tem0*1.6022d-9/rme) !(cm/sec)the  central 
                                         !thermal velocity
      cdvt=clight/vt0
      vmaxdt=dsqrt(rtem0)
      call ec_condh(n,Y,nll,vmaxdt/cdvt,ires,v0dc,vmax1dc,vpar0dc,
     +vpar1dc,vpar2dc,thet1,thet2)
    
c     write(8,*)'theta integration boundaries: thet1,thet2',
c     thet1,thet2      
    
c-----Subdivide theta range of integration
      dth=(thet2-thet1)/(jmax-1)
      do j=1,jmax
         xint=thet1+(j-1)*dth
         cs=dcos(xint)
         sn=dsin(xint)
         p_par=(vpar0dc-dabs(v0dc)*cs) !vper/Vt
         p_perp=vmax1dc*sn             !vpar/Vt
       
         if(p_perp.lt.1.d-3) p_perp=1.d-3 !to avoid singularity  

         call g_n_emis_theta(p_perp,p_par,y,nll,np,theta,n,i_fkin,
     &   r,z,phi,g_n,g_fluctcur_n)
          
         if((j.ne.jmax).and.(j.ne.1)) then
           sn=dsin(xint+dth)
           p_perp_pl=vmax1dc*sn                !vper/c        
           sn=dsin(xint-dth)
           p_perp_min=vmax1dc*sn   
           h=0.5d0*(dabs(p_perp_pl-p_perp)+dabs(p_perp-p_perp_min))
         else
           if(j.eq.1) then
             sn=dsin(xint+dth)
             p_perp_pl=vmax1dc*sn              !vper/c
             h=0.5d0*dabs(p_perp_pl-p_perp)
           else
             !j=jmax
             sn=dsin(xint-dth)
             p_perp_min=vmax1dc*sn             !vper/c
             h=0.5d0*dabs(p_perp-p_perp_min)
           endif 
         endif

         integral(1,1)=integral(1,1)+h*g_n(1,1)
         integral(1,2)=integral(1,2)+h*g_n(1,2)
         integral(1,3)=integral(1,3)+h*g_n(1,3)
         integral(2,2)=integral(2,2)+h*g_n(2,2)
         integral(2,3)=integral(2,3)+h*g_n(2,3)
         integral(3,3)=integral(3,3)+h*g_n(3,3)
    

         fluctcur_n(1,1)=fluctcur_n(1,1)+h*g_fluctcur_n(1,1)
         fluctcur_n(1,2)=fluctcur_n(1,2)+h*g_fluctcur_n(1,2)
         fluctcur_n(1,3)=fluctcur_n(1,3)+h*g_fluctcur_n(1,3)
         fluctcur_n(2,2)=fluctcur_n(2,2)+h*g_fluctcur_n(2,2)
         fluctcur_n(2,3)=fluctcur_n(2,3)+h*g_fluctcur_n(2,3)
         fluctcur_n(3,3)=fluctcur_n(3,3)+h*g_fluctcur_n(3,3)
      enddo 
      goto 20
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     end
c--------------------------------------------------------------------

 2    continue
c --------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     begin
c-------------------------------------------------------------------- 
cSm060306  
      do j=1,jmax
        p_perp=h*(j-0.5)

        call g_n_emis(p_perp,y,nll,np,theta,n,i_fkin,r,z,phi,
     .  g_n,g_fluctcur_n)

        integral(1,1)=integral(1,1)+g_n(1,1)
        integral(1,2)=integral(1,2)+g_n(1,2)
        integral(1,3)=integral(1,3)+g_n(1,3)
        integral(2,2)=integral(2,2)+g_n(2,2)
        integral(2,3)=integral(2,3)+g_n(2,3)
        integral(3,3)=integral(3,3)+g_n(3,3)

        fluctcur_n(1,1)=fluctcur_n(1,1)+g_fluctcur_n(1,1)
        fluctcur_n(1,2)=fluctcur_n(1,2)+g_fluctcur_n(1,2)
        fluctcur_n(1,3)=fluctcur_n(1,3)+g_fluctcur_n(1,3)
        fluctcur_n(2,2)=fluctcur_n(2,2)+g_fluctcur_n(2,2)
        fluctcur_n(2,3)=fluctcur_n(2,3)+g_fluctcur_n(2,3)
        fluctcur_n(3,3)=fluctcur_n(3,3)+g_fluctcur_n(3,3)
      enddo 

      integral(1,1)=integral(1,1)*h
      integral(1,2)=integral(1,2)*h
      integral(1,3)=integral(1,3)*h
      integral(2,2)=integral(2,2)*h
      integral(2,3)=integral(2,3)*h
      integral(3,3)=integral(3,3)*h  

      fluctcur_n(1,1)=fluctcur_n(1,1)*h
      fluctcur_n(1,2)=fluctcur_n(1,2)*h
      fluctcur_n(1,3)=fluctcur_n(1,3)*h
      fluctcur_n(2,2)=fluctcur_n(2,2)*h
      fluctcur_n(2,3)=fluctcur_n(2,3)*h
      fluctcur_n(3,3)=fluctcur_n(3,3)*h
      goto 20
c --------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     end
c--------------------------------------------------------------------

 3    continue
cSm060306
c --------------------------------------------------------------------
c     integration along the resonance curve using trapezoidal
c     formula
c     begin
c-------------------------------------------------------------------- 
      do j=0,jmax
         p_int=1.d0
         if((j.eq.1).or.(j.eq.jmax)) p_int=0.5d0
         p_perp=h*j
         if(p_perp.lt.1.d-3) p_perp=1.d-3
         if(j.eq.jmax) p_perp=p_perp-1.d-3

         call g_n_emis(p_perp,y,nll,np,theta,n,i_fkin,r,z,phi,
     .   g_n,g_fluctcur_n)
                 
         integral(1,1)=integral(1,1)+p_int*g_n(1,1)
         integral(1,2)=integral(1,2)+p_int*g_n(1,2)
         integral(1,3)=integral(1,3)+p_int*g_n(1,3)
         integral(2,2)=integral(2,2)+p_int*g_n(2,2)
         integral(2,3)=integral(2,3)+p_int*g_n(2,3)
         integral(3,3)=integral(3,3)+p_int*g_n(3,3)
    

         fluctcur_n(1,1)=fluctcur_n(1,1)+p_int*g_fluctcur_n(1,1)
         fluctcur_n(1,2)=fluctcur_n(1,2)+p_int*g_fluctcur_n(1,2)
         fluctcur_n(1,3)=fluctcur_n(1,3)+p_int*g_fluctcur_n(1,3)
         fluctcur_n(2,2)=fluctcur_n(2,2)+p_int*g_fluctcur_n(2,2)
         fluctcur_n(2,3)=fluctcur_n(2,3)+p_int*g_fluctcur_n(2,3)
         fluctcur_n(3,3)=fluctcur_n(3,3)+p_int*g_fluctcur_n(3,3)
      enddo 
 
      integral(1,1)=integral(1,1)*h
      integral(1,2)=integral(1,2)*h
      integral(1,3)=integral(1,3)*h
      integral(2,2)=integral(2,2)*h
      integral(2,3)=integral(2,3)*h
      integral(3,3)=integral(3,3)*h

      fluctcur_n(1,1)=fluctcur_n(1,1)*h
      fluctcur_n(1,2)=fluctcur_n(1,2)*h
      fluctcur_n(1,3)=fluctcur_n(1,3)*h
      fluctcur_n(2,2)=fluctcur_n(2,2)*h
      fluctcur_n(2,3)=fluctcur_n(2,3)*h
      fluctcur_n(3,3)=fluctcur_n(3,3)*h

      goto20
c --------------------------------------------------------------------
c     end integration along the resonance curve using trapezoidal
c     formula
c     end
c--------------------------------------------------------------------

 4    continue
cSm070417
c -------------------------------------------------------------------
c     integration along the resonance curve using 
c     the adaptive Simpson function
c--------------------------------------------------------------------

c      write(*,*)'emission.f intgr_emis before adaptive_simpson'

      call  calc_emission_integral_array_by_adaptive_simpson
     &(y,nll,np,theta,n,i_fkin,r,z,phi,p_permax,epsi,
     &integral,fluctcur_n)
c      write(*,*)'intgr_emis after adaptive_simpson'
c      write(*,*)'integral',integral
c      write(*,*)'fluctcur_n',fluctcur_n
c--------------------------------------------------------------------
 20   continue

      integral(2,1)=-integral(2,1)
      integral(3,1)= integral(3,1)
      integral(3,2)=-integral(3,2)
   
      fluctcur_n(2,1)=-fluctcur_n(2,1)
      fluctcur_n(3,1)= fluctcur_n(3,1)
      fluctcur_n(3,2)=-fluctcur_n(3,2)
  
 10   continue
      return
      end   
  

      subroutine g_n_emis(p_perp,y,nll,np,theta,n,i_fkin,r,z,phi,
     .g_n,g_fluctcur_n)
c     calculates under integral complex marix functions
c     g_n(3,3)=sum{k}g_nk and g_fluctcur_n(3,3)=sum{k}g_fluctcur_nk

c     input
c       p_perp  - momentum divided by (mc)
c       y       = omega_ce/omega for the electron rest mass    
c       nll     - N_parallel
c       np      - N_perpendicular 
c       theta   - mc**2/T
c       n       - n number of the given resonance harmonic
c       i_fkin =0 the usage of the analytical relativistic distributin
c              =1 the usage of the numerical 3D distribution from diskf 
c                 written be CQL3D code
c       r       - the major radius (it is used for i_fkin=1)
c       z       - the vertical coordinate  (it is used for i_fkin=1)
c       phi     - toroidal angle(it is used for i_fkin=1)
c     output
c       g_n(3,3) under integral matrix double complex function
c       g_fluctcur_n(3,3)
c-------------------------------------------------------
      implicit none
c-----input        
      double precision p_perp,y,nll,np,theta,r,z,phi
      integer n,i_fkin
c-----external root_res, s_calcn, zeroK, unde_intgr
      double precision u_n,u_flcur_n    
c-----output
      double complex g_n(3,3), g_fluctcur_n(3,3)
      
c-----local
      double complex sn_k(3,3),g_nk(3,3),g_fluctcur_nk(3,3)
      double precision gamma, p_par_rl(2),coeff,eps,ps,ps_max,pi,
     +dgdp_par,u_nk,u_flcur_nk,coeff_fc
      integer k,kmax,i,j
c     kmax - total number of the resonace condition root p_perp_k(p_perp)
      
      pi=4.d0*datan(1.d0)
c-----calculations ot the roots p_par_rl of the resonace condition
c-----gamma=N_par*p_par_rl+nY

      call root_res(p_perp,nll,n,Y,kmax,p_par_rl)
c      write(*,*)'g_n_emis number of roots p_par_rl(p_perp) kmax',kmax
c-----initialize g_n
      call zeroK(g_n)
      call zeroK(g_fluctcur_n)

      if (kmax.eq.0) goto 20
  
c      eps=1.d-9 ! the min value of Maxwell exponent
c      ps_max=(1.d0-dlog(eps)/theta)**2-1.d0

      
      do k=1,kmax  
         ps=p_perp*p_perp+p_par_rl(k)*p_par_rl(k)
         gamma=dsqrt(1.d0+ps)

c         if(ps.gt.ps_max) then
c           write(*,*)'emission.f g_n_emis k,ps,ps_max',k,ps,ps_max
c           goto 10
c         endif

c         write(*,*)'g_n_emis k,p_perp,p_par_rl(k),ps',
c     .   k,p_perp,p_par_rl(k),ps

         call s_calcn(n,p_perp,p_par_rl(k),np,y,sn_k)
c         write(*,*)'g_n_emis k,sn_k',k,sn_k

c--------resonance condition uses the delta function with argument
c        g_delta=gamma-nll*p_par-n*y
c        the derivative from this argument d(g_delta)/dp_par=dgdp_par

         dgdp_par=dabs((p_par_rl(k)-nll*gamma)/gamma)

cSm030415
c         if(dgdp_par.lt.1.d-3)then
c             write(*,*)'dgdp_par k,p_par_rl(k),nll,p_perp,gamma',
c     +        k,p_par_rl(k),nll,p_perp,gamma
c         endif
    

         call  unde_intgr(p_perp,p_par_rl(k),y,nll,theta,n,i_fkin,
     +   r,z,phi,u_nk,u_flcur_nk)

c         write(*,*)'in g_n_emis after unde_integr u_nk,u_flcur_nk',
c     +   u_nk,u_flcur_nk

         coeff=2.d0*pi*u_nk/dgdp_par
         coeff_fc=2.d0*pi*u_flcur_nk/dgdp_par
c         write(*,*)'in g_n_emis u_nk,u_flcur_nk',u_nk,u_flcur_nk
c         write(*,*)'in g_n_emis coeff,coeff_fc',coeff,coeff_fc

         g_nk(1,1)=coeff*sn_k(1,1)  
         g_nk(1,2)=coeff*sn_k(1,2)
         g_nk(1,3)=coeff*sn_k(1,3)
         g_nk(2,2)=coeff*sn_k(2,2)
         g_nk(2,3)=coeff*sn_k(2,3)
         g_nk(3,3)=coeff*sn_k(3,3)

         g_nk(2,1)=-g_nk(1,2)
         g_nk(3,1)= g_nk(1,3)
         g_nk(3,2)=-g_nk(2,3)      


         g_fluctcur_nk(1,1)=coeff_fc*sn_k(1,1)  
         g_fluctcur_nk(1,2)=coeff_fc*sn_k(1,2)
         g_fluctcur_nk(1,3)=coeff_fc*sn_k(1,3)
         g_fluctcur_nk(2,2)=coeff_fc*sn_k(2,2)
         g_fluctcur_nk(2,3)=coeff_fc*sn_k(2,3)
         g_fluctcur_nk(3,3)=coeff_fc*sn_k(3,3)


         do i=1,3
            do j=1,3
               g_n(i,j)=g_n(i,j)+g_nk(i,j)
               g_fluctcur_n(i,j)=g_fluctcur_n(i,j)+g_fluctcur_nk(i,j)
            enddo
         enddo
  
 10      continue
      enddo !kmax

 20   continue

      return
      end

      subroutine unde_intgr(p_perp,p_par,y,nll,
     +theta,n,i_fkin,r,z,phi,u_n,u_flcur_n)
    

c     calculates under integral functions
c     u_n=U_n(np=p_per,p_perp=p_parallel)=(1/gamma)*
c     (n*Y*df/dp^_perpendicular+N_parallel*p^_perpendicular*df/dp^_parallel)
c     u_n is used in under integral complex marix function G_nk(3,3)
c
c     u_flcur_n=f*p^perp/gamma
c     u_flcur_n is used in under integral complex matrix fluct_car_nk(3,3)
c     Here p^=p/mc
c     input
c       p_perp  - perpendicular momentum divided by (mc)
c       p_par     parallel momentum/mc
c       y       = omega_ce/omega for the electron rest mass    
c       nll     - N_parallel 
c       theta   - mc**2/T
c       n       - n number of the given resonance harmonic
c       i_fkin =0 the usage of the analytical relativistic Maxwellian distributin
c              =1 the usage of the numerical 3D distribution from diskf or
c                 netcdfnm.nc files or the mesh distribution given analytically
c                 written be CQL3D code
c       r       - the major radius (it is used for i_fkin=1)
c       z       - the vertical coordinate  (it is used for i_fkin=1)
c       phi     - toroidal angle
c-------------------------------------------------------
      implicit none
c     input        
      double precision p_perp,p_par,y,nll,theta
      double precision r,z,phi
      integer n,i_fkin

c     external besk2as, root_res,psif,rhopsi,bmin_psi,b
      double precision psif,rhopsi,bmin_psi,b,fdens_fdist  

c     output
      double precision u_n,u_flcur_n

c_for test
      double precision vnormloc,massloc
      COMMON /dskin1/vnormloc,massloc  !v (cm/sec),m( g)
cendtest

c     local
      
      double precision gamma, bk2,f_maxw,ps,p,ps_max,eps,pi,d_maxw_dp,
     + m_e,clight,energy,psi,rho,pitch,pitch0,bmin,btotal,
     + c_pitch,s_pitch,s_pitch0,c_pitch0,
     + dptc0dpt,                          !d(pitch0)/d(pitch)
     + dens,u_np,u_flcur_np
     
c     distributin function and its derivatives
      double precision fdist0,dfdx,dfdpitch,dfdpitc0,dfdp
      integer initial     
c-----!
      pi=4*datan(1.d0)
      ps=p_perp*p_perp+p_par*p_par
      p=dsqrt(ps)
      gamma=dsqrt(1.d0+ps)

      if (i_fkin.eq.0) then
c--------usage of the analytical relativistic Maxwellian distribution

c        calculation the Mackdonalds function bk2=K_2(theta)*EXP(theta)	
         call besk2as(theta,bk2)
  
         eps=1.d-9 ! the min value of Maxwell exponent
         ps_max=(1.d0-dlog(eps)/theta)**2-1.d0
      
c         if(ps.gt.ps_max) then
c           u_n=0.d0
c           goto 10
c         endif
           
          f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)
          u_n=-p_perp*theta*f_maxw*(n*y+nll*p_par)/(gamma*gamma)

          u_flcur_n=p_perp*f_maxw/gamma
cSm070126  
c          write(*,*)' emission unde_int'
c          write(*,*)' old u_n, old f_maxw,u_n,u_flcur_n',
c     &     f_maxw,u_n,u_flcur_n
          call analytical_distrib(theta,p,r,z,phi,
     &    f_maxw,d_maxw_dp)

          u_n=p_perp*d_maxw_dp*(n*y+nll*p_par)/gamma/p

          u_flcur_n=p_perp*f_maxw/gamma
c          write(*,*)'unde_intgr new f_maxw,u_n u_flcur_n',
c     &     f_maxw,u_n,u_flcur_n
         goto10 
      endif

      if (i_fkin.eq.1)then
c--------usage of the numerical 3D relativistic distribution from diskf file
c        written by CQL3D code
         initial=0
c        energy=p**2/(2m) ps=(p/mc)**2. energy should be in KeV
c        pitch
c        rho
         m_e=9.1094d-28                !g
         clight=2.99792458d10          !cm/sec
         energy=0.5d0*ps*m_e*clight**2 !erg 
         energy=energy/1.6022d-9         !KeV
 
         psi=psif(z,r)                 !poloidal flux
         rho=rhopsi(psi)!small radius

         if (p.ne.0.d0) then
             pitch=dacos(p_par/p)
         else
             pitch=0.d0
         endif
         if (pitch.gt.pi) pitch=pi
         if (pitch.lt.-pi) pitch=-pi
          

c------- pitch0 is pitch angle for the point with the minimal
c------- value of the magnetic field at the give flux surface with the poloidal flux=psi  
         pitch0=pitch
         dptc0dpt=1.d0
cSm060227
c         goto 100

         bmin=bmin_psi(psi)
         btotal=b(z,r,phi)
         
         if (p.ne.0.d0) then
             c_pitch=p_par/p     !cos(pitch)
             s_pitch=p_perp/p    !sin(pitch)
         else
             c_pitch=1.d0     !cos(pitch)
             s_pitch=0.d0     !sin(pitch)
         endif
         
c--------s_pitch0=sin and c_pitch0=cos of the pitch angle for the minimal b
         s_pitch0=s_pitch*dsqrt(bmin/btotal) !sin(pitch0) 
c         write(*,*)'bmin,btotal,s_pitch,s_pitch0',
c     +   bmin,btotal,s_pitch,s_pitch0

         if (s_pitch0.gt.1.d0) then
            s_pitch0=1.d0
            c_pitch0=0.d0
         else    
            c_pitch0=dsqrt(1.d0-s_pitch0**2) !dabs(cos(pitch0))
         endif
   
         if (c_pitch.lt.0.d0) then
            c_pitch0=-c_pitch0
         endif

         pitch0=dacos(c_pitch0)
         if (pitch0.gt.pi) pitch0=pi
         if (pitch0.lt.-pi) pitch0=-pi
         
c         write(*,*)'emission in u_n pitch,pitch0',pitch,pitch0

c--------d(pitch0)/d(pitch)
c         dptc0dpt=(c_pitch/c_pitch0)*(s_pitch0/s_pitch)

cSm051206
         if((dabs(c_pitch0).lt.1.d-12).or.(dabs(s_pitch).lt.1.d-12))
     &   then
            dptc0dpt=1.d0
         else
            dptc0dpt=(c_pitch/c_pitch0)*(s_pitch0/s_pitch)
         endif

c        write(*,*)'c_pitch,c_pitch0',c_pitch,c_pitch0
c         write(*,*)'s_pitch,s_pitch0',s_pitch,s_pitch0

 100     continue

c--------calculations of the distribution function fdist0 
c        and its derivatives df/dx,df/dpitch0
c        x = momentum-per-mass(nomalized to maximum 1.)
c        vel=dsqrt(2.d0*energy*1.d3*1.6022d-12/fmass(1))
c        at the point (rho,energy,pitch0). Here pitch0=pitch0(r,z,pitch)

c         write(*,*)'unde_intgr initial,energy,pitch0,rho',
c     +                         initial,energy,pitch0,rho
         call dskin(initial,energy,pitch0,rho,fdist0,dfdx,dfdpitc0,
     +              dfdp,2)

c         write(*,*)'unde_intgr after dskin fdist0,dfdx,dfdpitc0,dfdp',
c     +                                     fdist0,dfdx,dfdpitc0,dfdp

         dfdpitch=dfdpitc0*dptc0dpt
c         write(*,*)'dfdpitc0,dptc0dpt,dfdp,dfdx,fdist0',
c     +   dfdpitc0,dptc0dpt,dfdp,dfdx,fdist0

c         u_n=1.d0/gamma*(
c     +   n*y*(dfdp*p_perp/p-dfdpitch*p_par/p**2)+
c     +   nll*p_perp*(dfdp*p_par/p+dfdpitch*p_perp/p**2))*
c     +   (clight/vnormloc)**3

         u_n=1.d0/gamma*(
     +   n*y*(dfdp*p_perp/p+dfdpitch*p_par/p**2)+
     +   nll*p_perp*(dfdp*p_par/p-dfdpitch*p_perp/p**2))
     +   *(clight/vnormloc)**3

c         write(*,*)'unde_intgr u_n,gamma,p,vnormloc,n,y',
c     +   u_n,gamma,p,vnormloc,n,y
c         write(*,*)'dfdp,p_perp,dfdpitch,p_par',
c     +   dfdp,p_perp,dfdpitch,p_par
c         write(*,*)'nll,clight',nll,clight

         u_flcur_n=(p_perp*fdist0/gamma)*(clight/vnormloc)**3
c         u_flcur_n=(p_perp*fdist0/gamma)

c         write(*,*)'u_n,u_flcur_n',u_n,u_flcur_n                       

ctest_begin(comparison of fdist0 and the analytical maxwellian diatribution and its derivatives
c--------calculation the Mackdonalds function bk2=K_2(theta)*EXP(theta)	
c          call besk2as(theta,bk2)
c          f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)
c          f_maxw=f_maxw*(vnormloc/clight)**3
c          write(*,*)'in u_n_emis vnormloc,theta,gamma,f_maxw,fdist0',
c     +    vnormloc,theta,gamma,f_maxw,fdist0
c          u_np=-p_perp*theta*f_maxw*(n*y+nll*p_par)/(gamma*gamma)
c          write(*,*)'unde_intgr u_np,u_n',u_np,u_n
c          u_flcur_np=(p_perp*fdist0/gamma)
c          write(*,*)'unde_intgr u_flcur_np u_flcur_n',
c     +    u_flcur_np,u_flcur_n
c          u_n=u_np
c          u_flcur_n=u_flcur_np
cend_test
      endif !i_fkin.eq.1

 10   continue

      return
      end    

      subroutine calc_emis_coef(xe,ye,T_kev,cnpar,cnper,cnray,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +i_resonance_curve_integration_method,epsi,
     +i_fkin,r,z,phi,
     +cex,cey,cez,cflown,vgroup,frqncy,
     +iabsorp_collisional,coll_mult,tempe,dense,zeff,
     +al_emis,j_emis,temprad)
c---- calcultes the emission absorptivity al_emis (1/cm) 
c     emissivity j_emis (erg*sec/cm**3)  and
c     radiation  temperature temprad

      implicit none

c-----input
      double precision xe,ye,T_kev,cnpar,cnper,cnray,z,r,phi,epsi
c     cnray  is the ray refractive index N_ray
      integer n_relt_harm1, n_relt_harm2,n_relt_intgr,i_fkin,
     &i_resonance_curve_integration_method
c----------------------------------------------------------------------
c     i_resonance_curve_interation_method
c
c     i_resonance_curve_integration_method=1 angle integration
c     now it works for ellipse case only                             c
c
c     i_resonance_curve_integration_method=2 rectangle formula
c     for p_perp integration
c     
c     i_resonance_curve_integration_method=3 trapezoidal formula
c     for p_perp integration
c
c     i_resonance_curve_integration_method=4 !adaptive Simpson 
c     for p_perp integration
c
c     epsi is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------ 
      double complex cex,cey,cez
c     wave energy flux for the normalized field |E|=1
c     cflown=(vectorB*vector_B^*+vector_E*d(omega*tensor_eps_herm)/dw*vector_E)
      double complex cflown  
      double precision vgroup  !in clight
      double precision frqncy   !GhZ
c     the data for collisional absorption
      integer iabsorp_collisional ! =0 do not calculate collisional absorption
                                  ! =1 to calculte collisonal absorption
      double precision coll_mult  ! multiplier for coll absorption
      double precision tempe,     !electron temperature  (keV)  
     & dense,                     !electron density   (10*13/cm**3	 
     & zeff                       !plasma charge
c-----external
      double precision rrind     
c-----output  
      double precision al_emis,j_emis,temprad
c-----local        
      double complex aK(3,3),fluctcur(3,3),ce(3)
      double precision v_gr_cm_sec !cm/sec
      double precision omega       !2*pi*f, here f in [HZ]
      double precision omegpedce,omegdce!omega_pe/omega_ce and omega/omega_ce
      double precision cn          !refructive index
      double precision clight,pi      

      write(*,*)'calc_emis_coef'
      write(*,*)'i_resonance_curve_integration_method,epsi',
     &i_resonance_curve_integration_method,epsi

      call emis_tens(xe,ye,T_kev,cnpar,cnper,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +i_resonance_curve_integration_method,epsi,
     &i_fkin,r,z,phi,
     +aK,fluctcur)

      ce(1)=cex
      ce(2)=cey
      ce(3)=cez
      clight=2.99792458d10          !cm/sec
      v_gr_cm_sec=clight*vgroup     !cm/sec
      pi=4.d0*datan(1.d0)
      omega=2.d0*pi*frqncy*1.d9     !1/esc
           
c      write(*,*)'calc_coef omega,frqncy,cn,cnray,cflown',
c     +omega,frqncy,cn,cnray,cflown
     
c      write(*,*)'in calc_emis_coef aK',aK
c      write(*,*)'in calc_emis_coef fluctcur',fluctcur
c      write(*,*)'in calc_emis_coef cflown',cflown
c      write(*,*)'in calc_emis_coef ce',ce

      call emiss_coef(aK,fluctcur,ce,cflown,v_gr_cm_sec,omega,cnray,
     +iabsorp_collisional,coll_mult,tempe,dense,zeff,
     +al_emis,j_emis)  
     
c      write(*,*)'in calc_emis_coef cnray,al_emis,j_emis',
c     &cnray,al_emis,j_emis

c-----radiation temperature      
      if(al_emis.ne.0.d0)  then
         if(al_emis.gt.1.d-12) then
c----------we choosed this number 1.d-12 to avoid the division
c          the small number for the other small number.
           temprad=(2.d0*pi)**3*(clight)**2/(omega*cnray)**2
     +             *j_emis/al_emis/1.6022d-9
         else
           temprad=0.d0
         endif
c         write(*,*)'in calc_emis_coef T_kev,temprad=',T_kev,temprad
      else
         temprad=0.d0
      endif

      return
      end 


      subroutine emis_bin_sn(n,al_emis_n,al_emis_n1,j_emis_n,j_emis_n1,
     +nr_n,nr_n1,s_n,s_n1,in_sn)  
c-----calculates the emission at the bin boundary s=s_n from one bin s_n<s<s_n+1
      implicit none
c-----input
      integer n ! the number of ray bin, s_n(at n=1)=0,It should be: 1=<n<nrelt-1 
      double precision al_em is_n,j_emis_n,    !absorption and emission at s=s_n
     +                 al_emis_n1,j_emis_n1,  !absorption and emission at s=s_(n+1)
     +                 nr_n,nr_n1,            !N_ray rayrefractive index at s=s_n and s_(n+1)
     +                 s_n,s_n1               !the lengths along the ray (cm) s_n1=s_n+1
c     al_emis [1/cm]
c     j_emis  [egr/(cm**3*sec)]
c-----output
      double precision in_sn !emission I_n at the detector side of nth bin at s=s_n
      
c-----local
      double precision j_emis_np, !j_emis(s_n+0.5)
     +al_emis_np,                 !al_emis(s_n+0.5)
     +nrs_np,                     !N_r**2(s_n+0.5)
     +p,p1

        
c-----calculation of I_n=I_n(s=s_n) at the detector side of nth bin at s=s_n
      j_emis_np=0.5d0*(j_emis_n+j_emis_n1)
      al_emis_np=0.5d0*(al_emis_n+al_emis_n1)
      nrs_np=(0.5d0*(nr_n+nr_n1))**2

      p=nrs_np*al_emis_np
      p1=j_emis_np*nr_n**2
      
      if ((p1.eq.0.d0).or.(p.eq.0.d0)) then
         in_sn=0.d0
      else
         in_sn=p1/p*(1.d0-dexp(-al_emis_np*(s_n1-s_n)))
cSm011227-to conform with Bob
c          in_sn=j_emis_n*(1.d0-dexp(-al_emis_n*(s_n1-s_n)))
      endif

      return
      end  

      subroutine emis_bin_s1(n,nrayelt,wsn,wal_emis,wnray,in_sn,
     +tau_n,in_0,nrayelt_o_cutoff,transm_ox)  
c-----calculates the emission at the plasma boundary s=s_1
c     from one nth bin s_n<s<s_n+1
      implicit none
c-----input
      integer n, !is the number of ray bin, s_n(at n=1)=0,
                 !It should be: 1=<n<nrelt-1 
     +nrayelt,   ! the total number of the output points along the ray
     &nrayelt_o_cutoff !the number of the output point 
                       !wsn(nrayelt_o_cutoff)=s_O_mode 
                       !along the ray (at O-mode part)
                       !where the OX transmission
                       !coefficient transm_ox has the maximal value
      double precision wal_emis(*),    !absorption at s=s_n
     +wsn(*),  !the lengths s=sn  along the ray (cm) 
     +wnray(*),!N_ray ray refractive index along the ray 
     +in_sn,   !emission I_n at the detector side of nth bin at s=s_n
     *transm_ox !transmission coefficient for OX conversion
c     wal_emis [1/cm]

c-----output
      double precision in_0 !emission I_n at the plasma boundary from nth bin at s=s_1=0
      double precision tau_n                  !tau_n from n_th bin
c-----local
      double precision 
     +al_emis_np                 !al_emis(s_n+0.5)                  
      integer m

      if (n.lt.1) then
         write(*,*)'in emis_bin_s1 n<1, It should be 1=< n =but n=',n
         stop
      endif

      if (n.gt.(nrayelt-1)) then
         write(*,*)'in emis_bin_s1 n.gt.nrayelt-1.
     +   It should be n=<nrayelt-1'
         write(*,*)'in emis_bin_s1 n=',n,'nrayelt=',nrayelt
         stop
      endif

c-----tau_n=integral{0,s_n}(al_emis(s)ds)       

      tau_n=0.d0
      if (n.gt.1) then
         do m=2,n
            if(m.ne.nrayelt_o_cutoff+1) then
              al_emis_np=0.5d0*(wal_emis(m)+wal_emis(m-1))
c              tau_n=tau_n+(wsn(m)-wsn(m-1))*al_emis_np
cSm011227-to conform with Bob
             tau_n=tau_n+(wsn(m)-wsn(m-1))*wal_emis(m)
            endif

         enddo
      endif

c-----calculation of I_n=I_n(s=s(1)=0), at the plasma boundary
cSm070120
c      in_0=(wnray(1)/wnray(n))**2*in_sn*dexp(-tau_n)
      in_0=(1.d0/wnray(n))**2*in_sn*dexp(-tau_n)

c      write(*,*)'emission.f in emis_bin_s1 n,wnray(1),wnray(n),in_sn',
c     &n,wnray(1),wnray(n),in_sn
c      write(*,*)'tau_n,in_0',tau_n,in_0

      if ((n.ge.nrayelt_o_cutoff+1).and.(nrayelt_o_cutoff.gt.1)) then
         in_0=in_0*transm_ox
      endif  

c      write(*,*)'nrayelt_o_cutoff,transm_ox,in_0',
c     &nrayelt_o_cutoff,transm_ox,in_0 
    
      return
      end  



      subroutine emission(tol_emis,nrayelt,nrelta,wsn,wal_emis,wj_emis,
     &wnray,in_sn,in_0,i_0,i_0sn,wtaun_em,nrayelt_emis,
     &nrayelt_o_cutoff,transm_ox,nrayelt_o_cutoff_emis)  
c-----calculates the emission at the plasma boundary s=s_1 from 
c     one the ray s(1) < s < s(nrelt)
      implicit none
c-----input 
      integer nrayelt        ! the total number of the output points
                             ! along the ray
      integer nrelta  !the max number of the output points 
      integer nrayelt_o_cutoff !the number of the output point 
                               !wsn(nrayelt_o_cutoff)=s_O_mode 
                               !along the ray (at O-mode part)
                               !where the OX transmission
                               !coefficient transm_ox has the maximal value
      double precision
     +tol_emis,    ! tolerance parameter to add the mesh point s_n 
                   ! if in_0(n) > tol_emis*i_0  
     +wal_emis(*),wj_emis(*),!absorption and emission at s=s_n
     +wsn(*),    !the length s=s_n along the ray (cm) 
     +wnray(*),  !N_ray along the ray
     &transm_ox  !transmission coefficient for OX conversion
     
c     wal_emis [1/cm]
c     wj_emis  [egr/(cm**3*sec)]

c-----output
      double precision
     +in_sn(*), ! emission I_n at the detector side of nth bin at s=s_n
     +in_0(*),  ! emission at the plasma boundary s=s_1 from one
                ! nth bin s_n<s<s_n+1
     +i_0, !emission I_0 at the plasma boundary from the ray at s=s(1)=0
     +i_0sn(*), ! sum{k=1,n}[in_0(k)]! emission at the plasma boundary
                ! from the parts 0<s<sn of the ray
     +wtaun_em(*)! tau_n from n-th bin for plotting 
      integer nrayelt_emis ! number of the emission output points after
                           ! the addition of new points
      integer nrayelt_o_cutoff_emis !the number of the output point
                                    !(after addition new points)  
                                    !wsn(nrayelt_o_cutoff_emis)=s_O_mode 
                                    !along the ray (at O-mode part)
                                    !where the OX transmission
                                    !coefficient transm_ox has the maximal value
c-----local
      integer n,n_midway,isplit,
     +n_checked,nrayelt_emis_new,k,i_ox_jump
        
      nrayelt_emis=nrayelt
      nrayelt_o_cutoff_emis=nrayelt_o_cutoff

cSAP080422
      if (nrayelt.lt.2) then
         write(*,*)'in emission.f subroutine emission'
         write(*,*)'nrayelt=',nrayelt
         write(*,*)'nrayelt.lt.2, The code will not calculate emission'
         return 
      endif

      if (nrayelt.gt.nrelta) then
         write(*,*)'in emission nrayelt.gt.nrelta nrayelt,nrelta',
     .   nrayelt,nrelta
         write(*,*)'it should be nrayelt=<nrelta'
         write(*,*)'increas nrelta in param.i'
         stop
      endif  

c------------------------------------------------------------------
c     initialization 
      call bcast(in_0,0.d0,nrelta)
c-------------------------------------------------------------------

      if((nrayelt_o_cutoff.gt.1).and.(nrayelt_o_cutoff.lt.nrayelt))then
        i_ox_jump=1 !to use OX conversion jump
      else
        i_ox_jump=0 !do not use OX conversion jump
      endif
   
      write(*,*)'emission i_ox_jump',i_ox_jump 

      do n=1,nrayelt-1
c         write(*,*)'emssion n=',n
         if((i_ox_jump.eq.1).and.(n.eq.nrayelt_o_cutoff)) then
           in_sn(n)=0.d0 !in OX jump bin area
           in_0(n)=0.d0  !from OX jump bin area
         else  
c----------calculate in_sn emission I_n at the detector side
c          of n th bin at s=s_n
c          write(*,*)'in emission before emis_bin_sn' 

           call emis_bin_sn(n,wal_emis(n),wal_emis(n+1),wj_emis(n),
     +     wj_emis(n+1),wnray(n),wnray(n+1),wsn(n),wsn(n+1),in_sn(n))
 
c        write(*,*)'in emission after emis_bin_sn n,in_sn(n)',n,in_sn(n)  

c----------calculate in_0 the emission at the plasma boundary s=s_1
c          from one nth bin s_n<s<s_n+1
           call emis_bin_s1(n,nrayelt,wsn,wal_emis,wnray,in_sn(n),
     +     wtaun_em(n),in_0(n),nrayelt_o_cutoff,transm_ox)      
 
c           write(*,*)'in emission after emis_bin_s1 n,in_0(n),i_0sn(n)',
c     &     n,in_0(n),i_0sn(n)
        endif
      enddo
     
c      write(*,*)'emssion before add point nrayelt ',nrayelt
         
c      do n=1,nrayelt
c         write(*,*)'n,wsn(n)',n,wsn(n)
c      enddo


      n_checked=1
        
 10   continue
c-----calculate i_0 the total emission at the plasma boundary from one ray 
      i_0=0.d0
c      do n=1,nrayelt_emis
      do n=nrayelt_emis,1,-1
         i_0=i_0+in_0(n)
         i_0sn(n)=i_0
c         write(*,*)'n,in_0(n)',n,in_0(n)
      enddo  
         
c      write(*,*)'before add point i_0',i_0     

ctest
c      do n=1,nrayelt_emis
c         write(*,*)'emission n,wal_emis(n),wj_emis(n)',
c     .   n,wal_emis(n),wj_emis(n)
c      enddo
cendtest

cSm070110
      if (nrayelt_emis .lt.4) then
c-------------------------------------------------------------
c         Too small ray elements for spline approximation
c         The code will not check condition in_0(n)>tol_emis*i_0 
c         and will not add points along the ray
c------------------------------------------------------------------
         goto 30
      endif

c-----check the condition in_0(n)>tol_emis*i_0 

      do n=1,nrayelt_emis-1
c         write(*,*)'emission n=',n,'in_0(n),i_0,tol_emis*i_0',
c     &                              in_0(n),i_0,tol_emis*i_0
         if (in_0(n).gt.tol_emis*i_0) then
c-----------add the new midway point
            isplit=1
            nrayelt_emis_new=nrayelt_emis+1

            if((i_ox_jump.eq.1).and.(n.ge.nrayelt_o_cutoff_emis+1)) then
              nrayelt_o_cutoff_emis=nrayelt_o_cutoff_emis+1          
            endif
           
            if (nrayelt_emis_new.gt.nrelta) then
                write(*,*)'in emission nrayelt_emis_new.gt.nrelta 
     .          nrayelt_emis_new,nrelta', nrayelt_emis_new,nrelta
                write(*,*)'it should be nrayelt_emis_new=<nrelta'
                write(*,*)'increas nrelta in param.i'
                stop
            endif
  
            n_midway=n+1
c            write(*,*)'emission before add_ray_point n',n
            call add_ray_point(n_midway,nrayelt_emis)! It writes nrayelt_emis
                                                     ! in write.i
c            write(*,*)'emission after add_ray_point n',n

            do k=n,n+1
c-----------------------------------------------------------------------
c              calculate the emission In_sn(n) at the boundaries
c              of two new bins 
c              (s_n,s_n_1=n_midway) and (s_n+1=n_midway,s_n+2=n_midway+1)
c------------------------------------------------------------------------
               call emis_bin_sn(k,wal_emis(k),wal_emis(k+1),wj_emis(k),
     +         wj_emis(k+1),wnray(k),wnray(k+1),wsn(k),wsn(k+1),
     +         in_sn(k))
c               write(*,*)'after emis_bin_sn k',k
c               write(*,*) 'wal_emis(k),wal_emis(k+1)',
c     &         wal_emis(k),wal_emis(k+1)
c               write(*,*)'wj_emis(k),wj_emis(k+1)',
c     &         wj_emis(k),wj_emis(k+1)
c               write(*,*)'wnray(k),wnray(k+1)',wnray(k),wnray(k+1)
c               write(*,*)'wsn(k),wsn(k+1),in_sn(k)',
c     &         wsn(k),wsn(k+1),in_sn(k)
            enddo !k

            do k=n,nrayelt_emis_new-1
c----------------------------------------------------------------------
c              calculate in_0 the emission at the plasma boundary s=s_1
c              from one nth bin s_n<s<s_n+1
c----------------------------------------------------------------------
               call emis_bin_s1(k,nrayelt_emis_new,wsn,wal_emis,wnray,
     +                          in_sn(k),wtaun_em(k),in_0(k),
     +                          nrayelt_o_cutoff_emis,transm_ox)
c               write(*,*)'after add n,in_sn(k),wtaun_em(k),in_0(k)',
c     &         n,in_sn(k),wtaun_em(k),in_0(k)          
            enddo  !k
            goto 20
         else
            isplit=0
         endif
      enddo   

 20   continue
      if (isplit.eq.1) then
          n_checked=n_midway
          nrayelt_emis=nrayelt_emis_new
          goto10
      endif
ctest
c      do n=1,nrayelt_emis
c         write(*,*)'end emission n,wal_emis(n),wj_emis(n),in_0(n)',
c     &    n,wal_emis(n),wj_emis(n),in_0(n)
c      enddo
cendtest

 30   continue

      return
      end

      real*8 function rrind_new_100923(npar,nper,omega,omegpe)
c     input npar,nper,omega=omega/omega_ce,omegape=omegape/omegac_e
c     In this new version (100923) the multipiler st was removed
c     from the ratio
c     rrind=n2*st*droot*ddenom
c     Here ddenom =st*ddenom.
c     Now this version uses
c     rri2=n2*droot/ddenom1

      implicit none
c-----input 
      real*8 npar,nper,omega,omegpe 
c-----locals
      real*8 npar2,nper2,n2,eps1,eps2,eps3,
     &rn,st,ct,st2,ct2,var1,var2,var3,eps22,
     &dnumer,ddenom,dndt,droot,deriv,rri2,
     &ddenom1,dndt1,dnumer1
c---- external for terst only
      real*8 rrind
c
c     calculate ray refractive index for cold electron plasma
c

      eps1=1.d0-omegpe*omegpe/(omega*omega-1.d0)
      eps2=-omegpe*omegpe/(omega*(omega*omega-1.d0))
      eps3=1.d0-omegpe*omegpe/(omega*omega)
c
      npar2=npar*npar
c
      nper2=nper*nper
      n2=nper2+npar2

c
c     calculate ray refractive index squared
c
      rn=dsqrt(n2)
      st=nper/rn
      ct=npar/rn
      st2=st*st
      ct2=ct*ct
      var1=eps1-n2
      var2=eps3-nper2
      var3=eps1-npar2
      eps22=eps2*eps2

c      dnumer=-st*ct*var1*var2+st*ct*var3*var1-eps22*st*ct
c     1       +n2*var1*st*ct*(ct2-st2)

c      write(*,*)'rrind new 1 dnumer', dnumer

cSAP100921
      dnumer1=(-ct*var1*var2+ct*var3*var1-eps22*ct
     1       +n2*var1*ct*(ct2-st2))
      dnumer=st*dnumer1

c      write(*,*)'rrind new  dnumer1', dnumer1
c      write(*,*)'rrind new 2 dnumer', dnumer

      ddenom=eps22*st2+n2*n2*st2*ct2-var1*var2*ct2-var3*var2
     1       -st2*var3*var1-2.0*n2*st2*ct2*var1
       
c      write(*,*)'rrind new ddenom',ddenom
c      dndt=rn*dnumer/ddenom
c      write(*,*)'rrind new 1 dndt',dndt
cSAP100921
      dndt1=rn*dnumer1/ddenom
      dndt=st*dndt1

c      write(*,*)'rrind new 2 dndt',dndt
c      write(*,*)'rrind new dndt1',dndt1

      droot=dsqrt(1.d0+(dndt/rn)*(dndt/rn))

c      write(*,*)'rrind new droot',droot
c
      deriv=( var1*var2*(st2-ct2)+2.d0*rn*dndt*var2*st*ct
     1       +2.d0*rn*dndt*st2*st*ct*var1+2.0*st2*ct2*var1*n2
     2       +(ct2-st2)*var3*var1-
     1       2.d0*rn*dndt*st*ct*(ct2*var1+var3)
     3       +2.d0*st2*ct2*n2*var1
     4       +eps22*(st2-ct2)
     5       +2.d0*rn*dndt*st*ct*(ct2-st2)*(var1-n2)
     6       +n2*var1*(ct2*ct2-3.d0*ct2*st2+st2*st2-
     1       3.d0*st2*ct2) )/ddenom

c      write(*,*)'rrind new 1 deriv',deriv

      deriv=deriv - (dnumer/(ddenom*ddenom))*(2.d0*eps22*st*ct
     1      +n2*n2*2.d0*st*ct*(ct2-st2)+4.d0*n2*rn*dndt*ct2*st2
     2      +2.d0*ct*st*var1*var2+
     1      2.d0*rn*dndt*ct2*(var2+var1*st2)
     3      +n2*2.d0*st*ct*ct2*var1+2.d0*rn*dndt*
     1      (ct2*var2+st2*var3)
     4      +2.d0*st*ct*(var3-var2)*n2
     5      -2.d0*ct*st*(var3*var1+n2*st2*var1)
     6      +2.d0*rn*dndt*st2*(ct2*var1+var3)
     7      -4.d0*rn*dndt*st2*ct2*(var1-n2)
     8      -4.d0*st*ct*n2*var1*(ct2-st2) )

c      write(*,*)'rrind new 2 deriv',deriv

c      ddenom=-st/droot+ct*(dndt/rn)/droot+st*deriv/droot
c     1 - (ct+(dndt/rn)*st)*(dndt/rn)*deriv/(droot*droot*droot)

c      write(*,*)'rrind new 1 ddenom',ddenom

cSAP100921
      ddenom1=-1.d0/droot+ct*(dndt1/rn)/droot+deriv/droot
     1 - (ct+(dndt/rn)*st)*(dndt1/rn)*deriv/(droot*droot*droot)
      ddenom=st*ddenom1

c      write(*,*)'rrind new 2 ddenom',ddenom
c      write(*,*)'rrind new ddenom1',ddenom1

c      rri2=n2*st*droot/ddenom
cSAP100921
      rri2=n2*droot/ddenom1

c      write(*,*)'rrind new rri2',rri2

      rri2=dabs(rri2)
c
cSAP1009223 for test only 
      write(*,*)'old rrind=',rrind(npar,nper,omega,omegpe)
c_endtest
c
      rrind_new_100923=dsqrt(rri2)

cSm070121
c      write(*,*)'dsqrt(n2),st,droot,ddenom',dsqrt(n2),st,droot,ddenom
       write(*,*)'new rrind rrind_new_100923',rrind_new_100923
c      write(*,*)'st,ct,dndt/rn,deriv',st,ct,dndt/rn,deriv
c      write(*,*)'-st/droot+ct*(dndt/rn)/droot+st*deriv/droot',
c     &           -st/droot+ct*(dndt/rn)/droot+st*deriv/droot
c      write(*,*)'(ct+(dndt/rn)*st)*(dndt/rn)*deriv/(droot**3)',
c     &           (ct+(dndt/rn)*st)*(dndt/rn)*deriv/(droot**3)

      return
      end

      real*8 function rrind(npar,nper,omega,omegpe)
c     input npar,nper,omega=omega/omega_ce,omegape=omegape/omegac_e
c     In this version the variable st is used.
c     It has the problem for st=0 case (nperp=0)      
    
      !implicit double precision (a-h,o-z)
c      double precision npar,npar2,nper,nper2,n2
c-----input 
      real*8 npar,nper,omega,omegpe 
c-----locals
      real*8 npar2,nper2,n2,eps1,eps2,eps3,
     &rn,st,ct,st2,ct2,var1,var2,var3,eps22,
     &dnumer,ddenom,dndt,droot,deriv,rri2,
     &ddenom1,dndt1,dnumer1
c
c     calculate ray refractive index for cold electron plasma
c
      eps1=1.d0-omegpe*omegpe/(omega*omega-1.d0)
      eps2=-omegpe*omegpe/(omega*(omega*omega-1.d0))
      eps3=1.d0-omegpe*omegpe/(omega*omega)
c
      npar2=npar*npar
c
      nper2=nper*nper
      n2=nper2+npar2
c
c     calculate ray refractive index squared
c
      rn=dsqrt(n2)
      st=nper/rn
cSAP100923
c      if (dabs(st).lt.1.d-11) st=st+1.d-10
      if (dabs(st).lt.1.d-9) then
         st=st+1.d-8
c         write(*,*)'in rrind st',st
      endif

      ct=npar/rn
      st2=st*st
      ct2=ct*ct
      var1=eps1-n2
      var2=eps3-nper2
      var3=eps1-npar2
      eps22=eps2*eps2

      dnumer=-st*ct*var1*var2+st*ct*var3*var1-eps22*st*ct
     1       +n2*var1*st*ct*(ct2-st2)

c      write(*,*)'rrind old 1 dnumer', dnumer

      ddenom=eps22*st2+n2*n2*st2*ct2-var1*var2*ct2-var3*var2
     1       -st2*var3*var1-2.0*n2*st2*ct2*var1

c      write(*,*)'rrind old ddenom',ddenom

      dndt=rn*dnumer/ddenom

c      write(*,*)'rrind old 1 dndt',dndt

      droot=dsqrt(1.d0+(dndt/rn)*(dndt/rn))

c      write(*,*)'rrind old droot',droot
c
      deriv=( var1*var2*(st2-ct2)+2.d0*rn*dndt*var2*st*ct
     1       +2.d0*rn*dndt*st2*st*ct*var1+2.0*st2*ct2*var1*n2
     2       +(ct2-st2)*var3*var1-
     1       2.d0*rn*dndt*st*ct*(ct2*var1+var3)
     3       +2.d0*st2*ct2*n2*var1
     4       +eps22*(st2-ct2)
     5       +2.d0*rn*dndt*st*ct*(ct2-st2)*(var1-n2)
     6       +n2*var1*(ct2*ct2-3.d0*ct2*st2+st2*st2-
     1       3.d0*st2*ct2) )/ddenom

c      write(*,*)'rrind old 1 deriv',deriv
c
      deriv=deriv - (dnumer/(ddenom*ddenom))*(2.d0*eps22*st*ct
     1      +n2*n2*2.d0*st*ct*(ct2-st2)+4.d0*n2*rn*dndt*ct2*st2
     2      +2.d0*ct*st*var1*var2+
     1      2.d0*rn*dndt*ct2*(var2+var1*st2)
     3      +n2*2.d0*st*ct*ct2*var1+2.d0*rn*dndt*
     1      (ct2*var2+st2*var3)
     4      +2.d0*st*ct*(var3-var2)*n2
     5      -2.d0*ct*st*(var3*var1+n2*st2*var1)
     6      +2.d0*rn*dndt*st2*(ct2*var1+var3)
     7      -4.d0*rn*dndt*st2*ct2*(var1-n2)
     8      -4.d0*st*ct*n2*var1*(ct2-st2) )
c
c      write(*,*)'rrind old 2 deriv',deriv

      ddenom=-st/droot+ct*(dndt/rn)/droot+st*deriv/droot
     1 - (ct+(dndt/rn)*st)*(dndt/rn)*deriv/(droot*droot*droot)

c      write(*,*)'rrind old 1 ddenom',ddenom

      rri2=n2*st*droot/ddenom
c      write(*,*)'rrind old rri2',rri2

      rri2=dabs(rri2)
c
c      rrind=dsqrt(rri2)
      rrind=dsqrt(rri2)
cSm070121
c      write(*,*)'dsqrt(n2),st,droot,ddenom',dsqrt(n2),st,droot,ddenom
c       write(*,*)'rrind',rrind
c      write(*,*)'st,ct,dndt/rn,deriv',st,ct,dndt/rn,deriv
c      write(*,*)'-st/droot+ct*(dndt/rn)/droot+st*deriv/droot',
c     &           -st/droot+ct*(dndt/rn)/droot+st*deriv/droot
c      write(*,*)'(ct+(dndt/rn)*st)*(dndt/rn)*deriv/(droot**3)',
c     &           (ct+(dndt/rn)*st)*(dndt/rn)*deriv/(droot**3)

c
      return
      end

      double precision function rrind_hotpl(nbulk,mass_ar,x_ar,y_ar,
     .t_av_ar,tpop_ar,vflow_ar,nll_in,np_in)
c
c     calculate ray refructive index N_ray for hot non-relativistic plasma
c     n_ray**2=abs|n**2*sin(theta)*sqrt(1+(1/n*dn/d(heta))**2)/
c                  d/dtheta{[cos(theta)+(1/n*dn/d(theta))*sin(theta)]/
c                           [1+(1/n*dn/d(theta))**2]}|
c     G.Bekefi, Radiation Processes in Plasmas,(John Wiley and Sonc Inc.) 1966
c     page 32, formula (1.121)
c
c     input
c     INPUTS:
c      nbulk the total number of plasma species
c      mass_ar(*) - the masses  of the plasma species (in electron mass) 
c      x_ar(*) = (fpe/f)**2
c      y_ar(*) = fce/f   It is the algebraic number
c      t_av_ar(*)=average temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar(*) = T_perp/T_parallel  
c      vflow_ar(*)    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c      k_sum(3,3) -  the nine components of the dielectric tensor
c                    evaluated at (mass,X,Y,Te,tpop,vflow,nll,np)
      
      implicit none
cSm030226
      include 'param.i'      
c-----input
      integer nbulk
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double precision np_in,nll_in

c-----externals
      double complex dhot_sum

c-----locals
c     derivatives from the hermitian dispersion function D
c     dd(1,js)=dD/dX_js,dd(2,js)=dD/dY_js,dd(3,js)=dD/dT_av_js
c     dd(4,js)=dD/dtpop_js,dd(5,js)=dD/dVflow_js
c     ddnp_h=dD/N_perp,ddnll_h=dD/dN_parallel
c     derivartive from the complete dispersion function 
c     ddnp=dD/N_perp
c      integer nbulka
cSm030227
c      parameter (nbulka=5)
c      double complex dd(5,nbulka) !dD/d(X_s,Y_s,Tav_s,tpop_as,V_s) work array
      double complex dd(5,nbulka) !dD/d(X_s,Y_s,Tav_s,tpop_as,V_s) work array

c-----d=D, ddnll_h=dD_herm/dN_par, ddnp_h=dD_herm/dN_perp, ddnp=dD_full/dNperp
      double complex d,ddnll_h,ddnp,ddnp_h,
     . ddnll_plus,ddnll_min,ddnp_plus,ddnp_min

c     second derivatives
      double precision  d2d_dnp2,d2d_dnlldnp,d2d_dnll2,d2d_dnpdnll,
     . d_dddnpar_dt,d_dddnper_dt

      double precision cnparp,cnperp,cn,cost,sint,denumer,ddenom,
     .dndt_devn,
     .step,cnper_plus,cnper_min,cnpar_plus,cnpar_min,
     .d_denumer_dt,d_ddenom_dt,d_dndtdevn_dt,
     .dnpardt,dnperdt,droot,ddenum1
      double complex reps(3,3)
      integer i
      
      if (nbulk.gt.nbulka) then
         write(*,*)'emission.f in rrind_hotpl nbulk.gt.nbulka'
         write(*,*)'but it should be nbulk.le.nbulka'
         write(*,*)'nbulka=',nbulka,'nbulk=',nbulk
         write(*,*)'change nbulka in param.i and recompile'
         stop 
      endif

c      write(*,*)'rrind_hotpl nbulk,nll_in,np_in',nbulk,nll_in,np_in
c      do i=1,nbulk
c        write(*,*)'i,mass_ar(i),x_ar(i),y_ar(i),
c     . t_av_ar(i),tpop_ar(i),vflow_ar(i)',
c     . i,mass_ar(i),x_ar(i),y_ar(i),t_av_ar(i),tpop_ar(i),vflow_ar(i)
c      enddo
c


      call npnllmin(nll_in,np_in,cnparp,cnperp)

c      write(*,*)'cnparp,cnperp',cnparp,cnperp

c-----Now, compute the derivative of the D(w,k) with respesct all 
c     for Hermitian tensor

c-----calculate dielectric hot tensor: reps
      d=dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .           vflow_ar,cnparp,cnperp,2,reps)
     
c-----calculate the derivatives from the hot hermition dispersion function
c     ddnll_hd=D_hermitian/dN_par
c     ddnp_hd=D_hermitian/dN_perp

      call Ddhot(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar,
     .           cnparp,cnperp,reps,dd,ddnp_h,ddnll_h,ddnp)
     
c      write(*,*)'ddnp_h,ddnll_h',ddnp_h,ddnll_h

      cn=dsqrt(cnparp**2+cnperp**2)
      sint=cnperp/cn !sin(theta)
      cost=cnparp/cn !cos(theta)

      denumer=ddnll_h*sint-ddnp_h*cost
      ddenom= ddnll_h*cost+ddnp_h*sint

c       write(*,*)'denumer,ddenom', denumer,ddenom

c-----[dN/d(theta)]/N
      dndt_devn=denumer/ddenom
c      write(*,*)' dndt_devn', dndt_devn
c-----dNpar/d(theta), dNper/d(theta)     
      dnpardt= dndt_devn*cn*cost-cn*sint
      dnperdt= dndt_devn*cn*sint+cn*cost
c       write(*,*)'dnpardt,dnperdt',dnpardt,dnperdt
c-----calculate the second derivatives 
c     dD^2/dN_par^2, dD^2/dN_per^2,dD^2/(dN_par dN_per)
c     using the numerical derivatives from the first derivatives
c     calculated analytically
      step=1.d-7

      cnper_plus=cnperp*(1.d0+step)    
      d=dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .           vflow_ar,cnparp,cnper_plus,2,reps)
     
      call Ddhot(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar,
     .     cnparp,cnper_plus,reps,dd,ddnp_plus,ddnll_plus,ddnp)

c      write(*,*)'ddnp_plus,ddnll_plus',ddnp_plus,ddnll_plus

      cnper_min =cnperp*(1.d0-step)
      d=dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .           vflow_ar,cnparp,cnper_min,2,reps)
     
      call Ddhot(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar,
     .     cnparp,cnper_min,reps,dd,ddnp_min,ddnll_min,ddnp)
c      write(*,*)'ddnp_min,ddnll_min',ddnp_min,ddnll_min

      d2d_dnp2=(ddnp_plus-ddnp_min)/(2.d0*step*cnperp)      !dD^2/dN_perp^2
      d2d_dnlldnp=(ddnll_plus-ddnll_min)/(2.d0*step*cnperp) !dD^2/dNpar*dN_perp

c      write(*,*)'d2d_dnp2,d2d_dnlldnp',d2d_dnp2,d2d_dnlldnp

      cnpar_plus=cnparp*(1.d0+step)    
      d=dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .           vflow_ar,cnpar_plus,cnperp,2,reps)
     
      call Ddhot(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar,
     .     cnpar_plus,cnperp,reps,dd,ddnp_plus,ddnll_plus,ddnp)
c      write(*,*)' 2 ddnp_plus,ddnll_plus',ddnp_plus,ddnll_plus

      cnpar_min =cnparp*(1.d0-step)
      d=dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .           vflow_ar,cnpar_min,cnperp,2,reps)
     
      call Ddhot(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar,
     .     cnpar_min,cnperp,reps,dd,ddnp_min,ddnll_min,ddnp)
c      write(*,*)'ddnp_min,ddnll_min',ddnp_min,ddnll_min

      d2d_dnll2=(ddnll_plus-ddnll_min)/(2.d0*step*cnparp)   !dD^2/dN_par^2
      d2d_dnpdnll=(ddnp_plus-ddnp_min)/(2.d0*step*cnparp) !dD^2/dNper*dN_par
c      write(*,*)'d2d_dnll2,d2d_dnpdnll',d2d_dnll2,d2d_dnpdnll

c-----calculate d(dD/dN_par)/d(theta), d(dD/dN_per)/d(theta)
      d_dddnpar_dt=d2d_dnll2*dnpardt + d2d_dnlldnp*dnperdt
      d_dddnper_dt=d2d_dnlldnp*dnpardt + d2d_dnp2*dnperdt
c      write(*,*)'d_dddnpar_dt,d_ddnper_dt',d_dddnpar_dt,d_dddnper_dt
c-----d(denumer)/d(theta)
      d_denumer_dt=d_dddnpar_dt*sint+cost*ddnll_h-
     .             d_dddnper_dt*cost+sint*ddnp_h
c      write(*,*)'d_denumer_dt',d_denumer_dt
c-----d(ddenom)/d(theta)
      d_ddenom_dt=d_dddnpar_dt*cost-sint*ddnll_h+
     .            d_dddnper_dt*sint+cost*ddnp_h
c      write(*,*)'d_ddenom_dt',d_ddenom_dt
c-----d[1/n*dN/d(theta)]/d(theta)
      d_dndtdevn_dt=d_denumer_dt/ddenom-denumer*d_ddenom_dt/ddenom**2
c      write(*,*)'d_dndtdevn_dt',d_dndtdevn_dt
c-----ddenum1=d/dt{[cos(theta)+(1/n*dN/d(theta))*sin(theta)]/
c               sqrt(1+(1/n*dN/d(theta))**2)

      droot=dsqrt(1.d0+dndt_devn**2)
      ddenum1=(-sint+cost*dndt_devn+sint*d_dndtdevn_dt)/droot-
     .(cost+dndt_devn*sint)/droot**3*dndt_devn*d_dndtdevn_dt
c      write(*,*)'sint,droot,ddenum1',sint,droot,ddenum1

c     calculate ray refractive index
c
     
      rrind_hotpl=cn*dsqrt(dabs(sint*droot/ddenum1))

cSm070121      
c      write(*,*)' rrind_hotpl', rrind_hotpl
c      write(*,*)'cn,sint,droot,ddenum1',cn,sint,droot,ddenum1
c      write(*,*)'sint,cost,dndt_devn,d_dndtdevn_dt',
c     &           sint,cost,dndt_devn,d_dndtdevn_dt
c      write(*,*)'(-sint+cost*dndt_devn+sint*d_dndtdevn_dt)/droot',
c     &(-sint+cost*dndt_devn+sint*d_dndtdevn_dt)/droot
c      write(*,*)'(cost+dndt_devn*sint)/droot**3*dndt_devn*d_dndtdevn_dt'
c     &,(cost+dndt_devn*sint)/droot**3*dndt_devn*d_dndtdevn_dt

c  
      return
      end


      subroutine spl_ray     
c----------------------------------------------------------------------
c     creation of arrays of the spline coefficients                  
c     for radius and refractive vectors coordinates along the ray 	      
c     tr_r(nrelta4),cx_r(ndens4)                  		      
c     tr_z(nrelta4),cx_z(ndens4)
c     tr_phi(nrelta4),cx_phi(ndens4)                  		      
c     tr_cnr(nrelta4),cx_cnr(ndens4)
c     tr_cnz(nrelta4),cx_cnz(ndens4)                  		      
c     tr_cm(nrelta4),cx_cm(ndens4)
c
c     thise program uses the following subroutines:  		      
c								      
c     iac1r(rhom,ip,ip4,denm,lx,mxa,fxa,mxb,fxb,trdens,cxdens)	      
c     calculates the spline coefficients for 1d function	
c------------------------------------------------------------------
c     input data are from param.i, emission.i
c     output data are  in emission.i     	
c-----------------------------------------------------------------
      implicit none
      !implicit double precision (a-h,o-z)
      include 'param.i'
      include 'emissa.i'
      include 'write.i'
      include 'oxb.i'
      include 'one.i'     
c-----input in write.i
c      integer nrayelt! the number of points along the ray. The first point is=1
c-----locals
      real*8 fxa,fxb
      integer ip,ip4,j,lx,mxa,mxb
c-----external
      real*8 ias1r
      real*8 r_ray,z_ray,phi_ray,cnr_ray,cnz_ray,cm_ray

      write(*,*)'in spl_ray nrayelt,nrelta',nrayelt,nrelta

      if (nrayelt.gt.nrelta) then
        write(*,*)'in emission.df in spl_ray nrayelt>nrelta'
        write(*,*)'nrelta=',nrelta,'nrayelt=',nrayelt
        write(*,*)'increase nrelta in param.i'
	stop
      endif

c-----calculation of the spline coefficients for radius and refractive vectors
      write(*,*)'emission.f in spl_raynrayelt_o_cutoff',nrayelt_o_cutoff 
      if (nrayelt_o_cutoff.le.1) then
c-------one ray, no O-X jump for OX conversion
        ip=nrayelt  
        ip4=nrayelt+4
        lx=1
        mxa=0
        fxa=0.0
        mxb=0
        fxb=0.0
 
        call iac1r(wsn,ip,ip4,wz,lx,mxa,fxa,mxb,fxb,tr_z,cx_z)
        call iac1r(wsn,ip,ip4,wr,lx,mxa,fxa,mxb,fxb,tr_r,cx_r)
        call iac1r(wsn,ip,ip4,wphi,lx,mxa,fxa,mxb,fxb,tr_phi,cx_phi)
    
        call iac1r(wsn,ip,ip4,wcnz_em,lx,mxa,fxa,mxb,fxb,tr_cnz,cx_cnz)
        call iac1r(wsn,ip,ip4,wcnr_em,lx,mxa,fxa,mxb,fxb,tr_cnr,cx_cnr)
        call iac1r(wsn,ip,ip4,wcm_em,lx,mxa,fxa,mxb,fxb,tr_cm,cx_cm)
	
c for test
c        write(*,*)'spl_ray testno O_X jump  nrayelt',nrayelt

c        do j=1,nrayelt
c	  s=wsn(j)
c	  write(*,*)'emission test j,s',j,s
c          pz=z_ray(s)
c          pr=r_ray(s)
c          pphi=phi_ray(s)

c          pcnz=cnz_ray(s)
c          pcnr=cnr_ray(s)
c          pcm=cm_ray(s)

c          write(*,*)'pz,wz(j)',pz,wz(j)
c          write(*,*)'pr,wr(j)',pr,wr(j)
c          write(*,*)'pphi,wphi(j)',pphi,wphi(j)

c          write(*,*)'pcnz,wcnz_em(j)',pcnz,wcnz_em(j)
c          write(*,*)'pcnr,wcnr_em(j)',pcnr,wcnr_em(j)
c          write(*,*)'pcm,wcm_em(j)',pcm,wcm_em(j)
c        enddo

c        do j=1,nrayelt-1
c          s=0.5d0*(wsn(j)+wsn(j+1))
c	  write(*,*)'emission test j,s',j,s
c          pz=z_ray(s)
c          pr=r_ray(s)
c          pphi=phi_ray(s)

c          pcnz=cnz_ray(s)
c          pcnr=cnr_ray(s)
c          pcm=cm_ray(s)

c          write(*,*)'wz(j),pz,wz(j+0.5),wz(j+1)',
c     +    wz(j),pz,0.5d0*(wz(j)+wz(j+1)),wz(j+1)
c          write(*,*)'wr(j),pr,wr(j+0.5),wr(j+1)',
c     +    wr(j),pr,0.5d0*(wr(j)+wr(j+1)),wr(j+1)
c          write(*,*)'wphi(j),pphi,wphi(j+0.5),wphi(j+1)',
c     +    wphi(j),pphi,0.5d0*(wphi(j)+wphi(j+1)),wphi(j+1)

c          write(*,*)'wcnz_em(j),pcnz,wcnz_em(j+0.5),wcnz_em(j+1)',
c     +    wcnz_em(j),pcnz,0.5d0*(wcnz_em(j)+wcnz_em(j+1)),wcnz_em(j+1)
c          write(*,*)'wcnr_em(j),pcnr,wcnr_em(j+0.5),wcnr_em(j+1)',
c     +    wcnr_em(j),pcnr,0.5d0*(wcnr_em(j)+wcnr_em(j+1)),wcnr_em(j+1)
c          write(*,*)'wcm_em(j),pcm,wcm_em(j+0.5),wcm_em(j+1)',
c     +    wcm_em(j),pcm,0.5d0*(wcm_em(j)+wcm_em(j+1)),wcm_em(j+1)
        
c        enddo
cendtest
      else
c--------------------------------------------------------
c       ray has O-X jump for OX conversion,
c       this point divides the ray for two parts
c       O mode ray: is [1:nrayelt_o_cutoff]
c       X_EBW mode ray: is [nrayelt_o_cutoff+1,nrayelt]
c-------------------------------------------------------
c       O - ray
        do j=1,nrayelt
c          write(*,*)'j,wsn(j),wr(j)',j,wsn(j),wr(j)
        enddo

        do j=1,nrayelt_o_cutoff
          wsn_o(j)=wsn(j)
          wz_o(j)=wz(j)
          wr_o(j)=wr(j)
          wphi_o(j)=wphi(j)
          wcnz_em_o(j)=wcnz_em(j)
          wcnr_em_o(j)=wcnr_em(j)
          wcm_em_o(j)=wcm_em(j)
        enddo

        ip=nrayelt_o_cutoff
        ip4=ip+4
        lx=1
        mxa=0
        fxa=0.0
        mxb=0
        fxb=0.0
 
        call iac1r(wsn_o,ip,ip4,wz_o,lx,mxa,fxa,mxb,fxb,tr_z_o,cx_z_o)
        call iac1r(wsn_o,ip,ip4,wr_o,lx,mxa,fxa,mxb,fxb,tr_r_o,cx_r_o)
        call iac1r(wsn_o,ip,ip4,wphi_o,lx,mxa,fxa,mxb,fxb,tr_phi_o,
     &  cx_phi_o)
    
        call iac1r(wsn_o,ip,ip4,wcnz_em_o,lx,mxa,fxa,mxb,fxb,tr_cnz_o,
     &  cx_cnz_o)
        call iac1r(wsn_o,ip,ip4,wcnr_em_o,lx,mxa,fxa,mxb,fxb,tr_cnr_o,
     &  cx_cnr_o)
        call iac1r(wsn_o,ip,ip4,wcm_em_o,lx,mxa,fxa,mxb,fxb,tr_cm_o,
     &  cx_cm_o)
	
c for test
c        write(*,*)'spl_ray O_mode test nrayelt_o_cutoff',
c     &             nrayelt_o_cutoff

c        p_norm_r=0.d0
c        p_norm_z=0.d0
c        p_norm_phi=0.d0
c        p_norm_nr=0.d0
c        p_norm_nz=0.d0
c        p_norm_m=0.d0

c        do j=1,nrayelt_o_cutoff
c         write(*,*)'j,wsn(j),wr_o(j)',j,wsn(j),wr_o(j)
c        enddo


c        do j=1,nrayelt_o_cutoff
c	  s=wsn_o(j)
c	  write(*,*)'emission o-mode test j,s',j,s
c          pz=z_ray(s)
c          pr=r_ray(s) 
c          pr=ias1r(tr_r_o,nrayelt_o_cutoff,nrayelt_o_cutoff+4,
c     &             cx_r_o,0,s)
c          pphi=phi_ray(s)

c          pcnz=cnz_ray(s)
c          pcnr=cnr_ray(s)
c          pcm=cm_ray(s)

c          write(*,*)'pz,wz_o(j)',pz,wz_o(j)
c          write(*,*)'pr,wr_o(j)',pr,wr_o(j)
c          write(*,*)'pphi,wphi_o(j)',pphi,wphi_o(j)

c          write(*,*)'pcnz,wcnz_em_o(j)',pcnz,wcnz_em_o(j)
c          write(*,*)'pcnr,wcnr_em_o(j)',pcnr,wcnr_em_o(j)
c          write(*,*)'pcm,wcm_em_o(j)',pcm,wcm_em_o(j)
c          p_norm_r=p_norm_r+(pr-wr_o(j))**2
c          p_norm_z=p_norm_z+(pz-wz_o(j))**2
c          write(*,*)'j,pz,wz_o(j),p_norm_z',j,pz,wz_o(j),p_norm_z
c          p_norm_phi=p_norm_phi+(pphi-wphi_o(j))**2
c          p_norm_nr=p_norm_nr+(pcnr-wcnr_em_o(j))**2
c          p_norm_nz=p_norm_nz+(pcnz-wcnz_em_o(j))**2
c          p_norm_m=p_norm_m+(pcm-wcm_em_o(j))**2

c        enddo
c        write(*,*)'p_norm_r,p_norm_z,p_norm_phi',
c     &  p_norm_r,p_norm_z,p_norm_phi
c        write(*,*)'p_norm_r',dsqrt(p_norm_r/nrayelt_o_cutoff)
c        write(*,*)'p_norm_z',dsqrt(p_norm_z/nrayelt_o_cutoff)
c        write(*,*)'p_norm_phi',dsqrt(p_norm_phi/nrayelt_o_cutoff)
c        write(*,*)'p_norm_nr',dsqrt(p_norm_nr/nrayelt_o_cutoff)
c        write(*,*)'p_norm_nz',dsqrt(p_norm_nz/nrayelt_o_cutoff)
c        write(*,*)'p_norm_m',dsqrt(p_norm_m/nrayelt_o_cutoff)



c        do j=1,nrayelt_o_cutoff-1
c          s=0.5d0*(wsn_o(j)+wsn_o(j+1))
c          write(*,*)'emission o-mode test j,s',j,s
c          pz=z_ray(s)
c          pr=r_ray(s)
c          pphi=phi_ray(s)

c          pcnz=cnz_ray(s,nrayelt)
c          pcnr=cnr_ray(s,nrayelt)
c          pcm=cm_ray(s)

c         write(*,*)'wz(j),pz,wz(j+0.5),wz(j+1)',
c     +   wz(j),pz,0.5d0*(wz(j)+wz(j+1)),wz(j+1)
c         write(*,*)'wr(j),pr,wr(j+0.5),wr(j+1)',
c     +   wr(j),pr,0.5d0*(wr(j)+wr(j+1)),wr(j+1)
c         write(*,*)'wphi(j),pphi,wphi(j+0.5),wphi(j+1)',
c     +   wphi(j),pphi,0.5d0*(wphi(j)+wphi(j+1)),wphi(j+1)

c         write(*,*)'wcnz_em(j),pcnz,wcnz_em(j+0.5),wcnz_em(j+1)',
c     +   wcnz_em(j),pcnz,0.5d0*(wcnz_em(j)+wcnz_em(j+1)),wcnz_em(j+1)
c         write(*,*)'wcnr_em(j),pcnr,wcnr_em(j+0.5),wcnr_em(j+1)',
c     +   wcnr_em(j),pcnr,0.5d0*(wcnr_em(j)+wcnr_em(j+1)),wcnr_em(j+1)
c         write(*,*)'wcm_em(j),pcm,wcm_em(j+0.5),wcm_em(j+1)',
c     +   wcm_em(j),pcm,0.5d0*(wcm_em(j)+wcm_em(j+1)),wcm_em(j+1)
        
c       enddo
cendtest     ip=nrayelt 

c       X_EBW - ray

        do j=1,nrayelt-nrayelt_o_cutoff
          wsn_x(j)=wsn(j+nrayelt_o_cutoff)
          wz_x(j)=wz(j+nrayelt_o_cutoff)
          wr_x(j)=wr(j+nrayelt_o_cutoff)
          wphi_x(j)=wphi(j+nrayelt_o_cutoff)
          wcnz_em_x(j)=wcnz_em(j+nrayelt_o_cutoff)
          wcnr_em_x(j)=wcnr_em(j+nrayelt_o_cutoff)
          wcm_em_x(j)=wcm_em(j+nrayelt_o_cutoff)
        enddo

        ip=nrayelt-nrayelt_o_cutoff
        ip4=ip+4
        lx=1
        mxa=0
        fxa=0.0
        mxb=0
        fxb=0.0
 
        call iac1r(wsn_x,ip,ip4,wz_x,lx,mxa,fxa,mxb,fxb,tr_z_x,cx_z_x)
        call iac1r(wsn_x,ip,ip4,wr_x,lx,mxa,fxa,mxb,fxb,tr_r_x,cx_r_x)
        call iac1r(wsn_x,ip,ip4,wphi_x,lx,mxa,fxa,mxb,fxb,tr_phi_x,
     &             cx_phi_x)
    
        call iac1r(wsn_x,ip,ip4,wcnz_em_x,lx,mxa,fxa,mxb,fxb,tr_cnz_x,
     &             cx_cnz_x)
        call iac1r(wsn_x,ip,ip4,wcnr_em_x,lx,mxa,fxa,mxb,fxb,tr_cnr_x,
     &             cx_cnr_x)
        call iac1r(wsn_x,ip,ip4,wcm_em_x,lx,mxa,fxa,mxb,fxb,tr_cm_x,
     &             cx_cm_x)
	
c for test
c       write(*,*)'spl_ray x-mode test nrayelt_o_cutoff',nrayelt_o_cutoff

c       do j=nrayelt_o_cutoff+1,nrayelt
c	 s=wsn(j)
c	 write(*,*)'emission test j,s',j,s
c         pz=z_ray(s)
c         pr=r_ray(s)
c         pphi=phi_ray(s)

c         pcnz=cnz_ray(s)
c         pcnr=cnr_ray(s)
c         pcm=cm_ray(s)

c         write(*,*)'pz,wz(j)',pz,wz(j)
c         write(*,*)'pr,wr(j)',pr,wr(j)
c         write(*,*)'pphi,wphi(j)',pphi,wphi(j)

c         write(*,*)'pcnz,wcnz_em(j)',pcnz,wcnz_em(j)
c         write(*,*)'pcnr,wcnr_em(j)',pcnr,wcnr_em(j)
c         write(*,*)'pcm,wcm_em(j)',pcm,wcm_em(j)
c       enddo

c       do j=nrayelt_o_cutoff+1,nrayelt-1
c         s=0.5d0*(wsn(j)+wsn(j+1))
c	 write(*,*)'emission test j,s',j,s
c         pz=z_ray(s)
c         pr=r_ray(s)
c         pphi=phi_ray(s)

c         pcnz=cnz_ray(s)
c         pcnr=cnr_ray(s)
c         pcm=cm_ray(s)

c         write(*,*)'wz(j),pz,wz(j+0.5),wz(j+1)',
c     +   wz(j),pz,0.5d0*(wz(j)+wz(j+1)),wz(j+1)
c         write(*,*)'wr(j),pr,wr(j+0.5),wr(j+1)',
c     +   wr(j),pr,0.5d0*(wr(j)+wr(j+1)),wr(j+1)
c         write(*,*)'wphi(j),pphi,wphi(j+0.5),wphi(j+1)',
c     +   wphi(j),pphi,0.5d0*(wphi(j)+wphi(j+1)),wphi(j+1)

c         write(*,*)'wcnz_em(j),pcnz,wcnz_em(j+0.5),wcnz_em(j+1)',
c     +   wcnz_em(j),pcnz,0.5d0*(wcnz_em(j)+wcnz_em(j+1)),wcnz_em(j+1)
c         write(*,*)'wcnr_em(j),pcnr,wcnr_em(j+0.5),wcnr_em(j+1)',
c     +   wcnr_em(j),pcnr,0.5d0*(wcnr_em(j)+wcnr_em(j+1)),wcnr_em(j+1)
c         write(*,*)'wcm_em(j),pcm,wcm_em(j+0.5),wcm_em(j+1)',
c     +   wcm_em(j),pcm,0.5d0*(wcm_em(j)+wcm_em(j+1)),wcm_em(j+1)
        
c       enddo
cendtest
      endif

c      stop 'spl_ray'

      return
      end


      real*8 FUNCTION r_ray(s)
c     double precision FUNCTION r_ray(s,nrayelt)
c-------------------------------------------------------
c     major radius r on the length s (cm) along the ray (from spline)
c-------------------------------------------------------
      !implicit double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'emissa.i'
      include 'oxb.i'
      include 'write.i'
c-----input
c      integer nrayelt
      real*8 s
c-----external
      real*8 ias1r
c-----locals
      integer idx

      idx=0
      
      if (nrayelt_o_cutoff.le.1) then !no O-X jump
         r_ray=ias1r(tr_r,nrayelt,nrayelt+4,cx_r,idx,s)
      else
         !It was O_X jump
         if (s.le.wsn_o(nrayelt_o_cutoff)) then !o mode ray
           r_ray=ias1r(tr_r_o,nrayelt_o_cutoff,nrayelt_o_cutoff+4,
     &                 cx_r_o,idx,s)
         else
           if(s.ge.wsn_x(1)) then !x-ebw mode ray
              r_ray=ias1r(tr_r_x,nrayelt-nrayelt_o_cutoff,
     &                    nrayelt-nrayelt_o_cutoff+4,
     &                    cx_r_x,idx,s)
           else
             ! s is inside jump area
             write(*,*)'emission.f in r_ray O-X jump was'
             write(*,*)'nrayelt_o_cutoff ', nrayelt_o_cutoff
             write(*,*)'s is inside jump area'
             write(*,*)'s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)'
     &       ,s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)
           endif
         endif
      endif

      return
      end



      real*8 FUNCTION z_ray(s)
c-------------------------------------------------------
c     vertical coordinate z  on the length s (cm) along the ray (from spline)
c-------------------------------------------------------
      implicit none
      include 'param.i'
      include 'emissa.i'
      include 'oxb.i'
      include 'write.i'
c-----input
      real*8 s
c-----external
      real*8 ias1r
c-----locals
      integer idx

      idx=0
 
      if (nrayelt_o_cutoff.le.1) then !no O-X jump
         z_ray=ias1r(tr_z,nrayelt,nrayelt+4,cx_z,idx,s)
      else
         !It was O_X jump
         if (s.le.wsn_o(nrayelt_o_cutoff)) then !o mode ray
           z_ray=ias1r(tr_z_o,nrayelt_o_cutoff,nrayelt_o_cutoff+4,
     &                 cx_z_o,idx,s)
         else
           if(s.ge.wsn_x(1)) then !x-ebw mode ray
              z_ray=ias1r(tr_z_x,nrayelt-nrayelt_o_cutoff,
     &                    nrayelt-nrayelt_o_cutoff+4,
     &                    cx_z_x,idx,s)
           else
             ! s is inside jump area
             write(*,*)'emission.f in z_ray O-X jump was'
             write(*,*)'nrayelt_o_cutoff ', nrayelt_o_cutoff
             write(*,*)'nrayelt_o_cutoff ', nrayelt_o_cutoff
             write(*,*)'s is inside jump area'
             write(*,*)'s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)'
     &       ,s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)
           endif
         endif
      endif

      return
      end


      real*8 FUNCTION phi_ray(s)
c-------------------------------------------------------
c     toroidakl angle phi (radian) on the length s (cm) along the ray (from spline)
c-------------------------------------------------------
      implicit none
      include 'param.i'
      include 'emissa.i'
      include 'oxb.i'
      include 'write.i'
c-----input
      real*8 s
c-----external
      real*8 ias1r
c-----locals
      integer idx

      idx=0

      if (nrayelt_o_cutoff.le.1) then !no O-X jump
         phi_ray=ias1r(tr_phi,nrayelt,nrayelt+4,cx_phi,idx,s)
      else
         !It was O_X jump
         if (s.le.wsn_o(nrayelt_o_cutoff)) then !o mode ray
           phi_ray=ias1r(tr_phi_o,nrayelt_o_cutoff,nrayelt_o_cutoff+4,
     &                 cx_phi_o,idx,s)
         else
           if(s.ge.wsn_x(1)) then !x-ebw mode ray
              phi_ray=ias1r(tr_phi_x,nrayelt-nrayelt_o_cutoff,
     &                    nrayelt-nrayelt_o_cutoff+4,
     &                    cx_phi_x,idx,s)
           else
             ! s is inside jump area
             write(*,*)'emission.f in phi_ray O-X jump was'
             write(*,*)'nrayelt_o_cutoff ', nrayelt_o_cutoff
             write(*,*)'nrayelt_o_cutoff ', nrayelt_o_cutoff
             write(*,*)'s is inside jump area'
             write(*,*)'s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)'
     &       ,s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)
           endif
         endif
      endif

      return
      end


      real*8 FUNCTION cnr_ray(s)
c-------------------------------------------------------
c     radial refructive index N_r on the length s (sm) 
c     along the ray (from spline)
c-------------------------------------------------------
      implicit none
      include 'param.i'
      include 'emissa.i'
      include 'oxb.i'
      include 'write.i'

c-----input
      real*8 s
c-----external
      real*8 ias1r
c-----locals
      integer idx

      idx=0

      if (nrayelt_o_cutoff.le.1) then !no O-X jump
          cnr_ray=ias1r(tr_cnr,nrayelt,nrayelt+4,cx_cnr,idx,s)
      else
         !It was O_X jump
         if (s.le.wsn_o(nrayelt_o_cutoff)) then !o mode ray
           cnr_ray=ias1r(tr_cnr_o,nrayelt_o_cutoff,nrayelt_o_cutoff+4,
     &                 cx_cnr_o,idx,s)
         else
           if(s.ge.wsn_x(1)) then !x-ebw mode ray
              cnr_ray=ias1r(tr_cnr_x,nrayelt-nrayelt_o_cutoff,
     &                    nrayelt-nrayelt_o_cutoff+4,
     &                    cx_cnr_x,idx,s)
           else
             ! s is inside jump area
             write(*,*)'emission.f in cnr_ray O-X jump was'
             write(*,*)'nrayelt_o_cutoff ', nrayelt_o_cutoff
             write(*,*)'s is inside jump area'
             write(*,*)'s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)'
     &       ,s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)
           endif
         endif
      endif

      return
      end


      real*8 FUNCTION cnz_ray(s)
c-------------------------------------------------------
c     vertical refractive index N_z radius on the length s (sm)
c     along the ray (from spline)
c-------------------------------------------------------
      implicit none
      include 'param.i'
      include 'emissa.i'
      include 'oxb.i'
      include 'write.i'

c-----input
      real*8 s
c-----external
      real*8 ias1r
c-----locals
      integer idx

      idx=0

      if (nrayelt_o_cutoff.le.1) then !no O-X jump
          cnz_ray=ias1r(tr_cnz,nrayelt,nrayelt+4,cx_cnz,idx,s)
      else
         !It was O_X jump
         if (s.le.wsn_o(nrayelt_o_cutoff)) then !o mode ray
           cnz_ray=ias1r(tr_cnz_o,nrayelt_o_cutoff,nrayelt_o_cutoff+4,
     &                 cx_cnz_o,idx,s)
         else
           if(s.ge.wsn_x(1)) then !x-ebw mode ray
              cnz_ray=ias1r(tr_cnz_x,nrayelt-nrayelt_o_cutoff,
     &                    nrayelt-nrayelt_o_cutoff+4,
     &                    cx_cnz_x,idx,s)
           else
             ! s is inside jump area
             write(*,*)'emission.f in cnz_ray O-X jump was'
             write(*,*)'nrayelt_o_cutoff ', nrayelt_o_cutoff
             write(*,*)'s is inside jump area'
             write(*,*)'s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)'
     &       ,s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)
           endif
         endif
      endif

      return
      end


      real*8 FUNCTION cm_ray(s)
c-------------------------------------------------------
c     toroidal refractive index M=r*N_phi on the length s (sm)
c     along the ray (from spline)
c-------------------------------------------------------
      implicit none
      include 'param.i'
      include 'emissa.i' 
      include 'oxb.i'
      include 'write.i'

c-----input
      real*8 s
c-----external
      real*8 ias1r
c-----locals
      integer idx

      idx=0 

c      write(*,*)'cm_ray s,nrayelt_o_cutoff',s,nrayelt_o_cutoff

      if (nrayelt_o_cutoff.le.1) then !no O-X jump
         cm_ray=ias1r(tr_cm,nrayelt,nrayelt+4,cx_cm,idx,s)
      else
         !It was O_X jump
         if (s.le.wsn_o(nrayelt_o_cutoff)) then !o mode ray
           cm_ray=ias1r(tr_cm_o,nrayelt_o_cutoff,nrayelt_o_cutoff+4,
     &                 cx_cm_o,idx,s)
         else
           if(s.ge.wsn_x(1)) then !x-ebw mode ray
              cm_ray=ias1r(tr_cm_x,nrayelt-nrayelt_o_cutoff,
     &                    nrayelt-nrayelt_o_cutoff+4,
     &                    cx_cm_x,idx,s)
           else
             ! s is inside jump area
             write(*,*)'emission.f in cm_ray O-X jump was'
             write(*,*)'nrayelt_o_cutoff ', nrayelt_o_cutoff
             write(*,*)'s is inside jump area'
             write(*,*)'s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)'
     &       ,s,wsn(nrayelt_o_cutoff),wsn(nrayelt_o_cutoff+1)
           endif
         endif
      endif
      return
      end

      subroutine shift(n_midway,nrayelt_emis,nrelta,ar)
c-----shift of the elements of array ar(k) to ar(k+1) k=n_midway,...,nrayelt_emis
      implicit none
c-----input
      integer n_midway,nrayelt_emis
      integer nrelta ! max number of elements in ar(brelta)
      double precision ar(*) ! the input and output array
c-----output ar(*)

c-----local
      integer k
      double precision p
      if ((nrayelt_emis+1).gt.nrelta) then
         write(*,*)'in emmision.f in shift (nrayelt_emis+1).gt.nrelta'
         write(*,*)'it should be: (nrayelt_emis+1).le.nrelta'
         write(*,*)'increase nrelta in param.i'
         stop
      endif

      do k=nrayelt_emis,n_midway,-1
         ar(k+1)=ar(k)
      enddo 

     

      return
      end

      subroutine add_ray_point(n_midway,nrayelt_emis_l)
c-----add the midway ray point along the ray for the emission calcultions
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'emissa.i'
      include 'write.i'
      include 'ions.i'

c-----input
      integer n_midway ! the number of the additional midway point
      integer nrayelt_emis_l ! the input nuber of the points along the ray
c-----output
c     nrayelt_emis (in common/write/ in write.i file
c
c-----external
      double precision r_ray,z_ray,phi_ray,cnz_ray,cnr_ray,cm_ray,
     +rrind,x,y,b,tempe
c-----local
      double precision u(6),deru(6),bf(3),vgr(3)
      double complex cex,cey,cez,cflown
      integer i_fkin

      double precision x_ar(nbulka),y_ar(nbulka),
     .t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka)
cfor test
      integer itest
      double precision tal,tj,temp_rad

      if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c--------usage of the analytical relativistic function and its derivatives 
         i_fkin=0
      else 
c--------usage of the mech relativistic function and its derivatives
         i_fkin=1
      endif
      
      nrayelt_emis=nrayelt_emis_l

c-----shift of the arrays points ar(k) to ar(k+1) k=n_midway,nrelt_emis 

c-----wsn the length along the ray         
      call shift(n_midway,nrayelt_emis,nrelta,wsn) 
c      do k=1,nrayelt_emis+1
c        write(*,*)'add_ray_point after shift k,`wsn(k)',k,wsn(k)
c      enddo 
c-----radius vector cordinates along the ray
      call shift(n_midway,nrayelt_emis,nrelta,wr_em)
      call shift(n_midway,nrayelt_emis,nrelta,wz_em)
      call shift(n_midway,nrayelt_emis,nrelta,wphi_em)
c-----refractive index cordinates along the ray
      call shift(n_midway,nrayelt_emis,nrelta,wcnz_em)
      call shift(n_midway,nrayelt_emis,nrelta,wcnr_em)
      call shift(n_midway,nrayelt_emis,nrelta,wcm_em)
c-----ray refractive index
      call shift(n_midway,nrayelt_emis,nrelta,wnray)
c-----absorption and emission coefficients
      call shift(n_midway,nrayelt_emis,nrelta,wal_emis)
      call shift(n_midway,nrayelt_emis,nrelta,wj_emis)
c-----emission In from one bib s_n<s<s+n+1 at s=s_n   
      call shift(n_midway,nrayelt_emis,nrelta,win_sn)
c-----plasma and radiation temperature
      call shift(n_midway,nrayelt_emis,nrelta,wtemp_em)
      call shift(n_midway,nrayelt_emis,nrelta,wtemp_rad_em)
c-----tau_n frm the n-th bin
      call shift(n_midway,nrayelt_emis,nrelta,wtaun_em)

c*****the additional midway point
c      write(*,*)'add_ray_poin n_midway,wsn(n_midway-1),wsn(n_midway+1)',
c     +n_midway,wsn(n_midway-1),wsn(n_midway+1)
     
      wsn(n_midway)=0.5d0*(wsn(n_midway-1)+wsn(n_midway+1)) ! midway length
     
c      write(*,*)'add_ray_point new point:wsn(n_midway)',wsn(n_midway)
    
c      do k=1,nrayelt_emis+1
c        write(*,*)'add_ray_point !!! with new point k,wsn(k),wr_em(k)',
c     +  k,wsn(k),wr_em(k)
c      enddo 
c-----midway point cordinates
     
c      for test
c       itest=1
c 100   continue
c       if (itest.eq.2) goto 300
c      end test

      wr_em(n_midway)=r_ray(wsn(n_midway))/(100.d0*r0x)
      wz_em(n_midway)=z_ray(wsn(n_midway))/(100.d0*r0x)
      wphi_em(n_midway)=phi_ray(wsn(n_midway))
      wcnr_em(n_midway)=cnr_ray(wsn(n_midway))
      wcnz_em(n_midway)=cnz_ray(wsn(n_midway))
      wcm_em(n_midway)=cm_ray(wsn(n_midway))

cfor test
c 300  continue
cendtest

c      write(*,*)'wr_em(n_midway-1),wr_em(n_midway),wr_em(n_midway+1)',
c     +n_midway,wr_em(n_midway-1),wr_em(n_midway),wr_em(n_midway+1)
c      write(*,*)'wz_em(n_midway-1),wz_em(n_midway),wz_em(n_midway+1)',
c     +wz_em(n_midway-1),wz_em(n_midway),wz_em(n_midway+1)
c      write(*,*)'wcnr_em(n_midway-1),wcnr_em(n_midway),
c     +wcnr_em(n_midway+1)',
c     +wcnr_em(n_midway-1),wcnr_em(n_midway),wcnr_em(n_midway+1)
c      write(*,*)'wcnz_em(n_midway-1),wcnz_em(n_midway),
c     +wcnz_em(n_midway+1)',
c     +wcnz_em(n_midway-1),wcnz_em(n_midway),wcnz_em(n_midway+1)
c      write(*,*)'wcm_em(n_midway-1),wcm_em(n_midway),
c     +wcm_em(n_midway+1)',
c     +wcm_em(n_midway-1),wcm_em(n_midway),wcm_em(n_midway+1)

c-----normalised r and z 
      r=wr_em(n_midway)
      z=wz_em(n_midway)
      phi=wphi_em(n_midway)
c      write(*,*)'wr_em(n_midway),r',wr_em(n_midway),r
c      write(*,*)'wz_em(n_midway),z',wz_em(n_midway),z
c      write(*,*)'wphi_em(n_midway),phi',wphi_em(n_midway),phi
      cnz=wcnz_em(n_midway)
      cnr=wcnr_em(n_midway)
      cm=wcm_em(n_midway)

c-----emission and absorption in the midway point
 
c-----calculate the data for emission
      bmod=b(z,r,phi)
c      write(*,*)'rho',rho
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)    !it will be used as negative for electron
      T_kev=tempe(z,r,phi,1)  
      wtemp_em(n_midway)=T_kev ! for plotting

      u(1)=z
      u(2)=r
      u(3)=phi
      u(4)=cnz
      u(5)=cnr
      u(6)=cm

      bf(1)=bz
      bf(2)=br
      bf(3)=bphi

      gam=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam)
      dc=dcos(gam)
      cnt=cn(r,cnz,cnr,cm)
      cnpar=cnt*dc
      cnper=cnt*ds

      i_geom_optic_loc=i_geom_optic
      i_geom_optic=1 !to get group velocity in deru 
      call rside1(0.d0,u,deru)
      i_geom_optic=i_geom_optic_loc

c      write(*,*)'after rside u',u
c      write(*,*)'after rside deru',deru

      vgr(1)=deru(1)
      vgr(2)=deru(2)
      vgr(3)=r*deru(3)
      vgrmods=vgr(1)**2+vgr(2)**2+vgr(3)**2
      vgroup=dsqrt(vgrmods) !in clightc
c      write(*,*)'vgroup',vgroup
c-----calculate cnray
      do i=1,nbulk
         x_ar(i)=x(z,r,phi,i)
	 y_ar(i)=y(z,r,phi,i)
         if(i.eq.1) y_ar(1)=-y_ar(1)
	 te=tempe(z,r,phi,i) ! kev
	 t_av_ar(i)=te*1000.d0      ! ev 
         tpop_ar(i)=tpop_zrp(z,r,phi,i) !YuP[2024-08-14] was tpoprho(rho,i)
         vflow_ar(i)=vflowrho(rho,i)
      enddo  
c      write(*,*)'emssion in  add_ray_point before rayrind rho',rho  
      call rayrind(i_rrind,nbulk,dmas,x_ar,y_ar,
     .             t_av_ar,tpop_ar,vflow_ar,cnpar,cnper,cnray)
c      write(*,*)'emssion in  add_ray_point after rayrind'  

c      omegpedce=dsqrt(xe)/dabs(ye) !omega_pe/omega_ce
c      omegdce=1.d0/dabs(ye)        !omega/omega_ce
c      cnray=rrind(cnpar,cnper,omegdce,omegpedce)
c      cnn=dsqrt(cnpar**2+cnper**2)
c      write(*,*)'add_ray_point cnpar,cnper,cnt,cnn',
c     +           cnpar,cnper,cnt,cnn
cTest 
c      cnray=cnn
cendtest
 
      wnray(n_midway)=cnray
    
c-----calculate electric field (cex,cey,cez) and flux cflown  
    
      call elf_emis(t,u,deru,cex,cey,cez,cflown)
   
c         write(*,*)'emis xe,ye,T_kev',xe,ye,T_kev
c         write(*,*)'cnpar,cnper,cnray',
c     +   cnpar,cnper,cnray
c         write(*,*)'n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin',
c     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin
c         write(*,*)'r,z,phi',r,z,phi
c         write(*,*)'cex,cey,cez',cex,cey,cez
c         write(*,*)'cflown,vgroup,frqncy',cflown,vgroup,frqncy

c-----the data for collisional absorption
      iabsorp_collisional_l=iabsorp_collisional
      coll_mult_l=coll_mult
      tempe_l=tempe(z,r,phi,1) ! kev
      dense_l=dense(z,r,phi,1) 
cSAP110316
c      zeff_l=zeff(z,r,phi,1) 
      zeff_l=zeff(z,r,phi) 
ctest
c      if (itest.eq.2) then
cSm060314
c         call calc_emis_coef(xe,-ye,T_kev,cnpar,cnper,cnray,
c         call calc_emis_coef(xe,ye,T_kev,cnpar,cnper,cnray,
c     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,
c     +   phi,
c     +   cex,cey,cez,cflown,vgroup,frqncy,
c     +   iabsorp_collisional_l,coll_mult_l,tempe_l,dense_l,zeff_l,
c     +   tal,tj,temp_rad)
c         write(*,*)'in midway (n_midawy-1),tal,tj,temp_rad',
c     +   (n_midway-1),tal,tj,temp_rad
c         goto 200
c      endif
    
c      write(*,*)'xe,-ye,T_kev,cnpar,cnper,cnray,
c     +n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,
c     +phi,
c     +cex,cey,cez,cflown,vgroup,frqncy',
c     +xe,-ye,T_kev,cnpar,cnper,cnray,
c     +n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,
c     +phi,
c     +cex,cey,cez,cflown,vgroup,frqncy

cSm060314
c      write(*,*)'emssion in  add_ray_point before calc_emis_coef'
  
      call calc_emis_coef(xe,-ye,T_kev,cnpar,cnper,cnray,
c      call calc_emis_coef(xe,ye,T_kev,cnpar,cnper,cnray,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +i_resonance_curve_integration_method,epsi,
     +i_fkin,r,z,
     +phi,
     +cex,cey,cez,cflown,vgroup,frqncy,
     +iabsorp_collisional,coll_mult,tempe_l,dense_l,zeff_l,
     +wal_emis(n_midway),wj_emis(n_midway),wtemp_rad_em(n_midway))
       
c      write(*,*)'emssion in  add_ray_point after calc_emis_coef'

c      write(*,*)'in midway wj_emis(n_m-1),wj_emis(n_m),wj_emis(n_m+1)',
c     +n_midway,wj_emis(n_midway-1),wj_emis(n_midway),wj_emis(n_midway+1)

ctest           
c      write(*,*)'! n_midway wr_em(n_midway),wr_em(n_midway+1)',
c     +n_midway,wr_em(n_midway),wr_em(n_midway+1)

c      if (itest.eq.1) then
c         n_midway=n_midway+1
c         itest=2
c         goto 100
c      endif
c 200  continue
cendtest

c      do k=1,nrayelt_emis+1
c        write(*,*)'add_ray_point !!!! with new point k,wsn(k),wr_em(k)',
c     +  k,wsn(k),wr_em(k)
c      enddo 

      return
      end  


      subroutine elf_emis(t,u,deru,cex_em,cey_em,cez_em,cflown)
c---------------------------------------------------------------
c     calculate the electric field: cex_em,cey_em,cez_em 
c     and the flux: cflown 
c     for the emission subroutines
c---------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'write.i'
      include 'cefield.i'
      include 'cone.i'
      include 'grill.i'
      include 'eps.i'

      dimension u(*),deru(*),vgr(3),bf(3)
      dimension tempiar(nbulka)
      double complex cnx
      double complex dhot,dhot_rlt,dhot_sum
      double precision x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     .tpop_ar(nbulka),vflow_ar(nbulka)

      double complex fr_func_noeps,roots(2),dispfun
      integer info(2),ier,nsig

!      external dhot
CENM For the iabsorp=12 damping option to work (find full complex solution of
C    D(nperp)=0), the muller root-finding subroutine is called, and either
C    Ram's relativistic function is used (for id=14) or Nelson-Melby's (for id=11)
      external fr_func_noeps  ! for dispersion function root finder (id=14)     
      external dispfun  ! for dispersion function root finder (id=11)

c      integer nbulkaa
cSm030226  
c      parameter (nbulkaa=5)
c      parameter (nbulkaa=nbulka)

      double complex K(3,3),dK(3,3,7),dd(5),d,ddn,aK(3,3)
c      double complex dd5(5,nbulkaa),ddnp_h,ddnll_h,ddnp
      double complex dd5(5,nbulka),ddnp_h,ddnll_h,ddnp
      double complex image
c     for cold plasma+relativistic tensor
      double complex disp_func,dcold_rlt
     
c     for test the flux from grpde2 
cSAP110211
c      complex exde,eyde,ezde
c      real rnpar,rnper,romega,romegpe,rvgrpdc,redenfac
      complex*16 exde,eyde,ezde
      real*8 rnpar,rnper,romega,romegpe,rvgrpdc,redenfac
c-----output
      double complex  cex_em,cey_em,cez_em,cflown  

cSAP120605
c----- For Nicola Bertelli scattering. It was added output argument
c     dwepsdw(3,3) to the subroutine flown (z,r,phi,cnz,cnr,cm,cflown,dwepsdw)
      complex*16 dwepsdw(3,3) 
c--------------------------------------------------



      pi=4.d0*datan(1.d0)

c      if(nbulkaa.lt.nbulk) then
      if(nbulka.lt.nbulk) then
c        write(*,*)'in emission.f nbalkaa.lt.nbulk'
c        write(*,*) 'nbulkaa=',nbulkaa,'nbulk=',nbulk
c        write(*,*)'nbulkaa=nbulka'
c        write(*,*)'change nbulka in param.i'
        write(*,*)'in emission.f nbulka.lt.nbulk'
        write(*,*) 'nbulka=',nbulka,'nbulk=',nbulk
        write(*,*)'change nbulka in param.i'
        stop
      endif
c----------------------------------------
c     cvac (cm/sec)
      cvac=2.997930D+10
c----------------------------------------
c     cld (cm),frgncy(GHz)
      cld=cvac/(2.d0*pi*frqncy*1.0d+09)
c----------------------------------------
c     now proposed that r0x=1 m
      r00=100.d0*r0x
      t00=cvac/(2.d0*pi*frqncy*r00)
c-----------------------------------------
c     z,r (m)
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)
c---------------------------------------
c     bmod,bz,br,bphi (Tl)
      bmod=b(z,r,phi)
      bf(1)=bz
      bf(2)=br
      bf(3)=bphi

      gam=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam)
      dc=dcos(gam)
      cnt=cn(r,cnz,cnr,cm)
      cnpar=cnt*dc
      cnper=cnt*ds
       
c-----electric field calculations----
cSm061127
      if (((iabsorp.eq.3).or.(iabsorp.eq.2))
     &    .or.((iabsorp.eq.91).or.(iabsorp.eq.92))) then
c------------------------------------------------------------
c       electric field using the cold plasma dielectric tensor
        call tensrcld(u(1),u(2),u(3))
        cnx=dcmplx(cnper,0.d0)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
      endif

      if (iabsorp.eq.1) then
   
c-------EC wave case.The complex electric field calculations using
c       hermitian or full mazzucato tensor (with antihermitian part)
c       reps(3,3), and N_per=N_per(N_par) from mazzucato solver
        ihermloc=iherm        
c-------EC absorption (from Mazzucato code)
        ioptmaz=1 ! estimation of cnper1 from coldm
        cnper1=cnper
c        write(*,*)'prep3d  bef cnperm ioptmaz=1 cnper1,cnper',
c     *  cnper1,cnper
c        write(*,*)'perp3d cnper will calculate cnper1,cnprim' 
c        write(*,*)'using the estimation of nperp from cold plasma'
        ihermloc=2
        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)
     
        ioptmaz=2 ! estimation of cnper1 from input parameter cnper
c        cnper1=cnper
c        ihermloc=2
c        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)
c        write(*,*)'ioptmaz=2 cnper1,cnper,cnprim',cnper1,cnper,cnprim

c-------electric field for mazzucato tensor

        ihermloc=iherm
        ihermloc=2 !full hermition + antihermition Mazzucato tens.
        
        call hamiltmuz_ful(cnpar,ihermloc,z,r,phi,cnper1,
     .  cnprim,hamiltmz)

        cnx=dcmplx(cnper1,cnprim)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)

      endif

      goto30 
      if (iabsorp.eq.1) then
c       EC wave case.The complex electric field calculations
c       using hermition or full mazzucato tensor (with antihermition part)
c       reps(3,3), and N_per=N_per(N_par) from mazzucato solver
        ihermloc=iherm
c-------EC absorption (from Mazzucato code)
        ioptmaz=1 ! estimation of cnper1 from coldm
        cnper1=cnper
        ihermloc=2
        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)
    
c        ioptmaz=2 ! estimation of cnper1 from input parameter cnper
c        cnper1=cnper
c        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)

c-------electric field for mazzucato
        ihermloc=iherm
        ihermloc=2 !full hermition + antihermition Mazzucato tens.

        call hamiltmuz_ful(cnpar,ihermloc,z,r,phi,cnper1,
     .  cnprim,hamiltmz) 
 
        cnx=dcmplx(cnper1,cnprim)        
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)

c       write(*,*)'in emission.f efild mazz ex,ey,ez,cnper1,cnpar',
c     1	ex,ey,ez,cnper1,cnpar
      endif

 30   continue
      goto 40  
      if (iabsorp.eq.1) then   
c-------EC wave case.The complex electric field calculations using
c       hermitian or full mazzucato tensor (with antihermitian part)
c       reps(3,3), and N_per=N_per(N_par) from mazzucato solver
        ihermloc=iherm
c-------EC absorption (from Mazzucato code)
        ioptmaz=1 ! estimation of cnper1 from coldm
        cnper1=cnper
c        write(*,*)'perp3d cnper will calculate cnper1,cnprim' 
c        write(*,*)'using the estimation of nperp from cold plasma'
        ihermloc=2
        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)
        write(*,*)'ioptmaz=1 ihermloc,new cnper1,old cnper,new cnprim'
     *, ihermloc,cnper1,cnper,cnprim

        ioptmaz=2 ! estimation of cnper1 from input parameter cnper
c        cnper1=cnper
c        ihermloc=2
c        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)
c        write(*,*)'ioptmaz=2 cnper1,cnper,cnprim',cnper1,cnper,cnprim

        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0

c-------electric field for mazzucato tensor
        ihermloc=iherm
        ihermloc=2 !full hermition + antihermition Mazzucato tens.        
        call hamiltmuz_ful(cnpar,ihermloc,z,r,phi,cnper1,
     .  cnprim,hamiltmz)

        cnx=dcmplx(cnper1,cnprim)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
c       write(*,*)'in prep3d efield mazz ex,ey,ez,enp,cnper1,cnpar',
c     1	ex,ey,ez,enp,cnper1,cnpar
      endif !if(iabsorp.eq.1)
 40   continue

      if(iabsorp.eq.4) then
c       EC wave case.The complex electric field calculations
c       using Forest tensor (with antihermition part)
c       reps(3,3), and N_per=N_per(N_par) from forest solver

        cnparp=cnpar
        cnperp=cnper
        
        do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpop_zrp(z,r,phi,i) !YuP[2024-08-14] was tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo

c-------electric field for Forest tensor
        cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)

        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .  vflow_ar,cnparp,cnperp,1,reps)

        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
       
c       write(*,*)'in prep3d efild forest  ex,ey,ez,cnper,cnpar',
c     1	ex,ey,ez,cnper,cnpar

      endif ! iabsorp.eq.4

c------------------------------------------------------------------
      if(iabsorp.eq.5) then
c-------EC wave case.
c       EC absorption (from Shkarofsky code)
c       The complex dielectric field calculations
c       using Shkarofsky tensor 
c----------------------------------------
	xe=x(z,r,phi,1)
        ye=y(z,r,phi,1)
        te=tempe(z,r,phi,1)
        friq=frqncy
        call absorpsh(cnpar,cnper,xe,ye,te_kev,friq,cnprim)

c-------electric field for Shkarofsky tensor
        call efiedshk (cnpar,cnper,xe,ye,te_kev,friq,cex,cey,cez,
     *  ex,ey,ez)      
      
c       write(*,*)'in prep3d efild sharofsky ex,ey,ez,cnper1,cnpar',
c     1	ex,ey,ez,cnper1,cnpar

      endif ! iabsorp.eq.5

      if(iabsorp.eq.6) then
c       EC wave case.The complex electric field calculations
c       using Forest tensor (with antihermition part)
c       reps(3,3), and N_per=N_per(N_par) from forest solver
c-------EC relativistic absorption 
c       dielectric tensor=hermitian part(from Forest code)+
c                         anti-hermitian part(full relativistic) 
	
        cnparp=cnpar
        cnperp=cnper
        
c-------Hermitian non-relativistic tensor reps        

        do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpop_zrp(z,r,phi,i) !YuP[2024-08-14] was tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo
 
cSm060315
        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .  vflow_ar,cnparp,cnperp,1,reps)
        
c-------anti-hermitian relativistic tensor aK   
        
        if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c----------usage of the analytical relativistic function and its derivatives 
           i_fkin=0
        else 
c----------usage of the mesh relativistic function and its derivatives
           i_fkin=1
        endif 
cSm060315                         
        call anth_rlt(x_ar(1),y_ar(1),t_av_ar(1)*1.d-3,cnparp,cnperp,
     +  n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +  i_resonance_curve_integration_method,epsi,
     +  i_fkin,r,z,phi,
     +  aK)
        image=dcmplx(0.d0,1.d0)
        do i=1,3
        do j=1,3
         reps(i,j)=reps(i,j)+image*aK(i,j)
        enddo
        enddo 
 
c        cnprim=0.d0 
c        d=dhot_rlt(reps,aK,cnparp,cnperp,cnprim)
c        dham=dreal(d)
      
c	call Ddhot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
c     . ,cnparp,cnperp,reps,dd5,ddnp_h,ddnll_h,ddnp)
        
c-------electric field for Forest tensor
        cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)
cSm060315
c        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c     .  vflow_ar,cnparp,cnperp,1,reps)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        
c       write(*,*)'in prep3d efild forest  ex,ey,ez,cnper,cnpar',
c     1	ex,ey,ez,cnper,cnpar

      endif ! iabsorp.eq.6

      if(iabsorp.eq.7) then
c       EC wave case.The complex electric field calculations
c       using Cold plasma tensor +antihermition relativistic tensor
c-------EC relativistic absorption 
c       dielectric tensor=hermitian part(cold plasma)+
c                         anti-hermitian part(full relativistic) 
	
        cnparp=cnpar
        cnperp=cnper
        
c-------calculate Hermitian cold plasma complex tensor reps. 
c       It will be in eps.i        

        call  tensrcld(z,r,phi)
    
        do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
cSm060314
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpop_zrp(z,r,phi,i) !YuP[2024-08-14] was tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo
         
c-------anti-hermitian relativistic tensor aK   
c      n_relt_harm1  - min number of used cyclotron jn harmonics,
c      n_relt_harm2  - max number of used cyclotron jn harmonics,
c                      n_relt_harm1 <= jn <= n_relt_harm2 
c                      jn is number of EC harmonics in relativistic tensor aK

c       n_relt_intgr is the number of points for the numrerical integration
c       over p_perp for the calculations of anti-hermitian
c       relativistic tensor aK.
        
        if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c----------usage of the analytical relativistic function and its derivatives 
           i_fkin=0
        else 
c----------usage of the numerical relativistic function and its derivatives
           i_fkin=1
        endif

c        write(*,*)'elf_emis before anth_rlt x_ar(1),y_ar(1)',
c     &  x_ar(1),y_ar(1)
c        write(*,*)'t_av_ar(1)*1.d-3,cnparp,cnperp',
c     &  t_av_ar(1)*1.d-3,cnparp,cnperp
c        write(*,*)'n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin',
c     &  n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin
c        write(*,*)'r,z,phi',r,z,phi

        call anth_rlt(x_ar(1),y_ar(1),t_av_ar(1)*1.d-3,cnparp,cnperp,
     +  n_relt_harm1,n_relt_harm2,n_relt_intgr,
     &  i_resonance_curve_integration_method,epsi,
     +  i_fkin,r,z,phi,
     +  aK)

c        write(*,*)'elf_emis after anth_rlt aK',aK

c------ complex dispersion function calculated from the sum of
c       of the cold electron plasma dielectric tensor eps_h
c       and the relativistic electron anti-hermition dielectric tensor eps_a
        disp_func=dcold_rlt(reps,aK,cnparp,cnperp)
       
c-------calculate the derivative d(D_hermitian)/d(ReN_perp)
c       from the electron cold plasma dispersion function D
        ddnp=dDcold(reps,cnpar,cnper)
        
        cnprim = dabs(DIMAG(disp_func) / DREAL(ddnp))
	
        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0
     
c-------electric field for cold plasma
        cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)
        
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
        
      endif ! iabsorp.eq.7

      if(iabsorp.eq.8) then
c-------The absorption is calculated for Westerhof- Tokman dispersion
c       It works for id =10,12,13,15 cases.
        cnprim_in=0.d0
        call Im_nperp_Westerhof_Tokman(z,r,phi,
     &  cnpar,cnper,cnprim_in,id,cnprim)

        write(*,*)'after Im_nperp_Westerhof_Tokman cnper,cnprim',
     &  cnper,cnprim

        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0

c-------electric field using the dielectric tensor acording id=10,12 or 13
        cnx=dcmplx(cnper,cnprim)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c       electric field parallel to wave
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
        
      endif !iabsorp=8	

      if(iabsorp.eq.9) then
c--------------------------------------------------------------
c        The absorption is calculated for hot dispersion
c        using the formula from Stix book p.74 (17,18,21)
c        Im(k_perp)= 0.5*Power_abs/(P^+T^) 
c
c        Here 
c    
c        Power_abs=omega/(8pi)[ E~(i) . (eps_a_herm(i,j) . E(j)]
c
c        P^ = (c/16pi)[E~*B+E*B~] is Poining vector,calculated
c             using hot plasma complex dieletric tensor.
c
c        T^ = -omega/(16pi)[ E~(i) . d/dk^(eps_herm(i,j) . E(j)]
c             Is a flux of nonelectromagnetic energy
c----------------------------------------------------------------
         call absorp_hot(z,r,phi,cnpar,cnper,nbulk,
     &   cnprim_e,cnprim_i,cnprim_s)
     
         cnprim=cnprim_e+cnprim_i
         cnprim_cl=0.d0
        
         write(*,*)'after absorp_hot cnprim_e,cnprim_i,cnprim_s',
     &   cnprim_e,cnprim_i,cnprim_s
c----------------------------------------------------------------
c        electric field polarization (cex,cey,cez)
c        calculted  using the hot plasma dielectric tensor
c        will be common /efield.i/
c---------------------------------------------------------------
c        hot plasma (full, non-hermitian) dielectric tensor reps(3,3)
c        will be in common /eps.i/
c--------------------------------------------------------------- 
         cnx=dcmplx(cnper,cnprim)
         call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
      endif !iabsorp=9	

c----------------------------------------------------------------
      if(iabsorp.eq.10) then
c--------------------------------------------------------------
c        The absorption is calculated for relativistic dispersion
c        (combined E. Nelson-Melby  and  A.Ram)
c        using the formula from Stix book p.74 (17,18,21)
c        Im(k_perp)= 0.5*Power_abs/(P^+T^) 
c
c        Here 
c    
c        Power_abs=omega/(8pi)[ E~(i) . (eps_a_herm(i,j) . E(j)]
c
c        P^ = (c/16pi)[E~*B+E*B~] is Poining vector,calculated
c             using hot plasma complex dieletric tensor.
c
c        T^ = -omega/(16pi)[ E~(i) . d/dk^(eps_herm(i,j) . E(j)]
c             Is a flux of nonelectromagnetic energy
c----------------------------------------------------------------
         call absorp_relativist_disp_combined(z,r,phi,cnpar,cnper,
     &   cnprim_e)
 
         cnprim_i=0.d0     
         cnprim=cnprim_e+cnprim_i
         cnprim_cl=0.d0
        
         write(*,*)'after absorp_relativist_disp_combined'
         write(*,*)'cnprim_e',cnprim_e

c----------------------------------------------------------------
c        electric field polarization (cex,cey,cez)
c        calculated using the relativistic dielectric tensor
c        will be common /efield.i/
c---------------------------------------------------------------
c        relativistic (full, non-hermitian) dielectric tensor reps(3,3)
c        will be in common /eps.i/
c---------------------------------------------------------------        
        cnx=dcmplx(cnper,cnprim)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
         
      endif !iabsorp=10	

      if(iabsorp.eq.11) then
c       EC wave case.The complex electric field calculations
c       using relativisic tensor (Ram or Nelson-Melby) (with antihermition part)
c       reps(3,3), and N_per=N_per(N_par) from projection 
c       method for D=Determinant with Ram or Nelson-Melby (depending on id)
c       tensor.


        cnparp=cnpar
        cnperp=cnper
                
        x_e=x(z,r,phi,1)
	y_e=y(z,r,phi,1)
	t_e=tempe(z,r,phi,1) ! kev
        call relativist_absorp_projection_method_det
     &  (t_e,cnpar,x_e,y_e,cnper,reps,D,cnprim)

c-------electric field for Ram tensor
        cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)

        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
       write(*,*)'reps',reps
       write(*,*)'in emission efield Ram  ex,ey,ez,cnper,cnpar',
     1	ex,ey,ez,cnper,cnpar
       

      endif ! iabsorp.eq.11
C---------------------------------------------
      if (iabsorp.eq.12) then
C------------------------------------------
C       Test basic elements needed for iabsorp

C----- use Muller algorithm to find dispersion relation root given
C----- n_parallel and frequency and plasma parameters.
         
C************ Here are some hard-coded numbers for the accuracy with
C************ which to search for the root for the iabsorp.eq.12 option
C************ Since it should be sticking very close, usually the iterations
C************ won't be much of an issue.
         errabs=1.d-6
         nsig=6
         nknown=0
         nrts=1                 ! just look for one root
         nguess=nrts 
         nnew=nrts
         itmax=50
C******* For initial guess, just search with the real part from before and just
C******* 0 imaginary part. Usually, if the ray is propagating, the damping
C******* is low anyway, so it shouldn't be too large imaginary part.
         roots(1)=dcmplx(cnper,0.0d0)
      write(*,*)'()()(()()(()()initial data'
      write(*,*)'z=',z,'r=',r,'phi=',phi
      write(*,*)'cnz=',cnz,'cnr=',cnr,'cm=',cm
      write(*,*)'rho=',rho
c      print *,'====== is=',is
      print *,'cnprim: ',cnprim,' cnpar: ',cnpar,' cnper: ',cnper
c      cn2=cnz**2+cnr**2+cnphi**2
      cn2=cnz**2+cnr**2+(cm/r)**2
      print *,'cn2=',cn2,' cnper(calculated) ',dsqrt(cn2-cnpar**2)
c         if (is.lt.1) then
c            print *,'****** is=',is
c         else
c            print *,'nper before: ',wnper(is-1)
c         endif
C id.eq.11 or id.eq.12, call root-finding using Nelson-Melby dielectric function
         if (id.eq.11 .or. id.eq.12) then
           call muller(dispfun,errabs,nsig,nknown,nguess,nnew,roots,
     +     itmax,info,ier)
         else
C id.eq.14 or 15, call fr_func_noeps, using Ram's dielectric function
           call muller(fr_func_noeps,errabs,nsig,nknown,nguess,nnew,
     +     roots,itmax,info,ier)
         endif
         cnprim=abs(imag(roots(1)))
         cnper=abs(dble(roots(1)))

         cnprim_cl=0.d0
         cnprim_e=cnprim
         cnprim_i=0.d0

C******* To be consistent with all other methods of calculating
C******* force cnprim and cnper to be positive.
         print *,'&&&&&&&&&&&& cnprim: ',cnprim,' cnper: ',cnper

         cnx=dcmplx(cnper,cnprim)
         call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
      endif   !iabsorp=12

C----------------------------- END OF IABSORP MODULES -----------------

c--------------------------------------------------------
c     cflown - dimensionless for E_x/E,E_y/E,E_z/E
c  !!!!now flown is calculated using Mazzucato dielectric tensor
c  !!!!and electric field was calculated using Mazzucato tensor
c  !!!!it is only for EC wave .For LH and FW it is necessery
c !!!!!to create new subroutine flown ?what tensor?
c      call flown (u(1),u(2),u(3),u(4),u(5),u(6),cflown,dwepsdw)
c       id_old=id
c      if (iabsorp.eq.4) then
c         id=6 !hot plasma dispersion
c         call flown (u(1),u(2),u(3),u(4),u(5),u(6),cflown,dwepsdw)
c         id=id_old
c      endif

c      if (iabsorp.eq. 1) then
c         id=4
c         call flown (u(1),u(2),u(3),u(4),u(5),u(6),cflown,dwepsdw)
c         id=id_old
c      else !uses the dispersion with the given 'id' 
c         call flown (u(1),u(2),u(3),u(4),u(5),u(6),
c     1                   cflown,dwepsdw)
c      endif

c-----flux from the cold plasma    
c      xe=x(z,r,phi,1)
c      ye=y(z,r,phi,1)
c      rnpar=cnpar
c      rnper=cnper
c      romega=1/ye
c      romegpe=dsqrt(xe)/ye
c      write(*,*)'emission.f prep3d xe,ye',xe,ye
c      nsigma=ioxm
c      nsigma=1
c      nsigma=-1
c      call grpde2 (rnpar, rnper, romega, romegpe, nsigma,
c     .                   rvgrpdc, redenfac, exde, eyde, ezde)
c      write(*,*)'+ redenfac',redenfac
c      cflown=2.d0*redenfac
 
      if (iflux.eq.1) then
cSAP120605
          call flown (u(1),u(2),u(3),u(4),u(5),u(6),cflown,dwepsdw)
      endif

      if (iflux.eq.2) then
c------- flux from the cold plasma    
         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
         rnpar=cnpar
         rnper=cnper
         romega=1/ye
         romegpe=dsqrt(xe)/ye
         nsigma=ioxm
c        nsigma=1
c        nsigma=-1
         call grpde2 (rnpar, rnper, romega, romegpe, nsigma,
     .                   rvgrpdc, redenfac, exde, eyde, ezde)

c        write(*,*)'+ redenfac',redenfac
         cflown=2.d0*redenfac
c         write(*,*)'prep3d cold from grpde2 cflown',cflown    
      endif !cold electron plasma flux
    
c-------------------------
      cex_em=cex
      cey_em=cey
      cez_em=cez
    
      return
      END    


      subroutine intinit(omega,n,npar,cdvt,vmax,v0,vpar0,vpar1,
     & vpar2,thet1,thet2,ires)
      implicit double precision (a-h,o-z)
c
c  For integration we require that at least one of the vpar1,vpar2
c  (vperp=0. ) points of the resonance curve lie
c  within +/- vmax.
c  Integration is carried out in the clockwise direction along only
c  that portion of the resonance curve which
c  lies inside vmax and includes that
c  vpar1 or vpar2 which has the smallest absolute value laying
c  inside vmax.
c     Normalized velocities are used.
c-----input omega=emega/omegac=1/y
c           n    is the number of EC harmonic
c           npar is the parallel refructive index
c           cdvt -is the velocity of the normalization (cm/sec)
c           vmax  is max

      double precision npar

c-----output 
c     resonace ellipse N_par**2<1
c     (vpar-vpar0)**2/v0**2+vperp**2/vmax*2=1
c     (v0/c)**2=((n*Y)**2-(1-npar**2))/(1-npar**2)**2
c     vpar0/c=Y*npar/(1=npar**2)    
c     vmax=v0*dsqrt(1-npar**2)
c     vmax is the max V_per at the ellipse
c     v0   is the 
c     vper1=vpar0-dsign(vpar0)*v0
c     vpar2=vpar0+dsign(vpar0)*v0
c
c     ires indicates the subroutine result:
c     ires= 1,  OK.  v0,vpar0,vpar1,vpar2,thet1,thet2 are found.
c           2,  npar**2.ge.1,  no results
c           3,  imaginary resonance curve.
c           4,  integration path not on grid.
c
c     Determine intersections of resonance curve with v=vmax and
c     work out the limits on the angles in the integration
c     theta1,theta 
c  
c---- local
      double precision npar2   !N_par**2
      double precision npar2m1 !N_par**1-1

c 
      pi=4.d0*datan(1.d0) !YuP[2018-10-13][2020-01-27] BUG was: pi=4.d0*dtan(1.d0)

      ires=1
      vpar1=0.d0
      vpar2=0.d0
      cdvt2=cdvt*cdvt
      npar2=npar*npar
      npar2m1=npar2-1.d0
c
      if(npar2m1.ge.0.d0)  ires=2 ! hyperbole case
      if(ires.ne.1)  return
c
      v02=(((n/omega)**2+npar2m1)/npar2m1**2)*cdvt*cdvt
      if(v02.le.0d0)  ires=3
c  If no real resonance
      if(ires.ne.1)  return
c
      v0=dsqrt(v02)
      vpar0=-n/omega*npar/npar2m1*cdvt
      s=dsign(1.d0,vpar0)
      vpar1=vpar0-s*v0
      vpar2=vpar0+s*v0
c
c  Discover relationship of resonance ellipse to v=vmax.
      if(dabs(vpar1).ge.vmax)  ires=4
c  If resonance curve at vperp=0. doesn't lie on grid:
      if(ires.ne.1)  return
c
c  Determine intersections of resonance curve with v=vmax and
c  work out the limits on the angles in the integration
c
      thet1=0.d0
      thet2=pi
      if(dabs(vpar2).lt.vmax) return
      gammax=dsqrt(1.d0+vmax*vmax/(cdvt*cdvt))
      vparr=cdvt*(gammax-n/omega)/npar
      vperr=dsqrt(vmax**2-vparr**2)
c
      if(s.gt.0.d0) then
        thet1=0.d0
        thet2=datan2(vperr/dsqrt(-npar2m1),vpar0-vparr)
      else
        thet2=pi
        thet1=datan2(vperr/dsqrt(-npar2m1),vpar0-vparr)
      endif
c
      if(thet1.lt.0.d0) stop 111
      if(thet1.gt.pi) stop 112
      if(thet2.lt.0.d0) stop 113
      if(thet2.gt.pi) stop 114
      if(thet2.le.thet1) stop 114
c
      return
      end






      subroutine unde_intgrh(vperdvt,vpardvt,cdvt,y,kpar,
     +theta,n,i_fkin,r,z,phi,u_n,u_flcur_n)
    

c     calculates under integral functions
c     u_n=U_n(np=p_per,p_perp=p_parallel)=(1/gamma)*
c     (n*Y*df/dp^_perpendicular+N_parallel*p^_perpendicular*df/dp^_parallel)
c
c     u_n=n/dabs(kpar)*df/dvperdvt+dsign(kpar)*vperdvt*df/dvpardvt
c     kpar =K_par*Vt/omegac
c     u_n is used in under integral complex marix function G_nk(3,3)
c
c     u_flcur_n=f*p^perp/gamma
c     u_flcur_n is used in under integral complex matrix fluct_car_nk(3,3)
c     Here p^=p/mc
c
c     input
c       vperdvt  =perpendicular momentum divided by (mVt)
c       vpardvt  =parallel momentum/mVt
c       cdvt     =c/Vt
c       y       = omega_ce/omega for the electron rest mass    
c       kpar    = =K_par*Vt/omegac
c       theta   - mc**2/T
c       n       - n number of the given resonance harmonic
c       i_fkin =0 the usage of the analytical relativistic Maxwellian distributin
c              =1 the usage of the numerical 3D distribution from diskf or
c                 netcdfnm.nc files or the mesh distribution given analytically
c                 written be CQL3D code
c       r       - the major radius (it is used for i_fkin=1)
c       z       - the vertical coordinate  (it is used for i_fkin=1)
c       phi     - toroidal angle
c-------------------------------------------------------
      implicit none
c     input        
      double precision vperdvt,vpardvt,cdvt,y,kpar,theta
      double precision r,z,phi
      integer n,i_fkin

c     external besk2as, root_res,psif,rhopsi,bmin_psi,b
      double precision psif,rhopsi,bmin_psi,b,fdens_fdist  

c     output
      double precision u_n,u_flcur_n

c_for test
      double precision vnormloc,massloc
      COMMON /dskin1/vnormloc,massloc  !v (cm/sec),m( g)
cendtest

c     local
      
      double precision gamma, bk2,f_maxw,ps,p,ps_max,eps,pi,
     + m_e,clight,energy,psi,rho,pitch,pitch0,bmin,btotal,
     + c_pitch,s_pitch,s_pitch0,c_pitch0,
     + dptc0dpt,                          !d(pitch0)/d(pitch)
     + dens,u_np,
     + p_perp,p_par,  !momentum components/mc
     + nll,            !n_par
     + dfdvpert,dfdvpart !df/d(v_per/Vt)
c     distributin function and its derivatives
      double precision fdist0,dfdx,dfdpitch,dfdpitc0,dfdp
      integer initial     
      
      pi=4*datan(1.d0)
      p_perp=vperdvt/cdvt
      p_par=vpardvt/cdvt
      ps=p_perp*p_perp+p_par*p_par
      p=dsqrt(ps)
      gamma=dsqrt(1.d0+ps)

cSm050927
      nll=kpar*y*cdvt

c      write(*,*)'unde_intgrlh i_fkin,theta,gamma,',i_fkin,theta,gamma
      if (i_fkin.eq.0) then
c--------usage of the analytical relativistic Maxwellian distribution

c        calculation the Mackdonalds  function bk2=K_2(theta)*EXP(theta)	
         call besk2as(theta,bk2)              

c--------the following function is normalized as
c        integral{f_maxw*(v/c)**2*d(v/c)*2*pi}
           
         f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)

c--------the following function is normalized as
c        integral{f_maxw*(v/Vt)**2*d(v/Vt)*2*pi}
         f_maxw=f_maxw/(cdvt)**3 
                        
c         u_n=n/dabs(kpar)*df/d(vperdvt)+dsign(kpar)*vperdvt*df/d(vpardvt)
c         kpar =K_par*Vt/omegac
c         df/d(vperdvt)=df/d(gamma)*d(gamma)/d(vperdvt)
c         d(gamma)/d(vperdvt)=vperdvt/gamma/(cdvt)**2
c         d(gamma)/d(vpardvt)=vpardvt/gamma/(cdvt)**2
c         df/d(gamma)=-f_maxw*theta

c          u_n=-p_perp*theta*f_maxw*(n*y+nll*p_par)/(gamma*gamma)
          
c          u_n=-vperdvt*theta*f_maxw*
c     +        (n/dabs(kpar)+dsign(1.d0,kpar)*vperdvt)! like in Horace

          dfdvpart=-theta*f_maxw*vpardvt/(cdvt**2*gamma)
          dfdvpert=-theta*f_maxw*vperdvt/(cdvt**2*gamma)
          u_n=n/dabs(kpar)*dfdvpert+dsign(1.d0,kpar)*vperdvt*dfdvpart
          write(*,*)'kpar,dfdvpert,dfdvpart,vperdvt,f_maxw,u_n',
     +    kpar,dfdvpert,dfdvpart,vperdvt,f_maxw,u_n
 
c          u_flcur_n=p_perp*f_maxw/gamma
c           u_flcur_n=vperdvt*f_maxw/gamma
          u_flcur_n=f_maxw

           write(*,*)'u_nh vperdvt,n,y,kpar,f_maxw,u_n,u_flcur_n',
     +     vperdvt,n,y,kpar,f_maxw,u_n,u_flcur_n

         goto10 
      endif

      if (i_fkin.eq.1)then
c--------usage of the numerical 3D relativistic distribution from diskf file
c        written by CQL3D code
         initial=0
c        energy=p**2/(2m) ps=(p/mc)**2. energy should be in KeV
c        pitch
c        rho
         m_e=9.1094d-28                !g
         clight=2.99792458d10          !cm/sec
         energy=0.5d0*ps*m_e*clight**2 !erg 
         energy=energy/1.6022d-9         !KeV
 
         psi=psif(z,r)                 !poloidal flux
         rho=rhopsi(psi)!small radius

         if (p.ne.0.d0) then
             pitch=dacos(p_par/p)
         else
             pitch=0.d0
         endif
         if (pitch.gt.pi) pitch=pi
         if (pitch.lt.-pi) pitch=-pi
          

c------- pitch0 is pitch angle for the point with the minimal
c------- value of the magnetic field at the give flux surface with the poloidal flux=psi  
         pitch0=pitch
         dptc0dpt=1.d0
         goto 100

         bmin=bmin_psi(psi)
         btotal=b(z,r,phi)
         
         if (p.ne.0.d0) then
             c_pitch=p_par/p     !cos(pitch)
             s_pitch=p_perp/p    !sin(pitch)
         else
             c_pitch=1.d0     !cos(pitch)
             s_pitch=0.d0     !sin(pitch)
         endif
         
c--------s_pitch0=sin and c_pitch0=cos of the pitch angle for the minimal b
         s_pitch0=s_pitch*dsqrt(bmin/btotal) !sin(pitch0) 
c         write(*,*)'bmin,btotal,s_pitch,s_pitch0',
c     +   bmin,btotal,s_pitch,s_pitch0

         if (s_pitch0.gt.1.d0) then
            s_pitch0=1.d0
            c_pitch0=0.d0
         else    
            c_pitch0=dsqrt(1.d0-s_pitch0**2) !dabs(cos(pitch0))
         endif
   
         if (c_pitch.lt.0.d0) then
            c_pitch0=-c_pitch0
         endif

         pitch0=dacos(c_pitch0)
         if (pitch0.gt.pi) pitch0=pi
         if (pitch0.lt.-pi) pitch0=-pi
         
c         write(*,*)'emission in u_n pitch,pitch0',pitch,pitch0

c--------d(pitch0)/d(pitch)
         dptc0dpt=(c_pitch/c_pitch0)*(s_pitch0/s_pitch)
c         write(*,*)'c_pitch,c_pitch0',c_pitch,c_pitch0
c         write(*,*)'s_pitch,s_pitch0',s_pitch,s_pitch0

 100     continue

c--------calculations of the distribution function fdist0 
c        and its derivatives df/dx,df/dpitch0
c        x = momentum-per-mass(nomalized to maximum 1.)
c        vel=dsqrt(2.d0*energy*1.d3*1.6022d-12/fmass(1))
c        at the point (rho,energy,pitch0). Here pitch0=pitch0(r,z,pitch)

         call dskin(initial,energy,pitch0,rho,fdist0,dfdx,dfdpitc0,
     +              dfdp,2)
         
         dfdpitch=dfdpitc0*dptc0dpt
c         write(*,*)'dfdpitc0,dptc0dpt',dfdpitc0,dptc0dpt

         u_n=1.d0/gamma*(
     +   n*y*(dfdp*p_perp/p-dfdpitch*p_par/p**2)+
     +   nll*p_perp*(dfdp*p_par/p+dfdpitch*p_perp/p**2))*
     +   (clight/vnormloc)**3

         u_flcur_n=(p_perp*fdist0/gamma)*(clight/vnormloc)**3
        
ctest_begin(comparison of fdist0 and the analytical maxwellian diatribution and its derivatives
c--------calculation the Mackdonalds function bk2=K_2(theta)*EXP(theta)	
c         call besk2as(theta,bk2)
c         f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)
c         f_maxw=f_maxw*(vnormloc/clight)**3
c         write(*,*)'in u_n_emis vnormloc,theta,gamma,f_maxw,fdist0',
c     +   vnormloc,theta,gamma,f_maxw,fdist0
c        u_np=-p_perp*theta*f_maxw*(n*y+nll*p_par)/(gamma*gamma)
c        write(*,*)'unde_intgr u_np,u_n',u_np,u_n
cend_test
      endif !i_fkin.eq.1

 10   continue

      return
      end    


      subroutine ec_condh(n,Y,nll,vmaxdc,ires,v0dc,vmax1dc,vpar0dc,
     +vpar1dc,vpar2dc,thet1,thet2)
c     controls EC resonance condition gamma=N_par*p_par+nY
c     and calcultes the max value p_perp0 of the perpendicular moment
c     for the ellipse case
c     input 
c      n      - the number of EC harmonic
c      Y      - omega_ce/omega for the electron rest mass
c      nll    - N_parallel to the magnetic field
c      vmaxdc - vmax/dc the max value of the momentum/mc
c     output
c      ires =0 no EC resonance
c           =1 resonance ellipse
c           =2 resonance parabola
c           =3 resonance hyperbole
c      for the elipse case
c      (vpar-vpar0)**2/v0**2+vperp**2/vmax1*2=1
c      (v0/c)**2=((n*Y)**2-(1-npar**2))/(1-npar**2)**2
c      vpar0/c=Y*npar/(1-npar**2)    
c      vmax1=v0*dsqrt(1-npar**2)
c
c      v0dc      can be determine if (v0/c)**2.gt.0 . In opposite  case ires=0 no resonance
c      vmax1dc - is the max perpendicular to the magnetic field momentum 
c                divided by (mc) for the ellipse case;
c                m is the electron rest mass; c is the light speed. 
c      vpar0dc   vpar0/clight is the center of the ellipse
c      vpar1dc   vpar0dc-dsign(vpar0)*v0dc
c      vpar2dc   vpar0dc+dsign(vpar0)*v0dc
c
c      thet1     theta  (radian) is the along the ellipse:
c      thet2     vpar=vpar0-v0*cos(theta),vperp=vmax1*sin(theta)
c                thet1=< theta =<thet2 are the limits of theta in the  integration,
c                determined by the ellipse intersection with v=vmax 
 
      implicit none
c-----input
      integer n
      double precision Y,nll,vmaxdc

c-----output
      integer ires
      double precision v0dc,vmax1dc,vpar0dc,
     +vpar1dc,vpar2dc,thet1,thet2

c-----local
      double precision nlls,nY,nYs
      double precision npar2m1,v0dcs,s,pi,gammax,vparr,vperr,      
     &thet_1,thet_2     

      pi=4.d0*datan(1.d0)
      nY=n*Y 
      nYs=nY*nY     
      nlls=nll*nll
     
      npar2m1=nlls-1.d0
      vpar1dc=0.d0
      vpar2dc=0.d0
      ires=0
     
      if(nlls.lt.1.d0)then
c-------resonace ellipse N_par**2<1
c       v=p/me
c       (vpar-vpar0)**2/v0**2+vperp**2/vmax1*2=1
c       (v0/c)**2=((n*Y)**2-(1-npar**2))/(1-npar**2)**2
c       vpar0/c=Y*npar/(1-npar**2)    
c       vmax1=v0*dsqrt(1-npar**2)
c       vmax1 is the max value of V_perp at the ellipse
c       v0   is the value of |vpar-vpar0|=v0 at vperp=0
c       vpar1=vpar0-dsign(vpar0)*v0
c       vpar2=vpar0+dsign(vpar0)*v0
c-------------------------------------------------------
        v0dcs=(nYs-(1.d0-nlls))/(1.d0-nlls)**2 !(v0/clight)**2
        
        if (v0dcs.gt.0.d0)then
           ires=1
           v0dc=dsqrt(v0dcs)
           vmax1dc=v0dc*dsqrt(1.d0-nlls)!max value of vperp/clight
           vpar0dc=nY*nll/(1.d0-nlls)   !vpar0/clight is the center of the ellipse
           s=dsign(1.d0,vpar0dc)
           vpar1dc=vpar0dc-s*v0dc       !the values of v_par at v_per=0
           vpar2dc=vpar0dc+s*v0dc       ! 
c
c          Determine intersections of resonance curve with v=vmax and
c          work out the limits on the angles in the integration
c
           thet1=0.d0
           thet2=pi 
          
            if(dabs(vpar2dc).lt.vmaxdc) return 
           gammax=dsqrt(1.d0+vmaxdc*vmaxdc)
           vparr=(gammax-nY)/nll

cSm060418
           if((vmaxdc**2-vparr**2).ge.0.d0) then
c------------intersection of the ellipse with the circle maxdc
             vperr=dsqrt(vmaxdc**2-vparr**2)
           else
             !no intersection
           endif
c           write(*,*)'ec_condh vmaxdc,vparr,vperr',vmaxdc,vparr,vperr
c           write(*,*)'s',s
           if(s.gt.0.d0) then
             thet1=0.d0
             thet2=datan2(vperr/dsqrt(-npar2m1),vpar0dc-vparr)
           else
c             write(*,*)'ec_condh vperr,vpar0dc,vparr,npar2m1',
c     &                  vperr,vpar0dc,vparr,npar2m1
             thet2=pi
             thet1=datan2(vperr/dsqrt(-npar2m1),vpar0dc-vparr)
c              write(*,*)'ec_condh thet2,thet1',thet2,thet1
           endif
c
           if(thet1.lt.0.d0) then
              write(*,*)'emission.f thet1.lt.0'
              stop
           endif

           if(thet1.gt.pi) then
             write(*,*)'emission.f thet1.gt.pi'
             stop
           endif
             
           if(thet2.lt.0.d0) then
              write(*,*)'emission.f thet2.lt.0'
              stop
           endif

           if(thet2.gt.pi) then
              write(*,*)'emission.f thet2.gt.pi'
              stop
           endif

           if(thet2.le.thet1) then
              write(*,*)'thet2.le.thet1'
              stop 
           endif
c
        else
c----------no resonace
           ires=0
           goto 10
        endif    
                             

c        if(((nlls+nYs).lt.1.d0).or.(nY.le.0.d0)) then
cc         no resonace
c          ires=0
c        else    
c          ires=1   
c          p_perp0=dsqrt((nlls+nYs-1.d0)/(1.d0-nlls))
c        endif
c-------------------------------------------------------
      else  
        if(nlls.gt.1.d0)then  
c         resonance hyperbole
c-------------------------------------------------------
          ires=3
c-------------------------------------------------------
        else 
c         nlls.eq.1.d0
c         resonace parabola 
c-------------------------------------------------------
          if(nY.le.0.d0) then
c           no resonace 
            ires=0
          else
            ires=2 
          endif
c--------------------------------------------------------
        endif 
      endif

 10   continue
      thet_2=thet2
      thet_1=thet1
      thet1=dmin1(thet_1,thet_2)
      thet2=dmax1(thet_1,thet_2)
      return
      end

      subroutine intgr_emish(n,nll,kper,Y,theta,cdvt,vmaxdt,
     +n_relt_intgr,
     +i_fkin,r,z,phi,
     +integral,fluctcur_n,ires,v0dc,vmax1dc,vpar0dc,vpar1dc,vpar2dc)
c-----------------------------------------------------------------
c     calculates the matrix: double complex integrals
c     for the relativistic electron plasma 
c     I(n)^=integral(0<=p_perp<=p_perp0)g_n, g_n={sum_k=1,2,G_nk(p_perp)
c     I(n)^=integral(thet1<thet<thet2)g_n
c     fluctcar(n)^=integral(0<=p_perp<=p_perp0)g_fluctcur_n.
c               g_fluctcur__n={sum_k=1,2,g_fluctcur_nk(p_perp)
c     fluctcar(n)^=integral(thet1<thet<thet2)g_fluctcur_n.
c               g_fluctcur__n={sum_k=1,2,g_fluctcur_nk(p_perp)
c     for the EC harmonic with number 'n'
c----------------------------------------------------------------- 
c     input
c       n        - EC harmonic number
c       nll      = N_par
c       kper     =K_perp*Vt/omegac  
c       Y        = omega_ce/omega for the electron rest mass
c       theta    = mc**2/T
c       cdvt     =clight/Vt
c       vmaxdt   = max velocity on the greed for distribution function  Vmax/Vt
c       n_relt_intgr - the number of points for the integration over p_perp
c       p_perp0  - max value of the perpendicular momentum divided by mc
c                  on the resonanse ellipse
c       i_fkin   =0 the usage of the analytical relativistic distributin
c                =1 the usage of the numerical 3D distribution from diskf 
c                   written be CQL3D code
c       r        - the major radius (it is used for i_fkin=1)
c       z        - the vertical coordinate  (it is used for i_fkin=1)
c       phi        - toroidal angle (it is used for i_fkin=1)
c     output
c       integral(3,3) double complex integral from G_nk over 0<=p_perp=<p_perpmax
c       fluctcur(3,3) double complec integral from fluct_cur_nk over 0<=p_perp=<p_perpmax
c       ires,v0dc,vmax1dc,vpar0dc,vpar1dc,vpar2dc
c       ires     =0 no resonance,=1 ellipse,=2 parabola,=3 hyperbole 
      implicit none
c     input
      double precision nll,kper,Y,theta,vmaxdt,cdvt
      integer n,n_relt_intgr,i_fkin 
      double precision r,z,phi
      double precision vnormloc,massloc ! cm/sec, g
      COMMON /dskin1/vnormloc,massloc

      integer ninta
      parameter(ninta=201)      ! max number of the integration points
      
c     output
      double complex integral(3,3)
      double complex fluctcur_n(3,3)
      double precision v0dc,vmax1dc,vpar0dc,vpar1dc,vpar2dc
      integer ires 

c     local
      double precision kpar,eps,p_permax,h,p_perp,p,p_t,clight
      double complex i,g_n(3,3),g_fluctcur_n(3,3)
      integer j,j1,j2,nint

      double precision dth,thet1,thet2,dvper
      double precision xint(ninta),bes2(ninta),besp2(ninta),cs(ninta),
     &     sn(ninta),vper(ninta),vpar(ninta),gamma(ninta),b(ninta),
     &     ff(ninta),u(ninta),edsde(ninta),denom(ninta)
c     den(ird)
      double complex g_n_ar(3,3,ninta),g_flcrn_ar(3,3,ninta)
      integer ir

c     external
c     ec_cond, zeroK, g_n 
      
      if(n_relt_intgr.gt.ninta) then
         write(*,*)'intgr_emish n_relt_intgr.gt.ininta'
         write(*,*)'reduce n_relt_intgr in genray.in'
         write(*,*)'or increase ninta in intgr_emish'
         stop
      endif
      
      i = ( 0.0d0,1.0d0)        !imaginary number     
      
      kpar=nll/(y*cdvt)         !K_par*Vt/omegac
      write(*,*)'intgr_emish cdvt',cdvt
c-----determine the resonace curves and the angles thet1,thet2 for the
c     integration along the resonce courve angle (thet)
      call ec_condh(n,Y,nll,vmaxdt/cdvt,ires,v0dc,vmax1dc,vpar0dc,
     +vpar1dc,vpar2dc,thet1,thet2)
      write(*,*)'intgr_emis after ec_condh ires,v0dc,vmax1dc,vpar0dc',
     + ires,v0dc,vmax1dc,vpar0dc
      write(*,*)'vpar1dc,vpar2dc,thet1,thet2',
     +vpar1dc,vpar2dc,thet1,thet2
      write(*,*)'1 cdvt',cdvt

      write(*,*)'kpar,kper',kpar,kper

      if (ires.eq.0) goto 10 !no resonace

      if (i_fkin.eq.0) then
c--------i_fkin =0 the usage of the analytical relativistic Maxwellian distributin
c        the accuracy for the min value of Maxwellian exp(-theta(gamma-1))>eps
         eps=1.d-7                
         p_t=10.d0*dsqrt((1.d0/theta+1.d0)**2-1.d0) ! 10*thermal momentum
         p_permax=p_t                    
      endif !i_fkin=0

      write(*,*)'2 cdvt',cdvt 

      if (i_fkin.eq.1) then
c         i_fkin=1 the usage of the numerical 3D distribution from diskf
c         or netcdfnm.nc or the mesh distributiorn given analytically 
          clight=2.99792458d10          !cm/sec
c         
      endif

      write(*,*)'3 cdvt',cdvt 

c     calculations of the integrals over p_perp

      call zeroK(integral)

      write(*,*)'4 cdvt',cdvt 
      call zeroK(fluctcur_n)
      write(*,*)'5 cdvt',cdvt 
c     Subdivide theta range of integration
      dth=(thet2-thet1)/(n_relt_intgr-1)
      write(*,*)'thet1,thet2,n_relt_intgr,dth',
     +thet1,thet2,n_relt_intgr,dth
      write(*,*)'6 cdvt',cdvt
      do 100  ir=1,n_relt_intgr
         xint(ir)=thet1+(ir-1)*dth
         cs(ir)=dcos(xint(ir))
         sn(ir)=dsin(xint(ir))
         vpar(ir)=(vpar0dc*cdvt-v0dc*cdvt*cs(ir)) !vper/Vt
         vper(ir)=vmax1dc*cdvt*sn(ir)             !vpar/Vt
c         b(ir)=kper*vper(ir)
        
c         write(*,*)'100 ir,xint(ir),cs(i),sn(ir)',
c     +   ir,xint(ir),cs(ir),sn(ir)
         
c         write(*,*)'cdvt,vpar0dc*cdvt,v0dc*cdvt,vmax1dc*cdvt,vpar(ir),
c     +   vper(ir)',
c     +   cdvt,vpar0dc*cdvt,v0dc*cdvt,vmax1dc*cdvt,vpar(ir),vper(ir)
            
         gamma(ir)=dsqrt(1.d0+(vpar(ir)**2+vper(ir)**2)/cdvt**2)
         denom(ir)=dabs(gamma(ir)*kpar-vpar(ir)/(cdvt**2*Y))

c         write(*,*)'gamma(ir),kpar,vpar(ir),cdvt2,omega=1/Y,denom(ir)',
c     +   gamma(ir),kpar,vpar(ir),cdvt**2,1/Y,denom(ir)
  
       call g_n_emish(vper(ir),vpar(ir),y,kper,kpar,cdvt,
     .   theta,n,i_fkin,r,z,phi,
     .   g_n,g_fluctcur_n)
c      write(*,*)'after g_n_emish ir g_n'
         do j1=1,3
           do  j2=1,3
              g_n_ar(j1,j2,ir)=g_n(j1,j2)
              g_flcrn_ar(j1,j2,ir)=g_fluctcur_n(j1,j2)
c              write(*,*)'ir,j1,j2,g_n_ar(j1,j2,ir)',
c     +        ir,j1,j2,g_n_ar(j1,j2,ir)
           enddo
         enddo
        
 100  continue
        
c     Integrating
      do j1=1,3
         do j2=1,3
            integral(j1,j2)=0.d0
            fluctcur_n(j1,j2)=0.d0
         enddo
      enddo

      nint=n_relt_intgr
      do 200  ir=1,nint,nint-1
         dvper=0.5d0*dabs(vper(2)-vper(1))
         if(ir.ne.1)  dvper=0.5d0*dabs(vper(nint)-vper(nint-1))
c         e1=e1+dvper*u(ir)*edsde(ir)/denom(ir)
c         g1=g1+dvper*vper(ir)*edsde(ir)*ff(ir)/denom(ir)

         do j1=1,3
            do j2=1,3
               integral(j1,j2)=integral(j1,j2)+
     +         dvper*g_n_ar(j1,j2,ir)/denom(ir)
               fluctcur_n(j1,j2)=fluctcur_n(j1,j2)+
     +         dvper*vper(ir)*g_flcrn_ar(j1,j2,ir)/denom(ir)
            enddo
         enddo
c         write(*,*)'200 ir,g_n_ar(1,1,ir),g_n_ar(1,2,ir)',
c     +    ir,g_n_ar(1,1,ir),g_n_ar(1,2,ir)
c         write(*,*)'200 integral(1,1),integral(1,2)',
c     +   integral(1,1),integral(1,2)

 200  continue
     

      do 201  ir=2,nint-1
         dvper=0.5d0*dabs(vper(ir+1)-vper(ir-1))
c         e1=e1+dvper*u(ir)*edsde(ir)/denom(ir)
c         g1=g1+dvper*vper(ir)*edsde(ir)*ff(ir)/denom(ir)
         do j1=1,3
            do j2=1,3
               integral(j1,j2)=integral(j1,j2)+
     +         dvper*g_n_ar(j1,j2,ir)/denom(ir)
               fluctcur_n(j1,j2)=fluctcur_n(j1,j2)+
     +         dvper*vper(ir)*g_flcrn_ar(j1,j2,ir)/denom(ir)
            enddo
         enddo
         
c         write(*,*)'201 ir,g_n_ar(1,1,ir),g_n_ar(1,2,ir)',
c     +    ir,g_n_ar(1,1,ir),g_n_ar(1,2,ir)
c         write(*,*)'201 integral(1,1),integral(1,2)',
c     +   integral(1,1),integral(1,2)
c         write(*,*)'201 fluctcur_n(1,1),fluctcur_n(1,2)',
c     +   fluctcur_n(1,1),fluctcur_n(1,2)

 201  continue   
     
 10   continue
      return
      end   
  
      subroutine g_n_emish(vperdvt,vpardvt,y,kper,kpar,cdvt,
     .theta,n,i_fkin,
     .r,z,phi,
     .g_n,g_fluctcur_n)
c     calculates under integral complex matrix functions
c     g_n(3,3)=sum{k}g_nk and g_fluctcur_n(3,3)=sum{k}g_fluctcur_nk

c     input
c       vperdvt =parallel momentum divided by (mVt)
c       vpardvt =perpendicular momentum/mVt, Vt is the thermal velocity?
c       y       =omega_ce/omega for the electron rest mass    
c       kpar    =K_par*Vt/omegac  
c       kper    =K_perp*Vt/omegac 
c       cdvt   = clight/Vt
c       theta   = mc**2/T
c       n         n number of the given resonance harmonic
c       i_fkin =0 the usage of the analytical relativistic distributin
c              =1 the usage of the numerical 3D distribution from diskf 
c                 written be CQL3D code
c       r       - the major radius (it is used for i_fkin=1)
c       z       - the vertical coordinate  (it is used for i_fkin=1)
c       phi     - toroidal angle(it is used for i_fkin=1)
c     output
c       g_n(3,3) under integral matrix double complex function
c       g_fluctcur_n(3,3)
c-------------------------------------------------------
      implicit none
c-----input        
      double precision vperdvt,vpardvt,y,kpar,kper,cdvt,theta,r,z,phi
      integer n,i_fkin
c-----external s_calcnh, zeroK, unde_intgrh
      
c-----output
      double complex g_n(3,3), g_fluctcur_n(3,3)
      
c-----local
      double complex sn(3,3)
      double precision v2,cdvt2,gamma,omega,coeff,eps,ps,ps_max,pi,
     +dgdp_par,u_n,u_flcur_n,coeff_fc,
     +denom
      integer i,j

      pi=4.d0*datan(1.d0)
     
c-----initialize g_n
      call zeroK(g_n)
      call zeroK(g_fluctcur_n)
  
c      eps=1.d-9 ! the min value of Maxwell exponent
c      ps_max=(1.d0-dlog(eps)/theta)**2-1.d0
      
      v2=vperdvt**2+vpardvt**2
      cdvt2=cdvt*cdvt
      gamma=dsqrt(1.d0+v2/cdvt2)
      omega=1.d0/y

c      write(*,*)'g_n_emis bef s_calch n,vperdvt,vpardvt,kper',
c     +n,vperdvt,vpardvt,kper

      call s_calcnh(n,vperdvt,vpardvt,kper,sn)
c     write(*,*)'g_n_emish ,sn',sn
c         do i=1,3
c           do  j=1,3
c            write(*,*)'i,j,g_n(i,j)',i,j,sn(i,j)
c           enddo
c         enddo
c-----resonance condition uses the delta function with argument
c     g_delta=gamma-nll*p_par-n*y
c     the derivative from this argument d(g_delta)/dp_par=denom
      
c      denom=dabs(gamma*kpar-vpardvt/cdvt2*omega)
      

      call  unde_intgrh(vperdvt,vpardvt,cdvt,y,kpar,theta,n,i_fkin,
     +   r,z,phi,u_n,u_flcur_n)

c      coeff=u_n/denom
c      coeff_fc=u_flcur_n/denom
       coeff=u_n
       coeff_fc=u_flcur_n

       g_n(1,1)=coeff*sn(1,1)  
       g_n(1,2)=coeff*sn(1,2)
       g_n(1,3)=coeff*sn(1,3)
       g_n(2,2)=coeff*sn(2,2)
       g_n(2,3)=coeff*sn(2,3)
       g_n(3,3)=coeff*sn(3,3)

       g_n(2,1)=-g_n(1,2)
       g_n(3,1)= g_n(1,3)
       g_n(3,2)=-g_n(2,3)      


       g_fluctcur_n(1,1)=coeff_fc*sn(1,1)  
       g_fluctcur_n(1,2)=coeff_fc*sn(1,2)
       g_fluctcur_n(1,3)=coeff_fc*sn(1,3)
       g_fluctcur_n(2,2)=coeff_fc*sn(2,2)
       g_fluctcur_n(2,3)=coeff_fc*sn(2,3)
       g_fluctcur_n(3,3)=coeff_fc*sn(3,3)
      
      g_fluctcur_n(2,1)=-g_fluctcur_n(1,2)
      g_fluctcur_n(3,1)=-g_fluctcur_n(1,3)
      g_fluctcur_n(3,2)=-g_fluctcur_n(2,3)
      
      return
      end


      subroutine s_calcnh(n,vperdvt,vpardvt,kper,sn)
c     calculates tensor sn=p_perp*S(n)
c     input
c      n      - the number of the given EC harmonic
c      vperdvt,vpardvt are the components of momentum divided by m*Vt
c      kper     =k_perpendicular*vt/omegac (algebraic)
c     output
c      sn(3,3)- double complex  sn=p_perp*S(n)

      implicit none
c     input
      integer n
      double precision vperdvt,vpardvt,kper
c     output  
      double complex sn(3,3)
c     external
c     DBESJ
     
c     locals ATTENTION!!!: the local variables were real and complex for besj
      double complex i    
      double precision  bj,bj1,bj_prim,b,d,b_abs
      double precision  bjs
      integer n_abs,ier,k,ir,jr
      double complex snc(3,3)

      integer nz ! should be =0, If nz=1 the bessel function(dBESJ) will be J=0

      
      i=dcmplx(0.0d0,1.0d0)
c-----Calculation of the Bessel function bj= J_n(b) and 
c     its derivative bj_prim=dJ_n(b)/db.
c-----d is the relative error (input),ier -error switch,ier=0 is OK

c      write(*,*)'s_calcnh n,vperdvt,vpardvt,kper',
c     +n,vperdvt,vpardvt,kper
        
      n_abs=iabs(n)            ! the order of the Bessel function

      k=1                      ! coefficient =1 for J_n with not-negative n      
      if(n.lt.0) k=(-1)**n_abs ! coefficient for Bessel function J_n with n<0

      b=kper*vperdvt           ! Bessel function argument
      b_abs=dabs(b)
      
c      write(*,*)'s_calcnh b_abs',b_abs
       goto10

c      d=1.e-3                  ! relative accuracy for Bessel function besj
c      d=1.e-5
      
c       call besj(b_abs,n_abs,bj,d,ier)

       call DBESJ(b_abs,dble(n_abs),1,bj,nz)

c      write(*,*)'s_calcnh b_abs,n_anbs,bj,nz',b_abs,n_abs,bj,nz
      if(b.lt.0.d0) bj=bj*(-1)**n_abs   !negative argument b for J_(n_abs)(b)
 
c      if(ier.ne.0) then
c       write(*,*)'s_calcnh in besj ier should  =0 but ier=',ier
c      write(*,*)'b,n,n_abs,bj',b,n,n_abs,bj 
c      endif

c      call  besj(b_abs,n_abs+1,bj1,d,ier) 

      call DBESJ(b_abs,n_abs+1.d0,1,bj1,nz)

      if(b.lt.0.d0) bj1=bj1*(-1)**(n_abs+1) !negative argument b for J_(n_abs+1)(b)
 
c      if(ier.ne.0) then 
c      write(*,*)'s_calcnh in besj ier should =0 but ier=',ier
c      write(*,*)'b,n,n_abs+1,bj1',b,n,n_abs+1,bj1 
c      endif
     
      if(n.eq.0) then 
        bj_prim=-bj1
      else
c-------n_abs.ne.0
        if(b.eq.0.d0)then 
          if(n_abs.eq.1)then
            bj_prim=0.5d0
          else
c-----------n_abs.ge.2
            bj_prim=0.0d0     
          endif
        else 
c---------argument b.ne.0
          bj_prim=-bj1+n_abs*bj/b
        endif
      endif
      bj=k*bj               ! for negative n
      bj_prim=k*bj_prim     ! for negative n
      write(*,*)'in s_calcnh bj,bj_prim',bj,bj_prim
 10   continue
      call bes_calc(b,n,bj,bj_prim)
c      write(*,*)'relat_tens new n_abs,bj,bj_prim',n_abs,bj,bj_prim
      
      bjs=bj*bj

c-----complex 
c      write(*,*)'s_calch kper',kper
       snc(1,1)=n*n*bj*bj/(kper*kper)
c      write(*,*)'snc(1,1)',snc(1,1)
      snc(1,2)=-i*n*bj*bj_prim*vperdvt/kper
c      write(*,*)'snc(1,2)',snc(1,2)
cS     snc(1,2)=i*n*bj*bj_prime*vperdvt/kper

      snc(2,1)=-snc(1,2)
      snc(2,2)=vperdvt*vperdvt*bj_prim*bj_prim
      snc(1,3)=n*vpardvt*bj*bj/kper
      snc(3,1)=snc(1,3)
      snc(2,3)=i*vpardvt*vperdvt*bj*bj_prim
cS     snc(2,3)=-i*vpardvt*vperdvt*bj*bj_prim

      snc(3,2)=-snc(2,3)
      snc(3,3)=vpardvt*vpardvt*bj*bj
     
      do ir=1,3
        do jr=1,3
          sn(ir,jr)=snc(ir,jr)
        enddo
      enddo
      
            
      return
      end

      subroutine emis_tensh(X,Y,T_kev,nll_in,np_in,cdvt,vmaxdt,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +i_fkin,r,z,phi,
     +aK,fluctcur)

c-----calulates the anti-hermitian relativistic dielectric tensor aK
c     for electron plasma
c     and the correlation tensor fluctcur=G  for the fluctuating current density
c
c     see the article
c     R.W.Harvey at all Phys.Fluids B 5(2),Febrary 1993, page 446
c     G is the correlation tensor    formula(7)
c     aK  formula (6)
c
c     INPUTS:
c
c      X = (fpe/f)**2
c      Y = fce/f
c      T_kev  - electron temperature
c      nll_in - parallel index of refraction N.
c      np_in  - perpendicular refractive index N
c      cdvt   =clight/Vt
c      vmaxdt =Vmax/Vt
c      n_relt_harm1  - min number of used cyclotron jn harmonics,
c      n_relt_harm2  - max number of used cyclotron jn harmonics,
c                      n_relt_harm1 <= jn <= n_relt_harm2 

c      n_relt_intgr - the number of points for the integration over p_perp
c      i_fkin =0 the usage of the analytical relativistic distributin
c             =1 the usage of the numerical 3D distribution from diskf or netcdfnm.nc 
c                written be CQL3D code or created analytically at mesh points
c      r      - the major radius (it is used for i_fkin=1)
c      z      - the vertical coordinate  (it is used for i_fkin=1)
c      phi    - toroidal angle
c     OUTPUT:
c      aK(3,3):  the nine components of the anti-hermition
c                part dielectric relativistic tensor
c                evaluated at (X,Y,Te,nll,np,n)
c      
c      fluctcur(3,3) the correlation tensor fluctcur=G (egr)  
c                    for the fluctuating current density

      implicit none
c     input 
      double precision X,Y,T_kev,nll_in,np_in,cdvt,vmaxdt
      integer n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin
      double precision r,z,phi

c     output
      double complex aK(3,3),fluctcur(3,3)
    
c     local
      double precision c,mass_e,k_1_kev,nll,np,nlls,theta,pi,kper,
     .kpar,dens,xpidens,xpidensmcs,p
      double complex integral(3,3)
      double complex fluctcur_n(3,3)
      integer jn
  
      integer ires
      double precision v0dc,vmax1dc,vpar0dc,vpar1dc,vpar2dc
c      ires  =0 no resonaces,=1 ellipse,=2 parabola,=3 hyperbole 
c      v0dc,vmax1dc,vpar0dc,vpar1dc,vpar2dc the parameters of the resonace courve

c     external zeroK, npnllmin, ec_cond, intgr_rl,fdens_fdist
      double precision fdens_fdist
 
      c =2.99792458d10          !speed of light     (cm/sec)
      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)      

      theta=mass_e*c**2/(k_1_kev*T_kev)
      pi=4.d0*datan(1.d0)

      call npnllmin(nll_in,np_in,nll,np)
      nlls = nll**2
c-----normalized k
      kper=np/(Y*cdvt)
      kpar=nll/(Y*cdvt)
c-----initialization aK=0
      call zeroK(aK)
      call zeroK(fluctcur)      
c-----the loop over the cyclotron harmonics
      kper=np/(Y*cdvt)
      kpar=nll/(Y*cdvt)


      do jn=n_relt_harm1,n_relt_harm2
       
       
c        write(*,*)'emis_tensh jn,np,Y,cdvt,kper,nll',
c     +  jn,np,Y,cdvt,kper,nll

cSm060315 -jn
        call intgr_emish(-jn,nll,kper,Y,theta,cdvt,vmaxdt,
     +  n_relt_intgr,
     +  i_fkin,r,z,phi,
     +  integral,fluctcur_n,ires,v0dc,vmax1dc,vpar0dc,vpar1dc,vpar2dc)

c        write(*,*)'emis_tensh after end intgrl_emish'
c        write(*,*)' integral(1,1),integral(1,2)',
c     +  integral(1,1),integral(1,2)

        aK(1,1)=aK(1,1)+integral(1,1)
        aK(1,2)=aK(1,2)+integral(1,2)
        aK(1,3)=aK(1,3)+integral(1,3)
        aK(2,2)=aK(2,2)+integral(2,2)
        aK(2,3)=aK(2,3)+integral(2,3)
        aK(3,3)=aK(3,3)+integral(3,3)

c         write(*,*)'emis_tensh aK(1,1),aK(1,2)',aK(1,1),aK(1,2)

        fluctcur(1,1)=fluctcur(1,1)+fluctcur_n(1,1)
        fluctcur(1,2)=fluctcur(1,2)+fluctcur_n(1,2)
        fluctcur(1,3)=fluctcur(1,3)+fluctcur_n(1,3)
        fluctcur(2,2)=fluctcur(2,2)+fluctcur_n(2,2)
        fluctcur(2,3)=fluctcur(2,3)+fluctcur_n(2,3)
        fluctcur(3,3)=fluctcur(3,3)+fluctcur_n(3,3)


10      continue      
      enddo   !jn

      if (i_fkin.eq.0) then       
        dens=1.d0
      else      
        dens=fdens_fdist(r,z,phi)
      endif
      
c-----normalization for the unit electron density

c      xpidens=X*pi/dens 
      xpidens=X*2.d0*pi**2*dabs(kpar)/dens

      aK(1,1)=-xpidens*aK(1,1)
      aK(1,2)=-xpidens*aK(1,2)
      aK(1,3)=-xpidens*aK(1,3)
      aK(2,2)=-xpidens*aK(2,2)
      aK(2,3)=-xpidens*aK(2,3)
      aK(3,3)=-xpidens*aK(3,3)
c      write(*,*)'emis_tensh xpi aK(1,1),aK(1,2)',aK(1,1),aK(1,2)
      aK(2,1)=-aK(1,2)
      aK(3,1)= aK(1,3)
      aK(3,2)=-aK(2,3)
       
c      xpidensmcs=X*pi/(2.d0*pi)**5*(mass_e*c**2)

      xpidensmcs=1.d0/(2.d0*pi)**3*0.5d0*mass_e*(c/cdvt)**2*
     +X*Y/dens

      fluctcur(1,1)=fluctcur(1,1)*xpidensmcs
      fluctcur(1,2)=fluctcur(1,2)*xpidensmcs
      fluctcur(1,3)=fluctcur(1,3)*xpidensmcs
      fluctcur(2,2)=fluctcur(2,2)*xpidensmcs
      fluctcur(2,3)=fluctcur(2,3)*xpidensmcs
      fluctcur(3,3)=fluctcur(3,3)*xpidensmcs
      
      fluctcur(2,1)=-fluctcur(1,2)
      fluctcur(3,1)= fluctcur(1,3)
      fluctcur(3,2)=-fluctcur(2,3)
      

      return
      end


      subroutine calc_emis_coefh(xe,ye,T_kev,cnpar,cnper,cnray,
     +cdvt,vmaxdt,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,phi,
     +cex,cey,cez,cflown,vgroup,frqncy,
     +iabsorp_collisional,coll_mult,tempe,dense,zeff,
     +al_emis,j_emis,temprad)
c---- calcultes the emission absorptivity al_emis (1/cm) 
c     emissivity j_emis (erg*sec/cm**3)  and
c     radiation  temperature temprad

      implicit none

c-----input
      double precision xe,ye,T_kev,cnpar,cnper,cnray,z,r,phi
c      cdvt   =clight/Vt
c      vmaxdt =Vmax/Vt
      double precision cdvt,vmaxdt
c     cnray  is the ray refractive index N_ray
      integer n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin
      double complex cex,cey,cez
c     wave energy flux for the normalized field |E|=1
c     cflown=(vectorB*vector_B^*+vector_E*d(omega*tensor_eps_herm)/dw*vector_E)
      double complex cflown  
      double precision vgroup  !in clight
      double precision frqncy   !GhZ
c-----data for absorption calculations
      integer iabsorp_collisional ! =0 do not calculate collisional absorption
                                  ! =1 to calculte collisonal absorption
      double precision coll_mult  ! multiplier for coll absorption
      double precision tempe,     !electron temperature  (keV)  
     & dense,                     !electron density   (10*13/cm**3	 
     & zeff                       !plasma charge
c-----external
      double precision rrind     
c-----output  
      double precision al_emis,j_emis,temprad
c-----local        
      double complex aK(3,3),fluctcur(3,3),ce(3)
      double precision v_gr_cm_sec !cm/sec
      double precision omega       !2*pi*f, here f in [HZ]
      double precision omegpedce,omegdce!omega_pe/omega_ce and omega/omega_ce
      double precision cn          !refructive index
      double precision clight,pi      
   
      call emis_tensh(xe,ye,T_kev,cnpar,cnper,cdvt,vmaxdt,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,phi,
     +aK,fluctcur)

      ce(1)=cex
      ce(2)=cey
      ce(3)=cez
      clight=2.99792458d10          !cm/sec
      v_gr_cm_sec=clight*vgroup     !cm/sec
      pi=4.d0*datan(1.d0)
      omega=2.d0*pi*frqncy*1.d9     !1/esc
           
c      write(*,*)'calc_coef omega,frqncy,cn,cnray,cflown',
c     +omega,frqncy,cn,cnray,cflown
     
c      write(*,*)'in calc_emis_coef aK',aK
c      write(*,*)'in calc_emis_coef fluctcur',fluctcur
c      write(*,*)'in calc_emis_coef cflown',cflown
c      write(*,*)'in calc_emis_coef ce',ce

      call emiss_coef(aK,fluctcur,ce,cflown,v_gr_cm_sec,omega,cnray,
     +iabsorp_collisional,coll_mult,tempe,dense,zeff,
     +al_emis,j_emis)  

c     write(*,*)'in calc_emis_coef al_emis,j_emis',al_emis,j_emis

c-----radiation temperature      
      if(al_emis.ne.0.d0)  then
         temprad=(2.d0*pi)**3*(clight)**2/(omega*cnray)**2
     +             *j_emis/al_emis/1.6022d-9
c         write(*,*)'in calc_emis_coef T_kev,temprad=',T_kev,temprad
      else
         temprad=0.d0
      endif

      return
      end 

      subroutine maxw_test
c-----check the density from maxwellian distribution
      implicit none
      integer i,imax
      double precision theta,bk2,pi,gamma,h,v,den,f_maxw,clight
      theta=120.d0

      call besk2as(theta,bk2)
      pi=4.d0*datan(1.d0) 
      den=0.d0 
      imax=100
      clight=3.0d10
cSm0305125
c      h=0.2d0/dreal(imax)
      h=0.2d0/dfloat(imax)

      do i=1,imax
         v=i*h !v/c
         gamma=dsqrt(1+v**2)
         f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)
         den=den+f_maxw*v**2*h
      enddo
      den=den*4.d0*pi
      write(*,*)' test den',den
      return
      end
         
      subroutine rayrind(i_rrind,nbulk,mass_ar,x_ar,y_ar,
     .t_av_ar,tpop_ar,vflow_ar,nll_in,np_in,n_ray)
c-----calculte ray refructive inedx for emission
c     G.Bekefi, Radiation Processes in Plasmas,(John Wiley and Sonc Inc.)
c     1966 page 32, formula (1.121)
c-----input
c     i_rrind chooses the subroutine to calculate N_ray 
c     i_rrind=0  N_ray=N
c     i_rrind =1 from cold electron plasma using rrind   
c     i-rrind =2 from hot non-relativistic dispersion relation  
c     
c     the parameters for rrind_hotpl
c      nbulk the total number of plasma species
c      mass_ar(*) - the masses  of the plasma species (in electron mass) 
c      x_ar(*) = (fpe/f)**2
c      y_ar(*) = fce/f   It is the algebraic number
c      t_av_ar(*)=average temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar(*) = T_perp/T_parallel  
c      vflow_ar(*)    cm/sec
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c      k_sum(3,3) -  the nine components of the dielectric tensor
c                    evaluated at (mass,X,Y,Te,tpop,vflow,nll,np)           
c-----output 
c     n_ray      is the ray refractive index
      implicit none
c-----input

      integer i_rrind,nbulk
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double precision np_in,nll_in

c-----externals
      double precision rrind,rrind_hotpl

c-----output
      double precision n_ray
c-----locals
      double precision omegpedce,omegdce          

c-----the ray refructive index equals to the wave refractive index cnray=N
      if (i_rrind.eq.0) then
          n_ray=dsqrt(nll_in**2+np_in**2)
          goto10
      endif

c-----ray_refractive index for the cold electron plasma using
c     double precision function rrind
      if (i_rrind.eq.1) then
         omegpedce=dsqrt(x_ar(1))/dabs(y_ar(1)) !omega_pe/omega_ce
         omegdce=1.d0/dabs(y_ar(1))             !omega/omega_ce
         n_ray=rrind(nll_in,np_in,omegdce,omegpedce)
         goto10
      endif

c-----ray_refractive index for the hot non-relativistic plasma using
c     double precision function rrind_hotpl
      if (i_rrind.eq.2) then
         n_ray=rrind_hotpl(nbulk,mass_ar,x_ar,y_ar,
     .          t_av_ar,tpop_ar,vflow_ar,nll_in,np_in)
         goto10
      endif

 10   continue
      return
      end

      
      double precision function f_y(r)
c-----calculate Ye(r,z=0)-0.5.
c     This function will be used to find the major radius of
c     the second EC harmonic reconance point.
c     f_y(r)=omega_ce(r,z=0)/omega-0.5d0=y(r,z=0,phi=0,1)-0.5
      !implicit none
      implicit double precision (a-h,o-z) 
      include 'param.i'
      include 'one.i' !
c-----input
      double precision r !major radius
c-----external
      double precision b,y
c-----local
      double precision yloc
c      write(*,*)'emission.f in f_y r',r
      bmod=b(0.d0,r,0.d0)
c      write(*,*)'emission.f in f_y bmod',bmod
c      yloc=y(0.d0,r,0.d0,1)
c      write(*,*)'f_y yloc',yloc
      f_y=y(0.d0,r,0.d0,1)-0.5d0
c      write(*,*)'emission.f in f_y',f_y
      return
      end

      double precision function r_2nd_harm(t)
c     input parameter can be arbitray. It is not used inside the function. 
c-----calculate the major radius of the EC second harmonic
c     resonance point omega_ce(r,z=0)/omega
c
      implicit none !double precision (a-h,o-z) 
      include 'param.i'
      include 'five.i' ! gives rmax and rmin of the tokamak
cSm050929
      include 'three.i' ! gives xma
      include 'one.i' 
      
      real*8 t !INPUT
c-----output
c     r_2nd_harm ! the major radius where w_ce(r,z=0)/w=0.5
c-----external
      double precision f_y     !,rtbis [not used]

      integer JMAX,J
c     PARAMETER (JMAX=40)
      parameter(JMAX=100)      
      double precision FMID,F,X1,X2,XACC,XMID,DX
      real*8 b,y ! external
      real*8 phi0,y_e !local

      xacc=1.d-5 
      !write(*,*)'in r_2nd_harm bef rtbis rmin,rmax',rmin,rmax 
      !write(*,*)'in r_2nd_harm bef rtbis rmin,rmax',rmin,rmax
c      r_2nd_harm=rtbis(f_y,rmin+0.1d0,rmax-0.05d0,xacc)
     
      
c      X1=  rmin+0.1d0
c      X2=  rmax-0.05d0
      X1=  rmin
      X2=  rmax
      write(*,*)'emission r_2nd_harm X1,X2,XACC',X1,X2,XACC

      FMID=f_y(X2)
      write(*,*)'emission r_2nd_harm  X2,FMID',X2,FMID  
      F=f_y(X1)
      write(*,*)'emission r_2nd_harmX1, F',X1,F  

c      IF(F*FMID.GE.0.d0) PAUSE 
c     +'(emission.f r_2nd_harm: Root must be bracketed for bisection.'
      IF(F*FMID.GE.0.d0) then
         write(*,*)'WARNING: in emission.f r_2nd_harm '//
     &'Root must be bracketed for bisection '//
     &'the code set r_2d_harm =rmax*(1.'
         write(*,*)'f_y(rmin)',f_y(rmin)
         write(*,*)'f_y(rmax)',f_y(rmax)
         phi0=0.d0
         bmod=b(yma,xma,phi0)
         y_e=y(yma,xma,phi0,1)
         r_2nd_harm=2.d0*y_e*xma
         write(*,*)'xma,r_2nd_harm',xma,r_2nd_harm
         return
      ENDIF

      IF(F.LT.0.d0)THEN
        r_2nd_harm=X1
        DX=X2-X1
      ELSE
        r_2nd_harm=X2
        DX=X1-X2
      ENDIF
       DO 11 J=1,JMAX
         DX=DX*.5d0
         XMID= r_2nd_harm+DX
         FMID=f_y(XMID)
         write(*,*)'emission r_2nd_harm J,DX,XMID,FMID',J,DX,XMID,FMID  
         IF(FMID.LE.0.d0) r_2nd_harm=XMID
         IF(dabs(DX).LT.XACC .OR. FMID.EQ.0.d0) RETURN
          IF(dabs(DX).LT.XACC .OR. dabs(FMID).LT.XACC) RETURN
 11    CONTINUE 
       PAUSE 'r_2nd_harm:too many bisections,increase parameter JMAX i'

      return
      end function r_2nd_harm


      subroutine g_n_emis_theta(p_perp,p_par,y,nll,np,theta,n,i_fkin,
     &r,z,phi,g_n,g_fluctcur_n)
c     calculates under integral complex marix functions
c     g_n(3,3)=sum{k}g_nk and g_fluctcur_n(3,3)=sum{k}g_fluctcur_nk

c     input
c       p_perp  - perpendicular momentum divided by (mc)
c       p_par     paralell momentum along resonance curve
c       y       = omega_ce/omega for the electron rest mass    
c       nll     - N_parallel
c       np      - N_perpendicular 
c       theta   - mc**2/T
c       n       - n number of the given resonance harmonic
c       i_fkin =0 the usage of the analytical relativistic distributin
c              =1 the usage of the numerical 3D distribution from diskf 
c                 written be CQL3D code
c       r       - the major radius (it is used for i_fkin=1)
c       z       - the vertical coordinate  (it is used for i_fkin=1)
c       phi     - toroidal angle(it is used for i_fkin=1)
c     output
c       g_n(3,3) under integral matrix double complex function
c       g_fluctcur_n(3,3)
c-------------------------------------------------------
      implicit none
c-----input        
      double precision p_perp,p_par,y,nll,np,theta,r,z,phi
      integer n,i_fkin
c-----external root_res, s_calcn, zeroK, unde_intgr
      double precision u_n,u_flcur_n    
c-----output
      double complex g_n(3,3), g_fluctcur_n(3,3)
      
c-----local
      double complex sn_k(3,3),g_nk(3,3),g_fluctcur_nk(3,3)
      double precision gamma, p_par_rl(2),coeff,eps,ps,ps_max,pi,
     +dgdp_par,u_nk,u_flcur_nk,coeff_fc
      integer k,kmax,i,j
c     kmax - total number of the resonace condition root p_perp_k(p_perp)
      
c for test
      double precision resonan

      pi=4.d0*datan(1.d0)
c-----calculations ot the roots p_par_rl of the resonace condition
c-----gamma=N_par*p_par_rl+nY

      call root_res(p_perp,nll,n,Y,kmax,p_par_rl)

c      write(*,*)'g_n_emis number of roots p_par_rl(p_perp) kmax',kmax

c-----test resonance condition
c      resonan=dsqrt(1.d0+p_perp**2+p_par_rl(2)**2)-nll*p_par_rl(2)-n*Y
c      write(*,*)'1 resonan', resonan
c      resonan=dsqrt(1.d0+p_perp**2+p_par_rl(1)**2)-nll*p_par_rl(1)-n*Y
c      write(*,*)'2 resonan', resonan 
c      resonan=dsqrt(1.d0+p_perp**2+p_par**2)-nll*p_par-n*Y
c      write(*,*)'3 resonan', resonan
c      write(*,*)'g_n_emis_theta number of roots p_par_rl(p_perp) kmax',
c     &kmax,' p_par_rl ',p_par_rl,'perp,p_par',p_perp,p_par
    
c-----initialize g_n
      call zeroK(g_n)
      call zeroK(g_fluctcur_n)

      if (kmax.eq.0) goto 20
  
c      eps=1.d-9 ! the min value of Maxwell exponent
c      ps_max=(1.d0-dlog(eps)/theta)**2-1.d0
      
      ps=p_perp*p_perp+p_par*p_par
      gamma=dsqrt(1.d0+ps)
   
c     if(ps.gt.ps_max) then
c        write(*,*)'emission.f g_n_emis ps,ps_max',ps,ps_max
c        goto 10
c      endif

c         write(*,*)'g_n_emis k,p_perp,p_par_rl(k),ps',
c     .   k,p_perp,p_par_rl(k),ps

c      write(*,*)'before s_calcn n,p_perp,p_par,np,y',
c     & n,p_perp,p_par,np,y
     
      call s_calcn(n,p_perp,p_par,np,y,sn_k)
c      write(*,*)'g_n_emis sn_k',sn_k

c-----resonance condition uses the delta function with argument
c     g_delta=gamma-nll*p_par-n*y
c     the derivative from this argument d(g_delta)/dp_par=dgdp_par

      dgdp_par=dabs((p_par-nll*gamma)/gamma)

c      write(*,*)'p_perp,p_par,nll,gamma',p_perp,p_par,nll,gamma
c      write(*,*)'dgdp_par',dgdp_par

cSm030415
c      if(dgdp_par.lt.1.d-3)then
c       write(*,*)'dgdp_par p_par,nll,p_perp,gamma',
c     +            p_par,nll,p_perp,gamma
c      endif
    

      call  unde_intgr(p_perp,p_par,y,nll,theta,n,i_fkin,
     +r,z,phi,u_nk,u_flcur_nk)

c      write(*,*)'in g_n_emis after unde_integr u_nk,u_flcur_nk',
c     +u_nk,u_flcur_nk

      coeff=2.d0*pi*u_nk/dgdp_par
      coeff_fc=2.d0*pi*u_flcur_nk/dgdp_par
c      write(*,*)'in g_n_emis u_nk,u_flcur_nk',u_nk,u_flcur_nk
c      write(*,*)'in g_n_emis coeff,coeff_fc',coeff,coeff_fc

      g_nk(1,1)=coeff*sn_k(1,1)  
      g_nk(1,2)=coeff*sn_k(1,2)
      g_nk(1,3)=coeff*sn_k(1,3)
      g_nk(2,2)=coeff*sn_k(2,2)
      g_nk(2,3)=coeff*sn_k(2,3)
      g_nk(3,3)=coeff*sn_k(3,3)

      g_nk(2,1)=-g_nk(1,2)
      g_nk(3,1)= g_nk(1,3)
      g_nk(3,2)=-g_nk(2,3)      


      g_fluctcur_nk(1,1)=coeff_fc*sn_k(1,1)  
      g_fluctcur_nk(1,2)=coeff_fc*sn_k(1,2)
      g_fluctcur_nk(1,3)=coeff_fc*sn_k(1,3)
      g_fluctcur_nk(2,2)=coeff_fc*sn_k(2,2)
      g_fluctcur_nk(2,3)=coeff_fc*sn_k(2,3)
      g_fluctcur_nk(3,3)=coeff_fc*sn_k(3,3)


      do i=1,3
         do j=1,3
            g_n(i,j)=g_n(i,j)+g_nk(i,j)
            g_fluctcur_n(i,j)=g_fluctcur_n(i,j)+g_fluctcur_nk(i,j)
         enddo
      enddo
  
 10   continue
      
 20   continue

      return
      end

      subroutine g_n_calc_theta(p_perp,p_par,y,nll,np,theta,n,i_fkin,
     &r,z,phi,g_n)
c     under integral complex marix function G_n(3,3)
c     input
c       p_perp  - perpendicular momentum divided by (mc)
c       p_par     parallel momentum divided by (mc) along resonance curve
c       y       = omega_ce/omega for the electron rest mass    
c       nll     - N_parallel
c       np      - N_perpendicular 
c       theta   - mc**2/T
c       n       - n number of the given resonance harmonic
c       i_fkin =0 the usage of the analytical relativistic distributin
c              =1 the usage of the numerical 3D distribution from diskf 
c                 written be CQL3D code
c       r       - the major radius (it is used for i_fkin=1)
c       z       - the vertical coordinate  (it is used for i_fkin=1)
c       phi     - toroidal angle(it is used for i_fkin=1)
c     output
c       g_n(3,3) under integral matrix double complex function
c-------------------------------------------------------
      implicit none
c-----input        
      double precision p_perp,y,nll,np,theta,r,z,phi,p_par
      integer n,i_fkin
c-----external root_res, s_calcn, zeroK, u_n
      double precision u_n    

c-----output
      double complex g_n(3,3) 

c-----local
      double complex sn_k(3,3),g_nk(3,3)
      double precision gamma, p_par_rl(2),coeff,eps,ps,ps_max,pi,
     +dgdp_par,u_nk
      integer k,kmax,i,j
c     kmax - total number of the resonace condition root p_perp_k(p_perp)
      
    

      pi=4.d0*datan(1.d0)
c-----calculations ot the roots p_par_rl of the resonance condition
c-----gamma=N_par*p_par_rl+nY

      call root_res(p_perp,nll,n,Y,kmax,p_par_rl)
c      write(*,*)'g_n_calc_theta number of roots p_par_rl(p_perp) kmax',
c     &kmax,' p_par_rl ',p_par_rl,'perp,p_par',p_perp,p_par
c      write(*,*)'in g_n_calc: p_perp, nll,n,Y,kmax,p_par_rl'
c     &,p_perp, nll,n,Y,kmax,p_par_rl
c-----initialize g_n
      call zeroK(g_n)

      if (kmax.eq.0) goto 20
  
c     eps=1.d-9 ! the min value of Maxwell exponent
c     p s_max=(1.d0-dlog(eps)/theta)**2-1.d0
c     write(*,*)'g_n_calc kmax',kmax      
        
      ps=p_perp*p_perp+p_par*p_par
      gamma=dsqrt(1.d0+ps)

c     if(ps.gt.ps_max) then
c        write(*,*)'relat_tens g_n k,ps,ps_max',k,ps,ps_max
c        goto 10
c     endif
   
      call s_calcn(n,p_perp,p_par,np,y,sn_k)
c     write(*,*)'g_n_calc after s_calcn: sn_k',sn_k  

c-----resonance condition uses the delta function with argument
c     g_delta=gamma-nll*p_par-n*y
c     the derivative from this argument d(g_delta)/dp_par=dgdp_par

c     dgdp_par=(p_par_rl-nll*y)/gamma
      dgdp_par=dabs((p_par-nll*gamma)/gamma)
      u_nk=u_n(p_perp,p_par,y,nll,theta,n,i_fkin,r,z,phi)

cSm0960306
c     write(*,*)'g_n   dgdp_par',  dgdp_par

      coeff=2.d0*pi*u_nk/dgdp_par

c     write(*,*)'g_n k,p_perp,p_par_rl(k)',k,p_perp,p_par_rl(k)
c     write(*,*)'g_n dsqrt(ps),u_nk',dsqrt(ps),u_nk
c     write(*,*)'g_n u_nk,dgdp_par,coeff',u_nk,dgdp_par,coeff 
c     write(*,*)'coeff',coeff           

      g_nk(1,1)=coeff*sn_k(1,1)  
      g_nk(1,2)=coeff*sn_k(1,2)
      g_nk(1,3)=coeff*sn_k(1,3)
      g_nk(2,2)=coeff*sn_k(2,2)
      g_nk(2,3)=coeff*sn_k(2,3)
      g_nk(3,3)=coeff*sn_k(3,3)

      g_nk(2,1)=-g_nk(1,2)
      g_nk(3,1)= g_nk(1,3)
      g_nk(3,2)=-g_nk(2,3)        

      do i=1,3
         do j=1,3
            g_n(i,j)=g_n(i,j)+g_nk(i,j)
         enddo
      enddo

 10   continue
       
 20   continue

      return
      end

      subroutine get_rtem0_from_one(rtem0_loc)
      implicit none       
      include 'param.i'
      include 'one.i'
      double precision  rtem0_loc
      rtem0_loc=rtem0
      return
      end
      
      subroutine calc_emission_integral_array_by_adaptive_simpson
     &(y,nll,np,theta,n,i_fkin,r,z,phi,p_permax,epsi,
     &integral,fluctcur_n)
c------------------------------------------------------------------------
c     alculates the matrix: double complex integrals
c     for the relativistic electron plasma 
c     I(n)^=integral(0<=p_perp<=p_perp0)g_n, g_n={sum_k=1,2,G_nk(p_perp)
c     fluctcar(n)^=integral(0<=p_perp<=p_perp0)g_fluctcur_n.
c               g_fluctcur_n={sum_k=1,2,g_fluctcur_nk(p_perp)
c     for the EC harmonic with number 'n'
c-------------------------------------------------------------------------
       implicit none
c------------------------------------------------------------------------
c      input
c------------------------------------------------------------------------
      integer n,i_fkin
      real*8 y,nll,np,theta,z,r,phi,p_permax ,epsi
c     n        is EC harmonic number
c     i_fkin   =0 the usage of the analytical relativistic distributin
c              =1 the usage of the numerical 3D distribution from diskf 
c                   written be CQL3D code
c     nll      is N_parallel
c     theta    = mc**2/T
c     r        is the major radius (it is used for i_fkin=1)
c     z        is  the vertical coordinate  (it is used for i_fkin=1)
c     phi      is  the toroidal angle (it is used for i_fkin=1)
c     p_permax is  the maximal boundary of integration
c     epsi is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------------------
c     output
c       integral(3,3) double complex integral from G_nk over 0<=p_perp=<p_perpmax
c       fluctcur(3,3) double complec integral from fluct_cur_nk over 0<=p_perp=<p_perpmax
c------------------------------------------------------------------------
      double complex integral(3,3),fluctcur_n(3,3) 

c-------------------------------------------------------------------------
c     locals
c--------------------------------------------------------------------------
      integer k_max !max number of elements
      parameter (k_max=12)
      integer IST,IFS,k   
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
    
      integer n_l,i_fkin_l
      real*8 y_l,nll_l,np_l,theta_l,z_l,r_l,phi_l 
      common /fcn_input/y_l,nll_l,np_l,theta_l,n_l,i_fkin_l,
     &z_l,r_l,phi_l
 
      EXTERNAL FCN_emission_vector

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
c     set common /fcn_input/
c-----------------------------------------------------------------
      y_l=y
      nll_l=nll
      np_l=np
      theta_l=theta
      n_l=n
      i_fkin_l=i_fkin
      z_l=z
      r_l=r
      phi_l=phi
c-----------------------------------------------------------------
c     integration boundaries
c------------------------------------------------------------------
c      A = 0.0d0
      A = 1.d-5
c      A = 1.d-3

      B = p_permax-1.d-5
     
      CALL ASMP(FCN_emission_vector,k_max,fcn,
     &A,B,EPSI,LVMAX,STACK,IST,
     &FSTACK,IFS,SUM,IFLAG) 

c      do k=1,k_max 
c        write(*,*)'k=',k 
c        PRINT 3,SUM(k)
c      enddo  
     

      IF(IFLAG .EQ. 0) THEN 
c        PRINT 4   
        write(*,*) '          WITH NO BAD SUBINTERVALS:'
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
      integral(1,1)=dcmplx(Sum(1),0.d0) 
      integral(1,2)=dcmplx(0.d0,Sum(2))
      integral(1,3)=dcmplx(Sum(3),0.d0) 
      integral(2,2)=dcmplx(Sum(4),0.d0)
      integral(2,3)=dcmplx(0.d0,Sum(5))
      integral(3,3)=dcmplx(Sum(6),0.d0) 

      fluctcur_n(1,1)=dcmplx(Sum(7),0.d0) 
      fluctcur_n(1,2)=dcmplx(0.d0,Sum(8))
      fluctcur_n(1,3)=dcmplx(Sum(9),0.d0) 
      fluctcur_n(2,2)=dcmplx(Sum(10),0.d0)
      fluctcur_n(2,3)=dcmplx(0.d0,Sum(11))
      fluctcur_n(3,3)=dcmplx(Sum(12),0.d0)
 
      return
      END
 
      subroutine fcn_emission_vector(p_perp,fcn)
c---------------------------------------------------------------------------
c     culculates under integral function for
c     absorption  anti-hermitian dielectric tesor
c---------------------------------------------------------------------------

      implicit none

      integer k_max !max number of elements
      parameter (k_max=12) 
c----------------------------------------------------------------
c     input
c----------------------------------------------------------------- 
      real*8 p_perp !under integral function argument
c-----------------------------------------------------------------
c     output
c-----------------------------------------------------------------
      real*8 fcn(k_max) !array of under integral fuctions values
c-----------------------------------------------------------------
      integer n_l,i_fkin_l
      real*8 y_l,nll_l,np_l,theta_l,z_l,r_l,phi_l 
      
      common /fcn_input/y_l,nll_l,np_l,theta_l,n_l,i_fkin_l,
     &z_l,r_l,phi_l 
c-----------------------------------------------------------------
c     locals
c----------------------------------------------------------------
      double complex g_n(3,3),fluctcur_n(3,3) !complex under integral functions

      call g_n_emis(p_perp,y_l,nll_l,np_l,theta_l,n_l,i_fkin_l,
     &r_l,z_l,phi_l,g_n,fluctcur_n)

      FCN(1) =dreal(g_n(1,1)) 
      FCN(2) =dimag(g_n(1,2))
      FCN(3) =dreal(g_n(1,3))
      FCN(4) =dreal(g_n(2,2))
      FCN(5) =dimag(g_n(2,3))
      FCN(6) =dreal(g_n(3,3))

      FCN(7) =dreal(fluctcur_n(1,1)) 
      FCN(8) =dimag(fluctcur_n(1,2))
      FCN(9) =dreal(fluctcur_n(1,3))
      FCN(10) =dreal(fluctcur_n(2,2))
      FCN(11) =dimag(fluctcur_n(2,3))
      FCN(12) =dreal(fluctcur_n(3,3))

      return
      end
