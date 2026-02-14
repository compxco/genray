c       calculate the spectrum of the emission flux
c

      subroutine emission_spectrum(ifreq,iray_start_point,
     &nrayelt_o_cutoff_emis)
c---------------------------------------------------------------------
c     calculate the energy spectrum of the emission 
c---------------------------------------------------------------------
      implicit none
      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'
      include 'emissa.i'
      include 'write.i'
      include 'ions.i'
      include 'oxb.i'  ! gives nrayelt_o_cutoff
c-----input
      integer ifreq,
     & nrayelt_o_cutoff_emis !the OX cutoff point at the ray 
      integer iray_start_point !number of ray launching point              
              !it will put ray trajectories data for central ray only
              !at iray_start_point=1, 
              !for iray_start_point > 1 the rays will not be put in 
              !netcdf file
c-----output


c-----external
      real*8   x,y,tempe,b,gamma1,cn,tpoprho,dense,zeff,vflowrho, zeffi
      external x,y,tempe,b,gamma1,cn,tpoprho,dense,zeff,vflowrho, zeffi
      real*8 tpop_zrp   !external func
      external tpop_zrp !external func

c-----local
      integer i,j,iray_loc,i_geom_optic_loc,n_elt
      integer i_fkin,iabsorp_collisional_l

      real*8 z,r,phi,cnz,cnr,cm,xe,ye,T_kev,
     &ds,dc,cnt,cnpar,cnper,
     &vgrmods,vgroup,
     &te,u(6),deru(6),bf(3),vgr(3),
     &cnray,
     &coll_mult_l,tempe_l,dense_l,zeff_l
cSAP090808 nbulk -> nbulk
c      reaL*8 x_ar(nbulka),y_ar(nbulka),
c     &t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka)
      real*8 x_ar(nbulk), y_ar(nbulk), 
     +      t_av_ar(nbulk), tpop_ar(nbulk), vflow_ar(nbulk)

      complex*16 cex,cey,cez,cflown
cSAP090808  jx_kin_a-> jx_kin
c      real*8 j_emis_x(jx_kin_a),
      real*8 j_emis_x(jx_kin)
      real*8 al_emis_loc,j_emis_loc,temp_rad_em_loc

cSAP090808 jx_kin_a -> jx_kin
c      real*8 wj_emis_x(nrelta,jx_kin_a), 
c     &work_j_emis_x_n(jx_kin_a),work_j_emis_x_np1(jx_kin_a),  
c     &work_in_sn_x_n(jx_kin_a),work_in_0_x_n(jx_kin_a), 
c     &in_sn_x(nrelta,jx_kin_a), ! emission I_n at the detector side
c                                ! of nth bin at s=s_n
c     &in_0_x(nrelta,jx_kin_a),  ! emission at the plasma boundary s=s_1 from one
c                                ! nth bin s_n<s<s_n+1
c     &i_0_x(jx_kin_a),          ! emission I_0 at the plasma boundary from the
c                                ! ray at s=s(1)=0
c     &i_0sn_x(nrelta,jx_kin_a)  ! sum{k=1,n}[in_0(k)]! emission at the plasma boundary
     
      real*8 work_j_emis_x_n(jx_kin), work_j_emis_x_np1(jx_kin),
     &work_in_sn_x_n(jx_kin), work_in_0_x_n(jx_kin),
     &i_0_x(jx_kin)    ! emission I_0 at the plasma boundary from the
                       ! ray at s=s(1)=0

      integer jx_kin_a ! YuP: added for 4 arrays below.
      parameter (jx_kin_a=200) ! if jx_kin used instead of jx_kin_a,
                               ! the code overflows (in IVF-10 on PC)
      real*8 wj_emis_x(nrelta,jx_kin_a),
     &in_sn_x(nrelta,jx_kin_a), ! emission I_n at the detector side
               ! of nth bin at s=s_n
     &in_0_x(nrelta,jx_kin_a),  ! emission at the plasma boundary s=s_1 from one
               ! nth bin s_n<s<s_n+1
     &i_0sn_x(nrelta,jx_kin_a)  ! sum{k=1,n}[in_0(k)]! emission at the plasma boundary
     
c---- for testr
      real*8 test_i_0
      integer iray_loc_

c      write(*,*)'emission_spectrum: 1st exec statemnt bombs w viz:ifort'
c      write(*,*)'emission_spectrum:i_emission,i_emission_spectrum=',
c     +     i_emission,i_emission_spectrum,nrayelt_o_cutoff_emis

      if (jx_kin_a .lt. jx_kin)  then  ! Added by YuP
       print*,'Need to increase jx_kin_a parameter in emission_spectrum'
       print*,'set jx_kin_a equal to jx_kin in subr.emission_spectrum()'
       stop
      endif  
      
      if (i_emission.ne.1) go to 10 !nothing to do
      if (i_emission_spectrum.ne.1) go to 10 !nothing to do

c----------------------------------------------------------------------
c     create kinetic energy grid: kinetic_energy_kev(jx_kin)
c----------------------------------------------------------------------
cSAP090808 jx_kin_a -> jx_kin,
c      call kinetic_energy_grid_spectrum(jx_kin_a,jx_kin,    
      call kinetic_energy_grid_spectrum(jx_kin,jx_kin,
     &max_kin_energy_kev,kinetic_energy_kev)
c      write(*,*)'in emission_spectrum.f' 
c      write(*,*)'jx_kin_a,jx_kin,max_kin_energy_kev',
c     &jx_kin,jx_kin,max_kin_energy_kev
c      write(*,*)'kinetic_energy_kev',kinetic_energy_kev

c------------------------------------------------------------------------ 
      if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c--------usage of the analytical relativistic function and its derivatives 
         i_fkin=0
      else 
c--------usage of the mech relativistic function and its derivatives
         i_fkin=1
      endif

      write(*,*)'emission_spectrum i_diskf,i_fkin',i_diskf,i_fkin
      if( iray_start_point.eq.1) then 
c---------------------------------------------------------------
c       central ray data will be saved for writing in netcdf file
c-----------------------------------------------------------------
c        nrayelt_emis=nrayelt_emis_nc(ifreq)

        iray_loc=1
c        write(*,*)'nrayelt_emis',nrayelt_emis
        do n_elt=1,nrayelt_emis

           write(*,*)' n_elt=', n_elt

c          wsn_nc(n_elt,ifreq)=wsn(n_elt)                  !ray distance
c          wr_em_nc(n_elt,ifreq)=wr_em(n_elt)              !r
c          wz_em_nc(n_elt,ifreq)=wz_em(n_elt)              !z
c          wphi_em_nc(n_elt,ifreq)=wphi_em(n_elt)          !phi
c          wal_emis_nc(n_elt,ifreq)=wal_emis(n_elt)        !alpha_emis
c          wj_emis_nc(n_elt,ifreq)=wj_emis(n_elt)          !j_emis
c          wnray_nc(n_elt,ifreq)=wnray(n_elt)              !n_ray 
        
          r=wr_em(n_elt)   
          z=wz_em(n_elt)  
          phi=wphi_em(n_elt)
          cnz=wcnz_em(n_elt)
          cnr=wcnr_em(n_elt)
          cm=wcm_em(n_elt)

c---------calculate the data for emission
          bmod=b(z,r,phi)
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)    !it will be used as negative for electron
          T_kev=tempe(z,r,phi,1)  
          
          u(1)=z
          u(2)=r
          u(3)=phi        
c      nrayelt_emis=nrayelt
c      nrayelt_o_cutoff_emis=nrayelt_o_cutoff

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

c          write(*,*)'before rside1'

          call rside1(0.d0,u,deru)

c          write(*,*)'after rside1'

          i_geom_optic=i_geom_optic_loc

          vgr(1)=deru(1)
          vgr(2)=deru(2)
          vgr(3)=r*deru(3)
          vgrmods=vgr(1)**2+vgr(2)**2+vgr(3)**2
          vgroup=dsqrt(vgrmods) !in clightc

c---------calculate cnray
          do i=1,nbulk
             x_ar(i)=x(z,r,phi,i)
	     y_ar(i)=y(z,r,phi,i)
             if(i.eq.1) y_ar(1)=-y_ar(1)
	     te=tempe(z,r,phi,i) ! kev
	     t_av_ar(i)=te*1000.d0      ! ev 
             tpop_ar(i)=tpop_zrp(z,r,phi,i) !YuP[2024-08-14] was tpoprho(rho,i)
             vflow_ar(i)=vflowrho(rho,i)
          enddo        
 

c          call rayrind(i_rrind,nbulk,dmas,x_ar,y_ar,
c     .             t_av_ar,tpop_ar,vflow_ar,cnpar,cnper,cnray)

c          write(*,*)'cnray from rayrind ',cnray
                  
          cnray=wnray(n_elt)    

          write(*,*)'cnray from wbray(n_elt) ',cnray

c---------calculate electric field (cex,cey,cez) and flux cflown     
c          write(*,*)'before elf_emis' 
          call elf_emis(0.d0,u,deru,cex,cey,cez,cflown)
c          write(*,*)'after elf_emis' 
      
c          write(*,*)'z,r,phi',z,r,phi
c          write(*,*)'cex,cey,cez',cex,cey,cez
c          write(*,*)'cflown',cflown,'vgroup',vgroup
c          write(*,*)'cnpar,cnper,cnray',cnpar,cnper,cnray

c---------the data for collisional absorption
          iabsorp_collisional_l=iabsorp_collisional
          coll_mult_l=coll_mult
          tempe_l=tempe(z,r,phi,1) ! kev
          dense_l=dense(z,r,phi,1) 
          zeff_l=zeff(z,r,phi) 

c          if (n_elt.ne.nrayelt_emis) then
c            win_sn_nc(n_elt,ifreq)=win_sn(n_elt)         !In_sn
c            win_0_nc(n_elt,ifreq)=win_0(in_elt)           !In_0
c            wi_0sn_nc(n_elt,ifreq)=wi_0sn(n_elt)         |I_0_sn  
c            wtaun_em_nc(n_elt,ifreq)=wtaun_em(n_elt)
c          endif

c---------calculate emission coefficient
c         emission absorptivity al_emis (1/cm) 
c         emissivity j_emis (erg*sec/cm**3)  and
c         radiation  temperature temprad

c---------the data for collisional absorption
          iabsorp_collisional_l=iabsorp_collisional
          coll_mult_l=coll_mult
          tempe_l=tempe(z,r,phi,1) ! kev
          dense_l=dense(z,r,phi,1) 
          zeff_l=zeff(z,r,phi) 
c          write(*,*)'before calc_emis_coef_spectrum'
          call calc_emis_coef_spectrum(xe,-ye,T_kev,cnpar,cnper,cnray,
     +    n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,
     +    phi,
     +    cex,cey,cez,cflown,vgroup,frqncy,
     +    iabsorp_collisional,coll_mult,tempe_l,dense_l,zeff_l,
     +    al_emis_loc,j_emis_loc,temp_rad_em_loc,
     +    jx_kin,kinetic_energy_kev,j_emis_x)
        
c          write(*,*)'after calc_emis_coef_spectrum'  
       
          do j=1,jx_kin
             wj_emis_x_nc(n_elt,j,ifreq)=j_emis_x(j)          
           enddo
           
        enddo !n_elt

c-------input data for calc_emission_spectrum
        do n_elt=1,nrayelt_emis
          do j=1,jx_kin
             wj_emis_x(n_elt,j)=wj_emis_x_nc(n_elt,j,ifreq)
          enddo
        enddo

c        write(*,*)'before calc_emission_spectrum'
     
        call calc_emission_spectrum(nrelta,wsn,
     &  wal_emis,wj_emis,
     &  wnray,win_sn,win_0, wi_0(iray_loc,ifreq),
        !YuP[2026-01-14] Corrected: input argument i_0 must be a scalar. 
     &  wi_0sn,wtaun_em,nrayelt_emis,
     &  transm_ox,nrayelt_o_cutoff_emis,
     &  jx_kin,kinetic_energy_kev,
     &  wj_emis_x,in_sn_x,in_0_x,i_0_x,i_0sn_x,
     &  work_j_emis_x_n,work_j_emis_x_np1,
     &  work_in_sn_x_n,work_in_0_x_n) 

c        write(*,*)'after calc_emission_spectrum'

c-------write output spectrums to netcdf file for centrum rays
        do j=1,jx_kin
          do n_elt=1,nrayelt_emis
            win_sn_x_nc(n_elt,j,ifreq)=in_sn_x(n_elt,j)
            win_0_x_nc(n_elt,j,ifreq) =in_0_x(n_elt,j)
            wi_0sn_x_nc(n_elt,j,ifreq) =i_0sn_x(n_elt,j)
          enddo !n_elt
          wi_0_x_nc(j,ifreq) =i_0_x(j)
        enddo !j

c-------test
        test_i_0=0.d0
        do j=1,jx_kin
          test_i_0=test_i_0+i_0_x(j)
        enddo
        iray_loc_=1
        write(*,*)'ifreq,wi_0(iray_loc_,ifreq),test_i_0',
     &             ifreq,wi_0(iray_loc_,ifreq),test_i_0

      endif ! iray_start_point.eq.1

 10   continue 
      return
      end
      

c-----------------------------------------------------------------------------
      subroutine calc_emis_coef_spectrum(xe,ye,T_kev,cnpar,
     +cnper,cnray,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,phi,
     +cex,cey,cez,cflown,vgroup,frqncy,
     +iabsorp_collisional,coll_mult,tempe,dense,zeff,
     +al_emis,j_emis,temprad,
     +jx_kin,kinetic_energy_kev,j_emis_x)
c---- calcultes the emission absorptivity al_emis (1/cm) 
c     emissivity j_emis (erg*sec/cm**3)  and
c     radiation  temperature temprad

      implicit none
cSm061213
      include 'param.i'
c-----input
      real*8 xe,ye,T_kev,cnpar,cnper,cnray,z,r,phi
c     cnray  is the ray refractive index N_ray
      integer n_relt_harm1, n_relt_harm2,n_relt_intgr,i_fkin
      complex*16 cex,cey,cez
c     wave energy flux for the normalized field |E|=1
c     cflown=(vectorB*vector_B^*+vector_E*d(omega*tensor_eps_herm)/dw*vector_E)
      complex*16 cflown  
      
      real*8 frqncy   !GhZ
c     the data for collisional absorption
      integer iabsorp_collisional ! =0 do not calculate collisional absorption
                                  ! =1 to calculte collisonal absorption
      real*8
     & vgroup,     !in clight
     & coll_mult,  !multiplier for coll absorption
     & tempe,      !electron temperature  (keV)  
     & dense,      !electron density   (10*13/cm**3	 
     & zeff        !plasma charge
cSm061213
      integer jx_kin              !number of grid points in momentum  
      real*8 kinetic_energy_kev(jx_kin) ! kinetic energy grid used
c                                      ! for emission spectrum 
c-----external
      real*8 rrind     
c-----output  
      real*8 al_emis,j_emis,temprad
cSm061213
      real*8 j_emis_x(jx_kin)  !(erg*sec/sm**3)
c-----local        
      complex*16 aK(3,3),fluctcur(3,3),ce(3)
      real*8
     & v_gr_cm_sec,       !cm/sec
     & omega,             !2*pi*f, here f in [HZ]
     & omegpedce,omegdce, !omega_pe/omega_ce and omega/omega_ce
     & cn,                !refructive index
     & clight,pi  
cSm061213
cSAP090808 jx_kin_a -> jx_kin    
c      complex*16 aK_x(3,3,jx_kin_a),fluctcur_x(3,3,jx_kin_a)
      complex*16, dimension(1:3,1:3,1:jx_kin) :: K_x,fluctcur_x
ctest
      integer j
 
c      write(*,*)'calc_emis_coef_spectrum bef emis_tens'

      call emis_tens_spectrum(xe,ye,T_kev,cnpar,cnper,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,phi,
     +aK,fluctcur,
cSm051213
     +jx_kin,kinetic_energy_kev,fluctcur_x)

      ce(1)=cex
      ce(2)=cey
      ce(3)=cez
      clight=2.99792458d10          !cm/sec
      v_gr_cm_sec=clight*vgroup     !cm/sec
      pi=4.d0*datan(1.d0)
      omega=2.d0*pi*frqncy*1.d9     !1/esc

c      write(*,*)'in calc_emis_coef_spectrum'

      call emiss_coef_spectrum(aK,fluctcur,ce,cflown,
     +v_gr_cm_sec,omega,cnray,
     +iabsorp_collisional,coll_mult,tempe,dense,zeff,
     +al_emis,j_emis,
cSm061213
     +fluctcur_x,
     +jx_kin,kinetic_energy_kev,j_emis_x)  
      
c      write(*,*)'after calc_emis_coef_spectrum'
c      do j=1,jx_kin
c        write(*,*)'j,j_emis_x(j)',j,j_emis_x(j) 
c      enddo

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

      subroutine emiss_coef_spectrum(aK,fluctcur,ce,flux,v_group,omega,
     +nr,
     +iabsorp_collisional,coll_mult,tempe,dense,zeff,
     +al_emis,j_emis,
cSm061213
     +fluctcur_x, 
     +jx_kin,kinetic_energy_kev,j_emis_x)

c     al_emis (1/cm)
c     j_emis  (erg*sec/cm**3)
c      
c-----calculates the emission coefficient j_emis(erg*sec/cm**3)
c     energy spectrum of emission coefficient j_emis_x(erg*sec/cm**3)
c     and absorption al_emis (1/cm) 
c     
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
cSm061213
      include 'param.i'
c-----input
      complex*16
     & aK(3,3),        ! anti-hermitian relativistic tensor
     & fluctcur(3,3),  ! correltion tensor G for the fluctuating current(erg)
     & ce(3),          ! wave (normalized) electric field in Stix frame
c     wave energy flux for the normalized field |E|=1
c     flux=(vectorB*vector_B^*+vector_E*d(omega*tensor_eps_herm)/dw*vector_E)
     & flux            ! wave energy flux for the normalized field |E|=1

      real*8
     & omega,      ! wave angle frequency (1/sec) =2*pi*f
     & v_group,    ! wave group velocity [cm/sec]
     & nr          ! ray refructive index

      integer iabsorp_collisional ! =0 do not calculate collisional absorption
                                  ! =1 to calculte collisonal absorption
      real*8
     & coll_mult,                 ! multiplier for coll absorption
     & tempe,                     !electron temperature  (keV)  
     & dense,                     !electron density   (10*13/cm**3	 
     & zeff                       !plasma charge
cSm061213
      integer jx_kin              !number of kinetic energy grid points        
      complex*16 fluctcur_x(3,3,jx_kin)  !energy spectrum of fluctcur
      real*8 kinetic_energy_kev(jx_kin)  ! kinetic energy grid used
                                         ! for emission spectrum
c-----output
      real*8
     & al_emis,  !(1/sm)
     & j_emis,   !(erg*sec/sm**3)
cSm061213
     & j_emis_x(jx_kin)  !(erg*sec/sm**3)
c-----local
      complex*16 cec(3),cal_emis,cj_emis
      integer i,j,j1
      real*8 s,pi,clight,absflux,p,
     &arg,cln,r_nu_ei,al_emis_collisional

cSm061213
      complex*16, dimension(1:jx_kin) :: cj_emis_x
c-----for test
      real*8 test_j_emis
           
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

cSm061213
      do j1=1,jx_kin
        cj_emis_x(j1)=dcmplx(0.0d00,0.0d00)
      enddo     


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
cSm061213
            do j1=1,jx_kin            
              cj_emis_x(j1) =cj_emis_x(j1)+
     &                      cec(i)*fluctcur_x(i,j,j1)*ce(j)
            enddo
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
cSm061213
      do j1=1,jx_kin
         cj_emis_x(j1) =cj_emis_x(j1)*pi*(nr*omega/clight)**2/s
      enddo

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

cSm061213
      do j1=1,jx_kin
        p=cj_emis_x(j1)*dconjg(cj_emis_x(j1))
        j_emis_x(j1)=dsqrt(p)*omega
        j_emis_x(j1)=j_emis_x(j1)*4.d0*pi      
      enddo

ctest----------------------------------
c      write(*,*)'in  emiss_coef_spectrum'
      test_j_emis=0.d0
      do j1=1,jx_kin
         test_j_emis=test_j_emis+j_emis_x(j1)
c         write(*,*)'j1,j_emis_x(j1)',j1,j_emis_x(j1)
      enddo
c      write(*,*)'test_j_emis, j_emis',test_j_emis, j_emis
cend test ---------------------------

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



      subroutine emis_tens_spectrum(X,Y,T_kev,nll_in,np_in,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +i_fkin,r,z,phi,
     +aK,fluctcur,
cSm051213
     +jx_kin,kinetic_energy_kev,fluctcur_x)
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
c     OUTPUT:
c      aK(3,3):  the nine components of the anti-hermition
c                part dielectric relativistic tensor
c                evaluated at (X,Y,Te,nll,np_out,n)
c      
c      fluctcur(3,3) the correlation tensor fluctcur=G (egr)  
c                    for the fluctuating current density
c      
c     jx_kin !the number of kinetic energy grid  points
c     kinetic_energy_kev(jx_kin) the kinetic energy grid used
c                                  for emission spectrum 
      implicit none
cSm061213
      include 'param.i'

c     input 
      real*8 X,Y,T_kev,nll_in,np_in
      integer n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin
      real*8 r,z,phi
cSm061213
      integer jx_kin
      real*8 kinetic_energy_kev(jx_kin) 

      complex*16 aK(3,3),fluctcur(3,3)

cSm061213
      complex*16 fluctcur_x(3,3,jx_kin) !velosity spectrum of fluctcur

c     local
      real*8 c,mass_e,k_1_kev,nll,np_out,nlls,p_perp0,
     &theta,pi,
     *dens,xpidens,xpidensmcs,p,psi,rho
      complex*16 integral(3,3),fluctcur_n(3,3)  

c061213
      complex*16, dimension(1:3,1:3,1:jx_kin) :: fluctcur_n_x

      integer jn,ires,jn_negative
cSm061213
      integer j,k1,k2

c     external zeroK, npnllmin, ec_cond, intgr_rl,fdens_fdist
      real*8 fdens_fdist,fpsi,rhopsi,densrho

c-----test
      complex*16 test_fluctcur(3,3),test_fluctcur_sum(3,3),
     &test_fluctcur_n(3,3)



      call zeroK(test_fluctcur)    
c--------------------------------
 
      c =2.99792458d10          !speed of light     (cm/sec)
      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)      

      theta=mass_e*c**2/(k_1_kev*T_kev)
      pi=4.d0*datan(1.d0)

      call npnllmin(nll_in,np_in,nll,np_out)
      nlls = nll**2

c      write(*,*)'in emis_tens_spectrum'

c-----initialization aK=0
      call zeroK(aK)
      call zeroK(fluctcur)    

c for test
      call zeroK(test_fluctcur_sum)
cendtet
    
cSm061213
      do k1=1,3
        do k2=1,3 
          do j=1,jx_kin
             fluctcur_x(k1,k2,j)=dcmplx(0.d0,0.d0)
          enddo
        enddo
      enddo

        
c-----the loop over the cyclotron harmonics
      
      do jn=n_relt_harm1,n_relt_harm2               
        jn_negative=-jn
c-------control the EC resonance condition
cSm060315         
c        call ec_cond(jn,Y,nll,ires,p_perp0)
        call ec_cond(jn_negative,Y,nll,ires,p_perp0)

c        write(*,*)'emis_tens_spectrum jn,ires,p_perp0',jn,ires,p_perp0
        
        if(ires.eq.0) then
c---------no resonace
          goto 10
        endif
 
        call intgr_emis_spectrum(jn_negative,nll,np_out,Y,theta,ires,
     +  p_perp0,
     +  n_relt_intgr,
     +  i_fkin,r,z,phi,
     +  integral,fluctcur_n,
cSm061213
     +   jx_kin,kinetic_energy_kev,fluctcur_n_x)

c       write(*,*)'emis_tens_spectrum after integr_emis integral',
c     &  integral
c        write(*,*)'emis_tens_spectrum after integr_emis fluctcur_n',
c     &  fluctcur_n

c-------test
c        call zeroK(test_fluctcur_n)
c      write(*,*)'in emis_tens_spectrum afterintgr_emis_spectrum jn',jn
c        do j=1,jx_kin
c          do k1=1,3
c            do k2=1,3
c               test_fluctcur_n(k1,k2)=test_fluctcur_n(k1,k2)+
c     &                                fluctcur_n_x(k1,k2,j)
c            enddo
c           enddo
c        enddo

c        do k1=1,3
c           do k2=1,3
c            write(*,*)'k1,k2',k1,k2
c            write(*,*)'fluctcur_n(k1,k2),test_fluctcur_n(k1,k2)',
c     &                 fluctcur_n(k1,k2),test_fluctcur_n(k1,k2)     
c           enddo
c        enddo
c-------end_test    



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

cSm061213
        do j=1,jx_kin 
          fluctcur_x(1,1,j)=fluctcur_x(1,1,j)+fluctcur_n_x(1,1,j)
          fluctcur_x(1,2,j)=fluctcur_x(1,2,j)+fluctcur_n_x(1,2,j)
          fluctcur_x(1,3,j)=fluctcur_x(1,3,j)+fluctcur_n_x(1,3,j)
          fluctcur_x(2,2,j)=fluctcur_x(2,2,j)+fluctcur_n_x(2,2,j)
          fluctcur_x(2,3,j)=fluctcur_x(2,3,j)+fluctcur_n_x(2,3,j)
          fluctcur_x(3,3,j)=fluctcur_x(3,3,j)+fluctcur_n_x(3,3,j)
        enddo

ctest   
c        call zeroK(test_fluctcur_sum)    
c        do j=1,jx_kin
c         do k1=1,3
c            do k2=1,3
c               test_fluctcur_sum(k1,k2)=test_fluctcur_sum(k1,k2)+
c     &                                 fluctcur_n_x(k1,k2,j)
c               test_fluctcur_sum(k1,k2)=test_fluctcur_sum(k1,k2)+
c     &                                 fluctcur_x(k1,k2,j)
c            enddo
c         enddo
c        enddo

c        do k1=1,3
c          do k2=1,3
c            write(*,*)'k1,k2,test_fluctcur_sum(k1,k2),fluctcur(k1,k2)'
c     &                ,k1,k2,test_fluctcur_sum(k1,k2),fluctcur(k1,k2)
c          enddo
c        enddo
cendtest

10      continue      
      enddo   !jn


ctest   
c        call zeroK(test_fluctcur)    
c        do j=1,jx_kin
c         do k1=1,3
c            do k2=1,3
c               test_fluctcur(k1,k2)=test_fluctcur(k1,k2)+
c     &                                 fluctcur_n_x(k1,k2,j)
c            enddo
c         enddo
c        enddo

c        do k1=1,3
c          do k2=1,3
c            write(*,*)'k1,k2,test_fluctcur(k1,k2),fluctcur(k1,k2)'
c     &                ,k1,k2,test_fluctcur(k1,k2),fluctcur(k1,k2)
c          enddo
c        enddo
cendtest


      if (i_fkin.eq.0) then       
        dens=1.d0
      else      
        dens=fdens_fdist(r,z,phi)
c        write(*,*)'fdens_fdis dens',dens
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


      xpidensmcs=X*pi/(2.d0*pi)**5*(mass_e*c**2)/dens

cSm061213
      do j=1,jx_kin           
        fluctcur_x(1,1,j)=fluctcur_x(1,1,j)*xpidensmcs
        fluctcur_x(1,2,j)=fluctcur_x(1,2,j)*xpidensmcs
        fluctcur_x(1,3,j)=fluctcur_x(1,3,j)*xpidensmcs
        fluctcur_x(2,2,j)=fluctcur_x(2,2,j)*xpidensmcs
        fluctcur_x(2,3,j)=fluctcur_x(2,3,j)*xpidensmcs
        fluctcur_x(3,3,j)=fluctcur_x(3,3,j)*xpidensmcs
      
        fluctcur_x(2,1,j)=-fluctcur_x(1,2,j)
        fluctcur_x(3,1,j)= fluctcur_x(1,3,j)
        fluctcur_x(3,2,j)=-fluctcur_x(2,3,j)
      enddo
    

      fluctcur(1,1)=fluctcur(1,1)*xpidensmcs
      fluctcur(1,2)=fluctcur(1,2)*xpidensmcs
      fluctcur(1,3)=fluctcur(1,3)*xpidensmcs
      fluctcur(2,2)=fluctcur(2,2)*xpidensmcs
      fluctcur(2,3)=fluctcur(2,3)*xpidensmcs
      fluctcur(3,3)=fluctcur(3,3)*xpidensmcs
      
      fluctcur(2,1)=-fluctcur(1,2)
      fluctcur(3,1)= fluctcur(1,3)
      fluctcur(3,2)=-fluctcur(2,3)
      
c-----test
c      do j=1,jx_kin
c        do k1=1,3
c          do k2=1,3 
c            test_fluctcur(k1,k2)=test_fluctcur(k1,k2)+
c     &                             fluctcur_x(k1,k2,j)
c          enddo
c        enddo         
c      enddo

c      write(*,*)'in emis_tens_spectrum'
c      do k1=1,3
c         do k2=1,3
c            write(*,*)'k1,k2,test_fluctcur(k1,k2),fluctcur(k1,k2)',
c     &                 k1,k2,test_fluctcur(k1,k2),fluctcur(k1,k2)
c         enddo
c      enddo
c-----endtest

      return
      end



      subroutine intgr_emis_spectrum(n,nll,nper,Y,theta,ires,p_perp0,
     +n_relt_intgr,
     +i_fkin,r,z,phi,
     +integral,fluctcur_n,
cSm061213
     + jx_kin,kinetic_energy_kev,fluctcur_n_x)
c-----------------------------------------------------------------
c     calculates the matrix: double complex integrals
c     for the relativistic electron plasma 
c     I(n)^=integral(0<=p_perp<=p_perp0)g_n, g_n={sum_k=1,2,G_nk(p_perp)
c     fluctcar(n)^=integral(0<=p_perp<=p_perp0)g_fluctcur_n.
c               g_fluctcur__n={sum_k=1,2,g_fluctcur_nk(p_perp)
c     for the EC harmonic with number 'n'
c----------------------------------------------------------------- 
c     input
c       n        - EC harmonic number
c       nll      - N_parallel
c       nper       - N_perpendicular
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
c       jx_kin       the number of the kinetic energy grid points 
c       kinetic_energy_kev(jx_kin)! kinetic energy grid used
                                    ! for emission spectrum
c     output
c       integral(3,3) double complex integral from G_nk over 0<=p_perp=<p_perpmax
c       fluctcur(3,3) double complec integral from fluct_cur_nk over 0<=p_perp=<p_perpmax
c
c       fluctcur_n_x(3,3,jx_kin_a) is the velosity spectrum of fluctcur
c----------------------------------------------------------------    
c      The integration method is specified by the variable
c      i_resonance_curve_integration_method.
c      Now this variable is set inside this subroutine.
c----------------------------------------------------------------
      implicit none
cSm061213
      include 'param.i'
c     input
      real*8 nll,nper,Y,theta,p_perp0
      integer n,ires,n_relt_intgr,i_fkin 
      real*8 r,z,phi,
     &vnormloc,massloc ! cm/sec, g
      COMMON /dskin1/vnormloc,massloc
cSm061213
      integer jx_kin
      real*8 kinetic_energy_kev(jx_kin)  ! kinetic energy grid used
                                                   ! for emission spectrum
c-----output
      complex*16 integral(3,3),
     &fluctcur_n(3,3),
cSm061213
     &fluctcur_n_x(3,3,jx_kin)
cSm061213
c     local
cSm060725
      real*8 eps,p_permax,h,p_perp,p,p_t,clight,p_int,
c     double precision eps,p_permax,h,p_perp,p,p_t,jmax,clight,p_int,
     &vmax_d_vt
      complex*16 i,g_n(3,3),g_fluctcur_n(3,3)
cSm060725
      integer j,jmax
cSm061213
      integer k1,k2

c-----for integration along ellipse by angle
      real*8 rme,rtem0,vper,xint,cs,sn,thet1,thet2,p_par,
     & p_par_min,p_par_pl,p_perp_min,p_perp_pl,v0dc,vmax1dc,vpar0dc,
     +vpar1dc,vpar2dc,dth,vt0,cdvt,tem0,dvper,vmaxdt

      integer i_resonance_curve_integration_method

c     external
c     ec_cond, zeroK, g_n 
      real*8 temperho

      integer k_root,k,jx
      real*8 kinetic_energy_2(2),beta,beta_p
      complex*16 g_fluctcur_n_2(3,3,2)
                                         

c_for test
      complex*16 test_fluctcur_n(3,3)

      i = ( 0.0d0,1.0d0)        !imaginary number     
      jmax=n_relt_intgr

cSm060310
      if(ires.eq.0)goto10 ! no resonance
    
      vmax_d_vt=10.d0
      call p_perp_max_calc(i_fkin,theta,n,Y,nll,vmax_d_vt,
     &vnormloc,p_permax,ires)

c      write(*,*)'intgr_emis_spectrum ires,p_permax,vmax_d_vt,vnormloc',
c     &ires,p_permax,vmax_d_vt,vnormloc

      if(ires.eq.4) then
c--------the resonance curve is outside the grid
         call zeroK(integral)
         call zeroK(fluctcur_n)

         do j=1,jx_kin
          do k1=1,3
            do k2=1,3
               fluctcur_n_x(k1,k2,j)=dcmplx(0.d0,0.d0)
            enddo
          enddo
         enddo

         goto 10
      else
         h=p_permax/(n_relt_intgr)   ! step of integration over p_perp
      endif

c     calculations of the integrals over p_perp
      call zeroK(integral)
      call zeroK(fluctcur_n)

cSm061213
      do j=1,jx_kin
         do k1=1,3
            do k2=1,3
               fluctcur_n_x(k1,k2,j)=dcmplx(0.d0,0.d0)
            enddo
         enddo
      enddo


c------------------
c     i_resonance_curve_interation_method
c
c     i_resonance_curve_interation_method=1 angle integration
c     now it works for ellipse case only                             c
c     i_resonance_curve_interation_method=2 rectangle formula
c     for p_perp integraion
c     
c     i_resonance_curve_interation_method=3 trapezoidal formula
c     for p_perp integraion
c     Adaptive Simpson method can not be used for spectrum calculations 
c------------------------------------------------------------ 
c      i_resonance_curve_integration_method=1
c      i_resonance_curve_integration_method=2
      i_resonance_curve_integration_method=3
     
      goto (1,2,3) i_resonance_curve_integration_method

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

         call g_n_emis_theta_spectrum(p_perp,p_par,y,nll,nper,theta,n,
     &   i_fkin,
     &   r,z,phi,g_n,g_fluctcur_n)
c     &   k_root,kinetic_energy_2,g_fluctcur_n_2)
          
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
cSm061213
         fluctcur_n_x(1,1,j)=h*g_fluctcur_n(1,1)
         fluctcur_n_x(1,2,j)=h*g_fluctcur_n(1,2)
         fluctcur_n_x(1,3,j)=h*g_fluctcur_n(1,3)
         fluctcur_n_x(2,2,j)=h*g_fluctcur_n(2,2)
         fluctcur_n_x(2,3,j)=h*g_fluctcur_n(2,3)
         fluctcur_n_x(3,3,j)=h*g_fluctcur_n(3,3)

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

        call g_n_emis_spectrum(p_perp,y,nll,nper,theta,n,
     &  i_fkin,r,z,phi,
     &  g_n,g_fluctcur_n,
     &  k_root,kinetic_energy_2,g_fluctcur_n_2)

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

cSm061213

        fluctcur_n_x(1,1,j)=g_fluctcur_n(1,1)*h
        fluctcur_n_x(1,2,j)=g_fluctcur_n(1,2)*h
        fluctcur_n_x(1,3,j)=g_fluctcur_n(1,3)*h
        fluctcur_n_x(2,2,j)=g_fluctcur_n(2,2)*h
        fluctcur_n_x(2,3,j)=g_fluctcur_n(2,3)*h
        fluctcur_n_x(3,3,j)=g_fluctcur_n(3,3)*h
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
c      write(*,*)'intgr_emis_spectrum 3: jmax=',jmax
      do j=0,jmax

c         write(*,*)'intgr_emis_spectrum: j',j

         p_int=1.d0
         if((j.eq.1).or.(j.eq.jmax)) p_int=0.5d0
         p_perp=h*j
         if(p_perp.lt.1.d-3) p_perp=1.d-3
         if(j.eq.jmax) p_perp=p_perp-1.d-3


c         write(*,*)'intgr_emis_spectrum before g_n_emis_spectrum p_perp'
c     &   ,p_perp

         call g_n_emis_spectrum(p_perp,y,nll,nper,theta,n,i_fkin,
     &   r,z,phi,
     &   g_n,g_fluctcur_n,
     &   k_root,kinetic_energy_2,g_fluctcur_n_2)

ctest                 
c         write(*,*)'intgr_emis_spectrum after g_n_emis_spectrum k_root'
c     &   ,k_root

c         do k1=1,3
c            do k2=1,3
c               write(*,*)'k1,k2,g_fluctcur_n(k1,k2)',
c     &                    k1,k2,g_fluctcur_n(k1,k2)
c               do k=1,k_root
c                  write(*,*)'k,g_fluctcur_n_2(k1,k2,k)',
c     &                       k,g_fluctcur_n_2(k1,k2,k)
c               enddo
c            enddo
c         enddo
cendtest

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

cSm061213      
c------------------------------------------------------------
c        find the point at the kinetic enegry grid jx
c        and put p_int*g_fluctcur_n*h to   fluctcur_n_x(jx)
c---------------------------------------------------------------
c         write(*,*)'intgr_emis_spectrum: j,k_root',j,k_root

         do k=1,k_root 

cSm070105
           call search_grid_point(kinetic_energy_kev,jx_kin,
     &     kinetic_energy_2(k),jx)


c           write(*,*)'k,kinetic_energy_2(k),jx,kinetic_energy_kev(jx)',
c     &                k,kinetic_energy_2(k),jx,kinetic_energy_kev(jx)       

           fluctcur_n_x(1,1,jx)=fluctcur_n_x(1,1,jx)+
     &                          p_int*g_fluctcur_n_2(1,1,k)*h
           fluctcur_n_x(1,2,jx)=fluctcur_n_x(1,2,jx)+
     &                          p_int*g_fluctcur_n_2(1,2,k)*h
           fluctcur_n_x(1,3,jx)=fluctcur_n_x(1,3,jx)+
     &                          p_int*g_fluctcur_n_2(1,3,k)*h
           fluctcur_n_x(2,2,jx)=fluctcur_n_x(2,2,jx)+
     &                          p_int*g_fluctcur_n_2(2,2,k)*h
           fluctcur_n_x(2,3,jx)=fluctcur_n_x(2,3,jx)+
     &                          p_int*g_fluctcur_n_2(2,3,k)*h
           fluctcur_n_x(3,3,jx)=fluctcur_n_x(3,3,jx)+
     &                          p_int*g_fluctcur_n_2(3,3,k)*h
 
c------------------------------------------------------------------------
c           to get more smooth energy spectrum
c 
c           call search_grid_point_smooth(kinetic_energy_kev,jx_kin,
c     &     kinetic_energy_2(k),jx,beta,beta_p)
           
c           fluctcur_n_x(1,1,jx)=fluctcur_n_x(1,1,jx)+
c     &                          p_int*g_fluctcur_n_2(1,1,k)*h*beta
c           fluctcur_n_x(1,2,jx)=fluctcur_n_x(1,2,jx)+
c     &                          p_int*g_fluctcur_n_2(1,2,k)*h*beta
c           fluctcur_n_x(1,3,jx)=fluctcur_n_x(1,3,jx)+
c     &                          p_int*g_fluctcur_n_2(1,3,k)*h*beta
c           fluctcur_n_x(2,2,jx)=fluctcur_n_x(2,2,jx)+
c     &                          p_int*g_fluctcur_n_2(2,2,k)*h*beta
c           fluctcur_n_x(2,3,jx)=fluctcur_n_x(2,3,jx)+
c     &                          p_int*g_fluctcur_n_2(2,3,k)*h*beta
c           fluctcur_n_x(3,3,jx)=fluctcur_n_x(3,3,jx)+
c     &                          p_int*g_fluctcur_n_2(3,3,k)*h*beta

c           if (jx.lt.jx_kin) then
c              fluctcur_n_x(1,1,jx+1)=fluctcur_n_x(1,1,jx+1)+
c     &                          p_int*g_fluctcur_n_2(1,1,k)*h*beta_p
c              fluctcur_n_x(1,2,jx+1)=fluctcur_n_x(1,2,jx+1)+
c     &                          p_int*g_fluctcur_n_2(1,2,k)*h*beta_p
c              fluctcur_n_x(1,3,jx+1)=fluctcur_n_x(1,3,jx+1)+
c     &                          p_int*g_fluctcur_n_2(1,3,k)*h*beta_p
c              fluctcur_n_x(2,2,jx+1)=fluctcur_n_x(2,2,jx+1)+
c     &                          p_int*g_fluctcur_n_2(2,2,k)*h*beta_p
c              fluctcur_n_x(2,3,jx+1)=fluctcur_n_x(2,3,jx+1)+
c     &                          p_int*g_fluctcur_n_2(2,3,k)*h*beta_p
c              fluctcur_n_x(3,3,jx+1)=fluctcur_n_x(3,3,jx+1)+
c     &                          p_int*g_fluctcur_n_2(3,3,k)*h*beta_p
c           endif
c----------------------------------------------------------------------------
         enddo !k
       
c         do k1=1,3
c           do k2=k1,3
c             write(*,*)'k1,k2,g_fluctcur_n(k1,k2)',
c     &                  k1,k2,g_fluctcur_n(k1,k2)
c             do k=1,k_root
c               write(*,*)'k,g_fluctcur_n_2(k1,k2,k)',
c     &                    k,g_fluctcur_n_2(k1,k2,k)
c             enddo
c           enddo
c         enddo

      enddo !j 
 
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
c     integration along the resonance curve using trapezoidal
c     formula
c     end
c--------------------------------------------------------------------
 20   continue

      integral(2,1)=-integral(2,1)
      integral(3,1)= integral(3,1)
      integral(3,2)=-integral(3,2)
   
      fluctcur_n(2,1)=-fluctcur_n(2,1)
      fluctcur_n(3,1)= fluctcur_n(3,1)
      fluctcur_n(3,2)=-fluctcur_n(3,2)
cSm061213
      do j=1,jx_kin
        fluctcur_n_x(2,1,j)=-fluctcur_n_x(2,1,j)
        fluctcur_n_x(3,1,j)= fluctcur_n_x(3,1,j)
        fluctcur_n_x(3,2,j)=-fluctcur_n_x(3,2,j)
      enddo

c-----test
c      call zeroK(test_fluctcur_n)

c      do j=1,jx_kin
c         do k1=1,3
c            do k2=1,3
c               test_fluctcur_n(k1,k2)=test_fluctcur_n(k1,k2)+
c     &                                fluctcur_n_x(k1,k2,j)
c            enddo
c         enddo
c      enddo

c      do k1=1,3
c         do k2=1,3
c            write(*,*)'k1,k2',k1,k2
c            write(*,*)'fluctcur_n(k1,k2),test_fluctcur_n(k1,k2)',
c     &                 fluctcur_n(k1,k2),test_fluctcur_n(k1,k2)     
c         enddo
c      enddo
c-----end_test    

  
 10   continue
      return
      end   
  

      subroutine emis_bin_sn_spectrum(n,al_emis_n,al_emis_n1,
     +j_emis_n,j_emis_n1,
     +nr_n,nr_n1,s_n,s_n1,in_sn,
cSm061215
     +j_emis_n_x,j_emis_n1_x,
     +jx_kin,in_sn_x)  
c-----calculates the emission at the bin boundary s=s_n from one bin s_n<s<s_n+1
      implicit none
cSm061215
      include 'param.i'
c-----input
      integer n ! the number of ray bin, s_n(at n=1)=0,It should be: 1=<n<nrelt-1 
      real*8 al_em is_n,j_emis_n,    !absorption and emission at s=s_n
     +al_emis_n1,j_emis_n1,  !absorption and emission at s=s_(n+1)
     +nr_n,nr_n1,            !N_ray rayrefractive index at s=s_n and s_(n+1)
     +s_n,s_n1               !the lengths along the ray (cm) s_n1=s_n+1
cSm061215
      integer jx_kin         !number of momentum grid points    
      real*8 j_emis_n_x(jx_kin),j_emis_n1_x(jx_kin)

c     al_emis [1/cm]
c     j_emis  [egr/(cm**3*sec)]
c-----output
      real*8 in_sn !emission I_n at the detector side of nth bin at s=s_n

cSm061215
      real*8 in_sn_x(jx_kin) 
      
c-----local
      real*8 j_emis_np,           !j_emis(s_n+0.5)
     +al_emis_np,                 !al_emis(s_n+0.5)
     +nrs_np,                     !N_r**2(s_n+0.5)
     +p,p1

cSm061215
      integer j
  
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
      endif

cSm061215
      do j=1,jx_kin
         j_emis_np=0.5d0*(j_emis_n_x(j)+j_emis_n1_x(j)) 
         p1=j_emis_np*nr_n**2
         if ((p1.eq.0.d0).or.(p.eq.0.d0)) then
           in_sn_x(j)=0.d0
         else
           in_sn_x(j)=p1/p*(1.d0-dexp(-al_emis_np*(s_n1-s_n)))
         endif
      enddo

      return
      end  

      subroutine emis_bin_s1_spectrum(n,nrayelt,wsn,wal_emis,wnray,
     +in_sn,
     +tau_n,in_0,nrayelt_o_cutoff,transm_ox,
cSm061215
     +jx_kin,in_sn_x,in_0_x)  
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
      real*8 wal_emis(*),    !absorption at s=s_n
     +wsn(*),  !the lengths s=sn  along the ray (cm) 
     +wnray(*),!N_ray ray refractive index along the ray 
     +in_sn,   !emission I_n at the detector side of nth bin at s=s_n
     *transm_ox !transmission coefficient for OX conversion
c     wal_emis [1/cm]
cSm061215
      include 'param.i'
      integer jx_kin   ! the number of momentum grid points
      real*8 in_sn_x(jx_kin) !momentum spectrum of in_sn
c-----output
      real*8 in_0, !emission I_n at the plasma boundary from nth bin at s=s_1=0
     &tau_n,       !tau_n from n_th bin
cSm061215
     &in_0_x(jx_kin) !momentum spectrum of in_0
c-----local
      real*8 
     +al_emis_np                 !al_emis(s_n+0.5)                  
      integer m

cSm061215
      integer j


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
    
cSm061215
      do j=1,jx_kin
cSm070120
c         in_0_x(j)=(wnray(1)/wnray(n))**2*in_sn_x(j)*dexp(-tau_n)
         in_0_x(j)=(1.d0/wnray(n))**2*in_sn_x(j)*dexp(-tau_n)
         if ((n.ge.nrayelt_o_cutoff+1).and.(nrayelt_o_cutoff.gt.1)) then
          in_0_x(j)=in_0_x(j)*transm_ox
         endif    
      enddo

      return
      end  



      subroutine calc_emission_spectrum(nrelta,wsn,
     &wal_emis,wj_emis,
     &wnray,in_sn,in_0,i_0,i_0sn,wtaun_em,nrayelt_emis,
     &transm_ox,nrayelt_o_cutoff_emis,
cSm061215
cSAP090808 deklete jx_kin_a -> jx_kin
     &jx_kin,kinetic_energy_kev,
     &wj_emis_x,in_sn_x,in_0_x,i_0_x,i_0sn_x,
     &work_j_emis_x_n,work_j_emis_x_np1,
     &work_in_sn_x_n,work_in_0_x_n)        
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
      real*8
     +wal_emis(*),wj_emis(*),!absorption and emission at s=s_n
     +wsn(*),    !the length s=s_n along the ray (cm) 
     +wnray(*),  !N_ray along the ray
     &transm_ox  !transmission coefficient for OX conversion
     
c     wal_emis [1/cm]
c     wj_emis  [egr/(cm**3*sec)]

cSm061206
      integer jx_kin  !the number of the kinetic energy grid points
c      integer jx_kin_a!max number of the kinetic energy grid points
      real*8 kinetic_energy_kev(jx_kin)! kinetic energy grid used
                                       ! for emission spectrum
cSAP090808 jx_kin_a -> jx_kin
      real*8 wj_emis_x(nrelta,jx_kin),
     &work_j_emis_x_n(jx_kin),work_j_emis_x_np1(jx_kin),
     &work_in_sn_x_n(jx_kin),work_in_0_x_n(jx_kin) 

      integer nrayelt_emis ! number of the emission output points after
                           ! the addition of new points
      integer nrayelt_o_cutoff_emis !the number of the output point
                                    !(after addition new points)  
                                    !wsn(nrayelt_o_cutoff_emis)=s_O_mode 
                                    !along the ray (at O-mode part)
                                    !where the OX transmission
                                    !coefficient transm_ox has the maximal value

c-----output
      real*8
     +in_sn(*), ! emission I_n at the detector side of nth bin at s=s_n
     +in_0(*),  ! emission at the plasma boundary s=s_1 from one
                ! nth bin s_n<s<s_n+1
     +i_0, !emission I_0 at the plasma boundary from the ray at s=s(1)=0
     +i_0sn(*), ! sum{k=1,n}[in_0(k)]! emission at the plasma boundary
                ! from the parts 0<s<sn of the ray
     +wtaun_em(*)! tau_n from n-th bin for plotting 

cSm061215
      ! bellow are momentum spectrums
cSAP090808 jx_kin_a -> jx_kin
      real*8 in_sn_x(nrelta,jx_kin), ! emission I_n at the detector side
                                           ! of nth bin at s=s_n
     +in_0_x(nrelta,jx_kin),  ! emission at the plasma boundary s=s_1 from one
                              ! nth bin s_n<s<s_n+1
     +i_0_x(jx_kin),          ! emission I_0 at the plasma boundary from the
                              ! ray at s=s(1)=0
     +i_0sn_x(nrelta,jx_kin)  ! sum{k=1,n}[in_0(k)]! emission at the plasma boundary
                              ! from the parts 0<s<sn of the ray
     

c-----local
      integer n,n_midway,isplit,
     +n_checked,nrayelt_emis_new,k,i_ox_jump
     
cSm061215
      integer k1,j
     

cSm061219
      nrayelt=nrayelt_emis
      nrayelt_o_cutoff=nrayelt_o_cutoff_emis

c      write(*,*)'in calc_emission_spectrum nrayelt=', nrayelt

      if (nrayelt.gt.nrelta) then
         write(*,*)'in calc_emission_spectrum'
         write(*,*)'nrayelt.gt.nrelta nrayelt,nrelta',
     .   nrayelt,nrelta
         write(*,*)'it should be nrayelt=<nrelta'
         write(*,*)'increas nrelta in param.i'
         stop
      endif  

c------------------------------------------------------------------
c     initialization 
      call bcast(in_0,0.d0,nrelta)

cSm061215
      do k1=1,nrelta
         do j=1,jx_kin
           in_0_x(k1,j)=0.d0
         enddo
      enddo
c-------------------------------------------------------------------

      if((nrayelt_o_cutoff.gt.1).and.(nrayelt_o_cutoff.lt.nrayelt))then
        i_ox_jump=1 !to use OX conversion jump
      else
        i_ox_jump=0 !do not use OX conversion jump
      endif
   
      write(*,*)'emission i_ox_jump',i_ox_jump
      do n=1,nrayelt-1
         if((i_ox_jump.eq.1).and.(n.eq.nrayelt_o_cutoff)) then
           in_sn(n)=0.d0 !in OX jump bin area
           in_0(n)=0.d0  !from OX jump bin area

cSm061215
           do j=1,jx_kin
             in_sn_x(n,j)=0.d0 !in OX jump bin area
             in_0_x(n,j)=0.d0  !from OX jump bin area
           enddo

         else  
c----------calculate in_sn emission I_n at the detector side
c          of n th bin at s=s_n
 
cSm061215
           do j=1,jx_kin
             work_j_emis_x_n(j)=wj_emis_x(n,j)
             work_j_emis_x_np1(j)=wj_emis_x(n+1,j)
           enddo

           call emis_bin_sn_spectrum(n,wal_emis(n),wal_emis(n+1),
     &     wj_emis(n),
     +     wj_emis(n+1),wnray(n),wnray(n+1),wsn(n),wsn(n+1),in_sn(n),
cSm061215
     +     work_j_emis_x_n,work_j_emis_x_np1,
     +     jx_kin,work_in_sn_x_n) 
 
cSm061215
           do j=1,jx_kin
             in_sn_x(n,j)=work_in_sn_x_n(j)
           enddo

c         write(*,*)'calc_emission_spectrum aft emis_bin_sn n,in_sn(n)',
c     &     n,in_sn(n)  

c----------calculate in_0 the emission at the plasma boundary s=s_1
c          from one nth bin s_n<s<s_n+1
           call emis_bin_s1_spectrum(n,nrayelt,wsn,wal_emis,wnray,
     &     in_sn(n),
     +     wtaun_em(n),in_0(n),nrayelt_o_cutoff,transm_ox,
cSm061215
     &     jx_kin,work_in_sn_x_n,work_in_0_x_n) 
           
           do j=1,jx_kin
             in_sn_x(n,j)=work_in_sn_x_n(j)
             in_0_x(n,j)= work_in_0_x_n(j) 
           enddo
    
c      write(*,*)'calc_emission_spect emis_bin_s1 n,in_0(n),i_0sn(n)'
c     &     ,n,in_0(n),i_0sn(n)
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
cSm061215
      do j=1,jx_kin
        i_0_x(j)=0.d0  
      enddo

      do n=nrayelt_emis,1,-1
         i_0=i_0+in_0(n)
         i_0sn(n)=i_0
c         write(*,*)'n,in_0(n)',n,in_0(n)
cSm061215
         do j=1,jx_kin
           i_0_x(j)=i_0_x(j)+in_0_x(n,j)
           i_0sn_x(n,j)=i_0_x(j)
         enddo
      enddo  
         

ctest
c      do n=1,nrayelt_emis
c         write(*,*)'emission n,wal_emis(n),wj_emis(n)',
c     .   n,wal_emis(n),wj_emis(n)
c      enddo
cendtest


      return
      end

     
      subroutine g_n_emis_spectrum(p_perp,y,nll,np,theta,n,
     &i_fkin,r,z,phi,
     &g_n,g_fluctcur_n,
     &kmax,kinetic_energy_2,g_fluctcur_n_2)
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
c       kmax    - number of the roots N_par(N_perp) at the resonance
c                 curve. It can be kmax=0,1,2
c       kinetic_energy_2(kmax) The kinetic energy values of the point at
c                            at the resonance curve for the given p_perp [KeV]
c       g_fluctcur_n_2(3,3,kmax) g_fluctcur_n at calculated kinetic energys
c                                kinetic_energy_2(kmax)
c-------------------------------------------------------
      implicit none
c-----input        
      double precision p_perp,y,nll,np,theta,r,z,phi
      integer n,i_fkin
c-----external root_res, s_calcn, zeroK, unde_intgr
      double precision u_n,u_flcur_n    
c-----output
      double complex g_n(3,3),g_fluctcur_n(3,3),
     &g_fluctcur_n_2(3,3,2) 
      double precision kinetic_energy_2(2) 
      integer kmax !total number of the resonance condition roots p_perp_k(p_perp)
c-----local
      double complex sn_k(3,3),g_nk(3,3),g_fluctcur_nk(3,3)
      double precision gamma, p_par_rl(2),coeff,eps,ps,ps_max,pi,
     +dgdp_par,u_nk,u_flcur_nk,coeff_fc,
     &mass_e,clight,k_1_kev,m_cs_d_2

      integer k,i,j

      clight =2.99792458d10     !speed of light     (cm/sec)
      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)      
      m_cs_d_2=0.5d0*mass_e*clight**2/k_1_kev  !    [KeV] 

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

c      write(*,*)'in g_n_emis_spectrum: kmax',kmax

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

cSm061222
         g_fluctcur_nk(2,1)=-g_fluctcur_nk(1,2)
         g_fluctcur_nk(3,1)= g_fluctcur_nk(1,3)
         g_fluctcur_nk(3,2)=-g_fluctcur_nk(2,3)



c         kinetic_energy_2(k)= m_cs_d_2*ps
cSm070103
         kinetic_energy_2(k)=(gamma-1.d0)*2.d0*m_cs_d_2 !(gamma-1)*m_e*clight**2
           
         do i=1,3
            do j=1,3
               g_n(i,j)=g_n(i,j)+g_nk(i,j)
               g_fluctcur_n(i,j)=g_fluctcur_n(i,j)+g_fluctcur_nk(i,j)
               g_fluctcur_n_2(i,j,k)=g_fluctcur_nk(i,j)
            enddo
         enddo
          
 10      continue


      enddo !kmax


 20   continue

c      do i=1,3
c        do j=1,3
c           write(*,*)'i,j,g_fluctcur_n(i,j)',i,j,g_fluctcur_n(i,j)
c           do k=1,kmax
c             write(*,*)'k,g_fluctcur_n_2(i,j,k)',
c     &                  k,g_fluctcur_n_2(i,j,k)
c           enddo
c        enddo
c      enddo

      return
      end


      subroutine kinetic_energy_grid_spectrum(jx_kin_a,jx_kin,
     &max_kin_energy_kev,kinetic_energy_kev)

c-----fit the kinetic energy grid for the emission spectrums
      implicit none
c-----input
      integer jx_kin, ! the number of grid poits
     &jx_kin_a        ! max number of grid points 
      double precision max_kin_energy_kev ! max grid kinetic energy [KeV] 
c-----output
      double precision kinetic_energy_kev(jx_kin_a) !kinetic energy grid [KeV]
c---- locals
      integer j
      double precision step

      if( jx_kin.gt.jx_kin_a) then
        write(*,*)'in kinetic_energy_grid_spectrum'
        write(*,*)'jx_kin.gt.jx_kin_a'
        write(*,*)'jx_kin,jx_kin_a',jx_kin,jx_kin_a
        write(*,*)'Please change jx_kin in genray.dat'
        write(*,*)'or jx_kin_a in parem.i and recompile the code' 
        stop 'kx_kin'
      endif

      step= max_kin_energy_kev/(jx_kin-1)

      do j=1,jx_kin
         kinetic_energy_kev(j)=step*(j-1)
      enddo

      return
      end

      subroutine search_grid_point(x,n,a,na)
c-----search the number na
c     x(i) input grid i=1,...,n      
c
c     x(1)--xc(1)--x(2)--xc(2)--x(3)--xc(i-1)--x(i)--xc(i)--x(i+1)--xc(n-1)--x(n)
c
c     if a< x(1) na=1
c     if a> x(n) na=n
c
c     middle points: xc(i)=0,5(x(i)+x(i+1)) i=1,...,n-1
c
c     if      a < xc(1) na=1
c     if      a > xc(n-1) na=n
c     if      xc(i-1) =< a < xc(i) na=i      
c     
      implicit none
c-----input
      double precision x(*), !array x(i=1,n)
     &a                      !input argument
      integer n              !number of array points
c-----output
      integer na
c-----local
      integer i
      double precision xc_i

      if (a.lt.0.5d0*(x(1)+x(2))) then
        na =1 
        goto 10
      endif

      if (a.gt.0.5d0*(x(n-1)+x(n))) then
        na =n 
        goto 10
      endif

      do i=2,n !
        xc_i=0.5d0*(x(i)+x(i+1))
        if (a.lt.xc_i) then
           na=i
           goto 10  
        endif
      enddo

10    continue
      return
      end  

      subroutine search_grid_point_smooth(x,n,a,na,b,bp)
c-----search the number na
c     x(i) input grid i=1,...,n      
c
c     x(1)----x(2)----x(3)----x(i)----x(i+1)----x(n)
c
c     if a< x(1) na=1
c     if a> x(n) na=n
c
c     if      a < x(1) na=1, b=1, bp=0
c     if      a > x(n) na=n, b=0, bp=1
c     if      x(i) =< a < x(i+1) na=i, b=(a-x(na))/(x(na+1)-x(na)), bp=1-b           
c     
      implicit none
c-----input
      double precision x(*), !array x(i=1,n)
     &a                      !input argument
      integer n              !number of array points
c-----output
      integer na
      double precision b,bp
c-----local
      integer i

      if (a.lt.x(1)) then
        na =1
        b=1.d0 
        goto 10
      endif

      if (a.gt.x(n)) then
        na =n 
        b=0.d0
        goto 10
      endif

      do i=2,n 
        if (a.lt.x(i)) then
           na=i
           b=(a-x(na))/(x(na+1)-x(na))     
           goto 10  
        endif
      enddo

10    continue
      bp=1-b     

      return
      end  

      subroutine g_n_emis_theta_spectrum(p_perp,p_par,y,nll,np,theta,n,
     &i_fkin,
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
