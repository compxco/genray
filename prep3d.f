


c        ********************** prep3d***********************
c        *                      -----                       *
c        *  prep3d -subroutine to prepare the parameters    *
c        *  for	 output files for FP code  (e.g., CQL3D)    *
c        ****************************************************
c         input parameters: iray -number of the ray from antenna  
c                                 is in common/cone/ 
c output parameter: iraystop (if power in the ray channel 
c   delpwr(is).lt. delpwrmn*delpwr(1), then iraystop=1). 
c Also, stop ray if fluxn(is).lt.0.0, iraystop=1, and set fluxn 
c  to previous value fluxn(is-1).   [RWH:030428]
c-----------------------------------------------------------------------
      subroutine prep3d(t,u,deru,iraystop)   !ends at line 2384
      implicit real*8 (a-h,o-z)
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
            
      include 'one.i'
      include 'ions.i'
      include 'write.i'
      include 'cefield.i'
      include 'cone.i'
      include 'grill.i'
      include 'eps.i'
      include 'fourb.i'
      include 'oxb.i'
      include 'output.i'
      include 'three.i'
cSm070128 for test
      include 'six.i'
c-------------------
      include 'lsc_approach.i'
c------prmt(1;9)
      include 'rkutta.i'
      dimension u(*),deru(*),vgr(3),bf(3)
c      dimension u(6),deru(6),vgr(3),bf(3)

      
      dimension cnprim_s(nbulka)  !Added for indiv ion contrib[BH041009]
      dimension ckvipl_s(nbulka)  !Added for indiv ion contrib[BH041009]
      
      dimension tempiar(nbulka)
C------ 2 is maximum number of roots to look for in the muller root finding
      double complex cflown,cnx,fr_func_noeps,dfrfunc,roots(2)
      integer info(2),ier,nsig
      double complex dhot,dhot_rlt,dhot_sum,dispfun

      double precision x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     .tpop_ar(nbulka),vflow_ar(nbulka),image_d
      
      double complex hotnp
      double complex dhot_sumt,dhot_sum_e,dhot_sum_i
cfor_test
      double complex k_aherm(3,3)
      double complex eps1(3,3),cez1,eps2(3,3),cez2
      double complex dhot_sum_c,cd 
      external dhot_sum_c ! calculates complex hot non-relativistic 
                          ! dispersion  function
CENM For the iabsorp=12 damping option to work (find full complex solution of
C    D(nperp)=0), the muller root-finding subroutine is called, and either
C    Ram's relativistic function is used (for id=14) or Nelson-Melby's (for id=11)
      complex*16 cnper_c,cnpar_c !to test id=12 lambda
      

      external fr_func_noeps  ! for dispersion function root finder (id=14)     
      external dispfun  ! for dispersion function root finder (id=11)

cend_for_test

      double complex integral(3,3),fluctcur_n(3,3)
cfor emission test

cfor test ono tensor
      double complex K_dx_ar(3,3,nbulka),dK_dy_ar(3,3,nbulka),
     &dK_dt_ar(3,3,nbulka),dK_dnper(3,3),dK_dnpar(3,3)
      double precision mass_ar(nbulka)
cendtest
c     to calculate data for ploting the resonance ellipse boundary
      double precision p_par_rl(2)

!      external dhot

cSm030226      
c      integer nbulkaa
c      parameter (nbulkaa=5)
      double complex K(3,3),dK(3,3,7),dd(5),d,ddn,aK(3,3)
c      double complex dd5(5,nbulkaa),ddnp_h,ddnll_h,ddnp
      double complex dd5(5,nbulka),ddnp_h,ddnll_h,ddnp
 
c-----------------------------------------------------------------------
c     for test of relativistic absorption
      double complex eps_test(3,3)
c-----------------------------------------------------------------------
c     for test grpde2
cSAP110211
c      complex exde,eyde,ezde
c      real rnpar,rnper,romega,romegpe,rvgrpdc,redenfac
      complex*16 exde,eyde,ezde
      real*8 rnpar,rnper,romega,romegpe,rvgrpdc,redenfac

c     for cold plasma+relativistic tensor
      double complex disp_func,dcold_rlt

c-----for tokman flux
      real*8 flux_tokman_ar(3)

c-----for cpu_time
      real*4 time_prep3d_1,time_prep3d_2,time_prep3d
     &time_prep3d_emis_1,time_prep3d_emis_2,time_prep3d_emis

c-----to check eigenvalues
      complex*16   eigenvalue(3)

c-----to check hermitian part of K
      complex*16 K_herm(3,3)

      double precision optical_depth
c-----real*8 p_flux !for flux calculations


c-----for reflection lost 
      integer irefl_old
      real*8  tot_pow_absorb_at_refl
 
c-----for test N perpendicular coinside at the given ray pint
c     with the dispersion relation solution
      real*8
     &cnper_p,cnper_m      

cSAP120605
c-----for Nicola Bertelli scattering. It was added output argument
c     dwepsdw(3,3) to the subroutine flown (z,r,phi,cnz,cnr,cm,cflown, dwepsdw)
      complex*16 dwepsdw(3,3) 

c-----external for lsc LH absorption
cBH110715      real*8   d_f_e_maxw_d_v_lsc,d_f_e_nonmaxw_d_v_lsc
cBH110715      external d_f_e_maxw_d_v_lsc,d_f_e_nonmaxw_d_v_lsc
      real*8   d_f_e_nonmaxw_d_v_lsc
      external d_f_e_nonmaxw_d_v_lsc
c----------------------------------
      data irefl_old /0/
      data tot_pow_absorb_at_refl/0.d0/
c-------------------------------------------
      data optical_depth/0.d0/

      data time_prep3d /0./
      data time_prep3d_em /0./

      save cnprim_old
      save time_prep3d,time_prep3d_em
      save  optical_depth
      save irefl_old
      save tot_pow_absorb_at_refl

cSAP081111
      save ckvipold

c      write(*,*)'prep3d t',t
c      write(*,*)'prep3d (u(i),i=1,6)', (u(i),i=1,6)
c      write(*,*)'prep3d begin xma,yma',xma,yma
      call cpu_time(time_prep3d_1)

      nrayelt=nrayelt+1
      if (nrayelt.gt.nrelta) then
         write(*,*)'in prep3d nrayelt.gt.nrelta nrayelt,nrelta',
     .   nrayelt,nrelta
         write(*,*)'it should be nrayelt=<nrelta'
         write(*,*)'increase nrelta in param.i'
         stop
      endif  
      is=nrayelt

cSAP100507
      wt(is)=t 

      iraystop=0
      pi=4.d0*datan(1.d0)

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
cSAP080906
c----------------------------------------------------
c     poloidal angle 0< theta_pol [degree] <360
c-----------------------------------------------------
c      write(*,*)'prerp3d r,z,xma,yma', r,z,xma,yma
      call theta_xy((r-xma),(z-yma),wtheta_pol(is))
      wtheta_pol(is)=(wtheta_pol(is))*180d0/pi
      
      if (is.gt.1) then
 
c         write(*,*)'wtheta_pol(is-1),wtheta_pol(is)',
c     &              wtheta_pol(is-1),wtheta_pol(is)

         if ((wtheta_pol(is-1).le.90d0).and.
     &       (wtheta_pol(is).ge.180.d0)) then
            wtheta_pol(is)=wtheta_pol(is)-2*180.d0
         endif
      endif
c      write(*,*)'prep3d is,wtheta_pol(is)',is,wtheta_pol(is)
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

     
ctest
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
c      write(*,*)'prep3d u',u(1),u(2),u(3),u(4),u(5),u(6)
      write(*,*)'prep3d cnpar,cnper,dsqrt(cnpar**2+cnper**2) ',
     &                  cnpar,cnper,dsqrt(cnpar**2+cnper**2)
c       write(*,*)'prep3d iabsorp',iabsorp
      endif ! outprint
cendtest

cSAP100210 these lines were put here  for LCS absorption calcullations
      if(istart.eq.1) then
c        electron cyclotron waves
         wdnpar(is)=dabs(0.05d0*wnpar(1))
c        write(*,*)'in prep3d old wdnpar(is)',wdnpar(is)
      else
cSmirnov961205
c        grill conditions for the waves (LH or FW)
         wdnpar(is)=wdnpar0(iray)
c        write(*,*)'in prep3d new wdnpar(is)',wdnpar(is)
      endif

c-----waves absorption and electric field calculations----
c     Absorption if clauses end at line 1682.
      if ((iabsorp.eq.3).or.(iabsorp.eq.2)) then  !ends at line 590
c	absorption for LH and FW
c------------------------------------------------------------
c       electric field using the cold plasma dielectric tensor
        call tensrcld(u(1),u(2),u(3))
        cnx=dcmplx(cnper,0.d0)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c       electric field parallel to wave
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
        
c-------put cold plasma tensor to w_ceps array
        do i=1,3
           do j=1,3
              w_ceps(i,j,is)=reps(i,j) !from  eps.i
           enddo
        enddo
           
c-------------------------------------------------------------
        temp_e=tempe(z,r,phi,1)
        do i=2,nbulk
          tempiar(i)=tempe(z,r,phi,i)  
	enddo
       
        dens_e=dense(z,r,phi,1)
        z_eff=zeff(z,r,phi)

        if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
        write(*,*)'rho,temp_e,dens_e,z_eff',rho,temp_e,dens_e,z_eff
        do i=2,nbulk
        write(*,*)'i,tempiar(i)',i,tempiar(i)
        enddo
        endif ! outprint

        if (iabsorp.eq.3) then   !ends at line 520
c----------FW absorption
           cnprim_cl=0.d0
c----------absorpfd uses complex function modified bessel zfunc(argz,zf,zfp) 
c           call absorpf1(z,r,phi,cnpar,cnper,temp_e,dens_e,tempiar,
c----------absorpfd uses double complex function modified bessel 
c          czeta(argz,zf,zfp,ierror) and calculates the dielectric tensor
c          reps )(in common eps.i)
c          for the electron plasma with the hot correction using
c          Chiu et al, Nucl.Fus Vol. 29, No.12(1989) p.2175
c          formula (2),(3),(4),(5),(6) and (7)

cSAP091016  using  beslci1 for calulations I_n(x),
c           sometimes beslci1 could not calculate I_n(x),
c           call absorpfd(z,r,phi,cnpar,cnper,temp_e,dens_e,tempiar,
c     1                   nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
c           write(*,*)'after absorpfd cnprim_s',cnprim_s
c           write(*,*)'aft absorpfd:cnprim_e,cnprim_i',cnprim_e,cnprim_i

cSAP091016 using ibess0, ibessn for calulations I_n(x)*exp(-x)
c          It was set inside this subroutine the maximal number of harmonius
c          nb=200 
c          (nb=maximal number order of modified Bessel functions 
c           I_n) calcul

           call absorpfd_091016(z,r,phi,cnpar,cnper,temp_e,dens_e,
     &                   tempiar,
     1                   nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
           if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
           write(*,*)'after absorpfd_091016 cnprim_s',cnprim_s
           write(*,*)'aft absorpfd:cnprim_e,cnprim_i',cnprim_e,cnprim_i
           endif ! outprint
cBH050222:  Following two subroutine calls are testing of Pinsker-
cBH050222:  Porkolab expressions and generalizations, by Smirnov
cBh050222:  while at CompX, Jan-Feb, 2005.
c           call absorpfw_pinsker(z,r,phi,cnpar,cnper,temp_e,dens_e,
c     &     tempiar,nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
c           write(*,*)'after absorpfw_pisker cnprim_s',cnprim_s

c
c           call absorpfw_pinsker_1(z,r,phi,cnpar,cnper,temp_e,dens_e,
c     &     tempiar,frqncy,nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
c           write(*,*)'after absorpfw_pisker_1 cnprim_s',cnprim_s

c           call absorpfw_pinsker_2(z,r,phi,cnpar,cnper,temp_e,dens_e,
c     &     tempiar,nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
c           write(*,*)'after absorpfw_pisker_2 cnprim_e,cnprim_i',
c     &     cnprim_e,cnprim_i
c           write(*,*)'after absorpfw_pisker_2 cnprim_s',cnprim_s

c    c       call absorp_hot(z,r,phi,cnpar,cnper,nbulk,
c    c &     cnprim_e,cnprim_i,cnprim_s)
c    c       write(*,*)'after absorp_hot cnprim_s',cnprim_s

c          call absorpfw(z,r,phi,cnpar,cnper,temp_e,dens_e,tempiar,
cBH041009     1                   nbulk,bmod,cnprim_e,cnprim_i)
cBH041009: Adding cnprim_s()=individual contr to ion absorp.
     1                   
cBH040914: Adding cnprim=...  else, cnprim is from a previous call.
c           cnprim_ip=0.d0
c           do kk=2,nbulk
c              cnprim_ip=cnprim_ip+cnprim_s(kk)
c           enddo
c           write(*,*)'prep3d, after absorpfd: cnprim_i,cnprim_ip=',
c     1               cnprim_i,cnprim_ip


cSAP080303
c           if (ion_absorption.eq.'enabled') cnprim=cnprim_e+cnprim_i
c           write(*,*)'in prep3d  after absorpfd'

cSAP081111
           if (ion_absorption.eq.'enabled') then 
              cnprim=cnprim_e+cnprim_i
           else
              cnprim=cnprim_e !only electron absorption
           endif 

           if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
           write(*,*)'aft absorpfd:cnprim_e,cnprim_i,cnprim_cl,cnprim',
     1               cnprim_e,cnprim_i,cnprim_cl,cnprim
           endif ! outprint
c           write(*,*)'cnprim_s(1:nbulka)',cnprim_s
            if (cnprim.lt.0.d0) then
             if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
             write(*,*)'WARNING: cnprim <0 cnprim=',cnprim
             write(*,*)'The code will use abs(cnprim)'
             endif ! outprint
             cnprim=dabs(cnprim) 
             cnprim_e=dabs(cnprim_e)
             cnprim_i=dabs(cnprim_i)
             do i=1,nbulk
                cnprim_s(i)=dabs(cnprim_s(i))
             enddo
           endif
c----------electric field calculations using the dielectric tensor
c          from absorpfd (the electron plasma with the thermal correction)
c
c          electric field using the cold plasma dielectric tensor
cSm030512
c           call tensrcld(u(1),u(2),u(3))
c
           cnx=dcmplx(cnper,0.d0)
cSAP080916  
           cnx=dcmplx(cnper,(cnprim_e+cnprim_i))
c           cnx=dcmplx(cnper,(cnprim_e+cnprim_i))


           do i=1,3
              do j=1,3
                 w_ceps(i,j,is)=reps(i,j) !cold plasma with thermal correction
              enddo
           enddo
         
           call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c          electric field parallel to wave
           enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)

ctest of Ono dielectric tensor
cSm030512
           goto 11

           do i=1,nbulk
             x_ar(i)=x(z,r,phi,i)
	     y_ar(i)=y(z,r,phi,i)
	     te=tempe(z,r,phi,i) ! kev
	     t_av_ar(i)=te*1000.d0      ! ev 
             mass_ar(i)=dmas(i)

c             write(*,*)'prep3d i',i
c             write(*,*)'x_ar(i)',x_ar(i)             
c            write(*,*)'y_ar(i)',y_ar(i)
c            write(*,*)'t_av_ar(i)[eV]',t_av_ar(i) 

           enddo
           call ono_tens(mass_ar,x_ar,y_ar,t_av_ar,nbulk,
     &     cnpar,cnpern,1,0,
     &     eps1,d,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)
           do i=1,3
              do j=1,3
                 w_ceps(i,j,is)=eps1(i,j) !Ono tensor
              enddo
           enddo
           
           if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
        write(*,*)'w_ceps1',w_ceps(1,1,is),w_ceps(1,2,is),w_ceps(1,3,is)
        write(*,*)'w_ceps2',w_ceps(2,1,is),w_ceps(2,2,is),w_ceps(2,3,is)
        write(*,*)'w_ceps3',w_ceps(3,1,is),w_ceps(3,2,is),w_ceps(3,3,is)
           endif ! outprint
           
           k4 = 4.19d7                   !constant for elec therm vel
           c = 3.0d10 
           vt=k4*dsqrt(2.0d0*t_av_ar(1)/dmas(1)) !vt= sqrt(2.0*kT[eV]/mass)
                                                 !    cm/sec 
           if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
           write(*,*)'d=',d
           write(*,*)'k4,t_av_ar(1),dmas(1),cnpar,c,vt',
     &     k4,t_av_ar(1),dmas(1),cnpar,c,vt
           write(*,*)'y0',c/(cnpar*vt)
           endif ! outprint
           
 11        continue
cendtest_ono

cSAP081216----------------------------------------------
c           if (is.eq.81) then
           if (is.eq.1000000) then
c------------ test N_perp and polarization in one given ray point WF case
              if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
              write(*,*)'is=81'
              write(*,*)'cnpar,cnper',cnpar,cnper
              write(*,*)'cex,cey,cez',cex,cey,cez
              endif ! outprint
              call npernpar(z,r,phi,cnpar,cnper2p,cnper2m)
              if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
              write(*,*)'cnper2p,cnper2m',cnper2p,cnper2m
              endif ! outprint
              if (cnper2p.ge.0.d0) then
                 cnper_p=dsqrt(cnper2p)
                 if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
                 write(*,*)'cnper_p,cnper',cnper_p,cnper
                 endif ! outprint
              endif
cSAP090122
c              if (cnperp2m.ge.0.d0) then
              if (cnper2m.ge.0.d0) then
                 cnper_m=dsqrt(cnper2m)
                 if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
                 write(*,*)'cnper_m,cnper',cnper_m,cnper
                 endif ! outprint
              endif


              psi_loc=psi_rho(rho)
              r_m=r
              z_m=z  
              rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)
              
              if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc=Inf
              cos_theta_pol=1.d0
              sin_theta_pol=0.d0
              else
              cos_theta_pol=(r_m-xma)/rho_loc
              sin_theta_pol=(z_m-yma)/rho_loc
              endif
              
              if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
              if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
              if (sin_theta_pol.ge.0.d0) then
                 theta_pol=+dacos(cos_theta_pol)
              else  
                 theta_pol=-dacos(cos_theta_pol)
              endif

              if (theta_pol.lt.0.d0) then
                 theta_pol=theta_pol+2.d0*pi
              endif

              if (theta_pol.gt.2.d0*pi) then
                 n_theta_pol=theta_pol/(2.d0*pi)
                 theta_pol=theta_pol-2.d0*pi*n_theta_pol
              endif

c-------------calculate the hot plasma full dielectric tensor and the 
c             electric field
c              do i=1,nbulk
                 x_ar(i)=x(z,r,phi,i)
c	         y_ar(i)=y(z,r,phi,i)
cSm000324
c                 if(i.eq.1) y_ar(1)=-y_ar(1)
c	         te=tempe(z,r,phi,i) ! kev
c	         t_av_ar(i)=te*1000.d0      ! ev 
c                 tpop_ar(i)=tpoprho(rho,i)
c                 vflow_ar(i)=vflowrho(rho,i)
c              enddo
c              d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c     &                   vflow_ar,cnparp,cnperp,2,reps)
c
c              call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c-------------end calculation of the hot plasma full dielectric tensor 
c             and the electric field

cSAP110323
c              call CD_adj_efficiency_test(cnpar,cnper,
c     &        psi_loc,theta_pol,cey,cez,
c     &        eff(is))
              call CD_adj_efficiency_test(cnpar,cnper,
     &        psi_loc,theta_pol,cex,cey,cez,
     &        eff(is))
c              stop 'prep3d.f is=81'  
            endif !is.eq.81)
          
c-----------end test

        endif ! iabsorp=3

	if (iabsorp.eq.2) then
c----------LH absorption
cyup           write(*,*)'prep3d.f i_lsc_approach',i_lsc_approach

           if (i_lsc_approach.eq.1) then
cyup             write(*,*)'prep3d before LH_absorption_LSC'
c     &       ,'z,r,phi,cnpar,cnper,cnz,cnr,cm',
c     &       z,r,phi,cnpar,cnper,cnz,cnr,cm
              
              if (rho.lt.1.d0) then                 
                if (is.eq.1)then    
                   call LH_absorption_LSC(z,r,phi,cnpar,
     &             cnper,cnz,cnr,cm,reps, d_f_e_nonmaxw_d_v_lsc,
     &             wdnpar(is),cnprim_e)
                else
                   call LH_absorption_LSC(wz(is-1)/r00,wr(is-1)/r00,
     &             wphi(is-1),wnpar(is-1),
     &             wnper(is-1),wn_z(is-1),wn_r(is-1),wmtor(is-1),
     &             reps, d_f_e_nonmaxw_d_v_lsc,
     &             wdnpar(is-1),cnprim_e)
                endif
cyup                write(*,*)'LH LSC new nonmaxw cnprim_e', cnprim_e  
  
c                call LH_absorption_LSC_maxwellian(z,r,phi,cnpar,
c     &                  cnper,cnz,cnr,cm,reps,cnprim_e)

c                write(*,*)'LH maxw cnprim_e', cnprim_e 
                 cnprim_i=0.d0
                 cnprim_cl=0.d0
              else
                if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
                write(*,*)'WARNING: in prep3d rho>1 rho=',rho 
                write(*,*)'in this case LH_absorption_LSC can not work'
                write(*,*)'and the subroutine absorplh will be used'
                write(*,*)'for absorption calculations'
                endif ! outprint
                call absorplh(u,cnpar,cnper,temp_e,dens_e,tempiar
     1                  ,bz,br,bphi,nbulk,bmod,frqncy,z_eff,
     1                   cnprim_e,cnprim_i,cnprim_cl)
                cnprim_cl=coll_mult*cnprim_cl  !BH191207
              endif
           else ! i.e., i_lsc_approach.ne.1
              call absorplh(u,cnpar,cnper,temp_e,dens_e,tempiar
     1                  ,bz,br,bphi,nbulk,bmod,frqncy,z_eff,
     1                   cnprim_e,cnprim_i,cnprim_cl)
              cnprim_cl=coll_mult*cnprim_cl  !BH191207
cyup              write(*,*)'LH iabsorp_2 cnprim_e', cnprim_e
           endif  ! On i_lsc_approach.eq.1

c           call absorplh(u,cnpar,cnper,temp_e,dens_e,tempiar
c     1                  ,bz,br,bphi,nbulk,bmod,frqncy,z_eff,
c     1                   cnprim_e,cnprim_i,cnprim_cl)
	endif !iabsorp=2

cSAP080303
c        if (ion_absorption.eq.'enabled') cnprim=cnprim_e+cnprim_i   
cSAP081111
         if (ion_absorption.eq.'enabled') then
            cnprim=cnprim_e+cnprim_i
         else
            cnprim=cnprim_e
         endif
         
             cnprim=cnprim+cnprim_cl

c        cnprim=(cnprim_e+cnprim_i+cnprim_cl)
cyup        write(*,*)'cnprim_e,cnprim_i,cnprim_cl,cnprim',
cyup     1             cnprim_e,cnprim_i,cnprim_cl,cnprim
cSAP090403
c	cnprim_i=cnprim_i+cnprim_cl


      endif !iabsorp=2 or =3  Begins at line 255

      if (iabsorp.eq.1) then  !Ends at line 645
   
c-------EC wave case.The complex electric field calculations using
c       hermitian or full mazzucato tensor (with antihermitian part)
c       reps(3,3), and N_per=N_per(N_par) from mazzucato solver
        ihermloc=iherm
c       ihermloc=1
        ihermloc=2
c-------EC absorption (from Mazzucato code)
        ioptmaz=1 ! estimation of cnper1 from coldm
        cnper1=cnper
c        write(*,*)'prep3d  bef cnperm ioptmaz=1 cnper1,cnper',
c     *  cnper1,cnper
c        write(*,*)'perp3d cnper will calculate cnper1,cnprim' 
c        write(*,*)'using the estimation of nperp from cold plasma'
        ihermloc=2
        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)
cyup        write(*,*)'ioptmaz=1 ihermloc,new cnper1,old cnper,new cnprim'
cyup     *, ihermloc,cnper1,cnper,cnprim

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
c        ihermloc=1
        ihermloc=2 !full hermition + antihermition Mazzucato tens.

c        call hamiltmuz(cnpar,ihermloc,z,r,phi,cnper1,hamiltmz)
        
        call hamiltmuz_ful(cnpar,ihermloc,z,r,phi,cnper1,
     .  cnprim,hamiltmz)

c-------put Mazzucato tensor from complex N perpendicular to w_ceps array
        do i=1,3
           do j=1,3
              w_ceps(i,j,is)=reps(i,j) !from  eps.i
           enddo
        enddo

        cnx=dcmplx(cnper1,cnprim)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
c       write(*,*)'in prep3d efield mazz ex,ey,ez,enp,cnper1,cnpar',
c     1	ex,ey,ez,enp,cnper1,cnpar
      endif  !End of iabsorp.eq.1, line 592

      if(iabsorp.eq.4) then  !Ends at line 965
c       EC wave case.The complex electric field calculations
c       using Forest tensor (with antihermition part)
c       reps(3,3), and N_per=N_per(N_par) from forest solver
c-------EC absorption (from Forest code)

        cnparp=cnpar
        cnperp=cnper
     
        if (((id.eq.1).or.(id.eq.2)).or.(id.eq.3)) then
cc----------calculate N_perp(N_par) from hot plasma disp.relation   
cc          calculates two roots from the cold plasma as the initial
cc          approximation 
c           call npernpar(z,r,phi,cnpar,cnper2p,cnper2m)

c           write(*,*)'prep3d cnpar,cnper,cnper2p,cnper2m',
c     +     cnpar,cnper,cnper2p,cnper2m
cc----------set the data to common npercom.i for hotnp function
            call set_nperpcom(cnpar,nbulk,z,r,phi,dmas)
cc          calculates Nper(Npar) from hot plasma with the initial 
cc          iteration cnper from cold plasma

c           hotnperp=hotnp(nbulk,ibw,cnperp,cnper2p,cnper2m,K,iraystop)
c           write(*,*)'prep3d cnper,hotnperp',cnper,hotnperp
c           cnperp=hotnperp
        endif      
 
        do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
cSm000324
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo

        if(i_im_nperp.eq.1) then
c---------calculate ImN_perp using the formula
c         ImN_perp=abs(ImD_full/dD_hermitian/dReN_perp))
          d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .    vflow_ar,cnparp,cnperp,2,reps)

          dham=dreal(d)

	  call Ddhot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
     .    ,cnparp,cnperp,reps,dd5,ddnp_h,ddnll_h,ddnp)
c          write(*,*)'prep3d d,ddnp',d,ddnp
          cnprim = dabs(DIMAG(D) / DREAL(ddnp))
cyup          write(*,*)'cnprim=(ImD/dD/dn_perp)= ',cnprim

          call hot_nperp_muller(nbulk,dmas,x_ar,y_ar,t_av_ar,
     &    tpop_ar,vflow_ar,cnparp,cnperp,cnprim)
cyup          write(*,*)'muller cnprim',cnprim 
          
          do i=1,3
             do j=1,3
                w_ceps(i,j,is)=reps(i,j) !hot plasma tensor
             enddo
          enddo

c---------dD/domega
c         wf=frqncy
c         call hotdervs(nbulk,u,wf,dddcnz,dddcnr,dddcm,
c     .   dddz,dddr,dddph,dddw)
c         vgr_perp_w=dreal(ddnp)/dddw
c         wi=-dimag(d)/dddw
c         write(*,*)'dimagD,dddw,wi',dimag(D),dddw,wi
c-------------------------------
          cnper_new=cnperp !it will be used in the electric field calculations
        endif  ! i_im_nperp.eq.1  

        if(i_im_nperp.eq.2) then 
c---------find (Im_N_perp,ReN_perp) the root of the complex dispersion relation
c         using the Newton method with numerical derivatives (the chord method)
          iter_max=100 !max number of the iterations
          iter_max=3 !max number of the iterations
c---------initial values of Im_N_perp=cnprim Re_N_perp=cnperp
          cnprim=0.d0   
c          cnper_new=cnperp
c          cnprim=cnprim_old

          if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
          write(*,*)'prep3d before call solv_nperp_hot cnprim=',cnprim
          write(*,*)'nbulk,dmas,x_ar,y_ar',
     &               nbulk,dmas,x_ar,y_ar
          write(*,*)'t_av_ar',t_av_ar
          write(*,*)'tpop_ar',tpop_ar      
          write(*,*)'vflow_ar',vflow_ar
          write(*,*)'cnparp,cnperp,iter_max,cnper_new,cnprim',
     &               cnparp,cnperp,iter_max,cnper_new,cnprim
          endif ! outprint

          call solv_nperp_hot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .    vflow_ar,cnparp,cnperp,iter_max,cnper_new,cnprim)
cyup          write(*,*)'Newton cnperp,cnper_new,cnprim',
cyup     .    cnperp,cnper_new,cnprim
          cnprim_old=cnprim
         
ctest
c          temp_e=tempe(z,r,phi,1)
c          do i=2,nbulk
c             tempiar(i)=tempe(z,r,phi,i)
c	  enddo
c
c          dens_e=dense(z,r,phi,1)
c          z_eff=zeff(z,r,phi)
c
c           call absorpfd(z,r,phi,cnpar,cnper_new,temp_e,dens_e,tempiar,
c     1                   nbulk,bmod,cnprim_e,cnprim_i)
c           cnprim=cnprim_e+cnprim_i
c
cendtest
          cnprim=dabs(cnprim)

ctest
c        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c     .  vflow_ar,cnparp,cnper_new,2,reps)
c        cnprim = dabs(DIMAG(D) / DREAL(ddnp))
c        write(*,*)'nperp from Newton, cnprim=(ImD/dD/dn_perp)= ',cnprim
cendtes
        endif !i_im_nperp.eq.2

c       gradient method
c       write(*,*)'before gradient method'
c       call solv_nperp_hot_grad(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c     . vflow_ar,cnparp,cnprim,cnperp,cnper_new,cnprim_new)
c       write(*,*)'gradient method cnperp,cnper_new,cnprim_new',
c     . cnperp,cnper_new,cnprim_new
c       cnprim_i=0.d0



cSm020312 
c       call dhot_sum_ei(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c     . vflow_ar,cnparp,cnperp,reps,dhot_sumt,dhot_sum_e,dhot_sum_i)
c       write(*,*)'prep3d dhot_sumt',dhot_sumt
c       write(*,*)'prep3d dhot_sum_e',dhot_sum_e
c       write(*,*)'prep3d dhot_sum_i',dhot_sum_i
c       cnprim_e=dabs(dimag(dhot_sum_e)/dreal(ddnp))
c       cnprim_i=dabs(dimag(dhot_sum_i)/dreal(ddnp))
       

cSm020312 numerical derivatives
c       step=1.d-7 !for n_perp
c
c       cnper_pl=cnperp+step
c       d_pl=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c    .  vflow_ar,cnparp,cnper_pl,2,reps)
c
c       cnper_min=cnperp-step
c       d_min=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c    .  vflow_ar,cnparp,cnper_min,2,reps)
c        
c       ddnp=(d_pl-d_min)/(2.d0*step)
c       write(*,*)'prep3d numerical deriv. ddnp',ddnp
c              
c       cnprim1 = dabs(DIMAG(D) / DREAL(ddnp))
c       write(*,*)'prep3d Ddhot cnprim1=',cnprim1
cend_numerical deriv.

c       write(*,*)'prep3d te_ev'1,te_ev
c	write(*,*)'prep3d forest cnper,cnprim',cnper,cnprim
        cnprim_cl=0.d0
        cnprim_e=dabs(cnprim)
        cnprim_i=0.d0

c-------the  calculation of the contours 
c       of the dispersion function D at the plane (ImN_perp,ReN_perp)
c       go to 30

c-------set the map(Re_Nperp,ImN_perp) dimensions          


c        m_r_nperp=10
c        m_i_nperp=10
      
c        dmax_i_nperp=2.5d0*dabs(cnprim)
c        dmin_i_nperp=-2.5d0*dabs(cnprim)
 
cSm021031

c        dmax_r_nperp=cnper_new*1.5d0
c        dmin_r_nperp=cnper_new*0.5d0
 
                                     
c        if (nrayelt.lt.30) goto 30

        if (i_plot_d.eq.1) then
c          plotting dispersion function contours D(ImN_perp,ReN_perp)
c          to plot.ps file using PGplot
           m_r_nperp=number_map_points_real_nperp 
           m_i_nperp=number_map_points_image_nperp 
           dmax_r_nperp=cnper_new*ratio_max_r_nperp
           dmin_r_nperp=cnper_new*ratio_min_r_nperp

           if (dabs(cnprim).gt.1.d0) then        
             dmax_i_nperp=dabs(cnprim)*ratio_min_i_nperp
             dmin_i_nperp=dabs(cnprim)*ratio_max_i_nperp
           else
             dmax_i_nperp= 1.d0
             dmin_i_nperp=0.d0
           endif
cyup           write(*,*)'prep3d.f before map_dhot_nper'
cyup           write(*,*)'ratio_min_r_nperp,ratio_max_r_nperp,cnper',
cyup     &     ratio_min_r_nperp,ratio_max_r_nperp,cnper
cyup           write(*,*)'number_map_points_real_nperp',
cyup     &     number_map_points_real_nperp
cyup           write(*,*)'number_map_points_image_nperp',
cyup     &     number_map_points_image_nperp

           call map_dhot_nper(m_r_nperp,m_i_nperp, 
     .     dmax_r_nperp,dmin_r_nperp,dmax_i_nperp,dmin_i_nperp,
     .     n_contour_plot_disp,
     .     nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .     vflow_ar,cnparp,cnper_new,cnprim,iabsorp,
     .     n_relt_harm1,n_relt_harm2,n_relt_intgr,
     .     i_resonance_curve_integration_method,epsi,
     .     i_diskf,r,z,phi )
           
 
           m_y=100
           dmax_r_nperp=200.d0
           dmin_r_nperp=30.d0
           dmin_y=0.45
           dmax_y=0.95d0
c          plotting dispersion funxtion contours D(ReN_perp,Y)
c           call map_dhot_nper_omeg(m_r_nperp,m_y, 
c     .     dmax_r_nperp,dmin_r_nperp,dmax_y,dmin_y,
c     .     nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c     .     vflow_ar,cnparp)

        endif
 30     continue
        cnprim=dabs(cnprim)

c-------electric field for Forest tensor
c       cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)

c-------Hermitian tensor
c       d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c     . vflow_ar,cnparp,cnperp,1,reps)

c-------full tensor with cnperp calculated by Runge-Kutta 
c       along the ray trajectrory

c       d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c     . vflow_ar,cnparp,cnperp,2,reps)      
c       cnx=dcmplx(cnper,cnprim)

c-------full tensor with new n_perp calculated by solver
        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
cSAP090813 cnper_new gave jumps in reps
c     .  vflow_ar,cnparp,cnper_new,2,reps)     
     .  vflow_ar,cnparp,cnper,2,reps)      
        cnx=dcmplx(cnper_new,cnprim)
        cnx=dcmplx(cnper,0.d0)
 

c for checking of the dielectric tensor
c       do i=1,3
c          do j=1,3
c              w_ceps(i,j,is)=reps(i,j)!hot plasma
c           enddo
c        enddo
cendtes
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
c        write(*,*)'prep3d hot ez',ez
         
ctest of the cold plasma polarization
c        do i=1,3
c           do j=1,3
c              eps1(i,j)=reps(i,j) !put hot plasma tensor
c           enddo
c        enddo
c        cez1=cez !hot
c        ez1=ez   !hot

c        call tensrcld(u(1),u(2),u(3))
c        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c        write(*,*)'prep3d cold ez',ez
c
c        do i=1,3
c           do j=1,3
c              eps2(i,j)=reps(i,j) !put cold plasma tensor
c           enddo
c        enddo
c        cez2=cez !cold
c        ez2=ez   !cold

c-------test of cold plasma polarization with the thermal correction 
c        temp_e=tempe(z,r,phi,1)
c        do i=2,nbulk
c          tempiar(i)=tempe(z,r,phi,i)
c	enddo

c        dens_e=dense(z,r,phi,1)       
c        call absorpfd(z,r,phi,cnparp,cnperp,temp_e,dens_e,tempiar,
c     1                   nbulk,bmod,cnprim_e,cnprim_i)
c        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c        write(*,*)'prep3d cold_correction ez',ez

c        write(*,*)'i,j,hot eps1(i,j),cold eps2,col_correct_reps'
c        do i=1,3
c           do j=1,3
c              write(*,*)i,j
c              write(*,*)eps1(i,j),eps2(i,j),reps(i,j)
c           enddo
c        enddo
c        write(*,*)'cez1,cez2,cez',cez1,cez2,cez

c        write(*,*)'ez1,ez2,ez',ez1,ez2,ez   
c endtest cold polarization
  
c       write(*,*)'in prep3d efield forest  ex,ey,ez,enp,cnper,cnpar',
c     1	ex,ey,ez,enp,cnper,cnpar

      endif ! iabsorp.eq.4, Begins line 647

c------------------------------------------------------------------
      if(iabsorp.eq.5) then  !Ends line 991
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
c	write(*,*)'prep3d shkarofsky cnper,cnprim',cnper,cnprim
        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0

c-------electric field for Shkarofsky tensor
        call efiedshk (cnpar,cnper,xe,ye,te_kev,friq,cex,cey,cez,
     *  ex,ey,ez)      
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
c       write(*,*)'in prep3d efield sharofsky ex,ey,ez,enp,cnper1,cnpar',
c     1	ex,ey,ez,enp,cnper1,cnpar

      endif ! iabsorp.eq.5    !  Begins on line 968

      if(iabsorp.eq.6) then  !Ends at line 1179
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
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo
 
c        write(*,*)'prep3d.f before dhot_sum'        
c       calculate the Stix hot plsams disp function and tensor
c       with forest.f:
        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .   vflow_ar,cnparp,cnperp,1,reps)

c        write(*,*)'prep3d.f after dhot_sum '        

c-------anti-hermitian relativistic tensor aK   
c       n_relt_harm is the total number of EC harmonics in relativistic 
c              anti-hermitian dielectric tensor calculations
c   n_relt_harm1 min number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations
c   n_relt_harm2 max number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations
c               n_relt_harm1 =< n= <n_relt_harm2

c       n_relt_intgr is the number of points for the numerical integration
c       over p_perp for the calculations of anti-hermitian
c       relativistic tensor aK.
        
        if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c----------usage of the analytical relativistic function and its derivatives 
           i_fkin=0
        else 
c----------usage of the mech relativistic function and its derivatives
           i_fkin=1
        endif
cSm060314 
c       ye=dabs(y_ar(1)) !positive ye for positive harmonics
c       call anth_rlt(x_ar(1),ye,t_av_ar(1)*1.d-3,cnparp,cnperp,

        call anth_rlt(x_ar(1),y_ar(1),t_av_ar(1)*1.d-3,cnparp,cnperp,
     +  n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +  i_resonance_curve_integration_method,epsi,
     +  i_fkin,r,z,phi,
     +  aK)
         
        cnprimp=0.d0
        d=dhot_rlt(reps,aK,cnparp,cnperp,cnprimp)
        dham=dreal(d)
cyup        write(*,*)'reps',reps
cyup        write(*,*)'aK',aK
cyup        write(*,*)'d',d
	call Ddhot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
     . ,cnparp,cnperp,reps,dd5,ddnp_h,ddnll_h,ddnp)
        
        cnprim = dabs(DIMAG(D) / DREAL(ddnp))
c        cnprim = -DIMAG(D) / DREAL(ddnp)

c	write(*,*)'prep3d forest cnper, relativistic cnprim',cnper,cnprim
cyup        write(*,*)'dimag(d)',dimag(d)
cyup        write(*,*)'dreal(ddnp),cnprim',dreal(ddnp),cnprim

        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0

        goto 50

c        if (dabs(cnprim).lt.1.d-30) cnprim=1.d-1 
 
c        the number of mesh points for the map D(ReN_perp,Im_N_perp)        
c           m_r_nperp=100
c           m_i_nperp=100
c           m_r_nperp=10
c           m_i_nperp=10

c           dmax_r_nperp=cnperp*1.5d0
c           dmin_r_nperp=cnperp*0.5d0
c           dmax_i_nperp=2.d0*dabs(cnprim)
c           dmin_i_nperp=-2.d0*dabs(cnprim)

c           write(*,*)'dmax_i_nperp,dmin_i_nperp',
c     .     dmax_i_nperp,dmin_i_nperp
 
c           if(dmax_i_nperp.gt.100.d0) then
c             dmax_i_nperp=100.d0
c             dmin_i_nperp=-100.d0
c           endif


           if (nrayelt.lt.60) goto 50       
           if(i_plot_d.eq.1) then
c-------------plotting dispersion function contours D(ImN_perp,ReN_perp)
c-------------to plot.ps file using PGplot
              m_r_nperp=number_map_points_real_nperp 
              m_i_nperp=number_map_points_image_nperp 

              dmax_r_nperp=cnperp*ratio_max_r_nperp
              dmin_r_nperp=cnperp*ratio_min_r_nperp

             if (dabs(cnprim).gt.1.d0) then        
                dmax_i_nperp=dabs(cnprim)*ratio_min_i_nperp
                dmin_i_nperp=dabs(cnprim)*ratio_max_i_nperp
             else
                dmax_i_nperp= 1.d0
                dmin_i_nperp=0.d0
             endif

              call map_dhot_nper(m_r_nperp,m_i_nperp, 
     .        dmax_r_nperp,dmin_r_nperp,dmax_i_nperp,dmin_i_nperp,
     .        n_contour_plot_disp,
     .        nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .        vflow_ar,cnparp,cnperp,cnprim,iabsorp,
     .        n_relt_harm1,n_relt_harm2,n_relt_intgr,
     .        i_resonance_curve_integration_method,epsi,
     .        i_diskf,r,z,phi)
           endif

         
     
 50        continue

        goto 40
ctest Mazzucato
c       EC wave case.The complex electric field calculations
c       using hermition or full mazzucato tensor (with antihermition part)
c       reps(3,3), and N_per=N_per(N_par) from mazzucato solver
        ihermloc=iherm
c       ihermloc=1
        ihermloc=2
c-------EC absorption (from Mazzucato code)
        ioptmaz=1 ! estimation of cnper1 from coldm
        cnper1=cnper
c        write(*,*)'prep3d  bef cnperm ioptmaz=1 cnper1,cnper',
c     *  cnper1,cnper
c        write(*,*)'perp3d cnper will calculate cnper1,cnprim' 
c        write(*,*)'using the estimation of nperp from cold plasma'
        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)

c        write(*,*)'prep3d ioptmaz=1 ihermloc,cnper1,cnper,cnprim'
c     *, ihermloc,cnper1,cnper,cnprim

        ioptmaz=2 ! estimation of cnper1 from input parameter cnper
        cnper1=cnper
        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)
cyup        write(*,*)'prep3d ioptmaz=2 cnper1,cnper,cnprim',
cyup     +  cnper1,cnper,cnprim

        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0

        call aherm(reps,k_aherm)
cyup        write(*,*)'prep3d ani_herm Muzzucato',k_aherm
cend_test_Mazzucato
 40     continue

c-------electric field for Forest tensor
        cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)

c        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
c     .  vflow_ar,cnparp,cnperp,1,reps)

        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .  vflow_ar,cnparp,cnperp,2,reps)

        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
c       write(*,*)'in prep3d efield forest  ex,ey,ez,enp,cnper,cnpar',
c     1	ex,ey,ez,enp,cnper,cnpar

      endif ! iabsorp.eq.6, Begins line 993
c-------------------------------------------------------------
      if(iabsorp.eq.7) then  !Ends line 1284
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
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo
         
ctest   relativistic Maxwell normalization
c        if (is.eq.25) call test_density_relat_maxwell
cendtest

c-------anti-hermitian relativistic tensor aK   
c       n_relat_harm is the total number of EC harmonics in relativistic tensor aK
c       n_relt_harm1 min number of EC harmonics used in anti-hermitian
c                 dielectric tensor calculations
c       n_relt_harm2 max number of EC harmonics used in anti-hermitian
c                 dielectric tensor calculations
c                 n_relt_harm1 =< n= <n_relt_harm2

c       n_relt_intgr is the number of points for the numrerical integration
c                 over p_perp for the calculations of anti-hermitian
c                 relativistic tensor aK.
        
        if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c----------use of the analytical relativistic function and its derivatives 
           i_fkin=0
        else 
c----------use of the numerical relativistic function and its derivatives
           i_fkin=1
        endif

c        write(*,*)'prep3d before anth_rlt x_ar(1),y_ar(1)',
c     &  x_ar(1),y_ar(1)
c        write(*,*)'t_av_ar(1)*1.d-3,cnparp,cnperp',
c     &  t_av_ar(1)*1.d-3,cnparp,cnperp
c        write(*,*)'n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin',
c     &  n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin
c        write(*,*)'r,z,phi,rho',r,z,phi,rho

cSm060314 
c        ye=dabs(y_ar(1)) !positive ye for positive harmonics
        
         !write(*,*)'prep3d iabsorp==7 before anth_rlt'

         call anth_rlt(x_ar(1),y_ar(1),t_av_ar(1)*1.d-3,cnparp,cnperp,
     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +   i_resonance_curve_integration_method,epsi,
     +   i_fkin,r,z,phi,
     +  aK)
 
         !write(*,*)'prep3d iabsorp==7 after anth_rlt aK',aK

c------ complex dispersion function calculated from the sum of
c       of the cold electron plasma dielectric tensor eps_h
c       and the relativistic electron anti-hermition dielectric tensor eps_a

cyup        write(*,*)'prep3d reps',reps

        disp_func=dcold_rlt(reps,aK,cnparp,cnperp)

c        write(*,*)'reps',reps
c        write(*,*)'aK',aK
c        write(*,*)'  disp_func',disp_func

c-------calculate the derivative d(D_hermitian)/d(ReN_perp)
c       from the electron cold plasma dispersion function D
        ddnp=dDcold(reps,cnpar,cnper)
c        write(*,*)'ddnp',ddnp
        
        cnprim = dabs(DIMAG(disp_func) / DREAL(ddnp))
	
c        write(*,*)'iabsorp=7 DIMAG(disp_func)', DIMAG(disp_func)
cyup        write(*,*)'DREAL(ddnp),cnprim',DREAL(ddnp),cnprim

        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0
     
c-------electric field for cold plasma
        cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)
        
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)

      endif ! iabsorp.eq.7, Begins line 1181


      if(iabsorp.eq.8) then !Ends line 1309
c-------The absorption is calculated for Westerhof- Tokman dispersion
c       It works for id =10,12,13,15 cases.
        cnprim_in=0.d0
        call Im_nperp_Westerhof_Tokman(z,r,phi,
     &  cnpar,cnper,cnprim_in,id,cnprim)

cyup        write(*,*)'after Im_nperp_Westerhof_Tokman cnper,cnprim',
cyup     &  cnper,cnprim

        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0

c-------electric field using the dielectric tensor acording id=10,12 or 13
        cnx=dcmplx(cnper,cnprim)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c       electric field parallel to wave
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
        
      endif !iabsorp=8, begins line 1287	
c--------------------------------------------------------------
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
    
cSAP080303
c         if (ion_absorption.eq.'enabled') cnprim=cnprim_e+cnprim_i
cSAP08111
         if (ion_absorption.eq.'enabled') then
            cnprim=cnprim_e+cnprim_i
         else
            cnprim=cnprim_e
         endif        
c         cnprim=cnprim_e+cnprim_i
         cnprim_cl=0.d0
        
cyup         write(*,*)'after absorp_hot cnprim_e,cnprim_i,cnprim_s',
cyup     &   cnprim_e,cnprim_i,cnprim_s
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
        
cyup         write(*,*)'after absorp_relativist_disp_combined'
cyup         write(*,*)'cnprim_e',cnprim_e

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

c----------------------------------------------------------------
      if(iabsorp.eq.11) then
c--------------------------------------------------------------
c        The absorption is calculated for relativistic dispersion
c        D=determinant using relativistic tensor (A.Ram if id=14 or
c        Nelson-Melby if id=11) and projection method
         X_e=x(z,r,phi,1)
         Y_e=y(z,r,phi,1)
         T_e=tempe(z,r,phi,1)  

         call relativist_absorp_projection_method_det
     &   (T_e,cnpar,X_e,Y_e,cnper,reps,D,cnprim_e)

         cnprim_i=0.d0     
         cnprim=cnprim_e+cnprim_i
         cnprim_cl=0.d0
        
cyup         write(*,*)'aft absorp_relativist iabsorp=11 cnprim_e',cnprim_e	
c----------------------------------------------------------------
c        electric field polarization (cex,cey,cez)
c        calculated using the relativistic dielectric tensor
c        will be common /efield.i/
c---------------------------------------------------------------
c        relativistic A.Ram(full, non-hermitian) dielectric tensor reps(3,3)
c        will be in common /eps.i/
c--------------------------------------------------------------- 
          cnx=dcmplx(cnper,cnprim)
          call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        
         endif !iabsorp=11

C---------------------------------------------
      if (iabsorp.eq.12) then   !Ends line 1559
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
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
          write(*,*)'()()(()()(()()initial data'
          write(*,*)'z=',z,'r=',r,'phi=',phi
          write(*,*)'cnz=',cnz,'cnr=',cnr,'cm=',cm
          write(*,*)'rho=',rho
          print *,'====== is=',is
          print *,'cnprim: ',cnprim,' cnpar: ',cnpar,' cnper: ',cnper
         endif ! outprint
c      cn2=cnz**2+cnr**2+cnphi**2
      cn2=cnz**2+cnr**2+(cm/r)**2

cyup      print *,'cn2=',cn2,' cnper(calculated) ',dsqrt(cn2-cnpar**2)
cSAP091004
c        if (is.lt.1) then
         if (is.le.1) then
            print *,'****** is=',is
         else
            print *,'nper before: ',wnper(is-1)
         endif

C id.eq.11 or id.eq.12, call root-finding using Nelson-Melby dielectric function
         if (id.eq.11 .or. id.eq.12) then
cSAP091104
cyup           write(*,*)'prep3d.f before call muller'
c SAP091216 
           X_e=x(z,r,phi,1)
           Y_e=y(z,r,phi,1)
           T_e=tempe(z,r,phi,1)  
           cnx=dcmplx(cnper,0.d0)
           call Disp_Nelson_Melby(T_e,cnpar,X_e,Y_e,
     &   cnx,reps,D)

           call muller(dispfun,errabs,nsig,nknown,nguess,nnew,roots,
     +     itmax,info,ier)

cSAP091104
cyup           write(*,*)'prep3d.f after call muller'

         else
C id.eq.14 or 15, call fr_func_noeps, using Ram's dielectric function
 
c-----------print relativistic tensor for testing
            X_e=x(z,r,phi,1)
            Y_e=y(z,r,phi,1)
            T_e=tempe(z,r,phi,1)  
            cnx=dcmplx(cnper,0.d0)
            if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
            write(*,*)'Xe,Y_e ',Xe,Y_e
            write(*,*)'cnpar',cnpar
            write(*,*)'cnx',cnx
            endif ! outprint
            call Disp_Ram(T_e,cnpar,X_e,Y_e,cnx,K,d) !K is in z-y plane
                               !in Stix coordinates
cyup            write(*,*)'prep3d K',K 
            call herm(K,K_herm)
cyup            write(*,*)'prep3d K_herm',K_herm
cyup            write(*,*)'D ',D 
c           end print relativistic tensor for testing
c-------------------------------------------------------------

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
cyup         print *,'&&&&&&&&&&&& cnprim: ',cnprim,' cnper: ',cnper

         cnx=dcmplx(cnper,cnprim)

cSAP091104
         X_e=x(z,r,phi,1)
         Y_e=y(z,r,phi,1)
         T_e=tempe(z,r,phi,1)  
c-----------print relativistic tensor for testing            
c            write(*,*)'X_e,Y_e',X_e,Y_e
c            write(*,*)'T_e ', T_e
c            write(*,*)'cnpar',cnpar
c            write(*,*)'cnx',cnx
c            call Disp_Ram(T_e,cnpar,X_e,Y_e,cnx,K,d) !K is in z-y plane
                               !in Stix coordinates
c            write(*,*)'prep3d K',K 
c            call herm(K,K_herm)
c            write(*,*)'prep3d K_herm',K_herm
c            write(*,*)'D ',D 
c-----------end print relativistic tensor for testing

c-------------------------------------------------------------
cSm060719
c--------calculate dielectric tensor reps for electric field calculations
cSAP091104
c         write(*,*)'prep3d.f before Disp_combined'
c         write(*,*)'T_e,cnpar,X_e,Y_e,cnx',T_e,cnpar,X_e,Y_e,cnx

         call Disp_combined(T_e,cnpar,X_e,Y_e,cnx,reps,d)

cSAP091104
c         write(*,*)'prep3d.f Disp_combined'

cSAP091216 
         cnper_c=cmplx(cnper,cnprim)
         cnpar_c=cmplx(cnpar,0.d0)    
         call lambda(cnper_c,cnpar_c,reps,eigenvalue)
  
cyup         write(*,*)'eigenvalue',eigenvalue

c         eps=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))
c         write(*,*)'in prep3d  epshamilt1=',eps
        
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
      endif   !iabsorp=12
      if(iabsorp.eq.91) then
c--------------------------------------------------------------
c        The FW absorption is calculated using
c        formula from using formula:   	                            
c        1) for ion absorption from                      
c           R.i.Pinsker,M.Porkolab                          
c           e-mail: pinsker@fusion.gat.com                  
c           05/02/11                                         
c        2) for electron absorption from                 
c          S.C.Chiu,V.S.Chan,R.W.Harvey,M.Porkolab	 
c          Theory of fast wave current drive for tokamak   
c          plasma,					 
c          Nuclear fusion,1989,vol.29,No.12, p.2175-2186   
c       
c        It calculates the dielectric tensor reps
c        using S.C.Chiu article
c----------------------------------------------------------------
         temp_e=tempe(z,r,phi,1)
         do i=2,nbulk
            tempiar(i)=tempe(z,r,phi,i)
	 enddo  
         dens_e=dense(z,r,phi,1)
         z_eff=zeff(z,r,phi)

         call absorpfw_pinsker_1(z,r,phi,cnpar,cnper,temp_e,dens_e,
     &   tempiar,frqncy,nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
cyup         write(*,*)'after absorpfw_pisker_1 cnprim_s',cnprim_s

cSAP080303
c        if (ion_absorption.eq.'enabled') cnprim=cnprim_e+cnprim_i
cSAP08111
        if (ion_absorption.eq.'enabled') then
           cnprim=cnprim_e+cnprim_i
        else
           cnprim=cnprim_e
        endif

c         cnprim=cnprim_e+cnprim_i
         cnprim_cl=0.d0
        
cyup         write(*,*)'after absorpfw_pinsker_1 cnprim_e,cnprim_i,cnprim_s'
cyup     &   ,cnprim_e,cnprim_i,cnprim_s
c----------------------------------------------------------------
c        Calulate the electric field polarization (cex,cey,cez) 
c        using dielectric tensor reps from S.C.Chiu article 
c        will be common /efield.i/
c---------------------------------------------------------------
c        dielectric tensor reps(3,3)
c        will be in common /eps.i/
c--------------------------------------------------------------- 
c        cnx=dcmplx(cnper,cnprim)
c        electric field using the cold plasma dielectric tensor
         call tensrcld(u(1),u(2),u(3))
         cnx=dcmplx(cnper,0.d0)
         call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c        electric field parallel to wave
         enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
      endif !iabsorp=91


	
      if(iabsorp.eq.92) then
c--------------------------------------------------------------
c        The FW absorption is calculated using
c        formula from using formula:   	                            
c        for ion and electron absorption from                      
c           R.i.Pinsker,M.Porkolab                          
c           e-mail: pinsker@fusion.gat.com                  
c           05/02/11                                         
c       
c        It calculates the dielectric tensor reps
c        using the article
c        S.C.Chiu,V.S.Chan,R.W.Harvey,M.Porkolab	 
c        Theory of fast wave current drive for tokamak   
c        plasma,					 
c        Nuclear fusion,1989,vol.29,No.12, p.2175-2186   
c----------------------------------------------------------------
         temp_e=tempe(z,r,phi,1)
         do i=2,nbulk
            tempiar(i)=tempe(z,r,phi,i)
	 enddo  
         dens_e=dense(z,r,phi,1)
         z_eff=zeff(z,r,phi)

         call absorpfw_pinsker_2(z,r,phi,cnpar,cnper,temp_e,dens_e,
     &   tempiar,nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
cyup         write(*,*)'after absorpfw_pisker_1 cnprim_s',cnprim_s

cSAP080303
c        if (ion_absorption.eq.'enabled') cnprim=cnprim_e+cnprim_i 
cSAP080303
        if (ion_absorption.eq.'enabled') then
           cnprim=cnprim_e+cnprim_i
        else
           cnprim=cnprim_e
        endif

c         cnprim=cnprim_e+cnprim_i
         cnprim_cl=0.d0
        
cyup         write(*,*)'after absorpfw_pinsker_1 cnprim_e,cnprim_i,cnprim_s'
cyup     &   ,cnprim_e,cnprim_i,cnprim_s
c----------------------------------------------------------------
c        Calulate the electric field polarization (cex,cey,cez) 
c        using dielectric tensor reps from S.C.Chiu article 
c        will be common /efield.i/
c---------------------------------------------------------------
c        dielectric tensor reps(3,3)
c        will be in common /eps.i/
c--------------------------------------------------------------- 
c        cnx=dcmplx(cnper,cnprim)
c        electric field using the cold plasma dielectric tensor
         call tensrcld(u(1),u(2),u(3))
         cnx=dcmplx(cnper,0.d0)
         call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
         enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
      endif !iabsorp=92	


C======================================================================
C----------------------------- END OF IABSORP MODULES -----------------
C======================================================================


c--------------------------------------------------------------------
c        put dielectric tensor reps into w_ceps 
c--------------------------------------------------------------------
         do i=1,3
           do j=1,3
               w_ceps(i,j,is)=reps(i,j) !Put the tensor reps for writing
           enddo
         enddo

c      nrayelt=nrayelt+1

c      if (nrayelt.gt.nrelta) then
c         write(*,*)'in prep3d nrayelt.gt.nrelta nrayelt,nrelta',
c     .   nrayelt,nrelta
c         write(*,*)'it should be nrayelt=<nrelta'
c         write(*,*)'increase nrelta in param.i'
c         stop
c      endif  

    
      is=nrayelt
c        write(*,*)'in prep3d is',is
      wye(is)=y(z,r,phi,1)
      if(nbulk.gt.1) wyi(is)=y(z,r,phi,2)

      wxe(is)=x(z,r,phi,1)
      if(nbulk.gt.1) wxi(is)=x(z,r,phi,2)

cBH050310
cBH050310c-----for plotting  ion-ion resonance by draw     
cBH050310      if (nbulk.gt.2) then
cBH050310        wyi2(is)=y(z,r,phi,4)
cBH050310        wxi2(is)=x(z,r,phi,4)
cBH050310      endif
               if (nbulk.gt.3) then
                 wyi2(is)=y(z,r,phi,4)
                 wxi2(is)=x(z,r,phi,4)
               endif
c---------------------------------------------------------
c     ws (cm)
      if (is.eq.1) then
c         write(*,*)'in prep3d is=1 rho',rho
         ws(1)=0.d0  !poloidal ray distance 
         wsn(1)=0.d0 !total ray distance for emission
         rhoold=rho
         zold=0.d0
         rold=0.d0 
         phiold=0.d0
         i_ox_conversion=0
c         write(*,*)'prep3d is=1 wsn(1) ',wsn(1)
      else
c         write(*,*)'prep3d is,z,r,phi', is,z,r,phi
c         write(*,*)'prep3d zold,rold,phiold',zold,rold,phiold
c         write(*,*)'dsqrt((z-zold)**2+(r-rold)**2)*r0x*100.d0',
c     &              dsqrt((z-zold)**2+(r-rold)**2)*r0x*100.d0

         ws(is)=ws(is-1)+dsqrt((z-zold)**2+(r-rold)**2)*r0x*100.d0

c         write(*,*)'prep3d ws(is-1),ws(is)',ws(is-1),ws(is)
c         wsn(is)=ws(is-1)+dsqrt((z-zold)**2+(r-rold)**2+
c     +           (0.5d0*(r+rold)*(phi-phiold))**2)*r0x*100.d0
         wsn(is)=wsn(is-1)+dsqrt((z-zold)**2+(r-rold)**2+
     +           (0.5d0*(r+rold)*(phi-phiold))**2)*r0x*100.d0
c         write(*,*)'ws(is-1),wsn(is-1)',ws(is-1),wsn(is-1)
         delws=dsqrt((z-zold)**2+(r-rold)**2)*r0x*100.d0
c         write(*,*)'delws',delws
         delwsn=dsqrt((z-zold)**2+(r-rold)**2+
     +           (0.5d0*(r+rold)*(phi-phiold))**2)*r0x*100.d0
c         write(*,*)'delwsn',delwsn
      end if
cyup      write(*,*)'is,ws(is)',is,ws(is)
cyup      write(*,*)'irefl',irefl
c      if(is.gt.1) write(*,*)'is,wsn(is-1),wsn(is)',is,wsn(is-1),wsn(is)
cyup      write(*,*)'rho',rho

      psi_s=psi_rho(rho)
      q_s=qsafety_psi(psi_s)
c      write(*,*)'psi_s,q_s',psi_s,q_s
      wphi(is)=u(3)

c--------------------------------------------------------
c      call prepripl(t,u,is)
      
      call prepebw(t,u,is)
           
      zold=z
      rold=r
      
      phiold=phi
c--------------------------------------------------------
c     cflown - dimensionless for E_x/E,E_y/E,E_z/E
c  !!!!now flown is calculated using Mazzucato dielectric tensor
c  !!!!and electric field was calculated using Mazzucato tensor
c  !!!!it is only for EC wave .For LH and FW it is necessery
c !!!!!to create new subroutine flown ?what tensor?
c      call flown (u(1),u(2),u(3),u(4),u(5),u(6),
c     1                   cflown,dwepsdw)


c       id_old=id
c       if (iabsorp.eq.4) then
c          id=6 !hot plasma dispersion
c          call flown (u(1),u(2),u(3),u(4),u(5),u(6),cflown,dwepsdw)
c          id=id_old
c      endif

c      if (iabsorp.eq. 1) then
c         id=4
c         call flown (u(1),u(2),u(3),u(4),u(5),u(6),cflown,dwepsdw)
c         id=id_old
c      else !uses the dispersion with the given 'id' 
c         call flown (u(1),u(2),u(3),u(4),u(5),u(6),
c     1                   cflown,dwepsdw)
c      endif

c      write(*,*)'prep3d before flown iflux',iflux
           
      if (iflux.eq.1) then
cSAP20605
          call flown (u(1),u(2),u(3),u(4),u(5),u(6),cflown,dwepsdw)
c          write(*,*)'prep3d iflux=1 ater call flown cflown',cflown
      endif

      if (iflux.eq.2) then
c------- flux from the cold plasma    
         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
         rnpar=cnpar
         rnper=cnper
         romega=1.d0/ye
         romegpe=dsqrt(xe)/ye
cyup        write(*,*)'prep3d before grpde xe,ye',xe,ye
         nsigma=ioxm
c        nsigma=1
c        nsigma=-1
         call grpde2 (rnpar, rnper, romega, romegpe, nsigma,
     .                   rvgrpdc, redenfac, exde, eyde, ezde)
c        write(*,*)'+ redenfac',redenfac
         cflown=2.d0*redenfac
         if(rvgrpdc.le.0.d0)then
         !YuP[2018-05-23]added: stop the ray when vgrpdc<0 from grpde2()
            iraystop=1
            iray_status_one_ray=11 ! means: ray is stopped by fluxn<0
         endif
c         write(*,*)'prep3d cold from grpde2 cflown',cflown    
      endif !cold electron plasma flux
    
c-------------------------------------------------------
c     deru(i) -dimensionless
c      write(*,*)'in prep3d before rside1 nrayelt',nrayelt

c-----calculate tokman flux
      if (((id.eq.10).or.(id.eq.12).or.(id.eq.13)).or.(id.eq.15)) then
        call flux_tokman(z,r,phi,cnz,cnr,cm,flux_tokman_ar)
        wflux_tokman_z(is)=flux_tokman_ar(1)
        wflux_tokman_r(is)=flux_tokman_ar(2)
        wflux_tokman_phi(is)=flux_tokman_ar(3)
      endif


c-----calculate the group velocity
      id_old=id

c      if (iabsorp.eq.1) then
c         id=4 
c---------the group velocity from Mazzucato dispersion
c         call rside1(t,u,deru)
c
c         write(*,*)'id=4,dsqrt(deru(1)**2+deru(2)**2+(r*deru(3))**2)',
c     .   dsqrt(deru(1)**2+deru(2)**2+(r*deru(3))**2)
c
c         id=id_old
c      else
c         call rside1(t,u,deru)
c      endif


      i_geom_optic_loc=i_geom_optic
      i_geom_optic=1 !to get group velocity in deru 
      call rside1(t,u,deru)
      i_geom_optic=i_geom_optic_loc

      vgr(1)=deru(1)
      vgr(2)=deru(2)
      vgr(3)=r*deru(3)
      
      vgrmods=vgr(1)**2+vgr(2)**2+vgr(3)**2
      if(vgrmods.gt.1.01d0)then !YuP[2020-08-17] Changed gt.1 to gt.1.01
         !otherwise too much of printout; 
         ! vgroup=1.00005 (slightly larger than 1) might happen in near-vacuum
         write(*,*)
         write(*,*) '*************************************************'
       write(*,*)'prep3d: WARNING: vgroup>1, abs(vgroup)=',sqrt(vgrmods)
         write(*,*) '*************************************************'
         write(*,*)
      endif

c      write(*,*)'prep3d vgr',vgr(1),vgr(2),vgr(3)
c      write(*,*)'prep3d vgrmods',vgrmods

c-----the data for mnemonic.nc output file
c     the group velocity normalized to c
      wvgr_z(is)=vgr(1)
      wvgr_r(is)=vgr(2)
      wvgr_phi(is)=vgr(3)
c-----refractive index
      wn_z(is)=cnz
      wn_r(is)=cnr
      wn_phi(is)=cm/r
      
c-----------------------------------------------------------
c     vdotb -projection of the group velocity  on the
c            magnetic field multiplited by bmod
c-----------------------------------------------------------
      vdotb=vgr(1)*bz+vgr(2)*br+vgr(3)*bphi
c-----------------------------------------------------------
c     vgperps -perpendicular (to magnetic field)
c              component of the group velocity
c              in the second degree
c-----------------------------------------------------------
      vgperps=0.0d0
      do i=1,3
        vgperps=vgperps+(vgr(i)-vdotb*bf(i)/bmod**2)**2
      enddo
cSm030514
c      write(*,*)'dsqrt(vgrperps)',dsqrt(vgperps)
c----------------------------------------------------------
c     collisional damping 
c----------------------------------------------------------
      if( iabsorp_collisional.eq.1) then
         temp_e=tempe(z,r,phi,1)
         dens_e=dense(z,r,phi,1)
         z_eff=zeff(z,r,phi)
         frqncy_l=frqncy
         v_gr_perp= dsqrt(vgperps)
         call absorp_collisional(temp_e,dens_e,frqncy_l,z_eff,
     &   v_gr_perp,coll_mult,
     &   cnprim_cl)
cyup         write(*,*)'prep3d cnprim,cnprim_cl',cnprim,cnprim_cl
 2       cnprim=cnprim+cnprim_cl
cSAP090403
c         cnprim_e=cnprim_e+cnprim_cl
      endif !iabsorp_collisional=1 
c----------------------------------------------------------
      vgrs=vgr(1)**2+vgr(2)**2+vgr(3)**2
      ckvi=dsqrt(vgperps/vgrmods)*cnprim/cld
      vgrpls=vgr(1)**2+vgr(2)**2
      vgrpol=dsqrt(vgrpls)
      wf=frqncy
c----------------------------------------------------------
c     ckvipol  (1/cm)
      vratio=dsqrt(vgperps/vgrpls)
      ckvipol=vratio*cnprim/cld
c      write(*,*)'prep3d, ws,ckvipol with grp vel:',ws(is),ckvipol
cSmirnov970101 beg
      ckvipl_e=vratio*cnprim_e/cld
      ckvipl_i=vratio*cnprim_i/cld
      ckvipl_cl=vratio*cnprim_cl/cld
 
cyup      write(*,*)'prep3d: ckvipl_e,ckvipl_i,ckvipl_cl',
cyup     &                   ckvipl_e,ckvipl_i,ckvipl_cl


cBH041009  Only germaine if iabsorp.eq.3:
      do kk=2,nbulk
         ckvipl_s(kk)=vratio*cnprim_s(kk)/cld
      enddo

cSm040911
c     for id=10,12,13 vgperps aand  vgrpls have not sense
c     we will use Tokman flux
      
      if (((id.eq.10).or.(id.eq.12).or.(id.eq.13)).or.(id.eq.15)) then  
c--------------------------------------------------------------
c        flux_tokman_b -projection of the tokman flux on the
c                    magnetic field multiplited by bmod
c--------------------------------------------------------------
         flux_tokman_b=flux_tokman_ar(1)*bz+flux_tokman_ar(2)*br+
     &                 flux_tokman_ar(3)*bphi

c-----------------------------------------------------------
c        flux_tokman_perp_s -perpendicular (to magnetic field)
c              component of the tokman flux
c              in the second degree
c-----------------------------------------------------------
         flux_tokman_perp_s=0.d0
         do i=1,3
            flux_tokman_perp_s=flux_tokman_perp_s+
     &      (flux_tokman_ar(i)-flux_tokman_b*bf(i)/bmod**2)**2
         enddo

         flux_tokman_s=flux_tokman_ar(1)**2+flux_tokman_ar(2)**2+
     &                 flux_tokman_ar(3)**2

         flux_tokman_mod=sqrt(flux_tokman_s)

         ckvi=dsqrt(flux_tokman_perp_s/flux_tokman_s)*cnprim/cld

         flux_tokman_pol_s=flux_tokman_ar(1)**2+flux_tokman_ar(2)**2
         flux_tokman_pol=dsqrt(flux_tokman_pol_s)
         ratio_perp_pol=dsqrt(flux_tokman_perp_s/flux_tokman_pol_s)

         ckvipl_e=ratio_perp_pol*  cnprim_e/cld
         ckvipl_i=ratio_perp_pol*  cnprim_i/cld
         ckvipl_cl=ratio_perp_pol*  cnprim_cl/cld
         ckvipol=ratio_perp_pol*  cnprim/cld
cBH041009  Only germaine if iabsorp.eq.3:
         do kk=2,nbulk
            ckvipl_s(kk)=ratio_perp_pol*cnprim_s(kk)/cld
         enddo

cyup         write(*,*)'prep3d, ckvipol with tokman flux ratio,ckvipl_e',
cyup     &   ckvipl_e
         
      endif

c      write(*,*)'prep3d: ckvipl_s=',ckvipl_s
c---------------------------------------------------------
c     to check eigenvalues
c      call eigenvalue_hot(z,r,phi,cnpar,cnper,cnprim,
c     &eigenvalue)
c      do i=1,3
c        write(*,*)'i, eigenvalue(i)',i, eigenvalue(i)
c      enddo

c--------------------------------------------------------
cSmirnov970101 end
c----------------------------------------------------------
      seikon(is)=0.d0
      spsi(is)=rho/a
c---------------------------------------------------------
c     wr and wz (cm),wphi (radian)
      wr(is)=r*r00
      wphi(is)=phi
      wz(is)=z*r00

      wnpar(is)=cnpar
c      write(*,*)'prep3d wnpar(is)',wnpar(is)
      wnper(is)=cnper
c      write(*,*)'prep3d is,cnper,wnper(is)',is,cnper,wnper(is)
      wmtor(is)=u(6)

c-------------------------------------------------------------------
c     Here delpwr(is) is the power(erg/c) in the ray
c     channel.It is equal powini_iray at antenna.
      if (is.eq.1) then    !ends line 2255
c        powj(iray) (erg/c) was calculated in cone_ec
c        powinilh(iray) (erg/c) was calculated in grill_lh
c        powini=powj or powinilh
         p=powini
cyup         write(*,*)'prep3d is=1, powini=',powini
cSAP081202
c        delpwr(is)=dexp(-2.d0*ckvipol*(ws(is)))*p
         delpwr(is)=p
c         write(*,*)'is=1 delpwr(1)',delpwr(1)
         optical_depth=0.d0
cSm021101
         t_old=0.d0
         iflref_old=0

cSAP081111
         ckvipold=ckvipol
cRWH140909  Stop ray immediately, if has no starting power
cRWH160227  except for i_ox=1 calculation for rays with optimal OX 
cRWH160227  mode converison.
         if (dabs(powini).lt.1.d-100 .and. i_ox.ne.1) then
            iraystop=1
            iray_status_one_ray=14 !means: abs(powini).lt.1.d-100
         endif

      else  !On is.eq.1, that is, else is.gt.1, begins line 2023
                                               ! ends line 2255

        if (iabsorp_ql.eq.0) then  !ends line 2210
c----------do not use QL flux for absorption 
cSmirnov970105
c          argexp=(-2.d0*ckvipol*(ws(is)-ws(is-1)))
cSAP081111
c	   ckvipold=0.5d0*(sdpwr(is-1)+salphal(is-1)+salphac(is-1))

c           write(*,*)'ckvipold',ckvipold
c           write(*,*)'0.5d0*(sdpwr(is-1)+salphal(is-1)+salphac(is-1))',
c     &               0.5d0*(sdpwr(is-1)+salphal(is-1)+salphac(is-1))

           argexp= -(ckvipol+ckvipold)*(ws(is)-ws(is-1)) 

c           write(*,*)'ws(is)-ws(is-1)',ws(is)-ws(is-1)
c           write(*,*)'(ckvipol+ckvipold),argexp',
c     &                (ckvipol+ckvipold),argexp

cSAP100202         
cyup           write(*,*)'prep3d.f dabs(1.d0-dexp(argexp)),eps_delta_pow',
cyup     &                         dabs(1.d0-dexp(argexp)),eps_delta_pow

           i_go_to_previous_out_step=0   

           if ((is.gt.1).and.
cSAP100416
c     &         (delpwr(is-1).gt.delpwr(1)*eps_delta_pow)) then
     &         (delpwr(is-1).gt.delpwr(1)*eps_delta_pow*1.d-2)) then

              if (dabs(1.d0-dexp(argexp)).gt.eps_delta_pow) then
c--------------------------------------------------------------
c                If power variation for one output step is bigger
c                than eps_delta_pow then the code will reduce 
c                output step prmt(6) and will return to the previous 
c                output point
c--------------------------------------------------------------
                 i_go_to_previous_out_step=1

c                 write(*,*)'prep3d.f is,ws(is),ws(is-1)',
c     &                            is,ws(is),ws(is-1)

c                 write(*,*)'prep3d.f argexp,(ws(is)-ws(is-1))',
c     &                            argexp,(ws(is)-ws(is-1))

c                 write(*,*)'ckvipol+ckvipold',ckvipol+ckvipold

                 prmt6_new=prmt(6)/2.d0
                 zold=wz(is-1)*0.01d0/r0x
                 rold=wr(is-1)*0.01d0/r0x
                 phiold=wphi(is-1)
cyup                 write(*,*)'prep3d.f argexp,prmt6_new',argexp,prmt6_new
cyup                 write(*,*)'prep3d.f
cyup     &                      prmt6_new*(ckvipol+ckvipold)*r0x*1.d2'
cyup     &                     ,prmt6_new*(ckvipol+ckvipold)*r0x*1.d2
                 return
              else
c                 if (dabs(1.d0-dexp(argexp)).lt.0.5*eps_delta_pow) then
c--------------------It was too small step prmt(6).                  
c                     i_go_to_previous_out_step=0  
c                     write(*,*)'prep3d.f small prmt(6),prmt6',
c     &                                         prmt(6),prmt6
c                     prmt6_new=min(prmt(6)*2.d0,prmt6) 
c                     write(*,*)'prep3d prt6_new',prt6_new
c                     prmt(6)=prmt6_new
c                     write(*,*)'prep3d.f small step, new prmt(6)',
c     &                          prmt(6)
c                 else    
                     i_go_to_previous_out_step=0
c                 endif               
              endif !On dabs(1.d0-dexp(argexp)).gt.eps_delta_pow
           else  !On is
             prmt(6)=prmt6
             i_go_to_previous_out_step=0
           endif !On (is.gt.1).and.(delpwr(is-1).gt.delpwr(1)...[l 2050]
                 
cSAPO81111
           ckvipold=ckvipol
cSmirnov970105
c          write(*,*)'in prep3d ckvipol,argexp',ckvipol,argexp

c	   if (dabs(argexp).gt.90.d0) then
c	     write(*,*)'in prep3d argexp.gt.90',argexp
c             delpwr(is)=0.d0
c	   else 
c             pwexp=dexp(argexp)
c	     if (pwexp.lt.1.d-50) then
c	         write(*,*)'in prep3d pwexp.lt.1.d-50',pwexp
c                 delpwr(is)=0.d0
c	      endif
c	   endif

           optical_depth=optical_depth+dabs(argexp)
c          write(*,*)'optical_depth',optical_depth
                  
cSm050225
c-YuP-130605: changed ray-stopping criterion from argexp>0.d0 to this:
           if(argexp.gt. 1.d-30)then
              WRITE(*,*)'**********************************************'
              WRITE(*,*)'WARNING: in prep3d.f argexp>0 :' , argexp
              WRITE(*,*)'It would give growing ray power-> stopping ray'
              WRITE(*,*)'ckvipol,ckvipold=',ckvipol,ckvipold
              WRITE(*,*)'**********************************************'
              argexp=0.d0
              iraystop=1
              iray_status_one_ray=10 !means: argexp>1d-30 (growing ray power)
           endif
c-YuP-130605: From print-out: 
c Even though sometimes ckvipol and ckvipold both are zero,
c yet the value of  argexp= -(ckvipol+ckvipold)*(ws(is)-ws(is-1))
c is not zero (but rather a small number ~ 1e-321).
c Because of this seemingly insignificant rounding error,
c the rays were stopped prematurely.
c It only happens on Hopper/PGI, not on IntelVisualFortran.
cSm021101
c          write(*,*)'delpwr(is-1),argexp',delpwr(is-1),argexp
cSAP100227
cSAP110114
           if(i_lsc_approach.eq.1) then
cBH150117    Somehow, delpwr(is-1) in following appears to have
cBH150117    attained a zero value, giving a divide by zero.
cBH150117    Error located by Rob Andre; only known to occur for
cBH150117    the PPPL pathscale compiler.
cBH150117    if (argexp.lt.dlog(delpwrmn*delpwr(1)/delpwr(is-1))) then
cBH150117        delpwr(is)=0.d0         
cBH150117    else
cBH150117        delpwr(is)=delpwr(is-1)*dexp(argexp)
cBH150117    endif
             if (delpwr(is-1).gt.0d0) then
               if (argexp.lt.dlog(delpwrmn*delpwr(1)/delpwr(is-1))) then   
                 delpwr(is)=0.d0
               else
                 delpwr(is)=delpwr(is-1)*dexp(argexp)
               endif
             else
               delpwr(is)=0.d0
             endif
           else  !On i_lsc_approach
             delpwr(is)=delpwr(is-1)*dexp(argexp)
           endif  !On i_lsc_approach

c          delpwr(is)=delpwr(is-1)*dexp(-2.d0*dabs(wi)*
c     &               (ws(is)-ws(is-1))/(dsqrt(vgrs)*cvac)*    
c     &              (2.d0*pi*frqncy*1.d+9))

         else  !On if(iabsorp_ql.eq.0), line 2050
c----------to use QL flux for absorption calculations
c          iabsorp_ql=1
cyup           write(*,*)'prep3d before absorbed_power_using_ql_flux'
c           call tensrcld(u(1),u(2),u(3))
c           cnx=dcmplx(cnper,0.d0)
c           call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)

           call absorbed_power_using_ql_flux(wnpar(is-1),wnper(is-1),
     &     wz(is-1),wr(is-1),wphi(is-1),
     &     fluxn(is-1),delpwr(is-1),(ws(is)-ws(is-1)),
     &     absorbed_power_ql)
           delpwr(is)=delpwr(is-1)-absorbed_power_ql
           
           if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
           write(*,*)'QL absorption'
           write(*,*)'delpwr(is-1),absorbed_power_ql,delpwr(is)',
     &                delpwr(is-1),absorbed_power_ql,delpwr(is)
           endif ! outprint

         endif !iabsorp_ql, begins at line 2050
 
c-------------------------------------------------------------------
c        reflection lost at the plasma edge
c------------------------------------------------------------------
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'prep3d refl_loss,irefl,irefl_old',
     &             refl_loss,irefl,irefl_old
         endif ! outprint
         
         tot_pow_absorb_at_refl=tot_pow_absorb_at_refl+
     &           delpwr(is)*refl_loss*(irefl-irefl_old)

cyup         write(*,*)'prep3d before refl_looss delpwr(i)',delpwr(is)

         delpwr(is)=delpwr(is)*(1.d0-refl_loss*(irefl-irefl_old))
         irefl_old=irefl 
         w_tot_pow_absorb_at_refl=tot_pow_absorb_at_refl
c-------------------------------------------------------------------    
          
cSAP090603
c         write(*,*)'prep3d delpwr(is),delpwr(is-1)-delpwr(is)',
c     &                     delpwr(is),delpwr(is-1)-delpwr(is)

c         write(*,*)'delpwr(1)',delpwr(1)

cSAP090603
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         if((delpwr(is).gt.1.d-200).and.(delpwr(1).gt.1.d-200))then
           if(i_ox.ne.1) write(*,*)'-dlog(delpwr(is)/delpwr(1))',
     &             -dlog(delpwr(is)/delpwr(1))
         endif
         endif ! outprint
 
c         if(argexp.gt.0.d0)then
c           write(*,*)'******************************************'
c           write(*,*)'WARNING: in prep3d.f argexp>0' 
c           write(*,*)'It will give the growing ray power'
c          write(*,*)'******************************************'
c         endif

         if(i_lsc_approach.eq.0) then        
          if(delpwr(is).lt.delpwrmn*delpwr(1))then
            WRITE(*,*)'**prep3d delpwr(is)<delpwrmn*delpwr(1). iray,is='
     +      ,iray,is     
c           stop ray_iray calculations
            iraystop=1
            iray_status_one_ray=2 ! delpwr(is).lt.delpwrmn*delpwr(1)
          endif
         endif 
      end if  !On is.eq.1 [l 2009], else [l 2026]

      if(i_ox.eq.2) then
c---------------------------------------------------------------
c       It works for i_ox=2 case after OX mode conversion point,
c       where i_ox_conversion=1
c       It will reduce the power from O mode to X mode using 
c       transmission coefficient transm_ox

cyup        write(*,*)'i_call_prep3d_in_output_at_i_ox_conversion_eq_1',
cyup     &  i_call_prep3d_in_output_at_i_ox_conversion_eq_1

        if (i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.0)goto 20

c        write(*,*)'prep3d is,delpwr(is),delpwr_o',
c     &  is,delpwr(is),delpwr_o

        if (i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.1) then
          nrayelt=nrayelt-1
          is=is-1
          nrayelt_o_cutoff=is  !the number of ray point where O
                               ! cutoff was found
cyup          write(*,*)'prep3d: nrayelt_o_cutoff',nrayelt_o_cutoff
          goto 20 !the first call of prep3d in output after OX conversion 
        endif

        if (i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.2) then        
          transm_ox_loc=transm_ox
c          write(*,*)'prep3d OX transm i_ox_conversion,delpwr(is),
c     &    transm_ox',i_ox_conversion,delpwr(is),transm_ox
          delpwr_o=delpwr(is-1) !o-mode before jump
c          call OX_power_transmission(is,i_ox,i_ox_conversion,
c     &    delpwr_o,transm_ox_loc,delpwr_x)
          delpwr_x=delpwr_o*transm_ox
c          write(*,*)'prep3d delpwr_o,transm_ox_loc,delpwr_x',
c     &                      delpwr_o,transm_ox_loc,delpwr_x
          delpwr(is)=delpwr_x   !x-mode after jump
        endif

 20   continue
  
cyup      write(*,*)'prep3d after transm delpwr(is)',delpwr(is)
      endif ! i_ox.eq.2

c     sdpwr(is)=0.d0
c----------------------------------------------------------------------
cSAP100210 wdnpar(is) was set in prep3d.f. See Line 220
cSmirnov961122
c      if(istart.eq.1) then
c        electron cyclotron waves
c         wdnpar(is)=dabs(0.05d0*wnpar(1))
c        write(*,*)'in prep3d old wdnpar(is)',wdnpar(is)
c      else
cSmirnov961205
c        grill conditions for the waves (LH or FW)
c         wdnpar(is)=wdnpar0(iray)
c        write(*,*)'in prep3d new wdnpar(is)',wdnpar(is)
c      endif

      cwexde(is)=cex
      cweyde(is)=cey
      cwezde(is)=cez
c---------------------------------------------------------
c     vgrpol(cm/sec)=vgrpol*r00/t00
c     vrgpol/c=vgrpol/wf
c---------------------------------------

c     for test
c      ppp=vgrpol/wf
c      ppp1=dsqrt(vgrmods)/wf
c       write(*,*)'vgrpol=',ppp,'dsqrt(vgrmods)',ppp1,'wf',wf
c     for test	end
c---------------------------------------
      wf=frqncy
c      fluxn(is)=dreal(cflown)*vgrpol*0.5d0/wf
c020312
c      fluxn(is)=dreal(cflown)*vgrpol*0.5d0
CSAP081122 
      p_flux=cflown*dconjg(cflown)    
      fluxn(is)=dsqrt(p_flux)*vgrpol*0.5d0

      if (((id.eq.10).or.(id.eq.12).or.(id.eq.13)).or.(id.eq.15)) then
        fluxn(is)=dsqrt(flux_tokman_ar(1)**2+
     &   flux_tokman_ar(2)**2)
      endif

c       write(*,*)'prep3d dreal(cflown)',dreal(cflown)
c       write(*,*)'prep3d.f vgrpol,wf', vgrpol,wf
c       write(*,*)'vgrpol',vgrpol,'fluxn(is)',fluxn(is)

c----------------------------------------------------------------------
c   Stop ray if flux.lt.0.0 (set fluxn previous value)  [RWH:030428]
c----------------------------------------------------------------------
         if(fluxn(is).lt.0.0)then
            write(*,*)
            write(*,*) '***********************************************'
            write(*,*) '***in prep3d fluxn(is).lt.0.0. Set iraystop=1'
            write(*,*) '***in prep3d Set fluxn(is)=fluxn(is-1)'
            write(*,*) '***********************************************'
            write(*,*)
            if(is.gt.1) then
               fluxn(is)=fluxn(is-1)
            else
               fluxn(is)=1.
            endif
c           stop ray_iray calculations
cSm051130
            iraystop=1
cSAP100514
            iray_status_one_ray=11 ! means: fluxn<0
         endif
     
c-----------------------
cBH040816      sbtot(is)=bmod*10000.d0*b0
      one=1.d0
      sbtot(is)=bmod*10000.d0*b0*dsign(one,feqd(1))
cBH040915:  Magnetic field components
      sb_z(is)=1.d4*bz
      sb_r(is)=1.d4*br
      sb_phi(is)=1.d4*bphi
c-----------------------------------------------------------
c     dense - dimensionless
      sene(is)=dense(z,r,phi,1)*1.0d+13
      ste(is)=tempe(z,r,phi,1)
      szeff(is)=zeff(z,r,phi)
c     salphac(is)=0.0d0
c     salphal(is)=2.d0*ckvipol
c Smirnov970105 beg
c BH991017   sdpwr(is)=2.d0*ckvipl_e   ! electron damping coefficient
c BH991017   salphal(is)=2.d0*ckvipl_i  ! ion damping coefficient
      salphac(is)=2.d0*ckvipl_cl  ! collisional damping coefficient
c Smirnov970105 end
cHarvey991017 beg
      sdpwr(is)=  2.d0*ckvipl_i     ! ion damping coefficient
      salphal(is)=2.d0*ckvipl_e     ! electron damping coefficient

c      write(*,*)'prep3d salphal(is)',salphal(is)

cBH041009
      do kk=2,nbulk
         salphas(is,kk)=2.d0*ckvipl_s(kk)
c         write(*,*)'is,kk,salphas(is,kk)',is,kk,salphas(is,kk)
      enddo
cHarvey991017 end
c------------------------------------------------------------
c     for xdraw plotter
      xarr(is)=r*dcos(phi)*r00
      yarr(is)=r*dsin(phi)*r00
      rez(is)=cdabs(cez)
c     if (istart.eq.2) then
c        start point is inside the plasma (for LH and FW)
c        using cold plasma dielectric tensor
c        rez(is)=ezcold
c     end if
c------------------------------------------------------------

      if (i_emission.eq.1) then
c--------calculate the data for emission

         if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c----------usage of the analytical relativistic function and its derivatives 
           i_fkin=0
         else 
c----------usage of the mech relativistic function and its derivatives
           i_fkin=1
         endif
                          
         wr_em(is)=wr(is)/r00     !for emission
         wz_em(is)=wz(is)/r00     !for emission
         wphi_em(is)=wphi(is) !for emission
        
         xe=x(z,r,wphi,1)
         ye=y(z,r,phi,1) !it will be used as negative for electron
         T_kev=tempe(z,r,phi,1)         
         wtemp_em(is)=T_kev  
         vgroup=dsqrt(vgrmods) !in clightc

c--------calculate cnray

c         omegpedce=dsqrt(xe)/dabs(ye) !omega_pe/omega_ce
c         omegdce=1.d0/dabs(ye)        !omega/omega_ce
c         cnray=rrind(wnpar(is),wnper(is),omegdce,omegpedce)
c
c         cnray=dsqrt(wnpar(ia)**2+wnper(is)**2)

c--------initialize arrays
         do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
         enddo

c         cnray1=rrind_hotpl(nbulk,dmas,x_ar,y_ar,
c     .   t_av_ar,tpop_ar,vflow_ar,cnpar,cnper)

c         write(*,*)'prep3d before rayrind'

         call rayrind(i_rrind,nbulk,dmas,x_ar,y_ar,
     .   t_av_ar,tpop_ar,vflow_ar,cnpar,cnper,cnray)
         
c         write(*,*)'prep3d after rayrind cnray',cnray

         cnn=dsqrt(cnpar**2+cnper**2)
           
         wnray(is)=cnray
         if (is.eq.1) then
c           write(*,*)'is,wnray(is)',is,wnray(is)
         endif

         wcnz_em(is)=cnz 
         wcnr_em(is)=cnr
         wcm_em(is)=cm

c         write(*,*)'emis is,xe,ye,T_kev',is,xe,ye,T_kev
c         write(*,*)'wnpar(is),wnper(is),cnray',
c     +   wnpar(is),wnper(is),cnray
c         write(*,*)'n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin',
c     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin 
   
c         write(*,*)'r,z,phi',r,z,phi
c         write(*,*)'cex,cey,cez',cex,cey,cez
c         write(*,*)'cflown,vgroup,frqncy',cflown,vgroup,frqncy

c-------- electric field using the cold plasma dielectric tensor
c         call tensrcld(u(1),u(2),u(3))
c         cnx=dcmplx(cnper,0.d0)
c         call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c         cwexde(is)=cex
c         cweyde(is)=cey
c         cwezde(is)=cez
c         rez(is)=cdabs(cez)
c-------------------------
c         write(*,*)'prep3d before calc_emis_coef cex,cey,cez',
c     &   cex,cey,cez

c         write(*,*)'xe,-ye,T_kev,wnpar(is),wnper(is),cnray,
c     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,phi,
c     +   cex,cey,cez,cflown,   data tot_pow_absorb_at_refl/0.d0/vgroup,frqncy',
c     +   xe,-ye,T_kev,wnpar(is),wnper(is),cnray,
c     +   n_relt_harm1,n_relt_harm2,,n_relt_intgr,i_fkin,r,z,phi,
c     +   cex,cey,cez,cflown,vgroup,frqncy

c--------the data for collisional absorption
         iabsorp_collisional_l=iabsorp_collisional
         coll_mult_l=coll_mult
         tempe_l=tempe(z,r,phi,1) ! kev
         dense_l=dense(z,r,phi,1) ! kev
cSm060725
c         zeff_l=zeff(z,r,phi,1) ! kev
cBH060726         zeff_l=zeff(z,r,phi,1) ! kev
         zeff_l=zeff(z,r,phi) ! kev

         call cpu_time(time_prep3d_emis_1)

cSm060314
         call calc_emis_coef(xe,-ye,T_kev,wnpar(is),wnper(is),cnray,
c         call calc_emis_coef(xe,ye,T_kev,wnpar(is),wnper(is),cnray,
     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,  
     +   i_resonance_curve_integration_method,epsi,
     +   i_fkin,r,z,phi,
     +   cex,cey,cez,cflown,vgroup,frqncy,
     +   iabsorp_collisional_l,coll_mult_l,tempe_l,dense_l,zeff_l,
     +   wal_emis(is),wj_emis(is),wtemp_rad_em(is))
c         write(*,*)'prep3d after calc_emis_coef'
cyup         write(*,*)'is wal_emis,wj_emis',is,wal_emis(is),wj_emis(is)
cyup         write(*,*)'wtemp_em(is),wtemp_rad_em(is)',
cyup     +   wtemp_em(is),wtemp_rad_em(is)

         call cpu_time(time_prep3d_emis_2)

c----------------------------------------------------
c         calculate data to plot the resonance ellipse boundary
c         for several harmonics 
         y_e=y(z,r,phi,1)
         clight=2.99792458d10
         rme=9.1094d-28 !electron mass [gram]
         temp_e=tempe(z,r,phi,1) !kev
         vt=dsqrt(temp_e*1.6022d-9/rme) !(cm/sec) thermal velocity 
c          write(*,*)'is,r,y_e,temp_e,cnpar',is,r,y_e,temp_e,cnpar
c         write(*,*)'is,vt,clight,vt/clight,',is,vt,clight,vt/clight
         do n=n_relt_harm1,n_relt_harm2 
             n_negative=-n  
cSm060315
c             call ec_cond(n,y_e,cnpar,ires,p_perp0)
             call ec_cond(n_negative,-y_e,cnpar,ires,p_perp0)
             w_ires(is,n)=ires
c             call root_res(0.d0,cnpar,n,y_e,kmax,p_par_rl)
cSm060315
c             call root_res_test(0.d0,cnpar,n,y_e,kmax,p_par_rl,
            call root_res_test(0.d0,cnpar,n_negative,-y_e,kmax,p_par_rl,
     &       p0,del_p)

cSm060303
             if(ires.eq.1)then
              wp_0_dmvt(is,n)=p0*clight/vt
              wdel_p_dmvt(is,n)=del_p*clight/vt
             else
              wp_0_dmvt(is,n)=0.d0
              wdel_p_dmvt(is,n)=0.d0
             endif

             if (ires.eq.1) then
               wp_perpmax_dmvt(is,n)=p_perp0
     &         *clight/vt
             else
               wp_perpmax_dmvt(is,n)=0.d0
             endif
 
             if(kmax.gt.0) then
               wp_parmin_dmvt(is,n)=dmin1(p_par_rl(1),p_par_rl(2))
     &         *clight/vt
               wp_parmax_dmvt(is,n)=dmax1(p_par_rl(1),p_par_rl(2))
     &         *clight/vt
             else
               wp_parmin_dmvt(is,n)=0.d0
               wp_parmax_dmvt(is,n)=0.d0
             endif
             
c            write(*,*)'vt,clight,p_perp0,p_par_rl',
c     &      vt,clight,p_perp0,p_par_rl

c            write(*,*)'n,ires',n,ires
c            write(*,*)'clight/vt',clight/vt
c            write(*,*)'cnpar,y_e,n*y_e',cnpar,y_e,n*y_e
c            write(*,*)'wp_0_dmvt(is,n)',wp_0_dmvt(is,n)
c            write(*,*)'wdel_p_dmvt(is,n)',wdel_p_dmvt(is,n)

c            write(*,*)'wp_perpmax_dmvt(is,n)', wp_perpmax_dmvt(is,n)
c            write(*,*)'wp_parmin_dmvt(is,n),wp_parmax_dmvt(is,n)',
c     &      wp_parmin_dmvt(is,n),wp_parmax_dmvt(is,n)

             wp_perpmax_dmc(is,n)=wp_perpmax_dmvt(is,n)*vt/clight
             wp_parmin_dmc(is,n)= wp_parmin_dmvt(is,n)*vt/clight
             wp_parmax_dmc(is,n)= wp_parmax_dmvt(is,n)*vt/clight
         enddo

  
c         write(*,*)'wp_parmin_dmc(is,2),wp_parmin_dmvt(is,2)',
c     &   wp_parmin_dmc(is,2),wp_parmin_dmvt(is,2)


c--------data for test like horace emission
c         clight=2.99792458d10
c         rme=9.1094d-28 
c         tem0=5.d0 !keV
c         tem0=temperho(0.d0,1)
c         vt0=dsqrt(2.d0*tem0*1.6022d-9/rme) !(cm/sec)the  central thermal velocity 
c         cdvt=clight/vt0
c         vmaxdt=9.5d0
c         vmaxdt=100.5d0
         
c         call calc_emis_coefh(xe,ye,T_kev,wnpar(is),wnper(is),cnray,
c     +   cdvt,vmaxdt,
c     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,phi,
c     +   cex,cey,cez,cflown,vgroup,frqncy,
c     +   iabsorp_collisional_l,coll_mult_l,tempe_l,dense_l,zeff_l,
c     +   wal_emis(is),wj_emis(is),wtemp_rad_em(is))

c         write(*,*)'wal_emis,wj_emis',is,wal_emis(is),wj_emis(is)
c         write(*,*)'wtemp_em(is),wtemp_rad_em(is)',
c     +   wtemp_em(is),wtemp_rad_em(is)

      goto 15
      if (is.eq.1)then !test genray and like horace emission coef.
           call maxw_test !check the density from maxwellian distribution
c---------------------------------------------------------          
c           variant close to ns=18 in horace
c           ye=0.992653766d0
c           xe=0.199603072d0          
c           cnll=0.147982438d0
c           wnpar(1)=cnll
c           wnper(1)=0.884755135d0
c           ckper=0.124d0
c           n=1
c           T_kev=4.99854639
c           cex=dcmplx(-0.132492995d0,0.00110550533d0)
c           cey=dcmplx(-0.0147405718d0,0.0840407114d0)
c           cez=dcmplx(0.987504055d0,0.d0)
c----------------------------------------------------------
c           variant close to ns=23 in horace
c           ye=1.00802004d0
c           xe=0.199594029d0     
c           cnll=0.150273259d0
c           wnpar(1)=cnll
c           wnper(1)=0.884450246d0
c           ckper=0.124d0
c           n=1
c           T_kev=4.99831992d0
c           cex=dcmplx(-0.131414611d0,-0.000980708504d0)
c           cey=dcmplx(-0.0460538893d0,0.101026237d0)
c           cez=dcmplx(0.985089832d0,0.)
c           write(*,*)'cflown,vgroup',cflown,vgroup
c           cflown=2.d0
c           vgroup=0.896719192d0
c----------------------------------------------------------
         call calc_emis_coefh(xe,ye,T_kev,wnpar(is),wnper(is),cnray,
     +   cdvt,vmaxdt,
     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,phi,
     +   cex,cey,cez,cflown,vgroup,frqncy,
     +   iabsorp_collisional_l,coll_mult_l,tempe_l,dense_l,zeff_l,
     +   wal_emis(is),wj_emis(is),wtemp_rad_em(is))

c         call calc_emis_coef(xe,-ye,T_kev,wnpar(is),wnper(is),cnray,
c     +   n_relt_harm1,n_relt_harm2,n_relt_intgr, 
c     +   i_resonance_curve_integration_method,
c     +   i_fkin,r,z,phi,
c     +   cex,cey,cez,cflown,vgroup,frqncy,
c     +   iabsorp_collisional_l,coll_mult_l,tempe_l,dense_l,zeff_l,
c     +   wal_emis(is),wj_emis(is),wtemp_rad_em(is))

        stop
      endif
cendtes
 15   continue

c         call calc_emis_coefh(xe,ye,T_kev,wnpar(is),wnper(is),cnray,
c     +   cdvt,vmaxdt,
c     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,r,z,phi,
c     +   cex,cey,cez,cflown,vgroup,frqncy,
c     +   iabsorp_collisional_l,coll_mult_l,tempe_l,dense_l,zeff_l,
c     +   wal_emis(is),wj_emis(is),wtemp_rad_em(is))

c!!!!!for test
c       itest=1    

c-----normalised r and z
c      n_midway=is 
c      r=wr_em(n_midway)
c      z=wz_em(n_midway)
c      phi=wphi_em(n_midway)
c      write(*,*)'wr_em(n_midway),r',wr_em(n_midway),r
c      write(*,*)'wz_em(n_midway),z',wz_em(n_midway),z
c      write(*,*)'wphi_em(n_midway),phi',wphi_em(n_midway),phi
c      cnz=wcnz_em(n_midway)
c      cnr=wcnr_em(n_midway)
c      cm=wcm_em(n_midway)

c-----emission and absorption in the midway point
 
c-----calculate the data for emission
c      bmod=b(z,r,phi)
c      write(*,*)'bmod',bmod
c      xe=x(z,r,phi,1)
c      ye=y(z,r,phi,1)    !it will be used as negative for electron
c      write(*,*)'emission.f in add_ray_point xe,ye',xe,ye
c      T_kev=tempe(z,r,phi,1)  
c      wtemp_em(n_midway)=T_kev ! for plotting

c      u(1)=z
c      u(2)=r
c      u(3)=phi
c      u(4)=cnz
c      u(5)=cnr
c      u(6)=cm

c      bf(1)=bz
c      bf(2)=br
c      bf(3)=bphi

c      gam=gamma1(z,r,phi,cnz,cnr,cm)
c      ds=dsin(gam)
c      dc=dcos(gam)
c      cnt=cn(r,cnz,cnr,cm)
c      cnpar=cnt*dc
c      cnper=cnt*ds
     
c      call rside1(0.d0,u,deru)

c      vgr(1)=deru(1)
c      vgr(2)=deru(2)
c      vgr(3)=r*deru(3)
c      vgrmods=vgr(1)**2+vgr(2)**2+vgr(3)**2
c      vgroup=dsqrt(vgrmods) !in clightc

c     calculate cnray
c      omegpedce=dsqrt(xe)/dabs(ye) !omega_pe/omega_ce
c      omegdce=1.d0/dabs(ye)        !omega/omega_ce
c      cnray=rrind(cnpar,cnper,omegdce,omegpedce)
c      cnn=dsqrt(cnpar**2+cnper**2)
c      write(*,*)'cnpar,cnper,cnt,cnn,cnray',
c     +           cnpar,cnper,cnt,cnn,cnray
c      wnray(n_midway)=cnray

c-----calculate electric field (cex,cey,cez) and flux cflown
c      call elf_emis(t,u,deru,cex,cey,cez,cflown)
      
c      call calc_emis_coef(xe,-ye,T_kev,cnpar,cnper,cnray,
c     +n_relt_harm1,n_relt_harm2,n_relt_intgr, 
c     +i_resonance_curve_integration_method,
c     +i_fkin,r,z,
c     +phi,
c     +cex,cey,cez,cflown,vgroup,frqncy,
c     +iabsorp_collisional_l,coll_mult_l,tempe_l,dense_l,zeff_l,
c     +wal_emis(n_midway),wj_emis(n_midway),wtemp_rad_em(n_midway))

c      write(*,*)'in midway wal_emis(n_midway),wj_emis(n_mindway)',
c     +n_midway,wal_emis(n_midway),wj_emis(n_midway)
c      write(*,*)'wtemp_em(n_midway),wtemp_rad_em(n_midway)',
c     +   wtemp_em(n_midway),wtemp_rad_em(n_midway)

cendtest

      endif ! i_emission=1

c-----------------------------------------------------------
      
cc      write(*,*)'prep3d is=',is
cc      write(*,*)'ws=',ws(is)
cc      write(*,*)'seikon=',seikon(is)
cc      write(*,*)'spsi=',spsi(is)
cc      write(*,*)'wr=',wr(is)
cc      write(*,*)'wphi=',wphi(is)
cc      write(*,*)'wz=',wz(is)
cc      write(*,*)'wnpar=',wnpar(is)
cc      write(*,*)'delpwr=',delpwr(is)
cc      write(*,*)'sdpwr=',sdpwr(is)
cc      write(*,*)'wdnpar=',wdnpar(is)
cc      write(*,*)'cwexde=',cwexde(is)
cc      write(*,*)'cweyde=',cweyde(is)
cc      write(*,*)'cwezde=',cwezde(is)
cc      write(*,*)'fluxn=',fluxn(is)
cc      write(*,*)'sbtot=',sbtot(is)
cc      write(*,*)'sene=',sene(is)
cc      write(*,*)'salphac=',salphac(is)
cc      write(*,*)'salphal=',salphal(is)
c-----------------------------------------------------------
c       data for onetwo
c-----------------------------------------------------------
c      write(*,*)'in prep3d before ionetwo'
c      write(*,*)'pre3d is,rhoold,rho,cnpar',is,rhoold,rho,cnpar
c      if (rho.gt.1.d0) then
c         rho=1.d0
c      endif

      if(ionetwo.eq.1) then
        if(is.gt.1) then       
           call check_monotonic_radius_in_ray_element(is,
     &     i_non_monotonic,z_center,r_center,rho_center)
        endif
c        write(*,*)'is,i_non_monotonic',is,i_non_monotonic
c        write(*,*)'prep3d before p_c_prof cex,cey,cez',
c     &                                    cex,cey,cez

cSAP_110503
        rhonew=rho
c        write(*,*)'prep3d before p_c_prof,rhoold,rhonew',rhoold,rhonew
        call p_c_prof(is,rhoold,rhonew,cnpar,cnper,cex,cey,cez)
c        call p_c_prof(is,rhoold,rho,cnpar,cnper,cex,cey,cez)
      endif

cSAP_110503		 		       
c      rhoold=rho
      rhoold=rhonew
c      write(*,*)'in prep3d after p_c_profm,rhoold',rhoold

c      call cpu_time(time_prep3d_2)
c      time_prep3d=time_prep3d+(time_prep3d_2-time_prep3d_1)
c      time_prep3d_emis=time_prep3d_emis+(time_prep3d_emis_2-
c     &                                    time_prep3d_emis_1)

c      write(*,*)'************* time from one call of prep3d ********'
c      write(*,*)'time_prep3d_2-time_prep3d_1',
c     &            time_prep3d_2-time_prep3d_1
c      write(*,*)'time_prep3d_emis_2-time_prep3d_emis_1',
c     &           time_prep3d_emis_2-time_prep3d_emis_1
c
 
c     write(*,*)'*** total time of all prep3d calls ******'

c      write(*,*)'time_prep3d,time_prep3d_emis',
c     &time_prep3d,time_prep3d_emis

c      write(*,*)'prep3d end'

      return
      END subroutine prep3d


	  double precision
     1    FUNCTION u_res(jwave,cnpar,temp,ye)
c resonanse velosity (for nonrelativistic case)
c----------------------------------------------------------------------
c      input parameters: jwave - wave harmonic number
c                        cnpar - refractive index parallel to 
c                                magnetic field
c                        temp -  temperature in kev
c			 ye -    (omega_Be/omega)
c      output parameter: u_res in thermal velosity v/ve, ve=sqrt(2Te/me)
c-----------------------------------------------------------------------

	  IMPLICIT double precision (a-h,o-z)

c				        ve in (m/sec)
	  ve=1.87d7*dsqrt(temp)
c                               c - light velocity in (m/sec)
	  c=3.0d8
	  u_res=c/ve*(1.0d0-dble(jwave)*ye)/cnpar
c	  write(*,*)'in u_res ye,cnpar,u_res',ye,cnpar,u_res
	  return
	  END function u_res


      double precision
     1FUNCTION efficien(z_eff,u1,jwave,temp,den)
c------------------------------------------------------------------
c     RF current drive efficiency(asymptotic formulas
c     ,nonrelativistic case)
c------------------------------------------------------------------
c     input parameters: z_eff-effective charge
c               	u1-resonance velocity  v/ve,ve=sqrt(2T_e/m_e)
c                       jwave -wave harmonic number
c                       temp - temperature in kev
c                       den - density in 10**19/m**3
c     output parameter:  efficiency in (A/cm**2)/(erg/(sec*cm**3))
c------------------------------------------------------------------
c YuP: this function is supposed to be called for ieffic=1 only,
c but can be called when ieffic>1, for a comparison/test purposes.

      IMPLICIT double precision (a-h,o-z)
cyup      write(*,*)'efficien: jwave=',jwave 
      ! jwave= ! can be (0 - LH wave, -1 AW, 1 - EC wave) 
      if(jwave.eq.1) then
c        EC wave first harmonic
         efficien=1.5d0/(5.d0+z_eff)*
     1   (u1*u1+(2.0d0+3.0d0/2.0d0/(3.d0+z_eff)))
      elseif(jwave.eq.0) then
c        LH wave
         efficien=2.d0/(5.d0+z_eff)*
     1   (u1*u1+(7.0d0/4.0d0+9.0d0/4.0d0/(3.d0+z_eff)))
      else ! YuP[07-2017] Added: other values of jwave
cyup         write(*,*)'WARNING: func.efficien() : wrong value of jwave !!!'
         !efficien=0.d0
         ! Or maybe set to be the same as for jwave=1?
         ! Here, same as for jwave=1:
         efficien=1.5d0/(5.d0+z_eff)*
     1   (u1*u1+(2.0d0+3.0d0/2.0d0/(3.d0+z_eff)))
      endif
cyup      write(*,*)'efficien: nondim asimptotic efficiency',efficien
c-----------------------------------------------------------
c     efficiency  in (A/cm**2)/(erg/(sec*cm**3))
c     temperature temp in kev
c     density     den in 10**19/m**3
c-----------------------------------------------------------
c     cln -Coulomb Ln
      cln=17.d0
      arg1=1.d3/temp*dsqrt(10.d0*den)	!temp KeV, den 10**13/cm**3
      cln=24.d0-dlog(arg1)
      efficien=efficien*temp/den*(17.d0/cln)*9.d0*1.d-6
      if (u1.ge.0.d0) then
         s=1.d0
      else
         s=-1.d0
      endif
cSm050329
c      efficien=efficien*s
      efficien=-efficien*s

c      write(*,*)'dim asimptotic efficien',efficien
      return
      END function efficien

c-----------------------------------------------------------------
c     this  subroutine calculates current drive efficiency
c     for the toroidal  plasma using CURBA code
c     ATTENTION:all parameters and variables input to curba
c       are real(no double precision)
c     Called only for ieffic=3, else might need to adjust lh=
c       calculation.
c-----------------------------------------------------------------
      subroutine effcurb(z,r,yma,xma,r0x,z_eff,temp,den,jwave,cnpar,ye,
     1                   efficient)
c     input parameters: z,r -coordinates of the ray point(cm)

c                       yma,zma -coordinates of magnetic axis(cm)
c                       r0x character length (m)
c                       z_eff-effective charge
c                       temp-temperature (keV)
c                       den-density(10**13 cm**-3)
c                       jwave -wave harmonic number
c                       cnpar -parallel to magnetic field refractive
c                              index
c                     	ye-omega_Be/omega_wave
c     output parameter: J/P=efficient  in (A/m**2)/(W/m**3)
c------------------------------------------------------------------
c     these double precision parameters are to call function: 
c        psif(z,r) and  subroutine: zr_psith(psi,theta,z,r)
c-------------------------------------------------------------------
      implicit none
c-----input
      real*8 z,r,yma,xma,z_eff,temp,den,cnpar,ye
      real*8 r0x
      integer jwave
c-----output
      real*8 efficient

c-----locals
      real*8 psid,zd,rd,rmaxpsid,rminpsid,pid,
     &zbmin,zbmax,zb
      real*8 tc,tol,elomom,prho,ctheta,stheta,pi,theta,thtc,enpar,denom,
     *aspct,zacosarg,rjpd,rjpd0,ratjpd,efficien
      integer ig,n0,model,lh
c-----external 
      real*8 psif,bmin_psi,bmax_psi,rmax_psi,rmin_psi,b,zfac_f

      
c-------------------------------------------------------------------
c     for curba_GA
      real*8 rjpd_d,rjpd0_d,ratjpd_d,denom_d,aspct_d,enpar_d, 
     &tc_d, thtc_d, theta_d, elomom_d,zeff_d,tol_d
c-----external
c-----for efficiency transformation like in currn in toray
      real*8 phi,zfac
c-------------------------------------------------------------------
cyup      write(*,*)'effcurb z,r,yma,xma,r0x,z_eff',
cyup     &z,r,yma,xma,r0x,z_eff
cSAP080617
cyup      write(*,*)'temp,den,jwave,cnpar,ye',temp,den,jwave,cnpar,ye

c     ig=+1 selects new	Green's fcn in curba.-1 selects old
c        +2 or +3 are older models
c---------------------------------------------------------------------
      ig=1
c--------------------------------------------------------------------
c     n0 is number of mesh points for gaussian integration of Green's
c        function; n0=2,4,6,8,10,12,16,20,24,32,48 or 64. n0=128 gives
c        64 points over resonanse if resonanse has only single passing
c        particle segment; 64 points on each if two segments.20,24,32,
c        48 or 64.
c---------------------------------------------------------------------
c      n0=4
       n0=32
c      n0=20

c--------------------------------------------------------------------
c    tc is a bulk electron temperature in keV
c--------------------------------------------------------------------
      tc=temp
      
c--------------------------------------------------------------------
c     tol is absolute toleranse for D02GBF integrator; starting in
c       version	1.1 the variables integrated by D02GBF are normalized
c       to be O(1), so tol can be viewed as relative error tolerense.
c--------------------------------------------------------------------
c      tol=2.e-3
      tol=1.d-2
c--------------------------------------------------------------------
c     model Absolute value of model selects collisional model: 1 for
c     full bounce av, 2 for square well(numerical solution), 3 for
c     analytic solution to square well.negative model does parallel
c     heating (lower hybrid, fast wave)
c-------------------------------------------------------------------
      model=3
c-------------------------------------------------------------------
c     z_eff is ion effective charge
c-------------------------------------------------------------------
c     lh ABS(lh) is power of e_perp in diffusion coeff; sign governs
c       power of p_parallel in diffusion coeff: + gives p_par^0;
c       - gives p_par^2. For ECRH, |lh|	is harmonic number (lh=1
c       for fundamental);+ for E_- contribution, - for E_parallel;
c        + with lh-->lh+2 for E_+ component of electric field E.
c       For now ,there is noprovision for the p_par^2 option with lh=0,
c       as the compiler doesn't know the difference between +0 and -0.
c       If there is any interest, sent a message to use 313 and
c       revision including this option will be created.
c--------------------------------------------------------------------

cBH151018 For ECCD (ieffic=3 or 4) choosing the harmonic for the
cBH151018 efficiency calculation to correspond to harmonic of
cBH151018 maximum EC power absorption.
cBH151018 Formula in accord with TORBEAM (as pointed out by Nicola
cBH151018 Bertelli).  This enables multiple harmonic CD, for example,
cBH151018 in DIII-D, 3rd harmonic near outer edge, 2nd near plasma
cBH151018 center.
cBH151018 Justification for the formula can be seen by considering,
cBH151018 e.g., multiharmonic case along with formulas 3.37 and 
cBH151018 Fig. 3.3 of http://www.compxco.com/cql3d_manual_150122.pdf .

      if (jwave.ne.0) then
         lh=jwave
      else ! jwave=0 (or -1 or 2?)
         !YuP lh=int(sqrt(1d0-cnpar**2)/abs(ye))+1
         !write(*,*)' effcurb: lh,ye,cnpar=',lh,ye,cnpar
         !YuP[2020-08-04] BUG: when cnpar>1, it results in 
         ! large negative values for lh, which gives INF 
         ! in other subroutines.  
         ! Adjusted to sqrt(max(0.d0, 1d0-cnpar**2))
         lh= idint( dsqrt(dmax1(0.d0,1d0-cnpar**2))/dabs(ye) ) +1
      endif
c     lh=2
c     lh=0

c--------------------------------------------------------------------
c     elomom is ABS (lh)*cyclotron frec over wave frec for Y<0,
c       interpret as evaluated at poloidal angle where electrons are
c       first trapped in bucket.
c--------------------------------------------------------------------
      elomom=lh*dabs(ye)
c--------------------------------------------------------------------
c     theta is a poloidal angle at which power is absorbed (measured
c       from outside); for thtc<0,theta is poloidal angle at which
c       electrons become trapped in bucket; used to calculate resonant
c       energy given elomom and enpar. Note if calling programm fixes
c       minimum resonant energy and calculates enpar, then theta has no
c       significance to the physics; rjpd depends only on
c       B_min*Y/B(theta) ,not on theta or Y alone. But it is useful to
c       be able to specify theta and Y separately in order to keep
c       track of what higher harmonic are doing.
c--------------------------------------------------------------------
      prho=dsqrt((z-yma)**2+(r-xma)**2) !YuP: it can get to 0.d0
      !write(*,*)' effcurb: prho=',prho
      if(prho.gt.0.d0)then  !YuP[2020-08-04] Added check of prho>0
        ctheta=(r-xma)/prho
        stheta=(z-yma)/prho
      else
        ctheta=1.d0
        stheta=0.d0
      endif
      pi=4.d0*datan(1.d0)
c      write(*,*)' effcurb ctheta_geom,stheta',ctheta,stheta
      if(stheta.ge.0.0d0) then
         theta=dacos(ctheta)
      else
         theta=2.0d0*pi-dacos(ctheta)
      end if
c      write(*,*)' effcurba acos(ctheta),theta',acos(ctheta),theta
c--------------------------------------------------------------------
c     thtc is ratio of temperature along characteristic to that of
c       bulk, except for thtc <0, -thtc is energy (in keV) of bucket
c       rise;
c--------------------------------------------------------------------
      thtc=1.0d0
c--------------------------------------------------------------------
c     enpar is k_parallel*c/wave frec. Note enpar**2 <1-elomom**2
c       implies no resonant particles; rjpd set to zero and rjpd0
c       determined for vparallel given by nonrelativistic resonance
c       condition
c--------------------------------------------------------------------
      enpar=cnpar
c--------------------------------------------------------------------
c     aspct is inverse aspct ratio
c--------------------------------------------------------------------
c bmax,bmin  in toray
c     if(dabs(bmax-bmin).lt.1.0d-8) bmax=bmin+1.0d-8
c 	aspct=(bmax-bmin)/(bmax+bmin)
c--------------------------------------------------------------------
c     Here :aspct is calculaed using rmax_psi and rmin_psi
c  conversion from real to double precision
      pid=4.d0*datan(1.d0)
c conversion from non-dimensional to m
c      write(*,*)'z,r',z,r
      zd=dble(z)*0.01d0/r0x
      rd=dble(r)*0.01d0/r0x
      psid=psif(zd,rd)
c      write(*,*)'zd,rd,psid',zd,rd,psid
cSm040426
cc     rmax_psid and rmin_psid are the largest and and the least
cc     values of the radius (r) on the given flux surface
cc     psid=psid(zd,rd)
c------------------------------
c      call zr_psith(psid,0.d0,zd,rmaxpsid)
c      call zr_psith(psid,pid,zd,rminpsid)
      
c      if(abs(rmaxpsi-rminpsi).lt.1.0e-8) rminpsi=rmaxpsi-1.0e-8
c      aspct=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
      
c----------------------------------
c     Here :aspct is calculating using double prcision functions
c           rmax_psi(psid) and rmin_psi(psid)
c----------------------------------------
c     rmaxpsi and rminpsi are the largest and and the least
c     values of the radius (rmaxpsi and rminpsi in (m))
c     on the given flux surface psi=psi(zd,rd)
c--------------------------------------------
c      rmaxpsid=rmax_psi(psid)
c      rminpsid=rmin_psi(psid)
c      if (rmaxpsid.lt.rd)rmaxpsid=rd
c      if (rminpsid.gt.rd)rminpsid=rd
c      if(dabs(rmaxpsid-rminpsid).lt.1.0d-8) rminpsid=rmaxpsid-1.0d-8
c      aspct=(rmaxpsid-rminpsid)/(rmaxpsid+rminpsid)
cSm040426
       zbmin=bmin_psi(psid)    
       zbmax=bmax_psi(psid)
       aspct=((zbmax-zbmin)/(zbmax+zbmin))
cSm060201
       if (aspct.lt.1.d-8) aspct=1.d-8 

c       write(*,*)'psid,zbmax,zbmin,aspct',psid,zbmax,zbmin,aspct
c      write(*,*)'aspct from functionsrmin_psi, rmax_psi',aspct
c      write(*,*)'aspct',aspct,'enpar',enpar,'tc',tc,'thtc',thtc
c      write(*,*)'theta',theta,'elomom',elomom,'lh',lh
c      write(*,*)'z_eff',z_eff,'model',model,'n0',n0,'ig',ig

cSm040503
c     function b() uses the arguments z and r in [m]
      zb=b(.01d0*dble(z)/r0x,.01d0*dble(r)/r0x,0.d0)!it changes 
                                         !bz,br,bphi and rho in one.i
      zacosarg=(zbmax + zbmin - 2.d0*zbmax*zbmin/zb)/
     &         (zbmin - zbmax)

c      write(*,*)'zbmax,zb,zbmin,zacosarg',zbmax,zb,zbmin,zacosarg
      if (dABS (zacosarg) .gt. 1.0d0) zacosarg = 
     1                               DSIGN (1.0d0, dble(zacosarg))
      theta =dACOS (zacosarg) 
      denom=1.d0
c      write(*,*)'effcurb befre curba denom,aspct,enpar,tc,thtc',
c     &denom,aspct,enpar,tc,thtc
c      write(*,*) 'theta,elomom,lh,z_eff,model,tol,n0,ig',
c     &            theta,elomom,lh,z_eff,model,tol,n0,ig

c-----using  TorGAcurba with real*8 arguments
      denom_d= denom
      aspct_d=aspct
      enpar_d=enpar 
      tc_d=tc
      thtc_d=thtc
      theta_d=theta
      elomom_d=elomom
      zeff_d=z_eff
      tol_d=tol
cSAP080617
      !write(*,*)'effcurb: before TorGA_curba lh=',lh

      call TorGA_curba(rjpd_d,rjpd0_d,ratjpd_d,denom_d,aspct_d,enpar_d,
     & tc_d, thtc_d, theta_d, elomom_d, lh, zeff_d, model, tol_d, n0,ig)
      rjpd=rjpd_d
      rjpd0=rjpd0_d
      ratjpd=ratjpd_d
      
!      write(*,*)'effcurb: aft TorGA_curba rjpd,rjpd0,ratjpd,denom_d,lh',
!     & rjpd,rjpd0,ratjpd,denom_d,lh

c----------------------------------------------------------
c     efficiency  in (A/cm**2)/(erg/(sec*cm**3))
c     temperature temp in kev
c     density     den in 10**19/m**3
c-----------------------------------------------------------
c-----genray transformation of CD efficiency 
c     cln -Coulomb Ln
c     it is necessary to write the formula to calculate cln

c      cln=17.0
c      arg1=1.d3/temp*dsqrt(10.d0*den) !temp KeV, den 10**13/cm**3
c      cln=24.d0-alog(arg1)
       efficien=rjpd 
c      efficient=efficien*temp/den*(17.0/cln)*4.5*1.e-6
c      write(*,*)'efficiecy coefficient=temp/den*(17.0/cln)*4.5*1.e-6',
c     1temp/den*(17.0/cln)*4.5*1.e-6
c      write(*,*)'in effic rjpd,rjpd0,ratjpd,denom',
c     1           rjpd,rjpd0,ratjpd,denom
c      write(*,*)'end effic curba efficient',efficient

c-----the efficiency transformation from torayGA
c      write(*,*)'before zfac_f temp',temp
      zfac=zfac_f(dble(z*1.d-2),dble(r*1.d-2),0.d0,dble(temp))
      zfac=zfac*1.d-7 ! watt => egr/sec
      zfac=zfac*1.d-4 ! m**2 => cm**2
c      write(*,*)'zfac',zfac
      efficient=efficien*zfac
cSm050923
      efficient=-efficient 

      return
      end subroutine effcurb

      double precision function zfac_f(z,r,phi,temp_kev)
c-----It calculates the coefficient
c     that transforms the CD efficiency from
c     curba variables to A/cm**2/erg/sec*cm**3))
c     It uses the formula from currn in toray.
c-----input  z(m),r(m),phi(radians)
c            temp_kev is the electron temperature in keV    
c     It uses function x() that uses the small radius rho.
c     rho is calculated inside function b() and b() puts rho to common/one/
c     So, wee need to call function b() before using zfac_f 


      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      double precision ld
c-----externals: x
      data zelchg/4.8032d-10/, zconst/10.83259d0/,zconvr/33.33333d0/
      data cvac/2.9979d10/
c       write(*,*)'in zfac_f z,r,phi,temp_kev', z,r,phi,temp_kev
      alfa=x(z,r,phi,1) !(omega_pe/omega)**2
      zralfa=dsqrt(alfa)

c     Subroutine to interface between TORAY and R. Cohen's calculation
c     of current.
c
c     zelchg is the electron charge in statcoulombs.
c     zconst is the ln of Boltzmann's constant (1.3807e-16) times the
c            factor to convert temperature from eV to degrees Kelvin
c            (1.1605e4) divided by the product of cvac times Planck's
c            constant (1.0546e-27).
c     cvac   is the light speed sm/sec
c     ld     is cvac divided by omega (2*pi*f).
c     zconvr is the constant that converts from (statamps/cm**2) divided
c            by (ergs/sec) to (amps/m**2) divided by watts.
      
      pi=4.d0*datan(1.d0)
      ld=cvac/(2.d0*pi*frqncy*1.d9) !c/omega cm
      zlds=ld*ld
      zfac1=zconst+dlog(ld)
      zte=temp_kev*1.d3 !eV       
      alfa=x(z,r,phi,1) !(omega_pe/omega)**2
      zralfa=dsqrt(alfa)
      zfac2=(zlds /(zelchg*alfa))*zte*2.d0/511.0d3
      zfac3=zfac1+ dlog (zte/zralfa)
      zfac_f=zconvr*zfac2/zfac3

      return
      end
      
      


      SUBROUTINE p_c_prof(is,rhobegin,rhoend,cnpar,cnper,
     &cefldx,cefldy,cefldz)
c----------------------------------------------------------------
c     this subroutine calculates absorted power profiles
c     calculates RF current  profiles	 (from deposited power)
c-----------------------------------------------------------------
      !implicit double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      INCLUDE 'three.i'
      INCLUDE 'onetwo.i'
      INCLUDE 'write.i'
      INCLUDE 'gr.i'
      INCLUDE 'rho.i'
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
     
c-----input
      real*8 rhobegin,rhoend,cnpar,cnper
      integer is
      complex*16 cefldx,cefldy,cefldz !polarization

c-----locals
      real*8 delpow_s,ppow_s
      dimension delpow_s(nbulka) !Added for indiv ion contrib[BH041009]
      dimension ppow_s(nbulka)  !Added for indiv ion contrib[BH041009]

      real*8 z_r,r_r,ymar,xmar,z_effr,tempr,denr,cnparr,yer,effic_r

      real*8 r_m,z_m,temp,ye,u1,den,z_eff,hro,rholeft,rhoright,
     &rhomax,rhomin,eff_rho_max,eff_rho_min,hrho,delpower,delpow_i,
     &delpow_e,delpow_cl,del,r0,z0,r0m,z0m,psiloc,zfacgeom,aspct,
     &delcurr,rho0,poloidlen,rho0_pol,delrho,
     &ppow,pcur,ppow_e,ppow_i,ppow_cl,
     &r_bmin_right,r_bmin_left,z_bmin_right,z_bmin_left,theta,psi

      integer i,kk,jbinmin,jbinmax,j,k
cSAP070831
      real*8 psi_loc,cos_theta_pol,sin_theta_pol,rho_loc,theta_pol
      integer n_theta_pol,ir
c-----external
      real*8  tempe,y,u_res,dense,zeff,psi_rho,b,efficien,rhov,rhos,
     &qsafety_psi,bmin_psi,bmax_psi,dvol_dpsi,rho_lrho,
     &b_average

      pi=4.d0*datan(1.d0)

c     wr and wz are in (cm) 
      r_m=wr(is)*0.01d0/r0x ! normalization
      z_m=wz(is)*0.01d0/r0x ! normalization

c       write(*,*)'p_c_prof NR,is,ieffic',NR,is,ieffic

c  radius rho is in common/one/,it was calculated in subroutine: b
c  rho is used in functions:tempe,dense,z_eff
c  The magnetic field bmode is in common/one/ .bmode is used
c  in function y(z,r,phi)


c begin if is=1      

c      write(*,*)'prep3d.f p_c_prof is,rhobegin,rhoend,rho,cnpar,cnper',
c     &                             is,rhobegin,rhoend,rho,cnpar,cnper

      if(is.eq.1) then

c         write(*,*)'prep3d.f p_c_prof is=1 rho=',rho
cSAP080731
        if (rho.ge.1.d0) then
c---------------------------------------------------------------
c          zero CD efficiency outside the plasma
c----------------------------------------------------------------
           eff(is)=0.d0
           goto 100
        endif
  


c        calculation of efficiency on the first step
         temp=tempe(z_m,r_m,wphi(is),1)
         ye=y(z_m,r_m,wphi(is),1)
         u1=u_res(jwave,cnpar,temp,ye)
         den=dense(z_m,r_m,wphi(is),1)
         z_eff=zeff(z_m,r_m,wphi(is))
         if (ieffic.eq.1) then
c -------------------------------------------------------------------
c     calculation of current drive efficiency using asymptotic formula
c--------------------------------------------------------------------
          eff(is)=efficien(z_eff,u1,jwave,temp,den) 

c         write(*,*)'prep3d asymptotic eff',eff(is)
         endif

         if (ieffic.eq.2) then
c -------------------------------------------------------------------
c     calculation of current drive efficiency using Ehst-Karney formula
c--------------------------------------------------------------------
c          write(*,*)'prep3d bef Karn z_m,r_m,r0x,z_eff,temp,den,jwave',
c     &    z_m,r_m,r0x,z_eff,temp,den,jwave

cfor adj
c          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)
c          cos_theta_pol=(r_m-xma)/rho_loc
c          sin_theta_pol=(z_m-yma)/rho_loc
c          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
c          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
c          if (sin_theta_pol.ge.0.d0) then
c             theta_pol=+dacos(cos_theta_pol)
c          else  
c             theta_pol=-dacos(cos_theta_pol)
c          endif

c          if (theta_pol.lt.0.d0) then
c             theta_pol=theta_pol+2.d0*pi
c          endif

c          if (theta_pol.gt.2.d0*pi) then
c             n_theta_pol=theta_pol/(2.d0*pi)
c             theta_pol=theta_pol-2.d0*pi*n_theta_pol
c          endif
         
  
c          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof ,cefldy,cefldz',cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
c          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
c     &    psi_loc,theta_pol

c          call CD_adj_LH_efficiency(cnpar,cnper,
c     &    psi_loc,theta_pol,cefldy,cefldz,
c     &    eff(is))
c          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
c-------------------------------------------------------------------
c           write(*,*)'prep3d.f  in p_c_prof before'//
c     &               'asimptotic_power_CD_density_LSC'
 
c           call asimptotic_power_CD_density_LSC(
c     &     z_m,r_m,0.d0,cnpar,cnper)
c     &     CD_dens,power_dens,CD_efficiency)
c           write(*,*)'prep3d.f  in p_c_prof after'//
c     &               'asimptotic_power_CD_density_LSC'

c           write(*,*)'prep3d z_eff',z_eff
           call efKarney(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
     +                 z_eff,temp,den,jwave,
     1                 cnpar,eff(is))
c           write(*,*)'prep3d Karney is,cnpar,eff(is)',is,cnpar,eff(is)
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coincided exactly with the efficiency obtained 
c          calculated bu subroutine call efKarney
c
c-------------------------------------------------------
c          call efKarney_Bonoli(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
c     +                 z_eff,temp,den,jwave,
c     1                 cnpar,eff(is))
c          write(*,*)'prep3d Karney_Bonoli is,cnpar,eff(is)',
c     &                is,cnpar,eff(is)
c------------------------------------------------------------
	 endif   ! ieffic.eq.2

         if (ieffic.eq.3) then
c -------------------------------------------------------------------
c          calculation of current drive efficiency using curba
c--------------------------------------------------------------------
           z_r=(wz(is))
           r_r=(wr(is))
           ymar=(yma*100.d0*r0x) !cm
           xmar=(xma*100.d0*r0x) !cm
           z_effr=(z_eff)
           tempr=(temp)
           denr=(den)
           cnparr=(cnpar)
           yer=(ye)

c          write(*,*)'in p_c_prof before effcurb is=1 z_r,r_r,ymar,xmar',
c     .    z_r,r_r,ymar,xmar
c          write(*,*)'z_effr,tempr,denr,cnparr,yer',
c     .    z_effr,tempr,denr,cnparr,yer
c          write(*,*)'in p_c_prof before effcurb'

      call effcurb(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,jwave,cnparr,
     +             yer,effic_r)
cyup           write(*,*)'in p_c_prof after effcurb efffic_r',effic_r
           eff(is)=(effic_r)
         endif     ! ieffic.eq.3

         if (ieffic.eq.4) then
c -------------------------------------------------------------------
c          calculation of current drive efficiency using Lin_liu
c          TorGA_curgap
c--------------------------------------------------------------------
           z_r=(wz(is))
           r_r=(wr(is))
           ymar=(yma*100.d0*r0x) !cm
           xmar=(xma*100.d0*r0x) !cm
           z_effr=(z_eff)
           tempr=(temp)
           denr=(den)
           cnparr=(cnpar)
           yer=(ye)

           !write(*,*)'in p_c_prof befor eff_Lin_Liu efffic_r',effic_r
           call eff_Lin_Liu(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,
     &                 jwave,
     &                 cnparr,
     &                 cnper,ioxm,ye,
     &                 cefldx,cefldy,cefldz,
     &                 effic_r)

           !write(*,*)'in p_c_prof after eff_Lin_Liu efffic_r',effic_r

           eff(is)=(effic_r)
         endif    ! ieffic.eq.4

         if (ieffic.eq.5) then
c -------------------------------------------------------------------
c        calculation of current drive efficiency using adj function 
c        for n_harmonic=0 case (LH wave)
c--------------------------------------------------------------------
          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)

          if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc=Inf
          cos_theta_pol=1.d0
          sin_theta_pol=0.d0
          else
          cos_theta_pol=(r_m-xma)/rho_loc
          sin_theta_pol=(z_m-yma)/rho_loc
          endif
          
          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
          if (sin_theta_pol.ge.0.d0) then
             theta_pol=+dacos(cos_theta_pol)
          else  
             theta_pol=-dacos(cos_theta_pol)
          endif

          if (theta_pol.lt.0.d0) then
             theta_pol=theta_pol+2.d0*pi
          endif

          if (theta_pol.gt.2.d0*pi) then
             n_theta_pol=theta_pol/(2.d0*pi)
             theta_pol=theta_pol-2.d0*pi*n_theta_pol
          endif
           
          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof ,cefldy,cefldz',cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
c          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
c     &    psi_loc,theta_pol

          call CD_adj_LH_efficiency(cnpar,cnper,
     &    psi_loc,theta_pol,cefldy,cefldz,
     &    eff(is))
cyup          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
         endif !ieffic.eq.5


         if (ieffic.eq.6) then
c -------------------------------------------------------------------
c        calculation of current drive efficiency using adj function 
c        for all harmoniucs general case
c--------------------------------------------------------------------
          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)

          if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc=Inf
          cos_theta_pol=1.d0
          sin_theta_pol=0.d0
          else
          cos_theta_pol=(r_m-xma)/rho_loc
          sin_theta_pol=(z_m-yma)/rho_loc
          endif

          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
          if (sin_theta_pol.ge.0.d0) then
             theta_pol=+dacos(cos_theta_pol)
          else  
             theta_pol=-dacos(cos_theta_pol)
          endif

          if (theta_pol.lt.0.d0) then
             theta_pol=theta_pol+2.d0*pi
          endif

          if (theta_pol.gt.2.d0*pi) then
             n_theta_pol=theta_pol/(2.d0*pi)
             theta_pol=theta_pol-2.d0*pi*n_theta_pol
          endif
           
          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof cefldx,cefldy,cefldz',
c     &                                    cefldx,cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
cyup          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
cyup     &    psi_loc,theta_pol

      
          call CD_adj_efficiency(cnpar,cnper,
     &    psi_loc,theta_pol,cefldx,cefldy,cefldz,
     &    eff(is))

cyup          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
         endif !ieffic.eq.6
         
          
c--------------------------------------------------------------------
c        toroidal and poloidal current drives efficiencies for is=1
 100     continue
         bmod=b(z_m,r_m,0.d0)        
c--------------------------------------------------------------------
         allpower=0.0d0
         allpw_e=0.0d0
         allpw_i=0.0d0
         allpw_cl=0.0d0

         allcur=0.0d0
        
c------- initialization arrays
c        for power and current
         do i=1,NR
	    power(i)=0.0d0
	    current(i)=0.0d0           
	    power_e(i)=0.0d0
	    power_i(i)=0.0d0
	    power_cl(i)=0.0d0
            cur_den_parallel(i)=0.d0 
	 enddo

c--------binvol(NR) calculations
         theta=0.d0 
         hro=1.d0/dble(NR-1)
         do i=1,NR-1 !for onetwo.i in [cm**3]
            rholeft=hro*(i-1)
            rhoright=rholeft+hro
            binvol(i)=voltot*(rhov(rhoright)**2-rhov(rholeft)**2)*1.d6
c            write(*,*)'in prep3d.f'
c            write(*,*)'i,rholeft,rhoright',rholeft,rhoright
c            write(*,*)'rhov(rholeft),rhov(rhoright)',
c     &                 rhov(rholeft),rhov(rhoright)
c            write(*,*)'NR,i,binvol(i)',NR,i,binvol(i)

           binarea(i)=areatot*(rhos(rhoright)**2-rhos(rholeft)**2)*1.d4

            psi=psi_rho(rholeft)
            call zr_psith(psi,theta,z_bmin_left,r_bmin_left)
            psi=psi_rho(rhoright)
            call zr_psith(psi,theta,z_bmin_left,r_bmin_right)

            binarea_pol(i)=pi*(r_bmin_right-r_bmin_left)*
     &                        (r_bmin_right+r_bmin_left)*1.d4

         enddo

c         stop 'prep3d.f  binvol'
 
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               do i=1,NR
                  power_s(i,kk)=0.0d0
               enddo
            enddo
         endif

	 goto 30
      endif
c end if is=1

c      write(*,*)' is,p_c_prof: rhobegin,rhoend',is,rhobegin,rhoend

      if(rhoend.gt.rhobegin) then
         rhomax=rhoend
         rhomin=rhobegin
      else
	 rhomin=rhoend
	 rhomax=rhobegin
      endif
cSAP080831
      if(rhomax.gt.1.d0) rhomax=1.d0
      if(rhomin.gt.1.d0) rhomin=1.d0-1.d-10    

      hrho=1.0d0/(NR-1)
      jbinmin=1
      jbinmax=NR-1

c      write(*,*)'NR,hrho,rhomin,rhomax',NR,hrho,rhomin,rhomax

      do j=1,NR-1
         if(rhomin.lt.(hrho*j)) then
           jbinmin=j
	    goto 10
	 endif
      enddo
10    continue
      do j=jbinmin,NR-1
         if(rhomax.lt.(hrho*j)) then
            jbinmax=j
	    goto 20
	 endif
      enddo
20    continue

c      write(*,*)'p_c_prof is,jbinmax&min',is,jbinmax,jbinmin

      if (jbinmin.eq.NR) jbinmin=NR-1

c      write(*,*)'1 p_c_prof is,jbinmax&min',is,jbinmax,jbinmin
     
c-----------------------------------------------------------
c     here delpower and allpower are in (erg/sec)
c------------------------------------------------------------
      delpower=delpwr(is-1)-delpwr(is)

cSmirnov970105 beg
cBH001017      delpow_e=0.5d0*(ws(is)-ws(is-1))*
cBH001017     1(delpwr(is-1)*sdpwr(is-1)+delpwr(is)*sdpwr(is))
cBH001017      delpow_i=0.5d0*(ws(is)-ws(is-1))*
cBH001017     1(delpwr(is-1)*salphal(is-1)+delpwr(is)*salphal(is))
cHarvey991017 beg
      delpow_i=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*sdpwr(is-1)+delpwr(is)*sdpwr(is))
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            delpow_s(kk)=0.5d0*(ws(is)-ws(is-1))*
     1        (delpwr(is-1)*salphas(is-1,kk)+delpwr(is)*salphas(is,kk))
         enddo
      endif
      
      delpow_e=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*salphal(is-1)+delpwr(is)*salphal(is))  

c      write(*,*)'p_c_prof is,delpwr(is-1),delpwr(is),delpower,delpow_e'
c     &,is,delpwr(is-1),delpwr(is),delpower,delpow_e

cHarvey991017 end
      delpow_cl=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*salphac(is-1)+delpwr(is)*salphac(is))

      del=delpow_i+delpow_e+delpow_cl
cHarvey970111 beg
cSmirnov050215      if (del.ne.0.d0) then
cSmirnov050215      Change 0.d0 to a small number
      if (del.gt.1.d-100) then
         delpow_e=delpow_e*delpower/del
         delpow_i=delpow_i*delpower/del
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               delpow_s(kk)=delpow_s(kk)*delpower/del
            enddo
         endif
         delpow_cl=delpow_cl*delpower/del

      endif
cSmirnov970105 end
      allpower=allpower+delpower
cSAP091021
c      write(*,*)'delpower,delpow_i,delpow_e',delpower,delpow_i,delpow_e
c      write(*,*)'delpow_s(k),k=2,nbulk',(delpow_s(k),k=2,nbulk)
c      write(*,*)'allpower,delpower',allpower,delpower
cSmirnov970106 beg
      allpw_e=allpw_e+delpow_e
      allpw_i=allpw_i+delpow_i
      allpw_cl=allpw_cl+delpow_cl
      
cSmirnov970106 end
c     write(*,*)' allpower,allpw_e,allpw_i,allpw_cl',
c    1 allpower,allpw_e,allpw_i,allpw_cl
c     write(*,*)'delpower,delpow_e,delpow_cl',
c    1delpower,delpow_e,delpow_cl
c-----------------------------------------------------------
c     r0 (cm)
      r0=0.5d0*(wr(is)+wr(is-1))

cSAP080831
c      write(*,*)'p_c_prof is>1  rho',rho
      if (rho.ge.1.d0) then
c---------------------------------------------------------------
c        zero CD efficiency outside the plasma
c----------------------------------------------------------------
         eff(is)=0.d0
         goto 110
      endif
c---- calculation of the efficiency 
      temp=tempe(z_m,r_m,wphi(is),1)
      ye=y(z_m,r_m,wphi(is),1)
      u1=u_res(jwave,cnpar,temp,ye)
      den=dense(z_m,r_m,wphi(is),1)
      z_eff=zeff(z_m,r_m,wphi(is))
   
      if (ieffic.eq.1) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using asymptotic formula
c--------------------------------------------------------------------
        eff(is)=efficien(z_eff,u1,jwave,temp,den)
c       write(*,*)'prep3d asymptotic eff',eff(is)
      endif

      if (ieffic.eq.2) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using Ehst-Karney formula
c--------------------------------------------------------------------
c      write(*,*)'prep3d bef Karn z_m,r_m,r0x,z_eff,temp,den,jwave',
c     &    z_m,r_m,r0x,z_eff,temp,den,jwave

cfor adj
c          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)
c          cos_theta_pol=(r_m-xma)/rho_loc
c          sin_theta_pol=(z_m-yma)/rho_loc
c          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
c          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
c          if (sin_theta_pol.ge.0.d0) then
c             theta_pol=+dacos(cos_theta_pol)
c          else  
c             theta_pol=-dacos(cos_theta_pol)
c          endif

c          if (theta_pol.lt.0.d0) then
c             theta_pol=theta_pol+2.d0*pi
c          endif

c          if (theta_pol.gt.2.d0*pi) then
c             n_theta_pol=theta_pol/(2.d0*pi)
c             theta_pol=theta_pol-2.d0*pi*n_theta_pol
c          endif
         
c          psi_loc=psi_rho(rho)
c          write(*,*)'prep3d.f in p_c_prof cefldy,cefldz',cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
c          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
c     &    psi_loc,theta_pol

c          call CD_adj_LH_efficiency(cnpar,cnper,
c     &    psi_loc,theta_pol,cefldy,cefldz,
c     &    eff(is))
c          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
c-------------------------------------------------------------------
c      write(*,*)'prep3d z_eff',z_eff   
c           write(*,*)'prep3d.f  in p_c_prof before'//
c     &               'asimptotic_power_CD_density_LSC'
 
c           call asimptotic_power_CD_density_LSC(
c     &     z_m,r_m,0.d0,cnpar,cnper)
c     &     CD_dens,power_dens,CD_efficiency)
c           write(*,*)'prep3d.f  in p_c_prof after'//
c     &               'asimptotic_power_CD_density_LSC'

      call efKarney(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
     & z_eff,temp,den,jwave,cnpar,
     &              eff(is))
c      write(*,*)'prep3d Karney is,cnpar,eff(is)',is,cnpar,eff(is)
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coinsided exectly wiht the efficiency obtained 
c          calculated bu subroutine call efKarney
c
c-------------------------------------------------------
c      call efKarney_Bonoli(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
c     & z_eff,temp,den,jwave,cnpar,
c     &              eff(is))
c      write(*,*)'prep3d Karney_Bonoli is,cnpar,eff(is)',is,cnpar,eff(is)

      endif  !ieffic.eq.2

      if (ieffic.eq.3) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using curba
c-------------------------------------------------------------------
        z_r=(wz(is))
        r_r=(wr(is))
        ymar=(yma*100.d0*r0x)
        xmar=(xma*100.d0*r0x)
        z_effr=(z_eff)
        tempr=(temp)
        denr=(den)
        cnparr=(cnpar)
        yer=real(ye)


c         write(*,*)'in p_c_prof before effcurbis.ne.1 z_r,r_r,
c     .   ymar,xmar,r0x,jwave',
c     .   z_r,r_r,ymar,xmar,r0x,jwave
c         write(*,*)'z_effr,tempr,denr,cnparr,yer',
c     .   z_effr,tempr,denr,cnparr,yer
c         write(*,*)'in p_c_prof before effcurb' 

ctest
        u1=u_res(jwave,cnparr,tempr,yer)
cyup        write(*,*)'jwave,cnpar,tempr,yer,u1',jwave,cnparr,tempr,yer,u1
        effic_r=efficien(z_effr,u1,jwave,tempr,denr) ! just for comparison
cyup        write(*,*)'asimptotic: z_effr,denr,effic_r',z_effr,denr,effic_r
cendtest
        call effcurb(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,jwave,
     +  cnparr,yer,effic_r)

cyup        write(*,*)'in p_c_prof after effcurb is,effic_r',is,effic_r
 
        eff(is)=effic_r
      endif !ieffic.eq.3

      if (ieffic.eq.4) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using eff_Lin_Liu
c-------------------------------------------------------------------
        z_r=(wz(is))
        r_r=(wr(is))
        ymar=(yma*100.d0*r0x)
        xmar=(xma*100.d0*r0x)
        z_effr=(z_eff)
        tempr=(temp)
        denr=(den)
        cnparr=(cnpar)
        yer=real(ye)


c         write(*,*)'in p_c_prof before effcurbis.ne.1 z_r,r_r,
c     .   ymar,xmar,r0x,jwave',
c     .   z_r,r_r,ymar,xmar,r0x,jwave
c         write(*,*)'z_effr,tempr,denr,cnparr,yer',
c     .   z_effr,tempr,denr,cnparr,yer
        !write(*,*)'p_c_prof before effcurb'

        call effcurb(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,jwave,
     +  cnparr,yer,effic_r)

        !write(*,*)'p_c_prof after effcurb is,effic_r',is,effic_r

        call eff_Lin_Liu(z_r,r_r,ymar,xmar,r0x,z_effr,tempr,denr,
     &                 jwave,
     &                 cnparr,
     &                 cnper,ioxm,ye,
     &                 cefldx,cefldy,cefldz,
     &                 effic_r)
        !write(*,*)'p_c_prof after eff_Lin_Liu effic_r',effic_r
   
        eff(is)=effic_r
      endif !4

      if (ieffic.eq.5) then
c -------------------------------------------------------------------
c        calculation of current drive efficiency using adj function 
c        for n_harmonic=0 case (LH wave)
c--------------------------------------------------------------------
          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)

          if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc=Inf
          cos_theta_pol=1.d0
          sin_theta_pol=0.d0
          else
          cos_theta_pol=(r_m-xma)/rho_loc
          sin_theta_pol=(z_m-yma)/rho_loc
          endif

          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
          if (sin_theta_pol.ge.0.d0) then
             theta_pol=+dacos(cos_theta_pol)
          else  
             theta_pol=-dacos(cos_theta_pol)
          endif

          if (theta_pol.lt.0.d0) then
             theta_pol=theta_pol+2.d0*pi
          endif

          if (theta_pol.gt.2.d0*pi) then
             n_theta_pol=theta_pol/(2.d0*pi)
             theta_pol=theta_pol-2.d0*pi*n_theta_pol
          endif
           
          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof ,cefldy,cefldz',cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
cyup          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
cyup     &    psi_loc,theta_pol

          call CD_adj_LH_efficiency(cnpar,cnper,
     &    psi_loc,theta_pol,cefldy,cefldz,
     &    eff(is))
cyup          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)


c         if (is.ge.65) call CD_adj_LH_efficiency_test(cnpar,cnper,
c     &      psi_loc,theta_pol,cefldy,cefldz,
c     &      eff(is))

      endif ! ieffic.eq.5


      if (ieffic.eq.6) then
c -------------------------------------------------------------------
c        calculation of current drive efficiency using adj function 
c        for n_harmonic=0 case (LH wave)
c--------------------------------------------------------------------
          rho_loc=dsqrt((r_m-xma)**2+(z_m-yma)**2)

          if(rho_loc.eq.0.d0)then ! YuP [06,2015] To avoid 1/rho_loc=Inf
          cos_theta_pol=1.d0
          sin_theta_pol=0.d0
          else
          cos_theta_pol=(r_m-xma)/rho_loc
          sin_theta_pol=(z_m-yma)/rho_loc
          endif

          if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
          if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
          if (sin_theta_pol.ge.0.d0) then
             theta_pol=+dacos(cos_theta_pol)
          else  
             theta_pol=-dacos(cos_theta_pol)
          endif

          if (theta_pol.lt.0.d0) then
             theta_pol=theta_pol+2.d0*pi
          endif

          if (theta_pol.gt.2.d0*pi) then
             n_theta_pol=theta_pol/(2.d0*pi)
             theta_pol=theta_pol-2.d0*pi*n_theta_pol
          endif
           
          psi_loc=psi_rho(rho)

c          write(*,*)'prep3d.f in p_c_prof ,cefldx,cefldy,cefldz',
c     &                                     cefldx,cefldy,cefldz
c          write(*,*)'prep3d.f in p_c_prof cnpar,cnper',cnpar,cnper
cyup          write(*,*)'prep3d.f in p_c_prof psi_loc,theta_pol',
cyup     &    psi_loc,theta_pol

ctest
        u1=u_res(jwave,cnpar,temp,ye)
cyup        write(*,*)'jwave,cnpar,temp,ye,u1',jwave,cnpar,temp,ye,u1
        effic_r=efficien(z_eff,u1,jwave,temp,den) ! just for comparison
cyup        write(*,*)'asimptotic: z_eff,den,effic_r',z_eff,den,effic_r
cendtest
          call CD_adj_efficiency(cnpar,cnper,
     &    psi_loc,theta_pol,cefldx,cefldy,cefldz,
     &    eff(is))

c          if ((is.ge.262).and.(is.le.290)) then 
c            call CD_adj_efficiency_test(cnpar,cnper,
c    &       psi_loc,theta_pol,cefldx,cefldy,cefldz,
c    &       eff(is))
c          endif

cyup          write(*,*)'prep3d ADJ is,cnpar,eff(is)',is,cnpar,eff(is)
      endif ! ieffic.eq.6
 110  continue

c      write(*,*)'p_c_prof after 110' 

c--------------------------------------------------------------------
c     delpower(erg/sec),delcurr(Ampere),r0(cm)
c-----------------------------------------------------------
cSmirnov970105 beg

c-----calculate parallel CD using the geometric factor 1/(2pi*r0) 
c      delcurr=delpow_e*0.5d0*(eff(is-1)+eff(is))/(2*pi*r0)
c      write(*,*)'1/(2*pi*r0)',1.d0/(2.d0*pi*r0)

c------geometric factor 1/(r*pi*R_mag_axis)
c      zfacgeom = 1.0d0 / (2.0d0 * pi * xma)/100.d0
c      write(*,*)'1/(2*pi*xma)/100.d0',1.d0/(2.d0*pi*xma)/100.d0

c-----calculate toroidal CD
      z0=0.5d0*(wz(is)+wz(is-1))
      r0m=r0*1.d-2
      z0m=z0*1.d-2      
      bmod=b(z0m,r0m,0.d0)
c      write(*,*)'z0m,r0m,bmod,rho',z0m,r0m,bmod,rho
      psiloc=psi_rho(rho)
c-----geometric factor ~ 1/b_averaged 
      zfacgeom = -2.d0 * pi * qsafety_psi(psiloc)/
     &     (dvol_dpsi(psiloc)*dpsimax*b_average(psiloc))/100.d0

c      write(*,*)'qsafety_psi(psiloc),b_average(psiloc)',
c     &           qsafety_psi(psiloc),b_average(psiloc)

c      write(*,*)'dvol_dpsi(psiloc),dpsimax',dvol_dpsi(psiloc),dpsimax
c      write(*,*)'zfacgeom',zfacgeom

c-----geometric factor ~ 1/b_ min/(1+aspct) old from Toray
c      zfacgeom = -2.d0 * pi * qsafety_psi(psiloc)/
c     &     (dvol_dpsi(psiloc)*dpsimax*bmin_psi(psiloc))/100.d0
c      aspct=(bmax_psi(psiloc)-bmin_psi(psiloc))/
c     &(bmax_psi(psiloc)+bmin_psi(psiloc))
c      zfacgeom = zfacgeom / (1.d0 + aspct)
c      write(*,*)'old zfacgeom',zfacgeom

      delcurr=delpow_e*0.5d0*(eff(is-1)+eff(is))*zfacgeom !toroidal current
                                                          !created by delpow_e
c      write(*,*)'prep3d zfacgeom', zfacgeom 

c      write(*,*)' prep3d p_c_prof delpow_e,eff(is-1),eff(is),delcurr',
c     +delpow_e,eff(is-1),eff(is),delcurr

cSAP080902 to plot delta power and delta current along the ray
      delpow_e_ar(is)=delpow_e
      delcur_par_ar(is)=delpow_e*0.5d0*(eff(is-1)+eff(is))
      
c-----toroidal and poloidal CD from old genray version      
      rho0=0.5*(spsi(is)+spsi(is-1)) ! the small radius
      rho0_pol=rho_lrho(rho0)
      poloidlen=rho0_pol*totlength*100.d0     ! poloidal length cm

c      write(*,*)'p_c_prof rho0_pol,totlength,r0,poloidlen',
c     +rho0_pol,totlength,r0,poloidlen
c      write(*,*)'p_c_prof r0,2*pi*r0,eff(is)',r0,2*pi*r0
c      write(*,*)'p_c_prof eff(is)',
c     *eff(is)
c      write(*,*)'bphi,bpol,bmod',bphi,dsqrt(bz**2+br**2),bmod
c      write(*,*)'p_c_prof rho0,rho0_pol,poloidlen(cm)',
c     * rho0,rho0_pol,poloidlen
            
cSmirnov970105 end
      allcur=allcur+delcurr                     !total toroidal current
      
c      write(*,*)'allcur',allcur
c      write(*,*)'delcurr',delcurr

c      write(*,*)'p_c_prof rhoend,rhobegin',rhoend,rhobegin
c      write(*,*)'is,rho0',is,rho0

      if(rhoend.gt.rhobegin) then
          eff_rho_max=eff(is)
          eff_rho_min=eff(is-1)
      else
          eff_rho_max=eff(is-1)
          eff_rho_min=eff(is)
      endif

c      write(*,*)'jbinmin,jbinmax',jbinmin,jbinmax     
       
      if (jbinmin.eq.jbinmax) then

c         write(*,*)'jbinmon=jbinmax,delpower,power(jbinmin)',
c     .              jbinmin,delpower,power(jbinmin)

         power(jbinmin)=power(jbinmin)+delpower

cSAP080731
c         write(*,*)'power(jbinmin)',power(jbinmin)
 
         cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
     &   delpow_e*0.5d0*(eff_rho_min+eff_rho_max)/binvol(jbinmin)

cSAP080731
c        write(*,*)'delpow_e,eff_rho_min,eff_rho_max,binvol(jbinmin)',
c     &             delpow_e,eff_rho_min,eff_rho_max,binvol(jbinmin)

c        write(*,*)'cur_den_parallel(jbinmin)',cur_den_parallel(jbinmin)

         current(jbinmin)=current(jbinmin)+delcurr   !toroidal current from bin

c         write(*,*)'cur_den_parallel(jbinmin)*binarea(jbinmin)',
c     &   cur_den_parallel(jbinmin)*binarea(jbinmin)

c         write(*,*)'current(jbinmin)',current(jbinmin)
c         write(*,*)'zfacgeom',zfacgeom
c         write(*,*)'0.5d0*(eff(is-1)+eff(is))',
c     &   0.5d0*(eff(is-1)+eff(is))
c         write(*,*)'is,eff(is-1),eff(is)',is,eff(is-1),eff(is)
c         write(*,*)'binarea(jbinmin)/binvol(jbinmin)',
c     &   binarea(jbinmin)/binvol(jbinmin)
c         write(*,*)'0.5d0*(eff_rho_min+eff_rho_max)',
c     &   0.5d0*(eff_rho_min+eff_rho_max)
c         write(*,*)'is,eff(is-1),eff(is)',is,eff(is-1),eff(is)
c         write(*,*)'eff_rho_min,eff_rho_max',eff_rho_min,eff_rho_max
      
         power_e(jbinmin)=power_e(jbinmin)+delpow_e
         power_i(jbinmin)=power_i(jbinmin)+delpow_i

c         write(*,*)'jbinmin,power_e(jbinmin),power_i(jbinmin)',
c     &              jbinmin, power_e(jbinmin),power_i(jbinmin)

         if (iabsorp.eq.3) then
            do kk=2,nbulk
               power_s(jbinmin,kk)=power_s(jbinmin,kk)+delpow_s(kk)
            enddo
         endif
       
         power_cl(jbinmin)=power_cl(jbinmin)+delpow_cl     

         goto 30
      endif

c      write(*,*)'prep3d  after goto30 delpower rhomax,rhomin',
c     . delpower,rhomax,rhomin    
       
c      delrho=dmax1((rhomax-rhomin),hrho) 
cSm040426
      delrho=rhomax-rhomin

      ppow=delpower/delrho 

      pcur=delcurr/delrho
      
cSmirnov970106 beg
      ppow_e=delpow_e/delrho
      ppow_i=delpow_i/delrho
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            ppow_s(kk)=delpow_s(kk)/delrho
         enddo
      endif
      
      ppow_cl=delpow_cl/delrho
cSmirnov970106 end

c------------------------------------------------------------
c     power (erg/sec), current(A)
c------------------------------------------------------------

c      write(*,*)'jbinmin,power(jbinmin),ppow,(hrho*jbinmin-rhomin)',
c     &jbinmin,power(jbinmin),ppow,(hrho*jbinmin-rhomin)

      power(jbinmin)=power(jbinmin)+ppow*(hrho*jbinmin-rhomin)     

c      write(*,*)'power(jbinmin)',power(jbinmin)

cSmirnov970106 beg
      power_e(jbinmin)=power_e(jbinmin)+ppow_e*(hrho*jbinmin-rhomin)
      power_i(jbinmin)=power_i(jbinmin)+ppow_i*(hrho*jbinmin-rhomin)
      power_cl(jbinmin)=power_cl(jbinmin)+ppow_cl*(hrho*jbinmin-rhomin)
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            power_s(jbinmin,kk)=
     &           power_s(jbinmin,kk)+ppow_s(kk)*(hrho*jbinmin-rhomin)
         enddo
      endif

cSmirnov970106 end
c      write(*,*)'p_c_prof delpow_e',delpow_e

      cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
     &(delpow_e/delrho)*(hrho*jbinmin-rhomin)/binvol(jbinmin)*
     &0.5d0*(eff_rho_min+
     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
     
     &                   (hrho*jbinmin-rhomin)/(delrho)))

      current(jbinmin)=
     1 current(jbinmin)+pcur*(hrho*jbinmin-rhomin)  !toroidal current from bin

c     write(*,*)'p_c_prof jbinmin'
c     write(*,*)'current(jbinmin)',current(jbinmin)
c     write(*,*)'cur_den_parallel(jbinmin)*binarea(jbinmin)',
c    &cur_den_parallel(jbinmin)*binarea(jbinmin)

   
c      write(*,*)'jbinmax,power(jbinmax),ppow,(hrho*jbinmax-rhomin)',
c     &jbinmax,power(jbinmax),ppow,(hrho*jbinmax-rhomin)
 
      power(jbinmax)=power(jbinmax)+ppow*(rhomax-hrho*(jbinmax-1))

c      write(*,*)' power(jbinmax)', power(jbinmax)
cSmirnov970106 beg
      power_e(jbinmax)=power_e(jbinmax)+ppow_e*(rhomax-hrho*(jbinmax-1))
      power_i(jbinmax)=power_i(jbinmax)+ppow_i*(rhomax-hrho*(jbinmax-1))
      power_cl(jbinmax)=power_cl(jbinmax)+
     1                  ppow_cl*(rhomax-hrho*(jbinmax-1))
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            power_s(jbinmax,kk)=
     &          power_s(jbinmax,kk)+ppow_s(kk)*(rhomax-hrho*(jbinmax-1))
         enddo
      endif
cSmirnov970106 end

      cur_den_parallel(jbinmax)=cur_den_parallel(jbinmax)+
     &(delpow_e/delrho)*(rhomax-hrho*(jbinmax-1))/binvol(jbinmax)*
     &0.5d0*(eff_rho_max+
     &   (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &   (rhomax-hrho*(jbinmax-1))/(delrho)))
     
      current(jbinmax)=
     1	  current(jbinmax)+pcur*(rhomax-hrho*(jbinmax-1)) !toroidal current from bin

c     write(*,*)'p_c_prof jbinmax'
c     write(*,*)'current(jbinmax)',current(jbinmax)
c      write(*,*)'cur_den_parallel(jbinmax)*binarea(jbinmax)',
c     &cur_den_parallel(jbinmax)*binarea(jbinmax)

c      write(*,*)'jbinmin,jbinmax',jbinmin,jbinmax
c      write(*,*)'power(jbinmin),power(jbinmax)',
c     &           power(jbinmin),power(jbinmax)

c      write(*,*)'power_e(jbinmin)',
c     1           power_e(jbinmin)
c      write(*,*)'power_e(jbinmax)',
c     1           power_e(jbinmax)


      if(jbinmax.gt.(jbinmin+1)) then
         do j=(jbinmin+1),(jbinmax-1)
c            write(*,*)'j,power(j),ppow,hrho',j,power(j),ppow,hrho
            power(j)=power(j)+ppow*hrho
c            write(*,*)'power(j)',power(j)
cSmirnov970106 beg
            power_e(j)=power_e(j)+ppow_e*hrho
            power_i(j)=power_i(j)+ppow_i*hrho
            power_cl(j)=power_cl(j)+ppow_cl*hrho
            if (iabsorp.eq.3) then
               do kk=2,nbulk
                  power_s(j,kk)=power_s(j,kk)+ppow_s(kk)*hrho
               enddo
            endif

cSmirnov970106 end
            cur_den_parallel(j)=cur_den_parallel(j)+
     &      (delpow_e/delrho)*hrho/binvol(j)*
     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &      (hrho*(j-0.5d0)-rhomin)/delrho)
     
            current(j)=current(j)+pcur*hrho             !toroidal current from bins
         
c            write(*,*)'j,current(j),cur_den_parallel(j)*binarea(j)',
c     &      j,current(j),cur_den_parallel(j)*binarea(j)

c            write(*,*)'j,power(j)',j,power(j)

c            write(*,*)'j,power_e(j)',
c     .                 j,power_e(j)

         enddo
      endif

30    continue
      
      
c*******test
c       if (is.gt.1)then
c          write(*,*)'is,jbinmin,jbinmax.is',is,jbinmin,jbinmax
c          do ir=jbinmin,jbinmax
c            write(*,*)'ir,power(ir)',ir,power(ir)
c          enddo
c       endif
        
c      write(*,*)'allpower',allpower
c       write(*,*)'p_c_prof allcur', allcur
c      do ir=1,NR
c        write(*,*)'ir,power(ir),current(ir),temparr(ir),
c     1   spower(ir),scurrent(ir)',
c     2   ir,power(ir),current(ir),temparr(ir),
c     1   spower(ir),scurrent(ir)
c      enddo
c       do ir=1,NR
c        write(*,*)'ir,power(ir)',ir,power(ir)
c       enddo
c       spower_tot=0.d0
c       power_tot=0.d0
c       do ir=1,NR-1
c         spower_tot=spower_tot+spower(ir)
c         power_tot=power_tot+power(ir)
c       enddo
c       write(*,*)'delpwr(1)-delpwr(is)',delpwr(1)-delpwr(is)
c       write(*,*)'allpower,spower_tot,power_tot',
c     &            allpower,spower_tot,power_tot
c       write(*,*)'allpower,power_tot',
c     &            allpower,power_tot
c       write(*,*)'in p_c_prof before 99 return'
c       do ir=1,NR-1
c       write(*,*)'ir,cur_den_parallel(ir)',ir,cur_den_parallel(ir)
c       enddo
     
c99    return
      END SUBROUTINE p_c_prof


c-----------------------------------------------------------------------
c     this subroutine SONETWO is called after each ray finished
c     it calculates power (spower(i)) and current
c     (scurrent(i)) profiles
c     as sum the same profiles for all rays.
c-----------------------------------------------------------------------
c     input parameter:iray -number of the ray is in common/cone/
c-----------------------------------------------------------------------
      SUBROUTINE sonetwo
      IMPLICIT double precision (a-h,o-z)
      include 'param.i'
      include 'onetwo.i'
      INCLUDE 'cone.i'
c for test
      INCLUDE 'write.i'
      INCLUDE 'one.i'     !BH041020, for iabsorp,nbulk

c      write(*,*)'****************in sonetwo***************'
c      do i=1,NR
c        write(*,*)'i,power(i),spower(i),current(i),scurrent(i)',
c     1	 i,power(i),spower(i),current(i),scurrent(i)
c      end do
c-----------------------------------------------------------------------
c     spower (erg/sec), scurrent(A)
c-----------------------------------------------------------------------
      do i=1,NR-1
         spower(i)=spower(i)+power(i)
cSmirnov970106 beg
         spower_e(i)=spower_e(i)+power_e(i)
         spower_i(i)=spower_i(i)+power_i(i)
         spower_cl(i)=spower_cl(i)+power_cl(i)
cSmirnov970106 end
         scurrent(i)=scurrent(i)+current(i)                !toroidal current

c         write(*,*)' sonetwo i,cur_den_parallel(i)',
c     &                       i,cur_den_parallel(i)
     
         s_cur_den_parallel(i)=s_cur_den_parallel(i)+cur_den_parallel(i)

c         write(*,*)' sonetwo i,s_cur_den_parallel(i)',
c     &                       i,s_cur_den_parallel(i)
      enddo
      if (iabsorp.eq.3) then
         do kk=1,nbulk
            do i=1,NR-1
               spower_s(i,kk)=spower_s(i,kk)+power_s(i,kk)
            enddo
         enddo
      endif

      do i=1,NR
c       write(*,*)'i,power_e(i),spower_e(i)',
c     1   i,power_e(i),spower_e(i)
      enddo
c----------------------------------------------------------------------
c   test
c   absorbed power on the ray abspwer
      abspwer=0.d0
cSmirnov970106 beg
      abspw_e=0.d0
      abspw_i=0.d0
      abspw_cl=0.d0
cSmirnov970106 end
      curntray=0.d0
      curtrray=0.d0     !toroidal current
      curplray=0.d0     !poloidal current
      do i=1,NR
         abspwer=abspwer+power(i)
cSmirnov970106 beg
         abspw_e=abspw_e+power_e(i)
         abspw_i=abspw_i+power_i(i)
         abspw_cl=abspw_cl+power_cl(i)
cSmirnov970106 end
         curntray=curntray+current(i)    !toroidal current from one ray        
      end do
cSm040426
c      write(*,*)'abspwer',abspwer,'curentray',curntray
c      write(*,*)'abspw_e,abspw_i,abspw_cl',abspw_e,abspw_i,abspw_cl
c      absdpwr=0.0d0
c     write(*,*)'in sonetwo nrayelt',nrayelt
c      do i=1,nrayelt
c        write(*,*)'i,delpwr(i)',
c     1	 i,delpwr(i)
c      end do
cSm040426
c      absdpwr=delpwr(1)-delpwr(nrayelt)
c      write(*,*)'absdpwr',absdpwr
      return
      end


c-----------------------------------------------------------------------
c     this subroutine DNONETWO calculates power and current density
c     profiles on rho. Arrays: powden(NR) (erg/(cm**3*sec)
c                           and currden(NR).
c-----------------------------------------------------------------------
      SUBROUTINE dnonetwo
      !IMPLICIT double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'onetwo.i'
      include 'gr.i'
      include 'rho.i'
      include 'three.i'
      include 'five.i'
      include 'one.i'
      include 'lsc_approach.i'

c-----externals
      real*8 rhov,rhos,rho_lrho,psi_rho,ias1r,b

c-----locals
      real*8 pwtot_s
      dimension pwtot_s(nbulka)  !Added for indiv ion contrib[BH041009]
      real*8 hro,rholeft,rhoright,binplarea,h_rho,psi,f_eqd,
     &b_av,bs_av,r_av,dr_av,drs_av,theta,r_bmin,z_bmin,b_pol_bmin,
     & powertot,rho0,rho0_pol,poloidlen,
     &powrtot1,currtot1,pwtot_e,pwtot_i,pwtot_cl,
     &rho_l,psi_l,
     &total_power  
      integer i,kk,k,idx,nx4,j
     
      pi=4.0d0*datan(1.0d0)
c-----------------------------------------------------------------------
cSAP080731
cyup      write(*,*)'in dnonetwo'
c     write(*,*)'spower'
c     write(*,*)(spower(i), i=1,NR-1)
c     write(*,*)'spower_e'
c     write(*,*)(spower_e(i), i=1,NR-1)
c     write(*,*)'spower_i'
c     write(*,*)(spower_i(i), i=1,NR-1)
c     write(*,*)'spower_cl'
c     write(*,*)(spower_cl(i), i=1,NR-1)
c      write(*,*)'spower_s'
      
c      do i=1,NR-1
c        do kk=1,nbulk
c         write(*,*)'kk,i,spower_s(i,kk)',kk,i,spower_s(i,kk)
c        enddo
c      enddo

c      write(*,*)((spower_s(i,kk), i=1,NR-1),kk=1,nbulk)
c      write(*,*)'scurrent'
c      write(*,*)(scurrent(i), i=1,NR-1)      
c     write(*,*)'volume and area total:',voltot,areatot

c----------------------------------------------------------------------
c     spower (erg/sec), scurrent(A), binvol(cm**3),binarea(cm**2)
c     powden(erg/(sec*cm**3)),currden(A/cm**2)
c     voltot (m**3), areatot (m**2)
c----------------------------------------------------------------------
      hro=1.d0/dble(NR-1)
      do i=1,NR-1
         rholeft=hro*(i-1)
         rhoright=rholeft+hro

c        write(*,*)'RR,i,rholeft,rhoright',NR,i,rholeft,rhoright

         binvol(i)=voltot*(rhov(rhoright)**2-rhov(rholeft)**2)
     1          *1000000.d0
c         write(*,*)'rhov(rhoright)',rhov(rhoright)
c         write(*,*)'rhov(rholeft)',rhov(rholeft)
c         write(*,*)'voltot(m**3),binvol(cm**3',voltot,binvol(i)

         binarea(i)=areatot*(rhos(rhoright)**2-rhos(rholeft)**2)
     1            *10000.d0
c         write(*,*)'rhos(rhoright)',rhos(rhoright)
c         write(*,*)'rhos(rholeft)',rhos(rholeft)
c         write(*,*)'areatot(m**2),binarea(cm**2)',areatot,binarea(i)

         pollen(i)=
     1        totlength*50.d0*(rho_lrho(rhoright)+rho_lrho(rholeft))
c         write(*,*)'rho_lrho(rhoright)',rho_lrho(rhoright)
c         write(*,*)'prep3d totlength(m),pollen(cm)',
c     1   totlength,,pollen(i)
         binplarea=binvol(i)/pollen(i)
c	 write(*,*)'i',i, 'binvol!!(cm**3)',binvol(i)
c         write(*,*)'binarea cm**2',binarea(i)
c         write(*,*)'poloidal pollen (cm)',pollen(i)
c         write(*,*)'binplarea cm**2',binplarea

	 powden(i)=spower(i)/binvol(i)
cSmirnov970106 beg
	 powden_e(i)=spower_e(i)/binvol(i)
	 powden_i(i)=spower_i(i)/binvol(i)


         if(iabsorp.eq.3) then
            do kk=2,nbulk
               powden_s(i,kk)=spower_s(i,kk)/binvol(i)
            enddo
         endif

	 powden_cl(i)=spower_cl(i)/binvol(i)

cSmirnov970106 end
	 currden(i)=scurrent(i)/binarea(i)    !toroidal current density      
      enddo  ! i=1,NR-1

c      do kk=2, nbulk
c         write(*,*)'prep3d: first powden_s(i,kk)=', 
c     1                   (powden_s(i,kk),i=1,NR-1)
c      enddo
         
c       write(*,*)'spower'
c       write(*,*)(spower(i), i=1,NR-1)
c       write(*,*)'scurrent'
c       write(*,*)(scurrent(i), i=1,NR-1)

c$$$      do i=2,NR-1
c$$$cSm040426
c$$$         powden_s(i)=0.5d0*(powden(i-1)+powden(i))
c$$$
c$$$         powden_e_s(i)=0.5d0*(powden_e(i-1)+powden_e(i))
c$$$         powden_i_s(i)=0.5d0*(powden_i(i-1)+powden_i(i))
c$$$         powden_cl_s(i)=0.5d0*(powden_cl(i-1)+powden_cl(i))
c$$$
c$$$         currden_s(i)=0.5d0*(currden(i-1)+currden(i))
c$$$      end do
c$$$
c$$$      powden_s(1)=powden(1)
c$$$      powden_e_s(1)=powden_e(1)
c$$$      powden_i_s(1)=powden_i(1)
c$$$      powden_cl_s(1)=powden_cl(1)
c$$$      currden_s(1)=currden(1)
c$$$
c$$$      powden_s(NR)=powden(NR-1)
c$$$      powden_e_s(NR)=powden_e(NR-1)
c$$$      powden_i_s(NR)=powden_i(NR-1)
c$$$      powden_cl_s(NR)=powden_cl(NR-1)
c$$$      currden_s(NR)=currden(NR-1)

c      write(*,*)' powden_e',(powden_e(i), i=1,NR)

c      do k=2,nbulk
c        write(*,*)'k=',k,'powden_s'
c        write(*,*)(powden_s(i,k), i=1,NR)
c      enddo

c      write(*,*)'curdens'
c      write(*,*)(currden_s(i), i=1,NR)

c-----------------------------------------------------------
c     CD calculation using GA memo
c-----------------------------------------------------------
c     cur_den_onetwo=<j_parallel.B>/B_0=<j_parallel><B**2>/<B>/B_0
c----------------------------------------------------------
      h_rho=1.d0/(NR-1)
      idx=0 
      nx4=nx+4
 
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,1010)'cur_den_tor=p_tor*cur_den_par, '
     &   // 'p_tor = drs_av*f_eqd/(b_av*dr_av)'
      write(*,1010)'cur_den_pol=p_pol*cur_den_par, '
     &   // 'p_pol = b_pol_bmin/b_av'

      write(*,1010)'i rho_bin_center powden_e   cur_den_par  p_tor'
     & // '      cur_den_tor'
     & // '   p_pol    cur_den_pol'
      endif ! outprint
    
 1010    format(/,1x,a)
 1011    format(i3,7(1pe12.4))

      do i=1,NR-1
        rho_l=h_rho*(i-0.5d0)    
        psi_l=psi_rho(rho_l)
c        write(*,*)'i,rho_l,psi_l',i,rho_l,psi_l
cSAP080321 argument psi in f_eqd was changed from psi to psi_l
        f_eqd=ias1r(txf,nx,nx4,cx,idx,psi_l)

        call average_variables(psi_l,b_av,bs_av,r_av,dr_av,drs_av)

c        write(*,*)'after average_variables psi_l',psi_l
        
c        write(*,*)'b_av,bs_av,r_av,dr_av,drs_av',
c     &             b_av,bs_av,r_av,dr_av,drs_av

c        write(*,*)'s_cur_den_parallel(i)',s_cur_den_parallel(i)

        s_cur_den_onetwo(i)   = s_cur_den_parallel(i)*bs_av/(b_av*beqd)       
        s_cur_den_toroidal(i) = s_cur_den_parallel(i)*drs_av*f_eqd/
     &                          (b_av*dr_av)

c        write(*,*)'i,rho_l,s_cur_den_parallel(i)',
c     &  i,rho_l,s_cur_den_parallel(i)
c        write(*,*)'drs_av,f_eqd,b_av,dr_av,s_cur_den_toroidal(i)', 
c     &             drs_av,f_eqd,b_av,dr_av,s_cur_den_toroidal(i)

        theta=0.d0

c        write(*,*)'before zr_psith psi_l,theta ',psi_l,theta
 
        call zr_psith(psi_l,theta,z_bmin,r_bmin)

c        write(*,*)'psi_l,theta,z_bmin,r_bmin ',
c     &  psi_l,theta,z_bmin,r_bmin

        bmod=b(z_bmin,r_bmin,0.d0)
        b_pol_bmin=dsqrt(bz**2+br**2) ! poloidal B at the point
                                      ! with minimal B
                                      ! with theta poloidal =0

c        write(*,*)'bz,br,b_pol_bmin ',bz,br,b_pol_bmin 

        s_cur_den_poloidal(i) = s_cur_den_parallel(i)*b_pol_bmin/
     &                          b_av
c        write(*,*)'s_cur_den_onetwo(i)', s_cur_den_onetwo(i)
c        write(*,*)'s_cur_den_toroidal(i)',s_cur_den_toroidal(i)
c        write(*,*)'s_cur_den_poloidal(i)',s_cur_den_poloidal(i) 

c        write(*,*)'b_pol_bmin,b_av,b_pol_bmin/b_av ',
c     &             b_pol_bmin,b_av,b_pol_bmin/b_av
c        write(*,*)'drs_av,f_eqd,b_av,dr_av ',drs_av,f_eqd,b_av,dr_av 
c        write(*,*)'drs_av*f_eqd/(b_av*dr_av)',drs_av*f_eqd/(b_av*dr_av)
c----------------------------------------------------------------
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,1011)i,rho_l,powden_e(i),s_cur_den_parallel(i),
     &   drs_av*f_eqd/(b_av*dr_av),s_cur_den_toroidal(i),
     &   b_pol_bmin/b_av,s_cur_den_poloidal(i)
         endif ! outprint
      enddo 
c-----------------------------------------------------------
c     powertot (erg/sec), currtot(A)
c-----------------------------------------------------------
     
cSm040426begin
      powertot=0.0d0
      powtot_e=0.0
      powtot_i=0.0d0
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            powtot_s(kk)=0.0d0
         enddo
      endif
      powtot_cl=0.0d0
      currtot=0.0d0
   
      do i=1,NR-1
c     integration formulas ****INT1**
         powertot=powertot+powden(i)*
     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2	          voltot*1000000.0d0
         powtot_e=powtot_e+powden_e(i)*
     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2	          voltot*1000000.0d0
         powtot_i=powtot_i+powden_i(i)*
     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2	          voltot*1000000.0d0
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               powtot_s(kk)=powtot_s(kk)+powden_s(i,kk)*
     1              (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2              voltot*1000000.0d0
            enddo
         endif
         powtot_cl=powtot_cl+powden_cl(i)*
     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2	          voltot*1000000.0d0

	 currtot=currtot+currden(i)*
cSm030428     1	         (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     1	         (rhos(hro*i)**2-rhos(hro*(i-1))**2)*
     2	          areatot*10000.0d0

c--------poloidal total current
         rho0=hro*(i+0.5d0)                 !  small radius
         rho0_pol=rho_lrho(rho0)
c        poloidlen=rho0_pol*totlength*100.d0*r0x! poloidal length cm
cSm040728         poloidlen=rho0_pol*totlength*100.d0 ! poloidal length cm
         poloidlen=totlength*50.d0*(rho_lrho(hro*i)+rho_lrho(hro*(i-1)))
	 binvol(i)=voltot*(rhov(hro*i)**2-rhov(hro*(i-1))**2)
     1          *1000000.d0
          	          
      enddo
     
      parallel_cur_total=0.d0
      toroidal_cur_total=0.d0
      poloidal_cur_total=0.d0
      do j=1,NR-1
        parallel_cur_total=parallel_cur_total+
     &                     s_cur_den_parallel(j)*binarea(j)
        toroidal_cur_total=toroidal_cur_total+
     &                     s_cur_den_toroidal(j)*binarea(j)
        poloidal_cur_total=poloidal_cur_total+
     &                     s_cur_den_poloidal(j)*binarea_pol(j)
      enddo

      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'parallel_cur_total, toroidal_cur_total ',
     &parallel_cur_total, toroidal_cur_total 
      write(*,*)'poloidal_cur_total ',poloidal_cur_total

      write(*,*)'INT1 testing DNONETWO. powertot=erg/sec',powertot,
     1' currtot=A',currtot
 
      write(*,*)'testing 1 DNONETWO powtot_e,powtot_i,powtot_cl',
     &     powtot_e,powtot_i,powtot_cl
      write(*,*)'powtot_s(1:nbulk)=',(powtot_s(kk),kk=1,nbulk)
      endif ! outprint
cSm040426end

c$$$      powertot=0.0d0
c$$$cSmirnov970106 beg
c$$$      powtot_e=0.0
c$$$      powtot_i=0.0d0
c$$$      if (iabsorp.eq.3) then
c$$$      do kk=2,nbulk
c$$$         powtot_s(kk)=0.0d0
c$$$      enddo
c$$$      endif
c$$$      powtot_cl=0.0d0
c$$$cSmirnov970106 end
c$$$      currtot=0.0d0
c$$$cSm040426begin
c$$$c      do i=1,NR-1
c$$$c     integration formulas ****INT**
c$$$c         powertot=powertot+(powden(i)+powden(i+1))*
c$$$c     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)/4.0d0*
c$$$c     2	          2.0d0*voltot*1000000.0d0
c$$$cSmirnov970106 beg
c$$$c         powtot_e=powtot_e+(powden_e(i)+powden_e(i+1))*
c$$$c     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)/4.0d0*
c$$$c     2	          2.0d0*voltot*1000000.0d0
c$$$c         powtot_i=powtot_i+(powden_i(i)+powden_i(i+1))*
c$$$c     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)/4.0d0*
c$$$c     2	          2.0d0*voltot*1000000.0d0
c$$$c         powtot_cl=powtot_cl+(powden_cl(i)+powden_cl(i+1))*
c$$$c     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)/4.0d0*
c$$$c     2	          2.0d0*voltot*1000000.0d0
c$$$cSmirnov970106 end
c$$$c	 currtot=currtot+(currden(i)+currden(i+1))*
c$$$c     1	         (rhov(hro*i)**2-rhov(hro*(i-1))**2)/2.0d0*
c$$$c     2	          areatot*10000.0d0

c$$$c--------poloidal total current
c$$$c         rho0=hro*(i+0.5d0)                 !  small radius
c$$$c         rho0_pol=rho_lrho(rho0)
c$$$cc        poloidlen=rho0_pol*totlength*100.d0*r0x! poloidal length cm
c$$$c         poloidlen=rho0_pol*totlength*100.d0 ! poloidal length cm
c$$$cc         write(*,*)'prep3d totlength,poloidlen',totlength,poloidlen
c$$$
c$$$cSm990901
c$$$c	 binvol(i)=voltot*(rhov(hro*i)**2-rhov(hro*(i-1))**2)
c$$$c     1          *1000000.d0
c$$$cc        write(*,*)'prep3d voltot,rhov(hro*i),rhov(hro*(i-1)),binvol',
c$$$cc     +  voltot,rhov(hro*i),rhov(hro*(i-1)),binvol(i)
c$$$

c$$$cc       write(*,*)'prep3d i ,binvol,poloidlen',
c$$$cc     + i,binvol(i),poloidlen
c$$$     	          
c$$$c      enddo
c$$$      write(*,*)'INT testing DNONETWO. powertot=erg/sec',powertot,
c$$$     1' currtot=A',currtot
c$$$      write(*,*)'testing 2 DNONETWO powtot_e,powtot_i,powtot_cl',
c$$$     1 powtot_e,powtot_i,powtot_cl

      powrtot1=0.d0
      currtot1=0.d0
      
cSmirnov970106 beg
      pwtot_e=0.d0
      pwtot_i=0.d0
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            pwtot_s(kk)=0.d0
         enddo
      endif
      pwtot_cl=0.d0
cSmirnov970106 end
      do i=1,NR-1
         powrtot1=powrtot1+spower(i)
cSmirnov970106 beg
         pwtot_e=pwtot_e+spower_e(i)
         pwtot_i=pwtot_i+spower_i(i)
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               pwtot_s(kk)=pwtot_s(kk)+spower_s(i,kk)
            enddo
         endif
              
         pwtot_cl=pwtot_cl+spower_cl(i)
         
cSmirnov970106 end
         currtot1=currtot1+scurrent(i)             !totaql toroidal current        
      enddo
      
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added      
      write(*,*)'INT2 testing DNONETWO. powrtot1=erg/sec',powrtot1,
     1' currtot1=A',currtot1
      
      write(*,*)'testing 2 DNONETWO pwtot_e,pwtot_i,pwtot_cl',
     &     pwtot_e,pwtot_i,pwtot_cl
      write(*,*)'pwtot_s(1:nbulk)=',(pwtot_s(kk),kk=1,nbulk)

      write(*,*)'powden: ',powden
      write(*,*)'powden_e: ',powden_e
      write(*,*)'powden_i: ',powden_i
      do kk=2,nbulk
         write(*,*)'kk,powden_s(,kk): ',kk,(powden_s(i,kk),i=1,NR-1)
      enddo
      endif ! outprint
         

c -----------------------------------------------------------------
c     normalization of density profiles
c     it dependents on numerical integration formulas
c     In the case  ***INT**** we have the following normalization:

cyup      write(*,*)'NR',NR

cSAP091021
c      do i=1,NR
      do i=1,NR-1
cHarvey970112
         if (powertot.ne.0.d0) then
             powden(i)=powden(i)*powrtot1/powertot
	 endif
cSmirnov970106 beg
         if (powtot_e.ne.0.d0) then
	   powden_e(i)=powden_e(i)*pwtot_e/powtot_e
	 endif
         if (powtot_i.ne.0.d0) then
           powden_i(i)=powden_i(i)*pwtot_i/powtot_i
	 endif
         if (iabsorp.eq.3) then
         do kk=2,nbulk
            if (powtot_s(kk).ne.0.d0) then
               powden_s(i,kk)=powden_s(i,kk)*pwtot_s(kk)/powtot_s(kk)
            endif
         enddo
         endif
         if (powtot_cl.ne.0.d0) then

           powden_cl(i)=powden_cl(i)*pwtot_cl/powtot_cl

	 endif
cSmirnov970106 end
         if (currtot.ne.0.d0) then
            currden(i)=currden(i)*currtot1/currtot
         endif
        
      enddo
c------------------------------------------------------------------

      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added

      do kk=2, nbulk
         write(*,*)'prep3d: second powden_s(i,kk)=', 
     1        (powden_s(i,kk),i=1,NR-1)
      enddo


      write(*,998)
 998  format(//)

cSAP090306
c      write(*,999)powtott
c 999  format(' total injected power(erg/sec)   =',1pe14.6)
      write(*,999)powtott*1.d-7
 999  format(' total injected power(watt)   =',1pe14.6)
      write(*,1005)i_total_bad_initial_conditions
cSAP090306  
c      write(*,1006)power_launched
c      write(*,1000)powrtot1
c 1000 format(' total absorbed power(erg/sec)   =',1pe14.6)
c      write(*,1001)pwtot_e
c 1001 format(' total absorbed power_e(erg/sec) =',1pe14.6)
c      write(*,1002)pwtot_i
c 1002 format(' total absorbed power_i(erg/sec) =',1pe14.6)
      write(*,1006)power_launched*1.d-7
      write(*,1000)powrtot1*1.d-7
 1000 format(' total absorbed power(watt)   =',1pe14.6)
      write(*,1001)pwtot_e*1.d-7
 1001 format(' total absorbed power_e(watt) =',1pe14.6)
      write(*,1002)pwtot_i*1.d-7
 1002 format(' total absorbed power_i(watt) =',1pe14.6)
      if (iabsorp.eq.3) then
cSAP090306
c         write(*,10021)(kk-1,pwtot_s(kk),kk=2,nbulk)
c10021    format(' total absrbd power_s(erg/sec) for ion species',
c     1       i2,' =',1pe14.6)
        write(*,10021)(kk-1,pwtot_s(kk)*1.d-7,kk=2,nbulk)
10021    format(' total absrbd power_s(power) for ion species',
     1       i2,' =',1pe14.6)
      endif
cSAP090306
c      write(*,1003)pwtot_cl
c 1003 format(' total absorbed power_cl(erg/sec)=',1pe14.6)
      write(*,1003)pwtot_cl*1.d-7
 1003 format(' total absorbed power_cl(power)=',1pe14.6)

cSAP080928
c      write(*,1004)currtot1
 1004 format(' total tor curr drive per LL memo[Cohen/(1+eps)] (A)=',
     1     1pe14.6)
 1005 format(' number of rays with bad initial conditions ='1i6)
cSAP090306
c 1006 format(' total launched power with good initial conditions',/,
c     1 '                        (egs/sec)=',1pe14.6)
 1006 format(' total launched power with good initial conditions',/,
     1 '                        (watt)=',1pe14.6)
ctest      
c      write(*,*)' testing DNONETWO. allpower=',allpower

cSAP090306
      write(*,1007)parallel_cur_total
      write(*,1008)toroidal_cur_total
      write(*,1009)poloidal_cur_total
 1007 format('total parallel current parallel_cur_total (A)=',
     &    1pe14.6)
 1008 format('total toroidal current toroidal_cur_total (A)=',
     &    1pe14.6)
 1009 format('total poloidal current poloidal_cur_total (A)=',
     &    1pe14.6)


        if(i_lsc_approach.eq.1)then
         write(*,*)'NR-1,n_psi_TSC',NR-1,n_psi_TSC

         write(*,1030)'i rho_bin_center powden currden '

         total_power=0.d0
         do i=1,NR-1
            write(*,1031)i,rho_bin_center(i),powden(i),currden(i)

            total_power=total_power+
     &               powden_e(i)*
     &               binvol(i)*1.d-6*1.d-1
         enddo
         write(*,*)'in dnonetwo total_power [watt] = ',total_power

          write(*,1030)'i,rho_bin_center_lsc(i) '
     &      //'power_dens_watt_m3_TSC_1D(i)'
     &      //'CD_dens_no_Edc_a_m2_TSC_1D(i)'

c         do i=1,n_psi_TSC
c            write(*,1031)i,rho_bin_center_lsc(i),
c     &      power_dens_watt_m3_TSC_1D(i),
c     &      CD_dens_no_Edc_a_m2_TSC_1D(i)
c         enddo

1030    format(/,1x,a)
1031    format(i3,3(1pe12.4))

        endif ! i_lsc_approach
      
      endif ! outprint

      return
      END




  
      subroutine grpde2_old_real4_110211                     !Not called
     &(npar, nper, omega, omegpe, nsigma,                    !Not called
     .                   vgrpdc, edenfac, exde, eyde, ezde)  !Not called
c
      complex          rootm1,exde,eyde,ezde,d11,d12,d13,d22,d33,denom
      real             npar,npar2,nper,nper2,n2,nperol
****  double precision deps1,deps2,deps3,domega,domegpe,
****  double precision dquad
      real             dquad
      real deps1,deps2,deps3,domega,domegpe,aquad,bquad,cquad
c
c     calculate group velocity, wave polarizations, and energy density factor,
c     using cold plasma relations. (M O'Brien's version)
c
      nperol = nper
      rootm1 = (0.0, 1.0)
c
c     DOUBLE PRECISION is used since a problem was encountered in the
c                 real is used since a problem was encountered in the
c     following subtraction within the formation of rootqd.
c
c --- double precision mentioned above was disabled 21 Aug 94 by Joe Freeman
c --- not needed since the WHOLE PROGRAM is done with 64-bit real arithmetic
c --- on CRAY this is automatic, elsewhere a compilation switch ensures this
c
      domega  = omega
      domegpe = omegpe
c      write(*,*)'in grpde2 npar,domega,domegpe,nsigma',
c     1                     npar,domega,domegpe,nsigma
      deps1   = 1.0 - domegpe*domegpe/        (domega*domega-1.0)
      deps2   =      -domegpe*domegpe/(domega*(domega*domega-1.0))
      deps3   = 1.0 - domegpe*domegpe/        (domega*domega)
      eps1    = deps1
      eps2    = deps2
      eps3    = deps3
c
      npar2 = npar*npar
      aquad = deps1
      bquad = - ( (deps1-npar2)*(deps1+deps3) - deps2*deps2 )
      cquad = deps3 * ( (deps1-npar2)*(deps1-npar2) - deps2*deps2 )
      dquad = bquad*bquad - 4.0*aquad*cquad
      
      if (dquad .gt. 0.0) then
        rootqd = SQRT (dquad)
      else
        rootqd = 0.0
      end if
      nper2=(-bquad+nsigma*SIGN(1.0,omega-1.0)*rootqd)/(2.0*aquad)

c
c     Let nper2 become (only) slightly negative:
c
      if (nper2 .le. 0.0) then
        if (nper2 .le. -1.0e-2) then
          STOP 'subroutine GRPDE2_old: nper2 is too negative'
        else
          nper2=1.0e-10
          nper=1.0e-5
        end if
      else
        nper = SQRT (nper2)
      end if
c
      n2=nper2+npar2

      d11=eps1-npar2
      d12=-rootm1*eps2
      d13=npar*nper
      d22=eps1-n2
      d33=eps3-nper2
c 
      
      exde=d22/d12
      ezde=-d13*exde/d33
      denom = SQRT (exde*conjg(exde)+1.0+ezde*conjg(ezde))
      exde=exde/denom
      eyde=1.0/denom
      ezde=ezde/denom
c      write(*,*)'in grpde2  exde,eyde,ezde',exde,eyde,ezde

c
c calculate derivatives of eps's wrt omega
c
      dep1do=omegpe*omegpe*2.0*omega/
     .        ((omega*omega-1.0)*(omega*omega-1.0))
      dep2do=-eps2/omega-eps2*2.0*omega/(omega*omega-1.0)
      dep3do=omegpe*omegpe*2.0/(omega**3)
c
c find derivatives of D the dispersion reln wrt omega, kper and kpar
c
      dddnpp=-(2.0*nper)*( (eps1-npar2)*(eps1-n2+eps3-nper2)
     .                   - eps2*eps2 + npar2*(eps1-n2-nper2) )
      dddnpl=-(2.0*npar)*( (eps3-nper2)*(eps1-npar2+eps1-n2)
     .                   + nper2*(eps1-n2-npar2) )
c
      dddome=dddnpp*(-nper/omega) + dddnpl*(-npar/omega)
     .     + dep1do*( (eps3-nper2)*(eps1-n2+eps1-npar2) - nper2*npar2 )
     .     + dep2do*( -2.0*eps2*(eps3-nper2) )
     .     + dep3do*( (eps1-npar2)*(eps1-n2) - eps2*eps2 )
c
c find vgroup/c
c
      vgrpdc = SQRT (dddnpp*dddnpp+dddnpl*dddnpl) / ABS (dddome*omega)
      
c     energy density factor
c
c      1  | ~ | 2    1  ~*      ~       ~     ~
c   =  -  | B |    + -  E .tens.E  with B and E normalized to |E|
c      2  |   |      2       _      _
c                         d |     h  |    h
c    and tens the tensor  - |  w e   |   e  is the hermitian component
c                         dw|_      _!
c
c    of the dielectric tensor
c
      bmod2 = REAL ( n2*eyde*conjg(eyde)
     .    + (npar*exde-nper*ezde)*conjg(npar*exde-nper*ezde) )
c      write(*,*)'in grpde2  bmod2',bmod2
      tens11=eps1+omega*dep1do
      tens12=-(eps2+omega*dep2do)
      tens21=-tens12
      tens22=eps1+omega*dep1do
      tens33=eps3+omega*dep3do
      emod2 = REAL (conjg(exde)*exde*tens11
     .            + conjg(exde)*eyde*rootm1*tens12
     .            + conjg(eyde)*exde*rootm1*tens21
     .            + conjg(eyde)*eyde*tens22
     .            + conjg(ezde)*ezde*tens33)
c      write(*,*)'in grpde2  emod2',emod2
c
      edenfac=0.5*(bmod2+emod2)
c      write(*,*)'in grpde2  edenfac',edenfac

c
c     nper=nperol
      return
c
      end


c       *******************testel*******************
c       test of polarization and fluxn calculations
c----------------------------------------------------
c       input: cnpar,cnpert -refractive index components
c              xe,ye
c       output:ex,ey,ez,flux
c----------------------------------------------------
      subroutine testel(cnpar,cnper,xe,ye,ex,ey,ez,flux)
      implicit double precision (a-h,o-z)
c------------------------------------------------------
c			| s , -id,  0|
c     dielectric tensor=| id,	s,  0|
c			| 0 ,	0,  p|
c-----------------------------------------------------
      cn2=cnper*cnper+cnpar*cnpar
      cs=cnpar/dsqrt(cn2)
      cos2=cs*cs
      sn=cnper/dsqrt(cn2)
      sin2=sn*sn
      cn2=cnper*cnper+cnpar*cnpar
      ye2=ye*ye
      s=1.d0-xe/(1.d0-ye2)
      d=-xe*ye/(1.d0-ye2)
      p=1.d0-xe
      write(*,*)'in testel xe,ye,cos2,sin2,cn2'
      write(*,*)xe,ye,cos2,sin2,cn2
      write(*,*)'in testel s,d,p'
      write(*,*)s,d,p
c------------------------------------------------------
      ez=(s-cn2)/d*cn2*sn*cs/(p-cn2*sin2)
      ex=-(s-cn2)/d
      emod2=1.d0+ez*ez+ex*ex
      ey2=1.d0/emod2
c----------let .ge.0
      if(ez.lt.0.d0) then
        ey=-dsqrt(ey2)
      else
        ey=-dsqrt(ey2)
      endif
c-----------------------------
      ex=ex*ey
      ez=ez*ey
      write(*,*)'in testel real_ex imag_ey real_ez'
      write(*,*)ex,ey,ez
c-----------------------------
      dwsdw=s+2.d0*xe/((1.d0-ye2)*(1.d0-ye2))
      dwddw=d+xe*ye*(3.d0-ye2)/((1.d0-ye2)*(1.d0-ye2))
      write(*,*)'1 dwddw',dwddw,'dwsdw',dwsdw
      dwddw=2.d0*xe*ye/((1.d0-ye2)*(1.d0-ye2))
      write(*,*)'2 dwddw',dwddw
      dwpdw=1.d0+xe
c------------------------------------------------------
      fluxel1=(ex*ex+ey*ey)*dwsdw
      fluxel2=ez*ez*dwpdw
      fluxel3=2.d0*ex*ey*dwddw
      write(*,*)'ex,ey,ex*ex+ey*ey',ex,ey,ex*ex+ey*ey
      write(*,*)'fluxel1',fluxel1
      write(*,*)'fluxel21',fluxel2
      write(*,*)'fluxel3',fluxel3
      fluxel=fluxel1+fluxel2+fluxel3
      write(*,*)'fluxel',fluxel
c------------------------------------------------------
      cbx=-cnpar*ey
      cby=-cnper*ez+cnpar*ex
      cbz=cnper*ey
      write(*,*)'cbx,cby,cbz',cbx,cby,cbz
      bmodb=cbx*cbx+cby*cby+cbz*cbz
      write(*,*)'in testel bmodb',bmodb
      flux=bmodb+fluxel
      write(*,*)'in testel flux',flux
      return
      end
      
      
      
      subroutine efKarney(z,r,r0x,z_eff,temp,den,jwave,cnpar,effKarn)
c------------------------------------------------------------------
c  RF current drive efficiency(asymptotic formula, nonrelativistic case)  
c  D.A. Ehst and C.F.F.Karney Nucl.Fus. Vol.31, No.10 (1991), p 1933-1938
c  A subroutine to calculate the local current density driven by RF waves
c--Uses CD efficiency empirical formula based on numerical Fokker-
c  Planck bounce-averaged calculations
c  This formula is for
c  1)Landau damping of lower hybrid (LH) slow waves resonant at parallel
c    velocities	above the electron thermal velocity
c  2)Slow frequency fast (compressional Alfven) wave (AW) may resonant with
c    low phase velocity electrons via combined Landay damping and transit
c    time magnetic damping  
c------------------------------------------------------------------
c     input parameters: z,r -coordinates of the ray point(cm!!!!!)ATTENTION
c                       z_eff-effective charge
c                       temp-temperature (keV)
c                       den-density(10**13 cm**-3)
c                       jwave (=islofa)=-1 Alfven wave , 0- Landau damp.
c                       cnpar -paralell to magnetic field refractive
c                              index
c                       r0x character length
c     computed:  J/P=efficient in (A/m**2)/(W/m**3), converted to 
c     output parameter: efficien J/P in (A/cm**2)/(erg/(sec*cm**3))
!     Conversion is 1.0(A/m^2)/(W/m^3) = 1e-5*(A/cm**2)/(erg/(sec*cm**3))
c------------------------------------------------------------------
c     It uses:
c     double precision function: psif(z,r) and
c     subroutine: zr_psith(psi,theta,z,r)
c     double precision functions from zr_psith.f: rmax_psi(psi),
c     rmin_psi(psi), bmax_psi(psi)
c------------------------------------------------------------------
      IMPLICIT double precision (a-h,o-z)
      islofa=jwave
      zeff=z_eff

c      write(*,*)'in efKarney jwave,z_eff',jwave,z_eff
      if(zeff.le.0.d0) then
         WRITE(*,*)'efKarney: illegal INPUT z_eff=',zeff
         STOP
      endif
      if(r0x.le.0.d0) then
         WRITE(*,*)'efKarney: illegal INPUT r0x=',r0x
         STOP
      endif
      if(temp.le.0.d0) then
         WRITE(*,*)'efKarney: illegal INPUT temp=',temp
         STOP
      endif
      if(den.le.0.d0) then
         WRITE(*,*)'efKarney: illegal INPUT den=',den
         STOP
      endif
      if(cnpar.eq.0.d0) then
         WRITE(*,*)'efKarney: illegal INPUT cnpar=',cnpar
         STOP
      endif

      if(islofa.eq.-1) then	!Alfven damping  !Note: islofa=jwave
        akcd = 11.91d0/(0.678d0+zeff)
        c0cd = 4.13d0/zeff**0.707d0
        amcd = 2.48d0
        ccd = 0.0987d0
        acd = 12.3d0
      elseif(islofa.eq.0)then !Landau damping: jwave=0 (n=0 case only)
        !YuP[2021-02-11] Added islofa.eq.0, and warning for any other jwave
        akcd = 3.d0/zeff
        c0cd = 3.83d0/zeff**0.707d0
        amcd = 1.38d0
        ccd = 0.389d0
        acd = 0.d0
      else !YuP[2021-02-11] 
        !A warning message is printed in the beginning of run.
        !Use settings as in jwave=0 case (although not physically meaningful)
        akcd = 3.d0/zeff
        c0cd = 3.83d0/zeff**0.707d0
        amcd = 1.38d0
        ccd = 0.389d0
        acd = 0.d0
      endif
c      write(*,*)'in efKarney (K,D,m,c,a) kcd,c0cd,amcd,ccd,acd',
c     &                       akcd,c0cd,amcd,ccd,acd
c--------------------------------------------------------------------
c     aspct is inverse aspct ratio
c--------------------------------------------------------------------
c     bmax,bmin  in toray
c     if(dabs(bmax-bmin).lt.1.0d-8) bmax=bmin+1.0d-8
c     aspct=(bmax-bmin)/(bmax+bmin)
c--------------------------------------------------------------------
c     Here :aspct is calculating using rmax_psi and rmin_psi
c     conversion from cm to m
      zd=z*0.01d0/r0x
      rd=r*0.01d0/r0x
c----------------------------------------
c      write(*,*)'in efKarney zd,rd',zd,rd
      psid=psif(zd,rd)
c      write(*,*)'efKarney psid',psid
c----- rmaxpsi and rminpsi are the largest and and the least
c     values of the radius (rmaxpsi and rminpsi in (m))
c     on the given flux surface psi=psi(zd,rd)
      rmaxpsi=rmax_psi(psid)
      rminpsi=rmin_psi(psid)
      if (rmaxpsi.lt.rd)rmaxpsi=rd
      if (rminpsi.gt.rd)rminpsi=rd
      if(dabs(rmaxpsi-rminpsi).lt.1.0d-8) rminpsi=rmaxpsi-1.0d-8
      aspct=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
      epsil=aspct
 
c      write(*,*)'efKarney eplis',epsil
c----------------------------------------
      phi=0.d0 !   ????
      bmod=b(zd,rd,phi)
      bmax=bmax_psi(psid)
      if(bmax.eq.0.d0) then
         WRITE(*,*)'efKarney: illegal value bmax=',bmax
         STOP
      endif
      if (bmod.gt.bmax) bmod=bmax
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(Te/me)),    ve in (cm/sec),Te(keV)     !article p.1935
      ve=dsqrt(temp*1.6022d-9/9.1094d-28) !(sqrt(T_e/m_e) electron thermal
                                          !velocity [cm/sec]
c      write(*,*)'Ehst Karney test ve',ve
c      ve=1.32d9*dsqrt(temp)	
c                      !temp keV
c      write(*,*)'Ehst Karney ve',ve

c     normalized resonance parallel velocity:
c     wte=u_res=v_par/v_e=cvac/(ve*cnpar)
      wte=cvac/(ve*cnpar)			      !article p.1934

c      write(*,*)'cvac,temp,ve,cnpar,wte',cvac,temp,ve,cnpar,wte
c      write(*,*)'Ehst Karney wte',wte
      u1=wte/dsqrt(2.d0) 
cSAP080905
      wte=dabs(wte)
c----------------------------------------
      alt2=1.d0-bmod/bmax			      !article p.1935 
      alt2=dabs(alt2)
      alt1 =dsqrt(alt2)
c     write(*,*)'efKarney (labmda_t) alt1',alt1
c----------------------------------------
      if (alt2.ne.0.d0) then
         ytt = (1.d0-alt2)*wte**2/alt2		      !article (11)
         rprof = 1.0d0-(epsil**0.77d0*dsqrt(12.25d0+wte**2))/
     1  (3.5d0*epsil**0.77d0+wte)					      !article (7)

c         write(*,*)'efKarney epsil,wte,epsil**0.77d0',
c     &                       epsil,wte,epsil**0.77d0
c         write(*,*)'efKarney epsil**0.77d0*dsqrt(12.25d0+wte**2))',
c     &                       epsil**0.77d0*dsqrt(12.25d0+wte**2)
c         write(*,*)'efKarney (3.5d0*epsil**0.77d0+wte)',
c     &                       (3.5d0*epsil**0.77d0+wte)
c         write(*,*)'efKarney in R second term',
c     &  (epsil**0.77d0*dsqrt(12.25d0+wte**2))/
c     &  (3.5d0*epsil**0.77d0+wte)

         arg = (ccd*ytt)**amcd
         cprof = 1.d0-dexp(-arg)	      	       !article (9)
         amprof = 1.d0+acd*(alt1/wte)**3	       !article (10)
      else
        rprof=1.d0
        cprof=1.d0
        amprof=1.d0
      endif
c      write(*,*)'efKarney (R) rprof',rprof
      eta0 = akcd/(dabs(wte))+c0cd+4.d0*wte**2/(5.d0+zeff) !article (14)
      eta=cprof*amprof*eta0*rprof 		      !article (13)
c      write(*,*)'4.d0*wte**2/(5.d0+zeff)',4.d0*wte**2/(5.d0+zeff) 
c      write(*,*)'eta,cprof,amprof,eta0,rprof ',
c     &           eta,cprof,amprof,eta0,rprof 
c       write(*,*)'eta0,eta',eta0,eta
c-----------------------------------------------------------
c     efficiency  in (A/m**2)/(W/m**3)
c     temperature temp in kev
c     density     dens in 10**19 /m**3
c-----------------------------------------------------------
c     cln -Coulomb Ln
      arg1 = 1.d3/temp*dsqrt(10.d0*den)	!temp KeV, den 10**13/cm**3
      cln=24.d0-dlog(arg1)
      effKarn=eta*3.84d0*temp/(cln*den)	!(A/m**2)/(W/m**3)
c      write(*,*)'temp,den,cln,eta,effKarn',temp,den,cln,eta,effKarn

c-----------------------------------------------------------
c     efficiency will be converted to (A/cm**2)/(erg/(sec*cm**3))
c     temperature temp in kev
c     density     den in 10**19/m**3
c-----------------------------------------------------------
ctest beg
c     LH wave
c      efficien=2.d0/(5.d0+z_eff)*
c     1         (u1*u1+(7.0d0/4.0d0+9.0d0/4.0d0/(3.d0+z_eff)))
c      efficien=efficien*temp/den*(17.d0/cln)*9.d0*1.d-6
c      write(*,*)'u1',u1  
c     write(*,*)'asimptotic efficiency (A/cm**2)/(erg/(sec*cm**3))',efficien
c      efficien=8.d0/(5.d0+z_eff)*u1*u1 !test
c      write(*,*)'test efficien ',efficien
c
c      efficien=efficien*temp/den*(17.d0/cln)*9.d0*1.d-6 
c     write(*,*)'efficiency (A/cm**2)/(erg/(sec*cm**3))',efficien
c      efficien=efficien*3.84d0*temp/den/cln !test only
c      write(*,*)'asymptotic efficiency(A/m**2)/(W/m**3))',efficien
ctest end    
c-------------------------------------------------------------
c     determination the of the current drive sign
      if (u1.ge.0.d0) then
         s=1.d0
      else
         s=-1.d0
      endif
      efficien=efficien*s
c-------------------------------------------------------------
cSm050923
c      effKarn=effKarn*s*1.d-5  !  (A/cm**2)/(erg/(sec*cm**3))
      effKarn=-effKarn*s*1.d-5  !  (A/cm**2)/(erg/(sec*cm**3))
c      write(*,*)'effcient A/cm**2/(egr/(sec*cm**3))',efficien
c      write(*,*)'effKarn  A/cm**2/(erg/(sec*cm**3))',effKarn
c     stop
      return
      END subroutine efKarney



       subroutine map_dhot_nper(m_r_nperp,m_i_nperp, 
     .max_r_nperp,min_r_nperp,max_i_nperp,min_i_nperp,
     .nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .n_contour_plot_disp,
     .vflow_ar,nll_in,re_nperp0,im_nperp0,iabsorp,
     .n_relt_harm1,n_relt_harm2,n_relt_intgr,
     .i_resonance_curve_integration_method,epsi,
     .i_diskf,r,z,phi)
     
c     calculate the map for the complex hot plasma dispersion function
c     on the 2D mesh (ReN_perp,ImN_perp)
c
      implicit none
c     input u
      integer m_r_nperp,m_i_nperp, !the number of points of the 2D-mesh
                                   ! for (ReN_perp,ImN_perp)
     .n_contour_plot_disp          ! the number of contours

      double precision max_r_nperp,min_r_nperp!max and min ReN_perp
      double precision max_i_nperp,min_i_nperp !max and min ImN_perp
c      nbulk - the total number of the plasma species 
c      mass_ar - the masses of the plasma  species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=average temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar = T_perp/T_parallel  
c      vflow_ar    cm/sec
c      nll_in - parallel index of refraction n.
c      re_nperp0   Re(N_perp) the data for N_perp
c      im_nperp0    Im(N_perp) calculated along the ray (for information)
c      iabsor =4 hot plasma dispersion
c             =6 for Hermitian hot plasma tensor +
c                anti-Hhrmihian relativistic electron tensor
c     the data for relativistic anti Hermitian tensor
c     n_relt_harm1,n_relt_harm2
c     n_relat_intgr
c     i_diskf
c     z,r,phi the space coordinates
c-------------------------------------------------------------------
c       i_resonance_curve_integration_method=1 !rectangle 
c                                              !integration other angle
c                                              !for ellipse case only
c       i_resonance_curve_integration_method=2 !rectangle integration
c                                              !other p_perp integration
c                                              !Ellipse and hyperbola
c       i_resonance_curve_integration_method=3 !trapezoidal integration
c                                              !other p_perp integration
c                                              !Ellipse and hyperbola
c       i_resonance_curve_integration_method=4 !adaptive Simpson 
c                                              !other p_perp integration
c                                              !Ellipse and hyperbola
c
c       i_resonance_curve_integration_method is used in subroutine intgr_rl
c       to choose the numerical integration method 
c       for anti-hermitian relativistic dielectric tensor calculation
c-----------------------------------------------------------------------
c       epsi  absolute accuracy used in adaptive Simpson 
c----------------------------------------------------------------------
      integer nbulk,iabsorp,i_diskf,n_relt_harm1,n_relt_harm2,
     &n_relt_intgr,i_resonance_curve_integration_method
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double precision nll_in,z,r,phi
      double precision re_nperp0,im_nperp0,epsi
     
c-----external
      double complex dhot_sum     !from RealN_perp)
      double complex dhot_sum_c   !from complex N_perp=(nperp,nprim)
c     dhot_rlt complex dispersion function from Hermitian hot plasma +
c     anti-Hermitin electron tensor and comple N_perp=(nperp,nprim)
      double complex dhot_rlt  
c-----local
      integer i,j,i_fkin
      integer iherm   !=1 hermitian dielectric tensor, 2-full
      double precision step_real,step_imag,nperp,nprim,ye
      double complex K_sum(3,3),d_herm,d_full,d_compl
      double complex aK(3,3) !relativistic anti Hermitian electron tensor 
      integer n_param
      parameter(n_param=6)
      character*(10)name_param(n_param) !names of the input parameters
      real param(n_param)               !values of the input parameters
      
c-----for pgplot

c-----parameter n_contour_plot_disp is the number of contours           
c      parameter(m_r_nperp_a=10,m_i_nperp_a=10,n_contour=10)
      include 'param.i'

      real x_axis(m_r_nperp_a),y_axis(m_i_nperp_a),
     .f2d_r(m_r_nperp_a,m_i_nperp_a),f2d_i(m_r_nperp_a,m_i_nperp_a),
     .f2d_abs(m_r_nperp_a,m_i_nperp_a),
     .contour(n_contour_plot_disp_a)

cSm060919
      if (m_r_nperp_a.lt.m_r_nperp) then 
         write(*,*)'prep3d.f in map_dhot_nper m_r_nperp_a.lt.m_r_nperp'
         write(*,*)'m_r_nperp_a=',m_r_nperp_a
         write(*,*)'m_r_nperp=',m_r_nperp
         write(*,*)'change m_r_nperp_a in prep3d.f in' 
         write(*,*)'subroutine map_dhot_nper and recompile the code'
         stop 'prep3d.f in  subroutine map_dhot_nper'
      endif

      if (m_i_nperp_a.lt.m_i_nperp) then 
         write(*,*)'prep3d.f in map_dhot_nper m_i_nperp_a.lt.m_i_nperp'
         write(*,*)'m_i_nperp_a=',m_i_nperp_a
         write(*,*)'m_i_nperp=',m_i_nperp
         write(*,*)'change m_i_nperp_a in prep3d.f in' 
         write(*,*)'subroutine map_dhot_nper and recompile the code'
         stop 'prep3d.f in  subroutine map_dhot_nper'
      endif



c      write(*,*)'in map_dhot_nperp:  nbulk,nll_in',nbulk,nll_in
c      write(*,*)'x_ar',(x_ar(i), i=1,nbulk)
c      write(*,*)'y_ar',(y_ar(i), i=1,nbulk)
c      write(*,*)'t_av_ar',(t_av_ar(i), i=1,nbulk)

      step_real=(max_r_nperp-min_r_nperp)/dble(m_r_nperp-1) !step ReN_perp
      step_imag=(max_i_nperp-min_i_nperp)/dble(m_i_nperp-1) !step ImN_perp
c      write(*,*)'before open 15'
c      open(15,file='dhot_1D.dat')
c 1    format(5(1pe16.7))
c      open(16,file='dhot_2D.dat')
c 2    format(5(1pe16.7))
c      write(15,*)'nperp  ReD_herm ImD_herm ReD_full ImD_full'
c      write(16,*)'ReNperp ImNperp ReDc ImDc absDc'
      do i=1,m_r_nperp  !Ends line 5619
         nperp=min_r_nperp+step_real*(i-1) !Real N_perp
         x_axis(i)=nperp
c--------Hermitian hot plasma dispersion function
         iherm=1          
         d_herm=dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .   vflow_ar,nll_in,nperp,iherm,K_sum)

         if (iabsorp.eq.4) then ! full hot plasma dispersion function ImN_per=0
            iherm=2
            d_full=dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .      vflow_ar,nll_in,nperp,iherm,K_sum)
         endif

         if (iabsorp.eq.6) then !Ends line 
c-----------the dispersion function for ImN_per=0
c           for Hermitian hot plasma tensor + relativistic
c           anti Hermitian electron tensor
c
c           Hermitian hot plasma tensor K_sum was calculated previously 
c           before if(iabsorp.eq.4) 
c
c-----------anti-hermitian relativistic tensor aK   
c         
c           n_relt_harm1 min number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations
c           n_relt_harm2 max number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations
c               n_relt_harm1 =< n= <n_relt_harm2
c           n_relt_intgr is the number of points for the numrerical integration
c           over p_perp for the calculations of anti-hermitian
c           relativistic tensor aK.
        
           if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c--------------usage of the analytical relativistic function and its derivatives 
           i_fkin=0
               else 
c--------------usage of the mesh relativistic function and its derivatives
               i_fkin=1
           endif  !On iabsorp.eq.6
cSm060314 
c           ye=dabs(y_ar(1)) !positive ye for positive harmonics 
c           call anth_rlt(x_ar(1),ye,t_av_ar(1)*1.d-3,nll_in,nperp,
           call anth_rlt(x_ar(1),y_ar(1),t_av_ar(1)*1.d-3,nll_in,nperp,
     +     n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +     i_resonance_curve_integration_method,epsi,
     +     i_fkin,r,z,phi,aK)

           nprim=0.d0
           d_full=dhot_rlt(K_sum,aK,nll_in,nperp,nprim)
         endif

c        write(15,1)nperp,dreal(d_herm),dimag(d_herm),
c     .  dreal(d_full),dimag(d_full)

         do j=1,m_i_nperp
            nprim=min_i_nperp+step_imag*(j-1)
            y_axis(j)=nprim

            if (iabsorp.eq.4) then 
c--------------full hot plasma dispersion function for Im_N_per.ne.0
               iherm=2
               d_compl=dhot_sum_c(nbulk,mass_ar,x_ar,y_ar,t_av_ar,
     .         tpop_ar,vflow_ar,nll_in,nperp,nprim,iherm,K_sum)
            endif

            
            if (iabsorp.eq.6) then 
c--------------the dispersion function for ImN_per .ne. 0
c              for Hermitian hot plasma tensor + relativistic
c              anti Hermitian electron tensor
               d_compl=dhot_rlt(K_sum,aK,nll_in,nperp,nprim)
            endif

c            write(16,2)nperp,nprim,dreal(d_compl),dimag(d_compl),
c     .      dsqrt((dreal(d_compl))**2+(dimag(d_compl))**2)
            f2d_r(i,j)=real(dreal(d_compl))
            f2d_i(i,j)=real(dimag(d_compl))
            f2d_abs(i,j)=sqrt(f2d_r(i,j)**2+f2d_i(i,j)**2)
         enddo
      enddo  ! Begins line 5539
c      close(15)
c      close(16)

c      call plotinit
      name_param(1)='Y'
      name_param(2)='X'
      name_param(3)='N_parallel'
      name_param(4)='T'
      name_param(5)='ReNperp_0'
      name_param(6)='ImNperp_0'
      param(1)=y_ar(1)
      param(2)=x_ar(1)
      param(3)=nll_in
      param(4)=t_av_ar(1)
      param(5)=re_nperp0
      param(6)=im_nperp0
cSm030519
c      call contour2d(m_r_nperp,m_i_nperp,f2d_r,x_axis,y_axis,
c     .'real_D', 'ReN_perp','ImN_perp',
c     .n_contour_plot_disp,contour,
c     .name_param,param,n_param)

c      call contour2d(m_r_nperp,m_i_nperp,f2d_i,x_axis,y_axis,
c     .'ImageD', 'ReN_perp','ImN_perp',
c     .n_contour_plot_disp,contour,
c     .name_param,param,n_param)

c      call contour2d(m_r_nperp,m_i_nperp,f2d_abs,x_axis,y_axis,
c     .'absD', 'ReN_perp','ImN_perp',
c     .n_contour_plot_disp,contour,
c     .name_param,param,n_param)
cSm060920
      write(*,*)'prep3d.f in map_dhot_nper n_contour_plot_disp',
     &n_contour_plot_disp

      call contour2d_S(m_r_nperp_a,m_i_nperp_a,n_contour_plot_disp_a,
     .m_r_nperp,m_i_nperp,f2d_r,x_axis,y_axis,
     .'real_D', 'ReN_perp','ImN_perp',
     .n_contour_plot_disp,contour,
     .name_param,param,n_param)

      call contour2d_S(m_r_nperp_a,m_i_nperp_a,n_contour_plot_disp_a,
     .m_r_nperp,m_i_nperp,f2d_i,x_axis,y_axis,
     .'ImageD', 'ReN_perp','ImN_perp',
     .n_contour_plot_disp,contour,
     .name_param,param,n_param)

      call contour2d_S(m_r_nperp_a,m_i_nperp_a,n_contour_plot_disp_a,
     .m_r_nperp,m_i_nperp,f2d_abs,x_axis,y_axis,
     .'absD', 'ReN_perp','ImN_perp',
     .n_contour_plot_disp,contour,
     .name_param,param,n_param)

      return
      end


c!!!!!!!!!!!!!!!!!!!!
      subroutine map_dhot_nper_omeg(m_r_nperp,m_y, 
     .max_r_nperp,min_r_nperp,max_y,min_y,
     .nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .vflow_ar,nll_in)
     
c     calculate the map for the complex hot plasma dispersion function
c     on the 2D mesh (ReN_perp,Y)
c
      implicit none
c     input 
      integer m_r_nperp,m_y !the number of points of the 2D-mesh for (ReN_perp,Y)
      double precision max_r_nperp,min_r_nperp!max and min ReN_perp
      double precision max_y,min_y !max and min Y
c      nbulk - the total number of the plasma species 
c      mass_ar - the masses of the plasma  species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=average temperature=Te(1+2*tpop)/3 in ev   
c      tpop_ar = T_perp/T_parallel  
c      vflow_ar    cm/sec
c      nll_in - parallel index of refraction n


      integer nbulk   
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double precision nll_in,z,r,phi
     
c-----external
      double complex dhot_sum     !from RealN_perp)
      double complex dhot_sum_c   !from complex N_perp=(nperp,nprim) 
c-----local
      integer i,j,i_fkin
      integer iherm   !=1 hermitian dielectric tensor, 2-full
      double precision step_nperp,step_y,nperp,nprim,y_loc
      double complex K_sum(3,3),d_herm,d_full 

c     for pgplot
      integer m_r_nperp_a,m_y_a,n_contour
      parameter(m_r_nperp_a=100,m_y_a=100,n_contour=20)
       
      real x_axis(m_r_nperp_a),y_axis(m_y_a),
     .f2d_r(m_r_nperp_a,m_y_a),f2d_i(m_r_nperp_a,m_y_a),
     .f2d_abs(m_r_nperp_a,m_y_a),
     .contour(n_contour)

      integer n_param
      parameter(n_param=4)
      character*(10)name_param(n_param)
      real param(n_param)
     
      if (m_r_nperp_a.lt.m_r_nperp) then 
         write(*,*)'prep3d.f in map_dhot_nper m_r_nperp_a.lt.m_r_nperp'
         write(*,*)'m_r_nperp_a=',m_r_nperp_a
         write(*,*)'m_r_nperp=',m_r_nperp
         write(*,*)'change m_r_nperp_a in prep3d.f in' 
        write(*,*)'subroutine map_dhot_nper_omeg and recompile the code'
         stop 'prep3d.f in  map_dhot_nper_omeg'
      endif

      if (m_y_a.lt.m_y) then 
         write(*,*)'prep3d.f in map_dhot_nper m_y_a.lt.m_y'
         write(*,*)'m_y_a=',m_y_a
         write(*,*)'m_y=',m_y
         write(*,*)'change m_r_a in prep3d.f in' 
        write(*,*)'subroutine map_dhot_nper_omeg and recompile the code'
         stop 'prep3d.f in  map_dhot_nper_omeg'
      endif
    


      step_nperp=(max_r_nperp-min_r_nperp)/dble(m_r_nperp-1)
      step_y=(max_y-min_y)/dble(m_y-1)
      
      do i=1,m_r_nperp
         nperp=min_r_nperp+step_nperp*(i-1)
         x_axis(i)=nperp                
         do j=1,m_y
            y_loc=min_y+step_y*(j-1)
            y_axis(j)=1.d0/y_loc
c-----------Hermitian hot plasma dispersion function
            iherm=1        
            y_ar(1)=y_loc
            d_herm=dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .      vflow_ar,nll_in,nperp,iherm,K_sum)
            f2d_r(i,j)=real(dreal(d_herm))
c            f2d_i(i,j)=real(dimag(d_herm))
         enddo
      enddo

      name_param(1)='Y'
      name_param(2)='X'
      name_param(3)='N_parallel'
      name_param(4)='T'
      param(1)=y_ar(1)
      param(2)=x_ar(1)
      param(3)=nll_in
      param(4)=t_av_ar(1)

cSm030519      
      call contour2d(m_r_nperp,m_y,f2d_r,x_axis,y_axis,
     .'real_D', 'ReN_perp','1/Y',
     .n_contour,contour,
     .name_param,param,n_param)

      return
      end


      double complex function dcold_rlt(K,aK,nll,np)
c     calculates complex dispersion function dcold_rlt
c     using the cold plasma hermition tesor K and 
c     relativistic anti-hermition dielectric tensor aK
c     for electron plasma
c     input
c       K(3,3)   complex hermition tensor
c       aK(3,3)  complec anti-hermition relativistic tensor
c       nll      N_parallel
c       np       N_perpendicular
      implicit none
c     input
      double complex K(3,3), aK(3,3)
      double precision nll,np
c     local
      double complex Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz
      double precision nlls,nps
      double complex i

c      write(*,*)'dcold_rlt nll,np',nll,np		
c      write(*,*)'dcold_rlt K',K
c      write(*,*)'dcold_rlt aK',aK
      i=dcmplx(0.d0,1.d0)

      Kxx=K(1,1)+aK(1,1)*i
      Kxy=K(1,2)+aK(1,2)*i
      Kxz=K(1,3)+aK(1,3)*i
      Kyx=K(2,1)+aK(2,1)*i
      Kyy=K(2,2)+aK(2,2)*i
      Kyz=K(2,3)+aK(2,3)*i
      Kzx=K(3,1)+aK(3,1)*i
      Kzy=K(3,2)+aK(3,2)*i
      Kzz=K(3,3)+aK(3,3)*i

      nlls=nll*nll
      nps=np*np

c      write(*,*)'nlls,nps',nlls,nps 
c      write(*,*)'Kxx,Kxy,Kxz',Kxx,Kxy,Kxz 
c      write(*,*)'Kyx,Kyy,Kyz',Kyx,Kyy,Kyz 
c      write(*,*)'Kzx,Kzy,KXz',Kzx,Kzy,Kzz 
      dcold_rlt=(Kxx-nlls) * (Kyy-nlls-nps) * (Kzz-nps)
     .+  Kxy * Kyz * (Kzx+np*nll) 
     .+ (Kxz+np*nll) * Kyx * Kzy
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
     .-  Kzy * Kyz * (Kxx-nlls)
     .- (Kzz - nps) * Kyx * Kxy

c      write(*,*)'dcold_rlt',dcold_rlt

      return
      end

      double precision function dDcold(K,nll,np)
c     calculates the derivative from dD/d(N_perp) the electron cold plasma
c     dispersion function D
c
c     input
c       K(3,3)   complex cold plasma tensor
c       nll      N_parallel
c       np       N_perpendicular
      implicit none
c     input
      double complex K(3,3)
      double precision nll,np
c     local
      double complex Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz
      double precision nlls,nps

      nlls=nll*nll
      nps=np*np

c-----cold plasma dispersion function 
c      dcold_rlt=(Kxx-nlls) * (Kyy-nlls-nps) * (Kzz-nps) 
c     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
c     .- (Kzz - nps) * Kyx * Kxy

      Kxx=K(1,1)
      Kxy=K(1,2)
      Kxz=K(1,3)
      Kyx=K(2,1)
      Kyy=K(2,2)
      Kyz=K(3,3)
      Kzx=K(3,1)
      Kzy=K(3,2)
      Kzz=K(3,3)
 
      dDcold=(Kxx-nlls) * (-2.d0*np) * ((Kzz-nps)+(Kyy-nlls-nps)) 
     .- (nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
     .- (Kzx+np*nll) * (-2.d0*np) * (Kxz+np*nll)
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (nll)
     .- (-2.d0*np) * Kyx * Kxy
      
      return
      end


      subroutine check_monotonic_radius_in_ray_element(is,
     &i_non_monotonic,z_center,r_center,rho_center)
c-----------------------------------------------------------
c     check if the small radius is a monotonic function along 
c     the ray element [M_(is-1),M_is]
c     Here: M_is={z(is),r(is),phi(is)}
c-----------------------------------------------------------
      implicit none
      include 'param.i'
      include 'write.i'
c-----input
      integer is  !  >1 the number of ray element
c-----output
      integer i_non_monotonic   ! =0 the maximal and minimal rho are
                                !    at the boundaries of the ray element  
                                ! =1 the maximal or minimal rho are
                                !    inside the ray element
      real*8  z_center,r_center ! The coordinates of the extremum rho
                                ! at  i_non_monotonic=1 case
      real*8 rho_center         !extremal rho at i_non_monotonic=1 case
c-----locals
      real*8 psi_left,psi_right,rho_left,rho_right,h_z,h_r,
     &z_i,r_i,psi_i,rho_i,rho_min,rho_max

      integer n_steps,i_max,i_min,i
c-----externals
      real*8 psif, rhopsi

      psi_left=psif(wz(is-1),wr(is-1))
      psi_right=psif(wz(is),wr(is))

      rho_left=rhopsi(psi_left)
      rho_right=rhopsi(psi_right)

      n_steps=100
      h_z=(wz(is)-wz(is-1))/(n_steps-1)
      h_r=(wr(is)-wr(is-1))/(n_steps-1)

      rho_min = rho_left
      rho_max = rho_left
      i_max=0
      i_min=0
      do i=1,n_steps
        z_i=wz(is-1)+h_z*i
        r_i=wr(is-1)+h_r*i
        psi_i=psif(z_i,r_i)
        rho_i=rhopsi(psi_i)
       
        if (rho_i.gt.rho_max) then
           i_max=i
           rho_max=rho_i
        endif

        if (rho_i.lt.rho_max) then
           i_min=i
           rho_min=rho_i
        endif
      enddo

      if((i_min.gt.1).and.(i_min.lt.n_steps)) then
c--------the minimal small radius value is inside the ray element
         i_non_monotonic=1
         rho_center=rho_min
         z_center=wz(is-1)+h_z*i_min
         r_center=wr(is-1)+h_r*i_min
         goto 10
      endif

      if((i_max.gt.1).and.(i_max.lt.n_steps)) then
c--------the maximal small radius value is inside the ray element
         i_non_monotonic=1
         rho_center=rho_max
         z_center=wz(is-1)+h_z*i_max
         r_center=wr(is-1)+h_r*i_max
         goto 10
      endif

      i_non_monotonic=0

 10   continue

      return
      end     


c     this  subroutine calculates current drive efficiency
c     for the toroidal  plasma using TorGA_curgap codo written by Lin-Liu
c     ATTENTION:all parameters and variables input to curba
c       are real(no double precision)
c     Called only for ieffic=4, else might need to adjust lh=
c       calculation.
c-----------------------------------------------------------------
      subroutine eff_Lin_Liu(z,r,yma,xma,r0x,z_eff,temp,den,jwave,cnpar,
     &                   cnper,ioxm,ye,
     &                   cefldx,cefldy,cefldz,
     &                   efficient)

c     input parameters: z,r -coordinates of the ray point(cm)

c                       yma,zma -coordinares of magnetic axis(cm)
c                       r0x character length (m)
c                       z_eff-effective charge
c                       temp-temperature (keV)
c                       den-density(10**13 cm**-3)
c                       jwave -wave harmonic number
c                       cnpar -parallel to magnetic field refractive
c                              index
c                       cnper Re(N_perp)    
c                       ioxm= +1 O mode,
c                             -1 X mode
c                     	ye-omega_Be/omega_wave
c     cefldx = the x-component of the wave electric field (COMPLEX)
c     cefldy = the y-component                            (COMPLEX)
c     cefldy = the z-component                            (COMPLEX)
c     output parameter: J/P=efficient  in (A/m**2)/(Wt/m**3)
c------------------------------------------------------------------
c     these double precision parameters are to call function: 
c        psif(z,r) and  subroutine: zr_psith(psi,theta,z,r)
c-------------------------------------------------------------------
      implicit none
c-----input
      real*8 z,r,yma,xma,z_eff,temp,den,cnpar,ye,cnper
      real*8 r0x
      complex*16 cefldx,cefldy,cefldz
      integer jwave,ioxm
c-----output
      real*8 efficient

c-----locals
      real*8 psid,zd,rd,rmaxpsid,rminpsid,pid,
     &zbmin,zbmax,zb
      real*8 tc,tol,elomom,prho,ctheta,stheta,pi,theta,thtc,enpar,denom,
     *aspct,zacosarg,rjpd,rjpd0,ratjpd,efficien
      integer ig,n0,model,lh
c-----external 
      real*8 psif,bmin_psi,bmax_psi,rmax_psi,rmin_psi,b,zfac_f

      
c-------------------------------------------------------------------
c     for TorGA_curgap
      real*8 rjpd_d,rjpd0_d,ratjpd_d,denom_d,aspct_d,enpar_d, 
     &tc_d, thtc_d, theta_d, elomom_d,zeff_d,tol_d
c-----external
c-----for efficiency transformation like in currn in toray
      real*8 phi,zfac
c-------------------------------------------------------------------
c      write(*,*)'eff_lin_liu z,r,yma,xma,r0x,z_eff',
c     &z,r,yma,xma,r0x,z_eff
c----------------------------------------------------------------
c     ig     = 1 : the relativistic Green's function
c            = 0 : the non-relativistic approximation
c---------------------------------------------------------------------
      ig=1
c--------------------------------------------------------------------
c     n0 number of points in the Gaussian quadrature (ngauss = 64
c              is recommended.)
c---------------------------------------------------------------------
      n0=64
c--------------------------------------------------------------------
c     tc is a bulk electron temperature in keV
c--------------------------------------------------------------------
      tc=temp
c--------------------------------------------------------------------
c     tol = The relative error tolerence in numerical integration is set
c           to be MAX (tol, 1.0E-6).
c--------------------------------------------------------------------
c     tol=2.e-3
      tol=1.d-2
c--------------------------------------------------------------------
c     model Absolute value of model selectes collisional model: 1 for
c     full bounce av, 2 for square well(numerical solution), 3 for
c     analytic solution to square well.negative model does parallel
c     heating (lower hybrid, fast wave)
c
c     model < 5  gives rjpd in CURGAC
c           = 5  gives rjpd using the exact polarization-dependent
c                rf diffusion operator
c           > 5  gives rjpd using the polarization-dependent rf diffusion
c                operator with but small gyro-radius expansion 
c-------------------------------------------------------------------
      model=3
      model=5
c-------------------------------------------------------------------
c     z_eff is ion effective charge
c-------------------------------------------------------------------
c     lh = the cyclotron harmonic number
c--------------------------------------------------------------------

cBH151018 For ECCD (ieffic=3 or 4) chosing the harmonic for the
cBH151018 efficiency calculation to correspond to harmonic of
cBH151018 maximum EC power absorption.
cBH151018 Formula in accord with TORBEAM (as pointed out by Nicola
cBH151018 Bertelli).  This enables multiple harmonic CD, for example,
cBH151018 in DIII-D, 3rd harmonic near outer edge, 2nd near plasma
cBH151018 center.
cBH151018 Justification for the formula can be seen by considering,
cBH151018 e.g., multiharmonic case along with formulas 3.37 and 
cBH151018 Fig. 3.3 of http://www.compxco.com/cql3d_manual_150122.pdf .


      if (jwave.ne.0) then
         lh=jwave
      else
         !YuP lh=int(sqrt(1d0-cnpar**2)/abs(ye))+1
         !YuP[2020-08-04] BUG: when cnpar>1, it results in 
         ! large negative values for lh, which gives INF 
         ! in other subroutines.  
         ! Adjusted to sqrt(max(0.d0, 1d0-cnpar**2))
         lh= idint( dsqrt(dmax1(0.d0,1d0-cnpar**2))/dabs(ye) ) +1
      endif
c     lh=2
c     lh=0
c--------------------------------------------------------------------
c     elomom is ABS (lh)*cyclotron frec over wave frec for Y<0,
c       interpret as evaluated at poloidal angle where electrons are
c       first trapped in bucket.
c     elomom=yy = lh*omega_c/omega (y in Refs [1] and [2])
c--------------------------------------------------------------------
      elomom=lh*dabs(ye)
c--------------------------------------------------------------------
c     theta  = poloidal angle at which power is absorbed in radians
c              0: outborad; pi: inboard
c--------------------------------------------------------------------
      prho=dsqrt((z-yma)**2+(r-xma)**2)
      ctheta=(r-xma)/prho
      stheta=(z-yma)/prho
      pi=4.d0*datan(1.d0)
c      write(*,*)' effcurba ctheta_geom,stheta',ctheta,stheta
      if(stheta.ge.0.0d0) then
         theta=dacos(ctheta)
      else
         theta=2.0d0*pi-dacos(ctheta)
      end if
c      write(*,*)' effcurba acos(ctheta),theta',acos(ctheta),theta
c--------------------------------------------------------------------
c     thtc is ratio of temperature along characteristic to that of
c          bulk, except for thtc <0, -thtc is energy (in keV) of bucket
c          rise; 
c          an obsolete variable
c--------------------------------------------------------------------
      thtc=1.0d0
c--------------------------------------------------------------------
c     enpar is k_parallel*c/wave frec. Note enpar**2 <1-elomom**2
c       implies no resonant particles; rjpd set to zero and rjpd0
c       determined for vparallel given by nonrelativistic resonance
c       condition
c--------------------------------------------------------------------
      enpar=cnpar
c--------------------------------------------------------------------
c     aspct is inverse aspct ratio
c--------------------------------------------------------------------
c bmax,bmin  in toray
c     if(dabs(bmax-bmin).lt.1.0d-8) bmax=bmin+1.0d-8
c 	aspct=(bmax-bmin)/(bmax+bmin)
c--------------------------------------------------------------------
c     Here :aspct is calculaed using rmax_psi and rmin_psi
c  conversion from real to double precision
      pid=4.d0*datan(1.d0)
c conversion from non-dimensional to m
c      write(*,*)'z,r',z,r
      zd=dble(z)*0.01d0/r0x
      rd=dble(r)*0.01d0/r0x
      psid=psif(zd,rd)
c      write(*,*)'zd,rd,psid',zd,rd,psid
cSm040426
cc     rmax_psid and rmin_psid are the largest and and the least
cc     values of the radius (r) on the given flux surface
cc     psid=psid(zd,rd)
c------------------------------
c      call zr_psith(psid,0.d0,zd,rmaxpsid)
c      call zr_psith(psid,pid,zd,rminpsid)
      
c      if(abs(rmaxpsi-rminpsi).lt.1.0e-8) rminpsi=rmaxpsi-1.0e-8
c      aspct=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
      
c----------------------------------
c     Here :aspct is calculating using double precision functions
c           rmax_psi(psid) and rmin_psi(psid)
c----------------------------------------
c     rmaxpsi and rminpsi are the largest and and the least
c     values of the radius (rmaxpsi and rminpsi in (m))
c     on the given flux surface psi=psi(zd,rd)
c--------------------------------------------
c      rmaxpsid=rmax_psi(psid)
c      rminpsid=rmin_psi(psid)
c      if (rmaxpsid.lt.rd)rmaxpsid=rd
c      if (rminpsid.gt.rd)rminpsid=rd
c      if(dabs(rmaxpsid-rminpsid).lt.1.0d-8) rminpsid=rmaxpsid-1.0d-8
c      aspct=(rmaxpsid-rminpsid)/(rmaxpsid+rminpsid)
cSm040426
       zbmin=bmin_psi(psid)    
       zbmax=bmax_psi(psid)
       aspct=((zbmax-zbmin)/(zbmax+zbmin))
cSm060201
       if (aspct.lt.1.d-8) aspct=1.d-8 

c       write(*,*)'psid,zbmax,zbmin,aspct',psid,zbmax,zbmin,aspct
c      write(*,*)'aspct from functionsrmin_psi, rmax_psi',aspct
c      write(*,*)'aspct',aspct,'enpar',enpar,'tc',tc,'thtc',thtc
c      write(*,*)'theta',theta,'elomom',elomom,'lh',lh
c      write(*,*)'z_eff',z_eff,'model',model,'n0',n0,'ig',ig

cSm040503
c     function b() uses the arguments z and r in [m]
      zb=b(.01d0*dble(z)/r0x,.01d0*dble(r)/r0x,0.d0)!it changes 
                                         !bz,br,bphi and rho in one.i
      zacosarg=(zbmax + zbmin - 2.d0*zbmax*zbmin/zb)/
     &         (zbmin - zbmax)

c      write(*,*)'zbmax,zb,zbmin,zacosarg',zbmax,zb,zbmin,zacosarg
      if (dABS (zacosarg) .gt. 1.0d0) zacosarg = 
     1                               DSIGN (1.0d0, dble(zacosarg))
      theta =dACOS (zacosarg) 
      denom=1.d0
c      write(*,*)'effcurb befre curba denom,aspct,enpar,tc,thtc',
c     &denom,aspct,enpar,tc,thtc
c      write(*,*) 'theta,elomom,lh,z_eff,model,tol,n0,ig',
c     &            theta,elomom,lh,z_eff,model,tol,n0,ig

c-----using  TorGAcurba with real*8 arguments
      denom_d= denom
      aspct_d=aspct
      enpar_d=enpar 
      tc_d=tc
      thtc_d=thtc
      theta_d=theta
      elomom_d=elomom
      zeff_d=z_eff
      tol_d=tol

      !write(*,*)'prep3d/LinLiu befor TorGA_curgap' 
      call TorGA_curgap(rjpd_d,rjpd0_d,ratjpd_d,denom_d,
     &aspct_d,enpar_d,
     &cnper, ioxm,cefldx,cefldy,cefldz,
     &tc_d,thtc_d,theta_d,elomom_d,lh,zeff_d,model,tol_d,n0,ig)

      rjpd=rjpd_d
      rjpd0=rjpd0_d
      ratjpd=ratjpd_d
      !write(*,*)'prep3d/LinLiu after TorGA_curgap' 
      !write(*,*)'rjpd,rjpd0,ratjpd,lh',rjpd,rjpd0,ratjpd,lh

c----------------------------------------------------------
c     efficiency  in (A/cm**2)/(erg/(sec*cm**3))
c     temperature temp in kev
c     density     den in 10**19/m**3
c-----------------------------------------------------------
c-----genray transformation of CD efficiency 
c     cln -Coulomb Ln
c     it is necessary to write the formula to calculate cln

c      cln=17.0
c      arg1=1.d3/temp*dsqrt(10.d0*den) !temp KeV, den 10**13/cm**3
c      cln=24.d0-alog(arg1)
       efficien=rjpd 
c      efficient=efficien*temp/den*(17.0/cln)*4.5*1.e-6
c      write(*,*)'efficiecy coefficient=temp/den*(17.0/cln)*4.5*1.e-6',
c     1temp/den*(17.0/cln)*4.5*1.e-6
c      write(*,*)'in effic rjpd,rjpd0,ratjpd,denom',
c     1           rjpd,rjpd0,ratjpd,denom
c      write(*,*)'end effic curba efficient',efficient

c-----the efficiency transformation from torayGA
c      write(*,*)'before zfac_f temp',temp
      zfac=zfac_f(dble(z*1.d-2),dble(r*1.d-2),0.d0,dble(temp))
      zfac=zfac*1.d-7 ! watt => egr/sec
      zfac=zfac*1.d-4 ! m**2 => cm**2
c      write(*,*)'zfac',zfac
      efficient=efficien*zfac
cSm080629
c      efficient=-efficient 

      return
      end subroutine eff_Lin_Liu



       subroutine etajrf(z,eps,w,rto,eta,eta0,irfq) 
! +-------------------------------------------------------------------+
! | Last modified: July 3,1997                                        |
! | Written by P. Bonoli & J. Wright                                  |
! +-------------------------------------------------------------------+
! |  Evaluates the normalized current-drive efficiency using          |
! |  the parametrization of Ehst and Karney                           |
! |  (Argonne National Laboratory Report ANL/FPP/TM-247 (1990)).      |
! |  IRFQ = 1 indicates alfven wave type damping                      |
! |  IRFQ = 2 indicates landau wave type damping                      |
! +-------------------------------------------------------------------+
!   w=wte
!   rto=bmod/bmax
!   irfq=jwave+2 = For Alfven +1
!                = for Ladaw damp +2
      implicit none

      integer, intent(in) :: irfq

      real*8, intent(in) :: z, eps, w, rto
      real*8, intent(out) :: eta, eta0
      real*8 :: xr, en, ek, a, fk, fc, ec, em, r, flamt, fm,
     &            yt, ytt, gc, rtp

c      write(*,*)'in etajrf z,eps,w,rto,irfq',z,eps,w,rto,irfq 

      xr = 3.5d0
      en = 0.77d0
      ek = 3.d0

      if(irfq.eq.1)  then
         a = 12.3d0
         fk = 11.91d0/(0.678d0 + z)
         fc = 4.13d0/z**0.707d0
         ec = 0.0987d0
         em = 2.48d0
      endif

      if(irfq .eq. 2) then
         a = 0.d0
         fk = 3.d0/z
         fc = 3.83d0/z**0.707d0
         ec = 0.389d0
         em = 1.38d0
      endif

c      write(*,*)'in etajrf a,fk,fc,ec ,em',a,fk,fc,ec ,em

      eta0 = (fk/w + fc + 4.d0*w*w/(5.d0 + z))
      
c      write(*,*)'in etajrf eta0',eta0
     
      r = 1.d0 - eps**en*sqrt(xr*xr + w*w)/(eps**en*xr + w)
c      R=r=rprof

       rtp = min(rto,0.9999999d0)
       flamt = sqrt(1.d0 - rtp)
c------flamt=alt1=sqrt(1-bmod/bmx) => rtp=bmod/bmax 
       fm = 1.d0 + a*(flamt/w)**ek

       yt = (1.d0 - flamt*flamt)/(flamt*flamt)*w*w
c-----      Y_t=yt=ytt
c-----      flamt=alt1
       ytt = (ec*yt)**em
       ytt = min(500.d0,ytt)
       gc = 1.d0 - exp(-ytt)
       eta = gc*fm*eta0*r !normalized eta

       return
       end subroutine etajrf


       Subroutine efKarney_Bonoli(z,r,r0x,z_eff,temp,den,jwave,cnpar,
     & effKarn)
c-------RF current drive efficiency(asymptotic formula, nonrelativistic case)  
c  D.A. Ehst and C.F.F.Karney Nucl.Fus. Vol.31, No.10 (1991), p 1933-1938
c  A subroutine to calculate the local current density driven by RF waves
c--Uses CD efficiency empirical formula based on numerical Fokker-
c  Planck bounce-averaged calculations
c  This formula is for
c  1)Landau damping of lower hybrid (LH) slow waves resonant  at parallel
c    velocities	above the electron thermal velocity
c  2)Slow frequency fast (compressional Alfven) wave (AW) may resonant with
c    low phase velocity electrons via combined Landay damping and transit
c    time magnetic damping  
c
c It uses Bonoli subroutine etajrf
c     input parameters: z,r -coordinates of the ray point(cm!!!!!)ATTENTION
c                       z_eff-effective charge
c                       temp-temperature (keV)
c                       den-density(10**13 cm**-3)
c                       jwave (=islofa)=-1 ALfen wave , 0- Landau damp.
c                       cnpar -paralell to magnetic field refractive
c                              index
c                       r0x character length
c     computed:  J/P=efficiency  in (A/m**2)/(W/m**3)
c     converted to 
c     output parameter: effKarn in (A/cm**2)/(erg/(sec*cm**3))
c------------------------------------------------------------------
c     It uses:
c     double precision function: psif(z,r) and
c     subroutine: zr_psith(psi,theta,z,r)
c     double precision functions from zr_psith.f: rmax_psi(psi),
c     rmin_psi(psi), bmax_psi(psi)

      implicit none
c-----input
      real*8 z,r, !cm
     &r0x,        !normalization lemgth
     &z_eff,      !effective charge
     &temp,       !temperatura KeV
     &den,        !density 10**13
     &cnpar       !parallel to magnetic field refractive index

      integer jwave 
c-----output
      real*8 ffKarn    !efficiency  in (A/cm**2)/(erg/(c*cm**3))


c-----locals
      real*8
     &w,           !  w=wte=u_res=v_par/v_e=cvac/(ve*cnpar)
     &cvac,        !  light speed [cm/sec]
     &ve,          !  ve=sqrt(T_e/m_e)
     &zd,rd,       !  coordinates in [m]
     &phi,psid,rmaxpsi,rminpsi,aspct,eps,bmod,bmax,
     &rto,         !  b/bmax
     &eta0,eta,    !  normalized efficiency from Bonoli subrourine
     &arg1,cln,s
      integer irfq !=jwave+2
       
c-----outpt
      real*8 effKarn !efficiency  in (A/cm**2)/(erg/(sec*cm**3)

c-----externals
      real*8 psif,b,bmax_psi,rmin_psi,rmax_psi

      cvac=2.997930d+10
      ve=1.32d9*dsqrt(temp)	 
      w=cvac/(ve*cnpar)
      w=dabs(w)     !w=|wte|
      write(*,*)'Bon w',w
c-------------------------------------------
c     conversion from cm to m
      zd=z*0.01d0/r0x
      rd=r*0.01d0/r0x
      phi=0.d0       !arbitrary vatoroidal angle
      write(*,*)'Bon z,r,cnpar',z,r,cnpar 
c----------------------------------------------
c     calculate epsilon
c-------------------------------------------
      psid=psif(zd,rd) !poloidal flux at z,r
      rmaxpsi=rmax_psi(psid)
      rminpsi=rmin_psi(psid)
      if (rmaxpsi.lt.rd)rmaxpsi=rd
      if (rminpsi.gt.rd)rminpsi=rd 
      if(dabs(rmaxpsi-rminpsi).lt.1.0d-8) rminpsi=rmaxpsi-1.0d-8
      aspct=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
      eps=aspct 
      write(*,*)'Bon eps',eps
c------------------------------------------------
c     calculate rto=b/bmax
c--------------------------------------------------
      bmod=b(zd,rd,phi)
      bmax=bmax_psi(psid)
      if (bmod.gt.bmax) bmod=bmax 
      rto=bmod/bmax
      write(*,*)'Bon rto',rto
c--------------------------------------------------
      irfq=jwave+2
      write(*,*)'Bon irfq',irfq
c----------------------------------------------------
c     calculate normalized efficiency eta 
c---------------------------------------------------
      write(*,*)'Bon before etajrf z_eff',z_eff,eps,w,rto,irfq
      call etajrf(z_eff,eps,w,rto,eta,eta0,irfq)
      write(*,*)'Bon after etajr rta0,eta',eta0,eta
c-------------------------------------------------- 
c     efficiency  in (A/m**2)/(W/m**3)
c     temperature temp in kev
c     density     dens in 10**19 /m**3
c-----------------------------------------------------------
c     cln -Coulomb Ln
      arg1 = 1.d3/temp*dsqrt(10.d0*den)	!temp KeV, den 10**13/cm**3
      cln=24.d0-dlog(arg1)

      effKarn=eta*3.84d0*temp/(cln*den)	!(A/m**2)/(W/m**3)
      write(*,*)'Bon temp,den,cln,eta,effKarn',temp,den,cln,eta,effKarn
c--------------------------------------------------------------  
c     determination the of the current drive sign
c------------------------------------------------------
      if (cnpar.ge.0.d0) then
         s=1.d0
      else
         s=-1.d0
      endif
c------------------------------------------------------------
      effKarn=-effKarn*s*1.d-5  !(A/cm**2)/(erg/(sec*cm**3)) converted

      return
      end subroutine efKarney_Bonoli ! Not used; for tests only

      subroutine plot_1_Karney
c----------------------------------------------------
c     creates data (arrays) for plot like Fig. 1
c     at Karney article
c     Nuclear Fusion 1991 p. 1934
c--------------------------------------------------- 
      implicit none
      include 'param.i'
      include 'three.i'

      integer n_points_w,  !number of plot points in w direction
     &n_epsilon,  !number of epsilons
     &n_theta_pol !number of polidal angles
      parameter (n_points_w=11)
c      parameter (n_points_w=3)
      parameter (n_epsilon=3)
      parameter (n_theta_pol=4)

      integer i,j,k

      real*8 
     &epsilon,        !inverse aspect ratio
     &epsilon_ar(n_epsilon),
     &theta_pol,      !poloidal angle [radians]
     &theta_pol_ar(n_theta_pol),
     &w_ar(n_points_w),w,    !=clight/(N_parallel*V_te)
                      !here m_e*V_te**2=T_e
     &w_min,w_max,step,
     &pi,
     &psi,            !poloidal flux(epsilon)
     &z_eff,
     &eta,eta0,       !CD efficiency
     &eta_ar_fw(n_points_w,n_theta_pol,n_epsilon),
     &eta0_ar_fw(n_points_w,n_theta_pol,n_epsilon),
     &eta_ar_lh(n_points_w,n_theta_pol,n_epsilon),
     &eta0_ar_lh(n_points_w,n_theta_pol,n_epsilon),
     &z,r,phi,        !space coordinates
     &rto,            !rto=b(z,r,phi)/bmax
     &bmod,bmax,
     &accuracy       !accuracy to solve the equation
                     ! epsilon_psi(psi)=epsilon
      real*8 rmaxpsi,rminpsi,rmax_psi,rmin_psi

      integer irfq !=1 indicates alfven wave type damping                      
                    != 2 indicates landau wave type damping       

      character*16 file_nm ! name of output netcdf file
                   
c-----externals
      real*8 psi_epsilon,bmax_psi,b,rhopsi


      pi=4.d0*datan(1.d0)

      theta_pol_ar(1)=0.d0
      theta_pol_ar(2)=pi/2.d0
      theta_pol_ar(3)=3.d0*pi/4.d0
      theta_pol_ar(4)=pi

      epsilon_ar(1)=0.d0
      epsilon_ar(2)=0.03
      epsilon_ar(3)=0.1d0

      z_eff=1.d0
      irfq=1
      phi=0.d0
      w_min=0.1d0
      w_max=10.d0
c      step=(w_max-w_min)/(n_points_w-1)
      step=(dlog(w_max/w_min)/dlog(10.d0))/(n_points_w-1 )    
      write(*,*)'in prep3d subroutine plot_1_Karney'
      accuracy=1.d-7
      do i=1,n_epsilon
         epsilon=0.d0
         epsilon=epsilon_ar(i)
         write(*,*)'i,epsilon',i,epsilon

         psi=psi_epsilon(epsilon,psimag,psilim, accuracy)
         write(*,*)'after psi_epsilon psi=',psi
         write(*,*)'rhopsi(psi)',rhopsi(psi)
c--------check epsilon
         rmaxpsi=rmax_psi(psi)
         rminpsi=rmin_psi(psi)
     
         if(dabs(rmaxpsi-rminpsi).lt.1.0d-8) rminpsi=rmaxpsi-1.0d-8
         epsilon=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
         write(*,*)'check epsilon',epsilon

         do j=1,n_theta_pol
            theta_pol=theta_pol_ar(j)
            call zr_psith(psi,theta_pol,z,r)
            bmod=b(z,r,phi)
            bmax=bmax_psi(psi)
            if (bmod.gt.bmax) bmod=bmax 
            rto=bmod/bmax
            write(*,*)'i,j,rto',i,j,rto
            do k=1,n_points_w
               w_ar(k)=w_min+step*(k-1)
c               x=dlog(w_min)/dlog(10.d0)+step*(k-1)
               w_ar(k)=dexp(dlog(10.d0)*
     &                      (dlog(w_min)/dlog(10.d0)+step*(k-1)))
               w=w_ar(k)
               write(*,*)'k, w_min,w',k, w_min,w
               irfq=1 !fw
               call etajrf(z_eff,epsilon,w,rto,eta,eta0,irfq)
               eta_ar_fw(k,j,i)=eta 
               eta0_ar_fw(k,j,i)=eta0 
               write(*,*)'fw eta0,eta',eta0,eta
               irfq=2 !lh
               call etajrf(z_eff,epsilon,w,rto,eta,eta0,irfq)
               eta_ar_lh(k,j,i)=eta 
               eta0_ar_lh(k,j,i)=eta0
               write(*,*)'lh eta0,eta',eta0,eta
            enddo
        enddo
      enddo

c--------------------------------------------------------------
c     write CD efficiency eta and w for fig. 1 to netcdf file
c------------------------------------------------------------
      file_nm='Ehst-Karney_plot' ! set name of output netcdf file

      write(*,*)'before netcdf_Karne !sqrt(T/m)y file_nm=',file_nm

cSAP081110
c      if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
c        file_nm='Ehst-Karney_plot' ! set name of output netcdf file
c        write(*,*)'before netcdf_Karne !sqrt(T/m)y file_nm=',file_nm
c      call netcdf_Karney(file_nm,n_points_w,w_ar,
c     &n_theta_pol,theta_pol_ar,n_epsilon,epsilon_ar,
c     &eta0_ar_fw,eta_ar_fw,eta0_ar_lh,eta_ar_lh)
c        write(*,*)'after netcdf_Karney '
c      endif


      return
      end

      real*8 function epsilon_psi(psi)
c-------------------------------------------------------------
c     calculates inverse aspect ratio 
c     epsilon=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
c     at given ploidal flux psi
c---------------------------------------------------------------
      implicit none
     
c-----input
      real*8 psi ! poloidal flux

c-----externals rmax_psi,rmin_psi
      real*8 rmax_psi,rmin_psi


c-----locals
      real*8 rmaxpsi,rminpsi

      rmaxpsi=rmax_psi(psi)
      rminpsi=rmin_psi(psi)

      if(dabs(rmaxpsi-rminpsi).lt.1.0d-8) rminpsi=rmaxpsi-1.0d-8

      epsilon_psi=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)

      return
      end

      real*8 function epsilon_psi_minus_epsilon(psi)
c------------------------------------------------------------------
c     Calculates (epsilon_psi(psi)-epsilon)
c     It used to find the root of the equation
c     epsilon_psi(psi)-epsilon=0
c-----------------------------------------------------------------
      implicit none
c-----input
      real*8 psi !poloidal flux

      real*8 epsilon_in
      common /epsilon_psi_minos_epsilon/ epsilon_in     
c-----externals
      real*8 epsilon_psi

      epsilon_psi_minus_epsilon=epsilon_psi(psi)-epsilon_in
       
      return
      end

      real*8 function psi_epsilon(epsilon,psimag,psilim, accuracy)
c-----------------------------------------------------------
c     calculates poloidal flux at given inverse aspect ratio epsioln
c     epsilon=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
c
c     It uses bisection method at psimag < psi < psilim
c-------------------------------------------------------------
      implicit none
     
c-----input
      real*8 epsilon, !inverse aspect ratio
     &psimag,         !polidal flux at the magnetic axis 
     &psilim,         !poloidal flux at the last closed flux surface
     &accuracy        !bisection method accuracy 
                      !It uses the condition |difference of root| < accuracy
      real*8 epsilon_in
      common /epsilon_psi_minos_epsilon/ epsilon_in ! will be set equal to
                                                    ! the input epsilon
c-----externals
      real*8  epsilon_psi_minus_epsilon,rtbis
      external  epsilon_psi_minus_epsilon

c-----initialization of the common block
      epsilon_in =epsilon

      if(epsilon.lt.1.d-10) then
        psi_epsilon=psimag
      else

        write(*,*)'in psi_epsilon psimag,psilim,accuracy',
     &                            psimag,psilim,accuracy

        psi_epsilon=RTBIS(epsilon_psi_minus_epsilon,
     &psimag,psilim,accuracy)

      endif

      return
      end
     
      subroutine netcdf_Karney(file_nm,n_points_w,w_ar,
     &n_theta_pol,theta_pol_ar,n_epsilon,epsilon_ar,
     &eta0_ar_fw,eta_ar_fw,eta0_ar_lh,eta_ar_lh)
c-------------------------------------------------------------
c     writes data for fig1 Ehst-Karney  to the netcdf file
c     file_nm.nc
c--------------------------------------------------------------
      implicit none
      include 'netcdf.inc'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

c-----input
      character*(*) file_nm        !name of output nc file
      integer n_points_w,          !number of points
     &n_theta_pol,                 !number of poloidal angle values
     &n_epsilon                    !number of epsilon values
      real*8,  dimension(1:n_points_w)  :: w_ar         ! c/(N_par*v_te)
      real*8,  dimension(1:n_theta_pol) :: theta_pol_ar !ploidal angle [radian]
      real*8,  dimension(1:n_epsilon)   :: epsilon_ar   !inverse aspect ratio
      real*8,  dimension(1:n_points_w, 1:n_theta_pol, 1:n_epsilon) 
     &                                  :: eta_ar_fw,eta0_ar_fw, !FW efficiency
     &                                     eta_ar_lh,eta0_ar_lh  !LH efficiency
c-----locals----------------------------------------  
      integer i,j,k
c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,n_points_w_id,n_theta_pol_id,
     &n_epsilon_id,
     &start(3),starts(3),eta_dims(3),eta_count(3),
     +nccre2,ncvdef2,ncdid2,ncddef2

      character ltitle*256,filenc*128

c-----Storage tem1 is used in netcdf write.   
      real*8, dimension (1: n_points_w*n_theta_pol) :: tem1
      real*8, dimension (1: n_points_w,1:n_theta_pol) :: tem2
c-----externals
      integer length_char

      data start/1,1,1/

      if(myrank.ne.0) return
     
      write(*,*)'in netcdf_Karney begin'

      eta_count(1)=n_points_w
      eta_count(2)=n_theta_pol
      eta_count(3)=1

      write(*,*)'in netcdf_Karney n_points_w,n_theta_pol,n_epsilon',
     &                            n_points_w,n_theta_pol,n_epsilon
      write(*,*)'in netcdf_Karney file_nm=',trim(file_nm)

      if( length_char(file_nm).gt.124)
     &   stop 'Adjust file_nm in netcdf_Karney'

      write(filenc,1002) trim(file_nm)
 1002 format(a,".nc")
      write(*,*)'in netcdf_Karney after length_char(file_nm)'
c-------------------------------------------------------------
c      create net CDF file  define dimensions,variables
c          and attributes
c------------------------------------------------------------
      ncid=nccre2(trim(filenc),NCCLOB,istatus)
      call check_err(istatus)

c     Brief description added to file:
      ltitle='netCDF file of data for fig.1 WEhst-Karney article'
      if( length_char(ltitle).gt.256 ) 
     &   stop 'Adjust ltitle in netcdf_Karney'
      call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,length_char(ltitle),
     &     ltitle,istatus)

c.......................................................................
cl    1.1.2 define dimensions
c.......................................................................
      n_epsilon_id=ncddef2(ncid,'n_epsilon',n_epsilon,istatus)
      n_theta_pol_id=ncddef2(ncid,'n_theta_pol',n_theta_pol,istatus)
      n_points_w_id=ncddef2(ncid,'n_points_w',n_points_w,istatus)
     
c.......................................................................
cl    1.1.3 define variables
c.......................................................................
      vid=ncvdef2(ncid,'n_points_w',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27,
     +            'Number of argument points_w',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'n_epsilon',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of epsilons',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'n_theta_pol',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +            'Number of poloidal angles',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_ar',NCDOUBLE,1,n_points_w_id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +     'array of phase velocities w=c/(N_par*v_e)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'theta_pol_ar',NCDOUBLE,1,n_theta_pol_id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +     'array of poloidal angles',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'radians',istatus)
      call check_err(istatus)  

      vid=ncvdef2(ncid,'epsilon_ar',NCDOUBLE,1,n_epsilon_id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,16,
     +     'array of epsilon',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      eta_dims(1)=n_points_w_id
      eta_dims(2)=n_theta_pol_id
      eta_dims(3)=n_epsilon_id

      write(*,*)'ncid=',ncid
      write(*,*)'eta_dims',eta_dims
      vid=ncvdef2(ncid,'eta_ar_fw',NCDOUBLE,3,eta_dims,istatus)
      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'FW CD efficiency eta',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10,
     +           'normalized',istatus)

      write(*,*)'ncid=',ncid
      write(*,*)'eta_dims',eta_dims
      vid=ncvdef2(ncid,'eta0_ar_fw',NCDOUBLE,3,eta_dims,istatus)
      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +           'FW CD efficiency eta0',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10,
     +           'normalized',istatus)

      write(*,*)'ncid=',ncid
      write(*,*)'eta_dims',eta_dims
      vid=ncvdef2(ncid,'eta_ar_lh',NCDOUBLE,3,eta_dims,istatus)
      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'LH CD efficiency eta',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10,
     +           'normalized',istatus)

      write(*,*)'ncid=',ncid
      write(*,*)'eta_dims',eta_dims
      vid=ncvdef2(ncid,'eta0_ar_lh',NCDOUBLE,3,eta_dims,istatus)
      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +           'LH CD efficiency eta0',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10,
     +           'normalized',istatus)

c------------------------------------------------------------
cl    1.1.4 end the define-mode and start the data-mode
c     p. 51-2 of manual
c----------------------------------------------------------
      call ncendf2(ncid,istatus)
      call check_err(istatus)

      write(*,*)'end initialization'
c----------------------------------------------------------------
c     write data
c-----------------------------------------------------------------
      write(*,*)'before vid=ncvid(ncid,n_epsilon'
      call ncvid2(vid,ncid,'n_epsilon',istatus)
      call ncvpt_int2(ncid,vid,1,1,n_epsilon,istatus)

      call ncvid2(vid,ncid,'n_theta_pol',istatus)
      call ncvpt_int2(ncid,vid,1,1,n_theta_pol,istatus)

      call ncvid2(vid,ncid,'n_points_w',istatus)
      call ncvpt_int2(ncid,vid,1,1,n_points_w,istatus)

      call ncvid2(vid,ncid,'epsilon_ar',istatus)  
      call ncvpt_doubl2(ncid,vid,1,n_epsilon,epsilon_ar,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'theta_pol_ar',istatus)  
      call ncvpt_doubl2(ncid,vid,1,n_theta_pol,theta_pol_ar,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_ar',istatus)  
      call ncvpt_doubl2(ncid,vid,1,n_points_w,w_ar,istatus)
      call check_err(istatus)
 
c-----FW    
      do i=1,n_epsilon
c        write(*,*)'eta i=',i
         do j=1,n_theta_pol
c            write(*,*)'eta j',j
            do k=1,n_points_w
               tem2(k,j)=eta_ar_fw(k,j,i) 
c               write(*,*)'k,eta_ar_fw(k,j,i)',k,eta_ar_fw(k,j,i) 
            enddo
         enddo
          
         call pack21(tem2,1,n_points_w,1,n_theta_pol,
     &               tem1,n_points_w,n_theta_pol)

         starts(1)=start(1)
         starts(2)=start(2)
         starts(3)=i  ! number of epsilon put in .nc file

         call ncvid2(vid,ncid,'eta_ar_fw',istatus)
         call ncvpt_doubl2(ncid,vid,starts,eta_count,tem1,istatus)
      enddo

      do i=1,n_epsilon
         do j=1,n_theta_pol
            do k=1,n_points_w
               tem2(k,j)=eta0_ar_fw(k,j,i)  
            enddo
         enddo
          
         call pack21(tem2,1,n_points_w,1,n_theta_pol,
     &               tem1,n_points_w,n_theta_pol)

         starts(1)=start(1)
         starts(2)=start(2)
         starts(3)=i  ! number of epsilon put in .nc file

         call ncvid2(vid,ncid,'eta0_ar_fw',istatus)
         call ncvpt_doubl2(ncid,vid,starts,eta_count,tem1,istatus)
      enddo

c-----LH
      do i=1,n_epsilon
c        write(*,*)'eta i=',i
         do j=1,n_theta_pol
c            write(*,*)'eta j',j
            do k=1,n_points_w             
                tem2(k,j)=eta_ar_lh(k,j,i)
c               write(*,*)'k,eta_ar_lh(k,j,i)',k,eta_ar_lh(k,j,i) 
            enddo
         enddo
          
         call pack21(tem2,1,n_points_w,1,n_theta_pol,
     &               tem1,n_points_w,n_theta_pol)

         starts(1)=start(1)
         starts(2)=start(2)
         starts(3)=i  ! number of epsilon put in .nc file

         call ncvid2(vid,ncid,'eta_ar_lh',istatus)
         call ncvpt_doubl2(ncid,vid,starts,eta_count,tem1,istatus)
      enddo

      do i=1,n_epsilon
         do j=1,n_theta_pol
            do k=1,n_points_w
               tem2(k,j)=eta0_ar_lh(k,j,i)  
            enddo
         enddo
          
         call pack21(tem2,1,n_points_w,1,n_theta_pol,
     &               tem1,n_points_w,n_theta_pol)

         starts(1)=start(1)
         starts(2)=start(2)
         starts(3)=i  ! number of epsilon put in .nc file

         call ncvid2(vid,ncid,'eta0_ar_lh',istatus)
         call ncvpt_doubl2(ncid,vid,starts,eta_count,tem1,istatus)
      enddo
C-----------------------------------------------------------------------
c
cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)
      return
      end

      
      subroutine ADJ_eff_at_w(z,r,phi,theta_pol,w_t,cd_efficiency)  
c-----calculate CD efficiency cd_efficiency using ADJ functuion
c     at the give space point  (z,r,phi,theta_pol)
c     and given w_t=c/(N_parallel*v_te)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'cefield.i'
cfor test
      include 'three.i'
c-----input
      real*8 
     &z,r,phi,                    ! space coordinates
     &w_t,                        ! c/(N_parallel*v_te),   v_te=sqrt(T_e/m_e)
     &theta_pol                   ! poloidal angle [radian]
   
c-----output
      real*8  cd_efficiency
c-----locals
      real*8  v_te,clight,cnpar,cnper,ex,ey,ez,eplus,eminus,epar,
     & rho_loc,psi_loc,rme,
     & T_e_kev,                   ! electron temperature [KeV]
     & dense,                     !electron density in 10**13/cm**2
     & arg1,cln
      complex*16 cnx,             ! complex N_perp
     & cex_loc,cey_loc,cez_loc    ! electric field polarization

      integer id_loc              ! choose the cold plasma dispesrsion relation
                                  ! id_loc=1,2 or 3
      integer iraystop
c-----externals
      real*8 b,temperho,psi_rho,densrho

c-----calculate small radius rho inside b. rho will be in common /one/
      bmod=b(z,r,phi)

c-----calculate electron temperaure T_e_kev
      rho_loc=rho
      T_e_kev=temperho(rho_loc,1)
      dense=densrho(rho_loc,1)

      write(*,*)'ADJ_eff_at_w rho,T_e_kev,dense,w_t',
     &rho,T_e_kev,dense,w_t

      rme=9.1094d-28                    ! electron mass [g]
      v_te=4.19d7*dsqrt(T_e_kev*1.d3)   ! sqrt(T_e/m_e) thermal velocity [cm/sec]
      write(*,*)'1 v_te',v_te
      v_te=dsqrt(T_e_kev*1.6022d-9/rme) !(sqrt(T_e/m_e) thermal velocity [cm/sec]
      write(*,*)'2 v_te',v_te

      clight=2.99792458d10              ! light speed [cm.sec]
     
      cnpar=clight/(v_te*w_t)             ! N_parallel
c--------------------------------------------------------------
c     It solves the cold plasma id=1,2 or
c     Appleton-Harty id=3 dispersion relation N_per=N_per(n_par)
c     for given ioxm_n_npar                              
c     Then subroutine calculates the initial components  
c     of the refractive index  cnz,cnr,cm                
c     ******************************************************
      write(*,*)'cnpar',cnpar
      id_loc=id
      call  nper_npar_ioxm_n_npar(id_loc,z,r,phi,cnpar,
     &cnper,iraystop)
      write(*,*)'iraystop,cnper',iraystop,cnper

c-----calculate  the cold plasma dielectric tensor to 'eps.i'
      call tensrcld(z,r,phi) 
c-----calculate electric field (cex,cey,cez) using the cold plasma dielectric
c     tensor to 'cefield.i'
c
      cnx=dcmplx(cnper,0.d0)
      call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)

c-------------------------------------------------------------------------
c     calculate CD efficiency using adj chi function
c     cd_efficiency is in (A/cm**2)/(erg/(sec*cm**3))
c-------------------------------------------------------------------------
       psi_loc=psi_rho(rho)
       write(*,*)'psimag,psi_loc,psilim',psimag,psi_loc,psilim
       cex_loc=cex
       cey_loc=cey
       cez_loc=cez

       write(*,*)'cex_loc,cey_loc,cez_loc',cex_loc,cey_loc,cez_loc

      call CD_adj_efficiency(cnpar,cnper,
     &psi_loc,theta_pol,cex_loc,cey_loc,cez_loc,
     &cd_efficiency)
      write(*,*)'1 cd_efficiency',cd_efficiency
c--------------------------------------------------------------------------
c     normalized efficiency
c----------------------------------------------------------------------------
c-----cln -Coulomb Ln
      arg1 = 1.d3/t_e_kev*dsqrt(10.d0*dense)	!temp KeV, dense 10**13/cm**3
      cln=24.d0-dlog(arg1)
      
c---------------------------------------------------------------------------
      cd_efficiency=cd_efficiency*1.d5         ! (A/m**2)/(W/m**3)
 
      cd_efficiency=cd_efficiency*(cln*dense)/(3.84d0*t_e_kev)

      write(*,*)'2 cd_efficiency',cd_efficiency
      return
      end

     
      subroutine plot_1_ADJ_Karney
c----------------------------------------------------
c     creates data (arrays) for plot like Fig. 1
c     at Karney article Nuclear Fusion 1991 p. 1934
c     using ADJ efficiency
c--------------------------------------------------- 
      implicit none
      include 'param.i'
      include 'three.i'

      integer n_points_w,  !number of plot points in w direction
     &n_epsilon,  !number of epsilons
     &n_theta_pol !number of polidal angles
      parameter (n_points_w=11)
c      parameter (n_points_w=3)
      parameter (n_epsilon=3)
      parameter (n_theta_pol=4)

      integer i,j,k

      real*8 
     &epsilon,        !inverse aspect ratio
     &epsilon_ar(n_epsilon),
     &theta_pol,      !poloidal angle [radians]
     &theta_pol_ar(n_theta_pol),
     &w_ar(n_points_w),w,    !=clight/(N_parallel*V_te)
                      !here m_e*V_te**2=T_e
     &w_min,w_max,step,
     &pi,
     &psi,            !poloidal flux(epsilon)
     &z_eff,
     &eta,eta0,       !CD efficiency
     &eta_ar(n_points_w,n_theta_pol,n_epsilon),
     &eta0_ar(n_points_w,n_theta_pol,n_epsilon),
     &z,r,phi,        !space coordinates
     &rto,            !rto=b(z,r,phi)/bmax
     &bmod,bmax,
     &accuracy       !accuracy to solve the equation
                     ! epsilon_psi(psi)=epsilon
      real*8 rmaxpsi,rminpsi,rmax_psi,rmin_psi

      integer irfq !=1 indicates alfven wave type damping                      
                    != 2 indicates landau wave type damping       

      character*16 file_nm ! name of output netcdf file
                   
c-----externals
      real*8 psi_epsilon,bmax_psi,b,rhopsi


      pi=4.d0*datan(1.d0)

      theta_pol_ar(1)=0.d0
      theta_pol_ar(2)=pi/2.d0
      theta_pol_ar(3)=3.d0*pi/4.d0
      theta_pol_ar(4)=pi

      epsilon_ar(1)=0.d0
      epsilon_ar(2)=0.03
      epsilon_ar(3)=0.1d0

      z_eff=1.d0
      irfq=1
      phi=0.d0
      w_min=0.1d0
      w_max=10.d0
c      step=(w_max-w_min)/(n_points_w-1)
      step=(dlog(w_max/w_min)/dlog(10.d0))/(n_points_w-1 )    
      write(*,*)'in prep3d subroutine plot_1_Karney'
      accuracy=1.d-7
      do i=1,n_epsilon
         epsilon=0.d0
         epsilon=epsilon_ar(i)
         write(*,*)'i,epsilon',i,epsilon

         psi=psi_epsilon(epsilon,psimag,psilim, accuracy)
         write(*,*)'after psi_epsilon psi=',psi
         write(*,*)'rhopsi(psi)',rhopsi(psi)
c--------check epsilon
         rmaxpsi=rmax_psi(psi)
         rminpsi=rmin_psi(psi)
     
         if(dabs(rmaxpsi-rminpsi).lt.1.0d-8) rminpsi=rmaxpsi-1.0d-8
         epsilon=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
         write(*,*)'check epsilon',epsilon

         do j=1,n_theta_pol
            theta_pol=theta_pol_ar(j)
            call zr_psith(psi,theta_pol,z,r)
            bmod=b(z,r,phi)
            bmax=bmax_psi(psi)
            if (bmod.gt.bmax) bmod=bmax 
            rto=bmod/bmax
            write(*,*)'i,j,rto',i,j,rto
            do k=1,n_points_w
               w_ar(k)=w_min+step*(k-1)
c               x=dlog(w_min)/dlog(10.d0)+step*(k-1)
               w_ar(k)=dexp(dlog(10.d0)*
     &                      (dlog(w_min)/dlog(10.d0)+step*(k-1)))
               w=w_ar(k)
               write(*,*)'k, w_min,w',k, w_min,w
c               call etajrf(z_eff,epsilon,w,rto,eta,eta0,irfq)           
               call ADJ_eff_at_w(z,r,phi,theta_pol,w,eta)  
               eta0=0.d0
               eta_ar(k,j,i)=eta 
               eta0_ar(k,j,i)=eta0 
               write(*,*)'eta0,eta',eta0,eta
            enddo
        enddo
      enddo

c--------------------------------------------------------------
c     write CD efficiency eta and w for fig. 1 to netcdf file
c------------------------------------------------------------
c      if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
c        file_nm='Karney_ADJ_plot1' ! name of output netcdf file
c        write(*,*)'before netcdf_Karne !sqrt(T/m)y file_nm=',file_nm
cSAP081110
c      call netcdf_Karney(file_nm,n_points_w,w_ar,
c     &n_theta_pol,theta_pol_ar,n_epsilon,epsilon_ar,
c     &eta0_ar,eta_ar)
c        write(*,*)'after netcdf_Karney '
c      endif

      return
      end

 

 
      subroutine grpde2 (npar, nper, omega, omegpe, nsigma,
     .                   vgrpdc, edenfac, exde, eyde, ezde)
      implicit none
 
c-----input
      integer nsigma
      real*8 npar, nper, omega, omegpe
c-----output
      real*8  vgrpdc, edenfac
      complex*16  exde, eyde, ezde

c-----locals
      complex*16       rootm1,d11,d12,d13,d22,d33,denom

      real*8           npar2,nper2,n2,nperol
     

****  double precision deps1,deps2,deps3,domega,domegpe,
****  double precision dquad
      real*8             dquad
      real*8 deps1,deps2,deps3,domega,domegpe,aquad,bquad,cquad,
     &tens11,tens12,tens21,tens22,tens33,
     &bmod2,eps1,eps2,eps3,dddnpp,dddnpl,dddome,
     &dep1do,dep2do,dep3do,rootqd,emod2                    
           
c     calculate group velocity, wave polarizations, and energy density factor,
c     using cold plasma relations. (M O'Brien's version)
c
      nperol = nper
      rootm1 = (0.0d0, 1.0d0)
c
c     DOUBLE PRECISION is used since a problem was encountered in the
c                 real is used since a problem was encountered in the
c     following subtraction within the formation of rootqd.
c
c --- double precision mentioned above was disabled 21 Aug 94 by Joe Freeman
c --- not needed since the WHOLE PROGRAM is done with 64-bit real arithmetic
c --- on CRAY this is automatic, elsewhere a compilation switch ensures this
c
      domega  = omega
      domegpe = omegpe
c      write(*,*)'in grpde2 npar,domega,domegpe,nsigma',
c     1                     npar,domega,domegpe,nsigma
      deps1   = 1.0d0 - domegpe*domegpe/        (domega*domega-1.0d0)
      deps2   =      -domegpe*domegpe/(domega*(domega*domega-1.0d0))
      deps3   = 1.0d0 - domegpe*domegpe/        (domega*domega)
      eps1    = deps1
      eps2    = deps2
      eps3    = deps3
c
      npar2 = npar*npar
      aquad = deps1
      bquad = - ( (deps1-npar2)*(deps1+deps3) - deps2*deps2 )
      cquad = deps3 * ( (deps1-npar2)*(deps1-npar2) - deps2*deps2 )
      dquad = bquad*bquad - 4.0*aquad*cquad
      
      if (dquad .gt. 0.0d0) then
        rootqd = DSQRT (dquad)
      else
        rootqd = 0.0d0
      end if
      nper2=(-bquad+nsigma*DSIGN(1.0d0,omega-1.0d0)*rootqd)/
     &(2.0d0*aquad)

c
c     Let nper2 become (only) slightly negative:
c
      if (nper2 .le. 0.0d0) then
        if (nper2 .le. -1.0d-2) then
          !YuP/was: STOP 'subroutine GRPDE2: nper2 is too negative'
          !YuP[2018-05-23] Changed to
          WRITE(*,*)'sub.GRPDE2: nper2 is too negative. nperol^2, nper2'
     +      ,nperol**2, nper2
          vgrpdc=0.d0 !it will stop the ray (but would not halt the run)
          edenfac=0.d0
          return
        else
          nper2=1.0d-10
          nper=1.0d-5
        end if
      else
        nper = DSQRT (nper2)
      end if
c
      n2=nper2+npar2

      d11=eps1-npar2
      d12=-rootm1*eps2
      d13=npar*nper
      d22=eps1-n2
      d33=eps3-nper2
c 
      
      exde=d22/d12
      ezde=-d13*exde/d33
      denom = DSQRT (dble(exde*dconjg(exde)+1.0d0+ezde*dconjg(ezde)))
      exde=exde/denom
      eyde=1.0d0/denom
      ezde=ezde/denom
c      write(*,*)'in grpde2  exde,eyde,ezde',exde,eyde,ezde

c
c calculate derivatives of eps's wrt omega
c
      dep1do=omegpe*omegpe*2.0d0*omega/
     .        ((omega*omega-1.0d0)*(omega*omega-1.0d0))
      dep2do=-eps2/omega-eps2*2.0d0*omega/(omega*omega-1.0d0)
      dep3do=omegpe*omegpe*2.0d0/(omega**3)
c
c find derivatives of D the dispersion reln wrt omega, kper and kpar
c
      dddnpp=-(2.0d0*nper)*( (eps1-npar2)*(eps1-n2+eps3-nper2)
     .                   - eps2*eps2 + npar2*(eps1-n2-nper2) )
      dddnpl=-(2.0d0*npar)*( (eps3-nper2)*(eps1-npar2+eps1-n2)
     .                   + nper2*(eps1-n2-npar2) )
c
      dddome=dddnpp*(-nper/omega) + dddnpl*(-npar/omega)
     .     + dep1do*( (eps3-nper2)*(eps1-n2+eps1-npar2) - nper2*npar2 )
     .     + dep2do*( -2.0d0*eps2*(eps3-nper2) )
     .     + dep3do*( (eps1-npar2)*(eps1-n2) - eps2*eps2 )
c
c find vgroup/c
c
      vgrpdc = DSQRT (dddnpp*dddnpp+dddnpl*dddnpl) / DABS (dddome*omega)
      
c     energy density factor
c
c      1  | ~ | 2    1  ~*      ~       ~     ~
c   =  -  | B |    + -  E .tens.E  with B and E normalized to |E|
c      2  |   |      2       _      _
c                         d |     h  |    h
c    and tens the tensor  - |  w e   |   e  is the hermitian component
c                         dw|_      _!
c
c    of the dielectric tensor
c
      bmod2 = DBLE ( n2*eyde*dconjg(eyde)
     .    + (npar*exde-nper*ezde)*dconjg(npar*exde-nper*ezde) )
c      write(*,*)'in grpde2  bmod2',bmod2
      tens11=eps1+omega*dep1do
      tens12=-(eps2+omega*dep2do)
      tens21=-tens12
      tens22=eps1+omega*dep1do
      tens33=eps3+omega*dep3do
      emod2 = DBLE (dconjg(exde)*exde*tens11
     .            + dconjg(exde)*eyde*rootm1*tens12
     .            + dconjg(eyde)*exde*rootm1*tens21
     .            + dconjg(eyde)*eyde*tens22
     .            + dconjg(ezde)*ezde*tens33)
c      write(*,*)'in grpde2  emod2',emod2
c
      edenfac=0.5d0*(bmod2+emod2)
c      write(*,*)'in grpde2  edenfac',edenfac

c
c     nper=nperol
      return
c
      end


c      
