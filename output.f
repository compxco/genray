
c        ********************** outpt***********************
c        *                      -----                       *
c        *   output is used by drkgs as	an		    *
c        *   output subroutine. it has not to change the    *
c        *   values of its formal input  parameters.  if    *
c        *   prmt(5) is not equal to zero,  the  control    *
c        *   will be transferred to the main program.       *
c        *       its formal parameters are:                 *
c        *         x, y, dery, ihlf, ndim, prmt             *
c        *                                                  *
c        ****************************************************
c        !   this code prints: values of u(ndim);	    !
c        !   creates the data for 3D code; 		    !
c        !   it writes array u(6) in file( the name of      !
c        !   file is given as  the first parameter in	    !
c        !   genray.in file),                   	    !
c        !   eps- value of Hamiltonian. 		    !
c        !   It controls the conditions: if the ray point   !
c        !   is inside the plasma or not .		    !
c        !   It controls the conditions: if the refractive  !
c        !   index less the cnmax or not .		    !
c        !   It writes the output date in mnemonic.txt file	    !
c        !   It corrects the trajectory for Hamiltonian     !
c        !   conservation 
c        !   It scatters the perpendicular refractive index
c----------------------------------------------------------------------
c         input parameters:,t,u,deru,ihlf,ndim,prmt,iflagh,ht
c         output parameters:
c            if iraystop=1  then end of the j_ray calculation
c            if iraystop=0  then continuation of the j_ray calculations
c            iflagh=(2 rays outside  the plasma after the correction;
c                    1 ray is near the plasma boundary after Runge-Kutta
c                      procedure (it is after reflection)
c                    3 ordinary situation ray is inside the plasma
c                      after correction procedure)
c            u(i) after correction and reflection
c            it prepares and writes  parameters for 3d and onetwo
c            codes
c            it changes the u():
c            1)due to correction procedure that conserves Hamiltonian D=0
c            2)due to boundary reflection in reflection points
c            3)due to n_perp scattering after reflection
c            4)due to n_perp scattering in the given point rhoscat(1:nscat_n)
c
c
c           For i_ox.eq.1 it calculates antenna vertex coordinates
c           z_st_ox,r_st_ox,phi_st_ox,alpha_st_ox,beta_st_ox
c           and puts these coordinates into  cone_ec
c           For i_ox.eq.2 it creates OX mode conversion jump
c           in the small radial direction and calculates the 
c           OX transmission coefficient: transm_ox
c-----------------------------------------------------------------*
c           this program uses the following functions and	  *
c           subroutins b,gamma1,hamilt1,  bound,prep3d,write3d    *
c                      scatperp                                   *
c-----------------------------------------------------------------*
      subroutine outpt(t,u,deru,ihlf,ndim,prmt,iflagh,iraystop,
cSAP120517
c     &i_go_to_previous_output_step)
     &i_go_to_previous_output_step,s_total,s_total_old )
cSAP1200603
c      implicit double precision (a-h,o-z)
      implicit none
      include 'param.i' 
      include 'one.i'
      include 'write.i'
      include 'cefield.i'
      include 'cone.i'
c-----the following variable is for N_perp scatterin procedure   
      include 'scatnper.i'
c-----the following to get xma,yma for ox_conversion
      include 'three.i'


c-----the following include is for hot dispersion with electron+ion, 
c     it has the arrays for the correction subroutines 
      include 'nperpcom.i'
      include 'ions.i'
      
c for test only
      double complex chamilt
cend test

      include 'output.i'
      include 'oxb.i'
    
c-----input
cSAP120603
      real*8
     & t, !time or length
     & s_total,     !total ray length used for Nicola scatterin
     & s_total_old  !total length from previous time step

      real*8
     & u(6),    ! ray 6D coordinates
     & deru(6), ! ray coordinates derivatives du/dt
     & prmt(9) ! numerical parameters for RK subroutine

      integer
     &ihlf,      !do not used
     &ndim,
     &iflagh,
     &i_go_to_previous_output_step
c-----output
      integer
     &iraystop
c-----locals      
      real*8
     &pu(6),
     & r_old,       !major radius at previous point
     &t_old,        !t=poloidal length at previus point
     &ham_old,       ! Hamiltonian at previous step
     &prmt_7_old,
     &uout1,uout2,uout3,uout4,uout5,uout6,  
     &z1,r1,phi1,cnz1,cnr1,cm1,cnphi1,
     &cnz,
     &dc2,cos_phi,ds2,sin_phi,
     &dpsi_dz,dpsi_dr,
     & accurcy,eps1,dlambd,dip,eps,epscor,eps_xe_loc,uh,
     &hamilt_old,
     &ad,bd,cd,cn2,cnz_t,cnr_t,cm_t,cnznew,cnrnew,
     &cnpar2, cnper2,cnpar,cnper,cnper_test, cnpar_test,cn0_s,
     &cnmode,cnmax,cnz0,cnr0,cnphi0,cnprim,
     &b_z0,b_r0,b_phi0,cm0,
     &z_ref,r_ref,phi_ref,cnzref,cnrref,cmref,
     &z0,r0,phi0,
     &rnew,
     &clight,vt_e,rme,temp_e,
     &del_uh,del_y0,
     &r_st,phi_st,z_st,alpha_st,beta_st,
     &r_x,z_x,phi_x,cnr_x,cnz_x,cm_x,cnteta,cny,cnx,
     &cnper_max_ebw,
     &dif_phi,sum,vgrmods,yj,xi,xe,yi,ye,yma_loc,xma_loc,
     &y0_cr, Y_abs,s_total_new


      integer
     &iter,i,j,k,n_gam,n_nperp,itermax, naccurc,
     &i_ox_conversion_loc,iflref,iboundc,ibound,icor,ihermloc,
     &ioxm_loc,id_initial,i_output_data,icone,
     &m_limiter,inside_limiter,nearest_lim_boundary
     
      integer ibmx !local

c-----externals
      real*8 temperho,y,x,b,gamma1,tempe,cn,tpoprho,vflowrho,
     &hamilt1


c      save id_initial !to change Ono dispersion
c      logical first
c      data first /.true./
    
      save r_old,t_old
    
c      if (first) then
c        id_initial=id
c        first=.false.
c        r_old=u(2)
c      endif

c      write(*,*)'outpt t',t
c      write(*,*)'nstep_rk',nstep_rk
c      write(*,*)'outpt prmt(6),prmt(7)', prmt(6),prmt(7)
c-----------------------------------------------------------------------
c     test 1
c      if (nstep_rk.eq.50) then
c        z_in=u(1)
c        r_in=u(2)
c        phi_in=u(3) 
c        cnz_in=u(4)
c        cnr_in=u(5)
c        cnphi_in=u(6)/u(2) 
c        bmod=b(z_in,r_in,phi_in)
cc-------trasnform the vector coordinates Nr_in,Nphi_in to Nx_in,Ny_in
c        call transform_N_rphiz_to_xyz(cnr_in,cnphi_in,phi_in,
c     &   cnx_in,cny_in)
cc-------creates the tangent to magnetic surface
c        !components  of the refractive index cntheta,cnphi  
c        call ninit_ec(z_in,r_in,phi_in,cnx_in,cny_in,cnz_in,
c     &  cntheta,cnphi)

cc-------calculations of the parallel (to the magnetic field)
c        !refractive index component cnpar
c        !N_par=(B_phi*N_phi+N_theta*e_theta*{e_z*B_z+e_r*B_r})/bmod
c        !e_theta=(e_z*dpsi/dr-e_r*dpsi/dz)/abs(grad(psi))
c        ppp=(dpdzd*dpdzd+dpdrd*dpdrd)
c        gradpsi=dsqrt(dpdzd*dpdzd+dpdrd*dpdrd)
c        cnpar=(cnphi*bphi+cntheta*(bz*dpdrd-br*dpdzd)/gradpsi)/
c     &  bmod

cc        id_loc=id
cc        ioxm_loc=ioxm 
cc        id=2

c        cn2=cnpar**2+cnper2
c        cnrho2=cn2-cntheta**2-cnphi**2
c	cnrho=dsqrt(cnrho2)
c        cnperp=cnrho
          
c        call cnzcnr(z_in,r_in,phi_in,cntheta,cnphi,cnrho,cnz,cnr,cm)
cc        id=ioxm_loc


c        write(*,*)'output u(4),u(5),u(6)', u(4),u(5),u(6)
c        write(*,*)'output cnz,cnr,cm',cnz,cnr,cm
c        stop 'output'
c      endif
cc-----end_test_1

c      write(*,30)r0x*u(2),r0x*u(1),u(3),u(4),u(5),u(5)
      iraystop=0
      iflagh=3
c-----------------------------------------------------------

      z1=u(1)
      r1=u(2)
      phi1=u(3)
      cnz1=u(4)
      cnr1=u(5)
      cm1=u(6)
      iter=0

      cnphi1=cm1/r1
      cnpar=(bz*cnz1+br*cnr1+bphi*cnphi1)/bmod
     
      cnpar2=cnpar*cnpar
      cn2=cnz1*cnz1+cnr1*cnr1+cnphi1*cnphi1
      cnper=dsqrt(dabs(cn2-cnpar2))
c---------------------------------------------------------------
cSAP090413
c     If the ray intersected the limiter at the
c     toroidal limiter boundary  phi_limiter(j,m)
c     then stopping this ray.
c-------------------------------------------------------------------
c      write(*,*)'output.f max_limiters',max_limiters

      if (max_limiters.ge.1) then

cc         write(*,*)'before limiter z1,r1,phi1',z1,r1,phi1

         call check_if_point_inside_limiter(z1,r1,phi1,
     &   inside_limiter,m_limiter,dif_phi,nearest_lim_boundary)

cc         write(*,*)'after limiter z1,r1,phi1',z1,r1,phi1

         if (inside_limiter.eq.1) then
            write(*,*)'******ray is at toroidal limiter boundary******'
            iraystop=1
cSAP100514
            iray_status_one_ray=4
            return
         endif
      endif 
c------------------------------------------------
     

c------------------------------------------------------------
c     if the number of the time step nstep_rk is bigger that maxsteps_rk
c     then stop ray calculations
      if(nstep_rk.gt.maxsteps_rk) then
        write(*,*)'**********nstep_rk.gt.maxsteps_rk****************'
!nme        if (iout3d.eq.'enable') then
c----------writing data to output mnemonic.txt file
!nme	   call write3d
!nme        endif
	iraystop=1
cSAP100514
        iray_status_one_ray=5
	return
      end if
      nstep_rk=nstep_rk+1
c--------------------------------------------------------------------
c     if ray is close to resonance point then stop ray calculations
      cnmode=cn(r1,cnz1,cnr1,cm1)
      cnmax=10000.d0
      cnmax=1.d+6

      if(cnmode.gt.cnmax) then
	write(*,*)'***************nmode.gt.cnmax***********'
!nme        if (iout3d.eq.'enable') then
c----------writing data to output mnemonic.txt file
!nme	   call write3d
!nme        endif

	iraystop=1
cSAP100514
        iray_status_one_ray=6
	return
      end if

c------------------------------------------------------------------
c      check that cold plasma X mode is close to UH resonance
c      ann can be EBW conversion 
c-------------------------------------------------------------
cSAP070807
c      cnmax_cold=40.d0
c      if(((id.eq.1).or.(id.eq.2)).and.(cnmode.gt.cnmax_cold))then
c	write(*,*)'**************nmode.gt.cnmax_cold***********'
c          iraystop=1
c	return
c      end if

      xe=x(u(1),u(2),u(3),1)
      ye=y(u(1),u(2),u(3),1)
      uh=dsqrt(xe+ye*ye)
      uh=xe+ye*ye
      del_uh=1.d-2
      temp_e=temperho(rho,1) 
      rme=9.1094d-28 
      vt_e=dsqrt(2.d0*temp_e*1.6022d-9/rme) !(cm/sec)
                                             ! thermal velocity
      clight=2.99792458d10
      cnper_max_ebw=ye*(clight/vt_e)
      
      if((((id.eq.1).or.(id.eq.2)).and.(uh.gt.1.d0).and.
     & ((uh-1.d0).lt.del_uh))
     &.and.(cnper.gt.cnper_max_ebw)) then
        write(*,*)'********************************************'
        write(*,*)'For cold plasma dispersion relation id=',id
        write(*,*)'the ray is close to upper hybrid resonance'
        write(*,*)' (uh-1.d0).lt.del_uh, uh=',uh,'del_uh=',del_uh
        write(*,*)'N>1, N=',dsqrt(cn2),'cnper=',cnper
        write(*,*)'EBW condition: cnper*(vt_e/clight)/ye=',
     &  cnper*(vt_e/clight)/ye
        write(*,*)'It can be Xmode to close to the UH resonance'
        write(*,*)'Xmode can be transformed to EBW'
        write(*,*)'For Xmode- EBW convertion hot plasma dispersion'
        write(*,*)'can be used id=6,9'
        write(*,*)'The ray calculation stopped' 
        write(*,*)'********************************************'
        iraystop=1
cSAP100514
        iray_status_one_ray=7
	return
      endif

      if (t.gt.poldist_mx) then
	 write(*,*)'***************t.gt.poldist_mx***********'
!nme         if(iout3d.eq.'enable') then
c----------writing data to output mnemonic.txt file
!nme           write(*,*)'in outpt: t.gt.poldit_mx before write3d'
!nme	   call write3d
!nme	   write(*,*)'in output after write3d'
!nme         endif
      	 iraystop=1
cSAP100514
         iray_status_one_ray=1
	 return
      end if

c     if Group Velocity, normalized to c, is greater than 1.1
c     (give 10 percent grace) then stop the ray.  [RWH: 030427].

      vgrmods=deru(1)**2+deru(2)**2+(u(2)*deru(3))**2

      if ((vgrmods.gt.1.1).and.
     1    ((id.ne.10).and.(id.ne.12).and.(id.ne.13)))
     1   then
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*) 'outpt: Stopping ray.'
         write(*,*)' vgroup>1.1,   abs(vgroup) = ',dsqrt(vgrmods)  
         write(*,*) '*************************************************'
         write(*,*)

!nme         if(iout3d.eq.'enable') then
c-----------writing data to output mnemonic.txt file
!nme            write(*,*)'in outpt: t.gt.poldit_mx before write3d'
!nme            call write3d
!nme            write(*,*)'in output after write3d'
!nme         endif!nme        
!nme         endif
cSm060830
c            if(iout3d.eq.'enable') then
c------------writing data to output mnemonic.txt file
c             write(*,*)'in outpt: vgrmods.gt.1 before write3d'
c             call write3d
c             iraystop=1
c             write(*,*)'in output after write3d'
c            endif! 

CENM 1Sep05 -- best to really stop the ray when using the relativistic
C    dispersion relation, otherwise it goes on and on without much progress.
C    vgrmods.gt.1.1 usually when it has nearly reached full depletion of power
C    in the ray anyway.
      	 if (id.eq.14) iraystop=1
cSAP100514
         iray_status_one_ray=8 
cBH131024c      	 iraystop=1
         iraystop=1
cSAP090724
cBH131024c	 return
         return
      end if

c---------------------------------------------------------------------
c     change of the dispersion function near the cyclotron resonance 
      bmod=b(z1,r1,phi1)
      
ctest
c      do i=1,nbulk
c        yj=y(z1,r1,phi1,i)
c        write(*,*)'output i,1/y',i,1.d0/yj
c      enddo
cendtest

      yj=y(z1,r1,phi1,jy_d)
c      write(*,*)'output 1/yj',1.d0/yj
cSm070730
      if(iswitch.eq.1) then
cNazikian case 070802: comment previous line and uncomment folowing lones
c       cnt=cn(r1,cnz1,cnr1,cm1)
c      if (cnt.gt.5) then 
c        i_geom_optic=2
c        ray_direction=-1
c      endif 
c      write(*,*)'output id',id,'i_geom_optic',i_geom_optic
c      write(*,*)'ray_direction',ray_direction
c
c      if((iswitch.eq.1).and.(cnt.lt.3.d0)) then
c_end Nazikian case

        if (id_initial.eq.8) then !Ono dispersion
          del_y0=0.03d0
          y0_cr=0.95d0
          y0_cr=0.94d0
          del_y0=0.04d0
          idswitch=2
          call switch_d_ono(r1,z1,phi1,cnz1,cnr1,cm1,del_y0,y0_cr,
     &    id,idswitch)
        else 
c-------  change of the dispersion function and the absorption subroutines
c         near the cyclotron resonance points       
c         write(*,*)'yj',yj,',id',id
          call switch_da(yj,del_y,id,iabsorp,idswitch,iabswitch)
c         write(*,*)'output: id,iabsorp',id,iabsorp
        endif
      endif
      
c	write(*,*)'output before corrections'
c      write(*,*)'in output befe cor u(1-3)',u(1),u(2),u(3)
c      write(*,*)'in output befe cor u(4-6)',u(4),u(5),u(6)
c---------------------------------------------------------
c     correction
c---------------------------------------------------------
c     the  switch off the Hamiltonian correction procedure
      if(icorrect.eq.0) goto 11
c---------------------------------------------------------
      epscor=prmt(4)
      bmod=b(z1,r1,phi1)
      gam=gamma1(z1,r1,phi1,cnz1,cnr1,cm1)
      eps=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))
c      write(*,*)'output after hamilt1 eps',eps
      eps1=eps
  
      cnphi1=cm1/r1
      cnpar=(bz*cnz1+br*cnr1+bphi*cnphi1)/bmod
     
      cnpar2=cnpar*cnpar
      cn2=cnz1*cnz1+cnr1*cnr1+cnphi1*cnphi1
      cnper=dsqrt(dabs(cn2-cnpar2))

c      write(*,*)'output eps,epscor',eps,epscor

      if(dabs(eps).lt.epscor) goto 53

c      write(*,*)'ouput after goto 53'

c      if ((cnper.ge.1000.d0).and.
c      write(*,*)'output cnperid',cnper,id
cSAP080829
c      goto 60

      if ((cnper.ge.0.d0).and.
     *(((id.eq.4).or.(id.eq.5)).or.((id.eq.6).or.(id.eq.7)))) then
c         write(*,*)'output before call correct2'
c--------correction from the solution n_perp=n_perp(n_parallel)
c        Now it is for id=4,5,6,7
         
        do j=1,nbulk
          massc(j)=dmas(j)
          xc(j)=x(z1,r1,phi1,j)
          yc(j)=y(z1,r1,phi1,j)
cSm000324
          if(j.eq.1) yc(1)=-yc(1)
          tec(j)=tempe(z1,r1,phi1,j)*1.d+3 !(eV) averaged temperature
          tpopc(j)=tpoprho(rho,j)
          vflowc(j)=vflowrho(rho,j)         
        enddo
 
        cnpar=(bz*cnz1+br*cnr1+bphi*cm1/r1)/bmod
        accurcy=epscor 
        naccurc=5
        cnper=dsqrt(dabs((cnz1**2+cnr1**2+(cm1/r1)**2)-cnpar**2))
c        write(*,*)'output before solvnperp cnpar,cnper',cnpar,cnper
        ihermloc=iherm
        
c        if(dabs(yc(1)).lt.1.d0) then
c          write(*,*)'in output in correction before solvnperp y<1'  
c          goto 11
c        endif

        call solvnperp(nbulk,massc,xc,yc,tec,tpopc,vflowc,cnpar,id,
     *  ihermloc,accurcy,naccurc,cnprim,cnper)
c        write(*,*)'output after solvnperp cnpar,cnper',cnpar,cnper
c        write(*,*)'output before correct2 cnz1,cnr1',cnz1,cnr1
        
        call  correct2(cnpar,cnper,
     *  cnz1,cnr1,cm1,r1,bz,br,bphi,bmod,cnznew,cnrnew)
        goto 21

        eps=1.d-8
        itermax=20
c        write(*,*)'output before correct3 z1,r1,phi1',z1,r1,phi1
        call correct3(cnpar,cnper,cnz1,cnr1,cm1,z1,r1,phi1,
     .  eps,itermax,cnznew,cnrnew,rnew)
c        write(*,*)'output after correct3 rnew,cnznew,cnrnew',
c     .  rnew,cnznew,cnrnew        
         u(2)=rnew

 21     continue
c        write(*,*)'output after correct2 cnznew,cnrnew',cnznew,cnrnew        
        u(4)=cnznew
        u(5)=cnrnew

        goto 11
      endif 
cSAP080829
 60   continue

c--------correction from the problem
c        sum{delt(x_i)**2}=min
c        disp_func(x+delt(x))=0 
	 iter=1
	 icor=0 
c----------------------------------------------------------------------  
ctest
         bmod=b(z1,r1,phi1)
         gam=gamma1(z1,r1,phi1,cnz1,cnr1,cm1)
         eps=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))
c         write(*,*)'output test before correction_lagrange_1 eps',eps 
c         write(*,*)'bmod,gam',bmod,gam
c         write(*,*)'u(1),u(2),u(3),u(4),u(5),u(6)',
c     &              u(1),u(2),u(3),u(4),u(5),u(6)
cendtest
     
c         iter_max=20
c         call correction_lagrange_1(u(1),u(2),u(3),u(4),u(5),u(6),
c     &                              iter_max)
c         goto 53
c---------------------------------------------------------------------
c         write(*,*)'output before 52 eps,eps1',eps,eps1
 52	 continue
         dip=1.d0
c	 write(*,*)'output before dddrz1'
         call dddrz1(t,u,deru)
c         write(*,*)'output after dddrz1'
cSm070116
c	 sum=deru(1)**2+deru(2)**2+deru(4)**2+deru(5)**2
	 sum=deru(1)**2+deru(2)**2
         hamilt_old=eps1
 12      dlambd=2.d0*eps/sum*dip
c-----------------------------------------------------------
c         write(*,*)'12 dip=',dip
	 do 50 i=1,2
cSm070116
c 50	   pu(i)=u(i)+0.5*dlambd*deru(i+3)
 50	   pu(i)=u(i)

         write(*,*)'dlambd',dlambd

	 do 51 i=4,5
 51      pu(i)=u(i)-0.5*dlambd*deru(i-3)

         pu(3)=u(3)
         pu(6)=u(6)
c-----------------------------------------------------------
         z1=pu(1)
         r1=pu(2)
         phi1=pu(3)
         cnz1=pu(4)
         cnr1=pu(5)
         cm1=pu(6)
c         write(*,*)'output  before boundc r1,z1',r1,z1
c         write(*,*)'output  cnz1,cnr1,cm1', cnz1,cnr1,cm1
         call boundc(z1,r1,iboundc)

c         write(*,*)'output after boundc iboundc',iboundc

	 if (iboundc.eq.1) then
c           write(*,*)'in output iboundc=1 lflagh=2'
c	   write(*,*)'outside the plasma after the correction'
c---------------------------------------------------------------------
c          corrected ray point is outside the plasma
c          it is necessary to reduce the time step in the Runge -Kutta
c          procedure  and to recalculate the array u(i)
c----------------------------------------------------------------------
           iflagh=2
c           goto120

c          after the correction step ray is out of plasma in up()
c          stop the correction procedure at the last correction step u()
           iflagh=3
           goto 11
	 endif

         bmod=b(z1,r1,phi1)
         gam=gamma1(z1,r1,phi1,cnz1,cnr1,cm1)

cSAP080802
         hamilt_old=eps !hamiltonian at previous iteration step

c         write(*,*)'output pu(1),pu(2),pu(3),pu(4),pu(5),pu(6)',
c     &                     pu(1),pu(2),pu(3),pu(4),pu(5),pu(6)

         eps=hamilt1(pu(1),pu(2),pu(3),pu(4),pu(5),pu(6))
c        epsn=hamilt1(pu(1),pu(2),pu(3),pu(4),pu(5),pu(6))
c         write(*,*)'in output after correction 1 hamilt_old,epsham=',
c     &                                           hamilt_old,eps
c         write(*,*)'eps1',eps1
c-----------------------------------------------------------
cSAP080802
c         if(dabs(eps) .gt. abs(hamilt_old))then
c            write(*,*)'in output eps.gt.hamilt_old',eps,hamilt_old
c           
cc            write(*,*)'cnz1,cnr1',cnz1,cnr1
c            m_nz=100
c            m_nr=100
c            n_contour=20
c            n_contour=10
c            cmax_nz=cnz1*1.5d0
c            cmin_nz=cnz1*0.5d0
c            cmax_nr=cnr1*1.5d0
c            cmin_nr=cnr1*0.5d0
c            cmax_nz=cnz1*1.05d0
c            cmin_nz=cnz1*0.95d0
c            cmax_nr=cnr1*1.01d0
c            cmin_nr=cnr1*0.95d0

c            call map_dcold_nz_nr(m_nz,m_nr,
c     &                          cmax_nz,cmin_nz,cmax_nr,cmin_nr,
c     &                          n_contour,z1,r1,phi1,cnz1,cnr1,cm1)
c            iraystop=1
c            write(*,*)'**********dabs(eps).gt.abs(hamilt_old**********'
c	    return

c            goto 53 !to end of corrections

cyup         write(*,*)'output iter,eps',iter,eps
cSAP090321
         if(eps.gt.1.d120) then
c         if(eps.gt.1.d20) then
           write(*,*)'in output eps.gt.1.d20 iter=',iter
           write(*,*)'Hamiltonian correction procedure stoped'
           write(*,*)'Hamiltonian=',eps
           goto 53
	 end if
  
	 if (dabs(eps).gt.dabs(eps1)) then
cyup	     write(*,*)'in output eps.gt.eps1,dip',eps,eps1,dip             
	     dip=0.5d0*dip
	     goto 12
	 end if
         
c-------------------------------------------------------------
	 eps1=eps
c-------------------------------------------------------------
	 do 14 i=1,ndim
 14      u(i)=pu(i)

	 iter=iter+1
	 if(dabs(eps).lt.epscor) goto 53
         
	 if(iter.gt.20)then
	   write(*,*)'Hamiltonian correction procedure
     1	   made 20 iterations and stopped , Hamiltonian=',eps
	   goto 53
	 end if

	 goto 52
 53      continue

c         write(*,*)'output after iterations, label 53'
c         write(*,*)'output u(1),u(2),u(3),u(4),u(5),u(6)',
c     &                     u(1),u(2),u(3),u(4),u(5),u(6)

c--------------------------------------------------------------
c     end of correction
c--------------------------------------------------------------
 11   continue
c--------------------------------------------------------------
c     measure error in the dispersion realation
c     If D/(N|gradD|) > toll_hamilt stop ray calculation
c--------------------------------------------------------------

c      write(*,*)'output.f before refractive_index_relative_error'

      call refractive_index_relative_error(u(1),u(2),u(3),u(4),u(5),u(6)
     &,iraystop)

c      write(*,*)'output.f after refractive_index_relative_error iraystop
c     &',iraystop

      if(iraystop.eq.1) then
        write(*,*)
         write(*,*) '*************************************************'
         write(*,*) 'outpt: Stopping ray.'
         write(*,*) ' D/(N|gradD|) > toll_hamilt'  
         write(*,*) '*************************************************'
         write(*,*)
cSAP100514
         iray_status_one_ray=9 
         return
      endif
c------------------------------------------------------------------
c     Creates the jump of the ray point throw the OX mode conversion
c     area where X_e=1 (V_perp=0)
c-------------------------------------------------------------------
      if (i_ox.eq.2) then
        write(*,*)'was_not_ox_conversion',was_not_ox_conversion
        i_call_prep3d_in_output_at_i_ox_conversion_eq_1=0
        if (was_not_ox_conversion) then
          xma_loc=xma
          yma_loc=yma
          eps_xe_loc=eps_xe
          
c---------calculate the coordinates of X mode cutoff point
c---------r_x,z_x,phi_x,cnr_x,cnz_x,cm_x 
cyup          write(*,*)'output'
c          write(*,*)'before ox_conversion u(1),u(2),u(3),u(4),u(5),u(6)'
c     &    ,u(1),u(2),u(3),u(4),u(5),u(6)

ctest
c          call  transmit_coef_ox(u(1),u(2),u(3),u(4),u(5),u(6),
c     &                           transm_ox_l)
c          write(*,*)'before ox_conversion transm_ox_l',transm_ox_l
cendtest

          call  ox_conversion(u(2),u(1),u(3),u(5),u(4),u(6),
     &    xma_loc,yma_loc, !temporally
     &    eps_xe_loc,
     &    r_x,z_x,phi_x,cnr_x,cnz_x,cm_x,i_ox_conversion_loc)
          i_ox_conversion=i_ox_conversion_loc

cyup          write(*,*)'i_ox_conversion=',i_ox_conversion
c          write(*,*)'after ox_conversion u(1),u(2),u(3),u(4),u(5),u(6)'
c     &    ,u(1),u(2),u(3),u(4),u(5),u(6)

ctest
c          call  transmit_coef_ox(u(1),u(2),u(3),u(4),u(5),u(6),
c     &                           transm_ox_l)
c          write(*,*)'after ox_conversion transm_ox_l',transm_ox_l
cendtest
          
          if (i_ox_conversion.eq.1) then
c           calculate transmission coefficient for OX mode
c           conversion:transm_ox

c--------------------------------------------------------
c           The following data are for write.i to prepare 
c           the output data for mnemonic.nc file            
              
            Y_abs=dabs(y(u(1),u(2),u(3),1))
            cn_par_optimal=dsqrt(Y_abs/(Y_abs+1)) !optimal N parallel
                                                  !for OX conversion
            bmod=b(u(1),u(2),u(3))
            cnpar_ox=(bz*u(4)+br*u(5)+bphi*u(6)/u(2))/bmod ! N_parallel
                                                           ! before OX
                                                           ! conversion
            cn_b_gradpsi=(u(5)*bphi*dpdzd-
     &      (u(6)/u(2))*(br*dpdzd-bz*dpdrd)-u(4)*bphi*dpdrd)/
     &   dsqrt((bphi*dpdzd)**2+(br*dpdzd-bz*dpdrd)**2+(bphi*dpdrd)**2)


c-------------------------------------------------------
   
c-----------transmission coefficient in the point before the jump 
            call  transmit_coef_ox(u(1),u(2),u(3),u(4),u(5),u(6),
     &                           transm_ox)
cyup            write(*,*)'output OX mode conversion coefficient transm_ox'
cyup     &      ,transm_ox

c--------------------------------------------------------------------
c           write the data fo O mode  before jump
            i_call_prep3d_in_output_at_i_ox_conversion_eq_1=1
c            write(*,*)'output before ox jump before prep3d'
c            write(*,*)'u(1),u(2)',u(1),u(2)
            call prep3d(t,u,deru,iraystop)
c--------------------------------------------------------------------

            u(1)=z_x
            u(2)=r_x
            u(3)=phi_x
            u(4)=cnz_x
            u(5)=cnr_x
            u(6)=cm_x
cyup            write(*,*)'output after jump before prep3d u= ',
cyup     &     u(1),u(2),u(3),u(4),u(5),u(6)

c-----------transmission coefficient in the point after  the jump
c           call  transmit_coef_ox(u(1),u(2),u(3),u(4),u(5),u(6),
c     &                           transm_ox)
c            write(*,*)'output OX mode conversion coefficient transm_ox'
c            write(*,*)'after jump',transm_ox
c------------------------------------------------------------------------
c           write the data for x mode  after jump
            i_call_prep3d_in_output_at_i_ox_conversion_eq_1=2
            call prep3d(t,u,deru,iraystop)
c--------------------------------------------------------------------
            was_not_ox_conversion= .false.
          endif ! i_ox_conversion.eq.1
        endif   ! was_not_ox_conversion 
      endif     ! i_ox.eq.2
c-------------------------------------------------------------------


c       write(*,*)'output iter=',iter
c      write(*,*)'outpt after cor u(2),u(1)',u(2),u(1)
c      write(*,*)'in output after cor deru(1),deru(2)',deru(1),deru(2)
c-----------------------------------------------------------
cSAP090515
c       write(*,*)'090515 output t=',t,'prmt(7)',prmt(7)
     
c       call  transmit_coef_ox(u(1),u(2),u(3),u(4),u(5),u(6),
c     &                           transm_ox)
c       write(*,*)'transm_ox',transm_ox

c      write(*,*)'output n_plot_disp',n_plot_disp
c      write(*,*)'output point_plot_disp ',point_plot_disp

c--------------------------------------------------------------------
c     plot hot plasma dispersion function contours D(ReN_perp,Im_N_perp)
c     in given poloidal distance or major radius
c---------------------------------------------------------------------
      if( point_plot_disp.eq.'rz_nparallel_points') goto40 !do not plot here

      if(n_plot_disp.gt.0) then
        !plot dispersion function contours D(ReN_perp,Im_N_perp)
        do k=1,n_plot_disp

          if (point_plot_disp.eq.'poloidl_dist')then

c            write(*,*)'t,t_old,k,s_poloid_plot_disp(k)',
c     &                 t,t_old,k,s_poloid_plot_disp(k)

            if((t-s_poloid_plot_disp(k))*(t_old-s_poloid_plot_disp(k)).
     &        lt.0.d0)then
cyup            write(*,*)'in output plot contourD(Re_N_perp,Im_N_perp) t',t
     
              call output_map_disp_nperp(u)    
    
              goto 40
            endif
          endif

          if (point_plot_disp.eq.'major_radius')then
c            write(*,*)'output.f k,u(2),r_old,r_plot_disp(k)',
c     &                        k,u(2),r_old,r_plot_disp(k)
            if((u(2)-r_plot_disp(k))*(r_old-r_plot_disp(k)).lt.0.d0)then
              write(*,*)'**************************************'
            write(*,*)'in output plot contourD(Re_N_perp,Im_N_perp) t',t
              write(*,*)'major radius k,u(2),r_plot_disp(k),r_old',
     &                              k,u(2),r_plot_disp(k),r_old
              call output_map_disp_nperp(u)         
              goto 40
            endif
          endif

        enddo !k
      endif !n_plot_disp.gt.0
 40   continue
     
c-----------------------------------------------------------------------------
c     end of plot hot plasma dispersion function contours D(ReN_perp,Im_N_perp)
c-----------------------------------------------------------------------------

c--------------------------------------------------------------------
c     plot cold dispersion function D(ReN_perp)
c     in given poloidal distance or major radius
c---------------------------------------------------------------------
      
      if(n_plot_disp_cold.gt.0) then
        !plot dispersion function D_cold(ReN_perp)

        n_nperp=1000     ! set the number of nperp points at the plot D(N_perp)
        n_gam=500        ! set the number of angle points at the interval (0,2*pi) 
                         ! at the wave normal surface plot
       

        do k=1,n_plot_disp_cold
       
          if (point_plot_disp_cold.eq.'poloidl_dist')then

c            write(*,*)'t,t_old,k,s_poloid_plot_disp_cold(k)',
c     &                 t,t_old,k,s_poloid_plot_disp_cold(k)

            if((t-s_poloid_plot_disp_cold(k))*
     &         (t_old-s_poloid_plot_disp_cold(k)).lt.0.d0)then


               call plot_disp_cold_output(u,t,n_nperp)

              if (i_plot_wave_normal_cold.eq.1) then
c---------------plot wave normal surfaces
                bmod=b(u(1),u(2),u(3))

                cnpar=(bz*u(4)+br*u(5)+bphi*u(6)/u(2))/bmod
                cn2=u(4)**2+u(5)**2+(u(6)/u(2))**2

                if((cn2-cnpar**2).ge.0.d0)then
                  cnper=dsqrt(dabs(cn2-cnpar**2))
                else
                  cnper=-1.d0
                endif

cyup                write(*,*)'cn2,cnpar',cn2,cnpar

cyup                write(*,*)'in output before  wave_normal_surface'
cyup                write(*,*)'cnpar,cnper,dsqrt(cnpar**2+cnper**2) ',
cyup     &                     cnpar,cnper,dsqrt(cnpar**2+cnper**2) 

                ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
                call wave_normal_surface(u(1),u(2),u(3),
     &               cnpar,cnper,t,n_gam,ibmx)
 
                call wave_ray_normal_surface(u(1),u(2),u(3),
     &               cnpar,cnper,t,n_gam,ibmx)
               endif

               goto 41
            endif
          endif

          if (point_plot_disp_cold.eq.'major_radius')then
cyup            write(*,*)'output.f k,u(2),r_old,r_plot_disp_cold(k)',
cyup     &                        k,u(2),r_old,r_plot_disp_cold(k)
            if((u(2)-r_plot_disp_cold(k))*(r_old-r_plot_disp_cold(k)).
     &         lt.0.d0)then
 
c              write(*,*)'**************************************'
c               write(*,*)'in output plot D_cold(Re_N_perp) t',t
c               write(*,*)'major r k,u(2),r_plot_disp_cold(k),r_old',
c     &                            k,u(2),r_plot_disp_cold(k),r_old

               call plot_disp_cold_output(u,t,n_nperp)

               if (i_plot_wave_normal_cold.eq.1) then
c---------------plot wave normal surfaces 
                bmod=b(u(1),u(2),u(3))

                cnpar=(bz*u(4)+br*u(5)+bphi*u(6)/u(2))/bmod
                cn2=u(4)**2+u(5)**5+(u(6)/u(2))**2

                if((cn2-cnpar**2).ge.0.d0)then
                  cnper=dsqrt(dabs(cn2-cnpar**2))
                else
                  cnper=-1.d0
                endif

                ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
                call wave_normal_surface(u(1),u(2),u(3),
     &               cnpar,cnper,t,n_gam,ibmx) 

                call wave_ray_normal_surface(u(1),u(2),u(3),
     &               cnpar,cnper,t,n_gam,ibmx)
               endif

              goto 41
            endif
          endif

        enddo !k
      endif !n_plot_disp_cold.gt.0
 41   continue
     
c-----------------------------------------------------------------------------
c     end of plot cold plasma dispersion function D(ReN_perp)
c-----------------------------------------------------------------------------
      r_old=u(2)
      t_old=t
      i_output_data=0 !this time step is not the output step

cSAP090515
c      write(*,*)'output.f t,prmt(7)',t,prmt(7) 
     
      if(t.lt.prmt(7))goto 20
      i_output_data=1 !this time step is the output step
c-------------------------------------------------------
cyup      write(*,*)'tau=',t,' 1/y(1)=',1.d0/y(z1,r1,phi1,1),
cyup     &'1/y(2)=',1.d0/y(z1,r1,phi1,2)
      
      
c     The work on the data preparation for the output
c     at given time steps
c---------------------------------------------------------
c     denormalization of ray point coordinates
c---------------------------------------------------------
      uout1=r0x*u(1)
      uout2=r0x*u(2)
      uout3=u(3)
      uout4=u(4)
      uout5=u(5)
      uout6=u(6)*r0x
c-----------------------------------------------------------
c      write(*,30)uout2,uout1,uout3,uout4,uout5,uout6
cyup      write(*,*)uout2,uout1,uout3,uout4,uout5,uout6

c      write(*,130)uout2,uout1,uout3,uout4,uout5,uout6
130   format(3x,6(' ',e16.9))

c      write(i1_,110)uout2,uout1,uout3,uout4,uout5,uout6
      
c      prmt(7)=prmt(7)+prmt(6)
c
cSAP100202
      prmt_7_old=prmt(7)
c      write(*,*)'prmt_7_old',prmt_7_old

      prmt(7)=(dint(t/prmt(6))+1)*prmt(6)

cSAP090515
cyup      write(*,*)'output.f 1 t,prmt(6),prmt(7)',t,prmt(6),prmt(7)

 110  format(3x,6(' ',e13.6))
c***************************************************
      bmod=b(z1,r1,phi1)

c      write(*,*)'output bef eps z1,r1,phi1',z1,r1,phi1
c      write(*,*)'output bef eps bz,br,bphi,bmod',bz,br,bphi,bmod
c      write(*,*)'output bef eps cnz1,cnr1,cm1',cnz1,crn1,cm1
c       write(*,*)'output rho',rho
c     gam=gamma1(z1,r1,phi1,cnz1,cnr1,cm1)
      gam=gamma1(u(1),u(2),u(3),u(4),u(5),u(6))
      ds2=dsin(gam)**2
      dc2=dcos(gam)**2
c      write(*,*)'output bef eps gam,ds,dc',gam,dsin(gam),dcos(gam)    
      call abc(u(1),u(2),u(3),ds2,dc2,ad,bd,cd)

c      write(*,*)'output bef eps u(1),u(2),u(3)',u(1),u(2),u(3)
c      write(*,*)'output bef eps u(4),u(5),u(6)',u(4),u(5),u(6)
c      write(*,*)'output bef eps bz,br,bphi,bmod',bz,br,bphi,bmod
      
      eps=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))
cyup      write(*,*)'in output epshamilt1=',eps

      cn2=cnz1*cnz1+cnr1*cnr1+(cm1/r1)**2

c testc      write(*,*)'output bef eps z1,r1,phi1',z1,r1,phi1
c      write(*,*)'output bef eps bz,br,bphi,bmod',bz,br,bphi,bmod
c      write(*,*)'output bef eps cnz1,cnr1,cm1',cnz1,crn1,cm1
c       write(*,*)'output rho',rho
      cnpar=(cnz1*bz+cnr1*br+(cm1/r1)*bphi)/bmod
      cnper2=cn2*ds2
      cnper=dsqrt(cnper2)
      xe=x(u(1),u(2),u(3),1)
      ye=y(u(1),u(2),u(3),1)
      uh=dsqrt(xe+ye*ye)
      uh=xe+ye*ye
cyup      write(*,*)'output xe,ye,uh',xe,ye,uh

      if (i_uh_switch.eq.1) then
         !change the output step prmt(6) for prmt6_uh_switch 
         if(uh.lt.uh_switch) then
           prmt(6)=prmt6_uh_switch
         else
         ! regular output step
           prmt(6)=prmt6
         endif
      endif


      if (i_power_switch_resonance.eq.1) then
         !change the output step prmt(6) for prmt6_power_switch_resonance
         do k=1,n_power_switch_resonance 
           if(dabs(ye-y_power_switch_resonance(k)).lt.
     &       del_y_power_switch_resonance) then
             prmt(6)=prmt6_power_switch_resonance
           else
           ! regular output step
             prmt(6)=prmt6
           endif
         enddo
      endif  


c      write(*,*)'output cnper,cnpar',cnper,cnpar
      do i=2,nbulk
         xi=x(u(1),u(2),u(3),i)
         yi=y(u(1),u(2),u(3),i)
c         write(*,*)'output i,xi,yi',i,xi,yi
      enddo      
c end test
 568  format(3x,'eps=',d13.6)
     
cyup      write(*,*)'output.f before prep3d t',t

      call prep3d(t,u,deru,iraystop)
     
cyup      write(*,*)'output aft prep3d iraystop,t',iraystop,t
cSAP100202      
c      write(*,*)'output i_go_to_previous_out_step',
c     &                  i_go_to_previous_out_step
c      write(*,*)'0 output.f prmt(7),prmt(6)', 
c     &                      prmt(7),prmt(6)
      if (i_go_to_previous_out_step.eq.1) then
c--------return to the previous output point an will use new prmt6_new
cyup         write(*,*)'in output.f i_go_to_previous_out_step.eq.1'
cyup         write(*,*)'in output.f prmt_7_old,prmt(6),prmt6,prmt6_new',
cyup     &                          prmt_7_old,prmt(6),prmt6,prmt6_new
         prmt(7)=prmt_7_old-prmt(6)+prmt6_new
         prmt(6)=prmt6_new
c         write(*,*)'output.f nrayelt prmt(6),prmt(6)', 
c     &                       nrayelt,prmt(6),prmt(6)

         u(1)=wz(nrayelt-1)*0.01d0/r0x 
         u(2)=wr(nrayelt-1)*0.01d0/r0x 
         u(3)=wphi(nrayelt-1)
         u(4)=wn_z(nrayelt-1)
         u(5)=wn_r(nrayelt-1)
         u(6)=wn_phi(nrayelt-1)*u(2)
         t=wt(nrayelt-1)
        
cyup         write(*,*)'output.f prmt(7),prmt(6)',prmt(7),prmt(6)
cyup         write(*,*)'output.f u(1),u(2),u(3),u(4),u(5),u(6)',
cyup     &              u(1),u(2),u(3),u(4),u(5),u(6)
cyup         write(*,*)'output.f t',t
         nrayelt=nrayelt-1
         i_go_to_previous_output_step=i_go_to_previous_out_step

         return
      else
         i_go_to_previous_output_step=0
      endif
     


      if (iraystop.eq.1) then

!nme         if(iout3d.eq.'enable') then
c-----------writing data to output mnemonic.txt file
!nme	    call write3d
!nme         endif
cSm060630
c         if(iout3d.eq.'enable') then
c-----------writing data to output mnemonic.txt file
c	    call write3d
c         endif

	 return
      end if

cyup      write(*,*)'!!!!!output nrayelt,nrelt',nrayelt,nrelt
cBH131024: Somehow, a return from this subroutine before
cBH131024: reaching the following test, is causing the test
cBH131024: to not be reached.  Catch it on the next try.
cBH131024      if (nrayelt.eq.nrelt) then
      if (nrayelt.ge.nrelt) then
	 write(*,*)'***************nrayelt=nrelt***********'

!nme         if(iout3d.eq.'enable') then
c----------writing data to output mnemonic.txt file
c	   write(*,*)'in output before write3d
c     1	   nrayelt,nrelt',nrayelt,nrelt
!nme	   call write3d
c	   write(*,*)'in output after write3d'
!nme         endif

      	 iraystop=1
cSAP100514
         iray_status_one_ray =12
	 return
      end if

  20  continue
c---------------------------------------------------
c     end of output data preparation at the given time steps
c--------------------------------------------------
  30  format(1x,6(1x,e11.4))
c-----------------------------------------------------------
c     control of the reflection moment
c-----------------------------------------------------------
c     this call is for the calculation of
c     dzdt=deru(1) and drdt=deru(2).
c     subroutine bound uses deru(1) and deru(2)	to detemine
c     if the ray point goes into or out the plasma.
c      write(*,*)'in output.f before  rside1 t',t

      call rside1(t,u,deru)

c      write(*,*)'in output.f after  rside1 t',t
c-----------------------------------------------------------
c      write(*,*)' 20 output before bound u(2),u(1)',u(2),u(1)
      call bound(u(1),u(2),u(3),u(4),u(5),u(6),iflref,
     &      z_ref,r_ref,phi_ref,cnzref,cnrref,cmref,
     &      ibound,deru(1),deru(2))

c      write(*,*)'20 output after bound ibound,iflref,u(2),u(1)',
c     &           ibound,iflref,u(2),u(1)

      if (iflref.eq.1) then

        write(*,*)'the data in reflection point before reflection'
        write(*,*)'z=',u(1),'r=',u(2),'phi=',u(3)
        write(*,*)'cnz=',u(4),'cnr=',u(5),'cm=',u(6)
        bmod=b(u(1),u(2),u(3))    
        write(*,*)'cnpar=',(u(4)*bz+u(5)*br+u(6)*bphi/u(2))/bmod
        write(*,*)'cn2=',u(4)**2+u(5)**2+(u(6)/u(2))**2
        write(*,*)'cnper=',dsqrt((u(4)**2+u(5)**2+(u(6)/u(2))**2)-
     &                 ((u(4)*bz+u(5)*br+u(6)*bphi/u(2))/bmod)**2)
        write(*,*)'after reflection'

        u(1)= z_ref
        u(2)= r_ref
        u(3)= phi_ref

        write(*,*)'z_ref=',z_ref,'r_ref=',z_ref,'pfi_ref=',phi_ref
        write(*,*)'cnzref=',cnzref,'cnrref=',cnrref,'cmref=',cmref
        write(*,*)'cnpar=',(cnzref*bz+cnrref*br+u(6)*bphi/u(2))/bmod
        write(*,*)'cn2=',cnzref**2+cnrref**2+(u(6)/u(2))**2
        write(*,*)'cn2=',cnzref**2+cnrref**2+(cmref/u(2))**2
        write(*,*)'cnper=',dsqrt((cnzref**2+cnrref**2+(u(6)/u(2))**2)-
     &                 ((cnzref*bz+cnrref*br+u(6)*bphi/u(2))/bmod)**2)
        write(*,*)'cnper=',dsqrt((cnzref**2+cnrref**2+(cmref/u(2))**2)-
     &                 ((cnzref*bz+cnrref*br+cmref*bphi/u(2))/bmod)**2)

      endif

      if ((i_ox.eq.1).and.(iflref.eq.1)) then
c--------calculate the EC cone vertex coordinates
c        for the optimal OX mode conversion
           bmod=b(u(1),u(2),u(3))     
           b_z0=bz
           b_r0=br
           b_phi0=bphi
           dpsi_dz=dpdzd
           dpsi_dr=dpdrd
           z0=u(1)
           r0=u(2)
           phi0=u(3)
           cnz0=u(4)
           cnr0=u(5)
           cm0=u(6)
           cnphi0=cm0/r0              
cyup           write(*,*)'output i_ox=1'
cyup           write(*,*)'z0,r0,phi0',z0,r0,phi0
cyup           write(*,*)'cnz0,cnr0,cm0,cnphi0',cnz0,cnr0,cm0,cnphi0
           
c------test
c           write(*,*)'bz,br,bphi',bz,br,bphi
cyup            write(*,*)'output.f before call antenna_vertex'
           cnpar_test=(cnz0*bz+cnr0*br+cnphi0*bphi)/bmod
cyup           write(*,*)'output cnpar_test',cnpar_test
           cn0_s=cnz0**2+cnr0**2+cnphi0**2
           cnper_test=dsqrt(cn0_s-cnpar_test**2)
c           write(*,*)'dsqrt(cn0_s)',dsqrt(cn0_s)
cyup           write(*,*)'output cnper_test',cnper_test
           sin_phi=dsin(phi0)
           cos_phi=dcos(phi0)
           cnx=cnr0*cos_phi-cnphi0*sin_phi
           cny=cnr0*sin_phi+cnphi0*cos_phi
           cnz=cnz0  
cyup           write(*,*)'cnx,cny,cnz',cnx,cny,cnz
           call ninit_ec(z0,r0,phi0,cnx,cny,cnz,cnteta,cnphi0)
cyup           write(*,*)'cnteta,cnphi0',cnteta,cnphi0
cyup           write(*,*)'cnrho',dsqrt(cn0_s-cnteta**2-cnphi0**2) 
           ioxm_loc=ioxm
           ioxm=1  
           call cninit12(z0,r0,phi0,cnpar_test,cnteta,cnphi0,
     1                  cnz_t,cnr_t,cm_t,iraystop)
           ioxm=ioxm_loc 
cyup           write(*,*)'output.f after cninit12 cnz_t,cnr_t,cm_t',
cyup     &     cnz_t,cnr_t,cm_t
c----------endtest
           icone=icone_loc
           call antenna_vertex(u(2),u(3),u(1),-u(5),-u(6)/u(2),-u(4),
     &     b_r0,b_phi0,b_z0,dpsi_dz,dpsi_dr,
     &     r_st,phi_st,z_st,alpha_st,beta_st,icone)

c----------put the data into cone.i
           r_st_ox=r_st
           phi_st_ox=phi_st
           z_st_ox=z_st
           alpha_st_ox=alpha_st
           beta_st_ox=beta_st

cyup           write(*,*)'oxb r_st,phi_st,z_st,alpha_st,beta_st',
cyup     &     r_st,phi_st,z_st,alpha_st,beta_st
 
      endif ! ((i_ox.eq.1).and.(iflref.eq.1)  

      u(4)=cnzref
      u(5)=cnrref
      u(6)=cmref

cSAP081010
      if (iflref.eq.1) call prep3d(t,u,deru,iraystop)
c----------------------------------------------------------
c     call b() to calculate the small radius rho (inside b())
c     the resulting will be in common block  one.i
c-----------------------------------------------------------
      bmod=b(u(1),u(2),u(3)) 
c      write(*,*)'output iflref=',iflref
      if (iflref.eq.1) then         
c-----------------------------------------------------------
c        Reflecting procedure have worked
c        The scattering of N_perp after the reflection

c         write(*,*)'output.f iscat=',iscat,'scatd(0)=',scatd(0)

         if ((iscat.eq.1).and.(scatd(0).ne.0.d0)) then
             bmod=b(u(1),u(2),u(3)) 
c             write(*,*)'output u(4-6) before scattering '
c             write(*,*)u(4),u(5),u(6)
             call scatperp(bz,br,bphi,u(2),u(4),u(5),u(6),scatd(0))
c             write(*,*)'output u(4-6) after scattering'
c             write(*,*)u(4),u(5),u(6)
         endif
c-----------------------------------------------------------
         iflagh=1
      else
c-----------------------------------------------------------
c         The N_perp scattering inside  the plasma
c      	 write(*,*)'output iscat=',iscat
         if (iscat.eq.1) then
           do i=1,nscat_n                 
c             write(*,*)'rhooldsc,rho,i,rhoscat(i)'  
c             write(*,*)rhooldsc,rho,i,rhoscat(i)  
             if ((rhooldsc-rhoscat(i))*(rho-rhoscat(i)).lt.0.d0) then
c               write(*,*)'rho,rhooldsc,i,rhoscat(i)',
c     1         rho,rhooldsc,i,rhoscat(i)        
c               write(*,*)'output u(4-6) before interior scattering '
c               write(*,*)u(4),u(5),u(6)
               call scatperp(bz,br,bphi,u(2),u(4),u(5),u(6),scatd(i))
c               write(*,*)'output u(4-6) after interior scattering'
c               write(*,*)u(4),u(5),u(6)
             endif
	   enddo
         endif

c---------------------------------------------------------------
      end if

      rhooldsc=rho
cSAP101217
       if(iscat_lh_Nicola.eq.1) then
c--------scattering using Nicola subroutine
         s_total_new=s_total

cSAPO120603
cyup         write(*,*)'output.f bef sc u',u
cyup         write(*,*)'s_total_new, s_total_old',s_total_new,s_total_old

         call lhscattering(u, s_total_new, s_total_old)

cSAP120603
cyup         write(*,*)'output.f aft sc u',u
cyup         write(*,*)'s_total_new, s_total_old',s_total_new,s_total_old

       endif

c      write(*,*)'output irefl,ireflm',irefl,ireflm 
      if (irefl.ge.ireflm) then
cyup         write(*,*)'****irefl.ge.ireflm u(2)',u(2)
cyup         write(*,*)'i_output_data ',i_output_data 
         if (i_output_data.eq.0) then
cyup             write(*,*)'output.f i_output_data.eq.0 befor prep3d'
             call prep3d(t,u,deru,iraystop)
cyup           write(*,*)'output.f i_output_data.eq.0 aft prep3d iraystop',
cyup     &                iraystop
         endif
!nme         if(iout3d.eq.'enable') then
c----------writing data to output mnemonic.txt file
!nme	   call write3d
!nme         endif

 	 iraystop=1     
         iray_status_one_ray=3  
      else
 	 iraystop=0

      end if

cyup      write(*,*)'output.f after irefl.ge.ireflm iraystop,t',iraystop,t

 120  continue
c----------------------------------------------------
c     for control phi_deviation for one time step
c     write(*,*)'in output phiold,u(3)',phiold,u(3)
      phiold=u(3)  ! phiold is in common write.i
c----------------------------------------------------
c      write(*,*)'end outpt iraystop',iraystop
      return
      end

      subroutine correct2(cnpar,cnper,
     *cnz,cnr,cm,r,bz,br,bphi,bmod,cnz1,cnr1)

c      subroutine correct2(xe,ye,te_kev,tpope,vflowe,cnpar,id,
c     *ihermloc,accurcy,naccurc,cnprim,cnper,
c     *cnz,cnr,cm,r,bz,br,bphi,bmod,cnz1,cnr1)
c-----------------------------------------------------------
c     it calculates new values of N_perp coordinate
c     to get the Hamiltonian conservation
c     from the solution Nperp(N_parallel)
      implicit none
c-----input
      double precision cnpar,cnper,
     *cnz,cnr,cm,r,bz,br,bphi,bmod

c-----output
      double precision cnr1,cnz1
c-----local
      double precision cnper2,cnpar2,d1,d2,a0,a1,a2,det,
     *delnm,delnp,cnz1p,cnz1m,cnr1p,cnr1m,cnperold
      integer inewn
c     inewn=0 ! There were used the old values(cnz,cnr) inside 
c               this subbroutine
c          =1 ! There were used the new corrected
c               values(cnz1,cnr1) inside this subroutine

cfor test
      double precision oldnpol,newnpol,cn2new,bpol,bplnpl,
     *cosnb_pl
      
c------------------------------------------------------------  
c     vector{cnper}={nz-(bz/bmod)N_par,
c                    nr-(br/bmod)N_par
c                    nphi-(bphi/bmod)N_par}
c     we will determine the new (cnz1,cnr1) from the following system
c     cnr1*br+cnz1*bz=cnpar*bmod-(cm/r)*bphi=d1
c     cnr1**2+cnz1**2=cnpar**2+cnper**2-(cm/r)**2=d2
c-----------------------------------------------------------------
      cnper2=cnper*cnper
      cnpar2=cnpar*cnpar    
      d1=cnpar*bmod-(cm/r)*bphi
      d2=cnpar2+cnper2-(cm/r)**2

ctest
c      write(*,*)'output_correct2 input cnpar',cnpar
c      write(*,*)'output_correct2 input cnz,cnr,cm,r',cnz,cnr,cm,r
c      write(*,*)'output_correct2 input bz,br,bphi,bmod',
c     1bz,br,bphi,bmod
c      write(*,*)'output_correct2 input n*b/bmode',
c     1(cnz*bz+cnr*br+(cm/r)*bphi)/bmod
ctest_end

ctest
      oldnpol=dsqrt(cnz*cnz+cnr*cnr)
      cn2new=cnper2+cnpar2
      newnpol=dsqrt(cn2new-(cm/r)**2)
      bpol=dsqrt(bz*bz+br*br)/bmod
      bplnpl=cnpar-(cm/r)*(bphi/bmod)
c      write(*,*)'correct2 oldnpol,newnpol',oldnpol,newnpol
      cosnb_pl=bplnpl/newnpol
c      write(*,*)'correct2 bpol,cosnb_pl,bplnpl',
c     *bpol,cosnb_pl,bplnpl
ctestend


      inewn=0
      if (d2.lt.0.d0) then
         write(*,*)'correct2 d2=cnpar2+cnper2-(cm/r)**2<0
     *   correction was not made'       
         goto 10
      else
         if ((br*br+bz*bz).eq.0.d0) then
c           write(*,*)'correct2 (br*br+bz*bz).eq.0.d0'
            cnperold=dsqrt(cnz*cnz+cnr*cnr)
c           write(*,*)'coorect2 cnper,cnperold',cnper,cnperold
            if(cnperold.eq.0.d0) then
               cnz1=cnper/dsqrt(2.d0)
               cnr1=cnper/dsqrt(2.d0)
            else 
               cnz1=cnz*cnper/cnperold
               cnr1=cnr*cnper/cnperold
            endif
            inewn=1
c            write(*,*)'correct2 cnz1,cnr1',cnz1,cnr1
            goto 10
         endif

         if(br.ne.0.d0) then
c----------equation a2*cnz1**2-2*a1*cnz1+a0=0

           a2=1.d0+(bz/br)**2
           a1=(bz*d1)/(br*br)
           a0=(d1/br)**2-d2 
           det=a1*a1-a2*a0

           if(det.lt.0.d0) then
             write(*,*)'correct2 br.ne.0 det<0
     *       correction in correct2 was not made'
             goto 11
           else
             cnz1m=(a1-dsqrt(det))/a2
             cnz1p=(a1+dsqrt(det))/a2
             cnr1m=(d1-bz*cnz1m)/br
             cnr1p=(d1-bz*cnz1p)/br
             inewn=1        
           endif
         endif
 11      continue

         if(bz.ne.0.d0) then
c----------equation a2*cnr1**2-2*a1*cnr1+a0=0
           a2=1.d0+(br/bz)**2
           a1=(br*d1)/(bz*bz)
           a0=(d1/bz)**2-d2 
           det=a1*a1-a2*a0

           if(det.lt.0.d0) then
             write(*,*)'correct2 bz.ne.0 det<0
     *       correction in correct2 was not made'
             goto 10
           else
             cnr1m=(a1-dsqrt(det))/a2
             cnr1p=(a1+dsqrt(det))/a2
             cnz1m=(d1-br*cnr1m)/bz
             cnz1p=(d1-br*cnr1p)/bz
             inewn=1	     
           endif
         endif

      endif

ctest
c      write(*,*)'output_correct2 cnr1p*br+cnz1p*bz',cnr1p*br+cnz1p*bz
c      write(*,*)'output_correct2 cnr1m*br+cnz1m*bz',cnr1m*br+cnz1m*bz
c      write(*,*)'output_correct2 d1',d1
c      write(*,*)'output_correct2 cnr1p**2+cnz1p**2,d2',
c     1cnr1p**2+cnz1p**2,d2  
c      write(*,*)'output_correct2 input cnpar',cnpar
c      write(*,*)'output_correct2 output cnz1p,cnr1p,cm,r',
c     1cnz1p,cnr1p,cm,r
c      write(*,*)'output_correct2 output cnz1m,cnr1m,cm,r',
c     1cnz1m,cnr1m,cm,r
c      write(*,*)'output_correct2 output bz,br,bphi,bmod',
c     1bz,br,bphi,bmod
c      write(*,*)'output_correct2 output n1p*b/bmode',
c     1(cnz1p*bz+cnr1p*br+(cm/r)*bphi)/bmod
c      write(*,*)'output_correct2 output n1m*b/bmode',
c     1(cnz1m*bz+cnr1m*br+(cm/r)*bphi)/bmod
ctest_end

      delnm=(cnz-cnz1m)**2+(cnr-cnr1m)**2
      delnp=(cnz-cnz1p)**2+(cnr-cnr1p)**2
c      write(*,*)'output_correct2 delnm,delnp',delnm,delnp
      if (delnm.le.delnp) then
        cnz1=cnz1m
        cnr1=cnr1m
      else
        cnz1=cnz1p
        cnr1=cnr1p
      endif

 10   continue

      if (inewn.eq.0) then
c-------the correction was not made
        cnz1=cnz
        cnr1=cnr
      endif

c      write(*,*)'solvnper end inewn cnz1,cnr1',inewn,cnz1,cnr1     
ctest
c      write(*,*)'output_correct2 input cnpar',cnpar
c      write(*,*)'output_correct2 output cnz1,cnr1,cm,r',cnz1,cnr1,cm,r
c      write(*,*)'output_correct2 output bz,br,bphi,bmod',
c     1bz,br,bphi,bmod
c      write(*,*)'output_correct2 output n1*b/bmod',
c     1(cnz1*bz+cnr1*br+(cm/r)*bphi)/bmod
ctest_end
c

      return
      end


c-------------------------------------------------------------
c     calculation n_perp(n_par) for Mazzucato and Forest codes
c     It uses the input cnper=Re(n_perp) for the first iteration
c     It calculates cnprim=Im(N_perp) for forest case (id=6)
      subroutine solvnperp(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .vflow_ar,
     .cnpar,id,
     .ihermloc,accurcy,naccurc,cnprim,cnper)
      implicit none
c-----input
      integer nbulk
      double precision dmas(*),x_ar(*),y_ar(*),te_ev_ar(*),tpop_ar(*),
     .vflow_ar(*),
     .cnpar,cnper,cnprim,accurcy
      integer id,ihermloc,naccurc
c-----output 
c     cnper,cnprim   
c-----local
      double precision step,hamr,hamrp,hamrm,dhamrdnr
      double precision cnperp,cnperm,dcnper,hamold,cnpernew,cnpernw1
      double complex chamilt
      integer iter
c     calculation cnper from the solution of the equation
c     dhamr/d_nper*delnper=-hamr

      
      cnprim=0.d0
      step=1.d-7
      iter=0
c      write(*,*)'solvnper xe,ye,te_kev',xe,ye,te_kev
c      write(*,*)'solvnper cnpar,cnper,id,naccurc',cnpar,cnper,id,naccurc
      cnpernew=cnper !initial value
      cnpernw1=cnpernew

 10   continue

c      write(*,*)'1 solvnper cnpar,cnper,id',cnpar,cnper,id
c      write(*,*)'solvnperp bef dispfun1 cnpar,iter,cnpernew',
c     *cnpar,iter,cnpernew

      call dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     &vflow_ar,
     *cnpar,cnpernew,cnprim,id,ihermloc,
     *chamilt)
      hamr=dreal(chamilt)

c      write(*,*)'solvnperp after dispfun1 iter,hamr,cnpernw1,cnpernew'
c     *,iter,hamr,cnpernw1,cnpernew

      if (iter.eq.0) then 
        hamold=hamr
      else        
        if (dabs(hamold).le.dabs(hamr)) then
          write(*,*)'solvnperp dabs(hamold)<dabs(hamr) iter=',iter
          write(*,*)'solution Nperp not found hamold,hamr',hamold,hamr
          write(*,*)'cnper=cnpernw1',cnpernw1
          cnper=cnpernw1     
          goto 20
        else
          cnper=cnpernew
          cnpernw1=cnpernew
        endif
      endif

      if (dabs(hamr).le.accurcy) goto 20

      cnperp=cnper+step*cnper
      call dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     &vflow_ar,
     *cnpar,cnperp,cnprim,id,ihermloc,chamilt)
      hamrp=dreal(chamilt)
      
      cnperm=cnper-step*cnper
      call dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     &vflow_ar,
     *cnpar,cnperm,cnprim,id,ihermloc,chamilt)
      hamrm=dreal(chamilt)
         
      dhamrdnr=(hamrp-hamrm)/(2.d0*step*cnper)
       
c      write(*,*)'dhamrdnr,hamr',dhamrdnr,hamr
      if (dhamrdnr.eq.0.d0) then
        write(*,*)'solvnperp dhamdnnr=0'
        goto 20
      else    
        dcnper=-hamr/dhamrdnr
        cnpernew=cnper+dcnper
        iter=iter+1

        if(cnpernew.lt.0.d0) then
          write(*,*)'output.f in solvnper cnpernew<0'
          cnpernew=cnper
          goto 20
        endif

        if (iter.gt.naccurc) then
           write(*,*)'solvnper iter>naccurc'
           cnper=cnpernew
           goto 20
        endif

      endif
      hamold=hamr
      goto 10

 20   continue
c      write(*,*)'solvnper bef end cnper,iter',cnper,iter
      return
      end

      subroutine dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .vflow_ar,
     .cnpar,cnper,cnprim,id,
     .ihermloc,chamilt)
c     It calculates the complex dispersion function chamilt
c     for different tensors (NOw only for id=4,5 and 6)
c     id=1
c     id=2
c     id=3 Appleton-Hartry
c     id=4 Mazzucato (hermitian part)
c     id=5 Mazzucato total
c     id=6 Forest code
c-----it calculates Imaginary part cnprim=Im(N_perp) for Forest code
      implicit none
c-----input
      integer nbulk
      double precision dmas(*),x_ar(*),y_ar(*),te_ev_ar(*),
     .tpop_ar(*),vflow_ar(*)      
      double precision xe,ye,te_kev,tpope,vflowe
      double precision cnpar,cnper,cnprim
      integer id,ihermloc
c-----output
      double complex chamilt
c     cnprim : is a solution for forest code
c-----external dhot,complx1
      double complex dhot_sum       
c-----local
      double precision hamr,hami,mode
      double precision mox2de,mu
      double precision cnper2,te_ev,cnparp,cnperp
      double complex cpp,sol
      double complex ceps(3,3)
   
      if (((id.eq.1).or.(id.eq.2)).or.(id.eq.3)) then
c        call tensrcld(z,r,phi)
         write(*,*)'dispfun id=',id,'now id can be id=4,=5 or 6'
         stop         
      endif
      
      if ((id.eq.4). or.(id.eq.5)) then
c--------Mazzucato dispersion function
         xe=x_ar(1)
         ye=y_ar(1)
         te_kev=te_ev_ar(1)*1.d-3
         mode=dfloat(1) ! it not used
         mu=512.44d00/te_kev ! te (kev)
         cpp=dcmplx(cnper,0.0d0)
         cnper2=cnper*cnper
         sol=cpp*cpp
c         write(*,*)'dispfun_output mode,xe,ye',mode,xe,ye
c         write(*,*)'dispfun_output mu,cnpar,cnper2',mu,cnpar,cnper2
c         write(*,*)'dispfun_output sol,ihermloc',sol,ihermloc
         
         call complx1(mode,xe,ye,mu,cnpar,cnper2,sol,ihermloc,chamilt)
c         write(*,*)'in dispfun chamilt',chamilt
         goto 100
      endif ! id =4,5 Mazzucato function

      if (id.eq.6) then
c--------Hot dispersion function
c        electron + ions non-relativistic plasma
                 
         cnparp=cnpar
         cnperp=cnper

         hamr=dreal(dhot_sum(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .   vflow_ar,cnpar,cnper,1,ceps))
         hami=0.d0
         chamilt=dcmplx(hamr,hami)
         goto 100
      endif !id=6 Hot non-relativistic dispersion

 100  continue
      return
      end

      subroutine rdispfun(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .vflow_ar,
     .cnpar,cnper,cnprim,id,
     .ihermloc,hamilt)

c-----it calculates Imaginary part cnprim=Im(N_perp) for Forest code
      implicit none
c-----input     
      integer nbulk
      double precision dmas(*),x_ar(*),y_ar(*),te_ev_ar(*),
     .tpop_ar(*),vflow_ar(*) !cm/sec
      double precision cnpar,cnper,cnprim
      integer id,ihermloc
      external dispfun_output
      
c-----output 
      double precision hamilt
c     as output cnprim is a solution for the forest code
        
c-----local
      double complex chamilt
     
      call dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     &vflow_ar,
     .cnpar,cnper,cnprim,id,
     .ihermloc,chamilt)

      hamilt=dreal(chamilt)

      return
      end


      subroutine mphynper(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .vflow_ar,cnpar,cnper,cnprim,id,ihermloc,
     .y_min,y_max,m_y,nper_min,nper_max,m_nperp)
c-----this subroutine calculates the table Real part of dispersion function
c     Re(hamiltonian)(nperp,Y_k) for the given nperp,Yc     
c     k=1: Y_1 =Y_e
c     Subroutine creates the table (Y_j1,N_perp_j2,hamitonian(Y_j1,N_perp_j2))
c     j1=1,...m_y, j2=1,...m_nperp
c     y_min=< Y_j1 =<y_max, nper_min=< N_perp_j2 =<nper_max 
c     The table for Real(hamiltonian) will bi in hamr_ynp.dat file
c     The table for Image(hamiltonian) will bi in hami_ynp.dat file
     
      implicit none
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c-----input
      integer nbulk
      double precision dmas(*),x_ar(*),y_ar(*),te_ev_ar(*),
     .tpop_ar(*),vflow_ar(*)      
      double precision xe,ye,te_kev,tpope,vflowe
      double precision cnpar,cnper,cnprim
      integer id,ihermloc
c-----the data for the map mesh
      integer m_y,m_nperp
      double precision y_min,y_max,nper_min,nper_max
c-----locals 
      integer j1,j2
      double precision hy,hnperp,yloc,nperploc
      double complex chamilt

      if(myrank.ne.0) return
      
      open(10,file='hamrhynp.dat') !hot plasma
      open(20,file='hamihinp.dat') 
      write(10,*)'y nper hamr' 
      write(20,*)'y nper hami'
      
      open(50,file='hamrmynp.dat') !mazzucato plasma
      open(60,file='hamiminp.dat')
      write(50,*)'y nper hamr' 
      write(60,*)'y nper hami'

 1    format (3(' ',1pe15.7))

      hy=(y_max-y_min)/(m_y-1.d0)
      hnperp=(nper_max-nper_min)/(m_nperp-1.d0)

      do j1=1,m_y
        yloc=y_min+hy*(j1-1)
        y_ar(1)=yloc
        do j2=1,m_nperp
          nperploc=nper_min+hnperp*(j2-1)

c---------hot plasma
cSAP110316
c          call dispfun(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
          call dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .    vflow_ar,
     .    cnpar,nperploc,cnprim,id,
     .    ihermloc,chamilt)         
          write(10,1)yloc,nperploc,dreal(chamilt)
          write(20,1)yloc,nperploc,dimag(chamilt)


c---------mazzucato plasma
          nperploc=nper_min+hnperp*(j2-1)
          y_ar(1)=-yloc
cSAP110316
c          call dispfun(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,    
          call dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .    vflow_ar,
     .    cnpar,nperploc,cnprim,4,
     .    ihermloc,chamilt)         
          write(50,1)yloc,nperploc,dreal(chamilt)
          write(60,1)yloc,nperploc,dimag(chamilt)


        enddo
      enddo

      close(10)
      close(20)

      close(50)
      close(60)

      return
      end


      subroutine switch_da(y,del_y,id,iabsorp,idswitch,iabswitch)
c-----switches on(off) the dispersion function and the absorption near the cyclotron
c     resonance points 
c     input:
c     y    (omega_c/omega) algebraic
c     del_y the distance in Y from the resonace to swith on dispersion and absorption
c     id        the dispersion far from reconace
c     iabsorp   the absorpsion far from the resonance
c     idswitch  the dispersion near resonance
c     iabswitch  the absorpsion near resonance
c     output: as a local save variable
c     iswitch !=0 before the first change of  dispersion
c             !=1 after first changing of the dispersion
      implicit none
c-----input
      double precision y,del_y !
      integer id,iabsorp,idswitch,iabswitch     
c-----local
      integer iy
      integer iswitch,id_old,iabs_old 
      save iswitch,id_old,iabs_old

      iy=nint(1.d0/y)
      write(*,*)'1 iswitch',iswitch,'y',y,'iy',iy,'del_y',del_y
c      if((dabs(y-1.d0/dfloat(iy)).le.del_y).and.(iswitch.ne.1)) then
      if((dabs(1.d0/y-dfloat(iy)).le.del_y).and.(iswitch.ne.1)) then

c------ the change of the dispersion function and absorption near the cyclotron resonance
        iswitch=1  
        id_old=id
        id=idswitch       ! new
        iabs_old=iabsorp
        iabsorp=iabswitch ! new 
        write(*,*)'in output:  switch_da'
        write(*,*)'id=',id,'iabsorp=',iabsorp
      endif
      write(*,*)'2 iswitch',iswitch,'y',y
c      if ((dabs(y-1.d0/dfloat(iy)).ge.del_y).and.(iswitch.eq.1)) then
      if ((dabs(1.d0/y-dfloat(iy)).ge.del_y).and.(iswitch.eq.1)) then

c------ the back switch of the dispersion function after ec resonance 
        iswitch=0
        id=id_old
        iabsorp=iabs_old
        write(*,*)'in output: switch_da back switch of dispersion'
        write(*,*)'id=',id,'iabsorp=',iabsorp
      endif         
      write(*,*)'output.f 3 iswitch,id',iswitch,id
      return
      end

      subroutine switch_d_ono(r,z,phi,cnz,cnr,cm,del_y0,y0_cr,
     &id,idswitch)
c-----switches on(off) the Ono dispersion function (id=8) 
c     to cold plasma dispersion (id=2) near the critical point, 
c     where Z function has a maximum at y0=y0_c: 
c     |y0-y0_c| < del_y0 . Here y0_cr/(N_parallel*V_te),
c     y0_c is about 0.94,del_y is about 0.01  
c     input:
c     z,r,phi,cnz,cnr,cm -coordinates
c     y0_cr  the critical value of y0 where Z haz a maximum
c     del_y0 the distance in y0 from the critical point to swith on dispersion
c     id        the dispersion far from critical point
c     idswitch  the dispersion near resonance
c     output: as a local save variable
c     iswitch !=0 before the first change of  dispersion
c             !=1 after first changing of the dispersion
c     id      new or old id
      implicit none
c-----input
      double precision r,z,phi,cnz,cnr,cm,y0_cr,del_y0 
      integer id,idswitch     
c-----local
      integer iy
      integer iswitch,id_old,sign 
      double precision bmod,cnt,gam,cnpar,t_eV,vt,c,k4,y0
c-----external
      double precision b,gamma1,cn,tempe

      save iswitch,id_old

      bmod=b(z,r,phi)
      cnt=cn(r,cnz,cnr,cm)
      gam=gamma1(z,r,phi,cnz,cnr,cm)
      cnpar=cnt*dcos(gam)
c      write(*,*)'switch_d_ono cnt,gam,cnpar',cnt,gam,cnpar

      c = 3.0d10                    !speed of light [cm/sec]
      k4 = 4.19d7                   !constant for electron thermal velocity
      t_eV=tempe(z,r,phi,1)*1.d3    !electron temperature in eV  
      vt =k4*dsqrt(2.0d0*t_eV) !elevctron vt= sqrt(2.0*kT_e[eV]/mass_e) cm/sec 

      y0=c/(cnpar*vt)

      write(*,*)'1 iswitch Ono ',iswitch,' y0_cr,y0,del_y0',
     &y0_cr,y0,del_y0

      if (y0.gt.0.d0)then
         sign=1.d0
      else
         sign=-1.d0
      endif

      write(*,*)'y0,sign,dabs(y0_cr),dabs(y0-sign*dabs(y0_cr)),del_y0',
     &y0,sign,dabs(y0_cr),dabs(y0-sign*dabs(y0_cr)),del_y0
      if((dabs(y0-sign*dabs(y0_cr)).le.del_y0).and.(iswitch.ne.1)) then
c------ the change of the dispersion function near the critical point
        iswitch=1  
        id_old=id
        id=idswitch       ! new

        write(*,*)'in output:  switch_d_ono'
        write(*,*)'id=',id
      endif
      write(*,*)'2 iswitch',iswitch,'y0',y0

      if((dabs(y0-sign*dabs(y0_cr)).ge.del_y0).and.(iswitch.eq.1)) then
c------ the back switch of the dispersion function after the critical point
        iswitch=0
        id=id_old
        write(*,*)'in output: switch_d_ono back switch of dispersion'
        write(*,*)'id=',id
      endif         
      write(*,*)'output.f 3 Ono iswitch,id',iswitch,id
      return
      end

      subroutine correct3(cnpar,cnper,cnz,cnr,cm,z,r,phi,
     .eps,itermax,cnz1,cnr1,r1)
c     calculates the corrected values 
c     cnz1=cnz+delnz,cnr1=cnr+delnr,r1=r+delr from the problem
c     J=delnz**2+delnr**2_delr**2
c     min(J)
c     N_par(Nz+delnz,Nr+delnz,r+delr)=cnpar
c     N_perp(Nz+delnz,Nr+delnr,r+delnr)=cnperp

c-----input
      double precision cnpar,cnper,cnz,cnr,cm,z,r,phi
      double precision eps      !the accuracy of iterations
      integer itermax           !max number of steps in iterations
c-----output
      double precision cnz1,cnr1,r1
c-----locals
      double precision dr,dnz,dnr
      double precision f1,df1ddnz,df1ddnr,df1ddr,
     .f2,df2ddnz,df2ddnr,df2ddr       
      integer iter
      double precision derg(3,3),g123ar(3),det,det1,det2,det3,norma
      double precision bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl

      iter=0
      cnz1=cnz
      cnr1=cnr
      r1=r
      dnz=0.d0
      dnz=0.d0
      dr=0.d0
      
 10   continue
      write(*,*)'output correct3 z,r,phi',z,r,phi
      write(*,*)'output correct3 cnz,cnr,cm',cnz,cnr,cm

      call bcomp(z,r,phi,bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl)
  
c      write(*,*)'correct3 bzl,brl,bphil,bmodl',bzl,brl,bphil,bmodl     
c      write(*,*)'correct3 dbzdzl,dbzdrl,dbrdzl,dbrdrl',
c     .dbzdzl,dbzdrl,dbrdzl,dbrdrl
c      write(*,*)'correct3dbphdzl,dbphdrl,dbmdzl,dbmdrl',
c     .dbphdzl,dbphdrl,dbmdzl,dbmdrl
      call g123(cnz1,cnr1,cm,z,r1,phi,dnr,dnz,dr,cnpar,cnper,
     .bzl,brl,bphil,bmodl,dbzdrl,dbrdrl,dbphdrl,dbmdrl,
     .g123ar)
c      write(*,*)'output g123ar',g123ar
      norma=g123ar(1)**2+g123ar(1)**2+g123ar(1)**2
      write(*,*)'output correct3 norma',norma
c      if(norma.lt.eps) goto 10
      if(norma.lt.eps) goto 20

      call drvg123(cnz1,cnr1,cm,z,r1,phi,dnr,dnz,dr,cnpar,cnper,
     .bzl,brl,bphil,bmodl,dbzdrl,dbrdrl,dbphdrl,dbmdrl,
     .derg)
      
      det=derg(1,1)*(derg(2,2)*derg(3,3)-derg(2,3)*derg(3,2))-
     .derg(1,2)*(derg(2,1)*derg(3,3)-derg(2,3)*derg(3,1))+
     .derg(1,3)*(derg(2,1)*derg(3,2)-derg(2,2)*derg(3,1))
      
      det1=-g123ar(1)*(derg(2,2)*derg(3,3)-derg(2,3)*derg(3,2))+
     .g123ar(2)*(derg(1,2)*derg(3,3)-derg(1,3)*derg(3,2))-
     .g123ar(3)*(derg(1,2)*derg(2,3)-derg(1,3)*derg(2,2))

      det2=g123ar(1)*(derg(2,1)*derg(3,3)-derg(2,3)*derg(3,1))-
     .g123ar(2)*(derg(1,1)*derg(3,3)-derg(1,3)*derg(3,1))+
     .g123ar(3)*(derg(1,1)*derg(2,3)-derg(2,1)*derg(1,3))

      det3=-g123ar(1)*(derg(2,1)*derg(3,2)-derg(2,2)*derg(3,1))+
     .g123ar(2)*(derg(1,1)*derg(3,2)-derg(1,2)*derg(3,1))-
     .g123ar(3)*(derg(1,1)*derg(2,2)-derg(1,2)*derg(2,1))

      if (det.eq.0.d0) then 
        write(*,*)' in output.d correct3 det=0'
        stop
      endif

      dnz=det1/det
      dnr=det2/det
      dr=det3/det
      iter=iter+1

      if(iter.gt.itermax) then     
        write(*,*)'in correct3 iter>itermax'
        goto 20
      endif

      goto 10
      
 20   continue

      cnz1=cnz+dnz
      cnr1=cnr+dnr
      r1=r+dr

      return
      end

      double precision function f1(cnz,cnr,cm,r,bz,br,bphi,
     .bmod,cnpar)
      implicit none
      double precision cnz,cnr,cm,r,bz,br,bphi,bmod,cnpar
      f1=cnz*bz+cnr*br+cm*bphi/r-cnpar*bmod
      return
      end

      double precision function f2(cnz,cnr,cm,r,cn)
      implicit none
c-----cn**2=Npar**2+Nper**2  
      double precision cnz,cnr,cm,r,cn
      f2=cnz**2+cnr**2+(cm/r)**2-cn**2
      return
      end

      double precision function df1dnz(bz)
      implicit none
      double precision bz
      df1dnz=bz
      return
      end

      double precision function df1dnr(br)
      implicit none
      double precision br
      df1dnr=br
      return
      end

      double precision function df1dr(cnz,cnr,cm,r,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)
      implicit none
      double precision cnz,cnr,cm,r,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr 
      df1dr=cnz*dbzdr+cnr*dbrdr-cm*bphi/r**2+
     .cm*dbphdr/r-cnpar*dbmdr
      return
      end

      double precision function df2dnz(cnz)
      implicit none
      double precision cnz,cnpar,bz,bmod
      df2dnz=2*cnz
      return
      end

      double precision function df2dnr(cnr)
      implicit none
      double precision cnr,cnpar,br,bmod
      df2dnr=2*cnr
      return
      end

      double precision function df2dr(cm,r)
      implicit none
      double precision cm,r
      df2dr=-cm*cm/r**3
      return
      end
    
      subroutine al12(dnz,dnr,df1dnz,df1dnr,df2dnz,df2dnr,al1,al2)
c     calculates al1 and al2
      implicit none
c-----input
      double precision dnz,dnr,df1dnz,df1dnr,df2dnz,df2dnr
c-----output
      double precision al1,al2
c-----local 
      double precision det,det1,det2

      det=df1dnz*df2dnr-df1dnr*df2dnz
      det1=2.d0*(dnz*df2dnr-dnr*df2dnz)
      det2=2.d0*(-dnz*df1dnr+dnr*df1dnz)

      if (det.eq.0.d0) then
         write(*,*)'output.f in al12 det=0 stop'
         stop
      else
         al1=det1/det
         al2=det2/det
      endif

      return
      end

      double precision function dlagrdr(dr,al1,al2,df1dr,df2dr)
c     calculates derivative d(Lagrange function)/d(deltar) 
      implicit none
c-----input
      double precision dr,al1,al2,df1dr,df2dr 
      dlagrdr=2.d0*dr+al1*df1dr+al2*df2dr
      return
      end
  
      subroutine g123(cnz,cnr,cm,z,r,phi,dnr,dnz,dr,cnpar,cnper,
     .bz,br,bphi,bmod,dbzdr,dbrdr,dbphdr,dbmdr,
     .g123ar)
c     input magnetic field should be done in
c     (r+delr) point
c     calculates
c     g123ar(1)=f1(Nz+delNz,Nr+delNr,r+delr)
c     g123ar(2)=f2(Nz+delNz,Nr+delNr,r+delr)
c     g12ar(3)=2*delr+al1(delNz,delNr,delr)*df1(delNz,delNr,delr)/d(delr)+
c               al2(delNz,delNr,delr)*df2(delNz,delNr,delr)/d(delr)

      implicit none
c-----input 
      double precision cnz,cnr,cm,z,r,phi,dnz,dnr,dr,cnpar,cnper,
     .bz,br,bphi,bmod,dbzdr,dbrdr,dbphdr,dbmdr
c-----external
      double precision dlagrdr,df1dr,df2dr,f1,f2,
     .df1dnz,df1dnr,df2dnz,df2dnr
c-----output
      double precision g123ar(3)
c-----local
      double precision al1,al2,df1drl,df2drl,cn,cnzl,cnrl,rl,
     .df1dnzl,df1dnrl,df2dnzl,df2dnrl
c      write(*,*)'output g123 cnz,cnr,r',cnz,cnr,r
c      write(*,*)'output g123 dnz,dnr,dr',dnz,dnr,dr
      cnzl=cnz+dnz
      cnrl=cnr+dnr
      rl=r+dr
c      write(*,*)'output g123 cnzl,cnrl,rl',cnzl,cnrl,rl
c      write(*,*)'output g123 cm,bphi,bmod,cnpar',cm,bphi,bmod,cnpar
c      write(*,*)'output g123 dbzdr,dbrdr,dbphdr,dbmdr',
c     .dbzdr,dbrdr,dbphdr,dbmdr
      df1drl=df1dr(cnzl,cnrl,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)
c      write(*,*)'output g123 df1drl',df1drl

      df2drl=df2dr(cm,rl)
c      write(*,*)'output g123 df2drl',df2drl

      cn=dsqrt(cnper**2+cnpar**2)

      g123ar(1)=f1(cnzl,cnrl,cm,rl,bz,br,bphi,bmod,cnpar)
      g123ar(2)=f2(cnzl,cnrl,cm,rl,cn)
     
c      write(*,*)'output g123 g123ar(1),g123ar(1)',g123ar(1),g123ar(2)

      df1dnzl = df1dnz(bz)
      df1dnrl = df1dnr(br)
      df2dnzl = df2dnz(cnz)
      df2dnrl = df2dnr(cnr)
c      write(*,*)'output g123  df1dnzl, df1dnrl,df2dnzl,df2dnrl',
c     .df1dnzl, df1dnrl,df2dnzl,df2dnrl

      call al12(dnz,dnr,df1dnzl,df1dnrl,df2dnzl,df2dnrl,al1,al2)
c      write(*,*)'output al1,al2',al1,al2

      g123ar(3)=dlagrdr(dr,al1,al2,df1drl,df2drl)
c       write(*,*)'output g123 g123ar(3)',g123ar(3)
      return
      end

      subroutine drvg123(cnz,cnr,cm,z,r,phi,dnr,dnz,dr,cnpar,cnper,
     .bz,br,bphi,bmod,dbzdr,dbrdr,dbphdr,dbmdr,
     .derg)
c     calculates derivatives derg from g123
      implicit none
c-----input 
      double precision cnz,cnr,cm,z,r,phi,dnz,dnr,dr,cnpar,cnper,
     .bz,br,bphi,bmod,dbzdr,dbrdr,dbphdr,dbmdr
c-----external
      double precision dlagrdr,df1dr,df2dr,f1,f2,
     .df1dnz,df1dnr,df2dnz,df2dnr
c-----output
      double precision derg(3,3)
c-----local
      double precision al1,al2,df1drl,df2drl, 
     .cnzl,cnrl,rl,dnzl,dnrl,drl,
     .df1dnzl,df1dnrl,df2dnzl,df2dnrl,
     .df1dnzp,df1dnrp,df2dnzp,df2dnrp,
     .df1dnzm,df1dnrm,df2dnzm,df2dnrm,
     .df1drp,df1drm,df2drp,df2drm,
     .step,cnzp,cnzm,cnrp,cnrm,rp,rm,g3p,g3m,
     .dnzp,dnzm,dnrp,dnrm,drp,drm,
     .bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl

      cnzl=cnz+dnz
      cnrl=cnr+dnr
      rl=r+dr


      df1drl=df1dr(cnzl,cnrl,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2drl=df2dr(cm,r)

      derg(1,1)=df1dnz(bz)                                  !dg1/ddnz
      derg(1,2)=df1dnr(br)                                  !dg1/ddnr
      derg(1,3)=df1dr(cnzl,cnrl,cm,rl,bphi,bmod,cnpar,      !dg1/ddr
     .dbzdr,dbrdr,dbphdr,dbmdr)

      derg(2,1)=df2dnz(cnzl)                                !dg2/ddnz
      derg(2,2)=df2dnr(cnrl)                                !dg2/ddnz
      derg(2,3)=df2dr(cm,rl)                                !dg2/ddr


 
      step=1.d-7   
c-----d(g3)/d(dnz)
      dnzp=dnz+step
      cnzp=cnz+dnzp

      df1drp=df1dr(cnzp,cnrl,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2dnzp = df2dnz(cnzp)

      call al12(dnzp,dnr,df1dnzl,df1dnrl,df2dnzp,df2dnrl,al1,al2)

      g3p=dlagrdr(dr,al1,al2,df1drp,df2drl)
      
      dnzm=dnz-step
      cnzm=cnz+dnzm

      df1drm=df1dr(cnzm,cnrl,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2dnzm = df2dnz(cnzm)

      call al12(dnzm,dnr,df1dnzl,df1dnrl,df2dnzm,df2dnrl,al1,al2)

      g3m=dlagrdr(dr,al1,al2,df1drm,df2drl)
      
      derg(3,1)=(g3p-g3m)/(2.d0*step)      !dg3ddnz

c-----d(g3)/d(dnr)
      dnrp=dnr+step
      cnrp=cnr+dnr

      df1drp=df1dr(cnzl,cnrp,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2dnrp = df2dnr(cnrp)

      call al12(dnz,dnrp,df1dnzl,df1dnrl,df2dnzl,df2dnrp,al1,al2)

      g3p=dlagrdr(dr,al1,al2,df1drp,df2drl)
      
      dnrm=dnr-step
      cnrm=cnr+dnrm

      df1drm=df1dr(cnzl,cnrm,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2dnrm = df2dnr(cnrm)

      call al12(dnz,dnrm,df1dnzl,df1dnrl,df2dnzl,df2dnrm,al1,al2)

      g3m=dlagrdr(dr,al1,al2,df1drm,df2drl)
      
      derg(3,2)=(g3p-g3m)/(2.d0*step)         !dg3ddnr

c-----dgr/d(dr)
      drp=dr+step
      rp=r+drp
  
      call bcomp(z,rp,phi,bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl)

      df1drp=df1dr(cnzl,cnrl,cm,rp,bphil,bmodl,cnpar,
     .dbzdrl,dbrdrl,dbphdrl,dbmdrl)

      df2drp=df2dr(cm,rp)

      df1dnzp = df1dnz(bzl)
      df1dnrp = df1dnr(brl)

      call al12(dnz,dnr,df1dnzp,df1dnrp,df2dnzl,df2dnrl,al1,al2)

      g3p=dlagrdr(drp,al1,al2,df1drp,df2drp)

      drm=dr-step
      rm=r+drm
  
      call bcomp(z,rm,phi,bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl)

      df1drm=df1dr(cnzl,cnrl,cm,rm,bphil,bmodl,cnpar,
     .dbzdrl,dbrdrl,dbphdrl,dbmdrl)

      df2drm=df2dr(cm,rm)

      df1dnzm = df1dnz(bzl)
      df1dnrm = df1dnr(brl)

      call al12(dnz,dnr,df1dnzm,df1dnrm,df2dnzp,df2dnrl,al1,al2)

      g3m=dlagrdr(drm,al1,al2,df1drp,df2drp)

      derg(3,3)=(g3p-g3m)/(2.d0*step)         !dg3ddr

      return
      end
    
      subroutine bcomp(z,r,phi,bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl)
c     calculates the magnetic field and its derivatives 
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      double precision b      
c      write(*,*)'bcomp z,r,phi',z,r,phi
      bmodl=b(z,r,phi)
c      write(*,*)'bcomp bz,br,bphi',bz,br,bphi       
      bzl=bz
      brl=br
      bphil=bphi
      dbzdzl=dbzdz
      dbzdrl=dbzdr
      dbrdzl=dbrdz
      dbrdrl=dbrdr
      dbphdzl=dbphdz
      dbphdrl=dbphdr
      dbmdzl=dbmdz
      dbmdrl=dbmdr

      return
      end


      subroutine set_output
c-----initialize output.i
      include 'param.i'
      include 'output.i'
        
      first=.true.
      was_not_ox_conversion=.true.

      return
      end
       


 
      subroutine map_disp_nper(id_plot_disp,m_r_nperp,m_i_nperp, 
     .max_r_nperp,min_r_nperp,max_i_nperp,min_i_nperp,
     .n_contour_plot_disp,
     .nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .vflow_ar,nll_in,re_nperp0,im_nperp0,iabsorp,
     .n_relt_harm1,n_relt_harm2,n_relt_intgr,
     .i_diskf,r,z,phi)
     
c     calculate the map for the complex dispersion function
c     on the 2D mesh (ReN_perp,ImN_perp)
c
      implicit none
c     input 
      integer id_plot_disp !sets the plotted dispesion function type
      integer m_r_nperp,m_i_nperp ! the number of points of the 2D-mesh
                                  ! for (ReN_perp,ImN_perp)
      integer n_contour_plot_disp !is the number of contours
   
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
c      iabsorp =4 hot plasma dispersion
c             =6 for Hermitian hot plasma tensor +
c                anti-Hhrmihian relativistic electron tensor
c     the data for relativistic anti Hermitian tensor
c     n_relt_harm1,n_relt_harm2
c     n_relat_intgr
c     i_diskf
c     z,r,phi the space coordinates
cid_plot_disp.eq.6
c
      integer nbulk,iabsorp,i_diskf,n_relt_harm1,n_relt_harm2,
     &n_relt_intgr   
      double precision mass_ar(*)
      double precision x_ar(*),y_ar(*),t_av_ar(*),
     .tpop_ar(*),vflow_ar(*)
      double precision nll_in,z,r,phi
      double precision re_nperp0,im_nperp0
     
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
      double complex aK(3,3), !relativistic anti Hermitian electron tensor 
     &compl_nper
      integer n_param
      parameter(n_param=7)
c      parameter(n_param=9)
      character*(15)name_param(n_param) !names of the input parameters
      real param(n_param)               !values of the input parameters
     
c     for pgplot
      
c-----parameter n_contour_plot_disp_a is the number of contours          
      include 'param.i' 
      
      real x_axis(m_r_nperp_a),y_axis(m_i_nperp_a),
     .f2d_r(m_r_nperp_a,m_i_nperp_a),f2d_i(m_r_nperp_a,m_i_nperp_a),
     .f2d_abs(m_r_nperp_a,m_i_nperp_a),
     .contour(n_contour_plot_disp_a)


c     write(*,*)'map_disp_nper numbers of points at ReN_perp,ImN_perp'
c     write(*,*)'m_r_nperp,m_i_nperp',m_r_nperp,m_i_nperp

      if (m_r_nperp_a.lt.m_r_nperp) then 
         write(*,*)'prep3d.f in map_dhot_nper m_r_nperp_a.lt.m_r_nperp'
         write(*,*)'m_r_nperp_a=',m_r_nperp_a
         write(*,*)'m_r_nperp=',m_r_nperp
         write(*,*)'change m_r_nperp_a in map_disp_nper'
         write(*,*)'and recompile the code'
         stop 'output.f subroutine map_disp_nper'
      endif

      if (m_i_nperp_a.lt.m_i_nperp) then 
         write(*,*)'prep3d.f in map_dhot_nper m_i_nperp_a.lt.m_i_nperp'
         write(*,*)'m_i_nperp_a=',m_i_nperp_a
         write(*,*)'m_i_nperp=',m_i_nperp
         write(*,*)'change m__i_nperp_a in map_disp_nper'
         write(*,*)'and recompile the code'
         stop 'output.f subroutine map_disp_nper'
      endif

c      write(*,*)'in map_disp_nperp:  nbulk,nll_in',nbulk,nll_in
c      write(*,*)'x_ar',(x_ar(i), i=1,nbulk)
c      write(*,*)'y_ar',(y_ar(i), i=1,nbulk)
c      write(*,*)'t_av_ar',(t_av_ar(i), i=1,nbulk)

      step_real=(max_r_nperp-min_r_nperp)/dble(m_r_nperp-1) !step ReN_perp
      step_imag=(max_i_nperp-min_i_nperp)/dble(m_i_nperp-1) !step ImN_perp

c      write(*,*)'max_r_nperp,min_r_nperp,step_real',
c     &max_r_nperp,min_r_nperp,step_real
c      write(*,*)'max_i_nperp,min_i_nperp,step_imag',
c     &max_i_nperp,min_i_nperp,step_imag
      
c      write(*,*)'output.f map_disp_nperp m_r_nperp,m_i_nperp',
c     &                                   m_r_nperp,m_i_nperp

c      write(*,*)'id_plot_disp',id_plot_disp
      if((id_plot_disp.ne.6).and.(id_plot_disp.ne.11).and.
     &   (id_plot_disp.ne.14)) then
        write(*,*)'**********************************************'
        write(*,*)'output.f in subroutine map_disp_nper'
        write(*,*)'id_plot_disp=',id_plot_disp
        write(*,*)'for this id_plot_disp'
        write(*,*)'map_disp_nper can not plot contour'
        write(*,*)'**********************************************'
        return
      endif

      do i=1,m_r_nperp
         nperp=min_r_nperp+step_real*(i-1) !Real N_perp
         x_axis(i)=nperp

         do j=1,m_i_nperp
            nprim=min_i_nperp+step_imag*(j-1)
            y_axis(j)=nprim

            if (id_plot_disp.eq.6) then 
c--------------full hot plasma dispersion function for Im_N_per.ne.0
               iherm=2
               d_compl=dhot_sum_c(nbulk,mass_ar,x_ar,y_ar,t_av_ar,
     .         tpop_ar,vflow_ar,nll_in,nperp,nprim,iherm,K_sum)
            endif !id_plot_disp.eq.6

            if (id_plot_disp.eq.11) then 
c--------------relativistic dispersion function
c              from Eric Nelson_Melby dielectric tensor
               compl_nper=dcmplx(nperp,nprim)
               call Disp_Nelson_Melby(t_av_ar(1),nll_in,x_ar(1),y_ar(1),
     +              compl_nper,K_sum, d_compl) ! K_sum is in z-y plane
                                               ! in Stix coordinate system              
            endif !id_plot_disp.eq.11)

            if (id_plot_disp.eq.14) then 
c--------------relativistic dispersion function
c             is just the Ram tensor
               compl_nper=dcmplx(nperp,nprim)   
c               write(*,*)'output i,j,t_av_ar(1),nll_in',
c     +         i,j,t_av_ar(1),nll_in
c               write(*,*)'x_ar(1),y_ar(1),compl_nper',
c     +         x_ar(1),y_ar(1),compl_nper
               call Disp_Ram(t_av_ar(1),nll_in,x_ar(1),y_ar(1),
     +              compl_nper,K_sum,d_compl)  !K is in z-y plane
                               !in Stix coordinates
c               write(*,*)'d_compl',d_compl

            endif !id_plot_disp.eq.14)

c            write(16,2)nperp,nprim,dreal(d_compl),dimag(d_compl),
c     .      dsqrt((dreal(d_compl))**2+(dimag(d_compl))**2)

            f2d_r(i,j)=real(dreal(d_compl))
            f2d_i(i,j)=real(dimag(d_compl))
            f2d_abs(i,j)=sqrt(f2d_r(i,j)**2+f2d_i(i,j)**2)
         enddo
      enddo

cfor test
c      write(*,*)'in map_disp_nper'
c      do i=1,m_r_nperp
c        do j=1,m_i_nperp
c           write(*,*)'i,j,x_axis(i),y_axis(j)',i,j,x_axis(i),y_axis(j)
c           write(*,*)'f2d_i(i,j)',f2d_i(i,j)
c        enddo
c      enddo
cend test   


c      call plotinit
      name_param(1)='Y'
      name_param(2)='X'
      name_param(3)='N_parallel'
      name_param(4)='T'
c      name_param(5)='ReNperp_0'
c      name_param(6)='ImNperp_0'
      name_param(5)='r' 
      name_param(6)='z'
      name_param(7)='id_plot_disp'
      
      param(1)=y_ar(1)
      param(2)=x_ar(1)
      param(3)=nll_in
      param(4)=t_av_ar(1)
c      param(5)=re_nperp0
c      param(6)=im_nperp0
      param(5)=r  
      param(6)=z
      param(7)=id_plot_disp

cSm030519
c      call contour2d(m_r_nperp,m_i_nperp,f2d_r,x_axis,y_axis,
c     .'real_D', 'ReN_perp','ImN_perp',
c     .n_contour,contour,
c     .name_param,param,n_param)


c      call contour2d(m_r_nperp,m_i_nperp,f2d_i,x_axis,y_axis,
c     .'ImageD', 'ReN_perp','ImN_perp',
c     .n_contour,contour,
c     .name_param,param,n_param)

c      call contour2d(m_r_nperp,m_i_nperp,f2d_abs,x_axis,y_axis,
c     .'absD', 'ReN_perp','ImN_perp',
c     .n_contour,contour,
c     .name_param,param,n_param)

c      write(*,*)'output.f in map_disp_nper n_contour_plot_disp',
c     &n_contour_plot_disp

      call contour2d_S(m_r_nperp_a,m_i_nperp_a,n_contour_plot_disp_a,
     &m_r_nperp,m_i_nperp,f2d_r,x_axis,y_axis,
     .'real_D', 'ReN_perp','ImN_perp',
     .n_contour_plot_disp,contour,
     .name_param,param,n_param)

      call contour2d_S(m_r_nperp_a,m_i_nperp_a,n_contour_plot_disp_a,
     &m_r_nperp,m_i_nperp,f2d_i,x_axis,y_axis,
     .'ImageD', 'ReN_perp','ImN_perp',
     .n_contour_plot_disp,contour,
     .name_param,param,n_param)

      call contour2d_S(m_r_nperp_a,m_i_nperp_a,n_contour_plot_disp_a,
     &m_r_nperp,m_i_nperp,f2d_abs,x_axis,y_axis,

     .'absD', 'ReN_perp','ImN_perp',
     .n_contour_plot_disp,contour,
     .name_param,param,n_param)


      return
      end

      subroutine output_map_disp_nperp(u)
c-------------------------------------------------------------
c     plots dispersion function contours D(ReN_perp,Im_N_perp)
c     inside outpt subroutine
c-------------------------------------------------------------
      implicit none

      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'eps.i'
      include 'output.i'

      real*8 u(*)
      integer 
     &m_r_nperp, !number of map points in ReN_perp dimension
     &m_i_nperp, !number of map points in InN_perp dimension
     &iter_max   !max number of the iterations in Newton method

      integer i

      real*8 
     &max_r_nperp,min_r_nperp, !max and min ReN_perp
     &max_i_nperp,min_i_nperp, !max and min ImN_perp
     &bf(3),                   !magnetic field: bz,br,bphi
     &cnper,                   !ReN_perpendicular
     &cnpar,                   !ReN_parallel   
     &cnprim,                  !ImN_perpendicular
     &z,r,phi,cnz,cnr,cm,cnphi,
     &dham,
     &cnperp,cnparp,cnt,ds,dc,te,cnper_new,cnprim_old

      real*8 x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     &tpop_ar(nbulka),vflow_ar(nbulka)

      complex*16 dd5(5,nbulka),ddnp_h,ddnll_h,ddnp,
     &K(3,3),dK(3,3,7),dd(5),d,ddn,aK(3,3)
c-----------------------------------------------------------------
c     external
c-------------------------------------------------------------------
      real*8 b,gamma1,cn,tempe,x,y,tpoprho,vflowrho
      complex*16 dhot_sum

c      write(*,*)'output.f output_map_disp_nperp'
c      write(*,*)'ratio_min_r_nperp,ratio_max_r_nperp,cnper',
c     &ratio_min_r_nperp,ratio_max_r_nperp,cnper
c      write(*,*)'=number_map_points_real_nperp',
c     &number_map_points_real_nperp
c      write(*,*)'=number_map_points_image_nperp',
c     &number_map_points_image_nperp
 
      m_r_nperp=number_map_points_real_nperp 
      m_i_nperp=number_map_points_image_nperp 

      z=u(1)
      r=u(2)
      phi=u(3) 
      cnz=u(4)
      cnr=u(5)
      cm=u(6) 
      cnphi=u(6)/u(2) 
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

c      max_r_nperp=cnper*1.5d0
c      min_r_nperp=cnper*0.5d0  
      max_r_nperp=cnper*ratio_max_r_nperp
      min_r_nperp=cnper*ratio_min_r_nperp


      do i=1,nbulk
         x_ar(i)=x(z,r,phi,i)
	 y_ar(i)=y(z,r,phi,i)
         if(i.eq.1) y_ar(1)=-y_ar(1)
	 te=tempe(z,r,phi,i) ! kev
	 t_av_ar(i)=te*1000.d0      ! ev 
         tpop_ar(i)=tpoprho(rho,i)
         vflow_ar(i)=vflowrho(rho,i)
      enddo

c      write(*,*)'output_map_disp_nperp(u) i_im_nperp',i_im_nperp
c---------------------------------------------------------------
c     calculate Im N_perp using hot plasma dispersion function
c---------------------------------------------------------------
      if(i_im_nperp.eq.1) then
c-------calculate ImN_perp using the formula
c       ImN_perp=abs(ImD_full/dD_hermitian/dReN_perp))
        cnparp=cnpar
        cnperp=cnper
 
        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &    vflow_ar,cnparp,cnperp,2,reps)

        dham=dreal(d)

	call Ddhot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
     &  ,cnparp,cnperp,reps,dd5,ddnp_h,ddnll_h,ddnp)
c       write(*,*)' output_map_disp_nperp(u) d,ddnp',d,ddnp
        cnprim = dabs(DIMAG(D) / DREAL(ddnp))
        write(*,*)'output cnprim=(ImD/dD/dn_perp)= ',cnprim

        cnper_new=cnperp !it will be used in the electric field calculations
        write(*,*)'cnper,cnper_new',cnper,cnper_new
      endif  ! i_im_nperp.eq.1  

      if(i_im_nperp.eq.2) then 
c---------find (Im_N_perp,ReN_perp) the root of the complex dispersion relation
c         using the Newton method with numerical derivatives (the chord method)
          iter_max=100 !max number of the iterations
          iter_max=3 !max number of the iterations
c---------initial values of Im_N_perp=cnprim Re_N_perp=cnperp
          cnprim=0.d0   
          write(*,*)'prep3d before call solv_nperp_hot cnprim=',cnprim
          call solv_nperp_hot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .    vflow_ar,cnparp,cnperp,iter_max,cnper_new,cnprim)
          write(*,*)'Newton cnperp,cnper_new,cnprim',
     .    cnperp,cnper_new,cnprim
          cnprim_old=cnprim
          cnprim=dabs(cnprim)
      endif !i_im_nperp.eq.2

      if (dabs(cnprim).gt.1.d0) then        
c        max_i_nperp= 2.5d0*dabs(cnprim)
c        min_i_nperp=-2.5d0*dabs(cnprim)
        max_i_nperp=dabs(cnprim)*ratio_min_i_nperp
        min_i_nperp=dabs(cnprim)*ratio_max_i_nperp
      else
        max_i_nperp= 1.d0
        min_i_nperp=0.d0
      endif

c----------------------------------------------------------------
c     plotting dispersion function contours D(ImN_perp,ReN_perp)
c     to plot.ps file using PGplot
c----------------------------------------------------------------

      call map_disp_nper(id_plot_disp,m_r_nperp,m_i_nperp, 
     .max_r_nperp,min_r_nperp,max_i_nperp,min_i_nperp,
     .n_contour_plot_disp,
     .nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .vflow_ar,cnpar,cnperp,cnprim,iabsorp,
     .n_relt_harm1,n_relt_harm2,n_relt_intgr,
     .i_diskf,r,z,phi)
    
      return
      end

      subroutine plot_disp_cold_output(u,t,n_nperp)
c---------------------------------------------------------------------------
c     Creates plots of the cold disperion function and writes
c     these plots to plot.ps file
c     This subroutine is called from output
c----------------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'

c-----input
      
      real*8 u(*),! ray coordinates
     &t           ! poloidal distance or time along the ray
      integer n_nperp !number of map points at n_perpendicular interval 

c-----parameters
      integer n_param
      parameter(n_param=2)
     
c-----externals
      real*8 b

c-----local
      real*8 z,r,phi,cnz,cnr,cm,cnpar,cnper2m,cnper2p,
     &cnper_max,cnper_min 

      character*(20)name_param(n_param) !names  of the input parameters
      real param(n_param)               !values of the input parameters
      integer i_n_per_min_negative,i_n_per_max_negative 

c      write(*,*)'plot_disp_cold_output n_param,n_nperp',n_param,n_nperp

      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6) 
      bmod=b(z,r,phi)
      cnpar=(cnz*bz+cnr*br+(cm/r)*bphi)/bmod  
               
      call npernpar(z,r,phi,cnpar,cnper2p,cnper2m)

      cnper_min=dmin1(cnper2p,cnper2m)
      cnper_max=dmax1(cnper2p,cnper2m)

      if(cnper_min.gt.0.d0) then 
        cnper_min=dsqrt(cnper_min) !minimal value of nperp
        i_n_per_min_negative=0 
      else
        cnper_min=0.d0
        i_n_per_min_negative=1
      endif

      if(cnper_max.gt.0.d0) then 
        cnper_max=dsqrt(cnper_max) !maximal value of nperp
        i_n_per_max_negative=0
      else
        cnper_max=1.d0
        i_n_per_max_negative=1
      endif

      name_param(1)='r'
      name_param(2)='poloidal distance'
      param(1)=r
      param(2)=t

      call map_d_cold(z,r,phi,cnpar,n_nperp,
     &cnper_min*0.9d0,cnper_max*1.1d0,name_param,param,n_param)

      if (  i_n_per_min_negative.eq.0) then
         call map_d_cold(z,r,phi,cnpar,n_nperp,
     &   cnper_min*0.9d0,cnper_min*1.1d0,name_param,param,n_param)
      endif

      if (  i_n_per_max_negative.eq.0) then
         call map_d_cold(z,r,phi,cnpar,n_nperp,
     &   cnper_max*0.9d0,cnper_max*1.1d0,name_param,param,n_param)
      endif

      return
      end

      subroutine plot_disp_cold_grill(z,r,phi,cnpar,t,n_nperp)
c---------------------------------------------------------------------------
c     Creates plots of the cold disperion function and writes
c     these plots to plot.ps file
c     This subroutine is called from grill_lh.f
c----------------------------------------------------------------------------
      implicit none
     
c-----input
      
      real*8 z,r,phi, ! ray coordinates
     &cnpar,          ! n_parallel
     &t               ! poloidal distance or time along the ray
     
      integer n_nperp !number of map points at n_perpendicular interval 

c-----parameters
      integer n_param
      parameter(n_param=2)
c-----local
      real*8 cnper2m,cnper2p,cnper_max,cnper_min,rho
      integer i,i_n_per_min_negative,i_n_per_max_negative 

      character*(20)name_param(n_param) !names  of the input parameters
      real param(n_param)               !values of the input parameters



      write(*,*)'plot_disp_cold_grill n_param,n_nperp',n_param,n_nperp
                     
      call npernpar(z,r,phi,cnpar,cnper2p,cnper2m)

      cnper_min=dmin1(cnper2p,cnper2m)
      cnper_max=dmax1(cnper2p,cnper2m)

      if(cnper_min.gt.0.d0) then 
        cnper_min=dsqrt(cnper_min) !minimal value of nperp
        i_n_per_min_negative=0 
      else
        cnper_min=0.d0
        i_n_per_min_negative=1
      endif

      if(cnper_max.gt.0.d0) then 
        cnper_max=dsqrt(cnper_max) !maximal value of nperp
        i_n_per_max_negative=0
      else
        cnper_max=1.d0
        i_n_per_max_negative=1 
      endif

      name_param(1)='r'
      name_param(2)='poloidal distance'
      param(1)=r
      param(2)=t


      write(*,*)'plot_disp_cold_grill'
      write(*,*)'cnpar=',cnpar,'n_nperp=',n_nperp
      write(*,*)'cnper_min,cnper_max',cnper_min,cnper_max
      write(*,*)'n_param=',n_param
      write(*,*)'name_param=',(name_param(i),i=1,n_param)
      write(*,*)'param',(param(i),i=1,n_param)
      write(*,*)'before map_d_cold' 

      call map_d_cold(z,r,phi,cnpar,n_nperp,
     &cnper_min*0.9d0,cnper_max*1.1d0,name_param,param,n_param)

      if (  i_n_per_min_negative.eq.0) then
         call map_d_cold(z,r,phi,cnpar,n_nperp,
     &   cnper_min*0.9d0,cnper_min*1.1d0,name_param,param,n_param)
      endif
     
      if (  i_n_per_max_negative.eq.0) then
         call map_d_cold(z,r,phi,cnpar,n_nperp,
     &   cnper_max*0.9d0,cnper_max*1.1d0,name_param,param,n_param)
      endif

      return
      end


      subroutine refractive_index_relative_error(z,r,phi,cnz,cnr,cm,
     &iraystop)
c-------------------------------------------------------------------
c     Calculate the relative error of delta_N/N=delta_n_devide_n
c     of the dispersion relation D(cnz,cnr,cm)=0 at the given point
c     (z,z,phi,cnz,cnr,cm)
c     delta_N/N=(D/N)/|gradD|
c     Here
c     D=D(z,z,phi,cnz,cnr,cm)
c     N=|N|=sqrt(Nz**2+Nr**2+(cm/r)**2)
c          N_phi=cm/r
c     |gradD|=sqrt((dD/dN_z)**2+(dD/dN_r)**2+(dD/dN_phi)**2)=
c            =sqrt((dD/dN_z)**2+(dD/dN_r)**2+(dD/dcm)**2*r**2)
c
c     If delta_N/N=(D/N)/|gradD| > toll_hamilt it will put iraystop=1
c     else         iraystop=0
c     Variable toll_hamilt is set in one.i
c------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 z,r,phi, !space coordinates
     *cnz,cnr,cm      !refractive index coordinates
c-----output
      integer iraystop 
        
c-----externals dddrz1
      real*8 b,gamma1,hamilt1

c-----locals
      real*8 u(6),deru(6),
     &grad_d,  ! |gradD in N space|
     &cn,d,     
     &delta_n_devide_n ! delta_N/N=(D/N)/|gradD|

      u(1) = z
      u(2) = r                                                    
      u(3) = phi                                                  
      u(4) = cnz                                                  
      u(5) = cnr                                                  
      u(6) = cm
c      write(*,*)'in  refractive_index_relative_error'
      bmod=b(z,r,phi)
      gam=gamma1(z,r,phi,cnz,cnr,cm)

c      write(*,*)'in  refractive_index_relative_error gam',gam

      call  dddrz1(0.d0,u,deru)

c      write(*,*)'in refractive_index_relative_error after dddrz1'

      d=hamilt1(z,r,phi,cnz,cnr,cm)
 
c      write(*,*)'in  refractive_index_relative_error d=',d


      grad_d=dsqrt(deru(1)**2+deru(2)**2+deru(3)**2*r**2)
      cn=dsqrt(cnz**2+cnr**2+(cm/r)**2)
       
      delta_n_devide_n=d/(cn*grad_d)

c      write(*,*)'delta_n_devide_n=',delta_n_devide_n
c      write(*,*)'toll_hamilt=',toll_hamilt

      if( delta_n_devide_n.gt.toll_hamilt) then 
         iraystop=1
      else
         iraystop=0
      endif

      if(iraystop.eq.1) then
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*) 'outpt: in  refractive_index_relative_error'
         write(*,*) ' D/(N|gradD|) > toll_hamilt'
         write(*,*) 'D=',d,'N=',cn,'grad_d=',grad_d
         write(*,*) 'D/(N|gradD|)=delta_n_devide_n',delta_n_devide_n
         write(*,*) 'toll_hamilt=',toll_hamilt
         write(*,*) '*************************************************'
         write(*,*)
      endif
   
      return
      end

      subroutine map_dcold_nz_nr(m_nz,m_nr, 
     &max_nz,min_nz,max_nr,min_nr,n_contour,
     &z,r,phi,cnz,cnr,cm)
     
c     calculate the map for the cold plasma dispersion function
c     on the 2D mesh (N_z,N_r)
c
      implicit none

      include 'param.i'
      include 'one.i'

c     input 
      integer m_nz,m_nr, !the number of points of the 2D-mesh
                         ! for (N_z,N_r)
     .n_contour          ! the number of contours

      real*8
     &max_nz,min_nz, !max and min N_z
     &max_nr,min_nr, !max and min N_r
     &z,r,phi,      !the space coordinates
     &cnz,cnr,cm    !refractive index components
     
c-----locals                
      integer i,j
    
      real*8 step_nz,step_nr,n_r,n_z
    
      integer n_param
      parameter(n_param=4)
      character*(10)name_param(n_param) !names of the input parameters
      real*4 param(n_param)               !values of the input parameters
c-----for contours of the coefficient a of cold dispersion relation
      real*8 ds,dc,ds2,dc2,cnt,ad,bd,cd,
     &cn2,cn4,d4,d2,d0,hamilt_id1,hamilt_id2
c-----for pgplot

c-----n_contour_plot_disp is the number of contours           

      real*4, dimension (1:m_nr) ::  x_axis
      real*4, dimension (1:m_nz) ::  y_axis
      real*4, dimension (1:m_nz, 1:m_nz) :: f2d_a_r,f2d_b_r,f2d_c_r
      real*4, dimension (1:n_contour) :: contour

c-----externals
      real*8 b,gamma1,hamilt1,cn

      bmod=b(z,r,phi)

      step_nr = (max_nr-min_nr)/dble(m_nr-1)                !step Nr
      step_nz = (max_nz-min_nz)/dble(m_nz-1)                !step Nz
 
      write(*,*)'************map_dcol_nz_nr********************'
      write(*,*)'max_nr,min_nr, step_nr', max_nr,min_nr, step_nr
      write(*,*)'max_nz,min_nz, step_nz', max_nz,min_nz, step_nz

      do i=1,m_nr
         n_r=min_nr+step_nr*(i-1) !N_r
         x_axis(i)=n_r

         do j=1,m_nz
            n_z=min_nz+step_nz*(j-1)
            y_axis(j)=n_z
 
            gam=gamma1(z,r,phi,n_z,n_r,cm)
            f2d_a_r(i,j)=hamilt1(z,r,phi,n_z,n_r,cm)
c---------------------------------------------------------
c           write(*,*)'i,j,x_axis(i),y_axis(j),f2d_a_r(i,j)',
c     &                i,j,x_axis(i),y_axis(j),f2d_a_r(i,j)
c            ds=dsin(gam)
c            dc=dcos(gam)
c            ds2=ds*ds
c            dc2=dc*dc 

c            cnt=cn(r,n_z,n_r,cm)
c            cn2=cnt*cnt
c       	    cn4=cn2*cn2
c            call abc(z,r,phi,ds2,dc2,ad,bd,cd)
c            d4=ad
c       	    d2=bd
c       	    d0=cd
c            hamilt_id1=d4*cn4+d2*cn2+d0
c            write(*,*)'hamilt_id1',hamilt_id1
c            hamilt_id2=cn2+(d2-ioxm*dsqrt(d2*d2-4.d0*d4*d0))/
c     *                (2.d0*d4)
c            write(*,*)'hamilt_id2',hamilt_id2                                     
         enddo
      enddo

c      write(*,*)'x_axis',x_axis
c      write(*,*)'y_axis',y_axis
      write(*,*)'*********the map of the hamitonian******************'     
c      do i=1,m_nr
c         do j=1,m_nz          
c            write(*,*)'i,j,x_axis(i),y_axis(j),f2d_a_r(i,j)',
c     &                 i,j,x_axis(i),y_axis(j),f2d_a_r(i,j)
c         enddo
c      enddo

      name_param(1)='z'
      name_param(2)='r'
      name_param(3)='N_z'
      name_param(4)='N_r'
c      name_param(5)='ReNperp_0'
c      name_param(6)='ImNperp_0'

      param(1)=real(z)
      param(2)=real(r)
      param(3)=real(cnz)
      param(4)=real(cnr)

c      param(5)=re_nperp0
c      param(6)=im_nperp0

c      call contour2d_S(m_nr,m_nz,n_contour,
c     .m_nr,m_nz,f2d_r,x_axis,y_axis,
c     .'D', 'N_z','N_r',
c     .n_contour,contour,
c     .name_param,param,n_param)

      call contour2d_1(m_nr,m_nz,f2d_a_r,x_axis,y_axis,
     .'D', 'N_r','N_z',
     .n_contour,contour,
     .name_param,param,n_param)

c--------------------------------------------------------------------
c     contours coefficient in cold plasma dispersion relation
c     D=AN**4+B**N**2+C     
      do i=1,m_nr
         n_r=x_axis(i)
         do j=1,m_nz
            n_z=y_axis(j)
            gam=gamma1(z,r,phi,n_z,n_r,cm)
            ds=dsin(gam)
            dc=dcos(gam)
            ds2=ds*ds
            dc2=dc*dc
     
            cnt=cn(r,cnz,cnr,cm)
            call abc(z,r,phi,ds2,dc2,ad,bd,cd)  
          
c            write(*,*)'i,j,n_r,n_z,ds2,dc2,ad,bd,cd',
c     &                 i,j,n_r,n_z,ds2,dc2,ad,bd,cd
c            write(*,*)'i,j',i,j
c           write(*,*)'cnt**2,(bd-ioxm*dsqrt(bd**2-4.d0*ad*cd)/(2*ad))',
c     &                cnt**2,(bd-ioxm*dsqrt(bd**2-4.d0*ad*cd)/(2*ad))
c           write(*,*)'cnt**2+(bd-ioxm*dsqrt(bd**2-4.d0*ad*cd)/(2*ad))',
c     &                cnt**2+(bd-ioxm*dsqrt(bd**2-4.d0*ad*cd)/(2*ad))

            f2d_a_r(i,j)=ad
            f2d_b_r(i,j)=bd      
            f2d_c_r(i,j)=cd    
         enddo
      enddo

      write(*,*)'****** map of the cold dispersion coefficients***'
c      do i=1,m_nr
c         do j=1,m_nz
c            write(*,*)'i,j,x_axis(i),y_axis(j),f2d_a_r(i,j)',
c     &                 i,j,x_axis(i),y_axis(j),f2d_a_r(i,j)
c            write(*,*)'f2d_b_r(i,j), f2d_c_r(i,j)',
c     &                 f2d_b_r(i,j), f2d_c_r(i,j)
c         enddo
c      enddo

      call contour2d_1(m_nr,m_nz,f2d_a_r,x_axis,y_axis,
     .'A', 'N_r','N_z',
     .n_contour,contour,
     .name_param,param,n_param)

      call contour2d_1(m_nr,m_nz,f2d_b_r,x_axis,y_axis,
     .'B', 'N_r','N_z',
     .n_contour,contour,
     .name_param,param,n_param)

      call contour2d_1(m_nr,m_nz,f2d_c_r,x_axis,y_axis,
     .'C', 'N_r','N_z',
     .n_contour,contour,
     .name_param,param,n_param)

      return
      end



      subroutine jump_through_resonace(us,u_in,h,u_out,prmt7)
c-----jump in RK through gyro resonace point 
c     input:
c     y_in    (omega_c/omega) algebraic
c     del_y the distance in Y from the resonace to jump in one.i
c     u_in(6)  ray coordinates befor resonance 
c     h RK time step
c     us distance along the ray at the in point
c-----output:
c     u_out(6) ray coordinates after resonance
c     us
c     prmt7=us
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 u_in(6),h,us

c-----output
      real*8
     &u_out(6),prmt7
    
c-----local
      integer iy,i
      real*8 deru(6),y_in,y_out,h_out,
     &dif_y_in,dif_y_out,
     &dels2,dels1

c-----externals
      real*8 y,b

c      write(*,*)'jump_through_resonace'

      do i=1,6
         u_out(i)=u_in(i)
c         write(*,*)'i,u_out(i),u_in(i)',i,u_out(i),u_in(i)
      enddo

      bmod=b(u_in(1),u_in(2),u_in(3))
c      write(*,*)'bmod',bmod
      y_in=y(u_in(1),u_in(2),u_in(3),1)

      iy=nint(1.d0/y_in)
      h_out=h

      dif_y_in=1.d0/y_in-dfloat(iy)
       
c      write(*,*)'jump_through_resonace h,y_in,dif_y_in,del_y',
c     &h,y_in,dif_y_in,del_y

      if(dabs(dif_y_in).le.del_y) then
c------ the jump
        write(*,*)'dabs(dif_y_in).le.del_y'        
        call rside1(0.d0,u_in,deru)
 10     continue
        do i=1,6
           u_out(i)=u_in(i)+h_out*deru(i)
        enddo

        bmod=b(u_out(1),u_out(2),u_out(3))
        y_out=y(u_out(1),u_out(2),u_out(3),1)
        dif_y_out=1.d0/y_out-dfloat(iy)
c        write(*,*)'h_out, y_out,dif_y_out,del_y',
c     &             h_out, y_out,dif_y_out,del_y
        if ((dif_y_in*dif_y_out.gt.0.d0).or.
     &      (dabs(dif_y_out).le.del_y)) then
c----------in and out points are at the same side from the resonance point
           h_out=h_out*2.d0
           go to 10  
        else
c----------in and out points are at different sides from the resonance point
           if (i_output.eq.1) then
c              poloidal distance
               dels2=(u_out(1)-u_in(1))**2+(u_out(2)-u_in(2))**2
           endif

           if (i_output.eq.2) then
c              total distance
               dels2=(u_out(1)-u_in(1))**2+(u_out(2)-u_in(2))**2+
     &         (u_out(2)*u_out(3)-u_in(2)*u_in(3))**2
           endif

           dels1=dsqrt(dels2)

           prmt7=us+h_out
           us=prmt7
        endif
      endif
  
      return
      end 
c=====================================================================
c     functions to plot D(ImNperp,ReN_perp)
c     at given RZ N_parallel points
c=====================================================================
      subroutine plot_disp_D_at_RZ_n_parallel_points
c-------------------------------------------------------------
c     plot hot D(ReN_perp,ImN_perp) 
c     dispersion relation contous  at given RZ_N_parallel points 
c------------------------------------------------------------------------------
c id_plot_disp  determines the dispersion function D type
c          used for contours plots
c          Dispersion function contours are ploted
c          only at   id_plot_disp=6    - full hot plasma tensor
c                    id_plot_disp=11   - relativistic dispersion function
c                                        from Eric Nelson_Melby dielectric
c                                        tensor
c                    id_plot_disp=14   - relativistic dispersion function
c                                       just the Ram tensor
c
c at id_plot_disp=8 and point_plot_disp = 'rz_nparallel' 
c                   instead of D contours the code will plot Real part
c                   of the Ono dispersion function D_Ono(N_perpendicular)
c                   versur Real(N_perpendicular) at given N_parallel
c                   In this case at each RZ_N_parallel points with number
c                   i=1:n_plot_disp the code
c                   will plot 1D curve Re(D(N_perp)) at the N_perp interval
c                   [min_i_nperp_plot_disp(i)=<N_perp=<max_i_nperp_plot_disp(i)]
c                   using uniform n_perp mesh 1:number_map_points_real_nperp.
c                   About number_map_points_real_nperp see bellow.
c                   
c-------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'  
      include 'output.i'
c-----input    
      integer ::
     &m_r_nperp, !number of map points in ReN_perp dimension
     &m_i_nperp  !number of map points in InN_perp dimension 
      real*8 ::
     &max_r_nperp,min_r_nperp, !max and min ReN_perp
     &max_i_nperp,min_i_nperp !max and min ImN_perp

c-----locals
      integer :: i 
      real*8 ::
     &r_l,z_l,     !rz point coordinates
     &phi_l,
     &cnpar_l    !N parallel
      write(*,*)'in subroutine plot_dispd_D_at_RZ_n_parallel_points'

      do i=1,n_plot_disp
        r_l=     r_plot_disp(i)
        z_l=     z_plot_disp(i)  
        phi_l=   0.d0
        cnpar_l= n_parallel_plot_disp(i)

        write(*,*)'i,r_l,z_l,phi_l,cnpar_l=',
     &             i,r_l,z_l,phi_l,cnpar_l

        max_r_nperp=max_r_nperp_plot_disp(i)
        min_r_nperp=min_r_nperp_plot_disp(i)
        max_i_nperp=max_i_nperp_plot_disp(i)
        min_i_nperp=min_i_nperp_plot_disp(i)

        m_r_nperp=number_map_points_real_nperp 
        m_i_nperp=number_map_points_image_nperp 

        if (id_plot_disp.eq.8) then 
          call map_d_Ono(z_l,r_l,phi_l,cnpar_l,
     &                   max_r_nperp,min_r_nperp)
        else   
          call output_map_disp_nperp_RZ_N_parallel(r_l,z_l,cnpar_l,
     &    phi_l,
     &    m_r_nperp,m_i_nperp,
     &    max_r_nperp,min_r_nperp,max_i_nperp,min_i_nperp)
        endif

      enddo

      return
      end

      subroutine output_map_disp_nperp_RZ_N_parallel(r,z,cnpar,phi,
     & m_r_nperp,m_i_nperp,
     &max_r_nperp,min_r_nperp,max_i_nperp,min_i_nperp)
c-------------------------------------------------------------
c     plots dispersion function contours D(ReN_perp,Im_N_perp)
c     inside outpt subroutine
c--------------max ReN_perp boundary for plot-----------------------------------------------
      implicit none

      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'eps.i'
      include 'output.i'

c-----input
      real*8 ::
     &r,z,     !rz point coordinates
     &phi,
     &cnpar    !N parallel

      integer ::
     &m_r_nperp, !number of map points in ReN_perp dimension
     &m_i_nperp  !number of map points in InN_perp dimension 
      real*8 ::
     &max_r_nperp,min_r_nperp, !max and min ReN_perp
     &max_i_nperp,min_i_nperp !max and min ImN_perp

c-----locals
      real*8 :: te , cnperp, cnprim
      integer i

      real*8, dimension (1: nbulk) :: x_ar,y_ar,
     &t_av_ar,tpop_ar,vflow_ar

c-----------------------------------------------------------------
c     external
c-------------------------------------------------------------------
      real*8 b,tempe,x,y,tpoprho,vflowrho
      
c      m_r_nperp=number_map_points_real_nperp 
c      m_i_nperp=number_map_points_image_nperp 
      write(*,*)'in subroutine output_map_disp_nperp_RZ_N_parallel'
      write(*,*)'z,r,phi',z,r,phi
      bmod=b(z,r,phi) 
      write(*,*)'rho',rho
c      max_r_nperp=0.d0
c      min_r_nperp=10.d0
c      max_i_nperp=0.d0
c      min_i_nperp=10.d0
      

      do i=1,nbulk
         x_ar(i)=x(z,r,phi,i)
	 y_ar(i)=y(z,r,phi,i)
         if(i.eq.1) y_ar(1)=-y_ar(1)
	 te=tempe(z,r,phi,i) ! kev
	 t_av_ar(i)=te*1000.d0      ! ev 
         tpop_ar(i)=tpoprho(rho,i)
         vflow_ar(i)=vflowrho(rho,i)
        write(*,*)'i,x_ar(i),y_ar(i),t_av_ar(i),tpop_ar(i),vflow_ar(i)'
     &            ,i,x_ar(i),y_ar(i),t_av_ar(i),tpop_ar(i),vflow_ar(i)
      enddo

c----------------------------------------------------------------
c     plotting dispersion function contours D(ImN_perp,ReN_perp)
c     to plot.ps file using PGplot
c----------------------------------------------------------------
      cnperp=0.d0  !this argument is not used 
      cnprim=0.d0  !this argument is not used

      write(*,*)'id_plot_disp,m_r_nperp,m_i_nperp',
     &           id_plot_disp,m_r_nperp,m_i_nperp
      write(*,*)'max_r_nperp,min_r_nperp,max_i_nperp,min_i_nperp',
     &           max_r_nperp,min_r_nperp,max_i_nperp,min_i_nperp
      write(*,*)'n_contour_plot_disp',n_contour_plot_disp
      write(*,*)'iabsorp,n_relt_harm1,n_relt_harm2,n_relt_intgr',
     &           iabsorp,n_relt_harm1,n_relt_harm2,n_relt_intgr
      write(*,*)'i_diskf,r,z,phi',i_diskf,r,z,phi

      call map_disp_nper(id_plot_disp,m_r_nperp,m_i_nperp, 
     .max_r_nperp,min_r_nperp,max_i_nperp,min_i_nperp,
     .n_contour_plot_disp,
     .nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .vflow_ar,cnpar,cnperp,cnprim,iabsorp,
     .n_relt_harm1,n_relt_harm2,n_relt_intgr,
     .i_diskf,r,z,phi)
    
      return
      end

      subroutine map_d_Ono(z,r,phi,npar,nperp_min,nperp_max)
c-----creates----------------------------------------------------------
c     1D array D(nperp) for plotting of the Ono plasma dispersion function
c     D(nperp) versus perpendicular refractive index N_perpendicular=nperp
c     at the given point (r,z,phi) and given N_parallel=npar
c
c     Plots D_Ono(ReN_perp) to plot.ps file
c---------------------------------------------------------------------------- 
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'eps.i'
      include 'output.i'

c-----input:
      integer n_nperp                         !number of points in nperp mesh
      real*8 nperp_min,nperp_max,   !minimal and maximal nperp values
     &z,r,phi, ! space point coordinates
     &npar     ! parallel refractive index


c-----locals
      integer n_param                    !number of input parameters
      parameter (n_param=7)
      character*(15) name_param(n_param)  !names  of the input parameters
      real param(n_param)                 !values of the input parameters   

      real*8, dimension(1:number_map_points_real_nperp) :: 
     &d_ar,     !array of dispersion function values
     &nperp_ar  !mesh of nperp

      real*8, dimension (1: nbulk) :: x_ar,y_ar,
     &t_av_ar,tpop_ar,vflow_ar,mass_ar

      real*8 :: te,step,nperp
      integer i,iherm_l,i_deriv_l

      complex*16, dimension(1:3,1:3) :: K,dK_dnpar,dK_dnper
      complex*16, dimension(1:3,1:3,1:nbulk) :: 
     & dK_dx_ar,dK_dy_ar,dK_dt_ar
     
      complex*16 :: d
c-----externals
      real*8 b,tempe,x,y,tpoprho,vflowrho

      write(*,*)'in subroutine map_d_Ono z,r,phi',z,r,phi
      write(*,*)'npar,nperp_min,nperp_max',
     &npar,nperp_min,nperp_max

      bmod=b(z,r,phi) 
  
      step=(nperp_max-nperp_min)/dfloat(number_map_points_real_nperp-1)

      do i=1,nbulk
         x_ar(i)=x(z,r,phi,i)
	 y_ar(i)=y(z,r,phi,i)
         if(i.eq.1) y_ar(1)=-y_ar(1)
	 te=tempe(z,r,phi,i) ! kev
	 t_av_ar(i)=te*1000.d0      ! ev 
         tpop_ar(i)=tpoprho(rho,i)
         vflow_ar(i)=vflowrho(rho,i)
         mass_ar(i)=dmas(i)
        write(*,*)'i,x_ar(i),y_ar(i),t_av_ar(i),tpop_ar(i),vflow_ar(i)'
     &            ,i,x_ar(i),y_ar(i),t_av_ar(i),tpop_ar(i),vflow_ar(i)
      enddo

      iherm_l=2   !gives D fumction for complete tensor
      i_deriv_l=0 !no derivatives calculations

      write(*,*)'number_map_points_real_nperp',
     &           number_map_points_real_nperp

      do i=1,number_map_points_real_nperp
        nperp=nperp_min+(i-1)*step           
c        write(*,*)'i,nperp',i,nperp    
        
        call ono_tens(mass_ar,X_ar,Y_ar,t_av_ar,nbulk,
     &  npar,nperp,iherm_l,i_deriv_l,
     &  K,d,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

        nperp_ar(i)=nperp
        d_ar(i)=dreal(d)
        write(*,*)'i,nperp_ar(i),d_ar(i)',
     &             i,nperp_ar(i),d_ar(i)
      enddo

      name_param(1)='Y'
      name_param(2)='X'
      name_param(3)='N_parallel'
      name_param(4)='T'
      name_param(5)='r' 
      name_param(6)='z'
      name_param(7)='id_plot_disp'

      param(1)=y_ar(1)
      param(2)=x_ar(1)
      param(3)=npar
      param(4)=t_av_ar(1)
      param(5)=r  
      param(6)=z
      param(7)=id_plot_disp

      call plot1dt_param(nperp_ar,d_ar,0,0,number_map_points_real_nperp
     &,1,'linlin',0.d0,0.d0,
     &'Ono dispersion function',
     &'n_perpendicular','dispersion function',n_param,name_param,
     &param)

      return
      end
