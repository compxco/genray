c        **********************cninit**************************
c        *                        -                           *
c        * It solves the dispersion relation N=N(n_par)       *
c        * Then subroutine calculates the initial components  *
c        * of the refractive index N: cnz=N_z,cn=N_r,cm=M     *
c        ******************************************************
c
c------------------------------------------------------------------
c					                          !
c        input parameters				          !
c        z,r,phi,cnpar,cnteta,cnphi 		       	          !
c        if i_n_poloidal=4 cnpar=N_parallel will be calculated
c        inside this subroutine using cnteta,cnphi
c        istart from common block/ one/                           !
c        output parameters				          !
c        u(4)=cnz,u(5)=cnr,u(6)=cm 			          !
c        iraystop=1 is switch to stop the ray calculation         !
c------------------------------------------------------------------
c        it uses the following functions and subroutines          !
c        ias1r,b   ,y,x,gamma1,s,abc,hamilt1,                     !
c------------------------------------------------------------------
      subroutine cninit(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'      
      include 'grill.i'
      
      complex*16 hamilt_c !complex hamiltonian used for id=10
      save
      double complex hotnperp,cmplnper,relativistic_nperp
      double precision b
      external hotnperp,relativistic_nperp,b      

c-----for hot plasma roots
      double precision 
     &N_perp_root_ar(n_hot_roots_a)     !hot plasma roots

      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnteta*cnteta+cnphi*cnphi
      bmod=b(z,r,phi)

      if(istart.eq.2) then !grill conditions
        if (i_n_poloidal.eq.3) then
c----------refractive index is specified by N_parallel 
c          and ksi_nperp the angle between grad(psi) and_ N_perpendicular
c
c          set zero values for following variables
c          to get the solution of the dispersion function N_perp(N_parallel)
c          In this case N_perp=N

           cnteta=0.d0
           cnphi=0.d0
           cntang2=0.d0
        endif

        if (i_n_poloidal.eq.4) then
c---------refractive index is specified by N_toroidal and N_poloidal
c         calculation of the parallel refractive index component
c         N_parallel=cnpar
          gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
          b_teta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field
          cnpar=(cnphi*bphi+cnteta*b_teta)/bmod !N_parallel
          write(*,*)'cninit.f cnteta,cnphi,cnpar',cnteta,cnphi,cnpar
        endif

      endif !istart.eq.2
c-----------------------------------------------------------------
      iroot=0
c------------------------------------------------------------------
c     epsmode is the accuracy  for mode selection (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
      epsmode=0.0000001d0
c--------------------------------------------------------------------
c     the initial condition in plasma boundary
c     for cn2 from cnpar by using dispersion relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c------------------------------------------------------------------
c     Appleton - Hartree dispersion relation
c
c if 1
      if (id.eq.3) then

         call cninit3(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
         cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar**2)
         goto 111
      end if !id=3
c end if 1
c-----------------------------------------------------------------
c     initial condition in plasma boundary
c     for cn2 from cnpar by using dispersin relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
c if 0
      if (((id.eq.1).or.(id.eq.2)).or.(id.eq.16)) then

cyup         write(*,*)'cninit before cninit12 z,r,phi,cnpar,cnteta,cnphi',
cyup     &   z,r,phi,cnpar,cnteta,cnphi
cSAP111115          
         if(id.eq.16)then
           id_loc=id
           id=2
         else
           id_loc=id
         endif

         call cninit12(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)

cSAP111115          
         if (id_loc.eq.16)then
            id=id_loc
         endif


         if (iraystop.eq.1)then
           write(*,*)'the given conditions did not give root' 
           return
         endif
         cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar**2)
c         write(*,*)'cninit after call cnini12 cnpar,cnper,cn2',
c     &   cnpar,cnper,cnpar**2+cnper**2
c         write(*,*)'cninit cnz,cnr,cm',cnz,cnr,cm
c_test_begin
c         bmod=b(z,r,phi)
c         gam=gamma1(z,r,phi,cnz,cnr,cm)
c         epsham=hamilt1(z,r,phi,cnz,cnr,cm)
c         write(*,*)'epsham',epsham
c_test_end
         goto 111
      end if !id=1 or id=2
c end if 0

c if 4    
      if (((id.eq.4).or.(id.eq.5)).or.(id.eq.7))then
c-----------------------------------------------------------------
c       initial condition on plasma boundary
c       for cnper2 from cnpar by using dispersion relation
c       which was used in "Mazzucato's" code
c       Phys.Fluids 30(12),December 1987 p.3745
c-----------------------------------------------------------------
        ihermloc=iherm
        
        write(*,*)'in cninit ihermloc=',ihermloc,'cnpar=',cnpar
        ioptmaz=1 ! calculation of the estimation cnper from the cold plasma
        write(*,*)'in cninit before cnpermuz cnper,cnprim',cnper,cnprim
        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper,cnprim,ioptmaz)
        write(*,*)'in cninit aft cnpermuz cnper,cnprim',cnper,cnprim

      endif !id=4,5,7

      if(((id.eq.6).or.(id.eq.8)).or.(id.eq.9))then
c-----------------------------------------------------------------
c       initial condition on plasma boundary
c       for cnper2 from cnpar by using the hot dispersin relation
c-----------------------------------------------------------------
        if (i_look_roots.eq.2) then     
          call calculate_hot_nperp_roots(z,r,phi,cnpar,
     &    n_hot_roots,N_perp_root_ar)
          cnper=N_perp_root_ar(k_hot_root)
          cnprim=0.d0
        else
          ihermloc=iherm
          write(*,*)'cninit.f before hotnperp'
          cmplnper=hotnperp(z,r,phi,cnpar,cnteta,cnphi,K,iraystop)
          write(*,*)'cninit.f after hotnperp cmplnper',cmplnper
          if (iraystop.eq.1)then
           write(*,*)'the given conditions did not give hot root' 
          return
          endif
          cnper=dreal(cmplnper)
          cnprim=dimag(cmplnper)            
          write(*,*)'in cninit aft hotnperp cnper,cnprim',cnper,cnprim
        endif

      endif !id=5,8,9,16

c Smirnov 040304    
      if (id.eq.10) then
c-----------------------------------------------------------------
c       initial condition on plasma boundary
c       for cnper2 from cnpar by using dispersion relation
c       which was used in "Mazzucato's" code
c       Phys.Fluids 30(12),December 1987 p.3745
c       These initial condition will be used to choose the eigenvalue
c       of the dispersion tensor k_root
c       for the dispersion function D=Re(eigenvalue(k_root))
c-----------------------------------------------------------------

c-------calculate (cnper,cnprim) using mazzucato solver 
        ihermloc=iherm        

        write(*,*)'in cninit ihermloc=',ihermloc,'cnpar=',cnpar
        ioptmaz=1 ! calculation of the estimation cnper from the cold plasma
        write(*,*)'in cninit before cnpermuz cnper,cnprim',cnper,cnprim
        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper,cnprim,ioptmaz)
        write(*,*)'in cninit aft cnpermuz cnper,cnprim',cnper,cnprim

c-------determine k_root the number of the dispersion tensor eigenvalue
        hamilt_min=10.d15
        do k=1,3
 
          call hamilt_eigenvalue_muz(z,r,phi,cnpar,cnper,cnprim,k,
     &    hamilt,hamilt_c)

          if (dabs(hamilt_min).gt.dabs(hamilt)) then
            k_root=k      
            hamilt_min=hamilt
          endif
          write(*,*)'k,hamilt',k,hamilt
        enddo ! k      
        write(*,*)'cninit k_root',k_root 
c      stop 'cninit k_root'
      endif !id=10
c-----------------------------------------------------------------------    

      if ((id.eq.11).or.(id.eq.14)) then
c--------------------------------------------------------------------
c         id =11 Eric Nelson-Melby relativistic tensor
c         Dispersion function = Re(Det) 
c         id=14 Abhay Ram relativistic electron dispersion function
c----------------------------------------------------------------------
c       initial condition by using dispersion relation
c       which was used in "Mazzucato's" code 
c       and double complex function relativistic_nperp
c----------------------------------------------------------------------
        cmplnper=relativistic_nperp(z,r,phi,cnpar,cnteta,cnphi,K,
     &                              iraystop)
        cnper=dreal(cmplnper)
        cnprim=dimag(cmplnper)

        write(*,*)'cninit id14 z,r,phi,cnpar,cnteta,cnphi,cnper,cnprim',
     &  z,r,phi,cnpar,cnteta,cnphi,cnper,cnprim
      endif ! (id.eq.11).or.(id.eq.14)
c----------------------------------------------------------------        

      if (id.eq.12 .or. id.eq.15) then
c-----------------------------------------------------------------
c       initial condition on plasma boundary
c       for cnper2 from cnpar by using dispersion relation
c       which was used in "Mazzucato's" code
c       Phys.Fluids 30(12),December 1987 p.3745
c       These initial condition will be used to choose the eigenvalue
c       of the Eric Nelson Melby dispersion tensor k_root
c       for the dispersion function D=Re(eigenvalue(k_root))
c-----------------------------------------------------------------

c-------calculate (cnper,cnprim) using mazzucato solver 

c        ihermloc=iherm        
c        write(*,*)'in cninit ihermloc=',ihermloc,'cnpar=',cnpar
c        ioptmaz=1 ! calculation of the estimation cnper from the cold plasma
c        write(*,*)'in cninit befor cnpermuz cnper,cnprim',cnper,cnprim
c        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper,cnprim,ioptmaz)
c        write(*,*)'in cninit aft cnpermuz cnper,cnprim',cnper,cnprim
  
        cmplnper=relativistic_nperp(z,r,phi,cnpar,cnteta,cnphi,K,
     &                              iraystop)
        cnper=dreal(cmplnper)
        cnprim=dimag(cmplnper)

c-------determine k_root the number of the dispersion tensor eigenvalue
        hamilt_min=10.d15
        do k=1,3
 
c          call hamilt_eigenvalue_eric(z,r,phi,cnpar,cnper,cnprim,k,
c     &    hamilt,hamilt_c)  !old version

          call hamilt_eigenvalue_combined(z,r,phi,cnpar,cnper,cnprim,k,
     &    hamilt,hamilt_c)

          if (dabs(hamilt_min).gt.dabs(hamilt)) then
            k_root=k      
            hamilt_min=hamilt
          endif
          write(*,*)'k,hamilt',k,hamilt
        enddo ! k      
        write(*,*)'cninit k_root',k_root 
c      stop 'cninit k_root'
      endif !id=12
c--------------------------------------------------------------------
c Smirnov 040906    
      if (id.eq.13) then
c-----------------------------------------------------------------
c       initial condition on plasma boundary
c       for cnper2 from cnpar by using hot dispersin relation 
c       These initial condition will be used to choose the eigenvalue
c       of the Hot plasma dispersion tensor k_root
c       for the dispersion function D=Re(eigenvalue(k_root))
c-----------------------------------------------------------------

c-------calculate (cnper,cnprim) using hotnperp solver 
        ihermloc=2 
        cmplnper=hotnperp(z,r,phi,cnpar,cnteta,cnphi,K,iraystop)
c        write(*,*)'cninit.f after hotnperp cmplnper',cmplnper
        cnper=dreal(cmplnper)
        cnprim=dimag(cmplnper)
      
c-------determine k_root the number of the dispersion tensor eigenvalue
        hamilt_min=10.d15
        do k=1,3
 
          call hamilt_eigenvalue_hot(z,r,phi,cnpar,cnper,cnprim,k,
     &    hamilt,hamilt_c)

          if (dabs(hamilt_min).gt.dabs(hamilt)) then
            k_root=k      
            hamilt_min=hamilt
          endif
          write(*,*)'cninit k,hamilt',k,hamilt
        enddo ! k      
        write(*,*)'cninit k_root',k_root 
c      stop 'cninit k_root'
      endif !id=13
c--------------------------------------------------------------------

c-----calculations of initial values cnz,cnr,cm from cnper for id=4,5,6,7    
      cn2=cnper**2+cnpar**2   
      cnrho2=cn2-cnteta**2-cnphi**2
      if (cnrho2.lt.0.d0) cnrho2=1.d-16
      cnrho=dsqrt(cnrho2)

      write(*,*)'in cninit before cnzcnr cnper,cnpar,cnrho',
     & cnper,cnpar,cnrho
      write(*,*) 'cn=',dsqrt(cn2)
c---------------------------------------------------------------------
      call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)
c---------------------------------------------------------------------
cyup      write(*,*)'in cninit after cnzcnr cnz,cnr,cm',cnz,cnr,cm
      cm=cnphi*r
      gam=gamma1(z,r,phi,cnz,cnr,cm)
cyup      write(*,*)'cninit.f before ddd=hamilt1'
      ddd=hamilt1(z,r,phi,cnz,cnr,cm)
cyup      write(*,*)'in cninit ddd=',ddd


  111 continue

c----------------------------------------------------------------
c     for the grill conditions it will use the different grill types
c---------------------------------------------------------------
cyup      write(*,*)'cninit istart,i_n_poloidal',istart,i_n_poloidal
      if(istart.eq.2) then !grill conditions
         if (i_n_poloidal.eq.2) then !input N_parallel, N_poloidal
c            calculate N_phi,N_theta,N_rho
c-----------------------------------------------------------------------
c            dpdrd=dpsidr dpdzr=d psidr were calculated by b(z,r,phi)
c            They are in  common/one/
c---------------------------------------------------------------------
             gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
             b_teta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field
             cnteta_full=n_theta_pol
             write(*,*)'cninit i_n_poloidal,n_theta_pol',
     &       i_n_poloidal,n_theta_pol
             write(*,*)'cnteta_full',cnteta_full

             alpha_teta=(cnteta_full-cnpar*b_teta/bmod)/cnper
             write(*,*)'cnpar,b_teta,bmod,cnper,alpha_teta',
     &       cnpar,b_teta,bmod,cnper,alpha_teta
     
             cnphi_full=cnpar*bphi/bmod-alpha_teta*cnper*b_teta/bphi
             write(*,*)'cnphi_full',cnphi_full
c_test
             cnteta_full=cnpar*b_teta/bmod+cnper*alpha_teta
             write(*,*)'2 cnteta_full=',cnteta_full

             arg= 1.d0-(alpha_teta*bmod/bphi)**2
             if (arg.lt.0.d0) then
               WRITE(*,*)'cninit.f 1.d0-(alpha_teta*bmod/bphi)**2<0'
               WRITE(*,*)'cninit.f change n_theta_pol'
               STOP
             endif

             cnrho_full=cnper*dsqrt(arg)
             write(*,*)'cnteta_full,cnphi_full,cnrho_full,cn**2',
     &       cnteta_full,cnphi_full,cnrho_full,
     &       cnteta_full**2+cnphi_full**2+cnrho_full**2

         endif !i_n_poloidal=2
cSAP090504
cSAP091026
          if (i_n_poloidal.eq.3) then !input N_parallel, ksi_nperp
c------------calculate N_phi,N_theta,N_rho
             gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
             b_teta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field
             rad_ksi_nperp=ksi_nperp*pi/180.d0 !transfrm degrees to radians
             cnteta_full=(cnpar*b_teta+cnper*bphi*dsin(rad_ksi_nperp))
     &                    /bmod
             cnphi_full=(cnpar*bphi-cnper*b_teta*dsin(rad_ksi_nperp))
     &                    /bmod   
             cnrho_full=-cnper*dcos(rad_ksi_nperp)
             write(*,*)'cninit.f i_n_poloidal=3 cnrho_full',cnrho_full
          endif ! i_n_pol=3

         if (i_n_poloidal.eq.4) then !input N_toroidal and N_poloidal
c-------------calculate N_phi,N_theta,N_rho
              cnteta_full=cnteta
              cnphi_full=cnphi
              write(*,*)'cninit i_n_poloidal=4 cnteta,cnphi,cnper,cnpar'
     &                   ,cnteta,cnphi,cnper,cnpar
     
              if (cnper.ne.0.d0) then
                 sin_ksi=(cnphi*b_teta-cnteta*bphi)/(bmod*cnper)
                 if (dabs(sin_ksi).le.1.d0) then
                    cos_ksi=dsqrt(1.d0-sin_ksi**2)
                    cnrho_full=cnper*cos_ksi*i_vgr_ini
                 else
                    write(*,1010)
 1010               format('in cninit.f i_n_poloidal=4 case',/,
     &              'dabs(sin_ksi)>1,it is imposible to find N_rho')
                    iraystop=1
                    return                
                 endif
              else
                 if ((cnphi/bphi).ne.(cnteta/b_teta)) then
                    write(*,1000)
 1000               format('in cninit.f i_n_poloidal=4 case',/,
     &             'cnper=0 and it is impossible to find n_parallel')
                   iraystop=1
                   return
                 else
                   cnrho_full=0.d0
                 endif
              endif
              
              write(*,*)'cninit.f i_n_poloidal=4 cnrho_full',cnrho_full
         endif ! i_n_pol=4

cSAP090504
cSAP091026 
         if(i_n_poloidal.ne.1) then
c         if((i_n_poloidal.ne.1).and.(i_n_poloidal.ne.3)) then
             call cnzcnr(z,r,phi,cnteta_full,cnphi_full,cnrho_full,
     &       cnz,cnr,cm)
             if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
             write(*,*)'cninit end grill condition'
             write(*,*)'cnteta_full,cnphi_full,cnrho_full,cn**2',
     &       cnteta_full,cnphi_full,cnrho_full,
     &       cnteta_full**2,cnphi_full**2+cnrho_full**2
             write(*,*)'cnz**2+cnr**2+(cm/r)**2',
     &       cnz**2+cnr**2+(cm/r)**2
             write(*,*)'cnz,cnr,cm',cnz,cnr,cm
             write(*,*)'(cnz*bz+cnr*br+cm*bphi/r)/bmod',
     &        (cnz*bz+cnr*br+cm*bphi/r)/bmod
             endif ! outprint
c_test_begin
c         bmod=b(z,r,phi)
c         gam=gamma1(z,r,phi,cnz,cnr,cm)
c         epsham=hamilt1(z,r,phi,cnz,cnr,cm)
c         write(*,*)'epsham',epsham
c_test_end
         endif

      endif !grill conditions


      return
      end
c_____________________________________________________________
c        **********************cnzcnr**************************
c        This subroutine calculates the initial value		 *
c        cnz,cnr,cm					         *
c        It directs the wave into or out the plasma              *
c        the input parameter i_vgr_ini is in common /one/        *
c        i_vgr_ini =+1 the wave is directed into the plasma      *
c                      (in the initial point)                    *
c                  =-1 the wave is directed out the plasma       *
c-----------------------------------------------------------------
c      * input parameters:z,r,phi ,cnteta,cnphi,cnrho	      	 *
c      *                  cnrho is directed inside the plasma	 *
c________________________________________________________________
c     * output parameters: cnz,cnr,cm -refractive index components*
c-----------------------------------------------------------------
      subroutine cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)
      implicit double precision (a-h,o-z)
c-----------------------------------------------------------------
c	           z       r      phi
c       e_phi= {   0   ,   0    ,   1    }
c       e_psi= { dpsidz, dpsidr ,   0    }/mod(grad(psi))
c       e_teta={ dpsidr,-dpsidz ,   0    }/mod(grad(psi))
c
c-----------------------------------------------------------------
      include 'param.i'
      include 'one.i'
      include 'grill.i'
      dimension u(6),deru(6)
c---------------------------------------------------------------------
c     dpdrd=dpsidr dpdzr=dpsidr were calculated by b(z,r,phi)
c     They are in  common/one/
c---------------------------------------------------------------------
      gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
c-------------------------------------------------------------------
c     the initialization of cnrho 
      cirho=1.d0
c------------------------------------------------------------------
10    continue
cyup      write(*,*)'cninit.f in cnzcnr cnteta,cirho,cnrho',
cyup     &cnteta,cirho,cnrho
c     write(*,*)'cninit.f in cnzcnr dpdzd,dpdrd,gradpsi',
c    &dpdzd,dpdrd,gradpsi
      cnz=(cnteta*dpdrd-cirho*cnrho*dpdzd)/gradpsi
      cnr=(-cnteta*dpdzd-cirho*cnrho*dpdrd)/gradpsi
      cm=cnphi*r
cyup      write(*,*)'cninit.f in cnzcnr cnz,cnr,cm',cnz,cnr,cm
cSm030515
c      if((i_n_poloidal.eq.1).or.(i_n_poloidal.eq.2)) then
cSm040415
c     I tried it for id=4 ECR case when rside1 used the derivatives
c     by the trajectory d/dl (not by time d/dt).
c     In some case d/dl runs the ray to the plasma edge.
c     As I understand this changing of the cnrho sign is not correct
c     fo EC launch. It changes the initial NR and NZ signes. 
      if (id.eq.10) return

      if(((i_n_poloidal.eq.1).or.(i_n_poloidal.eq.2)).
     &    or.(istart.eq.1)) then 

  
c*****************************************************************
c        determination cirho to obey the situation
c        when groop velocity is directed inside or outside  the plasma
c        in the initial point
c-----------------------------------------------------------------
         u(1)=z
         u(2)=r
         u(3)=phi
         u(4)=cnz
         u(5)=cnr
         u(6)=cm
cyup        write(*,*)'cninit.f in cnzcnr before rside1'    
         call rside1(0.d0,u,deru) 
cyup        write(*,*)'cninit.f in cnzcnr after rside1 deru',deru        
c----------------------------------------------------------------
c        cmultpl is the scalar multiplication V_groop*grad(psi)
c        gradient(psi) is directed outside the plasma
c----------------------------------------------------------------
         cmultpl=dpdzd*deru(1)+dpdrd*deru(2)
cyup        write(*,*)'cninit in cnzcnr cmultpl,i_vgr_ini',cmultpl,i_vgr_ini

         if ((cmultpl*i_vgr_ini).gt.0.d0) then
c-----------the poloidal direction of the group velocity is opposite
c           to the direction determined by the parameter i_vgr_ini     
c           We change the sign of the poloidal group velocity 
            cirho=-1.d0
            go to 10
         end if
c****************************************************************
c         end irho determination
      endif !i_n_poloidal =1 or =2
c----------------------------------------------------------------
cyup      write(*,*)'cninit.f at end, cnzcnr cnz,cnr,cm:',cnz,cnr,cm
      return
      end ! cnzcnr

c     **********************nsolv**************************
c     *                        -                           *
c     * It solves the dispersion relation N=N(n_par)       *
c     ******************************************************
c-----------------------------------------------------------
c     							   !
c     input parameters					   !
c     z,r,phi,cnz,cnr,cm      				   !
c     output parameters					   !
c     cn2p,cn2m - roots of the dispersion relation =N**2   !
c-----------------------------------------------------------
c     it uses the following functions and subroutines      !
c     b ,y,x,s,abc,hamilt1                                 !
c-----------------------------------------------------------
      subroutine nsolv(z,r,phi,cnz,cnr,cm,cn2p,cn2m)
c      implicit double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
c-----input
      real*8 z,r,phi,cnz,cnr,cm
c-----output
      real*8 cn2p,cn2m

c-----locals
      real*8 cnpar,cnper,cn2,xi,yi,pyp,
     &f1e,f0e,g1e,g0e,w1e,w0e,dele,fd,gd,wd,detin,delib,
     &f0b,f1b,g0b,g1b,w1b,w0b,
     &xe,xb,ye,yb,
     &s1,s2,s3,s4,s6,s7,cnprim,cnpar2,
     &px,px2,pypb,pyme,pyme2,py2,py4,pype
     
      real*8 sign_del
      integer ibmx !local

      integer ioxmold,ioptmaz,iraystop        

c-----externals
      real*8 b,x,y

      pi=4*datan(1.d0)
      bmod=b(z,r,phi)
      cnpar=(cnz*bz+cnr*br+cm*bphi/r)/bmod
      cn2=cnz*cnz+cnr*cnr+cm*cm/(r*r)
c     write(*,*)'in nsolv cn2',cn2,'cnpar',cnpar
      cnpar2=cnpar*cnpar
c     solution of the dispersion relation
c     for cn2 from cnpar by using dispersin relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)

      ! Note: a=A*delta, b=B*delta, c=C*delta (where delta=1-Y)
      ! and sqrt(det)= sqrt(B*B-4*A*C)*|delta|
      ! in   N^2 = (-b +ioxm*sqrt(b*b-4*a*c))/(2a)
      
      ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
      delib= 1.d0-y(z,r,phi,ibmx) ! 1-Y (electrons or ions)
      sign_del=1.d0 !sign(1.d0,delib) ! in nsolv

c------------------------------------------------------------------
c     Applton - Hartry dispersion relation
c
c if 1
      if (id.eq.3) then
         xi=x(z,r,phi,1)
         yi=y(z,r,phi,1)
         py2=yi*yi
         py4=py2*py2
         px=1.-xi
         px2=px*px
c------------------------------------------------------------------
         pyp=xi/(1.+yi)
         f1e=1.d0
         f0e=-pyp
         g1e=cnpar2*(-xi)+pyp+xi-2.
         g0e=cnpar2*pyp+xi*(1.-pyp)+pyp*(1.-xi)
         w1e=cnpar2*(xi-pyp)+(1.-pyp)*(1.-xi)
         w0e=cnpar2*(pyp*(1-xi)-xi*(1.-pyp))-(1.-xi)*(1-pyp)*xi
         dele=1.-yi
         fd=f1e*dele+f0e
         gd=g1e*dele+g0e
         wd=w1e*dele+w0e  !==(1-Y)*w
         detin=gd**2-4.*fd*wd
         if (detin.lt.0d0) then
            write(*,*)' 1 in nsolv detin  less then zero '
            return
	   end if
	   cn2p=(-gd+sign_del*dsqrt(detin))/(2.*fd) ! N^2(Npar) using (f,g,w)
	   cn2m=(-gd-sign_del*dsqrt(detin))/(2.*fd) ! N^2(Npar) using (f,g,w)
	   WRITE(*,*)'Apl nsolv cn2p,cn2m',cn2p,cn2m
           if((cn2m.lt.0.d0).and.(cn2p.lt.0d0)) then
            write(*,*)'in cninit2 two roots of the dispersion < 0'
            write(*,*)'cn2m,cn2p',cn2m,cn2p
            write(*,*)'the given wave can not exist in plasma'
	    iraystop=1
	    return
	   endif
c-------------------------------------------------------------------
      end if
c end if 1
c-----------------------------------------------------------------
c     solution of the dispersion relation
c     for cn2 from cnpar by using the dispersion relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
c if 0
      if ((id.eq.1).or.(id.eq.2)) then
        call s(z,r,phi,s1,s2,s3,s4,s6,s7)

c       ib=1 electron resonance condition may be present in plasma
c  if 2
        if (ib.eq.1) then
c          write(*,*)'cold plasma ib=1 '
           xe=x(z,r,phi,1)
           ye=y(z,r,phi,1)
	     pyp=xe/(1.+ye)
	     dele=1.-ye
	     f1e=s7
	     f0e=-pyp
	     g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
	     g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4
	     w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4
             w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe
	     fd=f1e*dele+f0e
	     gd=g1e*dele+g0e
	     wd=w1e*dele+w0e  !==(1-Y)*w
	     detin=gd**2-4.*fd*wd
	     if (detin.lt.0d0) then
	        write(*,*)' 2 in nsolv detin  less then zero '
	        return
	     end if
	     cn2p=(-gd+sign_del*dsqrt(detin))/(2.*fd)
	     cn2m=(-gd-sign_del*dsqrt(detin))/(2.*fd)

             if((cn2m.lt.0.d0).and.(cn2p.lt.0d0)) then
               write(*,*)'in cninit2 two roots of the dispersion < 0'
               write(*,*)'cn2m,cn2p',cn2m,cn2p
               write(*,*)'the given wave can not exist in plasma'
	       iraystop=1
	       return
	     endif

	  endif ! ib=1
c end if 2
c
c       ib.gt.1 ion (species i=ib) resonance condition may be present in plasma
c  if 3
        if (ib.gt.1) then
c          write(*,*)'cold plasma ib .gt.1  '
           ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
           xb=x(z,r,phi,ibmx) 
           yb=y(z,r,phi,ibmx) ! i=ib species 
           xe=x(z,r,phi,1)
           ye=y(z,r,phi,1)
	   pype=xe/(1.+ye)
	   pyme=xe/(1.-ye)
	   pyme2=pype/(1.-ye)
	   pypb=xb/(1.+yb)
	   delib=1.d0-yb
	   f1b=(s1-pyme2)
	   f0b=-pypb

	   g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
	   g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

	   w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1	       s4*(s2-pype)*(s3-pyme)
           w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
	   fd=f1b*delib+f0b
	   gd=g1b*delib+g0b
	   wd=w1b*delib+w0b
	   detin=gd**2-4.*fd*wd
	   if (detin.lt.0d0) then
	      write(*,*)' 3 in nsolv detin  less then zero '
	      return
	   end if
	   cn2p=(-gd+sign_del*dsqrt(detin))/(2.*fd)
	   cn2m=(-gd-sign_del*dsqrt(detin))/(2.*fd)
           if((cn2m.lt.0.d0).and.(cn2p.lt.0d0)) then
            write(*,*)'in cninit2 two roots of the dispersion < 0'
            write(*,*)'cn2m,cn2p',cn2m,cn2p
            write(*,*)'the given wave can not exist in plasma'
	    iraystop=1
	    return
	   endif
        end if ! ib>1
c end if 3
      end if
c end if 0
c if 4
      if (id.eq.4) then
c-----------------------------------------------------------------
c       initial condition in plasma boundary
c       for cnper2 from cnpar by using dispersion relation
c       which was used in "Mazzucato's" code
c       Phys.Fluids 30(12),December 1987 p.3745
c-----------------------------------------------------------------
        write(*,*)'in cninit nsolv iherm=',iherm,'cnpar=',cnpar
        ioxmold=ioxm
        ioxm=1
        ioptmaz=1 ! estimation of cnper from the cold plasma
        call cnpermuz(cnpar,iherm,z,r,phi,cnper,cnprim,ioptmaz)
        write(*,*)'in cninit nsolv 1 cnper,cnprim',cnper,cnprim
        cn2p=cnper*cnper+cnpar2
        ioxm=-1
        ioptmaz=1 !estimation of cnper from the cold plasma
        call cnpermuz(cnpar,iherm,z,r,phi,cnper,cnprim,ioptmaz)
        write(*,*)'in cninit nsolv 2 cnper,cnprim',cnper,cnprim
        cn2m=cnper*cnper+cnpar2
        ioxm=ioxmold
	go to 111
      end if
c end if 4
  111 continue
      return
      end ! nsolv

c        **********************npernpar************************
c        *                        -                           *
c        * It solves the dispersion relation Nper=N(n_par)    *
c        * Then subroutine calculates the initial components  *
c        * of the refractive index  cnz,cnr,cm                *
c        ******************************************************
c
c---------------------------------------------------
c     							   !
c     input parameters					   !
c     z,r,phi,cnpar 	                                   !
c     output parameters:cnper2p,cnper2m	                   !
c---------------------------------------------------
c     it uses the following functions and subroutines      !
c     ias1r,b   ,y,x,gamma1,s,abc,hamilt1,                 !
c-----------------------------------------------------------
      subroutine npernpar(z,r,phi,cnpar,cnper2p,cnper2m)
c      implicit double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
c-----input
      real*8 z,r,phi,cnpar
c-----output
      real*8 cnper2p,cnper2m

c-----locals
      real*8 cnper,cnpar2,cnpar4,xi,yi,py2,py4,px,px2,pyp,
     &f1e,f0e,g1e,g0e,w1e,w0e,dele,fd,gd,wd,gnew,wnew,
     &detin,delib,
     &s1,s2,s3,s4,s6,s7,xe,ye,xb,yb,
     &pype,pyme,pyme2,pypb,
     &f1b,f0b,g1b,g0b,w1b,w0b,
     &cnprim
      real*8 sign_del
      integer ibmx !local

      integer
     &ioxmold,ioptmaz
c-----externals
      real*8 b,x,y

      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cnpar4=cnpar2*cnpar2
      bmod=b(z,r,phi)
c-------------------------------------------------------------
c     calculations of  cnrer2p and cnper2m
c     from cnpar2 by using the dispersin relation
c     f*cnper**4 +(2f*cnpar**2+g)*cnper**2+(f*cnpar**4+g*cnpar**2+w)=0
c     f*cnper**4 +gnew*cnper**2+wnew=0
c     gnew=2fcnpar**2+g, wnew=f*cnpar**4+g*cnpar**2+w
c     in the following form cnper**2=(-gnew+,-*sqrt(gnew**2-4*f*wnew))/(2*f)

      ! Note: a=A*delta, b=B*delta, c=C*delta (where delta=1-Y)
      ! and sqrt(det)= sqrt(B*B-4*A*C)*|delta|
      ! in   N^2 = (-b +ioxm*sqrt(b*b-4*a*c))/(2a)
      ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
      delib= 1.d0-y(z,r,phi,ibmx) ! 1-Y (electrons or ions)
      sign_del=1.d0 !sign(1.d0,delib) ! in npernpar
c------------------------------------------------------------------
c     Appleton - Hartry dispersion relation
c
c if 1
      if (id.eq.3) then
         xi=x(z,r,phi,1)
         yi=y(z,r,phi,1)
         py2=yi*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
c------------------------------------------------------------------
         pyp=xi/(1.d0+yi)
         f1e=1.d0
         f0e=-pyp
         g1e=cnpar2*(-xi)+pyp+xi-2.
         g0e=cnpar2*pyp+xi*(1.d0-pyp)+pyp*(1.d0-xi)
         w1e=cnpar2*(xi-pyp)+(1.-pyp)*(1.-xi)
         w0e=cnpar2*(pyp*(1-xi)-xi*(1.-pyp))-(1.d0-xi)*(1.d0-pyp)*xi
         dele=1.d0-yi
         fd=f1e*dele+f0e
         gd=g1e*dele+g0e
         wd=w1e*dele+w0e  !==(1-Y)*w
c new coefficients
         gnew=gd+2.d0*fd*cnpar2
         wnew=wd+gd*cnpar2+fd*cnpar4
         gd=gnew
         wd=wnew

         detin=gd**2-4.*fd*wd
         if (detin.lt.0d0) then
            write(*,*)' 1 in npernpar detin  less then zero '
            return
         endif
         cnper2p=(-gd+sign_del*dsqrt(detin))/(2.d0*fd)
         cnper2m=(-gd-sign_del*dsqrt(detin))/(2.d0*fd)
c         WRITE(*,*)'Aplt cnpernpar cnper2p,cnper2m',cnper2p,cnper2m
      end if
c end if 1
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
c if 0
      if ((id.eq.1).or.(id.eq.2)) then
        call s(z,r,phi,s1,s2,s3,s4,s6,s7)

c       ib=1 electron resonance condition may be
c  if 2
        if (ib.eq.1) then
c         write(*,*)'cnint cold plasma ib=1 '
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)
          pyp=xe/(1.d0+ye)
          dele=1.d0-ye
          f1e=s7
          f0e=-pyp
          g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
          g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4
          w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4
          w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe
          fd=f1e*dele+f0e
          gd=g1e*dele+g0e
          wd=w1e*dele+w0e  !==(1-Y)*w
c new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew

          detin=gd**2-4.d0*fd*wd
          if (detin.lt.0d0) then
             write(*,*)' 2 in npenpar detin  less then zero'
             cnper2p=-1.d0
             cnper2m=-1.d0
             return
	  end if
          cnper2p=(-gd+sign_del*dsqrt(detin))/(2.d0*fd)
          cnper2m=(-gd-sign_del*dsqrt(detin))/(2.d0*fd)
c          write(*,*)'in cninit.f  1 npernpar cnper2p,cnper2m',
c     .    cnper2p,cnper2m
          goto 111
        end if
c end if 2
c
c     ib.gt.1 iones resonance condition may be
c  if 3
        if (ib.gt.1) then
c          write(*,*)'cold plasma ib .gt.1  '
           ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
           xb=x(z,r,phi,ibmx)
           yb=y(z,r,phi,ibmx)
           xe=x(z,r,phi,1)
           ye=y(z,r,phi,1)
           pype=xe/(1.d0+ye)
           pyme=xe/(1.d0-ye)
           pyme2=pype/(1.d0-ye)
           pypb=xb/(1.d0+yb)
           delib=1.d0-yb
	   f1b=(s1-pyme2)
	   f0b=-pypb

	   g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
           g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

	   w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1	     s4*(s2-pype)*(s3-pyme)
           w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
           fd=f1b*delib+f0b ! see below
           gd=g1b*delib+g0b ! see below
           wd=w1b*delib+w0b ! see below
c new coefficients
           gnew=gd+2.d0*fd*cnpar2
           wnew=wd+gd*cnpar2+fd*cnpar4
           gd=gnew
           wd=wnew

           detin=gd**2-4.d0*fd*wd
           if (detin.lt.0d0) then
              write(*,*)' 3 in dinit detin  less then zero '
              return
           end if
           cnper2p=(-gd+sign_del*dsqrt(detin))/(2.d0*fd)
           cnper2m=(-gd-sign_del*dsqrt(detin))/(2.d0*fd)
c        write(*,*)'in npernpar ib.qt.1 cnper2p,cnper2m',cnper2p,cnper2m
           goto 111
        end if
c end if 3
      end if
c end if 0
c if 4
      if ((id.eq.4).or.(id.eq.5)) then
c-----------------------------------------------------------------
c        for cnper2 from cnpar by using dispersin relation
c        which was used in "Mazzucato's" code
c        Phys.Fluids 30(12),December 1987 p.3745
c-----------------------------------------------------------------
         ioxmold=ioxm
         ioxm=1
         ioptmaz=1 !estimation of cnper from the cold plasma
         call cnpermuz(cnpar,iherm,z,r,phi,cnper,cnprim,ioptmaz)
         cnper2p=cnper*cnper
         ioxm=-1
         ioptmaz=1 !estimation of  cnper from the cold plasma
         call cnpermuz(cnpar,iherm,z,r,phi,cnper,cnprim,ioptmaz)
         cnper2m=cnper*cnper
c         write(*,*)'in cninit.f npernpar cnper2p,cnper2m',
c     +   cnper2p,cnper2m
         ioxm=ioxmold
         go to 111
      end if
c end if 4
  111 continue

      return
      end ! npernpar
      
c        **********************cninit3*************************
c        *                        -                           *
c        * It solves the dispersion relation N=N(n_par)       *
c        * for Appleton-Hartree disperstion relation
c        * Then subroutine calculates the initial components  *
c        * of the refractive index  cnz,cnr,cm        	      *
c        ******************************************************
c
c------------------------------------------------------------------
c					                          !
c        input parameters				          !
c        z,r,phi,cnpar,cnpar,cnphi 		       	          !
c        istart from common block/ one/                           !
c        output parameters				          !
c        u(4)=cnz,u(5)=cnr,u(6)=cm 			          !
c        iraystop=1 end ray calculation                           !
c------------------------------------------------------------------
c        it uses the following functions and subroutines          !
c        ias1r,b   ,y,x,gamma1,s,abc,hamilt1,                     !
c------------------------------------------------------------------
      subroutine cninit3(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      save
      double complex cmplnper      

      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnteta*cnteta+cnphi*cnphi
      bmod=b(z,r,phi)

c-----------------------------------------------------------------
      iroot=0
c------------------------------------------------------------------
c     epsmode is the accuracy  for selection mode (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
      epsmode=0.0000001d0
c--------------------------------------------------------------------
c     initial condition in plasma boundary
c     for cn2 from cnpar by using dispersion relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c------------------------------------------------------------------
c     Appleton - Hartree dispersion relation
c
      xi=x(z,r,phi,1)
      yi=y(z,r,phi,1)
      py2=yi*yi
      py4=py2*py2
      px=1.d0-xi
      px2=px*px
c------------------------------------------------------------------
      pyp=xi/(1.d0+yi)
      f1e=1.d0
      f0e=-pyp
      g1e=cnpar2*(-xi)+pyp+xi-2.d0
      g0e=cnpar2*pyp+xi*(1.-pyp)+pyp*(1.d0-xi)
      w1e=cnpar2*(xi-pyp)+(1.d0-pyp)*(1.d0-xi)
      w0e=cnpar2*(pyp*(1.d0-xi)-xi*(1.d0-pyp))-(1.d0-xi)*(1-pyp)*xi
      dele=1.d0-yi
      fd=f1e*dele+f0e
      gd=g1e*dele+g0e
      wd=w1e*dele+w0e  !==(1-Y)*w
      detin=gd**2-4.d0*fd*wd
      if (detin.lt.0d0) then
	 write(*,*)' 3 in cninit detin  less then zero '
	 iraystop=1
	 return
      end if

      ! Note: a=A*delta, b=B*delta, c=C*delta (where delta=1-Y)
      ! and sqrt(det)= sqrt(B*B-4*A*C) * |delta|
      delib= dele  ! 1-Y (here electrons)
      sign_del=1.d0 !sign(1.d0,delib) ! in cninit3

      cn2p=(-gd+sign_del*dsqrt(detin))/(2.d0*fd) ! corr to ioxm=+1
      cn2m=(-gd-sign_del*dsqrt(detin))/(2.d0*fd) ! corr to ioxm=-1
      WRITE(*,*)'apl cninit cn2p,cn2m',cn2p,cn2m
      if((cn2m.lt.0.d0).and.(cn2p.lt.0d0)) then
            write(*,*)'in cninit2 two roots of the dispersion < 0'
            write(*,*)'cn2m,cn2p',cn2m,cn2p
            write(*,*)'the given wave can not exist in plasma'
	    iraystop=1
	    return
      endif
10    iroot=iroot+1
      if(iroot.eq.1) then
        cn2=cn2p
	if(cn2.lt.cntang2)then
           write(*,*)'in cninit2 cn2p.lt.cntang2'
           go to 10
        end if
      else
        cn2=cn2m
        if(cn2.lt.cntang2)then
           write(*,*)'in cninit2 cn2m.lt.cntang2'
           write(*,*)'the given mode can not exist in plasma'
           iraystop=1
           return
        end if
      endif
c     write(*,*)'in cninit cn2=',cn2
      cnrho2=cn2-cnteta**2-cnphi**2
      cnrho=dsqrt(cnrho2)
      cnperp=cnrho
cSm050826
      cnperp=dsqrt(cn2-cnpar2)
c------------------------------------------------------------------
      call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)
c-------------------------------------------------------------------
      cm=cnphi*r
      gam=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam)
      dc=dcos(gam)
      ds2=ds*ds
      dc2=dc*dc
      ds4=ds2*ds2
c--------------------------------------------------------------------
c     controle that cn2 and gam are the solution of the dispersion
c     relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c--------------------------------------------------------------------
      sqrdet=dsqrt(py4*ds4+4.d0*py2*px2*dc2)
      pz= 2.d0*px -py2*ds2 +ioxm*sqrdet
      cn2new=1.d0-2.d0*xi*px/pz
c     write(*,*)'cn2new=',cn2new
      if (iroot.eq.1) then
	dnp=dabs(cn2-cn2new)
	if (dnp.gt.epsmode)then
	   goto 10
	end if
      else
	dnm=dabs(cn2-cn2new)
	if (dnm.gt.epsmode)then
           write(*,*)'the given mode can not exist in plasma'
	   iraystop=1
	   return
	end if
      end if
      goto 111
      
  111 continue
      return
      end !  cninit3

c        **********************cninit12_n_gam******************
c        *                        -                           *
c        * It solves the dispersion relation 
c        & N=N(n_par)=N(gam,ioxm)                             *
c        * for cold plasma for given ioxm                     *
c        * Then subroutine calculates the initial components  *
c        * of the refractive index  cnz,cnr,cm         	      *
c        ******************************************************
c
c------------------------------------------------------------------
c					                          !
c        input parameters				          !
c        z,r,phi,cnpar,cnteta,cnphi 		       	          !
c        istart from common block/ one/                           !
c        output parameters				          !
c        u(4)=cnz,u(5)=cnr,u(6)=cm 			          !
c        iraystop=1 end ray calculation                           !
c------------------------------------------------------------------
c        it uses the following functions and subroutines          !
c        ias1r,b   ,y,x,gamma1,s,abc,hamilt1,                     !
c------------------------------------------------------------------
      subroutine cninit12_n_gam(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      save
      double complex cmplnper      
      
      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnteta*cnteta+cnphi*cnphi
      bmod=b(z,r,phi)
      
      ! Note: a=A*delta, b=B*delta, c=C*delta (where delta=1-Y)
      ! and sqrt(det)= sqrt(B*B-4*A*C)*|delta|
      !in    N^2 = (-b +ioxm*sqrt(b*b-4*a*c))/(2a)
      ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
      delib= 1.d0-y(z,r,phi,ibmx) ! 1-Y (electrons or ions)
      sign_del=1.d0 !sign(1.d0,delib) ! in cninit12_n_gam
      
c      write(*,*)'cninit12_n_gam ioxm=',ioxm
c      write(*,*)'z,r,phi,cnpar,cnteta,cnphi',z,r,phi,cnpar,cnteta,cnphi
c-----------------------------------------------------------------
      iroot=0
c------------------------------------------------------------------
c     epsmode is the accuracy  for selection mode (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
      epsmode=0.0000001d0
      epsmode=1.d-7
c--------------------------------------------------------------------
c     initial condition in plasma boundary
c     for cn2 from cnpar by using dispersion relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
c if 0
c      write(*,*)'cninit12 z,r,phi',z,r,phi
      if ((id.eq.1).or.(id.eq.2)) then
        call s(z,r,phi,s1,s2,s3,s4,s6,s7)

c       ib=1 electron resonance condition may be
c  if 2
        if (ib.eq.1) then
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)
c          write(*,*)'cninit12 xe,ye',xe,ye
	    pyp=xe/(1.d0+ye)
	    dele=1.d0-ye
	    f1e=s7  
	    f0e=-pyp 
	    g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
	    g0e=cnpar2*pyp +xe*(s6-pyp)+pyp*s4
            w1e=cnpar2*(-s7*s4+(s6-pyp)*s3) +(s6-pyp)*s3*s4
            w0e=cnpar2*(pyp*s4-xe*(s6-pyp)) -s4*(s6-pyp)*xe
	    fd=f1e*dele+f0e !==(1-Ye)*eps1  
	    gd=g1e*dele+g0e !==
	    wd=w1e*dele+w0e !==
	    detin=gd**2-4.d0*fd*wd
c            write(*,*)'cninit12 fd,gd,wd,detin',fd,gd,wd,detin
	    if (detin.lt.0d0) then
	       write(*,*)' 2 in dinit detin  less then zero '
	       iraystop=1
	       return
	    end if
	    cn2p=(-gd+sign_del*dsqrt(detin))/(2.d0*fd)
	    cn2m=(-gd-sign_del*dsqrt(detin))/(2.d0*fd)
c	    write(*,*)'in cninit fd,gd,wd,detin'
	    !write(*,*)fd,gd,wd,detin
cyup	    write(*,*)'cninit12_n_gam cn2p,cn2m',cn2p,cn2m
cyup            write(*,*)'cninit12_n_gam cn2p-cnpar2,cn2m-cnpar2',
cyup     .      (cn2p-cnpar2),(cn2m-cnpar2)

ctest angle
            if (cn2p.gt.0.d0)then
               if(cnpar2.le.cn2p) then
                 dc_l=cnpar/dsqrt(cn2p)
                 gam_l=dacos(dc_l)
cyup                 write(*,*)'in cninit12_n_gam cn2p,gam_l',cn2p,gam_l
                else
cyup                 write(*,*)'in cninit12_n_gam cnpar2>cn2p',cnpar2,cn2p
                endif
            endif  

            if (cn2m.gt.0.d0)then
              if(cnpar2.le.cn2m) then
                dc_l=cnpar/dsqrt(cn2m)
                gam_l=dacos(dc_l)
cyup                write(*,*)'in cninit12_n_gam cn2m,gam_l',cn2m,gam_l
              else
cyup                write(*,*)'in cninit12_n_gam cnpar2>cn2m',cnpar2,cn2m
              endif
            endif


            if((cn2m.lt.0.d0).and.(cn2p.lt.0d0)) then
cyup             write(*,*)'cninit12_n_gam two roots of the dispersion < 0'
cyup             write(*,*)'cn2m,cn2p=',cn2m,cn2p
cyup             write(*,*)'the given wave can not exist in plasma'
	     iraystop=1
	     return
	    endif

20	    iroot=iroot+1
cyup            write(*,*)'cninit12_n_gam iroot,cntang2',iroot,cntang2

	    if(iroot.eq.1) then
              cn2=cn2p
	      if(cn2.lt.cntang2)then
cyup	        write(*,*)'cninit12_n_gam cn2p.lt.cntang2'
	        go to 20
	      end if
	    else
              cn2=cn2m
	      if(cn2.lt.cntang2)then
	        if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
	        write(*,*)'cninit12_n_gam cn2m.lt.cntang2'
	        write(*,*)'the given wave can not exist in plasma'
	        endif ! outprint
	        iraystop=1
	        return
	      end if
	    end if
	    
cyup            write(*,*)'cninit12_n_gam cn2=',cn2
	    cnrho2=cn2-cnteta**2-cnphi**2
	    cnrho=dsqrt(cnrho2)
            cnperp=cnrho
cSm060906
            cnperp=dsqrt(cn2-cnpar2)
cyup            write(*,*)'cnpar,dsqrt(cn2-cnpar2)',cnpar,dsqrt(cn2-cnpar2)
c-------------------------------------------------------------------
cyup          write(*,*)'cninit12_n_gam cn2',cn2
          call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)

cyup          write(*,*)'cninit12_n_gam  after cnzcnr z,r,phi', z,r,phi
cyup          write(*,*)'cnteta,cnphi,cnrho',cnteta,cnphi,cnrho
cyup          write(*,*)'cnz,cnr,cm',cnz,cnr,cm
c-------------------------------------------------------------------
          cm=cnphi*r
	  gam=gamma1(z,r,phi,cnz,cnr,cm)
          ds=dsin(gam)
          dc=dcos(gam)
          ds2=ds*ds ! sin^2(gam)
          dc2=dc*dc ! cos^2(gam)
          ds4=ds2*ds2 ! not used here
c---------------------------------------------------------------------
c         verify that cn2 and gam are the solution of the dispersion
c         relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c---------------------------------------------------------------------
          call abc(z,r,phi,ds2,dc2,ad,bd,cd) ! here: cninit12_n_gam
          
      ! Note: a=A*delta, b=B*delta, c=C*delta (where delta=1-Y)
      ! and sqrt(det)= sqrt(B*B-4*A*C)*|delta|
      !in    N^2 = (-b +ioxm*sqrt(b*b-4*a*c))/(2a)


          d4=ad
	  d2=bd
	  d0=cd
          det=d2*d2-4.d0*d4*d0
cyup          write(*,*)'--- cninit12_n_gam, for abc: gam===',gam
c-------------------------------------------------------------
          cn1=dsqrt(d2*d2-4.d0*d4*d0)
          cn2p_ib1=(-d2+sign_del*cn1)/(2.d0*d4)
          cn2m_ib1=(-d2-sign_del*cn1)/(2.d0*d4)
cyup          write(*,*)'from abc: (-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
cyup          write(*,*)'from abc: (-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)
cyup          write(*,*)'cninit12_n_gam ib=1 from abc: cn2p,cn2m',cn2p_ib1,cn2m_ib1
          
          cn2new=(-d2+ioxm*sign_del*cn1)/(2.d0*d4)
cyup       write(*,*)'cninit12_n_gam cn2new=',cn2new,'iroot=',iroot
	  if(iroot.eq.1) then
	      dnp=dabs(cn2-cn2new)
c              write(*,*)'cninit12 dnp,epsmode',dnp,epsmode
	      if (dnp.gt.epsmode)then
	         goto 20
	      end if
	  else
	      dnm=dabs(cn2-cn2new)
	      if (dnm.gt.epsmode)then
                write(*,*)'the given mode can not exist in plasma'
	        iraystop=1
	        return
	      end if
	  end if
          goto 111
        end if ! ib=1
c end if 2
c
c       ib.gt.1 iones resonance condition may be
c  if 3
        if (ib.gt.1) then
c         write(*,*)'in cninit: cold plasma ib .gt.1 x,r.phi,ib',
c     +   z,r,phi,ib
          ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
          xb=x(z,r,phi,ibmx)
          yb=y(z,r,phi,ibmx)
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)

c          write(*,*)'in cninit12_n_gam z,r,ib,xb,yb',z,r,ib,xb,yb
c          write(*,*)'in cninit12_n_gam cnpar',cnpar
c          write(*,*)'in cninit12_n_gam s1,s2,s3,s4',s1,s2,s3,s4

	  pype=xe/(1.d0+ye)
	  pyme=xe/(1.d0-ye)
	  pyme2=pype/(1.d0-ye)
	  pypb=xb/(1.d0+yb)
	  delib=1.d0-yb
	  f1b=(s1-pyme2)
	  f0b=-pypb

	  g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
	  g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

	  w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1            s4*(s2-pype)*(s3-pyme)
          w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
	  fd=f1b*delib+f0b ! see below
	  gd=g1b*delib+g0b ! see below
	  wd=w1b*delib+w0b ! see bolow
	  detin=gd**2-4.d0*fd*wd

c       write(*,*)'in cninit12_n_gam f1b,delib,f0b,fd',f1b,delib,f0b,fd
c       write(*,*)'in cninit12_n_gam g1b,delib,g0b,gd',g1b,delib,g0b,gd
c       write(*,*)'in cninit12_n_gam w1b,delib,w0b,wd',w1b,delib,w0b,wd

c          write(*,*)'in cninit12_n_gam gd,fd,wd,detin',gd,fd,wd,detin

	  if (detin.lt.0d0) then
	      write(*,*)' 3 cninit12_n_gam detin  less then zero '
	      iraystop=1
	      return
	  end if
c          write(*,*)'cninit: gd,fd,wd,detin',gd,fd,wd,detin
	  cn2p=(-gd+sign_del*dsqrt(detin))/(2.d0*fd)
	  cn2m=(-gd-sign_del*dsqrt(detin))/(2.d0*fd)
         
          if((cn2m.lt.0.d0).and.(cn2p.lt.0.d0)) then
            if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
            write(*,*)'cninit12_n_gam two roots of the dispersion < 0'
            write(*,*)'cn2m,cn2p',cn2m,cn2p
            write(*,*)'the given wave can not exist in plasma'
            endif ! outprint
	    iraystop=1
	    return
	   endif

cyup	  write(*,*)'cninit12_n_gam ib.qt.1 cn2p,cn2m',cn2p,cn2m
c	  write(*,*)'in cninit fd,gd,wd',fd,gd,wd
30	  iroot=iroot+1

cyup          write(*,*)'cninit12_n_gam iroot,cntang2',iroot,cntang2

	  if(iroot.eq.1) then
             cn2=cn2p
	     if(cn2.lt.cntang2)then
cyup	        write(*,*)'cninit12_n_gam cn2p.lt.cntang2'
	        go to 30
	     end if
	  else
             cn2=cn2m
	     if(cn2.lt.cntang2)then
	        if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
	        write(*,*)'cninit12_n_gam cn2m.lt.cntang2'
	        write(*,*)'the given wave can not exist in plasma'
	        endif ! outprint
	        iraystop=1
	        return
	     end if
	  end if

cyup          write(*,*)'cninit12_n_gam cn2=',cn2
	  cnrho2=cn2-cnteta**2-cnphi**2
	  cnrho=dsqrt(cnrho2)   
          cnperp=cnrho
cyup          write(*,*)'cnperp=cnrho',cnperp
cSm050826
          cnperp=dsqrt(cn2-cnpar2)
cyup          write(*,*)'cnpar,dsqrt(cn2-cnpar2)',cnpar,dsqrt(cn2-cnpar2)
c------------------------------------------------------------------
c          write(*,*)'cnini12 ib>1 before call cnzcnr'
c          write(*,*)'z,r,phi,cnteta,cnphi,cnrho',
c     &    z,r,phi,cnteta,cnphi,cnrho 
          call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)
c          write(*,*)'cninit12 after cnzcnr cnz,cnr,cm',cnz,cnr,cm

c---------test begin
cSm050906
cyup          cnpar_test=(bz*cnz+br*cnr+bphi*cm/r)/bmod
cyup          write(*,*)'cninit12_n_gam cnpar,cnpar_test', cnpar,cnpar_test
cyup          write(*,*)'cn2,cnz**2+cnr**2+(cm/r)**2',
cyup     &    cn2,cnz**2+cnr**2+(cm/r)**2
c---------test_end

c---------------------------------------------------------------------
          cm=cnphi*r
	  gam=gamma1(z,r,phi,cnz,cnr,cm)
          !write(*,*)'cninit12 from gamma1 gam',gam 

cSm050825
c---------test begin
c          gam=dacos(cnpar/dsqrt(cn2))
c---------test_end

          ds=dsin(gam)
          dc=dcos(gam)
          ds2=ds*ds
          dc2=dc*dc
          ds4=ds2*ds2 ! not used here
c--------------------------------------------------------------------
c test hamilt
          !write(*,*)'cninit12_n_gam cnpar2,cn2*dc2',cnpar2,cn2*dc2
c--------------------------------------------------------------------
c         verify that cn2 and gam are the solution of the dispersion
c         relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c---------------------------------------------------------------------
          call abc(z,r,phi,ds2,dc2,ad,bd,cd)
          d4=ad
	  d2=bd
	  d0=cd
          det=d2*d2-4.d0*d4*d0
          !write(*,*)'gam,d4,d2,d0,det',gam,d4,d2,d0,det
c---------------------------------------------------
c test hamilt=0?
          hamtest=fd*cn2*cn2+gd*cn2+wd

         !write(*,*)'cninit12_n_gam fd,gd,wd,cn2*cn2,cn2,hamtest'
         !write(*,*) fd,gd,wd,cn2*cn2,cn2,hamtest

          hamtest=d4*cn2*cn2+d2*cn2+d0
          pt4=d4*cn2*cn2
          pt2=d2*cn2
          pt=pt4+pt2+d0
          ptm=pt4-pt2+d0
c          write(*,*)'cninit12 pt4,pt2,d0',pt4,pt2,d0
c          write(*,*)'pt,ptm',pt,ptm
          !write(*,*)'cninit12_n_gam d4,d2,d0,cn4,cn2,hamtest'
          !write(*,*)d4,d2,d0,cn2*cn2,cn2,hamtest
c         write(*,*)'(fd-d4),(gd-d2),(wd-d0)'
c         write(*,*)(fd-d4),(gd-d2),(wd-d0)
          gdt=dc2*cn2*((1.d0-xe-xb)
     1       -(1.d0-xb/(1.d0-yb*yb)-xe/(1.d0-ye*ye)))
     1       +(-(1.d0-xb/(1.d0-yb)-xe/(1.d0+ye))*
     1       (1.d0-xb/(1.d0+yb)-xe/(1.d0-ye))-
     1       (1.d0-xb/(1.d0-yb*yb)-xe/(1.d0-xe*xe))*(1.d0-xe-xb))
          gdt=gdt*delib ! for test only
c         write(*,*)'gd,gdt',gd,gdt
          btest=-(1.d0-xb/(1.d0-yb*yb)-
     1         xe/(1.d0-ye*ye))*(1.d0-xe-xb)*
     1         (1.d0+dc2)-(1.d0-xb/(1.d0-yb)-xe/(1.d0+ye))*
     1         (1.d0-xb/(1.d0+yb)-xe/(1.d0-ye))*ds2
c         write(*,*)'bd,btest',bd,btest
          ptt=(xe*ye*ye*delib/(1.d0-ye*ye)+xb*yb*yb/(1.d0+yb)) ! for test only
          fmina=-dc2*ptt
c         write(*,*)'xe,ye,xb,yb',xe,ye,xb,yb
          ptt1=(-(1.d0-xb/(1.d0-yb)-xe/(1.d0+ye))*
     1         (1.d0-xb/(1.d0+yb)-xe/(1.d0-ye))+
     1    (1.d0-xe-xb)*(1.d0-xe/(1.d0-ye*ye)-xb/(1.d0-yb*yb)))*delib ! for test only
          gminb=dc2*(-cn2*ptt-ptt1)
          wminc=-dc2*cn2*ptt1
c         write(*,*)'fmina,gminb,wminc'
c         write(*,*)fmina,gminb,wminc
c-------------------------------------------------------------
          cn1=dsqrt(d2*d2-4.d0*d4*d0)
cyup          write(*,*)'(-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
cyup          write(*,*)'(-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)
          cn2new=(-d2+ioxm*sign_del*cn1)/(2.d0*d4)
cyup	  write(*,*)'cn2new=',cn2new
	  if(iroot.eq.1) then
	     dnp=dabs(cn2-cn2new)
	     if (dnp.gt.epsmode)then
	        goto 30
	     end if
	  else
	     dnm=dabs(cn2-cn2new)
	     if (dnm.gt.epsmode)then
	     if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
           write(*,*)'cninit12_n_gam given mode cannot exist in plasma'
           endif ! outprint
	        iraystop=1
	        return
	     end if
	  end if
          goto 111
        end if
c end if 3
      end if
c end if 0
 
  111 continue
cyup      write(*,*)'cninit12_n_gam FOUND: cn2=',cn2,cnz**2+cnr**2+(cm/r)**2
                                                              
      return
      end ! cninit12_n_gam



      subroutine n_cold_gam(z,r,phi,ib,gam,cn_p,cn_m,iraystop_p,
     + iraystop_m)
c----------------------------------------------------------------
c     Used by cninit12/ioxm_n_npar section to find a matching ioxm
c     that satisfies N(Npar,ioxm_n_npar) = N(gam,ioxm)
c     Calculates cold plasma dispersion relation roots
c     cn=N(gam,ioxm) in the given space point (z,r,phi).
c     as a root of the equation a*N**4 + b*N**2 + c =0,
c     cn=(-b+ioxm*sqrt(b**2-4ac))/2a
c     cn_p=(-b+sqrt(b**2-4ac))/2a for ioxm=+1
c     cn_m=(-b-sqrt(b**2-4ac))/2a for ioxm=-1
c
c     Here: 
c     coefficients a=A*delta, b=B*delta, c=C*delta,
c
c     See book  Krall,Trivelpiece
c     A= eps_1*sin(gam)**2+eps_3*cos(gam)**2
c     B=-eps_1*eps_2(1+cos(gam)**2)-(eps_1**2-eps_2**2)sin(gam)**2
c     C= eps_3*(eps_1**2-eps_2**2)
c
c     delta=1-y_e for electrons ib=1 (ib is used in subroutine
c                                     abc) 
c          =1-y_i for ions      ib=i>1 (ib=i means a resonance 
c                               with this ion species is expected)
c----------------------------------------------------------------     
      implicit none
      
c-----input
      real*8 z,r,phi, ! space coordinates of the given point
     *gam             ! the angle [radians] between the wave vector
                      ! and the magnetic field                      
      ! ib  YuP[2018-05-24] Added ib as an input, see above
      
c-----output
      real*8 cn_p,cn_m    !roots cn_p=N(gam,ioxm=+1), cn_m=N(gam,ioxm=-1),  
      integer iraystop_p, !=0 the root with ioxm=+1 found
                          !=1 no root  with ioxm=+1 did not find
     &        iraystop_m  !=0 the root with ioxm=-1 found
                          !=1 no root  with ioxm=-1 did not find
 
c-----locals
      real*8 ds,dc,ds2,dc2,
     &a,b,c,det,cn2_p,cn2_m
     
      real*8 delib, sign_del, y
      integer ibmx !local
      integer ib

      iraystop_p=0
      iraystop_m=0
      
      ds=dsin(gam)
      dc=dcos(gam)
      ds2=ds*ds
      dc2=dc*dc

      call abc(z,r,phi,ds2,dc2,a,b,c)

      det=b**2-4.d0*a*c
      if(det.lt.0.d0)then
        write(*,*)'in subroutine n_cold_gam det<0'
        write(*,*)'the cold plasma dispersion eq.'
        write(*,*)'has no root N_perp(gam)'
        write(*,*)'for the given angle gam=',gam
        iraystop_p=1
        iraystop_m=1
        return
      endif

      ! Note: a=A*delta, b=B*delta, c=C*delta (where delta=1-Y)
      ! and sqrt(det)= sqrt(B*B-4*A*C) * |delta|
      ibmx=ib !min(ib,nbulk) ! safety check: not to exceed nbulk
      delib= 1.d0-y(z,r,phi,ibmx) ! 1-Y (can be electrons or ions)
      sign_del=1.d0 !sign(1.d0,delib) ! in n_cold_gam
      
      cn2_p=(-b+sign_del*dsqrt(det))/(2.d0*a) !N**2, Corresp to ioxm=+1
      if(cn2_p.gt.0.d0)then
        cn_p=dsqrt(cn2_p) 
      else
        write(*,*)'in subr. n_cold_gam: Negative N**2<0 (cn2_p<0)'
        write(*,*)'The cold plasma dispersion eq.'
        write(*,*)'has no positive root N^2(gam,ioxm)'
        write(*,*)'for the given angle gam=',gam,' and ioxm=+1'
        iraystop_p=1
      endif

      cn2_m=(-b-sign_del*dsqrt(det))/(2.d0*a) !N**2  Corresp to ioxm=-1
      if(cn2_m.gt.0.d0)then
        cn_m=dsqrt(cn2_m) 
      else
        write(*,*)'in subr. n_cold_gam: Negative N**2<0 (cn2_m<0)'
        write(*,*)'The cold plasma dispersion eq.'
        write(*,*)'has not positive root N^2(gam,ioxm)'
        write(*,*)'for the given angle gam=',gam,' and ioxm=-1'
        iraystop_p=1
      endif

      return
      end !  n_cold_gam

      

c        **********************nper_npar_ioxm_n_npar**********
c                                -                            
c         It solves the cold plasma id=1,2 or
c         Appleton-Harty id=3 dispersion relation N_per=N_per(N_par)
c         for given ioxm_n_npar                              
c         Then subroutine calculates the initial components  
c         of the refractive index  cnz,cnr,cm                
c        ******************************************************
c
c---------------------------------------------------
c     							   !
c     input parameters					   !
c     z,r,phi,cnpar 	                                   !
c     ioxm_n_npar is inside one.i common block             !
c     It works for cold plasma dispesion functions         !
c     id_loc=1,2,3                                         !
c     output parameters:N_per(N_par)=cnper                 !
c     iraystop=0 the root was found                        !
c     iraystop=1 root was not found                        !
c----------------------------------------------------------
c     it uses the following functions and subroutines      !
c     ias1r,b,y,x,gamma1,s                            !
c-----------------------------------------------------------

      subroutine nper_npar_ioxm_n_npar(id_loc,z,r,phi,cnpar,
     &cnper,iraystop)

      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      iraystop=0

      if(((id_loc.ne.1).or.(id_loc.ne.2)).or.(id_loc.ne.3))then
c        write(*,*)'nper_npar_ioxm_n_npar id_loc.ne.1, 2 or 3'
c        iraystop=1
c        return
      endif

      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar ! from INPUT
      cnpar4=cnpar2*cnpar2
      bmod=b(z,r,phi)
c-------------------------------------------------------------
c     calculations of  cnrer2p and cnper2m
c     from cnpar2 by using the dispersion relation
c     f*cnper**4 +(2f*cnpar**2+g)*cnper**2+(f*cnpar**4+g*cnpar**2+w)=0
c     f*cnper**4 +gnew*cnper**2+wnew=0
c     gnew=2fcnpar**2+g, wnew=f*cnpar**4+g*cnpar**2+w
c     in the following form cnper**2=(-gnew+,-*sqrt(gnew**2-4*f*wnew))/(2*f)
      ! Note: a=A*delta, b=B*delta, c=C*delta (where delta=1-Y)
      ! and sqrt(det)= sqrt(B*B-4*A*C) * |delta|
      ! in   N^2 = (-b +ioxm*sqrt(b*b-4*a*c))/(2a)
      ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
      delib= 1.d0-y(z,r,phi,ibmx) ! 1-Y (electrons or ions)
      sign_del=1.d0 !sign(1.d0,delib) ! in nper_npar_ioxm_n_npar
c------------------------------------------------------------------
c     Appleton - Hartry dispersion relation
      if (id_loc.eq.3) then ! id.eq.3
         xi=x(z,r,phi,1)
         yi=y(z,r,phi,1)
         py2=yi*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
c------------------------------------------------------------------
         pyp=xi/(1.d0+yi)
         f1e=1.d0
         f0e=-pyp
         g1e=cnpar2*(-xi)+pyp+xi-2.d0
         g0e=cnpar2*pyp+xi*(1.d0-pyp)+pyp*(1.d0-xi)
         w1e=cnpar2*(xi-pyp)+(1.d0-pyp)*(1.d0-xi)
         w0e=cnpar2*(pyp*(1.d0-xi)-xi*(1.d0-pyp))-
     &       (1.d0-xi)*(1.d0-pyp)*xi
         dele=1.-yi
         fd=f1e*dele+f0e
         gd=g1e*dele+g0e
         wd=w1e*dele+w0e  !==(1-Y)*w
c new coefficients
         gnew=gd+2.d0*fd*cnpar2
         wnew=wd+gd*cnpar2+fd*cnpar4
         gd=gnew
         wd=wnew

         detin=gd**2-4.d0*fd*wd
         if(detin.lt.0.d0) then
           write(*,*)'1. nper_npar_ioxm_n_nparn detin less then zero'
           write(*,*)'nper_npar_ioxm_n_npar no roots'
           iraystop=1
           return
         endif

         cnper2=(-gd+ioxm_n_npar*sign_del*dsqrt(detin))/(2.d0*fd)
         if(cnper2.lt.0d0) then
            write(*,*)'2. nper_npar_ioxm_n_nparn cnper2<0'
            write(*,*)'N_perp^2 is negative for the given ioxm_n_npar=',
     +                 ioxm_n_npar
            write(*,*)'nper_npar_ioxm_n_npar no positive root.'
            write(*,*)'FYI: the root for opposite ioxm_n_npar is',
     +                 (-gd -ioxm_n_npar*dsqrt(detin))/(2.d0*fd)
            iraystop=1
            return
         else
            cnper=dsqrt(cnper2)
         endif           
      end if !id_loc.eq.3
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
      if ((id_loc.eq.1).or.(id_loc.eq.2)) then  ! id.eq.1 or id.eq.2
        call s(z,r,phi,s1,s2,s3,s4,s6,s7)

        if (ib.eq.1) then
c---------ib=1 electron resonance condition may be
c         write(*,*)'nper_npar_ioxm_n_npar cold plasma ib=1 '
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)
          pyp=xe/(1.d0+ye)
          dele=1.d0-ye
          f1e=s7
          f0e=-pyp
          g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
          g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4
          w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4
          w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe
          fd=f1e*dele+f0e !== (1-Ye)*f == (1-Ye)*eps1 
          gd=g1e*dele+g0e !== (1-Ye)*g
          !gd= (1-Ye)*[Npar^2 *(eps3-eps1) +eps2^2 -eps1^2 -eps1*eps3]
          wd=w1e*dele+w0e !== (1-Ye)*w
          !wd= (1-Ye)*[Npar^2 *(eps1^2-eps2^2-eps1*eps3) +eps3*(eps1^2-eps2^2)]
c new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew

          detin=gd**2-4.d0*fd*wd
          if(detin.lt.0.d0) then           
            write(*,*)'2 nper_npar_ioxm_n_npar detin less then zero'
            write(*,*)'nper_npar_ioxm_n_npar no roots'
            iraystop=1
            return
          end if

          
          cnper2=(-gd+ioxm_n_npar*sign_del*dsqrt(detin))/(2.d0*fd) !Nper^2(Npar)
          if(cnper2.lt.0.d0)then
            write(*,*)'2, nper_npar_ioxm_n_nparn cnper2<0'
            write(*,*)'N_perp^2 is negative for the given ioxm_n_npar=',
     +                 ioxm_n_npar
            write(*,*)'nper_npar_ioxm_n_npar no positive root.'
            write(*,*)'FYI: the root for opposite ioxm_n_npar is',
     +                 (-gd -ioxm_n_npar*dsqrt(detin))/(2.d0*fd)
            iraystop=1
            return
          else
            cnper=dsqrt(cnper2)
          endif           
          goto 111
        end if !ib.eq.1


        if(ib.gt.1) then
c---------ib.gt.1 iones resonance condition may be
c         write(*,*)'cold plasma ib .gt.1  '
          ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
          xb=x(z,r,phi,ibmx)
          yb=y(z,r,phi,ibmx)
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)
          pype=xe/(1.d0+ye)
          pyme=xe/(1.d0-ye)
          pyme2=pype/(1.d0-ye)
          pypb=xb/(1.d0+yb)
          delib=1.d0-yb
	  f1b=(s1-pyme2)
	  f0b=-pypb

	  g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
          g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

	  w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1	     s4*(s2-pype)*(s3-pyme)
          w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
          fd=f1b*delib+f0b ! see below
          gd=g1b*delib+g0b ! see below
          wd=w1b*delib+w0b ! see below
c new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew

          detin=gd**2-4.d0*fd*wd
          if(detin.lt.0.d0) then           
            write(*,*)'2 nper_npar_ioxm_n_npar detin less then zero'
            write(*,*)'nper_npar_ioxm_n_npar no roots'
            iraystop=1
            return
	  end if

          cnper2=(-gd+ioxm_n_npar*sign_del*dsqrt(detin))/(2.d0*fd) !Nper^2(Npar)
          if(cnper2.lt.0.d0) then
            write(*,*)'3. nper_npar_ioxm_n_nparn cnper2<0'
            write(*,*)'N_perp^2 is negative for the given ioxm_n_npar=',
     +                 ioxm_n_npar
            write(*,*)'nper_npar_ioxm_n_npar no positive root.'
            write(*,*)'FYI: the root for opposite ioxm_n_npar is',
     +                 (-gd -ioxm_n_npar*dsqrt(detin))/(2.d0*fd)
            iraystop=1
            return
          else
            cnper=dsqrt(cnper2)
          endif           
          goto 111
        end if !ib>1
      end if !(id_loc.eq.1).or.(id_loc.eq.2)
 
  111 continue

      return
      end ! nper_npar_ioxm_n_npar


c        **********************cninit12***************************
c        *                        -                               *
c        * It solves the dispersion relation N=N(n_par)           *
c        * for cold plasma                                        *
c        * Then subroutine calculates the initial components      *
c        * of the refractive index  cnz,cnr,cm         	          *
c        * If ioxm_n_npar=0 it uses root N(gam,ioxm)              *
c        * If ioxm_n_npar=1 or -1 it uses root N(N_par,ioxm_n_npar)*
c        *                then calculates angle gam(N_par,N_per)  *
c        *                then finds ioxm which gives the root    *
c        *                N(gam,ioxm)=N(N_par,ioxm_n_npar)        *
c        
c        *********************************************************
c
c------------------------------------------------------------------
c					                          !
c     input parameters		         		          !
c     z,r,phi,cnpar,cnteta,cnphi 		       	          !
c                                                                 !
c     output parameters				                  !
c     cnz,cnr,cm are the components of the refractive index       !        
c     iraystop=1 end ray calculation                              !
c------------------------------------------------------------------
c     it uses the following functions                             !
c     b,gamma1                                                    !
c------------------------------------------------------------------
      subroutine cninit12(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
cSAP090504
      include 'grill.i'

      save
      double complex cmplnper      
      
      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnteta*cnteta+cnphi*cnphi
      bmod=b(z,r,phi)
cyup      write(*,*)'cninit12 ioxm_n_npar=',ioxm_n_npar
      
      if(ioxm_n_npar.eq.0) then ! ioxm_n_npar=0  Use ioxm instead
c-------------------------------------------------------------
c       calculates the root N(Npar)=N(gam,ioxm) for givem ioxm
c-------------------------------------------------------------
        call cninit12_n_gam(z,r,phi,cnpar,cnteta,cnphi,
     &                      cnz,cnr,cm,iraystop)
      else ! ioxm_n_npar.eq.1 or -1
c-------------------------------------------------------------------
c       calculates the root N(N_par,ioxm_n_npar) for given ioxm_n_npar,
c       then finds matching ioxm to get  N(N_par,ioxm_n_npar)=N(gam,ioxm)
c-------------------------------------------------------------------
cyup        write(*,*)'cninit12/ioxm_n_npar: z,r,phi,cnpar,cnteta,cnphi',
cyup     &             z,r,phi,cnpar,cnteta,cnphi
c-----------------------------------------------------------------
        iroot=0 ! will be 1, or 2(if two matching ioxm are found)
c------------------------------------------------------------------
c       epsmode is the accuracy  for selection mode (ioxm) on
c             the plasma boundary
c------------------------------------------------------------------
        epsmode=1.d-8

        id_loc=2
c------------------------------------------------------------------
c       solve the cold plasma id=1,2  
c       dispersion relation cnpar=N_per(n_par)
c------------------------------------------------------------------        
        call nper_npar_ioxm_n_npar(id_loc,z,r,phi,cnpar,
     &  cnper,iraystop) !ioxm_n_npar was set in one.i
cyup        write(*,*)'cninit12 after nper_npar_ioxm_n_npar'
cyup        write(*,*)'ioxm_n_npar,cnper,iraystop',
cyup     &             ioxm_n_npar,cnper,iraystop
c-------------------------------------------------------------------
        if(iraystop.eq.1) then
          write(*,*)'cninit12/ioxm_n_npar: no root N(Npar,ioxm_n_npar)'
          return
        else 
cSAP090504  
           if (i_n_poloidal.eq.3) then !input N_parallel, ksi_nperp
c-------------calculate N_phi,N_theta,N_rho

              gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
              b_teta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field

c              write(*,*)'cninit.f i_n_poloidal=3 ksi_nperp',ksi_nperp

              rad_ksi_nperp=ksi_nperp*pi/180.d0 !transfrm degrees to radians

c              write(*,*)'rad_ksi_nperp',rad_ksi_nperp

              cnteta=(cnpar*b_teta+cnper*bphi*dsin(rad_ksi_nperp))
     &                    /bmod
              cnphi=(cnpar*bphi-cnper*b_teta*dsin(rad_ksi_nperp))
     &                    /bmod   

cSAP091026 It was that for ksi=0 the refractive vector
c          was directed opposite to grad(psi)
c              cnrho=cnper*dcos(rad_ksi_nperp)
              cnrho=-cnper*dcos(rad_ksi_nperp)
              if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
              write(*,*)'cninit12/i_n_poloidal=3 cnteta,cnphi,cnrho',
     &                                           cnteta,cnphi,cnrho
              write(*,*)'cninit12/ioxm_n_npar:  i_n_poloidal=3', 
     &        'rad_ksi_nperp,dcos(rad_ksi_nperp),cnper',
     &        rad_ksi_nperp,dcos(rad_ksi_nperp),cnper
              endif ! outprint

c              write(*,*)'cnteta**2+cnrho**2',cnteta**2+cnrho**2

c              write(*,*)'cninit.f i_n_poloidal=3 cn2=cnpar**2+cnper**2',
c     &                                           cnpar**2+cnper**2

c              write(*,*)'cn2=cnteta**2+cnphi**2+cnrho**2',
c     &                       cnteta**2+cnphi**2+cnrho**2

c              write(*,*)'ksi_nperp,rad_ksi_nperp*180.d0/pi',
c     &                   ksi_nperp,rad_ksi_nperp*180.d0/pi

              cntang2=cnteta*cnteta+cnphi*cnphi              
              
           endif ! i_n_pol=3

          
          cn2_npar=cnper**2+cnpar**2
cyup          write(*,*)'cninit12/ioxm_n_npar: cn2_npar,cntang2',
cyup     +     cn2_npar,cntang2

          if(cn2_npar.lt.cntang2)then
            if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
            write(*,*)'cninit12/ioxm_n_npar: (cn2_npar.lt.cntang2)'
            write(*,*)'   no root N(Npar,ioxm_n_npar) > N_tang'
            endif ! outprint
            iraystop=1
            return
          else
cSAP091026  At (i_n_poloidal.eq.3) cnrho was calculated before this line, 
c           It can be negative (for i_n_poloidal.eq.3)
c           if the refractive index is directed otside the plasma
            if (i_n_poloidal.ne.3) then
               cnrho2=cn2_npar-cnteta**2-cnphi**2
	       cnrho=dsqrt(cnrho2)
            endif 
c--------------------------------------------------------------
c           calculate refractive index components cnz,cnr,cm
c           for given z,r,phi,cntheta,cnphi,cnrho
c-------------------------------------------------------------- 
cyup            write(*,*)'cninit12/ioxm_n_npar: before cnzcnr'
cyup            write(*,*)'r,phi,cnteta,cnphi,cnrho',
cyup     &                 r,phi,cnteta,cnphi,cnrho

cSAP090601
c           In this case ioxm_n_npar.ne.0 and
c           ioxm is not determined
c           Subroutine cnzcnr uses subroutine rside1, 
c           which for id=2 case needs the determined value ioxm 
c           To avoid this problem we will use id=1 in cnzcnr
c           May be signs of the radial group velocity
c           for id=1 and id=2 can be different?
c           If yes, then it can creat the new problem. 

            if (id.eq.2) then
              id_loc_2=id
              id=1
            endif 

            call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)
cSAP090601          
            id=id_loc_2 

cyup            write(*,*)'cninit12/ioxm_n_npar:  after cnzcnr'
cyup            write(*,*)'cnz,cnr,cm',cnz,cnr,cm

ccSAP091024------------------------------------------------test-------
c            if (i_n_poloidal.eq.3) then !input N_parallel, ksi_nperp
cc--------------test ksi_nperp
cc               write(*,*)'cnteta**2+cnrho**2',cnteta**2+cnrho**2
cc               write(*,*)'cnz**2+cnr**2',cnz**2+cnr**2

c               bmod=b(z,r,phi)

c               cnpar=(cnz*bz+cnr*br+(cm/r)*bphi)/bmod

c               cnpar_z=cnpar*bz/bmod
c               cnpar_r=cnpar*br/bmod
c               cnpar_phi=cnpar*bphi/bmod

cc               write(*,*)'cnpar,cnpar_z,cnpar_r,cnpar_phi',
cc     &                    cnpar,cnpar_z,cnpar_r,cnpar_phi

c               cnper_z=cnz-cnpar_z
c               cnper_r=cnr-cnpar_r
c               cnper_phi=(cm/r)-cnpar_phi

cc               write(*,*)'cnper_phi,cnpar_phi,cm/r',
cc     &                    cnper_phi,cnpar_phi,cm/r
             
cc               write(*,*)'dsqrt(cnper_z**2+cnper_r**2+cnper_phi**2)',
cc     &                    dsqrt(cnper_z**2+cnper_r**2+cnper_phi**2)
cc               write(*,*)'cnrho',cnrho,'cnper',cnper

cc               write(*,*)'cnper_z,cnper_r,cnper_phi',
cc     &                    cnper_z,cnper_r,cnper_phi

cc               write(*,*)'dpdzd,dpdrd,dsqrt(dpdrd*2+dpdzd**2)',
cc     &         dpdzd,dpdrd,dsqrt(dpdrd*2+dpdzd**2)

c               cnper_test=dsqrt(cnper_r**2+cnper_z**2+cnper_phi**2)
 
cc              write(*,*)'cnpar,cnper,cnper_test',cnpar,cnper,cnper_test

c               write(*,*)' cninit================='
cc               write(*,*)'cnper_z,cnper_r',cnper_z,cnper_r
cc               write(*,*)'dpdzd,dpdrd',dpdzd,dpdrd
cc               write(*,*)'dsqrt(dpdzd**2+dpdrd**2)',
cc     &                    dsqrt(dpdzd**2+dpdrd**2)
cc               write(*,*)'cnper_test',cnper_test

c               cos_ksi_test=(cnper_z*dpdzd+cnper_r*dpdrd)/
c     &                      (dsqrt(dpdzd**2+dpdrd**2)*cnper_test)
    

c               write(*,*)'cos_ksi_test,dcos(rad_ksi_nperp)',
c     &                    cos_ksi_test,dcos(rad_ksi_nperp)
c            endif !test
c--------------------------------------------------------------
c           calculate the angle gam between the refractive index 
c           and the magnetic field
c-------------------------------------------------------------
            gam=gamma1(z,r,phi,cnz,cnr,cm)
c-------------------------------------------------------------
c           calculate cold plasma roots roots 
c           cn_p=N(gam,ioxm=+1) and cn_m=N(gam,ioxm=-1)
c           Trying to find a matching ioxm 
c             for a found solution N(Npar,ioxm_n_npar)
c-------------------------------------------------------------  
            ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
            call n_cold_gam(z,r,phi,ibmx,gam,cn_p,cn_m,
     &                      iraystop_p,iraystop_m)
c-------------------------------------------------------------
c           choose ioxm for which
c           N(N_par,ioxm_n_npar)=N(gam,ioxm)
c-------------------------------------------------------------
c           check ioxm=+1 root
c-------------------------------------------------------------
            if(iraystop_p.eq.0) then
              delta=dabs(cn2_npar-cn_p**2)
              if(delta.gt.epsmode)then
c---------------root N(N_par,ioxm_n_npar).ne.N(gam,ioxm=+1)
              else
c---------------root N(N_par,ioxm_n_npar)=N(gam,ioxm=+1)
                ioxm=+1
                iroot=iroot+1
                write(*,*)'cninit12/ioxm_n_npar: found ioxm=',ioxm 
c               goto 10 
              endif                
            else
c-------------no root for N(gam,ioxm=+1)
            endif !iraystop_p.eq.0

c-------------------------------------------------------------
c           check ioxm=-1 root
c-------------------------------------------------------------
            if(iraystop_m.eq.0) then
              delta=dabs(cn2_npar-cn_m**2)
              if(delta.gt.epsmode)then
c---------------root N(N_par,ioxm_n_npar).ne.N(gam,ioxm=-1)
              else
c---------------root N(N_par,ioxm_n_npar)=N(gam,ioxm=-1)
                ioxm=-1                  
                write(*,*)'cninit12/ioxm_n_npar: found ioxm=',ioxm 
                iroot=iroot+1
c               goto 10 
              endif                
            else
c-------------no root for N(gam,ioxm=-1)
            endif !iraystop_m.eq.0

c------------------------------------------------------------
            if(iroot.eq.0) then
              write(*,*)'cninit12/ioxm_n_npar: no matching ioxm '
              write(*,*)'to get N(N_par,ioxm_n_npar)=N(gam,ioxm)'
              iraystop=1
              return 
            else
              if(iroot.eq.2) then
                write(*,*)'*******WARNING******************'
                write(*,*)'cninit12/ioxm_n_npar: two ioxm was found'
                write(*,*)'to get N(N_par,ioxm_n_npar)=N(gam,ioxm)'
              endif !iroot.eq.2
            endif !iroot.eq.0
          endif !cn2_npar.lt.cntang2
        endif !iraystop.eq.1
        
      endif !ioxm_n_npar=0        ioxm_n_npar.eq.0

c  10 continue
cyup      write(*,*)'end of subroutine cninit12 iraystop',iraystop
      return
      end ! cninit12


c        **********************npernpar_test******************
c        *                        -                           *
c        * It solves the dispersion relation Nper=N(n_par)    *
c        * Then subroutine calculates the initial components  *
c        * of the refractive index  cnz,cnr,cm                *
c        ******************************************************
c
c---------------------------------------------------
c     							   !
c     input parameters					   !
c     z,r,phi,cnpar 	                                   !
c     output parameters:cnper2p,cnper2m	                   !
c---------------------------------------------------
c     it uses the following functions and subroutines      !
c     ias1r,b   ,y,x,gamma1,s,abc,hamilt1,                 !
c-----------------------------------------------------------
      subroutine npernpar_test(z,r,phi,cnpar,cnper2p,cnper2m)
c      implicit double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
c-----input
      real*8 z,r,phi,cnpar
c-----output
      real*8 cnper2p,cnper2m

c-----locals
      real*8 cnper,cnpar2,cnpar4,xi,yi,py2,py4,px,px2,pyp,
     &f1e,f0e,g1e,g0e,w1e,w0e,dele,fd,gd,wd,gnew,wnew,
     &detin,delib,
     &s1,s2,s3,s4,s6,s7,xe,ye,xb,yb,
     &pype,pyme,pyme2,pypb,
     &f1b,f0b,g1b,g0b,w1b,w0b,
     &cnprim
      real*8 sign_del
      integer ibmx !local

      integer
     &ioxmold,ioptmaz
c-----externals
      real*8 b,x,y

      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cnpar4=cnpar2*cnpar2
      bmod=b(z,r,phi)
c-------------------------------------------------------------
c     calculations of  cnrer2p and cnper2m
c     from cnpar2 by using the dispersin relation
c     f*cnper**4 +(2f*cnpar**2+g)*cnper**2+(f*cnpar**4+g*cnpar**2+w)=0
c     f*cnper**4 +gnew*cnper**2+wnew=0
c     gnew=2fcnpar**2+g, wnew=f*cnpar**4+g*cnpar**2+w
c     in the following form cnper**2=(-gnew+,-*sqrt(gnew**2-4*f*wnew))/(2*f)
      ! Note: a=A*delta, b=B*delta, c=C*delta (where delta=1-Y)
      ! and sqrt(det)= sqrt(B*B-4*A*C) * |delta|
      ! in   N^2 = (-b +ioxm*sqrt(b*b-4*a*c))/(2a)
      ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
      delib= 1.d0-y(z,r,phi,ibmx) ! 1-Y (electrons or ions)
      sign_del=1.d0 !sign(1.d0,delib) ! in npernpar_test
c------------------------------------------------------------------
c     Appleton - Hartry dispersion relation
c
c if 1
      if (id.eq.3) then
         xi=x(z,r,phi,1)
         yi=y(z,r,phi,1)
         py2=yi*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
c------------------------------------------------------------------
         pyp=xi/(1.d0+yi)
         f1e=1.d0
         f0e=-pyp
         g1e=cnpar2*(-xi)+pyp+xi-2.
         g0e=cnpar2*pyp+xi*(1.d0-pyp)+pyp*(1.d0-xi)
         w1e=cnpar2*(xi-pyp)+(1.-pyp)*(1.-xi)
         w0e=cnpar2*(pyp*(1-xi)-xi*(1.-pyp))-(1.d0-xi)*(1.d0-pyp)*xi
         dele=1.d0-yi
         fd=f1e*dele+f0e
         gd=g1e*dele+g0e
         wd=w1e*dele+w0e  !==(1-Y)*w
c new coefficients
         gnew=gd+2.d0*fd*cnpar2
         wnew=wd+gd*cnpar2+fd*cnpar4
         gd=gnew
         wd=wnew

         detin=gd**2-4.*fd*wd
         if (detin.lt.0d0) then
            write(*,*)' 1 in npernpar detin  less then zero '
            return
         endif
         cnper2p=(-gd+sign_del*dsqrt(detin))/(2.d0*fd)
         cnper2m=(-gd-sign_del*dsqrt(detin))/(2.d0*fd)
c         WRITE(*,*)'Aplt cnpernpar cnper2p,cnper2m',cnper2p,cnper2m
      end if
c end if 1
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
c if 0
      if ((id.eq.1).or.(id.eq.2)) then
        call s(z,r,phi,s1,s2,s3,s4,s6,s7)

c       ib=1 electron resonance condition may be
c  if 2
        if (ib.eq.1) then
c         write(*,*)'cnint cold plasma ib=1 '
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)
          pyp=xe/(1.d0+ye)
          dele=1.d0-ye
          f1e=s7
          f0e=-pyp
          g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
          g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4
          w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4
          w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe
          fd=f1e*dele+f0e
          gd=g1e*dele+g0e
          wd=w1e*dele+w0e  !==(1-Y)*w
c new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew

          detin=gd**2-4.d0*fd*wd
          if (detin.lt.0d0) then
             write(*,*)' 2 in npenpar detin  less then zero'
             cnper2p=-1.d0
             cnper2m=-1.d0
             return
	  end if
          cnper2p=(-gd+sign_del*dsqrt(detin))/(2.d0*fd)
          cnper2m=(-gd-sign_del*dsqrt(detin))/(2.d0*fd)
c          write(*,*)'in cninit.f  1 npernpar cnper2p,cnper2m',
c     .    cnper2p,cnper2m
          goto 111
        end if
c end if 2
c
c     ib.gt.1 ions resonance condition may be
c  if 3
        if (ib.gt.1) then
c          write(*,*)'cold plasma ib .gt.1  '
           ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
           xb=x(z,r,phi,ibmx)
           yb=y(z,r,phi,ibmx)
           xe=x(z,r,phi,1)
           ye=y(z,r,phi,1)
           pype=xe/(1.d0+ye)
           pyme=xe/(1.d0-ye)
           pyme2=pype/(1.d0-ye)
           pypb=xb/(1.d0+yb)
           delib=1.d0-yb
	   f1b=(s1-pyme2)
	   f0b=-pypb

	   g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
           g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

	   w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1	     s4*(s2-pype)*(s3-pyme)
           w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
           fd=f1b*delib+f0b
           gd=g1b*delib+g0b
           wd=w1b*delib+w0b
c new coefficients
           gnew=gd+2.d0*fd*cnpar2
           wnew=wd+gd*cnpar2+fd*cnpar4
           gd=gnew
           wd=wnew

           detin=gd**2-4.d0*fd*wd
           if (detin.lt.0d0) then
              write(*,*)' 3 in dinit detin  less then zero '
              return
           end if
           cnper2p=(-gd+sign_del*dsqrt(detin))/(2.d0*fd)
           cnper2m=(-gd-sign_del*dsqrt(detin))/(2.d0*fd)
c        write(*,*)'in npernpar ib.qt.1 cnper2p,cnper2m',cnper2p,cnper2m
           goto 111
        end if
c end if 3
      end if
c end if 0
c if 4
      if ((id.eq.4).or.(id.eq.5)) then
c-----------------------------------------------------------------
c        for cnper2 from cnpar by using dispersin relation
c        which was used in "Mazzucato's" code
c        Phys.Fluids 30(12),December 1987 p.3745
c-----------------------------------------------------------------     if ((id.eq.4).or.(id.eq.5)) then
c-----------------------------------------------------------------
c        for cnper2 from cnpar by using dispersin relation
c        which was used in "Mazzucato's" code
c        Phys.Fluids 30(12),December 1987 p.3745
c-----------------------------------------------------------------
         ioxmold=ioxm
         ioxm=1
         ioptmaz=1 !estimation of cnper from the cold plasma
         call cnpermuz(cnpar,iherm,z,r,phi,cnper,cnprim,ioptmaz)
         cnper2p=cnper*cnper
         ioxm=-1
         ioptmaz=1 !estimation of  cnper from the cold plasma
         call cnpermuz(cnpar,iherm,z,r,phi,cnper,cnprim,ioptmaz)
         cnper2m=cnper*cnper
c         write(*,*)'in cninit.f npernpar cnper2p,cnper2m',
c     +   cnper2p,cnper2m
         ioxm=ioxmold
         go to 111
      end if
c end if 4
  111 continue
      write(*,*)'in npernpar ib.qt.1 cnper2p,cnper2m',cnper2p,cnper2m
      return
      end ! npernpar_test
