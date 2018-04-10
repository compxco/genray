c
c        *********************dense_no_RZ_spline***************
c        *                        -                           *
c        * this function calculates the density profile       *
c        * as a function  from psi(z,r,phi) 		      *
c        *  rho=psi-psimag 				      *
c        *  a=psilim-psimag inside the function b()           *
c        *  these variables are in common 'one'
c        * This function has not chamber wall effects.
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c       rho is from common/one/                                    !
c------------------------------------------------------------------
c      uses
c         constants dense0,denseb,rn1de,rn2de,idens are in common 'one'
c         they were read in subroutine dinit from the file genray.in
c         rho is  in common 'one'
c      functions
c         _densrho_(rho)    finds dense from spline
c------------------------------------------------------------------
      real*8 FUNCTION dense_no_RZ_spline(z,r,phi,i)
c      IMPLICIT double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 z,r,phi!space coordinates
c     rho is in common  /one/
      integer i !plasma species number
c-----local
      real*8  theta_pol, !0=< theta_pol <2pi
     &rho_small,dens_rho_theta,
     &d_dens_rho_theta_d_rho,d_dens_rho_theta_d_theta
c-----externals
      real*8 thetapol,   ! -pi <thetapol =<pi
     &vardens,densrho
c      write(*,*)'in dense z,r,phi,i,nbulk',z,r,phi,i,nbulk
      if(i.gt.nbulk)then
        write(*,*)'in dense i.gt.nbulk: i,nbulk',i,nbulk
	stop
      endif
c      write(*,*)'in dense rho',rho
      if (rho.gt.1.d0-1.d-10) then
cSAP090206
c         dense=densrho(rho,i)

         theta_pol=thetapol(z,r) ! -pi <thetapol =<pi
         if (theta_pol.lt.0d0) then
            theta_pol=theta_pol+2*pi !pi< theta_pol<2pi
         endif

c--------calculate density outside LCFS
c         write(*,*)'dense.f in function dense before'
c         write(*,*)'dens_rho_theta_LCFS i',i
         rho_small=rho

c         write(*,*)'dense.f z,r,phi,i,rho_small',z,r,phi,i,rho_small

         call dens_rho_theta_LCFS(rho_small,theta_pol,i,
     &        dens_rho_theta,d_dens_rho_theta_d_rho,
     &        d_dens_rho_theta_d_theta)
c         write(*,*)'dense.f in function dense after'
c         write(*,*)'dens_rho_theta_LCFS i',i
          dense_no_RZ_spline = dens_rho_theta 

cSAP090222
c         write(*,*)'in dense i,rho,theta_pol,dense',
c     &                       i,rho,theta_pol,dense

      else
c      write(*,*)'in dense0 dense',dense
ctest beg
c         write(*,*)'in dense rho,i',rho,i
c         denst=densrho(rho,i)
c         vardent=vardens(z,r,phi)
c	 write(*,*)'in dense denst,vardent',denst,vardent
c tes end
          dense_no_RZ_spline=densrho(rho,i)*(1.0d0+vardens(z,r,phi))
c         write(*,*)'in dense0 dense',dense
      endif
c         write(*,*)'in dense0 dense',dense
      return
      END



c        *********************tempe ***************************
c        *                        -                           * 
c        * this function calculates temperature        
c        * as a function  from psi(z,r,phi) 		       	
c        *  a=psilim-psimag inside the function b()           
c        *  these variables are in common 'one'
c        *
c        * For n_wall=0
c        * temperho=dmax1(temperho1,temp_min_edge)
c        * Here
c        * temperho1=tempedg*dexp(-(rho_small-1.d0)/sigmedgt)     
c        *
c        * For n_wall > 0
c        * this function uses 2D spline at RZ mesh outside LCFS
c        * which has temperture fall otutside LCFS.      
c        
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c       rho is from common/one/
c------------------------------------------------------------------
      real*8 function tempe(z,r,phi,i)
      implicit none

      include 'param.i'
      include 'one.i'
      include 'edge_prof_nml.i'

c-----input
      real*8 z,r,phi !sapce coordinates
      integer i     !number of plasma specie

c-----externals
      real*8  temperho,temperature_r_z_i

c-----locals

 
      if(i.gt.nbulk)then
        write(*,*)'in tempe i.gt.nbulk: i,nbulk',i,nbulk
	stop
      endif

cBH151016: As discoverd by Syun'ichi Shiraiwa, this bug precludes use of
cBH151016: the temperature profiles outside the LCFS
c      if ((rho.gt.(1.d0-1.d-10)).and.(i_edge_dens_rz_mesh.gt.2)) then
      if ((rho.gt.(1.d0-1.d-10)).and.(i_edge_dens_rz_mesh.ge.2)) then
c----------------------------------------------------------
c        temperature outside LCFS 
c        using 2D spline at RZ mesh
c-----------------------------------------------------------
c           write(*,*)'in tempe temperaure_r_z_i'
            tempe=temperature_r_z_i(z,r,phi,i)
c           write(*,*)'in tempe after temperature_r_z_i '
c         endif
      else
	 tempe=temperho(rho,i)
      endif
cc      write(*,*)'in dens.for rho,i,tempe=',rho,i,tempe

      return
      end

      real*8 function temperature_r_z_i(z,r,phi,i)
c----------------------------------------------------------------------
c     calculate temperature using spline coefficients for temperature_r_z
c---------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'

c-----input
      real*8 z,r,phi !space coordinates [m]/r0x
      integer i      !the number of plasma specie 
  

c-----locals
      integer idm

c-----externals
      real*8 terp2p_Sm

c      write(*,*)'function temperature_r_z_i'
            
      idm=nxeqd_add_a
c      write(*,*)'before terp2p_Sm'       
      temperature_r_z_i=terp2p_Sm(r,z,
     &                    nxeqd_add,rr_add,
     &                    nyeqd_add,zz_add,
     &            temperature_r_z(1,1,i),temperature_r_z_rr(1,1,i),
     &            temperature_r_z_zz(1,1,i),
     &            temperature_r_z_rrzz(1,1,i),idm,0,0,1,nxeqd_add_a) 
c       write(*,*)'after terp2p_Sm'      
      return
      end


c        *********************zeff  ***************************
c        *                        -                           *
c        * this function calculates the Z_eff profile         *
c        * as a function  from psi(z,r,phi) 		      *
c        *  rho=psi-psimag 				      *
c        *  a=psilim-psimag inside the function b()           *
c        *  these variables are in common 'one'
c        ******************************************************
c--------------------------------- ---------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c       rho is from common/one/
c------------------------------------------------------------------
      double precision FUNCTION zeff(z,r,phi)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----externals
      real*8  zeffrho,zeff_r_z_f
c-----input
      real*8 z,r,phi !space coordinates [m]/r0x

      if ((rho.gt.(1.d0-1.d-10)).and.(i_edge_dens_rz_mesh.eq.2)) then
c----------------------------------------------------------
c        zeff outside LCFS 
c        using 2D spline at RZ mesh
c-----------------------------------------------------------
c           write(*,*)'in zeff zeff_r_z '
            zeff=zeff_r_z_f(z,r,phi)
c           write(*,*)'in zeff after zeff_r_z '
      else
	  zeff=zeffrho(rho)
      endif

cc      write(*,*)'in dense.for zeff=',zeff
      return
      END


 

      real*8 function zeff_r_z_f(z,r,phi)
c----------------------------------------------------------------------
c     calculate zeff using spline coefficients for zeff_r_z
c---------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'

c-----input
      real*8 z,r,phi !space coordinates [m]/r0x

c-----locals
      integer idm

c-----externals
      real*8 terp2p_Sm

c      write(*,*)'function zeff_r_z'
            
      idm=nxeqd_add_a
c      write(*,*)'before terp2p_Sm'       
      zeff_r_z_f=terp2p_Sm(r,z,
     &                    nxeqd_add,rr_add,
     &                    nyeqd_add,zz_add,
     &            zeff_r_z, zeff_r_z_rr,
     &            zeff_r_z_zz,
     &            zeff_r_z_rrzz,idm,0,0,1,nxeqd_add_a) 
c       write(*,*)'after terp2p_Sm'      
      return
      end



      subroutine splcoef_temperature_r_z
c----------------------------------------------------------------------
c     calculate spline coefficients for temperature_r_z
c---------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input

c-----locals
      integer dens_nwka  

c     in general,need dens_nwka=1+3*max(nxeqd_add_a,nyeqd_add_a)
      parameter (dens_nwka=1+3*nxeqd_add_a)   
      
      real*8      
     &wk_dens(dens_nwka),dens_norm,dens_spline 

      integer ibd(4),idm
      integer i,j,k

c-----externals
      real*8 terp2p_Sm

      write(*,*)'splcoef_temperature_r_z'
      write(*,*)'nxeqd_add_a,nyeqd_add_a,dens_nwka',
     &           nxeqd_add_a,nyeqd_add_a,dens_nwka
 
c------------------------------------------------------------------
c     check parameter dens_nwka value 
c------------------------------------------------------------------
      if (dens_nwka.lt.(1+3*nxeqd_add_a)) then
         write(*,*)'**************************************' 
         write(*,*)'chamber_wall,f in  splcoef_temperatur_r_z'
         write(*,*)'dens_nwka.lt.(1+3*nxeqd_add_a)'
         write(*,*)'it shoud be'
         write(*,*)'dens_nwka=1+3*max(nxeqd_add_a,nyeqd_add_a)'
         write(*,*)'dens_nwka,nxeqd_add_a,nyeqd_add_a',
     &              dens_nwka,nxeqd_add_a,nxeqd_add_a
         write(*,*)'Please change parameter in splcoef_temperature_r_z'
         write(*,*)'for parameter (dens_nwka=1+3*nxeqd_add_a)'
         write(*,*)'and recompile code'  
         write(*,*)'**************************************' 
         stop
      endif


      if (dens_nwka.lt.(1+3*nyeqd_add_a)) then
         write(*,*)'**************************************' 
         write(*,*)'chamber_wall,f in  splcoef_temperature_r_z'
         write(*,*)'dens_nwka.lt.(1+3*nyeqd_add_a)'
         write(*,*)'it shoud be'
         write(*,*)'dens_nwka=1+3*max(nxeqd_add_a,nyeqd_add_a)'
         write(*,*)'dens_nwka,nxeqd_add_a,nxyqd_add_a',
     &              dens_nwka,nxeqd_add_a,nxeqd_add_a
         write(*,*)'Please change parameter in splcoef_temperature_r_z'
         write(*,*)'for parameter (dens_nwka=1+3*nyeqd_add_a)'
         write(*,*)'and recompile code'  
         write(*,*)'**************************************' 
         stop
      endif
      
c-----creates 2D spline coefficients for temperature_r_z functions
      ibd(1)=2
      ibd(2)=2
      ibd(3)=2
      ibd(4)=2

      idm=nxeqd_add_a

      write(*,*)'chamber_wall.f in subroutine splcoef_temperature_r_z'
      write(*,*)'nbulk=',nbulk
      write(*,*)'nxeqd_add,nyeqd_add',nxeqd_add,nyeqd_add
      write(*,*)'nxeqd_add_a,nyeqd_add_a',nxeqd_add_a,nyeqd_add_a
      

      do k=1,nbulk

c        write(*,*)'k',k
       
         do j=1,nyeqd_add_a
            temperature_r_z_rr(1,j,k)=0.d0
            temperature_r_z_rr(nxeqd_add_a,j,k)=0.d0
         enddo    

          do i=1,nxeqd_add_a
            temperature_r_z_zz(i,1,k)=0.d0
            temperature_r_z_zz(i,nyeqd_add_a,k)=0.d0
          enddo

c          write(*,*)'before coeff2_Sm(nxeqd_add,rr_'

          call coeff2_Sm(nxeqd_add,rr_add,nyeqd_add,zz_add,
     &      temperature_r_z(1,1,k),
     &      temperature_r_z_rr(1,1,k),temperature_r_z_zz(1,1,k),
     &      temperature_r_z_rrzz(1,1,k),
     &      idm,ibd,wk_dens,nxeqd_add_a)

c          write(*,*)'after coeff2_Sm(nxeqd_add,rr_'
       enddo !k
   

c-test
c      do k=1,nbulk    

c         write(*,*)'test k',k

c         tepm_norm=0.d0    
c         do i=1,nxeqd_add
c            do j=1,nyeqd_add
c               write(*,*)'before terp2p k,j,i,temperature_r_z(i,j,k)'
c     &                                 ,k,j,i,temperature_r_z(i,j,k)

c              write(*,*)'before terp2p i,j',i,j
c              temp_spline=terp2p_Sm(rr_add(i),zz_add(j),
c     &                    nxeqd_add,rr_add,
c     &                    nyeqd_add,zz_add,
c     &             temperature_r_z(1,1,k),temperature_r_z_rr(1,1,k),
c     &             temperature_r_z_zz(1,1,k),
c     &             temperature_r_z_rrzz(1,1,k),idm,0,0,1,nxeqd_add_a) 
c              write(*,*)'after terp2p temp_spline',temp_spline
         
c              temp_norm=temp_norm+dabs(temp_spline-temperature_r_z(i,j,k))

c            enddo
c         enddo 
c         write(*,*)'k,temp_norm',k,temp_norm
c      enddo      
c-endtest

      return
      end

      subroutine splcoef_zeff_r_z
c----------------------------------------------------------------------
c     calculate spline coefficients for zeff_r_z
c---------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input

c-----locals
      integer dens_nwka  

c     in general,need dens_nwka=1+3*max(nxeqd_add_a,nyeqd_add_a)
      parameter (dens_nwka=1+3*nxeqd_add_a)   
      
      real*8      
     &wk_dens(dens_nwka),dens_norm,dens_spline 

      integer ibd(4),idm
      integer i,j

c-----externals
      real*8 terp2p_Sm

      write(*,*)'splcoef_zeff_r_z'
      write(*,*)'nxeqd_add_a,nyeqd_add_a,dens_nwka',
     &           nxeqd_add_a,nyeqd_add_a,dens_nwka
 
c------------------------------------------------------------------
c     check parameter dens_nwka value 
c------------------------------------------------------------------
      if (dens_nwka.lt.(1+3*nxeqd_add_a)) then
         write(*,*)'**************************************' 
         write(*,*)'chamber_wall,f in  splcoef_zeff_r_z'
         write(*,*)'dens_nwka.lt.(1+3*nxeqd_add_a)'
         write(*,*)'it shoud be'
         write(*,*)'dens_nwka=1+3*max(nxeqd_add_a,nyeqd_add_a)'
         write(*,*)'dens_nwka,nxeqd_add_a,nyeqd_add_a',
     &              dens_nwka,nxeqd_add_a,nxeqd_add_a
         write(*,*)'Please change parameter in splcoef_zefft_r_z'
         write(*,*)'for parameter (dens_nwka=1+3*nxeqd_add_a)'
         write(*,*)'and recompile code'  
         write(*,*)'**************************************' 
         stop
      endif


      if (dens_nwka.lt.(1+3*nyeqd_add_a)) then
         write(*,*)'**************************************' 
         write(*,*)'chamber_wall,f in  splcoef_zeff_r_z'
         write(*,*)'dens_nwka.lt.(1+3*nyeqd_add_a)'
         write(*,*)'it shoud be'
         write(*,*)'dens_nwka=1+3*max(nxeqd_add_a,nyeqd_add_a)'
         write(*,*)'dens_nwka,nxeqd_add_a,nxyqd_add_a',
     &              dens_nwka,nxeqd_add_a,nxeqd_add_a
         write(*,*)'Please change parameter in splcoef_zeff_r_z'
         write(*,*)'for parameter (dens_nwka=1+3*nyeqd_add_a)'
         write(*,*)'and recompile code'  
         write(*,*)'**************************************' 
         stop
      endif
      
c-----creates 2D spline coefficients for zeff_r_z functions
      ibd(1)=2
      ibd(2)=2
      ibd(3)=2
      ibd(4)=2

      idm=nxeqd_add_a

      write(*,*)'chamber_wall.f in subroutine splcoef_zeff_r_z'
      write(*,*)'nbulk=',nbulk
      write(*,*)'nxeqd_add,nyeqd_add',nxeqd_add,nyeqd_add
      write(*,*)'nxeqd_add_a,nyeqd_add_a',nxeqd_add_a,nyeqd_add_a
      

    
       
      do j=1,nyeqd_add_a
         zeff_r_z_rr(1,j)=0.d0
         zeff_r_z_rr(nxeqd_add_a,j)=0.d0
      enddo    

      do i=1,nxeqd_add_a
         zeff_r_z_zz(i,1)=0.d0
         zeff_r_z_zz(i,nyeqd_add_a)=0.d0
      enddo

c      write(*,*)'in splcoef_zeff before coeff2_Sm(nxeqd_add,rr_'

c      do j=1,nyeqd_add 
c         do i=1,nxeqd_add_a
c            write(*,*)'j,i,zeff_r_z(i,j)',j,i,zeff_r_z(i,j)
c         enddo
c      enddo    

c      stop 'in dense.f splcoef_zeff'



      call coeff2_Sm(nxeqd_add,rr_add,nyeqd_add,zz_add,
     &      zeff_r_z,
     &      zeff_r_z_rr,zeff_r_z_zz,
     &      zeff_r_z_rrzz,
     &      idm,ibd,wk_dens,nxeqd_add_a)

c     write(*,*)'after coeff2_Sm(nxeqd_add,rr_'
     
   

c-test

c         zeff_norm=0.d0    
c         do i=1,nxeqd_add
c            do j=1,nyeqd_add
c               write(*,*)'before terp2p k,j,i,zeff_r_z(i,j)'
c     &                                 ,k,j,i,zeff_r_z(i,j)

c              write(*,*)'before terp2p i,j',i,j
c              temp_spline=terp2p_Sm(rr_add(i),zz_add(j),
c     &                    nxeqd_add,rr_add,
c     &                    nyeqd_add,zz_add,
c     &             zeff_r_z,zeff_r_z_rr,
c     &             zeff_r_z_zz,
c     &             zeff_r_z_rrzz,idm,0,0,1,nxeqd_add_a) 
c              write(*,*)'after terp2p zeff_spline',zeff_spline
         
c              zeff_norm=zeff_norm+dabs(zeff_spline-zeff_r_z(i,j))

c            enddo
c         enddo 
c         write(*,*)'zeff_norm',zeff_norm     
c-endtest

      return
      end
