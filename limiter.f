
      subroutine rz_lim_on_rho_geom_theta(rho_geom,theta_pol,
     &r_mag,z_mag,z,r)
c-----calculate (r,z) coordinate at the point with coordinates 
c     (rho_geom,theta_pol)
  
      implicit none
c-----input
      real*8
     &r_mag,z_mag,        !magnetic axis coordinates
     &rho_geom,           !geometrical radius from the magnetic axis
     &theta_pol           !poloidal angle [radians], measured
                          !anticlockwise from the meridional plane
c-----output
      real*8 r,z          !space coordinates

c      write(*,*)'in rz_lim_on_rho_geom_theta,rho_geom,theta_pol,
c     &r_mag,z_mag',rho_geom,theta_pol,r_mag,z_mag

      z = z_mag + rho_geom*dsin(theta_pol)
      r = r_mag + rho_geom*dcos(theta_pol)

c      write(*,*)'in rz_lim_on_rho_geom_theta 
c     &dsin(theta_pol),dcos(theta_pol),z,r',
c     &dsin(theta_pol),dcos(theta_pol),z,r

      return
      end

      subroutine RHS_limiter_line(theta_pol,u,deru)
c-----------------------------------------------------------
c     calculate RHS of the equation for field line at LCFS
c
c     d_rho_geom/d_theta_pol = RHS_limiter_lne=
c     =rho*[sin(theta_pol)d_psi_dr-cos(theta_pol)*dpsi_dz]/
c      [cos(theta_pol)d_psi_dr+sin(theta_pol)*dpsi_dz]
c---------------------------------------------------
      implicit none

      real*8 
     &z_mag_l,r_mag_l  
      common /RHS_limiter_line_common/
     &z_mag_l,r_mag_l                     !magnetix axis cordinates
c-----input
      real*8 theta_pol, !poloidal angle along field line
     & u(1),deru(1) 

c-----local 
      real*8
     &r,z,             ! coordinates at the fild line
     &dpsidz,dpsidr,   !derivatives from the polidal flux
     &phi,rho_geom_u

      rho_geom_u=u(1)
c-----calculate z,r coordinates at given rho_geom_u,theta_pol,r_mag,z_mag
       
      call rz_lim_on_rho_geom_theta(rho_geom_u,theta_pol,
     &r_mag_l,z_mag_l,z,r)

c      write(*,*)'RHS rho_geom_u,theta_pol,z,r',rho_geom_u,theta_pol,z,r

c-----calculate derivatives from the poloidal flux at LCFS
c     at given r,z
      phi=0.d0                       !toroidal angle can be arbitrary
      call dpsidzdr(z,r,phi,dpsidz,dpsidr)
      
      deru(1) = rho_geom_u*
     &(dsin(theta_pol)*dpsidr - dcos(theta_pol)*dpsidz)/
     &(dcos(theta_pol)*dpsidr + dsin(theta_pol)*dpsidz)

c      write(*,*)'dpsidr,dpsidz,deru(1)', dpsidr,dpsidz,deru(1)

      return
      return
      end
 
      subroutine calc_limiter_field_line(psilim,nlimit,
     &z_mag,r_mag,rlimit,zlimit,theta_pol_limit,rho_geom_limit)

c-----------------------------------------------------
c     calculate limiter field line at LCFS at the poloidal flux 
c     psi=psilim
c     It will create limiter arrays: rlimit(i),zlimit(i)
c     at the poloidal uniform mesh 0=<theta_pol_lim(i)<2*pi
c     i=1,nlimit,   with step h_theta=2pi/(nlimit-1)
c
c     Binary mejhod using
c
c     Not works
c     Runge-Kutta method  of  4th-order                  !
c     with the constant step integration.   
c     subroutine drkgs(prmt,u,deru,ndim,ihlf,fct,outp,aux)

      implicit none
c-----input
      real*8 psilim,     !poloidal flux at limiter LCFS
     &z_mag,r_mag        !magnetic axis coordinates
      integer nlimit     !number of limiter arrays points 
   
c-----output
      real*8
     & rlimit(nlimit),  !limiter points coordinates
     & zlimit(nlimit),
     & theta_pol_limit(nlimit),   
     & rho_geom_limit(nlimit)

c-----external
      real*8 fpsi
      external RHS_limiter_line
      
c-----local
      integer i

      real*8 r_0,step_r,r_l,r_r,r_c,psi1,psi2,
     &epspsi              !accuray

      real*8 step_rho,rho_l,rho_r,z_r,z_l,rho_c,z_c

c-----for RK solver
      real*8
     &rho_lim(1),d_rho_d_theta(1),aux_lim(8,1),
     &pi,step_theta_pol
      integer ndim 
c-------------------------------------------------------
      real*8   
     &z_mag_l,r_mag_l  
      common /RHS_limiter_line_common/
     &z_mag_l,r_mag_l                  !magnetix axis cordinates
c----------------------------------------------
c     additional RK steps
      integer nlimit_add,k
      real*8 step_add,theta_pol_add
c----------------------------------------------
      nlimit_add=10

      z_mag_l=z_mag
      r_mag_l=r_mag

      pi=4.d0*datan(1.d0)

      epspsi=1.d-6        !set accuracy

c-----------------------------------------------------------
c     calculate the major radius of the initial point
c     at the field line 
c     r_0 at LCFS at z=yma (at theta_pol=0) 
c
c     To find r_0 we will solve the equation
c     F(r)=psi(r,z=z_mag)-psilim=0
c    
c-----------------------------------------------------------
      step_r=r_mag/100.d0
      r_r=r_mag

 10   r_l=r_r
      r_r=r_r+step_r
      if (fpsi(r_r,z_mag).lt.psilim) then 
        goto 10        
      else
c---------------------------------------------------
c       Boundaries r_l,r_r were found
c       fpsi(r_l,z_mag) < psilim <fpsi(r_r,z_mag)
c       
c       Solution of the equation F(r)=psi(r,z=z_mag)-psilim=0 at 
c       using the binary method
c--------------------------------------
 
        do while ((r_r-r_l).gt.epspsi)
             r_c=r_l+(r_r-r_l)*0.5d0          
             psi1=fpsi(r_c,z_mag)-psilim
             psi2=fpsi(r_r,z_mag)-psilim
             if ((psi1*psi2).gt.0) then
               r_r=r_c
             else
               r_l=r_c
             end if
        end do

        write(*,*)'initial point at the fiel line r_c',r_c,
     &  'fpsi(r_c,z_mag)-psilim',fpsi(r_c,z_mag)-psilim
        
c       -----------------------------------------------
c       end of the binary method gives r_0=r_c
c       -----------------------------------------------
      endif

c-----calculate limiter points at the line field at LCFS 
c     initial point at the field line: (r_c,z_mag) with theta_pol=0

      ndim=1
      step_theta_pol=2.d0*pi/(nlimit-1)
      rho_lim(1)=(r_c-r_mag) !at the initial =point
      rlimit(1)=r_c
      zlimit(1)=z_mag
      theta_pol_limit(1)=0.d0 ![radians]
      rho_geom_limit(1)=rho_lim(1)

      do i=2,nlimit 
         
         theta_pol_limit(i)=theta_pol_limit(i-1)+step_theta_pol

c-----------------------------------------------------------
c        calculate the point  r_0,z_0 at LCFS  
c        r=r_mag+rho_geom*cos(theta_pol_limit)
c        z=z_mag+rho_geom*sin(theta_pol_limit)
c
c        from the solution of the field line equation 
c        using RK solver
c        
c-----------------------------------------------------------
      
c         write(*,*)'i,theta_pol_limit(i)',i,theta_pol_limit(i)

cSAP091004
         if (nlimit_add.ge.1)then
     
             step_add=step_theta_pol/nlimit_add
             theta_pol_add=theta_pol_limit(i-1)

             do k=1,nlimit_add
c                step_add=step_theta_pol/nlimit_add
c---------------additional RK steps for better accuracy
cSAP091004

                call drkgs0_limiter( theta_pol_add,step_add,
     &                         rho_lim,d_rho_d_theta,ndim,
     &                         RHS_limiter_line,aux_lim)

                theta_pol_add=theta_pol_add+step_add

c                write(*,*)'after drkgs0 k, rho_lim',i, rho_lim
             enddo
         else
             call drkgs0_limiter(theta_pol_limit(i-1),step_theta_pol,
     &               rho_lim,d_rho_d_theta,ndim,
     &               RHS_limiter_line,aux_lim)  
         endif

c        write(*,*)'after drkgs0 i, rho_lim',i, rho_lim
         call rz_lim_on_rho_geom_theta(rho_lim(1),theta_pol_limit(i),
     &   r_mag,z_mag,zlimit(i),rlimit(i))

         rho_geom_limit(i)=rho_lim(1)
        
      enddo

      write(*,*)'in limiter.f in calc_limiter_field_line nlimit',nlimit

      do i=1,nlimit
        write(*,*)'i,theta_pol_limit(i),rho_geom_limit(i)',
     &             i,theta_pol_limit(i),rho_geom_limit(i)
        write(*,*)'zlimit(i),rlimit(i)',zlimit(i),rlimit(i)
        write(*,*)'psilim,fpsi(rlimit(i),zlimit(i)',
     &             psilim,fpsi(rlimit(i),zlimit(i))
        write(*,*)'psilim-fpsi(rlimit(i),zlimit(i)',
     &             psilim-fpsi(rlimit(i),zlimit(i))
      enddo

      return
      end




c        ********************** drkgs0_limiter *******************
c         4_th order Runge-Kutta subroutine with constant time 
c         step finds the solution of
c         the system of the ordinary differential equations
c         only for one time step         
c        ****************************************************
c         It is used for solution of the field linwe equation
c         d_rho_geomer/df_theta=fct  at ndim=1
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c          Runge-Kutta method  of  4th-order                     !
c          with the constant step integration
c          It makes only one step!
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       2. input parameters                                      !
c                                                                !
c          h     - time step                                     ! 
c                                                                !
c          u     - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  after one time stept                         !
c                                                                !
c          deru  - the  vector of derivatives of u at point t	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array.                          !
c                                                                !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system. it has not to change the  values  of    !
c                t and u. its formal parameters are:             !
c                t, u, deru.                                     !
c                                                                !
c       			  		 !
c    u(1)=z,u(2)=r,u(3)=phi,u(4)=nz,u(5)=nr,u(6)=m               !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs0_limiter(t,h,u,deru,ndim,fct,aux)
      !implicit double precision (a-h,o-z)
      implicit none !YuP[2020-01-14]
      integer i,ndim !YuP[2020-01-14]
      real*8 t,h,u,deru,aux, up !YuP[2020-01-14]
      dimension u(1),deru(1),aux(8,1),up(1)

c      write(*,*)'in Runge -Kutta drkgs0'
c      write(*,*)'t,h',t,h

      call fct(t,u,deru)
c      write(*,*)'t,u,deru',t,u,deru
      do i=1,ndim
         aux(1,i)=h*deru(i)
	 aux(5,i)=u(i)+0.5d0*aux(1,i)
 	 up(i)=aux(5,i)
      enddo
c------------------------------------------------------------------
      call fct(t+0.5d0*h,up,deru)
      do i=1,ndim
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5d0*aux(2,i)
	 up(i)=aux(5,i)
      enddo
c--------------------------------------------------------------------
      call fct(t+0.5d0*h,up,deru)
      do i=1,ndim
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
 	 up(i)=aux(5,i)
      enddo
c--------------------------------------------------------------------
      call fct(t+h,up,deru)
      do i=1,ndim
         aux(4,i)=h*deru(i)
      enddo
c--------------------------------------------------------------------
      do i=1,ndim
         up(i)=u(i)+1.d0/6.d0*(aux(1,i)+2.d0*aux(2,i)+2.d0*aux(3,i)
     *          +aux(4,i))
      enddo

      do i=1,ndim
         u(i)=up(i)
      enddo

      return
      end

 
      subroutine theta_pol_at_angles(theta_right_top,
     &theta_left_top,theta_left_bottom,theta_right_bottom)
c----------------------------------------------------- 
c     calculate four po loidal angles 
c     at the corners of the eqdsk rectangle
c-----------------------------------------------------       

      implicit none

      include 'param.i'
      include 'three.i'
      include 'fourb.i'
c-----output
      real*8
     &theta_right_top,  !corner's angles [radians] 0=< angle < 2pi
     &theta_left_top,
     &theta_left_bottom,
     &theta_right_bottom

c-----local
      real*8 rho_geom,pi

      pi=4.d0*atan(1.d0)

c-----at right top corner
      rho_geom=dsqrt( (rr(nxeqd)-xma)**2 + (zz(nyeqd)-yma)**2)
      theta_right_top=dacos((rr(nxeqd)-xma)/rho_geom) 

c-----at left top corner
      rho_geom=dsqrt((rr(1)-xma)**2 + (zz(nyeqd)-yma)**2)
      theta_left_top=dacos((rr(1)-xma)/rho_geom) 

c-----at left bottom corner
      rho_geom=dsqrt((rr(1)-xma)**2 + (zz(1)-yma)**2)
      theta_left_bottom=2.d0*pi - dacos((rr(1)-xma)/rho_geom) 

c-----at right bootom corner
      rho_geom=dsqrt((rr(nxeqd)-xma)**2 + (zz(1)-yma)**2)
      theta_right_bottom=2.d0*pi - dacos((rr(nxeqd)-xma)/rho_geom) 

      return
      end

      real*8  function max_rho_geom(theta_pol)
c-----calculate maximal geometrical radius at the eqsdk rectangle 
c     at given poloidal angle theta_pol [radian]
c     geometrical radius is a distance from the magnetic axis(xma,yma)
c----------------------------------------------------------------------
      implicit none

      include 'param.i'
      include 'three.i'
      include 'fourb.i'
c-----input
      real*8 theta_pol ! poloidal angle measued anrti-clockwise
                       ! from the equatorial plane (z=yma)
                       ! the_pol=0 at the z=yma at outward tokomak side  

c------locals
      real*8
     &theta_right_top,  !corner's angles [radians] 0=< angle < 2pi
     &theta_left_top,
     &theta_left_bottom,
     &theta_right_bottom,
     &pi,theta_pol_l,ratio,z,r

      integer k

      pi=4.d0*datan(1.d0)

c      write(*,*)'in max_rho_geom theta_pol',theta_pol

c----------------------------------------------------- 
c     calculate four poloidal angles 
c     at corners of the eqdsk rectangle
c-----------------------------------------------------       
      call theta_pol_at_angles(theta_right_top,
     &theta_left_top,theta_left_bottom,theta_right_bottom)

c      write(*,*)'theta_right_top,
c     &theta_left_top,theta_left_bottom,theta_right_bottom',
c     &theta_right_top,
c     &theta_left_top,theta_left_bottom,theta_right_bottom

c------------------------------------------------------
c     transform the angle teta_pol to the region 
c     0=< teta_pol_l <2*pi
c-------------------------------------------------------
      ratio=theta_pol/(2.d0*pi)
      theta_pol_l=theta_pol

      if (ratio.gt.1) then
          k=int(ratio)
          theta_pol_l=theta_pol-2.d0*pi*k
      endif
     
      if (ratio.lt.0) then
          k=int(ratio)-1
          theta_pol_l=theta_pol-2.d0*pi*k
      endif

c      write(*,*)'in max_rho_geom theta_pol_l',theta_pol_l

c-------------------------------------------------------
      if((0.d0.le.theta_pol_l).and.(theta_pol_l.le.theta_right_top))then
         r=rr(nxeqd)-xma
         max_rho_geom=r/dcos(theta_pol_l)
      elseif(theta_pol_l.le.theta_left_top) then
cSAP101004
c         z=zz(nxeqd)-yma
         z=zz(nyeqd)-yma
         max_rho_geom=z/dsin(theta_pol_l)

c        write(*,*)'in max_rho_geom theta_pol_l.le.theta_left_top'
c        write(*,*)'nxeqd,zz(nxeqd),yma',nxeqd,zz(nxeqd),yma
c        write(*,*)'z,dsin(theta_pol_l),max_rho_geom',
c     &             z,dsin(theta_pol_l),max_rho_geom

      elseif(theta_pol_l.le.theta_left_bottom) then
         r=rr(1)-xma
         max_rho_geom=dabs(r/dcos(theta_pol_l))
      elseif(theta_pol_l.le.theta_right_bottom) then
cSAP101004
c         z=zz(nxeqd)-yma
         z=zz(nyeqd)-yma
         max_rho_geom=dabs(z/dsin(theta_pol_l))
      else
         r=rr(nxeqd)-xma
         max_rho_geom=dabs(r/dcos(theta_pol_l))
      endif  

c      write(*,*)'in max_rho_geo max_rho_geom',max_rho_geom

      return
      end

      subroutine calc_limiter_binary (psilim,nlimit,
     &z_mag,r_mag,rlimit,zlimit,theta_pol_limit,rho_geom_limit)

c-----------------------------------------------------
c     calculate limiter field line at LCFS at the poloidal flux 
c     psi=psilim
c     It will create limiter arrays: rlimit(i),zlimit(i)
c     at the poloidal uniform mesh 0=<theta_pol_lim(i)<2*pi
c     i=1,nlimit,   with step h_theta=2pi/(nlimit-1)
c
c     Binary method is used. 

      implicit none

c-----input
      real*8 psilim,     !poloidal flux at limiter LCFS
     &z_mag,r_mag        !magnetic axis coordinates
      integer nlimit     !number of limiter arrays points 
   
c-----output
      real*8
     & rlimit(nlimit),  !limiter points coordinates
     & zlimit(nlimit),
     & theta_pol_limit(nlimit),   
     & rho_geom_limit(nlimit)

c-----external
      real*8 fpsi,max_rho_geom
      
c-----local
      integer i

      real*8 
     &step_r,r_l,r_r,r_c,psi1,psi2,
     &epspsi,                                                   !accuray
     &pi,rho_max_theta,rho_l,rho_r,step_rho,step_theta_pol,
     &tr,tl,t,z,r,rho_c_m,rho_c_p,rho_c,rho_m,rho_p,
     &del_psi_m,del_psi_c,del_psi_p,del_psi_c_p,del_psi_c_m

      integer n_rho_steps

c----------------------------------------------
      n_rho_steps=100  !number of steps in rho

c      write(*,*)'calc_limiter_binary n_rho_steps',n_rho_steps

      pi=4.d0*datan(1.d0)

      epspsi=1.d-6        !set accuracy

c-----------------------------------------------------------
c     calculate the major radius of the initial point
c     at the field line 
c     r_0 at LCFS at z=yma (at theta_pol=0) 
c
c     To find r_0 we will solve the equation
c     F(r)=psi(r,z=z_mag)-psilim=0
c    
c-----------------------------------------------------------
      step_r=r_mag/100.d0
      r_r=r_mag

 10   continue
      r_l=r_r
      r_r=r_r+step_r
      if (fpsi(r_r,z_mag).lt.psilim) then 
        goto 10        
      else
c---------------------------------------------------
c       Boundaries r_l,r_r were found
c       fpsi(r_l,z_mag) < psilim <fpsi(r_r,z_mag)
c       
c       Solution of the equation F(r)=psi(r,z=z_mag)-psilim=0 at 
c       using the binary method
c--------------------------------------
 
        do while ((r_r-r_l).gt.epspsi)
             r_c=r_l+(r_r-r_l)*0.5d0          
             psi1=fpsi(r_c,z_mag)-psilim
             psi2=fpsi(r_r,z_mag)-psilim
             if ((psi1*psi2).gt.0) then
               r_r=r_c
             else
               r_l=r_c
             end if
        end do

        write(*,*)'initial point r_c',r_c,
     &  'fpsi(r_c,z_mag)-psilim',fpsi(r_c,z_mag)-psilim
        
c       -----------------------------------------------
c       end of the binary method gives r_0=r_c
c       -----------------------------------------------
      endif

c-----calculate limiter points at the line field at LCFS 
c     initial point at the field line: (r_c,z_mag) with theta_pol=0

      step_theta_pol=2.d0*pi/(nlimit-1)
      rho_geom_limit(1)=(r_c-r_mag) !at the initial =point
      rlimit(1)=r_c
      zlimit(1)=z_mag
      theta_pol_limit(1)=0.d0 ![radians]
      

      do i=2,nlimit 
         
         theta_pol_limit(i)=theta_pol_limit(i-1)+step_theta_pol

c-----------------------------------------------------------
c        calculate the point  r_0,z_0 at LCFS  
c        r=r_mag+rho_geom*cos(theta_pol_limit)
c        z=z_mag+rho_geom*sin(theta_pol_limit)
c
c        using bynary solver
c        
c-----------------------------------------------------------      
         write(*,*)'i,theta_pol_limit(i)',i,theta_pol_limit(i)

         rho_max_theta=max_rho_geom(theta_pol_limit(i))
 
         write(*,*)'rho_max_theta, step_rho',rho_max_theta, step_rho

         step_rho=rho_max_theta/dble(n_rho_steps)

         rho_m=0.d0
         rho_c=rho_m+step_rho 

 15      continue

         call rz_lim_on_rho_geom_theta(rho_m,
     &        theta_pol_limit(i),r_mag,z_mag,z,r)
         del_psi_m=fpsi(r,z)-psilim

c         write(*,*)'rho_m,z,r,del_psi_m',rho_m,z,r,del_psi_m    

         call rz_lim_on_rho_geom_theta(rho_c,
     &        theta_pol_limit(i),r_mag,z_mag,z,r)
         del_psi_c=fpsi(r,z)-psilim

c         write(*,*)'rho_c,z,r,del_psi_c',rho_c,z,r,del_psi_c

         if (del_psi_m*del_psi_c.lt.0d0) then
c-----------LCFS point at theta_pol_limit(i) is between rho_m and rho_c
            tl=rho_m
            tr=rho_c            
            do while ((tr-tl).gt.epspsi)

               t=tl+(tr-tl)*0.5d0
               call rz_lim_on_rho_geom_theta(t,
     &         theta_pol_limit(i),r_mag,z_mag,z,r)
               psi1=fpsi(r,z)-psilim

               call rz_lim_on_rho_geom_theta(tr,
     &         theta_pol_limit(i),r_mag,z_mag,z,r)
               psi2=fpsi(r,z)-psilim
                            
               if ((psi1*psi2).gt.0) then
                  tr=t
               else
                  tl=t
               end if
            end do

            rlimit(i)=r
            zlimit(i)=z

            write(*,*)'bin:i,theta_pol_limit(i)',i,theta_pol_limit(i)
            write(*,*)'rlimit(i),zlimit(i)', rlimit(i),zlimit(i)
            write(*,*)'fpsi(rlimit(i),zlimit(i))-psilim',
     &                 fpsi(rlimit(i),zlimit(i))-psilim
            goto 20

         endif

         rho_p=rho_c+step_rho
           
         call rz_lim_on_rho_geom_theta(rho_p,
     &        theta_pol_limit(i),r_mag,z_mag,z,r)
         del_psi_p=fpsi(r,z)-psilim

         write(*,*)'rho_p,z,r,del_psi_p',rho_p,z,r,del_psi_p
         if (del_psi_c*del_psi_p.lt.0d0) then
c-----------LCFS point at theta_pol_limit(i) is between rho_c and rho_p
            tl=rho_c
            tr=rho_p            
            do while ((tr-tl).gt.epspsi)

               t=tl+(tr-tl)*0.5d0
               call rz_lim_on_rho_geom_theta(t,
     &         theta_pol_limit(i),r_mag,z_mag,z,r)
               psi1=fpsi(r,z)-psilim

               call rz_lim_on_rho_geom_theta(tr,
     &         theta_pol_limit(i),r_mag,z_mag,z,r)
               psi2=fpsi(r,z)-psilim
                            
               if ((psi1*psi2).gt.0) then
                  tr=t
               else
                  tl=t
               end if
            end do

            rlimit(i)=r
            zlimit(i)=z

            write(*,*)'bin:i,theta_pol_limit(i)',i,theta_pol_limit(i)
            write(*,*)'rlimit(i),zlimit(i)', rlimit(i),zlimit(i)
            write(*,*)'fpsi(rlimit(i),zlimit(i))-psilim',
     &                 fpsi(rlimit(i),zlimit(i))-psilim
            goto 20

         endif

         if((dabs(del_psi_c).lt.dabs(del_psi_m)).and.
     &      (dabs(del_psi_c).lt.dabs(del_psi_p)))then

c-----------LCFS point at theta_pol_limit(i) is between rho_m and rho_p
c           poloidal flux has maximal value betweeen rho_m and rho_p
c------------------------------------------------------------------------
 30         continue

            if ((rho_p - rho_m).lt.epspsi) then
c--------------LCFS point is found
               call rz_lim_on_rho_geom_theta(rho_c,
     &              theta_pol_limit(i),r_mag,z_mag,z,r)
               rlimit(i)=r
               zlimit(i)=z
              
               write(*,*)'max:i,theta_pol_limit(i)',i,theta_pol_limit(i)
               write(*,*)'rlimit(i),zlimit(i)', rlimit(i),zlimit(i)
               write(*,*)'fpsi(rlimit(i),zlimit(i))-psilim',
     &                    fpsi(rlimit(i),zlimit(i))-psilim
               goto 20
            endif

           
            rho_c_m=0.5d0*(rho_m + rho_c)
            call rz_lim_on_rho_geom_theta(rho_c_m,
     &           theta_pol_limit(i),r_mag,z_mag,z,r)
            del_psi_c_m=fpsi(r,z)-psilim

            if((dabs(del_psi_c_m).lt.dabs(del_psi_m)).and.
     &         (dabs(del_psi_c_m).lt.dabs(del_psi_c)))then
c--------------LCFS point at theta_pol_limit(i) is between rho_m and rho_c
c              poloidal flux has maximal value betweeen rho_m and rho_c
               rho_p=rho_c
               del_psi_p= del_psi_c
               rho_c=rho_c_m
               del_psi_c=del_psi_c_m
               goto 30
            endif

            rho_c_p=0.5d0*(rho_c + rho_p)
            call rz_lim_on_rho_geom_theta(rho_c_p,
     &           theta_pol_limit(i),r_mag,z_mag,z,r)
            del_psi_c_p=fpsi(r,z)-psilim

            if((dabs(del_psi_c_p).lt.dabs(del_psi_p)).and.
     &         (dabs(del_psi_c_p).lt.dabs(del_psi_c)))then
c--------------LCFS point at theta_pol_limit(i) is between rho_c and rho_p
c              poloidal flux has maximal value betweeen rho_c and rho_p
               rho_m=rho_c
               del_psi_m=del_psi_c
               rho_c=rho_c_p
               del_psi_c=del_psi_c_p
               goto 30
            endif

c-----------LCFS point at theta_pol_limit(i) is between rho_c_m and rho_c_p
c           poloidal flux has maximal value betweeen rho_c_m and rho_c_p
            rho_m=rho_c_m
            del_psi_m= del_psi_c_m
            rho_p=rho_c_p
            del_psi_p=del_psi_c_p
            goto 30
         endif

         rho_m=rho_c 
         rho_c=rho_p 
          
         goto 15

 20      continue
      enddo

      write(*,*)'in limiter.f in calc_limiter_field_line nlimit',nlimit

      do i=1,nlimit
        write(*,*)'i,theta_pol_limit(i),rho_geom_limit(i)',
     &             i,theta_pol_limit(i),rho_geom_limit(i)
        write(*,*)'zlimit(i),rlimit(i)',zlimit(i),rlimit(i)
        write(*,*)'psilim,fpsi(rlimit(i),zlimit(i)',
     &             psilim,fpsi(rlimit(i),zlimit(i))
        write(*,*)'psilim-fpsi(rlimit(i),zlimit(i)',
     &             psilim-fpsi(rlimit(i),zlimit(i))
      enddo

      return
      end

      subroutine calc_r_z_psi_theta_binary(psi,
     &z_mag,r_mag,r,z,theta_pol)

c-----------------------------------------------------
c     calculate z,r on (psi,theta_pol)
c     It will create: r,z
c     at 0=<theta_pol<2*pi
c
c     Binary method is used. 

      implicit none

c-----input
      real*8 psi,        !poloidal flux
     &z_mag,r_mag,        !magnetic axis coordinates
     & theta_pol
   
c-----output
      real*8
     & r,z 
 
c-----external
      real*8 fpsi,max_rho_geom
      
c-----local
      integer i

      real*8 
     &step_r,r_l,r_r,r_c,psi1,psi2,
     &epspsi,                                                   !accuray
     &pi,rho_max_theta,rho_l,rho_r,step_rho,step_theta_pol,
     &tr,tl,t,rho_c_m,rho_c_p,rho_c,rho_m,rho_p,
     &del_psi_m,del_psi_c,del_psi_p,del_psi_c_p,del_psi_c_m,
     &rho_geom         
      integer n_rho_steps

c----------------------------------------------
      n_rho_steps=100  !number of steps in rho'

c      write(*,*)'calc_r_z_psi_theta_binary psi,z_mag,r_mag,theta_pol',
c     &psi,z_mag,r_mag,theta_pol
c      write(*,*)'calc_r_z_psi_theta_binary n_rho_steps',n_rho_steps

      pi=4.d0*datan(1.d0)

      epspsi=1.d-6        !set accuracy

c-----------------------------------------------------------
c     calculate the point  r_0,z_0 at psi,theta_pol  
c        r=r_mag+rho_geom*cos(theta_pol)
c        z=z_mag+rho_geom*sin(theta_pol)
c
c     using binary solver
c        
c-----------------------------------------------------------      

      rho_max_theta=max_rho_geom(theta_pol)
    
      step_rho=rho_max_theta/dble(n_rho_steps)

c      write(*,*)'rho_max_theta, step_rho',rho_max_theta, step_rho

      rho_m=0.d0
      rho_c=rho_m+step_rho 
       
 15   continue
 
c      write(*,*)'before rz_lim_on_rho_geom_theta rho_m,rho_c', 
c     & rho_m,rho_c

      call rz_lim_on_rho_geom_theta(rho_m,
     &        theta_pol,r_mag,z_mag,z,r)
      del_psi_m=fpsi(r,z)-psi

c      write(*,*)'rho_m,z,r,del_psi_m',rho_m,z,r,del_psi_m    

      call rz_lim_on_rho_geom_theta(rho_c,
     &        theta_pol,r_mag,z_mag,z,r)
         del_psi_c=fpsi(r,z)-psi

c      write(*,*)'rho_c,z,r,del_psi_c',rho_c,z,r,del_psi_c

      if (del_psi_m*del_psi_c.lt.0d0) then
c--------point at psi, theta_pol is between rho_m and rho_c

c         write(*,*)'calc_r_z_psi_theta_binary 1'

         tl=rho_m
         tr=rho_c            
         do while ((tr-tl).gt.epspsi)

            t=tl+(tr-tl)*0.5d0
            call rz_lim_on_rho_geom_theta(t,
     &         theta_pol,r_mag,z_mag,z,r)
            psi1=fpsi(r,z)-psi

            call rz_lim_on_rho_geom_theta(tr,
     &         theta_pol,r_mag,z_mag,z,r)
            psi2=fpsi(r,z)-psi
                            
            if ((psi1*psi2).gt.0) then
               tr=t
            else
               tl=t
            end if
         end do


c         write(*,*)'bin:theta_pol',theta_pol
c         write(*,*)'r,z', r,z
c         write(*,*)'fpsi(r,z)-psi',
c     &              fpsi(r,z)-psi
         goto 20

      endif

      rho_p=rho_c+step_rho
           
      call rz_lim_on_rho_geom_theta(rho_p,
     &        theta_pol,r_mag,z_mag,z,r)
      del_psi_p=fpsi(r,z)-psi

c      write(*,*)'rho_p,z,r,del_psi_p',rho_p,z,r,del_psi_p

      if (del_psi_c*del_psi_p.lt.0d0) then
c--------point at psi,theta_pol is between rho_c and rho_p

c         write(*,*)'calc_r_z_psi_theta_binary 2'

         tl=rho_c
         tr=rho_p            
         do while ((tr-tl).gt.epspsi)

            t=tl+(tr-tl)*0.5d0
            call rz_lim_on_rho_geom_theta(t,
     &         theta_pol,r_mag,z_mag,z,r)
            psi1=fpsi(r,z)-psi

            call rz_lim_on_rho_geom_theta(tr,
     &         theta_pol,r_mag,z_mag,z,r)
            psi2=fpsi(r,z)-psi
                            
            if ((psi1*psi2).gt.0) then
               tr=t
            else
               tl=t
            end if
         end do

c         write(*,*)'bin:theta_pol',theta_pol
c         write(*,*)'r,z', r,z
c         write(*,*)'fpsi(r,z)-psi',
c     &              fpsi(r,z)-psi
         goto 20

      endif

      if((dabs(del_psi_c).lt.dabs(del_psi_m)).and.
     &   (dabs(del_psi_c).lt.dabs(del_psi_p)))then

c--------point at psi,theta_pol is between rho_m and rho_p
c        poloidal flux has maximal value betweeen rho_m and rho_p
c------------------------------------------------------------------------
 30      continue

         if ((rho_p - rho_m).lt.epspsi) then
c-----------point is found

c            write(*,*)'calc_r_z_psi_theta_binary 3'

            call rz_lim_on_rho_geom_theta(rho_c,
     &              theta_pol,r_mag,z_mag,z,r)
                           
c            write(*,*)'max:theta_pol',theta_pol
c            write(*,*)'r,z', r,z
c            write(*,*)'fpsi(r,z)-psi',
c     &                 fpsi(r,z)-psi
            goto 20
         endif

           
         rho_c_m=0.5d0*(rho_m + rho_c)
         call rz_lim_on_rho_geom_theta(rho_c_m,
     &           theta_pol,r_mag,z_mag,z,r)
         del_psi_c_m=fpsi(r,z)-psi

         if((dabs(del_psi_c_m).lt.dabs(del_psi_m)).and.
     &      (dabs(del_psi_c_m).lt.dabs(del_psi_c)))then
c-----------point at psi,theta_pol is between rho_m and rho_c
c           poloidal flux has maximal value betweeen rho_m and rho_c

c            write(*,*)'calc_r_z_psi_theta_binary 4'

            rho_p=rho_c
            del_psi_p= del_psi_c
            rho_c=rho_c_m
            del_psi_c=del_psi_c_m
            goto 30
         endif

         rho_c_p=0.5d0*(rho_c + rho_p)
         call rz_lim_on_rho_geom_theta(rho_c_p,
     &           theta_pol,r_mag,z_mag,z,r)
         del_psi_c_p=fpsi(r,z)-psi

         if((dabs(del_psi_c_p).lt.dabs(del_psi_p)).and.
     &      (dabs(del_psi_c_p).lt.dabs(del_psi_c)))then
c-----------point at psi,theta_pol is between rho_c and rho_p
c           poloidal flux has maximal value betweeen rho_c and rho_p

c            write(*,*)'calc_r_z_psi_theta_binary 5'

            rho_m=rho_c
            del_psi_m=del_psi_c
            rho_c=rho_c_p
            del_psi_c=del_psi_c_p
            goto 30
         endif

c--------point at psi,theta_pol is between rho_c_m and rho_c_p
c        poloidal flux has maximal value betweeen rho_c_m and rho_c_p

c         write(*,*)'calc_r_z_psi_theta_binary 6'

         rho_m=rho_c_m
         del_psi_m= del_psi_c_m
         rho_p=rho_c_p
         del_psi_p=del_psi_c_p
         goto 30
      endif

      rho_m=rho_c 
      rho_c=rho_p 
          
      goto 15

 20   continue
      
c      write(*,*)'theta_pol',theta_pol
c      write(*,*)'z,r',z,r
c      write(*,*)'psi,fpsi(r,z',psi,fpsi(r,z)
c      write(*,*)'psi-fpsi(r,z',psi-fpsi(r,z)
     

      return
      end



       
c_for test only
      subroutine calc_r_z_psi_theta_binary_test(psi,
     &z_mag,r_mag,r,z,theta_pol)

c-----------------------------------------------------
c     calculate z,r on (psi,theta_pol)
c     It will create: r,z
c     at 0=<theta_pol<2*pi
c
c     Binary method is used. 

      implicit none

c-----input
      real*8 psi,        !poloidal flux
     &z_mag,r_mag,        !magnetic axis coordinates
     & theta_pol
   
c-----output
      real*8
     & r,z 
 
c-----external
      real*8 fpsi,max_rho_geom
      
c-----local
      integer i

      real*8 
     &step_r,r_l,r_r,r_c,psi1,psi2,
     &epspsi,                                                   !accuray
     &pi,rho_max_theta,rho_l,rho_r,step_rho,step_theta_pol,
     &tr,tl,t,rho_c_m,rho_c_p,rho_c,rho_m,rho_p,
     &del_psi_m,del_psi_c,del_psi_p,del_psi_c_p,del_psi_c_m,
     &rho_geom         
      integer n_rho_steps

c----------------------------------------------
      n_rho_steps=100  !number of steps in rho
      write(*,*)'calc_r_z_psi_theta_binary_t psi,z_mag,r_mag,theta_pol',
     &psi,z_mag,r_mag,theta_pol
      write(*,*)'calc_r_z_psi_theta_binary_test n_rho_steps',n_rho_steps

      pi=4.d0*datan(1.d0)

      epspsi=1.d-6        !set accuracy

c-----------------------------------------------------------
c     calculate the point  r_0,z_0 at psi,theta_pol  
c        r=r_mag+rho_geom*cos(theta_pol)
c        z=z_mag+rho_geom*sin(theta_pol)
c
c     using binary solver
c        
c-----------------------------------------------------------      

      rho_max_theta=max_rho_geom(theta_pol)
 
      write(*,*)'rho_max_theta, step_rho',rho_max_theta, step_rho

      step_rho=rho_max_theta/dble(n_rho_steps)

      rho_m=0.d0
      rho_c=rho_m+step_rho 

 15   continue
       write(*,*)'before rz_lim_on_rho_geom_theta'
      call rz_lim_on_rho_geom_theta(rho_m,
     &        theta_pol,r_mag,z_mag,z,r)
      del_psi_m=fpsi(r,z)-psi

      write(*,*)'rho_m,z,r,del_psi_m',rho_m,z,r,del_psi_m    

      call rz_lim_on_rho_geom_theta(rho_c,
     &        theta_pol,r_mag,z_mag,z,r)
         del_psi_c=fpsi(r,z)-psi

      write(*,*)'rho_c,z,r,del_psi_c',rho_c,z,r,del_psi_c

      if (del_psi_m*del_psi_c.lt.0d0) then
c--------point at psi, theta_pol is between rho_m and rho_c
         tl=rho_m
         tr=rho_c            
         do while ((tr-tl).gt.epspsi)

            t=tl+(tr-tl)*0.5d0
            call rz_lim_on_rho_geom_theta(t,
     &         theta_pol,r_mag,z_mag,z,r)
            psi1=fpsi(r,z)-psi

            call rz_lim_on_rho_geom_theta(tr,
     &         theta_pol,r_mag,z_mag,z,r)
            psi2=fpsi(r,z)-psi
                            
            if ((psi1*psi2).gt.0) then
               tr=t
            else
               tl=t
            end if
         end do


         write(*,*)'bin:theta_pol',theta_pol
         write(*,*)'r,z', r,z
         write(*,*)'fpsi(r,z)-psi',
     &              fpsi(r,z)-psi
         goto 20

      endif

      rho_p=rho_c+step_rho
           
      call rz_lim_on_rho_geom_theta(rho_p,
     &        theta_pol,r_mag,z_mag,z,r)
      del_psi_p=fpsi(r,z)-psi

      write(*,*)'rho_p,z,r,del_psi_p',rho_p,z,r,del_psi_p
      if (del_psi_c*del_psi_p.lt.0d0) then
c--------point at psi,theta_pol is between rho_c and rho_p
         tl=rho_c
         tr=rho_p            
         do while ((tr-tl).gt.epspsi)

            t=tl+(tr-tl)*0.5d0
            call rz_lim_on_rho_geom_theta(t,
     &         theta_pol,r_mag,z_mag,z,r)
            psi1=fpsi(r,z)-psi

            call rz_lim_on_rho_geom_theta(tr,
     &         theta_pol,r_mag,z_mag,z,r)
            psi2=fpsi(r,z)-psi
                            
            if ((psi1*psi2).gt.0) then
               tr=t
            else
               tl=t
            end if
         end do

         write(*,*)'bin:theta_pol',theta_pol
         write(*,*)'r,z', r,z
         write(*,*)'fpsi(r,z)-psi',
     &              fpsi(r,z)-psi
         goto 20

      endif

      if((dabs(del_psi_c).lt.dabs(del_psi_m)).and.
     &   (dabs(del_psi_c).lt.dabs(del_psi_p)))then

c--------point at psi,theta_pol is between rho_m and rho_p
c        poloidal flux has maximal value betweeen rho_m and rho_p
c------------------------------------------------------------------------
 30      continue

         if ((rho_p - rho_m).lt.epspsi) then
c-----------point is found
            call rz_lim_on_rho_geom_theta(rho_c,
     &              theta_pol,r_mag,z_mag,z,r)
                           
            write(*,*)'max:theta_pol',theta_pol
            write(*,*)'r,z', r,z
            write(*,*)'fpsi(r,z)-psi',
     &                 fpsi(r,z)-psi
            goto 20
         endif

           
         rho_c_m=0.5d0*(rho_m + rho_c)
         call rz_lim_on_rho_geom_theta(rho_c_m,
     &           theta_pol,r_mag,z_mag,z,r)
         del_psi_c_m=fpsi(r,z)-psi

         if((dabs(del_psi_c_m).lt.dabs(del_psi_m)).and.
     &      (dabs(del_psi_c_m).lt.dabs(del_psi_c)))then
c-----------point at psi,theta_pol is between rho_m and rho_c
c           poloidal flux has maximal value betweeen rho_m and rho_c
            rho_p=rho_c
            del_psi_p= del_psi_c
            rho_c=rho_c_m
            del_psi_c=del_psi_c_m
            goto 30
         endif

         rho_c_p=0.5d0*(rho_c + rho_p)
         call rz_lim_on_rho_geom_theta(rho_c_p,
     &           theta_pol,r_mag,z_mag,z,r)
         del_psi_c_p=fpsi(r,z)-psi

         if((dabs(del_psi_c_p).lt.dabs(del_psi_p)).and.
     &      (dabs(del_psi_c_p).lt.dabs(del_psi_c)))then
c-----------point at psi,theta_pol is between rho_c and rho_p
c           poloidal flux has maximal value betweeen rho_c and rho_p
            rho_m=rho_c
            del_psi_m=del_psi_c
            rho_c=rho_c_p
            del_psi_c=del_psi_c_p
            goto 30
         endif

c--------point at psi,theta_pol is between rho_c_m and rho_c_p
c        poloidal flux has maximal value betweeen rho_c_m and rho_c_p
         rho_m=rho_c_m
         del_psi_m= del_psi_c_m
         rho_p=rho_c_p
         del_psi_p=del_psi_c_p
         goto 30
      endif

      rho_m=rho_c 
      rho_c=rho_p 
          
      goto 15

 20   continue
     
      write(*,*)'theta_pol',theta_pol
      write(*,*)'z,r',z,r
      write(*,*)'psi,fpsi(r,z',psi,fpsi(r,z)
      write(*,*)'psi-fpsi(r,z',psi-fpsi(r,z)
     

      return
      end



       
c_for test only

       
