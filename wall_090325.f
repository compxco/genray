      subroutine wall_limiter_reflection(z_n,r_n,phi_n)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'fourb.i'

c-----input
      real*8 z_n,r_n,phi_n !ray coordinates

c-----locals
      real*8
     &z_old,r_old,phi_old, !ray coordinates at previous point
     &thetapol_n,          !poloidal angle at given (z_n,r_n) point
     &x_l,y_l,z_b,r_b,l_z_b_r_b

      integer i,i_pol_wall,i_pol_limiter,
     &i_pol_wall_count,
     &i_pol_wall_ar(10)
c-----externals
      real*8 thetapol

      save z_old,r_old,phi_old

      if(rho.le.1.d0) goto 10 !nohting to do

      x_l=r_n-xma 
      y_l=z_n-yma
c
      call theta_xy(x_l,y_l,thetapol_n) !poloidal angle at given (z_n,r_n) point

c      thetapol_n=thetapol(z_n,r_n) !poloidal angle at given (z_n,r_n) point

      i_pol_wall_count=0
c------------------------------------------------------------------------
c     using the ray poloidal angle thetapol_n determine 
c     numbers of wall points i_pol_wall_ar(1:i_pol_wall_count):
c
c     thetapol_wall(i_pol_wal-1) =< thetapol_n < thetapol_wall(i_pol_wal)
c
c-------------------------------------------------------------------------
      do i=2,n_wall
         if ((thetapol_n.lt.thetapol_wall(i)).and.
     &       (thetapol_n.ge.thetapol_wall(i-1))) then
             i_pol_wall_count=i_pol_wall_count+1
             i_pol_wall=i
             i_pol_wall_ar(i)=i
c             goto 20 
         endif
      enddo
 20   continue

c------------------------------------------------------------------------
c     using the ray poloidal angle thetapol_n determine 
c     number  of limiter point i_pol_limiter:
c
c     thetapol_limiter(i_pol_limiter-1)=< thetapol_n < thetapol_limirer(i_pol_limiter)
c
c-------------------------------------------------------------------------
      do i=2,n_limiter
         if ((thetapol_n.lt.thetapol_limiter(i)).and.
     &       (thetapol_n.ge.thetapol_limiter(i-1))) then
             i_pol_limiter=i
             goto 30 
         endif
      enddo

 30   continue

 10   continue

      return
      end

      subroutine wall_limiter_theta_pol_rho
c-----calculate poloidal angles [at radian] and small radius
c     of wall and limiter points
c     thetapol_wall(i=1,..,n_wall)
c     thetapol_limiter(i=1,...,n_limiter(j),j=1,max_limiters)
c     rho_wall(i=1,..,n_wall)
c     rho_limiter(i=1,..,n_limiter(j),j=1,max_limiters))
      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'fourb.i'

c-----input data are from common /one_nml/:
c     n_wal,n_limiter,
c     r_wall,z_wall,r_limiter,z_limiter

c-----output data
c     thetapol_wall(n_wall),thetapol_limiter(n_limiter)
c     rho_wall(n_wall),rho_limiter(n_limiter)
c     will be in common /fourb/

c-----locals
      integer i
      real*8 x_l,y_l,l_zr,l_z_b_r_b,z_b,r_b
c-----externals
      real*8 thetapol

      do i=1,n_wall 
         x_l=r_wall(i)-xma 
         y_l=z_wall(i)-yma
c--------calculate the wall poloidal angle: thetapol_wall
         call theta_xy(x_l,y_l,thetapol_wall(i))
c        thetapol_wall(i)=thetapol(z_wall,r_wall)

c--------calculate coordinates x_b,z_b at the LCFS surface psi=psilim
c        and at the poloidal angle thetapol_wall
         call zr_psith(psilim,thetapol_wall(i),z_b,r_b)

         l_zr=dsqrt(x_l**2+y_l**2)                   !distance between the magnetic axis 
                                                     !and the wall point
         l_z_b_r_b=dsqrt((r_b-xma)**2+(z_b-yma)**2)  !distance between the magnetic axis and
                                                     !LCFS at thetapol_wall(i)
c--------calculate the wall small radius
         rho_wall(i) = l_zr / l_z_b_r_b

      enddo
  
      do j=1,max_limiters
        do i=1,n_limiter
           x_l=r_limiter(i,j)-xma 
           y_l=z_limiter(i,j)-yma
c----------calculate limiters poloidal angles: thetapol_limiter(i)
           call theta_xy(x_l,y_l,thetapol_limiter(i,j))
c          thetapol_limiter(i)=thetapol(z_limiter,r_limiter)
c----------calculate coordinates x_b,z_b at the LCFS surface psi=psilim
c          and the poloidal angle thetapol_wall

           call zr_psith(psilim,thetapol_limiter(i,j),z_b,r_b)

           l_zr=dsqrt(x_l**2+y_l**2)                   !distance between the magnetic axis 
                                                       !and the limite point
           l_z_b_r_b=dsqrt((r_b-xma)**2+(z_b-yma)**2)  !distance between the magnetic axis and
                                                       !LCFS at thetapol_wall(i)
           rho_limiter(i,j) = l_zr / l_z_b_r_b
        enddo
      enddo

      return
      end


      subroutine ray_wall_limiter_intersection(z_old,r_old,phi_old,
     &z_n,r_n,phi_n,i,i_pol_wall_count,i_pol_wall_ar)
     
c-----check if the ray element [(z_old,r_old,phi_old).(z_n,r_n,phi_n)]
c     intersects wall or limiter

      implicit none

c-----input
      real*8
     &z_old,r_old,phi_old, !coordinates old ray point 
     &z_n,r_n,phi_n         !coordinates new ray point

      integer
     &i_pol_wall_count, ! the number of possible ray-wall intersections
     &i_pol_wall_ar(*)  ! numbers of poloidal angles where the ray 
                        ! can intersect the wall


c-----locals 
      integer i

c      z_p
c      do i=1,i_pol_wall_count
        
c      enddo


      return
      end

      subroutine ray_boundary_intersection(z_old,r_old,phi_old,
     &z_n,r_n,phi_n,z_m,r_m,z_p,r_p,
     &i_intersection,r_i,z_i)
     
c-----check if the ray element [(z_old,r_old,phi_old).(z_n,r_n,phi_n)]
c     intersects with the straight line
c     connected points at the boundary (z_m,r_m),(z_p,r_p)


      implicit none

c-----input
      real*8
     &z_old,r_old,phi_old, !coordinates old ray point 
     &z_n,r_n,phi_n,       !coordinates new ray point
     &z_m,r_m,             !coordinates of points 'm' at the boundary
     &z_p,r_p              !coordinates of points 'p' at the boundary
c-----output
      integer
     &i_intersection       !=0 no intersection
                           !=1 it has intersection
      real*8
     &r_i,z_i,phi_i        !intersection point coordinates    
c---------------------------------------     
c     The equation of the srtaihgt line connecting boundary points
c
c     r=r_m + (r_p-r_m)*q, z=z_m + (z_p-z_m)*q                    eq.(1)
c
c     The equation of the srtaihgt line connecting ray points
c
c     x=x_old+(x_n-x_old)*t, y=y_old+(y_n-y_old)*t,               eq.(2)
c     z=z_old+(x_n-z_old)*t 
c
c     In cylindrical coordinates x=rsos(phi), y=r*sin(phi)
c     this equation connecting ray points has the form
c
c     r^2=r_old^2 + 2*r_old*t*[r_n*cos(pfi_n-phi_old)-r_old] +   eq.(3)
c     + t^2*[r_n^2-2*r_old*r_n*cos(phi_n-phi_old)+r_old^2 ] 
c
c     The intersection point (z_i,r_i) equations are
c
c     r_i^2=r_old^2 + 2*r_old*t_i*[r_n*cos(pfi_n-phi_old)-r_old]+
c         +t_i^2*[r_n^2-2*r_old*r_n*cos(phi_n-phi_old)+r_old^2]=  eq.(4)
c         =[r_m + (r_p-r_m)*q_i]^2  
c
c     z_i=z_old+(x_n-z_old)*t_i=z_m + (z_p-z_m)*q_i               eq.(5)     
c
c     at (z_p.ne.z_m)
c        for non-vertical boundary eq.(5) gives
c        q_i=[(z_old-z_m) + (z_n-z_old)*t_i]/(z_p-z_m)            eq.(6)
c
c        Using eq.(6) in eq. (4) gives Eq. for t_i
c
c       a*t_i^2 +2*b*t_i+c=0                                      eq.(7)
c
c     a=[r_n^2-2*r_old*r_n*cos(phi_n-phi_old)+r_old^2]-
c       -[(r_p-r_m)*(z_n-z_old)/(z_p-z_m)]^2                      eq.(8a)
c
c     b=r_old*{r_n*cos(phoi_n-phi_old)-r_old]-
c       -[r_m+(r_p-r_m)(z_old-z_m)/(z_p_z_m)]*                    eq.(8b)
c        *[(r_p-r_m)(z_n-z_old)/(z_p-z_m)]
c    
c     c=r_old^2-[r_m+(r_p-r_m)(z_old-z_m)/(z_p-z_m)]^2            eq.(8c)
c
c-----locals 
      real*8 cos_phi,a,b,c,det,t_i,t_i_p,t_i_m,q_i,
     &l_boundary,l_boundary_p,l_boundary_m,pi

      integer n_phi ![phi_n/(2*pi)]

      pi=4.d0*datan(1.d0)
      n_phi=phi_n/(2.d0*pi)

      i_intersection=0

      cos_phi=dcos(phi_n-phi_old)
      if (z_p.ne.z_m) then 
c--------Using eq.(6) in eq. (4) gives Eq. for t_i
         a=(r_n**2-2.d0*r_old*r_n*cos_phi+r_old**2)+
     &     ((r_p-r_m)*(z_n-z_old)/(z_p-z_m))**2
         
         b=r_old*(r_n*cos_phi-r_old)-
     &     (r_m+(r_p-r_m)*(z_old-z_m)/(z_p-z_m))*
     &     ((r_p-r_m)*(z_n-z_old)/(z_p-z_m))

         c=r_old**2-(r_m+(r_p-r_m)*(z_old-z_m)/(z_p-z_m))**2
         
         det=b**2-a*c


         if(det.lt.0.d0) then
c----------no roots
           goto10
         else
           if (det.eq.0.d0) then
c------------one root
             t_i=b/a
             q_i=((z_old-z_m)+(z_n-z_old)*t_i)/(z_p-z_m)
             z_i=z_old+(z_n-z_old)*q_i
             r_i=r_m+(r_p-r_m)*q_i
             l_boundary=dsqrt((z_p-z_p)**2+(r_p-r_m)**2)
             l_boundary_p=dsqrt((z_p-z_i)**2+(r_p-r_i)**2)
             l_boundary_m=dsqrt((z_m-z_i)**2+(r_m-r_i)**2)
             if((l_boundary_p.lt.l_boundary).and.
     &          (l_boundary_m.lt.l_boundary)) then
                 i_intersection=1
             endif
           else
c------------two roots
             t_i_p=b+dsqrt(det)/a
             t_i_m=b-dsqrt(det)/a

c------------check first root
             t_i=t_i_p
             q_i=((z_old-z_m)+(z_n-z_old)*t_i)/(z_p-z_m)
             z_i=z_old+(z_n-z_old)*q_i
             r_i=r_m+(r_p-r_m)*q_i
             l_boundary=dsqrt((z_p-z_p)**2+(r_p-r_m)**2)
             l_boundary_p=dsqrt((z_p-z_i)**2+(r_p-r_i)**2)
             l_boundary_m=dsqrt((z_m-z_i)**2+(r_m-r_i)**2)
             if((l_boundary_p.lt.l_boundary).and.
     &          (l_boundary_m.lt.l_boundary)) then
                 i_intersection=1
                 goto 10
             endif
c------------check secon root
             t_i=t_i_m
             q_i=((z_old-z_m)+(z_n-z_old)*t_i)/(z_p-z_m)
             z_i=z_old+(z_n-z_old)*q_i
             r_i=r_m+(r_p-r_m)*q_i
             l_boundary=dsqrt((z_p-z_p)**2+(r_p-r_m)**2)
             l_boundary_p=dsqrt((z_p-z_i)**2+(r_p-r_i)**2)
             l_boundary_m=dsqrt((z_m-z_i)**2+(r_m-r_i)**2)
             if((l_boundary_p.lt.l_boundary).and.
     &          (l_boundary_m.lt.l_boundary)) then
                 i_intersection=1
                 goto 10
             endif
           
           endif
         endif
      else  
c-------z_p.ne.z_m the verticle boundary
c       t_i*(z_n-z_old)=(z_m-z_0)        
        if (z_n.ne.z_old) then
           t_i=(z_m-z_old)/(z_n-z_old)
           z_i=z_m
           r_i=dsqrt(r_old**2+2.d0*r_old*t_i*(r_n*cos_phi-r_old)+
     &               t_i**2*(r_n**2-2.d0*r_old*r_n*cos_phi+r_old**2))
        else
           !the ray is parallel to verrtical boundary
           !no intersection
        endif
      endif
   
 10   continue
      return
      end

     
      

      
      subroutine wall_limiter_rho(theta,phi,n_rho_wall,rho_ar)
c--------------------------------------------------------------
c     calculates several n_rho_wall small radii values rho_ar(1:n_rho_wall) 
c     of the chamber boundary (wall plus limiter) at given
c     poloidal theta and and toroidal phi angles [radians]
c  
c     the poloidal angle is measured counter-clockwise
c                        from the outer part the equatorial plane  
c
c     the toroidal angle is measured counter-clockwise
c                  around the vertical Z axis
c
c---------------------------------------------------------------
      implicit none
    
      include 'param.i'
c      include 'three.i'
      include 'fourb.i'

c-----input
      real*8 theta,phi !poloidal and toroidal angles [radians]

c-----output
      integer n_rho_wall !number of small radius (chamber points) at given
                    !poloidal and toroidal angles 
      real*8 rho_ar(n_rho_wall_a) !wall plus limiter small radii

c-----locals
      integer i
      return
      end
