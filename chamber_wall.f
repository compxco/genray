      subroutine wall_limiter_reflection_point(z_n,r_n,phi_n,
     &cnz_n,cnr_n,cm_n,i_reflection,
     &cnz_refl,cnr_refl,cm_refl,z_refl,r_refl,phi_refl)
c---------------------------------------------------------------------------
c     calculate wall or limiter reflection point coordinates:
c     z_refl,r_refl,phi_refl,
c     and the refractive index after reflection
c     cnz_refl,cnr_refl,cm_refl
c     output: 
c     index i_reflection=1 if the given point z_n,r_n,phi_n
c                          and the previous ray 'old' point 
c                          are at the different sides of the wall.
c                       =0 z_n,r_n,phi_n point does not give the reflection
c                          
c     z_refl,r_refl,phi_refl  - reflection point coordinates:
c     cnz_old,cnr_old,cm_old, - refractive index coordinates at
c                               the privious ray point
c 
c--------------------------------------------------------------------
c     the wall element M_p(z_p,r_p),M_m(z_m,r_m) near the reflection 
c     point is a straight line
c
c     The unit vector parallel to the wall element lw(lw_z.lw_r):
c
c     lw_z=(z_p-z_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     lw_r=(r_p-r_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     
c     The unit vector perpendicular to the wall element Nw(nw_z,nw_r)
c     is the solution of the following equation: 
c     lw_z*nm_z+lw_r*nm_r=0
c
c     If (lw_z.ne.0) then nm_z=-(lw_r/lw_z)*nm_r
c     else (lw_z.ne.0) then nm_r=-(lw_z/lw_r)*nm_z
c
c     Unit Nw: nm_z/sqrt(mn_z**2+mnr**2),nm_r/sqrt(mn_z**2+mnr**2)
c
c     The refractive vector perpendicular to the wall element :
c     N_perp^=(cnz_perp,cnr_perp)
c
c     (N^.Nw^)=cnz*nm_z+cnr*nm_r
c
c     cnr_perp= (N^.Nw^)*nm_r,  cnr_perp= (N^.Nw^)*nm_r
c
c     The reflected refactive index N_refl^=N^-2*N_perp^
c-------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'fourb.i'

c-----input
      real*8 z_n,r_n,phi_n, !ray coordinates
     &cnz_n,cnr_n,cm_n      !ray refractive index coordinates
c-----output
      integer i_reflection !=1 at reflection point
                           !=0 no reflection
      real*8 z_refl,r_refl,phi_refl, !reflection point coordinates
     &cnz_refl,cnr_refl,cm_refl      !refractive index in the ray point
                                     !after reflection
c-----locals
      real*8
     &z_old,r_old,phi_old, !ray coordinates at previous point
     &thetapol_n,          !poloidal angle at given (z_n,r_n) point
     &x_l,y_l,z_b,r_b,l_zr,l_z_b_r_b,rho_n,rho_old,
     &rho_wall_ar(n_rho_wall_a), !wall small radii at given thetapol_n
     &dif_rho_min,               !min(abs(rho_n-rho_wall_arr))   
     &dif_rho,                   !rho_n-rho_wall_ar at minimal
                                 ! abs(rho_n-rho_wall_arr))   
     &dif_rho_old,
     &cnz_old,cnr_old,cm_old,        !refractive index at the previous point
cfor test
     &u_l(6),hamilt                                     
cend for test 
      integer i,
     &n_rho_wall,                 !number of wall small radii at 
                                  !given thetapol_n
     &n_rho,                      !number of point at which wall small radius
                                  !is close to ray small radius 
                                  !at given polidal angle in reflection point
     &i_wall_ar_m(n_rho_wall_a),  !numbers of wall points around the straight 
     &i_wall_ar_p(n_rho_wall_a)   !line from the magnetic directed at the 
                                  !poloidal angle theta
     
      logical first

      data first/.true./
      save z_old,r_old,phi_old,rho_old,dif_rho_old,
     &cnz_old,cnr_old,cm_old 

      write(*,*)'wall_limiter_reflection_point rho',rho
      i_reflection=0

      if(rho.le.1.d0) goto 10 !nohting to do

cSAP090403
c      write(*,*)'n_wall',n_wall
      if (n_wall.eq.0) then
c--------There ais no wall 
         goto 10
      endif


      x_l=r_n-xma 
      y_l=z_n-yma
c
      call theta_xy(x_l,y_l,thetapol_n) !poloidal angle at given (z_n,r_n) point

c      write(*,*)'z_n,r_n,thetapol_n', z_n,r_n,thetapol_n

c
c-----calculate coordinates x_b,z_b at the LCFS surface psi=psilim
c     and the poloidal angle thetapol_wall
      call zr_psith(psilim,thetapol_n,z_b,r_b)

      l_zr=dsqrt(x_l**2+y_l**2)                   !distance between the magnetic axis 
                                                  !and the ray point
      l_z_b_r_b=dsqrt((r_b-xma)**2+(z_b-yma)**2)  !distance between the magnetic axis and
                                                  !LCFS at thetapol_wall_n 
      rho_n = l_zr / l_z_b_r_b                    !small radius at given z_n,r_n

c      write(*,*)'rho_n',rho_n
c-----------------------------------------------------------------------------
c     calculates wall small radius at the given poloidal angle theta [radian]
c     It is proposed that between two wall neighbor points
c     wall small radius is a linear function on poloidal angle
c-----------------------------------------------------------------------------
      call rho_wall_at_theta(z_n,r_n,rho_wall_ar,n_rho_wall,
     &i_wall_ar_m,i_wall_ar_p)

      if (n_rho_wall.eq.0) then
c--------There is no ray-wall intersection 
         goto 10
      endif

c      write(*,*)'after rho_wall_at_theta_2 n_rho',n_rho_wall
c      write(*,*)'rho_wall_ar(i)',(rho_wall_ar(i),i=1,n_rho_wall)
c-------------------------------------------------------------------
c     choose small wall radius more closer to rho_n
c--------------------------------------------------------------------
c      write(*,*)'rho_n',rho_n
      n_rho=0

      dif_rho_min=1.d9
      do i=1,n_rho_wall
        if ( dif_rho_min .gt. dabs(rho_n-rho_wall_ar(i))) then
           dif_rho_min=dabs(rho_n-rho_wall_ar(i))
           dif_rho=rho_n-rho_wall_ar(i)
           n_rho=i
        endif 
      enddo
c      write(*,*)'n_rho',n_rho

c      write(*,*)'first',first
       
      if (first) then
          z_old=z_n
          r_old=r_n
          phi_old=phi_n
          cnz_old=cnz_n
          cnr_old=cnr_n 
          cm_old= cm_n 
          dif_rho_old=-1.d0 
          first=.false.
      endif

c      write(*,*)'dif_rho_old,dif_rho',dif_rho_old,dif_rho
c      write(*,*)'cnz_n,cnr_n,cm_n',cnz_n,cnr_n,cm_n

c      write(*,*)'cnz_old,cnr_old,cm_old',
c     &           cnz_old,cnr_old,cm_old

      if (dif_rho_old*dif_rho.lt.0.d0) then
c--------ray-wall intersection point is between old and new points
c        The simple solution is to use the old point as a reflection point. 
         i_reflection=1
         z_refl=z_old
         r_refl=r_old
         phi_refl=phi_old  
         first=.true.
c         first=.false.

c         write(*,*)'i_wall_ar_m(n_rho),i_wall_ar_p(n_rho)',
c     &              i_wall_ar_m(n_rho),i_wall_ar_p(n_rho)
c         write(*,*)'r_wall(i_wall_ar_m(n_rho))',
c     &              r_wall(i_wall_ar_m(n_rho))
c         write(*,*)'z_wall(i_wall_ar_m(n_rho))',
c     &              z_wall(i_wall_ar_m(n_rho))

c         write(*,*)'r_wall(i_wall_ar_p(n_rho))',
c     &              r_wall(i_wall_ar_p(n_rho))
c         write(*,*)'z_wall(i_wall_ar_p(n_rho))',
c     &              z_wall(i_wall_ar_p(n_rho))

c         write(*,*)'r_old,z_old',r_old,z_old
c--------check hamiltonian value before reflection
         u_l(1)=z_old
         u_l(2)=r_old 
         u_l(3)=phi_old
         u_l(4)=cnz_old
         u_l(5)=cnr_old 
         u_l(6)=cm_old

         call callc_eps_ham(u_l,hamilt)
         write(*,*)'before reflection: hamilt=',hamilt
         write(*,*)'before reflection: h_l=',u_l
c--------calculate refractive index after reflection
         cm_refl=cm_old         
         call wall_limiter_N_reflection(cnz_old,cnr_old,n_rho,
     &       i_wall_ar_p,i_wall_ar_m,cnz_refl,cnr_refl)

         write(*,*)'wall_limiter_reflection_poi z_refl,r_refl,phi_refl'
     &              ,z_refl,r_refl,phi_refl 
         write(*,*)'cnz_refl,cnr_refl',cnz_refl,cnr_refl

c         z_old=z_n
c         r_old=r_n     
c         phi_old=phi_n
c         rho_old=rho_n
c         dif_rho_old=dif_rho 

         cnz_old=cnz_refl
         cnr_old=cnr_refl
c         cm_old=cm_n

c---------check hamiltonian value after reflection
          u_l(1)=z_refl
          u_l(2)=r_refl 
          u_l(3)=phi_refl
          u_l(4)=cnz_refl
          u_l(5)=cnr_refl
          u_l(6)=cm_refl

          call callc_eps_ham(u_l,hamilt)
          write(*,*)'after reflection: hamilt=',hamilt
          write(*,*)'after reflection: h_l=',u_l

          goto 10
      endif

c-----set old point equal to given point
      z_old=z_n
      r_old=r_n     
      phi_old=phi_n
      rho_old=rho_n
      dif_rho_old=dif_rho 

      cnz_old=cnz_n
      cnr_old=cnr_n
      cm_old=cm_n


c      write(*,*)'cnz_old,cnr_old,cm_old',
c     &           cnz_old,cnr_old,cm_old

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
c  
c     These arrays will be in common /fourb/ in file fourb.i
c
      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'fourb.i'

c-----input data are from common /one_nml/:
c     n_wal,n_limiter,
c     r_wall,z_wall,r_limiter,z_limiter

c-----output data
c     thetapol_wall(n_wall),thetapol_limiter(n_limiter_a,max_limiters_a),
c     rho_wall(n_wall),rho_limiter(n_limiter_a,max_limiters_a),
c     will be in common /fourb/

c-----locals
      integer i,i_theta_pol_min,n
      real*8 x_l,y_l,l_zr,l_z_b_r_b,z_b,r_b,
     &dif_m,dif_p,theta_pol_min,
     &temp(n_wall_a)       !work array
c-----externals
      real*8 thetapol

      write(*,*)'in wall_limiter_theta_pol_rho'

      do i=1,n_wall 
         x_l=r_wall(i)-xma 
         y_l=z_wall(i)-yma
c--------calculate the wall poloidal angle: thetapol_wall
         call theta_xy(x_l,y_l,thetapol_wall(i)) !0 =< theta_l < 2*pi

c        thetapol_wall(i)=thetapol(z_wall,r_wall)

c--------calculate coordinates x_b,z_b at the LCFS surface psi=psilim
c        and the poloidal angle thetapol_wall
         call zr_psith(psilim,thetapol_wall(i),z_b,r_b)

         l_zr=dsqrt(x_l**2+y_l**2)                   !distance between 
                                                     ! the magnetic axis 
                                                     !and the wall point
         l_z_b_r_b=dsqrt((r_b-xma)**2+(z_b-yma)**2)  !distance between
                                                     ! the magnetic axis and
                                                     !LCFS at thetapol_wall(i)
c--------calculate the wall small radius
         rho_wall(i) = l_zr / l_z_b_r_b
c         write(*,*)'i,r_wall(i),z_wall(i),xma,yma', 
c     &              i,r_wall(i),z_wall(i),xma,yma
c         write(*,*)'z_b,r_b, l_zr, l_z_b_r_b,rho_wall(i)',
c     &              z_b,r_b, l_zr, l_z_b_r_b,rho_wall(i)
          
      enddo

c-----renumeration of wall points accoding to poloidal angles
c     the first wall point will have the minimal poloidal angle 

      theta_pol_min=2.d0*pi
      do i=1,n_wall
        if (theta_pol_min.gt.thetapol_wall(i)) then
          theta_pol_min=thetapol_wall(i)
          i_theta_pol_min=i
        endif
      enddo

ctest081105
      if (i_theta_pol_min.eq.1) goto 20
cendtest

      dif_p=dabs(thetapol_wall(i_theta_pol_min+1)-
     &           thetapol_wall(i_theta_pol_min))

      dif_m=dabs(thetapol_wall(i_theta_pol_min-1)-
     &           thetapol_wall(i_theta_pol_min))

      call  renumeration_wall(n_wall,i_theta_pol_min,r_wall,temp,
     &dif_p,dif_m)
      call  renumeration_wall(n_wall,i_theta_pol_min,z_wall,temp,
     &dif_p,dif_m)
      call  renumeration_wall(n_wall,i_theta_pol_min,rho_wall,
     &temp,dif_p,dif_m)
c      write(*,*)'after renumeration rho_wall',rho_wall 
      call  renumeration_wall(n_wall,i_theta_pol_min,thetapol_wall,
     &temp,dif_p,dif_m)
 
      r_wall(n_wall)=r_wall(1)
      z_wall(n_wall)=z_wall(1)
      rho_wall(n_wall)=rho_wall(1)
      thetapol_wall(n_wall)=thetapol_wall(1)

c      write(*,*)'n_wall',n_wall
c      do i=1,n_wall
c        write(*,*)'i,thetapol_wall(i),rho_wall(i),r_wall(i),z_wall(i)',
c     &             i,thetapol_wall(i),rho_wall(i),r_wall(i),z_wall(i)       
c      enddo
c      !stop 'wall_limiter_theta_pol_rho'

 20   continue
c      goto 10
c------------------------------------------------------------------------
c     limiter points
c------------------------------------------------------------------------
c     calculte limiter points poloidal angles and small radii

c      write(*,*)'max_limiters',max_limiters
   
      do n=1, max_limiters

c         write(*,*)'n',n,'n_limiter(n)',n_limiter(n)

        do i=1,n_limiter(n)
          write(*,*)'i',i
          x_l=r_limiter(i,n)-xma 
          y_l=z_limiter(i,n)-yma

c          call theta_xy(x_l,y_l,thetapol_limiter(i,n)) !0 =< theta_pol < 2*pi

          thetapol_limiter(i,n)=                         !-pi=<theta_pol<pi 
     &                          thetapol(z_limiter(i,n),r_limiter(i,n))

c---------calculate coordinates x_b,z_b at the LCFS surface psi=psilim
c         and the poloidal angle thetapol_limiter

          call zr_psith(psilim,thetapol_limiter(i,n),z_b,r_b)

          l_zr=dsqrt(x_l**2+y_l**2)        !distance between the magnetic axis 
                                           !and the limiter point
          l_z_b_r_b=dsqrt((r_b-xma)**2+(z_b-yma)**2)  !distance between
                                                      !the magnetic axis and
                                                      !LCFS at thetapol_wall(i)
          rho_limiter(i,n) = l_zr / l_z_b_r_b
 
c          write(*,*)'i,n,thetapol_limiter(i,n)',
c     &               i,n,thetapol_limiter(i,n)

        enddo
      enddo 

      do n=1,max_limiters    
         theta_pol_min=2.d0*pi

c         write(*,*)'1 n=',n,'n_limiter(n)',n_limiter(n)

         do i=1,n_limiter(n)
c           write(*,*)'i',i
           if (theta_pol_min.gt.thetapol_limiter(i,n)) then

c              write(*,*)'1 theta_pol_min,thetapol_limiter(i,n)',
c     &                   theta_pol_min,thetapol_limiter(i,n)

              theta_pol_min=thetapol_limiter(i,n)
              i_theta_pol_min=i

c              write(*,*)'theta_pol_min, i_theta_pol_min',
c     &                   theta_pol_min, i_theta_pol_min
           endif
            
c           write(*,*)'i_theta_pol_min',i_theta_pol_min

        enddo !i

        dif_p=dabs(thetapol_limiter(i_theta_pol_min+1,n)-
     &             thetapol_limiter(i_theta_pol_min,n))

        dif_m=dabs(thetapol_limiter(i_theta_pol_min-1,n)-
     &             thetapol_limiter(i_theta_pol_min,n))

        call  renumeration_limiter(n_limiter(n),i_theta_pol_min,
     &  r_limiter(1,n),temp,dif_p,dif_m)
        call  renumeration_limiter(n_limiter(n),i_theta_pol_min,
     &  z_limiter(1,n),temp,dif_p,dif_m)
        call  renumeration_limiter(n_limiter(n),i_theta_pol_min,
     &  rho_limiter(1,n),temp,dif_p,dif_m)
        write(*,*)'after renumeration_limiter rho_wall',rho_wall 
        call  renumeration_limiter(n_limiter(n),i_theta_pol_min,
     &  thetapol_limiter(1,n),temp,dif_p,dif_m)
 
c-test---------------------------------------------------------
c        do i=1,n_limiter(n)
c           write(*,*)'i,r_limiter(i,n),z_limiter(i,n)',
c     &                i,r_limiter(i,n),z_limiter(i,n)
c           write(*,*)'rho_limiter(i,n),thetapol_limiter(i,n)',
c     &                rho_limiter(i,n),thetapol_limiter(i,n)
c        enddo
c-end test------------------------------------------------------

      enddo !n

 10   continue
c      !stop ' wall_limiter_theta_pol_rho'
      return
      end
  

      subroutine renumeration_wall(n_max,i_min,ar,temp,difp,difm)
c-----renumeration !it will set ar(n_max)=ar(1)
      implicit none
      
c-----input      
      integer n_max,  !the number of points in array ar(1:n_max)
     &i_min           !the number of point in the original arrray ="ar"
                      ! which will be renumerated to number =1
      real*8 ar(*),   !array which should be renumerated
     &temp(*)         !work array with dimension(tepm) >= dimension(ar) 
      real*8 difp,difm  !if difm<difp

c-----output
c     ar(*)             after renumeration 
                      
c-----locals
c      real*8,  dimension(1,n_max) :: temp
      integer i,k

      do i=1,n_max
        temp(i)=ar(i)
      enddo

      if(difm.lt.difp) then
         do i=i_min,1,-1
            k=i_min-i+1  !k=1,...,i_min
            ar(k)=temp(i)
         enddo

         do i=n_max-1,i_min+1,-1
            k=n_max-i+i_min !k=i_min+1,...,n_max-1
            ar(k)=temp(i)
         enddo
      else
c--------difm.ge.difp
         do i=i_min,n_max-1
            k=i-i_min+1     !k=1,...,(n_max-i_min)
            ar(k)=temp(i)
         enddo

         do i=1,i_min-1
            k=n_max-i_min+i !k=(n_max-i_min+1),...,n_max-1
            ar(k)=temp(i)
         enddo
        
      endif

      ar(n_max)=ar(1)

      return
      end      

      subroutine wall_limiter_N_reflection(cnz,cnr,n_rho,
     &i_wall_ar_p,i_wall_ar_m,cnzrefl,cnrrefl)
c---------------------------------------------------------------------------
c     calculate wall refactive index after reflection:
c     cnzrefl,cnrrefl
c--------------------------------------------------------------------
c     the wall element M_p(z_p,r_p),M_m(z_m,r_m) near the reflection 
c     point is a straight line
c
c     The unit vector parallel to the wall element lw(lw_z.lw_r):
c
c     lw_z=(z_p-z_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     lw_r=(r_p-r_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     
c     The unit vector perpendicular to the wall element Nw(nw_z,nw_r)
c     is the solution of the following equation: 
c     lw_z*nm_z+lw_r*nm_r=0
c
c     If (lw_z.ne.0) then nm_z=-(lw_r/lw_z)*nm_r
c     else (lw_z.ne.0) then nm_r=-(lw_z/lw_r)*nm_z
c
c     Unit Nw: nm_z/sqrt(mn_z**2+mnr**2),nm_r/sqrt(mn_z**2+mnr**2)
c
c     The refractive vector perpendicular to the wall element :
c     N_perp^=(cnz_perp,cnr_perp)
c
c     (N^.Nw^)=cnz*nm_z+cnr*nm_r
c
c     cnr_perp= (N^.Nw^)*nm_r,  cnr_perp= (N^.Nw^)*nm_r
c
c     The reflected refactive index N_refl^=N^-2*N_perp^

      implicit none
      include 'param.i'
      include 'one_nml.i'
      include 'fourb.i'
c-----input 
      integer 
     &n_rho,                      !number of point at which wall small radius
                                  !is close to ray small radius 
                                  !at given polidal angle in reflection point
     &i_wall_ar_m(n_rho_wall_a),  !numbers of wall points around the straight 
     &i_wall_ar_p(n_rho_wall_a)   !line from the magnetic axis directed at the 
                                  !poloidal angle theta    

      real*8 cnz,cnr !refractive index before reflection
c-----output
      real*8
     &cnzrefl,         ! refractive index after reflection            
     &cnrrefl          ! at wall element 
c-----locals
      real*8
     &z_p,r_p,z_m,r_m, !cordinates of points M_p,M_m
     &lw_z,lw_r,       !unit vector along wall element (M_m,M_p)
     &nw_z,nw_r,       !unit vector normal to wall element (M_m,M_p)
     &cnz_perp_w,      !The refractive vector perpendicular 
     &cnr_perp_w,      !to the wall element 
     &p
c-------------------------------------------------------------------
c     the wall element M_p(z_p,r_p),M_m(z_m,r_m) near the reflection 
c     point 
c------------------------------------------------------------------
      z_p=z_wall(i_wall_ar_p(n_rho))
      r_p=r_wall(i_wall_ar_p(n_rho))
      z_m=z_wall(i_wall_ar_m(n_rho))
      r_m=r_wall(i_wall_ar_m(n_rho))

c------------------------------------------------------------------
c     The unit vector parallel to the wall element lw(lw_z.lw_r):
c
c     lw_z=(z_p-z_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     lw_r=(r_p-r_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c------------------------------------------------------------------
      lw_z=z_p-z_m
      lw_r=r_p-r_m
      p=1.d0/dsqrt(lw_z**2+lw_r**2)
      lw_z=lw_z*p
      lw_r=lw_r*p
      write(*,*)'lw_z,lw_r',lw_z,lw_r
c-----------------------------------------------------------------------
c     The unit vector perpendicular to the wall element Nw(nw_z,nw_r)
c     is the solution of the following equation: 
c     lw_z*nw_z+lw_r*nw_r=0
c-----------------------------------------------------------------------
      
      if (lw_z.ne.0) then
         nw_r=1.d0
         nw_z=-(lw_r/lw_z)*nw_r
         p=1.d0/dsqrt(nw_z**2+nw_r**2)
         nw_z=  nw_z*p
         nw_r=  nw_r*p
      else
c--------lw_r.ne.0) 
         nw_z=1.d0
         nw_r=-(lw_z/lw_r)*nw_z
         p=1.d0/dsqrt(nw_z**2+nw_r**2)
         nw_z=  nw_z*p
         nw_r=  nw_r*p
      endif

      write(*,*)'nw_z,nw_r',nw_z,nw_r
      write(*,*)'lw_z*nw_z+lw_r*nc',lw_z*nw_z+lw_r*nw_r
c-----------------------------------------------------------------------
c     The refractive vector perpendicular to the wall element :
c     N_perp^=(cnz_perp_w,cnr_perp_w)
c
c     (N^.Nw^)=cnz*nw_z+cnr*nw_r
c
c     cnr_perp_w= (N^.Nw^)*nm_r,  cnr_perp_w= (N^.Nw^)*nm_r
c-----------------------------------------------------------------------
      p=cnz*nw_z+cnr*nw_r !(N^.Nw^)

      cnz_perp_w= nw_z*p
      cnr_perp_w= nw_r*p
      write(*,*)'(N^.Nw^)',p
c------------------------------------------------------------------------
c     The refactive index after reflection N_refl^=N^-2*N_perp^
c-----------------------------------------------------------------------
      cnzrefl=cnz-2.d0*cnz_perp_w
      write(*,*)'cnz,cnz_perp_w,cnzrefl',cnz,cnz_perp_w,cnzrefl
      cnrrefl=cnr-2.d0*cnr_perp_w
      write(*,*)'cnr,cnr_perp_w,cnrrefl',cnr,cnr_perp_w,cnrrefl
      return
      end


      subroutine renumeration_limiter(n_max,i_min,ar,temp,difp,difm)
c-----renumeration 
      implicit none
      
c-----input      
      integer n_max,  !the number of points in array ar(1:n_max)
     &i_min           !the number of point in the original arrray ="ar"
                      ! which will be renumerated to number =1
      real*8 ar(*),   !array which should be renumerated
     &temp(*)         !work array with dimension(tepm) >= dimension(ar) 
      real*8 difp,difm  !if difm<difp

c-----output
c     ar(*)             after renumeration 
                      
c-----locals
c      real*8,  dimension(1,n_max) :: temp
      integer i,k

      do i=1,n_max
        temp(i)=ar(i)
      enddo

      if(difm.lt.difp) then
         do i=i_min,1,-1
            k=i_min-i+1  !k=1,...,i_min
            ar(k)=temp(i)
         enddo

         do i=n_max,i_min+1,-1
            k=n_max-i+i_min+1 !k=i_min+1,...,n_max
            ar(k)=temp(i)
         enddo
      else
c--------difm.ge.difp
         do i=i_min,n_max
            k=i-i_min+1     !k=1,...,(n_max-i_min+1)
            ar(k)=temp(i)
         enddo

         do i=1,i_min-1
            k=n_max-i_min+i+1 !k=(n_max-i_min+2),...,n_max
            ar(k)=temp(i)
         enddo
      endif

      return
      end      

      subroutine rho_wall_at_theta(z_q,r_q,rho_wall_ar,n_rho_wall,
     &i_wall_ar_m,i_wall_ar_p)
c-----calculates wall small radius at poloidal angle theta [radian]
c     of the given ray point (z_q,r_q)
c     It is proposed that between two wall neighbor points
c     wall small radius is a linear function on poloidal angle
c
c     inside this subroutine 0=<theta=thetapol<2*pi

      implicit none
      include 'param.i'
      include 'one_nml.i'
      include 'three.i'
      include 'fourb.i'
c-----input
      real*8
     &z_q,r_q  !space coordinates of the ray point

c       theta !poloidal angle [radians]

c-----output
      integer  n_rho_wall,        !number of small radii at the given poloidal angle
     &i_wall_ar_m(n_rho_wall_a),  !numbers of wall points around the straight 
     &i_wall_ar_p(n_rho_wall_a)   !line from the magnetic directed at the 
                                  !poloidal angle theta
                                
      real*8 rho_wall_ar(n_rho_wall_a) ! small radii at given poloidal angle
    
c-----locals
      integer i,j
      real*8     
     &theta,                      !poloidal angle [radians] of the ray point(z_q,r_q)
     & pi,r_p,z_p,r_m,z_m,
c     &cos_theta,sin_theta,
     &a,b,a_q,b_q,
     &r_b,z_b,l_z_b_r_b,r_k,z_k,
     &x_l,y_l        

      write(*,*)'rho_wall_at_theta_2 z_q,r_q',z_q,r_q

      pi=4*datan(1.d0)


      x_l=r_q-xma 
      y_l=z_q-yma
      call theta_xy(x_l,y_l,theta)  !0=<theta<2*pi
      
      n_rho_wall=0
c------------------------------------------------------------------------
c     determine 
c     numbers of wall points:
c     i_wall_ar_m(1:n_rho__wall),i_wall_ar_p(1:n_rho__wall),
c
c     thetapol_wall(i-1) =< theta < thetapol_wall(i)
c     or
c     thetapol_wall(i) < theta =< thetapol_wall(i-1)
c     
c-------------------------------------------------------------------------
      do i=1,n_wall
c        write(*,*)'i,theta,thetapol_wall(i)',
c     &             i,theta,thetapol_wall(i)
        if (i.eq.1) then
           if((theta.le.thetapol_wall(1)).and.(theta.ge.0.d0)) then
              n_rho_wall=n_rho_wall+1
cSAP090321
c              i_wall_ar_m(n_rho_wall)=0
              i_wall_ar_m(n_rho_wall)=n_wall-1
              i_wall_ar_p(n_rho_wall)=1
           endif
           goto 10
        endif

        if (i.eq.n_wall) then
          if((theta.gt.thetapol_wall(n_wall-1)).and.
     &      (theta.le.2.d0*pi)) then
          n_rho_wall=n_rho_wall+1
          i_wall_ar_m(n_rho_wall)=n_wall-1
          i_wall_ar_p(n_rho_wall)=n_wall
          endif
          goto 10
        endif

        if(((theta.gt.thetapol_wall(i-1)).and.
     &      (theta.le.thetapol_wall(i))).or.
     &     ((theta.gt.thetapol_wall(i)).and.
     &      (theta.le.thetapol_wall(i-1)))) then
          n_rho_wall=n_rho_wall+1
          i_wall_ar_m(n_rho_wall)=i-1
          i_wall_ar_p(n_rho_wall)=i
          write(*,*)'n_rho_wall',n_rho_wall
        endif

 10   continue

      enddo
c--------------------------------------------------------------------------
c
c     equation for the wall small radius versus poloidal angle theta
c
c     rho=rho_wall_m+(theta-thetapol_wall_m)/(thetapol_wall_p-thetapol_wall_m)    
c      write(*,*)'in rho_wall_at_theta theta,n_rho_wall',theta,n_rho_wall
      
c-----calculate coordinates x_b,z_b at the LCFS surface psi=psilim
c     and the poloidal angle theta
      call zr_psith(psilim,theta,z_b,r_b)

c      write(*,*)'Q point theta,z_b,r_b',theta,z_b,r_b

      l_z_b_r_b=dsqrt((r_b-xma)**2+(z_b-yma)**2)  ! the distance between
                                                  ! the magnetic axis and
                                                  ! LCFS at theta
c      write(*,*)'Q point dsqrt((r_q-xma)**2+(z_q-yma)**2)',
c     &                      dsqrt((r_q-xma)**2+(z_q-yma)**2)
c      write(*,*)'Q point l_z_b_r_b', l_z_b_r_b
c      write(*,*)'Q point rho',dsqrt((r_q-xma)**2+(z_q-yma)**2)/l_z_b_r_b

      do j=1, n_rho_wall

c         write(*,*)'j',j
c         write(*,*)'i_wall_ar_m(j)',i_wall_ar_m(j)
c         write(*,*)'r_wall(i_wall_ar_m(j)),z_wall(i_wall_ar_m(j))',
c     &              r_wall(i_wall_ar_m(j)),z_wall(i_wall_ar_m(j))
c         write(*,*)'rho_wall(i_wall_ar_m(j))',rho_wall(i_wall_ar_m(j))
c         write(*,*)'theta,thetapol_wall(i_wall_ar_m(j))',
c     &              theta,thetapol_wall(i_wall_ar_m(j))

c         write(*,*)'i_wall_ar_p(j)',i_wall_ar_p(j)
c         write(*,*)'r_wall(i_wall_ar_p(j)),z_wall(i_wall_ar_p(j))',
c     &              r_wall(i_wall_ar_p(j)),z_wall(i_wall_ar_p(j))
c         write(*,*)'rho_wall(i_wall_ar_p(j))',rho_wall(i_wall_ar_p(j))
c         write(*,*)'theta,thetapol_wall(i_wall_ar_p(j))',
c     &              theta,thetapol_wall(i_wall_ar_p(j))

         r_p=r_wall(i_wall_ar_p(j))
         z_p=z_wall(i_wall_ar_p(j))
         r_m=r_wall(i_wall_ar_m(j))
         z_m=z_wall(i_wall_ar_m(j))

c         write(*,*)'r_p,r_m,r_q',r_p,r_m,r_q
c         write(*,*)'z_p,z_m,z_q',z_p,z_m,z_q

         if (r_p.ne.r_m) then
            b=(z_p-z_m)/(r_p-r_m)
            a=z_m-b*r_m

            if (r_q.ne.xma) then
              b_q=(z_q-yma)/(r_q-xma)
              a_q=z_q-b_q*r_q
c-------------intersection point K
              r_k=-(a-a_q)/(b-b_q)
              z_k=a+b*r_k
            else
c-------------r_q.eq.xma
              r_k=xma
              z_k=a+b*r_k
            endif 

         else
c-----------r_m=r_p
            if (r_q.ne.xma)then 
               b_q=(z_q-yma)/(r_q-xma)
               a_q=z_q-b_q*r_q

               r_k=r_m
               z_k=a_q+b_q*r_k
            else
c-------------r_q.eq.xma
              r_k=r_m
              WRITE(*,*)'r_m.eq.r_p, r_1.eq_r_m'
              STOP ' rho_wall_at_theta_2'
            endif

         endif
  
c         write(*,*)'r_p,r_m,r_q,r_k',r_p,r_m,r_q,r_k
c         write(*,*)'z_p,z_m,z_q,z_k',z_p,z_m,z_q,z_k
c         write(*,*)'j,i_wall_ar_m(j),i_wall_ar_p(j)', 
c     &              j,i_wall_ar_m(j),i_wall_ar_p(j)
c         write(*,*)'rho_wall(i_wall_ar_m(j)),rho_wall(i_wall_ar_p(j))',
c     &              rho_wall(i_wall_ar_m(j)),rho_wall(i_wall_ar_p(j))

c-----calculate coordinates x_b,z_b at the LCFS surface psi=psilim
c     and the poloidal angle theta at the K point
         x_l=r_k-xma 
         y_l=z_k-yma
         call theta_xy(x_l,y_l,theta)
         call zr_psith(psilim,theta,z_b,r_b)   
c         write(*,*)'r_k,z_k',r_k,z_k
c         write(*,*)'K point theta,z_b,r_b',theta,z_b,r_b
         l_z_b_r_b=dsqrt((r_b-xma)**2+(z_b-yma)**2)  ! the distance between
                                                     ! the magnetic axis and
                                                     ! LCFS at theta
c         write(*,*)'K point dsqrt((r_k-xma)**2+(z_k-yma)**2)',
c     &                      dsqrt((r_k-xma)**2+(z_k-yma)**2)

c         write(*,*)'K point l_z_b_r_b',l_z_b_r_b

         rho_wall_ar(j)=dsqrt((r_k-xma)**2+(z_k-yma)**2)/l_z_b_r_b

   

c         write(*,*)'rho_wall_ar(j)',rho_wall_ar(j)

      enddo


      return
      end


      subroutine create_fine_mesh_for_chamber_wall_limiter_coordinates 
c     It uses input wall arrays r_wall(n_wall_a),z_wall(n_wall_a)
c     from one_nml.i
c
c     Using linear approximation it creates new arrays 
c
c     r_wall_add(n_wall_add(m),max_limiters),z_wall_add(n_wall_add(max_limiters)
c
c     At m=0 it creates mesh for the chamber wall
c     At m>0 it creates mesh for the limiter with niumber 1=< m =< max_limiters
c
c     The created meshes ywill be in  fourb.i
c
c     The distance between additional point and old mesh points
c     is  h_add_wall*k 


      implicit none
      include 'param.i'
      include 'one.i'
      include 'fourb.i'

c-----locals
      integer i,j,k,ki,m,n,
     &n_wall_all   !number of points in r_wall_add and z_wall_add arrays

      real*8 
     &length_wall,               !poloidal wall length    
     &delta_r_wall,delta_z_wall, !r,z lengths of the wall element
     &delta_length_wall

      write(*,*)'create_fine_mesh_for_chamber_wall_limiter_coordinates' 
      
c-------------------------------------------------------------------------------
c     add chamber wall points at m=0
c-------------------------------------------------------------------------------
      length_wall=0.d0
      do i=1,n_wall-1
         length_wall=length_wall+dsqrt((r_wall(i+1)-r_wall(i))**2+
     &                         (z_wall(i+1)-z_wall(i))**2)
      enddo
      
      j=0
      m=0
      do i=1,n_wall-1
         delta_r_wall=r_wall(i+1)-r_wall(i)
         delta_z_wall=z_wall(i+1)-z_wall(i)
         delta_length_wall=dsqrt(delta_r_wall**2+delta_z_wall**2)

         j=j+1
         r_wall_add(j,m)=r_wall(i)
         z_wall_add(j,m)=z_wall(i)

         if (delta_length_wall.gt.h_add_wall) then
c-----------add additional wall points
          

            ki= delta_length_wall/h_add_wall ! number of additional wall points
                                             ! at the interval(i,i+1)
            do k=1,ki
               j=j+1
               
               if (j.gt.n_wall_add_a) then
                  WRITE(*,*)'in create_fine_wall_mesh'
                  WRITE(*,*)'j.gt.n_wall_add_a'
                  WRITE(*,*)'j,n_wall_ad_a',j,n_wall_add_a
                  WRITE(*,*)'Please increase n_wall_add_a'
                  WRITE(*,*)'in param.i and recompile the code' 
                  STOP 'create_seldom_mesh_for_chamber_wall_coordinates'
               endif 
               r_wall_add(j,m)=r_wall(i)+
     &                 delta_r_wall*h_add_wall*k/delta_length_wall
               z_wall_add(j,m)=z_wall(i)+
     &                 delta_z_wall*h_add_wall*k/delta_length_wall
            enddo
          endif

      enddo

      n_wall_add(0)=j

c-------------------------------------------------------------------------------
c     add limiter points at m>0
c-------------------------------------------------------------------------------
      do m=1,max_limiters 
         j=0         
         do i=1,n_limiter(m)-1             
           delta_r_wall=r_limiter(i+1,m)-r_limiter(i,m)
           delta_z_wall=z_limiter(i+1,m)-z_limiter(i,m)  
           delta_length_wall=dsqrt(delta_r_wall**2+delta_z_wall**2)
           j=j+1
           r_wall_add(j,m)=r_limiter(i,m)
           z_wall_add(j,m)=z_limiter(i,m)  
           if (delta_length_wall.gt.h_add_wall) then
c-------------add additional limiter points
              ki= delta_length_wall/h_add_wall ! number of additional limiter points
                                               ! at the interval(i,i+1)
                                   ! at the interval(i,i+1)
              do k=1,ki
                 j=j+1
               
                 if (j.gt.n_wall_add_a) then
                    WRITE(*,1000)
            STOP 'create_fine_mesh_for_chamber_wall_limiter_coordinates'
                 endif
 
                 r_wall_add(j,m)=r_limiter(i,m)+
     &                 delta_r_wall*h_add_wall*k/delta_length_wall
                 z_wall_add(j,m)=z_limiter(i,m)+
     &                 delta_z_wall*h_add_wall*k/delta_length_wall
              enddo !k
           endif
         enddo !i=0,n_limiter(m)   
         n_wall_add(m)=j
      enddo ! m=1,max_limiters
 1000 format('in create_fine_mesh_for_chamber_wall_limiter_coordinates'
     &,/,'at limiter calculations j.gt.n_wall_add_a',/,
     &'Please increase n_wall_add_a in param.i and recompile the code')


      write(*,*)'max_limiters,n_wall_add',max_limiters,n_wall_add
c      !stop 'create_fine_mesh_for_chamber_wall_limiter_coordinates' 
      return
      end
       
      subroutine distance_for_wall
c-----calculate distances  of RZ points determined  by  
c     arrays rr_add(nxeqd_add),zz_add(nyeqd_add)
c     1)from the chamber wall distance_to_wall(i,j,m=0)
c        at m=0
c     2)from {the chamber wall plus one limiter with number m} distance_to_wall(i,j,m>0)
c        at 1=<m=<max_limiters
c    

      implicit none

      include 'param.i'
      include 'fourb.i'
      include 'edge_prof_nml.i'

      include 'one_nml.i'
      integer i,j,k,m
      real*8  rho_l

c      write(*,*)'in distance_for_wall'
c      write(*,*)'n_wall_add',n_wall_add  

c------------------------------------------------------------------
       do i=1,nxeqd_add
         do j=1,nyeqd_add          

            do m=0,max_limiters     
c-------------------------------------------------------------------
c              m=0 distance to chamber wall,
c              m>1 distance to limiter
c------------------------------------------------------------------
               distance_to_wall(i,j,m)=1.d10 !initialization
               if (n_wall_add(m).eq.0) goto 10 ! gives huge distance_to_wall(i,j,m)=1.d10 

               do k=1,n_wall_add(m)
c                 write(*,*)'i,j,k,m',i,j,k,m
c                 write(*,*)'rr_add(i),r_wall_add(k,m)',
c     &                    rr_add(i),r_wall_add(k,m)
c                 write(*,*)'zz_add(i),z_wall_add(k,m)',
c     &                    zz_add(i),z_wall_add(k,m)

                  rho_l=dsqrt((rr_add(i)-r_wall_add(k,m))**2+
     &                     (zz_add(j)-z_wall_add(k,m))**2)

c                 write(*,*)'rho_l', rho_l

                  if (distance_to_wall(i,j,m).gt.rho_l) then
                     distance_to_wall(i,j,m)=rho_l
                  endif

c                  if(m.gt.0) then
cc-------------------distance to {wall plus limiter} 
c                    if(distance_to_wall(i,j,m).gt.
c     &                 distance_to_wall(i,j,0))
c     &                 distance_to_wall(i,j,m)=distance_to_wall(i,j,0)
c                  endif !m>0

               enddo !k

               if(m.gt.0) then
c-----------------distance to {wall plus limiter} 
                  if(distance_to_wall(i,j,m).gt.
     &               distance_to_wall(i,j,0))
     &               distance_to_wall(i,j,m)=distance_to_wall(i,j,0)
               endif !m>0

ctest-------------------------------------------------
c               if((m.ge.1).and.
c     &       (distance_to_wall(i,j,m).lt.distance_to_wall(i,j,0))) then
c               if(distance_to_wall(i,j,m).lt.(1.d-3)) then
c               write(*,*)'i,j,m,distance_to_wall(i,j,m),
c     &         distance_to_wall(i,j,0)',
c     &         i,j,m,distance_to_wall(i,j,m),distance_to_wall(i,j,0)
c               endif
c               endif
cend test----------------------------------------------

 10            continue

            enddo !m
          enddo !j
       enddo !i

c-test
c       do m=1,max_limiters
c        do i=1,nxeqd_add
c         do j=1,nyeqd_add
c           if(distance_to_wall(i,j,m).lt.2.d-3) then
c           write(*,*)'i,j,rr_add(i),zz_add(j),distance_to_wall(i,j,m)',
c     &               i,j,rr_add(i),zz_add(j),distance_to_wall(i,j,m)
c           endif
c          enddo !j
c        enddo !i
c       enddo !m  
c-end test

      return
      end

      subroutine density_at_zr_plane
c-----calculate density array density_r_z(nxeqd_add,nyeqd_add,nbulk,max_limiters)
c     at (rr_add,zz_add) mesh.
c
c     density_r_z(nxeqd_add,nyeqd_add,nbulk,m=0) is density with chamber wall effect only
c
c     density_r_z(nxeqd_add,nyeqd_add,nbulk,m>=1) is density with chamber wall effect plus
c                   limiter effect for limmiter with number m=1..,max_limiters
c                                     
c     It creates density fall near the chamber wall and limiters.
c     The code calculates factor
c     factor=1.d0-dexp(-(distance_to_wall(i,j)/sigma_wall_n)**2) or
c     factor=1.d0-dexp(-(distance_to_limiters(i,j,m)/sigma_wall_n)**2) or
c     which is equal to zero at the wall and it is equal to unit
c     far from the wall.
c
c     The final density is a product of this factor and
c     the density calculated without density fall. 

      implicit none

      include 'param.i'  
      include 'one.i'
      include 'fourb.i'
      include 'edge_prof_nml.i'

      integer i,j,k,m
      real*8 phi,dens_loc,factor,factor_limiter
c-----output
c     real*8  density_r_z(i,j,k,m)
c     This array is a density of k specie
c     at R(i)Z(j) mesh after product with calculted "factor".
c-----external
      real*8 rho_zr,b,dense_no_RZ_spline

      phi=0.d0

c-------------------------------------------------------------------------------
c     density with wall chamber effect
c------------------------------------------------------------------------------
      write(*,*)'in chamber_wall.f'
      write(*,*)'density_at_zr_plane n_wall=',n_wall
 
      do i=1,nxeqd_add ! over r
         do j=1,nyeqd_add !over z
          bmod= b(zz_add(j),rr_add(i),phi) !calculates small radius rho

c          write(*,*)'i,j,rr_add(i),zz_add(j),rho',
c     &               i,j,rr_add(i),zz_add(j),rho
          factor=1.d0 ! it will be equal to unit inside LCFS

          if(n_wall.gt.1) then
            if(rho.gt.1.d0) then
c-------------chamber wall effect
              factor=1.d0-
     &               dexp(-(distance_to_wall(i,j,0)/sigma_wall_n)**2)
        
c          write(*,*)'i,j,rho,distance_to_wall(i,j),sigma_wall_n,factor',
c     &               i,j,rho,distance_to_wall(i,j),sigma_wall_n,factor
            endif             
          endif

          do k=1,nbulk

            dens_loc=dense_no_RZ_spline(zz_add(j),rr_add(i),phi,k)

            do m=0,max_limiters
c----------------------------------------------------------------------
c               m=0 chamber wall effect only 
c               m=1,...max_limiters chamber wall plus limiter effects
c----------------------------------------------------------------------
               factor_limiter=1.d0 !far from limiter
               if (m.gt.0) then
c---------------- limiter effect

                  factor_limiter=1.d0-
     &            dexp(-(distance_to_wall(i,j,m)/sigma_wall_n)**2)

                  if(factor_limiter.lt.factor) factor=factor_limiter

               endif
 

               if (rho.gt.1.d0) then                    
                  density_r_z(i,j,k,m)= dens_loc*factor
               else
                  density_r_z(i,j,k,m)= dens_loc
               endif
            enddo !m
c            write(*,*)'k,dens_loc,density_r_z(i,j,k,m)',
c     &                 k,dens_loc,density_r_z(i,j,k,m)
c            write(*,*)'k,j,i,m,dens_loc,density_r_z(i,j,k,m)',
c     &                 k,j,i,m,dens_loc,density_r_z(i,j,k,m)
          enddo !k

        enddo !j
      enddo !i

      return
      end
     
      subroutine splcoef_density_r_z
c----------------------------------------------------------------------
c     calculate spline coefficients for density_r_z
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
c      parameter (dens_nwka=1+3*nxeqd_add_a)   
      
c      real*8      
c     &wk_dens(dens_nwka),dens_norm,dens_spline 

      real*8, dimension(1:3*max0(nxeqd_add,nyeqd_add)+1):: wk_dens
      real*8 dens_norm,dens_spline 

      integer ibd(4),idm
      integer i,j,k,m

c-----externals
      real*8 terp2p_Sm

     
      dens_nwka=1+3*max0(nxeqd_add,nyeqd_add)
      write(*,*)'splcoef_density_r_z'
      write(*,*)'nxeqd_add,nyeqd_add,dens_nwka',
     &           nxeqd_add,nyeqd_add,dens_nwka
c------------------------------------------------------------------
c     check parameter dens_nwka value 
c------------------------------------------------------------------
      if (dens_nwka.lt.(1+3*nxeqd_add)) then
         write(*,*)'**************************************' 
         WRITE(*,*)'chamber_wall,f in  splcoef_density_r_z'
         WRITE(*,*)'dens_nwka.lt.(1+3*nxeqd_add)'
         WRITE(*,*)'it shoud be'
         WRITE(*,*)'dens_nwka=1+3*max(nxeqd_add,nyeqd_add)'
         WRITE(*,*)'dens_nwka,nxeqd_add,nyeqd_add',
     &              dens_nwka,nxeqd_add,nxeqd_add
c         write(*,*)'Please change parameter in splcoef_density_r_z'
c         write(*,*)'for parameter (dens_nwka=1+3*nxeqd_add_a)'
c         write(*,*)'and recompile code'  
         write(*,*)'**************************************' 
         STOP
      endif

      if (dens_nwka.lt.(1+3*nyeqd_add)) then
         write(*,*)'**************************************' 
         WRITE(*,*)'chamber_wall,f in  splcoef_density_r_z'
         WRITE(*,*)'dens_nwka.lt.(1+3*nyeqd_add)'
         WRITE(*,*)'it shoud be'
         WRITE(*,*)'dens_nwka=1+3*max(nxeqd_add,nyeqd_add)'
         WRITE(*,*)'dens_nwka,nxeqd_add,nxyqd_add',
     &              dens_nwka,nxeqd_add,nxeqd_add
c         write(*,*)'Please change parameter in splcoef_density_r_z'
c         write(*,*)'for parameter (dens_nwka=1+3*nyeqd_add_a)'
c         write(*,*)'and recompile code'  
         write(*,*)'**************************************' 
         STOP
      endif
      
c-----creates 2D spline coefficients for density_r_z functions
      ibd(1)=2
      ibd(2)=2
      ibd(3)=2
      ibd(4)=2

      idm=nxeqd_add

      write(*,*)'chamber_wall.f in subroutine splcoef_density_r_z'
      write(*,*)'nbulk=',nbulk
      write(*,*)'nxeqd_add,nyeqd_add',nxeqd_add,nyeqd_add
c      write(*,*)'nxeqd_add_a,nyeqd_add_a',nxeqd_add_a,nyeqd_add_a
      do m=0,max_limiters

        do k=1,nbulk

          write(*,*)'k=',k,'m=',m
       
          do j=1,nyeqd_add
            density_r_z_rr(1,j,k,m)=0.d0
            density_r_z_rr(nxeqd_add,j,k,m)=0.d0
          enddo

          do i=1,nxeqd_add
            density_r_z_zz(i,1,k,m)=0.d0
            density_r_z_zz(i,nyeqd_add,k,m)=0.d0
          enddo

         write(*,*)'in subroutine splcoef_density_r_z before
     &              coeff2_Sm(nxeqd_add,rr_'

          call coeff2_Sm(nxeqd_add,rr_add,nyeqd_add,zz_add,
     &      density_r_z(1,1,k,m),
     &      density_r_z_rr(1,1,k,m),density_r_z_zz(1,1,k,m),
     &      density_r_z_rrzz(1,1,k,m),
     &      idm,ibd,wk_dens,nxeqd_add)

          write(*,*)'after coeff2_Sm(nxeqd_add,rr_'

        enddo !k
      enddo !m

c-test
c      do m=0,max_limiters
c      do k=1,nbulk    

c         write(*,*)'test k',k

c         dens_norm=0.d0    
c         do i=1,nxeqd_add
c            do j=1,nyeqd_add
c               write(*,*)'before terp2p k,j,i,density_r_z(i,j,k,0)'
c     &                                 ,k,j,i,density_r_z(i,j,k,0)

c              write(*,*)'before terp2p i,j',i,j
c              dens_spline=terp2p_Sm(rr_add(i),zz_add(j),
c     &                    nxeqd_add,rr_add,
c     &                    nyeqd_add,zz_add,
c     &             density_r_z(1,1,k,m),density_r_z_rr(1,1,k,m),
c     &             density_r_z_zz(1,1,k,m),
c     &             density_r_z_rrzz(1,1,k.m),idm,0,0,1,nxeqd_add) 
c              write(*,*)'after terp2p dens_spline',dens_spline
         
c              dens_norm=dens_norm+dabs(dens_spline-density_r_z(i,j,k,m))

c            enddo
c         enddo 
c         write(*,*)'k,dens_norm',k,dens_norm
c      enddo
c      enddo      
c-endtest

      return
      end


      real*8 function density_r_z_i_m(z,r,i,m)
c----------------------------------------------------------------------
c     calculate density using spline coefficients for density_r_z
c---------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input
      real*8 z,r !space coordinates [m]/r0x
      integer i,  !the number of plasma specie 
     &m           !number of limiter
                  !at m=0 no limiters 
      real*8      
     &dens_spline 

      integer idm

c-----externals
      real*8 terp2p_Sm

c      write(*,*)'function density_r_z_i'
            
      idm=nxeqd_add
c      write(*,*)'before terp2p_Sm'       
      density_r_z_i_m=terp2p_Sm(r,z,
     &                    nxeqd_add,rr_add,
cSAP110403
c  &                    nxeqd_add,zz_add
     &                  nyeqd_add,zz_add,
     &             density_r_z(1,1,i,m),density_r_z_rr(1,1,i,m),
     &             density_r_z_zz(1,1,i,m),
     &             density_r_z_rrzz(1,1,i,m),idm,0,0,1,nxeqd_add) 
c       write(*,*)'after terp2p_Sm'      
      return
      end

      real*8 function d_density_r_z_i_d_r_m(z,r,i,m)
c----------------------------------------------------------------------
c     calculate derivative d(density)/d(r)
c     using spline coefficients for density_r_z
c---------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input
      real*8 z,r !space coordinates [m]/r0x
      integer i,  !the number of plasma specie 
     &m           !number of limiter
                  !at m=0 no limiters 
      real*8      
     &dens_spline 

      integer idm

c-----externals
      real*8 terp2p_Sm

c      write(*,*)'function d_density_r_z_i_d_r i,m',i,m
           
      idm=nxeqd_add
          
      d_density_r_z_i_d_r_m=terp2p_Sm(r,z,
     &                    nxeqd_add,rr_add,
cSAP110304
c     &                    nxeqd_add,zz_add,
     &                    nyeqd_add,zz_add,
     &             density_r_z(1,1,i,m),density_r_z_rr(1,1,i,m),
     &             density_r_z_zz(1,1,i,m),
     &             density_r_z_rrzz(1,1,i,m),idm,1,0,1,nxeqd_add) 

      return
      end

 
      real*8 function d_density_r_z_i_d_z_m(z,r,i,m)
c----------------------------------------------------------------------
c     calculate derivative d(density)/d(z)
c     using spline coefficients for density_r_z
c---------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input
      real*8 z,r !space coordinates [m]/r0x
      integer i,  !the number of plasma specie 
     &m           !number of limiter
                  !at m=0 no limiters 
      real*8      
     &dens_spline 

      integer idm

c-----externals
      real*8 terp2p_Sm

c      write(*,*)'function d_density_r_z_i_d_z'
           
      idm=nxeqd_add
    
      d_density_r_z_i_d_z_m=terp2p_Sm(r,z,
     &                    nxeqd_add,rr_add,
cSAP110304
c     &                    nxeqd_add,zz_add,   
     &                    nyeqd_add,zz_add,
     &             density_r_z(1,1,i,m),density_r_z_rr(1,1,i,m),
     &             density_r_z_zz(1,1,i,m),
cSAP110322
c     &             density_r_z_rrzz(1,1,i,m),idm,0,1,1,nxeqd_add_a,m) 
     &             density_r_z_rrzz(1,1,i,m),idm,0,1,1,nxeqd_add) 
      return
      end

      real*8 function density_r_z_i(z,r,phi,i)
c----------------------------------------------------------------------
c     calculate density using spline coefficients for density_r_z
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input
      real*8 z,r,phi !space coordinates [m]/r0x
      integer i      !the number of plasma specie 
c-----external
      real*8 density_r_z_i_m
c-----local
      integer k, 
     &m,              !number of limiter
                      !at m=0 no limiters  
     &m_limiter,      !=0  the point is outside all limiters
                      !number of limiter (like m)
     &inside_limiter, !=1 the point is inside limiter with number 'm_limiter'
                      !=0  the point is utside all limiters
     &nearest_lim_boundary    !=1 or 2 number of nearest to the ray point
                              ! toroidal limiter boudndary
                              ! phi_limiter(nearest_lim_boundary,m_limiter)

      real*8 phi_loc, !0 =< phi_loc <2*pi
     &p,
     &dif_phi,dens_lim,d_dens_lim_d_r,
     &d_dens_lim_d_z,d_dens_lim_d_phi

c-----transform the toroidal angle phi
c     to the toroidal angle phi_loc: 0 =< phi_loc <2pi
      call phi_0_2pi(phi,phi_loc)

      p=pi/180.d0

      m=0

c------------------------------------------------------------------
c     Check if the toroidal andgle phi_loc
c     is inside location of any limiter
c     phi_limiter(1,m) =< phi_loc < phi_limiter(2,m)
c
c     If m=0 phi_loc is outside the toroidal location of any limiter
c---------------------------------------------------------------------
      do k=1,max_limiters

c         write(*,*)'density_r_z_i k,phi_limiter(1,k),phi_limiter(2,k)',
c     &                            k,phi_limiter(1,k),phi_limiter(2,k)
c         write(*,*)'phi,phi_loc', phi,phi_loc

         if ((phi_limiter(1,k)*p.le.phi_loc).and.
     &      (phi_loc.lt.phi_limiter(2,k)*p))then
            m=k
         endif
      enddo  

c      write(*,*)'density_r_z_i,z,r,phi,phi_loc,m',
c     &                         z,r,phi,phi_loc,m

cSAP090425
      inside_limiter=0
      if (max_limiters.gt.0) then
         call check_if_point_inside_limiter(z,r,phi,
     &   inside_limiter,m_limiter,dif_phi,nearest_lim_boundary)
      endif

      if (inside_limiter.eq.1)then
c--------the point is inside limiter

         call density_fall_at_toridal_limiter_boundary(z,r,phi,i,
     &   m_limiter,dens_lim,d_dens_lim_d_r,
     &   d_dens_lim_d_z,d_dens_lim_d_phi)

         density_r_z_i=dens_lim
      else 
c--------the point outside is limiter         
         density_r_z_i=density_r_z_i_m(z,r,i,m)
      endif
 
      return
      end

      real*8 function d_density_r_z_i_d_z(z,r,phi,i)
c----------------------------------------------------------------------
c     calculate density derivative using spline coefficients for density_r_z
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input
      real*8 z,r,phi !space coordinates [m]/r0x
      integer i      !the number of plasma specie 
c-----external
      real*8 d_density_r_z_i_d_z_m
c-----local
      integer k, 
     &m,           !number of limiter
                  !at m=0 no limiters 
     &m_limiter,      !=0  the point is outside all limiters
                      !number of limiter (like m)
     &inside_limiter, !=1 the point is inside limiter with number 'm_limiter'
                      !=0  the point is utside all limiters
     &nearest_lim_boundary    !=1 or 2 number of nearest to the ray point
                              ! toroidal limiter boudndary
                              ! phi_limiter(nearest_lim_boundary,m_limiter)
      real*8 phi_loc, !0 =< phi_loc <2*pi
     &p,
     &dif_phi,dens_lim,d_dens_lim_d_r,
     &d_dens_lim_d_z,d_dens_lim_d_phi

      call phi_0_2pi(phi,phi_loc)

      p=pi/180.d0

      m=0
      do k=1,max_limiters
         if ((phi_limiter(1,k)*p.le.phi_loc).and.
     &      (phi_loc.lt.phi_limiter(2,k)*p))then
            m=k
         endif        
      enddo   

cSAP090425
      inside_limiter=0
      if (max_limiters.gt.0) then
         call check_if_point_inside_limiter(z,r,phi,
     &   inside_limiter,m_limiter,dif_phi,nearest_lim_boundary)
      endif

      if (inside_limiter.eq.1)then
c--------the point is inside limiter

         call density_fall_at_toridal_limiter_boundary(z,r,phi,i,
     &   m_limiter,dens_lim,d_dens_lim_d_r,
     &   d_dens_lim_d_z,d_dens_lim_d_phi)

         d_density_r_z_i_d_z=d_dens_lim_d_z
      else 
c--------the point is outside limiter         
         d_density_r_z_i_d_z=d_density_r_z_i_d_z_m(z,r,i,m)
      endif

      return
      end     
 
      real*8 function d_density_r_z_i_d_r(z,r,phi,i)
c----------------------------------------------------------------------
c     calculate density derivative using spline coefficients for density_r_z
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input
      real*8 z,r,phi !space coordinates [m]/r0x
      integer i      !the number of plasma specie 
c-----external
      real*8 d_density_r_z_i_d_r_m
c-----local
      integer k, 
     &m,           !number of limiter
                  !at m=0 no limiters 
     &m_limiter,      !=0  the point is outside all limiters
                      !number of limiter (like m)
     &inside_limiter, !=1 the point is inside limiter with number 'm_limiter'
                      !=0  the point is utside all limiters
     &nearest_lim_boundary    !=1 or 2 number of nearest to the ray point
                              ! toroidal limiter boudndary
                              ! phi_limiter(nearest_lim_boundary,m_limiter)

      real*8 phi_loc, !0 =< phi_loc <2*pi
     &p,
     &dif_phi,dens_lim,d_dens_lim_d_r,
     &d_dens_lim_d_z,d_dens_lim_d_phi

      call phi_0_2pi(phi,phi_loc)

      p=pi/180.d0

      m=0
      do k=1,max_limiters
         if ((phi_limiter(1,k)*p.le.phi_loc).and.
     &      (phi_loc.lt.phi_limiter(2,k)*p))then
            m=k
         endif        
      enddo   

cSAP090425
      inside_limiter=0
      if (max_limiters.gt.0) then
         call check_if_point_inside_limiter(z,r,phi,
     &   inside_limiter,m_limiter,dif_phi,nearest_lim_boundary)
      endif

      if (inside_limiter.eq.1)then
c--------the point is inside limiter

         call density_fall_at_toridal_limiter_boundary(z,r,phi,i,
     &   m_limiter,dens_lim,d_dens_lim_d_r,
     &   d_dens_lim_d_z,d_dens_lim_d_phi)

         d_density_r_z_i_d_r=d_dens_lim_d_r
      else 
c--------the point is outside limiter         
         d_density_r_z_i_d_r=d_density_r_z_i_d_r_m(z,r,i,m)
      endif

      return
      end     


      real*8 function d_density_r_z_i_d_phi(z,r,phi,i)
c----------------------------------------------------------------------
c     calculate density derivative using spline coefficients for density_r_z
      implicit none
      include 'param.i'
      include 'fourb.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input
      real*8 z,r,phi !space coordinates [m]/r0x
      integer i      !the number of plasma specie 
c-----external
      real*8 d_density_r_z_i_d_r_m
c-----local
      integer k, 
     &m,           !number of limiter
                  !at m=0 no limiters 
     &m_limiter,      !=0  the point is outside all limiters
                      !number of limiter (like m)
     &inside_limiter, !=1 the point is inside limiter with number 'm_limiter'
                      !=0  the point is utside all limiters
     &nearest_lim_boundary    !=1 or 2 number of nearest to the ray point
                              ! toroidal limiter boudndary
                              ! phi_limiter(nearest_lim_boundary,m_limiter)

      real*8 phi_loc, !0 =< phi_loc <2*pi
     &p,
     &dif_phi,dens_lim,d_dens_lim_d_r,
     &d_dens_lim_d_z,d_dens_lim_d_phi

      call phi_0_2pi(phi,phi_loc)

      p=pi/180.d0

      m=0
      do k=1,max_limiters
         if ((phi_limiter(1,k)*p.le.phi_loc).and.
     &      (phi_loc.lt.phi_limiter(2,k)*p))then
            m=k
         endif        
      enddo   

cSAP090425
      inside_limiter=0
      if (max_limiters.gt.0) then
         call check_if_point_inside_limiter(z,r,phi,
     &   inside_limiter,m_limiter,dif_phi,nearest_lim_boundary)
      endif

      if (inside_limiter.eq.1)then
c--------the point is inside limiter

         call density_fall_at_toridal_limiter_boundary(z,r,phi,i,
     &   m_limiter,dens_lim,d_dens_lim_d_r,
     &   d_dens_lim_d_z,d_dens_lim_d_phi)

         d_density_r_z_i_d_phi=d_dens_lim_d_phi
      else 
c--------the point is outside limiter         
         d_density_r_z_i_d_phi=0.d0
      endif

      return
      end     
 
 
      subroutine creat_fine_2D_poloidal_mesh 
c-----creates at the  tokamak poloidal plane
c     arrays of (r,z) coordinates rr_add(nxeqd_add),zz_add(nyeqd_add)
c     with small distance betwen points
c     It is for the creation the density fall near the wall
c 

      implicit none

      include 'param.i'
      include 'three.i'
      include 'fourb.i'
      include 'edge_prof_nml.i'
      integer i 
      real*8 dstep

      dstep=xdimeqd/(nxeqd_add-1)
      do i=1,nxeqd_add
        rr_add(i)=redeqd+dstep*(i-1)
      enddo

      dstep=ydimeqd/(nyeqd_add-1)
      do i=1,nyeqd_add
        zz_add(i)=ymideqd+dstep*(i-1)-ydimeqd*0.5d0
      enddo 

      write(*,*)' creat_fine_2D_poloidal_mesh '
      write(*,*)'nxeqd_add,nyeqd_add',nxeqd_add,nyeqd_add
c      write(*,*)'rr',rr
c      write(*,*)'rr_add',rr_add
c      write(*,*)'zz',zz
c      write(*,*)'zz_add',zz_add
      return
      end
       
      subroutine phi_0_2pi(phi_in,phi_out)
c-----calculates phi_out  0 =< phi_out < 2*pi
      implicit none
c-----input
      real*8 phi_in  ! radian
c-----output      
      real*8 phi_out ! 0 =< phi_out < 2*pi
       
c-----locals 
      real*8 two_pi,phi_p
      integer k

      two_pi=8.d0*datan(1.d0)
 
      if (phi_in.ge.0) then 
         k=phi_in/two_pi 
         phi_out=phi_in-k*two_pi
      else
         phi_p=-phi_in
         k=phi_p/two_pi
         phi_out=two_pi-(phi_p-k*two_pi)
      endif   

      return
      end


      subroutine add_horizontal_limiter_walls
c---------------------------------------------------------------------
c     Add two horizontal limiter walls for each given limiters
c     These horizontal line will connect the top and the bottom points
c     with the right RZ mesh boundary  rr(nxeqd)
c-------------------------------------------------------------------
      implicit none 

      include 'param.i'
      include 'one.i'
      include 'fourb.i'
      include 'three.i'

c-----locals
      integer m,i,j,k,ki,j1
      real*8  step,delta,temp_r,temp_z
      real*8, dimension(1:n_wall_add_a) ::  temp_r_ar,temp_z_ar
      real*8, dimension(1:n_wall_add_a) ::  tem_r_ar,tem_z_ar
         

c      write(*,*)'add_horizontal_limiter_walls n_limiter_a',n_limiter_a

      do m=1,max_limiters       
c------------------------------------------------------------------------------
c       add additional wall points before the first given limiter point
c-----------------------------------------------------------------------------   
c        write(*,*)'m,n_wall_add(m)',m,n_wall_add(m)
c        do j=1,n_wall_add(m)
c           write(*,*)'1 j,r_wall_add(j,m),z_wall_add(j,m)',
c     &                  j,r_wall_add(j,m),z_wall_add(j,m)
c        enddo

        delta=rr(nxeqd)-r_wall_add(1,m)
        temp_r=rr(nxeqd)
        temp_z=z_wall_add(1,m)

        tem_r_ar(1)= temp_r
        tem_z_ar(1)= temp_z
        j=0
        if (delta.gt.h_add_wall) then
           ki= delta/h_add_wall ! number of additional wall points
                                ! at the old interval(1,2)
c           write(*,*)'ki',ki

           do k=1,ki
              j=j+1
c              write(*,*)'j',j
              tem_r_ar(j)=temp_r-h_add_wall*k
              tem_z_ar(j)=temp_z
           enddo
        endif
        j1=j      !number of added points

c        write(*,*)'j1',j1
c        do j=1,j1
c          write(*,*)'j,tem_r_ar(j),tem_z_ar(j)',
c     &               j,tem_r_ar(j),tem_z_ar(j)
c        enddo
c-------------------------------------------------------------------------------
c       put found additional points into arrays r_wall_add(*.m),r_wall_add(*,m)
c-------------------------------------------------------------------------------
        do j=1,n_wall_add(m)
c----------shift limiter points index at j1
           
           temp_r_ar(j+j1)=r_wall_add(j,m)
           temp_z_ar(j+j1)=z_wall_add(j,m)
        enddo

        do j=1,n_wall_add(m)+j1
           if(j.le.j1)then
              r_wall_add(j,m)=tem_r_ar(j)
              z_wall_add(j,m)=tem_z_ar(j)
           else
              r_wall_add(j,m)=temp_r_ar(j)
              z_wall_add(j,m)=temp_z_ar(j)
           endif
        enddo
        n_wall_add(m)=n_wall_add(m)+j1

c        write(*,*)'n_wall_add(m)',n_wall_add(m)
c        do j=1,n_wall_add(m)
c          write(*,*)'2 j,r_wall_add(j,m),r_wall_add(j,m)',
c     &                 j,r_wall_add(j,m),z_wall_add(j,m)
c        enddo
c---------------------------------------------------------------------------
c       add additional wall points after last given limiter point
c---------------------------------------------------------------------------
        delta=rr(nxeqd)-r_wall_add(n_wall_add(m),m)
        temp_r=r_wall_add(n_wall_add(m),m)
        temp_z=z_wall_add(n_wall_add(m),m)

c        write(*,*)'r_wall_add(n_wall_add(m),m)',
c     &             r_wall_add(n_wall_add(m),m)
c        write(*,*)'z_wall_add(n_wall_add(m),m)',
c     &             z_wall_add(n_wall_add(m),m)

        j=0
        if (delta.gt.h_add_wall) then
           ki= delta/h_add_wall ! number of additional wall points
                                ! at the old interval(1,2)
           do k=1,ki

              j=j+1
              tem_r_ar(j)=temp_r+h_add_wall*k
              tem_z_ar(j)=temp_z
           enddo
        endif
        j1=j      !number of added points

c        write(*,*)'j1',j1
c        do j=1,j1
c          write(*,*)'j,tem_r_ar(j),tem_z_ar(j)',
c     &               j,tem_r_ar(j),tem_z_ar(j)
c        enddo
c------------------------------------------------------------------------------
c       put found additional points into arrays r_wall_add(*.m),r_wall_add(*,m)
c------------------------------------------------------------------------------

        do j=n_wall_add(m)+1,n_wall_add(m)+j1
           r_wall_add(j,m)=tem_r_ar(j-n_wall_add(m))
           z_wall_add(j,m)=temp_z

c           write(*,*)'j,r_wall_add(j,m),z_wall_add(j,m)',
c     &                j,r_wall_add(j,m),z_wall_add(j,m)

        enddo
        n_wall_add(m)=n_wall_add(m)+j1+1
        r_wall_add(n_wall_add(m),m)=rr(nxeqd)
        z_wall_add(n_wall_add(m),m)=temp_z

c        write(*,*)'n_wall_add(m)',n_wall_add(m)
c        do j=1,n_wall_add(m)
c          write(*,*)'3j,r_wall_add(j,m),z_wall_add(j,m)',
c     &                j,r_wall_add(j,m),z_wall_add(j,m)
c        enddo

      enddo !m
          
c      !stop 'add_horizontal_limiter_walls'

      return
      end


      subroutine rho_limiter_at_theta(z_q,r_q,m,rho_limiter_theta_q,
     &n_rho_limiter)

c-----calculates limiter small radius at poloidal angle theta [radian]
c     of the given ray point (z_q,r_q)
c     It is proposed that between two wall neighbor points
c     wall small radius is a linear function on poloidal angle
c
c     inside this subroutine -pi=<theta=thetapol<pi

      implicit none
      include 'param.i'
      include 'one_nml.i'
      include 'three.i'
      include 'fourb.i'
c-----input
      real*8
     &z_q,r_q    !space coordinates of the ray point
      integer m  !the number of the given limiter 

c-----output
      real*8 
     &rho_limiter_theta_q !small radius at the poloidal angle theta(z_q,r_q)
      integer   n_rho_limiter !=0 intersection of the straight line OQ 
                              !   with the limiter was no found.
                              !   Here OA is the stright line connecting the magnetic axis
                              !   with the given Q point).
                              !   In this case small radius rho_limiter_theta_q
                              !   was not found.
                              !=1 the intersection and small radius rho_limiter_theta_q
                              !   were found    
c-----locals
      integer i,i_limiter_m,i_limiter_p
      real*8 r_p,z_p,r_m,z_m,thetapol_m,thetapol_p,theta_q

c-----externals
      real*8 thetapol

      theta_q=thetapol(z_q,r_q)                      !-pi=<theta_pol<pi 

c------------------------------------------------------------------------
c     determine 
c     n_rho_limiter
c     i_limiter_m=i-1
c     i_limiter_p=i     
c-------------------------------------------------------------------------
      n_rho_limiter=0

c      write(*,*)'rho_limiter_at_theta m,theta_q,n_limiter(m)',
c     &                                m,theta_q,n_limiter(m)

      do i=2,n_limiter(m)    

c         write(*,*)'i,thetapol_limiter(i-1,m),thetapol_limiter(i,m)',
c     &              i,thetapol_limiter(i-1,m),thetapol_limiter(i,m)

         if (((theta_q.gt.thetapol_limiter(i-1,m)).and.
     &       (theta_q.le.thetapol_limiter(i,m))).or.
     &       ((theta_q.ge.thetapol_limiter(i,m)).and.
     &       (theta_q.lt.thetapol_limiter(i-1,m)))) then
             n_rho_limiter=1
             i_limiter_m=i-1
             i_limiter_p=i
             goto 10
         endif   
      enddo

      if (n_rho_limiter.eq.0) then
c--------It is no interection of the stright line QO with limiter.
c        So, could not find  i_limiter_m and i_limiter_p
         goto 20
      endif

 10   continue
c------------------------------------------------------------------------
c      write(*,*)'in rho_limiter_at_theta n_rho_limiter',n_rho_limiter
c      write(*,*)'i_limiter_p,i_limiter_m,n_limiter,n_limiter_a',
c     &            i_limiter_p,i_limiter_m,n_limiter,n_limiter_a

      r_p=r_limiter(i_limiter_p,m)
      z_p=z_limiter(i_limiter_p,m)
      r_m=r_limiter(i_limiter_m,m)
      z_m=z_limiter(i_limiter_m,m)

      thetapol_m=thetapol(z_m,r_p)    
      thetapol_p=thetapol(z_p,r_p)    

      rho_limiter_theta_q=rho_limiter(i_limiter_m,m)+
     &    (rho_limiter(i_limiter_p,m)-rho_limiter(i_limiter_m,m))/
     &    (thetapol_p-thetapol_m)

 20   continue

      return
      end


      subroutine check_if_point_inside_limiter(z_q,r_q,phi_q,
     &inside_limiter,m_limiter,dif_phi,nearest_lim_boundary)
c-----It checks conditions that Q poit is inside limiter
     
      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'fourb.i'

c-----input
      real*8 z_q,r_q,phi_q ! coordinates of Q point
c-----output
      integer inside_limiter, !=1 the point is inside limiter with number 'm_limiter'
                              !=0  the point is outside all limiters
     &m_limiter,              !number of limiter
     &nearest_lim_boundary    !=1 or 2 number of nearest to the ray point
                              ! toroidal limiter boudndary
                              ! phi_limiter(nearest_lim_boundary,m_limiter)
      real*8 dif_phi
c-----locals
      integer m,n_rho_limiter,i_z_limiter,j

      real*8 rho_q,rho_limiter_theta_q,phi_q_loc,p

c-----externals: rho_limiter_at_theta, phi_0_2pi
      real*8 rho_zr

c----------------------------------------------------------------
c     calculate toroidal angle phi_q_loc:  0 = <phi_q_loc < 2pi using
c     the input toroidal angle phi_q
c-----------------------------------------------------------------------
      call phi_0_2pi(phi_q,phi_q_loc)
c-----------------------------------------------------------------------
c     calculate small radius of Q point
c----------------------------------------------------------------------
      rho_q= rho_zr(z_q,r_q)
      
      inside_limiter=0

      if(rho_q.le.1.d0) goto 20
      

      p=pi/180.d0
c      write(*,*)'check_if_point_inside_limiter z_q,r_q,phi_q',
c     &z_q,r_q,phi_q
c      write(*,*)'phi_q_loc,rho_q',phi_q_loc,rho_q

      do m=1,max_limiters

c          write(*,*)'m,phi_limiter(1,m)*p,phi_limiter(2,m)*p',
c     &               m,phi_limiter(1,m)*p,phi_limiter(2,m)*p

          if (((phi_limiter(1,m)*p.le.phi_q_loc).and.
     &         (phi_q_loc.lt.phi_limiter(2,m)*p)).or.
     &        ((phi_limiter(2,m)*p.le.phi_q_loc).and.
     &         (phi_q_loc.lt.phi_limiter(1,m)*p))) then
c------------the toroidal angle of Q point is inside 
c            the toroidal interval of 'm' limiter localization'
c
c            phi_q is inside (phi_limiter(m,1),phi_limiter(m,2))
             
c             write(*,*)'z_limiter(1,m),z_limiter(n_limiter(m),m),z_q',
c     &                  z_limiter(1,m),z_limiter(n_limiter(m),m),z_q

             i_z_limiter=0
             if (((z_limiter(1,m).lt.z_q).and.
     &            (z_q.lt.z_limiter(n_limiter(m),m))).or.
     &           ((z_limiter(n_limiter(m),m).lt.z_q).and.
     &            (z_q.lt.z_limiter(1,m)))) then
c----------------z_q is inside (z_limiter(1,m),z_limiter(n_limiter(m),m))
                 i_z_limiter=1
             else
                 i_z_limiter=0 ! outside z interval
                goto 10
             endif
c------------------------------------------------------------------
c            calculates limiter small radius at poloidal angle theta [radian]
c            of the given ray point (z_q,r_q)
c------------------------------------------------------------------            

c             write(*,*)'i_z_limiter',i_z_limiter

             if (i_z_limiter.eq.1) then
c---------------check that Q point is inside poloidal location of limiter
c--------------------------------------------------------------------------
c               calculate small radius rho_limiter_theta_q 
c               of the limiter intersection with the 
c               stright line connecting the magnetic axis
c               with Q(z_q,r_q) point
c----------------------------------------------------------------------------                  
                call rho_limiter_at_theta(z_q,r_q,m,rho_limiter_theta_q,
     &                                    n_rho_limiter)
c                write(*,*)'after rho_limiter_at_theta n_rho_limiter',
c     &                     n_rho_limiter

                if (n_rho_limiter.eq.1) then
c------------------the intersection and its small radius 
c                  rho_limiter_theta_q were found

c                   write(*,*)'rho_q,rho_limiter_theta_q',
c     &                        rho_q,rho_limiter_theta_q

                   if (rho_q.gt.rho_limiter_theta_q) then
c---------------------the Q point is inside limiter
                      m_limiter=m
                      inside_limiter=1
                      dif_phi=1.d9
                      do j=1,2
                        if (dabs(phi_q_loc-phi_limiter(j,m)*p).lt.
     &                     dif_phi) then
                           dif_phi=dabs(phi_q_loc-phi_limiter(j,m)*p)
                           nearest_lim_boundary=j
                        endif                         
                      enddo !j 

                      write(*,*)'nearest_lim_boundary',
     &                           nearest_lim_boundary  

                      write(*,*)'inside_limiter',inside_limiter

                      goto 20
                   endif
                endif
             endif
          endif
 10       continue
      enddo !m
 20   continue

      return
      end
      

      subroutine density_fall_at_toridal_limiter_boundary(z,r,phi,k,
     &m_limiter,dens_lim,d_dens_lim_d_r,d_dens_lim_d_z,d_dens_lim_d_phi)

c-----calculate density dens and its derivative
      implicit none
      include 'param.i'
      include 'one.i'
      include 'edge_prof_nml.i'
c-----input 
      real*8 
     &z,r,phi !space coordinates
     

      integer
     &k,         ! the number of plasma specie
     &m_limiter, !the number of the given limiter
     &nearest_lim_boundary    !=1 or 2 number of nearest to the ray point
                              ! toroidal limiter boudndary
                              ! phi_limiter(nearest_lim_boundary,m_limiter)
c-----output
      real*8 dens_lim,d_dens_lim_d_z,d_dens_lim_d_r,d_dens_lim_d_phi
c-----external
      real*8 density_r_z_i_m,
     &d_density_r_z_i_d_z_m,d_density_r_z_i_d_r_m
c-----locals
      real*8 factor,dens,d_dens_d_z,d_dens_d_r,
     &sigma_lim_toroidal_radian
      integer m 

      sigma_lim_toroidal_radian=sigma_lim_toroidal_degree*pi/180.d0

      factor=dexp(-(phi-phi_limiter(nearest_lim_boundary,m_limiter)/
     &              sigma_lim_toroidal_radian)**2)

      m=0 
      dens=density_r_z_i_m(z,r,k,m)
      d_dens_d_z=d_density_r_z_i_d_z_m(z,r,k,m)
      d_dens_d_r=d_density_r_z_i_d_r_m(z,r,k,m)

      dens_lim=dens*factor
      d_dens_lim_d_z=d_dens_d_z*factor
      d_dens_lim_d_r=d_dens_d_r*factor
      d_dens_lim_d_phi=dens_lim*2.d0*
     &(phi_limiter(nearest_lim_boundary,m_limiter)-phi)/
     &sigma_lim_toroidal_radian**2
     
      return
      end
