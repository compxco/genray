c        ********************** cone_ec *********************
c        *                      -----                       *
c        * this subroutine calculates initial data          *
c        * for ecr power distributon in cones		    *
c        ****************************************************
c
!        Called for each icone (inside icone=1,ncone loop)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c      	input parameters					    !
c	alpha1 - the angular width of the power distribution(radian)!
c       na1    - number of cones                                    !
c       na2    - number of rays per cone                            !
c       alpha2 - starting angle on cone(radian)			    !
c       phist  - antenna toroidal angle (radian)        	    !
c       phin   -toroidal angle of the central ray (radian)from e_r  !
c       tetan  -poloidal angle of the central ray (radian)from e_z  !
c       powtot-total power at antenna(in MW)			    !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c      	output parameters					    !
c       nray=na1*na2+1 number of rays				    !
!       Aiming angles for each ray:
c       alphaj(nray)-array: toroidal angles(radian) of j_ray	    !
!          (toroidal angle measured from R-vector through source)
!
c       betaj(nray)-array: angles(radian) between the  horizontal   !
c                          plane and j_ray    
!
c       powj(nray)-array:power flowing in the ray channel at    !
c               antenna normalized Sum_{i=1,nray}powj(i)=powtot !
c               (erg/sec)				            !
c--------------------------------------------------------------------
        subroutine cone_ec(alpha1,na1,na2,alpha2,phist,
     1                     phin,tetan,powtot,nray,alphaj,betaj,powj)
        implicit double precision (a-h,o-z)
c        dimension betaj(1),alphaj(1),powj(1)
        dimension betaj(*),alphaj(*),powj(*)
        real*8 n_r,n_phi,n_z, na1r, na2r
        write(*,*)'in cone_ec'
        pi=datan(1.0d0)*4.0d0
	tet=pi/180.0d0
	tetanr=tetan
	phinr=phin ! toroidal angle of the central ray (radian)from e_r
	phistr=phist

	alpha2r=alpha2
	alpha1r=alpha1

        ctetan=dcos(tetanr)
        stetan=dsin(tetanr)
	cphin=dcos(phinr)
	sphin=dsin(phinr)
	cphist=dcos(phistr)
	sphist=dsin(phistr)
c--------------------------------------------------------------------
c       number of the rays
	nray=na1*na2+1
	write(*,*)'in cone_ec na1,na2',na1,na2
	write(*,*)'number of the rays in EC cone: nray=na1*na2+1=',nray
c-----------------------------------------------------------------
c       central refractive index vector at antenna
	n_r=stetan*cphin
	n_phi=stetan*sphin
	n_z=ctetan
c----------------------------------------------------------------
c       only central ray in antenna cone
	if (na1.eq.0) then
	   da1=alpha1
	   goto 1
	endif
c----------------------------------------------------------------
	na1r=dfloat(na1)
	na2r=dfloat(na2)
	da1=alpha1r/na1r
	da2=2.0d0*pi/na2r
c-----------------------------------------------------------------
c       central ray
1	  iray=1
	  alphaj(1)=phinr
	  betaj(1)=0.5d0*pi-tetanr
	  write(*,*)'in cone_ec 1central ray iray,alphaj(1),betaj(1)',
     1	  iray,alphaj(1),betaj(1)
c----------------------------------------------------------------
c       only central ray data in antenna cone
	if (na1.eq.0) then
	  powj(1)=1.d0
c	  write(*,*)'only central ray na1,powj(1)',
c     1	  na1,powineci(1)
	  goto 30
	else
	  powj(1)=2.d0*pi*(1.0d0-dcos(0.5d0*da1))
c	  write(*,*)'na1.ne.0 powj(1)',powj(1)
	end if
c-----------------------------------------------------------------
c       cone rays around the central ray
c----------------------------------------------------------------
      alpha_avg=0.d0 !YuP[2020-07-24]: To evaluate average toroidal aiming angle
	do 10 i=2,na1+1
c          cone angles
           alpha1i=da1*(i-1)
c	   write(*,*)'alpha1i,i',alpha1i,i
           do 20 j=1,na2
c             angles around the cone
c           write(*,*)'in loop 20 j',j
              alpha2j=alpha2r+da2*(j-1)
	      iray=1+(i-2)*na2+j
c	   write(*,*)'alpha2j,j',alpha2j,j,'iray=',iray
c-----------------------------------------------------------
c             ray direction in e_n, e_phi, e_ksi system
c-----------------------------------------------------------
              cn_en=dcos(alpha1i)
              cn_eph=dsin(alpha1i)*dcos(alpha2j)
              cn_eksi=dsin(alpha1i)*dsin(alpha2j)
c-----------------------------------------------------------
c             ray direction in r, phi, z system
c-----------------------------------------------------------
	      cn_r=cn_en*stetan*cphin
     1             -cn_eksi*ctetan*cphin
     2		   -cn_eph*sphin

	      cn_phi=cn_en*stetan*sphin
     1               -cn_eksi*ctetan*sphin
     2		     +cn_eph*cphin

	      cn_z=cn_en*ctetan+cn_eksi*stetan
c-----------------------------------------------------------
c             ray direction in x,y,z system
c-----------------------------------------------------------
c	      cnx=cn_r*cphist-cn_phi*sphist
c	      cny=cn_r*sphist+cn_phi*cphist
c	      cnz=cn_z
c----------------------------------------------------------
    	  sbetaj=cn_z
	  if(cn_z.gt.1.0d0) sbetaj=1.0d0
	  if(cn_z.lt.-1.0d0) sbetaj=-1.0d0
	  betaj(iray)=dasin(sbetaj)

	  cbetaj=dcos(betaj(iray))
	  if (dabs(cbetaj).eq.0.0d0) then
	     alphaj(iray)=0.0d0
	     goto 2
	  end if

	  calphaj=cn_r/cbetaj
	  if(calphaj.gt.1.0d0) calphaj=1.0d0
	  if(calphaj.lt.-1.0d0) calphaj=-1.0d0

	  salphaj=cn_phi/cbetaj
	  if (salphaj.ge.0.0d0) then
	     alphaj(iray)=dacos(calphaj)
	  else
	     alphaj(iray)=-dacos(calphaj)
	  end if
 2        continue

cTJenkins	  powj(iray)=dexp(-2.0d0*(alpha1i/alpha1r))*
	  powj(iray)=dexp(-2.0d0*(alpha1i/alpha1r)**2)*
     1    (dcos(alpha1i-0.5d0*da1)-dcos(alpha1i+0.5d0*da1))*da2
     
           alpha_avg=alpha_avg+alphaj(iray) !YuP[2020-07-24]aver tor aiming angle

 20        continue ! j=1,na2
 10     continue ! i=2,na1+1
 
        alpha_avg=alpha_avg/float(nray) !YuP[2020-07-24] average tor aiming angle
 
        !YuP[2020-07-24] Compare average tor. aiming angle (alpha_avg)
        !over rays on the cone with tor. aiming angle alphaj(1) of central ray.
        if( (alpha_avg.lt.0.d0).and.(alphaj(1).gt.0.d0) )then
          ! If they are opposite, redefine tor. aiming 
          ! angle of the central ray :
          alphaj(1)= alphaj(1) -2.d0*pi
          ! It is not important for computations,
          ! but it is good to do this for saving data --> plots.
        endif
        if( (alpha_avg.gt.0.d0).and.(alphaj(1).lt.0.d0) )then
          alphaj(1)= alphaj(1) +2.d0*pi
        endif
 
c--------------------------------------------------------------
c       normalization of the angle distribution	 of the power
c--------------------------------------------------------------
 30	  continue ! handle to skip the above, in case of central ray only.
 
        psum=powj(1)
        do iray=2,nray
           psum=psum+powj(iray)
        enddo
c-----------------------------------------------------------
c  powtot is in MW
c  The transformation of total antenna power from MW to erg/sec
c-----------------------------------------------------------------
cSm050902        powtot=powtot*1.0d+13
c----------------------------------------------------------------
      p=powtot*1.0d+13/psum
      do iray=1,nray
         powj(iray)=powj(iray)*p
         write(*,*)'cone_ec iray,alphaj(iray),powj(iray)',
     &   iray,alphaj(iray),powj(iray)
      enddo
      
      return
      end subroutine cone_ec





c        **********************raypat**************************
c        *                        -                           *
c        * this subroutine calculates initial data            *
c        * for ecr power distributon in cones,		      *
c        * setting ray directions according to toray system   *
c        ******************************************************
c
c------------------------------------------------------------------

      subroutine raypat(thetac, phic, div, cr, nray, 
     & gzone, mray, theta, phi)
cProlog

c     From TORAY, Apr 5, 2004, editted down to 
c       R. Prater's Gaussian ray pattern generator.
c     A Guassian pattern of power injection, 
c       p(alpha)=p_norm*exp(-2.*(alpha/div)**2) is obtained,
c       where alpha is angle measured from the central ray.
c       Each ray is to be equi-weighted in power: the solid angle
c       density of ray trajectories gives the Gaussian power 
c       density profile.
c       The Gaussian distn is truncated at 0.98 of total power,
c         [else, adjust sfracpwr below].
c
c     Input:
c       thetac:   polar injection angle of central ray (degrees)
c       phic:     azimuthal injection angle of central ray (degrees)
c       div:      half-width at half power of the beam (in degrees)
c       gzone:    if 0 then 48 ray case, as specified by mray() below
c                 if 1 then there can only be 1 ray, the central ray.
c                 if .gt.1 then describes number of elements in mray
c       mray(*):  if gzone .gt.0, use gaussian formulation with this number
c          of rays in corresponding annular zone, otherwise use the
c          usual 1,5,12,12,18  (48 ray) arrangement.
c          mray(1) is effectively 1.
c       cr(*):    azimuthal phase of ray pattern for each zone, in radians;
c         same size array as mray for gaussian formulation.
c         Standard setting for gzone=0 is 0.0,0.1,0.05,-0.05,0.05.
c     outputs:
c       nray: total number of rays, using mray(*) if gzone.ne.0.
c       theta(nray): polar injection angle of a given ray (degrees)
c       phi(nray):   toroidal injection angle of a given ray (degrees)
c------------------------------------------------------------------------
      implicit none
c
      integer m, n, kk
      integer nray, mray(*), nzone, gzone, numrays
      !doubleprecision to real*8 !YuP[2020-01-27]
      real*8 r0, r1, r2, r3, ar, am, pi, sfracpwr
      real*8 d, ppz, f, df, ang, x, maxzone
      real*8 thetac, phic, div
      real*8 theta(*), phi(*), cr(*)
c
      pi=acos(-1.0d0)
c
      theta(1)=thetac
      phi(1)=phic
      numrays=1
cBH050603      if ((nray .eq. 1) .or. (gzone .eq. 1)) then
c      write(*,*)'raypat gzone',gzone 
      if (gzone .eq. 1) then
         nray = 1
         return
      endif
c     ***
      d=div/sqrt(log(2.0d0)/2.0d0)
c     Neglect lost power (ie, truncate the gaussian at 1-f)
      sfracpwr=0.98d0
      f=1.0d0-sfracpwr
      df=0.0d0
c
      if (gzone .gt. 1) then
         nzone=gzone
         numrays=1
         do m=2, nzone
            numrays=numrays+mray(m)
         enddo
      else
         mray(2)=5
         mray(3)=12
         mray(4)=12
         mray(5)=18
         cr(1)=0.0d0
         cr(2)=0.1d0
         cr(3)=0.05d0
         cr(4)=-0.05d0
         cr(5)=0.05d0
         numrays=48
         nzone=5
      endif
      nray=numrays
c     
      kk=numrays
      do n=nzone,2,-1
         ppz=mray(n)*sfracpwr/numrays
         f=f+0.5d0*(ppz+df)
         df=ppz
         ang=d*sqrt(-log(f)/2.0d0)
         do m=mray(n),1,-1
            x=pi*(2.0d0*float(m)/float(mray(n)) + cr(n))
            theta(kk)=thetac-ang*cos(x)
            phi(kk)=phic-ang*sin(x)
            kk=kk-1
         enddo
      enddo
c     
      return
      end

c        *********radial_integral_at_uniform_mesh**************
c        *                        -                           *
c        * this subroutine calculates integral{0,x}f(x)xdx    *
c        * for ecr power distribution at the launching disk   *
c        * using uniform radial greed                         *
c        ******************************************************
c
c------------------------------------------------------------------
      subroutine radial_integral_at_uniform_mesh(func,b,n,
     &x,xc,ar_integral,total_integral)
c-----calculates integral{0,x}f(x)xdx at the uniform mesh
c     and puts it in ar_integral
c     0=<  x <b
c     using the formula
c     Sum{k=2,j-1}[f(x(k))*((0.5*(x(k+1)+x(k)))**2-(0.5*(x(k)+x(k-1)))**2)],
c     for  j=2,n-1
c
c     This subroutine will be used to calculate non-uniform
c     radial mesh at the launching disk for raypatt='diskdisk' case
c--------------------------------------------------------------------
c     input fync(x) !double precision function
c                   !It should be described as external in the program
c                   !using call at_uniform_mesh
c
c     0 =< b        !the boundary of integration 
c
c     n is number of x mesh points
c
c-----output:
c     
c     x mesh, the points where the integral is calculated 
c     step h=b/(n-1.5)
c
c     0 h/2    h         h        h   b     
c     |----|--------|--------|--------|        (example for n=5)
c     x(1) x(2)     x(3)     x(4)     x(5)
c
c     x(1)=0, 
c     x(j)=h/2+(j-2)*h     Here j=2,n,   x(1)=0,x(2)=h/2,x(n)=b
c
c     xc mesh, the points where the function is calculated
c
c     0    h       h         h     h/2        
c     |--------|--------|--------|----b        (example for n=5)
c     xc(1)   xc(2)    xc(3)    xc(4)     
c
c     xc(1)=0
c     xc(j)=h(j-1), j=2,n-1
c 
c 
c
c     I(x(1))=ar_integral(1)=0.d0
c     I(x(j))=ar_integral(j)=Sum{k=1,j-1}[f(xc(k))*(x**2(k+1)-x**2(k))/2], 
c             j=2,n
c
c     total_integral=ar_integral(n)

      implicit none

c-----input      
      real*8 b !the radius of the launching disk  
      integer n !number of radial points of the uniform mesh
      double precision func !this function is the radial power distribution 
                            !at the launching disk

c-----output
      double precision x(n),xc(n-1),ar_integral(n),total_integral
        
c-----locals
      integer j,k
      double precision step 

      step=b/(n-1.5d0)

      x(1)=0.d0
      do j=2,n         
        x(j)=(j-1.5d0)*step
      enddo

      xc(1)=0.d0
      do j=2,n-1         
        xc(j)=step*(j-1)
      enddo

      do j=1,n
        ar_integral(j)=0.d0
      enddo

      do j=2,n
        ar_integral(j)=ar_integral(j-1)+func(xc(j-1))*
     &                 0.5d0*(x(j)**2-x(j-1)**2)
      enddo

      total_integral=ar_integral(n)
c-----normalization to get the total integral =1
      do j=1,n
        ar_integral(j)=ar_integral(j)/total_integral
      enddo
      total_integral=ar_integral(n)

      return
      end


c        ******** create_equal_power_radial_mesh **************
c        *                        -                           *
c        * this subroutine creates the non-uniform radial mesh*
c        * ar the EC launching disk                           *
c        * to get the equal launched power for each ray       *
c        ******************************************************
c
c------------------------------------------------------------------

      subroutine create_equal_power_radial_mesh(b,func,
     &n_mesh_radial_bin,n_mesh_angle,x_mesh,x_mesh_c,f_equal)

c-----creates X mesh and array for the function func:
c     n_mesh_radial_bin is the number of radial bins 
c     n_mesh=n_mesh_radial_bin + 1 is the number of radial points in which
c              the integral is calculated     
c     array x_mesh_c(n_mesh_radial_bin)  (0<  x_mesh_c <b) ! Bin centers
c     array x_mesh(n_mesh) (0=<  x_mesh =<b)   Bin edges (bin boundaries)
c
c     array f_equal(n_mesh_radial_bin) of the function func
c                                        at the created mesh
c
c     This mesh will be used to calculate the integral
c     integral{0,b}f(x)xdx=sum{j=1,n_mesh}f_equal(j)*del(j)
c
c     x_mesh(1)=0
c     x_mesh(n_mesh)=b
c     0=<x_mesh(j)=<b   j=1,...,n_mesh
c     for j=2,....n_mesh
c     del(j)=0.5(x_mesh(j)**2-x_mesh(j-1)**2) j=1,...,n_mesh
c     
c     x_mesh_c(1)=0
c     x_mesh_c(j)=0.5(x_mesh(j)+x-mesh(j-1)) j=2,..., n_mesh_radial_bin
c
c     In each radial bin (j= n_mesh_radial_bin1,n_mesh-1) there will be
c                         n_mesh_angle(j)
c     rays with different angles. In the first radial bin there will be 
c     one ray n_mesh_angle(1)=1
c
c     The total number of rays n_disk_rays
c     n_disk_rays = Sum{j=1,n_mesh_radial_bin}
c
c     The mesh and the function array will be calculated to give
c
c     f_equal(j)*del(j)=Constant=total_integral/n_mesh
c
c
c     total_integral is calculated from the given function func 
c     at the equispaced mesh with n points (n-1 bins)
c
c     h=b/(n-1.5)
c     x(1)=0
c     x(j)=(j-1.5)*h          Here j=2,n, x(1)=0, x(n)=b 
c
c     xc(1)=0
c     xc(j)=0.5(x(j)+x(j+1))  Here j=2,n-1, xc(1)=0, xc(n-1)=b-h/2
c
c     ar_integral(1)=0.d0
c     ar_integral(j)=Sum{k=1,j-1}[f(xc(k))*0.5(x(k+1)**2-x(k))**2], j=2,n
c
c     total_integral=ar_integral(n)
c
c-----------------------------------------------------------------------

      implicit none
     
      integer n !the number of points for the equispaced mesh x(n+1)
      parameter (n=1000)

c-----input
      integer n_mesh_radial_bin !the number of radial bins
      real*8 b                  ! the boundary of the integration area
      integer n_mesh_angle(n_mesh_radial_bin) !the number of angles
                                              ! at each radius bin 
c-----externals
      real*8 func
      external func

c-----output
      real*8 x_mesh(n_mesh_radial_bin+1),x_mesh_c(n_mesh_radial_bin),
     &f_equal(n_mesh_radial_bin)
c-----locals
      integer n_mesh ! the number of new mesh points
      integer j,k,k0,n_disk_rays,j_rays,i
      real*8 integral_value
      real*8 x(n),xc(n-1),ar_integral(n),total_integral ! local work arrays

      real*8  total_integral_mesh,total_integral_per_one_ray,
     &integral_per_one_ray
      real*8 p_expected ! local, for printout
c------------------------------------------------------------
c     calculate the integral from the function func
c     ar_integral(j) j=1,...,n
c     at the equispaced mesh x()

      write(*,*)'create_equal_radial_mesh b,n',b,n

c------------------------------------------------------------
      call radial_integral_at_uniform_mesh(func,b,n,x,xc,
     &ar_integral,total_integral)
c------------------------------------------------------------

      !write(*,*)'create_equal_mesh: ar_integral',ar_integral
      !write(*,*)'create_equal_mesh: x', x

c------------------------------------------------------------
c     x_mesh and the function array will be calculated to give
c
c     f_equal(j)*del(j)=Constant=total_integral/n_mesh
c---------------------------------------------------------
      n_mesh=n_mesh_radial_bin+1
      k0=1
      x_mesh(1)=0.d0 
      x_mesh(n_mesh)=b
      x_mesh_c(1)=0.d0

c-----calculate the total number of rays n_disk_rays
      n_disk_rays=1
      do j=2,n_mesh_radial_bin
        n_disk_rays=n_disk_rays+n_mesh_angle(j)
c        write(*,*)'j,n_mesh_angle(j)',j,n_mesh_angle(j)
      enddo
c      write(*,*)'n_disk_rays',n_disk_rays

c      write(*,*)'total_integral',total_integral

      integral_per_one_ray=total_integral/n_disk_rays

c      write(*,*)'integral_per_one_ray',integral_per_one_ray
      
      if (n_mesh_radial_bin.eq.1)then
        j_rays=1
        x_mesh_c(1)=0.d0 
        f_equal(1)=1.d0
        x_mesh(1)=0.d0
      else
        j_rays=0
        do j=2,n_mesh_radial_bin
c---------number of rays
          do i=1,n_mesh_angle(j-1) ! from n_mesh_disk_angle_bin(j-1)
            j_rays=j_rays+1
          enddo

          integral_value=j_rays*integral_per_one_ray

c         write(*,*)'j,integral_value',j,integral_value
          do k=k0,n+1

c            write(*,*)'k0,k,integral_value,ar_integral(k)',
c    &                  k0,k,integral_value,ar_integral(k)

            if(integral_value.le.ar_integral(k)) then
c---------------------------------------------------------------------
c             The equation to find xmesh(j)
c           (integral_value-ar_integral(k-1))/(x_mesh(j)**2-x(k-1)**2)=
c             =(ar_integral(k)-ar_integral(k-1))/(x_(k)**2-x(k-1)**2)
c-----------------------------------------------------------------------

c               write(*,*)'k,x(k-1),x(k)',k,x(k-1),x(k)

            x_mesh(j)=dsqrt(x(k-1)**2+(integral_value-ar_integral(k-1))*
     &         (x(k)**2-x(k-1)**2)/(ar_integral(k)-ar_integral(k-1)))

c               write(*,*)'j,x_mesh(j)',j,x_mesh(j)
c               write(*,*)'x_mesh',x_mesh            
           
               f_equal(j-1)=integral_per_one_ray/
     &                      (0.5d0*(x_mesh(j)**2-x_mesh(j-1)**2))

c               write(*,*)'j,integral_per_one_ray',
c     &                    j,integral_per_one_ray
c               write(*,*)'x_mesh(j-1),x_mesh(j),f_equal(j-1)',
c     &                    x_mesh(j-1),x_mesh(j),f_equal(j-1)
               if(j.gt.2) x_mesh_c(j-1)=0.5d0*(x_mesh(j)+x_mesh(j-1))  
c               write(*,*)'j,x_mesh_c(j)',j,x_mesh_c(j)

               if(k.eq.n+1) then
c              if(k.eq.n) then
                  goto 10
               else
                  k0=k+1
                  goto 20
               endif
            endif
          enddo !k
 20       continue     
c          write(*,*)'after 20 j,x_mesh',j,x_mesh
        enddo !j=2,n_mesh_radial_bin  (radial bin centers)
        !Note: in the above, f_equal(j-1) was defined.
        ! Still need to set f_equal(n_mesh_radial_bin)
     
        x_mesh_c(n_mesh_radial_bin)=0.5d0*
     &                              (x_mesh(n_mesh)+x_mesh(n_mesh-1))
        f_equal(n_mesh_radial_bin)=(integral_per_one_ray)/
     &                (0.5d0*(x_mesh(n_mesh)**2-x_mesh(n_mesh-1)**2))
           
 10     continue
      endif
      write(*,*)'x_mesh (bin boundaries in disk)[m]',x_mesh
      write(*,*)'x_mesh_c  (bin centers in disk)[m]',x_mesh_c
      write(*,*)'f_equal',f_equal
      !YuP: The meaning of f_equal is that it is ~1/dArea(j)
      ! where dArea(j) is the area of j-th radial bin on the disk.
      ! The power density in each radial bin on disk can be calculated as
      ! f_equal(j)/(2*pi) *powtot(1)*n_mesh_disk_angle_bin(j)
      write(*,*)
     &  '  power den[MW/m^2] in each rad bin on disk (for Ptotal=1MW):'
      do j=1,n_mesh_radial_bin
        write(*,*) (f_equal(j)/6.2831853)*1.0*n_mesh_angle(j)
      enddo
      !Compare to expected power density distribution 
      ! p_expected= func(x_mesh_c(j)) !=disk_power_distribution(x_mesh_c(j))
      !which is normalized as 
      ! INTEGRAL[0;rho_launching_disk] {p_expected(x) xdx} =1.
      ! Effectively, p_expected includes 2*pi factor;
      ! More logically, it should be 
      ! INTEGRAL[0;rho_launching_disk] {p_expected_without2pi(x) 2*pi*xdx} =1.
      !So, when we printout, we take this 2*pi factor out:
      write(*,*)'  disk_power_distribution(x_mesh_c(j))/2pi :'
      do j=1,n_mesh_radial_bin
        p_expected= func(x_mesh_c(j)) !=disk_power_distribution(x_mesh_c(j))
        write(*,*) p_expected/6.2831853  ![MW/m^2, if mult-ed by 1MW]
      enddo
   
c--------------------------------------------------------------------
c     check the total integral on the non-uniform mesh: should be 1.0
      total_integral_mesh=0.d0
      do j=1,n_mesh_radial_bin
         total_integral_mesh=total_integral_mesh+
     &                       f_equal(j)*n_mesh_angle(j)*
     &                       0.5d0*(x_mesh(j+1)**2-x_mesh(j)**2)
        write(*,*)'j,f_equal(j),n_mesh_angle(j)',
     &             j,f_equal(j),n_mesh_angle(j)
        write(*,*)'0.5d0*(x_mesh(j+1)**2-x_mesh(j)**2) = ',
     &             0.5d0*(x_mesh(j+1)**2-x_mesh(j)**2)
        write(*,*)'f_equal(j)*n_mesh_angle(j)*
     &   0.5d0*(x_mesh(j+1)**2-x_mesh(j)**2) = ',
     &   f_equal(j)*n_mesh_angle(j)*
     &   0.5d0*(x_mesh(j+1)**2-x_mesh(j)**2)
      enddo
      write(*,*)'total_intergal_mesh',total_integral_mesh !should be 1.0
      !stop

      return
      end subroutine create_equal_power_radial_mesh
      
      

      real*8 function disk_power_distribution(x)
c     the radial power distribution at the launching disk
c     x is the radius at the launching disk
c     Gaussian radial distribution at the disk
      implicit none
      include 'param.i'
      include 'cone.i'
c-----input
      real*8 x
c           sigma_launching     !these data are from cone.i
c           rho_launching_disk  !rho_launching_disk should be calculated
c                               !before usage of this function
c-----local
      real*8 const

c      disk_power_distribution=1.d0

      const=2.d0/sigma_launching_disk**2/
c     &(1.d0-(rho_launching_disk/sigma_launching_disk)**2)
     &(1.d0-dexp(-(rho_launching_disk/sigma_launching_disk)**2))

      disk_power_distribution=const*dexp(-(x/sigma_launching_disk)**2)

      return
      end function disk_power_distribution 




      subroutine disk_to_disk_rays_initial_launching_data(nray)
c--------------------------------------------------
c     calculates the initial ray's data for all rays
c     launched from the first disk to the second disk
c---------------------------------------------------
c-----input
      implicit none
      include 'param.i'
      include 'cone.i'
c-----output
      integer nray                        !the total number of the rays
c      real*8 ar_z_disk_launch(nraymax),  !space coordinates of
c     & ar_r_disk_launch(nraymax),        !launch point at the
c     & ar_phi_disk_launch(nraymax),      !first disk
c     & ar_beta_disk_launch(nraymax),     !poloidal and toroidal angles
c     & ar_alpha_disk_launch(nraymax)     !give the ray direction
c-----locals
      integer j,i,iray_l
      real*8 pi,rho_launch,
     &x_launch_disk_center,y_launch_disk_center,z_launch_disk_center,
     &cosbet,sinbet,cnz,cnx,cny,
     &x_f_center,y_f_center,z_f_center,
     &rho_0_x,rho_0_y,rho_0_z,
     &rho_perp_x,rho_perp_y,rho_perp_z,rho_perp,
     &rho_disk_x,rho_disk_y,rho_disk_z,rho_disk,
     &eta_launch,coseta,sineta,
     &rho_focus_disk_x,rho_focus_disk_y,rho_focus_disk_z,
     &x_disk_launch,y_disk_launch,z_disk_launch,
     &cnx_launch,cny_launch,cnz_launch,cn_launch,
     &r_launch,cos_phi_launch,phi_launch,
     &beta_launch,gamma_launch,alpha_launch
      real*8 disk_power_distribution
      
      real*8 f_focus_launch, rho_focus_disk_c !local


      external  disk_power_distribution
      write(*,*)'in disk_to_disk_rays_initial_launching_data(nray)'
c---------------------------------------------------
c     calculate non-uniform radial mesh on the first disk
c     disk_rho_mesh_center,disk_rho_mesh_edge
c     to get the equal power at each ray
c     disk_rho_mesh_center
c--------------------------------------------------- 
      write(*,*)'in  disk_to_disk_rays_initial_launching_data'
      write(*,*)'before create_equal_power_radial_mesh'
      write(*,*)'rho_launching_disk',rho_launching_disk
      call create_equal_power_radial_mesh(rho_launching_disk,
     &disk_power_distribution,
     &n_mesh_disk_radial_bin,n_mesh_disk_angle_bin,
     &disk_rho_mesh_edge,
     &disk_rho_mesh_center,disk_power_equal)
      write(*,*)'disk_rho_mesh_edge',disk_rho_mesh_edge
      write(*,*)'disk_rho_mesh_center',disk_rho_mesh_center
      write(*,*)'disk_power_equal',disk_power_equal
      write(*,*)'power_per_ray',power_per_ray
      pi=4.d0*datan(1.d0)
     
c--------------------------------------------------------------
c     coordinates of the the launching disk center 
c     R_launch_disk_center_vector=
c     ={x_launch_disk_center,y_launch_disk_center,z_launch_disk_center}
c--------------------------------------------------------------
      z_launch_disk_center=zst(1)
      x_launch_disk_center=rst(1)*dcos(phist(1))
      y_launch_disk_center=rst(1)*dsin(phist(1))
      write(*,*)'x_launch_disk_center',x_launch_disk_center
      write(*,*)'y_launch_disk_center',y_launch_disk_center
      write(*,*)'z_launch_disk_center',z_launch_disk_center

c--------------------------------------------------------------
c     coordinates of the refracive vector 
c     which gives the cenral ray direction
c     N_st_vector={cnx,cny,cnz}
c-------------------------------------------------------------
      cosbet=dcos(betast(1))
      sinbet=dsin(betast(1))
      cnz=sinbet
      cnx=cosbet*dcos(alfast(1)+phist(1))
      cny=cosbet*dsin(alfast(1)+phist(1))
      write(*,*)'betast(1)),alfast(1),phist(1)',
     &betast(1)*180/pi,alfast(1)*180/pi,phist(1)*180./pi
c      write(*,*)'cosbet,sinbet',cosbet,sinbet
c      write(*,*)'dcos(alfast(1)+phist(1))',dcos(alfast(1)+phist(1))
c      write(*,*)'dsin(alfast(1)+phist(1))',dsin(alfast(1)+phist(1))
      write(*,*)'cnx,cny,cnz',cnx,cny,cnz
      write(*,*)'dsqrt(cnx**2+cny**2+cnz**2)',
     &dsqrt(cnx**2+cny**2+cnz**2)
c----------------------------------------------------------------
c     calculate the coordinate of the focus disk center
c     R_f_center_vector=R_disk_center_vector+d_disk*N_st_vector
c-----------------------------------------------------------------
      x_f_center=x_launch_disk_center+d_disk*cnx
      y_f_center=y_launch_disk_center+d_disk*cny
      z_f_center=z_launch_disk_center+d_disk*cnz
      write(*,*)'d_disk',d_disk
      write(*,*)'x_f_center,y_f_center,z_f_center',
     &           x_f_center,y_f_center,z_f_center
c-----------------------------------------------------------------
c     coordinates of the unit vertor lying on the launching disk
c     perpendicular vertical Z axis={ 0*e_x, 0*e_y, 1*e_z }:
c     rho_0_vector
c     and coordinates of the unit vectot lying on the launching 
c     disk perpendicular to rho_0_vector: rho_perp_vector
c-----------------------------------------------------------------
      if ((cnx**2+cny**2).ne.0.d0) then
c-------------------------------------------------------------
c       launching disk is not perpendicuar to Z axis
c------------------------------------------------------------
c       rho_0_vector=[e_z_vector x N_st_vector]=
c                   ={-cny,cnx,0}/sqrt(cnx**2+cny**2)
c-------------------------------------------------------------
        rho_0_x=-cny/dsqrt(cnx**2+cny**2)
        rho_0_y= cnx/dsqrt(cnx**2+cny**2)
        rho_0_z= 0.d0
        write(*,*)'rho_0_x,rho_0_y,rho_0_z',
     &             rho_0_x,rho_0_y,rho_0_z
        write(*,*)'dsqrt(rho_0_x**2+rho_0_y**2+rho_0_z**2)',
     &             dsqrt(rho_0_x**2+rho_0_y**2+rho_0_z**2)
c-------------------------------------------------------------
c       rho_perp_vector=[N_st_vector x rho_0_vector]
c------------------------------------------------------------
c        rho_perp_x = - cnx*rho_0_y
c        rho_perp_y =   cny*rho_0_x
c        rho_perp_z =   cnx*rho_0_y - cny*rho_0_x
        rho_perp_x = - cnz*rho_0_y
        rho_perp_y =   cnz*rho_0_x
        rho_perp_z =   cnx*rho_0_y - cny*rho_0_x

        rho_perp = dsqrt(rho_perp_x**2+rho_perp_y**2+rho_perp_z**2)
c-------rho_perp_vector normalization
        rho_perp_x = rho_perp_x/ rho_perp
        rho_perp_y = rho_perp_y/ rho_perp
        rho_perp_z = rho_perp_z/ rho_perp
        write(*,*)'rho_perp_x,rho_perp_y,rho_perp_z',
     &             rho_perp_x,rho_perp_y,rho_perp_z
        write(*,*)'dsqrt(rho_perp_x**2+rho_perp_y**2+rho_perp_z**2)',
     &             dsqrt(rho_perp_x**2+rho_perp_y**2+rho_perp_z**2)
      else
c-------------------------------------------------------------
c       launching disk is prependicuar to Z axis
c-------------------------------------------------------------
c       rho_0_vector=[e_z_vector x [x_sr,y_st.0]]=
c                   ={-y_st,x_st,0}/sqrt(x_st**2+y_st**2)
c-------------------------------------------------------------
        rho_0_x =-y_launch_disk_center/
     &  dsqrt(x_launch_disk_center**2+y_launch_disk_center**2)
        rho_0_y = x_launch_disk_center/
     &  dsqrt(x_launch_disk_center**2+y_launch_disk_center**2)
        rho_0_z = 0.d0
c-------------------------------------------------------------
c       rho_perp_vector=[N_st_vector x rho_0_vector]
c------------------------------------------------------------
        rho_perp_x = x_launch_disk_center/
     &  dsqrt(x_launch_disk_center**2+y_launch_disk_center**2)
        rho_perp_y = y_launch_disk_center/
     &  dsqrt(x_launch_disk_center**2+y_launch_disk_center**2)
        rho_perp_z = 0.d0
      endif

      iray_l=0
      write(*,*)',n_mesh_disk_angle_bin',n_mesh_disk_angle_bin
      do j=1,n_mesh_disk_radial_bin
         write(*,*)'j=',j
         write(*,*)'n_mesh_disk_angle_bin(j)',n_mesh_disk_angle_bin(j)
         do i=1,n_mesh_disk_angle_bin(j)
c-------------------------------------------------------------
c           calculate coordinates of the vector lying on
c           the launching disk and directed from the first disk center
c           to the launch point
c------------------------------------------------------------
            iray_l=iray_l+1
            write(*,*)'iray_l',iray_l
            write(*,*)'i=',i
            if(iray_l.gt.nraymax)then
               WRITE(*,*)'iray_l.gt.nraymax'
               WRITE(*,*)'please increase nraymax in param.i'
               WRITE(*,*)'and recompile the code'            
               WRITE(*,*)'nraymax = ',nraymax               
               WRITE(*,*)'i=',i,'iray_l',iray_l
               STOP 'incone_ec.f'
            endif
               
            eta_launch=initial_azimuth_angle_degree(j)*pi/180.d0+
     &      2.d0*pi*(i-1)/n_mesh_disk_angle_bin(j)
            rho_disk=disk_rho_mesh_center(j)
            write(*,*)'eta_launch(degree),rho_disk',
     &                 eta_launch*180.d0/pi,rho_disk
            coseta=dcos(eta_launch)
            sineta=dsin(eta_launch)

            rho_disk_x=rho_disk*(rho_0_x*coseta+rho_perp_x*sineta)
            rho_disk_y=rho_disk*(rho_0_y*coseta+rho_perp_y*sineta)
            rho_disk_z=rho_disk*(rho_0_z*coseta+rho_perp_z*sineta)
            write(*,*)'rho_disk_x,rho_disk_y,rho_disk_z',
     &                 rho_disk_x,rho_disk_y,rho_disk_z
!            write(*,*)'dsqrt(rho_disk_x**2+rho_disk_y**2+rho_disk_z**2)'
!     &                ,dsqrt(rho_disk_x**2+rho_disk_y**2+rho_disk_z**2)       
            write(*,*)'rho_disk (Launching: bin centers)[m]',rho_disk   
            
            !YuP[2020-08-07] Corrected error in procedure.
            !Note that rho_disk is the radius of bin CENTER
            !where the rays are launched (from launching_disk).
            !The largest rho_disk is NOT the same as rho_launching_disk,
            !because rho_launching_disk is the outer BOUNDARY
            !of the last bin, and not the center.
            !The value of rho_disk was found by 
            !subr.create_equal_power_radial_mesh().
            !Recall that in the focus disk, the input variable
            ! rho_focus_disk also points to the outer BOUNDARY
            !of the last bin. So, in the focus disk, we also need
            !to find the corresponding values of bin CENTERS.
            ! The easiest way is to scale coordinates by 
            ! factor of 
            f_focus_launch= rho_focus_disk/rho_launching_disk !YuP
            rho_focus_disk_c= rho_disk*f_focus_launch !bin CENTER !YuP
              
            write(*,*)'rho_focus_disk_c (bin centers)',rho_focus_disk_c
c------------------------------------------------------------------
c           calculate coordinates of the vector lying on
c           the focus disk directed from the focus disk center
c           to the ray point
c           rho_focus_disk_vector=rho_disk_vector*
c                                 rho_focus_disk_c/rho_disk
c------------------------------------------------------------------    
            if (iray_l.eq.1) then
c----------------------------------------------------------------------
c              central ray is directed from the center of the launching disk
c              to the center of the second ("focus") disk
c-----------------------------------------------------------------------
               rho_focus_disk_x=0.d0
               rho_focus_disk_y=0.d0
               rho_focus_disk_z=0.d0
            else 
            write(*,*)'Launching disk: rho_disk_x,rho_disk_y,rho_disk_z'
     &                    ,rho_disk_x,rho_disk_y,rho_disk_z
!           write(*,*)'dsqrt(rho_disk_x**2+rho_disk_y**2+rho_disk_z**2)'
!     &                ,dsqrt(rho_disk_x**2+rho_disk_y**2+rho_disk_z**2) !same as rho_disk
            write(*,*)'At bin centers: rho_disk, rho_focus_disk_c',
     &                 rho_disk, rho_focus_disk_c

               rho_focus_disk_x=rho_disk_x*f_focus_launch !rho_focus_disk/rho_disk !YuP corrected
               rho_focus_disk_y=rho_disk_y*f_focus_launch !rho_focus_disk/rho_disk
               rho_focus_disk_z=rho_disk_z*f_focus_launch !rho_focus_disk/rho_disk
!         write(*,*)'rho_focus_disk_x,rho_focus_disk_y,rho_focus_disk_z'
!     &        ,rho_focus_disk_x,rho_focus_disk_y,rho_focus_disk_z
!              write(*,*)'dsqrt(rho_focus_disk_x**2+rho_focus_disk_y**2+
!     &rho_focus_disk_z**2)',
!     &dsqrt(rho_focus_disk_x**2+rho_focus_disk_y**2+rho_focus_disk_z**2)
            endif
c----------------------------------------------------------------
c           calculate coordinates of the launching point M on the
c           first disk:
c           R_M_vector=R_disk_center_vector+rho_disk_vector
c----------------------------------------------------------------
            x_disk_launch=x_launch_disk_center+rho_disk_x
            y_disk_launch=y_launch_disk_center+rho_disk_y
            z_disk_launch=z_launch_disk_center+rho_disk_z
            write(*,*)'launching point at first disk'
            write(*,*)'z_disk_launch,x_disk_launch,y_disk_launch',
     &      z_disk_launch,x_disk_launch,y_disk_launch
            write(*,*)'r=dsqrt(x_disk_launch**2+y_disk_launch**2)',
     &      dsqrt(x_disk_launch**2+y_disk_launch**2)
c----------------------------------------------------------------
c           calculate coordinates of the ray direction
c           N_launching_vector is parallel to the vector
c           (R_P_vector-R_M_vector).

c           Ray starts from the launch point M
c           R_M_vector=R_st_vector+rho_disk_vector
c           and directed to the boundary point P of the focus disk 
c           R_P_vector=R_f_center_vector+rho_focus_disk_vector
c
c           (R_P_vector-R_M_vector)=
c           =R_f_center_vector+rho_focus_disk_vector-
c            -(R_st_vector+rho_disk_vector)=
c           =(R_f_center_vector-R_st_vector)+
c            +(rho_focus_disk_vector-rho_disk_vector)=
c           =d_disk*N_st_vector+
c            +rho_disk_vector*(rho_focus_disk/rho_launching_disk-1)
c------------------------------------------------------------------
            if(iray_l.eq.1) then
c-------------central ray direction
              cnx_launch=cnx 
              cny_launch=cny
              cnz_launch=cnz
            else
c-------------non-central rays
              cnx_launch=d_disk*cnx+
     &              rho_disk_x*(rho_focus_disk_c/rho_disk-1.d0) !YuP[2020-08-07]corrected
              cny_launch=d_disk*cny+
     &              rho_disk_y*(rho_focus_disk_c/rho_disk-1.d0) !YuP[2020-08-07]corrected
              cnz_launch=d_disk*cnz+
     &              rho_disk_z*(rho_focus_disk_c/rho_disk-1.d0) !YuP[2020-08-07]corrected
              !YuP[2020-08-07]corrected: Now rho_focus_disk_c and rho_disk
              ! refer to bin CENTERs (on focus and launching disks, accordingly)
              write(*,*)'d_disk,cnz,rho_disk_z',d_disk,cnz,rho_disk_z
              write(*,*)'rho_focus_disk,rho_launching_disk',
     &                   rho_focus_disk,rho_launching_disk ! At bin BOUNDARIES
               write(*,*)'non normalized cnz_launch',cnz_launch
            endif
c-----------normalization
            cn_launch=dsqrt(cnx_launch**2+cny_launch**2+cnz_launch**2)
            cnx_launch=cnx_launch/cn_launch
            cny_launch=cny_launch/cn_launch
            cnz_launch=cnz_launch/cn_launch
            write(*,*)'cnx_launch,cny_launch,cnz_launch',
     &                cnx_launch,cny_launch,cnz_launch
            write(*,*)'dsqrt(cnx_launch**2+cny_launch**2)',
     &      dsqrt(cnx_launch**2+cny_launch**2)
c-------------------------------------------------------------------
c           transformation from (x,y,z) to phi using
c           x=r*cos(phi), y=r*sin(phi), r=sqrt(x**2+y**2)
c
c           back transformation is:
c
c                +arccos(cos(phi)) for y>=0      
c           phi=        
c                -arccos(cos(phi)) for y<0
c--------------------------------------------------------------------
            r_launch=dsqrt(x_disk_launch**2+y_disk_launch**2)
            cos_phi_launch=x_disk_launch/r_launch
            write(*,*)'x_disk_launch,y_disk_launch,r_launch',
     &      x_disk_launch,y_disk_launch,r_launch
            write(*,*)'cos_phi_launch',cos_phi_launch
            if(y_disk_launch.ge.0d0) then
               phi_launch= dacos(cos_phi_launch)
            else
               phi_launch=-dacos(cos_phi_launch)
            endif
            write(*,*)'x_disk_launch,y_disk_launch',
     &      x_disk_launch,y_disk_launch
            write(*,*)'phi_launch*180/pi',phi_launch*180/pi
c--------------------------------------------------------------------
c           Transformation from (nc,ny,nz) to alpha and beta angles    
c           and from (nx,ny,nx) to alpha,beta
c           using the formula:
c           cosbet=cos(betast)
c           sinbet=sin(betast)
c           cosalf=cos(alfast)
c           nz=sin(betast)
c           nx=cosbet*cos(alfast+phist)
c           ny=cosbet*sin(alfast+phist)
c           nr=cosbet*cosalf
c           Back transformation:
c---------------------------------------------------------------------
c           Let gamma=alphast+phistc
c----------------------------------------------------------------------
c
c                   arcsin(nz)       for cos(betast) >= 0  (A)
c           betast=
c                   pi-arcsin(nz)    for cos(betast) < 0   (B)
c
c
c           (A)If we will use betast=arcsin(nz), cos(betast) >= 0 then
c               cos(betast)= + sqrt(nx**2+ny**2) >= 0
c               it gives equations for gamma
c               nx/sqrt(nx**2+ny**2)=cos(gamma)
c               ny/sqrt(nx**2+ny**2)=sin(gamma)
c
c                     arccos(nx/sqrt(nx**2+ny**2)) for ny >= 0
c               gamma=
c                    -arccos(nx/sqrt(nx**2+ny**2)) for ny < 0
c
c           (B)If we will use betast=pi-arcsin(nz), cos(betast) <0 then
c               cos(betast)= - sqrt(nx**2+ny**2) < 0
c               it gives equations for gamma
c               -nx/sqrt(nx**2+ny**2)=cos(gamma)
c               -ny/sqrt(nx**2+ny**2)=sin(gamma)
c
c                     arccos(-nx/sqrt(nx**2+ny**2)) for ny <= 0
c               gamma=
c                    -arccos(-nx/sqrt(nx**2+ny**2)) for ny > 0
c
c           We will use (A) case
c----------------------------------------------------------------------
            write(*,*)'cnx_launch,cny_launch,cnz_launch',
     &                 cnx_launch,cny_launch,cnz_launch

            beta_launch=dasin(cnz_launch)

            if (cny_launch.ge.0d0) then
              gamma_launch=dacos(cnx_launch/
c     &        dsqrt(cnx_launch**2+cny_launch**2))
     &          dcos(beta_launch))
            else
              gamma_launch=-dacos(cnx_launch/
c     &        dsqrt(cnx_launch**2+cny_launch**2))
     &         dcos(beta_launch))
            endif

            write(*,*)'gamma_launch*180/pi',gamma_launch*180/pi

            alpha_launch=gamma_launch-phi_launch

            write(*,*)'alpha_launch*180/pi',alpha_launch*180/pi
            ar_z_disk_launch(iray_l)=z_disk_launch
            ar_r_disk_launch(iray_l)=dsqrt(x_disk_launch**2+
     &                                 y_disk_launch**2)
            ar_phi_disk_launch(iray_l)=phi_launch !radians
            ar_beta_disk_launch(iray_l)=beta_launch
            ar_alpha_disk_launch(iray_l)=alpha_launch
            write(*,*)'ar_z_disk_launch(iray_l)',
     &                 ar_z_disk_launch(iray_l)
            write(*,*)'ar_r_disk_launch(iray_l)',
     &                 ar_r_disk_launch(iray_l)
            write(*,*)'ar_phi_disk_launch(iray_l)*180/pi',
     &      ar_phi_disk_launch(iray_l)*180/pi
            write(*,*)'ar_beta_disk_launch(iray_l)*180/pi',
     &      ar_beta_disk_launch(iray_l)*180/pi
            write(*,*)'ar_alpha_disk_launch(iray_l)*180/pi',
     &      ar_alpha_disk_launch(iray_l)*180/pi
         enddo !i
      enddo !j
      write(*,*)'power_per_ray',power_per_ray
      nray=iray_l

      write(*,*)'cone_ec.f nray=',nray
c      stop 'cone_ec disk'

      do iray_l=1,nray
        powj(iray_l)=powtot(1)*1.d13/nray
        zstj(iray_l)= ar_z_disk_launch(iray_l)
        rstj(iray_l)=  ar_r_disk_launch(iray_l)
        phistj(iray_l)= ar_phi_disk_launch(iray_l)
        betaj(iray_l)= ar_beta_disk_launch(iray_l)
        alphaj(iray_l)=ar_alpha_disk_launch(iray_l)
      enddo

      return
      end subroutine disk_to_disk_rays_initial_launching_data


      subroutine disk_beam_rays_initial_launching_data(nray)
c--------------------------------------------------
c     calculates the initial ray's data for all rays
c     launched from the disk parallel to the central ray
c---------------------------------------------------
c-----input
      implicit none
      include 'param.i'
      include 'cone.i'
c-----output
      integer nray                         !the total number of the rays 
c      real*8 ar_z_disk_launch(nraymax),    !space coordinates of
c     & ar_r_disk_launch(nraymax),          !launch point at the
c     & ar_phi_disk_launch(nraymax),        !first disk
c     & ar_beta_disk_launch(nraymax),       !poloidal and toroidal angles
c     & ar_alpha_disk_launch(nraymax)       !give the ray direction
c-----locals
      integer j,i,iray_l
      real*8 pi,rho_launch,
     &x_launch_disk_center,y_launch_disk_center,z_launch_disk_center,
     &cosbet,sinbet,cnz,cnx,cny,
     &rho_0_x,rho_0_y,rho_0_z,
     &rho_perp_x,rho_perp_y,rho_perp_z,rho_perp,
     &rho_disk_x,rho_disk_y,rho_disk_z,rho_disk,
     &eta_launch,coseta,sineta,
     &x_disk_launch,y_disk_launch,z_disk_launch,
     &cnx_launch,cny_launch,cnz_launch,cn_launch,
     &r_launch,cos_phi_launch,phi_launch,
     &beta_launch,gamma_launch,alpha_launch
      real*8 disk_power_distribution


      external  disk_power_distribution
      write(*,*)'in disk_beam_rays_initial_launching_data nray',nray
c---------------------------------------------------
c     calculate non-uniform radial mesh on the launching disk
c     disk_rho_mesh_center,disk_rho_mesh_edge
c     to get the equal power at each ray
c     disk_rho_mesh_center
c--------------------------------------------------- 
      write(*,*)'in  disk_beam_rays_initial_launching_data'
      write(*,*)'before create_equal_power_radial_mesh'
      write(*,*)'rho_launching_disk',rho_launching_disk
      !YuP[2020-08-06] Was bug/typo in 1st argument:
      call create_equal_power_radial_mesh(rho_launching_disk,
     &disk_power_distribution,
     &n_mesh_disk_radial_bin,n_mesh_disk_angle_bin,
     &disk_rho_mesh_edge,
     &disk_rho_mesh_center,disk_power_equal)
      write(*,*)'disk_rho_mesh_edge',disk_rho_mesh_edge
      write(*,*)'disk_rho_mesh_center',disk_rho_mesh_center
      write(*,*)'disk_power_equal',disk_power_equal
      write(*,*)'power_per_ray',power_per_ray
      pi=4.d0*datan(1.d0)
     
c--------------------------------------------------------------
c     coordinates of the the launching disk center 
c     R_launch_disk_center_vector=
c     ={x_launch_disk_center,y_launch_disk_center,z_launch_disk_center} 
c--------------------------------------------------------------
      z_launch_disk_center=zst(1)
      x_launch_disk_center=rst(1)*dcos(phist(1))
      y_launch_disk_center=rst(1)*dsin(phist(1))
      write(*,*)'x_launch_disk_center',x_launch_disk_center
      write(*,*)'y_launch_disk_center',y_launch_disk_center
      write(*,*)'z_launch_disk_center',z_launch_disk_center

c--------------------------------------------------------------
c     coordinates of the refracive vector 
c     of the cenral ray direction
c     N_st_vector={cnx,cny,cnz}
c-------------------------------------------------------------
      cosbet=dcos(betast(1))
      sinbet=dsin(betast(1))
      cnz=sinbet
      cnx=cosbet*dcos(alfast(1)+phist(1))
      cny=cosbet*dsin(alfast(1)+phist(1))
      write(*,*)'betast(1)),alfast(1),phist(1)',
     &betast(1)*180/pi,alfast(1)*180/pi,phist(1)*180./pi
c      write(*,*)'cosbet,sinbet',cosbet,sinbet
c      write(*,*)'dcos(alfast(1)+phist(1))',dcos(alfast(1)+phist(1))
c      write(*,*)'dsin(alfast(1)+phist(1))',dsin(alfast(1)+phist(1))
      write(*,*)'cnx,cny,cnz',cnx,cny,cnz
      write(*,*)'dsqrt(cnx**2+cny**2+cnz**2)',
     &dsqrt(cnx**2+cny**2+cnz**2)

c-----------------------------------------------------------------
c     coordinates of the unit vertor lying on the launching disk
c     perpendicular vertical Z axis={ 0*e_x, 0*e_y, 1*e_z }:
c     rho_0_vector
c     and coordinates of the unit vectot lying on the launching 
c     disk perpendicular to rho_0_vector: rho_perp_vector
c-----------------------------------------------------------------
      if ((cnx**2+cny**2).ne.0.d0) then
c-------------------------------------------------------------
c       launching disk is not perpendicuar to Z axis
c------------------------------------------------------------
c       rho_0_vector=[e_z_vector x N_st_vector]=
c                   ={-cny,cnx,0}/sqrt(cnx**2+cny**2)
c-------------------------------------------------------------
        rho_0_x=-cny/dsqrt(cnx**2+cny**2)
        rho_0_y= cnx/dsqrt(cnx**2+cny**2)
        rho_0_z= 0.d0
        write(*,*)'rho_0_x,rho_0_y,rho_0_z',
     &             rho_0_x,rho_0_y,rho_0_z
        write(*,*)'dsqrt(rho_0_x**2+rho_0_y**2+rho_0_z**2)',
     &             dsqrt(rho_0_x**2+rho_0_y**2+rho_0_z**2)
c-------------------------------------------------------------
c       rho_perp_vector=[N_st_vector x rho_0_vector]
c------------------------------------------------------------
        rho_perp_x = - cnz*rho_0_y
        rho_perp_y =   cnz*rho_0_x
        rho_perp_z =   cnx*rho_0_y - cny*rho_0_x

        rho_perp = dsqrt(rho_perp_x**2+rho_perp_y**2+rho_perp_z**2)
c-------rho_perp_vector normalization
        rho_perp_x = rho_perp_x/ rho_perp
        rho_perp_y = rho_perp_y/ rho_perp
        rho_perp_z = rho_perp_z/ rho_perp
        write(*,*)'rho_perp_x,rho_perp_y,rho_perp_z',
     &             rho_perp_x,rho_perp_y,rho_perp_z
        write(*,*)'dsqrt(rho_perp_x**2+rho_perp_y**2+rho_perp_z**2)',
     &             dsqrt(rho_perp_x**2+rho_perp_y**2+rho_perp_z**2)
      else
c-------------------------------------------------------------
c       launching disk is prependicuar to Z axis
c-------------------------------------------------------------
c       rho_0_vector=[e_z_vector x [x_sr,y_st.0]]=
c                   ={-y_st,x_st,0}/sqrt(x_st**2+y_st**2)
c-------------------------------------------------------------
        rho_0_x =-y_launch_disk_center/
     &  dsqrt(x_launch_disk_center**2+y_launch_disk_center**2)
        rho_0_y = x_launch_disk_center/
     &  dsqrt(x_launch_disk_center**2+y_launch_disk_center**2)
        rho_0_z = 0.d0
c-------------------------------------------------------------
c       rho_perp_vector=[N_st_vector x rho_0_vector]
c------------------------------------------------------------
        rho_perp_x = x_launch_disk_center/
     &  dsqrt(x_launch_disk_center**2+y_launch_disk_center**2)
        rho_perp_y = y_launch_disk_center/
     &  dsqrt(x_launch_disk_center**2+y_launch_disk_center**2)
        rho_perp_z = 0.d0
      endif

      iray_l=0
      write(*,*)',n_mesh_disk_angle_bin',n_mesh_disk_angle_bin
      do j=1,n_mesh_disk_radial_bin
         write(*,*)'j=',j
         write(*,*)'n_mesh_disk_angle_bin(j)',n_mesh_disk_angle_bin(j)
         do i=1,n_mesh_disk_angle_bin(j)
c-------------------------------------------------------------
c           calculate coordinates of the vector lying on
c           the launching disk and directed from the first disk center
c           to the launch point
c------------------------------------------------------------
            iray_l=iray_l+1
            write(*,*)'iray_l',iray_l
            write(*,*)'i=',i
            if(iray_l.gt.nraymax)then
               WRITE(*,*)'iray_l.gt.nraymax'
               WRITE(*,*)'please increase nraymax in param.i'
               WRITE(*,*)'and recompile the code'            
               WRITE(*,*)'nraymax = ',nraymax               
               WRITE(*,*)'i=',i,'iray_l',iray_l
               STOP 'incone_ec.f'
            endif
               
            eta_launch=initial_azimuth_angle_degree(j)*pi/180.d0+
     &      2.d0*pi*(i-1)/n_mesh_disk_angle_bin(j)
            rho_disk=disk_rho_mesh_center(j)
            write(*,*)'eta_launch(degree),rho_disk',
     &                 eta_launch*180.d0/pi,rho_disk
            coseta=dcos(eta_launch)
            sineta=dsin(eta_launch)

            rho_disk_x=rho_disk*(rho_0_x*coseta+rho_perp_x*sineta)
            rho_disk_y=rho_disk*(rho_0_y*coseta+rho_perp_y*sineta)
            rho_disk_z=rho_disk*(rho_0_z*coseta+rho_perp_z*sineta)
            write(*,*)'rho_disk_x,rho_disk_y,rho_disk_z',
     &                 rho_disk_x,rho_disk_y,rho_disk_z
            write(*,*)'dsqrt(rho_disk_x**2+rho_disk_y**2+rho_disk_z**2)'
     &                ,dsqrt(rho_disk_x**2+rho_disk_y**2+rho_disk_z**2)       
            write(*,*)'rho_disk',rho_disk     
        
c------------------------------------------------------------------    
c           calculate coordinates of the launchig point M on the
c           first disk:
c           R_M_vector=R_disk_center_vector+rho_disk_vector
c----------------------------------------------------------------
            x_disk_launch=x_launch_disk_center+rho_disk_x
            y_disk_launch=y_launch_disk_center+rho_disk_y
            z_disk_launch=z_launch_disk_center+rho_disk_z
            write(*,*)'launching point at first disk'
            write(*,*)'z_disk_launch,x_disk_launch,y_disk_launch',
     &      z_disk_launch,x_disk_launch,y_disk_launch
            write(*,*)'r=dsqrt(x_disk_launch**2+y_disk_launch**2)',
     &      dsqrt(x_disk_launch**2+y_disk_launch**2)
c----------------------------------------------------------------
c           calculate coordinates of the ray direction
c           N_launching_vector is parallel to the central ray
c           
c           Ray starts from the launch point M
c           R_M_vector=R_st_vector+rho_disk_vector
c           and directed parallel to the central ray    
c------------------------------------------------------------------
            
            cnx_launch=cnx 
            cny_launch=cny
            cnz_launch=cnz

c-----------normalization
            cn_launch=dsqrt(cnx_launch**2+cny_launch**2+cnz_launch**2)
            cnx_launch=cnx_launch/cn_launch
            cny_launch=cny_launch/cn_launch
            cnz_launch=cnz_launch/cn_launch
            write(*,*)'cnx_launch,cny_launch,cnz_launch',
     &                cnx_launch,cny_launch,cnz_launch
            write(*,*)'dsqrt(cnx_launch**2+cny_launch**2)',
     &      dsqrt(cnx_launch**2+cny_launch**2)
c-------------------------------------------------------------------
c           transformation from (x,y,z) to phi using
c           x=r*cos(phi), y=r*sin(phi), r=sqrt(x**2+y**2)
c
c           back transformation is:
c
c                +arccos(cos(phi)) for y>=0      
c           phi=        
c                -arccos(cos(phi)) for y<0
c--------------------------------------------------------------------
            r_launch=dsqrt(x_disk_launch**2+y_disk_launch**2)
            cos_phi_launch=x_disk_launch/r_launch
            write(*,*)'x_disk_launch,y_disk_launch,r_launch',
     &      x_disk_launch,y_disk_launch,r_launch
            write(*,*)'cos_phi_launch',cos_phi_launch
            if(y_disk_launch.ge.0d0) then
               phi_launch= dacos(cos_phi_launch)
            else
               phi_launch=-dacos(cos_phi_launch)
            endif
            write(*,*)'x_disk_launch,y_disk_launch',
     &      x_disk_launch,y_disk_launch
            write(*,*)'phi_launch*180/pi',phi_launch*180/pi
c--------------------------------------------------------------------
c           Transformation from (nc,ny,nz) to alpha and beta angles         
c           and from (nx,ny,nx) to alpha,beta
c           using the formula:
c           cosbet=cos(betast)
c           sinbet=sin(betast)
c           cosalf=cos(alfast)
c           nz=sin(betast)
c           nx=cosbet*cos(alfast+phist)
c           ny=cosbet*sin(alfast+phist)
c           nr=cosbet*cosalf
c           Back transformation:
c---------------------------------------------------------------------
c           Let gamma=alphast+phistc
c----------------------------------------------------------------------
c
c                   arcsin(nz)       for cos(betast) >= 0  (A)
c           betast=
c                   pi-arcsin(nz)    for cos(betast) < 0   (B)
c
c
c           (A)If we will use betast=arcsin(nz), cos(betast) >= 0 then
c               cos(betast)= + sqrt(nx**2+ny**2) >= 0
c               it gives equations for gamma
c               nx/sqrt(nx**2+ny**2)=cos(gamma)
c               ny/sqrt(nx**2+ny**2)=sin(gamma)
c
c                     arccos(nx/sqrt(nx**2+ny**2)) for ny >= 0
c               gamma=
c                    -arccos(nx/sqrt(nx**2+ny**2)) for ny < 0
c
c           (B)If we will use betast=pi-arcsin(nz), cos(betast) <0 then
c               cos(betast)= - sqrt(nx**2+ny**2) < 0
c               it gives equations for gamma
c               -nx/sqrt(nx**2+ny**2)=cos(gamma)
c               -ny/sqrt(nx**2+ny**2)=sin(gamma)
c
c                     arccos(-nx/sqrt(nx**2+ny**2)) for ny <= 0
c               gamma=
c                    -arccos(-nx/sqrt(nx**2+ny**2)) for ny > 0
c
c           We will use (A) case
c----------------------------------------------------------------------
            write(*,*)'cnx_launch,cny_launch,cnz_launch',
     &                 cnx_launch,cny_launch,cnz_launch

            beta_launch=dasin(cnz_launch)

            if (cny_launch.ge.0d0) then
              gamma_launch=dacos(cnx_launch/dcos(beta_launch))
            else
              gamma_launch=-dacos(cnx_launch/dcos(beta_launch))
            endif

            write(*,*)'gamma_launch*180/pi',gamma_launch*180/pi

            alpha_launch=gamma_launch-phi_launch

            write(*,*)'alpha_launch*180/pi',alpha_launch*180/pi
            ar_z_disk_launch(iray_l)=z_disk_launch
            ar_r_disk_launch(iray_l)=dsqrt(x_disk_launch**2+
     &                                 y_disk_launch**2)
            ar_phi_disk_launch(iray_l)=phi_launch !radians
            ar_beta_disk_launch(iray_l)=beta_launch
            ar_alpha_disk_launch(iray_l)=alpha_launch
            write(*,*)'ar_z_disk_launch(iray_l)',
     &                 ar_z_disk_launch(iray_l)
            write(*,*)'ar_r_disk_launch(iray_l)',
     &                 ar_r_disk_launch(iray_l)
            write(*,*)'ar_phi_disk_launch(iray_l)*180/pi',
     &      ar_phi_disk_launch(iray_l)*180/pi
            write(*,*)'ar_beta_disk_launch(iray_l)*180/pi',
     &      ar_beta_disk_launch(iray_l)*180/pi
            write(*,*)'ar_alpha_disk_launch(iray_l)*180/pi',
     &      ar_alpha_disk_launch(iray_l)*180/pi
         enddo !i
      enddo !j
      write(*,*)'power_per_ray',power_per_ray
      nray=iray_l

      write(*,*)'cone_ec.f nray=',nray
c      stop 'cone_ec disk_beam'

      do iray_l=1,nray
        powj(iray_l)=powtot(1)*1.d13/nray
        zstj(iray_l)= ar_z_disk_launch(iray_l)
        rstj(iray_l)=  ar_r_disk_launch(iray_l)
        phistj(iray_l)= ar_phi_disk_launch(iray_l)
        betaj(iray_l)= ar_beta_disk_launch(iray_l)
        alphaj(iray_l)=ar_alpha_disk_launch(iray_l)
        write(*,*) 'iray_l,betaj(iray_l),alphaj(iray_l)',
     &              iray_l,betaj(iray_l),alphaj(iray_l)
        write(*,*)'betast(1),alfast(1)',betast(1),alfast(1)
      enddo

      return
      end subroutine disk_beam_rays_initial_launching_data
