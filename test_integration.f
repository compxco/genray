      program test
      implicit none
      integer n,n_mesh_radial_bin,n_mesh
c     parameter (n=1000)
      parameter (n_mesh_radial_bin=3,n_mesh=n_mesh_radial_bin+1)

      double precision b
      integer n_mesh_angle(n_mesh_radial_bin)
      double precision x_mesh(n_mesh),x_mesh_c(n_mesh_radial_bin),
     &f_equal(n_mesh_radial_bin)


      double precision func
      external func

      integer i
      double precision analytic_integral,dif,sum_dif,
     &total_integral_analytical

      n_mesh_angle(1)=1
      n_mesh_angle(2)=3
      n_mesh_angle(3)=5 
      b=1.d0

c      call radial_integral_at_uniform_mesh(func,b,n,x,xc,
c     &ar_integral,total_integral)

c      write(*,*)'n',n
c      do i=1,n
c        write(*,*)'i,x(i)',i,x(i)
c      enddo

c      do i=1,n-1
c        write(*,*)'i,xc(i)',i,xc(i)
c      enddo

c      sum_dif=0.d0
c      do i=1,n
cc       analytic_integral=x(i)**2/2.d0 !func=1
cc        analytic_integral=x(i)**3/3.d0 !func=x
c        analytic_integral=x(i)**2/2.d0-x(i)**4/4.d0 !func=1-x**2
c        dif=analytic_integral-ar_integral(i)
c        write(*,*)'i,x(i),ar_integral(i),analytic_integral,dif',
c     &             i,x(i),ar_integral(i),analytic_integral,dif
c        sum_dif=sum_dif+dif
c      enddo
c      total_integral_analytical=0.5d0-0.25d0 !func=1-x**2
c      write(*,*)'total_integral,sum_dif',total_integral,sum_dif
c      write(*,*)'total_integral_analytical',total_integral_analytical
      
      write(*,*)'before  create_equal_radial_mesh'
      call create_equal_power_radial_mesh(b,func,n_mesh_radial_bin,
     &n_mesh_angle,x_mesh,x_mesh_c,f_equal)
      write(*,*)'after  create_equal_radial_mesh'

      return
      end


      subroutine radial_integral_at_uniform_mesh(func,b,n,
     &x,xc,ar_integral,total_integral)
c-----calculates integral{0,x}f(x)xdx and puts it in ar_integral
c     0=<  x <b
c     using the formula
c     Sum{k=2,j-1}[f(x(k))*((0.5*(x(k+1)+x(k)))**2-(0.5*(x(k)+x(k-1)))**2)],
c     for  j=2,n-1
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
      real*8 b 
      integer n 
      double precision func

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

c      write(*,*)'radial_integral_at_uniform_mesh n, x',n,x

      xc(1)=0.d0
      do j=2,n-1         
        xc(j)=step*(j-1)
      enddo
c      write(*,*)'radial_integral_at_uniform_mesh xc',xc

      do j=1,n
        ar_integral(j)=0.d0
      enddo

      do j=2,n
        ar_integral(j)=ar_integral(j-1)+func(xc(j-1))*
     &                 0.5d0*(x(j)**2-x(j-1)**2)
      enddo
      total_integral=ar_integral(n)

c      do j=2,n
c        write(*,*)'j,xc(j-1),x(j),ar_integral(j)',
c     *             j,xc(j-1),x(j),ar_integral(j)
c      enddo
c      write(*,*)'total_integral',total_integral

c-----normalization to get the total integral =1
      do j=1,n
        ar_integral(j)=ar_integral(j)/total_integral
      enddo
      total_integral=ar_integral(n)

      return
      end


      double precision function func(x)
      implicit none
c-----input
      double precision x

      func=1.d0
c      func=x
c      func=1.d0-x*x

      return
      end      

      subroutine create_equal_power_radial_mesh(b,func,
     &n_mesh_radial_bin,n_mesh_angle,x_mesh,x_mesh_c,f_equal)

c-----creates X mesh and array for the function func:
c     n_mesh_radial_bin is the number of radial bins 
c     n_mesh=n_mesh_radial_bin + 1 is the number of radius point in which
c              the integral is calculated     
c     array x_mesh_c(n_mesh_radial_bin)  (0<  x_mesh_c <b) 
c     array x_mesh(n_mesh) (0=<  x_mesh =<b)
c     array f_equal(n n_mesh_radial_bin) of the function func
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
c     In each radial bin (j= n_mesh_radial_bin1,n_mesh-1) there will be n_mesh_angle(j)
c     rays with different angles. In the first radial bin there will be 
c     one ray n_mesh_angle(1)=1
c
c     The total number of rays n_disk_rays
c     n_disk_rays = Sum{j=1.n_mesh_radial_bin}
c
c     The mesh and the function array will be calculate to give
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
      real*8 x(n),xc(n-1),ar_integral(n),total_integral

      real*8  total_integral_mesh,total_integral_per_one_ray,
     &integral_per_one_ray

c------------------------------------------------------------
c     calculate the integral from the function func
c     ar_integral(j) j=1,...,n
c     at the equispaced mesh x()

      write(*,*)'create_equal_radial_mesh b,n',b,n

c------------------------------------------------------------
      call radial_integral_at_uniform_mesh(func,b,n,x,xc,
     &ar_integral,total_integral)
c------------------------------------------------------------

      write(*,*)'create_equal_mesh: ar_integral',ar_integral
      write(*,*)'create_equal_mesh: x', x

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
        write(*,*)'j,n_mesh_angle(j)',j,n_mesh_angle(j)
      enddo
      write(*,*)'n_disk_rays',n_disk_rays

      write(*,*)'total_integral',total_integral

      integral_per_one_ray=total_integral/n_disk_rays

      write(*,*)'integral_per_one_ray',integral_per_one_ray
      
      j_rays=0
      do j=2,n_mesh_radial_bin

c-------number of rays
        do i=1,n_mesh_angle(j-1)
           j_rays=j_rays+1
        enddo

        integral_value=j_rays*integral_per_one_ray

        write(*,*)'j,integral_value',j,integral_value
        do k=k0,n+1

           write(*,*)'k0,k,integral_value,ar_integral(k)',
     &                k0,k,integral_value,ar_integral(k)

          if(integral_value.le.ar_integral(k)) then
c---------------------------------------------------------------------
c           The equation to find xmesh(j)
c           (integral_value-ar_integral(k-1))/(x_mesh(j)**2-x(k-1)**2)=
c           =(ar_integral(k)-ar_integral(k-1))/(x_(k)**2-x(k-1)**2)
c-----------------------------------------------------------------------

             write(*,*)'k,x(k-1),x(k)',k,x(k-1),x(k)

            x_mesh(j)=dsqrt(x(k-1)**2+(integral_value-ar_integral(k-1))*
     &      (x(k)**2-x(k-1)**2)/(ar_integral(k)-ar_integral(k-1)))

            write(*,*)'j,x_mesh(j)',j,x_mesh(j)
            write(*,*)'x_mesh',x_mesh            
           
            f_equal(j-1)=integral_per_one_ray/
     &                   (0.5d0*(x_mesh(j)**2-x_mesh(j-1)**2))

            write(*,*)'j,integral_per_one_ray',
     &                 j,integral_per_one_ray
            write(*,*)'x_mesh(j-1),x_mesh(j),f_equal(j-1)',
     &                 x_mesh(j-1),x_mesh(j),f_equal(j-1)
            if(j.gt.2) x_mesh_c(j-1)=0.5d0*(x_mesh(j)+x_mesh(j-1))  
            write(*,*)'j,x_mesh_c(j)',j,x_mesh_c(j)

            if(k.eq.n+1) then
c            if(k.eq.n) then
              goto 10
            else
              k0=k+1
              goto 20
            endif
          endif
        enddo !k
 20     continue     
      write(*,*)'after 20 j,x_mesh',j,x_mesh
      enddo !j
     
      x_mesh_c(n_mesh_radial_bin)=0.5d0*
     &                              (x_mesh(n_mesh)+x_mesh(n_mesh-1))
      f_equal(n_mesh_radial_bin)=(integral_per_one_ray)/
     &                (0.5d0*(x_mesh(n_mesh)**2-x_mesh(n_mesh-1)**2))
           
 10   continue

      write(*,*)'x_mesh',x_mesh
      write(*,*)'x_mesh_c',x_mesh_c
      write(*,*)'f_equal',f_equal
    
   
c--------------------------------------------------------------------
c     check the total integral on the non-uniform mesh
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
     &0.5d0*(x_mesh(j+1)**2-x_mesh(j)**2) = ',
     &f_equal(j)*n_mesh_angle(j)*
     &0.5d0*(x_mesh(j+1)**2-x_mesh(j)**2)

      enddo
      write(*,*)'total_intergal_mesh',total_integral_mesh

      return
      end

    
