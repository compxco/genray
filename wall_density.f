      subroutine create_wall_density
      implicit none

      return
      end

      subroutine create_fine_mesh_for_chamber_wall_coordinates 
c     It uses input wall arrays r_wall(n_wall_a),z_wall(n_wall_a)
c     from one_nml.i
c
c     Using linear approximation it creates new arrays 
c     r_wall_add(n_wall_add_a),z_wall_add(n_wall_add_a)
c     they will be in  fourb.i
c     The distance between additional point and old mesh points
c     is  h_add_wall*k


      implicit none
      include 'param.i'
      include 'one.i'
      include 'fourb.i'

c-----locals
      integer i,j,k,ki,
     &n_wall_all   !number of points in r_wall_add and z_wall_add arrays

      real*8 
     &length_wall,               !poloidal wall length    
     &delta_r_wall,delta_z_wall, !r,z lengths of the wall element
     &delta_length_wall

      length_wall=0.d0
      do i=1,n_wall-1
         length_wall=length_wall+dsqrt((r_wall(i+1)-r_wall(i))**2+
     &                         (z_wall(i+1)-z_wall(i))**2)
      enddo
      
      j=0
      do i=1,n_wall-1
         delta_r_wall=r_wall(i+1)-r_wall(i)
         delta_z_wall=r_wall(i+1)-r_wall(i)
         delta_length_wall=dsqrt(delta_r_wall**2+delta_z_wall**2)

         j=j+1
         r_wall_add(j)=r_wall(i)
         z_wall_add(j)=z_wall(i)

         if (delta_length_wall.gt.h_add_wall) then
c-----------add additional wall points
            j=j+1
          
            ki= delta_length_wall/h_add_wall ! number of additional wall points
                                             ! at the interval(i,i+1)
            do k=1,ki
               j=j+1
               if (j.gt.n_wall_add_a) then
                  write(*,*)'in create_fine_wall_mesh'
                  write(*,*)'j.gt.n_wall_add_a'
                  write(*,*)'j,n_wall_ad_a',j,n_wall_add_a
                  write(*,*)'Please increase n_wall_add_a'
                  write(*,*)'in param.i and recompile the code' 
                  stop 'create_seldom_mesh_for_chamber_wall_coordinates'
               endif 
               r_wall_add(j)=r_wall(i)+
     &                 delta_r_wall*h_add_wall*k/delta_length_wall
               z_wall_add(j)=z_wall(i)+
     &                 delta_z_wall*h_add_wall*k/delta_length_wall
            enddo
          endif

      enddo
      
      return
      end

    
