      subroutine plot_cold_n_perp_omega_npar(z,r,phi,cnpar,
     &ratio_freq_min,ratio_freq_max,n_plot_freq)
c     
      implicit none

      include 'param.i'
      include 'one.i'
c-----input
      real*8 
     &z,r,phi,         ! space coordinates
     & ratio_freq_min, ! frequency_min/frqncy
     & ratio_freq_max,  ! frequency_max/frqncy
     & cnpar           ! N_parallel

      integer n_plot_freq !number of points of frequency mesh
c-----locals
      real*8 freq_loc,freq_min,freq_max,step,df, cnpar_loc,
     &eps,g,eta,det,yi,xi,a2,a0

      real*8, dimension(1:nbulk) :: v_loc,w_loc
      real*8, dimension(1:n_plot_freq) :: freq_ar,n_perp2_p,n_perp2_m
      integer n,i
c-----external
      real*8 b,x,y

      freq_min=ratio_freq_min*frqncy
      freq_max=ratio_freq_max*frqncy

      step=(freq_max-freq_min)/(n_plot_freq-1)

      bmod=b(z,r,phi)
      do i=1,nbulk
	  v_loc(i)=v(i)
	  w_loc(i)=w(i)
      enddo

      do n=1,n_plot_freq

         freq_loc=  freq_min+step*(n-1)  
         freq_ar(n)= freq_loc
         df=frqncy/freq_loc 

         do i=1,nbulk
       	    v(i)=v_loc(i)*df* df
	    w(i)=w_loc(i)*df
         enddo !i
       
         cnpar_loc=cnpar*df

         eps=1.d0
         g=0.d0
         eta=1.d0

         do i=1,nbulk
            xi=x(z,r,phi,i)
            yi=y(z,r,phi,i)

       	    eps=eps-xi/(1.d0-yi**2)
            eta=eta-xi

            if (i.eq.1) then
               g=g+xi*yi/(1.d0-yi**2)
            else
	       g=g-xi*yi/(1.d0-yi**2)
            endif  
         enddo !i

c--------dispersion relation
c        eps*N_perp**4-a2*N_perp**2+a0=0

         a2=(eps+eta)*(eps-cnpar_loc**2)-g*g
         a0=eta*((eps-cnpar_loc**2)**2-g*g)

         det=a2*a2-4.d0*eps*a0
  
         if (det.ge.0) then
           n_perp2_p(n)=(a2+dsqrt(det))/(2.d0*eps)
           n_perp2_p(n)=(a2-dsqrt(det))/(2.d0*eps)
         else
           n_perp2_p(n)=0.d0
           n_perp2_p(n)=0.d0
         endif

      enddo !n

      do i=1,nbulk
       	 v(i)=v_loc(i)
	 w(i)=w_loc(i)
      enddo !i   

       call plot1dt(freq_ar,n_perp2_p,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','N_perp**2_p')

      call plot1dt(freq_ar,n_perp2_m,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','N_perp**2_m')

      return
      end

      
