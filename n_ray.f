      subroutine wave_ray_normal_surface(z,r,phi,cnpar_ray,cnper_ray,
     &t,n_gam,ib)
c-----calculates and creates plot of the wave normal surface 
c     and the wave normal surface for the ray refractive index (1/n_r)
c     for cold plasma dispersion
c     in z,r,phi point
c     
      implicit none 
c-----input
      real*8 z,r,phi, ! space coordinates 
     &cnpar_ray,cnper_ray, ! parallel and perpendicular refractive index
                           ! components along the ray
     &t               ! poloidal distance of the ray point
      integer n_gam,  ! the number of angle points
     &ib ! 
c-----local
      integer n_gam_a !maximal value of n_gam
      parameter (n_gam_a=500)
      real*8 gam,     !the angle between B_field and k_vector
     & gam_ar(n_gam_a),! array of angles
     & step,           ! angle step      
     & nperp_1,nperp_2, ! roots of the dispersion relation
     & ns_1_ar(n_gam_a),! first root array
     & ns_2_ar(n_gam_a),! second root array
     & x_ar_1(n_gam_a),z_ar_1(n_gam_a),!projections of the 
     & x_ar_2(n_gam_a),z_ar_2(n_gam_a),!the wave normal vector
     & x_ar_3(n_gam_a,2),z_ar_3(n_gam_a,2),!the wave normal vector
     & x_Nray_ar_1(n_gam_a),z_Nray_ar_1(n_gam_a),!projections of the
                                                 !the wave ray normal vector
     & x_Nray_ar_2(n_gam_a),z_Nray_ar_2(n_gam_a),!projections of the
                                                 !the wave ray normal vector
     & pi,
     & ad,bd,cd,         ! dispersion relation coefficients
     & d4,d2,d0,det,cn1,ds,dc,ds2,dc2,ds4,
     & xe,xi,ye,yb,dele,delib,dn,
     & omegpedce,omegdce,cnpar_loc,cnper_loc,
     & N_ray, !ray refractive index
     & cnpar_Nray,cnper_Nray,gam_Nray     

      integer n_param                    !number of parameters at figure
      character*(20) name_param(4)       !names  of parameters at figure
      real param(4)                      !values of parameters at figure 

c-----for Drawing Markers:
      integer n_marker        ! number of points to be marked,
      parameter  (n_marker=1)   
      real  x_marker(n_marker),y_marker(n_marker)
      integer n_symbol_marker
      real*8 gam_ray
c-----externals
      real*8 x,y,rrind

      integer i,i_sign    

      i_sign=1

      if (n_gam.gt.n_gam_a) then
         write(*,*)'in wave_normal_surface'
         write(*,*)'n_gam > n_gam_a'
         write(*,*)'n_gam=',n_gam,'n_gam_a=',n_gam_a
         write(*,*)'Please change parameter n_gam_a'
         write(*,*)' in wave_normal_surface and recompile the code'
         stop 'in wave_normal_surface'
      endif


      do i=1,n_gam
         ns_1_ar(i)=0.d0
         ns_2_ar(i)=0.d0
         x_ar_1(i)=0.d0
         x_ar_2(i)=0.d0
         z_ar_1(i)=0.d0
         z_ar_2(i)=0.d0
         x_Nray_ar_1(i)=0.d0
         x_Nray_ar_2(i)=0.d0
         z_Nray_ar_1(i)=0.d0
         z_Nray_ar_2(i)=0.d0
      enddo     

      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      omegpedce=dsqrt(xe/dabs(ye))      !omega_pe/omega_ce
      omegdce=1.d0/dabs(ye)             !omega/omega_ce

      pi=4.d0*datan(1.d0)

      step=2.d0*pi/dfloat(n_gam-1)
      
      do i=1,n_gam
        gam=step*(i-1)
        gam_ar(i)=gam
        ds=dsin(gam)
        dc=dcos(gam) 
        ds2=ds*ds
        dc2=dc*dc
        ds4=ds2*ds2
       
        call abc(z,r,phi,ds2,dc2,ad,bd,cd) 
c-------dispersion relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
        d4=ad
        d2=bd
        d0=cd
        det=d2*d2-4.d0*d4*d0
        cn1=dsqrt(d2*d2-4.d0*d4*d0)

        if (ib.eq.1) then
c----------ib.eq.1 electron resonance condition may be
           xe=x(z,r,phi,1)
           ye=y(z,r,phi,1)
c           write(*,*)'wave_normal_surface xe,ye',xe,ye	 
           dele=1.d0-ye
           if (dele.lt.0.d0) i_sign=-1
c---------------------------------------------------------------------
c          calculate roots of the dispersion
c          relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c---------------------------------------------------------------------
c           write(*,*)'(-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
c           write(*,*)'(-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)
           ns_1_ar(i)=(-d2+cn1)/(2.d0*d4)
           ns_2_ar(i)=(-d2-cn1)/(2.d0*d4)         
	  
        endif ! ib.eq.1

        if (ib.gt.1) then
c----------ib.gt.1 iones resonance condition may be
           yb=y(z,r,phi,ib)
	   delib=1.d0-yb     
           if (delib.lt.0.d0) i_sign=-1
c--------------------------------------------------------------------
c          calculate roots of the dispersion
c          relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c--------------------------------------------------------------------
           cn1=dsqrt(d2*d2-4.d0*d4*d0)
c           write(*,*)'(-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
c           write(*,*)'(-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)       
	   ns_1_ar(i)=(-d2+cn1)/(2.d0*d4)
           ns_2_ar(i)=(-d2-cn1)/(2.d0*d4)     
        endif ! ib.gt.1

c        write(*,*)'i, ns_1_ar(i)',i, ns_1_ar(i)
        if (ns_1_ar(i).gt.0.d0) then
           dn=1.d0/dsqrt(ns_1_ar(i))
           z_ar_1 (i)=dc*dn
           x_ar_1 (i)=ds*dn
c------------------------------------------------------------------------
c          calculate ray refractive index: cnray
c--------------------------------------------------------------
           cnpar_loc=dc/dn
           cnper_loc=ds/dn
           N_ray=rrind(cnpar_loc,cnper_loc,omegdce,omegpedce)
           z_Nray_ar_1 (i)=dc/N_ray
           x_Nray_ar_1 (i)=ds/N_ray
        endif 

c        write(*,*)'i,x_ar_1(i),z_ar_1(i)',i,x_ar_1(i),z_ar_1(i)

c        write(*,*)'i, ns_2_ar(i)',i, ns_2_ar(i)
        if(ns_2_ar(i).gt.0.d0) then
           dn=1.d0/dsqrt(ns_2_ar(i))
           z_ar_2 (i)=dc*dn
           x_ar_2 (i)=ds*dn
c------------------------------------------------------------------------
c          calculate ray refractive index: cnray
c--------------------------------------------------------------
           cnpar_loc=dc/dn
           cnper_loc=ds/dn
           N_ray=rrind(cnpar_loc,cnper_loc,omegdce,omegpedce)
           z_Nray_ar_2 (i)=dc/N_ray
           x_Nray_ar_2 (i)=ds/N_ray
        endif
c        write(*,*)'i,x_ar_2(i),z_ar_2(i)',i,x_ar_2(i),z_ar_2(i)

      enddo !i=1,n_gam
 
      n_param=4
      name_param(1)='r [m]'
      name_param(2)='t [m]'
c----------------------------------------------------------------------------
c     n**2 at the dipersion curve VS gam
c----------------------------------------------------------------------------
      name_param(3)='n**2 at ray'
      name_param(4)='gamma  at ray'
      param(1)=r
      param(2)=t
      goto10 
      param(3)=(cnpar_ray**2+cnper_ray**2)
      gam_ray=dacos(cnpar_ray/dsqrt(cnpar_ray**2+cnper_ray**2))
      param(4)=gam_ray

c-----with marker

      n_symbol_marker=9
      gam_ray=dacos(cnpar_ray/dsqrt(cnpar_ray**2+cnper_ray**2))
      x_marker(1)=gam_ray
      y_marker(1)=param(3)
        

      call plot1dt_marker_param(gam_ar,ns_1_ar,0,0,n_gam,2,
     & n_marker,x_marker,y_marker,n_symbol_marker,    
     &'linlin',0.d0,0.d0,
     &'n**2 1','gam','ns 1',
     &n_param,name_param,param)

      call plot1dt_marker_param(gam_ar,ns_2_ar,0,0,n_gam,2,
     & n_marker,x_marker,y_marker,n_symbol_marker,    
     &'linlin',0.d0,0.d0,
     &'n**2 2','gam','ns 2',
     &n_param,name_param,param)

c--------------------------------------------------------------------
c     wave normals
c--------------------------------------------------------------------
      name_param(3)='n_perp/n**2 at ray'
      name_param(4)='n_par/n**2  at ray'
      param(3)=cnper_ray/(cnpar_ray**2+cnper_ray**2)
      param(4)=cnpar_ray/(cnpar_ray**2+cnper_ray**2)

c-----with marker

      n_symbol_marker=9
      x_marker(1)=param(3)
      y_marker(1)=param(4)

      call plot1dt_marker_param(x_ar_1,z_ar_1,0,0,n_gam,2,
     & n_marker,x_marker,y_marker,n_symbol_marker,    
     &'linlin',0.d0,0.d0,
     &'wave normal n^/n**2 1 root','n_perp/n**2','n_par/n**2',
     &n_param,name_param,param)

      call plot1dt_marker_param(x_ar_2,z_ar_2,0,0,n_gam,2,
     & n_marker,x_marker,y_marker,n_symbol_marker,    
     &'linlin',0.d0,0.d0,
     &'wave normal n^/n**2 2 root','n_perp/n**2','n_par/n**2',
     &n_param,name_param,param)

c--------------------------------------------------------------------
c     wave ray normals
c--------------------------------------------------------------------
 10   continue
      name_param(3)='n_ray_perp/n_ray**2 at ray'
      name_param(4)='n_ray_par/n_ray**2  at ray'
      param(3)=cnper_ray/(cnpar_ray**2+cnper_ray**2)
      param(4)=cnpar_ray/(cnpar_ray**2+cnper_ray**2)
c      call plot1dt_param(x_Nray_ar_1,z_Nray_ar_1,0,0,n_gam,1,'linlin',
c     & 0.d0,0.d0,
c     &'wave ray normal n^/n**2 1 root','n_ray_perp/n_ray**2',
c     &'n_ray_par/n_ray**2',
c     &n_param,name_param,param)
                               
c      call plot1dt_param(x_Nray_ar_2,z_Nray_ar_2,0,0,n_gam,1,'linlin',
c     & 0.d0,0.d0,
c     &'wave ray normal n^/n**2 1 root','n_ray_perp/n_ray**2',
c     &'n_ray_par/n_ray**2',
c     &n_param,name_param,param)
c-----with marker

      n_symbol_marker=9
      x_marker(1)=param(3)
      y_marker(1)=param(4)
      N_ray=rrind(cnpar_ray,cnper_ray,omegdce,omegpedce)
      cnpar_Nray=dc/N_ray
      cnper_Nray=ds/N_ray
      gam_Nray=dacos(cnpar_Nray/dsqrt(cnpar_Nray**2+cnper_Nray**2))

      x_marker(1)=gam_Nray
      y_marker(1)=1.d0/(cnpar_Nray**2+cnper_Nray**2)
        
      call plot1dt_marker_param(x_Nray_ar_1,z_Nray_ar_1,0,0,n_gam,2,
     & n_marker,x_marker,y_marker,n_symbol_marker,    
     &'linlin',0.d0,0.d0,
     &'wave ray normal n_ray^/n_ray**2 1 root','n_ray_perp/n**2',
     &'n_ray_par/n**2',
     &n_param,name_param,param)

      call plot1dt_marker_param(x_Nray_ar_2,z_Nray_ar_2,0,0,n_gam,2,
     & n_marker,x_marker,y_marker,n_symbol_marker,    
     &'linlin',0.d0,0.d0,
     &'wave ray normal n_ray^/n_ray**2 2 root','n_ray_perp/n_ray**2',
     &'n_ray_par/n_ray**2',
     &n_param,name_param,param)

      return
      end
