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
      !IMPLICIT double precision (a-h,o-z)
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


c        *********************tempe_old ***********************
c        *                        -                           *
c        * this function calculates the electron temperature  *
c        * profile as a function  from rho       	      *
c        * rho=rhopsi(psi) 				      *
c        * it calculates in function b()		      *
c        * rho is  in common 'one'			      *
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
	double precision
     1function tempe_old(z,r,phi,i)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
c------------------------------------------------------------------
c     te0,teb,rn1te,rn2te  are in common 'one'
c     they were read in subroutine dinit from the file genray.in
c           rho is  in common 'one'
c------------------------------------------------------------------
cc     analytic form of the i component temperature profile
c     if(idens.eq.0) then
c	  x1=(rho)**rn1te(i)
c          tempe=(te0(i)-teb(i))*(1.d0-x1)**rn2te(i)+teb(i)
cc                         spline form of the electron temperature profile
c      else
	  tempe_old=temperho(rho,i)
c      endif
cc      write(*,*)'in dense.for rho,i,tempe=',rho,i,tempe

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
      double precision FUNCTION zeff_old(z,r,phi)
      IMPLICIT double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'

c------------------------------------------------------------------
c      zeff0,zeffb,rn1zeff,rn2zeff  are in common 'one'
c      they were read in subroutine dinit_mr from the file genray.in
c------------------------------------------------------------------
c                                analytic form of the Z_eff profile
c      if((idens.eq.0).and.(izeff.eq.1)) then
c	  x1=(rho)**rn1zeff
c          zeff=(zeff0-zeffb)*(1.d0-x1)**rn2zeff+zeffb
cc                                spline form of the Z_eff profile
c      else
	  zeff_old=zeffrho(rho)
c      endif
cc      write(*,*)'in dense.for zeff=',zeff
      return
      END


      subroutine sigma_edge_n_theta_pol(theta_pol_radian,
     &sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     &d2_sigma_edge_n_d2_theta_pol)
c-----calculates sigma_edge_n and its derivative by poloidal angle
c     theta_pol [radians]
c     
c     It will use the analytical formula or spline approximation
c     depends on the variable i_edge_dens_anal given in edge_prof_nml.i
c     i_edge_dens_anal, ! to choose analytic formula or spline:
c                       ! =0 sigmedgn=const which was set in 
c                       !    the input namelist/write/
c                       ! =1 analytic formula for sigmedgn(theta_pol)
c                       ! =2 table data for  sigmedgn(theta_pol)
      implicit none
      include 'param.i'
      include 'edge_prof.i'  
      include 'one.i'
c-----input
      real*8 
     &theta_pol_radian       !poloidal angle [radians]
c-----output
      real*8
     &sigma_edge_n,
     &d_sigma_edge_n_d_theta_pol,  !first derivate
     &d2_sigma_edge_n_d2_theta_pol !second derivative
c-----locals 
     
      real*8, dimension(3*n_pol_edge_dens+1) :: workk
      real*8 tabl(3),exp1,exp2,del1,del2,
     & d_arg_exp1_d_theta, d_arg_exp2_d_theta,
     & d2_arg_exp1_d2_theta,d2_arg_exp2_d2_theta

      integer i,iop(2),itabl(3)

      logical first 
      data  first /.true./
      save first
  
      pi=4.d0*datan(1.d0)
     
c      write(*,*)'sigma_edge_n_theta_pol theta_pol_radian',
c     &           theta_pol_radian

c      write(*,*)'i_edge_dens_anal',i_edge_dens_anal
c      write(*,*)'first',first

      if ( n_pol_edge_dens.eq.1) then 
         i_edge_dens_anal=1
          sigma_edgen_0= sigmedgn_ar(1)
      endif

cSAP090222
      if (i_edge_dens_anal.eq.0) then
c        sigma_edge_n=sigmedgn=const set in the input namelist /edge_prof_nml/
         sigma_edge_n=sigmedgn
         d_sigma_edge_n_d_theta_pol=0.d0
         d2_sigma_edge_n_d2_theta_pol=0.d0
         goto 10
      endif

      if (i_edge_dens_anal.eq.1) then
c--------analytical formula---------------
c        
c         sigma_edge_n=sigma_edgen_0
c         d_sigma_edge_n_d_theta_pol=0.d0
c         d2_sigma_edge_n_d2_theta_pol=0.d0

         if (first) then
c-----------transform from degree to radian
c           write(*,*)'dense.f pi',pi  
c           write(*,*)'dense.f sigma_theta_pol_edge_1_degree',
c     &                        sigma_theta_pol_edge_1_degree

           sigma_theta_pol_edge_1_radian=
     &           sigma_theta_pol_edge_1_degree*pi/180.d0    
 
c           write(*,*)'dense.f sigma_theta_pol_edge_1_radian',
c     &                        sigma_theta_pol_edge_1_radian

           sigma_theta_pol_edge_2_radian=
     &           sigma_theta_pol_edge_2_degree*pi/180.d0

           theta_pol_edge_1_radian=theta_pol_edge_1_degree*pi/180.d0
           theta_pol_edge_2_radian=theta_pol_edge_2_degree*pi/180.d0

c          write(*,*)'dense.f sigma_theta_pol_edge_1_radian',
c     &                        sigma_theta_pol_edge_1_radian

           first=.false.  
         endif

c         write(*,*)'dense.f sigma_edge_n_theta_pol'
c         write(*,*)'theta_pol_radian',theta_pol_radian
c         write(*,*)'theta_pol_edge_1_radian',theta_pol_edge_1_radian
c         write(*,*)'theta_pol_edge_2_radian',theta_pol_edge_2_radian
c         write(*,*)'sigma_theta_pol_edge_1_radian',
c     &              sigma_theta_pol_edge_1_radian
c         write(*,*)'sigma_theta_pol_edge_2_radian',
c     &              sigma_theta_pol_edge_2_radian

         exp1=dexp(-((theta_pol_radian-theta_pol_edge_1_radian)/
     &                sigma_theta_pol_edge_1_radian)**2)

         exp2=dexp(-((theta_pol_radian-theta_pol_edge_2_radian)/
     &                sigma_theta_pol_edge_2_radian)**2)
 
c        write(*,*)'exp1,exp2',exp1,exp2

c         write(*,*)'sigma_edgen_max_1',sigma_edgen_max_1
c         write(*,*)'sigma_edgen_max_2',sigma_edgen_max_2
c         write(*,*)'sigma_edgen_min',sigma_edgen_min

         del1=sigma_edgen_1-sigma_edgen_0

         del2=sigma_edgen_2-sigma_edgen_0

c         write(*,*)'del1,del2',del1,del2

         d_arg_exp1_d_theta=
     &         -2.d0*(theta_pol_radian-theta_pol_edge_1_radian)/
     &         sigma_theta_pol_edge_1_radian**2

         d_arg_exp2_d_theta=
     &         -2.d0*(theta_pol_radian-theta_pol_edge_2_radian)/
     &         sigma_theta_pol_edge_2_radian**2

         d2_arg_exp1_d2_theta=-2.d0/sigma_theta_pol_edge_1_radian**2
         d2_arg_exp2_d2_theta=-2.d0/sigma_theta_pol_edge_2_radian**2
       
c         write(*,*)'sigmedgn',sigmedgn

         sigma_edge_n=sigma_edgen_0+del1*exp1+del2*exp2

c         write(*,*)'sigma_edge_n',sigma_edge_n

         d_sigma_edge_n_d_theta_pol=
     &        +del1*exp1*d_arg_exp1_d_theta
     &        +del1*exp2*d_arg_exp2_d_theta

c         write(*,*)'d_sigma_edge_n_d_theta_pol',
c     &              d_sigma_edge_n_d_theta_pol

         d2_sigma_edge_n_d2_theta_pol=
     &        del1*exp1*(d_arg_exp1_d_theta**2+d2_arg_exp1_d2_theta)
     &       +del2*exp2*(d_arg_exp2_d_theta**2+d2_arg_exp2_d2_theta)
     
c        write(*,*)'d2_sigma_edge_n_d2_theta_pol',
c     &             d2_sigma_edge_n_d2_theta_pol

      endif

      if(i_edge_dens_anal.eq.2) then
c-------using spline ---------------
        if (first) then
c----------create spline coefficients
c
c----------check input data

c           write(*,*)'theta_pol_edge_dens_ar_degree(1)',
c     &     theta_pol_edge_dens_ar_degree(1)
          

           if (dabs(theta_pol_edge_dens_ar_degree(1)).gt.1.d-13) 
     &         then
               write(*,*)'in the input file genray.dat or genray.in'
               write(*,*)'in namelist edge_prof_nml'
               write(*,*)'theta_pol_edge_dens_ar_degree(1).ne.0'
               write(*,*)'theta_pol_edge_dens_ar_degree(1)=',
     &                    theta_pol_edge_dens_ar_degree(1)
               write(*,*)'Please set theta_pol_edge_dens_ar_degree(1)=0'
               write(*,*)'in the input file' 
               stop 'check input edge_prof_nml'
           endif

c           write(*,*)'n_pol_edge_dens',n_pol_edge_dens
c           write(*,*)'theta_pol_edge_dens_ar_degree(n_pol_edge_dens)',
c     &                theta_pol_edge_dens_ar_degree(n_pol_edge_dens)

           if(dabs(theta_pol_edge_dens_ar_degree(n_pol_edge_dens)-
     &              360.d0).gt.1.d-13) then
               write(*,*)'in the input file genray.dat or genray.in'
               write(*,*)'in namelist edge_prof_nml'
               write(*,*)'theta_pol_edge_dens_ar_degree(n_pol_edge_dens)
     &      .ne.360'
            write(*,*)'theta_pol_edge_dens_ar_degree(n_pol_edge_dens)=',
     &      theta_pol_edge_dens_ar_degree(n_pol_edge_dens)
               write(*,*)'Please set 
     &       theta_pol_edge_dens_ar_degree(n_pol_edge_dens)=360'
               write(*,*)'in the input file' 
               stop 'check input edge_prof_nml'
           endif

           do i=1,n_pol_edge_dens
             
             if(theta_pol_edge_dens_ar_degree(i).lt.0.d0) then
               write(*,*)'in the input file genray.dat or genray.in'
               write(*,*)'in namelist edge_prof_nml'
               write(*,*)'theta_pol_edge_dens_ar_degree(i).le.0'
               write(*,*)'at i=',i
               write(*,*)'theta_pol_edge_dens_ar_degree(i)=',
     &                    theta_pol_edge_dens_ar_degree(i)
              write(*,*)'Please change theta_pol_edge_dens_ar_degree(i)'
               write(*,*) 'in input file'
               stop 'check input edge_prof_nml'
             endif

             if(theta_pol_edge_dens_ar_degree(i).gt.360.d0) then
               write(*,*)'in the input file genray.dat or genray.in'
               write(*,*)'in namelist edge_prof_nml'
               write(*,*)'theta_pol_edge_dens_ar_degree(i).gt.360'
               write(*,*)'at i=',i
               write(*,*)'theta_pol_edge_dens_ar_degree(i)=',
     &                    theta_pol_edge_dens_ar_degree(i)
              write(*,*)'Please change theta_pol_edge_dens_ar_degree(i)'
               write(*,*) 'in the input file'
               stop 'check input edge_prof_nml'
             endif

             if (i.gt.1) then
c----------------check that the given poloidal angle array is increasing
                if(theta_pol_edge_dens_ar_degree(i-1).ge.
     &             theta_pol_edge_dens_ar_degree(i)) then
                   write(*,*)'in the input file genray.dat or genray.in'
                   write(*,*)'in namelist edge_prof_nml'
                   write(*,*)'it was given non-momotonic array'
                   write(*,*)'theta_pol_edge_dens_ar_degree(i))'
                   write(*,*)'at i=',i
              write(*,*)'Please change theta_pol_edge_dens_ar_degree(i)'
                   write(*,*) 'in the input file'
                   stop 'check input edge_prof_nml'
               endif
             endif !i.gt.1

             if (sigmedgn_ar(i).le.0.d0) then
                write(*,*)'in the input file genray.dat or genray.in'
                write(*,*)'in namelist edge_prof_nml'
                write(*,*)'sigmedgn_ar(i).le.0 at i=',i
                write(*,*)'sigmedgn_ar(i)=',sigmedgn_ar(i)
                write(*,*)'it shoud be positive'
                write(*,*)'Please change sigmedgn_ar(i)'
                write(*,*) 'in the input file'
                stop 'check input edge_prof_nml'
             endif

           enddo 

           do i=1,n_pol_edge_dens
               theta_pol_edge_dens_ar_radian(i)=
     &         theta_pol_edge_dens_ar_degree(i)*pi/180.d0
           enddo
c-------------------------------------------------------------------
c          creation of spline coefficients for sigmedgn
c-------------------------------------------------------------------     

           iop(1)=3 ! periodic spline boundary conditions
           iop(2)=3  

           call coeff1(n_pol_edge_dens,
     &                 theta_pol_edge_dens_ar_radian,
     &                 sigmedgn_ar,sigmedgn_deriv,
     &                 iop,1,workk)

c----test
c          itabl(1)=1
c          itabl(2)=1
c          itabl(3)=1
c          write(*,*)'theta_pol_edge_dens_ar_radian',
c     &               theta_pol_edge_dens_ar_radian
c          write(*,*)'sigmedgn_ar',sigmedgn_ar
c          write(*,*)'sigmedgn_deriv',sigmedgn_deriv

c          do i=1,n_pol_edge_dens
             
c             call terp1(n_pol_edge_dens,
c     &       theta_pol_edge_dens_ar_radian,
c     &       sigmedgn_ar,sigmedgn_deriv,
c     &       theta_pol_edge_dens_ar_radian(i),1,tabl,itabl)

c            sigma_edge_n=tabl(1)
c            d_sigma_edge_n_d_theta_pol=tabl(2)
c            d2_sigma_edge_n_d2_theta_pol=tabl(3)
c            write(*,*)'i,theta_pol_edge_dens_ar_radian(i)',
c     &                 i,theta_pol_edge_dens_ar_radian(i)
c            write(*,*)'sigma_edge_n',sigma_edge_n
c            write(*,*)'d_sigma_edge_n_d_theta_pol',
c     &                 d_sigma_edge_n_d_theta_pol
c            write(*,*)'d2_sigma_edge_n_d2_theta_pol',
c     &                 d2_sigma_edge_n_d2_theta_pol
c          enddo
c----end test

          first=.false. 
        endif ! if first: create spline coefficients
c
c-------calculate sigmedgen and its derivative by poloidal angle
c  

        itabl(1)=1
        itabl(2)=1
        itabl(3)=1

c        write(*,*)'theta_pol_edge_dens_ar_radian',
c     &             theta_pol_edge_dens_ar_radian
c        write(*,*)'sigmedgn_ar',sigmedgn_ar
c        write(*,*)'sigmedgn_deriv',sigmedgn_deriv
c        write(*,*)'theta_pol_radian',theta_pol_radian
 
        call terp1(n_pol_edge_dens,
     &  theta_pol_edge_dens_ar_radian,
     &  sigmedgn_ar,sigmedgn_deriv,
     &  theta_pol_radian,1,tabl,itabl)

c        write(*,*)'itabl',itabl
c        write(*,*)'tabl',tabl

        sigma_edge_n=tabl(1)
        d_sigma_edge_n_d_theta_pol=tabl(2)
        d2_sigma_edge_n_d2_theta_pol=tabl(3)  

      endif !i_edge_dens_anal.eq.2

 10   continue

      return
      end

      subroutine dens_rho_theta_LCFS(rho_small,theta_pol,i,
     &dens_rho_theta,d_dens_rho_theta_d_rho,
     &d_dens_rho_theta_d_theta)
c-----------------------------------------------------------------------
c     Calculate density and derivatives:d_dentsity/d_rho
c     d_dentsity/d_theta_poloidal
c     versus radius and poloidal angle outside LCFS 
c     1 < rho_small, 0=< theta_pol <2pi
c
c     for i plasma specie
c-----------------------------------------------------------------------
      implicit none

      include 'param.i'
      include 'one.i'
      include 'six.i'
      include 'edge_prof_nml.i'

c-----input
      real*8
     &rho_small, !small radius [non-dimensional]
     &theta_pol  !poloidal angle [radians] 0=< theta <2pi
                 ! =zero at outer side at the equatorial plane
                 ! =pi   at inner side at the equatorial plane
                 ! =pi/2 at the top of the polidal cross-section
!     &dens_min_edge   !minimal density it is set in 'edge_prof.nml'
!                 ! if (dens_rho_theta < dens_min_edge) then
!                 !    dens_rho_theta=dens_min_edge
!                 !    d_dens_rho_theta_d_theta=0

      integer i  ! number of plasma species i=1,..,nbulk
c-----output 
      real*8
     &dens_rho_theta,                !density
     &d_dens_rho_theta_d_rho,        !derivative d_density/d_rho
     &d_dens_rho_theta_d_theta       !derivative d_density/d_theta_poloidal

c-----locals
      real*8 rhoedge,densedge,d2_dens_drho(ndensa),
     &sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     &d2_sigma_edge_n_d2_theta_pol,tabl(3),
     &densedge_e,dens_ratio

      integer itabl(3),k

      integer iwarn
      data iwarn/0/
      save iwarn

      rhoedge=1.d0 !small radius at LCFS

      itabl(1)=1
      itabl(2)=0
      itabl(3)=0 

c      write(*,*)'dense.f dens_rho_theta_LCFS ndens,i,nbulk',
c     & ndens,i,nbulk

      do k=1,ndens
         densm(k)=dens1(k,i)
         d2_dens_drho(k)=d2_dens_drho1(k,i)
      enddo

      call terp1(ndens,rhom,densm,d2_dens_drho,rhoedge,1,tabl,itabl)

      densedge=tabl(1) 

      call sigma_edge_n_theta_pol(theta_pol,
     &sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     &d2_sigma_edge_n_d2_theta_pol)

      dens_rho_theta=densedge*dexp(-(rho_small-1.d0)/sigma_edge_n)

cSAP090222
c      write(*,*)'densedge,i,rho_small,sigma_edge_n,dens_rho_theta',
c     &           i,densedge,rho_small,sigma_edge_n,dens_rho_theta


      d_dens_rho_theta_d_theta = dens_rho_theta*
     &   ((rho_small-1.d0)/sigma_edge_n**2)*d_sigma_edge_n_d_theta_pol

      d_dens_rho_theta_d_rho =- dens_rho_theta/sigma_edge_n
     
     
cSAP090513
c      if (densedge.lt.dens_min_edge) then
c         write(*,*)'WARNINNG:in dens_rho_theta_LCFS'
c         write(*,*)'densedge.lt.dens_min_edge'
c         write(*,*)'It wil be dens_rho_theta=densedge'
c         dens_rho_theta=densedge 
c         d_dens_rho_theta_d_rho=0.d0
c         d_dens_rho_theta_d_theta=0.d0
c      else
c         if (dens_rho_theta.lt.dens_min_edge) then
c         write(*,*)'WARNINNG:in dens_rho_theta_LCFS'
c         write(*,*)'dens_rho_theta.lt.dens_min_edge'
c         write(*,*)'It wil be dens_rho_theta=dens_min_edge '
c            dens_rho_theta=dens_min_edge 
c            d_dens_rho_theta_d_theta=0.d0
c            d_dens_rho_theta_d_rho=0.d0
c         endif
c      endif

c------------------------------------------------------------
cSAP090526
      dens_ratio=1.d0

      if (i.gt.1) then
c----------------------------------------------------------
c        calculates the ratio density(i)/density_e at rho=1
c-----------------------------------------------------------
c        calculate the electron density at rho=1
         do k=1,ndens
            densm(k)=dens1(k,1)
            d2_dens_drho(k)=d2_dens_drho1(k,1)
         enddo
         call terp1(ndens,rhom,densm,d2_dens_drho,rhoedge,1,tabl,itabl)
         densedge_e=tabl(1) 
         dens_ratio=densedge/densedge_e
      endif

cSAP090526     
c      if (dens_rho_theta.lt.dens_min_edge) then
      if (dens_rho_theta.lt.dens_min_edge*dens_ratio) then
         if (iwarn.lt.50) then
         iwarn=iwarn+1   
         write(*,*)'WARNING: [dens_rho_theta_LCFS] rho_small=',rho_small
         !YuP: these messages could be printed too often, 
         !especially in LHCD runs, when rays travel at rho>1 frequently.
         write(*,*)'densedge.lt.dens_min_edge*dens_ratio'
         write(*,*)'It will be dens_rho_theta=dens_min_edge*dens_ratio'
         endif
cSAP090526
c         dens_rho_theta=dens_min_edge
         dens_rho_theta=dens_min_edge*dens_ratio
         d_dens_rho_theta_d_rho=0.d0
         d_dens_rho_theta_d_theta=0.d0
      endif      

      return
      end
!
c
c        *********************dense ***************************
c        *                        -                           
c        * this function calculates the density profile       
c        * as a function  from psi(z,r,phi) 		      
c        *  rho=psi-psimag 				      
c        *  a=psilim-psimag inside the function b()           
c        *  these variables are in common 'one'
c        *                                                    
c        * This function uses 2D spline at RZ mesh outside LCFS
c        * which has density fall near the chamber wall.      
c        *
c        * For n_wall=0 case this function should give density
c        * close to the density created by
*        * FUNCTION dense_no_RZ_spline(z,r,phi,i)
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
      real*8 FUNCTION dense(z,r,phi,i)
      !IMPLICIT double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 z,r,phi !space coordinates
c     rho is in common  /one/
      integer i !plasma species number
c-----local
      real*8  theta_pol, !0=< theta_pol <2pi
     &rho_small,dens_rho_theta,
     &d_dens_rho_theta_d_rho,d_dens_rho_theta_d_theta
c-----externals 
      real*8 density_r_z_i
      real*8 thetapol,   ! -pi <thetapol =<pi
     &vardens,densrho
     
      integer isp,itype !local
      real*8 val,dvalz,dvalr,dvalphi !local

c      write(*,*)'in dense z,r,phi,i,nbulk',z,r,phi,i,nbulk
c      write(*,*)'in dense rho',rho

      if((model_rho_dens.eq.7).and.(dens_read.eq.'enabled')) then 
         ![2024-08-14] added 
         !Interpolate from dengrid_zrp(iz,ir,iphi) points
         !that are adjacent to the given point (z,r,phi)
         itype=1 !=1 for density (and derivs); =2 for temperature
         isp=i !species number
         call interp_zrp(z,r,phi, val,dvalz,dvalr,dvalphi, isp,itype) !val=n()
         dense=val
         return
      endif ![2024-08-14]

      if(i.gt.nbulk)then
c        write(*,*)'in dense i.gt.nbulk: i,nbulk',i,nbulk
	stop
      endif
c      write(*,*)'in dense rho',rho
      if (rho.gt.1.d0-1.d-10) then
c----------------------------------------------------------
c         density outside LCFS
c----------------------------------------------------------
cSAP090403   
          if(n_wall.gt.1) then
c-----------------------------------------------------------
c           calculate density using 2D spline at RZ mesh
c-----------------------------------------------------------
           !write(*,*)'in dense before density_r_z_i'
            dense=density_r_z_i(z,r,phi,i)
           !write(*,*)'in dense after density_r_z_i '
          else
c----------------------------------------------------------
cSAP090206
c           dense=densrho(rho,i)
c----------------------------------------------------------------
c           calculate density using formula versus small radius
c           and poloidal angle
c----------------------------------------------------------------
            theta_pol=thetapol(z,r) ! -pi <thetapol =<pi
            if (theta_pol.lt.0d0) then
              theta_pol=theta_pol+2*pi !pi< theta_pol<2pi
            endif

c           write(*,*)'dense.f in function dense before'
c           write(*,*)'dens_rho_theta_LCFS i',i

            rho_small=rho
   
           !write(*,*)'dense.f z,r,phi,i,rho_small',z,r,phi,i,rho_small

            call dens_rho_theta_LCFS(rho_small,theta_pol,i,
     &        dens_rho_theta,d_dens_rho_theta_d_rho,
     &        d_dens_rho_theta_d_theta)
             !write(*,*)'dense.f in function dense after'
             !write(*,*)'dens_rho_theta_LCFS i',i
            dense = dens_rho_theta 

cSAP090222
c           write(*,*)'in dense i,rho,theta_pol,dense',
c     &                       i,rho,theta_pol,dense
         endif
      else
      !write(*,*)'in dense0 dense',dense
ctest beg
c         write(*,*)'in dense rho,i',rho,i
c         denst=densrho(rho,i)
c         vardent=vardens(z,r,phi)
c	 write(*,*)'in dense denst,vardent',denst,vardent
c tes end
         dense=densrho(rho,i)*(1.0d0+vardens(z,r,phi))
c         write(*,*)'in dense0 dense',dense
      endif
         !write(*,*)'in dense dense',dense
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
      real*8 z,r,phi !space coordinates
      integer i     !number of plasma specie

c-----externals
      real*8  temperho,temperature_r_z_i, rhof

c-----locals
      integer isp,itype !local
      real*8 val,dvalz,dvalr,dvalphi !local
 
      if(i.gt.nbulk)then
        write(*,*)'in tempe i.gt.nbulk: i,nbulk',i,nbulk
	stop
      endif

      if((model_rho_dens.eq.7).and.(temp_read.eq.'enabled')) then 
         ![2024-08-14] added flag/condition for temp_read
         !Interpolate from tempgrid_zrp(iz,ir,iphi) points
         !that are adjacent to the given point (z,r,phi)
         itype=2 !=1 for density (and derivs); =2 for temperature
         isp=i !species number
         call interp_zrp(z,r,phi, val,dvalz,dvalr,dvalphi, isp,itype) !val=T()
         tempe=val
         return
      endif !model_rho_dens.eq.7


      rho=rhof(z,r,phi) !YuP[2024-08-14] added: not always defined ?
      
cBH151016: As discovered by Syun'ichi Shiraiwa, this bug precludes use of
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
            
      idm=nxeqd_add
c      write(*,*)'before terp2p_Sm'       
      temperature_r_z_i=terp2p_Sm(r,z,
     &                    nxeqd_add,rr_add,
     &                    nyeqd_add,zz_add,
     &            temperature_r_z(1,1,i),temperature_r_z_rr(1,1,i),
     &            temperature_r_z_zz(1,1,i),
     &            temperature_r_z_rrzz(1,1,i),idm,0,0,1,nxeqd_add) 
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
            
      idm=nxeqd_add
c      write(*,*)'before terp2p_Sm'       
      zeff_r_z_f=terp2p_Sm(r,z,
     &                    nxeqd_add,rr_add,
     &                    nyeqd_add,zz_add,
     &            zeff_r_z, zeff_r_z_rr,
     &            zeff_r_z_zz,
     &            zeff_r_z_rrzz,idm,0,0,1,nxeqd_add) 
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
c      parameter (dens_nwka=1+3*nxeqd_add_a)   
      
c      real*8      
c     &wk_dens(dens_nwka),dens_norm,dens_spline 

      real*8, dimension(1:3*max0(nxeqd_add,nyeqd_add)):: wk_dens
      real*8 dens_norm,dens_spline 
      integer ibd(4),idm
      integer i,j,k

c-----externals
      real*8 terp2p_Sm

      dens_nwka=1+3*max0(nxeqd_add,nyeqd_add)
      write(*,*)'splcoef_temperature_r_z'
      write(*,*)'nxeqd_add,nyeqd_add,dens_nwka',
     &           nxeqd_add,nyeqd_add,dens_nwka
 
c------------------------------------------------------------------
c     check parameter dens_nwka value 
c------------------------------------------------------------------
      if (dens_nwka.lt.(1+3*nxeqd_add)) then
         write(*,*)'**************************************' 
         WRITE(*,*)'chamber_wall,f in  splcoef_temperatur_r_z'
         WRITE(*,*)'dens_nwka.lt.(1+3*nxeqd_add)'
         WRITE(*,*)'it shoud be'
         WRITE(*,*)'dens_nwka=1+3*max(nxeqd_add,nyeqd_add)'
         WRITE(*,*)'dens_nwka,nxeqd_add,nyeqd_add',
     &              dens_nwka,nxeqd_add,nxeqd_add
c         write(*,*)'Please change parameter in splcoef_temperature_r_z'
c         write(*,*)'for parameter (dens_nwka=1+3*nxeqd_add_a)'
c         write(*,*)'and recompile code'  
         write(*,*)'**************************************' 
         STOP
      endif


      if (dens_nwka.lt.(1+3*nyeqd_add)) then
         write(*,*)'**************************************' 
         WRITE(*,*)'chamber_wall,f in  splcoef_temperature_r_z'
         WRITE(*,*)'dens_nwka.lt.(1+3*nyeqd_add)'
         WRITE(*,*)'it shoud be'
         WRITE(*,*)'dens_nwka=1+3*max(nxeqd_add,nyeqd_add)'
         WRITE(*,*)'dens_nwka,nxeqd_add,nxyqd_add',
     &              dens_nwka,nxeqd_add,nxeqd_add
c         write(*,*)'Please change parameter in splcoef_temperature_r_z'
c         write(*,*)'for parameter (dens_nwka=1+3*nyeqd_add_a)'
c         write(*,*)'and recompile code'  
         write(*,*)'**************************************' 
         STOP
      endif
      
c-----creates 2D spline coefficients for temperature_r_z functions
      ibd(1)=2
      ibd(2)=2
      ibd(3)=2
      ibd(4)=2

      idm=nxeqd_add

      write(*,*)'chamber_wall.f in subroutine splcoef_temperature_r_z'
      write(*,*)'nbulk=',nbulk
      write(*,*)'nxeqd_add,nyeqd_add',nxeqd_add,nyeqd_add
c      write(*,*)'nxeqd_add_a,nyeqd_add_a',nxeqd_add_a,nyeqd_add_a
      

      do k=1,nbulk

c        write(*,*)'k',k
       
         do j=1,nyeqd_add
            temperature_r_z_rr(1,j,k)=0.d0
            temperature_r_z_rr(nxeqd_add,j,k)=0.d0
         enddo    

          do i=1,nxeqd_add
            temperature_r_z_zz(i,1,k)=0.d0
            temperature_r_z_zz(i,nyeqd_add,k)=0.d0
          enddo

c          write(*,*)'before coeff2_Sm(nxeqd_add,rr_'

          call coeff2_Sm(nxeqd_add,rr_add,nyeqd_add,zz_add,
     &      temperature_r_z(1,1,k),
     &      temperature_r_z_rr(1,1,k),temperature_r_z_zz(1,1,k),
     &      temperature_r_z_rrzz(1,1,k),
     &      idm,ibd,wk_dens,nxeqd_add)

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
c     &             temperature_r_z_rrzz(1,1,k),idm,0,0,1,nxeqd_add) 
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
c      parameter (dens_nwka=1+3*nxeqd_add_a)   
      
c      real*8      
c     &wk_dens(dens_nwka),dens_norm,dens_spline 
      real*8, dimension(1:3*max0(nxeqd_add,nyeqd_add)):: wk_dens
      real*8 dens_norm,dens_spline
 
      integer ibd(4),idm
      integer i,j

c-----externals
      real*8 terp2p_Sm

      dens_nwka=1+3*max0(nxeqd_add,nyeqd_add)
      write(*,*)'splcoef_zeff_r_z'
      write(*,*)'nxeqd_add,nyeqd_add,dens_nwka',
     &           nxeqd_add,nyeqd_add,dens_nwka
 
c------------------------------------------------------------------
c     check parameter dens_nwka value 
c------------------------------------------------------------------
      if (dens_nwka.lt.(1+3*nxeqd_add)) then
         write(*,*)'**************************************' 
         WRITE(*,*)'chamber_wall,f in  splcoef_zeff_r_z'
         WRITE(*,*)'dens_nwka.lt.(1+3*nxeqd_add)'
         WRITE(*,*)'it shoud be'
         WRITE(*,*)'dens_nwka=1+3*max(nxeqd_add,nyeqd_add)'
         WRITE(*,*)'dens_nwka,nxeqd_add,nyeqd_add',
     &              dens_nwka,nxeqd_add,nxeqd_add
c         write(*,*)'Please change parameter in splcoef_zefft_r_z'
c         write(*,*)'for parameter (dens_nwka=1+3*nxeqd_add_a)'
c         write(*,*)'and recompile code'  
c         write(*,*)'**************************************' 
         STOP
      endif


      if (dens_nwka.lt.(1+3*nyeqd_add)) then
         write(*,*)'**************************************' 
         WRITE(*,*)'chamber_wall,f in  splcoef_zeff_r_z'
         WRITE(*,*)'dens_nwka.lt.(1+3*nyeqd_add)'
         WRITE(*,*)'it shoud be'
         WRITE(*,*)'dens_nwka=1+3*max(nxeqd_add,nyeqd_add)'
         WRITE(*,*)'dens_nwka,nxeqd_add,nxyqd_add',
     &              dens_nwka,nxeqd_add,nxeqd_add
c         write(*,*)'Please change parameter in splcoef_zeff_r_z'
c         write(*,*)'for parameter (dens_nwka=1+3*nyeqd_add_a)'
c         write(*,*)'and recompile code'  
         write(*,*)'**************************************' 
         STOP
      endif
      
c-----creates 2D spline coefficients for zeff_r_z functions
      ibd(1)=2
      ibd(2)=2
      ibd(3)=2
      ibd(4)=2

      idm=nxeqd_add

      write(*,*)'chamber_wall.f in subroutine splcoef_zeff_r_z'
      write(*,*)'nbulk=',nbulk
      write(*,*)'nxeqd_add,nyeqd_add',nxeqd_add,nyeqd_add
c      write(*,*)'nxeqd_add_a,nyeqd_add_a',nxeqd_add_a,nyeqd_add_a
      

    
       
      do j=1,nyeqd_add
         zeff_r_z_rr(1,j)=0.d0
         zeff_r_z_rr(nxeqd_add,j)=0.d0
      enddo    

      do i=1,nxeqd_add
         zeff_r_z_zz(i,1)=0.d0
         zeff_r_z_zz(i,nyeqd_add)=0.d0
      enddo

c      write(*,*)'in splcoef_zeff before coeff2_Sm(nxeqd_add,rr_'

c      do j=1,nyeqd_add 
c         do i=1,nxeqd_add
c            write(*,*)'j,i,zeff_r_z(i,j)',j,i,zeff_r_z(i,j)
c         enddo
c      enddo    

c      stop 'in dense.f splcoef_zeff'



      call coeff2_Sm(nxeqd_add,rr_add,nyeqd_add,zz_add,
     &      zeff_r_z,
     &      zeff_r_z_rr,zeff_r_z_zz,
     &      zeff_r_z_rrzz,
     &      idm,ibd,wk_dens,nxeqd_add)

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
c     &             zeff_r_z_rrzz,idm,0,0,1,nxeqd_add) 
c              write(*,*)'after terp2p zeff_spline',zeff_spline
         
c              zeff_norm=zeff_norm+dabs(zeff_spline-zeff_r_z(i,j))

c            enddo
c         enddo 
c         write(*,*)'zeff_norm',zeff_norm     
c-endtest

      return
      end



      subroutine create_edge_prof_table
c-----calculate tables for edge_prof
c     theta_pol_edge_dens_ar_degree(i=1,n_pol_edge_dens)
c     sigmedgn_ar(i=1,n_pol_edge_dens)
c
c     from the input data for analytical profile
c     using analytical profile  like in
c     subroutine sigma_edge_n_theta_pol(theta_pol_radian
c
c     calculated profiles will be in edge_prof_nml.i

      implicit none

      include 'param.i'
      include 'edge_prof.i'

c-----locals
      integer i
      real*8 pi,step_degree,exp1,exp2,del1,del2

      pi=4.d0*datan(1.d0)
     
      sigma_theta_pol_edge_1_radian=
     &           sigma_theta_pol_edge_1_degree*pi/180.d0 

      sigma_theta_pol_edge_2_radian=
     &           sigma_theta_pol_edge_2_degree*pi/180.d0

      theta_pol_edge_1_radian=theta_pol_edge_1_degree*pi/180.d0

      theta_pol_edge_2_radian=theta_pol_edge_2_degree*pi/180.d0

      step_degree=360.d0/(n_pol_edge_dens-1)

      do i=1,n_pol_edge_dens
       
         theta_pol_edge_dens_ar_degree(i)=step_degree*(i-1)

         write(*,*)'i,theta_pol_edge_dens_ar_degree(i)',
     &              i,theta_pol_edge_dens_ar_degree(i)

         theta_pol_edge_dens_ar_radian(i)=
     &         theta_pol_edge_dens_ar_degree(i)*pi/180.d0

          exp1=dexp(-((theta_pol_edge_dens_ar_radian(i)-
     &               theta_pol_edge_1_radian)/
     &               sigma_theta_pol_edge_1_radian)**2)

          exp2=dexp(-((theta_pol_edge_dens_ar_radian(i)-
     &               theta_pol_edge_2_radian)/
     &               sigma_theta_pol_edge_2_radian)**2)

          del1=sigma_edgen_1-sigma_edgen_0
          del2=sigma_edgen_2-sigma_edgen_0

          sigmedgn_ar(i)=sigma_edgen_0+del1*exp1+del2*exp2

      enddo


      return
      end

c======================================================================
c======================================================================


c        *********************tpop_zrp ***********************
c        * this function calculates Tpop==T_par/T_perp  
c        * profile as a function  of (z,r,phi)       	 
c        **************************************************
c
c------------------------------------------------------------------
c								   
c        input parameters					   
c      								   
c      z,r,phi - coordinates of the  point    
c              where  Tpop  is calculated.      		         	    
c------------------------------------------------------------------
	real*8 function tpop_zrp(z,r,phi,i)
      implicit none !integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i' ! stores rho
      include 'three.i'
      include 'fourb.i'
      include 'five.i'
      !--------------
      real*8 z,r,phi !IN
      integer i !IN (species)
      integer itype,isp !local
      real*8 val,dvalx,dvaly,dvalz, den !local
      real*8 dvalr,dvalphi
      real*8 tpoprho !external func. tpoprho(rho,i)
      real*8 rhof !external func.
      real*8 rholoc !local rho

      if((model_rho_dens.eq.7).and.(tpop_read.eq.'enabled')) then 
         ![2024-08-14] added 
         !Interpolate from tpoprid_zrp(iz,ir,iphi) points
         !that are adjacent to the given point (z,r,phi)
         itype=3 !=1 for density (and derivs); =2 for T, =3 for Tpop
         isp=i !species number
         call interp_zrp(z,r,phi, val,dvalz,dvalr,dvalphi, isp,itype) !val=Tpop()
         tpop_zrp=val
         !Note: no rho is needed 
         return
      endif !model_rho_dens.eq.7
	
      !All other cases - based on rho 
      rholoc=rhof(z,r,phi)
      if(rholoc.gt.1.d0)then
      tpop_zrp=1.d0 !assume Tpar=Tper outside of LCFS
      else
	tpop_zrp=tpoprho(rholoc,i) !This is based on 1D spline over rho-grid
	endif
      !Includes case of (model_rho_dens.eq.7) but (tpop_read.eq.'disabled')

      return
      end function tpop_zrp
