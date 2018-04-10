c 02/11/14
c     The dielectric tensor, its derivatives for
c     fast wave conditions using Masayuki Ono approximation:
c     Masayuki Ono, High harmonics fast waves in high beta plasmas,
c     Phys.Plasmas 2 (11),November 1995, p.4075-4081
c
c     The signs of the gyro frequencies for electrons and ions are positive
c     So, Y_e and Y_i has the same sign.

      subroutine ono_tens(mass_ar,X_ar,Y_ar,T_ar,nbulk,
     & nll_in,np_in,iherm,i_deriv,
     & K,d,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)
c-----calculates dielectric tensor K ,the dispersion function d
c     and the derivatives from the dielectric tensor (for i_deriv=1)
c     dK_dx_ar(3,3,nbulk)  d_eps(i,j)/d_X(s)
c     dK_dy_ar(3,3,nbulk)  d_eps(i,j)/d_Y(s)
c     dK_dt_ar(3,3,nbulk)  d_eps(i,j)/d_T(s)  here T in [eV]
c     dK_dnpar(3,3)        d_eps(i,j)/d_N_perpendicular
c     dK_dnper(3,3)        d_eps(i,j)/d_N_parallel
c     for fast waves in Ono approximation      
c     Masayuki Ono, High harmonics fast waves in high beta plasmas
c     Phys.Plasmas 2 (11),November 1995, p.4075-4081
c  
c     iherm=1 gives all results for hermitian tensor
c          =2  for complete tensor
c------------------------------------------------------------
c     INPUTS:
c      nbulk is a number of the plasma cpecies
c      mass_ar(nbulk) - the mass of the given specie (in electron mass) 
c      X_ar(nbulk) = (fp/f)**2
c      Y_ar(nbulk) = fc/f ! for electron Ye should has the sign opposite to Yi
c      T_ar(nbulk)=temperature in eV  
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c      iherm =1 Hermitian dielectric tensor, 2-full tensor
c      i_deriv=0 no derivatives calculations
c             =1 calculate derivatives dd_dx_ar=dD/dX_ar, dd_dy_ar=Dd/d_Y_ar
c                dd_dt_ar=dD/dT_ar,dd_dnpar=dD/dN_parallel,dd_dnper=dD/dN_perpendicular     
c     OUTPUT:
c      K(3,3):  the nine complex components of the dielectric 
c               evaluated at (mass,X,Y,T,Nll,np)
c      d        dispersion function
c      dK_dx_ar(3,3,nbulk)  d_eps(i,j)/d_X(s)
c      dK_dy_ar(3,3,nbulk)  d_eps(i,j)/d_Y(s)
c      dK_dt_ar(3,3,nbulk)  d_eps(i,j)/d_T(s)  here T in [eV]
c      dK_dnpar(3,3)        d_eps(i,j)/d_N_perp
c      dK_dnper(3,3)        d_eps(i,j)/d_N_parallel
c----------------------------------------------
       
      implicit none
cSm030226
      include 'param.i'
c-----input
      integer nbulk
      double precision mass_ar(*),X_ar(*),Y_ar(*),T_ar(*),
     & nll_in,np_in
      integer iherm,i_deriv
c-----output
      double complex K(3,3),d,
     &dK_dx_ar(3,3,*),dK_dy_ar(3,3,*),dK_dt_ar(3,3,*),
     &dK_dnpar(3,3),dK_dnper(3,3)
c-----locals
      double complex ci
      double precision npar,nper,pi,c,npars,npers
cSm030519
c      double precision y_0,k4,vt,vt_dc,delta_m,delta_x,delta
      double precision y_0,k4,vt,vt_dc,delta_x,delta
      double complex delta_m
      integer*4 IERR
      double complex CZ0,CZ1,CZ2,cy_0,
     &K_xxc,K_xyc,K_yyc,K_xzc,K_yz,K_zz,K_zze,
     &a_coef,b_coef,c_coef,
     &Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz
      integer ncomp_l
cSm030226
c      parameter (ncomp_l=4)
      parameter (ncomp_l=nbulka)

      double complex
     &dK_dx_ar_p(3,3,ncomp_l),dK_dy_ar_p(3,3,ncomp_l),
     &dK_dt_ar_p(3,3,ncomp_l),
     &dK_dnpar_p(3,3),dK_dnper_p(3,3)
      integer i,j,l 

c_test
      double precision det,r_a,r_b,r_c,r_d,r_delta_m,r_K_zze
      double precision npers_p,npers_m 
c_endtest

cSm030512
c test ono using cold
      double precision xe,ye,epsper,g,epspar,xi,yi
c end test ono using cold
     
      if(ncomp_l.lt.nbulk) then
        write(*,*)'ono_disp.f ono_tens ncomp_l=',ncomp_l,'nbulk=',nbulk
        write(*,*)'it should be ncompl_l.ge.nbulk'
        write(*,*)'ncomp_l=nbulka'
        write(*,*)'change parameter nbulka in param.i and recompile'
        stop
      endif
       
      if ((iherm.lt.1).or.(iherm.gt.2)) then
         write(*,*)'ono_disp.f ono_tens iherm=',iherm
         write(*,*)'it should be iherm=1 or iherm=2'
         stop
      endif

      pi = 4*datan(1.d0)
      c = 3.0d10                    !speed of light
      k4 = 4.19d7                   !constant for elec therm vel
    
      ci = ( 0.0d0,1.0d0)
        
c-----if n<nmin then set n=nmin
      call npnllmin(nll_in,np_in,npar,nper)
cSm030512
      npar=nll_in
      nper=np_in

      npars = npar**2
      npers = nper**2

      delta=0.d0
      do i=1,nbulk
         vt =k4*dsqrt(2.0d0*t_ar(i)/mass_ar(i))   !vt= sqrt(2.0*kT[eV]/mass) cm/sec 
         vt_dc = vt / c                
         if (i.eq.1) then  !electron terms
            K_xxc=1.d0+X_ar(1)/Y_ar(1)**2 ! (3a)
            K_xyc=X_ar(1)/Y_ar(1)         ! (3b)
            y_0= 1.d0/(npar*vt_dc) 
            cy_0=dcmplx(y_0,0.d0)
            call CZETA0(npar,cy_0,CZ0,CZ1,IERR)
            delta_m=X_ar(1)/(Y_ar(1)**2*npar)*vt_dc*CZ0         !(3c)

c        write(*,*)'ono_tens X_ar(1),Y_ar(1),npar,CZ0,vt_dc,delta_m',
c     &      X_ar(1),Y_ar(1),npar,vt_dc,CZ0,delta_m

            delta_x=0.5d0*vt_dc**2*npar/y_ar(1)                 !(3e)
            K_zz=-x_ar(1)*CZ1/(npar*vt_dc)**2                   !(3f)
            K_zze=K_zz
            K_yz=-nper*delta_x*K_zz                             !(3e)
         endif ! electron terms

         if (i.gt.1) then !ions terms
            K_xxc=K_xxc-X_ar(i)/(1.d0-Y_ar(i)**2)               !(3a) 
            K_xyc=K_xyc+X_ar(i)*Y_ar(i)/(1.d0-Y_ar(i)**2)       !(3b)
            delta=delta-x_ar(i)*vt_dc**2/(1.d0-Y_ar(i)**2)**2   !(3d)
         endif   

      enddo
      
      K_yyc = K_xxc+npers*delta_m                                 !(3c)
      K_xzc = nper*npar*delta                                     !(3d)

      K(1,1)=K_xxc
      K(1,2)=-ci*K_xyc
      K(1,3)=K_xzc
      K(2,1)=-K(1,2)
      K(2,2)=K_yyc
      K(2,3)=ci*K_yz
      K(3,1)=K(1,3)
      K(3,2)=-K(2,3)
      K(3,3)=K_zz

c     Hermitian part
      if (iherm.eq.1) then
          Kxx = 0.5D0 * ( K(1,1) + dconjg(K(1,1) )) 
          Kyy = 0.5D0 * ( K(2,2) + dconjg(K(2,2) )) 
          Kzz = 0.5D0 * ( K(3,3) + dconjg(K(3,3) )) 
          Kxy = 0.5D0 * ( K(1,2) + dconjg(K(2,1) ))
          Kxz = 0.5D0 * ( K(1,3) + dconjg(K(3,1) )) 
          Kyz = 0.5D0 * ( K(2,3) + dconjg(K(3,2) )) 
          Kyx = 0.5D0 * ( K(2,1) + dconjg(K(1,2) )) 
          Kzx = 0.5D0 * ( K(3,1) + dconjg(K(1,3) )) 
          Kzy = 0.5D0 * ( K(3,2) + dconjg(K(2,3) ))
      endif

c     full tensor
      if (iherm.eq.2) then
          Kxx = K(1,1)
          Kyy = K(2,2) 
          Kzz = K(3,3) 
          Kxy = K(1,2) 
          Kxz = K(1,3) 
          Kyz = K(2,3) 
          Kyx = K(2,1) 
          Kzx = K(3,1) 
          Kzy = K(3,2)
      endif

      K(1,1) =  Kxx
      K(2,2) =  Kyy
      K(3,3) =  Kzz
      K(1,2) =  Kxy
      K(1,3) =  Kxz
      K(2,1) =  Kyx
      K(2,3) =  Kyz
      K(3,1) =  Kzx
      K(3,2) =  Kzy

c-----dispersion function d=a_coef*npers**2+b_coef*npers+c_coef
ctest delta_m=0
c      delta_m=0.d0
c      delta_x=0.d0
cendtest 
      a_coef=((npars-K_xxc)-npars*(1.d0+delta)**2)*(1.d0-delta_m)
      b_coef=-K_xyc**2-delta_x**2*(npars-K_xxc)*K_zze**2-
     &      (npars-K_xxc)*(1.d0-delta_m)*K_zze+
     &      2.d0*delta_x*npar*(1.d0+delta)*K_zze*K_xyc+(npars-K_xxc)**2  
     &      -npars*(1.d0+delta)**2*(npars-K_xxc)
      c_coef=(K_xyc**2-(npars-K_xxc)**2)*K_zze

c     dispersion function from the Ono article
      d=a_coef*npers**2+b_coef*npers+c_coef                !(5)
cSm030511
c      goto 10

c_test
      r_delta_m=dreal(delta_m)
     

      r_K_zze=dreal(K_zze)
      r_a=((npars-K_xxc)-npars*(1.d0+delta)**2)*(1.d0-r_delta_m)
      r_b=-K_xyc**2-delta_x**2*(npars-K_xxc)*r_K_zze**2-
     &     (npars-K_xxc)*(1.d0-r_delta_m)*r_K_zze+
     &     2.d0*delta_x*npar*(1.d0+delta)*r_K_zze*K_xyc+(npars-K_xxc)**2  
     &     -npars*(1.d0+delta)**2*(npars-K_xxc)
      r_c=(K_xyc**2-(npars-K_xxc)**2)*r_K_zze
      det=r_b**2-4.d0*r_a*r_c
      if (det.ge.0.d0) then
         npers_p=(-r_b+dsqrt(det))/(2.d0*r_a)
         npers_m=(-r_b-dsqrt(det))/(2.d0*r_a)
      endif
c      write(*,*)'ono_temns npers_p,npers_m',npers_p,npers_m
cendtest

c     dispersion function calculated as determinant
cSm030512
ctest ono using cold tensor
c         xe=x_ar(1)
c         ye=y_ar(1)
c	 epsper=1.d0-xe/(1.d0-ye*ye)
c         g=xe*ye/(1.d0-ye*ye)
c         epspar=1-xe

c	 do i=2,nbulk
c            xi=x_ar(i)
c            yi=y_ar(i)
c	    epsper=epsper-xi/(1.d0-yi*yi)
c	    epspar=epspar-xi
c            g=g-xi*yi/(1.d0-yi*yi)
c	 enddo
c	 K(1,1)=dcmplx(epsper,0.d0)
c	 K(1,2)=dcmplx(0.d0,g)
c	 K(2,1)=-K(1,2)
c	 K(2,2)=K(1,1)
c	 K(3,3)=dcmplx(epspar,0.d0)
c	 K(1,3)=dcmplx(0.d0,0.d0)
c	 K(3,1)=dcmplx(0.d0,0.d0)
c	 K(2,3)=dcmplx(0.d0,0.d0)
c	 K(3,2)=dcmplx(0.d0,0.d0)
cendtest ono using cold
      
      d =(Kxx-npars) * (Kyy-npars-npers) * (Kzz-npers)
     .+  Kxy * Kyz * (Kzx+nper*npar) 
     .+ (Kxz+nper*npar) * Kyx * Kzy
     .- (Kzx+nper*npar) * (Kyy-npars-npers) * (Kxz+nper*npar)
     .-  Kzy * Kyz * (Kxx-npars)
     .- (Kzz - npers) * Kyx * Kxy

c       d=npers-(-r_b+dsqrt(det))/(2.d0*r_a)

c     write(*,*)'in ono_disp ono_tens npar,nper',npar,nper
c      write(*,*)'Kxx,Kxy,Kxz',Kxx,Kxy,Kxz 
c      write(*,*)'Kyx,Kyy,Kyz',Kyx,Kyy,Kyz
c     write(*,*)'Kzx,Kzy,Kzz',Kzx,Kzy,Kzz 
c     write(*,*)'ono_disp.f in ono_tens d determinant=',d
 10   continue
      if (i_deriv.eq.1) then
c--------calculate derivatives
         do i=1,nbulk
            vt =k4*dsqrt(2.0d0*t_ar(i)/mass_ar(i))   !vt= sqrt(2.0* kT[eV]/mass) 
            vt_dc = vt / c                

            if (i.eq.1) then  !electron terms
          
               dK_dx_ar(1,1,1)=1.d0/Y_ar(1)**2            
               dK_dy_ar(1,1,1)=-2.d0*X_ar(1)/Y_ar(1)**3           
               dK_dt_ar(1,1,1)=0.d0
               dK_dnper(1,1)=0.d0
               dK_dnpar(1,1)=0.d0

               dK_dx_ar(1,2,1)=-ci/Y_ar(1)            
               dK_dy_ar(1,2,1)= ci*X_ar(1)/Y_ar(1)**2
               dK_dt_ar(1,2,1)=0.d0
               dK_dnper(1,2)=0.d0
               dK_dnpar(1,2)=0.d0

               dK_dx_ar(1,3,1)=0.d0
               dK_dy_ar(1,3,1)= 0.d0
               dK_dt_ar(1,3,1)=0.d0
               dK_dnper(1,3)=npar*delta
               dK_dnpar(1,3)=nper*delta
             
c-------------------------------------------------------------------
c              derivatives from K(2,2) by Xe,Ye,Te nper and npar 
               dK_dx_ar(2,2,1)=dK_dx_ar(1,1,1)+npers*
     &         delta_m/X_ar(1)
c      &        1.d0/(Y_ar(1)**2*npar)*vt_dc*CZ0

               dK_dy_ar(2,2,1)=dK_dy_ar(1,1,1)+npers*
     &         (-2.d0*X_ar(1)/(Y_ar(1)**3*npar)*vt_dc*CZ0)

               dK_dt_ar(2,2,1)=dK_dt_ar(1,1,1)+npers*X_ar(1)/Y_ar(1)**2
     &         *((vt_dc/(2.d0*t_ar(1)*npar))*CZ0-
     &         CZ1/(npars*2.d0*t_ar(1)))                    !1/[EV}

               dK_dnper(2,2)=dK_dnper(1,1)+2.d0*nper*delta_m

               dK_dnpar(2,2)=dK_dnpar(1,1)+npers*X_ar(1)/Y_ar(1)**2*
     &         (-vt_dc*CZ0/npars-CZ1/(npars*npar))

c-------------------------------------------------------------------
c              derivatives from K(3,3) by Xe,Ye,Te nper and npar 
c------------------------------------------------------------------
               CZ2=-2.d0*(CZ0+y_0*CZ1)              !d^2 Z_0(y_0)/d^2y_0
               dK_dx_ar(3,3,1)=-CZ1/(vt_dc**2*npars)
               dK_dy_ar(3,3,1)= 0.d0
               dK_dt_ar(3,3,1)=X_ar(1)/npars*
     &         (CZ1/(vt_dc**2*t_ar(1))+CZ2/(2.d0*t_ar(1)*vt_dc**3*npar))

c      write(*,*)'ono_tens X_ar(1),npars,CZ1,CZ2,vt_dc,dK_dt_ar(3,3,1)',
c     &         X_ar(1),npars,CZ1,CZ2,vt_dc,dK_dt_ar(3,3,1)

               dK_dnper(3,3)=0.d0  

               dK_dnpar(3,3)=2.d0*X_ar(1)*CZ1/(vt_dc**2*npars*npar)+
     &                       X_ar(1)*CZ2/(vt_dc**3*npars*npars)
c-------------------------------------------------------------------
c              derivatives from K(2,3) by Xe,Ye,Te nper and npar    
c--------------------------------------------------------------------
               dK_dx_ar(2,3,1)=-ci*nper*delta_x*K_zz/X_ar(1)
               dK_dy_ar(2,3,1)= ci*nper*delta_x/Y_ar(1)*K_zz 
               dK_dt_ar(2,3,1)=-ci*nper*(delta_x/t_ar(1)*K_zz+
     &                         delta_x* dK_dt_ar(3,3,1))

               dK_dnper(2,3)=-ci*(delta_x*K_zz+
     &                        nper*delta_x*dK_dnper(3,3))  

               dK_dnpar(2,3)=-ci*nper*(delta_x/npar*K_zz+
     &                       delta_x*dK_dnpar(3,3))
c-------------------------------------------------------------------
c              derivatives  from K(2,1)
c-----------------------------------------------------------------
               dK_dx_ar(2,1,1)= -dK_dx_ar(1,2,1)
               dK_dy_ar(2,1,1)= -dK_dy_ar(1,2,1)
               dK_dt_ar(2,1,1)= -dK_dt_ar(1,2,1)
               dK_dnper(2,1)  = -dK_dnper(1,2)
               dK_dnpar(2,1)  = -dK_dnpar(1,2)
c-----------------------------------------------------------------
c              derivatives from K(3,1)
c-----------------------------------------------------------------
               dK_dx_ar(3,1,1)=  dK_dx_ar(1,3,1)
               dK_dy_ar(3,1,1)=  dK_dy_ar(1,3,1)
               dK_dt_ar(3,1,1)=  dK_dt_ar(1,3,1)
               dK_dnper(3,1)  =  dK_dnper(1,3)
               dK_dnpar(3,1)  =  dK_dnpar(1,3)
c-----------------------------------------------------------------
c              derivatives from K(3,2)
c-----------------------------------------------------------------
               dK_dx_ar(3,2,1)= -dK_dx_ar(2,3,1)
               dK_dy_ar(3,2,1)= -dK_dy_ar(2,3,1)
               dK_dt_ar(3,2,1)= -dK_dt_ar(2,3,1)
               dK_dnper(3,2)  = -dK_dnper(2,3)
               dK_dnpar(3,2)  = -dK_dnpar(2,3)
c-----------------------------------------------------------------

            endif ! electron terms

            if (i.gt.1) then !ions terms

               dK_dx_ar(1,1,i)=-1.d0/(1.d0-Y_ar(i)**2)
               dK_dy_ar(1,1,i)=-2.d0*X_ar(i)*Y_ar(i)/
     &                          (1.d0-Y_ar(i)**2)**2
               dK_dt_ar(1,1,i)=0.d0

               dK_dx_ar(1,2,i)=-ci*Y_ar(i)/(1.d0-Y_ar(i)**2)
               dK_dy_ar(1,2,i)=-ci*X_ar(i)*(1.d0+Y_ar(i)**2)/
     &                         (1.d0-Y_ar(i)**2)**2
               dK_dt_ar(1,2,i)=0.d0

               dK_dx_ar(1,3,i)=-nper*npar*vt_dc**2/
     &                         (1.d0-Y_ar(i)**2)**2
               dK_dy_ar(1,3,i)=-nper*npar*4.d0*X_ar(i)*Y_ar(i)*vt_dc**2
     &                         /(1.d0-Y_ar(i)**2)**3
               dK_dt_ar(1,3,i)=-nper*npar*X_ar(i)*vt_dc**2/
     &                          ((1.d0-Y_ar(i)**2)**2*t_ar(i))
c-------------------------------------------------------------------
c              derivatives from K(2,2) by Xi,Yi,Ti  
               dK_dx_ar(2,2,i)=dK_dx_ar(1,1,i)
               dK_dy_ar(2,2,i)=dK_dy_ar(1,1,i)
               dK_dt_ar(2,2,i)=dK_dt_ar(1,1,i)           !1/[EV}

c-------------------------------------------------------------------
c              derivatives from K(3,3) by Xi,Yi,Ti  
c------------------------------------------------------------------
               dK_dx_ar(3,3,i)=0.d0
               dK_dy_ar(3,3,i)=0.d0
               dK_dt_ar(3,3,i)=0.d0

c-------------------------------------------------------------------
c              derivatives from K(2,3) by Xi,Yi,Ti    
c--------------------------------------------------------------------
               dK_dx_ar(2,3,i)=0.d0
               dK_dy_ar(2,3,i)=0.d0 
               dK_dt_ar(2,3,i)=0.d0

c-------------------------------------------------------------------
c              derivatives  from K(2,1)
c-----------------------------------------------------------------
               dK_dx_ar(2,1,i)= -dK_dx_ar(1,2,i)
               dK_dy_ar(2,1,i)= -dK_dy_ar(1,2,i)
               dK_dt_ar(2,1,i)= -dK_dt_ar(1,2,i)
c-----------------------------------------------------------------
c              derivatives from K(3,1)
c-----------------------------------------------------------------
               dK_dx_ar(3,1,i)=  dK_dx_ar(1,3,i)
               dK_dy_ar(3,1,i)=  dK_dy_ar(1,3,i)
               dK_dt_ar(3,1,i)=  dK_dt_ar(1,3,i)
c-----------------------------------------------------------------
c              derivatives from K(3,2)
c-----------------------------------------------------------------
               dK_dx_ar(3,2,i)= -dK_dx_ar(2,3,i)
               dK_dy_ar(3,2,i)= -dK_dy_ar(2,3,i)
               dK_dt_ar(3,2,i)= -dK_dt_ar(2,3,i)
c-----------------------------------------------------------------            
            endif   
         enddo
     
                
c--------Hermitian part of the derivatives
         if (iherm.eq.1) then
             do j=1,3
                do l=1,3
                   dK_dnper_p(j,l)=dK_dnper(j,l)
                   dK_dnpar_p(j,l)=dK_dnpar(j,l)
                   do i=1,nbulk
                      dK_dx_ar_p(j,l,i)=dK_dx_ar(j,l,i)
                      dK_dy_ar_p(j,l,i)=dK_dy_ar(j,l,i)
                      dK_dt_ar_p(j,l,i)=dK_dt_ar(j,l,i)
                   enddo
                enddo
             enddo

              do j=1,3
                do l=1,3
                   dK_dnper(j,l)=0.5d0*
     &                        (dK_dnper_p(j,l)+dconjg(dK_dnper_p(l,j)))
                   dK_dnpar(j,l)=0.5d0*
     &                        (dK_dnpar_p(j,l)+dconjg(dK_dnpar_p(l,j)))
                   do i=1,nbulk
                      dK_dx_ar(j,l,i)=0.5d0*
     &                (dK_dx_ar_p(j,l,i)+dconjg(dK_dx_ar_p(l,j,i)))
                      dK_dy_ar(j,l,i)=0.5d0*
     &                (dK_dy_ar_p(j,l,i)+dconjg(dK_dy_ar_p(l,j,i)))
                      dK_dt_ar(j,l,i)=0.5d0*
     &                (dK_dt_ar_p(j,l,i)+dconjg(dK_dt_ar_p(l,j,i)))
                   enddo
                enddo
             enddo
         endif !Hermition part of derivatives
     
      endif !ideriv.eq.1 derivatives calculations

c      do i=1,nbulk
c      do j=1,3
c         do l=1,3
c           write(*,*)'!!! end of ono_tens j,l,dK_dx_ar(j,l,i)',
c     &      j,l,dK_dx_ar(j,l,i)
c         enddo
c      enddo
c      enddo


      return
      end


      subroutine DdOno(nbulk,mass_ar,x_ar,y_ar,t_ar,
     .nll_in,np_in,dd_ono,dd_ono_np_h,dd_ono_nll_h,
     .i_dd_ono_np,dd_ono_np)
c
c     This function returns the derivatives of the Hermitian
c     part of the Ono dispersion function (used for fast waves)
c     with respect to X,Y,Te_av,nll_in, and np_in
c
c     INPUTS:
c      nbulk the total number of plasma species
c      mass_ar - the masses  of the plasma species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=average temperature in ev    
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c      i_dd_ono_np =0 no calcultions of dd_ono_np (dD/dN_perp for full tensor)
c                  =1 calculates        dd_ono_np (dD/dN_perp for full tensor)   
c     OUTPUT:
c      DD_ono(1-3,nbulk) derivatives from the dispersion function (for hot
c                  Hermitian dielectric tensor ) with respect
c                  X_js,Y_js,T_av_js (js=1,nbulk)
c      DD_ono(1,js)=dD/dX_js
c      DD_ono(2,js)=dD/dY_js
c      DD_ono(3,js)=dD/dT_av_js
c
c      derivatives from the dispersion function (for hot
c      Hermitian dielectric tensor ) with respect
c      dd_ono_np_h  = dD/DN_perpendicular
c      dd_ono_nll_h = dD/dN_parallel 
c
c      derivative from the dispersion function (for hot
c      complete dielectric tensor ) with respect
c      dd_ono_np = dD/DN_perpendicular
c-----------------------------------------------------------
      Implicit none
cSm030226
      include 'param.i'
c-----input
      integer nbulk
      double precision mass_ar(*)
      double precision X_ar(*),Y_ar(*),T_ar(*)
      double precision np_in,nll_in
      integer i_dd_ono_np
c-----external
c-----output
c     derivatives from the Hermitian dispersion function D
c     dd_ono(1,js)=dD/dX_js,dd_ono(2,js)=dD/dY_js,dd_ono(3,js)=dD/dT_av_js
c     dd_ono_np_h=dD/N_perp,dd_ono_nll_h=dD/dN_parallel
c     derivartive from the complete dispersion function 
c     dd_ono_np=dD/N_erp
      double complex dd_ono(3,*) !(3,nbulka) 
      double complex dd_ono_nll_h,dd_ono_np_h
      double complex dd_ono_np
 
c-----locals
      double precision nperp,nll,nlls,nps
      double complex k_herm(3,3),k_sum(3,3)

      integer ncomp_l
cSm030226
c      parameter (ncomp_l=4)
      parameter (ncomp_l=nbulka)

      double complex d_ono,d_ono_h,dd_ono_nll,
     &dK_dx_ar(3,3,ncomp_l),dK_dy_ar(3,3,ncomp_l),dK_dt_ar(3,3,ncomp_l),
     &dK_dnpar(3,3),dK_dnper(3,3),
     &dK_dx_ar_h(3,3,ncomp_l),dK_dy_ar_h(3,3,ncomp_l),
     &dK_dt_ar_h(3,3,ncomp_l),dK_dnpar_h(3,3),dK_dnper_h(3,3)

      double complex ddeps(3,3),ddeps_h(3,3)
      double complex dk(3,3,6)
      integer i_deriv,iherm
      integer js,j,i1,i2

      if(ncomp_l.lt.nbulk) then
        write(*,*)'ono_disp.f Dd_Ono ncomp_l=',ncomp_l,'nbulk=',nbulk
        write(*,*)'it should be ncompl_l.ge.nbulk'
        write(*,*)'ncomp_l=nbulka'
        write(*,*)'change parameter nbulka in param.i recompile'
        stop
      endif
       
c-----put the minimal values for nll,nperp
      call npnllmin(nll_in,np_in,nll,nperp)

      nlls=nll*nll
      nps=nperp*nperp

      i_deriv=1 ! Hermitian tensor and its derivatives
      iherm=1  
   
      call ono_tens(mass_ar,X_ar,Y_ar,T_ar,nbulk,
     & nll_in,np_in,iherm,i_deriv,
     & k_herm,d_ono_h,dK_dx_ar_h,dK_dy_ar_h,dK_dt_ar_h,
     & dK_dnper_h,dK_dnpar_h)

c       write(*,*)'DdOno after call ono_tens iherm=1 '
c       do js=1,nbulk
c          do i1=1,3
c            do i2=1,3
c               write(*,*)'js,i1,i2,dK_dx_ar_h(i1,i2,js)',
c     &         js,i1,i2,dK_dx_ar_h(i1,i2,js)
c            enddo
c          enddo
c       enddo

c      compute the derivatives from the disperstion function D
c     for the hermitian hot tensor k_herm
c     with respect dielectric tensor components: ddeps_h
c     with respect N_perpendicular ,N_parallel : ddnp_h,ddnll_h)

      call dddeps(k_herm,nll,nperp,ddeps_h,dd_ono_nll_h,dd_ono_np_h)

      if(i_dd_ono_np.eq.0) then
        dd_ono_np=dcmplx(0.d0,0.d0)
      endif
 
      if(i_dd_ono_np.eq.1) then
        i_deriv=1 ! complete tensor and its derivatives
        iherm=2   
                                                                      
        call ono_tens(mass_ar,X_ar,Y_ar,T_ar,nbulk,
     &  nll_in,np_in,iherm,i_deriv,
     &  K_sum,d_ono,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)
     
c       compute the derivatives from the disperstion function D
c       for the complete hot tensor k_sum
c       with respect dielectric tensor components: ddeps
c       with respect N_perpendicular ,N_parallel : ddnp,ddnll)

        call dddeps(k_sum,nll,nperp,ddeps,dd_ono_nll,dd_ono_np)
      endif

c-----initialization
     
      do j=1,3
         do js=1,nbulk
            dd_ono(j,js)=dcmplx(0.d0,0.d0)
         enddo
      enddo

c---------------------------------------------------------------------  
c     Now, compute the derivative of the D(w,k) with respesct
c     all 6 variables for Hermitian and complete tensor
c     Sum_j{i1=1,3, i2=1,3}[dD/deps(i1,i2)*deps_js(i1,i2)/d(X_js,Y_js,T_js,nper,npar]
c     the first 5 derivatives (j=1,..,5) are calculated from Hermitian D,
c     the derivative dd_ono_pn=dD/d_Nperp is calculated from the complete D
                   
      do i1=1,3
         do i2=1,3
            dd_ono_np_h =dd_ono_np_h +ddeps_h(i1,i2)*dK_dnper_h(i1,i2)      
            dd_ono_nll_h=dd_ono_nll_h+ddeps_h(i1,i2)*dK_dnpar_h(i1,i2)

            if(i_dd_ono_np.eq.1) then
c-------------complete tensor derivative dD/DN_perp 
              dd_ono_np   =dd_ono_np   +ddeps(i1,i2)*dK_dnper(i1,i2)
            endif

            do js=1,nbulk
               dd_ono(1,js)=dd_ono(1,js)+
     &                      ddeps_h(i1,i2)*dK_dx_ar_h(i1,i2,js)
               dd_ono(2,js)=dd_ono(2,js)+
     &                      ddeps_h(i1,i2)*dK_dy_ar_h(i1,i2,js)
               dd_ono(3,js)=dd_ono(3,js)+
     &                      ddeps_h(i1,i2)*dK_dt_ar_h(i1,i2,js)
            enddo !js

         enddo
      enddo

c      write(*,*)'dDOno dd_ono_np_h',dd_ono_np_h 
      return
      end

      subroutine DdOno1(nbulk,mass_ar,x_ar,y_ar,t_ar,
     .nll_in,np_in,dd_ono,dd_ono_np_h,dd_ono_nll_h,
     .i_dd_ono_np,dd_ono_np)
c
c     This function returns the derivatives of the Hermitian
c     part of the Ono dispersion function (used for fast waves)
c     with respect to X,Y,Te_av,nll_in, and np_in
c
c     INPUTS:
c      nbulk the total number of plasma species
c      mass_ar - the masses  of the plasma species (in electron mass) 
c      x_ar = (fpe/f)**2
c      y_ar = fce/f   It is the algebraic number
c      t_av_ar=average temperature in ev    
c      nll_in - parallel index of refraction n.
c      np_in - perpendicular index of refraction n.
c      i_dd_ono_np =0 no calcultions of dd_ono_np (dD/dN_perp for full tensor)
c                  =1 calculates        dd_ono_np (dD/dN_perp for full tensor)   
c     OUTPUT:
c      DD_ono(1-3,nbulk) derivatives from the dispersion function (for hot
c                  Hermitian dielectric tensor ) with respect
c                  X_js,Y_js,T_av_js (js=1,nbulk)
c      DD_ono(1,js)=dD/dX_js
c      DD_ono(2,js)=dD/dY_js
c      DD_ono(3,js)=dD/dT_av_js
c
c      derivatives from the dispersion function (for hot
c      Hermitian dielectric tensor ) with respect
c      dd_ono_np_h  = dD/DN_perpendicular
c      dd_ono_nll_h = dD/dN_parallel 
c
c      derivative from the dispersion function (for hot
c      complete dielectric tensor ) with respect
c      dd_ono_np = dD/DN_perpendicular
c-----------------------------------------------------------
      Implicit none
cSm030226
      include 'param.i'
c-----input
      integer nbulk
      double precision mass_ar(*)
      double precision X_ar(*),Y_ar(*),T_ar(*)
      double precision np_in,nll_in
      integer i_dd_ono_np
c-----external
c-----output
c     derivatives from the Hermitian dispersion function D
c     dd_ono(1,js)=dD/dX_js,dd_ono(2,js)=dD/dY_js,dd_ono(3,js)=dD/dT_av_js
c     dd_ono_np_h=dD/N_perp,dd_ono_nll_h=dD/dN_parallel
c     derivartive from the complete dispersion function 
c     dd_ono_np=dD/N_erp
      double complex dd_ono(3,*) !(3,nbulka) 
      double complex dd_ono_nll_h,dd_ono_np_h
      double complex dd_ono_np
 
c-----locals
      double precision nperp,nll,nlls,nps
      double complex k_herm(3,3),k_sum(3,3)
      integer ncomp_l
cSm030226
c      parameter (ncomp_l=4)
      parameter (ncomp_l=nbulka)
      double complex d_ono,d_ono_h,dd_ono_nll,
     &dK_dx_ar(3,3,ncomp_l),dK_dy_ar(3,3,ncomp_l),dK_dt_ar(3,3,ncomp_l),
     &dK_dnpar(3,3),dK_dnper(3,3),
     &dK_dx_ar_h(3,3,ncomp_l),dK_dy_ar_h(3,3,ncomp_l),
     &dK_dt_ar_h(3,3,ncomp_l),dK_dnpar_h(3,3),dK_dnper_h(3,3)

      double complex ddeps(3,3),ddeps_h(3,3)
      double complex dk(3,3,6)
      integer i_deriv,iherm
      double complex dp,dp1,dp2,dp6
      double complex ddp(6)
      integer js,j,i1,i2

      if(ncomp_l.lt.nbulk) then
        write(*,*)'ono_disp.f Dd_Ono ncomp_l=',ncomp_l,'nbulk=',nbulk
        write(*,*)'it should be ncompl_l.ge.nbulk'
        write(*,*)'ncomp_l=nbulka'
        write(*,*)'change parameter nbulka in param.i recompile'
        stop
      endif
       
c-----put the minimal values for nll,nperp
      call npnllmin(nll_in,np_in,nll,nperp)

      nlls=nll*nll
      nps=nperp*nperp

      i_deriv=2 ! Hermitian tensor and its derivatives
      iherm=1     
      call ono_tens(mass_ar,X_ar,Y_ar,T_ar,nbulk,
     & nll_in,np_in,iherm,i_deriv,
     & k_herm,d_ono_h,dK_dx_ar_h,dK_dy_ar_h,dK_dt_ar_h,
     & dK_dnper_h,dK_dnpar_h)

c      compute the derivatives from the disperstion function D
c     for the hermitian hot tensor k_herm
c     with respect dielectric tensor components: ddeps_h
c     with respect N_perpendicular ,N_parallel : ddnp_h,ddnll_h)

      call dddeps(k_herm,nll,nperp,ddeps_h,dd_ono_nll_h,dd_ono_np_h)

      if(i_dd_ono_np.eq.0) then
        dd_ono_np=dcmplx(0.d0,0.d0)
      endif
 
      if(i_dd_ono_np.eq.1) then
        i_deriv=2 ! complete tensor and its derivatives
        iherm=2     
        call ono_tens(mass_ar,X_ar,Y_ar,T_ar,nbulk,
     &  nll_in,np_in,iherm,i_deriv,
     &  K_sum,d_ono,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)
     
c       compute the derivatives from the disperstion function D
c       for the complete hot tensor k_sum
c       with respect dielectric tensor components: ddeps
c       with respect N_perpendicular ,N_parallel : ddnp,ddnll)

        call dddeps(k_sum,nll,nperp,ddeps,dd_ono_nll,dd_ono_np)
      endif

c-----initialization
      dp1=dcmplx(0.d0,0.d0)
      dp2=dcmplx(0.d0,0.d0)
      dp6=dcmplx(0.d0,0.d0)
 
                   
      do js=1,nbulk ! the loop over the plasma species
c--------calculations of the derivatives dK(3,3,6) from the
c        sucseptibilities K_js for js specie (i1=1,2,3, i2=1,2,3):
c        dK(i1,i2,1)_js=dK(i1,i2)_js/dN_perp       dK_nper
c        dK(i3,i2,2)_js=dK(i1,i2)_js/dN_parallel   dk_npar
c        dK(i1,i2,3)_js=dK(i1,i2)_js/dX_js
c        dK(i1,i2,4)_js=dK(i1,i2)_js/dY_js
c        dK(i3,i2,5)_js=dK(i1,i2)_js/dT_js
c        dK(i1,i2,3,jj=1-5) for Hermitian part of the tensor
c        dK(i1,i2.jj=6)=dD/dN_perp for the complete tensor

         do i1=1,3
           do i2=1,3
             dk(i1,i2,1)=dK_dnper_h(i1,i2)
             dk(i1,i2,2)=dK_dnpar_h(i1,i2)
             dk(i1,i2,3)=dK_dx_ar_h(i1,i2,js)
             dk(i1,i2,4)=dK_dy_ar_h(i1,i2,js)
             dk(i1,i2,5)=dK_dt_ar_h(i1,i2,js)
             dk(i1,i2,6)=dK_dnper(i1,i2)
           enddo
         enddo      

         do j=1,6
c----------Now, compute the derivative of the D(w,k) with respesct
c          all 6 variables  
c          for Hermitian and complete (j=8) tensor
c          Sum_j{i1=1,3, i2=1,3}[dD/deps(i1,i2)*deps_js(i1,i2)/d(variavble_j)]
c          variable_j=(j=1 n_perp),(j=2 n_par),(j=3 X_js),(j=4 Y_js),
c          (j=5 T_av_js), (j=6 n_perp)
c          the first 5 derivatives (j=1,..,5) are calculated from Hermitian D,
c          the last 6th derivative is calulated from the complete D

           dp=dcmplx(0.d0,0.d0) ! initialization     
                                 
	   do i1=1,3
	     do i2=1,3
	       if (j.ne.6) then
	          dp=dp+ddeps_h(i1,i2)*dk(i1,i2,j)
               else
	          dp=dp+ddeps(i1,i2)*dk(i1,i2,6)
	       endif
	     enddo
	   enddo
 
           if ((j.ge.3).and.(j.le.5)) dd_ono(j-2,js) = dp

c          summation by all js=1,...,nbulk species 
           if  (j.eq.1) dp1=dp1+dp 
	   if  (j.eq.2) dp2=dp2+dp 
	   if  (j.eq.6) dp6=dp6+dp 
         enddo ! j
      enddo  ! js

      dd_ono_np_h  = dd_ono_np_h  + dp1     ! DD_hermitian/DN_perp
      dd_ono_nll_h = dd_ono_nll_h + dp2     ! DD_hermitian/DN_parallel
      dd_ono_np    = dd_ono_np    + dp6     ! DD_complete/DN_perp works for i_dd_ono_np=1  
 
      return
      end



      subroutine ono_dervs(nbulk,u,wf,dddcnz,dddcnr,dddcm,
     .dddz,dddr,dddph,dddw)
c---------------------------------------------------------
c     analytical calculation of the derivatives for the ray-tracing
c     equations from the ono plasma
c     with ions and electrons for fast waves
c--------------------------------------------------------- 
      implicit none
cSm030226
      include 'param.i'
c-----input
      double precision u(6)  ! z,r,phi,Nz,Nr,M
      double precision wf    ! wave friquency
      integer nbulk          ! the total number of the plasma components
c-----output
      double precision dddcnz,dddcnr,dddcm,
     .dddz,dddr,dddph,dddw   !derivatives from D
c-----uses
      external ono_tens,Ddono,b,x,y,rhof,tempe,
     .dxdz,dxdr,dxdphi,dydz,dydr,dydphi,
     .temperho,dtempdz,dtempdr
      double precision b,x,y,rhof,tempe,
     .temperho,
     .dxdz,dxdr,dxdphi,dydz,dydr,dydphi,
     .dtempdz,dtempdr

c-----locals
      integer ncomp_l
cSm030226
c      parameter (ncomp_l=4)
      parameter (ncomp_l=nbulka)
      double complex dd_ono(3,ncomp_l) !dD/d(X_s,Y_s,T_s)

      double precision t_kev,nll,nperp
      double complex K_sum(3,3),dK(3,3,5),
     &dd_ono_np_h,dd_ono_nll_h,dd_ono_np,d
      double precision z,r,phi,cnz,cnr,cm,bmod,rholoc
      integer i,i_deriv
      double precision dxdze,dxdre,dxdphie,dydze,dydre,dydphie, 
     .dtempdze,dtempdre,
     .dtempdpe,
     .dxdwe,dydwe,
     .dnpdz,dnpdr,dnpdphi,dnpdcnz,dnpdcnr,dnpdcm,
     .dnlldz,dnlldr,dnlldphi,dnlldcnz,dnlldcnr,dnlldcm,
     .dnpdw,dnlldw
 
      double precision mass_ar(ncomp_l),x_ar(ncomp_l),y_ar(ncomp_l),
     .t_av_ar(ncomp_l)
      double complex dK_dx_ar(3,3,ncomp_l),dK_dy_ar(3,3,ncomp_l),
     &dK_dt_ar(3,3,ncomp_l),dK_dnper(3,3),dK_dnpar(3,3)

cfor test only
      double precision step,zp,zm,rholocp,xep,yep,te_keVp,tep,
     .tpopep,vflowep,rholocm,xem,yem,te_keVm,tem,
     .tpopem,vflowem,rp,rm,phip,phim,cnzp,cnzm,cnrp,cnrm,cmp,cmm,
     .nllp,nllm,cnpp,cnpm
      double complex dp,dm,ddnum(5),deld
      double precision rdp,rdm,rdeld,rddnum(5)
      integer j,iherm
      double precision ddtest(5),ktest(3,3),ddntest
     
      double precision x_ar_i_old,y_ar_i_old,t_ar_i_old,np_old,nll_old
      double complex K_sum_p(3,3),K_sum_m(3,3),deriv
      integer i1,i2
      
      step=1.d-7
      
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)

      bmod=b(z,r,phi)
      rholoc=rhof(z,r,phi)

      if(nbulk.gt.ncomp_l) then
	   write(*,*)'in ono_disp.f in ono_dervs nbulk.gt.ncomp_l'
	   write(*,*)'put the value of parameter ncomp_l.ge.nbulk'
	   write(*,*)'parameter ncomp=nbulka is given in param.i '
           write(*,*)'change nbulka in param.i and recompile'
           write(*,*)'ncomp_l,nbulk',ncomp_l,nbulk           
	   stop
      endif

c-----The initilization mass_ar
      call put_mass(mass_ar,nbulk)
      
c-----nll=(cnz*bz+cnr*br+cm*bphi/r)/bmod
c     nperp=dsqrt(cnz**2+cnr**2+(cm/r)**2-nll**2)
c     calculation nll=N_parallel and nperp=N_perp
c     for given z,r,phi,cnz,cnr,cm,
        
      call nllcnp_s(z,r,phi,cnz,cnr,cm,nll,nperp)
     
       do i=1,nbulk
          x_ar(i)=x(z,r,phi,i)
          y_ar(i)=y(z,r,phi,i)
          t_keV=tempe(z,r,phi,i)   !keV averaged temperature
          t_av_ar(i)=t_keV*1.d3    !eV
       enddo
	
c------numerical derivatives deps/dx
c       write(*,*)'ono_dervs'      
c       do i=1,nbulk
c        x_ar_i_old=x_ar(i)
c        X_ar(i)=x_ar_i_old*(1.d0+step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_p,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        X_ar(i)=x_ar_i_old*(1.d0-step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_m,dm,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        X_ar(i)=x_ar_i_old
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,1,
c     &  K_sum,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        do i1=1,3
c          do i2=1,3
c          deriv=(K_sum_p(i1,i2)-K_sum_m(i1,i2))/(2.d0*step*x_ar_i_old)
c             write(*,*)'i,i1,i2,dK_dx_ar(i1,i2,i),deriv',
c     &       i,i1,i2,dK_dx_ar(i1,i2,i),deriv            
c          enddo  
c        enddo
c       enddo !i

c------numerical derivatives deps/dy      
c       do i=1,nbulk
c        y_ar_i_old=y_ar(i)
c        y_ar(i)=y_ar_i_old*(1.d0+step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_p,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        y_ar(i)=y_ar_i_old*(1.d0-step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_m,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        Y_ar(i)=y_ar_i_old
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,1,
c     &  K_sum,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)


c        write(*,*)'x_ar(i),y_ar(i),step',x_ar(i),y_ar(i),step
c        do i1=1,3
c          do i2=1,3
c          deriv=(K_sum_p(i1,i2)-K_sum_m(i1,i2))/(2.d0*step*y_ar_i_old)
c             write(*,*)'i,i1,i2,dK_dy_ar(i1,i2,i),deriv',
c     &       i,i1,i2,dK_dy_ar(i1,i2,i),deriv            
c          enddo  
c        enddo
c      enddo !i

c------numerical derivatives deps/dT      
c       do i=1,nbulk
c        t_ar_i_old=T_av_ar(i)
c        T_av_ar(i)=t_ar_i_old*(1.d0+step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_p,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        T_av_ar(i)=t_ar_i_old*(1.d0-step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_m,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        T_av_ar(i)=y_ar_i_old
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,1,
c     &  K_sum,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)


c       write(*,*)'x_ar(i),y_ar(i),t_av_ar(i)',x_ar(i),y_ar(i),t_av_ar(i)
c        do i1=1,3
c          do i2=1,3
c          deriv=(K_sum_p(i1,i2)-K_sum_m(i1,i2))/(2.d0*step*t_ar_i_old)
c             write(*,*)'i,i1,i2,dK_dt_ar(i1,i2,i),deriv',
c     &       i,i1,i2,dK_dt_ar(i1,i2,i),deriv            
c          enddo  
c        enddo
c      enddo !i

c------numerical derivatives deps/dN_par      
c        nll_old=nll
c        write(*,*)'ono_dervs nll',nll
c        nll=nll_old*(1.d0+step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_p,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)
        
c        nll=nll_old*(1.d0-step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_m,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        nll=nll_old
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,1,
c     &  K_sum,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        write(*,*)'K_sum_p(1,3),K_sum_m(1,3)',K_sum_p(1,3),K_sum_m(1,3)
c        write(*,*)'step,nll_old',step,nll_old
c        do i1=1,3
c          do i2=1,3
c          deriv=(K_sum_p(i1,i2)-K_sum_m(i1,i2))/(2.d0*step*nll_old)
c             write(*,*)'i1,i2,dK_dnpar(i1,i2),deriv',
c     &       i1,i2,dK_dnpar(i1,i2),deriv            
c          enddo  
c        enddo

c------numerical derivatives deps/dN_par      
c        np_old=nperp
c        nperp=np_old*(1.d0+step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_p,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        nperp=np_old*(1.d0-step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_m,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        nperp=np_old
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,1,
c     &  K_sum,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        do i1=1,3
c          do i2=1,3
c          deriv=(K_sum_p(i1,i2)-K_sum_m(i1,i2))/(2.d0*step*np_old)
c             write(*,*)'i1,i2,dK_dnper(i1,i2),deriv',
c     &       i1,i2,dK_dnper(i1,i2),deriv            
c          enddo  
c        enddo

ctestend

	do i=1,nbulk                  
          x_ar(i)=x(z,r,phi,i)
          y_ar(i)=y(z,r,phi,i)
          t_keV=tempe(z,r,phi,i)   !keV averaged temperature
          t_av_ar(i)=t_keV*1.d3    !eV	
	enddo

c        write(*,*)'ono_dervs before call ono_tens'
        
        call ono_tens(mass_ar,X_ar,Y_ar,t_av_ar,nbulk,
     &  nll,nperp,1,0,
     &  K_sum,d,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        write(*,*)'ono_dervs after ono_tens'
c	d= dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
c     .vflow_ar,nll,nperp,1,K_sum)

        call DdOno(nbulk,mass_ar,x_ar,y_ar,t_av_ar,
     .  nll,nperp,dd_ono,dd_ono_np_h,dd_ono_nll_h,
     .  0,dd_ono_np)

ctest
c------numerical derivatives dD/dX
c       do i=1,nbulk      
c        x_ar_i_old=x_ar(i)
c        X_ar(i)=x_ar_i_old*(1.d0+step)

c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_p,dp,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        X_ar(i)=x_ar_i_old*(1.d0-step)
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,0,
c     &  K_sum_m,dm,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        X_ar(i)=x_ar_i_old
c        call ono_tens(mass_ar,X_ar,Y_ar,T_av_ar,nbulk,
c     &  nll,nperp,1,1,
c     &  K_sum,d,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)

c        deriv=(dp-dm)/(2.d0*step*x_ar_i_old)
c        write(*,*)'i,dd_ono(1,i),deriv',i,dd_ono(1,i),deriv            
c       enddo !i
cendtest

      call dnd(z,r,phi,cnz,cnr,cm,
     .dnpdz,dnpdr,dnpdphi,dnpdcnz,dnpdcnr,dnpdcm,
     .dnlldz,dnlldr,dnlldphi,dnlldcnz,dnlldcnr,dnlldcm,
     .dnpdw,dnlldw)     
      dddz = dd_ono_np_h*dnpdz+dd_ono_nll_h*dnlldz     
      dddr = dd_ono_np_h*dnpdr+dd_ono_nll_h*dnlldr
      dddph= dd_ono_np_h*dnpdphi+dd_ono_nll_h*dnlldphi
      dddcnz=dd_ono_np_h*dnpdcnz+dd_ono_nll_h*dnlldcnz
      dddcnr=dd_ono_np_h*dnpdcnr+dd_ono_nll_h*dnlldcnr
      dddcm =dd_ono_np_h*dnpdcm +dd_ono_nll_h*dnlldcm
      
      dddw=dd_ono_np_h*dnpdw+dd_ono_nll_h*dnlldw

      do i=1,nbulk

        dxdze=dxdz(z,r,phi,i)
        dxdre=dxdr(z,r,phi,i)
        dxdphie=dxdphi(z,r,phi,i)

        dydze=dydz(z,r,phi,i)
        dydre=dydr(z,r,phi,i)
        dydphie=dydphi(z,r,phi,i)
      
        dtempdze=dtempdz(z,r,phi,i)
        dtempdre=dtempdr(z,r,phi,i)
        dtempdpe=0.d0			      !d(T_average)/d(phi)

        dddz = dddz+
     &  dd_ono(1,i)*dxdze+dd_ono(2,i)*dydze+dd_ono(3,i)*dtempdze*1.d3
        
        dddr = dddr+
     &  dd_ono(1,i)*dxdre+dd_ono(2,i)*dydre+dd_ono(3,i)*dtempdre*1.d3
        dddph= dddph+
     &  dd_ono(1,i)*dxdphie+dd_ono(2,i)*dydphie+dd_ono(3,i)*
     &  dtempdpe*1.d3
      
	dxdwe=-2.d0*x_ar(i)/wf
        dydwe=-y_ar(i)/wf
        dddw=dddw+dd_ono(1,i)*dxdwe+dd_ono(2,i)*dydwe

      enddo
     
      return
      end

 
