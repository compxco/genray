c     
      subroutine hamilt_eigenvalue_muz(z,r,phi,cnpar,cnper,cnprim,
     &k_root,D,Dc)   
c-----calculate the Hamiltonian from mazzucato tensor
c     using D=real part(eigenvalue(k) of the dispersion tensor)
c     Dc=eigenvalue(k_root), k_root=1,2,3

      implicit none
      include 'eps.i' !reps(3,3) complex dielectric tensor
c-----input    
      real*8 z,r,phi,! space coordinates. 
c         Radius rho should be calculated before call hamilt_eigenvalue
c         in subroutine b()    
     &cnpar,cnper,cnprim!Re(N_parallel),Re(N_perpendicular),
                        !Im(N_perpendicular)
      integer k_root !the number of the eigenvalue, k_rrot=1,2,3
c-----output
      real*8 D      !dispersion function=Re(eigenvalue(k))
      complex*16 Dc !Dc=eigenvalue(k_root), k_k_root=1,2,3

c-----locals
      integer ihermloc
      real*8 hamiltmz
      complex*16 cnper_c,cnpar_c,eigenvalue(3)
      logical check 
      real*8 min            
      integer k

      check=.FALSE.
      check=(((k_root.eq.1).or.(k_root.eq.2)).or.(k_root.eq.3))
c      write(*,*)'hamilt_eigenvalue k_root ',k_root,'check ',check

      if (.not.check) then
         write(*,*)'in subroutine hamilt_eigenvalue'
         write(*,*)'kroot=',k_root
         write(*,*)'But it should be k_root=1,2,3'
         stop 'in subroutine hamilt_eigenvalue'
      endif

c---- calculate the mazzucato complex dielectric tensor reps(3,3)
c     This tensor will be in include file 'eps.i'.

c      write(*,*)'hamilt_eigenvalue z,r,phi,cnpar,cnper,cnprim,k_root',
c     & z,r,phi,cnpar,cnper,cnprim,k_root
        
      ihermloc=2 !for non-hermitian tensor
      call hamiltmuz_ful(cnpar,ihermloc,z,r,phi,cnper,cnprim,
     &hamiltmz)
 
c     write(*,*)'hamilt_eigenvalue after hamilt_muz_ful hamiltmz'
c     &,hamiltmz

c-----calculate the complex eigenvalues of the dispertion tensor
      cnper_c=cmplx(cnper,cnprim)
      cnpar_c=cmplx(cnpar,0.d0)         

      call lambda(cnper_c,cnpar_c,reps,eigenvalue)

c      write(*,*)'in hamilt_eigenvalue k_root',k_root,
c     &' eigenvalue',eigenvalue
  
c-----calculate Hamiltonian: D=Real(eigenvalue(k_root))
      Dc=eigenvalue(k_root)
      D=dreal(Dc)
      D=dreal(eigenvalue(1))*dreal(eigenvalue(2))*dreal(eigenvalue(3))
c-----choose eigenvalue(k) that gives min=dabs(dreal(eigenvalue(k)))
      min=1.d10
      do k=1,3
       if(dabs(dreal(eigenvalue(k))).lt.min) then
         min=dabs(dreal(eigenvalue(k)))
         Dc=eigenvalue(k)
         D=dreal(Dc)
       endif
      enddo 
 10   continue
c      call d0_lambda(cnper_c,cnpar_c,reps,Dc)
c      D=dreal(Dc)

      return
      end



      subroutine lambda(cnper,cnpar,eps,eigenvalue)
c-----Calculates three eigenvalues of the 
c     dispersion tensor eigenvalue

      implicit none

c-----input
      complex*16 eps(3,3),  ! dielectric tensor
     &cnper,                ! N_parallel
     &cnpar                 ! N_perpendicular 

c-----output
      complex*16 eigenvalue(3) !eigenvalues of the dispersion tensor

c-----local
      complex*16 a0,b0,c0,d0,! cubic equation coefficients
     &x1,x2,x3               ! three roots:the eigenvalue lambda
      
c-----calculate cubic equation coefficients
c     a0*lambda**3+b0*lambda**2+c0*lambda+d0=0

      call set_cubic_eq_coef(eps,cnpar,cnper,a0,b0,c0,d0)
c-----solve the cubic equation. The roots are x1,x2,x3 

c      write(*,*)'lambda a0,b0,c0,d0', a0,b0,c0,d0

      call cubic_solver(a0,b0,c0,d0,x1,x2,x3)
      eigenvalue(1)=x1
      eigenvalue(2)=x2
      eigenvalue(3)=x3
 
c      write(*,*)'lambda x1,x2,x3',x1,x2,x3
c      write(*,*)'f(x1)',a0*x1**3+b0*x1**2+c0*x1+d0
c      write(*,*)'f(x2)',a0*x2**3+b0*x2**2+c0*x2+d0
c      write(*,*)'f(x3)',a0*x3**3+b0*x3**2+c0*x3+d0
      return
      end

      subroutine cubic_solver(a0,b0,c0,d0,x1,x2,x3)
c-----solves cubic equation
c     a0*x**3+b0*x*2+c0*x+d0=0   (1)
      implicit none
c-----input
      complex*16 a0,b0,c0,d0
c-----output
      complex*16 x1,x2,x3 !roots

c-----externals
      complex*16 cdsqrt_Sm

c-----locals
      complex*16 i_c
      complex*16 a,b,c,p,q,Qc,sqrtQc,Ac3,Bc3,
     &Ac_array(3),Bc_array(3),Ac,Bc,
     &y_1,y_2,y_3
      integer k,l,k0,l0
      real*8 three,pnorm

      complex*16 zero

      zero=(0.0,0.0)
      i_c=(0.0,1.0)
      three=3.0
 
      if (a0.eq.zero) then
        write(*,*)'in cubic_root a0=0'
        goto 10
      endif
 
c-----solve the cubic equaton
c     x**3+a*x**2+b*x+c=0      (2)
      a=b0/a0
      b=c0/a0
      c=d0/a0
c      write(*,*)'cubic solver a,b,c',a,b,c

c-----using new variabe x=y-a/3 in (2) we have
c     y**3+p*y+q=0  
      p=-a**2/3+b
      q=2*(a/3)**3-a*b/3+c
      
      Qc=(p/3)**3+(q/2)**2
      
c      sqrtQc=cdsqrt(Qc)
      sqrtQc=cdsqrt_Sm(Qc)

c      write(*,*)'q/2,Qc,sqrtQc',q/2,Qc,sqrtQc
      Ac3=-q/2+sqrtQc  
      Bc3=-q/2-sqrtQc
           
c-----calculate three complex third order roots 
c     Ac_array=(Ac3)**1/3,      Bc_array=(Bc3)**1/3 
      call qubic_complex_roots(Ac3,Ac_array)
c      write(*,*)'Ac3,Ac_array',Ac3,Ac_array
      call qubic_complex_roots(Bc3,Bc_array)
c      write(*,*)'Bc3,Bc_array',Bc3,Bc_array


c-----choose of one complex root from the condition
c     Ac_array(k)*Bc_array(l)=-p/3
      pnorm=10.d0**10
      do k=0,2
         do l=0,2
c          write(*,*)'k,l,Ac_array(k+1)*Bc_array(l+1)+p/3',
c     &    k,l,Ac_array(k+1)*Bc_array(l+1)+p/3.0

c          if(abs(Ac_array(k+1)*Bc_array(l+1)+p/3.0).lt.1.d-13)then
c          if(abs(Ac_array(k+1)*Bc_array(l+1)+p/3.0).lt.1.d-12)then
          if(abs(Ac_array(k+1)*Bc_array(l+1)+p/3.0).lt.
     &        (abs(Ac_array(k+1))*1.d-11)) then
           k0=k
           l0=l
           goto 20
          endif

          if(abs(Ac_array(k+1)*Bc_array(l+1)+p/3.0).lt.pnorm)then
             pnorm=abs(Ac_array(k+1)*Bc_array(l+1)+p/3.0)
             k0=k
             l0=l             
          endif  
         enddo
      enddo
 20   continue
c      write(*,*)'k0,l0,abs(Ac_array(k0+1)*Bc_array(l0+1)+p/3.0)',
c     &k0,l0,abs(Ac_array(k0+1)*Bc_array(l0+1)+p/3.0)
      Ac=Ac_array(k0+1)
      Bc=Bc_array(l0+1)

      y_1=Ac+Bc
      y_2=-y_1/2.0+i_c*(Ac-Bc)*dsqrt(three)/2.0
      y_3=-y_1/2.0-i_c*(Ac-Bc)*dsqrt(three)/2.0

      x1=y_1-a/3.0
      x2=y_2-a/3.0
      x3=y_3-a/3.0

 10   continue
      return
      end

      subroutine qubic_complex_roots(x_cube,x)
c-----calculate complex cubic roots from the complex argument x_cubic
c     x(k+1)=x_cube_mod**1/3*(cos(phi+2*k*pi)/3)+i*sin((phi+2*k*pi)/3)
c     k=0,1,2
c     x_cube_mode=sqrt(Re(|x_cube|)**2+Im(|x_cube|**2))
      implicit none
c-----input
      complex*16 x_cube
c-----output
      complex*16 x(3)
c-----locals
      real*8 x_cube_mod,   !x_cube module
     &phi,                 !x_cube argument
     &pi,one
      complex*16 i_c
      integer k

      one=1.d0
      pi=4.d0*atan(one)
      i_c=(0.0,1.0)

      x_cube_mod=abs(x_cube)

      if (x_cube_mod.eq.0.) then
         do k=1,3
            x(k)=(0.d0,0.d0)
         enddo
         goto 10
      else
cBH060730         phi=atan2(imag(x_cube),dreal(x_cube))
         phi=atan2(aimag(x_cube),dreal(x_cube))
      endif

      x_cube_mod=x_cube_mod**(1.0d0/3.0d0) 
 
      do k=0,2
        x(k+1)=x_cube_mod*(cos((phi+2*k*pi)/3.0)
     &                +i_c*sin((phi+2*k*pi)/3.0))  
      enddo
      
 10   continue
      return
      end
       

      subroutine set_cubic_eq_coef(eps,cnpar,cnper,a0,b0,c0,d0)
c-----calculate the coefficient of the cubic equation
c     Det(D_mp-lambda(delat_mp)=a0*lambda**3+b0*labda**2+c0*lambda+d0=0
c     Here D_mp is the dispersion tensor
c          lambda is the dispersion tensor eigenvalue
c     D_mp-labmda*delta_mp=
c     |eps_xx-N_par*2-lambda,  eps_xy,             eps_xz+N_perp*N_par     |
c     |epx_yx,                 eps_yy-N**2-lambda, eps_yz                  |
c     |eps_zx+N_perp*N_par,    eps_zx,             eps_zz-N_perp**2-lambda |
c     
c     eps_mp is the dielectric tensor
c
c
c     Det=A+B+C
c     A=(eps_xx-N_par*2-lambda)*
c       [(eps_yy-N**2-lambda)*(eps_zz-N_perp**2-lambda)-eps_yz *eps_zx]
c
c     B=-eps_xy*[epx_yx*(eps_zz-N_perp**2-lambda)-eps_yz*(eps_zx+N_perp*N_par]
c
c     C=(eps_xz+N_perp*N_par)*
c       [epx_yx*eps_zx-(eps_yy-N**2-lambda)*(eps_zx+N_perp*N_par)]
c
c-----input
      complex*16 eps(3,3),  !dielectric tensor
     &cnper,                ! N_parallel
     &cnpar                 ! N_perpendicular
c-----output
      complex*16 a0,b0,c0,d0 !coefficients of the cubic equation     

c-----local
      complex*16 cnpar2,cnper2,cn2

      cnpar2=cnpar*cnpar
      cnper2=cnper*cnper
      cn2=cnper2+cnpar2

c      write(*,*)'cnper,cnpar,cnper2,cnpar2,cn2',
c     &cnper,cnpar,cnper2,cnpar2,cn2
c      write(*,*)'eps',eps

      a0=-1.                                                !A
        
      b0=(eps(1,1)-cnpar2)+(eps(2,2)-cn2+eps(3,3)-cnper2)   !A
      

      c0=-(eps(1,1)-cnpar2)*(eps(2,2)-cn2+eps(3,3)-cnper2)
     &   -(eps(2,2)-cn2)*(eps(3,3)-cnper2)                  !A
     &   +eps(2,3)*eps(3,2)
c-----------------------------------------------------------
     &   +eps(1,2)*eps(2,1)                                 !B
     &   +(eps(1,3)+cnpar*cnper)*(eps(3,1)+cnpar*cnper)     !B

      d0=(eps(1,1)-cnpar2)*(eps(2,2)-cn2)*(eps(3,3)-cnper2) !A
     &   -(eps(1,1)-cnpar2)*eps(2,3)*eps(3,2)
c----------------------------------------------------------
     &   -eps(1,2)*                                         !B
     &   (eps(2,1)*(eps(3,3)-cnper2)-eps(2,3)*(eps(3,1)+cnpar*cnper))

     &   +(eps(1,3)+cnpar*cnper)*                           !C
     &  (eps(2,1)*eps(3,2)-(eps(2,2)-cn2)*(eps(3,1)+cnpar*cnper))

      return
      end

      complex*16 function cdsqrt_Sm(z)
c-----sq_root from complex z=rho*exp(i*phi)
c     cdxqrt_Sm=dsqrt(rho)*eps(i*phi/2)

      implicit none
c-----input
      complex*16 z
     
c-----local
      real*8 Re_z,Im_z,phi,abs_z,pi,p

      pi=4.d0*datan(1.d0)
            
      Re_z=dreal(z)
      Im_z=dimag(z)
       
      abs_z=dsqrt(Re_z**2+Im_z**2)

      if (abs_z.eq.0.d0) then
        cdsqrt_Sm=dcmplx(0.d0,0.d0)
        goto 10
      endif 

      if (Im_z.ge.0.d0) then
c--------0 =< phi =< pi

         if(Re_z.ge.0.d0) then
c---------- 0 =< phi =< pi/2
            p=Im_z/abs_z
            if(p.le.1.d0) then
               phi=dasin(p)
            else
               phi=0.5d0*pi
            endif
         else
c---------- pi/2 =< phi =< pi
            p=Im_z/abs_z
            if(p.le.1.d0) then
               phi=pi-dasin(p)
            else
               phi=0.5d0*pi
            endif
         endif

      else
c--------pi =< phi =< 2*pi
         if(Re_z.ge.0.d0) then
c---------- 3*pi/2 =< phi < 2*pi
            p=Im_z/abs_z
            if(p.gt.-1.d0) then
               phi=2.d0*pi+dasin(p)
            else
               phi=1.5d0*pi
            endif
         else
c---------- pi < phi =< 3*pi/2
            p=Im_z/abs_z
            if(p.gt.-1.d0) then
               phi=pi-dasin(p)
            else
               phi=1.5d0*pi
            endif
         endif
      endif

      cdsqrt_Sm=dsqrt(abs_z)*dcmplx(dcos(0.5d0*phi),dsin(0.5d0*phi))
 10   continue

      return
      end 

c


      subroutine d0_lambda(cnper,cnpar,eps,d0)
c-----Calculates the coefficient d0
c     from the cubic equation  
c     a0*lambda**3+b0*lambda**2+c0*lambda+d0=0 
c     for the dispersion tensor eigenvalue

      implicit none

c-----input
      complex*16 eps(3,3),  ! dielectric tensor
     &cnper,                ! N_parallel
     &cnpar                 ! N_perpendicular 

c-----output
      complex*16 d0         !cubic equation coefficient

c-----local
      complex*16 a0,b0,c0   ! cubic equation coefficients
      
c-----calculate cubic equation coefficients
c     a0*lambda**3+b0*lambda**2+c0*lambda+d0=0

      call set_cubic_eq_coef(eps,cnpar,cnper,a0,b0,c0,d0)

c      write(*,*)'lambda d0',d0
      return
      end

      subroutine hamilt_eigenvalue_hot(z,r,phi,cnpar,cnper,cnprim,
     &k_root_loc,D,Dc)
   
c-----calculate the Hamiltonian from  Stix hot plasma tensor
c     using D=real part(eigenvalue(k) of the dispersion tensor)
c     Dc=eigenvalue(k_root_loc), k_root=1,2,3
c     complex*16 dielectric tensor reps(3,3) wil be in eps.i 
  
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'eps.i'
  
c-----input    
      real*8 z,r,phi,! space coordinates. 
c     Radius rho should be calculated before call hamilt_eigenvalue
c     in subroutine b()    
     &cnpar,cnper,cnprim!Re(N_parallel),Re(N_perpendicular),
                        !Im(N_perpendicular)
      integer k_root_loc !the number of the eigenvalue, k_root=1,2,3
c     nbulk is the number of plasma components. It was set in one.i
c-----output
      real*8 D      !dispersion function=Re(eigenvalue(k))
      complex*16 Dc !Dc=eigenvalue(k_root), k_k_root=1,2,3

c-----externals
      real*8 x,y,tempe,tpoprho,vflowrho,b
      complex*16 dhot_sum

c-----locals
      integer ihermloc,i
      complex*16 hamiltc
      complex*16 cnper_c,cnpar_c,eigenvalue(3)
             
      logical check 
      real*8 min            
      integer k

      real*8 x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     & tpop_ar(nbulka),vflow_ar(nbulka),te
      
      complex*16 compl_nper


      check=.FALSE.
      check=(((k_root_loc.eq.1).or.(k_root_loc.eq.2)).or.
     &       (k_root_loc.eq.3))
c      write(*,*)'hamilt_eigenvalue_hot k_root_loc ',k_root_loc,
c     &'check ',check

      if (.not.check) then
         write(*,*)'in subroutine hamilt_eigenvalue_hot'
         write(*,*)'kroot_loc=',k_root_loc
         write(*,*)'But it should be k_root_loc=1,2,3'
         stop 'in subroutine hamilt_eigenvalue_hot'
      endif

c---- calculate hot  non-relativistic plasma 
c     complex dielectric tensor reps(3,3)

c      write(*,*)'hamilt_eigenvalue_hot'
c      write(*,*)'z,r,phi,cnpar,cnper,cnprim,k_root_loc',
c     &z,r,phi,cnpar,cnper,cnprim,k_root_loc

      bmod=b(z,r,phi) !calculate small radius rho
                       !and put rho into one.i 

      do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpop_zrp(z,r,phi,i) !YuP[2024-08-14] was tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
      enddo

      compl_nper=dcmplx(cnper,cnprim)
      ihermloc=2 !for non-hermitian tensor

c-----calculate the hot plasma non-Hermitian tensor 
      hamiltc=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &   vflow_ar,cnpar,cnper,ihermloc,reps)
    
c      write(*,*)'hamilt_eigenvalue_hot after dhot_sum'
c      write(*,*)'hamilt_hot',hamiltc
      
c-----calculate the complex eigenvalues of the dispersion tensor
      cnper_c=cmplx(cnper,cnprim)
      cnpar_c=cmplx(cnpar,0.d0)         

      call lambda(cnper_c,cnpar_c,reps,eigenvalue)

c      write(*,*)'in hamilt_eigenvalue_hot k_root_loc',k_root_loc,
c     &' eigenvalue',eigenvalue
  
c-----calculate Hamiltonian: D=Real(eigenvalue(k_root_loc))
      Dc=eigenvalue(k_root_loc)
      D=dreal(Dc)
      D=dreal(eigenvalue(1))*dreal(eigenvalue(2))*dreal(eigenvalue(3))
c-----choose eigenvalue(k) that gives min=dabs(dreal(eigenvalue(k)))
      min=1.d10
      do k=1,3
       if(dabs(dreal(eigenvalue(k))).lt.min) then
         min=dabs(dreal(eigenvalue(k)))
         Dc=eigenvalue(k)
         D=dreal(Dc)
       endif
      enddo 
 10   continue

      return
      end

      subroutine hamilt_eigenvalue_combined (z,r,phi,cnpar,cnper,cnprim,
     &k_root,D,Dc)
   
c-----calculate the Hamiltonian from either
c     Eric Nelson Melby (id=12) or Abhay Ram (id=15) relativistic tensor
c     using D=real part(eigenvalue(k) of the dispersion tensor)
c     Dc=eigenvalue(k_root), k_root=1,2,3
c     complex*16 dielectric tensor reps(3,3) wil be in eps.i 
  
      implicit none
      include 'eps.i' !reps(3,3) complex dielectric tensor

c-----input    
      real*8 z,r,phi,! space coordinates. 
c     Radius rho should be calculated before call hamilt_eigenvalue
c     in subroutine b()    
     &cnpar,cnper,cnprim!Re(N_parallel),Re(N_perpendicular),
                        !Im(N_perpendicular)
      integer k_root !the number of the eigenvalue, k_rrot=1,2,3

c-----output
      real*8 D      !dispersion function=Re(eigenvalue(k))
      complex*16 Dc !Dc=eigenvalue(k_root), k_k_root=1,2,3

c-----externals
      real*8 x,y,tempe

c-----locals
      real*8 hamilt_combined
      complex*16 cnper_c,cnpar_c,eigenvalue(3)
    
      logical check 
      real*8 min            
      integer k

      real*8 X_e,Y_e,T_e
      complex*16 compl_nper

cSAP091216
c-----for muller solver
      real*8 errabs 
      integer nsig,nknown,nrts,nguess, nnew,itmax,ier,info(1)
      complex*16 roots(1)
c-----external
      external dispfun
      complex*16 dispfun
c----------------
      check=.FALSE.
      check=(((k_root.eq.1).or.(k_root.eq.2)).or.(k_root.eq.3))
c      write(*,*)'hamilt_eigenvalue_combined k_root ',k_root,
c      'check ',check

      if (.not.check) then
         write(*,*)'in subroutine hamilt_eigenvalue_combined'
         write(*,*)'kroot=',k_root
         write(*,*)'But it should be k_root=1,2,3'
         stop 'in subroutine hamilt_eigenvalue_combined'
      endif

c-----calculate combined  Eric Nelson Melby-Abhay Ram non-Hermitian tensor
c     reps(3,3)  

cSAP091104
c      write(*,*)'hamilt_eigenvalue_combined'
c      write(*,*)'z,r,phi,cnpar,cnper,cnprim,k_root',
c     &z,r,phi,cnpar,cnper,cnprim,k_root

      X_e=x(z,r,phi,1)
      Y_e=y(z,r,phi,1)    !question for electron -y? 
      T_e=tempe(z,r,phi,1) !kev
      compl_nper=dcmplx(cnper,cnprim)

cSAP091104
c      write(*,*) 'hamilt_eigenvalue_combined X_e,T_e.Y_e',
c     &X_e,T_e,Y_e

      call  Disp_combined(T_e,cnpar,X_e,Y_e,
     & compl_nper,reps,hamilt_combined) ! reps is in Stix system

cc      write(*,*) 'hamilt_eigenvalue_combined compl_nper',compl_nper
ccSAP091216 calculate Im(Nperp)c
c      errabs=1.d-6
c      nsig=6
c      nknown=0
c      nrts=1                 ! just look for one root
c      nguess=nrts 
c      nnew=nrts
c      itmax=50 
c      roots(1)=dcmplx(cnper,0.0d0)
c       write(*,*) 'hamilt_eigenvalue_combine before muller'
c      call muller(dispfun,errabs,nsig,nknown,nguess,nnew,roots,
c     +     itmax,info,ier)
cc      write(*,*) 'hamilt_eigenvalue_combine after muller'
c      cnprim=abs(imag(roots(1)))
c      cnper=abs(dble(roots(1)))  
c      compl_nper=dcmplx(cnper,cnprim)
cc      write(*,*)'hamilt_eigenvalue_comb muller compl_nper',compl_nper

c-----calculate the complex eigenvalues of the dispersion tensor
      cnper_c=cmplx(cnper,cnprim)
      cnpar_c=cmplx(cnpar,0.d0)         

      call lambda(cnper_c,cnpar_c,reps,eigenvalue)

c      write(*,*)'in hamilt_eigenvalue k_root',k_root,
c     &' eigenvalue',eigenvalue
  
c-----calculate Hamiltonian: D=Real(eigenvalue(k_root))
      Dc=eigenvalue(k_root)
      D=dreal(Dc)

cSAP091218
      goto 10

      D=dreal(eigenvalue(1))*dreal(eigenvalue(2))*dreal(eigenvalue(3))

cSAP091218
c      goto 10

c-----choose eigenvalue(k) that gives min=dabs(dreal(eigenvalue(k)))
      min=1.d10
      do k=1,3
       if(dabs(dreal(eigenvalue(k))).lt.min) then
         min=dabs(dreal(eigenvalue(k)))
         Dc=eigenvalue(k)
         D=dreal(Dc)
       endif
      enddo 
 10   continue

      return
      end


      subroutine eigenvalue_hot(z,r,phi,cnpar,cnper,cnprim,
     &eigenvalue)
   
c-----calculate three complex eigenvalues of the despersion tensor 
c     for the hot plasma dielectric tensor
 
    
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'eps.i'
  
c-----input    
      real*8 z,r,phi,! space coordinates.     
     &cnpar,cnper,cnprim!Re(N_parallel),Re(N_perpendicular),
                        !Im(N_perpendicular)      
c     nbulk is the number of plasma components. It was set in one.i
c-----output
      complex*16 eigenvalue(3)

c-----externals
      real*8 x,y,tempe,tpoprho,vflowrho,b
      complex*16 dhot_sum

c-----locals
      integer ihermloc,i
      complex*16 hamiltc
      complex*16 cnper_c,cnpar_c
             
      real*8 x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     & tpop_ar(nbulka),vflow_ar(nbulka),te
      
      complex*16 compl_nper

c---- calculate hot  non-relativistic plasma 
c     complex dielectric tensor reps(3,3)

      bmod=b(z,r,phi) !calculate small radius rho
                      !and put rho into one.i 

      do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpop_zrp(z,r,phi,i) !YuP[2024-08-14] was tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
      enddo

      compl_nper=dcmplx(cnper,cnprim)
      ihermloc=2 !for non-hermitian tensor

c-----calculate the hot plasma non-Hermitian tensor 
      hamiltc=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &   vflow_ar,cnpar,cnper,ihermloc,reps)
          
c-----calculate the complex eigenvalues of the dispersion tensor
      cnper_c=cmplx(cnper,cnprim)
      cnpar_c=cmplx(cnpar,0.d0)         

      call lambda(cnper_c,cnpar_c,reps,eigenvalue)

c     write(*,*)'hot eigenvalue',eigenvalue
  
      return
      end
