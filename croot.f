      subroutine croot(cz0,czf)
      implicit double complex(c)
      implicit double precision(a-b,d-h,j-z)
      parameter(ilw=40,iliw=2,im=2,in=2)
      dimension xc(in),iw(iliw),w(ilw)

cSm010111
c     in is a dimension of the output vector xc(in)
      integer LIV,LV
	parameter (LIV=59)
	parameter (LV=77+in*(in+17)/2)
	
	integer IV(LIV), UIPARM(1)
	double precision D(in), URPARM(1)
	double precision  V(LIV) !is an estimation of the relative error
	
	EXTERNAL CALCF, UFPARM
c     CALCF calcultes the function.It should be equal the sum of the 
c     of the squares of the real and imaginary parts of the dispersion determinant,
c     as a function of real and imaginary k.
cSm_end

c
c     see nag p.2 of e04fde for details
c
      xc(1)=dreal(cz0)
      ifail=1
      xc(2)=dimag(cz0)
c     The following routine is from Naglib.  It finds the minimum of
c     m functions of n variables.  For the case here, it minimizes
c     the sum of the squares of the real and imaginary parts of the
c     dispersion determinant, as a function of real and imaginary k.
cSm001226
c     call e04fde(im,in,xc,err,iw,iliw,w,ilw,ifail)
c	call e04fdf(im,in,xc,err,iw,iliw,w,ilw,ifail)

cSm010111
      small=1.d-7
      if (dabs(xc(1)).gt.small) then
          D(1)=1/xc(1)
	else
	    D(1)=1/small
      endif
      
	if (dabs(xc(2)).gt.small) then
          D(2)=1/xc(2)
	else
	    D(2)=1/small
      endif
      write(*,*)'croot LIV',LIV,'IV',IV  
      call DMNF(in,D,xc,CALCF, IV, LIV, LV, V, UIPARM, URPARM, UFPARM)       

C  ***  MINIMIZE GENERAL UNCONSTRAINED OBJECTIVE FUNCTION USING
C  ***  FINITE-DIFFERENCE GRADIENTS AND SECANT HESSIAN APPROXIMATIONS.
      czf=dcmplx(xc(1),xc(2))
      return
      end

cSm010112
	subroutine CALCF (N, X, NF, F, UIPARM, URPARM, UFPARM)
c     calulates the funtion that should be minimized F=(ReD)^2+(ImD)^2
c     I used the text from the subroutine lsfun1,that calcultes ReD and ReD.   
      implicit none
c-----input
      integer N !N=2
	double precision X(N)
	double precision UIPARM(*), URPARM(*)
      integer NF
      external UFPARM

c-----ouput
	double precision F

c-----locals
      double precision pf

	double complex ckr1f,ckr2f,ckr3f
	integer imode
      common/wdata/ckr1f,ckr2f,ckr3f,imode
      double complex cz,cd
      
	cz=dcmplx(X(1),X(2))
c      call cwarm(cz,cd)
      call  hot_disp(cz,cd)


c**      if(imode.eq.1) go to 100
c**c
c**c      the following procedure is designed to prevent finding a root
c**c      that has been previously obtained, by dividing the factor
c**c      associated with the root out of the dispersion determinant.
c**c     It also avoids root jumping (most of the time).
c**      cd=cd/(cz-ckr1f)
c**      if(imode.eq.2) go to 101
c**      cd=cd/(cz-ckr2f)
c**      if(imode.eq.3) go to 102
c**      cd=cd/(cz-ckr3f)
c**c
c**  102   continue
c**  101   continue
c**  100    continue
	      
	F=dreal(cd)**2+dimag(cd)**2

	return
	end    

	subroutine UFPARM
c      dummy subroutine
	return
	end

       subroutine grad_min(cz0,czf)
c------find the 2D point czf giving the minimal value to F=(ReD)**2+(ImD)**2
       implicit none
c------input
       double complex cz0 !initial value of complex N_perp
c------output
       double complex czf !root of the dispersion functoion=complex N_perp
c------local
       double complex cd !complex dispertion function
       double precision alpha ! uteration parameter
       double precision step,f_old,f_new,df_drealz, df_dimagz
       double precision re_z,imag_z,re_zp,imag_zp,re_zm,imag_zm,fp,fm
       double complex cz_p,cz_m,cdp,cdm,cd_old,cd_new 
       double precision eps !accuracy      
       integer i,iter_max

       iter_max=100
       re_z=dreal(cz0)
       imag_z=dimag(cz0)
c       write(*,*)'in croot grad_min cz0',cz0
       call  hot_disp(cz0,cd_old)
       f_old=dreal(cd_old)**2+dimag(cd_old)**2
c       f_old=dsqrt(f_old)
c       write(*,*)'in croot grad_min begin  f_old',f_old
       step=1.d-6
       eps=1.d-7
       alpha=1.d-12
       i=1

 10    continue

       re_zp=re_z+step
       imag_zp=imag_z+step
       re_zm=re_z-step
       imag_zm=imag_z-step   
    
       cz_p=dcmplx(re_zp,imag_z)
       call  hot_disp(cz_p,cdp)
       fp=dreal(cdp)**2+dimag(cdp)**2
c       fp=dsqrt(fp)
       cz_m=dcmplx(re_zm,imag_z)
       call  hot_disp(cz_m,cdm)
       fm=dreal(cdm)**2+dimag(cdm)**2
c       fm=dsqrt(fm)
       df_drealz=(fp-fm)/(2.d0*step)
c       write(*,*)'fp,fm,step,df_drealz',fp,fm,step,df_drealz

       cz_p=dcmplx(re_z,imag_zp)
       call  hot_disp(cz_p,cdp)
       fp=dreal(cdp)**2+dimag(cdp)**2
c       fp=dsqrt(fp)
       cz_m=dcmplx(re_z,imag_zm)
       call  hot_disp(cz_m,cdm)
       fm=dreal(cdm)**2+dimag(cdm)**2
c       fm=dsqrt(fm)
       df_dimagz=(fp-fm)/(2.d0*step)

c       write(*,*)'fp,fm,df_dimagz',fp,fm,df_dimagz

 30    continue
       czf=cz0-alpha*dcmplx(df_drealz,df_dimagz)
c       write(*,*)'cz0,alpha,czf',cz0,alpha,czf
       call  hot_disp(czf,cd_new)
       f_new=dreal(cd_new)**2+dimag(cd_new)**2
       f_new=dsqrt(f_new)
       i=i+1
c       write(*,*)'i,f_old,f_new',i,f_old,f_new
       if (i.gt.iter_max) goto 20
       
       if(f_new.lt.eps) goto 20
       if(f_new.gt.f_old) then
         alpha=0.5*alpha
c         write(*,*)'alpha',alpha
         goto 30
       else
         f_old=f_new
         cz0=czf 
         goto 10
       endif

 20    continue
       return
       end







