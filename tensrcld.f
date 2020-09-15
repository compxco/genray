c        ********************tensrcld************************
c        * 						      *
c        * this subroutine  calculates the components of      *
c        * the complex dielectric tensor reps(3,3)            *
c        * for cold plasma                                    *
c        *****************************************************
c-------------------------------------------------------------
c       input parameters				     !
c     	rz,r,phi coordinates  of ray point		     !
c	.....						     !
c	rho-small radius was calculated in subroutine: b     !
c            rho is in common one    			     !
c       output parameter				     !
c       reps(3,3)-complex dielectric tensor in common /eps/  !
c------------------------------------------------------------
      subroutine tensrcld(z,r,phi)
      implicit none !integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'eps.i'
      real*8 z,r,phi !INPUT
      real*8 x,y !external functions
      real*8 xe,ye,xi,yi,epsper,epspar,g !local
      integer i
c	 write(*,*)'in tensrcld '
       xe=x(z,r,phi,1)
       ye=y(z,r,phi,1)
       epsper=1.d0-xe/(1.d0-ye*ye)
       g=xe*ye/(1.d0-ye*ye)

       epspar=1.d0-xe

	 do i=2,nbulk
          xi=x(z,r,phi,i)
          yi=y(z,r,phi,i)
          epsper=epsper-xi/(1.d0-yi*yi)
          epspar=epspar-xi
          g=g-xi*yi/(1.d0-yi*yi)
	 enddo
c-------------------------------------------------
	 reps(1,1)=dcmplx(epsper,0.d0)
	 reps(1,2)=dcmplx(0.d0,g)
	 reps(2,1)=-reps(1,2)
	 reps(2,2)=reps(1,1)
	 reps(3,3)=dcmplx(epspar,0.d0)
	 reps(1,3)=dcmplx(0.d0,0.d0)
	 reps(3,1)=dcmplx(0.d0,0.d0)
	 reps(2,3)=dcmplx(0.d0,0.d0)
	 reps(3,2)=dcmplx(0.d0,0.d0)
c-------------------------------------------------
c       for test
c        gel=xe*ye/(1-ye*ye)
c	gion=0.d0
c	do i=2,nbulk
c            xi=x(z,r,phi,i)
c            yi=y(z,r,phi,i)
c	    gion=gion+xi*yi/(1.d0-yi*yi)
c	enddo
c	gtest=gel-gion
c	write(*,*)'in tensrcld gel,gion,gtest,reps(1,2)'
c	write(*,*)gel,gion,gtest,reps(1,2)
c-------------------------------------------------

	 return
	 end

