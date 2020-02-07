c        ********************** arrays ***********************
c        *                      ------                       *
c        * this subroutine initializes the arrays  deru  and *
c        * prmt that drkgs  (runge-kutta  subroutine)  later *
c        * uses.					     *
c        *****************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      ndim - number of differential equations			   !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c								   !
c        output parameters					   !
c								   !
c      deru - vector of weight coefficients of solution error.	   !
c             later it becomes the vector of right hand side of    !
c	      the geometrical optics equations.			   !
c      prmt - vector of parameters for drkgs (initial step of	   !
c             integration, accuracy of solution etc. see also	   !
c             drkgs.for).					   !
c------------------------------------------------------------------
      subroutine arrays(ndim,deru,prmt,ihlf)
      !implicit double precision (a-h,o-z)
      implicit none !YuP[2020-01-14]
      integer ndim,ihlf !YuP[2020-01-14]
      integer i  !YuP[2020-01-14]
      real*8 deru,prmt !YuP[2020-01-14]
      real*8 ttt  !YuP[2020-01-14]

      include 'param.i'
      include 'one.i'
      dimension deru(*),prmt(*)

      ttt=1.d0/ndim !YuP[2020-01-14] Not good to use 't'. 
      !It is in name_adj.i and adj_nml.i
      !YuP: changed 't' to 'ttt', just in case.
      do 100 i=1,ndim
         deru(i)=ttt
 100  continue
      prmt(7)=prmt(1)+prmt(6)
      prmt(8)=prmt(1)+prmt(3)
 
c      read(1,20)prmt(1),prmt(2),prmt(3)
c      read(1,25)prmt(4),prmt(6)
c      read(1,30)ihlf
c      ihlf=prmt(9)
c      prmt(9)=ihlf

      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'in arrays prmt(1)=',prmt(1)
      write(*,*)'in arrays prmt(2)=',prmt(2)
      write(*,*)'in arrays tau=prmt(3)=',prmt(3)
      write(*,*)'in arrays acccuracy=prmt(4)=',prmt(4)
      write(*,*)'in arrays hprint=prmt(6)=',prmt(6)
      write(*,*)'in arrays hamiltonian acccuracy=prmt(9)=',prmt(9)
c     write(*,*)'in arrays ihlf=',ihlf
      endif ! outprint
      
 20   format(f6.3/f8.3/d7.1)
 25   format(d7.1/d7.1)
 30   format(i3)

      return
      end

