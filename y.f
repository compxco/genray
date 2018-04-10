c        **********************   y   ***********************
c        *                        -                         *
c        * this function calculates the ratio of  cyclotron *
c        * frequency of species i and wave frequency        *
c        * y=omega_cyclotron_i/omega                        *
c        ****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of the  point  where  the  ratio  is  !
c                 calculated.      		         	    !
c      small radius rho is in common one
c-------------------------------------------------------------------
      double precision
     1function y(z,r,phi,i)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'one.i'
      y=w(i)*bmod
c      write(*,*)'in y i,w(i),bmod,rho,y',i,y,w(i),bmod,rho,y
      if(dabs(y).lt.1.d-8) y=1.d-8
      return
      end

      double precision
     1function y_test(z,r,phi,i)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'one.i'
      y=w(i)*bmod

c      write(*,*)'in y_test z,r,phi,ib,bmod,w(i),y',
c     &z,r,phi,ib,bmod,w(i),y

c      write(*,*)'in y i,w(i),bmod,rho,y',y,w(i),bmod,rho,y
      if(dabs(y).lt.1.d-8) y=1.d-8

c      write(*,*) 'in y_test y',y
cSAP090122
      y_test=y

      return
      end
