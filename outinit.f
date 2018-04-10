c        ********************** outinit***********************
c        *                      -----                       *
c        *   outinit is used in dinit to prepare output	    *
c        *   data for 3D  F-P code       		    *
c        *   outinit subroutine.                            *
c        *                                                  *
c        ****************************************************
      subroutine outinit(u)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'write.i'
      dimension u(*)

      z1=u(1)
      r1=u(2)
      phi1=u(3)
      cnz1=u(4)
      cnr1=u(5)
      cm1=u(6)

      nrayelt=0
      zold=u(1)
      rold=u(2)

      return
      end

