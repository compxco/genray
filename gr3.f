c*****************************************************
c----------------subroutine GR3---------------------*
c Creates circles on the horizontal plane (X,Y) with*
c radii RMIN,RMAX,RA (common block 'five')       *
c Output: arrays XT and YT into common block gr     *
c         (file gr.cb)                              *
c*****************************************************


      subroutine GR3

      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'five.i'
      include 'three.i'
      include 'one.i'
c      double precision RMIN,RMAX,RA,PI
      double precision RA

      double precision PHI
      INTEGER I
      include 'gr.i'

      pi=4.d0*datan(1.d0)
      RA=xma

      DO 10 I=0,NP
       PHI=2.0D0*PI*I/DBLE(NP)
       XT(1,I+1)=RMAX*DCOS(PHI)*100.d0*r0x
       YT(1,I+1)=RMAX*DSIN(PHI)*100.d0*r0x

       XT(2,I+1)=RMIN*DCOS(PHI)*100.d0*r0x
       YT(2,I+1)=RMIN*DSIN(PHI)*100.d0*r0x

       XT(3,I+1)=RA*DCOS(PHI)*100.d0*r0x
       YT(3,I+1)=RA*DSIN(PHI)*100.d0*r0x

10    CONTINUE

      END




