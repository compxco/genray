
c     The parameters used in mathcurba subroutine
c --- rrange.i
c
c --- to specify the range of real numbers
c
      real        rsmall,              rbig
c      parameter (rsmall = 1.797e-99,   rbig = 2.225e99)
      parameter (rsmall = 1.e-38,   rbig = 1.e38)
      real       srsmall,             srbig
c      parameter (srsmall = 3.434e-50, srbig = 4.717e49)
      parameter (srsmall = 1.e-19, srbig = 1.e19)
c
c      srsmall = rsmall**(1/2)     srbig = rbig**(1/2)
c
      real       trsmall,             trbig
c      parameter (trsmall = 1.216e-33, trbig = 1.305e33)
      parameter (trsmall = 2.154433e-13, trbig = 4.641593e12)
c
c      trsmall = rsmall**(1/3)     trbig = rbig**(1/3)
c
      real       ttrsmall,             ttrbig
c      parameter (ttrsmall = 1.478e-66, ttrbig = 1.704e66)
      parameter (ttrsmall = 4.641581e-26, ttrbig = 2.154438e25)

c
c     ttrsmall = rsmall**(2/3)    ttrbig = rbig**(2/3)
c
