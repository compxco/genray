c
cDefines some REAL*4 constants used when calling pgplot library routines
c
      REAL*4 :: R40=0.,R4P1=.1,R4P2=.2,R4P3=.3,R4P4=.4,R4P5=.5
      REAL*4 :: R4P8=.8,R4P95=.95,R41P5=1.5
      REAL*4 :: R41=1.,R42=2.,R43=3.,R44=4.,R45=5.,R46=6.,R47=7.,R48=8.
      REAL*4 :: R43P8=3.8,R431=31.,R470=70.,R4110=110.,R41P4=1.4
      REAL*4 :: R445=45.,R4165=165.,R44P75=4.75,R410=10.

!BH201124: For some reason, use of the following leads to a
!  "multiple definition of `pgconst_'" compiler error, when
!  both contour.f and mk_graph.f are loaded together.
!     Even number of cnsts for common, to maintain alignment
!      common /pgconst/R40,R4P1,R4P2,R4P3,R4P4,R4P5,R4P8,R4P95,R41P5,
!     +                R41,R42,R43,R44,R45,R46,R47,R48,
!     +                R43P8,R431,R470,R4110,R41P4,
!     +                R445,R4165,R44P75,R410 

