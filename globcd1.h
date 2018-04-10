!globcd1.h
!-------------------------------------------------------------------------
!
      REAL*8 enperp, omodev 
              !enperp=nperp          set in TorGA_curgap
              !omodev=omode          set in TorGA_curgap
      COMPLEX*16 cez, ceplus, ceminus
              ! cez=cefldz                            set in TorGA_curgap
              ! ceplus=cefldx+(0.d0, 1.d0)*cefldy  set in TorGA_curgap
              ! ceminus=cefldx-(0.d0, 1.d0)*cefldy set in TorGA_curgap
      COMMON /wavccd1/cez,ceplus,ceminus,omodev, enperp
