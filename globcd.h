!     globcd.h
!------------------------------------------------------------------------
      REAL*8 small
      parameter (small=1.d-3)
      REAL*8 tmass
      parameter (tmass=511.d0) 
!
      REAL*8 pi      !  pi=acos(-1.d0)   set in TorGA_curgap
      common/cnstcd/pi
!
      REAL*8 tolval  ! =max(tol,1.d-6) set  in TorGA_curgap
      common /sftwcd/tolval  
!
      integer modelv   !  modelv=model      set  in TorGA_curgap
      integer igv      !  igv=ig            set  in TorGA_curgap
      common /mdlpcd/modelv,igv
!
      REAL*8 yval, enz, enzsq, sgnnz ! yval=yy,             set  in TorGA_curgap
                                   ! enz=npara,               set  in TorGA_curgap
                                   ! enzsq=enz*enz            set  in TorGA_curgap
                                   ! sgnnz=sign(1.0d0,npara) set  in TorGA_curgap
      integer nharm                ! nharm=lh                 set  in TorGA_curgap
      common /wavccd/yval,enz,enzsq,sgnnz,nharm
!
      REAL*8 hloc, href, hav, hsqav, chrtav 
              ! set in TorGa_ceqmdl 
              ! hloc=(1.d0-eps)/(1.d0+eps*cos(thetap))
              ! href=(1.d0-eps)/(1.d0+eps)
              ! hav=1.d0-eps
              ! hsqav=(1.d0-eps)**2/sqrt(1.d0-eps**2)
              ! chrtav=(sqrt(2.d0*eps*(1.d0-eps))+(1.d0+eps)*   
              ! &        asin(sqrt(2.d0*eps/(1.d0+eps))))/pi
!
      REAL*8 cxi2, cxi4, fc, ft
              ! set in TorGA_curgac 
              ! cxi2,cxi4,
              ! set in TorGA_curgac or in TorGa_getftrap
              ! fc,ft
      common /geoqcd/hloc,href,hav,hsqav,chrtav,cxi2,cxi4,fc,ft
!
      REAL*8 zeff, zrat, tau, etcutoff, gammin, gammax
              ! zeff=zeffin            set in TorGA_curgac
              ! zrat=(zeff+1.d0)/fc   set in TorGA_curgac
              ! tau=tebulk/tmass       set in TorGA_curgac
              ! etcutoff=20.d0        set in TorGA_curgac
              ! gammin                 set in TorGa_getlims
              ! gammax                 set in TorGa_getlims
      common /ztmpcd/zeff,zrat,tau,etcutoff,gammin,gammax
!
      REAL*8  eta    ! set in different places of toray.f90
      common /collcd/eta
