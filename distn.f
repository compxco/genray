c
c
c
      subroutine distn(nv,nt,v0,thet0,ff,dfdx,dfdy,denn,den,ird,cdvt2)
      implicit double precision (a-h,o-z)
c      include 'clich1.h'
      dimension den(ird)
c     Could speed up by eliminating repeated searches of grid
c     in terp2, for x and y.
c
c interpolate (using density) to get distribution functions
c  Location alogrithm for interpolation points assumes that density
c  is a decreasing function of radius.
c
      do 1 i=1,ird
      ir=i
      if(denn.gt.den(ir)) goto 2
1     continue
2     idm=nva
c
c
      ff=(denn/den(ir))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,ir),fxx(1,1,ir),fyy(1,1,ir),fxxyy(1,1,ir),idm,0,0,1)
      dfdx=(denn/den(ir))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,ir),fxx(1,1,ir),fyy(1,1,ir),fxxyy(1,1,ir),idm,1,0,0)
      dfdy=(denn/den(ir))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,ir),fxx(1,1,ir),fyy(1,1,ir),fxxyy(1,1,ir),idm,0,1,0)
c
      if(ir.eq.1.or.(ir.eq.ird.and.denn.le.den(ird))) return
      im=ir-1
      ffm=(denn/den(im))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,im),fxx(1,1,im),fyy(1,1,im),fxxyy(1,1,im),idm,0,0,0)
      dfdxm=(denn/den(im))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,im),fxx(1,1,im),fyy(1,1,im),fxxyy(1,1,im),idm,1,0,0)
      dfdym=(denn/den(im))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,im),fxx(1,1,im),fyy(1,1,im),fxxyy(1,1,im),idm,0,1,0)
c
      ddm=denn-den(im)
      ddp=den(ir)-denn
      dd=ddm+ddp
      ff=( ffm*ddp + ff*ddm ) / dd
      dfdx=( dfdxm*ddp + dfdx*ddm ) / dd
      dfdy=( dfdym*ddp + dfdy*ddm ) / dd
c
      return
      end


c!

      subroutine distr(nv,nt,v0,thet0,ff,dfdx,dfdy,rho,rovera,ird,cdvt2)
      implicit double precision (a-h,o-z)
      include 'clich1.h'
      dimension den(ird)
c     Could speed up by eliminating repeated searches of grid
c     in terp2, for x and y.
c
c interpolate (using density) to get distribution functions
c  Location alogrithm for interpolation points assumes that density
c  is a decreasing function of radius.
c
      do 1 i=1,ird
      ir=i
      if(denn.gt.den(ir)) goto 2
1     continue
2     idm=nva
c
c
      ff=(denn/den(ir))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,ir),fxx(1,1,ir),fyy(1,1,ir),fxxyy(1,1,ir),idm,0,0,1)
      dfdx=(denn/den(ir))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,ir),fxx(1,1,ir),fyy(1,1,ir),fxxyy(1,1,ir),idm,1,0,0)
      dfdy=(denn/den(ir))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,ir),fxx(1,1,ir),fyy(1,1,ir),fxxyy(1,1,ir),idm,0,1,0)
c
      if(ir.eq.1.or.(ir.eq.ird.and.denn.le.den(ird))) return
      im=ir-1
      ffm=(denn/den(im))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,im),fxx(1,1,im),fyy(1,1,im),fxxyy(1,1,im),idm,0,0,0)
      dfdxm=(denn/den(im))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,im),fxx(1,1,im),fyy(1,1,im),fxxyy(1,1,im),idm,1,0,0)
      dfdym=(denn/den(im))*terp2p(v0,thet0,nv,vs,nt,thets,
     1      f(1,1,im),fxx(1,1,im),fyy(1,1,im),fxxyy(1,1,im),idm,0,1,0)
c
      ddm=denn-den(im)
      ddp=den(ir)-denn
      dd=ddm+ddp
      ff=( ffm*ddp + ff*ddm ) / dd
      dfdx=( dfdxm*ddp + dfdx*ddm ) / dd
      dfdy=( dfdym*ddp + dfdy*ddm ) / dd
c
      return
      end


