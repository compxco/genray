c-----------------------------------------------------
c      the data for the spline coefficients for the distribution function
c      f(y,x,rho,ngena)
c
c      rho(ir) is the normalized small radius ir=1,lrz
c      y(i,ir) is the pitch angle i=1,...iy_(ir)
c      x(j) is the normalized momentun x=P/(mVnorm) j=1,jx 
c      ngena is the max number of plasma species 
             
c      double precision fxx,fyy,fxxyy,wk
c      common/spline_distrib/ fxx(iya,jxa,lrza,ngena),
c     &fyy(iya,jxa,lrza,ngena),fxxyy(iya,jxa,lrza,ngena),
c     &ibd(4),wk(nwka)

      real*8, pointer ::
     &fxx(:,:,:,:),             !(iya,jxa,lrza,ngena)  
     &fyy(:,:,:,:),             !(iya,jxa,lrza,ngena)  
     &fxxyy(:,:,:,:),           !(iya,jxa,lrza,ngena)
     &wk(:)                     !(nwka=1+3*max(jxa,iya)),

      integer ibd(4),nwk

      common/spline_distrib/
     &fxx,             !(iya,jxa,lrza,ngena)  
     &fyy,             !(iya,jxa,lrza,ngena)  
     &fxxyy,           !(iya,jxa,lrza,ngena)  
     &wk,              !(nwka),
     &ibd,             !(4),
     &nwk              ! nwka
