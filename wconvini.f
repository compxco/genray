      subroutine owconvr (theta,x0,rhoconv,zconv,rconv)
c     It detemines the point (zconv,rconv) 
c     and rhoconv where xe=x0,For x0=1 it gives point of O_X conversion. 
c     This point is on the surface rhoconvr and has the
c     given poloidal angle theta (degree) 
c     if theta=0  point is on the outer part of the magnetic surface
c     if theta=90 point is at the top of the magnetic surface
c     input data:
c          xe(rho=1) < x0 < xe(rho=0)  a given value xe=(omega_pe/omega)**2
c          theta is a poloidal angle (degree)
c     output:
c          the point (zconv,rconv) and the small radius rhoconv where xe=xe0
c----------------------------------------------------------------------
c     It solves the equation x(zconv(rhoconv),rconv(rhoconv),phix,1)-xe0=0
      implicit none
      double precision theta,x0 
      double precision phix,racc 
      double precision rhoconv,zconv,rconv
      double precision thetax,xe0,rholeft,rhoright,psiconv,xpp 
      double precision rtbis,psi_rho,xe_eq 
      common /convers/ thetax,phix,xe0
      external xe_eq

      thetax=theta*datan(1.d0)*4.d0/180.d0 ! polidal angle in radian
      xe0=x0
      racc = 1.d-7
c      write(*,*)'wconvini thetax, xe0',thetax,xe0
      rholeft=0.d0
      xpp=xe_eq(rholeft)
c      write(*,*)'xe_eq(0)',xpp
      rhoright=1.d0
      xpp=xe_eq(rhoright)
c      write(*,*)'xe_eq(1)',xpp
c      write(*,*)'wconvini owconvr before rtbis xe0,thetax',xe0,thetax
c      write(*,*)'owconvr bef rtbis rholeft,rhoright',
c     *rholeft,rhoright
  
      rhoconv=rtbis(xe_eq,rholeft,rhoright,racc)
c      write(*,*)'wconvini rhoconv',rhoconv
      psiconv=psi_rho(rhoconv)
c      write(*,*)'wconvini psiconv,thetax',psiconv,thetax
      call zr_psith(psiconv,thetax,zconv,rconv)
c      write(*,*)'wconvin rhoconv,zconv,rconv',rhoconv,zconv,rconv
      return
      end

      double precision FUNCTION rtbis(FUNC,X1,X2,XACC)
c     bisection method of the solution of the equation func(rtbis)=0 
      implicit none
      integer JMAX,J
c      PARAMETER (JMAX=40)
      parameter(JMAX=100)      
      double precision FUNC,FMID,F,X1,X2,XACC,XMID,DX
c      write(*,*)'wconvini rtbis X1,X2,XACC',X1,X2,XACC  
      FMID=FUNC(X2)
c      write(*,*)'wconvini X2,FMID',X2,FMID  
      F=FUNC(X1)
c      write(*,*)'wconvini X1, F',X1,F  
      IF(F*FMID.GE.0.d0) PAUSE 'Root must be bracketed for bisection.'
      IF(F.LT.0.d0)THEN
        rtbis=X1
        DX=X2-X1
      ELSE
        rtbis=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5d0
        XMID=rtbis+DX
        FMID=FUNC(XMID)
c        write(*,*)'wconvini J,DX,XMID,FMID',J,DX,XMID,FMID  
        IF(FMID.LE.0.d0)rtbis=XMID
c        IF(dabs(DX).LT.XACC .OR. FMID.EQ.0.d0) RETURN
         IF(dabs(DX).LT.XACC .OR. dabs(FMID).LT.XACC) RETURN
11    CONTINUE 
      PAUSE 'too many bisections, increase parameter JMAX in rtbis'
      END

      double precision function	xe_eq(rhox)
c-----The function xe_eq=x(z(rhox),r(rhox),phix,1)-xe0
c     input:
c     rhox  -small radius
c     thetax poloidal angle 
c     phix   toroidal angle
c     xe0     the given xe value
c     output:
c     xe_eq=x(z(rhox),r(rhox),phix,1)-xe0
c--------------------------------------
      implicit none
      double precision rhox
      double precision thetax,phix,xe0
      common /convers/thetax,phix,xe0
      double precision psix,zx,rx,xp,bmode
      double precision psi_rho,x,b
c      write(*,*)'wconvini xe_eq rhox,thetax',rhox,thetax
      psix=psi_rho(rhox)
c      write(*,*)'wconvini xe_eq psix',psix
      call zr_psith(psix,thetax,zx,rx)
c      write(*,*)'wcon xe_eq rhox psix thetax',rhox,psix,thetax
c      write(*,*)'wcon xe_eq zx,rx',zx,rx
      bmode=b(zx,rx,phix)
c      write(*,*)'wcon zx,rx,phix,bmode',zx,rx,phix,bmode
      xp=x(zx,rx,phix,1)
c      write(*,*)'wcon xp',xp
      xe_eq=xp-xe0
c      write(*,*)'wconv xe_eq xe0,xp,xe_eq',xe0,xp,xe_eq

      return
      end
