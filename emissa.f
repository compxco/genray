c
c
c
      subroutine emissa(nv,nt,psi,bes,besp,nbes,bmax,nint,omega,omegpe,
     & omegce,nsigma,nwarm,npar,nper,rmedmi,zeff,cdvt,n,rri,alpha,emis,
     & denn,den,ird,temp,vmax1,vpar1,vpar2,vper1,polm,polp,
     & polx,poly,polz,crt)
c
c     emissa calculates the absorption and emission coefficients for
c     ece using the prescription of Bornatici
c
c     Form tensor elements for epsilon-a and G.
c     Then multiply by normalized electric field elements, and
c     divide by normalized energy density * group velocity.
c     Output is alpha and emis (the j-coefficient) in cgs.
c
      implicit double precision (a-h,o-z)
      
      parameter(ninta=201)
      double complex rootm1,exde,eyde,ezde,cedsde,crt,
     1     	cexde,ceyde,cezde,
     &        ccexde,cceyde,ccezde
      dimension bes(nbes),besp(nbes)
      dimension xint(ninta),bes2(ninta),besp2(ninta),cs(ninta),
     & sn(ninta),vper(ninta),vpar(ninta),gamma(ninta),b(ninta),
     & ff(ninta),u(ninta),edsde(ninta),denom(ninta),den(ird)
cSm
      double precision s11p(ninta),s12p(ninta),s13p(ninta),
     & s21p(ninta),s22p(ninta),s23p(ninta),
     & s31p(ninta),s32p(ninta),s33p(ninta)
 
      double complex epsa11,epsa12,epsa13,
     & epsa21,epsa22,epsa23,
     & epsa31,epsa32,epsa33

      
       
      double precision kpar,kper,npar,nper,me
c
      common /pi/ pi
c
      if(nint.gt.ninta)  stop 6
c
      rootm1=(0.d0,1.d0)
      root2=dsqrt(2.d0)
      c=2.9979d10
      me=9.1095d-28
      cdvt2=cdvt**2
      alpha=0.d0
      emis=0.d0
      crt=(0.d0,0.d0)
c
c     Find range of integration
      call intinit(omega,n,npar,cdvt,vmax1,v0,vpar0,vpar1,
     & vpar2,thet1,thet2,ires)
      write(*,*)'emissa vmax1,vpar0,vpar1', vmax1,vpar0,vpar1
      write(*,*)'emissa npar,nper,1/omega,(omegpe/omega)**2',
     &npar,nper,1.d0/omega,(omegpe/omega)**2
      
c
      if(ires.gt.1)  return
c
c     Obtain e-field polarizations and get correct value of nper
c     for this npar
c


c      call grpde2(npar,nper,omega,omegpe,nsigma,
c     &            vgrpdc,edenfac,exde,eyde,ezde)

c zero out ion effects in cold dielectric tensor
      icold=1

c      call grpden(npar,nper,omega,omegpe,1.d0,zeff,
c    & vgrpdc0,eplusde,eminde,edenfac0,exder,eyder,ezder,icold)

c
      polm=cdabs(ccexde-rootm1*cceyde)/root2
      polp=cdabs(ccexde+rootm1*cceyde)/root2
      polx=cdabs(ccexde)
      poly=cdabs(cceyde)
      polz=cdabs(ccezde)
c
c
c     normalized kpar,kper
      kpar=npar*omega/cdvt
      kper=nper*omega/cdvt
      dbes=bmax/(nbes-1)
      nbesm1=nbes-1
c     Subdivide theta range of integration
      dth=(thet2-thet1)/(nint-1)
      vper1=dsqrt(1.d0-npar**2)*v0
c     Prepare for integration
      do 100  i=1,nint
      xint(i)=thet1+(i-1)*dth
      cs(i)=dcos(xint(i))
      sn(i)=dsin(xint(i))
      vpar(i)=vpar0-v0*cs(i)
      vper(i)=vper1*sn(i)
      b(i)=kper*vper(i)
      ii=1.0+b(i)/dbes
      if(ii.le.nbesm1)  go to 110
      bes2(i)=0.d0
      besp2(i)=0.d0
      stop 11
      go to 120
 110  bes2(i)=bes(ii)+(b(i)-(ii-1)*dbes)*
     & (bes(ii+1)-bes(ii))/(dbes)
      besp2(i)=besp(ii)+(b(i)-(ii-1)*dbes)*
     & (besp(ii+1)-besp(ii))/(dbes)
c
c bes2 = bessel fn, besp2 = deriv of bessel fn. (don't square them)
c
c     bes2(i)=bes2(i)**2
c     besp2(i)=besp2(i)**2
c
 120  continue
      v2=(vpar(i)**2+vper(i)**2)
      gamma(i)=dsqrt(1.d0+v2/cdvt2)
      v=dsqrt(v2)
c     Find equatorial plane pitch angle for calc of distn fctn:
      rootpsi=dsqrt(psi)
      vpeout=vper(i)/rootpsi
      vpaout=dsign(1.d0,vpar(i))*dsqrt(v2-vpeout*vpeout)
      t0=datan2(vpeout,vpaout)
      call distn(nv,nt,v,t0,ff(i),dfdx,dfdy,denn,den,ird,cdvt2)
c      Transform derivatives to local poloidal angle:
      dfdy=vpar(i)/vpaout/rootpsi*dfdy
      ct=vpar(i)/v
      st=vper(i)/v
  130 continue
      dfdvpar=ct*dfdx-st*dfdy/v
      dfdvper=st*dfdx+ct*dfdy/v
      u(i)=n/dabs(kpar)*dfdvper+dsign(1.d0,kpar)*vper(i)*dfdvpar
c
c     E.S.E/E**2
c
      s11=n*n*bes2(i)*bes2(i)/(kper*kper)
      s12=-n*bes2(i)*besp2(i)*vper(i)/kper
cS      s12=n*bes2(i)*besp2(i)*vper(i)/kper

      s21=-s12
      s22=vper(i)*vper(i)*besp2(i)*besp2(i)
      s13=n*vpar(i)*bes2(i)*bes2(i)/kper
      s31=s13
      s23=vpar(i)*vper(i)*bes2(i)*besp2(i)
cS      s23=-vpar(i)*vper(i)*bes2(i)*besp2(i)

      s32=-s23
      s33=vpar(i)*vpar(i)*bes2(i)*bes2(i)
c
      cedsde= dconjg(ccexde)*ccexde*s11
     1      + dconjg(ccexde)*cceyde*rootm1*s12
     2      + dconjg(cceyde)*ccexde*rootm1*s21
     3      + dconjg(cceyde)*cceyde*s22
     4      + dconjg(ccexde)*ccezde*s13
     5      + dconjg(ccezde)*ccexde*s31
     6      + dconjg(cceyde)*ccezde*rootm1*s23
     7      + dconjg(ccezde)*cceyde*rootm1*s32
     8      + dconjg(ccezde)*ccezde*s33
      edsde(i)=dreal(cedsde)
c
c     denom(i)=gamma(i)*abs(gamma(i)*kpar-vpar(i)/cdvt2*omega)
      denom(i)=dabs(gamma(i)*kpar-vpar(i)/cdvt2*omega)
cSm
      s11p(i)=s11
      s12p(i)=s12
      s13p(i)=s13
      s21p(i)=s21
      s22p(i)=s22
      s23p(i)=s23
      s31p(i)=s31
      s32p(i)=s32
      s33p(i)=s33
      
      write(*,*)'i s11p',i,s11p(i),s12p(i),s13p(i),
     &s21p(i),s22p(i),s23p(i),
     &s31p(i),s32p(i),s33p(i)
      write(*,*)'edsde(i)',edsde(i)

 100  continue
c
c
c     Integrating
      e1=0.d0
      g1=0.d0

cSm
      epsa11=dcmplx(0.d0,0.d0)
      epsa12=0.d0
      epsa13=0.d0
      epsa22=0.d0
      epsa23=0.d0
      epsa33=0.d0


c     endpoints
c     do 200  i=1,nint,nint
      do 200  i=1,nint,nint-1
      dvper=0.5d0*dabs(vper(2)-vper(1))
      if(i.ne.1)  dvper=0.5d0*dabs(vper(nint)-vper(nint-1))
      e1=e1+dvper*u(i)*edsde(i)/denom(i)
c     g1=g1+gamma(i)*dvper*vper(i)*edsde(i)*ff(i)/denom(i)
      g1=g1+dvper*vper(i)*edsde(i)*ff(i)/denom(i)

cSm
      epsa11=epsa11+s11p(i)*dvper*u(i)/denom(i)
      epsa12=epsa12+s12p(i)*dvper*u(i)/denom(i)
      epsa13=epsa13+s13p(i)*dvper*u(i)/denom(i)
      epsa22=epsa22+s22p(i)*dvper*u(i)/denom(i)
      epsa23=epsa23+s23p(i)*dvper*u(i)/denom(i)
      epsa33=epsa33+s33p(i)*dvper*u(i)/denom(i)
 200  continue
c
      do 201  i=2,nint-1
      dvper=0.5d0*dabs(vper(i+1)-vper(i-1))
      e1=e1+dvper*u(i)*edsde(i)/denom(i)
c     g1=g1+gamma(i)*dvper*vper(i)*edsde(i)*ff(i)/denom(i)
      g1=g1+dvper*vper(i)*edsde(i)*ff(i)/denom(i)
cSm
      epsa11=epsa11+s11p(i)*dvper*u(i)/denom(i)
      epsa12=epsa12+s12p(i)*dvper*u(i)/denom(i)
      epsa13=epsa13+s13p(i)*dvper*u(i)/denom(i)
      epsa22=epsa22+s22p(i)*dvper*u(i)/denom(i)
      epsa23=epsa23+s23p(i)*dvper*u(i)/denom(i)
      epsa33=epsa33+s33p(i)*dvper*u(i)/denom(i)
 201  continue
c
      e1=-2.d0*pi**2*(omegpe/omega)**2*dabs(kpar)*e1/denn
      g1=1.d0/(2.d0*pi)**3*0.5d0*me*(c/cdvt)**2*omegpe**2*omegce*g1/denn
c

cSm
      epsa11=-2.d0*pi**2*(omegpe/omega)**2*dabs(kpar)/denn*epsa11
      epsa12=-2.d0*pi**2*(omegpe/omega)**2*dabs(kpar)/denn*epsa12
      epsa13=-2.d0*pi**2*(omegpe/omega)**2*dabs(kpar)/denn*epsa13
      epsa22=-2.d0*pi**2*(omegpe/omega)**2*dabs(kpar)/denn*epsa22
      epsa23=-2.d0*pi**2*(omegpe/omega)**2*dabs(kpar)/denn*epsa23
      epsa33=-2.d0*pi**2*(omegpe/omega)**2*dabs(kpar)/denn*epsa33     
      epsa21=-epsa12
cSm030515
c      epsa13=epsa13
      epsa31=epsa13
      epsa32=-epsa23

      write(*,*)'emissa alpha 00 e1',e1
      write(*,*)'emissa omega*omegce,vgrpdc*c*edenfac',
     +omega*omegce,vgrpdc*c*edenfac
      WRITE(*,*)'(omega*omegce)/(4.d0*pi)/(vgrpdc*c*edenfac)',
     +(omega*omegce)/(4.d0*pi)/(vgrpdc*c*edenfac)
      write(*,*)'vgrpdc,c,edenfac',vgrpdc,c,edenfac
      write(*,*)'epsa',epsa11,epsa12,epsa13
      write(*,*)epsa21,epsa22,epsa23
      write(*,*)epsa31,epsa32,epsa33

      
      alpha=(omega*omegce)/(4.d0*pi)*e1/(vgrpdc*c*edenfac)
      write(*,*)'alpha',alpha
      emis=pi*rri**2*(omega*omegce/c)**2*g1/(vgrpdc*c*edenfac)
      

c
c multiply by 4*pi since in cgs Poynting vector has a 4*pi in it
c
      alpha=alpha*4.d0*pi
      write(*,*)'emiisa alpha*4*pi',alpha
      emis=emis*4.d0*pi
c     omega=omegaold
c
      return
      end
