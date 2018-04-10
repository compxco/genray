c************************************zfunc***********************
c     calculates complex dispersion plasma function  and derivatives
c     for complex argument. 
c     This routine is from Brambilla LH/FW code
c-----------------------------------------------------------------
      subroutine zfunc(w,z,zp)
      complex z,zp,w
c     calculates plasma dispersion function z = (zr,zi) and its
c     derivative zprime = (zpr,zpi) for the point u = (x,y)
      data sqrtpi/1.7724538509055/
      x=real(w)
      y=aimag(w)
c      write(*,*)'in zfunc before pfunc x,y',x,y
      call pfunc(x,y,zr,zi)
c      write(*,*)'in zfunc after pfunc zr,zi',zr,zi
      xx=zr
      zr=-sqrtpi*zi
      zi=sqrtpi*xx
      zpr=-2*(1.+x*zr-y*zi)
      zpi=-2*(y*zr+x*zi)
      z=cmplx(zr,zi)
      zp=cmplx(zpr,zpi)
      return
      end
      subroutine pfunc(x,y,u,v)
c
c     pfunc calculates z/(sqrt(pie)*i) = (u,v) at point (x,y)
      save et2
      dimension w287(4),w283(4),et2(63)
      data pie/3.14159265358979/,sqrtpi/1.7724538509055/
      data osqpi/.56418958354776/
c     8 point hermite gaussian quadrature  (+,-) x(i)
      data w283/.381186990207322,1.157193712446780,1.981656756695843,
     1          2.930637420257244/
c     8 point hermite gaussian quadrature w(i)
      data w287/.6611470125582,.2078023258149,1.707798300741e-2,
     1          1.996040722114e-4/
      data itest/0/
      if(itest.gt.0) go to 5
c     set up array of e(-t*t)
c     values of exp(-t*t) for values of .08 le x le 5.04 by steps of .08
      do 3 i=1,63
      t=float(i)*0.08
      et2(i)=exp(-t*t)
    3 continue
      itest=1
c     find which quadrant point u = (x,y) is in
    5 ii=1
      assign 244 to j
      c5=x
      c6=y
      if(c5.lt.0.0)go to 8
      if(c6.lt.0.0)go to 287
      go to 11
    8 if(c6.ge.0.0)go to 14
      assign 245 to i
      go to 20
   11 assign 257 to i
      go to 46
   14 assign 255 to i
      go to 46
c     point has y lt 0, calculate analytic continuation part of
c     expression = 2*exp(-u*u)/sqrt(pie)
   20 z=c6*c6-c5*c5
      z=amin1(z,80.)
      z=amax1(z,-80.)
      co=exp(z)
      c7=co+co
      co=c5*c6
      c9=co+co
      c8=-c7*sin(c9)
      c7=c7*cos(c9)
   46 c5=abs(c5)
      c6=abs(c6)
c     test for magnitude of y to pick step size to do
c     integral 1./(pie*i)*exp(-t*t)/(t-u)
      if(c5.ge.6.0)go to 219
   50 if(c6.le.0.5)go to 65
      if(c6.le.3.0)go to 61
      if(c6.gt.6.0)go to 219
      ns=6
      go to 73
   61 if(c6.le.1.5)go to 71
      ns=3
      go to 73
c     will do integral for y = 0.5 even if y is actually lt 0.5
c     will use taylor series to get final answer.
   65 c10=c6
      c6=0.5
      assign 128 to j
   71 ns=1
c     evaluate integrand at zero
   73 c62=c6*c6
      c19=c5*c5
      c20=c62+c19
      c19=1./c20
      c17=c19*c6
      c18=c19*c5
      c9=0.08
c     do trapezoidal integration from -5 lt x lt 5
      do 100 m=ns,63,ns
      ws=float(m)*c9
      c1m=c5-ws
      c1p=c5+ws
      c20p=et2(m)/(c62+c1p*c1p)
      c20m=et2(m)/(c62+c1m*c1m)
      c17=c6*(c20p+c20m)+c17
      c18=c1m*c20m+c1p*c20p+c18
  100 continue
      ws=c9/pie*float(ns)
      c17=c17*ws
      c18=c18*ws
      go to j,(128,244)
c     do taylor series expansion for case y lt 0.5 to get correct answer
  128 c11=c17
      c12=c18
      c9=2.0
      c6=c10-0.5
      c6=c6+c6
      c10=c11/2.0
      c13=(c5*c12+c10-osqpi)*c6
      c10=c12/2.0
      c14=(-c5*c11+c10)*c6
      c17=c11+c13
      c18=c12+c14
  165 c10=c6/c9
      c19=c13/2.0
      c19=c5*c14+c19
      c15=(c6/2.0*c11+c19)*c10
      c17=c15+c17
      t1=c5*c13
      c19=(c6*c12+c14)/2.0
      c16=(-t1+c19)*c10
      c18=c16+c18
      t1=c17+c15
      if(abs(t1-c17).gt.1.e-7*abs(t1+c17)) go to 207
      t1=c18+c16
      if(abs(t1-c18).le.1.e-7*abs(t1+c18)) go to 244
  207 c11=c13
      c12=c14
      c13=c15
      c14=c16
      c9=c9+1.0
      go to 165
c     do 8 point hermite gaussian quadrature for values of x or y gt 6
  219 c17=0.0
      c18=0.0
      c62=c6*c6
      do 230 m=1,4
      c1m=c5-w283(m)
      c1p=c5+w283(m)
      c2m=w287(m)/(c62+c1m*c1m)
      c2p=w287(m)/(c62+c1p*c1p)
      c17=c6*(c2m+c2p)+c17
      c18=c1m*c2m+c1p*c2p+c18
  230 continue
      c17=c17/pie
      c18=c18/pie
c     reset value of z for correct quadrant
  244 go to i,(245,249,255,257)
  245 c8=-c8
      c18=-c18
  249 c17=c7-c17
      c18=c8-c18
  255 c18=-c18
  257 u=c17
      v=c18
c      write(*,*)'in pfunc before return u,v',u,v
      return
  287 c5=-c5
      assign 249 to i
      go to 20
      end
 






