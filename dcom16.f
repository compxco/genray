
************************************************************************
*               this subroutine is for the components                  *
*               of the dielectric tensor in a hot plasma.              *
*               (part of Mazzucato dispersion relation).               *
************************************************************************
        subroutine dcom16(np,y,mu,nz)
        implicit double precision (a-h,o-z)
	double precision mu,nz,mmdei
        common /comp/ a11p,a11m,a13p,a13m,a33p,a33m,a130,a330,
     1                b11p,b11m,b13p,b13m,b33p,b33m,b110,b130,b330,
     2                c11p,c11m,c13p,c13m,c33p,c33m,
     3                a11,a13,a33,b11,b13,b33,c11,c13,c33,
     4                d11,d13,d33,f11,f13,f33,g11,g13,g33,ep0
        common /step/ du
************************************************************************
*               calculate hermitian components                         *
************************************************************************
c	""""""""""""""""""""""""""""""
c	write(*,*)'beginning of dcom16 np=',np,'y=',y,'mu=',mu,'nz=',nz
c	"""""""""""""""""""""""""""""""
        a11p=0.0d0
        a11m=0.0d0
        a13p=0.0d0
        a13m=0.0d0
        a33p=0.0d0
        a33m=0.0d0
        a130=0.0d0
        a330=0.0d0
        b11p=0.0d0
        b11m=0.0d0
        b13p=0.0d0
        b13m=0.0d0
        b33p=0.0d0
        b33m=0.0d0
        b110=0.0d0
        b130=0.0d0
        b330=0.0d0
        c11p=0.0d0
        c11m=0.0d0
        c13p=0.0d0
        c13m=0.0d0
        c33p=0.0d0
        c33m=0.0d0
        a11=0.0d0
        a13=0.0d0
        a33=0.0d0
        b11=0.0d0
        b13=0.0d0
        b33=0.0d0
        c11=0.0d0
        c13=0.0d0
        c33=0.0d0
        d11=0.0d0
        d13=0.0d0
        d33=0.0d0
        f11=0.0d0
        f13=0.0d0
        f33=0.0d0
        g11=0.0d0
        g13=0.0d0
        g33=0.0d0
        ep0=0.0d0
        iopt=3
c*               uo from mu*(gamma1-1)=umax                              *
        umax=100.0d0
cc	write(*,*)'mu=',mu,' np=',np
        uo=dsqrt(umax*(umax+2.0d0*mu))/mu
        du=2.0d0*uo/np
        do 1000 i=1,np
        u=-uo+du*i
        gamma1=dsqrt(1.0d0+u*u)
        gex=-mu*(gamma1-1.0d0)
        x0=-mu*(gamma1-nz*u)
        xp=-mu*(gamma1-y-nz*u)
        xm=-mu*(gamma1+y-nz*u)
        x2p=-mu*(gamma1-2.0d0*y-nz*u)
        x2m=-mu*(gamma1+2.0d0*y-nz*u)
        if(x0.eq.0.0d0) go to 100
        ei0=mmdei(iopt,dble(x0),ier)
        go to 200
100     ei0=0.0d0
200     if(xp.eq.0.0d0) go to 300
        eip=mmdei(iopt,dble(xp),ier)
        go to 400
300     eip=0.0d0
400     if(xm.eq.0.0d0) go to 500
        eim=mmdei(iopt,dble(xm),ier)
        go to 600
500     eim=0.0d0
600     if(x2p.eq.0.0d0) go to 700
        ei2p=mmdei(iopt,dble(x2p),ier)
        go to 800
700     ei2p=0.0d0
800     if(x2m.eq.0.0d0) go to 900
        ei2m=mmdei(iopt,dble(x2m),ier)
        go to 950
900     ei2m=0.0d0
950     ep0=ep0+u*u*ei0*dexp(gex)
        del=(1.0d0+(2.0d0*mu*gamma1+xp)*(1.0d0-xp*eip))*dexp(gex)
        a11p=a11p+del
        a13p=a13p+u*del
        a33p=a33p+u*u*del
c       print 567,a11p,a13p,a33p
c 567     format(' a11p'1pd12.2,'  a13p'1pd12.2,'  a33p'1pd12.2)
        del=(1.0d0+(2.0d0*mu*gamma1+xm)*(1.0d0-xm*eim))*dexp(gex)
        a11m=a11m+del
        a13m=a13m+u*del
        a33m=a33m+u*u*del
        del=(1.0d0+(2.0d0*mu*gamma1+x0)*(1.0d0-x0*ei0))*dexp(gex)
        a130=a130+u*del
        a330=a330+u*u*del
        del=(6.0d0+2.0d0*(4.0d0*mu*gamma1+xp)+(4.0d0*mu**2*gamma1**2+
     1      4.d0*mu*xp*gamma1+xp**2)*(1.d0+xp*(1.d0-xp*eip)))*dexp(gex)
        b11p=b11p+del
        b13p=b13p+u*del
        b33p=b33p+u*u*del
        del=(6.0d0+2.0d0*(4.0d0*mu*gamma1+xm)+(4.0d0*mu**2*gamma1**2+
     1      4.d0*mu*xm*gamma1+xm**2)*(1.d0+xm*(1.d0-xm*eim)))*dexp(gex)
        b11m=b11m+del
        b13m=b13m+u*del
        b33m=b33m+u*u*del
        del=(6.0d0+2.0d0*(4.0d0*mu*gamma1+x0)+(4.0d0*mu**2*gamma1**2+
     1      4.d0*mu*x0*gamma1+x0**2)*(1.d0+x0*(1.d0-x0*ei0)))*dexp(gex)
        b110=b110+del
        b130=b130+u*del
        b330=b330+u*u*del
        del=(6.0d0+2.0d0*(4.0d0*mu*gamma1+x2p)+(4.0d0*mu**2*gamma1**2+
     1      4.d0*mu*x2p*gamma1+x2p**2)*(1.d0+x2p*(1.d0-x2p*ei2p)))
     *      *dexp(gex)
        c11p=c11p+del
        c13p=c13p+u*del
        c33p=c33p+u*u*del
        del=(6.0d0+2.0d0*(4.0d0*mu*gamma1+x2m)+(4.0d0*mu**2*gamma1**2+
     1      4.d0*mu*x2m*gamma1+x2m**2)*(1.d0+x2m*(1.d0-x2m*ei2m)))
     *      *dexp(gex)
        c11m=c11m+del
        c13m=c13m+u*del
        c33m=c33m+u*u*del
1000    continue
cc       write(*,*)'in dcom16 after 1000 continue'
************************************************************************
*               calculate anti-hermitian components                    *
************************************************************************
        r2=y**2-1.d0+nz**2
        if(r2.eq.0.d0.or.r2.lt.0.d0) go to 10000
        r=dsqrt(r2)
        vp=(nz*y+r)/(1.d0-nz*nz)
        vm=(nz*y-r)/(1.d0-nz*nz)
        ppp1=1.d0+vp*vp
        gp=dsqrt(1.d0+vp*vp)
        gm=dsqrt(1.d0+vm*vm)
        q=1.d0/nz/mu
*               q is q/csi of fidone
        csi=(1.d0-nz**2)/nz/mu/r
*               csi is 1/csi of fidone
        ep=-mu*(gp-1.d0)
        em=-mu*(gm-1.d0)
        a11=(1.d0+csi)*dexp(ep)+
     1      (1.d0-csi)*dexp(em)
        a13=(vp+vp*csi+2.d0*q+3.d0*q*csi)*dexp(ep)+
     1      (vm-vm*csi+2.d0*q-3.d0*q*csi)*dexp(em)
        a33=(vp**2+vp*(vp*csi+4.d0*q)+6.d0*q*(vp*csi+q)+
     *      12.d0*q**2*csi)*
     1      dexp(ep)+
     2      (vm**2-vm*(vm*csi-4.d0*q)-6.d0*q*(vm*csi-q)-
     *      12.d0*q**2*csi)*
     3      dexp(em)
        b11=(1.d0+3.d0*csi+3.d0*csi**2)*dexp(ep)-
     1      (1.d0-3.d0*csi+3.d0*csi**2)*dexp(em)
        b13=(vp+3.d0*(vp*csi+q)+3.d0*(vp*csi**2+4.d0*q*csi)+
     *      15.d0*q*csi**2)*
     1      dexp(ep)-
     2      (vm-3.d0*(vm*csi-q)+3.d0*(vm*csi**2-4.d0*q*csi)+
     *      15.d0*q*csi**2)*
     3      dexp(em)
        b33=(vp**2+3.d0*vp*(vp*csi+2.d0*q)+3.d0*(vp**2*csi**2+
     *       8.d0*q*vp*csi+
     1      4.d0*q**2)+30.d0*q*csi*(vp*csi+2.d0*q)+
     *      90.d0*q**2*csi**2)*dexp(ep)-
     2      (vm**2-3.d0*vm*(vm*csi-2.d0*q)+3.d0*(vm**2*csi**2-
     *       8.d0*q*vm*csi+
     3      4.d0*q**2)+30.d0*q*csi*(vm*csi-2.d0*q)+
     *      90.d0*q**2*csi**2)*dexp(em)
        d11=(1.d0+6.d0*csi+15.d0*csi**2+15.d0*csi**3)*dexp(ep)+
     1      (1.d0-6.d0*csi+15.d0*csi**2-15.d0*csi**3)*dexp(em)
        d13=(vp*(1.d0+6.d0*csi+15.d0*csi**2+15.d0*csi**3)+
     1      2.d0*q*(2.d0+15.d0*csi+45.d0*csi**2+
     *      105.d0/2.d0*csi**3))*dexp(ep)+
     2      (vm*(1.d0-6.d0*csi+15.d0*csi**2-15.d0*csi**3)+
     3      2.d0*q*(2.d0-15.d0*csi+45.d0*csi**2-105.d0/2.d0*csi**3))*
     *      dexp(em)
        d33=(vp**2*(1.d0+6.d0*csi+15.d0*csi**2+15.d0*csi**3)+
     1      4.d0*q*vp*(2.d0+15.d0*csi+45.d0*csi**2+105.d0/2.d0*csi**3)+
     2      10.d0*q**2*(2.d0+18.d0*csi+63.d0*csi**2+84.d0*csi**3))*
     3      dexp(ep)+
     4      (vm**2*(1.d0-6.d0*csi+15.d0*csi**2-15.d0*csi**3)+
     5      4.d0*q*vm*(2.d0-15.d0*csi+45.d0*csi**2-105.d0/2.d0*csi**3)+
     6      10.d0*q**2*(2.d0-18.d0*csi+63.d0*csi**2-84.d0*csi**3))*
     7      dexp(em)
10000   continue
        r2=(2.d0*y)**2-1.d0+nz**2
        if(r2.eq.0.d0.or.r2.lt.0.d0) go to 10001
        r=dsqrt(r2)
	pppp=1.d0-nz*nz
        vp=(nz*(2.d0*y)+r)/(1.d0-nz*nz)
        vm=(nz*(2.d0*y)-r)/(1.d0-nz*nz)
        gp=dsqrt(1.d0+vp*vp)
        gm=dsqrt(1.d0+vm*vm)
        q=1.d0/nz/mu
*               q is q/csi of fidone
        csi=(1.d0-nz**2)/nz/mu/r
*               csi is 1/csi of fidone
        ep=-mu*(gp-1.d0)
        em=-mu*(gm-1.d0)
        c11=(1.d0+3.d0*csi+3.d0*csi**2)*dexp(ep)-
     1      (1.d0-3.d0*csi+3.d0*csi**2)*dexp(em)
        c13=(vp+3.d0*(vp*csi+q)+3.d0*(vp*csi**2+4.d0*q*csi)+
     *      15.d0*q*csi**2)*
     1      dexp(ep)-
     2      (vm-3.d0*(vm*csi-q)+3.d0*(vm*csi**2-4.d0*q*csi)+
     *      15.d0*q*csi**2)*
     3      dexp(em)
        c33=(vp**2+3.d0*vp*(vp*csi+2.d0*q)+3.d0*(vp**2*csi**2+
     *       8.d0*q*vp*csi+
     1      4.d0*q**2)+30.d0*q*csi*(vp*csi+2.d0*q)+90.d0*q**2*csi**2)*
     *      dexp(ep)-
     2      (vm**2-3.d0*vm*(vm*csi-2.d0*q)+3.d0*(vm**2*csi**2-
     *      8.d0*q*vm*csi+
     3      4.d0*q**2)+30.d0*q*csi*(vm*csi-2.d0*q)+90.d0*q**2*csi**2)
     *      *dexp(em)
        f11=(1.d0+6.d0*csi+15.d0*csi**2+15.d0*csi**3)*dexp(ep)+
     1      (1.d0-6.d0*csi+15.d0*csi**2-15.d0*csi**3)*dexp(em)
        f13=(vp*(1.d0+6.d0*csi+15.d0*csi**2+15.d0*csi**3)+
     1      2.d0*q*(2.d0+15.d0*csi+45.d0*csi**2+105.d0/2.d0*csi**3))*
     *      dexp(ep)+
     2      (vm*(1.d0-6.d0*csi+15.d0*csi**2-15.d0*csi**3)+
     3      2.d0*q*(2.d0-15.d0*csi+45.d0*csi**2-105.d0/2.d0*csi**3))*
     *      dexp(em)
        f33=(vp**2*(1.d0+6.d0*csi+15.d0*csi**2+15.d0*csi**3)+
     1      4.d0*q*vp*(2.d0+15.d0*csi+45.d0*csi**2+105.d0/2.d0*csi**3)+
     2      10.d0*q**2*(2.d0+18.d0*csi+63.d0*csi**2+84.d0*csi**3))*
     3      dexp(ep)+
     4      (vm**2*(1.d0-6.d0*csi+15.d0*csi**2-15.d0*csi**3)+
     5      4.d0*q*vm*(2.d0-15.d0*csi+45.d0*csi**2-105.d0/2.d0*csi**3)+
     6      10.d0*q**2*(2.d0-18.d0*csi+63.d0*csi**2-84.d0*csi**3))*
     7      dexp(em)
10001   continue
        r2=(3.d0*y)**2-1.d0+nz**2
        if(r2.eq.0.d0.or.r2.lt.0.d0) go to 10002
        r=dsqrt(r2)
        vp=(nz*(3.d0*y)+r)/(1.d0-nz*nz)
        vm=(nz*(3.d0*y)-r)/(1.d0-nz*nz)
        gp=dsqrt(1.d0+vp*vp)
        gm=dsqrt(1.d0+vm*vm)
        q=1.d0/nz/mu
*               q is q/csi of fidone
        csi=(1.d0-nz**2)/nz/mu/r
*               csi is 1/csi of fidone
        ep=-mu*(gp-1.d0)
        em=-mu*(gm-1.d0)
        g11=(1.d0+6.d0*csi+15.d0*csi**2+15.d0*csi**3)*dexp(ep)+
     1      (1.d0-6.d0*csi+15.d0*csi**2-15.d0*csi**3)*dexp(em)
        g13=(vp*(1.d0+6.d0*csi+15.d0*csi**2+15.d0*csi**3)+
     1      2.d0*q*(2.d0+15.d0*csi+45.d0*csi**2+105.d0/2.d0*csi**3))*
     *      dexp(ep)+
     2      (vm*(1.d0-6.d0*csi+15.d0*csi**2-15.d0*csi**3)+
     3      2.d0*q*(2.d0-15.d0*csi+45.d0*csi**2-105.d0/2.d0*csi**3))*
     *      dexp(em)
        g33=(vp**2*(1.d0+6.d0*csi+15.d0*csi**2+15.d0*csi**3)+
     1      4.d0*q*vp*(2.d0+15.d0*csi+45.d0*csi**2+105.d0/2.d0*csi**3)+
     2      10.d0*q**2*(2.d0+18.d0*csi+63.d0*csi**2+84.d0*csi**3))*
     3      dexp(ep)+
     4      (vm**2*(1.d0-6.d0*csi+15.d0*csi**2-15.d0*csi**3)+
     5      4.d0*q*vm*(2.d0-15.d0*csi+45.d0*csi**2-105.d0/2.d0*csi**3)+
     6      10.d0*q**2*(2.d0-18.d0*csi+63.d0*csi**2-84.d0*csi**3))*
     7      dexp(em)
10002   continue

c The following were not properly commented out until Nov. '93.
c Originally were in code as part of checking out additional effects
c of terms added by Mazzucato.  (BobH).
c        f11=0.d0
c        f13=0.d0
c        f33=0.d0
c	"""""""""""""""""""""""""""""
c        write(*,*)'end of dcom16'
c        write(*,*)'a11,a13,a33',a11,a13,a33
c        write(*,*)'b11,b13,b33',b11,b13,b33
c        write(*,*)'c11,c13,c33',c11,c13,c33
c        write(*,*)'d11,d13,d33',d11,d13,d33
c        write(*,*)'f11,f13,f33',f11,f13,f33
c        write(*,*)'g11,g13,g33',g11,g13,g33
c	"""""""""""""""""""""""""""""

        return
        end

