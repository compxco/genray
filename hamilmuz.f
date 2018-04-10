c        ********************hamiltmuz************************
c        * 						      *
c        * this subroutine  calculates: hamiltmz -real part of*
c        *                              Hamiltonian, and      *
c        * reps(3,3) -components of complex dielectric tensor *
c        *            for Mazzucato code                      *
c        *****************************************************
c-------------------------------------------------------------
c       input parameters				     !
c       cnpar - parallel (to magnetic field) component	     !
c               of the refractive index			     !
c       cnper - perpendicular component of the refractive    !
c               index                    		     !
c	ihermloc=1 Hermitian dielectric tensor       	     !
c               =2 full      dielectric tensor               !
c       Attention: iherm is in common one.i		     !
c     	rz,r,phi coordinates  of ray point		     !
c	.....						     !
c	rho-small radius was calculated in subroutine: b     !
c            rho is in common one    			     !
c       mode = +1., O-mode   is in common one		     !
c              -1., X-mode				     !
c       output parameter				     !
c       hamiltmz-Hamiltonian				     !
c       reps(3,3)-complex dielectric tensor in common /eps/  !
c------------------------------------------------------------
        subroutine hamiltmuz(cnpar,ihermloc,z,r,phi,cnper,hamiltmz)
	implicit double precision (a-h,o-z)
	double precision mode,mu,nz,np2
	double precision ncld(2)
        include 'param.i'
        include 'one.i'
	dimension an2(2),icuto(2)
	double complex sol,cpp
	double complex chamilmz
	 mode=dfloat(ioxm)
c	 write(*,*)'in hamilmuz mode=',mode
         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
	 te=tempe(z,r,phi,1)
	 mu=512.44d00/te
c	 write(*,*)'in hamiltmuz xe,ye,te',xe,ye,te
         if(dabs(cnpar).lt.1.d-4) then
           if (cnpar.ge.0.d0) then
             sign=1.d0
           else 
             sign=-1.d0
           endif
           cnpar=sign*1.d-4
         endif

	 nz=cnpar
cSm060225
         call n_par_min_maz(cnpar,nz)

c	 write(*,*)'in hamilmuz cnpar,cnper,ihermloc',
c     *  cnpar,cnper,ihermloc

	 cpp=dcmplx(cnper,0.0d0)
         np2=cnper*cnper
         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
	 sol=cpp*cpp

	 call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)
c	 write(*,*)'hamiltmuz chamilmz',chamilmz
	 hamiltmz=dreal(chamilmz)
c	 write(*,*)'in hamiltmuz hamiltmz=',hamiltmz
	 return
	 end
c---------------------------------------------------------------------
*****************************complx1***********************************
*               this subroutine calculates: chamilmz-the complex       *
*               hamiltonian  of the hot dispersion relation;	       *
*               sol-complex N_perp.;                                   *
*               reps(3,3)-complex dielectric tensor (in common/eps/)   *
c----------------------------------------------------------------------
*               It uses dcom161.sub and dten161.sub.                   *
************************************************************************
        subroutine complx1(mode,x,y,mu,nz,np2,sol,iherm,chamilmz)
        implicit double precision (a-h,o-z)
        double precision mode,mu,nz,np2
c
c
c  input:
c     mode = +1., O-mode
c            -1., X-mode
c     x    = (omegape/omega)**2
c     y    = omegace/omega
c     mu   = c**2/(Te/m)
c     nz   = parallel refractive index
c     np2  = estimate of perp ref. index (typically from coldm)
c     sol  = dcomplx (cnperp,0)**2
c     iherm=1 Hermitian dielectric tensor
c          =2 full      dielectric tensor
c  output:
c     chamilmz
c     The routine is described by E.Mazzucato, I.Fidone, and G.Granata,
c     Phys. Fl., vol. 30, p. 3745-3751 (1987), and references therein.
c     Published results obtained with the aid of this module should
c     contain a reference to is source.
c     (I obtained this code from Ernesto Mazzucato in 1987, and modified
c      it slightly to put it into the present subroutine form.
c      Bob Harvey,1990).
c
c
c
      include 'param.i'
      include 'eps.i'
        double complex eps1(3,3),eps2(3,3),eps3(3,3),eps(3,3),ceps(3,3)
        double complex c1,c2,c3,so1,so2,solu,det,sss,sol,cpp2,cpp4
	double complex chamilmz
csmirno961217test
c      double complex  dlamb(3,3),ff,dd,p,af,ad,p2
c	""""""""""""""""""""""""""""""""""
c	write(*,*)'beg.  complx1'
c	""""""""""""""""""""""""""""""""""

        call dcom161(1000,y,mu,nz,iherm)
c        call dcom161(5000,y,mu,nz,iherm)
        call dten161(x,y,mu,nz,eps0,eps1,eps2,eps3,iherm)
        solu=sol
c	write(*,*)'in complx1 sol=',sol
c	write(*,*)'in complx1 solu=',solu
        do 2000 i=1,3
        do 2000 j=1,3
        eps(i,j)=eps1(i,j)+solu*eps2(i,j)+solu**2*eps3(i,j)
2000    continue
        c1=eps(1,1)+2.d0*nz*eps(1,3)+eps(1,3)**2+nz**2*eps(3,3)-
     1          eps(1,1)*eps(3,3)
c**	write(*,*)'c1=',c1
        c2=eps(1,1)*(nz**2-eps(2,2))-eps(1,2)**2+eps0*(nz**2-eps(1,1))+
     1     2.d0*nz**3*eps(1,3)-2.d0*nz*eps(2,2)*eps(1,3)+
     2     2.d0*nz*eps(1,2)*
     2     eps(2,3)+eps(1,1)*eps(2,3)**2-eps(2,2)*eps(1,3)**2+
     2     2.d0*eps(1,2)*eps(1,3)*eps(2,3)+
     3     nz**2*(eps(1,3)**2-eps(2,3)**2)+nz**2*eps(3,3)*(nz**2-
     4     eps(1,1)-eps(2,2))+eps(3,3)*(eps(1,2)**2+eps(1,1)*eps(2,2))
c**	write(*,*)'c2=',c2
        c3=eps0*(nz**2*(nz**2-eps(1,1)-eps(2,2))+eps(1,2)**2+
     1     eps(1,1)*eps(2,2))
c**	write(*,*)'c3=',c3
        sol=cdsqrt(solu)
	cpp2=solu
	cpp4=solu*solu
        chamilmz=cpp4*c1+cpp2*c2+c3
c        write(*,*)'hamilmuz solu',solu
c        write(*,*)'hamilmuz cpp4',cpp4
c        write(*,*)'hamilmuz c1',c1
c        write(*,*)'hamilmuz c2',c2
c        write(*,*)'hamilmuz c3',c3
c        write(*,*)'hamilmuz',chamilmz
c------------------------------------------------
      ceps(1,1)=eps(1,1)
      ceps(1,2)=eps(1,2)
      ceps(1,3)=sol*eps(1,3)
      ceps(2,2)=eps(2,2)
      ceps(2,3)=sol*eps(2,3)
      ceps(3,3)=eps0+solu*eps(3,3)
      ceps(2,1)=-ceps(1,2)
      ceps(3,1)=ceps(1,3)
      ceps(3,2)=-ceps(2,3)
c      write(*,*)'in hamilmuz chamiltmz',chamilmz
c      write(*,*)'in hamilmuz tes chamiltmz',
c     1(nz**2-ceps(1,1))*((nz**2+cpp2-ceps(2,2))*(cpp2-ceps(3,3))
c     1 -(-ceps(2,3))*(-ceps(3,2)))-
c     2(-ceps(1,2))*((-ceps(2,1))*(cpp2-ceps(3,3))-
c     2(-ceps(2,3))*(-nz*sol-ceps(3,1)))+
c     3(-nz*sol-ceps(1,3))*((-ceps(2,1))*(-ceps(3,2))-
c     3(nz*nz+cpp2-ceps(2,2))*(-nz*sol-ceps(3,1)))
c      write(*,*)'eps(1,1)=',ceps(1,1)
c      write(*,*)'eps(1,3)=',ceps(1,3)
c      write(*,*)'eps(3,3)=',ceps(3,3)
c      write(*,*)'eps(2,2)=',ceps(2,2)
c      write(*,*)'eps(1,2)=',ceps(1,2)
c      write(*,*)'eps(2,3)=',ceps(2,3)
c-----------------------------------------------------------
c	""""""""""""""""""""""""""""""""""
c	write(*,*)'end of complx1'
c	""""""""""""""""""""""""""""""""""
       do 1 i=1,3
         do 1 j=1,3
 1	  reps(i,j)=ceps(i,j)
cSmirnov961217 test beg
c      dlamb(1,1)=nz*nz-reps(1,1)
c      dlamb(1,2)=-reps(1,2)
c      dlamb(1,3)=-nz*sol-reps(1,3)
c      dlamb(2,1)=-reps(2,1)
c      dlamb(2,2)=nz*nz+cpp2-reps(2,2)
c      dlamb(2,3)=-reps(2,3)
c      dlamb(3,1)=-nz*sol-reps(3,1)
c      dlamb(3,2)=-reps(3,2)
c      dlamb(3,3)=cpp2-reps(3,3)
c      ff=dlamb(1,1)*dlamb(2,2)+dlamb(1,1)*dlamb(3,3)+
c     1 dlamb(2,2)*dlamb(3,3)
c     1 -dlamb(1,2)*dlamb(2,1)-dlamb(1,3)*dlamb(3,1)-
c     1 dlamb(2,3)*dlamb(3,2)
c      ss=dlamb(1,1)+dlamb(2,2)+dlamb(3,3)
c      dd=dlamb(1,1)*(dlamb(2,2)*dlamb(3,3)-dlamb(2,3)*dlamb(3,2))
c     1  -dlamb(1,2)*(dlamb(2,1)*dlamb(3,3)-dlamb(3,1)*dlamb(2,3))
c     2  +dlamb(1,3)*(dlamb(2,1)*dlamb(3,2)-dlamb(2,2)*dlamb(3,1))
c3     format(2d15.3)
c      write(*,*)'in complx1 ff,dd,af,ad'
c      write(*,3)dreal(ff),dimag(ff)
c      write(*,*)'in control ss'
c      write(*,3)dreal(ss),dimag(ss)
c      write(*,3)dreal(dd),dimag(dd)
c      p=dcmplx(1.d0,0.d0)/(ss*ss-2.d0*ff)
c      af=ff*p
c      p2=dd**(2.d0/3.d0)
c      ad=3*p2*p
c      write(*,*)'in control p,p2,af,ad'
c      write(*,3)dreal(p),dimag(p)
c      write(*,3)dreal(p2),dimag(p2)
c      write(*,3)dreal(af),dimag(af)
c      write(*,*)'in complx1 chamilmz'
c      write(*,3)dreal(chamilmz),dimag(chamilmz)
cSmirnov961217 test end
        return	   
        end
c----------------------------------------------------------------------
*               this subroutine is for the components                  *
*               of the dielectric tensor in a hot plasma.              *
************************************************************************
        subroutine dcom161(np,y,mu,nz,iherm)
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
c	write(*,*)'beginning of dcom161 np=',np,'y=',y,'mu=',mu,'nz=',nz
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
       if (iherm.eq.1) goto 161
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
161     continue
c	"""""""""""""""""""""""""""""
c	write(*,*)'end of dcom161'
c        write(*,*)'a11,a13,a33',a11,a13,a33
c        write(*,*)'b11,b13,b33',b11,b13,b33
c        write(*,*)'c11,c13,c33',c11,c13,c33
c        write(*,*)'d11,d13,d33',d11,d13,d33
c        write(*,*)'f11,f13,f33',f11,f13,f33
c        write(*,*)'g11,g13,g33',g11,g13,g33
c	"""""""""""""""""""""""""""""
        return
        end
c------------------------------------------------------------------

************************************************************************
*               this subroutine is for calculating                     *
*               the dielectric tensor in a hot plasma.                 *
************************************************************************
        subroutine dten161(x,y,mu,nz,eps0,eps1,eps2,eps3,iherm)
c	:::::::::::::Harvey:::::::::::::::::::::::::::::::::::;
c        implicit real (a-h,o-z)
c        implicit double precision (a-h,o-z)
c       real nz,mu
c        complex eps1(3,3),eps2(3,3),eps3(3,3)
c        real bes(1),mmu,order
c	::::::::::::::::::::::::::::::::::::::::::::::::

c	:::::::::::::Larisa:::::::::::::::::::::::::::::::::::;
        implicit double precision (a-h,o-z)
        double precision nz,mu,mmu
        double complex eps1(3,3),eps2(3,3),eps3(3,3)
	double precision bes(1),order
c	::::::::::::::::::::::::::::::::::::::::::::::::


        common /comp/ a11p,a11m,a13p,a13m,a33p,a33m,a130,a330,
     1                b11p,b11m,b13p,b13m,b33p,b33m,b110,b130,b330,
     2                c11p,c11m,c13p,c13m,c33p,c33m,
     3                a11,a13,a33,b11,b13,b33,c11,c13,c33,
     4                d11,d13,d33,f11,f13,f33,g11,g13,g33,ep0
        common /step/ du
        common /print/ it
************************************************************************
*               calculate hermitian components                         *
************************************************************************
c	""""""""""""""""""""""""""
c        write(*,*)'beg.of dten161'
c	""""""""""""""""""""""""""

        it=0
      pi=datan2(0.d0,-1.d0)
        order=2.d0
        n=1
        mmu=mu

c     mmbskr(arg,order,n,bk,ier)  is a imsl subroutine, computing
c     the modified Bessel function K_(order+n)(arg), i=0,n-1,
c     scaled by exp(arg),
c     for real positive arguments and nonnegative real fractional
c     order order.  The results are given by bk(i),i=1,n.
cc        call mmbskr(mmu,order,n,bes,ier)
cc        bk=bes(1)
c------------------------------------------------------------
      if (mmu.lt.50) then
        call besk(mmu,2,bk,ier)

c     Nov. 12, 1993.   Bob H.
c     The following line pointed out as necessary by Gary Smith.
c     (an error must have crept in in conversion from mmbskr, somewhere...)
      bk=bk*dexp(mmu)
c--------------------------------------------------------------
c     the nfollowing lines were made by A.Smirnov
      else
      call besk2as(mmu,bk)
      end if
c-------------------------------------------------------------
c       write(*,*)'bk=',bk,' y=',y
c       z=16.0
c       bk=sqrt(pi/(2.0*mu))*(1.0+(z-1.0)/(8.*mu)+
c     1          (z-1.)*(z-9.)/(128.*z**2)+(z-1.)*(z-9.)*(z-25.)/3072./
c     2          mu**3)
        eps0=1.d0+mu*mu*x*ep0/2.d0/bk*du
c      write(*,*)'eps0=',eps0
        c=-x/8.d0/bk*du
        aa11p=a11p*c
        aa11m=a11m*c
        c=c/y
        aa13p=a13p*c
        aa13m=a13m*c
        aa130=a130*c
        c=c/y
        aa33p=a33p*c
        aa33m=a33m*c
        aa330=a330*c
        c=x/y**2/32.d0/mu**2/bk*du
        bb11p=b11p*c
        bb11m=b11m*c
        bb110=b110*c
        c=c/y
        bb13p=b13p*c
        bb13m=b13m*c
        bb130=b130*c
        c=c/y
        bb33p=b33p*c
        bb33m=b33m*c
        bb330=b330*c
        c=-x/y**2/32.d0/mu**2/bk*du
        cc11p=c11p*c
        cc11m=c11m*c
        c=c/(2.d0*y)
        cc13p=c13p*c
        cc13m=c13m*c
        c=c/(2.d0*y)
        cc33p=c33p*c
        cc33m=c33m*c
c--------------------------------------------------------------------
c the following	 lines were made by A.Smirnov
	if (iherm.eq.1) then
	aa11=0.0d0
	aa13=0.0d0
	aa33=0.0d0
	bb11=0.0d0
	bb13=0.0d0
	bb33=0.0d0
	dd11=0.0d0
	dd13=0.0d0
	dd33=0.0d0
	cc11=0.0d0
	cc13=0.0d0
	cc33=0.0d0
	ff11=0.0d0
	ff13=0.0d0
	ff33=0.0d0
	gg11=0.0d0
	gg13=0.0d0
	gg33=0.0d0
	goto 161
	end if
************************************************************************
*               calculate anti-hermitian components                    *
************************************************************************
        r2=y*y-1.d0+nz*nz
        if(r2.lt.0.d0) r2=0.d0
        r=dsqrt(r2)
        a=pi*x/4.d0/nz**2/bk
        a0=a*r
        aa11=a11*a0
c       print 665,aa11
c 665     format(' aa11',1pd10.2)
        aa13=a13*a0/y
        aa33=a33*a0/y**2
        b0=a/y**2*r**2/nz/mu
        bb11=b11*b0
c       print 666,x,eps0,aa11,bb11
c 666     format('x',1pd10.2,'eps0',1pd10.2,'a11',1pd10.2,'b11',1pd10.2)
        bb13=b13*b0/y
        bb33=b33*b0/y**2
        d0=5.d0/8.d0/y**4*a*r**3/(nz*mu)**2
        dd11=d11*d0
        dd13=d13*d0/y
        dd33=d33*d0/y**2
10000   continue
        r2=(2.d0*y)**2-1.d0+nz*nz
        if(r2.lt.0.d0) r2=0.d0
        r=dsqrt(r2)
        c0=-a/y**2*r**2/nz/mu
        cc11=c11*c0
        cc13=c13*c0/(2.d0*y)
        cc33=c33*c0/(2.d0*y)**2
        f0=-a/y**4*r**3/(nz*mu)**2
        ff11=f11*f0
        ff13=f13*f0/(2.d0*y)
        ff33=f33*f0/(2.d0*y)**2
        r2=(3.d0*y)**2-1.d0+nz*nz
        if(r2.lt.0.d0) r2=0.d0
        r=dsqrt(r2)
        g0=3.d0/8.d0/y**4*a*r**3/(nz*mu)**2
        gg11=g11*g0
        gg13=g13*g0/(3.d0*y)
        gg33=g33*g0/(3.d0*y)**2
10001   continue
161     continue
************************************************************************
*               calculate the dielectric tensor                        *
************************************************************************
c        write(*,*)'hamilmuz.for dten161'
c        write(*,*)'aa11,aa13,aa33',aa11,aa13,aa33
c        write(*,*)'aa11,aa13,aa33',aa11,aa13,aa33
c        write(*,*)'bb11,bb13,bb33',bb11,bb13,bb33
c        write(*,*)'cc11,cc13,cc33',cc11,cc13,c33
c        write(*,*)'dd11,dd13,dd33',dd11,dd13,dd33
c        write(*,*)'ff11,ff13,ff33',ff11,ff13,ff33
c        write(*,*)'gg11,gg13,gg33',gg11,gg13,gg33
        a=1.d0+aa11p+aa11m
        b=aa11
        if(it.gt.0) print 24999
24999   format(' -------------------------------------')
        if(it.gt.0) print 25000,a,b
25000   format(' herm1=',1pd12.2,' anth1=',1pd12.2)
        eps1(1,1)=dcmplx(a,b)
        a=aa11
        b=-aa11p+aa11m
        if(it.gt.0) print 25000,a,b
        eps1(1,2)=dcmplx(a,b)
        a=aa13p-aa13m
        b=aa13
        if(it.gt.0) print 25000,a,b
        eps1(1,3)=dcmplx(a,b)
        a=1.d0+aa11p+aa11m
        b=aa11
        if(it.gt.0) print 25000,a,b
        eps1(2,2)=dcmplx(a,b)
        a=-aa13
        b=aa13p+aa13m-2.d0*aa130
        if(it.gt.0) print 25000,a,b
        eps1(2,3)=dcmplx(a,b)
        a=aa33p+aa33m-2.d0*aa330
        b=aa33
        if(it.gt.0) print 25000,a,b
        eps1(3,3)=dcmplx(a,b)
        if(it.gt.0) print 24999
        a=bb11p+cc11p+bb11m+cc11m
        b=bb11+cc11
        if(it.gt.0) print 25000,a,b
c 25001   format(' herm2=',1pd12.2,' anth2=',1pd12.2)
        eps2(1,1)=dcmplx(a,b)
        a=2.d0*bb11+cc11
        b=-2.d0*bb11p-cc11p+2.d0*bb11m+cc11m
        if(it.gt.0) print 25000,a,b
        eps2(1,2)=dcmplx(a,b)
        a=bb13p+cc13p-bb13m-cc13m
        b=bb13+cc13
        if(it.gt.0) print 25000,a,b
        eps2(1,3)=dcmplx(a,b)
        a=3.d0*bb11p+cc11p+3.d0*bb11m+cc11m-4.d0*bb110
        b=3.d0*bb11+cc11
        if(it.gt.0) print 25000,a,b
        eps2(2,2)=dcmplx(a,b)
        a=-2.d0*bb13-cc13
        b=2.d0*bb13p+cc13p+2.d0*bb13m+cc13m-3.d0*bb130
        if(it.gt.0) print 25000,a,b
        eps2(2,3)=dcmplx(a,b)
        a=bb33p+cc33p+bb33m+cc33m-1.5d0*bb330
        b=bb33+cc33
        if(it.gt.0) print 25000,a,b
        eps2(3,3)=dcmplx(a,b)
        if(it.gt.0) print 24999
        a=0.0d0
        b=dd11+ff11+gg11
        eps3(1,1)=dcmplx(a,b)
        a=3.d0*dd11+1.5d0*ff11+gg11
        b=0.0d0
        eps3(1,2)=dcmplx(a,b)
        a=0.0d0
        b=dd13+ff13+gg13
        eps3(1,3)=dcmplx(a,b)
        a=0.0d0
        b=37.d0/5.d0*dd11+2.d0*ff11+gg11
        eps3(2,2)=dcmplx(a,b)
        a=-(3.d0*dd13+1.5d0*ff13+gg13)
        b=0.0d0
        eps3(2,3)=dcmplx(a,b)
        a=0.0d0
        b=dd33+ff33+gg33
        eps3(3,3)=dcmplx(a,b)
cc fortest beg
c	write(*,*)'hamilmuz.for dten161'
c	do i=1,3
c	  do j=1,3
c	  write(*,*)'i,j,eps1(i,j),eps2(i,j),eps3(i,j)',
c     1	             i,j,eps1(i,j),eps2(i,j),eps3(i,j)
c	  enddo
c	enddo
cc fortest end

c	""""""""""""""""""""""""""
c	write(*,*)'end of dten161'
c	""""""""""""""""""""""""""

        return
        end

c        ********************hamiltmc************************
c        * 						      *
c        * this subroutine  calculates: hamiltmz -real part of*
c        * Hamiltonian and                                    *
c        * reps(3,3) -components of complex dielectric tensor *
c        * for Mazzucato code                                 *
c        *****************************************************
c-------------------------------------------------------------
c       input parameters				     !
c       cnpar - paralel (to magnetic field) component	     !
c               of the refractive index			     !
c       cnper - perpendicular component of the refractive    !
c               index                    		     !
c	iherm=2 Full full dielectric tenzor		     !
c     	rz,r,phi coordinates  of ray point		     !
c	.....						     !
c	rho-small radius was calculated in subroutine: b     !
c            rho is in common one    			     !
c       mode = +1., O-mode   is in common one		     !
c              -1., X-mode				     !
c       output parameter				     !
c       hamiltmp-real(Hamiltonian)			     !
c       hamiltmq-image(Hamiltonian)     		     !
c       reps(3,3)-complex dielectric tensor in common /eps/  !
c------------------------------------------------------------
      subroutine hamiltmc(cnz,cnr,cm,z,r,phi,hamiltmp,hamiltmq)
      implicit double precision (a-h,o-z)
      double precision mode,mu,nz,np2
      double precision ncld(2)
      include 'param.i'
      include'one.i'	   
      double complex sol,cpp
      double complex chamilmz
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      te=tempe(z,r,phi,1)
      mu=512.44d00/te
      cnpar=(cnz*bz+cnr*br+cm/r*bphi)/bmod
      cn2=cnz*cnz+cnr*cnr+cm*cm/(r*r)
      cnper=dsqrt(cn2-cnpar*cnpar)
      nz=cnpar
      cpp=dcmplx(cnper,0.0d0)
      np2=cnper*cnper
      sol=cpp*cpp
      call complx1(mode,xe,ye,mu,nz,np2,sol,2,chamilmz)
      hamiltmp=dreal(chamilmz)
      hamiltmq=dimag(chamilmz)
      return
      end


      double precision function fomega_maz(cnz,cnr,cm,z,r,phi)
c-----calculate the dispersion function and complex dielectric tensor reps
c     (it will be in eps.i) for the full (hermitian + anty-hermitian)
c     Mazzucato dielectric tensor using Berstein-Frieldnd formula
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'	   
      double precision mu,nz,np2
      double complex sol,cpp
      double complex chamilmz
      double precision wp(nbulka),vp(nbulka)
      step=0.0000001d0
      hw=step*1.0d0
      hfrqnc=hw
      hfrqnc=hfrqnc*frqncy
      bmod=b(z,r,phi)
      do 11 i=1,nbulk
        vp(i)=v(i)
        wp(i)=w(i)
 11   continue
      frqncpl=frqncy+hfrqnc
      df=frqncy/frqncpl
      do 12 i=1,nbulk
        v(i)=vp(i)*df* df
        w(i)=wp(i)*df
 12   continue
c************************************************
      cnrplus=cnr*df
      cnzplus=cnz*df
      cmplus=cm*df
      call hamiltmc(cnzplus,cnrplus,cmplus,z,r,phi,hamiltmp,hamiltmq)
      fplus=hamiltmp*hamiltmp+hamiltmq*hamiltmq
c*************************************************
c----------------------------------------------------------
      frqncmn=frqncy-hfrqnc
      df=frqncy/frqncmn
      do 15 i=1,nbulk
        v(i)=vp(i)*df*df
        w(i)=wp(i)*df
 15   continue
c************************************************
      cnrminus=cnr*df
      cnzminus=cnz*df
      cmminus=cm*df
      call hamiltmc(cnzminus,cnrminus,cmminus,z,r,phi,hamiltmp,hamiltmq)
      fminus=hamiltmp*hamiltmp+hamiltmq*hamiltmq
c*************************************************
      fomega_maz=(fplus-fminus)/(2.0d0*hw)
c-----------------------------------------------------------
      do 14 i=1,nbulk
        v(i)=vp(i)
        w(i)=wp(i)
 14   continue

cSm030426
c     calculate dielectric tensor reps. It will be in eps.i. 
      call hamiltmc(cnz,cnr,cm,z,r,phi,hamiltmp,hamiltmq)

      return
      end

c this subroutine for control of dispersion relation  
c input :double complex reps(3,3), double ccnz=npar, cnx=nper
c output:double complex
c f=eps11*eps22+eps11*eps33+eps22*eps33
c   -eps12*eps21-epps13*eps31-eps23*eps32
c s=eps11+eps22+eps33
c d=determinant=eps11*(eps22*eps33-eps23*eps32)
c              -eps12*(eps21*eps33-eps31*eps23)
c              +eps13*(eps21*eps32-eps22*eps31)
c af=f/(s*s-2.d0*f)
c ad=3*d**(2.d0/3.d0)/(s*s-2.d0*f)
      subroutine control(cnpar,cnper)
      implicit double precision (a-h,o-z)
      include 'eps.i'
      double precision cnpar,cnper
      double complex dlamb(3,3)
      double complex ff,ss,dd,af,ad,p
      double complex cnx,ccnz,p2
c      write(*,*)'in control cnpar,cnper,reps(2,2)',cnpar,cnper,reps(2,2)
c      p2=reps(2,2)
c      write(*,*)'in control p2=reps(2,2)',p2
      ccnz=dcmplx(cnpar,0.d0)
      cnx=dcmplx(cnper,0.d0)
c      write(*,*)'in control ccnz',ccnz
c      write(*,*)'in control cnx',cnx
      dlamb(1,1)=ccnz*ccnz-reps(1,1)
      dlamb(1,2)=-reps(1,2)
      dlamb(1,3)=-ccnz*cnx-reps(1,3)
      dlamb(2,1)=-reps(2,1)
      dlamb(2,2)=ccnz*ccnz+cnx*cnx-reps(2,2)
      dlamb(2,3)=-reps(2,3)
      dlamb(3,1)=-ccnz*cnx-reps(3,1)
      dlamb(3,2)=-reps(3,2)
      dlamb(3,3)=cnx*cnx-reps(3,3)
c------------------------------------------------
2     format(2i2,2d15.3)
c      do i=1,3
c       do j=1,3
c	 write(*,*)'in control i,j,reps(i,j)'
c	 write(*,2)i,j,dreal(reps(i,j)),dimag(reps(i,j))
c       enddo
c      enddo
c      do i=1,3
c       do j=1,3
c	 write(*,*)'in control i,j,dlamb(i,j)'
c	 write(*,2)i,j,dreal(dlamb(i,j)),dimag(dlamb(i,j))
c       enddo
c      enddo
      ff=dlamb(1,1)*dlamb(2,2)+dlamb(1,1)*dlamb(3,3)+
     1 dlamb(2,2)*dlamb(3,3)
     1 -dlamb(1,2)*dlamb(2,1)-dlamb(1,3)*dlamb(3,1)-
     1 dlamb(2,3)*dlamb(3,2)
      ss=dlamb(1,1)+dlamb(2,2)+dlamb(3,3)
      dd=dlamb(1,1)*(dlamb(2,2)*dlamb(3,3)-dlamb(2,3)*dlamb(3,2))
     1  -dlamb(1,2)*(dlamb(2,1)*dlamb(3,3)-dlamb(3,1)*dlamb(2,3))
     2  +dlamb(1,3)*(dlamb(2,1)*dlamb(3,2)-dlamb(2,2)*dlamb(3,1))
3     format(2d15.3)
      write(*,*)'in control ff,dd'
      write(*,3)dreal(ff),dimag(ff)
c      write(*,*)'in control ss'
c      write(*,3)dreal(ss),dimag(ss)
c      write(*,*)'in control dd'
      write(*,3)dreal(dd),dimag(dd)
      p=dcmplx(1.d0,0.d0)/(ss*ss-2.d0*ff)
      af=ff*p
      p2=dd**(2.d0/3.d0)
      ad=3*p2*p
c      write(*,*)'in control p,p2,af,ad'
c      write(*,3)dreal(p),dimag(p)
c      write(*,3)dreal(p2),dimag(p2)
      write(*,*)'in control af,ad'
      write(*,3)dreal(af),dimag(af)
      write(*,3)dreal(ad),dimag(ad)
      return
      end 

c        ********************dispmuz*************************
c        * this function  calculates: real part of            *
c        * Re(Mazzucato Hamiltonian) and as output            *
c        * parameter hamilti Im(Mazzucato Hamiltonian)        *
c        * reps(3,3) -components of complex dielectric tensor *
c        * for Mazzucato code in eps.i                        *
c        *****************************************************
c-------------------------------------------------------------
c       input parameters				     !
c       cnpar - paralel (to magnetic field) component	     !
c               of the refractive index			     !
c       cnper - Re(perpendicular component of the refractive !
c               index)                    		     !
c       cnperi -Im( perpendicular component of the refractive!
c               index)                    		     !
c	iherm=1 Hermitian dielectric tenzor		     !
c	iherm=2 Full full dielectric tenzor		     !
c       output parameter				     !
c       dispmaz-Re(Hamiltonian)			             !
c       hamilti-Im(Hamiltonian)          		     !
c       reps(3,3)-complex dielectric tensor in common /eps/  !
c------------------------------------------------------------
      double precision function
     *dispmuz(xe,ye,te,iherm,cnpar,cnper,cnperi,hamilti)
      implicit double precision (a-h,o-z)
      double precision mode,mu,nz,np2
      double precision ncld(2)	   
      double complex sol,cpp
      double complex chamilmz

      mode=dfloat(1) ! is not used here
      mu=512.44d00/te
      nz=cnpar
      cpp=dcmplx(cnper,cnperi)
      np2=cnper*cnper
      sol=cpp*cpp
      ihermloc=iherm
      call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)
c      write(*,*)'chamilmz',chamilmz
      hamiltr=dreal(chamilmz)
      hamilti=dimag(chamilmz)
      dispmuz=hamiltr
      return
      end


c        ********************hamiltmuz_f***********************
c        * 						      *
c        * this subroutine  calculates: hamiltmz -real part of*
c        *                              Hamiltonian, and      *
c        * reps(3,3) -components of complex dielectric tensor *
c        *            for Mazzucato code                      *
c        * for complex n_perp
c        *****************************************************
c-------------------------------------------------------------
c       input parameters				     !
c       cnpar - parallel (to magnetic field) component	     !
c               of the refractive index			     !
c       cnper - real part of the perpendicular component of
c                the refractive index
c       cnperim -Imaginary part N_perp                        !
c	ihermloc=1 Hermitian dielectric tensor       	     !
c               =2 full      dielectric tensor               !
c       Attention: iherm is in common one.i		     !
c     	rz,r,phi coordinates  of ray point		     !
c	.....						     !
c	rho-small radius was calculated in subroutine: b     !
c            rho is in common one    			     !
c       mode = +1., O-mode   is in common one		     !
c              -1., X-mode				     !
c       output parameter				     !
c       hamiltmz-Hamiltonian				     !
c       reps(3,3)-complex dielectric tensor in common /eps/  !
c------------------------------------------------------------
        subroutine hamiltmuz_ful(cnpar,ihermloc,z,r,phi,cnper,
     .  cnperim,hamiltmz)
	implicit double precision (a-h,o-z)
	double precision mode,mu,nz,np2
	double precision ncld(2)
        include 'param.i'
        include 'one.i'
	dimension an2(2),icuto(2)
	double complex sol,cpp
	double complex chamilmz
	 mode=dfloat(ioxm)
c	 write(*,*)'in hamilmuz mode=',mode

c         write(*,*)' hamiltmuz_ful cnpar,ihermloc,z,r,phi,cnper,cnprim',
c     &   cnpar,ihermloc,z,r,phi,cnper,cnprim

         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
	 te=tempe(z,r,phi,1)
	 mu=512.44d00/te
c	 write(*,*)'in hamiltmuz_ful xe,ye,te',xe,ye,te
         if(dabs(cnpar).lt.1.d-4) then
           if (cnpar.ge.0.d0) then
             sign=1.d0
           else 
             sign=-1.d0
           endif
           cnpar=sign*1.d-4
         endif
	 nz=cnpar

c	 write(*,*)'in hamilmuz cnpar,cnper,ihermloc',
c     *   cnpar,cnper,ihermloc

c	 cpp=dcmplx(cnper,0.0d0)
cSm011206
         cpp=dcmplx(cnper,cnperim)
         np2=cnper*cnper
         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
	 sol=cpp*cpp

c         write(*,*)'hamiltmuz_ful cnper,cnperim,cpp',
c     &   cnper,cnperim,cpp
c         write(*,*)'hamiltmuz_ful mode,xe,ye,mu,nz,np2,sol,ihermloc',
c     &   mode,xe,ye,mu,nz,np2,sol,ihermloc

	 call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)

c	 write(*,*)'hamiltmuz_ful chamilmz',chamilmz
	 hamiltmz=dreal(chamilmz)
c	 write(*,*)'in hamiltmuz_ful hamiltmz=',hamiltmz
	 return
	 end
c---------------------------------------------------------------------
      

