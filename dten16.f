
************************************************************************
*               this subroutine is for calculating                     *
*               the dielectric tensor in a hot relativistic 
c               mazzucato plasma.                                      *
************************************************************************
        subroutine dten16(x,y,mu,nz,eps0,eps1,eps2,eps3)
c	:::::::::::::Harvey:::::::::::::::::::::::::::::::::::;
c        implicit real (a-h,o-z)
c        implicit integer (i-n), real*8 (a-h,o-z)
c       real nz,mu
c        complex eps1(3,3),eps2(3,3),eps3(3,3)
c        real bes(1),mmu,order
c	::::::::::::::::::::::::::::::::::::::::::::::::
        implicit integer (i-n), real*8 (a-h,o-z)
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
c        write(*,*)'beg.of dten16'
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
************************************************************************
*               calculate anti-hermitian components                    *
************************************************************************
        r2=y*y-1.d0+nz*nz
        if(r2.lt.0.d0) r2=0.d0
c        write(*,*)'dten16 r2=',r2
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
************************************************************************
*               calculate the dielectric tensor                        *
************************************************************************
c        write(*,*)'dten16'
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
cc	write(*,*)'in dten16 a,b,eps2(1,3)',a,b,eps2(1,3)
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
c	write(*,*)'in dten16'
c	do i=1,3
c	  do j=1,3
c	  write(*,*)'i,j,eps1(i,j),eps2(i,j),eps3(i,j)',
c     1	             i,j,eps1(i,j),eps2(i,j),eps3(i,j)
c	  enddo
c	enddo
cc fortest end
c	""""""""""""""""""""""""""
c	write(*,*)'end of dten16'
c	""""""""""""""""""""""""""

        return
        end


      double precision function cosd(a)
	implicit double precision(a-h,o-z)
      pi=datan2(0.d0,-1.d0)
      cosd=dcos(a*pi/180.d0)
      return
      end


      double precision function sind(a)
	implicit double precision(a-h,o-z)
      pi=datan2(0.d0,-1.d0)
      sind=dsin(a*pi/180.d0)
      return
      end

