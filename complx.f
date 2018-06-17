
************************************************************************
*               this subroutine calculates the solution                *
*               of the hot dispersion relation. it uses                *
*               dcom16.sub and dten16.sub.                             *
************************************************************************
        subroutine complx(mode,x,y,mu,nz,np2,accrcy,naccrcy,ceps,sol)
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
c     accrcy=
c     accrcy=iterate for root until successive approximations differ
c           by less than accrcy
c     naccrcy=max number of iterations.
c  output:
c     ceps(3,3)= complex array of dielectric tensor elements.
c                Wave in x-z plane.
c     sol   = complex perpendicular refractive index.
c
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
c        real ddd
        double complex eps1(3,3),eps2(3,3),eps3(3,3),eps(3,3),ceps(3,3)
        double complex c1,c2,c3,so1,so2,solu,det,sss,sol
ccc	double complex cp1,cp2,cp3,cso1,cso2,csolu,csolu2,
ccc     1	cdet,csol,cpp1,cpp2,cpp3,cpeps
c	""""""""""""""""""""""""""""""""""
c	write(*,*)'beg.  complx'
c	""""""""""""""""""""""""""""""""""

        call dcom16(1000,y,mu,nz)
        call dten16(x,y,mu,nz,eps0,eps1,eps2,eps3)
        solu=np2
        do 2010 iter=1,naccrcy
c	write(*,*)'complx naccrcy,iter',naccrcy,iter
10      do 2000 i=1,3
        do 2000 j=1,3
c------
ccc	write(*,*)'i,j,eps1,eps2,eps3,solu,solu**2',
ccc     1  i,j,eps1(i,j),eps2(i,j),eps3(i,j),solu,solu**2
ccc	cp1=eps1(i,j)
ccc	cp2=eps2(i,j)
ccc	cp3=eps3(i,j)
ccc	csolu=solu
ccc	csolu2=solu**2
ccc	cpeps=cp1+csolu*cp2+csolu2*cp3
ccc	write(*,*)'cp1,cp2,cp3,csolu,csolu2',
ccc     1	           cp1,cp2,cp3,csolu,csolu2
ccc        write(*,*)'cp1',cp1
ccc        write(*,*)'csolu*cp2',csolu*cp2
ccc        write(*,*)'csolu2*cp3',csolu2*cp3
ccc        write(*,*)'cpeps',cpeps
ccc	cpp2=csolu*cp2
ccc	cpp3=csolu2*cp3
ccc	cpp1=cp1+cpp2+cpp3
ccc	rcp1=dreal(cp1)
ccc	cicp1=dimag(cp1)
ccc	rcp2=dreal(cp2)
ccc	cicp2=dimag(cp2)
ccc	rcp3=dreal(cp3)
ccc	cicp3=dimag(cp3)
ccc	rcsolu=dreal(csolu)
ccc	cicsolu=dimag(csolu)
ccc	rcsolu2=dreal(csolu2)
ccc	cicsolu2=dimag(csolu2)
ccc	r2=rcsolu*rcp2-cicsolu*cicp2
ccc	ci2=cicsolu*rcp2+rcsolu*cicp2
ccc	r3=rcsolu2*rcp3-cicsolu2*cicp3
ccc	ci3=cicsolu2*rcp3+rcsolu2*rcp3
ccc	reps=rcp1+r2+r3
ccc	cieps=cicp1+ci2+ci3
ccc	write(*,*)'cpp2,cpp3,cpp1',cpp2,cpp3,cpp1
ccc	write(*,*)'rcp1,cicp1',rcp1,cicp1
ccc	write(*,*)'rcp2,cicp2',rcp2,cicp2
ccc	write(*,*)'rcp3,cicp3',rcp3,cicp3
ccc	write(*,*)'r2,ci2',r2,ci2
ccc	write(*,*)'r3,ci3',r3,ci3
ccc	write(*,*)'reps,cieps',reps,cieps
c-------
        eps(i,j)=eps1(i,j)+solu*eps2(i,j)+solu**2*eps3(i,j)
ccc	write(*,*)'eps(i,j)',eps(i,j)
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
        ddd=dreal(solu)
        sol1=ddd
        det=c2**2-4.d0*c3*c1
        det=cdsqrt(det)
        so1=-c2/(2.d0*c1)
c**	write(*,*)'so1=',so1
        so2=det/(2.d0*c1)
c**	write(*,*)'so2=',so2
        ddd=dreal(so2)
        so2r=ddd
c**	write(*,*)'ddd=',ddd
        if(y.gt.1.d0) go to 2400
        c=mode
        go to 2401
2400    c=(-1.d0)*mode
        if(so2r.lt.0.d0) c=mode
2401    solu=so1+c*so2
cc	write(*,*)'solu=',solu
        ddd=dreal(solu)
        sol2=ddd
        delta=(sol2-sol1)/sol2
c	write(*,*)'complx delta=',delta,'accrcy=',accrcy
        delta=dabs(delta)
        if(delta.lt.accrcy) go to 3000
2010    continue
3000    continue !NME 21.03.2005
        if(iter.eq.10) print 566
566     format('  complx/hot index needs more than 10 iterations')
        if(sol2.lt.0.d0) print 567
567     format(' complx/warning: hot refractive index is negative!!')
        sol=cdsqrt(solu)
ccc	write(*,*)'complx sol=',sol,'solu=',solu
ctest
c        write(*,*)'comlx sol',sol
c        write(*,*)'comlx solu',solu
c	csolu2=solu*solu
c	cpeps=csolu2*c1+solu*c2+c3
c        write(*,*)'complx csolu2',csolu2
c        write(*,*)'complx c1',c1
c        write(*,*)'complx c2',c2
c        write(*,*)'complx c3',c3
c        write(*,*)'complx cpeps',cpeps
cendtest
      ceps(1,1)=eps(1,1)
      ceps(1,2)=eps(1,2)
      ceps(1,3)=sol*eps(1,3)
      ceps(2,2)=eps(2,2)
      ceps(2,3)=sol*eps(2,3)
      ceps(3,3)=eps0+solu*eps(3,3)
      ceps(2,1)=-ceps(1,2)
      ceps(3,1)=ceps(1,3)
      ceps(3,2)=-ceps(2,3)
c      write(*,*)'in copmlx'
c      write(*,*)'eps(1,1)=',ceps(1,1)
c      write(*,*)'eps(1,2)=',ceps(1,2)
c      write(*,*)'eps(1,3)=',ceps(1,3)
c      write(*,*)'eps(2,1)=',ceps(2,1)
c      write(*,*)'eps(2,2)=',ceps(2,2)
c      write(*,*)'eps(2,3)=',ceps(2,3)
c      write(*,*)'eps(3,1)=',ceps(3,1)
c      write(*,*)'eps(3,2)=',ceps(3,2)
c      write(*,*)'eps(3,3)=',ceps(3,3)

c	""""""""""""""""""""""""""""""""""
c	write(*,*)'end of complx'
c	""""""""""""""""""""""""""""""""""

        return
        end

