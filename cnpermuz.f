
c        ********************cnpermuz**************************
c        * 						      *
c        * this subroutine finds the perpendicular component  *
c        * of the refractive index from the solution of the   *
c        * dispersion relation D(cnpar,cnper)=0	   using      *
c        * Mazzucato code, using full dielectric tensor with  *
c        * the hermitian and anti-hermitian parts 	      *		
c        *****************************************************
c-------------------------------------------------------------
c       input parameters				     !
c       cnpar - paralel (to magnetic field) component	     !
c               of the refractive index			     !
c       ioxm =1 O mode, =-1 X mode			     !
c	ihermloc=1 Hermit dielectric tenzor,	 	     !
c               =2 full	dielectric tenzor	 	     !
c       Attentin! iherm is in common one.i		     !
c     	rz,r,phi coordinates of ray point (from common one.i)!
c	rho      the small radius (from common one.i) 	     !
c       ioptmaz  is the option for the estimation
c                of the perpendicular refractive index       !
c       ioptmaz=1 np2  = estimation  of perp ref. index      ! 
c                (from coldm)                                !
c       ioptmaz=2   np2  = estimateion of perp ref. index    !
c                   (from input parameter cnper)             !
c       cnper - real part of the perpendicular component of  !
c               the refractive index           		     !
c       output parameter  !
c       cnper - real part of the perpendicular component of  !
c               the refractive index           		     !
c       cnprim- imaginal part of the perpendicular component !
c               of the refractive index        		     !
c       reps(3,3)-complex dielectric tensor		     !
c                 (with antihermitian part) in commmon/eps/  !
c------------------------------------------------------------
        subroutine cnpermuz(cnpar,ihermloc,z,r,phi,cnper,cnprim,ioptmaz)
	implicit double precision (a-h,o-z)
	double precision mode,mu,nz,np2
	double precision ncld(2)
        include 'param.i'
        include 'one.i'
        include 'eps.i'
	dimension an2(2),icuto(2)
	double complex ceps(3,3)
	double complex cpp,sol,chamilmz

        mode=dfloat(ioxm)
        xe=x(z,r,phi,1)         
        ye=y(z,r,phi,1)


        te=tempe(z,r,phi,1)
        mu=512.44d00/te

        nz=cnpar
cSm060225
        call n_par_min_maz(cnpar,nz)

	az2=nz*nz
	accrcy=1.d-06
	naccrcy=500

        if(ioptmaz.eq.1) then 
c---------calculation of nperp for X and O modes
c         estimation of nperp from the cold plasma 
	  call coldm(nz,az2,xe,ye,an2,ncld,icuto)
	  if (ioxm.eq.1) then
	     np2=an2(2) !omode
	  end if
	  if (ioxm.eq.-1) then
	     np2=an2(1) !xmode
	  end if

c-------------------------------------------------------------------
c         it calcultes the solution of dispersion relation 
c         d (n_perp_real,d_nperp_im)=0 using
c         if ihermloc=2 the total Mazzucato tensor
c         if ihermloc=1 the hermitian part of the total Mazzucato tensor       
c----------------------------------------------------------------------
cSm030515
cyup	  write(*,*)'cnpermuz cold an2(1),an2(2),ioxm,np2',
cyup     +    an2(1),an2(2),ioxm,np2

          call complx2(mode,xe,ye,mu,nz,np2,accrcy,naccrcy,ceps,sol,
     +    ihermloc)

c          do i=1,3
c	    do j=1,3
c	       reps(i,j)=ceps(i,j)
c	       write(*,*)'i,j,reps(i,j)',i,j,reps(i,j)
c	    enddo
c          enddo

cyup          write(*,*)'cold zero iteration cnpermuz sol',sol
c          write(*,*)'cnpermuz ceps',ceps

          a1=dreal(sol)
	  cnper=a1
	  a2=dimag(sol)
	  cnprim=a2
c	  write(*,*)'Re(n)=',a1,'Im(n)=',a2
        endif !ioptmaz=1
c------------------------------------------------------------------- 
        if(ioptmaz.eq.2) then 
c---------calculation of the imaginary part cnprim of N_perpendicular
c         from the total dispersion function 
c         using the given real part cnper 
          
c---------calculation of the complex value of the dispersion function
c         from the total hermitian and anti-hermitian parts of the
c         dielectric tensor, using the real N_perp refractive index
c         cnper

          ihermlc1=2
          cpp=dcmplx(cnper,0.0d0)
          np2=cnper*cnper
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermlc1,chamilmz)
c          write(*,*)'cnpermuz 0 ihermloc=2 cpp=cnper',cpp
c          write(*,10)dreal(chamilmz),dimag(chamilmz)
 10       format('total chamilmz iherm=2 initial cpp=cnper'/
     *    2(' ',1pe15.7))
          hamimag=dimag(chamilmz)
          hamr=dreal(chamilmz)

c---------calculation of the numerical derivative from the
c         real part of the hermitian dispersion function
c         d(dispersion_function)/d(N_perp_real)
          ihermlc1=1
          step=1.0d-7
          if(dabs(cnper).gt.1.d0) step=cnper*step
            
          cnperp=cnper+step
          cpp=dcmplx(cnperp,0.0d0)
          np2=cnperp*cnperp
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermlc1,chamilmz)
          hamrp=dreal(chamilmz)
          

          cnperm=cnper-step
          cpp=dcmplx(cnperm,0.0d0)
          np2=cnperm*cnperm
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermlc1,chamilmz)
          hamrm=dreal(chamilmz)
          

          dhamrdnr=(hamrp-hamrm)/(2.d0*step)
          cnprim=-hamimag/dhamrdnr !imaginary part of N_perp

cyup          write(*,*)'cnpermuz.f ipotmaz=2 hamimag,dhamrdnr,cnprim',
cyup     .    hamimag,dhamrdnr,cnprim

ctest
c          ihermlc1=1
c          cpp=dcmplx(cnper,0.0d0)
c          np2=cnper*cnper
c          sol=cpp*cpp
cc          write(*,*)'cnpermuz 1 xe,ye,mu',xe,ye,mu
cc          write(*,*)'cnpermuz 1 nz,np2,ihermlc1',nz,np2,ihermlc1
cc          write(*,*)'cnpermuz 1 mode,sol',mode,sol
c          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermlc1,chamilmz)
cc          write(*,*)'cnpermuz 1 ihermloc=1 cpp=cnper',cpp
cc          write(*,12)dreal(chamilmz),dimag(chamilmz)
 11       format ('chamilmz',2(' ',1pe15.7))
 12       format ('hermitian chamilmz ihermloc=1 for initial cpp=cnper'/
     *    2(' ',1pe15.7))

          ihermloc=2
          cpp=dcmplx(cnper,cnprim)
          np2=cnper*cnper
          sol=cpp*cpp
c---------this call will create dielectric tensor reps and will put this tensor into eps.i 
          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermlc1,chamilmz)
c          write(*,*)'cnpermuz 2 iherml=2 cpp=(cnper,cnprim)',cpp
c          write(*,13)dreal(chamilmz),dimag(chamilmz)
 13       format('total chamilmz iherm=2 cpp=(initial_cnper,cnprim)'/
     *    2(' ',1pe15.7))
 21       format('solutition of Mazzucato solver for full dispesion'/
     *    'using the initial complexNperp=ReN_perp+i*ImN_perp'/
     *    'Im_Nperp=Im(D_full(N_par,ReN_perp)/dReD_full/dReN_perp'/
     *    'sol=',2(' ',1pe15.7))         
cyup           write(*,21)dreal(sol),dimag(sol) 
           goto 300           

c-------- calculation cnper,cnprim from the solution of the system
c         dhamr/dnr*delnr+dhamr/dni*delni=-hamr
c         dhami/dnr*delnr+dhami/dni*delni=-hami
          ihermloc=2
          step=1.0d-7
          cnprim=0.d0
          epsnprim=1.0d-2
          inprmax=1
          inprim=0

 100      continue

          cnperp=cnper+step*cnper
          cpp=dcmplx(cnperp,cnprim)
          np2=cnperp*cnperp
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)
          hamrp=dreal(chamilmz)
          hamip=dimag(chamilmz)
c          write(*,*)'chamilmz',chamilmz
c          write(*,*)'hamrp,hamip',hamrp,hamip
          cnperm=cnper-step*cnper
          cpp=dcmplx(cnperm,cnprim)
          np2=cnperm*cnperm
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)
          hamrm=dreal(chamilmz)
          hamim=dimag(chamilmz) 
c          write(*,*)'chamilmz,2*step*cnper',chamilmz,2*step*cnper
c          write(*,*)'hamrm,hamim',hamrm,hamim 
          dhamrdnr=(hamrp-hamrm)/(2.d0*step*cnper)
          dhamidnr=(hamip-hamim)/(2.d0*step*cnper)
c          write(*,*)'dhamrdnr,dhamidnr',dhamrdnr,dhamidnr
          
          cnprimp=cnprim+step
          cpp=dcmplx(cnper,cnprimp)
          np2=cnper*cnper
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)
          hamrp=dreal(chamilmz)
          hamip=dimag(chamilmz)
c          write(*,*)'chamilmz',chamilmz
c          write(*,*)'hamrp,hamip',hamrp,hamip
          cnprimm=cnprim-step
          cpp=dcmplx(cnper,cnprimm)
          np2=cnper*cnper
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)
          hamrm=dreal(chamilmz)
          hamim=dimag(chamilmz)
c          write(*,*)'chamilmz',chamilmz
c          write(*,*)'hamrm,hamim',hamrm,hamim 
          dhamrdni=(hamrp-hamrm)/(2.d0*step)
          dhamidni=(hamip-hamim)/(2.d0*step)
c          write(*,*)'dhamrdni,dhamidni',dhamrdni,dhamidni
          delt=dhamrdnr*dhamidni-dhamrdni*dhamidnr
          deltr=-(hamr*dhamidni-hamimag*dhamrdni)
          delti=-(dhamrdnr*hamimag-dhamidnr*hamr)
c         write(*,*)'delt,deltr,delti',delt,deltr,delti
          dcnper=deltr/delt
          cnprim=delti/delt
          cnper=cnper+dcnper
c          write(*,*)'3 dcnper,cnper,cnprim',dcnper,cnper,cnprim

          cpp=dcmplx(cnper,cnprim)
          np2=cnper*cnper
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)
c          write(*,*)'cnpermuz 3 ihermloc=2 cpp=(cnper,cnprim)',cpp
c          write(*,14)dreal(chamilmz),dimag(chamilmz)
 14       format('total chamilmz iherm=2 cpp=(cnper,cnprim)'/
     *    2(' ',1pe15.7))
          inprim=inprim+1
          if(inprim.gt.inprmax) then
c             write(*,*)'Warning: cnpermuz inprim > inprmax'
             goto 200
          endif

          if ( dsqrt(dreal(chamilmz)**2+dimag(chamilmz)**2)
     *        .gt.epsnprim) goto 100
 200      continue
         
cendtest
        endif !iopmuz=2
c-------------------------------------------------------------------
      
c        do i=1,3
c	    do j=1,3
c	       reps(i,j)=ceps(i,j)
c	       write(*,*)'i,j,reps(i,j)',i,j,reps(i,j)
c	    enddo
c	enddo
c        write(*,*)'reps_cnpermuz',reps
 300    continue
	return
	end

        subroutine complx2(mode,x,y,mu,nz,np2,accrcy,naccrcy,ceps,sol,
     +  iherm)
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
c     iherm =1 will use the hermitian tensor
c           =2 total tensor
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

        double complex eps1(3,3),eps2(3,3),eps3(3,3),eps(3,3),ceps(3,3)
        double complex c1,c2,c3,so1,so2,solu,det,sss,sol
        double complex csolu2,csolu
c	""""""""""""""""""""""""""""""""""
c	write(*,*)'beg.  complx2 iherm=',iherm
c	""""""""""""""""""""""""""""""""""

        call dcom161(1000,y,mu,nz,iherm)
        call dten161(x,y,mu,nz,eps0,eps1,eps2,eps3,iherm)
c        call dcom16(1000,y,mu,nz)
c        call dten16(x,y,mu,nz,eps0,eps1,eps2,eps3)
        solu=np2
c        write(*,*)'complx2 np2,solu',np2,solu
        do 2010 iter=1,naccrcy
c	write(*,*)'complx2 naccrcy,iter',naccrcy,iter
10      do 2000 i=1,3
        do 2000 j=1,3
c------
        eps(i,j)=eps1(i,j)+solu*eps2(i,j)+solu**2*eps3(i,j)
c        write(*,*)'complx2 nz,eps0',nnz,eps0
c	write(*,*)'complx2 eps(i,j)',eps(i,j)
2000    continue
        c1=eps(1,1)+2.d0*nz*eps(1,3)+eps(1,3)**2+nz**2*eps(3,3)-
     1          eps(1,1)*eps(3,3)
c	write(*,*)'c1=',c1
        c2=eps(1,1)*(nz**2-eps(2,2))-eps(1,2)**2+eps0*(nz**2-eps(1,1))+
     1     2.d0*nz**3*eps(1,3)-2.d0*nz*eps(2,2)*eps(1,3)+
     2     2.d0*nz*eps(1,2)*
     2     eps(2,3)+eps(1,1)*eps(2,3)**2-eps(2,2)*eps(1,3)**2+
     2     2.d0*eps(1,2)*eps(1,3)*eps(2,3)+
     3     nz**2*(eps(1,3)**2-eps(2,3)**2)+nz**2*eps(3,3)*(nz**2-
     4     eps(1,1)-eps(2,2))+eps(3,3)*(eps(1,2)**2+eps(1,1)*eps(2,2))
c	write(*,*)'c2=',c2
        c3=eps0*(nz**2*(nz**2-eps(1,1)-eps(2,2))+eps(1,2)**2+
     1     eps(1,1)*eps(2,2))
c	write(*,*)'c3=',c3

        ddd=dreal(solu)
        sol1=ddd
        det=c2**2-4.d0*c3*c1
c        write(*,*)'complx2 det',det
        det=cdsqrt(det)
        so1=-c2/(2.d0*c1)
c**     write(*,*)'so1=',so1
        so2=det/(2.d0*c1)
c**     write(*,*)'so2=',so2
        ddd=dreal(so2)
        so2r=ddd
c**     write(*,*)'ddd=',ddd
        if(y.gt.1.d0) go to 2400
        c=mode
        go to 2401
2400    c=(-1.d0)*mode
        if(so2r.lt.0.d0) c=mode
2401    solu=so1+c*so2
cc	 write(*,*)'solu=',solu
        ddd=dreal(solu)
        sol2=ddd
        delta=(sol2-sol1)/sol2       
c        write(*,*)'complx2 iter,delta',iter,delta
        delta=dabs(delta)
        if(delta.lt.accrcy) go to 3000
2010    continue
3000    continue !NME 21.03.2005
        if(iter.eq.10) print 566
566     format('  complx2/hot index needs more than 10 iterations')
        if(sol2.lt.0.d0) print 567
567     format(' complx2/ WARNING: hot refractive index is negative!!')
        sol=cdsqrt(solu)
cyup	write(*,*)'comlpx2 sol=',sol,'solu=',solu
      ceps(1,1)=eps(1,1)
      ceps(1,2)=eps(1,2)
      ceps(1,3)=sol*eps(1,3)
      ceps(2,2)=eps(2,2)
      ceps(2,3)=sol*eps(2,3)
      ceps(3,3)=eps0+solu*eps(3,3)
      ceps(2,1)=-ceps(1,2)
      ceps(3,1)=ceps(1,3)
      ceps(3,2)=-ceps(2,3)
c      write(*,*)'in copmlx2'
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
c	write(*,*)'end of complx2'
c	""""""""""""""""""""""""""""""""""

c       write(*,*)'cnpermuz.f complx2 eps',eps
c       write(*,*)'cnpermuz.f complx2 ceps',ceps

        return
        end

      subroutine maphmnri(ye,xe,te,cnpar,cnper,delnper,delnperim,
     * ihermloc)
c     creates the array hamiltar(xe_i,nrprp_j)
c     at the given Y,X,N_par, Te(kev)
c     It uses Mazzucato tensor.
c     It writes the output hamiltonian values into the files:
c     hamr.dat    < (+,-)dlog(dreal(1(+,-)chamilmz)) for given ihermloc
c     hami.dat    < (+,-)dlog(dimag(1(+,-)chamilmz)) for given ihermloc
c     ham0r.dat   < (+,-)dlog(dreal(1(+,-)chamilmz)) for ihermlc1=1 Hermitian
c     ham0i.dat   < (+,-)dlog(dimag(1(+,-)chamilmz)) for ihermlc1=1
c     ham.dat     <  dlog(hamr(i,j)**2+hami(i,j)**2+1)
            
      implicit none 
c-----input
      double precision cnpar,xe,ye,te,cnper,
     * delnper,delnperim  !the segments for nper_real and nper_image
      integer ihermloc !=1 for hermitian =2 for total Mazzucato tensor
      external complx1

      integer nr,ni
      parameter (nr=40,ni=40) !number of points for real and image Nperp 
c      parameter (nr=1,ni=1) !number of points for real and image Nperp 
    
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

c-----output
c     real and image parts of the Hamiltoniuan
      double precision ham0r(nr,ni),ham0i(nr,ni),hamr(nr,ni),hami(nr,ni) 
c-----local
      double precision hnper,hnprim,p
      double precision cnper1,cnprim1
      integer i,j,ihermlc1
      double precision mode,mu,nz,np2
      double complex cpp,sol,chamilmz
      
      if(myrank.ne.0) return

      open(10,file='hamr.dat')
      open(20,file='hami.dat')
      write(10,*)'cnper cnprim hamr'
      write(20,*)'cnper cnprim hami'
      open(30,file='ham0r.dat')
      open(40,file='ham0i.dat')
      open(50,file='ham.dat')
      write(30,*)'cnper cnprim ham0r'
      write(40,*)'cnper cnprim ham0i'
      write(50,*)'cnper cnprim ham'
 1    format (3(' ',1pe15.7))

      mode=1   
      mu=512.44d00/te !te kev
      nz=cnpar
      
      hnper=delnper/nr
      hnprim=delnperim/ni
      do i=1,nr
         cnper1=cnper+hnper*(i-nr/2)
c         cnper1=cnper
         do j=1,ni
            cnprim1=hnprim*(j-1-ni/2)
c            cnprim1=0.d0
            cpp=dcmplx(cnper1,cnprim1)
            np2=cnper1*cnper1
            sol=cpp*cpp
            call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)
c           write(*,*)'chamilmz',chamilmz
            hami(i,j)=dimag(chamilmz)
            hamr(i,j)=dreal(chamilmz)

            ihermlc1=1
          write(*,*)'cnpermuz xe,ye,mu',xe,ye,mu
          write(*,*)'cnpermuz nz,np2,ihermlc1',nz,np2,ihermlc1
          write(*,*)'cnpermuz sol',sol
            cpp=dcmplx(cnper1,cnprim1)
            np2=cnper1*cnper1
            sol=cpp*cpp
            call complx1(mode,xe,ye,mu,nz,np2,sol,ihermlc1,chamilmz)
          write(*,*)'cnpermuz chamilmz',chamilmz
            ham0r(i,j)=dreal(chamilmz)
            ham0i(i,j)=dimag(chamilmz)
               
        write(*,*)'cnper1,cnprim1,hamr(i,j)',cnper1,cnprim1,hamr(i,j)
        write(*,*)'cnper1,cnprim1,hami(i,j)',cnper1,cnprim1,hami(i,j)
        write(*,*)'cnper1,cnprim1,ham0r(i,j)',cnper1,cnprim1,ham0r(i,j)
        write(*,*)'cnper1,cnprim1,ham0i(i,j)',cnper1,cnprim1,ham0i(i,j)

            p=hamr(i,j)
            if (p.ge.0.d0) p=dlog(p+1.d0)
            if (p.lt.0.d0) p=-dlog(1.d0-p)
            write(10,1)cnper1,cnprim1,p
      
            p=hami(i,j)
            if (p.ge.0.d0) p=dlog(p+1.d0)
            if (p.lt.0.d0) p=-dlog(1.d0-p)
            write(20,1)cnper1,cnprim1,p

            p=ham0r(i,j)
            if (p.ge.0.d0) p=dlog(p+1.d0)
            if (p.lt.0.d0) p=-dlog(1.d0-p)
            write(30,1)cnper1,cnprim1,p
      
            p=ham0i(i,j)
            if (p.ge.0.d0) p=dlog(p+1.d0)
            if (p.lt.0.d0) p=-dlog(1.d0-p)
            write(40,1)cnper1,cnprim1,p

            p=hamr(i,j)**2+hami(i,j)**2
            p=dlog(p+1.d0)
            write(50,1)cnper1,cnprim1,p
         enddo
      enddo

      close(10)
      close(20)
      close(30)
      close(40)

      return
      end



      subroutine n_par_min_maz(nll_in,nll_cor)
c     gives the minimal values for n_parallel and N_near the zero value
      implicit none
c     input      
      double precision nll_in !parallel N components 
c     output
      double precision nll_cor!corrected value of N_parrallel 
c     local parameters
      double precision nll_min

c      nll_min=1.d-4   ! should be positive
      nll_min=1.d-3   ! should be positive
     
c      write(*,*)'n_par_min_maz ,nll_min',nll_min
c      write(*,*)'n_par_min_maz nll_in',nll_in
      nll_cor=nll_in
     
      if (dabs(nll_in).le.dabs(nll_min)) then
         if (nll_in.ge.0.d0) then
            nll_cor=dabs(nll_min)
         else 
            nll_cor=-dabs(nll_min)
         endif 
      endif
c      write(*,*)''n_par_min_maz nll_cor',nll_cor
      return
      end
