c
      subroutine dskin(initial,energy,pitch,rholoc,
     &fdist,dfdx,dfdpitch,dfdp,i_f)
c..................................................................
c
c     This f77 subroutine reads an ASCII file "diskf" produced by 
c     and described in the subroutine dskout in the 
c     CQL3D Fokker-Planck code.
c     if i_f=1  It calculates the distribution function f
c     if i_f=2  It calculates f and its derivatives dfdx ,dfdpitch and 
c     dfdp=dfdx*c/vnorm (p=momentum/cligth/rest_mass) 
c     in the arbitrary point
c     with the coordinates (energy,pitch,rholoc), x=momentum/rest_mass/Vnorm
c
c     The parameters (ending in "a" below, and dimensions of the 
c     variables to be read into must be set
c     in the file param.i (for dskin.i),
c     in accord with the dimensions (iy,jx,lrz,lrzmax) in "diskf":
c
c     dskin.i:
c     iya=iy_(1), jxa=jx, lrza=lrz, ngena=ngen
c    
c     Bob Harvey, June, 2000
c     Smirnov 23/11/2000; 27/12/2006
c
c     input      initial, energy(KeV), pitch(radian), rholoc(normalized to 1)
c
c     i_diskf is in common one.i
c     i_diskf=0 no usage of the non-maxwellian electron distribution 
c     i_diskf=1 reading the file diskf
c     i_diskf=2 readinf the file netcdfnm.nc 
c     i_diskf=3 analytical calculation of the non-maxwellian disribution
c     i_diskf=4 spline analytic calculation of continuous non-Maxwellian
c               distribution with three temperatures in  three energy ranges
c     i_diskf=5 pure analytic calculation of continuous non-Maxwellian
c               distribution with three temperatures in  three energy ranges
c
c     initial=1 it reads diskf file
c     initial.ne.1 it calculates the distribution function and its derivatives
c
c                At the beginning (for the initialization) it is necessary
c                to call dskin with initial=1    
c
c     output     fdist,dfdx,dfdpitch
c..................................................................
      implicit none
c      implicit integer (i-n), double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'dskin.i'
c      save
c-----input
      integer 
     &initial,             !initial =1 read data (at i_diskf=1,2) or
                           !           create arrays (at i_diskf=3,4)using analytical formulars
                           !initial.ne.1  calculate the distribution function and its derivatives
                           !
     &i_f                  !i_f=1  calculate the distribution function f
                           !   =2  calculate f and its derivatives dfdx ,dfdpitch and
                           !      dfdp=dfdx*c/vnorm (p=momentum/cligth/rest_mass) 
      real*8
     &energy,          !(KeV),
     &pitch,           !(radian),
     &rholoc           !(normalized to 1)

c-----output
      real*8
     &fdist,dfdx,dfdpitch, !distribution function and its derivatives
     &dfdp
c-----locals
      integer initialized
      data initialized/0/
     

      real*8 vnormloc,massloc !(cm/sec),(g) 
      COMMON /dskin1/vnormloc,massloc 
      character*8  idskf
      character*80 line

      integer l,ll,iy_l,k,kk,jx1,n,j_dfdy,j_dfdx,jxx,jxnew,j,i,
     &lpitch,lrr,iyy,i_dfdx,i_dfdy,ir,istat,itrap,iynew

      real*8 clight,del_dfdx,del_dfdy,dens,dfdxnorm,eps_maxw,
     &xdweighl,xdweighu,rweightu,rweightl,weightl,weightu,
     &dfdynorm,dfdy,
     &t_kev,thet0,theta,v0,vel,xvel,pdweighl,pdweighu,
     &pweightu,pweightl,rk_1_kev,rmass_e,rr,ff,xweightu,xweightl
c-----external
      integer luf
      reaL*8 temperho,fdens_fkin,densrho,fdens_fkin1

      save  initialized

cSAP090720
c      real*8 den(lrza,ngena)
      real*8, pointer :: den(:,:) !(lrza,ngena)
     
c      write(*,*)'in dskin initial=',initial,' i_diskf = ',i_diskf
    
       if ((i_diskf.eq.0).or.(i_diskf.eq.5)) return ! in this case the code will use
                                                    ! the analytical maxwellian distribution


       
      IF (initial .eq. 1) THEN

         initialized=1

         if (i_diskf.eq.1) then !reading diskf

c            write(*,*) 'dskin: initial = ',initial
            idskf="diskf"
            open(unit=4,file=idskf,status='old')
c..................................................................
c           The input namelist file has been  transcribed onto the beginning
c           of file idskf and we space over it to the FP data.
c..................................................................

 1          read(4,1003) line
            if (line(1:27).eq."***Begin computed output***") go to 2
            write(*,*)line
            go to 1
 2          continue 

c..................................................................
c
c           CQL3D has a facility to computer FP solutions on a subset of
c           the full radial mesh on which plasma parameters are specified.
c      
c           In the following, we are only interested in the given data
c           which is at flux surfaces numbered ll=1:lrz
c           (thus we can consider lrindx(ll)=ll).
c
c
c           FROM CQL3D:dskout.f
c           In the following disk write to file named idskf:
c           This subroutine is called to write data for each FP'd
c           flux surface.
c           ll=  FP flux surface number 
c           (ll=1:lrz, lrz.le.lrzmax, see cqlinput_help))
c           (lrindx(ll) gives flux surface number on the full CQL3D
c                       radial mesh. (lrindx(ll)=ll if lrz=lrzmax,
c                       and using cql3d mode(cqlpmod="disabled"))).
c           iy_(ll),jx= dimensions in theta and u(momentum/mass)
c                      (In the case where iy_ varies with ll, iy_(1) will be greatest.)
c           lrz= number of flux surfaces FP'd.
c           lrzmax= number of flux surfaces, including any not FP'd.
c           x = momentum-per-mass(nomalized to maximum 1. by vnorm)
c               (same for each flux surface lrindx(ll))
c           y = theta(radians) mesh at  each flux surface lrindx(ll)
c           rovera= normalized radius (ll)
c                   (rho, see Hinton and Haseltine for non-circ).
c                   (generally ~sqrt(tor. flux), but other coords are
c                   available for setup in the CQL3D input file.)
c           elecfld = toroidal electric field (volts/cm)
c           bthr - the poloidal magnetic field at theta-poloidal = pi/2.
c           btoru - the toroidal magnetic field at the same position. 
c           bthr0, btor0, are the poloidal and  toroidal
c                  magnetic fields at the outer midplane of the flux surface. 
c           reden= electron density at minimum B point on flux surface.
c           temp= initial electron temperature (keV)
c           radmin= plasma minor radius (cms).
c           vnorm= normalization momentum-per-mass (maximum on grid) (cm/sec)
c           vmaxdvt= vnorm/(temp/mass)**0.5
c           eovedd= electric  field, normalized to Driecer field 
c                   (calc'd in sub restvty).
c
c           distribution function normalised so that
c           integral( (dx)**3 f) = density at minumum B point.
c
c..................................................................
c
            l=1
c           Loop back to 10 reading data for each flux surfaces

cSAP090720
c 10         read(4,1004)  ll, iy_(l),jx,lrz,lrzmax,ngen
 10         read(4,1004)  ll, iy_l,jx,lrz,lrzmax,ngen

cSAP090720
c            l=ll+1            
c..................................................................
c           Check dimensions from dskin.i
c..................................................................
cSAP090720
c            write(*,*)'iy_l,iya,jx,jxa,lrz,lrza,ngen,ngena',
c     +      iy_(1),iya,jx,jxa,lrz,lrza,ngen,ngena
cSAP090808  delete ngena,jxa,lrza
c            write(*,*)'iy_l,iya,jx,jxa,lrz,lrza,ngen,ngena',
c     +      iy_l,iya,jx,jxa,lrz,lrza,ngen,ngena
            write(*,*)'iy_l,iya,jx,lrz,ngen',
     +      iy_l,iya,jx,lrz,ngen

cSAP090720
c            if (iy_(1).ne.iya .or. jx.ne.jxa .or. lrz.ne.lrza
c     +      .or. ngen.ne.ngena) stop 'Check dskin.i parameters'
             if (iy_l.gt.iya) then 
                write(*,*)'dskin.f iy_l.gt.iya'
                write(*,*)'iy_l',iy_l,'iya',iya
                stop 'Check dskin.i parameters'
             endif

            if (ll.eq.1) then
c-------------allocate pointers in dskin.i file 
              call ainalloc_dskin_i

              if (lrzmax.gt.lrz) then
                 write(*,*)'in dskin.f before dskinainalloc_dskin_i'
                 write(*,*)'lrzmax > lrz'
                 write(*,*)'lrzmax=',lrzmax,'lrz=',lrz
                 write(*,*)'set small radial dimension in pointers'
                 write(*,*)'as lrz'
                 write(*,*)'it can create problem with pointers' 
              endif
            endif
c----------------------------------------------------
            iy_(l)=iy_l 
            l=ll+1            
c------------------------------------------------
            write(*,*)'ll, iy_(ll),jx,lrz,lrzmax,ngen'
            write(*,1004)ll, iy_(ll),jx,lrz,lrzmax,ngen


            read(4,1004)  itl(ll),itu(ll)
            write(*,*)'itl(ll),itu(ll)'
            write(*,1004) itl(ll),itu(ll)

            read(4,1005)  (x(j),j=1,jx)

            write(*,*)'x(j) j=1,jx'
            write(*,1005)  (x(j),j=1,jx)

            read(4,1005)  (y(i,ll),i=1,iy_(ll))

            write(*,*)'iy_(ll)',iy_(ll)
            write(*,*)'y(i,ll) i=1,jy_(ll)'
            write(*,1005)  (y(i,ll),i=1,iy_(ll))

            do 20 k=1,ngen
               write(*,*)'20 k=',k
               read(4,1005)  bnumb(k),fmass(k)
               write(*,*) 'bnumb(k),fmass(k)'
               write(*,1005) bnumb(k),fmass(k)

               read(4,1005)  rovera(ll),elecfld(ll),bthr(ll),btoru(ll)
               write(*,*)  'rovera(ll),elecfld(ll),bthr(ll),btoru(ll)'
               write(*,1005) rovera(ll),elecfld(ll),bthr(ll),btoru(ll)

               read(4,1005)  bthr0(ll),btor0(ll),reden(ll,k),temp(ll,k)
               write(*,*)' bthr0(ll),btor0(ll),reden(ll,k),temp(ll,k)'
               write(*,1005)bthr0(ll),btor0(ll),reden(ll,k),temp(ll,k)

               read(4,1005)  radmin,vnorm,vmaxdvt,eovedd
               write(*,*)'  radmin,vnorm,vmaxdvt,eovedd'
               vnormloc=vnorm   
               kk=1 !for ions (now for electron) !????
c              kk=3 !for electrons
               massloc=fmass(kk)
               write(*,1005)  radmin,vnorm,vmaxdvt,eovedd
               write(*,*)'dskin vnormloc',vnormloc
                read(4,1005)  ((f(i,j,ll,k),i=1,iy_(ll)),j=1,jx)
 20         continue

            if (ll.lt.lrz) go to 10      
      
 1003 format(a80)
 1004 format(16i5)
 1005 format(5e16.8)
 1006 format(a27)
            close(unit=4)
         endif! (i_diskf.eq.1) reading diskf

         if (i_diskf.eq.2) then ! reading netcfnm.nc file
cSm030519 the next line will be commented only for test under windows
cSAP090720
c            call netcdfr3d(netcdfnm,iymax,iy_,jx,lrz,y,x,rovera,vnorm,f) 
c--------------------------------------------------
c          allocate pointers in dskin.i
c--------------------------------------------------
c            call ainalloc_dskin_i
c            write(*,*)'in dskin.f after call ainalloc_dskin_i'
c--------------------------------------------------------
            write(*,*)'in dskin.f before call netcdfr3df'
cSAP110316
c           call netcdfr3d(netcdfnm)
            call netcdfr3d

            write(*,*)'in dskin.f after call netcdfr3df'
            vnormloc=vnorm
            write(*,*)'dskin vnormloc',vnormloc
            rmass_e=9.1094d-28         !electron rest mass (g)
            massloc= rmass_e
            lrzmax=lrz
 
c            do ll=1,lrz 
c              do j=1,jx 
c                do i=1,iy_(ll)
c                    write(*,*)'i,j,ll,f(i,j,ll,1)',i,j,ll,f(i,j,ll,1)
c                enddo
c              enddo
c            enddo

         endif !(i_diskf.eq.2) reading netcdfnm.nc file

         if (i_diskf.eq.3) then ! creation of the  non-Maxwellian distribution
c--------------------------------------------------
c           allocate pointers in dskin.i
c--------------------------------------------------
            call ainalloc_dskin_i
c--------------------------------------------------
            call fokker(rtem0,
     &             rtail,r1t,r2t,ttail,rhot,r1h,r2h,thotpar,thotper,
     &             hotexp,hotmnpar,hotmxpar,hotmnper,hotmxper,
     &             rbeam,r1b,r2b,tbeam,ebeam,thbeam,
     &             jx,iym,lrz,ngen,
cSAP090720
c     &             jxa,iya,lrza,ngena,
     &             jx,iya,lrz,ngen,
     &             f,x,iy_,y,rovera,vnorm)
               
            vnormloc=vnorm
            rmass_e=9.1094d-28         !electron rest mass (g)
            massloc= rmass_e
            write(*,*)'after fokker vnormloc,massloc',vnormloc,massloc

            kk=1
            lrzmax=lrz
ctest comparison with analytical maxwellian distrib
c          clight =2.99792458d10
c          pi=4.d0*datan(1.d0)
c          dk_1_kev=1.6022d-9   
c          do ir=1,lrz  
c             rho=rovera(ir)             
c             T_kev=temperho(rho,kk)
c             theta=massloc*clight**2/(dk_1_kev*T_kev)
c             write(*,*)'dskin ir,rho,temperho(rho,kk),theta',
c     .        ir,rho,temperho(rho,kk),theta
c             call besk2as(theta,bk2)
c             write(*,*)'ir,theta',ir,theta
c             do i=1,iy_(ir)
c                do j=1,jx
c                  p=x(j)*vnorm/clight
c                  ps=p*p
c                  gamma=dsqrt(1.d0+ps)                 
c                  f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)
c      write(*,*)'ir,i,j,x(j),gamma,f_maxw(vnorm/clight)**3,f(i,j,ir,kk)'
c     .      ,ir,i,j,x(j),gamma,f_maxw*(vnorm/clight)**3,f(i,j,ir,kk)
c                enddo
c             enddo
c          enddo
c          stop
cendtes

         endif !(i_diskf.eq.3)
      
         if (i_diskf.eq.4) then !analytic calculation of continuous non-Maxwellian
                               !distribution with three temperatures in three
                               !energy ranges creation of the non-Maxwellian distribution

c--------------------------------------------------
c          allocate pointers in dskin.i
c--------------------------------------------------
           call ainalloc_dskin_i
c--------------------------------------------------
           call fokker_3_temperature(rtem0,
     &     rtemp1, rtemp2, rtemp3,
     &     rvtovte1,rvtovte2,
cSAP0907230
c     &     jx,iym,lrz,jxa,iya,lrza,ngena,
     &     jx,iym,lrz,jx,iya,lrz,ngen,
     &     f,x,iy_,y,rovera,vnorm)
               
            vnormloc=vnorm
            rmass_e=9.1094d-28         !electron rest mass (g)
            massloc= rmass_e
            write(*,*)'after fokker_3_temperature vnormloc,massloc',
     &                vnormloc,massloc

            kk=1
            lrzmax=lrz
ctest comparison with analytical maxwellian distrib
c          clight =2.99792458d10
c          pi=4.d0*datan(1.d0)
c          dk_1_kev=1.6022d-9   
c          do ir=1,lrz  
c             rho=rovera(ir)             
c             T_kev=temperho(rho,kk)
c             theta=massloc*clight**2/(dk_1_kev*T_kev)
c             write(*,*)'dskin ir,rho,temperho(rho,kk),theta',
c     .        ir,rho,temperho(rho,kk),theta
c             call besk2as(theta,bk2)
c             den1=densrho(rho,1)

c             write(*,*)'ir,theta,den',ir,theta,den
c             do i=1,iy_(ir)
c                do j=1,jx
c                  p=x(j)*vnorm/clight
c                  ps=p*p
c                  gamma=dsqrt(1.d0+ps)                 
c               f_maxw=den1*theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)
c      write(*,*)'ir,i,j,x(j),gamma,f_maxw(vnorm/clight)**3,f(i,j,ir,kk)'
c     .      ,ir,i,j,x(j),gamma,f_maxw*(vnorm/clight)**3,f(i,j,ir,kk)
c                enddo
c             enddo
c          enddo
c          stop
cendtes

         endif !(i_diskf.eq.4)
      

c--------test calculations of the density from the distribution function
c        at all radial points
c--------allocate den
         allocate( den(1:lrz,1:ngen),STAT=istat)
         call bcast(den,0.d0,SIZE(den))
 
c--------calculate the density (for test) from the distribution: f
c        array den(*,*) will be used in spline distribution

c         write(*,*)'dskin before dens_fkin'  
c         do ll=1,lrz 
c            do j=1,jx 
c               do i=1,iy_(ll)
c                  write(*,*)'i,j,ll,f(i,j,ll,1)',i,j,ll,f(i,j,ll,1)
c               enddo
c            enddo
c         enddo

         call dens_fkin(f,iy_,jx,lrz,ngen,iya,jx,lrz,ngen,
     +   x,y,den)
         
         do n=1,ngen
            do k=1,lrz
               write(*,*)'n,k,den(k,n)',n,k,den(k,n)
            enddo
         enddo  

c         write(*,*)'dskin after dens_fkin'  
c         do ll=1,lrz 
c            do j=1,jx 
c               do i=1,iy_(ll)
c                  write(*,*)'i,j,ll,f(i,j,ll,1)',i,j,ll,f(i,j,ll,1)
c               enddo
c            enddo
c         enddo

         deallocate( den,STAT=istat)

         goto 30
         

c-----------creation of the test analytical distribution f(i,j,k,n)       
c           usage of the maxwellian electron distribution 
c-----------calculations of normloc; setting massloc=mass_e
            clight =2.99792458d10      !speed of light     (cm/sec)
            rmass_e=9.1094d-28         !electron rest mass (g)
            massloc=rmass_e
            rk_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)
            T_kev=temperho(0.d0,1)
            theta=rmass_e*clight**2/(rk_1_kev*T_kev)
c            eps_maxw=1.d-20
            eps_maxw=1.d-40
            write(*,*)'dskin theta',theta
            vnormloc=clight*dsqrt(-dlog(eps_maxw)/theta)
            write(*,*)'dskin 0 vnormloc cm/sec',vnormloc
c           vnormloc=clight*dsqrt((1.d0-dlog(eps_maxw)/theta)**2-1.d0)
c           write(*,*)'dskin 1 vnormloc cm/sec',vnormloc
            vnorm=vnormloc
c           write(*,*)'in dskin vnormloc/clight',vnormloc/clight
c           write(*,*)'in dskin vnormloc=',vnormloc
            fmass(1)=rmass_e 
                    
            call set_f(iya,jx,lrz,ngen,iy_,jx,lrz,ngen,y,x,rovera,f)
            write(*,*)'test analytical distribution created'
   
c           call dens_fkin(f,iy_,jx,lrz,ngen,iya,jxa,lrza,ngena,
c     +     x,y,den)
cSAP090720
c            dens=fdens_fkin(f,iy_,jx,lrz,1,14,iya,jxa,lrz,ngen,
c     +      x,y)
            dens=fdens_fkin(f,iy_,jx,lrz,1,14,iya,jx,lrz,ngen,
     +      x,y)
            write(*,*)'dskin.f dens=fdens_fkin at k=14,n=1',dens
          dens=fdens_fkin1(i_diskf,f,iy_,jx,lrz,1,14,iya,jx,lrz,ngen,
     +    x,y,rovera)
            write(*,*)'dskin.f dens=fdens_fkin1 at k=14,n=1',dens
c-----------calculate the density (for test) from the distribution: f 
            call dens_fkin(f,iy_,jx,lrz,ngen,iya,jx,lrz,ngen,
     +      x,y,den)
c-------------------------------------------------------    

 30      continue
c.............................................
c        the calculation of the numerical derivatives from 
c        the distribition function f: df_dpitch,df_dx in the
c        intermedium points.
c        df_dpitch(i,j,k,n)=df/dpitch(pitch_i-0.5,x_j,rho_k,n) 1<i<iy_(k) 1=<j=<jx
c        df_dpitch(1,j,k,n)=df_dpitch(pitch_1)=0 1=<j=<jx
c        df_dpitch(i,j,k,n)=df_dpitch(pitch_iy)=0 1=<j=<jx
c        df_dx(i,j,k,n)=df/dx(pitch_i,x_j-0.5,rho_k,n)) 1<j<jx 1=<i=<iy_(k)
c        df_dx(i,1,k,n)=df/dx(pitch_i,x_1,rho_k,n))=0 1=<i=<iy_(k)
c        df_dx(i,jx,k,n)=df/dx(pitch_i,x_jx,rho_k,n))=0 1=<i=<iy_(k)
c............................................. 
         
c         call derivative(f,iy_,jx,lrz,ngen,iya,jxa,lrza,ngena,
c     +   iya1,jxa1,x,y,df_dpitch,df_dx)

c.............................................
c        calculate the spline coeficients for the distribution finction
c        of kk specie
         kk=1
         write(*,*)'kk',kk

cSm060228
         goto 40 
         do ir=1,lrz
            do j=1,jx
c               do i=1,iya
               do i=1,iy_(kk)
	          if (f(i,j,ir,kk).gt.1.d-50) then
                     f(i,j,ir,kk)=dlog(f(i,j,ir,kk))
                  else
                     f(i,j,ir,kk)=-50.d0*dlog(10.d0)
                  endif
               enddo
            enddo
         enddo
 40      continue

ctest calculation of the derivatives
c         pause 'testing calculations of the derivatives in dskin'

c         write(*,*)'dskin before derivatives'  
c         do ll=1,lrz 
c            do j=1,jx 
c               do i=1,iy_(ll)
c                  write(*,*)'i,j,ll,f(i,j,ll,1)',i,j,ll,f(i,j,ll,1)
c               enddo
c            enddo
c         enddo

         jx1=jx+1
         call derivative(f,iy_,jx,lrz,ngen,iya,jx,lrz,ngen,
     +   iya1,jx1,x,y,df_dpitch,df_dx)
cendtest
c         write(*,*)'dskin before splcoef_distrib'  
c         do ll=1,lrz 
c            do j=1,jx 
c               do i=1,iy_(ll)
c                  write(*,*)'i,j,ll,f(i,j,ll,1)',i,j,ll,f(i,j,ll,1)
c               enddo
c            enddo
c         enddo
         call splcoef_distrib(f,iy_,jx,lrz,kk,y,x)
         write(*,*)'spline coefficients for distribution were created'
         write(*,*)'dskin rovera',rovera

ctest spline
cSm060228
        ngen=1
        do k=1,lrz
c           denn=den(k,kk)
           rho=rovera(k)
           dfdxnorm=0.d0
           dfdynorm=0.d0
           del_dfdx=0.d0
           del_dfdy=0.d0
           i_dfdx=0
           j_dfdx=0
           i_dfdy=0
           j_dfdy=0
           do j=1,jx
              do i=1,iy_(k)
                 thet0=y(i,k)
                 v0=x(j)
c        call distn(lrz,iy_,jx,kk,thet0,v0,y,x,f,ff,dfdx,dfdy,denn,den)
     
         call distr(lrz,iy_,jx,kk,thet0,v0,y,x,f,ff,dfdy,dfdx,
     +   rho,rovera)

c        write(*,*)'v0,thet0',v0,thet0  
c        write(*,*)'dskin,k,j,i,f(i,j,k,kk),ff',k,j,i,f(i,j,k,ngen),ff
c        write(*,*)'dskin dfdy,dfdx',dfdy,dfdx

        if (j.gt.1) v0=x(j)-(x(j)-x(j-1))*0.5d0
        if (i.gt.1) thet0=y(i,k)-(y(i,k)-y(i-1,k))*0.5

c        write(*,*)'v0,thet0',v0,thet0

        call distr(lrz,iy_,jx,kk,thet0,v0,y,x,f,ff,dfdy,dfdx,
     +  rho,rovera)
        dfdxnorm=dfdxnorm+(df_dx(i,j,k,kk)-dfdx)**2
        dfdynorm=dfdynorm+(df_dpitch(i,j,k,kk)-dfdy)**2
        if (del_dfdx.lt.dabs(df_dx(i,j,k,kk)-dfdx)) then
            i_dfdx=i
            j_dfdx=j
        endif

        if (del_dfdy.lt.dabs(df_dpitch(i,j,k,kk)-dfdy)) then
            i_dfdy=i
            j_dfdy=j
        endif
         
c        write(*,*)'df_dx(i,j,k,kk),dfdx',df_dx(i,j,k,kk),dfdx
c        write(*,*)'df_dpitch(i,j,k,kk),dfdy',df_dpitch(i,j,k,kk),dfdy

               enddo
            enddo
            dfdxnorm=dsqrt(dfdxnorm)/dfloat(jx*iy_(k))
            dfdynorm=dsqrt(dfdynorm)/dfloat(jx*iy_(k))

c            write(*,*)'k,dfdxnorm,dfdynorm',k,dfdxnorm,dfdynorm
c            write(*,*)'i_dfdx,j_dfdx',i_dfdx,j_dfdx
c            write(*,*)'i_dfdy,j_dfdy',i_dfdy,j_dfdy

c        stop 'stop in dskin test spline k=1' 
        enddo
c        stop 'stop in dskin test spline'
cendtest spline

      ELSE     !  END OF INITIALIZATION

c..................................................................
c        BEGIN: Interpolate f for value at energy,pitch,rholoc
c        Single species, although easy to generalize
c..................................................................

c        write(*,*) 'dskin: initial = ',initial
         if (initialized.ne.1) stop 'Stop in dskin:Data not initialized'
         kk=1
c         kk=1 !for ions (now for electron)
c        kk=3 !for electrons
         clight=2.99792458d10
c        write(*,*)'in dskin energy,kk,fmass(kk)',energy,kk,fmass(kk)
c         vel=dsqrt(2.d0*energy*1.6022d-9/fmass(kk))
        
         vel=dsqrt(2.d0*energy*1.6022d-9/massloc)
         xvel=vel/vnormloc
         goto 100 !to spline

c------- First, determine if the particle at (energy,pitch,radius)
c        is trapped.
c        Use lookup to find weights for distribution at neighboring
c        tabulated points.  Use closest pitch angle mesh to
c        determine whether trapped.

         if (rholoc .le. rovera(1)) then
            ll=1
         elseif (rholoc .ge. rovera(lrz)) then
            ll=lrz
         else
            call lookup(rholoc,rovera,lrzmax,weightu,weightl,ll)
            if (weightl.ge.weightu) ll=ll-1
         endif

         lpitch=luf(pitch,y(1,ll),iy_(ll))
         if (lpitch.gt.(iy_(ll)/2)) lpitch=lpitch-1 ! Gives symmetry.
         itrap=0                                   ! Trapping flag.

c        if (lpitch.gt.itl(ll) .and. lpitch.lt.itu(ll)) itrap=1

c        Shift radial location of trapped particles by a half banana width
c        This will be an outward/inward shift for pos/neg bthr*bnumb
c        (bnumb is charge number, including sign).
c        Transiting particles are unshifted.

         rr=rholoc
c         charge=4.8032d-10
c         if (itrap.eq.1) then
c            qb_mc=bnumb(kk)*charge*bthr(ll)/(fmass(kk)*clight)
c            rban=xvel*cos(pitch)*vnorm/qb_mc

c           Shift radial location of trapped particles by a half banana width
c           This will be an outward/inward shift for pos/neg bthr*bnumb
c           (bnumb is charge number, including sign).
c            rr=rholoc+rban/radmin
c            if (rr.lt. 0.d0) rr=0.d0
c            if (rr.gt. 1.d0)  rr=1.d0
c         endif

c..................................................................
c        Use lookup to find weights for distribution at neighboring
c        tabulated points.  Interpolate to velocity and pitch angle
c        to the specified (xvel,pitch) values.
c..................................................................

c        Stay on the grid:

c        radius
         call lookup(rr,rovera,lrz,rweightu,rweightl,lrr)

c         write(*,*)'rr,lrr,rovera(1),rovera(2),rovera(lrr)',
c     .   rr,lrr,rovera(1),rovera(2),rovera(lrr)
c         write(*,*)'0 rweightu,rweightl',rweightu,rweightl

         if (lrr.gt.lrz) then
            lrr=lrz
            rweightl=0.d0
            rweightu=1.d0
         elseif (lrr.le.1) then
            lrr=2
            rweightl=1.d0
            rweightu=0.d0
         endif
c         write(*,*)'1 rweightu,rweightul,lrr',rweightu,rweightl,lrr

c        pitch angle
         call lookup(pitch,y(1,ll),iy_(ll),pweightu,pweightl,iyy)

         if (i_f.eq.2) then ! for the derivatives calculations
            call lookupd (pitch,y(1,ll),iy_(ll),iyy,iynew,pdweighu,
     +                    pdweighl)
         endif

         if (iyy.ge.iy_(ll)) then
            iyy=iy_(ll)
            pweightu=1.d0
            pweightl=0.d0
            pdweighu=1.d0
            pdweighl=0.d0
         elseif (iyy.le.1) then
            iyy=2
            pweightu=0.d0
            pweightl=1.d0
            pdweighu=0.d0
            pdweighl=1.d0
         endif

c        momemtum u=p/m
          call lookup(xvel,x,jx,xweightu,xweightl,jxx)
         if (i_f.eq.2) then ! for the derivatives calculations
            call lookupd (xvel,x,j,jxx,jxnew,xdweighu,xdweighl)
         endif

         if (jxx.ge.jx) then
            jxx=jx
            xweightu=1.d0
            xweightl=0.d0
            xdweighu=1.d0
            xdweighl=0.d0
         elseif (jxx.le.1) then
            jxx=2
            xweightu=0.d0
            xweightl=1.d0
            xdweighu=0.d0
            xdweighl=1.d0
         endif

         
c        Tri-linear interpolate (Eliminate f<0 values).
         
         fdist=rweightl*
     +             (pweightl*xweightl*dmax1(0d0,f(iyy-1,jxx-1,lrr-1,1))
     +               +pweightu*xweightl*dmax1(0d0,f(iyy,jxx-1,lrr-1,1))
     +               +pweightl*xweightu*dmax1(0d0,f(iyy-1,jxx,lrr-1,1))
     +               +pweightu*xweightu*dmax1(0d0,f(iyy,jxx,lrr-1,1)))
     +     +rweightu*(pweightl*xweightl*dmax1(0d0,f(iyy-1,jxx-1,lrr,1))
     +               +pweightu*xweightl*dmax1(0d0,f(iyy,jxx-1,lrr,1))
     +               +pweightl*xweightu*dmax1(0d0,f(iyy-1,jxx,lrr,1))
     +               +pweightu*xweightu*dmax1(0d0,f(iyy,jxx,lrr,1)))

         write(*,*)'the interpolation fdist',fdist
 
c        Following write statments for code checkout.
c         write(*,*)'iyy,jxx,lrr',iyy,jxx,lrr
c         write(*,*)'pitch,xvel,r r',pitch,xvel,rr
c         write(*,*) 'dskin: fdist, f.....................:',fdist
c        fdista=distrib1(pitch,xvel,rr)
c        write(*,*) 'test analytical fdista',fdista

c        write(*,*) 'dskin: rweightl,pweightl,xweightl'
c        write(*,*) '           ',rweightl,pweightl,xweightl
c        write(*,*) 'f(iyy-1,jxx-1,lrr-1,1),f(iyy,jxx-1,lrr-1,1)',
c     +  f(iyy-1,jxx-1,lrr-1,1),f(iyy,jxx-1,lrr-1,1)
c        write(*,*) 'f(iyy-1,jxx,lrr-1,1),f(iyy,jxx,lrr-1,1)',
c     +  f(iyy-1,jxx,lrr-1,1),f(iyy,jxx,lrr-1,1)
c        write(*,*) 'dskin: rweightu', rweightu
c        write(*,*) 'f(iyy-1,jxx-1,lrr,1),f(iyy,jxx-1,lrr,1)',
c     +  f(iyy-1,jxx-1,lrr,1),f(iyy,jxx-1,lrr,1)
c        write(*,*) 'f(iyy-1,jxx,lrr,1),f(iyy,jxx,lrr,1)',
c     +  f(iyy-1,jxx,lrr,1),f(iyy,jxx,lrr,1)
c        write(*,*) 'dskin..................................'

c        END of evaluation of f(energy,pitch,rholoc)

         if (i_f.eq.2) then ! the derivatives calculations

c---------- Tri-linear interpolate (of the derivatives fd/dpitch df/dmomentum

            dfdpitch=rweightl*
     +             (pdweighl*xweightl*df_dpitch(iynew-1,jxx-1,lrr-1,1)
     +               +pdweighu*xweightl*df_dpitch(iynew,jxx-1,lrr-1,1)
     +               +pdweighl*xweightu*df_dpitch(iynew-1,jxx,lrr-1,1)
     +               +pdweighu*xweightu*df_dpitch(iynew,jxx,lrr-1,1))
     +     +rweightu*(pdweighl*xweightl*df_dpitch(iynew-1,jxx-1,lrr,1)
     +               +pdweighu*xweightl*df_dpitch(iynew,jxx-1,lrr,1)
     +               +pdweighl*xweightu*df_dpitch(iynew-1,jxx,lrr,1)
     +               +pdweighu*xweightu*df_dpitch(iynew,jxx,lrr,1))

c           write(*,*)'rweightl,rweightu',rweightl,rweightu
c           write(*,*)'pdweighl,pdweighu',pdweighl,pdweighu

            dfdx=rweightl*
     +           (pweightl*xdweighl*df_dx(iyy-1,jxnew-1,lrr-1,1)
     +               +pweightu*xdweighl*df_dx(iyy,jxnew-1,lrr-1,1)
     +               +pweightl*xdweighu*df_dx(iyy-1,jxnew,lrr-1,1)
     +               +pweightu*xdweighu*df_dx(iyy,jxnew,lrr-1,1))
     +      +rweightu*(pweightl*xdweighl*df_dx(iyy-1,jxnew-1,lrr,1)
     +               +pweightu*xdweighl*df_dx(iyy,jxnew-1,lrr,1)
     +               +pweightl*xdweighu*df_dx(iyy-1,jxnew,lrr,1)
     +               +pweightu*xdweighu*df_dx(iyy,jxnew,lrr,1))

            dfdp=dfdx*clight/vnorm

c-----------checkout the interpolation of the derivatives from the distribution 
c         write(*,*)'the interpolation of the deriv dfdpitch,dfdx,dfdp',
c     +     dfdpitch,dfdx,dfdp 
        
c-----------test analytical derivatives
c           call ddistrib(pitch,xvel,rr,dfdpitcha,dfdua)
c           write(*,*) 'test analytical derivatives dfdpitcha,dfdua',
c     +     dfdpitcha,dfdua
      
         endif ! of evaluation of derivatifes dfdpitch,dfdu,dfdp
c------------------------------------------------------      
 100     continue !spline distribution

cSm060310
         if(xvel.ge.1.d0)then
            fdist=0.d0
            dfdp=0.d0
            dfdpitch=0.d0
            goto 110    
         endif

         denn=densrho(rholoc,kk)
c        rite(*,*)'denn',denn             
c        write(*,*)'fdist,dfdx,dfdy',fdist,dfdx,dfdy

c        call distn(lrz,iy_,jx,kk,pitch,xvel,y,x,f,ff,dfdy,dfdx,
c     +  denn,den)

c         write(*,*)'dskin lrz,jx,kk,pitch,xvel,rholoc',
c     +   lrz,jx,kk,pitch,xvel,rholoc
         call distr(lrz,iy_,jx,kk,pitch,xvel,y,x,f,ff,dfdy,dfdx,
     +   rholoc,rovera)
c        write(*,*)'dskin after distr ff,dfdy,dfdx',ff,dfdy,dfdx

         fdist=ff
c        write(*,*)'1 fdist,dfdx,dfdy',fdist,dfdx,dfdy
         dfdp=dfdx*clight/vnorm
         dfdpitch=dfdy
c        write(*,*)'spline fdist',fdist      
c        write(*,*)'spline 1 of the deriv dfdpitch,dfdx,dfdp',
c     +  dfdpitch,dfdx,dfdp 
             
cend spline------------------------------------------
cSm060228
c         fdist=dexp(fdist)
c         dfdp=dfdp*fdist 
c         dfdpitch=dfdpitch*fdist 

cSm060228   
         if (fdist.lt.0.d0)then
            fdist=0.d0
            dfdp=0.d0
            dfdpitch=0.d0    
         endif

110      continue
c         write(*,*)'spline 2 fdist,dfdpitch,dfdp',fdist,dfdpitch,dfdp
      ENDIF          
    
      return
      end



      subroutine derivative(f,iy_,jx,lrz,ngen,iya,jxa,lrza,ngena,
     +iya1,jxa1,x,y,df_dpitch,df_dx)
c.....calculates the numerical derivatives df_dpitch and df_dx
c     from the distribution function f
c     by pitch angle and momentum x=u/vnorm=p/m/vnorm:  
      implicit none
c.....input
      integer iya,iya1,jxa,jxa1,lrza,ngena! iya1=iya+1,jxa1=jxa+1
      integer jx   ! the number of mesh points in momentum mesh
      integer lrz  ! the number of mesh points in radial mesh
      integer ngen ! the total number of the plasma species 
      real*8  f(iya,jxa,lrza,ngena) ! the distribution function
      integer iy_(lrza) ! the number of mesh points in pitch angle mesh  
      real*8 x(jxa) ! the momentum mesh
      real*8 y(iya,lrza) !the pitch angle mesh
     
c.....output
      real*8  df_dpitch(iya1,jxa,lrza,ngena)
      real*8  df_dx(iya,jxa1,lrza,ngena)
c.....local
      integer i,j,k,n
      
c.....the derivatives df/dpitch angle
       write(*,*)'dskin in derivative jx,lrz,ngen,iya,jxa,lrza,ngena',
     &                                jx,lrz,ngen,iya,jxa,lrza,ngena

c       do k=1,lrz 
c          do j=1,jx 
c             do i=1,iy_(k)
c                write(*,*)'i,j,k,f(i,j,k,1)',i,j,k,f(i,j,k,1)
c             enddo
c          enddo
c      enddo

      do n=1,ngen
        do k=1,lrz
          do j=1,jx
             do i=1,iy_(k)+1
              if ((i.eq.1).or.(i.eq.iy_(k)+1)) then
                 df_dpitch(i,j,k,n)=0.d0
              else 
c----------------calculation of the derivative df_dpitch at y=y_(i-0.5)
                 df_dpitch(i,j,k,n)=
     +           (dmax1(0d0,f(i,j,k,n))-dmax1(0d0,f(i-1,j,k,n)))/
     +                              (y(i,k)-y(i-1,k))
              endif
            enddo
          enddo
        enddo
      enddo
      
c.....the derivatives df/d_momentum (u/m)

      do n=1,ngen
        do k=1,lrz
          do i=1,iy_(k)
            do j=1,jx+1
              if ((j.eq.1).or.(j.eq.jx+1)) then
                 df_dx(i,j,k,n)=0.d0
              else
c----------------calculation of the derivative df_x  at x=x_(j-0.5)
                 df_dx(i,j,k,n)=
     +           (dmax1(0.d0,f(i,j,k,n))-dmax1(0.d0,f(i,j-1,k,n)))/
     +                              (x(j)-x(j-1))
              endif
            enddo
          enddo
        enddo
      enddo

      return
      end  


      subroutine lookupd (x,xarray,length,iement,
     + iement_new,dweightu,dweightl)
c..................................................................
c     This interpolates the derivatives
c     to gain a bit of accuracy.
c     x is the argument; xarray is the monotonic array; length
c     is the length of the array. lement is the first index such
c     that xarray(lement).gt.x.
c..................................................................
      implicit none
      save
c     input
      double precision x,xarray
      integer length,iement
      dimension xarray(*)
c     output
      double precision dweightu,dweightl
      integer iement_new ! the number of the element from the aray of derivative 
c     local
      double precision xc,x_imin,x_ipl 
       
      if  (iement.eq.1) then
          dweightl=1.d0
          iement_new=1
      else
          xc=0.5d0*(xarray(iement)+xarray(iement-1)) !x_(iement-0.5)         
          if  (x.gt.xc) then
c------------ x_(iement-0.5) < x 
              x_imin=xc             !x_(iement-0.5)
              if  (iement.eq.length) then
                  iement_new=length+1  
                  if  (x.gt.xarray(length)) then
                      dweightl=0.d0
                  else
                      dweightl=(xarray(length)-x)/
     +                         (xarray(length)-x_imin)
                  endif
              else    
                  x_ipl=0.5d0*(xarray(iement)+xarray(iement+1)) !x_(iement+0.5)
                  dweightl=(x_ipl-x)/(x_ipl-x_imin)
                  iement_new=iement+1
              endif 
          else
c-------------x_(i-1)=< x = < x_(iement-0.5) 
              x_ipl=xc      
              if  (iement.eq.2) then
                  dweightl=(x_ipl-x)/(x_ipl-xarray(1))
                  iement_new=2
              else 
                  x_imin=0.5d0*(xarray(iement-1)+xarray(iement-2)) !x_(i-1.5)
                  dweightl=(x_ipl-x)/(x_ipl-x_imin)
                  iement_new=iement
              endif   
          endif
      endif
               
      dweightu=1.d0-dweightl

      return
      end

      subroutine set_f(iya,jxa,lrza,ngena,iy_,jx,lrz,ngen,y,x,rovera,f)
c     put the test distribution to array f
      implicit integer (i-n), double precision (a-h,o-z)
      integer iya,jxa,lrza,ngena,jx,lrz,ngen
      double precision f(iya,jxa,lrza,ngena),
     +y(iya,lrza),x(jxa),rovera(lrza)
      integer iy_(lrza)
c     local
      integer n,k,j,i

c      write(*,*)'in set_f ngen,lrza,lrz,jxa,jx',ngen,lrza,lrz,jxa,jx

      do n=1,ngen
c         write(*,*)'n=',n
         do k=1,lrz
c            write(*,*)'k=',k,'iy_(k)',iy_(k)  
            do j=1,jx
c               write(*,*)'j=',j
               do i=1,iy_(k)
c                  write(*,*)'i,j,k',i,j,k
c                  write(*,*)'y(i,k),x(j),rovera(k)',
c     +                       y(i,k),x(j),rovera(k)
c                  f(i,j,k,n)=distrib(y(i,k),x(j),rovera(k))
                   f(i,j,k,n)=distrib1(y(i,k),x(j),rovera(k))
c                   write(*,*)'f',f(i,j,k,n)
               enddo
            enddo
         enddo         
      enddo   

      return
      end

   
      double precision function distrib(pitch,u,rho)
c     test distribution function 
      implicit none
c     input
      double precision pitch,u,rho
c     local
      double precision p

      p=5.d0 

      distrib=0.1d0+0.9d0*(1-rho**2)*
     +        (1.d0+(dcos(pitch))**2)*dexp(-(p*u)**2)
      return
      end

      double precision function distrib1(pitch,u,rho)
c     test relativistic maxwellian electron distribution function 
      implicit none
c     input
      double precision pitch,u,rho !u=p/(m*vnorm)
      double precision vnormloc,massloc
      COMMON /dskin1/vnormloc,massloc  !v (cm/sec),m( g)

c     external
      double precision densrho,temperho      
c     local
      double precision dense,T_kev,theta,clight,mass_e,k_1_kev,
     +gamma,ps,bk2,pi,f_maxw 
  
      dense=densrho(rho,1) !electron density  non-dimensional in 10**19 (m**-3)
      T_kev=temperho(rho,1)!electron temperature KeV
      pi=4.d0*datan(1.d0)
      clight =2.99792458d10     !speed of light     (cm/sec)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)
      theta=massloc*clight**2/(k_1_kev*T_kev)       
      ps=(u*vnormloc/clight)**2
      gamma=dsqrt(1.d0+ps)
    
 
c     calculation the Macdonalds function bk2=K_2(theta)*EXP(theta)	
      call besk2as(theta,bk2)

      f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)

      distrib1=dense*f_maxw*(vnormloc/clight)**3

cfor_test dense=1
c      distrib1=f_maxw*(vnormloc/clight)**3
c      write(*,*)'distrib1 (vnormloc/clight)**3',(vnormloc/clight)**3
c      distrib1=f_maxw
    
      return
      end



      subroutine ddistrib(pitch,u,rho,dfdpitch,dfdu)
c     analytical derivatives from the test analytical distribution function
c     dfdpitch, dfdu
      implicit none
c     input 
      double precision pitch,u,rho
c     output
      double precision dfdpitch,dfdu
c     local
      double precision p
      p=5.d0 

      dfdpitch=0.9d0*(1-rho**2)*
     +            (-2.d0*dcos(pitch)*dsin(pitch))*dexp(-(p*u)**2)
      dfdu=0.9d0*(1-rho**2)*
     +            (1.d0+(dcos(pitch))**2)*dexp(-(p*u)**2)*(-2.d0*p*p*u)
      return
      end


      subroutine dervtest(i_diskf,energy,pitch,rho,dfdpitch,dfdenergy)
c     calculates the derivatives from the distribution f in the given poit
c     as the central derivatives using dskin
      implicit none
c     input 
      double precision energy,pitch,rho !eV,radian,normalized rasdius
      integer i_diskf
c     i_diskf=0 no usage of the non-maxwellian electron distribution 
c     i_diskf=1 reading the file diskf
c     i_diskf=2 readinf the file netcdfnm.nc 
c     i_diskf=3 analytical calculation of the non-maxwellian disribution
c     i_diskf=4 analytic calculation of continuous non-Maxwellian distribution
c            with three temperatures in  three energy ranges
c     output 
      double precision dfdpitch,dfdenergy
c     local
      double precision step,energyp,energym,pitchp,pitchm,fdistp,fdistm,
     +dfdxl,dfdpitchl,dfdp 
      integer initial

      initial=0
      step=1.d-7       

      energyp=energy*(1.d0+step)
      energym=energy*(1.d0-step)
      write(*,*)'in dervtest energyp,pitch,rho',energyp,pitch,rho
      call dskin(initial,energyp,pitch,rho,
     +fdistp,dfdpitchl,dfdxl,dfdp,1)
      write(*,*)'in dervtest fdistp',fdistp

      write(*,*)'in dervtest energym,pitch,rho',energym,pitch,rho
      call dskin(initial,energym,pitch,rho,
     +fdistm,dfdpitchl,dfdxl,dfdp,1)
      write(*,*)'in dervtest fdistm',fdistm

      dfdenergy=(fdistp-fdistm)/(2*energy*step)
      write(*,*)'dfdenergy',dfdenergy

      pitchp=pitch*(1.d0+step)
      pitchm=pitch*(1.d0-step)
      call dskin(initial,energy,pitchp,rho,
     +fdistp,dfdpitchl,dfdxl,dfdp,1)
      call dskin(initial,energy,pitchm,rho,
     +fdistm,dfdpitchl,dfdxl,dfdp,1)
      dfdpitch=(fdistp-fdistm)/(2*pitch*step)
      write(*,*)'dfdpitch',dfdpitch
      return
      end
 
      real*8 function 
     +fdens_fkin(f,iy_,jx,lrz,n,k,iya,jxa,lrza,ngena,
     +x,y)
c     calculates the density from the distribution f (given as array)
c     in the momentum space x=p/m/vnorm
c     in one radial mesh point with the number k 
c     for the  plasma specie n
c     fdens_fkin=2*pi*integral{0,pi}sin(pitch)*dpitch[ integral{0,x(jx)}x**2*dx ]
c 
      implicit none
c     input 
      integer iya,jxa,lrza,ngena
      integer jx   ! the number of mesh points in momentum mesh
      integer lrz  ! the number of mesh points in radial mesh
      integer n    ! the number of the plasma species
      integer k    ! the number of the radial point
      real*8  f(iya,jxa,lrza,ngena) ! the distribution function
      integer iy_(lrza)   ! the number of mesh points in pitch angle mesh  
      real*8  x(jxa)      ! the momentum mesh
      real*8  y(iya,lrza) !the pitch angle mesh. It depends on the radial point. 
      
c     output 
      
c     local
      real*8  dens ! the integral at the given radius 
      real*8  coefx,coefpitch,pi 
      integer j,i

      pi=4.d0*datan(1.d0)
      
      dens=0.d0
      do j=1,jx

          if ((j.ne.1).and.(j.ne.jx)) then
              coefx=x(j)**2*(x(j+1)-x(j-1))/2.d0
          else            
               if (j.eq.1) then
                  coefx=1.d0*x(2)**3/24.d0
              else
                  !j=jx
                  coefx=(1.d0/3.d0)*
     +                  (x(jx)**3-(0.5d0*(x(jx)+x(jx-1)))**3)
              endif
          endif                      

          do i=1,iy_(k)

              if ((i.ne.1).and.(i.ne.iy_(k))) then
                  coefpitch=dsin(y(i,k))*0.5d0*(y(i+1,k)-y(i-1,k))
              else
                  if (i.eq.1) then
                      coefpitch=1.d0-dcos(0.5d0*(y(2,k)+y(1,k)))
                  else
                      ! i = iy_(k)
                      coefpitch=1.d0+
     +                dcos(0.5d0*(y(iy_(k),k)+y(iy_(k)-1,k)))
                  endif
              endif
  
              dens=dens+f(i,j,k,n)*coefx*coefpitch
                    
          enddo !i
      enddo     !j

      dens=2.d0*pi*dens
  
      fdens_fkin=dens

      return
      end






      subroutine dens_fkin(f,iy_,jx,lrz,ngen,iya,jxa,lrza,ngena,
     +x,y,den)
c     calculates the density from the distribution f(given as array f(iya,jxa,lrza,ngena))
c     in the momentum space u=p/m
c     for all radial mesh points
      implicit none
c-----input 
      integer iya,jxa,lrza,ngena
      integer jx   ! the number of mesh points in momentum mesh
      integer lrz  ! the number of mesh points in radial mesh
      integer ngen ! the total number of the plasma species 
      double precision  f(iya,jxa,lrza,ngena) ! the distribution function
      integer iy_(lrza) !t the number of mesh points in pitch angle mesh  
      double precision x(jxa) ! the momentum mesh
      double precision y(iya,lrza) !the pitch angle mesh
c-----external
      double precision fdens_fkin       

c-----output 
      double precision den(lrz,ngen) !density distribution for all species  
c-----local
      double precision dens_kin ! the integral at the given radius
 
      integer n,k,j,i

      do n=1,ngen
        write(*,*)'in densfkin n=',n
        do k=1,lrz          
           dens_kin=fdens_fkin(f,iy_,jx,lrz,n,k,iya,jxa,lrza,ngena,
     +                         x,y)
           den(k,n)=dens_kin
           write(*,*)'n,k,dens_kin',n,k,dens_kin
        enddo      !k
      enddo        !n

      return
      end

c
      double precision function 
     +fdens_fkin1(i_diskf,f,iy_,jx,lrz,n,k,iya,jxa,lrza,ngena,
     +x,y,rovera)
c     calculates the density from the distribution f (given as array)
c     in the momentum space u=p/m/vnorm
c     in one radial mesh point with the number k 
c     for the  plasma specie n
      implicit none

c     input 

      integer i_diskf
c     i_diskf=0 no usage of the non-maxwellian electron distribution 
c     i_diskf=1 reading the file diskf
c     i_diskf=2 readinf the file netcdfnm.nc 
c     i_diskf=3 analytical calculation of the non-maxwellian disribution
c     i_diskf=4 analytic calculation of continuous non-Maxwellian distribution
c               with three temperatures in  three energy ranges
      integer iya,jxa,lrza,ngena
      integer jx   ! the number of mesh points in momentum mesh
      integer lrz  ! the number of mesh points in radial mesh
      integer n    ! the number of the plasma species
      integer k    ! the number of the radial point
      double precision  f(iya,jxa,lrza,ngena) ! the distribution function
      integer iy_(lrza) ! the number of mesh points in pitch angle mesh  
      double precision x(jxa) ! the momentum mesh
      double precision y(iya,lrza) !the pitch angle mesh
      double precision rovera(lrza)
c     output 
      
c     local
      double precision dens ! the integral at the given radius 
      double precision coefx,coefpitch,pi 
      integer j,i

      double precision dens1,energy,k_1_kev,fdist,dfdx,dfdpitch,dfdp
      double precision vnormloc,massloc !velocity(cm/sec),  mass(g)
      COMMON /dskin1/vnormloc,massloc
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)

      pi=4.d0*datan(1.d0)

      dens=0.d0
      dens1=0.d0
      write(*,*)'in fdens_fkin1:k,rovera(k)',k,rovera(k)

      do j=1,jx

          if ((j.ne.1).and.(j.ne.jx)) then
              coefx=x(j)**2*(x(j+1)-x(j-1))/2.d0
          else            
               if (j.eq.1) then
                  coefx=1.d0*x(2)**3/24.d0
              else
                  !j=jx
                  coefx=(1.d0/3.d0)*
     +                  (x(jx)**3-(0.5d0*(x(jx)+x(jx-1)))**3)
              endif
          endif                      
          energy=0.5d0*massloc*(vnormloc*x(j))**2/k_1_kev
          do i=1,iy_(k)

              if ((i.ne.1).and.(i.ne.iy_(k))) then
                  coefpitch=dsin(y(i,k))*0.5d0*(y(i+1,k)-y(i-1,k))
              else
                  if (i.eq.1) then
                      coefpitch=1.d0-dcos(0.5d0*(y(2,k)+y(1,k)))
                  else
                      ! i = iy_(k)
                      coefpitch=1.d0+
     +                dcos(0.5d0*(y(iy_(k),k)+y(iy_(k)-1,k)))
                  endif
              endif
  
              dens=dens+f(i,j,k,n)*coefx*coefpitch

c              call dskin(0,energy,y(i,k),rovera(k),
c     +                   fdist,dfdx,dfdpitch,dfdp,1)
c              dens1=dens1+fdist*coefx*coefpitch    
 
          enddo !i
      enddo     !j

      dens=2.d0*pi*dens
c      dens1=dens1*2.d0*pi
c      write(*,*)'in fdens_fkin1: dens,dens1',dens,dens1
      fdens_fkin1=dens

      return
      end


      subroutine test_fdist
c-----calculates the distribution in the 3d mesh points using dskin
c     and compares the results with the array f            
      implicit integer (i-n), double precision (a-h,o-z)
      include 'param.i'
      include 'dskin.i'
      include 'one.i'
      double precision vnormloc,massloc !(cm/sec),(g) 
      COMMON /dskin1/vnormloc,massloc 

      i_f=1
      initial=0
       
      rmass_e=9.1094d-28         !electron rest mass (g)
      massloc= rmass_e

      ir=1
      i=1
      j=10

      write(*,*)'!!!!!!!!!!!!!!!! test fdist'
      do ir=1,lrz
         rholoc=rovera(ir)
         do i=1,iy_(ir)
            pitch=y(i,ir)
            do j=1,jx
               energy=0.5d0*massloc*(x(j)*vnormloc)**2/
     .                (1.6d-12*1.d+3) !KeV
               write(*,*)'dskin test_fdist ir,i,j',ir,i,j
c               write(*,*)'dskin test_fdist rholoc,pitch,x(j),x*vnorm',
c     .         rholoc,pitch,x(j),x(j)*vnormloc
c               write(*,*)'dskin test_fdist energy',energy
               call dskin(initial,energy,pitch,rholoc,
     +fdist,dfdx,dfdpitch,dfdp,i_f)
               write(*,*)'dskin test_fdist fdist,f(i,j,ir,1)',
     +         fdist,f(i,j,ir,1)
            enddo
         enddo
      enddo

      return
      end


      subroutine fmaxw_en(energy,rho,fmaxw,dfdx)
c-----calculate relativistic maxwellian function on energy (KeV)
c     and rho (0,1)
      implicit none
c-----input
      double precision energy ! (KEV)
      double precision rho    ! small radius (0,1)     

      double precision vnormloc,massloc !(cm/sec),(g) 
      COMMON /dskin1/vnormloc,massloc 
c-----external      
      double precision distrib1
c-----output
      double precision fmaxw,dfdx
       
c-----local
      double precision u ,pitch
      pitch=0.d0
      
      u=dsqrt(2.d0*energy/massloc*1.6002d-9)/vnormloc
c      fmaxw=distrib1(pitch,u,rho)
      call ddistrib1(pitch,u,rho,fmaxw,dfdx)    
      return
      end

      subroutine ddistrib1(pitch,u,rho,distrib,dfdu)
c     test relativistic maxwellian electron distribution function 
      implicit none
c     input
      double precision pitch,u,rho !u=p/(m*vnorm)
      double precision vnormloc,massloc
      COMMON /dskin1/vnormloc,massloc  !v (cm/sec),m( g)
c     output
      double precision distrib,dfdu
c     external
      double precision densrho,temperho      
c     local
      double precision dense,T_kev,theta,clight,mass_e,k_1_kev,
     +gamma,ps,bk2,pi,f_maxw 
  

      dense=densrho(rho,1) !electron density  non-dimensional in 10**19 (m**-3)
      T_kev=temperho(rho,1)!electron temperature KeV
      pi=4.d0*datan(1.d0)
      clight =2.99792458d10     !speed of light     (cm/sec)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)
      theta=massloc*clight**2/(k_1_kev*T_kev)       
      ps=(u*vnormloc/clight)**2
      gamma=dsqrt(1.d0+ps)
     
c     calculation the Mackdonalds function bk2=K_2(theta)*EXP(theta)	
      call besk2as(theta,bk2)

      f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)

      distrib=dense*f_maxw*(vnormloc/clight)**3
      dfdu=distrib*(-theta*dsqrt(ps)/gamma)*(vnormloc/clight)
      return
      end


      subroutine splcoef_distrib(f,iy_,jx,lrz,ngen,y,x)
c-----calculate the set of 2D spline coefficients (in the velocity space)
c     for the table distribution function f(y=pitch,x=p/(m*vnorm),rho)
c     for all radial points.
c     rho=rho(ir) ir=1,lrz
c     x=x(j)      j=1,jx
c     y=y(i,ir)      i=1,iy_(ir)
c
c     input: 
c     f(iya,jxa,lrza,ngena) is the distribution function
c
c     ngena                is the maximal number of the plasma species
c                          given in param.i
c     iya,jxa,lrza         are the maximal dimentions given in param.i
c
c     iy_(lrxa)             (integer*4)is the number of the pitch angle
c                          mesh points at the each small radius
c     jx                   is the number of the velocity=p/(m*vnorm) 
c                          mesh points
c     lrz                  is the number of the small radius 
c                          mesh points
c     ngen                 is the number of used plasma specie
c
c     y(iya,lrz)           is the mesh of the pitch angle at 
c                          different radial points
c
c     x(jx)                is the mesh of momentum x=p/(m*vnorm)   
c
c     output:
c     second derivativers
c     fxx(iya,jxa,lrza,ngena) d^2f/dpitch
c     fyy(iya,jxa,lrza,ngena) d^2f/dx
c   
c     fourth derivatives                         
c     fxxyy(iya,jxa,lrza,ngenan) d^4f/(d/dpitch)**2*(d/dx)**2*

      implicit none
c      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'spline_distrib.i'
c-----input
      integer jx,lrz,ngen
cSAP090808 delete ngega,jxa,lrza  
c      integer iy_(lrza)
c      real*8 f(iya,jxa,lrza,ngena),y(iya,lrza),x(jxa)
      integer iy_(lrz)
      real*8 f(iya,jx,lrz,ngen),y(iya,lrz),x(jx)
c-----local 
      integer ir,j,idm,i
      real*8 ff,df_dx,df_dy,theta0,v0 
c-----externals
      real*8 terp2p_Sm

      SAVE
      
c***********************************************************************
c     Distribution function related declarations:
c      common/spline_distrib/ fxx(iya,jxa,lrza,ngena),
c     &fyy(iya,jxa,lrza,ngena),fxxyy(iya,jxa,lrza,ngenan),
c     &ftemp(iya,jxa),ibd(4),wk(nwka)
c      here x direction is the pitch angle
c           y direction is the velocity=x=p/(m*vnorm)
c***********************************************************************

c now set up the splines used to interpolate for the distn
c function in the calculation of emissivity and absorptivity
c

c-----2D  spline boundary conditions
cSm030519
c      idm=iy_(ir)
cSm040111      
c      idm=iya


      ibd(1)=2
      ibd(2)=2
c      ibd(1)=4
c      ibd(2)=4
      ibd(3)=4
      ibd(4)=4

c      write(*,*)'splcoef_distrib jx,lrz,ngen',jx,lrz,ngen
c      write(*,*)'iy_',(iy_(i),i=1,lrz)
c      write(*,*)'x',(x(j),j=1,jx)
c      do ir=1,lrz  
c      write(*,*)'ir=',ir,'y',(y(i,ir),i=1,iy_(ir))
c      enddo
c      do ir=1,lrz 
c        do j=1,jx 
c          do i=1,iy_(ir)
c            write(*,*)'i,j,ir,f(i,j,ir,1)',i,j,ir,f(i,j,ir,1)
c          enddo
c       enddo
c      enddo
         
      do ir=1,lrz
cSm040111
c         do j=1,jxa
         do j=1,jx
            fxx(1,j,ir,ngen)=0.d0
            fxx(iy_(ir),j,ir,ngen)=0.d0
         enddo

cSm040110         
          idm=iya
 
c         call coeff2(iy_(ir),y(1,ir),jx,x,f(1,1,ir,ngen),
c     1   fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),fxxyy(1,1,ir,ngen),
c     1   idm,ibd,wk)

c        idm=iy_(ir)
         idm=iya
         call coeff2_Sm(iy_(ir),y(1,ir),jx,x,f(1,1,ir,ngen),
     1   fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),fxxyy(1,1,ir,ngen),
     1   idm,ibd,wk,iya)
c-------------------------------------------------------        
c         write(*,*)'ir,iy_(ir),jx',ir,iy_(ir),jx
c         do j=1,jx
c         write(*,*)'j,f(i,j,ir,ngen)',j
          
c         write(*,*)(f(i,j,ir,ngen),i=1,iy_(ir))
c         write(*,*)'fxx(i,j,ir,ngen)'
c         write(*,*)(fxx(i,j,ir,ngen),i=1,iy_(ir))
c         write(*,*)'fxxyy(i,j,ir,ngen)'
c         write(*,*)(fxxyy(i,j,ir,ngen),i=1,iy_(ir))
c         write(*,*)'fxx(i,j,ir,ngen)'
c         write(*,*)(fxx(i,j,ir,ngen),i=1,iy_(ir))
c        enddo   
c---------------------------------------------------------

ctest
c      idm=iya
cc      idm=iy_(ir) !this case works wrong 
c      do j=1,jx
c         do i=1,iy_(ir)
c            v0=x(j)
c            theta0=y(i,ir)
cc       write(*,*)'ir,j,i,v0,theta0',ir,j,i,v0,theta0
cc      write(*,*)'splcoef_distrib ir,ngen,fxx(1,1,ir,ngen)',
cc     +ir,ngen,fxx(1,1,ir,ngen)
c      ff=terp2p(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1          f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1          fxxyy(1,1,ir,ngen),idm,0,0,1)
cSm040111
c      ff=terp2p_Sm(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1          f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1          fxxyy(1,1,ir,ngen),idm,0,0,1,iya)
c       write(*,*)'splcoef_distrib ir,j,i,v0,theta0,ff,f(i,j,ir,ngen)',
c     1 ir,j,i,v0,theta0,ff,f(i,j,ir,ngen)
c       df_dx=terp2p(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1          f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1          fxxyy(1,1,ir,ngen),idm,0,1,1)
cSm040111
c       df_dx=terp2p_Sm(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1          f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1          fxxyy(1,1,ir,ngen),idm,0,1,1,iya)
c       df_dy=terp2p(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1          f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1          fxxyy(1,1,ir,ngen),idm,1,0,1)
cSm04011
c       df_dy=terp2p_Sm(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1          f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1          fxxyy(1,1,ir,ngen),idm,1,0,1,iya)
c       write(*,*)'df_dx(vel),df_dy',df_dx,df_dy
c         enddo       
c      enddo 
c      stop 'splcoef_distrib'
cendtest      
      enddo
c      stop 'splcoef_distrib'
      return
      end


      subroutine distn(lrz,iy_,jx,ngen,
     .theta0,v0,y,x,f,ff,dfdx,dfdy,denn,den)

c-----calculate distribution function and its derivatives using 2D
c     cubic spline in the velocity space and linear interpolation
c     in rho direction using the density
c   
c     input
c     iya,jxa,lrza  are the maximal numbers of meshs points (seted in param.i)
c     ngena         is the maximal number of plasma species (seted in param.i)
c
c     iy_(lrza)      is the number of pitch pints at each radius 
c     jx            the number of the momentum x=p/(m*vnorm) mesh points  
c     lrz           is the number of the radial points
c     ngen          is the number of plasma specie
c
c     theta0,v0      are the point in velocity space at which interpolation 
c                   is requied
c
c     y(iya,lrza)   is the pitch angle mesh
c     x(jxa)        is the momentum mesh =P/(m*vnorm)
c
c     f(iya,jxa,lrza,ngena) is the distribution function
c     
c     den(lrza,ngena)is the density at the radial mesh
c 
c     denn           is the density in the given radial point
c
c-----output
c     ff,dfdy,dfdx function and its derivatives

      implicit double precision (a-h,o-z)

      include 'param.i'
      include 'spline_distrib.i'
c     spline_distrib.i contains the following common block
c     common/spline_distrib/ fxx(iya,jxa,lrza,ngena),
c     &fyy(iya,jxa,lrza,ngena),fxxyy(iya,jxa,lrza,ngena),
c     &ftemp(iya,jxa),ibd(4),wk(nwka)

c-----input
cSAP090808 delete ngena,jxa,lrza
c      integer iy_(lrza),jx,lrz,ngen
c      integer iy_(lrza),jx,lrz,ngen
c      real*8 theta0,v0,y(iya,lrza),x(jxa),
c     & f(iya,jxa,lrza,ngena),den(lrza,ngena) 
      integer iy_(lrz),jx,lrz,ngen 
      real*8 theta0,v0,y(iya,lrz),x(jx),
     & f(iya,jx,lrz,ngen),den(lrz,ngen)
c-----local
      integer i,ir,im,idm 
c-----output
      real*8 ff,dfdx,dfdy
c-----external
      real*8 terp2p,terp2p_Sm
     
c     Could speed up by eliminating repeated searches of grid
c     in terp2, for x and y.
c
c  interpolate (using density) to get distribution functions
c  Location alogrithm for interpolation points assumes that density
c  is a decreasing function of radius.
c
      do i=1,lrz
         ir=i
         if(denn.gt.den(ir,ngen)) goto 2                
      enddo
2     idm=iya
        
c     
c      write(*,*)'distn ir denn,den(ir,ngen)',ir,denn,den(ir,ngen)

      ff=(denn/den(ir,ngen))*terp2p(theta0,v0,iy_(ir),y(1,ir),jx,x,
     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
     1      fxxyy(1,1,ir,ngen),idm,0,0,1)
      dfdx=(denn/den(ir,ngen))*terp2p(theta0,v0,iy_(ir),y(1,ir),jx,x,
     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
     1      fxxyy(1,1,ir,ngen),idm,1,0,0)
      dfdy=(denn/den(ir,ngen))*terp2p(theta0,v0,iy_(ir),y(1,ir),jx,x,
     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
     1      fxxyy(1,1,ir,ngen),idm,0,1,0)

cSm040111 use spline with maximal length iya>=iy_(ir) direction
c      idm=iy_(ir)
c      ff=(denn/den(ir,ngen))*terp2p_Sm(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1      fxxyy(1,1,ir,ngen),idm,0,0,1,iya)
c      dfdx=(denn/den(ir,ngen))*terp2p_Sm(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1      fxxyy(1,1,ir,ngen),idm,1,0,0,iya)
c      dfdy=(denn/den(ir,ngen))*terp2p_Sm(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1      fxxyy(1,1,ir,ngen),idm,0,1,0,iya)

c
      if(ir.eq.1.or.(ir.eq.lrz.and.denn.le.den(lrz,ngen))) return
      im=ir-1
      ffm=(denn/den(im,ngen))*terp2p(theta0,v0,iy_(im),y(1,im),jx,x,
     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
     1      fxxyy(1,1,im,ngen),idm,0,0,0)
      dfdxm=(denn/den(im,ngen))*terp2p(theta0,v0,iy_(im),y(1,im),jx,x,
     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
     1      fxxyy(1,1,im,ngen),idm,1,0,0)
      dfdym=(denn/den(im,ngen))*terp2p(theta0,v0,iy_(im),y(1,im),jx,x,
     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
     1      fxxyy(1,1,im,ngen),idm,0,1,0)

cSm040111 iuse soline with maximal length iya>iy_(ir) direction
c      idm=iy_(im)
c      ffm=(denn/den(im,ngen))*terp2p_Sm(theta0,v0,iy_(im),y(1,im),jx,x,
c     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
c     1      fxxyy(1,1,im,ngen),idm,0,0,0,iya)
c      dfdxm=(denn/den(im,ngen))*terp2p_Sm(theta0,v0,iy_(im),y(1,im),jx,x,
c     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
c     1      fxxyy(1,1,im,ngen),idm,1,0,0,iya)
c      dfdym=(denn/den(im,ngen))*terp2p_Sm(theta0,v0,iy_(im),y(1,im),jx,x,
c     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
c     1      fxxyy(1,1,im,ngen),idm,0,1,0,iya)

c
      ddm=denn-den(im,ngen)
      ddp=den(ir,ngen)-denn
      dd=ddm+ddp
      ff=( ffm*ddp + ff*ddm ) / dd
      dfdx=( dfdxm*ddp + dfdx*ddm ) / dd
      dfdy=( dfdym*ddp + dfdy*ddm ) / dd
c  
      return
      end

c!!!
      subroutine distr(lrz,iy_,jx,ngen,
     .theta0,v0,y,x,f,ff,dfdx,dfdy,rho,rovera)

c-----calculate distribution function and its derivatives using 2D
c     cubic spline in the velocity space and linear interpolation
c     in rho direction using the density
c   
c     input
c     iya,jxa,lrza  are the maximal numbers of meshs points (seted in param.i)
c     ngena         is the maximal number of plasma species (seted in param.i)
c
c     iy_(lrza)      is the number of pitch pints at each radius 
c     jx            the number of the momentum x=p/(m*vnorm) mesh points  
c     lrz           is the number of the radial points
c     ngen          is the number of plasma specie
c
c     theta0,v0      are the point in velocity space at which interpolation 
c                   is requied
c
c     y(iya,lrza)   is the pitch angle mesh
c     x(jxa)        is the momentum mesh =P/(m*vnorm)
c
c     f(iya,jxa,lrza,ngena) is the distribution function
c
c     rovera(lrza)   is the radial mesh
c     
c     rho            is the normaized small radius    
c

c-----output
c     ff,dfdy,dfdx function and its derivatives

      implicit double precision (a-h,o-z)

      include 'param.i'
      include 'spline_distrib.i'
c     spline_distrib.i contains the following common block
c     common/spline_distrib/ fxx(iya,jxa,lrza,ngena),
c     &fyy(iya,jxa,lrza,ngena),fxxyy(iya,jxa,lrza,ngena),
c     &ftemp(iya,jxa),ibd(4),wk(nwka)

c-----input
cSAP090808 delete ngena,jxa,lrza
c      integer iy_(lrza),jx,lrz,ngen
c      real*8 theta0,v0,y(iya,lrza),x(jxa),
c     & f(iya,jxa,lrza,ngena),rho,rovera(lrza)
      integer jx,lrz,ngen,iy_(lrz)
      real*8 theta0,v0,y(iya,lrz),x(jx),
     & f(iya,jx,lrz,ngen),rho,rovera(lrz)
c-----local
      integer i,ir,im,idm 
c-----output
      real*8 ff,dfdx,dfdy
c-----external
      real*8 terp2p,terp2p_Sm
     
c     Could speed up by eliminating repeated searches of grid
c     in terp2, for x and y.

c      write(*,*)'distr v0,theta0,rho',v0,theta0,rho
      if (v0.gt.1d0) then
c         write(*,*)'distr v0.gt.1 v0',v0
         ff=-50.d0*dlog(10.d0)
         dfdx=0.d0
         dfdy=0.d0
         return
      endif

      do i=1,lrz
         ir=i 
         if(rho.lt.rovera(ir)) goto 2  
      enddo
2     idm=iya
c      write(*,*)'distr ir,rho,rovera(i)',ir,rho,rovera(i)
cSm040110
c      idm=iy_(ir)      
c      ff=terp2p(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1      fxxyy(1,1,ir,ngen),idm,0,0,1)    
c      dfdx=terp2p(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1      fxxyy(1,1,ir,ngen),idm,1,0,1)
c      dfdy=terp2p(theta0,v0,iy_(ir),y(1,ir),jx,x,
c     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
c     1      fxxyy(1,1,ir,ngen),idm,0,1,1)

cSm040111 use spline with maximal length iya>iy_(ir) in y direction
c       idm=iy_(ir)
       idm=iya
      ff=terp2p_Sm(theta0,v0,iy_(ir),y(1,ir),jx,x,
     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
     1      fxxyy(1,1,ir,ngen),idm,0,0,1,iya)    
      dfdx=terp2p_Sm(theta0,v0,iy_(ir),y(1,ir),jx,x,
     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
     1      fxxyy(1,1,ir,ngen),idm,1,0,1,iya)
      dfdy=terp2p_Sm(theta0,v0,iy_(ir),y(1,ir),jx,x,
     1      f(1,1,ir,ngen),fxx(1,1,ir,ngen),fyy(1,1,ir,ngen),
     1      fxxyy(1,1,ir,ngen),idm,0,1,1,iya)

c      write(*,*)'distr 1 ff,dfdx,dfdy',ff,dfdx,dfdy
      if(rho.gt.rovera(lrz)) return
      if(ir.eq.1) return 

      im=ir-1
cSm040110
c      ffm=terp2p(theta0,v0,iy_(im),y(1,im),jx,x,
c     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
c     1      fxxyy(1,1,im,ngen),idm,0,0,1)
c      dfdxm=terp2p(theta0,v0,iy_(im),y(1,im),jx,x,
c     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
c     1      fxxyy(1,1,im,ngen),idm,1,0,1)
c      dfdym=terp2p(theta0,v0,iy_(im),y(1,im),jx,x,
c     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
c     1      fxxyy(1,1,im,ngen),idm,0,1,1)
cSm040111 use spline with maximal length iya>=iy_(ir) in y direction
c      idm=iy_(im)
       idm=iya 
      ffm=terp2p_Sm(theta0,v0,iy_(im),y(1,im),jx,x,
     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
     1      fxxyy(1,1,im,ngen),idm,0,0,1,iya)
      dfdxm=terp2p_Sm(theta0,v0,iy_(im),y(1,im),jx,x,
     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
     1      fxxyy(1,1,im,ngen),idm,1,0,1,iya)
      dfdym=terp2p_Sm(theta0,v0,iy_(im),y(1,im),jx,x,
     1      f(1,1,im,ngen),fxx(1,1,im,ngen),fyy(1,1,im,ngen),
     1      fxxyy(1,1,im,ngen),idm,0,1,1,iya)

c      write(*,*)'distr im,ffm,dfdxm,dfdym',im,ffm,dfdxm,dfdym

      d_rho=rho-rovera(im)
      d_rovera=rovera(ir)-rovera(im)
      a=d_rho/d_rovera
      ff=(ff-ffm)*a+ffm
      dfdx=(dfdx-dfdxm)*a+dfdxm
      dfdy=(dfdy-dfdym)*a+dfdym
c      write(*,*)'distr a,ff,dfdx,dfdy',a,ff,dfdx,dfdy
      return
      end



