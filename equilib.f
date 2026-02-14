
c        *********************EQUILIB************************
c        *  EQUILIB is the code for reading data for       *
c        *  toroidally symmetric plasma equilibrium given   *
c        *  the EQDSK (ref. Lang Lao) format.              *
c        *  It transforms these data for input data for    *
c        *  GENRAY code.				   *
c        *  RAYINP normalizes the eqdsk parameters by	   *
c        *  r0x(in m) and b0(in Tl) and			   *
c        *  creates the arrays for spline codes which      *
c        *  calculates: the (psi) and (feqd) magnetic field   *
c        *  functions; zlimit upper(r) and zlimit under(r) *
c        ***************************************************
c
c-------------------------------------------------------------------!
c        input data are read from:				    !
c        equilib.dat=       eqdsk  file    			    !
c        genray.in=        a namelist file with the normalization  !
c                           parameters r0x and b0                   !
c-------------------------------------------------------------------!
c********************************************************************
c  This program uses following external files:
c  input , dinitr , limitr 
c*********************************************************************
      subroutine equilib
      !implicit double precision (a-h,o-z)
      implicit none 

      call input

      call dinitr
      
      return
      end


c        ********************INPUT**************************
c        *   subroutine INPUT reads the eqdsk data , 	   *
c        *   normalizes the eqdsk parameters by 	   *
c        *   r0x(in m) and b0(in Tl)                       *
c        ***************************************************
c
c-------------------------------------------------------------------!
c        input data are read from:				    !
c        eqdsk_cq           file                 		    !
c        genray.in         it is a  file to read  the normalization!
c                           parameters r0x and b0                   !
c        output data  for subroutine limitr() are writen in:	    !
c        common 'three' and 'fourb'                                 !
c-------------------------------------------------------------------!
      subroutine input

      !implicit double precision (a-h,o-z)
      implicit none 

      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'three.i'
      include 'fourb.i'
      include 'one.i'
      include 'onetwo.i'

* nxeqd,nyeqd      are the numbers of points in the horizontal (r)
*                  and vertical (z) directions.
* xdimeqd,ydimeqd  are the horizontal and vertical full-widths
*                  of the  rectangle [meters]
* reqd             is the nominal major radius of the torus.
* redeqd           is the major radius of the inner edge
*                  of rectangular grid.
* ymideqd          is the vertical shift of the rectangular box
*                  up-down simmetry plane.
* xma,yma          are the major and vertical height of magnetic axis.
* psimag,psilim    are the poloidal flux function values at the
*                  magnetic axis and the last closed flux surface
*                   (touching the limiter or the separatrix).
* beqd             is the toroidal magnetic field at reqd.
* toteqd           is the toroidal current.
* psimx1,psixm2,xax1,xax2,zax1,zax2,psisep,xsep,ysep - OBSOLETE.
*  ATTENTION:
* psimx1=psimag
* xax1=xma
* zax1=zma
* psisep=psilim
c       feqd(nxeqda),pres(nxeqda)
* feqd = r * B_phi  at nxeqd equispaced points in psi
*                   from psimag to psilim.
* p                 are the pressure values at the same points.
c       ffpeqd(nxeqda),ppeqd(nxeqda)
* ffpeqd           (=  feqd_prime) at the same points.
* ppeqd            (=  p_prime) at the same points.
c       peqd(nxeqda,nyeqda)
* peqd(nxeqd,nyeqd)  are the psi values on the nxeqd * nyeqd
*                     equispaced grid.
* qpsi(nxeqda)        are q characteristic.
* nnlim,nnves        are the numbers of point at limiters and
*                    vacuum vessel wall.
c       rlimit(nlimit),zlimit(nlimit)
* rlimit,zlimit      is the r,z location of the limiters wall.
c       rves(nves),zves(nves)
* rves,zves          is the r,z location of the limiters
*                    vacuum vessel wall.


      character*12 namfil
      
c      namelist /genr/ r0x,b0,outdat,stat,mnemonic,rayop,
c     +partner,outnetcdf,dielectric_op
c      namelist /tokamak/ indexrho,ipsi,ionetwo,ieffic,psifactr,
c     +deltripl,nloop,i_ripple,eqdskin,NR
cSm070426
      include 'name_genr.i'
      include 'name_tokamak.i'
      integer kode,i,j
      integer io1,io2
      integer ierr !YuP[2020-03-09] for read(,,iostat=ierr)
      logical file_exists !YuP[2020-03-09] for INQUIRE (FILE=fname,EXIST=exists)
      real*8  dpsiar,dflux,pr,pb,dpsi,dstep

      external length_char
      integer length_char
cSAP090923     
c     Data from /genr/ and /tokamak/ namelists were read
c     in genray.f using read_all_namelists
      goto 10
c-----------------------------------------------------------
c     input /genr/ and /tokamak/ namelists
c-----------------------------------------------------------
cSAP100914
c       open(1,file='genray.in',status='old',iostat=kode)
      open(1,file='genray.dat',status='old',iostat=kode,
     +     delim='apostrophe')
      if(myrank.eq.0)then
         WRITE(*,*)'equiib.f after open(1,file=genray.dat kode',kode 
         if(kode.eq.0) WRITE(*,*)'======== USING equilib.dat =========='
      endif

      if (kode.ne.0) then
         open(1,file='genray.in',status='old',iostat=kode,
     +        delim='apostrophe')
         if (kode.ne.0) then
            if(myrank.eq.0)then
            WRITE(*,*)'dinit:Neither genray.in or genray.dat r present'
            endif
            stop ! stop in all cores (ranks)
         else
            if(myrank.eq.0)then
            WRITE(*,*)'======== USING equilib.in ==========='
            endif
         endif
      endif

cyup      write(*,*)'equilib.f before read genr'
      rewind(unit=1)
      read(1,genr,iostat=kode)
cyup      write(*,*)'equilib.f after read genr kode=',kode
      write(*,genr)
      call check_read(kode,'genr')

c      if (partner.ne.'disabled') write(*,*)'partner = ',partner
cyup      write(*,*)'partner = ',partner
c     Check partner values to see if valid (and not old-style):
      if (partner.ne.'disabled'
     1    .and. partner.ne.'genray_profs_in.txt'
     1    .and. partner.ne.'genray_profs_in.nc'
     1    .and.partner.ne.'genray_profs_out.nc') then
         write(*,*)'equilib.f:  Problem with value of partner'
         stop
      endif
 
      rewind(unit=1)
      read(1,tokamak,iostat=kode)
      write(*,*)'equilib.f tokamak kode=',kode
      write(*,tokamak)
      call check_read(kode,'tokamak')

c-----check the input dtata in namelist /tokamak/
cSm050304
      if ((indexrho.lt.1).or.(indexrho.gt.6)) then
         write(*,*)'in equilib.f in reading namelist /tokamak/'
         write(*,*)'It should be 0<indexrho<7, but indexrho =',indexrho
         write(*,*)'Change indexrho in genray.in file'
         stop
      endif

      if ((ipsi.lt.0).or.(ipsi.gt.1)) then
         write(*,*)'in equilib.f in reading namelist /tokamak/'
         write(*,*)'It should be -1<ipsi<2, but ipsi =',ipsi
         write(*,*)'Change ipsi in genray.in file'
         stop
      endif

      if ((ionetwo.lt.0).or.(ionetwo.gt.1)) then
         write(*,*)'in equilib.f in reading namelist /tokamak/'
         write(*,*)'It should be -1<ionetwo<2, but ionetwo =',ionetwo
         write(*,*)'Change ionetwo in genray.in file'
         stop
      endif

      if ((ieffic.lt.1).or.(ieffic.gt.6)) then
         write(*,*)'in equilib.f in reading namelist /tokamak/'
         write(*,*)'It should be 0<ieffic<7, but ieffic =',ieffic
         write(*,*)'Change ieffic in genray.in file'
         stop
      endif

      if ((psifactr.le.0).or.(psifactr.gt.1)) then
         write(*,*)'in equilib.f in reading namelitsst /tokamak/'
         write(*,*)'It should be 0<psifactr=<1, but psifactr=',psifactr
         write(*,*)'Change psifactr in genray.in file'
         stop
      endif

      write(*,*)'deltripl,nloop,i_ripple',deltripl,nloop,i_ripple
      if((i_ripple.ne.1).and.(i_ripple.ne.2))then
        write(*,*)'in dinit.f i_riplle should be equal 1 or 2'
        write(*,*)'i_ripple=',i_ripple
        write(*,*)'Please change i_ripple in /tokamak/'
        write(*,*)'in genray.in.file'
        stop
      endif

      if (length_char(trim(eqdskin)).gt.256) then
         WRITE(*,1001)
 1001    format('STOP: eqdskin spec too long')
         STOP
      endif
    
      if (NR.gt.NRA) then
        WRITE(*,*)'NR > NRA'
        WRITE(*,*)'it should be NR.le.NRA'
        WRITE(*,*)'please change NR in genray.dat or NRA in param.i' 
        STOP 
      endif

      close(1)

c-----end check of namelist /tokamak/

cyup      write(*,*)'in equilib r0x,b0',r0x,b0
cyup      write(*,*)'ipsi',ipsi

cyup      write(*,*) ' End of input'

cSAP090923
 10   continue
c-----------------------------------------------------------
c     Read the EQDSK file
c-----------------------------------------------------------
      write(*,*)'eqdskin= ',trim(eqdskin)
      
      !YuP[2020-03-09] Added message+stopping if no eqdsk file
      !INQUIRE about file's existence:
      INQUIRE(FILE=trim(eqdskin), EXIST=file_exists)
      write(*,*)'file_exists=',file_exists
      if (.NOT. file_exists) then
        WRITE(*,'(2A/)') 'equilib: Cannot find file ', trim(eqdskin)
        STOP
      endif

      open(30,file=trim(eqdskin),iostat=ierr)
cBH020822  Adding nveqd (.le. nxeqd) for different number of
cBH020822  flux surfaces on which p,feqd,p',feqd',q are tabulated.
cBH020822  The standard EQDSK does not incorporate this feature,
cBH020822  although is is available in cql3d.
cBH020822  Bonoli uses it for ACCOME eqdsk output.

c      read(30,2) nxeqd,nyeqd
      nveqd=0
cyup      write(*,*)'equilib.f before read (30,2)'
      read(30,2,iostat=io1) nxeqd,nyeqd,nveqd
      if (io1.ne.0) then
         rewind 30
         read(30,2,iostat=io2) nxeqd,nyeqd
         nveqd=0
      endif

cyup      write(*,*)'nxeqd,nyeqd,nveqd',nxeqd,nyeqd,nveqd
2     format(52x,3i4)

      nreqd= nxeqd  ! YuP: R-grid size
      nzeqd= nyeqd  ! YuP: Z-grid size
      if (nveqd.gt.nxeqd) stop 'nveqd.gt.nxeqd NOT ENABLED'
      if (nveqd.eq.0) nveqd=nxeqd

      call check_param(1)
c      if ((nxeqd.gt.nxeqda).or.(nyeqd.gt.nyeqda)) then
c        write(*,1000) nxeqd,nxeqda,nyeqd,nyeqda
c 1000   format('in equilib.dat in input',/,
c     .  'the dimensions of eqdsk (in eqilib.dat) nxeqd or nyeqd',/,
c     .  'are bigger than the parameters nxeqda or nyeqda in param.i'
c     .  ,/'nxeqd=',I5,'nxeqda=',I5
c     .  ,/'nyeqd=',I5,'nyeqda=',I5
c     .  ,/'Please change nxeqda or nyeqda in param.i')
c        stop
c      endif

      read(30,3) xdimeqd,ydimeqd,reqd,redeqd,ymideqd
cyup      write(*,*)'xdimeqd,ydimeqd,reqd,redeqd,ymideqd'
cyup      write(*,*)xdimeqd,ydimeqd,reqd,redeqd,ymideqd
3     format(5e16.9)

!YuP[2020-07-23] added, for border control
      rdimeqd=xdimeqd
      zdimeqd=ydimeqd
      zmideqd=ymideqd
      xeqmax= redeqd+rdimeqd ! YuP max for x-grid of cartesian grid
      xeqmin=-xeqmax
      yeqmax= redeqd+rdimeqd ! YuP max for y-grid of cartesian grid
      yeqmin=-yeqmax
      zeqmax= zmideqd + 0.5d0*zdimeqd
      zeqmin= zmideqd - 0.5d0*zdimeqd
      write(*,*)'xeqmin,xeqmax=',xeqmin,xeqmax
      write(*,*)'yeqmin,yeqmax=',yeqmin,yeqmax
      write(*,*)'zeqmin,zeqmax=',zeqmin,zeqmax

      read(30,3) xma,yma,psimag,psilim,beqd
      write(*,*)'equilib/reading magn.axis coord xma[m],yma[m]=',xma,yma
!yup      write(*,*)xma,yma,psimag,psilim,beqd
      if(psimag.gt.psilim) then
        dpsimax=-1
      else
        dpsimax=1
      endif
      psimag=dpsimax*psimag
      psilim=dpsimax*psilim
cyup      write(*,*)'in equilib psimag,psilim',psimag,psilim
      read(30,3) toteqd,psimx1,psimx2,xax1,xax2
      read(30,3) zax1,zax2,psisep,xsep,ysep
      read(30,3) (feqd(i),i=1,nveqd)
c********************************************
c      do i=1,nveqd
c        feqd(i)=beqd*xma
c         feqd(i)=feqd(i)*2.1d0/1.13d0
c      end do
c********************************************
      read(30,3) (pres(i),i=1,nveqd)
      read(30,3) (ffpeqd(i),i=1,nveqd)
      read(30,3) (ppeqd(i),i=1,nveqd)
      read(30,3) ((peqd(i,j),i=1,nxeqd),j=1,nyeqd)
c      write(*,3)((peqd(i,j),i=1,nxeqd),j=1,nyeqd)
c------------------------------------------------------------
c     creation of the poloidal flux peqd(i,j) with the minimum value
c     on the magnetic axis (psimag<psilim)
      do i=1,nxeqd
        do j=1,nyeqd
	   peqd(i,j)=dpsimax*peqd(i,j)
	enddo
      enddo
cc      write(*,*)'in equilib after change of the peqd sign'
c------------------------------------------------------------
      read(30,3) (qpsi(i),i=1,nveqd)
      write(*,*)'in equilib after reading qpsi'
      
      read(30,4,iostat=ierr) nnlim,nnves !YuP[2020-03-09] Why was it skipped?
      if(ierr.eq.0)then !YuP[2020-03-09]
        continue
      else ! lines with 'nnlim,nnves' are not present in file
        nnlim=0
        nnves=0
        !Note: Declared rlimit(nlimit),zlimit(nlimit); nlimit in param.i
        !Note: Declared rves(nves),zves(nves), where nves is in param.i
      endif !YuP[2020-03-09]
      write(*,*)'in equilib after reading nnlim,nnves=',nnlim,nnves
      write(*,*)'(if nnlim=0 and nnves=0 then - no data in eqdsk)'
 4    format(2i5)
     
c-----------------------------
!YuP[2020-03-09] Why skipping?  Not skipping - affects results of a run.
      nnlim=0 !=0 means skip the data from eqdskin ! affects results !
      !nnves=0 !=0 means skipping data from eqdskin ! no effect on results
c-----------------------------
      if (nnlim.ne.0) then 
         if(nlimit.lt.nnlim) then !YuP[2020-03-09] adjusted
           write(*,*)'nlimit=',nlimit,' is smaller than nnlim=',nnlim
           write(*,*)' see  common five.i and fourb.i'
           WRITE(*,*)'equilib.f: increase nlimit in param.i'
           stop
         endif
         read(30,3) (rlimit(i),zlimit(i),i=1,nnlim)
      else
         !eqdsk data w/o limiter points, or skipped because of nnlim=0 above
         write(*,*)'nnlim=0: eqdsk data w/o limiter points, or skipped'
         psisep=psilim
         write(*,*)'psisep,psilim',psisep,psilim
      endif
      
      if (nnves.ne.0) then
         if(nves.lt.nnves) then !YuP[2020-03-09] adjusted
           write(*,*)'nves=',nves,' is smaller than nnves=',nnves
           write(*,*)' see  common five.i and fourb.i'
           WRITE(*,*)'equilib.f: increase nves in param.i'
           stop
         endif
         read(30,3) (rves(i),zves(i),i=1,nnves)
      else
         write(*,*)'nnves=0 : eqdsk data without vessel points'
      endif

      close(30) !file=eqdskin ---------------------------------
      
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)' nxeqd,nyeqd,nveqd'
      write(*,2) nxeqd,nyeqd,nveqd
      !pause !!!
cc      write(*,*)'xdimeqd,ydimeqd,reqd,redeqd,ymideqd'
cc      write(*,3) xdimeqd,ydimeqd,reqd,redeqd,ymideqd
cc      write(*,*)'xma,yma,psimag,psilim,beqd'
cc      write(*,3) xma,yma,psimag,psilim,beqd
cc      write(*,*)'in equilib psimag,psilim',psimag,psilim
cc      write(*,*)'toteqd,psimx1,psixm2,xax1,xax2'
cc      write(*,3) toteqd,psimx1,psixm2,xax1,xax2
cc      write(*,*)'zax1,zax2,psisep,xsep,ysep'
cc      write(*,3) zax1,zax2,psisep,xsep,ysep
      write(*,*) 'feqd(i),i=1,nveqd'
      write(*,3) (feqd(i),i=1,nveqd)
      write(*,*) 'pres(i),i=1,nveqd'
      write(*,3) (pres(i),i=1,nveqd)
      write(*,*) 'ffpeqd(i),i=1,nveqd'
      write(*,3) (ffpeqd(i),i=1,nveqd)
      write(*,*) 'ppeqd(i),i=1,nveqd'
      write(*,3) (ppeqd(i),i=1,nveqd)
      !pause !!!
c      do i=1,nxeqd
c         do j=1,nyeqd
c            write(*,*)'i,j,peqd(i,j)',i,j,peqd(i,j)
c         enddo
c      enddo
c      write(*,*) '(peqd(i,j),i=1,nxeqd),j=1,nyeqd)'
cc      write(*,3) ((peqd(i,j),i=1,nxeqd),j=1,nyeqd)
      write(*,*) 'qpsi(i),i=1,nveqd'
      write(*,3) (qpsi(i),i=1,nveqd)
      write(*,*)'nnlim,nnves'
      write(*,4) nnlim,nnves
cc      write(*,*) 'rlimit(i),zlimit(i),i=1,nnlim'
cc      write(*,3) (rlimit(i),zlimit(i),i=1,nnlim)
cc      write(*,*) 'rves(i),zves(i),i=1,nnves'
cc      write(*,3) (rves(i),zves(i),i=1,nnves)
      endif ! outprint




c-----------------------------------------------------------
c-----------------------------------------------------------

      if (nveqd.ne.nxeqd) then

c-----------------------------------------------------------
c     Interpolate feqd,pres,ffpeqd,ppeqd,qpsi to 
c     nxeqd equispaced points.
c     The nveqd option has been enabled in CQL3D, and we
c     enable it also for GENRAY.  It is not a part of the
c     standard EQDSK  (BobH, 020722).
c
c     We use the standard, well documented spline package,
c     in zcunix.f, rather than use the uncommented subroutines
c     used in the majority of GENRAY cases (BobH, 020722).
c-----------------------------------------------------------

cAdding new arrays: psiar,d2feqd,d2pres,d2ffpeqd,d2ppeqd,d2qpsi,
c                   workk,r8temp,tabl,i1p,itabl
c     Create cubic spline arrays
c
         dpsiar=(psilim-psimag)/(nveqd-1)
         do i=1,nveqd
            psiar(i)=psimag+(i-1)*dpsiar
         enddo
         write(*,*) 'psiar(i),i=1,nveqd'
         write(*,3) (psiar(i),i=1,nveqd)
c     Use cubic splines at end points to obtain coefficients:
         i1p(1)=4
         i1p(2)=4
         call coeff1(nveqd,psiar,feqd,d2feqd,i1p,1,workk)
         call coeff1(nveqd,psiar,pres,d2pres,i1p,1,workk)
         call coeff1(nveqd,psiar,ffpeqd,d2ffpeqd,i1p,1,workk)
         call coeff1(nveqd,psiar,ppeqd,d2ppeqd,i1p,1,workk)
         call coeff1(nveqd,psiar,qpsi,d2qpsi,i1p,1,workk)
         
c     Create nxeqd long array and spline onto it.
         dflux=(psilim-psimag)/(nxeqd-1)
         do i=1,nxeqd
            flux(i)=psimag+(i-1)*dflux
         enddo
         
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0
         do i=1,nxeqd
            call terp1(nveqd,psiar,feqd,d2feqd,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nxeqd
            feqd(i)=r8temp(i)
         enddo
         
         do i=1,nxeqd
            call terp1(nveqd,psiar,pres,d2pres,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nxeqd
            pres(i)=r8temp(i)
         enddo
         
         do i=1,nxeqd
            call terp1(nveqd,psiar,ffpeqd,d2ffpeqd,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nxeqd
            ffpeqd(i)=r8temp(i)
         enddo
         
         do i=1,nxeqd
            call terp1(nveqd,psiar,ppeqd,d2ppeqd,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nxeqd
            ppeqd(i)=r8temp(i)
         enddo
         
         do i=1,nxeqd
            call terp1(nveqd,psiar,qpsi,d2qpsi,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nxeqd
            qpsi(i)=r8temp(i)
         enddo
         
         write(*,*) 'flux(i),i=1,nxeqd'
         write(*,3) (flux(i),i=1,nxeqd)
         write(*,*) 'feqd(i),i=1,nxeqd'
         write(*,3) (feqd(i),i=1,nxeqd)
         write(*,*) 'pres(i),i=1,nxeqd'
         write(*,3) (pres(i),i=1,nxeqd)
         write(*,*) 'ffpeqd(i),i=1,nxeqd'
         write(*,3) (ffpeqd(i),i=1,nxeqd)
         write(*,*) 'ppeqd(i),i=1,nxeqd'
         write(*,3) (ppeqd(i),i=1,nxeqd)
         write(*,*) 'qpsi(i),i=1,nxeqd'
         write(*,3) (qpsi(i),i=1,nxeqd)

         
c-----------------------------------------------------------
c-----------------------------------------------------------
         
      endif
      

c---------------------------------------------------------
c  normalization of the eqdsk data
c--------------------------------------------------------
      pr=1.d0/r0x
      pb=1.d0/b0
      xdimeqd=xdimeqd*pr
      ydimeqd=ydimeqd*pr
      reqd=reqd*pr
      redeqd=redeqd*pr
      ymideqd=ymideqd*pr
      xma=xma*pr
      yma=yma*pr
      xeqmin=xeqmin*pr
      xeqmax=xeqmax*pr
      yeqmin=yeqmin*pr
      yeqmax=yeqmax*pr
      zeqmin=zeqmin*pr
      zeqmax=zeqmax*pr
      psimag=psimag*pr*pr*pb
      psilim=psilim*pr*pr*pb
      write(*,*)'equilib psimag,psilim,pr,pb',psimag,psilim,pr,pb

      beqd=beqd*pb
      psimx1=psimx1*pr*pr*pb
      psimx2=psimx2*pr*pr*pb
      xax1=xax1*pb
      xax2=xax2*pb
      zax1=zax1*pr
      zax2=zax2*pr
      psisep=psisep*pr*pr*pb
      xsep=xsep*pr
      ysep=ysep*pr
      write(*,*)'in equilib psimag,psilim',psimag,psilim
      write(*,*)'nveqd,nxeqda',nveqd,nxeqda
      
c----- GRIDS/extra ---------------------------------------------
      dstep= rdimeqd/(nreqd-1) !same as xdimeqd/(nxeqd-1)
      write(*,*) 'rdimeqd, redeqd=', rdimeqd, redeqd
      do i=1,nreqd
        req(i)= redeqd+dstep*(i-1) !same as rr(i)
      enddo
 
      dstep= (zeqmax-zeqmin)/(nzeqd-1)
      !where zeqmin= zmideqd - 0.5d0*zdimeqd
      !and   zeqmax= zmideqd + 0.5d0*zdimeqd
      do i=1,nzeqd
        zeq(i)= zeqmin+dstep*(i-1) ! YuP [zeqmin; zeqmax]
      enddo

      !Only for subroutine wrtnetcdf_plasma_prof
      dstep= (xeqmax-xeqmin)/(2*nxeqd-1)
      !where xeqmax= redeqd+rdimeqd
      !and   xeqmin=-xeqmax
      do i=1,2*nxeqd
        xeq(i)= xeqmin+dstep*(i-1) ! YuP [xeqmin; xeqmax] X-axis in top view
      enddo
      yeq(:)=xeq(:) !Y-axis in top view
c----- GRIDS/extra (done)------------------------------------------

      do 30 i=1,nveqd 
        feqd(i)=feqd(i)*pr*pb
	pres(i)=pres(i)*pr*pr*pb
 30   continue

      do 31 i=1,nxeqd
        do 32 j=1, nyeqd
 32       peqd(i,j)=peqd(i,j)*pr*pr*pb
 31   continue
      if (nnlim.gt.0) then
        do 33 i=1,nnlim
          rlimit(i)=rlimit(i)*pr
 33	  zlimit(i)=zlimit(i)*pr
      endif
c------------------------------------------------------------------
cBH020822      dpsi=(psilim-psimag)/(nxmax-1)
      dpsi=(psilim-psimag)/(nxeqd-1)
      do 15 i=1,nxeqd
        flux(i)=psimag+dpsi*(i-1)
15     continue

      dstep=xdimeqd/(nxeqd-1)
      do 100 i=1,nxeqd
        rr(i)=redeqd+dstep*(i-1)
 100  continue
      dstep=ydimeqd/(nyeqd-1)
      do 200 i=1,nyeqd
        zz(i)=ymideqd+dstep*(i-1)-ydimeqd*0.5d0
 200  continue

c-------------------------------------------------------------
c     test map psi(i,j) on (R,Z) plane
c------------------------------------
c      open(40,file='psirz.dat')
c41     format(3(' ',e16.9))
c      write(40,*)' r z psi'         
c      do i=1,nxeqd
c         do j=1,nyeqd
c         write(40,41)rr(i),zz(j) ,peqd(i,j)
c         enddo
c      enddo
c      close(40)
c      stop
c     end test map
c---------------------------------------------------

      return
      end


c*************************DINITR***************************************
c  creation of arrays for eqdsk spline coefficients		      *
c  for psieqd:tx(nx4a),ty(ny4a),cxy(nx4a,ny4a)			      *
c  for feqd     :txf(nx4a),cx(nx4a)				      *
c  for pres  :tpres(nx4a),cpres(nx4a)	        		      *
c  creation of arrays for limiter spline coefficients		      *
c  for zlimmiter plus(r)  :trlimp(nlim4),cxlimp(nlim4)		      *
c  for zlimmiter minus(r) :trlimm(nlim4),cxlimm(nlim4)		      *
c								      *
c  these program uses the following subroutines:  		      *
c								      *
c  limitr(zzp,rrp,zzm,rrm,ip,im,				      *
c         rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin)		      *
c  for creation arrays:{zzp(ip),rrp(ip)} for zlimmiter upper (r),     *
c                      {zzm(im),rrm(im)} for zlimmiter under (r)      *
c  from zlimit,rlimit						      *
c								      *
c  iac1r(rlimr,ip,ip4,zlimr,lx,mxa,fxa,mxb,fxb,trlimp,cxlimp)	      *
c  -calculates the spline coefficients for 1d function		      *
c								      *
c  IAC2R(rrr,nx,zzr,ny,peqdr,nfx,nfy,lx,ly,mxa,arfxa,mxb,arfxb,	      *
c  mya,arfya,myb,arfyb,tx,ty,cxy,ncx,ncy,nry,cy)		      *
c  - calculates the spline coefficients for 2d function		      *
c**********************************************************************

      subroutine dinitr
      implicit double precision (a-h,o-z)
      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'three.i'
      include 'fourb.i'
      include 'five.i'
      include 'gr.i'
cSAP091031 
c     include 'limiter.i'

      dimension zzp(nlimit),rrp(nlimit),zzm(nlimit),rrm(nlimit)
      double precision
     1     arfxa(nyeqda),arfxb(nyeqda),arfya(nxeqda),arfyb(nxeqda),
     1	   fxa,fxb,fr(nxeqda),fluxr(nxeqda),rrr(nxeqda),zzr(nyeqda),
     2	   peqdr(nxeqda,nyeqda),
     3	   zlimr(nlimit),rlimr(nlimit)
      double precision	ias1r,ias2r,ias2r_Sm
      character*12 namfil

c--------------------------------------------------------------------
c     determination of rmax,rmin,zmax,zmin
c--------------------------------------------------------------------
      write(*,*)'in dinitr nlimit=',nlimit

      if (nnlim.gt.0)  then

         call limitr(zzp,rrp,zzm,rrm,ip,im,
     1   rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin)
c         write(*,*)'in dinitr after call limitr
c     1   ip,im,rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin',
c     1   ip,im,rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin
c      else
c         call makelacn(zmax1,rmax1,zmin1,rmin1,
c     1   rzmax1,zrmax1,rzmin1,zrmin1)
c         rmax=rmax1
c         zmax=zmax1
c         rmin=rmin1
c         zmin=zmin1
c         zrmax=zrmax1
c         rzmax=rzmax1
c         zrmin=zrmin1
c         rzmin=rzmin1

      endif
c--------------------------------------------------------------------
c     spline coefficients for psi, feqd and pres creation 
c--------------------------------------------------------------------
cSm030224
      nx=nxeqd
      ny=nyeqd
      nx4=nx+4
      ny4=ny+4

      lx=1
      mxa=0
      fxa=0.d0
      mxb=1
      fxb=0.d0
      do 1111 i=1,nx
        fluxr(i)=flux(i)
	fr(i)=feqd(i)
 1111	rrr(i)=rr(i)
      do 1112 j=1,ny
 	zzr(j)=zz(j)
 1112   continue

      do 1113 i=1,nx
      do 1113 j=1,ny
 1113   peqdr(i,j)=peqd(i,j)
      			    
      call iac1r(fluxr,nx,nx4,fr,lx,mxa,fxa,mxb,fxb,txf,cx)
c-------------------------------------------------------------------	
      call iac1r(fluxr,nx,nx4,pres,lx,mxa,fxa,mxb,fxb,tpres,cpres)
c test pres beg
c      do i=1,nx 
c         psi=flux(i)
c         pressure=prespsi(psi)
c	 write(*,*)'in equilib i,psi,pres(i),pressure'
c	 write(*,*) i,psi,pres(i),pressure
c      enddo !test
c      stop
c test pres end
c--------------------------------------------------------------------
      ny=nyeqd
      lx=1
      ly=1
      mxa=0
      do 91 i=1,ny
 91      arfxa(i)=0.
      mxb=0
      do 92 i=1,ny
 92     arfxb(i)=0.
      mya=0
      do 93 i=1,nx
 93     arfya(i)=0.
      myb=0
      do 94 i=1,nx
 94     arfyb(i)=0.
      nfx=nx
      nfy=ny
      ncx=nx+4
      ncy=ny+4

cSm030224
      nry=max(nx,ny)+4
c      call iac2r(rrr,nx,zzr,ny,peqdr,nfx,nfy,lx,ly,mxa,arfxa,mxb,arfxb,
c     1                 mya,arfya,myb,arfyb,tx,ty,cxy,ncx,ncy,nry,cy)
           
      call iac2r_Sm(rrr,nx,zzr,ny,peqdr,nfx,nfy,lx,ly,mxa,arfxa,mxb,
     &                 arfxb,
     &                 mya,arfya,myb,arfyb,tx,ty,cxy,ncx,ncy,nry,cy,
     &                 nxeqda,nx4a)
      write(*,*)
     1     'spline coefficients for psi, feqd, and pres were created'


ctest     
      cnorm=0.d0     
      do i=1,nx
         do j=1,ny
c           psi_t=fpsi(rrr(i),zzr(j))
           psi_t=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,rrr(i),zzr(j),
     &     nx4a)
c           write(*,*)'i,j,rrr(i),zzr(j),peqdr(i,j),psi_t',
c     &     i,j,rrr(i),zzr(j),peqdr(i,j),psi_t
           cnorm=cnorm+(peqdr(i,j)-psi_t)**2
         enddo
      enddo
      cnorm=dsqrt(cnorm)/dfloat(nx*ny)
      write(*,*)'equilib cnorm of psi',cnorm
cendtes     
c--------------------------------------------------------------------
c  The calculations of the contours coordinates for flux functions
c          r(psi,teta)   z(psi,teta)
c          arrays rpsi(j,i) zpsi(j,i)
c          j=1,npsi(number of counturs poloidal flux=constant)
c          i=1,nteta+1(number of point in poloidal angle)
c  and creates arrays ar(nl,nteta+1),az(nl,nteta+1) for nl contours
c  plotting

cSAP091030
      write(*,*)'dinitr before gr2new psilim',psilim,'nnlim',nnlim
cSAP091214
c      if (nnlim.eq.0) then
cc------- calculate limitter


cc         call calc_limiter_field_line(psilim,nlimit,
cc     &   yma,xma,rlimit,zlimit,theta_pol_limit,rho_geom_limit)

c          call calc_limiter_binary(psilim,nlimit,
c     &    yma,xma,rlimit,zlimit,theta_pol_limit,rho_geom_limit)

c         call limitr(zzp,rrp,zzm,rrm,ip,im,
c     1   rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin)
         
c      endif

      call gr2new
cSAP091030
      write(*,*)'dinitr after gr2new psilim',psilim
c--------------------------------------------------------------------
c     creation of spline coefficients for limiter
c--------------------------------------------------------------------
      if (nnlim.eq.0) then
c        ------------------------------------------------
c        creation of arrays zlimit(nlimit), rlimit(nlimit)
c        They set equal to rpsi(npsi,i) zpsi(npsi,i) (the
c        coordinates of the magnetic surface psi(r,s)=psilim.
c        In this case it must be nlimit=nteta1
c      	 -------------------------------------------------

cSAP091208
         if (nlimit.ne.nteta1) then
            write(*,*)'equilib.f in dinitr at nnlim=0'
            write(*,*)'nlimit.ne.nteta1'
            write(*,*)'nlimit,nteta,nteta1',nlimit,nteta,nteta1
            write(*,*)'it should be nlimit=nteta1'
            write(*,*)'Please change nlimit or nteta in param.i'
            write(*,*)'and recompile the code'
            stop 'in equilib.f in dinitr nlimit.ne.nteta1'
         endif

	 do i=1,nlimit
cyup           write(*,*)'npsi,i,nteta1', npsi,i,nteta1
	   zlimit(i)=zpsi(npsi,i)
	   rlimit(i)=rpsi(npsi,i)
cyup	   write(*,*)'nlimit,npsi',nlimit,npsi
cyup	   write(*,*)'i,zlimit(i),rlimit(i)',i,zlimit(i),rlimit(i)
	 enddo
cSAP091030
cyup         write(*,*)'dinitr before limitr psilim',psilim

         call limitr(zzp,rrp,zzm,rrm,ip,im,
     1   rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin)
c         write(*,*)'in dinitr  after call limitr limitr
c     1   ip,im,rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin',
c     1   ip,im,rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin

cSAP091030
cyup         write(*,*)'dinitr after limitr psilim',psilim

      endif
      ip4=ip+4
      lx=1
      mxa=0
      fxa=0.0d0
      mxb=0
      fxb=0.0d0
      if (rrp(1).lt.rrp(ip)) then
        do 111 i=1,ip
          zlimr(i)=zzp(i)
 111      rlimr(i)=rrp(i)
      else
        do 2112 i=1,ip
          zlimr(i)=zzp(ip-i+1)
 2112     rlimr(i)=rrp(ip-i+1)
      end if
 17   format(2e16.9)

cyup      write(*,*)'in equilib.for rlimr,zlimr,i=1,ip'
c      read(*,*)
cyup      write(*,17)(rlimr(i),zlimr(i),i=1,ip)

      call iac1r(rlimr,ip,ip4,zlimr,lx,mxa,fxa,mxb,fxb,trlimp,cxlimp)
c---------------------------------------------------------------------
c test begin
c---------------------------------------------------------------------
      do i=1,ip-1
        idx=0
        ipx=ip
        ipx4=ip+4
	 r=rlimr(i)
	 zrrr=zlimr(i)
         zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,r)
	 zp=zzrp
         call zpzmlim(r,zplin,zmlin)
c	 write(*,*)'ip,i,rlimr(i)',ip,i,rlimr(i)
c	 write(*,*)'zrr(array),zp(spline),zplin',zrrr,zp,zplin
	 r=0.5d0*(rlimr(i)+rlimr(i+1))
         zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,r)
	 zp=zzrp
         call zpzmlim(r,zplin,zmlin)
c	 write(*,*)'r,zp(spline),zplin',r,zp,zplin
       enddo
c---------------------------------------------------------------------
c test end
c---------------------------------------------------------------------
      im4=im+4
      if (rrm(1).lt.rrm(im)) then
        do 112 i=1,im
          zlimr(i)=zzm(i)
 112      rlimr(i)=rrm(i)
      else
        do 1122 i=1,im
          zlimr(i)=zzm(im-i+1)
 1122     rlimr(i)=rrm(im-i+1)
      end if

cyup      write(*,*)'in equilib.for rlimr,zlimr,i=1,im'
c      read(*,*)
cyup      write(*,17)(rlimr(i),zlimr(i),i=1,im)

      call iac1r(rlimr,im,im4,zlimr,lx,mxa,fxa,mxb,fxb,trlimm,cxlimm)
c---------------------------------------------------------------------
c test begin
c---------------------------------------------------------------------
c      write(*,*)'in equilib test zmin'
      do i=1,im-1
	 r=rlimr(i)
	 zrr=zlimr(i)
         idx=0
         ipx=im
         ipx4=im+4
         zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,r)
         zm=zzrm
         call zpzmlim(r,zplin,zmlin)
c 	 read(*,*)
c	 write(*,*)'im,i,rlimr(i)',im,i,rlimr(i)
c	 write(*,*)'zrr(array),zm(spline),zmlin',zrr,zm,zmlin
	 r=0.5d0*(rlimr(i)+rlimr(i+1))
         zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,r)
	 zm=zzrm
         call zpzmlim(r,zplin,zmlin)
c	 read(*,*)
c	 write(*,*)'im,i,rlimr(i)',im,i,rlimr(i)
c	 write(*,*)'r,zm(spline),zmlin',r,zm,zmlin
      enddo
c      stop
c---------------------------------------------------------------------
c test end
c---------------------------------------------------------------------

        return
        end


c*************LIMITR**************************************************
c   this subroutine creates the arrays {zzp(ip),rrp(ip)},
c   {zzm(im),rrp(im)} for function zlimit upper(r) and zlimit under(r)
c   and determinates the coordinates of:
c      top limitter point  (zmax,rzmax),
c      botom limitter point(zmin,rzmin),
c      inner limitter point(zrmin,rmin),
c      outer limitter point(zrmax,rmax)
c**********************************************************************
c   input data are in common 'three' and 'fourb'
c   output data are all parameters of the limitr()
c----------------------------------------------------------------------
      subroutine limitr(zzp,rrp,zzm,rrm,ip,im,
     1 rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin)
      implicit double precision (a-h,o-z)
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'three.i'
      include 'fourb.i'
      include 'limit.i'
      dimension zzp(nlimit),rrp(nlimit),zzm(nlimit),rrm(nlimit)
      double precision zlimr(nlimit),rlimr(nlimit)
c------------------------------------------------------------------
c   determination rmin and rmax  of limitter and
c   coordinates of these points: (zrmax,rmax),(zrmin,rmin)
c------------------------------------------------------------------
cSmirnov961226 test beg
cyup      write(*,*)'in limitr nlimit',nlimit
cyup      do i=1,nlimit
cyup        write(*,*)'i,zlimit(i),rlimit(i)',i,zlimit(i),rlimit(i)
cyup      enddo
cSmirnov961226 test end

      irmax=1
      irmin=1
      rmax=rlimit(1)
      rmin=rlimit(1)
      do 10 i=2,nlimit

c        write(*,*)'in limitr i rmax,rlimit(i)',i,rmax,rlimit(i)

        if (rmax.lt.rlimit(i)) then
          rmax=rlimit(i)
          irmax=i
          goto 10
        end if

c        write(*,*)'in limitr i rmin,rlimit(i)',i,rmin,rlimit(i)

        if (rmin.gt.rlimit(i)) then
          rmin=rlimit(i)
          irmin=i
        end if
 10   continue

      zrmin=zlimit(irmin)
      zrmax=zlimit(irmax)
cyup        write(*,*)'irmin=',irmin,'rmin=',rmin,'zrmin=',zrmin
cyup        write(*,*)'irmax=',irmax,'rmax=',rmax,'zrmax=',zrmax
c---------------------------------------------------------------------
c  begin of arrays rrp,zzp  and  rrm,zzm  creation
c---------------------------------------------------------------------
      if (irmin.gt.irmax) then
	 do 20 i=1,nlimit-irmin
	     rrp(i)=rlimit(irmin+i-1)
 	     zzp(i)=zlimit(irmin+i-1)
 20      continue

	 do 21 i=nlimit-irmin+1,nlimit-1
	     rrp(i)=rlimit(i-nlimit+irmin)
 	     zzp(i)=zlimit(i-nlimit+irmin)
 21	 continue
      else
c--------- if irmin.le.irmax then--------------------------
	 do 22 i=1,nlimit-irmin
	     rrp(i)=rlimit(i+irmin-1)
 	     zzp(i)=zlimit(i+irmin-1)
 22      continue
	 do 23 i=nlimit-irmin+1,nlimit-1
	     rrp(i)=rlimit(i-nlimit+irmin)
 	     zzp(i)=zlimit(i-nlimit+irmin)
 23      continue
      end if
c------------------------------------------------------------------
      do 24 i=1, nlimit-1
	 rlimit(i)=rrp(i)
         zlimit(i)=zzp(i)
cSAP091208
c         write(*,*)'equilib.f renumerated i,rlimit(i),zlimit(i)',
c     &                                    i,rlimit(i),zlimit(i)
 24   continue
c------------------------------------------------------------------
c     determination rmin and rmax ,zmin and zmax of limitter and
c     coordinates of these points:(zrmax,rmax),(zrmin,rmin),
c     (zmax,rzmax),(zmin,rzmin)
c------------------------------------------------------------------
      irmax=1
      irmin=1
      izmax=1
      izmin=1
      rmax=rlimit(1)
      rmin=rlimit(1)
      zmax=zlimit(1)
      zmin=zlimit(1)


      do 30 i=1,nlimit-1

c          write(*,*)'30 i,rmin,rmax,rlimit(i)',i,rmin,rmax,rlimit(i)
c          write(*,*)'zmin,zmax,zlimit(i)',zmin,zmax,zlimit(i)
c          write(*,*)'irmin,irmax,izmin,izmax',irmin,irmax,izmin,izmax
cSAP091208 for test the difference irmax with gr090903_091030
c	  if (rmax.lt.rlimit(i)) then
	  if (rmax.lt.rlimit(i)-1.d-15) then
	    rmax=rlimit(i)
	    irmax=i
	    zrmax=zlimit(irmax)
	  end if

	  if (rmin.gt.rlimit(i)) then
	    rmin=rlimit(i)
	    irmin=i
	    zrmin=zlimit(irmin)
	  end if

	  if (zmax.lt.zlimit(i)) then
	    zmax=zlimit(i)
	    izmax=i
	    rzmax=rlimit(izmax)
	  end if

	  if (zmin.gt.zlimit(i)) then
	    zmin=zlimit(i)
	    izmin=i
	    rzmin=rlimit(izmin)
	  end if

c          write(*,*)'30 1 i,rmin,rmax,rlimit(i)',i,rmin,rmax,rlimit(i)
c          write(*,*)'zmin,zmax,zlimit(i)',zmin,zmax,zlimit(i)
c          write(*,*)'irmin,irmax,izmin,izmax',irmin,irmax,izmin,izmax


 30   continue

c      write(*,*)'equilib,f limitr'
c      write(*,*)'irmin=',irmin,'rmin=',rmin,'zrmin=',zrmin
c      write(*,*)'irmax=',irmax,'rmax=',rmax,'zrmax=',zrmax
c      write(*,*)'izmin=',izmin,'zmin=',zmin,'rzmin=',rzmin
c      write(*,*)'izmax=',izmax,'zmax=',zmax,'rzmax=',rzmax
c------------------------------------------------------------------
c     rrp(1) must be less then rrp(ip)
c     rrm(1) must be less then rrm(im)
c     following part of program creates arrays {rrp,zzp} and
c     {rrm,zzm}  in  which rrp(1).lt.rrp(ip),
c                          rrm(1).lt.rrm(im)
c------------------------------------------------------------------
      if ((izmax.ge.1).and.(izmax.le.irmax)) then
c-----------------------------------------------------------------
c rrp(irmax),zzp(irmax),rrm(nlimit-irmax+1),zzm(nlimit-irmax+1)
c -----------------------------------------------------------------

	   ip=irmax
	   im=nlimit-irmax+1
	   do 40 i=1,ip
	     rrp(i)=rlimit(i)
  40	     zzp(i)=zlimit(i)
           do 41 i=1,im-1
	     rrm(i)=rlimit(i+irmax-1)
  41	     zzm(i)=zlimit(i+irmax-1)
             rrm(im)=rrp(1)
	     zzm(im)=zzp(1)
      else
c
c          rrm(irmax),zzm(irmax),rrp(nlimit-irmax+1),zzp(nlimit-irmax+1)
c
	   im=irmax
	   ip=nlimit-irmax+1
	   do 42 i=1,im
	     rrm(i)=rlimit(i)
  42	     zzm(i)=zlimit(i)
           do 43 i=1,ip-1
	     rrp(i)=rlimit(i+irmax-1)
  43	     zzp(i)=zlimit(i+irmax-1)
             rrp(ip)=rrm(1)
	     zzp(ip)=zzm(1)
      end if
c-----------------------------------------------------------------
      if (rrp(1).gt.rrp(ip)) then
        do i=1,ip
          zlimr(i)=zzp(i)
          rlimr(i)=rrp(i)
	enddo
        do i=1,ip
          zzp(i)=zlimr(ip-i+1)
          rrp(i)=rlimr(ip-i+1)
	enddo
      end if

      if (rrm(1).gt.rrm(im)) then
        do i=1,im
          zlimr(i)=zzm(i)
          rlimr(i)=rrm(i)
	enddo
        do i=1,im
          zzm(i)=zlimr(im-i+1)
          rrm(i)=rlimr(im-i+1)
	enddo
      end if
c-----------------------------------------------------------------
c  creation of the monotonic limiter array zzp  near rmin and rmax
      izmax=1
      zmax=zzp(1)
      do i=2,ip
        if (zmax.lt.zzp(i)) then
          zmax=zzp(i)
          izmax=i
	endif
      enddo

c      write(*,*)'ip,izmax,zzp(izmax),rrp(izmax)',
c     1 ip,izmax,zzp(izmax),rrp(izmax)

c--------near rmin --------------------------------
c     the correction for the second point (begin)
c      write(*,*)'before correction rrp(2),rrp(1)',rrp(2),rrp(1)
c      write(*,*)'before correction zzp(2),zzp(1)',zzp(2),zzp(1)

c      write(*,*)'equlib.f in limitr izmax',izmax
cSAP091102
      
      if (rrp(2).lt.rrp(1)+1.d-5) then

cSAP091102
         j0=2

         do j=3,izmax
c            write(*,*)'equlib.f in limitr j',j
            if (rrp(j).gt.rrp(2)) then
c               write(*,*)'equlib.f in limitr 1 j',j
	       j0=j
	       goto 45
	    endif
	 enddo
 45      continue
c         write(*,*)'equilib.f in limitr nlimit j0',j0
         if(zzp(j0).eq.zzp(1)) then
	    rrp(2)=rrp(1)+1.d-5
	 else
cSmirnov970104 beg
	    if(dabs(rrp(1)-rrp(2)).lt.1.d-6) then
	       rrp(2)=0.5d0*(rrp(j0)+rrp(1))
	    else
	    rrp(2)=rrp(1)+(zzp(2)-zzp(1))*(rrp(j0)-rrp(1))/
     1	       	      (zzp(j0)-zzp(1))
	    endif
cSmirnov961226 end
	 endif
      endif
c      write(*,*)'after correction rrp(2),rrp(1)',rrp(2),rrp(1)
c      write(*,*)'after correction zzp(2),zzp(1)',zzp(2),zzp(1)
c     the correction for the second point (end)
c     --------------------------------------
      do i=3,izmax-1
c         if(rrp(i).le.rrp(i-1)) then
         if(rrp(i).le.(rrp(i-1)+5.d-3)) then
	   j0=izmax
	   do j=i+1,izmax
	     if (rrp(j).gt.rrp(i-1)) then
	       j0=j
	       goto 50
	     endif
	   enddo
 50        continue
c        -----------------------------------------------------
c        calculation of the new value rrp(i) using the equation:
c        (zzp(i)-zzp(i-1))/(rrp(i)-rrp(i-1))=
c        (zzp(j0)-zzp(i-1))/(rrp(i0)-rrp(i-1))
c        -----------------------------------------------------
            if(zzp(j0).eq.zzp(i-1)) then
	       rrp(i)=rrp(i)+1.d-5
	    else
	       rrp(i)=rrp(i-1)+(zzp(i)-zzp(i-1))*(rrp(j0)-rrp(i-1))/
     1	       	      (zzp(j0)-zzp(i-1))
	    endif
	 endif
      enddo
c--------near rmax
c     ------------------------------------
c     the correction for the ip-1 point (begin)

c      write(*,*)'before correction ip',ip
c      write(*,*)'before correction rrp(ip-1),rrp(ip)',rrp(ip-1),rrp(ip)
c      write(*,*)'before correction zzp(ip-1),zzp(ip)',zzp(ip-1),zzp(ip)

      if (rrp(ip-1).gt.rrp(ip)-1.d-5) then

cSAP091103
         j0=ip-1

         do j=ip-2,izmax,-1
            if (rrp(j).lt.rrp(ip-1)) then
	       j0=j
	       goto 55
	    endif
	 enddo
 55      continue
c         write(*,*)'j0,zzp(j0),zzp(ip)',j0,zzp(j0),zzp(ip)
         if(zzp(j0).eq.zzp(ip)) then
	    rrp(ip-1)=rrp(ip)-1.d-5
c	    write(*,*)'1 cor rrp(ip-1)',rrp(ip-1)
	 else
cSmirnov961226 beg
	    if(dabs(rrp(ip-1)-rrp(ip)).lt.1.d-6) then
	       rrp(ip-1)=0.5d0*(rrp(j0)+rrp(ip))
	    else
	       rrp(ip-1)=rrp(ip)+(zzp(ip-1)-zzp(ip))*(rrp(j0)-rrp(ip))/
     1	       	      (zzp(j0)-zzp(ip))
	    endif
cSmirnov961226 end
c	    write(*,*)'2 cor rrp(ip-1)',rrp(ip-1)
	 endif
      endif
c      write(*,*)'after correction rrp(ip-1),rrp(ip)',rrp(ip-1),rrp(ip)
c      write(*,*)'after correction zzp(ip-1),zzp(ip)',zzp(ip-1),zzp(ip)
c     the correction for the ip-1 point (end)
c     ------------------------------------

      do i=ip-2,izmax+1,-1
c         if(rrp(i).ge.rrp(i+1)) then
         if(rrp(i).ge.(rrp(i+1)-5.d-3)) then
	   j0=izmax
	   do j=i-1,izmax,-1
	     if (rrp(j).lt.rrp(i+1)) then
	       j0=j
	       goto 60
	     endif
	   enddo
 60        continue
c        -----------------------------------------------------
c        calculation of the new value rrp(i) using the equation:
c        (zzp(i)-zzp(i+1))/(rrp(i)-rrp(i+1))=
c        (zzp(j0)-zzp(i+1))/(rrp(j0)-rrp(i+1))
c        -----------------------------------------------------
            if(zzp(j0).eq.zzp(i+1)) then
	       rrp(i)=rrp(i)-1.d-5
	    else
	       rrp(i)=rrp(i+1)+(zzp(i)-zzp(i+1))*(rrp(j0)-rrp(i+1))/
     1	       	      (zzp(j0)-zzp(i+1))
	    endif
	 endif
      enddo

c-----------------------------------------------------------------
c  creation of monotonic limiter array zzm  near rmin and rmax

      izmin=1
      zmin=zzm(1)
      do i=2,im
        if (zmin.gt.zzm(i)) then
          zmin=zzm(i)
          izmin=i
	endif
      enddo
c      write(*,*)'im,izmin,zzm(izmin),rrm(izmin)',
c     1 im,izmin,zzm(izmin),rrm(izmin)

c--------near rmin
c     the correction for the second point (begin)
c      write(*,*)'before correction rrm(2),rrm(1)',rrm(2),rrm(1)
c      write(*,*)'before correction zzm(2),zzm(1)',zzm(2),zzm(1)
      if (rrm(2).lt.rrm(1)+1.d-5) then
         do j=3,izmax

cSAP091102
            j0=2

            if (rrm(j).gt.rrm(2)) then
	       j0=j
	       goto 65
	    endif
	 enddo
 65      continue
         if(zzm(j0).eq.zzm(1)) then
	    rrm(2)=rrm(1)+1.d-5
	 else
cSmirnov961226 beg
	    if(dabs(rrm(1)-rrm(2)).lt.1.d-6) then
	       rrm(2)=0.5d0*(rrm(j0)+rrp(1))
	    else
	       rrm(2)=rrm(1)+(zzm(2)-zzm(1))*(rrm(j0)-rrm(1))/
     1	       	      (zzm(j0)-zzm(1))
	    endif
cSmirnov961226 end
	 endif
      endif
c      write(*,*)'after correction rrm(2),rrm(1)',rrm(2),rrm(1)
c      write(*,*)'after correction zzm(2),zzm(1)',zzm(2),zzm(1)
c     the correction for the second point (end)

      do i=3,izmin-1
c         if(rrm(i).le.rrm(i-1)) then
         if(rrm(i).le.(rrm(i-1)+5.d-3)) then
	   j0=izmin
	   do j=i+1,izmin
	     if (rrm(j).gt.rrm(i-1)) then
	       j0=j
	       goto 70
	     endif
	   enddo
 70        continue
c        -----------------------------------------------------
c        calculation of the new value rrp(i) using the equation:
c        (zzm(i)-zzm(i-1))/(rrm(i)-rrm(i-1))=
c        (zzm(j0)-zzm(i-1))/(rrm(i0)-rrm(i-1))
c        -----------------------------------------------------
            if(zzm(j0).eq.zzm(i-1)) then
	       rrm(i)=rrm(i)+1.d-5
	    else
	       rrm(i)=rrm(i-1)+(zzm(i)-zzm(i-1))*(rrm(j0)-rrm(i-1))/
     1	       	      (zzm(j0)-zzm(i-1))
	    endif
	 endif
      enddo
c--------near rmax
c     ------------------------------------
c     the correction for the im-1 point (begin)
c      write(*,*)'before correction im',im
c      write(*,*)'before correction rrm(im-1),rrm(im)',rrm(im-1),rrm(im)
c      write(*,*)'before correction zzm(im-1),zzm(im)',zzm(im-1),zzm(im)
cSAP0-91103
c      if (rrm(ip-1).gt.rrm(im)-1.d-5) then 
      if (rrm(im-1).gt.rrm(im)-1.d-5) then

cSAP091103
         j0=im-1

         do j=im-2,izmax,-1
            if (rrm(j).lt.rrm(im-1)) then
	       j0=j
	       goto 75
	    endif
	 enddo
 75      continue
         if(zzm(j0).eq.zzm(im)) then
	    rrm(im-1)=rrm(im)-1.d-5
	 else
cSmirnov961226 beg
	    if(dabs(rrm(im-1)-rrm(im)).lt.1.d-6) then
	       rrm(im-1)=0.5d0*(rrm(j0)+rrm(im))
	    else
	       rrm(im-1)=rrm(im)+(zzm(im-1)-zzm(im))*(rrm(j0)-rrm(im))/
     1	       	      (zzm(j0)-zzm(im))
	    endif
cSmirnov961226 end
	 endif
      endif
c      write(*,*)'after correction rrm(im-1),rrm(im)',rrm(im-1),rrm(im)
c      write(*,*)'after correction zzm(im-1),zzm(im)',zzm(im-1),zzm(im)
c     the correction for the im-1 point (end)
c     ------------------------------------
      do i=im-2,izmin+1,-1
c         if(rrm(i).ge.rrm(i+1)) then
         if(rrm(i).ge.(rrm(i+1)-5.d-3)) then
	   j0=izmin
	   do j=i-1,izmin,-1
	     if (rrm(j).lt.rrm(i+1)) then
	       j0=j
	       goto 80
	     endif
	   enddo
 80        continue
c        -----------------------------------------------------
c        calculation of the new value rrm(i) using the equation:
c        (zzm(i)-zzm(i+1))/(rrm(i)-rrm(i+1))=
c        (zzm(j0)-zzm(i+1))/(rrm(j0)-rrm(i+1))
c        -----------------------------------------------------
            if(zzm(j0).eq.zzm(i+1)) then
	       rrm(i)=rrm(i)-1.d-5
	    else
	       rrm(i)=rrm(i+1)+(zzm(i)-zzm(i+1))*(rrm(j0)-rrm(i+1))/
     1	       	      (zzm(j0)-zzm(i+1))
	    endif
	 endif
      enddo

c-----------------------------------------------------------------
c  end of arrays rrp,zzp  and  rrm,zzm  creation
c------------------------------------------------------------------
c      write(*,*)'in limitr after monotonization ip,im',ip,im
      do i=1,ip
        rpl(i)=rrp(i)
        zpl(i)=zzp(i)
c        write(*,*)'i,rrp(i),zzp(i)',i,rrp(i),zzp(i)
      enddo
      do i=1,im
        rml(i)=rrm(i)
        zml(i)=zzm(i)
c        write(*,*)'in i,rrm(i),zzm(i)',i,rrm(i),zzm(i)
      enddo
      return
      end


c        ************************FPSI**************************
c        *                        -                           *
c        * this function calculates psi                       *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r  the coordinates
c------------------------------------------------------------------
      double precision function fpsi(r,z)
      !implicit double precision (a-h,o-z)
      implicit none !YuP[2020-01-14]
      integer nx4,ny4 !YuP[2020-01-14]
      real*8 z,r !YuP[2020-01-14]

      include 'param.i'
      include 'fourb.i'
      include 'five.i'
       double precision ias2r,ias2r_Sm
              nx4=nx+4
              ny4=ny+4
	      ncx=nx4
	      ncy=ny4
cSm030224
c              fpsi=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
              fpsi=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nx4a)

      return
      end


c
c***********************MAKELAC0****************************************
c      it calculates the parameters
c      zmax1,rmaz1,zmin1,rmin1 for the eqdsk without
c      limiter data (nnlim=0)
c      zmax1,rmaz1,zmin1,rmin1 are closed but inside the
c      boundaries of the Lackner rectangle
c             ----------------------
c            |		. zmax1	    |
c            |			    |
c            |	   the points       |
c            |	     of the	    |
c      rmin1 | .     last	  . |rmax1
c            |	     closed	    |
c            |	     flux	    |
c            |	    surface	    | Lackner Rectangle
c            |			    |
c            |	       . zmin1	    |
c             ----------------------
c
c     (these points are inside the  separatrix line)
c----------------------------------------------------------------------
c      input data:
c        peqd(nxeqda,nyeqda) is the psi function ,in common/fourb/
c        psisep,psilim     in common/three/
c        xma,yma           coordinates of magnetic axis in common/three/
c        nxeqd,nyeqd       the number of points in the horizontal(r)
c                          and vertical(z) directions,  in common/three/
c        rr(nxeqda),zz(nyeqda) eqdsk mesh points, in common/fourb/
c-----------------------------------------------------------------------
c      output data:
c         zmax1,rmax1,zmin1,rmin1 -boundaries of Lackner rectangle
c         rzmax1,zrmax1,rzmin1,zrmin1)
c-----------------------------------------------------------------------
      subroutine makelac0(zmax1,rmax1,zmin1,rmin1,
     1 rzmax1,zrmax1,rzmin1,zrmin1)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'three.i'
      include 'fourb.i'
      dimension ix(nxeqda)
      dimension rminj(nyeqda),rmaxj(nyeqda)
      dimension irminj(nyeqda),irmaxj(nyeqda)
      dimension psiminj(nyeqda),rpsiminj(nyeqda)
      write(*,*)'in equilib in makelac0'
c------------------------------------------------------------
c     determination of eqdsk mesh point close to the magnetic axis
c     zz(imag) close to yma, rr(jmag) close to xma
c------------------------------------------------------------
      dif=10.d15
      imag=1
      do i=1,nxeqd
	difnew=dabs(rr(i)-xma)
	if (difnew.lt.dif) then
	   dif=difnew
	   imag=i
	endif
      enddo !i

      dif=10.d15
      jmag=1
      do j=1,nyeqd
	difnew=dabs(zz(j)-yma)
	if (difnew.lt.dif) then
	   dif=difnew
	   jmag=j
	endif
      enddo !i

c      psiminj(jmag)=peqd(imag,jmag)
c      rpsiminj(jmag)=rr(imag)
      write(*,*)'imag,rr(imag),xma',imag,rr(imag),xma
      write(*,*)'jmag,zz(jmag),yma',jmag,zz(jmag),yma
      write(*,*)'peqd(imag,jmag),psimag',peqd(imag,jmag),psimag
c------------------------------------------------------------
c     Determination of rminj(nyeqd),rmaxj(nyeqd) -min and max major
c     radius(at given j) of the eqdsk points which are outside the plasma
c     and have the minimum distance from the plasma.
c     Determination of irminj(nyeqd),irmaxj(nyeqd) -the numbers in array r(i)
c     rminj(j)=r(irminj(j)), rmaxj(j)=r(irmax(j))
c     Determination of psiminj(nyeqd) and rpsimin(nyeqd)
c     psiminj(j) is the minimum value of psi along the major radius
c     from r=rminj(j) to r =rmaxj(j).
c     psiminj(j)=min{i=irminj(j)+1,irmaxj(j)-1} peqd(i,j)
c     rpsimin(j) is the value of major radius in which peqd(i,j)=psiminj(j)
c
c           {z(j),r(irminj(j)}             	{z(j),r(irmaxj(j)}
c     z(j) : .  .  .  *|  .  .   .   .   .  .  .  | *  .  .
c	        rminj(j)         plasma	            rmaxj(j)
c	       irminj(j)         	            irmaxj(j)
c------------------------------------------------------------
      jzmin=nyeqd
      jzmax=1
      do j=jmag,nyeqd
c     -------------------------------------------------------
c     for the eqdsk points with z.ge.yma
c     -------------------------------------------------------
	 ki=0
c        -------------------------------------------------------
c        ix(ki) is the array of the numbers i of eqdsk mesh rr(i) for which
c        peqd(ix,j) .lt. psilim	, ki is the quantity of the such points
c        -------------------------------------------------------
         do i=1,nxeqd
	    if (peqd(i,j).le.psilim) then
	       ki=ki+1
	       ix(ki)=i
	    endif
	 enddo   !i

	 if (ki.ge.1) then
c           ---------------------------------------------------
c           in this case there are exist the points in which
c           peqd(i,j).le.psilim for the given j
c           determination of the point of
c           min{i=irmin,irmax}peqd(i,j)  	for the given j
c           ---------------------------------------------------
	    irmin=ix(1)
	    irmax=ix(ki)
	    jzmax=j
	    psiminj(j)=peqd(irmin,j)
	    rpsiminj(j)=rr(irmin)
	    do i=1,ki
	       ip=ix(i)
	       if(peqd(ip,j).lt.psiminj(j)) then
	          psiminj(j)=peqd(ip,j)
	          rpsiminj(j)=rr(ip)
	       endif
	    enddo !i

	    rminj(j)=rr(irmin)
	    irminj(j)=irmin
	    rmaxj(j)=rr(irmax)
	    irmaxj(j)=irmax
	 else
c           ---------------------------------------------------
c           in this case there are no points in which
c           peqd(i,j).le.psilim for the given j
c           ---------------------------------------------------
	    rzmax1=rpsiminj(jzmax)
            goto 10
         endif

cc         write(*,*)'j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)',
cc     1              j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)
      enddo !j
 10   continue
      do j=jmag-1,1,-1
c     -------------------------------------------------------
c     for the eqdsk points with z.lt.yma
c     -------------------------------------------------------
         ki=0
c        -------------------------------------------------------
c        ix(ki) is the array of the numbers i of eqdsk mesh rr(i) for which
c        peqd(ix,j) .lt. psilim	, ki is the quantity of the such points
c        -------------------------------------------------------
         do i=1,nxeqd
	    if (peqd(i,j).le.psilim) then
	       ki=ki+1
	       ix(ki)=i
	    endif
	 enddo   !i

	 if (ki.ge.1) then
c           ---------------------------------------------------
c           in this case there are exist points where
c           peqd(i,j).le.psilim for the given j	,
c           determination of the point of
c           min{i=irmin,irmax}peqd(i,j)  	for the given j
c           ---------------------------------------------------
	    irmin=ix(1)
	    irmax=ix(ki)
	    jzmin=j
	    psiminj(j)=peqd(irmin,j)
	    rpsiminj(j)=rr(irmin)
	    do i=1,ki
	       ip=ix(i)
	       if(peqd(ip,j).lt.psiminj(j)) then
	          psiminj(j)=peqd(ip,j)
	          rpsiminj(j)=rr(ip)
	       endif
	    enddo !i

	    rminj(j)=rr(irmin)
	    irminj(j)=irmin
	    rmaxj(j)=rr(irmax)
	    irmaxj(j)=irmax
	 else
c           ---------------------------------------------------
c           in this case there are no points in which
c           peqd(i,j).le.psilim for the given j
c           ---------------------------------------------------
cc	      write(*,*)'out of plasmaj,zz(j)',j,zz(j)
	    rzmin1=rpsiminj(jzmin)
	    goto 20
         endif

cc         write(*,*)'j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)',
cc     1              j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)
      enddo !j
 20   continue
c-----------------------------------------------------------------------
c     determination of zmax1 and zmin1
c-----------------------------------------------------------------------
cc      write(*,*)'jzmin,jzmax,nyeqd',jzmin,jzmax,nyeqd
      if(jzmax.lt.nyeqd) then
         zmax1=zz(jzmax)
      else
c        ---------------------------------------------------
c        it is necessary to calculate the jzmax1
c        (in this case the separatrix area exist close to plasma)
c        ---------------------------------------------------
         do j=jmag+1,nyeqd
cc	    write(*,*)'j,psiminj(j),psiminj(j-1)'
cc	    write(*,*)j,psiminj(j),psiminj(j-1)
            if(psiminj(j).gt.psiminj(j-1))then
	      zmax1=zz(j)
	      jzmax=j
cc	      write(*,*)'j,zmax1,jzmax',j,zmax1,jzmax
	    else
	       goto30
            endif
         enddo
 30      continue
cc	 jzmax=jzmax-1
	 jzmax=jzmax
	 zmax1=zz(jzmax)
	 rzmax1=rpsiminj(jzmax)
cc         write(*,*)'jzmax,zmax1',jzmax,zmax1
      endif

      if(jzmin.gt.1) then
         zmin1=zz(jzmin)
      else
c        ---------------------------------------------------
c        it is necessary to calculate the jzmin1
c        (in this case the separatrix area exist close to plasma)
c        ---------------------------------------------------
         do j=jmag-1,1,-1
cc	    write(*,*)'j,psiminj(j),psiminj(j+1)'
cc	    write(*,*)j,psiminj(j),psiminj(j+1)
            if(psiminj(j).gt.psiminj(j+1))then
	       zmin1=zz(j)
	       jzmin=j
cc	       write(*,*)'j,zmin1,jzmin',j,zmin1,jzmin
	    else
	       goto40
            endif
         enddo
 40      continue
cc	 jzmin=jzmin+1
	 jzmin=jzmin
	 zmin1=zz(jzmin)
	 rzmin1=rpsiminj(jzmin)
      endif
cc      write(*,*)'jzmin,zmin1',jzmin,zmin1
c---------------------------------------------------------------
      rmax1=rr(1)
      rmin1=rr(nxeqd)
c-------------------------------------------
c      determination of rmax1 and rmin1
c-------------------------------------------
cc      write(*,*)'jzmin,jzmax,rmin1,rmax1',jzmin,jzmax,rmin1,rmax1
      do j=jzmin,jzmax
cc	 write(*,*)'j, rminj(j),rmin1',j,rminj(j),rmin1
	 if (rminj(j).lt.rmin1) then
cc	 write(*,*)'rminj(j),rmin1',rminj(j),rmin1
	    rmin1=rminj(j)
	    zrmin1=zz(j)
	 endif
cc	 write(*,*)'rmin1',rmin1
cc	 write(*,*)'j, rmaxj(j),rmax1',j,rmaxj(j),rmax1
	 if (rmaxj(j).gt.rmax1) then
	    rmax1=rmaxj(j)
	    zrmax1=zz(j)
	 endif
cc	 write(*,*)'rmax1',rmax1
      enddo !j

cc      write(*,*)'rmax1,rmin1,zmax1,zmin1',rmax1,rmin1,zmax1,zmin1
      return
      end


c***********************MAKELAC1****************************************
c      it calculates the parameters
c      zmax1,rmaz1,zmin1,rmin1 for the eqdsk without
c      limiter data (nnlim=0)
c      zmax1,rmaz1,zmin1,rmin1 are closed but outside the
c      boundaries of the Lackner rectangle
c            		. zmax1
c             ----------------------
c            |			    |
c            |	   the points       |
c            |	     of the	    |
c    rmin1 . |      last	    |. rmax1
c            |	     closed	    |
c            |	     flux	    |
c            |	    surface	    | Lackner Rectangle
c            |			    |
c            |	         	    |
c             ----------------------
c		       . zmin1
c     (these points are outside the  separatrix line)
c----------------------------------------------------------------------
c      input data:
c        peqd(nxeqda,nyeqda) is the psi function ,in common/fourb/
c        psisep,psilim     in common/three/
c        xma,yma           coordinates of magnetic axis in common/three/
c        nxeqd,nyeqd       the number of points in the horizontal(r)
c                          and vertical(z) directions,  in common/three/
c        rr(nxeqda),zz(nyeqda) eqdsk mesh points, in common/fourb/
c-----------------------------------------------------------------------
c      output data:
c         zmax1,rmax1,zmin1,rmin1 -boundaries of Lackner rectangle
c-----------------------------------------------------------------------
      subroutine makelac1(zmax1,rmax1,zmin1,rmin1,
     1 rzmax1,zrmax1,rzmin1,zrmin1)
c it is necessary to add the calculation of the last 4 data
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'three.i'
      include 'fourb.i'
      dimension ix(nxeqda)
      dimension rminj(nyeqda),rmaxj(nxeqda)
      dimension psiminj(nyeqda),rpsiminj(nyeqda)
c------------------------------------------------------------
c     determination of eqdsk mesh point close to magnetic axis
c     zz(imag) close to yma, rr(jmag) close to xma
c------------------------------------------------------------
      dif=10.d15
      imag=1
      do i=1,nxeqd
	difnew=dabs(rr(i)-xma)
	if (difnew.lt.dif) then
	   dif=difnew
	   imag=i
	endif
      enddo !i

      dif=10.d15
      jmag=1
      do j=1,nyeqd
	difnew=dabs(zz(j)-yma)
	if (difnew.lt.dif) then
	   dif=difnew
	   jmag=j
	endif
      enddo !i

      psiminj(jmag)=peqd(imag,jmag)
      rpsiminj(jmag)=rr(imag)
      write(*,*)'imag,rr(imag),xma',imag,rr(imag),xma
      write(*,*)'jmag,zz(jmag),yma',jmag,zz(jmag),yma
      write(*,*)'peqd(imag,jmag),psimag',peqd(imag,jmag),psimag
c------------------------------------------------------------
c     determination of rminj(nyeqd),rmaxj(nyeqd),
c                      psiminj(nyeqd),rpsimin(nyeqd)
c------------------------------------------------------------
      jzmin=nyeqd
      jzmax=1
      do j=jmag,nyeqd
c     -------------------------------------------------------
c     for the eqdsk points with z.ge.yma
c     -------------------------------------------------------
	 ki=0
c        -------------------------------------------------------
c        ix(ki) is the array of the numbers i of eqdsk mesh rr(i) for which
c        peqd(ix,j) .lt. psilim	, ki is the quantity of the such points
c        -------------------------------------------------------
         do i=1,nxeqd
	     write(*,*)'1makelacn ,i,nxeqd',i,nxeqd
	    if (peqd(i,j).le.psilim) then
	       ki=ki+1
	       ix(ki)=i
	     write(*,*)'2makelacn ,ki,i',ki,i
	    endif
	 enddo   !i

	 if (ki.ge.1) then
c           ---------------------------------------------------
c           in this case there are exist the points in which
c           peqd(i,j).le.psilim for the given j
c           determination of the point of
c           min{i=irmin,irmax}peqd(i,j)  	for the given j
c           ---------------------------------------------------
	    irmin=ix(1)-1
	    irmax=ix(ki)+1
	    jzmax=j
	    psiminj(j)=peqd(irmin,j)
	    rpsiminj(j)=rr(irmin)
	    do i=1,ki
	       ip=ix(i)
	     write(*,*)'3makelacn ,i,ki,ip',i,ki,ip
	       if(peqd(ip,j).lt.psiminj(j)) then
	     write(*,*)'4makelacn ,j,ip,rr(ip)',j,i,rr(ip)
	          psiminj(j)=peqd(ip,j)
	          rpsiminj(j)=rr(ip)
	       endif
	    enddo !i

	     write(*,*)'5rmakelacn,j,irmax,irmin,rr(irmax),rr(irmin)',
     1	     j,irmax,irmin,rr(irmax),rr(irmin)
	    rminj(j)=rr(irmin)
	    rmaxj(j)=rr(irmax)
	     write(*,*)'6makelacn,j,rminj(j),rmaxj(j)',
     1	     j,rminj(j),rmaxj(j)
	 else
c           ---------------------------------------------------
c           in this case there are no points in which
c           peqd(i,j).le.psilim for the given j
c           ---------------------------------------------------
	    jzmax=j
            goto 10
         endif

cc         write(*,*)'j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)',
cc     1              j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)
      enddo !j
 10   continue
      do j=jmag-1,1,-1
c     -------------------------------------------------------
c     for the eqdsk points with z.lt.yma
c     -------------------------------------------------------
         ki=0
c        -------------------------------------------------------
c        ix(ki) is the array of the numbers i of eqdsk mesh rr(i) for which
c        peqd(ix,j) .lt. psilim	, ki is the quantity of the such points
c        -------------------------------------------------------
         do i=1,nxeqd
	    if (peqd(i,j).le.psilim) then
	       ki=ki+1
	       ix(ki)=i
	    endif
	 enddo   !i

	 if (ki.ge.1) then
c           ---------------------------------------------------
c           in this case there are exist points where
c           peqd(i,j).le.psilim for the given j	,
c           determination of the point of
c           min{i=irmin,irmax}peqd(i,j)  	for the given j
c           ---------------------------------------------------
	    irmin=ix(1)-1
	    irmax=ix(ki)+1
	    jzmin=j
	    psiminj(j)=peqd(irmin,j)
	    rpsiminj(j)=rr(irmin)
	    do i=1,ki
	       ip=ix(i)
	       if(peqd(ip,j).lt.psiminj(j)) then
	          psiminj(j)=peqd(ip,j)
	          rpsiminj(j)=rr(ip)
	       endif
	    enddo !i

	    rminj(j)=rr(irmin)
	    rmaxj(j)=rr(irmax)
	 else
c           ---------------------------------------------------
c           in this case there are no points in which
c           peqd(i,j).le.psilim for the given j
c           ---------------------------------------------------
cc	      write(*,*)'out of plasmaj,zz(j)',j,zz(j)
	    jzmin=j
	    goto 20
         endif

cc         write(*,*)'j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)',
cc     1              j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)
      enddo !j
 20   continue
c-----------------------------------------------------------------------
c     determination of zmax1 and zmin1
c-----------------------------------------------------------------------
cc      write(*,*)'jzmin,jzmax,nyeqd',jzmin,jzmax,nyeqd
      if(jzmax.lt.nyeqd) then
         zmax1=zz(jzmax)
      else
c        ---------------------------------------------------
c        it is necessary to calculate the jzmax1
c        (in this case the separatrix area exist close to plasma)
c        ---------------------------------------------------
         do j=jmag+1,nyeqd
cc	    write(*,*)'j,psiminj(j),psiminj(j-1)'
cc	    write(*,*)j,psiminj(j),psiminj(j-1)
            if(psiminj(j).gt.psiminj(j-1))then
	      jzmax=j
	      zmax1=zz(jzmax)
cc	      write(*,*)'j,zmax1,jzmax',j,zmax1,jzmax
	    else
	       goto30
            endif
         enddo
 30      continue
	 jzmax=jzmax+1
	 zmax1=zz(jzmax)
cc         write(*,*)'jzmax,zmax1',jzmax,zmax1
      endif

      if(jzmin.gt.1) then
         zmin1=zz(jzmin)
      else
c        ---------------------------------------------------
c        it is necessary to calculate the jzmin1
c        (in this case the separatrix area exist close to plasma)
c        ---------------------------------------------------
         do j=jmag-1,1,-1
cc	    write(*,*)'j,psiminj(j),psiminj(j+1)'
cc	    write(*,*)j,psiminj(j),psiminj(j+1)
            if(psiminj(j).gt.psiminj(j+1))then
	       zmin1=zz(j)
	       jzmin=j
cc	       write(*,*)'j,zmin1,jzmin',j,zmin1,jzmin
	    else
	       goto40
            endif
         enddo
 40      continue
	 jzmin=jzmin-1
	 zmin1=zz(jzmin)
      endif
cc      write(*,*)'jzmin,zmin1',jzmin,zmin1
c---------------------------------------------------------------
      rmax1=rr(1)
      rmin1=rr(nxeqd)
c-------------------------------------------
c      determination of rmax1 and rmin1
c-------------------------------------------
cc      write(*,*)'jzmin,jzmax,rmin1,rmax1',jzmin,jzmax,rmin1,rmax1
      do j=jzmin+1,jzmax-1
cc	 write(*,*)'j, rminj(j),rmin1',j,rminj(j),rmin1
	 if (rminj(j).lt.rmin1) then
cc	 write(*,*)'rminj(j),rmin1',rminj(j),rmin1
	    rmin1=rminj(j)
	 endif
cc	 write(*,*)'rmin1',rmin1
cc	 write(*,*)'j, rmaxj(j),rmax1',j,rmaxj(j),rmax1
	 if (rmaxj(j).gt.rmax1) then
	    rmax1=rmaxj(j)
	 endif
cc	 write(*,*)'rmax1',rmax1
      enddo !j

cc      write(*,*)'rmax1,rmin1,zmax1,zmin1',rmax1,rmin1,zmax1,zmin1
      return
      end


c***********************MAKELACN****************************************
c      it calculates the parameters
c      zmax1,rmaz1,zmin1,rmin1 for the eqdsk without
c      limiter data (nnlim=0)
c      zmax1,rmaz1,zmin1,rmin1 are closed but inside the
c      boundaries of the Lackner rectangle
c             ----------------------
c            |		. zmax1	    |
c            |			    |
c            |	   the points       |
c            |	     of the	    |
c     rmin1 .|      last	    |.rmax1
c            |	     closed	    |
c            |	     flux	    |
c            |	    surface	    | Lackner Rectangle
c            |			    |
c             ----------------------
c            	       . zmin1
c
c     (these points are inside the  separatrix line)
c----------------------------------------------------------------------
c      input data:
c        peqd(nxeqda,nyeqda) is the psi function ,in common/fourb/
c        psisep,psilim     in common/three/
c        xma,yma           coordinates of magnetic axis in common/three/
c        nxeqd,nyeqd       the number of points in the horizontal(r)
c                          and vertical(z) directions,  in common/three/
c        rr(nxeqda),zz(nyeqda) eqdsk mesh points, in common/fourb/
c-----------------------------------------------------------------------
c      output data:
c         zmax1,rmax1,zmin1,rmin1 -boundaries of Lackner rectangle
c         rzmax1,zrmax1,rzmin1,zrmin1)
c-----------------------------------------------------------------------
      subroutine makelacn(zmax1,rmax1,zmin1,rmin1,
     1 rzmax1,zrmax1,rzmin1,zrmin1)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'three.i'
      include 'fourb.i'
      dimension ix(nxeqda)
      dimension rminj(nyeqda),rmaxj(nyeqda)
      dimension irminj(nyeqda),irmaxj(nyeqda)
      dimension psiminj(nyeqda),rpsiminj(nyeqda)
c      write(*,*)'in equilib in makelacn'
c------------------------------------------------------------
c     determination of eqdsk mesh point close to magnetic axis
c     zz(imag) close to yma, rr(jmag) close to xma
c------------------------------------------------------------
c      write(*,*)'in makelacn nxeqd,nyeqd',nxeqd,nyeqd
      dif=10.d15
      imag=1
      do i=1,nxeqd
	difnew=dabs(rr(i)-xma)
	if (difnew.lt.dif) then
	   dif=difnew
	   imag=i
	endif
      enddo !i

      dif=10.d15
      jmag=1
      do j=1,nyeqd
	difnew=dabs(zz(j)-yma)
	if (difnew.lt.dif) then
	   dif=difnew
	   jmag=j
	endif
      enddo !i

      psiminj(jmag)=peqd(imag,jmag)
      rpsiminj(jmag)=rr(imag)
c      write(*,*)'imag,rr(imag),xma',imag,rr(imag),xma
c      write(*,*)'jmag,zz(jmag),yma',jmag,zz(jmag),yma
c      write(*,*)'peqd(imag,jmag),psimag',peqd(imag,jmag),psimag
c------------------------------------------------------------
c     Determination of rminj(nyeqd),rmaxj(nyeqd) -min and max major
c     radius(at given j) of the eqdsk points which are outside the plasma
c     and have the minimum distance from the plasma.
c     Determination of irminj(nyeqd),irmaxj(nyeqd) -the numbers in array r(i)
c     rminj(j)=r(irminj(j)), rmaxj(j)=r(irmax(j))
c     Determination of psiminj(nyeqd) and rpsimin(nyeqd)
c     psiminj(j) is the minimum value of psi along the major radius
c     from r=rminj(j) to r =rmaxj(j).
c     psiminj(j)=min{i=irminj(j)+1,irmaxj(j)-1} peqd(i,j)
c     rpsimin(j) is the value of major radius in which peqd(i,j)=psiminj(j)
c
c           {z(j),r(irminj(j)}             	{z(j),r(irmaxj(j)}
c     z(j) : .  .  .  *|  .  .   .   .   .  .  .  | *  .  .
c	        rminj(j)         plasma	            rmaxj(j)
c	       irminj(j)         	            irmaxj(j)
c------------------------------------------------------------
      jzmin=nyeqd
      jzmax=1
      do j=jmag,nyeqd
c     -------------------------------------------------------
c     for the eqdsk points with z.ge.yma
c     -------------------------------------------------------
	 ki=0
c        -------------------------------------------------------
c        ix(ki) is the array of the numbers i of eqdsk mesh rr(i) for which
c        peqd(ix,j) .lt. psilim	, ki is the quantity of the such points
c        -------------------------------------------------------
         do i=1,nxeqd
	    if (peqd(i,j).le.psilim) then
	       ki=ki+1
	       ix(ki)=i
	    endif
	 enddo   !i

	 if (ki.ge.1) then
c           ---------------------------------------------------
c           in this case there are exist the points in which
c           peqd(i,j).le.psilim for the given j
c           determination of the point of
c           min{i=irmin,irmax}peqd(i,j)  	for the given j
c           ---------------------------------------------------
	    irmin=ix(1)-1
	    irmax=ix(ki)+1
	    if (irmin.lt.1) irmin=1
	    if (irmax.gt.nxeqd) irmax=nxeqd
	    jzmax=j
	    psiminj(j)=peqd(irmin,j)
	    rpsiminj(j)=rr(irmin)
	    do i=1,ki
	       ip=ix(i)
	       if(peqd(ip,j).lt.psiminj(j)) then
	          psiminj(j)=peqd(ip,j)
	          rpsiminj(j)=rr(ip)
	       endif
	    enddo !i
c	    write(*,*)'in equilib irmin',irmin
	    rminj(j)=rr(irmin)
	    irminj(j)=irmin
c	    write(*,*)'in equilib irmax',irmax
	    rmaxj(j)=rr(irmax)
	    irmaxj(j)=irmax
c         write(*,*)'j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)',
c     1              j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)
	 else
c           ---------------------------------------------------
c           in this case there are no points in which
c           peqd(i,j).le.psilim for the given j
c           ---------------------------------------------------
	    rminj(j)=rr(nxeqd)
	    irminj(j)=nxeqd
	    rmaxj(j)=rr(1)
	    irmaxj(j)=1

	    jzmax=j
	    rzmax1=rpsiminj(jzmax)
            goto 10
         endif

      enddo !j
 10   continue
      do j=jmag-1,1,-1
c     -------------------------------------------------------
c     for the eqdsk points with z.lt.yma
c     -------------------------------------------------------
         ki=0
c        -------------------------------------------------------
c        ix(ki) is the array of the numbers i of eqdsk mesh rr(i) for which
c        peqd(ix,j) .lt. psilim	, ki is the quantity of the such points
c        -------------------------------------------------------
         do i=1,nxeqd
	    if (peqd(i,j).le.psilim) then
	       ki=ki+1
	       ix(ki)=i
	    endif
	 enddo   !i

	 if (ki.ge.1) then
c           ---------------------------------------------------
c           in this case there are exist points where
c           peqd(i,j).le.psilim for the given j	,
c           determination of the point of
c           min{i=irmin,irmax}peqd(i,j)  	for the given j
c           ---------------------------------------------------
	    irmin=ix(1)-1
	    irmax=ix(ki)+1
	    if (irmin.lt.1) irmin=1
	    if (irmax.gt.nxeqd) irmax=nxeqd
	    jzmin=j
	    psiminj(j)=peqd(irmin,j)
	    rpsiminj(j)=rr(irmin)
	    do i=1,ki
	       ip=ix(i)
	       if(peqd(ip,j).lt.psiminj(j)) then
	          psiminj(j)=peqd(ip,j)
	          rpsiminj(j)=rr(ip)
	       endif
	    enddo !i

	    rminj(j)=rr(irmin)
	    irminj(j)=irmin
	    rmaxj(j)=rr(irmax)
	    irmaxj(j)=irmax
c         write(*,*)'j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)',
c     1              j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)
	 else
c           ---------------------------------------------------
c           in this case there are no points in which
c           peqd(i,j).le.psilim for the given j
c           ---------------------------------------------------
cc	      write(*,*)'out of plasmaj,zz(j)',j,zz(j)
	    rminj(j)=rr(nxeqd)
	    irminj(j)=nxeqd
	    rmaxj(j)=rr(1)
	    irmaxj(j)=1

	    jzmin=j
	    rzmin1=rpsiminj(jzmin)
	    goto 20
         endif

      enddo !j
 20   continue
c-----------------------------------------------------------------------
c     determination of zmax1 and zmin1
c-----------------------------------------------------------------------
cc      write(*,*)'jzmin,jzmax,nyeqd',jzmin,jzmax,nyeqd
      if(jzmax.lt.nyeqd) then
         zmax1=zz(jzmax)
	 rzmax1=rpsiminj(jzmax-1)
      else
c        ---------------------------------------------------
c        it is necessary to calculate the jzmax1
c        (in this case the separatrix area exist close to plasma)
c        ---------------------------------------------------
         do j=jmag+1,nyeqd
cc	    write(*,*)'j,psiminj(j),psiminj(j-1)'
cc	    write(*,*)j,psiminj(j),psiminj(j-1)
            if(psiminj(j).gt.psiminj(j-1))then
	      zmax1=zz(j)
	      jzmax=j
cc	      write(*,*)'j,zmax1,jzmax',j,zmax1,jzmax
	    else
	       goto30
            endif
         enddo
 30      continue
cc	 jzmax=jzmax-1
	 jzmax=jzmax
	 zmax1=zz(jzmax)
	 rzmax1=rpsiminj(jzmax)
cc         write(*,*)'jzmax,zmax1',jzmax,zmax1
      endif

      if(jzmin.gt.1) then
         zmin1=zz(jzmin)
	 rzmin1=rpsiminj(jzmin+1)
      else
c        ---------------------------------------------------
c        it is necessary to calculate the jzmin1
c        (in this case the separatrix area exist close to plasma)
c        ---------------------------------------------------
         do j=jmag-1,1,-1
cc	    write(*,*)'j,psiminj(j),psiminj(j+1)'
cc	    write(*,*)j,psiminj(j),psiminj(j+1)
            if(psiminj(j).gt.psiminj(j+1))then
	       zmin1=zz(j)
	       jzmin=j
cc	       write(*,*)'j,zmin1,jzmin',j,zmin1,jzmin
	    else
	       goto40
            endif
         enddo
 40      continue
cc	 jzmin=jzmin+1
	 jzmin=jzmin
	 zmin1=zz(jzmin)
	 rzmin1=rpsiminj(jzmin)
      endif
cc      write(*,*)'jzmin,zmin1',jzmin,zmin1
c---------------------------------------------------------------
      rmax1=rr(1)
      rmin1=rr(nxeqd)
c-------------------------------------------
c      determination of rmax1 and rmin1
c-------------------------------------------
c      write(*,*)'jzmin,jzmax,rmin1,rmax1',jzmin,jzmax,rmin1,rmax1
      do j=jzmin,jzmax
c	 write(*,*)'j, rminj(j),rmin1',j,rminj(j),rmin1
	 if (rminj(j).lt.rmin1) then
c	 write(*,*)'rminj(j),rmin1',rminj(j),rmin1
	    rmin1=rminj(j)
	    zrmin1=zz(j)
	 endif
c	 write(*,*)'rmin1',rmin1
c	 write(*,*)'j, rmaxj(j),rmax1',j,rmaxj(j),rmax1
	 if (rmaxj(j).gt.rmax1) then
	    rmax1=rmaxj(j)
	    zrmax1=zz(j)
	 endif
c	 write(*,*)'rmax1',rmax1
      enddo !j

c      write(*,*)'rmax1,rmin1,zmax1,zmin1',rmax1,rmin1,zmax1,zmin1
      return
      end


c     ****************PRESPSI********************************
c     * this function calculates plasma pressure on psi	     
c     *******************************************************
c     ! input parameter: psi				    !
      double precision
     1function prespsi(psi)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'five.i'
      double precision ias1r
c     idx derivativs order 0.ge.idx.le.3
      idx=0
cSm030224
      nx4=nx+4
      prespsi=ias1r(tpres,nx,nx4,cpres,idx,psi)
      return
      end

           
      subroutine wrtnetcdf_wall_limiter_data(netcdfnml)
 
c-------------------------------------------------------------
c     writes wall and limiter coordinates to the existing
c     netcdf file   netcdfnml
c
c     The wall and limiters data are in /one.nml/ and /fourb/ 

c     n_wall is a number of wall points
c
c     max_limiters is a number of limiters
c
c     wall and limiter (r,z) coordinates in [m]
c     r_wall(1:n_wall)   
c     z_wall(1:n_wall)  
c
c     n_limiter(1:max_limiters) number of limiter points
c
c     r_limiter(1:n_limiter_a,1:max_limiters_a)   
c     z_limiter(1:n_limiter_a,1:max_limiters_a)
c--------------------------------------------------------------
      implicit none
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'fourb.i'

      include 'netcdf.inc'
       

c-----input
      character*(*) netcdfnml      !name of output nc file
    
c-----locals----------------------------------------  
      integer i,j,n_limiter_max

c-----Storage tem1 is used in netcdf writes
      integer n1n2
      parameter (n1n2=n_limiter_a*max_limiters_a)
      real*8 tem1(n1n2)
      real*8 tem12(2*max_limiters_a)

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,n_wall_id,max_limiters_id,ncvdef2,ncdid2,
     +ncddef2
      integer n_limiter_max_id,i_two_id,two

      integer limiter_dims(2),start(2),limiter_count(2),   
     &phi_limiter_dims(2),phi_limiter_count(2)
    
      data start/1,1/

      if(myrank.ne.0) return

c-----calculate maximal value of limiter points n_limiter_max
c     at all limiters (from 1 to max_limiters)
      n_limiter_max=0 
      n_limiter_max_id=0
      max_limiters_id=0
      
      write(*,*)'in sub wrtnetcdf_wall_limiter_data'
 
      do i=1,max_limiters
          n_limiter_max=max(n_limiter_max,n_limiter(i))
         write(*,*)'wrtnetcdf_wall_limiter_data'
         write(*,*)'i,n_limiter(i),n_limiter_max',
     &              i,n_limiter(i),n_limiter_max
      enddo
      write(*,*)'n_limiter_max',n_limiter_max


      limiter_count(1)=n_limiter_max
      limiter_count(2)=max_limiters
      write(*,*)'limiter_count',limiter_count
      two=2
      phi_limiter_count(1)=two
      phi_limiter_count(2)=max_limiters
c.......................................................................
cl    1. Initialize part
c
c-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c-------------------------------------------------------------
c      create net CDF file  define dimensions, variables
c      and attributes
c------------------------------------------------------------

c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

      write(*,*)'in netcdf_wall_limiter_data netcdfnml=',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file

      write(*,*)'in netcdf_wall_limiter_data after ncopn(n netcdfnml'

      call check_err(istatus)
c,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,      
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf2(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf2(integer ncid,integer error_code)
     
c      call ncredf2(ncid,error_code)

      write(*,*)'before ncredf2(ncid,istatus)'

      call ncredf2(ncid,istatus)

      write(*,*)'before after(ncid,istatus)'

      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c
c.......................................................................
cl    1.1.2 define dimensions
c.......................................................................
c      write(*,*)'before ncddef(ncid,n_wall nwall=',n_wall

      if (n_wall.gt.0)
     &n_wall_id=ncddef2(ncid,'n_wall',n_wall,istatus)
      call check_err(istatus)

c      write(*,*)'after ncddef(ncid,n_wall,n_wall'
c      write(*,*)'ncid',ncid,'max_limiters',max_limiters

      if (max_limiters.gt.0)
     &max_limiters_id=ncddef2(ncid,'max_limiters',max_limiters,istatus)

c      write(*,*)'after ncddef(ncid,max_limiters'

      call check_err(istatus)

c      write(*,*)'brefore ncddef(ncid, n_limiter_max',n_limiter_max

      if (n_limiter_max.gt.0)
     &n_limiter_max_id=ncddef2(ncid,'n_limiter_max',n_limiter_max,
     &                        istatus)

      write(*,*)'netcdf_wall_limiter_data: n_limiter_max=',
     + n_limiter_max
      write(*,*)'n_wall',n_wall
      write(*,*)'max_limiters',max_limiters

      call check_err(istatus)

      i_two_id=ncddef2(ncid,'i_two_id',two,istatus)
      call check_err(istatus)

      limiter_dims(1)=n_limiter_max_id
      limiter_dims(2)=max_limiters_id
      write(*,*)'limiter_dims',limiter_dims
      phi_limiter_dims(1)=i_two_id
      phi_limiter_dims(2)=max_limiters_id
c.......................................................................
cl    1.1.3 define variables
c.......................................................................
      vid=ncvdef2(ncid,'n_wall',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'Number of wall points',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'max_limiters',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of limiters',istatus)
      call check_err(istatus)

     
      if (n_wall.gt.0) then

         vid=ncvdef2(ncid,'r_wall',NCDOUBLE,1,n_wall_id,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,17,
     &               'wall r coordinate',istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &              'm',istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'z_wall',NCDOUBLE,1,n_wall_id,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,17,
     &               'wall z coordinate',istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &              'm',istatus)
         call check_err(istatus)
      endif

      if(max_limiters.gt.0)then

         vid=ncvdef2(ncid,'n_limiter',NCLONG,1,max_limiters_id,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +            'Numbers of limiter points',istatus)
         call check_err(istatus)

        vid=ncvdef2(ncid,'r_limiter',NCDOUBLE,2,limiter_dims,istatus)
        call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'limiter r coordinate',istatus)
        call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'm',istatus)
        call check_err(istatus)

        vid=ncvdef2(ncid,'z_limiter',NCDOUBLE,2,limiter_dims,istatus)
        call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'limiter z coordinate',istatus)
        call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'm',istatus)
        call check_err(istatus)

       write(*,*)'before vid=ncvdef2(phi_limiter'
        vid=ncvdef2(ncid,'phi_limiter',NCDOUBLE,2,phi_limiter_dims,
     &            istatus)
        call ncaptc2(ncid,vid,'long_name',NCCHAR,31,
     +           'toroidal bounadries of limiters',istatus)
        call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'degree',istatus)
        call check_err(istatus)
        write(*,*)'after vid=ncvdef2(phi_limiter'
      endif

cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf2(ncid,istatus)
      call check_err(istatus)
          
c----------------------------------------------------------------
c     write data
c-----------------------------------------------------------------
      call ncvid2(vid,ncid,'n_wall',istatus)
      call ncvpt_int2(ncid,vid,1,1,n_wall,istatus)

      call ncvid2(vid,ncid,'max_limiters',istatus)
      call ncvpt_int2(ncid,vid,1,1,max_limiters,istatus)

      write(*,*)'n_wall',n_wall

      if (n_wall.gt.0) then
         do i=1,n_wall
            write(*,*)'i,r_wall(i)',i,r_wall(i)
         enddo

         write(*,*)'ncid',ncid 
         write(*,*)'before vid=ncvid(ncid,r_wall'
         call ncvid2(vid,ncid,'r_wall',istatus)
         call ncvpt_doubl2(ncid,vid,1,n_wall,r_wall,istatus)
         call check_err(istatus)

          write(*,*)'before vid=ncvid(ncid,z_wall'
         call ncvid2(vid,ncid,'z_wall',istatus)
         call ncvpt_doubl2(ncid,vid,1,n_wall,z_wall,istatus) 
         call check_err(istatus)
      endif

       write(*,*)'max_limiters',max_limiters

      if (max_limiters.gt.0) then

         write(*,*)'before vid=ncvid(ncid,n_limiter'
         call ncvid2(vid,ncid,'n_limiter',istatus)

         call ncvpt_int2(ncid,vid,1,max_limiters,n_limiter,istatus)
         call check_err(istatus)

         write(*,*)'before pack21r_limiter'

         call pack21(r_limiter,1,n_limiter_a,1,max_limiters_a,tem1,
     &               n_limiter_max,max_limiters)

         write(*,*)'after pack21 r_limiter'
         call ncvid2(vid,ncid,'r_limiter',istatus)
         call ncvpt_doubl2(ncid,vid,start,limiter_count,tem1,istatus)
         call check_err(istatus)

         call pack21(z_limiter,1,n_limiter_a,1,max_limiters_a,tem1,
     &               n_limiter_max,max_limiters)
         call ncvid2(vid,ncid,'z_limiter',istatus)
         call ncvpt_doubl2(ncid,vid,start,limiter_count,tem1,istatus)
         call check_err(istatus)

         write(*,*)'before pac21 phi_lim'
         call pack21(phi_limiter,1,2,1,max_limiters_a,tem12,
     &               2,max_limiters)
         write(*,*)'after pac21 phi_lim'
         call ncvid2(vid,ncid,'phi_limiter',istatus)
         call ncvpt_doubl2(ncid,vid,start,phi_limiter_count,tem12,
     +                     istatus)
         call check_err(istatus)
      endif


c---------------------------------------------------
c         4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

      return
      end       

      

     
      subroutine netcdf_eqdsk_data(netcdfnml)
c-------------------------------------------------------------
c     writes eqdsk poloidal flux array 
c     (peqd(i,j),i=1,nxeqd),j=1,nyeqd)
c     and rr(nxeqda),zz(nyeqda), mesh to the existing
c     netcdf file   netcdfnml
c
c--------------------------------------------------------------
      implicit none
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'fourb.i'
      include 'netcdf.inc'
      include 'three.i'
c-----input
      character*(*) netcdfnml      !name of output nc file
    
c-----locals----------------------------------------  
      integer i
c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,nxeqd_id,nyeqd_id,ncvdef2,ncdid2,ncddef2
      integer psi_dims(2),start(2),psi_count(3)

c-----Storage tem1 is used in netcdf writes.
      integer nxny
      parameter (nxny=nxeqda*nyeqda)
      real*8 tem1(nxny)
 
      data start/1,1/
c.......................................................................
cl    1. Initialize part
c
c-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c-------------------------------------------------------------
c      create net CDF file  define dimensions, variables
c      and attributes
c------------------------------------------------------------

c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

      if(myrank.ne.0) return

      write(*,*)'in  netcdf_eqdsk_data=',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)
c,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,      
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf2(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf2(integer ncid,integer error_code)
     
c      call ncredf2(ncid,error_code)
      call ncredf2(ncid,istatus)
      call check_err(istatus)
c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c
c.......................................................................
cl    1.1.2 define dimensions
c.......................................................................
      nxeqd_id=ncddef2(ncid,'nxeqd',nxeqd,istatus)
      nyeqd_id=ncddef2(ncid,'nyeqd',nyeqd,istatus)
      

      write(*,*)'nxeqd',nxeqd
      write(*,*)'nyeqd',nyeqd

c.......................................................................
cl    1.1.3 define variables
c.......................................................................
      vid=ncvdef2(ncid,'nxeqd',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of r points',istatus)
      call check_err(istatus)

      write(*,*)'vid=ncvdef2(ncid,nxeqd=',vid

      vid=ncvdef2(ncid,'nyeqd',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of z points',istatus)
      call check_err(istatus)

      write(*,*)'vid=ncvdef2(ncid,nyeqd=',vid

      vid=ncvdef2(ncid,'eqdsk_r',NCDOUBLE,1,nxeqd_id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     &            'eqdsk r array ',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &           'm',istatus)
      call check_err(istatus)

      write(*,*)'vid=ncvdef2(ncid,eqdsk_r=',vid

      vid=ncvdef2(ncid,'eqdsk_z',NCDOUBLE,1,nyeqd_id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     &            'eqdsk z array ',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &           'm',istatus)
      call check_err(istatus)

      write(*,*)'vid=ncvdef2(ncid,eqdsk_z=',vid

c      vid=ncvdef2(ncid,'eqdsk_psi',NCDOUBLE,nxeqd_id,nyeqd_id,istatus)
c      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
c     &            'eqdsk poloidal flux',istatus)
c      call check_err(istatus)
c      call ncaptc2(ncid,vid,'units',NCCHAR,1,
c     &           ' ',istatus)
c      call check_err(istatus)

c      write(*,*)'vid=ncvdef2(ncid,eqdsk_psi=',vid

      psi_dims(1)=nxeqd_id
      psi_dims(2)=nyeqd_id

      psi_count(1)=nxeqd
      psi_count(2)=nyeqd

      write(*,*)'psi_dims',psi_dims
      write(*,*)'psi_count',psi_count

      vid=ncvdef2(ncid,'eqdsk_psi',NCDOUBLE,2,psi_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     +         'eqdsk poloidal flux',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &           ' ',istatus)
      call check_err(istatus)

      write(*,*)'vid=ncvdef2(ncid,eqdsk_psi=',vid
      
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf2(ncid,istatus)
      call check_err(istatus)

c----------------------------------------------------------------
c     write data
c-----------------------------------------------------------------
      write(*,*)'ncid=',ncid
      call ncvid2(vid,ncid,'nxeqd',istatus)
      call ncvpt_int2(ncid,vid,1,1,nxeqd,istatus)

      call ncvid2(vid,ncid,'nyeqd',istatus)
      call ncvpt_int2(ncid,vid,1,1,nyeqd,istatus)
     
      call ncvid2(vid,ncid,'eqdsk_r',istatus)
      call ncvpt_doubl2(ncid,vid,1,nxeqd,rr,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'eqdsk_z',istatus)
      call ncvpt_doubl2(ncid,vid,1,nyeqd,zz,istatus)
      call check_err(istatus)

    
      call pack21(peqd,1,nxeqda,1,nyeqda,tem1,nxeqd,nyeqd)
      call ncvid2(vid,ncid,'eqdsk_psi',istatus)
      call ncvpt_doubl2(ncid,vid,start,psi_count,tem1,istatus)



c---------------------------------------------------
c         4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

      return
      end       

      subroutine read_density_temperature_at_rz_mesh_namelist
c      allocate(dens_r_z_in(1:nxeqd_add,1:nyeqd_add,1:nbulk),STAT=istat)
c      call bcast(dens_r_z_in,zero,SIZE( dens_r_z_in))

c      allocate(temperature_r_z_in(1:nxeqd_add,1:nyeqd_add,1:nbulk),
c     &         STAT=istat)
c      call bcast(temperature_r_z_in,zero,SIZE(temperature_r_z_in))

c      allocate(zeff_r_z_in(1:nxeqd_add,1:nyeqd_add),STAT=istat)
c      call bcast(zeff_r_z_in,1.d0,SIZE(zeff_r_z_in))


      return
      end
 
      subroutine put_density_temperature_zeff_radial_profiles_to_rz_mesh
c--------------------------------------------------------------------
c     put density ,temperature and zeff from radial profiles to arrays
c     dens_r_z_in(i,j,k),temperature_r_z_in(i,j,k),zeff_r_z_in(i,j) 
c     inside LCFS
c--------------------------------------------------------------------

      implicit none

      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'     
      include 'edge_prof_nml.i'
      include 'ions.i'
      include 'fourb.i'   

      real*8 phi,z,r

      integer i,j,k
c-----externals
      real*8 b,densrho, temperho,zeffrho 
c-----locals 
      integer
     &istat
   
      real*8 zero

c-----allocate arrays dens_r_z_in,temperature_r_z_in,zeff_r_z_in

      zero=0.d0

      allocate(dens_r_z_in(1:nxeqd_add,1:nyeqd_add,1:nbulk),STAT=istat)
      call bcast(dens_r_z_in,zero,SIZE( dens_r_z_in))

      allocate(temperature_r_z_in(1:nxeqd_add,1:nyeqd_add,1:nbulk),
     &         STAT=istat)
      call bcast(temperature_r_z_in,zero,SIZE(temperature_r_z_in))

      allocate(zeff_r_z_in(1:nxeqd_add,1:nyeqd_add),STAT=istat)
      call bcast(zeff_r_z_in,1.d0,SIZE(zeff_r_z_in))

      phi=0.d0
      do j=1,nyeqd_add           
         do i=1,nxeqd_add             
             z=zz_add(j)
             r=rr_add(i)    
             bmod=b(z,r,phi) !calculate rho
             if (rho.lt.1.d0-1.d-10) then
                do k=1,nbulk           
                   dens_r_z_in(i,j,k)= densrho(rho,k)
                   temperature_r_z_in(i,j,k)= temperho(rho,k)
                   zeff_r_z_in(i,j)= zeffrho(rho)
                enddo !k
             endif 
         enddo !i
      enddo !j
      return
      end

      subroutine
     & put_dens_temp_rz_into_dens_temp_r_z
c-----put density and temperature given at rz mesh 
c     from dens_rz_in and temperature_rz_in
c     to arrays: density_r_z, temperature_r_z
c
c     It changes number of limiters given in the namelist /tokamak/
c     for  max_limiters=0
c
c     and it
c     deallocate(dens_rz_in)
c
      implicit none

      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'     
      include 'edge_prof_nml.i'
      include 'ions.i'
      include 'fourb.i'

      integer i,j,k,m,istat

      max_limiters=0
      
      do m=0,max_limiters
         do k=1,nbulk     
            do j=1,nyeqd_add           
               do i=1,nxeqd_add  
                   density_r_z(i,j,k,m)= dens_r_z_in(i,j,k)
                   temperature_r_z(i,j,k)= temperature_r_z_in(i,j,k)
               enddo !i
            enddo!j
         enddo !k
      enddo !m

      deallocate(dens_r_z_in,STAT=istat)
      deallocate(temperature_r_z_in,STAT=istat)

      return
      end
    
      subroutine read_density_temperature_at_rz_mesh_txt
c-----It reads from input file dens_rz_in.txt 
c
c     nxeqd_add,nyeqd_add   dimensions of rz mesh
c     nbulk                 number of plasma species
c
c     density [MKS] and temperature [KeV] (i,j,k) at r(i),z(k) mesh
c                              for all plasma species k=1,nbulk
c     dens_rz_in(i,j,k)   
c     temperature_rz_in(i,j,k) 
c
c
c     Uniform rz mesh rr_add(1:nxeqd_add),zz_add(1:nyeqd_add) is determined by
c     the input parameters:
c     nxeqd_add, nxeqd_add - are numbers of mesh points. they are set in
c                input genray.dat (or genray.in) in the namelist /edge_prof_nml/
c     the lengths of the meshs are coinsided with eqdsk lengths   
c
c     k=1,...,nbulk the number of plasma species
c     Dimensions of (r,z) mesh are given by the input data
c     nxeqd_add and nyeqd_add which are set in the namelist /edge_prof_nml/  
c
c     This mesh was created in subroutine creat_fine_2D_poloidal_mesh 
c     dstep=xdimeqd/(nxeqd_add-1) 
c     do i=1,nxeqd_add
c        rr_add(i)=redeqd+dstep*(i-1)
c      enddo
c
c     dstep=ydimeqd/(nyeqd_add-1)
c     do i=1,nyeqd_add
c        zz_add(i)=ymideqd+dstep*(i-1)-ydimeqd*0.5d0
c     enddo 
c
c     It works for izeff=0 case only
c     izeff =0 zeff will be calculated using the given ions;
c             electron density will be calculated using ions;
c     So k=2,nbulk
c
c     The code will use  dens_rz_in  outside LCFS only
c     Inside LCFS density profiles versus small radius should be given and
c     density at rz mesh dens_rz_in inside LCFS should coinside with
c     given radial density profile at 
c     rr_add(1:nxeqd_add),zz_add(1:nyeqd_add) mesh points.
c
c     The code can set data for dens_rz_in(i,j,k) 
c     inside LCFS using the given radial density profiles.
c
c   
c     Using ion's density profiles at rz mesh the code will calculate
c     electron density dens_rz_in(i,j,k=1) at r-z mesh and
c     Zeff at r-z mesh zeff_rz 
      
c     Structure of the data for density in the input file.
c     The data are in 1D array prof(nxeqd_adda*nyeqd_adda*nbulka).
c     Nonzero data are in posicions from 1 to nxeqd_add*nyeqd_add*(nbulk-1)
c     There are nbulk  serias of length 1:nxeqd_add*nyeqd_add data 
c     Each seria the density for k=2,..., nbulk specie.
c
c     Inside a seria density is given at points in following order 
c     rr_add(1),rr_add(2), .., rr_add(nxeqd_add);  at zz_add(1),
c     rr_add(1),rr_add(2), .., rr_add(nxeqd_add);  at zz_add(2)
c........................................................
c     rr_add(1),rr_add(2), .., rr_add(nxeqd_add);  at zz_add((nyeqd_add)

      implicit none

      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'     
      include 'edge_prof_nml.i'
      include 'ions.i'
      include 'fourb.i'

c-----output 
   
c-----locals
      
      real*8 zero
      integer i,j,k,istat,i_unit_loc,
     &nxeqd_add_file,     !dimension R mesh used in input file
     &nyeqd_add_file,   !dimensions Z meshused in input file
     &nbulk_file          !number of plasma species used in input file

      integer kode          !  < 0   the end of input file was detected
                            !  =0     reading has complited succefully
                            !  >0     an error has occurred 

      zero=0.d0

      allocate(dens_r_z_in(1:nxeqd_add,1:nyeqd_add,1:nbulk),STAT=istat)
      call bcast(dens_r_z_in,zero,SIZE( dens_r_z_in))

      allocate(temperature_r_z_in(1:nxeqd_add,1:nyeqd_add,1:nbulk),
     &         STAT=istat)
      call bcast(temperature_r_z_in,zero,SIZE(temperature_r_z_in))

      allocate(zeff_r_z_in(1:nxeqd_add,1:nyeqd_add),STAT=istat)
      call bcast(zeff_r_z_in,1.d0,SIZE(zeff_r_z_in))
     
      i_unit_loc=2 
      rewind i_unit_loc

      open(i_unit_loc,file='dens_temp_rz_in.txt',status='old',
     &     iostat=kode)
  
      if (kode.ne.0) then
         write(*,*)'dens_temp_rz_in.txt is not present'
         stop
      endif

      write(*,*)'read_density_temperature_at_rz_mesh_txt'
      write(*,*)'before read (i_unit_loc,2)'

2     format(52x,3i5) 
      read(i_unit_loc,2) nxeqd_add_file,nyeqd_add_file,nbulk_file  !mesh dimensions

      if(nbulk_file.ne.nbulk) then
        write(*,*)'in dens_temp_rz_in.txt nbulk_file.ne.nbulk'
        write(*,*)'nbulk_file=',nbulk_file,'nbulk=',nbulk
        write(*,*)'Please change nbulk in genray.in or genray.dat'
        stop 'read_density_temperature_at_rz_mesh_txt'
      endif

      if(nxeqd_add_file.ne.nxeqd_add) then
        write(*,*)'in dens_temp_rz_in.txt nxeqd_add_file.ne.nxeqd_add'
        write(*,*)'nxeqd_add_file=',nxeqd_add_file
        write(*,*)'nxeqd_add=',nxeqd_add
        write(*,*)'Please change nxeqd_add in genray.in or genray.dat'
        stop 'read_density_temperature_at_rz_mesh_txt'
      endif

      if(nyeqd_add_file.ne.nyeqd_add) then
        write(*,*)'in dens_temp_rz_in.txt nyeqd_add_file.ne.nyeqd_add'
        write(*,*)'nyeqd_add_file=',nyeqd_add_file
        write(*,*)'nyeqd_add=',nyeqd_add
        write(*,*)'Please change nyeqd_add in genray.in or genray.dat'
        stop 'read_density_temperature_at_rz_mesh_txt'
      endif

3     format(5e16.9)

c-----electron density can be arbitry, it will be recalculated
c     bellow

      read(i_unit_loc,3)(((dens_r_z_in(i,j,k),
     &                   i=1,nxeqd_add),j=1,nyeqd_add),k=1,nbulk)

      read(i_unit_loc,3)(((temperature_r_z_in(i,j,k),
     &                   i=1,nxeqd_add),j=1,nyeqd_add),k=1,nbulk)
    
      
      close(i_unit_loc)

      do k=1,nbulk
        do j=1,nyeqd_add           
           do i=1,nxeqd_add             
             dens_r_z_in(i,j,k)=dens_r_z_in(i,j,k)*1.d-19   !in 1.d13/cm**3

c             write(*,*)'read k,j,i,dens_r_z_in(i,j,k)',
c     &                       k,j,i,dens_r_z_in(i,j,k)               
           enddo
        enddo
      enddo
c---------------------------------------------------------------------
c     calculate electron density and zeff from ions densities
c---------------------------------------------------------------------
      do j=1,nyeqd_add           
         do i=1,nxeqd_add
             dens_r_z_in(i,j,1)=0.d0   
             zeff_r_z_in(i,j)=0.d0          
             do k=2,nbulk           
                dens_r_z_in(i,j,1)= dens_r_z_in(i,j,1)+
     &                          charge(k)*dens_r_z_in(i,j,k) 
                zeff_r_z_in(i,j)= zeff_r_z_in(i,j)+
     &                          charge(k)**2*dens_r_z_in(i,j,k) 
             enddo
 
             zeff_r_z_in(i,j)=zeff_r_z_in(i,j)/ dens_r_z_in(i,j,1)
         enddo !i
      enddo !j

c      k=1
      do j=1,nyeqd_add           
         do i=1,nxeqd_add
            zeff_r_z(i,j)=zeff_r_z_in(i,j)
c            write(*,*)'read k,j,i,dens_r_z_in(i,j,1),zeff_r_z_in(i,j)',
c     &                      k,j,i,dens_r_z_in(i,j,1),zeff_r_z_in(i,j)
         enddo !i
      enddo !j

c      stop 'read_write*t.f in  read_density_temperature_at_rz_mesh'

      return
      end

      subroutine writes_dens_temp_psi_rho_rz_out_txt
c-----writes
c
c     to the filedens_temp_rz_out.txt
c     nxeqd_add,nyeqd_add  dimensions of rz mesh
c     nbulk                 number of plasma species
c     density and temperature (i,j,k) at r(i),z(k) mesh
c                              for all plasma species k=1,nbulk
c                                  
c     
c     to the file psi_rho_rz_out.txt:
c
c     nxeqd_add,nyeqd_add, !dimensions of rz mesh
c
c     some data from eqdsk file
c     xdimeqd,ydimeqd,     !horizontal and vertical full-widths 
c                          !of the rectangle [meters] 
c     reqd                 !the nominal major radius of the torus.
c     redeqd,              !major radius of the inner edge of rectangular grid.
c     ymideqd,             !vertical shift of the rectangular box
c
c     r mesh               !r(1,nxeqd_add)
c     z-mesh               !z(1,nyeqd_add
c     poloidal flux psi(i,j) at r(i),z(j) nmesh
c     small radius rho (i,j) at r(i),z(j) mesh  
c
      implicit none
      include 'param.i'   
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'    
      include 'edge_prof_nml.i' 
      include 'three.i'
      include 'fourb.i'
c-----locals
 
      integer i,j,k,istat,i_unit,kode,nbulk_loc

      real*8 zero

      real*8, dimension (1:nxeqd_add,1:nyeqd_add,1:nbulk) ::
     &                   dens_r_z_in_MKS,temp_prof
      real*8, dimension (1:nxeqd_add,1:nyeqd_add) ::
     &                   psi_r_z,rho_r_z

c-----external    
      real*8 psif,rhof,b,temperho



c     &nxeqd_add,nyeqd_add, !dimensions of rz mesh
c     &xdimeqd,ydimeqd,     !horizontal and vertical full-widths
c                           ! of the rectangle [meters]
c     reqd                  !the nominal major radius of the torus.
c     &redeqd,              !major radius of the inner edge of rectangular grid.
c     &ymideqd,             !vertical shift of the rectangular box
c                           ! up-down simmetry plane
   

c      write(*,*)'in writes_dens_temp_psi_rho_rz_out_txt'
c      do k=1,nbulk
c         do j=1,nyeqd_add
c            do i=1,nxeqd_add
c               write(*,*)'01k,j,i,density_r_z(i,j,k,0)',
c     &                      k,j,i,density_r_z(i,j,k,0)
c            enddo
c         enddo
c      enddo

      do k=1,nbulk
        do j=1,nyeqd_add           
           do i=1,nxeqd_add
cSAP110905
             bmod=b(zz_add(j),rr_add(i),0.d0) !calculate small radius
             
             dens_r_z_in_MKS(i,j,k) = density_r_z(i,j,k,0)*1.d19 !MKS units [1/m**3]                 
c             write(*,*)'wr k,j,i,density_r_z(i,j,k,0)*1.d19',
c     &                     k,j,i,density_r_z(i,j,k,0)*1.d19
             
             temp_prof(i,j,k) = temperho(rho,k) !Kev
c             write(*,*)'r.j,i,rho,prof_temp_prof(i,j,k)',
c     &                  k,j,i,rho,prof_temp_prof(i,j,k)
           enddo
        enddo
c        stop 'read_write writes_dens_temp_psi_rho_rz_out_dat'
      enddo
c----------------------------------------------------------------------------
c       put temperature poloidal flux and small radius rho at rz mesh 
c-----------------------------------------------------------------------------   
c      write(*,*)'read_write_input.f writes_dens_temp_psi_rho_rz_out_dat'
      do j=1,nyeqd_add      
c        write(*,*)'j',j     
        do i=1,nxeqd_add
         
c           write(*,*)'i,zz_add(i),rr_add(j)',
c     &                i,,zz_add(i),rr_add(j)
cSAP110905
           bmod=b(zz_add(j),rr_add(i),0.d0) !calculate small radius rho
c           write(*,*)'bmod',bmod
cSAP110905
c           prof_psi(iprof) = psif(zz_add(i),rr_add(j))
            psi_r_z(i,j) = psif(zz_add(j),rr_add(i))
c           write(*,*)' prof_psi(iprof)',prof_psi(iprof)

            rho_r_z(i,j) = rho
c           write(*,*)'from i,j,b rho',i,j,rho
           
        enddo !i
      enddo !j
      
      if(myrank.ne.0) return     

c      stop 'in writes_dens_temp_psi_rho_rz_out_txt'

c      write(*,dens_temp_rz)

      i_unit=1
      rewind i_unit

      open(i_unit,file = 'dens_temp_rz_out.txt',iostat=kode) !MKS [1/m**3]
 2    format(52x,3i5) 
      write(i_unit,2) nxeqd_add,nyeqd_add,nbulk  !mesh dimensions
      call check_param(i_unit)

 3    format(5e16.9)
      write(i_unit,3)(((dens_r_z_in_MKS(i,j,k),
     &                i=1,nxeqd_add),j=1,nyeqd_add),k=1,nbulk)     !MKS
      write(i_unit,3)(((temp_prof(i,j,k),
     &                i=1,nxeqd_add),j=1,nyeqd_add),k=1,nbulk)     !KeV
      close(i_unit)     
    
      
      i_unit=2

      open(i_unit,file = 'psi_rho_rz_out.txt',iostat=kode) 
 4    format(52x,2i5) 
      write(i_unit,4) nxeqd_add,nyeqd_add  !mesh dimensions
      write(i_unit,3) xdimeqd,ydimeqd,reqd,redeqd,ymideqd
      write(i_unit,3)(rr_add(i),i=1,nxeqd_add)
      write(i_unit,3)(zz_add(j),j=1,nyeqd_add)
      write(i_unit,3)((rho_r_z(i,j),i=1,nxeqd_add),j=1,nyeqd_add)      
      write(i_unit,3)((psi_r_z(i,j),i=1,nxeqd_add),j=1,nyeqd_add)
       
      close(i_unit)

      return
      end


