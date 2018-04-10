      subroutine adjint(ientry)
      use kind_spec

      implicit none

c-----input
      integer ientry
      include 'param.i'
      include 'adj.i'
c-----locals
      integer thpp0,icurdr,iswchi
      integer :: nw, nthp, ns, js, jp, jps ,nthpp0
cSm070625
      integer iadj
      real(kind=dp) :: zef1
      integer nmax,imax,lmax

      select case (ientry)
      case default

c----------------
c  assume nthp0,npsi0,psimx,psimn already known
c
      icurdr=1
  
c      write(*,*)'adj_control.f in  adjint ientry=',ientry

c      write(*,*)'adj_control.f before op_file 0 icurdr=',icurdr

      call op_file(0,icurdr)
  
c      write(*,*)'adj_control.f adjint before op_file 3 icurdr=',icurdr
  
      call op_file(3,icurdr)
c      write(*,*)'adj_control.f adjint before op_file 4 icurdr=',icurdr
      call op_file(4,icurdr)
c      write(*,*)'adj_control.f adjint before op_file 5 icurdr=',icurdr
      call op_file(5,icurdr)
c      write(*,*)'adj_control.f adjint before op_file 6 icurdr=',icurdr
      call op_file(6,icurdr)
c      call op_file(7,icurdr)
c      call op_file(8,icurdr)

      write(*,*)'adj_control.f in  sub adjint'
      write(*,*)'iout3,iout4,iout5,iout6,iout7,iout8',
     &iout3,iout4,iout5,iout6,iout7,iout8

      nthpp0=nthp0
c--------------------------------------------------------------
c     calculate
c     &thetas(nthpp0+1),   !poloidal angle mesh for b line integration 
c     &dla(0:nthpp0-1),    !length along b line 
c     &ba(0:nthpp0-1),     !b/b0   along b line
c     &thetae(0:nthpp0-1)  !poloidal angle mesh for b line integration 
c---------------------------------------------------------------------------   
cSm070625
c      write(*,*)'adj_control.f adjint before bavorbt nthp0,npsi0',
c     & nthp0,npsi0

      call bavorbt( nthp0,npsi0)

c      write(*,*)'adj_control.f adjint after bavorbt'

c      write(*,*)'adj_control.f adjint after bavorbt1'
c      write(*,*)'in sub adjint thetas',thetas
c      write(*,*)'in sub adjint thetae',thetae

      write(*,*)'in adj_control in adjint before end case 1' 
      write(*,*)'adj_control.f in  sub adjint'
      write(*,*)'iout3,iout4,iout5,iout6,iout7,iout8',
     &iout3,iout4,iout5,iout6,iout7,iout8
     

      return
      case (2)
c----------------
c     get Spitzer-function table
c     chi(theta,ua)=chi(0:imax1 - 1,0:nmax - 1) and output to adjout
c     At each radial point (np=1,npsi0) write chi and other variables 
c     to the file with file number: iout5 and file name: 'adjout 
        
        iswchi=0

        nmax=nmax_chi
        imax=imax_chi
        lmax=lmax_chi 

        write(*,*)'adj_control.f in adjint subadj nmax,imax,lmax',
     &             nmax,imax,lmax
        write(*,*)'umax,nthp0,dt',umax,nthp0,dt
        write(*,*)'tmax,alpha,ze,t,rho_',tmax,alpha,ze,t,rho_
        write(*,*)'iswchi,aerrmx,rerrmx',iswchi,aerrmx,rerrmx
        write(*,*)' iout3, iout5, iout4, iout7, iout8',
     &              iout3, iout5, iout4, iout7, iout8

        write(*,*)'adj_control.f in sub adjint before subadj'

        call subadj(nmax,umax,imax,nthp0,dt,tmax,alpha,ze,t,
     &   rho_, lmax, iout3, iout5, iout4, iout7, iout8, iswchi, aerrmx
     &   , rerrmx,
     &  npsi0,dens_averaged_ar)

        write(*,*)'adj_control.f adjint after subadj'
        write(*,*)'dens_averaged_ar',dens_averaged_ar


        return
      case (3)

c  get j/p table
c
cSm070625
c         iadj = 1
c         call curadj (nthp0, npsi0, psimx, psimn, nmax, umax1, imax,
c     1      lmax, psis(1), eden(1), etem(1), zefi(1), sjop(1), nwp, wpmx
c     2      , wpmn, dwp, wp(1), nthp1, thetp(1), iout4, iout5, iadj,
c     3      gwk2, gwk3, npnta)
c         write (iout6, 1000) jopid
c         write (iout6, 1010) nwp, nthp0, npsi0, nthp1
c         write (iout6, 1020) psimx, psimn, wpmx, wpmn, dwp
c         write (iout6, 1030) (wp(nw),nw=1,nwp)
c         write (iout6, 1030) (thetp(nthp),nthp=1,nthp1)
c         write (iout6, 1030) (psis(ns),ns=1,npsi0)
c         write (iout6, 1030) (eden(ns),ns=1,npsi0)
c         write (iout6, 1030) (etem(ns),ns=1,npsi0)
c         write (iout6, 1030) (zefi(ns),ns=1,npsi0)
c         do js = 1, npsi0
c            do jp = 1, nthp1
c               jps = (js - 1)*nthp1*nwp + (jp - 1)*nwp + 1
c               write (iout6, 1030) (sjop(nw),nw=jps,jps + nwp - 1)
c            end do
c         end do
c--------------------------------
 1000    format(1x,1a40)
 1010    format(1x,4i8)
 1020    format(5(1x,e18.10))
 1030    format(5(1x,e18.10))
c------------------------
         return
      end select
      end subroutine adjint

      subroutine adj_chi_function
c------------------------------------------------------------------
c     calculate adj chi function using Karney subroutins from Curray code
c------------------------------------------------------------------- 
c     open output files, calculate
c     &thetas(nthpp0+1),   !poloidal angle mesh for b line integration 
c     &dla(0:nthpp0-1),    !length along b line 
c     &ba(0:nthpp0-1),     !b/b0   along b line
c     &thetae(0:nthpp0-1)  !poloidal angle mesh for b line integration 
c---------------------------------------------------------------------------
c      write(*,*)'adj_chi_function before call adjint(1)'
      call adjint(1)
c      write(*,*)'adj_chi_function after  call adjint(1)'
c----------------------------------------------------------------------------
c     get Spitzer-function table
c     chi(theta,ua)=chi(0:imax1 - 1,0:nmax - 1) and output to adjout
c     At each radial point (np=1,npsi0) write chi and other variables 
c     to the file with file number: iout5 and file name: 'adjout
c----------------------------------------------------------------------------
c      write(*,*)'adj_chi_function before call adjint(2)'
      call adjint(2)
c      write(*,*)'adj_chi_function after  call adjint(2)'
c--------------------------------------------------------------------------
      return
      end

     
