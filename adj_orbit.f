

c     This module contains orbit-averaging routines
c--------------------------------------

c      npsi0 input 


      subroutine bavorb1(thetas, dla, ba, thetae,nthpp0)
     
      implicit none
      include 'param.i'
      include 'three.i'
      include 'one.i'
      include 'adj.i' ! has npsi0, nthp0
c-----input
      integer nthpp0   !the number of poloidal points for the integration

      real*8, dimension(nthpp0 + 1) :: thetas
      real*8, dimension(0:nthpp0 - 1) :: dla, ba, thetae
      
c      real*8 
c     &thetas(nthpp0+1),   !poloidal angle mesh for b line integration 
c     &dla(0:nthpp0-1),    !length along b line 
c     &ba(0:nthpp0-1),     !b/b0   along b line
c     &thetae(0:nthpp0-1)  !poloidal angle mesh for b line integration 
   
c--------------------------------------------------------------
c     Externals
c--------------------------------------------------------------     
      real*8 b,bmax_psi,bmin_psi, 
     &rhopsi,densrho,temperho,zeffrho,psi_rho

C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer nth, n, n_psi, imn, imx, nt, it, ier

      real*8 dthet,beq0,alenb,
     &bmax,bmin,
     &dpsis,psi,theta,phi,btot,bpol,ymx,ymn,z,r,rho_psis,
     &rho_in,
cSAP080204
     &ba_min,ba_max,drho,rho_max
         
cSAP080211 for comparison of DC conductivity with CQL
      integer npsi0_loc
      parameter (npsi0_loc=11)   
      real(8) rya(npsi0_loc)
c---------------------------------------------------------
c     output for writing
      real*8   psimx, psimn, 
     &dene,teme,zi

c-----for nonuniform ADJ radial mesh
      real*8 rho_L,rho_R,drho_L,drho_R

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
        
c-----small radius array from CQL for the comparison of DC conductivity
      data rya / 
     &0.0454545468091965, 0.136363640156659, 0.227272732962262, 
     &0.318181825226005, 0.409090916947889, 0.500000008127912, 
     &0.590909098766075, 0.681818188862378, 0.772727278416821, 
     &0.863636367429404, 0.954545455900127 /

c  given nthp0,npsi0,psimx,psimn
c  find dl(1:nthp0)=dl/L line elements along field line
c    and bb(1:nthp0)=B(l)/B(0)
c  (for the bounce average integrals in adj)
c  typically npsi0=26,nthp0=60
c
      write(*,*)'in bavorb1 nthpp0',nthpp0

      pi = 4.d0*datan(1.0d0)
      psimx=psilim   ! GR from three.i
      psimn=psimag   ! GR from three.i

c      write(*,*)'from CQL rya',rya 
cSAP080205  for special case: circular eqdsk rho_max=R_axis
      psimx=(psimx-psimn)*0.8d0 +psimn  !here was used 0.8 to reduse max small radius
                                        !The code had not convergency at eps>0.8
                                        !for test circular eqdsk which gave 0<eps<1     
c      rho_max=0.8d0
      rho_max=1.d0                      !for test case it was rho_max=0.8 from 1.
      drho=rho_max/(npsi0 - 1)          !uniform small radii mesh
      write(*,*)'0 drho,npsi0',drho,npsi0
      dpsis = (psimx - psimn)/(npsi0 - 1) !npsi0  uniform poloidal flux mesh
                                          !The original Karney code used uniform
                                          !poloidal flux mesh  

       dthet = 2*pi/nthp0

      write(*,*)'bavorb1  nthp0, npsi0', nthp0, npsi0

      if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
      write (iout3, 2000) nthp0,npsi0 !in bavorb1 ! need opened(iout3,
      ! out3 is for 'adjinp' file
      endif  !On myrank=0   ! myrank=0

      write(*,*)'bavorb1  psimx,c psimn', psimx, psimn

      if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
      write (iout3, 2010) psimx, psimn !in bavorb1
      endif  !On myrank=0   ! myrank=0

      if(nthpp0<nthp0) then
        if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
        WRITE(*,*)' bavorb1: nthpp0<nthp0' 
        WRITE(*,*)'it should be nthpp0 >= nthp0'
        WRITE(*,*)'nthpp0, nthp0',nthpp0, nthp0
        endif  !On myrank=0   ! myrank=0
        STOP
      endif

      do nth = 1, nthp0 + 1
         thetas(nth) = (nth - 1.25d0)*dthet
c         write(*,*)'nth,thetas(nth)',nth,thetas(nth)
      end do

      do nth = 0, nthp0 - 1
         thetae(nth) = (nth + 0.25d0)*dthet
c         write(*,*)'nth,thetae(nth)',nth,thetae(nth)
      end do

cSm070830
c      thetae_1(0:nthp0_a-1)=thetae(0:nthp0_a-1)
      do nth = 0, nthp0 - 1
         thetae_1(nth)=thetae(nth)
      enddo

      do nth = 1, nthp0 + 1
         write(*,*)'nth,thetas(nth)',nth,thetas(nth)
      enddo

      if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
      write (iout3, 2020) (thetas(n),n=1,nthp0 + 1) !in bavorb1
      endif  !On myrank=0   ! myrank=0

      do nth = 0, nthp0 - 1
         write(*,*)'nth,thetae(nth)',nth,thetae(nth)
      enddo

      if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
      write (iout3, 2020) (thetae(n),n=0,nthp0 - 1) !in bavorb1
      endif  !On myrank=0   ! myrank=0


c      write(*,*)' bavorb1 npsi0', npsi0

      do n_psi = 1, npsi0
        if (npsi0.eq.1) then
            rho_in=0.5d0                !small radius
c            temp_kev=temperho(rho_in,1) !electron temperatura
c            unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
c            c2=(2.99792458d10/unorm)**2 !(clight/unorm)**2
c            vnorm_cql=umax*unorm
c            eps=rho_in
            psis(n_psi)=psi_rho(rho_in)
        else
cSAP080206 
            psis(n_psi) = psimn + dpsis*(n_psi - 1) !uniform poloidal flux mesh used
                                                     !in the original Karney code
            rho_in=(n_psi-1)*drho                    !uniform small radial mesk
c------------------------------------------------------------------------
cSAP081114   for radial mesh localized at rho intervalxdraw =(0., 0.25)
c            if(n_psi.eq.npsi0) then
c               rho_in=(n_psi-1)*drho 
c            else  
c               rho_in=(n_psi-1)*drho/4.d0
c            endif   
 
cSAP0801208  
c            rho_R=0.08d0
c            rho_L=0.03d0
c            drho_R=(1.d0-rho_R)/10.d0
c            drho_L=(rho_R-rho_L)/(npsi0-10-1)
c            if (n_psi.le.npsi0-10) then
cc              rho_in=rho_L+(n_psi-1)*drho_L
cc-------------at npsi=npsi0-10
cc             rho_in=rho_L+(npsi0-10-1)*(rho_R-rho_L)/(npsi0-10-1)=rho_R
c            else
cc-------------npsi=npsi0-10+1,...,npsi0
c              rho_in=rho_R+(n_psi-(npsi0-10))*drho_R
cc------------at npsi=npsi0-10+1      
cc             rho_in=rho_R+((npsi0-10+1) -(npsi0-10))*drho_R=rho_R+drho_R
cc------------at npsi=npsi0
cc             rho_in=rho_R+(npsi0 -(npsi0-10))*drho_R=rho_R+(10*drho_R)=rho_R+10(1-rho_R)/10=1
c            endif
c            drho=0.04d0/(npsi0-1)
c            write(*,*)'n_psi,drho',n_psi,drho
c            if(n_psi.eq.npsi0) then
c               rho_in=1.d0-1.d-6
c              write(*,*)'n_psi.eq.npsi0 n_psi,rho_in', n_psi,rho_in
c            else  
c               rho_in=0.03d0+(n_psi-1)*drho
c            endif     
cSAP090307    for FW NSTX case CD was located at rho 0.02-0.08
c            and npsi0=10
c            if (n_psi.eq.2) rho_in=0.02d0
c            if ((n_psi.ge.3).and.(n_psi.le.8))then
c               rho_in=0.02d0+(n_psi-2)*0.01d0
c            endif
c            if (n_psi.gt.8) rho_in=0.08d0+(n_psi-8)*(1.d0-0.08d0)*0.5d0   
    
            if (rho_in.gt.1.d0) then
               if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
               WRITE(*,*)'adj_orbit.f  subroutine bavorb1 rho_in.gt.1'
               WRITE(*,*)'n_psi=',n_psi,'rho_in=',rho_in
               WRITE(*,*)'it should be for adj calculations rho_in.le.1'
               WRITE(*,*)'Please change setting of small radius of '
               WRITE(*,*)'of the adj flux surface rho_in in'
               WRITE(*,*)' subroutine bavorb1 and recompile the code'
               endif  !On myrank=0   ! myrank=0
               STOP 'in adj_orbit.f'
            endif

c-----------------------------------------------------------------------
            psis(n_psi) =psi_rho(rho_in)
            write(*,*)'adj_orbit n_psi,rho_in,psis(n_psi)',
     &                           n_psi,rho_in,psis(n_psi)
        endif
        psi = psis(n_psi)

c        write(*,*)'adj_orbit.f in  bavorb1 n_psi,psis(n_psi)',
c     &            n_psi,psis(n_psi)

c----------------------------------------------------------------
c       calculate density temperature,zef at given psi
c       call subpar (psi, eden(np), etem(np), zefi(np))
c        if (psi.eq.1.d0) psi = 1. - 1.Ed-5
        if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
        write (iout3, 2010) psi !in bavorb1
        endif  !On myrank=0   ! myrank=0

        write(*,*)'n_psi,psi',n_psi,psi
        write(*,*)' bavorb1 nthp0',nthp0

        Q_safety_adj(n_psi)=0.d0
    
        do nth = 1, nthp0 + 1
           theta = thetas(nth)
c----------------------------------------------------------------
c          calculate Z,R coordinates versus poloidal angle 
c          at the given flux surface
c----------------------------------------------------------------
           call zr_psith(psi,theta,z,r)

c           write(*,*)'n_psi,nth,psi,theta,z,r ',
c     &                n_psi,nth,psi,theta,z,r 

           cxc(nth,n_psi) = r  !array of major radius at given small radius
           czc(nth,n_psi) = z  !array of Z at given small radius

c           write(*,*)'nth,n_psi,cxc(nth,n_psi),czc(nth,n_psi)',
c     &                nth,n_psi,cxc(nth,n_psi),czc(nth,n_psi)

        end do  !nth 

        theta = 0.d0
        call zr_psith(psi,theta,z,r) !calculate (r,z) at given (psi,theta=0)
        btot=b(z,r,phi)
        beq0 = dabs(btot)            !B_total at given (psi,theta=0)
  
c        write(*,*)'B_total at given (psi,theta=0) beq0',beq0

c--------------------------------------------------------------------------
c       calculate arrays ba=B_total/B_tolal(theta=0) 
c       dla=steps of of B_field line length
c       at the poloidal mesh
c       at the given poloidal flux (small radius)
c----------------------------------------------------------------------- 
        Q_safety_adj(1)=1.d0
        do nth = 0, nthp0 - 1
            theta = thetae(nth)
            call zr_psith(psi,theta,z,r)
            btot=b(z,r,phi)
cSm080204
c            ba(nth) = abs(btot/beq0)
            ba(nth) = abs(btot)          !total B non divided by b_min                    

            bpol=dsqrt(br**2+bz**2)
            dla(nth) = sqrt((cxc(nth+2,n_psi)-cxc(nth+1,n_psi))**2+
     1                (czc(nth+2,n_psi)-czc(nth+1,n_psi))**2)*btot/bpol
          
c           write(*,*)'n_psi,nth,theta,dla(nth)',n_psi,nth,theta,dla(nth)
c           write(*,*)'btot,bpol,btot/bpol ',btot,bpol,btot/bpol

c            write(*,*)'del_pol', 
c     &                dsqrt((cxc(nth+2,n_psi)-cxc(nth+1,n_psi))**2+
c     &                (czc(nth+2,n_psi)-czc(nth+1,n_psi))**2)     
            if(n_psi.gt.1) then
              Q_safety_adj(n_psi)=Q_safety_adj(n_psi)+(bphi/btot)/
     &        (pi*(cxc(nth+2,n_psi)+cxc(nth+1,n_psi)))*dla(nth) 
            endif      
        end do
c-----------------------------------------------------------------------
c       calculate 
c       alenb =of B_field line length
c       at the given poloidal flux (small radius)
c----------------------------------------------------------------------
        alenb = 0.d0
        do nth=0,nthp0-1
           alenb =alenb+dla(nth)
        enddo
c        write(*,*)'n_psi,alenb',n_psi,alenb
c-----------------------------------------------------------------------
c       normalization of  
c       dla/alenb=(steps of of B_field line length/length of B_filed line)
c       at the poloidal mesh
c       at the given poloidal flux (small radius)
c-----------------------------------------------------------------------
        do nth=0,nthp0-1
           dla(nth)=dla(nth)/alenb
c           write(*,*)'n_psi,nth,dla(nth)',n_psi,nth,dla(nth)
        enddo
        
C
C       determine bmax,bmin at given poloidal flux psi 
C       using functions bmax_psi and bmin_psi
C       
        bmax= bmax_psi(psi)
        bmin= bmin_psi(psi)	 ! (in Tl)  
c        write(*,*)'n_psi,psi,bmax,bmin',n_psi,psi,bmax,bmin

c---------------------------------------------------------------------------
c       normalization bmax and bmin by beq0(B_total at given (psi,theta=0)
c---------------------------------------------------------------------------
        ymx=bmax/beq0            
        ymn=bmin/beq0

c---------------------------------------------------------------------------
c       ba=b_total(length)/b_min
c--------------------------------------------------------------------------
cSAP 080204
c-------calculate max and min total magnetic field at the flux surface
        ba_min=bmax          !obtained from function bmax_psi
        ba_max=bmin          !obtained from function bmin_psi

c-------compare bm_min and_ba_max with total B values in array ba()
c       and recalculation ba_min and ba_max 
        do nth=0,nthp0-1
          if (ba_min.gt.ba(nth)) ba_min=ba(nth)
          if (ba_max.lt.ba(nth)) ba_max=ba(nth)
        enddo 
        if (bmax.lt.ba_max) then
           write(*,*)'adj_orbit.f in sub bavorb1 bmax.lt.ba_max'
           write(*,*)'bmax,ba_max',bmax,ba_max
           bmax=ba_max+1.d-6
        endif
        if (bmin.gt.ba_min) then
           write(*,*)'adj_orbit.f in sub bavorb1 bmin.gt.ba_min'
           write(*,*)'bmin,ba_min',bmin,ba_min
           bmin=ba_min-1.d-6
        endif
 
c-------calculate ba=b/bmin
        do nth=0,nthp0-1
          ba(nth) = ba(nth)/bmin 
          write(*,*)'n_psi,nth,ba(nth),bmax/bmin',
     &               n_psi,nth,ba(nth),bmax/bmin

        enddo 

c        do nth=0,nthp0-1
c          ba(nth) = ba(nth)/ymn 
c          write(*,*)'n_psi,nth,ba(nth),bmax/bmin',
c     &               n_psi,nth,ba(nth),bmax/bmin       
c        enddo   

c---------------------------------------------------------------------------
c       integral{0,l}]((b_total/bmin)*dl)
c----------------------------------------------------------------------------
        bbar(n_psi)=0.d0
        do nth=0,nthp0-1
            bbar(n_psi) = bbar(n_psi) + ba(nth)*dla(nth)
c            write(*,*)'n_psi,nth,bbar(n_psi)',n_psi,nth,bbar(n_psi)
        enddo

        if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
        write (iout3, 2020) (dla(n),n=0,nthp0 - 1) !in bavorb1
        endif  !On myrank=0   ! myrank=0

        do n=0,nthp0-1
           write(*,*)'n,dla(n)',n,dla(n)
        enddo

        if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
        write (iout3, 2020) (ba(n),n=0,nthp0 - 1) !in bavorb1
        endif  !On myrank=0   ! myrank=0

        do n=0,nthp0-1
           write(*,*)'n,ba(n)',n,ba(n)
        enddo
cSm070830
c       ba_2d(0:nthp0_a-1,n_psi)=ba(0:nthp0-1)
        do nth=0,nthp0-1
           ba_2d(nth,n_psi)=ba(nth)
        enddo
c        write(*,*)'adj_orbit.f in bavorb1n_psi=',n_psi

c        do nth=0,nthp0-1
c          write(*,*)'nth,ba_2d(nth,n_psi)',nth,ba_2d(nth,n_psi)
c        enddo 
  
        if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
        write (iout3, 2010) alenb, bmax, bmin !in bavorb1
        endif  !On myrank=0   ! myrank=0

c        write(*,*)'alenb, bmax, bmin',alenb, bmax, bmin

        rho_psis=rhopsi(psis(n_psi)) !   normalized small radius at given psi GR
        dene=densrho(rho_psis,1)  !   electron density in 10**13 cm**-3
        teme=temperho(rho_psis,1) !   electron temperature in KeV
        zi=zeffrho(rho_psis)      !   zeff 

      enddo !n_psi
cSAP 070727
      if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
       close(iout3) !in bavorb1
      endif  !On myrank=0   ! myrank=0
c
c---------------------------
 2000 format(1x,2i8)
 2010 format(4(2x,e20.12))
 2020 format(5(2x,e20.12))
c 2010 format(4(2x,e24.16))
c 2020 format(5(2x,e24.16))

c      'end of  bavorb1'
      return
      end !subroutine bavorb1


      subroutine op_file(i_unit,icurdr)
      
      implicit none 

c-----input
      integer i_unit,icurdr
      include 'param.i'
      include 'adj_no_nml.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      integer iout,iounit1,iout0,iout2,
     & iout1

      write(*,*)'in op_file i_unit,icurdr', i_unit,icurdr
c
c
      select case (i_unit)
      case default
      !YuP/note: the following is done when i_unit=0
      iout = 6
      iounit1 = 7
      iout0 = 9
      iout2 = 12
      iout3 = 13
      iout4 = 14
      iout5 = 15
      iout6 = 16
      iout7 = 17
      iout8 = 18
      iout1 = 21

cSm070626
      iout = iout+200
      iounit1 =  iounit1+200  
      iout0 = iout0 +200
      iout2 = iout2 +200
      iout3 = iout3 +200
      iout4 = iout4 +200
      iout5 = iout5 +200
      iout6 = iout6 +200
      iout7 = iout7 +200
      iout8 = iout8 +200
      iout1 = iout1 +200

      write(*,*)'adj_orbit.f in sub op_file'
      write(*,*)'iout',iout
      write(*,*)'iounit1',iounit1  
      write(*,*)'iout0' , iout0 
      write(*,*)'iout2' , iout2 
      write(*,*)'iout3' , iout3 
      write(*,*)'iout4' , iout4 
      write(*,*)'iout5' , iout5 
      write(*,*)'iout6' , iout6 
      write(*,*)'iout7' , iout7 
      write(*,*)'iout8' , iout8 
      write(*,*)'iout1' , iout1
c

      case(1)
c         iout = iounit(i_unit) ! op_file is never called with (1)
         if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         open(iout, file='currayout', status='unknown')
         endif  !On myrank=0   ! myrank=0
c         write(*,*)'open file currayout'
      case(2)
c         iout2 = iounit(i_unit) ! op_file is never called with (2)    
         if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         open(iout2, file='rayop', status='unknown')
         endif  !On myrank=0   ! myrank=0
c         write(*,*)'open file rayop'
      case(3)
c         iout3 = iounit(i_unit)
         if (icurdr.eq.1) then
            if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
            open(iout3, file='adjinp', status='unknown') ! op_file
            endif  !On myrank=0   ! myrank=0
         else
            if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
            open(iout3, file='adjinp', status='old') ! op_file
            endif  !On myrank=0   ! myrank=0
         endif
c         write(*,*)'open file adjinp'
      case(4)
c         iout4 = iounit(i_unit)
         if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         open(iout4, file='scrach', status='unknown')
         endif  !On myrank=0   ! myrank=0
c         write(*,*)'open file scrach'
      case(5)
c         iout5 = iounit(i_unit)
         if (icurdr.eq.1) then
            if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
            open(iout5, file='adjout', status='unknown') ! in op_file
            endif  !On myrank=0   ! myrank=0
         else
            if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
            open(iout5, file='adjout', status='old') ! in op_file
            endif  !On myrank=0   ! myrank=0
         endif
c         write(*,*)'open file adjout'
      case(6)
c         iout6 = iounit(i_unit)
         if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
         open(iout6, file='joptab', status='unknown')
         endif  !On myrank=0   ! myrank=0
c         write(*,*)'open file joptab' 
      end select
      return
      end subroutine op_file

c      subroutine bavorbt( nthp00,nthp0,npsis,npsi0,iout4)
      subroutine bavorbt(nthp0,npsi0)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  11:01:19   6/14/01
C...Switches: -p4 -yb
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------

      use kind_spec

      implicit none

c-------------------------------------------------
c     input
c      integer nthp00,nthp0,npsis,npsi0,iout4
      integer nthp0,npsi0
C-----------------------------------------------
C   G l o b a l   P a r a m e t e r s
C-----------------------------------------------
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer:: pthets,pdla,pba,pthete,ndum
      integer:: pwk,pwk1,pwkt,ii,ierr,i,nthpp0
      real(kind=dp), allocatable, dimension(:) :: dum
      real(kind=dp) :: stat
C-----------------------------------------------
c
c  parameters
c  Nspect-->nspectm, Nion-->nions,  Nua-->nuam,  Nta-->ntam, Lrz-->lrza,
c  Nthp11-->nthp11a,  Npnt-->npnta,   Mpnt-->mpnta
c
c----------------------------
c       common block Comiorlh(1-8)
c
c  ipropt         integer   #option for n-T profile calculation: 1 if
c                        # calculated from prontp; else =0.


c--------SCC  5/23/95
c*******  TKM  06/97, 2/00
c*************************
c********    TKM    11/98
c*******************************
c  #   SCC 3/27/92   SCC 4/13/95
c---------------------------
c   common block for Comraylh
c-------------------------
c
c
c$$$$$$$$     TKM    11/98
c*********    TKM   2/97, 6/99, 9/99
c*********    TKM   2/97
c*******   TK Mau   8/98
c*************************
c  SCC  3/7/95

c  Storage related to EQDSK:
c  Units used in the EQDSK are entirely MKSA.
c  Comeqdsk-->ceqdsk.i

c
c   Common block for Comdynray
c
c
c  given nthp0,npsi0,psimx,psimn
c  find dl(1:nthp0)=dl/L line elements along field line
c    and bb(1:nthp0)=B(l)/B(0)
c  (for the bounce average integrals in adj)
c  typically npsi0=26,nthp0=60
c
c     pointer(dumptr,dum(1))
c     pointer(pthets,thetas(1:nthp0+1))
c     pointer(pdla,dla(0:nthp0-1))
c     pointer(pba,ba(0:nthp0-1))
c     pointer(pthete,thetae(0:nthp0-1))
c     pointer(pwk,wk(1:nthp0+4))
c     pointer(pwk1,wk1(1:nthp0+4,3))
c     pointer(pwkt,wkt(1:nthp0+4))
c
c  check if nthp0 is smaller than nthp00 and npsi0 is smaller than npsis
c

c      if (npsi0>npsis .or. nthp0>nthp00-1) then
c      if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
c         write (iout4, *) 'npsi0 or nthp0 is too large in bavorbt'
c         write (iout4, *) npsi0, nthp0
c         write (iout4, *) npsis, nthp00
c      endif  !On myrank=0   ! myrank=0
c         !stop
c      endif
      ii = 1
      pthets = ii
      ii = ii + nthp0 + 1
      pdla = ii
      ii = ii + nthp0
      pba = ii
      ii = ii + nthp0
      pthete = ii
c      ii = ii + nthp0
c      pwk = ii
c      ii = ii + nthp0 + 4
c      pwk1 = ii
c      ii = ii + (nthp0 + 4)*3
c      pwkt = ii
      ndum = ii + nthp0 + 4 - 1

c      call hpalloc(dumptr,ndum,ierr,1)
      if (.not.allocated(dum)) then
         ALLOCATE(dum(ndum),stat=ierr)
         if (ierr /= 0) then
            print *, 'dum could not be allocated in bavorbt'
            return
         endif
      endif
      dum(:ndum) = 0.
      nthpp0 = nthp0

      write(*,*)'in bavorb nthp0,nthpp0,npsi0',nthp0,nthpp0,npsi0

      call bavorb1 (dum(pthets), dum(pdla), dum(pba), dum(pthete),
c     1                 dum(pwk), dum(pwk1), dum(pwkt), nthpp0)
     &                 nthpp0)
c      call hpdeallc(dumptr,ierr,0)
      DEALLOCATE(dum)
      return
      end subroutine bavorbt
