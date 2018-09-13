c*************************gr2new***************************************
c  It calculates coordinates of the contours for the flux functions:
c          r(psi,teta)   z(psi,teta)     			      *
c          arrays rpsi(j,i) zpsi(j,i)  				      *
c          j=1,npsi(number of counturs poloidal flux=constant)	      *
c          i=1,nteta+1(number of point in poloidal angle)    	      *
c  and creates arrays ar(nl,nteta+1),az(nl,nteta+1) for nl contours   *
c  plotting       						      *
c---------------------------------------------------------------------
c  input data are in common blocks /three/,/five/,/gr.cb/	      *
c  nl-is the number of flux contours for plotting
c---------------------------------------------------------------------*
c* Output: arrays AR,AZ, zpsi,rpsi are into common block 'gr          *
c*         (file gr.cb)                                               *
c**********************************************************************
c
      subroutine gr2new0 ! not called
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'gr.i'
      include 'three.i'
      include 'five.i'
      include 'one.i'
      integer*2 it1,it2,it3,it4,itf1,itf2,itf3,itf4

      pi=4.d0*datan(1.d0)
      hteta=2.d0*pi/dble(nteta)
      hpsi=(psilim-psimag)/dble(npsi-1)
c      write(*,*)'gr2new.f: psilim,psimag,npsi,hpsi',
c     1           psilim,psimag,npsi,hpsi
      do 1 j=1,npsi
	 arpsi(j)=psimag+hpsi*(j-1)
c         wite(*,*)'gr2new j,arpsi(j)',j,arpsi(j)
1     continue
      do 2 i=1,nteta1
 	 arteta(i)=hteta*(dble(i)-0.5d0)
2     continue
c----------------------------------------------------------
c      if ipsi=1 then continue,ipsi=0 then read file psi.bin
cbegin
c       ipsi=1 !  to calculate contours
c       ipsi=0 !  to read contours data from file:psi.bin

       if (ipsi.eq.0) then
        open(1,file='psi.bin',form='unformatted',status='old')
	do i=1,npsi
	  do j=1,nteta1
             read(1)zpsi(i,j),rpsi(i,j)
	  end do
	end do
	close(1)
	go to 200
       endif
c------------------------------------------------------------------
       do 10 i=1,nteta
         teta=arteta(i)
         sintet=dsin(teta)
         costet=dcos(teta)
         call zrlacner(sintet,costet,yma,xma,zmax,rmax,rmin,zmin,tm)
         tr=tm
         do 20 j=npsi,2,-1
           tl=0.d0
	   c=arpsi(j)
c          -----------------------------------------------
c          determination of z(j,i)=z(arpsi(j),arteta(i))
c          and r(j,i)=r(arpsi(j),arteta(i)) using the binary method
c          -----------------------------------------------
           do while ((tr-tl).gt.epspsi)
             t=tl+(tr-tl)*0.5d0
             r=xma+t*costet
             z=yma+t*sintet
             psi1=fpsi(r,z)-c
             rtr=xma+tr*costet
             ztr=yma+tr*sintet
             psi2=fpsi(rtr,ztr)-c
             if ((psi1*psi2).gt.0) then
               tr=t
             else
               tl=t
             end if
           end do
c          -----------------------------------------------
c          end of the binary method
c          -----------------------------------------------
c	   write(*,*)'j,i',j,i
c	   write(*,*)'tr-tl,t',tr-tl,t
c	   write(*,*)'psi1,fpsi(r,z),c',psi1,fpsi(r,z),c
c	   write(*,*)'r,z',r,z
           rpsi(j,i)=r
           zpsi(j,i)=z
20       continue
10     continue
       do 30 i=1,nteta
         zpsi(1,i)=yma
         rpsi(1,i)=xma
30     continue

       do 40 j=1,npsi
           zpsi(j,nteta+1)=zpsi(j,1)
           rpsi(j,nteta+1)=rpsi(j,1)
40     continue
c----------------------------------------------------------
c      if ipsi=1 then continue,write file psi.bin
        open(1,file='psi.bin',form='unformatted')
	do i=1,npsi
	  do j=1,nteta1
             write(1)zpsi(i,j),rpsi(i,j)
	  end do
	end do
	close(1)
c------------------------------------------------------------
c      if ipsi=1 then continue,ipsi=0 then read file psi.bin
cend
200    continue
c-------------------------------------------------------------
c       do j=1,npsi
c         do  i=1,nteta1
c 	   zpsi1(j,i)=zpsi(j,i)
c 	   rpsi1(j,i)=rpsi(j,i)
c	 enddo
c       enddo
c       do 100 j=1,npsi
c         do 100 i=1,nteta1
c 	  write(*,*)'j,arpsi(j),i,arteta(i)',j,arpsi(j),i,arteta(i)
c 	  write(*,*)'zpsi,zpsi1,rpsi,rpsi1'
c 	  write(*,*)zpsi(j,i),zpsi1(j,i),rpsi(j,i),rpsi1(j,i)
c	  write(*,*)'j,i,zpsi,rpsi',j,i,zpsi(j,i),rpsi(j,i)
c 100   continue
c----------------------------------------------------
c  the arrays for contours ploting
       jstep=npsi/NL

       do 50 j=1,NL

         do 60 i=1,nteta+1
c   multiplication to get length in cm and denormalization
c	   AZ(j,i)=zpsi(jstep*j,i)*100.d0
c	   AR(j,i)=rpsi(jstep*j,i)*100.d0
c           write(*,*)'j,i AR(j,i),AZ(j,i)',j,i,AR(j,i),AZ(j,i)
	   AZ(j,i)=zpsi(jstep*j,i)*100.d0*r0x
	   AR(j,i)=rpsi(jstep*j,i)*100.d0*r0x
60	 continue
50     continue

      return
      end

c--------------------------------------------------------------------
c***********************zrlacner*****************************************
c      it calculates the paramerar t for which
c        (z(t),r(t)) are the coordinates of the point in which
c      line: z=zmag+sin(teta)*t
c            z=rmag+cos(teta)*t
c      intersects the Lackner rectangle
c----------------------------------------------------------------------
c      input data:sintet=sin(teta),costet=cos(teta)
c                       teta-poloidal angle
c                       It must be not equal 0;0.5pi;pi;1.5pi;2pi
c                  zmag,rmag -magnetic axis coordinates
c                  zmax,rmax,zmin,rmin -boundaries of Lackner rectangle
c-----------------------------------------------------------------------
c      output data:t -parameter (length along the line)
c      these data are used in     to find the counturs
c-----------------------------------------------------------------------
      subroutine zrlacner(sintet,costet,zmag,rmag,zmax,rmax,rmin,zmin,
     1                    t)
      implicit double precision (a-h,o-z)
c-------------------------------------------
      zmaxd=zmax
      zmind=zmin
      rmaxd=rmax
      rmind=rmin
      zmax=zmaxd+0.01d0
      rmax=rmaxd+0.01d0
      zmin=zmind-0.01d0
      rmin=rmind-0.01d0

      write(*,*)'in zrlacner zmax,rmax,zmin,rmin',zmax,rmax,zmin,rmin

      write(*,*)'in zrlacner sintet,costet,zmag,rmag',
     &                       sintet,costet,zmag,rmag

c----------------------------------------------------------
      if ((costet.gt.0.d0).and.(sintet.gt.0.d0)) then
c     1 quarter of the (r,z) plate
        tz=(zmax-zmag)/sintet
        tr=(rmax-rmag)/costet
	go to 1
      end if
c----------------------------------------------------------
      if ((costet.lt.0.d0).and.(sintet.gt.0.d0)) then
c     2 quarter of the (r,z) plate
        tz=(zmax-zmag)/sintet
        tr=(rmin-rmag)/costet
	go to 1
      end if
c----------------------------------------------------------
      if ((costet.lt.0.d0).and.(sintet.lt.0.d0)) then
c     3 quarter of the (r,z) plate
        tz=(zmin-zmag)/sintet
        tr=(rmin-rmag)/costet
	go to 1
      end if
c----------------------------------------------------------
      if ((costet.gt.0.d0).and.(sintet.lt.0.d0)) then
c     4 quarter of the (r,z) plate
        tz=(zmin-zmag)/sintet
        tr=(rmax-rmag)/costet
	go to 1
      end if
c----------------------------------------------------------
1       continue
        t=dmin1(tz,tr)
c----------------------------------------------------------
      rmax=rmaxd
      rmin=rmind
      zmax=zmaxd
      zmin=zmind

      write(*,*)'in zrlacner tz,tr,tm',tz,tr,tm

      return
      end

c*************************gr2new***************************************
c  It calculates contours coordinates for the flux functions:		
c          r(psi,teta)   z(psi,teta)     			      *
c          arrays rpsi(j,i) zpsi(j,i)  				      *
c          j=1,npsi(the number of the countours poloidal_flux=constant)	      
c          i=1,nteta+1(number of the point in the poloidal angle)    	      
c  and creates arrays ar(nl,nteta+1),az(nl,nteta+1) for nl contours   *
c  plotting       						      *
c---------------------------------------------------------------------
c  input data are in common blocks /three/,/five/,/gr.cb/	      *
c  nl-the number of flux contours for plotting
c---------------------------------------------------------------------*
c  Output: arrays AR,AZ, zpsi,rpsi into common block 'gr              *
c          (file gr.cb)                                               *
c**********************************************************************
c
      subroutine gr2new
      implicit none
c      implicit real*8 (a-h,o-z)
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'gr.i'
      include 'three.i'
      include 'five.i'
      include 'one.i' 

      integer*2 it1,it2,it3,it4,itf1,itf2,itf3,itf4

c       call gettim(it1,it2,it3,it4)
c       tstart=it1*3600.00d0+it2*60.00d0+it3*1.00d0+it4*0.01d0
c-----local 
      real*8 r_psi_1d(nteta1),z_psi_1d(nteta1),psi_loc,   
     &theta_pol_limit(nteta1),rho_geom_limit(nteta1)

c for test of field line and binary methods
      real*8 rpsi_field(npsi,nteta1),zpsi_field(npsi,nteta1),p,
     & diff_r_psi(npsi),diff_z_psi(npsi),diff_r_total,diff_z_total,
     & costet,sintet,hteta,hpsi,htini,htmin,epspsit,pm,pp,psi0,
     & rmaxp,rminp,zmaxp,zminp,rt0,zt0,t0,teta,tini,tm

      integer i,j,ierr,maxiter,jstep

c-----externals
      real*8 fpsi

      pi=4.d0*datan(1.d0)
      hteta=2.d0*pi/dble(nteta)
      if (nnlim.eq.0) then
c        In this case:
c        2)we will create the limiter points using the close flux
c        surface psi(r,z)=psilim*psifactr, here
c	 psifactr is a parameter (it must be .le.1) to avoide the
c        problems with the nonmonotonic psi function near the separatrix.
c        psifactr is given in genray.in file (it is in common/one/
c        3) the number of limiter points nlimit let equal to the number
c        of points on the magnetic surface nlimit=nteta	(see common
c        blocks: gr.cb, fourb.cb  and five
c        ------------------------------------
         write(*,*)'gr2new nnlim=0 eqdsk data without limiter points'
cSAP091030
         write(*,*)'gr2new psilim,psifactr,psimag',
     &                     psilim,psifactr,psimag

	 psilim=psimag+(psilim-psimag)*psifactr
c	 psilim=psisep
         write(*,*)'psisep,psilim,psifactr',psisep,psilim,psifactr
      endif
c-------------------------------------------------------
      hpsi=(psilim-psimag)/dble(npsi-1)

      write(*,*)'gr2new.f psilim,psimag',psilim,psimag

      do j=1,npsi
	 arpsi(j)=psimag+hpsi*(j-1)
c	 write(*,*)'in gr2new j,arpsi(j)',j,arpsi(j)
      enddo
      do i=1,nteta1
 	 arteta(i)=hteta*(dble(i)-0.5d0)
c	 write(*,*)'in gr2new i,arteta(i)',i,arteta(i)
      enddo
c----------------------------------------------------------
c     if ipsi=1 then continue,ipsi=0 then read file psi.bin
cbegin
      write(*,*)'in gr2new ipsi',ipsi
c        ipsi=1 !  to calculate contours
c        ipsi=0 !  to read contours data from file:psi.bin
      if (ipsi.eq.0) then
         !if(myrank.eq.0) then  ! MPI ! YuP[2018-09-10] added
           open(1,file='psi.bin',form='unformatted',status='old') ! to read
           do i=1,npsi
           do j=1,nteta1
            read(1)zpsi(i,j),rpsi(i,j)            
            write(*,*)'read i,j',i,j
            write(*,*)'zpsi(i,j),rpsi(i,j)',zpsi(i,j),rpsi(i,j)
           end do
           end do
	     close(1)
	   !endif ! myrank=0
         go to 200
      endif ! ipsi
c------------------------------------------------------------------      
      do  i=1,nteta1
         zpsi(1,i)=yma
         rpsi(1,i)=xma
      end do

cSAP091123
      goto 5

cSAP091104
c----------------------------------------------------------------
c     Calculate arrays rpsi(j,i),zpsi(j,i) at given (psi,theta)
c     using field lines.
c     It works for all nnlim values
c     (for eqdsk with and without limiter points). 
c     It was introduced in the version from gr091116
c---------------------------------------------------------------
      do j=2,npsi
c         write(*,*)'gr2new before calc_z_r_at_psi j',j
         psi_loc=arpsi(j)
         call calc_z_r_at_psi_at_field_line(xma,yma,psi_loc,arteta,
     &   r_psi_1d,z_psi_1d)
c         write(*,*)'gr2new after calc_z_r_at_psi j',j

         do i=1,nteta
            rpsi(j,i)=r_psi_1d(i)
            zpsi(j,i)=z_psi_1d(i)
         enddo
         rpsi(j,nteta1)= rpsi(j,1)           
         zpsi(j,nteta1)= zpsi(j,1)
      enddo

c      goto 70
c-----for  test two methods
      do j=1,npsi
           do i=1,nteta1
               rpsi_field(j,i)=rpsi(j,i)
               zpsi_field(j,i)=zpsi(j,i)
           enddo
      enddo

c-------------------------------------------------------------
c     Calculate arrays rpsi(j,i),zpsi(j,i) at given (psi,theta)
c     using the binary method for the solution of the equation
c     psi(r(t),z(t))=psi0 at given poloidal angle
c
c     This approach was used in the old version (before 091106)
c     Now following lines (until lable 70)
c     can be used for comparison of two methods:
c     1) version 091106: field line soltion
c     2) old version   :  binary method
c       
c------------------------------------------------------- 
cSAP091123
 5    continue
cc-----z,r verus poloidal angle at LCFS
c      psi_loc=arpsi(npsi)
c      j=npsi
c      call calc_z_r_at_psi_at_field_line(xma,yma,psi0,arteta,
c     &r_psi_1d,z_psi_1d)
c      do i=1,nteta
c         rpsi(j,i)=r_psi_1d(i)
c         zpsi(j,i)=z_psi_1d(i)
c      enddo
c      rpsi(j,nteta1)= rpsi(j,1)           
c      zpsi(j,nteta1)= zpsi(j,1)

cSAP101004
c      write(*,*)'gr2new before do 10 nteta',nteta
         
      do 10 i=1,nteta

c        teta=hteta*(dble(i)-0.5d0)
         teta=arteta(i)

         sintet=dsin(teta)
         costet=dcos(teta)
	 pp=1.d0
	 pm=1.d0
	 rmaxp=rmax*pp
	 rminp=rmin*pm
	 zmaxp=zmax*pp
	 zminp=zmin*pm
cSAP091030
c         write(*,*)'gr2new i,teta.sintet,costet',i,teta,sintet,costet
c         write(*,*)'gr2new zmaxp,rmaxp,rminp,zminp,tm',
c     &                     zmaxp,rmaxp,rminp,zminp,tm

cSAP091214
c         call zrlacner(sintet,costet,yma,xma,
c     1	 zmaxp,rmaxp,rminp,zminp,tm)

c         write(*,*)'gr2new after zrlacner zmaxp,rmaxp,rminp,zminp,tm',
c     &                     zmaxp,rmaxp,rminp,zminp,tm

c	 tini=0.d0
c	 htini=tm*0.02d0
c         maxiter=10

	 epspsit=epspsi

c	 htmin=1.d-5
cSAP091127
         do 20 j=2,npsi 
c         do 20 j=2,npsi-1
	   psi0=arpsi(j)

c 	   write(*,*)'gr2new bef zrcontor i,teta,j,psi0',i,teta,j,psi0
c           write(*,*)'gr2newj,psi0,psilim,psi0',j,psi0,psilim
c           write(*,*)'gr2new i,teta.sintet,costet',i,teta,sintet,costet
         
c------------------------------------------------
c          Newton method for solution of equation psi(r(t),z(t))=psi0
cc         call zrcontor(tini,psi0,costet,sintet,tm,htini,
cc     1   epspsit,maxiter,htmin,ierr,zt0,rt0,t0)
c------------------------------------------------
c          binary method for solution of equation psi(r(t),z(t))=psi0
c           write(*,*)'in gr2new before zrcntrbin tm,tini',tm,tini
c           call zrcntrbin(tini,psi0,costet,sintet,tm,htini,
c     1     ierr,zt0,rt0,t0)

cSAP101004
c          write(*,*)'i,j',i,j
c          if((j.eq.19).and.(i.eq.20)) then
cfor test only
c            write(*,*)'before calc_r_z_psi_theta_binary_test'
c            write(*,*)'psi0,yma,xma,teta',psi0,yma,xma,teta
c            call calc_r_z_psi_theta_binary_test(psi0,
c     &      yma,xma,rt0,zt0,teta)
c          endif

           call calc_r_z_psi_theta_binary(psi0,
     &     yma,xma,rt0,zt0,teta)
 
c          write(*,*)'in gr2new after zrcntrbin zt0,rt0,t0',zt0,rt0,t0
 
           zpsi(j,i)=zt0
           rpsi(j,i)=rt0
          
20       continue
10    continue
cSAP091030
c       stop 'in gr2new after 10'

      do 40 j=1,npsi
           zpsi(j,nteta+1)=zpsi(j,1)
           rpsi(j,nteta+1)=rpsi(j,1)
40    continue

cSAP091123
      goto 6
c-----test methods
      
      diff_r_total=0.d0
      diff_z_total=0.d0

      do j=1,npsi
          diff_z_psi(j)=0.d0
          diff_r_psi(j)=0.d0

          do i=1,nteta1

            p=dabs(fpsi(rpsi_field(j,i),zpsi_field(j,i))-arpsi(j)) 
            write(*,*)'j,i,abs(fpsi(r(j,i),z(j,i))-arpsi(j))',j,i,p

            p=dabs(rpsi_field(j,i)-rpsi(j,i))

            if (diff_r_psi(j).lt.p) then
                diff_r_psi(j)=p
            endif

            p=dabs(zpsi_field(j,i)-zpsi(j,i))

            if (diff_z_psi(j).lt.p) then
                diff_z_psi(j)=p
            endif
          enddo

          if (diff_r_total.lt. diff_r_psi(j)) then
             diff_r_total=diff_r_psi(j) 
          endif

          if (diff_z_total.lt. diff_z_psi(j)) then
             diff_z_total=diff_z_psi(j) 
          endif

          write(*,*)'j,diff_r_psi(j),diff_z_psi(j)',
     &               j,diff_r_psi(j),diff_z_psi(j)
      enddo

      write(*,*)'diff_r_total,diff_z_total',
     &            diff_r_total,diff_z_total

c      stop 'gr2new test two methods'
c     end of the calculations rpsi(j,i),zpsi(j,i) usin binary method          
c---------------------------------------------------------------
cSAP091123
 6    continue

cSAP091104
 70   continue
c      do j=1,npsi
c       do j=npsi,npsi
c	  do i=1,nteta1
c             zpsi(j,i)=zpsi_field(j,i)
c             rpsi(j,i)=rpsi_field(j,i)
c	  end do
c      end do
c----------------------------------------------------------
c     if ipsi=1 then continue,write file psi.bin

      if(myrank.eq.0)then
        open(1,file='psi.bin',form='unformatted') ! to write
        do i=1,npsi
        do j=1,nteta1
             write(1)zpsi(i,j),rpsi(i,j)
        end do
        end do
        close(1)
      endif ! myrank=0
c------------------------------------------------------------
c     if ipsi=1 then continue,ipsi=0 then read file psi.bin
200   continue
c-------------------------------------------------------------
c     do 100 j=1,npsi
c	j=npsi
c	write(*,*)'in gr2new test'
c        do 100 i=1,nteta1
c 	  write(*,*)'j,arpsi(j),i,arteta(i)',j,arpsi(j),i,arteta(i)
c 	  write(*,*)'zpsi,rpsi',zpsi(j,i),rpsi(j,i)
c	  read(*,*)
 100  continue
c----------------------------------------------------
c  the arrays for contours ploting
      jstep=npsi/NL
      do 50 j=1,NL
         do 60 i=1,nteta+1
c   multiplication to get length in cm and denormalization
c	   AZ(j,i)=zpsi(jstep*j,i)*100.d0
c	   AR(j,i)=rpsi(jstep*j,i)*100.d0
	   AZ(j,i)=zpsi(jstep*j,i)*100.d0*r0x
	   AR(j,i)=rpsi(jstep*j,i)*100.d0*r0x
60	 continue
50    continue
c------------------------------------------------------
c       call gettim(itf1,itf2,itf3,itf4)
c       tfinish=itf1*3600.00d0+itf2*60.00d0+itf3*1.00d0+itf4*0.01d0
c       timen=tfinish-tstart
c       write(*,*)'time cpu for gr2new contour calculations',timen
      return
      end


c*************************zrcontor*********************************** *
c  It calculates  contours coordinates of the countor point 	       *
c   r(psi,teta)   z(psi,teta) for the given:			       *
c   flux function psi(z,r)=psi0 and poloidal angle teta0(in radians)   *     			      *
c---------------------------------------------------------------------
c  input data are in common blocks /three/,/five/,/gr.cb/	       *
c  tini -initial (left) value of the t( ray parameter )		       *
c  psi0 -given value of the poloidal flux psi			       *
c  costet0,sintet0 -for the given poloidal angle		       *
c  tm   -maximal value of t (right)					       *
c  htini   ininial step of the t 	 			       *
c  maxiter - maximal number of iteratios for the equation solution     *
c  epspsit-accuracy of equation solution (=max(abs( t_n-t_n+1))        *
c  epspsi-accuracy of equation solution (=max(abs( psi_n-psi_n+1) in param.i)     *
c  htmin - the minimal parameter t step
c----------------------------------------------------------------------*
c  Output: ,zt0(psi0,teta0),rt0(psi0,teta0) and parameter t0           *
c           rt0=xma+t0*costet0 ,  zt0=yma+t0*sintet0		       *
c  ierr -index if the solution was obtained =1	else=0		       *
c**********************************************************************
      subroutine zrcontor(tini,psi0,costet0,sintet0,tm,htini,
     1 epspsit,maxiter,htmin,ierr,zt0,rt0,t0)
      implicit real*8 (a-h,o-z)
      include 'param.i'
      include 'three.i'
      ierr=1
      ht=htini
c      write(*,*)'in zrcontor psi0,costet0,sintet0',
c     1 psi0,costet0,sintet0
c      write(*,*)'tini,tm,htini,epspsit,epspsi'
c     1 ,tini,tm,htini,epspsit,epspsi
c      write(*,*)'maxiter,htmin',maxiter,htmin
 10   continue
      t1=tini
      rt1=xma+t1*costet0
      zt1=yma+t1*sintet0
      psi1=fpsi(rt1,zt1)
c      write(*,*)'10 t1,rt1,zt1,psi1',t1,rt1,zt1,psi1
 20   t2=t1+ht
c      write(*,*)'20 t2,t1,rt1,zt1,psi1',t2,t1,rt1,zt1,psi1
      t2=dmin1(t2,tm)
c     write(*,*)'t2',t2
      if((t2.eq.tm.and.t1.eq.tm).or.
     1   (t2.eq.0.d0.and.t1.eq.0.d0)) then
	 ierr=0
c        error exit
         goto50
      endif

      rt2=xma+t2*costet0
      zt2=yma+t2*sintet0
      psi2=fpsi(rt2,zt2)
      dpsi=(psi0-psi1)*(psi0-psi2)
c      write(*,*)'20 t2,rt2,zt2,psi2,dpsi',t2,rt2,zt2,psi2,dpsi
      if(dpsi.le.0.d0)go to 30
      t1=t2
      psi1=psi2
      goto 20
c-----------------------------------------------------------
c     psi0 is between psi1 and psi2
c     iteration methode
c-----------------------------------------------------------
 30   numbiter=0
      tn=t1+ht*0.5d0
      rtn=xma+tn*costet0
      ztn=yma+tn*sintet0
c      write(*,*)'30 tn,rtn,ztn',tn,rtn,ztn
 40   continue
      psin=fpsi(rtn,ztn)
      call dpsidzdr(z,r,phi,dpsidz,dpsidr)
c      write(*,*)'40 tn,rtn,ztn,psin',tn,rtn,ztn,psin
      dzdt=sintet0
      drdt=costet0
      dfpsidt=dpsidz*dzdt+dpsidr*drdt
c      write(*,*)'40 dpsidz,dpsidr,dfdpidt',dpsidz,dpsidr,dfpsidt
      if (dfpsidt.eq.0.d0) then
c         write(*,*)'in zrcontur dfdpsidt=0 psi0=',psi0,
c     1    'costet0=',costet0
         stop
      endif
      dpsi=psin-psi0
      delt=-dpsi/dfpsidt
      tn=tn+delt
      rtn=xma+tn*costet0
      ztn=yma+tn*sintet0
cc only for test      psin=fpsi(rtn,ztn)
      if(dabs(delt).lt.epspsit) then
c        write(*,*)'1 delt,epspsit',delt,epspsit
c	   write(*,*)'delt,tn',delt,tn
c	   write(*,*)'psin-psi0,psin,psi0',psin-psi0,psin,psi0
c	   write(*,*)'rtn,ztn',rtn,ztn
        go to 50
      endif
      if(psi0.eq.0.and.dabs(dpsi).lt.epspsi) then
c      	write(*,*)'2 psi0,dpsi,epspsi',psi0,dpsi,epspsi
        go to 50
      endif
c      if(dabs(psin-psi0).lt.epspsi)then
c        write(*,*)'4 psin-psi0',psin-psi0
c	   write(*,*)'delt,tn',delt,tn
c	   write(*,*)'psin-psi0,psin,psi0',psin-psi0,psin,psi0
c	   write(*,*)'rtn,ztn',rtn,ztn
c        go to 50
c      endif
      if(psi0.ne.0.)  then
        if(dabs(dpsi/psi0).lt.epspsi)then
c	  write(*,*)'3 dpsi/psi0,epspsi',dpsi/psi0,epspsi
c          write(*,*)'psin-psi0',psin-psi0
	  go to 50
	endif
      endif
c------------------------------------------------
      numbiter=numbiter+1
      if(numbiter.ge.maxiter)then
         if(ht.le.htmin) stop 'zrcontur'
	 ht=0.5d0*ht
	 goto 10
      endif
      goto 40
c------------------------------------------------
c     end of iterations
c------------------------------------------------
 50   continue
      t0=tn
      zt0=ztn
      rt0=rtn
      return
      end

c*************************zrcntrbin*********************************** *
c  It calculates  the contours coordinates of the countour points      *
c   r(psi,teta)   z(psi,teta) for the given:			       *
c   flux function psi(z,r)=psi0 and poloidal angle teta0(in radians)   *
c   It using the binary method                                        *
c---------------------------------------------------------------------
c  input data are in common blocks /three/,/five/,/gr.cb/	       *
c  tini -initial (left) value of the t( ray parameter )		       *
c  psi0 -given value of the poloidal flux psi			       *
c  costet0,sintet0 -for the given poloidal angle		       *
c  tm   -maximal value of t (right)				       *
c  htini   ininial step of the t 	 			       *
c  epspsi-accuracy of equation solution (=max(abs( psi_n-psi_n+1))(in param.i)
c----------------------------------------------------------------------*
c  Output: ,zt0(psi0,teta0),rt0(psi0,teta0) and parameter t0           *
c           rt0=xma+t0*costet0 ,  zt0=yma+t0*sintet0		       *
c  ierr -index if the solution was obtained =1	else=0		       *
c**********************************************************************
      subroutine zrcntrbin(tini,psi0,costet0,sintet0,tm,htini,
     1 ierr,zt0,rt0,t0)
      implicit real*8 (a-h,o-z)
      include 'param.i'
      include 'three.i'
      ierr=1
      ht=htini
      write(*,*)'in zrcontor !1 psi0,costet0,sintet0',
     1 psi0,costet0,sintet0
      write(*,*)'tini,tm,htini,epspsi'
     1 ,tini,tm,htini,epspsi

      t1=tini
      rt1=xma+t1*costet0
      zt1=yma+t1*sintet0
      psi1=fpsi(rt1,zt1)
 10   t2=t1+ht

      write(*,*)'10 t2,t1,rt1,zt1,psi1',t2,t1,rt1,zt1,psi1

      t2=dmin1(t2,tm)

cSAP091030
      write(*,*)'t2.tm',t2,tm

      if((t2.eq.tm.and.t1.eq.tm).or.
     1   (t2.eq.0.d0.and.t1.eq.0.d0)) then
	 ierr=0
	 write(*,*)'in zrcntrbin error t2,tm,t1',t2,tm,t1
	 write(*,*)'in zrcntrbin error exit ierr',ierr
c        error exit
         goto30
      endif

      rt2=xma+t2*costet0
      zt2=yma+t2*sintet0
      psi2=fpsi(rt2,zt2)
      dpsi=(psi0-psi1)*(psi0-psi2)
c      write(*,*)'10 t2,rt2,zt2,psi2,dpsi',t2,rt2,zt2,psi2,dpsi
      if(dpsi.le.0.d0)go to 20
      t1=t2
      psi1=psi2
      goto 10
c-----------------------------------------------------------
c     psi0 is between psi1 and psi2
c     binary iteration method
c-----------------------------------------------------------
 20   continue
      tr=t2
      tl=t1
      do while ((tr-tl).gt.epspsi)
         t=tl+(tr-tl)*0.5d0
         r=xma+t*costet0
         z=yma+t*sintet0
         psi1=fpsi(r,z)-psi0
         rtr=xma+tr*costet0
         ztr=yma+tr*sintet0
         psi2=fpsi(rtr,ztr)-psi0
         if ((psi1*psi2).gt.0) then
            tr=t
         else
            tl=t
         end if
      end do
c     -----------------------------------------------
c          end of the binary method
c     -----------------------------------------------
      t0=t
      zt0=z
      rt0=r
 30   continue
      write(*,*)'the end of zrcotrbin t0',t0
      return
      end


c        **********************dpsidzdr***********************
c        *                        -                           *
c        * this subroutine calculates the derivatives         *
c        * output:dpsidz and dpsidr           		      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c------------------------------------------------------------------
      subroutine dpsidzdr(z,r,phi,dpsidz,dpsidr)
      implicit double precision (a-h,o-z)

      include 'param.i'
      include 'five.i'

c-----externals      
      real*8 ias1r,ias2r,ias2r_Sm
c
c nx4 , ny4 and nry  are given in common five as parameters
c
cSm030224
      nx4=nx+4
      ny4=ny+4
      ncx=nx4
      ncy=ny4
c      psi=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
c      dpsidr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z)
c      dpsidz=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z)
      dpsidr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z,nx4a)
      dpsidz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z,nx4a)
      return
      end

              
      subroutine calc_z_r_at_psi_at_field_line(r_mag,z_mag,psi,arteta,
     &r_psi_1d,z_psi_1d)
c----------------------------------------------------------------
c     calculate coordinates z_psi,r_psi at the flux surface psi
c     at poloidal angles arteta(1:nteta1) 
c     The point with artheta(1) consides with the point
c     at arteta(nteta1)
c     nteta and nteta1=nteta+1 were set in the param.i file
c-----------------------------------------------------------
      implicit none
      include 'param.i'
c-----input
      real*8
     & z_mag,r_mag,   !magnetic axis coordinates
     & psi,           !poloidal flux
     & arteta(nteta1) !poloidal angle mesh
c-----output
      real*8 r_psi_1d(nteta1),z_psi_1d(nteta1) !points at psi flux surface
                                         !at given poloidal angles 
                                         !arteta(nteta1)

c-----local
      real*8 pi,hteta,hpsi,
     &step_r,r_l,r_r,r_c,psi1,psi2,
     &eps_psi,theta_pol,step_pol,r_psi_l,z_psi_l,p

      integer i,j,i_count
c-----for RK solver
      integer ndim,n_psi_add,
     &n_psi_add_max  !maximal number of n_psi_add,

      real*8
     &rho_geom(1),d_rho_d_theta(1),aux_lim(8,1),
     &rho_geom_add(1),rho_geom_ini(1),theta_pol_ini
 
      real*8 
     &z_mag_l,r_mag_l  
      common /RHS_limiter_line_common/
     &z_mag_l,r_mag_l    
  
c-----additional RK steps
      integer n_add,k
      real*8 step_add,theta_pol_add

c-----external
      real*8 fpsi
      external RHS_limiter_line

      z_mag_l=z_mag
      r_mag_l=r_mag

cSAP091108     
      n_psi_add_max=100

      eps_psi=1.d-6    !set accuracy

      pi=4.d0*datan(1.d0)
      hteta=arteta(2)-arteta(1)


c     arpsi(j) j=1,npsi
c
c     arpsi(npsi)=psilim
c     arpsi(1)=psimag
c
c     hpsi=(psilim-psimag)/(npsi-1)
c-------------------------------------------------
c     arteta(i) i=1,nteta1, nteta1=nteta+1
c
c     arteta(i)=hteta*(i-0.5)
c     hteta = 2.d0*pi/nteta
c
c     arteta(1) = hteta*0.5 
c     arteta(nteta1) =[2.d0*pi/nteta](nteta1-0.5)=
c                    =[2.d0*pi/nteta](nteta+0.5)=
c                    =2*pi(1 + 0.5/nteta)=2pi+0.5hteta
c    It means that the point with arteta(nteta1) coinsides
c    with the point arteta(nteta1)
c
      
      z_psi_1d(1)=z_mag
      r_psi_1d(1)=r_mag

c----------------------------------------------------------- 
c     calculate the major radius of the initial point
c     at the field line 
c     r_0 at LCFS at z=yma (at theta_pol=0) 
c
c     To find r_0 we will solve the equation
c     F(r)=psi(r,z=z_mag)-psilim=0
c-----------------------------------------------------------
      step_r=r_mag/100.d0
      r_r=r_mag

c      write(*,*)'r_mag,z_mag',r_mag,z_mag

c      write(*,*)'gr2ne.f in calc_z_r_at_psi psi',psi
   
      i_count=1


 10   r_l=r_r
      r_r=r_r+step_r
      
      if (fpsi(r_r,z_mag).lt.psi) then 
         goto 10        
      else
c---------------------------------------------------
c        Boundaries r_l,r_r were found
c        fpsi(r_l,z_mag) < psi <fpsi(r_r,z_mag)
c        
c        Solution of the equation F(r)=psi(r,z=z_mag)-psi=0 at 
c        using the binary method
c--------------------------------------
         do while (dabs(r_r-r_l).gt.eps_psi*1.d-1)

             i_count=i_count+1
c             write(*,*)'psi,i_count r_r-r_l',
c     &                  psi,i_count,r_r-r_l

             r_c=r_l+(r_r-r_l)*0.5d0          
             psi1=fpsi(r_c,z_mag)-psi
             psi2=fpsi(r_r,z_mag)-psi
             if ((psi1*psi2).gt.0) then
               r_r=r_c
             else
               r_l=r_c
             end if
         end do
c        -----------------------------------------------
c        end of the binary method gives r_0=r_c
c        -----------------------------------------------
      endif

      write(*,*)'gr2new in calc_z_r_at_psi r_c',r_c
      write(*,*)'initial point at the fiel line r_c',r_c,
     &  'fpsi(r_c,z_mag)-psi',fpsi(r_c,z_mag)-psi

c-----calculate limiter points at the line field with given psi
c     initial point at the field line: (r_c,z_mag) with theta_pol=0

      ndim=1
      
      rho_geom(1)=(r_c-r_mag) !at the initial =point
 
c      write(*,*)'gr2new.f nteta',nteta
 
      rho_geom_ini(1)=rho_geom(1)

      do i=1,nteta
 
c        write(*,*)'loop nteta i,psi',i,psi

c--------loop over poloidal points   

c-----------------------------------------------------------
c        calculate the point  r_0,z_0 at the flux surface psi  
c        r=r_mag+rho_geom*cos(theta_pol)
c        z=z_mag+rho_geom*sin(theta_pol)
c
c        from the solution of the field line equation 
c        using RK solver
c-----------------------------------------------------------
         if (i.eq.1) then
            theta_pol=0.d0
            step_pol=arteta(1)
         else
            theta_pol=arteta(i-1)
            step_pol=hteta
         endif

c         write(*,*)'gr2new.f i,theta_pol,step_pol',i,theta_pol,step_pol

c 20      continue
        
c         write(*,*)'gr2new.f  calc_z_r_at_psi i,arteta(i),n_psi_add',
c     &   i,arteta(i),n_psi_add

c         write(*,*)'gr2new.f after 20 i,theta_pol,step_pol',
c     &                                i,theta_pol,step_pol

         n_psi_add=1 !initialization, no additional steps

         theta_pol_ini=theta_pol
         rho_geom_ini(1)=rho_geom(1)

 20      continue

c        if (n_psi_add.ge.1) then
c--------using additional poloidal points for better RK solver accuracy
         step_add=step_pol/n_psi_add
         theta_pol_add=theta_pol_ini
         rho_geom(1)=rho_geom_ini(1)

c            write(*,*)'theta_pol_add',theta_pol_add

         do  k=1,n_psi_add

c            write(*,*)'k,theta_pol_add,step_add',
c     &                 k,theta_pol_add,step_add

c------------additional RK steps for better accuracy at k>1
                
             call drkgs0_limiter(theta_pol_add,step_add,
     &                            rho_geom,d_rho_d_theta,ndim,
     &                             RHS_limiter_line,aux_lim)

             theta_pol_add=theta_pol_add+step_add                
c             write(*,*)'gr2new k,rho_geom_add(1)',k,rho_geom_add(1)
         enddo

c            write(*,*)'gr2new rho_geom_add(1),rho_geom(1)',
c     &                        rho_geom_add(1),rho_geom(1)
c         else
c            call drkgs0_limiter(theta_pol,step_pol,
c     &               rho_geom,d_rho_d_theta,ndim,
c     &               RHS_limiter_line,aux_lim) 
c         endif 

         call rz_lim_on_rho_geom_theta(rho_geom,arteta(i),
     &   r_mag,z_mag,z_psi_l,r_psi_l)

c        write(*,*)'after rz_lim_on_rho_geom_theta rho_geom(1)',
c    &              rho_geom(1)

c--------fit the accuracy at the calculated point
c        write(*,*)'gr2new r_psi_l,z_psi_l,psi,fpsi(r_psi_l,z_psi_l)',
c    &                     r_psi_l,z_psi_l,psi,fpsi(r_psi_l,z_psi_l)

         p=dabs(psi-fpsi(r_psi_l,z_psi_l))
c        write(*,*)'dabs(fpsi(r_psi_l,z_psi_l)-psi)',p

         if (dabs(fpsi(r_psi_l,z_psi_l)-psi).gt.eps_psi) then
            n_psi_add=n_psi_add*2
    
c            write(*,*)'dabs(fpsi(r_psi_l,z_psi_l)-psi)',
c     &      dabs(fpsi(r_psi_l,z_psi_l)-psi)
c            write(*,*)'gr2new.f psi,i,n_psi_add',psi,i,n_psi_add
            
            if (n_psi_add.lt.n_psi_add_max) then
               goto 20
            else
               write(*,*)'WARNINGin gr2new.f'
               write(*,*)'in subroutine calc_z_r_at_psi_at_field_line'
               write(*,*)'it was made more than n_psi_add_max=',
     &                    n_psi_add_max
               write(*,*)'additional RK steps but'
               write(*,*)'dabs(fpsi(r_psi_l,z_psi_l)-psi) > eps_psi'
               write(*,*)'dabs(fpsi(r_psi_l,z_psi_l)-psi),eps_psi',
     &                    dabs(fpsi(r_psi_l,z_psi_l)-psi),eps_psi
               if (dabs(fpsi(r_psi_l,z_psi_l)-psi).gt.1.d-3) then
                  write(*,*)'The field line calculations'
                  write(*,*)'gave bad accuracy'
                  write(*,*)'In the point with cordinates psi,arteta(i)'
     &            ,psi,arteta(i)
                  write(*,*)'coordinates at field line r_psi_l,z_psi_l',
     &            r_psi_l,z_psi_l
                  write(*,*)'dabs(fpsi(r_psi_l,z_psi_l)-psi)=',
     &                       dabs(fpsi(r_psi_l,z_psi_l)-psi)
                  write(*,*)'Please set smaller psifacr in the input'
                  write(*,*)'genray.in or genray.dat' 
                  stop 'in gr2new chage psifactr ' 
              endif
            endif
         endif

c         write(*,*)'gr2new r_psi_l,z_psi_l,psi',r_psi_l,z_psi_l,psi

         p=dabs(psi-fpsi(r_psi_l,z_psi_l))
c         write(*,*)'psi,i,n_psi_add,dabs(psi-fpsi(r_psi_l,z_psi_l))',
c     &              psi,i,n_psi_add,p
         n_psi_add=1
         z_psi_1d(i)=z_psi_l
         r_psi_1d(i)=r_psi_l

c         write(*,*)'end of loop nteta psi',psi

      enddo !nteta

      z_psi_1d(nteta+1)=z_psi_1d(1)
      r_psi_1d(nteta+1)=r_psi_1d(1)

      return
      end
     
