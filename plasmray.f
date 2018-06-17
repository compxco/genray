c        *********************plasmaray************************
c        *                        -                           *
c        * this subroutine  detemines the point where	      *
c        * the EC ray intersects the plasma                   *
c        * boundary(zu0,ru0,phiu0)                            *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c        rst, zst (m),phist(radian) -start point coordinates	   !
c        alfast,betast-tor. and poloid.angl. of the start	   !
c        refractive index in the vacuum(radian)                    !
c        output parameters					   !
c        zu0,ru0,phiu0                                  	   !
c        iraystop -parameter for stoppage the ray calculation 	   !
c                  if iraystop =1                                  !    
c------------------------------------------------------------------
c        this program uses the following functions and subroutines !
c        integer function iregion(z,r)				   !
c        double precision function proot(p2,p1,p0)                 !
c        subroutine vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)	   !
c        double precision function psif(z,r)  		           !
c------------------------------------------------------------------
cSAP110316
c      subroutine plasmray(zst,rst,phist,alfast,betast,
c     1                    zu0,ru0,phiu0,iraystop,
c     1                    raypatt,xst,yst,nxst,nyst,nzst)
      subroutine plasmray(zst,rst,phist,alfast,betast,
     1                    zu0,ru0,phiu0,iraystop)

      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
c      character*8 raypatt
      iraystop=0

c-------------------------------------------
      zmaxd=zmax
      zmind=zmin
      rmaxd=rmax
      rmind=rmin
cBH040404      zmax=zmaxd+0.01d0
cBH040404      rmax=rmaxd+0.01d0
cBH040404      zmin=zmind-0.01d0
cBH040404      rmin=rmind-0.01d0
c      zmax=zmaxd+0.01d0
c      rmax=rmaxd+0.01d0
c      zmin=zmind-0.01d0
c      rmin=rmind-0.01d0

cBH040404:  Ensuring initial box is outside of source posn zst,rst.      

      zmax=max(zmaxd,zst)+0.01d0
      zmin=min(zmind,zst)-0.01d0
      rmax=max(rmaxd,rst)+0.01d0
      rmin=min(rmind,rst)-0.01d0
      if (rmin.lt. 0.d0) rmin=0.d0
c-------------------------------------------------------------
c     conditions that start point is inside the plasma
      ireg=iregion(zst,rst)
      
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'in plasmaray rmax,rmaxd',rmax,rmaxd
      write(*,*)'in plasmaray rmin,rmind',rmin,rmind
      write(*,*)'in plasmaray zmax,zmaxd',zmax,zmaxd
      write(*,*)'in plasmaray zmin,zmind',zmin,zmind
      write(*,*)'in plasmray zst,rst,phist,alfast,betast',
     1zst,rst,phist,alfast,betast
      write(*,*)'plasmray ireg=',ireg
      endif ! outprint
      

      if (ireg.eq.9) then
         psip=psif(zst,rst)
         if (psip.lt.psilim) then
cyup        write(*,*)'start point is inside the plasma'
	   zu0=rst
	   ru0=rst
           rmax=rmaxd
           rmin=rmind
           zmax=zmaxd
           zmin=zmind
	   return
	 end if
       end if
c---------------------------------------------------------------
       cosbet=dcos(betast)
       sinbet=dsin(betast)
       cosalf=dcos(alfast)
       cnz=sinbet
       cnx=cosbet*dcos(alfast+phist)
       cny=cosbet*dsin(alfast+phist)
       cnr=cosbet*cosalf
c       write(*,*)'cnxst=',cnx,'cnyst=',cny,'cnzst=',cnz,'cnr=',cnr
c       write(*,*)'sinbet=',sinbet,'cosbet',cosbet,'cosalf',cosalf
c       write(*,*)'phist',phist,'alfast+phist',alfast+phist
c-----------------------------------------------------------
c      vacume ray trajectory
c      X=X_st+p*cnx
c      Y=Y_st+p*cny
c      Z=Z_st+p*cnz
c      R**2=R_st**2+p**2*cosbeta**2+2*p*R_st*cos(betast)*cos(alfast)
c      p>0 -parameter
c-----------------------------------------------------------
       z0=zst
       r0=rst
c----- determination of the region index for start point
10     continue
cyup       write(*,*)'10 z0,r0',z0,r0
       ireg=iregion(z0,r0)
cyup       write(*,*)'ireg=',ireg
c-----------------------------------------------------------
       goto (1,2,3,4,5,6,7,8,9) ireg

1      if ((cnz.ge.0.d0).or.(cnr.gt.0d0)) then
         iraystop=1
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'1 ray can not go inside the plasma'
         endif ! outprint
         rmax=rmaxd
         rmin=rmind
         zmax=zmaxd
         zmin=zmind

	 return
       end if
       if (dabs(sinbet).lt.1.d-10) then
         r0=rmax-1.d-10
         go to 10
       end if
       a2=cosbet**2
       if (a2.lt.1.d-10) then
         z0=zmax-1.d-10
         go to 10
       end if
       a1=r0*cosbet*cosalf
       a0=r0**2-rmax**2
       p_rmax=proot(a2,a1,a0)
       p_zmax=(zmax-z0)/sinbet
       p=dmin1(p_rmax,p_zmax)
       call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
       z0=zp-1.d-10
       r0=rp-1.d-10
       go to 10
c---------------------------------------------------------------------
3      if ((cnz.le.0.d0).or.(cnr.gt.0d0)) then
         iraystop=1
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'3 ray can not go inside the plasma'
         endif ! outprint
         rmax=rmaxd
         rmin=rmind
         zmax=zmaxd
         zmin=zmind
         return
       end if

       if (dabs(sinbet).lt.1.d-10) then
         r0=rmax-1.d-10
         go to 10
       end if
       a2=cosbet**2
       if (a2.lt.1.d-10) then
         z0=zmin+1.d-10
       go to 10
       end if
       a1=r0*cosbet*cosalf
       a0=r0**2-rmax**2
       p_rmax=proot(a2,a1,a0)
       p_zmin=(zmin-z0)/sinbet
       p=dmin1(p_rmax,p_zmin)
       call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
       z0=zp+1.d-10
       r0=rp-1.d-10
       go to 10
c---------------------------------------------------------------------
5      if ((cnz.le.0.d0).or.(cnr.lt.0d0)) then
         iraystop=1
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'5 ray can not go inside the plasma'
         endif ! outprint
         rmax=rmaxd
         rmin=rmind
         zmax=zmaxd
         zmin=zmind
         return
       end if
       if (dabs(sinbet).lt.1.d-10) then
         r0=rmin+1.d-10
         go to 10
       end if
       a2=cosbet**2
       if (a2.lt.1.d-10) then
         z0=zmin+1.d-10
       go to 10
       end if
       a1=r0*cosbet*cosalf
       a0=r0**2-rmin**2
       p_rmin=proot(a2,a1,a0)
       p_zmin=(zmin-z0)/sinbet
       p=dmin1(p_rmin,p_zmin)
       call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
       z0=zp+1.d-10
       r0=rp+1.d-10
       go to 10
c-------------------------------------------------------------------
7      if ((cnz.ge.0.d0).or.(cnr.lt.0d0)) then
         iraystop=1
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'7 ray can not go inside the plasma'
         endif ! outprint
         rmax=rmaxd
         rmin=rmind
         zmax=zmaxd
         zmin=zmind
         return
       end if
       if (dabs(sinbet).lt.1.d-10) then
         r0=rmin+1.d-10
         go to 10
       end if
       a2=cosbet**2
       if (a2.lt.1.d-10) then
         z0=zmax-1.d-10
       go to 10
       end if
       a1=r0*cosbet*cosalf
       a0=r0**2-rmin**2
       p_rmin=proot(a2,a1,a0)
       p_zmax=(zmax-z0)/sinbet
       p=dmin1(p_rmin,p_zmax)
c      write(*,*)'in plasmray7 p_rmin=',p_rmin,'p_zmax=',p_zmax
       call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
       z0=zp-1.d-10
       r0=rp+1.d-10
       go to 10
c-------------------------------------------------------------------
2      if (cnr.gt.0d0) then
         iraystop=1
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'2 ray can not go inside the plasma'
         endif ! outprint
         rmax=rmaxd
         rmin=rmind
         zmax=zmaxd
         zmin=zmind
         return
       end if
cyup       write(*,*)'in 2  sinbet',sinbet
       if (dabs(sinbet).lt.1.d-10) then
         r0=rmax-1.d-10
         go to 10
       end if
       a2=cosbet**2
       if (a2.lt.1.d-10) then
         if (cnz.gt.0) then
           z0=zmax+1.d-10
	 else
           z0=zmin-1.d-10
	 end if
         go to 10
       end if
       a1=r0*cosbet*cosalf
       a0=r0**2-rmax**2
       p_rmax=proot(a2,a1,a0)
       if (cnz.gt.0) then
         p_zmax=(zmax-z0)/sinbet
         p=dmin1(p_rmax,p_zmax)
       else
         p_zmin=(zmin-z0)/sinbet
         p=dmin1(p_rmax,p_zmin)
       end if
c       write(*,*)'plasmray 2 before vacray'
       call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
c       write(*,*)'plasmray 2 after vacray'
       if (cnz.gt.0) then
         z0=zp+1.d-10
       else
         z0=zp-1.d-10
       end if
       r0=rp-1.d-10
cyup       write(*,*)'in 2 z0,r0,zrmax,rmax',z0,r0,zrmax,rmax
       go to 10
c-------------------------------------------------------------------
4      if (cnz.lt.0d0) then
         iraystop=1
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'4 ray can not go inside the plasma'
         endif ! outprint
         rmax=rmaxd
         rmin=rmind
         zmax=zmaxd
         zmin=zmind
         return
       end if
       if (dabs(sinbet).lt.1.d-10) then
         if (cnr.gt.0) then
           r0=rmax+1.d-10
	 else
           r0=rmin-1.d-10
	 end if
         go to 10
       end if
       a2=cosbet**2
       if (a2.lt.1.d-10) then
         z0=zmin+1.d-10
         go to 10
       end if
       p_zmin=(zmin-z0)/sinbet
       a1=r0*cosbet*cosalf
       a0=r0**2-rmax**2
         if (cnr.gt.0) then
           a0=r0**2-rmax**2
           p_rmax=proot(a2,a1,a0)
           p=dmin1(p_rmax,p_zmin)
	 else
           a0=r0**2-rmin**2
           p_rmin=proot(a2,a1,a0)
           p=dmin1(p_rmin,p_zmin)
	 end if
       call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
       if (cnr.gt.0) then
         r0=rp+1.d-10
       else
         r0=rp-1.d-10
       end if
       z0=zp+1.d-10
       go to 10
c-------------------------------------------------------------------
6      if (cnr.lt.0d0) then
         iraystop=1
         write(*,*)'6 ray can not go inside the plasma'
         rmax=rmaxd
         rmin=rmind
         zmax=zmaxd
         zmin=zmind
         return
       end if
       if (dabs(sinbet).lt.1.d-10) then
         r0=rmin+1.d-10
         go to 10
       end if
       a2=cosbet**2
       if (a2.lt.1.d-10) then
         if (cnz.gt.0) then
           z0=zmax+1.d-10
	 else
           z0=zmin-1.d-10
	 end if
         go to 10
       end if
       a1=r0*cosbet*cosalf
       a0=r0**2-rmin**2
       p_rmin=proot(a2,a1,a0)
       if (cnz.gt.0) then
         p_zmax=(zmax-z0)/sinbet
         p=dmin1(p_rmin,p_zmax)
       else
         p_zmin=(zmin-z0)/sinbet
         p=dmin1(p_rmin,p_zmin)
       end if
       call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
       if (cnz.gt.0) then
         z0=zp+1.d-10
       else
         z0=zp-1.d-10
       end if
       r0=rp+1.d-10
       go to 10
c-------------------------------------------------------------------
8      if (cnz.gt.0d0) then
         iraystop=1
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'8 ray can not go inside the plasma'
         endif ! outprint
         rmax=rmaxd
         rmin=rmind
         zmax=zmaxd
         zmin=zmind
         return
       end if
       if (dabs(sinbet).lt.1.d-10) then
         if (cnr.gt.0) then
           r0=rmax+1.d-10
	 else
           r0=rmin-1.d-10
	 end if
         go to 10
       end if
       a2=cosbet**2
       if (a2.lt.1.d-10) then
         z0=zmax-1.d-10
         go to 10
       end if
       p_zmax=(zmax-z0)/sinbet
       a1=r0*cosbet*cosalf
       if (cnr.gt.0) then
         a0=r0**2-rmax**2
         p_rmax=proot(a2,a1,a0)
         p=dmin1(p_rmax,p_zmax)
       else
         a0=r0**2-rmin**2
         p_rmin=proot(a2,a1,a0)
         p=dmin1(p_rmin,p_zmax)
       end if
       call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
       if (cnr.gt.0) then
         r0=rp+1.d-10
       else
         r0=rp-1.d-10
       end if
       z0=zp-1.d-10
       go to 10
c-------------------------------------------------------------------
9      continue
c-------------------------------------------------------------------
c     ray point is inside the region(9).
c     This region's boundaries are tangent lines to the limitter
c     flux surface
c------------------------------------------------------------------
c     pstep(part of r0x) is a parameter along the ray
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'ray is inside 9'
         endif ! outprint
         pstep=0.0020d0
	 p=0.0d0
	 call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
cyup         write(*,*)'in 9  after vacray z0,r0,p,zp,rp',z0,r0,p,zp,rp
         psip=psif(zp,rp)
	 epslim=0.001d0
         epslim=1.d-7
         epslim=1.d-8
cyup	 write(*,*)'in plasmaray 9 psip=',psip,'psilim=',psilim,
cyup     &   'epslim',epslim
c	 if (dabs(psip-psilim).lt.epslim) then
	 if ((psip.lt.psilim).and.(psip.gt.(psilim-epslim))) then
c---------------------------------------------------------------
c          the point where the ray intersects the plasma boundary
           ru0=r0
           zu0=z0
cyup	   write(*,*)' 1initial point r0=',r0,'z0=',z0
           rmax=rmaxd
           rmin=rmind
           zmax=zmaxd
           zmin=zmind
	   return
	 end if
c----------------------------------------------------------------
c        solution of the equation for intersection of ray and limiter
c----------------------------------------------------------------
         isp=0
20       pold=p
         p1=pold
         f1=psilim-psip
	 p=p+pstep
c         write(*,*)'20 before vacray p,p1 f1=psilim-psip',p,p1,f1,
c     1   'r0,z0,sinbet,cosbet,cosalf,p',
c     1   r0,z0,sinbet,cosbet,cosalf,p

	 call vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)

c	 write(*,*)'in 9 zmin,zp,zmax',zmin,zp,zmax,
c     1		    'rmin,rp,rmax',rmin,rp,rmax

	 if(((rp.lt.rmin).or.(rp.gt.rmax)).or.
     1	    ((zp.gt.zmax).or.(zp.lt.zmin))) then
           if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
           write(*,*)'9 ray can not go inside the plasma'
           endif ! outprint
	   iraystop=1
           rmax=rmaxd
           rmin=rmind
           zmax=zmaxd
           zmin=zmind
	   return
	 end if
         psip=psif(zp,rp)
	 p2=p
         f2=psilim-psip
	 f3=f2
c         write(*,*)'psip=psif(zp,rp)',psip,'f2=psilim-psip',f2,
c     1             'p,pold',p,pold
	 if (psilim.lt.psip) then
	 go to 20
	 end if
c---------------------------------------------------------------
c        p -parameter for ray point inside plazma
c        pold - parameter for ray point outside plasma
c---------------------------------------------------------------
c30	 if (dabs(f3).lt.epslim) then
cyup         write(*,*)'f3,epslim',f3,epslim
 30   continue
c      write(*,*)'30'
      if ((dabs(f3).lt.epslim).and.(f3.gt.0.d0)) then
c---------------------------------------------------------------
c        the point where ray intersects the plasma boundary
         ru0=rp
	 zu0=zp
cyup	 write(*,*)'2 initial point ru0=',ru0,'zu0=',zu0
c-----------------------------------------------------------------
         if (cnz.ne.0.d0) then
            p0=(zp-zst)/cnz
         else
            if (cosalf.lt.0.d0) then
             p0=-rst*cosalf-dsqrt(rst*rst*cosalf*cosalf-(rst**2-rp**2))
	    else
	     p0=-rst*cosalf+dsqrt(rst*rst*cosalf*cosalf-(rst**2-rp**2))
	    end if
         end if
c******************************control
         xp=rst*dcos(phist)+p0*cnx
         yp=rst*dsin(phist)+p0*cny
         rpp=dsqrt(xp*xp+yp*yp)
c        rpp1=dsqrt(rst**2+p0**2*cnr**2+2.d0*p0*rst*cosbet*cosalf)
c        zpp=zst+p0*cnz
         if (yp.gt.0d0) then
            phiu0=dacos(xp/rpp)
         else
            phiu0=-dacos(xp/rpp)
         end if
	 call vacray(zst,rst,p0,sinbet,cosbet,cosalf,zpp2,rpp2)
c        write(*,*)'in plasmaray control rpp=',rpp,'ru0=',ru0,
c     1  'phiu0=',phiu0,'zp=',zp,'zpp=',zpp,'rpp1=',rpp1,
c     2  'zpp2',zpp2,'rpp2',rpp2,
c     3  'xp=',xp,'yp=',yp

c----------------------------------------------------------------------
         rmax=rmaxd
         rmin=rmind
         zmax=zmaxd
         zmin=zmind
	 return
      end if
c---------------------------------------------------------------
      isp=isp+1
cyup      write(*,*)'plasmray.f isp= ',isp 
c      if (isp.eq.10) then
      if (isp.eq.100) then
        write(*,*)'plasmray.f isp.eq.100 stopping'   
        stop
      endif
      p3=(p1+p2)*0.5d0
c      write(*,*)'before vacray p1,p3,p2',p1,p3,p2
      call vacray(z0,r0,p3,sinbet,cosbet,cosalf,zp,rp)
c      write(*,*)'after vacray z0,r0,zp,rp',z0,r0,zp,rp
      psip=psif(zp,rp)
      f3=psilim-psip
c      write(*,*)'in plsamaray p3=',p3,'f3=',f3,'f1',f1
      if (f1*f3.lt.0.d0) then
          p1=p1
          p2=p3
          goto 30
      else
          p1=p3
	  p2=p2
	  goto 30
      end if
c--------------------------------------------------------------
      zu0=zp
      ru0=rp
      rmax=rmaxd
      rmin=rmind
      zmax=zmaxd
      zmin=zmind
cyup      write(*,*)'end of plasmaray'
      return
      end
c******************************************************************
      integer function iregion(z,r)
c----------------------------------------------------------------------
c     this function finds the number of the region
c     (outside or inside the plasma) in which	the point(r,z)
c     is located
c----------------------------------------------------------------------
c   z	       regions
c   !
c   !         !		 !
c   !	7     !	   8	 !      1
c   !	      !	 	 !
c   !	---------zmaz------------------
c   !	      !  plasma	 !
c   !	6    rmin  9	rmax    2
c   !         !	 cord	 !
c   !   ---------zmin------------------
c   !	      !		 !
c   !	5     !    4	 !      3
c   !------------------------------------------------r
c   input parameters:zmax,zmin,rmax,rmin -are in common/five/
c--------------------------------------------------------------
        implicit double precision (a-h,o-z)
        include 'param.i'
      include'five.i'
c	write(*,*)'in ireg r=',r,'z=',z
c	write(*,*)'rmax=',rmax,'rmin=',rmin,'zmax=',zmax,'zmin=',zmin
	if((r.gt.rmax).and.(z.gt.zmax)) then
	 iregion=1
c	write(*,*)'in iregion=',iregion
	 return
	end if
	if((r.gt.rmax).and.((z.le.zmax).and.(z.ge.zmin) )) then
	 iregion=2
c	write(*,*)'in iregion=',iregion
	 return
	end if
	if((r.gt.rmax).and.(z.lt.zmin)) then
	 iregion=3
c	write(*,*)'in iregion=',iregion
	 return
	end if
	if(((r.ge.rmin).and.(r.le.rmax)).and.(z.lt.zmin)) then
	 iregion=4
c	write(*,*)'in iregion=',iregion
	 return
	end if
	if((r.lt.rmin).and.(z.lt.zmin)) then
	 iregion=5
c	write(*,*)'in iregion=',iregion
	 return
	end if
	if((r.lt.rmin).and.((z.le.zmax).and.(z.ge.zmin))) then
	 iregion=6
c	write(*,*)'in iregion=',iregion
	 return
	end if
	if((r.lt.rmin).and.(z.gt.zmax)) then
	 iregion=7
c	write(*,*)'in iregion=',iregion
	 return
	end if
	if(((r.ge.rmin).and.(r.le.rmax)).and.(z.gt.zmax)) then
	 iregion=8
c	write(*,*)'in iregion=',iregion
	 return
	end if
c	write(*,*)'in irigion beg 9 rmin,r,rmax',rmin,r,rmax,
c     1	'zmin,z,zmax',zmin,z,zmax
	if(((r.ge.rmin).and.(r.le.rmax)).and.
     1      ((z.gt.zmin).and.(z.lt.zmax))) then
        iregion=9
c	write(*,*)'in iregion=',iregion
	end if
	return
	end
c******************************************************************
	double precision function proot(p2,p1,p0)
c-------------------------------------------------------------------
c       smaller positive root of the equation p2*X**2+2*p1*X+p0=0
c       for vac.ray
c-------------------------------------------------------------------
        implicit double precision (a-h,o-z)
	det=p1*p1-p0*p2
	if (det.lt.0.d0) then
	  iraystop=1
	  write(*,*)'in plasmray  in proot det<0 '
	end if
	prootm=(-p1-dsqrt(det))/p2
	prootp=(-p1+dsqrt(det))/p2
c	write(*,*)'in proot p2,p1,p0',p2,p1,p0,
c     1  'det',det,'prootm',prootm,'prootp',prootp
	if(prootm*prootp.gt.0.d0) then
	  if (prootm.lt.0.d0) then
c           write(*,*)'in proot two roots .lt.0'
	    stop
	  else
            proot=dmin1(prootm,prootp)
	  end if
	else
	  if (prootm.lt.0.d0) then
	    proot=prootp
	  else
	    proot=prootm
	  end if
	end if
	return
	end
c**********************************************************************
	subroutine vacray(z0,r0,p,sinbet,cosbet,cosalf,zp,rp)
c----------------------------------------------------------------------
c       this subroutine determines the vacuum ray point(zp,rp)
c      	input:z0,r0,p,sinbet,cosbet,cosalf
c       output:zp,rp
c---------------------------------------------------------------------
        implicit double precision (a-h,o-z)
	zp=z0+p*sinbet
	rp=dsqrt(r0*r0+p*p*cosbet*cosbet
     1  	 +2.d0*p*r0*cosbet*cosalf)
c	write(*,*)'in vacray z0,r0,p',z0,r0,p,
c     1	'sinbet,cosbet,cosalf',sinbet,cosbet,cosalf,'zp=',zp,'rp=',rp
	return
	end
c*********************************************************************
      double precision function psif(z,r)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'five.i'
      double precision ias2r,ias2r_Sm
c
c nx4a , ny4a and nrya  are given in param.i as parameters
c
cSm030224
      nx4=nx+4
      ny4=ny+4

      ncx=nx4
      ncy=ny4
c      psif=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
      psif=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nx4a)

      return
      end

c        *********************edgcor   ************************
c        *                        -                           *
c        * this subroutine  shifts the point(where	      *
c        * the EC ray intersects the plasma)inside the plasma * 
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c        r1, z1 coordinates	before the shift                   !
c        output parameters					   !
c        zu0,ru0,phiu0                                  	   !
c        iraystop -parameter for stoppage the ray calculation 	   !
c                  if iraystop =1                                  !            !
c------------------------------------------------------------------
      subroutine edgcor(z,r,zu0,ru0)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      double precision
     1 ias1r,ias2r,ias2r_Sm
      epsbnd=1.d-7 !it must be equal  ebsbnd in boundc

      delcor=0.9995d0
      delcor=1.d0-epsbnd

cyup      write(*,*)' in edgcor initial :z,r',z,r
1     continue
cc      write(*,*)'in edgcor after 1 continue'
      iedg=0
      if (r.lt.rmin+epsbnd) then
cc        write(*,*)'in edgcor z,rmin+epsbnd,r,rmax-epsbnd',
cc     1  z,rmin+epsbnd,r,rmax-epsbnd
	iedg=1
        goto 10
      endif

      if (r.gt.rmax-epsbnd) then
cc        write(*,*)'in edgcor z,rmin+epsbnd,r,rmax-epsbnd',
cc     1  z,rmin+epsbnd,r,rmax-epsbnd
	iedg=1
      end if
 10   continue

      rrr=r
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      idx=0
      ipx=ip
      ipx4=ip+4
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
      zm=zzrm
      if ((z.gt.zp-epsbnd).or.(z.lt.zm+epsbnd)) then
cc        write(*,*)'in edgcor rrr,zm,z,zp',rrr,zm,z,zp
cc        write(*,*)'in edgcor zm+epsbnd,z,zp-epsbnd',
cc     1	zm+epsbnd,z,zp-epsbnd
	iedg=1
      end if
      if (iedg.eq.1) then
c       shift of the point (z,r) inside the plasma
        delz=z-yma
        delr=r-xma
        r=xma+delcor*delr
        z=yma+delcor*delz
        iedg=0
cc	write(*,*)'in engcor new z,r',z,r
	goto 1
      endif
      zu0=z
      ru0=r
cyup      write(*,*)'edgcor new z,r',z,r

cSAP091127
      delz=z-yma
      delr=r-xma
      r=xma+delr*shift_rho_geom
      z=yma+delz*shift_rho_geom 

      zu0=z
      ru0=r
cyup      write(*,*)'edgcor after shift_rho_geom z,r',
cyup     &shift_rho_geom,z,r

      return
      end





