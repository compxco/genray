c        ************************ b ***************************     
c        this function calculates the components  of  the   
c        magnetic field bz,br,bphi, 
c        the absolute value of the total magnetic field bmod
c        and the derivatives of bz,br,bphi,bmod by rz,r,phi 
c        It calculates the rho,rho2,a for the programms     
c        which use the density profile :dense,dxdr,dxdz     
c	 						      
c        It controls if the ray point is inside the plasma  
c        if the ray point is outside the plasma then the    
c        calculations will be finished 		      
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c------------------------------------------------------------------
      double precision function b(z,r,phi)
      implicit none !double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'six.i'
      real*8 z,r,phi !INPUT
      real*8 ias1r,ias2r_Sm ! external
cSm070324
      real*8 thetapol, l_zr, l_z_b_r_b !local
      real*8 epsbnd , x_l,y_l,theta_l, z_b,r_b !local
      integer idx,ipx,ipx4,nx4,ny4
      real*8 res,rrr,zm,zp,zzrm,zzrp,ffr,ffd,dffr,dres,zer,pp,pb,ppb
      real*8 dpsidr,dpsidz,rhopsi,spol
      real*8 dpdzr,dpdrr,dpdzz
      real*8 bph1,dbph1dph,dbph1dz,dbph1dr,
     +     bz1,dbz1dph,dbz1dz,dbz1dr,
     +     br1,dbr1dph,dbr1dz,dbr1dr

      
      integer ifirstc
      data ifirstc/1/
      save ifirstc
c--------------------------------------------------------------
c     control that the point is inside the plasma
      epsbnd=1.d-10
      iboundb=-1
c-----------------------
c      write(*,*)'b.f'
cSAP080817
c      write(*,*)'b.f z,r,xma,yma',z,r,xma,yma
      x_l=r-xma
      y_l=z-yma
c      write(*,*)'b.f x_l,y_l,theta_l', x_l,y_l,theta_l
      if ((x_l**2+y_l**2).lt.1.d-10) then
         goto 10
      endif

c-----calculate poloidal angle  thetapol
      call theta_xy(x_l,y_l,theta_l) !0 =< theta_l < 2*pi

c      write(*,*)'b.f x_l,y_l,theta_l', x_l,y_l,theta_l

c      write(*,*)'b.f theta_l',theta_l
c-----calculate coordinates x_b,z_b at the limiter surface psi=psilim
c     and the poloidal angle thetapol
 
c      write(*,*)'in b.f before zr_psith(psilim,theta_l,z_b,r_b)'
c      write(*,*)'in b.f psilim,theta_l',psilim,theta_l

      call zr_psith(psilim,theta_l,z_b,r_b)
 
      !write(*,*)'in b.f after zr_psith(psilim,theta_l,z_b,r_b)'
c      write(*,*)'in b.f z_b,r_b',z_b,r_b
      
c      write(*,*)'b.f r_b,z_b,fpsi(r_b,z_b),psilim',
c     &r_b,z_b,fpsi(r_b,z_b),psilim

      l_zr=dsqrt(x_l**2+y_l**2)
      l_z_b_r_b=dsqrt((r_b-xma)**2+(z_b-yma)**2)

      !write(*,*)'b.f l_zr,l_z_b_r_b',l_zr,l_z_b_r_b

      if (l_zr.ge.l_z_b_r_b) then
         rho = l_zr / l_z_b_r_b         
c        drhodzb=(z-yma)/(((r_b-xma)**2+(z_b-yma)**2)*rho)

         iboundb=3
      endif
c      write(*,*)'b.f rho',rho
      goto 10
c----------------------------------------------------
      
      if ((r.le.rmin+epsbnd).or.(r.ge.rmax-epsbnd)) then
c         write(*,*)'in b rmin+epsbnd,r,rmax-epsbnd',
c     1   rmin+epsbnd,r,rmax-epsbnd
         iboundb=2
c	 delr=0.5d0*(rmax-rmin)
c         if (r.le.rmin+epsbnd) rho=1.d0+(rmin-r+epsbnd)/delr
c         if (r.ge.rmax-epsbnd) rho=1.d0+(r-rmax+epsbnd)/delr

cSm070324 
c         write(*,*)'b.f 2 rho',rho
         theta_l=thetapol(z,r) 

c         write(*,*)'in b.f 2 before zr_psith(psilim,theta_l,z_b,r_b)'
c         write(*,*)'in b.f psilim,theta_l',psilim,theta_l

         call zr_psith(psilim,theta_l,z_b,r_b)

c         write(*,*)'in b.f 2 before zr_psith(psilim,theta_l,z_b,r_b)'
c         write(*,*)'in b.f z_b,r_b',z_b,r_b

         rho=dsqrt(((r-xma)**2+(z-yma)**2)/((r_b-xma)**2+(z_b-yma)**2))

c         write(*,*)'b.f rho',rho

         goto 10
      end if
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      rrr=r
      idx=0
      ipx=ip
      ipx4=ip+4
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
      zm=zzrm
cSAP080816
      if ((z.ge.zp-epsbnd).or.(z.le.zm+epsbnd)) then 
c      if ((z.gt.zp).or.(z.lt.zm)) then
c         write(*,*)'in b rrr,zm,z,zp',rrr,zm,z,zp
c         write(*,*)'in b zm+epsbnd,z,zp-epsbnd',
c     1	 zm+epsbnd,z,zp-epsbnd
c	 delz=0.5d0*(zmax-zmin)
c         if (z.ge.zp-epsbnd) then
c	    rho=1.d0+(z-zp+epsbnd)/delz
c	    drhodzb=1/delz
c	 else
c            rho=1.d0+(zm-z+epsbnd)/delz
c	    drhodzb=-1/delz
c	 endif
         iboundb=1
cSm070324 
        
         theta_l=thetapol(z,r) 
         call zr_psith(psilim,theta_l,z_b,r_b)
         rho=dsqrt(((r-xma)**2+(z-yma)**2)/((r_b-xma)**2+(z_b-yma)**2))
         drhodzb=(z-yma)/(((r_b-xma)**2+(z_b-yma)**2)*rho)

c         write(*,*)'in b.f iboundb=1 z,r,z_b,r_b,rho',
c     &                               z,r,z_b,r_b,rho

         goto 10
      end if
c---------------------------------------------------------------------
 10     continue
c end of the control that the point is inside the plasma
c--------------------------------------------------------------
c
c nx4a , ny4a and nrya  are given in common five as parameters
c nxeqda,nyeqda,nlimit are given in common four as parameters
c

cSm_030224
c      ncx=nx4
c      ncy=ny4
c      res=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
c      dpsidr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z)
c      dpsidz=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z)
      nx4=nx+4
      ny4=ny+4
      ncx=nx4
      ncy=ny4
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nx4a)
      dpsidr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z,nx4a)
      dpsidz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z,nx4a)
   
c--------------------------------------------------------------------
c     idx derivativs order 0.ge.idx.le.3
c--------------------------------------------------------------------
      nx4=nx+4  
      if (iboundb.ge.1) then
c        the point is outside the plasma
         idx=0
         ffr=ias1r(txf,nx,nx4,cx,idx,psilim)
         ffd=ffr
         dffr=0.d0
         dres=dffr
      else
c        the point is inside the plasma
         idx=0
         ffr=ias1r(txf,nx,nx4,cx,idx,res)
         ffd=ffr
         idx=1
         dffr=ias1r(txf,nx,nx4,cx,idx,res)
         dres=dffr
      endif
      psid=res
      dpdzd=dpsidz
      dpdrd=dpsidr
c------------------------------------------------------------------
c this  part of the programm calculates the rho,rho2 and a
c for the functions dense,dxdr,dxdz

       a=1.d0       
       if (iboundb.eq.-1) rho=rhopsi(res)

c       write(*,*)'in b.f: r,z,iboundb,psi,rho', r,z,iboundb,res,rho

       zer=0.d0
       rho2=rho*rho
c------------------------------------------------------------------
c     \vector{B}=B_phi+\vector{B_pol}
c     B_phi=(Feqd(psi)/R)\hat{phi}
c     \vector{B_pol}=(1/R)*grad(psi) x \hat{phi}
c     here x is the vector product
c     then  B_z=(1/R)*dpsi/dr and B_r=-(1/R)*dpsi/dz  (*).
c     If the total current is not given, then the sign of the B_pol
c     determined fro psi, using the previous formula (*).
c     If the current 'toteqd'is given, then the sign of it will indicate 
c     the direction of B_pol, regardless of the convention for psi(R,Z)
c     toteqd is positive if it is directed along \hat{phi}.
c-------------------------------------------------------------------
c     The poloidal flux	array peqd(i,j) was transformed in equilib.f
c     in subroutine input. It was multiplied by dpsimax.
c     If (psimag.gt.psilim) then dpsimax=-1.d0 else  dpsimax=1.d0
c     So after this multiplication poloidal flux peqd(i,j)
c     always has the minimum on tme magnetic axis.
c-------------------------------------------------------------------
      pp=1.d0/r
      bz=-dpdrd*pp
      br=dpdzd*pp
      bphi=ffd*pp
 
	if(toteqd.ge.0.d0) spol=-1.d0
	if(toteqd.lt.0.d0) spol=1.d0
	spol=-spol
	if(toteqd.eq.0.d0) then
           if (ifirstc.eq.1) then
              write(*,*)'current is not given in eqdsk (toteqd=0)' 
              ifirstc=0
           endif
c          The total current is not given. The sign of the B_pol will be
c          determined from psi, using the original poloidal flux from eqdsk.
c          In this case the poloidal flux had max on the magnetic axis
c          (and dpsimax=-1d0), or it had min (and dpsimax=1.d0).
           spol=dpsimax
	endif
      
	bz=bz*spol
	br=br*spol
cSm030224
c              dpdzr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,1,1,r,z)
c              dpdrr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,2,0,r,z)
c              dpdzz=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,2,r,z)
              dpdzr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,1,r,z,nx4a)
              dpdrr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,2,0,r,z,nx4a)
              dpdzz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,2,r,z,nx4a)
                            !(tx,nx,ty,ny,cxy,ncx,ncy,idx,idy,xx,yy,nx4a)
                            !In func: tx(*),   ty(*),   cxy(nx4a,*)
                         ! in five.i: tx(nx4a),ty(ny4a),cxy(nx4a,ny4a)
                         !               nx4a=nxeqda+4, ny4a=nyeqda+4
		dpdzrd=dpdzr
		dpdrrd=dpdrr
		dpdzzd=dpdzz
c------------------------------------------------------------------
      dbzdz=-pp*dpdzrd
      dbzdr=pp*pp*dpdrd-pp*dpdrrd
      dbrdz=pp*dpdzzd
      dbrdr=-pp*pp*dpdzd+pp*dpdzrd
	dbzdz=dbzdz*spol
	dbzdr=dbzdr*spol
	dbrdz=dbrdz*spol
	dbrdr=dbrdr*spol
      dbzdph=0.d0
      dbrdph=0.d0
      dbpdph=0.d0
      dbphdz=pp*dres*dpdzd
      dbphdr=-pp*bphi+pp*dres*dpdrd
c--------------------------------------------------------------------
c      deltripl=0.00d0
      if(deltripl.gt.0.d0)then ! YuP[2020-08-19] Added: skip these two subr.
                               ! when not needed (saves cpu time)
        if (i_ripple.eq.1)then
          !-------- the ripple model approximating the DIII-D field
           call briplmd4(z,r,phi,
     +     bph1,dbph1dph,dbph1dz,dbph1dr,
     +     bz1,dbz1dph,dbz1dz,dbz1dr,
     +     br1,dbr1dph,dbr1dz,dbr1dr)
        endif

        if (i_ripple.eq.2)then
           !--------the ripple model using modified Bessel function I_0
           !        see the GENRAY description
           call briplmd5(z,r,phi,
     +     bph1,dbph1dph,dbph1dz,dbph1dr,
     +     bz1,dbz1dph,dbz1dz,dbz1dr,
     +     br1,dbr1dph,dbr1dz,dbr1dr)
        endif
 
        bz=bz+bz1
        br=br+br1
        bphi=bphi+bph1

        dbzdz=dbzdz+dbz1dz
        dbzdr=dbzdr+dbz1dr
        dbzdph=dbzdph+dbz1dph

        dbrdz=dbrdz+dbr1dz
        dbrdr=dbrdr+dbr1dr
        dbrdph=dbrdph+dbr1dph

        dbpdph=dbpdph+dbph1dph
        dbphdz=dbphdz+dbph1dz
        dbphdr=dbphdr+dbph1dr
      endif !(deltripl.gt.0.d0)then ! YuP[2020-08-19]
      
c-----------------------------------------------------------------
c     calculation of the model stellarator W7 magnetic field
c     bmagloc=2.3d0 !Tl
c      call bstelw7(z,r,bmagloc,
c     1bphi,dbpdph,dbphdz,dbphdr,
c     2bz,dbzdph,dbzdz,dbzdr,
c     3br,dbrdph,dbrdz,dbrdr)

c-----------------------------------------------------------------

c------------------------------------------------------------------
      pb=dsqrt(bphi*bphi+bz*bz+br*br)
c------------------------------------------------------------------
	ppb=1.d0/pb
	dbmdz=ppb*(bz*dbzdz+br*dbrdz+bphi*dbphdz)
	dbmdr=ppb*(bz*dbzdr+br*dbrdr+bphi*dbphdr)
	dbmdph=ppb*(bz*dbzdph+br*dbrdph+bphi*dbpdph)

c---------------------------------------------------------------------
        b=pb
c        write(*,*)'in b: z,r,phi,bz,br,bphi,rho,b',
c     &                   z,r,phi,bz,br,bphi,rho,b
      
      return
      end function b


c*************************briplmod*********************************
c  F=sin(N_loop*phi)*g(r,z)   - is a vacume ripple field potential *
c  bripl_phi(z,r,phi)=(dF/dphi)/r=cos(N_loop*phi)*g(r,z)*N_loop/r  *
c  bripl_z(z,r,phi) =(dF/dz)    =sin(N_loop*phi)*(dg/dz) 	   *
c  bripl_r(z,r,phi) =(dF/dr)    =sin(N_loop*phi)*(dg/dr) 	   *
c  model for function g:					   *
c  g_model=beqd*deltripl*exp(N*((rho(r,z)-1)))			   *
c  N is a number of the toroidal field coils			   *
c------------------------------------------------------------------
c        input parameters					   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c  nloop       a number of the toroidal field coils  (in one.i)    !
c  deltripl	   is an amplitude of the ripple field at the	   !
c                  last flux surface (at rho=1)(in one.i)	   !
c      bmag      - is a magnetic field on the magnetic axis (Tl)   !
c------------------------------------------------------------------
c        output parameters:					   |
c        components of the ripple magnetic field, and their	   |
c        derivatives                                               |
c        bph1,dbph1dph,dbph1dz,dbph1dr	-toroidal  (phi) components|
c        bz1,dbz1dph,dbz1dz,dbz1dr	-Z               components|
c        br1,dbr1dph,dbr1dz,dbr1dr	-R               components|
c------------------------------------------------------------------
      subroutine briplmod(z,r,phi,
     1bph1,dbph1dph,dbph1dz,dbph1dr,
     2bz1,dbz1dph,dbz1dz,dbz1dr,
     3br1,dbr1dph,dbr1dz,dbr1dr)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
cc      write(*,*)'in briplmod'
      rnloop=dfloat(nloop)
      cosphi=dcos(phi*rnloop)
      sinphi=dsin(phi*rnloop)
cc      write(*,*)'cosphi,sinphi',cosphi,sinphi
c------------------------------------------------------------------
      g=beqd*deltripl*dexp((rho-1.d0)*rnloop)
      dgdrho=g*rnloop
      dgdrhrh=dgdrho*rnloop
cc      write(*,*)'beqd,deltripl,rho,rnloop',beqd,deltripl,rho,rnloop
cc      write(*,*)'g,dgdrho,dgdrhrh',g,dgdrho,dgdrhrh

c ??  drhodz, drhodr
      psi=psid
      dro_dpsi=drhopsi(psi)
      drhodz=dro_dpsi*dpdzd
      drhodr=dro_dpsi*dpdrd
cc      write(*,*)'psi,dpdzd,dpdrd',psi,dpdzd,dpdrd
cc      write(*,*)'dro_dpsi,drhodz,drhodr',dro_dpsi,drhodz,drhodr

c ??  drhodzdz, drhodrdr,drhodzdr
      dropsps=drhpsps(psi)
      drhodzdz=dropsps*dpdzd*dpdzd+dro_dpsi*dpdzzd
      drhodrdr=dropsps*dpdrd*dpdrd+dro_dpsi*dpdrrd
      drhodzdr=dropsps*dpdzd*dpdrd+dro_dpsi*dpdzrd

      dgdz=dgdrho*drhodz
      dgdr=dgdrho*drhodr

      dgdzdz=dgdrhrh*drhodz*drhodz+dgdrho*drhodzdz
      dgdrdr=dgdrhrh*drhodr*drhodr+dgdrho*drhodrdr
      dgdzdr=dgdrhrh*drhodz*drhodr+dgdrho*drhodzdr
c------------------------------------------------------------------
      bz1=sinphi*dgdz
      br1=sinphi*dgdr
      bph1=cosphi*g*rnloop/r
cc      write(*,*)'bz1,br1,bph1',bz1,br1,bph1
c------------------------------------------------------------------
      dbz1dz=sinphi*dgdzdz
      dbz1dr=sinphi*dgdzdr
      dbz1dph=cosphi*rnloop*dgdz
cc      write(*,*)'dbz1dz,dbz1dr,dbz1dph',dbz1dz,dbz1dr,dbz1dph

      dbr1dz=sinphi*dgdzdr
      dbr1dr=sinphi*dgdrdr
      dbr1dph=cosphi*rnloop*dgdr
cc     write(*,*)'dbr1dz,dbr1dr,dbr1dph',dbr1dz,dbr1dr,dbr1dph

      dbph1dz=cosphi*rnloop*dgdz/r
      dbph1dr=cosphi*rnloop*(dgdr/r-g/(r*r))
      dbph1dph=-sinphi*rnloop*g/r
cc      write(*,*)'dbph1dz,dbph1dr,dbph1dph',dbph1dz,dbph1dr,dbph1dph

      return
      end
      subroutine briplmd1(z,r,phi,
     1bph1,dbph1dph,dbph1dz,dbph1dr,
     2bz1,dbz1dph,dbz1dz,dbz1dr,
     3br1,dbr1dph,dbr1dz,dbr1dr)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
cc      write(*,*)'in briplmod1'
      rnloop=dfloat(nloop)
      cosphi=dcos(phi*rnloop)
      sinphi=dsin(phi*rnloop)
cc      write(*,*)'cosphi,sinphi',cosphi,sinphi
c------------------------------------------------------------------
      sinr=dsin(r-xma)
      cosr=dcos(r-xma)
      sinz=dsin(z)
      cosz=dcos(z)

      g=beqd*deltripl*sinr*sinz
c------------------------------------------------------------------
      bz1=g*sinphi*rnloop
      br1=g*sinphi*0.9d0*rnloop
      bph1=cosphi*g
cc      write(*,*)'bz1,br1,bph1',bz1,br1,bph1
c------------------------------------------------------------------
      dbz1dz=beqd*deltripl*sinr*cosz*sinphi*rnloop
      dbz1dr=beqd*deltripl*cosr*sinz*sinphi*rnloop
      dbz1dph=g*cosphi*rnloop*rnloop
cc      write(*,*)'dbz1dz,dbz1dr,dbz1dph',dbz1dz,dbz1dr,dbz1dph

      dbr1dz=0.9*beqd*deltripl*sinr*cosz*sinphi*rnloop
      dbr1dr=0.9*beqd*deltripl*cosr*sinz*sinphi*rnloop
      dbr1dph=0.9*g*cosphi*rnloop*rnloop
cc     write(*,*)'dbr1dz,dbr1dr,dbr1dph',dbr1dz,dbr1dr,dbr1dph

      dbph1dz=beqd*deltripl*sinr*cosz*cosphi
      dbph1dr=beqd*deltripl*cosr*sinz*cosphi
      dbph1dph=-sinphi*g*rnloop
cc      write(*,*)'dbph1dz,dbph1dr,dbph1dph',dbph1dz,dbph1dr,dbph1dph

      return
      end




c*************************briplmd2********************************
c  F=sin(N_loop*phi)*g(r,z)   - is a vacume ripple field potential *
c  bripl_phi(z,r,phi)=(dF/dphi)/r=cos(N_loop*phi)*g(r,z)*N_loop/r  *
c  bripl_z(z,r,phi) =(dF/dz)    =sin(N_loop*phi)*(dg/dz) 	   *
c  bripl_r(z,r,phi) =(dF/dr)    =sin(N_loop*phi)*(dg/dr) 	   *
c  model for function g:					   *
c  g_model=beqd*reqd*p1*(1+(N*rho(r,z)/2)**2)
c          p1=deltripl/(1+0.25*N**2)/N
c  N is a number of the toroidal field coils			   *
c------------------------------------------------------------------
c        input parameters					   !
c      z, r, phi - point where the components of magnetic field    !
c                  are calculated in.				   !
c  nloop           is number of  toroidal field coils(in one.i)    !
c  deltripl	   is an amplitude of the ripple field at the	   !
c                  last flux surface (at rho=1)(in one.i)			   !
c      bmag      - is a magnetic field on the magnetic axis (Tl)   !
c------------------------------------------------------------------
c        output parameters:					   |
c        components of the ripple magnetic field, and their	   |
c        derivatives                                               |
c        bph1,dbph1dph,dbph1dz,dbph1dr	-toroidal  (phi) components|
c        bz1,dbz1dph,dbz1dz,dbz1dr	-Z               components|
c        br1,dbr1dph,dbr1dz,dbr1dr	-R               components|
c------------------------------------------------------------------
      subroutine briplmd2(z,r,phi,
     1bph1,dbph1dph,dbph1dz,dbph1dr,
     2bz1,dbz1dph,dbz1dz,dbz1dr,
     3br1,dbr1dph,dbr1dz,dbr1dr)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
cc      write(*,*)'in briplmod2'
      rnloop=dfloat(nloop)
      cosphi=dcos(phi*rnloop)
      sinphi=dsin(phi*rnloop)
cc      write(*,*)'cosphi,sinphi',cosphi,sinphi
c------------------------------------------------------------------
c  g_model=beqd*p1*(1+(N*(rho(r,z)/2)**2)
      p1=deltripl/(1.d0+0.25d0*rnloop*rnloop)/rnloop
      p=rnloop*rho*0.5d0
      p2=p*p
      g=beqd*reqd*p1*(1.d0+p2)
      dgdrho=beqd*reqd*p1*0.5d0*rnloop*rnloop*rho
      dgdrhrh=beqd*reqd*p1*0.5d0*rnloop*rnloop
cc      write(*,*)'beqd,deltripl,rho,rnloop',beqd,deltripl,rho,rnloop
cc      write(*,*)'g,dgdrho,dgdrhrh',g,dgdrho,dgdrhrh

      psi=psid
      dro_dpsi=drhopsi(psi)
      drhodz=dro_dpsi*dpdzd
      drhodr=dro_dpsi*dpdrd
cc      write(*,*)'psi,dpdzd,dpdrd',psi,dpdzd,dpdrd
cc      write(*,*)'dro_dpsi,drhodz,drhodr',dro_dpsi,drhodz,drhodr

      dropsps=drhpsps(psi)
      drhodzdz=dropsps*dpdzd*dpdzd+dro_dpsi*dpdzzd
      drhodrdr=dropsps*dpdrd*dpdrd+dro_dpsi*dpdrrd
      drhodzdr=dropsps*dpdzd*dpdrd+dro_dpsi*dpdzrd

      dgdz=dgdrho*drhodz
      dgdr=dgdrho*drhodr

      dgdzdz=dgdrhrh*drhodz*drhodz+dgdrho*drhodzdz
      dgdrdr=dgdrhrh*drhodr*drhodr+dgdrho*drhodrdr
      dgdzdr=dgdrhrh*drhodz*drhodr+dgdrho*drhodzdr
c------------------------------------------------------------------
      bz1=sinphi*dgdz
      br1=sinphi*dgdr
      bph1=cosphi*g*rnloop/r
cc      bz1=0.d0
cc      br1=0.d0
cc      write(*,*)'bz1,br1,bph1',bz1,br1,bph1
c test
cc      ppp=beqd*reqd/r
cc      write(*,*)'test r,rho,phi,bphi,bph1/bphi',r,rho,phi,ppp,bph1/ppp
c------------------------------------------------------------------
      dbz1dz=sinphi*dgdzdz
      dbz1dr=sinphi*dgdzdr
      dbz1dph=cosphi*rnloop*dgdz
cc      dbz1dz=0.d0
cc      dbz1dr=0.d0
cc      dbz1dph=0.d0
cc      write(*,*)'dbz1dz,dbz1dr,dbz1dph',dbz1dz,dbz1dr,dbz1dph

      dbr1dz=sinphi*dgdzdr
      dbr1dr=sinphi*dgdrdr
      dbr1dph=cosphi*rnloop*dgdr
cc      dbr1dz=0.d0
cc      dbr1dr=0.d0
cc      dbr1dph=0.d0
cc      write(*,*)'dbr1dz,dbr1dr,dbr1dph',dbr1dz,dbr1dr,dbr1dph

      dbph1dz=cosphi*rnloop*dgdz/r
      dbph1dr=cosphi*rnloop*(dgdr/r-g/(r*r))
      dbph1dph=-sinphi*rnloop*rnloop*g/r
cc      write(*,*)'dbph1dz,dbph1dr,dbph1dph',dbph1dz,dbph1dr,dbph1dph
      return
      end
c*************************briplmd3********************************
c  only for irho=4,rho=sqrt((psi-psimag)/(psilim-psimag))
c  F=sin(N_loop*phi)*g(r,z)   - is a vacume ripple field potential *
c  bripl_phi(z,r,phi)=(dF/dphi)/r=cos(N_loop*phi)*g(r,z)*N_loop/r  *
c  bripl_z(z,r,phi) =(dF/dz)    =sin(N_loop*phi)*(dg/dz) 	   *
c  bripl_r(z,r,phi) =(dF/dr)    =sin(N_loop*phi)*(dg/dr) 	   *
c  model for function g:					   *
c  g_model=beqd*p1*(1+(N*(rho(r,z)/2))**2)  			   *
c          p1=deltripl/(1+0.25*N**2)				   *
c  N is a number of the toroidal field coils			   *
c  bripl_phi=cos(N_loop*phi)*N_loop/r*{beqd*p1*			   *
c             {1+0.25*N*N(psi-psimag)/(psilim-psimag)}}             *
c  bripl_z=sin(N_loop*phi)*{beqd*p1*0.25*N*N/(psilim-psimag)*dpdzd}*
c  bripl_r=sin(N_loop*phi)*{beqd*p1*0.25*N*N/(psilim-psimag)*dpdrd}*
c------------------------------------------------------------------
c        input parameters					   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c  nloop           is number of toroidal field coils (in one.i)    !
c  deltripl	   is an amplitude of the ripple field at the	   !
c                  last flux surface (at rho=1) (in one.i)	   !
c      bmag      - is a magnetic field on the magnetic axis (Tl)   !
c------------------------------------------------------------------
c        output parameters:					   |
c        components of the ripple magnetic field, and their	   |
c        derivatives                                               |
c        bph1,dbph1dph,dbph1dz,dbph1dr	-toroidal  (phi) components|
c        bz1,dbz1dph,dbz1dz,dbz1dr	-Z               components|
c        br1,dbr1dph,dbr1dz,dbr1dr	-R               components|
c------------------------------------------------------------------
      subroutine briplmd3(z,r,phi,
     1bph1,dbph1dph,dbph1dz,dbph1dr,
     2bz1,dbz1dph,dbz1dz,dbz1dr,
     3br1,dbr1dph,dbr1dz,dbr1dr)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
cc      write(*,*)'in briplmod'
      rnloop=dfloat(nloop)
      cosphi=dcos(phi*rnloop)
      sinphi=dsin(phi*rnloop)
cc      write(*,*)'cosphi,sinphi',cosphi,sinphi
c------------------------------------------------------------------
      p1=deltripl/(1.d0+0.25d0*rnloop*rnloop)
c------------------------------------------------------------------
      bz1=sinphi*beqd*p1*0.25d0*rnloop*rnloop/(psilim-psimag)*dpdzd
      br1=sinphi*beqd*p1*0.25d0*rnloop*rnloop/(psilim-psimag)*dpdrd
      bph1=cosphi*rnloop/r*beqd*p1*
     1     (1.d0+0.25d0*rnloop*rnloop*(psid-psimag)/(psilim-psimag))
cc      write(*,*)'bz1,br1,bph1',bz1,br1,bph1
c------------------------------------------------------------------
      dbz1dz=sinphi*beqd*p1*0.25d0*rnloop*rnloop/(psilim-psimag)*
     1       dpdzzd
      dbz1dr=sinphi*beqd*p1*0.25d0*rnloop*rnloop/(psilim-psimag)*
     1       dpdzrd
      dbz1dph=cosphi*rnloop*
     1beqd*p1*0.25d0*rnloop*rnloop/(psilim-psimag)*dpdzd
cc      write(*,*)'dbz1dz,dbz1dr,dbz1dph',dbz1dz,dbz1dr,dbz1dph

      dbr1dz=sinphi*beqd*p1*0.25d0*rnloop*rnloop/(psilim-psimag)*
     1       dpdzrd
      dbr1dr=sinphi*beqd*p1*0.25d0*rnloop*rnloop/(psilim-psimag)*
     1       dpdrrd
      dbr1dph=cosphi*rnloop*
     1beqd*p1*0.25d0*rnloop*rnloop/(psilim-psimag)*dpdrd
cc     write(*,*)'dbr1dz,dbr1dr,dbr1dph',dbr1dz,dbr1dr,dbr1dph

      dbph1dz=cosphi*rnloop/r*beqd*p1*
     1     0.25d0*rnloop*rnloop/(psilim-psimag)*dpdzd
      dbph1dr=cosphi*rnloop*beqd*p1*(
     1 0.25d0*rnloop*rnloop/(psilim-psimag)*dpdrd/r-
     1(1.d0+0.25d0*rnloop*rnloop*(psid-psimag)/(psilim-psimag))/(r*r))
      dbph1dph=-sinphi*rnloop*rnloop/r*beqd*p1*
     1     (1.d0+0.25d0*rnloop*rnloop*(psid-psimag)/(psilim-psimag))
cc      write(*,*)'dbph1dz,dbph1dr,dbph1dph',dbph1dz,dbph1dr,dbph1dph

      return
      end


c*************************briplmd4********************************
c  F=sin(N_loop*phi)*g(r,z)   - is a vacume ripple field potential *
c  bripl_phi(z,r,phi)=(dF/dphi)/r=cos(N_loop*phi)*g(r,z)*N_loop/r  *
c  bripl_z(z,r,phi) =(dF/dz)    =sin(N_loop*phi)*(dg/dz) 	   *
c  bripl_r(z,r,phi) =(dF/dr)    =sin(N_loop*phi)*(dg/dr) 	   *
c  model for function g:					   *
c  g_model=beqd*reqd*p1*(r/rmax)**N    			           *
c          p1=deltripl/N					   *
c  N is a number of the toroidal field coils			   *
c  rmax - the major radius of the outboard side of plasma	   *
c------------------------------------------------------------------
c        input parameters					   !
c      z, r, phi - point where the components of magnetic field    !
c                  are calculated in.				   !
c  nloop           is number of toroidal field coils(in one.i)     !
c  deltripl	   is an amplitude of the ripple field at the	   !
c                  last flux surface (at rho=1)	(in one.i)		   !
c      bmag      - is a magnetic field on the magnetic axis (Tl)   !
c------------------------------------------------------------------
c        output parameters:					   |
c        components of the ripple magnetic field, and their	   |
c        derivatives                                               |
c        bph1,dbph1dph,dbph1dz,dbph1dr	-toroidal  (phi) components|
c        bz1,dbz1dph,dbz1dz,dbz1dr	-Z               components|
c        br1,dbr1dph,dbr1dz,dbr1dr	-R               components|
c------------------------------------------------------------------
      subroutine briplmd4(z,r,phi,
     1bph1,dbph1dph,dbph1dz,dbph1dr,
     2bz1,dbz1dph,dbz1dz,dbz1dr,
     3br1,dbr1dph,dbr1dz,dbr1dr)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

c      write(*,*)'in briplmod4 deltripl,nloop',deltripl,nloop

      rnloop=dfloat(nloop)
      cosphi=dcos(phi*rnloop)
      sinphi=dsin(phi*rnloop)
     
c------------------------------------------------------------------
c  g_model=beqd*reqd*p1*(r/rmax)**N
c          p1=deltripl/N
      p1=deltripl/rnloop
      g=beqd*reqd*p1*(r/rmax)**nloop
      dgdr=g*(rnloop)/r
      dgdz=0.d0

      dgdrdr=dgdr*(rnloop-1.d0)/r
      dgdzdz=0.0d0
      dgdzdr=0.0d0
c------------------------------------------------------------------
      bz1=sinphi*dgdz
      br1=sinphi*dgdr
      bph1=cosphi*g*rnloop/r

c------------------------------------------------------------------
      dbz1dz=sinphi*dgdzdz
      dbz1dr=sinphi*dgdzdr
      dbz1dph=cosphi*rnloop*dgdz
cc      dbz1dz=0.d0
cc      dbz1dr=0.d0
cc      dbz1dph=0.d0
c      write(*,*)'briplmod4,dbz1dz,dbz1dr,dbz1dph',dbz1dz,dbz1dr,dbz1dph

      dbr1dz=sinphi*dgdzdr
      dbr1dr=sinphi*dgdrdr
      dbr1dph=cosphi*rnloop*dgdr
cc      dbr1dz=0.d0
cc      dbr1dr=0.d0
cc      dbr1dph=0.d0
c      write(*,*)'dbr1dz,dbr1dr,dbr1dph',dbr1dz,dbr1dr,dbr1dph

      dbph1dz=cosphi*rnloop*dgdz/r
      dbph1dr=cosphi*rnloop*(dgdr/r-g/(r*r))
      dbph1dph=-sinphi*rnloop*rnloop*g/r
c      write(*,*)'dbph1dz,dbph1dr,dbph1dph',dbph1dz,dbph1dr,dbph1dph
      return
      end


c*************************briplmd5********************************
c  F=sin(N_loop*phi)*g(r,z)   - is a vacume ripple field potential *
c  bripl_phi(z,r,phi)=(dF/dphi)/r=cos(N_loop*phi)*g(r,z)*N_loop/r  *
c  bripl_z(z,r,phi) =(dF/dz)    =sin(N_loop*phi)*(dg/dz) 	   *
c  bripl_r(z,r,phi) =(dF/dr)    =sin(N_loop*phi)*(dg/dr) 	   *
c  model for function g:					   *
c  g_model=beqd*reqd*p1*I_0(N_loop*rho(z,r))	   		   *
c          p1=deltripl/(n_loop*I_0(n_loop))			   *
c  N is a number of the toroidal field coils			   *
c  I_0 is the zero order modifed Bessel function 		   *
c------------------------------------------------------------------
c        input parameters					   !
c      z, r, phi - point where the components of magnetic field    !
c                  are calculated in.				   !
c  nloop           is number of toroidal field coils(in one.i)     !
c  deltripl	   is an amplitude of the ripple field at the	   !
c                  last flux surface (at rho=1)(in one.i)	   !
c      bmag      - is a magnetic field on the magnetic axis (Tl)   !
c------------------------------------------------------------------
c        output parameters:					   |
c        components of the ripple magnetic field, and their	   |
c        derivatives                                               |
c        bph1,dbph1dph,dbph1dz,dbph1dr	-toroidal  (phi) components|
c        bz1,dbz1dph,dbz1dz,dbz1dr	-Z               components|
c        br1,dbr1dph,dbr1dz,dbr1dr	-R               components|
c------------------------------------------------------------------
      subroutine briplmd5(z,r,phi,
     1bph1,dbph1dph,dbph1dz,dbph1dr,
     2bz1,dbz1dph,dbz1dz,dbz1dr,
     3br1,dbr1dph,dbr1dz,dbr1dr)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      real bslr,bsli,rarg,rarg1
      dimension bslr(3),bsli(3)
cc      write(*,*)'in briplmod5'
      rnloop=dfloat(nloop)
      cosphi=dcos(phi*rnloop)
      sinphi=dsin(phi*rnloop)
cc      write(*,*)'cosphi,sinphi',cosphi,sinphi
c------------------------------------------------------------------
c  g_model=beqd*reqd*p1*I_0(N*rho(r,z))
c     izi=1 to calculate modified Bessel functions
c     BESSEL FUNCTIONS AND ARGUMENTS ARE REAL
      rarg1=rnloop
      ize=1
      nb=1
      call beslci(rarg1,0.,nb,ize,bslr,bsli,ncalc)
c      write(*,*)'rarg1,bslr(1),bslr(2),bslr(3)',
c     1 rarg1,bslr(1),bslr(2),bslr(3)
      cI_0=bslr(1)
      cI_0edg=cI_0
c      write(*,*)'cI_0,cI_0edg',cI_0,cI_0edg
      p1=deltripl/(rnloop*cI_0)
c      write(*,*)'p1',p1
c      write(*,*)'in briplmd5 phi,z,r,rho',phi,z,r,rho
      rarg=rnloop*rho
      nb=3
      call beslci(rarg,0.,nb,ize,bslr,bsli,ncalc)
c      write(*,*)'rarg,bslr(1),bslr(2),bslr(3)',
c     1 rarg,bslr(1),bslr(2),bslr(3)
      cI_0=bslr(1) ! modifed bessel function 0 order
      cI_1=bslr(2) ! modifed bessel function 1 order
      cI_2=bslr(3) ! modifed bessel function 2 order
      delcI_0=cI_0/cI_0edg
c      write(*,*)'rarg1,rarg',rarg1,rarg
c      write(*,*)'cI_0edg,cI_0',cI_0edg,cI_0,'delcI_0',delcI_0
      g=beqd*reqd*p1*cI_0
      dgdrho=beqd*reqd*p1*cI_1*rnloop
      dgdrhrh=beqd*reqd*p1*0.5d0*(cI_0+cI_2)*rnloop*rnloop
cc      write(*,*)'beqd,deltripl,rho,rnloop',beqd,deltripl,rho,rnloop
cc      write(*,*)'g,dgdrho,dgdrhrh',g,dgdrho,dgdrhrh

      psi=psid
      dro_dpsi=drhopsi(psi)
      drhodz=dro_dpsi*dpdzd
      drhodr=dro_dpsi*dpdrd
cc      write(*,*)'psi,dpdzd,dpdrd',psi,dpdzd,dpdrd
cc      write(*,*)'dro_dpsi,drhodz,drhodr',dro_dpsi,drhodz,drhodr

      dropsps=drhpsps(psi)
      drhodzdz=dropsps*dpdzd*dpdzd+dro_dpsi*dpdzzd
      drhodrdr=dropsps*dpdrd*dpdrd+dro_dpsi*dpdrrd
      drhodzdr=dropsps*dpdzd*dpdrd+dro_dpsi*dpdzrd

      dgdz=dgdrho*drhodz
      dgdr=dgdrho*drhodr

      dgdzdz=dgdrhrh*drhodz*drhodz+dgdrho*drhodzdz
      dgdrdr=dgdrhrh*drhodr*drhodr+dgdrho*drhodrdr
      dgdzdr=dgdrhrh*drhodz*drhodr+dgdrho*drhodzdr
c------------------------------------------------------------------
      bz1=sinphi*dgdz
      br1=sinphi*dgdr
      bph1=cosphi*g*rnloop/r
cc      bz1=0.d0
cc      br1=0.d0
c      write(*,*)'bz1,br1,bph1',bz1,br1,bph1
c       write(*,*)'cosphi,g,bph1',cosphi,g,bph1
       ppp=beqd*reqd/r
c       write(*,*)'beqd*reqd/r',ppp,'bph1/bphi',bph1/ppp
c       write(*,*)'cosphi,deltripl,delcI_0',cosphi,deltripl,delcI_0
c       write(*,*)'delta test=',cosphi*deltripl*delcI_0
c------------------------------------------------------------------
      dbz1dz=sinphi*dgdzdz
      dbz1dr=sinphi*dgdzdr
      dbz1dph=cosphi*rnloop*dgdz
cc      dbz1dz=0.d0
cc      dbz1dr=0.d0
cc      dbz1dph=0.d0
c      write(*,*)'dbz1dz,dbz1dr,dbz1dph',dbz1dz,dbz1dr,dbz1dph

      dbr1dz=sinphi*dgdzdr
      dbr1dr=sinphi*dgdrdr
      dbr1dph=cosphi*rnloop*dgdr
cc      dbr1dz=0.d0
cc      dbr1dr=0.d0
cc      dbr1dph=0.d0
c      write(*,*)'dbr1dz,dbr1dr,dbr1dph',dbr1dz,dbr1dr,dbr1dph

      dbph1dz=cosphi*rnloop*dgdz/r
      dbph1dr=cosphi*rnloop*(dgdr/r-g/(r*r))
      dbph1dph=-sinphi*rnloop*rnloop*g/r
c      write(*,*)'dbph1dz,dbph1dr,dbph1dph',dbph1dz,dbph1dr,dbph1dph
      return
      end

c*************************bstelw7********************************
c     calculates  the model magnetic field for Stellarator W7-AS

cc------------------------------------------------------------------
c        input parameters				          !
c      z, r,     - point where the components of magnetic field   !
c                  are calculated in.			          !
c      bmagloc   - is a magnetic field on the magnetic axis (Tl)  !
c------------------------------------------------------------------
c        output parameters:			                  |
c        components of the stellarator magnetic field, and their  |
c        derivatives                                              |
c        bph1,dbph1dph,dbph1dz,dbph1dr	-toroidal  (phi) components|
c        bz1,dbz1dph,dbz1dz,dbz1dr	-Z               components|
c        br1,dbr1dph,dbr1dz,dbr1dr	-R               components|
c------------------------------------------------------------------
      subroutine bstelw7(z,r,bmagloc,
     1bph1,dbph1dph,dbph1dz,dbph1dr,
     2bz1,dbz1dph,dbz1dz,dbz1dr,
     3br1,dbr1dph,dbr1dz,dbr1dr)
      implicit double precision (a-h,o-z)
      include 'param.i'
c      include 'one.i'
      include 'three.i'
      include 'five.i'
c-----input
      double precision z,r,bmagloc
c-----output
      double precision bph1,dbph1dph,dbph1dz,dbph1dr,
     2 bz1,dbz1dph,dbz1dz,dbz1dr,
     3 br1,dbr1dph,dbr1dz,dbr1dr
c-----local
      double precision  abig,asmall,reff

c     write(*,*)'in bstelw7'
      abig=10.5 d0
      asmall=rmax-xma
      reff=(r-xma)
c      write(*,*)'in bstelw7 z,r',z,r
c      write(*,*)'in bstelw7 asmal,reff',asmall,reff

      bz1=0.d0
      br1=0.d0
      bph1=bmagloc*abig/(abig+reff/asmall)
c      write(*,*)'bz1,br1,bph1',bz1,br1,bph1
c------------------------------------------------------------------
      dbz1dz=0.d0
      dbz1dr=0.d0
      dbz1dph=0.d0
c      write(*,*)'dbz1dz,dbz1dr,dbz1dph',dbz1dz,dbz1dr,dbz1dph

      dbr1dz=0.d0
      dbr1dr=0.d0
      dbr1dph=0.d0
c      write(*,*)'dbr1dz,dbr1dr,dbr1dph',dbr1dz,dbr1dr,dbr1dph

      dbph1dz=0.d0
      dbph1dr=-bph1/(abig+reff/asmall)/asmall
      dbph1dph=0.d0

c      write(*,*)'dbph1dz,dbph1dr,dbph1dph',dbph1dz,dbph1dr,dbph1dph
      return
      end


c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        ************************ b ***************************     
c        this function calculates the components  of  the   
c        magnetic field bz,br,bphi, 
c        the absolute value of the total magnetic field bmod
c        and the derivatives of bz,br,bphi,bmod by rz,r,phi 
c        It calculates the rho,rho2,a for the programms     
c        which use the density profile :dense,dxdr,dxdz     
c	 						      
c        It controls if the ray point is inside the plasma  
c        if the ray point is outside the plasma then the    
c        calculations will be finished 
c        It uses spline functions from zcuinx coeff2 coeff1
c        terp2p, terp1		      
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c------------------------------------------------------------------
      real*8
     &function b1test(z,r,phi) !YuP[2020-08-20] renamed, to avoid conflict with arrays
      implicit none 
      !implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'six.i'
      include 'fourb.i'
c-----input
      real*8 z,r,phi !sapce coordinates
            
c-----external 
      real*8
     1 ias1r,ias2r_Sm

c-----locals
   
cSm030120
      logical first
      real*8 f_pp(nxeqda),tab(3),wrk_f(3*nxeqda+1)
      integer iop(2),itab(3),ibd(4),idm
c-------------------------------------------------------
cSAP090328     
      integer nx4,ny4,idx,ipx4,ipx,ifirstc

      real*8
     &bph1,dbph1dph,dbph1dz,dbph1dr,
     &bz1,dbz1dph,dbz1dz,dbz1dr,
     &br1,dbr1dph,dbr1dz,dbr1dr,
     &ffd,ffr,dffr,dres,dpdzr,dpdrr,dpsidz,dpsidr,dpdzz,
     &epsbnd,l_zr,l_z_b_r_b,x_l,y_l,theta_l,z_b,r_b,
     &pb,ppb,res,rrr,spol,zer,pp,zzrm,zzrp,zm,zp
c-----external
      real*8 thetapol,rhopsi
c----------------------------------------------------------------
      real*8 peqd_rr(nxeqda,nyeqda),peqd_zz(nxeqda,nyeqda),
     &peqd_rrzz(nxeqda,nyeqda),wrk_peqd(3*nxeqda+1)

      real*8 terp2p,terp2p_Sm !external functions
      
      data first /.true./
      data iop /4,4/
      data itab/1,1,1/
      data ibd/4,4,4,4/
       
cend

      data ifirstc/1/
      save ifirstc


cSm030120
      save f_pp,peqd_rr,peqd_zz,peqd_rrzz,idm
c      write(*,*)'b1 first',first
cSm030220
       nx=nxeqd
       ny=nyeqd
       nx4=nx+4
       ny4=ny+4
       if (first) then
        first=.false.
c-------feqd coefficient calculations
cSm030221
c        call coeff1(nxmax,flux,feqd,f_pp,iop,1,wrk_f)
c        write(*,*)'b1 before coeff1'
        call coeff1(nx,flux,feqd,f_pp,iop,1,wrk_f)
c-------psi coefficient calculations
cSm030221
c        idm=nxmax
         idm=nx
cSm040115
         idm=nxeqda

c        write(*,*)' b before coeff2'
c        write(*,*)'rr',rr
c        write(*,*)'zz',zz

cSm030221
c        call coeff2(nxmax,rr,nymax,zz,peqd,peqd_rr,peqd_zz,peqd_rrzz,
c     &               idm,ibd,wrk_peqd)
        call coeff2_Sm(nx,rr,ny,zz,peqd,peqd_rr,peqd_zz,peqd_rrzz,
     &               idm,ibd,wrk_peqd,nxeqda)

c        write(*,*)' b after coeff2 nx,ny,idm',nx,ny,idm
c        write(*,*)'peqd_rr',peqd_rr
c        write(*,*)'peqd_zz',peqd_zz
c        write(*,*)'peqd_rrzz',peqd_rrzz
c        write(*,*)'ibd',ibd
ctest  
c       do i=1,nx
c          r=rr(i)
c          do j=1,ny
c            z=zz(j)
c            res=terp2p_Sm(r,z,nx,rr,ny,zz,peqd,peqd_rr,peqd_zz,
c     &      peqd_rrzz,idm,0,0,1,nxeqda)
c            res=terp2p_Sm(r,z,nx,rr,ny,zz,peqd,peqd_rr,peqd_zz,
c     &      peqd_rrzz,idm,0,0,1,nxeqda)
c            write(*,*)'i,j,r,z,peqd(i,j),res',
c     &      i,j,r,z,peqd(i,j),res
c          enddo
c        enddo
c        !stop 'b.f in b1'

       endif
c--------------------------------------------------------------
c     control that the point is inside the plasma
      epsbnd=1.d-10
      iboundb=-1
c-----------------------   
      x_l=r-xma
      y_l=z-yma
c      write(*,*)'b1.f x_l,y_l,theta_l', x_l,y_l,theta_l
      if ((x_l**2+y_l**2).lt.1.d-10) then
         goto 10
      endif

c-----calculate poloidal angle  thetapol
      call theta_xy(x_l,y_l,theta_l) 

c      write(*,*)'b1.f x_l,y_l,theta_l', x_l,y_l,theta_l

c      write(*,*)'b1.f theta_l',theta_l
c-----calculate coordinates x_b,z_b at the limiter surface psi=psilim
c     and the poloidal angle thetapol
      call zr_psith(psilim,theta_l,z_b,r_b)

c      write(*,*)'b1.f theta_l',theta_l
c-----calculate coordinates x_b,z_b at the limiter surface psi=psilim
c     and the poloidal angle thetapol
      call zr_psith(psilim,theta_l,z_b,r_b)
 
      
c      write(*,*)'b1.f r_b,z_b,fpsi(r_b,z_b),psilim',
c     &r_b,z_b,fpsi(r_b,z_b),psilim

      l_zr=dsqrt(x_l**2+y_l**2)
      l_z_b_r_b=dsqrt((r_b-xma)**2+(z_b-yma)**2)
c      write(*,*)'b1.f l_zr,l_z_b_r_b',l_zr,l_z_b_r_b
      if (l_zr.ge.l_z_b_r_b) then
         rho = l_zr / l_z_b_r_b
c        drhodzb=(z-yma)/(((r_b-xma)**2+(z_b-yma)**2)*rho)

         iboundb=3
      endif
      goto 10
c----------------------------------------------------

      if ((r.le.rmin+epsbnd).or.(r.ge.rmax-epsbnd)) then
c         write(*,*)'in b rmin+epsbnd,r,rmax-epsbnd',
c     1   rmin+epsbnd,r,rmax-epsbnd
         iboundb=2
c	 delr=0.5d0*(rmax-rmin)
c         if (r.le.rmin+epsbnd) rho=1.d0+(rmin-r+epsbnd)/delr
c         if (r.ge.rmax-epsbnd) rho=1.d0+(r-rmax+epsbnd)/delr
cSm070324 
c         write(*,*)'b.f rho',rho
         theta_l=thetapol(z,r) 
         call zr_psith(psilim,theta_l,z_b,r_b)
         rho=dsqrt(((r-xma)**2+(z-yma)**2)/((r_b-xma)**2+(z_b-yma)**2))
         goto 10
      end if
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      rrr=r
      idx=0
      ipx=ip
      ipx4=ip+4
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
      zm=zzrm
      if ((z.ge.zp-epsbnd).or.(z.le.zm+epsbnd)) then
c         write(*,*)'in b rrr,zm,z,zp',rrr,zm,z,zp
c         write(*,*)'in b zm+epsbnd,z,zp-epsbnd',
c     1	 zm+epsbnd,z,zp-epsbnd
c	 delz=0.5d0*(zmax-zmin)
c         if (z.ge.zp-epsbnd) then
c	    rho=1.d0+(z-zp+epsbnd)/delz
c	    drhodzb=1/delz
c	 else
c            rho=1.d0+(zm-z+epsbnd)/delz
c	    drhodzb=-1/delz
c	 endif
         iboundb=1
cSm070324 
c         write(*,*)'b.f rho',rho
         theta_l=thetapol(z,r) 
         call zr_psith(psilim,theta_l,z_b,r_b)
         rho=dsqrt(((r-xma)**2+(z-yma)**2)/((r_b-xma)**2+(z_b-yma)**2))
         drhodzb=(z-yma)/(((r_b-xma)**2+(z_b-yma)**2)*rho)

         goto 10
      end if
c---------------------------------------------------------------------
 10     continue
c end of the control that the point is inside the plasma
c--------------------------------------------------------------
c
c nx4a , ny4a and nrya  are given in common five as parameters
c nxeqda,nyeqda,nlimit are given in common four as parameters
cSm030224
      nx4=nx+4
      ny4=ny+4

      ncx=nx4
      ncy=ny4

c_old spline
c      res=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
c      dpsidr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z)
c      dpsidz=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z)

cSm030120
c new spline
      res=terp2p_Sm(r,z,nx,rr,ny,zz,peqd,peqd_rr,peqd_zz,
     &      peqd_rrzz,idm,0,0,1,nxeqda)
      dpsidr=terp2p_Sm(r,z,nx,rr,ny,zz,peqd,peqd_rr,peqd_zz,
     &      peqd_rrzz,idm,1,0,0,nxeqda)
      dpsidz=terp2p_Sm(r,z,nx,rr,ny,zz,peqd,peqd_rr,peqd_zz,
     &      peqd_rrzz,idm,0,1,0,nxeqda)
c--------------------------------------------------------------------
c     idx derivativs order 0.ge.idx.le.3
c--------------------------------------------------------------------
     
      if (iboundb.ge.1) then
c        the point is outside the plasma
         idx=0
         ffr=ias1r(txf,nx,nx4,cx,idx,psilim)         
         ffd=ffr
         dffr=0.d0
         dres=dffr
      else
c        the point is inside the plasma

c old spline
c         idx=0
c         ffr=ias1r(txf,nx,nx4,cx,idx,res)

c---------------
cSm030120
c new spline
         call terp1(nx,flux,feqd,f_pp,res,1,tab,itab)
         ffr=tab(1)
c--------------         

         ffd=ffr
       
c old spline
c        idx=1
c         dffr=ias1r(txf,nx,nx4,cx,idx,res)

c--------------
c new spline
         dffr=tab(2)         
c--------------

         dres=dffr         
      endif

      psid=res
      dpdzd=dpsidz
      dpdrd=dpsidr
c------------------------------------------------------------------
c this  part of the programm calculates the rho,rho2 and a
c for the functions dense,dxdr,dxdz

       a=1.d0
       if (iboundb.eq.-1) rho=rhopsi(res)
       zer=0.d0
       rho2=rho*rho
c------------------------------------------------------------------
c     \vector{B}=B_phi+\vector{B_pol}
c     B_phi=(Feqd(psi)/R)\hat{phi}
c     \vector{B_pol}=(1/R)*grad(psi) x \hat{phi}
c     here x is the vector product
c     then  B_z=(1/R)*dpsi/dr and B_r=-(1/R)*dpsi/dz  (*).
c     If the total current is not given, then the sign of the B_pol
c     determined fro psi, using the previous formula (*).
c     If the current 'toteqd'is given, then the sign of it will indicate 
c     the direction of B_pol, regardless of the convention for psi(R,Z)
c     toteqd is positive if it is directed along \hat{phi}.
c-------------------------------------------------------------------
c     The poloidal flux	array peqd(i,j) was transformed in equilib.f
c     in subroutine input. It was multiplied by dpsimax.
c     If (psimag.gt.psilim) then dpsimax=-1.d0 else  dpsimax=1.d0
c     So after this multiplication poloidal flux peqd(i,j)
c     always has the minimum on tme magnetic axis.
c-------------------------------------------------------------------
      pp=1.d0/r
      bz=-dpdrd*pp
      br=dpdzd*pp
      bphi=ffd*pp
      
	if(toteqd.ge.0.d0) spol=-1.d0
	if(toteqd.lt.0.d0) spol=1.d0
	spol=-spol
	if(toteqd.eq.0.d0) then
           if (ifirstc.eq.1) then
              write(*,*)'current is not given in eqdsk (toteqd=0)' 
              ifirstc=0
           endif
c          The total current is not given. The sign of the B_pol will be
c          determined from psi, using the original poloidal flux from eqdsk.
c          In this case the poloidal flux had max on the magnetic axis
c          (and dpsimax=-1d0), or it had min (and dpsimax=1.d0).
           spol=dpsimax
	endif
      
	bz=bz*spol
	br=br*spol

c old spline
c              dpdzr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,1,1,r,z)
c              dpdrr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,2,0,r,z)
c              dpdzz=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,2,r,z)
           
c--------------------
cSm030120
c            new spline
            dpdzr=terp2p_Sm(r,z,nx,rr,ny,zz,peqd,peqd_rr,peqd_zz,
     &      peqd_rrzz,idm,1,1,1,nxeqda)
c           write(*,*)'before 2 terp2'
            dpdrr=terp2p_Sm(r,z,nx,rr,ny,zz,peqd,peqd_rr,peqd_zz,
     &      peqd_rrzz,idm,2,0,1,nxeqda)
            dpdzz=terp2p_Sm(r,z,nx,rr,ny,zz,peqd,peqd_rr,peqd_zz,
     &      peqd_rrzz,idm,0,2,1,nxeqda)
c--------------------

		dpdzrd=dpdzr
		dpdrrd=dpdrr
		dpdzzd=dpdzz
c------------------------------------------------------------------
      dbzdz=-pp*dpdzrd
      dbzdr=pp*pp*dpdrd-pp*dpdrrd
      dbrdz=pp*dpdzzd
      dbrdr=-pp*pp*dpdzd+pp*dpdzrd
	dbzdz=dbzdz*spol
	dbzdr=dbzdr*spol
	dbrdz=dbrdz*spol
	dbrdr=dbrdr*spol
      dbzdph=0.d0
      dbrdph=0.d0
      dbpdph=0.d0
      dbphdz=pp*dres*dpdzd
      dbphdr=-pp*bphi+pp*dres*dpdrd
c--------------------------------------------------------------------
c      deltripl=0.00d0

      if (i_ripple.eq.1)then
c--------the ripple model approximating the DIII-D field
         call briplmd4(z,r,phi,
     +   bph1,dbph1dph,dbph1dz,dbph1dr,
     +   bz1,dbz1dph,dbz1dz,dbz1dr,
     +   br1,dbr1dph,dbr1dz,dbr1dr)
      endif

      if (i_ripple.eq.2)then
c--------the ripple model using modified Bessel function I_0
c        see the GENRAY description
         call briplmd5(z,r,phi,
     +   bph1,dbph1dph,dbph1dz,dbph1dr,
     +   bz1,dbz1dph,dbz1dz,dbz1dr,
     +   br1,dbr1dph,dbr1dz,dbr1dr)
      endif
 
      bz=bz+bz1
      br=br+br1
      bphi=bphi+bph1

      dbzdz=dbzdz+dbz1dz
      dbzdr=dbzdr+dbz1dr
      dbzdph=dbzdph+dbz1dph

      dbrdz=dbrdz+dbr1dz
      dbrdr=dbrdr+dbr1dr
      dbrdph=dbrdph+dbr1dph

      dbpdph=dbpdph+dbph1dph
      dbphdz=dbphdz+dbph1dz
      dbphdr=dbphdr+dbph1dr
c-----------------------------------------------------------------
c     calculation of the model stellarator W7 magnetic field
c     bmagloc=2.3d0 !Tl

c      call bstelw7(z,r,bmagloc,
c     1bphi,dbpdph,dbphdz,dbphdr,
c     2bz,dbzdph,dbzdz,dbzdr,
c     3br,dbrdph,dbrdz,dbrdr)

c-----------------------------------------------------------------

c------------------------------------------------------------------
      pb=dsqrt(bphi*bphi+bz*bz+br*br)
c------------------------------------------------------------------
	ppb=1.d0/pb
	dbmdz=ppb*(bz*dbzdz+br*dbrdz+bphi*dbphdz)
	dbmdr=ppb*(bz*dbzdr+br*dbrdr+bphi*dbphdr)
	dbmdph=ppb*(bz*dbzdph+br*dbrdph+bphi*dbpdph)
c---------------------------------------------------------------------
      b1test=pb !YuP[2020-08-20] renamed, to avoid conflict with arrays
      return
      end function b1test

      subroutine test_b_field
c-----compare the calculation of the magneic field 
c     obtained by the different spilnes
      !implicit double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'fourb.i'
      include 'five.i'

      integer i,j 
      real*8 r,z,bmod1,bmod2
     
      real*8
     &bz1,br1,bphi1,psid1,dpdzd1,dpdrd1,dpdzzd1,dpdrrd1,dpdzrd1,
     3            dbzdz1,dbzdr1,dbzdph1,dbrdz1,dbrdr1,dbrdph1,
     4            dbphdz1,dbphdr1,dbpdph1,
     5            dqdrho1,dbmdz1,dbmdr1,dbmdph1,
     &bz2,br2,bphi2,psid2,dpdzd2,dpdrd2,dpdzzd2,dpdrrd2,dpdzrd2,
     3            dbzdz2,dbzdr2,dbzdph2,dbrdz2,dbrdr2,dbrdph2,
     4            dbphdz2,dbphdr2,dbpdph2,
     5            dqdrho2,dbmdz2,dbmdr2,dbmdph2,
     &bz3,br3,bphi3,psid3,dpdzd3,dpdrd3,dpdzzd3,dpdrrd3,dpdzrd3,
     3            dbzdz3,dbzdr3,dbzdph3,dbrdz3,dbrdr3,dbrdph3,
     4            dbphdz3,dbphdr3,dbpdph3,
     5            dqdrho3,dbmdz3,dbmdr3,dbmdph3
c-----externals
      real*8 b,b1test

            bz3=0.d0
            br3=0.d0
            bphi3=0.d0
            psid3=0.d0
            dpdzd3=0.d0
            dpdrd3=0.d0
            dpdzzd3=0.d0
            dpdrrd3=0.d0
            dpdzrd3=0.d0
            dbzdz3=0.d0
            dbzdr3=0.d0
            dbzdph3=0.d0
            dbrdz3=0.d0
            dbrdr3=0.d0
            dbrdph3=0.d0
            dbphdz3=0.d0
            dbphdr3=0.d0
            dbpdph3=0.d0
            dqdrho3=0.d0
            dbmdz3=0.d0
            dbmdr3=0.d0
            dbmdph3=0.d0

      do i=1,nx
         r=rr(i)
         do j=1,ny
            z=zz(j)
            bmod1=b(z,r,0.d0)
            bz1=bz
            br1=br
            bphi1=bphi
            psid1=psid
            dpdzd1=dpdzd
            dpdrd1=dpdrd
            dpdzzd1=dpdzzd
            dpdrrd1=dpdrrd
            dpdzrd1=dpdzrd
            dbzdz1=dbzdz
            dbzdr1=dbzdr
            dbzdph1=dbzdph
            dbrdz1=dbrdz
            dbrdr1=dbrdr
            dbrdph1=dbrdph
            dbphdz1=dbphdz
            dbphdr1=dbphdr
            dbpdph1=dbpdph
            dqdrho1=dqdrho
            dbmdz1= dbmdz
            dbmdr1=dbmdr
            dbmdph1=dbmdph

            bmod2=b1test(z,r,0.d0) !YuP[2020-08-20] renamed, to avoid conflict with arrays
            bz2=bz
            br2=br
            bphi2=bphi
            psid2=psid
            dpdzd2=dpdzd
            dpdrd2=dpdrd
            dpdzzd2=dpdzzd
            dpdrrd2=dpdrrd
            dpdzrd2=dpdzrd
            dbzdz2=dbzdz
            dbzdr2=dbzdr
            dbzdph2=dbzdph
            dbrdz2=dbrdz
            dbrdr2=dbrdr
            dbrdph2=dbrdph
            dbphdz2=dbphdz
            dbphdr2=dbphdr
            dbpdph2=dbpdph
            dqdrho2=dqdrho
            dbmdz2= dbmdz
            dbmdr2=dbmdr
            dbmdph2=dbmdph


            bz3=bz3+abs(bz2-bz1)
            br3=br3+abs(br2-br1)
            bphi3=bphi3+abs(bphi2-bphi1)
            psid3=psid3+abs(psid2-psid1)
            dpdzd3=dpdzd3+abs(dpdzd2-dpdzd1)
            dpdrd3=dpdrd3+abs(dpdrd2-dpdrd1)
            dpdzzd3=dpdzzd3+abs(dpdzzd2-dpdzzd1)
            dpdrrd3=dpdrrd3+abs(dpdrrd2-dpdrrd1)
            dpdzrd3=dpdzrd3+abs(dpdzrd2-dpdzrd1)
            dbzdz3=dbzdz3+abs(dbzdz2-dbzdz1)
            dbzdr3=dbzdr3+abs(dbzdr2-dbzdr1)
            dbzdph3=dbzdph3+abs(dbzdph2-dbzdph1)
            dbrdz3=dbrdz3+abs(dbrdz2-dbrdz1)
            dbrdr3=dbrdr3+abs(dbrdr2-dbrdr1)
            dbrdph3=dbrdph3+abs(dbrdph2-dbrdph1)
            dbphdz3=dbphdz3+abs(dbphdz2-dbphdz1)
            dbphdr3=dbphdr3+abs(dbphdr2-dbphdr1)
            dbpdph3=dbpdph3+abs(dbpdph2-dbpdph1)
            dqdrho3=dqdrho3+abs(dqdrho2-dqdrho1)
            dbmdz3=dbmdz3+abs(dbmdz2-dbmdz1)
            dbmdr3=dbmdr3+abs(dbmdr2-dbmdr1)
            dbmdph3=dbmdph3+abs(dbmdph2-dbmdph1)
            write(*,*)'i,j,r,z,peqd(i,j),psid1,psid2',
     &      i,j,r,z,peqd(i,j),psid1,psid2
            write(*,*)'bz3,bz1,bz2',bz3,bz1,bz2
            write(*,*)'br3,br1,br2',br3,br1,br2
            write(*,*)'bphi3,bphi1,bphi2',bphi3,bphi1,bphi2
            write(*,*)'psid3,psid1,psid2',psid3,psid1,psid2
            write(*,*)'dpdzd3,dpdzd1,dpdzd2',dpdzd3,dpdzd1,dpdzd2
            write(*,*)'dpdrd3,dpdrd1,dpdrd2',dpdrd3,dpdrd1,dpdrd2
            write(*,*)'dpdzzd3,dpdzzd1,dpdzzd2',dpdzzd3,dpdzzd1,dpdzzd2
            write(*,*)'dpdrrd3,dpdrrd1,dpdrrd2',dpdrrd3,dpdrrd1,dpdrrd2
            write(*,*)'dpdzrd3,dpdzrd1,dpdzrd2',dpdzrd3,dpdzrd1,dpdzrd2
            write(*,*)'dbzdz3,dbzdz1,dbzdz2',dbzdz3,dbzdz1,dbzdz2
            write(*,*)'dbzdr3,dbzdr1,dbzdr2',dbzdr3,dbzdr1,dbzdr2
            write(*,*)'dbzdph3,dbzdph1,dbzdph2',dbzdph3,dbzdph1,dbzdph2
            write(*,*)'dbrdz3,dbrdz1,dbrdz2',dbrdz3,dbrdz1,dbrdz2
            write(*,*)'dbrdr3,dbrdr1,dbrdr2',dbrdr3,dbrdr1,dbrdr2
            write(*,*)'dbrdph3,dbrdph1,dbrdph2',dbrdph3,dbrdph1,dbrdph2
            write(*,*)'dbphdz3,dbphdz1,dbphdz2',dbphdz3,dbphdz1,dbphdz2
            write(*,*)'dbphdr3,dbphdr1,dbphdr2',dbphdr3,dbphdr1,dbphdr2
            write(*,*)'dbpdph3,dbpdph1,dbpdph2',dbpdph3,dbpdph1,dbpdph2
            write(*,*)'dqdrho3,dqdrho1,dqdrho2',dqdrho3
            write(*,*)'dbmdz3,dbmdz1,dbmdz2',dbmdz3,dbmdz1,dbmdz2
            write(*,*)'dbmdr3,dbmdr1,dbmdr2',dbmdr3,dbmdr1,dbmdr2
            write(*,*)'dbmdph3,dbmdph1,dbmdph2',dbmdph3,dbmdph1,dbmdph2
          enddo
      enddo
            write(*,*)'bz3',bz3
            write(*,*)'br3',br3
            write(*,*)'bphi3',bphi3
            write(*,*)'psid3',psid3
            write(*,*)'dpdzd3',dpdzd3
            write(*,*)'dpdrd3',dpdrd3
            write(*,*)'dpdzzd3',dpdzzd3
            write(*,*)'dpdrrd3',dpdrrd3
            write(*,*)'dpdzrd3',dpdzrd3
            write(*,*)'dbzdz3=',dbzdz3
            write(*,*)'dbzdr3=',dbzdr3
            write(*,*)'dbzdph3',dbzdph3
            write(*,*)'dbrdz3=',dbrdz3
            write(*,*)'dbrdr3=',dbrdr3
            write(*,*)'dbrdph3',dbrdph3
            write(*,*)'dbphdz3',dbphdz3
            write(*,*)'dbphdr3',dbphdr3
            write(*,*)'dbpdph3',dbpdph3
            write(*,*)'dqdrho3',dqdrho3
            write(*,*)'dbmdz3',dbmdz3
            write(*,*)'dbmdr3',dbmdr3
            write(*,*)'dbmdph3',dbmdph3
      stop 'b.f in test_b_field'

      return
      end
! 
    

