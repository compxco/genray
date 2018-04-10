c        ********************** dddrz1***********************
c        *                      -----                       *
c        * this subroutine calculates the derivatives
c        * from Hamiltonian  				    *
c        * for 6 ray equations                              *
c        * The Hamiltonian derivatives  are calculated	    *
c        * analytically(idif=1) or numerically (idif=2)	    *
c        * This subroutine is used for the Hamiltonian
c        * correction precedure
c        ****************************************************
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      t - parameter of trajectory                                  !
c                                  		         	    !
c      u - solution of geometrical optics equations at the point t  !
c	u(1) = z                                                    !
c	u(2) = r                                                    !
c	u(3) = phi                                                  !
c	u(4) = n_z                                                  !
c	u(5) = n_r                                                  !
c	u(6) = r*n_phi  					    !
c        output parameters					    !
c                                                                   !
c     deru(i) are right  hand  sides  of  geometrical  optics       !
c     equations							    !
c      n_z     - deru(1)                                            !
c      n_r     - deru(2)                                            !
c      r*n_phi - deru(3)                                            !
c      z       - deru(4)                                            !
c      r       - deru(5)                                            !
c      phi     - deru(6)                                            !
c     in fact deru are right  hand  sides  of  geometrical  optics  !
c     equation                                                      !
c-------------------------------------------------------------------
c     this program uses following functions and subroutines	    !
c     b,gamma1,dxdz,dxdr,dxdphi,dydz,dydr,dydphi,x,y,s,tensor,	    !
c     hamilt1, etc.    						    !
c-------------------------------------------------------------------
      subroutine dddrz1(t,u,deru)
      implicit none
c      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
c this line only for call control
      include 'eps.i'

c-----input
      real*8 t,
     &u(*)    !ray coordinates
c-----output
      real*8 
     &deru(*) !derivatives d_D/d_

c-----local
      real*8 vp(nbulka),wp(nbulka),
     &z,r,phi,cnz,cnr,cm,wf,
     &step,hz,hr,hphi,hnz,hnr,hms,hfrqnc,
     &cnzplus,cnzmins,hp,hm,dddcnz,
     &cnrplus,cnrmins,dddcnr,
     &cmplus,cmminus,hmin,dddcm,
     &zplus,zminus,dddz,
     &rplus,rminus,dddr,
     &phipls,phimin,dddphi,
     &frqncpl,df,
     &frqncmn,cnrminus,cnzminus,
     &dddw,p,
     &ds,dc,ds2,dc2,
     &dxidz,dxidr,dxidph,dyidz,dyidr,dyidphi,
     &xi,yi,py2,py3,py4,px,px2,ds4,det,zer,sqrdet,p1,
     &ddetdx,ddetdy,pz,pz2,dpzdx,dpzdy,dddx,dddy,dddc2,
     &pc,s1,s2,s3,s4,s6,s7,xib,yib,delib,xe,ye,
     &peym,peyp,a0e,a1e,b0e,b1e,c0e,c1e,dele,ad,bd,cd,
     &da1edc,da0edc,db1edc,db0edc,dc1edc,dc0edc,
     &dadc2,dbdc2,dcdc2,
     &ppe,peym2,
     &ppb,pbym,pbyp,pbym2,
     &a0b,a1b,b0b,b1b,c0b,c1b,
     &da1bdc,da0bdc,db1bdc,db0bdc,dc1bdc,
     &dadz,dadr,dadphi,
     &dbdz,dbdr,dbdphi,
     &dcdz,dcdr,dcdphi,
     &dadw,dbdw,dcdwm,
     &pyim,pyip,pyim2,
     &dc1bye, db0byi,da1eye,da1exi,da1byb,da0bxi,da0bxb,
     &da0bxe,da0exe,da0byi,da0byb,da0bye,da1bxb,
     &da0eye,da0exi,da0eyi,
     &da1bxe,da1bxi,da1exe,da1bye,da1byi, 
     &dadnz,da1eyi,dadm,dadnr,ds3dyi,dbdyi,db0bxb,dadxi,
     &dadyi,db0bxi,db0bxe,db0bye,db0byb,dbdm,db1eyi,db0exi, 
     &db0exe, db1exe,db1byb,db1bxe,db0eye,db0eyi,db1bxb,db1bxi,
     &db1bye,db1byi,db1exi,db1eye,dbdnr,ds4dxe,dc0bxb,dbdnz,
     &dbdxi,dc0bdc,dc1exe,dc0eye,dc0bxe,dc0exi,dc0byi,dc0byb,
     &dc0bxi,dc0bye,dc0exe,dc0eyi,dc1byb,dc1bxe,dc1bxb,dc1bxi,
     &dc1byi,dcdnr,dc1exi,dc1eyi,dc1eye,dcdm,
     &dcn4,dcdnz,dcn2,dcdw,dcdxi,dcdyi,dcn,dddph,dddn,
     &ds1dye,ds1dxb,ds1dyb,ds1dxe,ds1dxi,ds1dyi,
     &ds6dy,ds2dxb,ds4dyb,ds2dye,ds2dyb,ds2dxe,ds2dxi,ds3dxb,
     &ds2dyi,ds3dxe,ds3dyb,ds3dxi,ds3dye,ds4dxb,ds4dxi,ds6dxb,
     &ds4dyi,ds6dye,ds6dxe,ds4dye,ds6dxi,ds6dyb,ds6dyi,ds7dyb,
     &ds7dxb,ds7dxe,ds7dxi,ds7dye,ds7dyi,dyidph,hw,
     &p7,p6,p4,p5

      integer i
c-----externals
      real*8 b,gamma1,hamilt1,dxdz,dxdr,dxdphi,dydz,dydr,dydphi,x,y,
     *cn 
         z=u(1)
         r=u(2)
         phi=u(3)
         cnz=u(4)
         cnr=u(5)
         cm=u(6)

c         write(*,*)'in dddrz1 z,r,phi,cnz,cnr,cm',
c     1                        z,r,phi,cnz,cnr,cm
c         write(*,*)'dddrz1 rho',rho,'idif',idif
c-------------------------------------------------------
          wf=frqncy
c-------------------for analytical derivatives-----------
c	 idif=1
c-------------------for numerical derivatives------------
c	 idif=2
c------------------------------------------------------------------
      if (idif.eq.1) goto 1955

c----------------------------------------------------------------
c numerical hamiltonian derivative calculations
c----------------------------------------------------------------
           step=1.d-7
c           step=1.d-8
c     This quantity was adjusted to 1.d-3, 040807, to give
c     a smooth fluxn calculation with the Westerhof-Tokman
c     integration scheme, id=10, per Smirnov.
           step=1.d-3
         

	   hz=step
	   hr=step
	   hphi=step
	   hnz=step
	   hnr=step
	   hms=step
	   hw=step*1.0d0
c	   pi=3.1415926d0
c	   hfrqnc=hw/(2.0d0*pi)
           hfrqnc=hw
	   hfrqnc=hfrqnc*frqncy
c----------------------------------------------------------------
	   cnzplus=cnz+hnz
	   cnzmins=cnz-hnz
	   bmod=b(z,r,phi)
	   gam=gamma1(z,r,phi,cnzplus,cnr,cm)
	   hp=hamilt1(z,r,phi,cnzplus,cnr,cm)
	   gam=gamma1(z,r,phi,cnzmins,cnr,cm)
	   hm=hamilt1(z,r,phi,cnzmins,cnr,cm)  
	   dddcnz=(hp-
     1  	   hm)/(2.d0*hnz)
           deru(1)=dddcnz
     

	   cnrplus=cnr+hnr
	   cnrmins=cnr-hnr
	   gam=gamma1(z,r,phi,cnz,cnrplus,cm)
	   hp=hamilt1(z,r,phi,cnz,cnrplus,cm)
	   gam=gamma1(z,r,phi,cnz,cnrmins,cm)
	   hm=hamilt1(z,r,phi,cnz,cnrmins,cm)
	   dddcnr=(hp-
     2  	   hm)/(2.d0*hnr)
           deru(2)=dddcnr


	   cmplus=cm+hms
	   cmminus=cm-hms
	   gam=gamma1(z,r,phi,cnz,cnr,cmplus)
	   hp=hamilt1(z,r,phi,cnz,cnr,cmplus)
	   gam=gamma1(z,r,phi,cnz,cnr,cmminus)
	   hmin=hamilt1(z,r,phi,cnz,cnr,cmminus)
	   dddcm=(hp-
     3  	   hmin)/(2.d0*hms)
           deru(3)=dddcm
c           write(*,*)'rside1 cmplus,hp',cnrplus,hp
c           write(*,*)'rside1 cmminus,hm',cmminus,hm
c           write(*,*)'rside1 dddcm,deru(3)',dddcm,deru(3)  


	   zplus=z+hz
	   zminus=z-hz
	   bmod=b(zplus,r,phi)
	   gam=gamma1(zplus,r,phi,cnz,cnr,cm)
	   hp=hamilt1(zplus,r,phi,cnz,cnr,cm)
	   bmod=b(zminus,r,phi)
	   gam=gamma1(zminus,r,phi,cnz,cnr,cm)
	   hm=hamilt1(zminus,r,phi,cnz,cnr,cm)
	   dddz=(hp-
     4  	 hm)/(2.d0*hz)
           deru(4)=-dddz

   
	   rplus=r+hr
	   rminus=r-hr
	   bmod=b(z,rplus,phi)
	   gam=gamma1(z,rplus,phi,cnz,cnr,cm)
	   hp=hamilt1(z,rplus,phi,cnz,cnr,cm)
	   bmod=b(z,rminus,phi)
	   gam=gamma1(z,rminus,phi,cnz,cnr,cm)
	   hm=hamilt1(z,rminus,phi,cnz,cnr,cm)
	   dddr=(hp-
     5  	   hm)/(2.d0*hr)
           deru(5)=-dddr


	   phipls=phi+hphi
	   phimin=phi-hphi
	   bmod=b(z,r,phipls)
	   gam=gamma1(z,r,phipls,cnz,cnr,cm)
	   hp=hamilt1(z,r,phipls,cnz,cnr,cm)
	   bmod=b(z,r,phimin)
	   gam=gamma1(z,r,phimin,cnz,cnr,cm)
	   hm=hamilt1(z,r,phimin,cnz,cnr,cm)
	   dddphi=(hp-
     6  	   hm)/(2.d0*hphi)
           deru(6)=-dddphi
    
c-----------------------------------------------------------
	   bmod=b(z,r,phi)
	   gam=gamma1(z,r,phi,cnz,cnr,cm)
          do 11 i=1,nbulk
	  vp(i)=v(i)
	  wp(i)=w(i)
 11       continue
	   frqncpl=frqncy+hfrqnc
	   df=frqncy/frqncpl
          do 12 i=1,nbulk
	  v(i)=vp(i)*df* df
	  w(i)=wp(i)*df
 12       continue
c************************************************
	   cnrplus=cnr*df
	   cnzplus=cnz*df
	   cmplus=cm*df
	   hp=hamilt1(z,r,phi,cnzplus,cnrplus,cmplus)
c----------------------------------------------------------
	   frqncmn=frqncy-hfrqnc
	   df=frqncy/frqncmn
          do 15 i=1,nbulk
	  v(i)=vp(i)*df*df
	  w(i)=wp(i)*df
 15       continue
c************************************************
	   cnrminus=cnr*df
	   cnzminus=cnz*df
	   cmminus=cm*df
	   hm=hamilt1(z,r,phi,cnzminus,cnrminus,cmminus)

c*************************************************
	   dddw=(hp-hm)/(2.0d0*hw)
c-----------------------------------------------------------
          do 14 i=1,nbulk
	  v(i)=vp(i)
	  w(i)=wp(i)
 14       continue
          p=-1.d0/dddw
c-----------------------------------------------------------
	  go to 1953
c----------------------------------------------------------------
c  end of numerical calculation of hamiltonian derivatives
c------------------------------------------------------------------
 1955	  continue
c-----------------------------------------------------------------
c analytical calculations of hamiltonian derivatives
c-----------------------------------------------------------------
           if (((id.eq.4).or.(id.eq.5)).or.(id.eq.7)) then
             write(*,*)'dddrz1: the given dispersion relation id=',id
             write(*,*)'does not work with the analytical derivatives'
             stop
           endif

	   bmod=b(z,r,phi)
	   gam=gamma1(z,r,phi,cnz,cnr,cm)
           ds=dsin(gam)
           dc=dcos(gam)
             ds2=ds*ds
             dc2=dc*dc
             
          dddphi=0.d0 ! YuP: will be found below             
             
          if (id.eq.3) then
c
c         Appleton-Hartry dispersion relation
c
          dxidz=dxdz(z,r,phi,1)
          dxidr=dxdr(z,r,phi,1)
          dxidph=dxdphi(z,r,phi,1)
          dyidz=dydz(z,r,phi,1)
          dyidr=dydr(z,r,phi,1)
          dyidph=dydphi(z,r,phi,1)

          xi=x(z,r,phi,1)
          yi=y(z,r,phi,1)

            py2=yi*yi
            py3=py2*yi
            py4=py2*py2
            px=1.d0-xi
            px2=px*px
            ds4=ds2*ds2

              det=py4*ds4+4.*py2*px2*dc2
              zer=0.d0
              if (det.lt.zer) then
               write(*,*)'det in rside less then 0 det=',det
	       stop
              end if
              sqrdet=dsqrt(det)
              p1=0.5d0/sqrdet
                ddetdx=-8.d0*py2*px*dc2
                ddetdy=4.d0*py3*ds4+8.d0*yi*px2*dc2
                pz=2.d0*px-py2*ds2+ioxm*sqrdet
                pz2=1.d0/(pz*pz)
                  dpzdx=-2.d0+ioxm*p1*ddetdx
                  dpzdy=-2.d0*yi*ds2+ioxm*p1*ddetdy
                    dddx=-2.d0*((2.d0*xi-1.d0)*pz+xi*px*dpzdx)*pz2
                    dddy=-2.d0*xi*px*dpzdy*pz2

              dddc2=-2.d0*xi*px*pz2*(py2+ioxm*(-2.d0*py4*(1.d0-dc2)+
     1                    4.d0*py2*px2)*p1)
c-----------------------------------------------------------------
                      dddz=dddx*dxidz+dddy*dyidz+dddc2*dc2dz
                      dddr=dddx*dxidr+dddy*dyidr+dddc2*dc2dr-
     1                     2*cm*cm/(r**3)
                      dddphi=dddx*dxidph+dddy*dyidph+dddc2*dc2dph

                      dddcnz=2.d0*cnz+dddc2*dc2dnz
                      dddcnr=2.d0*cnr+dddc2*dc2dnr
                      dddcm=2.d0*cm/(r*r)+dddc2*dc2dm
c--------------------------------------------------------------------
            goto 50
           end if
c   end if Appleton - Hartry
c-----------------------------------------------------------------
           if ((id.eq.1).or.(id.eq.2)) then
c  ---------------------------------------------------------------
c  cold plasma dispersion relation with electrons and ions
c  x(i=1),y(i=1) electrons component
c  x(i may be=2,nbulk),y(i may be=2,nbulk) ions components
c  ib number of component for which delib=1-yib may be equal
c  zero inside the plasma,ib  may be=from 1 to nbulk
c  dispersion relation is multiplied by delib
c

	   pc=1.d0+dc2
         call s(z,r,phi,s1,s2,s3,s4,s6,s7)
         xib=x(z,r,phi,ib)
         yib=y(z,r,phi,ib)
           delib=1.d0-yib

c  ib =1 (cyclotron resonance conditions dele=0 may be
c         realised in plasma, dispersion relation is
c         multiplied by dele )
c         ad,bd,cd calculations in dispersion relation
c       d=a*n**4+b*n**2+c
c------------------------------------------------------------------
c  if 2 begin
      if (ib.eq.1) then
       xe=xib
       ye=yib
         !peym=xe/(1.d0-ye)!YuP[07-2017] commented: not used in this case; can be 1/0
	 peyp=xe/(1.d0+ye)
	   a0e=-peyp*ds2
	   a1e=s7*ds2+s4*dc2
	   b0e=s4*peyp*pc+xe*(s6-peyp)*ds2
	   b1e=-s4*s7*pc-s3*(s6-peyp)*ds2
	   c0e=-xe*s4*(s6-peyp)
	   c1e=s4*s3*(s6-peyp)
	   dele=1.d0-ye
	     ad=(dele*a1e+a0e)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
	     bd=(dele*b1e+b0e)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
	     cd=(dele*c1e+c0e)*sign(1.d0,dele) !YuP[07-2017] corrected by sign() 
c----------------------------------------------------------------------
	     da1edc=-s7+s4
	     da0edc=peyp
	     db1edc=-s7*s4+s3*(s6-peyp)
	     db0edc=s4*peyp-xe*(s6-peyp)
	     dc1edc=0.d0
	     dc0edc=0.d0
	     dadc2=(dele*da1edc+da0edc)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
	     dbdc2=(dele*db1edc+db0edc)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
	     dcdc2=(dele*dc1edc+dc0edc)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
c-------------------------------------------------------------------
      end if
c  if 2 end
c-----------------------------------------------------------------
c  ib .gt.1 (cyclotron resonance conditions delib=0 may be
c            realised in plasma, dispersion relation is
c            multiplied by delib )
c            ad,bd,cd calculations
c	     d=a*n**4+b*n**2+c
c
c  if 3 begin
      if (ib.gt.1) then
       xe=x(z,r,phi,1)
       ye=y(z,r,phi,1)
         ppe=1.d0/(1.d0-ye)
         peym=xe/(1.d0-ye)
	 peyp=xe/(1.d0+ye)
	 peym2=peyp*ppe

         !ppb=1.d0/(1.d0-yib)
         !pbym=xib/(1.d0-yib)!YuP but this can be 1/0 ? Not used 
	 pbyp=xib/(1.d0+yib)
	 !pbym2=pbyp/(1.d0-yib) !YuP but this can be 1/0 ? Not used

	   a0b=-pbyp*ds2
	   a1b=(s1-peym2)*ds2+s4*dc2
	   b0b=s4*pbyp*pc+xib*(s3-peym)*ds2
	   b1b=-s4*(s1-peym2)*pc-(s2-peyp)*(s3-peym)*ds2
	   c0b=-xib*s4*(s3-peym)
	   c1b=s4*(s2-peyp)*(s3-peym)
	   delib=1.d0-yib
	     ad=(delib*a1b+a0b)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
	     bd=(delib*b1b+b0b)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
	     cd=(delib*c1b+c0b)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
c--------------------------------------------------------------------
	     da1bdc=-(s1-peym2)+s4
	     da0bdc=pbyp
	     db1bdc=-s4*(s1-peym2)+(s2-peyp)*(s3-peym)
	     db0bdc=s4*pbyp-xib*(s3-peym)
	     dc1bdc=0.d0
	     dc0bdc=0.d0

	     dadc2=(delib*da1bdc+da0bdc)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
	     dbdc2=(delib*db1bdc+db0bdc)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
 	     dcdc2=(delib*dc1bdc+dc0bdc)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
c-------------------------------------------- -----------0-----------
      end if
c  if 3 end
c-----------------------------------------------------------------

c  derivatives calculations
c   dadz,dadr,dadphi
c   dbdz,dbdr,dbdphi
c   dcdz,dcdr,dcdphi
c   dadw,dbdw,dcdw
c----------------------------------------------------------------
       dadz=0.d0
       dadr=0.d0
       dadphi=0.d0
       dbdz=0.d0
       dbdr=0.d0
       dbdphi=0.d0
       dcdz=0.d0
       dcdr=0.d0
       dcdphi=0.d0
       dadw=0.0d0
       dbdw=0.0d0
       dcdw=0.0d0

c  if nbulk.gt.1  (ions components are in plasma )
c
c  if 4 begin
      if (nbulk.gt.1) then
       do 10 i=2,nbulk
c
c  derivatives calculations
c   dadx,dady
c   dbdx,dbdy
c   dcdx,dcdj

         dxidz=dxdz(z,r,phi,i)
         dxidr=dxdr(z,r,phi,i)
         dxidph=dxdphi(z,r,phi,i)
         dyidz=dydz(z,r,phi,i)
         dyidr=dydr(z,r,phi,i)
         dyidph=dydphi(z,r,phi,i)
c--------------------------------------------
cc test
cc	 write(*,*)'10 in rzide1 dxidz,dxidr,dxidph',
cc     1	 dxidz,dxidr,dxidph
cc	 write(*,*)'in rzide1 dyidz,dyidr,dyidph',
cc     1	 dyidz,dyidr,dyidph
cc test
c--------------------------------------------

         if ( i.eq.ib) goto 20 ! resonant ions
c
c  i.ne.ib ! Non-resonant ions
c
	 xi=x(z,r,phi,i)
         yi=y(z,r,phi,i)
	   pyim=1.d0/(1.d0-yi)
	   pyip=1.d0/(1.d0+yi)
	   pyim2=pyim*pyim

c if 5 begin
	     if (ib.gt.1) then
	       ds1dxi=pyim*pyip
	       ds2dxi=pyim
	       ds1dyi=2.d0*xi*yi*ds1dxi*ds1dxi
	       ds2dyi=xi*pyim2
	     end if
c if 5 end
             ds3dxi=pyip
	     ds4dxi=1.d0
             ds3dyi=-xi*pyip*pyip
             ds4dyi=0.d0
c if 6 begin
	     if (ib.eq.1) then
	       ds6dxi=pyim
	       ds7dxi=pyim*pyip
	       ds6dyi=xi*pyim2
	       ds7dyi=2.d0*xi*yi*ds7dxi*ds7dxi
	     end if
c if 6 end

c if 7 begin
	     if (ib.eq.1) then

               da1exi=-ds7dxi*ds2-ds4dxi*dc2
	       da0exi=0.d0
	       db1exi=(ds4dxi*s7+ds7dxi*s4)*pc+
     1		      (ds3dxi*(s6-xe/(1.d0+ye))+ds6dxi*s3)*ds2
               db0exi=-ds4dxi*xe/(1.d0+ye)*pc-
     1                 xe*ds6dxi*ds2
               dc1exi=-ds4dxi*s3*(s6-xe/(1.d0+ye))-
     1                 s4*ds3dxi*(s6-xe/(1.d0+ye))-s4*s3*ds6dxi
               dc0exi=xe*ds4dxi*(s6-xe/(1.d0+ye))+xe*s4*ds6dxi

	       da1eyi=-ds7dyi*ds2-ds4dyi*dc2
	       da0eyi=0.d0
	       db1eyi=(ds4dyi*s7+ds7dyi*s4)*pc+
     1                (ds3dyi*(s6-xe/(1.d0+ye))+ds6dyi*s3)*ds2
               db0eyi=-ds4dyi*xe/(1.d0+ye)*pc-xe*ds6dyi*ds2
	       dc1eyi=-ds4dyi*s3*(s6-xe/(1.d0+ye))-
     1    	       ds3dyi*s4*(s6-xe/(1.d0+ye))-ds6dyi*s4*s3
               dc0eyi=xe*ds4dyi*(s6-xe/(1.d0+ye))+ds6dyi*xe*s4

                 dadxi=(dele*da1exi+da0exi)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dbdxi=(dele*db1exi+db0exi)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dcdxi=(dele*dc1exi+dc0exi)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dadyi=(dele*da1eyi+da0eyi)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dbdyi=(dele*db1eyi+db0eyi)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dcdyi=(dele*dc1eyi+dc0eyi)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()

	       goto 30
	     end if
c if 7 end

c if 8 begin
	     if (ib.gt.1) then

               da1bxi=-ds1dxi*ds2-ds4dxi*dc2
	       da0bxi=0.d0
	       db1bxi=(ds4dxi*(s1-xe/(1.d0-ye*ye))+ds1dxi*s4)*pc+
     1		  (ds2dxi*(s3-xe/(1.d0-ye))+ds3dxi*(s2-xe/(1.d0+ye)))*
     2                 ds2
               db0bxi=-ds4dxi*xib/(1.d0+yib)*pc-
     1                 xib*ds3dxi*ds2
               dc1bxi=-ds4dxi*(s2-xe/(1.d0+ye))*(s3-xe/(1.d0-ye))-
     1                 ds2dxi*s4*(s3-xe/(1.d0-ye))-
     2                 ds3dxi*s4*(s2-xe/(1.d0+ye))
               dc0bxi=xib*ds4dxi*(s3-xe/(1.d0-ye))+xib*s4*ds3dxi

	       da1byi=-ds1dyi*ds2-ds4dyi*dc2
	       da0byi=0.d0
	       db1byi=(ds4dyi*(s1-xe/(1.d0-ye*ye))+s4*ds1dyi)*pc+
     1          (ds2dyi*(s3-xe/(1.d0-ye))+ds3dyi*(s2-xe/(1.d0+ye)))*
     2                 ds2
               db0byi=-ds4dyi*xib/(1.d0+yib)*pc-xib*ds3dyi*ds2
	       dc1byi=s4*(-ds2dyi*(s3-xe/(1.d0-ye))-
     1    	       ds3dyi*(s2-xe/(1.d0+ye)))
               dc0byi=xib*s4*ds3dyi

                 dadxi=(delib*da1bxi+da0bxi)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dbdxi=(delib*db1bxi+db0bxi)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dcdxi=(delib*dc1bxi+dc0bxi)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dadyi=(delib*da1byi+da0byi)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dbdyi=(delib*db1byi+db0byi)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dcdyi=(delib*dc1byi+dc0byi)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()

	       goto 30
	     end if
c if 8 end
 20     continue
c
c         i=ib, i.ne.1
c

c if 15 begin
	     if (ib.gt.1) then
	       ds1dxb=0.d0
	       ds2dxb=0.d0
	       ds1dyb=0.d0
	       ds2dyb=0.d0
	     end if
c if 15 end
             ds3dxb=1.d0/(1.d0+yib)
	     ds4dxb=1.d0
             ds3dyb=-xib/(1.d0+yib)**2
             ds4dyb=0.d0
c if 16 begin
	     if (ib.eq.1) then
	       ds6dxb=0.d0
	       ds7dxb=0.d0
	       ds6dyb=0.d0
	       ds7dyb=0.d0
	     end if
c if 16 end

c if 9 begin
	     if (ib.gt.1) then

               da1bxb=-ds1dxb*ds2-ds4dxb*dc2
               da0bxb=-1.d0/(1.d0+yib)*ds2
               db1bxb=(ds4dxb*(s1-xe/(1.d0-ye*ye))+ds1dxb*s4)*pc+
     1            (ds2dxb*(s3-xe/(1.d0-ye))+ds3dxb*(s2-xe/(1.d0+ye)))*
     2                 ds2
               db0bxb=(-ds4dxb*xib/(1.d0+yib)+s4/(1.d0+yib))*pc+
     1                ((s3-xe/(1.d0-ye))-xib*ds3dxb)*ds2
               dc1bxb=-ds4dxb*(s2-xe/(1.d0+ye))*(s3-xe/(1.d0-ye))-
     1                 ds2dxb*s4*(s3-xe/(1.d0-ye))-
     2                 ds3dxb*s4*(s2-xe/(1.d0+ye))
            dc0bxb=-s4*(s3-xe/(1.d0-ye))+xib*ds4dxb*(s3-xe/(1.d0-ye))+
     +                xib*s4*ds3dxb

               da1byb=0.d0
               da0byb=xib/(1.d0+yib)**2*ds2
               db1byb=(ds4dyb*(s1-xe/(1.d0-ye*ye))+s4*ds1dyb)*pc+
     1          (ds2dyb*(s3-xe/(1.d0-ye))+ds3dyb*(s2-xe/(1.d0+ye)))*
     2                 ds2
               db0byb=(-ds4dyb*xib/(1.d0+yib)-s4*xib/(1.d0+yib)**2)*pc-
     -                xib*ds3dyb*ds2
               dc1byb=s4*(-ds2dyb*(s3-xe/(1.d0-ye))-
     1                 ds3dyb*(s2-xe/(1.d0+ye)))
               dc0byb=xib*s4*ds3dyb

                 dadxi=(delib*da1bxb+da0bxb)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dbdxi=(delib*db1bxb+db0bxb)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dcdxi=(delib*dc1bxb+dc0bxb)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dadyi=(delib*da1byb+da0byb-a1b)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dbdyi=(delib*db1byb+db0byb-b1b)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dcdyi=(delib*dc1byb+dc0byb-c1b)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()

	     end if
c if 9 end

 30       continue

          dadz=dadz+dadxi*dxidz+dadyi*dyidz
          dadr=dadr+dadxi*dxidr+dadyi*dyidr
          dadphi=dadphi+dadxi*dxidph+dadyi*dyidph

          dbdz=dbdz+dbdxi*dxidz+dbdyi*dyidz
          dbdr=dbdr+dbdxi*dxidr+dbdyi*dyidr
          dbdphi=dbdphi+dbdxi*dxidph+dbdyi*dyidph

          dcdz=dcdz+dcdxi*dxidz+dcdyi*dyidz
          dcdr=dcdr+dcdxi*dxidr+dcdyi*dyidr
          dcdphi=dcdphi+dcdxi*dxidph+dcdyi*dyidph

c new-------------------------------------
	 xi=x(z,r,phi,i)
         yi=y(z,r,phi,i)
c*************************************c
	  dadw=dadw-2*dadxi*xi-dadyi*yi
	  dbdw=dbdw-2*dbdxi*xi-dbdyi*yi
	  dcdw=dcdw-2*dcdxi*xi-dcdyi*yi
c new--------------------------------------

 10       continue
         end if
c if 4 nbulk > 1 end


c
c     i = 1
c
             ds1dxe=0.d0
             ds2dxe=0.d0
             ds3dxe=0.d0
             ds4dxe=1.d0
             ds6dxe=0.d0
             ds7dxe=0.d0

             ds1dye=0.d0
             ds2dye=0.d0
             ds3dye=0.d0
             ds4dye=0.d0
             ds6dye=0.d0
             ds7dye=0.d0

c if 10 begin
c
c   ib > 1, i = 1
c
	     if (ib.gt.1) then
               da1bxe=-1.d0/(1.d0-ye*ye)*ds2-ds4dxe*dc2
               da0bxe=0.d0
               db1bxe=(ds4dxe*(s1-xe/(1.d0-ye*ye))+s4/(1.d0-ye*ye))*pc+
     1                (1.d0/(1.d0+ye)*(s3-xe/(1.d0-ye))+
     +                 1.d0/(1.d0-ye)*(s2-xe/(1.d0+ye)))*ds2
               db0bxe=-ds4dxe*xib/(1.d0+yib)*pc-xib/(1.d0-ye)*ds2
               dc1bxe=s4*(-1.d0/(1.d0+ye)*(s3-xe/(1.d0-ye))-
     -                     1.d0/(1.d0-ye)*(s2-xe/(1.d0+ye)))-
     -                 ds4dxe*((s2-xe/(1.d0+ye))*(s3-xe/(1.d0-ye)))
               dc0bxe=xib*s4/(1.d0-ye)+xib*ds4dxe*(s3-xe/(1.d0-ye))

               da1bye=-2.d0*xe*ye/(1.d0-ye*ye)**2*ds2-ds4dye*dc2
               da0bye=0.d0
               db1bye=s4*2.*xe*ye/(1.d0-ye*ye)**2*pc+
     +                (-xe/(1.d0+ye)**2*(s3-xe/(1.d0-ye))+
     +                  xe/(1.d0-ye)**2*(s2-xe/(1.d0+ye)))*ds2
               db0bye=-xib*xe/(1.d0-ye)**2*ds2
               dc1bye=s4*(xe/(1.d0+ye)**2*(s3-xe/(1.d0-ye))-
     -                    xe/(1.d0-ye)**2*(s2-xe/(1.d0+ye)))
               dc0bye=xib*s4*xe/(1.d0-ye)**2

                 dadxi=(delib*da1bxe+da0bxe)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dbdxi=(delib*db1bxe+db0bxe)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dcdxi=(delib*dc1bxe+dc0bxe)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dadyi=(delib*da1bye+da0bye)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dbdyi=(delib*db1bye+db0bye)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()
                 dcdyi=(delib*dc1bye+dc0bye)*sign(1.d0,delib) !YuP[07-2017] corrected by sign()

               goto 40
	     end if
c if 10 end

c if 11 begin
c
c   ib = 1, i = 1
c
             if (ib.eq.1) then

               da1exe=-ds4dxe*dc2
               da0exe=-1.d0/(1.d0+ye)*ds2
               db1exe=ds4dxe*s7*pc+s3/(1.d0+ye)*ds2
               db0exe=(-ds4dxe*xe/(1.d0+ye)+s4/(1.d0+ye))*pc+
     +                 (s6-2.d0*xe/(1.d0+ye))*ds2
               dc1exe=-ds4dxe*s3*(s6-xe/(1.d0+ye))-
     -                s4*s3/(1.d0+ye)
               dc0exe=-s4*(s6-2.d0*xe/(1.d0+ye))+xe*(s6-xe/(1.d0+ye))

               da1eye=0.d0
               da0eye=xe/(1.d0+ye)**2*ds2
               db1eye=-s3*xe/(1.d0+ye)**2*ds2
               db0eye=-s4*xe/(1.d0+ye)**2*pc+xe*xe/(1.d0+ye)**2*ds2
               dc1eye=-ds4dye*s3*(s6-xe/(1.d0+ye))+s4*s3*xe/(1.d0+ye)**2
               dc0eye=-xe*xe*s4/(1.d0+ye)**2

                 dadxi=(dele*da1exe+da0exe)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dbdxi=(dele*db1exe+db0exe)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dcdxi=(dele*dc1exe+dc0exe)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dadyi=(dele*da1eye+da0eye-a1e)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dbdyi=(dele*db1eye+db0eye-b1e)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()
                 dcdyi=(dele*dc1eye+dc0eye-c1e)*sign(1.d0,dele) !YuP[07-2017] corrected by sign()


	     end if
c if 11 end

 40     continue
         dxidz=dxdz(z,r,phi,1)
         dxidr=dxdr(z,r,phi,1)
         dxidph=dxdphi(z,r,phi,1)
         dyidz=dydz(z,r,phi,1)
         dyidr=dydr(z,r,phi,1)
         dyidph=dydphi(z,r,phi,1)
c------------------------------------------------------------------
        dadz=dadz+dadxi*dxidz+dadyi*dyidz+dadc2*dc2dz
        dadr=dadr+dadxi*dxidr+dadyi*dyidr+dadc2*dc2dr
        dadphi=dadphi+dadxi*dxidph+dadyi*dyidph+dadc2*dc2dph

        dbdz=dbdz+dbdxi*dxidz+dbdyi*dyidz+dbdc2*dc2dz
        dbdr=dbdr+dbdxi*dxidr+dbdyi*dyidr+dbdc2*dc2dr
        dbdphi=dbdphi+dbdxi*dxidph+dbdyi*dyidph+dbdc2*dc2dph

        dcdz=dcdz+dcdxi*dxidz+dcdyi*dyidz+dcdc2*dc2dz
        dcdr=dcdr+dcdxi*dxidr+dcdyi*dyidr+dcdc2*dc2dr
        dcdphi=dcdphi+dcdxi*dxidph+dcdyi*dyidph+dcdc2*dc2dph
c--------------------------------------------
c new-------------------------------------
	 xi=x(z,r,phi,1)
         yi=y(z,r,phi,1)
	  dadw=dadw-2*dadxi*xi-dadyi*yi
	  dbdw=dbdw-2*dbdxi*xi-dbdyi*yi
	  dcdw=dcdw-2*dcdxi*xi-dcdyi*yi

c new--------------------------------------

        dcn=cn(r,cnz,cnr,cm)
        dcn2=dcn*dcn
        dcn4=dcn2*dcn2

	dadnz=dadc2*dc2dnz
	dbdnz=dbdc2*dc2dnz
	dcdnz=dcdc2*dc2dnz

	dadnr=dadc2*dc2dnr
	dbdnr=dbdc2*dc2dnr
	dcdnr=dcdc2*dc2dnr

	dadm=dadc2*dc2dm
	dbdm=dbdc2*dc2dm
	dcdm=dcdc2*dc2dm
c------------------------------------------------------------------

        if (id.eq.1) then

c
c  dispersion relation a*n**4+b*n**2+c=0
c

          dddz=dcn4*dadz+dcn2*dbdz+dcdz
          dddr=dcn4*dadr+dcn2*dbdr+dcdr-
     1	       cm*cm*(4.d0*dcn2*ad+2.d0*bd)/(r**3)
          dddphi=dcn4*dadphi+dcn2*dbdphi+dcdphi

          dddn=4.d0*ad*dcn2+2.d0*bd
          dddcnz=dddn*cnz+dcn4*dadnz+dcn2*dbdnz+dcdnz
          dddcnr=dddn*cnr+dcn4*dadnr+dcn2*dbdnr+dcdnr
          dddcm=dddn*cm/(r*r)+dcn4*dadm+dcn2*dbdm+dcdm

        end if

        if (id.eq.2) then

c
c  dispersion relation n**2-(-b+iom*sqrt(b*2-4*a*c))/(2*a)=0
c

         det=dsqrt(bd*bd-4.d0*ad*cd)
         p4=0.5d0/(ad*ad)
         p5=1.d0/det
         p6=-bd+ioxm*det
c------------------------------------------------------------------
         p7=ioxm*(bd*dbdz-2.d0*dadz*cd-2.d0*ad*dcdz)*p5
         dddz=-((-dbdz+p7)*ad-dadz*p6)*p4
         p7=ioxm*(bd*dbdr-2.d0*dadr*cd-2.d0*ad*dcdr)*p5
         dddr=-((-dbdr+p7)*ad-dadr*p6)*p4-
     1		2.d0*cm**2/(r**3)
         p7=ioxm*(bd*dbdphi-2.d0*dadphi*cd-2.d0*ad*dcdphi)*p5
         dddphi=-((-dbdphi+p7)*ad-dadphi*p6)*p4

         p7=ioxm*(bd*dbdnz-2.d0*dadnz*cd-2.d0*ad*dcdnz)*p5
          dddcnz=2.d0*cnz-
     1      	 ((-dbdnz+p7)*ad-dadnz*p6)*p4
         p7=ioxm*(bd*dbdnr-2.d0*dadnr*cd-2.d0*ad*dcdnr)*p5
          dddcnr=2.d0*cnr-
     1      	 ((-dbdnr+p7)*ad-dadnr*p6)*p4
         p7=ioxm*(bd*dbdm-2.d0*dadm*cd-2.d0*ad*dcdm)*p5
          dddcm=2.d0*cm/(r*r)-
     1      	 ((-dbdm+p7)*ad-dadm*p6)*p4
c--------------------------------------------------------------------
        end if
c--------------------------------------------------------------------
        goto 50
        end if
c   end of cold plasma  dispersion
c-----------------------------------------------------------------
        if (id.eq.6)then
c-----------------------------------------------------------------
c       hot non-relativistic plasma dispersion from Forest code
c-----------------------------------------------------------------
            dddph=dddphi ! YuP[07-2017] added: dddph was not defined.
            !(may not be important - over-written in hotdervs)
c            call hotderiv(u,frqncy,dddcnz,dddcnr,dddcm,
c     .                  dddz,dddr,dddph,dddw)
            call hotdervs(nbulk,u,frqncy,dddcnz,dddcnr,dddcm,
     .                  dddz,dddr,dddph,dddw)
          goto 50
        end if
c       end of hot non-relativistic plasma dispersion from Forest code
c-----------------------------------------------------------------
        if (id.eq.8)then
c-----------------------------------------------------------------
c           Ono dispersion for fast waves 
c-----------------------------------------------------------------
          wf=frqncy
            dddph=dddphi ! YuP[07-2017] added: dddph was not defined
            !(may not be important - over-written in ono_dervs)
          call  ono_dervs(nbulk,u,wf,dddcnz,dddcnr,dddcm,
     .                  dddz,dddr,dddph,dddw)

          goto 50
        end if
c------------------------------------------------------------
       if (id.eq.16)then
c-----------------------------------------------------------------
c           Bonoli dispersion for LH waves 
c-----------------------------------------------------------------
          wf=frqncy
c          write(*,*)'in dddrz1.f can not work for i=16 idif=1'
c          stop 'in dddrz1.f'
        
           call bonoli_dervs(u,wf,
     &     dddz,dddr,dddphi,
     &     dddcnz,dddcnr,dddcm,dddw)   
          goto 50
        end if


c       end of hot non-relativistic plasma dispersion from Forest code
c-----------------------------------------------------------------
 50     continue
           deru(1)=dddcnz
           deru(2)=dddcnr
           deru(3)=dddcm
           deru(4)=-dddz
           deru(5)=-dddr
           deru(6)=-dddphi


ctest
c           sum=0.d0
c           do i=1,6
c              sum=sum+deru(i)**2
c           enddo
c           sum=dsqrt(sum)
c           do i=1,6
c              deru(i)=deru(i)/sum
c           enddo
cendtest
            

 1953   continue

	  return
          end


