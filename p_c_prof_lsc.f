      SUBROUTINE p_c_prof_lsc(is,rhobegin,rhoend,cnpar,cnper,
     &cefldx,cefldy,cefldz)
c----------------------------------------------------------------
c     this subroutine calculates absorted power profiles
c     calculates RF current  profiles	 (from deposited power)
c-----------------------------------------------------------------
c      IMPLICIT double precision (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      INCLUDE 'three.i'
      INCLUDE 'onetwo.i'
      INCLUDE 'write.i'
      INCLUDE 'gr.i'
      INCLUDE 'rho.i'
     
c-----input
      real*8 rhobegin,rhoend,cnpar,cnper
      integer is
      complex*16 cefldx,cefldy,cefldz !polarization

c-----locals
      real*8 delpow_s,ppow_s
      dimension delpow_s(nbulka) !Added for indiv ion contrib[BH041009]
      dimension ppow_s(nbulka)  !Added for indiv ion contrib[BH041009]

      real*8 z_r,r_r,ymar,xmar,z_effr,tempr,denr,cnparr,yer,effic_r

      real*8 r_m,z_m,temp,ye,u1,den,z_eff,hro,rholeft,rhoright,
     &rhomax,rhomin,eff_rho_max,eff_rho_min,hrho,delpower,delpow_i,
     &delpow_e,delpow_cl,del,r0,z0,r0m,z0m,psiloc,zfacgeom,aspct,
     &delcurr,rho0,poloidlen,rho0_pol,delrho,
     &ppow,pcur,ppow_e,ppow_i,ppow_cl,
     &r_bmin_right,r_bmin_left,z_bmin_right,z_bmin_left,theta,psi

      integer i,kk,jbinmin,jbinmax,j,k
cSAP070831
      real*8 psi_loc,cos_theta_pol,sin_theta_pol,rho_loc,theta_pol
      integer n_theta_pol,ir
c-----external
      real*8  tempe,y,u_res,dense,zeff,psi_rho,b,efficien,rhov,rhos,
     &qsafety_psi,bmin_psi,bmax_psi,dvol_dpsi,rho_lrho,
     &b_average

      pi=4.d0*datan(1.d0)

c     wr and wz are in (cm) 
      r_m=wr(is)*0.01d0/r0x ! normalization
      z_m=wz(is)*0.01d0/r0x ! normalization

c       write(*,*)'p_c_prof NR,is,ieffic',NR,is,ieffic

c  radius rho is in common/one/,it was calculated in subroutine: b
c  rho is used in functions:tempe,dense,z_eff
c  The magnetic field bmode is in common/one/ .bmode is used
c  in function y(z,r,phi)


c begin if is=1      

      if(is.eq.1) then

         write(*,*)'prep3d.f p_c_prof is=1 rho=',rho
cSAP080731
        if (rho.ge.1.d0) then
c---------------------------------------------------------------
c          zero CD efficiency outside the plasma
c----------------------------------------------------------------
           eff(is)=0.d0
           goto 100
        endif
  


c        calculation of efficiency on the first step
         temp=tempe(z_m,r_m,wphi(is),1)
         ye=y(z_m,r_m,wphi(is),1)
         u1=u_res(jwave,cnpar,temp,ye)
         den=dense(z_m,r_m,wphi(is),1)
         z_eff=zeff(z_m,r_m,wphi(is))
 

         if (ieffic.eq.2) then

c-------------------------------------------------------------------
           write(*,*)'lsc_approach.f  in p_c_prof_lsc before'//
     &               'asimptotic_power_CD_density_LSC'
 
           call asimptotic_power_CD_density_LSC(
     &     z_m,r_m,0.d0,cnpar,cnper)
c     &     CD_dens,power_dens,CD_efficiency)
           write(*,*)'lsc_approach.f  in p_c_prof_lsc after'//
     &               'asimptotic_power_CD_density_LSC'

           write(*,*)'prep3d z_eff',z_eff
           call efKarney(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
     +                 z_eff,temp,den,jwave,
     1                 cnpar,eff(is))
           write(*,*)'prep3d Karney is,cnpar,eff(is)',is,cnpar,eff(is)
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coincided exactly with the efficiency obtained 
c          calculated bu subroutine call efKarney
c
c-------------------------------------------------------
c          call efKarney_Bonoli(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
c     +                 z_eff,temp,den,jwave,
c     1                 cnpar,eff(is))
c          write(*,*)'prep3d Karney_Bonoli is,cnpar,eff(is)',
c     &                is,cnpar,eff(is)
c------------------------------------------------------------
	 endif   ! ieffic.eq.2
          
c--------------------------------------------------------------------
c        toroidal and poloidal current drives efficiencies for is=1
 100     continue
         bmod=b(z_m,r_m,0.d0)        
c--------------------------------------------------------------------
         allpower=0.0d0
         allpw_e=0.0d0
         allpw_i=0.0d0
         allpw_cl=0.0d0

         allcur=0.0d0
        
c------- initialization arrays
c        for power and current
         do i=1,NR
	    power(i)=0.0d0
	    current(i)=0.0d0           
	    power_e(i)=0.0d0
	    power_i(i)=0.0d0
	    power_cl(i)=0.0d0
            cur_den_parallel(i)=0.d0 
	 enddo

c--------binvol(NR) calculations
         theta=0.d0 
         hro=1.d0/dble(NR-1)
         do i=1,NR-1 !for onetwo.i in [cm**3]
            rholeft=hro*(i-1)
            rhoright=rholeft+hro
            binvol(i)=voltot*(rhov(rhoright)**2-rhov(rholeft)**2)*1.d6

c            write(*,*)'NR,i,binvol(i)',NR,i,binvol(i)

           binarea(i)=areatot*(rhos(rhoright)**2-rhos(rholeft)**2)*1.d4

            psi=psi_rho(rholeft)
            call zr_psith(psi,theta,z_bmin_left,r_bmin_left)
            psi=psi_rho(rhoright)
            call zr_psith(psi,theta,z_bmin_left,r_bmin_right)

            binarea_pol(i)=pi*(r_bmin_right-r_bmin_left)*
     &                        (r_bmin_right+r_bmin_left)*1.d4

         enddo
 
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               do i=1,NR
                  power_s(i,kk)=0.0d0
               enddo
            enddo
         endif

	 goto 30
      endif
c end if is=1

c      write(*,*)'rhobegin,rhoend',rhobegin,rhoend

      if(rhoend.gt.rhobegin) then
         rhomax=rhoend
         rhomin=rhobegin
      else
	 rhomin=rhoend
	 rhomax=rhobegin
      endif
cSAP080831
      if(rhomax.gt.1.d0) rhomax=1.d0
      if(rhomin.gt.1.d0) rhomin=1.d0-1.d-10    

      hrho=1.0d0/(NR-1)
      jbinmin=1
      jbinmax=NR-1

c      write(*,*)'NR,hrho,rhomin,rhomax',NR,hrho,rhomin,rhomax

      do j=1,NR-1
         if(rhomin.lt.(hrho*j)) then
           jbinmin=j
	    goto 10
	 endif
      enddo
10    continue
      do j=jbinmin,NR-1
         if(rhomax.lt.(hrho*j)) then
            jbinmax=j
	    goto 20
	 endif
      enddo
20    continue

      if (jbinmin.eq.NR) jbinmin=NR-1

c      write(*,*)'jbinmax&min',jbinmax,jbinmin
     
c-----------------------------------------------------------
c     here delpower and allpower are in (erg/sec)
c------------------------------------------------------------
      delpower=delpwr(is-1)-delpwr(is)
cSmirnov970105 beg
cBH001017      delpow_e=0.5d0*(ws(is)-ws(is-1))*
cBH001017     1(delpwr(is-1)*sdpwr(is-1)+delpwr(is)*sdpwr(is))
cBH001017      delpow_i=0.5d0*(ws(is)-ws(is-1))*
cBH001017     1(delpwr(is-1)*salphal(is-1)+delpwr(is)*salphal(is))
cHarvey991017 beg
      delpow_i=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*sdpwr(is-1)+delpwr(is)*sdpwr(is))
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            delpow_s(kk)=0.5d0*(ws(is)-ws(is-1))*
     1        (delpwr(is-1)*salphas(is-1,kk)+delpwr(is)*salphas(is,kk))
         enddo
      endif
      
      delpow_e=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*salphal(is-1)+delpwr(is)*salphal(is))

cSAP100223 !for test only electron absorption
      delpow_e=delpower 
      write(*,*)'p_c_prof_lsc delpow_e', delpow_e   
cHarvey991017 end
      delpow_cl=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*salphac(is-1)+delpwr(is)*salphac(is))

      del=delpow_i+delpow_e+delpow_cl
cHarvey970111 beg
cSmirnov050215      if (del.ne.0.d0) then
cSmirnov050215      Change 0.d0 to a small number
      if (del.gt.1.d-100) then
         delpow_e=delpow_e*delpower/del
         delpow_i=delpow_i*delpower/del
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               delpow_s(kk)=delpow_s(kk)*delpower/del
            enddo
         endif
         delpow_cl=delpow_cl*delpower/del

      endif
cSmirnov970105 end
      allpower=allpower+delpower
cSAP091021
      write(*,*)'delpower,delpow_i,delpow_e',delpower,delpow_i,delpow_e
      write(*,*)'delpow_s(k),k=2,nbulk',(delpow_s(k),k=2,nbulk)
c      write(*,*)'allpower,delpower',allpower,delpower
cSmirnov970106 beg
      allpw_e=allpw_e+delpow_e
      allpw_i=allpw_i+delpow_i
      allpw_cl=allpw_cl+delpow_cl
      
cSmirnov970106 end
c     write(*,*)' allpower,allpw_e,allpw_i,allpw_cl',
c    1 allpower,allpw_e,allpw_i,allpw_cl
c     write(*,*)'delpower,delpow_e,delpow_cl',
c    1delpower,delpow_e,delpow_cl
c-----------------------------------------------------------
c     r0 (cm)
      r0=0.5d0*(wr(is)+wr(is-1))

cSAP080831
c      write(*,*)'p_c_prof is>1  rho',rho
      if (rho.ge.1.d0) then
c---------------------------------------------------------------
c        zero CD efficiency outside the plasma
c----------------------------------------------------------------
         eff(is)=0.d0
         goto 110
      endif
c---- calculation of the efficiency 
      temp=tempe(z_m,r_m,wphi(is),1)
      ye=y(z_m,r_m,wphi(is),1)
      u1=u_res(jwave,cnpar,temp,ye)
      den=dense(z_m,r_m,wphi(is),1)
      z_eff=zeff(z_m,r_m,wphi(is))
   
 
      if (ieffic.eq.2) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using Ehst-Karney formula
c--------------------------------------------------------------------
           write(*,*)'prep3d.f  in p_c_prof before'//
     &               'asimptotic_power_CD_density_LSC'
 
           call asimptotic_power_CD_density_LSC(
     &     z_m,r_m,0.d0,cnpar,cnper)
c     &     CD_dens,power_dens,CD_efficiency)
           write(*,*)'prep3d.f  in p_c_prof after'//
     &               'asimptotic_power_CD_density_LSC'

      call efKarney(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
     & z_eff,temp,den,jwave,cnpar,
     &              eff(is))
      write(*,*)'prep3d Karney is,cnpar,eff(is)',is,cnpar,eff(is)
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coinsided exectly wiht the efficiency obtained 
c          calculated bu subroutine call efKarney
c
c-------------------------------------------------------
c      call efKarney_Bonoli(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
c     & z_eff,temp,den,jwave,cnpar,
c     &              eff(is))
c      write(*,*)'prep3d Karney_Bonoli is,cnpar,eff(is)',is,cnpar,eff(is)

      endif  !ieffic.eq.2

 110  continue

c      write(*,*)'p_c_prof after 110' 

c--------------------------------------------------------------------
c     delpower(erg/sec),delcurr(Ampere),r0(cm)
c-----------------------------------------------------------
cSmirnov970105 beg

c-----calculate parallel CD using the geometric factor 1/(2pi*r0) 
c      delcurr=delpow_e*0.5d0*(eff(is-1)+eff(is))/(2*pi*r0)
c      write(*,*)'1/(2*pi*r0)',1.d0/(2.d0*pi*r0)

c------geometric factor 1/(r*pi*R_mag_axis)
c      zfacgeom = 1.0d0 / (2.0d0 * pi * xma)/100.d0
c      write(*,*)'1/(2*pi*xma)/100.d0',1.d0/(2.d0*pi*xma)/100.d0

c-----calculate toroidal CD
      z0=0.5d0*(wz(is)+wz(is-1))
      r0m=r0*1.d-2
      z0m=z0*1.d-2      
      bmod=b(z0m,r0m,0.d0)
c      write(*,*)'z0m,r0m,bmod,rho',z0m,r0m,bmod,rho
      psiloc=psi_rho(rho)
c-----geometric factor ~ 1/b_averaged 
      zfacgeom = -2.d0 * pi * qsafety_psi(psiloc)/
     &     (dvol_dpsi(psiloc)*dpsimax*b_average(psiloc))/100.d0

c      write(*,*)'qsafety_psi(psiloc),b_average(psiloc)',
c     &           qsafety_psi(psiloc),b_average(psiloc)

c      write(*,*)'dvol_dpsi(psiloc),dpsimax',dvol_dpsi(psiloc),dpsimax
c      write(*,*)'zfacgeom',zfacgeom

c-----geometric factor ~ 1/b_ min/(1+aspct) old from Toray
c      zfacgeom = -2.d0 * pi * qsafety_psi(psiloc)/
c     &     (dvol_dpsi(psiloc)*dpsimax*bmin_psi(psiloc))/100.d0
c      aspct=(bmax_psi(psiloc)-bmin_psi(psiloc))/
c     &(bmax_psi(psiloc)+bmin_psi(psiloc))
c      zfacgeom = zfacgeom / (1.d0 + aspct)
c      write(*,*)'old zfacgeom',zfacgeom

      delcurr=delpow_e*0.5d0*(eff(is-1)+eff(is))*zfacgeom !toroidal current
                                                          !created by delpow_e
c      write(*,*)'prep3d zfacgeom', zfacgeom 
c      write(*,*)' prep3d p_c_prof delpow_e,eff(is-1),eff(is),delcurr',
c     +delpow_e,eff(is-1),eff(is),delcurr
cSAP080902 to plot delta power and delta current along the ray
      delpow_e_ar(is)=delpow_e
      delcur_par_ar(is)=delpow_e*0.5d0*(eff(is-1)+eff(is))
      
c-----toroidal and poloidal CD from old genray version      
      rho0=0.5*(spsi(is)+spsi(is-1)) ! the small radius
      rho0_pol=rho_lrho(rho0)
      poloidlen=rho0_pol*totlength*100.d0     ! poloidal length cm

cSmirnov970105 end
      allcur=allcur+delcurr                     !total toroidal current
      
c      write(*,*)'allcur',allcur
c      write(*,*)'delcurr',delcurr

c      write(*,*)'p_c_prof rhoend,rhobegin',rhoend,rhobegin

      if(rhoend.gt.rhobegin) then
          eff_rho_max=eff(is)
          eff_rho_min=eff(is-1)
      else
          eff_rho_max=eff(is-1)
          eff_rho_min=eff(is)
      endif

      write(*,*)'jbinmin,jbinmax',jbinmin,jbinmax     
       
      if (jbinmin.eq.jbinmax) then

c         write(*,*)'jbinmon=jbinmax,delpower,power(jbinmin)',
c     .              jbinmin,delpower,power(jbinmin)

         power(jbinmin)=power(jbinmin)+delpower

cSAP080731
        write(*,*)'power(jbinmin)',power(jbinmin)
 
         cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
     &   delpow_e*0.5d0*(eff_rho_min+eff_rho_max)/binvol(jbinmin)

cSAP080731
        write(*,*)'delpow_e,eff_rho_min,eff_rho_max,binvol(jbinmin)',
     &             delpow_e,eff_rho_min,eff_rho_max,binvol(jbinmin)

        write(*,*)'cur_den_parallel(jbinmin)',cur_den_parallel(jbinmin)

         current(jbinmin)=current(jbinmin)+delcurr   !toroidal current from bin

c         write(*,*)'cur_den_parallel(jbinmin)*binarea(jbinmin)',
c     &   cur_den_parallel(jbinmin)*binarea(jbinmin)

c         write(*,*)'current(jbinmin)',current(jbinmin)
c         write(*,*)'zfacgeom',zfacgeom
c         write(*,*)'0.5d0*(eff(is-1)+eff(is))',
c     &   0.5d0*(eff(is-1)+eff(is))
c         write(*,*)'is,eff(is-1),eff(is)',is,eff(is-1),eff(is)
c         write(*,*)'binarea(jbinmin)/binvol(jbinmin)',
c     &   binarea(jbinmin)/binvol(jbinmin)
c         write(*,*)'0.5d0*(eff_rho_min+eff_rho_max)',
c     &   0.5d0*(eff_rho_min+eff_rho_max)
c         write(*,*)'is,eff(is-1),eff(is)',is,eff(is-1),eff(is)
c         write(*,*)'eff_rho_min,eff_rho_max',eff_rho_min,eff_rho_max
      
         power_e(jbinmin)=power_e(jbinmin)+delpow_e
         power_i(jbinmin)=power_i(jbinmin)+delpow_i
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               power_s(jbinmin,kk)=power_s(jbinmin,kk)+delpow_s(kk)
            enddo
         endif
       
         power_cl(jbinmin)=power_cl(jbinmin)+delpow_cl     

         goto 30
      endif

c      write(*,*)'prep3d  after goto30 delpower rhomax,rhomin',
c     . delpower,rhomax,rhomin    
       
c      delrho=dmax1((rhomax-rhomin),hrho) 
cSm040426
      delrho=rhomax-rhomin

      ppow=delpower/delrho 

      pcur=delcurr/delrho
      
cSmirnov970106 beg
      ppow_e=delpow_e/delrho
      ppow_i=delpow_i/delrho
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            ppow_s(kk)=delpow_s(kk)/delrho
         enddo
      endif
      
      ppow_cl=delpow_cl/delrho
cSmirnov970106 end

c------------------------------------------------------------
c     power (erg/sec), current(A)
c------------------------------------------------------------

c      write(*,*)'jbinmin,power(jbinmin),ppow,(hrho*jbinmin-rhomin)',
c     &jbinmin,power(jbinmin),ppow,(hrho*jbinmin-rhomin)

      power(jbinmin)=power(jbinmin)+ppow*(hrho*jbinmin-rhomin)     

c      write(*,*)'power(jbinmin)',power(jbinmin)

cSmirnov970106 beg
      power_e(jbinmin)=power_e(jbinmin)+ppow_e*(hrho*jbinmin-rhomin)
      power_i(jbinmin)=power_i(jbinmin)+ppow_i*(hrho*jbinmin-rhomin)
      power_cl(jbinmin)=power_cl(jbinmin)+ppow_cl*(hrho*jbinmin-rhomin)
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            power_s(jbinmin,kk)=
     &           power_s(jbinmin,kk)+ppow_s(kk)*(hrho*jbinmin-rhomin)
         enddo
      endif

cSmirnov970106 end
c      write(*,*)'p_c_prof delpow_e',delpow_e

      cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
     &(delpow_e/delrho)*(hrho*jbinmin-rhomin)/binvol(jbinmin)*
     &0.5d0*(eff_rho_min+
     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &                   (hrho*jbinmin-rhomin)/(delrho)))

      current(jbinmin)=
     1 current(jbinmin)+pcur*(hrho*jbinmin-rhomin)  !toroidal current from bin

c     write(*,*)'p_c_prof jbinmin'
c     write(*,*)'current(jbinmin)',current(jbinmin)
c     write(*,*)'cur_den_parallel(jbinmin)*binarea(jbinmin)',
c    &cur_den_parallel(jbinmin)*binarea(jbinmin)

   
c      write(*,*)'jbinmax,power(jbinmax),ppow,(hrho*jbinmax-rhomin)',
c     &jbinmax,power(jbinmax),ppow,(hrho*jbinmax-rhomin)
 
      power(jbinmax)=power(jbinmax)+ppow*(rhomax-hrho*(jbinmax-1))

c     write(*,*)' power(jbinmax)', power(jbinmax)
cSmirnov970106 beg
      power_e(jbinmax)=power_e(jbinmax)+ppow_e*(rhomax-hrho*(jbinmax-1))
      power_i(jbinmax)=power_i(jbinmax)+ppow_i*(rhomax-hrho*(jbinmax-1))
      power_cl(jbinmax)=power_cl(jbinmax)+
     1                  ppow_cl*(rhomax-hrho*(jbinmax-1))
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            power_s(jbinmax,kk)=
     &          power_s(jbinmax,kk)+ppow_s(kk)*(rhomax-hrho*(jbinmax-1))
         enddo
      endif
cSmirnov970106 end

      cur_den_parallel(jbinmax)=cur_den_parallel(jbinmax)+
     &(delpow_e/delrho)*(rhomax-hrho*(jbinmax-1))/binvol(jbinmax)*
     &0.5d0*(eff_rho_max+
     &   (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &   (rhomax-hrho*(jbinmax-1))/(delrho)))
     
      current(jbinmax)=
     1	  current(jbinmax)+pcur*(rhomax-hrho*(jbinmax-1)) !toroidal current from bin

c     write(*,*)'p_c_prof jbinmax'
c     write(*,*)'current(jbinmax)',current(jbinmax)
c      write(*,*)'cur_den_parallel(jbinmax)*binarea(jbinmax)',
c     &cur_den_parallel(jbinmax)*binarea(jbinmax)

c      write(*,*)'jbinmin,jbinmax',jbinmin,jbinmax
c     write(*,*)'power(jbinmin),power(jbinmax)',
c    &power(jbinmin),power(jbinmax)

c      write(*,*)'power_e(jbinmin)',
c     1           power_e(jbinmin)
c      write(*,*)'power_e(jbinmax)',
c     1           power_e(jbinmax)


      if(jbinmax.gt.(jbinmin+1)) then
         do j=(jbinmin+1),(jbinmax-1)
c            write(*,*)'j,power(j),ppow,hrho',j,power(j),ppow,hrho
            power(j)=power(j)+ppow*hrho
c            write(*,*)'power(j)',power(j)
cSmirnov970106 beg
            power_e(j)=power_e(j)+ppow_e*hrho
            write(*,*)'j,power_e(j),ppow_e,hrho',
     &                 j,power_e(j),ppow_e,hrho

            power(j)=power(j)+ppow*hrho
            power_i(j)=power_i(j)+ppow_i*hrho
            power_cl(j)=power_cl(j)+ppow_cl*hrho
            if (iabsorp.eq.3) then
               do kk=2,nbulk
                  power_s(j,kk)=power_s(j,kk)+ppow_s(kk)*hrho
               enddo
            endif

cSmirnov970106 end
            cur_den_parallel(j)=cur_den_parallel(j)+
     &      (delpow_e/delrho)*hrho/binvol(j)*
     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &      (hrho*(j-0.5d0)-rhomin)/delrho)
     
            current(j)=current(j)+pcur*hrho             !toroidal current from bins
         
c            write(*,*)'j,current(j),cur_den_parallel(j)*binarea(j)',
c     &      j,current(j),cur_den_parallel(j)*binarea(j)

c            write(*,*)'j,power(j)',j,power(j)

c            write(*,*)'j,power_e(j)',
c     .                 j,power_e(j)

         enddo
      endif

30    continue
      
      
c*******test
c      write(*,*)'allpower',allpower
c       write(*,*)'p_c_prof allcur', allcur
c      do ir=1,NR
c        write(*,*)'ir,power(ir),current(ir),temparr(ir),
c     1   spower(ir),scurrent(ir)',
c     2   ir,power(ir),current(ir),temparr(ir),
c     1   spower(ir),scurrent(ir)
c      enddo
c       do ir=1,NR
c        write(*,*)'ir,power(ir)',ir,power(ir)
c       enddo
c       spower_tot=0.d0
c       power_tot=0.d0
c       do ir=1,NR-1
c         spower_tot=spower_tot+spower(ir)
c         power_tot=power_tot+power(ir)
c       enddo
c       write(*,*)'delpwr(1)-delpwr(is)',delpwr(1)-delpwr(is)
c       write(*,*)'allpower,spower_tot,power_tot',
c     &            allpower,spower_tot,power_tot
c       write(*,*)'allpower,power_tot',
c     &            allpower,power_tot
c       write(*,*)'in p_c_prof before 99 return'
c       do ir=1,NR-1
c       write(*,*)'ir,cur_den_parallel(ir)',ir,cur_den_parallel(ir)
c       enddo
     
c99    return
      END

