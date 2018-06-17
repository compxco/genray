


c        **********************_dinit_mr************************
c        *                        -                            *
c        * this subroutine reads the data from genray.in       *
c        *******************************************************
c
c-------------------------dinit_mr---------------------------------
c        It creates data for multiple ray case. 
c        for spline approximations                       	   
c        density,temperature and Z_effective profiles.             
c        Calculates initial data for all rays.                     
c        It uses input data from genray.in or genray.dat 
c        file for multiple ray case.                               
c        Data from genray.in or genray.dat file were read          
c        previously in genray.f file   
c        It reads or creates the non-maxwellian electron distribution
c
c        it uses the following functions and subroutines           
c        bmod,spldens1                                             

c      	 This program reads data from genray.in file for	   !
c        multiple ray case.                                        !
c        It creates the files for spline approximations 	   !
c        density,temperature and Z_effective profiles.             !
c        it uses the following functions and subroutines           !
c        bmod,spldens1                                             !
c------------------------------------------------------------------
c        output parameters:					   !
c                          ndim-number of the ray-tracing equations!
c                          nray-number of the rays at antenna      !
c------------------------------------------------------------------

      subroutine dinit_mr(ndim,nray)

      implicit none

      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'three.i'
      include 'five.i'
      include 'cone.i'
      include 'grill.i'
      include 'rkutta.i'
      include 'six.i'
      include 'scatnper.i'
      include 'write.i'
      include 'writencdf.i'
      include 'onetwo.i'
      include 'output.i'
      include 'emissa.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

c-----output
      integer
     &ndim,  !number of the ray-tracing equations
     &nray   !number of the rays at antenna
c..............................................................
c     these two arrays are for namelist work only


c      dimension prof2(nbulka,ndensa)
c      dimension prof(nbulka*ndensa)
c     Following are working arrays for the namelist input of density 
c     and temperature,... at nonuniform radii grids
      real*8 prof2_nonuniform(nbulka,ndensa),prof_radii(nbulka,ndensa)
c      integer nj_tab(nbulka) !the number of profile points for each species
c..............................................................
      include 'dinit_nml.i'
c..............................................................
      integer nbulk1,i,j,k,iraystop,i1,j1,ifreq,ii,initial,nray1,icone,
     &imax

      real*8 trnspi,h,v0,w0,hfreq,delt,zefftest,zion,psi,denstot,
     &szini,szi2ni,stini,pressure,den_test,prestest,tem_test,
     &energy,pitch,fdist,dfdx,dfdpitch,dfdp,psi_mag,xi,yi,tetan,x0,
     &zconv,rconv,theta,rhoconv,phiconv,xconv,yconv,cnparopt,
     &h_rho,dens,temp,psi_loc,
     &tpop,
cSAP090311 for test
     &zeff_loc,vflow_loc

c-----externals
      real*8 b,psi_rho,prespsi,densrho,temperho,psif,x,y,tpoprho,
cSAP090311
     &zeffrho,vflowrho
cSAP080711BH080714
c-----genray input file: genray.in or genray.dat
      character*10 genray_in_dat

      iraystop=0
      pi=4*datan(1.d0)
      trnspi=pi/180.d0

c---------------------------------------------------------------
c     read all namelists from input genray.in or genray.dat file
c---------------------------------------------------------------
cSAP080731
c      call read_all_namelists(genray_in_dat,ndim,nray)
c      write(*,*)'in dinit_mr, genray_in_dat=', genray_in_dat
c-----If input file is genray.in, then it will change
c     input data to genray.dat file format 
c      if (genray_in_dat.eq.'genray.in')  then
c         write(*,*)'genray.f before transform_genray_in_to_dat'
c         call transform_genray_in_to_dat
c      endif
c---------------------------------------------------------------------
      if (n_wall.gt.0) then    
c-------calculate poloidal angles [at radian] and small radius
c       of wall and limiter points
c       thetapol_wall(i=1,..,n_wall)
c       thetapol_limiter(i=1,...,n_limiter(j),j=1,max_limiters)
c       rho_wall(i=1,..,n_wall)
c     rho_limiter(i=1,..,n_limiter(j),j=1,max_limiters))

c  
c       These arrays will be in common /fourb/ in file fourb.i
        call wall_limiter_theta_pol_rho

      endif

c------------------------------------------------------------------------
cyup      write(*,*)'dinit_mr: Absorption'
c------------------------------------------------------------------------
c     iabsorp=Imag(N_perp)
c     (imaginary part of the perpendicular refructive index)
c-------------------------------------------------------------------------
cyup      write(*,*)'dinit_mr: iabsorp=',iabsorp
c     iabsorp=1 for EC waves from Mazzucato solver
c     iabsorp=2 for LH waves
c     iabsorp=3 for FW waves
c-----------------------------------------------------
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'dinit_mr:  n_relt_harm1',n_relt_harm1
      write(*,*)'dinit_mr:  n_relt_harm',n_relt_harm
      write(*,*)'dinit_mr:  n_relt_harma',n_relt_harma
      endif ! outprint
      
      if (n_relt_harm1.eq.9999)then
        n_relt_harm1=-n_relt_harm
        n_relt_harm2= n_relt_harm
      else
        n_relt_harm2=n_relt_harm1+n_relt_harm-1
      endif
       
cSm060313
      if(n_relt_harm1.lt.n_relt_harm1a) then
         write(*,*)'n_relt_harm1<n_relt_harm1a'
         write(*,*)'it should be n_relt_harm1=>n_relt_harm1a'
         write(*,*)'please change n_relt_harm1a'        
         write(*,*)'in param.i and recompile the code'
         write(*,*)'n_relt_harm1,n_relt_harm1a',
     &               n_relt_harm1,n_relt_harm1a
         stop 'in dinit_mr' 
      endif
 
      if(n_relt_harm2.gt.n_relt_harma) then
         write(*,*)'n_relt_harm2>n_relt_harma'
         write(*,*)'it should be n_relt_harm2=<n_relt_harma'
         write(*,*)'please change n_relt_harma'
         write(*,*)'in param.i and recompile the code'
         write(*,*)'n_relt_harm2,n_relt_harm2a',
     &               n_relt_harm2,n_relt_harm2a
         stop 'in dinit_mr' 
      endif
      
cyup      write(*,*)'dinit_mr:  1 nbulk=',nbulk 
      do i=1,nbulk 
        if ((i_salphal(i).ne.1).and.(i_salphal(i).ne.0)) then
          write(*,*)'(i_salphal(i).ne.1).or.(i_salphal(i).ne.0)'
          write(*,*)'It should i_salphal(i) =0 or =1'
          write(*,*)'i=',i,'i_salphal(i)',i_salphal(i)
          write(*,*)'Please chagne i_saplhal(i)'
          write(*,*)'in input genray.in or genray.dat file'
          stop 'in dinit_mr:  i_salphal problem'
        endif
      enddo 
           
      
      if (n_contour_plot_disp.gt.n_contour_plot_disp_a)then
         write(*,*)'in dinit.f'
         write(*,*)'n_contour_plot_disp.gt.n_contour_plot_disp_a'
         write(*,*)'n_contour_plot_disp',n_contour_plot_disp
         write(*,*)'n_contour_plot_disp_a',n_contour_plot_disp_a
         write(*,*)'Please change n_contour_plot_disp_a in param.i'
         write(*,*)'and recompile the code'
         write(*,*)'or change n_contour_plot_disp in genray.dat'
         stop 'dinit.f dinit_mr'
      endif


c$$$      endif
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'dinit_mr:  Plasma parameters'
      write(*,*)'dinit_mr:  izeff=',izeff
      do i=1,nbulk
         write(*,*)'dinit_mr: i,temp_scale(i),den_scale(i)',
     .   i,temp_scale(i),den_scale(i)
      enddo 
      do i=1,nbulk
         write(*,*)'dinit_mr: i, dmas(i)',i,dmas(i) 
         ! dmas(i)=Mass(i)/Mass_electron
      enddo
      endif ! outprint

      do i=1,nbulk
         te0(i)=ate0(i)
         teb(i)=ateb(i)
      enddo
c-----------------------------------------------------
c      /species/
c   plasma component charges charge(i)=mod(charge(i)/charge_electron)
c----------------------------------------------------

c-----------------------------------------------------
c  plasma components mass dmas(i)=Mass(i)/Mass_electron
c-----------------------------------------------------
      do i=1,nbulk
         write(*,*)'dinit_mr: i, dmas(i)',i,dmas(i)
      enddo
c--------------------------------------------------------------

c     calculation of nbulk1
      if(((izeff.eq.0).or.(izeff.eq.2)).or.(izeff.eq.3)) then
c        izeff=0, zeff will be calculated using the given ions;
c                 the electron density will be calculated using ion's densities;
c             =1  ion's densities(i) i=nbulk and i=nbulk-1 will be calculated  using
c                 Zeff, the electon density and ion's densities(i), i=2,nbulk-1;
c        izeff=2, zeff will not coincide with the plasma components
c             =3  it uses eqdsk pres (pressure) and ions densities_i
c                 for i=2,... nbulk
c                 Let temperature T_E=T_i
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate Zeff(rho),
c                 dens_electron(rho) and T_e(rho)=T_i(rho)
c             =4  it uses eqdsk pres (pressure), zeff,ions densities
c                 for i=2,... nbulk-2 (nbulk>2) and T_e(rho),T_i(rho)
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate dens_electron(rho) and
c                 ion densities for i=nbulk and i=nbulk-1
         nbulk1=nbulk         
      else
c        (izeff=1 or izeff=4), zeff is given, the ions component will be calculated
         if (nbulk.le.2) nbulk1=nbulk
         if (nbulk.eq.2) then
	    write(*,*)'dinit_mr:  nbulk=2, Zeff must be equal charge(2)'
	    write(*,*)'Please check it or use the option izeff=0'
cSAP090801
c	    stop
	 endif
         if (nbulk.gt.2) nbulk1=nbulk-2
      endif !izeff
cyup      write(*,*)'nbulk1=',nbulk1
c------------------------------------------------------------------
      h=1.d0/(ndens-1)
      do i=1,ndens
        rhom(i)=h*(i-1)
      enddo
c------------------------------------------------------------------
c     The parameters for the density fluctuations
c     /varden/
c------------------------------------------------------------------
      if(idens.eq.0) then      
      
       if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'dinit_mr:  Analytical radial profiles'
c-dense_(i)=(dense0(i)-denseb(i))*(1-rho**rn1de(i))**rn2de(i)+denseb(i)
c------------------------------------------------------------------\     
	 do i=1,nbulk1
	   write(*,*)'dinit_mr: i, dense0(i)',i,dense0(i)
	 enddo

	 do i=1,nbulk1
	   write(*,*)'dinit_mr: i, denseb(i)',i,denseb(i)
	 enddo

	 do i=1,nbulk1
	   write(*,*)'dinit_mr: i, rn1de(i)',i,rn1de(i)
	 enddo

	 do i=1,nbulk1
	   write(*,*)'dinit_mr: i, rn2de(i)',i,rn2de(i)
	 enddo

	 do i=1,nbulk1
	   write(*,*)'dinit_mr: i, dense0(i)',i,dense0(i)
	 enddo
	 
	 endif ! outprint

c--------creation the array dens1(ndensa,nbulka)
cSm080118
         do i=1,nbulk1
            do k=1,ndens
               rho=rhom(k)
	       dens1(k,i)=(dense0(i)-denseb(i))*
     1	                  (1-rho**rn1de(i))**rn2de(i)+denseb(i)
            enddo
         enddo

	 do i=nbulk1,1,-1
	    do k=1,ndens
	       rho=rhom(k)

	       if (((izeff.eq.0).or.(izeff.eq.3)).and.(i.eq.1)) then
	          dens1(k,1)=0.d0
                do j=2,nbulk
                dens1(k,1)=dens1(k,1)+charge(j)*dens1(k,j)
                enddo
	       else
	          dens1(k,i)=(dense0(i)-denseb(i))*
     1	                     (1-rho**rn1de(i))**rn2de(i)+denseb(i)
cSm070426
c-----------------multiply density profiles by den_scale
                dens1(k,i)=dens1(k,i)*den_scale(i)
                if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
                  write(*,*)'dinit_mr:  i,k,dens1(k,i)',i,k,dens1(k,i)
                  write(*,*)'dense0(i),denseb(i),rho,rn1de(i),rn2de(i)',
     *            dense0(i),denseb(i),rho,rn1de(i),rn2de(i)
                endif ! outprint
               endif
            enddo
         enddo
cyup	 write(*,*)'dinit_mr: end of analytical density profiles input'

c------------------------------------------------------------------         
c         /tpopprof/
c         /vflprof/
c--------------------------------------------------------------------       
c        creation the arrays for analytical profoliles
c        tpop1(ndensa,nbulka), vflow1(ndensa,nbulka)
	 do i=1,nbulk
	    do k=1,ndens
	       rho=rhom(k)	
	       tpop1(k,i)=(tp0(i)-tpb(i))*
     1                    (1-rho**rn1tp(i))**rn2tp(i)+tpb(i)

               vflow1(k,i)=(vfl0(i)-vflb(i))*
     1                    (1-rho**rn1vfl(i))**rn2vfl(i)+vflb(i)
	    enddo
	 enddo

c---------------------------------------------------------------------
c        /zprof/          
	 if (izeff.eq.3) goto 10
c        /tprof/
c----------------------------------------------------------
cTemperature
ctempe_(i)=(te0(i)-teb(i))*(1-rho**rn1te(i))**rn2te(i)+teb(i)
c-----------------------------------------------------------
       if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         do i=1,nbulk
          write(*,*)'dinit_mr: i, te0(i)',i,te0(i)	
          write(*,*)'dinit_mr: i, teb(i)',i,teb(i)
         enddo
         do i=1,nbulk
          write(*,*)'dinit_mr: i, rn1te(i)',i,rn1te(i)
         enddo
         do i=1,nbulk
          write(*,*)'dinit_mr: i, rn2te(i)',i,rn2te(i)
         enddo
         do i=1,nbulk
          write(*,*)'dinit_mr: i,tp0(i),tpb(i)',i,tp0(i),tpb(i)
          write(*,*)'dinit_mr: rn1tp(i),rn2tp(i)',rn1tp(i),rn2tp(i)
         enddo
       endif ! outprint

c------- creation of array temp1(ndensa,nbulka)
	 do i=1,nbulk
	    do k=1,ndens
	       rho=rhom(k)
	       temp1(k,i)=(te0(i)-teb(i))*
     1                    (1-rho**rn1te(i))**rn2te(i)+teb(i)
cSm070426
c--------------multiply temperature profiles by temp_scale
               temp1(k,i)=temp1(k,i)*temp_scale(i)
	    enddo
	 enddo
10	 continue         


cyup         write(*,*)'dinit_mr: zeff0,zeffb,rn1zeff,rn2zeff'
cyup         write(*,*)zeff0,zeffb,rn1zeff,rn2zeff
	 if(((izeff.eq.1).or.(izeff.eq.2)).or.(izeff.eq.4)) then
c           the given analytical Zeff profile
c-------------------------------------------
c           zeff=(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff+zeffb
c-------------------------------------------
c           the creation of array zeff1(ndens)
	    do k=1,ndens
	       rho=rhom(k)
               zeff1(k)=(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff+zeffb
	    enddo
	 endif

	 if((izeff.eq.0).or.(izeff.eq.3)) then
	    if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
	    write(*,*)'dinit_mr: izeff=1 zeff will be calculated using 
     1 the given ions densities'
          endif ! outprint
	 endif
      endif ! idens analytical
cyup      write(*,*)'dinit_mr:  partner=',partner
c--------------------------------------------------------------------
      if (partner.eq.'genray_profs_in.txt' .or. 
     1    partner.eq.'genray_profs_in.nc') then
c--------------------------------------------------------------
c        read density ,temperature and zeff  profiles from  
c        the file created by ONETWO : genray_profs_in
c
c        Read in from external file
c        1) the density, temperature and zeff profiles
c          at genray small radial mesh:
c          dens1,temp1,zeff1,eqdskin
c        2) nbulk - the number of plasma species 
c        3) set izeff =2 (for nbulk.ge.1)
c-------------------------------------------------------------       
c        call read_transport_prof(ndens,charge,dmas,
c     &  nbulk,dens1,temp1,tpop1,vflow1,zeff1,izeff,eqdskin,partner)
c        write(*,*)'*******************************************'
c        write(*,*)'sub read_transport_prof changes izeff, eqdskin'
c        write(*,*)'izeff=',izeff
c        write(*,*)'eqdskin=',eqdskin
cSm070203
c        write(*,*)'dinit.f after read_transport_prof'
c        write(*,*)'nbulk ',nbulk
c        write(*,*)'charge',charge
c        write(*,*)'dmas',dmas
c        dmas(1)=1.d0
c--------multiply density profiles by den_scale(i)
         do i=1,nbulk
            if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
            write(*,*)'dinit_mr: i,temp_scale(i)',i,temp_scale(i)
            endif ! outprint
            do k=1,ndens               
               dens1(k,i)=dens1(k,i)*den_scale(i)
               temp1(k,i)=temp1(k,i)*temp_scale(i)
            enddo
        enddo
         
cyup        write(*,*)'*******************************************'
        go to 20   !   <<<<<==========
      endif
c------------------------------------

      if (idens.eq.1) then
c        spline approximation of the density, temperature, zeff,
c        tpop=T_perp/T_parallel, vflow
c        radial profiles
c        input of the arrays on the radial mesh	from dtzprof.dat
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'dinit_mr: idens=1: Before mult by den_scale'
         write(*,*)'nbulka,ndensa,nbulk,ndens',nbulka,ndensa,nbulk,ndens
         do i=1,nbulk
           write(*,*)'i',i,'dens1(k,i)'
           write(*,'(a,2i4,e12.4)')'dens1=',i,1,    dens1(1,i)
           write(*,'(a,2i4,e12.4)')'dens1=',i,ndens,dens1(ndens,i)
         enddo
         endif ! outprint
         !pause !!!
cSm070426 
c---------multiply density profiles by den_scale(i)
cBh080102          do i=i1,nbulk          
          do i=1,nbulk
	    do k=1,ndens
               dens1(k,i)=den_scale(i)*dens1(k,i)
            enddo
         enddo
        
c 21      format(5e16.9)
cyup	 write(*,*)'dinit_mr: nbulk1',nbulk1,'ndens',ndens
	 if ((izeff.eq.0).or.(izeff.eq.3)) then
	   i1=2
	 else
	   i1=1
	 endif

       if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'dinit_mr: idens=1:  After mult by den_scale'
         do i=i1,nbulk1
	    write(*,*)'i',i,'dens1(k,i)'
	    write(*,*)(dens1(k,i),k=1,ndens)
         enddo
       endif ! outprint
         
c----------------------------------------------------------
c        calculation of the electron density from
c        the charge neutrality 
c----------------------------------------------------------
	 if ((izeff.eq.0).or.(izeff.eq.3)) then
cSmirnov/00/05/26           
           do k=1,ndens
	     dens1(k,1)=0.d0
	     do j=2,nbulk
	       dens1(k,1)=dens1(k,1)+charge(j)*dens1(k,j)
	     enddo
           enddo
         endif

         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'dinit_mr:  dens1(k,1)'
         write(*,*)(dens1(k,1),k=1,ndens)
         endif ! outprint

cSm070426 
c--------multiply temperature profiles by temp_scale

         do i=1,nbulk
	    do k=1,ndens
c               write(*,*)'i,k,temp1(k,i),temp_scale(i)',
c     &                    i,k,temp1(k,i),temp_scale(i)
               temp1(k,i)=temp1(k,i)*temp_scale(i)
c               write(*,*)'temp1(k,i)',temp1(k,i)
            enddo
         enddo

         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         do i=1,nbulk
            write(*,*)'i',i,'temp1(k,i)'
            write(*,*) (temp1(k,i),k=1,ndens)
         enddo
c---------tpoptab
         do i=1,nbulk
            write(*,*)'i',i,'tpop1(k,i)'
            write(*,*) (tpop1(k,i),k=1,ndens)
         enddo
c--------vflowtab
         do i=1,nbulk
            write(*,*)'i',i,'vflow1(k,i)'
            write(*,*) (vflow1(k,i),k=1,ndens)
         enddo
         write(*,*)'dinit_mr: idens=1:  izeff=',izeff
         if(((izeff.eq.1).or.(izeff.eq.2)).or.(izeff.eq.4)) then
c           the given Zeff profile is in the table form
            write(*,*)'dinit_mr: idens=1:  zeff1(k)'
         else
            write(*,*)'dinit_mr: idens=1:  uniform zeff1',zeff1
         endif !izeff              
         endif ! outprint
      endif ! idens=1
      
c-----------------------------------------------------------
 20   continue !from:  if (partner.eq.
cyup      write(*,*)'dinit_mr:  after    20   continue '
c------------------------------------------------------------
c     read the data for emission calculations
c     /emission/
cyup       write(*,*)'dinit_mr: i_emission',i_emission
      if (i_emission.eq.1) then
c       if (nfreqa.lt.nfreq) then
c        write(*,*)'dinit_mr nfreqa<nfreq'
c        write(*,*)'it should be nfreqa.ge.nfreq'
c        write(*,*)'change nfreqa in param.i and recompile the code
c     +  or reduce nfreq in genray.in'
c        stop
c       endif
cSAP090808 delete jx_kin_a
c       if ((i_emission_spectrum.eq.1).and.(jx_kin_a.lt.jx_kin)) then
c        write(*,*)'dinit_mr jx_kina<jx_kin'
c        write(*,*)'it should be jx_kin_a.ge.jx'
c        write(*,*)'change jx_kin_a in param.i and recompile the code
c     +  or reduce jx_kin in genray.in'
c        stop
c       endif
 
       if ((wallr.lt.0.d0).or.(wallr.gt.1.d0)) then
        WRITE(*,*)'dinit_mr:  it should be {0=< wallr =<1}'
        WRITE(*,*)'but wallr=',wallr
        WRITE(*,*)'change wallr in genray.in'
        STOP
       endif         
 
cyup       write(*,*)'dinit_mr:  i_rrind',i_rrind
      endif ! (i_emission.eq.1) 
c---------------------------------------------------------
c     read the data for for EC cone vertex coordinates calculations
c     This case is used for the optimal OX mode conversion.
      
c      /ox/


      if(i_ox.eq.1) then
        istart=3
        prmt(3)=-prmt3 !to create the negative time
        i_vgr_ini=+1
        ireflm=1   
cyup        write(*,*)'dinit_mr:  i_ox.eq.1 prmt(3)',prmt(3)
      endif

      if(((i_ox.ne.0).and.(i_ox.ne.1)).and.(i_ox.ne.2)) then
         write(*,*)'dinit_mr:  i_ox can  =0 or =1 or =2'
         write(*,*)'in namelist /ox/ i_ox=',i_ox
         write(*,*)'please change i_ox in input file'
         stop 'in dinit /ox/'
      endif  
         
    

c------------------------------------------------------------
      if(i_emission.eq.0) then
c       no emission calculations, set one frequency
        nfreq=1
      endif
c-----------------------------------------------
c     for the emission multi frequency case
c     calculate the electron gyro-frequency freqncy0 at the plasma center
c---------------------------------------------
      bmod=b(yma,xma,0.d0) !TL
      freqncy0=28.0*b0*bmod !GHZ

      if ((i_emission.eq.0).or.(nfreq.eq.1)) then
c--------no emission or only one frequency in the emission calculations   
         v0=806.2/(frqncy*frqncy)
         w0=28.0*b0/frqncy 
c-----------------------------------------------------
c        set the arrays for v and w
         v(1)=v0
         w(1)=w0      
         do i=2,nbulk
            v(i)=v0*charge(i)**2/dmas(i)
            w(i)=w0*charge(i)/dmas(i)
cyup            write(*,*)'dinit_mr:  i charge(i),dmas(i),v(i),w(i)'
cyup     +      ,i,charge(i),dmas(i),v(i),w(i)
         enddo
      else
c----------------------------------------------
c        for the emission multi frequency case
c        calculate the electron gyro-frequency freqncy0 at the plasma center
c        and create the array of the frequencies wfreq()
c---------------------------------------------
c         bmod=b(yma,xma,0.d0) !TL
c         freqncy0=28.0*b0*bmod !GHZ
c         hfreq=(freq01-freq00)/dfloat(nfreq-1)
c         write(*,*)'bmod,b0,freq01,freq00,nfreq,hfreq',
c     +   bmod,b0,freq01,freq00,nfreq,hfreq
c         write(*,*)'freqncy0',freqncy0

c         do ifreq=1,nfreq
c-----------set the array for the frequenciesc
c            wfreq(ifreq)=freqncy0*(freq00+hfreq*(ifreq-1))  
c         enddo

c         write(*,*)'dinit.f nfreq',nfreq
c--------choose the frequency wfreq(ifre0) most close to the central
c        second harmonic 2*freqncy0
c         ifreq0=1
c         delt=dabs(wfreq(1)-2.d0*freqncy0)  
c         do ifreq=1,nfreq
c            if (delt.gt.dabs(wfreq(ifreq)-2.d0*freqncy0)) then
c               delt=dabs(wfreq(ifreq)-2.d0*freqncy0)
c               ifreq0=ifreq
c            endif
c            write(*,*)'dinit.f ifreq,wfreq(ifreq)',ifreq,wfreq(ifreq)
c         enddo
c         write(*,*)'dinit.f ifreq0,wfreq(ifreq0)',ifreq0,wfreq(ifreq0)
      endif
     
c---------------------------------------------------------
      if (partner.eq.'genray_profs_in.txt' .or. 
     1    partner.eq.'genray_profs_in.nc') goto 30
c-------------------------------------------------------------
      if ((izeff.eq.0).or.(izeff.eq.3))then
c---------------------------------------------------------
c        calculation of the table for the radial profile zeff1(ndens)
         call zeffcalc
c        zeff1(ndens) is in common six.i
cyup         do i=1,ndens
cyup         write(*,*)'i',i,'zeff1(i)',zeff1(i)
cyup         enddo
      endif !izeff=0
c---------------------------------------------------------
      if(izeff.eq.1) then
c---------------------------------------------------------
c        calculation of the table for the ion densities profiles
c---------------------------------------------------------
         if (nbulk.lt.3) then
            WRITE(*,*)'dinit_mr: nbulk.lt.3, Zeff must be equal 
     1	charge(2), control it and use the option izeff=0'
            STOP
         else
c           nbulk.ge.3
c           calculation of the tables for the radial profile
c           dens1(ndens,nbulk) and dens1(ndens,nbulk-1)
	    if( charge(nbulk).eq.charge(nbulk-1)) then
	      WRITE(*,*)'Warning in dinit_mr: nbulk(.ge.3)=',nbulk
	      WRITE(*,*)'in dinit: charge(nbulk)=charge(nbulk-1)'
	      WRITE(*,*)'it is impossible to find the ions densities'
	      WRITE(*,*)'change charge(nulk) or charge(nbulk-1)'
	      WRITE(*,*)'it should be charge(nulk)>charge(nbulk-1)'
	      WRITE(*,*)'or use the option izeff=0'
	      STOP
	    endif

            call denscalc
            !pause !!!
c------------------------------------------------
c for test
cyup            do i1=1,nbulk
cyup               write(*,*)'dinit_mr: after call denscalc i1=',i1
cyup               do j1=1,ndens
cyup	             write(*,*)'j1=',j1,'dens1(j1,i1)',dens1(j1,i1)
cyup               enddo
cyup            enddo
            
cyup	    do j1=1,ndens
cyup	       zefftest=0.d0
cyup	       zion=0.d0
cyup	       do i1=2,nbulk
cyup	          if(dens1(j1,1).ne.0.d0) then
cyup	             zefftest=zefftest+(dens1(j1,i1)/dens1(j1,1))*
cyup     1                        charge(i1)*charge(i1)
cyup		     zion=zion+charge(i1)*dens1(j1,i1)/dens1(j1,1)
cyup		  else
cyup		     write(*,*)'dinit_mr: dens1(j1,1)=0'
cyup		     zefftest=zefftest+1.d0*charge(i1)*charge(i1)
cyup		  endif
cyup	       enddo
cyup	       write(*,*)'j1',j1,'zefftest',zefftest,'zion',zion
cyup	    enddo
c end test
c------------------------------------------------

	 endif ! nbulk
      endif ! izeff=1
      
      if (izeff.eq.3) then
	 do i=1,nbulk
	    do k=1,ndens
	       rho=rhom(k)
	       psi=psi_rho(rho)
	       denstot=dens1(k,1)
	       do ii=2,nbulk
             denstot=denstot+dens1(k,ii)
	       enddo
	       temp1(k,i)=prespsi(psi)/denstot/(1.6d3)
	       if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
	       write(*,*)'dinit_mr:  izeff=3,k,rho',k,rho
	       write(*,*)'prespsi(psi),denstot',prespsi(psi),denstot
	       write(*,*)'in dinit i,temp1(k,i)',i,temp1(k,i)
	       write(*,*)'in dinit i,dens1(k,i)',i,dens1(k,i)
	       write(*,*)'in dinit i,dens1(k,1)',i,dens1(k,1)
	       endif ! outprint
	    enddo
	 enddo
      endif !izeff=3

      if(izeff.eq.4) then

cSAP090801 
c        if(nbulk.lt.3) then
c	  write(*,*)'in dinit: izeff=4 nbulk=',nbulk
c	  write(*,*)'in dinit: it should be nbulk>2, use another '
c	  write(*,*)'in dinit: izeff option '
c	  stop
c	else

c         nbulk.ge.3
c         calculation of the tables for the radial profiles
c         dens1(ndens,nbulk),dens1(ndens,nbulk-1) and
c         dens1(ndens,1)
	  call denscalp
	  !pause !!!
c test izeff=4	beg
cSAP090801
        if (nbulk.ge.3) then
	  do j=1,ndens
	    szini=0.d0	   !sum(i=2,nbulk){charge(i)*dens1(j,i)}
	    szi2ni=0.d0	   !sum(i=2,nbulk){charge(i)**2*dens1(j,i)}
	    stini=0.d0	   !sum(i=2,nbulk){temp1(j,i)*dens1(j,i})
	    rho=rhom(j)
	    psi=psi_rho(rho)
	    pressure=prespsi(psi)/1.6d3
	    do i=2,nbulk
	       szini=szini+charge(i)*dens1(j,i)
	       szi2ni=szi2ni+charge(i)*charge(i)*dens1(j,i)
	       stini=stini+temp1(j,i)*dens1(j,i)
	    enddo
	    prestest=temp1(j,1)*dens1(j,1)+stini
	    zefftest=szi2ni/dens1(j,1)
	    if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
	    write(*,*)'dinit_mr:  izeff=4,j,rho',j,rho
	    write(*,*)'pressure,prestest',pressure,prestest
	    write(*,*)'zeff,zefftest',zeff1(j),zefftest
	    write(*,*)'dens1(j,1),szini',dens1(j,1),szini
	    endif ! outprint
	  enddo !j
c test izeff=4	end
	endif !nbulk.ge.3
      endif !izeff=4

 30   continue ! if (partner.eq. ....) goto20

c     creation of the density,temperature,zeff and
c     tpop, vflow
c     spline coefficients
      call spldens1
cyup      write(*,*)'dinit_mr:  after spldens1'
c-----test printing density,temperature,zeff
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      do j=1,ndens
         write(*,*)'j,rhom(j)',j,rhom(j)
         do i=1,nbulk
           den_test=densrho(rhom(j),i)
           tem_test=temperho(rhom(j),i)
           write(*,*)'i,dens(i,rhom(j)),temp(i,rhom(j))',
     &     i,den_test,tem_test
         enddo
      enddo
      endif ! outprint
c------------------------------------------------------------
c     Reading or creation the non-maxwellian electron distribution
c     for the calculation the anti-hermitian relativistic tensor
c------------------------------------------------------------
c     i_diskf=0 no usage of the non-maxwellian electron distribution 
c     i_diskf=1 reading the file diskf
c     i_diskf=2 readinf the file netcdfnm.nc 
c     i_diskf=3 analitical calculation of the non-maxwellian disribution
c     i_diskf=4 analytic calculation of continuous non-Maxwellian distribution
c              with three temperatures in  three energy ranges 
c              (uses splines or tabulated distributions)
c     i_diskf=5 fully analytic calculation of continuous non-Maxwellian 
c               distribution with three temperatures in  three energy ranges.
c               Generally much faster than i_diskf=4.
c----------------------------------------------------
      initial=1
cyup      write(*,*)'dinit_mr:  before dskin i_diskf=',i_diskf
      call dskin(initial,energy,pitch,rho,fdist,dfdx,dfdpitch,
     .           dfdp,1) 
cyup      write(*,*)'dinit_mr aft dskin non-Maxwellian distribution was set'

c---------------------------------------------------------
      if (istart.eq.1) then
c-----------EC wave-----------
           if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
           write(*,*)'dinit_mr: ncone=',ncone
           write(*,*)'zst=',(zst(i),i=1,ncone),
     1        'rst=',(rst(i),i=1,ncone),
     1        'phist=',(phist(i),i=1,ncone)
           write(*,*)'dinit_mr: betast=',(betast(i),i=1,ncone),
     1          'alfast=',(alfast(i),i=1,ncone)
           write(*,*)'dinit_mr: na1=', na1,'na2=',na2,
     1           'powtot=',(powtot(i),i=1,ncone)
           write(*,*)'dinit_mr: alpha1=', alpha1,
     1           'alpha2=',(alpha2(i),i=1,ncone)
           endif ! outprint
      endif
      if (istart.eq.2) then
c-----------LH and FW---------
c           write(*,*)'zst=',zst,'rst=',rst,'phist=',phist
c           write(*,*)'cnphi=',cnphi,'cnteta=',cnteta
      end if
      do i=1,ncone
      phist(i)=phist(i)*trnspi
      betast(i)=betast(i)*trnspi
      alfast(i)=alfast(i)*trnspi
      enddo
c--------------------------------------------------------------
c  normalization of the start parameters
      do i=1,ncone
      zst(i)=zst(i)/r0x
      rst(i)=rst(i)/r0x
      enddo
c  normaliszation 'diskdisk' parameters  
      d_disk=d_disk/r0x
      rho_launching_disk=rho_launching_disk/r0x
      rho_focus_disk=rho_focus_disk/r0x
      sigma_launching_disk=sigma_launching_disk/r0x
c--------------------------------------------------------------
c      xst=rst*dcos(phist)
c      yst=rst*dsin(phist)
c-----------------------------------------------------------
c   printing of total magnetic field on the magnetic axis (for control)
cyup      write(*,*)'dinit_mr:  yma,xma',yma,xma
      bmod=b(yma,xma,0.d0)
cyup      write(*,*)'dinit_mr: magn. field on the magnetic axis bmag=',bmod
      psi_mag=psif(yma,xma)
cyup      write(*,*)'dinit_mr: psi on the magnetic axis psi_mag=',psi_mag
cyup       write(*,*)'dinit_mr:  istart,nray=',istart, nray
ctest
      if ((i_emission.eq.0).or.(nfreq.eq.1)) then
c--------no emission or only one frequency in the emission calculations   
         xi=x(yma,xma,0.d0,1)
         yi=y(yma,xma,0.d0,1)
cyup         write(*,*)'dinit_mr:  at magnetic axis Xe,Ye ',xi,yi
c---------------------------------------------------------------
c     the creation the data for the contours X_e=const,Y_e=const
c     B_tot,B_tor, B_pol on the plate (rho,theta) 
c     These data will have the form of the tables in 
c     xybrhoth.dat: rho(i),theta(j),xe,ye,(xe+ye*ye),bmod,bphi,
c     *              dsqrt(bz**2+br**2)
c      nrhomap=30
c      nthetmap=100
c      call mapxyb(nrhomap,nthetmap)
c      write(*,*)'mapxyb was created'
c      stop
      endif
c---------------------------------------------------------------
c     if ray start point is outside the plasma (ECR -case)
c     then:
c     determination of arrays :1)for ray coordinates on the ECR cone
c     alphaj(nray),betaj(nray)-spherical and cylindrical angles
c                              2)for wave power	angle distribution
c     powj(nray) -power flowing in the ray chanel at antenna
c     with normalized  Sum(i=1,nray)delpw0(i)=1
c---------------------------------------------------------------
      if (istart.eq.1) then
c---------EC waves---------------------------------
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)'dinit_mr:  istart=',istart
         write(*,*)'dinit_mr:  raypatt=',raypatt
         write(*,*)'dinit_mr:  ncone=',ncone
         endif ! outprint

         nray1=1                ! Counter for position in ray arrays
         do icone=1,ncone
            if (raypatt.ne.'toray') then  !.eq.genray -
            !  defined below as tetan=0.5d0*pi-betast(icone)
               tetan=0.5d0*pi-betast(icone) ! here - just for printout
            else ! .eq.toray - defined below as tetan=betast(icone)/trnspi
               tetan=betast(icone) ! just for printout
            endif
cyup            write(*,*)'alpha1,na1,na2,alpha2,phist,alfast,tetan',
cyup     1           alpha1(icone),na1,na2,alpha2(icone),phist(icone),
cyup     1           alfast(icone),tetan
     
c            if (raypatt.ne.'toray') then
            if (raypatt.eq.'genray') then
               alpha1(icone)=alpha1(icone)*trnspi
               alpha2(icone)=alpha2(icone)*trnspi
               tetan=0.5d0*pi-betast(icone) !Polar angle (radians)
cyup               write(*,*)'in dinit_mr before cone_ec'
               call cone_ec(alpha1(icone),na1,na2,alpha2(icone),
     1              phist(icone),alfast(icone),tetan,powtot(icone),nray,
     1              alphaj(nray1),betaj(nray1),powj(nray1))
               do i=nray1,(nray1-1)+nray
                  zstj(i)=zst(icone)
                  rstj(i)=rst(icone)
                  phistj(i)=phist(icone)
               enddo
cyup               write(*,*)'in dinit_mr after cone_ec: nray',nray
cyup               do i=nray1,(nray1-1)+nray
cyup                  write(*,*)' i,powj(i)',i,powj(i)
cyup                  write(*,*)' i,betaj(i),alphaj(i)',i,betaj(i),alphaj(i)
cyup               enddo
            endif
            
            if(raypatt.eq.'toray') then
               alfast(icone)=alfast(icone)/trnspi
               tetan=betast(icone)/trnspi
cyup               write(*,*)'dinit_mr:  bef raypat: tetan,alfast',
cyup     1              tetan,alfast(icone)
               call raypat(tetan,alfast(icone),alpha1(icone),cr,nray_in,
     1              gzone,mray,betaj(nray1),alphaj(nray1))
               
               nray=nray_in
               if(nray.gt.nraymax)then
                  if(myrank.eq.0)then
                  WRITE(*,*)'dinit_mr: nray>nraymax. '
                  WRITE(*,*)'dinit_mr: Increase nraymax in param.i'
                  endif
                  STOP
               endif
               do i=nray1,(nray1-1)+nray
                  zstj(i)=zst(icone)
                  rstj(i)=rst(icone)
                  phistj(i)=phist(icone)
c                  write(*,*)'zstj(i),rstj(i)=',
c     +            zstj(i),rstj(i)
               enddo
               
cyup               do i=nray1,(nray1-1)+nray
cyup                  write(*,*)'Raypat ray starting angles (degrees):'
cyup                  write(*,*)' i,betaj(i),alphaj(i)',i,betaj(i),alphaj(i)
cyup               enddo
               do i=nray1,(nray1-1)+nray
                  powj(i)=(powtot(icone)/nray)*1.e13
                  alphaj(i)=alphaj(i)*trnspi
                  betaj(i)=(90.-betaj(i))*trnspi
               enddo
cyup               write(*,*)'in dinit_mr after raypat powj(i)'
cyup               do i=nray1,(nray1-1)+nray
cyup                  write(*,*)' i,powj(i)',i,powj(i)
cyup                  write(*,*)' i,betaj(i),alphaj(i)',i,betaj(i),alphaj(i)
cyup               enddo
            endif               !(On raypatt)
             
            if (raypatt.eq.'diskdisk') then
               call disk_to_disk_rays_initial_launching_data(nray)
            endif !diskdisk
c           stop 'dinit.f after disk_to_disk_rays_initial_launching_dat'

            if (raypatt.eq.'diskbeam') then

cyup               write(*,*)'dinit_mr:  before'
cyup               write(*,*)'disk_beam_rays_initial_launching_data'

               call disk_beam_rays_initial_launching_data(nray)

cyup               write(*,*)'dinit_mr:  after'
cyup               write(*,*)'disk_beam_rays_initial_launching_data'

            endif !diskdisk
c            stop 'dinit.f after disk_beam_rays_initial_launching_dat'           

cSAP050831            powtott=powtott+powtot(icone)
            powtott=powtott+powtot(icone)*1.e13   !(erg/sec)
            nray1=nray1+nray
         enddo                  !(On icone)
         nray=nray1-1
cyup        write(*,*)'dinit_mr: after EC starting conditions, total nray=',
cyup     1        nray
      endif                     !(On istart.eq.1, EC)

      if (istart.eq.2) then
c---------LH or WF----------
cyup	  write(*,*)'dinit_mr: before grill_lh '
cyup	  write(*,*)'ngrilla,ngrill,xma,yma',ngrilla,ngrill,xma,yma
	  !pause
c--------------------------------------------------------------
          call grill_lh(rhopsi0,ngrilla,ngrill,thgrill,phigrill,
     1    height,nthin,nthinmax,
     1    anmin,anmax,nnkpar,powers,powtott,
     1    antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1    n_theta_pol,
     1    xma,yma,psimag,
cSm050309
     1    fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,nnkprmax,
     1    anztorin,anzpolin,pwcpl_tp,nnktormax,nnkpolmax,
     1    anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1    nray,arzu0,arru0,arphiu0,arntheta,arnphi,powinilh,
     1    nraymaxl,wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1    ilaunch,r0launch,phi0launch,z0launch,
cSAP111030
     1    rlaunch,zlaunch,weight_power_launch,
     1    i_grill_pol_mesh,
     1    i_grill_npar_ntor_npol_mesh)
cSm050309
c     1    i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh)

cyup	   write(*,*)'dinit_mr: after grill nray=',nray

c          do iray=1,nray         
c             write(*,*)'iray,arzu0(iray),arru0(iray),arphiu0(iray)',
c     &                  iray,arzu0(iray),arru0(iray),arphiu0(iray)
c          enddo

      endif !LH or FW

      if (istart.eq.3) then
c---------ECR O_X_EBW mode conversion case----------
c         It uses the wave input data from the grill form
c         It sets i_n_poloidal=1
c	  It calculates: 
c	  1) the space point M_conv=(rhoconv,theta)=(zconv,rconv)
c            at which x_e=x0=1. The value of the poloidal
c            angle theta (degree) is given bellow
c            if theta=0  point is on the outer part of the magnetic surface
c            if theta=90 point is on the top of the magnetic surface
c         2) the optimal value N_parallel_optimal=cnparopt
c            for O_X mode conversion
c         3) gives
c            anmin(1)=cnparopt-0.01d0 
c            anmax(1)=cnparopt+0.01d0 
c         4) call grill_lh to calculate the initial data
c         5) gives the coordinates of the initial O_X mode conversion poin 
c            arzu0(1)=zconv
c            arru0(1)=rconv
c-------------------------------------------------------------
          i_n_poloidal=1
c         these lines are to create the start point inside the 
c         plasma for the given point xe=x0
cSAP050510	  x0=1.d0
          x0=1.002d0

c          theta=0.d0   !poloidal angle  (degree)
c          theta=-30.d0
          theta=thgrill(1)
cyup          write(*,*)'dinit before owconvr theta=',theta
          call owconvr (theta,x0,rhoconv,zconv,rconv)
cyup    	  write(*,*)'dinit_mr: rhoconv,zconv,rconv',rhoconv,zconv,rconv
	  rhopsi0(1)=rhoconv
          phiconv=0.d0
          bmod=b(zconv,rconv,phiconv)
          xconv=x(zconv,rconv,0.d0,1)
          yconv=y(zconv,rconv,0.d0,1)
c---------calculation of the optimal value N_parallel_optimal
c         for O_X mode conversion
	  cnparopt=dsqrt(yconv/(1.d0+yconv))
cyup    	  write(*,*)'dinit_mr: xconv,yconv,cnparopt',xconv,yconv,cnparopt
c         write(*,*)'dinit old value of rhopsi0(1)',rhopsi0(1) 
	  rhopsi0(1)=rhoconv
cyup          write(*,*)'dinit_mr: new rhopsi0(1)',rhopsi0(1) 
cyup          write(*,*)'dinit_mr: old anmin(1),anmax(1)',anmin(1),anmax(1)
          anmin(1)=cnparopt-0.01d0 
          anmax(1)=cnparopt+0.01d0 
cyup          write(*,*)'dinit_mr: new anmin(1),anmax(1)',anmin(1),anmax(1)
c---------------------------------------------------------------
cyup          write(*,*)'dinit_mr: before grill_lh ngrilla,ngrill',
cyup     1    ngrilla,ngrill
          call grill_lh(rhopsi0,ngrilla,ngrill,thgrill,phigrill,
     1    height,nthin,nthinmax,
     1    anmin,anmax,nnkpar,powers,powtott,
     1    antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1    n_theta_pol,
     1    xma,yma,psimag,
cSm050309
     1    fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,nnkprmax,
     1    anztorin,anzpolin,pwcpl_tp,nnktormax,nnkpolmax,
     1    anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1    nray,arzu0,arru0,arphiu0,arntheta,arnphi,powinilh,
     1    nraymaxl,wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1    ilaunch,r0launch,phi0launch,z0launch,
cYuP+BH120424_cSAP111030
     1    rlaunch,zlaunch,weight_power_launch,
     1    i_grill_pol_mesh,
     1    i_grill_npar_ntor_npol_mesh)
cSm050309
c     1    i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh)


c for ecr internal case 
          arzu0(1)=zconv
          arru0(1)=rconv
c--end ecr internal case
     
	   !write(*,*)'in dinitmr after grill'
      endif ! O_X_EBW  mode conversion case

ctest density profile
      imax=101
      h_rho=1.d0/dfloat(imax-1)
      if(myrank.eq.0) open(88,file='dens.bin',form='unformatted')
c      write(*,*)'imax,h_rho',imax,h_rho
 200  format(2(1pe10.3))

      do i=0,imax-1
         rho=h_rho*i
         dens=densrho(rho,1)
         temp=temperho(rho,1)
         tpop=tpoprho(rho,1)
         !write(*,*)'i,rho,dens,temp,tpop',i,rho,dens,temp,tpop
cSAP090311
         zeff_loc=zeffrho(rho)
         vflow_loc=vflowrho(rho,1) 
c         write(*,*)'zeff_loc,vflow_loc',zeff_loc,vflow_loc
c         write(*,*)'i,rho,dens,temp',i,rho,dens,temp
         psi_loc=psi_rho(rho)
         if(myrank.eq.0) write(88)real(rho),real(dens),real(temp)
      enddo
      
      if(myrank.eq.0)then
         write(88) 
         close(88)
      endif ! myrank=0
      
ctest temperature profiles
c       rho_small=0.9904922704086799d0
c       do j=1,nbulk
c          temp_test=temperho(rho_small,j)
c          write(*,*)'dinit j,rho_small,temp_test',j,rho_small,temp_test
c       enddo 
c      imax=10001
c      h_rho=1.d0/dfloat(imax-1)
c      do j=1,nbulk
c      write(*,*)'plasma component number j ',j 
c       do i=0,imax
c         rho=h_rho*i
c         dens=densrho(rho,j)
c         temp=temperho(rho,j)
c         write(*,*)'j,i,rho,dens,temp',j,i,rho,dens,temp
c       enddo !i
c      enddo  !j
c      stop 'dinit_mr test temperature profiles'
c----endtest
cendtest density
c      stop
c      write(*,*)'in dinit_mr i_diskf',i_diskf

c-----allocate pointers at writencdf.i and write_i
cyup      write(*,*)'dinit_mr:  before ainalloc_writencdf nray',nray
      !if(outnetcdf.eq.'enabled')then !YuP[2018-01-17] Added
      call ainalloc_writencdf_i(nray)
      !endif ! outnetcdf
cBH130508      call ainalloc_write_i(nray)
 
      if ((i_emission.gt.0).and.(nfreq.gt.1)) then
c----------------------------------------------
c        for the emission multi frequency case
c        calculate the electron gyro-frequency freqncy0 at the plasma center
c        and create the array of the frequencies wfreq()
c---------------------------------------------
         bmod=b(yma,xma,0.d0) !TL
         freqncy0=28.0*b0*bmod !GHZ
         hfreq=(freq01-freq00)/dfloat(nfreq-1)
cyup         write(*,*)'dinit_mr: bmod,b0,freq01,freq00,nfreq,hfreq',
cyup     +   bmod,b0,freq01,freq00,nfreq,hfreq
cyup         write(*,*)'dinit_mr: freqncy0',freqncy0

         do ifreq=1,nfreq
c-----------set the array for the frequencies
            wfreq(ifreq)=freqncy0*(freq00+hfreq*(ifreq-1))  
         enddo

cyup         write(*,*)'dinit_mr: nfreq',nfreq
c--------choose the frequency wfreq(ifre0) most close to the central
c        second harmonic 2*freqncy0
         ifreq0=1
         delt=dabs(wfreq(1)-2.d0*freqncy0)  
         do ifreq=1,nfreq
            if (delt.gt.dabs(wfreq(ifreq)-2.d0*freqncy0)) then
               delt=dabs(wfreq(ifreq)-2.d0*freqncy0)
               ifreq0=ifreq
            endif
cyup            write(*,*)'dinit_mr:  ifreq,wfreq(ifreq)',ifreq,wfreq(ifreq)
         enddo
cyup         write(*,*)'dinit_mr: ifreq0,wfreq(ifreq0)',ifreq0,wfreq(ifreq0)
      endif

      return ! dinit_mr: 
      end  !==============================================================

c        **********************dinit_1ray**********************
c        * this subroutine                                    *
c        * 1)if ray was launched outside the plasma           *
c        *   (istart=1, ECR wave)			      *
c        *   it detemines the point where the ray inter-    *
c        *   sects the plasma boundary(zu0,ru0,phiu0),and the *
c        *   tangent components of the refractive index	      *
c        *   cnteta,cnphi} .				      *
c        *   \hat{theta}=grad(psi) x \hat{phi}/abs(grad(psi)).*
c        *   Here \hat{theta} is the poloidal unit vector.    *
c        *   Here x is the vector product.		      *
c        *   If istart.ne.1 it launch the ray                 *
c        *   in the given point inside the plasma .	      *
c        * 2)	Then it calculates the parallel to magnetic   *
c        *   field component of the refractive index  cnpar.  *
c        *      Then it determinates the normal   to   magne- *
c        *   tic surface component of refractive index cnrho. *
c        *   It is directed inside the plasma		      *
c        *      Then it determinates the initial components   *
c        *   of the refractive index:			      *
c        *              cnz,cnr,cm=cnphi*r		      *
c	 *   Then it calculate and prints the value of        *
c        *   dispersion function d=eps                        *
c        *   The result is  the initial condition  for the    *
c        *   ray equations				      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c        zst,rst,phist                                            !
c	 EC wave:						   !
c        alfast- toroidal angle(radian) of the ray at antenna      !
c        betast- angle(radian) between the  horizontal             !
c                          plane and the ray at antenna            !
c        LH and FW waves:					   !
c        cnteta,cnphi						   !
c        output parameters					   !
c        u(1)=zu0,u(2)=ru0,u(3)=phiu0,u(4)=cnz,u(5)=cnr,u(6)=cm	   !
c        iraystop- index to stop the ray determination(=1)	   !
c       	   or make the ray determination (=0)		   !
c------------------------------------------------------------------
c        it uses the following functions and subroutines           !
c        ias1r,bmod,y,x,gamma1,s,abc,hamilt1,plasmaray,ninit_ec     !
c        cinit                                                     !
c------------------------------------------------------------------
      subroutine dinit_1ray(zst,rst,phist,alfast,betast,
     1                      cnteta,cnphi,u,iraystop)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'three.i'
      include 'five.i'
      include 'grill.i'
      include 'write.i'
      include 'rkutta.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      dimension u(6),deru(6)

cfor test_to plot D(ReN_perp,ImN_perp)
      integer n_param
      parameter (n_param=3)
      character*(15) name_param(n_param) !names of the input parameters
      real param(n_param)               !values of the input parameters     
cendtest

c-----for hot plasma roots
      double precision 
     &N_perp_root_ar(n_hot_roots_a)     !hot plasma roots
c---------------------------------------------------------------
c------for  rho_ini_LHFW
      integer  i_rho_ini_LHWH_found


      iraystop=0
      pi=4*datan(1.d0)
      trnspi=pi/180.d0
c---------------------------------------------------------------
c     if ray start point is outside the plasma (ECR -case)
c     then:
c          determine:    1) where the vacuum ray
c                           intersects the plasma boundary
c                        2) cnteta,cnphi-toroidal and poloidal
c                           components of the vacuum
c                           refractive index
c---------------------------------------------------------------
      call set_output ! initialize output.i
                      ! for each new ray
      
c      if (i_ox.eq.1) then
          call set_oxb          ! initialize oxb.i
                                ! for each new ray
c      endif
      if (istart.eq.1) then
c--------EC wave
cyup         write(*,*)' bef plasmray zst,rst,phist,alfast,betast',
cyup     1	 zst,rst,phist,alfast,betast
         call plasmray(zst,rst,phist,alfast,betast,
     1                  zu0,ru0,phiu0,iraystop)
cyup         write(*,*)'in dinit_1ray after plasmaray zu0=',zu0,'ru0=',ru0,
cyup     1   'phiu0=',phiu0,'iraystop=',iraystop
         if (iraystop.eq.1) then
	    return
         end if
c	 -------------------------------
c        the shift of the initial point inside the plasma from the boundary
	 z=zu0
	 r=ru0
cSAP091127
         phi=0.d0
         bmod=b(z,r,phi)
cyup         write(*,*)'befor edgcor z,r,rho',z,r,rho
         call edgcor(z,r,zu0,ru0)
c        end of the shift
cyup	 write(*,*)'in dinit_1ray after initial point shift'
c	 write(*,*)'ru0,zu0',ru0,zu0
cSAP091127
cRobtAndre140412         bmod=b(z0,r0,phi)  BH: but no effective change.
         bmod=b(zu0,ru0,phi)
cRobtAndre140412         write(*,*)'after edgcor z0,r0,rho',z0,r0,rho
cyup         write(*,*)'after edgcor z0,r0,rho',zu0,ru0,rho
c	 -------------------------------
c        nx,ny,nz in start point
         cnzst=dsin(betast)
         cnxst=dcos(betast)*dcos(alfast+phist)
         cnyst=dcos(betast)*dsin(alfast+phist)

cyup         write(*,*)'dinit betast,alfast+phist',betast,alfast+phist
c         write(*,*)'cnxst=',cnxst,'cnyst=',cnyst,'cnzst=',cnzst
c----------------------------------------------------------------
	 bmod=b(zu0,ru0,phiu0)
c----------------------------------------------------------------
c        ninit_ec creates the tangent to magnetic surface
c        components  of the refractive index cnteta,cnp                 
c        in the initial point (zu0,ru0,phiu0) for ECR wave
c-----------------------------------------------------------------
c         write(*,*)'in dinit_1ray before ninit_ec bz,br,bphi,bmod'
c         write(*,*)bz,br,bphi,bmod
cyup         write(*,*)'dinit_1ray cnxst,cnyst,cnzst',cnxst,cnyst,cnzst
cyup         write(*,*)'dinit_1ray dsqrt(cnxst**2+cnyst**2+cnzst**2)',
cyup     &dsqrt(cnxst**2+cnyst**2+cnzst**2)
         call ninit_ec(zu0,ru0,phiu0,cnxst,cnyst,cnzst,cnteta,cnphi)
cyup	 write(*,*)' after ninit_ec cnteta,cnphi',cnteta,cnphi

      endif

      if ((istart.eq.2).or.(istart.eq.3)) then
c--------LH and FW wave, OX-conversion pt.
         zu0=zst
         ru0=rst
         phiu0=phist
         bmod=b(zu0,ru0,phiu0)
         xe=x(zu0,ru0,phiu0,1)
         ye=y(zu0,ru0,phiu0,1)
         if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
         write(*,*)' in dinit_1ray istart=',istart
         write(*,*)'magnetic field in initial point'
         write(*,*)'bmod =',bmod,'zu0,ru0,phiu0',zu0,ru0,phiu0
         write(*,*)'in dinit_1ray ksi_nperp', ksi_nperp
ctest for initial conditions
         write(*,*)'dinit_1ray initia values xe,ye,uh'
         write(*,*)xe,ye,dsqrt(xe+ye*2)
         write(*,*)'dinit_1ray: istart,ioxm',  istart,ioxm
         endif ! outprint
cendtest
      end if
c--------------------------------------------------------------
      z=zu0
      r=ru0
      phi=phiu0
      cosphi=dcos(phiu0)
      sinphi=dsin(phiu0)

c---------------------------------------------------------------
c     calculations of the parallel (to the magnetic field)
c     refractive index component cnpar
c     N_par=(B_phi*N_phi+N_theta*e_theta*{e_z*B_z+e_r*B_r})/bmod
c     e_theta=(e_z*dpsi/dr-e_r*dpsi/dz)/abs(grad(psi))
c--------------------------------------------------------------
      ppp=(dpdzd*dpdzd+dpdrd*dpdrd)
      gradpsi=dsqrt(dpdzd*dpdzd+dpdrd*dpdrd)
      cnpar1=(cnphi*bphi+cnteta*(bz*dpdrd-br*dpdzd)/gradpsi)/bmod
      cnpar2=cnpar1**2
c Smirnov 961210 beg
c     test of the initial conditions
c      btheta0=(bz*dpdrd-br*dpdzd)/gradpsi
c      write(*,*)'in dinit1_ray btheta0',btheta0,'bphi',bphi,'bmod',bmod
c      write(*,*)'in dinit1_ray cnpat1,cnpar2',cnpar1,cnpar2
c Smirnov 961210 end
c-------------------------------------------------------------------
       if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
       write(*,*)'in 1ray cnpar1',cnpar1,'cnteta',cnteta,'cnphi',cnphi
       write(*,*)'i_rho_find_hot_nperp_roots',i_rho_find_hot_nperp_roots
       endif ! outprint
c---------------------------------------------------------------
      if (i_rho_find_hot_nperp_roots.eq.1) then
c-----------------------------------------------------------------
c       finds the small radius rho_ini > rho_min_find_hot_nperp_roots
c       at the vector rho where
c       hot plasma dispersdion function D_hot(npar) has three roots.
c       The vector rho^ is starting at the edge point (r_edge,z_edge,phi_edge),
c       and directed to the magnetic axis O(xma,yma,phi_edge)
c-------------------------------------------------------------------
cyup       write(*,*)'dinit.f r,z,phi,cnpar1', r,z,phi,cnpar1
 
        call rho_ini_hot_nperp_roots(r,z,phi,cnpar1)     
c     &  rho_ini,z_ini,r_ini)
        stop 'dinit.f after call rho_ini_hot_nperp_roots'
      endif
c------------------------------------------------------------------    

c-------------------------------------------------------------
c     fit the initial value of rho_ini fort LH or FW cutoff
c--------------------------------------------------------------
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'dinit.f i_rho_cutoff ',i_rho_cutoff
      write(*,*)'cnteta,cnphi',cnteta,cnphi  
      write(*,*)'dinit.f before i_rho_cutoff=1 z,r',z,r
      endif ! outprint
      if (i_rho_cutoff.eq.1) then
         costheta=(r-xma)/dsqrt((z-yma)**2+(r-xma)**2)
         sintheta=(z-yma)/dsqrt((z-yma)**2+(r-xma)**2)
         if(sintheta.gt.0d0) then
            theta=dacos(costheta)
         else
            theta=2*pi-dacos(costheta)
         endif  

c-------rho_ini_LHFW changes input arguments n_theta_pol,n_toroidal
c       New names n_theta_pol=cnteta,n_toroidal=cnphi
c       were used as arguments of rho_ini_LHFW to conserve 
c       the values of the variables cnteta and cnphi
c       after  call rho_ini_LHFW
cSAP091017
         n_theta_pol=cnteta
cSAP091026
         n_toroidal=cnphi

cyup         write(*,*)'dinit_1ray before rho_ini i_n_poloidal,n_theta_pol'
cyup     &   ,i_n_poloidal,n_theta_pol
cyup         write(*,*)'cnteta,cnphi',cnteta,cnphi

         call rho_ini_LHFW(theta,phi,cnpar1,
cSAP091026
c     &   i_n_poloidal,n_theta_pol,cnphi,
     &   i_n_poloidal,n_theta_pol,n_toroidal,
     &   rho_ini,z_ini,r_ini,cntheta_ini,cnphi_ini,
     &   i_rho_ini_LHFW_found)

cyup         write(*,*)'dinit.f after rho_ini_LHFW i_rho_ini_LHFW_found=',
cyup     &              i_rho_ini_LHFW_found

         if (i_rho_ini_LHFW_found.eq.1) then
c-----------cutoff point with new z,r coordinates was found
            z=z_ini 
            r=r_ini
cyup            write(*,*)'dinit.f after  rho_ini_LHFW z,r',z,r
         else
c-----------cutoff point was not found
            write(*,*)'cutoff point was not found'
            iraystop=1
            return
         endif
      endif

c      write(*,*)'dinit.f after rho_ini_LHFW,cntheta_ini,cnphi_ini',
c     &cntheta_ini,cnphi_ini
               
      if (i_rho_cutoff.eq.0) then
         cntheta_ini=cnteta
         cnphi_ini=cnphi
      endif

cyup      write(*,*)'dinit.f after rho_ini_LHFW,cntheta_ini,cnphi_ini',
cyup     &cntheta_ini,cnphi_ini

c--------------------------------------------------------------------
c     cninit solves the dispersion relation N=N(n_par)
c     Then subroutine calculates the initial components
c     of the refractive index  cnz,cnr,cm
c---------------------------------------------------------
cyup      write(*,*)'dinit z,r,phi,cnpar1,cntheta_ini,cnphi_ini',
cyup     &z,r,phi,cnpar1,cntheta_ini,cnphi_ini

cBH070123 start
      psi=fpsi(r,z)
      rho=rhopsi(psi)
      bmod=b(z,r,phi)
!      Rho= 9.535660582E-01  Saveliev starting condition (temporary)
cyup      write(*,*)'dinit r,z,rho,dens,temp ',r,z,rho,densrho(rho,1),
cyup     +          temperho(rho,1)

cSAP090518
      do i=1,nbulk
cyup        write(*,*)'dinit i',i
        dens_i=dense(z,r,phi,i)
cyup        write(*,*)'dinit i,dens_i',i,dens_i
      enddo
      x_e=x(z,r,phi,1)
      y_e=y(z,r,phi,1)
cyup      write(*,*)'x_e,y_e',x_e,y_e

      if(nbulk.ge.2) then
        x_2=x(z,r,phi,2)
        y_2=y(z,r,phi,2)
cyup        write(*,*)'x_2,y_2',x_2,y_2
      endif

      if (nbulk.ge.3)then
        x_3=x(z,r,phi,3) 
        y_3=y(z,r,phi,3)
cyup        write(*,*)'x_3,y_3',x_3,y_3
      endif

cyup      write(*,*)'cnpar1',cnpar1
c      w_cut_d_w_p=y_2-x_e/(y_e*(cnpar1**2-1))
c      w_cut_d_w_m=-y_2+x_e/(y_e*(cnpar1**2-1))
c      write(*,*)'w_cut_d_w_p',w_cut_d_w_p
c      write(*,*)'w_cut_d_w_m',w_cut_d_w_m
c      eps_m_g=1.d0-(x_e/y_e)/(1-y_2)
c      eps_p_g=1.d0+(x_e/y_e)/(1+y_2)
c      write(*,*)'eps_m_g,eps_p_g',eps_m_g,eps_p_g
      if(nbulk.eq.2) then
         w_lh_d_w=dsqrt(x_2)
cyup         write(*,*)'w_lh_d_w',w_lh_d_w
       endif
       if(nbulk.eq.3) then
         w_lh_d_w=dsqrt(x_2+x_3)
cyup         write(*,*)'w_lh_d_w',w_lh_d_w
       endif
cBH070123 end

ctest
      ioxm_n_npar_loc= ioxm_n_npar
      ioxm_n_npar=1
       call nper_npar_ioxm_n_npar(2,z,r,phi,cnpar1,
     & cnper,iraystop) ! ioxm_n_npar will be set in one.i
cyup      write(*,*)'iraystop,ioxm_n_npar,cnper',iraystop,ioxm_n_npar,cnper

cyup      write(*,*)'dinit.f subroutine dinit_1ray ifreq_write',ifreq_write

      if(iraystop.eq.0) then
         wn_perp_ioxm_p(ifreq_write)=cnper
      else
         wn_perp_ioxm_p(ifreq_write)=-1.d0
      endif

       ioxm_n_npar=-1
       call nper_npar_ioxm_n_npar(2,z,r,phi,cnpar1,
     & cnper,iraystop) ! ioxm_n_npar will be set in one.i
cyup      write(*,*)'iraystop,ioxm_n_npar,cnper',iraystop,ioxm_n_npar,cnper

cyup      write(*,*)'dinit.f subroutine dinit_1ray nfreq',nfreq

cyup      write(*,*)'dinit.f subroutine dinit_1ray ifreq_write',ifreq_write

      if(iraystop.eq.0) then
         wn_perp_ioxm_m(ifreq_write)=cnper
      else
         wn_perp_ioxm_m(ifreq_write)=-1.d0
      endif
      ioxm_n_npar= ioxm_n_npar_loc

cyup      write(*,*)'dinit.f i_look_roots',i_look_roots
      if (i_look_roots.eq.1)then   
c-----------------------------------------------------------------
c       plot ReD_hot(nperp) at given npar
c-----------------------------------------------------------------
        name_param(1)='xe'
        name_param(2)='ye'
        name_param(3)='npar'
        param(1)=xe
        param(2)=ye
        param(3)=cnpar1
        write(*,*)'dinit.f before map_d_hot n_param', n_param
        write(*,*)'name_param ',name_param 
        write(*,*)'param ', param
        call map_d_cold(z,r,phi,cnpar1,n_nperp_plot,cnperp_plot_min,
     &  cnperp_plot_max,
     &  name_param,param,n_param)   
        call  map_d_hot(z,r,phi,cnpar1,n_nperp_plot,cnperp_plot_min,
     &  cnperp_plot_max,
     &  name_param,param,n_param)
        call  map_d_hot(z,r,phi,cnpar1,n_nperp_plot,cnperp_plot_min,
     &  cnperp_plot_max,
     &  name_param,param,n_param)
        write(*,*)'dinit.f after map_d_hot'
        call rho_ini_hot_nperp_roots(r,z,phi,cnpar1)  
c       stop 'after map_d_hot'
c--------------------------------------------------------------
c       calculate all roots N_perpendicular of the hot plasma
c       dispersion function D_hot(N_perp=0) at the interval
c       0 < N_perpendicular < cN_perp_root_max
c----------------------------------------------------------
        call calculate_hot_nperp_roots(z,r,phi,cnpar1,
     &  n_hot_roots,N_perp_root_ar)

c        stop 'dinit.f after calculate_hot_nperp_roots'
        iraystop=1
        return
      endif !   i_look_roots=1    
c-----------------------------------------------------------------
      if ((id.eq.1).or.(id.eq.2)) then
c-------calculate two cold plasma roots for different ioxm
        ioxm_loc=ioxm
        ioxm=+1
        call cninit(z,r,phi,cnpar1,cntheta_ini,cnphi_ini,
     1            cnz,cnr,cm,iraystop)

        if (iraystop.eq.0) then
           cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar1**2)   
        else
           cnper =-1.d0
        endif
        wn_perp_ioxm_p(ifreq_write)=cnper

        ioxm=-1
        call cninit(z,r,phi,cnpar1,cntheta_ini,cnphi_ini,
     1            cnz,cnr,cm,iraystop)

        if (iraystop.eq.0) then
           cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar1**2)   
        else
           cnper =-1.d0
        endif
        wn_perp_ioxm_m(ifreq_write)=cnper

        ioxm= ioxm_loc
        wye_0(ifreq_write)=y(z,r,phi,1)
        wxe_0(ifreq_write)=x(z,r,phi,1)
cend_test
      endif !id=1,2

c      write(*,*)'dinit_1ray before cninit ksi_nperp',ksi_nperp

c      call plot_cold_n_perp_omega_npar(z,r,phi,cnpar1,
c     &0.1d0,3.d0,1000)

cBH130506      call plot_cold_n_perp_omega_npar_iray(z,r,phi,cnpar1,
cBH130506      &0.1d0,3.d0,1000)
c     &0.5d0,40.d0,1000)

      call cninit(z,r,phi,cnpar1,cntheta_ini,cnphi_ini,
     1            cnz,cnr,cm,iraystop)

c      stop' dinit_1ray after cninit '

      if (i_plot_wave_normal_cold.eq.1) then
c------------------------------------------------------
c       plot wave normal surfaces at initial ray point 
c------------------------------------------------------
        n_gam=500 ! set the number of angle points at the interval (0,2*pi) 
                  ! at the wave normal surface plot       

        if (iraystop.eq.0) then
           cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar1**2)   
        else
           cnper =-1.d0
        endif

cyup        write(*,*)'dinit.f  before wave_normal_surface'
       
        ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
        call wave_normal_surface(z,r,phi,cnpar1,cnper,0.d0,n_gam,ibmx)

cyup        write(*,*)'dinit.f  before wave_ray_normal_surface'

      call wave_ray_normal_surface(z,r,phi,cnpar1,cnper,0.d0,n_gam,ibmx)
               
      endif

c--------------------------------------------------------
c      call cninit(z,r,phi,cnpar1,cnteta,cnphi,
c     1            cnz,cnr,cm,iraystop)

      if (iraystop.eq.1) then
         return
      end if

cyup      write(*,*)'dinit after cninit cn**2',cnz**2+cnr**2+(cm/r)**2,
cyup     &'cn=',dsqrt(cnz**2+cnr**2+(cm/r)**2)

      cn2p=cnz**2+cnr**2+cnphi**2
cSm030513
c      cn2=cnz**2+cnr**2+cnphi**2
      cnphi_loc=cm/r
      cn2=cnz**2+cnr**2+cnphi_loc**2
      cnper=dsqrt(cn2-cnpar1**2)
cyup      write(*,*)'in dinit cn2',cn2,'cnper,cnpar1',cnper,cnpar1

ctest_launch
c       z= 0.1364802434874849 !theta= 58.62661361694336
c       r= 0.9810885719351724
c       phi= 0.000000000000000E+00
c       cnz= -7.936947824097243E-02
c       cnr= 5.822046050918139E-02
c       cm= 0.5224828576538614
c       prmt(3)=-prmt(3) !to create the negative time
c       i_vgr_ini=+1
       
c      z= 0.1815012404341051 !bound point
c      r= 1.065852982228670
c      phi= 9.379259749792192E-02
c      cnz= 0.3146127910809676
c      cnr= 0.7224110148199956
c      cm= 0.5224828576538614

c       z= 0.1815526690186028 !bound point
c       r= 1.065977980682529 
c       phi= 9.386458441905626E-02
c       cnz= 0.3147511593634075
c       cnr= 0.7227255362319519
c       cm= 0.5224828576538614

c        z= 0.1815401084275075
c        r= 1.065986991882685
c        phi= 8.866550623716300E-02
c        cnz= 0.3147767582210387
c        cnr= 0.7224196375916611
c        cm= 0.5226713433206168
        
cend_test

      bmod=b(z,r,phi)
      gam=gamma1(z,r,phi,cnz,cnr,cm)
cyup      write(*,*)'1ray before d=hamilt1'
cyup      write(*,*)'z,r,phi before d=hamilt1 z,r,phi',z,r,phi
cyup      write(*,*)'z,r,phi before d=hamilt1 cnz,cnr,cm',cnz,cnr,cm

      dh=hamilt1(z,r,phi,cnz,cnr,cm)

cyup      write(*,*)'dinit_1ray after d=hamilt1 dh=',dh
c The check of the Hamiltonian value for the found initial conditions.
      epshamin=1.d-6
      epshamin=1.d-2
      epshamin=1.d-1
     
c      if (dabs(dh).gt.epshamin) then
c         write(*,*)'in dinit_1ray initial hamiltonian dh>epsham'
c	 write(*,*)'in dinit1ray dh,epshamin',dh,epshamin
c	 write(*,*)'in dinit1ray correct epshamin or initial conditions'
c         iraystop=1
c         return
c      endif
c-------------------------------------------------------------------

cyup      write(*,*)'before outinit'
      call outinit(u)
cyup      write(*,*)'in dinit_1ray after call outini nrayelt= ',nrayelt
      irefl=0
c      cnz=-0.5166646619861239d0 
c      cnr=-0.7701078676733518d0 
      u(1)=z
      u(2)=r
      u(3)=phi
      u(4)=cnz
      u(5)=cnr
      u(6)=cm
      
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'dinit_1ray before prep3d powini',powini
      write(*,*)'dinit_1ray before prep3d u',u
      write(*,*)'dinit_1ray before prep3d xma,yma',xma,yma
      endif ! outprint
  
      call prep3d(0.0,u,deru,iraystop)

      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'initial data'
      write(*,*)'z=',z,'r=',r,'phi=',phi
      write(*,*)'cnz=',cnz,'cnr=',cnr,'cm=',cm
      write(*,*)'rho=',rho
      write(*,*)'x_e=',x(z,r,phi,1),'y_e=',y(z,r,phi,1)
      endif ! outprint
c-----------------------------------------------------------
c     check the group velocity
c-----------------------------------------------
      bmod=b(z,r,phi)        

      call rside1(0.d0,u,deru)

      cnpar1=(cnz*bz+cnr*br+(cm/r)*bphi)/bmod
   
      v_gr=dsqrt(deru(1)**2+deru(2)**2+deru(3)**2)
      v_gr_rho=(dpdzd*deru(1)+dpdrd*deru(2))/
     &              dsqrt(dpdzd**2+dpdrd**2)

cyup      write(*,*)'v_gr,v_gr_rho',v_gr,v_gr_rho

c-----the angle between grad_psi and n_perp
c     vector_n_perp=vector_n-vector_n_par
c     vector_n_par=cnpar*vector_b/mod(b)
 
      cnpar_r=cnpar1*br/bmod
      cnpar_z=cnpar1*bz/bmod
      cnpar_phi=cnpar1*bphi/bmod

c      write(*,*)'cnpar1,cnpar_z,cnpar_r,cnpar_phi',
c     &           cnpar1,cnpar_z,cnpar_r,cnpar_phi

      cnper_r=cnr-cnpar_r
      cnper_z=cnz-cnpar_z
      cnper_phi=cm/r-cnpar_phi

c      write(*,*)'cnper_z,cnper_r,cnper_phi',
c     &           cnper_z,cnper_r,cnper_phi

c      write(*,*)'dpdzd,dpdrd,dsqrt(dpdrd*2+dpdzd**2)',
c     &           dpdzd,dpdrd,dsqrt(dpdrd*2+dpdzd**2)

      cnper_test=dsqrt(cnper_r**2+cnper_z**2+cnper_phi**2)


c      write(*,*)'dinit================'
c      write(*,*)'cnper_z,cnper_r',cnper_z,cnper_r
c      write(*,*)'dpdzd,dpdrd',dpdzd,dpdrd
c      write(*,*)'dsqrt(dpdzd**2+dpdrd**2)',
c     &           dsqrt(dpdzd**2+dpdrd**2)
c      write(*,*)'cnper_test',cnper_test

      cos_ksi_test=(cnper_z*dpdzd+cnper_r*dpdrd)/
     &                      (dsqrt(dpdzd**2+dpdrd**2)*cnper_test)
   
      cos_ksi_vg=v_gr_rho/v_gr
      
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'cos_ksi_test,dcos(ksi_nperp*pi/180.d0)',
     &           cos_ksi_test,dcos(ksi_nperp*pi/180.d0)
      write(*,*)'cos_ksi_vg',cos_ksi_vg
      write(*,*)'*******ksi_vg is the angle between group velocity'
      write(*,*)'and grad(psi) [degree].'
      write(*,*)'If angle less 90 degree then'
      write(*,*)'the group velocity is directed outside the plasma.'
      write(*,*)'ksi_vg=',dacos(cos_ksi_vg)*180.d0/pi

cRobtAndre140421      write(*,*)'ksi_nperp,dacos(cos_ksi_test)',
cRobtAndre140421     &           ksi_nperp,dacos(cos_ksi_test)
      write(*,*)'ksi_nperp,dacos(cos_ksi_test)',
     &           ksi_nperp,dacos(min(1.d0,cos_ksi_test))
      endif ! outprint

c-----safety factor calculations
      psi_initial=psi_rho(rho)
      q_initial=qsafety_psi(psi_initial)
      b_av=b_average(psi_initial)
      
      if(outprint.eq.'enabled')then !YuP[2018-01-17] Added
      write(*,*)'psi_initial,q_initial',psi_initial,q_initial
      write(*,*)'b(z,r,phi)',b(z,r,phi)
      write(*,*)'b_average(psi_initial)',b_av
      write(*,*)'end dinit_1ray'
      endif ! outprint
      
      return
      end




c        **********************zeffcalc************************
c        *                        -                           *
c        * this subroutine calculates zeff(ndens)             *
c        * using the ions densities 			      *
c        ******************************************************
c
c------------------------------------------------------------------
c        output parameters: array zeff(ndens)
c------------------------------------------------------------------

      subroutine zeffcalc
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'six.i'
c---------------------------------------------------------
c     calculation of the table for the radial profile zeff1(ndens)
c---------------------------------------------------------
      if (nbulk.eq.1) then
         write(*,*)'Warning dinit: nbulk=1, it will be created zeff=1'
      endif

      do j=1,ndens
         if (nbulk.eq.1) then
            zeff1(j)=1.d0
            goto 10
         endif
c---------------------------------------------------------
c        electron density : dens_el
         dens_el=dens1(j,1)
         zeff1(j)=0.d0
	 do i=2,nbulk
           dens_ion=dens1(j,i)
           if (dens_el.ne.0.d0) then
              denside=dens_ion/dens_el
           else
              if (dens_ion.eq.0.d0) then
                 denside=1.d0
              else
                 write(*,*)'in zeffcalc j',j
	         write(*,*)'densel=0 but dens_ion(i).ne.0'
	         write(*,*)'i number of plasma component=',i
	     write(*,*)'change the parameters of the density profiles'
	         stop
	      endif
	   endif

	   zeff1(j)=zeff1(j)+charge(i)*charge(i)*denside
         enddo !i=2,nbulk
 10      continue
c---------------------------------------------------------
      enddo !j=1,ndens
c---------------------------------------------------------
c     end of calculation of the table for the radial profile
      return
      end
c        **********************zeffcal1************************
c        *                        -                           *
c        * this subroutine calculates zeff(ndens)             *
c        * using the ions densities. 			      *
c        *It is the old version FOR the charge neutrality control
c        ******************************************************
c
c------------------------------------------------------------------
c        output parameters: array zeff(ndens)
c------------------------------------------------------------------

      subroutine zeffcal1
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'six.i'
c---------------------------------------------------------
c     calculation of the table for the radial profile zeff1(ndens)
c---------------------------------------------------------
      h=1.d0/(ndens-1.d0)
      if (nbulk.eq.1) then
         write(*,*)'Worning nbulk=1, it will be created zeff=1'
      endif

      do j=1,ndens
         rho=h*(j-1)
         if (nbulk.eq.1) then
            zeff1(j)=1.d0
            goto 10
         endif
c---------------------------------------------------------
c        electron density : dens_el
	 if (idens.eq.0) then
c           analytical representation of the electron density profile
	    dens_el=(dense0(1)-denseb(1))*
     1              (1-rho**rn1de(1))**rn2de(1)+denseb(1)
	 else
c           table representation of the electron density profile
	    dens_el=dens1(j,1)
	 endif
c---------------------------------------------------------
         zeff1(j)=0.d0
	 qusneutr=0.d0 ! for the control of the charge-neutrality
         do i=2,nbulk
	    if (idens.eq.0) then
c              analytical representation of the ion density profiles
	       dens_ion=(dense0(i)-denseb(i))*
     1                  (1-rho**rn1de(i))**rn2de(i)+denseb(i)
	    else
c              table representation of the ion density profiles
	       dens_ion=dens1(j,i)
	    endif

            if (dens_el.ne.0.d0) then
	       denside=dens_ion/dens_el
	    else
	       if (dens_ion.eq.0.d0) then
	         denside=1.d0
	       else
	         write(*,*)'in zeffcalc rho=',rho
		 write(*,*)'densel=0 but dens_ion(i).ne.0'
		 write(*,*)'i number of plasma component=',i
	     write(*,*)'change the parameters of the density profiles'
		 stop
	       endif
	    endif

c---------------------------------------------------------
c   control of the charge-neutrality condition in the point rho(i) 
c---------------------------------------------------------
	    qusneutr=qusneutr+charge(i)*denside
	    zeff1(j)=zeff1(j)+charge(i)*charge(i)*denside
	 enddo !i=2,nbulk
c---------------------------------------------------------
         if (qusneutr.ne.1.d0) then
            write(*,*)'in zeffcalc the bad quasi-neutrality'
	    write(*,*)'in rho(j)=',rho,'j=',j
            write(*,*)'sum({i=2,nbulk}(charge(i)*dens(i(i)/dens_e)',
     1	    qusneutr
            write(*,*)'change the charge(j) or the densities profiles'
    	    stop
         endif
 10      continue
      enddo !j=1,ndens
c---------------------------------------------------------
c     end of calculation of the table for the radial profile
      return
      end


c     **********************denscalc************************
c     *                        -                           *
c     * this subroutine calculates ions densities          *
c     * dens1(j,nbulk) and dens1(j,nbulk-1)                *
c     * (here j=1,ndens),				   *
c     * using the zeff and dens1(ndens,i) i=2,nbulk-2      *
c     ******************************************************
c
c------------------------------------------------------------------
c     output parameters: array dens1(ndens,nbulk)
c------------------------------------------------------------------

      subroutine denscalc
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'six.i'
c---------------------------------------------------------
c     calculation of the table for the radial profile
c     dens1(ndens,nbulk), dens1(ndens,nbulk-1)
c---------------------------------------------------------
cc      write(*,*)'in denscalc ndens,nbulk,idens,izeff'
cc      write(*,*)ndens,nbulk,idens,izeff
      if (nbulk.le.2) then
         write(*,*)'in denscalc nbulk.le.2'
         write(*,*)'zeff must be equal charge(2)'
         write(*,*)'and dense(2)=dense(1)/charge(2)'
         write(*,*)'use the option izeff=0 and these parameres'
	 stop
      else
c        nbulk.ge.3
         if( charge(nbulk).eq.charge(nbulk-1)) then
	    write(*,*)'Warning in denscalc: nbulk(.ge.3)=',nbulk
	    write(*,*)'in denscalc: charge(nbulk)=charge(nbulk-1)'
	    write(*,*)'it is impossible to find the ions denscalc'
	    write(*,*)'change charge(nulk) or charge(nbulk-1)'
	    write(*,*)'it should be charge(nulk)>charge(nbulk-1)'
	    write(*,*)'or use the option izeff=0'
	    stop
	 endif
      endif

      h=1.d0/(ndens-1.d0)
      do j=1,ndens
         rho=h*(j-1)
cSAP091015
         write(*,*)'in denscalc j,rho',j,rho
c---------------------------------------------------------
c        electron density : dens_el
	 if (idens.eq.0) then
c           analytical representation of the electron density profile
	    dens_el=(dense0(1)-denseb(1))*
     1              (1-rho**rn1de(1))**rn2de(1)+denseb(1)
	 else
c           table representation of the electron density profile
	    dens_el=dens1(j,1)
	 endif
c---------------------------------------------------------
         sum1=0.d0
         sum2=0.d0

cSAP091015
	 write(*,*)'in denscalc j=',j,'rho',rho,'dens_el',dens_el
         write(*,*)'nbulk',nbulk

	 if(nbulk.ge.4) then
c           sum1=sum{i=2,nbulk-2}(density(i)/electron_density*
c                charge(i)*(charge(nbulk)-charge(i))
c           sum2=sum{i=2,nbulk-2}(density(i)/electron_density*
c                charge(i)*(charge(nbulk-1)-charge(i))
            do i=2,nbulk-2
	       if ((charge(nbulk)-charge(i)).lt.0) then
	         write(*,*)'in denscalc charge(nbulk).lt.charge(i)'
		 write(*,*)'change array charge in genray.in'
		 write(*,*)'in array charge(i+1) must be >charge(i)'
		 write(*,*)'change the array charge(i) in genray.in'
		 stop
	       endif
	       if ((charge(nbulk-1)-charge(i)).lt.0) then
	         write(*,*)'in denscalc charge(nbulk-1).lt.charge(i)'
		 write(*,*)'change array charge in genray.in'
		 write(*,*)'in array charge(i+1) must be >charge(i)'
		 write(*,*)'change the array charge(i) in genray.in'
		 stop
	       endif
	       if (idens.eq.0) then
c                 analytical representation of the ion density profiles
	          dens_ion=(dense0(i)-denseb(i))*
     1                     (1-rho**rn1de(i))**rn2de(i)+denseb(i)
	       else
c                 table representation of the ion density profiles
	          dens_ion=dens1(j,i)
	       endif

               if (dens_el.ne.0.d0) then
	          denside=dens_ion/dens_el
	       else
	          if(dens_ion.eq.0.d0) then
	            denside=1.d0
		  else
		    write(*,*)'in denscalc dens_el=0,
     1		    but dens_ion.ne.0, j=',j,'i=',i
                    write(*,*)'change density profiles parameters'
		    stop
	          endif
	       endif
	       sum1=sum1+denside*charge(i)*(charge(nbulk)-charge(i))
	       sum2=sum2+denside*charge(i)*(charge(nbulk-1)-charge(i))
	    enddo !nbulk
	 endif ! nbulk.ge.4

cSAP091015
        write(*,*)'charge(nbulk),charge(nbulk-1)',
     &             charge(nbulk),charge(nbulk-1)

	 p3=charge(nbulk)-charge(nbulk-1)
	 if (p3.lt.0.d0) then
	    write(*,*)'in denscalc charge(nbulk).lt.charge(nbulk-1)'
	    write(*,*)'change array charge in genray.in'
	    write(*,*)'in array charge(nbulk) must be >charge(nbulk-1)'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif

cSAP091015
         write(*,*)'zeff1(j),sum1',zeff1(j),sum1
	
         p1=charge(nbulk)-zeff1(j)-sum1
	 if (p1.lt.0.d0) then
	    write(*,*)'in denscalc charge(nbulk)-zeff1(j)-sum1.lt.0'
	    write(*,*)'j=',j,'it is impossible to calculate
     1	    the density(j,nbulk-1) using zeff profile'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
	 p2=zeff1(j)-charge(nbulk-1)+sum2
	 if (p2.lt.0.d0) then
	    write(*,*)'in denscalc charge(nbulk)-zeff1(j)-sum1.lt.0'
	    write(*,*)'j=',j,'it is impossible to calculate
     1	    the density(j,nbulk) using zeff profile'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
	 dens1(j,nbulk-1)=p1/(charge(nbulk-1)*p3)*dens1(j,1)
	 dens1(j,nbulk)=p2/(charge(nbulk)*p3)*dens1(j,1)
cSAP091015
         write(*,*)'in denscalc dens1(j,nbulk-1),dens1(j,nbulk)'
         write(*,*)dens1(j,nbulk-1),dens1(j,nbulk)
      enddo ! j
 10   continue
      return
      end
c     **********************denscalp************************
c     * this subroutine calculates electron and		   *
c     * ions densities profiles: dense(j,1),               *
c     * dens1(j,nbulk) and dens1(j,nbulk-1)                *
c     * (here j=1,ndens),				   *
c     * using the zeff and tempe1(ndens,i) i=2,nbulk       *
c     * and eqdsk pres (pressure)
c     ******************************************************
c
c------------------------------------------------------------------
c     output parameters: array dens1(ndens,nbulk)
c------------------------------------------------------------------

      subroutine denscalp
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'six.i'
c---------------------------------------------------------
c     calculation of the table for the radial profile
c     dens1(ndens,nbulk), dens1(ndens,nbulk-1)
c---------------------------------------------------------
cc      write(*,*)'in denscalp ndens,nbulk,idens,izeff'
cc      write(*,*)ndens,nbulk,idens,izeff

      if (nbulk.eq.2) then
         write(*,*)'in denscalp nbulk.eq.2'
         write(*,*)'use another option izeff or another parameters'
	 stop
      else
        if (nbulk.ge.3) then
c        nbulk.ge.3
         if( charge(nbulk).eq.charge(nbulk-1)) then
	    write(*,*)'Warning in denscalc: nbulk(.ge.3)=',nbulk
	    write(*,*)'in denscalp: charge(nbulk)=charge(nbulk-1)'
	    write(*,*)'it is impossible to find the ions denscalp'
	    write(*,*)'change charge(nulk) or charge(nbulk-1)'
	    write(*,*)'it should be charge(nulk)>charge(nbulk-1)'
	    write(*,*)'or use another option izeff'
	    stop
	 endif
        endif
      endif

      h=1.d0/(ndens-1.d0)
      do j=1,ndens
         rho=h*(j-1)
	 psi=psi_rho(rho)
	 pressure=prespsi(psi)/1.6d3
cc         write(*,*)'in denscalp j,rho',j,rho

cSAP090801
         if(nbulk.eq.1) then
c----------electron density n=p/2T)
           dens1(j,1)=pressure/(2.d0*temp1(j,1))
           goto 10
         endif
c---------------------------------------------------------
         sum1=0.d0
         sum2=0.d0
         sum4=0.d0
cc	 write(*,*)'in denscalc j=',j,'rho',rho,'dens_el',dens_el
	 if(nbulk.ge.4) then
c           sum1=sum{i=2,nbulk-2}(density(i)*charge(i)*
c                                 (charge(nbulk)-charge(i))
c           sum2=sum{i=2,nbulk-2}(density(i)*charge(i)*
c                                 (charge(nbulk-1)-charge(i))
c           sum4=sum{i=2,nbulk-2}(density(i)*temp(i))
            do i=2,nbulk-2
	       if ((charge(nbulk)-charge(i)).lt.0) then
	         write(*,*)'in denscalp charge(nbulk).lt.charge(i)'
		 write(*,*)'change array charge in genray.in'
		 write(*,*)'in array charge(i+1) must be >charge(i)'
		 write(*,*)'change the array charge(i) in genray.in'
		 stop
	       endif
	       if ((charge(nbulk-1)-charge(i)).lt.0) then
	         write(*,*)'in denscalp charge(nbulk-1).lt.charge(i)'
		 write(*,*)'change array charge in genray.in'
		 write(*,*)'in array charge(i+1) must be >charge(i)'
		 write(*,*)'change the array charge(i) in genray.in'
		 stop
	       endif
	       if (idens.eq.0) then
c                 analytical representation of the ion density profiles
	          dens_ion=(dense0(i)-denseb(i))*
     1                     (1-rho**rn1de(i))**rn2de(i)+denseb(i)
	       else
c                 table representation of the ion density profiles
	          dens_ion=dens1(j,i)
	       endif

	       sum1=sum1+dens_ion*charge(i)*(charge(nbulk)-charge(i))
	       sum2=sum2+dens_ion*charge(i)*(charge(nbulk-1)-charge(i))
	       sum4=sum4+dens_ion*temp1(j,i)
	    enddo !nbulk
	 endif ! nbulk.ge.4
	 p3=charge(nbulk)-charge(nbulk-1)
	 if (p3.lt.0.d0) then
	    write(*,*)'in denscalp charge(nbulk).lt.charge(nbulk-1)'
	    write(*,*)'change array charge in genray.in'
	    write(*,*)'in array charge(nbulk) must be >charge(nbulk-1)'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
c----------------------------------------------------------------
c the formula for dense(nbulk-1) and dense(nbulk)
c dense(nbulk-1)=((charge(nbulk)-zeff)*dense(1)-sum1))/
c                (charge(nbulk-1)*(charge(nbulk)-charge(nbulk-1)))
c dense(nbulk)  =((zeff-charge(nbulk-1))*dense(1)+sum2))/
c                (charge(nbulk)*(charge(nbulk)-charge(nbulk-1)))
c----------------------------------------------------------------
c calculation of the electron density dens1(j,1)
c from the the plasma pressure:
c pressure=dens1(j,1)*temp1(j,1)+sum(i=2,nbulk1)(dens1(j,i)*temp1(j,i))
c    +dens1(j,nbulk-1)*temp1(j,nbulk-1)+dens1(j,nbulk)*temp1(j,nbulk)
         p4=temp1(j,1)+
     1   temp1(j,nbulk-1)*(charge(nbulk)-zeff1(j))/(charge(nbulk-1)*p3)
     1   +temp1(j,nbulk)*(zeff1(j)-charge(nbulk-1))/(charge(nbulk)*p3)
         p5=pressure-sum4+
     1   temp1(j,nbulk-1)*sum1/(charge(nbulk-1)*p3)-
     1   temp1(j,nbulk)*sum2/(charge(nbulk)*p3)
         if(p4.eq.0.d0) then
	   write(*,*)'in denscalp:izeff=4 p4=0 bad conditions j=',j
	   stop
	 else
	   dens1(j,1)=p5/p4
	   if(dens1(j,1).lt.0.d0) then
	     write(*,*)'in denscalp:izeff=4, des1(j,1)<0,j=',j
	     write(*,*)'bad conditions'
	     stop
	   endif
	 endif
         p1=dens1(j,1)*(charge(nbulk)-zeff1(j))-sum1
	 if (p1.lt.0.d0) then
	    write(*,*)'in denscalp j=',j
	    write(*,*)'dens1(j,1)*(charge(nbulk)-zeff1(j))-sum1.lt.0'
	    write(*,*)'j=',j,'it is impossible to calculate
     1	    the density(j,nbulk-1) using zeff profile'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
	 p2=dens1(j,1)*(zeff1(j)-charge(nbulk-1))+sum2
	 if (p2.lt.0.d0) then
	    write(*,*)'in denscalp j=',j
	    write(*,*)'dens1(j,1)*(charge(nbulk)-zeff1(j))-sum1.lt.0'
	    write(*,*)'j=',j,'it is impossible to calculate
     1	    the density(j,nbulk) using zeff profile'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
	 dens1(j,nbulk-1)=p1/(charge(nbulk-1)*p3)
	 dens1(j,nbulk)=p2/(charge(nbulk)*p3)
         write(*,*)'in denscalp: j=',j
	 write(*,*)'dens1(j,1),dens1(j,nbulk-1),dens1(j,nbulk)'
	 write(*,*)dens1(j,1),dens1(j,nbulk-1),dens1(j,nbulk)

 10      continue

      enddo ! j
      !pause !!!
      return
      end


c-----following functions are for the root (N_perp) determinations

      double complex function chammuz(cnper)
c-----calculates the complex Mazzucato hamiltonian for given Re(N_perp)
c     the part of the input data are in common
c     common/ham_h/mode,xe,ye,te_kev,tpope,vflowe,
c     .cnpar,cnprim,ihermloc,iham_loc
      implicit none 
c-----input
      double precision cnper
      double precision mode,xe,ye,te_kev,tpope,vflowe,cnpar,cnprim
      integer ihermloc,iham_loc
      common/ham_h/mode,xe,ye,te_kev,tpope,vflowe,
     .cnpar,cnprim,ihermloc,iham_loc
      external complx1
c-----locals
      double precision mu,nz,np2
      double complex cpp,sol,chamilmz
      write(*,*)' dinit chammuz mode,xe,ye,te_kev,cnpar,cnprim,
     *ihermloc,iham_loc'
      write(*,*)mode,xe,ye,te_kev,cnpar,cnprim,ihermloc,iham_loc
      mu=512.44d00/te_kev
      nz=cnpar
      cpp=dcmplx(cnper,0.0d0)
      sol=cpp*cpp
      np2=cnper*cnper
      write(*,*)'dinit in chammuz cnper',cnper
      write(*,*)'dinit in chammuz mode,xe,ye,mu,nz,np2,sol,ihermloc'
      write(*,*)mode,xe,ye,mu,nz,np2,sol,ihermloc
      call complx1(mode,xe,ye,mu,nz,np2,sol,ihermloc,chamilmz)
      chammuz=chamilmz
      write(*,*)'dinit in chammuz=',chammuz,'chamilmz',chamilmz

      return
      end             


      double precision function	rhammuz (cnper)
c-----calculates the real part Mazzucato hamiltonian for given Re(N_perp)
c     the part of the input data are in common /ham_h/
      implicit none 
c-----input	      
      double precision cnper
      double precision mode,xe,ye,te_kev,tpope,vflowe,cnpar,cnprim
      integer ihermloc,iham_loc
      common/ham_h/mode,xe,ye,te_kev,tpope,vflowe,
     .cnpar,cnprim,ihermloc,iham_loc
      double complex chammuz
      external chammuz
      write(*,*)'dinit rhammuz cnper',cnper
      rhammuz=dreal(chammuz(cnper))
      return
      end






c--------------mapxyb---------------------------------
      subroutine mapxyb(nrhomap,nthetmap)
   
c     (1) It creates the data for plots
c     X_e, Y_e B_tot, B_ptor, B_pol (rho,theta)
c     the results are in the output files
c     xybrhoth.dat: rho(i),theta(j),xe,ye,(xe+ye*ye),bmod,bphi,
c     *              dsqrt(bz**2+br**2)
c     Here (0<rho<1,0<theta<pi) are the points of mesh

c      implicit none
      implicit double precision (a-h,o-z)
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'

c-----input
      integer nrhomap,nthetmap  ! The number of mesh points (rho,theta) 
     
      double precision psi_rho,b,x,y
      external psi_rho  ! calculates psi(rho)
      external zr_psith ! calculates (z,r) at given (psi,theta_poloidal)
      external b,x,y
      
c-----local
      double precision hthetmap,hrhomap,
     * thetmap,rhomap,psix,
     * rmap,zmap,phi,xe,ye
      integer i,j

      if(myrank.ne.0) return
      
      write(*,*)'dinit begin mapxyb nrhomap,nthetmap',nrhomap,nthetmap

c-----steps of  mesh      
      hrhomap=(1.d0-1.d-4)/(nrhomap-1)
      write(*,*)'mapxyb hrhomap',hrhomap
        
      pi=4.d0*datan(1.d0) 
      hthetmap=pi/(nthetmap-1)
      write(*,*)'mapxyb thetmapn',hthetmap

      open(1,file='xybrhoth.dat')
     
1     format (8(' ',1pe11.4))      
      write(1,*)' rho theta xe ye uh bmod bphi bpol '
      phi=0.d0
      do i=1,nrhomap
        rhomap=hrhomap*(i-1)
        write(*,*)'i,rhomap',i,rhomap
        psix=psi_rho(rhomap)

        do j=1,nthetmap
           thetmap=hthetmap*(j-1)
           write(*,*)'j,thetmap',j,thetmap
           call zr_psith(psix,thetmap,zmap,rmap)
           bmod=b(zmap,rmap,phi)
           write(*,*)'bmod',bmod
           xe=x(zmap,rmap,phi,1)
           ye=y(zmap,rmap,phi,1)
           write(1,1)rhomap,thetmap,xe,ye,(xe+ye*ye),bmod,bphi,
     *              dsqrt(bz**2+br**2)
        enddo
      enddo
                    
      close(1)

      return
      end

      subroutine set_freq(ifreq)
c-----calculate the wave frequency frqncy,w0,v0,v(nbulk),w(nbulk)
c     and put it to common/one.i/
c     It used for the emission calculations with nfreq.ne.1
c     If (nfreq.eq.1) it will not do anything

      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'write.i'
c-----input
      integer ifreq !is the number of frequencies in wfreq(nfreqa)

      if(nfreq.eq.1) goto 10

      frqncy=wfreq(ifreq)
      v0=806.2/(frqncy*frqncy)
      w0=28.0*b0/frqncy 
      v(1)=v0
      w(1)=w0

      do  i=2,nbulk
         v(i)=v0*charge(i)**2/dmas(i)
         w(i)=w0*charge(i)/dmas(i)
cyup         write(*,*)'dinit.f set_freq: i charge(i),dmas(i),v(i),w(i)'
cyup     +   ,i,charge(i),dmas(i),v(i),w(i)
      enddo

 10   continue
cyup      write(*,*)'in set_freq ifreq,frqncy',ifreq,frqncy

      return
      end





      subroutine poloidal_vector(z,r,phi,poloid_z,poloid_r)
c--------------------------------------------------------
c     calculate the unit poloidal vector coordinates (poloid_z,poloid_r)
c     directed clockwise  
c--------------------------------------------------------------
c                                         |e_z         e_r         e_phi | 
c     poloidal vector =[grad_psi * e_phi]=|grad_psi_z  grad_psi_r    0   |=
c                                         |0           0             1   |
c
c     =e_z*grad_psi_r - e_r*grad_psi_z
c
c     popoidal_z =   grad_psi_r/dsqrt(grad_psi_z**2+grad_psi_r**2)
c     popoidal_z = - grad_psi_z/dsqrt(grad_psi_z**2+grad_psi_r**2)
c
c-------------------------------------------------------------------------
      implicit none  
      include 'param.i'
      include 'one.i'
c-----input
      real*8 z,r,phi
c-----output
      real*8 poloid_z,poloid_r !poloidal vector coordinates 
c-----externals
      real*8 b
c-----locals
      real*8 grad_psid

      bmod=b(z,r,phi)
      grad_psid=dsqrt(dpdzd**2+dpdrd**2)
      poloid_z=  dpdrd/grad_psid
      poloid_r= -dpdzd/grad_psid

      return
      end


      subroutine ainalloc_onetwo_no_nml_i
c-----allocate pointers in onetwo_no_nml.i
      implicit none
      include 'param.i'
      include 'one.i'
      include 'onetwo.i'
      real*8  zero
      integer istat

      zero=0.d0

      allocate(  powtot_s(1:nbulk),STAT=istat)
      call bcast(powtot_s,zero,SIZE(powtot_s))

      allocate(  power(1:NR),STAT=istat)
      call bcast(power,zero,SIZE(power))

      allocate(  current(1:NR),STAT=istat)
      call bcast(current,zero,SIZE(current))

      allocate(  temparr(1:NR),STAT=istat)
      call bcast(temparr,zero,SIZE( temparr))

      allocate(  zeffarr(1:NR),STAT=istat)
      call bcast(zeffarr,zero,SIZE(zeffarr))

      allocate(  spower(1:NR),STAT=istat)
      call bcast(spower,zero,SIZE(spower))

      allocate(  scurrent(1:NR),STAT=istat)
      call bcast(scurrent,zero,SIZE(scurrent))

      allocate(  powden(1:NR),STAT=istat)
      call bcast(powden,zero,SIZE(powden))

      allocate(  currden(1:NR),STAT=istat)
      call bcast(currden,zero,SIZE(currden))

      allocate(  power_e(1:NR),STAT=istat)
      call bcast(power_e,zero,SIZE(power_e))

      allocate(  power_i(1:NR),STAT=istat)
      call bcast(power_i,zero,SIZE(power_i))

      allocate(  power_s(1:NR,1:nbulk),STAT=istat)
      call bcast(power_s,zero,SIZE(power_s))

      allocate(  power_cl(1:NR),STAT=istat)
      call bcast(power_cl,zero,SIZE(power_cl))

      allocate(  spower_e(1:NR),STAT=istat)
      call bcast(spower_e,zero,SIZE(spower_e))

      allocate(  spower_i(1:NR),STAT=istat)
      call bcast(spower_i,zero,SIZE(spower_i))

      allocate(  spower_s(1:NR,1:nbulk),STAT=istat)
      call bcast(spower_s,zero,SIZE(spower_s))

      allocate(  spower_cl(1:NR),STAT=istat)
      call bcast(spower_cl,zero,SIZE(spower_cl))

      allocate(  powden_e(1:NR),STAT=istat)
      call bcast(powden_e,zero,SIZE(powden_e))

      allocate( powden_i(1:NR),STAT=istat)
      call bcast(powden_i,zero,SIZE(powden_i))

      allocate(  powden_s(1:NR,1:nbulk),STAT=istat)
      call bcast(powden_s,zero,SIZE( powden_s))

      allocate( powden_cl(1:NR),STAT=istat)
      call bcast(powden_cl,zero,SIZE(powden_cl))

      allocate(  currden_s(1:NR),STAT=istat)
      call bcast( currden_s,zero,SIZE( currden_s))

      allocate(  powden_e_s(1:NR),STAT=istat)
      call bcast(powden_e_s,zero,SIZE(powden_e_s))

      allocate(  powden_i_s(1:NR),STAT=istat)
      call bcast(powden_i_s,zero,SIZE(powden_i_s))

      allocate(  powden_cl_s(1:NR),STAT=istat)
      call bcast(powden_cl_s,zero,SIZE(powden_cl_s))

      allocate(  densprof(1:NR,1:nbulk),STAT=istat)
      call bcast(densprof,zero,SIZE(densprof))

      allocate(  temprof(1:NR,1:nbulk),STAT=istat)
      call bcast(temprof,zero,SIZE(temprof))

      allocate(  zefprof(1:NR),STAT=istat)
      call bcast(zefprof,zero,SIZE(zefprof))

      allocate(  rho_bin(1:NR),STAT=istat)
      call bcast(rho_bin,zero,SIZE(rho_bin))

      allocate(  rho_bin_center(1:NR-1),STAT=istat)
      call bcast(rho_bin_center,zero,SIZE(rho_bin_center))

      allocate(  binvol(1:NR-1),STAT=istat)
      call bcast(binvol,zero,SIZE(binvol))

      allocate(   binarea(1:NR-1),STAT=istat)
      call bcast( binarea,zero,SIZE( binarea))

      allocate( binarea_pol(1:NR-1),STAT=istat)
      call bcast(binarea_pol,zero,SIZE(binarea_pol))

      allocate( pollen(1:NR-1),STAT=istat)
      call bcast(pollen,zero,SIZE( pollen))

      allocate(  cur_den_parallel(1:NR),STAT=istat)
      call bcast(cur_den_parallel,zero,SIZE(cur_den_parallel))

      allocate( s_cur_den_parallel(1:NR-1),STAT=istat)
      call bcast(s_cur_den_parallel,zero,SIZE(s_cur_den_parallel))

      allocate( s_cur_den_onetwo(1:NR-1),STAT=istat)
      call bcast(s_cur_den_onetwo,zero,SIZE(s_cur_den_onetwo))

      allocate(  s_cur_den_toroidal(1:NR),STAT=istat)
      call bcast( s_cur_den_toroidal,zero,SIZE( s_cur_den_toroidal))

      allocate(  s_cur_den_poloidal(1:NR),STAT=istat)
      call bcast( s_cur_den_poloidal,zero,SIZE( s_cur_den_poloidal))

      return
      end

      subroutine ainalloc_emissa_no_nml_i
c-----allocate pointers in emissa_no_nml.i
      implicit none
      include 'param.i'
      include 'one.i'
      include 'emissa.i'
      
      real*8  zero
      integer istat

      zero=0.d0

cyup      write(*,*)'in ainalloc_emissa_no_nml_i nrelta4)',nrelta4

      allocate( cx_z(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate cx_z istat',istat
      call bcast(cx_z,zero,SIZE(cx_z))

      allocate( tr_z(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate tr_z istat',istat
      call bcast(tr_z,zero,SIZE(tr_z))

      allocate( cx_r(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate cx_r istat',istat
      call bcast(cx_r,zero,SIZE(cx_r))

      allocate( tr_r(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate tr_r istat',istat
      call bcast(tr_r,zero,SIZE(tr_r))

      allocate( cx_phi(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate cx_phi istat',istat
      call bcast(cx_phi,zero,SIZE(cx_phi))

      allocate( tr_phi(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate tr_phi istat',istat
      call bcast(tr_phi,zero,SIZE(tr_phi))

      allocate( cx_cnz(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate cx_cnz istat',istat
      call bcast(cx_cnz,zero,SIZE(cx_cnz))

      allocate( tr_cnz(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  tr_cnzc istat',istat
      call bcast(tr_cnz,zero,SIZE(tr_cnz))

      allocate( cx_cnr(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate   cx_cnr istat',istat
      call bcast(cx_cnr,zero,SIZE(cx_cnr))

      allocate( tr_cnr(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate    tr_cnr istat',istat
      call bcast(tr_cnr,zero,SIZE(tr_cnr))

      allocate( cx_cm(1:nrelta4),STAT=istat)
      write(*,*)'after allocate   cx_cm istat',istat
      call bcast(cx_cm,zero,SIZE(cx_cm))

      allocate( tr_cm(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate   tr_cm istat',istat
      call bcast(tr_cm,zero,SIZE(tr_cm))

c-----for ebw-x mode ray part 

      allocate( cx_z_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  cx_z_x  istat',istat
      call bcast(cx_z_x,zero,SIZE(cx_z_x))

      allocate( tr_z_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate tr_z_x  istat',istat
      call bcast(tr_z_x,zero,SIZE(tr_z_x))

      allocate( cx_r_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate cx_r_x  istat',istat
      call bcast(cx_r_x,zero,SIZE(cx_r_x))

      allocate( tr_r_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  tr_r_x istat',istat
      call bcast(tr_r_x,zero,SIZE(tr_r_x))

      allocate( cx_phi_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate cx_phi_x istat',istat
      call bcast(cx_phi_x,zero,SIZE(cx_phi_x))

      allocate( tr_phi_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate tr_phi_x istat',istat
      call bcast(tr_phi_x,zero,SIZE(tr_phi_x))

      allocate( cx_cnz_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate cx_cnz_x istat',istat
      call bcast(cx_cnz_x,zero,SIZE(cx_cnz_x))

      allocate( tr_cnz_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  tr_cnz_x istat',istat
      call bcast(tr_cnz_x,zero,SIZE(tr_cnz_x))

      allocate( cx_cnr_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  cx_cnr_x istat',istat
      call bcast(cx_cnr_x,zero,SIZE(cx_cnr_x))

      allocate( tr_cnr_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  tr_cnr_x istat',istat
      call bcast(tr_cnr_x,zero,SIZE(tr_cnr_x))

      allocate( cx_cm_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  cx_cm_x  istat',istat
      call bcast(cx_cm_x,zero,SIZE(cx_cm_x))

      allocate( tr_cm_x(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate tr_cm_x  istat',istat
      call bcast(tr_cm_x,zero,SIZE(tr_cm_x))

c-----for o mode ray part     

      allocate( cx_z_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  cx_z_o  istat',istat
      call bcast(cx_z_o,zero,SIZE(cx_z_o))

      allocate( tr_z_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate   tr_z_o  istat',istat
      call bcast(tr_z_o,zero,SIZE(tr_z_o))

      allocate( cx_r_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  cx_r_o  istat',istat
      call bcast(cx_r_o,zero,SIZE(cx_r_o))

      allocate( tr_r_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  tr_r_o  istat',istat
      call bcast(tr_r_o,zero,SIZE(tr_r_o))

      allocate( cx_phi_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  cx_phi_o istat',istat
      call bcast(cx_phi_o,zero,SIZE(cx_phi_o))

      allocate( tr_phi_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  tr_phi_o  istat',istat
      call bcast(tr_phi_o,zero,SIZE(tr_phi_o))

      allocate( cx_cnz_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  cx_cnz_o  istat',istat
      call bcast(cx_cnz_o,zero,SIZE(cx_cnz_o))

      allocate( tr_cnz_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  tr_cnz_o istat',istat
      call bcast(tr_cnz_o,zero,SIZE(tr_cnz_o))

      allocate( cx_cnr_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  cx_cnr_o istat',istat
      call bcast(cx_cnr_o,zero,SIZE(cx_cnr_o))

      allocate( tr_cnr_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  tr_cnr_o istat',istat
      call bcast(tr_cnr_o,zero,SIZE(tr_cnr_o))

      allocate( cx_cm_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  cx_cm_o istat',istat
      call bcast(cx_cm_o,zero,SIZE(cx_cm_o))

      allocate( tr_cm_o(1:nrelta4),STAT=istat)
cyup      write(*,*)'after allocate  tr_cm_o istat',istat
      call bcast(tr_cm_o,zero,SIZE(tr_cm_o))
 
c----- 
cyup      write(*,*)'jx_kin',jx_kin
      allocate( kinetic_energy_kev(1:jx_kin),STAT=istat)
cyup      write(*,*)'after allocate  kinetic_energy_kev istat',istat
      call bcast( kinetic_energy_kev,zero,SIZE( kinetic_energy_kev))

      return
      end


      subroutine ainalloc_six_no_nml_i
c-----allocate pointers in six_no_nml.i
      implicit none
      include 'param.i'
      include 'one.i'
      include 'six.i'
      
      real*8  zero
      integer istat,ndens4

      ndens4=ndens+4
    
      zero=0.d0

      allocate( cxdens1(1:ndens4,1:nbulk),STAT=istat)
      call bcast(cxdens1,zero,SIZE(cxdens1))

      allocate(  cxtemp1(1:ndens4,1:nbulk),STAT=istat)
      call bcast( cxtemp1,zero,SIZE( cxtemp1))

      allocate(  trdens1(1:ndens4,1:nbulk),STAT=istat)
      call bcast( trdens1,zero,SIZE( trdens1))

      allocate(   trtemp1(1:ndens4,1:nbulk),STAT=istat)
      call bcast(  trtemp1,zero,SIZE(  trtemp1))

      allocate(  trvflow1(1:ndens4,1:nbulk),STAT=istat)
      call bcast( trvflow1,zero,SIZE(trvflow1))

      allocate(  cxvflow1(1:ndens4,1:nbulk),STAT=istat)
      call bcast( cxvflow1,zero,SIZE(cxvflow1))

      allocate(  trvflow(1:ndens4),STAT=istat)
      call bcast( trvflow,zero,SIZE(trvflow))

      allocate(   cxvflow(1:ndens4),STAT=istat)
      call bcast(  cxvflow,zero,SIZE( cxvflow))

      allocate(  trtpop1(1:ndens4,1:nbulk),STAT=istat)
      call bcast( trtpop1,zero,SIZE(trtpop1))

      allocate(  cxtpop1(1:ndens4,1:nbulk),STAT=istat)
      call bcast( cxtpop1,zero,SIZE(cxtpop1))

      allocate(  trtpop(1:ndens4),STAT=istat)
      call bcast( trtpop,zero,SIZE( trtpop))

      allocate(  cxtpop(1:ndens4),STAT=istat)
      call bcast(cxtpop,zero,SIZE( cxtpop))

      allocate(  trdens(1:ndens4),STAT=istat)
      call bcast(trdens,zero,SIZE( trdens))

      allocate(  rhom(1:ndens4),STAT=istat)
      call bcast(rhom,zero,SIZE( rhom))

      allocate( densm(1:ndens4),STAT=istat)
      call bcast(densm,zero,SIZE(densm))

      allocate(  cxdens(1:ndens4),STAT=istat)
      call bcast(cxdens,zero,SIZE(cxdens))

      allocate(  trtempe(1:ndens4),STAT=istat)
      call bcast(trtempe,zero,SIZE(trtempe))

      allocate(  trzeff(1:ndens4),STAT=istat)
      call bcast(trzeff,zero,SIZE(trzeff))

      allocate(  cxtempe(1:ndens4),STAT=istat)
      call bcast(cxtempe,zero,SIZE(cxtempe))

      allocate(  cxzeff(1:ndens4),STAT=istat)
      call bcast(cxzeff,zero,SIZE( cxzeff))

      allocate(  d2_dens_drho1(1:ndens4,1:nbulk),STAT=istat)
      call bcast(d2_dens_drho1,zero,SIZE(d2_dens_drho1))

      allocate(  d2_temp_drho1(1:ndens4,1:nbulk),STAT=istat)
      call bcast(d2_temp_drho1,zero,SIZE(d2_temp_drho1))

      allocate(  d2_tpop_drho1(1:ndens4,1:nbulk),STAT=istat)
      call bcast(d2_tpop_drho1,zero,SIZE(d2_tpop_drho1))

      allocate( d2_vflow_drho1(1:ndens4,1:nbulk),STAT=istat)
      call bcast(d2_vflow_drho1,zero,SIZE(d2_vflow_drho1))

      allocate( d2_zeff_drho1(1:ndens4),STAT=istat)
      call bcast(d2_zeff_drho1,zero,SIZE(d2_zeff_drho1))

      return
      end



      subroutine ainalloc_writencdf_i(nray)
c-----allocate pointers in writencdf.i
      implicit none
      include 'param.i'
      include 'one.i'
      include 'onetwo.i'
      include 'writencdf.i'
      include 'emissa.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c-----input 
      integer nray !number of rays

c-----locals
      real*8  zero
      complex*16 compl_zero
      integer istat

      if(myrank.eq.0)WRITE(*,*)
     + 'in ainalloc_writencdf_i nray,nrelt,nfreq,nbulk',
     &                          nray,nrelt,nfreq,nbulk
     
      zero=0.d0
      compl_zero=dcmplx(zero,zero)
   
      allocate( ws_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate ws_nc istat',istat
      call bcast(ws_nc,zero,SIZE(ws_nc))

      allocate( seikon_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate seikon_nc istat',istat
      call bcast(seikon_nc,zero,SIZE(seikon_nc))

      allocate( spsi_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate spsi_nc istat',istat
      call bcast(spsi_nc,zero,SIZE(spsi_nc))

      allocate( wr_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  wr_nc istat',istat
      call bcast(wr_nc,zero,SIZE(wr_nc))

      allocate( wphi_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate wphi_nc istat',istat
      call bcast(wphi_nc,zero,SIZE(wphi_nc))

      allocate( wz_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate wz_nc istat',istat
      call bcast(wz_nc,zero,SIZE(wz_nc))

      allocate( wnpar_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  wnpar_nc istat',istat
      call bcast(wnpar_nc,zero,SIZE(wnpar_nc))

      allocate( wnper_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  wnper_nc istat',istat
      call bcast(wnper_nc,zero,SIZE(wnper_nc))

      allocate( delpwr_nc(1:nrelta,1:nray*nfreq),STAT=istat)
ccyup      write(*,*)'after allocate  delpwr_nc istat',istat
      call bcast(delpwr_nc,zero,SIZE(delpwr_nc))

      allocate( sdpwr_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  sdpwr_nc istat',istat
      call bcast(sdpwr_nc,zero,SIZE(sdpwr_nc))

      allocate( wdnpar_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  wdnpar_nc istat',istat
      call bcast(wdnpar_nc,zero,SIZE(wdnpar_nc))

      allocate( fluxn_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  fluxn_nc istat',istat
      call bcast(fluxn_nc,zero,SIZE(fluxn_nc))

      allocate( sbtot_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate sbtot_nc istat',istat
      call bcast(sbtot_nc,zero,SIZE(sbtot_nc))

      allocate( sb_z_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  sb_z_nc istat',istat
      call bcast(sb_z_nc,zero,SIZE(sb_z_nc))

      allocate( sb_r_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  sb_r_nc istat',istat
      call bcast(sb_r_nc,zero,SIZE(sb_r_nc))

      allocate( sb_phi_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  sb_phi_nc istat',istat
      call bcast(sb_phi_nc,zero,SIZE(sb_phi_nc))
      
      allocate( sene_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate sene_nc istat',istat
      call bcast(sene_nc,zero,SIZE(sene_nc))

      allocate( ste_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  ste_nc istat',istat
      call bcast(ste_nc,zero,SIZE(ste_nc))
      
      allocate( szeff_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate szeff_nc istat',istat
      call bcast(szeff_nc,zero,SIZE(szeff_nc))

      allocate( salphac_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  salphac_nc istat',istat
      call bcast(salphac_nc,zero,SIZE(salphac_nc))

      allocate( salphal_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  salphal_nc istat',istat
      call bcast(salphal_nc,zero,SIZE(salphal_nc))
      
      allocate( salphas_nc(1:nrelta,1:nray*nfreq,1:nbulk),STAT=istat)
cyup      write(*,*)'after allocate  salphas_nc istat',istat
      call bcast(salphas_nc,zero,SIZE(salphas_nc))
      
      allocate( vgr_z_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate vgr_z_nc istat',istat
      call bcast(vgr_z_nc,zero,SIZE(vgr_z_nc))
                     
      allocate( vgr_r_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate vgr_r_nc istat',istat
      call bcast(vgr_r_nc,zero,SIZE(vgr_r_nc))
            
      allocate( vgr_phi_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate vgr_phi_nc istat',istat
      call bcast(vgr_phi_nc,zero,SIZE(vgr_phi_nc))
            
      allocate( flux_z_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate flux_z_nc istat',istat
      call bcast(flux_z_nc,zero,SIZE(flux_z_nc))
            
      allocate( flux_r_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate flux_r_nc istat',istat
      call bcast(flux_r_nc,zero,SIZE(flux_r_nc))
       
      allocate( flux_phi_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate flux_phi_nc istat',istat
      call bcast(flux_phi_nc,zero,SIZE(flux_phi_nc))
            
      allocate( wn_r_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  wn_r_nc istat',istat
      call bcast(wn_r_nc,zero,SIZE(wn_r_nc))
      
      allocate( wn_z_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  wn_z_nc istat',istat
      call bcast(wn_z_nc,zero,SIZE(wn_z_nc))
            
      allocate( wn_phi_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  wn_phi_nc istat',istat
      call bcast(wn_phi_nc,zero,SIZE(wn_phi_nc))
      
      allocate( transm_ox_nc(1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  transm_ox_nc istat',istat
      call bcast(transm_ox_nc,zero,SIZE(transm_ox_nc))
      
      allocate( cn_par_optimal_nc(1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  cn_par_optimal_nc istat',istat
      call bcast(cn_par_optimal_nc,zero,SIZE(cn_par_optimal_nc))
      
      allocate( cnpar_ox_nc(1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  cnpar_ox_nc istat',istat
      call bcast(cnpar_ox_nc,zero,SIZE(cnpar_ox_nc))
      
      allocate( cn_b_gradpsi_nc(1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  cn_b_gradpsi_nc istat',istat
      call bcast(cn_b_gradpsi_nc,zero,SIZE(cn_b_gradpsi_nc))

      if (i_emission.eq.1) then
c--------emission data--------------------------      
 
         allocate( wsn_nc(1:nrelta,1:nfreq),STAT=istat)   
         write(*,*)'after allocate wsn_nc istat',istat
         call bcast(wsn_nc,zero,SIZE(wsn_nc))
 
         allocate( wcnpar_em_nc(1:nrelta,1:nfreq),STAT=istat) 
         write(*,*)'after allocate  wcnpar_em_n istat',istat
         call bcast(wcnpar_em_nc,zero,SIZE(wcnpar_em_nc))
 
         allocate( wcnper_em_nc(1:nrelta,1:nfreq),STAT=istat) 
         write(*,*)'after allocate  wcnper_em_nc istat',istat
         call bcast(wcnper_em_nc,zero,SIZE(wcnper_em_nc))
 
         allocate( wz_em_nc(1:nrelta,1:nfreq),STAT=istat)     
         write(*,*)'after allocate  wz_em_nc istat',istat
         call bcast(wz_em_nc,zero,SIZE(wz_em_nc))
 
         allocate( wr_em_nc(1:nrelta,1:nfreq),STAT=istat)     
         write(*,*)'after allocate  wr_em_nc istat',istat
         call bcast(wr_em_nc,zero,SIZE(wr_em_nc))
 
         allocate( wphi_em_nc(1:nrelta,1:nfreq),STAT=istat)    
         write(*,*)'after allocate wphi_em_nc  istat',istat
         call bcast(wphi_em_nc,zero,SIZE(wphi_em_nc))
 
         allocate( wal_emis_nc(1:nrelta,1:nfreq),STAT=istat)   
         write(*,*)'after allocate wal_emis_nc  istat',istat
         call bcast(wal_emis_nc,zero,SIZE(wal_emis_nc))
 
         allocate( wj_emis_nc(1:nrelta,1:nfreq),STAT=istat)    
         write(*,*)'after allocate wj_emis_nc  istat',istat
         call bcast(wj_emis_nc,zero,SIZE(wj_emis_nc))
 
         allocate( wnray_nc(1:nrelta,1:nfreq),STAT=istat)       
         write(*,*)'after allocate  wnray_nc istat',istat
         call bcast(wnray_nc,zero,SIZE(wnray_nc))
 
         allocate( win_sn_nc(1:nrelta,1:nfreq),STAT=istat) 
         write(*,*)'after allocate  win_sn_nc istat',istat
         call bcast(win_sn_nc,zero,SIZE(win_sn_nc))
 
         allocate( win_0_nc(1:nrelta,1:nfreq),STAT=istat)
         write(*,*)'after allocate win_0_nc istat',istat
         call bcast(win_0_nc,zero,SIZE(win_0_nc))
 
         allocate( w_specific_intensity_nc(1:nrelta,1:nfreq),STAT=istat)
         write(*,*)'after allocate w_specific_intensity_nc istat',istat
         call bcast(w_specific_intensity_nc,zero,
     &              SIZE(w_specific_intensity_nc))
 
         allocate( wi_0_nc(1:nrelta,1:nfreq),STAT=istat) 
         write(*,*)'after allocate wi_0_nc istat',istat
         call bcast(wi_0_nc,zero,SIZE(wi_0_nc))
 
         allocate( wi_0sn_nc(1:nrelta,1:nfreq),STAT=istat) 
         write(*,*)'after allocate wi_0sn_nc istat',istat
         call bcast(wi_0sn_nc,zero,SIZE(wi_0sn_nc)) 

         allocate( wtemp_em_nc(1:nrelta,1:nfreq),STAT=istat)
         write(*,*)'after allocate wtemp_em_nc istat',istat
         call bcast(wtemp_em_nc,zero,SIZE(wtemp_em_nc)) 

         allocate( wtemp_rad_em_nc(1:nrelta,1:nfreq),STAT=istat)
         write(*,*)'after allocate  wtemp_rad_em_nc istat',istat
         call bcast(wtemp_rad_em_nc,zero,SIZE(wtemp_rad_em_nc)) 

         allocate( wtemp_rad_fr_wall_nc(1:nray,1:nfreq),STAT=istat)
         write(*,*)'after allocate wtemp_rad_fr_wall_nc istat',istat
         call bcast(wtemp_rad_fr_wall_nc,zero,
     &              SIZE(wtemp_rad_fr_wall_nc))

         allocate( waveraged_temp_rad_fr_nc(1:nfreq),STAT=istat)
         write(*,*)'after allocate waveraged_temp_rad_fr_nc istat',istat
         call bcast(waveraged_temp_rad_fr_nc,zero,
     &              SIZE(waveraged_temp_rad_fr_nc)) 

         allocate( waveraged_temp_rad_fr_wall_nc(1:nfreq),
     &             STAT=istat)
         write(*,*)'after allocate  waveraged_temp_rad_fr_wall_nc
     &              istat',istat
         call bcast(waveraged_temp_rad_fr_wall_nc,zero,
     &              SIZE(waveraged_temp_rad_fr_wall_nc)) 

         allocate( wtemp_rad_fr_nc(1:nray,1:nfreq),STAT=istat)
         write(*,*)'after allocate  wtemp_rad_fr_nc istat',istat
         call bcast(wtemp_rad_fr_nc,zero,SIZE(wtemp_rad_fr_nc))

         allocate( wtemp_pl_fr_nc(1:nray,1:nfreq),STAT=istat)
         write(*,*)'after allocate   wtemp_pl_fr_nc istat',istat
         call bcast(wtemp_pl_fr_nc,zero,SIZE(wtemp_pl_fr_nc))

         allocate( wr0_em_nc(1:nray,1:nfreq),STAT=istat)
         write(*,*)'after allocate  wr0_em_nc istat',istat
         call bcast(wr0_em_nc,zero,SIZE(wr0_em_nc))

         allocate( wz0_em_nc(1:nfreq),STAT=istat)
         write(*,*)'after allocate  wz0_em_nc istat',istat
         call bcast(wz0_em_nc,zero,SIZE(wz0_em_nc))

         allocate( wrho0_em_nc(1:nray,1:nfreq),STAT=istat)
         write(*,*)'after allocate  wrho0_em_nc istat',istat
         call bcast(wrho0_em_nc,zero,SIZE(wrho0_em_nc))

         write(*,*)'ainalloc_writencdf befor allocate(wtaun_em_nc'
     &   ,'nrelta,nfreq',nrelta,nfreq

         allocate( wtaun_em_nc(1:nrelta,1:nfreq),STAT=istat)
         write(*,*)'after allocate  wtaun_em_nc istat',istat
         call bcast(wtaun_em_nc,zero,SIZE(wtaun_em_nc))

         allocate( wfreq_nc(1:nfreq),STAT=istat)
         write(*,*)'after allocate  wfreq_nc istat',istat
         call bcast(wfreq_nc,zero,SIZE(wfreq_nc))

         allocate( wtau_em_nc(1:nray,1:nfreq),STAT=istat)
         write(*,*)'after allocate  wtau_em_nc istat',istat
         call bcast(wtau_em_nc,zero,SIZE(wtau_em_nc))

         allocate( wi_0t_nc(1:nray,1:nfreq),STAT=istat)
         write(*,*)'after allocate wi_0t_nc istat',istat
         call bcast(wi_0t_nc,zero,SIZE(wi_0t_nc))

         allocate( wr_2nd_harm_nc(1:nray,1:nfreq),STAT=istat)
         write(*,*)'after allocate  wr_2nd_harm_nc istat',istat
         call bcast(wr_2nd_harm_nc,zero,SIZE(wr_2nd_harm_nc))

         allocate( wtemp_2nd_harm_nc(1:nray,1:nfreq),STAT=istat)
         write(*,*)'after allocate  wtemp_2nd_harm_nc istat',istat
         call bcast(wtemp_2nd_harm_nc,zero,SIZE(wtemp_2nd_harm_nc))

         allocate(wr_emis_initial_mesh_nc(1:nrelta,1:nfreq),STAT=istat)
         write(*,*)'after allocate wr_emis_initial_mesh_nc istat',istat
         call bcast(wr_emis_initial_mesh_nc,zero,
     &              SIZE(wr_emis_initial_mesh_nc))

         allocate(wp_perpmax_dmc_nc(1:nrelta,n_relt_harm1:n_relt_harm2,
     &                              1:nfreq),STAT=istat)
         write(*,*)'after allocate(wp_perpmax_dmc_nc istat',istat
         call bcast(wp_perpmax_dmc_nc,zero,SIZE(wp_perpmax_dmc_nc))

         allocate(wp_parmin_dmc_nc(1:nrelta,n_relt_harm1:n_relt_harm2,
     &                              1:nfreq),STAT=istat)
         write(*,*)'after allocate wp_parmin_dmc_nc istat',istat
         call bcast(wp_parmin_dmc_nc,zero,SIZE(wp_parmin_dmc_nc))

         allocate(wp_parmax_dmc_nc(1:nrelta,n_relt_harm1:n_relt_harm2,
     &                              1:nfreq),STAT=istat)
         write(*,*)'after allocate wp_parmax_dmc_nc istat',istat
         call bcast(wp_parmax_dmc_nc,zero,SIZE(wp_parmax_dmc_nc))

         allocate( wj_emis_x_nc(1:nrelta,1:jx_kin,1:nfreq),STAT=istat)
         write(*,*)'after allocate wj_emis_x_nc istat',istat
         call bcast(wj_emis_x_nc,zero,SIZE(wj_emis_x_nc))

         allocate( win_sn_x_nc(1:nrelta,1:jx_kin,1:nfreq),STAT=istat)
         write(*,*)'after allocate  win_sn_x_nc istat',istat
         call bcast(win_sn_x_nc,zero,SIZE(win_sn_x_nc))

         allocate( win_0_x_nc(1:nrelta,1:jx_kin,1:nfreq),STAT=istat)
         write(*,*)'after allocate   win_0_x_nc istat',istat
         call bcast(win_0_x_nc,zero,SIZE(win_0_x_nc))

         allocate( wi_0sn_x_nc(1:nrelta,1:jx_kin,1:nfreq),STAT=istat)
         write(*,*)'after allocate  wi_0sn_x_nc istat',istat
         call bcast(wi_0sn_x_nc,zero,SIZE(wi_0sn_x_nc))

         allocate( wi_0_x_nc(1:jx_kin,1:nfreq),STAT=istat)
         write(*,*)'after allocate   wi_0_x_nc istat',istat
         call bcast(wi_0_x_nc,zero,SIZE(wi_0_x_nc))

         allocate( w_specific_intensity_x_nc(1:nrelta,1:jx_kin,1:nfreq),
     &             STAT=istat)
         write(*,*)'aft allocate w_specific_intensity_x_nc istat',istat
         call bcast(w_specific_intensity_x_nc,zero,
     &              SIZE(w_specific_intensity_x_nc))

         allocate( wdye_nc(1:nrelta,1:nfreq),STAT=istat)
         write(*,*)'after allocate wdye_nc istat',istat
         call bcast(wdye_nc,zero,SIZE(wdye_nc))

         allocate( nrayelt_emis_nc(1:nfreq),STAT=istat)
         write(*,*)'after allocate  nrayelt_emis_nc istat',istat
         call ibcast(nrayelt_emis_nc,0,SIZE(nrayelt_emis_nc))

         allocate( nrayelt_emis_initial_mesh_nc(1:nfreq),STAT=istat)
         write(*,*)'after allocate  nrayelt_emis_initial_mesh_nc istat',
     &                                                           istat
         call ibcast(nrayelt_emis_initial_mesh_nc,0,
     &              SIZE(nrayelt_emis_initial_mesh_nc))

      endif
c------------------plasma profiles VS r at z=0  ------------------    

      allocate( w_dens_vs_r_nc(1:NR,1:nbulk),STAT=istat)
cyup      write(*,*)'after allocate w_dens_vs_r_nc istat',istat
      call bcast( w_dens_vs_r_nc,zero,SIZE( w_dens_vs_r_nc))

      allocate( w_temp_vs_r_nc(1:NR,1:nbulk),STAT=istat)
cyup      write(*,*)'after allocate  w_temp_vs_r_nc istat',istat
      call bcast( w_temp_vs_r_nc,zero,SIZE( w_temp_vs_r_nc))

      allocate( w_zeff_vs_r_nc(1:NR),STAT=istat)
cyup      write(*,*)'after allocate   w_zeff_vs_r_nc istat',istat
      call bcast( w_zeff_vs_r_nc,zero,SIZE( w_zeff_vs_r_nc))

      allocate( w_r_densprof_nc(1:NR),STAT=istat)
cyup      write(*,*)'after allocate   w_r_densprof_nc istat',istat
      call bcast( w_r_densprof_nc,zero,SIZE( w_r_densprof_nc))

c--------------------------------------------------------------
      allocate( w_eff_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  w_eff_nc istat',istat
      call bcast( w_eff_nc,zero,SIZE( w_eff_nc))

      allocate( w_theta_pol_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate   w_theta_pol_nc istat',istat
      call bcast( w_theta_pol_nc,zero,SIZE( w_theta_pol_nc))

c-----electric field complex polarization

      allocate( cwexde_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  cwexde_nc istat',istat
      call ccast(cwexde_nc,zero,SIZE( cwexde_nc))

      allocate( cweyde_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  cweyde_nc istat',istat
      call ccast(cweyde_nc,zero,SIZE( cweyde_nc))

      allocate( cwezde_nc(1:nrelta,1:nray*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  cwezde_nc istat',istat
      call ccast(cwezde_nc,zero,SIZE( cwezde_nc))

      if (dielectric_op.eq.'enabled') then
c--------write dielectric tensor elements

         allocate( cweps11_nc(1:nrelta,1:nray*nfreq),STAT=istat)
         write(*,*)'after allocate  cweps11_nc istat',istat
         call ccast( cweps11_nc,zero,SIZE( cweps11_nc))

         allocate( cweps12_nc(1:nrelta,1:nray*nfreq),STAT=istat)
         write(*,*)'after allocate  cweps12_nc istat',istat
         call ccast( cweps12_nc,zero,SIZE( cweps12_nc))

         allocate( cweps13_nc(1:nrelta,1:nray*nfreq),STAT=istat)
         write(*,*)'after allocate  cweps13_nc istat',istat
         call ccast( cweps13_nc,zero,SIZE( cweps13_nc))

         allocate( cweps21_nc(1:nrelta,1:nray*nfreq),STAT=istat)
         write(*,*)'after allocate  cweps21_nc istat',istat
         call ccast( cweps21_nc,zero,SIZE( cweps21_nc))

         allocate( cweps22_nc(1:nrelta,1:nray*nfreq),STAT=istat)
         write(*,*)'after allocate  cweps22_nc istat',istat
         call ccast( cweps22_nc,zero,SIZE( cweps22_nc))

         allocate( cweps23_nc(1:nrelta,1:nray*nfreq),STAT=istat)
         write(*,*)'after allocate  cweps23_nc istat',istat
         call ccast( cweps23_nc,zero,SIZE( cweps23_nc))

         allocate( cweps31_nc(1:nrelta,1:nray*nfreq),STAT=istat)
         write(*,*)'after allocate  cweps31_nc istat',istat
         call ccast( cweps31_nc,zero,SIZE( cweps31_nc))

         allocate( cweps32_nc(1:nrelta,1:nray*nfreq),STAT=istat)
         write(*,*)'after allocate  cweps32_nc istat',istat
         call ccast( cweps32_nc,zero,SIZE( cweps32_nc))

         allocate( cweps33_nc(1:nrelta,1:nray*nfreq),STAT=istat)
         write(*,*)'after allocate  cweps33_nc istat',istat
         call ccast( cweps33_nc,zero,SIZE( cweps33_nc))        

      endif
c----------------------------      
      allocate( i_ox_conversion_nc(1:nray*nfreq),STAT=istat)
cyup      write(*,*) 'after allocate  i_ox_conversion_nc istat=',istat
      call ibcast(i_ox_conversion_nc,0,SIZE( i_ox_conversion_nc))

cyup      write(*,*) 'ainalloc_writencdf nray=',nray,'nfreq',nfreq
      allocate( nrayelt_nc(1:nray*nfreq),STAT=istat)
cyup      write(*,*) 'after allocate nrayelt_nc istat=',istat
      call ibcast(nrayelt_nc,0,SIZE( nrayelt_nc))
c-------------------------------------------------------------
      allocate(iray_status_nc(1:nray,1:nfreq),STAT=istat)
cyup      write(*,*) 'after allocate iray_status istat=',istat
      call ibcast(iray_status_nc,0,SIZE(iray_status_nc))

      return
      end


      subroutine ainalloc_dskin_i
c-----allocate pointers in dskin.i and in spline_distrib.i
      implicit none
      include 'param.i'
      include 'one.i'
      include 'dskin.i'
      
      real*8  zero
      integer istat

      zero=0.d0

      if (i_diskf.eq.1) then
c--------i_diskf=1 read in the file diskf
        
         if(lrzmax.gt.lrz) then
           write(*,*)'in subroutine ainalloc_dskin_i'
           write(*,*)'i_diskf=1, read in the file diskf'
           write(*,*)'lrzmax > lrz'
           write(*,*)'lrzmax=',lrzmax,'lrz=',lrz
           write(*,*)'the code sets small radial dimension in pointers'
           write(*,*)'as lrz'
           write(*,*)'it can create problem with pointers' 
         endif
      endif

c-----real*8
      allocate( x(1:jx),STAT=istat)
      call bcast(x,zero,SIZE(x))
     
      allocate( y(1:iya,1:lrz),STAT=istat)
      call bcast(y,zero,SIZE(y))

      allocate( rovera(1:lrz),STAT=istat)
      call bcast(rovera,zero,SIZE(rovera))

      allocate( rya(1:lrz),STAT=istat)
      call bcast(rya,zero,SIZE(rya))

      allocate( elecfld(1:lrz),STAT=istat)
      call bcast(elecfld,zero,SIZE(elecfld))

      allocate( bthr(1:lrz),STAT=istat)
      call bcast(bthr,zero,SIZE(bthr))

      allocate( btoru(1:lrz),STAT=istat)
      call bcast(btoru,zero,SIZE(btoru))

      allocate( bthr0(1:lrz),STAT=istat)
      call bcast(bthr0,zero,SIZE(bthr0))

      allocate( btor0(1:lrz),STAT=istat)
      call bcast(btor0,zero,SIZE(btor0))

      allocate( reden(1:lrz,1:ngen),STAT=istat)
      call bcast(reden,zero,SIZE(reden))

      allocate( temp(1:lrz,1:ngen),STAT=istat)
      call bcast(temp,zero,SIZE(temp))

      allocate( bnumb(1:ngen),STAT=istat)
      call bcast(bnumb,zero,SIZE(bnumb))

      allocate( fmass(1:ngen),STAT=istat)
      call bcast(fmass,zero,SIZE(fmass))

cyup      write(*,*)'before  allocate f iya,jx,lrz,ngen',iya,jx,lrz,ngen

      allocate( f(1:iya,1:jx,1:lrz,1:ngen),STAT=istat)
      call bcast(f,zero,SIZE(f))
      write(*,*)'allocate f istat=',istat
      
      allocate( df_dpitch(1:iya+1,1:jx,1:lrz,1:ngen),STAT=istat)
      call bcast(df_dpitch,zero,SIZE(df_dpitch))
  
      allocate( df_dx(1:iya,1:jx+1,1:lrz,1:ngen),STAT=istat)
      call bcast(df_dx,zero,SIZE(df_dx))
    
c-----integer*4
      allocate( iy_(1:lrz),STAT=istat)
      call ibcast(iy_,0,SIZE(iy_))

      allocate( itl(1:lrz),STAT=istat)
      call ibcast(itl,0,SIZE(itl))

      allocate( itu(1:lrz),STAT=istat)
      call ibcast(itu,0,SIZE(itu))

c-----allocate pointers in spline_distrib.i
      call ainalloc_spline_distrib_i

cyup      write(*,*)'dskin.i and spline_distrib.i was allocated'

      return
      end

      subroutine ainalloc_spline_distrib_i
c-----allocate pointers spline_distrib.i
      implicit none
      include 'param.i'
      include 'one.i'
c      include 'dskin.i'
      include 'spline_distrib.i'

      real*8  zero

      integer istat

      zero=0.d0

      nwk=1+3*max(jx,iya)

c-----real*8 
  
      allocate(fxx(1:iya,1:jx,1:lrz,1:ngen),STAT=istat)
      call bcast(fxx,zero,SIZE(fxx))

      allocate(fyy(1:iya,1:jx,1:lrz,1:ngen),STAT=istat)
      call bcast(fyy,zero,SIZE(fyy))

      allocate(fxxyy(1:iya,1:jx,1:lrz,1:ngen),STAT=istat)
      call bcast(fxxyy,zero,SIZE(fxxyy))

      allocate(wk(1:nwk),STAT=istat)
      call bcast(wk,zero,SIZE(wk))

      return
      end

      subroutine ainalloc
c-----allocate pointers at
c     onetwo_no_nml.i 
c     emissa_no_nml.i
c     onetwo_no_nml.i 


      implicit none
      include 'param.i'
      include 'one.i'
      include 'adj_nml.i'
c-----input
      integer nray !!total number of rays 

      
      call ainalloc_six_no_nml_i

      if (i_emission.eq.1) call ainalloc_emissa_no_nml_i

      if (ionetwo.eq.1) call ainalloc_onetwo_no_nml_i

      write(*,*)'in ainalloc before ainalloc_adj_no_nml_i'

      if (i_adj.eq.1) call ainalloc_adj_no_nml_i


      return
      end


      subroutine ainalloc_write_i(nray)
c-----allocate pointers in write.i
      implicit none
      include 'param.i'
      include 'one.i'
      include 'onetwo.i'
      include 'write.i'
      include 'emissa.i'

c-----input 
      integer nray !number of rays

c-----locals
      real*8  zero
      
      integer istat
    
      zero=0.d0

cBH130508   If i_emission.ne.1, the following arrays will
cBH130508   still be allocated, for simplicity of coding related
cBH130508   to their use in subroutines sendrecv_alloc, senddata,
cBH130508   and recvdata, which are inserted in the code during make
cBH130508   for MPI purposes (using file ./mpi/mpi.ins).  These MPI
cBH130508   related subroutines use common/write/... as in write.i.
cBH130508   The arrays are declared as pointer variables in write.i.
cBH130508   [Help in this MPI bug fix from YuP.]
cBH130508      if (i_emission.eq.1) then
c--------for emission calculations
         write(*,*)' ainalloc_write_i nray,nfreq',nray,nfreq


cBH130508         allocate(  wi_0(1:nray,1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate(wi_0 istat',istat
cBH130508         call bcast( wi_0,zero,SIZE(wi_0))

cBH130508         allocate(  wtemp_rad_fr_wall(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate(wtemp_rad_fr_wall istat',istat
cBH130508         call bcast( wtemp_rad_fr_wall,zero,SIZE(wtemp_rad_fr_wall))

cBH130508         allocate(  wtemp_rad_fr(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate(wtemp_rad_fr istat',istat
cBH130508         call bcast( wtemp_rad_fr,zero,SIZE(wtemp_rad_fr))

cBH130508         allocate(  wtemp_pl_fr(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate(wtemp_pl_fr istat',istat
cBH130508         call bcast( wtemp_pl_fr,zero,SIZE(wtemp_pl_fr))

cBH130508         allocate(  wr0_em(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate(wr0_em istat',istat
cBH130508         call bcast( wr0_em,zero,SIZE(wr0_em))

cBH130508         allocate(  wz0_em(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate(wz0_em istat',istat
cBH130508         call bcast( wz0_em,zero,SIZE(wz0_em))

cBH130508         allocate(  wrho0_em(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate(wrho0_em istat',istat
cBH130508         call bcast( wrho0_em,zero,SIZE(wrho0_em))

cBH130508         allocate(  wfreq(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate(wfreq istat',istat
cBH130508        call bcast( wfreq,zero,SIZE(wfreq))

cBH130508         allocate(  wtau_em(1:nray,1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wtau_em istat',istat
cBH130508         call bcast( wtau_em,zero,SIZE(wtau_em))

cBH130508         allocate(  wi_0t(1:nray,1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wi_0t istat',istat
cBH130508         call bcast( wi_0t,zero,SIZE(wi_0t))

cBH130508         allocate(  wr_2nd_harm(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wr_2nd_harm istat',istat
cBH130508         call bcast( wr_2nd_harm,zero,SIZE(wr_2nd_harm))

cBH130508         allocate(  wtemp_2nd_harm(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wtemp_2nd_harm istat',istat
cBH130508         call bcast( wtemp_2nd_harm,zero,SIZE(wtemp_2nd_harm))

c-------------------to plot cold N_perp in initial ray point for test

cBH130508         allocate(  wn_perp_ioxm_p(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wn_perp_ioxm_p istat',istat
cBH130508         call bcast( wn_perp_ioxm_p,zero,SIZE(wn_perp_ioxm_p))

cBH130508         allocate(  wn_perp_ioxm_m(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wn_perp_ioxm_m istat',istat
cBH130508         call bcast( wn_perp_ioxm_m,zero,SIZE(wn_perp_ioxm_m))

cBH130508         allocate(  wn_perp_ioxm_n_npar_p(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wn_perp_ioxm_n_npar_p istat',istat
cBH130508         call bcast(wn_perp_ioxm_n_npar_p,zero,
cBH130508     &              SIZE(wn_perp_ioxm_n_npar_p))

cBH130508         allocate(  wn_perp_ioxm_n_npar_m(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wn_perp_ioxm_n_npar_m istat',istat
cBH130508         call bcast(wn_perp_ioxm_n_npar_m,zero,
cBH130508     &              SIZE(wn_perp_ioxm_n_npar_m))

cBH130508         allocate(  wye_0(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wye_0 istat',istat
cBH130508         call bcast(wye_0,zero,SIZE(wye_0))

cBH130508         allocate(  wxe_0(1:nfreq),STAT=istat)
cBH130508         write(*,*)'after allocate( wxe_0 istat',istat
cBH130508         call bcast(wxe_0,zero,SIZE(wxe_0))

cBH130508      else
c---------these array are used for some plotting at i_emission=0 
cBH130508         allocate(  wn_perp_ioxm_p(1:1),STAT=istat)
cBH130508         write(*,*)'after allocate( wn_perp_ioxm_p',istat
cBH130508         call bcast( wn_perp_ioxm_p,zero,SIZE(wn_perp_ioxm_p))

cBH130508         allocate(  wn_perp_ioxm_m(1:1),STAT=istat)
cBH130508         write(*,*)'after allocate( wn_perp_ioxm_m',istat
cBH130508         call bcast( wn_perp_ioxm_m,zero,SIZE(wn_perp_ioxm_m))

cBH130508         allocate(  wye_0(1:1),STAT=istat)
cBH130508         write(*,*)'after allocate( wye_0',istat
cBH130508         call bcast(wye_0,zero,SIZE(wye_0))

cBH130508         allocate(  wxe_0(1:1),STAT=istat)
cBH130508         write(*,*)'after allocate( wxe_0',istat
cBH130508         call bcast(wxe_0,zero,SIZE(wxe_0))

cBH130508      endif !(i_emission.eq.1)

      return
      end

      subroutine ainalloc_adj_no_nml_i
c-----allocate pointers in adj_no_nml.i
      implicit none
      include 'adj.i'
    
c-----locals
      real*8  zero
      
      integer istat
 
      zero=0.d0

cyup      write(*,*)'in ainalloc_adj_no_nml_i npsi0',npsi0
     
      allocate(psis(1:npsi0),STAT=istat)
      call bcast(psis,zero,SIZE(psis))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after psis istat',istat

cyup      write(*,*)'in ainalloc_adj_no_nml_i npsi0,nthp0',npsi0,nthp0
      allocate(cxc(1:nthp0+1,1:npsi0),STAT=istat)
      call bcast(cxc,zero,SIZE(cxc))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after cxc istat',istat

      allocate(czc(1:nthp0+1,1:npsi0),STAT=istat)
      call bcast(czc,zero,SIZE(czc))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after czc istat',istat

      allocate(bbar(1:npsi0),STAT=istat)
      call bcast(bbar,zero,SIZE(bbar))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after bbar istat',istat

c      write(*,*)'in ainalloc_adj_no_nml_i imax_chi',imax_chi
c      allocate(th0a(1:imax_chi),STAT=istat)
c      call bcast(th0a,zero,SIZE(th0a))
c      write(*,*)'in ainalloc_adj_no_nml_i after th0a istat',istat

cyup      write(*,*)'in ainalloc_adj_no_nml_i nmax_chi',nmax_chi
      allocate(ua(1:nmax_chi),STAT=istat)
      call bcast(ua,zero,SIZE(ua))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after ua istat',istat

      allocate(dene_chi(1:npsi0),STAT=istat)
      call bcast(dene_chi,zero,SIZE(dene_chi))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after dene_ch istat',istat

      allocate(teme_chi(1:npsi0),STAT=istat)
      call bcast(teme_chi,zero,SIZE(teme_chi))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after teme_ch istat',istat

      allocate(zi_chi(1:npsi0),STAT=istat)
      call bcast(zi_chi,zero,SIZE(zi_chi))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after zi_ch istat',istat

      allocate(chi_3d(1:imax_chi,1:nmax_chi,1:npsi0),STAT=istat)
      call bcast(chi_3d,zero,SIZE(chi_3d))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after ch_3d istat',istat

      allocate(bmin_chi(1:npsi0),STAT=istat)
      call bcast(bmin_chi,zero,SIZE(bmin_chi))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after bmim_chi istat',istat

      allocate(bmax_chi(1:npsi0),STAT=istat)
      call bcast(bmax_chi,zero,SIZE(bmax_chi))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after bmax_chi istat',istat

c      allocate(umax_chi(1:npsi0),STAT=istat)
c      call bcast(umax_chi,zero,SIZE(umax_chi))
c      write(*,*)'in ainalloc_adj_no_nml_i after umax_chi istat',istat

      allocate(chi_tt(1:imax_chi,1:nmax_chi,1:npsi0),STAT=istat)
      call bcast(chi_tt,zero,SIZE(chi_tt))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after chi_tt istat',istat

      allocate(chi_uu(1:imax_chi,1:nmax_chi,1:npsi0),STAT=istat)
      call bcast(chi_uu,zero,SIZE(chi_uu))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after chi_uu istat',istat

      allocate(chi_ttuu(1:imax_chi,1:nmax_chi,1:npsi0),STAT=istat)
      call bcast(chi_ttuu,zero,SIZE(chi_ttuu))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after chi_ttuu istat',istat

c      allocate(th0max_chi(1:npsi0),STAT=istat)
c      call bcast(th0max_chi,zero,SIZE(th0max_chi))
ccyup      write(*,*)'in ainalloc_adj_no_nml_i after th0max_chi istat',istat

      allocate(th0a_chi(1:imax_chi,1:npsi0),STAT=istat)
      call bcast(th0a_chi,zero,SIZE(th0a_chi))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after th0a_chi istat',istat

      allocate(thetae_1(0:nthp0-1),STAT=istat)
      call bcast(thetae_1,zero,SIZE(thetae_1))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after thetae_1 istat',istat

      allocate(ba_2d(0:nthp0-1,1:npsi0),STAT=istat)
      call bcast(ba_2d,zero,SIZE(ba_2d))
cyup      write(*,*)'in ainalloc_adj_no_nml_i after ba_2 istat',istat

      allocate(d_chi_d_theta(1:imax_chi,1:nmax_chi,1:npsi0),STAT=istat)
      call bcast(d_chi_d_theta,zero,SIZE(d_chi_d_theta))
cyup      write(*,*)'ainalloc_adj_no_nml_i after d_chi_d_theta istat',istat

      allocate(d_chi_d_u(1:imax_chi,1:nmax_chi,1:npsi0),STAT=istat)
      call bcast(d_chi_d_u,zero,SIZE(d_chi_d_u))
cyup      write(*,*)'ainalloc_adj_no_nml_i after d_chi_d_u istat',istat

      allocate(Q_safety_adj(1:npsi0),STAT=istat)
      call bcast(Q_safety_adj,zero,SIZE(Q_safety_adj))
cyup      write(*,*)'ainalloc_adj_no_nml_i after Q_safety_adj istat',istat

      allocate(sigma_E_adj(1:npsi0),STAT=istat)
      call bcast(sigma_E_adj,zero,SIZE(sigma_E_adj))
cyup      write(*,*)'ainalloc_adj_no_nml_i after sigma_E_adj istat',istat

      allocate(dens_averaged_ar(1:npsi0),STAT=istat)
      call bcast(dens_averaged_ar,zero,SIZE(dens_averaged_ar))
cyup      write(*,*)'ainalloc_adj_no_nml_i aft dens_averaged_ar istat',istat

cyup      write(*,*)'npsi0',npsi0
      allocate(rho_adj(1:npsi0),STAT=istat)
      call bcast(rho_adj,zero,SIZE(rho_adj))
cyup      write(*,*)'ainalloc_adj_no_nml_i after rho_adj istat',istat

cyup      write(*,*)'in ainalloc_adj_no_nml_i before end'
      return
      end

      subroutine ainalloc_lsc_approach_no_nml_i
c-----allocate pointers in lsc_approach_no_nml
      implicit none
      include 'param.i'
      include 'one.i'
      include 'lsc_approach.i' 
      include 'onetwo.i'
      include 'writencdf.i'
c-----locals
      real*8  zero
     
      integer istat
c-----input 
      integer nray !number of rays

cyup      write(*,*)'in lsc_approach_no_nml nv_lsc_ql_lh', nv_lsc_ql_lh

      allocate(rho_bin_lsc(1:n_psi_TSC+1),STAT=istat)
      call bcast(rho_bin_lsc,zero,SIZE(rho_bin_lsc))
cyup      write(*,*)'in lsc_approach_no_nml aft (rho_bin_lsc istat',
cyup     &istat 
    
      allocate(rho_bin_center_lsc(1:n_psi_TSC),STAT=istat)
      call bcast(rho_bin_center_lsc,zero,SIZE(rho_bin_center_lsc))
cyup      write(*,*)'in lsc_approach_no_nml aft (rho_bin_center_lsc istat',
cyup     &istat 

      allocate(binvol_lsc(1:n_psi_TSC),STAT=istat)
      call bcast(binvol_lsc,zero,SIZE(binvol_lsc))
cyup      write(*,*)'in lsc_approach_no_nml aft binvol_lsc istat',
cyup     &istat 


      allocate(binarea_lsc(1:n_psi_TSC),STAT=istat)
      call bcast(binarea_lsc,zero,SIZE(binarea_lsc))
cyup      write(*,*)'in lsc_approach_no_nml aft binarea_lsc istat',
cyup     &istat 


cSAP100122 NR-1 ->  n_psi_TSC
c      allocate(v_te_bin_center_lsc(1:NR-1),STAT=istat)
      allocate(v_te_bin_center_lsc(1:n_psi_TSC),STAT=istat)
      call bcast(v_te_bin_center_lsc,zero,SIZE(v_te_bin_center_lsc))
cyup      write(*,*)'in lsc_approach_no_nml aft v_te_bin_center_lsc istat',
cyup     &istat 

c      allocate(tau_n_bin_center_lsc(1:NR-1),STAT=istat)
      allocate(tau_n_bin_center_lsc(1:n_psi_TSC),STAT=istat)
      call bcast(tau_n_bin_center_lsc,zero,SIZE(tau_n_bin_center_lsc))
cyup      write(*,*)'in lsc_approach_no_nml aft tau_n_bin_center_lsc istat',
cyup     &istat 

c      allocate(del_s_pol_bin_lsc(1:NR-1),STAT=istat
      allocate(del_s_pol_bin_lsc(1:n_psi_TSC),STAT=istat)
      call bcast(del_s_pol_bin_lsc,zero,SIZE(del_s_pol_bin_lsc))
cyup      write(*,*)'in lsc_approach_no_nml after del_s_pol_bin_lsc',
cyup     &istat 

cyup      write(*,*)'in lsc_approach_no_nml bef allcate v_par_mesh_lsc
cyup     &nv_lsc_ql_lh',nv_lsc_ql_lh
      allocate(v_par_mesh_lsc(-nv_lsc_ql_lh : nv_lsc_ql_lh),STAT=istat)
      call bcast(v_par_mesh_lsc,zero,SIZE(v_par_mesh_lsc))
cyup      write(*,*)'in lsc_approach_no_nml after v_par_mesh_lsc istat',
cyup     &istat 

c      allocate(d_ql_lsc_lh_ar(1:NR-1,-nv_lsc_ql_lh:nv_lsc_ql_lh),
      allocate(d_ql_lsc_lh_ar(1:n_psi_TSC,-nv_lsc_ql_lh:nv_lsc_ql_lh),
     &STAT=istat)
      call bcast(d_ql_lsc_lh_ar,zero,SIZE(d_ql_lsc_lh_ar))
cyup      write(*,*)'in lsc_approach_no_nml after d_ql_lsc_lh_ar istat',
cyup     &istat

      allocate(d_ql_lsc_lh_ar_loc(-nv_lsc_ql_lh:nv_lsc_ql_lh),
     &STAT=istat)
      call bcast(d_ql_lsc_lh_ar_loc,zero,SIZE(d_ql_lsc_lh_ar_loc))
cyup      write(*,*)'in lsc_approach_no_nml after d_ql_lsc_lh_ar_loc istat',
cyup     &istat

c      allocate(integral_fe_lsc(1:NR-1,-nv_lsc_ql_lh:nv_lsc_ql_lh),
      allocate(integral_fe_lsc(1:n_psi_TSC,-nv_lsc_ql_lh:nv_lsc_ql_lh)
     &,STAT=istat)
      call bcast(integral_fe_lsc,zero,SIZE(integral_fe_lsc))
cyup      write(*,*)'in lsc_approach_no_nml after integral_fe_lsc istat',
cyup     &istat

c      allocate(fe_lsc(1:NR-1,-nv_lsc_ql_lh:nv_lsc_ql_lh),
      allocate(fe_lsc(1:n_psi_TSC,-nv_lsc_ql_lh:nv_lsc_ql_lh),
     &STAT=istat)
      call bcast(fe_lsc,zero,SIZE(fe_lsc))
cyup      write(*,*)'in lsc_approach_no_nml after fe_lsc istat',
cyup     &istat

c       allocate(d_fe_dv_lsc(1:NR-1,-nv_lsc_ql_lh:nv_lsc_ql_lh),
      allocate(d_fe_dv_lsc(1:n_psi_TSC,-nv_lsc_ql_lh:nv_lsc_ql_lh),
     &STAT=istat)
      call bcast(d_fe_dv_lsc,zero,SIZE(d_fe_dv_lsc))
cyup      write(*,*)'in lsc_approach_no_nml after d_fe_dv_lsc istat',
cyup     &istat

c      allocate(d_ql_lsc_lh_ar_old(1:NR-1,-nv_lsc_ql_lh:nv_lsc_ql_lh),
      allocate(d_ql_lsc_lh_ar_old
     &(1:n_psi_TSC,-nv_lsc_ql_lh:nv_lsc_ql_lh),
     &STAT=istat)
      call bcast(d_ql_lsc_lh_ar_old,zero,SIZE(d_ql_lsc_lh_ar_old))
cyup      write(*,*)'in lsc_approach_no_nml after d_ql_lsc_lh_ar_old istat',
cyup     &istat

c      allocate(fe_lsc_old(1:NR-1,-nv_lsc_ql_lh:nv_lsc_ql_lh),
      allocate(fe_lsc_old(1:n_psi_TSC,-nv_lsc_ql_lh:nv_lsc_ql_lh),
     &STAT=istat)
      call bcast(fe_lsc_old,zero,SIZE(fe_lsc_old))
cyup      write(*,*)'in lsc_approach_no_nml after fe_lsc_old istat',
cyup     &istat

cyup      write(*,*)'in lsc_approach_no_nml bef allocate delpwr_nc_old_lsc,
cyup     &nrelta,nrayl,nfreq',nrelta,nrayl,nfreq

      allocate( delpwr_nc_old_lsc(1:nrelta,1:nrayl*nfreq),STAT=istat)
cyup      write(*,*)'after allocate  delpwr_nc_old_lsc istat',istat
      call bcast(delpwr_nc_old_lsc,zero,SIZE(delpwr_nc_old_lsc))
cyup      write(*,*)'in lsc_approach_no_nml before end'

c-----for E_DC profile calculations

c     EdcTSC_1D and JparTSC_1D are used in the namelists
c     So,they can be pointers and they are set in lsc_approach_nml.i
c
c      allocate(EdcTSC_1D(1:n_psi_TSC),STAT=istat)
c      call bcast(EdcTSC_1D,zero,SIZE(EdcTSC_1D))
c      write(*,*)'in lsc_approach_no_nml aft EdcTSC_1D istat',
c     &istat 

c      allocate(JparTSC_1D(1:n_psi_TSC),STAT=istat)
c      call bcast(JparTSC_1D,zero,SIZE(JparTSC_1D))
c      write(*,*)'in lsc_approach_no_nml aft JparTSC_1D istat',
c     &istat 

      allocate(j_rf_TSC_1D(1:n_psi_TSC),STAT=istat)
      call bcast(j_rf_TSC_1D,zero,SIZE(j_rf_TSC_1D))
cyup      write(*,*)'in lsc_approach_no_nml aft j_rf_TSC_1D istat',
cyup     &istat 

      allocate(d_j_rf_d_Edc_TSC_1D(1:n_psi_TSC),STAT=istat)
      call bcast(d_j_rf_d_Edc_TSC_1D,zero,SIZE(d_j_rf_d_Edc_TSC_1D))
cyup      write(*,*)'in lsc_approach_no_nml aft d_j_rf_d_Edc_TSC_1D istat',
cyup     &istat 

      allocate(d_ln_j_rf_d_ln_Edc_TSC_1D(1:n_psi_TSC),STAT=istat)
      call bcast(d_ln_j_rf_d_ln_Edc_TSC_1D,zero,
     &SIZE(d_ln_j_rf_d_ln_Edc_TSC_1D))
cyup      write(*,*)'in lsc_approach_no_nml aft d_ln_j_rf_d_ln_Edc_TSC_1D
cyup     & istat',istat

    
      allocate(power_dens_watt_m3_TSC_1D(1:n_psi_TSC),STAT=istat)
      call bcast(power_dens_watt_m3_TSC_1D,zero,
     &SIZE(power_dens_watt_m3_TSC_1D))
cyup      write(*,*)'in lsc_approach_no_nml aft power_densTSC_1D istat',
cyup     &istat

      allocate(CD_dens_no_Edc_a_m2_TSC_1D(1:n_psi_TSC),STAT=istat)
      call bcast(CD_dens_no_Edc_a_m2_TSC_1D,zero,
     &SIZE(CD_dens_no_Edc_a_m2_TSC_1D))
cyup      write(*,*)'in lsc_approach_no_nml aft CD_dens_no_Edc_a_m2_TSC_1D
cyup     &istat',istat

      allocate(CD_dens_small_Edc_a_m2_TSC_1D(1:n_psi_TSC),STAT=istat)
      call bcast(CD_dens_small_Edc_a_m2_TSC_1D,zero,
     &SIZE(CD_dens_small_Edc_a_m2_TSC_1D))
cyup      write(*,*)'in lsc_approach_no_nml aft CD_dens_small_Edc_a_m2_TSC_1D
cyup     &istat',istat

      return
      end

      subroutine ainalloc_fourb_i
c-----allocate pointers in fourb.i
      implicit none
      include 'param.i'
      include 'one_nml.i'
      include 'edge_prof_nml.i'
      include 'fourb.i'
c-----locals
      real*8  zero
      
      integer istat
    
      zero=0.d0

cyup      write(*,*)'in ainalloc_fourb_i'

      allocate(rr_add(1:nxeqd_add),STAT=istat)
      call bcast(rr_add,zero,SIZE(rr_add))
cyup      write(*,*)'in ainalloc_fourb_i after rr_add istat',istat

      allocate(zz_add(1:nyeqd_add),STAT=istat)
      call bcast(zz_add,zero,SIZE(zz_add))
cyup      write(*,*)'in ainalloc_fourb_i after zz_add istat',istat

      allocate(distance_to_wall(1:nxeqd_add,
     &1:nyeqd_add,0:max_limiters),STAT=istat)
      call bcast(distance_to_wall,zero,SIZE(distance_to_wall))
cyup      write(*,*)'in ainalloc_fourb_i after distance_to_wall istat',istat

      allocate(density_r_z(1:nxeqd_add,1:nyeqd_add,1:nbulk,
     &0:max_limiters),STAT=istat)
      call bcast(density_r_z,zero,SIZE(density_r_z))
cyup      write(*,*)'in ainalloc_fourb_i after density_r_z istat',istat

      allocate(density_r_z_rr(1:nxeqd_add,1:nyeqd_add,1:nbulk,
     &0:max_limiters),STAT=istat)
      call bcast(density_r_z_rr,zero,SIZE(density_r_z_rr))
cyup      write(*,*)'in ainalloc_fourb_i after &density_r_z_rr istat',istat

      allocate(density_r_z_zz(1:nxeqd_add,1:nyeqd_add,1:nbulk,
     &0:max_limiters),STAT=istat)
      call bcast(density_r_z_zz,zero,SIZE(density_r_z_zz))
cyup      write(*,*)'in ainalloc_fourb_i after &density_r_z_zz istat',istat

      allocate(density_r_z_rrzz(1:nxeqd_add,1:nyeqd_add,1:nbulk,
     &0:max_limiters),STAT=istat)
      call bcast(density_r_z_rrzz,zero,SIZE(density_r_z_rrzz))
cyup      write(*,*)'in ainalloc_fourb_i aft &density_r_z_rrzz istat',istat

      allocate(temperature_r_z(1:nxeqd_add,1:nyeqd_add,1:nbulk),
     &STAT=istat)
      call bcast(temperature_r_z,zero,SIZE(temperature_r_z))
cyup      write(*,*)'in ainalloc_fourb_i after temperature_r_z istat',istat

      allocate(temperature_r_z_rr(1:nxeqd_add,1:nyeqd_add,1:nbulk),
     &STAT=istat)
      call bcast(temperature_r_z_rr,zero,SIZE(temperature_r_z_rr))
cyup      write(*,*)'in ainalloc_fourb_i aft temperature_r_z_rr istat',istat

      allocate(temperature_r_z_zz(1:nxeqd_add,1:nyeqd_add,1:nbulk),
     &STAT=istat)
      call bcast(temperature_r_z_zz,zero,SIZE(temperature_r_z_zz))
cyup      write(*,*)'in ainalloc_fourb_i aft temperature_r_z_zz istat',istat

      allocate(temperature_r_z_rrzz(1:nxeqd_add,1:nyeqd_add,1:nbulk),
     &STAT=istat)
      call bcast(temperature_r_z_rrzz,zero,SIZE(temperature_r_z_rrzz))
cyup      write(*,*)'in ainalloc_fourb_i after temperature_r_z_rrzz istat',
cyup     &istat

      allocate(zeff_r_z(1:nxeqd_add,1:nyeqd_add),STAT=istat)
      call bcast(zeff_r_z,zero,SIZE(zeff_r_z))
cyup      write(*,*)'in ainalloc_fourb_i aft zeff_r_z istat',istat

      allocate(zeff_r_z_rr(1:nxeqd_add,1:nyeqd_add),STAT=istat)
      call bcast(zeff_r_z_rr,zero,SIZE(zeff_r_z_rr))
cyup      write(*,*)'in ainalloc_fourb_i aft zeff_r_z_rr istat',istat

      allocate(zeff_r_z_zz(1:nxeqd_add,1:nyeqd_add),STAT=istat)
      call bcast(zeff_r_z_zz,zero,SIZE(zeff_r_z_zz))
cyup      write(*,*)'in ainalloc_fourb_i aft zeff_r_z_zz istat',istat

      allocate(zeff_r_z_rrzz(1:nxeqd_add,1:nyeqd_add),STAT=istat)
      call bcast(zeff_r_z_rrzz,zero,SIZE(zeff_r_z_rrzz))
cyup      write(*,*)'in ainalloc_fourb_i aft zeff_r_z_rrzz istat',istat

      return
      end

