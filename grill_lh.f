
c     This subroutine calculates the initial conditions for      
c     the grill type wave launch
c     ***************************************************************
c       input parameters:				               
c
c       igrillpw: option specifying the N_parallel power spectra       
c                =1 > power=powers/nnkpar	(default)
c                =2 > power=sin**2(x)/x**2,			      
c                     x=2pi(n_par-0.5(anmax-anmin))/(anmax-anmin)
c                =3 > power=exp-((n_par-a1)/a2)**2,
c                           and normalized to powers().
c                           a1=anmin(1:ngrill), a2=anmax(1:ngrill)     
c       ngrill   is a number of the poloidal grill angles
c       ngrilla    is a maximal number of the poloidal grill angles
c       thgrill(ngrilla)    poloidal angle of grill,measured counter   
c                          clockwise from horizontal through the        
c                          magnetic axis (degrees)		       
c       height(ngrilla)     is a poloidal length (m) of grill 	       
c                          (giving poloidal power distribution	       
c                          of each grill).    			       
c       nthin(ngrilla)      is a number of ray near the each poloidal
c                          center, simulating a grill                  
c       phigrill(1;ngrilla) is a toroidal grill angle of grill (degrees)
c                n_parallel (ngrill=nspect)
c
c       igrilltw: option specifying poloidal distribution of power
c                =1 > spread with equal weight over height 
c                =2 > {cos(pi(theta_pol-thgrill(i))/(height/radius))}**2
c                     (default).
c                     -0.5height/radius<theta_pol-thgrill(i)<0.5height/radius
c
c       i_grill_pol_mesh: option specifying the poloidal mesh wtheta(j)
c                         near the central grill angle thgrill(i)
c                         =1 equispaced mesh 
c                            wtheta(j)-wtheta(j-1)=zdth=Const (default)
c                         =2 poloidal mesh will be chosen to get the equal
c			     power fpwth(j) for all rays near the central 
c                            grill angle fpwth(j)=1/nthini
c
c       i_grill_npar_ntor_npol_mesh: option specifying the refractive
c                         index meshes.
c
c                         For  i_n_poloidal=1,2,3 it is specifying
c                         n_parallel mesh anzin(n) for the power
c                         spectrum pwcpl(n) n=1,...,nnkpari
c                         =1 equispaced mesh 
c                            anzin(n)-anzin(n-1)=hnpar=Const (default)
c                         =2 n_parallel mesh will be chosen to get the equal
c			     power pwcpl(n) for all rays in the given power
c                            spectrum  pwcpl(n)=1.d0/nnkpari 
c                            pwcpl(n)=power_spectrum(anzin(n))*
c                                     delta_npar_bin(n)= 1.d0/nnkpari
c
c                            For  i_n_poloidal=4 it is specifying two meshes:
c                            a) n_toroidal mesh anztorin(ntor) and             
c                            b) n_poloidal mesh anzpolin(npolmesh) 
c                            for the power spectrum
c                            pwcpl_tp(1:nnktori,1:nnkpoli)=pwcpl_t*pwcpl_t 
c                         =1 equispaced meshes (default)
c                            anztorin(ntor)- anztorin(ntor-1)=hntor=Const 
c                            anzpolin(npol)- anzpolin(npol-1)=hnpol=Const
c                         =2 the meshes anztorin(1:nntori) anzpolin(1:nnkpoli)
c                            will be chosen to get the equal
c			     power pwcpl_tp(ntor,npol) for all rays in 
c                            the given power spectrum 
c                            pwcpl_tp(ntor,npol)=1.d0/(nnktori*nnkpoli) 
c                         
c
c       anmin(1:ngrilla)    position of the left bound of	       
c                          power spectrum P(n_parallel) (Can be neg).
c                          Or, as specified for igrillpw=3 above.
c                          It needs for i_n_poloidal=1,2,3  
c       anmax(1:ngrilla)    position of the right bounds  	       
c                          of power spectrum P(n_parallel)	       
c                          Or, as specified for igrillpw=3 above.
c                          It needs for i_n_poloidal=1,2,3   
c       nnkpar(1:ngrilla)   number of points  of power spectrum	       
c                          P(n_parallel)		
c                          It needs for i_n_poloidal=1,2,3 	       
c       powers(1:ngrilla)   power in one grill (MWatts)	    	       
c                          (total power of grill(in MWatts) will be    
c                          powtotlh=sum{powers}			       
c       rhopsi0(1:ngrilla)  initial psi for wave front (0<rhopsi0<1)
c	xma (m)	           major radius of magnetic axis	       
c       yma (m)            vertical shift of magnetic axis
c       psimag	       
c       nnkprmax=max{i=1,ngrill}nnkpar(i)
c
c       i_n_poloidal gives the type of the grill launch	
c                    =1,2,3 N_parallel and other variables will be given
c                    =4     N_toroidal and N_poloidal will be given
c
c       The following variables work for i_n_poloidal=4
c
c       antormin(1:ngrilla)    position of the left bound of	       
c                          power spectrum P(n_parallel) (Can be neg).
c                          Or, as specified for igrillpw=3 above.  
c       antormax(1:ngrilla)    position of the right bounds  	       
c                          of power spectrum P(n_parallel)	       
c                          Or, as specified for igrillpw=3 above.  
c       nnktor(1:ngrilla)   number of points  of power spectrum	       
c                          P(n_parallel)	
c       anpolmin(1:ngrilla)    position of the left bound of	       
c                          power spectrum P(n_parallel) (Can be neg).
c                          Or, as specified for igrillpw=3 above.  
c       anpolmax(1:ngrilla)    position of the right bounds  	       
c                          of power spectrum P(n_parallel)	       
c                          Or, as specified for igrillpw=3 above.  
c       nnkpol(1:ngrilla)   number of points  of power spectrum	       
c                          P(n_parallel)	
c       nnktormax=max{i=1,ngrill}nnktor(i)
c       nnkpolmax=max{i=1,ngrill}nnkpol(i)
c
c       n_theta_pol   the poloidal refractive index
c                     used in i_n_poloidal=2 case 
c--------------------------------------------------------------------
c       output parameters:					       
c       nray               total number of the rays 		       
c                          nray=Sum{i=1,ngrill}Sum_{j=1,nthin(i)}      
c                          Sum{n=1,nnkpar(i)}nnkpar(n)                 
c       arzu0(nray)   (m)  initial positions of	                       
c       arru0(nray)   (m)                  ray_iray                    
c       arphiu0(nray) (radian)             start points		       
c       arntheta(nray)     initial  poloidal refractive indexes	       
c       arnphi(nray)	   initial  toroidal refractive indexes        
c       powinilh(nray)    initial power flowing in the iray wave channel
c		           normalized Sum_{i=1,nray}powinilh(i)=powtotlh
c                          (erg/c)				       
c       fpwth(nthinmax)    spectrum poloidal distributions             
c                          normalized 1=Sum{j=1,nthin}fpwth(j)
c       wtheta(nthinmax)   poloidal angle mesh wtheta(j) radian        
c                          near the grill poloidal angles thgrill(i)   
c                          j=1,nthin(i)
c       wtheta_edge(nthinmax+1): wtheta_edge(j) j=1,,nthin=1
c                          wtheta(j)=0.5*
c                          (wtheta_edge(j)+wtheta_edge(j+1)),j=1,nthin
c		       
c       anzin(nnkprmax)    are the points of mesh n_parallel	       
c                          P(n_parallel)  anzin(j), j=1,nnkpar(i)
c
c       anzin_edge(nnkprmax+1) : anzin_edge(j), j=1,nnkpar(i)+1
c                         anzin(j)=0.5(anzin_edge(j)+anzin_edge(j+1))
c          	          j=1,nnkpar(i)
c
c       pwcpl(nnkprmax)    power spectrum pwcpl(n) on   	       
c                          n_parallel, normalization		       
c                          Sum{n=1,nnkpar(i)}pwspl(n)=powers(i)        
c       wdnpar0(nray)	   initial values of the n_parallel width of   
c                          the rays
c  
c       pwcpl_tp(nnktormax,nnkpolmax)    power spectrum pwcpl on
c                          (N_toroidal,N_poloidal) mesh parallel,
c                          It is for i_n_poloidal=4
c       normalization		       
c       Sum{ntor=1,nnktor(i),npol=1,nnkpol(i)}pwspl_tp(ntpr,npol)
c               =powers(i)  
c
c       powtotlh
c----------------------------------------------------------------
c       for i_n_poloidal=4 case
c       anztorin(nnktormax)   are the points of mesh n_toroidal       
c                          P(n_toroidal,n_ploidal), j=1:nnktori
c       anztpolin(nnkpolmax)   are the points of mesh n_poloidal
c                              j=1:nnkpoli
c       anztorin_edge(nnktormax+1) N_toroidal mesh, j=1:nnktori+1 
c                anztorin(j)=0.5(anztorin_edge(j+1)+anztorin_edge(j))
c
c       anzpolin_edge(nnkpolmax+1) N_poloidal mesh, j=1:nnkpoli+1
c                anztpolin(j)=0.5(anzpolin_edge(j+1)+anzpolin_edge(j))
c
c       pwcpl_t_ar(nnktormax), power spectrum on n_toroidal mesh 
c       pwcpl_p_ar(nnkpolmax), power spectrum on n_poloidal mesh 
c       pwcpl_tp(nnktormax,nnkpolmax)    power spectrum on   	       
c                          n_toroidal,n_poloidal mesh, 
c       normalization
c       Sum{ntor=1,nnktor(i),n_pol=1,nnkpol(i)}(pwspl_tp(ntor,npol)
c                                 =powers(i)
c     	       
c       ilaunch=1 gives explicit r0launch,phi0launch,z0launch launch
c                 location (meters and radians). This is for test cases.
cSAP111008
c       ilaunch=2   gives explicit rlaunch(nthin(igrill),igrill),
c                   zlaunch(nthin(igrill),igrill),
c                   phigrill(igrill),
c-----------------------------------------------------------------------
c  uses: double precision functions r(psi,theta),z(psi,theta),psiff(rho)
c-----------------------------------------------------------------------

      subroutine grill_lh(rhopsi0,ngrilla,ngrill,thgrill,phigrill,
     1 height,nthin,nthinmax,
     1 anmin,anmax,nnkpar,powers,powtotlh,
     1 antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1 n_theta_pol,
     1 xma,yma,psimag,
     1 fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,nnkprmax,
     1 anztorin,anzpolin,pwcpl_tp,nnktormax,nnkpolmax,
     1 anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1 nray,arzu0,arru0,arphiu0,arntheta,arnphi,powinilh,
     1 nraymaxl,wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1 ilaunch,r0launch,phi0launch,z0launch,
cSAP111030
     1 rlaunch,zlaunch,weight_power_launch,
     1 i_grill_pol_mesh,
     1 i_grill_npar_ntor_npol_mesh)

      implicit none
      integer nraymaxl,ngrilla,ngrill,nray,iray,nthini,nnkpari,
     &nthinmax,
     .nnkprmax,igrillpw,igrilltw,i_n_poloidal,ilaunch,
     .i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh,
     &nnktormax,nnkpolmax
      double precision thgrill(ngrilla),height(ngrilla),
     1          phigrill(ngrilla),rhopsi0(ngrilla),
     1          anmin(ngrilla),anmax(ngrilla),powers(ngrilla),
     1          antormin(ngrilla),antormax(ngrilla),
     1          anpolmin(ngrilla),anpolmax(ngrilla),
     1          n_theta_pol
      double precision r0launch,phi0launch,z0launch,
cSAP111030
     1          rlaunch(nthinmax,ngrilla),
     1          zlaunch(nthinmax,ngrilla),
     1          weight_power_launch(nthinmax,ngrilla)
 
      integer   nnkpar(ngrilla),nthin(ngrilla),
     & nnktor(ngrilla),nnkpol(ngrilla)
      double precision arzu0(*),arru0(*),arphiu0(*),arntheta(*),
     1       arnphi(*),powinilh(*)
      double precision fpwth(nthinmax),wtheta(nthinmax),
cSm050309
     1       wtheta_edge(nthinmax+1)
      double precision anzin(nnkprmax),anzin_edge(nnkprmax+1),
     &pwcpl(nnkprmax)

      double precision anztorin(nnktormax),anzpolin(nnkpolmax),
     &pwcpl_tp(nnktormax,nnkpolmax),
     &anztorin_edge(nnktormax+1),anzpolin_edge(nnkpolmax+1),
     &pwcpl_t_ar(nnktormax),pwcpl_p_ar(nnkpolmax)

      double precision xma,yma,psimag,powtotlh,pi,tet,psi0,z,r,
     .rhomag,thaper,zdth,delth,zcos,znorm,an0,difrnpar,andwdth,
     .hnpar,anorm,delnpar,x1,xx,phi0,theta0,z0,r0,cnphi,cnteta
      double precision psi_rho,powert
c-----for i_n_poloidal=4 case, input N_toroidal and N_poloidal
      double precision antor0,anpol0,difrntor,difrnpol,
     &andwdth_tor,andwdth_pol,hntor,hnpol,hnwidth,delntor,delnpol,
     &pwcpl_t,pwcpl_p,y1,yy,anpar
      integer nnktori,nnkpoli


      integer i,j,n,ntor,npol,k

      double precision wdnpar0(*)

cSm050309
      double precision f_pow_poloid_1,f_pow_poloid_0
      external f_pow_poloid_1,f_pow_poloid_0

c-----external for test
      real*8 rhopsi,fpsi,fpsi_

      double precision aperture
      common /aperture/ aperture
      double precision f_pow_npar_2,f_pow_npar_3
      external f_pow_npar_2,f_pow_npar_3
      
      if (ilaunch.eq.1) then
         write(*,*)
         write(*,*)'grill_lh: ilaunch=1; ngrill,nthin(1),'
         write(*,*)'grill_lh: nnkpar(1) have been set =1'
         write(*,*)
      endif

cSAP111030 
      if (ilaunch.eq.2) then
         write(*,*)
         write(*,*)'grill_lh: ilaunch=2:ngrill=',ngrill
         write(*,*)'nthinmax',nthinmax
         do i=1,ngrill
            write(*,*)'grill_lh:i,nthin(i)',i,nthin(i)
            if (nthin(i).gt.nthinmax) then
               write(*,*)'grill_lh: nthin(i).gt.nthinmax'
               write(*,*)'It should be nthin(i).le.nthinmax'
               write(*,*)'Please change in param.i file'
               write(*,*)'and recompile the code' 
               stop 'in grill_lh'
            endif
         enddo !i
      endif
     
      write(*,*)'in grill_lh ngrill,i_n_poloidal,ilaunch',
     &                       ngrill,i_n_poloidal,ilaunch

      if (ngrill.lt.1) then 
        write(*,*)'in ngrill_lh ngrill<1 but it should be .ge.1'
        write(*,*)'Please change ngrill in genray.in file' 
        stop
      endif

      if (ngrill.gt.ngrilla) then
        write(*,*)'in ngrill_lh ngrilla<ngrill but it should be .ge.'
        write(*,*)'please change ngrilla in param.i and recompile'
        stop
      endif
   
      do i=1,ngrill
cSAP111008 
         if(ilaunch.eq.0)then
           write(*,*)'i,thgrill(i),nthin(i)',i,thgrill(i),nthin(i)
           write(*,*)'height(i),phigrill(i),rhopsi0(i),powers(i)',
     1                height(i),phigrill(i),rhopsi0(i),powers(i)
         endif

         if(ilaunch.eq.2)then
           write(*,*)'grill number i=',i,'nthin(i)=',nthin(i)
           write(*,*)'phigrill(i)',phigrill(i)

           do n=1,nthin(i)  
              write(*,*)'i,n,rlaunch(n,i),zlaunch(n,i)',
     &                   i,n,rlaunch(n,i),zlaunch(n,i)
              write(*,*)'weight_power_launch(n,i)',
     &                   weight_power_launch(n,i)
           enddo
         endif
 
         if ((i_n_poloidal.ge.1).and.(i_n_poloidal.le.3)) then
           write(*,*)'anmin(i),anmax(i),nnkpar(i)',
     1                anmin(i),anmax(i),nnkpar(i)
         else !i_n_poloidal=4
           write(*,*)'antormin(i),antormax(i),nnktor(i)',
     1               antormin(i),antormax(i),nnktor(i)
           write(*,*)'anpolmin(i),anpolmax(i),nnkpol(i)',
     1                anpolmin(i),anpolmax(i),nnkpol(i)
         endif
      enddo !i
c-----------------------------------------------------------
c     powtotlh is in MWatts
c     Transformation of the total antenna power from MW to erg/sec
c-----------------------------------------------------------------
      powtotlh=0.d0
      do i=1,ngrill
CAYP201124 this was bug         powers(i)=powers(i)*1.d13   !erg/sec  

      !Originally, it was envisioned this subroutine would be called
      !   only once, at each execution. OXB interation (i_ox=1) of
      !   could lead to many calls (although power not relevant
      !   for this case).
      !YuP[2020-11-25] Causing overflow, if many calls?
      !YuP[2020-11-25] However, array powers() is used further 
      !in this subroutine, in erg/sec units.
      !If we don't rescale powers() here, 
      ! it should be done at all those lines below.
      !YuP adjusted all those lines below.
         powtotlh=powtotlh+powers(i)*1.d13
      enddo

      write(*,*)'powtotlh',powtotlh
c----------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      tet=pi/180.d0
c-------------------------------------------------------------------
cSAP111030
      if (ilaunch.eq.0) then
        do i=1,ngrill
           thgrill(i)=thgrill(i)*tet
           phigrill(i)=phigrill(i)*tet
        end do
      endif

      if (ilaunch.eq.2) then
        do i=1,ngrill
           phigrill(i)=phigrill(i)*tet
        end do
      endif

c----------------------------------------------------------------
c     For each i, calculation of
c     wtheta (j)  -poloidal mesh, j=1,nthin(i)
c     fpwth  (j)  -poloidal power distributions normalized to
c                    1=Sum{j=1,nthin(i)}fpwth(j)
c----------------------------------------------------------------

      iray=0      !Ray counter

      write(*,*)'grill_lh iray=0, ngrill',ngrill

      do i=1,ngrill
cSAP111030
         if (ilaunch.eq.0) then
            psi0=psi_rho(rhopsi0(i))
           write(*,*)'in grill_lh i,rhopsi0(i),psi0,thgrill(i)',
     *                            i,rhopsi0(i),psi0,thgrill(i)
         endif 
         nthini=nthin(i)
         write(*,*)'nthini',nthini

cSm050309
         if(nthini.eq.1) then
            wtheta(1)=thgrill(i)
            fpwth(1)=1.d0
	 else
c------------------------------------------------------------
c           calculations: r=r(psi0,thgrill(i)),
c                         z=z(psi0,thgrill(i))
c------------------------------------------------------------
cSAP111030
            if (ilaunch.eq.0) then
               call zr_psith(psi0,thgrill(i),z,r)
c	       write(*,*)'grill_lh after zr_psith z,r',z,r 
c------------------------------------------------------------
	       rhomag=dsqrt((r-xma)**2+(z-yma)**2)
	       thaper=height(i)/rhomag
cSm050311
               aperture=thaper !put the data to common /aperture/

cSm050309
               if((i_grill_pol_mesh.eq.1).or.(igrilltw.eq.1))then
c------------------------------------------------------------
c                equispaced poloidal mesh 
c                wtheta(j)-wtheta(j-1)=zdth=Const(default)
c------------------------------------------------------------
                 zdth=thaper/dfloat(nthini)

                 do j=1,nthini
	            delth=-.5d0*thaper+zdth*(j-0.5)
	            wtheta(j)=thgrill(i)+delth
	            zcos=dcos(pi*delth/thaper)
                
c ----------------------------------------------------
c                   fpwth(j) is the poloidal power distribution on
c                       poloidal angle mesh wtheta (j)
c                       near the grill poloidal angle thgrill(i)
c----------------------------------------------------               
                    if (igrilltw.eq.1) then
                       fpwth(j)=1.0
                    else
                       fpwth(j)=zcos*zcos
                    endif
c                   write(*,*)'fpwth(j)',fpwth(j)
                 enddo !j
               endif !i_grill_pol_mesh.eq.1
cSm050309
               if((i_grill_pol_mesh.eq.2).and.(igrilltw.ne.1))then
c-------------------------------------------------------------
c                poloidal mesh will be chosen to get the equal
c                power fpwth(j) for all rays near the central 
c                grill angle fpwth(j)=1/nthini
c-------------------------------------------------------------
        
c_old_method using   f_pow_poloid_1,        
c                call create_equal_mesh(-0.5*thaper,0.5*thaper,
c     &          f_pow_poloid_1,nthini,
c     &          wtheta_edge,wtheta,fpwth)

c                do k=1,nthini
c                   wtheta(k)= wtheta(k)+thgrill(i)
c                enddo ! k

c                do k=1,nthini+1
c                   wtheta_edge(k)= wtheta_edge(k)+thgrill(i)
c                enddo ! k

c                do k=1,nthini
c                   write(*,*)'!1 k,wtheta(k),wtheta_edge(k),fpwth(k)',
c     &                           k,wtheta(k),wtheta_edge(k),fpwth(k)
c                enddo
c                write(*,*)'1 wtheta_edge(nthini+1)',
c     &                       wtheta_edge(nthini+1)

c_new method using f_pow_poloid_0
                 call create_equal_mesh(-0.5*pi,0.5*pi,
     &           f_pow_poloid_0,nthini,
     &           wtheta_edge,wtheta,fpwth)

                 do k=1,nthini
                    wtheta(k)= wtheta(k)*aperture/pi+thgrill(i)
                 enddo ! k
            
                 do k=1,nthini+1
                   wtheta_edge(k)=wtheta_edge(k)*aperture/pi+thgrill(i)
                 enddo ! k

c                do k=1,nthini
c                   write(*,*)'!2 k,wtheta(k),wtheta_edge(k),fpwth(k)',
c     &                           k,wtheta(k),wtheta_edge(k),fpwth(k)
c                enddo
c                write(*,*)'!2 wtheta_edge(nthini+1)',
c    &                         wtheta_edge(nthini+1)

               endif !i_grill_pol_mesh.eq.2
            endif !ilaunch.eq.0

cSAP111030
            if (ilaunch.eq.2) then
              do j=1,nthini
                 fpwth(j)=weight_power_launch(j,i)
              enddo
            endif !ilaunch.eq.2

	 endif ! (nthini.eq.1)
	 
c--------------------------------------------
c        normalization of fpwth
	 znorm=0.d0
	 do j=1,nthini
	    znorm=znorm+fpwth(j)
	 enddo !j

	 do j=1,nthini
	    fpwth(j)=fpwth(j)/znorm
	    write(*,*)'j,fpwth(j),wtheta(j)',j,fpwth(j),wtheta(j)
	 enddo !j
        
         if((i_n_poloidal.ge.1).and.(i_n_poloidal.le.3)) then
c-----------------------------------------------------------------------
c          calculation of anzin(n), the n_parallel mesh,
c          and the corresponding 
c          power spectrum pwcpl(n).
c      	   Normalization:  Sum{n=1,nnkpar(i)}pwspl(n)=powers(i)
c          BH020826:  Added igrillpw.eq.3 option, and slightly
c                     modified n_par grid of points.
c-----------------------------------------------------------------------
          
           if(igrillpw.eq.1 .or. igrillpw.eq.2) then
             an0=0.5d0*(anmin(i)+anmax(i)) !central value of N_parallel
	     difrnpar=anmax(i)-anmin(i)
           elseif (igrillpw.eq.3) then
             an0=anmin(i)
	     difrnpar=2.*anmax(i)  ! power out to exp(-4.)
           endif
            
	   nnkpari=nnkpar(i)
	   write(*,*)'nnkpari',nnkpari
	   if(nnkpari.eq.1)then
	     anzin(1)=an0
	     pwcpl(1)=powers(i)*1.d13 !YuP[2020-11-25] to erg/sec
c	     write(*,*)'nnkpari',nnkpari,'anzin(1),pwcpl(1)',
c     1                                   anzin(1),pwcpl(1)
	     goto 10
	   else
ccBH020826    difrnpar=anmax(i)-anmin(i)
	     andwdth=2.d0*pi/difrnpar
            
             if(i_grill_npar_ntor_npol_mesh.eq.1) then
c---------------------------------------------------------
c              equispaced mesh 
c              anzin(n)-anzin(n-1)=hnpar=Const (default)
c--------------------------------------------------------
cBH020826      hnpar=difrnpar/dfloat(nnkpari+1)
	       hnpar=difrnpar/dfloat(nnkpari)

	       do n=1,nnkpari
cBH020826        anzin(n)=anmin(i)+n*hnpar
                 anzin(n)=an0-0.5*difrnpar+(n-0.5)*hnpar
	       enddo

	       anorm=0.d0
	       do n=1,nnkpari
                 delnpar=anzin(n)-an0
	         x1=andwdth*delnpar
	         xx=x1*x1
	         if(delnpar.ne.0.d0) then
                   if(igrillpw.eq.1) pwcpl(n)=1.d0/nnkpari
    	           if(igrillpw.eq.2) pwcpl(n)=dsin(x1)*dsin(x1)/xx
                  if(igrillpw.eq.3) pwcpl(n)=exp(-(delnpar/anmax(i))**2)
	         else
	           pwcpl(n)=1.d0
	         endif
	         anorm=anorm+pwcpl(n)
	       enddo !n
c              write(*,*)'anorm='anorm

	       do n=1,nnkpari
	         pwcpl(n)=pwcpl(n)/anorm*powers(i)*1.d13 !YuP[2020-11-25] to erg/sec
	         write(*,*)'n,pwcpl(n)',n,pwcpl(n)
	       enddo !n
             endif !i_grill_npar_ntor_npol_mesh.eq.1
            
             if(i_grill_npar_ntor_npol_mesh.eq.2) then
c--------------------------------------------------------------------
c               n_parallel mesh will be chosen to get the equal
c		power pwcpl(n) for all rays in the given power
c               spectrum  pwcpl(n)=1.d0/nnkpari 
c               pwcpl(n)=power_spectrum(anzin(n))*
c                        delta_npar_bin(n)= 1.d0/nnkpari
c---------------------------------------------------------------------
                if (igrillpw.eq.2) then
                  andwdth=2.d0*pi/difrnpar
                  call create_equal_mesh(-pi,pi,f_pow_npar_2,nnkpari,
     &            anzin_edge,anzin,pwcpl)

                  do k=1,nnkpari 
                    anzin(k)=anzin(k)/andwdth+an0 
                    write(*,*)'k, anzin(k)',k, anzin(k)                  
                  enddo 

                  do k=1,nnkpari+1                   
                    anzin_edge(k)=anzin_edge(k)/andwdth+an0
                    write(*,*)'k, anzin_edge(k)',k,anzin_edge(k)
                  enddo

ctest 
c                  do k=1,nnkpari                   
c               write(*,*)'k,anzin(k),0.5(anzin_edge(k)+anzin_edge(k+1))'
c     &             ,k,anzin(k),0.5d0*(anzin_edge(k)+anzin_edge(k+1))
c                  enddo                
cendtes          
                endif !igrillpw.eq.2

                if(igrillpw.eq.3) then 
                  
                  call create_equal_mesh(-1.d0,1.d0,f_pow_npar_3,
     &            nnkpari,anzin_edge,anzin,pwcpl)

                  do k=1,nnkpari 
                    anzin(k)=anzin(k)*0.5d0*difrnpar+an0 
                    write(*,*)'k,anzin(k)',k,anzin(k)
                  enddo
 
                  do k=1,nnkpari+1 
                    anzin_edge(k)=anzin_edge(k)*0.5d0*difrnpar+an0 
                    write(*,*)'k,anzin_edge(k)',k,anzin_edge(k)
                  enddo
 
                endif !igrillpw.eq.3

c------------------------------------------------------------------------
c               power normalization
                anorm=0.d0
	        do n=1,nnkpari
	           anorm=anorm+pwcpl(n)
	        enddo !n
c               write(*,*)'anorm='anorm

	        do n=1,nnkpari
	         pwcpl(n)=pwcpl(n)/anorm*powers(i)*1.d13 !YuP[2020-11-25] to erg/sec
	         write(*,*)'n,pwcpl(n)',n,pwcpl(n)
	        enddo !n
 
            endif !i_grill_npar_ntor_npol_mesh.eq.2

	   endif ! nnkpari.eq.1
10         continue
c--------------------------------------------------------------------
         endif !i_n_poloidal=1,2,3
         
         if(i_n_poloidal.eq.4) then
c-----------------------------------------------------------------------
c          calculation of anztorin(ntor), the N_toroidal mesh,
c          calculation of anzpolin(ntor), the N_poloidal mesh,
c
c          and the corresponding 
c          power spectrum pwcpl_tp(ntor,npol).
c      	   Normalization: 
c           Sum{ntor=1,nnktor(i),npol=1,nnkpol(i)}pwspl_tp(ntor,npol)
c                   =powers(i)
c          BH020826:  Added igrillpw.eq.3 option, and slightly
c                     modified n_par grid of points.
c-----------------------------------------------------------------------

           write(*,*)'i_n_poloidal.eq.4 igrillpw',igrillpw

           if((igrillpw.eq.1).or.(igrillpw.eq.2)) then
             antor0=0.5d0*(antormin(i)+antormax(i)) 
                                           !central value of N_toroidal
             anpol0=0.5d0*(anpolmin(i)+anpolmax(i)) 
                                           !central value of N_poloidal
	     difrntor=antormax(i)-antormin(i)
             difrnpol=anpolmax(i)-anpolmin(i)

             write(*,*)'i,antormax(i),antormin(i),antor0,difrntor',
     &                  i,antormax(i),antormin(i),antor0,difrntor
             write(*,*)'i,anpolmax(i),anpolmin(i),anpol0,difrnpol',
     &                  i,anpolmax(i),anpolmin(i),anpol0,difrnpol
             

           elseif (igrillpw.eq.3) then
             antor0=antormin(i)
             anpol0=anpolmin(i)
	     difrntor=2.*antormax(i)  ! power out to exp(-4.)
             difrnpol=2.*anpolmax(i)  ! power out to exp(-4.)
           endif

	   nnktori=nnktor(i)
           nnkpoli=nnkpol(i)       

	   write(*,*)'nnktori,nnkpoli',nnktori,nnkpoli
	   if((nnktori.eq.1).and.(nnkpoli.eq.1))then
	     anztorin(1)=antor0
             anzpolin(1)=anpol0
 	     pwcpl_tp(1,1)=powers(i)*1.d13 !YuP[2020-11-25] to erg/sec

	     write(*,*)'nnktori,nnkpoli',nnktori,nnkpoli

             write(*,*)'anztorin(1),anzpolin(1),pwcpl_tp(1,1)',
     &       anztorin(1),anzpolin(1),pwcpl_tp(1,1)

c            Save composite width for nominal npar_width, wdnpar:
             hnwidth=0.05*sqrt(antor0**2+anpol0**2)

	     goto 20
	   else

	     andwdth_tor=2.d0*pi/difrntor
             andwdth_pol=2.d0*pi/difrnpol
 
             if(i_grill_npar_ntor_npol_mesh.eq.1) then
c---------------------------------------------------------
c              equispaced meshs 
c              anztorin(n)-anztorin(n-1)=hntor=Const (default)
c              anzpolin(n)-anzpolin(n-1)=hnpol=Const (default)
c---------------------------------------------------------
	       hntor=difrntor/dfloat(nnktori)
               hnpol=difrnpol/dfloat(nnkpoli)
             
	       do ntor=1,nnktori
                 anztorin(ntor)=antor0-0.5*difrntor+(ntor-0.5)*hntor
	       enddo

	       do npol=1,nnkpoli
                 anzpolin(npol)=anpol0-0.5*difrnpol+(npol-0.5)*hnpol
	       enddo

c              Save composite width for nominal npar_width, wdnpar:
               hnwidth=sqrt(hntor**2+hnpol**2)

	       anorm=0.d0

	       do ntor=1,nnktori
                 delntor=anztorin(ntor)-antor0
	         x1=andwdth_tor*delntor
	         xx=x1*x1
                 do npol=1,nnkpoli
                   delnpol=anzpolin(npol)-anpol0
cShiraiwa15113	           y1=andwdth_tor*delntor
	           y1=andwdth_pol*delntor
	           yy=y1*y1
	           if(delntor.ne.0.d0) then
                     if(igrillpw.eq.1) pwcpl_t=1.d0/nnktori
    	             if(igrillpw.eq.2) pwcpl_t=dsin(x1)*dsin(x1)/xx
               if(igrillpw.eq.3) pwcpl_t=exp(-(delntor/antormax(i))**2)
	           else
	              pwcpl_t=1.d0
	           endif

	           if(delnpol.ne.0.d0) then
cShiraiwa15113                    if(igrillpw.eq.1) pwcpl_p=1.d0/nnktori
                     if(igrillpw.eq.1) pwcpl_p=1.d0/nnkpoli
    	             if(igrillpw.eq.2) pwcpl_p=dsin(y1)*dsin(y1)/yy
               if(igrillpw.eq.3) pwcpl_p=exp(-(delnpol/anpolmax(i))**2)
	           else
cShiraiwa151105	             pwcpl_t=1.d0
	             pwcpl_p=1.d0
	           endif
	           write(*,*)'pwcpl_p, pwcpl_t',
     &                        pwcpl_p, pwcpl_t 
cBH151110                   pwcpl_tp(ntor,npol)=pwcpl_t*pwcpl_t
                   pwcpl_tp(ntor,npol)=pwcpl_t*pwcpl_p

	           anorm=anorm+pwcpl_tp(ntor,npol)
	         enddo !npol
               enddo ! ntor
c              write(*,*)'anorm='anorm

	       do ntor=1,nnktori
                 do npol=1,nnkpoli
	         pwcpl_tp(ntor,npol)=pwcpl_tp(ntor,npol)/anorm*powers(i)*1.d13 !YuP[2020-11-25] to erg/sec
	           write(*,*)'ntor,npol,pwcpl_tp(ntor,npol)',
     &                    ntor,npol,pwcpl_tp(ntor,npol)
                  enddo !npol
	       enddo !ntor
             endif ! i_grill_npar_ntor_npol_mesh.eq.1

             if(i_grill_npar_ntor_npol_mesh.eq.2) then
c------------------------------------------------------------------
c               n_toroidal mesh anztorin() will be chosen to get the equal
c		power pwcpl_t_ar(n) for all rays in the given power
c
c               n_poloidal mesh anzpolin() will be chosen to get the equal
c		power pwcpl_p_ar(n) for all rays in the given power
c
c               two dimentional array pwcpl_tp(nnktormax,nnkpolmax)
c               for power distribution on (N_toroidal,N_poloidal)
c               In this case it will pwcpl_tp(i,j)=Const for all rays
c------------------------------------------------------------------
                if(igrillpw.eq.2) then

c-----------------N_toroidal mesh calculations

                  call create_equal_mesh(-pi,pi,f_pow_npar_2,nnktori,
     &            anztorin_edge,anztorin,pwcpl_t_ar)

                  write(*,*)'andwdth_tor,antor0',andwdth_tor,antor0
                  do k=1,nnktori 
                    anztorin(k)=anztorin(k)/andwdth_tor+antor0
                    write(*,*)'k,anztorin(k)',k,anztorin(k)                  
                  enddo 

                  do k=1,nnktori+1                   
                     anztorin_edge(k)=anztorin_edge(k)/andwdth_tor+
     &                                 antor0
                    write(*,*)'k,anztorin_edge(k)',k,anztorin_edge(k)
                  enddo

c-----------------N_poloidal mesh calculations

                  call create_equal_mesh(-pi,pi,f_pow_npar_2,nnkpoli,
     &            anzpolin_edge,anzpolin,pwcpl_p_ar)

                  write(*,*)'andwdth_pol,anpol0',andwdth_pol,anpol0
                  do k=1,nnkpoli 
                    anzpolin(k)=anzpolin(k)/andwdth_pol+anpol0
                    write(*,*)'k,anzpolin(k)',k,anzpolin(k)                  
                  enddo 

                  do k=1,nnkpoli+1                   
                     anzpolin_edge(k)=anzpolin_edge(k)/andwdth_pol+
     &                                anpol0
                    write(*,*)'k,anzpolin_edge(k)',k,anzpolin_edge(k)
                  enddo

                endif ! igrillpw.eq.2

                if(igrillpw.eq.3) then

c-----------------N_toroidal mesh calculations

                  call create_equal_mesh(-1.d0,1.d0,f_pow_npar_3,nnktori
     &            ,anztorin_edge,anztorin,pwcpl_t_ar)

                  do k=1,nnktori 
                    anztorin(k)=anztorin(k)*0.5d0*difrntor+antor0
                    write(*,*)'k,anztorin(k)',k,anztorin(k)                  
                  enddo 

                  do k=1,nnktori+1                   
                     anztorin_edge(k)=anztorin_edge(k)*0.5d0*difrntor+
     &                                antor0
                    write(*,*)'k,anztorin_edge(k)',k,anztorin_edge(k)
                  enddo

c-----------------N_poloidal mesh calculations

                  call create_equal_mesh(-1.d0,1.d0,f_pow_npar_3,nnkpoli
     &            ,anzpolin_edge,anzpolin,pwcpl_p_ar)

                  do k=1,nnkpoli 
                    anzpolin(k)=anzpolin(k)*0.5d0*difrnpol+anpol0
                    write(*,*)'k,anzpolin(k)',k,anzpolin(k)                  
                  enddo 

                  do k=1,nnkpoli+1                   
                     anzpolin_edge(k)=anzpolin_edge(k)*0.5d0*difrnpol+
     &                                anpol0
                    write(*,*)'k,anzpolin_edge(k)',k,anzpolin_edge(k)
                  enddo

                endif ! igrillpw.eq.3

                do ntor=1,nnktori
                   pwcpl_t_ar(ntor)=1.d0/nnktori
                enddo

                do npol=1,nnkpoli
                   pwcpl_p_ar(npol)=1.d0/nnkpoli
                enddo

                do ntor=1,nnktori
                  do npol=1,nnkpoli
                    pwcpl_tp(ntor,npol)=pwcpl_t_ar(ntor)*
     &                                  pwcpl_p_ar(npol)
                  enddo
                enddo

                anorm=0.d0
 
                do ntor=1,nnktori
                  do npol=1,nnkpoli
                    anorm=anorm+pwcpl_tp(ntor,npol)                    
                  enddo
                enddo
c               write(*,*)'anorm='anorm

                do ntor=1,nnktori
                  do npol=1,nnkpoli
	         pwcpl_tp(ntor,npol)=pwcpl_tp(ntor,npol)/anorm*powers(i)*1.d13 !YuP[2020-11-25] to erg/sec
	           write(*,*)'ntor,npol,pwcpl_tp(ntor,npol)',
     &                        ntor,npol,pwcpl_tp(ntor,npol)
                  enddo !npol
	        enddo !ntor

             endif ! i_grill_npar_ntor_npol_mesh.eq.2

	   endif
20         continue
         endif !i_n_poloidal=4

c        calculation of the intial data for rays_iray
	 phi0=phigrill(i)
c	 write(*,*)'phi0=',phi0
	 do j=1,nthini
cSAP111030
            if (ilaunch.eq.0) then
               theta0=wtheta(j)
               write(*,*)'j,psi0,theta0',j,psi0,theta0
c------------------------------------------------------------
c              calculates  z0=z(psi0,theta0) and r0
c------------------------------------------------------------
               call zr_psith(psi0,theta0,z0,r0)
              write(*,*)'grill_lh after call zr_psith z0,r0',z0,r0
c------------------------------------------------------------------
            endif
cSAP111030
            if (ilaunch.eq.1) then
               r0=r0launch
               z0=z0launch
               phigrill(1)=phi0launch*tet
            endif
cSAP111030
            if (ilaunch.eq.2) then
               r0=rlaunch(j,i)
               z0=zlaunch(j,i)
               phi0=phigrill(i)
            endif

            if((i_n_poloidal.ge.1).and.(i_n_poloidal.le.3)) then 
c-------------i_n_poloidal=1,2,3
	      do n=1,nnkpari
	        write(*,*)'grill_lh: 1 n=',n
	        iray=iray+1
                write(*,*)'grill_lh: iray',iray                
	        if (iray.gt.nraymaxl) then
	         write(*,*)'in grill_lh iray=',iray,'nraymaxl',nraymaxl
	         write(*,*)'in grill_lh the total number of the rays is'
		 write(*,*)'grater than number of elements in arrays'
		 write(*,*)'iray.gt.nraymaxl'
		 write(*,*)'please change parameter nraymaxl in param.i'
		 stop
	        endif

	        arzu0(iray)=z0
	        arru0(iray)=r0
	        arphiu0(iray)=phigrill(i)

c-----------------------------------------------------
cSm050826
               if((i_n_poloidal.eq.1).or.(i_n_poloidal.eq.3)) then
c----------------calculations of N_theta and N_phi from N_par
c
c                For (i_n_poloidal.eq.3) case cnteta and cnphi are 
c                the componetnts of the vector N_parallel
c                They are not the components of the total vector N.  
            
                 call nphiteta(z0,r0,phigrill(i),anzin(n),cnteta,cnphi)

               else
                 if(i_n_poloidal.eq.2) then
                  cnteta=n_theta_pol
c-----------------calculation on N_phi from N_par and N_theta 
                  call grill_i_n_poloidal_2(z0,r0,phigrill(i),
     &            anzin(n),cnteta,cnphi)
                 endif 
               endif                  
           write(*,*)'grill_lh after nphiteta iray,z0,r0,phigrill(i)',
     &             iray,z0,r0,phigrill(i)
                 write(*,*)'grill_lh after nphiteta cnteta,cnphi',
     &          cnteta,cnphi

	        arnphi(iray)=cnphi
	        arntheta(iray)=cnteta
c-----------------------------------------------------
	        powinilh(iray)=fpwth(j)*pwcpl(n)
	write(*,*)'iray,arzu0(iray),arru0(iray),arphiu0(iray)',
     1iray,arzu0(iray),arru0(iray),arphiu0(iray)
        write(*,*)'arnphi(iray),arntheta(iray),powinilh(iray)',
     1arnphi(iray),arntheta(iray),powinilh(iray)

c calculation of the initial values of the n_parallel width
	        if (nnkpari.eq.1) then
                   write(*,*)'grill_lh nnkpari.eq.1'
        	   wdnpar0(iray)=dabs(anzin(n))*0.05d0
c        	   wdnpar0(iray)=dabs(anzin(n))*0.1d0
	        else
        	  wdnpar0(iray)=dabs(difrnpar)/dfloat(nnkpari)
	        endif
                write(*,*)'grill_lh  i_n_poloidal,iray,wdnpar0(iray)',
     +                               i_n_poloidal,iray,wdnpar0(iray)

	      end do !n
            endif !i_n_poloidal=1,2,3

            if(i_n_poloidal.eq.4)then 
c-------------i_n_poloidal=4
	      do ntor=1,nnktori
                do npol=1,nnkpoli
c	           write(*,*)'1 ntor,npol',ntor,npol
	           iray=iray+1
	           if (iray.gt.nraymaxl) then
	         write(*,*)'in grill_lh iray=',iray,'nraymaxl',nraymaxl
	         write(*,*)'in grill_lh the total number of the rays is'
		 write(*,*)'greater than number of elements in arrays'
		 write(*,*)'iray.gt.nraymaxl'
		 write(*,*)'please change parameter nraymaxl in param.i'
		      stop
	           endif

	           arzu0(iray)=z0
	           arru0(iray)=r0
	           arphiu0(iray)=phigrill(i)

c-----------------------------------------------------
	           arnphi(iray)=anztorin(ntor)
	           arntheta(iray)=anzpolin(npol)
c-----------------------------------------------------
	           powinilh(iray)=fpwth(j)*pwcpl_tp(ntor,npol)
	         write(*,*)'iray,arzu0(iray),arru0(iray),arphiu0(iray)',
     1             iray,arzu0(iray),arru0(iray),arphiu0(iray)
                 write(*,*)'arnphi(iray),arntheta(iray),powinilh(iray)',
     1             arnphi(iray),arntheta(iray),powinilh(iray)

c Setting the initial values of the n_parallel width
                 wdnpar0(iray)=hnwidth

                 write(*,*)'grill_lh i_n_poloidal,iray,wdnpar0(iray)',
     +                               i_n_poloidal,iray,wdnpar0(iray)

	        end do !npol
               enddo !ntor
            endif !i_n_poloidal=4
              
	 end do !j
      end do !i=1,ngrill
      nray=iray		   

cSAP100503 for test one ray from all rays
c      powinilh(2)=0.d0
c      powinilh(1)=0.d0
c      powinilh(3)=0.d0
c      powinilh(4)=0.d0
c      powinilh(5)=0.d0


      powert=0.d0
      do iray=1,nray
         powert=powert+powinilh(iray)
         write(*,*)'iray,arzu0(iray),arru0(iray),arphiu0(iray)',
     &              iray,arzu0(iray),arru0(iray),arphiu0(iray)
      enddo
      write(*,*)'grill_lh.f nray,powert',nray,powert
      !pause
      
c--------------------------------------------------------------------

      return
      end subroutine grill_lh



c for test  only
      double precision
     1function psiff(rho,psimag)
      implicit double precision (a-h,o-z)
      psiff=psimag*(1.d0-rho**2)
      return
      end

c        **********************nphiteta************************
c        * This subroutine calculates		              *
c        * the components cnpar_phi and cnpar_theta.          *
c        * They are the components of the parallel
c        * refructive index  N_parallel                       *
c        * in the point(z,r,phi).                             *
c        * It uses the component cnpar , which is the	      *	
c        * parallel to the magnetic field B     	      *
c        * refructive index ( N_parallel=cnpar) 	      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c        z,r,phi,cnpar                                             !
c        output parameters					   !
c        cnpar_theta,cnpar_phi						   !
c------------------------------------------------------------------
c        it uses the following function:                           !
c        b(z,r,phi)                                                !
c------------------------------------------------------------------
cSAP
c      subroutine nphiteta(z,r,phi,cnpar,cnteta,cnphi)
      subroutine nphiteta(z,r,phi,cnpar,cnpar_theta,cnpar_phi)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
c---------------------------------------------------------
c     e_theta=\hat{theta}=grad(psi) x \hat{phi}/abs(grad(psi)).
c     Here e_theta=\hat{theta} is the poloidal unit vector.
c	           z       r      phi
c       e_phi= {   0   ,   0    ,   1    }
c       e_psi= { dpsidz, dpsidr ,   0    }/mod(grad(psi))
c	e_theta={ dpsidr,-dpsidz ,   0    }/mod(grad(psi))
c
c              	   x	                y                z
c       e_phi ={-sin(phiu0)      ,  cos(phiu0)     ,     0    }
c	e_theta={-cos(phiu0)dpsidz,-sin(phiu0)dpsidz,   dpsidr }/
c               abs(grad(psi))
c	e_r   ={cos(phiu0)       , sin(phiu0)      ,     0    }
c
c       N_theta=(vector{N_par}*e_teta) ,N_phi=(vector{N_par}*e_phi)
c       vector{N_par}=N_par*{e_z*B_z+e_r*B_r+e_phi*B_phi}/abs(B)
c                                             
c       vector{N_par}=N_par*{e_x*(B_r*cos(phi)-B_phi*sin(phi))+
c                            e_y*(B_r*sin(phi)+B_phi*cos(phi))+
c                            e_z*B_z}/abs(B)
c       Npar_theta=N_par*(-cos(phi)*dpsidz*(B_r*cos(phi)-B_phi*sin(phi))
c                     -sin(phi)*dpsidz*(B_r*sin(phi)+B_phi*cos(phi))
c                     +dpsidr*B_z)/(abs(B)*abs*grad(psi))
c       Npar_phi=N_par*(-sin(phi)*(B_r*cos(phi)-B_phi*sin(phi))
c                    +cos(phi)*(B_r*sin(phi)+B_phi*cos(phi)))/abs(B)
c----------------------------------------------------------------------
      bmod=b(z,r,phi)
      write(*,*)'in grill_lh:z r,phi,bmod',z,r,phi,bmod
      cnphi=cnpar*bphi/bmod
      bpol1=dsqrt(bz*bz+br*br)
      gradpsi=dsqrt(dpdzd*dpdzd+dpdrd*dpdrd)
      if(gradpsi.gt.0.d0)then !YuP[2020-08-17] Added check of gradpsi=0
        bpol=(bz*dpdrd-br*dpdzd)/gradpsi
      else
        bpol=0.d0
      endif

      write(*,*)'in nphiteta bpol,bpol1',bpol,bpol1
      write(*,*)'in nphiteta bz,br,bphi,bmod',bz,br,bphi,bmod
      write(*,*)'in nphiteta dpdrd,r,bz1',dpdrd,r,-dpdrd/r
      write(*,*)'in nphiteta dpdzd,r,br1',dpdzd,r,dpdzd/r
      cntheta=cnpar*bpol/bmod
      if(gradpsi.gt.0.d0)then !YuP[2020-08-17] Added check of gradpsi=0
        cnpar1=(cnphi*bphi+cntheta*(bz*dpdrd-br*dpdzd)/gradpsi)/bmod
      else
        cnpar1=cnphi*bphi/bmod
      endif
      write(*,*)'in nphiteta cnpar,cnpar1',cnpar,cnpar1
c--------------------------------------------------------------------
      write(*,*)'nphiteta cnpar1',cnpar1,'cntheta',cntheta,'cnphi',cnphi
c--------------------------------------------------------------------
      cs=dcos(phi)
      sn=dsin(phi)
cSAP090513
c     cnteta=cnpar*(-cs*dpdzd*(br*cs-bphi*sn)
      if(gradpsi.gt.0.d0)then !YuP[2020-08-17] Added check of gradpsi=0
        cnpar_theta=cnpar*(-cs*dpdzd*(br*cs-bphi*sn)
     1              -sn*dpdzd*(br*sn+bphi*cs)
     1              +dpdrd*bz)/(bmod*gradpsi)
      else
        cnpar_theta=0.d0
      endif
cSAP090513
c     cnphi=cnpar*(-sn*(br*cs-bphi*sn)+cs*(br*sn+bphi*cs))/bmod
      cnpar_phi=cnpar*(-sn*(br*cs-bphi*sn)+cs*(br*sn+bphi*cs))/bmod
c      write(*,*)'nphiteta 1 cntheta',cntheta,'cnphi',cnphi
      write(*,*)'nphiteta 1 cnpar_theta',cnpar_theta,
     &'cnpar_phi',cnpar_phi
c      stop
c--------------------------------------------------------------------
      return
      end

      subroutine rho_ini_LHFW(theta,phi,cnpar,
     &i_n_poloidal,n_theta_pol,n_toroidal,
     &rho_ini,z_ini,r_ini,cntheta_ini,cnphi_ini,
     & i_rho_ini_LHFW_found)
c-----finds the small radius rho_ini at the vector rho^ where
c     the wave specified by ioxm (LH or  FW) has the cutoff.
c     The vector rho^ starting at the magnetic axis O(xma,yma).
c     The poloidal angle of this vector is theta (radians).
c
c     If (this subroutine could find rho_ini) then
c        it will set: i_rho_ini_LHFW_found=1
c     else
c        it will set i_rho_ini_LHWF_found=0

      implicit none
      !implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
cSAP091027
      include 'three.i'          
c-----input
      real*8 theta, ! poloidal angle(radians) measured from the
                    ! midplane in clockwise direction 
     &cnpar,        !N_parallel
     & phi          !toroidal angle (radians)

      integer i_n_poloidal ! the grill case (see inigrill)
 
      real*8  n_theta_pol, !N  poloidal for i_n_poloidal=4 case
     &        n_toroidal   !N  toroidal  for i_n_poloidal=4 case
c-----output
      real*8  rho_ini, !normalized small radius in cutoff point
     & z_ini,r_ini,
     & cntheta_ini,   !N poloidal in rho_ini point
     & cnphi_ini       !N_poloidal in rho_ini point 

      integer i_rho_ini_LHFW_found

c-----external zr_psith

      real*8 psi_rho,b
c     psi_rho(rho),zr_psith(psi,theta,z,r))
c
c-----local

cSm061107
cSAP091027
      real*8 psi,rho_loc,z,r,
     & hstep,cnz,cnr,cm,cnphi,cnteta,cnparz,
     & b_teta,cnper,gradpsi,p, 
     & z_b,r_b,    !coordinates at LCFS at (psi=psilim,theta)
     & rho_b_geom, ! distance from the magnetic axis
                   ! until (z_b,r_b) point 
     & rho_loc_geom ! geometrical sdistace from the magnetic axis
                    ! radius at (rho_loc,theta)
     

      integer 
     &n_gam,       !=500 set the number of angle points at the interval
     &n_nperp,     !=1000 the number of nperp points at the plot D(N_perp) 
     &iraystop_loc

!      write(*,*)'grill_lh.f rho_ini_LHFW i_n_poloidal,n_theta_pol',
!     &i_n_poloidal,n_theta_pol

      hstep=rho_step_find_LHFW_cutoff 

      pi=4.d0*datan(1.d0)
       
cSm061107
c     rho_loc =1.d0-rho_step_find_LHFW_cutoff  !initialization
      rho_loc=rho_initial_find_LHFW_cutoff     !initialization

 10   continue

cSm061107
      psi=psi_rho(rho_loc)
      !write(*,*)'grill_lh in rho_ini_LHFW rho_loc= ',rho_loc
      !write(*,*)'grill_lh in rho_ini_LHFW psi= ',psi

cSAP0921027
      if (rho_loc.lt.1.d0) then
         call zr_psith(psi,theta,z,r)
      else
c--------calculate coordinates (r_b,z_b) at the LCFS(psi=psilim,theta)          
         call zr_psith(psilim,theta,z_b,r_b)
         rho_b_geom=dsqrt((z_b-yma)**2+(r_b-xma)**2)    ! distance from magnetic
                                                     ! axis
                                                     !until (z_b,r_b) point 
         rho_loc_geom=rho_b_geom*rho_loc
         z=yma+rho_loc_geom*dsin(theta)
         r=xma+rho_loc_geom*dcos(theta)

      endif

      !write(*,*)'grill_lh psi,theta,z,r',psi,theta,z,r
   
ctest cutoff3 frequency     
c      bmod=b(z,r,phi)
c      write(*,*)'bmod',bmod
c      rho=rho_loc
c      x_e=x(z,r,phi,1)
c      y_e=y(z,r,phi,1)
c      cutoff3_d_omega=dsqrt(x_e+0.25*y_e**2)-0.5*y_e
c      write(*,*)'x_e,y_e',x_e,y_e
c      write(*,*)'grill_lh ps cutoff3_d_omega',cutoff3_d_omega
cendtest cutoff3 frequency    

cSm050408   
      if(i_n_poloidal.eq.1) then
        call nphiteta(z,r,phi,cnpar,cnteta,cnphi)
cSAP0921027
c       call cninit12(z,r,phi,cnpar,cnteta,cnphi,cnz,cnr,cm,iraystop
        call cninit12(z,r,phi,cnpar,cnteta,cnphi,cnz,cnr,cm,
     &               iraystop_loc)
      else    
!        write(*,*)'grill_lh in rho_ini_LHFW n_theta_pol,n_toroidal',
!     & n_theta_pol,n_toroidal

cSm050826
        if(i_n_poloidal.eq.2) then
           call grill_i_n_poloidal_2(z,r,phi,
     &     cnpar,n_theta_pol,n_toroidal)
c           write(*,*)'grill_lh n_theta_pol,n_toroidal',
c     &     n_theta_pol,n_toroidal

        endif  !i_n_poloidal.eq.2

!        write(*,*)'ion rho_ini_LHFW before cninit',
!     &  'cnparz,r,phi,cnpar,n_theta_pol,n_toroidal',
!     &  cnparz,r,phi,cnpar,n_theta_pol,n_toroidal
 
cSAP091027      
        call cninit(z,r,phi,cnpar,n_theta_pol,n_toroidal,
c     &              cnz,cnr,cm,iraystop)
     &              cnz,cnr,cm,iraystop_loc)
c        write(*,*)'ion rho_ini_LHFW after cninit iraystop',iraystop

        bmod=b(z,r,phi) 
cSm070804 cnteta=(cnz*bz+cnr*br)/bmod
        cnteta=-(cnz*bz+cnr*br)/dsqrt(bz**2+br**2)        
        cnphi=(cm/r)

!        write(*,*)'in rho_ini_LHFW after cninit cnteta,cnphi,
!     & iraystop_loc',
!     &  cnteta,cnphi,iraystop_loc

      endif ! i_n_poloidal.eq.1

c---------------------------------------------------------------
c     create the plot:
c     cold plasma dispersion function on N_perp
c     at the initial ray point
c---------------------------------------------------------------
      if (i_plot_disp_cold.eq.1) then
        n_nperp=1000 ! the number of nperp points at the plot D(N_perp)
cSAP091027
c        if (iraystop.eq.0) then
        if (iraystop_loc.eq.0) then
           cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar**2)   
        else
           cnper =-1.d0
        endif

        call plot_disp_cold_grill(z,r,phi,cnpar,0.d0,n_nperp)
        n_gam=500        ! set the number of angle points at the interval
                         ! (0,2*pi) at the wave normal surface plot

      endif
      !write(*,*)'grill_lh.f in  rho_ini_LHFW after '
      !write(*,*)'plot_disp_cold_grill and wave_normal_surface'
c---------------------------------------------------------------
cSAP091027
c      if(iraystop.eq.1) then
      if(iraystop_loc.eq.1) then
c-------the given wave does not exist in this point
        !write(*,*)'in rho_ini_LHFW rho_loc,hstep',rho_loc,hstep

cSm061107
        rho_loc=rho_loc-hstep
        !write(*,*)'in rho_ini_LHFW new rho_loc',rho_loc
        if (rho_loc.lt.0d0) then 
           !write(*,*)'rho_ini_LFW did not find the rho point'
cSAP081028
           i_rho_ini_LHFW_found=0
           return
c           stop
        endif

        go to 10                   
      else
        cntheta_ini=cnteta
        cnphi_ini=cnphi
        cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar**2)
        !write(*,*)'grill_lh.f rho_ini_LHFW cnper',cnper
        if (i_n_poloidal.eq.2) then 
c---------check the condition
          gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
          if(gradpsi.gt.0.d0)then !YuP[2020-08-17] Added check of gradpsi=0
            b_teta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field
          else
            b_teta=0.d0
          endif
          p=cnper**2-((n_theta_pol-cnpar*b_teta/bmod)*(bmod/bphi))**2
          !write(*,*)'p ',p
          if (p.lt.0.d0) then
cSm061107
            rho_loc=rho_loc-hstep
            goto 10
          endif
        endif !i_n_poloidal.eq.2

cSAP081028
        i_rho_ini_LHFW_found=1        
!        write(*,*)'in sub rho_ini_LHFW i_rho_ini_LHFW_found=',
!     &            i_rho_ini_LHFW_found

cSm061107
        rho_ini=rho_loc

        r_ini=r
        z_ini=z
        cntheta_ini=cnteta
        cnphi_ini=cnphi
        cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar**2)
        !write(*,*)'in rho_ini_LHFW rho_ini,cnper',rho_ini,cnper

      endif !iraystop.eq.1

      return
      end


      subroutine  integral_1(func,a,b,n,x,ar_integral,total_integral)
c-----calculates integral{a,x}f(x)dx and puts it in ar_integral
c     a<  x <b
c     using the formula
c     ar_integral(j)=Sum{k=1,j-1}[f(xc(k))*(x(k+1)-x(k))], j=2,n+1
c
c
c     input fync(x) !double precision function
c                   !It should be described as external in the program
c                   !using call intergal_1  
c     a < b         !the boundaries 
c
c     number of x mesh points
c-----output:
c     x(j)=a+(j-1)*(b-a)/n    Here j=1,n+1, x(1)=a, x(n+1)=b 
c     xc(j)=0.5(x(j)+x(j+1))  Here j=1,n    xc(1)=a+0.5(x(2)-x(1))
c                                            xc(n)=b-0.5(x(n+1)-x(n))
c     ar_integral(1)=0.d0
c     ar_integral(j)=Sum{k=1,j-1}[f(xc(k))*(x(k+1)-x(k))], j=2,n+1
c
c     total_integral=ar_integral(n+1)

      implicit none

c-----input      
      real*8 a,b 
      integer n 
      double precision func

c-----output
      double precision x(n+1),xc(n),ar_integral(n+1),total_integral
  
c-----locals
      integer j,k
      double precision step 

      step=(b-a)/n

      do j=1,n+1         
        x(j)=a+(j-1)*step
      enddo
c      write(*,*)'integral_1 x',x

cSAP080726
c      do j=1,n+1      
      do j=1,n
        xc(j)=x(j)+0.5d0*step
      enddo
c      write(*,*)'integral_1 xc',xc


      do j=1,n+1
        ar_integral(j)=0.d0
      enddo

      do j=2,n+1
        ar_integral(j)=ar_integral(j-1)+step*func(xc(j-1))
c        write(*,*)'j,xc(j-1),func(xc(j-1)),ar_integral(j)',
c     &             j,xc(j-1),func(xc(j-1)),ar_integral(j)
      enddo
      total_integral=ar_integral(n+1)

c      do j=1,n
c        write(*,*)'j,xc(j),x(j+1),ar_integral(j+1)',
c     *             j,xc(j),x(j+1),ar_integral(j+1)
c      enddo
c      write(*,*)'total_integral',total_integral

      return
      end



      subroutine create_equal_mesh(a,b,func,n_mesh,
     &x_mesh,x_mesh_c,f_equal)

c-----creates X mesh and array for the function func: 
c     array x_mesh_c(n_mesh)  (a<  x_mesh_c <b)
c     array x_mesh(n_mesh+1) (a=<  x_mesh =<b)
c     array f_equal(n_mesh) of the function func at the created mesh
c
c     This mesh will be used to calculate the integral
c     integral{a,b}f(x)dx=sum{j=1,n_mesh}f_equal(j)*del(j)
c
c     x_mesh(1)=a
c     x_mesh(n_mesh+1)=b
c     a<x_mesh(j)<b   j=2,...,n_mesh
c
c     del(j)=x_mesh(j+1)-x_mesh(j) j=1,...,n_mesh
c     
c     x_mesh_c(j)=x_mesh(j)+0.5*del(j) j=1,...,n_mesh
c
c     The mesh and the function array will be calculate to give
c
c     f_equal(j)*del(j)=Constant=total_integral/n_mesh
c
c
c     total_integral is calculated from the given function func 
c     at the equispaced mesh with n points
c     x(j)=a+(j-1)*(b-a)/n    Here j=1,n+1, x(1)=a, x(n+1)=b 
c     xc(j)=0.5(x(j)+x(j+1))  Here j=1,n    xc(1)=a+0.5(x(2)-x(1))
c                                            xc(n)=b-0.5(x(n+1)-x(n))
c     ar_integral(1)=0.d0
c     ar_integral(j)=Sum{k=1,j-1}[f(xc(k))*(x(k+1)-x(k))], j=2,n+1
c
c     total_integral=ar_integral(n+1)
c
c-----------------------------------------------------------------------

      implicit none
     
      integer n !the number of points for the equispaced mesh x(n+1)
      parameter (n=1000)

c-----input
      integer n_mesh ! the number of new mesh points
      real*8 a,b     ! the boundaries of the integration area
      
c-----externals
      real*8 func
      external func

c-----output
      real*8 x_mesh(n_mesh+1),x_mesh_c(n_mesh),f_equal(n_mesh)

c-----locals
      integer j,k,k0
      real*8 integral_value
      real*8 x(n+1),xc(n),ar_integral(n+1),total_integral

c------------------------------------------------------------
c     calculate the integral from the function func
c     ar_integral(j) j=1,...,n+1
c     at the equispaced mesh x()

c      write(*,*)'grill_lh.f: create_equal_mesh a,b,n',a,b,n

c------------------------------------------------------------
      call integral_1(func,a,b,n,x,ar_integral,total_integral)
c------------------------------------------------------------

c      write(*,*)'create_equal_mesh: ar_integral',ar_integral
c      write(*,*)'create_equal_mesh: x', x

c------------------------------------------------------------
c     x_mesh and the function array will be calculated to give
c
c     f_equal(j)*del(j)=Constant=total_integral/n_mesh
c---------------------------------------------------------
      k0=1
      x_mesh(1)=a
      x_mesh(n_mesh+1)=b
      
      do j=2,n_mesh

        integral_value=(j-1)*total_integral/n_mesh
c        write(*,*)'j,integral_value',j,integral_value
        do k=k0,n+1

c           write(*,*)'k0,k,integral_value,ar_integral(k)',
c     &                k0,k,integral_value,ar_integral(k)

          if(integral_value.le.ar_integral(k)) then
c---------------------------------------------------------------------
c           The equation to find xmesh(j)
c           (integral_value-ar_integral(k-1))/(x_mesh(j)-x(k-1))=
c           =(ar_integral(k)-ar_integral(k-1))/(x_(k)-x(k-1))
c-----------------------------------------------------------------------

c             write(*,*)'k,x(k-1),x(k)',k,x(k-1),x(k)

            x_mesh(j)=x(k-1)+(integral_value-ar_integral(k-1))*
     &      (x(k)-x(k-1))/(ar_integral(k)-ar_integral(k-1))

c            write(*,*)'x_mesh(j)',x_mesh(j)

            f_equal(j-1)=(total_integral/n_mesh)/(x_mesh(j)-x_mesh(j-1))
            x_mesh_c(j-1)=0.5d0*(x_mesh(j)+x_mesh(j-1))
           
            if(k.eq.n+1) then
              goto 10
            else
              k0=k+1
              goto 20
            endif
          endif
        enddo !k
 20     continue     
      enddo !j
     
      x_mesh_c(n_mesh)=0.5d0*(x_mesh(n_mesh+1)+x_mesh(n_mesh))
      f_equal(n_mesh)=(total_integral/n_mesh)/
     &                (x_mesh(n_mesh+1)-x_mesh(n_mesh))
           
 10   continue

c      write(*,*)'x_mesh',x_mesh
c      write(*,*)'x_mesh_c',x_mesh_c
c      write(*,*)'f_equal',f_equal

      do j=1,n_mesh
         f_equal(j)= f_equal(j)*(x_mesh(j+1)-x_mesh(j))
c         write(*,*)'j,x_mesh_c(j), f_equal(j) ',j,x_mesh_c(j),f_equal(j) 
      enddo

      return
      end


       double precision function f_pow_poloid_0(x)
c-----------------------------------------------------
c     poloidal distribution of power spectrum
c------------------------------------------------------

      implicit none
c-----input:
      real*8 x ! normalized deviation of the poloidal angle
               ! near the central polidal angle of the grill
               ! x=pi*(del_theta/aperture)
               ! -pi/2< x < pi/2
               ! Here 
               ! del_theta=theta_pol-theta_pol_center
               ! aperture (rad) poloidal aperture of the grill
               ! -0.5*appertura < del_theta < 0.5*appertura

      f_pow_poloid_0=(dcos(x))**2/ 0.5d0 
    
      return
      end


      double precision function f_pow_poloid_1(del_theta)
c-----------------------------------------------------
c     poloidal distribution of power spectrum
c------------------------------------------------------
      implicit none
c-----input:
      real*8 aperture ! (rad) poloidal aperture of the grill
                       ! -0.5*appertura < del_theta < 0.5*appertura
      common /aperture/ aperture

      real*8 del_theta ! (rad) deviation of the poloidal angle
                       ! near the central polidal angle of the grill
                       ! del_theta=theta_pol-theta_pol_centr     
c-----local
      real*8 pi

      pi=4.d0*datan(1.d0)

c      write(*,*)'grill_lh.h in  f_pow_poloid_1:del_theta,aperture',
c     &del_theta,aperture

      f_pow_poloid_1=(dcos(pi*del_theta/aperture))**2
      f_pow_poloid_1=f_pow_poloid_1/0.5d0

      return
      end

      double precision function f_pow_npar_2(x)
c-----------------------------------------------------
c     n_parallel distribution of power spectrum: f_pow_npar_2=(sin(x)/x)**2  
c------------------------------------------------------
      implicit none
c-----input:
      real*8 x
c-----------------------------------------------------
c     it as assumed in grill_lh: dabs(x).le.pi
c-----------------------------------------------------
c      write(*,*)'grill_lh.h in   f_pow_npar_2: x,dabs(x)',x,dabs(x)

      
      if (dabs(x).le.1.d-100) then 
         f_pow_npar_2=1.d0
      else
         f_pow_npar_2=(dsin(x)/x)**2  
      endif

      return
      end

      
      double precision function f_pow_npar_3(x)
c-----------------------------------------------------
c     n_parallel distribution of power spectrum
c     f_pow_npar_3=exp(-x*x)
c------------------------------------------------------
      implicit none
c-----input:
      real*8 x
c-----------------------------------------------------
c     it as assumed in grill_lh: dabs(x).le.1.d0
c-----------------------------------------------------
c     write(*,*)'grill_lh.h in   f_pow_npar_3:xd',x
      
      f_pow_npar_3=dexp(-x*x)  
   
      return
      end


      subroutine grill_i_n_poloidal_2(z,r,phi,cnpar,cntheta,cnphi)
c-----------------------------------------------------      
c     calculates the toroidal refractive index (cnphi) from the given
c     parallel (cnpar) and poloidal (cntheta) refractive indexes             
c     in the given point P(z,r,phi).
c
c     It is used for i_n_poloidal=2 grill case
c
c     It works for cnpar.ne.0
c----------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
c-----input 
      real*8 z,r,phi, ! the coordinates of the given point P
     &cnpar,cntheta   ! parallel and poidal refractive indexes
      
c-----output
      real*8 cnphi   !toroidal refractive index  

c-----externals
      real*8 b

c-----locals
      real*8 
     & b_theta, ! poloidal magnetic field
     & gradpsi,cnpar_theta,cnpar_phi,cnper_theta,cnper_phi


      bmod=b(z,r,phi) ! put the magnetic filed componets into one.i

      gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
      if(gradpsi.gt.0.d0)then !YuP[2020-08-17] Added check of gradpsi=0
        b_theta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field
      else
        b_theta=0.d0
      endif

c-----the parallel refractive index components
      cnpar_phi = cnpar*bphi/bmod  
      cnpar_theta = cnpar* b_theta/bmod

      write(*,*)'grill_i_n_poloidal_2: cnpar_phi,cnpar_theta',
     &cnpar_phi,cnpar_theta

c-----perpendicular refractive index components 
      cnper_theta = cntheta - cnpar_theta
      write(*,*)'grill_i_n_poloidal_2: cntheta,cnper_theta',
     &cntheta,cnper_theta
cSm060830
cSm060830
c      cnper_phi=-cnper_theta*cnpar_theta/cnpar_phi
c      write(*,*)'grill_lh old formula cnper_phi',cnper_phi
      cnper_phi=-cnper_theta*b_theta/bphi

      write(*,*)'cnper_phi,cnper_theta,b_theta,bphi,bmod',
     &cnper_phi,cnper_theta,b_theta,bphi,bmod

      write(*,*)'grill_lh new formula cnper_phi',cnper_phi
c-----toroidal component od the total refractive index
      cnphi=cnpar_phi+cnper_phi

      return
      end


       
