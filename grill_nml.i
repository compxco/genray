
c     the variables for the grill parameters from namelists
c------------------------------------
c     nnkprmax=max{i=1,ngrill}nnkpar(i)
c------------------------------------
c     parameter (ngrilla) is set in param.i
      real*8
c------from namelist /grill/
     &  n_theta_pol,ksi_nperp,
     &  rhopsi0,thgrill,
     &  phigrill,height,
     &  anmin,anmax,
     &  powers,
     &  antormin,antormax,
     &  anpolmin,anpolmax,
     &  r0launch,z0launch,phi0launch,
     &  rlaunch,zlaunch,weight_power_launch
c----------------------------------------
      integer 
c------from namelist /grill/
     &  i_n_poloidal,i_rho_cutoff,ngrilld,ngrill,
     &  igrillpw,igrilltw,
     &  nthin,nnkpar,
     &  nnktor,nnkpol,
     &  ilaunch,i_grill_pol_mesh,
     &  i_grill_npar_ntor_npol_mesh

       common /grill_nmnl/
c------real*8
     &  n_theta_pol,ksi_nperp,
     &  rhopsi0(ngrilla),thgrill(ngrilla),
     &  phigrill(ngrilla),height(ngrilla),
     &  anmin(ngrilla),anmax(ngrilla),
     &  powers(ngrilla),
     &  antormin(ngrilla),antormax(ngrilla),
     &  anpolmin(ngrilla),anpolmax(ngrilla),
     &  rlaunch(nthinmax,ngrilla),
     &  zlaunch(nthinmax,ngrilla),
     &  weight_power_launch(nthinmax,ngrilla),
     &  r0launch,z0launch,phi0launch,

c-------integer 
     &  i_n_poloidal,i_rho_cutoff,ngrilld,ngrill,
     &  igrillpw,igrilltw,
     &  nthin(ngrilla),nnkpar(ngrilla),
     &  nnktor(ngrilla),nnkpol(ngrilla),
     &  ilaunch,i_grill_pol_mesh,
     &  i_grill_npar_ntor_npol_mesh
c--------------------------------------------------------------------      
c       ngrilla          is a max number of the poloidal grill angles	    
c       ngrill           is a number of the poloidal grill angles	    
c       thgrill(1:ngrilla) poloidal  angle of grill, measured counter
c                         clockwise from horizontal through the
c                         magnetic axis	(degrees).
c       height(1:ngrilla)  is a poloidal length (m) of grill 
c                         (giving poloidal power distribution of each grill).
c       nthin(1:ngrilla)   is a number of rays near the each poloidal center,
c                         simulating a grill   
c       phigrill(1;ngrilla)is a toroidal grill angle of grill (degrees)
c       anmin(1:ngrilla)   position of the left bound   
c                         of power spectrum P(n_parallel) (Can be neg).
c       anmax(1:ngrilla)   position of the right bounds  
c                         of power spectrum P(n_parallel)		    
c       nnkpar(1:ngrilla)  number of points  of power spectrum
c                         P(n_parallel)				    
c       powers(1:ngrilla)  power in one grill (MWatts)	    
c                         (total power of grill(in MWatts) will be   
c                         powtott=sum{powers}
c       rhopsi0(1:ngrilla) initial psi for wave front (0<rhopsi0<1)     
c       nnkprmax=max{i=1,ngrill}nnkpar(i)
c       nthinmax=max{i=1,ngrill}nthin(i)
c       wtheta(nthinmax)     poloidal angle mesh wtheta(j) radian      !
c                           near the grill poloidal angles thgrill(i)  !
c                           j=1,nthin(i)			       !
c       wtheta_edge(nthinmax+1) poloidal angle mesh [radian]           !
c                           near the grill poloidal angles thgrill(i)  !
c                           j=1,nthin(i)+1			       !       
c
c                           wtheta(j)=0.5*
c                           (wtheta_edge(j)+wtheta_edge(j+1)),j=1,nthin
c                            
c       anzin(nnkprmax)      are the points of mesh n_parallel	       
c                            P(n_parallel): anzin(j), j=1,nthin
c   
c       anzin_edge(nnkprmax+1) are the points of mesh n_parallel	       
c                            P(n_parallel):anzin_edge(j),j=1,=nnkpar(i)+1
c 
c                            anzin(j)=0.5(anzin_edge(j)+anzin_edge(j+1))
c                                     j=1,=nnkpar(i)
c
c       pwcpl(nnkprmax)      power spectrum pwcpl(n)    on	       !
c                             n_parallel, normalization		       !
c                             Sum{n=1,nnkpar(i)}pwspl(n)=powers(i)     !
c       wdnpar0(nraymaxl)    initial values of the n_parallel width    !
c                           nraymaxl mast be .ge. nray (the total      !
c                           number of rays)             	       !
c       igrillpw      specifies the form of N_parallel power spectra
c                    =1 power=powers/nnkpar, =2 power=sin**2x/x**2
c                    =3 exp
c       igrilltw      specifies the form poloidal variation of power,
c                    =1 uniform over height, =2 cos**2 variation.
c
c i_n_poloidal =1         The input parameter is N_parallel(from grill).
c by default =1           N_phi,N_theta are calculated from given N_parallel 
c                         N_rho=N_perpendicular(N_parallel) is determined 
c                         from the dispersion relation. It is directed
c                         along +,- gradient(psi) 
c
c i_n_poloidal =2         The input parameters: N_parallel(from grill)
c                         and  n_theta_pol. By default N_theta=0. 
c                         N_perpendicular(N_parallel) is determined 
c                         from the dispersion relation. 
c                         N_phi is calculated from N_parallel and N_theta
c                         N_rho is calculated form N_perpendicular, N_parallel
c                         and N_theta. 
c                         It is directed along +,- gradient(psi)
c
c i_n_poloidal=3          The given parameters: N_parallel and the angle
c                         0<<ksi_nperp<<180 between the vector N_perpendicular 
c                         and gradient(psi). By default ksi_nperp=0.
c                         N_perpendicular(N_parallel) is determined 
c                         from the dispersion relation.
c                         N_phi,N_theta and N_rho are calculated from
c                         N_parallel,N_perpendicular and ksi_nperp.
c
c n_theta_pol             The poloidal refractive index component
c                         It is used for i_n_poloidal =2     
c                         By_default n_theta=0.
c
c ksi_nperp               (degree) the angle 0<<ksi_nperp<<180
c                         between the vector N_perpendicular 
c                         and gradient(psi). By default ksi_nperp=0.
c-------------------------------------------------------------------
!  Calculations of the small radius value near the plasma edge
!  where LH or FW have cutoof.  
!  i_rho_cutoff=0 (default) no calculations
!              =1 use these calculations
!-----------------------------------------------
c       antormin(1:ngrilla)   position of the left bound   
c                of power spectrum P(n_toroidal) (Can be neg).
c       antormax(1:ngrilla)   position of the right bounds  
c                of power spectrum P(n_toroidal)		    
c       nnktor(1:ngrilla)  number of points  of power spectrum
c                P(n_toroidal in n_toroidal direction	
c       anpolmin(1:ngrilla)   position of the left bound   
c                of power spectrum P(n_poloidal) (Can be neg).
c       anpolmax(1:ngrilla)   position of the right bounds  
c                of power spectrum P(n_poloidal)		    
c       nnkpol(1:ngrilla)  number of points  of power spectrum
c                P(n_toroidal,n_poloidal) in n_poloidal direction	
c---------------------------------------------------------------------
c  ilaunch=1, to launch a single ray at r0launch,phi0launch,z0launch
c             in the plasma (meters and degs)
c         =0, no effect (the default)
c  This option is added for comparison with other codes.
c  r0launch is the major radius of the launch point [m]
c  z0launch is the vertical position of the launch pimt [m]
c  phi0launch is the toroidal angle of the launch point [degree] 
c---------------------------------------------------------------------

c       i_grill_pol_mesh: option specifying the poloidal mesh wtheta(j)
c                         near the central grill angle thgrill(i)
c                         =1 equaspaced mesh 
c                            wtheta(j)-wtheta(j-1)=zdth=Const(default)
c                         =2 poloidal mesh will be chosen to get the equal
c			     power fpwth(j) for all rays near the central 
c                            grill angle fpwth(j)=1/nthini
c  
c       i_grill_npar_ntor_npol_mesh: option specifying the refactive
c                         index meshes.
c
c                         For  i_n_poloidal=1,2,3 it is specifing
c                         n_parallel mesh anzin(n) for the power
c                         spectrum pwcpl(n) n=1,...,nnkpari
c                         =1 equaspaced mesh 
c                            anzin(n)-anzin(n-1)=hnpar=Const (default)
c                         =2 n_parallel mesh will be chosen to get the equal
c			     power pwcpl(n) for all rays in the given power
c                            spectrum  pwcpl(n)=1.d0/nnkpari 
c                            pwcpl(n)=power_spectrum(anzin(n))*
c                                     delta_npar_bin(n)= 1.d0/nnkpari
c
c                            For  i_n_poloidal=4 it is specifing two meshes:
c                            a) n_toroidal mesh anztorin(ntor) and             c                            b) n_poloidal mesh anzpolin(npolmesh) 
c                            for the power spectrum
c                            pwcpl_tp(1:nnktori,1:nnkpoli)=pwcpl_t*pwcpl_t 
c                         =1 equaspaced meshs (default)
c                            anztorin(ntor)- anztorin(ntor-1)=hntor=Const 
c                            anzpolin(npol)- anzpolin(npol-1)=hnpol=Const
c                         =2 the meshes anztorin(1:nntori) anzpolin(1:nnkpoli)
c                            will be chosen to get the equal
c			     power pwcpl_tp(ntor,npol) for all rays in 
c                            the given power spectrum 
c                            pwcpl_tp(ntor,npol)=1.d0/(nnktori*nnkpoli) 
c------------------------------------------------------------------------      c   pwcpl_t_ar(nnktormax) N toroidal power spectrum                
c   pwcpl_p_ar(nnkpolmax) N poloidal power spectrum 
c
c   anztorin_edge(nnktormax+1): N_toroidal mesh j=1:nnktori+1         
c   anzpolin_edge(nnkpolmax+1): N_poloidal mesh j=1:nnkpoli+1
c
c  rlaunch(nthinmax,ngrilla), launch points R coordinates [m]  
c  zlaunch(nthinmax,ngrilla), launch points Z coordinates [m]  
c  weight_power_launch(nthinmax,ngrilla) 
