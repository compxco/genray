


c************************************************************
      character version*64
      parameter(version="genray_v11.0_210212")
      ! previous: genray_v10.15_201206 genray_v10.14_200831
      ! genray_v10.13_200117 genray_v10.12_180912 
      ! genray_v10.12_180529 
c************************************************************

c************************************************************
c number of rays in common/cone/ EC-cone
       integer nraymax
       parameter(nraymax=240)
c nraymax must be greater than or equal to nray 
c (nray is calculated)

c maximum number of EC cones
       integer nconea
       parameter(nconea=20)
c maximal number of radial bins at the launching disk
       integer n_mesh_disk_radial_bin_a
       parameter (n_mesh_disk_radial_bin_a=5)
c************************************************************
c     for common/five/
      integer nxeqda,nyeqda,nlimit,nx4a,ny4a,nrya,nlim4

c      parameter (nxeqda=201,nyeqda=201,nlimit=101)
c      parameter (nxeqda=129,nyeqda=257,nlimit=101)
c      parameter (nxeqda=136,nyeqda=257,nlimit=101)
      parameter (nxeqda=257,nyeqda=257,nlimit=101)
      !parameter (nxeqda=257,nyeqda=257,nlimit=501) !YuP[2020-03] Increased nlimit
      !However, it slightly changes the results,
      !probably because of accuracy for tracing surfaces in subr.gr2new.
      !Need to keep it in mind when performing auto-tests 
      !(verification against data in genray.nc)
      
      parameter (nx4a=nxeqda+4,ny4a=nyeqda+4,nrya=ny4a)
      parameter (nlim4=nlimit+4)
c It should be nrya=max(nxeqda,nyeqda)+4
c************************************************************
c     for common/fourb/
      integer nves
      parameter (nves=202) !YuP[2020-03] Increased nves; was (nves=62)
      !parameter (nves=502) !YuP[2020-03] Increased nves; was (nves=62)
c************************************************************
c     for common gr.i
      INTEGER NL,NP,nteta,npsi,nteta1
      real*8 epspsi
      parameter (nteta=nlimit-1,npsi=50,epspsi=0.0001d0,nteta1=nteta+1)
      !YuP[2020-03] set nteta=nlimit-1, as required by code
cSAP091211
c      parameter (nteta=10000,npsi=50,epspsi=0.0001d0,nteta1=nteta+1)
cSAP080819
c      parameter (nteta=1010,npsi=1000,epspsi=0.0001d0,nteta1=nteta+1)
      PARAMETER (NL=5,NP=nteta)
c npsi is a number of contours
c nteta  is a number of the points along the each contours
c                       (in the poloidal direction)
c epspsi is the accuracy for the determination of the contours
c         poins coordinates	(zpsi,rpsi)
c NL is a number of contours for the ploting of the tokamak
c              cross section
c NP is a number of points along each contours in poloidal direction
c     for the plotting of the tokamak cross section
c************************************************************
c     for common/grill/
      integer ngrilla,nnkprmax,nraymaxl,nthinmax
c      parameter (ngrilla=10,nnkprmax=60,nraymaxl=1000,nthinmax=64)
      parameter (ngrilla=10,nnkprmax=60,nraymaxl=1000,nthinmax=160)
c      parameter (ngrilla=10,nnkprmax=60,nraymaxl=1000,nthinmax=40)
c     for N_toroidal N_poloidal initial condition i_n_poloidal=4 
      integer nnktormax,nnkpolmax
      parameter (nnktormax=60,nnkpolmax=60)

c------------------------------------
c nnkprmax =max{i=1,ngril1}nnkpar(i)
c nnktormax=max{i=1,ngril1}nnktor(i)
c nnkpolmax=max{i=1,ngril1}nnkpol(i)
c------------------------------------
c ngrilla    is a maximal number of the N_parallel spectra
c ngrill     is a number of the N_parallel spectra
c nnkprmax=max{i=1,ngrill}nnkpar(i)
c nthinmax=max{i=1,ngrill}nthin(i)
c nraymaxl is the max number of rays from the grill
c************************************************************
c     for common/ions/
      integer nbulka,nbulkma
c      parameter (nbulka=12)
c      parameter (nbulka=3)
c      parameter (nbulka=4)
c      parameter (nbulka=5)
cBH111106: For SWIM Plasma State, need nbulka=.ge.8
      parameter (nbulka=8)

      parameter (nbulkma=nbulka-1)
c ncomp=nbulk -the number of plasma species.ge.1 & .le.nbulka
c************************************************************
c     for common/loopb/
c      parameter (ntor=21,nloop=24,npol=50)
c ntor  is the number of points in the 'toroidal angle' mesh
c nloop is the number of the loops
c npol  is the number of points along the loop
c************************************************************
c     for common/onetwo/
      integer NRA
c      PARAMETER (NRA=21)
c      PARAMETER (NRA=51)
       PARAMETER (NRA=201)
c NRA is the maximal number of bin boundaries in the small radius direction 
c for the calculation of the power and current drive radial profiles.
c Power and current is tabulated at (NR-1) bin centers.
c************************************************************
c     for common/rho/
      integer npsi4,nteta14,nzy
      parameter (npsi4=npsi+4,nteta14=nteta1+4)
c npsi  and nteta1 are in common/gr/
c nzy=max(npsi,nteta1)+4
      parameter(nzy=nteta1+4)
c************************************************************
c     for common/six/ and /lsc_approach_nml/
      integer ndensa,ndens4a
      parameter (ndensa=201)
c      parameter (ndensa=51)
c      parameter (ndensa=101)
      parameter (ndens4a=ndensa+4)
c ndensa is the max number of points in arrays with the 
c plasma density, temperature, zeff. tpop, and vflow
c************************************************************
c     for common/write/
cSAP090809 delete nfreqa,nraya
cBH130508:  Restore nraya,nfreqa for limited purpose of MPI subs
cBH130508:  senddata/recvdata.  Compiler does not permit equivencing
cBH130508:  with variable dimension specification.
      integer nraya,nfreqa
      parameter (nraya=300)
      parameter (nfreqa=51)

      integer nrelta
c,nfreqa
      integer n_relt_harma,n_relt_harm1a,n_relt_harm2a
c      parameter (nrelta=5000) !BH080128
c      parameter (nrelta=1000)
c       parameter (nrelta=1200) !SAP 071203
c      parameter (nrelta=2000)

       parameter (nrelta=10000)
     
c nrelta is maximum value for nrelt
c nrelt is the max number of the ray elements along every ray

c nraya is the max value for nray (the number of the rays)


c nfreqa is the max value for nfreq 
c (nfreq is the number of the emission frequencies)
c n_relt_harma  is the number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations  n=<n_relt_harm 

c      parameter (n_relt_harma=3)
c      parameter (n_relt_harm1a=-3)
c      parameter (n_relt_harm2a=3)
      parameter (n_relt_harma=11)
      parameter (n_relt_harm1a=-5)
      parameter (n_relt_harm2a=5)
c************************************************************
c     for common/scatnper/ (the data for the n_perp scattering )
      integer nscat_n
      parameter (nscat_n=10)
c nscat is a namber of point in rhoscat(nscat_n) and scatd(0:nscat_n)
c*************************************************************
c     for common/dskincomm/ the data for reading dskin file
c     with 3D distribution function f(iya,jx,lrz,ngen) written by CQL3D code
c     iya is the number of the points in the pitch angle mesh
c     jx is the  number of the points in the momentum  mesh
c     lrz is the number of the points in the radial mesh
c     ngen is the number of the plasma species
cSAP090808 delete jxa,lrza,ngena,jxa1  
      integer iya
      parameter(iya=120)
      integer iya1          
      parameter(iya1=iya+1)
c***************************************************************
c     for common/emissa/
      integer nrelta4
      parameter (nrelta4=nrelta+4)
cSAP090808 delete jx_kin_a  
c**************************************************************
c     for distribution function spline approximation
c     in general,need nwka=1+3*max(jxa,iya)
c     for common/distrib/
cSAP090808 delete parameter nwka
c*************************************************************
c     for common /output/
      integer n_plot_dispa,n_plot_disp_colda,
     &m_r_nperp_a,m_i_nperp_a,n_contour_plot_disp_a
      parameter (n_plot_dispa=100)
      parameter (n_plot_disp_colda=100)
c      parameter(m_r_nperp_a=10,m_i_nperp_a=10,n_contour_plot_disp_a=10)
      parameter(m_r_nperp_a=20,m_i_nperp_a=20,n_contour_plot_disp_a=20)

c*************************************************************
c-----for hot plasma roots n_hot_roots_a max number of hot plasma root
      integer n_hot_roots_a
      parameter (n_hot_roots_a=4)

c******************************************************************
c     for small output step near EC resonance points for power calculation 
c     It will be used for data at namelist /numercl/ and common /one/
      integer n_power_switch_resonance_a
      parameter (n_power_switch_resonance_a=3)

c******************************************************************
c     for toray EC launch
      integer gzonemax
      parameter (gzonemax=20)
c******************************************************************
c     for adj.i
cSAP090828 delete npsi0_a,nthp0_a,nmax_chi_a,imax_chi_a  
c     integer npsi0_a,nthp0_a,nmax_chi_a,imax_chi_a           
c     npsi0_a  is a maximal number of radial points for chi function calculation
c     nthp0_a  is a maximal number of poloidal points for integration along field line
c     nmax_chi_a is a maximal number of grid points in u_0 direction (for adj mesh)
c     imax_chi_a is a maximal the number of grid points in pitch angle(for adj mesh) 

c********************************************************************
c     to read wall and limiter coordinates (r,z)
c     for writencdf.i:
      integer n_wall_a,n_limiter_a,max_limiters_a
c     n_wall_a,n_limiter_a are maximal numbes of wall and limiter points
c     max_limiters_a is a maximal number of limiters
      parameter (n_wall_a=200)
      parameter (n_limiter_a=200)
      parameter (max_limiters_a=1)

c     maximal number of wall points having the given poloidal angle theta
      integer n_rho_wall_a
      parameter (n_rho_wall_a=10)

c***********************************************************************
c     to read normalized exponential density falls at poloidal mesh
c     for edge_prof_nml,i
c     n_pol_edge_dens_a is a maximal number of mesh points
      integer  n_pol_edge_dens_a
      parameter (n_pol_edge_dens_a=100)
c********************************************************************
c     to create fine meshes of wall coordinates with additional points
      integer n_wall_add_a
      parameter (n_wall_add_a=10001) 
c*********************************************************************
c     for (R,Z) meshes rr_add, zz_add at the poloidal plane
c     for the density fall near the wall
c      nxeqd_add_a and nxeqd_add_a will be removed
c      integer
c     &nxeqd_add_a, !max number of r mesh points at the poloidal plane
c     &nyeqd_add_a  !max number of z mesh points at the poloidal plane

c      parameter (nxeqd_add_a=1000) 
c      parameter (nyeqd_add_a=1000)
c      parameter (nxeqd_add_a=10) 
c      parameter (nyeqd_add_a=10)
