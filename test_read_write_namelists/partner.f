      subroutine read_transport_prof(ndens,charge,dmass,
     &nbulk,dens1,temp1,tpop1,vflow1,zeff1,izeff,eqdskin,partner)
c--------------------------------------------------
c     This is an interface routine between the transport code [ONETWO]
c     and the all frequencies ray tracing code GENRAY.

c     It reads in plasma profile data from the file: genray_profs_in.txt
c     or genray_profs_in.nc.
c     [Below is subroutine write_transport_prof which writes
c      results of the genray run to genray_profs_out, for use
c      by the transport code.]
c     This data is read in after the code reads the genray
c      namelist data, and hence updates the plasma profiles.
c
c     If partner='genray_profs_in.nc':
c     the input plasma profiles are from a netCDF file genray_profs_in.nc
c     If partner='genray_profs_in.txt':
c     the input plasma profiles are from a text file genray_profs_in.txt
c     Both above cases also output genray calculated depostion and
c     current drive profiles in netCDF file 'genray_profs_out.nc'.
c     If partner='genray_profs_out.nc'
c     There is only 'genray_profs_out.nc' output of profiles, no input.
c     Presently, the profile data in only in bin-boundary format,
c     suitable for the ONETWO transport code.  Bin-centered profiles,
c     suitable for the Plasma State, out output into the general
c     genray netCDF output file.
c
c     nbulk from genray.dat is set to nspec, from genray_profs_in.
c     If nbulk from the genray.dat file is .eq.1, then it is reset
c     to 1 at the end of this subroutine so can do simple electron
c     calculation with genray of power depostion and CD.
c
c     Using this data, set the radial profiles on genray radial mesh:
c     density     dens1(ndensa,nbulka), filled in 1:ndens,1:nbulk
c     temperature temp1(ndensa,nbulka), filled in 1:ndens,1:nbulk
c     zeff        zeff1(ndensa), filled in 1:ndens
c     Also, set izeff=2, tpop1=1., vflow1=0.

c--------------------------------------------------     

c     call getioun(igenray,43)  ! local unit no. for files
                                ! genray_profs_in, genray_profs_out
      implicit none

      include 'param.i' ! gives nbulka,ndensa
      include 'transport_prof.i'
               !       in transport_prof.i':
               !nspeca is maximal number of transport code species,
               !      a parameter.
               !nja   is maximal number of points in the transport 
               !      code small radius mesh, a parameter.
               !r_transport(nja) is transport code small radius mesh [cm]

      integer iworkka
      parameter(iworkka=3*nja+1) ! is the length of the work array      

c-----input 
      integer 
     &ndens       !dimension of Genray small radial uniform mesh
c     &            !used for radial profiles
c     &ndensa,     !maximal value of ndens, a parameter
c     &            !set in param.i
c     &nbulka      !maximal number of plasma species (e+ions),   
c     &            !a parameter set in param.i
      real*8 charge(nbulka),dmass(nbulka) !species charge,weight numbers
c-----output
      integer nbulk  ! number of plasma species (e+ions)
      integer nbulk_save
      integer izeff  ! Set equal to 2 in this routine; gives method
                     ! of specification of the density profiles.
      double precision 
     &dens1(ndensa,nbulka),
     &temp1(ndensa,nbulka),
     &tpop1(ndensa,nbulka),
     &vflow1(ndensa,nbulka),
     &zeff1(ndensa)

      character (len=*) eqdskin,partner

c-----locals
      integer
c     & nj,          ! the number of radial points in transport profiles
c                    ! see include 'transport_prof.i' file
     & igenray,      ! the number (unit) of the open file
     & nspec,        ! the number of species
     & i,j
     
      real*8
     &r(nja),          ! transport code normalized small radius mesh 
c     &r_transport(nja),   ! transport code small radius mesh [cm]
c                       !    see include 'transport_prof.i' file
     &en(nja,nspeca),   !  density profiles [/cm**3]
     &temp(nja,nspeca), !  temperature profiles [keV]
     &zeff(nja)         !  effective charge profile

c-----function
      integer length_char1

c-----------------------------------------------------------------------
      integer i1p(2),itabl(3),k,i_total_n,i_total_t,i_fd,i_fdbeam
      double precision d2ene(nja),workk(iworkka), 
c     & d2work_1d(nja),work_1d(nja),tabl(3),rho_genray,h_rho
     & d2work_1d(nja),tabl(3),rho_genray,h_rho
   

c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus
      integer dims_id(2),start(2),count(2)
      integer nspec_id,nj_id,char256_id
      integer char256
      character*256 eqdsk_name
      character*128 name

      if (partner.eq.'genray_profs_in.txt') then
      
      igenray=2                            ! Choose a unit number
      open (unit=igenray,file='genray_profs_in.txt',status='UNKNOWN')
      read  (igenray,*)                    ! file ident data is skipped
      read  (igenray, 1003)  nspec,nj
 1003 format (10i5)
      
      elseif (partner.eq.'genray_profs_in.nc') then

      ncid = ncopn('genray_profs_in.nc',NCNOWRIT,istatus)
      write(*,*)'after ncopn ncid=',ncid,'istatus',istatus

c     read in dimension IDs

      nj_id = ncdid(ncid,'nj_dim',istatus)
      write(*,*)'partner: after ncdid nj_id',nj_id,'istatus',istatus
      nspec_id = ncdid(ncid,'nspecgr_dim',istatus)
      write(*,*)'partner: after ncdid nspec_id=',nspec_id,'istatus',
     .     istatus
      char256_id = ncdid(ncid,'char256dim',istatus)
      write(*,*)'partner: after ncdid char256_id=',char256_id,'istatus',
     .     istatus

c --- inquire about dimension sizes:#species,grid size---
      call ncdinq(ncid,nj_id,name,nj,istatus)
      write(*,*)'partner: after ncdinq, # of rad pts =',nj,'istatus=',
     .     istatus
      call ncdinq(ncid,nspec_id,name,nspec,istatus)
      write(*,*)'partner: after ncdinq, # of species =',nspec,' 
     .     istatus=',istatus
      call ncdinq(ncid,char256_id,name,char256,istatus)
      write(*,*)'partner: after ncdinq, length of eqdsk_name =',char256,
     .     ' istatus=',istatus

      start(1)=1
      start(2)=1
      count(1)=nj
      count(2)=nspec

      endif      !endif on partner

      nbulk_save=nbulk
      nbulk=nspec

      if (nspec.gt.nspeca) then
         write(*,*)'*************************************************'
         write(*,*)'in partner.f in read_transport_prof' 
         write(*,*)'nspec > nspeca. Need nspec.le.nspeca'
         write(*,*)'nspeca = ',nspeca,' nspec = ',nspec 
         write(*,*)'Please increase parameter nspeca in the genray'
         write(*,*)'file onewto_prof.i and recompile genray'
         write(*,*)'*************************************************'
         STOP 'in subroutine read_transport_prof'
      endif

      if (nspec.gt.nbulka) then
         write(*,*)'*************************************************'
         write(*,*)'in partner.f in read_transport_prof' 
         write(*,*)'nspec > nbulka. Need nspec.le.nbulka'
         write(*,*)'nbulka = ',nbulka,' nspec = ',nspec 
         write(*,*)'Please increase parameter nbulka in the genray'
         write(*,*)'file param.i and recompile genray'
         write(*,*)'*************************************************'
         STOP 'in subroutine read_transport_prof'
      endif

      if (nj.gt.nja) then
         write(*,*)'*************************************************'
         write(*,*)'in partner.f in read_transport_prof' 
         write(*,*)'nj > nja. but it should be nj.le.nja'
         write(*,*)'nja = ',nja,' nj = ',nj 
         write(*,*)'Please increase parameter nja in the file'
         write(*,*)'onewto_prof.i and recompile genray'
         write(*,*)'*************************************************'
         STOP 'in subroutine read_transport_prof'
      endif
      
      write(*,1004)  nj,nspec
 1004 format('partner: nj,nspec = ',2i5)


c-------------------------------------------------------------------------
c     Read rest of genray_profs_in
c-------------------------------------------------------------------------


      if (partner.eq.'genray_profs_in.txt') then
      
      read (igenray, 1001)  (r_transport(j), j=1,nj)     !cms

      read (igenray, 1001)  (zeff(j), j=1,nj)

      do i=1,nspec
         read (igenray, 1001) charge(i),dmass(i) !Units of electon charge,mass
         read (igenray, 1001)  (en(j,i), j=1,nj)  ! /cm**3
         read (igenray, 1001)  (temp(j,i), j=1,nj)  ! keV
      enddo

 1001 format (5(1pe16.9))
   
      close (unit =  igenray)


      elseif (partner.eq.'genray_profs_in.nc') then

      vid = ncvid(ncid,'eqdsk_name',istatus)
      call ncvgtc(ncid,vid,1,char256,eqdsk_name,char256,istatus)
      write(*,*)'partner: eqdsk_name = ',eqdsk_name
c     trimming length of eqdsk_name to first blank
      eqdsk_name=eqdsk_name(1:length_char1(eqdsk_name))
      write(*,*)'partner: trimmed eqdsk_name = ',eqdsk_name
      write(*,*)'partner: LEN_TRIM(eqdsk_name) = ',LEN_TRIM(eqdsk_name)
    
      vid = ncvid(ncid,'r',istatus)
      call ncvgt(ncid,vid,1,nj,r_transport,istatus) !normalized rho
      write(*,*)'partner: r_transport'
      write(*, 1001)  (r_transport(j), j=1,nj)

      vid = ncvid(ncid,'zeff',istatus)   
      call ncvgt(ncid,vid,1,nj,zeff,istatus) !normalized rho
      write(*,*)'partner: zeff'
      write(*, 1001)  (zeff(j), j=1,nj)

      vid = ncvid(ncid,'charge',istatus)   
      call ncvgt(ncid,vid,1,nspec,charge,istatus) !normalized rho
      vid = ncvid(ncid,'dmass',istatus)
      call ncvgt(ncid,vid,1,nspec,dmass,istatus) !normalized rho
      do i=1,nspec
         start(2)=i
         count(2)=1
         write(*,*)'partner: start=',start
         write(*,*)'partner: count=',count
         vid = ncvid(ncid,'en',istatus)   
         call ncvgt(ncid,vid,start,count,en(1,i),istatus) !normalized rho
         vid = ncvid(ncid,'temp',istatus)   
         call ncvgt(ncid,vid,start,count,temp(1,i),istatus) !normalized rho
      enddo
      start(2)=1
      count(2)=nspec

      call ncclos(ncid,istatus)

c-------------------------------------------------------------------------
c     overwrite eqdskin
c--------------------------------------------------------------------------

      eqdskin=eqdsk_name(1:LEN_TRIM(eqdsk_name))

      endif      !endif on second set of partner ifs


c--------------------------------------------------------------------------

      do i=1,nspec
         write(*,*)'partner: i,charge(i),dmass(i)'
         write(*,*)'partner: en(j,i), i= ',i,charge(i),dmass(i)
         write(*, 1001)  (en(j,i), j=1,nj)  !/cm**3
         write(*,*)'partner: temp(j,i), i= ',i
         write(*, 1001)  (temp(j,i), j=1,nj)  !keV
      enddo

c-------------------------------------------------------------------------
c     create normalized small radius: r
c--------------------------------------------------------------------------

      do j=1,nj
        r(j)= r_transport(j)/r_transport(nj)
      enddo
  
      write(*,*)' r',r

c------------------------------------------------------------------------
c     calculate the density and temperature profiles 
c     at the small radial Genray mesh
c------------------------------------------------------------------------

      i1p(1)=4
      i1p(2)=4
      itabl(1)=1
      itabl(2)=0
      itabl(3)=0

      do i=1,nspec

         call coeff1(nj,r,en(1,i),d2work_1d,i1p,1,workk)
         h_rho=1.d0/(ndens-1)        
         do k=1,ndens 
            rho_genray=h_rho*(k-1) ! GENRAY radial mesh points
            call terp1(nj,r,en(1,i),d2work_1d,rho_genray,1,tabl,itabl)
            dens1(k,i)=tabl(1)*1.d-13 !Convert to GENRAY units
         enddo                  !k
         
         call coeff1(nj,r,temp(1,i),d2work_1d,i1p,1,workk)
         do k=1,ndens 
            rho_genray=h_rho*(k-1) ! GENRAY radial mesh points
            call terp1(nj,r,temp(1,i),d2work_1d,rho_genray,1,tabl,itabl)
            temp1(k,i)=tabl(1)
         enddo                  !k

      enddo                     !i

c---------------------------------------------------------------------
c      calculate zeff1 profile at GENRAY small radius mesh
c---------------------------------------------------------------------
      
      izeff=2                   ! electron plasma for EC or EBW case,
                                ! zeff is given
      call coeff1(nj,r,zeff,d2work_1d,i1p,1,workk)
      do k=1,ndens 
         rho_genray=h_rho*(k-1) ! GENRAY radial mesh points
         call terp1(nj,r,zeff,d2work_1d,rho_genray,1,tabl,itabl)
         zeff1(k)=tabl(1)
      enddo                     !k
      
c---------------------------------------------------------------------
c     Set totp1 and vflow1 on the GENRAY small radius mesh
c---------------------------------------------------------------------

      do i=1,nbulk
         do k=1,ndens
            tpop1(k,i)=1.D0
         enddo
      enddo
      
      do i=1,nbulk
         do k=1,ndens
            vflow1(k,i)=0.D0
         enddo
      enddo

c---------------------------------------------------------------------
c     Reset nbulk=1 if nbulk_save=1 (for high frequency cases)
c---------------------------------------------------------------------

      if (nbulk_save.eq.1) nbulk=1      

      return
      end



      subroutine write_transport_prof(NR,nbulk,igenray,indexrho,
     &powden_e,powden_s,powden_i,currden,powtot_e,powtot_i,
     &powtot_s,powtott,currtot) 
c  -----------------------------------------------------------------------
c     Write data produced by ray-tracing code
c     in the file: genray_profs_out.nc and genray_profs_out.txt
c     using above subroutine inputs.
c  -----------------------------------------------------------------------
c
c     Variables  written to genray_profs_out include (nbulk.gt.1):
c     pgre_transport (j), j=1,nj    power density to electrons (W/cm**3)
c     pgri_transport (j), j=1,nj    power to (all) ion species (W/cm**3)
c     pgrc_transport (j), j=1,nj    RF current density (A/cm**2)
c     totgrp, total injected rf power (W)
c     totgrpe, totgrpi, totgrc, are respectively total power
c        to electron (W), to ions (W), and total RF current (Amps)
c     If (nbulk.eq.1) only powers to electrons and the currents are
c     written out.

c     Input variables ending in "s" are for specific ion species. Ion
c     results are in entries 2:nbulk.  The first entry in the
c     arrays is not used [dimensioning is for programming convenience].
c     Output variables to genray_profs_out(.nc) ending in "s" are
c     given for ions: i=1,nbulk-1 (for nbulk.gt.1).

c      implicit none

      include 'param.i'    !gives parameters nra,nbulka

      include 'transport_prof.i'
               !       in transport_prof.i':
               !nspeca is a maximal number of transport code species
               !nja   is a maximal number of points in 
               !      transport small radial mesh 
               !nj    is number of points in 
               !      transport small radial mesh 
               !r_transport(nj) is transport small radial mesh [cm]
    
c-----inputs to subroutine:
      integer
     &nbulk,        ! number of species,(electrons are the first) 
                    ! including possible hot beam species
     &igenray       ! the number (unit) of the open file

      real*8
     &powden_e(NRA),      !power density profile to electrons [erg/(sec*cm**3)]
     &powden_s(NRA,nbulka),!power density profile to each ion [erg/(sec*cm**3)]
     &powden_i(NRA),      !summed pwr density profile to ions [erg/(sec*cm**3)]
     &currden(NRA),       !current density profile (along B-field) [A/cm**2]
     &powtot_e,           !total power to electrons            [erg/sec]
     &powtot_i,           !total power to ions                 [erg/sec]
     &powtot_s(nbulka),  !total power to each ion component   [erg/sec]
     &powtott,            !total injected power
     &currtot             !total current along B field         [A]
c-----outputs for transport code via genray_profs_out:
      real*8
     &pgre_transport(nja),       ! power density to electrons    [W/cm**3] 
     &pgri_transport(nja),       ! power density to ions         [W/cm**3] 
     &pgrs_transport(nja,nbulkma),! power density to ion specie  [W/cm**3]
     &pgrc_transport(nja),       ! current density               [A/cm**2]
     &totgrpe,                   ! total power to electrons      [W]
     &totgrpi,                   ! total power to ions           [W]
     &totgrps(nbulkma),          ! total power to each ion specie[W]
     &totgrp,                    ! total injected rf power
     &totgrc                     ! total current                 [A]

c-----locals:
      !*******for the spline functions******
      integer iworkka 
      parameter(iworkka=3*NRA+1) ! is the length of the work array  
      integer i1p(2),itabl(3)
      double precision workk(iworkka), 
c     & d2work_1d(NRA),work_1d(NRA),tabl(3)
     & d2work_1d(NRA),tabl(3)

      integer n,j
      
      double precision
     &h_rho,                 ! genray grid step
     &r_genray(NRA),         ! GENRAY normalized rho grid used for
                             ! power profiles
     &rho_transport(nja)     ! normalized transport radial mesh

c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus
      character filenc*128,ltitle*128
      integer dims_id(2),start(2),dimsm_id(2),countm(2)
                                           !countm,dimsm_id for ions (no e)
      integer nbulk_id,nprofs_id,nbulkm_id !nbulkm_id for ions (no e) 


ctest transport mesh and nj to run this subroutine without given transport grid
c      nj=21
c      do j=1,nj
c         r_transport(j)=(j-1)/dfloat(nj-1)
c      enddo
cendtest
     
c      write(*,*)'before open (unit = igenray, file = genray_profs_out'

      open (unit = igenray,file='genray_profs_out.txt',status='UNKNOWN')

c     Set value of nj=nr, if this is a case where it is not been read 
c     in by read_transport_prof.
      if (nj.eq.0) then
         nj=nr
         dr_transp=1.0d0/(nj-1)
         do j=1,nj
            r_transport(j)=(j-1)*dr_transp
         enddo
      endif

      write   (igenray, '(3i5)')   nbulk,nj,indexrho  
                                   !#species,grid size,radial coord type


c-----------------------------------------------------------------
c     GENRAY small radial mesh
c------------------------------------------------------------------
      h_rho=1.d0/dfloat(NR-1)
      do j=1,NR-1
         r_genray(j)=h_rho*(j-0.5d0)
      enddo
c      write(*,*)'partner.f r_genray',r_genray 


c--------------------------------------------------------------------
c     calculate the power density to electrons
c     at the small radial transport mesh
c     Linear interpolation is used below, rather than splines
c     which can create problems (negative power) in radial
c     regions of rapid change.
c---------------------------------------------------------------------

      do j=1,nj
        rho_transport(j)=r_transport(j)/r_transport(nj)! normalized transport 
      enddo                                ! code radial mesh points
      
      write(*,*)'rho_transport',rho_transport
      write(*,*)'powden_e',powden_e
      write(*,*)'partnenr.f in write_transport_prof'
      write(*,*)'before lin_inter(ppowden_e) NR',NR,'nj',nj 
      call lin_interp(r_genray,powden_e,NR-1,rho_transport,
     +     pgre_transport,nj)
      write(*,*)'after lin_inter(ppowden_e)'

      do j=1,nj
         pgre_transport(j)=pgre_transport(j)*1.d-7 !W/cm**3
      enddo
      write(*,*)'pgre_transport',pgre_transport

c--------------------------------------------------------------------
c     calculate the power density to ions
c     on the minor radius radial transport mesh
c---------------------------------------------------------------------
      
c      write(*,*)'powden_i',powden_i
      if (nbulk.gt.1) then

         write(*,*)'nbulk,gt.1 before lin_inter(ppowden_i)NR',NR,'nj',nj 
         call lin_interp(r_genray,powden_i,NR-1,rho_transport,
     +        pgri_transport,nj)

         write(*,*)'after lin_inter(ppowden_i)'
 
         do j=1,nj
            pgri_transport(j)=pgri_transport(j)*1.d-7 !W/cm**3
         enddo
      endif
c
c--------------------------------------------------------------------
c     calculate the power density to the different ion components
c     on the small radial transport mesh
c---------------------------------------------------------------------
      if (nbulk.gt.1) then
         do k=2,nbulk
            write(*,*)'before lin_inter(ppowden_i_s) k',k,'nj',nj 
            call lin_interp(r_genray,powden_s(1,k),NR-1,rho_transport,
     +           pgrs_transport(1,k-1),nj)
            write(*,*)'after lin_inter(ppowden_s)'
 
            do j=1,nj
               pgrs_transport(j,k-1)=pgrs_transport(j,k-1)*1.d-7 !W/cm**3
            enddo
            
         enddo                  !k
      endif

c      do k=2,nbulk
c         write(*,*)'partner:k,pgrs_transport(j,k-1)= ',k,
c     .        (pgrs_transport(j,k-1),j=1,nj)
c      enddo

c--------------------------------------------------------------------
c     calculate the current density profile 
c     on the small radial transport mesh
c---------------------------------------------------------------------
      write(*,*)'before lin_interp currden NR',NR,'nj',nj 
      call lin_interp(r_genray,currden,NR-1,rho_transport,
     +     pgrc_transport,nj)

      write(*,*)'after lin_interp currden'
c--------------------------------------------------------------------
c     write the profiles in genray_prof_out 
c     on the small radial transport mesh
c---------------------------------------------------------------------
c
      write (igenray,   3010 )  (r_transport(j), j=1,nj)
      write (igenray,   3010 )  (pgre_transport(j), j=1,nj) ! power to e
                                                            ! W/cm**3 

c-----current  density (along B field) profile
      write (igenray,   3010 )  (pgrc_transport(j), j=1,nj) ! A/cm**2

      totgrp  = 1.d-7*powtott  ! total injected power        [W]
      totgrpe = 1.d-7*powtot_e ! total power to electrons    [W]
      totgrc  = currtot        ! total current               [A]
      write (igenray, 3020) totgrp, totgrpe, totgrc

      if (nbulk.gt.1) then
         do k=2,nbulk
            totgrps(k-1)= 1.d-7*powtot_s(k)! total power to ion species[W]
         enddo
         write (igenray, 3010 )  (pgri_transport(j), j=1,nj)
                                ! power to ions, W/cm**3 
         do k=2,nbulk
            write (igenray, 3010 ) (pgrs_transport(j,k-1),j=1,nj) !W/cm**3
         end do
         totgrpi = 1.d-7*powtot_i ! total power to ions         [W]
         write (igenray, 3010) totgrpi, totgrps(1:nbulk-1)
      endif
                                                            

 3010 format (5(1pe19.6))
 3020 format (4(1pe19.6))

      close (unit =  igenray)

c
c-----------------------------------------------------------------
c     Open netCDF file for output
c------------------------------------------------------------------
      filenc='genray_profs_out.nc'
      ncid=nccre(filenc,NCCLOB,istatus)
      call check_err(istatus)
      write(*,*)'In write_transport_prof after nccre, istatus=',istatus

c     Brief description to be added to file:
      ltitle='Power and current profile data passed from GENRAY'
      call ncaptc(ncid,NCGLOBAL,'title',NCCHAR,length_char(ltitle),
     +     ltitle,istatus)


c-----------------------------------------------------------------------
c     Set up and write data to netcdf file
c-----------------------------------------------------------------------

      start(1)=1
      start(2)=1

      countm(1)=nj
      countm(2)=nbulk-1

cl    Define dimensions for netcdf file
      nprofs_id=ncddef(ncid,'nprofs_dim',nj,istatus)
      nbulk_id=ncddef(ncid,'nbulk_dim',nbulk,istatus)
      if (nbulk.gt.1) 
     .     nbulkm_id=ncddef(ncid,'nbulkm_dim',nbulk-1,istatus)

c     Define vector of dimensions
      dims_id(1)=nprofs_id
      dims_id(2)=nbulk_id
      if (nbulk.gt.1) then
         dimsm_id(1)=nprofs_id
         dimsm_id(2)=nbulkm_id
      endif

c     Define variable names for netcdf file

c$$$Don't Need
c$$$      vid=ncvdef(ncid,'nj',NCLONG,0,0,istatus)
c$$$      call ncaptc(ncid,vid,'long_name',NCCHAR,24,
c$$$     +           'Dimension of radial mesh',istatus)
c$$$      call check_err(istatus)

      vid=ncvdef(ncid,'indexrho',NCLONG,0,0,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,35,
     +           'Radial coordinate type: should be 2',istatus)
      call check_err(istatus)

c$$$Don't Need
c$$$      vid=ncvdef(ncid,'nbulk',NCLONG,0,0,istatus)
c$$$      call ncaptc(ncid,vid,'long_name',NCCHAR,40,
c$$$     +           'Number of plasma species, incl electrons',istatus)
c$$$      call check_err(istatus)
c$$$


      vid=ncvdef(ncid,'rgenray',NCDOUBLE,1,nprofs_id,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,18,
     +           'genray radial mesh',istatus)
      call ncaptc(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'pgre',NCDOUBLE,1,nprofs_id,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,39,
     +           'FSA power density to electrons vs radius',istatus)
      call ncaptc(ncid,vid,'units',NCCHAR,11,
     +           'Watts/cm**3',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'totgrpe',NCDOUBLE,0,0,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,29,
     +           'Integrated power to electrons', istatus)
      call ncaptc(ncid,vid,'units',NCCHAR,5,
     +           'Watts',istatus)
      call check_err(istatus)

      if (nbulk.gt.1) then
         vid=ncvdef(ncid,'pgri',NCDOUBLE,2,dimsm_id,istatus)
         call ncaptc(ncid,vid,'long_name',NCCHAR,38,
     +        'FSA power density to each ion vs radius',istatus)
         call ncaptc(ncid,vid,'units',NCCHAR,11,
     +        'Watts/cm**3',istatus)
         call check_err(istatus)
         
         vid=ncvdef(ncid,'pgrit',NCDOUBLE,1,nprofs_id,istatus)
         call ncaptc(ncid,vid,'long_name',NCCHAR,39,
     +        'FSA power density to all ions vs radius',istatus)
         call ncaptc(ncid,vid,'units',NCCHAR,11,
     +        'Watts/cm**3',istatus)
         call check_err(istatus)
         
         vid=ncvdef(ncid,'totgrps',NCDOUBLE,1,nbulkm_id,istatus)
         call ncaptc(ncid,vid,'long_name',NCCHAR,36,
     +        'Integrated power to each ion species', istatus)
         call ncaptc(ncid,vid,'units',NCCHAR,5,
     +        'Watts',istatus)
         call check_err(istatus)
         
         vid=ncvdef(ncid,'totgrpi',NCDOUBLE,0,0,istatus)
         call ncaptc(ncid,vid,'long_name',NCCHAR,28,
     +        'Integrated power to all ions', istatus)
         call ncaptc(ncid,vid,'units',NCCHAR,5,
     +        'Watts',istatus)
         call check_err(istatus)
      endif

      vid=ncvdef(ncid,'totgrp',NCDOUBLE,0,0,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,37,
     +     'Total rf injected power in this mode', istatus)
      call ncaptc(ncid,vid,'units',NCCHAR,5,
     +     'Watts',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'pgrc',NCDOUBLE,1,nprofs_id,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,28,
     +           'Rf current density <j.B/B_0>',istatus)
      call ncaptc(ncid,vid,'units',NCCHAR,10,
     +           'Amps/cm**2',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'totgrc',NCDOUBLE,0,0,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,23,
     +           'Total rf driven current', istatus)
      call ncaptc(ncid,vid,'units',NCCHAR,4,
     +           'Amps',istatus)
      call check_err(istatus)

c     End the define-mode and start the data-mode

      call ncendf(ncid,istatus)
      call check_err(istatus)


c     Write data into netcdf file

      vid=ncvid(ncid,'indexrho',istatus)
      call ncvpt(ncid,vid,1,1,indexrho,istatus)

      vid=ncvid(ncid,'rgenray',istatus)
      call ncvpt(ncid,vid,start(1),countm(1),r_transport,istatus)

      vid=ncvid(ncid,'pgre',istatus)
      call ncvpt(ncid,vid,start(1),countm(1),pgre_transport,istatus)

      vid=ncvid(ncid,'totgrpe',istatus)
      call ncvpt(ncid,vid,1,1,totgrpe,istatus)

      if (nbulk.gt.1) then
         vid=ncvid(ncid,'pgri',istatus)
         do k=1,nbulk-1
            start(2)=k
            countm(2)=1
            call ncvpt(ncid,vid,start,countm,pgrs_transport(1,k),
     .           istatus)
         enddo
         start(2)=1
         countm(2)=nbulk-1

         vid=ncvid(ncid,'pgrit',istatus)
         call ncvpt(ncid,vid,start(1),countm(1),pgri_transport,istatus)
         
         vid=ncvid(ncid,'totgrps',istatus)
         call ncvpt(ncid,vid,start(2),countm(2),totgrps,istatus)
         
         vid=ncvid(ncid,'totgrpi',istatus)
         call ncvpt(ncid,vid,1,1,totgrpi,istatus)
      endif

      vid=ncvid(ncid,'totgrp',istatus)
      call ncvpt(ncid,vid,1,1,totgrp,istatus)

      vid=ncvid(ncid,'pgrc',istatus)
      call ncvpt(ncid,vid,start(1),countm(1),pgrc_transport,istatus)

      vid=ncvid(ncid,'totgrc',istatus)
      call ncvpt(ncid,vid,1,1,totgrc,istatus)


c     Close netcdf file

      call ncclos(ncid,istatus)
      call check_err(istatus)


      return
      end

 
