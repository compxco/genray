cSAP090720
c      subroutine netcdfr3d(netcdfnm,iy,iy_,jx,lrz,y,x,rya,vnorm,f)
      subroutine netcdfr3d

c-----It reads the distribution and the mesh from netcdfnm.nc file

      implicit none
      !implicit integer (i-n), real*8 (a-h,o-z)
c       save

      
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'
      
c-----determine the dimensions 
c     parameter (iya,jxa,lrza,ngena)
      include 'param.i'
      include 'one_nml.i' 
      include 'dskin.i'
c-----output
c     It will put output data to dskin.i and one_nml.i files
c      real*8 vnorm          ->  dskin.i
c      integer jx,lrz        ->  one_nml.i
c      integer iy_(iya)=iy_(iya) ->  dskin.i
c      real*8 x(jx)          ->  dskin.i
c      real*8 y(iya,lrz)     ->  dskin.i
c      real*8 rya(lrza)=rovera(lrz) ->dskin.i     !normalized small radius
c      real*8 f(iya,jxa,lrza,ngen) -> dskin.i
      
c-----local
      integer lrzmax,k
c --- some stuff for netCDF file ---
      character*128 name
      integer ncid,istatus,ncvdef2,ncdid2,ncddef2
      integer xdim,ydim,rdim,kdim,vid
      integer ntotal 
      integer ll,j,i
      integer start(3),count(3),start_y(2),count_y(2)

cSAP090720
      integer iylrz2a,ijl,iyjxlrz3a
c      parameter (iylrz2a=iya*lrza) 
c      real*8 tem1(iylrz2a)
c      parameter (iyjxlrz3a=iylrz2a*jxa) 
c      real*8 tem2(iyjxlrz3a)

      real*8, pointer :: tem1(:),tem2(:)
      integer istat
      
      data start/1,1,1/,start_y/1,1/
 
c     Open previous netCDF file
c      write(*,*)'netcdfr3d 1 iya,jxa,lrza',iya,jxa,lrza
c      write(*,*)'before ncopn netcdfnm=',netcdfnm
c      ncid = ncopn(netcdfnm,NCNOWRIT,istatus) 
c     Open previous netCDF file

c      write(*,*)'netcdfr3d 1 iya,jxa,lrza',iya,jxa,lrza
c      write(*,*)'before ncopn netcdfnm=',netcdfnm
      call ncopn2(netcdfnm,NCNOWRIT,ncid,istatus)
c      write(*,*)'after ncopn ncid=',ncid,'istatus',istatus
c      write(*,*)'after ncopn ncid=',ncid,'istatus',istatus
c.......................................................................
c     read in dimension IDs and sizes

c      write(*,*)'before ncdid xdim'
      xdim =ncdid2(ncid,'xdim',istatus)
c      write(*,*)'after ncdid xdim=',xdim,'istatus',istatus
      ydim =ncdid2(ncid,'ydim',istatus)
c      write(*,*)'after ncdid ydim=',ydim,'istatus',istatus
      rdim =ncdid2(ncid,'rdim',istatus)
c      write(*,*)'after ncdid rdim=',rdim,'istatus',istatus
      kdim=ncdid2(ncid,'species_dim',istatus)
c      write(*,*)'afterncid kim=',kdim,'istatus',istatus

c --- inquire about dimension sizes ---
c     ncdinq2(netCDF_id, dimension_id_from_ncdid, returned_dim_name,
c     returned_dim_size)
c     Note: for unlimited dimension, returned_dim_size=current maximum
c     which is the same as the maximum record number



c      write(*,*)'netcdfr3d.f before ncdinq2'
c      write(*,*)'jx,lrz,ntotal,ngen,iya',jx,lrz,ntotal,ngen,iya
         
      call ncdinq2(ncid,ydim,name,iy,istatus)
      call ncdinq2(ncid,xdim,name,jx,istatus)
      call ncdinq2(ncid,rdim,name,lrz,istatus)
      call ncdinq2(ncid,kdim,name,ntotal,istatus) 
        
      call ncvid2(vid,ncid,'ngen',istatus)   
      call ncvgt_int2(ncid,vid,1,1,ngen,istatus)
c      write(*,*)'after ncvgp ngen=',ngen
 
c      write(*,*)'after ncdinq2 iy,jx,lrz,ntotal,ngen',
c     &                        iy,jx,lrz,ntotal,ngen
cSAP090808 delete jxa,lrza
c      write(*,*)'netcdfr3d 2 iya,jxa,lrza',iya,jxa,lrza
c      write(*,*)'netcdfr3d 2 iya',iya
c-----check the dimensions
 

c      if (lrz.gt.lrza) then
c         write(*,*)'netcdfr3d lrz.gt.lrza' 
c         write(*,*)'lrz from netcdfnm.nc file =',lrz
c         write(*,*)'lrza from param.i file =',lrza
c         write(*,*)'Attention!!! It should be lrz=lrza'
c         write(*,*)'Please change lrza in param.i and recomplile codes'
c         stop
c      endif

c      write(*,*)'netcdfr3d 3 iya,jxa,lrza',iya,jxa,lrza
c      write(*,*)'netcdfr3d 3 iy,iy_,jx,lrz',iy,iy_,jx,lrz

c       if (jx.gt.jxa) then
c         write(*,*)'netcdfr3d jx.gt.jxa'
c         write(*,*)'jx from netcdfnm.nc file =',jx
c         write(*,*)'jxa from param.i file =',jxa
c         write(*,*)'Please change jxa in param.i'
c         stop
c      endif

cSAP090720
c      if (iy.gt.iya) then
      if (iy.gt.iya) then
         write(*,*)'netcdfr3d iy.gt.iya'
         write(*,*)'iy from netcdfnm.nc file =',iy
         write(*,*)'iya from param.i file =',iya
         write(*,*)'Please change iya in param.i'
         stop
      endif

c--------------------------------------------------
c     allocate pointers in dskin.i
c--------------------------------------------------
c      write(*,*)'netcdf3d before ainalloc_dskin_i'
      call ainalloc_dskin_i
c      write(*,*)'netcdf3d after ainalloc_dskin_i'
c---------------------------------------------------
c     allocate pointers tem1 and tem2
c--------------------------------------------------
      iylrz2a=iya*lrz 
      allocate( tem1(1:iylrz2a),STAT=istat)
      call bcast(tem1,0.d0,SIZE(tem1)) 

      iyjxlrz3a=iylrz2a*jx
      allocate( tem2(1:iyjxlrz3a),STAT=istat)
      call bcast(tem2,0.d0,SIZE(tem2)) 
c----------------------------------
      count(1)=iy        
      count(2)=jx
      count(3)=lrz

c-----normalized momentum x (momentum/mass/vnorm) variables

c     vnorm - character velocity (momentum-per-mass)[cms/sec]
      call ncvid2(vid,ncid,'vnorm',istatus)   
      call ncvgt_doubl2(ncid,vid,1,1,vnorm,istatus)
c      write(*,*)'after ncvgp vnorm=',vnorm 

      call ncvid2(vid,ncid,'x',istatus)
      call ncvgt_doubl2(ncid,vid,1,jx,x,istatus)
c      write(*,*)'netcdfr3d x',(x(j),j=1,jx)

c      call ncvid2(vid,ncid,'dx',istatus)
c      call ncvgt_doubl2(ncid,vid,1,jx,dx,istatus)     
c      write(*,*)'dx',(dx(j),j=1,jx)

c      call ncvid2(vid,ncid,'cint2',istatus)
c      call ncvgt_doubl2(ncid,vid,1,jx,cint2,istatus)     
c      write(*,*)'cint2',(cint2(j),j=1,jx)

c-----pitch angle variables y

cSAP090720
c      call ncvid2(vid,ncid,'iy_',istatus)
c      call ncvgt_doubl2(ncid,vid,1,lrz,iy_,istatus)
c      write(*,*)'netcdfr3d iy_',(iy_(ll),ll=1,lrz)

      call ncvid2(vid,ncid,'iy_',istatus)
      call ncvgt_int2(ncid,vid,1,lrz,iy_,istatus)
c      write(*,*)'netcdfr3d iy_',(iy_(ll),ll=1,lrz)

cSAP090720
      count_y(1)=iy
      count_y(2)=lrz 
c      write(*,*)'netcdfr3d iy,lrz',iy,lrz

      call ncvid2(vid,ncid,'y',istatus)
c      call ncvgt_doubl2(ncid,vid,start,count_y,y,istatus)
cSm040111
      call ncvgt_doubl2(ncid,vid,start,count_y,tem1,istatus)
 
c      do i=1,iylrz2a
c         write(*,*)'i,tem1(i)',i,tem1(i)
c      enddo

      j=0
      do ll=1,lrz
cSAP090720
c         do i=1,iy_(ll)
         do i=1,iy_(ll)
            j=j+1
            y(i,ll)=tem1(j)
         enddo
      enddo
 
c      do ll=1,lrza
c         write(*,*)'ll=',ll,'iy_(ll)=',iy_(ll)              
c         write(*,*)'y',(y(i,ll),i=1,iy_(ll))
c      enddo
c      stop 'netcdfr3d'

c      call ncvid2(vid,ncid,'dy',istatus)
c      call ncvgt_doubl2(ncid,vid,start,count_y,dy,istatus)
c      do ll=1,lrza
c         write(*,*)'ll=',ll,'iy_(ll)=',iy_(ll)              
c         write(*,*)'dy',(dy(i,ll),i=1,iy_(ll))
c      enddo

c      call ncvid2(vid,ncid,'dy',istatus)
c      call ncvgt_doubl2(ncid,vid,start,count_y,cynt2,istatus)
c      do ll=1,lrza
c         write(*,*)'ll=',ll,'iy_(ll)=',iy_(ll)              
c         write(*,*)'cynt2',(cynt2(i,ll),i=1,iy_(ll))
c      enddo

c-----normalized small radius
cSAP090720
      call ncvid2(vid,ncid,'rya',istatus)
      call ncvgt_doubl2(ncid,vid,1,lrz,rya,istatus)              
c      write(*,*)'rya',(rya(ll),ll=1,lrz)
      rovera=rya
c      write(*,*)'rovera',(rovera(ll),ll=1,lrz)
c-----distribution function f(i,j,ll) [vnorm**3/(cm**3*(cm/sec)**3)]
      call ncvid2(vid,ncid,'f',istatus)
c      call ncvgt_doubl2(ncid,vid,start,count,f,istatus)
cSm04011
c      write(*,*)'netcdr3d count',count
      call ncvgt_doubl2(ncid,vid,start,count,tem2,istatus)
c      write(*,*)'netcdr3d read f istatus',istatus
ctest
c      write(*,*)'iyjxlrz3a',iyjxlrz3a
c     ijl=0
c      do ll=1,lrz
c         do j=1,jx
c            do i=1,iy_(ll)
c               ijl=ijl+1
c               write(*,*)'ll,j,i,ijl,tem2(ijl)',ll,j,i,ijl,tem2(ijl)
c            enddo
c         enddo
c      enddo  
      
c      do ll=1,lrz
c         do j=1,jx
c            do i=1,iy_(ll)
c               write(*,*)'ll,j,i,f(i,j,ll,1)',ll,j,i,f(i,j,ll,1)
c            enddo
c         enddo
c      enddo  
c endtest
       ijl=0
       do ll=1,lrz
          do j=1,jx
             do i=1,iy_(ll)
                ijl=ijl+1
cSAP090720
c               f(i,j,ll)=tem2(ijl)
                 f(i,j,ll,1)=tem2(ijl)
c                write(*,*)'i,j,ll,ijl,f(i,j,ll,1)',
c     &                     i,j,ll,ijl,f(i,j,ll,1)
            enddo
         enddo
      enddo

c      do ll=1,lrz
c         write(*,*)' netcdfr3 ll=',ll
c         do j=1,jx
c            write(*,*)'ll=',ll,'j=',j,'f(i=1,...,iy_(ll))'
c            write(*,*)(f(i,j,ll,1),i=1,iy_(ll))   
c            do i=1,iy_(ll)
c              write(*,*)'i,j,ll,f(i,j,ll,1)',
c     &                   i,j,ll,f(i,j,ll,1)
c            enddo  
c         enddo
c      enddo
c      write(*,*)'netcdfr3d after read f'
c      stop 'netcdfr3d test f'
c-----Close netCDF file
      call ncclos2(ncid,istatus)
      call check_err(istatus)

c-----normalization of the distribution function 
c     for 10**13(1/cm**3) density
cSAP090720
      do k=1,lrz 
         do j=1,jx 
           do i=1,iy_(k)
c         do i=1,iya
c            do j=1,jx
              f(i,j,k,1)=f(i,j,k,1)*1.d-13
c              write(*,*)'i,j,k,f(i,j,k,1)',i,j,k,f(i,j,k,1)
            enddo
         enddo
      enddo

      deallocate (tem1,STAT=istat)
      deallocate (tem2,STAT=istat)

      return
      end
       
     
      subroutine wrtnetcdf(kopt)
      !implicit integer (i-n), double precision (a-h,o-z)
      implicit none
c
c     Write ray tracing data (as in mnemonic.txt) into a netCDF file.

cSm0940727
c     If the parameter ionetwo.eq.1 it will write the current
c     and power profiles into a netcdf file

      include 'param.i'
      include 'writencdf.i'
      include 'one.i'
      include 'ions.i'
      include 'adj.i'
      include 'cone_nml.i'     !nccone
      include 'grill_nml.i'    !ngrill
c--------------------------
cSm040727 to write the current and power profiles into netcdf file
c     Done in subroutine, wrtnetcdf_prof
c      include 'onetwo.i'
c--------------------------
c-----input
      integer
     & kopt  !1 create netCDF file and define dimensions,variables
c            !and attributes
             !0 Write data into created netcdf file
     
c-----local     
      integer n1n2,istat
c     Storage tem1 is used in netcdf writes, including complex numbers.
cSAP090905
c      parameter (n1n2=2*nrelta*nraya)
c      real*8 tem1(n1n2) 
       real*8, pointer :: tem1(:)
    
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,nccre2,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid,twoid,char64id,char128id,char256id,
cSAP
     &char8id
    
      integer nbulkmid  !ion species dimension
      integer nbulkid
      integer nrho_adjid

      integer ray_dims(3),start(3),ray_count(3)
      integer ray_dimss(3),starts(3),ray_counts(3)

      double complex cei

      character ltitle*256

      integer neltmax,iray,i,j,ii,ll
      integer length_char

      data start/1,1,1/,starts/1,1,1/
       
      save

cSAP090903
      n1n2=2*nrelta*nrayl
c      write(*,*)'nrelta,nrayl,n1n2',nrelta,nrayl,n1n2
c------------------------------------------
c     allocate pointers tem1
c-------------------------------------------
      allocate( tem1(1:n1n2),STAT=istat)
      call bcast(tem1,0.d0,SIZE(tem1)) 
c      write(*,*)'wrtnetcdf after allocate tem1 istat=',istat
c------------------------------------------
      cei=(0.d0,1.d0)

c     Maximum number of ray elements per ray:
      neltmax=0

c      write(*,*)'wrtnetcdf nray=',nrayl


      do iray=1,nrayl
         neltmax=max(neltmax,nrayelt_nc(iray))
c         write(*,*)'wrtnetcdf iray,nrayelt_nc(iray),neltmax',
c     &                        iray,nrayelt_nc(iray),neltmax
      enddo

cSm05038
      if (neltmax.eq.0) neltmax=1

      ray_count(1)=neltmax
      ray_count(2)=nrayl

      ray_count(3)=2
      ray_counts(1)=neltmax
      ray_counts(2)=nrayl
      ray_counts(3)=1

c      write(*,*)'ray_count',ray_count

c.......................................................................
cl    1. Initialize part, creating new netcdf file
c

c --- begin if ---
      if ( kopt.eq.1 ) then
c      write(*,*)'wrtnetcdf,  kopt=',kopt

C-----------------------------------------------------------------------
c
cl     1.1 create netCDF file and define dimensions,variables
c          and attributes
c

c.......................................................................
cl    1.1.1 create netCDF filename (Entering define mode.)
c     integer function nccre(filename,overwrite?,error_code)
c     Ref to page 46 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

c      write(*,*)'wrtnetcdf before  ncid=nccre(filenc,NCCLOB,istatus)'

      ncid=nccre2(trim(filenc),NCCLOB,istatus)
      call check_err(istatus)
c      write(*,*)'ncid',ncid      
c     Brief description added to file:
      ltitle='netCDF file of ray data from GENRAY version: '//version
      if( length_char(trim(ltitle)).gt.256 )
     +   stop 'Adjust ltitle in netcdfrw2'
      call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,
     +  length_char(trim(ltitle)), trim(ltitle), istatus)
        
c.......................................................................
cl    1.1.2 define dimensions
c     integer function ncddef(ncid,dim_name,dim_siz,error_code)
c       returns dimension id.
c     p. 67 of netcdf-2 manual

c     For ray data:
c      write(*,*)'wrtnetcdf before ncddef(neltmax) neltmax=',neltmax

      neltmaxid=ncddef2(ncid,'neltmax',neltmax,istatus)         
     
c      write(*,*)'wrtnetcdf before ncddef(ncid,nrays) nrayl=',nrayl

      nraysid=ncddef2(ncid,'nrays',nrayl,istatus)

c      write(*,*)'wrtnetcdf before ncddef(ncid,two,2,istatus)'
      
      twoid=ncddef2(ncid,'two',2,istatus)
    
c      write(*,*)'wrtnetcdf before ncddef(ncid,nbulk)'

      nbulkid=ncddef2(ncid,'nbulk',nbulk,istatus) 

cSAP080303
c      write(*,*)'wrtnetcdf before ncddef(char8dim)'
      char8id=ncddef2(ncid,'char8dim',8,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char64dim)'

      char64id=ncddef2(ncid,'char64dim',64,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char128dim)'

      char128id=ncddef2(ncid,'char128dim',128,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char5212dim)'

      char256id=ncddef2(ncid,'char256dim',256,istatus)

c      write(*,*)'neltmaxid',neltmaxid
c      write(*,*)'nraysid',nraysid
c      write(*,*)'twoid',twoid

      ray_dims(1)=neltmaxid
      ray_dims(2)=nraysid
      ray_dims(3)=twoid
    
      if (nbulk.gt.1) then

c         write(*,*)'wrtnetcdf before ncdef(nbulkm,nbulk-1 )'

         nbulkmid=ncddef2(ncid,'nbulkm',nbulk-1,istatus)

         ray_dimss(1)=ray_dims(1)
         ray_dimss(2)=ray_dims(2)
         ray_dimss(3)=nbulkmid
      endif
      
c.......................................................................
cl    1.1.3 define variables

c     Note, the variable IDs (denoted below as vid) are
c     not saved here in this subroutine; rather, the IDs
c     are retrieved from the netCDF data file, as needed,
c     by calling the netCDF routine ncvid.

c     netCDF variable_define:
c     integer function ncvdef(ncid,variable_name,variable_type,
c                number_of_dimensions,
c                vector_for_length_per_dimension,error_code)
c       returns varid.
c     Note: Unlimited dimension must be last
c     Refer to p. 77, netcdf-2 manual

c     NCDOUBLE for REAL*8.  This is independent of the
c       internal representation of data on a particular
c       machine, e.g., 32- or 64-bit based calculations.


c--------------------------
c     Genray version
c--------------------------

c      write(*,*)'before ncvdef2(version)'
      vid=ncvdef2(ncid,'version',NCCHAR,1,char64id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +          'GENRAY version number',istatus)

c--------------------------
c     Mnemonic for the run
c--------------------------
c      write(*,*)'before ncvdef2(mnemonic)'
      vid=ncvdef2(ncid,'mnemonic',NCCHAR,1,char128id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +          'Mnemonic run identifier',istatus)

c-------------------------------
c     Run Descriptive Parameters
c-------------------------------
c      write(*,*)'before ncvdef2(vid)'
      vid=ncvdef2(ncid,'id',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +            'Disp relation identifier',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(iabsorp)'
      vid=ncvdef2(ncid,'iabsorp',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Absorp calc identifier',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ieffic)'
      vid=ncvdef2(ncid,'ieffic',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29,
     +            'Current drive calc identifier',istatus)
      call check_err(istatus)

cSAP080303
c      write(*,*)'before ncvdef2(ion_absorption)'
      vid=ncvdef2(ncid,'ion_absorption',NCCHAR,1,char8id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +'Switch on/off ion absorption at iabsorp=3,9,91,92',istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(refl_loss)'
      vid=ncvdef2(ncid,'refl_loss',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +    'fraction of power loss at each reflection',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(iflux)'
      vid=ncvdef2(ncid,'iflux',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Flux calc, non-Westerhof-Tokman id',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ioxm)'
      vid=ncvdef2(ncid,'ioxm',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     +            'wave mode indicator (1 - om, -1 - xm )',istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(ioxm_n_npar)'
      vid=ncvdef2(ncid,'ioxm_n_npar',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,68,
     +'wave mode indicator: sign before square root to find '//
     +' N(N_parallel)',
     +istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(jwave)'
      vid=ncvdef2(ncid,'jwave',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37,
     +            'Wave harmonic, for CD efficiency calc',istatus)
      call check_err(istatus)

      if (iabsorp.eq.4) then
      vid=ncvdef2(ncid,'i_im_nperp',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     +       'iabsorp=4  nperp:1,ImD_full/dD/dnp; 2,Cmplx soln',istatus)
      call check_err(istatus)
      endif

c      write(*,*)'before ncvdef2(i_geom_optic)'
      vid=ncvdef2(ncid,'i_geom_optic',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Integrate rays wrt (1)time,(2)dist',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(istart)'
      vid=ncvdef2(ncid,'istart',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,50,
     +     'Ray launch type: 1,eccone; 2,grill, 3,OX in plasma',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ncone)'
      vid=ncvdef2(ncid,'ncone',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,39,
     +     'Number of rf source ray cones, istart=1',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,40,
     +     'Each cone has (nray/ncone) launched rays, istart=1',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ngrill)'
      vid=ncvdef2(ncid,'ngrill',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,41,
     +     'Number of rf source ray grills, istart=2,3',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,47,
     +     ',nray/ngrill rays launched per grill istart=2,3',istatus)
      call ncaptc2(ncid,vid,'long_name2',NCCHAR,55,
     +     ',Each grill has (nray/ngrill) launched rays, istart=2,3',
     +     istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ionetwo)'
      vid=ncvdef2(ncid,'ionetwo',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +            'if ionetwo=1 then calculate CD',istatus)
      call check_err(istatus)
c--------------------------
c     Ray data
c--------------------------
c      write(*,*)'before ncvdef2(nray)'
      vid=ncvdef2(ncid,'nray',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +            'Number of rays',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(nharm)'
      vid=ncvdef2(ncid,'nharm',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'First harmonic number',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(freqncy)'
      vid=ncvdef2(ncid,'freqcy',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +            'Wave frequency',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'Hz',istatus)

c      write(*,*)'before ncvdef2(i_emission)'
      vid=ncvdef2(ncid,'i_emission',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,63,
     +'emission switch: =0 now emission, =1 use emission calculations',
     +istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(i_emission_spectrum)'
      vid=ncvdef2(ncid,'i_emission_spectrum',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,91,
     +'emission spectrum switch: =0 no emission spectrum, '//
     +'  =1 use emission spectrum  calculations',
     +istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(nrayelt)'

      vid=ncvdef2(ncid,'nrayelt',NCLONG,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     +            'Number of ray elements for each ray',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ws)'

      vid=ncvdef2(ncid,'ws',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29,
     +           'poloidal distance along a ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

c      write(*,*)'before ncvdef2(seikon)'

      vid=ncvdef2(ncid,'seikon',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,7,
     +           'eikonal',istatus)


c      write(*,*)'before ncvdef2(spsi)'

      vid=ncvdef2(ncid,'spsi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,45,
     +         'normalized small radius=rho given by indexrho',istatus)

c      write(*,*)'before ncvdef2(wr)'

      vid=ncvdef2(ncid,'wr',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)


c      write(*,*)'before ncvdef2(wphi)'

      vid=ncvdef2(ncid,'wphi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'toroidal angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'radians',istatus)

c      write(*,*)'before ncvdef2(wz)'

      vid=ncvdef2(ncid,'wz',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,15,
     +           'vertical height',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

c      write(*,*)'before ncvdef2(w_theta_pol)'

      vid=ncvdef2(ncid,'w_theta_pol',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'poloidal angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'degree',istatus)

c      write(*,*)'before ncvdef2(wnpar)'

      vid=ncvdef2(ncid,'wnpar',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'parallel refractive index',istatus)

c      write(*,*)'before ncvdef2(wnper)'

      vid=ncvdef2(ncid,'wnper',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'perpendicular refractive index',istatus)

c      write(*,*)'before ncvdef2(delpwr)'

      vid=ncvdef2(ncid,'delpwr',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'power in ray channel',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)

c      write(*,*)'before ncvdef2(sdpwr)'


      vid=ncvdef2(ncid,'sdpwr',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,48,
     +      'Ion collisionless absorption coeff (all species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

cSm050411      if (iabsorp.eq.3.and.nbulk.gt.1) then
cSm061127      if ((iabsorp.eq.3.or.iabsorp.eq.9).and.nbulk.gt.1) then
      if (((iabsorp.eq.3.or.iabsorp.eq.9).or.
     &    (iabsorp.eq.91.or.iabsorp.eq.92))
     &     .and.nbulk.gt.1) then
      vid=ncvdef2(ncid,'salphas',NCDOUBLE,3,ray_dimss,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +     'Ion collisionless absorption coeff (each species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)
      endif

c      write(*,*)'before ncvdef2(wdnpar)'

      vid=ncvdef2(ncid,'wdnpar',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,6,
     +           'wdnpar',istatus)

c      write(*,*)'before ncvdef2(cwexde)'

c     Added 3rd dimension equal to 2 accomodates complex data.
c      write(*,*)'ncid=',ncid
c      write(*,*)'ray_dims',ray_dims
      vid=ncvdef2(ncid,'cwexde',NCDOUBLE,3,ray_dims,istatus)
c      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ex/E Polarization',istatus)

c      write(*,*)'before ncvdef2(cweyde)'

      vid=ncvdef2(ncid,'cweyde',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ey/E Polarization',istatus)

c      write(*,*)'before ncvdef2(cwezde)'

      vid=ncvdef2(ncid,'cwezde',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ez/E Polarization',istatus)

c      write(*,*)'before ncvdef2(fluxn)'

      vid=ncvdef2(ncid,'fluxn',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +           'fluxn, Stix norm, |E|=1',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13,
     +           'ergs/sec/cm^2',istatus)

c      write(*,*)'before ncvdef2(sbtot)'

      vid=ncvdef2(ncid,'sbtot',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +           'Magnetic field strength',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'gauss',istatus)

c      write(*,*)'before ncvdef2(cene)'

      vid=ncvdef2(ncid,'sene',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,17,
     +           'Density along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14,
     +           'particles/cm^3',istatus)
     
c      write(*,*)'before ncvdef2(ste)'

      vid=ncvdef2(ncid,'ste',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +           'Temperature along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)

      vid=ncvdef2(ncid,'szeff',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'Zeff along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13,
     +           'dimensionless',istatus)
c      write(*,*)'before ncvdef2(salphac)'

      vid=ncvdef2(ncid,'salphac',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'Collisional damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

c      write(*,*)'before ncvdef2(salphal)'

      vid=ncvdef2(ncid,'salphal',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Linear damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

c      write(*,*)'before ncvdef2(sb_r)'

      vid=ncvdef2(ncid,'sb_r',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_r magnetic field component',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'guass',istatus)

c      write(*,*)'before ncvdef2(sb_z)'

      vid=ncvdef2(ncid,'sb_z',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_z magnetic field component',istatus) 
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'guass',istatus)

c      write(*,*)'before ncvdef2(sb_phi)'

      vid=ncvdef2(ncid,'sb_phi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'B_phi magnetic field component',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'guass',istatus)

c      write(*,*)'before ncvdef2(wn_r)'

      vid=ncvdef2(ncid,'wn_r',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_r refractive index component',istatus)

c      write(*,*)'before ncvdef2(wn_z)'

      vid=ncvdef2(ncid,'wn_z',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_z refractive index component',istatus) 

c      write(*,*)'before ncvdef2(wn_phi)'

      vid=ncvdef2(ncid,'wn_phi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32,
     +           'N_phi refractive index component',istatus)
     
c      write(*,*)'before ncvdef2(wgr_r)'

      vid=ncvdef2(ncid,'vgr_r',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_r normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(wgr_z)'

      vid=ncvdef2(ncid,'vgr_z',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_z normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(wgr_hpi)'

      vid=ncvdef2(ncid,'vgr_phi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +           'vgroup_phi normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(flux_z)'

       
      vid=ncvdef2(ncid,'flux_z',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +           'flux_z normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(flux_r)'

       
      vid=ncvdef2(ncid,'flux_r',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +           'flux_r normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(flux_phi)'

       
      vid=ncvdef2(ncid,'flux_phi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'flux_phi normalized to c',istatus)

c--------------------------------------------------------
c      CD efficiency along rays
c--------------------------------------------------------
      if (ionetwo.eq.1)then
c         write(*,*)'before ncvdef2(w_eff_nc)'
         vid=ncvdef2(ncid,'w_eff_nc',NCDOUBLE,2,ray_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'CD efficiency along a ray',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,27,
     +           '(A/cm**2)/(erg/(sec*cm**3))',istatus)
       endif
c--------------------------------------------------------
c      DC electric field conductivity from adj 
c--------------------------------------------------------
       if (i_adj.eq.1)then

c         write(*,*)'before ncvdef2(i_adj)'
         vid=ncvdef2(ncid,'i_adj',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     &   '=0 no adj; =1 use adj calculations ',istatus)
         call check_err(istatus)

c         write(*,*)'before ncvdef2 npsi0)'
         vid=ncvdef2(ncid,'npsi0',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,53,
     &   'Number of flux surfaces where ADJ equation was solved',
     &   istatus)
         call check_err(istatus)


         nrho_adjid=ncddef2(ncid,'nrho_adj',npsi0,istatus)
         call check_err(istatus)

c         write(*,*)'sigma_E_adj',sigma_E_adj
c         write(*,*)'rho_adj',rho_adj
c         write(*,*)'before ncvdef2(sigma_E_adj)'
         vid=ncvdef2(ncid,'sigma_E_adj',NCDOUBLE,1,nrho_adjid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,74,
     &   'DC electric field conductivity at radii where '//
     &   'adj function is calculated ',istatus) 
         call ncaptc2(ncid,vid,'units',NCCHAR,19,
     +           'sigma/sigma_Spitzer',istatus)

c         write(*,*)'before ncvdef2(rho_adj)'
         vid=ncvdef2(ncid,'rho_adj',NCDOUBLE,1,nrho_adjid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,50,
     &   'small radii mesh where adj function is calculated ',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)

       endif
c--------------------------
c     determine OX conversion data for i_ox=2 case
c--------------------------
      if (i_ox.eq.2) then

         vid=ncvdef2(ncid,'i_ox_conversion',NCLONG,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     &            'Option equals 1 after OX conversion',istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'transm_ox',NCDOUBLE,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     &   'OX transmission coefficient Preinhaelter and Kopecky 1973',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cn_par_optimal',NCDOUBLE,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     &              'optimal N parallel for OX transmission',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cnpar_ox',NCDOUBLE,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,33,
     &              'N parallel before OX transmission',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cn_b_gradpsi',NCDOUBLE,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,45,
     &              'N along [B^*grad(psi)] before OX transmission',
     &              istatus)
         call check_err(istatus)
         
      endif !i_ox=2

c--------------------------
c     Plasma data
c--------------------------
     
c      write(*,*)'before ncvdef2(eqdskin)'

       
      vid=ncvdef2(ncid,'eqdskin',NCCHAR,1,char256id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39,
     +          'Name of input eqdsk, for eqsource=eqdsk',istatus)
     
c      write(*,*)'before ncvdef2(nbulk)'

       
      vid=ncvdef2(ncid,'nbulk',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,43,
     +          'Number of Maxl plasma cmpts, electrons+ions',istatus)
      call check_err(istatus)     

c      write(*,*)'before ncvdef2(dmas)'

       

      vid=ncvdef2(ncid,'dmas',NCDOUBLE,1,nbulkid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +           'plasma species mass: electrons, then ions',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,27,
     +           'Normalized to electron mass',istatus)
     
c      write(*,*)'before ncvdef2(charge)'


      vid=ncvdef2(ncid,'charge',NCDOUBLE,1,nbulkid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44,
     +           'plasma species charge: electrons, then ions',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,31,
     +           'Normalized to electronic charge',istatus)

c------------------------------------------------------------------
c     for total power [erg/sec] absorbed at all reflections at all rays
c------------------------------------------------------------------
c      write(*,*)'before ncvdef2(w_tot_pow_absorb_at_refl_nc)'
      vid=ncvdef2(ncid,'w_tot_pow_absorb_at_refl_nc',NCDOUBLE,
     &0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,47,
     +    'Total power absorbed at reflections of all rays',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
c.................................
cl    1.1.4 end the define-mode and start the data-mode
c     p. 51-2 of manual
c      write(*,*)'ncid',ncid
      call ncendf2(ncid,istatus)
      call check_err(istatus)

c      write(*,*)'end initialization'

      endif               ! End initialize




c.......................................................................
cl    1. Writing data
c

      if ( kopt.eq.0 ) then
c      write(*,*)'wrtnetcdf,  kopt=',kopt
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c--------------------------

      call ncvid2(vid,ncid,'version',istatus)
      ll=length_char(version)
      call ncvptc2(ncid,vid,1,ll,version,ll,istatus)
c      write(*,*)'wrtnetcdf after ncptc version istatus',istatus

      call ncvid2(vid,ncid,'mnemonic',istatus)
      ll=length_char(trim(mnemonic))
      call ncvptc2(ncid,vid,1,ll,trim(mnemonic),ll,istatus)
c      write(*,*)'wrtnetcdf after ncptc mnemonic istatus',istatus

      call ncvid2(vid,ncid,'eqdskin',istatus)
      ll=length_char(trim(eqdskin))
      call ncvptc2(ncid,vid,1,ll,trim(eqdskin),ll,istatus)
c      write(*,*)'wrtnetcdf after ncptc eqdskin istatus',istatus

      call ncvid2(vid,ncid,'dmas',istatus)
      call ncvpt_doubl2(ncid,vid,1,nbulk,dmas,istatus)
c      write(*,*)'wrtnetcdf after ncptc dmas istatus',istatus

      call ncvid2(vid,ncid,'charge',istatus)
      call ncvpt_doubl2(ncid,vid,1,nbulk,charge,istatus)
c      write(*,*)'wrtnetcdf after ncptc charge istatus',istatus

      call ncvid2(vid,ncid,'istart',istatus)
      call ncvpt_int2(ncid,vid,1,1,istart,istatus)
c      write(*,*)'wrtnetcdf after ncptc istart istatus',istatus
  
      call ncvid2(vid,ncid,'ncone',istatus)
      call ncvpt_int2(ncid,vid,1,1,ncone,istatus)
c      write(*,*)'wrtnetcdf after ncptc nccone istatus',istatus

      call ncvid2(vid,ncid,'ngrill',istatus)
      call ncvpt_int2(ncid,vid,1,1,ngrill,istatus)
c      write(*,*)'wrtnetcdf after ncptc ngrill istatus',istatus

      call ncvid2(vid,ncid,'ionetwo',istatus)
      call ncvpt_int2(ncid,vid,1,1,ionetwo,istatus)
c      write(*,*)'wrtnetcdf after ncptc ionetwo istatus',istatus

c--------------------------
c     Run specs
c--------------------------
      call ncvid2(vid,ncid,'id',istatus)
      call ncvpt_int2(ncid,vid,1,1,id,istatus)

c      write(*,*)'netcdfr3d.f refl_loss',refl_loss

      call ncvid2(vid,ncid,'refl_loss',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,refl_loss,istatus)
c      write(*,*)'wrtnetcdf after ncvpt refl_loss istatus',istatus

      call ncvid2(vid,ncid,'iabsorp',istatus)
      call ncvpt_int2(ncid,vid,1,1,iabsorp,istatus)
c      write(*,*)'wrtnetcdf after ncptc iabsorp istatus',istatus


      call ncvid2(vid,ncid,'ieffic',istatus)
      call ncvpt_int2(ncid,vid,1,1,ieffic,istatus)
c      write(*,*)'wrtnetcdf after ncptc ieffic istatus',istatus


cSAP080303
      call ncvid2(vid,ncid,'ion_absorption',istatus)
      ll=length_char(ion_absorption)
      call ncvptc2(ncid,vid,1,ll,ion_absorption,ll,istatus)
c      write(*,*)'wrtnetcdf after ncptc ion_absorption istatus',istatus

      call ncvid2(vid,ncid,'iflux',istatus)
      call ncvpt_int2(ncid,vid,1,1,iflux,istatus)
c      write(*,*)'wrtnetcdf after ncpvt iflux istatus',istatus

      call ncvid2(vid,ncid,'ioxm',istatus)
      call ncvpt_int2(ncid,vid,1,1,ioxm,istatus)
c      write(*,*)'wrtnetcdf after ncvpt ioxm istatus',istatus

      call ncvid2(vid,ncid,'ioxm_n_npar',istatus)
      call ncvpt_int2(ncid,vid,1,1,ioxm_n_npar,istatus)
c      write(*,*)'wrtnetcdf after ncvpt ioxm_n_npar istatus',istatus


      call ncvid2(vid,ncid,'jwave',istatus)
      call ncvpt_int2(ncid,vid,1,1,jwave,istatus)
c      write(*,*)'wrtnetcdf after ncvpt jwave istatus',istatus

      call ncvid2(vid,ncid,'i_geom_optic',istatus)
      call ncvpt_int2(ncid,vid,1,1,i_geom_optic,istatus)
c      write(*,*)'wrtnetcdf after ncvpt i_geom_optic istatus',istatus

      
      if (iabsorp.eq.4) then
         call ncvid2(vid,ncid,'i_im_nperp',istatus)
         call ncvpt_int2(ncid,vid,1,1,i_im_nperp,istatus)
      endif

      call ncvid2(vid,ncid,'nbulk',istatus)
      call ncvpt_int2(ncid,vid,1,1,nbulk,istatus)
c      write(*,*)'wrtnetcdf after ncvpt nbulk istatus',istatus

c--------------------------
c     Ray data
c--------------------------
      call ncvid2(vid,ncid,'nray',istatus)
      call ncvpt_int2(ncid,vid,1,1,nrayl,istatus)
c      write(*,*)'wrtnetcdf after ncvpt nray istatus',istatus


      call ncvid2(vid,ncid,'nharm',istatus)
      call ncvpt_int2(ncid,vid,1,1,nharm,istatus)
c      write(*,*)'wrtnetcdf after ncvpt nharm istatus',istatus

      call ncvid2(vid,ncid,'freqcy',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,freqcy,istatus)
c      write(*,*)'wrtnetcdf after ncvpt freqcy istatus',istatus

      call ncvid2(vid,ncid,'i_emission',istatus)
      call ncvpt_int2(ncid,vid,1,1,i_emission,istatus)      
c      write(*,*)'wrtnetcdf after ncvpt i_emission istatus',istatus

      
      call ncvid2(vid,ncid,'i_emission_spectrum',istatus)
      call ncvpt_int2(ncid,vid,1,1,i_emission_spectrum,istatus)
c      write(*,*)'wrtnetcdf after ncvpt i_emission_spectrum istatus',
c     &istatus

      call ncvid2(vid,ncid,'nrayelt',istatus)
      call ncvpt_int2(ncid,vid,1,nrayl,nrayelt_nc,istatus)
c      write(*,*)'wrtnetcdf after ncvpt nrayelt istatus',istatus

c      write(*,*)'in wrtnetcdf nrayl,neltmax',nrayl,neltmax
c      do i=1,nrayl
c         write(*,*)'number of rayi,nrayelt_nc(i)',i,nrayelt_nc(i)
c         do j=1,nrayelt_nc(i)
c            write(*,*)'j,ws_nc(j,i)',j,ws_nc(j,i)
c         enddo
c      enddo

cSAP090903
c      call pack21(ws_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
c      write(*,*)'wrtnetcdf before  pack21 ws_ns'
c      write(*,*)'nrelta,nrayl,neltmax',
c     &           nrelta,nrayl,neltmax

      call pack21(ws_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
c      write(*,*)'wrtnetcdf after  pack21 ws_ns'
c      do i=1,nrayl
c         write(*,*)'wrtnetcdf tem1 i=',i
c         do j=1,neltmax          
c            k=(i-1)*neltmax+j
c            write(*,*)'j,i,k tem1(k)',j,i,k,tem1(k)
c         enddo
c      enddo
c      write(*,*)'ws: vid,start,ray_count',vid,start,ray_count
      call ncvid2(vid,ncid,'ws',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(seikon_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(seikon_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'seikon',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(spsi_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(spsi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'spsi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(wr_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(wr_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(wphi_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(wphi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wphi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(wz_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(wz_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wz',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(w_theta_pol_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(w_theta_pol_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'w_theta_pol',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(wnpar_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(wnpar_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wnpar',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(wnper_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(wnper_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wnper',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(delpwr_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(delpwr_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'delpwr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(sdpwr_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(sdpwr_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sdpwr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSm050411      if (iabsorp.eq.3) then
cSm061127      if (iabsorp.eq.3.or.iabsorp.eq.9) then
      if ((iabsorp.eq.3.or.iabsorp.eq.9).or.
     &    (iabsorp.eq.91.or.iabsorp.eq.92)) then
      do i=2,nbulk
cSAP090903
c      call pack21(salphas_nc(1,1,i),1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(salphas_nc(1,1,i),1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      starts(1)=start(1)
      starts(2)=start(2)
      starts(3)=i-1  ! nbulk-1 ion species salphas are put in .nc file
      call ncvid2(vid,ncid,'salphas',istatus)
      call ncvpt_doubl2(ncid,vid,starts,ray_counts,tem1,istatus)
      enddo
      endif

cSAP090903
c      call pack21(wdnpar_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(wdnpar_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wdnpar',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      ii=0
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=0.5d0*(cwexde_nc(i,j)+dconjg(cwexde_nc(i,j)))
         enddo
      enddo
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=-cei*0.5d0*(cwexde_nc(i,j)-dconjg(cwexde_nc(i,j)))
         enddo
      enddo
      call ncvid2(vid,ncid,'cwexde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      ii=0
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=0.5d0*(cweyde_nc(i,j)+dconjg(cweyde_nc(i,j)))
         enddo
      enddo
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=-cei*0.5d0*(cweyde_nc(i,j)-dconjg(cweyde_nc(i,j)))
         enddo
      enddo
      call ncvid2(vid,ncid,'cweyde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      ii=0
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=0.5d0*(cwezde_nc(i,j)+dconjg(cwezde_nc(i,j)))
         enddo
      enddo
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=-cei*0.5d0*(cwezde_nc(i,j)-dconjg(cwezde_nc(i,j)))
         enddo
      enddo
      call ncvid2(vid,ncid,'cwezde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(fluxn_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(fluxn_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'fluxn',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(sbtot_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(sbtot_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sbtot',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(sene_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(sene_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sene',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(ste_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(ste_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'ste',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cBH111028
      call pack21(szeff_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'szeff',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(salphac_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(salphac_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'salphac',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(salphal_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(salphal_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'salphal',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(sb_r_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(sb_r_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sb_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(sb_z_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(sb_z_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sb_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(sb_phi_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(sb_phi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sb_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(wn_r_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(wn_r_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wn_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(wn_z_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(wn_z_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wn_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(wn_phi_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(wn_phi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wn_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(vgr_r_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(vgr_r_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'vgr_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(vgr_z_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(vgr_z_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'vgr_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(vgr_phi_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(vgr_phi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'vgr_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(flux_r_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(flux_r_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'flux_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(flux_z_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
      call pack21(flux_z_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'flux_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call pack21(flux_phi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call pack21(flux_phi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'flux_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c--------------------------------------------------------------
c      CD efficiency along rays
c--------------------------------------------------------
c      write(*,*)'netcdfr3d ionetwo',ionetwo
       if (ionetwo.eq.1) then            
cSAP090903
c         call pack21(w_eff_nc,1,nrelta,1,nraya,tem1,neltmax,nrayl)
         call pack21(w_eff_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
          call ncvid2(vid,ncid,'w_eff_nc',istatus)
          call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
       endif
c--------------------------------------------------------
c      DC electric field conductivity from adj 
c---------------------------`-----------------------------
c      write(*,*)'netcdfr3d i_adj',i_adj
       if (i_adj.eq.1)then  
c           write(*,*)'before vid=ncvid(i_adj'
           call ncvid2(vid,ncid,'i_adj',istatus)
           call ncvpt_int2(ncid,vid,1,1,i_adj,istatus)

c           write(*,*)'before vid=ncvid(npsi0i'
           call ncvid2(vid,ncid,'npsi0',istatus)
           call ncvpt_int2(ncid,vid,1,1,npsi0,istatus)

c            write(*,*)'before vid=ncvid(sigma_E_adj'

           call ncvid2(vid,ncid,'sigma_E_adj',istatus)
           call ncvpt_doubl2(ncid,vid,start,npsi0,sigma_E_adj,istatus)

c           write(*,*)'before vid=ncvid(rho_adj'

           call ncvid2(vid,ncid,'rho_adj',istatus)
           call ncvpt_doubl2(ncid,vid,start,npsi0,rho_adj,istatus)
c           write(*,*)'after vid=ncvpt(rho_adj'

       endif
c--------------------------
c     write OX conversion data for i_ox=2 case
c--------------------------
      if (i_ox.eq.2) then

         call ncvid2(vid,ncid,'i_ox_conversion',istatus)
         call ncvpt_int2(ncid,vid,1,nrayl,i_ox_conversion_nc,istatus)

         call ncvid2(vid,ncid,'transm_ox',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,transm_ox_nc,istatus)

         call ncvid2(vid,ncid,'cn_par_optimal',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_par_optimal_nc,istatus)

c         call ncvid2(vid,ncid,'cn_par_optimal',istatus)
c         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_par_optimal_nc,istatus)

         call ncvid2(vid,ncid,'cnpar_ox',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cnpar_ox_nc,istatus)

         call ncvid2(vid,ncid,'cn_b_gradpsi',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_b_gradpsi_nc,istatus)

      endif !i_ox=2 

c------------------------------------------------------------------
c     for total power absorbed at all reflections at all rays
c------------------------------------------------------------------
      call ncvid2(vid,ncid,'w_tot_pow_absorb_at_refl_nc',istatus)
      call ncvpt_doubl2(ncid,vid,1,1, w_tot_pow_absorb_at_refl_nc,
     +                  istatus)
      
C-----------------------------------------------------------------------
c
cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

      endif                    ! End write netcdf data

cSAP090903
      deallocate (tem1,STAT=istat)

      return
      end
c
c


      subroutine wrtnetcdf_prof(netcdfnml,kopt)
      implicit none
      !implicit integer (i-n), real*8 (a-h,o-z)

c     If the namelist ionetwo.eq.1 it will write the current
c     and power profiles into existing netcdf file: netcdfnml 

c-----profiles description:

c     normalized small radius bin boundaries rho_bin(i), i=1,NR
c     for small radius bin(i)=[rho_bin(i),rho+bin(i+1)], i=1,NR-1
c     rho_bin(1)=0, rho_bin(NR)=1

c     the centers of small radius bin(i) rho_bin_center(i), i=1,NR-1
c     rho_bin_center(i)=0.5*[rho_bin(i)+rho_bin(i+1)],
c     rho_bin_center(1)=0.5*rho_bin(2)
c     rho_bin_center(NR-1)=0.5*(rho_bin(NR-1)+1)

c     spower(i)   i=1,NR-1, power [erg/sec] in small radius bin(i)
c     powden(i)  i=1,NR-1, power density [erg/(cm**3*sec)]
c                           in small radius bin(i)
c
c     s_cur_den_parallel(i) i=1,NR-1 parallel averaged current density
c                           [A/cm**2] 
c               <j_parallel>=Integral{dl_poloidal*j_parallel/B_poloidal}
c     s_cur_den_onetwo(i) i=1,NR-1 ONETWO current density 
c                        <j^.B^>/B_0=<>j_parallel*modB>/B_0=
c                        =<j_parallel><B**2>/(B_0*<B>)
c     s_cur_den_toroidal(i) i=1,NR-1 toroidal current density
c                         j_phi=<j_parallel>f<1/R**2>/(<B><1/R>)
c     s_cur_den_poloidal(i) i=1,NR-1   poloidal current density
c                         i_poloidal=<j_parallel>B_poloidal/<B>
c                         B_poloidal is taken at theta_poloidal=0
c

c-----total: absorbed power and toroidal current

c     total power = power_total [erg/sec]
c     total toroidal current = tor_curr_total [A]

c     do i=1,NR-1
c        power_total=power_total+spower(i)
c        tor_curr_total=tor_curr_total+scurrent(i)
c     enndo

      include 'param.i'     
      include 'one.i'
      include 'onetwo.i'  !For profiles and total current.
      include 'rho.i'     !For areatot,voltot,torftot,totlength
cSm070201
      include 'writencdf.i'     !
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'
c-----input
      integer kopt !kopt.eq.1 initialize writing
                   !kopt.eq.0 write data
      character(*) netcdfnml ! input filename
c-----locals
      integer kk

      real*8 hrho,power_total,tor_curr_total,GA_tor_cur_total

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,error_code,ncvdef2,ncdid2,ncddef2
      integer nrhoid,nrhomid,nbulkid,nbulkmid  !ion species dimension
     &,nbulk_ion_id 
      integer nrho,i,rdims(2),rdimss(2)
      integer start(2),starts(2),count(2),counts(2)

      
c      real*8 rho_bin(NR) ! normalized small radius bin(i) boundaries
c      real*8 rho_bin_center(NR-1) ! centers of small radius bin(i) 

      real*8 tem1(NRA*nbulka)

c-----externals
      real*8 zeffrho,temperho,densrho

      data start/1,1/, starts/1,1/

      save

c      write(*,*)'netcdfr3d.f in wrtnetcdf_prof start=',start

c      write(*,*)'wrtnetcdf_prof,  ionetwo=',ionetwo

c      write(*,*)'wrtnetcdf_prof: powtot_e,powtot_i',powtot_e,powtot_i

      if (ionetwo.ne.1) return    !nothing to do
      
c.......................................................................
cl    1. Initialize part, creating new netcdf file
c

c --- begin if ---
      if ( kopt.eq.1 ) then
c      write(*,*)'wrtnetcdf_prof,  kopt=',kopt

C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c

c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus)    ! Open existing netCDF file
      call check_err(istatus)
      
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf2(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf2(integer ncid,integer error_code)

      call ncredf2(ncid,error_code)
      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 define dimensions
c     integer function ncddef(ncid,character(*) dim_name,
c                             integer cdim_siz, integer error_code)
c     returns dimension id.
c     p. 52 of netcdf manual

      
c-----For radial profiles data:

c-----define dimensions      
      nrhoid=ncddef2(ncid,'nrho',NR,istatus)
      call check_err(istatus)

      nrhomid=ncddef2(ncid,'nrhom',NR-1,istatus)
      call check_err(istatus)

      if (nbulk.gt.1) then
         nbulkmid=ncdid2(ncid,'nbulkm',istatus)
         call check_err(istatus)
         rdims(1)=nrhomid
         rdims(2)=nbulkmid
      endif

      nbulkid=ncdid2(ncid,'nbulk',istatus)
      call check_err(istatus)

      count(1)=NR-1
      count(2)=nbulk-1
      rdimss(1)=nrhoid
      rdimss(2)=nbulkid
      counts(1)=NR
      counts(2)=nbulk
     
c    
c     For toroidal current and power:
c-----define variables and attributes
      vid=ncvdef2(ncid,'parallel_cur_total',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Total parallel current',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'A',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'toroidal_cur_total',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Total toroidal current',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'A',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'poloidal_cur_total',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Total poloidal current',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'A',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'power_inj_total',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +            'Total injected power',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'power_total',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +            'Total absorbed power',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'powtot_e',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +            'Total power to electrons',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'powtot_i',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     +            'Total power to ions',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'powtot_cl',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +            'Collisional power absorbed',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)

      if (iabsorp.eq.3. and. nbulk.gt.1) then
         vid=ncvdef2(ncid,'powtot_s',NCDOUBLE,1,nbulkmid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +        'Power to individual ion species',istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +        'erg/sec',istatus)
         call check_err(istatus)
      endif         

c.......................................................................
cl    1.1.4 define variables

c     Note, the variable IDs (denoted below as vid) are
c     not saved here in this subroutine; rather, the IDs
c     are retrieved from the netCDF data file, as needed,
c     by calling the netCDF routine ncvid.

c     netCDF variable_define:
c     integer function ncvdef(ncid,variable_name,variable_type,
c                number_of_dimensions,
c                vector_for_length_per_dimension,error_code)
c       returns varid.
c     Note: Unlimited dimension must be last
c     Refer to p. 77, netcdf-2 manual

c     NCDOUBLE for REAL*8.  This is independent of the
c       internal representation of data on a particular
c       machine, e.g., 32- or 64-bit based calculations.


c--------------------------
c     current and power profiles data
c--------------------------
      vid=ncvdef2(ncid,'NR',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29,
     +            'Number of small radius points',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'voltot',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40,
     +            'Total volume in last closed flux surface',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm^3',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'areatot',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Total cross-sectional area of LCFS',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm^2',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'pollentot',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +            'Poloidal length of LCFS',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'cm',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'torftot',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +            'Toroidal flux through LCFS',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,9,
     +           'Tesla*m^2',istatus)
      call check_err(istatus)
 
      vid=ncvdef2(ncid,'indexrho',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40,
     +           'Radial coord type: 2 gives sqrt(tor flx)',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'psifactr',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +            'Reduces Psi-Value of LCFS',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'Should be .le.1',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'binvol',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Volumes of radial bins',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm^3',istatus)
      call check_err(istatus)
 
      vid=ncvdef2(ncid,'binarea',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +            'Areas of radial bins',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm^2',istatus)
      call check_err(istatus)
 
      vid=ncvdef2(ncid,'pollen',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31,
     +            'Poloidal lengths of radial bins',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'cm',istatus)
      call check_err(istatus)
 
      vid=ncvdef2(ncid,'rho_bin',NCDOUBLE,1,nrhoid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     +            'normalized small radius bin boundaries',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'rho_bin_center',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     +            'normalized small radius bin centers',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'densprof',NCDOUBLE,2,rdimss,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44,
     +     'plasma density at bin boundaries, e and ions',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14,
     +     'particles/cm^3',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'temprof',NCDOUBLE,2,rdimss,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +     'plasma temperatures at bin boundaries, e and ions',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +     'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'zefprof',NCDOUBLE,1,rdimss(1),istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +     'plasma Zeff at bin boundaries, e and ions',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8,
     +     'unitless',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'spower',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'power per bin profile',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)

c----------------------------------------------------------------
c     plasma profiles VS major radius
c---------------------------------------------------------------
      vid=ncvdef2(ncid,'w_r_densprof_nc',NCDOUBLE,1,nrhoid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39,
     +            'major radius r mesh for plasma profiles',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +     'm ',istatus)
      call check_err(istatus)
cSAP090731
c      vid=ncvdef2(ncid,'w_dens_vs_r_nc',NCDOUBLE,1,nrhoid,istatus)
      vid=ncvdef2(ncid,'w_dens_vs_r_nc',NCDOUBLE,2,rdimss,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     +     'plasma density at bin bndries, e and ions' //
     +' vs major radius',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13,
cSAP090731
c     +     'particles/cm^3',istatus)
     +     'particles/m^3',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_temp_vs_r_nc',NCDOUBLE,2,rdimss,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,62,
     +     'plasma temperatures at bin bndries, e and ions' //
     +' vs major radius',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +     'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_zeff_vs_r_nc',NCDOUBLE,1,nrhoid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,42,
     +     'plasma Zeff at bin bndries vs major radius ',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8,
     +     'unitless',istatus)
      call check_err(istatus)


c----------------------------------------------------------------
c     averaged current densities according to GA-memo
c----------------------------------------------------------------
      vid=ncvdef2(ncid,'GA_tor_cur_total',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +            'GA memo total toroidal current',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'A',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'s_cur_den_parallel',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     &            'averaged parallel current density profile',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     &           'A/cm**2',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'s_cur_den_onetwo',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40,
     &            'ONETWO current <j.B>/B_0 density profile',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     &           'A/cm**2',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'s_cur_den_toroidal',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,60,
     &'toroidal current <j_par>f<1/r**2>/(<B><1/r>) density profile',
     &istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'A/cm**2',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'s_cur_den_poloidal',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,62,
     &'poloidal current <j_par>B_pol(theta_pol=0)/<B> density profile',
     &istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'A/cm**2',istatus)
      call check_err(istatus)

c-------------------------------------------------------

      vid=ncvdef2(ncid,'powden',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'power density profile',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'erg/(cm**3*sec)',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'powden_e',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,48,
     +     'power density profile to electrons, bin centered',istatus)  
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,26,
     +           'From wave-particle damping',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'erg/(cm**3*sec)',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'powden_cl',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,48,
     +     'power density profile to electrons, bin centered',istatus)  
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,24,
     +           'From collisional damping',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'erg/(cm**3*sec)',istatus)
      call check_err(istatus)
      vid=ncvdef2(ncid,'powden_i',NCDOUBLE,1,nrhomid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,43,
     +     'power density profile to ions, bin centered',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,26,
     +           'From wave-particle damping',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'erg/(cm**3*sec)',istatus)
      call check_err(istatus)

      if (iabsorp.eq.3) then
         vid=ncvdef2(ncid,'powden_s',NCDOUBLE,2,rdims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     +      'power density profile to individual ion species, bin cent',
     +      istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +        'erg/(cm**3*sec)',istatus)
         call check_err(istatus)
      endif

    
c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf2(ncid,istatus)
      call check_err(istatus)

      endif               ! End initialize, kopt=1


c.......................................................................
cl    1. Writing data
c
      

      if ( kopt.eq.0 ) then
c      write(*,*)'wrtnetcdf_prof,  kopt,NR=',kopt,NR
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c--------------------------
      nrho=NR
c-----create:
c     small radius array rho(i) for bin boundaries
      hrho=1.d0/(NR-1)
      do i=1,NR
         rho_bin(i)=hrho*(i-1)
      enddo

c-----create:
c     small radius array rho_bin_center(i) for bin(i) centers
      do i=1,NR-1
         rho_bin_center(i)=0.5d0*(rho_bin(i)+rho_bin(i+1))
      enddo

c-----evaluate radial profiles of density, temperature, zeff
c-----to put with radial profile data.
c-----(densrho is in units 10**13/cm**3, temprof in keV.)
      do kk=1,nbulk
         do i=1,NR
            densprof(i,kk)=densrho(rho_bin(i),kk)*1.e13
            temprof(i,kk)=temperho(rho_bin(i),kk)
         enddo
      enddo
      do i=1,NR
         zefprof(i)=zeffrho(rho_bin(i))
      enddo
      

c-----calculate total toroidal current and total absorbed power
      power_total=0.d0
      tor_curr_total=0.d0
    
      do i=1,NR-1
         power_total=power_total+spower(i)              
      enddo


c--------------------------
c     write
c     the number of radial points NR for bin boundaries  and 
c     total toroidal current and total power data
c--------------------------     
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.1'
      call ncvid2(vid,ncid,'NR',istatus)
      call ncvpt_int2(ncid,vid,1,1,NR,istatus)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.2'
      call ncvid2(vid,ncid,'voltot',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,1.d6*voltot,istatus)

      call ncvid2(vid,ncid,'areatot',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,1.d4*areatot,istatus)

      call ncvid2(vid,ncid,'pollentot',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,1.d2*totlength,istatus)

      call ncvid2(vid,ncid,'torftot',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,torftot,istatus)

      call ncvid2(vid,ncid,'indexrho',istatus)
      call ncvpt_int2(ncid,vid,1,1,indexrho,istatus)

      call ncvid2(vid,ncid,'psifactr',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,psifactr,istatus)

      call ncvid2(vid,ncid,'binvol',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,binvol,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'binarea',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,binarea,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'pollen',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,pollen,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'parallel_cur_total',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,parallel_cur_total,istatus)

      call ncvid2(vid,ncid,'toroidal_cur_total',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,toroidal_cur_total,istatus)
   
      call ncvid2(vid,ncid,'poloidal_cur_total',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,poloidal_cur_total,istatus)

      call ncvid2(vid,ncid,'power_inj_total',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,powtott,istatus)

      call ncvid2(vid,ncid,'power_total',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,power_total,istatus)

      call ncvid2(vid,ncid,'powtot_e',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,powtot_e,istatus)

      call ncvid2(vid,ncid,'powtot_i',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,powtot_i,istatus)

      call ncvid2(vid,ncid,'powtot_cl',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,powtot_cl,istatus)

      if (iabsorp.eq.3) then
         call ncvid2(vid,ncid,'powtot_s',istatus)
         call ncvpt_doubl2(ncid,vid,1,nbulk-1,powtot_s(2),istatus)
         call check_err(istatus)
      endif
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.3'



c --------------------------
c     current and power profiles data
c--------------------------  
      
      call ncvid2(vid,ncid,'rho_bin',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR,rho_bin,istatus)
      call check_err(istatus)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.3.1'

      call ncvid2(vid,ncid,'rho_bin_center',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,rho_bin_center,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'densprof',istatus)
c      call pack21(densprof,1,NRA,1,nbulka,tem1,NR,nbulk)
      call pack21(densprof,1,NR,1,nbulk,tem1,NR,nbulk)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.3.3'
      call ncvpt_doubl2(ncid,vid,starts,counts,tem1,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'temprof',istatus)
c      call pack21(temprof,1,NRA,1,nbulka,tem1,NR,nbulk)
      call pack21(temprof,1,NR,1,nbulk,tem1,NR,nbulk)
      call ncvpt_doubl2(ncid,vid,starts,counts,tem1,istatus)
      call check_err(istatus)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.3.4'

      call ncvid2(vid,ncid,'zefprof',istatus)
      call ncvpt_doubl2(ncid,vid,starts,counts,zefprof,istatus)
      call check_err(istatus)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.3.5'

      call ncvid2(vid,ncid,'spower',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,spower,istatus)
      call check_err(istatus)
c----------------------------------------------------------------
c     plasma profiles VS major radius
c---------------------------------------------------------------
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.4'


c------------------------------------------------------
c     Creates 1D arrays of density,temperature,zeff
c     versus major radius r at given z
c-------------------------------------------------------
      call plasma_profiles_vs_r(0.d0)
c------------------------------------------------------
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.5'
     
      call ncvid2(vid,ncid,'w_r_densprof_nc',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR,w_r_densprof_nc,istatus)
      call check_err(istatus)
 
      call ncvid2(vid,ncid,'w_dens_vs_r_nc',istatus)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.5.1'
c      call pack21(w_dens_vs_r_nc,1,NRA,1,nbulka,tem1,NR,nbulk)
      call pack21(w_dens_vs_r_nc,1,NR,1,nbulk,tem1,NR,nbulk)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.5.2'
cSAP090731 transform density units from  1/cm**3 to 1/m**3
      !tem1=tem1*1d6 !YuP[07-2017] Moved all conversion factors to
      ! subroutine plasma_profiles_vs_r(z)
      ! where w_dens_vs_r_nc() is defined.
      call ncvpt_doubl2(ncid,vid,starts,counts,tem1,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_temp_vs_r_nc',istatus)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.5.3'
c      call pack21(w_temp_vs_r_nc,1,NRA,1,nbulka,tem1,NR,nbulk)
      call pack21(w_temp_vs_r_nc,1,NR,1,nbulk,tem1,NR,nbulk)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.5.4'
      call ncvpt_doubl2(ncid,vid,starts,counts,tem1,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_zeff_vs_r_nc',istatus)
      call ncvpt_doubl2(ncid,vid,1,NR,w_zeff_vs_r_nc,istatus)
      call check_err(istatus)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.6'
c----------------------------------------------------------------
c     averaged current densities according to GA-memo
c----------------------------------------------------------------
c-----calculate total GA memo toroidal current 
      GA_tor_cur_total=0.d0
c      write(*,*)'in wrtnetcdf_prof'
      do i=1,NR-1
         GA_tor_cur_total=GA_tor_cur_total+s_cur_den_toroidal(i)*
     &                     binarea(i)
c         GA_tor_cur_total=GA_tor_cur_total+s_cur_den_onetwo(i)*
c     &                     binarea(i)

      enddo

c      write(*,*)'GA_tor_cur_total ',GA_tor_cur_total
cSAP090306
      write(*,1002)GA_tor_cur_total     
 1002 format('total toroidal current GA_tor_cur_total (A)=',
     &   1pe14.6)

      call ncvid2(vid,ncid,'GA_tor_cur_total',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,GA_tor_cur_total,istatus)

      call ncvid2(vid,ncid,'s_cur_den_parallel',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,s_cur_den_parallel,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'s_cur_den_onetwo',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,s_cur_den_onetwo,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'s_cur_den_toroidal',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,s_cur_den_toroidal,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'s_cur_den_poloidal',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,s_cur_den_poloidal,istatus)
      call check_err(istatus)


c----------------------------------------------------------------
      call ncvid2(vid,ncid,'powden',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,powden,istatus)
      call check_err(istatus)
 
      call ncvid2(vid,ncid,'powden_e',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,powden_e,istatus)
      call check_err(istatus)
 
      call ncvid2(vid,ncid,'powden_cl',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,powden_cl,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'powden_i',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,powden_i,istatus)
      call check_err(istatus)
     
      if (iabsorp.eq.3 .and. nbulk.gt.1) then
c      call pack21(powden_s(1,2),1,NRA,1,nbulka-1,tem1,NR-1,nbulk-1)
      write(*,*)'before write powden_s 2'
      call pack21(powden_s(1,2),1,NR,1,nbulk,tem1,NR,nbulk-1)
       write(*,*)'before write powden_s 2'
      call ncvid2(vid,ncid,'powden_s',istatus)
      call ncvpt_doubl2(ncid,vid,start,count,tem1,istatus)
      call check_err(istatus)
c      write(*,*)
c      do kk=2,nbulk
c         do i=1,NR-1
c            write(*,*)' i,kk,powden_s(i,kk) :',i,kk,powden_s(i,kk)
c         enddo
c      enddo
    
      endif
 

c      do i=1,NR-1
c         write(*,*)'i,scurrent(i),currden(i)',i,scurrent(i),currden(i)
c         write(*,*)'i,spower(i),powden(i)',i,spower(i),powden(i)
c      enddo     
      write(*,1000)'i rho_bin_center scurrent    currden     spower    '
     1          //'powden'
      do i=1,NR-1
         write(*,1001)i,rho_bin_center(i),scurrent(i),currden(i),
     1        spower(i),powden(i)
      enddo
 1000    format(/,1x,a)
 1001    format(i3,5(1pe12.4))



C-----------------------------------------------------------------------
c
cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

      endif                    ! End write netcdf data, kopt=0

      return
      end
c
c

      subroutine wrtnetcdf_eps(netcdfnml)
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none !YuP[2020-01-14]
      integer error_code  !YuP[2020-01-14]
      integer neltmax,iray  !YuP[2020-01-14]

c     If the namelist variable dielectric_op="enabled", this routine
c     will write the complex dielectric tensor elements along the rays
c     into into existing netcdf file: netcdfnml 


      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'

      integer n1n2,istat
c     Storage tem1 is used in netcdf writes, including complex numbers.
cSAP090903
c      parameter (n1n2=2*nrelta*nraya)
c      real*8 tem1(n1n2) 
      real*8, pointer :: tem1(:)

c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid,twoid
      integer ray_dims(3),start(3),ray_count(3)

c----input
      character(*) netcdfnml ! input filename

      data start/1,1,1/

cSAP090903
      n1n2=2*nrelta*nrayl
c------------------------------------------
c     allocate pointers tem1
c-------------------------------------------
      allocate( tem1(1:n1n2),STAT=istat)
      call bcast(tem1,0.d0,SIZE(tem1))
c----------------------------------------------
 
c     Maximum number of ray elements per ray:
      neltmax=0
      do iray=1,nrayl
         neltmax=max(neltmax,nrayelt_nc(iray))
      enddo
      
      ray_count(1)=neltmax
      ray_count(2)=nrayl
      ray_count(3)=2

c.......................................................................
cl    1. Initialize part
c

C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c

c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.
      
      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)
      
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf2(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf2(integer ncid,integer error_code)
     
      call ncredf2(ncid,error_code)
      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c     dimension named 'lat'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c     returns dimension id.
c     p. 54 of netcdf manual



      neltmaxid=ncdid2(ncid,'neltmax',istatus)
      nraysid=ncdid2(ncid,'nrays',istatus)
      twoid=ncdid2(ncid,'two',istatus)
     
      ray_dims(1)=neltmaxid
      ray_dims(2)=nraysid
      ray_dims(3)=twoid

c     For ray data:
c-----define variables 

c-----Added 3rd dimension equal to 2 accomodates complex data.
  
      vid=ncvdef2(ncid,'cweps11',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps11',istatus)

      vid=ncvdef2(ncid,'cweps12',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps12',istatus)  

      vid=ncvdef2(ncid,'cweps13',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps13',istatus)

      vid=ncvdef2(ncid,'cweps21',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps21',istatus)

      vid=ncvdef2(ncid,'cweps22',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps22',istatus)  

      vid=ncvdef2(ncid,'cweps23',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps23',istatus)

      vid=ncvdef2(ncid,'cweps31',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps31',istatus)

      vid=ncvdef2(ncid,'cweps32',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps32',istatus)  

      vid=ncvdef2(ncid,'cweps33',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps33',istatus)
     
c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf2(ncid,istatus)
      call check_err(istatus)
   
c     End initialize

c.......................................................................
cl    1. Writing data
c
 
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c 
c--------------------------
c     write eps data
c--------------------------     
cSAP090903
c      call  storage_compl_2d(nrelta,nraya,nrayl,neltmax,
c     &cweps11_nc,tem1)       
      call  storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps11_nc,tem1)    
      call ncvid2(vid,ncid,'cweps11',istatus)     
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call  storage_compl_2d(nrelta,nraya,nrayl,neltmax,
c     &cweps12_nc,tem1)
      call  storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps12_nc,tem1)
      call ncvid2(vid,ncid,'cweps12',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call  storage_compl_2d(nrelta,nraya,nrayl,neltmax,
c     &cweps13_nc,tem1)
      call  storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps13_nc,tem1)
      call ncvid2(vid,ncid,'cweps13',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call  storage_compl_2d(nrelta,nraya,nrayl,neltmax,
c     &cweps21_nc,tem1)
      call  storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps21_nc,tem1)
      call ncvid2(vid,ncid,'cweps21',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call  storage_compl_2d(nrelta,nraya,nrayl,neltmax,
c     &cweps22_nc,tem1)
      call  storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps22_nc,tem1)
      call ncvid2(vid,ncid,'cweps22',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call  storage_compl_2d(nrelta,nraya,nrayl,neltmax,
c     &cweps23_nc,tem1)
      call  storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps23_nc,tem1)
      call ncvid2(vid,ncid,'cweps23',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call  storage_compl_2d(nrelta,nraya,nrayl,neltmax,
c     &cweps31_nc,tem1)
      call  storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps31_nc,tem1)
      call ncvid2(vid,ncid,'cweps31',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call  storage_compl_2d(nrelta,nraya,nrayl,neltmax,
c     &cweps32_nc,tem1)
      call  storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps32_nc,tem1)
      call ncvid2(vid,ncid,'cweps32',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090903
c      call  storage_compl_2d(nrelta,nraya,nrayl,neltmax,
c     &cweps33_nc,tem1)
      call  storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps33_nc,tem1)
      call ncvid2(vid,ncid,'cweps33',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

c     End write netcdf data

      deallocate( tem1,STAT=istat)

      return
      end
c
c


      subroutine storage_compl_2d(nrelta,nraya,nrayl,neltmax,
     &compl_2d,tem1)
c     put complex*16 2D array: compl_2d(neltmax,nrayl)
c     to real*8 1D array: tem1(2*nrelta*nraya)
c
c     It will be used in  netcdf subroutines

      
      implicit none

c-----input
      integer nrelta,nraya,nrayl,neltmax
      complex*16 compl_2d(nrelta,nraya)

c-----output
c     Storage tem1 is used in netcdf writes, including complex numbers.
      real*8 tem1(*) 

c-----locals    
      integer j,i,ii
      complex*16 cei

      cei=(0.d0,1.d0)

      ii=0
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=0.5d0*(compl_2d(i,j)+dconjg(compl_2d(i,j)))
         enddo
      enddo

      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=-cei*0.5d0*(compl_2d(i,j)-dconjg(compl_2d(i,j)))
         enddo
      enddo

      return
      end



      subroutine wrtnetcdf_grill_launch(netcdfnml)
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none !YuP[2020-01-14]
      integer error_code  !YuP[2020-01-14]

c     If the namelist variable istart=2 or =3  this routine
c     will write the ray starting coordinates 
c     into existing netcdf file: netcdfnml 



      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'
      include 'grill.i'
     
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid
      
c----input
      character(*) netcdfnml ! input filename

c.......................................................................
c     1. Initialize part
c
C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c
c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.
      
c      write(*,*)'netcdfr3d.f in wrtnetcdf_grill_launch'
c      write(*,*)'before ncid=ncopn(netcdfnml,NCWRITE,istatus)'
c      write(*,*)'netcdfnml = ',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)

c      write(*,*)'after ncid=ncopn2(netcdfnml,NCWRITE,) ncid= ',ncid
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf2(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf2(integer ncid,integer error_code)
     
      call ncredf2(ncid,error_code)
      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c     dimension named 'lat'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c     returns dimension id.
c     p. 54 of netcdf manual
     
      nraysid=ncdid2(ncid,'nrays',istatus)


c     For ray data:
c-----define variables 

c-----starting coordinate.
c      write(*,*)'before ncvdef2(z_starting) ncid,nraysid',ncid,nraysid
      vid=ncvdef2(ncid,'z_starting',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,10,
     +           'z_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'm',istatus)

c      write(*,*)'before ncvdef2(r_starting) ncid,nraysid',ncid,nraysid
      vid=ncvdef2(ncid,'r_starting',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,10,
     +           'r_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'm',istatus)

c      write(*,*)'before ncvdef2(phi_strting) ncid,nraysid',ncid,nraysid
      vid=ncvdef2(ncid,'phi_starting',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'phi_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'radians',istatus)

c      write(*,*)'before ncvdef2(N_toroidal_starting) ncid,nraysid',
c     +           ncid,nraysid
      vid=ncvdef2(ncid,'N_toroidal_starting',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     +           'N_toroidal_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)

c      write(*,*)'before ncvdef(N_pol_starting)ncid,nraysid',ncid,nraysid
      vid=ncvdef2(ncid,'N_poloidal_starting',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     +           'N_poloidal_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)

c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf2(ncid,istatus)
      call check_err(istatus)

c     End initialize

c.......................................................................
cl    1. Writing data
c
 
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c 
c--------------------------
c     write ray starting coordinates data for all rays
c--------------------------     
c      write(*,*)'before ncvid( z_starting ) ncid',ncid
      call ncvid2(vid,ncid,'z_starting',istatus)     
c      write(*,*)'after ncvid( z_starting vid',vid
      call ncvpt_doubl2(ncid,vid,1,nrayl,arzu0,istatus)

c      write(*,*)'before ncvid( r_starting ) ncid',ncid
      call ncvid2(vid,ncid,'r_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,arru0,istatus)
           
c      write(*,*)'before ncvid( phi_starting ) ncid',ncid
      call ncvid2(vid,ncid,'phi_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,arphiu0,istatus)
   
c      write(*,*)'before ncvid( N_toroidal_starting ) ncid',ncid
      call ncvid2(vid,ncid,'N_toroidal_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,arnphi,istatus)

c      write(*,*)'before ncvid( N_poloidal_starting ncid',ncid
      call ncvid2(vid,ncid,'N_poloidal_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,arntheta,istatus)

cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

c     End write netcdf data

      return
      end



      subroutine wrtnetcdf_EC_launch(netcdfnml)
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none !YuP[2020-01-14]
      integer error_code  !YuP[2020-01-14]
c     If the namelist variable istart=1  this routine
c     will write the ray starting coordinates 
c     into existing netcdf file: netcdfnml 



      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'
      include 'cone.i'
     
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid
      
c----input
      character(*) netcdfnml ! input filename

c.......................................................................
c     1. Initialize part
c
C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c
c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.
      
c      write(*,*)'netcdfr3d.f in wrtnetcdf_grill_launch'
c      write(*,*)'before ncid=ncopn2(netcdfnml,NCWRITE,istatus)'
c      write(*,*)'netcdfnml = ',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)

c      write(*,*)'after ncid=ncopn2(netcdfnml,NCWRITE,) ncid= ',ncid
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf2(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf2(integer ncid,integer error_code)
     
      call ncredf2(ncid,error_code)
      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c     dimension named 'lat'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c     returns dimension id.
c     p. 54 of netcdf manual

c      write(*,*)'before nraysid=ncdid(ncid,nrays,istatus)'

      nraysid=ncdid2(ncid,'nrays',istatus)

c      write(*,*)'after nraysid=ncdid(ncid,nrays,) nraysid ',nraysid

c     For ray data:
c-----define variables 

c-----starting coordinates.
      vid=ncvdef2(ncid,'z_starting',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,10,
     +           'z_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'm',istatus)

      vid=ncvdef2(ncid,'r_starting',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,10,
     +           'r_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'm',istatus)

      vid=ncvdef2(ncid,'phi_starting',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'phi_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'degrees',istatus)

      vid=ncvdef2(ncid,'alphast',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,52,
     &'toroidal angle measured from R-vector through source',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,'degree',istatus)

      vid=ncvdef2(ncid,'betast',NCDOUBLE,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,73,
     &'poloidal angle measured from z=constant plane,'//
     & ' pos above plane, neg below',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'degree',istatus)

c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf2(ncid,istatus)
      call check_err(istatus)

c     End initialize

c.......................................................................
cl    1. Writing data
c
 
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c--------------------------
c     write ray starting coordinates data for all rays
c--------------------------     
c      write(*,*)'before ncvid( z_starting)ncid ',ncid
      call ncvid2(vid,ncid,'z_starting',istatus)     
c      write(*,*)'after ncvid( z_starting vid',vid
      call ncvpt_doubl2(ncid,vid,1,nrayl,zstj,istatus)

c      write(*,*)'before ncvid( r_starting)ncid ',ncid
      call ncvid2(vid,ncid,'r_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,rstj,istatus)

c      write(*,*)'before ncvid( phi_starting) ncid ',ncid
      call ncvid2(vid,ncid,'phi_starting',istatus)     
cSAP080704      call ncvpt_doubl2(ncid,vid,1,nrayl,arphiu0,istatus)
      call ncvpt_doubl2(ncid,vid,1,nrayl,180.d0/pi*phistj,istatus)
      !YuP[2019-10-01] corrected phist --> phistj
   
c      write(*,*)'before ncvid( alphast ) ncid',ncid
      call ncvid2(vid,ncid,'alphast',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,180.d0/pi*alphaj,istatus)

c      write(*,*)'before ncvid( betast ) ncid',ncid
      call ncvid2(vid,ncid,'betast',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,180.d0/pi*betaj,istatus)

cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

c     End write netcdf data

      return
      end




      subroutine wrtnetcdf_emission(netcdfnml,nray_emission)
      implicit none   
      !implicit integer (i-n), real*8 (a-h,o-z)

c     If the namelist variable i_emission=1, this routine
c     will write the emiision data 
c     into into existing netcdf file: netcdfnml 

c     It write the emission data for several frequencies
c     but only for one launching point and one launching direction
c     The variable iray should be fixed

c     The data along the ray will be written for all used frequencies
c     for the central ray.
c     For non-central rays it will write only the data for the temperature
c     and data at the detector.
     
      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'
      include 'write.i'    
      include 'emissa.i'     

c      integer n1n2,nraya_nfreqa 
c      parameter (n1n2=2*nrelta*nfreqa)
c      parameter (nraya_nfreqa=nraya*nfreqa)

c      real*8 tem1(n1n2)
c      real*8 tem2(nraya_nfreqa) 
       real*8, pointer :: tem1(:) !(n1n2)
       real*8, pointer :: tem2(:) !(nraya_nfreqa) 
      integer n1n2,nraya_nfreqa 
    
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c     nfreq !(in one.i) is the number of frequencies
c     nfreq=1 gives detailed plot of emission from a single ray     
c     nfreq.gt.1 gives spectra covering the specified frequency
c     range
c   freq00=<freq01 (in one.i)are ratios of the minimal and maximal emission
c         frequencies to the central electron gyro-frequency
c         f_ce(rho=0)
c    wallr(in one.i) is the wall reflection coefficient, {0=< wallr =<1}
c    i_rrind  (in one.i)chooses the subroutine to calculate N_ray (default =0)
c    i_rrind=0  N_ray=N
c    i_rrind =1 from cold electron plasma using rrind   
c    i-rrind =2 from hot non-relativistic dispersion relation  

c    i_r_2nd_harm=1  (in one.i)
c    to calculate the major radius of the EC 2nd harmonic
c                =0 do not calculate (default =0) 
c              (to use drawemfr.in file i_r_2nd_harm=1)
c              (to use drawemf1.in file i_r_2nd_harm=0)
c  c nfreqa is the max value for nfreq (param.i)

c n_relt_harma (param.i) is the number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations  n=<n_relt_harm (one.i)
c n_relt_harm1a,n_relt_harm2a (param.i) are  the minimal and maximal
c               numbers of EC harmonics used in anti-hermitian
c               dielectric tensor calculations
c               n_relt_harm1=<n=<n_relt_harm2 (one.i)
c     from 'write.i'


c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,error_code,ncvdef2,ncdid2,ncddef2
 
      integer nfreqdim,neltmaxid,neltmax_initial_mesh_dim,
     & nraysid

      integer ray_dims(2),ray_dims_initial_mesh(2)
      integer start(2),ray_count(2)     !redefindable emission mesh 
                                        !with additional points
      integer ray_count_initial_mesh(2) !non-redefindable initial mesh,
                                        !no additional points
      integer temperature_dims(2)
      integer temperature_count(2)     
c-----input
      character(*) netcdfnml ! input filename
      integer nray_emission !number of rays at each frequency
      data start/1,1/
      save

c-----local for test
      integer nrelt_l,ifreq_l,i,
     &ifreq,iifreq,istat,neltmax_emis,nrayelt_l,
     &neltmax_emis_initial_mesh      
      real*8 ds,ds_old
c------------------------------------------------------
c     allocate tem1 and tem2
c----------------------------------------------------------
      n1n2=2*nrelta*nfreq
      allocate( tem1(1:n1n2),STAT=istat)
      call bcast(tem1,0.d0,SIZE(tem1))

cSAP090903
c      nraya_nfreqa=nraya*nfreq
      nraya_nfreqa=nrayl*nfreq
      allocate( tem2(1:nraya_nfreqa),STAT=istat)
      call bcast(tem2,0.d0,SIZE(tem2))
      
c......................................................................
cl    1. Initialize part, creating new netcdf file
c

C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c
c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

ctest
c      write(*,*)'netcdfr3d.f in wrtnetcdf_emission'
c      write(*,*)'nray_emission',nray_emission
c      do iray=1,nray_emission
c         do ifreq=1,nfreq
c            write(*,*)'iray,ifreq',iray,ifreq
c            write(*,*)'wi_0(iray,ifreq)',wi_0(iray,ifreq)
c            write(*,*)'wi_0_nc(iray,ifreq)',wi_0_nc(iray,ifreq)
c         enddo      
c      enddo
cendtest
c      write(*,*)'before ncid=ncopn2(netcdfnml,NCWRITE,istatus)'
c      write(*,*)'netcdfnml = ',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)

c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf2(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf2(integer ncid,integer error_code)
     
c      write(*,*)'before ncredf2'
      call ncredf2(ncid,error_code)
c      write(*,*)'after  ncredf2'
      call check_err(istatus)

c-------------------------------------------------------------------
c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c     dimension named 'lat'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c     returns dimension id.
c     p. 54 of netcdf manual
C-----------------------------------------------------------------------
c
c      cl     1.1 create netCDF file and define dimensions,variables
c          and attribute
c
c.......................................................................
cl    1.1.2 define dimensions
c     integer function ncddef(ncid,dim_name,dim_siz,error_code)
c       returns dimension id.
c     p. 67 of netcdf-2 manual
c----------------------------------------------------------------------
c     Maximum number of ray elements( with additinal points) per ray
c     for all frequencies :
c----------------------------------------------------------------------
      neltmax_emis=0

      do ifreq=1,nfreq         
         neltmax_emis=max(neltmax_emis,nrayelt_emis_nc(ifreq))        
      enddo

      if (neltmax_emis.eq.0) neltmax_emis=1

      ray_count(1)=neltmax_emis
      ray_count(2)=nfreq
c      write(*,*)'ray_count',ray_count
      neltmaxid=ncddef2(ncid,'neltmax_emis',neltmax_emis,istatus)
c      write(*,*)'neltmaxid',neltmaxid
      nfreqdim=ncddef2(ncid,'nfreq',nfreq,istatus)
c      write(*,*)'nfreqdim',nfreqdim
      nraysid=ncddef2(ncid,'nray_emission',nray_emission,istatus)
c      write(*,*)'nray_emission',nray_emission
c      write(*,*)'nraysid',nraysid

      ray_dims(1)=neltmaxid
      ray_dims(2)=nfreqdim
c      write(*,*)'ray_dims',ray_dims

c---------------------------------------------------------------------
c     the variavles to write the data at the end of the all rays
c     at the detector
c---------------------------------------------------------------------
      
      temperature_count(1)=nray_emission
      temperature_count(2)=nfreq
      
      temperature_dims(1)=nraysid
      temperature_dims(2)=nfreqdim

c----------------------------------------------------------------------
c     Maximum number of ray elements( with no additinal points) per ray
c     for all frequencies :
c----------------------------------------------------------------------
      neltmax_emis_initial_mesh=0

      do ifreq=1,nfreq
         neltmax_emis_initial_mesh=max(neltmax_emis_initial_mesh,
     &   nrayelt_emis_initial_mesh_nc(ifreq))
      enddo

      if (neltmax_emis_initial_mesh.eq.0) neltmax_emis_initial_mesh=1

      ray_count_initial_mesh(1)=neltmax_emis_initial_mesh
      ray_count_initial_mesh(2)=nfreq

      neltmax_initial_mesh_dim=ncddef2(ncid,'neltmax_emis_initial_mesh',
     &neltmax_emis_initial_mesh,istatus)
   
      ray_dims_initial_mesh(1)=neltmax_initial_mesh_dim
      ray_dims_initial_mesh(2)=nfreqdim

c-------------------------------------------------------------------------
c     Run Descriptive Parameters
c-------------------------------
c-----the emission input data set in the input files: genray.dat or genray.in
      
      vid=ncvdef2(ncid,'nfreq',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     +            'The number of the emission frequencies',istatus)
      call check_err(istatus)
     
      vid=ncvdef2(ncid,'nray_emission',NCLONG,0,0,istatus)      
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32,
     +            'The number of rays per frequency',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'freq00',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,84,
     +            'The ratio of the minimal emission frequencies to '//
     +'the central electron gyro-frequency',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'freq01',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,84,
     +            'The ratio of the maximal emission frequencies to '// 
     +'the central electron gyro-frequency',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wallr',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31,
     +            'The wall reflection coefficinet',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'tol_emis',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,46,
     +            'The tolerance parameter to add new mesh points',
     +            istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'i_rrind',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44,
     +            'It chooses the subroutine to calculate N_ray',
     +            istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'i_r_2nd_harm',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,68,
     +            'The parameter (=1) to calculate the major radius '//
     +'of EC 2nd harmonics',istatus)
      call check_err(istatus)
 
      vid=ncvdef2(ncid,'freqncy0',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     +            'The central electron gyro-frequency ',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,'GHZ',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'nrayelt_emis_nc',NCLONG,1,nfreqdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,65,
     +'the number of points along the central ray with ' //
     +'additional points',
     +           istatus)
      call check_err(istatus)
   
      vid=ncvdef2(ncid,'wsn_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,63,
     +'total distance along the central ray ' // 
     +'with the additional points',
     +           istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wz_em_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,64,
     +'vertical height along the central ray with the additional points'
     +           ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)
      call check_err(istatus)
     
      vid=ncvdef2(ncid,'wr_em_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,61,
     +'major radius along the central ray with the additional points',
     +           istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)
      call check_err(istatus)
    
      vid=ncvdef2(ncid,'wphi_em_nc',NCDOUBLE,2,ray_dims,istatus)     
      call ncaptc2(ncid,vid,'long_name',NCCHAR,63,
     +'toroidal angle along the central ray with the additional points',
     +           istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'radians',istatus)
      call check_err(istatus)
    
      vid=ncvdef2(ncid,'wcnpar_em_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,73,
     +'parallel refractive index along the cenral ray '//
     +'with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wcnper_em_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,80,
     +'perpendicular refractive index along the central ray '//
     +' with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wal_emis_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,72,
     +'absorption coefficient along the central ray '//
     +' with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           '(1/cm)',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wj_emis_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,70,
     +'emission coefficient along the central ray '//
     +' with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           '(erg*sec/cm**3)',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wnray_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,70,
     +'ray refractive index along the central ray '//
     +' with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'win_sn_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,102,
     +'emission I_n at the detector side of nth bin at s=s_n along'//
     +' the central ray with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'win_0_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,102,
     +'emission I_n_0 at the plasma edge of nth bin at s=s_n along'//
     +' the central ray with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_specific_intensity_nc',NCDOUBLE,2,ray_dims,
     +istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,96,
     +'specific intensity of the emission I_n_0/delta(s) '//
     +'from the central ray with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)
   
      vid=ncvdef2(ncid,'wi_0_nc',NCDOUBLE,2,temperature_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,56,
     +'emission I_0 at plasma boundary from each ray at s=s_1=0'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wi_0sn_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,92,
     +'sum{k=1,n}[in_0(k)] emission at the plasma boundary '//
     +'from the parts 0<s<sn of the central ray'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wtemp_em_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,68,
     +'plasma temperature along the central ray '//
     +' with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wtemp_rad_em_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,71,
     +'radiation temperature along the central ray '//
     +' with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wtemp_rad_fr_nc',NCDOUBLE,2,temperature_dims,
     +istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,88,
     +'emission temperature for multi-frequency case '//
     +'for each ray  =2*pi*c**2*I_0t/frequency**2'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wtemp_rad_fr_wall_nc',NCDOUBLE,2,
     +temperature_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,108,
     +'emission temperature for multi-frequency case with wall '//
     +'reflection for each ray =2*pi*c**2*I_0t/frequency**2'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'waveraged_temp_rad_fr_nc',NCDOUBLE,1,nfreqdim,
     +istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,136,
     +'averaged other all rays emission temperature '//
     +' for multi-frequency case without wall '//
     +'reflection for each ray =2*pi*c**2*I_0t/frequency**2'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'waveraged_temp_rad_fr_wall_nc',NCDOUBLE,1,
     +nfreqdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,133,
     +'averaged other all rays emission temperature '//
     +' for multi-frequency case with wall '//
     +'reflection for each ray =2*pi*c**2*I_0t/frequency**2'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wtemp_pl_fr_nc',NCDOUBLE,2,temperature_dims,
     &istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,81,
     +'plasma temperature for multi-frequency case'//
     +'in point of max emission for each rays'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wr0_em_nc',NCDOUBLE,2,temperature_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,67,
     +'major radius of max emission for multi-frequency case '//
     +' for each ray'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'cm',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'wrho0_em_nc',NCDOUBLE,2,temperature_dims,
     +            istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,67,
     +'small radius of max emission for multi-frequency case '//
     +' for each ray'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'normalized',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'wtaun_em_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,84,
     +'tau_n =integral{0,s_n}(al_emis(s)ds)'//  
     +'along the central ray with the additional points'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wfreq_nc',NCDOUBLE,1,nfreqdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,45,
     +'emission frequencies for multi-frequency case'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'GHZ',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'wtau_em_nc',NCDOUBLE,2,temperature_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,55,
     +'the optical depth for multi-frequency case for each ray' 
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wi_0t_nc',NCDOUBLE,2,temperature_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,103,
     +'emission at the boundary  with wall reflection'//
     +'=I_0/ (1.d0-wallr*dexp(-wtau_em(iray,ifreq)) for each ray'
     +          ,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'nrayelt_emis_initial_mesh_nc',
     +NCLONG,1,nfreqdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,90,
     +'the number of non-refindable output points along '//
     +'the central ray for different frequencies',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wr_emis_initial_mesh_nc',
     +NCDOUBLE,2,ray_dims_initial_mesh,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,92,
     +'major radius of non-refindable output points along '// 
     +'the central ray for different frequencies',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wp_perpmax_dmc_nc',NCDOUBLE,2,
     &ray_dims_initial_mesh,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,122,
     +'p_perpendicular_max/(rest_me*c_light) at the second harmonic'// 
     +'resonance line along the central ray without additional points',
     +istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wp_parmax_dmc_nc',NCDOUBLE,2,
     &ray_dims_initial_mesh,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,117,
     +'p_parallel_max/(rest_me*c_light) at the second harmonic'// 
     +'resonance line along the central ray without additional points',
     +istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,'non-dimensional',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wp_parmin_dmc_nc',NCDOUBLE,2,
     &ray_dims_initial_mesh,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,113,
     +'p_parallel_min/(rest_me*c_light) at the second harmonic'// 
     +'resonance line along the central without additional points',
     +istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,'non-dimensional',istatus)
      call check_err(istatus)


      if (i_r_2nd_harm.eq.1) then  
        vid=ncvdef2(ncid,'wr_2nd_harm_nc',NCDOUBLE,2,temperature_dims,
     +istatus)
        call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +  'major radius at second harmonic point for each ray'
     +          ,istatus)
        call check_err(istatus)
        call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'cm',istatus)
        call check_err(istatus)

        vid=ncvdef2(ncid,'wtemp_2nd_harm_nc',NCDOUBLE,2,
     +temperature_dims,istatus)
        call ncaptc2(ncid,vid,'long_name',NCCHAR,56,
     +  'plasma temperature at second harmonic point for each ray'
     +          ,istatus)
        call check_err(istatus)
        call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'KeV',istatus)
        call check_err(istatus)
      endif !  i_r_2nd_harm.eq.1

      vid=ncvdef2(ncid,'wdye_nc',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,107,
     +'ratio of the wave frequency to the electron gyro frequency ' //
     +'along the central ray with the additional points',
     +           istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'non-dimentional',istatus)
      call check_err(istatus)


c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual 

      call ncendf2(ncid,istatus)
      call check_err(istatus)

c     End initialize

c.......................................................................
cl    1. Writing data
c
 
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c-------------------------- 
      call ncvid2(vid,ncid,'nfreq',istatus)     
      call ncvpt_int2(ncid,vid,1,1,nfreq,istatus)

      call ncvid2(vid,ncid,'nray_emission',istatus)     
      call ncvpt_int2(ncid,vid,1,1,nray_emission,istatus)


      call ncvid2(vid,ncid,'freq00',istatus)     
      call ncvpt_doubl2(ncid,vid,1,1,freq00,istatus)

      call ncvid2(vid,ncid,'freq01',istatus) 
      call ncvpt_doubl2(ncid,vid,1,1,freq01,istatus)

      call ncvid2(vid,ncid,'wallr',istatus)     
      call ncvpt_doubl2(ncid,vid,1,1,wallr,istatus)

      call ncvid2(vid,ncid,'tol_emis',istatus)     
      call ncvpt_doubl2(ncid,vid,1,1,tol_emis,istatus)

      call ncvid2(vid,ncid,'i_rrind',istatus)     
      call ncvpt_int2(ncid,vid,1,1,i_rrind,istatus)

      call ncvid2(vid,ncid,'i_r_2nd_harm',istatus)     
      call ncvpt_int2(ncid,vid,1,1,i_r_2nd_harm,istatus)

      call ncvid2(vid,ncid,'freqncy0',istatus)       
      call ncvpt_doubl2(ncid,vid,1,1,freqncy0,istatus)

      call ncvid2(vid,ncid,'nrayelt_emis_nc',istatus)
      call ncvpt_int2(ncid,vid,1,nfreq,nrayelt_emis_nc,istatus)

      call pack21(wsn_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)

c_test_begin
c      write(*,*)'netcdfr3d.f neltmax_emis,nfreq',neltmax_emis,nfreq
c      do ifreq_l=1,nfreq
c         write(*,*)'ifreq_l',ifreq_l
c         do nrelt_l=1,neltmax_emis
c          write(*,*)'nrelt_l,wsn_nc(nrelt_l,ifreq_l)',
c     &    nrelt_l,wsn_nc(nrelt_l,ifreq_l)
c         enddo
c      enddo
c      do i=1,nfreq*neltmax_emis
c         write(*,*)'i,tem1(i)',i,tem1(i)
c      enddo
c_test_end
         
      call ncvid2(vid,ncid,'wsn_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wz_em_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'wz_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wr_em_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'wr_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wphi_em_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'wphi_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wcnpar_em_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,
     &nfreq)
      call ncvid2(vid,ncid,'wcnpar_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wcnper_em_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,
     &nfreq)
      call ncvid2(vid,ncid,'wcnper_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wal_emis_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'wal_emis_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wj_emis_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'wj_emis_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wnray_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'wnray_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(win_sn_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'win_sn_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(win_0_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'win_0_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

c----------------------------------------------------------------
c     specific intensity calculations
c----------------------------------------------------------------
c      write(*,*)'netcdfr3d.f nfreq',nfreq 
      do ifreq=1,nfreq

        do nrayelt_l=1,nrayelt_emis_nc(ifreq)-1
 
           ds=wsn_nc(nrayelt_l+1,ifreq)-wsn_nc(nrayelt_l,ifreq)
           if(ds.gt.1.d-18) then
             w_specific_intensity_nc(nrayelt_l,ifreq)=
     &       win_0_nc(nrayelt_l,ifreq)/ds
           else
             w_specific_intensity_nc(nrayelt_l,ifreq)=
     &       win_0_nc(nrayelt_l,ifreq)/ds_old
           endif
c
c           write(*,*)'nrayelt_l,wsn_nc(nrayelt_l+1,ifreq)',
c     &     nrayelt_l,wsn_nc( nrayelt_l+1,ifreq)
c           write(*,*)'wsn_nc(nrayelt_l,ifreq),ds',
c     &     wsn_nc( nrayelt_l,ifreq),ds

c           write(*,*)'win_0_nc(nrayelt_l,ifreq)',
c     &     win_0_nc(nrayelt_l,ifreq)

c           write(*,*)'w_specific_intensity_nc(nrayelt_l,ifreq)',
c     &     w_specific_intensity_nc(nrayelt_l,ifreq)

           ds_old=ds
        enddo
      enddo
c      write(*,*)'before pack21 w_specific_intensity_nc'

      call pack21(w_specific_intensity_nc,1,nrelta,1,nfreq,tem1,
     &neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'w_specific_intensity_nc',istatus)      
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

c      write(*,*)'after pack21 w_specific_intensity_nc'

      call pack21(wdye_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'wdye_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

c-------------------------------------------------------------------
      call pack21(wi_0sn_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'wi_0sn_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)  
 
      call pack21(wtemp_em_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,nfreq)
      call ncvid2(vid,ncid,'wtemp_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
      
      call pack21(wtemp_rad_em_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,
     &nfreq)
      call ncvid2(vid,ncid,'wtemp_rad_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
       
c-----the data for all rays
cSAP090903 
c      call pack21(wi_0_nc,1,nraya,1,nfreq,tem2,nray_emission,nfreq)
      call pack21(wi_0_nc,1,nrayl,1,nfreq,tem2,nray_emission,nfreq)
      call ncvid2(vid,ncid,'wi_0_nc',istatus) 
      call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,
     &istatus)

cSAP090903       
c      call pack21(wtemp_rad_fr_nc,1,nraya,1,nfreq,tem2,
c     &nray_emission,nfreq)
      call pack21(wtemp_rad_fr_nc,1,nrayl,1,nfreq,tem2,
     &nray_emission,nfreq)
      call ncvid2(vid,ncid,'wtemp_rad_fr_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,istatus)
        
c      do iray_l=1,nray_emission
c        do ifreq_l=1,nfreq
c        write(*,*)'iray_l,ifreq_l,wtemp_rad_fr_wall_nc(iray_l,ifreq_l)'
c     &            ,iray_l,ifreq_l,wtemp_rad_fr_wall_nc(iray_l,ifreq_l)     
c        enddo
c      enddo
cSAP090903   
c      call pack21(wtemp_rad_fr_wall_nc,1,nraya,1,nfreq,tem2,
c     &nray_emission,nfreq)
      call pack21(wtemp_rad_fr_wall_nc,1,nrayl,1,nfreq,tem2,
     &nray_emission,nfreq)
      call ncvid2(vid,ncid,'wtemp_rad_fr_wall_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,istatus)
  
      call ncvid2(vid,ncid,'waveraged_temp_rad_fr_nc',istatus)  
      call ncvpt_doubl2(ncid,vid,1,nfreq, waveraged_temp_rad_fr_nc,
     +                  istatus)    

      call ncvid2(vid,ncid,'waveraged_temp_rad_fr_wall_nc',istatus)
      call ncvpt_doubl2(ncid,vid,1,nfreq, waveraged_temp_rad_fr_wall_nc,
     +                  istatus)
  
c      do iray_l=1,nray_emission
c        do ifreq_l=1,nfreq
c        write(*,*)'iray_l,ifreq_l,wtemp_pl_fr_nc(iray_l,ifreq_l)'
c     &            ,iray_l,ifreq_l,wtemp_pl_fr_nc(iray_l,ifreq_l)    
c        enddo
c      enddo
cSAP090903   
c      call pack21(wtemp_pl_fr_nc,1,nraya,1,nfreq,tem2,
c     &nray_emission,nfreq)
      call pack21(wtemp_pl_fr_nc,1,nrayl,1,nfreq,tem2,
     &nray_emission,nfreq)
      call ncvid2(vid,ncid,'wtemp_pl_fr_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,istatus)

cSAP090903   
c      call pack21(wr0_em_nc,1,nray,1,nfreq,tem2,nray_emission,nfreq)
      call pack21(wr0_em_nc,1,nrayl,1,nfreq,tem2,nray_emission,nfreq)
      call ncvid2(vid,ncid,'wr0_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,istatus)

cSAP090903   
c      call pack21(wrho0_em_nc,1,nraya,1,nfreq,tem2,nray_emission,nfreq)
      call pack21(wrho0_em_nc,1,nrayl,1,nfreq,tem2,nray_emission,nfreq)
      call ncvid2(vid,ncid,'wrho0_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,istatus)
    
      call pack21(wtaun_em_nc,1,nrelta,1,nfreq,tem1,neltmax_emis,
     &nfreq)   
      call ncvid2(vid,ncid,'wtaun_em_nc',istatus)    
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call ncvid2(vid,ncid,'wfreq_nc',istatus)
      call ncvpt_doubl2(ncid,vid,1,nfreq,wfreq_nc,istatus)

cSAP090903   
c      call pack21(wtau_em_nc,1,nraya,1,nfreq,tem2,nray_emission,
c     &nfreq)
      call pack21(wtau_em_nc,1,nrayl,1,nfreq,tem2,nray_emission,
     &nfreq)
      call ncvid2(vid,ncid,'wtau_em_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,istatus)

cSAP090903   
c      call pack21(wi_0t_nc,1,nraya,1,nfreq,tem2,nray_emission,
c     &nfreq)
      call pack21(wi_0t_nc,1,nrayl,1,nfreq,tem2,nray_emission,
     &nfreq)
      call ncvid2(vid,ncid,'wi_0t_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,istatus)
c-----------------------------------------------------------------------------------
c     resonance curve:
c     nrayelt_emis_initial_mesh_nc is the number of initial (non-refindable) mesh points 
c     wr_emis_initial_mesh_nc  is major radius initial (non-refindable) mesh
c--------------------------------------------------------------------------------- 
      call ncvid2(vid,ncid,'nrayelt_emis_initial_mesh_nc',istatus)
      call ncvpt_int2(ncid,vid,1,nfreq,nrayelt_emis_initial_mesh_nc,
     +                istatus)


c      do ifreq=1,nfreq
c       write(*,*)'write emission nfreq:ifreq', ifreq
c       write(*,*)'write emission nrayelt_emis_initial_mesh_nc(ifreq)',
c     & nrayelt_emis_initial_mesh_nc(ifreq)
c       n_=nrayelt_emis_initial_mesh_nc(ifreq)
c       do i=1,n_
c       write(*,*)'write emission nfreq:i,wr_emis_initial_mesh_nc'
c     &,i,wr_emis_initial_mesh_nc(i,ifreq)
c       enddo
c      enddo

      call pack21(wr_emis_initial_mesh_nc,1,nrelta,1,nfreq,tem1,
     &neltmax_emis_initial_mesh,nfreq)
      call ncvid2(vid,ncid,'wr_emis_initial_mesh_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count_initial_mesh,tem1,
     +                  istatus)


       iifreq=0
      do ifreq=1,nfreq
        do i=1,neltmax_emis_initial_mesh     
          iifreq=iifreq+1
          if (i.le.nrayelt_emis_initial_mesh_nc(ifreq)) then   
            tem1(iifreq)=wp_perpmax_dmc_nc(i,2,ifreq)
          else
            tem1(iifreq)=0.d0
          endif
        enddo
      enddo 
 
      call ncvid2(vid,ncid,'wp_perpmax_dmc_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count_initial_mesh,tem1,
     +                  istatus)
    
      iifreq=0
      do ifreq=1,nfreq
        do i=1,neltmax_emis_initial_mesh     
          iifreq=iifreq+1
          if (i.le.nrayelt_emis_initial_mesh_nc(ifreq)) then   
            tem1(iifreq)=wp_parmax_dmc_nc(i,2,ifreq)
          else
            tem1(iifreq)=0.d0
          endif
        enddo
      enddo 
 
      call ncvid2(vid,ncid,'wp_parmax_dmc_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count_initial_mesh,tem1,
     +                  istatus)

  
      iifreq=0
      do ifreq=1,nfreq
        do i=1,neltmax_emis_initial_mesh     
          iifreq=iifreq+1
          if (i.le.nrayelt_emis_initial_mesh_nc(ifreq)) then   
            tem1(iifreq)=wp_parmin_dmc_nc(i,2,ifreq)
          else
            tem1(iifreq)=0.d0
          endif
        enddo
      enddo 
 
      call ncvid2(vid,ncid,'wp_parmin_dmc_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count_initial_mesh,tem1,
     +                  istatus)

      if (i_r_2nd_harm.eq.1) then
c         write(*,*)'in netcdfr3d.f nrayl,nray_emission',
c     &                            nrayl,nray_emission
c         write(*,*)'in netcdfr3d.f nfreq',nfreq

cSAP090903
c         call pack21(wr_2nd_harm_nc,1,nraya,1,nfreq,tem2,nray_emission,
c     &   nfreq)
         call pack21(wr_2nd_harm_nc,1,nrayl,1,nfreq,tem2,nray_emission,
     &   nfreq)
         call ncvid2(vid,ncid,'wr_2nd_harm_nc',istatus)     
         call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,
     +                     istatus)

cSAP090903
c         call pack21(wtemp_2nd_harm_nc,1,nraya,1,nfreq,tem2,
c     &   nray_emission,nfreq)
         call pack21(wtemp_2nd_harm_nc,1,nrayl,1,nfreq,tem2,
     &   nray_emission,nfreq)
         call ncvid2(vid,ncid,'wtemp_2nd_harm_nc',istatus)
         call ncvpt_doubl2(ncid,vid,start,temperature_count,tem2,
     +                     istatus)

      endif !i_r_2nd_harm.eq.1

cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

c     End write netcdf data

      deallocate( tem1,STAT=istat)
      deallocate( tem2,STAT=istat)

      return
      end 


      subroutine put_emission_in_writencdf_i(ifreq,iray_start_point)
c-----put emission data in writencdf.i 
      implicit none  
      !implicit integer (i-n), real*8 (a-h,o-z)

      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'
      include 'write.i'
      include 'emissa.i'
c-----input
      integer ifreq
c     &,iray
cSm060203
      integer iray_start_point !number of ray launching point
              !it will put ray trajectories data for central ray only
              !at iray_start_point=1, 
              !for iray_start_point > 1 the rays will not be put in 
              !netcdf file
c-----local
      integer i,n,iray_loc,istatus  
      real*8 cm,cnpar,cnper,cnr,cnt,cnz,ds,dc,phi,r,z

c-----external
      real*8 b,gamma1,cn,y

c      write(*,*)'put_emission_in_writencdf_i ifreq,nrayelt_emis',
c     &ifreq,nrayelt_emis  
c      write(*,*)'begin of put_emission_in_writencdf_i'     
c      write(*,*)'iray_start_point,ifreq',iray_start_point,ifreq
c      write(*,*)'wi_0(iray_start_point,ifreq)',
c     &           wi_0(iray_start_point,ifreq)
        
      if( iray_start_point.eq.1) then 
c---------------------------------------------------------------
c     central ray data will be saved for writing in netcdf file
c-----------------------------------------------------------------
      nrayelt_emis_nc(ifreq)= nrayelt_emis
      iray_loc=1 

      do i=1,nrayelt_emis
         wsn_nc(i,ifreq)=wsn(i)
         wr_em_nc(i,ifreq)=wr_em(i)
         wz_em_nc(i,ifreq)=wz_em(i)
         wphi_em_nc(i,ifreq)=wphi_em(i)
         wal_emis_nc(i,ifreq)=wal_emis(i)
         wj_emis_nc(i,ifreq)=wj_emis(i)
         wnray_nc(i,ifreq)=wnray(i)
         wtemp_em_nc(i,ifreq)=wtemp_em(i)
         wtemp_rad_em_nc(i,ifreq)=wtemp_rad_em(i)
         if (i.ne.nrayelt_emis) then
            win_sn_nc(i,ifreq)=win_sn(i)
            win_0_nc(i,ifreq)=win_0(i)
            wi_0sn_nc(i,ifreq)=wi_0sn(i)

c            write(*,*)'netcdfr3d.f put_emission_in_writencdf_',
c     &                 'i,ifreq',i,ifreq

            wtaun_em_nc(i,ifreq)=wtaun_em(i)
         endif
c-----------------------------------------------
c        calculate N_par and N_per
c-----------------------------------------------
         r=wr_em(i)
         z=wz_em(i)
         phi=wphi_em(i)
         cnz= wcnz_em(i)
         cnr=wcnz_em(i)
         cm=wcm_em(i)
         bmod=b(z,r,phi)
         gam=gamma1(z,r,phi,cnz,cnr,cm)
         ds=dsin(gam)
         dc=dcos(gam)
         cnt=cn(r,cnz,cnr,cm)
         cnpar=cnt*dc
         cnper=cnt*ds
c----------------------------------------------
         wcnpar_em_nc(i,ifreq)=cnpar
         wcnper_em_nc(i,ifreq)=cnper 

c         write(*,*)'put_emission_in_writencdf_i i,ifreq',i,ifreq
c         write(*,*)'wcnpar_em_nc(i,ifreq),wcnper_em_nc(i,ifreq)',
c     &   wcnpar_em_nc(i,ifreq),wcnper_em_nc(i,ifreq)
c         write(*,*)'wi_0sn_nc(i,ifreq)',wi_0sn_nc(i,ifreq)
c----------------------------------------------------
c         calculate the ratio omega/omega_ce along the ray
c----------------------------------------------------
          wdye_nc(i,ifreq)=1.d0/y(z,r,phi,1)

      enddo

      nrayelt_emis_nc(ifreq)=nrayelt_emis
      nrayelt_emis_initial_mesh_nc(ifreq)= nrayelt
      wfreq_nc(ifreq)=wfreq(ifreq)

cSm060203
c----------------------------------------------------------------     
c      wi_0_nc(ifreq) =wi_0(iray_loc,ifreq)
c      wtemp_rad_fr_nc(ifreq)=wtemp_rad_fr(ifreq)
c      wtemp_rad_fr_wall_nc(ifreq)=wtemp_rad_fr(ifreq)
c      wtemp_pl_fr_nc(ifreq)=wtemp_pl_fr(ifreq)
c      wr0_em_nc(ifreq)=wr0_em(ifreq)
c      wrho0_em_nc(ifreq)=wrho0_em(ifreq)
c      wfreq_nc(ifreq)=wfreq(ifreq)
c      wtau_em_nc(ifreq)=wtau_em(1,ifreq)
c      wi_0t_nc(ifreq)=wi_0t(1,ifreq)

c      if (i_r_2nd_harm.eq.1) then
c        wr_2nd_harm_nc(ifreq)=wr_2nd_harm(ifreq)
c        wtemp_2nd_harm_nc(ifreq)=wtemp_2nd_harm(ifreq)
c      endif 
 
c     write(*,*)'put_emission_in_writencdf_i ifreq,wi_0_nc(ifreq)',
c     &ifreq,wi_0_nc(ifreq)
    
c      nrayelt_emis_initial_mesh_nc(ifreq)= nrayelt
c------------------------------------------------------------------       
c      write(*,*)'put_emission_in_writencdf_i nrayelt',nrayelt    
      do i=1,nrayelt
cSAP090805
c       wr_emis_initial_mesh_nc(i,ifreq)=wr_em(i)
       wr_emis_initial_mesh_nc(i,ifreq)=wr(i)/(100.d0*r0x)
 
c       write(*,*)'put_emission_in_writencdf_i ifreq,i,wr(i),r0x',
c     & ifreq,i,wr(i),r0x

       do n=n_relt_harm1,n_relt_harm2     
         wp_perpmax_dmc_nc(i,n,ifreq)=wp_perpmax_dmc(i,n)
         wp_parmin_dmc_nc(i,n,ifreq)=wp_parmin_dmc(i,n)
         wp_parmax_dmc_nc(i,n,ifreq)=wp_parmax_dmc(i,n)
       enddo
      enddo

c      write(*,*)'put_emission_in_writencdf_i ifreq,nrayelt',
c     &ifreq,nrayelt
c      write(*,*)'put_emission:  ifreq wr_emis_initial_nc(i,ifreq)',
c     &(wr_emis_initial_mesh_nc(i,ifreq),i=1,nrayelt)

c-----end_if(iray_start_point.eq.1) case      
      endif

c---------------------------------------------------------------
c     some data for all rays will be saved for writing in netcdf file
c-----------------------------------------------------------------     
c     wi_0_nc(iray_start_point,ifreq) =wi_0(iray_loc,ifreq)
    
      wi_0_nc(iray_start_point,ifreq)=wi_0(iray_start_point,ifreq)
      
      wtemp_rad_fr_nc(iray_start_point,ifreq)=wtemp_rad_fr(ifreq)
      
      wtemp_rad_fr_wall_nc(iray_start_point,ifreq)=
     &wtemp_rad_fr_wall(ifreq)
      
      wtemp_pl_fr_nc(iray_start_point,ifreq)=wtemp_pl_fr(ifreq)
    
      wr0_em_nc(iray_start_point,ifreq)=wr0_em(ifreq)
    
      wrho0_em_nc(iray_start_point,ifreq)=wrho0_em(ifreq)
    
      wtau_em_nc(iray_start_point,ifreq)=wtau_em(1,ifreq)
     
      wi_0t_nc(iray_start_point,ifreq)=wi_0t(1,ifreq)
    
      if (i_r_2nd_harm.eq.1) then
         wr_2nd_harm_nc(iray_start_point,ifreq)=wr_2nd_harm(ifreq)
c         write(*,*)'wr_2nd_harm(ifreq)',wr_2nd_harm(ifreq)
         wtemp_2nd_harm_nc(iray_start_point,ifreq)=
     &                           wtemp_2nd_harm(ifreq)
c         write(*,*)'wtemp_2nd_harm(ifreq)', wtemp_2nd_harm(ifreq)
      endif 

      return
      end
   
      subroutine read_nc(netcdfnm_out)
c-----read output *.nc file to check the resonance curves

c     implicit none
      implicit none
      !implicit integer (i-n), real*8 (a-h,o-z)
      character(*) netcdfnm_out ! input filename
c      character*9 netcdfnm ! input filename

      include 'netcdf.inc'
      include 'param.i' 
      include 'one.i'
      include 'writencdf.i'
      include 'emissa.i'

c      integer n1n2
c      parameter (n1n2=2*nrelta*nfreqa)
c      real*8 tem1(n1n2) 
      integer n1n2
      real*8, pointer :: tem1(:)      !(n1n2)
 
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2

      integer nfreqdim,neltmaxid,neltmax_initial_mesh_dim

      integer ray_dims(2),ray_dims_initial_mesh(2)

      integer start(2),ray_count(2),ray_count_initial_mesh(2)

      character*128 name

      integer i,j,
     &istat,ifreq,neltmax,neltmax_id,nfreq_id,
     & neltmax_emis_initial_mesh
       
      data start/1,1/
       
c      write(*,*)'********************************************'
c      write(*,*)'in read_nc'
c      write(*,*)'********************************************'
c      write(*,*)'netcdfnm_out= ',netcdfnm_out
c-------------------------------------------------
c     Open previous netCDF file
c-------------------------------------------------
      call ncopn2(netcdfnm_out,NCNOWRIT,ncid,istatus)
      call check_err(istatus)
c      write(*,*)'ncid',ncid
c------------------------------------------------
c     Put Open NetCDF file into Define Mode
c------------------------------------------------
c      call ncredf2(ncid,error_code)
c      call check_err(istatus)
c----------------------------------------------------
c     determine dimension ID of  dimension named '...'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c     read ID name and sizes
c-----------------------------------------------------
      neltmax_id =ncdid2(ncid,'neltmax',istatus)
      call check_err(istatus)      
      call ncdinq2(ncid,neltmax_id,name,neltmax,istatus)
      call check_err(istatus)   
c      write(*,*)'neltmax_id,neltmax',neltmax_id,neltmax

      nfreq_id =ncdid2(ncid,'nfreq',istatus)
      call check_err(istatus)   
      call ncdinq2(ncid,nfreq_id,name,nfreq,istatus)
      call check_err(istatus)  
c      write(*,*)'nfreq_id,nfreq',nfreq_id,nfreq
     
c-----allocate tem1
      n1n2=nrelta*nfreq
      allocate( tem1(1:n1n2),STAT=istat)
      call bcast(tem1,0.d0,SIZE(tem1))
c--------------------------------------
      neltmax_initial_mesh_dim=ncdid2(ncid,'neltmax_emis_initial_mesh',
     &istatus)
      call check_err(istatus)
      call ncdinq2(ncid,neltmax_initial_mesh_dim,name,
     &neltmax_emis_initial_mesh,istatus)
      call check_err(istatus)   
c      write(*,*)'neltmax_emis_initial_mesh',neltmax_emis_initial_mesh

      ray_count_initial_mesh(1)=neltmax_emis_initial_mesh
      ray_count_initial_mesh(2)=nfreq

      ray_dims(1)=neltmax
      ray_dims(2)=nfreq

      ray_dims_initial_mesh(1)=neltmax_emis_initial_mesh
      ray_dims_initial_mesh(2)=nfreq
   
c-----------------------------------------------------------------------
c     read numbers of non-refindable mesh points
c     for all frequencies
c     nrayelt_emis_initial_mesh_nc(ifreq), ifreq=1,nfreq
c-----------------------------------------------------------------
      call ncvid2(vid,ncid,'nrayelt_emis_initial_mesh_nc',istatus)
      call ncvgt_int2(ncid,vid,1,nfreq, nrayelt_emis_initial_mesh_nc,
     +                istatus)

      do ifreq=1,nfreq
c        write(*,*)'ifreq,nrayelt_emis_initial_mesh_nc(ifreq)',
c     &             ifreq,nrayelt_emis_initial_mesh_nc(ifreq)
      enddo
c-----------------------------------------------------------------------
c     read the major radus initial (non-refindable) mesh
c     along rays for all frequencies 
c     wr_emis_initial_mesh_nc(i,2,ifreq)
c     Here
c     i=1,nrayelt_emis_initial_mesh_nc(ifreq)
c     ifreq=1,nfreq
c----------------------------------------------------------------------
      call ncvid2(vid,ncid,'wr_emis_initial_mesh_nc',istatus)
      call check_err(istatus)  
      call ncvgt_doubl2(ncid,vid,start, ray_count_initial_mesh,tem1,
     +                  istatus)
      call check_err(istatus)  
      
      j=0
      do ifreq=1,nfreq       
         do i=1,neltmax_emis_initial_mesh
            j=j+1            
            wr_emis_initial_mesh_nc(i,ifreq)=tem1(j)
c            write(*,*)'i,wr_emis_initial_mesh_nc(i,ifreq)',
c     &                 i,wr_emis_initial_mesh_nc(i,ifreq)
         enddo      
      enddo      
c-----------------------------------------------------------------------
c     read p_perpmax/me/clight at the initial (non-refindable) mesh
c     along rays for all frequencies for the second harmonic
c     wp_perpmax_dmc_nc(i,2,ifreq)
c     Here
c     i=1,nrayelt_emis_initial_mesh_nc(ifreq)
c     ifreq=1,nfreq
c----------------------------------------------------------------------
      call ncvid2(vid,ncid,'wp_perpmax_dmc_nc',istatus)
      call check_err(istatus)  
      call ncvgt_doubl2(ncid,vid,start, ray_count_initial_mesh,tem1,
     +                  istatus)
      call check_err(istatus)  
       
      j=0
      do ifreq=1,nfreq        
         do i=1,neltmax_emis_initial_mesh
            j=j+1
            wp_perpmax_dmc_nc(i,2,ifreq)=tem1(j)
c            write(*,*)'i,wp_perpmax_dmc_nc(i,2,ifreq)',
c     &                 i,wp_perpmax_dmc_nc(i,2,ifreq)
         enddo      
      enddo      

c-----------------------------------------------------------------------
c     read p_parmax/me/clight and p_parmin/me/cligh
c     at the initial (non-refindable) mesh
c     along rays for all frequencies for the second harmonic
c     wp_parpmax_dmc_nc(i,2,ifreq), wp_parpm_dmc_nc(i,2,ifreq)
c     Here
c     i=1,nrayelt_emis_initial_mesh_nc(ifreq)
c     ifreq=1,nfreq
c----------------------------------------------------------------------
      call ncvid2(vid,ncid,'wp_parmax_dmc_nc',istatus)
      call check_err(istatus)  
      call ncvgt_doubl2(ncid,vid,start,ray_count_initial_mesh,tem1,
     +                  istatus)
      call check_err(istatus)  

      j=0
      do ifreq=1,nfreq
         do i=1,neltmax_emis_initial_mesh
            j=j+1
            wp_parmax_dmc_nc(i,2,ifreq)=tem1(j)
c            write(*,*)'i,wp_parmax_dmc_nc(i,2,ifreq)',
c     &                 i,wp_parmax_dmc_nc(i,2,ifreq)
         enddo      
      enddo      

      call ncvid2(vid,ncid,'wp_parmin_dmc_nc',istatus)
      call check_err(istatus)  
      call ncvgt_doubl2(ncid,vid,start,ray_count_initial_mesh,tem1,
     +                  istatus)
      call check_err(istatus)  

      j=0
      do ifreq=1,nfreq
         do i=1,neltmax_emis_initial_mesh
            j=j+1
            wp_parmin_dmc_nc(i,2,ifreq)=tem1(j)
c            write(*,*)'i,wp_parmin_dmc_nc(i,2,ifreq)',
c     &                 i,wp_parmin_dmc_nc(i,2,ifreq)
         enddo      
      enddo      
      call ncclos2(ncid,istatus)
      call check_err(istatus)

c--------------------------------------
      deallocate( tem1,STAT=istat)
c---------------------------------
c      write(*,*)'********************************************'
c      write(*,*)'end of read_nc'
c      write(*,*)'********************************************'
      return
      end



     
      subroutine wrtnetcdf_emission_spectrum(netcdfnml)
c     ,nray_emission)
      implicit none
      !implicit integer (i-n), real*8 (a-h,o-z)

c     If the namelist variable i_emission=1 and i_emission_spectrum=1 than
c     this routine will write the emiision spectrum data 
c     into into existing netcdf file: netcdfnml 

c     It write the emission spectrum data for several frequencies
c     but only for one launching point and one launching direction
c     The variable iray should be fixed

c     The data along the ray will be written for all used frequencies
c     for the central ray.
     
      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'
      include 'write.i'   
      include 'emissa.i'      

c-----input
      character(*) netcdfnml ! input filename
c      integer nray_emission

c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'


c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2

      integer nfreq_dim,nelt_spectrum_dim,jx_kindim

      integer ray_dims(3),start(3),ray_counts(3)       
      integer flux_dims(2),start2(2),flux_counts(2)    

c-----locals
      integer neltmax_emis_spectrum,ifreq,error_code
      integer n1n2real, n_jxkin_nfreq     
c     Storage tem1r tem2r are used in netcdf writes. 
c      parameter (n1n2real=nrelta*jx_kin_a)
c      real*8 tem1r(n1n2real) 
c      parameter (n_jxkin_nfreq=jx_kin_a*nfreqa) 
c      real*8 tem2r(n_jxkin_nfreq)   
      real*8, pointer :: tem1r(:)      !(n1n2real) 
      real*8, pointer :: tem2r(:)      !(n_jxkin_nfreqa)

      integer istat

      real*8 ds,ds_old


      real*8 wye_nc
ctest
      integer n_elt,j

      data start/1,1,1/
      data start2/1,1/

      save

c-----allocate tem1r and tem2r
      n1n2real=nrelta*jx_kin
      allocate( tem1r(1:n1n2real),STAT=istat)
      call bcast(tem1r,0.d0,SIZE(tem1r))
      n_jxkin_nfreq=jx_kin*nfreq 
      allocate( tem2r(1:n_jxkin_nfreq),STAT=istat)
      call bcast(tem2r,0.d0,SIZE(tem2r))
c-----------------------------------------
c      write(*,*)'in wrtnetcdf_emission_spectrum netcdfnml=',netcdfnml

c      write(*,*)'i_emission,i_emission_spectrum',
c     &i_emission,i_emission_spectrum

      if (i_emission.eq.0) goto 10          ! nothing to do
      if (i_emission_spectrum.eq.0) goto 10 ! nothing to do
c......................................................................
cl    1. Initialize part
c
C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c
c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.
c.......................................................................
c      write(*,*)'before ncid=ncopn2(netcdfnml,NCWRITE,istatus)'
c      write(*,*)'netcdfnml = ',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)

c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf2(integer ncid,integer error_code)
c     p. 41 of netcdf manual
     
c      write(*,*)'before ncredf2'
      call ncredf2(ncid,error_code)
c      write(*,*)'after  ncredf2 ncid=',ncid
      call check_err(istatus)
c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c     dimension named 'lat'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c     returns dimension id.
c     p. 54 of netcdf manual
C-----------------------------------------------------------------------
c
c      1.1 define dimensions,variables and attribute
c
c.......................................................................
cl    1.1.2 define dimensions
c     integer function ncddef(ncid,dim_name,dim_siz,error_code)
c       returns dimension id.
c     p. 67 of netcdf-2 manual
c----------------------------------------------------------------------
c     Maximum number of ray elements( with additinal points) per ray
c     for all frequencies :
c----------------------------------------------------------------------
      neltmax_emis_spectrum=0

      do ifreq=1,nfreq         
         neltmax_emis_spectrum=max(neltmax_emis_spectrum,
     &                             nrayelt_emis_nc(ifreq))
c         write(*,*)'ifreq,nrayelt_emis_nc(ifreq)=',
c     &              ifreq,nrayelt_emis_nc(ifreq)
      enddo

      if (neltmax_emis_spectrum.eq.0) neltmax_emis_spectrum=1

c      write(*,*)'neltmax_emis_spectrum',neltmax_emis_spectrum

      ray_counts(1)=neltmax_emis_spectrum
      ray_counts(2)=jx_kin
      ray_counts(3)=1

      flux_counts(1)=jx_kin
      flux_counts(2)=nfreq

c      write(*,*)'ray_counts',ray_counts
      
c      write(*,*)'before neltmaxid =ncddef(ncid,neltmax_emis_spectrum)'
c      write(*,*)'neltmax_emis_spectrum=',neltmax_emis_spectrum

      nelt_spectrum_dim=ncddef2(ncid,'neltmax_emis_spectrum',
     &               neltmax_emis_spectrum,istatus)
      call check_err(istatus)
c      write(*,*)'after nelt_spectrum_dim=', nelt_spectrum_dim

      jx_kindim=ncddef2(ncid,'jx_kin',jx_kin,istatus)
      call check_err(istatus)
c      write(*,*)'after jx_kindim=', jx_kindim

      nfreq_dim=ncddef2(ncid,'nfreq_spectrum',nfreq,istatus)
      call check_err(istatus)
c      write(*,*)'after nfreq_dim=', nfreq_dim

c      nraysid=ncddef2(ncid,'nray_emission',nray_emission,istatus)

      ray_dims(1)=nelt_spectrum_dim
      ray_dims(2)=jx_kindim
      ray_dims(3)=nfreq_dim
c      write(*,*)'ray_dims',ray_dims

      flux_dims(1)=jx_kindim
      flux_dims(2)=nfreq_dim
c      write(*,*)'flux_dims',flux_dims
    
c-------------------------------------------------------------------------
c     Run Descriptive Parameters
c-------------------------------
c-----the emission input data set in the input files: genray.dat or genray.in

c-----kinetic energy grid descitption

c      write(*,*)'before =ncvdef2(ncid,jx_kin'
      vid=ncvdef2(ncid,'jx_kin',NCLONG,0,0,istatus)
c      write(*,*)'after =ncvdef2(ncid,jx_kin vid=',vid
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40,
     +            'number of points in kinetic energy grid ',istatus)
      call check_err(istatus)
c      write(*,*)'after =ncaptc2(ncid,jx_kin vid='

      vid=ncvdef2(ncid,'max_kin_energy_kev',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44,
     +          'maximal energy [KeV] at kinetic energy grid ',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'kinetic_energy_kev',NCDOUBLE,1,jx_kindim,
     + istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     +           'kinetic energy grid',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'KeV ',istatus)

c-----emission spectrum description
    
      vid=ncvdef2(ncid,'wj_emis_x_nc',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,73,
     +' emission coefficient spectrum  along the central ray '//
     +' with the additional points ' ,istatus)     
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           '(erg*sec/cm**3)',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ncid,w_specific_intensity_x_nc'
      vid=ncvdef2(ncid,'w_specific_intensity_x_nc',NCDOUBLE,3,ray_dims,
     &           istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,107,
     +' spectrum of the emission specific intensity I_n_0/del(s) '//
     +' from the central ray with the additional points ' ,istatus)     
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimentional',istatus)
      call check_err(istatus)

c      write(*,*)'after ncvdef2(ncid,w_specific_intensity_x_nc'


c      write(*,*)'flux_dims',flux_dims

      vid=ncvdef2(ncid,'wi_0_x_nc',NCDOUBLE,2,flux_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,97,
     +' emission flux spectrum I_0 at plasma boundary  '//
     +' from the central ray with the additional points ' ,istatus)     
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimentional',istatus)
      call check_err(istatus)



c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual 

      call ncendf2(ncid,istatus)
      call check_err(istatus)

c     End initialize

c.......................................................................
cl    1. Writing data
c
 
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c------------------------------------------------------------------------
      call ncvid2(vid,ncid,'jx_kin',istatus)
      call ncvpt_int2(ncid,vid,1,1,jx_kin,istatus)

      call ncvid2(vid,ncid,'max_kin_energy_kev',istatus)
      call ncvpt_int2(ncid,vid,1,1,max_kin_energy_kev,istatus)

      call ncvid2(vid,ncid,'kinetic_energy_kev',istatus)  
      call ncvpt_doubl2(ncid,vid,1,jx_kin,kinetic_energy_kev,istatus)
      call check_err(istatus)

c-----emissivity spectrum
      do ifreq=1,nfreq
ctest
c         write(*,*)'ifreq',ifreq
c         do n_elt=1,neltmax_emis_spectrum        !rayelt_emis_nc(ifreq)) 
c            do j=1,jx_kin
c               write(*,*)'n_elt,j,wj_emis_x_nc(n_elt,j,ifreq)',
c     &                    n_elt,j,wj_emis_x_nc(n_elt,j,ifreq)
c            enddo
c         enddo
cendtest
c         write(*,*)' neltmax_emis_spectrum',neltmax_emis_spectrum
         call pack21(wj_emis_x_nc(1,1,ifreq),1,nrelta,1,jx_kin,tem1r,
     &               neltmax_emis_spectrum,jx_kin)
c         write(*,*)'tem1r',tem1r         
         call ncvid2(vid,ncid,'wj_emis_x_nc',istatus)
         start(3)=ifreq
         call ncvpt_doubl2(ncid,vid,start,ray_counts,tem1r,istatus)
      enddo


c------------spectrum of the specific intensity
       do ifreq=1,nfreq   
c          write(*,*)'wrtnetcdf_emission_spectrum ifreq=',ifreq

          do n_elt=1,nrayelt_emis_nc(ifreq)-1 
            ds=wsn_nc(n_elt+1,ifreq)-wsn_nc(n_elt,ifreq)
            do j=1,jx_kin
               if(ds.gt.1.d-18) then
                 w_specific_intensity_x_nc(n_elt,j,ifreq)=
     &           win_0_x_nc(n_elt,j,ifreq)/ds
               else
                 w_specific_intensity_x_nc(n_elt,j,ifreq)=
     &           win_0_x_nc(n_elt,j,ifreq)/ds_old
               endif
            enddo
   
            ds_old=ds
          enddo

ctest
c         write(*,*)'ifreq',ifreq
c         do n_elt=1,neltmax_emis_spectrum        !rayelt_emis_nc(ifreq)) 
c          do j=1,jx_kin
c           write(*,*)'n_elt,j,w_specific_intensity_x_nc(n_elt,j,ifreq)'
c     &               ,n_elt,j,w_specific_intensity_x_nc(n_elt,j,ifreq)
c            enddo
c         enddo
cendtest

c         write(*,*)' neltmax_emis_spectrum',neltmax_emis_spectrum
         call pack21(w_specific_intensity_x_nc(1,1,ifreq),1,nrelta,
     &   1,jx_kin,tem1r,neltmax_emis_spectrum,jx_kin)
c         write(*,*)'tem1r',tem1r        
         call ncvid2(vid,ncid,'w_specific_intensity_x_nc',istatus)
         call check_err(istatus)
c         write(*,*)'before call ncvpt(ncid,vid,start,ray_counts,tem1r'
         start(3)=ifreq
         call ncvpt_doubl2(ncid,vid,start,ray_counts,tem1r,istatus)
         call check_err(istatus)
      enddo
      

c--------------i_0_x-----------------
c      do ifreq=1,nfreq
c         do j=1,jx_kin
c           write(*,*)'ifreq,j,wi_0_x_nc(j,ifreq)',
c     &                ifreq,j,wi_0_x_nc(j,ifreq)
c         enddo
c      enddo

      call pack21(wi_0_x_nc,1,jx_kin,1,nfreq,tem2r,jx_kin,nfreq)
c      write(*,*)'tem2r',tem2r
      call ncvid2(vid,ncid,'wi_0_x_nc',istatus)
      call ncvpt_doubl2(ncid,vid,start2,flux_counts,tem2r,istatus)  
 

c     4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

c     End write netcdf data

 10   continue

      deallocate( tem1r,STAT=istat)
      deallocate( tem2r,STAT=istat)

      return
      end


      subroutine plasma_profiles_vs_r(z)
c------------------------------------------------------
c     Creates 1D arrays of density,temperature,zeff
c     versus major radius r at given z
c-------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i' 
      include 'onetwo_nml.i'
      include 'five.i'
      include 'writencdf.i'
      
c-----input
      double precision z    ! vertical coordinate         
c-----externals
      double precision psi_rho,rhopsi,psif,zeffrho,densrho,temperho
c-----locals
      double precision 
     &accuracy,  ! accuracy of the solutions
                 ! r=r_min_z ans r=r_max_z
                 ! of the equation psi_limitter-psi(r,z)=0   
                 ! at the given z
     &psi_loc,h_r,r_loc,rho_loc,z_loc,rmin_psi_z,rmax_psi_z

      integer i,k

c     &nr                    ! number of points
c                            ! of profiles mesh is in onetwo_nml.i
cSAP090710
c      nr=NRA

      accuracy=1.d-4
      psi_loc=psi_rho(1.d0)
      z_loc=z
c      write(*,*)'psi_loc,z_loc,accuracy',psi_loc,z_loc,accuracy
c      write(*,*)'rmax,rmin',rmax,rmin
      call r_min_r_max_from_psi_z(psi_loc,z_loc,accuracy,
     &rmin_psi_z,rmax_psi_z)

      h_r=(rmax-rmin)/dfloat(nr-1)
      do i=1,nr
         r_loc=rmin+h_r*(i-1)
c         write(*,*)'i,rmin,r_loc',i,rmin,r_loc
         w_r_densprof_nc(i)=r_loc
         psi_loc=psif(z_loc,r_loc)
         rho_loc=rhopsi(psi_loc)
c         write(*,*)'i,z_loc,rho,psi_loc',i,z_loc,rho_loc,psi_loc 
         do k=1,nbulk    
           w_dens_vs_r_nc(i,k)=densrho(rho_loc,k)*1.d19
           !YuP[07-2017] corrected 2nd index. Was '1', now 'k'.
           !YuP[07-2017] w_dens_vs_r_nc() is saved in units ptcls/m^3
           ! densrho is in  10^13 ptcls/cm^3  == 10^19 ptcls/m^3
           !So, in the above, added factor 1.d19
           w_temp_vs_r_nc(i,k)=temperho(rho_loc,k) ! keV
c           write(*,*)'k,R_loc,Z_loc,w_dens_vs_r_nc(irad,k)=',
c     +                k,R_loc,Z_loc,w_dens_vs_r_nc(i,k)
         enddo
         w_zeff_vs_r_nc(i)=zeffrho(rho_loc)
      enddo
      
      return
      end

                           
      subroutine wrtnetcdf_one_ray_point(kopt,is0,nray0)            
      implicit none
c
c     Write ray tracing data into a netCDF file genray_one_ray_point.nc
c     in one ray point:
c     at the ray with number nray0 
c     at the pointwith number is0.

cSm0940727
c     If the parameter ionetwo.eq.1 it will write the current
c     and power profiles into a netcdf file

      include 'param.i'
      include 'writencdf.i'
      include 'one.i'
      include 'ions.i'
      include 'adj.i'
      include 'cone_nml.i'     !nccone
      include 'grill_nml.i'    !ngrill
cSAP090314
c      include 'writencdf_one_ray_point.i'
c--------------------------
cSm040727 to write the current and power profiles into netcdf file
c     Done in subroutine, wrtnetcdf_prof
c      include 'onetwo.i'
c--------------------------
c-----input
      integer kopt
      integer nray0, !the number of ray
     &is0            !the number of ray point 
 
c-----locals
      character filenc_one_ray_point* 128
cSAP090903
      integer n1n2,istat
c     Storage tem1 is used in netcdf writes, including complex numbers.
c      parameter (n1n2=2*nrelta*nraya)
c      real*8 tem1(n1n2) 
      real*8, pointer :: tem1(:)      !(n1n2)
 
cSAP090314
c     Storage to put array(nraya) to one point array(nraya1=1)
c     at one ray 
      integer nrelta1,nraya1
      parameter(nrelta1=1)
      parameter(nraya1=1)
      real*8 point_2D(nrelta1,nraya1),point_3D(nrelta1,nraya1,nbulka)
      integer  ipoint_1D(nraya1)
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,nccre2,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid,twoid,char64id,char128id,char256id,
cSAP
     &char8id

    
      integer nbulkmid  !ion species dimension
      integer nbulkid
      integer nrho_adjid

      integer ray_dims(3),start(3),ray_count(3)
      integer ray_dimss(3),starts(3),ray_counts(3)

      double complex cei

      character ltitle*256

      integer neltmax,iray,i,j,ii,ll,nrayl0
      integer length_char

      data start/1,1,1/,starts/1,1,1/


      save

c      write(*,*)'@@@ wrtnetcdf_one_ray_point begin'
      if(is0.gt.nrayelt_nc(nray0)) then
         write(*,*)'WARNING is0.gt.nrayelt_nc(nray0)'
         write(*,*)'is0,nray0,nrayelt_nc(nray0)',
     &              is0,nray0,nrayelt_nc(nray0)
         write(*,*)'the code will not create genray_one_ray_point.nc'
         return
      endif

c----------------------------------------------
c     allocate tem1
c----------------------------------------------
      n1n2=2*nrelta1*nraya1
      allocate( tem1(1:n1n2),STAT=istat)
      call bcast(tem1,0.d0,SIZE(tem1))
c----------------------------------------------

cSAP090309
      filenc_one_ray_point='genray_one_ray_point.nc'

      cei=(0.d0,1.d0)

c     Maximum number of ray elements per ray:
      neltmax=0
cSAP090309 for one point of the ray
      neltmax=1
cSAP090314 number of rays nrayl0
      nrayl0=1  
c      write(*,*)'wrtnetcdf nray=',nrayl

cSm05038
      if (neltmax.eq.0) neltmax=1

      ray_count(1)=neltmax
cSAP090314
c     ray_count(2)=nrayl
      ray_count(2)=nrayl0
      ray_count(3)=2
      ray_counts(1)=neltmax
cSAP090314
c      ray_counts(2)=nrayl
      ray_counts(2)=nrayl0
      ray_counts(3)=1

c      write(*,*)'ray_count',ray_count

c.......................................................................
cl    1. Initialize part, creating new netcdf file
c

c --- begin if ---
      if ( kopt.eq.1 ) then
c      write(*,*)'wrtnetcdf,  kopt=',kopt

C-----------------------------------------------------------------------
c
cl     1.1 create netCDF file and define dimensions,variables
c          and attributes
c

c.......................................................................
cl    1.1.1 create netCDF filename (Entering define mode.)
c     integer function nccre(filename,overwrite?,error_code)
c     Ref to page 46 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

c      write(*,*)'wrtnetcdf before 
c     &ncid=nccre(filenc_one_ray_point,NCCLOB,istatus)'

      ncid=nccre2(filenc_one_ray_point,NCCLOB,istatus)
      call check_err(istatus)
c      write(*,*)'ncid',ncid      
c     Brief description added to file:
      ltitle='netCDF file of ray data from GENRAY version: '//version
      if( length_char(trim(ltitle)).gt.256 ) 
     + stop 'Adjust ltitle in netcdfrw2'
      call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,
     +     length_char(trim(ltitle)), trim(ltitle), istatus)
        
c.......................................................................
cl    1.1.2 define dimensions
c     integer function ncddef(ncid,dim_name,dim_siz,error_code)
c       returns dimension id.
c     p. 67 of netcdf-2 manual

c     For ray data:
c      write(*,*)'wrtnetcdf before ncddef(neltmax) neltmax=',neltmax

      neltmaxid=ncddef2(ncid,'neltmax',neltmax,istatus)         
     
c      write(*,*)'wrtnetcdf before ncddef(ncid,nrays) nrayl0=',nrayl0
cSAP090314
c      nraysid=ncddef2(ncid,'nrays',nrayl,istatus)
      nraysid=ncddef2(ncid,'nrays',nrayl0,istatus)

c      write(*,*)'wrtnetcdf before ncddef(ncid,two,2,istatus)'
      
      twoid=ncddef2(ncid,'two',2,istatus)
    
c      write(*,*)'wrtnetcdf before ncddef(ncid,nbulk)'

      nbulkid=ncddef2(ncid,'nbulk',nbulk,istatus) 

cSAP080303
c      write(*,*)'wrtnetcdf before ncddef(char8dim)'
      char8id=ncddef2(ncid,'char8dim',8,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char64dim)'

      char64id=ncddef2(ncid,'char64dim',64,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char128dim)'

      char128id=ncddef2(ncid,'char128dim',128,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char5212dim)'

      char256id=ncddef2(ncid,'char256dim',256,istatus)

c      write(*,*)'neltmaxid',neltmaxid
c      write(*,*)'nraysid',nraysid
c      write(*,*)'twoid',twoid

      ray_dims(1)=neltmaxid
      ray_dims(2)=nraysid
      ray_dims(3)=twoid
    
      if (nbulk.gt.1) then

c         write(*,*)'wrtnetcdf before ncdef(nbulkm,nbulk-1 )'

         nbulkmid=ncddef2(ncid,'nbulkm',nbulk-1,istatus)

         ray_dimss(1)=ray_dims(1)
         ray_dimss(2)=ray_dims(2)
         ray_dimss(3)=nbulkmid
      endif
      
c.......................................................................
cl    1.1.3 define variables

c     Note, the variable IDs (denoted below as vid) are
c     not saved here in this subroutine; rather, the IDs
c     are retrieved from the netCDF data file, as needed,
c     by calling the netCDF routine ncvid.

c     netCDF variable_define:
c     integer function ncvdef2(ncid,variable_name,variable_type,
c                number_of_dimensions,
c                vector_for_length_per_dimension,error_code)
c       returns varid.
c     Note: Unlimited dimension must be last
c     Refer to p. 77, netcdf-2 manual

c     NCDOUBLE for REAL*8.  This is independent of the
c       internal representation of data on a particular
c       machine, e.g., 32- or 64-bit based calculations.


c--------------------------
c     Genray version
c--------------------------

c      write(*,*)'before ncvdef2(version)'
c      write(*,*)'version',version
      vid=ncvdef2(ncid,'version',NCCHAR,1,char64id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +          'GENRAY version number',istatus)

c--------------------------
c     Mnemonic for the run
c--------------------------
c      write(*,*)'before ncvdef2(mnemonic)'
      vid=ncvdef2(ncid,'mnemonic',NCCHAR,1,char128id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +          'Mnemonic run identifier',istatus)

c-------------------------------
c     Run Descriptive Parameters
c-------------------------------
c      write(*,*)'before ncvdef2(vid)'
      vid=ncvdef2(ncid,'id',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +            'Disp relation identifier',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(iabsorp)'
      vid=ncvdef2(ncid,'iabsorp',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Absorp calc identifier',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ieffic)'
      vid=ncvdef2(ncid,'ieffic',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29,
     +            'Current drive calc identifier',istatus)
      call check_err(istatus)

cSAP080303
c      write(*,*)'before ncvdef2(ion_absorption)'
      vid=ncvdef2(ncid,'ion_absorption',NCCHAR,1,char8id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +'Switch on/off ion absorption at iabsorp=3,9,91,92',istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(refl_loss)'
      vid=ncvdef2(ncid,'refl_loss',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +    'fraction of power loss at each reflection',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(iflux)'
      vid=ncvdef2(ncid,'iflux',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Flux calc, non-Westerhof-Tokman id',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ioxm)'
      vid=ncvdef2(ncid,'ioxm',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     +            'wave mode indicator (1 - om, -1 - xm )',istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(ioxm_n_npar)'
      vid=ncvdef2(ncid,'ioxm_n_npar',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,68,
     +'wave mode indicator: sign before square root to find '//
     +' N(N_parallel)',
     +istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(jwave)'
      vid=ncvdef2(ncid,'jwave',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37,
     +            'Wave harmonic, for CD efficiency calc',istatus)
      call check_err(istatus)

      if (iabsorp.eq.4) then
      vid=ncvdef2(ncid,'i_im_nperp',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     +       'iabsorp=4  nperp:1,ImD_full/dD/dnp; 2,Cmplx soln',istatus)
      call check_err(istatus)
      endif

c      write(*,*)'before ncvdef2(i_geom_optic)'
      vid=ncvdef2(ncid,'i_geom_optic',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Integrate rays wrt (1)time,(2)dist',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(istart)'
      vid=ncvdef2(ncid,'istart',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,50,
     +     'Ray launch type: 1,eccone; 2,grill, 3,OX in plasma',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ncone)'
      vid=ncvdef2(ncid,'ncone',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,39,
     +     'Number of rf source ray cones, istart=1',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,40,
     +     'Each cone has (nray/ncone) launched rays, istart=1',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ngrill)'
      vid=ncvdef2(ncid,'ngrill',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,41,
     +     'Number of rf source ray grills, istart=2,3',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,47,
     +     ',nray/ngrill rays launched per grill istart=2,3',istatus)
      call ncaptc2(ncid,vid,'long_name2',NCCHAR,55,
     +     ',Each grill has (nray/ngrill) launched rays, istart=2,3',
     +     istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ionetwo)'
      vid=ncvdef2(ncid,'ionetwo',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +            'if ionetwo=1 then calculate CD',istatus)
      call check_err(istatus)
c--------------------------
c     Ray data
c--------------------------
c      write(*,*)'before ncvdef2(nray)'
      vid=ncvdef2(ncid,'nray',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +            'Number of rays',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(nharm)'
      vid=ncvdef2(ncid,'nharm',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'First harmonic number',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(freqncy)'
      vid=ncvdef2(ncid,'freqcy',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +            'Wave frequency',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'Hz',istatus)

c      write(*,*)'before ncvdef2(i_emission)'
      vid=ncvdef2(ncid,'i_emission',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,63,
     +'emission switch: =0 now emission, =1 use emission calculations',
     +istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(i_emission_spectrum)'
      vid=ncvdef2(ncid,'i_emission_spectrum',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,91,
     +'emission spectrum switch: =0 no emission spectrum, '//
     +'  =1 use emission spectrum  calculations',
     +istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(nrayelt)'

      vid=ncvdef2(ncid,'nrayelt',NCLONG,1,nraysid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     +            'Number of ray elements for each ray',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ws)'

      vid=ncvdef2(ncid,'ws',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29,
     +           'poloidal distance along a ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

c      write(*,*)'before ncvdef2(seikon)'

      vid=ncvdef2(ncid,'seikon',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,7,
     +           'eikonal',istatus)


c      write(*,*)'before ncvdef2(spsi)'

      vid=ncvdef2(ncid,'spsi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,45,
     +         'normalized small radius=rho given by indexrho',istatus)

c      write(*,*)'before ncvdef2(wr)'

      vid=ncvdef2(ncid,'wr',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)


c      write(*,*)'before ncvdef2(wphi)'

      vid=ncvdef2(ncid,'wphi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'toroidal angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'radians',istatus)

c      write(*,*)'before ncvdef2(wz)'

      vid=ncvdef2(ncid,'wz',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,15,
     +           'vertical height',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

c      write(*,*)'before ncvdef2(w_theta_pol)'

      vid=ncvdef2(ncid,'w_theta_pol',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'poloidal angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'degree',istatus)

c      write(*,*)'before ncvdef2(wnpar)'

      vid=ncvdef2(ncid,'wnpar',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'parallel refractive index',istatus)

c      write(*,*)'before ncvdef2(wnper)'

      vid=ncvdef2(ncid,'wnper',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'perpendicular refractive index',istatus)

c      write(*,*)'before ncvdef2(delpwr)'

      vid=ncvdef2(ncid,'delpwr',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'power in ray channel',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)

c      write(*,*)'before ncvdef2(sdpwr)'


      vid=ncvdef2(ncid,'sdpwr',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,48,
     +      'Ion collisionless absorption coeff (all species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

cSm050411      if (iabsorp.eq.3.and.nbulk.gt.1) then
cSm061127      if ((iabsorp.eq.3.or.iabsorp.eq.9).and.nbulk.gt.1) then
      if (((iabsorp.eq.3.or.iabsorp.eq.9).or.
     &    (iabsorp.eq.91.or.iabsorp.eq.92))
     &     .and.nbulk.gt.1) then
      vid=ncvdef2(ncid,'salphas',NCDOUBLE,3,ray_dimss,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +     'Ion collisionless absorption coeff (each species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)
      endif

c      write(*,*)'before ncvdef2(wdnpar)'

      vid=ncvdef2(ncid,'wdnpar',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,6,
     +           'wdnpar',istatus)

c      write(*,*)'before ncvdef2(cwexde)'

c     Added 3rd dimension equal to 2 accomodates complex data.
c      write(*,*)'ncid=',ncid
c      write(*,*)'ray_dims',ray_dims
      vid=ncvdef2(ncid,'cwexde',NCDOUBLE,3,ray_dims,istatus)
c      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ex/E Polarization',istatus)

c      write(*,*)'before ncvdef2(cweyde)'

      vid=ncvdef2(ncid,'cweyde',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ey/E Polarization',istatus)

c      write(*,*)'before ncvdef2(cwezde)'

      vid=ncvdef2(ncid,'cwezde',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ez/E Polarization',istatus)

c      write(*,*)'before ncvdef2(fluxn)'

      vid=ncvdef2(ncid,'fluxn',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +           'fluxn, Stix norm, |E|=1',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13,
     +           'ergs/sec/cm^2',istatus)

c      write(*,*)'before ncvdef2(sbtot)'

      vid=ncvdef2(ncid,'sbtot',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +           'Magnetic field strength',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'guass',istatus)

c      write(*,*)'before ncvdef2(cene)'

      vid=ncvdef2(ncid,'sene',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,17,
     +           'Density along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14,
     +           'particles/cm^3',istatus)
     
c      write(*,*)'before ncvdef2(ste)'

      vid=ncvdef2(ncid,'ste',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +           'Temperature along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)

      vid=ncvdef2(ncid,'szeff',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'Zeff along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13,
     +           'dimensionless',istatus)

      vid=ncvdef2(ncid,'salphac',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'Collisional damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

      vid=ncvdef2(ncid,'salphal',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Linear damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

c      write(*,*)'before ncvdef2(sb_r)'

      vid=ncvdef2(ncid,'sb_r',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_r magnetic field component',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'guass',istatus)

c      write(*,*)'before ncvdef2(sb_z)'

      vid=ncvdef2(ncid,'sb_z',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_z magnetic field component',istatus) 
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'guass',istatus)

c      write(*,*)'before ncvdef2(sb_phi)'

      vid=ncvdef2(ncid,'sb_phi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'B_phi magnetic field component',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'guass',istatus)

c      write(*,*)'before ncvdef2(wn_r)'

      vid=ncvdef2(ncid,'wn_r',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_r refractive index component',istatus)

c      write(*,*)'before ncvdef2(wn_z)'

      vid=ncvdef2(ncid,'wn_z',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_z refractive index component',istatus) 

c      write(*,*)'before ncvdef2(wn_phi)'

      vid=ncvdef2(ncid,'wn_phi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32,
     +           'N_phi refractive index component',istatus)
     
c      write(*,*)'before ncvdef2(wgr_r)'

      vid=ncvdef2(ncid,'vgr_r',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_r normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(wgr_z)'

      vid=ncvdef2(ncid,'vgr_z',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_z normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(wgr_hpi)'

      vid=ncvdef2(ncid,'vgr_phi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +           'vgroup_phi normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(flux_z)'

       
      vid=ncvdef2(ncid,'flux_z',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +           'flux_z normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(flux_r)'

       
      vid=ncvdef2(ncid,'flux_r',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +           'flux_r normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(flux_phi)'

       
      vid=ncvdef2(ncid,'flux_phi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'flux_phi normalized to c',istatus)

c--------------------------------------------------------
c      CD efficiency along rays
c--------------------------------------------------------
      if (ionetwo.eq.1)then
c         write(*,*)'before ncvdef2(w_eff_nc)'
         vid=ncvdef2(ncid,'w_eff_nc',NCDOUBLE,2,ray_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'CD efficiency along a ray',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,27,
     +           '(A/cm**2)/(erg/(sec*cm**3))',istatus)
       endif
c--------------------------------------------------------
c      DC electric field conductivity from adj 
c--------------------------------------------------------
       if (i_adj.eq.1)then

c         write(*,*)'before ncvdef2(i_adj)'
         vid=ncvdef2(ncid,'i_adj',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     &   '=0 no adj; =1 use adj calculations ',istatus)
         call check_err(istatus)

c         write(*,*)'before ncvdef2 npsi0)'
         vid=ncvdef2(ncid,'npsi0',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,53,
     &   'Number of flux surfaces where ADJ equation was solved',
     &   istatus)
         call check_err(istatus)


         nrho_adjid=ncddef2(ncid,'nrho_adj',npsi0,istatus)
         call check_err(istatus)

c         write(*,*)'sigma_E_adj',sigma_E_adj
c         write(*,*)'rho_adj',rho_adj
c         write(*,*)'before ncvdef2(sigma_E_adj)'
         vid=ncvdef2(ncid,'sigma_E_adj',NCDOUBLE,1,nrho_adjid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,74,
     &   'DC electric field conductivity at radii where '//
     &   'adj function is calculated ',istatus) 
         call ncaptc2(ncid,vid,'units',NCCHAR,19,
     +           'sigma/sigma_Spitzer',istatus)

c         write(*,*)'before ncvdef2(rho_adj)'
         vid=ncvdef2(ncid,'rho_adj',NCDOUBLE,1,nrho_adjid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,50,
     &   'small radii mesh where adj function is calculated ',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)

       endif
c--------------------------
c     determine OX conversion data for i_ox=2 case
c--------------------------
      goto 10
      if (i_ox.eq.2) then

         vid=ncvdef2(ncid,'i_ox_conversion',NCLONG,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     &            'Option equals 1 after OX conversion',istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'transm_ox',NCDOUBLE,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     &   'OX transmission coefficient Preinhaelter and Kopecky 1973',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cn_par_optimal',NCDOUBLE,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     &              'optimal N parallel for OX transmission',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cnpar_ox',NCDOUBLE,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,33,
     &              'N parallel before OX transmission',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cn_b_gradpsi',NCDOUBLE,1,nraysid,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,45,
     &              'N along [B^*grad(psi)] before OX transmission',
     &              istatus)
         call check_err(istatus)
         
      endif !i_ox=2
 10   continue
c--------------------------
c     Plasma data
c--------------------------
     
c      write(*,*)'before ncvdef2(eqdskin)'

       
      vid=ncvdef2(ncid,'eqdskin',NCCHAR,1,char256id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39,
     +          'Name of input eqdsk, for eqsource=eqdsk',istatus)
     
c      write(*,*)'before ncvdef2(nbulk)'

       
      vid=ncvdef2(ncid,'nbulk',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,43,
     +          'Number of Maxl plasma cmpts, electrons+ions',istatus)
      call check_err(istatus)     

c      write(*,*)'before ncvdef2(dmas)'

       

      vid=ncvdef2(ncid,'dmas',NCDOUBLE,1,nbulkid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +           'plasma species mass: electrons, then ions',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,27,
     +           'Normalized to electron mass',istatus)
     
c      write(*,*)'before ncvdef2(charge)'


      vid=ncvdef2(ncid,'charge',NCDOUBLE,1,nbulkid,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44,
     +           'plasma species charge: electrons, then ions',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,31,
     +           'Normalized to electronic charge',istatus)

c------------------------------------------------------------------
c     for total power [erg/sec] absorbed at all reflections at all rays
c------------------------------------------------------------------
c      write(*,*)'before ncvdef2(w_tot_pow_absorb_at_refl_nc)'
      vid=ncvdef2(ncid,'w_tot_pow_absorb_at_refl_nc',NCDOUBLE,
     &0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,47,
     +    'Total power absorbed at reflections of all rays',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
c.................................
cl    1.1.4 end the define-mode and start the data-mode
c     p. 51-2 of manual
c      write(*,*)'ncid',ncid
      call ncendf2(ncid,istatus)
      call check_err(istatus)

c      write(*,*)'end initialization'

      endif               ! End initialize




c.......................................................................
cl    1. Writing data
c

      if ( kopt.eq.0 ) then
c      write(*,*)'wrtnetcdf,  kopt=',kopt
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c--------------------------
c      write(*,*)'one_ray_point version',version
c      write(*,*)'ncid',ncid
      call ncvid2(vid,ncid,'version',istatus)
c      write(*,*)'vid',vid 
      ll=length_char(version)
c      write(*,*)'ll',ll
      call ncvptc2(ncid,vid,1,ll,version,ll,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'mnemonic',istatus)
      ll=length_char(trim(mnemonic))
      call ncvptc2(ncid,vid,1,ll,trim(mnemonic),ll,istatus)

      call ncvid2(vid,ncid,'eqdskin',istatus)
      ll=length_char(trim(eqdskin))
      call ncvptc2(ncid,vid,1,ll,trim(eqdskin),ll,istatus)

      call ncvid2(vid,ncid,'dmas',istatus)
      call ncvpt_doubl2(ncid,vid,1,nbulk,dmas,istatus)

      call ncvid2(vid,ncid,'charge',istatus)
      call ncvpt_doubl2(ncid,vid,1,nbulk,charge,istatus)

      call ncvid2(vid,ncid,'istart',istatus)
      call ncvpt_int2(ncid,vid,1,1,istart,istatus)

      call ncvid2(vid,ncid,'ncone',istatus)
      call ncvpt_int2(ncid,vid,1,1,ncone,istatus)

      call ncvid2(vid,ncid,'ngrill',istatus)
      call ncvpt_int2(ncid,vid,1,1,ngrill,istatus)

      call ncvid2(vid,ncid,'ionetwo',istatus)
      call ncvpt_int2(ncid,vid,1,1,ionetwo,istatus)

c--------------------------
c     Run specs
c--------------------------  
      call ncvid2(vid,ncid,'id',istatus)
      call ncvpt_int2(ncid,vid,1,1,id,istatus)

      call ncvid2(vid,ncid,'refl_loss',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,refl_loss,istatus)

      call ncvid2(vid,ncid,'iabsorp',istatus)
      call ncvpt_int2(ncid,vid,1,1,iabsorp,istatus)

      call ncvid2(vid,ncid,'ieffic',istatus)
      call ncvpt_int2(ncid,vid,1,1,ieffic,istatus)

      call ncvid2(vid,ncid,'ion_absorption',istatus)
      ll=length_char(ion_absorption)
      call ncvptc2(ncid,vid,1,ll,ion_absorption,ll,istatus)

      call ncvid2(vid,ncid,'iflux',istatus)
      call ncvpt_int2(ncid,vid,1,1,iflux,istatus)

      call ncvid2(vid,ncid,'ioxm',istatus)
      call ncvpt_int2(ncid,vid,1,1,ioxm,istatus)

      call ncvid2(vid,ncid,'ioxm_n_npar',istatus)
      call ncvpt_int2(ncid,vid,1,1,ioxm_n_npar,istatus)

      call ncvid2(vid,ncid,'jwave',istatus)
      call ncvpt_int2(ncid,vid,1,1,jwave,istatus)

      call ncvid2(vid,ncid,'i_geom_optic',istatus)
      call ncvpt_int2(ncid,vid,1,1,i_geom_optic,istatus)

      
      if (iabsorp.eq.4) then
         call ncvid2(vid,ncid,'i_im_nperp',istatus)
         call ncvpt_int2(ncid,vid,1,1,i_im_nperp,istatus)
      endif

      call ncvid2(vid,ncid,'nbulk',istatus)
      call ncvpt_int2(ncid,vid,1,1,nbulk,istatus)

c--------------------------
c     Ray data
c--------------------------
      call ncvid2(vid,ncid,'nray',istatus)
      call ncvpt_int2(ncid,vid,1,1,nrayl0,istatus)

      call ncvid2(vid,ncid,'nharm',istatus)
      call ncvpt_int2(ncid,vid,1,1,nharm,istatus)

      call ncvid2(vid,ncid,'freqcy',istatus)
      call ncvpt_doubl2(ncid,vid,1,1,freqcy,istatus)
      
      call ncvid2(vid,ncid,'i_emission',istatus)

      call ncvpt_int2(ncid,vid,1,1,i_emission,istatus)      

      
      call ncvid2(vid,ncid,'i_emission_spectrum',istatus)
      call ncvpt_int2(ncid,vid,1,1,i_emission_spectrum,istatus)
      
      call ncvid2(vid,ncid,'nrayelt',istatus)
cSAP090314
c      call ncvpt_int2(ncid,vid,1,nrayl,nrayelt_nc,istatus)
      ipoint_1D(1)=1
      call ncvpt_int2(ncid,vid,1,nrayl0,ipoint_1D,istatus)

c      write(*,*)'in wrtnetcdf nrayl,neltmax',nrayl,neltmax
c      do i=1,nrayl
c         write(*,*)'number of rayi,nrayelt_nc(i)',i,nrayelt_nc(i)
c         do j=1,nrayelt_nc(i)
c            write(*,*)'j,ws_nc(j,i)',j,ws_nc(j,i)
c         enddo
c      enddo

cSAP090314
      point_2D(1,1)=ws_nc(is0,nray0)

cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
       call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
c      do i=1,nrayl
c         write(*,*)'wrtnetcdf tem1 i=',i
c         do j=1,neltmax          
c            k=(i-1)*neltmax+j
c            write(*,*)'j,i,k tem1(k)',j,i,k,tem1(k)
c         enddo
c      enddo
c      write(*,*)'ws: vid,start,ray_count',vid,start,ray_count
      call ncvid2(vid,ncid,'ws',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
c      write(*,*)'is0,nray0',is0,nray0

      point_2D(1,1)=seikon_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'seikon',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=spsi_nc(is0,nray0)
      
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'spsi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wr_nc(is0,nray0)
      
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314      
      point_2D(1,1)=wphi_nc(is0,nray0)
cSAP090903  
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
       call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wphi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wz_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wz',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=w_theta_pol_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'w_theta_pol',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wnpar_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wnpar',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wnper_nc(is0,nray0)
      
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
       call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wnper',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=delpwr_nc(is0,nray0)
cSAP090903  
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'delpwr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=sdpwr_nc(is0,nray0)
cSAP090903       
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sdpwr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)



cSm050411      if (iabsorp.eq.3) then
cSm061127      if (iabsorp.eq.3.or.iabsorp.eq.9) then
      if ((iabsorp.eq.3.or.iabsorp.eq.9).or.
     &    (iabsorp.eq.91.or.iabsorp.eq.92)) then

cSAP090314
      
c      do j=1,nbulka
      do j=1,nbulk
         point_3D(1,1,j)=salphas_nc(is0,nray0,j)
      enddo


      do i=2,nbulk
cSAP090903      
c        call pack21(point_3D(1,1,i),1,nrelta,1,nraya,tem1,neltmax,
c     &  nrayl0)
c        call pack21(point_3D(1,1,i),1,nrelta,1,nrayl,tem1,neltmax,
c     &  nrayl0)
        call pack21(point_3D(1,1,i),1,1,1,1,tem1,neltmax,
     &  nrayl0)
        starts(1)=start(1)
        starts(2)=start(2)
        starts(3)=i-1  ! nbulk-1 ion species salphas are put in .nc file
        call ncvid2(vid,ncid,'salphas',istatus)
        call ncvpt_doubl2(ncid,vid,starts,ray_counts,tem1,istatus)
      enddo
      endif

cSAP090314
      point_2D(1,1)=wdnpar_nc(is0,nray0)
cSAP090903 
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wdnpar',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c-----------------------------------------------------------------
      ii=0
cSAP090314
      j=nray0
      i=is0
      ii=ii+1
      tem1(ii)=0.5d0*(cwexde_nc(i,j)+dconjg(cwexde_nc(i,j)))
      
      ii=ii+1
      tem1(ii)=-cei*0.5d0*(cwexde_nc(i,j)-dconjg(cwexde_nc(i,j)))
      call ncvid2(vid,ncid,'cwexde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c-------------------------------------------------------------------
      ii=0
cSAP090314
      j=nray0
      i=is0
      ii=ii+1
      tem1(ii)=0.5d0*(cweyde_nc(i,j)+dconjg(cweyde_nc(i,j)))

      ii=ii+1
      tem1(ii)=-cei*0.5d0*(cweyde_nc(i,j)-dconjg(cweyde_nc(i,j)))
      call ncvid2(vid,ncid,'cweyde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

c-------------------------------------------------------------
      ii=0
cSAP090314
      j=nray0
      i=is0
      ii=ii+1
      tem1(ii)=0.5d0*(cwezde_nc(i,j)+dconjg(cwezde_nc(i,j)))

      ii=ii+1
      tem1(ii)=-cei*0.5d0*(cwezde_nc(i,j)-dconjg(cwezde_nc(i,j)))
      call ncvid2(vid,ncid,'cwezde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c---------------------------------------------------------------

cSAP090314
      point_2D(1,1)=fluxn_nc(is0,nray0)
cSAP090903 
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelt,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'fluxn',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=sbtot_nc(is0,nray0)
cSAP090903 
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sbtot',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=sene_nc(is0,nray0)
cSAP090903 
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelt,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sene',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=ste_nc(is0,nray0)
cSAP090903 
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'ste',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cBH111028
      point_2D(1,1)=szeff_nc(is0,nray0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'szeff',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=salphac_nc(is0,nray0)
cSAP090903 
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'salphac',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=salphal_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'salphal',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=sb_r_nc(is0,nray0)

cSAP090903 
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sb_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=sb_z_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sb_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=sb_phi_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sb_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wn_r_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wn_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wn_z_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wn_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wn_phi_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wn_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=vgr_r_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'vgr_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=vgr_z_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'vgr_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=vgr_phi_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'vgr_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=flux_r_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'flux_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=flux_z_nc(is0,nray0)
cSAP090903
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'flux_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=flux_phi_nc(is0,nray0)
cSAP090903  
c      call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'flux_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c--------------------------------------------------------------
c      CD efficiency along rays
c--------------------------------------------------------
c       write(*,*)'netcdfr3d ionetwo',ionetwo
cSAP090314
      point_2D(1,1)=w_eff_nc(is0,nray0)

       if (ionetwo.eq.1) then            
cSAP090903  
c          call pack21(point_2D,1,nrelta,1,nraya,tem1,neltmax,nrayl0)
c          call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
          call pack21(point_2D,1,1,1,1,tem1,neltmax,nrayl0)
          call ncvid2(vid,ncid,'w_eff_nc',istatus)
          call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
       endif
c--------------------------------------------------------
c      DC electric field conductivity from adj 
c---------------------------`-----------------------------
c       write(*,*)'netcdfr3d i_adj',i_adj
       if (i_adj.eq.1)then  
c           write(*,*)'before vid=ncvid(i_adj'
           call ncvid2(vid,ncid,'i_adj',istatus)
           call ncvpt_int2(ncid,vid,1,1,i_adj,istatus)

c           write(*,*)'before vid=ncvid(npsi0i'
           call ncvid2(vid,ncid,'npsi0',istatus)
           call ncvpt_int2(ncid,vid,1,1,npsi0,istatus)

c            write(*,*)'before vid=ncvid(sigma_E_adj'

           call ncvid2(vid,ncid,'sigma_E_adj',istatus)
           call ncvpt_doubl2(ncid,vid,start,npsi0,sigma_E_adj,istatus)

c           write(*,*)'before vid=ncvid(rho_adj'

           call ncvid2(vid,ncid,'rho_adj',istatus)
           call ncvpt_doubl2(ncid,vid,start,npsi0,rho_adj,istatus)
c           write(*,*)'after vid=ncvpt(rho_adj'

       endif
c--------------------------
c     write OX conversion data for i_ox=2 case
c--------------------------
      goto 20
      if (i_ox.eq.2) then

         call ncvid2(vid,ncid,'i_ox_conversion',istatus)
         call ncvpt_int2(ncid,vid,1,nrayl,i_ox_conversion_nc,istatus)

         call ncvid2(vid,ncid,'transm_ox',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,transm_ox_nc,istatus)

         call ncvid2(vid,ncid,'cn_par_optimal',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_par_optimal_nc,istatus)

c         call ncvid2(vid,ncid,'cn_par_optimal',istatus)
c         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_par_optimal_nc,istatus)

         call ncvid2(vid,ncid,'cnpar_ox',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cnpar_ox_nc,istatus)

         call ncvid2(vid,ncid,'cn_b_gradpsi',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_b_gradpsi_nc,istatus)

      endif !i_ox=2 
 20   continue
c------------------------------------------------------------------
c     for total power absorbed at all reflections at all rays
c------------------------------------------------------------------
      call ncvid2(vid,ncid,'w_tot_pow_absorb_at_refl_nc',istatus)
      call ncvpt_doubl2(ncid,vid,1,1, w_tot_pow_absorb_at_refl_nc,
     +                  istatus)
      
C-----------------------------------------------------------------------
c
cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos2(ncid,istatus)
      call check_err(istatus)

      endif                    ! End write netcdf data

      deallocate (tem1,STAT=istat)
      return
      end
c
c


      subroutine writencdf_iray_status(netcdfnml)
      implicit none
c-----writes iray_status(nray,nfreq) into existing netcdf file: netcdfnml

      include 'param.i'  
      include 'one.i'
      include 'writencdf.i'

 
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc' 
c-----input      
      character(*) netcdfnml ! input filename

c---- some stuff for netCDF file --- 
      integer ncid,vid,istatus,error_code,ncvdef2,ncdid2,ncddef2
      integer nraysid,nfreqdim
      integer ray_dims(2),start(2),ray_count(2)
c-----local
      integer, pointer :: item2(:) !(nraya_nfreqa)
c      real*8, pointer :: tem2(:) !(nraya_nfreqa)
c      real*8,  pointer :: r_iray_status_nc(:,:)

      integer nraya_nfreqa,iray,ifreq

      data start/1,1/

      write(*,*)'in writencdf_iray_status netcdfnml=',netcdfnml

c------------------------------------------------------
c     allocate item2
c----------------------------------------------------------
      write(*,*)'nrayl,nfreq',nrayl,nfreq

      nraya_nfreqa=nrayl*nfreq

      allocate( item2(1:nraya_nfreqa),STAT=istatus)
      call ibcast(item2,0,SIZE(item2))

c      allocate( tem2(1:nraya_nfreqa),STAT=istatus)
c      call bcast(tem2,0.d0,SIZE(tem2))

c      allocate(r_iray_status_nc(1:nrayl,1:nfreq),STAT=istatus)
c      call bcast(r_iray_status_nc,0.d0,SIZE(r_iray_status_nc))
      

      write(*,*)'after allocate item2 istatus=',istatus
c.......................................................................
cl    1. Initialize part
c

C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes    
      write(*,*)'before ncopn2(netcdfnml'

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)
      write(*,*)'after ncopn(netcdfnml istatus=',istatus
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf2(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf2(integer ncid,integer error_code)
      call ncredf2(ncid,error_code)
      call check_err(istatus)
c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c     dimension named 'lat'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c     returns dimension id.
c     p. 54 of netcdf manual
c      write(*,*)'before ncdid(ncid nrays'
      nraysid=ncdid2(ncid,'nrays',istatus)
c      write(*,*)'after ncdid(ncid nrays'
c      write(*,*)'nraysid=',nraysid,'istaus=',istatus

c      write(*,*)'i_emission',i_emission
      if (i_emission.eq.1) then
c         write(*,*)'before nfreqdim=ncdid(ncid,nfreq'
         nfreqdim=ncdid2(ncid,'nfreq',istatus)
         call check_err(istatus)
c         write(*,*)'after nfreqdim=ncdid(ncid,nfreq istatus',istatus
c         write(*,*)'nfreqdim',nfreqdim
      else
c         write(*,*)'before nfreqdim=ncddef(ncid,nfreq ncid=',ncid 
         nfreqdim=ncddef2(ncid,'nfreq',nfreq,istatus)
         call check_err(istatus)
c         write(*,*)'after nfreqdim=ncddef(ncid,nfreq istatus',istatus
c         write(*,*)'nfreq',nfreq
c         write(*,*)'nfreqdim',nfreqdim
c-------------------------------------------------------------------------
c        Run Descriptive Parameters
c------------------------------------------------------------------------- 
         vid=ncvdef2(ncid,'nfreq',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     &   'The number of the emission frequencies',istatus)
         call check_err(istatus)
      endif

      ray_count(1)=nrayl
      ray_count(2)=nfreq
      
      write(*,*)'ray_count',ray_count

      ray_dims(1)=nraysid
      ray_dims(2)=nfreqdim

      write(*,*)'ray_dims',ray_dims
c-------------------------------------------------------------------------
c     Run Descriptive Parameters
c------------------------------------------------------------------------- 
c      write(*,*)'before vid=ncvdef(ncid,iray_status_nc'
      vid=ncvdef2(ncid,'iray_status_nc',NCLONG,2,ray_dims,istatus)
      call check_err(istatus)
c      write(*,*)'after vid=ncvdef(ncid,iray_status_n istatus',istatus
c      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
c    &        'type of stopping of each ray:
      call ncaptc2(ncid,vid,'long_name',NCCHAR,54,
     &'Status code giving the reason each ray stops          ', istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,54,
     &'=-1 bad initial condition;                            ', istatus)
      call ncaptc2(ncid,vid,'long_name2',NCCHAR,54,
     &' =1  t.gt.poldist_mx;                                 ', istatus)
      call ncaptc2(ncid,vid,'long_name3',NCCHAR,54,
     &' =2 in prep3d delpwr(is).lt.delpwrmn*delpwr(1)        ', istatus)
      call ncaptc2(ncid,vid,'long_name4',NCCHAR,54,
     &' =3  irefl.ge.ireflm;                                 ', istatus)
      call ncaptc2(ncid,vid,'long_name5',NCCHAR,54,
     &' =4  ray is at tor limiter bndry,  max_limiters.ge.1; ', istatus)
      call ncaptc2(ncid,vid,'long_name6',NCCHAR,54,
     &' =5  nstep_rk.gt.maxsteps_rk;                         ', istatus)
      call ncaptc2(ncid,vid,'long_name7',NCCHAR,54,
     &' =6  nmode.gt.cnmax, cnmax set in output.f cnmax=10000', istatus)
      call ncaptc2(ncid,vid,'long_name8',NCCHAR,54,
     &' =7  for ((id.eq.1,2).and.(uh.gt.1.d0).and.           ', istatus)
      call ncaptc2(ncid,vid,'long_name9',NCCHAR,54,
     &' ((uh-1.d0).lt.del_uh)).and.(cnper.gt.cnper_max_ebw)).', istatus)
      call ncaptc2(ncid,vid,'long_name10',NCCHAR,54,
     &'  The ray is close to upper hybrid resonance.         ', istatus)
      call ncaptc2(ncid,vid,'long_name11',NCCHAR,54,
     &'  It can be Xmode to close to the UH resonance.       ', istatus)
      call ncaptc2(ncid,vid,'long_name12',NCCHAR,54,
     &' =8  if ((vgrmods.gt.1.1).and. ((id.ne.10).and.       ', istatus)
      call ncaptc2(ncid,vid,'long_name13',NCCHAR,54,
     &  '(id.ne.12).and.(id.ne.13))) and if (id.eq.14);      ', istatus)
      call ncaptc2(ncid,vid,'long_name14',NCCHAR,54,
     &' =9  D/(N|gradD|) > toll_hamilt;                      ', istatus)
      call ncaptc2(ncid,vid,'long_name15',NCCHAR,54,
     &' =10 in prep3d.f argexp>0;                            ', istatus)
      call ncaptc2(ncid,vid,'long_name16',NCCHAR,54,
     &' =11 in prep3d fluxn(is).lt.0.0;                      ', istatus)
      call ncaptc2(ncid,vid,'long_name17',NCCHAR,54,
     &' =12 nrayelt=nrelt;                                   ', istatus)
      call ncaptc2(ncid,vid,'long_name18',NCCHAR,54,
     &' =13 in RK drkgs2 small time step h.lt.1.d-11         ', istatus)
      call ncaptc2(ncid,vid,'long_name19',NCCHAR,54,
     &'  can not get the given accuracy prmt4                ', istatus)
      call ncaptc2(ncid,vid,'long_name20',NCCHAR,54,
     &'  =14, ray started with zero power (lt.1e-100)        ', istatus)

      call check_err(istatus)
c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual 
      call ncendf2(ncid,istatus)
      call check_err(istatus)
c     End initialize
     
c.......................................................................
cl    1. Writing data

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual

c------------------------- 
      write(*,*)'writing data i_emission',i_emission
      if (i_emission.eq.0) then
         write(*,*)'before ncvid(ncid,'
         call ncvid2(vid,ncid,'nfreq',istatus)
         write(*,*)'after ncvid(ncid istatus',istatus
         call ncvpt_int2(ncid,vid,1,1,nfreq,istatus)
         write(*,*)'after  ncvpt(ncid,vid,1,1,nfreq istatus',istatus
         do iray=1,nrayl
            item2(iray)=iray_status_nc(iray,1)
         enddo
      else
c         call pack21(iray_status_nc,1,nrayl,1,nfreq,tem2,nrayl,nfreq)
          call pack21_integer(iray_status_nc,1,nrayl,1,nfreq,item2,
     &    nrayl,nfreq)      
      endif

cPrint out to screen the defn of iray_status_nc, for convenience:

      write(*,*) \\
     &  "iray_status_nc: Status code giving the reason each ray stops"\
     &	"iray_status_nc:=-1 bad initial condition"\
     &  "iray_status_nc:= 1 t.gt.poldist_mx"\
     &  "iray_status_nc:= 2 in prep3d delpwr(is).lt.delpwrmn*delpwr(1)"\
     &  "iray_status_nc: =3 irefl.ge.ireflm"\
     &  "iray_status_nc: =4 ray at tor limiter bndry, max_limtrs.ge.1"\
     &  "iray_status_nc: =5 nstep_rk.gt.maxsteps_rk"\
     &  "iray_status_nc: =6 nmode.gt.cnmax, in output.f cnmax=10000"\
     &  "iray_status_nc: =7 for ((id.eq.1,2).and.(uh.gt.1.d0).and."\
     &  "      ((uh-1.d0).lt.del_uh)).and.(cnper.gt.cnper_max_ebw))."\
     &  "      The ray is close to upper hybrid resonance.         "\
     &  "      It can be Xmode to close to the UH resonance.       "\
     &  "iray_status_nc: =8 if ((vgrmods.gt.1.1).and. ((id.ne.10)"\
     &  "      .and.(id.ne.12).and.(id.ne.13))) and if (id.eq.14)  lo"\
     &  "iray_status_nc: =9 D/(N|gradD|) > toll_hamilt"\
     &  "iray_status_nc: =10 in prep3d.f argexp>0"\
     &  "iray_status_nc: =11 in prep3d fluxn(is).lt.0.0 "\ 
     &  "iray_status_nc: =12 nrayelt=nrelt "\
     &  "iray_status_nc: =13 in RK drkgs2 small time step h.lt.1.d-11"\
     &  "      can not get the given accuracy prmt4"\
     &  "iray_status_nc:  =14, ray started with zero power (lt.1e-100) "

      do iray=1,nrayl
        do ifreq=1,nfreq
           write(*,*)'iray,ifreq,iray_status_nc(iray,ifreq)',
     &               iray,ifreq,iray_status_nc(iray,ifreq)
        enddo
      enddo

c      do iray=1,nrayl
c        do ifreq=1,nfreq
c           r_iray_status_nc(iray,ifreq)=iray_status_nc(iray,ifreq)             
c        enddo
c      enddo
      

c      write(*,*)'before pack21(iray_status_nc nrayl,nfreq',nrayl,nfreq
c      call pack21_integer(r_iray_status_nc,1,nrayl,1,nfreq,tem2,
c     &nrayl,nfreq)       

c      call pack21_integer(iray_status_nc,1,nrayl,1,nfreq,item2,
c     &nrayl,nfreq)      

c      write(*,*)'after pack21(iray_status_nc tem2',tem2
c      do iray=1,nraya_nfreqa
c         item2(iray)=tem2(iray)
c      enddo
c      goto 10

      call ncvid2(vid,ncid,'iray_status_nc',istatus)

      write(*,*)'after ncvid(ncid,iray_status_nc istatus',istatus
      write(*,*)'ray_count',ray_count
      call ncvpt_int2(ncid,vid,start,ray_count,item2,istatus)
      write(*,*)'after ncvpt(ncid,vid,start, istatus',istatus

 10   continue
      deallocate (item2,STAT=istatus)
c      deallocate (tem2,STAT=istatus)
c      deallocate (r_iray_status_nc,STAT=istatus)

c    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
 

      call ncclos2(ncid,istatus)
      write(*,*)'in writencdf_iray_status after ncclos2'
      call check_err(istatus)

      return
      end





C=YuP=> ADDED: conversion from Netcdf-2 to Netcdf-3 or higher ==========
C These routines/function convert old routines:

      integer function ncvdef2(NCID,VARNAM,VARTYP,NDIMS,VDIMS,istatus)
      ! vid=ncvdef2() is renamed to vid=ncvdef2() in *.f files
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) VARNAM
      INTEGER istatus,vid,NCID, VARTYP,XTYPE, NDIMS,VDIMS(*)
      if (VARTYP.eq.NCFLOAT)  XTYPE=NF_FLOAT    ! 32 BITS
      if (VARTYP.eq.NCDOUBLE) XTYPE=NF_DOUBLE   ! 64 BITS
      if (VARTYP.eq.NCCHAR)   XTYPE=NF_CHAR
      if (VARTYP.eq.NCBYTE)   XTYPE=NF_BYTE
      if (VARTYP.eq.NCSHORT)  XTYPE=NF_SHORT    ! 16 BITS
      if (VARTYP.eq.NCLONG)   XTYPE=NF_INT      ! 32 BITS
      istatus = NF_DEF_VAR(NCID,VARNAM,XTYPE,NDIMS,VDIMS,vid)
      ncvdef2 = vid
      end

      integer function ncdid2(ncid,VARNAM,istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) VARNAM
      INTEGER istatus,ncid,ndim
c     Get dimension ID from dimension name
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c-YuP:      neltdim=ncdid(ncid,'neltmax',istatus)
      istatus= NF_INQ_DIMID(ncid,VARNAM,ndim) 
      ncdid2 = ndim
      end

      integer function ncddef2(ncid,VARNAM,LEN,istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) VARNAM
      INTEGER istatus,ncid,LEN,RECID
c     returns dimension id.
      if(LEN.eq.NCUNLIM) then
        ! unlimited dimension for time, dimension name= 'time'
        istatus= NF_DEF_DIM(ncid, VARNAM, NF_UNLIMITED, RECID)
      else
        istatus= NF_DEF_DIM(ncid, VARNAM, LEN, RECID)
      endif
      ncddef2 = RECID
      end

      integer function nccre2(filename,MODE,istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) filename
      INTEGER istatus,ncid,MODE
cl    create netCDF filename (Entering define mode.)
      if(MODE.eq.NCCLOB)then
        ! over-write existing file
        istatus= NF_CREATE(filename, NF_CLOBBER, ncid) 
      else
        ! Do not over-write existing file
        istatus= NF_CREATE(filename, NF_NOCLOBBER, ncid) 
      endif
      nccre2 = ncid
      end

      subroutine ncvgtc2(NCID, vid, START, COUNTS, TEXT, LEN, istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) TEXT
      INTEGER istatus,vid,NCID, LEN, START(*), COUNTS(*)
      istatus= NF_GET_VARA_TEXT (NCID, vid, START, COUNTS, TEXT)
      return 
      end

      subroutine ncvgt_int2(NCID, vid, START, COUNTS,  IVALS, istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      INTEGER IVALS(*)
      istatus= NF_GET_VARA_INT(NCID, vid, START, COUNTS, IVALS)
      return 
      end

      subroutine ncvgt_doubl2(NCID, vid, START, COUNTS,  DVALS, istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      REAL*8 DVALS(*)
      istatus= NF_GET_VARA_DOUBLE(NCID, vid, START, COUNTS, DVALS)
      return 
      end
      
      
      subroutine ncaptc2(NCID, vid, NAME, ATTYPE, LEN, TEXT, istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) NAME
      CHARACTER*(*) TEXT
      INTEGER istatus,vid,NCID, LEN, ATTYPE
      istatus= NF_PUT_ATT_TEXT(NCID, vid, NAME, LEN, TEXT)
      return 
      end
        
      subroutine ncvptc2(NCID, vid, START, COUNTS, TEXT, LEN, istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) TEXT
      INTEGER istatus,vid,NCID, LEN, START(*), COUNTS(*)
      istatus= NF_PUT_VARA_TEXT(NCID, vid, START, COUNTS, TEXT)
      return 
      end

      subroutine ncvpt_doubl2(NCID, vid, START, COUNTS,  DVALS, istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      REAL*8 DVALS(*)
      istatus= NF_PUT_VARA_DOUBLE(NCID, vid, START, COUNTS, DVALS)
      return 
      end

      subroutine ncvpt_int2(NCID, vid, START, COUNTS,  IVALS, istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      INTEGER IVALS(*)
      istatus= NF_PUT_VARA_INT(NCID, vid, START, COUNTS, IVALS)
      return 
      end
      
      subroutine ncopn2(filename,MODE,ncid,istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,MODE,ncid
      CHARACTER*(*) filename
      !integer function ncopn(filename,write,error_code) ! NetCDF-2
      !ncid=ncopn(netcdfnml,NCWRITE,istatus) !Open existing netCDF file
      if(MODE.eq.0) then
         istatus= NF_OPEN(filename, NF_NOWRITE, ncid) 
         else
         istatus= NF_OPEN(filename, NF_WRITE, ncid) 
      endif
      if (istatus .NE. NF_NOERR) then         
         write(*,*)'   ***   Problem opening .nc data file   ***'
         Stop
      endif  
      return                                
      end

      subroutine ncclos2(ncid,istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,ncid
      istatus= NF_CLOSE(ncid)
      return 
      end

      subroutine ncvid2(vid,ncid,NAME,istatus)
      ! vid=ncvid(ncid,NAME,istatus) ! NetCDF-2
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) NAME
      INTEGER istatus,ncid,vid
      istatus= NF_INQ_VARID(ncid,NAME,vid)
      return
      end

      subroutine ncdinq2(ncid,DIMID,NAME,LEN,istatus) 
c     Query netcdf file for dimensions:
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) NAME
      INTEGER istatus,ncid,DIMID,LEN
      istatus= NF_INQ_DIM(ncid,DIMID,NAME,LEN)
      return
      end

      subroutine ncredf2(ncid,istatus)
c     Put Open NetCDF file into Define Mode
      INCLUDE 'netcdf.inc'
      INTEGER istatus,ncid
      istatus= NF_REDEF(ncid) ! put into define mode
      return
      end

      subroutine ncendf2(ncid,istatus)
c     end the define-mode and start the data-mode
      INCLUDE 'netcdf.inc'
      INTEGER istatus,ncid
      istatus= NF_ENDDEF(ncid) 
      return
      end

         
