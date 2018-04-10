c        ********************** write3d1*********************
c        *                      -----                       *
c        *  write3d1 -subroutine writes ray data:           *
c        *  nray,nharm,freqcy    for	 FP code            *
c        ****************************************************
c        *  input parameters:nray -number of the rays at antenna  *
c        *  jwave -number of harmonic from common block one	  *
c	 *  frqncy-wave friquency from common block one		  *
c        *  i_ -number of output file(for 3D FP) from common one  *
c-----------------------------------------------------------------*
      subroutine write3d1(nray)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'writencdf.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c     open(io9,file=filetxt) in genray.f
c---------------------------------------------------------------------
c     nray- total number of rays
c     nharm-harmonic for ql diffusion (0,for Landau/TTMP;
c     ioxm=1 OM,ioxm=-1 XM for ECR)
c     freqcy-wave antenna frequency(Hz)
c---------------------------------------------------------------------
c     for ECR wave
cSmirnov961206 begin
c      nharm=ioxm
      nharm=jwave
cSmirnov961206 end
      freqcy=frqncy*1.d+9
c      write(*,*) 'in write3d1 i_',i_

      if(myrank.eq.0) then
      write(*,*) 'in write3d1 nray,nharm,freqcy',nray,nharm,freqcy
      if(rayop.eq."text" .or. rayop.eq."both") then
         write(i_,1) nray,nharm,freqcy
 1       format (2i5,1pe16.9)
      endif
      endif ! myrank=0
     
c-----initialization of the number of the ray 
      irayl=0
      nrayl=nray
      if(myrank.eq.0) write(*,*)'in sub write3d1 nrayl=',nrayl
      return
      end


c     ********************** write3d**********************
c     *                      -----                       *
c     *  write3d -subroutine writes ray data             *
c     *  for  3D code                         	         *
c     ****************************************************
c        ****************************************************
c        *  input parameters:                                     *
c        *  i_ -number of output file(for 3D FP) from common one  *
c        *  nrayelt -number of ray elements along ray iray from   *
c        *           common block write(was calculated in prep3d) *				  *
c-----------------------------------------------------------------*
      subroutine write3d
      implicit double precision (a-h,o-z)
      dimension u(6),deru(6)
      include 'param.i'
      include 'write.i'
      include 'one.i'
      include 'writencdf.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      if(myrank.ne.0) return

c	nrayelt- number of ray elements for ray iray
c       The next nine integers in the write are the final states of the mode,
c       number of refl., integ.steps, number of electron and ion power calcs
c       and the restart status of the ray. (For EC, can be set equel to 0)
c       (Also can be 0, if no restart needed).
c       For EC, next 4 real numbers may be zero).
c       sxxrt(iray,k) - final value of integ.variable
c       skpsi(iray,k) - cannonical n_psi
c       skth(iray,k)  -canonical n_th
c       skphi - canonical n_phi
c-----------------------------------------------------------
c       The next 8 integers =0 for test
      jslofas=0
      nurefls=0
      keiks=0
      jpes=0
      jpis=0
      istarts=0
      jprmt5=0
      jhlfs=0
c-------------------------------------------------------------
c     For EC
      sxxrt=0.d0
      skpsi=0.d0
      skth=0.d0
      skphi=0.d0
            
      write(*,*)'in write3d i_',i_
      write(*,*)'in write3d rayop=',rayop

      if(rayop.eq."text" .or. rayop.eq."both") then
      
c      write(*,*)'nrayelt,jslofas,nurefls,keiks,jpes,jpis,istarts,
c     +           jprmt5,jhlfs',
c     +           nrayelt,jslofas,nurefls,keiks,jpes,jpis,istarts,
c     +           jprmt5,jhlfs
      write(i_,4) nrayelt,jslofas,
     +           nurefls,keiks,jpes,
     +           jpis,istarts,
     +           jprmt5,jhlfs
c--------------------------------------------------------------
      write(i_,2)sxxrt,skpsi,skth,
     1           skphi
c--------------------------------------------------------------
      io9=i_
c     open(io9,file=filetxt) in genray.f
c     write(io9,4) nrayelt
      if (nrayelt.eq.0) then
c        ray iray can not go to plasma
         goto 9
      endif
      write(io9,2) (ws(is),is=1,nrayelt)
      write(io9,2) (seikon(is),is=1,nrayelt)
      write(io9,2) (spsi(is),is=1,nrayelt)
      write(io9,2) (wr(is),is=1,nrayelt)
      write(io9,2) (wphi(is),is=1,nrayelt)
      write(io9,2) (wz(is),is=1,nrayelt)
      write(io9,2) (wnpar(is),is=1,nrayelt)
      write(io9,2) (wnper(is),is=1,nrayelt)
      write(io9,2) (delpwr(is),is=1,nrayelt)
      write(io9,2) (sdpwr(is),is=1,nrayelt)
      write(io9,2) (wdnpar(is),is=1,nrayelt)
      write(io9,2) (cwexde(is),is=1,nrayelt)
      write(io9,2) (cweyde(is),is=1,nrayelt)
      write(io9,2) (cwezde(is),is=1,nrayelt)
      write(io9,2) (fluxn(is),is=1,nrayelt)
      write(io9,2) (sbtot(is),is=1,nrayelt)
      write(io9,2) (sene(is),is=1,nrayelt)
      write(io9,2) (salphac(is),is=1,nrayelt)
      write(io9,2) (salphal(is),is=1,nrayelt)
c     write(*,*)'cwezde',(cwezde(is),is=1,nrayelt)
c     write(*,*)'rez',(rez(is),is=1,nrayelt)
 1    format (2i5,1pe16.9)
 2    format (5(1pe16.9))
 3    format (9i5)
 4    format (9(' ',i4))

 9    continue

      endif  !On rayop

c-----the output ray data for mnemonic.nc netCDF file
      irayl=irayl+1 

      write(*,*)'write3d irayl,nrayelt',irayl,nrayelt

      nrayelt_nc(irayl)=nrayelt

      write(*,*)'write3d irayl,nrayelt',irayl,nrayelt

      if (nrayelt.eq.0) goto 10

      do is=1,nrayelt
         ws_nc(is,irayl)=ws(is)
         seikon_nc(is,irayl)=seikon(is)
         spsi_nc(is,irayl)=spsi(is)
         wr_nc(is,irayl)=wr(is)
         wphi_nc(is,irayl)=wphi(is)
         wz_nc(is,irayl)=wz(is)
         wnpar_nc(is,irayl)=wnpar(is)
         wnper_nc(is,irayl)=wnper(is)
         delpwr_nc(is,irayl)=delpwr(is)
         sdpwr_nc(is,irayl)=sdpwr(is)
         do i=2,nbulk
            salphas_nc(is,irayl,i)=salphas(is,i)
         enddo
         wdnpar_nc(is,irayl)=wdnpar(is)
         cwexde_nc(is,irayl)=cwexde(is)
         cweyde_nc(is,irayl)=cweyde(is)
         cwezde_nc(is,irayl)=cwezde(is)
         fluxn_nc(is,irayl)=fluxn(is)
         sbtot_nc(is,irayl)=sbtot(is)
         sene_nc(is,irayl)=sene(is)
         ste_nc(is,irayl)=ste(is)
         szeff_nc(is,irayl)=szeff(is)
         salphac_nc(is,irayl)=salphac(is)

c         write(*,*)'write3d.f i_salphal=',i_salphal
cSm060913
cSm061127         if((iabsorp.eq.3).or.(iabsorp.eq.9)) then
         if(((iabsorp.eq.3).or.(iabsorp.eq.9)).or.
     &      ((iabsorp.eq.91).or.(iabsorp.eq.92))) then
           salphal_nc(is,irayl)=0.d0
           do i=1,nbulk
c             write(*,*)'in write3d i,i_salphal(i)',i,i_salphal(i)
             if (i.eq.1)then
               if(i_salphal(1).eq.1)then 
                salphal_nc(is,irayl)=salphal_nc(is,irayl)+salphal(is)
               endif
             else  !i>1
               if(i_salphal(i).eq.1)then 
                salphal_nc(is,irayl)=salphal_nc(is,irayl)+salphas(is,i)
c          write(*,*)'in write3d is,salphas(is,i),salphal_nc(is,irayl)',
c     &                          is,salphas(is,i),salphal_nc(is,irayl)
               endif
             endif
           enddo
         else
           salphal_nc(is,irayl)=salphal(is)
         endif!(iabsorp.eq.3.or.iabsorp.eq.9)=91, =92 

c         salphal_nc(is,irayl)=salphal(is)

         sb_z_nc(is,irayl)=sb_z(is)
         sb_r_nc(is,irayl)=sb_r(is)
         sb_phi_nc(is,irayl)=sb_phi(is)
         wn_z_nc(is,irayl)=wn_z(is)
         wn_r_nc(is,irayl)=wn_r(is)
         wn_phi_nc(is,irayl)=wn_phi(is)
         vgr_z_nc(is,irayl)=wvgr_z(is)
         vgr_r_nc(is,irayl)=wvgr_r(is)
         vgr_phi_nc(is,irayl)=wvgr_phi(is)
c--------tokman flux
         flux_z_nc(is,irayl)=wflux_tokman_z(is)
         flux_r_nc(is,irayl)=wflux_tokman_r(is)
         flux_phi_nc(is,irayl)=wflux_tokman_phi(is)
c--------dielectric tensor
cSAP090710
        if (dielectric_op.eq.'enabled') then
           cweps11_nc(is,irayl)= w_ceps(1,1,is)
           cweps12_nc(is,irayl)= w_ceps(1,2,is)
           cweps13_nc(is,irayl)= w_ceps(1,3,is)
           cweps21_nc(is,irayl)= w_ceps(2,1,is)
           cweps22_nc(is,irayl)= w_ceps(2,2,is)
           cweps23_nc(is,irayl)= w_ceps(2,3,is)
           cweps31_nc(is,irayl)= w_ceps(3,1,is)
           cweps32_nc(is,irayl)= w_ceps(3,2,is)
           cweps33_nc(is,irayl)= w_ceps(3,3,is)
         endif
c--------CD efficiency
         if(ionetwo.eq.1) w_eff_nc(is,irayl)=eff(is)
         w_theta_pol_nc(is,irayl)=wtheta_pol(is)
      enddo

      if (i_ox.eq.2) then
c--------the data for OX conversion 
         i_ox_conversion_nc(irayl) = i_ox_conversion
         if (i_ox_conversion.eq.0) then
           transm_ox_nc(irayl) = 0.d0
           cn_par_optimal_nc(irayl) =0.d0 
           cnpar_ox_nc(irayl)= 0.d0
           cn_b_gradpsi_nc(irayl)= 0.d0
         else
           transm_ox_nc(irayl) = transm_ox
           cn_par_optimal_nc(irayl) = cn_par_optimal
           cnpar_ox_nc(irayl)= cnpar_ox
           cn_b_gradpsi_nc(irayl)=cn_b_gradpsi
         endif

         write(*,*)'write3d i_ox.eq.2'
         write(*,*)'irayl,i_ox_conversion_nc(irayl)',
     &              irayl,i_ox_conversion_nc(irayl)
         write(*,*)'transm_ox_nc(irayl)',transm_ox_nc(irayl)
         write(*,*)'cn_par_optimal_nc(irayl)',cn_par_optimal_nc(irayl)
         write(*,*)'cnpar_ox_nc(irayl)',cnpar_ox_nc(irayl)
         write(*,*)'cn_b_gradpsi_nc(irayl)',cn_b_gradpsi_nc(irayl)
      endif
      w_tot_pow_absorb_at_refl_nc=w_tot_pow_absorb_at_refl

 10   continue
      write(*,*)'write3d after 10 irayl,nrayelt',irayl,nrayelt
       
      return
      end

c     ********************** read3d**********************
c     *                      -----                      *
c     *  read3d -subroutine reads  ray data             *
c     *  from old_3d.dat                    	        *
c     ***************************************************
c        ****************************************************
c        *  input parameters:                                     *
c        *  i_ -number of output file(for 3D FP) from common one  *
c        *  nrayelt -number of ray elements along ray iray from   *
c        *           common block write(was calculated in prep3d)
c-----------------------------------------------------------------*
      subroutine read3d(io)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'write.i'
      include 'one.i'
c      io is the number of the open file.
c	nrayelt- number of ray elements for ray iray
c       The next nine integers in the write are the final states of the mode,
c       number of refl., integ.steps, number of electron and ion power calcs
c       and the restart status of the ray. (For EC, can be set equel to 0)
c       (Also can be 0, if no restart needed).
c       For EC, next 4 real numbers may be zero).
c       sxxrt(iray,k) - final value of integ.variable
c       skpsi(iray,k) - cannonical n_psi
c       skth(iray,k)  -canonical n_th
c       skphi - canonical n_phi
c-----------------------------------------------------------
c        read ray data:       nray,nharm,freqcy                
c        *  input parameters:nray -number of the rays at antenna  *
c        *  jwave -number of harmonic from common block one	  *
c	 *  frqncy-wave friquency from common block one		  *
c        *  i_ -number of output file(for 3D FP) from common one  *
c-----------------------------------------------------------------*
      open(io,file='old_3d.dat',status='old') 
c---------------------------------------------------------------------
c     nray- total number of rays
c     nharm-harmonic for ql diffusion (0,for Landau/TTMP;
c     ioxm=1 OM,ioxm=-1 XM for ECR)
c     freqcy-wave antenna frequency(Hz)
c--------------------------------------------------------------

      write(*,*) 'in read3d'
      read(io,1)nray,nharm,freqcy
 11    format (2i5,1pe16.9)
c      write(*,*)'nray,nharm,freqcy'
c      write(*,11)nray,nharm,freqcy

c--------------------------------------------------------------
      read(io,4) nrayelt,jslofas,
     +           nurefls,keiks,jpes,
     +           jpis,istarts,
     +           jprmt5,jhlfs
c      write(*,*)'nrayelt,jslofas,nurefls,keiks,jpes,jpis,istarts,
c     +           jprmt5,jhlfs'
c      write(*,4) nrayelt,jslofas,
c     +           nurefls,keiks,jpes,
c     +           jpis,istarts,
c     +           jprmt5,jhlfs

c--------------------------------------------------------------
      read(io,2)sxxrt,skpsi,skth,skphi
c      write(*,*)'sxxrt,skpsi,skth,skphi'
c      write(*,2)sxxrt,skpsi,skth,skphi

c-----------------------------------------------------------

      if (nrayelt.eq.0) then
c        ray iray can not go to plasma
         goto 10
      endif

      read(io,2) (ws(is),is=1,nrayelt)
      read(io,2) (seikon(is),is=1,nrayelt)
      read(io,2) (spsi(is),is=1,nrayelt)
      read(io,2) (wr(is),is=1,nrayelt)
      read(io,2) (wphi(is),is=1,nrayelt)
      read(io,2) (wz(is),is=1,nrayelt)
      read(io,2) (wnpar(is),is=1,nrayelt)
      read(io,2) (wnper(is),is=1,nrayelt)
      read(io,2) (delpwr(is),is=1,nrayelt)
      read(io,2) (sdpwr(is),is=1,nrayelt)
      read(io,2) (wdnpar(is),is=1,nrayelt)
      read(io,2) (cwexde(is),is=1,nrayelt)
      read(io,2) (cweyde(is),is=1,nrayelt)
      read(io,2) (cwezde(is),is=1,nrayelt)
      read(io,2) (fluxn(is),is=1,nrayelt)
      read(io,2) (sbtot(is),is=1,nrayelt)
      read(io,2) (sene(is),is=1,nrayelt)
      read(io,2) (salphac(is),is=1,nrayelt)
      read(io,2) (salphal(is),is=1,nrayelt)

 1    format (2i5,1pe16.9)
 2    format (5(1pe16.9))
 3    format (9i5)
 4    format (9(' ',i4))
 10   continue
   
      return
      end



      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
c      write(*,*) nf_strerror(iret)
         write(*,*) 'netCDF error'
      stop 'check_err:'
      endif
      end
c
c



cfor test
      subroutine pack21_test(a,ibot,itop,jbot,jtop,b,iy,jx)
      implicit integer (i-n), real*8 (a-h,o-z)

c.......................................................................
c     It sometimes becomes necessary to take a
c     2-D array dimensioned ibot:itop by jbot:jtop
c     and repack it as though it were
c     dimensioned 1:iy by 1:jx, starting at a(1,1).
c     This routine does this, transfering relevant data
c     from array a to b.
c.......................................................................

      save
      dimension a(ibot:itop,jbot:jtop)
      dimension b(iy*jx)
      write(*,*)'pack21 ibot,itop,jbot,jtop',ibot,itop,jbot,jtop
      write(*,*)'pack21 iy,jx',iy,jx

      do 1 j=1,jx
        i1=(j-1)*iy+1
c        call scopy(iy,a(1,j),1,b(i1),1)

c         write(*,*)'pack21 a(1,j)'
c         do i=1,iy
c           write(*,*)'i,a(i,j)',i,a(i,j)
c         enddo
        write(*,*)'j,i1',j,i1
        call dcopy_test(iy,a(1,j),1,b(i1) ,1)

       write(*,*)'pack21_test b(i1,j)'
c        do k=i1,i1+jx-1
        write(*,*)'i1,i1+jx-2',i1,i1+jx-2
        do k=i1,i1+jx-2
           write(*,*)'k,b(k)',k,b(k)
       enddo

 1    continue
      return
      end
c
      subroutine dcopy_test(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real*8 dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      write(*,*)'dcopy incx,incy',incx,incy

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      write(*,*)'dcopy_test n,m',n,m
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
        write(*,*)'i,dy(i),dx(i)',i,dy(i),dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      write(*,*)'dcopy_test mp1',mp1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
        write(*,*)'dx(i),dx(i+1),dx(i+2),dx(i+3)',
     &             dx(i),dx(i+1),dx(i+2),dx(i+3)
        write(*,*)'dx(i+4),dx(i+5),dx(i+6)',
     &             dx(i+4),dx(i+5),dx(i+6)
   50 continue
      return
      end

