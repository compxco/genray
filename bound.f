c        *********************BOUND****************************
c        *                        -                           *
c        * This function gives the ray position in            *
c        * the plasma (If the point of the ray is inside the  *
c        * limitter the BOUND=-1 , else BOUND=1),	      *
c        * and it reflects the ray from plasma boundary.      *
c        * It calculates the number of the reflections-irefl  *
c        * It gives the command to stop ray calculation	      *
c        * if irefl.ge.ireflm				      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        Input parameters					   !
c      								   !
c      Z, R, PHI - point which the components of magnetic field	   !
c                  are calculated in.
c
c     the input parameter epsbnd(the small distance from the boundary)
c     is set inside this subroutine
c------------------------------------------------------------------!
c     output parameters:iflref=1 after reflection =-1 before refl.
c			cnzref,cnrref,cmref= N_z,N_r,M after reflectin point,or
c				     =N_z,N_r,M before reflection point
c                       ibound=1 point is outside the plasma =-1 in plasma
c----------------------------------------------------------------------
      subroutine bound(z,r,phi,cnz,cnr,cm,iflref,
     &z_ref,r_ref,phi_ref,cnzref,cnrref,cmref,
     &ibound,dzdt,drdt)
      !implicit double precision (a-h,o-z)
      implicit none

c-----input
      real*8 
     &z,r,phi,     !space coordinates
     &cnz,cnr,cm,  !refractive index coordinates
     &dzdt,drdt    !wave group velocity in RHS of the ray-tracing equations
c-----output
      integer
     &iflref,     !iflref=1 after reflection =-1 before refl.
     &ibound      !ibound=1 point is outside the plasma =-1 in plasma

      real*8 
     &z_ref,r_ref,phi_ref, ! reflection point coordinates
     &cnzref,cnrref,cmref  ! N_z,N_r,M after reflectin point,or
                           ! =N_z,N_r,M before reflection point
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
cSm041006
      include 'rkutta.i'
      include 'oxb.i' ! YuP[2020-07-29] needed now for i_ox=1 run.
      !oxb.i now stores the value of rhoconv; 
      ! rhoconv is rho at Xe=1.0 conversion layer,
      ! calculated in sub.ox_conversion_grill_in_poloidal_point.
      ! It also stores:
      ! rho_plasm_vac is rho at plasma-vacuum border;
      ! normally, it would be at rho=1.0, but if edge density is high,
                    ! Xe=1 layer could be just outside of rho=1.0
      ! Then, we define rho_plasm_vac as max(rhoconv,1.0)+0.1
      !  (or slightly different, see sub.bound)
      ! The value of rho_plasm_vac is used in sub.bound in i_ox=1 run.
      
c-----externals
      real*8
     &ias1r,ias2r,ias2r_Sm,b,rhopsi

c-----locals
      real*8
     &epsbnd,   !accuracy of point near the plasma edge LCFS
     &zp,zm,zzrm,zzrp,res,rrr,dmrdt,dpsidr,dpsidz,
     &edg,edag,grad

      integer nx4,ny4,ipx,ipx4,idx,
     &i_reflection  !=1 it was wall limiter reflection
                    !=0 it was not wall limiter reflection
                    
      nx4=nx+4
      ny4=ny+4

      cnzref=cnz
      cnrref=cnr
      cmref=cm
      z_ref=z
      r_ref=r
      phi_ref=phi

      iflref=-1 !initialize;   iflref-->1 after reflection; 
      !epsbnd=1.d-8
cSm070121
      !epsbnd=1.d-7
cSm070613
      !epsbnd=1.d-3
cSm070613
      epsbnd=1.d-5 ![m] ! should be same in reflect, bound, boundc

      ibound=-1 ! initialize: -1 means inside plasma
      !write(*,*)'bound: rhoconv=',rhoconv

c-----------------------  YuP[2020-07-23] added more checks
      ! Check: if the ray is out of (R,Z)-grid:
      if ( (r.lt.xeqmin+epsbnd).or.(r.ge.xeqmax-epsbnd) .or.
     +     (z.lt.zeqmin+epsbnd).or.(z.ge.zeqmax-epsbnd)  ) then
         ibound=1 !here: outside of B-grid
         irefl=ireflm !  to stop the ray
         write(*,*)'bound: out of grid. z,r=',z,r
         goto 10 !-> procedure to make a reflection/stop
      end if
c--------------------------------------------------------------YuP
      
      !YuP: no_reflection=1 means artificial reflection is OFF.
      !The next section will check - is ray crossing the chamber?
      !It seems that no_reflection setting has nothing to do with this.
      !Maybe drop if(no_reflection.eq.1) condition ?
      !YuP: For i_ox.eq.1, skip this part anyway. 
      !The M0 point (see oxb_launch.pdf) is supposed to be close to LCFS,
      !and will be detected by rho>1.0 below.
      if(i_ox.ne.1)then !YuP[2020-07-29] Added condition
      if(no_reflection.eq.1)then !YuP: =0 means artificial reflection is ON
!          write(*,*)'before wall_limiter_reflection_point'
!          write(*,*)'z,r,phi,cnz,cnr,cm',
!     &               z,r,phi,cnz,cnr,cm
          call  wall_limiter_reflection_point(z,r,phi,cnz,cnr,cm,
     &    i_reflection, cnzref,cnrref,cmref,z_ref,r_ref,phi_ref)
          ! i_reflection becomes 1 only when old and new points
          ! of the ray are on opposite sides of chamber wall;
          ! otherwise i_reflection stays equal to 0.

!          write(*,*)'after wall_limiter_reflection_point i_reflection=',
!     &    i_reflection

          if (i_reflection.eq.1) then  
!             write(*,*)'z,r,phi,cnz,cnr,cm',
!     &                  z,r,phi,cnz,cnr,cm
!             write(*,*)'z_ref,r_ref,phi_ref',z_ref,r_ref,phi_ref
!             write(*,*)'cnzref,cnrref,cmref',
!     &                  cnzref,cnrref,cmref
             irefl=irefl+1
             ibound=1 !point is outside of plasma (here: outside of wall)
             iflref=1 ! iflref=1 after reflection; =-1 before reflection
          endif
          goto 30 !-> return/end
      endif ! no_reflection.eq.1  means artificial reflection is OFF
      endif !(i_ox.ne.1) !YuP[2020-07-29] Added condition
      
c-----------------------
      if ((r.le.rmin+epsbnd).or.(r.ge.rmax-epsbnd)) then
         write(*,*)'in bound rmin,r,rmax,epsbnd',rmin,r,rmax,epsbnd
         ibound=1 !here R < min(R_LCFS)  or  R > max(R_LCFS)
         goto 10 !-> procedure to make a reflection/stop
      end if
      
      rrr=r !=R major radius
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      !Check where the ray element is with resp. to the Z(R) for the LCFS
      idx=0
      ipx=ip
      ipx4=ip+4
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr) !Z(R) for the upper half of LCFS
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr) !Z(R) for the lower half of LCFS
      zm=zzrm
      ! Check that the ray is inside flux surface
      if ((z.ge.zp-epsbnd).or.(z.le.zm+epsbnd)) then
cyup         write(*,*)'in bound r,zm,z,zp,epsbnd',r,zm,z,zp,epsbnd
         ibound=1 !here Z>Z_LCFS_upper(R) or  Z<Z_LCFS_lower(R)
         goto 10 !-> procedure to make a reflection/stop
      endif
      
c-----------------------
      !Check where the ray element is with resp. to the psilim
      nx4=nx+4
      ny4=ny+4
      ncx=nx4
      ncy=ny4
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nx4a)
      rho=rhopsi(res) ! rho value at this surface
      if(i_ox.eq.1)then !YuP[2020-07-29] Modified procedure, in case of i_ox=1
        !Check rho at ray element: is it larger than rho_plasm_vac ?
        ! rho_plasm_vac is rho at plasma-vacuum border;
        ! normally, it would be at rho=1.0, but if edge density is high,
        ! Xe=1 conversion layer could be just outside of rho=1.0
        ! Then, we define rho_plasm_vac as max(rhoconv,1.0)+0.1
        !  (or slightly different, see sub.bound)
        ! where rhoconv is rho at Xe=1.0 conversion layer,
        ! calculated in sub.ox_conversion_grill_in_poloidal_point.
        if(rhoconv.lt. 1.d0)then
          rho_plasm_vac=1.d0
        else ! rhoconv=1.d0 or larger
          rho_plasm_vac= rhoconv+0.1d0
        endif
        !rho_plasm_vac= max(rhoconv,1.0)+0.5d0 ! Another attempt
        if (rho.ge. rho_plasm_vac-epsbnd) then
         ibound=1 !here ray element is outside of rho_plasm_vac
         goto 10 !-> make a "reflection" of N, 
         !(then shoot rays in vacuum, in other subr.)
        else
         !All checks are done, none gave ibound --> +1
         goto 30 !return/end  ray element is inside plasma
        end if
      else ! i_ox=2 (original procedure, it was used for i_ox=1 also)
        if (res.ge.(psilim*(1.d0-epsbnd))) then
         ibound=1 !here psi>psilim (LCFS)
         !write(*,*)'in bound res,psilim,epsbnd',res,psilim,epsbnd
         goto 10 !-> procedure to make a reflection/stop
         !YuP: looks like this check is not needed, since few lines below
         !a check is made for rho=1.0 crossing 
        end if
        !Check rho at ray element: is it larger than 1.0 ?
        !bmod=b(z,r,phi) !YuP: Why needed here?
        if (rho.ge.1.d0-epsbnd) then
         ibound=1 !here ray element is outside of rho=1.0
         !write(*,*)'in bound rho,epsbnd',rho,epsbnd
         goto 10 !-> procedure to make a reflection/stop
        else
         !All checks are done, none gave ibound --> +1
         goto 30 !return/end  ray element is inside plasma
        end if
      endif ! i_ox
c-----------------------

 10   continue 
      !-> proceed to make a reflection/stop

cyup      write(*,*)'in bound no_reflection,i',no_reflection,'ibound',ibound
      
      !YuP: no_reflection=1 means artificial reflection is OFF.
      !The next section will check - is ray crossing the chamber?
      !It seems that no_reflection setting has nothing to do with this.
      !Maybe drop (no_reflection.eq.1) from if()then ?
      !YuP: For i_ox.eq.1, skip this part anyway. 
      !The M0 point (see oxb_launch.pdf) is supposed to be close to LCFS,
      !and will be detected by rho>1.0 below.
      if(i_ox.ne.1)then !YuP[2020-07-29] Added condition
      if((no_reflection.eq.1).and.(ibound.eq.1))then 
          !YuP: no_reflection=0 means artificial reflection is ON
          ibound=-1 !checking wall_limiter   !initialization
cyup          write(*,*)'before wall_limiter_reflection_point'
cyup          write(*,*)'z,r,phi,cnz,cnr,cm',
cyup     &               z,r,phi,cnz,cnr,cm

          call  wall_limiter_reflection_point(z,r,phi,cnz,cnr,cm,
     &    i_reflection, cnzref,cnrref,cmref,z_ref,r_ref,phi_ref)
          ! i_reflection becomes 1 only when old and new points
          ! of the ray are on opposite sides of chamber wall;
          ! otherwise i_reflection stays equal to 0.

cyup          write(*,*)'after wall_limiter_reflection_point i_reflection=',
cyup     &    i_reflection

          if (i_reflection.eq.1) then ! Reflection happened
cyup             write(*,*)'z,r,phi,cnz,cnr,cm',
cyup     &                  z,r,phi,cnz,cnr,cm
cyup             write(*,*)'z_ref,r_ref,phi_ref',z_ref,r_ref,phi_ref
cyup             write(*,*)'cnzref,cnrref,cmref',
cyup     &                  cnzref,cnrref,cmref
             irefl=irefl+1
             ibound=1 !checked wall_limiter: reached 
             iflref=1 ! iflref=1 after reflection; =-1 before reflection
          endif
          goto 30 !return/end, skip the rest (skip call_reflect)
      endif ! (no_reflection.eq.1)&(ibound.eq.1)  
      endif !(i_ox.ne.1) !YuP[2020-07-29] Added condition

c-------------------------------------------------------
 100  continue ! YuP: no handle ?
       
      nx4=nx+4
      ncx=nx4
      ncy=ny4
      dpsidr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z,nx4a)
      dpsidz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z,nx4a)
      edg=dpsidz*dzdt+dpsidr*drdt
      edg= edg*prmt(3)/dabs(prmt(3))
c      write(*,*)'in bound dpsidr,dpsidz',dpsidr,dpsidz
c      write(*,*)'in bound drdt,dzdt',drdt,dzdt,'edg',edg
      grad=dsqrt(dpsidz*dpsidz+dpsidr*dpsidr)
      dmrdt=dsqrt(drdt*drdt+dzdt*dzdt)
      edag=edg/(grad*dmrdt)
      if (edg.gt.0d0) then !Pos. deriv. of PSI (normal, since psilim>psimag)
         irefl=irefl+1
         write(*,*)'in bound: edg>0. edg,grad,dmrdt',edg,grad,dmrdt
         write(*,*)'in bound r[m],z[m], irefl=', r,z,irefl
         call reflect(z,r,phi,cnz,cnr,cnzref,cnrref)
         iflref=1 ! Just after call_reflect: reflection is done
      else ! negative grad(PSI)  no reflection needed; continue
         grad=dsqrt(dpsidz*dpsidz+dpsidr*dpsidr)
         dmrdt=dsqrt(drdt*drdt+dzdt*dzdt)
         write(*,*)'in bound: edg<0. edg,grad,dmrdt=',edg,grad,dmrdt
      endif
      
 30   continue

      return
      end subroutine bound

c======================================================================
c======================================================================
c        *********************BOUNDC correction****************
c        *                        -                           *
c        * This function gives  the ray position in           *
c        * the plasma. If the point of the ray is inside the  *
c        * limiting flux surface, leave iBOUNDc=-1, 
c        * else set iBOUNDc=1                                 *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        Input parameters					   !
c      								   !
c      Z, R, - point which the components of magnetic field	   !
c                  are calculated in.
c
c     the input parameter epsbnd(the small distance from the boundary)
c     is set inside this subroutine
c------------------------------------------------------------------!
c        Output parameter:iboundc=1 point is outside the  plasma   !				   !
c                         iboundc=-1 point is inside the  plasma   !				   !
c------------------------------------------------------------------!
      subroutine boundc(z,r,iboundc)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'oxb.i' ! YuP[2020-07-29] needed now for i_ox=1 run.
      !oxb.i now stores the value of rhoconv; 
      ! rhoconv is rho at Xe=1.0 conversion layer,
      ! calculated in sub.ox_conversion_grill_in_poloidal_point.
      ! It also stores:
      ! rho_plasm_vac is rho at plasma-vacuum border;
      ! normally, it would be at rho=1.0, but if edge density is high,
                    ! Xe=1 layer could be just outside of rho=1.0
      ! Then, we define rho_plasm_vac as max(rhoconv,1.0)+0.1
      !  (or slightly different, see sub.bound)
      ! The value of rho_plasm_vac is used in sub.bound in i_ox=1 run.

      double precision
     1 ias1r,ias2r,ias2r_Sm
      !epsbnd=1.d-7
      !epsbnd=1.d-8
      !epsbnd=1.d-10  ![m]!
      epsbnd=1.d-5 ![m] it must be equal ebsbnd in bound and boundc
      iboundc=-1

c--------------------------------------YuP[2020-07-23] added more checks
      ! Check: if the ray is out of R,Z-grid:
      if ( (r.lt.xeqmin+epsbnd).or.(r.ge.xeqmax-epsbnd) .or.
     +     (z.lt.zeqmin+epsbnd).or.(z.ge.zeqmax-epsbnd) ) then
         write(*,*)'boundc: out of grid. r,z,epsbnd=',r,z,epsbnd
         iboundc=1
         goto 10
      end if
c--------------------------------------------------------------YuP added

cSAP080727
      !YuP: For i_ox.eq.1, skip the check related to no_reflection. 
      !The M0 point (see oxb_launch.pdf) is supposed to be close to LCFS,
      !and will be detected by rho>1.0 in subr.bound().
      if(i_ox.ne.1)then !YuP[2020-07-29] Added condition
        if (no_reflection.eq.1) goto 10 !->return/end
        !YuP: =1 means artificial reflection is OFF
      endif !(i_ox.ne.1) !YuP[2020-07-29]

c-----------------------
      if ((r.lt.rmin+epsbnd).or.(r.gt.rmax-epsbnd)) then
cyup         write(*,*)'in boundc rmin+epsbnd,r,rmax-epsbnd',
cyup     1   rmin+epsbnd,r,rmax-epsbnd
c	 read(*,*)
         iboundc=1
         goto 10
      end if
      rrr=r
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      idx=0
      ipx=ip
      ipx4=ip+4
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
      zm=zzrm
c      write(*,*)'in boundc rrr,zm,z,zp',rrr,zm,z,zp
      if ((z.gt.zp-epsbnd).or.(z.lt.zm+epsbnd)) then
cyup         write(*,*)'in boundc rrr,zm,z,zp',rrr,zm,z,zp
cyup         write(*,*)'in boundc zm+epsbnd,z,zp-epsbnd',
cyup     1	 zm+epsbnd,z,zp-epsbnd
         iboundc=1
         goto 10
      end if
c---------------------------------------------------------------------
cSm0302243
      nx4=nx+4
      ny4=ny+4
      ncx=nx4
      ncy=ny4
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nx4a)
      rho=rhopsi(res)
      
      if(i_ox.eq.1)then !YuP[2020-07-29] Modified procedure, in case of i_ox=1
      
        !Check rho at ray element: is it larger than rho_plasm_vac ?
        ! rho_plasm_vac is rho at plasma-vacuum border;
        ! normally, it would be at rho=1.0, but if edge density is high,
        ! Xe=1 conversion layer could be just outside of rho=1.0
        ! Then, we define rho_plasm_vac as max(rhoconv,1.0)+0.1
        !  (or slightly different, see sub.bound)
        ! where rhoconv is rho at Xe=1.0 conversion layer,
        ! calculated in sub.ox_conversion_grill_in_poloidal_point.
        !if(rhoconv.lt. 1.d0)then
        !  rho_plasm_vac=1.d0
        !else ! rhoconv=1.d0 or larger
        !  rho_plasm_vac= rhoconv+0.1d0 !already set in sub.bound
        !endif
        if (rho.ge. rho_plasm_vac-epsbnd) then
         ibound=1 !here ray element is outside of rho_plasm_vac
         goto 10 
        end if
        
      else ! i_ox=2 (original procedure, it was used for i_ox=1 also)

        if (res.gt.(psilim*(1.d0-epsbnd))) then
         !YuP: This check looks excessive, since we are also checking rho below.
         iboundc=1
         goto 10
        end if
        if (rho.gt.1.d0-epsbnd) then
         iboundc=1 ! rho>1.0
         goto 10
        end if
        
      endif ! i_ox
c---------------------------------------------------------------------
 10     continue

      return
      end subroutine boundc

c        *********************REFLECT**************************
c        *                        -                           *
c        * This subroutine reflects                           *
c        * the ray from plasma boundary	                      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        Input parameters					   !
c      								   !
c      Z, R, PHI - point which the components of magnetic field	   !
c                  are calculated in.				   !
c      cnz,cnr						   !
c------------------------------------------------------------------!
c     output parameters:					   !
c	       	cnzref,cnrref= N_z,N_r after reflectin point	   !
c----------------------------------------------------------------------
      subroutine reflect(z,r,phi,cnz,cnr,cnzref,cnrref)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      double precision
     1 ias1r,ias2r,ias2r_Sm

cSm030224
              nx4=nx+4
              ny4=ny+4
	      ncx=nx4
	      ncy=ny4
c              dpsidz=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z)
c              dpsidr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z)
              dpsidz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z,nx4a)
              dpsidr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z,nx4a)

	      dpsimd=dsqrt(dpsidz**2+dpsidr**2)
	      pp=1.0d0/dpsimd
	      cnpsi=-pp*(dpsidz*cnz+dpsidr*cnr)
	      cnteta= pp*(-dpsidz*cnr+dpsidr*cnz)
c	      write(*,*)'before reflect cnpsi,cnteta',cnpsi,cnteta
	      cnzref=pp*(cnpsi*dpsidz+cnteta*dpsidr)
	      cnrref=pp*(cnpsi*dpsidr-cnteta*dpsidz)
cSAP081010 
cyup              write(*,*)'reflect1 cnzref,cnrref',cnzref,cnrref
              grad_psi_z=dpsidz*pp
              grad_psi_r=dpsidr*pp
cyup              write(*,*)'grad_psi_z,grad_psi_r',grad_psi_z,grad_psi_r
cyup              write(*,*)'cnpsi,cnz,cnr',cnpsi,cnz,cnr
             
              cnz_psi=-grad_psi_z*cnpsi
              cnr_psi=-grad_psi_r*cnpsi
cyup              write(*,*)'cnz_psi,cnr_psi',cnz_psi,cnr_psi

              cnzref=cnz-2.d0*cnz_psi
              cnrref=cnr-2.d0*cnr_psi
cyup              write(*,*)'reflect2 cnzref,cnrref',cnzref,cnrref

	      cnpsi=-pp*(dpsidz*cnzref+dpsidr*cnrref)
	      cnteta=pp*(-dpsidz*cnrref+dpsidr*cnzref)
c	      write(*,*)'after	reflect	cnzref,cnrref,cnpsi,cnteta',
c     1	      cnzref,cnrref,cnpsi,cnteta
cyup          write(*,*)'reflection point'

      return
      end
c        *********************BOUNTest*************************
c        *                        -                           *
c        * This function controles the ray output from        *
c        * the plasma(for test)	        		      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        Input parameters					   !
c      Z, R, PHI - point which the components of magnetic field	   !
c                  are calculated in.				   !
c        Output parameters:it prints r,z,zm,zp 			   !
c------------------------------------------------------------------!
      subroutine boundt(z,r)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      double precision
     1 ias1r,ias2r,ias2r_Sm
c-----------------------
      write(*,*)'in boundt ,r,z',r,z
      write(*,*)'in boundt rmin,rmax',rmin,rmax
      read(*,*)
      rrr=r
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      idx=0
      ipx=ip
      ipx4=ip+4
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
      zm=zzrm
      write(*,*)'in boundt r,zm,z,zp',r,zm,z,zp
      read(*,*)
c-----------------------
cSm030224
      nx4=nx+4
      ny4=ny+4

      ncx=nx4
      ncy=ny4
c      res=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nx4a)
      write(*,*)'in boundt res,psilim',res,psilim
      read(*,*)
      rho=rhopsi(res)
c-----------------------
      bmod=b(z,r,phi)
      write(*,*)'in boundt rho',rho
      read(*,*)
      return
      end






c        *********************zpzmlim**************************
c        *                        -                           *
c        * This function calculates the limiter               *
c        * zp(r) and rm(r)				      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        Input parameters					   !
c      	 r,							   |
c        zml(im),rml(im),zpl(ip),rpl(ip)                               |
c------------------------------------------------------------------!
c     output parameters:zp,zm
c-------------------------------------------------------------------
      subroutine zpzmlim(r,zp,zm)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'five.i'
      include 'limit.i'

cyup      write(*,*)'in zpzmlim r,ip,im',r,ip,im

      imm=1
      ipp=1
      if(r.gt.rpl(ip)+1.d-13)then
        WRITE(*,*)'in zpzmlim r.gt.rpl(ip)'
        WRITE(*,*)'r,rpl(ip)',r,rpl(ip)
        STOP ' bound.f in zpzmlim'
      endif

      if(r.lt.rpl(1)-1.d-13)then
        WRITE(*,*)'in zpzmlim r.lt.rpl(1)'
        WRITE(*,*)'r,rpl(1)',r,rpl(1)
        STOP ' bound.f in zpzmlim'
      endif
      do i=1,ip-1
c         write(*,*)'i,r,rpl(i),rpl(i+1)',i,r,rpl(i),rpl(i+1)
         if ((rpl(i).le.r).and.(r.le.rpl(i+1))) then
            ipp=i
c	    write(*,*)'ipp',ipp
	    goto 20
         endif
      enddo
 20   continue

      do i=1,im-1
c         write(*,*)'i,r,rml(i),rml(i+1)',i,r,rml(i),rml(i+1)
         if ((rml(i).le.r).and.(r.le.rml(i+1))) then
            imm=i
c	    write(*,*)'imm',imm
	    goto 30
         endif
      enddo
 30   continue

c      write(*,*)'20,ipp,imm',ipp,imm
      delrp=rpl(ipp+1)-rpl(ipp)
      if (delrp.ne.0.d0) then
         zp=zpl(ipp)+(zpl(ipp+1)-zpl(ipp))*(r-rpl(ipp))/delrp
      else
	 zp=zpl(ipp)
      endif

      delrm=rml(imm+1)-rml(imm)
      if (delrm.ne.0.d0) then
c         write(*,*)'imm,zml(imm),zml(imm+1)',imm,zml(imm),zml(imm+1)
c	 write(*,*)'rml(imm),r,rml(imm+1),delrm',
c     1	 rml(imm),r,rml(imm+1),delrm
         zm=zml(imm)+(zml(imm+1)-zml(imm))*(r-rml(imm))/delrm
      else
	 zm=zml(imm)
      endif
c      write(*,*)'zp,zm',zp,zm
      return
      end






  
