c        ********************** efield1***********************
c        * this subroutine calculates the components of  the *
c        * wave electric field using complex dielectric      *
c        * tensor -  reps(3,3),    complex N_perp - cnperp   *
c        * and real N_par- cnpar                             *
c        *****************************************************
c
!--------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      cnpar     real*8 N_parallel                        !
c      cnperp    complex*16 N_perp      (YuP: renamed cnx to cnperp)                        !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c								    !
c        output parameters					    !
c     ex,ey,ez				      			    !
c     complex*16 cex,cey,cez in common/cefield/		    !
c     eplus, eminus, epar - wave polarization 			    !
!--------------------------------------------------------------------
      subroutine efield1(cnpar,cnperp,ex,ey,ez,eplus,eminus,epar)
      implicit none !double precision (a-h,o-z)
      include 'param.i'
      include 'one.i' ! contains ioxm, etc.
      include 'eps.i' ! contains double complex reps(3,3)
      include 'cefield.i' ! Store cex,cey,cez
      complex*16 dlamb(3,3),aa(3,3),ccnz,cnperp,ci,cp,det,ccn2
      complex*16 cdelz ! for tests only
      real*8 ex,ey,ez,cnpar,eplus,eminus,epar,cpa !YuP[2020-03-10] added
      real*8 cnper1,cnpar2,cnprim,cnper2,dcn2,dmaxaa,ppp,cmodee,dmode
      real*8 root2, wka(3,3), pp1,pp2,pp3, pp01,pp02
      integer i0,j0,i,j, j00(3), i01,i02,i03
  
c      write(*,*)'in efield1'
      
      ci=dcmplx(0.0d0,1.0d0)
      cnpar2=cnpar*cnpar
      cnper1=dreal(cnperp)
      cnprim=dimag(cnperp)
      cnper2=cnper1*cnper1+cnprim*cnprim
      dcn2=cnpar2+cnper2
c --------------------------------------------------------
      !write(*,*)'1--- efield cnpar,cnper1',cnpar,cnper1
      
c      write(*,'(a,6e12.3)')
c     +   'reps(1,1:3)=',reps(1,1),reps(1,2),reps(1,3)
c      write(*,'(a,6e12.3)')
c     +   'reps(2,1:3)=',reps(2,1),reps(2,2),reps(2,3)
c      write(*,'(a,6e12.3)')
c     +   'reps(3,1:3)=',reps(3,1),reps(3,2),reps(3,3)

      ccnz=dcmplx(cnpar,0.d0)
      ccn2= ccnz*ccnz+cnperp*cnperp ! =N^2

      ! For the 1st row of ((-D_tensor)).(E_vector)=0, see Eqn(1.6) in manual
      dlamb(1,1)=ccnz*ccnz-reps(1,1)     ! -eps_xx + Npar^2 
      dlamb(1,2)=-reps(1,2)              ! -eps_xy
      dlamb(1,3)=-ccnz*cnperp-reps(1,3)  ! -eps_xz - Npar*Nperp
      
      ! For the 2nd row of ((-D_tensor)).(E_vector)=0
      dlamb(2,1)=-reps(2,1)              ! -eps_yx
      dlamb(2,2)=ccn2-reps(2,2)          ! -eps_yy + N^2  !differ id=6 vs id=2
      dlamb(2,3)=-reps(2,3)              ! -eps_yz
      
      ! For the 3rd row of ((-D_tensor)).(E_vector)=0
      dlamb(3,1)=-ccnz*cnperp-reps(3,1)  ! -eps_zx - Npar*Nperp
      dlamb(3,2)=-reps(3,2)              ! -eps_zy
      dlamb(3,3)=cnperp*cnperp-reps(3,3) ! -eps_zz + Nperp^2 !differ id=6 vs id=2
      !write(*,*)'--- efield: cnperp=',cnperp
      !write(*,'(a,4e12.3)')'dlamb(22),dlamb(33)=',dlamb(2,2),dlamb(3,3)

!-------------------------------------------------------------
c     cofactors (see Genray manual, Ch.8)
      aa(1,1)= dlamb(2,2)*dlamb(3,3)-dlamb(2,3)*dlamb(3,2)
      aa(1,2)=-dlamb(2,1)*dlamb(3,3)+dlamb(2,3)*dlamb(3,1)
      aa(1,3)= dlamb(2,1)*dlamb(3,2)-dlamb(2,2)*dlamb(3,1)
      
      aa(2,1)=-dlamb(1,2)*dlamb(3,3)+dlamb(3,2)*dlamb(1,3)
      aa(2,2)= dlamb(1,1)*dlamb(3,3)-dlamb(1,3)*dlamb(3,1)
      aa(2,3)=-dlamb(1,1)*dlamb(3,2)+dlamb(1,2)*dlamb(3,1)

      aa(3,1)= dlamb(1,2)*dlamb(2,3)-dlamb(2,2)*dlamb(1,3)
      aa(3,2)=-dlamb(1,1)*dlamb(2,3)+dlamb(2,1)*dlamb(1,3)
      aa(3,3)= dlamb(1,1)*dlamb(2,2)-dlamb(1,2)*dlamb(2,1)
c      write(*,'(a,6e12.3)')
c     +   'aa(1,1:3)=',aa(1,1),aa(1,2),aa(1,3)
c      write(*,'(a,6e12.3)')
c     +   'aa(2,1:3)=',aa(2,1),aa(2,2),aa(2,3)
c      write(*,'(a,6e12.3)')
c     +   'aa(3,1:3)=',aa(3,1),aa(3,2),aa(3,3)

!-------------------------------------------------------------
c choice of max(abs(aa(i,j))
      !YuP[2020-08-21] Re-designed this part - on selection of row (i0)
      ! in aa(i,j) matrix that will be used 
      ! for calculation of Erf components.
      ! The problem occurs in a very-low density plasma.
      ! In this case, eps_xx~1, eps_yy~1, eps_zz~1, eps_xy~0,
      ! and then the three rows of ((D)).E=0 equations
      ! (see manual, Eqns. 1.4-1.6) become:
      !   (1-Npar^2)*Ex=0
      !   (1-N^2)*Ey=0
      !   (1-Nperp^2)*Ez=0
      ! If Npar<<Nperp, 
      ! then Ex=0, and we have to select which of the remaining
      ! two equations to use. Practically, the original procedure 
      ! of detecting the row i0 in aa(i,j) matrix 
      ! that contains the largest value of |aa| yields two rows 
      ! (2nd or 3rd) with approximately same value of |aa|.
      ! As a result, the solution jumps from 2nd to 3rd row and back,
      ! giving either predominantly Ey or predominantly Ez polarization.
      ! In view of this confusion, the procedure of selecting 
      ! the proper row in aa matrix is revised.
      ! Now we find TWO rows in aa matrix that contain TWO largest
      ! elements of |aa(i,j)|. Then, we check - if these two elements
      ! are about same then the selection of row is based on ioxm value.
      
      ! To get the original functionality, simply use a small number in
      ! in front of (pp01+pp02), (e.g., 0.001) in the line below:
      ! if( abs(pp01-pp02).lt. 0.001*(pp01+pp02) ) then
      
      wka(:,:)=abs(aa) !To be searched for largest element in each row.
      j00=MAXLOC(wka,DIM=2) ! j-index of max value in each row 
      pp1=wka(1,j00(1)) ! max value among abs(aa(1,:))  [1st row]
      pp2=wka(2,j00(2)) ! max value among abs(aa(2,:))  [2nd row]
      pp3=wka(3,j00(3)) ! max value among abs(aa(3,:))  [3rd row]
      ! Sort the three values in descending order,
      ! Record the i-indices (row number) of the two largest elements
      if(pp1.gt.pp2)then
        if(pp1.gt.pp3)then
          if(pp2.gt.pp3)then
            ! The order is {pp1,pp2,pp3}
            i01=1 ! 1st largest
            i02=2 ! 2nd largest
          else ! pp3>pp2
            ! The order is {pp1,pp3,pp2}
            i01=1 ! 1st largest
            i02=3 ! 2nd largest      
          endif
        else ! pp3>pp1
          ! The order is {pp3,pp1,pp2}
            i01=3 ! 1st largest
            i02=1 ! 2nd largest      
        endif
      else ! pp2>pp1
        if(pp2.gt.pp3)then
          if(pp1.gt.pp3)then
            ! The order is {pp2,pp1,pp3}
            i01=2 ! 1st largest
            i02=1 ! 2nd largest      
          else ! pp3>pp1
            ! The order is {pp2,pp3,pp1}
            i01=2 ! 1st largest
            i02=3 ! 2nd largest      
          endif
        else !pp3>pp2
          ! The order is {pp3,pp2,pp1}
            i01=3 ! 1st largest
            i02=2 ! 2nd largest      
        endif
      endif
      ! So, the two largest elements are 
      pp01= wka(i01,j00(i01))
      pp02= wka(i02,j00(i02))
      i0=i01 ! by default, select the row with the largest element.
      ! BUT: The two values (pp01 and pp02)
      ! can be very close to each other.
      ! Check - if they differ by ~~50% or less, 
      ! make selection of the row based on ioxm
      if( abs(pp01-pp02).lt. 0.001*(pp01+pp02) ) then
       !Note from test: 
       !  Using 0.5 factor in  if( abs(pp01-pp02).lt.0.5*(pp01+pp02) )
       !  still gives occasional jumps between different polarizations.
       !  Factor 1.0 seems to be optimal. 
       ![YuP: actually not; depends on particular run]
       if(ioxm.eq.1)then ! Select 3rd row - to yield Ez~1 polarization
         if(i01.eq.3) i0=i01
         if(i02.eq.3) i0=i02
       else ! ioxm=-1;   ! Select 2nd row - to yield Ey~1 polarization
         if(i01.eq.2) i0=i01
         if(i02.eq.2) i0=i02
       endif
      endif
      !--------- DONE: i0 row is selected ---------------------------
      !write(*,'(a,3e12.3)')'pp1,pp2,pp3=', pp1,pp2,pp3
      !write(*,'(a,2i5)')   'i01,i02=', i01,i02
      !write(*,'(a,2e12.3)')'wka for i01,i02 ', pp01,pp02

!YuP: old version for selection of i0      
!      i0=1
!      j0=1
!      dmaxaa=0.d00
!      do 30 i=1,3
!         do 40 j=1,3
!            ppp=abs(aa(i,j))
!            if(ppp.gt.dmaxaa) then
!              dmaxaa=ppp
!              i0=i
!              j0=j
!            end if
!40       continue
!30    continue
!YuP: old version end

!---------------------------------------------------------------
c      write(*,'(a,6e12.3)')
c     +   'dlamb(1,1:3)=',dlamb(1,1),dlamb(1,2),dlamb(1,3)
c      write(*,'(a,6e12.3)')
c     +   'dlamb(2,1:3)=',dlamb(2,1),dlamb(2,2),dlamb(2,3)
c      write(*,'(a,6e12.3)')
c     +   'dlamb(3,1:3)=',dlamb(3,1),dlamb(3,2),dlamb(3,3)
      
      cex=aa(i0,1)
      cey=aa(i0,2)
      cez=aa(i0,3) 
      
! cez==Ez is along B-field line.
c
c      p=1.d0/dsqrt(cdabs(cex)*cdabs(cex)+cdabs(cey)*cdabs(cey)+
c     1   cdabs(cez)*cdabs(cez))
c      cex=cex*p
c      cey=cey*p
c      cez=cez*p

c      cdelz=cdabs(cez)/cez

!Sm021021
c      thetz=datan2(dimag(cez),dabs(dreal(cez)))    
c      if (dreal(cez).ge.0.d0) then
c         p=1.d0
c      else
c         p=-1.d0
c      endif     
c      thetz=thetz*p
c      cdelz=dcmplx(dcos(thetz),-dsin(thetz))
         
!Sm021021
c      write(*,*)'cex',cex
c      write(*,*)'cey',cey
c      write(*,*)'cez',cez
c      write(*,*)'dimag(cez),dreal(cez)',dimag(cez),dreal(cez)
c       thetz=datan2(dimag(cez),dreal(cez))
c      write(*,*)'1 thetz',thetz
c      write(*,*)'dimag(cez),cdabs(cez)',dimag(cez),cdabs(cez)
c       thetz=dasin(dimag(cez)/cdabs(cez))
c      write(*,*)'2 thetz',thetz

c       cdelz=dcmplx(dcos(thetz),-dsin(thetz))

c       if (dreal(cez).eq.0.d0) goto 10
        if(abs(cez).gt.0.d0)then !YuP[2020-08-21] added this check/condition
          cp=dcmplx(dreal(cez),-dimag(cez))/abs(cez)
c        if (dreal(cez).lt.0d0) cp=-cp
c       cp=dcmplx(dreal(cez),-dimag(cez))/abs(cez)**2*
c     & dreal(cez)/abs(cez)
         cex=cex*cp
         cey=cey*cp
         cez=cez*cp
       else
         ! YuP: do nothing?
         WRITE(*,*)'efield1: abs(cez)=0.  cez=', cez
         !pause
       endif
c       write(*,*)'efield1 cez',cez
c       write(*,*)'efield1 cex',cex 
c       write(*,*)'efield1 cey',cey

 10    continue !YuP: not used now
c       cdelz=1.d0
c      if (j0.eq.1) cdelz=cdabs(cex)/cex
c      if (j0.eq.2) cdelz=cdabs(cey)/cey
c      if (j0.eq.3) cdelz=cdabs(cez)/cez
!-----------------------------------------------------------
c      cp=cdelz/dsqrt(cdabs(cex)*cdabs(cex)+cdabs(cey)*cdabs(cey)+
c     1   cdabs(cez)*cdabs(cez))

      !YuP[2020-03-10] Changed "cp" (which is a complex value) to "cpa" (real*8)
      cpa=1.d0/dsqrt(abs(cex)*abs(cex)+abs(cey)*abs(cey)+
     &   abs(cez)*abs(cez))    
      !write(*,*)'efield1: cex=',cex,'cey=',cey,'cez=',cez, ' cpa=',cpa
      !Now normalize to have Ex/|E|, etc.
      cex=cex*cpa
      cey=cey*cpa
      cez=cez*cpa
c      write(*,*)'efield1 cez',cez
c      write(*,*)'efield1 cex',cex 
c      write(*,*)'efield1 cey',cey
!      cez1=dsqrt(1.d0-cda bs(cex)**2-abs(cey)**2)
      ex=abs(cex)
      ey=abs(cey)
      ez=abs(cez)
      cmodee=ex*ex+ey*ey+ez*ez
      !write(*,'(a,6e12.3,i3)')'cex,cey,cez; i0=',cex,cey,cez,i0
c      write(*,*)'ex=',ex,'ey=',ey,'ez=',ez,'cmodee',cmodee
c      write(*,*)'thetz',thetz,'cdelz',cdelz,'cpa=',cpa
c      write(*,*)'efield1 cex=',cex,'cey=',cey,'cez=',cez,'cez1=',cez1
c      write(*,*)'efield1 cez,ez',cez,ez
!---------------------------------------------------------------
      root2=dsqrt(2.d0)
      eplus=abs(cex+ci*cey)/root2
      eminus=abs(cex-ci*cey)/root2
      dmode=eplus**2+eminus**2+ez**2
c      write(*,*)'eplus=',eplus,'eminus=',eminus,'epar=',epar
c     1 ,'dmode=',dmode
      epar=abs(cez*cnpar+cex*dreal(cnperp))
     &    /dsqrt(cnpar**2+dreal(cnperp)**2)
c      write(*,*)'efield1 cez,cex,cey',cez,cex,cey
c      write(*,*)'efield1 cnpar',cnpar
c      write(*,*)'efield1 cnperp',cnperp
c      write(*,*)'efield1 epar',epar

      return
      end subroutine efield1











