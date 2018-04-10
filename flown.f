c        ********************** flown  ***********************
c        *                      ------                       *
c        * this subroutine calculates the wave energy flow   *
c        * for normalized electric field modeE=1             *
c        *****************************************************
c
c--------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of the point where the electric field !
c                 is calculated.      		         	    !
c								    !
c      cnz, cnr, cm - n_z, n_r, r*n_phi - components of  wave  ref- !
c                     ractive index. 				    !
c       							    !
c      complex components of dielectric tensor reps(i,j) are 	    !
c      in common block 'eps'.These components are created in        !
c      subroutine hamilt1					    !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c								    !
c        output parameters					    !
c								    !
c        cflown                               			    !
c        cflown=B~(i)*B(i)+E~(i)*d(omega*eps(i,j))/domega*E(j)	    !
c        B~ =conjg(B), E~=conjg(E)				    !
c--------------------------------------------------------------------
      subroutine flown (z,r,phi,cnz,cnr,cm,
     1                   cflown, dwepsdw)
c     cex,cey,cez - complex electr.field polarization E_x/E,E_y/E,E_z/E
c     from common bloc /cefield/
c     cbx,cby,cbz - complex magnet.field polarization B=N*E modeE=1

      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'eps.i'
      include 'cefield.i'
      dimension vp(nbulka),wp(nbulka)
      double complex cbx,cby,cbz
      double complex cepsp,cepsm,dwepsdw,ce,cec,cb,cbc,cflown
      dimension cepsp(3,3),cepsm(3,3),dwepsdw(3,3)
      dimension ce(3),cec(3),cb(3),cbc(3)
cSAP081122
      complex*16 ceps_herm(3,3) !Hermitian part of reps
c--------------------------------------------------------------
c     for test comparison grpde2
c     bmod2,emod2
c--------------------------------------------------------------
c     ce - wave electric field
c     cec- comlex conjugate wave electric field
c     cb - wave magnetic field
c     cbc- comlex conjugate wave magnetic field
c---------------------------------------------------------------
      bmod=b(z,r,phi)
      gam1=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam1)
      dc=dcos(gam1)
      dcn=cn(r,cnz,cnr,cm)
      cnpar=dcn*dc
      cnper=dcn*ds
c      write(*,*)'in flown dcn,cnpar,cnper',dcn,cnpar,cnper
c      write(*,*)'cex,cey,cez',cex,cey,cez
c----------------------------------------------------------
c    complex components of the wave magnetic field

      cbx=-cnpar*cey
      cby=-cnper*cez+cnpar*cex
      cbz=cnper*cey
c      write(*,*)'cbx,cby,cbz',cbx,cby,cbz
c----------------------------------------------------------
      cb(1)=cbx
      cb(2)=cby
      cb(3)=cbz
      ce(1)=cex
      ce(2)=cey
      ce(3)=cez
c      write(*,*)'cb(i)',cb(1),cb(2),cb(3)
c      write(*,*)'ce(i)',ce(1),ce(2),ce(3)
c-----------------------------------------------------------
      step=0.0000001d0
      step=1.d-7
      hw=step*1.0d0
cSm990727
      hw=hw*frqncy
      pi=3.1415926d0

      hfrqnc=hw
      w0=frqncy
c-----------------------------------------------------------
      gam=gamma1(z,r,phi,cnz,cnr,cm)
      do 1 i=1,nbulk
         vp(i)=v(i)
         wp(i)=w(i)
 1    continue
      frqncpl=frqncy+hfrqnc
      df=frqncy/frqncpl
      do 2 i=1,nbulk
         v(i)=vp(i)*df*df
         w(i)=wp(i)*df
 2    continue

c************************************************
      cnrplus=cnr*df
      cnzplus=cnz*df
      cmplus=cm*df
c      write(*,*)'flown before hp=hamilt1'
      hp=hamilt1(z,r,phi,cnzplus,cnrplus,cmplus)
c      write(*,*)'flown after hp=hamilt1=',hp
c*************************************************
      w0p=w0+hw
cSAP081122
c      do 3 i=1,3
c        do 3 j=1,3
c           cepsp(i,j)=w0p*reps(i,j)
c 3    continue

cSAP081122 calculates hermitian part of the dielectric tensor reps
      call hermitian_part(reps,ceps_herm)
      do 3 i=1,3
        do 3 j=1,3
           cepsp(i,j)=w0p*ceps_herm(i,j)
 3    continue
c----------------------------------------------------------
      frqncmn=frqncy-hfrqnc
      df=frqncy/frqncmn
      do 4 i=1,nbulk
        v(i)=vp(i)*df*df
        w(i)=wp(i)*df
 4    continue
c************************************************
      cnrminus=cnr*df
      cnzminus=cnz*df
      cmminus=cm*df
      hm=hamilt1(z,r,phi,cnzminus,cnrminus,cmminus)
c*************************************************
      w0m=w0-hw
cSAP081122
c      do 5 i=1,3
c         do 5 j=1,3
c           cepsm(i,j)=w0m*reps(i,j)
c          write(*,*)'after hm i,j,cepsm(i,j)',i,j,cepsm(i,j)
c 5    continue
       call hermitian_part(reps,ceps_herm)

cSAP081122 calculates hermitian part of the dielectric tensor reps
      do 5 i=1,3
         do 5 j=1,3
           cepsm(i,j)=w0m*ceps_herm(i,j)
c          write(*,*)'after hm i,j,cepsm(i,j)',i,j,cepsm(i,j)
 5    continue
c-----------------------------------------------------------
      do 6 i=1,nbulk
        v(i)=vp(i)
        w(i)=wp(i)
 6    continue
c-----------------------------------------------------------
      do 7 i=1,3
         do 7 j=1,3
           dwepsdw(i,j)=(cepsp(i,j)-cepsm(i,j))/(2.d0*hw)
c          write(*,*)'i,j,dwepsdw(i,j)',i,j,dwepsdw(i,j)
 7    continue
      do 8 i=1,3
c        do 8 j=1,3
           cec(i)=dconjg(ce(i))
 	   cbc(i)=dconjg(cb(i))
c          write(*,*)'i,cec(i),cbc(i)',i,cec(i),cbc(i)
 8    continue

      cflown=dcmplx(0.0d00,0.0d00)
c----------------------------
c      bmod2=0.0d0
c----------------------------

      do 9 i=1,3
c         bmod2=bmod2+cb(i)*cbc(i)          
         cflown=cflown+cb(i)*cbc(i)
         do 9 j=1,3
 	    cflown=cflown+cec(i)*dwepsdw(i,j)*ce(j)
c          write(*,*)'j,ce(i),cec(j),dwepsdw(i,j)',
c     1               j,ce(i),cec(j),dwepsdw(i,j)
 9    continue
c      write(*,*)'end flown cflown=',cflown
c      write(*,*)'end flown'

      return
      end





c      if (id.eq.6) then !hot non-relativistic plasma
cc--------calculate the flux using the analytical derivatives 
cc        from the hermitian non-relativistic hot plasma dielectric tensor
cc        d(eps_herm^)/d(omega)=
cc        Sum{k=1,nbulk}[d(eps_herm^)/dXk * dXk/domega+d(eps_herm^)/dYk * dYk/domega]
cc        +d(eps_herm^)/dN_perp * dN_perp/domega 
cc        +d(eps_herm^)/dN_parallel * dN_parallel/domega      
cc          
c               
c         do i=1,3
c            do j=1,3
c               depsdomeg(i,j)=0.d0  
c            enddo
c         endoo

c         do i=1,3
c            do j=1,3
c               depsdomeg(i,j)=depsdomeg(i,j)
c            enddo
c         endoo


      subroutine flux_tokman(z,r,phi,cnz,cnr,cm,flux)
c-----calculates the wave energy flux using Tokman formula
c     S^=c/(8*pi)*dRe(lambda)/dN^ 
c
c     
c     flux^=(8*pi)S^/c !flux normalized at light speed
c
c     Now it works for id=10 case

c-----input    
c     z, r phi - coordinates of the point       
c								    
c     cnz, cnr, cm - n_z, n_r, r*n_phi - components of  wave 
c                    refractive index. 		

c-----output
c     flux(1)=flux_z,flux(2)=flux_r, flux(3)=flux_phi
c     
      implicit none
c-----input
      real*8 z,r,phi,cnz,cnr,cm
c-----output
      real*8 flux(3)

c-----local
      real*8 pi,p
      real*8 t,u(6),deru(6)

      pi=4.d0*datan(1.d0)
      p=1.d0/(8.d0*pi)
      p=0.5d0
      

      u(1)=z
      u(2)=r
      u(3)=phi
      u(4)=cnz
      u(5)=cnr
      u(6)=cm

c      deru(i) are right  hand  sides  of  geometrical  optics       !
c      equations							    !
c      n_z     - deru(1) =dD/dNz                                            !
c      n_r     - deru(2) =dD/dNr                                            !
c      r*n_phi - deru(3) =dD/dM                                           !
c      z       - deru(4) =dD/dz                                            !
c      r       - deru(5) =dD/dr                                           !
c      phi     - deru(6) =dD/dphi 

      call dddrz1(t,u,deru)

      flux(1)=p*deru(1) !z component
      flux(2)=p*deru(2) !r
      flux(3)=p*deru(3) !phi

      return
      end

      subroutine hermitian_part(K,K_herm)
c-----calculate Hermition part K_herm(3,3) of the complex tensor
c     K(3,3)

      implicit none

c-----input
      complex*16 K(3,3) !complex tensor
c-----output
      complex*16 K_herm(3,3) !hermitian part of K

       K_herm(1,1) = 0.5D0 * ( K(1,1) + dconjg(K(1,1) )) !kxx
       K_herm(2,2) = 0.5D0 * ( K(2,2) + dconjg(K(2,2) )) !kyy
       K_herm(3,3) = 0.5D0 * ( K(3,3) + dconjg(K(3,3) )) !kxx
       K_herm(1,2) = 0.5D0 * ( K(1,2) + dconjg(K(2,1) )) !kxy
       K_herm(1,3) = 0.5D0 * ( K(1,3) + dconjg(K(3,1) )) !kxz
       K_herm(2,3) = 0.5D0 * ( K(2,3) + dconjg(K(3,2) )) !kxz
       K_herm(2,1) = 0.5D0 * ( K(2,1) + dconjg(K(1,2) )) !kyx
       K_herm(3,1) = 0.5D0 * ( K(3,1) + dconjg(K(1,3) )) !kzx
       K_herm(3,2) = 0.5D0 * ( K(3,2) + dconjg(K(2,3) )) !kxz

       return
       end
