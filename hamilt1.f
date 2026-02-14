
c        ********************** hamilt1 **********************
c        *                      ------                      *
c        * this function calculates the hamiltonian of the  *
c        * the system of geometrical optics equations       *
c        * It will calculate the dielectric tensor
c        * reps(3,3) and will put this tensor to eps.i file
c        ****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of the point where the hamiltonian is !
c                 calculated.      		         	    !
c								    !
c      cnz, cnr, cm - n_z, n_r, r*n_phi - components of  wave  ref- !
c                     ractive index at this point.                  !
c      the angle gam between refractive index n and magnetic field  !
c      from common 'one'					    !
c-------------------------------------------------------------------
      double precision
     1function hamilt1(z,r,phi,cnz,cnr,cm)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'eps.i'
      include 'ions.i'
      double precision dshkarof
      double complex dhot,dhot_sum
      double complex ceps(3,3),hamiltc
      external x,y,tempe,dshkarof
      external dhot_sum

      double precision x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     .tpop_ar(nbulka),vflow_ar(nbulka)
      
      double complex K(3,3),dK_dx_ar(3,3,nbulka),dK_dy_ar(3,3,nbulka),
     &dK_dt_ar(3,3,nbulka),dK_dnper(3,3),dK_dnpar(3,3)


      double complex compl_nper !for Eric tensor
      double complex eps_weiss(3,3)

c-----for i_disp=16 Bonoli dispersion
c      real*8,  dimension(1:nbulk) :: d_D_0_d_x_ar
c      real*8,  dimension(1:nbulk) :: d_D_0_d_y_ar
c      real*8,  dimension(1:nbulk) :: d_D_0_d_v_thermal_ar
      real*8   D_0   
c      real*8   d_D_0_d_n_par,d_D_0_d_n_perp 
c-----external
      double complex det       

      bmod=b(z,r,phi)
           
      ds=dsin(gam)
      dc=dcos(gam)
      ds2=ds*ds
      dc2=dc*dc
      if (((id.eq.1).or.(id.eq.2)).or.(id.eq.3)) then
        call tensrcld(z,r,phi)
      end if
c---------------------------------------------------------
c
c     Appleton-Hartree dispersion relation
c---------------------------------------------------------
      if (id.eq.3) then
         ds4=ds2*ds2
         xi=x(z,r,phi,1)
         yi=y(z,r,phi,1)
         py2=yi*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
         sqrdet=dsqrt(py4*ds4+4.*py2*px2*dc2)
         pz=2.d0*px-py2*ds2+ioxm*sqrdet
         cnt=cn(r,cnz,cnr,cm)
         cn2=cnt*cnt
         hamilt1=cn2-(1.d0-2.d0*xi*px/pz)
         goto 10
      end if

c--------------------------------------------------------
      cnt=cn(r,cnz,cnr,cm)
      if ((id.eq.1).or.(id.eq.2)) then
c--------------------------------------------------------
c       id=1 or id=2 cold plasma dispersion relation
       	cn2=cnt*cnt
       	cn4=cn2*cn2
c--------------------------------------------------------
        call abc(z,r,phi,ds2,dc2,ad,bd,cd)
        d4=ad
       	d2=bd
       	d0=cd
c--------------------------------------------------------
        if (id.eq.1) then
           hamilt1=d4*cn4+d2*cn2+d0 ! a*N^4 +b*N^2 +c
        end if

        if (id.eq.2) then ! See Eq.(4.12) in Genray manual
         ! For id=1 or 2,   
         ! Note: a=A*delta, b=B*delta, c=C*delta (where delta=1-Y)
         ! and sqrt(det)= sqrt(B*B-4*A*C) * |delta|
         ! in   N^2 = (-b +ioxm*sqrt(b*b-4*a*c))/(2a)
           ibmx=min(ib,nbulk) ! safety check: not to exceed nbulk
           delib= 1.d0-y(z,r,phi,ibmx) ! 1-Y (electrons or ions)
           sign_del=1.d0 !sign(1.d0,delib)
           oxm= ioxm*sign_del ! replaced ioxm by oxm below
        
           hamilt1=cn2+(d2-oxm*dsqrt(d2*d2-4.d0*d4*d0))/
     *                (2.d00*d4)
c  YuP[07-2017] However, not always ioxm is used to select roots.
c  Sometimes ioxm_n_npar is used.
c  Need to fix?  Maybe not - the matching ioxm
c  [ such that N(N_par,ioxm_n_npar)=N(gam,ioxm) ]
c  should be found in cninit12/nper_npar_ioxm_n_npar
        end if
        go to 10
      end if ! id=2
c     end cold plasma dispersion relation
c-----------------------------------------------------------
c     det=dsqrt(d2*d2-4.d00*d4*d0)
c     cn2od=(-d2+det)/(2.d00*d4)
c     cn2ex=(-d2-det)/(2.d00*d4)
c     write(*,*)'cn2od=',cn2od,'cn2ex=',cn2ex
c     cn2c=cnz*cnz+cnr*cnr+cm*cm/(r*r)
c     write(*,*)'cn2c=',cn2c
c-----------------------------------------------------------
c     id=4 dispersion relation from mazzucato code
c     for hermitian part of the dielectric tensor(ihermloc=1)
      if (id.eq.4) then
         cnpar=cnt*dc
         cnper=cnt*ds
         ihermloc=1
         call hamiltmuz(cnpar,ihermloc,z,r,phi,cnper,hamltmuz)
         hamilt1=hamltmuz
         go to 10
      end if
c **  end if mazzucato
c-----------------------------------------------------------
c     id=5 dispersion relation from mazzucato code fomega_maz
c     for full the dielectric tensor(iherm=2)
      if (id.eq.5) then
         hamilt1=fomega_maz(cnz,cnr,cm,z,r,phi)
         go to 10
      end if
c **  end if id=5
c-------------------------------------------------------------
c     Hot non-relativistic plasma
      if (id.eq.6) then
         cnpar=cnt*dc
         cnper=cnt*ds
cSAP080102
c         write(*,*)'hamilt1.f id=6 cnpar,cnper',cnpar,cnper
c         write(*,*)'z,r,phi',z,r,phi

         do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpop_zrp(z,r,phi,i) !YuP[2024-08-14] was tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
         enddo

         hamiltc=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .   vflow_ar,cnpar,cnper,1,reps)

         hamilt1=dreal(hamiltc)
c         write(*,*)'hamilt1.f hamilt1=',hamilt1
         go to 10
      end if
c **  end if forest
c-------------------------------------------------------------
c     Shkarofsky code relativistic, electron plasma
      if (id.eq.7) then
         cnpar=cnt*dc
         cnper=cnt*ds
         xe=x(z,r,phi,1)
	 ye=y(z,r,phi,1)
	 te=tempe(z,r,phi,1) ! kev
         friq=frqncy
         hamilt1=dshkarof(xe,te,ye,cnpar,cnper,friq,hami)        

         go to 10
      end if
c **  end if shkarofsky
c-------------------------------------------------------------
c     Ono tensor for fast waves 
      if (id.eq.8) then
         cnt=cn(r,cnz,cnr,cm)
         cnpar=cnt*dc
         cnper=cnt*ds

         do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
         enddo

         call ono_tens(dmas,X_ar,Y_ar,T_av_ar,nbulk,
     &   cnpar,cnper,1,0,
cSm030426
c     &   K,hamiltc,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)
     &   reps,hamiltc,dK_dx_ar,dK_dy_ar,dK_dt_ar,dK_dnper,dK_dnpar)
        
         hamilt1=dreal(hamiltc)
         go to 10
      end if
c **  end if Ono dispersion fast waves 
c-------------------------------------------------------------
c     id=9 dispersion relation from full(non-hermitian) hot plasma tensor 
c     for full hot plasma the dielectric tensor(iherm=2)
      if (id.eq.9) then
c         write(*,*)'hamilt1 before fomega_hot'
         hamilt1=fomega_hot(cnz,cnr,cm,z,r,phi)
c         write(*,*)'hamilt1 after fomega_hot'
         go to 10
      end if
c **  end if id=9
c-----------------------------------------------------------

c     id=10 dispersion function
c     D=real part(eigenvalue(k) of the dispersion tensor)
c     Dispersion tensor uses the Mazzucato dielectric tensor.
      if (id.eq.10) then
c         write(*,*) 'hamilt1 id=10'    
         cnpar=cnt*dc
         cnper=cnt*ds
c--------k_root is the number of the cubic equation root: lambda.
c        k_root was set at the plasma edge.

c--------calculate the dispersion function: hamilt1=Real(hamiltc)
c        Here hamiltc=eigenvalue(k_root).
c        It will calculate the Mazzucato non-Hermitian 
c        complex dielectric tensor reps(3,3). 
c        This tensor will be in eps.i

         cnprim=0.d0
         call hamilt_eigenvalue_muz(z,r,phi,cnpar,cnper,cnprim,k_root,
     &   hamilt1,hamiltc) 

         go to 10
      end if
c **  end if id=10
c-----------------------------------------------------------

c     id=11 dispersion function
c     D=real part(D relativistic dispersion function)
c     from Eric Nelson_Melby dielectric tensor

      if (id.eq.11) then
c         write(*,*) 'hamilt1 id=11'    
         cnpar=cnt*dc
         cnper=cnt*ds
         x_ar(1)=x(z,r,phi,1)
         y_ar(1)=y(z,r,phi,1)    !question for electron -y? 
         te=tempe(z,r,phi,1) !kev

         compl_nper=dcmplx(cnper,0.d0)

c         call Disp_Nelson_Melby(te,cnpar,x_ar(1),y_ar(1),
c     +   compl_nper,eps_weiss,hamiltc) !eps_weiss is in Weiss system
c                                       !k is in z-y plane
c
c         call eps_weiss_to_stix(eps_weiss,K) !trasform eps_weiss to
                                             !K in Stix system
         call Disp_Nelson_Melby(te,cnpar,x_ar(1),y_ar(1),
     +   compl_nper,K,hamiltc) ! K is in z-y plane
                               ! in Stix coordinate system

         if (iherm.eq.1) then
c-----------calcualte hermitian part reps of the complex tensor K
            call herm(K,reps)
            hamiltc=det(reps,cnpar,compl_nper)   
         endif
         hamilt1=dreal(hamiltc)
         go to 10
      end if
c **  end if id=11
c-----------------------------------------------------------
c     id=12 or id=15 dispersion function
c     D=real part(eigenvalue(k) of the relativistic dispersion tensor)
c     for id=12: Dispersion tensor uses Eric Nelson_Melby dielectric tensor
c     for id=15: Dispersion tensor uses Abhay Ram dielectric tensor
c       these use the Westerhof-Tokman method with the eigenvalues
      if (id.eq.12 .or. id.eq.15) then

c         write(*,*) 'hamilt1 id=12'    

         cnpar=cnt*dc
         cnper=cnt*ds 
      
c--------k_root is the number of the cubic equation root: lambda.
c        k_root was set at the plasma edge.

c--------calculate the dispersion function: hamilt1=Real(hamiltc)
c        Here hamiltc=eigenvalue(k_root).
c        It will calculate the Eric Nelson_Melby non-Hermitian - old version
c        complex dielectric tensor reps(3,3). 
c        This tensor will be in eps.i

         cnprim=0.d0
c         call hamilt_eigenvalue_eric(z,r,phi,cnpar,cnper,cnprim,k_root,
c     &   hamilt1,hamiltc) 

c        It will calculate either the  Eric Nelson_Melby or Abhay Ram
c         non-Hermitian relativistic
c        complex dielectric tensor reps(3,3). 
c        This tensor will be in eps.i

cSAP091104
c         write(*,*)'hamilt1.f  hamilt1 before
c     & hamilt_eigenvalue_combined z,r,phi,cnpar,cnper,cnprim',
c     & z,r,phi,cnz,cnr,cm

         call hamilt_eigenvalue_combined(z,r,phi,cnpar,cnper,cnprim,
     &   k_root,hamilt1,hamiltc) 

cSAP091104
c        write(*,*)'hamilt1.f  hamilt1 after hamilt_eigenvalue_combined'

         go to 10
      end if
c **  end if id=12 or id=15
c---------------------

c     id=13 dispersion function
c     D=real part(eigenvalue(k) of the dispersion tensor.
c     Dispersion tensor uses Stix Hot-plasma dielectric tensor
 
      if (id.eq.13) then
c         write(*,*) 'hamilt1 id=13'    

         cnpar=cnt*dc
         cnper=cnt*ds 
      
c--------k_root is the number of the cubic equation root: lambda.
c        k_root was set at the plasma edge.

c--------calculate the dispersion function: hamilt1=Real(hamiltc)
c        Here hamiltc=eigenvalue(k_root).
c        It will calculate the Hot plasma non-Hermitian 
c        complex dielectric tensor reps(3,3). 
c        This tensor will be in eps.i

         cnprim=0.d0
         call hamilt_eigenvalue_hot(z,r,phi,cnpar,cnper,cnprim,k_root,
     &   hamilt1,hamiltc) 

         go to 10
      end if
c **  end if id=13
c-----------------------------------------------------------
c     id=14 dispersion function
c     D=real part(D relativistic dispersion function)
c     from Eric Nelson_Melby dielectric tensor ,if npar.le. 0.38D0
c     from Ram Abhay dielectric tensor ,if npar.ge. 0.38D0
C 1Sep2005 -- Now Disp_combined is entirely Abhay's tensor, which 
C works just as well or faster than Nelson-Melby's, when the resolution
C is lower than the maximum like it used to be.
C
C ENM 15Mar2006 -- After finding that the combined version which jumps from
C     one to another dispersion relation depending on n_parallel didn't work
C     very well, and finding that the Ram version works well, even for fairly
C     small n_parallel, as long as you adjust the resolution parameters (see
C     genray.dat template), now id=14 is just the Ram tensor (id=11 is just the
C     Nelson-Melby tensor)
C
      if (id.eq.14) then
c        write(*,*) 'hamilt1 id=14'    
         cnpar=cnt*dc
         cnper=cnt*ds
c        write(*,*)'hamilt1 cnt,dc,ds,cnpar,cnper',cnt,dc,ds,cnpar,cnper
     
         x_ar(1)=x(z,r,phi,1)
         y_ar(1)=y(z,r,phi,1)    !question for electron -y? 
         te=tempe(z,r,phi,1) !kev

         compl_nper=dcmplx(cnper,0.d0)

c         call Disp_combined(te,cnpar,x_ar(1),y_ar(1),
c     +   compl_nper,eps_weiss,hamiltc) !eps_weiss is in Weiss system
c                                       !k is in z-y plane
c
c         call eps_weiss_to_stix(eps_weiss,K) !trasform eps_weiss to
                                             !K in Stix system

         call Disp_Ram(te,cnpar,x_ar(1),y_ar(1),
     +   compl_nper,K,hamiltc) !K is in z-y plane
                               !in Stix coordinates

c         write(*,*)'hamilt1 id=14 te,cnpar,x_ar(1),y_ar(1),compl_nper',
c     &   te,cnpar,x_ar(1),y_ar(1),compl_nper
c         write(*,*)'hamilt1 id=14 hamiltc', hamiltc


         if (iherm.eq.1) then
c-----------calcualte hermitian part reps of the complex tensor K
            call herm(K,reps)
            hamiltc=det(reps,cnpar,compl_nper)   
         endif
         hamilt1=dreal(hamiltc)
         go to 10
      end if
c **  end if id=14
c------------------------------------------------------------
c      Paul T. Bonoli and Ronald C. Englade dispersion function
c     for LH wave: cold plasma plus thermal correction     
      if (id.eq.16) then

         cnt=cn(r,cnz,cnr,cm)
         cnpar=cnt*dc
         cnper=cnt*ds

         do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
         enddo

c        write(*,*)'hamilt1 before  bonoli_Dispersion'

         call  bonoli_Dispersion(
     &   X_ar,Y_ar,T_av_ar,dmas,nbulk,
     &   cnpar,cnper,D_0)

c        write(*,*)'hamilt1 after  bonoli_Dispersion'

c         call  bonoli_Dispersion(
c         call bonoli_D_der_d_x_d_y_d_v_t_d_n_par_d_n_perp(
c     &   dmas,X_ar,Y_ar,T_av_ar,nbulk,
c     &   cnpar,cnper,0,
c     &   D_0,
c     &   d_D_0_d_x_ar,d_D_0_d_y_ar,
c     &   d_D_0_d_v_thermal_ar,
c     &   d_D_0_d_n_par,
c     &   d_D_0_d_n_perp) 

         hamilt1 = D_0

      endif
c**  end if id=16

 10   continue
      return
      end
