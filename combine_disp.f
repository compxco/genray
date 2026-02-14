 
      Subroutine Disp_combined(T_e,N_parallel,X_e,Y_e,N_perp,eps,D)

      DOUBLE PRECISION T_e,N_parallel,X_e,Y_e   !T_e in[KeV]
      DOUBLE COMPLEX   N_perp,eps(1:3,1:3),D
      External Disp_Nelson_Melby, Disp_Ram
      include 'param.i'
      include 'one.i'

C***** Disp_Nelson_Melby uses transformed I. Weiss method, which is faster
C**** than Disp_Ram, which uses Trubnikov method, but not quite as accurate
C**** and can be more susceptible to numerical errors, especially when 
C**** Im N_perp is more than about half of Re N_perp.
C**** Trubnikov method speeds up around N_parallel of 0.3 or 0.4.
C
C**** Weiss method does not work for N_parallel > 1.0, so Trubnikov method
C**** is necessary for N_parallel > 1.0

C******** Found that Ram's method works actually better than Nelson-Melby
C****** when using lower resolution integration parameters
C**** Still nearly as accurate, and faster.  However, Nelson-Melby is still
C**** necessary when n_parallel is very small (near zero). 
C
c      write(*,*)'Disp_combined N_parallel',N_parallel
C      If (dabs(N_parallel) .le. 0.38D0) then
c         write(*,*)'Disp_combined before Nelson'
C        call Disp_Nelson_Melby(T_e,N_parallel,X_e,Y_e,N_perp,eps,D)
C      else
c         write(*,*)'Disp_combined before Ram'
C
C Rather than swapping dispersion relations for a given n_parallel (didn't work
C  too well at boundary, and for small n_parallel, usually Nelson-Melby version
C  works for entire run), call appropriate dispersion relation, according to
C  id:  id=11, |D^Hermitian| relativistic Nelson-Melby
C       id=12, eigenvalue of complex D, relativistic Nelson-Melby
C       id=14, |D^Hermitian| relativistic Ram
C       id=15, eigenvalue of complex D, relativistic Ram

      if (id.eq.11 .or. id.eq.12) then 

cSAP091104
c        write(*,*)'combine_disp.f in Disp_combined before
c     &   Disp_Nelson_Melby T_e,N_parallel,X_e,Y_e,N_perp,eps',
c     &   T_e,N_parallel,X_e,Y_e,N_perp,eps

        call Disp_Nelson_Melby(T_e,N_parallel,X_e,Y_e,N_perp,eps,D)
      else
        call Disp_Ram(T_e,N_parallel,X_e,Y_e,N_perp,eps,D)
      endif
C      endif

      end
C------------------------------------------------
c      include "Abhay_disp.f"
c      include "Eric_disp.f"



      double complex function relativistic_nperp(z,r,phi,
     &nll,cnteta,cnphi,K,iraystop)
c     calculates n_perp( and the relayvistic dielectric tensor K)
c     from the relativistic electron
c     dispresion function (combined Eric Nelson-Melby - Abhay Ram)
c     if ibw=0 it calculates the O and X modes , using the cold plasma roots
c              as the initial approximation
c     if ibw=1 it calculates BW root N_perp_bw> N_perp_O and N_perp_X
c--------------------------------------------------------------------
c     INPUTS:    
c     z,r,phi - space cordinates
c     nll     - parallel index of refraction n
c     cnteta  _ N component parallel to the poloidal direction
c     cnphi    _ N component parallel to the toroidal direction
c    
c     OUTPUTS:
c       function hotnperp=N_perpendiculae (double complex)
c       K double complex dielectric tensor
c       iraystop=0 the root have been found,
c       iraystop=1 the root was not found
c-----------------------------------------------------    
      !implicit double precision (a-h,o-z) 
cSm060725
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'nperpcom.i'

c-----input
c      integer nbulk 
      double precision z,r,phi
      double precision nll ! N_parallel
cSm060725
      double precision cnteta,cnphi
c-----external 
      double complex hotnp
      double precision b,x,y,tempe,tpoprho,vflowrho
      double precision ddwrap_relativistic, rtbis

cBH090202:  After pathscale compiler error:
      external ddwrap_relativistic
      external b,x,y,tempe,tpoprho,vflowrho
      real*8 tpop_zrp   !external func
      external tpop_zrp !external func
 
c-----local  
     
      double precision mass_ar(nbulka),x_ar(nbulka),y_ar(nbulka),
     .t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka)
      integer j,id_old

      double precision cnper,cnprim,
     &cnper_o,cnper_x,
cSm060725
     &cnpar2,cntang2,cm,cnz,cnr,cnper_02,cnper_0,xacc,cnperp_o,cnperp_x,
     &cnpermax,x1,x2
  
      integer ioxm_old,ihermloc,ioptmaz,
cSm060725
     &ntry
 
      double precision te_kev
cSm060725
      double complex compl_nper
      double complex d
      double complex eps_weiss(3,3) !eps_weiss is in Weiss system
                                    !k is in z-y plane
      
c-----output      
      double complex K(3,3)
      integer iraystop

      cnpar2=nll*nll
      cntang2=cnteta*cnteta+cnphi*cnphi
      
      if(ibw.eq.0) then ! X and O mode cnper calculation 
c-------calculation of perpendicular refractive index cnper
c       for the given ioxm=(+/-)1 mode  using Mazzucato solver
c       ioxm is in common /one/

        ihermloc=1
        ioptmaz=1
        call cnpermuz(nll,ihermloc,z,r,phi,cnper,cnprim,ioptmaz)
        write(*,*)'in function relativistic_nperp cnper =',cnper
        if(cnper.ge.0.d0) then
           iraystop=0
        else
           iraystop=1
           write(*,*)'in function relativistic_nperp cnpermuz could not
     &     find the root for mode ioxm=',ioxm
           id_old=id
           id=1            ! to use the cold plasma dispersion
   
           call cninit12(z,r,phi,nll,cnteta,cnphi,
     &                   cnz,cnr,cm,iraystop)
           id=id_old 
           cnper_02=cnz*cnz+cnr*cnr+(cm/r)**2-nll*nll
           write(*,*)'cold plasma cnper_02=',cnper_02
           cnper_0=dsqrt(cnper_02)
           WRITE(*,*)'cold plasma cnper_0=',cnper_0  
           STOP 'in function relativistic_nperp'
        endif
      else
        ! calculation of the EBW
c-------initialization common/nperpcom/, it is in nperpcom.i
        nllc=nll
        nbulkc=nbulk
        if(nbulka.lt.nbulk) then
          WRITE(*,*)'in forest.f in hotnperp nbulka.lt.nbulk'          
          WRITE(*,*)'nbulka=',nbulka,'nbulk=',nbulk
	  WRITE(*,*)'change parameter nbulka in file param.i'  
          STOP
        endif

c       write(*,*)'hotnperp before loop j nbulk=',nbuklk

        do j=1,nbulk
          massc(j)=dmas(j)
          xc(j)=x(z,r,phi,j)
          yc(j)=y(z,r,phi,j)
          if(j.eq.1) yc(1)=-yc(1) ! negative Y=(omega_ce/omega)
                                  ! for electrons
          tec(j)=tempe(z,r,phi,j)*1.d+3 !(eV) averaged temperature
          tpopc(j)=tpop_zrp(z,r,phi,j) !YuP[2024-08-14] wastpoprho(rho,j)
          vflowc(j)=vflowrho(rho,j)         
        enddo

        xacc=1.d-14  ! accuracy of the root calculations
        ntry=200


c-------calculate X and O modes using Mazucato solver 
c       to find the boundary for EBW root calculations 
        ihermloc=1
        ioptmaz=1
        ioxm_old=ioxm

c-------O mode calculations 
        ioxm=1
        call cnpermuz(nll,ihermloc,z,r,phi,cnper,cnprim,ioptmaz)
        write(*,*)'in function relativistic_nperp ioxm,cnper',ioxm,cnper
         
        if(cnper.gt.0.d0) then
           cnperp_o=cnper
        else !cnper =<0
           cnperp_o=0.d0
c----------check the cnper root using cold plasma: cnper_0
           write(*,*)'in function relativistic_nperp cnpermuz could not
     &     find the root for mode ioxm=',ioxm
           id_old=id
           id=1            ! to use the cold plasma dispersion        
           call cninit12(z,r,phi,nll,cnteta,cnphi,
     &                   cnz,cnr,cm,iraystop)
           id=id_old 
           cnper_02=cnz*cnz+cnr*cnr+(cm/r)**2-nll*nll
           write(*,*)'cold plasma cnper_02=',cnper_02
           cnper_0=dsqrt(cnper_02)
           write(*,*)'cold plasma cnper_0=',cnper_0  
        endif

c-------X mode calculations   
        ioxm=-1
        call cnpermuz(nll,ihermloc,z,r,phi,cnper,cnprim,ioptmaz)
        write(*,*)'in function relativistic_nperp ioxm,cnper',ioxm,cnper
         
        if(cnper.gt.0.d0) then
           cnperp_x=cnper
        else !cnper =<0
           cnperp_x=0.d0
c----------check the cnper root using cold plasma: cnper_0
           write(*,*)'in function relativistic_nperp cnpermuz could not
     &     find the root for mode ioxm=',ioxm
           id_old=id
           id=1            ! to use the cold plasma dispersion        
           call cninit12(z,r,phi,nll,cnteta,cnphi,
     &                   cnz,cnr,cm,iraystop)
           id=id_old 
           cnper_02=cnz*cnz+cnr*cnr+(cm/r)**2-nll*nll
           write(*,*)'cold plasma cnper_02=',cnper_02
           cnper_0=dsqrt(cnper_02)
           write(*,*)'cold plasma cnper_0=',cnper_0  
        endif  

        ioxm=ioxm_old

        cnpermax=dmax1(cnperp_x,cnperp_o)
        x1 = dlog(cnpermax+1.d-2)  !nperp
        x2 = x1+dlog(1.d0)         !nperp

c       expand the range of x2 until the first root is found
        iraystop=0
        do j = 1, ntry
          if (j.eq.ntry) then
            iraystop=1
            write(*,*)'relativistic_nperp could not find ebw root'
            return
          endif
ctest
c          p1=ddwrap_relativistic(x1)
c          p2=ddwrap_relativistic(x2)
c          write(*,*)'relativistic_nperp j,x2,x1',j,x2,x1
c          write(*,*)'ddwrap_relativistic(x1),ddwrap_relativistic(x2)',
c          p1,p2
cendtest
          if(ddwrap_relativistic(x2)*ddwrap_relativistic(x1).lt.0) then
            go to 12
          endif

          x1 = x2

          if(ddwrap_relativistic(x2)*ddwrap_relativistic(x1).ge.0) then
            x2 = x2 +dlog(2.0D0)
          endif   
        enddo
 12     continue

c        write(*,*)'relativistic_nperp before dexp(rtbis)'
        cnper=dexp(rtbis(ddwrap_relativistic,x1,x2,xacc))
        write(*,*)'in relativistic_nperp cnper=',cnper  

        compl_nper=dcmplx(cnper,0.d0)
        iraystop=0
        write(*,*)'in relativistic_nperp compl_nper=',compl_nper
        
        te_kev = tec(1)*1.d-3 !electron temperature [kev]
        write(*,*)'te_kev,nllc,xc(1),dabs(yc(1))',
     &  te_kev,nllc,xc(1),dabs(yc(1))

c        call Disp_combined(te_kev,nllc,xc(1),dabs(yc(1)),
c     +  compl_nper,eps_weiss,d) !eps_weiss is in Weiss system
c                                !k is in z-y plane        
c        call eps_weiss_to_stix(eps_weiss,K) !trasform eps_weiss to
c                                             !K in Stix system

        if (id.eq.11 .or. id.eq.12) then
          call Disp_Nelson_Melby(te_kev,nllc,xc(1),dabs(yc(1)),
     +  compl_nper,K,d) !K is in z-y plane 
                        ! in Stix coordinates
        else
          call Disp_Ram(te_kev,nllc,xc(1),dabs(yc(1)),
     +  compl_nper,K,d) !K is in z-y plane 
        endif

        write(*,*)'after Disp_combined d',d

      endif   !ibw=1 BW

      relativistic_nperp=cnper

      write(*,*)'in relativistic_nperp cnper=',cnper,
     &'iraystop=',iraystop

      return
      end

      double precision function ddwrap_relativistic(log_nper)
c-----calculates Real(relativistic combined
c     Eric Eric Nelson-Melby - Abhay Ram dispersion function)
c     INPUT:
c       log_nper=ln(Re(N_perpendicular))
c       the input data from common /nperpcom/
cSm060725
      Implicit none
    
c-----input
      double precision log_nper
c-----local
      double precision nperp
      double complex  K(3,3),d
      include 'param.i'
      include 'nperpcom.i'
CENM included one.i so that id would be visible
      include 'one.i'

      complex*16 compl_nper
      double precision te_kev

      nperp = dexp(log_nper)
      compl_nper=dcmplx(nperp,0.d0)

c-----data from include 'nperpcom.i'
      te_kev = tec(1)*1.d-3 !electron temperature [kev]
 
C      print *,'IN ddwrap_relativistic :, id=',id
      if (id.eq.11 .or. id.eq.12) then
        call Disp_Nelson_Melby(te_kev,nllc,xc(1),dabs(yc(1)),
     &compl_nper,K,d) !  K is in z-y plane
                      ! in Stix coordinates
      else
        call Disp_Ram(te_kev,nllc,xc(1),dabs(yc(1)),
     &compl_nper,K,d) !  K is in z-y plane
                      ! in Stix coordinates
      endif
      ddwrap_relativistic = dreal(d) 

      return
      end



      subroutine Im_nperp_Westerhof_Tokman(z,r,phi,
     &cnpar,cnper,cnprim_in,id,cnprim_out)
C*** For use with iabsorp=8 option *****
c-----calculate Im(N_perpendicular) and 
c     the electric field polarization E
c     using the Westerhof-Tokman form of dispersion function
c     for different tensors:
c     id=10 Mazzucato relativistic tensor
c     id=12 relativistic tensor: Eric Nelson-Melby (for n|| approx. 0)
c     id=13 hot plasma tensor
c     id=15 relativistic tensor: Abhay Ram
c
c-------------------------------------------------
c      input: 
c      z,r,phi are the space coordinates 
c      cnpar   is N_parallel. It is Real
c      cnper,cnprim_in are Re(N_perp) and Im(N_perp)  
c-------------------------------------------------
c      output:
c      cnprim_out  Im(N_perp)     
c-------------------------------------------------
      implicit none
      include 'eps.i' !contains reps(3,3) complex dielectric tensor

c-----input    
      real*8 z,r,phi,! space coordinates. 
c     Radius rho should be calculated before call hamilt_eigenvalue
c     in subroutine b()    
     &cnpar,cnper,cnprim_in!Re(N_parallel),Re(N_perpendicular),
                        !Im(N_perpendicular)
      integer id        !dispersion function type

c-----output
      real*8 cnprim_out  !Im(N_perp)
c-----externals
      real*8 b

c-----locals
      real*8 bmod,step,step_2,cnperp,cnperm
      real*8 hamilt,     !Re(eigenvalue)
     &hamilt_im,
     &hamiltp,hamiltm,d_hamilt_d_nperp

      complex*16 hamiltc !eigenvalue

      step=1.d-5 ! N_perp step for numerical derivatives

      cnprim_in=0.d0 ! Im(N_perp) initialization

c      if(i_im_nperp.eq.1) then
c--------calculate ImN_perp using the formula
c        cnprim_out=ImN_perp=abs(ImD_full/dD_hermitian/dReN_perp))
         
         bmod=b(z,r,phi)

cSm060719
c         call hamilt_nperp_npar_Westerhof_Tokman(z,r,phi,
c     &   cnpar,cnper,cnprim_in,id,hamilt,hamiltc)
        
c         hamilt_im=dimag(hamiltc)
         
         cnperp=cnper+step

         if (cnper.gt.step) then
            cnperm=cnper-step
            step_2=2.d0*step
         else
            cnperm=cnper
            step_2=step
         endif
         call hamilt_nperp_npar_Westerhof_Tokman(z,r,phi,
     &   cnpar,cnperp,cnprim_in,id,hamiltp,hamiltc)
        
         call hamilt_nperp_npar_Westerhof_Tokman(z,r,phi,
     &   cnpar,cnperm,cnprim_in,id,hamiltm,hamiltc)
       
         d_hamilt_d_nperp=(hamiltp-hamiltm)/step_2

cSm060719 calculate Im(D) and dielectric tensor reps (in eps.i)
         call hamilt_nperp_npar_Westerhof_Tokman(z,r,phi,
     &   cnpar,cnper,cnprim_in,id,hamilt,hamiltc)
        
         hamilt_im=dimag(hamiltc)
         
         cnprim_out=dabs(hamilt_im/ d_hamilt_d_nperp)
c      endif  

      return
      end


      subroutine hamilt_nperp_npar_Westerhof_Tokman(z,r,phi,
     & cnpar,cnper,cnprim,id,D,Dc)
c-----calculate hamiltonian as eigenvalue(N^) 
c     input: 
c     z,r,phi are the space coordinates 
c     cnpar   is N_parallel. It is Real
c     cnper,cnprim are Re(N_perp) and Im(N_perp)
c
c     id  : gives the dispersion function type 
c     id=10 Mazzucato relativistic tensor
c                   Westerhof - Tokman dispersion function
c     id=12 relativistic tensor: Eric Nelson-Melby (for n|| approx. 0)
c                   Westerhof - Tokman dispersion function
c     id=13 hot plasma tensor
c                   Westerhof - Tokman dispersion function
c     id=15 relativistic tensor: Abhay Ram
c                   Westerhof - Tokman dispersion function
c-----output:
c     D=real part(eigenvalue(k) of the dispersion tensor)
c     Dc=eigenvalue(k) of the dispersion tensor

      implicit none
      include 'eps.i' !contains reps(3,3) complex dielectric tensor

c-----input    
      real*8 z,r,phi,! space coordinates. 
c     Radius rho should be calculated before call hamilt_eigenvalue
c     in subroutine b()    
     &cnpar,cnper,cnprim!Re(N_parallel),Re(N_perpendicular),
                        !Im(N_perpendicular)
      integer id ! the dispersion function type 
                 ! used with Westehof Tokman procedure
                 ! id=10 Mazzucato relativistic tensor
                 ! id=12 relativistic tensor: Eric Nelson Melby (for small n||)
                 ! id=13 hot plasma tensor
                 ! id=15 relativistic tensor: Abhay Ram
c-----output
      real*8 D      !dispersion function=Re(eigenvalue(k))
      complex*16 Dc !Dc=eigenvalue(k_root), k_k_root=1,2,3

c-----locals
      integer k_root !the number of the eigenvalue, k_root=1,2,3

c-----check id value
      if ((id.ne.10).and.(id.ne.12).and.(id.ne.13).and.(id.ne.15)) then
        WRITE(*,*)'hamilt_nperp_npar_Westerhof_Tokman'
        WRITE(*,*)'it should be id=10,12,13,15'
        WRITE(*,*)'id=',id
        STOP 'in hamilt_nperp_npar_Westerhof_Tokman'
      endif

      k_root=1 !really k_root is not used in the following subroutines

      if (id.eq.10) then !Mazzucato tensor
         call hamilt_eigenvalue_muz(z,r,phi,cnpar,cnper,cnprim,
     &   k_root,D,Dc)
      endif

      if (id.eq.12 .or. id.eq.15) then !either relativistic
                         !Nelson Melby or Ram tensor
         call hamilt_eigenvalue_combined(z,r,phi,cnpar,cnper,cnprim,
     &   k_root,D,Dc)
      endif

      if (id.eq.13) then !hot plasma tensor
         call hamilt_eigenvalue_hot(z,r,phi,cnpar,cnper,cnprim,
     &   k_root,D,Dc)
      endif

      return
      end


c     *****  absorp_relativist_disp_combined  ***********
c 
c     *  this subroutine calculates imaginary part       
c     *  of refractive index using formula from:   	 
c     *  2k^_i_s*(P^+T^)=P_abs_s {Stix p.74 (17)}        
c     *  and relativistic tensor from Disp_combined         
c     *   05/03/18                                       
c     ****************************************************
c         input parameters: z,r,phi  (m,m,radian) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c
c                           It calculates the electric field polarization
c                           using the hot plasma tensor.The polarization
c                           (cex,cey,cez) will be in common /efield.i/ 
c------------------------------------------------------------------
      subroutine absorp_relativist_disp_combined(z,r,phi,cnpar,cnper,
     &cnprim_e)
      !implicit double precision (a-h,o-z)
cSm060725
      implicit none
      include 'param.i'
      include 'eps.i'
      include 'cefield.i'
      
c-----input 
      real*8  z,r,phi,cnper,cnpar
c-----output
      real*8 cnprim_e
    
c-----local
      real*8 X_e,Y_e,T_e,
     & cnprim,power_abs_e,omega,bmod,pi,
     & T_perp_r,d_reps_herm_d_n_perp,
     & ex,ey,ez,eplus,eminus,epar,pp,step,
     & Poynting_vector_perp_r

      complex*16 D_rel,! relativistic dispersion function Determinant
     & K_s(3,3),      ! susceptibilties of the given specie  
     & cnx,           ! complex n_perp
     & ci,
     & kappa_e(3,3), !anti-Hermitian part of eps for electron specie
     & T_perp_c,power_abs_e_c,
     & reps_p_c(3,3),reps_m_c(3,3),d_reps_herm_d_n_perp_c,
     & cbx,cby,cbz,ce(3),cec(3),cb(3),cbc(3),
     & Poynting_vector_x_c, c_cnper,cnper_p,cnper_m
C ENM 26jul06 -- added c_cnper to be a complex version of cnper (with 0 imaginary part)
C  because disp_combined expects a complex n_perp as input (Also moved cnper_p and cnper_m to the complex*16 list
      integer iherm_loc,i,j,j1,j2
 
c-----external
      real*8 b,x,y,tempe,tpoprho

      pi=4.d0*datan(1.d0) !YuP[2018-10-13][2020-01-27] BUG was: pi=4.d0*dtan(1.d0)
      c_cnper=dcmplx(cnper,0.0d0)
c----------------------------------------------
c     calculate relativisticplasma tensor reps(). 
c     reps() will be in common /eps/ in eps.i file.
c---------------------------------------------
      bmod=b(z,r,phi) ! small radius rho calculations.
                      ! rho will be in common /one.i/

     
      X_e=x(z,r,phi,1)
      Y_e=y(z,r,phi,1)
c     if(i.eq.1) y_ar(1)=-y_ar(1) ? question
      T_e=tempe(z,r,phi,1)        ! kev
             
      call Disp_combined(T_e,cnpar,X_e,Y_e,c_cnper,reps,D_rel)
    
c--------------------------------------------------------------
c     Calculate electric field polarization cex,cey,cez
c     using relativistic dielectric tensor reps()
c     and real N_perp=cnper.
c     The polarization will be in common /cefield/ in cefield.i file.
c-------------------------------------------------------------------   
      cnx=dcmplx(cnper,0.d0)
      call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
      ce(1)=cex
      ce(2)=cey   
      ce(3)=cez
c-------------------------------------------------------------------
c     complex components of the wave magnetic field
c
c     ce - wave electric field
c     cec- comlex conjugate wave electric field
c     cb - wave magnetic field
c     cbc- comlex conjugate wave magnetic field
c-------------------------------------------------------------------
      cbx=-cnpar*cey
      cby=-cnper*cez+cnpar*cex
      cbz=cnper*cey
 
      cb(1)=cbx
      cb(2)=cby
      cb(3)=cbz
      
      do j=1,3
         cec(j)=dconjg(ce(j))
 	 cbc(j)=dconjg(cb(j))
      enddo      
c------------------------------------------------------------------
c     Calculate apsorped power for all species
c    
c     Power=omega/(8Pi)[E~(i) . eps_a_herm(i,j) 'E~(j)] 
c------------------------------------------------------------------
c     The code calculates absorped powes for electron specie
c     using anti-hermitian part of the susceptibilty
c     The code uses normalized Power=Power/omega
c------------------------------------------------------------------
c     1. Calculate susceptibilities K_s() for all species.      
c     2. Calculate absorped power  power_abs_s(i) for
c        all "i=1,...,nbulk" species.      
c-----------------------------------------------------------------
      ci=dcmplx(0.d0,1.d0)

      cnprim=0.d0
c     omega=2.d0*pi*frqncy*1.d9       !frqncy in GHZ

     
      power_abs_e_c=0.d0

      do j1=1,3
         do j2=1,3
            kappa_e(j1,j2)=0.5d0*(reps(j1,j2)-dconjg(reps(j2,j1)))/ci
            power_abs_e_c=power_abs_e_c+cec(j1)*kappa_e(j1,j2)*ce(j2)
         enddo
      enddo

      power_abs_e_c=power_abs_e_c/(8.d0*pi)
      power_abs_e=dreal(power_abs_e_c)

      write(*,*)'absorp__relativist_disp_combined'
      write(*,*)'power_abs_e_c,power_abs_e',
     &          power_abs_e_c,power_abs_e
     
c------------------------------------------------------------------
c     calculations of the Poynting flux 
c
c                      c_light
c     Poynting_flux^=  ------- [ E^~ (vectr*) B^+ E^ (vectr*) B^~ ]
c                       16 Pi 
c-------------------------------------------------------------------
c     x direction of the Stix sytem is parallel to N_perp^.
c     So, x component of the Poynting flux is parallel to N_perp^.
c
c     Poynting_vector_x_c is complex x component of 
c     Poynting flux in the Stix coordinates. 
c------------------------------------------------------------
c     In the code it was used the normalized variable 
c     Poynting_vector_perp_r= Poynting /c_light
c-------------------------------------------------------------
      pp=1.d0/(16.d0*pi)
      Poynting_vector_x_c = pp*((cec(2)*cb(3)-cec(3)*cb(2))+
     &                          (ce(2)*cbc(3)-ce(3)*cbc(2)))     

      Poynting_vector_perp_r=dreal(Poynting_vector_x_c)
c-------------------------------------------------------------------
c     calculations T=T_vector the flux of nonelectromagnetic energy
c     Stix book p.74 (19)
c           omega        d                c_light       d
c     T^= - ----- E^~.----(eps_h^^).E^= - -------- E^~. ---(eps_h^^) .E^
c           16Pi         dk^              16Pi          dN^
c
c     For calculations of Im (K_perp) we need T^_perp.
c     T^perp it is parallel to Re( k^_perp)
c
c                c_light          d
c     T^_perp= - -------- E^~(i). -----------(eps_h(i,j) . E(j)
c                 16Pi            dRe(N_perp)              
c   
c-------------------------------------------------------------------
c     In the code it was used the normalized variable
c     T_perp_r = (T/c_light) 
c------------------------------------------------------------------
c      step=1.d-7
      step=1.d-3     

      cnper_p=c_cnper*(1.d0+step)
      cnper_m=c_cnper*(1.d0-step)
         
      call Disp_combined(T_e,cnpar,X_e,Y_e,cnper_p,reps_p_c,D_rel)

      call Disp_combined(T_e,cnpar,X_e,Y_e,cnper_m,reps_m_c,D_rel)
      
      T_perp_c=dcmplx(0.d0,0.d0)

      do j1=1,3
        do j2=1,3
          d_reps_herm_d_n_perp_c=(reps_p_c(j1,j2)-
     &                            reps_m_c(j1,j2))/(2.d0*step*cnper)  

          T_perp_c=T_perp_c+cec(j1)*d_reps_herm_d_n_perp_c*ce(j2)
          enddo
      enddo

      T_perp_r=-dreal(T_perp_c)/(16.d0*pi)
c      write(*,*)'T_perp_c,T_perp_r',T_perp_c,T_perp_r 
      
      cnprim_e=0.5d0*power_abs_e/(T_perp_r+Poynting_vector_perp_r)

      cnprim_e=dabs(cnprim_e)

      if(cnprim_e.lt.0.d0) then
        write(*,*)'!!!cnprim_e<0 cnprim_e', cnprim_e 
      endif
     
      write(*,*)'T_perp_r,Poynting_vector_perp_r,T/P',
     &T_perp_r,Poynting_vector_perp_r,T_perp_r/Poynting_vector_perp_r

      write(*,*)'cnprim_e',cnprim_e
         
c--------------------------------------------------------
c     calculated comples relativistic dielectric tensor reps
c     reps will be in common /eps.i/
c--------------------------------------------------------     
      call Disp_combined(T_e,cnpar,X_e,Y_e,c_cnper,reps,D_rel)
    

      return
      end
