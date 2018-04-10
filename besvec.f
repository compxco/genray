c
c
cend
ctop
c
c This is a package of subroutines needed for the calculation of ece
c Written by Bob Harvey and Martin O'Brien, summer 1986
c
      subroutine besvec(n,bes,besp,nbes,bmax,besdel)
c
c besvec calculates an array of Bessel Functions used for interpolation
c it calls besj
c
cSm001216
c      implicit double precision (a-h,o-z)
cSAP080422
      implicit none

c-----input
      integer n,   !  is the order of bessel function , n.ge.0
     &nbes         !  the number of elements of the calculated array
                   !  of bessel function values bes(1:nbes)=J_n(argument(i))
      real*8  bmax !  is the maximal value of the Bessel function argument                  
c-----output
      real*8 bes(*), !the array of Bessel function of order n values
                     !bes(i)=J_n(argument(i)) with i=1,...,nbes
     &besp(*),       !the array of Bessel function of order (n+1) values
                     !bes(i)=J_(n+1)(argument(i)) with i=1,...,nbes
     &besdel         !It is not used
c-----local
      
      integer nz,  ! should be =0, 
                   !If nz=1 the bessel function(dBESJ) will be J=0
     &np,          !=n+1
     &i

      real*8 dbes, !step of bessel function argument 
     &x            !argument of the Bessel function
   
      

c-----it shoud be n.ge.0 for DBESJ
c-----it shoud bmax n.ge.0 for DBESJ
cSAP080422
c      dimension bes(1),besp(1)

      if(n.lt.0) then
        write(*,*)'in sub besvec n<0 '
        write(*,*)'n should be .ge.0'
        stop 'in besvec'
      endif

      if(bmax.lt.0.d0) then
        write(*,*)'in sub besvec bmax<0 '
        write(*,*)'n should bmax .ge.0'
        stop 'in besvec'
      endif
     
      dbes=bmax/(nbes-1)

      do 10  i=1,nbes
         x=(i-1)*dbes
         call DBESJ(x,dble(n),1,bes(i),nz)     
         np=n+1
         call DBESJ(x,dble(np),1,besp(i),nz)
  
         if(n.eq.0) then        
           besp(i)=-besp(i)
         else
c----------n.ne.0 
           if(x.eq.0.d0)then 
              if(n.eq.1)then
                besp(i)=0.5d0
              else
c-------------- n.ge.2
                besp(i)=0.0  d0   
              endif
            else 
c-------------argument x.ne.0
              besp(i)=n*bes(i)/x-besp(i)  !derivative           
            endif
      endif

 10   continue
  
      return
      end

      subroutine bes_calc(z,n,bes_f,bes_prime)
c-----calculate the table for bessel function and its derivative
c     then calculate n_th order bessel function and its derivative
c     with argument z
      implicit none

      include 'param.i'      
      include 'one.i'
c-----input
      integer n !bessel function order
      double precision z !bessel function argument

c-----output
      double precision bes_f,bes_prime

c-----local
      integer nbesa,nbes
c      parameter (nbesa=201)
      parameter (nbesa=10001)


      double precision bmaxon,besdel,z_abs
      
      double precision bes_ar,besp_ar,bmax_ar,dbes,bes_1
      common /bes/ bes_ar(nbesa,n_relt_harm1a:n_relt_harm2a),
     &besp_ar(nbesa,n_relt_harm1a:n_relt_harm2a),
     &bmax_ar(n_relt_harm1a:n_relt_harm2a)

      integer isave,j,n_abs
      save isave
      integer nh,nb,inh,ii,nbesm1,k,nh_abs
      integer nz ! should be =0, 
                 !If nz=1 the bessel function(dBESJ) will be J=0
      data isave/0/

      nbes=nbesa
      besdel=1.d-4
      bmaxon=4.d0
      z_abs=dabs(z)
      n_abs=iabs(n)
    
      if (isave.eq.1) goto10
 
c-----initialization       
      inh=0
      do nh=n_relt_harm1,n_relt_harm2        
         bmax_ar(nh)=bmaxon*(abs(nh)+1)
c--------it will calculate bessel function and its derivative tables
c        for non-negative harmonics number abs(nh) 
         nh_abs=abs(nh)
         call besvec(nh_abs,bes_ar(1,nh),besp_ar(1,nh),nbes,bmax_ar(nh),
     &               besdel)     
      enddo

      isave=1

c      write(*,*)'bes_calc: table for bessel function created'    
c      do nh=n_relt_harm1,n_relt_harm2
c        write(*,*)'nh bmax_ar(nh)',nh,bmax_ar(nh)
c        do j=1,nbes
c          write(*,*)'j,bes_ar(j,nh),besp_ar(j,nh)',
c     &               j,bes_ar(j,nh),besp_ar(j,nh)
c        enddo
c      enddo

 10   continue

c-----calculation of bessel function's interpolation 

      dbes=bmax_ar(n_abs)/(nbes-1)
      nbesm1=nbes-1
      ii=1.d0+z_abs/dbes
      if(ii.le.nbesm1)  go to 110 
cSAP080107
c      bes_f=0.d0
c      bes_prime=0.d0
      call DBESJ(z_abs,dble(n_abs),1,bes_f,nz)
      call DBESJ(z_abs,dble(n_abs+1),1,bes_1,nz)

      if(n_abs.eq.0) then 
        bes_prime=-bes_1
      else
c-------n_abs.ne.0 
c        write(*,*)'in besvec.f before z_abs=0 z_abs',z_abs
        if(z_abs.eq.0.d0)then 
          if(n_abs.eq.1)then
            bes_prime=0.5d0
          else
c-----------n_abs.ge.2
            bes_prime=0.0d0   
          endif
        else 
c---------argument z_abs.ne.0
           bes_prime=-bes_1+n_abs*bes_f/z_abs
        endif
c         write(*,*)'besvec.f in  sub bes_calc bes_prime=',bes_prime
      endif

      goto 120

110   continue


      bes_f=bes_ar(ii,n_abs)+(z_abs-(ii-1)*dbes)*
     &(bes_ar(ii+1,n_abs)-bes_ar(ii,n_abs))/(dbes)
        
      bes_prime=besp_ar(ii,n_abs)+(z_abs-(ii-1)*dbes)*
     &(besp_ar(ii+1,n_abs)-besp_ar(ii,n_abs))/(dbes)

120   continue

      if(z.lt.0d0) then
        bes_f=bes_f*(-1)**n_abs  !negative argument z for J_n(z)
        bes_prime=bes_prime*(-1)**(n_abs+1)  !negative argument z for dJ_n/dz(z)
      endif
    
      if(n.lt.0)then 
         k=(-1)**n_abs ! coefficient for Bessel function J_n with n<0
         bes_f=k*bes_f               ! for negative n
         bes_prime=k*bes_prime       ! for negative n
      endif

      return
      end
