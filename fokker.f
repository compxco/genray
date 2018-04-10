c
      subroutine fokker(rtem0,
     &             rtail,r1t,r2t,ttail,rhot,r1h,r2h,thotpar,thotper,
     &             hotexp,hotmnpar,hotmxpar,hotmnper,hotmxper,
     &             rbeam,r1b,r2b,tbeam,ebeam,thbeam,
     &             jx,iym,lrz,ngen,
     &             jxa,iya,lrza,ngena,
     &             f,x,iy_,y,rovera,vnorm)

c
c     calculates the analytical 3D distribution function in the form of
c     two maxvellian Maxwellian, plus the hot tail, plus the beam.
c     f(v=p/m,pith,rho)=dens(rho)(1-rtatil-rhot-rbeam)*f_relativist_Maxw(T(rho))+
c                                 +rtail*f_tail+rhot*f_thot+rbeam*f_beam) 
c     here dens(rho) is the electron density normalized to 10**19 m**(-3)
c 
c     Input parameters dimensions(iya,jxa,lrza,ngena) are in the file param.i 

      implicit none

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c-----input
      integer 
     .jx, ! the number of used normalized momentum/ mesh points
     .lrz,! the number of used radial mesh points
     .iym,! the number of used pitch-angle mesh points (here the same at each radius)
     .ngen, ! the number of plasma species (here we use only electron specie with 
c            the number of specie k=1) 
     .jxa,iya,lrza,ngena ! the max values for jx,iym,lrz,ngen
      real*8 
     .rtem0, !relation tem0/electron_temperature(rho=0)
c            tem0 is the energy for the momentum normalization (KeV) 
c-----tail parameters
c     f_tail=H(rho,rt1,rt2)*f_rel_Maxw(ttail),
c     H(x,x1,x2) is the box function. H=1 for x1<x<x2 otherwise H=0  
     .r1t,r2t,                !small normalized radii for the tail localization  
     .rtail,                  !the relation the tail density to the total density
     .ttail,                  !tail temperature (KeV)!tail temperature (KeV)
c-----hot parameters
c     f_hot=H(rho,rh1,rh2)*H(epar,hotmnpar,hotmxpar)*H(eper,hotmnper,hotmxper)*
c     *(p_per/mc)**hotexp*exp{-mu(thoppar)(p_par/m_ec)**2-mu(thopper)(p_per/m_ec)**2}.
c     Here mu(T)=m_e*c**2/T  
     .r1h,r2h,                !small normalized radii for the hot localization
     .rhot,                   !the relation of hot hot density to the total density
     .thotpar,thotper,        !parallel and perpendicular hot temperatures (KeV)
     .hotmnpar,hotmxpar,      !the boundaries of the parallel energy box(KeV)
c                              hotmnpar < epar < hotmxpar  
     .hotmnper,hotmxper,      !the boundaries of the perpendicular energy box(KeV)
c                              hotmnper < eper < hotmxper  
     .hotexp,                  !the degree of the perpendicular momentum: (p_per/mc)**hotexp
c-----beam parameters
c     f_beam=H(rho,rb1,rb2)*exp{-0.5*mu(tbeam)*
c              [(p_par-p_beam_par)**2+(p_per-p_beam_per)**2]/(m_e*c)**2}
c     Here
c          (p_beam /m_e*c)**2=ebeam**2/(m_e**2*c**4)-1
c           p_beam_par=p_beam*cos(thbeam)
c           p_beam_per=p_beam*sin(thbeam)

     .r1b,r2b,                !small normalized radii for the beam localization
     .rbeam,                  !the relation of the beam density to the total density
     .ebeam,                  !beam energy (KeV)
     .thbeam,                 !beam pitch angle (0=<degree=<180) 
     .tbeam                   !beam temperature (KeV)

c-----external
      real*8 temperho,densrho

c-----output    
c     electron distribution function. It works for ngena=1 
      real*8 f(iya,jxa,lrza,ngena),
c     the meshs
     .x(jxa), ! momentum/electron_restmass/vnorm mesh: x(j), j=1,jx x(1)=0, x(jx)=1
     .rovera(lrza),!normalized radial mesh: rovera(ir), ir=1,lrz rovera(1)=0,rovera(lrz)=1
     .y(iya,lrza), !pitch angle mesh y(i,ir) ir=1,lrz, i=1,iy_(ir),
c                   y(1,ir)=0,  y(iy_(ir),ir)=2*pi
     .vnorm       ! momentum/electron_restmass for normalization (cm/sec)
      integer 
     .iy_(lrza) !the number of pitch angles at different radii: iy_(ir)=iym, ir=1,lrz

c-----local
      integer k
      integer i,j,ir,itail,ihot,ibeam
      real*8 step,pi,rho,rmu,gamma,fmax,ftail,tem0,
     .dv,vvol,vovc,vovc2,
     .vpar,vpar2,vper,vper2,gampar,gamper,epar,eper,
     .vparb,vperb,gammab,vbovc,fact,
     .dnormt,rmutail,
     .rmubeam,
     .dnormh,rmuhpar,rmuhper,dnormb
      real*8 dnorm(lrza),den,ftemp(iya,jxa)      
      real*8 thets(iya),coss(iya,lrza),sinn(iya,lrza)
c      double precision vs(jxa)

      real*8 mass_e,k_1_kev,clight,m_cs_d_2
ctest       
      real*8 bk2,plogf,penergy

      write(*,*)'begin of fokker.f'

      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !1 egr in 1 KEV
      clight =2.99792458d10     !sm/secC
c-----check the dimensions

      if(jx.gt.jxa) then 
        write(*,*)'it should be jx.le.jxa, but jx>jxa'
        write(*,*)'Please change jx in genray.in or jxa in param.i'
        stop
      endif

      if(lrz.gt.lrza) then 
        write(*,*)'it should be lrz.le.lrza, but lrz>jrza'
        write(*,*)'Please change lrz in genray.in or lrza in param.i'
        stop
      endif

      if(iym.gt.iya) then 
        write(*,*)'it should be iym.le.iya, but iym>iya'
        write(*,*)'Please change iym in genray.in or iya in param.i'
        stop
      endif


c-----creation the grids for the variables:
c     x(jx)=(momentum/restmass/vnorm) 0=<x(j)=<1
c     rovera(lrz)= normalized small radius 
c     y(iya,lrz)= pitch angle. At the given radius rho=rovera(ll)
c     i=1,..iy_(ir)
c
      step=1.d0/dble(jx-1)
      do j=1,jx
         x(j)=step*(j-1) 
      enddo

      step=1.d0/dble(lrz-1)       
      do ir=1,lrz
         rovera(ir)=step*dble(ir-1)
      enddo
      
      do ir=1,lrz
         iy_(ir)=iym
      enddo
     
      pi=4.d0*datan(1.d0) 
            
      do ir=1,lrz
         step=pi/dble(iy_(ir)-1)    
         do i=1,iy_(ir)
            y(i,ir)=step*dble(i-1)
         enddo
      enddo

c      do j=1,jx
c         vs(j)=x(j) 
c      enddo
      
      do ir=1,lrz    
         do i=1,iy_(ir)
            coss(i,ir)=dcos(y(i,ir))
            sinn(i,ir)=dsin(y(i,ir))
         enddo
      enddo

c
c set up electron distribution functions for ECE calcn (Maxwellians)
c
      k=1 ! for the electron distribution

      tem0=rtem0*temperho(0.d0,k)
      vovc2=2.d0*1.6022d-19*1.d3*tem0/9.1095d-31/2.9979d8**2
      vovc=dsqrt(vovc2) !v_t/c
   
      vnorm=vovc*2.9979d10    !cm/sec
      write(*,*)'in fokker.f rtem0= ',rtem0
      write(*,*)'in fokker.f temperho(0.d0,k)= ',temperho(0.d0,k)
      write(*,*)'fokker vovc',vovc
      pi=4.d0*datan(1.d0)

      do 1 ir=1,lrz

        rho=rovera(ir)
        dnorm(ir)=0.0d0
        rmu = 2.d0/(0.00391392d0*temperho(rho,k))
        rmu =mass_e*clight**2/(k_1_kev*temperho(rho,k))

        write(*,*)'fokker ir,rho,temperho(rho,k),rmu',
     .  ir,rho,temperho(rho,k),rmu

        call besk2as(rmu,bk2)

        den=densrho(rho,k)   !electron density 10**19 m**(-3)

cSmtest
c        den=1.d0

c        write(*,*)'fokker ir,den',ir,den
        if (rho.ge.r1t.and.rho.le.r2t.and.rtail.ne.0.d0) then
           itail=1
           dnormt=0.d0
           rmutail=2.d0/(0.00391392d0*ttail)
        else
           itail=0
           dnormt=1.d0
        endif
c        write(*,*)'fokker r1t,r2t,ttail,itail',r1t,r2t,ttail,itail

        if (rho.ge.r1h.and.rho.le.r2h.and.rhot.ne.0.d0) then
           ihot=1
        else
           ihot=0
        endif

        if (rho.ge.r1b.and.rho.le.r2b.and.rbeam.ne.0.d0) then
           ibeam=1
        else
           ibeam=0
        endif
        
        do i=1,iy_(ir)
           thets(i)=y(i,ir)
        enddo
        
        do j=1,jx !x
           gamma=dsqrt(1.d0+x(j)*x(j)*vovc2)
          
           fmax=dexp(-rmu*(gamma-1.d0))
           fmax=fmax*rmu/(4.d0*pi*bk2)
c           write(*,*)'fokker0 j,x(j),gamma,fmax',j,x(j),gamma,fmax

           if(itail.ne.0)  ftail=dexp(-rmutail*(gamma-1.d0))
c           write(*,*)'itail,ftail',itail,ftail
           dv=0.d0
           if(j.eq.1) dv=(0.5d0*(x(2)+x(1)))**3/3.d0
           if(j.gt.1.and.j.lt.jx) dv=0.5d0*(x(j+1)-x(j-1))
           vvol=4.d0*pi*x(j)*x(j)*dv
           dnorm(ir)=dnorm(ir)+vvol*fmax

           do i=1,iy_(ir)
              f(i,j,ir,k)=fmax
c              write(*,*)'fokker j,f(i,j,ir,k)=fmax',j,f(i,j,ir,k)
           enddo

c----------initilization
           do i=1,iy_(ir)
              ftemp(i,j)=0.d0
           enddo


           if (itail.ne.0)  then
              dnormt=dnormt+vvol*ftail
              do i=1,iy_(ir)
                 ftemp(i,j)=ftail
c                 write(*,*)'fokker tail i,ftepm(i,j)',i,ftemp(i,j)
              enddo
           endif
c           write(*,*)'fokker dnormt',dnormt
        enddo

c        write(*,*)'fokker ir,rmu,bk2,dnorm(ir)',ir,rmu,bk2,dnorm(ir)
        do j=1,jx
           do i=1,iy_(ir)
c           write(*,*)'fokker max i,j,ir,den,dnorm(ir),dnormt',
c     &     i,j,ir,den,dnorm(ir),dnormt
c           write(*,*)'fokkker 0 f(i,j,ir,1)',
c     &     f(i,j,ir,1)
c           write(*,*)'fokker itail,rtail,ihot,rhot,ibeam,rbeam',
c     &     itail,rtail,ihot,rhot,ibeam,rbeam
c           write(*,*)'fokker ftemp(i,j)',ftemp(i,j)
            f(i,j,ir,1)=den*
     &            ((1.d0-itail*rtail-ihot*rhot-ibeam*rbeam)
     &            *f(i,j,ir,1)/dnorm(ir)
     &            +itail*rtail*ftemp(i,j)/dnormt)
          

c           write(*,*)'fokkker f(i,j,ir,1)',
c     &     f(i,j,ir,1)

           enddo

        enddo
ctest   
c        write(*,*)'ir dnormt,dnorm(ir)',ir,dnormt,dnorm(ir)
c        dnormt=0.d0
c        do j=1,jx
c           dv=0.d0
c           if(j.eq.1) dv=(0.5d0*(x(2)+x(1)))**3/3.d0
c           if(j.gt.1.and.j.lt.jx) dv=0.5d0*(x(j+1)-x(j-1))
c           vvol=2.d0*pi*x(j)*x(j)*dv
c           do i=2,iy_(ir)
c              dnormt=dnormt+vvol*0.5d0*(f(i-1,j,ir,k)+f(i,j,ir,k))*
c     .        dsin(0.5d0*(thets(i-1)+thets(i)))*(thets(i)-thets(i-1))
c           enddo
c        enddo
c        write(*,*)'dnormt=',dnormt
cend test
      
c
c       Add in hot component
        if (ihot.eq.1)  then
           dnormh=0.d0
           rmuhpar=2.d0/(0.00391392d0*thotpar)
           rmuhper=2.d0/(0.00391392d0*thotper)

           do j=1,jx
              dv=0.d0
              if(j.eq.1) dv=(0.5d0*(x(2)+x(1)))**3/3.d0
              if(j.gt.1.and.j.lt.jx) dv=0.5d0*(x(j+1)-x(j-1))
c              write(*,*)'fokker dhot j,dv',j,dv
              do i=1,iy_(ir)
                 vpar2=0.5d0*(coss(i,ir)*x(j))**2*vovc2
                 vper2=0.5d0*(sinn(i,ir)*x(j))**2*vovc2
                 if (hotexp.eq.0.d0) then
                    fact=1.d0
                 else
                    fact=vper2**(hotexp/2.d0)
                 endif
                 ftemp(i,j)=fact*dexp(-rmuhpar*vpar2-rmuhper*vper2)
                 gampar=sqrt(1.d0+2.d0*vpar2)
                 gamper=sqrt(1.d0+2.d0*vper2)
                 epar=(gampar-1.d0)*511.d0
                 eper=(gamper-1.d0)*511.d0
                 if (epar.lt.hotmnpar.or.epar.gt.hotmxpar) 
     +              ftemp(i,j)=0.d0
                 if (eper.lt.hotmnper.or.eper.gt.hotmxper)
     +              ftemp(i,j)=0.d0
              enddo

              do i=2,iy_(ir)
                 dnormh=dnormh+2.d0*pi*x(j)**2*dv*0.5d0*(ftemp(i-1,j)+
     1           ftemp(i,j))
     &         *dsin(0.5d0*(thets(i-1)+thets(i)))*(thets(i)-thets(i-1))
              enddo

           enddo

           if (dnormh.eq.0.d0) then
              write(*,*) "Prblm with total hot tail density. Check nml."
              STOP
           endif

           do j=1,jx
              do i=1,iy_(ir)
                f(i,j,ir,k)=f(i,j,ir,k)+rhot*den*ftemp(i,j)/dnormh
c                write(*,*)'fokker hot j,i,ir,k,f',j,i,ir,k,f(i,j,ir,k)
              enddo
           enddo

        endif
c
c
c       Add in beam component
        if (ibeam.eq.1)  then

           dnormb=0.d0
           rmubeam=2.d0/(0.00391392d0*tbeam)
           gammab=ebeam/511.d0+1.d0
           vbovc=dsqrt(gammab**2-1.d0)
           vparb=vbovc*dcos(thbeam*pi/180.d0)
           vperb=vbovc*dsin(thbeam*pi/180.d0)

           do j=1,jx
              dv=0.d0
              if(j.eq.1) dv=(0.5d0*(x(2)+x(1)))**3/3.d0
              if(j.gt.1.and.j.lt.jx) dv=0.5d0*(x(j+1)-x(j-1))

              do i=1,iy_(ir)
                 vpar=coss(i,ir)*x(j)*vovc
                 vper=sinn(i,ir)*x(j)*vovc
                 ftemp(i,j)=dexp(-0.5d0*rmubeam*((vpar-vparb)**2+
     1           (vper-vperb)**2))
              enddo
 
              do i=2,iy_(ir)
                 dnormb=dnormb+2.d0*pi*x(j)**2*dv*0.5d0*(ftemp(i-1,j)+
     1           ftemp(i,j))
     &         *dsin(0.5d0*(thets(i-1)+thets(i)))*(thets(i)-thets(i-1))
             enddo
           enddo

           if (dnormb.eq.0.d0) then
              write(*,*) "Prblm with total beam, density. Check nml."
              STOP
           endif
               
           do j=1,jx
              do i=1,iy_(ir)
                f(i,j,ir,k)=f(i,j,ir,k)+rbeam*den*ftemp(i,j)/dnormb
              enddo 
           enddo

        endif

1     continue

ctest
c      write(*,*)'end of fokker.f'
c      write(*,*)'jx,lrz',jx,lrz
c      write(*,*)'x',x
c      write(*,*)'iy_',iy_
c      do ir=1,lrz
c        do i=1,iy_(ir)
c          write(*,*)'i,ir,y(i,ir)',i,ir,y(i,ir)
c        enddo
c      enddo
         
      m_cs_d_2=0.5d0*mass_e*clight**2/k_1_kev  !    [KeV] m*c**2/2

      open(17,file='fokker.bin',form='unformatted')
      do ir=1,lrz
        dnorm(ir)=0.0d0
        do j=1,jx !
           dv=0.d0
           if(j.eq.1) dv=(0.5d0*(x(2)+x(1)))**3/3.d0
           if(j.gt.1.and.j.lt.jx) dv=0.5d0*(x(j+1)-x(j-1))
           vvol=2.d0*pi*x(j)*x(j)*dv
           
           do i=2,iy_(ir)
              dnorm(ir)=dnorm(ir)+vvol*0.5d0*(f(i-1,j,ir,k)+f(i,j,ir,k))*
     .        *dsin(0.5d0*(thets(i-1)+thets(i)))*(thets(i)-thets(i-1))
c              write(*,*)'fokker ir,j,i,f(i,j,ir,1)',ir,j,i,f(i,j,ir,1)
           enddo
        enddo
c        write(*,*)'ir dnorm(ir)',ir,dnorm(ir)


        do j=1,jx

          if (f(1,j,ir,1).gt.1.d-300) then
             plogf=dlog(f(1,j,ir,1))
          else
             plogf=-728
          endif

c         penergy=(dsqrt(1.d0+(x(j)*vovc)**2)-1.d0)*(2.d0/0.00391392d0)
c--------------------------------------------------------------------
c        kinetic energy [KeV]
c--------------------------------------------------------------------
         penergy=(dsqrt(1.d0+(x(j)*vovc)**2)-1.d0)*(2.d0*m_cs_d_2) !kev

         if(myrank.eq.0) WRITE(17) REAL(penergy),REAL(plogf)
c          write(*,*) 'ir,j,x(j)**2,dlog(f(1,j,ir,1)),f(1,j,ir,1)',
c     .    ir,j,penergy,plogf,f(1,j,ir,1)
        enddo
        if(myrank.eq.0) WRITE(17)

      enddo

cend test
      return   
      end

      subroutine fokker_3_temperature(rtem0,
     & rtemp1, rtemp2, rtemp3,
     & rvtovte1,rvtovte2,
     & jx,iym,lrz,jxa,iya,lrza,ngena,
     & f,x,iy_,y,rovera,vnorm)

c
c     calculates the spline analytical 3D electron distribution
c     function in the form of
c     Maxwellian with 3 temperature bins and put this function to 3D grid 
c-----Three temperature case  (i_diskf=4)
c     rvtovte1,rvtovte2 = ratio of momentum-per-mass (electrons) to on-axis
c                       thermal velocity vte0= sqrt(Te/me), defining the
c                       three velocity ranges for the temperatures.
c                       defaults=1.e6,1.e6 [i.e., effectively infinity]
c     rtemp1, rtemp2, rtemp3 = ratios of temperatures in each of the
c                       three velocity (energy) bins to the radially
c                       local temperature.
c     In summary:  The three momentum-per-mass bins are [0.,rvtovte1*vte0],
c                [rvtovte1*vte0,rvtovte2*vte0], and [rvtovte2*vte0,infinity].
c                These bins are constant as a function of radius.
c                The temperatures in each bin are given by rtemp[1-3]
c                and vary as a function of radius as the bulk temperature.
c
c     f(rho,v=p/(m_e*c))=dens(rho)(C/4pi)exp[-mu_3_bin(rho,v)(gamma*-1)]
c             f(iya,jxa,lrza,ngena)=f                       
c     Here  
c
c     dens(rho) is the electron density normalized to 10**19 m**(-3)
c
c     mu_3_bin(rho,v)=m_e*clight**2/T_3_bin(rho,v)
c
c     T_3_bin(rho,v)=T_e(rho){rtemp1*H(0,v,rvt0te1*V_Te(rho=0))+
c                             rtemp2*H(rvt0te1*V_Te(rho=0)),v,rvt0te2*V_Te(rho=0))+
c                             rtemp3*H(rvt0te2*V_Te(rho=0)),v,infinity}
c      
c     V_Te(rho=0)=sqrt(T_e(rho=0)/m_e) on-axis thermal velocity
c
c     H(x1,x,x2)=1 at x1 =< x <x2
c               =0    otherwise   
c
c     Input parameters dimensions(iya,jxa,lrza,ngena) are set in the file param.i 

      implicit none

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c-----input
      integer 
     &jx, ! the number of used normalized momentum/ mesh points
     &lrz,! the number of used radial mesh points
     &iym,! the number of used pitch-angle mesh points (here the same at each radius)
     
     &jxa,iya,lrza,ngena ! the max values for jx,iym,lrz,ngen
      double precision rtem0, !relation tem0/electron_temperature(rho=0)
c                              tem0 is the energy for the momentum normalization (KeV) 
     &rtemp1, rtemp2, rtemp3, !relation bin temperature/electron_temperature(rho)
c-----momentum bins boundaries
     &rvtovte1,rvtovte2 != ratio of momentum-per-mass (electrons) to on-axis
                        !thermal velocity vte0= sqrt(Te/me), defining the
                        !three velocity ranges for the temperatures.
                        !defaults=1.e6,1.e6 [i.e., effectively infinity]
c-----external
      double precision temperho,densrho

c-----output    
c     electron distribution function. It works for ngen=1 
      double precision f(iya,jxa,lrza,ngena),
c     the mesh
     &x(jxa), ! momentum/electron_restmass/vnorm mesh: x(j), j=1,jx x(1)=0, x(jx)=1
     &rovera(lrza),!normalized radial mesh: rovera(ir), ir=1,lrz rovera(1)=0,rovera(lrz)=1
     &y(iya,lrza), !pitch angle mesh y(i,ir) ir=1,lrz, i=1,iy_(ir),
c                   y(1,ir)=0,  y(iy_(ir),ir)=2*pi
     &vnorm       ! momentum/electron_restmass for normalization (cm/sec)
      integer 
     &iy_(lrza) !the number of pitch angles at different radii: iy_(ir)=iym, ir=1,lrz

c-----local
      integer i,j,ir,jb1,jb2
      double precision step,pi,rho,rmu,gamma,f_maxwell,tem0,bin_factor,
     &dv,vvol,vovc,vovc2,dnormt
      double precision den 
      double precision thets(iya),coss(iya,lrza),sinn(iya,lrza)
c      double precision vs(jxa)

      double precision mass_e,k_1_kev,clight,m_cs_d_2,
     & const_bin_1,const_bin_2,const_bin_3,bin_const,gamma1,gamma2
ctest        
      double precision plogf,penergy
      double precision t_test

      write(*,*)'begin of fokker_3_temperature'

      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !1 egr in 1 KEV
      clight =2.99792458d10     !sm/secC
c-----check the dimensions

      if(jx.gt.jxa) then 
        write(*,*)'it should be jx.le.jxa, but jx>jxa'
        write(*,*)'Please change jx in genray.in or jxa in param.i'
        stop
      endif

      if(lrz.gt.lrza) then 
        write(*,*)'it should be lrz.le.lrza, but lrz>jrza'
        write(*,*)'Please change lrz in genray.in or lrza in param.i'
        stop
      endif

      if(iym.gt.iya) then 
        write(*,*)'it should be iym.le.iya, but iym>iya'
        write(*,*)'Please change iym in genray.in or iya in param.i'
        stop
      endif


c-----creation the grids for the variables:
c     x(jx)=(momentum/restmass/vnorm) 0=<x(j)=<1
c     rovera(lrz)= normalized small radius 
c     y(iya,lrz)= pitch angle. At the given radius rho=rovera(ll)
c     i=1,..iy_(ir)
c
      step=1.d0/dble(jx-1)
      do j=1,jx
         x(j)=step*(j-1) 
      enddo

      step=1.d0/dble(lrz-1)       
      do ir=1,lrz
         rovera(ir)=step*dble(ir-1)
      enddo
      
      do ir=1,lrz
         iy_(ir)=iym
      enddo
     
      pi=4.d0*datan(1.d0) 
            
      do ir=1,lrz
         step=pi/dble(iy_(ir)-1)    
         do i=1,iy_(ir)
            y(i,ir)=step*dble(i-1)
         enddo
      enddo

c      do j=1,jx
c         vs(j)=x(j) 
c      enddo
      
      
c
c-----set up electron distribution functions for ECE calcn (32 Maxwellians)
c

      tem0=rtem0*temperho(0.d0,1)
      vovc2=2.d0*1.6022d-19*1.d3*tem0/9.1095d-31/2.9979d8**2
      vovc=dsqrt(vovc2) !v_t/c
   
      vnorm=vovc*2.9979d10    !cm/sec
      write(*,*)'in fokker.f temperho(0.d0,1)= ',temperho(0.d0,1)
      write(*,*)'fokker vovc',vovc
ctest
      t_test=(vnorm**2*mass_e/2.d0)/k_1_kev  !Kev
      write(*,*)'t_test [kev]=',t_test 

      den=densrho(rho,1)   !electron density 10**19 m**(-3)
               
      do ir=1,lrz    
         do i=1,iy_(ir)
            coss(i,ir)=dcos(y(i,ir))
            sinn(i,ir)=dsin(y(i,ir))
         enddo
      enddo

c-------------------------------------------------------------------
      const_bin_1=0.d0
      jb1=jx
      jb2=jx
      do j=1,jx-1 
         if ((x(j).lt.rvtovte1/dsqrt(2.d0*rtem0)).and.
     &       (x(j+1).ge.rvtovte1/dsqrt(2.d0*rtem0))) jb1=j+1                   
         if ((x(j).lt.rvtovte2/dsqrt(2.d0*rtem0)).and.
     &       (x(j+1).ge.rvtovte2/dsqrt(2.d0*rtem0))) jb2=j+1
      enddo             
      write(*,*)'fokker.f jb1,jb2',jb1,jb2

      gamma1=dsqrt(1.d0+x(jb1)*x(jb1)*vovc2)
      gamma2=dsqrt(1.d0+x(jb2)*x(jb2)*vovc2)
      clight =2.99792458d10     !speed of light     (cm/sec)
      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)  
      write(*,*)'1: gamma1,gamma2',gamma1,gamma2      
c      gamma1=dsqrt(1.d0+rvtovte1*vovc2/(2.d0*rtem0))
c      gamma2=dsqrt(1.d0+rvtovte2*vovc2/(2.d0*rtem0))
c      write(*,*)'2: gamma1,gamma2',gamma1,gamma2 
c-------------------------------------------------------------------

      do ir=1,lrz

        rho=rovera(ir)
        den=densrho(rho,1)  

        rmu = 2.d0/(0.00391392d0*temperho(rho,1))
        rmu =mass_e*clight**2/(k_1_kev*temperho(rho,1))

        write(*,*)'ir,rho,temperho(rho,1),rmu',
     &             ir,rho,temperho(rho,1),rmu
                            
        do i=1,iy_(ir)
           thets(i)=y(i,ir)
        enddo
        write(*,*)'rvtovte1,rvtovte2',rvtovte1,rvtovte2
c---------------------------------------------------------------------
c
c       x = P/(m_e*vnorm)
c       vnorm = sqrt[2*T(rho=0)*rtem0/m_e]
c       rvtovt1 = (P/m_e)/sqrt[T(rho=0)/m_e]
c
c       At the velosity range boundary:
c      
c       P_boundary = x_boundary*m_e*vnorm = rvtovt1*m_e*sqrt[T(rho=0)/m_e]
c
c       x_boundary = rvtovt1 *sqrt[T(0)/m_e]/sqrt[2T(0)*rtem0/m_e] =
c                  = rvtovt1/sqrt[2*rtem0]
c      
c                        rvtovt1/sqrt(2rtem0)   rvtovt2/sqrt(2rtem0)
c                             |                      |
c      0--------------------------------------------------------------------> x
c      |    T=T(rho)*rtemp1   |    T=T(rho)*rtemp2   |     T=T(rho)*rtemp3
c
c-----------------------------------------------------------------------------
                
        const_bin_2=(gamma1-1.d0)*rmu*(1.d0/rtemp2-1.d0/rtemp1)     
        const_bin_3= const_bin_2+
     &             (gamma2-1.d0)*rmu*(1.d0/rtemp3-1.d0/rtemp2)
       
        do j=1,jx !x
cSm070109
c           if(x(j).lt.rvtovte1*vovc) then

            if(x(j).lt.rvtovte1/dsqrt(2.d0*rtem0)) then
c------------first velocity bin
             bin_factor=1.d0/rtemp1
             bin_const=const_bin_1
             write(*,*)'jb1,j',jb1,j
           else
             if (x(j).lt.rvtovte2/dsqrt(2.d0*rtem0)) then
c---------------second velosity bin
                bin_factor=1.d0/rtemp2
                bin_const=const_bin_2
                write(*,*)'jb2,j,bin_const',jb2,j
c                gamma1=dsqrt(1.d0+x(jb1)*x(jb1)*vovc2)
c                const_bin_2=(gamma1-1.d0)*rmu*(1.d0/rtemp2-1.d0/rtemp1) 
c                bin_const=const_bin_2
             else
c---------------third velocity bin
                bin_factor=1.d0/rtemp3
                bin_const=const_bin_3
                write(*,*)'jb2,j',jb2,j
             endif
           endif

c           if (j.lt.jb1) then
c              bin_factor=1.d0/rtemp1
c              bin_const=const_bin_1
c           else
c              if (j.lt.jb2) then
c                 bin_factor=1.d0/rtemp2
c                 bin_const=const_bin_2
c              else
c                bin_factor=1.d0/rtemp3
c                bin_const=const_bin_3
c              endif
c           endif
 
      write(*,*)'j,x(j),rvtovte1/dsqrt(2.d0*rtem0)',
     &           j,x(j),rvtovte1/dsqrt(2.d0*rtem0)
      write(*,*)'j,x(j),rvtovte2/dsqrt(2.d0*rtem0)',
     &           j,x(j),rvtovte2/dsqrt(2.d0*rtem0)

           write(*,*)'rtemp1,rtemp2,bin_factor',rtemp1,rtemp2,bin_factor


           gamma=dsqrt(1.d0+x(j)*x(j)*vovc2)
          
           f_maxwell=dexp(-rmu*bin_factor*(gamma-1.d0)+bin_const)
                       
           do i=1,iy_(ir)
             f(i,j,ir,1)=den*f_maxwell
           enddo !i
        enddo !j

      
        dnormt=0.d0
        do j=1,jx
           dv=0.d0
           if(j.eq.1) dv=(0.5d0*(x(2)+x(1)))**3/3.d0
           if(j.gt.1.and.j.lt.jx) dv=0.5d0*(x(j+1)-x(j-1))
           vvol=2.d0*pi*x(j)*x(j)*dv
           do i=2,iy_(ir)
              dnormt=dnormt+vvol*0.5d0*(f(i-1,j,ir,1)+f(i,j,ir,1))*
     &        dsin(0.5d0*(thets(i-1)+thets(i)))*(thets(i)-thets(i-1))
           enddo
        enddo
        write(*,*)'ir,dnormt',ir,dnormt

        do j=1,jx
           do i=1,iy_(ir)
             f(i,j,ir,1)=den*f(i,j,ir,1)/dnormt
           enddo !i
        enddo !j

      enddo !ir

      open(17,file='fokker.bin',form='unformatted')


      m_cs_d_2=0.5d0*mass_e*clight**2/k_1_kev  !    [KeV] m*c**2/2
ctest  
      do ir=1,lrz
        dnormt=0.0d0
        do j=1,jx !x
           dv=0.d0
           if(j.eq.1) dv=(0.5d0*(x(2)+x(1)))**3/3.d0
           if(j.gt.1.and.j.lt.jx) dv=0.5d0*(x(j+1)-x(j-1))
           vvol=2.d0*pi*x(j)*x(j)*dv
           
           do i=2,iy_(ir)
              dnormt=dnormt+vvol*0.5d0*(f(i-1,j,ir,1)+f(i,j,ir,1))*
     &        *dsin(0.5d0*(thets(i-1)+thets(i)))*(thets(i)-thets(i-1))
           enddo
        enddo

c_test
        rho=rovera(ir)      
        den=densrho(rho,1)  
        write(*,*)'ir dnormt,den',ir,dnormt,den
c_end_test

        do j=1,jx

          if (f(1,j,ir,1).gt.1.d-300) then
             plogf=dlog(f(1,j,ir,1))
          else
             plogf=-728
          endif

c         penergy=(dsqrt(1.d0+(x(j)*vovc)**2)-1.d0)*(2.d0/0.00391392d0)
c--------------------------------------------------------------------
c        kinetic energy [KeV]
c--------------------------------------------------------------------
         penergy=(dsqrt(1.d0+(x(j)*vovc)**2)-1.d0)*(2.d0*m_cs_d_2) !kev
          if(myrank.eq.0) WRITE(17) REAL(penergy),REAL(plogf)
c          write(*,*) 'ir,j,x(j)**2,dlog(f(1,j,ir,1)),f(1,j,ir,1)',
c     &    ir,j,penergy,plogf,f(1,j,ir,1)
        enddo
        if(myrank.eq.0) WRITE(17)

      enddo !ir

      return   
      end

      subroutine fokker_3_temperature_analytical(rtem0,
     & rtemp1, rtemp2, rtemp3,
     & rvtovte1,rvtovte2,
     & jx,jxa,x,vnorm,ndens,ndensa,rhom,
     & p_in,rho_in,
     & f_analytical,df_dp)

c
c     Calculates:
c     1) The pure analytical 3D electron distribution function
c         f_analytical
c        (for i_diskf=5 case)
c        in the form of Maxwellian with 3 temperature bins
c        at the given small radius rho_in an the normalized
c          momentum p=momentum/mc
c                                    x=P/(m_e*vnorm) 
c
c     2) The derivative from this distribution by x: df_dx
c
c     3) vnorm =sqrt(2T(rho=0)*rtem0/me) the velocity for normalization
c
c-----Three temperature case  (i_diskf=4)
c     rvtovte1,rvtovte2 = ratio of momentum-per-mass (electrons) to on-axis
c                       thermal velocity vte0= sqrt(Te/me), defining the
c                       three velocity ranges for the temperatures.
c                       defaults=1.e6,1.e6 [i.e., effectively infinity]
c     rtemp1, rtemp2, rtemp3 = ratios of temperatures in each of the
c                       three velocity (energy) bins to the radially
c                       local temperature.
c     In summary:  The three momentum-per-mass bins are [0.,rvtovte1*vte0],
c                [rvtovte1*vte0,rvtovte2*vte0], and [rvtovte2*vte0,infinity].
c                These bins are constant as a function of radius.
c                The temperatures in each bin are given by rtemp[1-3]
c                and vary as a function of radius as the bulk temperature.
c
c     f(rho,v=p/(m_e*c))=dens(rho)(C/4pi)exp[-mu_3_bin(rho,v)(gamma*-1)]
c             f(iya,jxa,lrza,ngena)=f       
c             It will be here dens(rho)=1 to get normalized disrtibution                 
c     Here  
c
c     dens(rho) is the electron density normalized to 10**19 m**(-3)
c
c     mu_3_bin(rho,v)=m_e*clight**2/T_3_bin(rho,v)
c
c     T_3_bin(rho,v)=T_e(rho){rtemp1*H(0,v,rvt0te1*V_Te(rho=0))+
c                             rtemp2*H(rvt0te1*V_Te(rho=0)),v,rvt0te2*V_Te(rho=0))+
c                             rtemp3*H(rvt0te2*V_Te(rho=0)),v,infinity}
c      
c     V_Te(rho=0)=sqrt(T_e(rho=0)/m_e) on-axis thermal velocity
c
c     H(x1,x,x2)=1 at x1 =< x <x2
c               =0    otherwise   
c
c     Input parameters dimensions (jxa,ndensa) are set in the file param.i 

      implicit none

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c-----input
      integer 
     &jx,      ! the number of used normalized momentum/ mesh points
     &jxa,     ! maximal value of jx
     &ndens,   ! number of radial points used for plotting
     &ndensa   ! maximal value of ndens
      double precision rhom(*),! rhom(ndens) the small radius array 
                               ! used for plotting
     *rtem0,                   ! relation tem0/electron_temperature(rho=0)
                               ! tem0 is the energy for the momentum
                               ! normalization (KeV) 
     &rtemp1, rtemp2, rtemp3,  ! relation bin temperature/electron_temperature(rho)
c-----momentum bins boundaries
     &rvtovte1,rvtovte2,!= ratio of momentum-per-mass (electrons) to on-axis
                        !thermal velocity vte0= sqrt(Te/me), defining the
                        !three velocity ranges for the temperatures.
                        !defaults=1.e6,1.e6 [i.e., effectively infinity]
     &p_in,             ! P/(m_e*clight) 
     &rho_in            ! normalized small radius 
c-----external
      double precision temperho,densrho

c-----output    
c     electron distribution function. It works for ngen=1 
      double precision 
c     the mesh
     &x(jxa),      ! momentum/electron_restmass/vnorm mesh:
                   !  x(j), j=1,jx x(1)=0, x(jx)=1
     &vnorm,       ! momentum/electron_restmass for normalization (cm/sec)
     &f_analytical,
     &df_dp
c---------------------------------------------------------------------
c
c       x = P/(m_e*vnorm) = p*clight/vnorm
c       vnorm = sqrt[2*T(rho=0)*rtem0/m_e]
c       rvtovt1 = (P/m_e)/sqrt[T(rho=0)/m_e]
c
c       At the velosity range boundary:
c      
c       P_boundary = x_boundary*m_e*vnorm = rvtovt1*m_e*sqrt[T(rho=0)/m_e]
c
c       x_boundary = rvtovt1 *sqrt[T(0)/m_e]/sqrt[2T(0)*rtem0/m_e] =
c                  = rvtovt1/sqrt[2*rtem0]
c      
c                        rvtovt1/sqrt(2rtem0)   rvtovt2/sqrt(2rtem0)
c                             |                      |
c      0--------------------------------------------------------------------> x
c      |    T=T(rho)*rtemp1   |    T=T(rho)*rtemp2   |     T=T(rho)*rtemp3
c-----------------------------------------------------------------------------

c-----local
      integer i,j,ir,jb1,jb2
      double precision step,pi,rho,f_3T_Maxwell,tem0,bin_factor,
     &dv,vvol,vovc,vovc2,dnormt,x_in,df_dx
      double precision den
c      double precision vs(jxa)

      double precision mass_e,k_1_kev,clight,m_cs_d_2

      integer ir_maxwell 
      parameter (ir_maxwell=100)
      double precision ratio_dens(ir_maxwell),rho_maxwell(ir_maxwell),
     & ratio_dens_spl


c-------------------------------------------------------------------
c     for ratio_dens spline 
c------------------------------------------------------------------
      integer i1p(2),
     &itabl(3) 
      double precision d2_ratio_dens_drho(ir_maxwell),
     &work(3*ir_maxwell+1),
     &tabl(3)
       

      integer  i_first
      data i_first/0/
      save i_first,
     &ratio_dens,rho_maxwell,d2_ratio_dens_drho

ctest 
      double precision theta_test
c-----external
c      double precision temperho

ctest        
      double precision plogf,penergy

c      write(*,*)'begin of fokker_3_temperature_analytical'

      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !1 egr in 1 KEV
      clight =2.99792458d10     !speed of light [sm/sec]
     
c-----check the dimensions

      if(jx.gt.jxa) then 
        write(*,*)'it should be jx.le.jxa, but jx>jxa'
        write(*,*)'Please change jx in genray.in or jxa in param.i'
        stop
      endif

      if(ndens.gt.ndensa) then 
        write(*,*)'it should be ndens.le.ndensa, but ndens>ndensa'
        write(*,*)'Please change ndens in genray.in or ensa in param.i'
        stop
      endif

c-----creation the grids for the variables:
c     x(jx)=(momentum/restmass/vnorm) 0=<x(j)=<1
c     r_ar(lrz)= normalized small radius 
c
      step=1.d0/dble(jx-1)
      do j=1,jx
         x(j)=step*(j-1) 
      enddo    
 
      step=1.d0/dble(ir_maxwell-1)       
      do ir=1,ir_maxwell
         rho_maxwell(ir)=step*dble(ir-1)
      enddo
     
      pi=4.d0*datan(1.d0) 
            
c     do j=1,jx
c        vs(j)=x(j) 
c     enddo

c-----set up electron distribution functions for ECE calcn (32 Maxwellians)

      tem0=rtem0*temperho(0.d0,1)
c      vovc2=2.d0*1.6022d-19*1.d3*tem0/9.1095d-31/2.9979d8**2
      vovc2=(2.d0*tem0*k_1_kev/mass_e)/clight**2
      vovc=dsqrt(vovc2) !v_t/c
   
      vnorm=vovc*clight         !cm/sec
c      write(*,*)'in fokker.f temperho(0.d0,1)= ',temperho(0.d0,1)
c      write(*,*)'fokker vovc,vnorm',vovc,vnorm
      x_in=p_in*clight/vnorm
c      write(*,*)'fokker p_in,x_in,i_first',p_in,x_in,i_first
c      write(*,*)'clight,vnorm',clight,vnorm
c      write(*,*)'clight/vnorm,vnorm/clight',clight/vnorm,vnorm/clight
c-------------------------------------------------------------------
      if (i_first.eq.0) then
c-------------------------------------------------------------------
c        initialization
c------------------------------------------------------------------
         i_first=1
            
c-------------------------------------------------------------------
         do ir=1,ir_maxwell
            rho=rho_maxwell(ir)
c            den=densrho(rho,1)  !electron density 10**19 m**(-3)
            den=1.d0             ! normalized electron density 10**19 m**(-3)
            dnormt=0.d0
            do j=1,jx
               call Three_Temp_Maxwell_1D(rtem0,rtemp1, rtemp2, rtemp3,
     &         rvtovte1,rvtovte2,x(j), rho,f_3T_Maxwell,df_dx)       
               dv=0.d0
               if(j.eq.1) dv=(0.5d0*(x(2)+x(1)))**3/3.d0
               if(j.gt.1.and.j.lt.jx) dv=0.5d0*(x(j+1)-x(j-1))
               vvol=4.d0*pi*x(j)*x(j)*dv
               vvol=vvol*(vnorm/clight)**3
               dnormt=dnormt+vvol*f_3T_Maxwell           
            enddo !j

            write(*,*)'ir,dnormt',ir,dnormt

            ratio_dens(ir)=den/dnormt

         enddo !ir
c--------------------------------------------------------
         i1p(1)=4
         i1p(2)=4
c--------------------------------------------------------
         call coeff1(ir_maxwell,rho_maxwell,ratio_dens,
     &               d2_ratio_dens_drho,i1p,1,work) 

c-----------------------------------------------------------
ctest  
        itabl(1)=1
        itabl(2)=0
        itabl(3)=0
 
        do ir=1,ir_maxwell
            rho=rho_maxwell(ir)
            den=densrho(rho,1)
            call terp1(ir_maxwell,rho_maxwell,ratio_dens,
     &               d2_ratio_dens_drho,rho,1,tabl,itabl)
            ratio_dens_spl=tabl(1)

            write(*,*)'ir,ratio_dens(ir),ratio_dens_spl',
     &                 ir,ratio_dens(ir),ratio_dens_spl
        enddo ! ir
c end test
c------------------------------------------------------------
c       Write the distribution to fokker.bin for xdraw plotting
c------------------------------------------------------------
        open(17,file='fokker.bin',form='unformatted')
        m_cs_d_2=0.5d0*mass_e*clight**2/k_1_kev  !    [KeV] m*c**2/2
  
        itabl(1)=1
        itabl(2)=0
        itabl(3)=0

        do ir=1,ndens
           rho_in=rhom(ir)        
           call terp1(ir_maxwell,rho_maxwell,ratio_dens,
     &               d2_ratio_dens_drho,rho_in,1,tabl,itabl)

           ratio_dens_spl=tabl(1)

           do j=1,jx
             call Three_Temp_Maxwell_1D(rtem0,rtemp1, rtemp2, rtemp3,
     &       rvtovte1,rvtovte2,x(j),rho_in,f_3T_Maxwell,df_dx)

             f_analytical=f_3T_Maxwell*ratio_dens_spl
             if (f_analytical.gt.1.d-300) then
                plogf=dlog(f_analytical)
             else
                plogf=-728
             endif
c--------------------------------------------------------------------
c            kinetic energy [KeV]
c--------------------------------------------------------------------
             penergy=(dsqrt(1.d0+(x(j)*vovc)**2)-1.d0)*
     &               (2.d0*m_cs_d_2)                        ![kev]
             if(myrank.eq.0) WRITE(17) REAL(penergy),REAL(plogf)
c            write(*,*) 'ir,j,x(j)**2,dlog(f_analytical),f_analytical)'
c     &               ,ir,j,penergy,dlog(f_analytical),f_analytical)
           enddo !jx

           if(myrank.eq.0) WRITE(17)
        enddo !ir
      endif !(i_first.eq.0) end of initialization

      if (i_first.eq.1) then
c-----------------------------------------------------------------------
c       calculate:
c       3 temperature Maxwellian distribution: f_analytical
c       and the derivative from this disribution by x=P/(m_e*vnorm): df_dx
c-------------------------------------------------------------------------  
        itabl(1)=1
        itabl(2)=0
        itabl(3)=0
 
        call terp1(ir_maxwell,rho_maxwell,ratio_dens,
     &               d2_ratio_dens_drho,rho_in,1,tabl,itabl)

        ratio_dens_spl=tabl(1)

        call Three_Temp_Maxwell_1D(rtem0,rtemp1, rtemp2, rtemp3,
     &  rvtovte1,rvtovte2,x_in,rho_in,f_3T_Maxwell,df_dx)
 
c        write(*,*)'x_in,f_3T_Maxwell,rho_in,df_dx,ratio_dens_spl',
c     &  x_in,f_3T_Maxwell,df_dx,ratio_dens_spl

        f_analytical=f_3T_Maxwell*ratio_dens_spl
        df_dx= df_dx*ratio_dens_spl

        df_dp= df_dx*clight/vnorm


c        write(*,*)'rho_in,temperho(rho_in,1)',rho_in,temperho(rho_in,1)
c        theta_test=mass_e*clight**2/(k_1_kev*temperho(rho_in,1))
c        write(*,*)'theta_test',theta_test
c        write(*,*)'f_analytical*theta_test*p_in/dsqrt(1.d0+p_in**2)',
c     &             (f_analytical*theta_test*p_in/dsqrt(1.d0+p_in**2))
c        write(*,*)'df_dp',df_dp

      endif

      return   
      end


      subroutine Three_Temp_Maxwell_1D(rtem0,
     &rtemp1, rtemp2, rtemp3,
     &rvtovte1,rvtovte2,x,rho,f_3T_Maxwell,df_dx)
c----------------------------------------------------------
      implicit none
c-----input
      double precision rtem0, !relation tem0/electron_temperature(rho=0)
                              !tem0 is the energy for the momentum normalization (KeV) 
     &rtemp1, rtemp2, rtemp3, !relation bin temperature/electron_temperature(rho)
c-----------------------momentum bins boundaries
     &rvtovte1,rvtovte2,!= ratio of momentum-per-mass (electrons) to on-axis
                        !thermal velocity vte0= sqrt(Te/me), defining the
                        !three velocity ranges for the temperatures.
                        !defaults=1.e6,1.e6 [i.e., effectively infinity]
     &x,                ! P/(m_e*vnorm) normalized momemtum per rest electron mass
     &rho               ! Normalized small radius
c-----output
      double precision f_3T_Maxwell, ! Three temperature relativistic Maxwellian 
                                     ! electron distribution. It has not the normalization
                                     ! at the density
     &df_dx                          ! The derivative from three temperature relativistic Maxwellian 
                                     ! electron distribution by the normalized momentum x=P/(me*vnorm).
c-----locals
      double precision pi,rmu,gamma,f_maxwell,tem0,bin_factor,
     &dv,vvol,vovc,vovc2,vnorm

      double precision den,dnormt,bk2 
      
      double precision mass_e,k_1_kev,clight,m_cs_d_2,
     & const_bin_1,const_bin_2,const_bin_3,bin_const,gamma1,gamma2
c-----external besk2as
      double precision temperho,densrho

      pi=4.0*datan(1.d0)

      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !1 egr in 1 KEV
      clight =2.99792458d10     !sm/secC
   
      tem0=rtem0*temperho(0.d0,1)
      vovc2=(2.d0*tem0*k_1_kev/mass_e)/clight**2
 
c      vovc2=2.d0*1.6022d-19*1.d3*tem0/9.1095d-31/2.9979d8**2
      vovc=dsqrt(vovc2)     !v_tem0/c
   
      vnorm=vovc*clight     !cm/sec  V_tem0 =sqrt(2*rtem0*T(rho=0)/m_e)

      gamma1=dsqrt(1.d0+rvtovte1**2/(2.d0*rtem0)*vovc2)
      gamma2=dsqrt(1.d0+rvtovte2**2/(2.d0*rtem0)*vovc2)

c      den=densrho(rho,1)  ! electron density 10**19 m**(-3)
      den=1.d0             ! the distribution will normalized by
                           ! unit density     
      rmu =mass_e*clight**2/(k_1_kev*temperho(rho,1))

c      write(*,*)'in Three_Temp_Maxwell rho,temperho(rho,1),rmu',
c     &rho,temperho(rho,1),rmu

      const_bin_1=0.d0
      const_bin_2=(gamma1-1.d0)*rmu*(1.d0/rtemp2-1.d0/rtemp1)     
      const_bin_3= const_bin_2+
     &             (gamma2-1.d0)*rmu*(1.d0/rtemp3-1.d0/rtemp2)

      if (x.lt.rvtovte1/dsqrt(2.d0*rtem0)) then
c--------first velocity bin
         bin_factor=1.d0/rtemp1
         bin_const=const_bin_1
      else
         if (x.lt.rvtovte2/dsqrt(2.d0*rtem0)) then
c-----------second velosity bin
            bin_factor=1.d0/rtemp2
            bin_const=const_bin_2
         else
c-----------third velocity bin
            bin_factor=1.d0/rtemp3
            bin_const=const_bin_3
         endif
      endif
      gamma=dsqrt(1.d0+x*x*vovc2)          

c      write(*,*)'clight,vnorm',clight,vnorm
c      write(*,*)'clight/vnorm,vnorm/clight',clight/vnorm,vnorm/clight
c      write(*,*)'x,dsqrt(vovc2),x*dsqrt(vovc2),gamma',
c     &x,dsqrt(vovc2),x*dsqrt(vovc2),gamma

c------------------------------------------------------------------------
c     calculation the Mackdonalds function bk2=K_2(theta)*EXP(theta)
c------------------------------------------------------------------------	
      call besk2as(rmu,bk2)

c      f_3T_Maxwell=dexp(-rmu*bin_factor*(gamma-1.d0)+bin_const)
      f_3T_Maxwell=rmu*dexp(-rmu*bin_factor*(gamma-1.d0)+bin_const)/
     &             (4.d0*pi*bk2)                       
      df_dx=f_3T_Maxwell*(-rmu*bin_factor*x/gamma)*vovc2                      
 
c      write(*,*)'in Three_Temp_Maxwell_1D,f_3T_Maxwell',f_3T_Maxwell
c      write(*,*)'rmu,bin_factor,x,gamma',rmu,bin_factor,x,gamma
c      write(*,*)'df_dx=',df_dx
c      write(*,*)'(rmu*x/gamma)*vovc2',(rmu*x/gamma)*vovc2
      return
      end


      subroutine analytical_distrib(theta,p,r,z,phi,
     &f_maxw,d_maxw_dp)
c------------------------------------------------------------   
c     calculates analytical ditribution function and its derivative
c         i_diskf=0  one temperature relativistic Maxwellian distributin
c         i_diskf=5  theree temperature relativistic Maxwellian distributin
c----------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'dskin.i' !      vnorm
      include 'six.i'   !      rhom
c-----input
      double precision theta, ! mc**2/T for 1 Temp maxwell
     &p,                      ! momentum/electron_restmass/clight 
     &r,z,phi                 ! space coordinates 
c-----output
      double precision 
     &f_maxw,     !relativstic maxwellian distribution
                  !nirmalized at density=1.d0 
     &d_maxw_dp   !derivative d_maxw/dp
c-----local

  
c-----expernal
      double precision b,
     &temperho !for test only


      if (i_diskf.eq.0) then
         
         call relativictic_maxwell_1_temperature_analytical(theta,
     &   p,f_maxw,d_maxw_dp)
      endif

      if (i_diskf.eq.5) then

         bmod=b(z,r,phi)!calculates rho
c         write(*,*)'fokker.f before fokker_3_temperature_analytical'
c         write(*,*)'theta',theta
c         write(*,*)'rtem0,rtemp1, rtemp2, rtemp3',
c     &              rtem0,rtemp1, rtemp2, rtemp3
c         write(*,*)'rvtovte1,rvtovte2',rvtovte1,rvtovte2
c         write(*,*)'jx,jxa,vnorm,ndens,ndensa',
c     &              jx,jxa,vnorm,ndens,ndensa
c        write(*,*)'rhom',rhom
c         write(*,*)'p,rho,tempergo(rho,1)',p,rho,temperho(rho,1)

cSAP090808 delete jxa -> jx     
         call fokker_3_temperature_analytical(rtem0,
     &   rtemp1, rtemp2, rtemp3,
     &   rvtovte1,rvtovte2,
cSAP090808 delete jxa -> jx
c    &   jx,jxa,x,vnorm,ndens,ndensa,rhom,
     &   jx,jx,x,vnorm,ndens,ndensa,rhom,
     &   p,rho,
     &   f_maxw,d_maxw_dp)
c         write(*,*)'fokker.f after fokker_3_temperature_analytical'
c         write(*,*)'f_maxw,d_maxw_dp', f_maxw,d_maxw_dp
      endif
          
      return
      end

    
      subroutine relativictic_maxwell_1_temperature_analytical(theta,
     &p,f_maxw,d_maxw_dp)
c------------------------------------------------------------------------
c     calculates relativistic Maxwell distribution function
c     normalized at unit
c     f_maxw=theta/{4pi[K_2(theta)exp(theta)]}exp(theta(1-gamma))
c            gamma=sqrt(1+p**2) , here p=(momentum/mc)
c     and its derivative
c     d_maxw_dp=d_f_maxw/dp
c------------------------------------------------------------------------
      implicit none
c-----input
      double precision theta,  !m_e*clight**2/T
     &p  ! momentum/electron_restmass/clight 
c-----output
      double precision 
     &f_maxw,     !relativstic maxwellian distribution
                  !nirmalized at density=1.d0 
     &d_maxw_dp   !derivative d_maxw/dp

c-----local
      double precision pi,gamma,bk2
 
c      vnorm,       ! momentum/electron_restmass for normalization (cm/sec)
c      x(jxa),      ! momentum/electron_restmass/vnorm mesh:
       
      pi=4*datan(1.d0)
      gamma=dsqrt(1.d0+p*p) !     p=(momentum/mc)=(x*vnorm/clight)
      call besk2as(theta,bk2)
           
      f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)
      d_maxw_dp=p*theta*f_maxw/gamma

      return
      end 


      subroutine test_density_relat_maxwell
      implicit none
      include 'param.i'
      include 'one.i'
      include 'dskin.i'
      include 'six.i'
c-----input
     


      integer ir_maxwell 
      parameter (ir_maxwell=100)
c-----external  besk2as
      double precision temperho,densrho
c-----local
      double precision rho_maxwell(ir_maxwell),ratio_dens(ir_maxwell)
      double precision bk2,step,dv,rho_loc,den,dnormt,vvol,f_maxw,
     *gamma,vovc2,tem0,theta,temp_kev,clight,mass_e,k_1_kev,vnorm_test             
      integer ir,j

      step=1.d0/dble(ir_maxwell-1)       
      do ir=1,ir_maxwell
         rho_maxwell(ir)=step*dble(ir-1)
      enddo

      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !1 egr in 1 KEV
      clight =2.99792458d10     !speed of light [sm/sec]
     
      pi=4.d0*datan(1.d0) 
      tem0=rtem0*temperho(0.d0,1)
      vovc2=(2.d0*tem0*k_1_kev/mass_e)/clight**2
      vnorm_test=dsqrt(2.d0*tem0*k_1_kev/mass_e)

      write(*,*)'fokker.f in test_density_relat_maxwell'
      write(*,*)'ir_maxwell,jx',ir_maxwell,jx
      write(*,*)'dsqrt(vovc2)',dsqrt(vovc2)
      write(*,*)'vnorm,vnorm_test',vnorm,vnorm_test
      write(*,*)'(vnorm/clight)',(vnorm/clight)

cfor test
      step=1.d0/dble(jx-1)
      do j=1,jx
        x(j)=step*(j-1)
      enddo
cend for test

      do ir=1,ir_maxwell
         rho_loc=rho_maxwell(ir)
         den=1.d0
         temp_kev=temperho(rho_loc,1)
         theta=mass_e*clight**2/(k_1_kev*temp_kev)

         write(*,*)'ir,rho_loc,temp_kev,theta',ir,rho_loc,temp_kev,theta

         call besk2as(theta,bk2)
         
         dnormt=0.d0
         do j=1,jx
            gamma=dsqrt(1.d0+x(j)**2*vovc2)
            f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)
            dv=0.d0
            if(j.eq.1) dv=(0.5d0*(x(2)+x(1)))**3/3.d0
            if(j.gt.1.and.j.lt.jx) dv=0.5d0*(x(j+1)-x(j-1))
            vvol=4.d0*pi*x(j)*x(j)*dv              
            vvol=vvol*(vnorm/clight)**3
            dnormt=dnormt+vvol*f_maxw           
        enddo !j

        ratio_dens(ir)=den/dnormt
        write(*,*)'ir,dnormt',ir,dnormt
         
      enddo !ir

      return
      end



