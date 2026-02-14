c        ********************** drkgs    *******************
c         4_th order Runge-Kutta subroutine with constant time 
c         step finds the solution of
c         the system of the ordinary differential equations         
c        ****************************************************
c
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c             Runge-Kutta method  of  4th-order                  !
c          with the constant step integration.                   !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       2. input parameters                                      !
c                                                                !
c          prmt  - the vector of input and output data:          !
c              prmt(1) - the lower boundary of interval of       !
c                        integration;                            !
c              prmt(2) - the upper boundary ...;                 !
c              prmt(3) - the initial step of integration;        !
c              prmt(4) - the accuracy of solution;               !
c              prmt(5) - if it is not equal to zero the control  !
c                        is transferred to the main program;     !
c                                                                !
c          y     - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  at the intermediate points x;                 !
c                                                                !
c          dery  - the  vector of derivatives of y at point x	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array.                          !
c                                                                !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system. it has not to change the  values  of    !
c                x and y. its formal parameters are:             !
c                x, y, dery.                                     !
c                                                                !
c          outp - output subroutine. it has not to change the    !
c                 values of all  its  formal  parameters.  if    !
c                 prmt(5) is not equal to zero,  the  control    !
c                 will be transferred to the main program.       !
c                 its formal parameters are:			 !
c                 x, y, dery, ihlf, ndim, prmt                   !
c          output data						 !
c    u(1)=z,u(2)=r,u(3)=phi,u(4)=nz,u(5)=nr,u(6)=m               !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs(prmt,u,deru,ndim,ihlf,fct,outp,aux)
      implicit none 
c-----input
      integer ndim,ihlf 
      real*8 u(*),deru(*),aux(8,*),prmt(*),t
c-----locals
      real*8 us,uz,ur,uphi,dels1,dels2,up(6),tend,h,dtstep

cSAP120517
      real*8 us_total,dels_total,us_total_old

      integer iraystop,iboundc,ibound,i_go_to_previous_output_step,
     &i,iflagh,ntau
c      integer*2 itb1,itb2,itb3,itb4,itf1,itf2,itf3,itf4,
c     1itd1,itd2,itd3,itd4

cyup      write(*,*)'in Runge -Kutta drkgs'
c----------------------------------------------------------------
c  the program for calculation of the computer time 
c----------------------------------------------------------------
c for sun but not working
c       call gettim(time1)
c---------------------------------------------------------------
c for ibm
c      call gettim(it1,it2,it3,it4)
c      tstart=it1*3600.00d0+it2*60.00d0+it3*1.00d0+it4*0.01d0
c----------------------------------------------------------------
      t=prmt(1)
cyup      write(*,*)'in rk_new t,prmt(1)',t,prmt(1)
      tend=prmt(2)
      h=prmt(3)
      prmt(5)=0.d0
      us=0.d0
cSAP120517
      us_total=0.d0 
      us_total_old=0.d0
c--------------------------------------------------------------
c    ntau is the number of runge-kutta time step. it may be deleted
c    if you do not want to find the calculation of the computer time
      ntau=0
      iflagh=3
c---------------------------------------------------------------
10    continue
cSAP090726
c-----jump in RK through gyro resonace point 
      call jump_through_resonace(us,u,h,up,prmt(7))
      do i=1,ndim
         u(i)=up(i)
      enddo
c--------------------------------------------------
      uz=u(1)
      ur=u(2)
      uphi=u(3)
c      write(*,*)' in rk0 t',t
c      write(*,*)'h=',h
      call fct(t,u,deru)
      do 1 i=1,ndim
         aux(1,i)=h*deru(i)
	 aux(5,i)=u(i)+0.5d0*aux(1,i)
 1	 up(i)=aux(5,i)
c        aux(1,i)=k1(i)
c        ayx(5,i)=y0(i)+0.5*k1(i)

c      write(*,*)'in rk1 u(1-6)'
c      write(*,*)u(1),u(2),u(3)
c      write(*,*)u(4),u(5),u(6)
c      write(*,*)'in rk1 deru(1-6)'
c      write(*,*)deru(1),deru(2),deru(3)
c      write(*,*)deru(4),deru(5),deru(6)
c      write(*,*)'in rk1 up(1-6)'
c      write(*,*)up(1),up(2),up(3)
c      write(*,*)up(4),up(5),up(6)
c------------------------------------------------------------------
c     control that point up(i) is inside the plasma
      call boundc(up(1),up(2),iboundc)
cc      write(*,*)'in rk 1 iboundc',iboundc
      ibound=1
      if (iboundc.eq.1) then
c        ray point is outside the plasma,
c        we change time step h
cc	 write(*,*)'1 u(i),i=1,ndim h=',h
cc	 write(*,*)(u(i),i=1,ndim)
cc	 call boundt(u(1),u(2))
cc	 write(*,*)'deru(i),i=1,ndim'
cc	 write(*,*)(deru(i),i=1,ndim)
cc	 write(*,*)'aux(1,i),i=1,ndim'
cc	 write(*,*)(aux(1,i),i=1,ndim)
cc	 write(*,*)'up(i),i=1,ndim'
cc	 write(*,*)(up(i),i=1,ndim)
cc	 read(*,*)
         iflagh=-1
         h=0.5d0*h
cyup         write(*,*)'in rk1 iboundc=1,h',h
         go to 10
      end if
c------------------------------------------------------------------
      call fct(t+0.5d0*h,up,deru)
      do 2 i=1,ndim
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5d0*aux(2,i)
 2	 up(i)=aux(5,i)
c        aux(2,i)=k2(i)
c        ayx(5,i)=y0(i)+0.5*k2(i)

c      write(*,*)'in rk2 u(1-6)'
c      write(*,*)u(1),u(2),u(3)
c      write(*,*)u(4),u(5),u(6)
c      write(*,*)'in rk2 deru(1-6)'
c      write(*,*)deru(1),deru(2),deru(3)
c      write(*,*)deru(4),deru(5),deru(6)
c      write(*,*)'in rk2 up(1-6)'
c      write(*,*)up(1),up(2),up(3)
c      write(*,*)up(4),up(5),up(6)
c------------------------------------------------------------------
c        control that point up(i) is inside the plasma
	 call boundc(up(1),up(2),iboundc)
      ibound=1
cc      write(*,*)'in rk 2 iboundc',iboundc
	 if (iboundc.eq.1) then
c          ray point is outside the plasma,
c          we change time step h
cc	 write(*,*)'2 u(i),i=1,ndim h=',h
cc	 write(*,*)(u(i),i=1,ndim)
cc	 call boundt(u(1),u(2))
cc	 write(*,*)'deru(i),i=1,ndim'
cc	 write(*,*)(deru(i),i=1,ndim)
cc	 write(*,*)'aux(2,i),i=1,ndim'
cc	 write(*,*)(aux(2,i),i=1,ndim)
cc	 write(*,*)'up(i),i=1,ndim'
cc	 write(*,*)(up(i),i=1,ndim)
cc	 read(*,*)
	   iflagh=-1
	   h=0.5d0*h
c	   write(*,*)'in rk2 iboundc=1,h',h
	   go to 10
	 end if
c--------------------------------------------------------------------
      call fct(t+0.5d0*h,up,deru)
      do 3 i=1,ndim
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
 3	 up(i)=aux(5,i)
c        aux(3,i)=k3(i)
c        ayx(5,i)=y0(i)+k3(i)

c      write(*,*)'in rk3 u(1-6)'
c      write(*,*)u(1),u(2),u(3)
c      write(*,*)u(4),u(5),u(6)
c      write(*,*)'in rk3 deru(1-6)'
c      write(*,*)deru(1),deru(2),deru(3)
c      write(*,*)deru(4),deru(5),deru(6)
c      write(*,*)'in rk3 up(1-6)'
c      write(*,*)up(1),up(2),up(3)
c      write(*,*)up(4),up(5),up(6)
c------------------------------------------------------------------
c        control that point up(i) is inside the plasma
	 call boundc(up(1),up(2),iboundc)
cc      write(*,*)'in rk 3 iboundc',iboundc
	 if (iboundc.eq.1) then
      ibound=1
c          ray point is outside the plasma,
c          we change time step h
cc	 write(*,*)'3 u(i),i=1,ndim h=',h
cc	 write(*,*)(u(i),i=1,ndim)
cc	 call boundt(u(1),u(2))
cc	 write(*,*)'deru(i),i=1,ndim'
cc	 write(*,*)(deru(i),i=1,ndim)
cc	 write(*,*)'aux(3,i),i=1,ndim'
cc	 write(*,*)(aux(3,i),i=1,ndim)
cc	 write(*,*)'up(i),i=1,ndim'
cc	 write(*,*)(up(i),i=1,ndim)
cc	 read(*,*)
	   iflagh=-1
	   h=0.5d0*h
c	   write(*,*)'in rk3 iboundc=1,h',h
	   go to 10
	 end if
c--------------------------------------------------------------------
      call fct(t+h,up,deru)
      do 4 i=1,ndim
 4       aux(4,i)=h*deru(i)
c        aux(4,i)=k4(i)
c--------------------------------------------------------------------
      do 5 i=1,ndim
 5       up(i)=u(i)+1.d0/6.d0*(aux(1,i)+2.d0*aux(2,i)+2.d0*aux(3,i)
     1          +aux(4,i))

c      write(*,*)'in rk5 u(1-6)'
c      write(*,*)u(1),u(2),u(3)
c      write(*,*)u(4),u(5),u(6)
c      write(*,*)'in rk5 deru(1-6)'
c      write(*,*)deru(1),deru(2),deru(3)
c      write(*,*)deru(4),deru(5),deru(6)
c      write(*,*)'in rk5 up(1-6)'
c      write(*,*)up(1),up(2),up(3)
c      write(*,*)up(4),up(5),up(6)
c------------------------------------------------------------------
c        control that point up(i) is inside the plasma
	 call boundc(up(1),up(2),iboundc)
cc      write(*,*)'in rk 4 iboundc',iboundc
	 if (iboundc.eq.1) then
c          ray point is outside the plasma,
c          we change time step h
cc	 write(*,*)'4 u(i),i=1,ndim h=',h
cc	 write(*,*)(u(i),i=1,ndim)
cc	 call boundt(u(1),u(2))
cc	 write(*,*)'deru(i),i=1,ndim'
cc	 write(*,*)(deru(i),i=1,ndim)
cc	 write(*,*)'aux(4,i),i=1,ndim'
cc	 write(*,*)(aux(4,i),i=1,ndim)
cc	 write(*,*)'up(i),i=1,ndim'
cc	 write(*,*)(up(i),i=1,ndim)
cc	 read(*,*)
	   iflagh=-1
	   h=0.5d0*h
	   write(*,*)'in rk4 iboundc=1,h',h
	   go to 10
	 end if
c--------------------------------------------------------------------
c      do 51 i=1,ndim
c 51       u(i)=up(i)
c--------------------------------------------------------------------
c      write(*,*)'in rk0 before t=t+h t,h',t,h
      t=t+h
c      write(*,*)'in rk0 after t=t+h t',t
c     h=prmt(3)
c--------------------------------------------------------------------
      dels2=(uz-up(1))*(uz-up(1))+(ur-up(2))*(ur-up(2))
      dels1=dsqrt(dels2)
      us=us+dels1
c      write(*,*)'in rk0 t,h',t,h
c      write(*,*)'in rk0 us,dels1',us,dels1
cSAP120517 for Nicola scattering
      dels_total=dsqrt((uz-u(1))**2+(ur-u(2))**2
     &                 +(ur*uphi-u(2)*u(3))**2)
      us_total_old=us_total
      us_total=us_total+dels_total
      call outp(us,up,deru,ihlf,ndim,prmt,iflagh,iraystop,
     &i_go_to_previous_output_step,us_total,us_total_old )
c--------------------------------------------------------------------
      if (iraystop.eq.1) then
        goto 100
      end if
c      write(*,*)'in rk0 after outpt  t',t
c     write(*,*)'in rk after outpt iflagh=',iflagh
      if (iflagh.eq.1) then
c if 1
c----------------------------------------------------------------
c        ray is near the plasma boundary after Runge-Kutta
c        procedure (it is after reflection)
c-----------------------------------------------------------------
         h=prmt(3)
	do 54 i=1,ndim
54       u(i)=up(i)
	write(*,*)'in rk after output  iflagh=',iflagh,'h=',h
      else
        if (iflagh.eq.2) then
c if 2
c-------------------------------------------------------------------
c         Ray is outside  the plasma after the correction.
c         It is necessary to reduce the time step in the Runge-Kutta
c         procedure and to recalculate the array u(i)
c----------------------------------------------------------------------
          t=t-h
          h=h*0.5d0
	  dtstep=prmt(3)/h
c------------------------------------------------------------------
	  if (dtstep.lt.10000.d0)then
c if 2.1
c            write(*,*)'in rk after output  iflagh, h=',iflagh,h
	    goto 10
          else
cyup            write(*,*)'dtstep.gt.10000.0'
c previous Runge-Kutta time step is the reflection point
c            call reflect(u(1),u(2),u(3),u(4),u(5),cnzref,cnrref)
c	    u(4)=cnzref
c	    u(5)=cnrref
c	    irefl=irefl+1
c	    write(*,*)'reflection point1 in rk'

	    goto 10
	  endif
c end if 2.1
c-------------------------------------------------------------------
        else
c-------------------------------------------------------------------
c     ordinary ray point	  iflagh=3
c     value of u(i) after correction procedure
            h=prmt(3)
           do 55 i=1,ndim
55            u(i)=up(i)
c-------------------------------------------------------------------
        endif
c end if 2
      endif
c end if 1
c--------------------------------------------------------------------
      do 51 i=1,ndim
51       u(i)=up(i)
c--------------------------------------------------------------------
      if (iraystop.eq.1) then
        goto 100
      end if
      if (prmt(5).gt.0) goto 100
c------------------------------------------------------------------
c  the program for timing determination
c      ntau=ntau+1
c      if (ntau.eq.10) then
c        call gettim(itf1,itf2,itf3,itf4)
c	itd1=itf1-it1
c	itd2=itf2-it2
c	itd3=itf3-it3
c	itd4=itf4-it4
c	tstart=it1*3600.00d0+it2*60.00d0+it3*1.00d0+it4*0.01d0
c	tfinish=itf1*3600.00d0+itf2*60.00d0+itf3*1.00d0+itf4*0.01d0
c	deltim=tfinish-tstart
c	write(*,*)'deltim=',deltim
c	write(*,*)'it1=',it1,'it2=',it2,'itd3=',it3,'itd4=',it4
c	write(*,*)'itf1=',itf1,'itf2=',itf2,'itf3=',itf3,'itf4=',itf4
c	write(*,*)'itd1=',itd1,'itd2=',itd2,'itd3=',itd3,'itd4=',itd4
c        call gettim(time2)
c        time=time2-time1
c        write(*,*) ' 2 - calculation time for 100 step:',time
c      end if
c--------------------------------------------------------
      goto 10
  100 continue
      return
      end

c        4_th order method Runge-Kutta with automatical step select 
c        The time step can be reduced only.
c
c        ********************** drkgs3     ******************
c        *                      -----                       
c        * this Runge-Kutta subroutine finds the solution of the system 
c        *      of ordinary differential equations          
c        ****************************************************
c
c
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c             runge-kutta method  of  4th-order                  !
c          with constant step inegration.                        !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       2. input parameters                                      !
c                                                                !
c          prmt  - the vector of input and output data:          !
c              prmt(1) - the lower boundary of interval of       !
c                        integration;                            !
c              prmt(2) - the upper boundary ...;                 !
c              prmt(3) - the initial step of integration;        !
c              prmt(4) - the accuracy of solution;               !
c              prmt(5) - if it is not equal to zero the control  !
c                        is transferred to the main program;     !
c                                                                !
c          y     - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  at the intermediate points x;                 !
c                                                                !
c          dery  - the  vector of derivatives of y at point x	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array.                          !
c                                                                !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system. it has not to change the  values  of    !
c                x and y. its formal parameters are:             !
c                x, y, dery.                                     !
c                                                                !
c          outp - output subroutine. it has not to change the    !
c                 values of all  its  formal  parameters.  if    !
c                 prmt(5) is not equal to zero,  the  control    !
c                 will be transferred to the main program.       !
c                 its formal parameters are:			 !
c                 x, y, dery, ihlf, ndim, prmt                   !
c          output date						 !
c    u(1)=z,u(2)=r,u(3)=phi,u(4)=nz,u(5)=nr,u(6)=m               !                                         !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs3(prmt,u,deru,ndim,ihlf,fct,outp,aux)
      dimension u(*),deru(*),aux(8,*),prmt(*),up(6),uu(6),startu(6)
      double precision prmt,u,deru,aux,t,tend,h,up,tt,hh,uu,startu,dd
c      integer*2 itb1,itb2,itb3,itb4,itf1,itf2,itf3,itf4,
c     1itd1,itd2,itd3,itd4
cSAP120517
      real*8 us_total,dels_total,us_total_old

cyup      write(*,*)'in Runge -Kutta !!!'
c----------------------------------------------------------------
c  the program for calculation time determination
c----------------------------------------------------------------
c for sun but not working
c       call gettim(time1)
c---------------------------------------------------------------
c for ibm
c      call gettim(it1,it2,it3,it4)
c      tstart=it1*3600.00d0+it2*60.00d0+it3*1.00d0+it4*0.01d0
c----------------------------------------------------------------
      t=prmt(1)
      tt=t
      tend=prmt(2)
      h=prmt(3)
      hh=h/2.0d0
      prmt(5)=0.d0
      eps=1.0d-10
cSAP120517
      us_total=0.d0 
      us_total_old=0.d0
c---------------------------------------------------------------
c    ntau is the number of runge-kutta time step.it may be deleted
c    if you do not want to find the calculation time
      ntau=0
c---------------------------------------------------------------
 10   continue
      do i=1,ndim
	uu(i)=u(i)
	startu(i)=u(i)
      enddo

      call fct(t,u,deru)
      do i=1,ndim
         aux(1,i)=h*deru(i)
	 aux(5,i)=u(i)+0.5*aux(1,i)
 	 up(i)=aux(5,i)
      enddo
c        aux(1,i)=k1(i)
c        ayx(5,i)=y0(i)+0.5*k1(i)
c------------------------------------------------------------------
      call fct(t+0.5*h,up,deru)
      do i=1,ndim
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5*aux(2,i)
 	 up(i)=aux(5,i)
      enddo
c        aux(2,i)=k2(i)
c        ayx(5,i)=y0(i)+0.5*k2(i)
c--------------------------------------------------------------------
      call fct(t+0.5*h,up,deru)
      do i=1,ndim
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
 	 up(i)=aux(5,i)
      enddo
c        aux(3,i)=k3(i)
c        ayx(5,i)=y0(i)+k3(i)
c--------------------------------------------------------------------
      call fct(t+h,up,deru)
      do i=1,ndim
        aux(4,i)=h*deru(i)
      enddo
c        aux(4,i)=phi3(i)
c--------------------------------------------------------------------
      do i=1,ndim
        u(i)=u(i)+1./6.*(aux(1,i)+2.*aux(2,i)+2.*aux(3,i)+aux(4,i))
      enddo
c--------------------------------------------------------------------
c=====================================================================
      do 320 j=1,2

      call fct(tt,uu,deru)
      do i=1,ndim
         aux(1,i)=hh*deru(i)
	 aux(5,i)=uu(i)+0.5*aux(1,i)
 	 up(i)=aux(5,i)
      enddo
c        aux(1,i)=k1(i)
c        ayx(5,i)=y0(i)+0.5*k1(i)
c------------------------------------------------------------------
      call fct(tt+0.5*hh,up,deru)
      do i=1,ndim
         aux(2,i)=hh*deru(i)
         aux(5,i)=uu(i)+0.5*aux(2,i)
 	 up(i)=aux(5,i)
      enddo
c        aux(2,i)=k2(i)
c        ayx(5,i)=y0(i)+0.5*k2(i)
c--------------------------------------------------------------------
      call fct(tt+0.5*hh,up,deru)
      do i=1,ndim
         aux(3,i)=hh*deru(i)
         aux(5,i)=uu(i)+aux(3,i)
 	 up(i)=aux(5,i)
      enddo
c        aux(3,i)=k3(i)
c        ayx(5,i)=y0(i)+k3(i)
c--------------------------------------------------------------------
      call fct(tt+hh,up,deru)
      do i=1,ndim
        aux(4,i)=hh*deru(i)
      enddo
c        aux(4,i)=phi3(i)
c--------------------------------------------------------------------
      do i=1,ndim
        uu(i)=uu(i)+1./6.*(aux(1,i)+2.*aux(2,i)+2.*aux(3,i)+aux(4,i))
      enddo
      tt=tt+hh
c--------------------------------------------------------------------
 320   continue
c=====================================================================
      dd=0.0d0
      do i=1,ndim
         dd=dd+(u(i)-uu(i))*(u(i)-uu(i))
      enddo
      dd=dsqrt(dd)
c      write(*,*)'test rk_new! eps=',dd
      if (dd .lt. eps) goto 191

      h=hh
      hh=h/2.0d0
      tt=t
      do i=1,ndim
        u(i)=startu(i)
      enddo
c      write(*,*)'!!!!!!!!!!!!!!!! back !!!!!!!!!!!!!!!'
      goto 10


c=====================================================================
191	continue
          t=t+h
	  uout1=u(1)
	  uout2=u(2)
	  uout3=u(3)
	  uout4=u(4)
	  uout5=u(5)
	  uout6=u(6)
c-----------------------------------------------------------
cSAP120517 for Nicola scattering
      dels_total=dsqrt((uz-u(1))**2+(ur-u(2))**2
     &                 +(ur*uphi-u(2)*u(3))**2)
      us_total_old=us_total
      us_total=us_total+dels_total

      call outp(t,u,deru,ihlf,ndim,prmt,iflagh,iraystop,
cSAP100202
c     & i_go_to_previous_output_step)
cSAP120517
     &i_go_to_previous_output_step,us_total,us_total_old )
      if (iraystop.eq.1) then
        goto 100
      end if
      if (prmt(5).gt.0) goto 100
  30      format(1x,6(1x,e11.4))
c------------------------------------------------------------------
c  the program for calculation time determination
      ntau=ntau+1
      if (ntau.eq.10) then
c        call gettim(itf1,itf2,itf3,itf4)
c	itd1=itf1-it1
c	itd2=itf2-it2
c	itd3=itf3-it3
c	itd4=itf4-it4
c	tstart=it1*3600.00d0+it2*60.00d0+it3*1.00d0+it4*0.01d0
c	tfinish=itf1*3600.00d0+itf2*60.00d0+itf3*1.00d0+itf4*0.01d0
c	deltim=tfinish-tstart
c	write(*,*)'deltim=',deltim
c	write(*,*)'it1=',it1,'it2=',it2,'itd3=',it3,'itd4=',it4
c	write(*,*)'itf1=',itf1,'itf2=',itf2,'itf3=',itf3,'itf4=',itf4
c	write(*,*)'itd1=',itd1,'itd2=',itd2,'itd3=',itd3,'itd4=',itd4
c        call gettim(time2)
c        time=time2-time1
c        write(*,*) ' 2 - calculation time for 100 step:',time
      end if
c--------------------------------------------------------
      goto 10
  100 continue
      return
      end



c        4_th order  Runge-Kutta method with automatic 
c        time step selection
c
c        **********************     drkgs2   ****************
c        *                      -----                       
c        * this Runge-Kutta subroutine finds the solution of the 
c        *      system of ordinary differential equations          
c        ****************************************************
c
c
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c             runge-kutta method  of  4th-order                  !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       2. input parameters                                      !
c                                                                !
c          prmt  - the vector of input and output data:          !
c              prmt(1) - the lower boundary of interval of       !
c                        integration;                            !
c              prmt(2) - the upper boundary ...;                 !
c              prmt(3) - the initial step of integration;        !
c              prmt(4) - the accuracy of solution;               !
c              prmt(5) - if it is not equal to zero the control  !
c                        is transferred to the main program;     !
c                                                                !
c          y     - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  at the intermediate points x;                 !
c                                                                !
c          dery  - the  vector of derivatives of y at point x	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array.                          !
c                                                                !
c          i_output =1 the output at the poloidal distance steps !
c                   =2 the output at the total distance steps    !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system. it has not to change the  values  of    !
c                x and y. its formal parameters are:             !
c                x, y, dery.                                     !
c                                                                !
c          outp - output subroutine. it has not to change the    !
c                 values of all  its  formal  parameters.  if    !
c                 prmt(5) is not equal to zero,  the  control    !
c                 will be transferred to the main program.       !
c                 its formal parameters are:			 !
c                 x, y, dery, ihlf, ndim, prmt                   !
c          output date						 !
c    u(1)=z,u(2)=r,u(3)=phi,u(4)=nz,u(5)=nr,u(6)=m               !                                         !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs2(prmt,u,deru,ndim,ihlf,fct,outp,aux,i_output,
     & nbulk,nharm_refined_step) ![2025-12-16] added nbulk,nharm_refined_step
      implicit none
cSAP100514
      include 'param.i'
      include 'write.i'
c-----input

      dimension u(6),deru(6),aux(8,6),prmt(9)
      real*8 u,deru,aux,prmt
      integer i_output,ndim,ihlf
    
c-----locals
      dimension up(6),uu(6),startu(6),uplus(6)
      real*8 up,uu,startu,uplus
      real*8 t,tend,h,tt,hh,
     & dd,eps_prmt4
      real*8 us,uz,ur,dels1,dels2,uphi,
     &uout1,uout2,uout3,uout4,uout5,uout6,
     &usold,dp,delh,ham,dtstep,eps_ham,hnew,hold
     
      real*8 yi
      real*8 y !external
      real*8 prmt4_orig, prmt4_adj, prmt6_orig, prmt6_adj !original and adjusted
      integer nbulk,nharm_refined_step !INPUT
      integer ic_near,nh !local

cSAP120517
      real*8 us_total,dels_total,us_total_old
      integer iflag,i,j,ntau,i_delh,i_delh_1,iboundc,iflagh,iraystop,
     &i_go_to_previous_output_step
cSm030515
       external fct
c      integer*2 itb1,itb2,itb3,itb4,itf1,itf2,itf3,itf4,
c     1itd1,itd2,itd3,itd4
cyup      write(*,*)'in Runge -Kutta drkgs2!!! prmt',prmt
c----------------------------------------------------------------
c  the program for calculation time determination
c----------------------------------------------------------------
c for sun but not working
c       call gettim(time1)
c---------------------------------------------------------------
c for ibm
c      call gettim(it1,it2,it3,it4)
c      tstart=it1*3600.00d0+it2*60.00d0+it3*1.00d0+it4*0.01d0
c----------------------------------------------------------------

      prmt4_orig= prmt(4) !Save original value
      prmt4_adj=  prmt(4)
      prmt6_orig= prmt(6) !Save original value
      prmt6_adj=  prmt(6)
      
      t=prmt(1)
      tt=t
      tend=prmt(2)
      h=prmt(3)
      hh=h/2.0
      prmt(5)=0.d0 !if =1.d0, it will stop the ray after 1 step
      eps_prmt4=prmt(4)
      iflag=0  ! control of the RK method accuracy
      iflagh=3 ! control of the plasma edg intersecsion
      us=0.d0  ! poloidal length (initial value)
cSAP120517
      us_total=0.d0 
      us_total_old=0.d0
c---------------------------------------------------------------
c    ntau is the number of runge-kutta time step.it may be delitted
c    if you do not want to find the calculation time
      ntau=0
c---------------------------------------------------------------
 10   continue ! handle for time advances (in a loop)

c      write(*,*)'R-K 10 h',h,'t',t
cSAP090726
c-----jump in RK through gyro resonace point 
c      call jump_through_resonace(us,u,h,up,prmt(7))   
c      do i=1,ndim
c         u(i)=up(i)
c      enddo
c      write(*,*)'after jump prmt(7)',prmt(7)
c-------------------------------------------------
cSAP090228
cSAP090502

      ![2025-12-15] Adjust integration step near resonance.
      ic_near=0 ! To be changed if ray element is near resonance.
      !Restore to original (from genray namelist):
      prmt4_adj=prmt4_orig 
      prmt6_adj=prmt6_orig 
      prmt(4)=prmt4_orig
      prmt(6)=prmt6_orig
      if(nharm_refined_step>0)then
      uz=u(1)
      ur=u(2)
      uphi=u(3)
      do i=1,nbulk !scan species
        !If at least one species has nearly res.,
        ! change ic_near to 1 
        if(ic_near.eq.0)then
          yi= dabs(y(uz,ur,uphi,i)) != |omega_ci/omega|
          do nh=1,nharm_refined_step !scan harmonics
            if( dabs(1.d0-nh*yi) .le. 0.1d0 ) then
              ic_near=1 !flag for proximity of resonance
            endif
          enddo !nh
        endif
      enddo !i=1,nbulk
      if(ic_near.eq.1)then !Near resonance. Reduce steps:
        prmt4_adj=prmt4_orig*0.1 !Reduce by 10x (hardwired, for now)
        prmt6_adj=prmt6_orig*0.1 !Reduce by 10x (hardwired, for now)
        prmt(4)=prmt4_adj
        prmt(6)=prmt6_adj
      endif 
      eps_prmt4=prmt(4)
      endif !(nharm_refined_step>0)
      ![2025-12-15] done

      if (dabs(h).lt.1.d-11) then
         write(*,*)'***** In Runge-Kutta subroutine drkgs2 *****'
         write(*,*)'***** time step is too small prmt3<1.e-11  *****' 
         write(*,*)'***** Cannot achieve given accuracy prmt4 or prmt9'
         write(*,*)'h=prmt3=',h
         write(*,*)'  iraystop->1'
         iraystop=1 ! Cannot achieve given accuracy 
         iray_status_one_ray=13 ! dabs(h).lt.1.d-11 !in subroutine drkgs2
         prmt(4)=prmt4_orig ![2025-12-15] reset to original (in case it was adjusted)
         prmt(6)=prmt6_orig ![2025-12-15] reset to original (in case it was adjusted)
         goto 100 ! exit
      endif    

c       write(*,*)'R-K 10 h',h,'us',us
ctest 
c      call calc_hamiltonian(u,hamiltonian)
c      write(*,*)'in rk_new 0 hamiltonian=',hamiltonian
c      write(*,*)'u',u
cedtest

      uz=u(1)
      ur=u(2)
      uphi=u(3)
      do i=1,ndim
	uu(i)=u(i)
	startu(i)=u(i)
      enddo
      tt=t
      hh=h/2.0d0
c      write(*,*)'R-K 1 u',u
      call fct(t,u,deru)
c      write(*,*)'rk after fct1 deru' 
c      write(*,*)(deru(i),i=1,ndim)
c      write(*,*)'rk after fct1 deru(1),deru(2)',deru(1),deru(2)
      do i=1,ndim	       
         aux(1,i)=h*deru(i)
	 aux(5,i)=u(i)+0.5*aux(1,i)
 	 up(i)=aux(5,i)
      enddo
c      write(*,*)'rk after fct1 up'
c      write(*,*)(up(i),i=1,ndim)
c        aux(1,i)=k1(i)
c        ayx(5,i)=y0(i)+0.5*k1(i)
c------------------------------------------------------------------
      call fct(t+0.5*h,up,deru)
      do i=1,ndim
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5*aux(2,i)
 	 up(i)=aux(5,i)
      enddo
c      write(*,*)'rk after fct2 up h=',h
c      write(*,*)(u(i),i=1,ndim)
c      write(*,*)(deru(i),i=1,ndim)
c      write(*,*)(up(i),i=1,ndim)
c        aux(2,i)=k2(i)
c        ayx(5,i)=y0(i)+0.5*k2(i)
c------------------------------------------------------------------
cc      write(*,*)'rk before fct3 u'
c      write(*,*)u

      call fct(t+0.5*h,up,deru)
c      write(*,*)'rk after fct3 deru'
cc      write(*,*)(deru(i),i=1,ndim)
      do i=1,ndim
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
 	 up(i)=aux(5,i)
      enddo
c      write(*,*)'rk after fct3 up'
c      write(*,*)(up(i),i=1,ndim)
c        aux(3,i)=k3(i)
c        ayx(5,i)=y0(i)+k3(i)
c------------------------------------------------------------------
c      write(*,*)'rk before fct4 u'
c       write(*,*)u
      call fct(t+h,up,deru)
c      write(*,*)'rk after fct4 deru'
cc      write(*,*)(deru(i),i=1,ndim)
      do i=1,ndim
        aux(4,i)=h*deru(i)
      enddo
c        aux(4,i)=phi3(i)
c--------------------------------------------------------------------
      do i=1,ndim
        up(i)=u(i)+1.d0/6.d0*
     1        (aux(1,i)+2.d0*aux(2,i)+2.d0*aux(3,i)+aux(4,i))
      enddo
c      write(*,*)'rk after fct4 up'
c      write(*,*)(up(i),i=1,ndim)
c--------------------------------------------------------------------
      do i=1,ndim
         u(i)=up(i)
      enddo

ctest 
c      call calc_hamiltonian(u,hamiltonian)
c      write(*,*)'in rk_new before 320 hamiltonian=',hamiltonian
c      write(*,*)'u',u  
cedtest

c=====================================================================
      !for accuracy control: two steps starting with old uu 
      !but using h/2
      do 320 j=1,2
      ! u->uu was saved before step 1
      call fct(tt,uu,deru)
c      write(*,*)'rk 320 after fct1 deru',deru

      do i=1,ndim
         aux(1,i)=hh*deru(i)
	 aux(5,i)=uu(i)+0.5*aux(1,i)
 	 up(i)=aux(5,i)
      enddo
c        aux(1,i)=k1(i)
c        ayx(5,i)=y0(i)+0.5*k1(i)
c--------------------------------------------------------------------
c      write(*,*)'bef fct2 up',up
      call fct(tt+0.5*hh,up,deru)
c      write(*,*)'rk 320 after fct2 deru',deru
      do i=1,ndim
         aux(2,i)=hh*deru(i)
         aux(5,i)=uu(i)+0.5*aux(2,i)
 	 up(i)=aux(5,i)
      enddo
c        aux(2,i)=k2(i)
c        ayx(5,i)=y0(i)+0.5*k2(i)
c--------------------------------------------------------------------
c      write(*,*)'bef fct3 up',up
      call fct(tt+0.5*hh,up,deru)
c      write(*,*)'rk 320 after fct3 deru',deru
      do i=1,ndim
         aux(3,i)=hh*deru(i)
         aux(5,i)=uu(i)+aux(3,i)
 	 up(i)=aux(5,i)
      enddo
c        aux(3,i)=k3(i)
c        ayx(5,i)=y0(i)+k3(i)
c--------------------------------------------------------------------
c      write(*,*)'bef fct3 up',up
      call fct(tt+hh,up,deru)
c      write(*,*)'rk 320 after fct4 deru',deru
      do i=1,ndim
        aux(4,i)=hh*deru(i)
      enddo
c        aux(4,i)=phi3(i)
c--------------------------------------------------------------------
      do i=1,ndim
        uu(i)=uu(i)+1./6.*(aux(1,i)+2.*aux(2,i)+2.*aux(3,i)+aux(4,i))
      enddo
c      write(*,*)'bef 320 continue uu',uu 
c--------------------------------------------------------------------
      tt=tt+hh
c--------------------------------------------------------------------
 320  continue
      ! Done: uu() after 2steps using h/2 is obtained
!So, (uz,ur,uphi) is before step, u(i) - after 1 step with h, 
!                               uu(i) - after 2 steps with h/2

c=====================================================================
      dd=0.0d0
      
      do 80 i=1,ndim

c        if(i.eq.3) goto 80
cc--------relative error
         dp=dabs(u(i))
         if (dp.gt.1.d0) then
            dp=1.d0/dp
         else
            dp=1.d0
         endif
         dd=dd+(u(i)-uu(i))*(u(i)-uu(i))*dp
cc--------absolute error
c         dd=dd+(u(i)-uu(i))*(u(i)-uu(i))

 80   continue  
      
      dd=dsqrt(dd)
c       write(*,*)'2 test rk eps_prmt4=',dd
c       write(*,*)'test rk iflag=',iflag

cSAP080807
      eps_ham =prmt(9)
      call callc_eps_ham(uu,ham)
c      write(*,*)'ham,eps_ham',ham,eps_ham
     
c      write(*,*)'dd,eps_prmt4',dd,eps_prmt4
c      write(*,*)'dabs(ham),eps_ham',dabs(ham),eps_ham
      if ((dd .lt. eps_prmt4).and.(dabs(ham).lt.eps_ham)) goto 189
c      if (dd .lt. eps_prmt4) goto 189
        
c      write(*,*)'iflag',iflag
 
      if (iflag .le. 0) then
        h=hh
        do i=1,ndim
           u(i)=startu(i)
        enddo
c         write(*,*)'!!!!!!!!!!!!!!!! back !!!!!!!!!!!!!!!'
        iflag=-1
        !write(*,*)'drkgs2_xyz: dd>eps_prmt4  goto 10',dd,eps_prmt4
        goto 10 !another step
      else
        h=hh
        do i=1,ndim
           u(i)=uplus(i)
        enddo
c         write(*,*)'!!!!!!!!!!!!!!!! done !!!!!!!!!!!!!!!'
        iflag=2
      endif

 189  continue
c       write(*,*)'after 189 iflagh',iflagh 
c      if (iflagh.eq.-1) goto 191
       if (iflagh.eq.2) goto 191
c       write(*,*)'test rk 189 iflag=',iflag
      if (iflag .lt. 0) goto 191

      if ((dd.lt.eps_prmt4/16.d0) .and. (dabs(ham).lt.prmt(9))) then
        ! The time step might be too small (dd can be 0)
        ! Double the time step and try again:
        h=h*2.0d0
        do i=1,ndim
           uplus(i)=u(i)
        enddo
        do i=1,ndim
           u(i)=startu(i)
        enddo
c         write(*,*)'!!!!!!!!!!!!!!!! back++ !!!!!!!!!!!!!!!'
        iflag=1
        !write(*,*)'dd<eps_prmt4/16  h->2h, goto 10'
        goto 10 !another step
      endif

c=====================================================================

191   continue

c      t=t+h
      iflag=0
      uout1=u(1)
      uout2=u(2)
      uout3=u(3)
      uout4=u(4)
      uout5=u(5)
      uout6=u(6)
c-----------------------------------------------------------
c      write(*,*)'in rk before outpt deru(1),deru(2)',deru(1),deru(2)

      if (i_output.eq.1) then
c        poloidal distance
         dels2=(uz-u(1))*(uz-u(1))+(ur-u(2))*(ur-u(2))
      endif

      if (i_output.eq.2) then
c       total distance
         dels2=(uz-u(1))**2+(ur-u(2))**2+(ur*uphi-u(2)*u(3))**2
      endif
    
      dels1=dsqrt(dels2)
      
      usold=us
 
      us=us+dels1 ! the distance after the Runge-Kutta with the variable time step  
cSAP120517 for Nicola scattering
      dels_total=dsqrt((uz-u(1))**2+(ur-u(2))**2
     &                 +(ur*uphi-u(2)*u(3))**2)
      us_total_old=us_total
      us_total=us_total+dels_total

cyup       write(*,*)'rk_new before outpt usold',usold,'us',us
c     &           'dels1',dels1
c      write(*,*)'rk prmt(6),delh',prmt(6),delh
c      write(*,*)'rk prmt(7),prmt(7)+delh*prmt(6)',
c     &prmt(7),prmt(7)+delh*prmt(6)

c      goto 60
c-----check that the ray point is close to the output point
      delh=1.d-1     
 
      i_delh_1=0
 20   continue
      i_delh=0

      if(us.gt.(prmt(7)+delh*prmt6_adj)) then
c-------trajectory jumped past the output point us > prmt(7)+delh*prmt(6)
c-------We will reduse the time step to get the point close to the given output point 
        i_delh=1
        hold=h 
 
        if (i_delh_1.ne.1) hnew=h*(prmt(7)+delh*prmt6_adj-usold)/(dels1)
        i_delh_1=1

c-------one step using the Runge-Kutta with the constant time step h=hnew

c       go back (to the previous time step) to start of the Runge-Kutta procedure 
     
        do i=1,ndim
	   u(i)=startu(i)
        enddo
c------ one step Runge-Kutta procedure    
cyup        write(*,*)'rk_new before drkgs0 hnew,u(1),u(2)',hnew,u(1),u(2)
        call drkgs0(hnew,u,deru,ndim,fct,aux)
c        write(*,*)'0 rk_new after drkgs0 u(1),u(2)',u(1),u(2)

        call boundc(u(1),u(2),iboundc )
       
c        write(*,*)'rk_new.f after drkgs0 iboundc',iboundc       
      endif
c-----------------------------------------------------------
c      write(*,*)'if i_output.eq.1 poloidal dist eq.2 total',i_output
      if (i_output.eq.1) then
c        poloidal distance after one step using drkgs0
         dels2=(uz-u(1))*(uz-u(1))+(ur-u(2))*(ur-u(2))
      endif

      if (i_output.eq.2) then
c        total distance
         dels2=(uz-u(1))**2+(ur-u(2))**2+(ur*uphi-u(2)*u(3))**2
      endif

      dels1=dsqrt(dels2)

      us=usold
      us=us+dels1 ! the distance after one step with h=hnew

cSAP120517 for Nicola scattering
      dels_total=dsqrt((uz-u(1))**2+(ur-u(2))**2
     &                 +(ur*uphi-u(2)*u(3))**2)
      us_total_old=us_total
      us_total=us_total+dels_total ! the total distance after
                                   ! one step with h=hnew
cyup      write(*,*)'i_delh,usold,dels1,us',i_delh,usold,dels1,us
 
c      write(*,*)'rk 2 bef call output us,dels1',us,dels1
c      write(*,*)'ur,u(2)',ur,u(2)
c      write(*,*)'uz,u(1)',uz,u(1)

      if (i_delh.eq.0) goto 60
      goto 60  
      if(us.lt.(prmt(7)+delh*prmt6_adj)) then 
c--------hnew is too short to jump past prmt(7)+delh*prmt(6).
c        Now we will increase hnew.
          hnew=0.5d0*(hold+hnew)
          goto 20      
       endif
 60   continue
cyup      write(*,*)'rk_new before outp us,u(1),u(2)',us,u(1),u(2)
      call outp(us,u,deru,ihlf,ndim,prmt,iflagh,iraystop,
     &i_go_to_previous_output_step,us_total,us_total_old )
cyup      write(*,*)'rk_new.f in drkgs2 after outp us',us

       if (i_go_to_previous_output_step.eq.1) then
cyup          write(*,*)'in rk_new.f i_go_to_previous_output_step',
cyup     &                           i_go_to_previous_output_step
cyup          write(*,*)'in rk_new.f us',us
          goto 10 !another step
       endif
c       write(*,*)'rk_new aft outp iraystop',iraystop
c-----------------------------------------------------------
c      write(*,*)'in rk after outpt,iflagh',iflagh
      if (iflagh.eq.1) then
c if 1
c----------------------------------------------------------------
c        ray is near the plasma boundary after Runge-Kutta
c        procedure (it is after reflection)
c-----------------------------------------------------------------
c        h=prmt(3)
c         write(*,*)'in rk after output 54 iflagh=',iflagh,'h=',h
c	 write(*,*)'ray is near the plasma boundary after reflection'
      else
        if (iflagh.eq.2) then
c if 2
c-------------------------------------------------------------------
c         Ray is outside  the plasma after the correction.
c         It is necessary to reduce the time step in the Runge-Kutta
c         procedure and to recalculate the array u(i)
c----------------------------------------------------------------------
          t=t-h
          h=h*0.5d0
	  dtstep=prmt(3)/h
cyup	  write(*,*)'ray is outside the plasma after correction'
c------------------------------------------------------------------
          if ((prmt(3)/h).lt.100.d0)then !if 2.1
             do i=1,ndim
                u(i)=startu(i)
             enddo
             goto 10 !another step
          else
             write(*,*)' Initial_step/adjusted_step: prmt(3)/h > 100'
             do i=1,ndim
                u(i)=startu(i)
             enddo
             goto 10 !another step
          endif !end if 2.1
c-------------------------------------------------------------------
        else
c-------------------------------------------------------------------
c     ordinary ray point	  iflagh=3
c     value of u(i) after correction procedure
c-------------------------------------------------------------------
        endif !end if 2
      endif !end if 1
c--------------------------------------------------------------------
      if (iraystop.eq.1) then
         prmt(4)=prmt4_orig ![2025-12-15] reset to original (in case it was adjusted)
         prmt(6)=prmt6_orig ![2025-12-15] reset to original (in case it was adjusted)
         goto 100 !-> finish and stop the ray
      end if     
      if (prmt(5).gt.0) then ! handle to exit integration along ray
         prmt(4)=prmt4_orig ![2025-12-15] reset to original (in case it was adjusted)
         prmt(6)=prmt6_orig ![2025-12-15] reset to original (in case it was adjusted)
         !write(*,*)'prmt(5)>0  goto 100 finish',prmt(5)
         goto 100 !-> finish here
      endif
  30  format(1x,6(1x,e11.4))
c------------------------------------------------------------------
c  the program for calculation time determination
      ntau=ntau+1
c      if (ntau.eq.10) then
c        call gettim(itf1,itf2,itf3,itf4)
c	itd1=itf1-it1
c	itd2=itf2-it2
c	itd3=itf3-it3
c	itd4=itf4-it4
c	tstart=it1*3600.00d0+it2*60.00d0+it3*1.00d0+it4*0.01d0
c	tfinish=itf1*3600.00d0+itf2*60.00d0+itf3*1.00d0+itf4*0.01d0
c	deltim=tfinish-tstart
c	write(*,*)'deltim=',deltim
c	write(*,*)'it1=',it1,'it2=',it2,'itd3=',it3,'itd4=',it4
c	write(*,*)'itf1=',itf1,'itf2=',itf2,'itf3=',itf3,'itf4=',itf4
c	write(*,*)'itd1=',itd1,'itd2=',itd2,'itd3=',itd3,'itd4=',itd4
c        call gettim(time2)
c        time=time2-time1
c        write(*,*) ' 2 - calculation time for 100 step:',time
c      end if
c--------------------------------------------------------
      goto 10 !---> go back for another step along ray.
  100 continue
      return
      end subroutine drkgs2
	 


c global const: NMAX  - max order of system
c               MAXSTEP  - max namber of iteration
c               TINY=1.e-30  - by not division zero
c type of variable:double precision
c*******************************************************************************              
c      subroutine rk_top(prmt,u,deru,ndim,ihlf,fct,outp,aux)
c         1. input parameters                                      !
c                                                                !
c          prmt  - the vector of input and output data:          !
c              prmt(1) - the lower boundary of interval of       !
c                        integration;                            !
c              prmt(2) - the upper boundary ...;                 !
c              prmt(3) - the initial step of integration;        !
c              prmt(4) - the accuracy of solution;               !
c              prmt(5) - if it is not equal to zero the control  !
c                        is transferred to the main program;     !
c                                                                !
c          y     - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  at the intermediate points x;                 !
c                                                                !
c          dery  - the  vector of derivatives of y at point x    !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array.                          !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system. it has not to change the  values  of    !
c                x and y. its formal parameters are:             !
c                x, y, dery.                                     !
c                                                                !
c          outp - output subroutine. it has not to change the    !
c                 values of all  its  formal  parameters.  if    !
c                 prmt(5) is not equal to zero,  the  control    !
c                 will be transferred to the main program.       !
c                 its formal parameters are:                     !
c                 x, y, dery, ihlf, ndim, prmt                   !
c          output date                                           !
c    u(1)=z,u(2)=r,u(3)=phi,u(4)=nz,u(5)=nr,u(6)=m               !                                         !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c****************************************************************************************      
c methode Runge-Kutta. 5 order with change step.  
c scale abs(y)+abs(h*dydx)+TINY    
      subroutine drkgs1(prmt,y0,dery,n,ihlf,func,outp,aux)
c     subroutine rkadap(prmt,y0,dery,n,ihlf,func,outp,aux)
      integer n,ihlf,MAXSTEP,NMAX    
c number  of variable,attribute of out      
      double precision TINY
      double precision y0(*),dery(*),aux(8,*),prmt(*) 
      external func,outp         
      
      parameter(NMAX=10)
      parameter(MAXSTEP=1000000)
      parameter(TINY=1.d-30) 
      double precision y(1:NMAX),yscal(1:NMAX)
      integer i,nstep,iflagh,iraystop,flag
      double precision h,hdid,hnext,x,hmax,us,uz,ur 
      
c      write(*,*)'ok1 rkadap'
      us=0
      uz=y0(1)
      ur=y0(2)
      x=prmt(1)                   
      h=sign(prmt(3),prmt(2)-prmt(1))
      hmax=10
      do 10 i=1,n
        y(i)=y0(i)
10    continue  
      do 100 nstep=1,MAXSTEP
        call func(x,y,dery)
        do 20 i=1,n                 
          yscal(i)=dabs(y(i))+dabs(h*dery(i))+TINY
20      continue          
        if((x+h-prmt(2))*(x+h-prmt(1)) .GE. 0.) h=prmt(2)-x 
        call rk5step(y,y0,n,x,h,prmt(4),yscal,hdid,hnext,func,hmax)
        if((x-prmt(2))*(prmt(2)-prmt(1)) .GE. 0.d0) then 
          flag=1
cyup          write(*,*)'ok end rkadap x,prmt(2)',x,prmt(2)
          return
        end if
        us=us+dsqrt((uz-y0(1))**2+(ur-y0(2))**2)
c        write(*,*)'rk1 uz,y0(1)',uz,y0(1)
c        write(*,*)'rk1 ur,y0(2)',ur,y0(2)
c        write(*,*)'rk1 dus',dsqrt((uz-y0(1))**2+(ur-y0(2))**2) 
        uz=y0(1)
        ur=y0(2)

cSAP120517 for Nicola scattering
        dels_total=dsqrt((uz-y0(1))**2+(ur-y0(2))**2
     &                 +(ur*uphi-y0(2)*y0(3))**2)
        us_total_old=us_total
        us_total=us_total+dels_total
c        write(*,*)'rk1 before output us,t',us,x
        call outp(us,y0,dery,ihlf,n,prmt,iflagh,iraystop,
cSAP100202
c     & i_go_to_previous_output_step)
cSAP120517
     &  i_go_to_previous_output_step,s_total,s_total_old )

c        write(*,*)'rk_new after outp iraystop',iraystop
        if (iraystop.eq.1) then
cyup          write(*,*)'rk_new iraystop=1 goto 200'
          goto 200
        endif

        if(iflagh .EQ. 2) then 
          x=x-hdid
          h=hdid/2
          hmax=hdid
        else  
          h=hnext
          hmax=10
          do 30 i=1,n
            y(i)=y0(i)
30        continue 
        endif
100   continue 
200   continue
c      write(*,*)'too many steps.rk5s'
      flag=0
c      write(*,*)'ok end bad rkadap'
      return
      end         
c***********************************************************************
c method of Runge-Kutta. 5 order.1 step with count next step      
      subroutine rk5step(y0,yout,n,x,htry,eps,yscal,hdid,hnext,func,
     +hmax)
      integer n,NMAX                   
c namber of variable
      double precision x,htry,hdid,hnext,eps,hmax
c start,step,step which did,next step,accurate      
      double precision y0(1:n),yscal(1:n),yout(1:n)    
c start value and result,scale error      
c      external func,rk
      external func
      parameter(NMAX=10)
      integer i
      double precision errmax,h,htmp,xnew,errcon,safety,pgrow,pshrnk
      double precision yerr(1:NMAX),ytmp(1:NMAX)
      parameter(safety=0.999d0,pgrow=-0.2d0,pshrnk=-0.25d0)
      
c      write(*,*)'ok1 rk5step'      
      errcon=(5.d0/safety)**(1/pgrow)         
      h=htry                
10    call rk51(y0,x,h,ytmp,n,func,yerr) 
      errmax=0.
      do 20 i=1,n
        errmax=dmax1(errmax,dabs(yerr(i)/yscal(i)))
20    continue              
      errmax=errmax/eps 
      if(errmax .GT. 1.d0)then  
        htmp=safety*h*(errmax**pshrnk)
        h=dsign(dmax1(dabs(htmp),0.1d0*dabs(h)),h)
        xnew=x+h
        if(xnew .eq. x) pause'step=0.rk5step'
        goto 10   
      else
        if(errmax .GT. errcon) then
          hnext=safety*h*(errmax**pgrow)
        else 
          hnext=(hmax-hdid)/2
          if(hmax .GT. 1.d0) then hnext=5.d0*h 
        end if 
30      hdid=h
        x=x+h 
        do 40 i=1,n
         yout(i)=ytmp(i)
40      continue      
c        write(*,*)'ok end rk5step'
        return
      end if
      end                                                                     
c******************************************************************************************      
c methode Runge-Kutta. 5 order.1 step.with error      
      subroutine rk51(y0,a,h,yout,n,func,yerr)
      integer n,NMAX
c number of variable    
      double precision a,h,y0(1:n),yout(1:n),yerr(1:n)
c start pset,step,start value,solution,error      
      external func  
c left part        
      parameter(NMAX=10)
      double precision q1(1:NMAX),q2(1:NMAX),q3(1:NMAX),q4(1:NMAX),
     +q5(1:NMAX),q6(1:NMAX)
      double precision a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,
     +b53,b54,b61,b62,b63,b64,b65,c1,c3,c4,c6,cc1,cc3,cc4,cc5,cc6         
      integer i
      parameter(b21=.2d0,b31=3.d0/40.d0,b32=9.d0/40.d0,
     +b41=.3d0,b42=-.9d0,b43=1.2d0,
     +b51=-11.d0/54,b52=2.5d0,b53=-70.d0/27.d0,b54=35.d0/27.d0,
     +b61=1631.d0/55296.d0,b62=175.d0/512.d0,b63=575.d0/13824.d0,
     +b64=44275.d0/110592.d0,
     +b65=253.d0/4096.d0,
     +a2=b21,a3=b31+b32,a4=b41+b42+b43,a5=b51+b52+b53+b54,
     +a6=b61+b62+b63+b64+b65,
     +c1=37.d0/378.d0,c3=250.d0/621.d0,c4=125.d0/594.d0,
     +c6=512.d0/1771.d0,
     +cc1=c1-2825.d0/27648.d0,cc3=c3-18575.d0/48384.d0,
     +cc4=c4-13525.d0/55296.d0,cc5=-277.d0/14336.d0,cc6=c6-.25d0) 
      
c c2 and c5 and cc2 = 0 then no count           
c      write(*,*)'ok1 rk51'  
      call func(a,y0,q1)  
      do 20 i=1,n  
        yerr(i)=y0(i)+h*b21*q1(i)
20    continue             
      call func(a+a2*h,yerr,q2) 
      do 30 i=1,n  
        yerr(i)=y0(i)+h*(b31*q1(i)+b32*q2(i))
30    continue      
      call func(a+a3*h,yerr,q3)
      do 40 i=1,n  
        yerr(i)=y0(i)+h*(b41*q1(i)+b42*q2(i)+b43*q3(i))
40    continue      
      call func(a+a4*h,yerr,q4)
      do 50 i=1,n  
        yerr(i)=y0(i)+h*(b51*q1(i)+b52*q2(i)+b53*q3(i)+b54*q4(i))
50    continue      
      call func(a+a5*h,yerr,q5) 
      do 60 i=1,n  
        yerr(i)=y0(i)+h*(b61*q1(i)+b62*q2(i)+b63*q3(i)+b64*q4(i)+
     +b65*q5(i))
60    continue      
      call func(a+a6*h,yerr,q6)
      do 70 i=1,n
        yout(i)=y0(i)+h*(c1*q1(i)+c3*q3(i)+c4*q4(i)+c6*q6(i))       
70    continue
      do 80 i=1,n
        yerr(i)=h*(cc1*q1(i)+cc3*q3(i)+cc4*q4(i)+cc5*q5(i)+cc6*q6(i))       
80    continue
c      write(*,*)'ok end rk51'             
      return
      end







c        ********************** drkgs0    *******************
c         4_th order Runge-Kutta subroutine with constant time 
c         step finds the solution of
c         the system of the ordinary differential equations
c         only for one time step         
c        ****************************************************
c
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c          Runge-Kutta method  of  4th-order                     !
c          with the constant step integration
c          It makes only one step!
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       2. input parameters                                      !
c                                                                !
c          h     - time step                                     ! 
c                                                                !
c          u     - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  after one time stept                         !
c                                                                !
c          deru  - the  vector of derivatives of u at point t	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array.                          !
c                                                                !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system. it has not to change the  values  of    !
c                t and u. its formal parameters are:             !
c                t, u, deru.                                     !
c                                                                !
c       			  		 !
c    u(1)=z,u(2)=r,u(3)=phi,u(4)=nz,u(5)=nr,u(6)=m               !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs0(h,u,deru,ndim,fct,aux)
      !implicit double precision (a-h,o-z)
      implicit none !YuP[2020-01-14]
      real*8 h,u,deru,aux,up, t
      integer ndim,i

      dimension u(*),deru(*),aux(8,*),up(6)

c      write(*,*)'in Runge -Kutta drkgs0'

      
      call fct(t,u,deru)
      do i=1,ndim
         aux(1,i)=h*deru(i)
	 aux(5,i)=u(i)+0.5d0*aux(1,i)
 	 up(i)=aux(5,i)
      enddo
c------------------------------------------------------------------
      call fct(t+0.5d0*h,up,deru)
      do i=1,ndim
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5d0*aux(2,i)
	 up(i)=aux(5,i)
      enddo
c--------------------------------------------------------------------
      call fct(t+0.5d0*h,up,deru)
      do i=1,ndim
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
 	 up(i)=aux(5,i)
      enddo
c--------------------------------------------------------------------
      call fct(t+h,up,deru)
      do i=1,ndim
         aux(4,i)=h*deru(i)
      enddo
c--------------------------------------------------------------------
      do i=1,ndim
         up(i)=u(i)+1.d0/6.d0*(aux(1,i)+2.d0*aux(2,i)+2.d0*aux(3,i)
     *          +aux(4,i))
      enddo

      do i=1,ndim
         u(i)=up(i)
      enddo

      return
      end


      subroutine callc_eps_ham(u,ham)
c-------------------------------------
c     calculates hamiltonian
c     to check rk accuracy
c-----------------------------
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 u(6)    !ray coordinates
c-----output
      real*8 ham
c-----externals
      real*8 b,gamma1,hamilt1

c-----local

      bmod=b(u(1),u(2),u(3))
      gam=gamma1(u(1),u(2),u(3),u(4),u(5),u(6))
      ham=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))

      return
      end

      subroutine calc_hamiltonian(u,hamiltonian)
c-----calculate hamitonian
      implicit none
      include 'param.i'
      include 'one.i'

c-----input
      real*8 u(*)
c-----output
      real*8 hamiltonian   

c-----externals
      real*8 b,gamma1,hamilt1,
     &cn
c-----local
      real*8 ds,dc,ds2,dc2,cnt,d4,d2,d0,hamilt_id1,hamilt_id2,
     &cn2,cn4,ad,bd,cd

      bmod=b(u(1),u(2),u(3))
      gam=gamma1(u(1),u(2),u(3),u(4),u(5),u(6))
      hamiltonian=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))
      
cyup      write(*,*)'hamiltonian',hamiltonian
c-----------------------------------
c      ds=dsin(gam)
c      dc=dcos(gam)
c      ds2=ds*ds
c      dc2=dc*dc 
c      cnt=cn(u(2),u(4),u(5),u(6))
c      cn2=cnt*cnt
c      cn4=cn2*cn2
c      call abc_test(u(1),u(2),u(3),ds2,dc2,ad,bd,cd)
            
c      d4=ad
c      d2=bd
c      d0=cd
c      hamilt_id1=d4*cn4+d2*cn2+d0
c      write(*,*)'hamilt_id1',hamilt_id1
c      hamilt_id2=cn2+(d2-ioxm*dsqrt(d2*d2-4.d0*d4*d0))/
c     &                (2.d0*d4)
c      write(*,*)'d4,d2,d0,cn2,hamilt_id2',
c     &d4,d2,d0,cn2,hamilt_id2

      return
      end
c======================================================================
c======================================================================

c        4_th order  Runge-Kutta method with automatic 
c        time step selection
c
c        **********************  drkgs_auto  ****************
c        *                      -----                       
c        * this Runge-Kutta subroutine finds the solution of the 
c        *      system of ordinary differential equations          
c        ****************************************************
c
c        drkgs_auto: called if(isolv=1 & irkmeth=3)
c
c The set of equations to be solved :
!        Integration in time (default):
!        ray-tracing equations right hand side=
!        dr^/dt_code = -(dD/dN^)/(dD/domega)
!        dN^/dt_code = +(dD/dr^)/(dD/domega)
! sqrt(deru(1)^2 + deru(2)^2 +(R*deru(3))^2) gives v_group (normalized).
c        Note: dt_code here is NOT an actual physical time, but
c        normalized as dt_code= dt_physical[sec]*c/omega  
c                                               
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c          runge-kutta method  of  4th-order                  
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       2. input parameters                                      !
c                         
c          dL_step
c          dN_step  
c     dL_step=1.d-3 ! [m]  max allowed change in cartesian coords.
c     dN_step=1.d-2 ! max allowed change in refraction index.
c     The code will set the time step h = dt_code for integration
c     in such a way that the change |dr| in configuration space
c     is not larger than dL_step, and also
c     the change in refr. index |N| is not larger than dN_step.
c          prmt(3)=initial step of integration (not really needed)  
c              Step h is set below, based on dL_step and dN_step.                                  
c          prmt(6)= distance step along ray [m] for saving data.
c                  Does NOT affect the step of integration 
!                  in this subroutine!
!     Other prmt* values are not needed.
!     This new option (irkmeth=3) is ~4x faster than irkmeth=2 (per ray element)
!     ! The advantage of using both drdt/dL_step and dNdt/dN_step is:
!     ! At plasma edge, for example, when starting as O-mode,
!     ! the value of N is ~1, change of N is also small, dN ~1 
!     ! (going down to 0 at O-X cutoff), while the group velocity 
!     ! is large Vgr~c. In this case the time step is determined by 
!     ! step_inv= drdt/dL_step, so that dt_code= dL_step/drdt.
!     ! In opposite case of slow electrostatic waves, 
!     ! the group velocity can be very low, so that drdt~0.
!     ! On the other hand, the change in N can be very steep,
!     ! dN>10 over a short travel distance.
!     ! In this case the time step is determined by dNdt/dN_step term
!     ! so that dt_code= dN_step/dNdt .
c
c          u()   - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  at the intermediate points u;                 !
c        
c          deru() - the  vector of derivatives at point z,r,phi	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array (working array).  
c                                                                !
c      
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system. It has not to change the  values  of    !
c                u().     Its formal parameters are: u, deru.                                     !
c                                                                !
c          outp - output subroutine. it should not change the
c                 values of all its  formal  parameters.  
c                 output data:                                   !
c          u(1)=z, u(2)=r, u(3)=phi, u(4)=Nz, u(5)=Nr, u(6)=cm     !                                         !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs_auto(prmt,u,deru,ndim,ihlf,fct,outp,
     +  aux,dL_step,dN_step,nbulk,nharm_refined_step)  
       ![2025-12-16] added nbulk,nharm_refined_step
      implicit none !integer (i-n), real*8 (a-h,o-z)
      integer ndim,ihlf !Input
      real*8 dL_step,dN_step !input
      dimension u(6),deru(6),aux(8,6),prmt(9),up(6),uu(6),startu(6)
      integer i_output
      dimension uplus(6)
      double precision prmt,u,deru,aux,t,tend,h,up,hh,uu,startu,
     1 dd,uplus
      double precision us, dels1,dels2, uz,ur,uphi
      integer iflag,iflagh, i, iraystop, i_go_to_previous_output_step !local
      real*8 tdum,dt_code, us_total, us_total_old, dndt,drdt,step_inv
      
      real*8 dL_step_orig, dL_step_adj !original and adjusted
      integer nbulk,nharm_refined_step !INPUT ![2025-12-16]added
      integer ic_near,nh ![2025-12-16]
      real*8 yi
      real*8 y !external
      
      external fct !=rside1() !YuP: Does not actually dep. on t
      external outp != outpt()
      
      dL_step_orig= dL_step ![2025-12-16] Save original value

      h=prmt(3) !step of integration; not really needed here:
      !h is set below, based on dL_step and dN_step.
      !(But reserve for future, maybe will be used?)
      prmt(5)=0.d0 ! YuP: if =1.d0, it will stop the ray after 1 step
      iflag=0  ! control of the RK method accuracy
      iflagh=3 !drkgs_auto: initialize (control of the plasma edge intersection)
      us=0.d0  !cumulative distance along ray (initial value)
      dels2= 0.d0 !  to initialize
      tdum=0.d0 !YuP: rside1 does not actually dep. on t, so we use tdum
      
c---------------------------------------------------------------
 10   continue ! handle for time-advancing (loop)

      ![2025-12-15] Adjust integration step near resonance.
      ic_near=0 ! To be changed if ray element is near resonance.
      !Restore to original (from genray namelist):
      dL_step_adj= dL_step_orig
      if(nharm_refined_step>0)then
      uz=u(1)
      ur=u(2)
      uphi=u(3)
      do i=1,nbulk !scan species
        !If at least one species has nearly res. - change ic_near to 1
        if(ic_near.eq.0)then
          yi= dabs(y(uz,ur,uphi,i)) != |omega_ci/omega|
          do nh=1,nharm_refined_step !scan harmonics
            if( dabs(1.d0-nh*yi) .le. 0.1d0 ) then
              ic_near=1 !flag for proximity of resonance
            endif
          enddo !nh
        endif
      enddo
      if(ic_near.eq.1)then !Near resonance. Reduce steps:
        dL_step_adj= dL_step_orig*0.1 !Reduce by 10x (hardwired, for now)
      endif 
      endif !(nharm_refined_step>0)
      ![2025-12-15] done

      ! Save values from previous time step:
      uz=u(1) !=Z
      ur=u(2) !=R
      uphi=u(3) !=phi
      do i=1,ndim
        uu(i)=u(i)
        startu(i)=u(i)
      enddo
      call fct(tdum,u,deru) ! rhs of ODE: get dr/dt and dN/dt
      !Note: deru(1)= dD/dN_z, deru(2)= dD/dN_r, deru(3)= dD/dCM
      !      deru(4)=-dD/dZ,   deru(5)=-dD/dR,   deru(6)=-dD/dphi
      
      !Set the time step.
      ! This is the speed of ray propagation (normalized Vgroup):
      drdt= sqrt(deru(1)*deru(1)+deru(2)*deru(2)+(ur*deru(3))**2)
      ! This is the normalized "speed" of change of refr. index:
      dNdt= sqrt(deru(4)*deru(4)+deru(5)*deru(5)+(deru(6)/ur)**2)
      
      ! consider the inverse time steps 
      ! 1/dt ~ (dr/dt)/dL_step
      ! and 
      ! 1/dt ~ (dN/dt)/dN_step
      ! Select the larger of the two - it will correspond to the smaller 
      ! value of dt
      step_inv= max( drdt/dL_step_adj , dNdt/dN_step )
      ! Why we need to consider the inverse time step:
      ! because one of these speeds, dr/dt or dN/dt, can be exactly 0
      ! or at least very close to 0.

      if(step_inv.gt.0.d0)then
        dt_code= 1.d0/step_inv !Normalized units: dt_code= dt[sec]*c/omega
        h= dt_code
      else ! both drdt and dNdt are 0 
           !(could be very small, below rounding error)
        write(*,*)'***** In Runge-Kutta subroutine drkgs_auto *****'
        write(*,*)'***** time step is too small.'
        write(*,*)'***** dr/dt, dN/dt=',drdt,dNdt
        write(*,*)'***** dL_step=',dL_step_adj
        write(*,*)'***** dN_step=',dN_step
        write(*,*)'  iraystop->1'
        iraystop=1 ! step_inv is too small (=0)
        !iray_status_one_ray=13 ! step_inv is too small !in subroutine drkgs_auto
        goto 100 ! exit
      endif
      ! The advantage of using both drdt/dL_step and dNdt/dN_step is:
      ! At plasma edge, for example, when starting as O-mode,
      ! the value of N is ~1, change of N is also small, dN ~1 
      ! (going down to 0 at O-X cutoff), while the group velocity 
      ! is large Vgr~c. In this case the time step is determined by 
      ! step_inv= drdt/dL_step, so that dt_code= dL_step/drdt.
      ! In opposite case of slow electrostatic waves, 
      ! the group velocity can be very low, so that drdt~0.
      ! On the other hand, the change in N can be very steep,
      ! dN>10 over a short travel distance.
      ! In this case the time step is determined by dNdt/dN_step term
      ! so that dt_code= dN_step/dNdt .
      
      do i=1,ndim	       
        aux(1,i)=h*deru(i)
        aux(5,i)=u(i)+0.5*aux(1,i)
        up(i)=aux(5,i)
      enddo

      call fct(tdum,up,deru)
      do i=1,ndim
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5*aux(2,i)
         up(i)=aux(5,i)
      enddo

      call fct(tdum,up,deru)
      do i=1,ndim
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
         up(i)=aux(5,i)
      enddo
      
      call fct(tdum,up,deru)
      do i=1,ndim
        aux(4,i)=h*deru(i)
      enddo
      
      do i=1,ndim
        up(i)=u(i)+1.d0/6.d0*
     1        (aux(1,i)+2.d0*aux(2,i)+2.d0*aux(3,i)+aux(4,i))
      enddo
      ! Done: new u(i) after one step is obtained.
      do i=1,ndim
        u(i)=up(i)
      enddo
c--------------------------------------------------------------------
c     This section can be used for hamiltonian (dispersion D=0) control,
c     but we skip it for now.
c      zu=  u(1)
c      ru=  u(2)
c      phiu=u(3)
c      cnzu=   u(4)    !=Nz refr.index
c      cnru=   u(5)    !=Nr
c      cmu=    u(6)
c      cnphiu= u(6)/ru !=Nphi
c      cn2u= cnzu*cnzu + cnru*cnru + cnphiu*cnphiu
c      bmod=b(zu,ru,phiu) !-> get b and derivs of b
c      gam=  gamma1(zu,ru,phiu,cnzu,cnru,cmu) ! uses bz,br,bphi,bmod
c      ham= hamilt1(zu,ru,phiu,cnzu,cnru,cmu) ! uses gam
c--------------------------------------------------------------------      

c-----------------------------------------------------------
c     total distance traveled in one step
      dels2= (uz-u(1))**2 + (ur-u(2))**2 + (ur*uphi-u(2)*u(3))**2
      dels1= dsqrt(dels2)
      us= us+dels1 ! accumulate the distance with each step. 
      iflagh=3 !drkgs_auto/before outp : initialize
      !========================================================
      call outp(us,u,deru,ihlf,ndim,prmt,iflagh,iraystop,
     & i_go_to_previous_output_step,us_total,us_total_old )
      !outpt, among other control, increments the time step
      !nstep_rk ; if nstep_rk is bigger than maxsteps_rk
      !it will stop the ray.
      !Among other staff, the outpt checks: if(us>prmt(7)),
      !which is initially prmt(7)=prmt(6), 
      !then save ray data and also increase prmt(7): add prmt(6), 
      !so the new value of prmt(7) will be used 
      !during next call of outpt.
      !This way, the data is saved each prmt(6) steps.
      !========================================================
      !Ray can also be stopped if Vgroup/c >2, or 
      ! D/(N|gradD|) > toll_hamilt
      ! Also, the value of iflagh can be changed 3->1
      ! if(iflref.eq.1)
      if(iraystop.eq.1) then
         goto 100 ! finish and stop the ray
      endif
      
      goto 10  !drkgs_auto: another time step
      
  100 continue !drkgs_auto: exit handle
      return
      end subroutine drkgs_auto
c======================================================================
c======================================================================
