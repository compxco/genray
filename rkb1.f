c        **********************  rkb1 ************************
c        *                      -----                        *
c        * this subroutine finds the solution of the system  *
c        * of ordinary differential equations    for isolv=2 *
c        *****************************************************
c
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c             runge-kutta method  of  4th-order  with  constant  !
c          integration step.                                     !
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
c          dery  - the input vector of weight coefficients of    !
c                  error. the sum of its components  must  be    !
c                  equal to 1 at the  beginning.  later  dery    !
c                  becomes the  vector of derivatives of y at 	 !
c                  point x;                                      !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          ihlf  - the number of divisions of initial step.      !
c                  if ihlf=11 the control is transferred to      !
c                  the main program.  the  message  ihlf=12      !
c                  is obtained when prmt(3)=0. and  ihlf=13      !
c                  if sign(prmt(3)).ne.sign(prmt(2)*prmt(1)).    !
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
c                                                                !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                 solution of 5 ray-tracing equations            *
c		  for u(i) i=1,ndim, i.ne.i0			 *
c		  and u(i0) is a  dispersion relation solution   *
c                 with given accuracy    			 *
c----------------------------------------------------------------!
      subroutine rkb1(prmt,u,deru,ndim,ihlf,fct,outp,aux)
      dimension u(*),deru(*),aux(8,*),prmt(*),up(6)
      double precision prmt,u,deru,aux,t,tend,h,up,epsten
c      integer*2 itb1,itb2,itb3,itb4,itf1,itf2,itf3,itf4,
c     1itd1,itd2,itd3,itd4
      write(*,*)'in rkb1 '
c----------------------------------------------------------------
c  the program for timing determination
c----------------------------------------------------------------
c for sun but not working
c       call gettim(time1)
c---------------------------------------------------------------
c for ibm
c      call gettim(it1,it2,it3,it4)
c      tstart=it1*3600.00d0+it2*60.00d0+it3*1.00d0+it4*0.01d0
c----------------------------------------------------------------
      t=prmt(1)
      tend=prmt(2)
      h=prmt(3)
      epsten=prmt(4)
      write(*,*)'accuracy of the solution hamiltonian equation for u(i0)
     1hamiltonian=epsten, epsten=',epsten
      prmt(5)=0.d0
c---------------------------------------------------------------
c    ntau is the number of runge-kutta time step.it may be deleted
c    if you do not want to find the calculation time
      ntau=0
c------------------------------------------------------------------
 10   continue
      write(*,*)'before number u',u(1),u(2),u(3),u(4),u(5),u(6)
      write(*,*)'deru',deru(1),deru(2),deru(3),deru(4),deru(5),deru(6)
      i0=number1(u,deru)
      write(*,*)'after number i0',i0
      call fct(t,u,deru,i0,epsten)
      write(*,*)'after fct1'
      do 1 i=1,ndim
         if (i.eq.i0) go to 1
         aux(1,i)=h*deru(i)
	 aux(5,i)=u(i)+0.5*aux(1,i)
 	 up(i)=aux(5,i)
 1    continue
	 up(i0)=u(i0)
c        aux(1,i)=k1(i)
c        ayx(5,i)=y0(i)+0.5*k1(i)
c------------------------------------------------------------------
      write(*,*)' rkb1 up',up
      write(*,*)'rkb1 before fct2'
      call fct(t+0.5*h,up,deru,i0,epsten)
      write(*,*)'after fct2'
      do 2 i=1,ndim
         if (i.eq.i0) go to 2
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5*aux(2,i)
 	 up(i)=aux(5,i)
 2     continue
	 up(i0)=u(i0)
c        aux(2,i)=k2(i)
c        ayx(5,i)=y0(i)+0.5*k2(i)
c--------------------------------------------------------------------
      call fct(t+0.5*h,up,deru,i0,epsten)
      write(*,*)'after fct3'
      do 3 i=1,ndim
         if (i.eq.i0) go to 3
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
 	 up(i)=aux(5,i)
 3     continue
	 up(i0)=u(i0)
c        aux(3,i)=k3(i)
c        ayx(5,i)=y0(i)+k3(i)
c--------------------------------------------------------------------
c        control that point up(i) is inside the plasma
         write(*,*)'before boundc'
	 call boundc(up(1),up(2),iboundc)
         write(*,*)'after boundc'
	 if (iboundc.eq.1) then
c          ray point is outside the plasma,
c          we change time step h
	   iflagh=-1
	   h=0.5d0*h
	   write(*,*)'in rk3 iboundc=1,h',h
	   go to 10
	 end if
      write(*,*)'before fct4'
      call fct(t+h,up,deru,i0,epsten)
      write(*,*)'after fct4'
      do 4 i=1,ndim
         if (i.eq.i0) go to 4
         aux(4,i)=h*deru(i)
 4     continue
c        aux(4,i)=phi3(i)
c--------------------------------------------------------------------
      do 5 i=1,ndim
         if (i.eq.i0) go to 5
         u(i)=u(i)+1./6.*(aux(1,i)+2.*aux(2,i)+2.*aux(3,i)+aux(4,i))
 5     continue
c--------------------------------------------------------------------
      t=t+h
c-----------------------------------------------------------
      write(*,*)'before outp'
      call outp(t,u,deru,ihlf,ndim,prmt,i0,iraystop)
      write(*,*)'after outp'
      if (iraystop.eq.1) goto 100
      if (prmt(5).gt.0) goto 100
  30      format(1x,6(1x,e11.4))
c----------------------------------------------------------
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
