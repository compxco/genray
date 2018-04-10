c        ********************** outputb1*********************
c        *                      -----                       *
c        *   output is used id rkb as			    *
c        *   output subroutine. it has not to change the    *
c        *   values of all  its  formal  parameters.  if    *
c        *   prmt(5) is not equal to zero,  the  control    *
c        *   will be transferred to the main program.       *
c        *       its formal parameters are:                 *
c        *         x, y, dery, ihlf, ndim, prmt             *
c        *                                                  *
c        ****************************************************
c        !   this code printes: values of u(ndim);	    !
c        !   components of the electric field; 		    !
c        !   it writes array u(6) in file( the name of      !
c        !   file is given as  the first parameter in	    !
c        !   genray.in file)                   	    !
c        !   eps- value of hamiltonian. 		    !
c        !   it controls the conditions: if the ray point   !
c        !   is inside the plasma or not .		    !
c        !   it controls the conditions: if the refractive  !
c        !   index less the cnmax or not .		    !
c        !   it writes the output date in mnemonic.txt file !
c        !   it calculates the sixth ray variable	    !
c------------------------------------------------------------------
c          output data						  *
c    u(1)=z,u(2)=r,u(3)=phi,u(4)=nz,u(5)=nr,u(6)=m,iraystop       *                                         !
c----------------------------------------------------------------------
c         output parameter:
c            if iraystop=1  then end of the j_ray calculation
c            if iraystop=0  then continuation of the j_ray calculations
c           it prepares and writes  parameters for 3d code
c----------------------------------------------------------------------
c---------------------------------------------------------------------*
c           this program uses the following functions and	      *
c           subroutines  b,gamma1,hamilt1,bound,sdr,prep3d,write3d    *
c           Is used for isolv=2 case as a rkb1 external subroutine.   *
c---------------------------------------------------------------------*
      subroutine outptb1(t,u,deru,ihlf,ndim,prmt,i0,iraystop)
      implicit double precision (a-h,o-z)
      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'write.i'
      include 'cefield.i'
      dimension u(*),deru(*),prmt(*),pu(6)
c	 write(*,*)'in outptb1'
	  iraystop=0
          z1=u(1)
          r1=u(2)
          phi1=u(3)
          cnz1=u(4)
          cnr1=u(5)
          cm1=u(6)
         if(t.lt.prmt(7))goto 20
          u(i0)=sdr1(u,deru,i0,prmt(4))
c          write(*,*)'tau=',t,' y=',y(z1,r1,phi1,1)
c--------------------------------------------------------------------
c        if ray is close to resonance point then stop ray calculations
          cnmode=cn(r1,cnz1,cnr1,cm1)
          cnmax=500.d0
	    if(cnmode.gt.cnmax) then
c	      write(*,*)'***************nmode.gt.cnmax***********'
              call prep3d(t,u,deru,iraystop)

!nme              if(iout3d.eq.'enable') then  [iout3d is DEFUNCT]
c----------------writing data to output 3d.dat file
!nme	         call write3d
!nme              endif

	      iraystop=1
	      return
	    end if
c---------------------------------------------------------
c denormalization of ray point coordinates
c---------------------------------------------------------
	  uout1=r0x*u(1)
	  uout2=r0x*u(2)
	  uout3=u(3)
	  uout4=u(4)
	  uout5=u(5)
	  uout6=u(6)
c-----------------------------------------------------------
c          write(*,30)uout2,uout1,uout3,uout4,uout5,uout6
c          write(*,*)uout2,uout1,uout3,uout4,uout5,uout6
        if(myrank.eq.0)then
           write(i1_,110)uout2,uout1,uout3,uout4,uout5,uout6
        endif ! myrank=0

c         do 200 i=1,ndim
c           write(2,110)u(i)
c 200     continue
          prmt(7)=prmt(7)+prmt(6)
 110     format(3x,6(' ',e13.6))
c-------------------------------------------------------------------
         bmod=b(u(1),u(2),u(3))
         gam=gamma1(u(1),u(2),u(3),u(4),u(5),u(6))
         ds2=dsin(gam)**2
         dc2=dcos(gam)**2
         eps=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))
c------------------------------------------------------------------
c         write(*,568)eps
 568     format(3x,'eps=',d13.6)
         call prep3d(t,u,deru,iraystop)
         
	 if (iraystop.eq.1) then

!nme           if(iout3d.eq.'enable') then
c-------------writing data to output 3d.dat file  [iout3d is DEFUNCT]
!nme	      call write3d
!nme           endif

	   return
	 end if
	 if (nrayelt.eq.nrelt) then
c	   write(*,*)'***************nrayelt=nrelt***********'

!nme           if(iout3d.eq.'enable') then
c-------------writing data to output 3d.dat file  [iout3d is DEFUNCT]
!nme	      call write3d
!nme           endif

	   iraystop=1
	   return
	 end if
c***************************************************
  20     continue
  30      format(1x,6(1x,e11.4))
c------------------------------------------------
c       call bound(u(1),u(2),u(3),u(4),u(5),iflref,cnzref,cnrref,ibound)
      call bound(u(1),u(2),u(3),u(4),u(5),u(6),iflref,
     &z_ref,r_ref,phi_ref,cnzref,cnrref,cmref,
     &ibound,deru(1),deru(2))
	 if (irefl.ge.ireflm) then
          if(myrank.eq.0)then
            write(*,*)'irelf.ge.ireflm'
          endif ! myrank=0
          call prep3d(t,u,deru,iraystop)

!nme           if(iout3d.eq.'enable') then
c------------writing data to output 3d.dat file  [iout3d is DEFUNCT]
!nme	     call write3d
!nme           endif

 	   iraystop=1
	 end if
  
         if (iflref.eq.1) then
           u(1)= z_ref
           u(2)= r_ref
           u(3)= phi_ref
           u(4)=cnzref
           u(5)=cnrref
           u(6)=cmref
         endif

         if(ibound.gt.0.)then
	  prmt(5)=1.
c	  write(*,*)'bound gt 0'
	 end if

       return
       end
