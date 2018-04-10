c     ******************************besk2as***********************
c     *  besk2as  calculates the Macdonalds function		 *
c     *  of the second order K_2(x) multiplied by EXP(x)	 *
c     *  fox x>50 by using assimptotical representation          *
c     *	 result: bk=exp(x)*K_2(x)				 *
c     ***********************************************************
c     * input parameter real X					 *
c     -----------------------------------------------------------
      subroutine besk2as(x,bk)
      implicit double precision(a-h,o-z)
c      write(*,*)'beg.of besk2as'
      pi=4.d0*datan(1.d0)
      m=20
      pp=dsqrt(pi/(2.d0*x))
       s=0.d0
       coef=1.d0
       do 10 j=1,m
        if (j.eq.1) then
	  coef=1.d0
	else
	  ppp=2.d0*(j-1)-1.d0
	  ppp=ppp**2
          ppp=16.0d0-ppp
	  ppp=ppp/(j-1)
          coef=coef*ppp
	end if
          ppp1=2.d0*(j-1)
	  ppp=2.d0**ppp1
        if (j.eq.1) then
	  ppp=ppp*1.d0
	else
	  ppp=ppp*(j-1)
	end if
          coef1=coef/ppp
	  s=s+coef1/(2.0d0*x)**(j-1)
10      continue
	bk=pp*s
       end
