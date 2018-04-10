c        **********************   se   **********************
c        *                        -                         *
c        * this subroutine calculates eight different  sums *
c        * that are used by efield in cold plasma           *
c        ****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of  the  point  where  the  angle  is !
c                 calculated.      		         	    !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c								    !
c        output parameters                                          !
c								    !
c     s1 = 1 - sum_(i.ne.ib) { x_i/[1-(y_i)**2] }, at    ib.ne.1    !
c								    !
c     s2 = 1 - sum_(i.ne.ib) { x_i/[1-y_i] },      at    ib.ne.1    !
c  								    !
c     s3 = 1 - sum_(i) { x_i/[1+y_i] },      		            !
c								    !
c     s4 = 1 - sum_(i) { x_i } - x_1 ,      		            !
c								    !
c     s5 = sum_(i) { x_i * y_i/[1-y_i*y_i] } ,     at  	 ib.eq.1    !
c		  						    !
c     s6 = 1 - sum_(i) {x_i/[1-y_i] ,	           at    ib.eq.1    !
c								    !
c     s7 = 1 - sum_(i) {x_i/[1-(y_i)**2] ,	   at    ib.eq.1    !
c								    !
c     s8 = sum_(i) (i.ne.ib) { x_i * y_i/[1-y_i*y_i] } ,            !
c								    !
c     i=2,nbulk; ib=1,nbulk; nbulk - number of bulk particles       !
c-------------------------------------------------------------------
      subroutine se(z,r,phi,s1,s2,s3,s4,s5,s6,s7,s8)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'one.i'
	ds1=0.d0
        ds2=0.d0
        ds3=0.d0
        ds4=0.d0
	ds5=0.d0
	ds6=0.d0
	ds7=0.d0
	ds8=0.d0

      if (nbulk.gt.1) then
c if 1 begin
        do 10 i=2,nbulk

          xi=x(z,r,phi,i)
          yi=y(z,r,phi,i)
          ds3=ds3+xi/(1.+yi)
          ds4=ds4+xi
c if 2 begin
            if (ib.eq.1) then
	      ds5=ds5+xi*yi/(1.-yi*yi)
	      ds6=ds6+xi/(1.-yi)
	      ds7=ds7+xi/(1.-yi*yi)
            end if
c if 2 end

          if ((i.eq.ib).or.(ib.eq.1)) goto 10
            ds1=ds1+xi/(1.-yi*yi)
            ds2=ds2+xi/(1.-yi)
	    ds8=ds8+xi*yi/(1.-yi*yi)
 10     continue

      end if
c if 1 end
c finish if ef nbulk.gt.1
        s1=1.-ds1
        s2=1.-ds2
        s3=1.-ds3
        pps=x(z,r,phi,1)
        s4=1.-ds4-pps
        s5=ds5
	s6=1.-ds6
	s7=1.-ds7
	s8=ds8
      return
      end
