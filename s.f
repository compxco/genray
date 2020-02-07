c        **********************   s   ************************
c        *                        -                          *
c        * this subroutine calculates  six  different  sums  *
c        * that are used by rside for cold plasma dispersion *
c        *****************************************************
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
c								    
c     s1 = 1 - sum_(i.ne.ib) { x_i/[1-(y_i)**2] }, sum is present only if ib>1    
c              (sum_i is present only if ib>1)								    
c     s2 = 1 - sum_(i.ne.ib) { x_i/[1-y_i] },      sum is present only if ib>1 
c              sum_(i.ne.ib) means the sum over non-resonant ion species   
c  								    
c     s3 = 1 - sum_(i) { x_i/[1+y_i] },      		            
c								    
c     s4 = 1 - sum_(i) { x_i } - x_1 ,      		            
c								    
c     s6 = 1 - sum_(i) {x_i/[1-y_i]} ,        sum_i is present only if ib=1
c								    
c     s7 = 1 - sum_(i) {x_i/[1-(y_i)**2]} ,   sum_i is present only if ib=1
c								    
c     sum_() is over i=2,nbulk (ion species); 
c     ib can be one of 1,nbulk;   nbulk - number of bulk particles
c-------------------------------------------------------------------
      subroutine s(z,r,phi,s1,s2,s3,s4,s6,s7)
      !implicit double precision (a-h,o-z)
      implicit none !YuP[2020-01-14]
      real*8 z,r, phi,s1,s2,s3,s4,s6,s7
      real*8 ds1,ds2,ds3,ds4,ds6,ds7,xi,x,yi,y,pps
      integer i
      include 'param.i'
      include 'one.i'
	ds1=0.d0
        ds2=0.d0
        ds3=0.d0
        ds4=0.d0
	ds6=0.d0
	ds7=0.d0

      if (nbulk.gt.1) then
c if 1 begin
        do 10 i=2,nbulk ! sum is over i=2,nbulk (ion species)

          xi=x(z,r,phi,i)
          yi=y(z,r,phi,i)
          ds3=ds3+xi/(1.d0+yi)
          ds4=ds4+xi
         
c if 2 begin
            if (ib.eq.1) then
	      ds6=ds6+xi/(1.d0-yi)
	      ds7=ds7+xi/(1.d0-yi*yi)             
            end if
c if 2 end

          if ((i.eq.ib).or.(ib.eq.1)) goto 10
            ds1=ds1+xi/(1.d0-yi*yi)
            ds2=ds2+xi/(1.d0-yi)             
 10     continue

      end if
c if 1 end
c finish if ef nbulk.gt.1
        s1=1.d0-ds1
        s2=1.d0-ds2
        s3=1.d0-ds3      
        pps=x(z,r,phi,1) ! Xe for electrons is added
        s4=1.d0-ds4-pps
	s6=1.d0-ds6
	s7=1.d0-ds7
       
      return
      end
