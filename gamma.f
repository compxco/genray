c        ********************** gamma1 ***********************
c        *                      -----                        *
c        * this function calculates the angle  between  the  *
c        * wave vector and the magnetic field and th         *
c	 * derivatives from cos(gamma1)by z,r,phi,cnz,cnr,cm *
c        *****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of  the  point  where  the  angle  is !
c                 calculated.      		         	    !
c								    !
c      cnz, cnr, cm - n_z, n_r, r*n_phi - components of  wave  ref- !
c                     ractive index.                                !
c-------------------------------------------------------------------
c     this program uses the function cn 			    !
c-------------------------------------------------------------------
	double precision
     1function gamma1(z,r,phi,cnz,cnr,cm)
      !implicit double precision (a-h,o-z)
      implicit none !YuP[2020-01-14]
      real*8 z,r,phi,cnz,cnr,cm  !YuP[2020-01-14]
      real*8 cn,arg,cnt,gg,dc  !YuP[2020-01-14]
      include 'param.i'
      include 'one.i'

	 cnt=cn(r,cnz,cnr,cm)
	gg=cnz*bz+cnr*br+cm*bphi/r

        if(((cnt.eq.0.d0).or.(bmod.eq.0.d0)).or.(r.eq.0.d0))then
          write(*,*)'gamma1,z,r,phi',z,r,phi
          write(*,*)'cnz,cnr,cm,cnt,bmod,gg',cnz,cnr,cm,cnt,bmod,gg
        endif

c        write(*,*)'gamma1 z,r,phi',gamma1 z,r,phi
c        write(*,*)'gamma1 cnz,cnr,cm',cnz,cnr,cm
      

        arg=gg/(cnt*bmod)  
c        write(*,*)'gg,cnt,bmod,arg',gg,cnt,bmod,arg

        if(dabs(arg).ge.1.d0) then
c        if(dabs(gg/(cnt*bmod)).gt.1.d0) then
c          write(*,*)'z,cnz,cnr,cm,cnt',z,cnz,cnr,cm,cnt
c          write(*,*)'bz,br,bphi,bmod',bz,br,bphi,bmod
          WRITE(*,*)'GAMMA: cnt,r,gg,cnz,cnr,cm',cnt,r,gg,cnz,cnr,cm
	  if  (abs(arg).gt.1.d+100) arg=sign(1.d+100,arg)
	  WRITE(*,*) 'GAMMA: ARG'
          if(arg.ge.1.d0)then 
            arg=1.d0
            dc=1.d0
            gamma1=0.d0
          else
            arg=-1.d0
            dc=-1.d0
            pi=4.d0*datan(1.d0)
            gamma1=pi
          endif
        else
          gamma1=dacos(arg)
          dc=dcos(gamma1)
        endif

c	gamma1=dacos(gg/(cnt*bmod))                
c       dc=dcos(gamma1)
       
	dc2dz=2.d0*dc*((cnz*dbzdz+cnr*dbrdz+cm*dbphdz/r)/(cnt*bmod)-
     1               dc*dbmdz/bmod)
	dc2dr=2.d0*dc*((cnz*dbzdr+cnr*dbrdr+cm*dbphdr/r-cm*bphi/(r*r))/
     1                (cnt*bmod)-
     2                dc*dbmdr/bmod+dc*cm*cm/(r**3*cnt**2))
	dc2dph=2.d0*dc*((cnz*dbzdph+cnr*dbrdph+cm*dbpdph/r)/(cnt*bmod)-
     1               dc*dbmdph/bmod)
        dc2dnz=2.d0*dc*(bz/(cnt*bmod)-dc*cnz/(cnt*cnt))
        dc2dnr=2.d0*dc*(br/(cnt*bmod)-dc*cnr/(cnt*cnt))
        dc2dm=2.d0*dc*(bphi/(r*cnt*bmod)-dc*cm/(cnt*cnt*r*r))


      return
      end
