c        ********************** efield ***********************
c        * this subroutine calculates the components of  the *
c        * wave electric field using dielectric tensor for   *
c        * cold plasma  or for Appleton-Hartry case           *
c        *********************w********************************
c--------------------------------------------------------------------
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of the point where the electric field !
c                 is calculated.      		         	    !
c								    !
c      cnz, cnr, cm - n_z, n_r, r*n_phi - components of  wave  ref- !
c                     ractive index.                                !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c        output parameters					    !
c     ex,ey,ez							    !
c     eplus, eminus, epar - wave polarization 			    !
c     cex,cey,cez in common/cefield/
c--------------------------------------------------------------------
      subroutine efield (z,r,phi,cnpar,cnper,cm,ex,ey,ez,
     1 eplus,eminus,epar)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include'one.i'
      include 'cefield.i'
      double complex ci
      double complex cepsd(3,3),cp
c      write(*,*)'in efield cnpar,cnper',cnpar,cnper
       ci=dcmplx(0.0d0,1.0d0)
       dcn=dsqrt(cnpar**2+cnper**2)
       ds=cnper/dcn
       dc=cnpar/dcn
       ds2=ds*ds
       dc2=dc*dc
c---------------------------------------------------------
c
c     Appleton-Hartry  dispersion relation
c
        if (id.eq.3) then
          cn2=dcn*dcn
          xi=x(z,r,phi,1)
          yi=y(z,r,phi,1)
           deta=1.-yi*yi
            ex=(cn2*ds2-1.+xi)/(cn2*ds*dc)
c-----------------------------------------------------------
c it us necessary to control this line g>-g (for cold plasma)
            ey=-(cn2*dc*ds*deta+ex*(deta*(1.-cn2*dc2)-xi))/(xi*yi)
c-----------------------------------------------------------
            ez=1.
            dmode=1./dsqrt(ex*ex+ey*ey+ez*ez)
              ex=ex*dmode
              ey=ey*dmode
              ez=ez*dmode
         eplus=(ex-ey)
         eminus=(ex+ey)
         eplus1=(ex-ey)/dsqrt(2.0d00)
         eminus1=(ex+ey)/dsqrt(2.0d00)
         epar=ez
	 go to 20
	 end if
c        end if Apleton-Hartry
c--------------------------------------------------------------
c     cold plasma
c
200   continue
      call se(z,r,phi,s1,s2,s3,s4,s5,s6,s7,s8)
      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
        p2=xe/(1.d0+ye)
        p3=p2*ye
        pc=1.d0+dc2

	if (ib.eq.1) then
        del=1.-ye
          epstg1=s7
          epstg0=-p2
          epsg1=-s5
          epsg0=+p3
	go to 10
	else
          xb=x(z,r,phi,ib)
          yb=y(z,r,phi,ib)
          del=1.d0-yb
          epstg1=s1-xe/(1.d0-ye*ye)
          epstg0=-xb/(1.d0+yb)
          epsg1=-s8+xe*ye/(1.d0-ye*ye)
          epsg0=-xb*yb/(1.d0+yb)
	end if
 10    continue
          epslon=s4
          cntg=dcn*ds
          cnlon=dcn*dc
c--------------------------------------------------------
	 cn2=dcn*dcn
	 cntg2=cntg*cntg
         cnlon2=cnlon*cnlon
c-------------------------------------------------------
       cepsd(1,1)=dcmplx(del*epstg1+epstg0,0d0)
       cepsd(1,2)=dcmplx(0.d0,del*epsg1+epsg0)
       cepsd(1,3)=dcmplx(0.d0,0.d0)
       cepsd(2,1)=-cepsd(1,2)
       cepsd(2,2)=cepsd(1,1)
       cepsd(2,3)=cepsd(1,3)
       cepsd(3,1)=cepsd(1,3)
       cepsd(3,2)=cepsd(1,3)
       cepsd(3,3)=dcmplx(del*epslon,0.d0)

c            write(*,*)'cnpar,cnper',cnpar,cnper,'cntg,cnlon',cntg,cnlon
       cp=dcmplx(1.d0/del,0.d0)
ccc       do i=1,3
ccc         do j=1,3
ccc	    cepsd(i,j)=cepsd(i,j)*cp
ccc	    write(*,*)'i,j,cepsd',i,j,cepsd(i,j)
ccc	 enddo
ccc      enddo
c-----------------------------------------------------------
       dlam33=cntg2-epslon
       dlam332=dlam33*dlam33
       dlam13=-cnlon*cntg
       dlam221=cn2-epstg1
       dlam220=-epstg0
       dlam121=-epsg1
       dlam120=-epsg0
       dlam22=del*dlam221+dlam220
       dlam12=del*dlam121+dlam120
       ex=dlam33*dlam22
c
c  ey is image,it must be multiplide by i
c
       ey=dlam33*dlam12
       ez=-dlam13*dlam22
       p=1.d00/dsqrt(ex*ex+ey*ey+ez*ez)
       ex=ex*p
       ey=ey*p
       ez=ez*p
       cex=dcmplx(ex,0.d0)
       cey=dcmplx(0.d0,ey)
       cez=dcmplx(ez,0.d0)
c       write(*,*)'cold cex=',cex,'cey=',cey,'cez=',cez
c------------------------------------------------------------
c  end cold plasma
       eplus=(ex-ey)/dsqrt(2.d0)
       eminus=(ex+ey)/dsqrt(2.d0)
       epar=ez
       dmode=dsqrt(eplus*eplus+eminus*eminus+epar*epar)
 20   continue
c       write(*,*)'eplus=',eplus,'eminus=',eminus,'epar=',epar
c     1 ,'dmode=',dmode
      return
      end











