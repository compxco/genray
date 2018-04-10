c        ********************** efield1***********************
c        * this subroutine calculates the components of  the *
c        * wave electric field using complex dielectric      *
c        * tensor -  reps(3,3),    complex N_perp - cnx	     *
c        * and real N_par- cnpar                             *
c        *****************************************************
c
c--------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      cnpar     double precision N_parallel                        !
c      cnx       double complex N_perp                              !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c								    !
c        output parameters					    !
c     ex,ey,ez				      			    !
c     double complex cex,cey,cez in common/cefield/		    !
c     eplus, eminus, epar - wave polarization 			    !
c--------------------------------------------------------------------
      subroutine efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'eps.i'
      include 'cefield.i'
      double complex dlamb(3,3),aa(3,3),
     1 ccnz,cnx,det,ci,cdelz,cp
  
       
c      write(*,*)'in efield1'
      ci=dcmplx(0.0d0,1.0d0)
      cnpar2=cnpar*cnpar
      cnper1=dreal(cnx)
      cnprim=dimag(cnx)
      cnper2=cnper1*cnper1+cnprim*cnprim
      dcn2=cnpar2+cnper2
c --------------------------------------------------------
c      write(*,*)'1 efield cnpar,cnper1,cnprim',cnpar,cnper1,cnprim
c      write(*,*)'reps',reps

      ccnz=dcmplx(cnpar,0.d0)

      dlamb(1,1)=ccnz*ccnz-reps(1,1)
      dlamb(1,2)=-reps(1,2)
      dlamb(1,3)=-ccnz*cnx-reps(1,3)
      dlamb(2,1)=-reps(2,1)
      dlamb(2,2)=ccnz*ccnz+cnx*cnx-reps(2,2)
      dlamb(2,3)=-reps(2,3)
      dlamb(3,1)=-ccnz*cnx-reps(3,1)
      dlamb(3,2)=-reps(3,2)
      dlamb(3,3)=cnx*cnx-reps(3,3)
c      do i=1,3
c        do j=1,3
c          write(*,*)'i,j,dlamb(i,j)',i,j,dlamb(i,j)
c        enddo
c      enddo
c-------------------------------------------------------------
c     cofactors
      aa(1,1)=dlamb(2,2)*dlamb(3,3)-dlamb(2,3)*dlamb(3,2)
      aa(1,2)=-dlamb(2,1)*dlamb(3,3)+dlamb(2,3)*dlamb(3,1)
      aa(1,3)=dlamb(2,1)*dlamb(3,2)-dlamb(2,2)*dlamb(3,1)
      aa(2,1)=-dlamb(1,2)*dlamb(3,3)+dlamb(3,2)*dlamb(1,3)
      aa(2,2)=dlamb(1,1)*dlamb(3,3)-dlamb(1,3)*dlamb(3,1)
      aa(2,3)=-dlamb(1,1)*dlamb(3,2)+dlamb(1,2)*dlamb(3,1)

      aa(3,1)=dlamb(1,2)*dlamb(2,3)-dlamb(2,2)*dlamb(1,3)
      aa(3,2)=-dlamb(1,1)*dlamb(2,3)+dlamb(2,1)*dlamb(1,3)
      aa(3,3)=dlamb(1,1)*dlamb(2,2)-dlamb(1,2)*dlamb(2,1)

c      do i=1,3
c         do j=1,3
c            write(*,*)'i,j,aa(i,j)',i,j,aa(i,j)
c         enddo
c      enddo
c-------------------------------------------------------------
c choice of max(abs(aa(i,j))

      i0=1
      j0=1
      dmaxaa=0.d00
      do 30 i=1,3
         do 40 j=1,3
            ppp=abs(aa(i,j))
            if(ppp.gt.dmaxaa) then
              dmaxaa=ppp
              i0=i
              j0=j
            end if
40       continue
30    continue
c---------------------------------------------------------------
c      write(*,*)'i0=',i0,'j0=',j0
      
      cex=aa(i0,1)
      cey=aa(i0,2)
      cez=aa(i0,3)

c
c      p=1.d0/dsqrt(cdabs(cex)*cdabs(cex)+cdabs(cey)*cdabs(cey)+
c     1   cdabs(cez)*cdabs(cez))
c      cex=cex*p
c      cey=cey*p
c      cez=cez*p

c      cdelz=cdabs(cez)/cez

cSm021021
c      thetz=datan2(dimag(cez),dabs(dreal(cez)))    
c      if (dreal(cez).ge.0.d0) then
c         p=1.d0
c      else
c         p=-1.d0
c      endif     
c      thetz=thetz*p
c      cdelz=dcmplx(dcos(thetz),-dsin(thetz))
         
cSm021021
c      write(*,*)'cex',cex
c      write(*,*)'cey',cey
c      write(*,*)'cez',cez
c      write(*,*)'dimag(cez),dreal(cez)',dimag(cez),dreal(cez)
c       thetz=datan2(dimag(cez),dreal(cez))
c      write(*,*)'1 thetz',thetz
c      write(*,*)'dimag(cez),cdabs(cez)',dimag(cez),cdabs(cez)
cSAP080228
c       thetz=dasin(dimag(cez)/cdabs(cez))
c      write(*,*)'2 thetz',thetz

c       cdelz=dcmplx(dcos(thetz),-dsin(thetz))

cSAP080422
c       if (dreal(cez).eq.0.d0) goto 10
        cp=dcmplx(dreal(cez),-dimag(cez))/abs(cez)
c        if (dreal(cez).lt.0d0) cp=-cp
        
cSAP080306   
c       cp=dcmplx(dreal(cez),-dimag(cez))/abs(cez)**2*
c     & dreal(cez)/abs(cez)

c       write(*,*)'cex,cey,cez',cex,cey,cez
c       write(*,*)' cp',cp

       cex=cex*cp
       cey=cey*cp
       cez=cez*cp
c       write(*,*)'efield1 cez',cez
c       write(*,*)'efield1 cex',cex 
c       write(*,*)'efield1 cey',cey

 10    continue
c       cdelz=1.d0
                            
c      if (j0.eq.1) cdelz=cdabs(cex)/cex
c      if (j0.eq.2) cdelz=cdabs(cey)/cey
c      if (j0.eq.3) cdelz=cdabs(cez)/cez
c-----------------------------------------------------------
cSAP080306
c      cp=cdelz/dsqrt(cdabs(cex)*cdabs(cex)+cdabs(cey)*cdabs(cey)+
c     1   cdabs(cez)*cdabs(cez))

      cp=1.d0/dsqrt(abs(cex)*abs(cex)+abs(cey)*abs(cey)+
     1   abs(cez)*abs(cez))    

      cex=cex*cp
      cey=cey*cp
      cez=cez*cp
c      write(*,*)'efield1 cez',cez
c      write(*,*)'efield1 cex',cex 
c      write(*,*)'efield1 cey',cey
      cez1=dsqrt(1.d0-cda bs(cex)**2-abs(cey)**2)
      ex=abs(cex)
      ey=abs(cey)
      ez=abs(cez)
      cmodee=ex*ex+ey*ey+ez*ez
c      write(*,*)'ex=',ex,'ey=',ey,'ez=',ez,'cmodee',cmodee
c      write(*,*)'thetz',thetz,'cdelz',cdelz,'cp=',cp
c      write(*,*)'efield1 cex=',cex,'cey=',cey,'cez=',cez,'cez1=',cez1
c      write(*,*)'efield1 cez,ez',cez,ez
c---------------------------------------------------------------
      root2=dsqrt(2.d0)
      eplus=abs(cex+ci*cey)/root2
      eminus=abs(cex-ci*cey)/root2
      dmode=eplus**2+eminus**2+ez**2
c      write(*,*)'eplus=',eplus,'eminus=',eminus,'epar=',epar
c     1 ,'dmode=',dmode
      epar=abs(cez*cnpar+cex*dreal(cnx))/dsqrt(cnpar**2+dreal(cnx)**2)
c      write(*,*)'efield1 cez,cex,cey',cez,cex,cey
c      write(*,*)'efield1 cnpar',cnpar
c      write(*,*)'efield1 cnx',cnx
c      write(*,*)'efield1 epar',epar

      return
      end











