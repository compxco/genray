c        ********************** abc   ***********************
c        * this subroutine calculates the coefficients      *
c        * a,b,c, for the cold dispersion relation          *
c        * d=a*n**4+b*n**2+c                                *
c        ****************************************************
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c      	input parameters					    !	
c	point coordinates: z,r,phi    				    !
c       ds2=sin(gam)**2,dc2=cos(gam)**2                             !
c       gam is the angle between magnetic field and refractiv index !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      subroutine abc(z,r,phi,ds2,dc2,ad,bd,cd)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
c  ---------------------------------------------------------------
c  cold plasma dispersion relation with electrons and ions
c  Xe==x(i=1), Ye==y(i=1) for electron component, and 
c  Xi==x(i may be=2,nbulk), Yi==y(i may be=2,...,nbulk) for ions.
c  ib is the species number for which delib=1-Yib is expected to be
c  zero within plasma; ib can be 1 (if ECR is expected) 
c  or other ib>1 (if ICR is expected).
c  The dispersion coefficients A,B,C 
c  (see Genray manual: A,B,C are given by (4.9-4.11))  
c  are multiplied by dele=1-Ye, in case of ib=1,
c  or by deli=1-Yi, in case of ib=i (ib>1, in case of ICR)
c  --------------------------------------------------------------
      pc=1.d0+dc2
      call s(z,r,phi,s1,s2,s3,s4,s6,s7)
      xib=x(z,r,phi,ib) !here ib can be 1 (electrons) or >1 (ions)
      yib=y(z,r,phi,ib)
      delib=1.d0-yib ! 1-Yib for a given ib, it is either dele or deli

c----------------------------------------------------------------
c  ib =1 (the cyclotron resonance conditions dele=0 may be
c         realised in plasma, dispersion relation is
c         multiplied by dele )
c         ad,bd,cd calculations
c------------------------------------------------------------------

c  if 2 begin
      if (ib.eq.1) then
         xe=xib
         ye=yib
         !peym=xe/(1.d0-ye) !YuP[07-2017] commented: not used in this case; can be 1/0
         peyp=xe/(1.d0+ye)
         a0e=-peyp*ds2
         a1e=s7*ds2+s4*dc2
         b0e=s4*peyp*pc+xe*(s6-peyp)*ds2
         b1e=-s4*s7*pc-s3*(s6-peyp)*ds2
         c0e=-xe*s4*(s6-peyp)
         c1e=s4*s3*(s6-peyp)
         dele=1.d0-ye ! or simply dele=delib, from the above
         ad=dele*a1e+a0e !ad=(1-Ye)*A,  A,B,C given by (4.9-4.11) in Genray manual
         bd=dele*b1e+b0e !bd=(1-Ye)*B
         cd=dele*c1e+c0e !cd=(1-Ye)*C
      end if
c  if 2 end
c-----------------------------------------------------------------
c  ib .gt.1 (the cyclotron resonance conditions delib=0 may be
c            realised in plasma, dispersion relation is
c            multiplied by delib )
c            ad,bd,cd calculations
c------------------------------------------------------------------
c  if 3 begin
      if (ib.gt.1) then ! one of ions can have ICR
         ! here deli=delib
         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
         ppe=1.d0/(1.d0-ye)
         peym=xe/(1.d0-ye)
         peyp=xe/(1.d0+ye)
         peym2=peyp*ppe ! == Xe/(1-Ye^2)

         !ppb=1.d0/(1.d0-yib) !YuP but this can be 1/0 ?
         !pbym=xib/(1.d0-yib) !YuP but this can be 1/0 ? Not used 
         pbyp=xib/(1.d0+yib)
         !pbym2=pbyp/(1.d0-yib) !YuP but this can be 1/0 ? Not used

         a0b=-pbyp*ds2              != -[Xib/(1+Yib)]*sin^2()
         a1b=(s1-peym2)*ds2+s4*dc2  != (s1 -Xe/(1-Ye^2))*sin^2() + s4*cos^2()
         b0b=s4*pbyp*pc+xib*(s3-peym)*ds2              
           != s4*(Xib/(1+Yib))*(1+cos^2) + Xib*(s3 -Xe/(1-Ye))*sin^2()
         b1b=-s4*(s1-peym2)*pc-(s2-peyp)*(s3-peym)*ds2 
           != -s4*(s1- Xe/(1-Ye^2))*(1+cos^2) - (s2 -Xe/(1+Ye))*(s3 -Xe/(1-Ye))*sin^2()
         c0b=-xib*s4*(s3-peym)      != -Xib*s4*(s3- Xe/(1-Ye))
         c1b=s4*(s2-peyp)*(s3-peym) != s4*(s2- Xe/(1+Ye))*(s3- Xe/(1-Ye))
         delib=1.d0-yib
         ad=delib*a1b+a0b !ad=(1-Ye)*A,  A,B,C given by (4.9-4.11) in Genray manual
         bd=delib*b1b+b0b
         cd=delib*c1b+c0b
      end if
c  if 3 end

!   A problem noticed by YuP[07-2017]:
!   In equation for the two roots,    
!	    N2p=(-B +sqrt(B^2-4AC))/(2A)      (4.12a) (O)
!	    N2m=(-B -sqrt(B^2-4AC))/(2A)      (4.12b) (X)
!   the sign in front of sqrt() [defined as ioxm] determines the mode: 
!   '+' is for O-mode, '-' is for X-mode .
!   However, after we multiplied A,B,C by (1-Yib) factor,
!   the result depends on the sign of delib=(1-Yib).
!   If (1-Yib) is positive, nothing changes 
!   in correspondence of ioxm=+/-1 and the two modes. 
!   But for the negative (1-Yib), e.g. (1-Yib)=-1, we get
!      (-b +ioxm*sqrt(b^2-4ac))/(2a) = [use a=(1-Yib)*A= -A, etc] 
!    = (+B +ioxm*sqrt(B^2-4AC))/(-2A)= 
!    = (-B -ioxm*sqrt(B^2-4AC))/(+2A)
!   Thus, for ioxm=+1, we are getting the branch (4.12b), 
!   which is the X mode. The meaning of modes is reversed !
!   To correct this problem, we simply need to further adjust 
!   the a,b,c cofficients:
!      sign_delib=sign(1.d0,delib) 
!      a=a*sign_delib 
!      b=b*sign_delib 
!      c=c*sign_delib 
!   So, effectively, we use a=|1-Yib|*A, etc.
!   Then, the meaning (mode type) defined through a,b,c
!   will remain the same as that defined through the original A,B,C.
      sign_delib=sign(1.d0,delib) ! YuP[07-2017]
      ad=ad*sign_delib ! YuP[07-2017] adjusted
      bd=bd*sign_delib ! YuP[07-2017] adjusted
      cd=cd*sign_delib ! YuP[07-2017] adjusted
!   Note that the solution of a*n**4 +b*n**2 +c= 0 equation
!   only depends on gamma angle, and not on value of Npar ! 
!   These coefficients (a,b,c)==(ad,bd,cd) are further used to find
!   the solution in form of N^2(gamma^2).
!   Important is to remember that the two solutions, 
!	    N2p=(-B +sqrt(B^2-4AC))/(2A)      (4.12a) (O)
!	    N2m=(-B -sqrt(B^2-4AC))/(2A)      (4.12b) (X)
!   will have different values of Npar (even by absolute value)

      return
      end


c        ********************** abc_test   ******************
c        * this subroutine calculates the coefficients      *
c        * a,b,c, for the cold dispersion relation          *
c        * d=a*n**4+b*n**2+c                   
c        * and prints some variables for test               *
c        ****************************************************
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c      	input parameters					    !	
c	point coordinates: z,r,phi    				    !
c       ds2=sin(gam)**2,dc2=cos(gam)**2                             !
c       gam is the angle between magnetic field and refractiv index !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      subroutine abc_test(z,r,phi,ds2,dc2,ad,bd,cd)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
c  ---------------------------------------------------------------
c  cold plasma dispersion relation with electrons and ions
c  x(i=1),y(i=1) for electron component
c  x(i may be=2,nbulk),y(i may be=2,...,nbulk) ions components.
c  ib is a number of component for which delib=1-yib may be equal
c  zero inside the plasma,ib may be=from 1 to nbulk
c  The dispersion relation is multiplied by delib
c  --------------------------------------------------------------

c test rho
c      bmod=b(z,r,phi)
c      write(*,*)'z,r,phi,rho',z,r,phi,rho
c end test

      pc=1.d0+dc2
      call s(z,r,phi,s1,s2,s3,s4,s6,s7)
      xib=x(z,r,phi,ib)
      yib=y_test(z,r,phi,ib)
      
      delib=1.d0-yib

c----------------------------------------------------------------
c  ib =1 (the cyclotron resonance conditions dele=0 may be
c         realised in plasma, dispersion relation is
c         multiplied by dele )
c         ad,bd,cd calculations
c------------------------------------------------------------------

c  if 2 begin
      if (ib.eq.1) then
         xe=xib
         ye=yib
         !peym=xe/(1.d0-ye)!YuP[07-2017] commented: not used in this case; can be 1/0
	 peyp=xe/(1.d0+ye)
	   a0e=-peyp*ds2
	   a1e=s7*ds2+s4*dc2
	   b0e=s4*peyp*pc+xe*(s6-peyp)*ds2
	   b1e=-s4*s7*pc-s3*(s6-peyp)*ds2
	   c0e=-xe*s4*(s6-peyp)
	   c1e=s4*s3*(s6-peyp)
	   dele=1.d0-ye
	     ad=dele*a1e+a0e
	     bd=dele*b1e+b0e
	     cd=dele*c1e+c0e
      end if
c  if 2 end
c-----------------------------------------------------------------
c  ib .gt.1 (the cyclotron resonance conditions delib=0 may be
c            realised in plasma, dispersion relation is
c            multiplied by delib )
c            ad,bd,cd calculations
c------------------------------------------------------------------
c  if 3 begin
      if (ib.gt.1) then
         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
         ppe=1.d0/(1.d0-ye)
         peym=xe/(1.d0-ye)
	 peyp=xe/(1.d0+ye)
	 peym2=peyp*ppe

         !ppb=1.d0/(1.d0-yib)
         !pbym=xib/(1.d0-yib)!YuP but this can be 1/0 ? Not used
	 pbyp=xib/(1.d0+yib)
         !pbym2=pbyp/(1.d0-yib) !YuP but this can be 1/0 ? Not used

	   a0b=-pbyp*ds2
	   a1b=(s1-peym2)*ds2+s4*dc2
	   b0b=s4*pbyp*pc+xib*(s3-peym)*ds2
	   b1b=-s4*(s1-peym2)*pc-(s2-peyp)*(s3-peym)*ds2
	   c0b=-xib*s4*(s3-peym)
	   c1b=s4*(s2-peyp)*(s3-peym)
	   delib=1.d0-yib
	     ad=delib*a1b+a0b
	     bd=delib*b1b+b0b
	     cd=delib*c1b+c0b

c             write(*,*)'ds2,dc2,delib,c1b,c0b,cd',
c     &                  ds2,dc2,delib,c1b,c0b,cd
c            write(*,*)'for c1b=s4*(s2-peyp)*(s3-peym)'
c            write(*,*)'s4,s2,peyp,s3,peym,xe,ye,yib',
c     &                 s4,s2,peyp,s3,peym,xe,ye,yib
c            write(*,*)'rho',rho,'densrho(rho,1)',densrho(rho,1)
c            write(*,*)'x_test(z,r,phi,1)',x_test(z,r,phi,1)
c            write(*,*)'(s2-peyp),(s3-peym),(s2-peyp)*(s3-peym)',
c     &                 (s2-peyp),(s3-peym),(s2-peyp)*(s3-peym)
      end if
c  if 3 end

      return
      end


