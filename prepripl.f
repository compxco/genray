c        ********************** prepripl***************************
c        *                      -----                             *
c        *  prepripl -subroutine to prepar the  data for styding  *
C        *  the influence of the ripple magnetic field on N_par   *
c        **********************************************************
      subroutine prepripl(t,u,is)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'write.i'
      include 'three.i'
      include 'n_parb.i'
      double precision modr,modrold
      dimension u(6)
c---------------------------------------------
c     z,r (m)
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)

c---------------------------------------
c     the creation of the array wthet(i): poloidal angle
c     along the ray trajectory (ir radians)
      zcomp=z-yma
      rcomp=r-xma
      modr=dsqrt(zcomp*zcomp+rcomp*rcomp)  !=|r|
cc    write(*,*)'in prepripl,is,zcomp,rcomp,z,yma,r,xma,modr',
cc     1	  is,zcomp,rcomp,z,yma,r,xma,modr
      btor(is)=bphi
      btot(is)=bmod
      gradpsi(is)=dsqrt(dpdzd**2+dpdrd**2)
      gradpdr(is)=gradpsi(is)/r
c--------------------------------------------------------------------
c     e_theta=e_psi x e_phi		   (x is a vector production)
c     e_psi=(e_z*dpsidz+e_r*dpsidr)/sqrt(dpsidz**2+dpsidr**2)
c     e_theta=(e_z*dpsidr-e_r*dpsidz)/sqrt(dpsidz**2+dpsidr**2)
c     Npol=N * e_theta
c--------------------------------------------------------------------
c      write(*,*)'bz,br,bphi',bz,br,bphi
      bpoloid(is)=(bz*dpdrd-br*dpdzd)/gradpsi(is)
c     write(*,*)'cnz,cnr',cnz,cnr
c     write(*,*)'dpdzd,dpdrd,gradpsi(is)',dpdzd,dpdrd,gradpsi(is)
      wnpol(is)=(cnz*dpdrd-cnr*dpdzd)/gradpsi(is)
c      write(*,*)'wnpol,bpoloid,btot',wnpol(is),bpoloid(is),btot(is)
c     read(*,*)
c--------------------------------------------------------------------
      if (is.eq.1) then
         if(zcomp.ge.0.d0) then
            wthet(1)=dacos(rcomp/modr)
         else
            wthet(1)=-dacos(rcomp/modr)
         endif
      else
cc       write(*,*)'in prepripl,zold,rold,xma,yma'
cc    1                        ,zold,rold,xma,yma
         zoldcomp=zold-yma
         roldcomp=rold-xma
	 modrold=dsqrt(zoldcomp*zoldcomp+roldcomp*roldcomp) !|r_old|
cc	 write(*,*)'in prepripl,is,zoldcomp,roldcomp,modrold',
cc     1 		    is,zoldcomp,roldcomp,modrold
c--------------------------------------------------------------------
c            [r*r_old]=|r|*|r_old|*sin(deltheta)
c            |e_z        e_r       e_phi |
c            |zcomp	 rcomp	   0	 |=e_phi*|r|*||r_old|*sin(deltheta)
c            |zoldcomp   roldcomp  0	 |
	 sindelth=(zcomp*roldcomp-rcomp*zoldcomp)/(modr*modrold)
c--------------------------------------------------------------------
c            (r*r_old)=|r|*|r_old|*cos(delttheta)
	 cosdelth=(zcomp*zoldcomp+rcomp*roldcomp)/(modr*modrold)
	 if(cosdelth.gt.1.d0-1.d-15) cosdelth=1.d0-1.d-15
	 if(cosdelth.lt.-1.d0+1.d-15) cosdelth=-1.d0+1.d-15
cc	 trigtest=sindelth**2+cosdelth**2
cc	 write(*,*)'sindelth,cosdelth,trigtest',
cc     1 		sindelth,cosdelth,trigtest
c--------------------------------------------------------------------
         if (sindelth.gt.0.d0) then
	    deltheta=dacos(cosdelth)
	 else
	    deltheta=-dacos(cosdelth)
	 endif
	 wthet(is)=wthet(is-1)+deltheta
      endif ! is
cc     write(*,*)'r=',r
      wmdevr(is)=cm/r*bphi/bmod
cc    write(*,*)'wthet(is)',wthet(is),'deltheta',deltheta
c--------------------------------------------------------------------
      if (istart.eq.1) goto 10
c     Determination of the npar boundary. Report of the APS Meeting 1994
c     LH ray tracing studies with reconstructed magnetic equilibria
c     in PBX-M. F.Paoletti, S,Bernabei, D.Ignat, R.Kaita, J.Kersner,
c     B.Le Blanc, F.Levinton, S.Luckhurdt
      xe=x(z,r,phi,1)
      wxe(is)=dsqrt(xe)
c     wnrho is the radial component of the refructive index.
c     It is positive if it is directed along the gradient 
c     of the flux surface. 
      wnrho(is)=(cnr*dpdrd+cnz*dpdzd)/gradpsi(is)
      gnpar(is)=bpoloid(is)/btot(is)
      wp(is)=gnpar(is)*wxe(is)
      p2=wp(is)*wp(is)
      pp=(1.d0+wp(is))*(1.d0-wp(is))
      det=1.d0-pp*(wnrho(is)/(cm/r)**2)/xe
      cnphi=cm/r
      if (det.gt.0.d0) then
         det=cnphi*wp(is)*dsqrt(det)
	 if (pp.ne.0.d0) then
	    wnparplb(is)=-(cnphi-det)/pp
	    wnparmnb(is)=-(cnphi+det)/pp
      write(*,*)'cnphi,det,pp',cnphi,det,pp
	 else			       
	    write(*,*)'in prepripl (1-wp)(1+wp)=0'
	    wnparplb(is)=10000.d0
	    wnparmnb(is)=10000.d0
	 endif  
      else
	 write(*,*)'in prepripl det.lt.0'
	 wnparplb(is)=-10000.d0
	 wnparmnb(is)=-10000.d0
      endif
      write(*,*)'wnrho(is)',wnrho(is),'cnphi',cm/r,'gnpar(is)',gnpar(is)
      write(*,*)'sqrt(xe)',wxe(is),'wp(is)',wp(is),'cnphi',cnphi
      write(*,*)'wnparplb(is)',wnparplb(is),
     1          'wnparmnb(is)',wnparmnb(is),
     1          'npar',(cnz*bz+cnr*br+cnphi*bphi)/bmod
      wnparplb(is)=-cnphi/(1.d0+wp(is))
      wnparmnb(is)=-cnphi/(1.d0-wp(is))
c      write(*,*)'nmin',-cnphi/(1.d0+wp(is)),'nmax',-cnphi/(1.d0-wp(is))
      write(*,*)'nmin',wnparplb(is),'nmax',wnparmnb(is)
c     control of the electrostatic limit for the dispersion Hamiltonian
      xi=x(z,r,phi,2)
      ye=y(z,r,phi,1)
      sh=1.d0-xi+xe/(xe*ye)
      ph=1.d0-xe
      cnparh=(cnz*bz+cnr*br+cnphi*bphi)/bmod
      cnparh2=cnparh*cnparh
      cnperh2=cnz*cnz+cnr*cnr+cnphi*cnphi-cnparh2
      hamelst=cnperh2*sh+cnparh2*ph
      write(*,*)'sh',sh,'ph',ph,'cnparh2',cnparh2,'cnperh2',cnperh2
      write(*,*)'in prepripl electrostaic limit hamiltonian',hamelst
c--------------------------------------------------------------------
 10   continue
      return
      END

c*****************************************************
c*----------------subroutine GRAPHnpr----------------*
c* Prepaires file drawgenr.in for xdraw              *
c* INPUT: from common blocks'gr.i''n_parb.i' and   *
c         'write'                                    *
c*****************************************************
      SUBROUTINE GRAPHnpr
      IMPLICIT double precision (a-h,o-z)
      INTEGER I
      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      include 'gr.i'
      include 'write.i'
      include 'n_parb.i'
      include 'three.i'

      if(myrank.ne.0) return
      
      OPEN(11,file='npar.bin',form='unformatted')

      DO 10 I=1,nrayelt
        WRITE(11) REAL(wthet(I)),REAL(wr(I)),REAL(wz(I)),
     1  REAL(ws(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(wmdevr(I)),
     3  REAL(wnpar(I)-wmdevr(I))
10    CONTINUE

      WRITE(11)

      CLOSE(11)
13	Format(16g15.5)
14	Format(14e11.4)
15	Format(12e12.4)

      OPEN(21,file='npar.sap')
      OPEN(23,file='nparadd.sap')
 1    format(' ',i3)
      write(21,*)'wthet wr wz ws psi npar nper npol
     1 bpol btor btot rex imey rez flux wmtor '
      write(23,*)'delpwr phi wx wy xmaga ymaga wnrho gnpar wxe wp 
     1wnparplb wnparmnb '
      DO 120 I=1,nrayelt
       WRITE(21,13) wthet(I),wr(I),wz(I),ws(I),spsi(I),
     2 wnpar(I),wnper(I),wnpol(I),
     3 bpoloid(I),btor(I),btot(I),
     4 dreal(cwexde(I)),dimag(cweyde(I)),dreal(cwezde(I)),fluxn(I),
     5 wmtor(I)
       WRITE(23,15) delpwr(I),wphi(I),wr(I)*cos(wphi(I)),
     1 wr(I)*sin(wphi(I)),100.d0*xma*dcos(wphi(I)),
     2 100.d0*xma*dsin(wphi(I)),
     3 wnrho(I),gnpar(I),wxe(I),wp(I),wnparplb(I),wnparmnb(I)
120    CONTINUE
      CLOSE(21)
      CLOSE(23)

c-------------------------------------------- 
c     5 contours psi=const
      OPEN(24,file='section.sap')
      write(24,*)
16    Format(10e12.4)
      do i=1,NP+1
        WRITE(24,16)AR(1,i),AZ(1,i),AR(2,i),AZ(2,i),AR(3,i),AZ(3,i),
     1              AR(4,i),AZ(4,i),AR(5,i),AZ(5,i)
      enddo
      CLOSE(24)
c-------------------------------------------- 
      RETURN
      END



c        ********************** prepebw***************************
c        *                      -----                             *
c        *  prepebw -subroutine to prepar the  data for studying   *
C        *  ebw                                                   *
c        **********************************************************
      subroutine prepebw(t,u,is)
      implicit double precision (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'write.i'
      include 'three.i'
      include 'n_parb.i'
      double precision modr,modrold
      dimension u(6)
c---------------------------------------------
c     z,r (m)
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)

c---------------------------------------
c     the creation of the array wthet(i): poloidal angle
c     along the ray trajectory (ir radians)
      zcomp=z-yma
      rcomp=r-xma
      modr=dsqrt(zcomp*zcomp+rcomp*rcomp)  !=|r|

      btor(is)=bphi
      btot(is)=bmod
      gradpsi(is)=dsqrt(dpdzd**2+dpdrd**2)
      gradpdr(is)=gradpsi(is)/r
c--------------------------------------------------------------------
c     e_theta=e_psi x e_phi		   (x is a vector production)
c     e_psi=(e_z*dpsidz+e_r*dpsidr)/sqrt(dpsidz**2+dpsidr**2)
c     e_theta=(e_z*dpsidr-e_r*dpsidz)/sqrt(dpsidz**2+dpsidr**2)
c     Npol=N * e_theta
c--------------------------------------------------------------------
c      write(*,*)'prepripl_ebw: bz,br,bphi',bz,br,bphi
      bpoloid(is)=(bz*dpdrd-br*dpdzd)/gradpsi(is)
c      write(*,*)'prepripl_ebw cnz,cnr,cnphi',cnz,cnr,cm/r
c      write(*,*)'prepripl_ebw dpdzd,dpdrd,gradpsi(is)',
c     1dpdzd,dpdrd,gradpsi(is)
      wnpol(is)=(cnz*dpdrd-cnr*dpdzd)/gradpsi(is)
c      write(*,*)'prepripl_ebw wnpol(is)',wnpol(is)
c      write(*,*)'prepripl_ebwwnpol,bpoloid,btot',
c     1wnpol(is),bpoloid(is),btot(is)
c--------------------------------------------------------------------
      if (is.eq.1) then
         if(zcomp.ge.0.d0) then
            wthet(1)=dacos(rcomp/modr)
         else
            wthet(1)=-dacos(rcomp/modr)
         endif
      else

         zoldcomp=zold-yma
         roldcomp=rold-xma
	 modrold=dsqrt(zoldcomp*zoldcomp+roldcomp*roldcomp) !|r_old|

c--------------------------------------------------------------------
c            [r*r_old]=|r|*|r_old|*sin(deltheta)
c            |e_z        e_r       e_phi |
c            |zcomp	 rcomp	   0	 |=e_phi*|r|*||r_old|*sin(deltheta)
c            |zoldcomp   roldcomp  0	 |
	 sindelth=(zcomp*roldcomp-rcomp*zoldcomp)/(modr*modrold)
c--------------------------------------------------------------------
c            (r*r_old)=|r|*|r_old|*cos(delttheta)
	 cosdelth=(zcomp*zoldcomp+rcomp*roldcomp)/(modr*modrold)
	 if(cosdelth.gt.1.d0-1.d-15) cosdelth=1.d0-1.d-15
	 if(cosdelth.lt.-1.d0+1.d-15) cosdelth=-1.d0+1.d-15

c--------------------------------------------------------------------
         if (sindelth.gt.0.d0) then
	    deltheta=dacos(cosdelth)
	 else
	    deltheta=-dacos(cosdelth)
	 endif
	 wthet(is)=wthet(is-1)+deltheta
      endif ! is

      wmdevr(is)=cm/r*bphi/bmod

c--------------------------------------------------------------------
      

      xe=x(z,r,phi,1)
      wxe(is)=dsqrt(xe)
c     wnrho is the radial component of the refructive index.
c     It is positive if it is directed along the gradient 
c     of the flux surface. 
      wnrho(is)=(cnr*dpdrd+cnz*dpdzd)/gradpsi(is)
      gnpar(is)=bpoloid(is)/btot(is)
      wp(is)=gnpar(is)*wxe(is)
     
      cnphi=cm/r
     
 10   continue
      return
      END

c*****************************************************
c*----------------subroutine GRAPHebw----------------*
c* Prepares for file drawgenr.in for xdraw           *
c* INPUT: from common blocks'gr.i''n_parb.i' and     *
c         'write'                                    *
c*****************************************************
      SUBROUTINE GRAPHebw
      IMPLICIT double precision (a-h,o-z)
      INTEGER I
      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'gr.i'
      include 'write.i'
      include 'n_parb.i'
      include 'three.i'
      
      if(myrank.ne.0) return

      OPEN(11,file='npar.bin',form='unformatted')

      DO 10 I=1,nrayelt
        WRITE(11) REAL(wthet(I)),REAL(wr(I)),REAL(wz(I)),
     1  REAL(ws(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(wmdevr(I)),
     3  REAL(wnpar(I)-wmdevr(I))
10    CONTINUE

      WRITE(11)

      CLOSE(11)
13	Format(16g15.5)
14	Format(14e11.4)
15	Format(12e12.4)

      OPEN(21,file='npar.sap')
      OPEN(23,file='nparadd.sap')
 1    format(' ',i3)
      write(21,*)'wthet wr wz ws psi npar nper npol
     1 bpol btor btot rex imey rez flux wmtor '
      write(23,*)'delpwr phi wx wy xmaga ymaga wnrho gnpar wxe wp 
     1wnparplb wnparmnb '
      DO 120 I=1,nrayelt
       WRITE(21,13) wthet(I),wr(I),wz(I),ws(I),spsi(I),
     2 wnpar(I),wnper(I),wnpol(I),
     3 bpoloid(I),btor(I),btot(I),
     4 dreal(cwexde(I)),dimag(cweyde(I)),dreal(cwezde(I)),fluxn(I),
     5 wmtor(I)
       WRITE(23,15) delpwr(I),wphi(I),wr(I)*cos(wphi(I)),
     1 wr(I)*sin(wphi(I)),100.d0*xma*dcos(wphi(I)),
     2 100.d0*xma*dsin(wphi(I)),
     3 wnrho(I),gnpar(I),wxe(I),wp(I),wnparplb(I),wnparmnb(I)
120    CONTINUE
      CLOSE(21)
      CLOSE(23)

c-------------------------------------------- 
c     5 contours psi=const
      OPEN(24,file='section.sap')
      write(24,*)
16    Format(10e12.4)
      do i=1,NP+1
        WRITE(24,16)AR(1,i),AZ(1,i),AR(2,i),AZ(2,i),AR(3,i),AZ(3,i),
     1              AR(4,i),AZ(4,i),AR(5,i),AZ(5,i)
      enddo
      CLOSE(24)
c-------------------------------------------- 
      RETURN
      END



