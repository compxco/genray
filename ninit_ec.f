c        **********************ninit_ec************************
c        * this subroutine creates                            
c        * tangent to the magnetic surface, the components
c        * of the refractive index cnteta,cnphi                               
c        * at the initial point (zu0,ru0,phiu0) for ECR wave  
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters:                                         !
c              zu0,ru0,phiu0,cnx,cny,cnz			   !
c        output parameters:					   !
c              cnteta,cnphi                                        !
c------------------------------------------------------------------
        subroutine ninit_ec(zu0,ru0,phiu0,cnx,cny,cnz,cnteta,cnphi)
        implicit double precision (a-h,o-z)
c------------------------------------------------------------------
c       e_phi=e_psi x e_phi,   x -is the vector product.
c	           z       r      phi
c       e_phi= {   0   ,   0    ,   1    }
c       e_psi= { dpsidz, dpsidr ,   0    }/mod(grad(psi))
c	e_teta={ dpsidr,-dpsidz ,   0    }/mod(grad(psi))
c
c              	   x	                y                z
c       e_phi= {-sin(phiu0)      ,  cos(phiu0)     ,     0    }
c	e_teta={-cos(phiu0)dpsidz,-sin(phiu0)dpsidz,   dpsidr }
c
c       N_teta=(vector{N}*e_teta) ,N_phi=(vector{N}*e_phi)
c----------------------------------------------------------------------
      include 'param.i'
      include 'one.i'
c---------------------------------------------------------------------
c     dpdrd=dpsidr dpdzr=dpsidr were calculated by b(z,r,phi)
c     They are in  common/one/
c---------------------------------------------------------------------
      gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
      cosphi=dcos(phiu0)
      sinphi=dsin(phiu0)
      costet=dpdzd/gradpsi
      sintet=dpdrd/gradpsi
      cntet=-cnx*costet*cosphi-cny*costet*sinphi+cnz*sintet

c      write(*,*)'cntet=',cntet
c      write(*,*)'ninit_ec.f dpdzd,dpdrd,gradpsi',
c     +dpdzd,dpdrd,gradpsi
c      write(*,*)'cnx,cny,cnz',cnx,cny,cnz
c      write(*,*)'cosphi,sinphi',cosphi,sinphi

      cnteta=(-dpdzd*(cnx*cosphi+cny*sinphi)+dpdrd*cnz)/gradpsi
      cnphi=-cnx*sinphi+cny*cosphi

c      write(*,*)'ninit_ec.f cnteta,cnphi',cnteta,cnphi

      return
      end






