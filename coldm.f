

      subroutine coldm(nz,az2,x,y,an2,ncld,icuto)

      implicit double precision (a-h,o-z)
c
c  input:
c     nz     = parallel refractive index
c     az2    = parallel refractive index squared
c     x      = (omegape/omega)**2
c     y      = (omegace/omega)
c  output:
c     an2(2) = X- and O-mode (perp. ref. index)**2
c     ncld(2)= X- and O-mode perp. ref. index, respectively.
c     icuto(2)=icuto(1)=1 if an2(1).lt.0., else 0
c              icuto(2)=2 if an2(2).lt.0., else 0
c
c
c       costruisce l'indice freddo per i due modi
c
      double precision ncld(2),nz,an2(2)
      dimension icuto(2)
c	"""""""""""""""""""""""""""
c	write(*,*)'beg. of coldm'
c	"""""""""""""""""""""""""""

      icuto(1)=0
      icuto(2)=0
      y2=y*y
      de=1.d0-x-y2
      d=y2*(1.d0-az2)**2+4.d0*az2*(1.d0-x)
      if(d.eq.0.0d0.or.d.gt.0.0d0) go to 1000
      icuto(1)=1
      icuto(2)=2
      go to 2000
1000  dd=dsqrt(d)
      b1=1.d0-az2-x-.5d0*x*y2*(1.d0+az2)/de
      b2=.5d0*x*y*dd/de
      an2(1)=b1-b2
      an2(2)=b1+b2
      if(an2(1).lt.0.d0) icuto(1)=1
      if(an2(2).lt.0.d0) icuto(2)=2
2000  if(icuto(1).eq.1) an2(1)=1.d0
      if(icuto(2).eq.2) an2(2)=1.d0
      ncld(1)=-1.d0
      ncld(2)=-1.d0
      if(an2(1).ge.0.)  ncld(1)=dsqrt(an2(1))
      if(an2(1).ge.0.)  ncld(2)=dsqrt(an2(2))
c      if (icuto(1).eq.1) write(*,*) 'cutoff modo  ',icuto(1)
c      if (icuto(2).eq.2) write(*,*) 'cutoff modo  ',icuto(2)

c	"""""""""""""""""""""""""""""""""""""
c	write(*,*)'end coldm'
c	"""""""""""""""""""""""""""""""""""""
      return
      end

