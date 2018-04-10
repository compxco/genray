
c the data for the contours flux_surface=psi(r,z)=const 
c and for plotting of the toroidal and poloidal tokamak cross sections

c      parameter (nteta,npsi,epspsi,nteta1=nteta+1)
c      PARAMETER (NL,NP=nteta)
      real*8 zpsi(npsi,nteta1),rpsi(npsi,nteta1),
     +                 zpsi1(npsi,nteta1),rpsi1(npsi,nteta1),
     +                 arpsi(npsi),arteta(nteta1),
cSAP091201 for spline approximation using coeff2
     &zpsi_psps(npsi,nteta1),zpsi_thth(npsi,nteta1),
     &zpsi_pspsthth(npsi,nteta1),
     &rpsi_psps(npsi,nteta1),rpsi_thth(npsi,nteta1),
     &rpsi_pspsthth(npsi,nteta1)

      real*8 AR(NL,NP+1),AZ(NL,NP+1)
      real*8 XT(3,NP+1),YT(3,NP+1)
      real*8 ar_min(npsi),ar_max(npsi),ab_max(npsi),
     +ab_min(npsi)

      common/gr/ AR,AZ,XT,YT,zpsi,rpsi,arpsi,arteta,
     +           ar_min,ar_max,ab_max,ab_min,
     &zpsi_psps,zpsi_thth,zpsi_pspsthth,
     &rpsi_psps,rpsi_thth,rpsi_pspsthth

c the arrays for contours of the flux surface: psi(r,z)=const
c npsi   is a number of contours
c nteta  is a number of the points along the each contours
c                       (in the poloidal direction)
c epspsi is the accuracy for the determination of the contours
c         poins coordinates	(zpsi,rpsi)
c zpsi(npsi,nteta1),rpsi(npsi,nteta1) arrays of the contour point's
c  		  coordinates
c arpsi(npsi) is the array of values the poloidal flux on the contours
c arteta(nteta1) is the array of the poloidal angles
c NL is a number of contours for ploting of the tokamak
c              cross section
c NP is a number of points along each contours in poloidal direction
c       for plotting of the tokamak cross section
c AR(NL,NP+1),AZ(NL,NP+1) arrays of the contour point's
c  		  coordinates	for the ploting of the tokamak cross section
c XT(3,NP+1),YT(3,NP+1) (x=r*cos(phi),y=r*sin(phi)) are the
c     coordinates for ploting of the middleplane tokamak cross section
c ar_min(npsi) is the array of the min value r(psi)
c ar_max(npsi) is the array of the max value r(psi)
c ab_max(npsi) is the array of the max value b_total(psi)
c ab_min(npsi) is the array of the min value b_total(psi)
