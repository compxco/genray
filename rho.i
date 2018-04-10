
c The arrays for the spline coefficients
c for the different forms of the small radius. 
c      parameter (npsi4=npsi+4,nteta14=nteta1+4)
c npsi  and nteta1 are in param.i
c      parameter(nzy=max0(npsi,nteta1)+4)
c      parameter(nzy=nteta1+4)
      real*8 
     1		  areatot,voltot,torftot,totlength,
     1            txz,tyz,czy,
     1            czxy,
     1            txr,tyr,
     1            crxy,
     1            trhos,crhspsi,
     1            trhov,crhvpsi,
     1            trhobt,crhbtpsi,
     1            trhos2,crhsps2,
     1            trhov2,crhvps2,
     1            trhobt2,crhbtps2,
     1            trhos_rh,crhsrho,
     1            trhov_rh,crhvrho,
     1		  tpsi,cpsirho,
     1		  trminpsi,crminpsi,
     1		  trmaxpsi,crmaxpsi,
     1		  tbmaxpsi,cbmaxpsi,
     1		  tbminpsi,cbminpsi,
     1            trholpol,crhlpsi,
     1            trhol_rh,crhlrho,
cSm050302
     1            trhor,crhrpsi,
     1            trhor2,crhrps2,
     1            trhor_rh,crhrrho,
cSAP080815 for zcunix spline 
     &            d2_rhos_psi,
     &            d2_rhotf_psi,
     &            d2_rhovl_psi,
     &            d2_rhol_psi,
     &            d2_rhor_psi,
     &            arrho_s,arrho_tf,arrho_vl,arrho_r,arrho_l,
cSAP0891122       for B_line length,
     &            length_b_ar,d2_length_b_psi

      common/rho/
     1		  areatot,voltot,torftot,totlength,
     1            txz(npsi4),tyz(nteta14),czy(nzy),
     1            czxy(npsi4,nteta14),
     1            txr(npsi4),tyr(nteta14),
     1            crxy(npsi4,nteta14),
     1            trhos(npsi4),crhspsi(npsi4),
     1            trhov(npsi4),crhvpsi(npsi4),
     1            trhobt(npsi4),crhbtpsi(npsi4),
     1            trhos2(npsi4),crhsps2(npsi4),
     1            trhov2(npsi4),crhvps2(npsi4),
     1            trhobt2(npsi4),crhbtps2(npsi4),
     1            trhos_rh(npsi4),crhsrho(npsi4),
     1            trhov_rh(npsi4),crhvrho(npsi4),
     1		  tpsi(npsi4),cpsirho(npsi4),
     1		  trminpsi(npsi4),crminpsi(npsi4),
     1		  trmaxpsi(npsi4),crmaxpsi(npsi4),
     1		  tbmaxpsi(npsi4),cbmaxpsi(npsi4),
     1		  tbminpsi(npsi4),cbminpsi(npsi4),
     1            trholpol(npsi4),crhlpsi(npsi4),
     1            trhol_rh(npsi4),crhlrho(npsi4),
cSm050302
     1            trhor(npsi4),crhrpsi(npsi4),
     1            trhor2(npsi4),crhrps2(npsi4),
     1            trhor_rh(npsi4),crhrrho(npsi4),
cSAP080815 for zcunix spline 
     &            d2_rhos_psi(npsi),
     &            d2_rhotf_psi(npsi),
     &            d2_rhovl_psi(npsi),
     &            d2_rhol_psi(npsi),
     &            d2_rhor_psi(npsi),
     &            arrho_s(npsi),
     &            arrho_tf(npsi),
     &            arrho_vl(npsi),
     &            arrho_r(npsi), 
     &            arrho_l(npsi),
cSAP081122 for B line length
     &            length_b_ar(npsi),d2_length_b_psi(npsi)

c   areatot -total area of the poloidal cross section
c   voltot  -total volume of the plasma cord
c   torftot -total toroidal flux in the plasma cord
c   totlength -the length of the intersection of the lust closed flux surface
c              with the poloidal verticle plane  

c  coefficients for the spline approximation of
c  z(psi,theta)        : txz,tyz,czxy  ,(czy-working array)
c  r(psi,theta)        : txr,tyrz,crxy ,(czy-working array)

c  rho_s(psi)          : tros,crhspsi
c  rho_v(psi)          : trhov,crhvpsi
c  rho_bt(psi)         : trhobt,crhbtpsi
c  (rho_s)**2(psi)     : trhos2,crhsps2
c  (rho_v)**2(psi)     : trhov2,crhvps2
c  (rho_bt)**2(psi)    : trhobt2,crhbtps2
c  rho_s(rho)          : trhos_rh,crhsrho
c  rho_v(rho)          : trhov_rh,crhvrho

c  rho_psi(psi)	       : tpsi,cpsirho

c  rmax_psi(psi)       : trmaxpsi(npsi4),crmaxpsi(npsi4),
c  rmin_psi(psi)       : trminsii(npsi4),crminpsi(npsi4),
c  bmax_psi(psi)       : tbmaxpsi(npsi4),cbmaxpsi(npsi4),
c  bmin_psi(psi)       : tbminpsi(npsi4),cbminpsi(npsi4),


c  rho_lpsi(psi)       : trholpol,crhlpsi
c  rho_lrho(rho)       : trhol_rh,crhlrho


c  rho_r(psi)          : trhor(npsi4),crhrpsi(npsi4)
c  (rho_r)**2(psi)     : trhor2(npsi4),crhrps2(npsi4)
c  rho_r(rho)          : trhor_rh(npsi4),crhrho(npsi4)

c to use spline subroutines
c coeff1 and terp1 from zcunix.f
c
c for rho_s:  d2_rhos_psi(npsi)
c for rho_bt: d2_rhotf_psi(npsi)
c for rho_v:  d2_rhovl_psi(npsi)
c for rho_l:  d2_rhol_psi(npsi)
c for rho_r:  d2_rhor_psi(npsi)

c  arrho_s(npsi), array rho from s (area)
c  arrho_tf(npsi), array rho from toroidal flux
c  arrho_vl(npsi), array rho from toroidal volume
c  arrho_r(npsi)   array rom from R major radius
cSAP081122 for B line length
c length_b_ar(npsi),d2_length_b_psi(npsi)
