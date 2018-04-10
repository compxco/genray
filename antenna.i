c     it is to specify the possible curves for EC cone vertex positions

      real*8 r_antenna_max,delta_z
      integer i_cone,i_ant
      common /antenna/ r_antenna_max, ! for i_ant=1 type of antenna curve
     &delta_z, ! vertical shift over (zmax,rzmax) point
     &i_cone,  !for i_ant=2 type of antenna curve i_ant=2
     &i_ant


c     r_antenna_max, ! for i_ant=1 type of antenna curve
c     &delta_z, ! [m] vertical shift over (zmax,rzmax) point
c     &i_cone,  !for i_ant=2 type of antenna curve i_ant=2
c     &i_ant    ! sets the form of the curve
c                 =1 circuler curve
c                 =2 ellipse curve 
