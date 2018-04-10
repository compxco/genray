
c      the data for the profiles of the absorbed power and current drive.
c      PARAMETER (NRA) is set in param.i
c      It has namelist data only

c-----from namelist /tokamak/
      integer
     &NR

      common/onetwo_nml/
     &         NR

c---------------------------------------------------------------------------------
c    NR is the maximal number of points in the small radius direction 
c      for the
c      calculations of the power and current drive radial profiles 
c      in each radial bin.
c    NRA= maximal value of NR. It is set in param.i
c---------------------------------------------------------------------------------
