
c     Data from  namelists for common /ions_nml/
c
c     the charges and masses of the plasma species.
c     parameter (nbulka)
c     nbulka>=nbulk

c-----from namelist /species/
      real*8
     &charge,dmas

      common/ions_nml/ charge(nbulka),dmas(nbulka)

