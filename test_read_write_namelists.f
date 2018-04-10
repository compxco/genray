       program test
c------reads genray.dat using read_all_namelists
c------writes new genray.dat using write_all_namelists
       call  prepare_genray_input
       stop
       end

      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
c      write(*,*) nf_strerror(iret)
         write(*,*) 'netCDF error'
      stop 'check_err:'
      endif
      end
