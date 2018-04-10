      module kind_spec
      integer, parameter :: rprec = selected_real_kind(12,100)
      integer, parameter :: iprec = selected_int_kind(8)
      integer, parameter :: dp = rprec
      real(kind=dp), parameter :: d1mach1 = epsilon(1.0_dp)
      real(kind=dp), parameter :: d1mach2 = huge(1.0_dp)
      real(kind=dp),parameter::d1mach3=epsilon(1.0_dp)/radix(1.0_dp)
      real(kind=dp), parameter :: d1mach4 = epsilon(1.0_dp)
      real(kind=dp), parameter :: d1machx = radix(1.0_dp)
! --- Note d1mach5 should be log10( FLT_RADIX) but we can not
!---  call the log function in an initialization
      real(kind=dp) :: d1mach5 =  0.301029995663981
      integer :: j1mach15 = minexponent(1.0_dp)
      integer :: j1mach14 = digits(1.0_dp)
      end module kind_spec
