c    The data for  subroutine in oxb.f 
      logical was_in_ox_vicinity
      real*8  transm_ox_old
      integer nrayelt_o_cutoff,
     &i_call_prep3d_in_output_at_i_ox_conversion_eq_1
      common /oxb/was_in_ox_vicinity,transm_ox_old,
     & nrayelt_o_cutoff,
     &i_call_prep3d_in_output_at_i_ox_conversion_eq_1
c     was_in_ox_vicinity !for OX jump 
c     transm_ox_old the previous value of the transmission coefficient
c
c     i_call_prep3d_in_output_at_i_ox_conversion_eq_1 - It is used
c     in prep3d.f to calculate the transmission coefficient 
c     when prep3d was called in output after OX conversion point
c     determination   
c
c      nrayelt_o_cutoff is the number of ray point where O cutoff was
c                       found
