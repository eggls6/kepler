module kind_m

!----------------------------------------------------------------
!contains definition of precision for real and integer variables
!
! written by Siegfried EGGL 20170201
!----------------------------------------------------------------

      PUBLIC
      !KIND                                                        
      integer,parameter::ip=selected_int_kind(8)            !kind for integer 
      integer,parameter::wp=selected_real_kind(15,307)       !kind for real number precision
      
      !integer, parameter :: sp = selected_real_kind(6, 37)
      !integer, parameter :: dp = selected_real_kind(15, 307)
      !integer, parameter :: qp = selected_real_kind(33, 4931)
end module
