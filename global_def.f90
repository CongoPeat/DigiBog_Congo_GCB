module global_definition
  implicit none
  !Global type definition.
  !Define a real kind type q with at least 8 decimal digits and an exponent
  !range from 10**30 to 10**(-30)
  integer, parameter :: q = SELECTED_REAL_KIND(P = 8, R = 30)
end module global_definition
