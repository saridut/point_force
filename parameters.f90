 module parameters

 implicit none

 integer,parameter    :: dbl= selected_real_kind(15,300)
 integer,parameter    :: lng = selected_int_kind(15)

 real(dbl),parameter  :: pi = 3.141592653589793_dbl
 real(dbl),parameter  :: sqrt2 = 1.414213562373095_dbl
 real(dbl),parameter  :: sqrt3 = 1.732050807568877_dbl
 real(dbl), parameter :: third = 0.333333333333333_dbl

 integer, parameter    :: ndim = 3
 integer, parameter    :: nbead = 2
 real (dbl), parameter :: a = 0.5_dbl

 end module parameters

!*******************************************************************************
