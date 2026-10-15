!**************************************************************
!   Copyright (C) 2008
!   FastOpt GmbH, Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
!
!    Purpose: Library for the replacement of non-differentiable
!             functions and constructs.
!
!    Usage: Refer to comments in interface declarations for usage.
!           All routines exist as subroutines and as functions.
!           If all arguments are scaler, the first (dimension) argument is optional.
!           The respective meaning of the tuning parameter is documented in the interface section
!           together with resulting maximal approximation error.
!
module substitutes
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! step functions !
!!!!!!!!!!!!!!!!!!
!
! logistic function: 
!
!    Regardless of parameter choice, the asymptotic approach to -1/1 is always faster than the atan version below.
!    This function is defined as 0.5*(1.+tanh(x/2))=1./(1.+exp(-x)).
!
!    tune is derivative at p==0
!    .75 reached at log(3)/4/tune ~= .2745/tune
!
! atan soft "step" function:
!
!    Slower asymptotic approch to -1/1 than logistic function. Probably faster convergence of grad.-based algor.
!
!    tune is derivative at p==0
!    .75 reached at p==1/pi/tune ~= .3183/tune
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! logistic function: result = logisticf ( integer :: n, real(n) :: argument, real :: tuning parameter ) 
  interface logisticf
     module procedure logisticf_s     ! real scalar
     module procedure logisticf_1     ! real 1d-array
  end interface

! logistic subroutine: call logistic ( integer :: n, real(n) :: argument, real :: tuning parameter, real(n) :: result ) 
  interface logistic
     module procedure logistic_s     ! real scalar
     module procedure logistic_1     ! real 1d-array
  end interface

! atan step function: result = atanstepf ( integer :: n, real(n) :: argument, real :: tuning parameter ) 
  interface atanstepf
     module procedure atanstepf_s     ! real scalar
     module procedure atanstepf_1     ! real 1d-array
  end interface

! atan step subroutine: call atanstep ( integer :: n, real(n) :: argument, real :: tuning parameter, real(n) :: result ) 
  interface atanstep
     module procedure atanstep_s     ! real scalar
     module procedure atanstep_1     ! real 1d-array
  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! min !
!!!!!!!
!
! min(a,b)-sqrt(tune)/2 <= fomin(a,b) <= min(a,b)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! differentiable min function: 
! result = fominf ( integer :: n, real(1|n) :: argument1, real(n|1) :: argument2, real :: tuning parameter [ , integer: order] ) 
  interface fominf
     module procedure fominf_ss     ! real scalar, scalar
     module procedure fominf_1s     ! real 1d-array, scalar
     module procedure fominf_s1     ! real scalar, 1d-array
     module procedure fominf_11     ! real 1d-array, 1d-array
     module procedure fominf_ss_o     ! real scalar, scalar
     module procedure fominf_1s_o     ! real 1d-array, scalar
     module procedure fominf_s1_o     ! real scalar, 1d-array
     module procedure fominf_11_o     ! real 1d-array, 1d-array
  end interface

! differentiable min subroutine: 
! call fomin ( integer :: n, real(1|n) :: argument1, real(n|1) :: argument2, real :: tuning parameter,
!                                                              [ , integer: order,]  real(n) :: result  )
  interface fomin
     module procedure fomin_ss     ! real scalar, scalar
     module procedure fomin_1s     ! real 1d-array, scalar
     module procedure fomin_s1     ! real scalar, 1d-array
     module procedure fomin_11     ! real 1d-array, 1d-array
     module procedure fomin_ss_o     ! real scalar, scalar
     module procedure fomin_1s_o     ! real 1d-array, scalar
     module procedure fomin_s1_o     ! real scalar, 1d-array
     module procedure fomin_11_o     ! real 1d-array, 1d-array
  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! max !
!!!!!!!
!
! max(a,b)+sqrt(tune)/2 >= fomax(a,b) >= max(a,b)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! differentiable max function: 
! result = fomaxf ( integer :: n, real(1|n) :: argument1, real(n|1) :: argument2, real :: tuning parameter [, integer: order] ) 
  interface fomaxf
     module procedure fomaxf_ss     ! real scalar, scalar
     module procedure fomaxf_1s     ! real 1d-array, scalar
     module procedure fomaxf_s1     ! real scalar, 1d-array
     module procedure fomaxf_11     ! real 1d-array, 1d-array
     module procedure fomaxf_ss_o     ! real scalar, scalar
     module procedure fomaxf_1s_o     ! real 1d-array, scalar
     module procedure fomaxf_s1_o     ! real scalar, 1d-array
     module procedure fomaxf_11_o     ! real 1d-array, 1d-array
  end interface

! differentiable max subroutine: 
! call fomax ( integer :: n, real(1|n) :: argument1, real(n|1) :: argument2, real :: tuning parameter, 
!                                                              [ integer: order,]  real(n) :: result  )
  interface fomax
     module procedure fomax_ss     ! real scalar, scalar
     module procedure fomax_1s     ! real 1d-array, scalar
     module procedure fomax_s1     ! real scalar, 1d-array
     module procedure fomax_11     ! real 1d-array, 1d-array
     module procedure fomax_ss_o     ! real scalar, scalar
     module procedure fomax_1s_o     ! real 1d-array, scalar
     module procedure fomax_s1_o     ! real scalar, 1d-array
     module procedure fomax_11_o     ! real 1d-array, 1d-array
  end interface

! differential max/min-routines. Use with constant as first argument, as in maxx(const.,x)
  interface minx
     module procedure minx_ss      ! real scalar, scalar
     module procedure minx_s1      ! real scalar, 1d-array
  end interface

  interface maxx
     module procedure maxx_ss      ! real scalar, scalar
     module procedure maxx_s1      ! real scalar,1d-array
  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! abs !
!!!!!!!
!
! abs(x) >= foabsf(x) >= abs(x) - .3... * sqrt(tune)
!
! maximum error = .3002831060007775 * sqrt( tune ) at x = sqrt( tune ) * sqrt( (sqrt(5)-1)/2 )
!
! or, in other words:
! (sqrt(sqrt(5) + 1)*(sqrt(2)* sqrt(5) - 3*sqrt(2)) + 2*sqrt(2) * sqrt(sqrt(5) - 1))/4 * sqrt(tune)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! differentiable abs function: 
! result = foabsf ( integer :: n, real(n) :: argument1, real :: tuning parameter [ , integer: order] ) 
  interface foabsf
     module procedure foabsf_s     ! real scalar
     module procedure foabsf_1     ! real 1d-array
! to be implemented     module procedure foabsf_s_o     ! real scalar
! to be implemented     module procedure foabsf_1_o     ! real 1d-array
  end interface

! differentiable abs subroutine: 
! call foabs ( integer :: n, real(n) :: argument1, real :: tuning parameter, [ integer: order,] real(n) :: result ) 
  interface foabs
     module procedure foabs_s     ! real scalar
     module procedure foabs_1     ! real 1d-array
! to be implemented     module procedure foabs_s_o     ! real scalar
! to be implemented     module procedure foabs_1_o     ! real 1d-array
  end interface

contains
!  function logisticf_s(ub,p,tune)
!    integer,optional :: ub
  function logisticf_s(p,tune)
    real :: p,logisticf_s
    real :: tune ! tune is derivative at p==0, .75 reached at log(3)/4/tune ~= .2745/tune
    logisticf_s = (1. + tanh( p * 2 * tune )) / 2
  end function logisticf_s

!  subroutine logistic_s(ub,p,tune,r)
!    integer,optional :: ub
  subroutine logistic_s(p,tune,r)
    real :: p,r
    real :: tune ! tune is derivative at p==0, .75 reached at log(3)/4/tune ~= .2745/tune
    r = (1. + tanh( p * 2 * tune )) / 2
  end subroutine logistic_s

  function logisticf_1(ub,p,tune)
    integer :: ub
    real :: p(ub)
    real :: logisticf_1(ub)
    real :: tune ! tune is derivative at p==0, .75 reached at log(3)/4/tune ~= .2745/tune
    logisticf_1(:) = (1. + tanh( p(:) * 2 * tune )) / 2
  end function logisticf_1

  subroutine logistic_1(ub,p,tune,r)
    integer :: ub
    real :: p(ub),r(ub)
    real :: tune ! tune is derivative at p==0, .75 reached at log(3)/4/tune ~= .2745/tune
    r(:) = (1. + tanh( p(:) * 2 * tune )) / 2
  end subroutine logistic_1

!  function atanstepf_s(ub,p,tune)
!    integer,optional :: ub
  function atanstepf_s(p,tune)
    real :: p,atanstepf_s
    real :: tune ! tune is derivative at p==0, .75 reached at p==1/pi/tune ~= .3183/tune
    real :: pi ! 3.1415...
    pi = 4*atan(1.) ! trusting the compiler to optimise this away
    atanstepf_s = 0.5 + atan( p * pi * tune ) / pi
  end function atanstepf_s

!  subroutine atanstep_s(ub,p,tune,r)
!    integer,optional :: ub
  subroutine atanstep_s(p,tune,r)
    real :: p,r
    real :: tune ! tune is derivative at p==0, .75 reached at p==1/pi/tune ~= .3183/tune
    real :: pi ! 3.1415...
    pi = 4*atan(1.) ! trusting the compiler to optimise this away
    r = 0.5 + atan( p * pi * tune ) / pi
  end subroutine atanstep_s

  function atanstepf_1(ub,p,tune)
    integer :: ub
    real :: p(ub),atanstepf_1(ub)
    real :: tune ! tune is derivative at p==0, .75 reached at p==1/pi/tune ~= .3183/tune
    real :: pi ! 3.1415...
    pi = 4*atan(1.) ! trusting the compiler to optimise this away
    atanstepf_1(:) = 0.5 + atan( p(:) * pi * tune ) / pi
  end function atanstepf_1

  subroutine atanstep_1(ub,p,tune,r)
    integer :: ub
    real :: p(ub),r(ub)
    real :: tune ! tune is derivative at p==0, .75 reached at p==1/pi/tune ~= .3183/tune
    real :: pi ! 3.1415...
    pi = 4*atan(1.) ! trusting the compiler to optimise this away
    r(:) = 0.5 + atan( p(:) * pi * tune ) / pi
  end subroutine atanstep_1

!  function fominf_ss(ub,p1,p2,tune)
!    integer,optional :: ub
  function fominf_ss(p1,p2,tune)
    real :: p1,p2,fominf_ss
    real :: tune ! sqrt(tune)/2 >= fomin - min(p1,p2)
    fominf_ss = sqrt( (p1-p2)**2 + tune )
    fominf_ss = ( p1 + p2 - fominf_ss ) / 2
  end function fominf_ss

!  subroutine fomin_ss(ub,p1,p2,tune,r)
!    integer,optional :: ub
  subroutine fomin_ss(p1,p2,tune,r)
    real :: p1,p2,r
    real :: tune ! sqrt(tune)/2 >= fomin - min(p1,p2)
    r = sqrt( (p1-p2)**2 + tune )
    r = ( p1 + p2 - r ) / 2
  end subroutine fomin_ss

  REAL FUNCTION minx_ss (x, y, x0)
    REAL :: x, y, x0
    IF (x.LE.y+x0) THEN
       minx_ss = x - x0*EXP((x-y)/x0-1.)
    ELSE
       minx_ss = y
    ENDIF
  END FUNCTION minx_ss

  FUNCTION minx_s1 (ub, x, y, x0)
    integer :: ub
    REAL :: x, x0
    real, dimension(ub) :: y,minx_s1
    where (x.LE.y+x0)
       minx_s1 = x - x0*EXP((x-y)/x0-1.)
    ELSEwhere
       minx_s1 = y
    END where
  END FUNCTION minx_s1

  REAL FUNCTION maxx_ss (x, y, x0)
    REAL :: x, y, x0
    IF (x.GE.y-x0) THEN
       maxx_ss = x + x0*EXP(-(x-y)/x0-1.)
    ELSE
       maxx_ss = y
    ENDIF
  END FUNCTION maxx_ss

  FUNCTION maxx_s1 (ub, x, y, x0)
    integer :: ub
    REAL :: x, x0
    real, dimension(ub) :: y, maxx_s1
    where (x.GE.y-x0)
       maxx_s1 = x + x0*EXP(-(x-y)/x0-1.)
    ELSEwhere
       maxx_s1 = y
    END where
  END FUNCTION maxx_s1

  function fominf_s1(ub,p1s,p2,tune)
    integer :: ub
    real :: p2(ub),fominf_s1(ub)
    real :: tune, p1s  ! sqrt(tune)/2 >= fomin - min(p1s,p2)
    fominf_s1(:) = sqrt( (p1s-p2(:))**2 + tune )
    fominf_s1(:) = ( p1s + p2(:) - fominf_s1(:) ) / 2
  end function fominf_s1

  subroutine fomin_s1(ub,p1s,p2,tune,r)
    integer :: ub
    real :: p2(ub),r(ub)
    real :: tune, p1s  ! sqrt(tune)/2 >= fomin - min(p1s,p2)
    r(:) = sqrt( (p1s-p2(:))**2 + tune )
    r(:) = ( p1s + p2(:) - r(:) ) / 2
  end subroutine fomin_s1

  function fominf_1s(ub,p1,p2s,tune)
    integer :: ub
    real :: p1(ub),fominf_1s(ub)
    real :: tune, p2s ! sqrt(tune)/2 >= fomin - min(p1,p2s)
    fominf_1s(:) = sqrt( (p1(:)-p2s)**2 + tune )
    fominf_1s(:) = ( p1(:) + p2s - fominf_1s(:) ) / 2
  end function fominf_1s

  subroutine fomin_1s(ub,p1,p2s,tune,r)
    integer :: ub
    real :: p1(ub),r(ub)
    real :: tune, p2s ! sqrt(tune)/2 >= fomin - min(p1,p2s)
    r(:) = sqrt( (p1(:)-p2s)**2 + tune )
    r(:) = ( p1(:) + p2s - r(:) ) / 2
  end subroutine fomin_1s

  function fominf_11(ub,p1,p2,tune)
    integer :: ub
    real :: p1(ub),p2(ub),fominf_11(ub)
    real :: tune ! sqrt(tune)/2 >= fomin - min(p1,p2)
    fominf_11(:) = sqrt( (p1(:)-p2(:))**2 + tune )
    fominf_11(:) = ( p1(:) + p2(:) - fominf_11(:) ) / 2
  end function fominf_11

  subroutine fomin_11(ub,p1,p2,tune,r)
    integer :: ub
    real :: p1(ub),p2(ub),r(ub)
    real :: tune ! sqrt(tune)/2 >= fomin - min(p1,p2)
    r(:) = sqrt( (p1(:)-p2(:))**2 + tune )
    r(:) = ( p1(:) + p2(:) - r(:) ) / 2
  end subroutine fomin_11

!  function fominf_ss(ub,p1,p2,tune)
!    integer,optional :: ub
  function fominf_ss_o(p1,p2,tune,ko)
    real :: p1,p2,fominf_ss_o
    real :: tune ! sqrt(tune)/2 >= fomin - min(p1,p2)
    integer :: ko
    fominf_ss_o = ( (p1-p2)**(ko*2) + tune**ko )
    fominf_ss_o = ( p1 + p2 - fominf_ss_o ) / 2
  end function fominf_ss_o

!  subroutine fomin_ss(ub,p1,p2,tune,ko,r)
!    integer,optional :: ub
  subroutine fomin_ss_o(p1,p2,tune,ko,r)
    real :: p1,p2,r
    real :: tune ! sqrt(tune)/2 >= fomin - min(p1,p2)
    integer :: ko
    r = sqrt( (p1-p2)**(ko*2)+ tune**ko )
    r = ( p1 + p2 - r ) / 2
  end subroutine fomin_ss_o

  function fominf_s1_o(ub,p1s,p2,tune,ko)
    integer :: ub
    real :: p2(ub),fominf_s1_o(ub)
    real :: tune, p1s  ! sqrt(tune)/2 >= fomin - min(p1s,p2)
    integer :: ko
    fominf_s1_o(:) = sqrt( (p1s-p2(:))**(ko*2)+ tune**ko )
    fominf_s1_o(:) = ( p1s + p2(:) - fominf_s1_o(:) ) / 2
  end function fominf_s1_o

  subroutine fomin_s1_o(ub,p1s,p2,tune,ko,r)
    integer :: ub
    real :: p2(ub),r(ub)
    real :: tune, p1s  ! sqrt(tune)/2 >= fomin - min(p1s,p2)
    integer :: ko
    r(:) = sqrt( (p1s-p2(:))**(ko*2)+ tune**ko )
    r(:) = ( p1s + p2(:) - r(:) ) / 2
  end subroutine fomin_s1_o

  function fominf_1s_o(ub,p1,p2s,tune,ko)
    integer :: ub
    real :: p1(ub),fominf_1s_o(ub)
    real :: tune, p2s ! sqrt(tune)/2 >= fomin - min(p1,p2s)
    integer :: ko
    fominf_1s_o(:) = sqrt( (p1(:)-p2s)**(ko*2)+ tune**ko )
    fominf_1s_o(:) = ( p1(:) + p2s - fominf_1s_o(:) ) / 2
  end function fominf_1s_o

  subroutine fomin_1s_o(ub,p1,p2s,tune,ko,r)
    integer :: ub
    real :: p1(ub),r(ub)
    real :: tune, p2s ! sqrt(tune)/2 >= fomin - min(p1,p2s)
    integer :: ko
    r(:) = sqrt( (p1(:)-p2s)**(ko*2)+ tune**ko )
    r(:) = ( p1(:) + p2s - r(:) ) / 2
  end subroutine fomin_1s_o

  function fominf_11_o(ub,p1,p2,tune,ko)
    integer :: ub
    real :: p1(ub),p2(ub),fominf_11_o(ub)
    real :: tune ! sqrt(tune)/2 >= fomin - min(p1,p2)
    integer :: ko
    fominf_11_o(:) = sqrt( (p1(:)-p2(:))**(ko*2)+ tune**ko )
    fominf_11_o(:) = ( p1(:) + p2(:) - fominf_11_o(:) ) / 2
  end function fominf_11_o

  subroutine fomin_11_o(ub,p1,p2,tune,ko,r)
    integer :: ub
    real :: p1(ub),p2(ub),r(ub)
    real :: tune ! sqrt(tune)/2 >= fomin - min(p1,p2)
    integer :: ko
    r(:) = sqrt( (p1(:)-p2(:))**(ko*2)+ tune**ko )
    r(:) = ( p1(:) + p2(:) - r(:) ) / 2
  end subroutine fomin_11_o

!  function fomaxf_ss(ub,p1,p2,tune)
!    integer,optional :: ub
  function fomaxf_ss(p1,p2,tune)
    real :: p1,p2,fomaxf_ss
    real :: tune ! sqrt(tune)/2 >= max(p1,p2) - fomax
    fomaxf_ss = sqrt( (p1-p2)**2 + tune )
    fomaxf_ss = ( p1 + p2 + fomaxf_ss ) / 2
  end function fomaxf_ss

!  subroutine fomax_ss(ub,p1,p2,tune,r)
!    integer,optional :: ub
  subroutine fomax_ss(p1,p2,tune,r)
    real :: p1,p2,r
    real :: tune ! sqrt(tune)/2 >= max(p1,p2) - fomax
    r = sqrt( (p1-p2)**2 + tune )
    r = ( p1 + p2 + r ) / 2
  end subroutine fomax_ss

  function fomaxf_s1(ub,p1s,p2,tune)
    integer :: ub
    real :: p2(ub),fomaxf_s1(ub)
    real :: tune, p1s ! sqrt(tune)/2 >= max(p1s,p2) - fomax
    fomaxf_s1(:) = sqrt( (p1s-p2(:))**2 + tune )
    fomaxf_s1(:) = ( p1s + p2(:) + fomaxf_s1(:) ) / 2
  end function fomaxf_s1

  subroutine fomax_s1(ub,p1s,p2,tune,r)
    integer :: ub
    real :: p2(ub),r(ub)
    real :: tune, p1s ! sqrt(tune)/2 >= max(p1s,p2) - fomax
    r(:) = sqrt( (p1s-p2(:))**2 + tune )
    r(:) = ( p1s + p2(:) + r(:) ) / 2
  end subroutine fomax_s1

  function fomaxf_1s(ub,p1,p2s,tune)
    integer :: ub
    real :: p1(ub),fomaxf_1s(ub)
    real :: tune, p2s ! sqrt(tune)/2 >= max(p1,p2s) - fomax
    fomaxf_1s(:) = sqrt( (p1(:)-p2s)**2 + tune )
    fomaxf_1s(:) = ( p1(:) + p2s + fomaxf_1s(:) ) / 2
  end function fomaxf_1s

  subroutine fomax_1s(ub,p1,p2s,tune,r)
    integer :: ub
    real :: p1(ub),r(ub)
    real :: tune, p2s ! sqrt(tune)/2 >= max(p1,p2s) - fomax
    r(:) = sqrt( (p1(:)-p2s)**2 + tune )
    r(:) = ( p1(:) + p2s + r(:) ) / 2
  end subroutine fomax_1s

  function fomaxf_11(ub,p1,p2,tune)
    integer :: ub
    real :: p1(ub),p2(ub),fomaxf_11(ub)
    real :: tune ! sqrt(tune)/2 >= max(p1,p2) - fomax
    fomaxf_11(:) = sqrt( (p1(:)-p2(:))**2 + tune )
    fomaxf_11(:) = ( p1(:) + p2(:) + fomaxf_11(:) ) / 2
  end function fomaxf_11

  subroutine fomax_11(ub,p1,p2,tune,r)
    integer :: ub
    real :: p1(ub),p2(ub),r(ub)
    real :: tune ! sqrt(tune/2) >= max(p1,p2) - fomax
    r(:) = sqrt( (p1(:)-p2(:))**2 + tune )
    r(:) = ( p1(:) + p2(:) + r(:) ) / 2
  end subroutine fomax_11

!  function fomaxf_ss(ub,p1,p2,tune)
!    integer,optional :: ub
  function fomaxf_ss_o(p1,p2,tune,ko)
    real :: p1,p2,fomaxf_ss_o
    real :: tune ! sqrt(tune)/2 >= fomax - min(p1,p2)
    integer :: ko
    fomaxf_ss_o = ( (p1-p2)**(ko*2) + tune**ko )
    fomaxf_ss_o = ( p1 + p2 + fomaxf_ss_o ) / 2
  end function fomaxf_ss_o

!  subroutine fomax_ss(ub,p1,p2,tune,ko,r)
!    integer,optional :: ub
  subroutine fomax_ss_o(p1,p2,tune,ko,r)
    real :: p1,p2,r
    real :: tune ! sqrt(tune)/2 >= fomax - min(p1,p2)
    integer :: ko
    r = sqrt( (p1-p2)**(ko*2)+ tune**ko )
    r = ( p1 + p2 + r ) / 2
  end subroutine fomax_ss_o

  function fomaxf_s1_o(ub,p1s,p2,tune,ko)
    integer :: ub
    real :: p2(ub),fomaxf_s1_o(ub)
    real :: tune, p1s  ! sqrt(tune)/2 >= fomax - min(p1s,p2)
    integer :: ko
    fomaxf_s1_o(:) = sqrt( (p1s-p2(:))**(ko*2)+ tune**ko )
    fomaxf_s1_o(:) = ( p1s + p2(:) + fomaxf_s1_o(:) ) / 2
  end function fomaxf_s1_o

  subroutine fomax_s1_o(ub,p1s,p2,tune,ko,r)
    integer :: ub
    real :: p2(ub),r(ub)
    real :: tune, p1s  ! sqrt(tune)/2 >= fomax - min(p1s,p2)
    integer :: ko
    r(:) = sqrt( (p1s-p2(:))**(ko*2)+ tune**ko )
    r(:) = ( p1s + p2(:) + r(:) ) / 2
  end subroutine fomax_s1_o

  function fomaxf_1s_o(ub,p1,p2s,tune,ko)
    integer :: ub
    real :: p1(ub),fomaxf_1s_o(ub)
    real :: tune, p2s ! sqrt(tune)/2 >= fomax - min(p1,p2s)
    integer :: ko
    fomaxf_1s_o(:) = sqrt( (p1(:)-p2s)**(ko*2)+ tune**ko )
    fomaxf_1s_o(:) = ( p1(:) + p2s + fomaxf_1s_o(:) ) / 2
  end function fomaxf_1s_o

  subroutine fomax_1s_o(ub,p1,p2s,tune,ko,r)
    integer :: ub
    real :: p1(ub),r(ub)
    real :: tune, p2s ! sqrt(tune)/2 >= fomax - min(p1,p2s)
    integer :: ko
    r(:) = sqrt( (p1(:)-p2s)**(ko*2)+ tune**ko )
    r(:) = ( p1(:) + p2s + r(:) ) / 2
  end subroutine fomax_1s_o

  function fomaxf_11_o(ub,p1,p2,tune,ko)
    integer :: ub
    real :: p1(ub),p2(ub),fomaxf_11_o(ub)
    real :: tune ! sqrt(tune)/2 >= fomax - min(p1,p2)
    integer :: ko
    fomaxf_11_o(:) = sqrt( (p1(:)-p2(:))**(ko*2)+ tune**ko )
    fomaxf_11_o(:) = ( p1(:) + p2(:) + fomaxf_11_o(:) ) / 2
  end function fomaxf_11_o

  subroutine fomax_11_o(ub,p1,p2,tune,ko,r)
    integer :: ub
    real :: p1(ub),p2(ub),r(ub)
    real :: tune ! sqrt(tune)/2 >= fomax - min(p1,p2)
    integer :: ko
    r(:) = sqrt( (p1(:)-p2(:))**(ko*2)+ tune**ko )
    r(:) = ( p1(:) + p2(:) + r(:) ) / 2
  end subroutine fomax_11_o

!  function foabsf_s(ub,p,tune)
!    integer,optional :: ub
  function foabsf_s(p,tune)
    real :: p,foabsf_s
    real :: tune ! second deriv. and curvature at p==0: 2/sqrt(tune)
    foabsf_s = sqrt( p**2 + tune )
    foabsf_s = foabsf_s - tune / foabsf_s
  end function foabsf_s

!  subroutine foabs_s(ub,p,tune,r)
!    integer,optional :: ub
  subroutine foabs_s(p,tune,r)
    real :: p,r
    real :: tune ! second deriv. and curvature at p==0: 2/sqrt(tune)
    r = sqrt( p**2 + tune )
    r = r - tune / r
  end subroutine foabs_s

  function foabsf_1(ub,p,tune)
    integer :: ub
    real :: p(ub),foabsf_1(ub)
    real :: tune ! second deriv. and curvature at p==0: 2/sqrt(tune)
    foabsf_1 = sqrt( p(:)**2 + tune )
    foabsf_1(:) = foabsf_1(:) - tune / foabsf_1(:)
  end function foabsf_1

  subroutine foabs_1(ub,p,tune,r)
    integer :: ub
    real :: p(ub),r(ub)
    real :: tune ! second deriv. and curvature at p==0: 2/sqrt(tune)
    r = sqrt( p(:)**2 + tune )
    r(:) = r(:) - tune / r(:)
  end subroutine foabs_1

end module substitutes
