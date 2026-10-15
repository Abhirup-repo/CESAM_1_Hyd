!**************************************************************
!   Copyright (C) 2008
!   FastOpt GmbH, Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
!
!    Purpose: Development aid for substmod.f90
!
!    Usage: Link with this module instead of substmod.f90 in the replacement
!           stage of your code development in order to check whether your changed
!           are consistent with your original code. substmod.f90 converges towards 
!           this module for extreme choices of tuning parameters.
!
module substitutes
  implicit none

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

  interface minx
     module procedure minx_ss      ! real scalar, scalar
     module procedure minx_s1      ! real scalar, 1d-array
  end interface

  interface maxx
     module procedure maxx_ss      ! real scalar, scalar
     module procedure maxx_s1      ! real scalar,1d-array
  end interface

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
    real :: tune
    if (p.gt.0) then
       logisticf_s = 1.
    else
       logisticf_s = 0.
    end if
  end function logisticf_s

!  subroutine logistic_s(ub,p,tune,r)
!    integer,optional :: ub
  subroutine logistic_s(p,tune,r)
    real :: p,r
    real :: tune
    if (p.gt.0) then
       r = 1.
    else
       r = 0.
    end if
  end subroutine logistic_s

  function logisticf_1(ub,p,tune)
    integer :: ub
    real :: p(ub)
    real :: logisticf_1(ub)
    real :: tune
    where (p(:).gt.0)
       logisticf_1 = 1.
    elsewhere
       logisticf_1 = 0.
    end where
  end function logisticf_1

  subroutine logistic_1(ub,p,tune,r)
    integer :: ub
    real :: p(ub),r(ub)
    real :: tune
    where (p(:).gt.0)
       r(:) = 1.
    elsewhere
       r(:) = 0.
    end where
  end subroutine logistic_1

!  function atanstepf_s(ub,p,tune)
!    integer,optional :: ub
  function atanstepf_s(p,tune)
    real :: p,atanstepf_s
    real :: tune
    if (p.gt.0) then
       atanstepf_s = 1.
    else
       atanstepf_s = 0.
    end if
  end function atanstepf_s

!  subroutine atanstep_s(ub,p,tune,r)
!    integer,optional :: ub
  subroutine atanstep_s(p,tune,r)
    real :: p,r
    real :: tune
    if (p.gt.0) then
       r = 1.
    else
       r = 0.
    end if
  end subroutine atanstep_s

  function atanstepf_1(ub,p,tune)
    integer :: ub
    real :: p(ub),atanstepf_1(ub)
    real :: tune
    where (p(:).gt.0)
       atanstepf_1(:) = 1.
    elsewhere
       atanstepf_1(:) = 0.
    end where
  end function atanstepf_1

  subroutine atanstep_1(ub,p,tune,r)
    integer :: ub
    real :: p(ub),r(ub)
    real :: tune
    where (p(:).gt.0)
       r(:) = 1.
    elsewhere
       r(:) = 0.
    end where
  end subroutine atanstep_1

!  function fominf_ss(ub,p1,p2,tune)
!    integer,optional :: ub
  function fominf_ss(p1,p2,tune)
    real :: p1,p2,fominf_ss
    real :: tune
    fominf_ss = min(p1,p2)
  end function fominf_ss

!  subroutine fomin_ss(ub,p1,p2,tune,r)
!    integer,optional :: ub
  subroutine fomin_ss(p1,p2,tune,r)
    real :: p1,p2,r
    real :: tune
    r = min(p1,p2)
  end subroutine fomin_ss

  function fominf_s1(ub,p1s,p2,tune)
    integer :: ub
    real :: p2(ub),fominf_s1(ub)
    real :: tune, p1s
    fominf_s1(:) = min(p1s,p2(:))
  end function fominf_s1

  subroutine fomin_s1(ub,p1s,p2,tune,r)
    integer :: ub
    real :: p2(ub),r(ub)
    real :: tune, p1s
    r(:) = min(p1s,p2(:))
  end subroutine fomin_s1

  REAL FUNCTION minx_ss (x, y, x0)
    REAL :: x, y, x0
    minx_ss = min(x,y)
  END FUNCTION minx_ss

  FUNCTION minx_s1 (ub, x, y, x0)
    integer :: ub
    REAL :: x, x0
    real, dimension(ub) :: y,minx_s1
    minx_s1 = min(x,y)
  END FUNCTION minx_s1

  REAL FUNCTION maxx_ss (x, y, x0)
    REAL :: x, y, x0
    maxx_ss = max(x,y)
  END FUNCTION maxx_ss

  FUNCTION maxx_s1 (ub, x, y, x0)
    integer :: ub
    REAL :: x, x0
    real, dimension(ub) :: y, maxx_s1
    maxx_s1 = max(x,y)
  END FUNCTION maxx_s1

  function fominf_1s(ub,p1,p2s,tune)
    integer :: ub
    real :: p1(ub),fominf_1s(ub)
    real :: tune, p2s
    fominf_1s(:) = min(p1(:),p2s)
  end function fominf_1s

  subroutine fomin_1s(ub,p1,p2s,tune,r)
    integer :: ub
    real :: p1(ub),r(ub)
    real :: tune, p2s
    r(:) = min(p1(:),p2s)
  end subroutine fomin_1s

  function fominf_11(ub,p1,p2,tune)
    integer :: ub
    real :: p1(ub),p2(ub),fominf_11(ub)
    real :: tune
    fominf_11(:) = min(p1(:),p2(:))
  end function fominf_11

  subroutine fomin_11(ub,p1,p2,tune,r)
    integer :: ub
    real :: p1(ub),p2(ub),r(ub)
    real :: tune
    r(:) = min(p1(:),p2(:))
  end subroutine fomin_11

!  function fominf_ss(ub,p1,p2,tune)
!    integer,optional :: ub
  function fominf_ss_o(p1,p2,tune,ko)
    real :: p1,p2,fominf_ss_o
    real :: tune
    integer :: ko
    fominf_ss_o = min(p1,p2)
  end function fominf_ss_o

!  subroutine fomin_ss(ub,p1,p2,tune,ko,r)
!    integer,optional :: ub
  subroutine fomin_ss_o(p1,p2,tune,ko,r)
    real :: p1,p2,r
    real :: tune
    integer :: ko
    r = min(p1,p2)
  end subroutine fomin_ss_o

  function fominf_s1_o(ub,p1s,p2,tune,ko)
    integer :: ub
    real :: p2(ub),fominf_s1_o(ub)
    real :: tune, p1s
    integer :: ko
    fominf_s1_o(:) = min(p1s,p2(:))
  end function fominf_s1_o

  subroutine fomin_s1_o(ub,p1s,p2,tune,ko,r)
    integer :: ub
    real :: p2(ub),r(ub)
    real :: tune, p1s
    integer :: ko
    r(:) = min(p1s,p2(:))
  end subroutine fomin_s1_o

  function fominf_1s_o(ub,p1,p2s,tune,ko)
    integer :: ub
    real :: p1(ub),fominf_1s_o(ub)
    real :: tune, p2s
    integer :: ko
    fominf_1s_o(:) = min(p1(:),p2s)
  end function fominf_1s_o

  subroutine fomin_1s_o(ub,p1,p2s,tune,ko,r)
    integer :: ub
    real :: p1(ub),r(ub)
    real :: tune, p2s
    integer :: ko
    r(:) = min(p1(:),p2s)
  end subroutine fomin_1s_o

  function fominf_11_o(ub,p1,p2,tune,ko)
    integer :: ub
    real :: p1(ub),p2(ub),fominf_11_o(ub)
    real :: tune
    integer :: ko
    fominf_11_o(:) = min(p1(:),p2(:))
  end function fominf_11_o

  subroutine fomin_11_o(ub,p1,p2,tune,ko,r)
    integer :: ub
    real :: p1(ub),p2(ub),r(ub)
    real :: tune
    integer :: ko
    r(:) = min(p1(:),p2(:))
  end subroutine fomin_11_o

!  function fomaxf_ss(ub,p1,p2,tune)
!    integer,optional :: ub
  function fomaxf_ss(p1,p2,tune)
    real :: p1,p2,fomaxf_ss
    real :: tune
    fomaxf_ss = max(p1,p2)
  end function fomaxf_ss

!  subroutine fomax_ss(ub,p1,p2,tune,r)
!    integer,optional :: ub
  subroutine fomax_ss(p1,p2,tune,r)
    real :: p1,p2,r
    real :: tune
    r = max(p1,p2)
  end subroutine fomax_ss

  function fomaxf_s1(ub,p1s,p2,tune)
    integer :: ub
    real :: p2(ub),fomaxf_s1(ub)
    real :: tune, p1s
    fomaxf_s1(:) = max(p1s,p2(:))
  end function fomaxf_s1

  subroutine fomax_s1(ub,p1s,p2,tune,r)
    integer :: ub
    real :: p2(ub),r(ub)
    real :: tune, p1s
    r(:) = max(p1s,p2(:))
  end subroutine fomax_s1

  function fomaxf_1s(ub,p1,p2s,tune)
    integer :: ub
    real :: p1(ub),fomaxf_1s(ub)
    real :: tune, p2s
    fomaxf_1s(:) = max(p1(:),p2s)
  end function fomaxf_1s

  subroutine fomax_1s(ub,p1,p2s,tune,r)
    integer :: ub
    real :: p1(ub),r(ub)
    real :: tune, p2s
    r(:) = max(p1(:),p2s)
  end subroutine fomax_1s

  function fomaxf_11(ub,p1,p2,tune)
    integer :: ub
    real :: p1(ub),p2(ub),fomaxf_11(ub)
    real :: tune
    fomaxf_11(:) = max(p1(:),p2(:))
  end function fomaxf_11

  subroutine fomax_11(ub,p1,p2,tune,r)
    integer :: ub
    real :: p1(ub),p2(ub),r(ub)
    real :: tune
    r(:) = max(p1(:),p2(:))
  end subroutine fomax_11

!  function fomaxf_ss(ub,p1,p2,tune)
!    integer,optional :: ub
  function fomaxf_ss_o(p1,p2,tune,ko)
    real :: p1,p2,fomaxf_ss_o
    real :: tune
    integer :: ko
    fomaxf_ss_o = max(p1,p2)
  end function fomaxf_ss_o

!  subroutine fomax_ss(ub,p1,p2,tune,ko,r)
!    integer,optional :: ub
  subroutine fomax_ss_o(p1,p2,tune,ko,r)
    real :: p1,p2,r
    real :: tune
    integer :: ko
    r = max(p1,p2)
  end subroutine fomax_ss_o

  function fomaxf_s1_o(ub,p1s,p2,tune,ko)
    integer :: ub
    real :: p2(ub),fomaxf_s1_o(ub)
    real :: tune, p1s
    integer :: ko
    fomaxf_s1_o(:) = max(p1s,p2(:))
  end function fomaxf_s1_o

  subroutine fomax_s1_o(ub,p1s,p2,tune,ko,r)
    integer :: ub
    real :: p2(ub),r(ub)
    real :: tune, p1s
    integer :: ko
    r(:) = max(p1s,p2(:))
  end subroutine fomax_s1_o

  function fomaxf_1s_o(ub,p1,p2s,tune,ko)
    integer :: ub
    real :: p1(ub),fomaxf_1s_o(ub)
    real :: tune, p2s
    integer :: ko
    fomaxf_1s_o(:) = max(p1(:),p2s)
  end function fomaxf_1s_o

  subroutine fomax_1s_o(ub,p1,p2s,tune,ko,r)
    integer :: ub
    real :: p1(ub),r(ub)
    real :: tune, p2s
    integer :: ko
    r(:) = max(p1(:),p2s)
  end subroutine fomax_1s_o

  function fomaxf_11_o(ub,p1,p2,tune,ko)
    integer :: ub
    real :: p1(ub),p2(ub),fomaxf_11_o(ub)
    real :: tune
    integer :: ko
    fomaxf_11_o(:) = max(p1(:),p2(:))
  end function fomaxf_11_o

  subroutine fomax_11_o(ub,p1,p2,tune,ko,r)
    integer :: ub
    real :: p1(ub),p2(ub),r(ub)
    real :: tune
    integer :: ko
    r(:) = max(p1(:),p2(:))
  end subroutine fomax_11_o

!  function foabsf_s(ub,p,tune)
!    integer,optional :: ub
  function foabsf_s(p,tune)
    real :: p,foabsf_s
    real :: tune
    foabsf_s = abs(p)
  end function foabsf_s

!  subroutine foabs_s(ub,p,tune,r)
!    integer,optional :: ub
  subroutine foabs_s(p,tune,r)
    real :: p,r
    real :: tune
    r = abs(p)
  end subroutine foabs_s

  function foabsf_1(ub,p,tune)
    integer :: ub
    real :: p(ub),foabsf_1(ub)
    real :: tune
    foabsf_1 = abs(p(:))
  end function foabsf_1

  subroutine foabs_1(ub,p,tune,r)
    integer :: ub
    real :: p(ub),r(ub)
    real :: tune
    r = abs(p(:))
  end subroutine foabs_1

end module substitutes
