!-*-f90-*-

! This file is cpp-preprocessed. Be sure to edit func.F90

!#define DEBUG
#define STATIC_MEMTAPE

subroutine func(n,x,m,y)
! author: FastOpt GmbH
! purpose: top level routine for TAF, defining the function to differentiate
! input:
!      integer :: n    ! dimension of control vector
!      real    :: x(n) ! control vector
!      integer :: m    ! dimension of function value
! output:
!      real    :: y(m) ! function value  (result vector)
!
  use pumamod
  use oceanmod, only : NLEV_OCE
  use landmod, only : NLSOIL
  use mo_mapping
  use observation, only: reinitobs, nobs
!********************************************************************
!
! input  : N      amount of parameters
!          X      array of parameters
!
! output : 
!          Y      value of cost function or result vector
!
!********************************************************************
  IMPLICIT NONE

  INTEGER :: n, m
  REAL    :: x(n), y(m)
#if defined STATIC_MEMTAPE
!$TAF INIT func1 = STATIC, 1
#else
        include 'mpif.h' 
!$TAF INIT func1 = "tapes/func1"
#endif
#if defined STATIC_MEMTAPE
!$TAF INIT master1 = STATIC, (nstepe-nstep0+1+nkits)
!$TAF INIT mastercp = STATIC, (nstepe-nstep0+1+nkits)
!$TAF INIT gridpointa = STATIC, (nstepe-nstep0+1+nkits)
!$TAF INIT calcgp1 = STATIC, (nstepe-nstep0+1+nkits)
!$TAF INIT calcgpn = STATIC, (NLEV-2)*(nstepe-nstep0+1+nkits)
!$TAF INIT inner1 = STATIC, (nstepe-nstep0+1+nkits)
!$TAF INIT innern = STATIC, NLEV*(nstepe-nstep0+1+nkits)
!$TAF INIT inner2n = STATIC, (2*NLEV+1)*(nstepe-nstep0+1+nkits)
!$TAF INIT innernn2 = STATIC, ((NLEV+1)*NLEV/2)*(nstepe-nstep0+1+nkits)
!$TAF INIT inner2nn = STATIC, (2*NLEV+1)*(NLEM-1)*(nstepe-nstep0+1+nkits)
!$TAF INIT innerh = STATIC, NHOR*(nstepe-nstep0+1+nkits)
!$TAF INIT inners = STATIC, NLSOIL*(nstepe-nstep0+1+nkits)
!$TAF INIT innerno = STATIC, (NLEV_OCE-2)*(nstepe-nstep0+1+nkits)
!$TAF INIT obs1 = STATIC, 1+(nstepe-nstep0)/ncost
!$TAF INIT obsn = STATIC, nobs*(1+(nstepe-nstep0)/ncost)
#else
!$TAF INIT master1 = "tapes/master1"
!$TAF INIT mastercp = "tapes/mastercp"
!$TAF INIT gridpointa = "tapes/gridpointa"
!$TAF INIT calcgp1 = "tapes/calcgp1"
!$TAF INIT calcgpn = "tapes/calcgpn"
!$TAF INIT inner1 = "tapes/inner1"
!$TAF INIT innern = "tapes/innern"
!$TAF INIT inner2n = "tapes/inner2n"
!$TAF INIT innernn2 = "tapes/innernn2"
!$TAF INIT inner2nn = "tapes/inner2nn"
!$TAF INIT innerh = "tapes/innerh"
!$TAF INIT inners = "tapes/inners"
!$TAF INIT innerno = "tapes/innerno"
!$TAF INIT obs1 = "tapes/obs1"
!$TAF INIT obsn = "tapes/obsn"
#endif
  !
  ! restore initial conditions (created during initfunc)
  !
  call read_internal_restart('internal_restart')
  !
  ! set iteration count
  !
  nstep1=0
  nstep2=0
  !
  ! unpack arguments for model
  !
  call cont2mod(n,x)
  !
  ! initialise cost function
  !
  call reinitobs
  !
  ! run model
  !
  call master
!$TAF STORE nstep0,nstepe = func1, REC = 1
  !
  ! pack results of model
  !
  call mod2y(m,y)
end subroutine func


