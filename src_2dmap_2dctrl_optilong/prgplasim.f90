!**************************************************************
!   Copyright (C) 2011
!   FastOpt GmbH, Hamburg, Germany
!
!   Contact: http://www.FastOpt.com
!
!   All rights reserved.
!**************************************************************

      program main
!*************************************************************
      implicit none
!=========================================
! declaration
!=========================================
      integer n, m

!==============================================================
! get the number of the independent and dependent variables
!==============================================================
      call setfunc( n, m )
!-----------------------------------------
! call the subroutine
!-----------------------------------------
      call prgfunction( n, m )

      end program main


      subroutine prgfunction( n, m )
        use pumamod, only : ncost
        use observation, only: postobs
!==============================================================
! declaration
!==============================================================
      implicit none
      
      integer n, m
      
      
      real  x(n)
      real  Y(m)
      integer nmax, mmax

      ncost = 0 ! do not evaluate cost function
!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================
      call initfunc( n, x )

!==============================================================
! call function to run the model
!==============================================================
      call func( n, x, m, y )

      call epilog
      call postobs
      call mpstop

      end subroutine prgfunction

