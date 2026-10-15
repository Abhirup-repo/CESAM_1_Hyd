! Dummy module replacement for guimod.f90
! Use this for environments with no X11

subroutine guistart
return
end subroutine guistart

subroutine guistep_puma
return
end subroutine guistep_puma

subroutine guistep_plasim
return
end subroutine guistep_plasim

subroutine guistop
return
end subroutine guistop

subroutine guihor(yn,f,klev,pm,pa) ! changed, FastOpt
  use pumamod
  implicit none ! FastOpt
  character (len=*) :: yn ! FastOpt
  integer :: klev ! FastOpt
  real :: f(NLON,NLPP,klev) ! FastOpt
  real :: pm, pa ! FastOpt
  return
end subroutine guihor

subroutine guigv(yn,f)
  use pumamod  ! FastOpt
  implicit none ! FastOpt
  character (len=*) :: yn
  real :: f(NLON,NLPP,NLEV) ! FastOpt
  return
end subroutine guigv

subroutine guips(fp)
  implicit none ! FastOpt
  real :: fp(*) ! FastOpt
  return
end subroutine guips

subroutine guiput(yn,f,k1,k2,k3)
  implicit none ! FastOpt
  character (len=*) :: yn ! FastOpt
  integer :: k1,k2,k3 ! FastOpt
  real (kind=4) :: f(k1,k2,k3) ! FastOpt
  return
end subroutine guiput

