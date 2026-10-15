

!-*-f90-*-
program writexb
  implicit none
  real, allocatable :: x(:), x1(:)
  integer :: n, i,j,iargc,getarg,nargs,iostat
  character(len=16) :: y
  external getarg, iargc
  nargs = iargc()
  if (nargs > 0) then
     i= getarg(1,y)
     read(y,'(i)',iostat=iostat) n
     if (iostat.ne.0 .or. nargs.ne.n+1 .and. nargs.ne.1 .and. nargs.ne.2 ) then
        stop "Usage: writexb [ n [ param_1 param_2 ... param_n ] ]"
     end if
     allocate(x(n))
     select case (nargs)
     case (1)
        x(:) = 0.
     case (2)
        i= getarg(2,y)
        read(y,*) x(1)
        x(:) = x(1)
     case default
        do j=1,n
           i= getarg(1+j,y)
!           print*,n,i,j,y
           read(y,*) x(j)
        end do
     end select
  else
     allocate(x(118806))
     allocate(x1(118806))
     open(10,file="x.b_opti",form="unformatted")
        read (10) x
     close (10)
     print*,'x_opti: ',x(1:6)
!     n = 122902 ! set by makefile macro
!     allocate(x(n))
!     x(1:6) = 0.2
!     x(7:122902)=0.
  end if
  x1(:)=x;
  x1(1)=0.;
  open(10,file="x.b_001",form="unformatted")
  write (10) x1
  close (10)
  print*,x1(1:6)

  x1(:)=x;
  x1(2)=0.;
  open(10,file="x.b_002",form="unformatted")
  write (10) x1
  close (10)
  print*,x1(1:6)

  x1(:)=x;
  x1(3)=0.;
  open(10,file="x.b_003",form="unformatted")
  write (10) x1
  close (10)
  print*,x1(1:6)

  x1(:)=x;
  x1(4)=0.;
  open(10,file="x.b_004",form="unformatted")
  write (10) x1
  close (10)
  print*,x1(1:6)

  x1(:)=x;
  x1(5)=0.;
  open(10,file="x.b_005",form="unformatted")
  write (10) x1
  close (10)
  print*,x1(1:6)

  x1(:)=x;
  x1(6)=0.;
  open(10,file="x.b_006",form="unformatted")
  write (10) x1
  close (10)
  print*,x1(1:6)

  x1(:)=x;
  x1(7:118806)=0.;
  open(10,file="x.b_007",form="unformatted")
  write (10) x1
  close (10)
  print*,x1(1:6)

  print*,"ok"


end program writexb
