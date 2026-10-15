      module restartmod
      integer, parameter :: nresdim  = 200   ! Max number of records
      integer, parameter :: nreaunit =  33   ! FORTRAN unit for reading
      integer, parameter :: nwriunit =  34   ! FORTRAN unit for writing
      integer            :: nexcheck =   1   ! Extended checks
      integer            :: nresnum  =   0   ! Actual number of records
      integer            :: nlastrec =   0   ! Last read record
      character (len=16) :: yresnam(nresdim) ! Array of record names
      end module restartmod

!     ======================
!     SUBROUTINE RESTART_INI
!     ======================

      subroutine restart_ini(lrestart,yrfile)
      use restartmod

      logical :: lrestart
      character (len=*)  :: yrfile
      character (len=16) :: yn ! variable name

      inquire(file=yrfile,exist=lrestart)
      if (lrestart) then
         open(nreaunit,file=yrfile,form='unformatted')
         iostat=0 ! FastOpt
! FastOpt         do
         do while (iostat==0) ! changed, FastOpt
            read (nreaunit,IOSTAT=iostat) yn
! FastOpt            if (iostat /= 0) exit
            if (iostat==0) then! FastOpt
               nresnum = nresnum + 1
               yresnam(nresnum) = yn
               read (nreaunit,IOSTAT=iostat)
! FastOpt            if (iostat /= 0) exit
               if (nresnum >= nresdim) then
                  print *,'Too many variables in restart file'
                  print *,'Increase dimension NRESDIM in module restartmod'
                  print *,'*** Error Stop ***'
                  stop
               endif
            endif ! FastOpt
         enddo
   
         if (yrfile=='puma_restart') then ! added, FastOpt
         write (*,'(a,i4,3a/)') 'Found ',nresnum, &
               ' variables in file <',yrfile,'>'
         do j = 1 , nresnum
            write(*,'(i4," : ",a)') j,yresnam(j)
         enddo
         end if  ! added, FastOpt
         nlastrec = nresnum
      endif ! (lrestart)

!     file must be left open for further access

      return
      end subroutine restart_ini     


!     ==========================
!     SUBROUTINE RESTART_PREPARE
!     ==========================

      subroutine restart_prepare(ywfile)
      use restartmod

      character (len=*) :: ywfile

      open(nwriunit,file=ywfile,form='unformatted')

      return
      end subroutine restart_prepare


!     =======================
!     SUBROUTINE RESTART_STOP
!     =======================

      subroutine restart_stop
      use restartmod

      close (nreaunit)
      close (nwriunit)

      return
      end subroutine restart_stop


!     ==============================
!     SUBROUTINE GET_RESTART_INTEGER
!     ==============================

      subroutine get_restart_integer(yn,kv)
      use restartmod

      character (len=*) :: yn
      integer :: kv

 ! FastOpt      do j = 1 , nresnum
      j=1 ! FastOpt
      do while (j<nresnum.and.(trim(yn) .ne. trim(yresnam(j)))) ! FastOpt
         j=j+1 ! FastOpt
      end do ! FastOpt
      if (trim(yn) == trim(yresnam(j))) then
         call reseek(yn,j)
         read (nreaunit) kv
         nlastrec = nlastrec + 1
 ! FastOpt            return
 ! FastOpt         endif
 ! FastOpt      enddo
      else ! FastOpt
         if (nexcheck == 1) then
            print *,'*** Error in get_restart_integer ***'
            print *,'Requested integer {',yn,'} was not found'
            stop
         endif
      end if ! FastOpt
      return
      end subroutine get_restart_integer


!     ============================
!     SUBROUTINE GET_RESTART_ARRAY
!     ============================

      subroutine get_restart_array(yn,pa,k1,k2,k3)
      use restartmod

      character (len=*) :: yn
      real :: pa(k2,k3)

 ! FastOpt      do j = 1 , nresnum
      j=1 ! FastOpt
      do while (j<nresnum.and.(trim(yn) .ne. trim(yresnam(j)))) ! FastOpt
         j=j+1 ! FastOpt
      end do ! FastOpt
      if(j.le.nresnum) then ! FastOpt
         if (trim(yn) == trim(yresnam(j))) then
            call reseek(yn,j)
            read (nreaunit) pa(1:k1,:)
            nlastrec = nlastrec + 1
 ! FastOpt            return
 ! FastOpt         endif
 ! FastOpt      enddo
         else ! FastOpt
            stop "logical error in get_restart_array" ! FastOpt
         end if ! FastOpt
      else ! FastOpt
         if (nexcheck == 1) then
            print *,'*** Error in get_restart_array ***'
            print *,'Requested array {',yn,'} was not found'
            stop
         endif
      end if ! FastOpt
      return
      end subroutine get_restart_array


!     ==============================
!     SUBROUTINE PUT_RESTART_INTEGER
!     ==============================

      subroutine put_restart_integer(yn,kv)
      use restartmod

      character (len=*)  :: yn
      character (len=16) :: yy
      integer :: kv

      yy = yn
      write (nwriunit) yy
      write (nwriunit) kv
      return
      end subroutine put_restart_integer


!     ============================
!     SUBROUTINE PUT_RESTART_ARRAY
!     ============================

      subroutine put_restart_array(yn,pa,k1,k2,k3)
      use restartmod

      character (len=*)  :: yn
      character (len=16) :: yy
      integer :: k1,k2,k3
      real :: pa(k2,k3)

      yy = yn
      write (nwriunit) yy
      write (nwriunit) pa(1:k1,1:k3)
      return
      end subroutine put_restart_array


!     =================
!     SUBROUTINE RESEEK
!     =================

      subroutine reseek(yn,k)
      use restartmod

      character (len=*)  :: yn
      character (len=16) :: yy
      logical :: l_seek
      
! FastOpt      write (*,*) 'Pos:',nlastrec,'   Want:',k
      if (k <= nlastrec) then
! FastOpt         write (*,*) 'Rewinding'
         rewind nreaunit
         nlastrec = 0
      endif

      iostat=0 ! FastOpt
      l_seek=.true. ! FastOpt
! FastOpt      do
      do while (iostat==0.and.l_seek) ! FastOpt
         read (nreaunit,iostat=iostat) yy
! FastOpt         if (iostat /= 0) exit
         if (iostat==0) then ! FastOpt
            if (trim(yn) == trim(yy)) then ! FastOpt
               l_seek=.false. ! FastOpt
            else ! FastOpt
! FastOpt            if (trim(yn) == trim(yy)) return ! success
               read (nreaunit,iostat=iostat)    ! skip data
! FastOpt         if (iostat /= 0) exit
               if (iostat==0) then ! FastOpt
                  nlastrec = nlastrec + 1
               end if ! FastOpt
            end if ! FastOpt
         end if ! FastOpt
      enddo
! FastOpt      write (*,*) 'Variable <',trim(yn),'> not in restart file'
      if (l_seek.or.iostat.ne.0) write (*,*) 'Variable <',trim(yn),'> not in restart file' ! FastOpt
      return
      end


!     =========================
!     SUBROUTINE CHECK_EQUALITY
!     =========================

      subroutine check_equality(yn,pa,pb,k1,k2)
      character (len=*) :: yn
      real :: pa(k1,k2)
      real :: pb(k1,k2)

      do j2 = 1 , k2
      do j1 = 1 , k1
         if (pa(j1,j2) /= pb(j1,j2)) then
            print *,'No Equality on ',yn,'(',j1,',',j2,')',pa(j1,j2),pb(j1,j2)
            return
         endif
      enddo
      enddo
      print *,'Array {',yn,'} is OK'
      return
      end

