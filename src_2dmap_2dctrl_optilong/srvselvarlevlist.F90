!-*-f90-*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                          !!
!! PURPOSE: Select only records with code and levels in a given list.       !!
!!          Reads and writes unformatted SERVICE data.                      !!
!!                                                                          !!
!! INPUT:   std.input:  binary data in SERVICE format                       !!
!!                                                                          !!
!! OUTPUT:  std.output: binary data in SERVICE format                       !!
!!                                                                          !!
!! COMPILE WITH: f90 srvselvarlevlist.F -Vaxlib  (intel)                    !!
!!               f90 srvselvarlevlist.F -lU77    (absoft)                   !!
!!               f90 srvselvarlevlist.F -lc      (some others)              !!
!!               pgf90 srvselvarlevlist.F        (portland group)           !!
!!                                                                          !!
!! USAGE:  srvselvarlevlist input.srv  output.srv  < list.txt               !!
!!                                                                          !!
!! AUTHOR: Simon Blessing (http://www.fastopt.com)                          !!
!!                                                                          !!
!! LAST MAINTENANCE: snb 2010-03-24                                         !!
!!                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program srvselvarlevlist
        implicit none
        type t_codeset
           integer :: code
           integer :: nlev
           integer, allocatable :: levlst(:)
        end type t_codeset
        integer, allocatable :: itime (:)
        integer :: ihtime, iarg, i, nsets
        character*8 :: ydate
        character(len=255) :: inputfile, outputfile
        integer(kind=4), dimension(8) :: ih
        real(kind=4), allocatable :: f(:)
        integer :: iargc,getarg,id,iostat
        logical :: lcopy
        type(t_codeset), allocatable :: codeset(:)
        external getarg,iargc

! parse command line
        select case (iargc())
        case (2)
           id=getarg(1,inputfile)
           id=getarg(2,outputfile)
        case default
           call usage
           stop
        end select

! read standard input
        read(*,*) nsets ! number of sets
        allocate(codeset(nsets))
        do i=1,nsets
           read(*,*) codeset(i)%code  ! code of set
           read(*,*) codeset(i)%nlev  ! number of levels of set
           allocate(codeset(i)%levlst(codeset(i)%nlev))
           read(*,*) codeset(i)%levlst(:) ! list of levels of set
        end do

! open files
        print'(a,a,a)',"Opening ",trim(inputfile)," for reading."
        open(10,file=inputfile,form="unformatted",status="old")
        
        print'(a,a,a)',"Opening ",trim(outputfile)," for writing."
        open(20,file=outputfile,form="unformatted",status="unknown")

        read(10,iostat=iostat) ih
        if (.not.allocated(f)) allocate(f(ih(5)*ih(6)))
        do while (iostat==0)
           lcopy = .false.
           do id=1,nsets
              if ((ih(1).eq.codeset(id)%code).and. &
                   any(ih(2).eq.codeset(id)%levlst)) lcopy = .true.
           end do

           if (lcopy) then
              read(10) f(:)
! ** hack **
              if (ih(1) == 31) ih(1) = 210 ! code translation
! /hack
              write(20) ih(:)
              write(20) f(:)
           else
              read(10)
           end if

        read(10,iostat=iostat) ih

        end do
           
        close(10)
        close(20)

        do i=1,nsets
           deallocate(codeset(i)%levlst)
        end do
        deallocate(codeset)

      end program srvselvarlevlist

      subroutine usage
      open(77,file="/dev/stderr",form="formatted")
      write(77,*)"Selects records from a SERVICE file with "
      write(77,*)"headers matching codes and levels in a list"
      write(77,*)"read from standard input."
      write(77,*)"USAGE: srvselvarlevlist input.srv output.srv < list.txt"
      end subroutine usage
