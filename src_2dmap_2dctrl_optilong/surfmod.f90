      module surfmod
      use pumamod

      character(len=80) :: version = '26.02.2002 by larry'

      logical            :: surfdebug = .false.    ! debug printout
      integer, parameter :: nsurdim  = 100         ! max # of variables
      integer, parameter :: nsurunit =  35         ! read unit
      integer            :: nsurnum  =   0         ! total # of records
      integer            :: nlastsu  =   0         ! Last read record
      character (len=16) :: ysurnam(nsurdim) = ' ' ! Name
      integer            :: nsurptr(nsurdim) =  -1 ! Record #
      integer            :: nsurcod(nsurdim) =  -1 ! Code
      integer            :: nsurmin(nsurdim) =  13 ! First month
      integer            :: nsurmax(nsurdim) =   0 ! Last month
      integer            :: nsurlev(nsurdim) =   0 ! Max level

!     namelist parameter

      integer :: nsurf = 1
      integer :: noromax = NTRU

      real :: oroscale = 1.   ! scaling factor for orography
      real :: doro(NHOR)      ! orography

      end module surfmod


!     ===================
!     SUBROUTINE SURFCODE
!     ===================

      subroutine surfcode(kcode,yn)
      use surfmod
      integer :: kcode
      character (len=*) :: yn

!     Build a table of codes and corresponding names
!     After that variables may be read from <surface.txt>
!     referencing their name

      nsurnum = nsurnum + 1
      if (nsurnum >= nsurdim) then
          call mpabort('NSURDIM in surfmod too low')
      endif
      nsurcod(nsurnum) = kcode ! Code
      ysurnam(nsurnum) = yn    ! Name
      return
      end subroutine surfcode


!     ==================
!     SUBROUTINE SURTELL
!     ==================

      subroutine surtell(kh,krec)
      use surfmod
      integer :: kh(8)
      integer :: imon

!     The header "kh" was read at record position "krec"
!     Use the header information to update surface tables

      icod = kh(1)
      ilev = kh(2)
      imon = kh(3)
      if (imon >= 100) imon = mod(imon/100,100)
      if (imon <    0) call mpabort('imon <  0 in surtell')
      if (imon >   13) call mpabort('imon > 13 in surtell')
      do j = 1 , nsurnum
         if (kh(1) == nsurcod(j)) then
            if (nsurptr(j)    < 0) nsurptr(j) = krec
            if (ilev > nsurlev(j)) nsurlev(j) = ilev
            if (imon < nsurmin(j)) nsurmin(j) = imon
            if (imon > nsurmax(j)) nsurmax(j) = imon
         endif
      enddo
      return
      end subroutine surtell


!     ==================
!     SUBROUTINE SURLIST
!     ==================

      subroutine surlist
      use surfmod

!     Print table of variables found in <surface.txt>

      write (*,'(/,"***************************************")')
      write (*,'("* Contents of <surface.txt>           *")')
      write (*,'("***************************************")')
      write (*,'("* Code Name          Level Months Rec *")')
      write (*,'("* ----------------------------------- *")')
      do j = 1 , nsurnum
         if (nsurptr(j) >= 0) then
            if (nsurmin(j) == 13 .or. nsurmax(j) == 0 .or. &
                nsurmin(j) == nsurmax(j)) then
            write (*,'("*",i5,1x,a,i3,"       ",i4," *")') &
               nsurcod(j),ysurnam(j),nsurlev(j),nsurptr(j)
            else
            write (*,'("*",i5,1x,a,i3,i4,"-",i2,i4," *")') &
               nsurcod(j),ysurnam(j),nsurlev(j),nsurmin(j) &
               ,nsurmax(j),nsurptr(j)
            endif
         endif
      enddo
      write (*,'("***************************************",/)')

!     Print table of variables that may be read from file
!     <surface.txt>, but which where not present in this run

      write (*,'("*************************")')
      write (*,'("* Not included          *")')
      write (*,'("*************************")')
      write (*,'("* Code Name             *")')
      write (*,'("* --------------------- *")')
      do j = 1 , nsurnum
         if (nsurptr(j) < 0) then
            write (*,'("*",i5,1x,a," *")') nsurcod(j),ysurnam(j)
         endif
      enddo
      write (*,'("*************************")')

      return
      end subroutine surlist


!     ==================
!     SUBROUTINE CHECKIO
!     ==================

      subroutine checkio(kstat,ym)
      integer :: kstat               ! status of last i/o
      character (len=* ) :: ym       ! calling routine
      character (len=60) :: ymessage ! complete abort message

      if (kstat /= 0) then
         write(ymessage,'(2A)') 'I/O error in ',ym
         call mpabort(ymessage)
      endif
      return
      end


!     =========================
!     SUBROUTINE GET_SURF_ARRAY
!     =========================

      subroutine get_surf_array(yn,pa,k1,k2,kread)
      use surfmod

      character (len=*) :: yn ! array name
      real :: pa(k1,k2)       ! array to receive values
      integer :: k1           ! 1. dim of array
      integer :: k2           ! 2. dim of array
      integer :: kread        ! flag to indicate success
      integer :: imon,j
      integer :: ih(8)

      if (surfdebug) write(*,*) 'get_surf_array(',yn,'(',k1,',',k2,'))'

      kread = 0 ! Initialize as "not read"

! FastOpt      do j = 1 , nsurnum
      j=1 ! FastOpt
      do while (j<nsurnum.and.(yn .ne. trim(ysurnam(j)) .or. nsurptr(j) < 0)) ! FastOpt
         j=j+1 ! FastOpt
      end do ! FastOpt
      if (yn == trim(ysurnam(j)) .and. nsurptr(j) >= 0) then
         if (surfdebug) write(*,*) 'Seek ',yn,' rec',nsurptr(j)
         call suseek(nsurptr(j)) ! Set file position
         do j2 = 1 , k2
            read (nsurunit,'(8I10)',IOSTAT=iostat) ih(:)
            call checkio(iostat,'get_surf_array 1')
            if (ih(1) /= nsurcod(j)) then
               write(*,*) 'Seek code ',nsurcod(j),'   found ',ih(1)
               call mpabort('seek code error in get_surf_array')
            endif
            read (nsurunit,'(8E12.6)',IOSTAT=iostat) pa(:,j2)
            call checkio(iostat,'get_surf_array 2')
            nlastsu = nlastsu + 1
         enddo
         if (surfdebug) write(*,*) 'Read ',yn,' rec',nsurptr(j)
         kread = 1
! FastOpt    return
! FastOpt         endif
! FastOpt      enddo
      else ! FastOpt
         write(*,*) 'Array {',yn,'} was not in <surface.txt>'
      end if ! FastOpt
      return
      end subroutine get_surf_array


!     ======================
!     SUBROUTINE GET_3D_YEAR
!     ======================

      subroutine get_3d_year(yn,pa,kdim,klev,kread)
      use surfmod

      character (len=*) :: yn
      real :: pa(kdim,klev,0:13)
      integer :: ih(8)

      if (surfdebug) write(*,*) 'get_3d_year(',yn,'(',kdim,',',klev,'))'

      kread = 0 ! Initialize as "not read"

! FastOpt      do j = 1 , nsurnum
      j=1 ! FastOpt
      do while (j<nsurnum.and.(yn .ne. trim(ysurnam(j)) .or. nsurptr(j) < 0)) ! FastOpt
         j=j+1 ! FastOpt
      end do ! FastOpt
      if (yn == trim(ysurnam(j)) .and. nsurptr(j) >= 0) then
         if (surfdebug) write(*,*) 'Seek ',yn,' rec',nsurptr(j)
         call suseek(nsurptr(j))
         ico = nsurcod(j)
         ima = nsurmin(j)
         imz = nsurmax(j)
         ilv = nsurlev(j)
         do jmon = ima , imz
            do jlev = 1 , ilv
               read (nsurunit,'(8I10)',IOSTAT=io) ih(:)
               call checkio(io,'get_3d_year 1')
               if (ico /= ih(1)) then
                  write (*,*) 'Expecting code',ico,' found',ih(1)
                  call mpabort('get_3d_year 2')
               endif
               if (jlev /= ih(2)) then
                  write (*,*) 'Expecting level',jlev,' found',ih(2)
                  call mpabort('get_3d_year 3')
               endif
               imon = ih(3)
               if (imon >= 100) imon = mod(imon/100,100)
               if (jmon /= imon) then
                  write (*,*) 'Expecting month',jmon,' found',imon
                  call mpabort('get_3d_year 4')
               endif
               read (nsurunit,'(8E12.6)',IOSTAT=io) pa(:,jlev,jmon)
               call checkio(io,'get_3d_year 5')
               nlastsu = nlastsu + 1
               if (surfdebug) write(*,*) 'Month',imon,' for {',yn,'}'
            enddo
         enddo
         if (ima == 1 .and. imz == 12) then
            pa(:,:, 0) = pa(:,:,12)
            pa(:,:,13) = pa(:,:, 1)
            write (*,'("Array {",a,"} expanded to 14 months")') yn
         endif
         kread = 1
! FastOpt   return
! FastOpt         endif ! (yn == trim(ysurnam(j)))
! FastOpt      enddo
      else  ! FastOpt
         write(*,*) 'Array {',yn,'} was not in <surface.txt>'
         kread = 0
      end if ! FastOpt
      return
      end subroutine get_3d_year


!     ========================
!     SUBROUTINE GET_SURF_YEAR
!     ========================

      subroutine get_surf_year(yn,pa,kdim,kread)
      use surfmod

      character (len=*) :: yn
      integer :: ih(8)
      real :: pa(kdim,0:13)

      if (surfdebug) write(*,*) 'get_surf_year(',yn,'(',kdim,'))'

      kread = 0 ! Initialize as "not read"

! FastOpt      do j = 1 , nsurnum
      j=1 ! FastOpt
      do while (j<nsurnum.and.(yn .ne. trim(ysurnam(j)) .or. nsurptr(j) < 0)) ! FastOpt
         j=j+1 ! FastOpt
      end do ! FastOpt
      if (yn == trim(ysurnam(j)) .and. nsurptr(j) >= 0) then
         if (surfdebug) write(*,*) 'Seek ',yn,' rec',nsurptr(j)
         call suseek(nsurptr(j))
         ico = nsurcod(j)
         ima = nsurmin(j)
         imz = nsurmax(j)
         do jmon = ima , imz
            read (nsurunit,'(8I10)',IOSTAT=io) ih(:)
            call checkio(io,'get_surf_year 1')
            if (ico /= ih(1)) then
               write (*,*) 'Expecting code',ico,' found',ih(1)
               call mpabort('get_surf_year 2')
            endif
            imon = ih(3)
            if (imon >= 100) imon = mod(imon/100,100)
            if (jmon /= imon) then
               write (*,*) 'Expecting month',jmon,' found',imon
               call mpabort('get_surf_year 3')
            endif
            read (nsurunit,'(8E12.6)',IOSTAT=io) pa(:,jmon)
            call checkio(io,'get_surf_year 4')
            nlastsu = nlastsu + 1
            if (surfdebug) write(*,*) 'Month',imon,' for {',yn,'}'
         enddo
         if (ima == 1 .and. imz == 12) then
            pa(:, 0) = pa(:,12)
            pa(:,13) = pa(:, 1)
            write (*,'(" Array {",a,"} expanded to 14 months")') yn
         endif
         kread = 1
! FastOpt            return
! FastOpt         endif ! (yn == trim(ysurnam(j)))
! FastOpt      enddo
         else ! FastOpt
            write(*,*) 'Array {',yn,'} was not in <surface.txt>'
            kread = 0
         end if ! FastOpt
      return
      end subroutine get_surf_year


!     ======================
!     SUBROUTINE SURFACE_INI
!     ======================

      subroutine surface_ini
      use surfmod

      integer :: ih(8)
      integer :: iptr
      integer :: iskipl
      character (len=16) :: yn ! variable name

!     Define code - arrayname relationship
 
      call surfcode( 129,'doro'    )
      call surfcode( 172,'dls'     )

!     miscmod arrays

      call surfcode( 130,'dtnudge' )
      call surfcode(  22,'dfnudge' )

!     oceanmod arrays

      call surfcode( 169,'yclsst'  )
      call surfcode( 172,'yls'     )  ! same as dls
      call surfcode( 903,'yfsst'   )

!     icemod arrays

      call surfcode( 169,'xclsst'  )  ! same as zclsst
      call surfcode( 172,'xls'     )  ! same as dls
      call surfcode( 210,'xclicec' )
      call surfcode( 211,'xcliced' )
      call surfcode( 709,'xflxice' )

!     simba arrays

      call surfcode( 304,'dcveg'   )
      call surfcode( 305,'dcsoil'  )
      call surfcode(1604,'pgrow'   )
      call surfcode(1605,'plai'    )
      call surfcode(1606,'pgs'     )
      call surfcode(1607,'pz0_max' )

!     landmod arrays

      call surfcode( 173,'dz0clim' )
      call surfcode(1730,'dz0climo')
      call surfcode( 174,'dalbclim')
      call surfcode( 232,'dglac'   )
      call surfcode( 212,'dforest' )
      call surfcode( 229,'dwmax'   )
      call surfcode( 209,'dtclsoil')
      call surfcode( 169,'dtcl'    ) ! same as xclsst
      call surfcode( 140,'dwcl'    )

!     radmod arrays

      call surfcode( 200,'dqo3cl'  )


!     Scan start_data for codes and store sequence number

      iptr = 0
      open(nsurunit,file='surface.txt',form='formatted')
      io=0 ! FastOpt
! FastOpt      do
      do while (io==0) ! FastOpt
         read (nsurunit,'(8I10)',IOSTAT=io) ih(:)
! FastOpt         if (io /= 0) exit ! end-of-file
         if (io==0) then  ! FastOpt
            call surtell(ih,iptr)
            iptr = iptr + 1
            iskipl = (ih(5) * ih(6)) / 8
            do j = 1 , iskipl
               read (nsurunit,'(1X)',IOSTAT=io)
               call checkio(io,'surface_ini')
            enddo
         end if ! FastOpt
      enddo
      rewind nsurunit
      nlastsu = 0
      call surlist

!     file must be left open for further access

      return
      end subroutine surface_ini


!     ==================
!     SUBROUTINE SURFINI
!     ==================

      subroutine surfini
      use surfmod
!
!     initialize surface parameter
!
      character (len=16) :: yn ! variable name

      namelist/surfpar/nsurf,noromax,oroscale,surfdebug

      if (mypid == NROOT) then
         read(11,surfpar)
         write (*,'(/," *****************************************")')
         write (*,'(" * SURFMOD ",a29," *")') trim(version)
         write (*,'(" *****************************************")')
         write (*,'(" * Namelist SURFPAR from <puma_namelist> *")')
         write (*,'(" *****************************************")')
         write(*,surfpar)
      endif

      call mpbci(nsurf)

      if (nrestart == 0) then ! need to read start data
         call mpsurfgp('doro',doro,NHOR,1)
         call mpsurfgp('dls' ,dls ,NHOR,1)
         doro(:) = doro(:) * oroscale  ! Scale orography

!        Compute spectral orography

         call gp2fc(doro,NLON,NLPP)
         call fc2sp(doro,so)
         call mpsum(so,1)
         if (mypid == NROOT) then
            so(:) = so(:) / (cv*cv)
            if (noromax < NTRU) then
             jr=-1
             do jm=0,NTRU
              do jn=jm,NTRU
               jr=jr+2
               ji=jr+1
               if(jn > noromax) then
                so(jr)=0.
                so(ji)=0.
               endif
              enddo
             enddo
            endif ! (noromax < NTRU)
         endif ! (mypid == NROOT)

!*       adjust land sea mask  (workaround)

         dls(:)=AMAX1(dls(:),0.)
         dls(:)=AMIN1(dls(:),1.)

      endif ! (nrestart == 0)

      call landini                      ! land module
      call vegini                       ! vegetation module
      if (n_sea_points > 0) call seaini ! sea module

      return
      end subroutine surfini

!     ===================
!     SUBROUTINE SURFSTEP
!     ===================

      subroutine surfstep
      use surfmod
!$taf incompl dq,dt,dqs,dts
      call landstep
!$taf incompl dq,dt,dqs,dts
      if (n_sea_points > 0) call seastep

      return
      end subroutine surfstep

!     ===================
!     SUBROUTINE SURFSTOP
!     ===================

      subroutine surfstop
      use surfmod

      call landstop
      call vegstop
      if (n_sea_points > 0) call seastop

      return
      end subroutine surfstop


!     =================
!     SUBROUTINE SUSEEK
!     =================

      subroutine suseek(k)
      use surfmod
      integer :: ih(8)

      if (surfdebug) write (*,*) 'Pos:',nlastsu,'   Want:',k
      if (k == nlastsu) return
      if (k <  nlastsu) then
         if (surfdebug) write (*,*) 'Rewinding surface.txt'
         rewind nsurunit
         nlastsu = 0
      endif

      do jrec = nlastsu+1 , k
         read (nsurunit,'(8I10)',IOSTAT=io) ih(:)
         call checkio(io,'suseek 1')
         iskipl = (ih(5) * ih(6)) / 8
         do j = 1 , iskipl
            read (nsurunit,'(1X)',IOSTAT=io)
            call checkio(io,'suseek 2')
         enddo
      enddo
      nlastsu = k
      return
      end

