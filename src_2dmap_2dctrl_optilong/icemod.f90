      module icemod

      use resmod
      use pumamod, only : nstep2 ! added, FastOpt
      use pumamod, only : nstep0 ! added, FastOpt
      use pumamod, only : nstepp => nstep ! added, FastOpt
      use pumamod, only : dstore1, dstore2, dstore3
!
!     version identifier (date)
!
      character(len=80) :: version = '13.10.2005 by Larry'
!
!
!     Parameter
!
      parameter(NPRO = NPRO_ATM)        ! Number of processes (resmod)
      parameter(NLAT = NLAT_ATM)        ! Number of latitudes (resmod)
      parameter(NLON = NLAT + NLAT)     ! Number of longitudes
      parameter(NLPP = NLAT / NPRO)     ! Latitudes per process
      parameter(NHOR = NLON * NLPP)     ! Horizontal part
      parameter(NROOT = 0)              ! Master node
!
      parameter(TFREEZE=271.25)         ! freezing temp. for sea ice at S=34.7
      parameter(TMELT=273.16)           ! melting temp. for snow (0 deg C)
                                        ! ALL DENSITIES IN (kg/m**3)
      parameter(CRHOS = 1030.)          ! DENSITY OF SEA WATER AT S=34.7
      parameter(CRHOI = 920.)           ! DENSITY OF ICE
      parameter(CRHOF = 1003.8)         ! DENSITY OF 'FRESH' WATER AT S=5
      parameter(CRHOSN = 330.)          ! DENSITY OF SNOW
      parameter(CPS = 4180.)            ! SPECIFIC HEAT OF SEA WATER (J/kg*K)
      parameter(CPI = 2070.)            ! SPECIFIC HEAT OF ICE (J/(kg*K))
      parameter(CPSN = 2090.)           ! SPECIFIC HEAT OF SNOW (J/(kg*K))
      parameter(CKAPI = 2.03)           ! HEAT CONDUCTIVITY IN ICE (W/(m*K))
      parameter(CKAPSN = 0.31)          ! HEAT CONDUCTIVITY IN SNOW (W/(m*K))
      parameter(CLFI  = 3.28E5)         ! HEAT OF FUSION OF ICE (J/kg)
      parameter(CLFSN = 3.32E5)         ! HEAT OF FUSION OF SNOW (J/kg)
!
!     namelist parameters
!
      integer :: nice       = 1    ! compute ice yes/no (1/0)
      integer :: newsurf    = 0    ! 1: read surface data after restart
      integer :: nsnow      = 1    ! allow snow on ice yes/no (1/0)
      integer :: ntskin     = 1    ! compute skin temperature (0=clim.)
      integer :: ntspd      = 32   ! number of time steps per day
      integer :: noutput    = 1    ! master switch for output
      integer :: nout       = 32   ! output each nout timesteps
      integer :: nfluko     = 0    ! switch for flux correction
                                   ! (0 = none, 1 = heat-budget, 2 = newtonian)
      integer :: ncpl_ice_ocean = 1! coupling intervall ice - ocean

      integer :: nperpetual_ice = 0! perpetual climate conditions (day)
      integer :: nprint = 0        ! debug print out
      integer :: nprhor = 0        ! gp for debug printout
!
      real :: taunc         =  0.  ! time scale for newtonian cooling
      real :: xmind         = 0.1  ! minimal ice thickness (m)
      real :: stotf         = 0.054! factor to convert sal.(psu) in tfreeze (C)
!
      character(len=80) :: icefile = 'surface.txt'         ! start file
      character(len=80) :: flxfile = 'ice_flux_correction' ! flx-correction
!
!     global integer
!
      integer :: nstep      = 0    ! time step
      integer :: nrestart   = 0    ! restart switch
!
!     global real
!
      real :: xdt           = 0.   ! timestep (sec.)
      real :: cicemin       = 0.5  ! minimum compactness to be ice
      real :: solar_day  = 86400.0 ! length of day [sec]
!
!     global arrays
!
      real :: xshfl(NHOR)     = 0.   ! surface sensible heat flx
      real :: xshdt(NHOR)     = 0.   ! derivative of shfl w.r.t. temp
      real :: xlhfl(NHOR)     = 0.   ! surface latent heat flx
      real :: xlhdt(NHOR)     = 0.   ! derivative of slfl w.r.t. temp
      real :: xswfl(NHOR)     = 0.   ! surface short wave radiation
      real :: xlwfl(NHOR)     = 0.   ! surface long wave radiation
!
      real :: xts(NHOR)       = 0.   ! surface temperature (K)
      real :: xsst(NHOR)      = 0.   ! sea surface temperature (K)
      real :: xsss(NHOR)      = 35.  ! sea surface salinity (psu/ppt)
      real :: xmld(NHOR)      = 0.   ! mixed layer depth (m) (from ocean)
      real :: xls(NHOR)       = 0.   ! land sea mask (1/0)
      real :: xiced(NHOR)     = 0.   ! ice thickness (m)
      real :: xicec(NHOR)     = 0.   ! ice cover (1/0)
      real :: xsnow(NHOR)     = 0.   ! snow depth (h2o equiv. m)
      real :: xsmelt(NHOR)    = 0.   ! snow melt  (m/s water equiv.)
      real :: xsmflx(NHOR)    = 0.   ! flux for snow melt  (w/m2)
      real :: ximelt(NHOR)    = 0.   ! flux for ice melt/freeze  (w/m2)
      real :: xsndch(NHOR)    = 0.   ! snow depth change (m/s water equiv.)
      real :: xqmelt(NHOR)    = 0.   ! res. heat flux for ice melt (W/m^2)
      real :: xtfreeze(NHOR)  = 0.   ! freezing T dependent on sss (K)
      real :: xps(NHOR)       = 0.   ! atm. surface pressure (Pa)
!
      real :: xheat(NHOR)     = 0.  ! heat flux from atmosphere (W/m^2)
      real :: xcflux(NHOR)    = 0.  ! conductive heatflux through ice (W/m^2)
      real :: xprs(NHOR)      = 0.  ! snow precipitation from atmosphere (m/s)
      real :: xfresh(NHOR)    = 0.  ! fresh water flux (m/s)
      real :: xroff(NHOR)     = 0.  ! runoff (m/s)
      real :: xtaux(NHOR)     = 0.  ! zonal wind stress (pa)
      real :: xtauy(NHOR)     = 0.  ! meridional wind stress (pa)
      real :: xust3(NHOR)     = 0.  ! ustar**3 (m**3/s**3)
      real :: xoflux(NHOR)    = 0.  ! (residual) heat flux from ocean (w/m2)
      real :: xtsflux(NHOR)   = 0.  ! flux to warm/cool ice/snow (w/2)
      real :: xhnew(NHOR)     = 0.  ! cond. heatflux after skintemp (w/m2)
!
!     Climatological fields
!
      real :: xclsst(NHOR,0:13) =-1.! climatological sst
      real :: xclicec(NHOR,0:13)=-1.! climatological ice cover
      real :: xcliced(NHOR,0:13)=-1.! climatological ice thickness
      real :: xflxice(NHOR,0:13)= 0.! flux correction (W/m^2)
      real :: xclsst2(NHOR)   = 0.  ! climatological sst
      real :: xclicec2(NHOR)  = 0.  ! climatological ice cover
      real :: xcliced2(NHOR)  = 0.  ! climatological ice thickness
      real :: xflxice2(NHOR)  = 0.  ! flux correction (W/m^2)
!
!     fluxes for coupling to the ocean (accumulated)
!
      integer :: naccuo       = 0   ! counter for accumulation
      real :: cheat(NHOR)     = 0.  ! heat flux to the ocean (w/m2)
      real :: cfresh(NHOR)    = 0.  ! fresh water flux (m/s)
      real :: ctaux(NHOR)     = 0.  ! zonal wind stress (pa)
      real :: ctauy(NHOR)     = 0.  ! meridional wind stress (pa)
      real :: cust3(NHOR)     = 0.  ! ustar**3 (m**3/s**)
!
!     accumulated diagnostics
!
      integer :: naccuout     = 0   ! counter for accumulation
      real :: xflxicea(NHOR)  = 0.  ! flux correction (w/m2)
      real :: xheata(NHOR)    = 0.  ! flux from the atmosphere (w/m2)
      real :: xofluxa(NHOR)   = 0.  ! flux from the ocean (w/m2)
      real :: xqmelta(NHOR)   = 0.  ! res flux to melt ice (w/m2)
      real :: xcfluxa(NHOR)   = 0.  ! flux to the ocean (w/m2)
      real :: xsmelta(NHOR)   = 0.  ! flux for snow melt (w/m2)
      real :: ximelta(NHOR)   = 0.  ! flux used for ice melt (w/m2)
      real :: xtsfluxa(NHOR)  = 0.  ! flux to warm/cool ice/snow (w/2)
      real :: xhnewa(NHOR)    = 0.  ! cond. heatflux after skintemp (w/m2)
! ???ION???
      real :: xswrfla(NHOR)    = 0.  ! swr flux from the atmosphere (w/m2)
!
!     Parallel Stuff
!
      integer :: mpinfo  = 0
      integer :: mypid   = 0
      integer :: myworld = 0
      integer :: nproc   = NPRO
!
!     coupling with mit-ocean
!
      integer :: nmit  = 0          ! mit coupling switch (namelist)
      real    :: xsaltflx(NHOR) = 0.! salt flux into ocean
      real    :: osst(NHOR) = 0.    ! mit-sst
      real    :: osss(NHOR) = 0.    ! mit-sss

      end module icemod


!     ===========================
!     SUBROUTINE READ_ICE_SURFACE
!     ===========================

      subroutine read_ice_surface
      use icemod

      call mpsurfgp('xls',xls,NHOR,1)

      call mpsurfyear('xclsst' ,xclsst ,NHOR)
      call mpsurfyear('xclicec',xclicec,NHOR)
      call mpsurfyear('xcliced',xcliced,NHOR)

!        make sure, that land sea mask values are 0 or 1

         where (xls(:) > 0.5)
            xls(:) = 1.0
         elsewhere
            xls(:) = 0.0
         endwhere

      call mpmaxval(xclicec,NHOR,14,zmax)
      if (zmax > 1.1) then
         xclicec(:,:) = xclicec(:,:) / 100.0
         if (mypid == NROOT) &
         write (*,*) 'ice cover {xclicec} converted from % to fraction'
      endif

      if (xclicec(1,1) < 0.0) then ! xclicec was not read
         if (mypid == NROOT) &
         where (xclsst(:,:) <= TFREEZE) xclicec(:,:) = 1.0
         write (*,*) 'ice cover {xclicec} constructed from SST'
      endif

      if (xcliced(1,1) < 0.0) then ! xcliced was not read
         call make_ice_thickness
         if (mypid == NROOT) &
         write (*,*) 'ice thickness {xcliced} computed from ice cover'
      endif

!     correct climatological ice with land-sea mask

      do jm = 0 , 13
         where(xls(:) >= 1.0) ! .or. xclsst(:,jm) > TFREEZE)
            xclicec(:,jm) = 0.0
            xcliced(:,jm) = 0.0
         endwhere
      enddo
      return
      end subroutine read_ice_surface


!     =================
!     SUBROUTINE ICEINI
!     =================

      subroutine iceini(kstep,krestart,koutput,pts,psst,pmld,picec,piced,psnow)
      use icemod
!
      real :: pts(NHOR)
      real :: psst(NHOR)
      real :: pmld(NHOR)
      real :: picec(NHOR)
      real :: piced(NHOR)
      real :: psnow(NHOR)
!
      real :: zls(NLON,NLAT)
!
      real,allocatable :: zprf1(:),zprf2(:)

      namelist/icepar/nout,nfluko,nperpetual_ice,ntspd,nprint,nprhor    &
     &               ,nice,nsnow,ntskin,ncpl_ice_ocean,taunc,xmind      &
     &               ,newsurf,icefile,flxfile,nmit,stotf
!
!     copy input parameter to icemod
!
      nstep    = kstep
      nrestart = krestart
      noutput  = koutput
!
!     get process id
!
      call mpi_info(nproc,mypid)
!
!     print version number and read namelist
!
      if (mypid == NROOT) then
         open(20,file='ice_namelist',form='formatted')
         read(20,icepar)
         write (*,'(/," ***************************************")')
         write (*,'(" * ICEMOD ",a28," *")') trim(version)
         write (*,'(" ***************************************")')
         write (*,'(" * Namelist ICEPAR from <ice_namelist> *")')
         write (*,'(" ***************************************")')
         write(*,icepar)
         close(20)
      endif

      call mpbci(nout)
      call mpbci(nfluko)
      call mpbci(nice)
      call mpbci(newsurf)
      call mpbci(nsnow)
      call mpbci(ntskin)
      call mpbci(ntspd)
      call mpbci(nperpetual_ice)
      call mpbci(ncpl_ice_ocean)
      call mpbci(nprint)
      call mpbci(nprhor)
      call mpbcr(taunc)
      call mpbcr(xmind)
      call mpbcr(stotf)
      call mpbci(nmit)
!
!     set time step
!
      xdt   = solar_day / real(ntspd)
      taunc = solar_day * taunc

      if (nrestart == 0) then ! read start file
       
         call read_ice_surface
!
!        initialize
!
         call iceget
         xiced(:)=xcliced2(:)

!!LARRY
         xiced(:)=0.
!!ENDLARRY

         where (xiced(:) > xmind .and. xls(:) < 1.0)
            xicec(:)=1.0
         elsewhere
            xicec(:)=0.0
            xiced(:)=0.0 ! guokun fix this bug
         endwhere

      else ! (nrestart /= 0)
!
!        restart from restart file
!
         if (mypid == NROOT) then
            call get_restart_integer('nstep',nstep)
            call get_restart_integer('naccuice',naccuout)
            call get_restart_integer('naccuo',naccuo)
         endif
  
         call mpbci(nstep)
         call mpbci(naccuout)
         call mpbci(naccuo)
  
         call mpgetgp('xls'     ,xls     ,NHOR, 1)
         call mpgetgp('xts'     ,xts     ,NHOR, 1)
         call mpgetgp('xicec'   ,xicec   ,NHOR, 1)
         call mpgetgp('xiced'   ,xiced   ,NHOR, 1)
         call mpgetgp('xsnow'   ,xsnow   ,NHOR, 1)
         call mpgetgp('cheat'   ,cheat   ,NHOR, 1)
         call mpgetgp('cfresh'  ,cfresh  ,NHOR, 1)
         call mpgetgp('ctaux'   ,ctaux   ,NHOR, 1)
         call mpgetgp('ctauy'   ,ctauy   ,NHOR, 1)
         call mpgetgp('cust3'   ,cust3   ,NHOR, 1)
         call mpgetgp('xflxicea',xflxicea,NHOR, 1)
         call mpgetgp('xheata'  ,xheata  ,NHOR, 1)
         call mpgetgp('xofluxa' ,xofluxa ,NHOR, 1)
         call mpgetgp('xqmelta' ,xqmelta ,NHOR, 1)
         call mpgetgp('xcfluxa' ,xcfluxa ,NHOR, 1)
         call mpgetgp('xsmelta' ,xsmelta ,NHOR, 1)
         call mpgetgp('ximelta' ,ximelta ,NHOR, 1)
         call mpgetgp('xtsfluxa',xtsfluxa,NHOR, 1)
         call mpgetgp('xhnewa'  ,xhnewa  ,NHOR, 1)
         call mpgetgp('osst'    ,osst    ,NHOR, 1)
         call mpgetgp('osss'    ,osss    ,NHOR, 1)

         if (newsurf == 1) then 
            call read_ice_surface
         else
            call mpgetgp('xclicec' ,xclicec ,NHOR,14)
            call mpgetgp('xcliced' ,xcliced ,NHOR,14)
            call mpgetgp('xclsst'  ,xclsst  ,NHOR,14)
         endif ! (newsurf == 1)
      endif ! (nrestart == 0)
!
!     set xtfreeze
!
      xtfreeze(:)=TFREEZE
!
!     read flux correction
!
      if (nfluko == 1) then
         call mpsurfyear('xflxice',xflxice,NHOR)
      endif
!
!     open output file
!
      if (mypid == NROOT) then
         open(71,file='ice_output',form='unformatted')
      endif ! (mypid == NROOT)
!
!     initialize ocean
!
      call oceanini(nstep,nrestart,noutput,xsst,xmld,xoflux)
!
!     mit-coupling
!
!
!     initialize mit-ocean coupling
!
      if(nmit > 0) then
       call mpgagp(zls,xls,1)
       if(mypid == NROOT) then
        call cplatmini(zls,nlon,nlat)
       endif
       if(nrestart /= 0) then
        where(xls(:) < 1.)
         xsst(:)=osst(:)
         xsss(:)=osss(:)
         xoflux(:)=0.
! ???ION??? fix freezing temperature -1.9 degC
         xtfreeze(:)=-1.*xsss(:)*stotf+TMELT
! ???ION???
        endwhere
       endif
      endif
!
!     initialize skintemperature
!
      if (nrestart == 0) xts(:)=xsst(:)
!
!     copy output from icemod
!
      pts(:)=xts(:)
      psst(:)=xsst(:)
      pmld(:)=xmld(:)
      picec(:)=xicec(:)
      piced(:)=xiced(:)
      psnow(:)=xsnow(:)
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       call mpgagp(zprf1,xts,1)
       call mpgagp(zprf2,xsst,1)
       if(mypid==NROOT) then
        print*,'In iceini:'
        print*,'sst and ts: ',zprf2(nprhor),zprf1(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
      endif
!
      return
      end subroutine iceini

!     =====================================================================
!     SUBROUTINE icestep
!     =====================================================================

      subroutine icestep(pheat,pshfl,pshdt,plhfl,plhdt,plwfl,pswfl      &
     &                  ,pfresh,proff,pprs,ptaux,ptauy,pust3,pps        &
     &                  ,pts,picec,piced,psnow,psmelt,psndch,psst,pmld)
      use icemod
      implicit none
!
      real :: pheat(NHOR),pshfl(NHOR),pshdt(NHOR),plhfl(NHOR),plhdt(NHOR)
      real :: plwfl(NHOR),pswfl(NHOR),pfresh(NHOR),proff(NHOR),pprs(NHOR)
      real :: ptaux(NHOR),ptauy(NHOR),pust3(NHOR),pts(NHOR),picec(NHOR)
      real :: piced(NHOR),psnow(NHOR),psmelt(NHOR),psndch(NHOR),psst(NHOR)
      real :: pmld(NHOR),pps(NHOR)
!
      real :: zsnowold(NHOR),zicedold(NHOR)
      integer :: jrec ! added, FastOpt
!
      real :: zheat(NLON,NLAT) = 0.   ! net total heat flux to mit-ocean
      real :: zswr(NLON,NLAT)  = 0.   ! net short wave radiation to mit
      real :: zfresh(NLON,NLAT)= 0.   ! net fresh water flux (w/o roff) to mit 
      real :: zroff(NLON,NLAT) = 0.   ! runoff to mit
      real :: ztaux(NLON,NLAT) = 0.   ! u-stress to mit
      real :: ztauy(NLON,NLAT) = 0.   ! v-stress to mit
      real :: zps(NLON,NLAT)   = 0.   ! surface pressure to mit
      real :: zicem(NLON,NLAT) = 0.   ! snow+ice mass to mit
      real :: zicemp(NHOR)     = 0.  
      real :: zswrp(NHOR)      = 0.  
      real :: zrhoilfdt
!
!     debug arrays
!
      real,allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:)
      real,allocatable :: zprf5(:),zprf6(:)
!
!     set some helpful bits
!
      zrhoilfdt=CRHOI*CLFI/xdt
!
!     reset diagnostic arrays
!
      ximelt(:)=0.
      xsmelt(:)=0.
      xsmflx(:)=0.
      xtsflux(:)=0.
      xhnew(:)=0.
!
!     copy input to icemod
!
      where(xls(:) < 1.) ! FastOpt
         xheat(:)=pheat(:)
         xshfl(:)=pshfl(:)
         xshdt(:)=pshdt(:)
         xlhfl(:)=plhfl(:)
         xlhdt(:)=plhdt(:)
         xlwfl(:)=plwfl(:)
         xswfl(:)=pswfl(:)
         xfresh(:)=pfresh(:)
         xroff(:)=proff(:)
         xprs(:)=pprs(:)
         xtaux(:)=ptaux(:)
         xtauy(:)=ptauy(:)
         xust3(:)=pust3(:)
         xps(:)=pps(:)
      end where ! FastOpt
!
!     get climatology
!
      call iceget
!$TAF STORE xclicec2,xcliced2,xsst = inner1, REC = 1 + nstep2
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       call mpgagp(zprf1,xheat,1)
       call mpgagp(zprf2,xclicec2,1)
       call mpgagp(zprf3,xcliced2,1)
       call mpgagp(zprf4,xsst,1)
       if(mypid==NROOT) then
        print*,'In icestep:'
        print*,'sst (old time step): ',zprf4(nprhor)
        print*,'heatflx from atm: ',zprf1(nprhor)
        print*,'clim. ice c and d: ',zprf2(nprhor),zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
      endif
!
      where(xls(:) < 1.)
       zsnowold(:)=xsnow(:)
       zicedold(:)=xiced(:)
      endwhere
!
!     add new snow
!
      if(nsnow==1) then
!$TAF INCOMPLETE xsnow
       where(xicec(:) > cicemin) xsnow(:)=xsnow(:)+xprs(:)*xdt
      endif
!
!     compute surface temperature
!
!$TAF STORE xsnow = inner1, REC = 1 + nstep2
      if(ntskin==1) then
!$TAF STORE xheat,xiced,xoflux,xsnow,xsst,xts = inner1, REC = 1 + nstep2
       call skintemp
      else
       where(xiced(:)+1.E3/CRHOSN*xsnow(:) >= 0.1)
        xts(:)=xclsst2(:)
       elsewhere
        xts(:)=xsst(:)
       endwhere
      endif
!$TAF STORE xqmelt,xsnow,xsst,xts = inner1, REC = 1 + nstep2
!
!     make conductive heatflux
!
      call mkcflux
!
!     make new snow (alter conductive heatflux if neccessary)
!
      if(nsnow==1) then
       call subsnow
      else
       xsnow(:)=0.
      endif
!
!     make new ice
!
!$taf store xicec,xiced = inner1, REC = 1 + nstep2
      if(nice == 1) then
!
!     ice model
!
       call mkice

      else
!
!     climatological sea ice
!
!$TAF INCOMPLETE xiced, ximelt
       where(xiced(:) > 0. .or. xcliced2(:) > 0.)
!
!     flux correction diagnostics
!
        xflxice2(:)=(xiced(:)-xcliced2(:))*zrhoilfdt-xcflux(:)
!
!     set ice
!
        ximelt(:)=(xiced(:)-xcliced2(:))*zrhoilfdt
        xiced(:)=xcliced2(:)
!
!     reset fluxes
!
        xcflux(:)=0.
!
       elsewhere
!
        xflxice2(:)=0.
!
       endwhere
!
      endif
!
!     set sea ice compactness according to thickness
!
       where(xiced(:) > xmind)
        xicec(:)=1.
       elsewhere
        xicec(:)=0.
!        xiced(:)=0. ! guokun fix this bug
       end where
!$TAF STORE xicec = inner1, REC = 1 + nstep2
!
!     correct snow with new ice mask (ice was melted below snow)
!
      where(xicec < cicemin)
       xsmelt(:)=xsmelt(:)+xsnow(:)/xdt
       xsnow(:)=0.
      end where
!
!     correct fresh water flux into ocean due to ice/snow changes
!
      where(xls(:) < 1.) xfresh(:)=xfresh(:)                            &
     &                           -(xsnow(:)-zsnowold(:))/xdt            &
     &                           -(xiced(:)-zicedold(:))/xdt*CRHOI/1.E3
!
!     correct wind stress and ustar**3 according to ice
!
      where(xicec > cicemin)
       xtaux(:)=0.
       xtauy(:)=0.
       xust3(:)=0.
      end where
!
! ???ION??? no swr under ice
       where(xls(:) < 1.)
         zswrp(:)=xswfl(:)
       endwhere
       where(xls(:) < 1. .and. xicec(:) > cicemin)
         zswrp(:)=0.
       endwhere
! ???ION???
!
!     output
!
!     accumulate
!
      xflxicea(:)=xflxicea(:)+xflxice2(:)
      xheata(:)=xheata(:)+xheat(:)
      xofluxa(:)=xofluxa(:)+xoflux(:)
      xqmelta(:)=xqmelta(:)+xqmelt(:)
      ximelta(:)=ximelta(:)+ximelt(:)
      xsmelta(:)=xsmelta(:)+xsmflx(:)
      xcfluxa(:)=xcfluxa(:)+xcflux(:)
      xtsfluxa(:)=xtsfluxa(:)+xtsflux(:)
      xhnewa(:)=xhnewa(:)+xhnew(:)
! ???ION???    no swr under ice
!      xswrfla(:)=xswrfla(:)+xswfl(:)
      xswrfla(:)=xswrfla(:)+zswrp(:)
! ???ION???
!
      naccuout=naccuout+1
!
!     write out
!
      if(nout > 0) then
      if(mod(nstep,nout) == 0) then
       xflxicea(:)=xflxicea(:)/REAL(naccuout)
       xheata(:)=xheata(:)/REAL(naccuout)
       xofluxa(:)=xofluxa(:)/REAL(naccuout)
       xqmelta(:)=xqmelta(:)/REAL(naccuout)
       xcfluxa(:)=xcfluxa(:)/REAL(naccuout)
       ximelta(:)=ximelta(:)/REAL(naccuout)
       xsmelta(:)=xsmelta(:)/REAL(naccuout)
       xtsfluxa(:)=xtsfluxa(:)/REAL(naccuout)
       xhnewa(:)=xhnewa(:)/REAL(naccuout)
! ???ION???
       xswrfla(:)=xswrfla(:)/REAL(naccuout)
! ???ION???
       if(noutput > 0) call iceout
       xflxicea(:)=0.
       xheata(:)=0.
       xofluxa(:)=0.
       xqmelta(:)=0.
       xcfluxa(:)=0.
       ximelta(:)=0.
       xsmelta(:)=0.
       xtsfluxa(:)=0.
       xhnewa(:)=0.
! ???ION???
       xswrfla(:)=0.
! ???ION???
       naccuout=0
      endif
      endif
!
!     ocean coupling
!
!     accumulate fluxes
!
      where(xls(:) < 1.)
       cheat(:)=cheat(:)+xcflux(:)
       cfresh(:)=cfresh(:)+xfresh(:)
       ctaux(:)=ctaux(:)+xtaux(:)
       ctauy(:)=ctauy(:)+xtauy(:)
       cust3(:)=cust3(:)+xust3(:)
      end where
      naccuo=naccuo+1
!$taf store nstep = innern, REC = 1 + nstep2
      if(mod(nstep,ncpl_ice_ocean) == 0) then
       where(xls(:) < 1.)
        cheat(:)=cheat(:)/real(naccuo)
        cfresh(:)=cfresh(:)/real(naccuo)
        ctaux(:)=ctaux(:)/real(naccuo)
        ctauy(:)=ctauy(:)/real(naccuo)
        cust3(:)=cust3(:)/real(naccuo)
       endwhere
!$TAF STORE cheat = inner1, REC = 1 + nstep2

!       call oceanstep(xicec,xiced,cheat,cfresh,ctaux,ctauy,cust3,xsst   &
!     &               ,xmld,xoflux)
!
       where(xls(:) < 1.)
        cheat(:)=0.
        cfresh(:)=0.
        ctaux(:)=0.
        ctauy(:)=0.
        cust3(:)=0.
       end where
       naccuo=0
!
      endif
!
!     mit-ocean coupling
!
      if(nmit > 0) then
!
!     a) put atm.-variables to coupler
!
       call mpgagp(zheat,xcflux,1)
! ???ION??? no swr under ice
!       call mpgagp(zswr,xswfl,1)
       call mpgagp(zswr,zswrp,1)
! ???ION???
       call mpgagp(zfresh,xfresh,1)
       call mpgagp(zroff,xroff,1)
!!     call mpgagp(zps,xps,1)
!!     call mpgagp(zicem,zicemp,1)
       call mpgagp(ztaux,xtaux,1)
       call mpgagp(ztauy,xtauy,1)
       if(mypid==NROOT) then
        call cputhfla(zheat)
        call cputswra(zswr)
        call cputpmea(zfresh)
        call cputroffa(zroff)
!!      call cputpsa(zps)
!!      call cputima(zicem)
        call cputtaua(ztaux,ztauy)
       endif
!
!     b) overwrite oceanic sst with mit-values
!        (the oceanic variables are read in *getmitocean*)
!
       where(xls(:) < 1.)
        xsst(:)=osst(:)
        xsss(:)=osss(:)
        xoflux(:)=0.
! ???ION??? fix freezing temperature -1.9 degC
        xtfreeze(:)=-1.*xsss(:)*stotf+TMELT
! ???ION???
       endwhere
      endif
!
!     snow diagnostics
!
      where(xls(:) < 1.) xsndch(:)=(xsnow(:)-zsnowold(:))/xdt
!
!     copy output from icemod
!
      pts(:)=xts(:)
      picec(:)=xicec(:)
      piced(:)=xiced(:)
      psnow(:)=xsnow(:)
      psmelt(:)=xsmelt(:)
      psndch(:)=xsndch(:)
      psst(:)=xsst(:)
      pmld(:)=xmld(:)
!
!     advance time step
!
      nstep=nstep+1
!
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       call mpgagp(zprf1,xoflux,1)
       call mpgagp(zprf2,xicec,1)
       call mpgagp(zprf3,xiced,1)
       call mpgagp(zprf4,xsst,1)
       call mpgagp(zprf5,xsss,1)
       call mpgagp(zprf6,xtfreeze,1)
       if(mypid==NROOT) then
        print*,'final ice c and d to atm and oce: ',zprf2(nprhor),zprf3(nprhor)
        print*,'heatflx and sst from ocean: ',zprf1(nprhor),zprf4(nprhor)
        print*,'sss and tfreeze from ocean: ',zprf5(nprhor),zprf6(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
      endif
!
      return
      end subroutine icestep

!     ===================================
!     SUBROUTINE ICESTOP
!     ===================================

      subroutine icestop
      use icemod
!
!     close output file
!
      if (mypid == NROOT) then
       close(71)
      endif
!
!     write restart file
!
      if (mypid == NROOT) then
         call put_restart_integer('naccuice',naccuout)
         call put_restart_integer('naccuo',naccuo)
      endif

      call mpputgp('xls'     ,xls     ,NHOR, 1)
      call mpputgp('xts'     ,xts     ,NHOR, 1)
      call mpputgp('xicec'   ,xicec   ,NHOR, 1)
      call mpputgp('xiced'   ,xiced   ,NHOR, 1)
      call mpputgp('xsnow'   ,xsnow   ,NHOR, 1)
      call mpputgp('xclicec' ,xclicec ,NHOR,14)
      call mpputgp('xcliced' ,xcliced ,NHOR,14)
      call mpputgp('xclsst'  ,xclsst  ,NHOR,14)
      call mpputgp('cheat'   ,cheat   ,NHOR, 1)
      call mpputgp('cfresh'  ,cfresh  ,NHOR, 1)
      call mpputgp('ctaux'   ,ctaux   ,NHOR, 1)
      call mpputgp('ctauy'   ,ctauy   ,NHOR, 1)
      call mpputgp('cust3'   ,cust3   ,NHOR, 1)
      call mpputgp('xflxicea',xflxicea,NHOR, 1)
      call mpputgp('xheata'  ,xheata  ,NHOR, 1)
      call mpputgp('xofluxa' ,xofluxa ,NHOR, 1)
      call mpputgp('xqmelta' ,xqmelta ,NHOR, 1)
      call mpputgp('xcfluxa' ,xcfluxa ,NHOR, 1)
      call mpputgp('xsmelta' ,xsmelta ,NHOR, 1)
      call mpputgp('ximelta' ,ximelta ,NHOR, 1)
      call mpputgp('xtsfluxa',xtsfluxa,NHOR, 1)
      call mpputgp('xhnewa'  ,xhnewa  ,NHOR, 1)
      call mpputgp('osst'    ,osst    ,NHOR, 1)
      call mpputgp('osss'    ,osss    ,NHOR, 1)
!
!     finalize ocean
!
      call oceanstop
!
      return
      end subroutine icestop

!     =====================================================================
!     SUBROUTINE mkice
!     =====================================================================

      subroutine mkice
      use icemod
!
      parameter(zdcool=2.) ! depth of the ocean ML intercting with the ice
!
      real :: zmelt(NHOR) = 0.
!
!     debug arrays
!
      real, allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:)
      real, allocatable :: zprf5(:),zprf6(:)
!
      zrhoilfdt=CRHOI*CLFI/xdt
!
!     make flux correction
!
      if(nfluko == 1) call getflx
      if(nfluko == 2) call mkflukoi
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       call mpgagp(zprf1,xcflux,1)
       call mpgagp(zprf2,xflxice2,1)
       if(mypid==NROOT) then
        print*,'In mkice:'
        print*,'conductive hf and fluko: ',zprf1(nprhor),zprf2(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
      endif
!
!     add flux correction to conductive heat flux
!
      where(xls(:) < 1.) xcflux(:)=xcflux(:)+xflxice2(:)
!
!     make new ice
!
      where((xicec(:) > cicemin) .or. (xsst(:) <= xtfreeze(:)))
!
       ximelt(:)=xcflux(:)
       xiced(:)=xiced(:)-ximelt(:)/zrhoilfdt
!
!      reset heat fluxes
!
       xcflux(:)=0.
!
      end where
!
!     melt sea ice if sst > tfreeze and take heat from ocean
!
      zmelt(:)=0.
      where((xicec(:) > cicemin) .and. (xsst(:) > xtfreeze(:)))
! ???ION???  BUGGY
       zmelt(:)=(xsst(:)-xtfreeze(:))*CRHOS*zdcool*CPS/xdt
!       zmelt(:)=(xsst(:)-xtfreeze(:))*CRHOS*zdcool/xdt
! ???ION???
       ximelt(:)=ximelt(:)+zmelt(:)
       xcflux(:)=xcflux(:)-zmelt(:)
       xiced(:)=xiced(:)-zmelt(:)/zrhoilfdt
      end where
!      dstore1(:)=xiced(:) 
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       call mpgagp(zprf1,xcflux,1)
       call mpgagp(zprf2,ximelt,1)
       call mpgagp(zprf3,xiced,1)
       call mpgagp(zprf4,zmelt,1)
       call mpgagp(zprf5,xtfreeze,1)
       if(mypid==NROOT) then
        print*,'tfreeze: ',zprf5(nprhor)
        print*,'flux from ocean because t > tf: ',zprf4(nprhor)
        print*,'new ice and flux for ice: ',zprf3(nprhor),zprf2(nprhor)
        print*,'new (residual) conductive hf: ',zprf1(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
      endif
!
!     correct negative ice thickness (warm ocean)
!
!$TAF STORE xiced = inner1, REC = 1 + nstep2
      where(xiced(:) <= 0.)
       xcflux(:)=xcflux(:)-xiced(:)*zrhoilfdt
       ximelt(:)=ximelt(:)+xiced(:)*zrhoilfdt
       xiced(:)=0.
      end where
      
!
!     set infinitisimal sea ice to zero
!
      where(ABS(xiced(:)) < 1.E-9) xiced(:)=0.
!
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       call mpgagp(zprf1,xcflux,1)
       call mpgagp(zprf2,ximelt,1)
       call mpgagp(zprf3,xiced,1)
       if(mypid==NROOT) then
        print*,'final ice: ',zprf3(nprhor)
        print*,'final flux for ice: ',zprf2(nprhor)
        print*,'final conductive hf: ',zprf1(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
      endif
!
      return
      end subroutine mkice

!     =====================================================================
!     SUBROUTINE subsnow
!     =====================================================================

      subroutine subsnow
      use icemod
!
      real :: zhice(NHOR) = 0.    ! new ice thickness
!
      zrhosls=CRHOSN*CLFSN
      zrhoilfdt=CRHOI*CLFI/xdt
!
!     snow melt
!
      do jhor=1,NHOR
       if(xicec(jhor) > cicemin) then
        if(xqmelt(jhor) > (1.E3/CRHOSN)*xsnow(jhor)*zrhosls/xdt) then
         xsmelt(jhor)=xsnow(jhor)/xdt
         xqmelt(jhor)=xqmelt(jhor)-(1.E3/CRHOSN)*xsnow(jhor)*zrhosls/xdt
         xsmflx(jhor)=(1.E3/CRHOSN)*xsnow(jhor)*zrhosls/xdt
         xsnow(jhor)=0.
        elseif(xqmelt(jhor) > 0.) then
         xsmelt(jhor)=CRHOSN/1.E3*xqmelt(jhor)/zrhosls
         xsnow(jhor)=xsnow(jhor)-xsmelt(jhor)*xdt
         xsmflx(jhor)=xqmelt(jhor)
         xqmelt(jhor)=0.
        endif
!
!     add residual flux to conductive heat flux and reset qmelt
!
        xcflux(jhor)=xcflux(jhor)+xqmelt(jhor)
        xqmelt(jhor)=0.
!
       else
!
!       snow falls into water
!
        xsmelt(jhor)=xsnow(jhor)/xdt
        xsnow(jhor)=0.
!
       endif
      enddo
!
!     convert snow to sea ice if snow/ice interface is below sea level:
!     (note: xsnow = h2o equiv.) (new ice is build by modified cflux)
!
      where(xls(:) < 1. .and. 1.E3*xsnow(:) > (CRHOS-CRHOI)*xiced(:))
       zhice(:)=(1.E3*xsnow(:)+CRHOI*xiced(:))/CRHOS
       xsnow(:)=(CRHOS-CRHOI)*zhice(:)/1.E3
       xcflux(:)=xcflux(:)-(zhice(:)-xiced(:))*zrhoilfdt
      endwhere
!
      return
      end subroutine subsnow

!     =====================================================================
!     SUBROUTINE MKCFLUX
!     =====================================================================

      subroutine mkcflux
      use icemod
!
      real :: zckap(NHOR)
      real :: zhsnow(NHOR)
!
!     debug arrays
!
      real, allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:)
!
!     compute conductive heat flux
!
      where(xls(:) < 1.) xcflux(:)=xheat(:)
      where(xicec(:) > cicemin)
       zhsnow(:)=1.E3/CRHOSN*xsnow(:)
       zckap(:)=(CKAPSN*zhsnow(:)+CKAPI*xiced(:))/(zhsnow(:)+xiced(:))
       xcflux(:)=zckap(:)*(xts(:)-xsst(:))/(zhsnow(:)+xiced(:))
      end where
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       call mpgagp(zprf1,zhsnow,1)
       call mpgagp(zprf2,zckap,1)
       call mpgagp(zprf3,xcflux,1)
       if(mypid==NROOT) then
        print*,'in mkcflux: '
        print*,'snow depth (in m snow): ',zprf1(nprhor)
        print*,'kappa and conductive hf: ',zprf2(nprhor),zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
      endif
!
      return
      end subroutine mkcflux

!     =====================================================================
!     SUBROUTINE iceout
!     =====================================================================

      subroutine iceout
      use icemod
!
      integer :: ih(8)
!
      real :: zsnow(NHOR) = 0.
!
      where (xls(:) < 0.5) zsnow(:) = 1.e3/CRHOSN *xsnow(:)
!
      call ntomin(nstep,nmin,nhour,nday,nmonth,nyear)
!
      ih(2) = 0
      ih(3) = nyear*10000 +nmonth*100 +nday
      ih(4) = nhour*100 +nmin
      ih(5) = NLON
      ih(6) = NLAT
      ih(7) = 0
      ih(8) = 0
!
      ih(1) = 701
      call mpwritegph(71,xheata,NHOR,1,ih)
      ih(1) = 702
      call mpwritegph(71,xhnewa,NHOR,1,ih)
      ih(1) = 703
      call mpwritegph(71,xtsfluxa,NHOR,1,ih)
      ih(1) = 704
      call mpwritegph(71,xsmelta,NHOR,1,ih)
      ih(1) = 705
      call mpwritegph(71,ximelta,NHOR,1,ih)
      ih(1) = 706
      call mpwritegph(71,xofluxa,NHOR,1,ih)
      ih(1) = 707
      call mpwritegph(71,xcfluxa,NHOR,1,ih)
      ih(1) = 708
      call mpwritegph(71,xqmelta,NHOR,1,ih)
      ih(1) = 709
      call mpwritegph(71,xflxicea,NHOR,1,ih)
      ih(1) = 710
      call mpwritegph(71,xicec,NHOR,1,ih)
      ih(1) = 711
      call mpwritegph(71,xiced,NHOR,1,ih)
      ih(1) = 739
      call mpwritegph(71,xts,NHOR,1,ih)
      ih(1) = 741
      call mpwritegph(71,zsnow,NHOR,1,ih)
      ih(1) = 769
      call mpwritegph(71,xsst,NHOR,1,ih)
      ih(1) = 772
      call mpwritegph(71,xls,NHOR,1,ih)
      ih(1) = 790
      call mpwritegph(71,xclicec2,NHOR,1,ih)
      ih(1) = 791
      call mpwritegph(71,xcliced2,NHOR,1,ih)
! ???ION??? new vars in output
      ih(1) = 798
      call mpwritegph(71,xsss,NHOR,1,ih)
      ih(1) = 799
      call mpwritegph(71,xswrfla,NHOR,1,ih) 
! ???ION???
!
      return
      end subroutine iceout

!====================================================================
!     SUBROUTINE MKFLUKOI
!====================================================================

      subroutine mkflukoi
      use icemod
!
      zrhoilf=CRHOI*CLFI
!
      xflxice2(:)=0.
      if(taunc > 0.) then
       where(xls(:) < 1.)
        xflxice2(:)=(xiced(:)-xcliced2(:))*zrhoilf/taunc
       end where
      else
       where(xls(:) < 1.)
        xflxice2(:)=(xiced(:)-xcliced2(:))*zrhoilf/xdt                  &
     &             -xcflux(:)
       end where
      endif
!
      return
      end subroutine mkflukoi

!====================================================================
!     SUBROUTINE SKINTEMP
!====================================================================

      subroutine skintemp
      use icemod
!
      parameter(stb=5.67E-8)     ! stefan boltzmann constant
      parameter(zlayer=0.1)      ! min ice+snow thickness to compute skint.
!
      real :: ztso(NHOR)         ! old ts for diagnostics
      real :: zhsnow(NHOR)       ! snow depth (m snow)
      real :: zckap_mean(NHOR)   ! kappa (ice+snow layer)
!
!     debug arrays
!
      real, allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:),zprf5(:)
      real, allocatable :: zprf6(:),zprf7(:),zprf8(:),zprf9(:),zprf10(:)
!
!     ONLY VALID FOR ICE COVER
!     SST GIVEN BY OCEAN MODEL
!
!     F.LUNKEIT   UNIHH     FEB-07
!
!     PURPOSE.
!     --------
!     CALCULATE ICE SKIN-TEMPERATURE AS A PROGNOSTIC VARIABLE
!
!     preset fluxes
!
      xhnew(:)=xheat(:)
      xtsflux(:)=0.
      xqmelt(:)=0.
!
!     save old ts
!
      ztso(:)=xts(:)
!
!     start calculations
!
      do jhor=1,NHOR
       zhsnow(jhor)=1.E3/CRHOSN*xsnow(jhor)
       if(xiced(jhor)+zhsnow(jhor) >= zlayer) then
        zlsnow=AMIN1(zhsnow(jhor),zlayer)
        zlice=zlayer-zlsnow
        zcpdt=(CRHOI*CPI*zlice+CRHOSN*CPSN*zlsnow)/xdt
        zckap_mean(jhor)=(CKAPSN*zhsnow(jhor)+CKAPI*xiced(jhor))        &
     &                  /(zhsnow(jhor)+xiced(jhor))
!
!     new skin temperature (implicit w.r.t. heat conduction)
!
        zkapz=zckap_mean(jhor)/(xiced(jhor)+zhsnow(jhor))
        zflx=xheat(jhor)+xoflux(jhor)+zkapz*xsst(jhor)
        xts(jhor)=(zcpdt*xts(jhor)+zflx)/(zcpdt+zkapz)
!
!     residual flux going into ice/snow melt
!
        if(xts(jhor) > TMELT) then
         xqmelt(jhor)=xheat(jhor)+xoflux(jhor)+zkapz*(xsst(jhor)-TMELT) &
     &               -zcpdt*(TMELT-ztso(jhor))
         xts(jhor)=TMELT
        endif
!
!     diagnose fluxes
!
        xhnew(jhor)=zkapz*(xts(jhor)-xsst(jhor))
        xtsflux(jhor)=zcpdt*(xts(jhor)-ztso(jhor))
!
       else
        xts(jhor)=xsst(jhor)
        if(xicec(jhor) > 0.) then
         xqmelt(jhor)=xheat(jhor)+xoflux(jhor)
        endif
       endif
      enddo

!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       allocate(zprf8(NLON*NLAT))
       allocate(zprf9(NLON*NLAT))
       call mpgagp(zprf1,zhsnow,1)
       call mpgagp(zprf2,zckap_mean,1)
       call mpgagp(zprf3,xheat,1)
       call mpgagp(zprf4,xoflux,1)
       call mpgagp(zprf5,ztso,1)
       call mpgagp(zprf6,xts,1)
       call mpgagp(zprf7,xhnew,1)
       call mpgagp(zprf8,xtsflux,1)
       call mpgagp(zprf9,xqmelt,1)
       if(mypid==NROOT) then
        print*,'in skintemp: '
        print*,'snow depth (in m snow): ',zprf1(nprhor)
        print*,'kappa: ',zprf2(nprhor)
        print*,'atm. heat flux: ',zprf3(nprhor)
        print*,'oce. (residual) heat flux: ',zprf4(nprhor)
        print*,'old and new Ts: ',zprf5(nprhor),zprf6(nprhor)
        print*,'conductive heat flux (w. new ts): ',zprf7(nprhor)
        print*,'heat flux used to warm/cool ice: ',zprf8(nprhor)
        print*,'residual flux to melt snow/ice: ',zprf9(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
       deallocate(zprf7)
       deallocate(zprf8)
       deallocate(zprf9)
      endif
!
      return
      end subroutine skintemp

!====================================================================
!      SUBROUTINE ICEGET
!====================================================================

      subroutine iceget
      use icemod
!
      jstep=nstep+1   !advance time step (ts = ts(t)+ dtsdt)
      call ntomin(jstep,nmin,nhour,nday,nmonth,nyear)
!
      if(nperpetual_ice > 0) then
       nmonth=(nperpetual_ice-1)/30.+1
       nday=nperpetual_ice-nmonth*30
      endif
!
      jm1=nmonth
      if (nday.gt.15) then
       jm2=jm1+1
      else
       jm2=jm1-1
      endif
!
!     interpolate
!
      zgw2=abs((nday-1)*1440+nhour*60+nmin-21600)/43200.
! ???ION???
      if (zgw2 > 1.0) zgw2 = 1.0
! ???ENDION???
      zgw1=1.-zgw2
!
      xclsst2(:)=zgw1*xclsst(:,jm1)+zgw2*xclsst(:,jm2)
      xcliced2(:)=zgw1*xcliced(:,jm1)+zgw2*xcliced(:,jm2)
      where(xcliced2(:) > xmind)
       xclicec2(:)=1.
      elsewhere
       xclicec2(:)=0.
      end where
!
      return
      end subroutine iceget

!====================================================================
!      SUBROUTINE getflx
!====================================================================

      subroutine getflx
      use icemod
!
      jstep=nstep+1   !advance time step (ts = ts(t)+ dtsdt)
      call ntomin(jstep,nmin,nhour,nday,nmonth,nyear)
!
      if(nperpetual_ice > 0) then
       nmonth=(nperpetual_ice-1)/30.+1
       nday=nperpetual_ice-nmonth*30
      endif
!
      jm1=nmonth
      if (nday.gt.15) then
       jm2=jm1+1
       if (jm1.eq.12) jm2=1
      else
       jm2=jm1-1
       if (jm2.eq.0) jm2=12
      endif
!
!     interpolate
!
      zgw2=abs((nday-1)*1440+nhour*60+nmin-21600)/43200.
! ???ION???
      if (zgw2 > 1.0) zgw2 = 1.0
! ???ENDION???
      zgw1=1.-zgw2
      xflxice2(:)=zgw1*xflxice(:,jm1) +zgw2*xflxice(:,jm2)
!
      return
      end subroutine getflx


!     =============================
!     SUBROUTINE MAKE_ICE_THICKNESS
!     =============================

      subroutine make_ice_thickness
      use icemod

      parameter(cminn=0.1 ,cmins=0.25)
      parameter(cmaxn=0.9 ,cmaxs=1.0 )
      parameter(hminn=0.25,hmins=0.25)
      parameter(hmaxn=3.0 ,hmaxs=0.50)

      real :: zc(NLON,NLAT,0:13)
      real :: zd(NLON,NLAT,0:13)

      real :: zhfac(0:13)=(/0.912,0.942,1.,1.058,1.124,1.161,1.175,1.058,0.931 &
                        ,0.883,0.88,0.876,0.912,0.942/)

!     convert ice compactness to thickness (see CCM3 report pp 127-129)

      call mpgagp(zc,xclicec,14)
      zd(:,:,:) = 0.0
      
      if (mypid == NROOT) then
         do jm = 0 , 13
   
      !     northern hemisphere
      
            do jlat = 1 , NLAT/2
               where (zc(:,jlat,jm) >= cmaxn)
                  zd(:,jlat,jm)=hmaxn*zhfac(jm)
               elsewhere (zc(:,jlat,jm) >= cminn)
                  zd(:,jlat,jm)=(hminn+(hmaxn-hminn)*(zc(:,jlat,jm)-cminn) &
                               /(cmaxn-cminn))*zhfac(jm)
               elsewhere (zc(:,jlat,jm) > 0.0)
                  zd(:,jlat,jm)=hminn*zhfac(jm)
               endwhere
            enddo
      
      !     southern hemisphere
      
            do jlat = NLAT/2+1 , NLAT
               where (zc(:,jlat,jm) >= cmins)
                  zd(:,jlat,jm)=hmaxs
               elsewhere (zc(:,jlat,jm) > 0.0)
                  zd(:,jlat,jm)=hmins+zc(:,jlat,jm)
               endwhere
            enddo
   
         enddo ! jm
      endif ! (mypid == NROOT)

      call mpscgp(zd,xcliced,14)

      return
      end subroutine make_ice_thickness


!     =========================
!     SUBROUTINE ICE_INI_PLANET
!     =========================

      subroutine ice_ini_planet(pspd)

!     The ice module defaults to initialization for planet Earth
!     For other planets (e.g. Mars or Venus) this routine must
!     be called directly after "iceini", in order to adjust
!     the constants to the other planet

      use icemod

      solar_day = pspd

      return
      end


