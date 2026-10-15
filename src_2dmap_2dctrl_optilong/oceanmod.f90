      module oceanmod

      use resmod
      use pumamod, only : nstep2
!
!     version identifier (date)
!
      character(len=80) :: version = '13.10.2005 by Larry'
!
!     Parameter
!

      parameter(NPRO = NPRO_ATM)        ! Number of processes (resmod)
      parameter(NLAT = NLAT_ATM)        ! Number of latitudes (resmod)

      parameter(NLEV_OCE = 1)           ! Number of Layers
      parameter(NLON = NLAT + NLAT)     ! Number of longitudes
      parameter(NLPP = NLAT / NPRO)     ! Latitudes per process
      parameter(NHOR = NLON * NLPP)     ! Horizontal part
      parameter(NROOT = 0)              ! Master node
      parameter(CRHOS=1030.)            ! Density of sea water (kg/m**3)
      parameter(CPS=4180.)              ! Specific heat of sea water (J/kg*K)
      parameter(TFREEZE=271.25)         ! Freezing point (K)
!
!     namelist parameter
!
      integer :: nfluko           = 0   ! switch for flux correction
                                        ! (0=none,1=heat-budget,2=newtonian)
      integer :: ndiag            = 480 ! diagnostics each ndiag timesteps
      integer :: noutput          = 1   ! master switch for output: 0= no output
      integer :: nout             = 32  ! afterburner output each nout timesteps
      integer :: nocean           = 1   ! compute ocean yes/no
      integer :: newsurf          = 0   ! update surface arrays at restart
      integer :: ntspd            = 32  ! ocean timesteps per day
      integer :: nperpetual_ocean = 0   ! perpetual climate conditions
      integer :: nprint           = 0   ! print debug information
      integer :: nprhor           = 0   ! gp to print debug information
!
      real :: dlayer(NLEV_OCE)      = 50.   ! layer depth (m)
      real :: taunc             =  0.   ! newtonian cooling timescale (d)
      real :: vdiffk            = 1.E-4 ! vertikal diffusion coeff. [m**2/s]
!
!     global integer
!
      integer :: nstep    = 0           ! time step
      integer :: naccuout = 0           ! counter for accumulated output
      integer :: nrestart = 0           ! switch for restart
!
!     global reals
!
      real :: dtmix                     ! time step (s)
      real :: solar_day    = 86400.0    ! 24 * 60 * 60 (for Earth)
!
      real :: gw(NHOR)                  ! gaussian weights
!
      real :: ysst(NHOR,NLEV_OCE) = 0.      ! temperature (K)
      real :: ymld(NHOR,NLEV_OCE) = 0.      ! layer depth (m)
      real :: yls(NHOR)   = 0.          ! land sea mask (frac.)
      real :: yicec(NHOR) = 0.          ! ice cover (frac.)
      real :: yiced(NHOR) = 0.          ! ice thickness (frac.)
      real :: yheat(NHOR) = 0.          ! heat flux from ice/atm. (W/m**2)
      real :: yfresh(NHOR)= 0.          ! fresh water flux (m/s)
      real :: ytaux(NHOR) = 0.          ! zonal wind stress (pa)
      real :: ytauy(NHOR) = 0.          ! meridional wind stress (pa)
      real :: yust3(NHOR) = 0.          ! friction velocity**3 (m**3/s**3)
      real :: yiflux(NHOR)= 0.          ! heat flux into ice (residual) (W/m2)
      real :: ydsst(NHOR) = 0.          ! heat flux from vdiff (w/m2)

      real :: yclsst(NHOR,0:13)         ! climatological sst (K)
      real :: yfsst(NHOR,0:13) = 0.     ! flux corr. sst (W/m**2)

      real :: yclsst2(NHOR) = 0.        ! climatological sst (K)
      real :: yfsst2(NHOR) = 0.         ! flux corr. sst (W/m**2)
!
!     diagnostics (accum.)
!
      real :: yheata(NHOR) = 0.        ! accum. heat flux
      real :: yfssta(NHOR) = 0.        ! accum. flux correction
      real :: yifluxa(NHOR)= 0.        ! accum. heat flux into ice
      real :: ydssta(NHOR)= 0.         ! accum. heat flux from vdiff
!
!     Parallel Stuff
!
      integer :: mpinfo  = 0
      integer :: mypid   = 0
      integer :: myworld = 0
      integer :: nproc   = NPRO

      end module oceanmod


!     ===================================
!     SUBROUTINE OCEANINI
!     ===================================
!
      subroutine oceanini(kstep,krestart,koutput,psst,pmld,piflux)
      use oceanmod
!
      real :: psst(NHOR),pmld(NHOR),piflux(NHOR)
      real :: zsi(NLAT),zgw(NLAT),zgw2(NLON,NLAT)
!
      namelist/oceanpar/ndiag,nout,nfluko,ntspd,nocean,nprint,nprhor    &
    &                  ,nperpetual_ocean,taunc,dlayer   &
    &                  ,vdiffk,newsurf

!
!     get process id
!
      call mpi_info(nproc,mypid)
!
!     compute gaussian weights
!
      if(mypid == NROOT) then
       call gauawo(zsi,zgw,NLAT)
       do jlat=1,NLAT
        zgw2(:,jlat)=zgw(jlat)
       enddo
      endif
      call mpscgp(zgw2,gw,1)
!
!     copy input to seamod
!
      nstep    = kstep
      nrestart = krestart
      noutput  = koutput
!
!     read and print namelist and distribute it
!
      if (mypid == NROOT) then
         open(12,file='ocean_namelist',form='formatted')
         read(12,oceanpar)
         write (*,'(/," *******************************************")')
         write (*,'(" * OCEANMOD ",a30," *")') trim(version)
         write (*,'(" *******************************************")')
         write (*,'(" * Namelist OCEANPAR from <ocean_namelist> *")')
         write (*,'(" *******************************************")')
         write(*,oceanpar)
         close(12)
      endif

      call mpbci(ndiag)
      call mpbci(nout)
      call mpbci(nfluko)
      call mpbci(nocean)
      call mpbci(newsurf)
      call mpbci(nprint)
      call mpbci(nprhor)
      call mpbci(ntspd)
      call mpbci(nperpetual_ocean)
      call mpbcr(taunc)
      call mpbcr(vdiffk)
      call mpbcrn(dlayer,NLEV_OCE)
!
      do jlev=1,NLEV_OCE
         ymld(:,jlev) = dlayer(jlev)
      enddo
!
!     compute time step
!
      dtmix = solar_day / real(ntspd)
!
!     compute taunc in s
!
      taunc = solar_day  * taunc
!
      if (nrestart == 0) then ! new start (read start file)
         call mpsurfgp('yls',yls,NHOR,1)
         call mpsurfyear('yclsst',yclsst,NHOR)

!        make sure, that land sea mask values are 0 or 1

         where (yls(:) > 0.5)
            yls(:) = 1.0
         elsewhere
            yls(:) = 0.0
         endwhere
!
!        make clsst >= tfreeze
!
         yclsst(:,:) = MAX(yclsst(:,:),TFREEZE)
!
!        initialize sst
!
         call oceanget

         do jlev=1,NLEV_OCE
            ysst(:,jlev) = yclsst2(:)
         enddo

      else ! restart from restart file <surface.txt>
         if (mypid == NROOT) then
            call get_restart_integer('nstep'   ,nstep   )
            call get_restart_integer('naccuoce',naccuout)
         endif
         call mpbci(nstep)
         call mpbci(naccuout)

         call mpgetgp('yls'    ,yls    ,NHOR,   1)
         call mpgetgp('ysst'   ,ysst   ,NHOR,NLEV_OCE)
         call mpgetgp('yiflux' ,yiflux ,NHOR,   1)
         call mpgetgp('yclsst' ,yclsst ,NHOR,  14)
         call mpgetgp('yheata' ,yheata ,NHOR,   1)
         call mpgetgp('yfssta' ,yfssta ,NHOR,   1)
         call mpgetgp('yifluxa',yifluxa,NHOR,   1)
         call mpgetgp('ydssta' ,ydssta ,NHOR,   1)

         if (newsurf == 1) then ! Read new surface data
            call mpsurfyear('yclsst',yclsst,NHOR)
            yclsst(:,:) = MAX(yclsst(:,:),TFREEZE)
         endif
      endif ! (nrestart == 0)
!
!     read flux correction
!
      if (nfluko == 1) then
         call mpsurfyear('yfsst',yfsst,NHOR)
      endif
!
!     open output file
!
      if (mypid == NROOT) then
       open(31,file='ocean_output',form='unformatted')
      endif
!
!     copy output from seamod
!
      psst(:)=ysst(:,1)
      pmld(:)=ymld(:,1)
      piflux(:)=yiflux(:)
!
      return
      end subroutine oceanini

!     ===================================================================
!     SUBROUTINE oceanstep
!     ===================================================================

      subroutine oceanstep(picec,piced,pheat,pfresh,ptaux,ptauy,pust3   &
     &                    ,psst,pmld,piflux)
      use oceanmod
!
      real :: picec(NHOR),piced(NHOR),pheat(NHOR),pfresh(NHOR)
      real :: ptaux(NHOR),ptauy(NHOR),pust3(NHOR)
      real :: psst(NHOR),pmld(NHOR),piflux(NHOR)
!
!     dbug arrays
!
      real,allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:),zprf5(:)
!
!     set helpful bits
!
      zcmdt=CRHOS*CPS/dtmix
!
!     copy input to seamod
!
      yicec(:)=picec(:)
      yiced(:)=piced(:)
      yheat(:)=pheat(:)
      yfresh(:)=pfresh(:)
      ytaux(:)=ptaux(:)
      ytauy(:)=ptauy(:)
      yust3(:)=pust3(:)
!
!     get annual cycle
!
      call oceanget
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       call mpgagp(zprf1,yicec,1)
       call mpgagp(zprf2,yiced,1)
       call mpgagp(zprf3,yheat,1)
       call mpgagp(zprf4,yclsst2,1)
       if(mypid==NROOT) then
        print*,'in oceanstep: nstep= ',nstep
        print*,' '
        print*,'seaice c and d from ice: ',zprf1(nprhor),zprf2(nprhor)
        print*,'heatflux from ice: ',zprf3(nprhor)
        print*,'climatological sst: ',zprf4(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
      endif
!
!     compute new sst
!
      call mksst
!$TAF STORE ysst = inner1, REC = 1 + nstep2
!
!     a) climatological ocean
!
      if(nocean == 0) then
!
       where(yls(:) < 1.)
!
!     flux correction diagnostics
!
        yfsst2(:)=zcmdt*(yclsst2(:)-ysst(:,1))*ymld(:,1)
!
!     set sst
!
        ysst(:,1)=yclsst2(:)
!
       endwhere
!
!     print dbug information if needed
!
       if (nprint==2) then
        allocate(zprf1(NLON*NLAT))
        call mpgagp(zprf1,yfsst2,1)
        if(mypid==NROOT) then
         print*,'computed flux correction: ',zprf1(nprhor)
        endif
        deallocate(zprf1)
       endif
!
      else
!
!     b) interactive ocean
!
!     get flux correction
!
       if(nfluko == 0) then
        yfsst2(:)=0.
       elseif(nfluko == 1) then
        call getflxco
        call addfc
       elseif(nfluko == 2) then
        call mkfc
        call addfc
       endif
!
      endif
!$TAF STORE yiflux = inner1, REC = 1 + nstep2
!
!     output
!
      yfssta(:)=yfssta(:)+yfsst2(:)
      yheata(:)=yheata(:)+yheat(:)
      yifluxa(:)=yifluxa(:)+yiflux(:)
      ydssta(:)=ydssta(:)+ydsst(:)
      naccuout=naccuout+1
      if(mod(nstep,nout) == 0) then
       yheata(:)=yheata/real(naccuout)
       yfssta(:)=yfssta/real(naccuout)
       yifluxa(:)=yifluxa/real(naccuout)
       ydssta(:)=ydssta/real(naccuout)
       if(noutput > 0) call oceanout
       yfssta(:)=0.
       yheata(:)=0.
       yifluxa(:)=0.
       ydssta(:)=0.
       naccuout=0
      endif
!
!     make flux into ice
!
      call mkiflx
!
!     copy output to seamod
!
      psst(:)=ysst(:,1)
      pmld(:)=ymld(:,1)
      piflux(:)=yiflux(:)
!
!     advance time step
!
      nstep=nstep+1
!
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       call mpgagp(zprf1,psst,1)
       call mpgagp(zprf2,piflux,1)
       if(mypid==NROOT) then
        print*,'final sst to ice: ',zprf1(nprhor)
        print*,'final heatflux to ice: ',zprf2(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
      endif
!
      return
      end subroutine oceanstep

!     ===================================
!     SUBROUTINE OCEANSTOP
!     ===================================

      subroutine oceanstop
      use oceanmod
!
!     write restart file
!
      if (mypid == NROOT) then
         call put_restart_integer('nlev_oce',NLEV_OCE)
         call put_restart_integer('naccuoce',naccuout)
      endif

      call mpputgp('yls'    ,yls    ,NHOR,   1)
      call mpputgp('ysst'   ,ysst   ,NHOR,NLEV_OCE)
      call mpputgp('yiflux' ,yiflux ,NHOR,   1)
      call mpputgp('yclsst' ,yclsst ,NHOR,  14)
      call mpputgp('yheata' ,yheata ,NHOR,   1)
      call mpputgp('yfssta' ,yfssta ,NHOR,   1)
      call mpputgp('yifluxa',yifluxa,NHOR,   1)
      call mpputgp('ydssta' ,ydssta ,NHOR,   1)
!
!     close output file
!
      if (mypid == NROOT) then
       close(31)
      endif
!
      return
      end subroutine oceanstop

!     =====================================================================
!     SUBROUTINE oceanout
!     =====================================================================

      subroutine oceanout
      use oceanmod

      integer :: ih(8)

      call ntomin(nstep,nmin,nhour,nday,nmonth,nyear)

      ih(2) = 0
      ih(3) = nyear*10000 +nmonth*100 +nday
      ih(4) = nhour*100 +nmin
      ih(5) = NLON
      ih(6) = NLAT
      ih(7) = 0
      ih(8) = 0
!
!     heat flux from atm/ice
!
      ih(1) = 901
      call mpwritegph(31,yheata,NHOR,1,ih)
!
!     heat flux into ice
!
      ih(1) = 902
      call mpwritegph(31,yifluxa,NHOR,1,ih)
!
!     flux correction
!
      ih(1) = 903
      call mpwritegph(31,yfssta,NHOR,1,ih)
!
!     flux used for dsst/dt
!
      ih(1) = 904
      call mpwritegph(31,ydssta,NHOR,1,ih)
!
!     sea ice cover
!
      ih(1) = 910
      call mpwritegph(31,yicec,NHOR,1,ih)
!
!     sst
!
      ih(1) = 939
      do jlev=1,NLEV_OCE
       ih(2)=jlev
       call mpwritegph(31,ysst(:,jlev),NHOR,1,ih)
      enddo
      ih(2)=0
!
!     land sea mask
!
      ih(1) = 972
      call mpwritegph(31,yls,NHOR,1,ih)
!
!     clim. sst
!
      ih(1) = 990
      call mpwritegph(31,yclsst2,NHOR,1,ih)
!
      return
      end subroutine oceanout


!====================================================================
!      SUBROUTINE OCEANGET
!====================================================================

      subroutine oceanget
      use oceanmod

!     *********************
!     * get  annual cycle *
!     *********************

!     modified for 14 months version of SST
!     jm2: 0=Dec, 1-12:Months, 13=Jan

      jstep=nstep+1    !advance time step (ts = ts(t)+ dtsdt)
      call ntomin(jstep,nmin,nhour,nday,nmonth,nyear)

      if (nperpetual_ocean > 0) then
         nmonth=(nperpetual_ocean-1)/30.+1
         nday=nperpetual_ocean-nmonth*30
      endif

      jm1 = nmonth
      if (nday > 15) then
         jm2=jm1+1
      else
         jm2=jm1-1
      endif

!     ****************************************
!*    interpolate
!     ****************************************

       zgw2 = abs((nday-1)*1440+nhour*60+nmin-21600)/43200.
       zgw1 = 1.0 - zgw2
       yclsst2(:)  = zgw1*yclsst(:,jm1)+zgw2*yclsst(:,jm2)

       return
       end subroutine oceanget


!======================================================================
!     SUBROUTINE GETFLXCO
!======================================================================

      subroutine getflxco
      use oceanmod
!
!     dbug arrays
!
      real,allocatable :: zprf1(:),zprf2(:),zprf3(:)
!
!     get  annual cycle
!
      jstep=nstep+1    !advance time step (ts = ts(t)+ dtsdt)
      call ntomin(jstep,nmin,nhour,nday,nmonth,nyear)
!
      if(nperpetual_ocean > 0) then
       nmonth=(nperpetual_ocean-1)/30.+1
       nday=nperpetual_ocean-nmonth*30
      endif
!
      jm1=nmonth
      if(nday.gt.15) then
       jm2=jm1+1
       if(jm1.eq.12) jm2=1
      else
       jm2=jm1-1
       if(jm1.eq.1) jm2=12
      endif
!
!     interpolate
!
      zgw2=abs((nday-1)*1440+nhour*60+nmin-21600)/43200.
      zgw1=1.-zgw2
      yfsst2(:)=zgw1*yfsst(:,jm1)+zgw2*yfsst(:,jm2)
!
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       call mpgagp(zprf1,yfsst2,1)
       call mpgagp(zprf2,yfsst(:,jm1),1)
       call mpgagp(zprf3,yfsst(:,jm2),1)
       if(mypid==NROOT) then
        print*,'interpolated flux correction : ',zprf1(nprhor)
        print*,'from: ',zprf2(nprhor),' and ',zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
      endif
!
      return
      end subroutine getflxco

!======================================================================
!     SUBROUTINE MKSST
!======================================================================
!
      subroutine mksst
      use oceanmod
!
!     compute new temperatures
!
      real :: zold(NHOR)
      real :: zflx(NHOR)
!
!     dbug arrays
!
      real,allocatable :: zprf1(:),zprf2(:),zprf3(:)
!
!     set some useful bits
!
      zcpsdt=dtmix/(CRHOS*CPS)
!
!     add surface heat flux
!
      where(yls(:) < 1.)
       ysst(:,1)=ysst(:,1)+yheat(:)*zcpsdt/ymld(:,1)
      endwhere
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       call mpgagp(zprf1,ysst(:,1),1)
       if(mypid==NROOT) then
        print*,'in mksst:'
        print*,' '
        print*,'new sst from heatflx: ',zprf1(nprhor)
       endif
       deallocate(zprf1)
      endif
!
!     comput residual flux going into sea ice (if sst < tfreeze)
!
      yiflux(:)=0.
      where(ysst(:,1) < TFREEZE .and. yls(:) < 1.)
       yiflux(:)=(ysst(:,1)-TFREEZE)*ymld(:,1)/zcpsdt
       ysst(:,1)=TFREEZE
      endwhere
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       call mpgagp(zprf1,ysst(:,1),1)
       call mpgagp(zprf2,yiflux,1)
       if(mypid==NROOT) then
        print*,'new sst from correction sst < tfreeze: ',zprf1(nprhor)
        print*,'residual heat flux into ice: ',zprf2(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
      endif
!
!     make vertical diffusion (new temperatures) if more than 1 layer
!
      where(yls(:) < 1.)
       zold(:)=ysst(:,1)
      endwhere
!
!     vertical diffusion if more tahn 1 layer
!
      if(NLEV_OCE > 1) then
       call vdiffo
!
!     diagnose vdiff flux used for dsst/dt
!
       where(yls(:) < 1.)
        ydsst(:)=(ysst(:,1)-zold(:))*ymld(:,1)/zcpsdt
       endwhere
!
!     use residual flux (into ice) to cool water again (if > tfreeze)
!
       where(ysst(:,1) > TFREEZE .and. yiflux(:) < 0. .and. yls(:) < 1.)
        zflx(:)=(TFREEZE-ysst(:,1))*ymld(:,1)/zcpsdt
        zflx(:)=AMAX1(zflx(:),yiflux(:))
        ysst(:,1)=ysst(:,1)+zflx(:)*zcpsdt/ymld(:,1)
        yiflux(:)=yiflux(:)-zflx(:)
       endwhere
      endif
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       call mpgagp(zprf1,ysst(:,1),1)
       call mpgagp(zprf2,ydsst,1)
       call mpgagp(zprf3,yiflux,1)
       if(mypid==NROOT) then
        print*,'new sst from diffusion: ',zprf1(nprhor)
        print*,'heat from diffusion: ',zprf2(nprhor)
        print*,'new residual heat flux: ',zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
      endif
!
      return
      end subroutine mksst

!======================================================================
!     SUBROUTINE ADDFC
!======================================================================

      subroutine addfc
      use oceanmod
!
!     dbug arrays
!
      real,allocatable :: zprf1(:),zprf2(:)
!
!     add flux correction
!
      zcpsdt=dtmix/(CRHOS*CPS)
!
      where(yls(:) < 1.)
       ysst(:,1)=ysst(:,1)+yfsst2(:)*zcpsdt/ymld(:,1)
      endwhere
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       call mpgagp(zprf1,ysst(:,1),1)
       if(mypid==NROOT) then
        print*,'in addfc: '
        print*,' '
        print*,'new sst from fluko: ',zprf1(nprhor)
       endif
       deallocate(zprf1)
      endif
!
!     comput residual flux going into sea ice
!
!     a) sst < tfreeze
!
      where(ysst(:,1) < TFREEZE .and. yls(:) < 1.)
       yiflux(:)=yiflux(:)+(ysst(:,1)-TFREEZE)*ymld(:,1)/zcpsdt
       ysst(:,1)=TFREEZE
      endwhere
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       call mpgagp(zprf1,ysst(:,1),1)
       call mpgagp(zprf2,yiflux,1)
       if(mypid==NROOT) then
        print*,'new sst from correction sst < tfreeze: ',zprf1(nprhor)
        print*,'residual heat flux into ice: ',zprf2(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
      endif
!
!     b) sst > t climate and ice > 0 (note: min(clsst) = TFREEZE)
!
      where(ysst(:,1) > yclsst2(:) .and. yicec(:) > 0.)
       yiflux(:)=yiflux(:)+(ysst(:,1)-yclsst2(:))*ymld(:,1)/zcpsdt
       ysst(:,1)=yclsst2(:)
      endwhere
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       call mpgagp(zprf1,ysst(:,1),1)
       call mpgagp(zprf2,yiflux,1)
       if(mypid==NROOT) then
        print*,'new sst from correction sst > tclim: ',zprf1(nprhor)
        print*,'residual heat flux into ice: ',zprf2(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
      endif
!
      return
      end subroutine addfc

!======================================================================
!     SUBROUTINE MKIFLX
!======================================================================
!
      subroutine mkiflx
      use oceanmod
      real :: zsum(2)
!
!     debug arrays
!
      real,allocatable :: zprf1(:)
!
!     since iflux could be locally to large a global avergage (ice only)
!     iflux is computed which goes than into the sea ice model
!
      zsum(1)=SUM(yiflux(:)*gw(:),MASK=(yicec(:) > 0.))
      zsum(2)=SUM(gw(:),MASK=(yicec(:) > 0.))
      call mpsumbcr(zsum,2)
      if(zsum(1) /= 0.) then
       where(yicec(:) > 0.) yiflux(:)=zsum(1)/zsum(2)
      endif
!
!     print dbug information if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       call mpgagp(zprf1,yiflux,1)
       if(mypid==NROOT) then
        print*,'in mkiflx: '
        print*,' '
        print*,' global sums 1 and 2: ',zsum(1:2)
        print*,'new residual heat flux into ice: ',zprf1(nprhor)
       endif
       deallocate(zprf1)
      endif
!
      return
      end subroutine mkiflx

!====================================================================
!     SUBROUTINE MKFC
!====================================================================

      subroutine mkfc
      use oceanmod
!
      zcm=CRHOS*CPS
!
      if (taunc > 0.) then
       where(yls(:) < 1.)
        yfsst2(:)=zcm*(yclsst2(:)-ysst(:,1))*ymld(:,1)/taunc
       endwhere
      else
       where(yls(:) < 1.)
        yfsst2(:)=zcm/dtmix*(yclsst2(:)-ysst(:,1))*ymld(:,1)-yheat(:)
       endwhere
      endif
!
      return
      end subroutine mkfc

!     =================
!     SUBROUTINE VDIFFO
!     =================

      subroutine vdiffo
!
!     compute new sst due to vertical diffusion
!
      use oceanmod
      parameter(NLEM=NLEV_OCE-1)
!     note: bounds of local arrays due to compiler, only 1:NLEM used
!
      real :: zk(NHOR,0:NLEV_OCE)  = 0.   ! modfied diffusion coefficient
      real :: ztn(NHOR,0:NLEV_OCE) = 0.   ! new temperature
      real :: zebs(NHOR,0:NLEV_OCE)= 0.   ! work array for back substitution
!
!     modified diffusion coeffizient
!
      do jlev=1,NLEM
       where(yls(:) < 1.)
        zk(:,jlev)=dtmix*vdiffk*2./(ymld(:,jlev+1)+ymld(:,jlev))
       endwhere
      enddo
!$TAF STORE zk = inner1, REC = 1 + nstep2
!
!     semi implizit scheme:
!
!     a) top layer elimination
!
      where(yls(:) < 1.)
       zebs(:,1)=zk(:,1)/(ymld(:,1)+zk(:,1))
       ztn(:,1)=ymld(:,1)*ysst(:,1)/(ymld(:,1)+zk(:,1))
      endwhere
!
!     b) middle layer elimination
!
      if(NLEV_OCE > 2) then
       do jlev=2,NLEV_OCE-1
        jlem=jlev-1
!$taf store zebs(:,jlem) = innerno, rec = (jlev - 1) + nstep2*(NLEV_OCE-2)
        where(yls(:) < 1.)
         zebs(:,jlev)=zk(:,jlev)/(ymld(:,jlev)+zk(:,jlev)               &
     &                           +zk(:,jlem)*(1.-zebs(:,jlem)))
         ztn(:,jlev)=(ysst(:,jlev)*ymld(:,jlev)+zk(:,jlem)*ztn(:,jlem)) &
     &              /(ymld(:,jlev)+zk(:,jlev)                           &
     &               +zk(:,jlem)*(1.-zebs(:,jlem)))
        endwhere
       enddo
      endif
!$TAF STORE zebs = inner1, REC = 1 + nstep2
!
!     c) bottom layer elimination (new temperature)
!
      where(yls(:) < 1.)
       ysst(:,NLEV_OCE)=(ysst(:,NLEV_OCE)*ymld(:,NLEV_OCE)+zk(:,NLEM)*ztn(:,NLEM))  &
     &             /(ymld(:,NLEV_OCE)+zk(:,NLEM)*(1.-zebs(:,NLEM)))
      endwhere
!
!     d) back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       where(yls(:) < 1.)
        ysst(:,jlev)=ztn(:,jlev)+zebs(:,jlev)*ysst(:,jlep)
       endwhere
      enddo
!
      return
      end

!     =================
!     SUBROUTINE GAUAWO
!     =================

      subroutine gauawo(psi,pgw,k)
      real psi(k),pgw(k)
      real (kind=8) c,eps
      parameter(c = (1.0 - (2.0/3.14159265358979)**2)/4.0)
      parameter(eps = 1.0e-14)
      real (kind=8) pa(k),pw(k)
      real (kind=8) pk,pkm1,pkm2,pkmrk,zsp,xz

      kk = k/2
      call bsslzr(pa,kk)

      do is=1,kk
         xz = cos(pa(is)/sqrt((k+0.5)**2+c))
         iter = 0
         zsp  = 1.0
         do while (abs(zsp) > eps .and. iter <= 10)
            pkm1 = xz
            pkm2 = 1.0
            do n = 2 , k
               pk   = ((n+n-1.0)*xz*pkm1 - (n-1.0) * pkm2) / float(n)
               pkm2 = pkm1
               pkm1 = pk
            enddo
            pkm1  = pkm2
            pkmrk = (k * (pkm1 - xz * pk)) / (1.0 - xz*xz)
            zsp   = pk / pkmrk
            xz    = xz - zsp
            iter  = iter + 1
         enddo ! while
         if (iter > 10) then
            print *,' *** no convergence in gauaw ***'
            stop
         endif

         pa(is) = xz
         pw(is) = (2.0 * (1.0 - xz*xz)) / (k*k*pkm1*pkm1)
         pa(k+1-is) = -pa(is)
         pw(k+1-is) =  pw(is)
      enddo

      psi = pa
      pgw = pw
      return
      end

!     ==================
!     SUBROUTINE BSSLZRO
!     ==================

      subroutine bsslzro(pbes,kdim)
      parameter(NDIM=50)
      real (kind=8) pbes(kdim)
      real (kind=8) zbes(NDIM)
      real (kind=8) PI
      parameter(PI = 3.14159265359)
      data zbes        /  2.4048255577,   5.5200781103,                 &
     &    8.6537279129,  11.7915344391,  14.9309177086,  18.0710639679, &
     &   21.2116366299,  24.3524715308,  27.4934791320,  30.6346064684, &
     &   33.7758202136,  36.9170983537,  40.0584257646,  43.1997917132, &
     &   46.3411883717,  49.4826098974,  52.6240518411,  55.7655107550, &
     &   58.9069839261,  62.0484691902,  65.1899648002,  68.3314693299, &
     &   71.4729816036,  74.6145006437,  77.7560256304,  80.8975558711, &
     &   84.0390907769,  87.1806298436,  90.3221726372,  93.4637187819, &
     &   96.6052679510,  99.7468198587, 102.8883742542, 106.0299309165, &
     &  109.1714896498, 112.3130502805, 115.4546126537, 118.5961766309, &
     &  121.7377420880, 124.8793089132, 128.0208770059, 131.1624462752, &
     &  134.3040166383, 137.4455880203, 140.5871603528, 143.7287335737, &
     &  146.8703076258, 150.0118824570, 153.1534580192, 156.2950342685/

      idim = min(kdim,NDIM)
      pbes(1:idim) = zbes(1:idim)
      do j = idim+1 , kdim
         pbes(j) = pbes(j-1) + PI
      enddo
      return
      end
