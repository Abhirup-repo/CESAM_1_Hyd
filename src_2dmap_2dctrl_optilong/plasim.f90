

!-*-f90-*-

! FastOpt:  This file is cpp-preprocessed. Be sure to edit plasim.F90
! FastOpt:  This source has been modified for TAF compliance
!
!     **********************************************************
!     * PUMA 2.0 * Portable University Model of the Atmosphere *
!     **********************************************************

!     *****************
!     * Frank Lunkeit *       University of Hamburg
!     * Edilbert Kirk *      Meteorological Institute
!     *****************   Dept.: Theoretical Meteorology
!                             Head: Klaus Fraedrich

!     **************
!     * Name rules *
!     **************

!     i - local integer
!     j - loop index
!     k - integer dummy parameter
!     N - integer constants

!     g - real gridpoint arrays
!     p - real dummy parameter
!     s - real spectral arrays
!     z - local real
!
! FastOpt: main program moved to external driver routine
! FastOpt 
! FastOpt       program plasim_main
! FastOpt       use pumamod
! FastOpt 
! FastOpt !     ******************************************************************
! FastOpt !     *         Portable University Model of the Atmosphere            *
! FastOpt !     ******************************************************************
! FastOpt 
! FastOpt !     Frank (Larry) Lunkeit & Edilbert Kirk
! FastOpt !     Meteorologisches Institut
! FastOpt !     Universitaet Hamburg
! FastOpt 
! FastOpt !     based on SGCM (Simple Global Circulation Model)
! FastOpt !     by Ian James and James Dodd (April 1993)
! FastOpt !     Dept of Meteorology, University of Reading, UK.
! FastOpt 
! FastOpt !     With relaxation to basic restoration temerature state
! FastOpt !     The state is set in SETZT(EX) for initial runs,
! FastOpt !     on a time scale which can vary with the level
! FastOpt !     and which is entered in the namelist in restim(NLEV).
! FastOpt !     There is a linear drag which can again vary with level,
! FastOpt !     the time scale is entered for each level in days in tfrc(NLEV).
! FastOpt !     the diffusion is del-ndel with a time scale for diffusion of
! FastOpt !     tdiss days on every variable at the truncation wavenumber.
! FastOpt !
! FastOpt !     set version identifier
! FastOpt !
! FastOpt 
! FastOpt       pumaversion='2.2.1 (04.08.2005 by Larry)'
! FastOpt 
! FastOpt       call mpstart
! FastOpt       call prolog
! FastOpt       call master
! FastOpt       call epilog
! FastOpt       call mpstop
! FastOpt 
! FastOpt       stop
! FastOpt       end

!     =================
!     SUBROUTINE PROLOG
!     =================

      subroutine prolog
      use pumamod
      use observation, only : nudgini
      logical :: lrestart
      integer :: itime(8)
      integer :: dummy_nstep
!
      pumaversion='mit_coupling-svn'//"_FastOpt" ! LARRY

!     ************************************************************
!     * Initializations that cannot be run on parallel processes *
!     ************************************************************

      if (mypid == NROOT) then
         call date_and_time(values=itime)
         tmstart = itime(3) * 86400.0 + itime(5) * 3600.0 &
                 + itime(6) * 60.0    + itime(7) + itime(8) * 0.001

         write(*,'(/," ****************************************************")')
         write(*,'(" * PLANET SIMULATOR ",a31," *")') trim(pumaversion)
         write(*,'(" ****************************************************")')
         write(*,'(" * NTRU =",i4,"  NLEV =",i4,"  NLON = ",i4,"   NLAT =",i4," *")') &
            NTRU,NLEV,NLON,NLAT
         write(*,'(" ****************************************************")')
         if (NPRO > 1) then
           write(*,'(/," ****************************************************")')
           do jpro = 1 , NPRO
              write (*,'(" * CPU",i4,1x,a40," *")') jpro-1,ympname(jpro)
           enddo
           write(*,'(" ****************************************************")')
         endif
         call restart_ini(lrestart,'plasim_restart')
         if (lrestart) then
            nrestart = 1
            nkits    = 0
         endif
         call surface_ini           ! Read <surface.txt>
         call gauaw(si,gw,NLAT)     ! Gaussian abscissas and weights
         call inilat                ! Set latitudinal arrays
         call readnl                ! Open and read <puma_namelist>
         call initpm                ! Several initializations
         call initsi                ! Initialize semi implicit scheme
         call guistart              ! Initialize GUI
         call outini                ! Open output file <puma_output>

      if (makenudge==1) then
         call nudgini
      endif
      endif ! (mypid == NROOT)

!     ***********************
!     * broadcast & scatter *
!     ***********************

      call mpscrn(si  ,NLPP)  ! sine of latitude
      call mpscrn(gw  ,NLPP)  ! gaussian weights
      call mpscrn(gaussw ,NLPP) ! gaussian weights w/o *=1/csr
      call mpscrn(csq ,NLPP)  ! cosine squared of latitude

      do jlat = 1 , NLPP
         cola(jlat) = sqrt(csq(jlat))
         rcsq(1+(jlat-1)*NLON:jlat*NLON) = 1.0 / csq(jlat)
      enddo

      call mpbci(mars    ) ! switch for Mars parameter
      call mpbci(ntpal   ) ! color pallet for T

      call mpbci(kick    ) ! add noise for kick > 0
      call mpbci(noutput ) ! write data interval
      call mpbci(nafter  ) ! write data interval
      call mpbci(nwpd    ) ! number of writes per day
      call mpbci(ncoeff  ) ! number of modes to print
      call mpbci(ndiag   ) ! write diagnostics interval
      call mpbci(nexp    ) ! experiment number
      call mpbci(nexper  ) ! 1: dipole experiment
      call mpbci(ngui    ) ! GUI on (1) or off (0)
      call mpbci(nkits   ) ! number of initial timesteps
      call mpbci(nrestart) ! 1: read restart file 0: initial run

      call mpbci(nstep   ) ! current timestep
      call mpbci(ntspd   ) ! number of timesteps per day

      call mpbci(mpstep)   ! minutes per timestep
      call mpbci(n_days_per_month)
      call mpbci(n_days_per_year)
      call mpbci(n_run_steps)
      call mpbci(n_run_days)
      call mpbci(n_run_months)
      call mpbci(n_run_years)
      call mpbci(n_start_step)
      call mpbci(n_start_year)
      call mpbci(n_start_month)
      call mpbci(N_RESET_YEAR)

      call mpbci(nflux   ) !
      call mpbci(nrad    ) !
      call mpbci(neqsig  ) ! switch for equidistant sigma levels
      call mpbci(nhdiff  ) ! critical wavenumber for horizonal diffusion

      call mpbci(nprint    ) ! switch for extensive diagnostic pintout (dbug)
      call mpbci(nprhor    ) ! grid point to be printed (dbug)
      call mpbci(ndiaggp   ) ! switch for franks grid point diagnostics
      call mpbci(ndiagsp   ) ! switch for franks spectral diagnostics
      call mpbci(ndiagcf   ) ! switch for cloud forcing diagnostics
      call mpbci(ndiaggp3d ) ! no of 3d gp diagnostic arrays
      call mpbci(ndiaggp2d ) ! no of 2d gp diagnostic arrays
      call mpbci(ndiagsp3d ) ! no of 3d sp diagnostic arrays
      call mpbci(ndiagsp2d ) ! no of 2d sp diagnostic arrays
      call mpbci(ntime     ) ! switch to activate time consuming estimate
      call mpbci(nperpetual) ! day of perpetual integration

      call mpbcr(acpd    )   ! Specific heat for dry air
      call mpbcr(adv     )
      call mpbcr(akap    )
      call mpbcr(alr     )
      call mpbcr(cv      )
      call mpbcr(ct      )
      call mpbcr(dtep    )
      call mpbcr(dtns    )
      call mpbcr(dtrop   )
      call mpbcr(dttrp   )
      call mpbcr(ga      )
      call mpbcr(gascon  )
      call mpbcr(tgr     )
      call mpbcr(plarad  )
      call mpbcr(pnu     )
      call mpbcr(pnu21   )
      call mpbcr(psurf   )
      call mpbcr(ra1     )
      call mpbcr(ra2     )
      call mpbcr(ra4     )
      call mpbcr(rdbrv   )
      call mpbcr(ww      )
      call mpbcr(solar_day)
      call mpbcr(siderial_day)

      call mpbcin(ndel  ,NLEV) ! ndel
      call mpbcin(ndl   ,NLEV)

      call mpbcrn(tdissd,NLEV)
      call mpbcrn(tdissz,NLEV)
      call mpbcrn(tdisst,NLEV)
      call mpbcrn(tdissq,NLEV)
      call mpbcrn(dnudg, NLEV)
      call mpbcrn(znudg, NLEV)
      call mpbcrn(qnudg, NLEV)
      call mpbcrn(tnudg, NLEV)
      call mpbcr(pnudg)
      call mpbcr(smnudg)
      call mpbcrn(damp  ,NLEV)
      call mpbcrn(dsigma,NLEV)
      call mpbcrn(rdsig ,NLEV)
      call mpbcrn(restim,NLEV)
      call mpbcrn(sigma ,NLEV)
      call mpbcrn(sigmah,NLEV)
      call mpbcrn(t0    ,NLEV)
      call mpbcrn(t01s2 ,NLEV)
      call mpbcrn(tfrc  ,NLEV)
      call mpbcrn(tkp   ,NLEV)

      call mpbcrn(c     ,NLSQ)
      call mpbcrn(g     ,NLSQ)
      call mpbcrn(tau   ,NLSQ)

      call mpscin(nindex,NSPP)
      call mpscsp(sak,sakpp,NLEV)
      call mpbcrn(sigh  ,NLEV)

      call mpbci(ncost)    ! fastopt cost function parameter
      call mpbci(n_spinup) ! fastopt cost function parameter

      call mpbci(makenudge) ! A.B switch for nudging 
      call mpbci(chyear)

!     Copy some calendar variables to calmod

      call calini(n_days_per_month,n_days_per_year,n_start_step,ntspd,solar_day)

!     allocate additional diagnostic arrays, if switched on
!

      if(ndiaggp2d > 0) then
       allocate(dgp2d(NHOR,ndiaggp2d))
       dgp2d(:,:)=0.
      end if
      if(ndiaggp3d > 0) then
       allocate(dgp3d(NHOR,NLEV,ndiaggp3d))
       dgp3d(:,:,:)=0.
      end if
      if(ndiagsp2d > 0) then
       allocate(dsp2d(NESP,ndiagsp2d))
       dsp2d(:,:)=0.
      end if
      if(ndiagsp3d > 0) then
       allocate(dsp3d(NESP,NLEV,ndiagsp3d))
       dsp3d(:,:,:)=0.
      end if
      if(ndiagcf > 0) then
       allocate(dclforc(NHOR,7))
       dclforc(:,:)=0.
      end if

      call legini

      if (nrestart > 0) then
         call read_atmos_restart
      else
         call initfd
      endif

      
!
!*    FastOpt: compute length of main loop (moved here from master)
!

!FastOpt!     Use either month countdown (n_run_years * 12 + n_run_months)
!FastOpt!     or step-countdown (for debugging purposes)
!FastOpt      
!FastOpt      mocd = n_run_months                     ! month countdown
!FastOpt      nscd = n_run_days * ntspd + n_run_steps ! step countdown

      if (nrestart == 0) nstep = n_start_step ! timestep since 01-01-0001
! guokun set the initial date to 2000.01.01.00
! this is hard code here
      if (N_RESET_YEAR == 2100 .and. n_days_per_year == 365) then
         nstep=23010300
         print*, 'simulation time is reset to ', nstep
       endif
      call updatim(nstep)  ! set date & time array ndatim
      if(mypid == NROOT) write(*,'(a,i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2," UTC")') &
           'simulations starting ',ndatim(1),ndatim(2),ndatim(3),ndatim(4),ndatim(5)
! FastOpt      do while (mocd > 0 .or. nscd > 0)  ! main loop
      nstep0 = nstep ! FastOpt, remember start timestep

      if (n_days_per_year == 365) then ! added, FastOpt
! cal2step routines robust against days > 30
         call cal2step(nstepe,ntspd,ndatim(1)+(ndatim(2)+n_run_months)/12,mod(ndatim(2)+n_run_months,12), &
              & ndatim(3)+n_run_days,ndatim(4),ndatim(5)) ! added, FastOpt
      else ! added, FastOpt
         call cal2step30(nstepe,ntspd,ndatim(1)+(ndatim(2)+n_run_months)/12,mod(ndatim(2)+n_run_months,12), &
              & ndatim(3)+n_run_days,ndatim(4),ndatim(5)) ! added, FastOpt
      end if ! added, FastOpt

      nstepe = nstepe + n_run_steps + n_spinup ! added, FastOpt

      call updatim(nstepe)  ! set date & time array ndatim
      if(mypid == NROOT) write(*,'(a,i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2," UTC")') &
           'simulations stopping ',ndatim(1),ndatim(2),ndatim(3),ndatim(4),ndatim(5)

      call updatim(nstep)  ! set date & time array ndatim

      dprl(:)=0. ! added, FastOpt
      dprc(:)=0. ! added, FastOpt
      dprs(:)=0. ! added, FastOpt
      dftu(:,:)  = 0.0 ! added, FastOpt
      if (ndiag==0) then ! FastOpt
         ndiag = nstepe
         print'("ndiag set to ",i8)',ndiag ! FastOpt
         call mpbci(ndiag   ) ! write diagnostics interval
      end if ! FastOpt

!
!*    initialize miscellaneous additional parameterization
!     which are included in *miscmod*
!

      call miscini

!
!*    initialize surface fluxes and vertical diffusion
!

      call fluxini

!
!*    initialize radiation
!

      call radini

!
!*    initialize convective and large scale rain and clouds
!

      call rainini

!
!*    initialize surface parameterizations or models
!

      call surfini

!
!*    reset psurf according to orography
!

      if (mypid == NROOT) then
         write (*,'(" Surface pressure with no topography = ",f10.2," [hPa]")')&
            psurf * 0.01
         zmeanoro = so(1) * cv * cv / sqrt(2.0)
         zdlnp    = zmeanoro / gascon / tgr
         psurf    = EXP(LOG(psurf)-zdlnp)
         write (*,'(" Mean of topographic height          = ",f10.2," [m]")')&
            zmeanoro / ga
         write (*,'(" Mean of surface pressure            = ",f10.2," [hPa]")')&
            psurf * 0.01
      end if
      call mpbcr(psurf)

!
!*    broadcast and scatter
!

      call mpbcrn(sp,NESP)
      call mpbcrn(sd,NESP*NLEV)
      call mpbcrn(st,NESP*NLEV)
      call mpbcrn(sz,NESP*NLEV)
      call mpbcrn(sq,NESP*NLEV)

      call mpscsp(sd,sdp,NLEV)
      call mpscsp(st,stp,NLEV)
      call mpscsp(sz,szp,NLEV)
      call mpscsp(sq,sqp,NLEV)
      call mpscsp(sr,srp,NLEV)
      call mpscsp(sp,spp,1)
      call mpscsp(so,sop,1)

!
!*    close the namelist file
!

      if(mypid==NROOT) close(11)

      if (mypid == NROOT) call print_planets
!
!*    start time consuming calculations
!

      if(ntime==1) then
       call mksecond(zsec,0.)
       time0=zsec
      endif

      if (mypid == NROOT) then
         csu  = 0.
         csv  = 0.
         cst  = 0.
         span = 0.
      end if

      return
      end


!     ========================
!     SUBROUTINE PRINT_PLANETS
!     ========================

      subroutine print_planets
      use pumamod

      write (*,4000)
      write (*,1000)
      write (*,1100) 'This simulatiomn is for ',yplanet
      write (*,1000)
      write (*,2000) 'Parameter','Units','Earth','Mars'
      write (*,1000)
      write (*,3000) 'Mass'             ,'[10^24 kg]' ,   5.9736,   0.6419
      write (*,3000) 'Volume'           ,'[10^10 km3]', 108.321 ,  16.318 
      write (*,3000) 'Equatorial radius','[km]'       ,6378.0   ,3393.0
      write (*,3000) 'Polar radius'     ,'[km]'       ,6356.0   ,3373.0
      write (*,3000) 'Mean radius'      ,'[km]'       ,6371.0   ,3390.0
      write (*,3000) 'Ellipticity'      ,' '        ,   0.0034,   0.0065
      write (*,3000) 'Mean density'     ,'[kg/m3]'    ,5520.0   ,3933.0   
      write (*,3000) 'Surface gravity'  ,'[m/s2]'     ,GA_EARTH ,GA_MARS
      write (*,3000) 'Bond albedo'      ,' '        ,   0.385 ,   0.16
      write (*,3000) 'Solar irradiance' ,'[W/m2]' ,SOLCON_EARTH,SOLCON_MARS
      write (*,3000) 'Black-body temperature','[K]'   , 247.3   , 216.6
      write (*,3000) 'Topographic range'     ,'[km]'  ,  20.0   ,  36.0
      write (*,3000) 'Sidereal orbit period' ,'[days]', 365.256 , 686.980
      write (*,3000) 'Sidereal rotation period','[hrs]',23.9345 , 24.6229
      write (*,3000) 'Equatorial inclination'  ,'[deg]',23.44   , 23.98  
      write (*,3000) 'Perihelion'       ,'[10^6 km]'  , 147.1   , 206.6
      write (*,3000) 'Aphelion'         ,'[10^6 km]'  , 152.1   , 249.2
      write (*,3000) 'Orbit eccentricity'    ,' '   , 0.0167  , 0.0934
      write (*,1000)
      write (*,4000)
      return

 1000 format(60('*'))
 1100 format('* ',a24,1x,a31,' *')
 2000 format('* ',a24,1x,a11,2a10,' *')
 3000 format('* ',a24,1x,a11,2f10.4,' *')
 4000 format(/)
      end

!     ======================
!     SUBROUTINE INIT_MASTER
!     ======================

      subroutine init_master
      use pumamod
      use icemod, only : cheat, xmld, xsndch, xflxice2,xheat,xicec,xiced,xoflux,xprs,xqmelt,xsnow,xsst,xts ! FastOpt
      use seamod, only : clhdta,clhfla,clwfla,cpmea,croffa,cshdta,cshfla,cswfla,ctauxa,ctauya,cust3a,cheata, &
           &           cprsa, cicec, ciced, cmld, csmelt, csndch, csnow, csst, cts, dqs  ! FastOpt
      use seamod, only : naccua ! FastOpt
      use seamod, only : dtt=>dts
      use landmod, only : nlsoil, driver, dsnowz, dsoilt, dtclim, dtsm, dwclim ! FastOpt
      use landmod, only : dtu=>dts,dqs_land=>dqs
      use vegmod, only : valb_max,vz0_min,dvalb,dvlai,dvrhs,dvsoil,dvveg,dvz0,nbiomass,zfw,zveg ! FastOpt
      use radmod, only : gmu0,dqo3,gdist2 ! FastOpt
      use rainmod, only : icclev, icctot ! FastOpt
      use oceanmod, only : ysst,nlev_oce ! FastOpt
      use oceanmod, only : nstep_ocean=>nstep ! FastOpt
      use icemod, only : naccuo ! FastOpt
      use observation, only : genav,addcost,writeobs,nobs,cost
      implicit none ! FastOpt
      integer :: ikits, jkits, nscd, iyea, imon

! FastOpt: all loop-length computations moved to initialisation

!     ***************************
!     * short initial timesteps *
!     ***************************

! do initialization only 

      ikits = nkits

      do jkits=1,ikits
         delt  = (TWOPI/ntspd) / (2**nkits)
         delt2 = delt + delt
         if (mypid == NROOT) then
            write (*,*) 'Initial timestep ',jkits,'   delt = ',delt
         endif
         call gridpointa
         call makebm
         call spectrala
         call gridpointd
         call spectrald
         nkits = nkits - 1
       enddo

!     ****************************************************************
!     * The scaling factor "ww" is derived from the rotation "omega" *
!     * with 1 planetary rotation per siderial day (2 Pi)            *
!     * In order to compute "ntspd" timesteps per solar day, we must *
!     * scale the timestep "delt" with (solar_day / siderial_day)    *
!     ****************************************************************

       delt  = (TWOPI * solar_day) / (ntspd * siderial_day)
       delt2 = delt + delt
       call makebm

     end subroutine init_master

     
!     =================
!     SUBROUTINE MASTER
!     =================

     subroutine master(ksteps,kstepo,ikey_dynamics)
     use pumamod

     use icemod, only: naccuo,xfresh,xswfl
     use icemod, only: xtaux,xtauy
     use icemod, only : cheat,xmld,xsndch,xflxice2,xheat,xicec,xiced,xoflux,xprs,xqmelt,xsnow,xsst,xts ! FastOpt
     use seamod, only : cheata,clhdta,clhfla,clwfla,cpmea,croffa,cshdta,cshfla,cswfla,ctauxa,ctauya,cust3a,cheata,&
          &           cprsa,cicec,ciced,cmld,csmelt,csndch,csnow,csst,cts  ! FastOpt
     use seamod, only : naccua ! FastOpt
     use seamod, only : dts_sea=>dts
     use seamod, only : dqs_sea=>dqs
     use landmod, only: devap,dalb
     use landmod, only : nlsoil,driver,dsnowz,dsoilt,dtclim,dtsm,dwclim ! FastOpt
     use landmod, only : dts_land=>dts,dqs_land=>dqs
      use vegmod, only : valb_max,vz0_min,dvalb,dvlai,dvrhs,dvsoil,dvveg,dvz0,nbiomass,zfw,zveg ! FastOpt
     use radmod, only : gmu0,dqo3,gdist2 ! FastOpt
     use rainmod, only : icclev,icctot ! FastOpt
     use oceanmod, only : ysst,nlev_oce ! FastOpt
     use oceanmod, only : nstep_ocean=>nstep ! FastOpt
     use icemod, only : naccuo,xtfreeze,xroff ! FastOpt
     use observation, only : genav,addcost,writeobs,nobs,cost ! FastOpt
     implicit none ! FastOpt
     integer :: ikits, jkits, nscd, iyea, imon
     integer :: ksteps,jstep   ! LARRY
     integer :: kstepo ! step of ocean loop, FastOpt
     integer :: nstep_inner, ikey_dynamics

! FastOpt: all loop-length computations moved to initialisation

     do nstep_inner = 0, ksteps -1 ! added, FastOpt - month/day countdown included
        nstep2 = (ikey_dynamics-1)*ksteps + nstep_inner ! tape counter
        nstep  = nstep0 + (kstepo-1)*ksteps + nstep_inner ! time counter
!        if (mypid == NROOT) then ! debug, FastOpt
!           call updatim(nstep) ! debug, FastOpt
!           if(ncost.ne.0) then
!              print '(a,6i8)','main loop exec no.',nstep_inner,nstep2,nstep-nstep0-n_spinup,nstep, &
!                   mod((nstep-n_spinup-nstep0+1),ncost), &
!                   mod(nstep,nafter) ! debug, FastOpt
!           else
!              print '(a,6i8)','main loop exec no.',nstep_inner,nstep2,nstep-nstep0-n_spinup,nstep, &
!                   mod(nstep,nafter) ! debug, FastOpt
!           end if
!           write(*,'(a,i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2," UTC")') &
!                'simulating ',ndatim(1),ndatim(2),ndatim(3),ndatim(4),ndatim(5) ! debug, FastOpt
!        end if ! debug, FastOpt

!$taf STORE &
!$taf& dcsoil,dcveg,dforest,dicec,dnogrow,dql,dqt,dres,drhs,dsnow,dsnowz,dsoilt,dtclim,dtsm,dvlai,dvrhs,dvsoil,  &
!$taf& dvveg,dvz0,dwatc,dwclim,dwmax,nbiomass,xicec,xiced,ysst,zfw,zveg ,xsst,dlai&
!$taf& = inner1, REC = 1 + nstep2
!        ************************************************************
!        * calculation of non-linear quantities in grid point space *
!        ************************************************************
        call gridpointa
!$TAF STORE sqt = inner1, REC = 1 + nstep2
!        ******************************
!        * adiabatic part of timestep *
!        ******************************
        call spectrala
!!$taf store sdp,sqp,stp,szp,naccuo,nstep_ocean,xtfreeze = inner1, REC = 1 + nstep2
!$taf store sd,sdm,sdp,sp,spm,spp,sq,sqm,sqp,st,stm,stp,sz,szm,szp,naccuo,nstep_ocean,xtfreeze = inner1, REC = 1 + nstep2

!        *****************************
!        * diabatic part of timestep *
!        *****************************
        call gridpointd
!$taf store dp,dq,dqt,dt,sdt,sqt,stt,szt = inner1, REC = 1 + nstep2

        if (mypid == NROOT) then
           if (mod(nstep,nafter) == 0 .and. noutput > 0) then
              call outsp
           endif
           if (mod(nstep,ndiag) == 0) then
              call diag
           elseif (ngui > 0) then
              if (mod(nstep,ngui) == 0) call diag
           endif
        endif
        
        if (ngui > 0) call guistep_plasim
        call spectrald
!$taf store sdp,szp,stp,sqp,dp = inner1, REC = 1 + nstep2
        if (ncost.ne.0) call genav ! FastOpt: preparing averages for cost function
                  
        if (mod(nstep,nafter) == 0) then
           call outaccu
           if(noutput > 0) then
              call outgp
              if(ndiaggp3d+ndiaggp2d+ndiagsp3d+ndiagsp2d+ndiagcf>0)&
                   call outdiag
           endif
        else
           call outaccu
        end if

        if (ncost.ne.0) then
           if (mod(nstep-n_spinup-nstep0+1,ncost) == 0 .and. (nstep-nstep0) .ge. n_spinup ) then

              call addcost((nstep-n_spinup-nstep0+1)/ncost) ! FastOpt: adding contribution to cost function



           end if
        end if
         
        if (mod(nstep,nafter) == 0) then
           call outreset ! moved here, FastOpt
        endif

        iyea  = ndatim(1)    ! current year
        imon  = ndatim(2)    ! current month
        call updatim(nstep+1)  ! set date & time array ndatim
        if (imon /= ndatim(2)) then
           mocd = mocd - 1 ! next month
           if (mypid == NROOT) then
              write (*,"('Completed month ',I2.2,'-',I4.4)") imon,iyea  
           endif
        endif
        call stability_check ! aborts if model tends to explode

! move to higher level         end if ! .not.lunstable, FastOpt
     enddo ! month countdown
     nstep  = nstep0 + kstepo*ksteps  ! added, FastOpt
       
     end subroutine master


!     =================
!     SUBROUTINE EPILOG
!     =================

      subroutine epilog
      use pumamod
      integer :: itime(8)

!
!     deallocate additional diagnostic arrays, if switched on
!

      if(ndiaggp2d > 0) deallocate(dgp2d)
      if(ndiagsp2d > 0) deallocate(dsp2d)
      if(ndiaggp3d > 0) deallocate(dgp3d)
      if(ndiagsp3d > 0) deallocate(dsp3d)
      if(ndiagcf   > 0) deallocate(dclforc)
!
!     close output file
!
      if (mypid == NROOT) close(40)
     if (mypid == NROOT) close(50) ! file for nudging
!
!     write restart file
!

      if (mypid == NROOT) then
         call restart_prepare('plasim_status')
         call put_restart_integer('nstep'   ,nstep   )
         call put_restart_integer('naccuout',naccuout)
         call put_restart_integer('nlat'    ,NLAT    )
         call put_restart_integer('nlon'    ,NLON    )
         call put_restart_integer('nlev'    ,NLEV    )
         call put_restart_integer('nrsp'    ,NRSP    )

         call put_restart_array('sz',sz,NRSP,NESP,NLEV)
         call put_restart_array('sd',sd,NRSP,NESP,NLEV)
         call put_restart_array('st',st,NRSP,NESP,NLEV)
         call put_restart_array('sq',sq,NRSP,NESP,NLEV)
         call put_restart_array('sr',sr,NRSP,NESP,NLEV)
         call put_restart_array('sp',sp,NRSP,NESP,   1)
         call put_restart_array('so',so,NRSP,NESP,   1)
      endif

      call mpputsp('szm',szm,NSPP,NLEV)
      call mpputsp('sdm',sdm,NSPP,NLEV)
      call mpputsp('stm',stm,NSPP,NLEV)
      call mpputsp('sqm',sqm,NSPP,NLEV)
      call mpputsp('spm',spm,NSPP,   1)
!
!     gridpoint restart
!
      call mpputgp('dls'    ,dls    ,NHOR,1)
      call mpputgp('drhs'   ,drhs   ,NHOR,1)
      call mpputgp('dalb'   ,dalb   ,NHOR,1)
      call mpputgp('dz0'    ,dz0    ,NHOR,1)
      call mpputgp('dtaux'  ,dtaux  ,NHOR,1) ! added
      call mpputgp('dtauy'  ,dtauy  ,NHOR,1) ! added
      call mpputgp('dicec'  ,dicec  ,NHOR,1)
      call mpputgp('diced'  ,diced  ,NHOR,1)
      call mpputgp('dwatc'  ,dwatc  ,NHOR,1)
      call mpputgp('drunoff',drunoff,NHOR,1)
      call mpputgp('dveg'   ,dveg   ,NHOR,1)
      call mpputgp('dforest',dforest,NHOR,1)
      call mpputgp('dust3'  ,dust3  ,NHOR,1)
      call mpputgp('dcc'    ,dcc    ,NHOR,NLEV)
      call mpputgp('dql'    ,dql    ,NHOR,NLEV)
      call mpputgp('dqsat'  ,dqsat  ,NHOR,NLEV)
      call mpputgp('dt'     ,dt(1,NLEP),NHOR,1)
      call mpputgp('dq'     ,dq(1,NLEP),NHOR,1)
!
!     accumulated diagnostics
!
      call mpputgp('aprl'  ,aprl  ,NHOR,1)
      call mpputgp('aprc'  ,aprc  ,NHOR,1)
      call mpputgp('aprs'  ,aprs  ,NHOR,1)
      call mpputgp('aevap' ,aevap ,NHOR,1)
      call mpputgp('ashfl' ,ashfl ,NHOR,1)
      call mpputgp('alhfl' ,alhfl ,NHOR,1)
      call mpputgp('aroff' ,aroff ,NHOR,1)
      call mpputgp('asmelt',asmelt,NHOR,1)
      call mpputgp('asndch',asndch,NHOR,1)
      call mpputgp('acc'   ,acc   ,NHOR,1)
      call mpputgp('assol' ,assol ,NHOR,1)
      call mpputgp('asthr' ,asthr ,NHOR,1)
      call mpputgp('atsol' ,atsol ,NHOR,1)
      call mpputgp('atthr' ,atthr ,NHOR,1)
      call mpputgp('ataux' ,ataux ,NHOR,1)
      call mpputgp('atauy' ,atauy ,NHOR,1)
      call mpputgp('atsolu',atsolu,NHOR,1)
      call mpputgp('assolu',assolu,NHOR,1)
      call mpputgp('asthru',asthru,NHOR,1)
      call mpputgp('aqvi'  ,aqvi  ,NHOR,1)
      call mpputgp('atsa'  ,atsa  ,NHOR,1)
      call mpputgp('ats0'  ,ats0  ,NHOR,1)

!
!*    finish graphical user interface
!

      call guistop

!
!*    finish miscellaneous additional parameterizations
!

      call miscstop

!
!*    finish surface fluxes and vertical diffusion
!

      call fluxstop

!
!*    finish radiation
!

      call radstop

!
!*    finish large scale and convective rain and clouds
!

      call rainstop

!
!*    finish surface parameterizations or models
!

      call surfstop


      if (mypid == NROOT) then
         call restart_stop

         print 1020
         print 1030,nstep
         print 1020
      endif

!
!     time consumption
!

      if (mypid == NROOT) then
         call date_and_time(values=itime)
         tmrun = itime(3) * 86400.0 + itime(5) * 3600.0 &
               + itime(6) * 60.0 + itime(7) + itime(8) * 0.001 - tmstart
         ims = nint(tmrun)
         if (ims > 0) then
! FastOpt            isteps = nstep - n_start_step
            isteps = max(1,nstep - n_start_step) ! FastOpt
            zsps = tmrun / isteps ! seconds per step
            zmpy = zsps * ntspd / 60.0 * n_days_per_year
            zypd = (24.0 * 60.0 / zmpy)
            write(*,'(/," **********************************")')
            write(*,'(" * Wallclock time    : ",i6," sec *")') ims
            write(*,'(" * Minutes per year  : ",i6,5x,"*")') nint(zmpy)
            write(*,'(" * Sim years per day : ",i6,5x,"*")') nint(zypd)
            write(*,'(" **********************************")')
         endif
      endif

!     if(ntime==1 .and. mypid==NROOT) then
!      call mksecond(zsec,time0)
!      print*,'****************************************************'
!      print*,' total CPU usage (ROOT process only): ',zsec,'s'
!      print*,'****************************************************'
!     endif
!
      return
 1030 format(' *   Run finished at step ',i16,' *')
 1020 format(' *******************************************')
      end


!     =============================
!     SUBROUTINE READ_ATMOS_RESTART
!     =============================

      subroutine read_atmos_restart
      use pumamod

!     read scalars and full spectral arrays

      if (mypid == NROOT) then
         call get_restart_integer('nstep'   ,nstep)
         call get_restart_integer('naccuout',naccuout)
         call get_restart_array('sz',sz,NRSP,NESP,NLEV)
         call get_restart_array('sd',sd,NRSP,NESP,NLEV)
         call get_restart_array('st',st,NRSP,NESP,NLEV)
         call get_restart_array('sq',sq,NRSP,NESP,NLEV)
         call get_restart_array('sr',sr,NRSP,NESP,NLEV)
         call get_restart_array('sp',sp,NRSP,NESP,   1)
         call get_restart_array('so',so,NRSP,NESP,   1)
      endif

      call mpbci(nstep)     ! broadcast current timestep
      call mpbci(naccuout)  ! broadcast accumulation timestep for diagnostics

!     read and scatter spectral arrays

      call mpgetsp('szm',szm,NSPP,NLEV)
      call mpgetsp('sdm',sdm,NSPP,NLEV)
      call mpgetsp('stm',stm,NSPP,NLEV)
      call mpgetsp('sqm',sqm,NSPP,NLEV)
      call mpgetsp('spm',spm,NSPP,   1)

!     read and scatter surface grids

      call mpgetgp('dls'    ,dls    ,NHOR,   1)
      call mpgetgp('drhs'   ,drhs   ,NHOR,   1)
      call mpgetgp('dalb'   ,dalb   ,NHOR,   1)
      call mpgetgp('dz0'    ,dz0    ,NHOR,   1)
      call mpgetgp('dtaux'  ,dtaux  ,NHOR,1) ! added
      call mpgetgp('dtauy'  ,dtauy  ,NHOR,1) ! added
      call mpgetgp('dicec'  ,dicec  ,NHOR,   1)
      call mpgetgp('diced'  ,diced  ,NHOR,   1)
      call mpgetgp('dwatc'  ,dwatc  ,NHOR,   1)
      call mpgetgp('drunoff',drunoff,NHOR,   1)
      call mpgetgp('dveg'   ,dveg   ,NHOR,   1)
      call mpgetgp('dforest',dforest,NHOR,   1)
      call mpgetgp('dust3'  ,dust3  ,NHOR,   1)
      call mpgetgp('dcc'    ,dcc    ,NHOR,NLEV)
      call mpgetgp('dql'    ,dql    ,NHOR,NLEV)
      call mpgetgp('dqsat'  ,dqsat  ,NHOR,NLEV)
      call mpgetgp('dt'     ,dt(1,NLEP),NHOR,1)
      call mpgetgp('dq'     ,dq(1,NLEP),NHOR,1)

!     read and scatter accumulated diagnostics

      call mpgetgp('aprl'  ,aprl  ,NHOR,1)
      call mpgetgp('aprc'  ,aprc  ,NHOR,1)
      call mpgetgp('aprs'  ,aprs  ,NHOR,1)
      call mpgetgp('aevap' ,aevap ,NHOR,1)
      call mpgetgp('ashfl' ,ashfl ,NHOR,1)
      call mpgetgp('alhfl' ,alhfl ,NHOR,1)
      call mpgetgp('aroff' ,aroff ,NHOR,1)
      call mpgetgp('asmelt',asmelt,NHOR,1)
      call mpgetgp('asndch',asndch,NHOR,1)
      call mpgetgp('acc'   ,acc   ,NHOR,1)
      call mpgetgp('assol' ,assol ,NHOR,1)
      call mpgetgp('asthr' ,asthr ,NHOR,1)
      call mpgetgp('atsol' ,atsol ,NHOR,1)
      call mpgetgp('atthr' ,atthr ,NHOR,1)
      call mpgetgp('ataux' ,ataux ,NHOR,1)
      call mpgetgp('atauy' ,atauy ,NHOR,1)
      call mpgetgp('atsolu',atsolu,NHOR,1)
      call mpgetgp('assolu',assolu,NHOR,1)
      call mpgetgp('asthru',asthru,NHOR,1)
      call mpgetgp('aqvi'  ,aqvi  ,NHOR,1)
      call mpgetgp('atsa'  ,atsa  ,NHOR,1)
      call mpgetgp('ats0'  ,ats0  ,NHOR,1)

      return
      end subroutine read_atmos_restart


!     ==========================
!     SUBROUTINE STABILITY_CHECK
!     ==========================

      subroutine stability_check
      use pumamod


!     Some operating systems ignore floating point exceptions
!     and continue to run the program after an explosion,
!     e.g. after wrong settings for timestep or other conditions.
!     This subroutine checks some variables for valid ranges
!     and aborts the program if it's obviously broken.

      if (mypid == NROOT) then
         if (gd(1,1) > 1000.0 .or. gd(1,1) < -1000.0 .or. &
             gz(1,1) > 1000.0 .or. gz(1,1) < -1000.0 .or. &
             gt(1,1) > 1000.0 .or. gt(1,1) < -1000.0) then
            open(44,file='Abort_Message')
            write(44,*) 'Planet Simulator aborted'
            write(44,*) 'gd(1,1) = ',gd(1,1)
            write(44,*) 'gz(1,1) = ',gz(1,1)
            write(44,*) 'gt(1,1) = ',gt(1,1)
            close(44)
   
            write(*,*) 'Planet Simulator aborted'
            write(*,*) 'gd(1,1) = ',gd(1,1)
            write(*,*) 'gz(1,1) = ',gz(1,1)
            write(*,*) 'gt(1,1) = ',gt(1,1)
   
! FastOpt            stop
            lunstable = .true. ! FastOpt
         endif
      endif
      
      end subroutine stability_check
            

!     =================
!     SUBROUTINE INITFD
!     =================

      subroutine initfd
      use pumamod

      if (nkits < 1) nkits = 1
!     ============================================================
!     next subroutine call to set inital temperature.
!     model started from rest with stratification given by
!     trs calculated in setzt(ex). a perturbation is added to
!     constant log(sp) by subroutine noise, called from setzt(ex).
!     ============================================================
!
      if (mypid == NROOT) then
       call setzt
      endif
      call mpscsp(sp,spm,1)
      if (mypid == NROOT) then
          st(1,:) = sr(1,:)
         stm(1,:) = sr(1,:)
          sz(3,:) = EZ
         szm(3,:) = EZ
      endif
      return
      end

!     =================
!     SUBROUTINE READNL
!     =================

      subroutine readnl
      use pumamod

      namelist /inp/ kick    , nwpd    , ncoeff  , ndel    , ndiag      &
                   , nexp    , nexper  , nkits   , nrestart, noutput    &
                   , nstep   , ndebug  , ntspd             , neqsig     &
                   , nprint  , nprhor  , npacksp , npackgp , ngui       &
                   , nrad    , nflux   , ntime   , nperpetual , mars    &
                   , n_start_year , n_start_month, n_run_steps          &
                   , n_run_years , n_run_months  , n_run_days, mpstep   &
                   , n_days_per_month, n_days_per_year                  &
                   , dtep    , dtns    , dtrop   , dttrp                &
                   , tdissd  , tdissz  , tdisst  , tdissq  , tgr        &
                   , psurf   , ndl     , nhdiff                         &
                   , restim  , t0      , tfrc    , ndiaggp , ndiagsp    &
                   , ndiaggp2d, ndiaggp3d, ndiagsp2d, ndiagsp3d, ndiagcf&
                   , sigh , tnudg, dnudg,znudg,qnudg,pnudg &
                   , ncost,n_spinup,nafter,N_RESET_YEAR,makenudge,smnudg ! added, FastOpt


      namelist /planet/ akap, alr, ga, gascon, plarad, pnu, ww          &
                      , solar_day, siderial_day &
                      , ra1, ra2, ra4, yplanet
!
!     preset namelist parameter according to model set up
!

      if(NLEV==10) then
       tfrc(1)=20.
       tfrc(2)=100.
       tfrc(3:NLEV)=0.
      endif
!
      if(NTRU==42) then
       nhdiff=16
       ndel(:)=4
       tdissq(:)=0.76
       tdisst(:)=0.76
       tdissz(:)=0.3
       tdissd(:)=0.06
      endif

!
!     read namelist
!

      open(11,file='puma_namelist',form='formatted')
      read (11,inp)

!     shortcut for Mars settings

      if (mars == 1) then
         yplanet      = "Mars" 
         ntpal        = 5        ! Color pallet for T
         n_days_per_month = 57   ! simple calendar
         n_days_per_year  = 12 * n_days_per_month
         solar_day    = SOL_DAY_MARS
         siderial_day = SID_DAY_MARS
         akap    = AKAP_MARS
         alr     = ALR_MARS
         ga      = GA_MARS
         gascon  = GASCON_MARS
         plarad  = PLARAD_MARS
         pnu     = PNU_MARS
         ww      = WW_MARS
         ra1     = RA1_ICE
         ra2     = RA2_ICE
         ra4     = RA4_ICE
         dtrop     = 40000.0
         dttrp     =  20.0
         tgr       = 214.0
         psurf     = 820.0
         dveg(:)   =   0.0       ! No vegetation on Mars
         tdisst(:) =   2.0
!        restim(:) =  15.0       ! Instable
!        t0(:)     = 120.0       ! Explodes after 5 days
         tfrc(1)   =   1.0
         tfrc(2)   =   5.0
         tfrc(3:)  =  10.0
      endif

      read (11,planet)

!     set planet dependent variables

      acpd  = gascon / akap      ! Specific heat for dry air
      adv   = ACPV / acpd -1.0   ! Often used
      cv    = plarad * ww        ! cv
      ct    = cv * cv / gascon   ! ct
      pnu21 = 1.0 - 2.0 * pnu    ! Time filter 2
      rdbrv = gascon / RV        ! rd / rv

      write (*,'(/," *************************************")')
      write (*,'(" * Namelist INP from <puma_namelist> *")')
      write (*,'(" *************************************")')
      write(*,inp)
      write (*,'(/," ****************************************")')
      write (*,'(" * Namelist PLANET from <puma_namelist> *")')
      write (*,'(" ****************************************")')
      write(*,planet)
!
!     calendar and time control
!     set simulation length by using the following parameters:
!     "n_run_years" and "n_run_months"

!     Make sure that (mpstep * 60) * ntspd = solar_day

      if (mpstep > 0) then             ! timestep given in [min]
         ntspd = nint(solar_day) / (mpstep * 60)
      endif
      mpstep = solar_day  / (ntspd * 60)
      nafter = ntspd
      if (nwpd > 0 .and. nwpd <= ntspd) then
         nafter = ntspd / nwpd
      endif
   

!     Convert start date to timesteps since 1-Jan-0000

      if (n_days_per_year == 365) then
         call cal2step(n_start_step,ntspd,n_start_year,n_start_month,   &
                       1,0,0)
      else
         n_start_step = ntspd * (n_start_year     * n_days_per_year     &
                              + (n_start_month-1) * n_days_per_month)   
      endif

!     Compute simulation time in [months]

      n_run_months = n_run_years * 12 + n_run_months
      iyea = n_run_months / 12
      imon = mod(n_run_months,12)

!     Print some values

      write (*,'(/," ********************************")')
      write (*,'(" * Solar    day   :",f8.1," [s] *")') solar_day
      write (*,'(" * Siderial day   :",f8.1," [s] *")') siderial_day
      write (*,'(" * Timestep       :",i6," [min] *")') mpstep
      write (*,'(" * Timesteps / day:",i6,"       *")') ntspd

           if (iyea  > 1 .and. imon == 0) then
         write (*,'(" * Simulation time:",i5,"  years *")') iyea
      else if (iyea == 1 .and. imon == 0) then
         write (*,'(" * Simulation time:   one year  *")')
      else if (n_run_months > 1) then
         write (*,'(" * Simulation time:",i5," months *")') n_run_months
      else if (n_run_months == 1) then
         write (*,'(" * Simulation time:  one month  *")')
      else if (n_run_days  > 1) then
         write (*,'(" * Simulation time:",i5,"   days *")') n_run_days
      else if (n_run_days == 1) then
         write (*,'(" * Simulation time:    one day  *")')
      else if (n_run_steps  > 1) then
         write (*,'(" * Simulation time:",i5,"  steps *")') n_run_steps
      else if (n_run_steps == 1) then
         write (*,'(" * Simulation time: single step *")')
      endif
      write (*,'(" ********************************")')

!     set franks diagnostics


      if(ndiaggp==1) then
       ndiaggp3d=19+ndiaggp3d
      end if
      if(ndiagsp==1) then
       ndiagsp3d=3+ndiagsp3d
      end if

      return
      end

!     =================
!     SUBROUTINE INITPM
!     =================

      subroutine initpm
      use pumamod

      real (kind=8) radea,zakk

!     *************************************************************
!     * carries out all initialisation of model prior to running. *
!     * major sections identified with comments.                  *
!     * this s/r sets the model parameters and all resolution     *
!     * dependent quantities.                                     *
!     *************************************************************

      radea = plarad

!     *********************
!     * set vertical grid *
!     *********************

      if(neqsig==-1) then
       sigmah(:)=sigh(:)
      elseif(neqsig==1) then
       do jlev = 1 , NLEV
        sigmah(jlev) = real(jlev) / NLEV
       enddo
      else
       do jlev=1,NLEV
        zsk=REAL(jlev)/REAL(NLEV)
        sigmah(jlev)=0.75*zsk+1.75*zsk**3-1.5*zsk**4
       enddo
      end if

      dsigma(1     ) = sigmah(1)
      dsigma(2:NLEV) = sigmah(2:NLEV) - sigmah(1:NLEV-1)

      rdsig = 0.5 / dsigma

      sigma(1     ) = 0.5 * sigmah(1)
      sigma(2:NLEV) = 0.5 * (sigmah(1:NLEV-1) + sigmah(2:NLEV))

!     dimensionless coefficient for newtonian cooling
!     friction and timestep. of course a day is 2*PI in non dimensional
!     units using omega as the unit of frquency.

      where (restim > 0.0)
         damp = 1.0 / (TWOPI * restim)
      elsewhere
         damp = 0.0
      endwhere

      where (tfrc > 0.0)
          tfrc = 1.0 / (TWOPI * tfrc)
      elsewhere
          tfrc = 0.0
      endwhere

!     compute internal diffusion parameter (LAUERSON)

      do jlev=1,NLEV
       jdel = ndel(jlev)
       if (tdissd(jlev) > 0.0) then
        tdissd(jlev) = 1./(TWOPI*tdissd(jlev))
       else
        tdissd(jlev)=0.
       endif
       if (tdissz(jlev) > 0.0) then
        tdissz(jlev) = 1./(TWOPI*tdissz(jlev))
       else
        tdissz(jlev)=0.
       endif
       if (tdisst(jlev) > 0.0) then
        tdisst(jlev) = 1./(TWOPI*tdisst(jlev))
       else
        tdisst(jlev) = 0.
       endif
       if (tdissq(jlev) > 0.0) then
        tdissq(jlev) = 1./(TWOPI*tdissq(jlev))
       else
        tdissq(jlev)=0.
       endif
       zakk=1./((NTRU-nhdiff)**jdel)
       jr=-1
       do jm=0,NTRU
         do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            zsq = (jn - nhdiff)
            if(jn >= nhdiff) then
             sak(jr,jlev) = zakk*zsq**jdel
            else
             sak(jr,jlev) = 0.
            endif
            sak(ji,jlev) = sak(jr,jlev)
         enddo
       enddo
      enddo

if (makenudge==1) then 
     do jlev=1,NLEV
       if (tnudg(jlev)  > 0.0) then
        tnudg(jlev)  = 1./(TWOPI*tnudg(jlev))
       else
        tnudg(jlev)  = 0.
       endif
       if (dnudg(jlev)  > 0.0) then
        dnudg(jlev)  = 1./(TWOPI*dnudg(jlev) )
       else
        dnudg(jlev)  = 0.
       endif
       if (znudg(jlev)  > 0.0) then
        znudg(jlev)  = 1./(TWOPI*znudg(jlev) )
       else
        znudg(jlev)  = 0.
       endif
       if (qnudg(jlev)  > 0.0) then
        qnudg(jlev)  = 1./(TWOPI*qnudg(jlev) )
       else
        qnudg(jlev)  = 0.
       endif
     enddo

       if (pnudg> 0.0) then
        pnudg = 1./(TWOPI*pnudg)
       else
        pnudg = 0.
       endif

       if (smnudg> 0.0) then
        smnudg = 1./(TWOPI*smnudg)
       else
        smnudg = 0.
       endif
endif

!     set coefficients which depend on wavenumber

      zrsq2 = 1.0 / sqrt(2.0)

      jr=-1
      do jm=0,NTRU
         do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            nindex(jr)=jn
            nindex(ji)=jn
            spnorm(jr)=zrsq2
            spnorm(ji)=zrsq2
         enddo
         zrsq2=-zrsq2
      enddo

! finally make temperatures dimensionless

      dtns  = dtns  / ct
      dtep  = dtep  / ct
      dttrp = dttrp / ct
      t0    = t0    / ct ! (NLEV)

!     print out

      zakk=tdisst(NLEV)/((NTRU-nhdiff)**ndel(NLEV))
      write(*,'(/," *************************************************")')
      if (zakk == 0.0) then
      write (*,'(" * No lateral dissipation *")')
      else
      write(*,'(" * Lateral dissipation",13x,"NDEL(",i3,") =",i2," *")')&
         NLEV,ndel(NLEV)
      write(*,'(" * Diffusion coefficient = ",g14.4," [m**",i1,"] *")')&
         zakk*ww*radea**ndel(NLEV),ndel(NLEV)
      write(*,'(" * e-folding time for smallest scale =",f5.1," days *")')&
         1.0/(TWOPI*tdisst(NLEV))
      endif
      write(*,'(" *************************************************")')
      return
      end


!     =================
!     SUBROUTINE MAKEBM
!     =================

      subroutine makebm
      use pumamod

      zdeltsq = delt * delt

      do jlev1 = 1 , NLEV
         do jlev2 = 1 , NLEV
            zaq = zdeltsq * (t0(jlev1) * dsigma(jlev2)                  &
     &          + dot_product(g(:,jlev1),tau(jlev2,:)))
            bm1(jlev2,jlev1,:) = zaq
         enddo
      enddo

      do jn=1,NTRU
         do jlev = 1 , NLEV
            bm1(jlev,jlev,jn) = bm1(jlev,jlev,jn) + 1.0 / (jn*(jn+1))
         enddo
         call minvers(bm1(1,1,jn),NLEV)
      enddo
      return
      end

!     =================
!     SUBROUTINE INITSI
!     =================

      subroutine initsi
      use pumamod
!===========================================================
! carries out all initialisation of model prior to running.
! major sections identified with comments.
! this s/r sets the variables and arrays associated with
! the semi-implicit scheme.
!===========================================================
!
      dimension zalp(NLEV),zh(NLEV)
!
! this value, used in setting alpha(1), is irrelevant in the
! angular momentum conserving ecmwf scheme

      tkp = akap * t0
      t01s2(1:NLEM) = t0(2:NLEV) - t0(1:NLEM)
      t01s2(  NLEV) = 0.0

      zalp(2:NLEV) = log(sigmah(2:NLEV)) - log(sigmah(1:NLEM))

      g      = 0.0
      g(1,1) = 1.0
      do jlev = 2 , NLEV
         g(jlev,jlev) = 1.0 - zalp(jlev)*sigmah(jlev-1)/dsigma(jlev)
         g(jlev,1:jlev-1) = zalp(jlev)
      enddo

      do jlev = 1 , NLEV
         c(jlev,:) = g(:,jlev) * (dsigma(jlev) / dsigma(:))
      enddo

      zt01s2   = t01s2(1)
      zsig     = sigmah(1)
      tau(1,1) = 0.5 * zt01s2 * (zsig - 1.0) + tkp(1)
      tau(2:NLEV,1) = 0.5 * zt01s2 * dsigma(2:NLEV)

      do 1410 jlev=2,NLEV
        zttm=zt01s2
        zsigm=zsig
        zt01s2=t01s2(jlev)
        zsig=sigmah(jlev)
        do 1420 jlev2=1,NLEV
          ztm=0.
          ztmm=0.
          if(jlev2.le.jlev) ztm=1
          if(jlev2.lt.jlev) ztmm=1
          ztau=zttm*(zsigm-ztmm)
          if(jlev.lt.NLEV) ztau=ztau+zt01s2*(zsig-ztm)
          ztau=ztau*rdsig(jlev)*dsigma(jlev2)
          if(jlev2.le.jlev) ztau=ztau+tkp(jlev)*c(jlev2,jlev)
          tau(jlev2,jlev)=ztau
 1420   continue
 1410 continue
!
      zfctr=0.001*CV*CV/ga
      do 1500 jlev = 1 , NLEV
         zh(jlev) = dot_product(g(:,jlev),t0) * zfctr
 1500 continue
!
!     **********************************
!     * write out vertical information *
!     **********************************
!
      write(*,9001)
      write(*,9003)
      write(*,9002)
      do jlev = 1 , NLEV
        write(*,9004) jlev,sigma(jlev),t0(jlev),zh(jlev)
      enddo
      write(*,9000)

      write(*,9012)
      write(*,9013) (jlev,jlev = 1 , 5)
      write(*,9012)
      do jlev = 1 , NLEV
        write(*,9014) jlev,(c(i,jlev),i=1,5)
      enddo
      write(*,9012)
      return
 9000 format(1x,33('*'),/)
 9001 format(/,1x,33('*'))
 9002 format(1x,33('*'))
 9003 format(' * Lv *    Sigma Basic-T  Height *')
 9004 format(' *',i3,' * ',3f8.3,' *')
 9012 format(1x,69('*'))
 9013 format(' * Lv * C',i11,4i12,' *')
 9014 format(' *',i3,' * ',5f12.8,' *')
      end

!     ==================
!     SUBROUTINE MINVERS
!     ==================

      subroutine minvers(a,n)
      dimension a(n,n),b(n,n),indx(n)

      b = 0.0
      do j = 1 , n
         b(j,j) = 1.0
      enddo
      call ludcmp(a,n,indx)
      do j = 1 , n
         call lubksb(a,n,indx,b(1,j))
      enddo
      a = b
      return
      end

!     =================
!     SUBROUTINE LUBKSB
!     =================

      subroutine lubksb(a,n,indx,b)
      dimension a(n,n),b(n),indx(n)
      k = 0
      do i = 1 , n
         l    = indx(i)
         sum  = b(l)
         b(l) = b(i)
         if (k > 0) then
            do j = k , i-1
               sum = sum - a(i,j) * b(j)
            enddo
         else if (sum /= 0.0) then
            k = i
         endif
         b(i) = sum
      enddo

      do i = n , 1 , -1
         sum = b(i)
         do j = i+1 , n
            sum = sum - a(i,j) * b(j)
         enddo
         b(i) = sum / a(i,i)
      enddo
      return
      end

!     =================
!     SUBROUTINE LUDCMP
!     =================

      subroutine ludcmp(a,n,indx)
      dimension a(n,n),indx(n),vv(n)

      d = 1.0
      vv = 1.0 / maxval(abs(a),2)

      do 19 j = 1 , n
         do i = 2 , j-1
            a(i,j) = a(i,j) - dot_product(a(i,1:i-1),a(1:i-1,j))
         enddo
         aamax = 0.0
         do i = j , n
            if (j > 1)                                                  &
     &      a(i,j) = a(i,j) - dot_product(a(i,1:j-1),a(1:j-1,j))
            dum = vv(i) * abs(a(i,j))
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            endif
         enddo
         if (j .ne. imax) then
            do 17 k = 1 , n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
   17       continue
            d = -d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
         if (a(j,j) == 0.0) a(j,j) = tiny(a(j,j))
         if (j < n) a(j+1:n,j) = a(j+1:n,j) / a(j,j)
   19 continue
      return
      end

!     =================
!     SUBROUTINE INILAT
!     =================

      subroutine inilat
      use pumamod
      character(len=1) ch
      
      gaussw(:) = gw(:) ! added, FastOpt
      ch = 'N'
      do jlat = 1 , NLAT
         csq(jlat)  = 1.0 - si(jlat) * si(jlat)
         gw(jlat)   = gw(jlat) / csq(jlat)
         rcs(jlat)  = 1.0 / sqrt(csq(jlat))
      enddo
      do jlat = 1 , NLAT/2
         ideg = nint(180.0/PI * asin(si(jlat)))
         write(chlat(jlat),'(i2,a1)') ideg,'N'
         write(chlat(NLAT+1-jlat),'(i2,a1)') ideg,'S'
      enddo
      return
      end

!     ================
!     SUBROUTINE NOISE
!     ================

      subroutine noise
      use pumamod

!     if kick is set to 1 or 2
!     adds white noise perturbation to ln(surface pressure)
!     balanced initial state at t=0.
!     for kick=2, the white noise pertubation is
!     symmetric to the equator
!     eps sets magnitude of the noise

      call random_seed()
      zeps=1.e-4
      zscale=zeps/sqrt(2.0)

      write(*,'(/," ******************************************")')
      if (kick == 1) then
         jsp1=2*NTP1+1
         do jsp=jsp1,NRSP
            call random_number(zrand)
            sp(jsp)=sp(jsp)+zscale*(zrand-0.5)
         enddo
         write(*,'(" *      White noise added (KICK = 1)      *")')
      elseif (kick == 2) then
         jr=2*NTP1-1
         do jm=1,NTRU
            do jn=jm,NTRU
               jr=jr+2
               ji=jr+1
               if (mod(jn+jm,2) == 0) then
                  call random_number(zrand)
                  sp(jr)=sp(jr)+zscale*(zrand-0.5)
                  sp(ji)=sp(ji)+zscale*(zrand-0.5)
               endif
            enddo
         enddo
         write(*,'(" * Symmetric white noise added (KICK = 2) *")')
         print *,'symmetric white noise added'
      elseif (kick == 3) then
         sp(2*NTP1+3) = zscale
         sp(2*NTP1+4) = zscale * 0.5
         write(*,'(" *   Mode sp(1,1) disturbed (KICK = 3)    *")')
      endif
      write(*,'(" ******************************************")')
      return
      end

!     ================
!     SUBROUTINE SETZT
!     ================
      subroutine setzt
      use pumamod
!
      dimension ztrs(NLEV)
      dimension zfac(NLEV)
!
!*********************************************************************
!  this s/r sets up restoration temp field.
! the temperature at sigma = 1 is tgr, entered in kelvin.
! a lapse rate of ALR k/m is assumed under the tropopause and zero
! above. the actual profile tends to this away from then
! tropopause, with smooth interpolation depending on dttrp
! at the model tropopause.the height of
! the tropopause is given as dtrop m.
!*********************************************************************
!
      sr(:,:) = 0.0 ! NESP,NLEV
!
      zdttrp=dttrp*ct
!
      zsigprev=1.
      ztprev=tgr
      zzprev=0.
      do 1100 jlev=NLEV,1,-1
        zzp=zzprev+(gascon*ztprev/ga)*log(zsigprev/sigma(jlev))
        ztp=tgr-dtrop*ALR
        ztp=ztp+sqrt((.5*ALR*(zzp-dtrop))**2+zdttrp**2)
        ztp=ztp-.5*ALR*(zzp-dtrop)
        ztpm=.5*(ztprev+ztp)
        zzpp=zzprev+(gascon*ztpm/ga)*log(zsigprev/sigma(jlev))
        ztpp=tgr-dtrop*ALR
        ztpp=ztpp+sqrt((.5*ALR*(zzpp-dtrop))**2+zdttrp**2)
        ztpp=ztpp-.5*ALR*(zzpp-dtrop)
        ztrs(jlev)=ztpp
        zzprev=zzprev                                                   &
     &        +(.5*(ztpp+ztprev)*gascon/ga)*log(zsigprev/sigma(jlev))
        ztprev=ztpp
        zsigprev=sigma(jlev)
1100  continue
!
!     **********************************
!     * write out vertical information *
!     **********************************
!
      write(*,9001)
      write(*,9003)
      write(*,9002)
!
      do 1200 jlev = 1 , NLEV
         write(*,9004) jlev,sigma(jlev),ztrs(jlev)
         ztrs(jlev)=ztrs(jlev)/ct
 1200 continue
!
      write(*,9002)
!
!******************************************************************
! loop to set array zfac - this controls temperature gradients as a
! function of sigma in tres. it is a sine wave from one at
! sigma = 1 to zero at stps (sigma at the tropopause) .
!******************************************************************
! first find sigma at dtrop
!
      zttrop=tgr-dtrop*ALR
      ztps=(zttrop/tgr)**(ga/(ALR*gascon))
!
! now the latitudinal variation in tres is set up ( this being in terms
! of a deviation from t0 which is usually constant with height)
!
      zsqrt2=sqrt(2.)
      zsqrt04=sqrt(0.4)
      zsqrt6=sqrt(6.)
      do 2100 jlev = 1 , NLEV
        zfac(jlev)=sin(0.5*PI*(sigma(jlev)-ztps)/(1.-ztps))
        if (zfac(jlev).lt.0.0) zfac(jlev)=0.0
        sr(1,jlev)=zsqrt2*(ztrs(jlev)-t0(jlev))
        sr(3,jlev)=(1./zsqrt6)*dtns*zfac(jlev)
        sr(5,jlev)=-2./3.*zsqrt04*dtep*zfac(jlev)
 2100 continue
!
      call noise
!
      return
 9001 format(/,1x,26('*'))
 9002 format(1x,26('*'))
 9003 format(' * Lv *    Sigma Restor-T *')
 9004 format(' *',i3,' * ',f8.3,f9.3,' *')
      end


!     ===============
!     SUBROUTINE DIAG
!     ===============

      subroutine diag
      use pumamod
      if (mod(nstep,ndiag) == 0) then
         if (ncoeff .gt. 0) call prisp
         call xsect
      endif
      call energy
      return
      end

!     ================
!     SUBROUTINE PRISP
!     ================

      subroutine prisp
      use pumamod

      character(len=30) title

      scale = 100.0
      title = 'Vorticity [10-2]'
      do 100 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(sz(1,jlev),jlev,title,scale)
  100 continue

      title = 'Divergence [10-2]'
      do 200 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(sd(1,jlev),jlev,title,scale)
  200 continue

      scale = 1000.0
      title = 'Temperature [10-3]'
      do 300 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(st(1,jlev),jlev,title,scale)
  300 continue

      scale = 1000.0
      title = 'Specific Humidity [10-3]'
      do 400 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(sq(1,jlev),jlev,title,scale)
  400 continue

      title = 'Pressure [10-3]'
      call wrspam(sp,0,title,scale)

      return
      end

!     ==============
!     FUNCTION RMSSP
!     ==============

      function rmssp(pf)
      use pumamod
      real pf(NESP,NLEV)

      zsum = 0.0
      do jlev = 1 , NLEV
         zsum = zsum + dsigma(jlev)                                     &
     &        * (dot_product(pf(1:NZOM,jlev),pf(1:NZOM,jlev)) * 0.5     &
     &        +  dot_product(pf(NZOM+1:NRSP,jlev),pf(NZOM+1:NRSP,jlev)))
      enddo
      rmssp = zsum
      return
      end

!     =================
!     SUBROUTINE ENERGY
!     =================

      subroutine energy
      use pumamod

      parameter (idim=5) ! Number of scalars for GUI timeseries
      real (kind=4) ziso(idim)

      ziso(1) = umax   ! maximum value of csu (zonal mean cross section)
      ziso(2) = t2mean - TMELT  ! mean of 2m temperature
      ziso(3) = precip * 1.0e9  ! mean precipitation
      ziso(4) = evap   * 1.0e9  ! mean evaporation
      ziso(5) = olr             ! OLR

      call guiput("SCALAR" // char(0) ,ziso,idim,1,1)

      return
      end

!     ==================
!     SUBROUTINE UPDATIM
!     ==================

      subroutine updatim(kstep)
      use pumamod
       
      interface
         integer function nweekday(kday)
            integer, intent(in) :: kday
         end function nweekday
      end interface

      if (n_days_per_year == 365) then
         call step2cal(kstep,ntspd, &
            ndatim(1),ndatim(2),ndatim(3),ndatim(4),ndatim(5))
         ndatim(6) = nweekday(kstep/ntspd)
      else
         call step2cal30(kstep,ntspd, &
            ndatim(1),ndatim(2),ndatim(3),ndatim(4),ndatim(5))
      endif
      return
      end

!     =================
!     SUBROUTINE WRSPAM
!     =================

      subroutine wrspam(ps,klev,title,scale)
      use pumamod
!
      dimension ps(NRSP)
      character(len=30) title
      character(len=18) datch


      call ntodat(nstep,datch)
      write(*,'(1x)')
      write(*,20000)
      write(*,20030) datch,title,klev
      write(*,20000)
      write(*,20020) (i,i=0,9)
      write(*,20000)
      write(*,20100) (cab(i),i=1,10)
      write(*,20200) (cab(i),i=NTRU+2,NTRU+10)
      write(*,20300) (cab(i),i=2*NTRU+2,2*NTRU+9)
      write(*,20400) (cab(i),i=3*NTRU+1,3*NTRU+7)
      write(*,20000)
      write(*,'(1x)')

      return

20000 format(1x,78('*'))
20020 format(' * n * ',10i7,' *')
20030 format(' *   * ',a18,2x,a30,'  Level ',i2,11x,'*')
20100 format(' * 0 *',f8.2,9f7.2,' *')
20200 format(' * 1 *',8x,9f7.2,' *')
20300 format(' * 2 *',15x,8f7.2,' *')
20400 format(' * 3 *',22x,7f7.2,' *')
      contains
      function cab(i)
      cab=real(scale*sqrt(ps(i+i-1)*ps(i+i-1)+ps(i+i)*ps(i+i)))
      end function cab
      end

!     ===============
!     SUBROUTINE WRZS
!     ===============

      subroutine wrzs(zs,title,scale)
      use pumamod
!
      dimension zs(NLAT,NLEV)
      character(len=30) title
      character(len=18) datch

      ip = NLAT / 16
      ia = ip/2
      ib = ia + 7 * ip
      id = NLAT + 1 - ia
      ic = id - 7 * ip

      call ntodat(nstep,datch)
      write(*,'(1x)')
      write(*,20000)
      write(*,20030) datch,title
      write(*,20000)
      write(*,20020) (chlat(i),i=ia,ib,ip),(chlat(j),j=ic,id,ip)
      write(*,20000)
      do 200 jlev = 1 , NLEV
         write(*,20100) jlev,((int(zs(i,jlev)*scale)),i=ia,ib,ip),      &
     &                       ((int(zs(j,jlev)*scale)),j=ic,id,ip),jlev
  200 continue
      write(*,20000)
      write(*,'(1x)')

20000 format(1x,78('*'))
20020 format(' * Lv * ',16(1x,a3),' * Lv *')
20030 format(' *    * ',a18,2x,a30,20x,'*')
20100 format(' * ',i2,' * ',16i4,' * ',i2,' *')
      end

!     ================
!     SUBROUTINE XSECT
!     ================

      subroutine xsect
      use pumamod
      character(len=30) title

      scale = 10.0
      title = 'Zonal Wind [0.1 m/s]'
      call wrzs(csu,title,scale)
      title = 'Meridional Wind [0.1 m/s]'
      call wrzs(csv,title,scale)
      scale = 1.0
      title = 'Temperature [C]'
      call wrzs(cst,title,scale)
      scale = 10000.0
      title = 'specific humidity [0.1g/Kg]'
      call wrzs(csm,title,scale)
      scale = 100.0
      title = 'cloud cover [%]'
      call wrzs(ccc,title,scale)
      return
      end


!     =====================
!     SUBROUTINE GRIDPOINTA
!     =====================

      subroutine gridpointa
      use pumamod
      implicit none ! FastOpt
!
!*    Adiabatic Gridpoint Calculations
!
      real gtn(NHOR,NLEV)
      real gut(NHOR,NLEV)
      real gvt(NHOR,NLEV)
      real gqn(NHOR,NLEV)
      real guq(NHOR,NLEV)
      real gvq(NHOR,NLEV)
      real guz(NHOR,NLEV)
      real gvz(NHOR,NLEV)
      real gke(NHOR,NLEV)
      real gphi(NHOR,NLEV)
      real gvpp(NHOR)
      real gpmt(NLON,NLPP)
      real sdf(NESP,NLEV)
      real stf(NESP,NLEV)
      real sqf(NESP,NLEV)
      real szf(NESP,NLEV)
      real spf(NESP)

      real zgp(NLON,NLAT)
      real zqout(NHOR,NLEV)
!
!     local variables
!
      integer :: jlev,jlat,sec,j1,j2,jlon ! FastOpt
!
!     inverse Legendre transformation (spectral to fourier domain)
!     st -> gt  sq -> gq  sd -> gd  sz -> gz
!     (sd,sz) -> (gu,gv)
!     sp -> gp  sp -> gpj (dlnps/dphi)
!

      call invlega

      if (ngui > 0 .or. mod(nstep,ndiag) == 0) then
        do jlev = 1 , NLEV
          do jlat = 1 , NLPP
            sec = CV / sqrt(csq(jlat))
            csu(jlat,jlev) =  gu(1+(jlat-1)*NLON,jlev) * sec
            csv(jlat,jlev) =  gv(1+(jlat-1)*NLON,jlev) * sec
            cst(jlat,jlev) =(gt(1+(jlat-1)*NLON,jlev) + t0(jlev))*ct-TMELT
            csm(jlat,jlev) = (gq(1+(jlat-1)*NLON,jlev))
            j1=(jlat-1)*NLON+1
            j2=jlat*NLON
            ccc(jlat,jlev) = SUM(dcc(j1:j2,jlev))/real(NLON)
          enddo
        enddo
        umax = maxval(csu)
      endif

      do jlat = 1 , NLPP
         do jlon = 1 , NLON , 2
           gpmt(jlon  ,jlat) = -gp(jlon+1+(jlat-1)*NLON) * ((jlon-1)/2)
           gpmt(jlon+1,jlat) =  gp(jlon  +(jlat-1)*NLON) * ((jlon-1)/2)
         enddo
      enddo

      call fc2gp(gu  ,NLON,NLPP*NLEV)
      call fc2gp(gv  ,NLON,NLPP*NLEV)
      call fc2gp(gt  ,NLON,NLPP*NLEV)
      call fc2gp(gq  ,NLON,NLPP*NLEV)
      call fc2gp(gd  ,NLON,NLPP*NLEV)
      call fc2gp(gz  ,NLON,NLPP*NLEV)
      call fc2gp(gpj ,NLON,NLPP)
      call fc2gp(gpmt,NLON,NLPP)
      call fc2gp(gp  ,NLON,NLPP)
!$TAF STORE gp = inner1, REC = 1 + nstep2
      gp = exp(gp)
!$TAF STORE gu,gv,gt,gq,gd,gz,gpj,gpmt,gp = inner1, REC = 1 + nstep2

      call calcgp(gtn,gqn,guz,gvz,gpmt,gvpp,gphi)

      gut = gu * gt
      gvt = gv * gt
      guq = gu * gq
      gvq = gv * gq
      gke = gu * gu + gv * gv
!
!     add non linear geopotential terms
!
      gke(:,:) = gke(:,:) + gphi(:,:)

!
!     fft
!

      call gp2fc(gtn ,NLON,NLPP*NLEV)
      call gp2fc(gqn ,NLON,NLPP*NLEV)

      call gp2fc(gut ,NLON,NLPP*NLEV)
      call gp2fc(gvt ,NLON,NLPP*NLEV)
      call gp2fc(guq ,NLON,NLPP*NLEV)
      call gp2fc(gvq ,NLON,NLPP*NLEV)
      call gp2fc(guz ,NLON,NLPP*NLEV)
      call gp2fc(gvz ,NLON,NLPP*NLEV)
      call gp2fc(gke ,NLON,NLPP*NLEV)
      call gp2fc(gvpp,NLON,NLPP     )

!
!     direct Legendre transformation (fourier domain to spectral domain)
!

      call dirlega(sdf,stf,sqf,szf,spf,gke,gtn,gqn,gut,guz,gvt,gvz      &
     &            ,guq,gvq,gvpp)                                     !NEU

      call mpsumsc(spf,spt,1)
      call mpsumsc(stf,stt,NLEV)
      call mpsumsc(sqf,sqt,NLEV)
      call mpsumsc(sdf,sdt,NLEV)
      call mpsumsc(szf,szt,NLEV)

!     Diagnostic output and GUI calls

      if (ngui>0 .or. mod(nstep,ndiag)==0 .or. mod(nstep,nafter)==0) then
         do jlev=1,NLEV
          zqout(:,jlev)=gq(:,jlev)/gp(:)
          call gp2fc(zqout(:,jlev),NLON,NLPP)
          call fc2sp(zqout(:,jlev),sqout(:,jlev))
         enddo
         call mpsum(sqout,NLEV)
         call mpgagp(zgp,gp,1);
         call guips(zgp)                ! Send Ps for Hovmoeller
         call guihor("DQVI"//char(0),dqvi,1,1.0,0.0)! Vertically integrated q
         call guigv("GU"  // char(0),gu)            ! Send u to GUI
         call guigv("GV"  // char(0),gv)            ! Send v to GUI
         call mpgagp(zgp,dtsa,1)
         t2mean = sum(zgp) / (NLON * NLAT) ! Mean of 2m temperature
         call mpgagp(zgp,dprc,1)        ! Convective precip
         precip = sum(zgp)
         call mpgagp(zgp,dprl,1)        ! Large scale precip
         precip = (precip + sum(zgp)) / (NLON * NLAT)
         call mpgagp(zgp,devap,1)       ! Evaporation
         evap = sum(zgp) / (NLON * NLAT)
         call mpgagp(zgp,dftu,1)        ! OLR = dftu level 1
         olr = -sum(zgp) / (NLON * NLAT)! make positive upwards
         gp(:) = gp(:) - 1.0
         call gp2fc(gp,NLON,NLPP)
         call fc2sp(gp,span)
         call mpsum(span,1)
 
         call mpgacs(csu)
         call mpgacs(csv)
         call mpgacs(cst)
         call mpgacs(csm)
         call mpgacs(ccc)
 
         if (mypid == NROOT) then
            call guiput("CSU"  // char(0) ,real(csu ,kind=4),NLAT, NLEV,1)
            call guiput("CSV"  // char(0) ,real(csv ,kind=4),NLAT, NLEV,1)
            call guiput("CST"  // char(0) ,real(cst ,kind=4),NLAT, NLEV,1)
            call guiput("SPAN" // char(0) ,real(span,kind=4),NCSP,-NTP1,1)
         endif
      endif
      return
      end

!     =================
!     SUBROUTINE CALCGP
!     =================

      subroutine calcgp(gtn,gqn,guz,gvz,gpm,gvp,gphi)

!     *****************************************************
!     * computes nonlinear tendencies in grid point space *
!     *****************************************************

      use pumamod

      real gtn(NHOR,NLEV)
      real gqn(NHOR,NLEV)                                                       !NEU
      real gphi(NHOR,NLEV)
      real guz(NHOR,NLEV), gvz(NHOR,NLEV)
      real gpm(NHOR)     , gvp(NHOR)
      real zsdotp(NHOR,NLEM),zsumd(NHOR)
      real ztpta(NHOR),ztptb(NHOR)
      real zvgpg(NHOR,NLEV)
      real gtd(NHOR,NLEM)
      real gqm(NHOR,NLEM)                                                       !NEU

      real gud(NHOR,NLEM)
      real gvd(NHOR,NLEM)

      real ztv1(NHOR,NLEV),ztv2(NHOR,NLEV),zq(NHOR)

      do jlev = 1 , NLEV
         zvgpg(:,jlev) = rcsq  * (gu(:,jlev) * gpm + gv(:,jlev) * gpj)
!
!     set pseudo temperatures to use virtual temperatures
!     (note: gq=ps*q)
!
         zq(:)=AMAX1(gq(:,jlev)/gp(:),0.)
         ztv1(:,jlev)=(gt(:,jlev)+t0(jlev))*(1.+(1./rdbrv-1.)*zq(:))    &
     &               -t0(jlev)
         ztv2(:,jlev)=(gt(:,jlev)+t0(jlev))*(1.+(1./rdbrv-1.)*zq(:))    &
     &               /(1.+ADV*zq(:))                                    &
     &               -t0(jlev)

      enddo

!     *******
!     * gvp *
!     *******

      zsumd = dsigma(1) * gd(:,1)
      gvp   = dsigma(1) * zvgpg(:,1)
      zsdotp(:,1) = zsumd + gvp

      do jlev = 2 , NLEV-1
         zsumd = zsumd + dsigma(jlev) * gd(:,jlev)
         gvp   = gvp   + dsigma(jlev) * zvgpg(:,jlev)
         zsdotp(:,jlev) = zsumd + gvp
      enddo

      zsumd = zsumd + dsigma(NLEV) * gd(:,NLEV)
      gvp   = gvp   + dsigma(NLEV) * zvgpg(:,NLEV)

!     **************
!     * loop  400: *
!     **************
!$taf incompl gtd,gqm,gud,gvd
      do jlev = 1 , NLEM
         zsdotp(:,jlev) = (sigmah(jlev) * (zsumd+gvp) - zsdotp(:,jlev))
         gtd(:,jlev) = zsdotp(:,jlev) * (gt(:,jlev+1) - gt(:,jlev))
         gqm(:,jlev) = zsdotp(:,jlev) * (gq(:,jlev+1) + gq(:,jlev))
         gud(:,jlev) = zsdotp(:,jlev) * (gu(:,jlev+1) - gu(:,jlev))
         gvd(:,jlev) = zsdotp(:,jlev) * (gv(:,jlev+1) - gv(:,jlev))
      enddo

!     *************
!     * top level *
!     *************

      zsumd = zvgpg(:,1) * dsigma(1)

      gtn(:,1) = gt(:,1) * gd(:,1) - akap * ztv2(:,1) * gd(:,1)         &
     &   - rdsig(1)*(gtd(:,1) + t01s2(1) * (sigmah(1)*gvp-zsumd))

      gqn(:,1) = - rdsig(1) * gqm(:,1)

      guz(:,1) =-gu(:,1) * gz(:,1) - gpj * ztv1(:,1) - rdsig(1)*gvd(:,1)
      gvz(:,1) = gv(:,1) * gz(:,1) - gpm * ztv1(:,1) - rdsig(1)*gud(:,1)

      dw(:,1)=gd(:,1)*gp(:)*psurf*ww

!     ****************
!     * inner levels *
!     ****************
!$taf incompl gtn,gqn,guz,gvz,dw
      do jlev = 2 , NLEV-1
         ztpta = c(1,jlev) *  zvgpg(:,1)
         ztptb = c(1,jlev) * (zvgpg(:,1) + gd(:,1))

         do jlej = 2 , jlev
            ztpta = ztpta + c(jlej,jlev) *  zvgpg(:,jlej)
            ztptb = ztptb + c(jlej,jlev) * (zvgpg(:,jlej) + gd(:,jlej))
         enddo
!$TAF STORE ztptb = calcgpn, REC = (jlev-1) + nstep2 * (NLEV-2)
         zsumd = zsumd + zvgpg(:,jlev) * dsigma(jlev)

         gtn(:,jlev) = gt(:,jlev) * gd(:,jlev)                          &
     &       + akap * ztv2(:,jlev) * (zvgpg(:,jlev) - ztptb)            &
     &       + tkp(jlev) * (zvgpg(:,jlev) - ztpta)                      &
     &       - rdsig(jlev) * (gtd(:,jlev) + gtd(:,jlev-1)               &
     &                       +gvp*(t01s2(jlev)*sigmah(jlev)             &
     &                            +t01s2(jlev-1)*sigmah(jlev-1))        &
     &                       -zsumd*(t01s2(jlev-1)+t01s2(jlev))         &
     &                       +zvgpg(:,jlev)*dsigma(jlev)*t01s2(jlev-1))

         gqn(:,jlev) = - rdsig(jlev)*(gqm(:,jlev) - gqm(:,jlev-1))

         guz(:,jlev) = - gu(:,jlev) * gz(:,jlev) - gpj * ztv1(:,jlev)   &
     &       - rdsig(jlev)*(gvd(:,jlev) + gvd(:,jlev-1))

         gvz(:,jlev) =   gv(:,jlev) * gz(:,jlev) - gpm * ztv1(:,jlev)   &
     &       - rdsig(jlev)*(gud(:,jlev) + gud(:,jlev-1))

         dw(:,jlev)=(zvgpg(:,jlev)-ztptb)*gp(:)*psurf*ww

      enddo

!     ****************
!     * bottom level *
!     ****************

      ztpta = c(1,NLEV) *  zvgpg(:,1)
      ztptb = c(1,NLEV) * (zvgpg(:,1) + gd(:,1))

      do jlej = 2 , NLEV
         ztpta = ztpta + c(jlej,NLEV) *  zvgpg(:,jlej)
         ztptb = ztptb + c(jlej,NLEV) * (zvgpg(:,jlej) + gd(:,jlej))
      enddo
!$TAF STORE ztptb = inner1, REC = 1 + nstep2
      gtn(:,NLEV) = gt(:,NLEV) * gd(:,NLEV)                             &
     &   + akap*ztv2(:,NLEV)*(zvgpg(:,NLEV)-ztptb)                      &
     &   + tkp(NLEV)*(zvgpg(:,NLEV)-ztpta)                              &
     &   - rdsig(NLEV)*(gtd(:,NLEM)                                     &
     &                 +t01s2(NLEV-1)*(sigmah(NLEV-1)*gvp-zsumd))

      gqn(:,NLEV) = rdsig(NLEV) * gqm(:,NLEM)

      guz(:,NLEV) = -gu(:,NLEV)*gz(:,NLEV) - gpj*ztv1(:,NLEV)           &
     &    - rdsig(NLEV) * gvd(:,NLEM)
      gvz(:,NLEV) =  gv(:,NLEV)*gz(:,NLEV) - gpm*ztv1(:,NLEV)           &
     &    - rdsig(NLEV) * gud(:,NLEM)

      dw(:,NLEV)=(zvgpg(:,NLEV)-ztptb)*gp(:)*psurf*ww

!
!     compute non linear geopotential terms
!

      ztv1(:,:)=ztv1(:,:)-gt(:,:)
      do jlev=1,NLEV
      do jhor=1,NHOR
       gphi(jhor,jlev)=dot_product(g(:,jlev),ztv1(jhor,:))              &
     &                *2./rcsq(jhor)
      enddo
      enddo

      return
      end

!     ====================
!     SUBROUTINE SPECTRALA
!     ====================

      subroutine spectrala
      use pumamod
!
!*    Add adiabatic and diabatic tendencies
!
!     the adiabatic tendencies are added using the semi implicit scheme
!     described in
!     Hoskins and Simmons 1975 (Q.J.R.Meteorol.Soc.,101,637-655) (HS75)
!     To compare the code directly with HS75 the following notes might be
!     helpful (in addition to the comments below):
!
!     - *x* means model variable x
!     - eq.x referres to equation x of HS75
!     - the script T of HS75 (e.g. 1st term rhs of eq.9) is stored in *stt*
!     - the script D of HS75 (e.g. 1st term rhs of eq.8) is stored in *sdt*
!     - the script P of HS75 (e.g. 1st term rhs of eq.10) is stored in -*spt*
!     - the surface geopotential is stored in *so*
!
!     - the vertical scheme has being changed to the ECMWF scheme
!       (see e.g. Simmons and Burridge 1981, Mon.Wea.Rev.,109,758-766).
!       in this scheme,  matrix g differs from that in HS75.
!
!     - in addition to the dry HS75 version also the tendencys of specific
!       humidity are processed
!
      real azm(NSPP,NLEV)
      real atm(NSPP,NLEV)
      real aqm(NSPP,NLEV)
      real adt(NSPP,NLEV)
      real adm(NSPP,NLEV)
      real zgt(NSPP,NLEV)
      real zgm(NSPP,NLEV)
      real apm(NSPP)
!
!     franks diagnostics
!
      real, allocatable :: ztt(:,:)
!
!*    0. save prognostic variables at (t-dt)
!        and the non-linear divergence tendency terms
!
      apm(:)   = spm(:)   ! log surface pressure
      adm(:,:) = sdm(:,:) ! divergence
      azm(:,:) = szm(:,:) ! (absolut) vorticity
      atm(:,:) = stm(:,:) ! temperature
      aqm(:,:) = sqm(:,:) ! spec.  humidity
      adt(:,:) = sdt(:,:) ! divergence tendency
!
!*    1. calculate divergence on timelevel t (solving eq.17),
!        which will replace the divergence tendency sdt
!
!     1.a precompute (g*script T; *zgt*) and (phi-phi*=g*T, eq.11; *zgm*)
!         needed for the rhs of eq.17.
!         (note that phi is needed in eq.17 and, therefor,
!         the surface geopotential phi* is added later (loop 1.b))
!
      do jlev = 1 , NLEV
         do jsp=1,NSPP
            zgt(jsp,jlev) = dot_product(g(:,jlev),stt(jsp,:))
            zgm(jsp,jlev) = dot_product(g(:,jlev),atm(jsp,:))
         enddo
      enddo
!
!     1.b compute divergence at time t (i.e. solve eq.17)
!         and overwrite *sdt* with the result
!         for comparison with HS75 note:
!
!         - *spt* contains -(!)script P
!         - *spm* contains log(ps)(t-dt)
!         - surface geopotential (*so*) needs to be added (see 1.a)
!         - bm1 is the invers of matrix (1/cn I+B dt**2) (lhs eq.17)
!         - zz is set to the rhs of eq.17
!         - zsum holds the result of the dot product of rhs(eq17) and bm1
!           (therefor jlev2-loop)
!

      do jlev = 1 , NLEV
        do jsp = 1 , NSPP
          jn = nindex(jsp)
          zsum = 0.0
          if (jn > 0) then
            zq = 1.0 / (jn * (jn+1))
            do jlev2 = 1 , NLEV
              z0 = t0(jlev2)
              zt = zgt(jsp,jlev2) - z0 * spt(jsp)
              zm = zgm(jsp,jlev2) + z0 * spm(jsp)
              za = adt(jsp,jlev2) * zq + sop(jsp)
              zb = adm(jsp,jlev2) * zq
              zz = zb + delt * (zm + za + delt * zt)
              zsum = zsum + zz * bm1(jlev2,jlev,jn)
            enddo
          endif
          sdt(jsp,jlev) = zsum
        enddo
      enddo
      if (mypid == NROOT) then
       sdt(1:2,:) = 0.0 ! first mode should be zero (errors due to numerics)
      endif
!
!*    2. calculate (-log) surface pressure tendency from eq.15 (pi=*dsigma*)
!
      do jlev = 1 , NLEV
         spt = spt + dsigma(jlev) * sdt(:,jlev)
      enddo
!
!*    3. calculate temperature tendency from eq.14
!
      do jlev = 1 , NLEV
       do jsp = 1 , NSPP
        stt(jsp,jlev)=stt(jsp,jlev)-dot_product(tau(:,jlev),sdt(jsp,:))
       enddo
      enddo
!
!*    4. time stepping and time filtering (1st part)
!        the time filtering is spited into two parts:
!
!        1st part: xf(prel)=x+eps*(xf(-)-2x)
!        2nd part: xf=xf(prel)+eps*x(+)
!
!        together the complete filter is xf=x+eps(xf(-)-2x+x(+))
!
!        with: x        = the prognostic variable to be filtered
!              xf       = the filtered variable (current time step)
!              xf(-)    = the filtered variable (old time step)
!              xf(prel) = a preliminary filterd variable
!              x(+)     = x at time (t+2dt)
!              eps      = filter constant = *pnu*
!
!     4.a 1st part of time filter
!
      if (nkits == 0) then
         spm = pnu21 * spp + pnu * apm
         sdm = pnu21 * sdp + pnu * adm
         szm = pnu21 * szp + pnu * azm
         stm = pnu21 * stp + pnu * atm
         sqm = pnu21 * sqp + pnu * aqm
      endif
!
!     4.b add tendencies
!
      spp = apm - delt2 * spt   ! log surface pressure (negative tendencies)
      sdp =   2.0 * sdt - adm   ! note that sdt=adm+delt*tendencies (see 1.)
      stp = delt2 * stt + atm   ! temperature
      sqp = delt2 * sqt + aqm   ! spec. humidity
      szp = delt2 * szt + azm   ! vorticity
!
!     conserve p
!
      if (mypid == NROOT) then
         spp(1) = 0.0
         spp(2) = 0.0
      endif
!
!*    now the diabatic time step is finished beside the time filtering.
!     2nd part of time filtering is done in subroutine spectrald
!
!     finaly the partial arrays are gathered from all processors (mpi)
!
      call mpgallsp(sd,sdp,NLEV)
      call mpgallsp(sz,szp,NLEV)
      call mpgallsp(st,stp,NLEV)
      call mpgallsp(sq,sqp,NLEV)
      call mpgallsp(sp,spp,   1)
!
!     franks diagnostic
!
      if(ndiagsp==1) then
       allocate(ztt(NESP,NLEV))
       call mpgasp(ztt,stt,NLEV)
       dsp3d(:,1:NLEV,3)=ztt(:,1:NLEV)*ct*ww
       deallocate(ztt)
      endif
!
      return
      end

!-Diabatic-subroutines  (Frank (Larry) 26-Nov-99)

!     =====================
!     SUBROUTINE GRIDPOINTD
!     =====================

      subroutine gridpointd
      use pumamod
      use calmod, only : monaccu
      use seamod, only : naccua

      use landmod, only : NLSOIL
      use oceanmod, only : NLEV_OCE

      implicit none   ! added, FastOpt
      integer :: jlev ! added, FastOpt
!
      real sdf(NESP,NLEV)
      real stf(NESP,NLEV)
      real sqf(NESP,NLEV)
      real szf(NESP,NLEV)
!
!*    Diabatic Gridpoint Calculations
!
      gudt(:,:)=0.
      gvdt(:,:)=0.
      gtdt(:,:)=0.
      gqdt(:,:)=0.
      dudt(:,:)=0.
      dvdt(:,:)=0.
      dtdt(:,:)=0.
      dqdt(:,:)=0.

!
!     transform to gridpoint domain
!
      call invlegd

      call fc2gp(gu  ,NLON,NLPP*NLEV)
      call fc2gp(gv  ,NLON,NLPP*NLEV)
      call fc2gp(gt  ,NLON,NLPP*NLEV)
      call fc2gp(gq  ,NLON,NLPP*NLEV)
      call fc2gp(gp  ,NLON,NLPP)
!$taf store gq,gp = inner1, REC = 1 + nstep2
!
!     dimensionalize
!
      dp(:)=psurf*exp(gp(:))
!$taf INCOMPL du,dv,dt,dq
      do jlev = 1 , NLEV
       du(:,jlev)=cv*gu(:,jlev)*SQRT(rcsq(:))
       dv(:,jlev)=cv*gv(:,jlev)*SQRT(rcsq(:))
       dt(:,jlev)=ct*(gt(:,jlev)+t0(jlev))
       dq(:,jlev)=gq(:,jlev)*psurf/dp(:)
      enddo
!$taf store du,dv,dt,dq,dp = inner1, REC = 1 + nstep2
!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        print*,'-------------- SUBROUTINE GRIDPOINTD ---------------'
        print*,'NSTEP= ',nstep
       endif
       call prdbug1
      endif

!
!     PARAMETERIZATION ROUTINES TO BE INCLUDED BELOW
!
!     a) miscellaneous parameterizations (like q-fixer)
!
!$taf incompl dtdt,dq
      call miscstep
!$taf store dq,drhs,dtsa,dz0 = inner1, REC = 1 + nstep2


!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        print*,'After MISC:'
       endif
       call prdbug2
      endif

!
!     b) surface fluxes and vertical diffusion
!

      call fluxstep
!$taf store devap,dlhfl,dqdt,dshfl,dtaux,dtauy,dudt,dvdt &
!$taf& = inner1, REC = 1 + nstep2

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        print*,'After FLUX:'
       endif
       call prdbug2
      endif

!
!     c) radiation
!

      if(nrad > 0) call radstep
!$taf store dtdt,dflux,dswfl = inner1, REC = 1 + nstep2

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        print*,'After RAD:'
       endif
       call prdbug2
      endif

!
!     d) convective and large scale rain and clouds
!

      call rainstep
!$taf store dp,dprc,dprl,dprs,dq,dqdt = inner1, REC = 1 + nstep2

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        print*,'After RAIN:'
       endif
       call prdbug2
      endif

!
!     e) surface parmeterizations or models
!

      call surfstep

!
!     END OF PARAMETERISATION ROUTINES
!

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        print*,'FINAL:'
       endif
       call prdbug2
      endif

!
!     franks diagnostic
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEV,1)=dtdt(:,1:NLEV)
       do jlev=1,NLEV
        dgp3d(:,jlev,10)=dqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000
       enddo
      endif

!
!     de-dimensionalize
!

      do jlev = 1 , NLEV
       gudt(:,jlev)=dudt(:,jlev)/SQRT(rcsq(:))/cv/ww                    &
     &             +gudt(:,jlev)
       gvdt(:,jlev)=dvdt(:,jlev)/SQRT(rcsq(:))/cv/ww                    &
     &             +gvdt(:,jlev)
       gtdt(:,jlev)=dtdt(:,jlev)/ct/ww                                  &
     &             +gtdt(:,jlev)
       gqdt(:,jlev)=dqdt(:,jlev)*dp(:)/ww/psurf                         &
     &             +gqdt(:,jlev)
      enddo

!
!     transform to spectral space
!
!$taf INCOMPL gtdt
      call gp2fc(gtdt,NLON,NLPP*NLEV)
!$taf INCOMPL gqdt
      call gp2fc(gqdt,NLON,NLPP*NLEV)
!$taf INCOMPL gudt
      call gp2fc(gudt,NLON,NLPP*NLEV)
!$taf INCOMPL gvdt
      call gp2fc(gvdt,NLON,NLPP*NLEV)

!
!     direct Legendre transformation (fourier domain to spectral domain)
!

      call dirlegd(stf,sqf,gtdt,gqdt)
      call uv2dv(gudt,gvdt,sdf,szf)

      call mpsumsc(stf,stt,NLEV)
      call mpsumsc(sqf,sqt,NLEV)
      call mpsumsc(sdf,sdt,NLEV)
      call mpsumsc(szf,szt,NLEV)

      return
      end

!     ==================
!     SUBROUTINE PRDBUG1
!     ==================

      subroutine prdbug1
      use pumamod

      integer kmaxl(1),kminl(1)

      real zfpr1(NLON*NLAT,NLEV)
      real zfpr2(NLON*NLAT,NLEV)
      real zfpr3(NLON*NLAT,NLEV)
      real zfpr4(NLON*NLAT,NLEV)
      real zfpr5(NLON*NLAT,NLEV)
      real zfpr6(NLON*NLAT,NLEV)

      call mpgagp(zfpr1,dt,NLEV)
      call mpgagp(zfpr2,dq,NLEV)
      call mpgagp(zfpr3,du,NLEV)
      call mpgagp(zfpr4,dv,NLEV)
      call mpgagp(zfpr5,dcc,NLEV)
      call mpgagp(zfpr6,dql,NLEV)

      if(mypid==NROOT) then
       if(nprint==1) then
        print*,'Global diagnostics: '
        do jlev=1,NLEV
         zmaxv=MAXVAL(zfpr1(:,jlev))
         kmaxl=MAXLOC(zfpr1(:,jlev))
         zminv=MINVAL(zfpr1(:,jlev))
         kminl=MINLOC(zfpr1(:,jlev))
         print*,'L= ',jlev,' MAX T= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN T= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr2(:,jlev))
         kmaxl=MAXLOC(zfpr2(:,jlev))
         zminv=MINVAL(zfpr2(:,jlev))
         kminl=MINLOC(zfpr2(:,jlev))
         print*,'L= ',jlev,' MAX Q= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN Q= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr3(:,jlev))
         kmaxl=MAXLOC(zfpr3(:,jlev))
         zminv=MINVAL(zfpr3(:,jlev))
         kminl=MINLOC(zfpr3(:,jlev))
         print*,'L= ',jlev,' MAX U= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN U= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr4(:,jlev))
         kmaxl=MAXLOC(zfpr4(:,jlev))
         zminv=MINVAL(zfpr4(:,jlev))
         kminl=MINLOC(zfpr4(:,jlev))
         print*,'L= ',jlev,' MAX V= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN V= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr5(:,jlev))
         kmaxl=MAXLOC(zfpr5(:,jlev))
         zminv=MINVAL(zfpr5(:,jlev))
         kminl=MINLOC(zfpr5(:,jlev))
         print*,'L= ',jlev,' MAX C= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN C= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr6(:,jlev))
         kmaxl=MAXLOC(zfpr6(:,jlev))
         zminv=MINVAL(zfpr6(:,jlev))
         kminl=MINLOC(zfpr6(:,jlev))
         print*,'L= ',jlev,' MAX L= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN L= ',zminv,' NHOR= ',kminl(1)
        enddo
       elseif(nprint==2) then
        print*,'Local diagnostics at nhor= ',nprhor,': '
        do jlev=1,NLEV
         print*,'L= ',jlev,' T= ',zfpr1(nprhor,jlev)                    &
     &                    ,' Q= ',zfpr2(nprhor,jlev)
         print*,'L= ',jlev,' U= ',zfpr3(nprhor,jlev)                    &
     &                    ,' V= ',zfpr4(nprhor,jlev)
         print*,'L= ',jlev,' C= ',zfpr5(nprhor,jlev)                    &
     &                    ,' L= ',zfpr6(nprhor,jlev)
        enddo
       endif
      endif

      call mpgagp(zfpr1(1,1),dls,1)
      call mpgagp(zfpr2(1,1),dp,1)
      call mpgagp(zfpr3(1,1),drhs,1)
      call mpgagp(zfpr4(1,1),dalb,1)
      call mpgagp(zfpr4(1,2),dz0,1)
      call mpgagp(zfpr5(1,1),dt(1,NLEP),1)
      call mpgagp(zfpr6(1,1),dq(1,NLEP),1)

      if(mypid==NROOT) then
       if(nprint==1) then
        print*,'Global diagnostic:'
        zmaxv=MAXVAL(zfpr5(:,1))
        kmaxl=MAXLOC(zfpr5(:,1))
        zminv=MINVAL(zfpr5(:,1))
        kminl=MINLOC(zfpr5(:,1))
        print*,'SURF MAX T= ',zmaxv,' NHOR= ',kmaxl(1)
        print*,'SURF MIN T= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr6(:,1))
        kmaxl=MAXLOC(zfpr6(:,1))
        zminv=MINVAL(zfpr6(:,1))
        kminl=MINLOC(zfpr6(:,1))
        print*,'SURF MAX Q= ',zmaxv,' NHOR= ',kmaxl(1)
        print*,'SURF MIN Q= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr2(:,1))
        kmaxl=MAXLOC(zfpr2(:,1))
        zminv=MINVAL(zfpr2(:,1))
        kminl=MINLOC(zfpr2(:,1))
        print*,'SURF MAX P= ',zmaxv,' NHOR= ',kmaxl(1)
        print*,'SURF MIN P= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr3(:,1))
        kmaxl=MAXLOC(zfpr3(:,1))
        zminv=MINVAL(zfpr3(:,1))
        kminl=MINLOC(zfpr3(:,1))
        print*,'SURF MAX RHS= ',zmaxv,' NHOR= ',kmaxl(1)
        print*,'SURF MIN RHS= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr4(:,1))
        kmaxl=MAXLOC(zfpr4(:,1))
        zminv=MINVAL(zfpr4(:,1))
        kminl=MINLOC(zfpr4(:,1))
        print*,'SURF MAX ALB= ',zmaxv,' NHOR= ',kmaxl(1)
        print*,'SURF MIN ALB= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr4(:,2))
        kmaxl=MAXLOC(zfpr4(:,2))
        zminv=MINVAL(zfpr4(:,2))
        kminl=MINLOC(zfpr4(:,2))
        print*,'SURF MAX Z0= ',zmaxv,' NHOR= ',kmaxl(1)
        print*,'SURF MIN Z0= ',zminv,' NHOR= ',kminl(1)
       elseif(nprint==2) then
        print*,'Local diagnostics at nhor= ',nprhor,': '
        print*,'SURF: LS= ',zfpr1(nprhor,1)                             &
     &             ,' PS= ',zfpr2(nprhor,1)
        print*,'SURF: RHS= ',zfpr3(nprhor,1)                            &
     &             ,' ALB= ',zfpr4(nprhor,1)
        print*,'SURF: TS= ',zfpr5(nprhor,1)                             &
     &             ,' QS= ',zfpr6(nprhor,1)
        print*,'SURF: Z0= ',zfpr4(nprhor,2)
       endif
      endif

      if(nprint==2) then
       call mpgagp(zfpr1(1,1),dicec,1)

       if(mypid==NROOT) then
        print*,'SURF: ICEC= ',zfpr1(nprhor,1)
       endif
      endif

      return
      end subroutine prdbug1

!     ==================
!     SUBROUTINE PRDBUG2
!     ==================

      subroutine prdbug2
      use pumamod

      integer kmaxl(1),kminl(1)

      real zfpr1(NLON*NLAT,NLEV)
      real zfpr2(NLON*NLAT,NLEV)
      real zfpr3(NLON*NLAT,NLEV)
      real zfpr4(NLON*NLAT,NLEV)

      call mpgagp(zfpr1,dtdt,NLEV)
      call mpgagp(zfpr2,dqdt,NLEV)
      call mpgagp(zfpr3,dudt,NLEV)
      call mpgagp(zfpr4,dvdt,NLEV)

      if(mypid==NROOT) then
       if(nprint==1) then
        print*,'Global diagnostics: '
        do jlev=1,NLEV
         zmaxv=MAXVAL(zfpr1(:,jlev))
         kmaxl=MAXLOC(zfpr1(:,jlev))
         zminv=MINVAL(zfpr1(:,jlev))
         kminl=MINLOC(zfpr1(:,jlev))
         print*,'L= ',jlev,' MAX DT= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN DT= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr2(:,jlev))
         kmaxl=MAXLOC(zfpr2(:,jlev))
         zminv=MINVAL(zfpr2(:,jlev))
         kminl=MINLOC(zfpr2(:,jlev))
         print*,'L= ',jlev,' MAX DQ= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN DQ= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr3(:,jlev))
         kmaxl=MAXLOC(zfpr3(:,jlev))
         zminv=MINVAL(zfpr3(:,jlev))
         kminl=MINLOC(zfpr3(:,jlev))
         print*,'L= ',jlev,' MAX DU= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN DU= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr4(:,jlev))
         kmaxl=MAXLOC(zfpr4(:,jlev))
         zminv=MINVAL(zfpr4(:,jlev))
         kminl=MINLOC(zfpr4(:,jlev))
         print*,'L= ',jlev,' MAX DV= ',zmaxv,' NHOR= ',kmaxl(1)
         print*,'L= ',jlev,' MIN DV= ',zminv,' NHOR= ',kminl(1)
        enddo
       elseif(nprint==2) then
        print*,'Local diagnostics at nhor= ',nprhor,': '
        do jlev=1,NLEV
         print*,'L= ',jlev,' DT= ',zfpr1(nprhor,jlev)                   &
     &                    ,' DQ= ',zfpr2(nprhor,jlev)
         print*,'L= ',jlev,' DU= ',zfpr3(nprhor,jlev)                   &
     &                    ,' DV= ',zfpr4(nprhor,jlev)
        enddo
       endif
      endif

      return
      end subroutine prdbug2

!     ====================
!     SUBROUTINE SPECTRALD
!     ====================

      subroutine spectrald
      use pumamod
      use observation, only : getnudg, getnudg2d
!
!     franks diagnostics
!
      real, allocatable :: ztt(:,:)
      real :: frd(NUGP,NLEV)   = 0.0     ! temperature 
      real :: frdt(NHOR,NLEV)   = 0.0     ! temperature
      real :: srdt(NESP,NLEV)   = 0.0     ! temperature
      real :: srdtp(NSPP,NLEV)  = 0.0  ! Spectral Temperature Partial

!
!     add tendencies from diabatic parameterizations
!

      szp = szp + delt2 * szt
      stp = stp + delt2 * stt
      sqp = sqp + delt2 * sqt
      sdp = sdp + delt2 * sdt

!
!     franks diagnostic
!

      if(ndiagsp==1) then
       allocate(ztt(NESP,NLEV))
       call mpgasp(ztt,stt,NLEV)
       dsp3d(:,1:NLEV,1)=ztt(:,1:NLEV)*ct*ww
       deallocate(ztt)
       call mkdqtgp
       do jlev=1,NLEV
        dgp3d(:,jlev,18)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      endif

!     calculates spectral tendencies from restoration (if included)
!     and biharmonic diffusion (both implicitely).
!     add newtonian cooling and drag.
!     add in biharmonic difusion if required


      do jlev = 1 , NLEV
       stt(:,jlev)=((-damp(jlev)-tdisst(jlev)*sakpp(1:NSPP,jlev))       &
     &              *stp(:,jlev)+damp(jlev)*srp(:,jlev))                &
     &          /(1.+delt2*(damp(jlev)+tdisst(jlev)*sakpp(1:NSPP,jlev)))
       sqt(:,jlev)=-tdissq(jlev)*sakpp(1:NSPP,jlev)*sqp(:,jlev)         &
     &            /(1.+delt2*tdissq(jlev)*sakpp(1:NSPP,jlev))
       sdt(:,jlev)=((-tfrc(jlev)-tdissd(jlev)*sakpp(1:NSPP,jlev))       &
     &              *sdp(:,jlev))                                       &
     &          /(1.+delt2*(tfrc(jlev)+tdissd(jlev)*sakpp(1:NSPP,jlev)))
       szt(:,jlev)=((-tfrc(jlev)-tdissz(jlev)*sakpp(1:NSPP,jlev))       &
     &              *szp(:,jlev))                                       &
     &          /(1.+delt2*(tfrc(jlev)+tdissz(jlev)*sakpp(1:NSPP,jlev)))
      enddo

!     no friction on planetary vorticity
!     no difusion on EZ (planetary vorticity)

      if (mypid == NROOT) then
         spp(1) = 0.0
         spp(2) = 0.0
         szt(3,:) = szt(3,:) + EZ * (tfrc(:) + tdissz(:)*sakpp(3,:))    &
     &                       /(1.+delt2*(tfrc(:)+tdissz(:)*sakpp(3,:)))
      endif

      szp = szp + delt2 * szt
      stp = stp + delt2 * stt
      sqp = sqp + delt2 * sqt
      sdp = sdp + delt2 * sdt
!$taf store sdp,szp,stp,sqp = inner1, REC = 1 + nstep2
      if (nkits == 0) then
         szm = szm + pnu * szp
         stm = stm + pnu * stp
         sqm = sqm + pnu * sqp
         sdm = sdm + pnu * sdp
         spm = spm + pnu * spp
      endif


if (makenudge==1) then
! nudging on divergence
     frd=0. ! reset frdt
     frdt=0. ! reset frdt
     srdt=0.
     srdtp=0.
     if(mypid==NROOT) then
      call getnudg(nstep,155,frd) ! nudging for former step:divergence
     endif
      call mpscgp(frd,frdt,NLEV)
      frdt=frdt/ww  ! nondimesionalize
      do jlev=1,NLEV
        call gp2fc(frdt(:,jlev),NLON,NLPP)  
        call fc2sp(frdt(:,jlev),srdt(:,jlev))
      enddo
       call mpsumsc(srdt,srdtp,NLEV)

     do jlev=1,NLEV 
! implicit calculated
       sdp(:,jlev)= sdp(:,jlev)+ delt2*dnudg(jlev) *(srdtp(:,jlev)-sdp(:,jlev))/(1.+dnudg(jlev) *delt2) 
     enddo 

!nudging on vorticity
     frd=0. ! reset frdt
     frdt=0. ! reset frdt
     srdt=0.
     srdtp=0.
     if(mypid==NROOT) then
     call getnudg(nstep,138,frd) ! nudging for former step:vorticity
     endif
     call mpscgp(frd,frdt,NLEV)
     frdt=frdt/ww    ! nondimensionalize      
     do jlev=1,NLEV
       call gp2fc(frdt(:,jlev),NLON,NLPP)  
       call fc2sp(frdt(:,jlev),srdt(:,jlev))
 !      srdt(3,jlev)=srdt(3,jlev)+EZ    
     enddo     
    call mpsumsc(srdt,srdtp,NLEV)
    if (mypid == NROOT) then
         srdtp(3,:) = srdtp(3,:) +EZ ! transform to relative vorticity.
    endif

    
     do jlev=1,NLEV
!implicit nudging
       szp(:,jlev)= szp(:,jlev)+delt2*znudg(jlev) *(srdtp(:,jlev)-szp(:,jlev))/(1.+znudg(jlev) *delt2)
     enddo

!! nudging on temperature
      frd=0. ! reset frdt
      frdt=0. ! reset frdt
      srdt=0.
      srdtp=0.
     if(mypid==NROOT) then
      call getnudg(nstep,130,frd) ! nudging for former step:temperature
     endif      
     call mpscgp(frd,frdt,NLEV)
      do jlev=1,NLEV
        frdt(:,jlev)=frdt(:,jlev)/ct-t0(jlev)  ! nondimesionalize
      end do
      do jlev=1,NLEV
        call gp2fc(frdt(:,jlev),NLON,NLPP)  
        call fc2sp(frdt(:,jlev),srdt(:,jlev))
      enddo      

      call mpsumsc(srdt,srdtp,NLEV)
      do jlev=1,NLEV
!implicit nudging
        stp(:,jlev)= stp(:,jlev)+delt2*tnudg(jlev) *(srdtp(:,jlev)-stp(:,jlev))/(1.+tnudg(jlev) *delt2)
      enddo

!nudging on moist 
      frd=0. ! reset frdt
      frdt=0. ! reset frdt
      srdt=0.
      srdtp=0.
     if(mypid==NROOT) then
      call getnudg(nstep,133,frd) ! nudging for former step:moist
      endif
        call mpscgp(frd,frdt,NLEV)

     do jlev=1,NLEV
       frdt(:,jlev)=frdt(:,jlev)*dp(:)/psurf  ! nondimesionalize
      ENDDO
      do jlev=1,NLEV
         call gp2fc(frdt(:,jlev),NLON,NLPP)  
         call fc2sp(frdt(:,jlev),srdt(:,jlev))
      enddo     
      call mpsumsc(srdt,srdtp,NLEV)

      do jlev=1,NLEV   
!implicit nudging
       sqp(:,jlev)= sqp(:,jlev)+delt2*qnudg(jlev) *(srdtp(:,jlev)-sqp(:,jlev))/(1.+qnudg(jlev) *delt2)
      enddo

! nudging on pressure
      frdt=0. ! reset frdt
      srdt=0.
      srdtp=0.
      if(mypid==NROOT) then
      call getnudg2d(nstep,134,frd(:,1)) ! nudging for former step:pressure
      endif
      call mpscgp(frd(:,1),frdt(:,1),1)
      frdt(:,1)=log(frdt(:,1)/psurf)
      call gp2fc(frdt(:,1) ,NLON,NLPP)  ! fft
      call fc2sp(frdt(:,1),srdt(:,1))
      
      call mpsumsc(srdt,srdtp,NLEV)
      spp(:)=spp(:)+delt2*pnudg *(srdtp(:,1)-spp(:))/(1.+pnudg *delt2)
      if (mypid == NROOT) then
         spp(1) = 0.0
         spp(2) = 0.0
      endif

endif

      call mpgallsp(sd,sdp,NLEV)
      call mpgallsp(sz,szp,NLEV)
      call mpgallsp(st,stp,NLEV)
      call mpgallsp(sq,sqp,NLEV)
      call mpgallsp(sp,spp,   1)

!
!     franks diagnostic
!

      if(ndiagsp==1) then
       allocate(ztt(NESP,NLEV))
       call mpgasp(ztt,stt,NLEV)
       dsp3d(:,1:NLEV,2)=ztt(:,1:NLEV)*ct*ww
       deallocate(ztt)
       call mkdqtgp
       do jlev=1,NLEV
        dgp3d(:,jlev,19)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      endif

      return
      end

!     ===============
!     SUBROUTINE MKDQ
!     ===============

      subroutine mkdq(psq,pgq,KLEV)
      use pumamod
!
      real psq(NESP,klev),pgq(NHOR,klev)
!
      real zps(NHOR)
!
      call sp2fl(psq,pgq,KLEV)
      call fc2gp(pgq,NLON,NLPP*KLEV)
      call sp2fl(sp,zps,1)
      call fc2gp(zps,NLON,NLPP)
!
      zps(:)=psurf*EXP(zps(:))
      do jlev=1,KLEV
       pgq(:,jlev)=pgq(:,jlev)*psurf/zps(:)
      enddo
!
      return
      end subroutine mkdq

!     ===================
!     SUBROUTINE MKSECOND
!     ===================

      subroutine mksecond(ps,poff)
!
      real ps,poff
      integer ivalues(8)
!
      call date_and_time(values=ivalues)
!
      zdays = ivalues(3) - 1
      zhour = 24.0 * zdays + ivalues(5)
      zmins = 60.0 * zhour + ivalues(6)
      zsecs = 60.0 * zmins + ivalues(7) + 0.001 * ivalues(8)
      ps = zsecs - poff
!
      return
      end
