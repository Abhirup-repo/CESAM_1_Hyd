      module rainmod
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: version = '22.11.2002 by Larry'
!
!     namelist parameter:
!

      integer :: kbetta = 1  ! switch for betta in kuo (1/0=yes/no)
      integer :: nprl   = 1  ! switch for large scale precip (1/0=yes/no)
      integer :: nprc   = 1  ! switch for large convective precip (1/0=yes/no)
      integer :: ndca   = 1  ! switch for dry convective adjustment (1/0=yes/no)
      integer :: ncsurf = 1  ! switch for conv. starts from surface (1/0=yes/no)
      integer :: nmoment= 0  ! switch for momentum mixing (1/0=yes/no)

      real :: clwcrit1 = -0.1 ! 1st critical vertical velocity for clouds
      real :: clwcrit2 =  0.0 ! 2nd critical vertical velocity for clouds
                              ! not active if clwcrit1 < clwcrit2

      real :: rcrit(NLEV)    ! critical relative hum. for non conv. clouds
!
!     global scalras
!

      real :: clwfac    = 0. ! smothing for cloud suppression
      real :: time4rain = 0. ! CPU time needed for rainmod
      real :: time4prl  = 0. ! CPU time spent for large scale precip
      real :: time4prc  = 0. ! CPU time spent for convective precip
      real :: time4dca  = 0. ! CPU time spent dry convection
      real :: time4cl   = 0. ! CPU time spent for cloud formation

!
!     global arrays
!

      integer :: icclev(NHOR,NLEV)  ! flag for convective cloud layer
      integer :: icctot(NHOR)       ! total convective cloud layers

      end module rainmod

!     ==================
!     SUBROUTINE RAININI
!     ==================

      subroutine rainini
      use rainmod
!
      integer :: im
      namelist/rainpar/kbetta,nprl,nprc,ndca,ncsurf,nmoment             &
     &                ,rcrit,clwcrit1,clwcrit2
!
!     set default rcrit
!
      rcrit(:)=MAX(0.85,MAX(sigma(:),1.-sigma(:)))
!    guokun initialize 2d varriable
     do im=1,NUGP
       rcritg(im,1:NLEV)=rcrit(1:NLEV)
     enddo


!
      if(mypid==NROOT) then
         read(11,rainpar)
         write (*,'(/," *****************************************")')
         write (*,'(" * RAINMOD ",a29," *")') trim(version)
         write (*,'(" *****************************************")')
         write (*,'(" * Namelist RAINPAR from <puma_namelist> *")')
         write (*,'(" *****************************************")')
         write(*,rainpar)
      endif
!
!     broadcast namelist parameter
!
      call mpbci(kbetta)
      call mpbci(nprl)
      call mpbci(nprc)
      call mpbci(ndca)
      call mpbci(ncsurf)
      call mpbci(nmoment)
      call mpbcr(clwcrit1)
      call mpbcr(clwcrit2)
      call mpbcrn(rcrit,NLEV)
      call mpscgp(rcritg,rcritm,NLEV)
!
!     smoothing factor for cloud suppression
!
      if(clwcrit1 > clwcrit2) then
       clwfac=1./(clwcrit1-clwcrit2)
      elseif(clwcrit1 < clwcrit2) then
       clwfac=-1.
      else
       clwfac=1.
      endif
!
      dprl(:)=0.
      dprc(:)=0.
      dprs(:)=0.

      return
      end subroutine rainini

!     ===================
!     SUBROUTINE RAINSTEP
!     ===================

      subroutine rainstep
      use rainmod
!
!**   performance estimate
!
      if(ntime==1) call mksecond(zsec,0.)
!
!
!**   0) preset precipitation
!

      dprl(:)=0.
      dprc(:)=0.
      dprs(:)=0.

!
!**   1) transfer adiabatic q-tendencies to grid point space
!        (needed for kuo scheme)
!

      if(nprc == 1) call mkdqtgp
!$TAF STORE dqt,dp  = inner1, REC = 1 + nstep2
!
!**   2) kuo convection
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(nprc == 1) call kuo
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4prc=time4prc+zsec1
      endif
!$TAF STORE dprc,dprl,dq,dt,dw,icclev,icctot  = inner1, REC = 1 + nstep2
!
!**   3) dry convective adjustment
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(ndca == 1) call mkdca
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4dca=time4dca+zsec1
      endif
!$TAF STORE dqdt,dtdt = inner1, REC = 1 + nstep2
!
!**   4) large scale precipitation
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(nprl == 1) call mklsp
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4prl=time4prl+zsec1
      endif
!$TAF STORE dqdt,dtdt = inner1, REC = 1 + nstep2
!
!**   5) diagnose cloud cover
!

      if(ntime ==1) call mksecond(zsec1,0.)
      call mkclouds
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4cl=time4cl+zsec1
       call mksecond(zsec,zsec)
       time4rain=time4rain+zsec
      endif
!
      return
      end subroutine rainstep

!     ===================
!     SUBROUTINE RAINSTOP
!     ===================

      subroutine rainstop
      use rainmod
!
      if(mypid == NROOT .and. ntime == 1) then
       print*,'*********************************************'
       print*,' CPU usage in RAINSTEP (ROOT process only):  '
       print*,'    All routines   : ',time4rain,' s'
       print*,'    Large scale p  : ',time4prl,' s'
       print*,'    Convective p   : ',time4prc,' s'
       print*,'    Cloud formation: ',time4cl,' s'
       print*,'    Dry convection : ',time4dca,' s'
       print*,'*********************************************'
      endif
!
      return
      end subroutine rainstop


!     ==================
!     SUBROUTINE MKDQTGP
!     ==================

      subroutine mkdqtgp
      use rainmod
!
      real zsqt(NESP,NLEV)
!
      call mpgallsp(zsqt,sqt,NLEV)
!$taf incompl dqt
      call sp2fl(zsqt,dqt,NLEV)
      call fc2gp(dqt,NLON,NLPP*NLEV)
!
      do jlev=1,NLEV
       dqt(:,jlev)=dqt(:,jlev)*psurf*ww/dp(:)
      enddo
!
      return
      end subroutine mkdqtgp

!     ================
!     SUBROUTINE MKLSP
!     ================

      subroutine mklsp
      use rainmod
!
!     local arrays
!

      real zq(NHOR)          ! preliminary q at t+dt
      real zt(NHOR)          ! preliminary t at t+dt
! FastOpt      real zqn(NHOR)         ! q after condensation
      real :: zqn(NHOR) = 0. ! FastOpt
! FastOpt      real ztn(NHOR)         ! t after condensation
      real :: ztn(NHOR) = 0. ! FastOpt
      real zdqdt(NHOR)       ! local q-tendency
      real zdtdt(NHOR)       ! local t-tendency
      real zqsat(NHOR)       ! saturation humidity
      real zlcpe(NHOR)       ! L/CP
      real zcor(NHOR)        ! correction term for saturation humidity
      logical :: l_cond(NHOR) ! FastOpt
!
!     allocatable arrys for dbug diagnostic
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)

!
!*    set some const.
!

      zdelt=delt2/ww

!
!*    make condensation
!
!$taf LOOP = PARALLEL
      do jlev=1,NLEV

!
!     preset tendencies
!

       zdqdt(:)=0.
       zdtdt(:)=0.

!
!     compute preliminary t and q at t+dt
!

       zq(:)=dq(:,jlev)+dqdt(:,jlev)*zdelt
       zt(:)=dt(:,jlev)+dtdt(:,jlev)*zdelt
!$taf store zq = innern, rec = jlev + nstep2*NLEV
!
!*    set l/cp depending on t
!

       where(zt(:)>TMELT)
        zlcpe(:)=ALV/(acpd*(1.+ADV*zq(:)))
       elsewhere
        zlcpe(:)=ALS/(acpd*(1.+ADV*zq(:)))
       endwhere

!
!     saturation humidity
!

       zqsat(:)=rdbrv*ra1*exp(ra2*(zt(:)-TMELT)/(zt(:)-ra4))            &
     &         /(dp(:)*sigma(jlev))
       zcor(:)=1./(1.-(1./rdbrv-1.)*zqsat(:))
       zqsat(:)=zqsat(:)*zcor(:)

!      Fix for problem with higher Mars atmosphere

       where (zqsat(:) < 0.0)
          zqsat(:) = 0.0
       endwhere

!
!     look for supersaturation
!
!$taf store zqsat,zcor,zlcpe = innern, rec = jlev + nstep2*NLEV
!$taf incompl dqdt,dtdt,dprl
! FastOpt       where(zq(:) > zqsat(:))
       l_cond(:) = (zq(:) > zqsat(:)) ! FastOpt
       where (l_cond(:)) ! FastOpt
!
!     make condensation (2 iterations)
!

        zdqdt(:)=(zqsat(:)-zq(:))                                       &
     &          /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                          &
     &           *zqsat(:)*zcor(:)/(zt(:)-ra4)**2)
        ztn(:)=zt(:)-zdqdt(:)*zlcpe(:)
        zqn(:)=zq(:)+zdqdt(:)
        zqsat(:)=rdbrv*ra1*exp(ra2*(ztn(:)-TMELT)/(ztn(:)-ra4))         &
     &          /(dp(:)*sigma(jlev))
        zcor(:)=1./(1.-(1./rdbrv-1.)*zqsat(:))
     end where ! FastOpt
!$taf store zqsat,zcor,zqn,ztn = innern, rec = jlev + nstep2*NLEV
     where (l_cond) ! Fastopt
        zqsat(:)=zqsat(:)*zcor(:)
        zdqdt(:)=(zqsat(:)-zqn(:))                                      &
     &          /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                          &
     &           *zqsat(:)*zcor(:)/(ztn(:)-ra4)**2)
        ztn(:)=ztn(:)-zdqdt(:)*zlcpe(:)
        zqn(:)=zqn(:)+zdqdt(:)

!
!     final tendencies
!

        zdqdt(:)=(zqn(:)-zq(:))/zdelt
        zdtdt(:)=(ztn(:)-zt(:))/zdelt

!
!     add tendencies
!

        dqdt(:,jlev)=dqdt(:,jlev)+zdqdt(:)
        dtdt(:,jlev)=dtdt(:,jlev)+zdtdt(:)

!
!     add precipitation
!

        dprl(:)=dprl(:)-AMIN1(zdqdt(:),0.0)*dsigma(jlev)

       endwhere

!
!     dbug print out
!

       if(nprint==2) then
        allocate(zprf1(NLON*NLAT))
        allocate(zprf2(NLON*NLAT))
        allocate(zprf3(NLON*NLAT))
        allocate(zprf4(NLON*NLAT))
        allocate(zprf5(NLON*NLAT))
        allocate(zprf6(NLON*NLAT))
        call mpgagp(zprf1,zq,1)
        call mpgagp(zprf2,zt,1)
        call mpgagp(zprf3,zqsat,1)
        call mpgagp(zprf4,zdqdt,1)
        call mpgagp(zprf5,zdtdt,1)
        call mpgagp(zprf6,dprl,1)
        if(mypid==NROOT) then
         if(jlev==1) print*,'In mklsp:'
         print*,'L= ',jlev,' zt= ',zprf2(nprhor)                        &
     &         ,' zq= ',zprf1(nprhor),' zqsat= ',zprf3(nprhor)
         print*,'L= ',jlev,' dtdt= ',zprf5(nprhor)                      &
     &         ,' dqdt= ',zprf4(nprhor)
         print*,'L= ',jlev,' prl (mm/day)= ',zprf6(nprhor)              &
     &                                      *0.5*zdelt*ntspd*1000.
        endif
        deallocate(zprf1)
        deallocate(zprf2)
        deallocate(zprf3)
        deallocate(zprf4)
        deallocate(zprf5)
        deallocate(zprf6)
       endif

!
!     franks dbug
!

       if(ndiaggp==1) then
        dgp3d(:,jlev,7)=zdtdt(:)
        dgp3d(:,jlev,14)=zdqdt(:)*dp(:)*dsigma(jlev)/ga/1000.
       endif

      enddo

!
!     convert precipitation into m/s and set snow fall
!
!$TAF STORE dprl(:) = inner1, REC = 1 + nstep2
      dprl(:)=dprl(:)*dp(:)/ga/1000.
      where(dt(:,NLEV) <= TMELT) dprs(:)=dprs(:)+dprl(:)

      return
      end subroutine mklsp

!     ==============
!     SUBROUTINE KUO
!     ==============

      subroutine kuo
      use rainmod

!
!     local arrays
!

      integer ilift(NHOR)    ! lifting level (= cloud base + 1)
      integer itop(NHOR)     ! cloud top level
      real zqe(NHOR,NLEV)    ! moisture of environmental air
      real zte(NHOR,NLEV)    ! temperatur of environmental air
      real zdqtot(NHOR,NLEV) ! total moisture accession in resp. layer
      real zdqdt(NHOR,NLEV)  ! moisture diff. (parcel-environment)
      real zdtdt(NHOR,NLEV)  ! temperature diff. (parcel-environment)
      real zalcpe(NHOR,NLEV) ! help array for factors
!FastOpt      real zqp(NHOR)         ! humidity of lifted parcel
      real :: zqp(NHOR) = 0.! added, FastOpt
!FastOpt      real ztp(NHOR)         ! temperature of lifted parcel
      real :: ztp(NHOR) = 0.! added, FastOpt
      real zi(NHOR)          ! moisture supply in the collumn
      real zpt(NHOR)         ! energy surplus of the cloud (t-part)
      real zpq(NHOR)         ! energy surplus of the cloud (q-part)
      real zbetta(NHOR)      ! betta factor (new kuo scheme)
      real zpbetta(NHOR)     ! p-scaling for betta factor
      real zat(NHOR)         ! distribution factor (temperature)
      real zaq(NHOR)         ! distribution factor (moisture)
!
!     array for momentum mixing
!
      real zue(NHOR,NLEV)    ! u of environmental air
      real zve(NHOR,NLEV)    ! v of environmental air
      real zdudt(NHOR,NLEV)  ! u-momentum difference (cloud-environment)
      real zdvdt(NHOR,NLEV)  ! v-momentum difference (cloud-environment)
      real zau(NHOR)         ! distribution factor (momentum)
      real zup(NHOR)         ! reference (cloud) u
      real zvp(NHOR)         ! reference (cloud) v
      real zsums(NHOR)       ! sum(dsigma) for mixed levels


!
!     allocatable arrays for dbug diagnostic
!

      real, allocatable :: zprf1(:,:)
      real, allocatable :: zprf2(:,:)
      real, allocatable :: zprf3(:,:)
      real, allocatable :: zprf4(:,:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:)
      real, allocatable :: zprf8(:)
      real, allocatable :: zprf9(:)
      real, allocatable :: zprf10(:)
      real, allocatable :: zprf11(:)
      real, allocatable :: zprf12(:)

      zdelt=delt2/ww         ! dimensionalized timestep

!
!*    0) initialize some fields
!

      ilift(:)=-1
      itop(:)=-1
      icctot(:)=0
      zi(:)=0.
      zpt(:)=0.
      zpq(:)=0.
      zat(:)=0.
      zaq(:)=0.
      zbetta(:)=0.
      zpbetta(:)=0.
      zdtdt(:,:)=0.
      zdqdt(:,:)=0.
      icclev(:,:)=0

!
!     moisture accession (level and column)
!

      zdqtot(:,:)=(dqt(:,1:NLEV)+dqdt(:,1:NLEV))*zdelt

!
!     for original kuo (whole column) use the next loop
!
!     do jlev=1,NLEV
!      zi(:)=zi(:)+zdqtot(:,jlev)*dsigma(jlev)
!     enddo

!
!     franks dbug
!

      if(ndiaggp==1) then
       do jlev=1,NLEV
        dgp3d(:,jlev,15)=-dqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
        dgp3d(:,jlev,17)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      end if

!
!*    1) update q and t
!

      zqe(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*zdelt
      zte(:,:)=dt(:,1:NLEV)+dtdt(:,1:NLEV)*zdelt

!
!     set L/CP depending on t
!
!$taf incompl zalcpe
      where(zte(:,:) < TMELT)
       zalcpe(:,:)=ALS/(acpd*(1.+ADV*zqe(:,:)))
      elsewhere
       zalcpe(:,:)=ALV/(acpd*(1.+ADV*zqe(:,:)))
      endwhere

!
!*    2) search for lifting level
!
      if(ncsurf == 1) then

!     a) start at surface (lifting level only if dry unstable)
!

      zfac=sigma(NLEV)**akap
!$TAF LOOP = PARALLEL
      do jhor=1,NHOR
       ztnew=dt(jhor,NLEP)*zfac
       if(zdqtot(jhor,NLEV) > 0. .and. ztnew > zte(jhor,NLEV)) then

!
!     surface air specific humidity
!

        zqsnl=rdbrv*ra1*exp(ra2*(zte(jhor,NLEV)-TMELT)                  &
     &                          /(zte(jhor,NLEV)-ra4))                  &
     &       /(dp(jhor)*sigma(NLEV))
! changed, Fastopt        zqsnl=zqsnl/(1.-(1./rdbrv-1.)*zqsnl)
        zqsnlh1=zqsnl/(1.-(1./rdbrv-1.)*zqsnl) ! changed, Fastopt
        zqnew=dq(jhor,NLEP)*AMIN1(1.,zqe(jhor,NLEV)/zqsnlh1)

!
!     saturation humidity:
!

! changed, Fastopt        zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))              &
! changed, Fastopt     &       /(dp(jhor)*sigma(NLEV))
        zqsnlh2=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))            & ! changed, Fastopt
     &       /(dp(jhor)*sigma(NLEV))                                      ! changed, Fastopt
! changed, Fastopt        zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
        zcor=1./(1.-(1./rdbrv-1.)*zqsnlh2) ! changed, Fastopt
! changed, Fastopt        zqsnl=zqsnl*zcor
        zqsnlh3=zqsnlh2*zcor ! changed, Fastopt

!
!     condensation if supersaturated
!

! changed, Fastopt        if(zqsnl < zqnew) then
        if(zqsnlh3 < zqnew) then ! changed, Fastopt

!
!     update constants (ice/water phase)
!

         if(ztnew < TMELT) then
          zlcp=ALS/(acpd*(1.+ADV*zqnew))
         else
          zlcp=ALV/(acpd*(1.+ADV*zqnew))
         endif

!
!     temperature correction term
!

! changed, Fastopt         zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
         zdqsdt=ra2*(TMELT-ra4)*zqsnlh3*zcor/((ztnew-ra4)*(ztnew-ra4)) ! changed, Fastopt


!
!     update t and q
!

! changed, Fastopt         zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
         zrain=(zqnew-zqsnlh3)/(1.+zlcp*zdqsdt) ! changed, Fastopt
! changed, FastOpt         ztnew=ztnew+zlcp*zrain
         ztnewh1=ztnew+zlcp*zrain ! changed, FastOpt
! changed, Fastopt         zqnew=zqnew-zrain
         zqnewh1=zqnew-zrain ! changed, Fastopt

!
!     second iteration:
!

! changed, FastOpt         zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))             &
! changed, FastOpt     &        /(dp(jhor)*sigma(NLEV))
         zqsnlh4=rdbrv*ra1*exp(ra2*(ztnewh1-TMELT)/(ztnewh1-ra4))         & ! changed, FastOpt
     &        /(dp(jhor)*sigma(NLEV))                                     ! changed, FastOpt
! changed, Fastopt         zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
         zcor=1./(1.-(1./rdbrv-1.)*zqsnlh4) ! changed, Fastopt
! changed, Fastopt         zqsnl=zqsnl*zcor
         zqsnlh5=zqsnlh4*zcor ! changed, Fastopt
! changed, FastOpt         zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
         zdqsdt=ra2*(TMELT-ra4)*zqsnlh5*zcor/((ztnewh1-ra4)*(ztnewh1-ra4)) ! changed, FastOpt
! changed, Fastopt         zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
         zrainh=(zqnewh1-zqsnlh5)/(1.+zlcp*zdqsdt) ! changed, Fastopt
! changed, FastOpt         ztnew=ztnew+zlcp*zrain
         ztnewh2=ztnewh1+zlcp*zrainh ! changed, FastOpt
! changed, Fastopt         zqnew=zqnew-zrain
         zqnewh2=zqnewh1-zrainh ! changed, Fastopt

!
!     if convection: set lifting level and parcel q,t (zpq,ztp)
!                    calculate prelim. dq,dt (zdqdt,zdtdt)
!                    update surplus of total energy (zp)
!

         ilift(jhor)=NLEP
         itop(jhor)=NLEV
! changed, FastOpt         ztp(jhor)=ztnew
         ztp(jhor)=ztnewh2 ! changed, FastOpt
! changed, Fastopt         zqp(jhor)=zqnew
         zqp(jhor)=zqnewh2 ! changed, Fastopt
         zqsa=rdbrv*ra1*exp(ra2*(zte(jhor,NLEV)-TMELT)                  &
     &                          /(zte(jhor,NLEV)-ra4))                  &
     &       /(dp(jhor)*sigma(NLEV))
! changed, Fastopt         zqsa=zqsa/(1.-(1./rdbrv-1.)*zqsa)
         zqsah=zqsa/(1.-(1./rdbrv-1.)*zqsa) ! changed, Fastopt
! changed, Fastopt         zdqdt(jhor,NLEV)=AMAX1(zqsa-zqe(jhor,NLEV),0.)
         zdqdt(jhor,NLEV)=AMAX1(zqsah-zqe(jhor,NLEV),0.) ! changed, Fastopt
         zdtdt(jhor,NLEV)=ztp(jhor)-zte(jhor,NLEV)
         zpt(jhor)=zdtdt(jhor,NLEV)*dsigma(NLEV)/zalcpe(jhor,NLEV)
         zpq(jhor)=zdqdt(jhor,NLEV)*dsigma(NLEV)
         if(kbetta > 0 ) then
! changed, Fastopt          zbetta(jhor)=AMIN1(1.,zqe(jhor,NLEV)/zqsa)*dsigma(NLEV)
          zbetta(jhor)=AMIN1(1.,zqe(jhor,NLEV)/zqsah)*dsigma(NLEV) ! changed, Fastopt
          zpbetta(jhor)=dsigma(NLEV)
         endif

!
!     moisture accession (delete for original kuo)
!

         zi(jhor)=zdqtot(jhor,NLEV)*dsigma(NLEV)
        endif
       endif
      enddo

      endif
!
!     b) continue search for free atmosphere
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       zfac=(sigma(jlev)/sigma(jlep))**akap
!$TAF STORE ilift = innern, REC = jlev + NLEM * nstep2
!$TAF LOOP = parallel
       do jhor=1,NHOR
        if(ilift(jhor) == -1                                            &
!    &    .and. zi(jhor) > 0.                   ! use for original kuo
     &    .and. zdqtot(jhor,jlep) > 0.) then

!
!     adiabatic uplift of parcel from below
!

         ztnew=zte(jhor,jlep)*zfac
         zqnew=zqe(jhor,jlep)

!
!     saturation humidity:
!

         zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))             &
     &        /(dp(jhor)*sigma(jlev))
         zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
! changed, FastOpt         zqsnl=zqsnl*zcor
         zqsnlh1=zqsnl*zcor ! changed, FastOpt
!
!     condensation if supersaturated
!

! changed, FastOpt         if(zqsnl < zqnew) then
         if(zqsnlh1 < zqnew) then ! changed, FastOpt
!
!     update constants (ice/water phase)
!

          if(ztnew < TMELT) then
           zlcp=ALS/(acpd*(1.+ADV*zqnew))
          else
           zlcp=ALV/(acpd*(1.+ADV*zqnew))
          endif

!
!     temperature correction term
!

! changed, FastOpt          zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
          zdqsdt=ra2*(TMELT-ra4)*zqsnlh1*zcor/((ztnew-ra4)*(ztnew-ra4)) ! changed, FastOpt


!
!     update t and q
!

! changed, FastOpt          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          zrain=(zqnew-zqsnlh1)/(1.+zlcp*zdqsdt) ! changed, FastOpt
! changed, FastOpt          ztnew=ztnew+zlcp*zrain
          ztnewh1=ztnew+zlcp*zrain ! changed, FastOpt
! changed, FastOpt          zqnew=zqnew-zrain
          zqnewh1=zqnew-zrain ! changed, FastOpt
!
!     second iteration:
!

! changed, FastOpt          zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))            &
! changed, FastOpt     &         /(dp(jhor)*sigma(jlev))
          zqsnlh2=rdbrv*ra1*exp(ra2*(ztnewh1-TMELT)/(ztnewh1-ra4))            &! changed, FastOpt
     &         /(dp(jhor)*sigma(jlev))                                   ! changed, FastOpt
! changed, FastOpt          zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
          zcorh1=1./(1.-(1./rdbrv-1.)*zqsnlh2) ! changed, FastOpt
! changed, FastOpt          zqsnl=zqsnl*zcor
          zqsnlh3=zqsnlh2*zcorh1 ! changed, FastOpt
! changed, FastOpt          zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
          zdqsdt=ra2*(TMELT-ra4)*zqsnlh3*zcorh1/((ztnewh1-ra4)*(ztnewh1-ra4)) ! changed, FastOpt
! changed, FastOpt          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          zrainh1=(zqnewh1-zqsnlh3)/(1.+zlcp*zdqsdt) ! changed, FastOpt
! changed, FastOpt          ztnew=ztnew+zlcp*zrain
          ztnewh2=ztnewh1+zlcp*zrainh1 ! changed, FastOpt
! changed, FastOpt          zqnew=zqnew-zrain
          zqnewh2=zqnewh1-zrainh1 ! changed, FastOpt

!
!     check for buyoncy
!

! changed, FastOpt          if(ztnew > zte(jhor,jlev)) then
          if(ztnewh2 > zte(jhor,jlev)) then ! changed, FastOpt
!
!     if convection: set lifting level and parcel q,t (zpq,ztp)
!                    calculate prelim. dq,dt (zdqdt,zdtdt)
!                    update surplus of total energy (zp)
!

           ilift(jhor)=jlep
           itop(jhor)=jlev
! changed, FastOpt           ztp(jhor)=ztnew
           ztp(jhor)=ztnewh2 ! changed, FastOpt
! changed, FastOpt           zqp(jhor)=zqnew
           zqp(jhor)=zqnewh2 ! changed, FastOpt
           zqsa=rdbrv*ra1*exp(ra2*(zte(jhor,jlev)-TMELT)                &
     &                           /(zte(jhor,jlev)-ra4))                 &
     &         /(dp(jhor)*sigma(jlev))
! changed, FastOpt           zqsa=zqsa/(1.-(1./rdbrv-1.)*zqsa)
           zqsah=zqsa/(1.-(1./rdbrv-1.)*zqsa) ! changed, FastOpt
! changed, FastOpt           zdqdt(jhor,jlev)=AMAX1(zqsa-zqe(jhor,jlev),0.)
           zdqdt(jhor,jlev)=AMAX1(zqsah-zqe(jhor,jlev),0.) ! changed, FastOpt
           zdtdt(jhor,jlev)=ztp(jhor)-zte(jhor,jlev)
           zpt(jhor)=zdtdt(jhor,jlev)*dsigma(jlev)/zalcpe(jhor,jlev)
           zpq(jhor)=zdqdt(jhor,jlev)*dsigma(jlev)
           if(kbetta > 0 ) then
! changed, FastOpt            zbetta(jhor)=AMIN1(1.,zqe(jhor,jlev)/zqsa)*dsigma(jlev)
            zbetta(jhor)=AMIN1(1.,zqe(jhor,jlev)/zqsah)*dsigma(jlev) ! changed, FastOpt
            zpbetta(jhor)=dsigma(jlev)
           endif

!
!          check if lifting level is supersaturated
!

! changed, FastOpt           ztnew=zte(jhor,jlep)
           ztnewh3=zte(jhor,jlep) ! changed, FastOpt
! changed, FastOpt           zqnew=zqe(jhor,jlep)
           zqnewh3=zqe(jhor,jlep) ! changed, FastOpt
! changed, FastOpt           zqsnl=rdbrv*ra1*exp(ra2*(zte(jhor,jlep)-TMELT)               &
! changed, FastOpt     &                            /(zte(jhor,jlep)-ra4))                &
! changed, FastOpt     &         /(dp(jhor)*sigma(jlep))
           zqsnlh4=rdbrv*ra1*exp(ra2*(zte(jhor,jlep)-TMELT)             & ! changed, FastOpt
     &                            /(zte(jhor,jlep)-ra4))                & ! changed, FastOpt
     &         /(dp(jhor)*sigma(jlep))                                    ! changed, FastOpt
! changed, FastOpt           zqsnl=zqsnl/(1.-(1./rdbrv-1.)*zqsnl)
           zqsnlh5=zqsnlh4/(1.-(1./rdbrv-1.)*zqsnlh4) ! changed, FastOpt
! changed, FastOpt           if(zqe(jhor,jlep) > zqsnl) then
           if(zqe(jhor,jlep) > zqsnlh5) then ! changed, FastOpt
!
!     update constants (ice/water phase)
!

! changed, FastOpt           if(ztnew < TMELT) then
           if(ztnewh3 < TMELT) then ! changed, FastOpt
! changed, FastOpt            zlcp=ALS/(acpd*(1.+ADV*zqnew))
            zlcp=ALS/(acpd*(1.+ADV*zqnewh3)) ! changed, FastOpt
           else
! changed, FastOpt            zlcp=ALV/(acpd*(1.+ADV*zqnew))
            zlcp=ALV/(acpd*(1.+ADV*zqnewh3)) ! changed, FastOpt
           endif

!
!     temperature correction term
!

! changed, FastOpt           zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
           zdqsdt=ra2*(TMELT-ra4)*zqsnlh5*zcorh1/((ztnewh3-ra4)*(ztnewh3-ra4)) ! changed, FastOpt

!
!     update t and q
!

! changed, FastOpt           zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
           zrainh2=(zqnewh3-zqsnlh5)/(1.+zlcp*zdqsdt) ! changed, FastOpt
! changed, FastOpt           ztnew=ztnew+zlcp*zrain
           ztnewh4=ztnewh3+zlcp*zrainh2 ! changed, FastOpt
! changed, FastOpt           zqnew=zqnew-zrain
           zqnewh4=zqnewh3-zrainh2 ! changed, FastOpt

!
!     second iteration:
!

! changed, FastOpt           zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))            &
! changed, FastOpt     &          /(dp(jhor)*sigma(jlep))
           zqsnlh6=rdbrv*ra1*exp(ra2*(ztnewh4-TMELT)/(ztnewh4-ra4))            & ! changed, FastOpt
     &          /(dp(jhor)*sigma(jlep))                                    ! changed, FastOpt
! changed, FastOpt           zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
           zcorh2=1./(1.-(1./rdbrv-1.)*zqsnlh6) ! changed, FastOpt
! changed, FastOpt           zqsnl=zqsnl*zcor
           zqsnlh7=zqsnlh6*zcorh2 ! changed, FastOpt
! changed, FastOpt           zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
           zdqsdt=ra2*(TMELT-ra4)*zqsnlh7*zcorh2/((ztnewh4-ra4)*(ztnewh4-ra4)) ! changed, FastOpt
! changed, FastOpt           zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
           zrainh3=(zqnewh4-zqsnlh7)/(1.+zlcp*zdqsdt) ! changed, FastOpt
! changed, FastOpt           ztnew=ztnew+zlcp*zrain
           ztnewh5=ztnewh4+zlcp*zrainh3 ! changed, FastOpt
! changed, FastOpt           zqnew=zqnew-zrain
           zqnewh5=zqnewh4-zrainh3 ! changed, FastOpt
! changed, FastOpt           zdqdt(jhor,jlep)=AMAX1(zqnew-zqe(jhor,jlep),0.)
           zdqdt(jhor,jlep)=AMAX1(zqnewh5-zqe(jhor,jlep),0.) ! changed, FastOpt
! changed, FastOpt           zdtdt(jhor,jlep)=ztnew-zte(jhor,jlep)
           zdtdt(jhor,jlep)=ztnewh5-zte(jhor,jlep) ! changed, FastOpt
           zpt(jhor)=zdtdt(jhor,jlep)*dsigma(jlep)/zalcpe(jhor,jlep)    &
     &              +zpt(jhor)
           zpq(jhor)=zdqdt(jhor,jlep)*dsigma(jlep)                      &
     &              +zpq(jhor)
           if(kbetta > 0 ) then
! changed, FastOpt            zbetta(jhor)=AMIN1(1.,zqe(jhor,jlep)/zqsnl)*dsigma(jlep)    &
! changed, FastOpt     &                  +zbetta(jhor)
            zbetta(jhor)=AMIN1(1.,zqe(jhor,jlep)/zqsnlh7)*dsigma(jlep)  & ! changed, FastOpt
     &                  +zbetta(jhor)                                     ! changed, FastOpt
            zpbetta(jhor)=dsigma(jlep)+zpbetta(jhor)
           endif
           endif

!
!     moisture accession (delete for original kuo)
!

           zi(jhor)=zdqtot(jhor,jlev)*dsigma(jlev)                      &
     &             +zdqtot(jhor,jlep)*dsigma(jlep)
          endif
         endif
        endif
       enddo
      enddo
!$TAF STORE ilift = inner1, REC = 1 + nstep2

!
!*    3) search for top level
!
!$TAF STORE zte = inner1, REC = 1 + nstep2
      do jlev=NLEM,1,-1
       jlep=jlev+1
       zfac=(sigma(jlev)/sigma(jlep))**akap
!$TAF STORE itop,zqp,ztp = innern, REC = jlev + nstep2 * NLEM
!$TAF LOOP = parallel
       do jhor=1,NHOR
          itoph=itop(jhor) ! added, FastOpt
          ztph=ztp(jhor) ! added, FastOpt
        if(jlep == itop(jhor)) then

!
!     3a) parcel uplift
!

         ztnew=ztp(jhor)*zfac
         zqnew=zqp(jhor)

!
!     saturation humidity
!

         zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))             &
     &        /(dp(jhor)*sigma(jlev))
         zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
! changed, FastOpt         zqsnl=zqsnl*zcor
         zqsnlh1=zqsnl*zcor ! changed, FastOpt

!
!     3b) more condensation
!
! changed, FastOpt         if(zqsnl < zqnew) then
         if(zqsnlh1 < zqnew) then ! changed, FastOpt
!
!     update constants (ice/water phase)
!

          if(ztnew < TMELT) then
           zlcp=ALS/(acpd*(1.+ADV*zqnew))
          else
           zlcp=ALV/(acpd*(1.+ADV*zqnew))
          endif

!
!     temperature correction term
!

! changed, FastOpt          zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
          zdqsdt=ra2*(TMELT-ra4)*zqsnlh1*zcor/((ztnew-ra4)*(ztnew-ra4)) ! changed, FastOpt

!
!     update t,q
!

! changed, FastOpt          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          zrain=(zqnew-zqsnlh1)/(1.+zlcp*zdqsdt) ! changed, FastOpt
! changed, FastOpt          ztnew=ztnew+zlcp*zrain
          ztnewh1 = ztnew + zlcp * zrain ! changed, FastOpt
! changed, FastOpt          zqnew=zqnew-zrain
          zqnewh1=zqnew-zrain ! changed, FastOpt

!
!     second iteration:
!
! changed, FastOpt          zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))            &
! changed, FastOpt     &         /(dp(jhor)*sigma(jlev))
          zqsnlh2=rdbrv*ra1*exp(ra2*(ztnewh1-TMELT)/(ztnewh1-ra4))        &! changed, FastOpt
     &         /(dp(jhor)*sigma(jlev))                                   ! changed, FastOpt
! changed, FastOpt          zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
          zcorh=1./(1.-(1./rdbrv-1.)*zqsnlh2) ! changed, FastOpt
! changed, FastOpt          zqsnl=zqsnl*zcor
          zqsnlh3=zqsnlh2*zcorh ! changed, FastOpt
! changed, FastOpt          zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
          zdqsdt=ra2*(TMELT-ra4)*zqsnlh3*zcorh/((ztnewh1-ra4)*(ztnewh1-ra4)) ! changed, FastOpt
! changed, FastOpt          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          zrainh=(zqnewh1-zqsnlh3)/(1.+zlcp*zdqsdt) ! changed, FastOpt
! changed, FastOpt          ztnew=ztnew+zlcp*zrain
          ztnewh2=ztnewh1+zlcp*zrainh ! changed, FastOpt
! changed, FastOpt          zqnew=zqnew-zrain
          zqnewh2=zqnewh1-zrainh ! changed, FastOpt

!
!     3c) check for buyoncy
!
! changed, FastOpt          if(ztnew > zte(jhor,jlev)) then
          if(ztnewh2 > zte(jhor,jlev)) then ! changed, FastOpt
!
!     3d) if further convection: update parcel q,t (zpq,ztp)
!                                calculate prelim. dq,dt (zdqdt,zdtdt)
!                                update surplus of total energy (zp)
!

! changed, FastOpt           itop(jhor)=jlev
           itoph = jlev ! changed, FastOpt
! changed, FastOpt           ztp(jhor)=ztnew
           ztph=ztnewh2 ! changed, FastOpt
! changed, FastOpt           zqp(jhor)=zqnew
           zqp(jhor)=zqnewh2 ! changed, FastOpt
           zqsa=rdbrv*ra1*exp(ra2*(zte(jhor,jlev)-TMELT)                &
     &                            /(zte(jhor,jlev)-ra4))                &
     &         /(dp(jhor)*sigma(jlev))
! changed, FastOpt           zqsa=zqsa/(1.-(1./rdbrv-1.)*zqsa)
           zqsah=zqsa/(1.-(1./rdbrv-1.)*zqsa) ! changed, FastOpt
! changed, FastOpt           zdqdt(jhor,jlev)=AMAX1(zqsa-zqe(jhor,jlev),0.)
           zdqdt(jhor,jlev)=AMAX1(zqsah-zqe(jhor,jlev),0.) ! changed, FastOpt
!           zdqdt(jhor,jlev)=sqrt(zdqdt(jhor,jlev)**2+1.e-3**2)-1.e-3 ! added, FastOpt
! changed, FastOpt           zdtdt(jhor,jlev)=ztp(jhor)-zte(jhor,jlev)
           zdtdt(jhor,jlev)=ztph-zte(jhor,jlev) ! FastOpt
           zpt(jhor)=zdtdt(jhor,jlev)*dsigma(jlev)/zalcpe(jhor,jlev)    &
     &              +zpt(jhor)
           zpq(jhor)=zdqdt(jhor,jlev)*dsigma(jlev)                      &
     &              +zpq(jhor)
           if(kbetta > 0 ) then
            zbetta(jhor)=AMIN1(1.,zqe(jhor,jlev)/zqsah)*dsigma(jlev)     &
     &                  +zbetta(jhor)
            zpbetta(jhor)=zpbetta(jhor)+dsigma(jlev)
           endif

!
!     moisture accession (delete for original kuo)
!

           zi(jhor)=zdqtot(jhor,jlev)*dsigma(jlev)                      &
     &             +zi(jhor)
          endif
         endif
        endif
        ztp(jhor) = ztph ! added, FastOpt
        itop(jhor) = itoph ! added, FastOpt
       enddo
      enddo
!$TAF STORE itop = inner1, REC = 1 + nstep2

!
!*    4) calc. moisture distribution factors
!
!$taf STORE zdqdt,zdtdt,zpbetta,zi,zpt,zpq = inner1, REC = 1+nstep2
      zi(:)=AMAX1(0.,zi(:))

      if(kbetta==0) then

!
!     kuo without betta
!

       where(zpt(:)+zpq(:) > 0.)
        zat(:)=zi(:)/(zpt(:)+zpq(:))
        zaq(:)=zat(:)
       end where
      else

!
!     kuo with betta
!
!$taf store zbetta = inner1, REC = 1+nstep2
       where(zpbetta(:) > 0.) zbetta(:)=(1.-zbetta(:)/zpbetta(:))**3
       where(zpt(:) > 0.) zat(:)=zi(:)*(1.-zbetta(:))/zpt(:)
       where(zpq(:) > 0.) zaq(:)=zi(:)*zbetta(:)/zpq(:)

      endif
!$taf STORE zat,zaq = inner1, REC = 1+nstep2
!
!*    5) update dqdt and dtdt, compute precip.
!

      do jlev=1,NLEV

!
!     5a) in cloud: add tendencies due to convection and sub. advection
!         lifting level: sub. advection. (dtdt != 0 only if level supersat)
!
!$taf store ilift,itop = innern, REC = jlev+nstep2*NLEV
       where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.)
         dqdt(:,jlev)=zaq(:)*zdqdt(:,jlev)/zdelt                        &
     &               -dqt(:,jlev)
         dtdt(:,jlev)=zat(:)*zdtdt(:,jlev)/zdelt                        &
     &               +dtdt(:,jlev)

!
!     5b) calc. precipitation
!

         dprc(:)=dprc(:)                                                &
     &          +zat(:)*zdtdt(:,jlev)*dsigma(jlev)/zalcpe(:,jlev)

       endwhere

!
!     5c) original kuo:
!         rest of column: subtrac all q-tendencies (i.e. reset dqdt)
!
!      where(zi(:) > 0 .and. (jlev >= ilift(:) .or. jlev < itop(:)))
!       dqdt(:,jlev)=-1.*dqt(:,jlev)
!      endwhere

!
!     5d) flag and counter for clouds:
!

       where(zat(:)*zdtdt(:,jlev) > 0.)
        icclev(:,jlev)=1
        icctot(:)=icctot(:)+1
       endwhere
       if(sigma(jlev) < 0.4) then
        where(icclev(:,jlev)==1 .and. itop(:)==jlev)
         icclev(:,jlev)=2
        endwhere
       endif

      enddo
!$taf STORE dprc = inner1, REC = 1 + nstep2

!
!*    6) final precip. in m/s and snow fall
!

      zfac1=1./ga/zdelt/1000.
      dprc(:)=dprc(:)*dp(:)*zfac1
      where(dt(:,NLEV) <= TMELT) dprs(:)=dprs(:)+dprc(:)

!
!*    7) momentum mixing
!

      if(nmoment == 1) then

!
!     7a) update u and v
!

       zue(:,:)=du(:,1:NLEV)+dudt(:,1:NLEV)*zdelt
       zve(:,:)=dv(:,1:NLEV)+dvdt(:,1:NLEV)*zdelt

!
!     7a) make distribution factor
!

       zau(:)=0.
       where(zpt(:)+zpq(:) > 0.)
        zau(:)=AMIN1(1.,zi(:)/(zpt(:)+zpq(:)))
       end where

!
!     7b) reference u and v
!

       zup(:)=0.
       zvp(:)=0.
       zsums(:)=0.
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.)
         zup(:)=zup(:)+zue(:,jlev)*dsigma(jlev)
         zvp(:)=zvp(:)+zve(:,jlev)*dsigma(jlev)
         zsums(:)=zsums(:)+dsigma(jlev)
        endwhere
       enddo
!$TAF STORE zsums = inner1, REC = 1 + nstep2
       where(zsums(:) > 0.)
        zup(:)=zup(:)/zsums(:)
        zvp(:)=zvp(:)/zsums(:)
       endwhere
!$TAF STORE zup, zvp = inner1, REC = 1 + nstep2
!
!     7c) make and add tendencies
!
!$taf loop = parallel
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.)
          zdudt(:,jlev)=zau(:)*(zup(:)-zue(:,jlev))
          zdvdt(:,jlev)=zau(:)*(zvp(:)-zve(:,jlev))
          dudt(:,jlev)=dudt(:,jlev)+zdudt(:,jlev)/zdelt
          dvdt(:,jlev)=dvdt(:,jlev)+zdvdt(:,jlev)/zdelt
        endwhere
       enddo

      endif
!
!*    dbug diagnostic for nprint==2 (see pumamod)
!

      if(nprint==2) then
       allocate(zprf1(NLON*NLAT,NLEV))
       allocate(zprf2(NLON*NLAT,NLEV))
       allocate(zprf3(NLON*NLAT,NLEV))
       allocate(zprf4(NLON*NLAT,NLEV))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       allocate(zprf8(NLON*NLAT))
       allocate(zprf9(NLON*NLAT))
       allocate(zprf10(NLON*NLAT))
       allocate(zprf11(NHOR))
       allocate(zprf12(NHOR))
       zprf11(:)=real(ilift(:))
       zprf12(:)=real(itop(:))
       call mpgagp(zprf1,zdqdt,NLEV)
       call mpgagp(zprf2,zdtdt,NLEV)
       call mpgagp(zprf3,zdqtot,NLEV)
       call mpgagp(zprf4,dqt,NLEV)
       call mpgagp(zprf5,zprf11,1)
       call mpgagp(zprf6,zprf12,1)
       call mpgagp(zprf7,dprc,1)
       call mpgagp(zprf8,zaq,1)
       call mpgagp(zprf9,zat,1)
       call mpgagp(zprf10,zi,1)
       deallocate(zprf11)
       deallocate(zprf12)
       if(mypid==NROOT) then
        print*,'In Kuo:'
        print*,'ilift= ',zprf5(nprhor),' itop= ',zprf6(nprhor)
        do jlev=1,NLEV
         print*,'L= ',jlev,' dqtot= ',zprf3(nprhor,jlev)                &
     &                    ,' dqt= ',zprf4(nprhor,jlev)
         if(zprf5(nprhor) >= jlev .and. zprf6(nprhor) <= jlev           &
     &     .and. zprf10(nprhor) > 0.) then
          zzdt=zprf9(nprhor)*zprf2(nprhor,jlev)/zdelt
          zzdq=zprf8(nprhor)*zprf1(nprhor,jlev)/zdelt
          print*,'L= ',jlev,' dtdt= ',zzdt,' dqdt= ',zzdq
         endif
        enddo
        print*,'prc (mm/day)= ',zprf7(nprhor)*1000.*ntspd*zdelt*0.5
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
       deallocate(zprf10)
      endif

!
!     franks dbug
!

      if(ndiaggp==1) then
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.)
         dgp3d(:,jlev,8)=zat(:)*zdtdt(:,jlev)/zdelt
!!         dgp3d(:,jlev,15)=(zaq(:)*zdqdt(:,jlev)/zdelt-dqt(:,jlev))      &
         dgp3d(:,jlev,15)=(zaq(:)*zdqdt(:,jlev)-zdqtot(:,jlev))/zdelt   &
     &                   *dp(:)*dsigma(jlev)/ga/1000.
        elsewhere
         dgp3d(:,jlev,8)=0.
         dgp3d(:,jlev,15)=0.
        endwhere
       enddo
      end if

      return
      end subroutine kuo

!     ===================
!     SUBROUTINE MKCLOUDS
!     ===================

      subroutine mkclouds
      use rainmod
!
!     compute cloud properties (cover and liquid water content)
!

      parameter(zclmax=1.)
      parameter(zcca=0.245,zccb=0.125)
      parameter(zccmax=0.8,zccmin=0.05)

      real zq(NHOR),zt(NHOR),zqsat(NHOR),zrh(NHOR)
      real zrcrit(NHOR),zccconv(NHOR),zcctot(NHOR)
      real zwfac(NHOR),zcc(NHOR)
      real zzf(NHOR,NLEV),zzh(NHOR),zdh(NHOR)
      real dcch(NHOR) ! FastOpt - aux.variable to make adjoint tape handling easier

!
!*    preset
!

      zdelt=delt2/ww
      dql(:,:)=0.
      dcc(:,1:NLEV)=0.
      zcctot(:) = 0. ! FastOpt
      zccconv(:) = 0. ! FastOpt
!
!*    convective cloud cover for random overlab:
!

      zrfac = solar_day * 1000.0 ! convert m/s into mm/day
      where(dprc(:) > 0.)
       zcctot(:)=zcca+zccb*log(dprc(:)*zrfac)
       zcctot(:)=AMIN1(zccmax,AMAX1(zccmin,zcctot(:)))
       zccconv(:)=1.-(1.-zcctot(:))**(1./real(icctot(:)))
      end where

!
!*    cloud cover for the levels
!
!$taf loop = parallel
      do jlev=1,NLEV
       zt(:)=dt(:,jlev)+dtdt(:,jlev)*zdelt
       zq(:)=dq(:,jlev)+dqdt(:,jlev)*zdelt
       zqsat(:)=rdbrv*ra1*exp(ra2*(zt(:)-TMELT)/(zt(:)-ra4))            &
     &         /(dp(:)*sigma(jlev))
       zqsat(:)=zqsat(:)/(1.-(1./rdbrv-1.)*zqsat(:))
       zrh(:)=zq(:)/zqsat(:)
!
!     convective clouds
!
!
!     a) cu/cb
!
       where(icclev(:,jlev) > 0)
        dcc(:,jlev)=zccconv(:)
       endwhere
!  !$taf store dcc(:,jlev) = innern, rec = jlev + NLEV*nstep2
!
!     b) ci
!
       where(icclev(:,jlev)==2 .and. zcctot(:) >= 0.3)
        dcc(:,jlev)=AMIN1(dcc(:,jlev)+2.*(zcctot(:)-0.3),zclmax)
       endwhere
!   !$taf store dcc(:,jlev) = innern, rec = jlev + NLEV*nstep2
!
!     reduce relative humidity
!
        zrh(:)=zrh(:)*(1.-dcc(:,jlev))
!
!     c) non convective clouds
!rcritg
       zrcrit(:)=rcritm(:,jlev)
       zwfac(:)=1.
       if(sigma(jlev) > 0.8 .and. clwfac > 0.)  then
        zwfac(:)=clwfac*AMAX1(0.,clwcrit1-dw(:,jlev))
        where(clwcrit2 > dw(:,jlev)) zwfac(:)=1.
       endif
       zcc(:)=zwfac(:)*AMAX1(0.,(zrh(:)-zrcrit(:))/(1.-zrcrit(:)))**2
       dcc(:,jlev)=AMIN1(dcc(:,jlev)+zcc(:),zclmax)

      enddo

!
!*    precipitable water (and z which is later used to compute dql)
!

      dqvi(:)=0.
      dqviraw(:) = 0. ! test, FastOpt
      zzh(:)=0.
      do jlev=NLEV,2,-1
       dqvi(:)=dq(:,jlev)*dsigma(jlev)+dqvi(:)
       dqviraw(:)=dqraw(:,jlev)*dsigma(jlev)+dqviraw(:) ! test, FastOpt
       zdh(:)=-dt(:,jlev)*gascon/ga*ALOG(sigmah(jlev-1)/sigmah(jlev))
       zzf(:,jlev)=zzh(:)+zdh(:)*0.5
       zzh(:)=zzh(:)+zdh(:)
      enddo
      dqvi(:)=dq(:,1)*dsigma(1)+dqvi(:)
      dqviraw(:)=dqraw(:,1)*dsigma(1)+dqviraw(:) ! test, FastOpt
      zdh(:)=-dt(:,1)*gascon/ga*ALOG(sigma(1)/sigmah(1))*0.5
      zzf(:,1)=zzh(:)+zdh(:)

!
!     cloud liquid water (see CCM3 description: Kiehl et al. 1996 pp49+50)
!     and total cloud cover for random overlap
!

!      dcc(:,NLEP)=1.
      dcch(:) = 1. ! FastOpt
      dqvi(:)=dqvi(:)*dp(:)/ga
      zzh(:)=700.*ALOG(1.+dqvi(:))
!$taf store dcc = inner1, rec = 1 + nstep2
      do jlev=1,NLEV
!$taf store dcch(:) = innern, rec = jlev + nstep2*NLEV
       where(zzh(:) > 0. .and. dcc(:,jlev) > 0.)
        dql(:,jlev)=0.00021*EXP(-zzf(:,jlev)/zzh(:))*gascon*dt(:,jlev)  &
     &             /(sigma(jlev)*dp(:))
        dql(:,jlev)=MAX(dql(:,jlev),1.E-9)
       endwhere
!       dcc(:,NLEP)=dcc(:,NLEP)*(1.-dcc(:,jlev))
       dcch(:)=dcch(:)*(1.-dcc(:,jlev))
      enddo
!      dcc(:,NLEP)=1.-dcc(:,NLEP)
      dcc(:,NLEP)=1.-dcch(:)
!
      dqvi(:) = dqviraw(:)*dp(:)/ga ! FastOpt: for output
      return
      end subroutine mkclouds

!     ==============
!     SUBROUTINE DCA
!     ==============

      subroutine mkdca
      use rainmod
!
!     calculate t and q tendencies due to dry convection
!
      integer itop(NHOR),ibase(NHOR),icon(NHOR)
      real zdtdt(NHOR,NLEV)
      real zdqdt(NHOR,NLEV)
      real zt(NHOR,NLEV),zth(NHOR,NLEV)
      real zq(NHOR,NLEV),zqn(NHOR,NLEV)
      real ztht(NHOR)
      real zqt(NHOR)
!FastOpt      real zsum1(NHOR),zsum2(NHOR)
      real :: zsum1(NHOR) = 0.! added, FastOpt
      real :: zsum2(NHOR) = 0.! added, FastOpt
!FastOpt      real zsumq1(NHOR),zsumq2(NHOR)
      real :: zsumq1(NHOR) = 0.! added, FastOpt
      real :: zsumq2(NHOR) = 0.! added, FastOpt
      real zskap(NLEV)
!
      zdelt=delt2/ww
!
      do jlev = 1,NLEV
       zskap(jlev)=sigma(jlev)**akap
       zt(:,jlev)=dt(:,jlev)+dtdt(:,jlev)*zdelt
       zth(:,jlev)=zt(:,jlev)/zskap(jlev)
       zq(:,jlev)=dq(:,jlev)+dqdt(:,jlev)*zdelt
       zqn(:,jlev)=zq(:,jlev)
      enddo
!$TAF STORE zskap,zth = inner1, REC = 1 + nstep2
!FastOpt      jiter=0
!FastOpt  1000 continue
       icon(:)=0
       ibase(:)=0
       itop(:)=NLEV
       ztht(:)=0.
       zqt(:)=0.
       do jlev=1,NLEM
        jlep=jlev+1
        where(zth(:,jlev) < zth(:,jlep))
         ibase(:)=jlep
         itop(:)=jlev-1
         zsum1(:)=zskap(jlep)*dsigma(jlep)                              &
     &              +zskap(jlev)*dsigma(jlev)
         zsum2(:)=zskap(jlep)*dsigma(jlep)*zth(:,jlep)                  &
     &           +zskap(jlev)*dsigma(jlev)*zth(:,jlev)
         zsumq1(:)=dsigma(jlep)+dsigma(jlev)
         zsumq2(:)=dsigma(jlep)*zqn(:,jlep)+dsigma(jlev)*zqn(:,jlev)
         ztht(:)=zsum2(:)/zsum1(:)
         zqt(:)=zsumq2(:)/zsumq1(:)
         icon(:)=1
        endwhere
       enddo
       do jiter=0,2*NLEV !added, FastOpt
!$TAF STORE icon = inner2n, REC = jiter+1+(2*NLEV+1)*nstep2
          if(SUM(icon(:)) > 0) then
!$TAF STORE ibase,itop,zsum1,zsum2,zsumq1,zth,ztht = inner2n, REC = jiter+1+(2*NLEV+1)*nstep2
             do jlev=NLEM-1,1,-1
                where(itop(:) == jlev .and. zth(:,jlev) < ztht(:))
                   itop(:)=jlev-1
                   zsum1(:)=zsum1(:)+zskap(jlev)*dsigma(jlev)
                   zsum2(:)=zsum2(:)+zskap(jlev)*dsigma(jlev)*zth(:,jlev)
                   zsumq1(:)=zsumq1(:)+dsigma(jlev)
                   zsumq2(:)=zsumq2(:)+dsigma(jlev)*zqn(:,jlev)
                   ztht(:)=zsum2(:)/zsum1(:)
                   zqt(:)=zsumq2(:)/zsumq1(:)
                endwhere
             enddo
!$TAF STORE itop,zsum,zsumq1,ztht = inner2n, REC = jiter+1+(2*NLEV+1)*nstep2
             do jlev=1,NLEV
                where(jlev > itop(:) .and. jlev <= ibase(:))
                   zth(:,jlev)=ztht(:)
                   zqn(:,jlev)=zqt(:)
                endwhere
             enddo
!$TAF STORE zth = inner2n, REC = jiter+1+(2*NLEV+1)*nstep2
!FastOpt        jiter=jiter+1
!FastOpt        if(jiter > 2*NLEV) goto 1001
!FastOpt        goto 1000
             icon(:)=0                                                              ! added, FastOpt
             ibase(:)=0                                                             ! added, FastOpt
             itop(:)=NLEV                                                           ! added, FastOpt
             ztht(:)=0.                                                             ! added, FastOpt
             zqt(:)=0.                                                              ! added, FastOpt
!$TAF STORE ibase,icon,itop,zsum1,zsum2,zsumq1,zth,ztht = inner2n, REC = jiter+1+(2*NLEV+1)*nstep2
             do jlev=1,NLEM                                                         ! added, FastOpt
                jlep=jlev+1                                                         ! added, FastOpt
                where(zth(:,jlev) < zth(:,jlep))                                    ! added, FastOpt
                   ibase(:)=jlep                                                    ! added, FastOpt
                   itop(:)=jlev-1                                                   ! added, FastOpt
                   zsum1(:)=zskap(jlep)*dsigma(jlep)                              & ! added, FastOpt
                        &              +zskap(jlev)*dsigma(jlev)                    ! added, FastOpt
                   zsum2(:)=zskap(jlep)*dsigma(jlep)*zth(:,jlep)                  & ! added, FastOpt
                        &           +zskap(jlev)*dsigma(jlev)*zth(:,jlev)           ! added, FastOpt
                   zsumq1(:)=dsigma(jlep)+dsigma(jlev)                              ! added, FastOpt
                   zsumq2(:)=dsigma(jlep)*zqn(:,jlep)+dsigma(jlev)*zqn(:,jlev)      ! added, FastOpt
                   ztht(:)=zsum2(:)/zsum1(:)                                        ! added, FastOpt
                   zqt(:)=zsumq2(:)/zsumq1(:)                                       ! added, FastOpt
                   icon(:)=1                                                        ! added, FastOpt
                endwhere                                                            ! added, FastOpt
             enddo                                                                  ! added, FastOpt
          endif
       end do !added, FastOpt

!$taf STORE icon = inner1, REC = 1+nstep2
!
       if(SUM(icon(:)) ==  0) then ! added, Fastopt
!
          do jlev=1,NLEV
             zdtdt(:,jlev)=(zth(:,jlev)*zskap(jlev)-zt(:,jlev))/zdelt
             zdqdt(:,jlev)=(zqn(:,jlev)-zq(:,jlev))/zdelt
             dtdt(:,jlev)=dtdt(:,jlev)+zdtdt(:,jlev)
             dqdt(:,jlev)=dqdt(:,jlev)+zdqdt(:,jlev)
          enddo
!
!     franks dbug
!
          if(ndiaggp==1) then
             dgp3d(:,1:NLEV,9)=zdtdt(:,1:NLEV)
             do jlev=1,NLEV
                dgp3d(:,jlev,16)=zdqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
             enddo
          endif
       end if ! added, FastOpt
!
      return
!FastOpt  1001 continue
      end subroutine mkdca
