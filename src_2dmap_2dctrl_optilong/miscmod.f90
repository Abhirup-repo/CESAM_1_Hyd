      module miscmod
      use pumamod
!
!     version identifier
!
      character(len=80) :: version = '03.03.2003 by Larry'
!
!     namelist parameters:
!
      integer :: nfixer = 1 ! switch for negative humidity fix (1/0 : on/off)
      integer :: nudge  = 0 ! switch for t-nudging (1/0 : on/off, 2: flxcorr)
      real    :: tnudget = 10. ! timescale for t-nudging (days)
!
!     global scalars
!
      real :: time4mis = 0.  ! CPU time needed for misc routines
      real :: time4fix = 0.  ! CPU time needed for fixer
!
!     global arrays
!
      real :: zgw(NHOR)          ! gaussian weights
      real :: dtnudge(NHOR,0:13) ! climatological temperature (annual cycle)
      real :: dfnudge(NHOR,0:13) ! climatological fluxes correction (annual cycle)
!
      end module miscmod

!     ==================
!     SUBROUTINE MISCINI
!     ==================

      subroutine miscini
      use miscmod
!
!     initialization
!
      namelist/miscpar/nfixer,nudge,tnudget
!
      if (mypid == NROOT) then
         read(11,miscpar)
         write (*,'(/," *****************************************")')
         write (*,'(" * MISCMOD ",a29," *")') trim(version)
         write (*,'(" *****************************************")')
         write (*,'(" * Namelist MISCPAR from <puma_namelist> *")')
         write (*,'(" *****************************************")')
         write(*,miscpar)
      endif

      call mpbci(nfixer)
      call mpbci(nudge)
      call mpbcr(tnudget)

      if (nudge == 1) then
         call mpsurfyear('dtnudge',dtnudge,NHOR)
      endif

      if (nudge == 2) then
         call mpsurfyear('dfnudge',dfnudge,NHOR)
      endif
!
!     gaussian weights
!
      do jlat=1,NLPP
       jhor1=(jlat-1)*NLON+1
       jhor2=jlat*NLON
       zgw(jhor1:jhor2)=gw(jlat)*csq(jlat)
      enddo
!
!     nudging time in seconds
!
      tnudget=TWOPI*tnudget/ww
!
      return
      end subroutine miscini

!     ===================
!     SUBROUTINE MISCSTEP
!     ===================

      subroutine miscstep
      use miscmod

      if(ntime == 1) then
       call mksecond(zsec,0.)
       call mksecond(zsec1,0.)
      endif
! test, FastOpt      call fixer
      call fixer
!      dqraw(:,1:NLEP) = dq(:,1:NLEP) ! test, FastOpt
!      do jlev=1,NLEP                                      ! test, FastOpt
!!         dq(:,jlev) = fomaxf(NHOR,dqraw(:,jlev),0.,1.e-5) ! test, FastOpt
!         dq(:,jlev) = max(dqraw(:,jlev),0.)               ! test, FastOpt
!      end do                                              ! test, FastOpt
      if(ntime == 1) then
       call mksecond(zsec1,zsec1)
       time4fix=time4fix+zsec1
      endif
!
      if(nudge == 1) then
       call mknudge
      elseif(nudge == 2) then
       call mkflcor
      endif
!
      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4mis=time4mis+zsec
      endif

      return
      end subroutine miscstep

!     ===================
!     SUBROUTINE MISCSTOP
!     ===================

      subroutine miscstop
      use miscmod

      if(mypid == NROOT .and. ntime == 1) then
       print*,'*********************************************'
       print*,' CPU usage in MISCSTEP (ROOT process only):  '
       print*,'    All routines : ',time4mis,' s'
       print*,'    Fixer        : ',time4fix,' s'
       print*,'*********************************************'
      endif

      return
      end subroutine miscstop


      subroutine fixer
      use miscmod
!
      real zqn(NHOR,NLEV)    ! humidity in and out
      real zneg(NHOR)        ! moisture needed
      real zpos(NHOR)        ! moisture available
!
      real zsum(2)           ! work array for mpi (needs to be an array)
!
!*    franks dbug
!
      if(ndiaggp > 0) then
       dgp3d(:,1:NLEV,11)=dq(:,1:NLEV)
      endif
!
      zdelt=delt2/ww
!
      zqn(:,:)=dq(:,1:NLEV)
      zneg(:)=0.
      zpos(:)=0.
      do jhor=1,NHOR
       zneg(jhor)=-SUM(zqn(jhor,:)*dsigma(:),MASK=(zqn(jhor,:) < 0.))   &
     &           *dp(jhor)*zgw(jhor)
       zpos(jhor)=SUM(zqn(jhor,:)*dsigma(:),MASK=(zqn(jhor,:) > 0.))    &
     &           *dp(jhor)*zgw(jhor)
      enddo
!
!*    fix within grid point collumn
!
      if(SUM(zneg(:)) > 0.) call fix0d(zqn,zneg,zpos) ! problematic code for derivative - FastOpt
!
!*    fix within latitude
!
      if(SUM(zneg(:)) > 0.) call fix1d(zqn,zneg,zpos) ! problematic code for derivative - FastOpt
!
!*    fix globally
!
      zsum(1)=SUM(zneg(:))
      zsum(2)=SUM(zpos(:))
      call mpsumbcr(zsum,2)
      zzneg=zsum(1)
      zzpos=zsum(2)
      if(zzneg > 0.) then
       zfac=(zzpos-zzneg)/zzpos
       zqn(:,:)=AMAX1(zqn(:,:)*zfac,0.)
      endif
!
!*    tendencies
!
      do jlev=1,NLEV
       gqdt(:,jlev)=gqdt(:,jlev)                                        &
     &             +dp(:)*(zqn(:,jlev)-dq(:,jlev))/delt2/psurf
      enddo
!
!*    set new dq
!
!$TAF incompl dq
      dq(:,1:NLEV)=zqn(:,:)
!
!*    franks dbug
!
      if(ndiaggp > 0) then
       do jlev=1,NLEV
        dgp3d(:,jlev,11)=(dq(:,jlev)-dgp3d(:,jlev,11))/zdelt*dp(:)      &
     &                  *dsigma(jlev)/ga/1000.
       enddo
      endif
!
      return
      end subroutine fixer

!     ================
!     SUBROUTINE FIX0D
!     ================

      subroutine fix0d(pqn,pneg,ppos)
      use miscmod
!
      real pqn(NHOR,NLEV)
      real pneg(NHOR),ppos(NHOR)
!
      do jhor=1,NHOR
       if (pneg(jhor) > 0. .and. ppos(jhor) >= pneg(jhor)) then
        zfac=(ppos(jhor)-pneg(jhor))/ppos(jhor)
        pqn(jhor,:)=AMAX1(pqn(jhor,:)*zfac,0.)
        pneg(jhor)=0.
        ppos(jhor)=ppos(jhor)*zfac
       endif
      enddo
!
      return
      end subroutine fix0d

!     ================
!     SUBROUTINE FIX1D
!     ================

      subroutine fix1d(pqn,pneg,ppos)
      use miscmod
!
      real pqn(NHOR,NLEV)
      real pneg(NHOR),ppos(NHOR)
!
      do jlat=1,NLPP
       jhor1=(jlat-1)*NLON+1
       jhor2=jlat*NLON
       zneg=SUM(pneg(jhor1:jhor2))
       if(zneg > 0.) then
        zpos=SUM(ppos(jhor1:jhor2))
        if(zpos >= zneg) then
          zfac=(zpos-zneg)/zpos
          pqn(jhor1:jhor2,:)=AMAX1(pqn(jhor1:jhor2,:)*zfac,0.)
          pneg(jhor1:jhor2)=0.
          ppos(jhor1:jhor2)=ppos(jhor1:jhor2)*zfac
        endif
       endif
      enddo
!
      return
      end subroutine fix1d

!     ==================
!     SUBROUTINE MKNUDGE
!     ==================

      subroutine mknudge
      use miscmod
!
!     make nudging of upper level temperature
!
      real ztcl(NHOR)
      real zdtdt(NHOR)

!
!     get climatology
!

      call getnudge(ztcl)

!
!     tendencies due to nudging
!

      zdtdt(:)=(ztcl(:)-dt(:,1))/tnudget

!
!     add tendencies
!

      dtdt(:,1)=dtdt(:,1)+zdtdt(:)

!
!     franks diagnostic
!

      if(ndiaggp==1) then
       dgp3d(:,2:NLEV,2)=0.
       dgp3d(:,1,2)=zdtdt(:)
      endif
!

      return
      end

!     ==================
!     SUBROUTINE MKFLCOR
!     ==================

      subroutine mkflcor
      use miscmod
!
!     make flux correction for upper level temperature
!
      real zfcl(NHOR)
      real zdtdt(NHOR)

!
!     get climatology
!

      call getflcor(zfcl)

!
!     tendencies due to flux correction
!

      zdtdt(:)=zfcl(:)

!
!     add tendencies
!

      dtdt(:,1)=dtdt(:,1)+zdtdt(:)
!
!     franks diagnostic
!

      if(ndiaggp==1) then
       dgp3d(:,2:NLEV,2)=0.
       dgp3d(:,1,2)=zdtdt(:)
      endif
!
      return
      end

!     ===================
!     SUBROUTINE GETNUDGE
!     ===================

      subroutine getnudge(ptcl)
      use miscmod
!
      real ptcl(NHOR)
!
!     get surface temperature annual cycle
!
      jstep=nstep+1    !advance time step (ts = ts(t)+ dtsdt)
      call ntomin(jstep,nmin,nhour,nday,nmonth,nyear)
!
      if(nperpetual > 0) then
       nmonth=(nperpetual-1)/30.+1
       nday=nperpetual-nmonth*30
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
!*    interpolate
!
      zgw2=abs(nday-15)/30.
      zgw1=1.-zgw2
      ptcl(:)=zgw1*dtnudge(:,jm1)+zgw2*dtnudge(:,jm2)
!
      return
      end subroutine getnudge

!     ===================
!     SUBROUTINE GETFLCOR
!     ===================

      subroutine getflcor(pfcl)
      use miscmod
!
      real pfcl(NHOR)
!
!     get flux correction annual cycle
!
      jstep=nstep+1    !advance time step (ts = ts(t)+ dtsdt)
      call ntomin(jstep,nmin,nhour,nday,nmonth,nyear)
!
      if(nperpetual > 0) then
       nmonth=(nperpetual-1)/30.+1
       nday=nperpetual-nmonth*30
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
!*    interpolate
!
      zgw2=abs(nday-15)/30.
      zgw1=1.-zgw2
      pfcl(:)=zgw1*dfnudge(:,jm1)+zgw2*dfnudge(:,jm2)
!
      return
      end subroutine getflcor
