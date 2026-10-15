      module fluxmod
      use pumamod
!
!     version identifier
!

      character(len=80) :: version = '22.11.2002 by Larry'

!
!     global parameter
!

      parameter(vonkarman=0.4)    ! von karman const.

!
!     namelist parameter:
!

      integer :: nvdiff  = 1      ! switch for vertical diffusion
      integer :: nshfl   = 1      ! switch for sensible heat flux
      integer :: nevap   = 1      ! switch for evaporation/latend heat flux
      integer :: nstress = 1      ! switch for wind stress
      integer :: ntsa    = 2      ! flag for near surface temp calc.

      real :: vdiff_lamm  = 160.  ! const. used in vdiff (see parameterization)
      real :: vdiff_b     = 5.    !        "
      real :: vdiff_c     = 5.    !        "
      real :: vdiff_d     = 5.    !        "

!
!     arrays
!

      real :: dtransh(NHOR)       ! transfere coeficient for heat
      real :: dtransm(NHOR)       ! transfere coeficient for momentum

!
!     real scalars
!

      real :: time4fl = 0.        ! CPU time needed for fluxmod
      real :: time4sf = 0.        ! CPU time needed for surface fluxes
      real :: time4tr = 0.        ! CPU tiem needed for transfere coeff.
      real :: time4sh = 0.        ! CPU time needed for sensible heat flux
      real :: time4ev = 0.        ! CPU time needed for evaporation
      real :: time4st = 0.        ! CPU time needed for wind stress
      real :: time4vd = 0.        ! CPU time needed for vertical diffusion

      end module fluxmod

!     ==================
!     SUBROUTINE FLUXINI
!     ==================

      subroutine fluxini
      use fluxmod
!
      namelist/fluxpar/nvdiff,nshfl,nevap,nstress,ntsa                  &
     &                ,vdiff_lamm,vdiff_b,vdiff_c,vdiff_d
!
      if(mypid==NROOT) then
         read(11,fluxpar)
         write (*,'(/," *****************************************")')
         write (*,'(" * FLUXMOD ",a29," *")') trim(version)
         write (*,'(" *****************************************")')
         write (*,'(" * Namelist FLUXPAR from <puma_namelist> *")')
         write (*,'(" *****************************************")')
         write(*,fluxpar)
      endif

! guokun added to initialize parameter matirx 21.Nov
      alpha1g(1:NUGP)=1.
      beta1g(1:NUGP)=1.
      gamma1g(1:NUGP)=1.
      call mpscgp(alpha1g,alpha1m,1)
      call mpscgp(beta1g,beta1m,1)
      call mpscgp(gamma1g,gamma1m,1)

      call mpbci(nvdiff)
      call mpbci(nshfl)
      call mpbci(nevap)
      call mpbci(nstress)
      call mpbci(ntsa)
      call mpbcr(vdiff_lamm)
      call mpbcr(vdiff_b)
      call mpbcr(vdiff_c)
      call mpbcr(vdiff_d)

      return
      end subroutine fluxini

!     ===================
!     SUBROUTINE FLUXSTEP
!     ===================

      subroutine fluxstep
      use fluxmod

      if(ntime == 1) call mksecond(zsec,0.)
      call surflx
      if(ntime == 1) then
       call mksecond(zsec1,zsec)
       time4sf=time4sf+zsec1
       call mksecond(zsec1,0.)
      endif
!$taf STORE dqdt,dtdt,dudt,dvdt = inner1, REC = 1+nstep2
      if(nvdiff > 0) call vdiff
      if(ntime == 1) then
       call mksecond(zsec1,zsec1)
       time4vd=time4vd+zsec1
       call mksecond(zsec,zsec)
       time4fl=time4fl+zsec
      endif

      return
      end subroutine fluxstep

!     ===================
!     SUBROUTINE FLUXSTOP
!     ===================

      subroutine fluxstop
      use fluxmod

      if(mypid == NROOT .and. ntime == 1) then
       print*,'*********************************************'
       print*,' CPU usage in FLUXSTEP (ROOT process only):  '
       print*,'    All routines      : ',time4fl,' s'
       print*,'    All surface fluxes: ',time4sf,' s'
       print*,'    Vertical diffusion: ',time4vd,' s'
       print*,'    Transfere coeff   : ',time4tr,' s'
       print*,'    Sensible heat flux: ',time4sh,' s'
       print*,'    Evaporation       : ',time4ev,' s'
       print*,'    Wind stress       : ',time4st,' s'
       print*,'*********************************************'
      endif

      return
      end subroutine fluxstop

!
!     =================
!     SUBROUTINE SURFLX
!     =================
!
      subroutine surflx
      use fluxmod
!
!     local parameters
!

      parameter(zumin=1.)  ! minimum wind speed

!
!     local arrays
!

      real zabsu2(NHOR)    ! squared wind speed (m2/s2)
      real znl(NHOR)       ! z of lowermost level (m)
      real zbz0(NHOR)      ! z/z0
      real zri(NHOR)       ! bulk richardson number
      real zrifh(NHOR)     ! factor for heat flux transfere coef.
      real zrifm(NHOR)     ! factor for momentum flux transfere coef.
!
!     allocatable arrays for dbug diagnostics
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:)

      if(ntime == 1) call mksecond(zsec,0.)

!
!*    set some const.:
!

      zexp=-akap
      zlnsig=ALOG(sigma(NLEV))

!
!*    surface air temperature
!

      if(ntsa==1) then
       dtsa(:)=dt(:,NLEV)*sigma(NLEV)**zexp
      elseif(ntsa==2) then
       dtsa(:)=dt(:,NLEV)*sigma(NLEV)**zexp                             &
     &        *(1.+(1./rdbrv-1.)*dq(:,NLEV))
      else
       write(*,*) '!ERROR: wrong ntsa in fluxmod (only 1 or 2 valied)!'
       stop
      end if

!
!*    calculate transfere coef:
!
!     windspeed (squared):
!

! test, FastOpt      zabsu2(:)=AMAX1(zumin,du(:,NLEV)*du(:,NLEV)+dv(:,NLEV)*dv(:,NLEV))
      zabsu2(:)=zumin + du(:,NLEV)*du(:,NLEV)+dv(:,NLEV)*dv(:,NLEV) ! test, FastOpt

!
!     z off lowermost layer
!

      znl(:)=-gascon*0.5*(dt(:,NLEV)+dtsa(:))*zlnsig/ga

!
!     z/z0
!

      zbz0(:)=znl(:)/dz0(:)

!
!     bulk richardson number
!

      if(ntsa==1) then
       zri(:)=ga*znl(:)*(dtsa(:)-dt(:,NLEP))/(zabsu2(:)*dtsa(:))
      elseif(ntsa==2) then
       zri(:)=ga*znl(:)/(zabsu2(:)*dtsa(:))                             &
     &       *(dtsa(:)-dt(:,NLEP)*(1.+(1./rdbrv-1.)*dq(:,NLEP)))
      endif
!$taf LOOP = parallel
      do jhor=1,NHOR

       zkblnz2=(vonkarman/ALOG(zbz0(jhor)+1.))**2

!
!     richardson number dependent factor
!     see ECHAM REPORT 218
!     (unstable zrifh over ocean according to Miller et al. 1992)
!
! FastOpt       if(zri(jhor) <= 0.) then
       if(zri(jhor) < 0.) then ! changed, FastOpt
        zdenom=1.+3.*vdiff_c*vdiff_b                                    &
     &        *sqrt(-1.*zri(jhor)*(zbz0(jhor)+1.))*zkblnz2
        zrifm(jhor)=1.-2.*vdiff_b*zri(jhor)/zdenom
        if(dls(jhor) < 1.) then
! FastOpt: zrifh is 1 for zri==0, but not differentiable(?).
         zdth=-(dtsa(jhor)                                              &
     &         -dt(jhor,NLEP)*(1.+(1./rdbrv-1.)*dq(jhor,NLEP)))
         zrifh(jhor)=(1.+(0.0016*zdth**(1./3.)/SQRT(zabsu2(jhor))       &
     &                  /zkblnz2)**1.25)**0.8
        else
         zrifh(jhor)=1.-3.*vdiff_b*zri(jhor)/zdenom
        endif
       else
        zdenom=SQRT(1+vdiff_d*zri(jhor))
        zrifh(jhor)=1./(1.+3.*vdiff_b*zri(jhor)*zdenom)
        zrifm(jhor)=1./(1.+2.*vdiff_b*zri(jhor)/zdenom)
       endif

!
!     transfere coeff.
!

       ztrans=SQRT(zabsu2(jhor))*zkblnz2
       dtransh(jhor)=ztrans*zrifh(jhor)
       dtransm(jhor)=ztrans*zrifm(jhor)

! guokun added  moved to subroutine below
!       dtransh(jhor)=alpha1m(jhor)*ztrans*zrifh(jhor)
!       dtransm(jhor)=beta1m(jhor)*ztrans*zrifm(jhor)
      enddo

!
!*    dbug diagnostic output if nprint=2 (pumamod)
!

      if(nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       call mpgagp(zprf1,zabsu2,1)
       call mpgagp(zprf2,znl,1)
       call mpgagp(zprf3,zri,1)
       call mpgagp(zprf4,zrifh,1)
       call mpgagp(zprf5,zrifm,1)
       call mpgagp(zprf6,dtransh,1)
       call mpgagp(zprf7,dtransm,1)
       if(mypid==NROOT) then
        print*,'In Surflx:'
        print*,'abs(u)= ',SQRT(zprf1(nprhor)),' znl= ', zprf2(nprhor)
        print*,'ri-number= ',zprf3(nprhor)
        print*,'rifh= ',zprf4(nprhor),' rifm= ', zprf5(nprhor)
        print*,'transh= ',zprf6(nprhor),' transm= ', zprf7(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
       deallocate(zprf7)
      endif

      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4tr=time4tr+zsec
       call mksecond(zsec,0.)
      endif

!
!*    calc. fluxes and tendencies
!

      if (nstress > 0 ) call mkstress   ! wind stress
      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4st=time4st+zsec
       call mksecond(zsec,0.)
      endif
      if (nshfl   > 0 ) call mkshfl     ! sensible heat flux
      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4sh=time4sh+zsec
       call mksecond(zsec,0.)
      endif
      if (nevap   > 0 ) call mkevap     ! evaporation/latent heat flux
      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4ev=time4ev+zsec
      endif

      return
      end subroutine surflx

!     ===================
!     SUBROUTINE MKSTRESS
!     ===================

      subroutine mkstress
      use fluxmod
!
!*    calculate u,v tendencies due to surface wind stress
!
!     implicit formulation is used:
!
!     (u(t+dt)-u(t))/dt=g/ps/dsigma* (cd*abs(u)*rho*u(t+dt))
!     => u(t+dt) => dudt => stress=cd*abs(u)*rho*u(t+dt)
!
!
!     local arrays
!

      real zun(NHOR),zvn(NHOR)         ! u,v at new time level
      real zkdiff(NHOR)                ! transfere coefficients

!
!     allocatable arrays for dbug diagnostics
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:)

!
!*    set some const.
!

      zdelt=delt2/ww
      zkonst1=ga*zdelt/(gascon*dsigma(NLEV))
      zkonst2=dsigma(NLEV)/zdelt/ga

!
!*    modified transfere coeficient (incl. some const.)
!

      zkdiff(:)=zkonst1*gamma1m(:)*dtransm(:)/dt(:,NLEP)

!
!*    wind at t+dt
!

      zun(:)=du(:,NLEV)/(1.+zkdiff(:))
      zvn(:)=dv(:,NLEV)/(1.+zkdiff(:))

!
!*    add tendencies
!

      dudt(:,NLEV)=dudt(:,NLEV)+(zun(:)-du(:,NLEV))/zdelt
      dvdt(:,NLEV)=dvdt(:,NLEV)+(zvn(:)-dv(:,NLEV))/zdelt

!
!*    wind stress
!
      dtaux(:)=dp(:)*zkonst2*zkdiff(:)*zun(:)
      dtauy(:)=dp(:)*zkonst2*zkdiff(:)*zvn(:)

!
!*    ustar**3
!

      dust3(:)=SQRT(dtaux(:)**2+dtauy(:)**2)*gascon*dt(:,NLEP)/dp(:)
      dust3(:)=SQRT(dust3(:))
      dust3(:)=dust3(:)**3

!
!*    dbug diagnostic output if nprint=2 (pumamod)
!

      if(nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       call mpgagp(zprf1,zun,1)
       call mpgagp(zprf2,zvn,1)
       call mpgagp(zprf3,du(1,NLEV),1)
       call mpgagp(zprf4,dv(1,NLEV),1)
       call mpgagp(zprf5,dtaux,1)
       call mpgagp(zprf6,dtauy,1)
       call mpgagp(zprf7,dust3,1)
       if(mypid==NROOT) then
        print*,'In mkstress:'
        print*,'u-new= ',zprf1(nprhor),' v-new= ', zprf2(nprhor)
        print*,'du/dt (m/s2)= ',(zprf1(nprhor)-zprf3(nprhor))/zdelt     &
     &        ,' dv/dt (m/s2)= ',(zprf2(nprhor)-zprf4(nprhor))/zdelt
        print*,'taux= ',zprf5(nprhor),' tauy= ', zprf6(nprhor)
        print*,'ust3= ',zprf7(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
       deallocate(zprf7)
      endif

      return
      end subroutine mkstress

!     =================
!     SUBROUTINE MKSHFL
!     =================

      subroutine mkshfl
      use fluxmod
!
!
!*    calculate t tendency due to surface sensible heat flux
!
!     implicit formulation is used:
!
!     (t(t+dt)-t(t))/dt=-g/cp/ps/dsigma*(cp*cd*abs(u)*rho*(t'(t+dt)-ts)
!     => t(t+dt) => dtdt => sensible heat flux=cp*cd*abs(u)*rho*(t'(t+dt)-ts)
!     =>d(sensible heat flux)/dts (needed for ts calculation)
!     t'= t*(1+(Rv/rd-1)*q)/sigma(NLEV)**akap
!
!
!     local arrays
!

      real ztn(NHOR)      ! new temperature
      real zkdiff(NHOR)   ! modified transfere coef.
      real zfac(NHOR)     ! factor to extrapolate temp. to surface
                          ! (needs to be independent of temp)
!
!     allocatable arrays for dbug diagnostics
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)

!
!*    set some const.
!

      zexp=-akap
      zdelt=delt2/ww
      zkonst1=ga*zdelt/(gascon*dsigma(NLEV))
      zkonst2=dsigma(NLEV)/zdelt/ga
      zfac(:)=sigma(NLEV)**(zexp)

!
!*    modified transfere coeficien (incl. some const.)
!

      zkdiff(:)=zkonst1*alpha1m(:)*dtransh(:)/dt(:,NLEP)

!
!*    temperatue at t+dt
!

      ztn(:)=(dt(:,NLEV)+zkdiff(:)*dt(:,NLEP))/(1.+zkdiff(:)*zfac(:))

!
!*    add tendency
!

      dtdt(:,NLEV)=dtdt(:,NLEV)+(ztn(:)-dt(:,NLEV))/zdelt

!
!*    sensible heat flux and dshfl/dsurface_temperature
!

      dshfl(:)=zkonst2*zkdiff(:)*dp(:)*acpd*AMAX1(1.,1.+ADV*dq(:,NLEV)) &
     &        *(ztn(:)*zfac(:)-dt(:,NLEP))
      dshdt(:)=-1.*zkonst2*zkdiff(:)*dp(:)                              &
     &        *acpd*AMAX1(1.,1.+ADV*dq(:,NLEV))
!
      dtsa(:)=ztn(:)*zfac(:)
!
!
!*    dbug diagnostic output if nprint=2 (pumamod)
!
      if(nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       call mpgagp(zprf1,ztn,1)
       call mpgagp(zprf2,dt(1,NLEV),1)
       call mpgagp(zprf3,zfac,1)
       call mpgagp(zprf4,dshfl,1)
       call mpgagp(zprf5,dshdt,1)
       if(mypid==NROOT) then
        zzts=zprf2(nprhor)*zprf3(nprhor)
        zztsn=zprf1(nprhor)*zprf3(nprhor)
        print*,'In mkshfl:'
        print*,'t-new= ',zprf1(nprhor),' tsa(old)= ',zzts,' tsa(new)= ' &
     &        ,zztsn
        print*,'dt/dt (K/s)= ',(zprf1(nprhor)-zprf2(nprhor))/zdelt
        print*,'shfl= ',zprf4(nprhor),' dshdt= ', zprf5(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
      endif

!
!     franks dbug
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEM,3)=0.
       dgp3d(:,NLEV,3)=(ztn(:)-dt(:,NLEV))/zdelt
      end if

      return
      end subroutine mkshfl

!     =================
!     SUBROUTINE MKEVAP
!     =================

      subroutine mkevap
      use fluxmod
!
!*    calculate q tendency due to surface evaporation
!
!     implicit formulation is used:
!
!     (q(t+dt)-q(t))/dt=-g/ps/dsigma*(rhs*cd*abs(u)*rho*(q(t+dt)-qs)
!     => q(t+dt) => dqdt => evaporation=rhs*cd*abs(u)*rho*(q(t+dt)-qs)
!     =>latent heat flux=>d(latent heat flux)/dts (needed for ts calculation)
!     rhs= soil wetness factor
!
!     local arrays
!
!
      real zqn(NHOR)      ! new q
      real zkdiff(NHOR)   ! modified transfere coef.

!
!     allocatable arrays for dbug diagnostics
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)

!
!*    set some const.
!

      zdelt=delt2/ww
      zkonst1=ga*zdelt/(gascon*dsigma(NLEV))
      zkonst2=dsigma(NLEV)/zdelt/ga

!
!*    modified transfere coeficient (incl. some const.)
!

      zkdiff(:)=drhs(:)*zkonst1*beta1m(:)*dtransh(:)/dt(:,NLEP)

!
!*    q at t+dt
!

      zqn(:)=AMAX1(dq(:,NLEV)                                           &
     &            ,(dq(:,NLEV)+zkdiff(:)*dq(:,NLEP))/(1.+zkdiff(:)))

!
!*    add tendency
!

      dqdt(:,NLEV)=dqdt(:,NLEV)+(zqn(:)-dq(:,NLEV))/zdelt

!
!*    evaporation, latent heat flux and dlhfl/dtemperature
!

      devap(:)=-dp(:)*zkonst2/1000.*(zqn(:)-dq(:,NLEV))

      where(dt(:,NLEP) > TMELT )
       dlhfl(:)=devap(:)*ALV*1000.
       dlhdt(:)=-1.*ALV*zkdiff(:)*zkonst2*dp(:)                         &
     &         *ra2*(TMELT-ra4)*dq(:,NLEP)/(dt(:,NLEP)-ra4)**2
      elsewhere
       dlhfl(:)=devap(:)*ALS*1000.
       dlhdt(:)=-1.*ALS*zkdiff(:)*zkonst2*dp(:)                         &
     &         *ra2*(TMELT-ra4)*dq(:,NLEP)/(dt(:,NLEP)-ra4)**2
      endwhere
      where(dlhfl(:) == 0.) dlhdt(:)=0.
!
!*    dbug diagnostic output if nprint=2 (pumamod)
!
      if(nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       call mpgagp(zprf1,zqn,1)
       call mpgagp(zprf2,dq(1,NLEV),1)
       call mpgagp(zprf3,devap,1)
       call mpgagp(zprf4,dlhfl,1)
       call mpgagp(zprf5,dlhdt,1)
       if(mypid==NROOT) then
        print*,'In mkevap:'
        print*,'q-new= ',zprf1(nprhor)                                  &
     &        ,' dq/dt (1/s)= ',(zprf1(nprhor)-zprf2(nprhor))/zdelt
        print*,'lhfl= ',zprf4(nprhor),' dlhdt= ', zprf5(nprhor)
        print*,'evap (mm/day)= ',1000.*ntspd*0.5*zdelt*zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
      endif

!
!     franks dbug
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEM,12)=0.
       dgp3d(:,NLEV,12)=(zqn(:)-dq(:,NLEV))/zdelt*dp(:)*dsigma(NLEV)    &
     &                 /ga/1000.
      end if

      return
      end subroutine mkevap

!     ================
!     SUBROUTINE VDIFF
!     ================

      subroutine vdiff
      use fluxmod
      use rainmod, only : nprl,nprc
!
!     calculate t,q,u,v tendencies due to vertical diffusion
!     using the ECHAM semi-implicit scheme
!
      parameter(ztscal=250.)
      parameter(zumin=1.)
!
      real zdtdt(NHOR,NLEV)
      real zdudt(NHOR,NLEV)
      real zdvdt(NHOR,NLEV)
      real zdqdt(NHOR,NLEV)
!
      real zkdiffm(NHOR,NLEM),zkdiffh(NHOR,NLEM)
      real zabsu(NHOR,NLEV)
      real ztn(NHOR,NLEV)
      real zun(NHOR,NLEV)
      real zvn(NHOR,NLEV)
      real zqn(NHOR,NLEV)
      real zt(NHOR,NLEV)
      real zu(NHOR,NLEV)
      real zv(NHOR,NLEV)
      real zq(NHOR,NLEV)
      real zebs(NHOR,NLEM)
      real zskap(NLEV),zskaph(NLEV)

      zlamh=vdiff_lamm*SQRT(3.*vdiff_d*0.5)
      zdelt=delt2/ww
      zkonst1=ga*zdelt/gascon
      zkonst2=zkonst1*ga/gascon
      zkonst3=zkonst2*ga/gascon
!
      zskap(:)=sigma(:)**akap
      zskaph(:)=sigmah(:)**akap
!
      zt(:,:)=dt(:,1:NLEV)+dtdt(:,1:NLEV)*zdelt
      zq(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*zdelt
      zu(:,:)=du(:,1:NLEV)+dudt(:,1:NLEV)*zdelt
      zv(:,:)=dv(:,1:NLEV)+dvdt(:,1:NLEV)*zdelt

!
!     0) modified diffusion coefficents
!
      zabsu(:,:)=SQRT(zu(:,1:NLEV)**2+zv(:,1:NLEV)**2)
!
      do jlev=1,NLEM
       jlp=jlev+1
       zzlev=-gascon*ztscal*ALOG(sigmah(jlev))/ga
       zdz=gascon*ztscal*ALOG(sigma(jlp)/sigma(jlev))/ga
       zmixh=zlamh*vonkarman*zzlev/(zlamh+vonkarman*zzlev)
       zmixm=vdiff_lamm*vonkarman*zzlev/(vdiff_lamm+vonkarman*zzlev)
       zmixh2=zmixh*zmixh
       zmixm2=zmixm*zmixm
       zrdsig=1./(sigma(jlp)-sigma(jlev))
       zfac=zkonst3*sigmah(jlev)*sigmah(jlev)*sigmah(jlev)              &
     &     *zrdsig
       zzkfac=SQRT((((((zzlev+zdz)/zzlev)**(1./3.)-1.)/zdz)**3)/zzlev)
!$taf LOOP = parallel
       do jhor=1,NHOR
        zrthl=(dsigma(jlp)+dsigma(jlev))                                &
     &       /(zt(jhor,jlev)*dsigma(jlp)+zt(jhor,jlp)*dsigma(jlev))
        zdvds=AMAX1(zumin,ABS(zabsu(jhor,jlev)-zabsu(jhor,jlp)))        &
     &       *zrdsig
        zqf1=1.+(1./rdbrv-1.)*zq(jhor,jlev)
        zqf2=1.+(1./rdbrv-1.)*zq(jhor,jlp)
        zdthds=(zt(jhor,jlev)*zqf1/zskap(jlev)                          &
     &         -zt(jhor,jlp)*zqf2/zskap(jlp))                           &
     &        *zrdsig
        zri=zskaph(jlev)*gascon*zdthds                                  &
     &     /(sigmah(jlev)*zdvds*zdvds)
!FastOpt        if(zri <= 0.) then
        if(zri < 0.) then ! changed, FastOpt
         zdenomh=1.+3.*vdiff_c*vdiff_b*sqrt(-1.*zri)*zzkfac*zmixh2
         zdenomm=1.+3.*vdiff_c*vdiff_b*sqrt(-1.*zri)*zzkfac*zmixm2
         zrifh=1.-3.*vdiff_b*zri/zdenomh
         zrifm=1.-2.*vdiff_b*zri/zdenomm
        else
         zdenom=SQRT(1+vdiff_d*zri)
         zrifh=1./(1.+3.*vdiff_b*zri*zdenom)
         zrifm=1./(1.+2.*vdiff_b*zri/zdenom)
        endif
        zkdiffh(jhor,jlev)=zfac*zmixh2*zrifh*zdvds*zrthl*zrthl*zrthl
        zkdiffm(jhor,jlev)=zfac*zmixm2*zrifm*zdvds*zrthl*zrthl*zrthl
       enddo
      enddo
!$TAF STORE zkdiffh(:,:),zkdiffm(:,:) = inner1, REC = 1 + nstep2

!
!     2. semi implicit scheme
!
!     2a momentum
!
!     top layer elimination
!
      zebs(:,1)=zkdiffm(:,1)/(dsigma(1)+zkdiffm(:,1))
      zun(:,1)=dsigma(1)*zu(:,1)/(dsigma(1)+zkdiffm(:,1))
      zvn(:,1)=dsigma(1)*zv(:,1)/(dsigma(1)+zkdiffm(:,1))
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       zebs(:,jlev)=zkdiffm(:,jlev)                                     &
     &             /(dsigma(jlev)+zkdiffm(:,jlev)                       &
     &              +zkdiffm(:,jlem)*(1.-zebs(:,jlem)))
       zun(:,jlev)=(zu(:,jlev)*dsigma(jlev)+zkdiffm(:,jlem)*zun(:,jlem))&
     &            /(dsigma(jlev)+zkdiffm(:,jlev)                        &
     &             +zkdiffm(:,jlem)*(1.-zebs(:,jlem)))
       zvn(:,jlev)=(zv(:,jlev)*dsigma(jlev)+zkdiffm(:,jlem)*zvn(:,jlem))&
     &            /(dsigma(jlev)+zkdiffm(:,jlev)                        &
     &             +zkdiffm(:,jlem)*(1.-zebs(:,jlem)))
      enddo
!$TAF STORE zebs(:,:),zun(:,1:NLEM),zvn(:,1:NLEM) = inner1, REC = 1 + nstep2
!
!     bottom layer elimination
!
      zun(:,NLEV)=(zu(:,NLEV)*dsigma(NLEV)+zkdiffm(:,NLEM)*zun(:,NLEM)) &
     &           /(dsigma(NLEV)+zkdiffm(:,NLEM)*(1.-zebs(:,NLEM)))
      zvn(:,NLEV)=(zv(:,NLEV)*dsigma(NLEV)+zkdiffm(:,NLEM)*zvn(:,NLEM)) &
     &           /(dsigma(NLEV)+zkdiffm(:,NLEM)*(1.-zebs(:,NLEM)))
!$TAF STORE zebs(:,:),zun(:,:),zvn(:,:) = inner1, REC = 1 + nstep2
!
!     back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       zun(:,jlev)=zun(:,jlev)+zebs(:,jlev)*zun(:,jlep)
       zvn(:,jlev)=zvn(:,jlev)+zebs(:,jlev)*zvn(:,jlep)
      enddo
!
!     tendencies
!
      zdudt(:,1:NLEV)=(zun(:,1:NLEV)-zu(:,1:NLEV))/zdelt
      zdvdt(:,1:NLEV)=(zvn(:,1:NLEV)-zv(:,1:NLEV))/zdelt
      dudt(:,1:NLEV)=dudt(:,1:NLEV)+zdudt(:,1:NLEV)
      dvdt(:,1:NLEV)=dvdt(:,1:NLEV)+zdvdt(:,1:NLEV)
   if (nevap.gt.0.or.nprl.gt.0.or.nprc.gt.0) then
!
!     2b moisture
!
!     top layer elimination
!
      zebs(:,1)=zkdiffh(:,1)/(dsigma(1)+zkdiffh(:,1))
      zqn(:,1)=dsigma(1)*zq(:,1)/(dsigma(1)+zkdiffh(:,1))
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       zebs(:,jlev)=zkdiffh(:,jlev)                                     &
     &             /(dsigma(jlev)+zkdiffh(:,jlev)                       &
     &              +zkdiffh(:,jlem)*(1.-zebs(:,jlem)))
       zqn(:,jlev)=(zq(:,jlev)*dsigma(jlev)+zkdiffh(:,jlem)*zqn(:,jlem))&
     &            /(dsigma(jlev)+zkdiffh(:,jlev)                        &
     &             +zkdiffh(:,jlem)*(1.-zebs(:,jlem)))
      enddo
!$TAF STORE zebs(:,:),zqn(:,1:NLEM) = inner1, REC = 1 + nstep2
!
!     bottom layer elimination
!
      zqn(:,NLEV)=(zq(:,NLEV)*dsigma(NLEV)+zkdiffh(:,NLEM)*zqn(:,NLEM)) &
     &           /(dsigma(NLEV)+zkdiffh(:,NLEM)*(1.-zebs(:,NLEM)))
!
!     back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       zqn(:,jlev)=zqn(:,jlev)+zebs(:,jlev)*zqn(:,jlep)
      enddo
!
!     tendencies
!
      zdqdt(:,1:NLEV)=(zqn(:,1:NLEV)-zq(:,1:NLEV))/zdelt
      do jlev=1,NLEV
       dqdt(:,jlev)=dqdt(:,jlev)+zdqdt(:,jlev)
      enddo
   end if
!
!     2c potential temperature
!
      do jlev=1,NLEM
       zkdiffh(:,jlev)=zkdiffh(:,jlev)*zskaph(jlev)
      enddo
!
!     semi implicit scheme
!
!     top layer elimination
!
      zebs(:,1)=zkdiffh(:,1)/(dsigma(1)+zkdiffh(:,1)/zskap(1))
      ztn(:,1)=dsigma(1)*zt(:,1)/(dsigma(1)+zkdiffh(:,1)/zskap(1))
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       zebs(:,jlev)=zkdiffh(:,jlev)                                     &
     &             /(dsigma(jlev)+(zkdiffh(:,jlev)                      &
     &              +zkdiffh(:,jlem)*(1.-zebs(:,jlem)/zskap(jlem)))     &
     &              /zskap(jlev))
       ztn(:,jlev)=(zt(:,jlev)*dsigma(jlev)                             &
     &             +zkdiffh(:,jlem)/zskap(jlem)*ztn(:,jlem))            &
     &            /(dsigma(jlev)+(zkdiffh(:,jlev)                       &
     &              +zkdiffh(:,jlem)*(1.-zebs(:,jlem)/zskap(jlem)))     &
     &              /zskap(jlev))
      enddo
!
!     bottom layer elimination
!
!$TAF STORE zebs,ztn(:,NLEM) = inner1, REC = 1 + nstep2
      ztn(:,NLEV)=(zt(:,NLEV)*dsigma(NLEV)                              &
     &            +zkdiffh(:,NLEM)*ztn(:,NLEM)/zskap(NLEM))             &
     &           /(dsigma(NLEV)+zkdiffh(:,NLEM)/zskap(NLEV)             &
     &                         *(1.-zebs(:,NLEM)/zskap(NLEM)))
!
!     back-substitution
!
!$taf store ztn = inner1, REC = 1 + nstep2
      do jlev=NLEM,1,-1
       jlep=jlev+1
       ztn(:,jlev)=ztn(:,jlev)+zebs(:,jlev)*ztn(:,jlep)/zskap(jlep)
      enddo
!
!     tendencies
!
      zdtdt(:,1:NLEV)=(ztn(:,1:NLEV)-zt(:,1:NLEV))/zdelt
      dtdt(:,1:NLEV)=dtdt(:,1:NLEV)+zdtdt(:,1:NLEV)

!
!     franks dbug
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEV,4)=zdtdt(:,1:NLEV)
       do jlev=1,NLEV
        dgp3d(:,jlev,13)=zdqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      end if

      return
      end subroutine vdiff
