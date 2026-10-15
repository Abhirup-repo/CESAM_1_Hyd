      subroutine plasim_ini(myProcId,myworld)
      integer :: myProcId
      integer :: myworld
!
      call mpstart(myworld)
      call prolog
      call init_master
!
      return
      end
!
      subroutine plasim_step(ksteps,kstepo,ikey_dynamics,myProcId)
        use cplmod ! for ad-storage
        implicit none
        integer :: ksteps, kstepo,ikey_dynamics
        integer :: myProcId
!
      call getmitocean
      call master(ksteps,kstepo,ikey_dynamics)
!$TAF STORE ahfla,apmea,aroffa,aswra,atauxa,atauya,nhfla,npmea,nroffa,nswra,ntaua = innern, REC = ikey_dynamics
!
      return
      end
!
      subroutine plasim_stop
!
      call epilog
!
      return
      end

!     ======================
!     SUBROUTINE GETMITOCEAN
!     ======================

      subroutine getmitocean
      use pumamod, only: NLON,NLAT,NHOR,dt,dq,dz0,dp,drhs,dtaux,dtauy,dalb,ga,gascon,NLEP,mypid,NROOT,ra1,ra2,ra4,rdbrv,tmelt
      use seamod, only: dls_sea=>dls,dsst,dts_sea=>dts,dicec,dqs_sea=>dqs,charnock,albice,albsea,drhsice,drhssea,dz0ice,dz0sea
      use icemod, only: nmit,osst,osss
      use fluxmod, only: nstress
      implicit none
!
!     get oceanic variables from the coupler (if mit-coupling is
!     switched on.
!
!
      real :: zz0(NHOR) = 0.
      real :: zsst(NLON,NLAT)
      real :: zsss(NLON,NLAT)
      real :: zts_sea(NHOR) ! FastOpt, workaround
!
!     get new mit variables
!
!$TAF INCOMPLETE dalb, drhs, dz0, dqs_sea, dt
      if(nmit > 0) then
       if(mypid == NROOT) then
!
!      get ocean state
!
        call cgetssta(zsst)
        call cgetsssa(zsss)
       
       endif
       call mpscgp(zsst,osst,1)
       call mpscgp(zsss,osss,1)
!
!      overwrite sst and tsurf
!
!$TAF incomplete dq,dt,drhs
       where (dls_sea(:) < 0.5)
        dsst(:)=osst(:)
        zts_sea(:) = dts_sea(:) ! FastOpt, workaround
 ! FastOpt, workaround       dts_sea(:)=dicec(:)*dts_sea(:)+(1.-dicec(:))*dsst(:)
        dts_sea(:)=dicec(:)*zts_sea(:)+(1.-dicec(:))*dsst(:) ! FastOpt, workaround
!
!     set dependent surface variables
!
        dt(:,NLEP)=dts_sea(:)
        dqs_sea(:)=rdbrv*ra1*EXP(ra2*(dt(:,NLEP)-TMELT)                     &
     &        /(dt(:,NLEP)-ra4))/dp(:)
        dqs_sea(:)=dqs_sea(:)/(1.-(1./rdbrv-1.)*dqs_sea(:))
        dq(:,NLEP)=dqs_sea(:)
        drhs(:)=1.
       endwhere
!
!     z0 after charnock
!
       zz0 = 0. ! FastOpt
       if (nstress>0) then ! added, FastOpt:  zz0/dz0 positive and defined
        where(dls_sea(:) < 0.5)
         zz0(:)=charnock*SQRT(dtaux(:)**2+dtauy(:)**2)*gascon*dt(:,NLEP) &
     &         /(ga*dp(:))
         zz0(:)=AMAX1(zz0(:),dz0sea)
        end where ! FastOpt
       else ! added, FastOpt
        where(dls_sea(:) < 0.5)
          zz0(:) = dz0sea ! added, FastOpt
        end where
       end if ! added, FastOpt
       where(dls_sea(:) < 0.5)
        drhs(:)=drhssea*(1.-dicec(:))+drhsice*dicec(:)
        dalb(:)=albsea*(1.-dicec(:))                                     &
     &         +dicec(:)*AMIN1(albice,0.5+0.025*(273.-dts_sea(:)))
        dz0(:)=zz0(:)*(1.-dicec(:))+dz0ice*dicec(:)
       endwhere
      endif
!
      return
      end subroutine getmitocean

