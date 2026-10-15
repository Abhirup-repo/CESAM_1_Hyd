!     ***************************************************
!     * GUIMOD - Graphical User Interface Routines      *
!     * 21-Sep-2006 - Edilbert Kirk                     *
!     ***************************************************
!     * This file contains all interface routines for   *
!     * communication between model (PUMA or PLASIM)    *
!     * and the GUI routines in file "pumax.c".         *
!     * This file is identical for both models, if you  *
!     * make changes, make sure, that these either      *
!     * affect both models in a proper way or use the   *
!     * if (model == PUMA)    ...  endif                *
!     * if (model == PLASIM)  ...  endif                *
!     * statements. After changes copy the new version, *
!     * so that ../puma/src/guimod.f90  and             *
!     *       ../plasim/src/guimod.f90 are identical.   *
!     ***************************************************

      subroutine guistart
      use pumamod
      integer :: idebug = 0

      if (ngui == 0) return

      call initgui(model,idebug)

      return
      end subroutine guistart

!     ==================
!     SUBROUTINE GUISTOP
!     ==================

      subroutine guistop
      use pumamod

      if (mypid == NROOT .and. ngui > 0) call guiclose

      return
      end subroutine guistop

!     =======================
!     SUBROUTINE GUISTEP_PUMA
!     =======================

      subroutine guistep_puma
      use pumamod

      interface 
         integer(kind=4) function iguistep(parc,idatim)
            real   (kind=4), intent(inout) :: parc(*)
            integer(kind=4), intent(in)    :: idatim(6)
         end function iguistep
      end interface

      integer (kind=4) idatim(6)

      if (ngui == 0) return

      if (mypid == NROOT) then
         parc(1)   = disp
         parc(2)   = dtep * CT;
         parc(3)   = dtns * CT;
         crap(:)   = parc(:)
         idatim(:) = ndatim(:)
         nshutdown = iguistep(parc,idatim)    ! GUI event handler
         if (parc(1) /= crap(1)) call change_disp(1)
         if (parc(2) /= crap(2)) call change_dtep(2)
         if (parc(3) /= crap(3)) call change_dtns(3)
      endif
   
      call mpbci(nshutdown)
      call mpbcr(disp)
   
      call mpbcl(ldtep)
      if (ldtep) then
         call mpbcr(dtep)            ! broadcast changed value
         call mpscsp(sr1,srp1,NLEV)  ! scatter   changed array
         ldtep = .false.
      endif 
   
      call mpbcl(ldtns)
      if (ldtns) then
         call mpbcr(dtns)            ! broadcast changed value
         call mpscsp(sr2,srp2,NLEV)  ! scatter   changed array
         ldtns = .false.
      endif 
   
      return
      end subroutine guistep_puma

!     =========================
!     SUBROUTINE GUISTEP_PLASIM
!     =========================

      subroutine guistep_plasim
      use radmod

      interface 
         integer(kind=4) function iguistep(p2,k3)
            real   (kind=4), intent(inout) :: p2(*)
            integer(kind=4), intent(in)    :: k3(6)
         end function iguistep
      end interface

      integer(kind=4) :: idatim(6)

      if (mypid == NROOT) then
         parc(1)   = co2 
         parc(2)   = gsol0
         parc(3)   = dawn
         crap(:)   = parc(:)
         idatim(:) = ndatim(:)
         nshutdown = iguistep(parc,idatim)    ! GUI event handler
         if (parc(1) /= crap(1)) call change_co2(1)
         if (parc(2) /= crap(2)) call change_gsol0(2)
         if (parc(3) /= crap(3)) call change_dawn(3)
      endif
   
      call mpbci(nshutdown)
      call mpbcr(co2)
      call mpbcr(gsol0)
      call mpbcr(dawn)
   
      return
      end subroutine guistep_plasim

!     ================
!     SUBROUTINE GUIPS
!     ================

      subroutine guips(f)
      use pumamod
      real :: f(NLON,NLAT)
      real (kind=4) :: x(NLON+1,NLAT)

      if (ngui == 0) return
      if (mypid == NROOT) then
         x(1:NLON,:) = f(:,:) * PSURF * 0.01 ! [hPa]
         x(NLON+1,:) = x(1,:)
         call guiput("GP" // char(0) ,x,NLON+1,NLAT,1)
      endif
      return
      end

!     =================
!     SUBROUTINE GUIHOR
!     =================

      subroutine guihor(yname,f,klev,pm,pa)
      use pumamod
      character (len=*) :: yname
      real :: f(NLON,NLPP,klev)
      real :: z(NLON,NLAT,klev)
      real (kind=4) :: x(NLON+1,NLAT,klev)

      if (ngui == 0) return

      call mpgagp(z,f,klev)
      if (mypid == NROOT) then
         x(1:NLON,:,:) = z(:,:,:) * pm + pa
         x(NLON+1,:,:) = x(1,:,:)
         call guiput(yname,x,NLON+1,NLAT,klev)
      endif
      return
      end

!     ================
!     SUBROUTINE GUIGV
!     ================

      subroutine guigv(yname,f)
      use pumamod
      character (len=*) :: yname
      real :: f(NLON,NLPP,NLEV)
      real :: z(NLON,NLAT,NLEV)
      real (kind=4) :: x(NLON+1,NLAT,NLEV)

      if (ngui == 0) return
      
      call mpgagp(z,f,NLEV)
      if (mypid == NROOT) then
         do jlat = 1 , NLAT
            x(1:NLON,jlat,:) = z(:,jlat,:) * CV * rcs(jlat)
         enddo
         x(NLON+1,:,:) = x(1,:,:)
         call guiput(yname,x,NLON+1,NLAT,NLEV)
      endif
      return
      end

!     ======================
!     SUBROUTINE CHANGE_DISP
!     ======================

      subroutine change_disp(k)
      use pumamod
      write (*,7000) nstep,'DISP  ',crap(k),parc(k)
      disp = parc(k)
      return
 7000 format('Step',i8,': User changed ',a,' from ',f6.2,' to ',f6.2)
      end

!     ======================
!     SUBROUTINE CHANGE_DTEP
!     ======================

      subroutine change_dtep(k)
      use pumamod
      write (*,7000) nstep,'DTEP  ',crap(k),parc(k)
      dtep = parc(k) / CT;
      zttrop = tgr-dtrop*ALR
      ztps   = (zttrop/tgr)**(GA/(ALR*GASCON))
      do jlev = 1 , NLEV
         zfac = sin(0.5*PI*(sigma(jlev)-ztps)/(1.-ztps))
         if (zfac < 0.0) zfac = 0.0
         sr1(5,jlev) = -2.0/3.0 * sqrt(0.4) * dtep * zfac
      enddo
      ldtep = .true.
      return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
      end

!     ======================
!     SUBROUTINE CHANGE_DTNS
!     ======================

      subroutine change_dtns(k)
      use pumamod
      write (*,7000) nstep,'DTNS  ',crap(k),parc(k)
      dtns = parc(k) / CT;
      zttrop = tgr-dtrop*ALR
      ztps   = (zttrop/tgr)**(GA/(ALR*GASCON))
      do jlev = 1 , NLEV
         zfac = sin(0.5*PI*(sigma(jlev)-ztps)/(1.-ztps))
         if (zfac < 0.0) zfac = 0.0
         sr2(3,jlev) = (1.0 / sqrt(6.0)) * dtns * zfac
      enddo
      ldtns = .true.
      return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
      end

!     =====================
!     SUBROUTINE CHANGE_CO2
!     =====================

      subroutine change_co2(k)
      use pumamod
      write (*,7000) nstep,'CO2   ',crap(k),parc(k)
      co2  = parc(k)
      return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
      end

!     =======================
!     SUBROUTINE CHANGE_GSOL0
!     =======================

      subroutine change_gsol0(k)
      use pumamod
      use radmod
      write (*,7000) nstep,'GSOL0 ',crap(k),parc(k)
      gsol0  = parc(k)
      return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
      end

!     ======================
!     SUBROUTINE CHANGE_DAWN
!     ======================

      subroutine change_dawn(k)
      use radmod
      write (*,7000) nstep,'DAWN ',crap(k),parc(k)
      dawn  = parc(k)
      return
 7000 format('Step',i8,': User changed ',a,' from ',f7.2,' to ',f7.2)
      end
