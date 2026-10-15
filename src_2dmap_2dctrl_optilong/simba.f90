! =======================
! VEGETATION MODULE SIMBA
! =======================

module vegmod
use landmod
implicit none

! parameter for soil component

real, parameter :: vws_crit  = 0.25             ! critical soil wetness

! parameter for vegetation component
! all values are in m-k-s units

! FastOpt real, parameter :: q10       = 2.0              ! q10 value for resp
real :: q10       = 2.0              ! changed, FastOpt

!     * parameters for land surface parameters

real, parameter :: vlai_max  =  6.0               ! max LAI
real, parameter :: vlai_min  =  0.1               ! min LAI for bare soil
real, parameter :: valb_min  =  0.12              ! min albedo for veg
real, parameter :: valb_max  =  0.30              ! max albedo for bare soil
real, parameter :: vsalb_min =  0.35              ! snow albedo for veg
real, parameter :: vwmax_max =  1.0               ! max SWHC for veg
real, parameter :: vwmax_min =  0.01              ! SWHC for bare soil
real, parameter :: vz0_max   =  2.0               ! max roughness for veg
real, parameter :: vz0_min   =  0.05              ! min roughness for bare soil
! FastOpt real, parameter :: cveg_a    =  3.0               ! forest cover conversion
real :: cveg_a    =  3.0 ! changed, FastOpt
! FastOpt real, parameter :: cveg_b    =  1.0               ! forest cover conversion
real :: cveg_b    =  1.0 ! changed, FastOpt
real, parameter :: cveg_c    =  1.570796326794897 ! forest cover conversion
real, parameter :: cveg_k    =  0.5               ! conversion LAI
!FastOpt real, parameter :: ct_crit   =  5.0               ! conversion LAI
real :: ct_crit   =  5.0               ! changed, FastOpt

!* beta factor defined as GPP  = (1 + beta ln(C/C0))
!  beta  = 0.3-0.6, C0 = 360ppm - Harvey 1989

real, parameter :: co2_ref   = 360.     ! co2 reference concentration in ppm
real, parameter :: co2_comp  = 0.       ! light compensation point in ppm
real, parameter :: co2_sens  = 0.1      ! beta factor

real               :: cveg_d, cveg_e         ! forest cover conversion

!     * dynamically calculated land surface parameters

real,dimension(NHOR) :: zforest=   0.0         ! forest cover = dforest(:)
real,dimension(NHOR) :: zwmax  =   0.0         ! bucket depth = dwmax(:)
real,dimension(NHOR) :: dvlai  =   0.          ! leaf area index LAI
real,dimension(NHOR) :: dvz0   =   vz0_min     ! roughness length
real,dimension(NHOR) :: dvveg  =   0.          ! vegetation cover
real,dimension(NHOR) :: dvalb  =   valb_max    ! albedo
real,dimension(NHOR) :: dvrhs  =   0.          ! surface conductance
real,dimension(NHOR) :: dvsoil =   0.          ! soil cover
real,dimension(NHOR) :: zalbsn =   0.          ! snow albedo

!     * output of cumulative fluxes

real,dimension(NHOR) :: zlaim                  ! maximum lai
real,dimension(NHOR) :: zvegm                  ! maximum veg cover
real,dimension(NHOR) :: zveg                   ! veg cover
real,dimension(NHOR) :: zvpd, zft, zfw
real,dimension(NHOR) :: zesat_srf
integer,dimension(NHOR) :: nbiomass
real                 :: zbeta

end module vegmod


!========!
! VEGINI !
!========!

subroutine vegini
use vegmod
implicit none

if (nrestart == 0) then
   where (dls(:) == 0) dforest(:) = 0.0
endif
zfw(:)   = 0. ! FastOpt
zft(:)   = 0. ! FastOpt
nbiomass = 0. ! FastOpt
dvsoil   = 0. ! FastOpt
end subroutine vegini


!=========!
! VEGSTOP !
!=========!

subroutine vegstop
use vegmod
implicit none

return
end subroutine vegstop


!=========!
! VEGSTEP !
!=========!

subroutine vegstep
use vegmod
implicit none

real :: zdelt ! Timestep in seconds

! Initialization and copy of puma fields needed

zdelt  = delt / ww
cveg_d = -atan(-cveg_a/cveg_b)/cveg_c
cveg_e = 1./(1. + cveg_d)
zbeta  = max(0.,1.+co2_sens*log((co2-co2_comp)/(co2_ref-co2_comp)))

! Following plasim arrays are used but not modified
! -------------------------------------------------
! dwatc(:) = soil wetness [m]
! dgppl(:) = GPP light limited
! dswfl(:) = short wave radiation [W/m2]
! devap(:) = surface evaporation (negative)

! Make local copies of some plasim arrays
! They are copied back if NBIOME is set to 1 (interactive)

zforest(:) = dforest(:) ! Forest cover (0.0 - 1.0)
zwmax(:)   = dwmax(:)   ! Bucket depth

dnpp(:)   = 0.0  ! Net   primary production [kg C m-2 s-1]
dgpp(:)   = 0.0  ! Gross primary production [kg C m-2 s-1]
dgppw(:)  = 0.0  ! GPP water limited

zesat_srf(:) = 0. ! added, FastOpt
dgppl(:) = 0. ! added, FastOpt
dlitter(:) = 0. ! added, FastOpt
zlaim(:) = 0. ! added, FastOpt
zvegm(:) = 0. ! added, FastOpt
zvpd(:) = 0. ! added, FastOpt
zft(:) = 0. ! added, FastOpt
dvalb(:) = 0. ! added, FastOpt
!$TAF INCOMPLETE zesat_srf,zvpd,zlaim,zvegm,zfw,zveg,zft,dgppw,dgpp,dnpp,dlitter,dres,dnogrow,zforest,dvsoil
where (dls(:) > 0.0 .and. dglac(:) == 0.0) ! land cell with no glacier

   ! calculate vapor pressure deficit

   zesat_srf(:) = ra1*exp(ra2*(dt(:,NLEP)-TMELT)/(dt(:,NLEP)-ra4))
   zvpd(:)      = max(10.0, zesat_srf(:) - dp(:) * dq(:, NLEV))

   ! calculate structural limitation to leaf area

   zlaim(:) = vlai_max * zforest(:) + vlai_min * (1.0 - zforest(:))
   zvegm(:) = 1.0 - exp(-cveg_k * zlaim(:))

   ! calculate water stress factor

   where (zwmax(:) > 0.0) zfw(:) = dwatc(:) / (zwmax(:) * vws_crit)
   where (zfw(:)   < 0.0) zfw(:) = 0.0
   where (zfw(:)   > 1.0) zfw(:) = 1.0

   ! calculate actual leaf cover from actual intercepted radiation

   zveg(:) = min(zfw(:),zvegm(:))

   ! calculate gross primary productivity

   zft(:) = (dt(:,NLEP) - TMELT) / ct_crit
   where (zft(:) < 0.0) zft(:) = 0.0
   where (zft(:) > 1.0) zft(:) = 1.0

   ! light limited gpp [kg C / m2 /s]

   dgppl(:) = rlue * zbeta * zft(:) * zveg(:) * dswfl(:,NLEP)

   ! water limited gpp [kg C / m2 /s]

   dgppw(:) = co2conv * zveg(:) * (-devap(:)) / zvpd(:) * (0.7 * co2)
   where (dgppw(:) < 0.0) dgppw(:) = 0.0

   ! combine light & water limitations

   dgpp(:) = min(dgppl(:),dgppw(:))

   ! net primary production [kg C / m2 /s]

   dnpp(:) = pgrow(:) * 0.5 * dgpp(:)

   ! litterfall [kg C / m2 /s]

   dlitter(:) = dcveg(:) / tau_veg

   ! heterotrophic respiration [kg C / m2 /s]

   dres(:) = q10**((dt(:,NLEP)-TMELT-10.)/10.) * dcsoil(:)/tau_soil

   ! no growth allocation [kg C / m2 /s]

   dnogrow(:)  = (1.0 - pgrow(:)) * 0.5 * dgpp(:)

   ! carbon stored in biomass

   dcveg(:) = dcveg(:) + (dnpp(:) - dlitter(:)) * zdelt * ncveg  !! FL

   ! carbon stored in soil

   dcsoil(:) = dcsoil(:) + (dlitter(:) - dres(:)) * zdelt * ncveg !! FL

   ! forest cover

   zforest(:) = atan((2.0*plai(:)*dcveg(:)-cveg_a) /cveg_b) /cveg_c +cveg_d
   zforest(:) = cveg_e * zforest(:)
   where (zforest(:) < 0.0) zforest(:) = 0.0
   where (zforest(:) > 1.0) zforest(:) = 1.0

   ! soil cover

   dvsoil(:) = atan((2.0*(1.0-plai(:))*dcveg(:)-cveg_a)/cveg_b)/cveg_c+cveg_d
   dvsoil(:) = cveg_e * dvsoil(:)
   where (dvsoil(:) < 0.0) dvsoil(:) = 0.0
   where (dvsoil(:) > 1.0) dvsoil(:) = 1.0

elsewhere
   dcveg(:)  = 0.0
   dcsoil(:) = 0.0
   dgppl(:)  = 0.0
end where

!     * discretization of vegetation state

!$TAF INCOMPLETE nbiomass,zforest,dvsoil
if (rnbiocats .ge. 2.0) then
  nbiomass(:)    = (zforest(:) * rnbiocats)
  zforest(:)       = MIN(1., real(nbiomass(:))/(rnbiocats-1.))
  nbiomass(:)    = (dvsoil(:) * rnbiocats)
  dvsoil(:)      = MIN(1., real(nbiomass(:))/(rnbiocats-1.))
end if

!     * derivation of land surface parameters

!$TAF INCOMPLETE dvveg,dvlai,dvz0,dvalb,zwmax,dvrhs
where (dls(:) > 0. .and. dglac(:) == 0.)

  dvveg(:)         = zveg(:)
  dvlai(:)         = -LOG(1. - dvveg(:))/cveg_k
  dvz0(:)          = pz0_max(:) * zforest(:)                        &
     &                   + vz0_min * (1. - zforest(:))
  dvalb(:)         = valb_min * dvveg(:)                          &
     &                   + valb_max * (1. - dvveg(:))
  zwmax(:)        = vwmax_max * dvsoil(:)                        &
     &                   + vwmax_min * (1. - dvsoil(:))
  dvrhs(:)         = pgs(:) * AMIN1(1., AMAX1(0., dwatc(:)        &
     &                   / (zwmax(:) * vws_crit) ))

end where

!     * modify albedo due to snow

!$TAF INCOMPLETE dvalb
where (dls(:) > 0. .and. dglac(:) == 0. .and. dsnow(:) > 0. )

!       * copy PUMA-2 formulation of snow albedo and mix it with forested
!         albedo

  zalbsn(:)        = (albsmax - albsmin) * (dt(:,NLEP)-263.16)    &
     &                   / (TMELT-263.16)
  zalbsn(:)        = MAX(albsmin,MIN(albsmax,albsmax-zalbsn(:)))
  dvalb(:)         = dvalb(:)                                     &
     &                   + (zalbsn(:) - dvalb(:)) * dsnow(:)            &
     &                   / (dsnow(:)+0.01)
  dvalb(:)         = vsalb_min * zforest(:)                         &
     &                   + dvalb(:) * (1.-zforest(:))
end where

!     * copy to puma/land variables

!
!     a) diagnostics only
!
 where (dls(:) > 0. .and. dglac(:) == 0.)
  dveg(:)          = dvveg(:)
  dlai(:)          = dvlai(:)
 elsewhere
  dveg(:)          = 0.
  dlai(:)          = 0.
 end where
!
!     b) interactive coupling
!
!$TAF INCOMPLETE dz0,dlai,dwmax,drhs,dalb,dforest
if(nbiome == 1) then

 where (dls(:) > 0. .and. dglac(:) == 0.)
  dz0(:)           = SQRT(dvz0(:)**2+dz0climo(:)**2)
  dlai(:)          = dvlai(:)
  dwmax(:)         = zwmax(:)
  drhs(:)          = dvrhs(:)
  dalb(:)          = dvalb(:)
  dforest(:)       = zforest(:)
 end where

end if

!     c) Send some arrays to GUI

call guihor("DCVEG"   // char(0),dcveg  ,1,1000.0,0.0)
call guihor("DVSOIL"  // char(0),dvsoil ,1,1000.0,0.0)
call guihor("ZFOREST" // char(0),zforest,1,1000.0,0.0)

return
end subroutine vegstep
