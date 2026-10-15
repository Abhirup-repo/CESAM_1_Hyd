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

      module pumamod
        use substitutes ! added, FastOpt
!     ****************************************************************
!     * The module resmod defines all resolution parameters          *
!     * NPRO_ATM, NLAT_ATM, NLEV_ATM                                 *
!     * If MPI (Message Passing Interface) is to be used it will     *
!     * also include a 'use mpi' statement                           *
!     ****************************************************************

      use resmod

!     ****************************************************************
!     * The number of processes for processing on parallel machines  *
!     * It must be a power of 2 e.g. (1, 2, 4, 8, ... NLAT)          *
!     ****************************************************************

      parameter(NPRO = NPRO_ATM)          ! Number of processes (CPUs)

!     ****************************************************************
!     * Set the horizontal resolution by specifying the number of    *
!     * latitudes. Use only values, listed below (FFT restrictions)  *
!     ****************************************************************

      parameter(NLAT = NLAT_ATM)                !  Number of latitudes

!     currently allowed values:  32,  48,  64, 128,  192,  256
!     corresponding truncation: T21, T31, T42, T85, T127, T170

!     ****************************************************************
!     * Set the vertical resolution, a minimum of 5 is recommended   *
!     ****************************************************************

      parameter(NLEV = NLEV_ATM)                    ! Number of levels

!     ****************************************************************
!     * Don't touch the following parameter definitions !            *
!     ****************************************************************

!     ********************
!     * Global Constants *
!     ********************
!
      parameter(NLON = NLAT + NLAT)        ! Number of longitudes
      parameter(NTRU = (NLON-1) / 3)       ! Triangular truncation
      parameter(NLPP = NLAT / NPRO)        ! Latitudes per process
      parameter(NHOR = NLON * NLPP)        ! Horizontal part
      parameter(NUGP = NLON * NLAT)        ! Number of gridpoints
      parameter(NPGP = NLON * NLAT / 2)    ! Dimension of packed fields
      parameter(NLEM = NLEV - 1)           ! Levels - 1
      parameter(NLEP = NLEV + 1)           ! Levels + 1
      parameter(NLSQ = NLEV * NLEV)        ! Levels squared
      parameter(NTP1 = NTRU + 1)           ! Truncation + 1
      parameter(NRSP =(NTRU+1)*(NTRU+2))   ! No of real global    modes
      parameter(NCSP = NRSP / 2)           ! No of complex global modes
      parameter(NSPP = (NRSP+NPRO-1)/NPRO) ! Modes per process
      parameter(NESP = NSPP * NPRO)        ! Dim of spectral fields
      parameter(NVCT = 2 * (NLEV+1))       ! Dim of Vert. Coord. Tab
      parameter(NZOM  = 2 * NTP1)          ! Dim for zonal mean diagnostics
      parameter(NROOT = 0)                 ! Master node

      parameter(EZ     = 1.63299310207D0)  ! ez = 1 / sqrt(3/8)
      parameter(PI     = 3.14159265359D0)  ! Pi
      parameter(TWOPI  = PI + PI)          ! 2 Pi
      parameter(RV     = 461.51)           ! Gas constant for water vapour
      parameter(ALV    = 2.5008E6)         ! Latent heat of vaporization
      parameter(ALS    = 2.8345E6)         ! Latent heat of sublimation
      parameter(ACPV   = 1870.)            ! Specific heat for water vapour
      parameter(TMELT  = 273.16)           ! Melting point (H2O)
      parameter(TMELT_CO2 = 148.0)         ! Melting point (CO2) - for Mars

!     ***************
!     * Date & Time *
!     ***************

      integer :: nstep           =       0 ! current timestep
      integer :: nstep0    ! start timestep, added, FastOpt
      integer :: nstep1    ! record counter for storage, added, FastOpt
      integer :: nstep2    ! record counter for storage, added, FastOpt
      integer :: nstepe    ! end timestep, added, FastOpt
      integer :: mocd            =       0 ! month countdown
      integer :: n_start_year    =    0001 ! start year
      integer :: n_reset_year    =    0001 ! start year
      integer :: n_start_month   =      01 ! start month
      integer :: n_start_step    =       0 ! start step since 1-Jan-0000
      integer :: n_days_per_month=      30 ! set to   0 for real calendar
      integer :: n_days_per_year =     360 ! set to 365 for real calendar
      integer :: n_run_years     =       0 ! years to run
      integer :: n_run_months    =       0 ! months to run
      integer :: n_run_days      =       0 ! days  to run (debugging)
      integer :: n_run_steps     =       0 ! steps to run (debugging)
      integer :: n_inner_loop    =      16 ! FastOpt: length of inner loop (for checkpointing)
      integer :: n_spinup        =       0 ! FastOpt: additional spinup timesteps
      integer :: mpstep          =       0 ! minutes / timestep = 1 day / ntspd
      integer :: ntspd           =      32 ! number of timesteps per day
      integer :: nwpd            =       1 ! number of writes per day
      integer :: ndatim(6)       =      -1 ! date & time array
      real    :: tmstart         =     0.0 ! start of run

      logical :: lunstable       = .false. ! added, FastOpt
      

!     **************************
!     * Global Integer Scalars *
!     **************************

      integer :: kick     =  1  ! add noise for kick > 0
      integer :: noutput  =  1  ! master switch for output: 0=no output
      integer :: nafter   = 30  ! write data interval
      integer :: ncost    = 30  ! Cost function evaluation interval; added, FastOpt
      integer :: ncoeff   =  0  ! number of modes to print
      integer :: ndiag    = 300 ! write diagnostics interval
      integer :: ndebug   =  0  ! activates additional printouts
      integer :: nexp     =  0  ! experiment number
      integer :: nexper   =  0  ! 1: dipole experiment
      integer :: ngui     =  0  ! 1: run with GUI
      integer :: nkits    =  3  ! number of initial timesteps
      integer :: nrestart =  0  ! 1 for true, 0 for false
      integer :: nrad     =  1  ! switches radiation off/on  1/0
      integer :: nflux    =  1  ! vertical diffusion 1/0
      integer :: neqsig   =  0  ! equidistant sigma levels (1/0)=(y/n)
      integer :: nprint   =  0  ! comprehensive print out (only for checks!)
      integer :: nprhor   =  0  ! grid point for print out (only for checks!)
      integer :: npacksp  =  1  ! pack spectral fields on output
      integer :: npackgp  =  1  ! pack gridpoint fields on output
      integer :: naccuout =  0  ! accumulation counter for diagnistics
      integer :: ndiaggp  =  0  ! switch for frank's gp-diagnostic arrays
      integer :: ndiagsp  =  0  ! switch for frank's sp-diagnostic arrays
      integer :: ndiagcf  =  0  ! switch for cloud forcing diagnostic
      integer :: ndiaggp2d=  0  ! number of additional 2-d gp-diagnostic arrays
      integer :: ndiaggp3d=  0  ! number of additional 3-d gp-diagnostic arrays
      integer :: ndiagsp2d=  0  ! number of additional 2-d sp-diagnostic arrays
      integer :: ndiagsp3d=  0  ! number of additional 3-d sp-diagnostic arrays
      integer :: nhdiff   = 15  ! critical wavenumber for horizontal diffusion
      integer :: ntime    =  0  ! switch for time use diagnostics
      integer :: nperpetual = 0 ! radiation day for perpetual integration
      integer :: n_sea_points = 0 ! number of sea points on grid

!     ***********************
!     * Global Real Scalars *
!     ***********************

      real :: dawn  =     0.0 ! angle threshold for solar radiation
      real :: delt            ! 2 pi / ntspd timestep interval
      real :: delt2           ! 2 * delt
      real :: dtep  =     0.0
      real :: dtns  =     0.0
      real :: dtrop = 12000.0
      real :: dttrp =     2.0
      real :: tgr   =   288.0  ! Temperature ground in mean profile
      real :: psurf =101100.0  ! global mean surface pressure
      real :: time0 =     0.0  ! start time (for performance estimates)
      real :: co2   =   360.0  ! atm. co2 concentration (ppmv)
      real :: umax  =     0.0  ! diagnostic U max
      real :: t2mean=     0.0  ! diagnostic T2m mean
      real :: precip=     0.0  ! diagnostic precipitation mean
      real :: evap  =     0.0  ! diagnostic evaporation
      real :: olr   =     0.0  ! Outgoing longwave radiation

!     **************************
!     * Global Spectral Arrays *
!     **************************
      real :: dstore1(nhor)=0.
      real :: dstore2(nhor)=0.
      real :: dstore3(nhor)=0.
      real ::  sd(NESP,NLEV) = 0.0 ! Spectral Divergence
      real ::  st(NESP,NLEV) = 0.0 ! Spectral Temperature
      real ::  sz(NESP,NLEV) = 0.0 ! Spectral Vorticity
      real ::  sq(NESP,NLEV) = 0.0 ! Spectral Specific Humidity
      real ::  sp(NESP)      = 0.0 ! Spectral Pressure (ln Ps)
      real ::  so(NESP)      = 0.0 ! Spectral Orography
      real ::  sr(NESP,NLEV) = 0.0 ! Spectral Restoration Temperature

      real :: sdp(NSPP,NLEV) = 0.0 ! Spectral Divergence  Partial
      real :: stp(NSPP,NLEV) = 0.0 ! Spectral Temperature Partial
      real :: szp(NSPP,NLEV) = 0.0 ! Spectral Vorticity   Partial
      real :: sqp(NSPP,NLEV) = 0.0 ! Spectral S.Humidity  Partial
      real :: spp(NSPP)      = 0.0 ! Spectral Pressure    Partial
      real :: sop(NSPP)      = 0.0 ! Spectral Orography   Partial
      real :: srp(NSPP,NLEV) = 0.0 ! Spectral Restoration Partial

      real :: sdt(NSPP,NLEV) = 0.0 ! Spectral Divergence  Tendency
      real :: stt(NSPP,NLEV) = 0.0 ! Spectral Temperature Tendency
      real :: szt(NSPP,NLEV) = 0.0 ! Spectral Vorticity   Tendency
      real :: sqt(NSPP,NLEV) = 0.0 ! Spectral S.Humidity  Tendency
      real :: spt(NSPP)      = 0.0 ! Spectral Pressure    Tendency

      real :: sdm(NSPP,NLEV) = 0.0 ! Spectral Divergence  Minus
      real :: stm(NSPP,NLEV) = 0.0 ! Spectral Temperature Minus
      real :: szm(NSPP,NLEV) = 0.0 ! Spectral Vorticity   Minus
      real :: sqm(NSPP,NLEV) = 0.0 ! Spectral S.Humidity  Minus
      real :: spm(NSPP)      = 0.0 ! Spectral Pressure    Minus

      real :: sak(NESP,NLEV)   = 0.0 ! horizontal diffusion
      real :: sakpp(NSPP,NLEV) = 0.0 ! horizontal diffusion partial
      real :: sqout(NESP,NLEV) = 0.0 ! specific humidity for output
      real :: spnorm(NESP)     = 0.0 ! Factors for output normalization

      integer :: nindex(NESP) = NTRU ! Holds wavenumber
      integer :: nscatsp(NPRO)= NSPP ! Used for reduce_scatter op
      integer :: ndel(NLEV)   =    2 ! ndel for horizontal diffusion

!     ************************************************
!     * Global Gridpoint Arrays (un-dimensionalized) *
!     ************************************************

      real :: gd(NHOR,NLEV)   = 0. ! divergence
      real :: gt(NHOR,NLEV)   = 0. ! temperature (-t0)
      real :: gz(NHOR,NLEV)   = 0. ! absolut vorticity
      real :: gq(NHOR,NLEV)   = 0. ! spec. humidity
      real :: gu(NHOR,NLEV)   = 0. ! zonal wind (*cos(phi))
      real :: gv(NHOR,NLEV)   = 0. ! meridional wind (*cos(phi))
      real :: gtdt(NHOR,NLEV) = 0. ! t-tendency
      real :: gqdt(NHOR,NLEV) = 0. ! q-tendency
      real :: gudt(NHOR,NLEV) = 0. ! u-tendency
      real :: gvdt(NHOR,NLEV) = 0. ! v-tendency
      real :: gp(NHOR)        = 0. ! surface pressure or ln(ps)
      real :: gpj(NHOR)       = 0. ! dln(ps)/dphi

      real :: rcsq(NHOR)      = 0. ! 1/cos(phi)**2

!     *********************************************
!     * Global Gridpoint Arrays (dimensionalized) *
!     *********************************************

      real :: dt(NHOR,NLEP)   = 0.     ! temperature
      real :: dq(NHOR,NLEP)   = 0.     ! spec. humidity
      real :: dqraw(NHOR,NLEP)= 0.     ! added, FastOpt: spec. humidity, may contain negative values
      real :: du(NHOR,NLEP)   = 0.     ! zonal wind (*cos(phi))
      real :: dv(NHOR,NLEP)   = 0.     ! meridional wind (*cos(phi))
      real :: dp(NHOR)        = 0.     ! surface pressure
      real :: dqsat(NHOR,NLEP)= 0.     ! saturation humidity
      real :: dqt(NHOR,NLEP)  = 0.     ! adiabatic q-tendencies (for eg kuo)
      real :: dcc(NHOR,NLEP)  = 0.     ! cloud cover
      real :: dql(NHOR,NLEP)  = 0.     ! Liquid water content
      real :: dw(NHOR,NLEV)   = 0.     ! vertical velocity (dp/dt)
      real :: dtdt(NHOR,NLEP) = 0.     ! t-tendency
      real :: dqdt(NHOR,NLEP) = 0.     ! q-tendency
      real :: dudt(NHOR,NLEP) = 0.     ! u-tendency
      real :: dvdt(NHOR,NLEP) = 0.     ! v-tendency

!     *************
!     * Radiation *
!     *************

      real :: dalb(NHOR)               ! albedo
      real :: dswfl(NHOR,NLEP)         ! net solar radiation
      real :: dlwfl(NHOR,NLEP)         ! net thermal radiation
      real :: dflux(NHOR,NLEP)         ! net radiation (SW + LW)
      real :: dfu(NHOR,NLEP)           ! solar radiation upward
      real :: dfd(NHOR,NLEP)           ! solar radiation downward
      real :: dftu(NHOR,NLEP)          ! thermal radiation upward
      real :: dftd(NHOR,NLEP)          ! thermal radiation downward
! define 2d parameter matrix
      real :: th2ocg(NUGP)=0.
      real :: acllwrg(NUGP)=0.
      real :: tswr1g(NUGP)=0.
      real :: tswr2g(NUGP)=0.
      real :: tswr3g(NUGP)=0.
      real :: tpofmtg(NUGP)=0.
      real :: th2ocm(NHOR)=0.
      real :: acllwrm(NHOR)=0. 
      real :: tswr1m(NHOR)=0.
      real :: tswr2m(NHOR)=0.
      real :: tswr3m(NHOR)=0.
      real :: tpofmtm(NHOR)=0.
      real :: rcritg(NUGP,NLEV)=0.  ! guokun included for 2d par
      real :: rcritm(NHOR,NLEV)=0.
      real :: albgmaxg(NUGP)=0.  ! guokun included for 2d par
      real :: albgmaxm(NHOR)=0.
      real :: alpha1g(NUGP)=0.
      real :: beta1g(NUGP)=0.
      real :: gamma1g(NUGP)=0.
      real :: alpha1m(NHOR)=0.
      real :: beta1m(NHOR)=0.
      real :: gamma1m(NHOR)=0.

! New parameter for land 
      real    :: ss_field_gp (NUGP) =0.
      real    :: ss_field_hor (NHOR) =0. 

      real    :: ss_flux_gp (NUGP) =0.
      real    :: ss_flux_hor (NHOR) =0. 



!     ***********
!     * SURFACE *
!     ***********

      real :: drhs(NHOR)  = 0.  ! surface wetness
      real :: dls(NHOR)   = 1.  ! land(1)/sea(0) mask
      real :: dz0(NHOR)   = 0.  ! rougthness length
      real :: diced(NHOR) = 0.  ! ice thickness
      real :: dicec(NHOR) = 0.  ! ice cover
      real :: dtaux(NHOR) = 0.  ! x-surface wind stress
      real :: dtauy(NHOR) = 0.  ! y-surface wind stress
      real :: dust3(NHOR) = 0.  ! u-star**3 (needed eg. for coupling)
      real :: dshfl(NHOR) = 0.  ! surface sensible heat flx
      real :: dlhfl(NHOR) = 0.  ! surface latent heat flx
      real :: devap(NHOR) = 0.  ! surface evaporation
      real :: dtsa(NHOR)  = 0.  ! surface air temperature
      real :: dmld(NHOR)  = 0.  ! mixed-layer depth (output from ocean)
!
      real dshdt(NHOR),dlhdt(NHOR)

!     *********
!     * WATER *
!     *********

      real dprc(NHOR)       ! Convective Precip   (m/s)
      real dprl(NHOR)       ! Large Scale Precip  (m/s)
      real dprs(NHOR)       ! Snow Fall           (m/s)
      real dqvi(NHOR)       ! vertical integrated specific humidity (kg/m**2)
      real :: dqviraw(NHOR)= 0.     ! added, FastOpt: vertically integrated spec. humidity, may contain negative values

!     *********
!     * BIOME *
!     *********

      real :: dveg(NHOR)    = 0.5  ! vegetation cover (fract.)
      real :: dforest(NHOR) = 0.5  ! forest cover (fract.)
      real :: dlai(NHOR)    = 0.0  ! leaf area index (integer value)
      real :: dwmax(NHOR)   = 0.0  ! field capacity (m)
      real :: dnpp(NHOR)    = 0.0  ! net primary production (kg C /m2 /s)
      real :: dgpp(NHOR)    = 0.0  ! gross primary production (kg C /m2 /s)
      real :: dgppl(NHOR)   = 0.0  ! light limited gpp (kg C /m2 /s)
      real :: dgppw(NHOR)   = 0.0  ! water limited gpp (kg C /m2 /s)
      real :: dres(NHOR)    = 0.0  ! heterotrophic respiration (kg C /m2 /s)
      real :: dlitter(NHOR) = 0.0  ! litterfall (kg C /m2 /s)
      real :: dnogrow(NHOR) = 0.0  ! no growth allocation (kg C /m2 /s)

      real :: dcveg(NHOR)    =0.0  ! organic carbon stored in biomass
      real :: dcsoil(NHOR)   =0.0  ! organic carbon stored in soil
      real :: dwloss(NHOR)   =0.0  ! water loss

!     ********
!     * SOIL *
!     ********

      real :: dwatc(NHOR)   = 0.  ! soil wetness (m)
      real :: drunoff(NHOR) = 0.  ! surface runoff (m/s)
      real :: dsnow(NHOR)   = 0.  ! snow depth (m)
      real :: dsmelt(NHOR)  = 0.  ! snow melt (m/s water eq.)
      real :: dsndch(NHOR)  = 0.  ! snow depth change (m/s water eq.)
      real :: dtsoil(NHOR)  = 0.  ! soil temperature uppermost level (K)
      real :: dtd2(NHOR)    = 0.  ! soil temperature level 2 (K)
      real :: dtd3(NHOR)    = 0.  ! soil temperature level 3 (K)
      real :: dtd4(NHOR)    = 0.  ! soil temperature level 4 (K)
      real :: dtd5(NHOR)    = 0.  ! soil temperature lowermost level (K)
      real :: dglac(NHOR)   = 0.  ! glacier mask (0.,1.)

!     *********************
!     * Diagnostic Arrays *
!     *********************

      integer :: ndl(NLEV) = 0

      real :: csu(NLAT,NLEV),csv(NLAT,NLEV),cst(NLAT,NLEV),csm(NLAT,NLEV)
      real :: ccc(NLAT,NLEV)
      real :: span(NESP)

      real, allocatable :: dgp2d(:,:),dsp2d(:,:)     ! 2-d diagnostics
      real, allocatable :: dgp3d(:,:,:),dsp3d(:,:,:) ! 3-d diagnostics
      real, allocatable :: dclforc(:,:)  ! cloud forcing diagnostics

!
!     accumulated output
!

      real :: aevap(NHOR) = 0. ! acculumated evaporation
      real :: aprl(NHOR)  = 0. ! acculumated lage scale precip.
      real :: aprc(NHOR)  = 0. ! acculumated convective precip.
      real :: aprs(NHOR)  = 0. ! acculumated snow fall
      real :: ashfl(NHOR) = 0. ! acculumated sensible heat flux
      real :: alhfl(NHOR) = 0. ! acculumated latent heat flux
      real :: aroff(NHOR) = 0. ! acculumated surface runoff
      real :: asmelt(NHOR)= 0. ! acculumated snow melt
      real :: asndch(NHOR)= 0. ! acculumated snow depth change
      real :: acc(NHOR)   = 0. ! acculumated total cloud cover
      real :: assol(NHOR) = 0. ! acculumated surface solar radiation
      real :: asthr(NHOR) = 0. ! acculumated surface thermal radiation
      real :: atsol(NHOR) = 0. ! acculumated top solar radiation
      real :: atthr(NHOR) = 0. ! acculumated top thermal radiation
      real :: assolu(NHOR)= 0. ! acculumated surface solar radiation upward
      real :: asthru(NHOR)= 0. ! acculumated surface thermal radiation upward
      real :: atsolu(NHOR)= 0. ! acculumated top solar radiation upward
      real :: ataux(NHOR) = 0. ! acculumated zonal wind stress
      real :: atauy(NHOR) = 0. ! acculumated meridional wind stress
      real :: aqvi(NHOR)  = 0. ! acculumated vertical integrated q
      real :: atsa(NHOR)  = 0. ! accumulated surface air temperature
      real :: ats0(NHOR)  = 0. ! accumulated surface temperature

      real :: anpp(NHOR)    = 0.0
      real :: agpp(NHOR)    = 0.0
      real :: agppl(NHOR)   = 0.0
      real :: agppw(NHOR)   = 0.0
      real :: anogrow(NHOR) = 0.0
      real :: aresh(NHOR)   = 0.0
      real :: alitter(NHOR) = 0.0
      real :: awloss(NHOR)  = 0.0

!     *******************
!     * Latitude Arrays *
!     *******************

      character(len=3) chlat(NLAT)
      real :: si(NLAT)        ! sin(phi)
      real :: gw(NLAT)        ! Gaussian weights
      real :: gaussw(NLAT)    ! Gaussian weights, FastOpt
      real :: csq(NLAT)       ! cos(phi)**2
      real :: cola(NLAT)      ! cos(phi)
      real :: rcs(NLAT) = 0.0 ! 1 / cos(phi)

!     ****************
!     * Level Arrays *
!     ****************

      real :: tdissd(NLEV)  =  0.20 ! diffusion time scale for divergence [days]
      real :: tdissz(NLEV)  =  1.10 ! diffusion time scale for vorticity [days]
      real :: tdisst(NLEV)  =  5.60 ! diffusion time scale for temperature [days]
      real :: tdissq(NLEV)  =  5.60 ! diffusion time scale for sp. humidity [days]
      real :: restim(NLEV) =  00.0
      real :: t0(NLEV)     = 250.0
      real :: tfrc(NLEV)   = (/ (0.0,i=1,NLEM), 0.0 /)
      real :: sigh(NLEV)   = (/(0.0,i=1,NLEV)/)
      real damp(NLEV)
      real dsigma(NLEV)
      real rdsig(NLEV)
      real sigma(NLEV)
      real sigmah(NLEV)
      real t01s2(NLEV)
      real tkp(NLEV)
      real c(NLEV,NLEV)
      real g(NLEV,NLEV)
      real tau(NLEV,NLEV)
      real bm1(NLEV,NLEV,NTRU)
! guokun added nudging coefficient
      real :: dnudg(NLEV)    =  0. ! nudging time scale for divergence  [days]
      real :: tnudg(NLEV)    =  0. ! nidging time scale for temperature [days
      real :: znudg(NLEV)    =  0. ! nudging time scale for vorticity [days
      real :: qnudg(NLEV)    =  0. ! nudging time scale for sp huminidity [days
      real :: pnudg          =  0. ! nudging time scale for surface pressure [days
      real :: smnudg         =  0. ! nudging time scale for soil moisture [days
      
      integer :: makenudge   =  0 ! switch to on/off atmospheric nudging !! A.B 27.04.2026
      
!     ******************
!     * Parallel Stuff *
!     ******************

      integer :: mpinfo  = 0
      integer :: mypid   = 0
      integer :: myworld = 0
      integer :: nproc   = NPRO
      character (80),allocatable :: ympname(:)
!!      character (len=40) :: ympname(NPRO)

!     **********************************************
!     * version identifier (to be set in puma.f90) *
!     * will be printed in puma sub *prolog*       *
!     **********************************************

      character(len=80) :: pumaversion

!     ******************************
!     * Planet dependent variables *
!     ******************************

      parameter(AKAP_EARTH   = 0.286 )      ! Kappa Earth
      parameter(AKAP_MARS    = 0.2273)      ! Kappa Mars
      parameter(ALR_EARTH    = 0.0065)      ! Lapse rate Earth
      parameter(ALR_MARS     = 0.0025)      ! Lapse rate Mars
      parameter(GA_EARTH     = 9.81)        ! Gravity Earth
      parameter(GA_MARS      = 3.74)        ! Gravity Mars
      parameter(GASCON_EARTH = 287.0)       ! Gas constant for dry air on Earth
      parameter(GASCON_MARS  = 188.9)       ! Gas constant for dry air on Mars 
      parameter(PLARAD_EARTH = 6371000.0)   ! Earth radius
      parameter(PLARAD_MARS  = 3397000.0)   ! Mars radius
      parameter(PNU_EARTH    = 0.10)        ! Time filter Earth
      parameter(PNU_MARS     = 0.15)        ! Time filter Mars
      parameter(SID_DAY_EARTH= 86164.)      ! Siderial day Earth 23h 56m 04s
      parameter(SID_DAY_MARS = 88642.)      ! Siderial day Mars  24h 37m 22s
      parameter(SOL_DAY_EARTH= 86400.)      ! Solar day Earth    24h 00m 00s
      parameter(SOL_DAY_MARS = 88775.)      ! Solar day Mars     24h 39m 35s
      parameter(WW_EARTH  = TWOPI / SID_DAY_EARTH) ! Scaling parameter
      parameter(WW_MARS   = TWOPI / SID_DAY_MARS ) ! Scaling parameter
      parameter(RA1_LIQUID   = 610.78)      ! Parameter for Magnus-Teten-Formula
      parameter(RA2_LIQUID   =  17.2693882) ! for saturation vapor pressure
      parameter(RA4_LIQUID   =  35.86)      ! over liquid water (used on Earth)
      parameter(RA1_ICE      = 610.66)      ! Parameter for Magnus-Teten-Formula
      parameter(RA2_ICE      =  21.875)     ! for saturation vapor pressure
      parameter(RA4_ICE      =   7.65)      ! over ice (used on Mars)
      parameter(SOLCON_EARTH = 1365.0)      ! Solar constant [W/m2]
      parameter(SOLCON_MARS  =  595.0)      ! Solar constant [W/m2]
      parameter(ECCEN_EARTH  = 0.016715)    ! Eccentricity (AMIP-II value)
      parameter(ECCEN_MARS   = 0.09341233)  ! Eccentricity Mars
      parameter(OBLIQ_EARTH  = 23.441)      ! Obliquity [deg] (AMIP-II)
      parameter(OBLIQ_MARS   = 25.19)       ! Obliquity [deg] Mars
! FastOpt      parameter(LOPER_EARTH  = 102.7)       ! Longitude of perihelion Earth
      real, parameter :: LOPER_EARTH  = 102.7       !  changed, FastOpt
! FastOpt      parameter(LOPER_MARS   = 336.04084)   ! Longitude of perihelion Mars
      real, parameter :: LOPER_MARS   = 336.04084   ! changed, FastOpt

      character(len=32) :: yplanet = "Earth"

      integer :: mars = 0              ! Global switch for planet Mars

      real :: akap   = AKAP_EARTH      ! Kappa
      real :: alr    = ALR_EARTH       ! Lapse rate
      real :: ga     = GA_EARTH        ! Gravity
      real :: gascon = GASCON_EARTH    ! Gas constant for dry air
      real :: plarad = PLARAD_EARTH    ! Planet radius
      real :: pnu    = PNU_EARTH       ! Time filter
      real :: siderial_day = SID_DAY_EARTH ! Length of siderial day [sec]
      real :: solar_day    = SOL_DAY_EARTH ! Length of solar day [sec]
      real :: ww     = WW_EARTH        ! Scaling using 2 Pi / siderial day
      real :: ra1    = RA1_LIQUID      !
      real :: ra2    = RA2_LIQUID      !
      real :: ra4    = RA4_LIQUID      !
      real :: acpd   = 0.0             ! acpd = gascon / akap ! Specific heat for dry air
      real :: adv    = 0.0             ! acpv / acpd - 1.0
      real :: cv     = 0.0             ! cv = plarad * ww
      real :: ct     = 0.0             ! ct = CV * CV / gascon
      real :: pnu21  = 0.0             ! pnu21 = 1.0 - 2.0 * pnu ! Time filter 2
      real :: rdbrv  = 0.0             ! rdbrv = gascon / RV  ! Rd/Rv

!     ******************************************
!     * GUI (Graphical User Interface for X11) *
!     ******************************************

      integer, parameter :: PUMA   = 0
      integer, parameter :: PLASIM = 2
      parameter (NPARCS = 4)          ! Number of GUI parameters
      character(6) :: yguinam(NPARCS) ! Variable names for GUI display
      integer :: nshutdown = 0        ! Flag for shutdown request
      integer :: model     = PLASIM
      integer :: ntpal     = 2        ! Color pallette for temperature
      real(kind=4) :: parc(NPARCS)            ! Values of GUI parameters
      real(kind=4) :: crap(NPARCS)            ! Backup of parc(NPARCS)
      real(kind=4) :: guimin(NPARCS)          ! lower limit
      real(kind=4) :: guimax(NPARCS)          ! upper limit
      real(kind=4) :: guiinc(NPARCS)          ! increment
      real, allocatable :: sr1(:,:)   ! dummy array for PUMA guimod compatibility
      real, allocatable :: sr2(:,:)   ! dummy array for PUMA guimod compatibility
      logical :: ldisp   = .FALSE.    ! DISP changed by GUI
      logical :: ldtep   = .FALSE.    ! DTEP changed by GUI
      logical :: ldtns   = .FALSE.    ! DTNS changed by GUI
      logical :: lrotspd = .FALSE.    ! rotspd changed by GUI

      end module pumamod
