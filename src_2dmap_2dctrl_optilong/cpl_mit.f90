      module cplmod
        implicit none
!
        integer, parameter :: nimax = 100  ! max. no points for interpol. 
        integer, parameter :: nxa = 64
        integer, parameter :: nya = 32
        integer, parameter :: nxo = 90
        integer, parameter :: nyo = 52
        integer, parameter :: radea = 6.371E6
!
      character(len=30) :: version='cplmod 11/09/12'
!
!     define rot process (for MPI)
!
      integer, parameter :: NROOT = 0
!
      real :: xa(nxa,nya)        ! atm. longitudes
      real :: xa1(nxa,nya)       ! atm. longitudes western boundary 
      real :: xa2(nxa,nya)       ! atm. longitudes eastern boundary
      real :: ya(nxa,nya)        ! atm. latitudes 
      real :: ya1(nxa,nya)       ! atm. latitudes boundary 1
      real :: ya2(nxa,nya)       ! atm. latitudes boundary 2
      real :: xos(nxo,nyo)       ! oce. longitudes (tracer)
      real :: xos1(nxo,nyo)      ! oce. longitudes western boundary (scalar)
      real :: xos2(nxo,nyo)      ! oce. longitudes eastern boundary (scalar)
      real :: xou(nxo,nyo)       ! oce. longitudes (u-velocity)
      real :: xou1(nxo,nyo)      ! oce. longitudes western boundary (u)
      real :: xou2(nxo,nyo)      ! oce. longitudes eastern boundary (u)
      real :: xov(nxo,nyo)       ! oce. longitudes (v-velocity)
      real :: xov1(nxo,nyo)      ! oce. longitudes western boundary (v)
      real :: xov2(nxo,nyo)      ! oce. longitudes eastern boundary (v)
      real :: yos(nxo,nyo)       ! oce. latitudes (tracer)
      real :: yos1(nxo,nyo)      ! oce. latitudes boundary 1 (tracer)
      real :: yos2(nxo,nyo)      ! oce. latitudes boundary 2 (tracer)
      real :: you(nxo,nyo)       ! oce. latitudes (u-velocity)
      real :: you1(nxo,nyo)      ! oce. latitudes boundary 1 (u)
      real :: you2(nxo,nyo)      ! oce. latitudes boundary 2 (u)
      real :: yov(nxo,nyo)       ! oce. latitudes (v-velocity)
      real :: yov1(nxo,nyo)      ! oce. latitudes boundary 1 (v)
      real :: yov2(nxo,nyo)      ! oce. latitudes boundary 2 (v)

!
      real :: agw(nxa,nya)       ! atm. area weights
      real :: ogws(nxo,nyo)      ! oce. area weights (tracer)
      real :: ogwu(nxo,nyo)      ! oce. area weights (u-velocity)
      real :: ogwv(nxo,nyo)      ! oce. area weights (v-velocity)
!
      real :: aslm(nxa,nya)   = 0. ! atm. land(0)-see(1)-mask
      real :: oslms(nxo,nyo)  = 0. ! oce. land(0)-see(1)-mask (scalar)
      real :: oslmu(nxo,nyo)  = 0. ! oce. land(0)-see(1)-mask (u-velocity)
      real :: oslmv(nxo,nyo)  = 0. ! oce. land(0)-see(1)-mask (v-velocity)

      real :: asst(nxa,nya)   = 0. ! atm. sst
      real :: ats(nxa,nya)    = 0. ! atm. tsurf
      real :: aicec(nxa,nya)  = 0. ! atm. ice concentration
      real :: aiced(nxa,nya)  = 0. ! atm. ice thickness
      real :: apme(nxa,nya)   = 0. ! atm. fresh water flux
      real :: aroff(nxa,nya)  = 0. ! atm. runoff
      real :: ataux(nxa,nya)  = 0. ! atm. u-stress
      real :: atauy(nxa,nya)  = 0. ! atm. v-stress
      real :: ahfl(nxa,nya)   = 0. ! atm. heat flux
      real :: aswr(nxa,nya)   = 0. ! atm. solar radiation
      real :: aps(nxa,nya)    = 0. ! atm. surf. pressure
      real :: aim(nxa,nya)    = 0. ! atm. (ice+snow) mass 
!
      real :: osst(nxo,nyo)   = 0. ! oce. sst
      real :: osss(nxo,nyo)   = 0. ! oce. sss
      real :: ots(nxo,nyo)    = 0. ! oce. tsurf
      real :: oicec(nxo,nyo)  = 0. ! oce. ice concentration
      real :: oiced(nxo,nyo)  = 0. ! oce. ice thickness
      real :: opme(nxo,nyo)   = 0. ! oce. fresh water flux
      real :: otaux(nxo,nyo)  = 0. ! oce. u-stress
      real :: otauy(nxo,nyo)  = 0. ! oce. v-stress
      real :: ohfl(nxo,nyo)   = 0. ! oce. heat flux
      real :: oswr(nxo,nyo)   = 0. ! oce. solar radiation
      real :: ops(nxo,nyo)    = 0. ! oce. surf. pressure
      real :: oim(nxo,nyo)    = 0. ! oce. (ice+snow) mass
!
      real :: assta(nxa,nya)  = 0. ! atm. sst (accumulated)
      real :: asss(nxa,nya)   = 0. ! atm. sss (accumulated)
      real :: atsa(nxa,nya)   = 0. ! atm. tsurf (accumulated)
      real :: aiceca(nxa,nya) = 0. ! atm. ice concentration (acumulated)
      real :: aiceda(nxa,nya) = 0. ! atm. ice thickness (accumulated)
      real :: apmea(nxa,nya)  = 0. ! atm. frsh water flux (accumulated)
      real :: aroffa(nxa,nya) = 0. ! atm. runoff (accumulated)
      real :: atauxa(nxa,nya) = 0. ! atm. u-stress (accumulated)
      real :: atauya(nxa,nya) = 0. ! atm. v-stress (accumulated)
      real :: ahfla(nxa,nya)  = 0. ! atm. heat flux (accumulated)
      real :: aswra(nxa,nya)  = 0. ! atm. solar radiation (accumulated)
      real :: apsa(nxa,nya)   = 0. ! atm. surf. pressure (accumulated)
      real :: aima(nxa,nya)   = 0. ! atm. (ice+snow) mass (accumulated)
!
      real :: ossta(nxo,nyo)  = 0. ! oce. sst (accumulated)
      real :: osssa(nxo,nyo)  = 0. ! oce. sss (accumulated)
      real :: otsa(nxo,nyo)   = 0. ! oce. tsurf (accumulated)
      real :: oiceca(nxo,nyo) = 0. ! oce. ice concentration (accumulated)
      real :: oiceda(nxo,nyo) = 0. ! oce. ice thickness(accumulated)
      real :: opmea(nxo,nyo)  = 0. ! oce. fresh water flux (accumulated)
      real :: otauxa(nxo,nyo) = 0. ! oce. u-stress (accumulated)
      real :: otauya(nxo,nyo) = 0. ! oce. v-stress (accumulated)
      real :: ohfla(nxo,nyo)  = 0. ! oce. heat flux (accumulated)
      real :: oswra(nxo,nyo)  = 0. ! oce. solar radiation (accumulated)
      real :: opsa(nxo,nyo)   = 0. ! oce. surface pressure (accumulated)
      real :: oima(nxo,nyo)   = 0. ! oce. (ice+snow) mass (accumulated)
!
!     indices and weights for bi-linear interpolation
!
      integer :: ios2a(nxa,nya,2) = 0 ! ocean scalar -> atm.
      integer :: ia2os(nxo,nyo,2) = 0 ! atm. -> ocean scalar
      integer :: ia2ou(nxo,nyo,2) = 0 ! atm. -> ocean u
      integer :: ia2ov(nxo,nyo,2) = 0 ! atm. -> ocean v
      real :: wos2a(nxa,nya,4) = 0.   ! ocean scalar -> atm. 
      real :: wa2os(nxo,nyo,4) = 0.   ! atm. -> ocean scalar 
      real :: wa2ou(nxo,nyo,4) = 0.   ! atm. -> ocean u 
      real :: wa2ov(nxo,nyo,4) = 0.   ! atm. -> ocean v 
!
!     indices and weights for area overlap interpolation
!
      integer :: inos2a(nxa,nya) = 0          ! no ocean scalar -> atm.
      integer :: ina2os(nxo,nyo) = 0          ! no atm. -> ocean scalar
      integer :: ina2ou(nxo,nyo) = 0          ! no atm. -> ocean u
      integer :: ina2ov(nxo,nyo) = 0          ! no atm. -> ocean v
      integer :: ios2a_a(nxa,nya,nimax,2) = 0 ! ocean scalar -> atm.
      integer :: ia2os_a(nxo,nyo,nimax,2) = 0 ! atm. -> ocean scalar
      integer :: ia2ou_a(nxo,nyo,nimax,2) = 0 ! atm. -> ocean u
      integer :: ia2ov_a(nxo,nyo,nimax,2) = 0 ! atm. -> ocean v
      real :: wos2a_a(nxa,nya,nimax) = 0.     ! ocean scalar -> atm.
      real :: wa2os_a(nxo,nyo,nimax) = 0.     ! atm. -> ocean scalar
      real :: wa2ou_a(nxo,nyo,nimax) = 0.     ! atm. -> ocean u
      real :: wa2ov_a(nxo,nyo,nimax) = 0.     ! atm. -> ocean v
!
      integer :: iroffi(nxa,nya) = 0  ! runoff mask input (nroffc=2,3)
      integer :: iroffo(nxa,nya) = 0  ! runoff mask output (nroffc=3)
!  
      integer :: ncoupling = 2 ! switch for coupling method 
      integer :: nroffc    = 0 ! switch for runoff correction
!
      integer :: nssta = 0 ! counter for atm. sst accumulation
      integer :: ntsa  = 0 ! counter for atm. tsurf accumulation
      integer :: nicea = 0 ! counter for atm. ice accumulation
      integer :: npmea = 0 ! counter for atm. fresh water flux accumulation
      integer :: nroffa= 0 ! counter for atm. runoff
      integer :: ntaua = 0 ! counter for atm. stress accumulation
      integer :: nhfla = 0 ! counter for atm. heat flux accumulation
      integer :: nswra = 0 ! counter for atm. solar radiation
      integer :: npsa  = 0 ! counter for atm. surf. pressure
      integer :: nima  = 0 ! counter for atm. (ice+snow) mass 
      integer :: nssto = 0 ! counter for oce. sst accumulation
      integer :: nssso = 0 ! counter for oce. sss accumulation
      integer :: ntso  = 0 ! counter for oce. tsurf accumulation
      integer :: niceo = 0 ! counter for oce. ice accumulation
      integer :: npmeo = 0 ! counter for oce. fresh water flux accumulation
      integer :: ntauo = 0 ! counter for oce. stress accumulation 
      integer :: nhflo = 0 ! counter for oce. heat flux accumulation
      integer :: nswro = 0 ! counter for oce. solar radiation
      integer :: npso  = 0 ! counter for oce. surface pressure
      integer :: nimo  = 0 ! counter for oce. (ice+snow) mass 
!
      integer :: naomod = 0 ! atm-oce coupling interval
!
      integer :: ndatim(6) = 0 ! array containing time/date information 
      integer :: nyear     = 0 ! actual year
      integer :: nmonth    = 0 ! actual month
      integer :: nday      = 0 ! actual day
      integer :: nhour     = 0 ! actual hour
      integer :: nmin      = 0 ! actual minute
      integer :: nwday     = 0 ! actual day of the week
!
!     namelist parameter
!
      integer :: nprint    = 0 ! extended diagnostics
      integer :: ndbox     = 1 ! ocean diagnostic gp (i)
      integer :: ndboy     = 1 ! ocean diagnostic gp (j)
      integer :: nissta    = 0 ! choise for ssta interpolation
      integer :: nitsa     = 0 ! choise for tsa interpolation
      integer :: niicea    = 0 ! choise for icea interpolation
      integer :: nipmea    = 0 ! choise for pmea interpolation
      integer :: niroffa   = 0 ! choise for roff interpolation
      integer :: nitaua    = 0 ! choise for taua interpolation
      integer :: nihfla    = 0 ! choise for hfla interpolation
      integer :: niswra    = 0 ! choise for swra interpolation
      integer :: nipsa     = 0 ! choise for psa  interpolation
      integer :: niima     = 0 ! choise for ima  interpolation
      integer :: nissto    = 0 ! choise for ssto interpolation
      integer :: nissso    = 0 ! choise for ssso interpolation
      integer :: nitso     = 0 ! choise for tso interpolation
      integer :: niiceo    = 0 ! choise for iceo interpolation
      integer :: nihflo    = 1 ! choise for hflo interpolation
      integer :: niswro    = 1 ! choise for swro interpolation
      integer :: nipso     = 1 ! choise for pso  interpolation
      integer :: niimo     = 1 ! choise for imo  interpolation
      integer :: ncssta    = 3 ! choise for ssta correction
      integer :: nctsa     = 3 ! choise for tsa correction
      integer :: ncicea    = 2 ! choise for icea correction
      integer :: ncpmea    = 3 ! choise for pmea correction
      integer :: ncroffa   = 1 ! choise for roff correction
      integer :: nctaua    = 3 ! choise for taua correction
      integer :: nchfla    = 3 ! choise for hfla correction
      integer :: ncswra    = 1 ! choise for swra correction
      integer :: ncpsa     = 4 ! choise for psa correction
      integer :: ncima     = 1 ! choise for ima correction
      integer :: ncssto    = 4 ! choise for ssto correction
      integer :: ncssso    = 4 ! choise for ssso correction
      integer :: nctso     = 4 ! choise for tso correction
      integer :: nciceo    = 2 ! choise for iceo correction
      integer :: nchflo    = 3 ! choise for hflo correction
      integer :: ncswro    = 1 ! choise for swro correction
      integer :: ncpso     = 4 ! choise for pso correction
      integer :: ncimo     = 1 ! choise for imo correction
      integer :: newmask   = 0 ! switch to create a new alsm
!
      real :: cfacsst   = 1. ! coupling factor sst
      real :: cfactau   = 1. ! coupling factor stress
      real :: cfacice   = 1. ! coupling factor ice
      real :: cfacfresh = 1. ! coupling factor fresh water
!
      character(len=80) :: runoffmap='runoffmap.txt' ! file for runoff map
!
      real :: radeo  ! erath radius from the ocean model
!
      end module cplmod
!
!     ====================
!     subroutine CPLATMINI
!     ====================
!
      subroutine cplatmini(pslm,kxa,kya)
      use cplmod
      implicit none
!
!     initialize Plasim and the coupler (atm. side)
!
      real :: pslm(nxa,nya)
      real :: zgr
      integer :: kdatim(6)
      integer :: kcpl  = 2
      integer :: kxa,kya
      integer :: i,j,jt,ntape,kmiss
      logical :: lopen
      character(len=10) :: chform
      real :: zpi,zgw,zgm,zgs
!
!     dbug
!
      real :: zduma(nxa,nya) = 1
      real :: zdumo(nxo,nyo) = 1
      real :: zslmaos(nxo,nyo) 
      real :: zslmaov(nxo,nyo)
      real :: zslmoa(nxa,nya)
!
      character(len=8) :: cform
!
      namelist/cplpar/ndbox,ndboy,nprint,nissta,niicea,nipmea,nitaua    &
     &               ,nihfla,niswra,nissto,nihflo,niswro                &
     &               ,niroffa,ncroffa,nipsa,ncpso,nima,nimo,ncima,ncimo &
     &               ,nitsa,nitso,nctsa,nctso                           &
     &               ,ncssta,ncicea,ncpmea,nctaua,nchfla                &
     &               ,ncswra,nchflo,ncssto,ncswro,ncpso                 &
     &               ,cfacsst,cfacice,cfactau,cfacfresh,nroffc,newmask  &
     &               ,runoffmap
!
!     check atm dimensions
!
      if(kxa /= nxa .or. kya /= nya) then
       print*,'!ERROR! wrong Plasim dimensions in cplmod'
       stop
      endif
!
!     read namelist
!
      do jt=30,99
       ntape=jt
       inquire(unit=ntape,opened=lopen)
       if(.not. lopen) exit
      enddo
      open(ntape,file='cpl_namelist',form='formatted')
       read(ntape,cplpar)
       write (*,'(/," ****************************************")')
       write (*,'(" * CPLMOD ",a30," *")') trim(version)
       write (*,'(" ******************************************")')
       write (*,'(" * Namelist CPLPAR from <cpl_namelist>    *")')
       write (*,'(" ******************************************")')
       write(*,cplpar)
      close(ntape)
!
!     initialize coupling method from kcpl (= nlsg in oceanmod)
!     ncoupling=1 (atmospheric sst and oceanic hfl)
!     ncoupling=2 (atmospheric hfl and oceanic sst)
!     
      ncoupling=kcpl
!
!     initialize runoff correction
!
      if(nroffc > 1) then
       write(chform,'(A1,I6.6,A3)') '(',nxa,'I1)'
       open(ntape,file=trim(runoffmap),form='formatted')
       read(ntape,chform,end=9000) iroffi(:,:)
       if(nroffc == 3) then
        read(ntape,*,end=9000)
        read(ntape,chform,end=9000) iroffo(:,:)
       endif
       close(ntape)
      endif
!
!     make plasim grid
!
      call cinigrida
!
!     put atm land sea mask to common
!
      call cputlsma(pslm,nxa,nya)
!
!     put atm time to common
!
!      call updatim(0)
!      kdatim = ndatim
      kdatim(:) = 0
      call cputtime(kdatim)
!
!     initialize/compute index-arrays for bi-linear interpolation  
!
!     a) ocean scalar -> atm. 
!
      call cinter2_ini(oslms,wos2a,ios2a,ogws,xos,xa                    &
     &                ,yos,ya,nxo,nxa,nyo,nya,nprint)
!
!     b) atm. -> ocean scalar
!
      call cinter2_ini(aslm,wa2os,ia2os,agw,xa,xos                      &
     &                ,ya,yos,nxa,nxo,nya,nyo,nprint)
!
!     c) atm. -> ocean u
!
      call cinter2_ini(aslm,wa2ou,ia2ou,agw,xa,xou                      &
     &                ,ya,you,nxa,nxo,nya,nyo,nprint)
!
!     d) atm. -> ocean v
!
      call cinter2_ini(aslm,wa2ov,ia2ov,agw,xa,xov                      &
     &                ,ya,yov,nxa,nxo,nya,nyo,nprint)
!
!     initialize/compute index-arrays for area overlap interpolation
!
!     a) ocean scalar -> atm.
!
      call cinter_ini(oslms,wos2a_a,inos2a,ios2a_a,xos1,xos2,xa1,xa2    &
     &               ,yos1,yos2,ya1,ya2,nxo,nxa,nyo,nya,nimax,nprint)
!
!     b) atm. -> ocean scalar
!
      call cinter_ini(aslm,wa2os_a,ina2os,ia2os_a,xa1,xa2,xos1,xos2     &
     &               ,ya1,ya2,yos1,yos2,nxa,nxo,nya,nyo,nimax,nprint)
!
!     c) atm. -> ocean u
!
      call cinter_ini(aslm,wa2ou_a,ina2ou,ia2ou_a,xa1,xa2,xou1,xou2     &
     &               ,ya1,ya2,you1,you2,nxa,nxo,nya,nyo,nimax,nprint)
!
!     d) atm. -> ocean v
!
      call cinter_ini(aslm,wa2ov_a,ina2ov,ia2ov_a,xa1,xa2,xov1,xov2     &
     &               ,ya1,ya2,yov1,yov2,nxa,nxo,nya,nyo,nimax,nprint)
!
!     dbug print out
!
      if(nprint > 0) then
       zgr=180./(2.*asin(1.))
       write(*,*) ' '
       write(*,*) 'End of CMITINI: '
!
       write(*,*) 'Atmosphere y-coordiantes (x=1)'
       do j=1,nya
        write(*,*) 'j=',j,' Y= ',ya(1,j)*zgr
       enddo
!
       write(*,*) 'Ocean y-coordiantes (x=1)'
       do j=1,nyo
        write(*,*) 'j=',j,' ys= ',yos(1,j)*zgr,' yu= ',you(1,j)*zgr     &
      &                  ,' yv= ',yov(1,j)*zgr
       enddo
!
       write(cform,'(A1,I2.2,A3)') '(',nxa,'I1)'
       write(*,*) 'Atmosphere Mask'
       write(*,cform) nint(aslm)
       write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
       write(*,*) 'Ocean Scalar Mask'
       write(*,cform) nint(oslms)
       write(*,*) 'Ocean u-velocity MasK'
       write(*,cform) nint(oslmu)
       write(*,*) 'Ocean v-velocity MasK'
       write(*,cform) nint(oslmv)
       write(*,*) ' '
       zpi=2.*asin(1.)
       call cglm(aslm,zduma,zgs,zgm,zgw,agw,nxa,nya)
       write(*,*) 'Atm. ocean area (sum/mean) + global area : '         &
     &            ,zgs*radea*radea*4.*zpi,zgm,zgw*radea*radea*4.*zpi
       call cglm(oslms,zdumo,zgs,zgm,zgw,ogws,nxo,nyo)
       write(*,*) 'Oce. ocean area scalar (sum/mean) + global area : '  &
     &            ,zgs*radeo*radeo*4.*zpi,zgm,zgw*radeo*radeo*4.*zpi
       call cglm(oslmu,zdumo,zgs,zgm,zgw,ogwu,nxo,nyo)
       write(*,*) 'Oce. ocean area u (sum/mean) + global area : '       &
     &            ,zgs*radeo*radeo*4.*zpi,zgm,zgw*radeo*radeo*4.*zpi
       call cglm(oslmv,zdumo,zgs,zgm,zgw,ogwv,nxo,nyo)
       write(*,*) 'Oce. ocean area v (sum/mean) + global area : '       &
     &            ,zgs*radeo*radeo*4.*zpi,zgm,zgw*radeo*radeo*4.*zpi

!
!     check mask overlaps for area overlap interpolation:
!
       write(*,*) 'Check for area overlap interpolation'
       write(*,*)
!
!     a) atm. -> ocean scalar
!
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(SUM(wa2os_a(i,j,:)) < 0.5 .and. oslms(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(*,*) 'a) a2os:'
       write(*,*) kmiss                                                 &
     &           ,' Oce Scalar points must be extrapolated'
       if(kmiss > 0) then
        write(*,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(*,*) 'i,j= ',i,j,' Lon,Lat= ',xos(i,j)*zgr,yos(i,j)*zgr
         endif
        enddo
        enddo
        write(*,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(*,cform) nint(zdumo)
       endif
!
!     b) atm. -> ocean u
!
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(SUM(wa2ou_a(i,j,:)) < 0.5 .and. oslmu(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(*,*) 'b) a2ou:'
       write(*,*) kmiss                                                 &
     &           ,' Oce U-points must be extrapolated:'
       if(kmiss > 0) then
        write(*,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(*,*) 'i,j= ',i,j,' Lon,Lat= ',xou(i,j)*zgr,you(i,j)*zgr
         endif
        enddo
        enddo
        write(*,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(*,cform) nint(zdumo)
       endif
!
!     c) atm. -> ocean v
!
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(SUM(wa2ov_a(i,j,:)) < 0.5 .and. oslmv(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(*,*) 'c) a2ov:'
       write(*,*) kmiss                                                 &
     &           ,' Oce V-points must be extrapolated:'
       if(kmiss > 0) then
        write(*,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(*,*) 'i,j= ',i,j,' Lon,Lat= ',xov(i,j)*zgr,yov(i,j)*zgr
         endif
        enddo
        enddo
        write(*,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(*,cform) nint(zdumo)
       endif
!
!     d) ocean scalar -> atm
!
       kmiss=0
       do j=1,nya
       do i=1,nxa
        zduma(i,j)=0.
        if(SUM(wos2a_a(i,j,:)) < 0.5 .and. aslm(i,j) > 0.5) then
         kmiss=kmiss+1
         zduma(i,j)=1.
        endif
       enddo
       enddo
       write(*,*) 'd) os2a:'
       write(*,*) kmiss                                                 &
     &           ,' Atm points must be extrapolated:'
       if(kmiss > 0) then
        write(*,*) 'These points are:'
        do j=1,nya
        do i=1,nxa
         if(zduma(i,j) > 0.5) then
          write(*,*) 'i,j= ',i,j,' Lon,Lat= ',xa(i,j)*zgr,ya(i,j)*zgr
         endif
        enddo
        enddo
        write(*,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxa,'I1)'
        write(*,cform) nint(zduma)
       endif
       write(*,*)
!
!     check mask overlaps for bi-linear interpolation:
!
       write(*,*) 'Check for bi-linear interpolation'
       write(*,*)
!
!     a) atm. -> ocean scalar
!
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(wa2os(i,j,1)+wa2os(i,j,2)+wa2os(i,j,3)+wa2os(i,j,4) < 0.5    &
     &    .and. oslms(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(*,*) 'a) a2os:'
       write(*,*) kmiss                                                 &
     &           ,' Oce Scalar points must be extrapolated'  
       if(kmiss > 0) then
        write(*,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(*,*) 'i,j= ',i,j,' Lon,Lat= ',xos(i,j)*zgr,yos(i,j)*zgr
         endif
        enddo
        enddo
        write(*,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(*,cform) nint(zdumo)
       endif
!
!     b) atm. -> ocean u
!
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(wa2ou(i,j,1)+wa2ou(i,j,2)+wa2ou(i,j,3)+wa2ou(i,j,4) < 0.5    &
     &    .and. oslmu(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(*,*) 'b) a2ou:'
       write(*,*) kmiss                                                 &
     &           ,' Oce u-points must be extrapolated'
       if(kmiss > 0) then
        write(*,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(*,*) 'i,j= ',i,j,' Lon,Lat= ',xou(i,j)*zgr,you(i,j)*zgr
         endif
        enddo
        enddo
        write(*,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(*,cform) nint(zdumo)
       endif
!
!     c) atm. -> ocean v
!
       kmiss=0
       do j=1,nyo
       do i=1,nxo
        zdumo(i,j)=0.
        if(wa2ov(i,j,1)+wa2ov(i,j,2)+wa2ov(i,j,3)+wa2ov(i,j,4) < 0.5    &
     &    .and. oslmv(i,j) > 0.5) then
         kmiss=kmiss+1
         zdumo(i,j)=1.
        endif
       enddo
       enddo
       write(*,*) 'c) a2ov:'
       write(*,*) kmiss                                                 &
     &           ,' Oce v-points must be extrapolated'
       if(kmiss > 0) then
        write(*,*) 'These points are:'
        do j=1,nyo
        do i=1,nxo
         if(zdumo(i,j) > 0.5) then
          write(*,*) 'i,j= ',i,j,' Lon,Lat= ',xov(i,j)*zgr,yov(i,j)*zgr
         endif
        enddo
        enddo
        write(*,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxo,'I1)'
        write(*,cform) nint(zdumo)
       endif
!
!     d) ocean scalar -> atm.
!
       kmiss=0
       do j=1,nya
       do i=1,nxa
        zduma(i,j)=0.
        if(wos2a(i,j,1)+wos2a(i,j,2)+wos2a(i,j,3)+wos2a(i,j,4) < 0.5    &
     &    .and. aslm(i,j) > 0.5) then
         kmiss=kmiss+1
         zduma(i,j)=1.
        endif
       enddo
       enddo
       write(*,*) 'd) os2a:'
       write(*,*) kmiss                                                 &
     &           ,' Atm points must be extrapolated'
       if(kmiss > 0) then
        write(*,*) 'These points are:'
        do j=1,nya
        do i=1,nxa
         if(zduma(i,j) > 0.5) then
          write(*,*) 'i,j= ',i,j,' Lon,Lat= ',xa(i,j)*zgr,ya(i,j)*zgr
         endif
        enddo
        enddo
        write(*,*) 'Map of missed points'
        write(cform,'(A1,I2.2,A3)') '(',nxa,'I1)'
        write(*,cform) nint(zduma)
       endif
       write(*,*)
      endif
!
!     if selected: compute a new aslm to ensure a better matching
!
      if(newmask > 0) call mknewmask
!
!     handle points to be extrapolated
!
!     1. area overlap
!
!     a) atm. -> ocean scalar
!
! 
      if(nprint > 0) then
       write(*,*) 'going for the extrapolation: '
       write(*,*) '1) area overlap'
       write(*,*) '1.a) a2os' 
      endif
      call cextra(wa2os_a,ina2os,ia2os_a,oslms,nxo,nyo,nimax,nprint) 
!
!     b) atm. -> ocean u
!
      if(nprint > 0) write(*,*) '1.b) a2ou'
      call cextra(wa2ou_a,ina2ou,ia2ou_a,oslmu,nxo,nyo,nimax,nprint)
!
!     c) atm. -> ocean v
!
      if(nprint > 0) write(*,*) '1.c) a2ov'
      call cextra(wa2ov_a,ina2ov,ia2ov_a,oslmv,nxo,nyo,nimax,nprint)
!
!     d) ocean scalar -> atm
!
      if(nprint > 0) write(*,*) '1.d) os2a'
      call cextra(wos2a_a,inos2a,ios2a_a,aslm,nxa,nya,nimax,nprint)
!
!     2. bi linear
!
!     a) atm. -> ocean scalar
!
      if(nprint > 0) then
       write(*,*) '2) bi-linear'
       write(*,*) '2.a) a2os'
      endif
      call cextra2(wa2os,ia2os,oslms,nxo,nyo,4,nprint)
!
!     b) atm. -> ocean u
!
      if(nprint > 0) write(*,*) '2.b) a2ou'
      call cextra2(wa2ou,ia2ou,oslmu,nxo,nyo,4,nprint)
!
!     c) atm. -> ocean v
!
      if(nprint > 0) write(*,*) '2.c) a2ov'
      call cextra2(wa2ov,ia2ov,oslmv,nxo,nyo,4,nprint)
!
!     d) ocean scalar -> atm
!
      if(nprint > 0) write(*,*) '2.d) os2a'
      call cextra2(wos2a,ios2a,aslm,nxa,nya,4,nprint)
!
      return
!
 9000 continue
      write(*,*) 'ERROR: nroffc set to ',nroffc,' but file '            &
     &           ,trim(runoffmap)                                       &
     &           ,'not correctly given!'
      stop
!
      end subroutine cplatmini
!
!     ====================
!     subroutine CINIGRIDA
!     ====================
!
      subroutine cinigrida
      use cplmod
      implicit none
!
      real :: zsi(nya),zgw(nya)
      real :: zpi,zdxa,zsin,zsinm
      integer :: i,j
!
      zpi=2.*asin(1.)
!
!     atmospheric grid
!
      call cgauaw(zsi,zgw,nya)
      zdxa=2.*zpi/real(nxa)
      do i=1,nxa
      do j=1,nya
       xa(i,j)=zdxa*(i-1)
       ya(i,j)=asin(zsi(j))
      enddo
      enddo
      do i=2,nxa
       xa1(i,:)=0.5*(xa(i,:)+xa(i-1,:))
      enddo
      xa1(1,:)=xa1(2,:)-zdxa
      do i=1,nxa-1
       xa2(i,:)=0.5*(xa(i,:)+xa(i+1,:))
      enddo
      xa2(nxa,:)=xa2(nxa-1,:)+zdxa
      zgw(:)=2.*zgw(:)/sum(zgw(:))
      zsinm=sin(zpi*0.5)
      ya1(:,1)=asin(zsinm)
      do j=2,nya
       zsin=zsinm-zgw(j-1)
       ya1(:,j)=asin(zsin)
       ya2(:,j-1)=asin(zsin)
       zsinm=zsin
      enddo
      ya2(:,nya)=-zpi*0.5
!
!     (normalized) weights
!
      do j=1,nya
      do i=1,nxa
       agw(i,j)=zdxa*abs(sin(ya1(i,j))-sin(ya2(i,j)))
      enddo
      enddo
      agw(:,:)=agw(:,:)/SUM(agw(:,:))
!
      return
      end
!
!     =====================
!     subroutine CINIGRIDO1
!     =====================
!
      subroutine cinigrido1(pxc,pyc,pxg,pyg,pas,pau,pav,pms,pmu,pmv     &
     &                     ,pradeo,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
      integer :: kx,ky,ki,kj,myprocid
      real(kind=8) :: pxc(kx,ky),pyc(kx,ky)    ! x,y at centre of the c-grid (t-point)
      real(kind=8) :: pxg(kx,ky),pyg(kx,ky)    ! x,y at s-w corner of the c-grid
      real(kind=8) :: pas(kx,ky),pau(kx,ky),pav(kx,ky)  ! areas of s,u,v-point
      real(kind=8) :: pms(kx,ky),pmu(kx,ky),pmv(kx,ky)  ! l-s mask s,u,v-point
      real(kind=8) :: pradeo
      real(kind=8) :: zpf
      integer :: kx1,kxn,ky1,kyn
!
!     return if not root 
!
      if(myprocid /= NROOT) return
!
      zpf=2.*asin(1.d+0)/180.
!
      radeo=pradeo
!
      kx1=(ki-1)*kx+1
      kxn=ki*kx
      ky1=(kj-1)*ky+1
      kyn=kj*ky
!
      xos(kx1:kxn,ky1:kyn)=pxc(:,:)*zpf
      yos(kx1:kxn,ky1:kyn)=pyc(:,:)*zpf
      xou(kx1:kxn,ky1:kyn)=pxg(:,:)*zpf
      you(kx1:kxn,ky1:kyn)=pyc(:,:)*zpf
      xov(kx1:kxn,ky1:kyn)=pxc(:,:)*zpf
      yov(kx1:kxn,ky1:kyn)=pyg(:,:)*zpf
!
      ogws(kx1:kxn,ky1:kyn)=pas(:,:)
      ogwu(kx1:kxn,ky1:kyn)=pau(:,:)
      ogwv(kx1:kxn,ky1:kyn)=pav(:,:)
!
      oslms(kx1:kxn,ky1:kyn)=pms(:,:)
      oslmu(kx1:kxn,ky1:kyn)=pmu(:,:)
      oslmv(kx1:kxn,ky1:kyn)=pmv(:,:)
!
      return
      end
!
!     =====================
!     subroutine CINIGRIDO2
!     =====================
!
      subroutine cinigrido2(myprocid)
      use cplmod
      implicit none
      real(kind=8) :: zpi
      integer :: myprocid
!
!     return if not root
!
      if(myprocid /= NROOT) return
!
      xos1(:,:)=xou(:,:)
      xos2(1:nxo-1,:)=xou(2:nxo,:)
      xos2(nxo,:)=xou(1,:)
      xou1(2:nxo,:)=xos(1:nxo-1,:)
      xou1(1,:)=xos(nxo,:)
      xou2(:,:)=xos(:,:)
      xov1(:,:)=xou(:,:)
      xov2(1:nxo-1,:)=xou(2:nxo,:)
      xov2(nxo,:)=xou(1,:)
      yos1(:,:)=yov(:,:)
      yos2(:,1:nyo-1)=yov(:,2:nyo)
      yos2(:,nyo)=-yos1(:,1)
      you1(:,:)=yov(:,:)
      you2(:,1:nyo-1)=yov(:,2:nyo)
      you2(:,nyo)=-you1(:,1)
      yov1(:,2:nyo)=yos(:,1:nyo-1)
      yov1(:,1)=yov(:,1)                ! closed boundary at y=1
      yov2(:,:)=yos(:,:)
!
      zpi=2.*asin(1.d+0)
      where(xos1(:,:) > xos(:,:)) xos1(:,:)=xos1(:,:)-2.*zpi
      where(xou1(:,:) > xou(:,:)) xou1(:,:)=xou1(:,:)-2.*zpi
      where(xov1(:,:) > xov(:,:)) xov1(:,:)=xov1(:,:)-2.*zpi
      where(xos2(:,:) < xos(:,:)) xos2(:,:)=2.*zpi+xos2(:,:)
      where(xou2(:,:) < xou(:,:)) xou2(:,:)=2.*zpi+xou2(:,:)
      where(xov2(:,:) < xov(:,:)) xov2(:,:)=2.*zpi+xov2(:,:)
!
      ogws(:,:)=ogws(:,:)/SUM(ogws(:,:))
      ogwu(:,:)=ogwu(:,:)/SUM(ogwu(:,:))
      ogwv(:,:)=ogwv(:,:)/SUM(ogwv(:,:))
!
      return
      end

!     =======================  
!     subroutine reset_cplmod
!     =======================

      subroutine reset_cplmod
        use cplmod
        implicit none
!
        asst(:,:)   = 0. ! atm. sst
        ats(:,:)    = 0. ! atm. tsurf
        aicec(:,:)  = 0. ! atm. ice concentration
        aiced(:,:)  = 0. ! atm. ice thickness
        apme(:,:)   = 0. ! atm. fresh water flux
        aroff(:,:)  = 0. ! atm. runoff
        ataux(:,:)  = 0. ! atm. u-stress
        atauy(:,:)  = 0. ! atm. v-stress
        ahfl(:,:)   = 0. ! atm. heat flux
        aswr(:,:)   = 0. ! atm. solar radiation
        aps(:,:)    = 0. ! atm. surf. pressure
        aim(:,:)    = 0. ! atm. (ice+snow) mass 
!
        osst(:,:)   = 0. ! oce. sst
        osss(:,:)   = 0. ! oce. sss
        ots(:,:)    = 0. ! oce. tsurf
        oicec(:,:)  = 0. ! oce. ice concentration
        oiced(:,:)  = 0. ! oce. ice thickness
        opme(:,:)   = 0. ! oce. fresh water flux
        otaux(:,:)  = 0. ! oce. u-stress
        otauy(:,:)  = 0. ! oce. v-stress
        ohfl(:,:)   = 0. ! oce. heat flux
        oswr(:,:)   = 0. ! oce. solar radiation
        ops(:,:)    = 0. ! oce. surf. pressure
        oim(:,:)    = 0. ! oce. (ice+snow) mass
!
        assta(:,:)  = 0. ! atm. sst (accumulated)
        asss(:,:)   = 0. ! atm. sss (accumulated)
        atsa(:,:)   = 0. ! atm. tsurf (accumulated)
        aiceca(:,:) = 0. ! atm. ice concentration (acumulated)
        aiceda(:,:) = 0. ! atm. ice thickness (accumulated)
        apmea(:,:)  = 0. ! atm. frsh water flux (accumulated)
        aroffa(:,:) = 0. ! atm. runoff (accumulated)
        atauxa(:,:) = 0. ! atm. u-stress (accumulated)
        atauya(:,:) = 0. ! atm. v-stress (accumulated)
        ahfla(:,:)  = 0. ! atm. heat flux (accumulated)
        aswra(:,:)  = 0. ! atm. solar radiation (accumulated)
        apsa(:,:)   = 0. ! atm. surf. pressure (accumulated)
        aima(:,:)   = 0. ! atm. (ice+snow) mass (accumulated)
!
        ossta(:,:)  = 0. ! oce. sst (accumulated)
        osssa(:,:)  = 0. ! oce. sss (accumulated)
        otsa(:,:)   = 0. ! oce. tsurf (accumulated)
        oiceca(:,:) = 0. ! oce. ice concentration (accumulated)
        oiceda(:,:) = 0. ! oce. ice thickness(accumulated)
        opmea(:,:)  = 0. ! oce. fresh water flux (accumulated)
        otauxa(:,:) = 0. ! oce. u-stress (accumulated)
        otauya(:,:) = 0. ! oce. v-stress (accumulated)
        ohfla(:,:)  = 0. ! oce. heat flux (accumulated)
        oswra(:,:)  = 0. ! oce. solar radiation (accumulated)
        opsa(:,:)   = 0. ! oce. surface pressure (accumulated)
        oima(:,:)   = 0. ! oce. (ice+snow) mass (accumulated)
!
        nssta = 0 ! counter for atm. sst accumulation
        ntsa  = 0 ! counter for atm. tsurf accumulation
        nicea = 0 ! counter for atm. ice accumulation
        npmea = 0 ! counter for atm. fresh water flux accumulation
        nroffa= 0 ! counter for atm. runoff
        ntaua = 0 ! counter for atm. stress accumulation
        nhfla = 0 ! counter for atm. heat flux accumulation
        nswra = 0 ! counter for atm. solar radiation
        npsa  = 0 ! counter for atm. surf. pressure
        nima  = 0 ! counter for atm. (ice+snow) mass 
        nssto = 0 ! counter for oce. sst accumulation
        nssso = 0 ! counter for oce. sss accumulation
        ntso  = 0 ! counter for oce. tsurf accumulation
        niceo = 0 ! counter for oce. ice accumulation
        npmeo = 0 ! counter for oce. fresh water flux accumulation
        ntauo = 0 ! counter for oce. stress accumulation 
        nhflo = 0 ! counter for oce. heat flux accumulation
        nswro = 0 ! counter for oce. solar radiation
        npso  = 0 ! counter for oce. surface pressure
        nimo  = 0 ! counter for oce. (ice+snow) mass 
!
      end subroutine reset_cplmod

!
!     ===================
!     subroutine CPUTTIME
!     ===================
!
      subroutine cputtime(kdatim)
      use cplmod
      implicit none
!
!     put time to common
!
      integer :: kdatim(6)
!
      ndatim(:)=kdatim(:)
      nyear=ndatim(1)
      nmonth=ndatim(2)
      nday=ndatim(3)
      nhour=ndatim(4)
      nmin=ndatim(5)
      nwday=ndatim(6)
!      
      return
      end
!
!     ===================
!     subroutine CGETTIME
!     ===================
!
      subroutine cgettime(kdatim)
      use cplmod
      implicit none
!
!     get time from common
!
      integer :: kdatim(6)
!
      kdatim(:)=ndatim(:)
!
      if(nprint > 0) then
       write(*,*) ' '
       write(*,*) 'In CGETTIME, datim= ',ndatim(:)
       write(*,*) ' '
      endif
!       
      return
      end
!
!     ===================
!     subroutine CPUTLSMA
!     ===================
!
      subroutine cputlsma(plsm,kx,ky)
      use cplmod
      implicit none
!
!     put atm land see mask to common
!
      integer, intent(in) :: kx,ky
      real, intent(in) :: plsm(nxa,nya)      ! land see mask (input)
!
!     check dimensions
!
      if(nxa /= kx .or. nya /= ky) then
       write(*,*) '!ERROR! inconsistent atm. dimensions in cpl'
       stop
      endif
!
!     copy land see mask and convert 0/1 to 1/0
!
      aslm(:,:)=real(nint(1.-plsm(:,:)))
!
      return
      end
!
!     ===================
!     subroutine CPUTLSMO
!     ===================
!
      subroutine cputlsmo(plsms,plsmv,kx,ky)
      use cplmod
      implicit none
!
!     put oce land see mask to common
!
      integer, intent(in) :: kx,ky
      real :: plsms(nxo,nyo)   ! land see mask scalar (input)
      real :: plsmv(nxo,nyo)   ! land see mask vector (input)
!
!     check dimensions
!
      if(nxo /= kx .or. nyo /= ky) then
       write(*,*) '!ERROR! inconsistent oce. dimensions in cpl'
       stop
      endif
!
!     copy land see mask and convert 0/1 to 1/0
!
      oslms(1:nxo,1:nyo)=real(nint(plsms(1:nxo,1:nyo)))
      oslmv(1:nxo,1:nyo)=real(nint(plsmv(1:nxo,1:nyo)))
!
      return
      end
!
!     ===================
!     subroutine CPUTPMEA
!     ===================
!
      subroutine cputpmea(ppme)
      use cplmod
      implicit none
!
!     put atm fresh water flux to common (accumulated)
!
      real :: ppme(nxa,nya) ! fresh water flux (input)
!
!     accumulate
!
      apmea(:,:)=apmea(:,:)+ppme(:,:)
      npmea=npmea+1
!
      return
      end
!
!     ===================
!     subroutine CPUTSWRA
!     ===================
!
      subroutine cputswra(pswr)
      use cplmod
      implicit none
!
!     put atm short wave radiation to common (accumulated)
!
      real :: pswr(nxa,nya) ! short wave radiaton (input)
!
!     accumulate
!
      aswra(:,:)=aswra(:,:)+pswr(:,:)
      nswra=nswra+1
!
      return
      end
!
!     ==================
!     subroutine CPUTPSA
!     ==================
!
      subroutine cputpsa(pps)
      use cplmod
      implicit none
!
!     put atm surf. pressure to common (accumulated)
!
      real :: pps(nxa,nya) ! surf. pressure (input)
!
!     accumulate
!
      apsa(:,:)=apsa(:,:)+pps(:,:)
      npsa=npsa+1
!
      return
      end
!
!     ==================
!     subroutine CPUTIMA
!     ==================
!
      subroutine cputima(pim)
      use cplmod
      implicit none
!
!     put atm surf. pressure to common (accumulated)
!
      real :: pim(nxa,nya) ! ice + snow mass (input)
!
!     accumulate
!
      aima(:,:)=aima(:,:)+pim(:,:)
      nima=nima+1
!
      return
      end
!
!     ===================
!     subroutine CPUTTAUA
!     ===================
!
      subroutine cputtaua(ptaux,ptauy)
      use cplmod
      implicit none
!
!     put atm wind stress to common (accumulated)
!
      real :: ptaux(nxa,nya) ! u-stress (input)
      real :: ptauy(nxa,nya) ! v-stress (input)
!
!     accumulate
!
      atauxa(:,:)=atauxa(:,:)+ptaux(:,:)
      atauya(:,:)=atauya(:,:)+ptauy(:,:)
      ntaua=ntaua+1
!
      return
      end
!
!     ===================
!     subroutine CPUTHFLA
!     ===================
!
      subroutine cputhfla(phfl)
      use cplmod
      implicit none
!
!     put atm heat flux to common (accumulated)
!
      real :: phfl(nxa,nya) ! heat flux (input)
!
!     accumulate
!
      ahfla(:,:)=ahfla(:,:)+phfl(:,:)
      nhfla=nhfla+1
!
      return
      end
!
!     ====================
!     subroutine CPUTROFFA
!     ====================
!
      subroutine cputroffa(proff)
      use cplmod
      implicit none
!
!     put atm runoff to common (accumulated)
!
      real :: proff(nxa,nya) ! runoff (input)
      real :: zroff(nxa,nya) ! runoff (local)
!
      zroff(:,:)=proff(:,:)
!
!     remap runoff (if needed)
!
      if(nroffc > 0) call croffc(zroff)
!
!     accumulate
!
      aroffa(:,:)=aroffa(:,:)+zroff(:,:)
      nroffa=nroffa+1
!
      return
      end

!
!     ===================
!     subroutine CPUTHFLO
!     ===================
!
      subroutine cputhflo(phfl)
      use cplmod
      implicit none
!
!     put oce heat flux to common (accumulated)
!
      real(kind=8) :: phfl(nxo,nyo) ! heat flux (input)
!
!     accumulate
!
      ohfla(:,:)=ohfla(:,:)+phfl(:,:)
      nhflo=nhflo+1
!
      return
      end
!
!     ===================
!     subroutine CPUTSSTO
!     ===================
!
      subroutine cputssto(psst,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
!     put oce sst to common (accumulated)
!
      real(kind=8) :: psst(kx,ky) ! ocean sst (input)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
      integer :: myprocid         ! process id
!
      integer :: kx1,kxn,ky1,kyn
!
!     return if not root
!
      if(myprocid /= NROOT) return
!
!     compute indexes within global ocean field
!
      kx1=(ki-1)*kx+1
      kxn=ki*kx
      ky1=(kj-1)*ky+1
      kyn=kj*ky
!
      ossta(kx1:kxn,ky1:kyn)=ossta(kx1:kxn,ky1:kyn)+psst(:,:)
      nssto=nssto+1
!
      return
      end
!
!     ===================
!     subroutine CPUTSSSO
!     ===================
!
      subroutine cputssso(psss,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
!     put oce sst to common (accumulated)
!
      real(kind=8) :: psss(kx,ky) ! ocean sss (input)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
      integer :: myprocid         ! process id
!
      integer :: kx1,kxn,ky1,kyn
!
!     return if not root
!
      if(myprocid /= NROOT) return
!
!     compute indexes within global ocean field
!
      kx1=(ki-1)*kx+1
      kxn=ki*kx
      ky1=(kj-1)*ky+1
      kyn=kj*ky
!
      osssa(kx1:kxn,ky1:kyn)=osssa(kx1:kxn,ky1:kyn)+psss(:,:)
      nssso=nssso+1
!
      return
      end
!
!     ===================
!     subroutine CPUTTSO
!     ===================
!
      subroutine cputtso(pts,kx,ky,ki,kj)
      use cplmod
      implicit none
!
!     put oce tsurf to common (accumulated)
!
      real(kind=8) :: pts(kx,ky)  ! ocean tsurf (input)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
!
      integer :: kx1,kxn,ky1,kyn
!
!     compute indexes within global ocean field
!
      kx1=(ki-1)*kx+1
      kxn=ki*kx
      ky1=(kj-1)*ky+1
      kyn=kj*ky
!
      otsa(kx1:kxn,ky1:kyn)=otsa(kx1:kxn,ky1:kyn)+pts(:,:)
      ntso=ntso+1
!
      return
      end
!
!     ===================
!     subroutine CPUTICEO
!     ===================
!
      subroutine cputiceo(picec,piced,kx,ky,ki,kj)
      use cplmod
      implicit none
!
!     put oce ice to common (accumulated)
!
      real(kind=8) :: picec(kx,ky)! ocean ice frac. (input)
      real(kind=8) :: piced(kx,ky)! ocean ice thick. (input)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
!
      integer :: kx1,kxn,ky1,kyn
!
!     compute indexes within global ocean field
!
      kx1=(ki-1)*kx+1
      kxn=ki*kx
      ky1=(kj-1)*ky+1
      kyn=kj*ky
!
      oiceca(kx1:kxn,ky1:kyn)=oiceca(kx1:kxn,ky1:kyn)+picec(:,:)
      oiceda(kx1:kxn,ky1:kyn)=oiceda(kx1:kxn,ky1:kyn)+piced(:,:)
      niceo=niceo+1
!
      return
      end

!
!     ===================
!     subroutine CGETSSTO
!     ===================
!
      subroutine cgetssto(psst)
      use cplmod
      implicit none
      real, parameter :: rkelvin=273.15
!
!     get interpolated sst from common 
!
      real :: psst(nxo,nyo) ! sst (output)
!
      real :: zgm,zgs,zgw
!
!     get accumulated atm sst
!
      if(nssta == 0) then
       write(*,*) '!ERROR! no accumulated atm sst found'
       stop
      endif
      asst(:,:)=assta(:,:)/real(nssta)
!
!     interpolate
!
      if(nissta== 1) then
       call cinter(asst,aslm,wa2os_a,ina2os,ia2os_a,osst,oslms          &
     &            ,agw,ogws,nxa,nxo,nya,nyo,ncssta,nimax,nprint)
      else
       call cinter2(asst,aslm,wa2os,ia2os,osst,oslms,agw,ogws           &
     &             ,nxa,nxo,nya,nyo,ncssta,nprint)
      endif
!
!     copy to output (convert into C)
!
      psst(:,:)=osst(:,:)-rkelvin
!
!     dbug printout
!
      if(nprint > 0) then
       write(*,*) ' '
       write(*,*) 'In CGETSSTO: '
       write(*,*) 'sst from ',nssta,' accumulated values'
       call cglm(asst,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(*,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(asst,MASK=(aslm > 0.5))               &
     &            ,MINVAL(asst,MASK=(aslm > 0.5))
       call cglm(osst,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
       write(*,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(osst,MASK=(oslms > 0.5))              &
     &            ,MINVAL(osst,MASK=(oslms > 0.5))
       write(*,*) ' '
      endif
!
!     reset accumulated field
!
      assta(:,:)=0.
      nssta=0
!
      return
      end
!
!     ===================
!     subroutine CGETPMEO
!     ===================
!
      subroutine cgetpmeo(ppme,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
!     get interpolated atm pme from common 
!
      real(kind=8) :: ppme(kx,ky) ! P-E (output)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
      integer :: myprocid         ! process id
!
      integer :: kx1,kxn,ky1,kyn
!
      real :: zpme(nxo,nyo) ! p-e (global field)
      real :: zgm,zgs,zgw
!
!     do the interpolation for root only
!
      if(myprocid == NROOT) then
!
!     compute indexes within global ocean field
!
       kx1=(ki-1)*kx+1
       kxn=ki*kx
       ky1=(kj-1)*ky+1
       kyn=kj*ky
!
!     get accumulated atm pme
!
       if(npmea == 0) then
        write(*,*) '!ERROR! no accumulated atm pme found'
        stop
       endif
       apme(:,:)=apmea(:,:)/real(npmea)
!
!     interpolate
!
       if(nipmea == 1) then
        call cinter(apme,aslm,wa2os_a,ina2os,ia2os_a,zpme,oslms         &
     &             ,agw,ogws,nxa,nxo,nya,nyo,ncpmea,nimax,nprint)
       else
        call cinter2(apme,aslm,wa2os,ia2os,zpme,oslms,agw,ogws          &
     &              ,nxa,nxo,nya,nyo,ncpmea,nprint)
       endif
!
!     copy to output
!
       ppme(:,:)=-1.*zpme(kx1:kxn,ky1:kyn)
!
!     dbug printout
!
       if(nprint > 0) then
        write(*,*) ' '
        write(*,*) 'In CGETPMEO: '
        write(*,*) 'pme from ',npmea,' accumulated values'
        call cglm(apme,aslm,zgs,zgm,zgw,agw,nxa,nya)
        write(*,*) 'Atm (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(apme,MASK=(aslm > 0.5))              &
     &             ,MINVAL(apme,MASK=(aslm > 0.5))
        call cglm(zpme,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
        write(*,*) 'Oce (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(zpme,MASK=(oslms > 0.5))             &
     &            ,MINVAL(zpme,MASK=(oslms > 0.5))
        write(*,*) ' '
       endif
!
!     reset accumulated field
!
       apmea(:,:)=0.
       npmea=0
      endif
!
!     broadcast to the rest of the world
!
      call mpbcrn(ppme,kx*ky)
!
      return
      end
!
!     ===================
!     subroutine CGETSWRO
!     ===================
!
      subroutine cgetswro(pswr,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
!     get interpolated atm swr from common 
!
      real(kind=8) :: pswr(kx,ky) ! swr (output)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
      integer :: myprocid         ! process id
!
      real :: zswr(nxo,nyo) ! swr (global field)
      real :: zgm,zgs,zgw
      integer :: kx1,kxn,ky1,kyn

!
!     do the interpolation for root only
!
      if(myprocid == NROOT) then
!
!     compute indexes within global ocean field
!
       kx1=(ki-1)*kx+1
       kxn=ki*kx
       ky1=(kj-1)*ky+1
       kyn=kj*ky
!
!     get accumulated atm pme
!
       if(nswra == 0) then
        write(*,*) '!ERROR! no accumulated atm swr found'
        stop
       endif
       aswr(:,:)=aswra(:,:)/real(nswra)
!
!     interpolate
!
       if(niswra == 1) then
        call cinter(aswr,aslm,wa2os_a,ina2os,ia2os_a,zswr,oslms         &
     &             ,agw,ogws,nxa,nxo,nya,nyo,ncswra,nimax,nprint)
       else
        call cinter2(aswr,aslm,wa2os,ia2os,zswr,oslms,agw,ogws          &
     &              ,nxa,nxo,nya,nyo,ncswra,nprint)
       endif
!
!     copy to output
!
       pswr(:,:)=-1.*zswr(kx1:kxn,ky1:kyn)
!
!     dbug printout
!
       if(nprint > 0) then
        write(*,*) ' '
        write(*,*) 'In CGETSWRO: '
        write(*,*) 'swr from ',nswra,' accumulated values'
        call cglm(aswr,aslm,zgs,zgm,zgw,agw,nxa,nya)
        write(*,*) 'Atm (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(aswr,MASK=(aslm > 0.5))              &
     &             ,MINVAL(aswr,MASK=(aslm > 0.5))
        call cglm(zswr,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
        write(*,*) 'Oce (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(zswr,MASK=(oslms > 0.5))             &
     &            ,MINVAL(zswr,MASK=(oslms > 0.5))
        write(*,*) ' '
       endif
!
!     reset accumulated field
!
       aswra(:,:)=0.
       nswra=0
      endif
!
!     broadcast to the rest of the world
!
      call mpbcrn(pswr,kx*ky)
!
      return
      end
!
!     ==================
!     subroutine CGETPSO
!     ==================
!
      subroutine cgetpso(pps,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
!     get interpolated atm surf. pressure from common 
!
      real(kind=8) :: pps(kx,ky)  ! surf. pressure (output)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
      integer :: myprocid         ! process id
!
      real :: zps(nxo,nyo) ! ps (global field)
      real :: zgm,zgs,zgw
      integer :: kx1,kxn,ky1,kyn
!
!     do the interpolation for root only
!
      if(myprocid == NROOT) then
!
!     compute indexes within global ocean field
!
       kx1=(ki-1)*kx+1
       kxn=ki*kx
       ky1=(kj-1)*ky+1
       kyn=kj*ky
!
!     get accumulated atm pme
!
       if(npsa == 0) then
        write(*,*) '!ERROR! no accumulated atm ps found'
        stop
       endif
       aps(:,:)=apsa(:,:)/real(npsa)
!
!     interpolate
!
       if(nipsa == 1) then
        call cinter(aps,aslm,wa2os_a,ina2os,ia2os_a,zps,oslms           &
     &             ,agw,ogws,nxa,nxo,nya,nyo,ncpsa,nimax,nprint)
       else
        call cinter2(aps,aslm,wa2os,ia2os,zps,oslms,agw,ogws            &
     &              ,nxa,nxo,nya,nyo,ncpsa,nprint)
       endif
!
!     copy to output
!
       pps(:,:)=zps(kx1:kxn,ky1:kyn)
!
!     dbug printout
!
       if(nprint > 0) then
        write(*,*) ' '
        write(*,*) 'In CGETPSO: '
        write(*,*) 'surf. pres. from ',npsa,' accumulated values'
        call cglm(aps,aslm,zgs,zgm,zgw,agw,nxa,nya)
        write(*,*) 'Atm (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(aps,MASK=(aslm > 0.5))               &
     &             ,MINVAL(aps,MASK=(aslm > 0.5))
        call cglm(zps,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
        write(*,*) 'Oce (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(zps,MASK=(oslms > 0.5))              &
     &             ,MINVAL(zps,MASK=(oslms > 0.5))
        write(*,*) ' '
       endif
!
!     reset accumulated field
!
       apsa(:,:)=0.
       npsa=0
      endif
!
!     broadcast to the rest of the world
!
      call mpbcrn(pps,kx*ky)
!
      return
      end
!
!     ==================
!     subroutine CGETIMO
!     ==================
!
      subroutine cgetimo(pim,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
!     get interpolated snow and ice mass from common 
!
      real(kind=8) :: pim(kx,ky)  ! mass (output)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
      integer :: myprocid         ! process id
!
      real :: zim(nxo,nyo) ! mass (global field)
      real :: zgm,zgs,zgw
      integer :: kx1,kxn,ky1,kyn
!
!     do the interpolation for root only
!
      if(myprocid == NROOT) then
!
!     compute indexes within global ocean field
!
       kx1=(ki-1)*kx+1
       kxn=ki*kx
       ky1=(kj-1)*ky+1
       kyn=kj*ky
!
!     get accumulated atm mass
!
       if(nima == 0) then
        write(*,*) '!ERROR! no accumulated mass found'
        stop
       endif
       aim(:,:)=aima(:,:)/real(nima)
!
!     interpolate
!
       if(niima == 1) then
        call cinter(aim,aslm,wa2os_a,ina2os,ia2os_a,zim,oslms           &
     &             ,agw,ogws,nxa,nxo,nya,nyo,ncima,nimax,nprint)
       else
        call cinter2(aim,aslm,wa2os,ia2os,zim,oslms,agw,ogws            &
     &             ,nxa,nxo,nya,nyo,ncima,nprint)
       endif
!
!     copy to output
!
       pim(:,:)=zim(kx1:kxn,ky1:kyn)
!
!     dbug printout
!
       if(nprint > 0) then
        write(*,*) ' '
        write(*,*) 'In CGETPSO: '
        write(*,*) 'mass from ',nima,' accumulated values'
        call cglm(aim,aslm,zgs,zgm,zgw,agw,nxa,nya)
        write(*,*) 'Atm (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(aim,MASK=(aslm > 0.5))               &
     &             ,MINVAL(aim,MASK=(aslm > 0.5))
        call cglm(zim,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
        write(*,*) 'Oce (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(zim,MASK=(oslms > 0.5))              &
     &             ,MINVAL(zim,MASK=(oslms > 0.5))
        write(*,*) ' '
       endif
!
!     reset accumulated field
!
       aima(:,:)=0.
       nima=0
      endif
!
!     broadcast to the rest of the world
!
      call mpbcrn(pim,kx*ky)
!
      return
      end
!
!     ====================
!     subroutine CGETROFFO
!     ====================
!
      subroutine cgetroffo(proff,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
!     get interpolated atm runoff from common 
!
      real(kind=8) :: proff(kx,ky)! runoff (output)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
      integer :: myprocid         ! process id
!
      real :: zroff(nxo,nyo) ! runoff (global field)
      real :: zgs,zgm,zgw
      integer :: kx1,kxn,ky1,kyn
!
!     do the interpolation for root only
!
      if(myprocid == NROOT) then
!
!     compute indexes within global ocean field
!
       kx1=(ki-1)*kx+1
       kxn=ki*kx
       ky1=(kj-1)*ky+1
       kyn=kj*ky
!
!     get accumulated atm pme
!
       if(nroffa == 0) then
        write(*,*) '!ERROR! no accumulated atm roff found'
        stop
       endif
       aroff(:,:)=aroffa(:,:)/real(nroffa)
!
!     interpolate
!
       if(niroffa == 1) then
        call cinter(aroff,aslm,wa2os_a,ina2os,ia2os_a,zroff,oslms       &
     &             ,agw,ogws,nxa,nxo,nya,nyo,ncroffa,nimax,nprint)
       else
        call cinter2(aroff,aslm,wa2os,ia2os,zroff,oslms,agw,ogws        &
     &              ,nxa,nxo,nya,nyo,ncroffa,nprint)
       endif
!
!     copy to output
!
       proff(:,:)=-1.*zroff(kx1:kxn,ky1:kyn)
!
!     dbug printout
!
       if(nprint > 0) then
        write(*,*) ' '
        write(*,*) 'In CGETROFFO: '
        write(*,*) 'runoff from ',nroffa,' accumulated values'
        call cglm(aroff,aslm,zgs,zgm,zgw,agw,nxa,nya)
        write(*,*) 'Atm (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(aroff,MASK=(aslm > 0.5))             &
     &            ,MINVAL(aroff,MASK=(aslm > 0.5))
        call cglm(zroff,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
        write(*,*) 'Oce (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(zroff,MASK=(oslms > 0.5))            &
     &            ,MINVAL(zroff,MASK=(oslms > 0.5))
        write(*,*) ' '
       endif
!
!     reset accumulated field
!
       aroffa(:,:)=0.
       nroffa=0
      endif 
!
!     broadcast to the rest of the world
!
      call mpbcrn(proff,kx*ky)
!
      return
      end
!
!     ===================
!     subroutine CGETTAUO
!     ===================
!
      subroutine cgettauo(ptaux,ptauy,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
!     get interpolated atm stress from common 
!
      real(kind=8) :: ptaux(kx,ky)! u-stress (output)
      real(kind=8) :: ptauy(kx,ky)! v-stress (output)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
      integer :: myprocid         ! process id
!
      real :: ztaux(nxo,nyo) ! u-stress (global field)
      real :: ztauy(nxo,nyo) ! v-stress (global field)
      real :: zgs,zgm,zgw
      integer :: kx1,kxn,ky1,kyn
!
!     do the interpolation for root only
!
      if(myprocid == NROOT) then
!
!     compute indexes within global ocean field
!
       kx1=(ki-1)*kx+1
       kxn=ki*kx
       ky1=(kj-1)*ky+1
       kyn=kj*ky
!
!     get accumulated atm stress
!
       if(ntaua == 0) then
        write(*,*) '!ERROR! no accumulated atm stress found'
        stop
       endif
       ataux(:,:)=atauxa(:,:)/real(ntaua)
       atauy(:,:)=atauya(:,:)/real(ntaua)
!
!     interpolate
!
       if(nitaua == 1) then
        call cinter(ataux,aslm,wa2ou_a,ina2ou,ia2ou_a,ztaux,oslmu       &
     &             ,agw,ogwu,nxa,nxo,nya,nyo,nctaua,nimax,nprint)
        call cinter(atauy,aslm,wa2ov_a,ina2ov,ia2ov_a,ztauy,oslmv       &
     &             ,agw,ogwv,nxa,nxo,nya,nyo,nctaua,nimax,nprint)
       else
        call cinter2(ataux,aslm,wa2ou,ia2ou,ztaux,oslmu,agw,ogwu        &
     &              ,nxa,nxo,nya,nyo,nctaua,nprint)
        call cinter2(atauy,aslm,wa2ov,ia2ov,ztauy,oslmv,agw,ogwv        &
     &              ,nxa,nxo,nya,nyo,nctaua,nprint)
       endif
!
!     copy to output
!
       ptaux(:,:)=ztaux(kx1:kxn,ky1:kyn)
       ptauy(:,:)=ztauy(kx1:kxn,ky1:kyn)
!
!     dbug printout
!
       if(nprint > 0) then
        write(*,*) ' '
        write(*,*) 'In CGETTAUO: '
        write(*,*) 'taux from ',ntaua,' accumulated values'
        call cglm(ataux,aslm,zgs,zgm,zgw,agw,nxa,nya)
        write(*,*) 'Atm (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(ataux,MASK=(aslm > 0.5))             &
     &             ,MINVAL(ataux,MASK=(aslm > 0.5))
        call cglm(ztaux,oslmu,zgs,zgm,zgw,ogwu,nxo,nyo)
        write(*,*) 'Oce (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(ztaux,MASK=(oslmu > 0.5))            &
     &            ,MINVAL(ztaux,MASK=(oslmu > 0.5))
        write(*,*) ' '
        write(*,*) 'tauy from ',ntaua,' accumulated values'
        call cglm(atauy,aslm,zgs,zgm,zgw,agw,nxa,nya)
        write(*,*) 'Atm (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(atauy,MASK=(aslm > 0.5))             &
     &             ,MINVAL(atauy,MASK=(aslm > 0.5))
        call cglm(ztauy,oslmv,zgs,zgm,zgw,ogwv,nxo,nyo)
        write(*,*) 'Oce (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(ztauy,MASK=(oslmv > 0.5))            &
     &             ,MINVAL(ztauy,MASK=(oslmv > 0.5))
        write(*,*) ' '
       endif
!
!     reset accumulated field
!
       atauxa(:,:)=0.
       atauya(:,:)=0.
       ntaua=0
      endif
!
!     broadcast to the rest of the world
!
      call mpbcrn(ptaux,kx*ky)
      call mpbcrn(ptauy,kx*ky)
!
      return
      end
!
!     ===================
!     subroutine CGETHFLO
!     ===================
!
      subroutine cgethflo(phfl,kx,ky,ki,kj,myprocid)
      use cplmod
      implicit none
!
!     get interpolated atm hfl from common 
!
      real(kind=8) :: phfl(kx,ky) ! hfl (output)
      integer :: kx               ! x-dim of input
      integer :: ky               ! y-dim of input
      integer :: ki               ! ocean tile i
      integer :: kj               ! ocean tile j
      integer :: myprocid         ! process id
!
      real :: zhfl(nxo,nyo) ! hfl (global field)
      real :: zgs,zgm,zgw
      integer :: kx1,kxn,ky1,kyn
!
!     do the interpolation for root only
!
      if(myprocid == NROOT) then
!
!     compute indexes within global ocean field
!
       kx1=(ki-1)*kx+1
       kxn=ki*kx
       ky1=(kj-1)*ky+1
       kyn=kj*ky
!
!     get accumulated atm hfl
!
       if(nhfla == 0) then
        write(*,*) '!ERROR! no accumulated atm hfl found'
        stop
       endif
       ahfl(:,:)=ahfla(:,:)/real(nhfla)
!
!     interpolate
!
       if(nihfla == 1) then
        call cinter(ahfl,aslm,wa2os_a,ina2os,ia2os_a,zhfl,oslms         &
     &             ,agw,ogws,nxa,nxo,nya,nyo,nchfla,nimax,nprint)
       else
        call cinter2(ahfl,aslm,wa2os,ia2os,zhfl,oslms,agw,ogws          &
     &              ,nxa,nxo,nya,nyo,nchfla,nprint)
       endif
!
!     copy to output
!     note: in MIT positive is upward
!
       phfl(:,:)=-1.*zhfl(kx1:kxn,ky1:kyn)
!
!     dbug printout
!
       if(nprint > 0) then
        write(*,*) ' '
        write(*,*) 'In CGETHFLO: '
        write(*,*) 'hfl from ',nhfla,' accumulated values'
        call cglm(ahfl,aslm,zgs,zgm,zgw,agw,nxa,nya)
        write(*,*) 'Atm (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(ahfl,MASK=(aslm > 0.5))              &
     &             ,MINVAL(ahfl,MASK=(aslm > 0.5))
        call cglm(zhfl,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
        write(*,*) 'Oce (sum,mean,max,min): '                           &
     &             ,zgs,zgm,MAXVAL(zhfl,MASK=(oslms > 0.5))             &
     &             ,MINVAL(zhfl,MASK=(oslms > 0.5))
        write(*,*) ' '
       endif
!
!     reset accumulated field
!
       ahfla(:,:)=0.
       nhfla=0
      endif
!
!     broadcast to the rest of the world
!
      call mpbcrn(phfl,kx*ky)
!
      return
      end
!
!     ===================
!     subroutine CGETHFLA
!     ===================
!
      subroutine cgethfla(phfl)
      use cplmod
      implicit none
!
!     get interpolated oce hfl from common 
!
      real :: phfl(nxa,nya) ! hfl (output)
!
      real :: zgm,zgs,zgw
!
!     get accumulated oce hfl
!
      if(nhflo == 0) then
       write(*,*) '!ERROR! no accumulated oce hfl found'
       stop
      endif
      ohfl(:,:)=ohfla(:,:)/real(nhflo)
!
!     interpolate
!
      if(nihflo == 1) then
       call cinter(ohfl,oslms,wos2a_a,inos2a,ios2a_a,ahfl,aslm          &
     &            ,ogws,agw,nxo,nxa,nyo,nya,nchflo,nimax,nprint)
      else
       call cinter2(ohfl,oslms,wos2a,ios2a,ahfl,aslm,ogws,agw           &
     &             ,nxo,nxa,nyo,nya,nchflo,nprint)
      endif
!
!     copy to output
!
      phfl(:,:)=ahfl(:,:)
!
!     dbug printout
!
      if(nprint > 0) then
       write(*,*) ' '
       write(*,*) 'In CGETHFLA: '
       write(*,*) 'hfl from ',nhflo,' accumulated values'
       call cglm(ohfl,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
       write(*,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(ohfl,MASK=(oslms > 0.5))              &
     &            ,MINVAL(ohfl,MASK=(oslms > 0.5))
       call cglm(ahfl,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(*,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(ahfl,MASK=(aslm > 0.5))               &
     &            ,MINVAL(ahfl,MASK=(aslm > 0.5))
       write(*,*) ' '
      endif
!
!     reset accumulated field
!
      ohfla(:,:)=0.
      nhflo=0
!
      return
      end
!
!     ===================
!     subroutine CGETSSTA
!     ===================
!
      subroutine cgetssta(psst)
      use cplmod
      implicit none
      real :: zgw,zgm,zgs
!
!     get interpolated oce sst from common 
!
      real :: psst(nxa,nya) ! sst (output)
!
!     get accumulated oce sst
!
      if(nssto == 0) then
       write(*,*) '!ERROR! no accumulated oce sst found'
       stop
      endif
      osst(:,:)=ossta(:,:)/real(nssto)
!
!     interpolate
!
      if (nissto == 1) then
       call cinter(osst,oslms,wos2a_a,inos2a,ios2a_a,asst,aslm          &
     &            ,ogws,agw,nxo,nxa,nyo,nya,ncssto,nimax,nprint)
      else
       call cinter2(osst,oslms,wos2a,ios2a,asst,aslm,ogws,agw           &
     &             ,nxo,nxa,nyo,nya,ncssto,nprint)

      endif
!
!     copy to output
!
      psst(:,:)=asst(:,:)
!
!     dbug printout
!
      if(nprint > 0) then
       write(*,*) ' '
       write(*,*) 'In CGETSSTA: '
       write(*,*) 'sst from ',nssto,' accumulated values'
       call cglm(osst,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
       write(*,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(osst,MASK=(oslms > 0.5))              &
     &            ,MINVAL(osst,MASK=(oslms > 0.5))
       call cglm(asst,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(*,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(asst,MASK=(aslm > 0.5))               &
     &            ,MINVAL(asst,MASK=(aslm > 0.5))
       write(*,*) ' '
      endif
!
!     reset accumulated field
!
      ossta(:,:)=0.
      nssto=0
!
      return
      end
!
!     ===================
!     subroutine CGETSSSA
!     ===================
!
      subroutine cgetsssa(psss)
      use cplmod
      implicit none
!
!     get interpolated oce sss from common 
!
      real, intent(out) :: psss(nxa,nya) ! sss (output)
      real :: zgw,zgm,zgs
!
!     get accumulated oce sss
!
      if(nssso == 0) then
       write(*,*) '!ERROR! no accumulated oce sss found'
       stop
      endif
      osss(:,:)=osssa(:,:)/real(nssso)
!
!     interpolate
!
      if (nissso == 1) then
       call cinter(osss,oslms,wos2a_a,inos2a,ios2a_a,asss,aslm          &
     &            ,ogws,agw,nxo,nxa,nyo,nya,ncssso,nimax,nprint)
      else
       call cinter2(osss,oslms,wos2a,ios2a,asss,aslm,ogws,agw           &
     &             ,nxo,nxa,nyo,nya,ncssso,nprint)
      endif
!
!     copy to output
!
      psss(:,:)=asss(:,:)
!
!     dbug printout
!
      if(nprint > 0) then
       write(*,*) ' '
       write(*,*) 'In CGETSSSA: '
       write(*,*) 'sss from ',nssso,' accumulated values'
       call cglm(osss,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
       write(*,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(osss,MASK=(oslms > 0.5))              &
     &            ,MINVAL(osss,MASK=(oslms > 0.5))
       call cglm(asss,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(*,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(asss,MASK=(aslm > 0.5))               &
     &            ,MINVAL(asss,MASK=(aslm > 0.5))
       write(*,*) ' '
      endif
!
!     reset accumulated field
!
      osssa(:,:)=0.
      nssso=0
!
      return
      end
!
!     ==================
!     subroutine CGETTSA
!     ==================
!
      subroutine cgettsa(pts)
      use cplmod
      implicit none
!
!     get interpolated oce tsurf from common 
!
      real :: pts(nxa,nya) ! tsurf (output)
!
      real :: zgs,zgm,zgw
!
!     get accumulated oce tsurf
!
      if(ntso == 0) then
       write(*,*) '!ERROR! no accumulated oce tsurf found'
       stop
      endif
      ots(:,:)=otsa(:,:)/real(ntso)
!
!     interpolate
!
      if (nitso == 1) then
       call cinter(ots,oslms,wos2a_a,inos2a,ios2a_a,ats,aslm            &
     &            ,ogws,agw,nxo,nxa,nyo,nya,nctso,nimax,nprint)
      else
       call cinter2(ots,oslms,wos2a,ios2a,ats,aslm,ogws,agw             &
     &             ,nxo,nxa,nyo,nya,nctso,nprint)
      endif
!
!     copy to output
!
      pts(:,:)=ats(:,:)
!
!     dbug printout
!
      if(nprint > 0) then
       write(*,*) ' '
       write(*,*) 'In CGETTSA: '
       write(*,*) 'tsurf from ',ntso,' accumulated values'
       call cglm(ots,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
       write(*,*) 'Oce (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(ots,MASK=(oslms > 0.5))               &
     &            ,MINVAL(ots,MASK=(oslms > 0.5))
       call cglm(ats,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(*,*) 'Atm (sum,mean,max,min): '                            &
     &            ,zgs,zgm,MAXVAL(ats,MASK=(aslm > 0.5))                &
     &            ,MINVAL(ats,MASK=(aslm > 0.5))
       write(*,*) ' '
      endif
!
!     reset accumulated field
!
      otsa(:,:)=0.
      ntso=0
!
      return
      end
!
!     ===================
!     subroutine CGETICEA
!     ===================
!
      subroutine cgeticea(picec,piced)
      use cplmod
      implicit none
!
!     get interpolated oce ice from common 
!
      real :: picec(nxa,nya) ! ice frac. (output)
      real :: piced(nxa,nya) ! ice thick. (output)
!
      real :: zgs,zgm,zgw
!
!     get accumulated oce sst
!
      if(niceo == 0) then
       write(*,*) '!ERROR! no accumulated oce ice found'
       stop
      endif
      oicec(:,:)=oiceca(:,:)/real(niceo)
      oiced(:,:)=oiceda(:,:)/real(niceo)
!
!     interpolate
!
      if (niiceo == 1) then
       call cinter(oicec,oslms,wos2a_a,inos2a,ios2a_a,aicec,aslm        &
     &            ,ogws,agw,nxo,nxa,nyo,nya,nciceo,nimax,nprint)
       call cinter(oiced,oslms,wos2a_a,inos2a,ios2a_a,aiced,aslm        &
     &            ,ogws,agw,nxo,nxa,nyo,nya,nciceo,nimax,nprint)
      else
       call cinter2(oicec,oslms,wos2a,ios2a,aicec,aslm,ogws,agw         &
     &             ,nxo,nxa,nyo,nya,nciceo,nprint)
       call cinter2(oiced,oslms,wos2a,ios2a,aiced,aslm,ogws,agw         &
     &             ,nxo,nxa,nyo,nya,nciceo,nprint)
      endif
!
!     copy to output
!
      picec(:,:)=aicec(:,:)
      piced(:,:)=aiced(:,:)
!
!     dbug printout
!
      if(nprint > 0) then
       write(*,*) ' '
       write(*,*) 'In CGETICEA: '
       write(*,*) 'ice from ',niceo,' accumulated values'
       call cglm(oicec,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
       write(*,*) 'Oce c. (sum,mean,max,min): '                         &
     &            ,zgs,zgm,MAXVAL(oicec,MASK=(oslms > 0.5))             &
     &            ,MINVAL(oicec,MASK=(oslms > 0.5))
       call cglm(aicec,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(*,*) 'Atm c. (sum,mean,max,min): '                         &
     &            ,zgs,zgm,MAXVAL(aicec,MASK=(aslm > 0.5))              &
     &            ,MINVAL(aicec,MASK=(aslm > 0.5))
       call cglm(oiced,oslms,zgs,zgm,zgw,ogws,nxo,nyo)
       write(*,*) 'Oce th. (sum,mean,max,min): '                        &
     &            ,zgs,zgm,MAXVAL(oiced,MASK=(oslms > 0.5))             &
     &            ,MINVAL(oiced,MASK=(oslms > 0.5))
       call cglm(aiced,aslm,zgs,zgm,zgw,agw,nxa,nya)
       write(*,*) 'Atm th. (sum,mean,max,min): '                        &
     &            ,zgs,zgm,MAXVAL(aiced,MASK=(aslm > 0.5))              &
     &            ,MINVAL(aiced,MASK=(aslm > 0.5))

       write(*,*) ' '
      endif
!
!     reset accumulated field
!
      oiceca(:,:)=0.
      oiceda(:,:)=0.
      niceo=0
!
      return
      end
!
!     =====================
!     subroutine CINTER_INI
!     =====================
!
      subroutine cinter_ini(pslmi,pfw,ifn,ifw,pxi1,pxi2,pxo1,pxo2       &
     &                     ,pyi1,pyi2,pyo1,pyo2,nxi,nxo,nyi,nyo,nwi,npr)
      implicit none
!
      real, parameter :: zerr=-9.E09
!
      integer, intent(in) :: nwi               ! max no of fields for inter.
      integer, intent(in) :: npr               ! print switch
      integer, intent(in) :: nxi               ! x dimension input
      integer, intent(in) :: nxo               ! x dimension output
      integer, intent(in) :: nyi               ! y dimension input
      integer, intent(in) :: nyo               ! y dimension input
      integer, intent(out):: ifn(nxo,nyo)      ! number of fields for inter.
      integer, intent(out):: ifw(nxo,nyo,nwi,2)! x and y indices for interpol.
!
      real, intent(in)  :: pslmi(nxi,nyi)      ! land sea mask of input
      real, intent(out) :: pfw(nxo,nyo,nwi)    ! weights
      real, intent(in)  :: pxi1(nxi,nyi)       ! 1st x border of input box 
      real, intent(in)  :: pxi2(nxi,nyi)       ! 2nd x border of input box
      real, intent(in)  :: pxo1(nxo,nyo)       ! 1st x border of output box
      real, intent(in)  :: pxo2(nxo,nyo)       ! 2nd x border of output box
      real, intent(in)  :: pyi1(nxi,nyi)       ! 1st y border of input box
      real, intent(in)  :: pyi2(nxi,nyi)       ! 2nd y border of input box
      real, intent(in)  :: pyo1(nxo,nyo)       ! 1st y border of output box
      real, intent(in)  :: pyo2(nxo,nyo)       ! 2nd y border of output box 
!
      real :: xi1(nxi,nyi)     ! left x's for input
      real :: xi2(nxi,nyi)     ! right x's for input
      real :: xo1(nxo,nyo)     ! left x's for output
      real :: xo2(nxo,nyo)     ! right x's for output
      real :: yi1(nxi,nyi)     ! southern y's for input
      real :: yi2(nxi,nyi)     ! northern y's for input
      real :: yo1(nxo,nyo)     ! southern y's for output
      real :: yo2(nxo,nyo)     ! northern y's for output
!
!
      integer :: ii,io,ji,jo,kmiss,iw,jw
      real :: pi
      real :: xil,xir,xol,xor
      real :: zgw,zweight
!
      pi=2.*asin(1.)
!
      if(pyi1(1,1) < pyi2(1,1)) then
       yi1(:,:)=pyi1(:,:)
       yi2(:,:)=pyi2(:,:)
      else
       yi1(:,:)=pyi2(:,:)
       yi2(:,:)=pyi1(:,:)
      endif
      xi1(:,:)=pxi1(:,:)
      xi2(:,:)=pxi2(:,:)
      if (pyo1(1,1) < pyo2(1,1)) then
       yo1(:,:)=pyo1(:,:)
       yo2(:,:)=pyo2(:,:)
      else
       yo1(:,:)=pyo2(:,:)
       yo2(:,:)=pyo1(:,:)
      endif
      xo1(:,:)=pxo1(:,:)
      xo2(:,:)=pxo2(:,:)
!
      pfw(:,:,:)=0.
      ifn(:,:)=0
      ifw(:,:,:,:)=0
      kmiss=0
      do jo=1,nyo
      do io=1,nxo
       iw=0
       xol=pxo1(io,jo)
       xor=pxo2(io,jo)
       do ji=1,nyi
       do ii=1,nxi
        if(pslmi(ii,ji) > 0.) then
         xil=pxi1(ii,ji)
         xir=pxi2(ii,ji)
         if(ABS(xol-xil) > pi) then
          if(xol > xil) then
           xor=xor-2.*pi
           xol=xol-2.*pi
          else
           xir=xir-2.*pi
           xil=xil-2.*pi
          endif
         endif
!
         if((xol <= xil .and. xor > xil)                                &
     &   .and. (yo1(io,jo) <= yi1(ii,ji) .and. yo2(io,jo) > yi1(ii,ji)))&
     &   then
          iw=iw+1
          if(iw > nwi) then
           print*,'ERROR: nwi too small!'
           stop
          endif
          zweight=(MIN(xor,xir)-xil)                                    &
     &           *ABS(sin(MIN(yo2(io,jo),yi2(ii,ji)))-sin(yi1(ii,ji)))
          pfw(io,jo,iw)=zweight
          ifw(io,jo,iw,1)=ii
          ifw(io,jo,iw,2)=ji
         endif
!
         if((xol > xil .and. xol < xir)                                 &
     &   .and. (yo1(io,jo) <= yi1(ii,ji) .and. yo2(io,jo) > yi1(ii,ji)))&
     &   then
          iw=iw+1
          if(iw > nwi) then
           print*,'ERROR: nwi too small!'
           stop
          endif
          zweight=(MIN(xor,xir)-xol)                                    &
     &           *ABS(sin(MIN(yo2(io,jo),yi2(ii,ji)))-sin(yi1(ii,ji)))
          pfw(io,jo,iw)=zweight
          ifw(io,jo,iw,1)=ii
          ifw(io,jo,iw,2)=ji
         endif
!
         if((xol <= xil .and. xor > xil)                                &
     &   .and. (yo1(io,jo) > yi1(ii,ji) .and. yo1(io,jo) < yi2(ii,ji))) &
     &   then
          iw=iw+1
          if(iw > nwi) then
           print*,'ERROR: nwi too small!'
           stop
          endif
          zweight=(MIN(xor,xir)-xil)                                    &
     &           *ABS(sin(MIN(yo2(io,jo),yi2(ii,ji)))-sin(yo1(io,jo)))
          pfw(io,jo,iw)=zweight
          ifw(io,jo,iw,1)=ii
          ifw(io,jo,iw,2)=ji
         endif
!
         if((xol > xil .and. xol < xir)                                 &
     &   .and. (yo1(io,jo) > yi1(ii,ji) .and. yo1(io,jo) < yi2(ii,ji))) &
     &   then
          iw=iw+1
          if(iw > nwi) then
           print*,'ERROR: nwi too small!'
           stop
          endif
          zweight=(MIN(xor,xir)-xol)                                    &
     &           *ABS(sin(MIN(yo2(io,jo),yi2(ii,ji)))-sin(yo1(io,jo)))
          pfw(io,jo,iw)=zweight 
          ifw(io,jo,iw,1)=ii
          ifw(io,jo,iw,2)=ji
         endif
        endif
       enddo
       enddo
       zgw=0.
       do jw=1,iw
        zgw=zgw+pfw(io,jo,jw)
       enddo
       if(zgw > 0.) then
        do jw=1,jw
         pfw(io,jo,jw)=pfw(io,jo,jw)/zgw  
        enddo
       endif
       ifn(io,jo)=iw
      enddo
      enddo
!
      return
      end subroutine cinter_ini

!
!     =================
!     subroutine CINTER
!     =================
!

      subroutine cinter(pfi,pslmi,pfw,ifn,ifw,pfo,pslmo,pgwi,pgwo       &
     &                  ,nxi,nxo,nyi,nyo,nc,nwi,npr)
!
      implicit none
      real, parameter :: zerr=-9.E09
!
      integer, intent(in) :: nc             ! switch for conservation method
      integer, intent(in) :: nwi            ! max. no of fields for inter.
      integer, intent(in) :: npr            ! switch for printout
      integer, intent(in) :: nxi            ! x-dimension input
      integer, intent(in) :: nxo            ! x-dimension output
      integer, intent(in) :: nyi            ! y-dimension input
      integer, intent(in) :: nyo            ! y-dimension output
!
      integer, intent(in) :: ifn(nxo,nyo)       ! no fields for inter. 
      integer, intent(in) :: ifw(nxo,nyo,nwi,2) ! x and y indices for interpol. 
      real, intent(in) :: pfi(nxi,nyi)          ! input field
      real, intent(in) :: pslmi(nxi,nyi)        ! input mask (0/1)
      real, intent(in) :: pfw(nxo,nyo,nwi)      ! interpolation weights
      real, intent(out):: pfo(nxo,nyo)          ! output field
      real, intent(in) :: pslmo(nxo,nyo)        ! output mask (0/1)
      real, intent(in) :: pgwi(nxi,nyi)         ! input weights
      real, intent(in) :: pgwo(nxo,nyo)         ! output weights
!
      real :: zgw
      real :: fo(nxo,nyo)         
!
      integer :: kmiss
      integer :: io,jo,ii,ji,jw
!
      kmiss=0
      do jo=1,nyo
      do io=1,nxo
       if(pslmo(io,jo) > 0.) then
        fo(io,jo)=0.
        zgw=0.
        do jw=1,ifn(io,jo)
         ii=ifw(io,jo,jw,1)
         ji=ifw(io,jo,jw,2)
         fo(io,jo)=fo(io,jo)+pfw(io,jo,jw)*pfi(ii,ji) 
         zgw=zgw+pfw(io,jo,jw)
        enddo
        if(zgw < 0.5) then
         fo(io,jo)=zerr
         kmiss=kmiss+1 
        endif
       else
        fo(io,jo)=0.
       endif
      enddo
      enddo
!
!!      if(kmiss > 0) then
!!         write(*,*) 'kmiss= ',kmiss
!!         stop "routine cfill not yet ready for adjoint" ! FastOpt, debug
! FastOpt, debug call cfill(fo,pslmo,nxo,nyo,zerr,npr)
!!      end if
!
      if(nc > 0) then
       call conserve(pfi,fo,pslmi,pslmo,pgwi,pgwo,nxi,nxo,nyi,nyo,nc    &
     &              ,npr)
      endif
!
      pfo(:,:)=fo(:,:)
!
      return
      end subroutine cinter

!
!     ======================
!     subroutine CINTER2_INI
!     ======================
!
      subroutine cinter2_ini(pslmi,pfw,ifw,pgwi,pxi,pxo                 &
     &                      ,pyi,pyo,nxi,nxo,nyi,nyo,npr)
!
      implicit none
      real :: pi
      integer, intent(in) :: npr            ! print switch
      integer, intent(in) :: nxi            ! x dimension input
      integer, intent(in) :: nxo            ! x dimension output
      integer, intent(in) :: nyi            ! y dimension input
      integer, intent(in) :: nyo            ! y dimension input
      integer, intent(out):: ifw(nxo,nyo,2) ! x and y indices for interpol.
!
      real, intent(in) :: pslmi(nxi,nyi)    ! input mask (0/1)
      real, intent(out):: pfw(nxo,nyo,4)    ! interpolation weights
      real, intent(in) :: pxi(nxi,nyi)      ! x's for input
      real, intent(in) :: pxo(nxo,nyo)      ! x's for output
      real, intent(in) :: pyi(nxi,nyi)      ! y's for input
      real, intent(in) :: pyo(nxo,nyo)      ! y's for output
      real, intent(in) :: pgwi(nxi,nyi)     ! input weights 
!
      real :: slmi(0:nxi,0:nyi+1)  ! input mask (0/1)
      real :: xi(0:nxi,0:nyi+1)    ! x's for input
      real :: yi(0:nxi,0:nyi+1)    ! y's for input
      real :: gwi(0:nxi,0:nyi+1)   ! input weights
!
      integer :: io,jo,ii,ji
      real :: zgw,zxo,zyo,zxi,zxil,zgwx1,zgwx11,zgwx21,zgwy,zgwy1,zgwy2
      real :: zgwx2,zgwx12,zgwx22
!
      pi=2.*asin(1.)
!
      slmi(1:nxi,1:nyi)=pslmi(:,:)
      slmi(0,1:nyi)=pslmi(nxi,:)
      slmi(:,0)=slmi(:,1)
      slmi(:,nyi+1)=slmi(:,nyi)
      gwi(1:nxi,1:nyi)=pgwi(:,:)
      gwi(0,1:nyi)=pgwi(nxi,:)
      gwi(:,0)=gwi(:,1)
      gwi(:,nyi+1)=gwi(:,nyi)
      yi(1:nxi,1:nyi)=pyi(:,:)
      yi(0,1:nyi)=pyi(nxi,:)
      if(pyi(1,1) < pyi(1,nyi)) then
       yi(:,0)=-pi*0.5
       yi(:,nyi+1)=pi*0.5
      else
       yi(:,0)=pi*0.5
       yi(:,nyi+1)=-pi*0.5
      endif
      xi(1:nxi,1:nyi)=pxi(:,:)
      xi(0,1:nyi)=pxi(nxi,:)-2.*pi
      xi(:,0)=xi(:,1)
      xi(:,nyi+1)=xi(:,nyi)
!
      pfw(:,:,:)=0.
!
      do jo=1,nyo
      do io=1,nxo
       zgw=0.
       zxo=pxo(io,jo)
       zyo=pyo(io,jo)
       do ji=1,nyi+1
       do ii=1,nxi
        zxi=xi(ii,ji)
        zxil=xi(ii-1,ji)
        if(zxil > zxi) zxil=zxil-2.*pi
        if(ABS(zxil-zxi) > pi) then
         if(zxil > zxi) then
          zxil=zxil-2.*pi
         else
          zxi=zxi-2.*pi
         endif
        endif
        if(ABS(zxo-zxi) > pi) then
         if(zxo > zxi) then
          zxo=zxo-2.*pi
         else
          zxi=zxi-2.*pi
          zxil=zxil-2.*pi
         endif
        endif
        if(((zyo <= yi(ii,ji) .and. zyo > yi(ii,ji-1))                  &
     &      .or.(zyo > yi(ii,ji) .and. zyo <= yi(ii,ji-1)))             &
     &     .and.(zxo <= zxi .and. zxo > zxil)) then
         zgwx11=abs(zxo-zxi)*slmi(ii-1,ji)
         zgwx21=abs(zxo-zxil)*slmi(ii,ji)
         zgwx1=zgwx11+zgwx21
         if(zgwx1 > 0.) then
          zgwx11=zgwx11/zgwx1
          zgwx21=zgwx21/zgwx1
          zgwy1=abs(zyo-yi(ii,ji-1))                                    &
     &         *(gwi(ii,ji)*slmi(ii,ji)+gwi(ii-1,ji)*slmi(ii-1,ji))
         else
          zgwy1=0.
         endif
         zgwx12=abs(zxo-zxi)*slmi(ii-1,ji-1)
         zgwx22=abs(zxo-zxil)*slmi(ii,ji-1)
         zgwx2=zgwx12+zgwx22
         if(zgwx2 > 0.) then
          zgwx12=zgwx12/zgwx2
          zgwx22=zgwx22/zgwx2
          zgwy2=abs(zyo-yi(ii,ji))                                      &
     &         *(gwi(ii,ji-1)*slmi(ii,ji-1)                             &
     &          +gwi(ii-1,ji-1)*slmi(ii-1,ji-1))
         else
          zgwy2=0.
         endif
         zgwy=zgwy1+zgwy2
         if(zgwy > 0.) then
          pfw(io,jo,1)=zgwy1*zgwx11/zgwy
          pfw(io,jo,2)=zgwy1*zgwx21/zgwy
          pfw(io,jo,3)=zgwy2*zgwx12/zgwy
          pfw(io,jo,4)=zgwy2*zgwx22/zgwy
         endif
         ifw(io,jo,1)=ii
         ifw(io,jo,2)=ji
        endif
       enddo
       enddo
      enddo
      enddo
!
      return
      end subroutine cinter2_ini
!
!     ==================
!     subroutine CINTER2
!     ==================
!
      subroutine cinter2(pfi,pslmi,pfw,ifw,pfo,pslmo,pgwi,pgwo          &
     &                  ,nxi,nxo,nyi,nyo,nc,npr)
!
      implicit none
      real, parameter :: zerr=-9.E09
!
      integer, intent(in) :: nc             ! switch for conservation method
      integer, intent(in) :: npr            ! switch for printout
      integer, intent(in) :: nxi            ! x-dimension input
      integer, intent(in) :: nxo            ! x-dimension output
      integer, intent(in) :: nyi            ! y-dimension input
      integer, intent(in) :: nyo            ! y-dimension output
!
      integer, intent(in) :: ifw(nxo,nyo,2) ! x and y indices for interpol. 
      real, intent(in) :: pfi(nxi,nyi)      ! input field
      real, intent(in) :: pslmi(nxi,nyi)    ! input mask (0/1)
      real, intent(in) :: pfw(nxo,nyo,4)    ! interpolation weights
      real, intent(out):: pfo(nxo,nyo)      ! output field
      real, intent(in) :: pslmo(nxo,nyo)    ! output mask (0/1)
      real, intent(in) :: pgwi(nxi,nyi)     ! input weights
      real, intent(in) :: pgwo(nxo,nyo)     ! output weights
!
      real :: fi(0:nxi,0:nyi+1)    ! input field
      real :: fo(nxo,nyo)          ! output field
!
      integer :: kmiss
      integer :: io,jo,ii,ji
!
      fi(1:nxi,1:nyi)=pfi(:,:)
      fi(0,1:nyi)=pfi(nxi,:)
      fi(:,0)=fi(:,1)
      fi(:,nyi+1)=fi(:,nyi)
!
      kmiss=0
      do jo=1,nyo
      do io=1,nxo
       if(pslmo(io,jo) > 0.) then
        ii=ifw(io,jo,1)
        ji=ifw(io,jo,2)
        fo(io,jo)=pfw(io,jo,1)*fi(ii-1,ji)+pfw(io,jo,2)*fi(ii,ji)       &
     &           +pfw(io,jo,3)*fi(ii-1,ji-1)+pfw(io,jo,4)*fi(ii,ji-1)
        if(pfw(io,jo,1)+pfw(io,jo,2)+pfw(io,jo,3)+pfw(io,jo,4)<0.5) then
         fo(io,jo)=zerr
         kmiss=kmiss+1 
        endif
       else
        fo(io,jo)=0.
       endif
      enddo
      enddo
!
!
!
!!      if(kmiss > 0) then
!!         write(*,*) 'kmiss= ',kmiss
!!         stop "routine cfill not yet ready for adjoint" ! FastOpt, debug
! FastOpt, debug call cfill(fo,pslmo,nxo,nyo,zerr,npr)
!!      end if
!
      if(nc > 0) then
       call conserve(pfi,fo,pslmi,pslmo,pgwi,pgwo,nxi,nxo,nyi,nyo,nc    &
     &              ,npr)
      endif
!
      pfo(:,:)=fo(:,:)
!
      return
      end subroutine cinter2

!     =================
!     subroutine CEXTRA
!     =================
!
      subroutine cextra(pw,kn,kw,pm,nx,ny,nw,npr)
!
!     complete the weights by a simple extrapolation
!     (area overlap)
!
      implicit none
      integer, intent(in) :: nx,ny,nw,npr
      real, intent(inout) :: pw(nx,ny,nw)
      integer, intent(inout) :: kn(nx,ny)
      integer, intent(inout) :: kw(nx,ny,nw,2)
      real, intent(in) :: pm(nx,ny)
      integer :: kmiss,kiter,kimax,i,j,ii,io,jj,jo,i1,i2,j1,j2
!
      kimax=nx*ny
!
      do j=1,ny
      do i=1,nx
       if(SUM(pw(i,j,:)) < 0.5 .and. pm(i,j) > 0.5) then
        kmiss=1
        kiter=1
        jo=-1
        io=0
        do while(kmiss > 0 .and. kiter < kimax)
         jo=jo+1
         io=io+1
         j1=max(1,j-jo)
         j2=min(j+jo,ny)
         i1=max(1,i-io)
         i2=min(i+io,nx) 
         do jj=j1,j2
         do ii=i1,i2
          if(sum(pw(ii,jj,:)) > 0.5 .and. kmiss > 0) then
           pw(i,j,:)=pw(ii,jj,:)
           kn(i,j)=kn(ii,jj)
           kw(i,j,:,:)=kw(ii,jj,:,:)
           kmiss=0
           if(npr > 0) then 
            write(*,*) 'found ii,jj= ',ii,jj,' for i,j= ',i,j           &
     &                ,' after iter= ',kiter
           endif
          endif
         enddo
         enddo
         kiter=kiter+1
        enddo
       endif
      enddo
      enddo  
!
!     check again
!
      kmiss=0
      do j=1,ny
      do i=1,nx
       if(SUM(pw(i,j,:)) < 0.5 .and. pm(i,j) > 0.5) kmiss=kmiss+1
      enddo
      enddo
      if(kmiss > 0) then
       write(*,*) 'In CEXTRA: ERROR! still ',kmiss,' points found'
       stop
      else
       if(npr > 0) then
        write(*,*) 'In CEXTRA: No (more) singular points found' 
       endif
      endif
!
      return
      end subroutine cextra

!     ==================
!     subroutine CEXTRA2
!     ==================
!
      subroutine cextra2(pw,kw,pm,nx,ny,nw,npr)
!
!     complete the weights by a simple extrapolation
!     (bi linear interpolation)
!
      implicit none
      integer, intent(in) :: nx,ny,nw,npr
      real, intent(inout) :: pw(nx,ny,nw)
      integer, intent(inout) :: kw(nx,ny,2)
      real, intent(in) :: pm(nx,ny)
      integer :: kmiss,kiter,kimax,i,j,ii,io,jj,jo,i1,i2,j1,j2
!
      kimax=nx*ny
!
      do j=1,ny
      do i=1,nx
       if(SUM(pw(i,j,:)) < 0.5 .and. pm(i,j) > 0.5) then
        kmiss=1
        kiter=1
        jo=-1
        io=0
        do while(kmiss > 0 .and. kiter < kimax)
         jo=jo+1
         io=io+1
         j1=max(1,j-jo)
         j2=min(j+jo,ny)
         i1=max(1,i-io)
         i2=min(i+io,nx)
         do jj=j1,j2
         do ii=i1,i2
          if(sum(pw(ii,jj,:)) > 0.5 .and. kmiss > 0) then
           pw(i,j,:)=pw(ii,jj,:)
           kw(i,j,:)=kw(ii,jj,:)
           kmiss=0
           if(npr > 0) then
            write(*,*) 'found ii,jj= ',ii,jj,' for i,j= ',i,j           &
     &                ,' after iter= ',kiter
           endif
          endif
         enddo
         enddo
         kiter=kiter+1
        enddo
       endif
      enddo
      enddo
!
!     check again
!
      kmiss=0
      do j=1,ny
      do i=1,nx
       if(SUM(pw(i,j,:)) < 0.5 .and. pm(i,j) > 0.5) kmiss=kmiss+1
      enddo
      enddo
      if(kmiss > 0) then
       write(*,*) 'In CEXTRA2: ERROR! still ',kmiss,' points found'
       stop
      else
       if(npr > 0) then
        write(*,*) 'In CEXTRA2: No (more) singular points found' 
       endif 
      endif
!
      return
      end subroutine cextra2
!
!     ===================
!     subroutine CONSERVE
!     ===================
!
      subroutine conserve(pfi,pfo,pmi,pmo,pgwi,pgwo                     &
     &                   ,nxi,nxo,nyi,nyo,nc,npr)
      implicit none
!
!     correct global mean (multiplicative or additive depending on nc)
!
!     nc=1 : multiply f to give sum(fi*gwi)=sum(fon*gwo) with fon=f*fo
!     nc=2 : multiply f to give sum(fi*gwi)/sum(gwi)=sum(fon*gwo)/sum(gwo)
!     nc=3 : add a to give sum(fi*gwi)=sum(fon*gwo) with fon=fo+a
!     nc=4 : add a to give sum(fi*gwi)/sum(gwi)=sum(fon*gwo)/sum(gwo)
!
      integer, intent(in) :: nxi,nxo,nyi,nyo,nc,npr
      real :: pfi(nxi,nyi)
      real :: pfo(nxo,nyo)
      real :: pmi(nxi,nyi)
      real :: pmo(nxo,nyo)
      real :: pgwi(nxi,nyi)
      real :: pgwo(nxo,nyo)
!
      integer :: i,j
      real :: zadd1,zadd2,zfac1,zfac2,zfi,zfo,zgw,zgwi,zgwo
!
      zgwi=0.
      zfi=0.
      do j=1,nyi
      do i=1,nxi
       zgw=pgwi(i,j)*pmi(i,j)
       zfi=zfi+pfi(i,j)*zgw
       zgwi=zgwi+zgw
      enddo
      enddo
      zgwo=0.
      zfo=0.
      do j=1,nyo
      do i=1,nxo
       zgw=pgwo(i,j)*pmo(i,j)
       zfo=zfo+pfo(i,j)*zgw
       zgwo=zgwo+zgw
      enddo
      enddo
!
      if(nc < 3) then
       if(zfo == 0.) then
        if(zfi /= 0.) write(*,*) 'WARNING (CPL): No Conservation (zfo=0.)'
        zfac1=1.
        zfac2=1.
       else
        zfac1=zfi/zfo
        zfac2=zfi/zfo*zgwo/zgwi
       endif
      endif
      zadd1=(zfi-zfo)/zgwo
      zadd2=zfi/zgwi-zfo/zgwo
!
      if(nc==1) pfo(:,:)=pfo(:,:)*zfac1
      if(nc==2) pfo(:,:)=pfo(:,:)*zfac2
      if(nc==3) pfo(:,:)=pfo(:,:)+zadd1*pmo(:,:)
      if(nc==4) pfo(:,:)=pfo(:,:)+zadd2*pmo(:,:)
!
      if(npr > 0) then
       if(nc==1) write(*,*) 'In CONSERVE: nc= ',nc,' Fac= ',zfac1
       if(nc==2) write(*,*) 'In CONSERVE: nc= ',nc,' Fac= ',zfac2
       if(nc==3) write(*,*) 'In CONSERVE: nc= ',nc,' Add= ',zadd1
       if(nc==4) write(*,*) 'In CONSERVE: nc= ',nc,' Add= ',zadd2
      endif
!
      return
      end
!
!     ===============
!     subroutine CGLM
!     ===============
!
      subroutine cglm(pf,pm,pgs,pgm,pgw,pgwi,nx,ny)
      implicit none
!
!     compute global means
!
      integer, intent(in) :: nx,ny
      real :: pf(nx,ny)
      real, intent(out) :: pgs
      real :: pgm
      real, intent(out) :: pgw
      real :: pm(nx,ny)
      real :: pgwi(nx,ny)
!
      real :: zgwd,zgsd,zgw
      integer :: i,j
!      
      zgwd=0.
      zgsd=0.
      do j=1,ny
      do i=1,nx
       zgw=pgwi(i,j)*pm(i,j)
       zgsd=zgsd+pf(i,j)*zgw
       zgwd=zgwd+zgw
      enddo
      enddo
      if(zgwd > 0.) then
       pgm=zgsd/zgwd
      else
       print*,'!WARNING! no sea points found in cglm2'
      endif
      pgs=zgsd
      pgw=zgwd
!
      return
      end subroutine cglm
!
!     ================
!     subroutine CPFAC
!     ================
!
      subroutine cpfac(pfsst,pftau,pffresh,pfice)
      use cplmod
      implicit none
!
      real :: pfsst,pftau,pffresh,pfice
!
      pfsst=cfacsst 
      pftau=cfactau
      pffresh=cfacfresh
      pfice=cfacice
!
      return
      end subroutine cpfac
!
!     =================
!     subroutine CROFFC
!     =================
!
      subroutine croffc(proff)
      use cplmod
      implicit none
!
      real :: proff(nxa,nya)
      real :: zroff(nxa,nya)
!
      integer :: ja,kroffa
      real :: zgm,zgs,zgw,zgwa,zroffa
!
      if(nroffc == 1) then
!
!     substitute lokal runoff by global mean
!
       call cglm(proff,aslm,zgs,zgm,zgw,agw,nxa,nya)
       proff(:,:)=zgm*aslm(:,:)
      endif
!
      if(nroffc == 2) then
!
!     substitute lokal runoff by area averages
!     where iroffi==0 local runoff is kept
!
!!     kroff=0
       kroffa=maxval(iroffi(:,:))
       do ja=1,kroffa
        zroffa=SUM(proff(:,:)*agw(:,:)*aslm(:,:)                        &
     &            ,mask=(ja == iroffi(:,:)))
        zgwa=SUM(agw(:,:)*aslm(:,:),mask=(ja == iroffi(:,:))) 
!!      kroff=kroff+SUM(nint(aslm(:,:)),mask=(ja == iroffi(:,:)))
        if(zgwa <= 0.) then
         zgwa=-1.
         zroffa=0.
        endif
        where(ja == iroffi(:,:))
         proff(:,:)=zroffa/zgwa*aslm(:,:)
        endwhere      
       enddo
!!     klsm=SUM(nint(aslm(:,:)))
!!     if(kroff /= klsm) then
!!      print*,'WARNING: runoff corrected only for ',kroff,' from '     &
!!   &        ,klsm,' gridpoints!'
!!     endif
      elseif(nroffc == 3) then
       zroff(:,:)=0.
!
!     redistribute lokal runoff using area averages
!     where iroffi==0 local runoff + area average is used
!
       where(iroffi(:,:) == 0) 
        zroff(:,:)=proff(:,:)
       endwhere
       kroffa=maxval(iroffi(:,:))
!!     kroffi=0
       do ja=1,kroffa
        zroffa=SUM(proff(:,:)*agw(:,:)*aslm(:,:)                        &
     &            ,mask=(ja == iroffi(:,:)))
        zgwa=SUM(agw(:,:)*aslm(:,:),mask=(ja == iroffo(:,:))) 
!!      kroffi=kroffi+SUM(nint(aslm(:,:)),mask=(ja == iroffi(:,:)))
        if(zgwa <= 0. .and. zroffa > 0.) then
         print*,'ERROR: zero area for runoff output no ',ja
         stop
        endif
        if(zgwa > 0.) then
         where(ja == iroffo(:,:))
          zroff(:,:)=zroffa/zgwa*aslm(:,:)+zroff(:,:)
         endwhere    
        endif  
       enddo
!!     klsm=SUM(nint(aslm(:,:)))
!!     if(kroffi /= klsm) then
!!      print*,'ERROR: not all runoff redistributed!',klsm,kroffi,kroffa 
!!      stop    
!!     endif
       proff(:,:)=zroff(:,:)
      endif
!
      return
      end subroutine croffc

!     =================
!     SUBROUTINE CGAUAW
!     =================

      subroutine cgauaw(psi,pgw,k)
      implicit none
      real psi(k),pgw(k)
      integer, intent(in) :: k
      real (kind=8) c,eps
      parameter(c = (1.0 - (2.0/3.14159265358979)**2)/4.0)
      parameter(eps = 1.0e-14)
      real (kind=8) pa(k),pw(k)
      real (kind=8) pk,pkm1,pkm2,pkmrk,zsp,xz
      integer :: is, iter, kk, n

      kk = k/2
      call cbsslzr(pa,kk)

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
      end subroutine cgauaw

!     ==================
!     SUBROUTINE CBSSLZR
!     ==================

      subroutine cbsslzr(pbes,kdim)
      implicit none
      integer, intent(in) :: kdim
      integer, parameter :: NDIM=50
      real (kind=8) pbes(kdim)
      real (kind=8) zbes(NDIM)
      real (kind=8) PI
      integer :: idim, j
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
      end subroutine cbsslzr
!
!     ====================
!     subroutine MKNEWMASK
!     ====================
!
      subroutine mknewmask
      use cplmod
      implicit none
!
      real :: zaslm(nxa,nya) = 1.
      real :: zoslm(nxo,nyo) = 1.
      real :: zaslmn(nxa,nya)= 0.
      real :: zaslmf(nxa,nya)= 0.
      real :: zduma(nxa,nya) = 1.
      real :: zdumo(nxo,nyo) = 1.

      integer :: kwa(nxa,nya,2) = 0 
      integer :: kwo(nxo,nyo,2) = 0
      integer :: kwo2(nxo,nyo,2) = 0
      integer :: kwna(nxa,nya) = 0
      integer :: kwno(nxo,nyo) = 0
      integer :: kwa_a(nxa,nya,nimax,2) =0
      integer :: kwo_a(nxo,nyo,nimax,2) =0

!
      real :: zwa(nxa,nya,4) = 0.  
      real :: zwo(nxo,nyo,4) = 0.
      real :: zwo2(nxo,nyo,4) = 0.
      real :: zwa_a(nxa,nya,nimax) = 0.
      real :: zwo_a(nxo,nyo,nimax) = 0.

      integer :: i,j,jt,ii,jj,iim,jjm,ntape,kmissa,kmisso
      real :: zth,zzth,zref,zdiff,zgs,zgm,zgw,zpi
      logical :: lopen
      character(len=8) :: cform
!
      zpi=2.*asin(1.)
!
      write(*,*) '--------------------------------------'
      write(*,*) 'Compute a new aslm for a better match:'
!
      call cinter2_ini(zaslm,zwo,kwo,agw,xa,xos                         &
     &                ,ya,yos,nxa,nxo,nya,nyo,0)
!
!     1) interpolate oslm to atm grid
!
      call cinter2_ini(zoslm,zwa,kwa,ogws,xos,xa                        &
     &                ,yos,ya,nxo,nxa,nyo,nya,0)
      call cinter2(oslms,zoslm,zwa,kwa,zaslmn,zaslm,ogws,agw            &
     &            ,nxo,nxa,nyo,nya,0,0)
!
!!
!!      call cinter_ini(zoslm,zwa_a,kwna,kwa_a,xos1,xos2,xa1,xa2          &
!!     &               ,yos1,yos2,ya1,ya2,nxo,nxa,nyo,nya,nimax,0)
!!      call cinter(oslms,zoslm,zwa_a,kwna,kwa_a,zaslmn,zaslm             &
!!     &           ,ogws,agw,nxo,nxa,nyo,nya,0,nimax,0)
!
      call cglm(oslms,zdumo,zgs,zgm,zgw,ogws,nxo,nyo)
      write(*,*) 'Oce: ocean area scalar (sum/mean) + global area : '   &
     &           ,zgs*radeo*radeo*4.*zpi,zgm,zgw*radeo*radeo*4.*zpi
      zref=zgs*radeo*radeo
!
      call cglm(aslm,zduma,zgs,zgm,zgw,agw,nxa,nya)
      write(*,*) 'Old atm ocean area (sum/mean) + global area :     '   &
     &           ,zgs*radea*radea*4.*zpi,zgm,zgw*radea*radea*4.*zpi
!
      kmisso=0
      do j=1,nyo
      do i=1,nxo
       if(SUM(wa2os(i,j,:)) < 0.5 .and. oslms(i,j) > 0.5) then
        kmisso=kmisso+1
       endif
      enddo
      enddo
      kmissa=0
      do j=1,nya
      do i=1,nxa
       if(SUM(wos2a(i,j,:)) < 0.5 .and. aslm(i,j) > 0.5) then
        kmissa=kmissa+1
       endif
      enddo
      enddo
      write(*,*) 'Old: ',kmisso,' ocean and ',kmissa                    &
     &          ,' atm points missed' 
!
!     2) fit the global ocean area
!
      zaslm(:,:)=0.
      where(zaslmn(:,:) > 0.) zaslm(:,:)=1.
      call cglm(zaslm,zduma,zgs,zgm,zgw,agw,nxa,nya)
      zdiff=abs(zref-zgs*radea*radea)
      zaslmf(:,:)=zaslm(:,:)
      zzth=0.
      do jt=1,9
       zth=0.1*real(jt)
       zaslm(:,:)=0.
       where(zaslmn(:,:) > zth) zaslm(:,:)=1.
       call cglm(zaslm,zduma,zgs,zgm,zgw,agw,nxa,nya)
       if(abs(zref-zgs*radea*radea) < zdiff) then
        zdiff=abs(zref-zgs*radea*radea)
        zaslmf(:,:)=zaslm(:,:)
        zzth=zth
       endif
      enddo
      write(*,*) 'Best fit with threshold : ',zzth
!
!     3) fix missing oceanic points (as many as possible)
!     a) scalar
!
      call cinter2_ini(zduma,zwo,kwo,agw,xa,xos                         &
     &                ,ya,yos,nxa,nxo,nya,nyo,0)
      call cinter2_ini(zaslmf,zwo2,kwo2,agw,xa,xos                      &
     &                ,ya,yos,nxa,nxo,nya,nyo,0)
      do j=1,nyo
      do i=1,nxo
       if(SUM(zwo2(i,j,:)) < 0.5 .and. oslms(i,j) > 0.5) then
        ii=kwo(i,j,1)
        jj=kwo(i,j,2)
        iim=ii-1
        if(iim < 1) iim=nxa
        jjm=jj-1
        if(jjm < 1) jjm=jj
        if(zaslmn(ii,jj) > 0.) then
         zaslmf(ii,jj)=1. 
        elseif(zaslmn(iim,jj) > 0.) then     
         zaslmf(iim,jj)=1.  
        elseif(zaslmn(ii,jjm) > 0.) then
         zaslmf(ii,jjm)=1.
        elseif(zaslmn(iim,jjm) > 0.) then
         zaslmf(iim,jjm)=1.
        else
         write(*,*) 'Can not find atm for s ',i,j
        endif
       endif
      enddo
      enddo
!
!     b) u
!
      call cinter2_ini(zduma,zwo,kwo,agw,xa,xou                         &
     &                ,ya,you,nxa,nxo,nya,nyo,0)
      call cinter2_ini(zaslmf,zwo2,kwo2,agw,xa,xou                      &
     &                ,ya,you,nxa,nxo,nya,nyo,0)
      do j=1,nyo
      do i=1,nxo
       if(SUM(zwo2(i,j,:)) < 0.5 .and. oslmu(i,j) > 0.5) then
        ii=kwo(i,j,1)
        jj=kwo(i,j,2)
        iim=ii-1
        if(iim < 1) iim=nxa
        jjm=jj-1
        if(jjm < 1) jjm=jj
        if(zaslmn(ii,jj) > 0.) then
         zaslmf(ii,jj)=1.
        elseif(zaslmn(iim,jj) > 0.) then
         zaslmf(iim,jj)=1.
        elseif(zaslmn(ii,jjm) > 0.) then
         zaslmf(ii,jjm)=1.
        elseif(zaslmn(iim,jjm) > 0.) then
         zaslmf(iim,jjm)=1.
        else
         write(*,*) 'Can not find atm for u',i,j
        endif
       endif
      enddo
      enddo
!
!     c) v
!
      call cinter2_ini(zduma,zwo,kwo,agw,xa,xov                         &
     &                ,ya,yov,nxa,nxo,nya,nyo,0)
      call cinter2_ini(zaslmf,zwo2,kwo2,agw,xa,xov                      &
     &                ,ya,yov,nxa,nxo,nya,nyo,0)
      do j=1,nyo
      do i=1,nxo
       if(SUM(zwo2(i,j,:)) < 0.5 .and. oslmv(i,j) > 0.5) then
        ii=kwo(i,j,1)
        jj=kwo(i,j,2)
        iim=ii-1
        if(iim < 1) iim=nxa
        jjm=jj-1
        if(jjm < 1) jjm=jj
        if(zaslmn(ii,jj) > 0.) then
         zaslmf(ii,jj)=1.
        elseif(zaslmn(iim,jj) > 0.) then
         zaslmf(iim,jj)=1.
        elseif(zaslmn(ii,jjm) > 0.) then
         zaslmf(ii,jjm)=1.
        elseif(zaslmn(iim,jjm) > 0.) then
         zaslmf(iim,jjm)=1.
        else
         write(*,*) 'Can not find atm for v ',i,j
        endif
       endif
      enddo
      enddo
!
      call cglm(zaslmf,zduma,zgs,zgm,zgw,agw,nxa,nya)
      write(*,*) 'New atm ocean area (sum/mean) + global area :     '   &
     &           ,zgs*radea*radea*4.*zpi,zgm,zgw*radea*radea*4.*zpi
!
      call cinter2_ini(zaslmf,zwo,kwo,agw,xa,xos                        &
     &                ,ya,yos,nxa,nxo,nya,nyo,0)
      kmisso=0
      do j=1,nyo
      do i=1,nxo
       if(SUM(zwo(i,j,:)) < 0.5 .and. oslms(i,j) > 0.5) then
        kmisso=kmisso+1
       endif
      enddo
      enddo
      kmissa=0
      do j=1,nya
      do i=1,nxa
       if(SUM(wos2a(i,j,:)) < 0.5 .and. zaslmf(i,j) > 0.5) then
        kmissa=kmissa+1
       endif
      enddo
      enddo
      write(*,*) 'New: ',kmisso,' ocean s and ',kmissa                  &
     &          ,' atm points missed'
!
      call cinter2_ini(zaslmf,zwo,kwo,agw,xa,xou                        &
     &                ,ya,you,nxa,nxo,nya,nyo,0)
      kmisso=0
      do j=1,nyo
      do i=1,nxo
       if(SUM(zwo(i,j,:)) < 0.5 .and. oslmu(i,j) > 0.5) then
        kmisso=kmisso+1
       endif
      enddo
      enddo
      call cinter2_ini(zaslmf,zwo,kwo,agw,xa,xov                        &
     &                ,ya,yov,nxa,nxo,nya,nyo,0)
      kmissa=0
      do j=1,nyo
      do i=1,nxo
       if(SUM(zwo(i,j,:)) < 0.5 .and. oslmv(i,j) > 0.5) then
        kmissa=kmissa+1
       endif
      enddo
      enddo
      write(*,*) 'New: ',kmisso,' ocean u and ',kmissa                  &
     &          ,' ocean v points missed'
!
      write(*,*) 'New alsm written to newmask.txt'
      do jt=30,99
       ntape=jt
       inquire(unit=ntape,opened=lopen)
       if(.not. lopen) exit
      enddo
      zaslmf(:,:)=real(nint(1.-zaslmf(:,:)))
      open(ntape,file='newmask.txt',form='formatted')
      write(ntape,'(8I10)') 172,0,0,0,nxa,nya,0,0
      write(ntape,'(8E12.6)') zaslmf(:,:)
      close(ntape)
      write(*,*) '--------------------------------------'
!
      return
      end subroutine mknewmask



