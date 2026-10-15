subroutine write_internal_restart(yfile)
! author: FastOpt GmbH
! modified by larry 03.01.2013
! purpose: write a snapshot of the model state to 'yfile'
! input:
!  character (len=*) :: yfile
!
  use restartmod
  use landmod, only : pgrow,plai,pgs,pz0_max
  use pumamod
  implicit none
  character (len=*) :: yfile
!  rewind (nwriunit)
! the following copied from s/r epilog
  if (mypid == NROOT) then
   call restart_prepare(yfile)
   call put_restart_integer('nstep'   ,nstep   )
   call put_restart_integer('naccuout',naccuout)
   call put_restart_integer('nlat'    ,NLAT    )
   call put_restart_integer('nlon'    ,NLON    )
   call put_restart_integer('nlev'    ,NLEV    )
   call put_restart_integer('nrsp'    ,NRSP    )

   call put_restart_integer('nkits'   ,NKITS   )

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
  !*    finish surface parameterizations or models
  !
 
  call surfstop ! writes further restart variables
  
! important parameter fields
  call mpputgp('pgrow'  ,pgrow  ,NHOR,1)! added, FastOpt
  call mpputgp('plai'   ,plai   ,NHOR,1)! added, FastOpt
  call mpputgp('pgs'    ,pgs    ,NHOR,1)! added, FastOpt
  call mpputgp('pz0_max',pz0_max,NHOR,1)! added, FastOpt
  
  call restart_stop
  return
end subroutine write_internal_restart

subroutine read_internal_restart(yfile)
! author: FastOpt GmbH
! modified by larry 03.01.2013
! purpose: read a snapshot of the model state from 'yfile'
! input:
!  character (len=*) :: yfile
!
  use restartmod
  use pumamod
  use radmod, only : dqo3,gdist2,gmu0
  use rainmod, only : icclev,icctot
  use surfmod, only : nlastsu
  use landmod, only : dtclim, dwater, dwclim
  use landmod, only : pgrow,plai,pgs,pz0_max
  use icemod, only:
  use oceanmod, only:
  use vegmod, only : dvalb,dvrhs,dvsoil,dvz0,zfw,zveg
  implicit none
  character (len=*) :: yfile
  logical :: lyfile, lopen
  lunstable = .false.
  lyfile=.false.
  if(mypid == NROOT) then
   inquire(unit=40,opened=lopen) ! puma output
   if (lopen) then
     close(40)
     call outini ! header must be re-written
   end if
   inquire(unit=31,opened=lopen) ! ocean output
   if (lopen) then
     rewind(31)
   else
     open(31,file='ocean_output',form='unformatted')
   endif
   inquire(unit=71,opened=lopen) ! ice output
   if (lopen) then
     rewind(71)
   else
     open(71,file='ice_output',form='unformatted')
   end if
   call restart_stop
   nresnum=0
   call restart_ini(lyfile,yfile)
   if (.not.lyfile) then
     print*,"Cannot find internal restart file ",trim(yfile)
     stop
   end if
  endif
! interface pre-initiatlisation  
  du = 0.  
  dv = 0.
  dtsa = 0.
  dsnow = 0.
  dqt = 0.
  dswfl(:,:) = 0.
  dlwfl(:,:) = 0.
  dflux(:,:) = 0.
!
  call read_atmos_restart
  if (mypid == NROOT) then
   call get_restart_integer('nkits'   ,NKITS   ) ! added
  endif
  call mpbci(nkits) 
! 
  dprl(:)=0. ! added
  dprc(:)=0. ! added
  dprs(:)=0. ! added

  dftu(:,:)  = 0.0 ! added

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
! reset tendencies
!
  szt(:,:) = 0.
  sdt(:,:) = 0.
  stt(:,:) = 0.
  sqt(:,:) = 0.
  spt(:) = 0.
!
! reset gridpoint fields
!
  gt(:,:) = 0.0
  gq(:,:) = 0.0
  gd(:,:) = 0.0
  gz(:,:) = 0.0
  gu(:,:) = 0.0
  gv(:,:) = 0.0
  gp(:)   = 0.0
  gpj(:)  = 0.0
!
! compute required fields
!
  call sp2gp(sp(:),gp(:),1)
  dp(:) = psurf*exp(gp(:))
! radmod pre-initialisation
  dqo3 = 0.
  gdist2 = 0.
  gmu0 = 0.
! rainmod pre-initialisation
  icclev = 0.
  icctot = 0.
!
!  surface restart
!
! surfmod pre-initialisation
  nlastsu=0.
! landmod pre-initialisation
  dtclim=0.
  dwater=0.
  dwclim=0.
  call read_internal_restart_landmod
! no restart info from soilini
! no restart info from roffini
! vegmod pre-initialisation
  dvalb = 0.
  dvrhs = 0.
  dvsoil = 0.
  dvz0 = 0.
  zfw = 0.
  zveg = 0.
! no restart info from vegini

  call read_internal_restart_seamod

! important parameter fields:
  call mpgetgp('pgrow'  ,pgrow  ,NHOR,1)! added, FastOpt
  call mpgetgp('plai'   ,plai   ,NHOR,1)! added, FastOpt
  call mpgetgp('pgs'    ,pgs    ,NHOR,1)! added, FastOpt
  call mpgetgp('pz0_max',pz0_max,NHOR,1)! added, FastOpt

  if(mypid==NROOT) then
   call restart_stop
  endif
  
  return
end subroutine read_internal_restart

subroutine read_internal_restart_landmod
! purpose: a version of landmod::landini which can be called more than once
!          to restore a model snapshot
!
  use restartmod
  use landmod
  implicit none
  integer :: ncountsea
  external ncountsea
!
! from landini (landmod)
!
  call mpgetgp('dtsl'    ,dts     ,NHOR,     1)
  call mpgetgp('dtsm'    ,dtsm    ,NHOR,     1)
  call mpgetgp('dqs'     ,dqs     ,NHOR,     1)
  call mpgetgp('driver'  ,driver  ,NHOR,     1)
  call mpgetgp('duroff'  ,duroff  ,NHOR,     1)
  call mpgetgp('dvroff'  ,dvroff  ,NHOR,     1)
  call mpgetgp('darea'   ,darea   ,NHOR,     1)
  call mpgetgp('dwmax'   ,dwmax   ,NHOR,     1)
  call mpgetgp('dtcl'    ,dtcl    ,NHOR,    14)
  call mpgetgp('dwcl'    ,dwcl    ,NHOR,    14)
  call mpgetgp('dsnowt'  ,dsnowt  ,NHOR,     1)
  call mpgetgp('dsnowz'  ,dsnowz  ,NHOR,     1)
  call mpgetgp('dsoilt'  ,dsoilt  ,NHOR,NLSOIL)
  call mpgetgp('dcveg'   ,dcveg   ,NHOR,     1)
  call mpgetgp('dcsoil'  ,dcsoil  ,NHOR,     1)
  call mpgetgp('dglac'   ,dglac   ,NHOR,     1)
  call mpgetgp('dz0clim' ,dz0clim ,NHOR,     1)
  call mpgetgp('dz0climo',dz0climo,NHOR,     1)
  call mpgetgp('dalbclim',dalbclim,NHOR,     1)
!
!     accumulated simba diagnostic
!
  call mpgetgp('agpp'   ,agpp   ,NHOR,1)
  call mpgetgp('anpp'   ,anpp   ,NHOR,1)
  call mpgetgp('agppl'  ,agppl  ,NHOR,1)
  call mpgetgp('agppw'  ,agppw  ,NHOR,1)
  call mpgetgp('anogrow',anogrow,NHOR,1)
  call mpgetgp('aresh'  ,aresh  ,NHOR,1)
  call mpgetgp('alitter',alitter,NHOR,1)
  
!      accumulated dune diagnostic

  call mpgetgp('awloss' ,awloss ,NHOR,1)
  
end subroutine read_internal_restart_landmod

subroutine read_internal_restart_icemod(pts,psst,pmld,picec,piced,psnow)
! purpose: a version of icemod::iceini which can be called more than once
!          to restore a model snapshot
!
  use restartmod
  use icemod
  implicit none
  real :: pts(NHOR)
  real :: psst(NHOR)
  real :: pmld(NHOR)
  real :: picec(NHOR)
  real :: piced(NHOR)
  real :: psnow(NHOR)
  ! icemod pre-initialisation
  xflxice2 = 0.
  xqmelt = 0.
!
! from iceini
!
  if (mypid == NROOT) then
     call get_restart_integer('nstep',nstep)
     call get_restart_integer('naccuice',naccuout)
     call get_restart_integer('naccuo',naccuo)
  endif
  
  call mpbci(nstep)
  call mpbci(naccuout)
  call mpbci(naccuo)
  
  call mpgetgp('xls'     ,xls     ,NHOR, 1)
  call mpgetgp('xts'     ,xts     ,NHOR, 1)
  call mpgetgp('xicec'   ,xicec   ,NHOR, 1)
  call mpgetgp('xiced'   ,xiced   ,NHOR, 1)
  call mpgetgp('xsnow'   ,xsnow   ,NHOR, 1)
  call mpgetgp('cheat'   ,cheat   ,NHOR, 1)
  call mpgetgp('cfresh'  ,cfresh  ,NHOR, 1)
  call mpgetgp('ctaux'   ,ctaux   ,NHOR, 1)
  call mpgetgp('ctauy'   ,ctauy   ,NHOR, 1)
  call mpgetgp('cust3'   ,cust3   ,NHOR, 1)
  call mpgetgp('xflxicea',xflxicea,NHOR, 1)
  call mpgetgp('xheata'  ,xheata  ,NHOR, 1)
  call mpgetgp('xofluxa' ,xofluxa ,NHOR, 1)
  call mpgetgp('xqmelta' ,xqmelta ,NHOR, 1)
  call mpgetgp('xcfluxa' ,xcfluxa ,NHOR, 1)
  call mpgetgp('xsmelta' ,xsmelta ,NHOR, 1)
  call mpgetgp('ximelta' ,ximelta ,NHOR, 1)
  call mpgetgp('xtsfluxa',xtsfluxa,NHOR, 1)
  call mpgetgp('xhnewa'  ,xhnewa  ,NHOR, 1)
  
  call mpgetgp('xclicec' ,xclicec ,NHOR,14)
  call mpgetgp('xcliced' ,xcliced ,NHOR,14)
  call mpgetgp('xclsst'  ,xclsst  ,NHOR,14)
  call mpgetgp( 'osst',osst,nhor,1 )
  call mpgetgp( 'osss',osss,nhor,1 )
  if (newsurf .eq. 1) then
     call read_ice_surface
  else
     call mpgetgp( 'xclicec',xclicec,nhor,14 )
     call mpgetgp( 'xcliced',xcliced,nhor,14 )
     call mpgetgp( 'xclsst',xclsst,nhor,14 )
  endif
  xtfreeze(:) = TFREEZE
!
!     read flux correction
!
  if (nfluko == 1) then
     call mpsurfyear( 'xflxice',xflxice,nhor )
  endif
!
!     restart ocean
!
      call read_internal_restart_oceanmod(xsst,xmld,xoflux)
!
!     mit-coupling
!
!
!     initialize mit-ocean coupling
!
      xtfreeze(:) = tfreeze
      if(nmit > 0) then
       if(nrestart /= 0) then
        where(xls(:) < 1.)
           xsst(:)=osst(:)
           xsss(:)=osss(:)
           xtfreeze(:)=-1.*xsss(:)*stotf+TMELT
        endwhere
       endif
      endif
!
!     copy output from icemod
!
      pts(:)=xts(:)
      psst(:)=xsst(:)
      pmld(:)=xmld(:)
      picec(:)=xicec(:)
      piced(:)=xiced(:)
      psnow(:)=xsnow(:)
!
end subroutine read_internal_restart_icemod

subroutine read_internal_restart_oceanmod(psst,pmld,piflux)
! purpose: a version of oceanmod::oceanini which can be called more than once
!          to restore a model snapshot
!
  use restartmod
  use oceanmod
  implicit none
  real :: psst(NHOR),pmld(NHOR),piflux(NHOR)
  !
  ! from oceanini (oceanmod)
  !
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
!
!     read flux correction
!
      if (nfluko == 1) then
         call mpsurfyear('yfsst',yfsst,NHOR)
      endif
!
!     copy output from seamod
!
      psst(:)=ysst(:,1)
      pmld(:)=ymld(:,1)
      piflux(:)=yiflux(:)
  
end subroutine read_internal_restart_oceanmod

subroutine read_internal_restart_seamod
! purpose: a version of seamod::seaini which can be called more than once
!          to restore a model snapshot
!
  use restartmod
  use seamod
  implicit none
  !
  call read_internal_restart_icemod(cts,csst,cmld,cicec,ciced,csnow)
  !
  !     set puma surface variables
  !

  where(dls(:) < 0.5)
!     dts(:)=cts(:)
!     dicec(:)=cicec(:)
!     diced(:)=ciced(:)
     dsnow(:)=csnow(:)
     dsst(:)=csst(:)
     dmld(:)=cmld(:)
!     dt(:,NLEP)=dts(:)
!     dqs(:)=dq(:,NLEP)
  endwhere
  !
  ! from seaini (seamod)
  !
  if (mypid == NROOT) then
     call get_restart_integer('naccua',naccua)
  endif
  call mpbci(naccua)
  
  call mpgetgp('dts'   ,dts   ,NHOR,1) ! surface temperature (K) 
  call mpgetgp('cheata',cheata,NHOR,1) ! total heat flux (w/m2)                        
  call mpgetgp('cpmea' ,cpmea ,NHOR,1) ! P-E (m/s)                                     
  call mpgetgp('cprsa' ,cprsa ,NHOR,1) ! snow precipitation (m/s)                      
  call mpgetgp('croffa',croffa,NHOR,1) ! Runoff (m/s) (puma -> ocean)                  
  call mpgetgp('ctauxa',ctauxa,NHOR,1) ! u-stress (pa)                                 
  call mpgetgp('ctauya',ctauya,NHOR,1) ! v-stress (pa)                                 
  call mpgetgp('cust3a',cust3a,NHOR,1) ! ustar**3 (m3/s3)                              
  call mpgetgp('cshfla',cshfla,NHOR,1) ! surface sens. heat flx (w/m2) (puma->ice)     
  call mpgetgp('cshdta',cshdta,NHOR,1) ! deriv. of dshfl (w/(m2 s))(puma->ice)         
  call mpgetgp('clhfla',clhfla,NHOR,1) ! surface latent heat flx (w/m2) (puma->ice)    
  call mpgetgp('clhdta',clhdta,NHOR,1) ! deriv. of dlhfl (w/(m2 s))(puma->ice)         
  call mpgetgp('cswfla',cswfla,NHOR,1) ! net solar radiation (w/m2) (puma -> ice)      
  call mpgetgp('clwfla',clwfla,NHOR,1) ! surface thermal radiation (w/m2) (puma -> ice)

!  where(dls(:) < 0.5)
!     dqs(:)  = rdbrv*ra1*EXP(ra2*(dt(:,NLEP)-TMELT)   &
!          &          /(dt(:,NLEP)-ra4))/psurf
!     dqs(:)  = dqs(:)/(1.-(1./rdbrv-1.)*dqs(:))
!     dq(:,NLEP) = dqs(:)
!     drhs(:)=drhssea*(1.-dicec(:))+drhsice*dicec(:)
!     dalb(:)=albsea*(1.-dicec(:))   &
!          &         +dicec(:)*AMIN1(albice,0.5+0.025*(273.-dts(:)))
!     dz0(:)=dz0sea*(1.-dicec(:))+dz0ice*dicec(:)
!  endwhere

end subroutine read_internal_restart_seamod
