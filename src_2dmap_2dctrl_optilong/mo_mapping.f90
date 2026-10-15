

!-*-f90-*-

! This file is cpp-preprocessed. Be sure to edit mo_mapping.F90

! The following two macros are mutually exclusive. First one wins, usually.



!#define NX_OCEAN 7920
!#define NX_OCEAN 237600 - last used by Guokun
!#define NX_OCEAN 6035040 - last used by Yulia for optiflux_and_TS_and_14atmpar

! 

MODULE mo_mapping 
! author: FastOpt GmbH
! purpose:
!------------------------------------------------------------------
! map vector of control parameters to physical model variables
!------------------------------------------------------------------
! contains:
!     subroutine cont2mod(n,xarg)
!     subroutine mod2cont(n,par)
!     subroutine addprior(n,x,m,y)
!     subroutine mod2y(m,y)
! public: everything
!
  use substitutes
  use pumamod, only : gt ,st, stp, stm, nroot, mypid,  nlev,nlep, nesp, ntru
  use pumamod, only : sz,szm,szp,sd,sdm,sdp,sq,sqm,sqp,sp,spp,spm,psurf
  use pumamod, only :gaussw,ct,t0,nugp,nlon,nlat,nesp,npackgp,spnorm,nhor,nlpp,nspp
  use pumamod, only :tdissd,tdissz,tdisst,tdissq,tfrc, dt,dls,ww
  use pumamod, only :th2ocg,tswr1g,tswr2g,tswr2g,tswr3g,tpofmtg,acllwrg,albgmaxg,rcritg,alpha1g,beta1g,gamma1g
  use pumamod, only :th2ocm,tswr1m,tswr2m,tswr2m,tswr3m,tpofmtm,acllwrm,albgmaxm,rcritm,alpha1m,beta1m,gamma1m
  use pumamod, only: ss_field_gp,ss_flux_gp, ss_field_hor,ss_flux_hor
  use calmod, only : n_days_per_year
!  use miscmod
  use fluxmod, only : vdiff_lamm,vdiff_d,vdiff_c,vdiff_b
  use radmod, only : tswr1,tswr2,tswr3,rcl1,rcl2,acllwr,acl2,tpofmt,clgray,th2oc
  use rainmod, only : rcrit
  use surfmod , only :
  use landmod, only : albsmin,albsmax,albsminf,albsmaxf,albland,albsmax
  use landmod, only : albgmin,albgmax,dz0land,dalbclim,drhsland,rlue,tau_soil
  use landmod, only : pgrow,plai,pgs,pz0_max
  use landmod, only : nlandt,nlandw,dtcl,dwcl,dts
  use vegmod, only : q10,cveg_a,cveg_b,ct_crit
  use seamod, only: dz0sea, dz0ice, charnock
  use icemod, only : xmind
  use oceanmod, only : vdiffk
  IMPLICIT NONE
  
  integer, parameter :: nx=222659872 !2 !3 ! dimensions of full control vector
  integer, parameter :: nx_ocean=222631200 ! dimensions of ocean control vector
  integer, parameter :: ny=1 !dimension of result vector
  integer :: iounit = 15






  real    :: prior(nx),std_prior(nx)
  integer :: ix ! index of x, marks last mapped control variable
  real    :: xsave(222659872)
CONTAINS

SUBROUTINE cont2mod(n,xarg)
! purpose: cont2mod maps natural params xarg onto model params
! input:
!    integer :: n     ! dimension of control vector
!    real :: xarg(222659872) ! control vector
!  contains
!      subroutine setmod(pmod,kix)
!      subroutine setmodarray(kmod,pmod,kix)
!
  IMPLICIT NONE
  integer :: n,ig
  real :: xarg(222659872)

  real :: zp



  real :: dtmp(NUGP)
  INTEGER k,lock_i,lock_j,lock_k
  
  xsave(:) = xarg(:) ! make a copy to module variable for later use in mo_mapping_ocean
  ix = 0 ! assumes that this is the first mapping routine


  if (mypid==nroot) then
 
! selected from the 36 parameters based on adjoint sensitivity
!     22 16   7  4  31 34  6 5 32 3 33
     call setmodarray_d01(NUGP,th2ocg(:),1,0.3,"th2oc")
     call setmodarray_d01(NUGP,acllwrg(:),NUGP+1,0.25,"acllwr")
     call setmodarray_d01(NUGP,tswr1g(:),2*NUGP+1,0.2,"tswr1")
     call setmodarray_d01(NUGP,tswr2g(:),3*NUGP+1,0.5,"tswr2")
     call setmodarray_d01(NUGP,tswr3g(:),4*NUGP+1,0.2,"tswr3")
     call setmodarray_d01(NUGP,tpofmtg(:),5*NUGP+1,1.1,"tpofmtg")
     call setmodarray_d01(NUGP,albgmaxg(:),6*NUGP+1,1.0,"albgmaxg")
     call setmodarray_d01(NUGP,rcritg(:,7),7*NUGP+1,1.0,"rcrit7")
     call setmodarray_d01(NUGP,rcritg(:,8),8*NUGP+1,1.0,"rcrit8")
     call setmodarray_d01(NUGP,rcritg(:,9),9*NUGP+1,1.0,"rcrit9")
     call setmodarray_d01(NUGP,rcritg(:,10),10*NUGP+1,1.0,"rcrit10")
     call setmodarray_d01(NUGP,alpha1g(:),11*NUGP+1,5.0,"alpha1")
!     call setmodarray_d01(NUGP,beta1g(:),12*NUGP+1,5.0,"beta1")
    !  call setmodarray_d01(NUGP,gamma1g(:),13*NUGP+1,5.0)
!A.B 
    call setmodarray_d01(NUGP,ss_flux_gp(:),12*NUGP+1,0.3,"sf")
    call setmodarray_d01(NUGP,ss_field_gp(:),13*NUGP+1,0.5,"ss")

! smooth paramters to make them smooth
     call smooth_data( th2ocg(:), 1 )
     call smooth_data( acllwrg(:), 1 )
     call smooth_data( tswr1g(:), 1 )
     call smooth_data( tswr2g(:), 1 )
     call smooth_data( tswr3g(:), 1 )
     call smooth_data( tpofmtg(:), 1 )
     call smooth_data( albgmaxg(:), 1 )
     call smooth_data( th2ocg(:), 1 )
     call smooth_data( acllwrg(:), 1 )
     call smooth_data( tswr1g(:), 1 )
     call smooth_data( tswr2g(:), 1 )
     call smooth_data( tswr3g(:), 1 )
     call smooth_data( tpofmtg(:), 1 )
     call smooth_data( albgmaxg(:), 1 )
     dtmp(:)= rcritg(:,7)
     call smooth_data( dtmp(:), 1 )
     rcritg(:,7) =dtmp(:)
     dtmp(:)= rcritg(:,8)
     call smooth_data( dtmp(:), 1 )
     rcritg(:,8) =dtmp(:)
     dtmp(:)= rcritg(:,9)
     call smooth_data( dtmp(:), 1 )
     rcritg(:,9) =dtmp(:)
     dtmp(:)= rcritg(:,10)
     call smooth_data( dtmp(:), 1 )
     rcritg(:,10) =dtmp(:)
     call smooth_data( alpha1g(:), 1 )
     call smooth_data( beta1g(:), 1 )
     call smooth_data( gamma1g(:), 1 )
     call smooth_data( alpha1g(:), 1 )
!     call smooth_data( beta1g(:), 1 )
!     call smooth_data( gamma1g(:), 1 )
! A.B 
     call smooth_data( ss_flux_gp(:), 1 ) 
     call smooth_data( ss_field_gp(:), 1 )
     

 end if

      call mpbci(ix)
      call mpscgp(th2ocg,th2ocm,1)
      call mpscgp(acllwrg,acllwrm,1)
      call mpscgp(tswr1g,tswr1m,1)
      call mpscgp(tswr2g,tswr2m,1)
      call mpscgp(tswr3g,tswr3m,1)
      call mpscgp(tpofmtg,tpofmtm,1)
      call mpscgp(albgmaxg,albgmaxm,1)
      call mpscgp(rcritg,rcritm,NLEV)
      call mpscgp(alpha1g,alpha1m,1)
!      call mpscgp(beta1g,beta1m,1)
!      call mpscgp(gamma1g,gamma1m,1)
! A.B 
 !     print*, "Min and Max of ss_field_hor",MINVAL(ss_field_hor(:)),maxval(ss_field_hor(:))
      call mpscgp(ss_flux_gp,ss_flux_hor,1)
      call mpscgp(ss_field_gp,ss_field_hor,1)
      

  contains

    subroutine setmod(pmod,kix)
! purpose: perturb model parameter
! input:
!      integer :: kix ! index
! output:
!      real :: pmod   ! perturbed model parameter value
!
      real :: pmod
      integer :: kix
      if (ix+1 .ne. kix) then
         print*,ix+1,kix
         stop "wrong index in setmod!"
      end if
      pmod = prior(kix) + std_prior(kix) * xsave(kix)
      print '(a,i5,a,g30.20)','setting parameter ',ix+1,' to ',pmod
      ix = kix
    end subroutine setmod

    subroutine setmod_log(pmod,kix)
      real :: pmod
      integer :: kix
      if (ix+1 .ne. kix) then
         print*,ix+1,kix
         stop "wrong index in setmod!"
      end if
        pmod =exp(prior(kix) + std_prior(kix) * xsave(kix))
      print '(a,i5,a,g30.20)','setting parameter ',ix+1,' to ',pmod
      ix = kix
    end subroutine setmod_log
    subroutine setmodarray_log(kmod,pmod,kix)
! purpose: perturb model parameter array
! input:
!      integer :: kix     ! index
!      integer :: kmod    ! dimension of pmod
! output:
!      real :: pmod(pmod) ! perturbed model parameter value
!
      integer :: kmod,kix
      real :: pmod(kmod)
      if (ix+1 .ne. kix) then
         print*,ix+1,kix
         stop "wrong index in setmodarray!"
      end if
      pmod(:) = exp(prior(kix:kix+kmod-1) + std_prior(kix:kix+kmod-1) * xsave(kix:kix+kmod-1))
      print '(a,i6,a,2g30.20)','setting parameter ',ix+1,' to (max/min)',maxval(pmod(:)),minval(pmod(:))
      ix = kix + kmod - 1
    end subroutine setmodarray_log

    subroutine setmod_d01(pmod,kix,parr)
! purpose: perturb model parameter
! input:
!      integer :: kix ! index
! output:
!      real :: pmod   ! perturbed model parameter value
!
      real :: pmod
      real :: parr  ! real value to control upper limit
      integer :: kix
      if (ix+1 .ne. kix) then
         print*,ix+1,kix
         stop "wrong index in setmod!"
      end if
      pmod = parr/(1+exp(-(prior(kix) + std_prior(kix) * xsave(kix)))) 
      print '(a,i6,a,g30.20)','setting parameter ',ix+1,' to ',pmod
      ix = kix
    end subroutine setmod_d01
    subroutine setmodarray_d01(kmod,pmod,kix,parr,varname)
! purpose: perturb model parameter array
! input:
!      integer :: kix     ! index
!      integer :: kmod    ! dimension of pmod
! output:
!      real :: pmod(pmod) ! perturbed model parameter value
!
      integer :: kmod,kix
      real :: pmod(kmod)
      real :: parr
      character(len=*), intent(in) :: varname
      
    ! Special treatment only for ss
      ! if (trim(varname) == 'ss' .and. &
      !   minval(xsave(kix:kix+kmod-1)) == 0.0 .and. &
      !   maxval(xsave(kix:kix+kmod-1)) == 0.0) then

      ! pmod(:) = 0.0

      ! else

      pmod(:) = parr/(1+exp(-(prior(kix:kix+kmod-1) + std_prior(kix:kix+kmod-1) * xsave(kix:kix+kmod-1)))) 
      print '(a,i6,a,2g30.20)','setting parameter ',ix+1,' to (max/min)',maxval(pmod(:)),minval(pmod(:))
      ix = kix + kmod - 1
!      endif
    end subroutine setmodarray_d01

    subroutine setmodarray(kmod,pmod,kix)
! purpose: perturb model parameter array
! input:
!      integer :: kix     ! index
!      integer :: kmod    ! dimension of pmod
! output:
!      real :: pmod(pmod) ! perturbed model parameter value
!
      integer :: kmod,kix
      real :: pmod(kmod)
      if (ix+1 .ne. kix) then
         print*,ix+1,kix
         stop "wrong index in setmodarray!"
      end if
      pmod(:) = prior(kix:kix+kmod-1) + std_prior(kix:kix+kmod-1) * xsave(kix:kix+kmod-1)
      print '(a,i3,a,2g30.20)','setting parameter ',ix+1,' to (max/min)',maxval(pmod(:)),minval(pmod(:))
      ix = kix + kmod - 1
    end subroutine setmodarray

END SUBROUTINE cont2mod

      subroutine enlarger( adj_var, a, size_z )

!      use pumamod, only : nhor,nlat, nlon
      implicit none
      integer :: size_z
      double precision :: a(nlon+2,nlat+2,size_z)
      double precision :: adj_var(nugp,size_z)
      double precision :: aux(nlon,nlat,size_z)
      integer :: loc_i
      integer :: loc_j
      integer :: loc_k

      loc_k = 0
      do loc_j = 1, nlat
        do loc_i = 1, nlon
          loc_k = loc_k+1
          aux(loc_i,loc_j,:) = adj_var(loc_k,:)
        end do
      end do
      a(:,:,:) = 0.
      a(2:nlon+1,2:nlat+1,:) = aux(1:nlon,1:nlat,:)
      a(1,2:nlat+1,:) = a(nlon+1,2:nlat+1,:)
      a(nlon+2,2:nlat+1,:) = a(2,2:nlat+1,:)
      a(:,1,:) = a(:,2,:)
      a(:,nlat+2,:) = a(:,nlat+1,:)
      end subroutine enlarger

      subroutine smooth_data( adj_var, size_z )

      use pumamod, only : nhor,nlat,nlon

      implicit none

      integer :: size_z
      double precision :: a(nlon+2,nlat+2,size_z)
      double precision :: as(nlon,nlat,size_z)
      double precision :: adj_var(nugp,size_z)


      double precision :: aux(nlon,nlat,size_z)
      integer :: loc_i
      integer :: loc_j
      integer :: loc_k
      call enlarger( adj_var, a, size_z )
      print '(a,2g30.20)','smooth_data: ', maxval(a),minval(a)
     do loc_k= 1, size_z
      do loc_j = 2, nlat+1
        do loc_i = 2, nlon+1
          as(loc_i-1,loc_j-1,loc_k) = sum(a(loc_i-1:loc_i+1,loc_j-1:loc_j+1,loc_k))/9.
        end do
      end do
     end do
      loc_k = 0
      do loc_j = 1, nlat
        do loc_i = 1, nlon
          loc_k = loc_k+1
          adj_var(loc_k,:) = as(loc_i,loc_j,:)
        end do
      end do

      end subroutine smooth_data

SUBROUTINE mod2cont(n,par)
! purpose: mod2cont  maps physical model variables to their natural control parameters 'par'
! input:
!    integer :: n     ! dimension of control vector
! output:
!    real :: xarg(222659872) ! control vector
!  contains
!        subroutine setcontrol(pmod)
!        subroutine setcontrolarray(kmod,pmod)
!
  IMPLICIT NONE
  integer :: n
  REAL, DIMENSION(n), INTENT(out) :: par ! control variables

  ix = 0


  
! selected from the 36 parameters based on adjoint sensitivity
     call setcontrolarray_d01(NUGP,th2ocg,0.3)
     call setcontrolarray_d01(NUGP,acllwrg,0.25)
     call setcontrolarray_d01(NUGP,tswr1g,0.2)
     call setcontrolarray_d01(NUGP,tswr2g,0.5)
     call setcontrolarray_d01(NUGP,tswr3g,0.2)
     call setcontrolarray_d01(NUGP,tpofmtg,1.1)
     call setcontrolarray_d01(NUGP,albgmaxg,1.0)
     call setcontrolarray_d01(NUGP,rcritg(:,7),1.0)
     call setcontrolarray_d01(NUGP,rcritg(:,8),1.0)
     call setcontrolarray_d01(NUGP,rcritg(:,9),1.0)
     call setcontrolarray_d01(NUGP,rcritg(:,10),1.0)
     call setcontrolarray_d01(NUGP,alpha1g,5.0)
!     call setcontrolarray_d01(NUGP,beta1g,5.0)
!     call setcontrolarray_d01(NUGP,gamma1g,5.0)

      call setcontrolarray_d01(NUGP,ss_flux_gp,0.3)
      call setcontrolarray_d01(NUGP,ss_field_gp,0.5)


!  par(:) = -.1 ! de-tune parameters for id-twin experiment
  par(:) = 0.
contains
  subroutine setcontrol_log(pmod)
    real :: pmod
    ! guokun use log transform to keep parameter positive
      prior(ix+1) = log(pmod)
      std_prior(ix+1) = prior(ix+1) ! makes x to be fractions of pmod
 if (mypid==nroot) then
      print '(a,i3,a,g30.20)','setting control ',ix+1,' to ',prior(ix+1)
 endif
    ix = ix + 1
  end subroutine setcontrol_log
  subroutine setcontrolarray_log(kmod,pmod)
! purpose: set first guess of control parameter vector section and its error
! input:
!    integer :: kmod    ! dimension of pmod
!    real :: pmod(kmod) ! model quantity
!
    integer :: kmod
    real :: pmod(kmod)
    prior(ix+1:ix+kmod) = log(pmod(:))
    std_prior(ix+1:ix+kmod) = prior(ix+1:ix+kmod) ! makes x to be fractions of pmod
    ix = ix + kmod
  end subroutine setcontrolarray_log
  subroutine setcontrol(pmod)
! purpose: set first guess of control parameter and its error
! input:
!    real :: pmod ! model quantity
!
    real :: pmod
    prior(ix+1) = pmod
    std_prior(ix+1) = prior(ix+1) ! makes x to be fractions of pmod
    ix = ix + 1
  end subroutine setcontrol
  subroutine setcontrol_d10(pmod)
! purpose: set first guess of control parameter and its error
! input:
!    real :: pmod ! model quantity
!
    real :: pmod
    prior(ix+1) = pmod
    std_prior(ix+1) = prior(ix+1)/10. ! makes x to be fractions of pmod
    ix = ix + 1
  end subroutine setcontrol_d10

  subroutine setcontrol_d01(pmod)
! purpose: set first guess of control parameter and its error
! input:
!    real :: pmod ! model quantity
!
    real :: pmod
    prior(ix+1) = -log(1./pmod-1.)
    std_prior(ix+1) = prior(ix+1) ! makes x to be fractions of pmod
    ix = ix + 1
  end subroutine setcontrol_d01
  subroutine setcontrolarray_d01(kmod,pmod,parr)
! purpose: set first guess of control parameter vector section and its error
! input:
!    integer :: kmod    ! dimension of pmod
!    real :: pmod(kmod) ! model quantity
!
    integer :: kmod
    real :: pmod(kmod)
    real :: parr  ! real value (0 1] to tune upper range
    prior(ix+1:ix+kmod) = -log(parr/pmod(:)-1.)
    std_prior(ix+1:ix+kmod) = prior(ix+1:ix+kmod) ! makes x to be fractions of pmod
    ix = ix + kmod
  end subroutine setcontrolarray_d01
  subroutine setcontrolarray(kmod,pmod)
! purpose: set first guess of control parameter vector section and its error
! input:
!    integer :: kmod    ! dimension of pmod
!    real :: pmod(kmod) ! model quantity
!
    integer :: kmod
    real :: pmod(kmod)
    prior(ix+1:ix+kmod) = pmod(:)
    std_prior(ix+1:ix+kmod) = prior(ix+1:ix+kmod) ! makes x to be fractions of pmod
    ix = ix + kmod
  end subroutine setcontrolarray

END SUBROUTINE mod2cont



subroutine addprior(n,x,m,y)
! purpose: add background term to y
! input:
!    integer :: n ! dimension of x
!    integer :: m ! dimension of y
!    real :: x(n) ! control vector
!    real :: y(m) ! cost function
! output:
!    real :: y(m) ! cost function
!
  integer :: n,m,i,j,k
  real :: x(n)
  real :: y(m)
  real :: bk_cost
  real :: par_mtx(NLAT,NLON)
  real :: par_mtxe(NLAT+4,NLON+4)
  real :: par_mtxs(NLAT,NLON)
  real :: flt(5,5)

    if (mypid==nroot) then

      call Gaussian_filter(flt,1.)
! enlarge matrix
    call enlarger_data( th2ocg(:)-th2oc, par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/th2oc)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo   
       print*,'bkcost is:' ,sum(par_mtxe),sum(th2ocg),sum(th2ocg(:)-th2oc)
    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( acllwrg(:)-acllwr, par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/acllwr)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( tswr1g(:)-tswr1, par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/tswr1)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 
    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( tswr2g(:)-tswr2, par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/tswr2)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( tswr3g(:)-tswr3, par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/tswr3)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( tpofmtg(:)-tpofmt, par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/tpofmt)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( albgmaxg(:)-albgmax, par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/albgmax)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( rcritg(:,7)-rcrit(7), par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/rcrit(7))**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( rcritg(:,8)-rcrit(8), par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/rcrit(8))**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( rcritg(:,9)-rcrit(9), par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/rcrit(9))**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( rcritg(:,10)-rcrit(10), par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/rcrit(10))**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( alpha1g(:)-1., par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/1.)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( beta1g(:)-1., par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/1.)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

    par_mtxe(:,:)=0.;
    par_mtxs(:,:)=0.;
! enlarge matrix
    call enlarger_data( gamma1g(:)-1., par_mtxe )
      do i=3,NLAT+2
         do j=3,NLON+2
           par_mtxs(j-2,i-2)=(sum(par_mtxe(i-2:i+2,j-2:j+2)*flt)/1.)**2;
         end do 
       enddo
        bk_cost=0.
       do i=1,NLAT
        do j=1,NLON
          bk_cost = bk_cost+ 0.5*par_mtxs(j,i)
        enddo
       enddo 

     print '(a,g20.12,g20.12)','adding prior to total : ',bk_cost*100, y(1)
 
  y(1) = y(1) + bk_cost*100.  ! weight bkcost by 2 times 
  endif ! if mypid==nroot
end subroutine addprior

subroutine Gaussian_filter(flt,radius)
   implicit none
   real :: flt(1:5,1:5)
   real:: radius,sumflt
   integer :: i,j
    do i=1,5
     do j=1,5
       flt(i,j)=(1./2./3.14/radius**2)*exp(-(i-3.)**2/2./radius**2-(j-3.)**2/2./radius**2);
     enddo
    enddo
  
    sumflt=sum(flt)
     do i=1,5
     do j=1,5
       flt(i,j)=flt(i,j)/sumflt;
     enddo
    enddo
  
end subroutine Gaussian_filter

      subroutine enlarger_data( adj_var, par_mtxe )
      use pumamod, only : nugp
      implicit none
      double precision :: par_mtxe(68,36)
      double precision :: adj_var(nhor)


      double precision :: aux(64,32)
      integer :: loc_i
      integer :: loc_j
      integer :: loc_k

      loc_k = 0
      do loc_j = 1, 32
        do loc_i = 1, 64
          loc_k = loc_k+1
          aux(loc_i,loc_j) = adj_var(loc_k)
        end do
      end do
      
      par_mtxe(:,:) = 0.
      par_mtxe(3:66,3:34) = aux(1:64,1:32)
      par_mtxe(1:2,3:34) = aux(63:64,1:32)
      par_mtxe(67:68,3:34) = aux(1:2,1:32)
      par_mtxe(:,1:2) = 0.
      par_mtxe(:,35:36) = 0.

      end subroutine enlarger_data

SUBROUTINE mod2y(m,y)
! purpose: map model result to y
! input:
!    integer :: m ! dimension of y
! output:
!    real :: y(m) ! cost function
!
  use observation, only : cost,yname
  IMPLICIT NONE
  real :: y(m)
  integer :: m


  y(1) = cost

END SUBROUTINE mod2y

END MODULE mo_mapping

SUBROUTINE init_global_mapping
! author: FastOpt GmbH
! purpose: initialise mapping of model and control variables
!
  USE mo_mapping
  use legmod
  use fftmod
  implicit none
  integer :: get_free_unit
! defaults
END SUBROUTINE init_global_mapping

integer function get_free_unit(kunit)
! author: FastOpt
! purpose: return an unused unit number
! input:
!     integer :: kunit ! unit number before trying higher ones
!
  implicit none
  integer :: kunit
  integer :: iunit
  logical :: l_opened
  l_opened=.true.
  iunit=max(6,kunit-1)
  do while (l_opened) ! find a free unit number
     iunit=iunit+1
     inquire(UNIT=iunit,OPENED=l_opened)
  end do
  get_free_unit=iunit
end function get_free_unit
