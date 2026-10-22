!-*-f90-*-
!======================================================================
!   file:               prgdfpmin.F90
!   purpose:            driver for optimisation using dpfmin
!   creation date:      09 05
!
!   prepared by:
!   FastOpt GmbH
!   Hamburg, Germany
!
!   Contact: http://www.FastOpt.com
!======================================================================

      module mo_n
      implicit none
      integer, parameter :: n = 263288032  ! set by preprocessor
      integer, parameter :: m = 1  ! set by preprocessor           
      end module mo_n

      module mo_h
      implicit none
      integer, parameter :: mh = 1
      integer, parameter :: mx = 263288032 ! set by preprocessor

      integer ih
      real(kind=8) hx(mx,mh)
      real(kind=8) hf(mh)
      real(kind=8) hg(mx,mh)
      real(kind=8) :: lx(mx) = 0.
      real(kind=8) :: mask(mx)
     end module mo_h

      program main
!*************************************************************
      use mo_n
      implicit none

!==============================================================
! get the number of control variables
!==============================================================
!      call setfunc(n,m)

!-----------------------------------------
! call the subroutine
!-----------------------------------------
      call opti

      end

      subroutine opti
      use mo_n
      implicit none

      REAL(KIND=8)    objf(1)

      REAL(KIND=8)   x(n)
      REAL(KIND=8)   adx(n)
      integer ifail

!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================
!      call initfunc( n, x )

      ifail=0
      OPEN (39,file='x.b',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'DPFMIN: parameter values set by initmod'
      else
         read(39) x
         print*,'DPFMIN: parameter values read from x.b'
         print*,x(1:8)
      endif
      CLOSE(39)
! reinitialized 
     
!       x(8)=x(8)+0.02
!       x(11)=x(11)-0.02
!==============================================================
! optimize the costfunction using control variables
!==============================================================
      call optimum( n, x, objf, adx )

!==============================================================
! start of postprocessor
!==============================================================
!      call postfunc( n, x, 1, objf )

! save optimum
      OPEN (39,file='xf.b',status='unknown',form='unformatted')
      WRITE(39) x
      CLOSE(39)

      end

      subroutine optimum( n, x, fc, adx )
      use mo_h
      implicit none

!==============================================================
! declare call parameters
!==============================================================
      integer n
      integer :: M = 4, NWORK
      REAL(KIND=8)    x(n)
      REAL(KIND=8)    fc(1)
      REAL(KIND=8)    adx(n)

!==============================================================
! control vector layout: EDIT HERE when grid or run length changes
!
!   x = [ atm : atm_fld x atm_gp ]              14 x 64*32
!       [ ocn : temp, salt, kapgm, kapredi ]     4 x ocn_gp
!       [ flx : hflux, sflux, fu, fv ]           4 x flx_gp
!==============================================================
      integer, parameter :: atm_fld = 14             ! number of atm parameters
      integer, parameter :: atm_gp  = 32*64          ! atm grid points
      integer, parameter :: ocn_gp  = 90*52*23       ! ocean points per parameter
      integer, parameter :: flx_gp  = 90*52*360*39   ! daily flux points per flux (360 d x 39 yr)
      integer, parameter :: i0_ocn  = atm_fld*atm_gp       ! last index before ocean block
      integer, parameter :: i0_flx  = i0_ocn + 4*ocn_gp    ! last index before flux block
      integer :: nopt     ! packed elements given to LBFGS (atm + wet ocean)

!==============================================================
! declare local variables
!==============================================================
      integer    ioptns, iter, i, unith, ifail,j
      parameter( ioptns = 7 )
      parameter( unith = 8 )
      REAL(KIND=8) epsg
      real,external ::  f
      namelist/ OPTIM/ epsg,iter
      integer :: IFLAG, ICALL, IPrint(2)
      real, allocatable :: W(:)
      real   ::fo,fn
      real, allocatable    :: DIAG(:)
      real    :: fc0, EPS, XTOL
      logical :: DiagCO
      real, external :: LB2
      logical :: lexist1,lexist2,lexist3
      integer ::  ncount, ncount1,nrec,ncountp
      REAL(KIND=8), allocatable ::   xp(:)
      REAL(KIND=8), allocatable ::   adxp(:)
! precondition factors
      real(kind=8):: wfac(1:22)
      REAL(KIND=8) :: mlim=30., abs_gm
      character(len=3) ::cname
      character(len=21) :: cmd
      character(len=51) :: cmd1
      INTEGER*4 system,status
      integer iosn

! check that n matches the layout above
      if (n .ne. i0_flx + 4*flx_gp) then
         print*,'ERROR: n =',n,' but layout gives ',i0_flx + 4*flx_gp
         stop
      endif

! the first test
!      wfac(1:5)= (/1.,1.,1.,1.,0.5/)
!      wfac(6:10)=(/5.,5.,5.,5.,0.5/)
!      wfac(11:14)=(/5.,1.,1.,2./)
!      wfac(15:16)= 50.
!      wfac(17:18) = (/4.,2./)

! the second test
!      wfac(1:5)= (/1.0,1.0,1.,1.,1./)
!      wfac(6:10)=(/5.0,1.0,1.0,1.0,0.2/)
!      wfac(11:14)=(/1.,1.,.3,.2/)
!      wfac(15:16)= (/1.,1./)
!      wfac(17:18) = (/.001,.001/)

! original from Guokun
!      wfac(1:5)= (/1.0,1.0,1.,1.,1./)
!      wfac(6:10)=(/1.0,1.0,1.0,1.0,1./)
!      wfac(11:14)=(/1.,1.,1.,1./)
!      wfac(17:18) = (/.1,.1/)
! Aug 2023
!     wfac(15:16)= (/1000,1000/) leads to too large salt incr
!     wfac(15:16)= (/100,100/) also leads to large salt incr
!     wfac(15:16)= (/10,10/) stops before 1st iter
      wfac(1:5)  =(/.1, 1., .1, .1, .1/)
      wfac(6:10) =(/1., .3, .3, .3, .1/)
      wfac(11:14)=(/.3, 1., .3, .3/)
! temp, salt, kapgm, kapredi
      wfac(15:18) = (/.2, .2, 1., 1./)
! hflux,sflux fu fv
      wfac(19:22)=(/100, 30, 100, 200/)
      
      epsg  = 1.E-2
      iter = -1
 
!======================================================
! read new optimisation switches from file
!======================================================
      open( ioptns, FILE   = 'dfpmin.par'&
                 , ACCESS = 'SEQUENTIAL'&
                 , FORM   = 'FORMATTED'  )
      read( ioptns, OPTIM, iostat = ifail )
      if (ifail==0) then
         print *, ' OPTIMI : options have been read'
      else if (ifail .lt. 0 ) then
         print *, ' OPTIMI : dfpmin.par empty?'
      else
         print *, ' OPTIMI : error reading from dfpmin.par.'
      end if
      close( ioptns )
      print OPTIM

!==============================================================
! perform optimisation
!==============================================================
      ih = 0 ! our own count of iterations
      hg = 0. 
      hf = 0.
      hx = 0.
! infomation needed for optimization
      IPrint(1) = 1
      IPrint(2) = 0

      DiagCO = .false.
      XTOL = 1.0D-20
      ICALL = 0
      IFLAG = 0     

! read in mask file
      inquire(file='mask.b',exist=lexist1)
      if(lexist1)then
      OPEN (39,file='mask.b',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'mask.b infor is wrong'
         stop
      else
         read(39) mask
      endif
      close(39)
      endif

! ncountp: all active elements; nopt: atm + wet ocean (what LBFGS optimises)
       ncountp=0
       nopt=0
       do i=1,n
         if(abs(mask(i)-1.) .lt. 0.01)then
         ncountp=ncountp+1  
         if(i .le. i0_flx) nopt=nopt+1
        endif
       enddo
       print*,'all net grid is: ', ncountp,n,maxval(mask)
       print*,'optimised elements (atm + ocean): ', nopt

 

      lexist1 = .false.
      lexist2 = .false.
      lexist3 = .false.

      inquire(file='fc.b_000',exist=lexist1)
      inquire(file='x.b_001',exist=lexist2)
      inquire(file='g.b_001',exist=lexist3)


      if(lexist1 .and. lexist2 .and. lexist3) then
      ifail=0
             ih = 1
      OPEN (39,file='x.b_001',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'x.b infor is wrong'
         stop
      else
         read(39) x
         print*,'x initial: ',x(1:8)
      endif
      CLOSE(39)
      ifail=0
      OPEN (39,file='g.b_001',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'g.b infor is wrong'
         stop
      else
         read(39) adx
         print*,'g initial: ',adx(1:8)
      endif
      CLOSE(39)

      ifail=0
      OPEN (39,file='fc.b_000',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'fc.b infor is wrong'
         stop
      else
         read(39) fc
         print*,'fc initial: ',fc(1)
      endif
      CLOSE(39)

      else
       ih = 0
!    calculate cost function F(x)
	fc(1)=f(x)

!    calculate gradient of F(x)  
 
      call df( x, adx )
      endif

      NWORK = ncountp*(2*M+1) + 2*M
      allocate(W(NWORK))    
      allocate(xp(ncountp))  
      allocate(adxp(ncountp))  
      allocate(DIAG(ncountp)) 
!      hx(1:8,1)=x(1:8)
!      hg(1:8,1)=adx(1:8)
!      hf(1)=fc(1)
!      print'(a,i5,a,1e26.16)','i = ', ih,' fc = ',fc(1)
!      print'(a,i5,a,8e26.16)','i = ', ih,' x = ',hx(1:8,1)
!      print'(a,i5,a,8e26.16)','i = ', ih,' g = ',hg(1:8,1)

       do ICALL = 1, iter
         fo=fc(1)
! 1)  filtering the sensivity
! a)  abs>global*50times
         do i=1,atm_fld
           abs_gm=sum(abs(adx((i-1)*atm_gp+1:i*atm_gp)))/real(atm_gp)
           do j=1,atm_gp
             if(abs(adx((i-1)*atm_gp+j))>=abs_gm*mlim)then
               adx((i-1)*atm_gp+j)=mlim*abs_gm*sign(1.0,adx((i-1)*atm_gp+j))
!               print*,'filter atm',i,j
              endif
            enddo
!           print*,'filter atm',i
         enddo

! ocean temp
        
        abs_gm=0.; 
        nrec=0; 
       do i=i0_ocn+1,i0_ocn+ocn_gp
         if(abs(mask(i)-1.) .lt. 0.01)then
         nrec=nrec+1  
         abs_gm=abs_gm+abs(adx(i))
        endif
        enddo
        abs_gm=abs_gm/real(nrec) 
         do i=i0_ocn+1,i0_ocn+ocn_gp  
           if(abs(mask(i)-1.) .lt. 0.01 .and. abs(adx(i)).ge. mlim*abs_gm)then
             adx(i)= mlim*abs_gm*sign(1.0,adx(i)) 
!               print*,'filter temi',i  
            endif  
         enddo
! ocean salt

        abs_gm=0.; 
        nrec=0; 
       do i=i0_ocn+1+ocn_gp,i0_ocn+ocn_gp*2
         if(abs(mask(i)-1.) .lt. 0.01)then
         nrec=nrec+1  
         abs_gm=abs_gm+abs(adx(i))
        endif
        enddo
        abs_gm=abs_gm/real(nrec) 
         do i=i0_ocn+1+ocn_gp,i0_ocn+ocn_gp*2
           if(abs(mask(i)-1.) .lt. 0.01 .and. abs(adx(i)).ge. mlim*abs_gm)then
             adx(i)= mlim*abs_gm*sign(1.0,adx(i)) 
!               print*,'filter temi',i  
            endif  
         enddo
!kapgm

        abs_gm=0.; 
        nrec=0; 
       do i=i0_ocn+1+ocn_gp*2,i0_ocn+ocn_gp*3
         if(abs(mask(i)-1.) .lt. 0.01)then
         nrec=nrec+1  
         abs_gm=abs_gm+abs(adx(i))
        endif
        enddo
        abs_gm=abs_gm/real(nrec) 
         print*,'abs_gm', abs_gm,nrec           
         do i=i0_ocn+1+ocn_gp*2,i0_ocn+ocn_gp*3
           if(abs(mask(i)-1.) .lt. 0.01 .and. abs(adx(i)).ge. mlim*abs_gm)then
             adx(i)= mlim*abs_gm*sign(1.0,adx(i)) 
!               print*,'filter temi',i  
            endif  
         enddo
        
!kapredi

        abs_gm=0.; 
        nrec=0; 
       do i=i0_ocn+1+ocn_gp*3,i0_ocn+ocn_gp*4
         if(abs(mask(i)-1.) .lt. 0.01)then
         nrec=nrec+1  
         abs_gm=abs_gm+abs(adx(i))
        endif
        enddo
        abs_gm=abs_gm/real(nrec) 

         do i=i0_ocn+1+ocn_gp*3,i0_ocn+ocn_gp*4
           if(abs(mask(i)-1.) .lt. 0.01 .and. abs(adx(i)).ge. mlim*abs_gm)then
             adx(i)= mlim*abs_gm*sign(1.0,adx(i)) 
!               print*,'filter temi',i  
            endif  
         enddo

!hflux

        abs_gm=0.; 
        nrec=0; 
        do i=i0_flx+1,i0_flx+flx_gp  
           if(abs(mask(i)-1.) .lt. 0.01)then
              nrec=nrec+1  
              abs_gm=abs_gm+abs(adx(i))
           endif
        enddo
        abs_gm=abs_gm/real(nrec) 

        do i=i0_flx+1,i0_flx+flx_gp 
           if(abs(mask(i)-1.) .lt. 0.01 .and. abs(adx(i)).ge. mlim*abs_gm)then
             adx(i)= mlim*abs_gm*sign(1.0,adx(i)) 
!               print*,'filter temi',i  
            endif  
         enddo

!sflux

        abs_gm=0.; 
        nrec=0; 
        do i=i0_flx+1+flx_gp,i0_flx+flx_gp*2
           if(abs(mask(i)-1.) .lt. 0.01)then
              nrec=nrec+1  
              abs_gm=abs_gm+abs(adx(i))
           endif
        enddo
        abs_gm=abs_gm/real(nrec) 

         do i=i0_flx+1+flx_gp,i0_flx+flx_gp*2    
           if(abs(mask(i)-1.) .lt. 0.01 .and. abs(adx(i)).ge. mlim*abs_gm)then
             adx(i)= mlim*abs_gm*sign(1.0,adx(i)) 
!               print*,'filter temi',i  
            endif  
         enddo
!fu

        abs_gm=0.; 
        nrec=0; 
         do i=i0_flx+1+flx_gp*2,i0_flx+flx_gp*3  
           if(abs(mask(i)-1.) .lt. 0.01)then
              nrec=nrec+1  
              abs_gm=abs_gm+abs(adx(i))
           endif
        enddo
        abs_gm=abs_gm/real(nrec) 

        do i=i0_flx+1+flx_gp*2,i0_flx+flx_gp*3
           if(abs(mask(i)-1.) .lt. 0.01 .and. abs(adx(i)).ge. mlim*abs_gm)then
             adx(i)= mlim*abs_gm*sign(1.0,adx(i)) 
!               print*,'filter temi',i  
            endif  
         enddo

!fv

        abs_gm=0.; 
        nrec=0; 
         do i=i0_flx+1+flx_gp*3,i0_flx+flx_gp*4
           if(abs(mask(i)-1.) .lt. 0.01)then
              nrec=nrec+1  
              abs_gm=abs_gm+abs(adx(i))
           endif
        enddo
        abs_gm=abs_gm/real(nrec) 

         do i=i0_flx+1+flx_gp*3,i0_flx+flx_gp*4
           if(abs(mask(i)-1.) .lt. 0.01 .and. abs(adx(i)).ge. mlim*abs_gm)then
             adx(i)= mlim*abs_gm*sign(1.0,adx(i)) 
!               print*,'filter temi',i  
            endif  
         enddo
         
! write out filterd sensivity
         OPEN (39,file='g.b_flt',status='unknown',form='unformatted')
         WRITE(39) adx(:)
         CLOSE(39)
         OPEN (39,file='maskn.b',status='unknown',form='unformatted')
         WRITE(39) mask(:)
         CLOSE(39)

         write(cname,'(I3.3)')ih     
         cmd='mv g.b_flt g.b_flt'//cname
         status = system(trim(cmd))
!         stop
! rescale x and adx for better preconditon
! also remove non-zero data.
! atmosphere part
! wfactor=[1 1 1 1 0.5 5 5 5 5 0.5 5 1 1 2];
         do i=1,atm_fld
           xp((i-1)*atm_gp+1:i*atm_gp)=x((i-1)*atm_gp+1:i*atm_gp)/wfac(i)
         end do
         ncount=i0_ocn
! ocean temp
         do i=i0_ocn+1,i0_ocn+ocn_gp        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               xp(ncount)=x(i)/wfac(15)
              endif
         enddo
! ocean salt
         do i=i0_ocn+1+ocn_gp,i0_ocn+ocn_gp*2        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               xp(ncount)=x(i)/wfac(16)
              endif
         enddo
!kapgm
         do i=i0_ocn+1+ocn_gp*2,i0_ocn+ocn_gp*3        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               xp(ncount)=x(i)/wfac(17)
              endif
         enddo
!kapredi
         do i=i0_ocn+1+ocn_gp*3,i0_ocn+ocn_gp*4        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               xp(ncount)=x(i)/wfac(18)
              endif
         enddo

!hflux
         do i=i0_flx+1,i0_flx+flx_gp       
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               xp(ncount)=x(i)/wfac(19)
              endif
         enddo
!sflux
          do i=i0_flx+1+flx_gp,i0_flx+flx_gp*2      
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               xp(ncount)=x(i)/wfac(20)
              endif
         enddo

!fu
          do i=i0_flx+1+flx_gp*2,i0_flx+flx_gp*3    
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               xp(ncount)=x(i)/wfac(21)
              endif
         enddo
!fv
          do i=i0_flx+1+flx_gp*3,i0_flx+flx_gp*4      
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               xp(ncount)=x(i)/wfac(22)
              endif
         enddo
! rescale sensivity
! atmosphere part

         do i=1,atm_fld
           adxp((i-1)*atm_gp+1:i*atm_gp)=adx((i-1)*atm_gp+1:i*atm_gp)*wfac(i)
         enddo

         ncount=i0_ocn


! ocean temp
         do i=i0_ocn+1,i0_ocn+ocn_gp        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               adxp(ncount)=adx(i)*wfac(15)
              endif
         enddo
! ocean salt
         do i=i0_ocn+1+ocn_gp,i0_ocn+ocn_gp*2        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               adxp(ncount)=adx(i)*wfac(16)
              endif
         enddo
!kapgm
         do i=i0_ocn+1+ocn_gp*2,i0_ocn+ocn_gp*3        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               adxp(ncount)=adx(i)*wfac(17)
              endif
         enddo
!kapredi
         do i=i0_ocn+1+ocn_gp*3,i0_ocn+ocn_gp*4        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               adxp(ncount)=adx(i)*wfac(18)
              endif
         enddo
!hflux
          do i=i0_flx+1,i0_flx+flx_gp        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               adxp(ncount)=adx(i)*wfac(19)
              endif
         enddo
!sflux
          do i=i0_flx+1+flx_gp,i0_flx+flx_gp*2   
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               adxp(ncount)=adx(i)*wfac(20)
              endif
         enddo
!fu
          do i=i0_flx+1+flx_gp*2,i0_flx+flx_gp*3       
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               adxp(ncount)=adx(i)*wfac(21)
              endif
         enddo
!fv
          do i=i0_flx+1+flx_gp*3,i0_flx+flx_gp*4       
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               adxp(ncount)=adx(i)*wfac(22)
              endif
           enddo
!A.B setting the last element to zero    
!         adxp(11*32*64+1:14*32*64)=0
!           adxp(1:12*32*64)=0

! Ocean parameters
!   adxp(i0_ocn:size(adxp))=0
! A.B Making the flux zero
   adxp(nopt+1:size(adxp))=0
	if ( ICALL==1 ) fc0 = fc(1)
	 print*,'OPTIM INFO: ' ,ICALL, fc(1), IFLAG, ih
	 print*,'SIZE INFO: ' ,size(xp),size(adxp),size(DIAG),ncountp
! limit the optimization to the atmospheric and ocean variables
          ncountp=nopt
         call LBFGS(ncountp,4,xp,fc(1),adxp,DIAGCO,DIAG,IPRINT,EPSG,xtol,W,IFLAG)
!  call LBFGS(ncountp,5,xp,fc(1),adxp,DIAGCO,DIAG,IPRINT,EPSG,xtol,W,IFLAG) 
 ! write out x(:) infor for next iteration
! put xp back to x

         print*,'size x1 is ', size(x)

        do i=1,atm_fld
          x((i-1)*atm_gp+1:i*atm_gp)=xp((i-1)*atm_gp+1:i*atm_gp)*wfac(i)
        enddo
         ncount=i0_ocn

! ocean temp
         do i=i0_ocn+1,i0_ocn+ocn_gp        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               x(i)=xp(ncount)*wfac(15)
              else
               x(i)=0.
              endif
         enddo
! ocean salt
         do i=i0_ocn+1+ocn_gp,i0_ocn+ocn_gp*2        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               x(i)=xp(ncount)*wfac(16)
              else
               x(i)=0.
              endif
         enddo
!kapgm
         do i=i0_ocn+1+ocn_gp*2,i0_ocn+ocn_gp*3        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               x(i)=xp(ncount)*wfac(17)
              else
               x(i)=0.
              endif
         enddo
!kapredi
         do i=i0_ocn+1+ocn_gp*3,i0_ocn+ocn_gp*4        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               x(i)=xp(ncount)*wfac(18)
              else
               x(i)=0.
              endif
         enddo
!hflux
          do i=i0_flx+1,i0_flx+flx_gp        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               x(i)=xp(ncount)*wfac(19)
              else
               x(i)=0.
              endif
         enddo
!sflux
          do i=i0_flx+1+flx_gp,i0_flx+flx_gp*2        
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               x(i)=xp(ncount)*wfac(20)
              else
               x(i)=0.
              endif
         enddo
 !fu
          do i=i0_flx+1+flx_gp*2,i0_flx+flx_gp*3     
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               x(i)=xp(ncount)*wfac(21)
              else
               x(i)=0.
              endif
         enddo
!fv
          do i=i0_flx+1+flx_gp*3,i0_flx+flx_gp*4     
             if(abs(mask(i)-1.) .lt. 0.01) then
               ncount=ncount+1
               x(i)=xp(ncount)*wfac(22)
              else
               x(i)=0.
              endif
         enddo
           
         print*,'size x2 is ', size(x),maxval(x),minval(x)

         OPEN (39,file='x.b',status='unknown',form='unformatted')
         WRITE(39) x(:)
         CLOSE(39)
          print*,'x write to "x.b" in opti: ',x(1:6)

         if ( IFLAG<=0 ) exit
          fn=f(x)
          if(fn>=fo)then
            print*,'Linesearchr ',fn
          elseif(fn<fo)then
           fc(1)=fn
           call df( x, adx )
! save state of iteration
      hx(1:8,1)=x(1:8)
      hg(1:8,1)=adx(1:8)
      hf(1)=fc(1)
      print'(a,i5,a,1e26.16)','i = ', ih,' fc = ',fc(1)
      print'(a,i5,a,8e26.16)','i = ', ih,' x = ',hx(1:8,1)
      print'(a,i5,a,8e26.16)','i = ', ih,' g = ',hg(1:8,1)

          endif



      enddo   
      DEALLOCATE(W)
      DEALLOCATE(adxp)
      DEALLOCATE(xp)
      DEALLOCATE(DIAG)
!==============================================================
! evaluate final gradient
!==============================================================

!      call df( x, adx )
      open (unith,file='f.dat',form='formatted')
      do i=1,1
         write(unith,'(i5,6e26.16)')i,hf(i)!,sum(hg(:,i))
      enddo
      close (unith)
      open (unith,file='x.dat',form='formatted')
      do i=1,1
         write(unith,'(i5,13e26.16)')i,sqrt(sum(hx(:,i)**2)),hx(1:min(12,n),i)
      enddo
      close (unith)
      open (unith,file='g.dat',form='formatted')
      do i=1,1
         write(unith,'(i5,13e26.16)')i,sqrt(sum(hg(:,i)**2)),hg(1:min(12,n),i)
      enddo
      close (unith)
      end

!===============================================
! this subroutine evaluates the function
!===============================================
      REAL(KIND=8) function f( x )
      use mo_n
      use mo_h
      implicit none
      REAL(KIND=8)    x(n)
      integer i,ifail
      REAL(KIND=8)    y(1)
      character(len=3) ::cname
      character(len=20) :: cmd
      character(len=50) :: cmd1
      INTEGER*4 system,status
!===============================================
! call adjoint model
!===============================================
      y = 0.
    
       status = system('mpirun -np 4 ./function>function.out')
       if ( status .ne. 0 ) stop 'system: error'
       status = system('rm -r mnc_test*')
!      call func( n, x, 1, y )
!      f = y(1)
      ifail=0
      OPEN (39,file='fc.b',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'cost value is not prepared'
      else
         read(39) y(1)
         print*,'cost is read from fc.b', y(1)
         f = y(1)
      endif
      CLOSE(39)

!      print'(a,g30.20)','fc in f:',f
      ! ATTENTION
      ! if model returns fc=NaN (e.g. if model is called with
      ! irregular parameters x) we want to return a value that
      ! is a bit less than "huge" to make numerical recipes routine
      ! lnsrch-routine continue
      !
      ! condition ".not.fc.lt.huge(fc)" becomes true for fc=NaN
      !
      if(.not.f.le.sqrt(huge(f))) then
         f = sqrt(huge(f))
         print*, 'PRGDFPMIN: avoiding function overflow'
         call flush(6)
      endif
      if(.not.f.ge.-sqrt(huge(f))) then
         f = sqrt(huge(f))
         print*, 'PRGDFPMIN: avoiding function overflow'
         call flush(6)
      endif
 !     print*, 'XXX Save previous point from line search '
 !     OPEN (39,file='lx.b',status='unknown',form='unformatted')
 !     WRITE(39) lx(1:n)
 !     CLOSE(39)
 !     lx(1:n) = x

!         OPEN (39,file='fc.b',status='unknown',form='unformatted')
!         WRITE(39) f
!         CLOSE(39)
         write(cname,'(I3.3)')ih
         cmd='mv fc.b fc.b_'//cname
         status = system(trim(cmd))
      end function f

!===============================================
! this subroutine evaluates the gradient only
! use taf -pure -reverse
! to generate subroutine admodel
!===============================================
      subroutine df( x, adx )
      use mo_n
      use mo_h
      implicit none
      REAL(KIND=8)    x(n)
      REAL(KIND=8)    y(1),  adx(n)

      REAL(KIND=8)    adfc(1)
      character(len=3) ::cname
      character(len=20) :: cmd
      character(len=50) :: cmd1
      INTEGER*4 system,status



      integer i,idx(1:1),ifail

!===============================================
! call adjoint model
!===============================================
      y = 0.

      adfc = 1.
      adx(:) = 0.
!      call func_ad( n, x, adx, 1, y(1), adfc )
       status = system('mpirun -np 4 ./tstadm>tstadm.out')
       if ( status .ne. 0 ) stop 'system: error'
       status = system('rm -r mnc_test*')
      ifail=0
      OPEN (39,file='g.b',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'adx is not prepared in g.b'
      else
         read(39) adx(:)
         print*,'adx is read from g.b'
         print*,adx(1:8)
      endif
      CLOSE(39)
! normalize parameters 
!       adx(1)=adx(1)/10.;adx(2)=adx(2)/200.;adx(3)=adx(3)/50.;adx(4)=adx(4)/100.;
!       adx(5)=adx(5)/100.;adx(6)=adx(6)/40.;
!        adx(8:10)=adx(8:10)/10.
      if (.not.all(adx(:).lt.huge(adx(1)))) then
         print*, 'PRGDFPMIN: avoiding derivative overflow for gradient'
         adx(:) = 0.
         idx = maxloc(abs(x))
         adx(idx) = sign(1.,x(idx))
      else if (.not.all(-adx(:).lt.huge(adx(1)))) then
         print*, 'PRGDFPMIN: avoiding derivative overflow for gradient'
         adx(:) = 0.
         idx = maxloc(abs(x))
         adx(idx) = sign(1.,x(idx))
      else
! save current x
         print*,'size of xb is: ',size(x)
         OPEN (39,file='x.b',status='unknown',form='unformatted')
         WRITE(39) x(:)
         CLOSE(39)
         OPEN (39,file='g.b',status='unknown',form='unformatted')
         WRITE(39) adx(:)
         CLOSE(39)
         ih=ih+1
         write(cname,'(I3.3)')ih
         
         cmd='mv x.b x.b_'//cname
         status = system(trim(cmd))
         cmd1='mv function.out function_'//cname
         status = system(trim(cmd1))
         cmd1='mv costfunction0000 costfunction00'//cname
         status = system(trim(cmd1))

         cmd='mv g.b g.b_'//cname
         status = system(trim(cmd))
         print*,'cmd is : ', cmd
         print*,'DPFMIN: current parameter values written to x.b'
 

     end if
       


      call flush(6)
      end subroutine df