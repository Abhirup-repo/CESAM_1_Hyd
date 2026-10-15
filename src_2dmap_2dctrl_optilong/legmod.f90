      module legmod
!     ************************
!     * Legendre Polynomials *
!     ************************

      use pumamod

!     polix : indirect transformation (spectral to fourier )
!     poldx :   direct transformation (fourier  to spectral)
!       x=p : polynomials
!       x=d : deviation polynomials

      real :: polip(NCSP,NLPP) ! indirect transform polynomials
      real :: polid(NCSP,NLPP) ! indirect transform deviations
      real :: plip(NLPP,NCSP),plid(NLPP,NCSP)
      real :: pliu(NCSP,NLPP),pliv(NCSP,NLPP)
      real :: pldp(NLPP,NCSP),pldd(NLPP,NCSP)
      real :: pldc(NLPP,NCSP),pldq(NLPP,NCSP)

      end module legmod

!     =================
!     SUBROUTINE LEGINI
!     =================


      subroutine legini
      use pumamod
      use legmod
!
!     initializing the legendre polynomials
!     (for a reference see eg Bourke 1988: spectral methods in global
!      climate and weather prediction models. in pysically-based
!      modelling and simulation of climate an climatic change part 1.
!      M.E.Schlesinger ed., NATO ASI series C vol.243. pp169-220.)
!

      do 500 jlat = 1 , NLPP

!        *************************
!        * set p(0,0) and p(0,1) *
!        *************************

         plip(jlat,1) = sqrt(0.5)
         f1m          = sqrt(1.5)
         plip(jlat,2) = f1m * si(jlat)
         plid(jlat,1) = 0.0
         lm           = 2

!        *************************
!        * loop over wavenumbers *
!        *************************

         do 200 m = 0 , NTRU
            e2  = sqrt(m+m+3.0)
            if (m .gt. 0) then
               f2m = -f1m * sqrt(csq(jlat)/(m+m))
               f1m =  f2m * e2
               lm  = lm + 1
               plip(jlat,lm) = f2m
               if (lm .lt. NCSP) then
                  lm = lm + 1
                  plip(jlat,lm  ) = f1m * si(jlat)
                  plid(jlat,lm-1) =  -m * plip(jlat,lm) / e2
               endif
            endif

            amsq = m * m
            do 100 n = m+2 , NTRU
               e1 = sqrt(((n-1)*(n-1) - amsq) / (4*(n-1)*(n-1)-1))
               e2 = sqrt((4*n*n-1) / (n*n-amsq))
               lm = lm + 1
               plip(jlat,lm  ) = e2 * (si(jlat)*plip(jlat,lm-1)         &
     &                         - e1 * plip(jlat,lm-2))
               plid(jlat,lm-1) = (1-n)*plip(jlat,lm)/e2+n*e1            &
     &                         * plip(jlat,lm-2)
  100       continue
            if(lm .lt. NCSP) then
             plid(jlat,lm) = -NTRU*si(jlat)*plip(jlat,lm)               &
     &                     + (NTRU+NTRU+1)*plip(jlat,lm-1)/e2
            else
             plid(jlat,lm) = -NTRU*si(jlat)*plip(jlat,lm)
            endif
  200    continue

         lm = 0
         do 400 m = 0 , NTRU
            do 300 n = m , NTRU
                 lm = lm + 1
                 znn1 = 0.0
                 if (n > 0) znn1 = 1.0 / (n*(n+1))
                 pliu(lm,jlat) = plip(jlat,lm) * znn1 * m
                 pliv(lm,jlat) = plid(jlat,lm) * znn1
                 pldp(jlat,lm) = plip(jlat,lm) * gw(jlat)
                 pldc(jlat,lm) = pldp(jlat,lm) * csq(jlat)
                 pldd(jlat,lm) = plid(jlat,lm) * gw(jlat)
                 pldq(jlat,lm) = pldp(jlat,lm) * n * (n+1) * 0.5
  300    continue
  400 continue
  500 continue

      do jlat = 1 , NLPP
         do jsp = 1 , NCSP
            polip(jsp,jlat) = plip(jlat,jsp)
            polid(jsp,jlat) = plid(jlat,jsp)
         enddo
      enddo

      return
      end

!     ================
!     SUBROUTINE GAUAW
!     ================

      subroutine gauaw(psi,pgw,k)
      real psi(k),pgw(k)
      real (kind=8) c,eps
      parameter(c = (1.0 - (2.0/3.14159265358979D0)**2)/4.0)
      parameter(eps = 1.0e-14)
      real (kind=8) pa(k),pw(k)
      real (kind=8) pk,pkm1,pkm2,pkmrk,zsp,xz

      kk = k/2
      call bsslzr(pa,kk)

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
      end

!     =================
!     SUBROUTINE BSSLZR
!     =================

      subroutine bsslzr(pbes,kdim)
      parameter(NDIM=50)
      real (kind=8) pbes(kdim)
      real (kind=8) zbes(NDIM)
      real (kind=8) PI
      parameter(PI = 3.14159265359D0)
      data zbes      /  2.4048255577D0, 5.5200781103D0,                 &
     &  8.6537279129D0,11.7915344391D0,14.9309177086D0,18.0710639679D0, &
     & 21.2116366299D0,24.3524715308D0,27.4934791320D0,30.6346064684D0, &
     & 33.7758202136D0,36.9170983537D0,40.0584257646D0,43.1997917132D0, &
     & 46.3411883717D0,49.4826098974D0,52.6240518411D0,55.7655107550D0, &
     & 58.9069839261D0,62.0484691902D0,65.1899648002D0,68.3314693299D0, &
     & 71.4729816036D0,74.6145006437D0,77.7560256304D0,80.8975558711D0, &
     & 84.0390907769D0,87.1806298436D0,90.3221726372D0,93.4637187819D0, &
     & 96.6052679510D0,99.7468198587D0,102.8883742542D0,                &
     & 106.0299309165D0,109.1714896498D0,112.3130502805D0,              &
     & 115.4546126537D0,118.5961766309D0,121.7377420880D0,              &
     & 124.8793089132D0,128.0208770059D0,131.1624462752D0,              &
     & 134.3040166383D0,137.4455880203D0,140.5871603528D0,              &
     & 143.7287335737D0,146.8703076258D0,150.0118824570D0,              &
     & 153.1534580192D0, 156.2950342685D0/

      idim = min(kdim,NDIM)
      pbes(1:idim) = zbes(1:idim)
      do j = idim+1 , kdim
         pbes(j) = pbes(j-1) + PI
      enddo
      return
      end

!     ================
!     SUBROUTINE SP2FL
!     ================

      subroutine sp2fl(psp,pfc,klev)
      use pumamod
      use legmod
      real psp(NESP,klev)
      real pfc(NHOR,klev)

      pfc(:,:) = 0.0

!$TAF LOOP = parallel
      do jlev = 1,klev
!$taf incompl pfc
         call sp2fc(psp(1,jlev),pfc(1,jlev),polip)
      enddo

      return
      end

!     ==================
!     SUBROUTINE INVLEGA
!     ==================

      subroutine invlega
      use pumamod
      use legmod

      gt(:,:) = 0.0
      gq(:,:) = 0.0
      gd(:,:) = 0.0
      gz(:,:) = 0.0
      gu(:,:) = 0.0
      gv(:,:) = 0.0
      gp(:)   = 0.0
      gpj(:)  = 0.0

!$taf LOOP = parallel
      do jlev = 1,NLEV
!$taf incompl gt
       call sp2fc(st(1,jlev),gt(1,jlev),polip)
!$taf incompl gq
       call sp2fc(sq(1,jlev),gq(1,jlev),polip)
!$taf incompl gd
       call sp2fc(sd(1,jlev),gd(1,jlev),polip)
!$taf incompl gz
       call sp2fc(sz(1,jlev),gz(1,jlev),polip)
       savez = sz(3,jlev)
       sz(3,jlev) = savez - EZ
!$taf incompl gu,gv
       call spufc(sd(1,jlev),sz(1,jlev),gu(1,jlev),gv(1,jlev),pliu,pliv)
       sz(3,jlev) = savez
      enddo

!$taf incompl gp
      call sp2fc(sp,gp ,polip)
!$taf incompl gpj
      call sp2fc(sp,gpj,polid)

      return
      end

!     ==================
!     SUBROUTINE INVLEGD
!     ==================

      subroutine invlegd
      use pumamod
      use legmod

      gt(:,:) = 0.0
      gq(:,:) = 0.0
      gu(:,:) = 0.0
      gv(:,:) = 0.0
      gp(:)   = 0.0

!$taf LOOP = parallel
      do jlev = 1,NLEV
!$taf incompl gt
       call sp2fc(st(1,jlev),gt(1,jlev),polip)
!$taf incompl gq
       call sp2fc(sq(1,jlev),gq(1,jlev),polip)
       savez = sz(3,jlev)
       sz(3,jlev) = savez - EZ
!$taf incompl gu,gv
       call spufc(sd(1,jlev),sz(1,jlev),gu(1,jlev),gv(1,jlev),pliu,pliv)
       sz(3,jlev) = savez
      enddo

!$taf incompl gp
      call sp2fc(sp,gp ,polip)

      return
      end

!     ================
!     SUBROUTINE SP2FC
!     ================

      subroutine sp2fc(psp,pfc,pol)
      use pumamod
      real psp(2,NESP/2)
      real pfc(2,NHOR/2)
      real pol(NCSP*NLPP)
!      integer :: iph,ish ! debug, FastOpt

!$taf INCOMPL pfc
!$taf loop = parallel
      do jlat =0,NLPP-1
! FastOpt         ip = 1 + jlat * NCSP    ! polynomial index
         im = 1 + jlat * NLON/2  ! fourier    index
! FastOpt         is = 1                  ! spectral   index
!$taf loop = parallel
         do jm = 0,NTRU
            jl = NTRU - jm
            ip = 1 + jlat * NCSP + (NTRU+1)*jm - (jm*(jm-1))/2 ! added, FastOpt
            is = 1 + (NTRU+1)*jm - (jm*(jm-1))/2 ! added, FastOpt
!            if (iph.ne.ip) print*,"error iph ",iph,ip,iph-ip ! debug, FastOpt
!            if (ish.ne.is) print*,"error ish ",ish,is,ish-is ! debug, FastOpt
            pfc(1,im+jm) = dot_product(pol(ip:ip+jl),psp(1,is:is+jl))
            pfc(2,im+jm) = dot_product(pol(ip:ip+jl),psp(2,is:is+jl))
! FastOpt            ip = ip + jl + 1
! FastOpt            is = is + jl + 1
         enddo ! jm
      enddo ! jlat
      return
      end

!     ================
!     SUBROUTINE DV2UV
!     ================

      subroutine dv2uv(pd,pz,pu,pv)
      use pumamod
      use legmod
      dimension pd(NESP,NLEV)
      dimension pz(NESP,NLEV)
      dimension pu(NHOR,NLEV)
      dimension pv(NHOR,NLEV)

      pu(:,:) = 0.0
      pv(:,:) = 0.0

      do jlev = 1,NLEV
       savez = pz(3,jlev)
       pz(3,jlev) = savez - EZ
!$taf incompl pu,pv
       call spufc(pd(1,jlev),pz(1,jlev),pu(1,jlev),pv(1,jlev),pliu,pliv)
       pz(3,jlev) = savez
      enddo ! jlev
      return
      end

!     ================
!     SUBROUTINE SPUFC
!     ================

      subroutine spufc(pd,pz,pu,pv,pa,pb)
      use pumamod
      real pz(2,NESP/2)   ! Vorticity
      real pd(2,NESP/2)   ! Divergence
      real pu(2,NHOR/2)   ! u-wind
      real pv(2,NHOR/2)   ! v-wind
      real pa(NCSP*NLPP)  ! Legendre polynomial pliu
      real pb(NCSP*NLPP)  ! Legendre polynomila pliv
!      integer :: iph,ish ! debug, FastOpt

!$taf incompl pu, pv
      do jlat =0,NLPP-1
! FastOpt         ip = 1 + jlat * NCSP    ! polynomial index
         im = 1 + jlat * NLON/2  ! fourier    index
! FastOpt         is = 1                  ! spectral   index
         do jm = 0,NTRU
            jl = NTRU - jm
            ip = 1 + jlat * NCSP + (NTRU+1)*jm - (jm*(jm-1))/2 ! added, FastOpt
            is = 1 + (NTRU+1)*jm - (jm*(jm-1))/2 ! added, FastOpt
!            if (iph.ne.ip) print*,"error iph ",iph,ip,iph-ip ! debug, FastOpt
!            if (ish.ne.is) print*,"error ish ",ish,is,ish-is ! debug, FastOpt
            pu(1,im+jm) = dot_product(pb(ip:ip+jl),pz(1,is:is+jl))      &
     &                   +dot_product(pa(ip:ip+jl),pd(2,is:is+jl))
            pu(2,im+jm) = dot_product(pb(ip:ip+jl),pz(2,is:is+jl))      &
     &                   -dot_product(pa(ip:ip+jl),pd(1,is:is+jl))
            pv(1,im+jm) = dot_product(pa(ip:ip+jl),pz(2,is:is+jl))      &
     &                   -dot_product(pb(ip:ip+jl),pd(1,is:is+jl))
            pv(2,im+jm) =-dot_product(pa(ip:ip+jl),pz(1,is:is+jl))      &
     &                   -dot_product(pb(ip:ip+jl),pd(2,is:is+jl))
! FastOpt            ip = ip + jl + 1
! FastOpt            is = is + jl + 1
         enddo ! jm
      enddo ! jlat
      return
      end

!     ================
!     SUBROUTINE FC2SP
!     ================

      subroutine fc2sp(fc,ps)
      use pumamod
      use legmod
      dimension ps(2,NESP/2)
      dimension fc(2,NLON/2,NLPP)

      jw = 1
      do m = 1 , NTP1
         do n = m , NTP1
            ps(1,jw) = dot_product(pldc(:,jw),fc(1,m,:))
            ps(2,jw) = dot_product(pldc(:,jw),fc(2,m,:))
            jw = jw + 1
         enddo
      enddo
      if (jw <= NESP/2) ps(:,jw:NESP/2) = 0.0
      return
      end

!     ==================
!     SUBROUTINE DIRLEGA
!     ==================

      subroutine dirlega(pdt,ptt,pqt,pzt,pft,gke,gtn,gqn,gut,guz,gvt,gvz&
     &                  ,guq,gvq,gvp)
      use pumamod
      use legmod

      dimension pdt(NESP*NLEV)
      dimension ptt(NESP*NLEV)
      dimension pzt(NESP*NLEV)
      dimension pqt(NESP*NLEV)
      dimension pft(NESP)

      dimension gke(NLON,NLPP,NLEV)
      dimension gtn(NLON,NLPP,NLEV)
      dimension gqn(NLON,NLPP,NLEV)
      dimension gut(NLON,NLPP,NLEV)
      dimension guq(NLON,NLPP,NLEV)
      dimension guz(NLON,NLPP,NLEV)
      dimension gvt(NLON,NLPP,NLEV)
      dimension gvq(NLON,NLPP,NLEV)
      dimension gvz(NLON,NLPP,NLEV)
      dimension gvp(NLON,NLPP)

      dimension zker(NLPP)
      dimension zkei(NLPP)
      dimension ztnr(NLPP)
      dimension ztni(NLPP)
      dimension zqnr(NLPP)
      dimension zqni(NLPP)
      dimension zutr(NLPP)
      dimension zuti(NLPP)
      dimension zuqr(NLPP)
      dimension zuqi(NLPP)
      dimension zuzr(NLPP)
      dimension zuzi(NLPP)
      dimension zvtr(NLPP)
      dimension zvti(NLPP)
      dimension zvqr(NLPP)
      dimension zvqi(NLPP)
      dimension zvzr(NLPP)
      dimension zvzi(NLPP)
!      integer :: jwh,jah,jbh,jrh,jih ! debug, FastOpt

      call fc2sp(gvp,pft)

      do jlev = 1 , NLEV
         jw = 1
!         jwh=1 ! debug, FastOpt
         jz = jlev * NESP
         ja = jz - NESP + 1
!         jah = jz - NESP + 1 ! debug, FastOpt
         jb = ja + 1
!         jbh = jah + 1 ! debug, FastOpt
         jr = 1
!         jrh = 1 ! debug, FastOpt
         ji = 2
!         jih = 2 ! debug, FastOpt
         do m = 0 , NTRU
            jr=2*m+1 ! FastOpt
            ji=2*m+2 ! FastOpt
!            if(jr.ne.jrh) print*,"error jr",jr,jrh,jrh-jr ! FastOpt
!            if(ji.ne.jih) print*,"error ji",ji,jih,jih-ji ! FastOpt
            zker = gke(jr,:,jlev)
            zkei = gke(ji,:,jlev)
            ztnr = gtn(jr,:,jlev)
            ztni = gtn(ji,:,jlev)
            zqnr = gqn(jr,:,jlev)
            zqni = gqn(ji,:,jlev)
            zutr = gut(jr,:,jlev)
            zuti = gut(ji,:,jlev)
            zuqr = guq(jr,:,jlev)
            zuqi = guq(ji,:,jlev)
            zuzr = guz(jr,:,jlev)
            zuzi = guz(ji,:,jlev)
            zvtr = gvt(jr,:,jlev)
            zvti = gvt(ji,:,jlev)
            zvqr = gvq(jr,:,jlev)
            zvqi = gvq(ji,:,jlev)
            zvzr = gvz(jr,:,jlev)
            zvzi = gvz(ji,:,jlev)
            do n = m , NTRU
               jw=n+m*ntru-((m-1)*m)/2+1 ! FastOpt
               ja=jz - NESP + 2*(n+m*ntru)-((m-1)*m)+1 ! FastOpt
               jb=ja+1 ! FastOpt
!               if(jw.ne.jwh) print*,"error jw",jw,jwh,jwh-jw ! FastOpt
!               if(ja.ne.jah) print*,"error ja",ja,jah,jah-ja ! FastOpt
!               if(jb.ne.jbh) print*,"error jb",jb,jbh,jbh-jb ! FastOpt
               zpd = pldd(1,jw)
               zpp = pldp(1,jw)
               zpq = pldq(1,jw)
               zduzr = zpd * zuzr(1)
               zduzi = zpd * zuzi(1)
               zdvzr = zpd * zvzr(1)
               zdvzi = zpd * zvzi(1)
               zpuzr = zpp * zuzr(1)
               zpuzi = zpp * zuzi(1)
               zpvzr = zpp * zvzr(1)
               zpvzi = zpp * zvzi(1)
               zqker = zpq * zker(1)
               zqkei = zpq * zkei(1)
               zdvtr = zpd * zvtr(1)
               zdvti = zpd * zvti(1)
               zdvqr = zpd * zvqr(1)
               zdvqi = zpd * zvqi(1)
               zputr = zpp * (ztnr(1)*csq(1)+m*zuti(1))
               zputi = zpp * (ztni(1)*csq(1)-m*zutr(1))
               zpuqr = zpp * (zqnr(1)*csq(1)+m*zuqi(1))
               zpuqi = zpp * (zqni(1)*csq(1)-m*zuqr(1))
               do jl = 2 , NLPP
                  zpd = pldd(jl,jw)
                  zpp = pldp(jl,jw)
                  zpq = pldq(jl,jw)
                  zduzr = zduzr + zpd * zuzr(jl)
                  zduzi = zduzi + zpd * zuzi(jl)
                  zdvzr = zdvzr + zpd * zvzr(jl)
                  zdvzi = zdvzi + zpd * zvzi(jl)
                  zpuzr = zpuzr + zpp * zuzr(jl)
                  zpuzi = zpuzi + zpp * zuzi(jl)
                  zpvzr = zpvzr + zpp * zvzr(jl)
                  zpvzi = zpvzi + zpp * zvzi(jl)
                  zqker = zqker + zpq * zker(jl)
                  zqkei = zqkei + zpq * zkei(jl)
                  zdvtr = zdvtr + zpd * zvtr(jl)
                  zdvti = zdvti + zpd * zvti(jl)
                  zdvqr = zdvqr + zpd * zvqr(jl)
                  zdvqi = zdvqi + zpd * zvqi(jl)
                  zputr = zputr + zpp * (ztnr(jl)*csq(jl)+m*zuti(jl))
                  zputi = zputi + zpp * (ztni(jl)*csq(jl)-m*zutr(jl))
                  zpuqr = zpuqr + zpp * (zqnr(jl)*csq(jl)+m*zuqi(jl))
                  zpuqi = zpuqi + zpp * (zqni(jl)*csq(jl)-m*zuqr(jl))
               enddo ! jl
               pdt(ja) = zqker - zduzr - m * zpvzi
               pdt(jb) = zqkei - zduzi + m * zpvzr
               pzt(ja) = zdvzr - m * zpuzi
               pzt(jb) = zdvzi + m * zpuzr
               ptt(ja) = zdvtr + zputr
               ptt(jb) = zdvti + zputi
               pqt(ja) = zdvqr + zpuqr
               pqt(jb) = zdvqi + zpuqi
!               jwh = jwh + 1
!               jah = jah + 2
!               jbh = jbh + 2
            enddo ! n
!            jrh = jrh + 2
!            jih = jih + 2
         enddo ! m
         ja=jlev*NESP - NESP + 2*(1+ntru+ntru*ntru)-((ntru-1)*ntru)+1 ! FastOpt
!         if(ja.ne.jah) print*,"error ja, pos.2",ja,jah,jah-ja ! FastOpt
         if (ja < jz) then
            pdt(ja:jz) = 0.0
            ptt(ja:jz) = 0.0
            pqt(ja:jz) = 0.0
            pzt(ja:jz) = 0.0
         endif
      enddo ! jlev

      return
      end

!     ==================
!     SUBROUTINE DIRLEGD
!     ==================

      subroutine dirlegd(ptf,pqf,ptdt,pqdt)
      use pumamod
      use legmod

      dimension ptf(NESP*NLEV)
      dimension pqf(NESP*NLEV)

      dimension ptdt(NLON,NLPP,NLEV)
      dimension pqdt(NLON,NLPP,NLEV)

      dimension ztdr(NLPP)
      dimension ztdi(NLPP)
      dimension zqdr(NLPP)
      dimension zqdi(NLPP)
!      integer :: jwh,jah,jbh,jrh,jih ! debug, FastOpt

      do jlev = 1 , NLEV
         jw = 1
!         jwh=1 ! debug, FastOpt
         jz = jlev * NESP
         ja = jz - NESP + 1
!         jah = jz - NESP + 1 ! debug, FastOpt
         jb = ja + 1
!         jbh = jah + 1 ! debug, FastOpt
         jr = 1
!         jrh = 1 ! debug, FastOpt
         ji = 2
!         jih = 2 ! debug, FastOpt
         do m = 0 , NTRU
            jr=2*m+1 ! FastOpt
            ji=2*m+2 ! FastOpt
!            if(jr.ne.jrh) print*,"error jr",jr,jrh,jrh-jr ! FastOpt
!            if(ji.ne.jih) print*,"error ji",ji,jih,jih-ji ! FastOpt
            ztdr = ptdt(jr,:,jlev)
            ztdi = ptdt(ji,:,jlev)
            zqdr = pqdt(jr,:,jlev)
            zqdi = pqdt(ji,:,jlev)
            do n = m , NTRU
               jw=n+m*ntru-((m-1)*m)/2+1 ! FastOpt
               ja=jz - NESP + 2*(n+m*ntru)-((m-1)*m)+1 ! FastOpt
               jb=ja+1 ! FastOpt
!               if(jw.ne.jwh) print*,"error jw",jw,jwh,jwh-jw ! FastOpt
!               if(ja.ne.jah) print*,"error ja",ja,jah,jah-ja ! FastOpt
!               if(jb.ne.jbh) print*,"error jb",jb,jbh,jbh-jb ! FastOpt
               zpc = pldc(1,jw)
               ztfr = zpc * ztdr(1)
               ztfi = zpc * ztdi(1)
               zqfr = zpc * zqdr(1)
               zqfi = zpc * zqdi(1)
               do jl = 2 , NLPP
                  zpc = pldc(jl,jw)
                  ztfr = ztfr + zpc * ztdr(jl)
                  ztfi = ztfi + zpc * ztdi(jl)
                  zqfr = zqfr + zpc * zqdr(jl)
                  zqfi = zqfi + zpc * zqdi(jl)
               enddo ! jl
               ptf(ja) = ztfr
               ptf(jb) = ztfi
               pqf(ja) = zqfr
               pqf(jb) = zqfi
!               jwh = jwh + 1
!               jah = jah + 2
!               jbh = jbh + 2
            enddo ! n
!            jrh = jrh + 2
!            jih = jih + 2
         enddo ! m
         ja=jlev*NESP - NESP + 2*(1+ntru+ntru*ntru)-((ntru-1)*ntru)+1 ! FastOpt
!         if(ja.ne.jah) print*,"error ja, pos.2",ja,jah,jah-ja ! FastOpt
         if (ja < jz) then
            ptf(ja:jz) = 0.0
            pqf(ja:jz) = 0.0
         endif
      enddo ! jlev

      return
      end


!     ================
!     SUBROUTINE FC3SP
!     ================

      subroutine fc3sp(fc,ps)
      use pumamod
      use legmod
      dimension ps(2,NESP/2,NLEV)
      dimension fc(2,NLON/2,NLPP,NLEV)

      do jlev = 1 , NLEV
       jw = 1
       do m = 0 , NTRU
         do n = m , NTRU
            ps(1,jw,jlev) = dot_product(pldc(:,jw),fc(1,m+1,:,jlev))
            ps(2,jw,jlev) = dot_product(pldc(:,jw),fc(2,m+1,:,jlev))
            jw = jw + 1
         enddo
       enddo
       if (jw <= NESP/2) then
        ps(:,jw:NESP/2,jlev) = 0.0
       endif
      enddo

      return
      end

!     ================
!     SUBROUTINE UV2DV
!     ================

      subroutine uv2dv(pu,pv,pd,pz)
      use pumamod
      use legmod
!
      dimension  pu(2,NLON/2,NLPP,NLEV)
      dimension  pv(2,NLON/2,NLPP,NLEV)
      dimension  pd(2,NESP/2,NLEV)
      dimension  pz(2,NESP/2,NLEV)
!
      do jlev = 1 , NLEV
         jw = 1
         do m = 0 , NTRU
           do n = m , NTRU
             pz(1,jw,jlev)=dot_product(pu(1,m+1,:,jlev),pldd(:,jw))     &
     &                    -m*dot_product(pv(2,m+1,:,jlev),pldp(:,jw))
             pz(2,jw,jlev)=dot_product(pu(2,m+1,:,jlev),pldd(:,jw))     &
     &                    +m*dot_product(pv(1,m+1,:,jlev),pldp(:,jw))
             pd(1,jw,jlev)=-dot_product(pv(1,m+1,:,jlev),pldd(:,jw))    &
     &                    -m*dot_product(pu(2,m+1,:,jlev),pldp(:,jw))
             pd(2,jw,jlev)=-dot_product(pv(2,m+1,:,jlev),pldd(:,jw))    &
     &                    +m*dot_product(pu(1,m+1,:,jlev),pldp(:,jw))
             jw = jw + 1
           enddo
        enddo
        if (jw <= NESP/2) then
         pz(:,jw:NESP/2,jlev) = 0.0
         pd(:,jw:NESP/2,jlev) = 0.0
        end if
      enddo

      return
      end

!     ================
!     SUBROUTINE SP2GP
!     ================

      subroutine sp2gp(psp,pgp,klev)
      use pumamod
      use legmod
      real psp(NESP,klev)
      real pgp(NHOR,klev)
!
      pgp(:,:) = 0.0
!
!$taf loop = parallel
      do jlev = 1,klev
!$taf incompl pgp
         call sp2fc(psp(1,jlev),pgp(1,jlev),polip)
      enddo
!
      call fc2gp(pgp,NLON,NLPP*klev)
!
      return
      end
