C This file is part of eOn.
C
C SPDX-License-Identifier: BSD-3-Clause
C
C Copyright (c) 2010--present, eOn Development Team
C All rights reserved.
C
C Repo:
C https://github.com/TheochemUI/eOn
C
!force
      subroutine FEFORCE(nm,x0,y0,z0,ispec,fx,fy,fz,pe,ax,ay,az)
      use vesin_fehe, only: vesin_fehe_build, vnl_row, vnl_jlist
          implicit  real*8 ( a-h,o-z )
c***************************************************
c
c  Ackland Many body potential for a-Fe ***********
c
      dimension x0(nm),y0(nm),z0(nm),ispec(nm)
      dimension dafrho(nm),dafrhos(nm)
      dimension fx(nm),fy(nm),fz(nm)
      common/potcut/rcutp_fe,rcutr_fe,rcut_he,rcut_her
      common/pothe/zeta1,zeta2,zeta3,zeta4,v1c,b_he(3),a_he(6),
     >             xhep(11)
      common/sband/dNs,psi,pi
      dimension etom(nm),emedtom(nm),emedtoms(nm)
c
c     write(6,796)
c 796 format(2x,'into force of Fe')
c
c      write(*,*) 'box: ',ax,ay,az
c      write(*,*) 'nm: ',nm
c      write(*,*) 'pos and input: '
c      do i = 1,nm
c        write(*,*) i,ispec(i),x0(i),y0(i),z0(i)
c      enddo
      point5 = 0.5d0
      pe = 0.0d0
      ev = 1.602177d-19
      one = 1.
      pi = 3.14159265358979324d0

      rcutp_fe = 5.3d0
      rcutr_fe = 4.2d0
      rcut_he = 3.9028d0
      rcut_her = 4.10d0
      rcut_hehe = 5.4d0

c the candidate list carries every species pair, so it spans the widest
c of the four Fe-Fe / Fe-He ranges and the He-He pair range.
      rcutmax = max(rcutp_fe,rcutr_fe,rcut_he,rcut_her,rcut_hehe)

C He-He pair potentials
      herm = 2.9683d0
      c6 = 1.35186623d0
      c8 = 0.41495143d0
      c10 = 0.17151143d0
      heaa = 186924.404d0
      hea = 10.5717543d0
      heb = -2.07758779d0
      hed = 1.438d0
      hee = 10.956*1.380658d-23/ev
c
      dNs = 20.0d0
      rb = 0.529177210818181818 ! Bohr radius
      psi = 1.5312426703208/rb

      xhep(1) = 559.804426025391
      xhep(2) = -45.916354995728
      xhep(3) = 35.550312671661
      xhep(4) = 164.319865173340
      xhep(5) = -1.727464405060
      xhep(6) = 0.106771826237
      xhep(7) = 0.073715269849
      xhep(8) = 0.038235287677
      xhep(9) = 0.220813420485
      xhep(10) = 1.367508764130
      xhep(11) = 3.382256025271

      r1 = 1.5440
      r2 = 1.6155
      r3 = 1.6896
      r4 = 1.8017
      r5 = 2.0482
      r6 = 2.3816
      r7 = 3.5067
      r8 = 3.9028
c Shifting atoms into -0.5 to 0.5
      do 102 i=1,nm
      etom(i) = 0.0d0
      fx(i) = 0.0d0
      fy(i) = 0.0d0
 102  fz(i) = 0.0d0

c
c Shifting atoms into -0.5 to 0.5
      do i = 1,nm
       x00 = mod((x0(i)/ax + 1000.5),one) - 0.5
       y00 = mod((y0(i)/ay + 1000.5),one) - 0.5
       z00 = mod((z0(i)/az + 1000.5),one) - 0.5
c changing to -0.5*ax and 0.5*ax
       x0(i) = x00*ax
       y0(i) = y00*ay
       z0(i) = z00*az
      end do
c
c      write(*,*) 'pos after -0.5 .. 0.5: '
c      do i = 1,nm
c        write(*,*) i,x0(i),y0(i),z0(i)
c      enddo
c
c vesin enumerates the candidate pairs (j > i) within rcutmax. The
c atoms sit in -0.5*ax .. 0.5*ax after the fold above, so the
c single-shot minimum image test in rhovlc and in the force loop below
c reaches the nearest image of every candidate.
c
      call vesin_fehe_build(nm,x0,y0,z0,ax,ay,az,rcutmax,iverr)
      if (iverr .ne. 0) then
        write(6,*) 'feforce: vesin neighbour list build failed'
        stop
      endif
c
      call rhovlc(nm,x0,y0,z0,ispec,ax,ay,az,
     >            dafrho,dafrhos,emedtom,emedtoms)
c
c
c half of the box in each direction, for the minimum image test
c
      axhalf = 0.5d0*ax
      ayhalf = 0.5d0*ay
      azhalf = 0.5d0*az
c
c now do the particle-particle interactions, one row of the candidate
c list per atom
c
      do 6001 ina = 1,nm
      rxi = x0(ina)
      ryi = y0(ina)
      rzi = z0(ina)
      dafim = dafrho(ina)
      dafims = dafrhos(ina)
      ssumx = 0.0d0
      ssumy = 0.0d0
      ssumz = 0.0d0
c
      do 6002 k = vnl_row(ina),vnl_row(ina+1)-1
      jna = vnl_jlist(k)
      dx=rxi-x0(jna)
      dy=ryi-y0(jna)
      dz=rzi-z0(jna)
c
c minimum image convention, rectangular cell
c
      if (dx .gt. axhalf) dx = dx - ax
      if (dx .lt.-axhalf) dx = dx + ax
      if (dy .gt. ayhalf) dy = dy - ay
      if (dy .lt.-ayhalf) dy = dy + ay
      if (dz .gt. azhalf) dz = dz - az
      if (dz .lt.-azhalf) dz = dz + az
C change rijx... into Amgstron *************
      rxij=dx
      ryij=dy
      rzij=dz
c ***** box vectors in Anmstrong ************
      rsq=rxij*rxij+ryij*ryij+rzij*rzij
C
      itype = ispec(ina) + ispec(jna)
      ipick1=itype+1
       r = dsqrt(rsq)
c
***********************************************************************
      vsum=0.0d0
      vsumd=0.0d0
      IF(ipick1.eq.1) THEN
      if (r.gt.rcutp_fe) goto 7776
c**********New ackland potential-2004*******************
       vsum =  (9.7342365892908E+03/r)
     #        *( 0.18180*exp(-2.8616724320005E+01*r)
     #          +0.50990*exp(-8.4267310396064E+00*r)
     #          +0.28020*exp(-3.6030244464156E+00*r)
     #          +0.02817*exp(-1.8028536321603E+00*r)
     #         )*HH(1.0000-r)
     #          +exp( 7.4122709384068E+00
     #               -6.4180690713367E-01*r
     #               -2.6043547961722E+00*r**2
     #               +6.2625393931230E-01*r**3
     #              )*HH(r-1.0000)*HH(2.0500-r)
     #             +(-2.7444805994228E+01*(HH(2.2-r)*(2.2-r)**3)
     #               +1.5738054058489E+01*(HH(2.3-r)*(2.3-r)**3)
     #               +2.2077118733936E+00*(HH(2.4-r)*(2.4-r)**3)
     #               -2.4989799053251E+00*(HH(2.5-r)*(2.5-r)**3)
     #               +4.2099676494795E+00*(HH(2.6-r)*(2.6-r)**3)
     #               -7.7361294129713E-01*(HH(2.7-r)*(2.7-r)**3)
     #               +8.0656414937789E-01*(HH(2.8-r)*(2.8-r)**3)
     #               -2.3194358924605E+00*(HH(3.0-r)*(3.0-r)**3)
     #               +2.6577406128280E+00*(HH(3.3-r)*(3.3-r)**3)
     #               -1.0260416933564E+00*(HH(3.7-r)*(3.7-r)**3)
     #               +3.5018615891957E-01*(HH(4.2-r)*(4.2-r)**3)
     #               -5.8531821042271E-02*(HH(4.7-r)*(4.7-r)**3)
     #               -3.0458824556234E-03*(HH(5.3-r)*(5.3-r)**3)
     #              )*HH(r-2.0500)
c                                                 Derivatives of pair potential
       vsumd = (-9.7342365892908E+03/r**2)
     #        *( 0.18180*exp(-2.8616724320005E+01*r)
     #          +0.50990*exp(-8.4267310396064E+00*r)
     #          +0.28020*exp(-3.6030244464156E+00*r)
     #          +0.02817*exp(-1.8028536321603E+00*r)
     #         )*HH(1.0000-r)
     #        +( 9.7342365892908E+03/r)
     #        *(-0.18180*exp(-2.8616724320005E+01*r)*2.8616724320005E+01
     #          -0.50990*exp(-8.4267310396064E+00*r)*8.4267310396064E+00
     #          -0.28020*exp(-3.6030244464156E+00*r)*3.6030244464156E+00
     #          -0.02817*exp(-1.8028536321603E+00*r)*1.8028536321603E+00
     #         )*HH(1.0000-r)
     #          +exp( 7.4122709384068E+00
     #               -6.4180690713367E-01*r
     #               -2.6043547961722E+00*r**2
     #               +6.2625393931230E-01*r**3
     #              )*HH(r-1.0000)*HH(2.0500-r)
     #               *(-6.4180690713367E-01
     #                 -2.6043547961722E+00*r*2.
     #                 +6.2625393931230E-01*r**2*3.
     #                )
     #             +(-2.7444805994228E+01*(HH(2.2-r)*(2.2-r)**2*3.)
     #               +1.5738054058489E+01*(HH(2.3-r)*(2.3-r)**2*3.)
     #               +2.2077118733936E+00*(HH(2.4-r)*(2.4-r)**2*3.)
     #               -2.4989799053251E+00*(HH(2.5-r)*(2.5-r)**2*3.)
     #               +4.2099676494795E+00*(HH(2.6-r)*(2.6-r)**2*3.)
     #               -7.7361294129713E-01*(HH(2.7-r)*(2.7-r)**2*3.)
     #               +8.0656414937789E-01*(HH(2.8-r)*(2.8-r)**2*3.)
     #               -2.3194358924605E+00*(HH(3.0-r)*(3.0-r)**2*3.)
     #               +2.6577406128280E+00*(HH(3.3-r)*(3.3-r)**2*3.)
     #               -1.0260416933564E+00*(HH(3.7-r)*(3.7-r)**2*3.)
     #               +3.5018615891957E-01*(HH(4.2-r)*(4.2-r)**2*3.)
     #               -5.8531821042271E-02*(HH(4.7-r)*(4.7-r)**2*3.)
     #               -3.0458824556234E-03*(HH(5.3-r)*(5.3-r)**2*3.)
     #              )*HH(r-2.0500)*(-1.)
c                                                         Derivatives of density
 7776 continue
        ELSEIF(ipick1.EQ.2) THEN
C Pair potential of He-Fe interaction
        if (r.gt.rcut_he) goto 7766
         vsum =     xhep(1)*(HH(r1-r)*(r1-r)**3)
     #             +xhep(2)*(HH(r2-r)*(r2-r)**3)
     #             +xhep(3)*(HH(r3-r)*(r3-r)**3)
     #             +xhep(4)*(HH(r4-r)*(r4-r)**3)
     #             +xhep(5)*(HH(r5-r)*(r5-r)**3)
     #             +xhep(6)*(HH(r6-r)*(r6-r)**3)
     #             +xhep(7)*(HH(r7-r)*(r7-r)**3)
     #             +xhep(8)*(HH(r8-r)*(r8-r)**3)

        vsumd =     (xhep(1)*(HH(r1-r)*(r1-r)**2*3.)
     #              +xhep(2)*(HH(r2-r)*(r2-r)**2*3.)
     #              +xhep(3)*(HH(r3-r)*(r3-r)**2*3.)
     #              +xhep(4)*(HH(r4-r)*(r4-r)**2*3.)
     #              +xhep(5)*(HH(r5-r)*(r5-r)**2*3.)
     #              +xhep(6)*(HH(r6-r)*(r6-r)**2*3.)
     #              +xhep(7)*(HH(r7-r)*(r7-r)**2*3.)
     #              +xhep(8)*(HH(r8-r)*(r8-r)**2*3.)
     #              )*(-1.)
 7766  continue
C - He-He interaction
       ELSEIF(ipick1.EQ.3) THEN
        if(r.gt.rcut_hehe) go to 7666 ! need to check
         xrm=r/herm
         sumc=c6/xrm**6.0+c8/xrm**8.0+c10/xrm**10.0
         term=heaa*exp(-hea*xrm+heb*xrm*xrm)
        if (xrm.lt.hed) then
         FFx=exp(-(hed/xrm-1.)**2.)
         vsum=hee*(term-FFx*sumc)
         vsumd=hee*(term*(-hea+2.*heb*xrm)-FFx*sumc*
     &       (2.*heD*(heD/xrm-1.)/(xrm*xrm))+FFx*(6.*c6/xrm**7.
     &       +8.*c8/xrm**9.+10.*c10/xrm**11.))/herm
        elseif (xrm.ge.hed)then
         FFx=1.0
         vsum=hee*(term-FFx*sumc)
         vsumd=hee*(term*(-hea+2.*heb*xrm)+FFx*
     &   (6.*c6/xrm**7.+8.*c8/xrm**9.+10.*c10/xrm**11.))/herm
        endif
 7666  continue
       ENDIF

        fpp=-vsumd/r
        pp=vsum
        etom(ina)=etom(ina)+pp/2.0d0
        etom(jna)=etom(jna)+pp/2.0d0
        psum=0.0d0
        psums = 0.0d0
c
c        pe = pe+pp
c
c*******************many body part***************
      if(ipick1.eq.1) then
       if (r.gt.rcutr_fe) goto 7777
        psum=1.1686859407970E+01*(HH(2.4-r)*(2.4-r)**2*3.)*(-1.)
     #        -1.4710740098830E-02*(HH(3.2-r)*(3.2-r)**2*3.)*(-1.)
     #        +4.7193527075943E-01*(HH(4.2-r)*(4.2-r)**2*3.)*(-1.)
 7777 continue
      elseif(ipick1.eq.2) then
       if(r.gt.rcut_her) goto 7781
        psums =dNs*r**2.*exp(-2.0*psi*r)*(3.-2.0*r*psi)
 7781  continue
      endif

        fcp = -psum/r
        fcps = -psums/r

       fcp = fcp * ( dafim + dafrho(jna))
       fcps = fcps * ( dafims +  dafrhos(jna))
c
c fp is (1/r)*(force on atom i from atom j)
c
      fp = fpp + fcp + fcps
      dfxk = -dx*fp
      dfyk = -dy*fp
      dfzk = -dz*fp
      ssumx = ssumx + dfxk
      ssumy = ssumy + dfyk
      ssumz = ssumz + dfzk
      fx(jna) = fx(jna) + dfxk
      fy(jna) = fy(jna) + dfyk
      fz(jna) = fz(jna) + dfzk
6002  continue
c
      fx(ina) = fx(ina) - ssumx
      fy(ina) = fy(ina) - ssumy
      fz(ina) = fz(ina) - ssumz
c
6001  continue

      do kl = 1, nm
         pe=pe+etom(kl)+emedtom(kl)+emedtoms(kl)
      end do

c      write(*,*) 'force:'
c      do i = 1, nm
c        write(*,*) i,fx(i),fy(i),fz(i)
c      enddo

      return
      end
c
      subroutine rhovlc(nm,x0,y0,z0,ispec,ax,ay,az,
     >                  dafrho,dafrhos,emedtom,emedtoms)
      use vesin_fehe, only: vnl_row, vnl_jlist
          implicit  real*8 ( a-h,o-z )
c     *****************
c
      dimension x0(nm),y0(nm),z0(nm),ispec(nm)
      dimension afrho(nm),dafrho(nm),rho(nm)
      dimension afrhos(nm),dafrhos(nm),rhos(nm)
      common/potcut/rcutp_fe,rcutr_fe,rcut_he,rcut_her
      common/pothe/zeta1,zeta2,zeta3,zeta4,v1c,b_he(3),a_he(6),
     >             xhep(11)
      common/sband/dNs,psi,pi
      dimension emedtom(nm),emedtoms(nm)
**************************************** afc nov 91 end *************
      point5 = 0.5d0
      ab2 = a/2.0d0

      do i=1,nm
       rho(i) = 0.0d0
       rhos(i) = 0.0d0
       emedtom(i)=0.0d0
       emedtoms(i)=0.0d0
       afrho(i) = 0.0d0
       dafrho(i) = 0.0d0
      enddo
c
c half of the box in each direction, for the minimum image test
c
      axhalf = 0.5d0*ax
      ayhalf = 0.5d0*ay
      azhalf = 0.5d0*az
c
c now do the particle-particle interactions, one row of the candidate
c list per atom. FEFORCE has already built the list.
c
      do 6001 ina = 1,nm
      rxi = x0(ina)
      ryi = y0(ina)
      rzi = z0(ina)
      ssumr = 0.0d0
      ssumrs = 0.0d0
c
      do 6002 k = vnl_row(ina),vnl_row(ina+1)-1
      jna = vnl_jlist(k)
      dx=rxi-x0(jna)
      dy=ryi-y0(jna)
      dz=rzi-z0(jna)
c
c minimum image convention, rectangular cell
c
      if (dx .gt. axhalf) dx = dx - ax
      if (dx .lt.-axhalf) dx = dx + ax
      if (dy .gt. ayhalf) dy = dy - ay
      if (dy .lt.-ayhalf) dy = dy + ay
      if (dz .gt. azhalf) dz = dz - az
      if (dz .lt.-azhalf) dz = dz + az
c convert to a
      rxij=dx
      ryij=dy
      rzij=dz
      rsq=rxij*rxij+ryij*ryij+rzij*rzij
      r = sqrt(rsq)
      itype = ispec(ina) + ispec(jna)
      ipick1=itype+1
c********new subroutine to calculate phi******************sjw may 91**
      phi=0.0d0
      phis = 0.0d0

      if(ipick1.eq.1) then
      if (r.gt.rcutr_fe) goto 7776
      phi = 1.1686859407970E+01*(HH(2.4-r)*(2.4-r)**3)
     #       -1.4710740098830E-02*(HH(3.2-r)*(3.2-r)**3)
     #       +4.7193527075943E-01*(HH(4.2-r)*(4.2-r)**3)
 7776 continue
      elseif(ipick1.eq.2) then
      if(r.gt.rcut_her) goto 7780
       phis = dNs*r**3.*exp(-2.0*psi*r)
 7780 continue
      endif

      ssumr = ssumr + phi
      ssumrs = ssumrs + phis
      rho(jna) = rho(jna) + phi
      rhos(jna) = rhos(jna) + phis

6002  continue
c
c rows are walked in ascending order and every neighbour has jna > ina,
c so every contribution to rho(ina) from an earlier row is already in.
c
      rho(ina) = rho(ina) + ssumr
      rhos(ina) = rhos(ina) + ssumrs
c
6001  continue

      do 110 i=1,nm
       rhoin = rho(i)
       rhosin = rhos(i)
       if(ispec(i).eq.0) then
        ipick=1
       elseif(ispec(i).eq.1) then
        ipick=2
       endif

        iw=1
        afrho(i) = embe(iw,ipick,rhoin)
        afrhos(i) = embes(iw,ipick,rhosin)
        emedtom(i) = afrho(i)
        emedtoms(i) = afrhos(i)
        iw=2
        dafrho(i) = embe(iw,ipick,rhoin)
        dafrhos(i) = embes(iw,ipick,rhosin)

      rho(i) = 0.0d0
      rhos(i) = 0.0d0
110   continue

      return
      end

c=======================================================================
c Embedding energy or its derivative
c      ro   -   density
c      iw=1 -   embe = energy;
c      iw=2 -   embe = derivative of energy
c=======================================================================
c
      real*8 function embe(iw,ipick,ro)
      implicit real*8 (a-h,o-z)
c
       if(ro.lt.1.0E-35) then
        embe = 0.0d0
         return
       endif

       if(ipick.eq.1) then

         if (iw.eq.1) embe = -ro**0.5
     #                       -6.7314115586063E-04*(ro**2)
     #                       +7.6514905604792E-08*(ro**4)
c                              Derivative of the embedding energy
c
         if (iw.eq.2) embe = -0.5/ro**0.5
     #                       -6.7314115586063E-04*(ro*2.)
     #                       +7.6514905604792E-08*(ro**3*4.)
       elseif(ipick.eq.2) then

         if (iw.eq.1) embe = 0.0d0
c
         if (iw.eq.2) embe = 0.0d0
       endif
c
       return
       end

      real*8 function embes(iw,ipick,ro)
      implicit real*8 (a-h,o-z)

      common/pothe/zeta1,zeta2,zeta3,zeta4,v1c,b_he(3),a_he(6),
     >             xhep(11)

       if(ro.lt.1.0E-35) then
         embes = 0.0d0
         return
       endif

       if(ipick.eq.1) then

         if (iw.eq.1) embes = xhep(9)*ro**0.5
     #                       +xhep(10)*(ro**2.)
     #                       +xhep(11)*(ro**4.)
c                              Derivative of the embedding energy
         if (iw.eq.2) embes = (xhep(9)*0.5)/ro**0.5
     #                        +xhep(10)*(ro*2.)
     #                        +xhep(11)*(ro**3.*4.)
       elseif(ipick.eq.2) then

         if (iw.eq.1) embes =  xhep(9)*ro**0.5
     #                        +xhep(10)*(ro**2.)
     #                        +xhep(11)*(ro**4.)
c                              Derivative of the embedding energy

         if (iw.eq.2) embes = (xhep(9)*0.5)/ro**0.5
     #                        +xhep(10)*(ro*2.)
     #                        +xhep(11)*(ro**3.*4.)
       endif

       return
       end
c
c-------------------------------------------------------------------
c Heaviside step function
c-------------------------------------------------------------------
c
      real*8 function HH (x)
      implicit real*8 (a-h,o-z)
c
       if (x.gt.0.) then
          HH = 1.d0
       else
          HH = 0.d0
       endif
       return
       end
