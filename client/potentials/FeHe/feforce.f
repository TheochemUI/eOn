!force
      subroutine FEFORCE(nm,x0,y0,z0,ispec,fx,fy,fz,pe,ax,ay,az)
          implicit  real*8 ( a-h,o-z )
          integer*4       tv0001
c***************************************************
c
      parameter (nmat=3,nnbrs=800)
c      parameter (nlcx=7,nlcy=7,nlcz=7)
c      parameter (nlcx=5,nlcy=5,nlcz=5)
      parameter (nlcx=3,nlcy=3,nlcz=3)
c      parameter (nlcx=1,nlcy=1,nlcz=1)
      parameter (nlc=nlcx*nlcy*nlcz)
c  Ackland Many body potential for a-Fe ***********
      integer*4 nix(27),niy(27),niz(27)
c Data block
      data nix  / 0,-1,-1,-1,0,0,-1,1,-1, 0, 1,-1,0,1,
     > 1, 1, 1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1/
      data niy  / 0, 0,-1, 1,1,0, 0,0,-1,-1,-1, 1,1,1,
     > 1, 0,-1,-1, 0, 0, 0,-1,-1,-1, 1, 1, 1/
      data niz  / 0, 0, 0, 0,0,1, 1,1, 1, 1, 1, 1,1,1,
     > 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/
c
      dimension x0(nm),y0(nm),z0(nm),ispec(nm)
      dimension afrho(nm),dafrho(nm),dafnbr(nnbrs),rho(nm),
     >       drho(nnbrs)
      dimension afrhos(nm),dafrhos(nm),dafnbrs(nnbrs),rhos(nm),
     >      drhos(nnbrs)
      dimension fx(nm),fy(nm),fz(nm),dfx(nnbrs),dfy(nnbrs),dfz(nnbrs)
     > ,fb(nmat,nmat)
      common/potcut/rcutp_fe,rcutr_fe,rcut_he,rcut_her
      common/pothe/zeta1,zeta2,zeta3,zeta4,v1c,b_he(3),a_he(6),
     >             xhep(11)
      common/sband/dNs,psi,pi
      dimension etom(nm),emedtom(nm),emedtoms(nm)
      dimension ltop(nlc),jaddr(nnbrs),link(nm)
      dimension xnbr(nnbrs),ynbr(nnbrs),znbr(nnbrs)
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
      call linvlc(nm,x0,y0,z0,ax,ay,az,ltop,link)
c
c
      call rhovlc(nm,x0,y0,z0,ispec,ax,ay,az,ltop,link,
     >            dafrho,dafrhos,emedtom,emedtoms)
c
c
c modifications for the vectorised link cell method follow.
c
      ix = 1
      iy = 1
      iz = 1
      s = 1.0d0
c
c primary loop over link cells
c
      do 600 ic = 1,nlc
      i = ltop(ic)
c
c bypass this cell if it is empty
c the goto 599 corrects a bug in the original heyes-smith formulation.
c the latter would have goto 600 thereby missing out the ix,iy,iz update
c
      if (i.eq.0) goto 599
      m = 0
99    m = m+1
c ******* jmh dec 88 : this statement moved from later ****************
c ******* point in routine - from alan foreman at harwell *************
      if (m.gt.nnbrs) goto 999
c *************** jmh dec 88 : end ************************************
      jaddr(m) = i
      xnbr(m) = x0(i)
      ynbr(m) = y0(i)
      znbr(m) = z0(i)
      dafnbr(m) = dafrho(i)
      dafnbrs(m) = dafrhos(i)
      i = link(i)
      if (i.gt.0) goto 99
      mstart = m
c
c secondary loop over neighbouring cells
c
      do 4001 kc = 2,14
      sx = 0.0d0
      sy = 0.0d0
      sz = 0.0d0
      jx = ix + nix(kc)
      jy = iy + niy(kc)
      jz = iz + niz(kc)
c
c minimum image conversion
c
      if ((ix.eq.nlcx).and.(jx.gt.ix)) then
      jx = 1
      sx = s
      elseif ((ix.eq.1).and.(jx.lt.ix)) then
      jx = nlcx
      sx = -s
      endif
      if ((iy.eq.nlcy).and.(jy.gt.iy)) then
      jy = 1
      sy = s
      elseif ((iy.eq.1).and.(jy.lt.iy)) then
      jy = nlcy
      sy = -s
      endif
      if ((iz.eq.nlcz).and.(jz.gt.iz)) then
      jz = 1
      sz = s
      elseif ((iz.eq.1).and.(jz.lt.iz)) then
      jz = nlcz
      sz = -s
      endif
c
c index of neighbouring cell
c
      jc = jx + nlcx*( (jy-1) + nlcy*(jz-1) )
      j  = ltop(jc)
c
c bypass this neighbouring cell if it is empty
c
      if (j.eq.0) goto 4001
199   m = m+1
c ******* jmh dec 88 : this statement moved from later ****************
c ******* point in routine - from alan foreman at harwell *************
      if (m.gt.nnbrs) goto 999
c *************** jmh dec 88 : end ************************************
      jaddr(m) = j
      dafnbr(m) = dafrho(j)
      dafnbrs(m) = dafrhos(j)
      xnbr(m) = x0(j) + sx*ax
      ynbr(m) = y0(j) + sy*ay
      znbr(m) = z0(j) + sz*az
      j = link(j)
      if (j.gt.0) goto 199
4001  continue
c
c we have now found all the neighbouring particles of cell ic.
c
      max = m
c
      tv0001  = mstart
      if (max.eq.mstart) tv0001  = mstart-1
      if (tv0001 .le.0) goto 599
c
c now do the particle-particle interactions
c
      do 6001 im = 1,tv0001
      ina= jaddr(im)
      rxi = xnbr(im)
      ryi = ynbr(im)
      rzi = znbr(im)
      dafim = dafnbr(im)
      dafims = dafnbrs(im)
      mm  = 0
c do the inner vectorised loop
c
      do 6002 m = im+1,max
      jna=jaddr(m)
      mm = mm+1
      dx=rxi-xnbr(m)
      dy=ryi-ynbr(m)
      dz=rzi-znbr(m)
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
        if(r.gt.5.4d0) go to 7666 ! need to check
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

       fcp = fcp * ( dafim + dafnbr(m))
       fcps = fcps * ( dafims +  dafnbrs(m))
c
c fp is (1/r)*(force on atom i from atom j)
c
      fp = fpp + fcp + fcps
      dfx(mm) = -dx*fp
      dfy(mm) = -dy*fp
      dfz(mm) = -dz*fp
6002  continue
c
      mmax = mm
c
c**************** jmh nov '88 - begin *************************
c******* replace cray intrinsic function ssum by : ************
          ssumx = 0.d0
          ssumy = 0.d0
          ssumz = 0.d0
          do 99000 is = 1,mmax
          ssumx = ssumx + dfx(is)
          ssumy = ssumy + dfy(is)
          ssumz = ssumz + dfz(is)
99000     continue
      i = jaddr(im)
c******* cray function ssum replaced in statements below : ****
      fx(i) = fx(i) - ssumx
      fy(i) = fy(i) - ssumy
      fz(i) = fz(i) - ssumz
c****************** jmh nov '88 - end *************************
c
c this loopis now vectorised.
c
*vocl loop,novrec
      do 3422 kk=1,mmax
      j = jaddr(im+kk)
      fx(j) = fx(j) + dfx(kk)
      fy(j) = fy(j) + dfy(kk)
      fz(j) = fz(j) + dfz(kk)
3422  continue
c
6001  continue
c
c primary cell index update
c
599   continue
      ix = ix+1
      if (ix.gt.nlcx) then
      ix = 1
      iy = iy+1
      if (iy.gt.nlcy) then
      iy = 1
      iz = iz+1
      endif
      endif
600   continue

      do kl = 1, nm
         pe=pe+etom(kl)+emedtom(kl)+emedtoms(kl)
      end do

c      write(*,*) 'force:'
c      do i = 1, nm
c        write(*,*) i,fx(i),fy(i),fz(i)
c      enddo

      return
999   write(idout,1004) m,nnbrs
1004  format(' sorry user, the run is stopping now in "kravlc"',/,
     >       ' because m=',i10,' which is greater than nnbrs=',i5,'.')
      stop
      end
c
      subroutine linvlc(nm,x0,y0,z0,ax,ay,az,ltop,link)
          implicit  real*8 ( a-h,o-z )
c     *****************
c
c sets up the link cell map
c
c      parameter (nlcx=7,nlcy=7,nlcz=7)
c      parameter (nlcx=5,nlcy=5,nlcz=5)
      parameter (nlcx=3,nlcy=3,nlcz=3)
c      parameter (nlcx=1,nlcy=1,nlcz=1)
      parameter (nlc=nlcx*nlcy*nlcz)
      parameter (hmeps=0.5d0-1d-9)
      dimension x0(nm),y0(nm),z0(nm)
      dimension ltop(nlc),link(nm)
c     write(6,300)
c300  format( '    enter subroutine linvlc to form linkage map')
      fnlcx = dfloat(nlcx)
      fnlcy = dfloat(nlcy)
      fnlcz = dfloat(nlcz)
      do 100 i=1,nlc
      ltop(i)=0
100   continue
c
c determine in which link cell atom i is situated.
c the parameter 0.5-epsilon is used instead of 0.5 (which was in the
c original program) for the assignment to link cells in case the int
c function operates on an atom exactly at the box boundary, when in
c the original version it would return an integer one greater than
c the correct one.
c
c      write(*,*) 'link cells:'
      do 110 i=1,nm
      xx0 = x0(i)/ax
      yy0 = y0(i)/ay
      zz0 = z0(i)/az
      ix = int( (xx0+hmeps)*fnlcx     ) + 1
      iy = int( (yy0+hmeps)*fnlcy     )
      iz = int( (zz0+hmeps)*fnlcz     )
c      write(*,*) i,':',ix,iy,iz
      ip = ix + nlcx*iy + nlcx*nlcy*iz
c
c assign atom i to link cell ip
c
      j = ltop(ip)
      ltop(ip) = i
      link(i) = j
110   continue
      return
      end
c
      subroutine rhovlc(nm,x0,y0,z0,ispec,ax,ay,az,ltop,link,
     >                  dafrho,dafrhos,emedtom,emedtoms)
          implicit  real*8 ( a-h,o-z )
          integer*4       tv0001
c     *****************
      parameter (nmat=3,nnbrs=800)
c      parameter (nlcx=7,nlcy=7,nlcz=7)
c      parameter (nlcx=5,nlcy=5,nlcz=5)
      parameter (nlcx=3,nlcy=3,nlcz=3)
c      parameter (nlcx=1,nlcy=1,nlcz=1)
      parameter (nlc=nlcx*nlcy*nlcz)
      integer*4 nix(27),niy(27),niz(27)
cdata block
      data nix  / 0,-1,-1,-1,0,0,-1,1,-1, 0, 1,-1,0,1,
     > 1, 1, 1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1/
      data niy  / 0, 0,-1, 1,1,0, 0,0,-1,-1,-1, 1,1,1,
     > 1, 0,-1,-1, 0, 0, 0,-1,-1,-1, 1, 1, 1/
      data niz  / 0, 0, 0, 0,0,1, 1,1, 1, 1, 1, 1,1,1,
     > 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/
c
      dimension x0(nm),y0(nm),z0(nm),ispec(nm)
      dimension afrho(nm),dafrho(nm),dafnbr(nnbrs),rho(nm),
     >             drho(nnbrs)
      dimension afrhos(nm),dafrhos(nm),dafnbrs(nnbrs),rhos(nm),
     >      drhos(nnbrs)
      common/potcut/rcutp_fe,rcutr_fe,rcut_he,rcut_her
      common/pothe/zeta1,zeta2,zeta3,zeta4,v1c,b_he(3),a_he(6),
     >             xhep(11)
      common/sband/dNs,psi,pi
      dimension emedtom(nm),emedtoms(nm)
      dimension ltop(nlc),jaddr(nnbrs),link(nm)
      dimension xnbr(nnbrs),ynbr(nnbrs),znbr(nnbrs)
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
c modifications for the vectorised link cell method follow.
      ix = 1
      iy = 1
      iz = 1
      s  = 1.0d0
c
c primary loop over link cells
c
      do 600 ic = 1,nlc
      i = ltop(ic)
c
c bypass this cell if it is empty
c
      if (i.eq.0) goto 599
      m = 0
99    m = m+1
c ******* jmh dec 88 : this statement moved from later ****************
c ******* point in routine - from alan foreman at harwell *************
cRES3  info from HFD via AFC at Liverpool                               RES3
cRES3 if (max.gt.nnbrs) goto 999                                        RES3
      if (m.gt.nnbrs) goto 999                                          RES3
c *************** jmh dec 88 : end ************************************
      jaddr(m) = i
      xnbr(m) = x0(i)
      ynbr(m) = y0(i)
      znbr(m) = z0(i)
      i = link(i)
      if (i.gt.0) goto 99
      mstart = m
c
c secondary loop over neighbouring cells
c
      do 4001 kc = 2,14
      sx = 0.0d0
      sy = 0.0d0
      sz = 0.0d0
      jx = ix + nix(kc)
      jy = iy + niy(kc)
      jz = iz + niz(kc)
c
c minimum image conversion
c
      if ((ix.eq.nlcx).and.(jx.gt.ix)) then
      jx = 1
      sx = s
      elseif ((ix.eq.1).and.(jx.lt.ix)) then
      jx = nlcx
      sx = -s
      endif
      if ((iy.eq.nlcy).and.(jy.gt.iy)) then
      jy = 1
      sy = s
      elseif ((iy.eq.1).and.(jy.lt.iy)) then
      jy = nlcy
      sy = -s
      endif
      if ((iz.eq.nlcz).and.(jz.gt.iz)) then
      jz = 1
      sz = s
      elseif ((iz.eq.1).and.(jz.lt.iz)) then
      jz = nlcz
      sz = -s
      endif
c
c index of neighbouring cell
c
      jc = jx + nlcx*( (jy-1) + nlcy*(jz-1) )
      j  = ltop(jc)
c
c bypass this neighbouring cell if it is empty
c
      if (j.eq.0) goto 4001
199   m = m+1
c ******* jmh dec 88 : this statement moved from later ****************
c ******* point in routine - from alan foreman at harwell *************
      if (m.gt.nnbrs) goto 999                                          RES3
c *************** jmh dec 88 : end ************************************
      jaddr(m) = j
      xnbr(m) = x0(j) + sx*ax
      ynbr(m) = y0(j) + sy*ay
      znbr(m) = z0(j) + sz*az
      j = link(j)
      if (j.gt.0) goto 199
4001  continue
c
c we have now found all the neighbouring particles of cell ic.
c
      max = m
c *** if (max.gt.nnbrs) goto 999 *** jmh dec 88 : too late! ***********
c *** so moved to earlier points in routine - from alan ***************
c *** foreman at harwell. jmh dec 88 : end **<<  SEE RES3  >>*********
      tv0001  = mstart
      if (max.eq.mstart) tv0001  = mstart-1
      if (tv0001 .le.0) goto 599
c
c now do the particle-particle interactions
c
      do 6001 im = 1,tv0001
      ina= jaddr(im)
      rxi = xnbr(im)
      ryi = ynbr(im)
      rzi = znbr(im)
      mm  = 0
c
c do the inner vectorised loop
c
      do 6002 m = im+1,max
      jna=jaddr(m)
      mm = mm+1
      dx=rxi-xnbr(m)
      dy=ryi-ynbr(m)
      dz=rzi-znbr(m)
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

      drho(mm) = phi
      drhos(mm) = phis

6002  continue
c
      mmax = mm
c**************** jmh nov '88 - begin *************************
c******* replace cray intrinsic function ssum by : ************
          ssumr = 0.d0
          ssumrs = 0.0d0
          do 99002 is = 1,mmax
          ssumr = ssumr + drho(is)
          ssumrs = ssumrs + drhos(is)
99002     continue
      i = jaddr(im)
c******* cray function ssum replaced in statement  below : ****
      rho(i) = rho(i) + ssumr
      rhos(i) = rhos(i) + ssumrs
c****************** jmh nov '88 - end *************************

c
c this loop is now vectorised! :
*vocl loop,novrec
      do 3422 kk=1,mmax
      j = jaddr(im+kk)
      rho(j) = rho(j) + drho(kk)
      rhos(j) = rhos(j) + drhos(kk)
3422  continue
c
6001  continue
c
c primary cell index update
c
599   continue
      ix = ix+1
      if (ix.gt.nlcx) then
      ix = 1
      iy = iy+1
      if (iy.gt.nlcy) then
      iy = 1
      iz = iz+1
      endif
      endif
600   continue

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
999   write(6,1004) m,nnbrs                                             RES3
1004  format(' sorry user, the run is stopping now in "rhovlc"',/,      RES3
     >       ' because m=',i10,' which is greater than nnbrs=',i5,'.',  RES3
     > /,' increase either nnbrs or nlcx in the parameter statements')  RES3
      stop
c ***************** jmh dec 88 : end ******************************
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
