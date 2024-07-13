! This file is part of eOn.
!
! SPDX-License-Identifier: BSD-3-Clause
!
! Copyright (c) 2010--present, eOn Development Team
! All rights reserved.
!
! Repo:
! https://github.com/TheochemUI/eOn

module eam_routines
  use eam_dat
  implicit none
  contains

      SUBROUTINE force_eam(natms,ndim,box,R,F,U)

        IMPLICIT NONE

        INTEGER natms(2),ndim,i
        REAL*8 U(1)
        REAL*8 R(ndim),F(ndim)
        REAL*8 box(3)

        nH=natms(2)
        nCuCl=natms(1)
        nCuQ=0
        nimpo=1
        nimrp=1
        ncoo=ndim
        natoms=natms(1)+natms(2)
        PL(1)=box(1)
        PL(2)=box(2)
        PL(3)=box(3)
        do i=1,natoms
           iatclass(i)=i
        end do
        CALL potinit()


        do i=1,ndim
          ra(i)=R(i)
       end do


        call eamh2cu()



        U(1)=utot
        do i=1,ndim
          F(i)=fa(i)
        end do



      end SUBROUTINE force_eam

      subroutine eamh2cu()

      IMPLICIT NONE


      REAL*8 rho(maxatoms),dFdrho(maxatoms),dr(3),dummy
      REAL*8 rhoij,f_phi,rhoji,f_rhoij,drhoji_r,time,rij
      REAL*8 phiij,del,dphiijdr_r,drhoij_r,f_rhoji
      REAL*8 dfdrhoi,rhoi,embedi,embed,phi,sum2,fact
      INTEGER jcuq,im,icuq,jh,ih,ixyz,jcoo,iat,icoo,jcl,jat
      INTEGER  jcucl,icucl,icl



      Utot = 0.
      do icoo = 1,ncoo
         fa(icoo) = 0.
      enddo


      fact = 1./float(nimpo*nimrp)
      phi = 0.
      embed = 0.
      do iat = 1,natoms
        rho(iat) = 0.
      enddo


      do iCuCl = 1,nCuCl       ! Loop over classical Cu atoms
        icl = iCuCl
        iat = iatclass(icl)         ! Atom corresponding to classical coord icl
        do jCuCl = iCuCl+1,nCuCl    ! Loop over classical Cu atoms
          jcl = jCuCl
          jat = iatclass(jcl)
          sum2 = 0.
          do ixyz = 1,3
            icoo = 3*(iat - 1) + ixyz
            jcoo = 3*(jat - 1) + ixyz
            del = ra(jcoo) - ra(icoo)
            if (ABS(del) .gt. 0.5*PL(ixyz))del = del*(1. - PL(ixyz)/ABS(del))
            sum2 = sum2 + del**2
            dr(ixyz) = del
         enddo

          if (sum2 .le. 0.)  goto 999
          if (sum2 .lt. rcut2) then
             rij = SQRT(sum2)

            CALL EAMPhiCuCu(rij,phiij,dphiijdr_r)
            phi = phi + phiij

            CALL EAMRhoCu (rij, rhoij, dummy, 0)  ! 0:  don't calculate derivative.

            rhoji = rhoij
            rho(iat) = rho(iat) + rhoij
            rho(jat) = rho(jat) + rhoji

            do ixyz = 1,3
              icoo = 3*(iat - 1) + ixyz
              jcoo = 3*(jat - 1) + ixyz
              f_phi = dr(ixyz)*dphiijdr_r
              fa(icoo) = fa(icoo) + f_phi
              fa(jcoo) = fa(jcoo) - f_phi
            enddo
          endif
        enddo

        do im = 1,nimpo*nimrp     ! Loop over all images in FPI chains
          do jH = 1,nH            ! Loop over H atoms
            jcl = nCuCl + nCuQ + jH
            jat = iatclass(jcl) + (im - 1)
            sum2 = 0.
            do ixyz = 1,3
              icoo = 3*(iat - 1) + ixyz
              jcoo = 3*(jat - 1) + ixyz
              del = ra(jcoo) - ra(icoo)
              if (ABS(del) .gt. 0.5*PL(ixyz))del = del*(1. - PL(ixyz)/ABS(del))
              sum2 = sum2 + del**2
              dr(ixyz) = del
           enddo

            if (sum2 .le. 0.)  goto 999
            if (sum2 .lt. rcut2) then
               rij = SQRT(sum2)

              CALL EAMPhiHCu(rij,phiij,dphiijdr_r)
              phi = phi + phiij*fact

              CALL EAMRhoH (rij, rhoij, dummy, 0)

              CALL EAMRhoCu (rij, rhoji, dummy, 0)

              rho(iat) = rho(iat) + rhoij*fact
              rho(jat) = rho(jat) + rhoji
              do ixyz = 1,3
                icoo = 3*(iat - 1) + ixyz
                jcoo = 3*(jat - 1) + ixyz
                f_phi = dr(ixyz)*dphiijdr_r*fact
                fa(icoo) = fa(icoo) + f_phi
                fa(jcoo) = fa(jcoo) - f_phi
              enddo
            endif
          enddo
        enddo      ! End loop over images in MAP chains
      enddo    ! End loop over classical Cu atoms



      do im = 1,nimpo*nimrp    ! Loop again over images
        do iH = 1,nH          ! Loop using H atoms as atom i
          icl = nCuCl + nCuQ + iH
          iat = iatclass(icl) + (im - 1)

          do jH = iH+1,nH        ! Other H atoms
            jcl = nCuCl + nCuQ + jH
            jat = iatclass(jcl) + (im - 1)
            sum2 = 0.
            do ixyz = 1,3
              icoo = 3*(iat - 1) + ixyz
              jcoo = 3*(jat - 1) + ixyz
              del = ra(jcoo) - ra(icoo)
              if (ABS(del) .gt. 0.5*PL(ixyz))del = del*(1. - PL(ixyz)/ABS(del))
              sum2 = sum2 + del**2
              dr(ixyz) = del
            enddo
            if (sum2 .le. 0.)  goto 999

            if (sum2 .lt. rcut2) then
              rij = SQRT(sum2)
              CALL EAMPhiHH(rij,phiij,dphiijdr_r)
              phi = phi + phiij*fact

              CALL EAMRhoH (rij, rhoij, dummy, 0)
              rhoji = rhoij
              rho(iat) = rho(iat) + rhoij
              rho(jat) = rho(jat) + rhoji

              do ixyz = 1,3
                icoo = 3*(iat - 1) + ixyz
                jcoo = 3*(jat - 1) + ixyz
                f_phi = dr(ixyz)*dphiijdr_r*fact
                fa(icoo) = fa(icoo) + f_phi
                fa(jcoo) = fa(jcoo) - f_phi
              enddo
            endif
          enddo
        enddo        ! End loop using H atoms as atom i
      enddo       ! End loop over images

      do iCuCl = 1,nCuCl
        icl = iCuCl
        iat = iatclass(icl)
        rhoi = rho(iat)
        CALL EAMFCu(rhoi,embedi,dFdrhoi)
        embed = embed + embedi
        dFdrho(iat) = dFdrhoi
      enddo

      do im = 1,nimpo*nimrp
        do iH = 1,nH
          icl = nCuCl + nCuQ + iH
          iat = iatclass(icl) + (im - 1)
          rhoi = rho(iat)
          CALL EAMFH(rhoi,embedi,dFdrhoi)
          embed = embed + embedi*fact
          dFdrho(iat) = dFdrhoi*fact
        enddo
      enddo

      do iCuCl = 1,nCuCl
        icl = iCuCl
        iat = iatclass(icl)

        do jCuCl = iCuCl+1,nCuCl    ! Loop over classical Cu atoms
          jcl = jCuCl
          jat = iatclass(jcl)
          sum2 = 0.
          do ixyz = 1,3
            icoo = 3*(iat - 1) + ixyz
            jcoo = 3*(jat - 1) + ixyz
            del = ra(jcoo) - ra(icoo)
            if (ABS(del) .gt. 0.5*PL(ixyz))del = del*(1. - PL(ixyz)/ABS(del))
            sum2 = sum2 + del**2
            dr(ixyz) = del
          enddo

          if (sum2 .lt. rcut2) then
            rij = SQRT(sum2)
            CALL EAMRhoCu (rij, dummy, drhoij_r, 1)   ! 1:  do calculate derivative
            drhoji_r = drhoij_r
            do ixyz = 1,3
              icoo = 3*(iat - 1) + ixyz
              jcoo = 3*(jat - 1) + ixyz
              f_rhoij = dr(ixyz)*drhoij_r*dFdrho(iat)
              f_rhoji = -dr(ixyz)*drhoji_r*dFdrho(jat)
              fa(icoo) = fa(icoo) + f_rhoij - f_rhoji
              fa(jcoo) = fa(jcoo) - f_rhoij + f_rhoji
            enddo
          endif
        enddo

        do im = 1,nimpo*nimrp
          do jH = 1,nH            ! Loop over H atoms
            jcl = nCuCl + nCuQ + jH
            jat = iatclass(jcl) + (im - 1)
            sum2 = 0.
            do ixyz = 1,3
              icoo = 3*(iat - 1) + ixyz
              jcoo = 3*(jat - 1) + ixyz
              del = ra(jcoo) - ra(icoo)
              if (ABS(del) .gt. 0.5*PL(ixyz))del = del*(1. - PL(ixyz)/ABS(del))
              sum2 = sum2 + del**2
              dr(ixyz) = del
            enddo

            if (sum2 .lt. rcut2) then
              rij = SQRT(sum2)
              CALL EAMRhoH (rij, dummy, drhoij_r, 1)
              CALL EAMRhoCu (rij, dummy, drhoji_r, 1)
              drhoij_r = drhoij_r*fact
              do ixyz = 1,3
                icoo = 3*(iat - 1) + ixyz
                jcoo = 3*(jat - 1) + ixyz
                f_rhoij = dr(ixyz)*drhoij_r*dFdrho(iat)
                f_rhoji = -dr(ixyz)*drhoji_r*dFdrho(jat)
                fa(icoo) = fa(icoo) + f_rhoij - f_rhoji
                fa(jcoo) = fa(jcoo) - f_rhoij + f_rhoji
              enddo
            endif
          enddo
        enddo
      enddo

      do im = 1,nimpo*nimrp
        do iH = 1,nH
          icl = nCuCl + nCuQ + iH
          iat = iatclass(icl) + (im - 1)

          do jH = iH+1,nH
            jcl = nCuCl + nCuQ + jH
            jat = iatclass(jcl) + (im - 1)
            sum2 = 0.
            do ixyz = 1,3
              icoo = 3*(iat - 1) + ixyz
              jcoo = 3*(jat - 1) + ixyz
              del = ra(jcoo) - ra(icoo)
              if (ABS(del) .gt. 0.5*PL(ixyz))del = del*(1. - PL(ixyz)/ABS(del))
              sum2 = sum2 + del**2
              dr(ixyz) = del
            enddo

            if (sum2 .lt. rcut2) then
              rij = SQRT(sum2)
              CALL EAMRhoH (rij, dummy, drhoij_r, 1)
              drhoji_r = drhoij_r
              do ixyz = 1,3
                icoo = 3*(iat - 1) + ixyz
                jcoo = 3*(jat - 1) + ixyz
                f_rhoij = dr(ixyz)*drhoij_r*dFdrho(iat)
                f_rhoji = -dr(ixyz)*drhoji_r*dFdrho(jat)
                fa(icoo) = fa(icoo) + f_rhoij - f_rhoji
                fa(jcoo) = fa(jcoo) - f_rhoij  + f_rhoji

210   format(a15,3i6)
211   format(a30,3g13.6)
              enddo
            endif
          enddo
        enddo
      enddo


      Utot = phi + embed

201   format(i5,3g13.6)

      return

999   write(6,899) 'Collision between atoms:  ',iat,jat,' at time:  ',time
899   format(a30,2i5,a12,g13.6)
      stop

    end subroutine eamh2cu
!:::::::::::::
!  EAMFCu.f
!:::::::::::::
      subroutine  EAMFCu (rhotot, embedE, dFdrho)

      IMPLICIT NONE
      REAL*8 rhotot2,rhotot3,rhotot4,rhotot5,rhotot6,rhotot7,rhotot8,embede,rhotot,dfdrho

      rhotot2 = rhotot**2
      rhotot3 = rhotot2*rhotot
      rhotot4 = rhotot3*rhotot
      rhotot5 = rhotot4*rhotot
      rhotot6 = rhotot5*rhotot
      rhotot7 = rhotot6*rhotot
      rhotot8 = rhotot7*rhotot

      embedE = FCu1*rhotot + FCu2*rhotot2 + FCu3*rhotot3 + FCu4*rhotot4 + FCu5*rhotot5 + FCu6*rhotot6 + FCu7*rhotot7 + FCu8*rhotot8
      dFdrho = FCu1 + 2.*FCu2*rhotot + 3.*FCu3*rhotot2 + 4.*FCu4*rhotot3 + 5.*FCu5*rhotot4 + 6.*FCu6*rhotot5 + 7.*FCu7*rhotot6 + &
           & 8. * FCu8*rhotot7

      return
    end subroutine EAMFCu
!:::::::::::::
!  EAMFH.f
!:::::::::::::
      subroutine  EAMFH (rhotot, embedE, dFdrho)

      IMPLICIT NONE
      REAL*8 rhotot2,rhotot3,rhotot4,rhotot5,rhotot,embede,dfdrho

      rhotot2 = rhotot**2
      rhotot3 = rhotot2*rhotot
      rhotot4 = rhotot3*rhotot
      rhotot5 = rhotot4*rhotot

      embedE = FH1*rhotot + FH2*rhotot2 + FH3*rhotot3 + FH4*rhotot4 + FH5*rhotot5
      dFdrho = FH1 + 2.*FH2*rhotot + 3.*FH3*rhotot2 + 4.*FH4*rhotot3 + 5.*FH5*rhotot4

      return
    end subroutine EAMFH
!:::::::::::::
!  EAMPhiCuCu.f
!:::::::::::::
      subroutine EAMPhiCuCu (R, phiij, dphidroverr)


      IMPLICIT NONE
      REAL*8 term1,term2,phiij,dphidroverr,R

      term1 = DACuCu*EXP(-alphaACuCu*R)
      term2 = DBCuCu*EXP(-alphaBCuCu*R)
      phiij = term1 + term2 - phicutCuCu
      dphidroverr = -(alphaACuCu*term1 + alphaBCuCu*term2)/R

      return
    end subroutine EAMPhiCuCu
    !:::::::::::::
    !  EAMPhiHCu.f
    !:::::::::::::
      subroutine EAMPhiHCu (R, phiij, dphidroverr)


      IMPLICIT NONE
      REAL*8 term1,term2,phiij,dphidroverr,R

      term1 = DAHCu*EXP(-alphaAHCu*R)
      term2 = DBHCu*EXP(-alphaBHCu*R)
      phiij = term1 + term2 - phicutHCu
      dphidroverr = -(alphaAHCu*term1 + alphaBHCu*term2)/R

      return
    end subroutine EAMPhiHCu
!:::::::::::::
!  EAMPhiHH.f
!:::::::::::::
      subroutine EAMPhiHH (R, phiij, dphidroverr)


      IMPLICIT NONE
      REAL*8 term1,term2,phiij,dphidroverr,R

      term1 = DAHH*EXP(-alphaAHH*R)
      term2 = DBHH*EXP(-alphaBHH*R)
      phiij = term1 + term2 - phicutHH
      dphidroverr = -(alphaAHH*term1 + alphaBHH*term2)/R

      return
    end subroutine EAMPhiHH
!:::::::::::::
!  EAMRhoCu.f
!:::::::::::::
      subroutine EAMRhoCu (R, rho, drhodroverr, indic)


      IMPLICIT NONE
      REAL*8 R,term1,rinv,term2,term3,rho,drhodroverr
      INTEGER indic

      Rinv = 1.0/R
      term1 = EXP(-betaACu*R)
      term2 = gammaCu*EXP(-betaBCu*R)
      term3 = scaleCu * R**netaCu

      rho = term3 * (term1 + term2)
      if (indic .ne. 0) drhodroverr = (netaCu*rho*Rinv - term3 * (betaACu*term1 + betaBCu*term2))*Rinv
      rho = rho - rhocutCu

      return
    end subroutine EAMRhoCu
!:::::::::::::
!  EAMRhoH.f
!:::::::::::::
      subroutine EAMRhoH (R, rho, drhodroverr, indic)


      IMPLICIT NONE
      REAL*8 R,term1,rinv,term2,term3,rho,drhodroverr
      INTEGER indic

      Rinv = 1.0/R
      term1 = EXP(-betaAH*R)
      term2 = gammaH*EXP(-betaBH*R)
      term3 = scaleH * R**netaH

      rho = term3 * (term1 + term2)
      if (indic .ne. 0) drhodroverr = (netaH*rho*Rinv - term3 * (betaAH*term1 + betaBH*term2))*Rinv
      rho = rho - rhocutH

      return
    end subroutine EAMRhoH

      SUBROUTINE potinit()


        IMPLICIT NONE
        REAL*8 phitemp,dphidr_r,rhotemp,drhodr_r,potpar(18,3)
        INTEGER i

        rcut = 6.10000
        rskin = 6.30000
        rcut2=rcut**2
        rskin2=rskin**2
        ! I don't want to schlep around pot.par
        ! So:
        ! aa = np.loadtxt("pot.par", skiprows=1)
        ! for idx, val in enumerate(aa, start=1):
        !   print(f"potpar({idx}, 1) = {val[0]}\npotpar({idx}, 2) = {val[1]}\npotpar({idx}, 3) = {val[2]}")
        potpar(1, 1) = 0.0
        potpar(1, 2) = 0.0
        potpar(1, 3) = 0.0
        potpar(2, 1) = 2862.3
        potpar(2, 2) = 79.5013
        potpar(2, 3) = 86.1495
        potpar(3, 1) = 3.51236
        potpar(3, 2) = 2.47961
        potpar(3, 3) = 4.21154
        potpar(4, 1) = -109.107
        potpar(4, 2) = -107.555
        potpar(4, 3) = 15355.5
        potpar(5, 1) = 1.75618
        potpar(5, 2) = 2.99918
        potpar(5, 3) = 6.07604
        potpar(6, 1) = 0.273072
        potpar(6, 2) = 2.14389
        potpar(6, 3) = 0.0
        potpar(7, 1) = 3.69051
        potpar(7, 2) = 3.77701
        potpar(7, 3) = 0.0
        potpar(8, 1) = 7.38101
        potpar(8, 2) = 0.0
        potpar(8, 3) = 0.0
        potpar(9, 1) = 6.0
        potpar(9, 2) = 0.0
        potpar(9, 3) = 0.0
        potpar(10, 1) = 512.0
        potpar(10, 2) = 0.0
        potpar(10, 3) = 0.0
        potpar(11, 1) = -112.945
        potpar(11, 2) = -81.7537
        potpar(11, 3) = 0.0
        potpar(12, 1) = 8510.04
        potpar(12, 2) = 838.668
        potpar(12, 3) = 0.0
        potpar(13, 1) = -261734.0
        potpar(13, 2) = -3952.35
        potpar(13, 3) = 0.0
        potpar(14, 1) = 4780090.0
        potpar(14, 2) = 8767.87
        potpar(14, 3) = 0.0
        potpar(15, 1) = -52341900.0
        potpar(15, 2) = -6599.37
        potpar(15, 3) = 0.0
        potpar(16, 1) = 339124000.0
        potpar(16, 2) = 0.0
        potpar(16, 3) = 0.0
        potpar(17, 1) = -1201500000.0
        potpar(17, 2) = 0.0
        potpar(17, 3) = 0.0
        potpar(18, 1) = 1796190000.0
        potpar(18, 2) = 0.0
        potpar(18, 3) = 0.0

        FCu1 = potpar(11,1)
        FCu2 = potpar(12,1)
        FCu3 = potpar(13,1)
        FCu4 = potpar(14,1)
        FCu5 = potpar(15,1)
        FCu6 = potpar(16,1)
        FCu7 = potpar(17,1)
        FCu8 = potpar(18,1)

        FH1 = potpar(11,2)
        FH2 = potpar(12,2)
        FH3 = potpar(13,2)
        FH4 = potpar(14,2)
        FH5 = potpar(15,2)

        DACuCu = potpar(2,1)
        alphaACuCu = potpar(3,1)
        DBCuCu = potpar(4,1)
        alphaBCuCu = potpar(5,1)

        DAHH = potpar(2,2)
        alphaAHH = potpar(3,2)
        DBHH = potpar(4,2)
        alphaBHH = potpar(5,2)

        DAHCu = potpar(2,3)
        alphaAHCu = potpar(3,3)
        DBHCu = potpar(4,3)
        alphaBHCu = potpar(5,3)

        scaleCu = potpar(6,1)
        betaACu = potpar(7,1)
        betaBCu = potpar(8,1)
        netaCu = potpar(9,1) + 0.01d0
        gammaCu = potpar(10,1)

        scaleH = potpar(6,2)
        betaAH = potpar(7,2)
        betaBH = potpar(8,2)
        netaH = potpar(9,2) + 0.01d0
        gammaH = potpar(10,2)

        phicutCuCu = 0.d0
        phicutHCu  = 0.d0
        phicutHH   = 0.d0
        rhocutCu   = 0.d0
        rhocutH    = 0.d0
        CALL EAMPhiCuCu(rcut,phitemp,dphidr_r)
        phicutCuCu = phitemp
        dphiCuCudr_r = dphidr_r
        CALL EAMPhiHCu(rcut,phitemp,dphidr_r)
        phicutHCu = phitemp
        dphiHCudr_r = dphidr_r
        CALL EAMPhiHH(rcut,phitemp,dphidr_r)
        phicutHH = phitemp
        dphiHHdr_r = dphidr_r
        CALL EAMRhoCu(rcut,rhotemp,drhodr_r,1)
        rhocutCu = rhotemp
        drhoCudr_r = drhodr_r
        CALL EAMRhoH(rcut,rhotemp,drhodr_r,1)
        rhocutH = rhotemp
        drhoHdr_r = drhodr_r

      END SUBROUTINE potinit

end module eam_routines
