C This file is part of eON.
C
C SPDX-License-Identifier: BSD-3-Clause
C
C Copyright (c) 2010--present, eON Development Team
C All rights reserved.
C
C Repo:
C https://github.com/TheochemUI/eON

c   Dblexp:    June 92
c  This routine calculates FOFRHO
c  Generalized to arbitrary order

      SUBROUTINE FINT(rho1,natm1,ipot,emb)

      implicit real*8 (a-h,o-z)

      include 'commonblks/parameters.cmn'
      include 'commonblks/comconf.cmn'
      include 'commonblks/compotent.cmn'

      dimension emb(MAXATOMS,3),natm1(3),rho1(MAXATOMS,3)

      do i = 1,natm1(ipot)
        sum = 0.0
        do io=1,iorderf
          sum=sum+potpar(indf1+io-1,ipot)*rho1(i,ipot)**io
        enddo
        emb(i,ipot) = sum + potpar(1,ipot)*rho1(i,ipot)
      end do
      RETURN
      END
