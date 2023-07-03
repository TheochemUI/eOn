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
