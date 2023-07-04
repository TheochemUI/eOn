c  Dblexp:  June 92
c  Generalized to arbitrary order   Feb 92
c  This routine evaluates teh derivative of the embedding function with respect
c  to its argument (the density).

        SUBROUTINE DFRHODRHO(rho1,natm1,ipot,fprho)

        implicit real*8 (a-h,o-z)

        include 'commonblks/parameters.cmn'
        include 'commonblks/comconf.cmn'
        include 'commonblks/compotent.cmn'

        dimension fprho(MAXATOMS,3),natm1(3),rho1(MAXATOMS,3)

c     ---------------------------------------------------------------
c   Interpret parameters in potpar:

c        gfact=potpar(1,ipot)

c    ----------------------------------------------------------------

        do i = 1,natm1(ipot)
          sum = potpar(indf1,ipot)
          do io=2,iorderf
              sum=sum+potpar(indf1+io-1,ipot)*io*rho1(i,ipot)**(io-1)
          enddo

          fprho(i,ipot) = sum + potpar(1,ipot)
c                      potpar(1,ipot) is the g-transformation parameter.

        end do

        RETURN
        END
