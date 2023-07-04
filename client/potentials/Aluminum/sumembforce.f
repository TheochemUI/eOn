      SUBROUTINE SUMFORCE(ipot,FA1,FA2)

      implicit real*8 (a-h,o-z)

      include 'commonblks/parameters.cmn'
      include 'commonblks/comintlis.cmn'

      common/things/ isame(3),ifisame(3),nshift(3)

      dimension FA1(MAXCOO),FA2(MAXCOO)

      do 160 i = 1,natm1(ipot)-ifisame(ipot)
        do ipr = iptpr1(i,ipot)+nshift(ipot),iptpr1(i+1,ipot)-1 +
     +           nshift(ipot)
          FA1(3*i-2) = FA1(3*i-2) + rhofij(3*ipr-2)
          FA1(3*i-1) = FA1(3*i-1) + rhofij(3*ipr-1)
          FA1(3*i) = FA1(3*i) + rhofij(3*ipr)
        end do
 160  continue

      do 170 j = 1 +nshift(ipot) ,nintp(ipot) + nshift(ipot)
        FA2(3*indpra(2,j)-2) = FA2(3*indpra(2,j)-2) - rhofij(3*j-2)
        FA2(3*indpra(2,j)-1) = FA2(3*indpra(2,j)-1) - rhofij(3*j-1)
        FA2(3*indpra(2,j)) = FA2(3*indpra(2,j)) - rhofij(3*j)
 170  continue

      return
      end
