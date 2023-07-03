c   version e93:   add gregs SPF routine for chain relaxation
c   EAM version b91:   add FPI forces and potential by calling FPIforce
c
C     THIS SUBROUTINE COMPUTES THE FORCES by calling GAGAFE for each type of
c     interaction.

      SUBROUTINE POTINIT()

      implicit real*8 (a-h,o-z)

      include 'commonblks/parameters.cmn'
      include 'commonblks/combaths.cmn'
      include 'commonblks/comgeom.cmn'
      include 'commonblks/comconf.cmn'
      include 'commonblks/comenergy.cmn'
      include 'commonblks/comenperat.cmn'
      include 'commonblks/compotent.cmn'
      include 'commonblks/comluns.cmn'
      include 'commonblks/comintlis.cmn'

      common /unscale/RALOCAL(MAXCOO)

c     -----------------------------------------------------------------------
c     Defining the parameters in potpar manually (for Al-Al only)

      initflag=1

      IPOT=1
      ISAME=1
      indf1=11
      iorderf=8
      NATYPE=1
      NATMS(2)=0
      nncount=0

      alpha=90
      beta=90
      gamma=90

      IQKMIN=.false.
      nFPI=0

c      Al Parameters

      rcut(1)=7.1
      rcut2(1)=rcut(1)**2
      rskin(1)=7.25
      rskin2(1)=rskin(1)**2

      potpar(1,1)=-0.60647886749203
      potpar(2,1)=2294.3609145535
      potpar(3,1)=3.0205380362464
      potpar(4,1)=-192.06894637533
      potpar(5,1)=1.5102690181232
      potpar(6,1)=0.87928364657088
      potpar(7,1)=3.3068828537934
      potpar(8,1)=6.6137657075868
      potpar(9,1)=6.0000000000000
      potpar(10,1)=512.00000000000
      potpar(11,1)=4.4572157051836
      potpar(12,1)=193.21775368064
      potpar(13,1)=-1173.6502684704
      potpar(14,1)=4203.3500116196
      potpar(15,1)=-8785.1680280827
      potpar(16,1)=10632.994102532
      potpar(17,1)=-6921.0410455328
      potpar(18,1)=1875.2698365752

      RETURN
      END
