c   EAM version e:
c   EAM version d:    Greg found bug in interactions between FPI chains.
c   EAM version c:    tags 100-199 for fixed-endpoint FPI
c                       -  200-299 for circular FPI with no constraints
c                       -  300-399  -     -      -   -   centroid fixed
c       Potential Dblexp:   June 92.   Similar to Generalized Voter pot., except
c                            phi and rho are a double exponentials,
c                            Also, F is evaluated with straight rho
c                           (including scale).
c   EAM version b92:
c                     Generalize potential, make order of F(rho) arbitrary
c                     add FPI.  Fix neighborlist for FPI images
c                     In neighbor list for each FPI image, all reference
c                     to other images of same FPI is excluded and
c                     reference to images of other FPI that do not correspond
c                     to the same imaginary time.
c                     Scale contributions of FPI images to pair pot and density
c                     by the number of images in each FPI chain (dispersion).
c
c      error
c
C     SUBROUTINE TO COMPUTE FORCE, ENERGY AND VIRIAL FOR A GROUP OF
C     ATOMS INTERACTING WITH ANOTHER GROUP OF ATOMS VIA L-J POTENTIAL.
C     THE GROUPS MAY BE THE SAME.
c     The force is scaled in such a way that the real force is ax*fa.
C
      SUBROUTINE GAGAFE(NATOM1,RA1,FA1,NATOM2,RA2,FA2,UTOT,virial,
     +                  IPOT,ISAME)

c-------------------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      logical IsatiFPI,IsatjFPI

      include 'commonblks/parameters.cmn'
!      include 'commonblks/comtime.cmn'
      include 'commonblks/comgeom.cmn'
      include 'commonblks/combaths.cmn'
      include 'commonblks/compotent.cmn'
      include 'commonblks/comconf.cmn'
      include 'commonblks/comenperat.cmn'
      include 'commonblks/comluns.cmn'
      include 'commonblks/comintlis.cmn'

      common /unscale/ RALOCAL(MAXCOO)
c
      DIMENSION DELTA(MAXCOO),phi(MAXPRS),phivirst(MAXPRS)
      DIMENSION RA1(MAXCOO),RA2(MAXCOO),FA1(MAXCOO),FA2(MAXCOO)

c ------------------------------------------------------------------------------

      lunout=6

      natm1(ipot) = natom1
      natm2(ipot) = natom2

      axhalf=0.5*ax
      ayhalf=0.5*ay
      azhalf=0.5*az

      sccut2=rcut2(ipot)
      sccut=rcut(ipot)
      skinr2=rskin2(ipot)

      iptcom(1)=1
      ishift=iptcom(ipot)-1

      ifisam=1
      if(indupd .ne. 1) go to 600

c      -------------------------------------------------------------------

c  Update the list of interacting pairs:
      if(ipot .eq. mxipot) indupd=0
c            if this is the last type of interaction then the list of
c            interactions will not have to be updated in the next call.
      INDX=0
      nintpr=0
c                 nintpr counts the total number of interacting pairs.
      DO 500 I=1,NATOM1
        iabs=i+iatshift1
        JSTRT=1
        IF(ISAME.EQ.1) JSTRT=I+1
        NJATMS=NATOM2-JSTRT+1
c
        DO 300 J=1,NJATMS
          DELTA(3*(J-1)+1)=RA1(INDX+1)-RA2(3*(J+JSTRT-2)+1)
          DELTA(3*(J-1)+2)=RA1(INDX+2)-RA2(3*(J+JSTRT-2)+2)
          DELTA(3*(J-1)+3)=RA1(INDX+3)-RA2(3*(J+JSTRT-2)+3)
300     CONTINUE

C     APPLY PERIODIC BOUNDARY CONDITIONS TO DELTA ARRAY  (only in x-y plane):

c  Rectangular simulation cell:
      DO 320 J=1,NJATMS
c   First the y coordinate:
          IF(DELTA(3*(j-1)+2) .GT. ayhalf)
     +       DELTA(3*(j-1)+2)=DELTA(3*(j-1)+2)-ay
          IF(DELTA(3*(j-1)+2) .LT.-ayhalf)
     +       DELTA(3*(j-1)+2)=DELTA(3*(j-1)+2)+ay
c   Then the x coordinate:
          IF(DELTA(3*(j-1)+1) .GT. axhalf)
     +       DELTA(3*(j-1)+1)=DELTA(3*(j-1)+1)-ax
          IF(DELTA(3*(j-1)+1) .LT.-axhalf)
     +       DELTA(3*(j-1)+1)=DELTA(3*(j-1)+1)+ax
c   Finally the z coordinate:
          IF(DELTA(3*(j-1)+3) .GT. azhalf)
     +       DELTA(3*(j-1)+3)=DELTA(3*(j-1)+3)-az
          IF(DELTA(3*(j-1)+3) .LT.-azhalf)
     +       DELTA(3*(j-1)+3)=DELTA(3*(j-1)+3)+az
320   CONTINUE
c
      do 371 J=1,NJATMS
        R2st(j)=
     +     DELTA(3*(J-1)+1)**2+DELTA(3*(J-1)+2)**2+DELTA(3*(J-1)+3)**2
371   continue
c
      nintwi=0
c                  nintwi counts the number of atoms interacting with atom i.

      do 372 j=1,njatms
        if(r2st(j) .gt. skinr2) go to 372
c
c   The maximum of R2 is RSKIN2 here to get effectively a neighbor list with
c   buffer region (rskin2 should be larger than rcut2)
c
c        add atom j to the neighbor list of atom i:
          nintwi=nintwi+1
          indpra(1,nintpr+ishift+nintwi)=i
          indpra(2,nintpr+ishift+nintwi)=j
372     continue

        do 3722 k=nintpr+ishift+1,nintpr+ishift+nintwi
          jd=indpra(2,k)
          r2pr(k)=r2st(jd)
          delpr(3*k-2)=delta(3*jd-2)
          delpr(3*k-1)=delta(3*jd-1)
          delpr(3*k)  =delta(3*jd)
3722    continue
        IF(ISAME.EQ.1) then
          do 3723 k=nintpr+ishift+1,nintpr+ishift+nintwi
            indpra(2,k)=indpra(2,k)+jstrt-1
3723      continue
        ENDIF

        if(nintwi .gt. maxneb) maxneb=nintwi
        iptpr1(i,ipot)=nintpr+1
        nintpr=nintpr+nintwi
        if(nintpr .gt. maxnpr) maxnpr=nintpr
        indx=indx+3
500   continue

      nintp(ipot) = nintpr
      iptpr1(natom1+1,ipot)=nintpr+1
      nintst(ipot)=nintpr

      if(nintpr .gt. maxprs) then
        write(lunout,221) nintpr,maxprs
221     format(/'   ERROR in gagafe:  nintpr > maxprs, = ',2i9)
        stop
      endif
      iptcom(ipot+1)=ishift+nintpr+1

      go to 601
c           ------------------------------------------------
600   continue
c
c   Since indupd was not 1, the list is not updated but only the distances:
c
      nintpr=nintst(ipot)
      ishift=iptcom(ipot)-1
      nint3=3*nintpr

      do 385 i=1,natom1-ifisam
        do 380 ipr=iptpr1(i,ipot)+ishift,iptpr1(i+1,ipot)-1+ishift
          delpr(3*ipr-2)=ra1(3*i-2)-ra2(3*indpra(2,ipr)-2)
          delpr(3*ipr-1)=ra1(3*i-1)-ra2(3*indpra(2,ipr)-1)
          delpr(3*ipr)  =ra1(3*i)  -ra2(3*indpra(2,ipr))
380     continue
385   continue
c
c   Periodic Boundary Conditions:
c
      do 381 ipr=1+ishift,nintpr+ishift
c   the x coordinates:
          IF(delpr(3*ipr-2) .GT. axhalf)
     +       delpr(3*ipr-2)=delpr(3*ipr-2)-ax
          IF(delpr(3*ipr-2) .LT.-axhalf)
     +       delpr(3*ipr-2)=delpr(3*ipr-2)+ax
c   the y coordinates:
          IF(delpr(3*ipr-1) .GT. ayhalf)
     +       delpr(3*ipr-1)=delpr(3*ipr-1)-ay
          IF(delpr(3*ipr-1) .LT.-ayhalf)
     +       delpr(3*ipr-1)=delpr(3*ipr-1)+ay
c   the z coordinates:
          IF(delpr(3*ipr) .GT. azhalf)
     +       delpr(3*ipr)=delpr(3*ipr)-az
          IF(delpr(3*ipr) .LT.-azhalf)
     +       delpr(3*ipr)=delpr(3*ipr)+az
381   continue
c
c   Find R2:
c
      do 382 ir=1+ishift,nintpr+ishift
        r2pr(ir)=delpr(3*ir-2)**2+delpr(3*ir-1)**2+delpr(3*ir)**2
382   continue
c
c           ------------------------------------------------

c     Potential specific part:
c
601   continue
c
c  Loop over all interacting pairs:
c
c   Ensure that the potential is effectively zero beond the cutoff:
c
c   Interpret the potential parameters:

      gtrans = potpar(1,ipot)
      preexpA = potpar(2,ipot)
      alphaA = potpar(3,ipot)
      preexpB = potpar(4,ipot)
      alphaB = potpar(5,ipot)
      scale = potpar(6,ipot)
      betaA = potpar(7,ipot)
      betaB = potpar(8,ipot)
      etap=potpar(9,ipot)
      netap=etap+0.01
      gammap=potpar(10,ipot)

      scale1=scale
      scale2=scale1
      betaA1=betaA
      betaA2=betaA1
      betaB1=betaB
      betaB2=betaB1
      eta1=etap
      eta2=eta1
      neta1=eta1+0.01
      neta2=neta1
      gamma1=gammap
      gamma2=gamma1

      sccut=rcut(ipot)
      cut1=sccut
      cut2=sccut

c   -------------------------------------------------------------------------
c  Evaluate pair potential and atomic electron density at the cutoff distance:
      rhocut1=1.0
      if(neta1 .ne. 0) rhocut1=cut1**neta1
      rhocut1 = rhocut1*(exp(-betaA1*cut1)+gamma1*exp(-betaB1*cut1))
      rhocut2=1.0
      if(neta2 .ne. 0) rhocut2=cut2**neta2
      rhocut2 = rhocut2*(exp(-betaA2*cut2)+gamma2*exp(-betaB2*cut2))
      phicut =preexpA*exp(-alphaA*sccut)+preexpB*exp(-alphaB*sccut)
      phicutst(ipot)=phicut
      rhocut1st(ipot)=rhocut1*scale1
      rhocut2st(ipot)=rhocut2*scale2

c   -------------------------------------------------------------------------

      UTOT=0.0
      virialold=0.0
      phitot = 0.0

      do 370 j=1+ishift,nintpr+ishift
        rpr = sqrt(r2pr(j))
        rinv = 1.0/rpr

c       adding a bunch of constants to speed up convergence

         eaA=exp(-alphaA*rpr)
         eaB=exp(-alphaB*rpr)
         ebA=exp(-betaA*rpr)
         ebB=exp(-betaB*rpr)
         fact1=ebA+gamma1*ebB
         fact2=betaA1*ebA+gamma1*betaB1*ebB

c       evaluate atomic electron density at atom i due to j and at j due to i:

        reta1=1.0
        if(neta1.ne.0) reta1=rpr**neta1
        rhoji(j)=(reta1*fact1-rhocut1)*scale1
        rhoij(j)=rhoji(j)

c      evaluate the derivative/r of the atomic electron density:

        if(neta1.ne.0) then
          rhoprji=eta1*rpr**(neta1-2)*fact1-rpr**(neta1-1)*fact2
        else
          rhoprji=-fact2*rinv
        endif
        rhoprji=rhoprji*scale1
        drhojir(j)=rhoprji*r2pr(j)
        rhoprij=rhoprji
        drhoijr(j)=drhojir(j)

c          The arrays drhoijr and drhojir now contain
c          the derivative of the atomic electron
c          density, drho(r), times the distance, r.

c                ---------------------------------------------
c      evaluate pair potential between atom i and j:

        phi(j)=preexpA*eaA+preexpB*eaB-phicut-2.*gtrans*rhoij(j)
        phitot=phitot+phi(j)

cPres:  the derivative of the pair potential divided by r:

        phipror=-(preexpA*alphaA*eaA+preexpB*alphaB*eaB)*rinv
     +          -2.0*gtrans*rhoprij
          phivirst(j)=phipror*r2pr(j)
          virialold=virialold+phivirst(j)

c         At this point the array delpr contains the x,y,z coordinates of the
c           separation vector between each interacting pair.

          j3=3*j
          rhofij(j3-2)=delpr(j3-2)*(-rhoprij)
          rhofij(j3-1)=delpr(j3-1)*(-rhoprij)
          rhofij(j3)  =delpr(j3)  *(-rhoprij)
          rhofji(j3-2)=rhofij(j3-2)
          rhofji(j3-1)=rhofij(j3-1)
          rhofji(j3)  =rhofij(j3)

c      Now use delpr temporarily to store the force components due to the
c            pair potential:
          delpr(j3-2)=delpr(j3-2)*(-phipror)
          delpr(j3-1)=delpr(j3-1)*(-phipror)
          delpr(j3)  =delpr(j3)  *(-phipror)

370       CONTINUE

c        --------------------------------------------------------

c  Fixes related to the cutoff distance:

cFixJune91:
c    Fix Force and Virial due to pairs beond the cutoff distance:

      sccutm=rcut2(ipot)*0.9999
      sccutm1=sccutm
      sccutm2=sccutm

      nprsbeond=0
      nprsbeond1=0
      nprsbeond2=0

      do 377 j=1+ishift,nintpr+ishift
        if(r2pr(j) .lt. sccutm) go to 377
           nprsbeond = nprsbeond + 1
           delpr(3*j-2)=0.0
           delpr(3*j-1)=0.0
           delpr(3*j)  =0.0
           phi(j)= 0.0
           phivirst(j)= 0.0
377   continue

      do 388 j=1+ishift,nintpr+ishift
        if(r2pr(j) .lt. sccutm2) go to 388
           nal = nal + 1
           rhofij(3*j-2)=0.0
           rhofij(3*j-1)=0.0
           rhofij(3*j)  =0.0
           rhoij(j) = 0.0
           drhoijr(j) = 0.0
388   continue

      do 399 j=1+ishift,nintpr+ishift
         if(r2pr(j) .lt. sccutm1) go to 399
           nni = nni + 1
           rhofji(3*j-2) = 0.0
           rhofji(3*j-1)=0.0
           rhofji(3*j)=0.0
           rhoji(j) = 0.0
           drhojir(j) = 0.0
399   continue
c                                          endfixJune91
c        --------------------------------------------------------

      totphi= 0.0
      totvirial= 0.0
      do j=1+ishift,nintpr+ishift
         totphi= totphi + phi(j)
         totvirial= totvirial + phivirst(j)
         Potpr(j)= phi(j)
      end do

      UTOT=totphi
      virial=totvirial

c     --------------------------------------------------------------------

c    Store the x,y and z components of the total force due to the pair
c       potential in the arrays fa1 and fa2:

      do 376 i=1,natom1-ifisam
        do 373 ipr=iptpr1(i,ipot)+ishift,iptpr1(i+1,ipot)-1+ishift
           potpat(i+iatshift1)=potpat(i+iatshift1)+potpr(ipr)
           fa1(3*i-2)=fa1(3*i-2)+delpr(3*ipr-2)
           fa1(3*i-1)=fa1(3*i-1)+delpr(3*ipr-1)
           fa1(3*i)  =fa1(3*i)  +delpr(3*ipr)
373     continue
376   continue

      do 375 j=1+ishift,nintpr+ishift
        potpat(indpra(2,j)+iatshift2)=potpat(indpra(2,j)+iatshift2)+
     +  potpr(j)
        fa2(3*indpra(2,j)-2)=fa2(3*indpra(2,j)-2)-delpr(3*j-2)
        fa2(3*indpra(2,j)-1)=fa2(3*indpra(2,j)-1)-delpr(3*j-1)
        fa2(3*indpra(2,j))  =fa2(3*indpra(2,j))  -delpr(3*j)
375   continue

      RETURN
      END
