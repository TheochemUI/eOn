c   version e93:   add gregs SPF routine for chain relaxation
c   EAM version b91:   add FPI forces and potential by calling FPIforce
c
C     THIS SUBROUTINE COMPUTES THE FORCES by calling GAGAFE for each type of
c     interaction.

      SUBROUTINE FORCE(N, R, F, U, BX, BY, BZ)

      implicit real * 8 (a - h, o - z)

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

c      dimension R(MAXCOO),F(MAXCOO)
      dimension R(3 * N), F(3 * N), U(1)

      AX = BX
      AY = BY
      AZ = BZ

      Rxyz(1) = AX
      Rxyz(2) = AY
      Rxyz(3) = AZ

      NATYPE = 1
      NATOMS = N
      NATMS(1) = N
      NATMS(2) = 0

      if(initflag.eq.1) then
         initflag = 2
         indupd = 1
      else
         do I = 1, NATOMS
            J = 3 * (I - 1)
            Rij(1) = Rold(J + 1) - R(J + 1)
            Rij(2) = Rold(J + 2) - R(J + 2)
            Rij(3) = Rold(J + 3) - R(J + 3)
c            write(*,*) Rij(1),Rij(2),Rij(3)
c            write(*,*)
            Rij(1) = Rij(1) / Rxyz(1)
            Rij(2) = Rij(2) / Rxyz(2)
            Rij(3) = Rij(3) / Rxyz(3)
            R1=1.0
            Rij(1) = mod(Rij(1) + 1000.5, R1) - 0.5
            Rij(2) = mod(Rij(2) + 1000.5, R1) - 0.5
            Rij(3) = mod(Rij(3) + 1000.5, R1) - 0.5
            Rij(1) = Rij(1) * Rxyz(1)
            Rij(2) = Rij(2) * Rxyz(2)
            Rij(3) = Rij(3) * Rxyz(3)
c            write(*,*) Rij(1),Rij(2),Rij(3)
c            write(*,*)
            dR = 0
            dR = dR + Rij(1) * Rij(1)
            dR = dR + Rij(2) * Rij(2)
            dR = dR + Rij(3) * Rij(3)
            dR = sqrt(dR)
            if(I.eq.1) then
               dR1 = dR
            elseif(I.eq.2) then
               if(dR1.gt.dR) then
                  dR2 = dR
               else
                  dR2 = dR1;
                  dR1 = dR;
               endif
            else
               if(dR1.lt.dR) then
                  dR2 = dR1
                  dR1 = dR
               elseif(dR2.lt.dR) then
                  dR2 = dR
               endif
            endif
         enddo
         if((dR1 + dR2).gt.(rskin(1) - rcut(1))) then
            indupd = 1
         else
            indupd = 0
         endif
      endif

      if(indupd.eq.1) then
c        nncount=nncount+1
c        write(*,*) 'NNUpdate: ',nncount
         do I = 1, 3 * NATOMS
            Rold(I) = R(I)
         enddo
      endif


      DO I = 1, 3 * NATOMS
         FA(I) = 0.
         RALOCAL(I) = R(I)
      enddo

      do i = 1, natoms
         Potenpat(i) = 0.
      enddo

C  INTERACTIONS BETWEEN ATOMS OF TYPE 'ITYPE' WITH THEMSELVES:
      ITPTR = 1
      DO 200 ITYPE = 1,NATYPE
         iatshift1 = (ITPTR-1)/3
         iatshift2 = iatshift1
         CALL GAGAFE(NATMS(ITYPE), RALOCAL(ITPTR), FA(ITPTR),
     X      NATMS(ITYPE), RALOCAL(ITPTR), FA(ITPTR),
     X      utotITYPE, virialITYPE, ITYPE,1)
         UT(ITYPE * (ITYPE + 1) / 2) = utotITYPE
         ITPTR = ITPTR + 3 * NATMS(ITYPE)
200   CONTINUE

      UTOT = 0.
      DO 600 I = 1, NATYPE * (NATYPE + 1) / 2
         UTOT = UTOT + UT(I)
c         write(*,*) UTOT
600   continue

      totpair = utot
c      write(*,*) totpair
      CALL EMBED(FRHOTOT, embvir)
c      write(*,*) FRHOTOT, UTOT
      UTOT = UTOT + FRHOTOT
      do i = 1, 3 * NATOMS
c      	 write (*,*) FA(i)
         F(i) = FA(i)
      end do
c      write(*,*) UTOT
      U(1)=UTOT
c      write(*,*) U


      RETURN
      END
