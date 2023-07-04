c  the case of exactly two components is hardwired into this routine(!!!)
C  This subroutine calculate the total embedding energy of the system and the
c  contribution of the embedding energy to the virial.


      SUBROUTINE EMBED(FRHOTOT,embvir)

      implicit real*8 (a-h,o-z)

      include 'commonblks/parameters.cmn'
      include 'commonblks/comconf.cmn'
      include 'commonblks/comluns.cmn'
      include 'commonblks/comintlis.cmn'

      common/things/ isame(3),ifisame(3),nshift(3)
      common /unscale/RALOCAL(MAXCOO)

      dimension rho1(MAXATOMS,3),rho2(MAXATOMS,3),fprho(MAXATOMS,3),
     +          emb(MAXATOMS,3),embvirst(MAXATOMS,3),
     +          drho1r(MAXATOMS,3),drho2r(MAXATOMS,3)

      do ipot = 1,3
        nshift(ipot)=iptcom(ipot)-1
      end do

      isame(1) = 1
      isame(2) = 1
      isame(3) = 0

      ifisame(1) = 1
      ifisame(2) = 1
      ifisame(3) = 0

c   Initialize the arrays that will have the total density of a
c   particular atom k.
      do 10 ipot = 1,3
         do k = 1,NATOMS
            rho1(k,ipot) = 0.0
            rho2(k,ipot) = 0.0
            drho1r(k,ipot) = 0.0
            drho2r(k,ipot) = 0.0
         end do
 10   continue

cw            write(lunout,*) '  Density contributions:  Just bef do 30'

        do 30 ipot = 1,3
c                loop over the three types of pairs, A-A, B-B, A-B.
c   Add up density contributions at atom i due to atom j of pair number
        do 40 i = 1,natm1(ipot) - ifisame(ipot)
cw            write(lunout,*) '  i=',i
          do ipr = iptpr1(i,ipot)+nshift(ipot),iptpr1(i+1,ipot)-1+
     +        nshift(ipot)
              rho1(i,ipot) = rho1(i,ipot) + rhoij(ipr)
              drho1r(i,ipot) = drho1r(i,ipot) + drhoijr(ipr)

cw              write(lunout,*) '   ipr = ',ipr,'  rhoij=',rhoij(ipr),
cw     +                        ' rhoji=',rhoji(ipr)

            end do
40          continue

          do 50 j = 1+nshift(ipot),nintp(ipot)+nshift(ipot)
            rho2(indpra(2,j),ipot) = rho2(indpra(2,j),ipot)+rhoji(j)
            drho2r(indpra(2,j),ipot) = drho2r(indpra(2,j),ipot)+
     +         drhojir(j)
cw            write(lunout,*) '   j = ',j,'  indpra(2,j)',indpra(2,j),
cw     +                      '  rhoji=',rhoji(j)
c           write(*,*)'rho2',i,ipot,rho2(i,ipot)
50          continue

30        continue

      do 60 ipot = 1,2
           do j = 1,natm1(ipot)
             rho1(j,ipot) = rho1(j,ipot) + rho2(j,ipot)
             drho1r(j,ipot) = drho1r(j,ipot) + drho2r(j,ipot)
         end do
60         continue

      do 70 ipot = 1,2
          do j = 1,natm2(ipot)
            rho2(j,ipot) = rho1(j,ipot)
            drho2r(j,ipot) = drho1r(j,ipot)
          end do
70        continue

      do j = 1,natm1(1)
          rho1(j,1) = rho1(j,1) + rho1(j,3)
          drho1r(j,1) = drho1r(j,1) + drho1r(j,3)
        end do

      do k = 1,natm1(2)
          rho1(k,2) = rho1(k,2) + rho2(k,3)
          drho1r(k,2) = drho1r(k,2) + drho2r(k,3)
      end do

c      At this point rho1(j,1) should contain the total density at the
c      location of atom j of component 1 and
c      rho1(j,2) should contain the total density at the
c      location of atom j of component 2.
c      write(*,*)'rhosum',rho1(1,1),rho1(2,1)
cw         write(lunout,*)
cw         write(lunout,*)' EMBED: Tot electron density at atom locations:'
cw          iatsh=0
cw          do 71 ic=1,2
cw            write(lunout,*) ' component ',ic,': atom:    rho1:     drho1r:'
cw            if(ipot .eq. 2) iatsh=Natms(1)
cw	    do k=1,natm1(ic)
cw	      write(lunout,*) '            ',k+iatsh,rho1(k,ic),drho1r(k,ic)
cw            enddo
cw71          continue

c         -------------------------------------------------------

c     Now find the force components due to the embedding function:
      do 80 ipot = 1,2
          do j = 1,NATOMS
            fprho(j,ipot) = 0.0
            emb(j,ipot) = 0.0
          end do
80        continue


      do 90 ipot = 1,2
        call dfrhodrho(rho1,natm1,ipot,fprho)
          call fint(rho1,natm1,ipot,emb)
c	  write(6,*)ipot,emb(1,ipot),fprho(1,ipot)
90        continue

        do 106 ipot = 1,2
          do i = 1,natm1(ipot)
            embvirst(i,ipot)=fprho(i,ipot)*drho1r(i,ipot)
          end do
 106    continue


cFPI:   scale embedding energy and virial contr. due to images in FPI chains:
        if(nFPI .gt. 0) then
          fact=1./float(nimFPI(1))
c             NB:  assume here that all FPI have the same number of images.
          iatsh=0
          do 66 ipot=1,2
            if(ipot .eq. 2) iatsh=Natms(1)
          do i= 1,natm1(ipot)
            if(itag(i+iatsh) .ge. 100 .and. itag(i+iatsh) .lt. 400) then

cw                  write(lunout,*) '  EMBED:  FPI corr.  i=',i,'  iatsh=',iatsh

                  emb(i,ipot)=emb(i,ipot)*fact
                  embvirst(i,ipot)=embvirst(i,ipot)*fact
              endif
          enddo
66          continue
        endif


c     Calculate the total embedding energy and contribution to virial:
      frhotot = 0.0
        do 100 ipot = 1,2
          do i = 1,natm1(ipot)
            frhotot = frhotot + emb(i,ipot)

cw	    write(lunout,*) '  i=',i,'  ipot=',ipot,'   emb = ',emb(i,ipot)

          end do
 100    continue

      embvir=0.0
        do 105 ipot = 1,2
          do i = 1,natm1(ipot)
            embvir= embvir + embvirst(i,ipot)
          end do
 105    continue
c             ----------------------------------------------

cw        write(lunout,*)
cw        write(lunout,*) '  The rhofij and rhofji vectors before do 110:'
cw	  write(lunout,*)
cw        do 160 ipot=1,3
cw	  write(lunout,*) ' ipot = ',ipot
cw	  write(lunout,*)'  ipair:   i:    j:   rhofij:                rhofji:'
cw	  do ip=1+nshift(ipot),nintp(ipot)+nshift(ipot)
cw	    write(lunout,161) ip,indpra(1,ip),indpra(2,ip),
cw     +                        rhofij(3*ip-2),rhofij(3*ip-1),rhofij(3*ip),
cw     +                        rhofji(3*ip-2),rhofji(3*ip-1),rhofji(3*ip)
cw161         format('   ',i4,i5,i5,3g12.4,'  ',3g12.4)
cw          enddo
cw160       continue


        do 110 ipot = 1,2
          do j = 1 + nshift(ipot),nintp(ipot)+nshift(ipot)
           rhofij(3*j-2) = rhofij(3*j-2)*fprho(indpra(1,j),ipot)
           rhofij(3*j-1) = rhofij(3*j-1)*fprho(indpra(1,j),ipot)
           rhofij(3*j) = rhofij(3*j)*fprho(indpra(1,j),ipot)
          end do
 110    continue

      do 120 ipot = 1,2
          do j = 1 + nshift(ipot),nintp(ipot) + nshift(ipot)
           rhofji(3*j-2) = rhofji(3*j-2)*fprho(indpra(2,j),ipot)
           rhofji(3*j-1) = rhofji(3*j-1)*fprho(indpra(2,j),ipot)
           rhofji(3*j) = rhofji(3*j)*fprho(indpra(2,j),ipot)
          end do
 120  continue

      do 130 j = 1 + nshift(3),nintp(3) +nshift(3)
           rhofij(3*j-2) = rhofij(3*j-2)*fprho(indpra(1,j),1)
           rhofij(3*j-1) = rhofij(3*j-1)*fprho(indpra(1,j),1)
           rhofij(3*j) = rhofij(3*j)*fprho(indpra(1,j),1)
 130  continue

      do 140 k = 1 + nshift(3),nintp(3) + nshift(3)
           rhofji(3*k-2) = rhofji(3*k-2)*fprho(indpra(2,k),2)
           rhofji(3*k-1) = rhofji(3*k-1)*fprho(indpra(2,k),2)
           rhofji(3*k) = rhofji(3*k)*fprho(indpra(2,k),2)
 140   continue

        ntot = nintp(1) + nintp(2) + nintp(3)

      do j = 1,ntot
           rhofij(3*j-2) = rhofij(3*j-2) + rhofji(3*j-2)
           rhofij(3*j-1) = rhofij(3*j-1) + rhofji(3*j-1)
           rhofij(3*j) = rhofij(3*j) + rhofji(3*j)
c	   write(*,*)rhofij(3*j-2),rhofij(3*j-1),rhofij(3*j)
      end do

c     At this point the array  rhofij(i) should contain
c     the force components due to embedding energy,
c          [F'k*rho'j + F'j*rho'k] \vec r_jk / r_jk
c     write out the info:
cw        write(lunout,*)
cw        write(lunout,*) '  The rhofij and rhofji vectors after do 140:'
cw	  write(lunout,*)
cw        do 150 ipot=1,3
cw	  write(lunout,*) ' ipot = ',ipot
cw	  write(lunout,*)'  ipair:   i:    j:   rhofij:                rhofji:'
cw	  do ip=1+nshift(ipot),nintp(ipot)+nshift(ipot)
cw	    write(lunout,151) ip,indpra(1,ip),indpra(2,ip),
cw     +                        rhofij(3*ip-2),rhofij(3*ip-1),rhofij(3*ip),
cw     +                        rhofji(3*ip-2),rhofji(3*ip-1),rhofji(3*ip)
cw151         format('   ',i4,i5,i5,3g12.4,'  ',3g12.4)
cw          enddo
cw150       continue


       call sumforce(1,FA(1),FA(1))

       istrt = 3*natm1(1) + 1
       jstrt = istrt
       call sumforce(2,FA(istrt),FA(jstrt))

       if (natm1(3).eq.0) go to 200

       mstrt = 1
       nstrt = 3*natm1(3) + 1

       call sumforce(3,FA(mstrt),FA(nstrt))

200    return
       end
