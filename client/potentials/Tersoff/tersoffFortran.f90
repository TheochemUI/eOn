!   for version EAMf of the MD code, Feb 1995.  (Modified by Graeme)
!   for EAM version E     Feb 1994
!   Tersoff potential:   Sept 14. 1993   based on EAMd version
!      two components
!      No neighborlist, no FPI chains

!     SUBROUTINE TO COMPUTE FORCE, ENERGY AND VIRIAL FOR A GROUP OF
!     ATOMS INTERACTING WITH ANOTHER GROUP OF ATOMS.
!     THE GROUPS MAY BE THE SAME.
!     The force is scaled in such a way that the real force is ax*fa.
!


!AA      SUBROUTINE tersoff(NATOM,RA,FA,UTOT,ax, ay, az, MAXCOO)
      SUBROUTINE tersoff(NATOM,RA,FA,UTOT,ax, ay, az)

      implicit real*8 (a-h,o-z)

      real*8 ax, ay, az
!AA      integer(4) NATOM, MAXCOO
      integer(4) NATOM

!AA      DIMENSION RA(MAXCOO),FA(MAXCOO)
      DIMENSION RA(3*NATOM),FA(3*NATOM)

      real*8 lamdaij,lamda,muij,kaij,nei,np,mu
      DIMENSION PL(3),drij(3),drik(3),XYZ(3,NATOM),FRC(3,NATOM),    &
     &          ftj(3,NATOM),ftk(3,NATOM),lid(NATOM)

!   ---------------------------------------------------------------------
!   Function definitions:

      Frfn(x)=Aij*exp(-lamdaij*x)
      Fafn(x)=-Bij*exp(-muij*x)
      Fcfn(x)=0.5+0.5*cos(pi*(x-RLD)/(Sij-RLD))
      Bfn(z)=kaij/(1.+(beta*z)**nei)**(0.5/nei)
      Gfn(t)=1.0+(ci/di)**2-ci**2/(di**2+(hi-t)**2)

      DFrfn(x)=-Aij*lamdaij*exp(-lamdaij*x)
      DFafn(x)=Bij*muij*exp(-muij*x)
      DFcfn(x)=-pi/(2.*(Sij-RLD))*sin(pi*(x-RLD)/(Sij-RLD))
      XSij(z)=-kaij*0.5*beta**nei*z**(nei-1)/(1.            &
     &     +(beta*z)**nei)**((0.5+nei)/nei)
      Qfn(t)=-2.*ci**2*(hi-t)/(di**2+(hi-t)**2)**2

      pi=4.0*atan(1.0)

! --------------------------------------------------------------
!  Potential parameters:  (for silicon)

      A=1.8308e+3
      B=4.7118e+2
      lamda=2.4799
      mu=1.7322
      bt=1.1000e-6
      np=7.8734e-1
      c=1.0039e+5
      d=1.6217e+1
      h=-5.9825e-1
      R=2.7
      S=3.0
      chipure=1.0

      PL(1)=ax
      PL(2)=ay
      PL(3)=az

      k=0
      do i=1,natom
       do j=1,3
        k=k+1
        XYZ(j,i)=RA(k)
       enddo
      enddo
      k=0

      do 99 i=1,natom
         do 99 l=1,3
            FRC(l,i)=0.0
 99   continue

      ENG=0.0

      do 100 i=1,natom
         ci=c
         di=d
         hi=h
         beta=bt
         nei=np
         do 101 j=1,natom
            if(j.eq.i) goto 101
            Aij=a
            Bij=b
            RLD=r
            Sij=s
            lamdaij=lamda
            muij=mu
            kaij=chipure

            rij2=0.0
            do 120 ii=1,3
               drij(ii)=XYZ(ii,i)-XYZ(ii,j)
               if(abs(drij(ii)).ge.Sij) then
                  drij(ii)=drij(ii)*(1.-PL(ii)/abs(drij(ii)))
                  if(abs(drij(ii)).ge.Sij) goto 101
               endif
               rij2=rij2+drij(ii)**2
 120        continue

            rij=sqrt(rij2)
            if(rij.ge.Sij) goto 101

            if(rij.le.RLD) then
               fcij=1.0
               dfcij=0.0
            else
               fcij=Fcfn(rij)
               dfcij=DFcfn(rij)
            endif

            zeta=0.0
            do 102 k=1,natom
               lid(k)=0
               if(k.eq.i.or.k.eq.j) goto 102
               RLD=r
               Sij=s
               rik2=0.0
               do 130 ii=1,3
                  drik(ii)=xyz(ii,i)-xyz(ii,k)
                  if(abs(drik(ii)).ge.Sij) then
                     drik(ii)=drik(ii)*(1.0-PL(ii)/abs(drik(ii)))
                     if(abs(drik(ii)).ge.Sij) goto 102
                  endif
                  rik2=rik2+drik(ii)**2
 130           continue

               rik=sqrt(rik2)
               if(rik.ge.Sij) goto 102
               lid(k)=1
               if(rik.le.RLD) then
                  fcik=1.0
                  dfcik=0.0
               else
                  fcik=Fcfn(rik)
                  dfcik=DFcfn(rik)
               endif

               tijk=0.0
               do 140 ii=1,3
                  tijk=tijk+drij(ii)*drik(ii)
 140           continue

               cosijk=tijk/(rij*rik)
               drijk=rij-rik
               zeta=zeta+fcik*Gfn(cosijk)

               cf1=fcik*Qfn(cosijk)

               do 104 l=1,3

                  ftj(l,k)=0.0
                  ftk(l,k)=0.0

                  ftk(l,k)=ftk(l,k)                  &
     &                 +drik(l)/rik*dfcik*Gfn(cosijk)

                  ftj(l,k)=ftj(l,k)+cf1*(drik(l)/rik       &
     &                 -cosijk*drij(l)/rij)/rij
                  ftk(l,k)=ftk(l,k)+cf1*(drij(l)/rij        &
     &                 -cosijk*drik(l)/rik)/rik

 104           continue
 102        continue

!         write(6,*) ' '
!         write(6,*) '  From Gagafe:  zeta = ',zeta
            if(zeta.lt.1.e-15) then
               cf=0.0
               cf1=1.0
            else
               cf=fcij*Fafn(rij)*XSij(zeta)
               Cf1=Bfn(zeta)
            endif
            Fcf=dfcij*(Frfn(rij)+Cf1*Fafn(rij))                 &
     &           +fcij*(DFrfn(rij)+Cf1*DFafn(rij))

            do 105 l=1,3
               fik=0.0
               fjk=0.0
               do 103 k=1,natom
                  if (lid(k).eq.0) goto 103

                  fik=fik+ftj(l,k)+ftk(l,k)
                  fjk=fjk+ftj(l,k)

                  FRC(l,k)=FRC(l,k)-cf*ftk(l,k)
 103           continue

               Cf2=drij(l)/rij

               FRC(l,i)=FRC(l,i)+cf2*Fcf+cf*fik
               FRC(l,j)=FRC(l,j)-cf2*Fcf-cf*fjk

 105        continue

            ENG=ENG+fcij*(Frfn(rij)+Cf1*Fafn(rij))

 101        continue

 100  continue

      do 200 i=1,natom
         do 200 l=1,3
            FRC(l,i)=FRC(l,i)/2.0
 200  continue

      ENG=ENG/2.0

!   Translate back:   Change sign to get force rather than gradient

      utot=eng
      do i=1,natom
         FA(3*(i-1)+1)=-FRC(1,i)
         FA(3*(i-1)+2)=-FRC(2,i)
         FA(3*(i-1)+3)=-FRC(3,i)
      enddo

      return
      end


!-----------------------------------------------------------------------------------!
!
!  SUBROUTINE potinit()
!     IMPLICIT NONE
!
!  END SUBROUTINE potinit
!
!-----------------------------------------------------------------------------------!
