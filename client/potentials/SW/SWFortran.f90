!     Last change:  AA    3 Apr 2002    2:31 pm
!-----------------------------------------------------------------------------------!

  SUBROUTINE sw(NATOMS,pos,F,energy,bx,by,bz)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: NATOMS
    REAL(8),intent(in),dimension(3*NATOMS) :: pos
    REAL(8),INTENT(IN) :: bx,by,bz
    REAL(8),intent(out) :: energy
    REAL(8),intent(out),dimension(3*NATOMS) :: F

    REAL(8),dimension(NATOMS) :: fx,fy,fz
    REAL(8),dimension(NATOMS) :: x,y,z
    integer :: i,j,i_id,j_id,ind_j,k,k_id,ind_k
    REAL(8) :: xi, yi, zi, xij, yij, zij, rij, rij2
    REAL(8) :: twobodyenergy, threebodyenergy
    REAL(8) :: cos_x_ij, cos_y_ij, cos_z_ij, invrij, rhoij
    REAL(8) :: one_o_a_ij, expo, gam_o_a_ij, exp_gam_ij, r_to_minusp
    REAL(8) :: one_o_a2, term1, term2,fact3_ij
    REAL(8) :: ffx,ffy,ffz, term_ij, term_ik
    REAL(8) :: dcosdxj,dcosdyj,dcosdzj,dcosdxk,dcosdyk,dcosdzk
    REAL(8) :: dhdcos_ij,dhdcos_ik,cos_p_1o3,cos_jik
    REAL(8) :: cos_x_ik, cos_y_ik, cos_z_ik, one_o_a_ik, gam_o_a_ik, fact
    REAL(8) :: xik,yik,zik,rik, rik2, invrik,rhoik,xxx

    INTEGER,PARAMETER :: MAXNATOMS=20000,                                      &
    &                    MAXNEI=20,                                            &
    &                    P=4
    REAL(8),PARAMETER :: SIGMA = 2.0951D0,                                     &
    &                    ALPHA = 1.8d0,                                        &
    &                    RCUT = (ALPHA * SIGMA),                               &
    &                    RSKIN = 4.0d0,                                        &
    &                    invsig=1.0d0/SIGMA,                                   &
    &                    EPSILON = 2.16826D0 * 1.06767385d0,                   & !* 1.07d0
    &                    LAMBDA = 21.0d0,                                      & !* 1.07d0
    &                    A =  7.049556277D0,                                   &
    &                    BETA = 0.6022245584D0,                                &
    &                    GAMMA = 1.2D0,                                        &
    &                    ONE_THIRD = (1.0d0/ 3.0d0),                           &
    &                    A_EPS = (A * EPSILON),                                &
    &                    PI = 3.14159265358979D0

    LOGICAL,SAVE :: initflag=.true. , indupd=.true.
    REAL(8),DIMENSION(MAXNATOMS),SAVE :: xold,yold,zold
    INTEGER,DIMENSION(MAXNATOMS),SAVE :: numnei
    INTEGER,DIMENSION(MAXNATOMS,MAXNEI),SAVE :: nei
    REAL(8),DIMENSION(2,2),SAVE :: S0,A0,rcut2

    INTERFACE
      SUBROUTINE check_neigbours(NATOMS,x,y,z,bx,by,bz,MAXNATOMS,xold,yold,zold,&
                                 initflag,indupd,RCUT,RSKIN)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NATOMS,MAXNATOMS
        REAL(8),INTENT(IN),DIMENSION(NATOMS) :: x,y,z
        REAL(8),INTENT(INOUT),DIMENSION(NATOMS) :: xold,yold,zold
        REAL(8),INTENT(IN) :: bx,by,bz,RCUT,RSKIN
        LOGICAL,INTENT(INOUT) :: initflag,indupd
      END SUBROUTINE check_neigbours

      !-----*****-----!

      SUBROUTINE neighbours(NATOMS,x,y,z,bx,by,bz,MAXNATOMS,MAXNEI,numnei,nei,RSKIN)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NATOMS,MAXNATOMS,MAXNEI
        REAL(8),intent(in),dimension(NATOMS) :: x,y,z
        REAL(8),INTENT(IN) :: bx,by,bz,RSKIN
        INTEGER,INTENT(INOUT),DIMENSION(MAXNATOMS) :: numnei
        INTEGER,INTENT(INOUT),DIMENSION(MAXNATOMS,MAXNEI) :: nei
      END SUBROUTINE
    END INTERFACE

    IF (initflag) THEN
      S0=0.0d0
      A0=0.0d0
      rcut2=0.0d0
      S0(1,1) = 3.6d0
      A0(1,1) = 0.0d0
      rcut2(1,1) = RCUT*RCUT
      xold=0.0d0
      yold=0.0d0
      zold=0.0d0
    END IF

    j=0
    DO i=1,3*NATOMS,3
      j=j+1
      x(j)=pos(i)
      y(j)=pos(i+1)
      z(j)=pos(i+2)
      fx(j)=0.0d0
      fy(j)=0.0d0
      fz(j)=0.0d0
    END DO

  !Generate a neighbour list
    CALL check_neigbours(NATOMS,x,y,z,bx,by,bz,MAXNATOMS,xold,yold,zold,initflag,&
                         indupd,RCUT,RSKIN)
    IF (indupd) CALL neighbours(NATOMS,x,y,z,bx,by,bz,MAXNATOMS,MAXNEI,numnei,nei,&
                                RSKIN)
!    CALL neighbours(NATOMS,x,y,z,bx,by,bz)

  ! We first set-up pointers for the x, y, z components in the force

    twobodyenergy = 0.0d0
    threebodyenergy = 0.0d0

    DO i=1, NATOMS
!      i_id = type(i)
      i_id=1
      xi = x(i)
      yi = y(i)
      zi = z(i)
      DO ind_j=1, numnei(i)
        j = nei(i,ind_j)
!        j_id = type(j)
        j_id=1

      ! Pair interactions of i and j
      ! Distance, with periodic boundary conditions
        xij=x(j)-xi
        yij=y(j)-yi
        zij=z(j)-zi
        xij=xij-bx*ANINT(xij/bx)
        yij=yij-by*ANINT(yij/by)
        zij=zij-bz*ANINT(zij/bz)

        rij2 = xij*xij + yij*yij + zij*zij

      ! Check the cut-off before proceeding
        IF( rij2 < rcut2(i_id,j_id) ) THEN
          rij = sqrt(rij2)
          invrij = 1.0d0 / rij
          rhoij = rij * invsig
          cos_x_ij = xij * invrij
          cos_y_ij = yij * invrij
          cos_z_ij = zij * invrij

        ! Some useful quantities
          one_o_a_ij =1.0/(rhoij-ALPHA)
          IF (one_o_a_ij > -300.0d0) THEN
            expo=exp(one_o_a_ij)
          ELSE
            expo=0.0d0
          END IF
          gam_o_a_ij=GAMMA*one_o_a_ij

          IF (gam_o_a_ij > -300.0d0) THEN
            exp_gam_ij=exp(gam_o_a_ij)
          ELSE
            exp_gam_ij=0.0d0
          END IF
          r_to_minusp=rhoij ** (-1*P)

        ! Two body energy and force
          term1=A_EPS*(BETA*r_to_minusp-1.0)*expo
          one_o_a2 = one_o_a_ij * one_o_a_ij
          term2=(one_o_a2*term1+A_EPS*P*BETA*r_to_minusp*expo/rhoij)*invsig

        ! Contribution to the binary repulsive term
          IF(rij <=  S0(i_id,j_id)) THEN
            term1=term1+A0(i_id,j_id)*(cos(PI*rij/S0(i_id,j_id))+1.0d0)
            term2=term2+A0(i_id,j_id)*PI/S0(i_id,j_id)* sin(PI*rij/S0(i_id,j_id))
          ENDIF
          twobodyenergy=twobodyenergy + 0.5d0*term1;

          fx(i)=fx(i)-term2*cos_x_ij;
          fy(i)=fy(i)-term2*cos_y_ij;
          fz(i)=fz(i)-term2*cos_z_ij;

        ! Prepare for the three body term
          fact3_ij=gam_o_a_ij*one_o_a_ij*invsig

          DO ind_k = ind_j+1, numnei(i)
          ! Triplet interaction with i in the middle; all interactions
            k = nei(i,ind_k)
!            k_id = type(k)
            k_id=1

          ! Distance, with periodic boundary conditions
            xik=x(k)-xi
            yik=y(k)-yi
            zik=z(k)-zi
            xik=xik-bx*ANINT(xik/bx)
            yik=yik-by*ANINT(yik/by)
            zik=zik-bz*ANINT(zik/bz)

            rik2 = xik*xik + yik*yik + zik*zik

          ! Check whether the distance is too large
            IF (rik2<rcut2(i_id,k_id))  THEN
              rik=sqrt(rik2)
              invrik=1.0/rik
              rhoik=rik*invsig
              cos_x_ik=xik*invrik
              cos_y_ik=yik*invrik
              cos_z_ik=zik*invrik

            ! Some useful quantities
              one_o_a_ik=1.0D0/(rhoik-ALPHA)
              gam_o_a_ik=GAMMA*one_o_a_ik

              IF (gam_o_a_ik > -300.0d0) THEN
                xxx=exp(gam_o_a_ik)
              ELSE
                xxx=0.0d0
              END IF
              fact=EPSILON*LAMBDA*xxx*exp_gam_ij

              cos_jik  =cos_x_ij*cos_x_ik+cos_y_ij*cos_y_ik+ cos_z_ij*cos_z_ik
              cos_p_1o3=cos_jik+ONE_THIRD

            ! Energy (added only to central atom)
              threebodyenergy=threebodyenergy+fact*cos_p_1o3*cos_p_1o3

            ! Force
              term_ij=fact*fact3_ij*cos_p_1o3*cos_p_1o3;
              dhdcos_ij=2*fact*cos_p_1o3;
              term_ik=fact*gam_o_a_ik*one_o_a_ik*cos_p_1o3*cos_p_1o3/SIGMA;
              dhdcos_ik=2*fact*cos_p_1o3;

              dcosdxj=(cos_x_ik-cos_jik*cos_x_ij)*invrij;
              dcosdyj=(cos_y_ik-cos_jik*cos_y_ij)*invrij;
              dcosdzj=(cos_z_ik-cos_jik*cos_z_ij)*invrij;

              dcosdxk=(cos_x_ij-cos_jik*cos_x_ik)*invrik;
              dcosdyk=(cos_y_ij-cos_jik*cos_y_ik)*invrik;
              dcosdzk=(cos_z_ij-cos_jik*cos_z_ik)*invrik;

              ffx=term_ij*cos_x_ij-dhdcos_ij*dcosdxj;
              ffy=term_ij*cos_y_ij-dhdcos_ij*dcosdyj;
              ffz=term_ij*cos_z_ij-dhdcos_ij*dcosdzj;
              fx(j)=fx(j)+ffx;
              fy(j)=fy(j)+ffy;
              fz(j)=fz(j)+ffz;
              fx(i)=fx(i)-ffx;
              fy(i)=fy(i)-ffy;
              fz(i)=fz(i)-ffz;
              ffx=term_ik*cos_x_ik-dhdcos_ik*dcosdxk;
              ffy=term_ik*cos_y_ik-dhdcos_ik*dcosdyk;
              ffz=term_ik*cos_z_ik-dhdcos_ik*dcosdzk;
              fx(k)=fx(k)+ffx;
              fy(k)=fy(k)+ffy;
              fz(k)=fz(k)+ffz;
              fx(i)=fx(i)-ffx;
              fy(i)=fy(i)-ffy;
              fz(i)=fz(i)-ffz;
            ENDIF
          END DO
        ENDIF
      END DO
    END DO

    j=0
    DO i=1,3*NATOMS,3
      j=j+1
      F(i)=fx(j)
      F(i+1)=fy(j)
      F(i+2)=fz(j)
    END DO

    energy= twobodyenergy+threebodyenergy

  END SUBROUTINE sw

!-----------------------------------------------------------------------------------!

  SUBROUTINE check_neigbours(NATOMS,x,y,z,bx,by,bz,MAXNATOMS,xold,yold,zold,initflag,&
                             indupd,RCUT,RSKIN)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: NATOMS,MAXNATOMS
    REAL(8),INTENT(IN),DIMENSION(NATOMS) :: x,y,z
    REAL(8),INTENT(INOUT),DIMENSION(NATOMS) :: xold,yold,zold
    REAL(8),INTENT(IN) :: bx,by,bz,RCUT,RSKIN
    LOGICAL,INTENT(INOUT) :: initflag,indupd

    REAL(8),DIMENSION(3) :: Rij
    REAL(8),PARAMETER :: R1=1.0d0
    REAL(8) :: dR0,dR1,dR2
    INTEGER :: i

    IF(initflag) THEN
      initflag=.false.
      indupd=.true.
    ELSE
      DO i=1,NATOMS
        Rij(1)=xold(i)-x(i)
        Rij(2)=yold(i)-y(i)
        Rij(3)=zold(i)-z(i)
        Rij(1)=Rij(1)-bx*ANINT(Rij(1)/bx)
        Rij(2)=Rij(2)-by*ANINT(Rij(2)/by)
        Rij(3)=Rij(3)-bz*ANINT(Rij(3)/bz)
        dR0=SQRT(DOT_PRODUCT(Rij,Rij))
        IF(i == 1) THEN
          dR1=dR0
        ELSE IF(i == 2) THEN
          IF(dR1 > dR0) THEN
            dR2=dR0
          ELSE
            dR2=dR1
            dR1=dR0
          END IF
        ELSE
          IF(dR1 < dR0) THEN
            dR2=dR1
            dR1=dR0
          ELSE IF(dR2 < dR0) THEN
            dR2=dR0
          END IF
        END IF
      END DO
      IF((dR1+dR2) > (RSKIN-RCUT)) THEN
        indupd=.true.
      ELSE
        indupd=.false.
      END IF
    END IF
    IF(indupd) THEN
      xold(1:NATOMS)=x
      yold(1:NATOMS)=y
      zold(1:NATOMS)=z
    END IF
  RETURN
  END SUBROUTINE check_neigbours

!-----------------------------------------------------------------------------------!

! Subroutine neighbours.f90
!
! This subroutine should depend only on the details of the physics and not
! on the ART algorithm per se. It should be called only from the force or
! energy routine
!
! Box should include all rescaling.
! pos is in box units, from -0.5 to 0.5
!
  SUBROUTINE neighbours(NATOMS,x,y,z,bx,by,bz,MAXNATOMS,MAXNEI,numnei,nei,RSKIN)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: NATOMS,MAXNATOMS,MAXNEI
    REAL(8),intent(in),dimension(NATOMS) :: x,y,z
    REAL(8),INTENT(IN) :: bx,by,bz,RSKIN
    INTEGER,INTENT(INOUT),DIMENSION(MAXNATOMS) :: numnei
    INTEGER,INTENT(INOUT),DIMENSION(MAXNATOMS,MAXNEI) :: nei

    integer i, i_ind, j, j_ind
    REAL(8) :: xi, yi, zi, xij, yij, zij, rij2

    numnei = 0  ! Vectorial assignment
    DO i=1, NATOMS
!      i_ind = type(i)
      i_ind=1
      xi = x(i)
      yi = y(i)
      zi = z(i)

      DO j=i+1, NATOMS
!        j_ind = type(j)
        j_ind=1

        xij=x(j)-xi
        yij=y(j)-yi
        zij=z(j)-zi
        xij=xij-bx*ANINT(xij/bx)
        yij=yij-by*ANINT(yij/by)
        zij=zij-bz*ANINT(zij/bz)

        rij2 = xij*xij + yij*yij + zij*zij

 !       IF (rij2 < rcut2(i_ind,j_ind) ) THEN
        IF (rij2 < RSKIN**2) THEN
          numnei(i) = numnei(i) + 1
          numnei(j) = numnei(j) + 1
          nei(i,numnei(i)) = j
          nei(j,numnei(j)) = i
        ENDIF
      END DO
    END DO

  END SUBROUTINE neighbours

!-----------------------------------------------------------------------------------!
!
!  SUBROUTINE potinit()
!     IMPLICIT NONE
!
!  END SUBROUTINE potinit
!
!-----------------------------------------------------------------------------------!

!END MODULE  potential
