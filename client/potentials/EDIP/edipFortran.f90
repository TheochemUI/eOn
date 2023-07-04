!  SUBROUTINE potinit()
!    IMPLICIT NONE
!  RETURN
!  END SUBROUTINE potinit

!-------------------------------------------------------------------------------!

!    subroutine bazant(nat,alat,rxyz0,fxyz,ener,coord,ener_var,coord_var,count)


    subroutine edip(nat,R,F,ener,ax,ay,az)

!     Evaluates the bazant silicon potential with linear scaling
!     Any publication describing research involving this software should
!     contain the following citations, at once and in this order,
!     1.  M. Z. Bazant and E. Kaxiras, Phys. Rev. Lett. 77, 4370 (1996).
!     2.  M. Z. Bazant, E. Kaxiras, J. F. Justo, Phys. Rev. B 56, 8542 (1997).
!     3.  J. F. Justo, M. Z. Bazant, E. Kaxiras, V. V. Bulatov, and S. Yip,
!           Phys. Rev. B 58, 2539 (1998).
!     4.  S. Goedecker, Comp. Phys. Commun. ??, ?? (2002) once
!         this article has been published or
!         S. Goedecker, cond-mat/0201475 (2002) before it has been published
!
!     Parallelized using OpenMP
! Good Compiler options (last option only if paralellization with OpenMp desired)\
! IBM Power3
!  xlf90_r -O2 -qarch=pwr3 -qtune=pwr3 -qmaxmem=-1 -qsmp=omp
! Dec Alpha
! f90 -arch ev67 -O2 -fast -omp
! Intel Pentium
!  ifc -w -xW -O2 -openmp

!  Copyright (C) 2001-2002 Stefan Goedecker, CEA Grenoble
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!
!     input: - "nat": number of atoms
!            - "alat": lattice constants of the orthorombic box
!               containing the particles
!            - "rxyz0": atomic positions in Angstroem.
!               If an atom is outside the box the program will bring it back
!               into the box by translations through alat
!     output:- "fxyz": forces in eV/A
!            - "ener": total energy in eV
!            - "coord": average coordination number
!            - "ener_var": variance of the energy/atom
!            - "coord_var": variance of the coordination number
!     I/Oput:- "count": is increased by one per call, has to be initialized
!               to 0.d0 before first call of bazant

        implicit real*8 (a-h,o-z)
        logical lgl

!$      interface
!$        integer ( kind=4 ) function omp_get_num_threads ( )
!$        end function omp_get_num_threads
!$      end interface
!$      interface
!$        integer ( kind=4 ) function omp_get_thread_num ( )
!$        end function omp_get_thread_num
!$      end interface

! AA 21-12-2002

        REAL*8,INTENT(IN),DIMENSION(3*nat) :: R
        REAL*8,INTENT(INOUT),DIMENSION(3*nat) :: F
        REAL*8,INTENT(IN) :: ax,ay,az
        INTEGER,INTENT(IN) :: nat
        REAL*8,INTENT(OUT) :: ener
        REAL*8,SAVE :: count=1.0d0

! END AA 21-12-2002

        dimension rxyz0(3,nat),fxyz(3,nat),alat(3)
        real*8, ALLOCATABLE, DIMENSION(:,:) :: rxyz
        integer, ALLOCATABLE, DIMENSION(:,:) :: lsta
        integer, ALLOCATABLE, DIMENSION(:) :: lstb
        integer, ALLOCATABLE, DIMENSION(:) :: lay
        integer, ALLOCATABLE, DIMENSION(:,:,:,:) :: icell
        real*8, ALLOCATABLE, DIMENSION(:,:) :: rel
        real*8, ALLOCATABLE, DIMENSION(:,:) :: txyz
        real*8, ALLOCATABLE, DIMENSION(:,:) :: s2,s3,sz
        integer, ALLOCATABLE, DIMENSION(:) :: num2,num3,numz


! AA 21-12-20002 (2)

        alat(1)=ax
        alat(2)=ay
        alat(3)=az

        j=1
        DO i=1,nat
          rxyz0(:,i)=R(j:j+2)
          j=j+3
        END DO

! END AA 21-12-2002 (2)

!        cut=par_a
        cut= 3.1213820d0 + 1.d-14

!        if (count.eq.0)  open(unit=10,file='bazant.mon',status='unknown')
!        count=count+1.d0

! linear scaling calculation of verlet list
        ll1=int(alat(1)/cut)
        if (ll1.lt.1) stop 'alat(1) too small'
        ll2=int(alat(2)/cut)
        if (ll2.lt.1) stop 'alat(2) too small'
        ll3=int(alat(3)/cut)
        if (ll3.lt.1) stop 'alat(3) too small'

! determine number of threads
        npr=1
!$omp parallel private(iam)  shared (npr)
!$       iam=omp_get_thread_num()
!$       if (iam.eq.0) npr=omp_get_num_threads()
!$omp end parallel

! linear scaling calculation of verlet list

     if (npr.le.1) then !serial if too few processors to gain by parallelizing

! set ncx for serial case, ncx for parallel case set below
        ncx=8
1234    ncx=ncx*2
        allocate(icell(0:ncx,-1:ll1,-1:ll2,-1:ll3))
        do 984,l3=-1,ll3
        do 984,l2=-1,ll2
        do 984,l1=-1,ll1
984     icell(0,l1,l2,l3)=0
        rlc1i=ll1/alat(1)
        rlc2i=ll2/alat(2)
        rlc3i=ll3/alat(3)

        do 983,iat=1,nat

        lgl=.false.
115    continue
        if (rxyz0(1,iat).ge.alat(1)) then
                  if (lgl) then
                    write(10,*) count,' bad x position ', iat,rxyz0(1,iat)
                    rxyz0(1,iat)=modulo(rxyz0(1,iat),alat(1))
                    goto 115
                  endif
            rxyz0(1,iat)=rxyz0(1,iat)-alat(1)
            lgl=.true.
            goto 115
        endif
        if (rxyz0(1,iat).lt.0.d0) then
                  if (lgl) then
                    write(10,*) count,' bad x position ', iat,rxyz0(1,iat)
                    rxyz0(1,iat)=modulo(rxyz0(1,iat),alat(1))
                    goto 115
                  endif
            rxyz0(1,iat)=rxyz0(1,iat)+alat(1)
            lgl=.true.
            goto 115
        endif
        l1=int(rxyz0(1,iat)*rlc1i)

        lgl=.false.
225    continue
        if (rxyz0(2,iat).ge.alat(2)) then
                  if (lgl) then
                    write(10,*) count,' bad y position ', iat,rxyz0(2,iat)
                    rxyz0(2,iat)=modulo(rxyz0(2,iat),alat(2))
                    goto 225
                  endif
            rxyz0(2,iat)=rxyz0(2,iat)-alat(2)
            lgl=.true.
            goto 225
        endif
        if (rxyz0(2,iat).lt.0.d0) then
                  if (lgl) then
                    write(10,*) count,' bad y position ', iat,rxyz0(2,iat)
                    rxyz0(2,iat)=modulo(rxyz0(2,iat),alat(2))
                    goto 225
                  endif
            rxyz0(2,iat)=rxyz0(2,iat)+alat(2)
            lgl=.true.
            goto 225
        endif
        l2=int(rxyz0(2,iat)*rlc2i)

        lgl=.false.
335    continue
        if (rxyz0(3,iat).ge.alat(3)) then
                  if (lgl) then
                    write(10,*) count,' bad z position ', iat,rxyz0(3,iat)
                    rxyz0(3,iat)=modulo(rxyz0(3,iat),alat(3))
                    goto 335
                  endif
            rxyz0(3,iat)=rxyz0(3,iat)-alat(3)
            lgl=.true.
            goto 335
        endif
        if (rxyz0(3,iat).lt.0.d0) then
                  if (lgl) then
                    write(10,*) count,' bad z position ', iat,rxyz0(3,iat)
                    rxyz0(3,iat)=modulo(rxyz0(3,iat),alat(3))
                    goto 225
                  endif
            rxyz0(3,iat)=rxyz0(3,iat)+alat(3)
            lgl=.true.
            goto 335
        endif
        l3=int(rxyz0(3,iat)*rlc3i)

        ii=icell(0,l1,l2,l3)
        ii=ii+1
        icell(0,l1,l2,l3)=ii
        if (ii.gt.ncx) then
        write(10,*) count,'NCX too small',ncx
        deallocate(icell)
        goto 1234
        endif
        icell(ii,l1,l2,l3)=iat
983     continue

     else  ! parallel case

! periodization of particles can be done in parallel
!$omp parallel do shared (alat,nat,rxyz0,count) private(lgl,iat)

        do 5983,iat=1,nat

        lgl=.false.
1155    continue
        if (rxyz0(1,iat).ge.alat(1)) then
                  if (lgl) then
                    write(10,*) count,' bad x position ', iat,rxyz0(1,iat)
                    rxyz0(1,iat)=modulo(rxyz0(1,iat),alat(1))
                    goto 1155
                  endif
            rxyz0(1,iat)=rxyz0(1,iat)-alat(1)
            lgl=.true.
            goto 1155
        endif
        if (rxyz0(1,iat).lt.0.d0) then
                  if (lgl) then
                    write(10,*) count,' bad x position ', iat,rxyz0(1,iat)
                    rxyz0(1,iat)=modulo(rxyz0(1,iat),alat(1))
                    goto 1155
                  endif
            rxyz0(1,iat)=rxyz0(1,iat)+alat(1)
            lgl=.true.
            goto 1155
        endif

        lgl=.false.
2255    continue
        if (rxyz0(2,iat).ge.alat(2)) then
                  if (lgl) then
                    write(10,*) count,' bad y position ', iat,rxyz0(2,iat)
                    rxyz0(2,iat)=modulo(rxyz0(2,iat),alat(2))
                    goto 2255
                  endif
            rxyz0(2,iat)=rxyz0(2,iat)-alat(2)
            lgl=.true.
            goto 2255
        endif
        if (rxyz0(2,iat).lt.0.d0) then
                  if (lgl) then
                    write(10,*) count,' bad y position ', iat,rxyz0(2,iat)
                    rxyz0(2,iat)=modulo(rxyz0(2,iat),alat(2))
                    goto 2255
                  endif
            rxyz0(2,iat)=rxyz0(2,iat)+alat(2)
            lgl=.true.
            goto 2255
        endif

        lgl=.false.
3355    continue
        if (rxyz0(3,iat).ge.alat(3)) then
                  if (lgl) then
                    write(10,*) count,' bad z position ', iat,rxyz0(3,iat)
                    rxyz0(3,iat)=modulo(rxyz0(3,iat),alat(3))
                    goto 3355
                  endif
            rxyz0(3,iat)=rxyz0(3,iat)-alat(3)
            lgl=.true.
            goto 3355
        endif
        if (rxyz0(3,iat).lt.0.d0) then
                  if (lgl) then
                    write(10,*) count,' bad z position ', iat,rxyz0(3,iat)
                    rxyz0(3,iat)=modulo(rxyz0(3,iat),alat(3))
                    goto 3355
                  endif
            rxyz0(3,iat)=rxyz0(3,iat)+alat(3)
            lgl=.true.
            goto 3355
        endif
5983        continue
!$omp end parallel do

! assignment to cell is done serially
! set ncx for parallel case, ncx for serial case set above
        ncx=8
4321    ncx=ncx*2
        allocate(icell(0:ncx,-1:ll1,-1:ll2,-1:ll3))
        do 3984,l3=-1,ll3
        do 3984,l2=-1,ll2
        do 3984,l1=-1,ll1
3984    icell(0,l1,l2,l3)=0

        rlc1i=ll1/alat(1)
        rlc2i=ll2/alat(2)
        rlc3i=ll3/alat(3)

        do 6983,iat=1,nat
        l1=int(rxyz0(1,iat)*rlc1i)
        l2=int(rxyz0(2,iat)*rlc2i)
        l3=int(rxyz0(3,iat)*rlc3i)
        ii=icell(0,l1,l2,l3)
        ii=ii+1
        icell(0,l1,l2,l3)=ii
        if (ii.gt.ncx) then
        write(10,*) count,'NCX too small',ncx
        deallocate(icell)
        goto 4321
        endif
        icell(ii,l1,l2,l3)=iat
6983    continue

    endif


! duplicate all atoms within boundary layer
        laymx=ncx*(2*ll1*ll2+2*ll1*ll3+2*ll2*ll3+4*ll1+4*ll2+4*ll3+8)
        nn=nat+laymx
        allocate(rxyz(3,nn),lay(nn))
        do  iat=1,nat
        lay(iat)=iat
        rxyz(1,iat)=rxyz0(1,iat)
        rxyz(2,iat)=rxyz0(2,iat)
        rxyz(3,iat)=rxyz0(3,iat)
        enddo
        il=nat
! xy plane
        do l2=0,ll2-1
        do l1=0,ll1-1

        in=icell(0,l1,l2,0)
        icell(0,l1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,l1,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,l2,ll3-1)
        icell(0,l1,l2,-1)=in
        do ii=1,in
        i=icell(ii,l1,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo
        enddo


! yz plane
        do l3=0,ll3-1
        do l2=0,ll2-1

        in=icell(0,0,l2,l3)
        icell(0,ll1,l2,l3)=in
        do ii=1,in
        i=icell(ii,0,l2,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,l2,l3)
        icell(0,-1,l2,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo
        enddo


! xz plane
        do l3=0,ll3-1
        do l1=0,ll1-1

        in=icell(0,l1,0,l3)
        icell(0,l1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,l1,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,l1,ll2-1,l3)
        icell(0,l1,-1,l3)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo
        enddo


! x axis
        do l1=0,ll1-1

        in=icell(0,l1,0,0)
        icell(0,l1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,l1,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,0,ll3-1)
        icell(0,l1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,l1,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,l1,ll2-1,0)
        icell(0,l1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,ll2-1,ll3-1)
        icell(0,l1,-1,-1)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo


! y axis
        do l2=0,ll2-1

        in=icell(0,0,l2,0)
        icell(0,ll1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,0,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,l2,ll3-1)
        icell(0,ll1,l2,-1)=in
        do ii=1,in
        i=icell(ii,0,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,l2,0)
        icell(0,-1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,l2,ll3-1)
        icell(0,-1,l2,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo


! z axis
        do l3=0,ll3-1

        in=icell(0,0,0,l3)
        icell(0,ll1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,0,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,0,l3)
        icell(0,-1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,0,ll2-1,l3)
        icell(0,ll1,-1,l3)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,ll2-1,l3)
        icell(0,-1,-1,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo


! corners
        in=icell(0,0,0,0)
        icell(0,ll1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,0,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,0,0)
        icell(0,-1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,ll2-1,0)
        icell(0,ll1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,ll2-1,0)
        icell(0,-1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,0,ll3-1)
        icell(0,ll1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,0,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,0,ll3-1)
        icell(0,-1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,0,ll2-1,ll3-1)
        icell(0,ll1,-1,-1)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,ll2-1,ll3-1)
        icell(0,-1,-1,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        allocate(lsta(2,nat))
        nnbrx=8
2345    nnbrx=3*nnbrx/2
        allocate(lstb(nnbrx*nat),rel(5,nnbrx*nat))

        indlstx=0

!$omp parallel  &
!$omp private(iat,cut2,iam,ii,indlst,l1,l2,l3,myspace,npr) &
!$omp shared (indlstx,nat,nn,nnbrx,ncx,ll1,ll2,ll3,icell,lsta,lstb,lay, &
!$omp rel,rxyz,cut,myspaceout)


        npr=1
!$       npr=omp_get_num_threads()
        iam=0
!$       iam=omp_get_thread_num()

        cut2=cut**2
! assign contiguous portions of the arrays lstb and rel to the threads
        myspace=(nat*nnbrx)/npr
        if (iam.eq.0) myspaceout=myspace
! Verlet list, relative positions
        indlst=0
      do 6000,l3=0,ll3-1
      do 6000,l2=0,ll2-1
      do 6000,l1=0,ll1-1
      do 6600,ii=1,icell(0,l1,l2,l3)
        iat=icell(ii,l1,l2,l3)
        if ( ((iat-1)*npr)/nat .eq. iam) then
!     write(6,*) 'sublstiat:iam,iat',iam,iat
        lsta(1,iat)=iam*myspace+indlst+1
        call sublstiatEDIP(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, &
             rxyz,icell,lstb(iam*myspace+1),lay,rel(1,iam*myspace+1),cut2,indlst)
        lsta(2,iat)=iam*myspace+indlst
!        write(6,'(a,4(x,i3),100(x,i2))') &
!               'iam,iat,lsta',iam,iat,lsta(1,iat),lsta(2,iat), &
!                    (lstb(j),j=lsta(1,iat),lsta(2,iat))
        endif

6600    continue
6000    continue
!$omp critical
        indlstx=max(indlstx,indlst)
!$omp end critical
!$omp end parallel

           if (indlstx.ge.myspaceout) then
               write(10,*) count,'NNBRX too  small', nnbrx
               deallocate(lstb,rel)
               goto 2345
           endif

!$omp parallel  &
!$omp private(iam,npr,iat,iat1,iat2,lot,istop,tcoord,tcoord2, &
!$omp tener,tener2,txyz,s2,s3,sz,num2,num3,numz,max_nbrs) &
!$omp shared (nat,nnbrx,lsta,lstb,rel,ener,ener2,fxyz,coord,coord2,istopg)

        npr=1
!$       npr=omp_get_num_threads()
        iam=0
!$       iam=omp_get_thread_num()

         max_nbrs=30
         istopg=0

        if (npr.ne.1) then
! PARALLEL CASE
! create temporary private scalars for reduction sum on energies and
!        temporary private array for reduction sum on forces
!$omp critical
        allocate(txyz(3,nat),s2(max_nbrs,8),s3(max_nbrs,7),sz(max_nbrs,6),  &
                 num2(max_nbrs),num3(max_nbrs),numz(max_nbrs))
!$omp end critical
        if (iam.eq.0) then
        ener=0.d0
        ener2=0.d0
        coord=0.d0
        coord2=0.d0
        do 121,iat=1,nat
        fxyz(1,iat)=0.d0
        fxyz(2,iat)=0.d0
121     fxyz(3,iat)=0.d0
        endif

! Each thread treats at most lot atoms
        lot=int(float(nat)/float(npr)+.999999999999d0)
        iat1=iam*lot+1
        iat2=min((iam+1)*lot,nat)
!       write(6,*) 'subfeniat:iat1,iat2,iam',iat1,iat2,iam
        call subfeniatEDIP(iat1,iat2,nat,lsta,lstb,rel,tener,tener2,  &
          tcoord,tcoord2,nnbrx,txyz,max_nbrs,istop,  &
          s2(1,1),s2(1,2),s2(1,3),s2(1,4),s2(1,5),s2(1,6),s2(1,7),s2(1,8),  &
          num2,s3(1,1),s3(1,2),s3(1,3),s3(1,4),s3(1,5),s3(1,6),s3(1,7),  &
          num3,sz(1,1),sz(1,2),sz(1,3),sz(1,4),sz(1,5),sz(1,6),numz)

!$omp critical
        ener=ener+tener
        ener2=ener2+tener2
        coord=coord+tcoord
        coord2=coord2+tcoord2
        istopg=istopg+istop
        do 8093,iat=1,nat
        fxyz(1,iat)=fxyz(1,iat)+txyz(1,iat)
        fxyz(2,iat)=fxyz(2,iat)+txyz(2,iat)
        fxyz(3,iat)=fxyz(3,iat)+txyz(3,iat)
8093    continue
        deallocate(txyz,s2,s3,sz,num2,num3,numz)
!$omp end critical

        else
! SERIAL CASE
        iat1=1
        iat2=nat
        allocate(s2(max_nbrs,8),s3(max_nbrs,7),sz(max_nbrs,6),  &
                 num2(max_nbrs),num3(max_nbrs),numz(max_nbrs))
        call subfeniatEDIP(iat1,iat2,nat,lsta,lstb,rel,ener,ener2,  &
          coord,coord2,nnbrx,fxyz,max_nbrs,istopg,  &
          s2(1,1),s2(1,2),s2(1,3),s2(1,4),s2(1,5),s2(1,6),s2(1,7),s2(1,8),  &
          num2,s3(1,1),s3(1,2),s3(1,3),s3(1,4),s3(1,5),s3(1,6),s3(1,7),  &
          num3,sz(1,1),sz(1,2),sz(1,3),sz(1,4),sz(1,5),sz(1,6),numz)
        deallocate(s2,s3,sz,num2,num3,numz)

        endif
!$omp end parallel

!         write(6,*) 'ener,norm force', &
!                    ener,DNRM2(3*nat,fxyz,1)
        if (istopg.gt.0) stop 'DIMENSION ERROR (see WARNING above)'
        ener_var=ener2/nat-(ener/nat)**2
        coord=coord/nat
        coord_var=coord2/nat-coord**2

        j=1
        DO i=1,nat
          F(j:j+2)=fxyz(:,i)
          j=j+3
        END DO

        deallocate(rxyz,icell,lay,lsta,lstb,rel)

        end subroutine edip

!----------------------------------------------------------------------------------!


        subroutine subfeniatEDIP(iat1,iat2,nat,lsta,lstb,rel,ener,ener2,  &
          coord,coord2,nnbrx,ff,max_nbrs,istop,  &
          s2_t0,s2_t1,s2_t2,s2_t3,s2_dx,s2_dy,s2_dz,s2_r,  &
          num2,s3_g,s3_dg,s3_rinv,s3_dx,s3_dy,s3_dz,s3_r,  &
          num3,sz_df,sz_sum,sz_dx,sz_dy,sz_dz,sz_r,numz)
! This subroutine is a modification of a subroutine that is available at
! http://www-math.mit.edu/~bazant/EDIP/ and for which Martin Z. Bazant
! and Harvard University have a 1997 copyright.
! The modifications were done by S. Goedecker on April 10, 2002.
! The routines are included with the permission of M. Bazant into this package.

        implicit none
!  ------------------------- VARIABLE DECLARATIONS -------------------------
          integer iat1,iat2,nat
          real*8 ener,ener2,coord,coord2
          real*8 xarg,coord_iat,ener_iat
          real*8 ff(3,nat)

        real*8 par_cap_A,par_cap_B,par_rh,par_a,par_sig,par_lam,par_gam, &
               par_b,par_c,par_delta,par_mu,par_Qo,par_palp, &
               par_bet,par_alp,par_bg,par_eta,u1,u2,u3,u4,u5

          integer nnbrx,max_nbrs,istop
          integer lsta(2,nat),lstb(nnbrx*nat)
          real*8 rel(5,nnbrx*nat)

          integer i,j,k,l,n
          real*8 dx,dy,dz,r
          real*8 rinv,rmainv,xinv,xinv3,den,Z,fZ
          real*8 dV2j,dV2ijx,dV2ijy,dV2ijz,pZ,dp
          real*8 temp0,temp1
          real*8 Qort,muhalf
          real*8 rmbinv,winv,dwinv,tau,dtau,lcos,x,H,dHdx,dhdl
          real*8 dV3rij,dV3rijx,dV3rijy,dV3rijz
          real*8 dV3rik,dV3rikx,dV3riky,dV3rikz
          real*8 dV3l,dV3ljx,dV3ljy,dV3ljz,dV3lkx,dV3lky,dV3lkz
          real*8 dV2dZ,dxdZ,dV3dZ
          real*8 dEdrl,dEdrlx,dEdrly,dEdrlz
          real*8 bmc,cmbinv
          real*8 fjx,fjy,fjz,fkx,fky,fkz

          real*8 s2_t0(max_nbrs)
          real*8 s2_t1(max_nbrs)
          real*8 s2_t2(max_nbrs)
          real*8 s2_t3(max_nbrs)
          real*8 s2_dx(max_nbrs)
          real*8 s2_dy(max_nbrs)
          real*8 s2_dz(max_nbrs)
          real*8 s2_r(max_nbrs)
          integer n2
!   size of s2[]
          integer num2(max_nbrs)
!   atom ID numbers for s2[]

          real*8 s3_g(max_nbrs)
          real*8 s3_dg(max_nbrs)
          real*8 s3_rinv(max_nbrs)
          real*8 s3_dx(max_nbrs)
          real*8 s3_dy(max_nbrs)
          real*8 s3_dz(max_nbrs)
          real*8 s3_r(max_nbrs)

          integer n3
!   size of s3[]
          integer num3(max_nbrs)
!   atom ID numbers for s3[]

          real*8 sz_df(max_nbrs)
          real*8 sz_sum(max_nbrs)
          real*8 sz_dx(max_nbrs)
          real*8 sz_dy(max_nbrs)
          real*8 sz_dz(max_nbrs)
          real*8 sz_r(max_nbrs)
          integer nz
!   size of sz[]
          integer numz(max_nbrs)
!   atom ID numbers for sz[]

          integer nj,nk,nl
!   indices for the store arrays

!   EDIP parameters
          par_cap_A = 5.6714030d0
          par_cap_B = 2.0002804d0
          par_rh = 1.2085196d0
          par_a = 3.1213820d0
          par_sig = 0.5774108d0
          par_lam = 1.4533108d0
          par_gam = 1.1247945d0
          par_b = 3.1213820d0
          par_c = 2.5609104d0
          par_delta = 78.7590539d0
          par_mu = 0.6966326d0
          par_Qo = 312.1341346d0
          par_palp = 1.4074424d0
          par_bet = 0.0070975d0
          par_alp = 3.1083847d0

          u1 = -0.165799d0
          u2 = 32.557d0
          u3 = 0.286198d0
          u4 = 0.66d0

          par_bg=par_a
          par_eta = par_delta/par_Qo

          do i=1, nat
            ff(1,i) = 0.0d0
            ff(2,i) = 0.0d0
            ff(3,i) = 0.0d0
          end do

          coord=0.d0
          coord2=0.d0
          ener=0.d0
          ener2=0.d0
          istop=0


!   COMBINE COEFFICIENTS

          Qort = sqrt(par_Qo)
          muhalf = par_mu*0.5D0
          u5 = u2*u4
          bmc = par_b-par_c
          cmbinv = 1.0D0/(par_c-par_b)



!  --- LEVEL 1: OUTER LOOP OVER ATOMS ---

          do 1000, i= iat1, iat2

!   RESET COORDINATION AND NEIGHBOR NUMBERS

            coord_iat=0.d0
            ener_iat=0.d0
            Z = 0.0d0
            n2 = 1
            n3 = 1
            nz = 1


!  --- LEVEL 2: LOOP PREPASS OVER PAIRS ---

            do n=lsta(1,i),lsta(2,i)
              j=lstb(n)


!   PARTS OF TWO-BODY INTERACTION r<par_a

                num2(n2) = j
                dx = -rel(1,n)
                dy = -rel(2,n)
                dz = -rel(3,n)
                r=rel(4,n)
                rinv=rel(5,n)
                rmainv = 1.d0/(r-par_a)
                s2_t0(n2) = par_cap_A*dexp(par_sig*rmainv)
                s2_t1(n2) = (par_cap_B*rinv)**par_rh
                s2_t2(n2) = par_rh*rinv
                s2_t3(n2) = par_sig*rmainv*rmainv
                s2_dx(n2) = dx
                s2_dy(n2) = dy
                s2_dz(n2) = dz
                 s2_r(n2) = r
                n2 = n2 + 1
                if (n2.gt.max_nbrs) then
                write(6,*) 'WARNING enlarge max_nbrs'
                istop=1
                return
                endif

! coordination number calculated with soft cutoff between first and
! second nearest neighbor
        if (r.le.2.36d0) then
        coord_iat=coord_iat+1.d0
        else if (r.ge.3.83d0) then
        else
        xarg=(r-2.36d0)*(1.d0/(3.83d0-2.36d0))
        coord_iat=coord_iat+(2*xarg+1.d0)*(xarg-1.d0)**2
        endif


!   RADIAL PARTS OF THREE-BODY INTERACTION r<par_b

                if(r .lt. par_bg)  then

                  num3(n3) = j
                  rmbinv = 1.d0/(r-par_bg)
                  temp1 = par_gam*rmbinv
                  temp0 = dexp(temp1)
                  s3_g(n3) = temp0
                  s3_dg(n3) = -rmbinv*temp1*temp0
                  s3_dx(n3) = dx
                  s3_dy(n3) = dy
                  s3_dz(n3) = dz
                  s3_rinv(n3) = rinv
                  s3_r(n3) = r
                  n3 = n3 + 1
                  if (n3.gt.max_nbrs) then
                  write(6,*) 'WARNING enlarge max_nbrs'
                  istop=1
                  return
                  endif


!   COORDINATION AND NEIGHBOR FUNCTION par_c<r<par_b

                  if(r .lt. par_b) then
                    if(r .lt. par_c) then
                    Z = Z + 1.d0
                   else
                    xinv = bmc/(r-par_c)
                    xinv3 = xinv*xinv*xinv
                    den = 1.d0/(1 - xinv3)
                    temp1 = par_alp*den
                    fZ = dexp(temp1)
                    Z = Z + fZ
                    numz(nz) = j
                    sz_df(nz) = fZ*temp1*den*3.d0*xinv3*xinv*cmbinv
!   df/dr
                    sz_dx(nz) = dx
                    sz_dy(nz) = dy
                    sz_dz(nz) = dz
                    sz_r(nz) = r
                    nz = nz + 1
                    if (nz.gt.max_nbrs) then
                    write(6,*) 'WARNING enlarge max_nbrs'
                    istop=1
                    return
                    endif
                   end if
!  r < par_C
                  end if
!  r < par_b
                  end if
!  r < par_bg
              end do


!   ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES

              do nl=1, nz-1
                sz_sum(nl)=0.d0
              end do


!   ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION

              temp0 = par_bet*Z
              pZ = par_palp*dexp(-temp0*Z)
!   bond order
              dp = -2.d0*temp0*pZ
!   derivative of bond order



!  --- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---

            do nj=1, n2-1

              temp0 = s2_t1(nj) - pZ


!   two-body energy V2(rij,Z)

              ener_iat = ener_iat + temp0*s2_t0(nj)

!   two-body forces

              dV2j = - s2_t0(nj) * (s2_t1(nj)*s2_t2(nj) + temp0 * s2_t3(nj))
!   dV2/dr
              dV2ijx = dV2j * s2_dx(nj)
              dV2ijy = dV2j * s2_dy(nj)
              dV2ijz = dV2j * s2_dz(nj)
              ff(1,i) = ff(1,i) + dV2ijx
              ff(2,i) = ff(2,i) + dV2ijy
              ff(3,i) = ff(3,i) + dV2ijz
              j = num2(nj)
              ff(1,j) = ff(1,j) - dV2ijx
              ff(2,j) = ff(2,j) - dV2ijy
              ff(3,j) = ff(3,j) - dV2ijz



!  --- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---

              dV2dZ = - dp * s2_t0(nj)
              do nl=1, nz-1
                 sz_sum(nl) =  sz_sum(nl) + dV2dZ
              end do

            end do


!   COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION

              winv = Qort*exp(-muhalf*Z)
!   inverse width of angular function
              dwinv = -muhalf*winv
!   its derivative
              temp0 = exp(-u4*Z)
              tau = u1+u2*temp0*(u3-temp0)
!   -cosine of angular minimum
              dtau = u5*temp0*(2*temp0-u3)
!   its derivative


!  --- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---

            do nj=1, n3-2

              j=num3(nj)


!  --- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---

              do nk=nj+1, n3-1

                k=num3(nk)


!   angular function h(l,Z)

                lcos=s3_dx(nj)*s3_dx(nk)+s3_dy(nj)*s3_dy(nk)+s3_dz(nj)*s3_dz(nk)
                x = (lcos + tau)*winv
                temp0 = exp(-x*x)

                H = par_lam*(1 - temp0 + par_eta*x*x)
                dHdx = 2*par_lam*x*(temp0 + par_eta)

                dhdl = dHdx*winv


!   three-body energy

                temp1 = s3_g(nj) * s3_g(nk)
                ener_iat = ener_iat + temp1*H


!   (-) radial force on atom j

                dV3rij = s3_dg(nj) * s3_g(nk) * H
                dV3rijx = dV3rij * s3_dx(nj)
                dV3rijy = dV3rij * s3_dy(nj)
                dV3rijz = dV3rij * s3_dz(nj)
                fjx = dV3rijx
                fjy = dV3rijy
                fjz = dV3rijz


!   (-) radial force on atom k

                dV3rik = s3_g(nj) * s3_dg(nk) * H
                dV3rikx = dV3rik * s3_dx(nk)
                dV3riky = dV3rik * s3_dy(nk)
                dV3rikz = dV3rik * s3_dz(nk)
                fkx = dV3rikx
                fky = dV3riky
                fkz = dV3rikz


!   (-) angular force on j

                dV3l = temp1*dhdl
                dV3ljx = dV3l * (s3_dx(nk) - lcos * s3_dx(nj)) * s3_rinv(nj)
                dV3ljy = dV3l * (s3_dy(nk) - lcos * s3_dy(nj)) * s3_rinv(nj)
                dV3ljz = dV3l * (s3_dz(nk) - lcos * s3_dz(nj)) * s3_rinv(nj)
                fjx = fjx + dV3ljx
                fjy = fjy + dV3ljy
                fjz = fjz + dV3ljz


!   (-) angular force on k

                dV3lkx = dV3l * (s3_dx(nj) - lcos * s3_dx(nk)) * s3_rinv(nk)
                dV3lky = dV3l * (s3_dy(nj) - lcos * s3_dy(nk)) * s3_rinv(nk)
                dV3lkz = dV3l * (s3_dz(nj) - lcos * s3_dz(nk)) * s3_rinv(nk)
                fkx = fkx + dV3lkx
                fky = fky + dV3lky
                fkz = fkz + dV3lkz


!   apply radial + angular forces to i, j, k

                ff(1,j) = ff(1,j) - fjx
                ff(2,j) = ff(2,j) - fjy
                ff(3,j) = ff(3,j) - fjz
                ff(1,k) = ff(1,k) - fkx
                ff(2,k) = ff(2,k) - fky
                ff(3,k) = ff(3,k) - fkz
                ff(1,i) = ff(1,i) + fjx + fkx
                ff(2,i) = ff(2,i) + fjy + fky
                ff(3,i) = ff(3,i) + fjz + fkz



!   prefactor for 4-body forces from coordination
                  dxdZ = dwinv*(lcos + tau) + winv*dtau
                  dV3dZ = temp1*dHdx*dxdZ


!  --- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---

                  do nl=1, nz-1
                    sz_sum(nl) = sz_sum(nl) + dV3dZ
                  end do
              end do
            end do


!  --- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---

            do nl=1, nz-1

                dEdrl = sz_sum(nl) * sz_df(nl)
                dEdrlx = dEdrl * sz_dx(nl)
                dEdrly = dEdrl * sz_dy(nl)
                dEdrlz = dEdrl * sz_dz(nl)
                ff(1,i) = ff(1,i) + dEdrlx
                ff(2,i) = ff(2,i) + dEdrly
                ff(3,i) = ff(3,i) + dEdrlz
                l = numz(nl)
                ff(1,l) = ff(1,l) - dEdrlx
                ff(2,l) = ff(2,l) - dEdrly
                ff(3,l) = ff(3,l) - dEdrlz


            end do

        coord=coord+coord_iat
        coord2=coord2+coord_iat**2
        ener = ener + ener_iat
        ener2 = ener2 + ener_iat**2

1000      continue


        return
        end

!---------------------------------------------------------------------------------------!

        subroutine sublstiatEDIP(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, &
                   rxyz,icell,lstb,lay,rel,cut2,indlst)
! finds the neighbours of atom iat (specified by lsta and lstb) and and
! the relative position rel of iat with respect to these neighbours
        implicit real*8 (a-h,o-z)
        dimension rxyz(3,nn),lay(nn),icell(0:ncx,-1:ll1,-1:ll2,-1:ll3), &
                  lstb(0:myspace-1),rel(5,0:myspace-1)

        do 6363,k3=l3-1,l3+1
        do 6363,k2=l2-1,l2+1
        do 6363,k1=l1-1,l1+1
        do 6363,jj=1,icell(0,k1,k2,k3)
          jat=icell(jj,k1,k2,k3)
          if (jat.eq.iat) goto 6363
          xrel= rxyz(1,iat)-rxyz(1,jat)
          yrel= rxyz(2,iat)-rxyz(2,jat)
          zrel= rxyz(3,iat)-rxyz(3,jat)
          rr2=xrel**2 + yrel**2 + zrel**2
          if ( rr2 .le. cut2 ) then
           indlst=min(indlst,myspace-1)
           lstb(indlst)=lay(jat)
!        write(6,*) 'iat,indlst,lay(jat)',iat,indlst,lay(jat)
           tt=sqrt(rr2)
           tti=1.d0/tt
           rel(1,indlst)=xrel*tti
           rel(2,indlst)=yrel*tti
           rel(3,indlst)=zrel*tti
           rel(4,indlst)=tt
           rel(5,indlst)=tti
           indlst= indlst+1
          endif
6363        continue

        return
        end
