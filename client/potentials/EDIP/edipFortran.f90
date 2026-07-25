! This file is part of eOn.
!
! SPDX-License-Identifier: BSD-3-Clause
!
! Copyright (c) 2010--present, eOn Development Team
! All rights reserved.
!
! Repo:
! https://github.com/TheochemUI/eOn

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

        use vesin, only: NeighborList

        implicit real*8 (a-h,o-z)

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

! END AA 21-12-2002

        dimension rxyz0(3,nat),fxyz(3,nat),alat(3),box(3,3)
        integer, ALLOCATABLE, DIMENSION(:,:) :: lsta
        integer, ALLOCATABLE, DIMENSION(:) :: lstb
        integer, ALLOCATABLE, DIMENSION(:) :: nnbr
        real*8, ALLOCATABLE, DIMENSION(:,:) :: rel
        real*8, ALLOCATABLE, DIMENSION(:,:) :: txyz
        real*8, ALLOCATABLE, DIMENSION(:,:) :: s2,s3,sz
        integer, ALLOCATABLE, DIMENSION(:) :: num2,num3,numz
        TYPE(NeighborList),SAVE :: nl
        LOGICAL,SAVE :: nl_ready=.false.


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

! the box has to hold at least one cutoff sphere per direction
        ll1=int(alat(1)/cut)
        if (ll1.lt.1) stop 'alat(1) too small'
        ll2=int(alat(2)/cut)
        if (ll2.lt.1) stop 'alat(2) too small'
        ll3=int(alat(3)/cut)
        if (ll3.lt.1) stop 'alat(3) too small'

! full neighbour list from vesin; the calculator lives across calls so that
! vesin reuses its internal buffers
        if (.not.nl_ready) then
          nl=NeighborList(cutoff=cut,full=.true.,sorted=.true., &
                          return_distances=.true.,return_vectors=.true.)
          nl_ready=.true.
        endif

        box=0.d0
        box(1,1)=alat(1)
        box(2,2)=alat(2)
        box(3,3)=alat(3)
        call nl%compute(rxyz0,box,periodic=[.true.,.true.,.true.], &
                        status=ivstat)
        if (ivstat.ne.0) then
          write(6,*) 'EDIP: vesin neighbour build failed: ',nl%errmsg
          stop 1
        endif

! per-atom neighbour counts give the contiguous lsta slices; nnbrx is the
! widest per-atom list, so lstb and rel keep the nnbrx*nat shape that
! subfeniatEDIP declares
        allocate(lsta(2,nat),nnbr(nat))
        nnbr(1:nat)=0
        do ip=1,int(nl%length)
          iat=int(nl%pairs(1,ip))+1
          if (iat.eq.int(nl%pairs(2,ip))+1) cycle
          if (nl%distances(ip).le.0.d0) cycle
          nnbr(iat)=nnbr(iat)+1
        enddo

        nnbrx=max(1,maxval(nnbr(1:nat)))
        allocate(lstb(nnbrx*nat),rel(5,nnbrx*nat))

        indlst=0
        do iat=1,nat
          lsta(1,iat)=indlst+1
          indlst=indlst+nnbr(iat)
          lsta(2,iat)=indlst
        enddo

! rel carries the unit vector of r_iat - r_jat, the distance and its
! inverse; vesin reports r_jat - r_iat with the periodic shift folded in
        nnbr(1:nat)=0
        do ip=1,int(nl%length)
          iat=int(nl%pairs(1,ip))+1
          jat=int(nl%pairs(2,ip))+1
          if (iat.eq.jat) cycle
          tt=nl%distances(ip)
          if (tt.le.0.d0) cycle
          nnbr(iat)=nnbr(iat)+1
          indlst=lsta(1,iat)+nnbr(iat)-1
          tti=1.d0/tt
          lstb(indlst)=jat
          rel(1,indlst)=-nl%vectors(1,ip)*tti
          rel(2,indlst)=-nl%vectors(2,ip)*tti
          rel(3,indlst)=-nl%vectors(3,ip)*tti
          rel(4,indlst)=tt
          rel(5,indlst)=tti
        enddo


!$omp parallel  &
!$omp private(iam,npr,iat,iat1,iat2,lot,istop,tcoord,tcoord2, &
!$omp tener,tener2,txyz,s2,s3,sz,num2,num3,numz,max_nbrs) &
!$omp shared (nat,nnbrx,lsta,lstb,rel,ener,ener2,fxyz,coord,coord2,istopg)

        npr=1
!$       npr=omp_get_num_threads()
        iam=0
!$       iam=omp_get_thread_num()

         max_nbrs=30

        if (npr.ne.1) then
! PARALLEL CASE
! create temporary private scalars for reduction sum on energies and
!        temporary private array for reduction sum on forces
!$omp critical
        allocate(txyz(3,nat),s2(max_nbrs,8),s3(max_nbrs,7),sz(max_nbrs,6),  &
                 num2(max_nbrs),num3(max_nbrs),numz(max_nbrs))
!$omp end critical
! Zero shared accumulators under single (implicit barrier on exit):
! every thread critical-adds into ener/fxyz below, so no thread may start
! accumulating before the clear completes.
!$omp single
        ener=0.d0
        ener2=0.d0
        coord=0.d0
        coord2=0.d0
        istopg=0
        do 121,iat=1,nat
        fxyz(1,iat)=0.d0
        fxyz(2,iat)=0.d0
121     fxyz(3,iat)=0.d0
!$omp end single

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
        istopg=0
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

        deallocate(lsta,lstb,rel,nnbr)

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
