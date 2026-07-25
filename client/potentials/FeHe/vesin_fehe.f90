! This file is part of eOn.
!
! SPDX-License-Identifier: BSD-3-Clause
!
! Copyright (c) 2010--present, eOn Development Team
! All rights reserved.
!
! Repo:
! https://github.com/TheochemUI/eOn

!> Half neighbour list for the Fe-He potential, built with vesin.
!!
!! The list is stored in CSR form over the first atom of each pair:
!! `vnl_jlist(vnl_row(i) : vnl_row(i+1)-1)` holds the neighbours `j > i` of
!! atom `i`, in ascending order. Iterating rows in order therefore visits
!! every pair exactly once, in the same (i, j) lexicographic order as a
!! nested `do i; do j > i` sweep.
!!
!! One list serves every species combination: the cutoff passed in is the
!! largest of the Fe-Fe, Fe-He and He-He ranges, and `rhovlc` / `FEFORCE`
!! apply the per-species cutoffs themselves. Those routines also derive the
!! separation vectors, so the minimum image convention and the `r_i - r_j`
!! sign convention stay with the caller.
module vesin_fehe
  use, intrinsic :: iso_c_binding, only: c_double
  use vesin, only: NeighborList
  implicit none

  private
  public :: vesin_fehe_build, vnl_row, vnl_jlist, vnl_npairs

  !> Row starts, `natoms+1` entries, one-based
  integer, allocatable :: vnl_row(:)
  !> Neighbour indices, grouped by row and ascending within a row
  integer, allocatable :: vnl_jlist(:)
  !> Number of pairs in `vnl_jlist`
  integer :: vnl_npairs = 0

  type(NeighborList), save :: nl
  logical, save :: nl_ready = .false.
  real(c_double), save :: nl_cutoff = -1.0_c_double

  real(c_double), allocatable, save :: points(:,:)
  integer, allocatable, save :: cnt(:), rstart(:), scratch(:)

contains

  !> Build the half list for `natoms` atoms within `cutoff`.
  !!
  !! Coordinates arrive as the three separate arrays that the fixed-form
  !! caller keeps. `ax`, `ay`, `az` are the orthorhombic cell lengths.
  !! `ierr` is 0 on success and 1 when vesin fails.
  subroutine vesin_fehe_build(natoms, x0, y0, z0, ax, ay, az, cutoff, ierr)
    integer, intent(in) :: natoms
    real(c_double), intent(in) :: x0(natoms), y0(natoms), z0(natoms)
    real(c_double), intent(in) :: ax, ay, az
    real(c_double), intent(in) :: cutoff
    integer, intent(out) :: ierr

    real(c_double) :: box(3,3)
    integer :: i, j, k, p, np, ip, jp, pos, ins, jv, prev, nout, vstatus

    ierr = 0
    vnl_npairs = 0

    if (.not. nl_ready .or. nl_cutoff /= cutoff) then
      nl = NeighborList(cutoff=cutoff, full=.false., sorted=.true.)
      nl_ready = .true.
      nl_cutoff = cutoff
    end if

    if (allocated(points)) then
      if (size(points, 2) /= natoms) deallocate(points)
    end if
    if (.not. allocated(points)) allocate(points(3, natoms))

    do i = 1, natoms
      points(1, i) = x0(i)
      points(2, i) = y0(i)
      points(3, i) = z0(i)
    end do

    box = 0.0_c_double
    box(1,1) = ax
    box(2,2) = ay
    box(3,3) = az

    call nl%compute(points, box, periodic=[.true., .true., .true.], &
                    status=vstatus)
    if (vstatus /= 0) then
      write(*,*) 'FeHe: vesin neighbour build failed: ', nl%errmsg
      ierr = 1
      return
    end if

    call grow_int(vnl_row, natoms + 1)
    call grow_int(cnt, natoms + 1)
    call grow_int(rstart, natoms + 1)

    np = int(nl%length)

    ! Bucket every pair under min(i,j); a pair repeated across periodic images
    ! lands in the same bucket and is collapsed after sorting.
    cnt(1:natoms) = 0
    do p = 1, np
      ip = int(nl%pairs(1, p)) + 1
      jp = int(nl%pairs(2, p)) + 1
      if (ip == jp) cycle
      i = min(ip, jp)
      cnt(i) = cnt(i) + 1
    end do

    rstart(1) = 1
    do i = 1, natoms
      rstart(i+1) = rstart(i) + cnt(i)
    end do
    call grow_int(scratch, max(rstart(natoms+1) - 1, 1))
    call grow_int(vnl_jlist, max(rstart(natoms+1) - 1, 1))

    cnt(1:natoms) = 0
    do p = 1, np
      ip = int(nl%pairs(1, p)) + 1
      jp = int(nl%pairs(2, p)) + 1
      if (ip == jp) cycle
      i = min(ip, jp)
      j = max(ip, jp)
      pos = rstart(i) + cnt(i)
      scratch(pos) = j
      cnt(i) = cnt(i) + 1
    end do

    nout = 0
    do i = 1, natoms
      vnl_row(i) = nout + 1

      do k = rstart(i) + 1, rstart(i) + cnt(i) - 1
        jv = scratch(k)
        ins = k - 1
        do while (ins >= rstart(i))
          if (scratch(ins) <= jv) exit
          scratch(ins+1) = scratch(ins)
          ins = ins - 1
        end do
        scratch(ins+1) = jv
      end do

      prev = 0
      do k = rstart(i), rstart(i) + cnt(i) - 1
        if (scratch(k) == prev) cycle
        prev = scratch(k)
        nout = nout + 1
        vnl_jlist(nout) = prev
      end do
    end do
    vnl_row(natoms+1) = nout + 1
    vnl_npairs = nout

  end subroutine vesin_fehe_build

  subroutine grow_int(arr, n)
    integer, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: n

    if (allocated(arr)) then
      if (size(arr) < n) deallocate(arr)
    end if
    if (.not. allocated(arr)) allocate(arr(n))
  end subroutine grow_int

end module vesin_fehe
