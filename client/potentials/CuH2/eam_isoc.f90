! This file is part of eOn.
!
! SPDX-License-Identifier: BSD-3-Clause
!
! Copyright (c) 2010--present, eOn Development Team
! All rights reserved.
!
! Repo:
! https://github.com/TheochemUI/eOn
  module eam_wrap
    use iso_c_binding, only: c_double, c_int
    use iso_fortran_env, only: real64
      implicit none
  contains

    subroutine c_force_eam (natms, ndim, box, R, F, U) bind(c)
      integer(c_int), intent(in) :: natms(2)
      integer(c_int), value, intent(in) :: ndim
      real(c_double), intent(in) :: box(3)
      real(c_double), intent(inout) :: U(1), R(ndim), F(ndim)

      ! write(*,*)"natms", natms
      ! write(*,*)"ndim", ndim
      ! write(*,*)"box", box
      ! write(*,*)"R", R
      call force_eam(natms, ndim, box, R, F, U)

    end subroutine c_force_eam



  end module eam_wrap
