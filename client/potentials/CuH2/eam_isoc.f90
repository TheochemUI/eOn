  module eam_wrap
    use iso_c_binding, only: c_double, c_int
    use iso_fortran_env, only: real64
      implicit none
  contains

    subroutine c_force_eam (natms, ndim, box, R, F, U) bind(c)
      integer, intent(in) :: natms(2), ndim
      real(c_double), intent(in) :: U(1), R(ndim), F(ndim), box(3)

      call force_eam(natms,ndim,box,R,F,U)

    end subroutine c_force_eam



  end module eam_wrap
