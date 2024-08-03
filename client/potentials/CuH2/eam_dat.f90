module eam_dat
    use iso_fortran_env, only: real64
    implicit none

    ! Parameters
    integer, parameter :: maxatoms = 500        ! Maximum number of atoms
    integer, parameter :: maxcoo = maxatoms*3   ! Maximum number of coordinates

    ! Variables from the COMMON blocks
    real(real64) :: PL(3)
    real(real64) :: Utot
    real(real64) :: ra(maxcoo)
    real(real64) :: fa(maxcoo)
    integer :: iatclass(maxatoms)
    integer :: ncoo, nim, natoms, nimrp, nimpo
    integer :: nCuCl, nCuQ, nH, netaCu, netaH

    ! EAM variables
    real(real64) :: FCu1, FCu2, FCu3, FCu4, FCu5, FCu6, FCu7, FCu8
    real(real64) :: FH1, FH2, FH3, FH4, FH5
    real(real64) :: DACuCu, alphaACuCu, DBCuCu, alphaBCuCu
    real(real64) :: DAHH, alphaAHH, DBHH, alphaBHH
    real(real64) :: DAHCu, alphaAHCu, DBHCu, alphaBHCu
    real(real64) :: scaleCu, betaACu, betaBCu, gammaCu
    real(real64) :: scaleH, betaAH, betaBH, gammaH
    real(real64) :: rcut, rskin, rcut2, rskin2
    real(real64) :: phicutCuCu, phicutHCu, phicutHH, rhocutCu, rhocutH
    real(real64) :: dphiCuCudr_r, dphiHCudr_r, dphiHHdr_r, drhoCudr_r, drhoHdr_r

end module eam_dat
