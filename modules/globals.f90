!-----------------------------------------------------------------------
!> @file    globals.f90
!! @brief   This module contains our project's global variables and some
!!          useful math and physical constants.
!! @author  bgeneto
!! @date    20150628
!-----------------------------------------------------------------------
module globals

    use working_precision

    implicit none
    public

    !-------------------------------------------------------------------
    ! Declaration of global namelists
    ! (mainly those contained in the input data file)
    !-------------------------------------------------------------------
    ! mesh-related variables and parameters
    integer  :: n          !< number of sites
    integer  :: nu         !< number of electrons up
    integer  :: nd         !< number of electrons down
    !> mesh-related namelist
    namelist /mesh/ n, nu, nd

    ! impurity-related variables and parameters
    integer  :: imp1       !< site number for the first impurity positioning
    integer  :: imp2       !< site number for the second impurity positioning
    integer  :: impon      !< impurity switch on time
    integer  :: impoff     !< impurity switch off time
    real(wp) :: timp1      !< first impurity strength (eV)
    real(wp) :: timp2      !< second impurity strength (eV)
    !> impurity-related namelist
    namelist /impurity/ imp1, imp2, timp1, timp2, impon, impoff

    ! iteration-related variables
    integer  :: niter      !< total number of iteration (total_time = 4*iteration*delt in femtosecond)
    real(wp) :: ierror     !< attained precision in each iteration
    real(wp) :: delt       !< integration time step
    ! iteration-related namelist
    namelist /iteration/ niter, ierror, delt

    ! electric field-related variables
    integer  :: ifieldon   !< electric field switch on time
    integer  :: ifieldoff  !< electric field switch off time
    real(wp) :: efield     !< electric field strength (in units of 1.3 v/Å or 0.01 = 1.3mv/Å; max 2.0 mv/Å)
    !> electric field-related namelist
    namelist /field/ ifieldon, ifieldoff, efield

    ! other physical parameters and variables
    real(wp) :: fcnst       !< force constant (K in eV/Å²)
    real(wp) :: alp         !< electron-phonon coupling constant (alpha in units of eV/Å)
    real(wp) :: t0          !< hopping term or transfer integral (eV)
    real(wp) :: u           !< onsite coulomb or onsite electron-electron interaction
    real(wp) :: v           !< intersite coulomb interaction or intersite electron-electron interaction
    real(wp) :: delta0      !< Brazovskii-Kirova symmetry breaking term. for cis (delta0 = 0.05), trans (delta0 = 0)
    real(wp) :: tau         !< switching smoothness parameter
    real(wp) :: tomega      !< omega/t0 (see ref ?)
    real(wp) :: gammal      !< damping term constant in fluctuation-dissipation theorem (viscosity term in Langevin's equation)
    real(wp) :: flut        !< fluctuation magnitude
    real(wp) :: temp        !< temperature in kelvins
    !> physical parameters namelist
    namelist /physical/ fcnst, alp, t0, u, v, delta0, tau, tomega, gammal, flut, temp

    !> output file names:
    character(len=255), target ::   exciton_polaron_fname,      &
                                    up_occ_fname,               &
                                    down_occ_fname,             &
                                    updown_occ_fname,           &
                                    lattice_displacement_fname, &
                                    energy_spectrum_up_fname,   &
                                    energy_spectrum_down_fname, &
                                    energies_fname,             &
                                    charge_density_fname,       &
                                    spin_density_fname,         &
                                    temperature_fname,          &
                                    site_velocity_fname
    !> filenames namelist
    namelist /filenames/            exciton_polaron_fname,      &
                                    up_occ_fname,               &
                                    down_occ_fname,             &
                                    updown_occ_fname,           &
                                    lattice_displacement_fname, &
                                    energy_spectrum_up_fname,   &
                                    energy_spectrum_down_fname, &
                                    energies_fname,             &
                                    charge_density_fname,       &
                                    spin_density_fname,         &
                                    temperature_fname,          &
                                    site_velocity_fname
    !====================================================================

    integer  :: ne              !< total number of electrons (nd + nu)
    real(wp) :: rlambda, alpha  !< auxiliary variables

    !-------------------------------------------------------------------
    ! Declaration of global constants (parameters)
    !-------------------------------------------------------------------
    ! useful numbers name
    real(wp), parameter :: zero  = 0.0_wp        !< floating point number zero
    real(wp), parameter :: third = 1.0_wp/3.0_wp !< one third
    real(wp), parameter :: half  = 0.5_wp        !< one and a half
    real(wp), parameter :: one   = 1.0_wp        !< floating point number one
    real(wp), parameter :: two   = 2.0_wp        !< floating point number two
    real(wp), parameter :: sqrt2 = sqrt(2.0_wp)  !< square root of two
    real(wp), parameter :: sqrt3 = sqrt(3.0_wp)  !< square root of three

    ! useful physical and mathematical constants and parameters
    real(wp), parameter :: pi = 4.0_wp*atan(one) !< mathematical PI constant

    ! useful output informational words
    !> default fatal error string
    character(len = *), parameter :: fatal_str = '***FATAL: '
    !< default error string
    character(len = *), parameter :: err_str   = '***ERROR: '
    !> default warning string
    character(len = *), parameter :: warn_str  = '***WARNING: '
    !> default information string
    character(len = *), parameter :: info_str  = '***INFO: '
    !> continuation string (for printing alignment purposes only)
    character(len = *), parameter :: cont_str  = '          '
    !====================================================================

    !-------------------------------------------------------------------
    ! Declaration of global allocatable arrays
    !-------------------------------------------------------------------
    real(wp), allocatable, dimension(:) :: ocd, ocu  !< occupation (up and down) distribution arrays

end module globals
