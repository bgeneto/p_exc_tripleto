! Linux build command line:
!
!-----------------------------------------------------------------------
!> @file    ssh.f90
!! @brief   Su-Schrieffer-Heeger model program.
!! @details Exciton dynamics.
!! @author  ribeirojr
!! @author  bgeneto
!! @author  fabinhofis
!! @date    20150623
!! @todo    Change some physical units to a more appropriated form.
!! @note    Compiled with: <tt>gfortran -O3 ssh.f90 -o ssh -lpthread -llapack95 -lopenblas -I/opt/OpenBLAS/include -L/opt/OpenBLAS/lib</tt>
!-----------------------------------------------------------------------
program SSH
    use utils              !< useful procedures
    use globals            !< global variables and constants module
    use omp_lib            !< OpenMP Fortran module
    use io_error_checking  !< file and array io checking module
    use working_precision  !< real and complex working precision module

    ! lapack95 interface library for the eigen problem solution
    !use f95_lapack, only: la_lamch, la_heevr   !< use these lapack95 subroutines interfaces

    ! every variable must be explicit declared!
    implicit none

    integer  :: tid = 0   !< OpenMP thread ID
    integer  :: nthreads  !< OpenMP number of threads available

    real(wp) :: ptime0, ptime   !< initial processing time in seconds

    !-------------------------------------------------------------------
    ! Program execution part
    !-------------------------------------------------------------------
    call cpu_time(ptime0)
    !!$ tid = omp_get_thread_num()
    !!$ nthreads = omp_get_num_threads()

    call read_input_data("input/parameters.txt")

    call init_unit_conversion()

    call occupation()
    call stationary()

    call cpu_time(ptime)
    write(*,'(a,f0.6,a)') "Total execution time: ", (ptime - ptime0), " s"

end program SSH

!> Initializes important variables and performs some physical units conversion.
!! @todo unit conversions
subroutine init_unit_conversion()
    use globals
    use utils
    use working_precision
    implicit none

    ne = nu + nd
    rlambda = two*alp/fcnst  !< dimensionless parameter
    alpha   = alp/t0         !< auxiliary constant

    if (imp1 > n .or. imp2 > n) then
        print*, fatal_str, "imp1 (=", string(imp1), ") or imp2 (=", string(imp2), &
                ") value cannot be larger then n (=", string(n), ")!"
    end if

end subroutine

!> Reads the input data file by using a Fortran feature called \c namelist.
!! @param[in] filename Input data filename.
subroutine read_input_data(filename)
    use globals
    use io_error_checking
    character(len=*), intent(in) :: filename
    unt = free_unit()
    open(unt, file=trim(filename), status='old', action='read', iostat=ios, err=100, iomsg=io_err_msg)
    read(unt, nml=mesh, iostat=ios, err=100)
    read(unt, nml=physical, iostat=ios, err=100)
    read(unt, nml=impurity, iostat=ios, err=100)
    read(unt, nml=iteration, iostat=ios, err=100)
    read(unt, nml=field, iostat=ios, err=100)
    read(unt, nml=filenames, iostat=ios, err=100)
    close(unt, iostat=ios, err=100)
    100 call io_error_check(unt, ios, io_err_msg)
end subroutine read_input_data


!> Build the electronic occupation arrays
!! @param[in] ocd Electron down occupation array
!! @param[in] ocu Electron up occupation array
!! @todo Allow user to enter the occupation scheme in an easy way
subroutine occupation()
    use utils
    use globals
    use io_error_checking
    use working_precision
    implicit none

    integer :: i
    character(len=*), parameter :: occfmt = '(4d15.7)' !< output format

    if ( .not.allocated(ocd) .and. .not.allocated(ocu) ) then
        allocate ( ocu(ne), ocd(ne), stat=ios )
        call io_error_check(ios, 'ocu, ocd')
    end if

    !$omp parallel workshare if (ne > 500) default(shared)
        ! electron up occupation rule
        ocu(1:nu) = one
        ocu(nu+1:ne) = zero
        ! electron down occupation rule
        ocd(1:nd) = one
        ocd(nd+1:ne) = zero
    !$omp end parallel workshare

    ! write (to file) up occupation scheme
    unt = free_unit()
    open(unt, file=trim(up_occ_fname), status='replace', action='write', iostat=ios, err=100, iomsg=io_err_msg)
    write(unt, '(i0,a,f0.7)', iostat=ios, err=100) (i, char(9), ocu(i), i = 1, n)
    close(unt, iostat=ios, err=100)

    call gplot(trim(up_occ_fname), title='occupation (up)', plotdir='plot')

    ! write (to file) down occupation scheme
    unt = free_unit()
    open(unt, file=trim(down_occ_fname), status='replace', action='write', iostat=ios, err=100, iomsg=io_err_msg)
    write(unt, occfmt, iostat=ios, err=100) ocd
    close(unt, iostat=ios, err=100 )

    ! write (to file) up+down occupation scheme
    unt = free_unit()
    open(unt, file=trim(updown_occ_fname), status='replace', action='write', iostat=ios, err=100, iomsg=io_err_msg)
    write(unt, occfmt, iostat=ios, err=100) ocu+ocd
    close(unt, iostat=ios, err=100)

    100 call io_error_check(unt, ios, io_err_msg)

end subroutine occupation

!> Stationary procedure
! @todo Comment every variable used in this procedure
subroutine stationary()
    use utils
    use globals
    use omp_lib
    use io_error_checking
    use working_precision
    implicit none
    real(wp), allocatable, dimension(:)   :: y1, ru, rd, ru1, rd1, tu, td, tu1, td1, lw, v2, eu, ed, y, yv, chden, sden
    real(wp), allocatable, dimension(:,:) :: w, aur, aui, adr, adi, evur, evdr, evui, evdi

    integer  :: i, j, ier, tid, iterationscf = 0
    real(wp) :: ccr, cct, ccy, ddr, ddt, ddy, ee, el, edy, et, eps, uterm, vterm, vterm1, vterm2, ytot, yerror=one, eerror=one


    allocate( y1(n), ru(n), rd(n), ru1(n), rd1(n), tu(n), td(n), tu1(n), &
              td1(n), lw(n), v2(n), eu(n), ed(n), y(n), yv(n), chden(n), sden(n), stat=ios )
    call io_error_check(ios, 'arrays(n)')
    allocate ( w(n,6), aur(n,n), aui(n,n), adr(n,n), adi(n,n), evur(n,n), evdr(n,n), evui(n,n), evdi(n,n), stat=ios )
    call io_error_check(ios, 'arrays(n,n)')

    unt = free_unit()
    open(unt, file=trim(exciton_polaron_fname), status='old', action='read', iostat=ios, err=100, iomsg=io_err_msg)
    read(unt, *, iostat=ios, err=100) ( y(i), i = 1, n )
    close(unt, status='keep', iostat=ios, err=100)

    ru = zero; rd = zero; tu = zero; td = zero

    !> beginning of scf iteration
    scf: do while (eerror + yerror > ierror)
        aur = zero; adr = zero; aui = zero; adi = zero

        forall (i=1:n-1)
            aur(i+1,i) = -(one + delta0*(-one)**i) + alpha*y(i) - v*tu(i)
            adr(i+1,i) = -(one + delta0*(-one)**i) + alpha*y(i) - v*td(i)
        end forall

        aur(n,1) = -(one + delta0*(-one)**n) + alpha*y(n) - v*tu(n)
        adr(n,1) = -(one + delta0*(-one)**n) + alpha*y(n) - v*td(n)

        i = 1
        vterm    = ru(i+1) + ru(n) + rd(i+1) + rd(n) - two
        aur(i,i) = u*(rd(i) - half) + v*vterm
        adr(i,i) = u*(ru(i) - half) + v*vterm

        do i = 2, n-1
            vterm    = ru(i+1) + ru(i-1) + rd(i+1) + rd(i-1) - two
            aur(i,i) = u*(rd(i) - half) + v*vterm
            adr(i,i) = u*(ru(i) - half) + v*vterm
        end do

        i=n
        vterm    = ru(1) + ru(n-1) + rd(1) + rd(n-1) - two
        aur(i,i) = u*(rd(i) - half) + v*vterm
        adr(i,i) = u*(ru(i) - half) + v*vterm

        !< impurity positioning
        aur(imp1,imp1) = aur(imp1,imp1) + timp1
        adr(imp1,imp1) = adr(imp1,imp1) + timp1
        aur(imp2,imp2) = aur(imp2,imp2) + timp2
        adr(imp2,imp2) = adr(imp2,imp2) + timp2

        if (iterationscf .eq. 5) write(10,'(20f8.4,a)') ( ( aur(i,j), j = 1, n), char(9), i = 1, n )
        if (iterationscf .eq. 5) write(11,'(20f8.4,a)') ( ( adr(i,j), j = 1, n), char(9), i = 1, n )

        eps = 1.0d-15
        call eigrs(aur,n,n,-n,n,eps,w,lw,eu,evur,ier)
        call eigrs(adr,n,n,-n,n,eps,w,lw,ed,evdr,ier)


        y1(:)    = zero
        chden(:) = zero
        sden(:)  = zero

        do i = 1, n-1
            do j = 1, ne
                y1(i) = y1(i) - evur(i,j)*evur(i+1,j)*ocu(j) - evdr(i,j)*evdr(i+1,j)*ocd(j)
            end do
        end do

        do j = 1, ne
            y1(n) = y1(n) - evur(n,j)*evur(1,j)*ocu(j) - evdr(n,j)*evdr(1,j)*ocd(j)
        end do

        ytot = sum(y1)/real(n)

        y1 = (y1 - ytot)*rlambda

        ru1 = zero; rd1 = zero; tu1 = zero; td1 = zero

        do i = 1, n
            do j = 1, ne
                ru1(i) = ru1(i) + evur(i,j)*evur(i,j)*ocu(j)
                rd1(i) = rd1(i) + evdr(i,j)*evdr(i,j)*ocd(j)
            end do
        end do

        do i = 1, n-1
            do j = 1, ne
                tu1(i) = tu1(i) + evur(i+1,j)*evur(i,j)*ocu(j)
                td1(i) = td1(i) + evdr(i+1,j)*evdr(i,j)*ocd(j)
            end do
        end do

        do j = 1, ne
            tu1(n) = tu1(n) + evur(1,j)*evur(n,j)*ocu(j)
            td1(n) = td1(n) + evdr(1,j)*evdr(n,j)*ocd(j)
        end do


        !verification of convergence/consistency (below)

        ccy = zero; ddy = zero

        do i = 1, n
            ccy = ccy + y1(i)**2
            ddy = ddy + (y(i) - y1(i))**2
        end do

        yerror  = ddy/ccy
        ccr = zero; ddr = zero

        do i = 1, n
            ccr = ccr + ru1(i)**2 + rd1(i)**2
            ddr = ddr + (ru(i) - ru1(i))**2 + (rd(i) - rd1(i))**2
        end do

        cct = zero; ddt = zero
        eerror = (ddr + ddt)/(ccr + cct)

        do i = 1, n
            cct = cct + tu1(i)**2 + td1(i)**2
            ddt = ddt + (tu(i) - tu1(i))**2 + (td(i) - td1(i))**2
        end do

        iterationscf = iterationscf + 1

        write(*,110) iterationscf, yerror+eerror, yerror, ccy, eerror

        ! update new values
        y  = y1
        ru = ru1
        rd = rd1
        tu = tu1
        td = td1

    end do scf ! end of SCF iteration

    evui = zero; evdi = zero
    yv = zero

    unt = free_unit()
    open(unt, file=trim(lattice_displacement_fname), status='replace', action='write', iostat=ios, err=100, iomsg=io_err_msg)
    write(unt, 120, iostat=ios, err=100)
    write(unt, 190, iostat=ios, err=100) ( y(i), i = 1, n )
    write(unt, 140, iostat=ios, err=100)
    write(unt, 190, iostat=ios, err=100) ( eu(i), i = 1, n )
    write(unt, 150, iostat=ios, err=100)
    write(unt, 190, iostat=ios, err=100) ( ed(i), i = 1, n )
    close(unt, iostat=ios, err=100)

    ! calculation of the energies
    ee = zero; el = zero; edy = zero

    do i = 1, n
        el = el + y(i)*y(i)
    end do
    el = el*fcnst/t0/two

    do i = 1, ne
        ee = ee + eu(i)*ocu(i) + ed(i)*ocd(i)
    end do

    uterm  = zero
    vterm1 = zero
    vterm2 = zero

    do i=1,n-1
        uterm  = uterm  + ru(i)*rd(i) - 0.25_wp
        vterm1 = vterm1 + ( ru(i) + rd(i) )*( ru(i+1) + rd(i+1) ) - one
        vterm2 = vterm2 + tu(i)*tu(i) + td(i)*td(i)
    end do
    uterm = uterm  + ru(n)*rd(n) - 0.25_wp
    vterm1= vterm1 + ( ru(n) + rd(n) )*( ru(1) + rd(1) ) - one
    vterm2= vterm2 + tu(n)*tu(n) + td(n)*td(n)
    ee = ee - u*uterm - v*vterm1 + v*vterm2
    et = ee + el + edy

    iterationscf = 0

    unt = free_unit()
    open(unt, file=trim(energies_fname), status='replace', action='write', iostat=ios, err=100, iomsg=io_err_msg)
    write(unt, 160, iostat=ios, err=100)
    write(unt,'(i7,4e18.10)', iostat=ios, err=100) iterationscf, ee, el, edy, et
    close(unt, iostat=ios, err=100)


    ! charge and spin density calculation

    chden = ru + rd
    sden  = half*( ru(i) - rd(i) )

    unt = free_unit()
    open(unt, file=trim(charge_density_fname), status='replace', action='write', iostat=ios, err=100, iomsg=io_err_msg)
    write(unt, 180, iostat=ios, err=100)
    write(unt, 190, iostat=ios, err=100) ( chden(i), i = 1, n)
    write(unt, 210, iostat=ios, err=100)
    write(unt, 190, iostat=ios, err=100) ( sden(i), i = 1, n)
    close(unt, iostat=ios, err=100)


    !calculating the temperature (t=mv2/kb) --> unit (m/kb)
    !set all the velocities equal to zero in the first step
    v2 = zero

    unt = free_unit()
    open(unt, file=trim(temperature_fname), status='replace', action='write', iostat=ios, err=100, iomsg=io_err_msg)
    write(unt, 300, iostat=ios, err=100)
    write(unt, 190, iostat=ios, err=100) ( v2(i), i = 1, n )
    write(unt, 310, iostat=ios, err=100)
    write(unt, 190, iostat=ios, err=100) ( yv(i), i = 1, n )
    close(unt, iostat=ios, err=100)

    return

    100 call io_error_check(unt, ios, io_err_msg)
    110   format (5x,i7,4(5x,d12.5))
    120   format (3x,'lattice displacement (bond) at 0-th step',/)
    140   format (3x,'electronic energy spectrum up at 0-th step',/)
    150   format (3x,'electronic energy spectrum down at 0-th step',/)
    160   format (/1x,'time step, electronic energy, lattice energy, kinetic energy, total energy',/)
    180   format (3x,'charge density at 0-th step',/)
    190   format (4d15.7)
    210   format (3x,'spin density at 0-th step',/)
    300   format (3x,'temperature v^2 at 0-th step',/)
    310   format (3x,'velocity of the sites at 0-th step',/)

end subroutine stationary

