!-----------------------------------------------------------------------
!> @file    io_error_checking.f90
!! @brief   IO Error checking module.
!! @details This module checks for IO errors when operating with arrays
!!          and files.
!! @author  bgeneto
!! @date    20150623
!-----------------------------------------------------------------------
module io_error_checking
    implicit none
    private
    public  :: io_error_check, ios, io_err_msg, unt, free_unit

    ! some useful global variables when dealing with IO operations
    integer   :: ios      !< IO status code
    integer   :: unt      !< a file unit number
    character(len=255) :: io_err_msg = ''   !< IO operation error message

    !> Checks if IO errors occurred while working with arrays or files.
    interface io_error_check
        module procedure io_check_file
        module procedure io_check_array
    end interface io_error_check

contains

    !> @brief     Checks if an error occurred while reading, writing, opening or
    !!            closing external files.
    !! @details   This routine tries to print (in the standard output unit) some
    !!            useful information relative to the triggered error, like the
    !!            error code, the error message, and the corresponding file name.
    !! @warning   This procedure immediately stops program execution.
    !! @param[in] unt File unit number.
    !! @param[in] ios IO status code.
    !! @param[in] io_err_msg IO error message.
    subroutine io_check_file(unt, ios, io_err_msg)
        integer, intent(in)          :: unt, ios
        character(len=*), optional   :: io_err_msg
        character(len=5)             :: sunt, sios
        character(len=90)            :: file_name
        character(len=*), parameter  :: cont_str = '          '
        logical                      :: is_named

        if (ios == 0) return ! no error. Nothing to do.

        if ( .not. present(io_err_msg) ) io_err_msg = ''
        write (sios,'(i0)') ios ! using internal write to convert from integer to string
        inquire (unit = unt, named = is_named)
        if (is_named) inquire (unit=unt, name=file_name)
        write(*,'(a)', advance='no') '***ERROR: '
        if ( len_trim(io_err_msg) > 0 ) write (*,'(a)') trim(io_err_msg)
        select case (ios)
        case(1:)
            if (is_named) then
                write (*,'(a)') "Error while reading from file named '" // trim(file_name) // "'."
                write (*,'(a)') cont_str // "Aborting program execution with read error no. " // trim(sios)
                stop
            else
                ! use internal write to convert from integer to string
                write (sunt,'(i0)') unt
                write (*,200) cont_str, trim(sunt)
                write (*,300) cont_str, trim(sios)
                200 format (a, 'Unable to open file unit = ', a, '. Please try again!')
                300 format (a,'Aborting program execution with open error no. ',a,'.')
                stop
            endif
        case(-1)
            write (*,'(a/a)')   "Unable to read from file named '" // trim(file_name) // "'.", &
                                cont_str // 'Possible cause: end of file (EOF) during read!'
            write (*,'(a)')     cont_str // 'Aborting program execution...'
            stop
        case(-2)
            write (*,'(a/a)')  "Unable to read from file named '" // trim(file_name) // "'.", &
                               cont_str // 'Possible cause: end of record during read!'
            write (*,'(a)')    cont_str // 'Aborting program execution...'
            stop
        end select
    end subroutine io_check_file


    !> Check dynamic allocation status.
    !! @param[in] ios IO status as returned by the (de)allocate command.
    !! @param[in] str (Optional) Variable name.
    !! @warning This subroutine stops program execution in case of
    !!          failure in allocation operation.
    subroutine io_check_array(ios, str)
        integer, intent(in)    :: ios
        character(len=*), optional :: str

        if ( .not. present(str) ) str = ''

        if (ios == 0) return ! no error, nothing to do.

        print *, '***FATAL: allocation error! Array ', str, ' cannot be allocated.'
        stop

    end subroutine io_check_array

    !> Finds the first available (not used) unit number from \c min value
    !! to a \c max value (if given).
    !! @details   Can be used with IO file operations where an unit number
    !!            is required. This way one can avoid the usage of an explicit
    !!            unit number while performing operations like open, read/write and close.
    !! @param[in] unt_min (Optional) Starts searching for a free unit number
    !!            beginning with this value (defaults to 20).
    !! @param[in] unt_max (Optional) Stops searching for free unit numbers at this value
    !!            (defaults to 90).
    !! @warning   Program execution is aborted if a (free) unit number cannot be found.
    !! @returns   An integer: the first (free) unit number found.
    function free_unit(unt_min, unt_max) result(unt)
        integer, optional :: unt_min, unt_max
        integer           :: unt, ini=20, fin=90
        logical           :: is_open

        if ( present(unt_min) ) ini = unt_min
        if ( present(unt_max) ) fin = unt_max

        search: do unt = ini, fin
            inquire (unit=unt, opened=is_open)
            if (.not. is_open) return
        end do search

        print'(2(a,i0))', "*** ERROR: Cannot find a (free) unit number. Searched from ", unt_min, " to ", unt_max
        stop
    end function free_unit

end module io_error_checking
