!-----------------------------------------------------------------------
!> @file    utils.f90
!! @brief   General purpose Fortran auxiliary functions and subroutines.
!! @details This module also contains useful mathematical and physical
!!          constants and other interesting parameters.
!! @author  bgeneto
!! @date    20150622
!-----------------------------------------------------------------------
module utils
    ! use iso kind (type) parameter values (this a Fortran 2008 feature!)
    use iso_fortran_env
    ! force declaration of all variables
    implicit none
    ! default access for this module is public
    public
    ! only identifiers below are private
    private :: int8_num_digits, int16_num_digits, int32_num_digits,    &
               int64_num_digits, real32_num_digits, real64_num_digits, &
               real128_num_digits


    interface num_digits
        module procedure int8_num_digits
        module procedure int16_num_digits
        module procedure int32_num_digits
        module procedure int64_num_digits
        module procedure real32_num_digits
        module procedure real64_num_digits
        module procedure real128_num_digits
    end interface num_digits

    interface string
        module procedure int8_to_str
        module procedure int16_to_str
        module procedure int32_to_str
        module procedure int64_to_str
        module procedure real32_to_str
        module procedure real64_to_str
        module procedure real128_to_str
    end interface string

contains

    !> Returns the number of digits of an \c integer (8-bit) number.
    !! @param[in] num The \c integer number.
    !! @returns   Number of digits contained in the given \c integer number \c num.
    pure function int8_num_digits(num) result(nchar)
        integer(int8), intent(in) :: num
        integer(int8) :: nchar
        character(40) :: tmp

        write(tmp, '(i0)') num
        nchar = len_trim(tmp, int8)
    end function int8_num_digits

    !> Returns the number of digits of an \c integer (16-bit) number.
    !! @param[in] num The \c integer number.
    !! @returns   Number of digits contained in the given \c integer number \c num.
    pure function int16_num_digits(num) result(nchar)
        integer(int16), intent(in) :: num
        integer(int16) :: nchar
        character(40)  :: tmp

        write(tmp, '(i0)') num
        nchar = len_trim(tmp, int16)
    end function int16_num_digits

    !> Returns the number of digits of an \c integer (32-bit) number.
    !! @param[in] num The \c integer number.
    !! @returns   Number of digits contained in the given \c integer number \c num.
    pure function int32_num_digits(num) result(nchar)
        integer(int32), intent(in) :: num
        integer(int32) :: nchar
        character(40)  :: tmp

        write(tmp, '(i0)') num
        nchar = len_trim(tmp, int32)
    end function int32_num_digits

    !> Returns the number of digits of an \c integer (64-bit) number.
    !! @param[in] num The \c integer number.
    !! @returns   Number of digits contained in the given \c integer number \c num.
    pure function int64_num_digits(num) result(nchar)
        integer(int64), intent(in) :: num
        integer(int64) :: nchar
        character(40)  :: tmp

        write(tmp, '(i0)') num
        nchar = len_trim(tmp, int64)
    end function int64_num_digits

    !> Returns the number of digits (including the decimal comma or dot)
    !! of a \c real (single precision) number.
    !! @param[in] num The \c real number.
    !! @returns   Number of digits contained in the given \c real number \c num.
    pure function real32_num_digits(num) result(nchar)
        real(real32), intent(in) :: num
        integer              :: nchar
        character(50)        :: tmp

        write(tmp, *) num
        nchar = len_trim(tmp)
    end function real32_num_digits

    !> Returns the number of digits (including the decimal comma or dot)
    !! of a \c real (double precision) number.
    !! @param[in] num The \c real number.
    !! @returns   Number of digits contained in the given \c real number \c num.
    pure function real64_num_digits(num) result(nchar)
        use iso_fortran_env
        real(real64), intent(in) :: num
        integer              :: nchar
        character(50)        :: tmp

        write(tmp, *) num
        nchar = len_trim(tmp)
    end function real64_num_digits

    !> Returns the number of digits (including the decimal comma or dot)
    !! of a \c real (quad precision) number.
    !! @param[in] num The \c real number.
    !! @returns   Number of digits contained in the given \c real number \c num.
    pure function real128_num_digits(num) result(nchar)
        use iso_fortran_env
        real(real128), intent(in) :: num
        integer              :: nchar
        character(50)        :: tmp

        write(tmp, *) num
        nchar = len_trim(tmp)
    end function real128_num_digits

    !> Converts an \c integer (8-bit) number to string.
    !! @param[in] num The \c integer number.
    !! @returns   The string representing the given \c integer number \c num.
    pure function int8_to_str(num) result(str)
        integer(int8), intent(in) :: num
        character( int8_num_digits(num) ) :: str

        write(str, '(i0)') num
    end function int8_to_str

    !> Converts an \c integer (16-bit) number to string.
    !! @param[in] num The \c integer number.
    !! @returns   The string representing the given \c integer number \c num.
    pure function int16_to_str(num) result(str)
        integer(int16), intent(in) :: num
        character( int16_num_digits(num) ) :: str

        write(str, '(i0)') num
    end function int16_to_str

    !> Converts an \c integer (32-bit) number to string.
    !! @param[in] num The \c integer number.
    !! @returns   The string representing the given \c integer number \c num.
    pure function int32_to_str(num) result(str)
        integer(int32), intent(in) :: num
        character( int32_num_digits(num) ) :: str

        write(str, '(i0)') num
    end function int32_to_str

    !> Converts an \c integer (64-bit) number to string.
    !! @param[in] num The \c integer number.
    !! @returns   The string representing the given \c integer number \c num.
    pure function int64_to_str(num) result(str)
        integer(int64), intent(in) :: num
        character( int64_num_digits(num) ) :: str

        write(str, '(i0)') num
    end function int64_to_str

    !> Converts a \c real (32-bit) number to string.
    !! @param[in] num The \c real number.
    !! @returns   The string representing the given \c real number \c num.
    pure function real32_to_str(num) result(str)
        real(real32), intent(in) :: num
        character( real32_num_digits(num) ) :: str

        write(str, *) num
    end function real32_to_str

    !> Converts a \c real (64-bit) number to string.
    !! @param[in] num The \c real number.
    !! @returns   The string representing the given \c real number \c num.
    pure function real64_to_str(num) result(str)
        real(real64), intent(in) :: num
        character( real64_num_digits(num) ) :: str

        write(str, *) num
    end function real64_to_str

    !> Converts a \c real (128-bit) number to string.
    !! @param[in] num The \c real number.
    !! @returns   The string representing the given \c real number \c num.
    pure function real128_to_str(num) result(str)
        real(real128), intent(in) :: num
        character( real128_num_digits(num) ) :: str

        write(str, *) num
    end function real128_to_str

    !> Returns the length the string \c filename without its last extension.
    !! @param[in] filename The file name.
    !! @returns   The length of the string \c filename without its last extension.
    pure function len_noext(filename) result(length)
        character(len=*), intent(in) :: filename
        integer :: length

        length = index(trim(filename), '.', back=.true.) - 1
        if (length < 1) length = len_trim(filename)
    end function len_noext

    !> Returns the filename string length without its path.
    !! @param[in] filename - nome do arquivo
    !! @returns   An \c integer: the filename string length without path.
    pure function len_nopath(filename) result(length)
        character(len=*), intent(in) :: filename
        integer :: length, fin

        fin = len_trim(filename)
        length = fin - index(trim(filename), '/', back=.true.)
        if (length == fin) length = fin - index(trim(filename), '\', back=.true.)
        if (length == fin) length = fin
    end function len_nopath

    !> Removes the filename extension or its last extension case more than
    !! one extension exists.
    !! @param[in] filename The file name.
    !! @returns   The string \c filename without its last extension.
    function noext(filename) result(str)
        character(len=*), intent(in)     :: filename
        character(len_noext(filename))  :: str
        integer :: n

        n   = len_noext(filename)
        str = filename(1:n)
    end function noext

    !> Removes the path from a file name.
    !! @param[in] filename The file name.
    !! @returns The string \c filename without its full path.
    function basename(filename) result(str)
        character(len = *), intent(in)     :: filename
        character(len_nopath(filename)) :: str
        integer :: ini, fin

        fin = len_trim(filename)
        ini = index(trim(filename), '/', back=.true.) + 1
        if (ini < 1) ini = index(trim(filename), '\', back=.true.) + 1
        if (ini < 1) ini = 1

        str = trim(filename(ini:fin))
    end function basename


    !> Plots the data stored in datafile and generates an output file in
    !! ps (postscript) format.
    !! @warning This routines requires the gnuplot software in your path
    !!          environment variable to work properly.
    !! @param   datafile Data file.
    !! @param   xlabel  (Optional) \a x axis label.
    !! @param   ylabel  (Optional) \a y axis label.
    !! @param   title   (Optional) A title for the plotted graphics.
    !! @param   offsets (Optional) Internal graphic margins.
    !! @param   nplots  (Optional) Number of curves to plot.
    !! @param   plotdir (Optional) Output directory to store the graphics in ps format.
    subroutine gplot(datafile,xlabel,ylabel,title,offsets,nplots,plotdir)
        use io_error_checking
        character(len=*), intent(in)           :: datafile
        integer,          intent(in), optional :: nplots
        character(len=*), intent(in), optional :: xlabel,   &
                                                  ylabel,   &
                                                  title,    &
                                                  offsets,  &
                                                  plotdir

        character(len_noext(trim(datafile))+3) :: cmdfile      !< gnuplot command file with .gp extension
        logical                                :: file_exists  !< used to check if datafile exists
        character(len=*), parameter            :: err_str = '***ERROR: '

        ! check if datafile exists
        file_exists = .false.
        inquire (file=trim(datafile), exist=file_exists)
        if ( .not. file_exists ) then
            write(*,'(a/a)')  err_str // 'GNUPLOT datafile ' //    &
                              trim(datafile) // ' does not exist!',      &
                             'Plot canceled!'
            return
        end if

        ! creates gnuplot command file if it does not already exist
        file_exists = .false.
        cmdfile = noext(trim(datafile)) // ".gp"  !< name of gnuplot command file
        inquire (file=trim(cmdfile), exist=file_exists)
        if ( .not. file_exists ) then
            unt = free_unit()
            open (unt, file=trim(cmdfile), action='write', status='new', iostat=ios, iomsg=io_err_msg)
            call io_error_check(unt, ios, io_err_msg)
            ! build command file
            write(unt,'(a)') 'set encoding iso_8859_1'     !< codepage de código
            write(unt,'(a)') 'set mxtics; set mytics'      !< x & y minor ticks
            write(unt,'(a)') 'unset key'                   !< removes legend
            write(unt,'(a)') 'set style data lines'        !< graph using lines
            write(unt,'(a)') 'set terminal pdf enhanced color'
            if ( present(plotdir) .and. len_trim(plotdir) > 0 ) then
                write(unt,'(a)') 'set output "'//plotdir//'/'//noext(basename(datafile))//'.pdf"'
            else
                write(unt,'(a)') 'set output "'//noext(trim(datafile))//'.pdf"'
            end if
            if ( present(title) )   &
                write(unt,'(a)') 'set title "' // trim(title) // '" font "Sans,24" enhanced'
            if ( present(xlabel) )  &
                write(unt,'(a)') 'set xlabel "' // trim(xlabel) // '" font "Sans,14" enhanced'
            if ( present(ylabel) )  &
                write(unt,'(a)') 'set ylabel "' // trim(ylabel) // '" font "Sans,14" enhanced'
            if ( present(offsets) ) &
                write(unt,'(2a)')'set offsets ', offsets
            if ( present(nplots) ) then
                select case (nplots)
                case (2)
                    write(unt,'(a)') 'plot "' // trim(datafile) // '" u 1:2,"' // trim(datafile) // '" u 1:3'
                case (3)
                    write(unt,'(a)') 'plot "' // trim(datafile) // '" u 1:2,"' // trim(datafile) // '" u 1:3,"' &
                                     // trim(datafile) // '" u 1:4'
                case default
                    write(unt,'(a)') 'plot "' // trim(datafile) // '"'
                end select
            else
                write(unt,'(a)') 'plot "' // trim(datafile) // '"'
            end if

            ! close the command file
            close(unt)
        end if

        ! runs gnuplot command
        call execute_command_line('gnuplot ' // trim(cmdfile), wait=.true., cmdstat=ios, cmdmsg=io_err_msg)

        ! checks execution status
        if ( ios /= 0 )  then
            write(*,'(a/a)') err_str // "The command 'gnuplot "// trim(cmdfile) // "' returned an error.", io_err_msg
        end if

    end subroutine gplot

    !> Skip comment lines at the beginning of an input file.
    !! @param[in] unt Unit number of the input file.
    !! @param[in] cchar (Optional) Comment character. (Default value is '!').
    !! @warning This routine \b stops the program execution in case of errors.
    subroutine skip_comments (unt, cchar)
        use io_error_checking, only: ios, io_error_check
        integer, intent(in) :: unt
        character, optional :: cchar
        character           :: ch, comment  !< read character and comment character, respectively.
        character(3)        :: rflag        !< read flag.
        character(len=90)   :: file_name    !< file name.
        logical             :: is_named, &  !< is named flag.
                               is_opened    !< is opened flag.

        comment = '!'
        if ( present(cchar) ) comment = cchar
        ! checks if the file is opened for reading
        inquire (unit = unt, opened = is_opened)
        if (is_opened) then
            inquire (unit = unt, named = is_named, name = file_name, read = rflag)
        else
            write(*,'(/a)')  '***ERROR while trying to skip comment lines ***'
            write(*,'(/a)')  'Possible cause: the input file is not opened!'
            write(*,'(/a/)') 'Aborting program execution...'
            stop
        end if

        if (rflag /= 'YES' .and. is_named) then
            write(*,'(/a)')  '***ERROR while trying to skip input file comment lines ***'
            write(*,'(/a,a)')'File name: ', trim(file_name)
            write(*,'(/a)')  'Possible cause: the file must be opened with read (or rw) action!'
            write(*,'(/a/)') 'Aborting program execution...'
            stop
        end if

        ! this actually skips lines beginning with a comment char.
        skip: do
          read(unt, '(a1)', iostat=ios, err=100) ch
          if (ch /= comment .and. ch /= '') exit skip
        end do skip
        ! rewinds the file pointer to the beginning of the last record to reread.
        backspace (unt, iostat=ios, err=100)
        100 call io_error_check(unt, ios)
    end subroutine skip_comments

end module utils
