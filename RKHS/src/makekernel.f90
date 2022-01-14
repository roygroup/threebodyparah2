subroutine parse_commandline(csvfile, kernelfile, power1, power2)
! =============================================================================
! Parses the command lines inputs
! csvfile    : the .csv file from which the kernel will be constructed
! kernelfile : the name of the .kernel file to be written
! power1     : the value of M (from 0 to 6 inclusive) for dimension 1
! power2     : the value of M (from 0 to 6 inclusive) for dimension 2
! =============================================================================
    implicit none
    character(len=256), intent(out) :: csvfile, kernelfile
    character(len=256)              :: strpower1, strpower2
    integer, intent(out)            :: power1, power2
    
    
    if (command_argument_count().ne.4) then
        write(*,*) "ERROR: requires 'csvfile', 'kernelfile', 'power1', 'power2'"
        stop
    endif
    
    ! read arguments from the command line
    call get_command_argument(1, csvfile)
    call get_command_argument(2, kernelfile)
    call get_command_argument(3, strpower1)
    call get_command_argument(4, strpower2)
    
    ! convert strpower to power
    read(strpower1, *) power1
    read(strpower2, *) power2
    
    if ( (power1.lt.0) .or. (power1.gt.6) ) then
        write(*,*) "ERROR: power1 must be between 0 and 6, inclusive"
        stop
    endif
    
    if ( (power2.lt.0) .or. (power2.gt.6) ) then
        write(*,*) "ERROR: power2 must be between 0 and 6, inclusive"
        stop
    endif
end subroutine

program example
use RKHS            ! This module needs to be used by your code
implicit none

type(kernel)       :: K                   ! The kernel type
character(len=256) :: csvfile, kernelfile ! names of files
integer            :: power1, power2      ! powers for the 1D kernels

call parse_commandline(csvfile, kernelfile, power1, power2) ! read command line inputs
call K%read_grid(trim(csvfile))                             ! read training data

!   choose one-dimensional kernel for dimension 1
if      (power1 == 0) then
    call K%k1d(1)%init(RECIPROCAL_POWER_N2_M0_KERNEL)
else if (power1 == 1) then
    call K%k1d(1)%init(RECIPROCAL_POWER_N2_M1_KERNEL)
else if (power1 == 2) then
    call K%k1d(1)%init(RECIPROCAL_POWER_N2_M2_KERNEL)
else if (power1 == 3) then
    call K%k1d(1)%init(RECIPROCAL_POWER_N2_M3_KERNEL)
else if (power1 == 4) then
    call K%k1d(1)%init(RECIPROCAL_POWER_N2_M4_KERNEL)
else if (power1 == 5) then
    call K%k1d(1)%init(RECIPROCAL_POWER_N2_M5_KERNEL)
else if (power1 == 6) then
    call K%k1d(1)%init(RECIPROCAL_POWER_N2_M6_KERNEL)
else
    write(*,*) "ERROR: something went wrong initializing the kernel for dimension 1"
    stop
endif

!   choose one-dimensional kernel for dimension 2
if      (power2 == 0) then
    call K%k1d(2)%init(RECIPROCAL_POWER_N2_M0_KERNEL)
else if (power2 == 1) then
    call K%k1d(2)%init(RECIPROCAL_POWER_N2_M1_KERNEL)
else if (power2 == 2) then
    call K%k1d(2)%init(RECIPROCAL_POWER_N2_M2_KERNEL)
else if (power2 == 3) then
    call K%k1d(2)%init(RECIPROCAL_POWER_N2_M3_KERNEL)
else if (power2 == 4) then
    call K%k1d(2)%init(RECIPROCAL_POWER_N2_M4_KERNEL)
else if (power2 == 5) then
    call K%k1d(2)%init(RECIPROCAL_POWER_N2_M5_KERNEL)
else if (power2 == 6) then
    call K%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)
else
    write(*,*) "ERROR: something went wrong initializing the kernel for dimension 2"
    stop
endif

! choose one-dimensional kernel for dimension 3
call K%k1d(3)%init(TAYLOR_SPLINE_N2_KERNEL)

call K%calculate_coefficients_fast()   ! calculate the coefficients (uses Tensor product)
call K%calculate_sums()                ! this call will calculate the complete lookup table
call K%save_to_file(trim(kernelfile))  ! save the kernel to a binary file
call K%free()                          ! free the kernel memory

end program example
