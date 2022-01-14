subroutine calculate_energies(kernelfile, sdata, rdata, udata, energies, n) bind(C, name="calculate_energies_")
    ! =============================================================================
    ! Calculate the total three-body potential for all beads between three atoms
    ! Reads in the Jacobi point coordinates for the atomic triplet
    ! =============================================================================
    
    ! use statements and parameter declarations
    use, intrinsic :: iso_c_binding
    use RKHS
    implicit none                                            ! no implicit variables
    integer :: i                                             ! indices
    integer :: j
    integer, parameter :: dp = selected_real_kind(14)        ! double precision
    integer, parameter :: D  = 3                             ! number of Jacobi dimensions
    character(len=256) :: fkernelfile

    ! function argument declarations
    integer, intent(in)                      :: n            ! size of sdata, rdata, udata
    character(kind=c_char, len=1), dimension(*), intent(in) :: kernelfile
    real(kind=dp), dimension(n), intent(in)  :: sdata        ! 1st Jacobi coordinate
    real(kind=dp), dimension(n), intent(in)  :: rdata        ! 2nd Jacobi coordinate
    real(kind=dp), dimension(n), intent(in)  :: udata        ! 3rd Jacobi coordinate
    real(kind=dp), dimension(n), intent(out) :: energies     ! output three-body energies

    ! variables needed to call the RKHS kernel functions and store energies
    ! energies calculated from Jacobi points
    real(kind=dp) :: x(D)   ! coordinate array
    real(kind=dp) :: f      ! the energy
    type(kernel)  :: K      ! instance that evaluates an RKHS model

    fkernelfile = " "
    loop_string: do j=1,256
        if ( kernelfile(j) == c_null_char ) then
            exit loop_string
        else
            fkernelfile(j:j) = kernelfile(j)
        end if
    end do loop_string
    
    ! get the kernel
    call K%load_from_file(trim(fkernelfile))
    
    ! calculate the three-body energies, save them in the 'energies' array
    do i = 1,n
        x = (/ sdata(i), rdata(i), udata(i) /)
        call K%evaluate_fast(x, f)
        energies(i) = f
    end do
    
    ! free the kernel resources
    call K%free()
    
end subroutine
