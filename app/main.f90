program fsieve
    use omp_lib
    implicit none

    ! Parameters for large range
    integer(kind=8), parameter :: lower_bound = 2_8  ! Starting exponent (adjust as needed)
    integer(kind=8), parameter :: upper_bound = 2147483647_8  ! Ending exponent (adjust as needed)
    integer(kind=8), parameter :: segment_size = 100000000  ! Size of each segment

    logical, allocatable :: is_prime(:)
    integer(kind=8) :: sieve_start, sieve_end, sqrt_upper, i, num_primes
    integer(kind=8), allocatable :: small_primes(:)
    integer(kind=8) :: small_primes_count
    real(kind=8) :: start_time, end_time
    integer :: unit_number

    ! Open the file for writing
    open (unit=unit_number, file='primes_si32.txt', status='replace')

    ! Start the timer
    start_time = omp_get_wtime()

    ! Precompute small primes up to sqrt(upper_bound)
    sqrt_upper = int(sqrt(real(upper_bound, kind=selected_real_kind(15, 307))), kind=8)
    call precompute_small_primes(sqrt_upper)

    ! Process large range in segments
    num_primes = 0
    do sieve_start = lower_bound, upper_bound, segment_size
        sieve_end = min(sieve_start + segment_size - 1, upper_bound)
        allocate (is_prime(sieve_end - sieve_start + 1))
        is_prime = .true.

        ! Sieve segment with small primes
        call sieve_segment(sieve_start, sieve_end, is_prime)

        ! Count and output primes in the segment
        do i = sieve_start, sieve_end
            if (is_prime(i - sieve_start + 1)) then
                num_primes = num_primes + 1
                write (unit_number, '(I0)') i  ! Write to file
                print *, i  ! Print to console
            end if
        end do

        deallocate (is_prime)
    end do

    ! End the timer
    end_time = omp_get_wtime()

    print *, "Total primes in range [", lower_bound, ",", upper_bound, "]: ", num_primes
    print *, "Execution time (seconds): ", end_time - start_time

    ! Close the file
    close (unit_number)

contains

    subroutine precompute_small_primes(limit)
        integer(kind=8), intent(in) :: limit
        logical, allocatable :: is_small_prime(:)
        integer(kind=8) :: i, j

        allocate (is_small_prime(limit))
        is_small_prime = .true.

        is_small_prime(1) = .false.  ! 1 is not a prime

        !$omp parallel do private(i, j) shared(is_small_prime)
        do i = 3, int(sqrt(real(limit))), 2  ! Skip even numbers, start from 3
            if (is_small_prime(i)) then
                do j = i*i, limit, i
                    is_small_prime(j) = .false.
                end do
            end if
        end do
        !$omp end parallel do

        ! Store the small primes in a global array
        allocate (small_primes(count(is_small_prime)))
        small_primes_count = 0

        !$omp parallel do private(i) shared(small_primes_count, small_primes)
        do i = 3, limit, 2  ! Skip even numbers, start from 3
            if (is_small_prime(i)) then
                !$omp critical
                small_primes_count = small_primes_count + 1
                small_primes(small_primes_count) = i
                !$omp end critical
            end if
        end do
        !$omp end parallel do

        ! Add 2 manually as it is the only even prime
        small_primes_count = small_primes_count + 1
        small_primes(small_primes_count) = 2

        deallocate (is_small_prime)
    end subroutine precompute_small_primes

    subroutine sieve_segment(start, end, is_prime)
        integer(kind=8), intent(in) :: start, end
        logical, intent(inout) :: is_prime(:)
        integer(kind=8) :: i, j, prime, offset, idx

        !$omp parallel do private(i, j, prime, offset, idx) schedule(dynamic)
        do i = 1, small_primes_count
            prime = small_primes(i)
            offset = max(prime*prime, ((start + prime - 1)/prime)*prime)

            ! Mark multiples of the prime as non-prime
            do j = offset, end, prime
                idx = j - start + 1
                is_prime(idx) = .false.
            end do
        end do
        !$omp end parallel do
    end subroutine sieve_segment

end program fsieve