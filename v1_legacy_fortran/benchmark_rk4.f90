! v1_legacy_fortran/benchmark_rk4.f90
program RK4_Benchmark
    implicit none
    real :: h, a, b, x0, y0, z0
    integer :: i, n
    real :: f, g
    real :: k1, k2, k3, k4, l1, l2, l3, l4, k, l
    real, ALLOCATABLE :: x(:), y(:), z(:)
    real :: start_time, end_time

    ! 1 Million Step Benchmark
    a = 0.0
    b = 10.0
    n = 1000000 
    h = (b - a) / real(n)

    ALLOCATE(x(n), y(n), z(n))

    ! Initial Conditions (Using your calculated superposition)
    x0 = a
    y0 = 0.0
    z0 = 0.2

    x(1) = x0
    y(1) = y0
    z(1) = z0

    print *, "Starting Fortran RK4 Benchmark (1,000,000 steps)..."
    CALL CPU_TIME(start_time)

    ! Core Math Loop
    do i = 1, n - 1
        k1 = h * f(x(i), y(i), z(i))
        l1 = h * g(x(i), y(i), z(i))

        k2 = h * f(x(i)+h/2.0, y(i)+k1/2.0, z(i)+l1/2.0)
        l2 = h * g(x(i)+h/2.0, y(i)+k1/2.0, z(i)+l1/2.0)

        k3 = h * f(x(i)+h/2.0, y(i)+k2/2.0, z(i)+l2/2.0)
        l3 = h * g(x(i)+h/2.0, y(i)+k2/2.0, z(i)+l2/2.0)

        k4 = h * f(x(i)+h, y(i)+k3, z(i)+l3)
        l4 = h * g(x(i)+h, y(i)+k3, z(i)+l3)

        k = (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
        l = (l1 + 2.0*l2 + 2.0*l3 + l4)/6.0
        
        x(i+1) = x(i) + h
        y(i+1) = y(i) + k
        z(i+1) = z(i) + l
    end do

    CALL CPU_TIME(end_time)
    
    print *, "Fortran Benchmark Complete!"
    print *, "Execution Time: ", end_time - start_time, " seconds"

    DEALLOCATE(x, y, z)
end program RK4_Benchmark

real function f(x, y, z)
    real :: x, y, z
    f = z
end function

real function g(x, y, z)
    real :: x, y, z
    real, PARAMETER :: L = 1.0, R = 10.0, C = 0.0025, w = 2.5
    real :: E_dot
    E_dot = 0.2 * sin(w * x)
    g = (E_dot/L) - (R/L)*z - (1.0/(L*C))*y
end function