! Solution of Differential form of LCR circuit using Runge-Kutta Order 4 method
! In an RLC circuit for given L, R & C and given boundary condition.
! ddi/dtt + R/L di/dt + 1/c i = E'/L, di/dt=z=f(t, i, z), dz/dt=ddi/dtt=E'-(R/L)z-(1/LC)i=g(t, i, z)  
! ddy/dxx + R/L dy/dx + 1/c y = E'/L, dy/dx=z=f(x, y, z), dz/dx=ddy/dxx=E'-(R/L)z-(1/LC)y=g(x, y, z) 
! Coded by: Ram Chandra Giri (Author) 
!Main Program
program RK4
    implicit none
    !Variables declaration
    real:: h, a, b, x0, y0, z0
    integer:: i, n
    !Functions declaration
    real:: f, g 
    !RK slope
    real:: k1, k2, k3, k4, l1, l2, l3, l4, k, l
    ! dynamic array declaration
    real, ALLOCATABLE :: x(:), y(:), z(:)
    
    !Varible Initialization
    a=0 !Initial time (sec)
    b=5 !Final Time    
    h=0.01 !Step-size (ms) depending on no. of points needed
    n=int((b-a)/h) !No. of steps

    !Array allocation
    ALLOCATE(x(n), y(n), z(n))

    !Creating a file to put data
    open(unit=10, file="data.dat")

    !Boundary condition (Initial Values)
    x0 = a !time (t)
    y0 = _ !current (i) or charge (q)
    z0 = _ !dummy variable, derivative of i/q

    !Array Initial values
    x(1) = x0 
    y(1) = y0
    z(1) = z0
    
    !RK4 part
    do i = 1, n, 1
        k1 = h* f(x(i), y(i), z(i))
        l1 = h* g(x(i), y(i), z(i))

        k2 = h* f(x(i)+h/2, y(i)+k1/2, z(i)+l1/2)
        l2 = h* g(x(i)+h/2, y(i)+k1/2, z(i)+l1/2)

        k3 = h* f(x(i)+h/2, y(i)+k2/2, z(i)+l2/2)
        l3 = h* g(x(i)+h/2, y(i)+k2/2, z(i)+l2/2)

        k4 = h* f(x(i)+h, y(i)+k3, z(i)+l3)
        l4 = h* g(x(i)+h, y(i)+k3, z(i)+l3)

        k = (k1 + 2*k2 + 2*k3 + k4)/6.0
        l = (l1 + 2*l2 + 2*l3 + l4)/6.0
        
        x(i+1) = x(i) + h
        y(i+1) = y(i) + k
        z(i+1) = z(i) + l        

        !Writing in a datafile: S. no., time, current/charge, derivative of current
        write(10,*) i, x(i), y(i), z(i)
    end do 
    !deallocate arrays
    DEALLOCATE(x, y, z)
    !Closing of File
    close(10)
end program RK4

!external function
real function f(x, y, z)
    real:: x, y, z
    f = z
end function
real function g(x, y, z)
    real:: x, y, z
    real, PARAMETER:: L=_ 
    real, PARAMETER:: R=_  
    real, PARAMETER:: C=_ 
    g = (E_dot/L)-(R/L)*z-(1/(L*C))*y
end function