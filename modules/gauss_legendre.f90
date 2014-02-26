module gauss_legendre
!Computes the Gauss-Legendre nodes/weights.

implicit none

    real(8), parameter, private :: pi = 3.141592653589793238462643d0

contains

subroutine gauss_legendre_quadrature(n,a,b,x,w)
!Computes the requested Gauss-Legendre rule.
!
!Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!Modified:
!
!    15 September 2010
!
!Author:
!
!    John Burkardt
!
!    Modified 5/14/2012 by Larz White, University of Idaho, Physics Department: Nuclear Theory Group.
!    I made it so when you call the program it outputs the arrays containing the nodes and weights;
!    instead of printing them to a file. I also got ride of all the extra stuff like outputting the
!    CPU time and string manipulation. Also, changed the function name from legendre_handle
!    to gauss_legendre_quadrature.
!
!    Modified 2/1/2013 by Larz White, University of Idaho, Physics Department: Nuclear Theory Group.
!    I went through every routine and made it more "Fortran 90" and less "Fortran 77". Also, the program
!    was using single precision intrinsic functions ( e.g. sin() ) instead of double precision intrinsic
!    functions ( e.g. dsin() ); I fixed that. Finally, I made it so it could be placed inside a module.

implicit none

    !n = Number of nodes/weights.
    !a = The lower integral bound.
    !b = The upper integral bound.
    !x = The nodes.
    !w = The weights.
    integer, intent(in) :: n
    real(8), intent(in) :: a,b
    real(8), dimension(n), intent(out) :: x,w

    call legendre_compute_glr(n,x,w)

    !Rescale the data.
    call rescale(n,a,b,x,w)

end subroutine

subroutine legendre_compute_glr(n,x,w)
!Legendre quadrature by the Glaser-Liu-Rokhlin method.
!
!Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!Modified:
!
!    20 October 2009
!
!Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
!    A fast algorithm for the calculation of the roots of special functions,
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.

implicit none

    !n = The order.
    !x,w = The abscissas/weights.
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: x,w

    real(8):: p,pp

    !Get the value and derivative of the N-th Legendre polynomial at 0.0
    call legendre_compute_glr0(n,p,pp)

    !If n is odd, then zero is a root.
    if( mod(n,2) == 1 ) then

        x((n+1)/2) = 0.0d0
        w((n+1)/2) = pp

        !If n is even, we have to compute a root.

    else

        call legendre_compute_glr2(p,n,x((n/2)+1),w((n/2)+1))

    end if

    !Get the complete set of roots and derivatives.
    call legendre_compute_glr1(n,x,w)

    !Compute w.
    w(1:n) = 2.0d0/(1.0d0 - x(1:n))/(1.0d0 + x(1:n))/w(1:n)/w(1:n)

    w(1:n) = ( 2.0d0*w(1:n) )/sum(w(1:n))

end subroutine

subroutine legendre_compute_glr0(n,p,pp)
!Gets a starting value for the fast algorithm.
!
!Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!Modified:
!
!    16 October 2009
!
!Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
!    A fast algorithm for the calculation of the roots of special functions,
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.

implicit none

    !n = The order of the Legendre polynomial.
    !p,pp = The value of the N-th Legendre polynomial and its derivative at 0.
    integer, intent(in) :: n
    real(8), intent(out) :: p,pp

    integer :: k
    real(8) :: pm1,pm2,ppm1,ppm2,rk

    !Compute coefficients of P_m(0), Pm'(0), m = 0,..,N
    pm2 = 0.0d0
    pm1 = 1.0d0
    ppm2 = 0.0d0
    ppm1 = 0.0d0

    do k = 0,n-1

        rk = dble(k)
        p = (-rk*pm2)/(rk + 1.0d0)
        pp = ( (2.0d0*rk + 1.0d0)*pm1 - rk*ppm2 )/(rk + 1.0d0)
        pm2 = pm1
        pm1 = p
        ppm2 = ppm1
        ppm1 = pp

    enddo

end subroutine

subroutine legendre_compute_glr1(n,x,w)
!Gets the complete set of Legendre points and weights.
!
!Discussion:
!
!    This routine requires that a starting estimate be provided for one
!    root and its derivative.  This information will be stored in entry
!    (n+1)/2 if n is odd, or n/2 if n is even, of x and w.
!
!Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!Modified:
!
!    15 November 2009
!
!Author:
!
!    Original C++ version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
!    A fast algorithm for the calculation of the roots of special functions,
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.

implicit none

    !n = The order of the Legendre polynomial.
    !x = On input, a starting value has been set in one entry.
    !    On output, the roots of the Legendre polynomial.
    !w = On input, a starting value has been set in one entry.
    !    On output, the derivatives of the Legendre polynomial at the zeros.
    integer, intent(in) :: n
    real(8), dimension(n), intent(inout) ::x,w

    !m = The number of terms in the Taylor expansion.
    integer, parameter :: m = 30
    integer :: j,k,l,n2,s
    real(8) :: dk,dn,h,xp
    real(8), dimension(m+1) :: up
    real(8), dimension(m+2) :: u


    if( mod(n,2) == 1 ) then

        n2 = ( (n - 1)/2 ) - 1
        s = 1

    else

        n2 = (n/2) - 1
        s = 0

    endif

    dn = dble(n)

    do j = n2+1,n-2

        xp = x(j+1)

        h = rk2_leg(pi/2.0d0,(-pi)/2.0d0,xp,n) - xp

        u(1) = 0.0d0
        u(2) = 0.0d0
        u(3) = w(j+1)

        up(1) = 0.0d0
        up(2) = u(3)

        do k = 0,m-2

            dk = dble(k)

            u(k+4) =(2.0d0*xp*(dk + 1.0d0)*u(k+3) + ( dk*(dk + 1.0d0) - dn*(dn + 1.0d0) )&
            &*u(k+2)/(dk + 1.0d0))/(1.0d0 - xp)/(1.0d0 + xp)/(dk + 2.0d0)

            up(k+3) = (dk + 2.0d0)*u(k+4)

        enddo

        do l = 0,4

            h = h - ( ts_mult(u,h,m)/ts_mult(up,h,m-1) )

        enddo

        x(j+2) = xp + h
        w(j+2) = ts_mult(up,h,m-1) 

    end do

    do k = 0,n2+s

        x(k+1) = -x(n-1-k+1)
        w(k+1) = w(n-1-k+1)

    enddo

end subroutine

subroutine legendre_compute_glr2(pn0,n,x1,d1)
!Finds the first real root.
!
!Discussion:
!
!    This routine is only called if n is even.
!
!Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!Modified:
!
!    20 October 2009
!
!Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
!    A fast algorithm for the calculation of the roots of special functions,
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.

implicit none

    !n = The order of the Legendre polynomial.
    !pn0 = The value of the N-th Legendre polynomial at 0.
    !x1 = The first real root.
    !d1 = The derivative at x1.
    integer, intent(in) :: n
    real(8), intent(in) :: pn0
    real(8), intent(out) :: x1,d1

    !m = The number of terms in the Taylor expansion.
    integer, parameter :: m = 30
    integer :: k,kk,l
    real(8) :: eps,rk,rn,scaling,step,theta
    real(8), dimension(m+1) :: u,up,x1k

    k = (n + 1)/2

    theta = pi*dble(4*k - 1)/dble(4*n + 2)

    x1 = (   1.0d0 - (dble(n - 1)/dble(8*n*n*n)) - (1.0d0/dble(384*n*n*n*n))&
    &*(  39.0d0 - ( 28.0d0/(dsin(theta)*dsin(theta)) )  )   )*dcos(theta)

    !Scaling.
    scaling = 1.0d0/x1

    !Recurrence relation for Legendre polynomials.
    u(1:m+1) = 0.0d0
    up(1:m+1) = 0.0d0

    rn = dble(n)

    u(1) = pn0

    do k = 0,m-2,2

        rk = dble(k)

        u(k+3) = ( rk*(rk + 1.0d0) - rn*(rn + 1.0d0) )*u(k+1)/(rk + 1.0d0)/(rk + 2.0d0)/scaling/scaling
        up(k+2) = (rk + 2.0d0)*u(k+3)*scaling

    enddo

    !Flip for more accuracy in inner product calculation
    u = u(m+1:1:-1)
    up = up(m+1:1:-1)

    x1k(1:m+1) = 1.0d0

    step = huge(0.0d0)
    l = 0

    !Newton iteration.
    eps = epsilon (0.0d0)

    do while ( eps < dabs(step) .and. l < 10 )

        l = l + 1
        step = dot_product(u(1:m+1),x1k(1:m+1))/dot_product(up(1:m+1),x1k(1:m+1))
        x1 = x1 - step
        x1k(1) = 1.0d0
        x1k(2) = scaling*x1

        do kk = 3,m+1
  
            x1k(kk) = x1k(kk-1)*scaling*x1

        enddo

        x1k(1:m+1) = x1k(m+1:1:-1)

    enddo

    d1 = dot_product(up(1:m+1),x1k(1:m+1))

end subroutine

subroutine rescale(n,a,b,x,w)
!Rescales a Legendre quadrature rule from (-1,1) to (a,b).
!
!Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!Modified:
!
!    18 October 2009
!
!Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
!    A fast algorithm for the calculation of the roots of special functions,
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.

implicit none

    !n = The order.
    !a,b = The endpoints of the new interval.
    !x = On input, the abscissas for (-1,1). On output, the abscissas for (a,b).
    !w = On input, the weights for (-1,1). On output, the weights for (a,b).
    integer, intent(in) :: n
    real(8), intent(in) :: a,b
    real(8), dimension(n), intent(inout) :: x,w

    x(1:n) = ( (a + b) + (b - a)*x(1:n) )/2.0d0
    w(1:n) = ( (b - a)*w(1:n) )/2.0d0

end subroutine

real(8) function rk2_leg(t1,t2,x,n)
!Advances the value of x(t) using a Runge-Kutta method. The output is the value of x at t2.
!
!Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!Modified:
!
!    15 November 2009
!
!Author:
!
!    Original C++ version by Nick Hale.
!    FORTRAN90 version by John Burkardt.

implicit none

    !n = The number of steps to take.
    !t1,t2 = The range of the integration interval.
    !x = The value of x at t1.
    integer, intent(in) :: n
    real(8), intent(in) :: t1,t2,x

    integer, parameter :: m = 10
    integer :: j
    real(8) :: f,h,k1,k2,snn1,t,x2

    x2 = x

    h = (t2 - t1)/dble(m)
    snn1 = dsqrt( dble(n*(n + 1)) )
    t = t1

    do j = 0,m-1

        f = (1.0d0 - x2)*(1.0d0 + x2)
        k1 = (-h*f)/( (snn1*dsqrt(f)) - (0.5d0*x2*dsin(2.0d0*t)) )
        x2 = x2 + k1

        t = t + h

        f = (1.0d0 - x2)*(1.0d0 + x2)
        k2 = (-h*f)/( (snn1*dsqrt(f)) - (0.5d0*x2*dsin(2.0d0*t)) )
        x2 = x2 + 0.5d0*(k2 - k1)

    enddo

    rk2_leg = x2

end function

real(8) function ts_mult(u,h,n)
!Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!Modified:
!
!    15 November 2009
!
!Author:
!
!    Original C++ version by Nick Hale.
!    FORTRAN90 version by John Burkardt.

implicit none

    integer, intent(in) :: n
    real(8), intent(in) :: h
    real(8), dimension(n+1), intent(in) :: u

    integer :: k
    real(8) :: hk,ts

    ts = 0.0d0
    hk = 1.0d0
    do k = 1,n

        ts = ts + ( u(k+1)*hk )
        hk = hk*h

    enddo

    ts_mult = ts

end function

end module gauss_legendre
