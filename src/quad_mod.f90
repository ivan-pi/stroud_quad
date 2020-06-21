!>
!   Gaussian quadrature routines
!
!### Author
! * Stroud, A. H. and Secrest, Don -- original authors (1966)
! * Pribec, Ivan -- refactoring and modernization (2019)

module quad_mod

    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    private

    public :: wp, pi
    public :: jacobi, hermit, laguer

    real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
    real(wp), parameter :: epsm = epsilon(1.0_wp)
    integer, parameter :: max_iter = 10

contains

!>
!   Calculates the zeros \(x_i\) of the `nn`-th order Jacobi
!   polynomial \(P^{(\alfa,\beta)}_n(x)\) for the segment \((-1,1)\)-
!
!   The largest zero will be stored in `x(1)`. Also calculates the
!   corresponding coefficients \(A_i\) of the `nn`-th order Gauss-Jacobi
!   quadrature formula of degree `2*nn-1`. This subroutine must be given
!   the coefficients
!   \[b_n = \frac{(\alfa + \beta)(\beta - \alfa)}{(\alfa+\beta+2n)(\alfa+\beta+2n-2)}\]
!   \[c_n = \frac{4(n-1)(\alfa+n-1)(\beta+n-1)(\alfa+\beta+n-1)}{(\alfa+\beta+2n-1)(\alfa+\beta+2n-2)^2 (\alfa+\beta+2n-3)}\]
!   in the recursion relation
!   \[P_n = (x - b_n)P_{n-1} - c_n P_{n-1}\]
!   for all \(n\) less than or equal to the highest degree `nn`.
!
    subroutine jacobi(nn,x,a,alf,bta,b,c,eps,csx,csa,tsx,tsa)
        integer, intent(in) :: nn           !! Degree of polynomial
        real(wp), intent(out) :: x(*)       !! Zeros of the polynomial
        real(wp), intent(out) :: a(*)       !! Quadrature weights
        real(wp), intent(in) :: alf         !! Exponent $\alfa$ of the Jacobi polynomial
        real(wp), intent(in) :: bta         !! Exponent $\beta$ of the Jacobi polynomial
        real(wp), intent(in) :: b(*), c(*)  !! Coefficients of the recurrence relation
        real(wp), intent(in), optional :: eps         !! Desired accuracy
        real(wp), intent(out) :: csx        !! Calculated sum of zeros \(x_i\)
        real(wp), intent(out) :: csa        !! Calculated sum of quadrature weights \(A_i\)
        real(wp), intent(out) :: tsx        !! True sum of zeros \(x_i\)
        real(wp), intent(out) :: tsa        !! True sum of quadrature weights \(A_i\)

        integer :: i, j
        real(wp) :: fn, beta, cc
        real(wp) :: an, bn, r1, r2, r3, ratio, xt, dpn, pn1
        real(wp) :: eps_

        eps_ = epsm
        if (present(eps)) eps_ = eps

        fn = nn
        csx = 0.0_wp
        csa = 0.0_wp
        beta = exp(log_gamma(alf+1.0_wp) + log_gamma(bta+1.0_wp) - log_gamma(alf+bta+2.0_wp))
        cc = 2.0_wp**(alf+bta+1.0_wp)*beta
        tsx = fn*(bta - alf)/(alf + bta + 2.0_wp*fn)
        tsa = cc
        do j=2,nn
          cc = cc*c(j)
        end do

        do i = 1, nn
            if (i == 1) then ! largest zero
                an = alf/fn
                bn = bta/fn
                r1 = (1.0_wp+alf)*(2.78_wp/(4.0_wp+fn*fn)+0.768_wp*an/fn)
                r2 = 1.0_wp+1.48_wp*an+0.96_wp*bn+0.452_wp*an*an+0.83_wp*an*bn
                xt = 1.0_wp - r1/r2
            else if (i == 2) then ! second zero
                r1 = (4.1_wp+alf)/((1.0_wp+alf)*(1.0_wp+0.156_wp*alf))
                r2 = 1.0_wp + 0.06_wp*(fn-8.0_wp)*(1.0_wp+0.12_wp*alf)/fn
                r3 = 1.0_wp + 0.012_wp*bta*(1.0_wp+0.25_wp*abs(alf))/fn
                ratio = r1*r2*r3
                xt = xt - ratio*(1.0_wp-xt)
            else if (i == 3) then
                r1 = (1.67_wp+0.28_wp*alf)/(1.0_wp+0.37_wp*alf)
                r2 = 1.0_wp + 0.22_wp*(fn-8.0_wp)/fn
                r3 = 1.0_wp + 8.0_wp*bta/((6.28_wp+bta)*fn*fn)
                ratio = r1*r2*r3
                xt = xt - ratio*(x(1) - xt)
            else if (i == nn - 1) then ! second last zero
                r1 = (1.0_wp + 0.235_wp*bta)/(0.766_wp+0.119_wp*bta)
                r2 = 1.0_wp/(1.0_wp + 0.639_wp*(fn-4.0_wp)/(1.0_wp+0.71_wp*(fn-4.0_wp)))
                r3 = 1.0_wp/(1.0_wp + 20.0_wp*alf/((7.5_wp+alf)*fn*fn) )
                ratio = r1*r2*r3
                xt = xt + ratio*(xt-x(i-2))
            else if (i == nn) then ! last zero
                r1 = (1.0_wp+0.37_wp*bta)/(1.67_wp+0.28_wp*bta)
                r2 = 1.0_wp/(1.0_wp + 0.22_wp*(fn-8.0_wp)/fn )
                r3 = 1.0_wp/(1.0_wp + 8.0_wp*alf/((6.28_wp+alf)*fn*fn) )
                ratio = r1*r2*r3
                xt = xt + ratio*(xt-x(i-2))
            else !middle zeros
                xt = 3.0_wp*x(i-1) - 3.0_wp*x(i-2) + x(i-3)
            end if

            call root(xt,nn,alf,bta,dpn,pn1,b,c,eps_)
            x(i) = xt
            a(i) = cc/(dpn*pn1)
            WRITE(*,'(2f6.2,2i3,2(1x,e16.10),2x,2(1x,e16.10))') alf,bta,nn,i,xt,a(i)
            csx = csx + xt
            csa = csa + a(i)
        end do
        write(*,'(2f6.2,2i3,2(1x,e16.10),2x,2(1x,e16.10))') alf,bta,nn,i,csx,csa,tsx,tsa
    end subroutine

!>
!   Improves the approximate root x of the Jacobi polynomial
!
!   In addition we also obtain the derivative of $P(n)$ at $x$
!   and the value of $P(n-1)$ at x
!
    subroutine root(x,nn,alf,bta,dpn,pn1,b,c,eps)
        real(wp), intent(inout) :: x        !! Approximate root value
        integer, intent(in) :: nn           !! Degree of polynomial
        real(wp), intent(in) :: alf         !! Exponent $\alfa$ of the Jacobi polynomial
        real(wp), intent(in) :: bta         !! Exponent $\beta$ of the Jacobi polynomial
        real(wp), intent(out) :: dpn        !! Derivative of P(n) at x
        real(wp), intent(out) :: pn1        !! Value of P(n-l) at x
        real(wp), intent(in) :: b(*), c(*)  !! Coefficients of the recurrence relation
        real(wp), intent(in) :: eps         !! Desired accuracy

        integer :: iter
        real(wp) :: d, p, dp

        do iter = 1, max_iter
            call recur(p,dp,pn1,x,nn,alf,bta,b,c)
            d = p/dp
            x = x - d
            if (abs(d) < eps) exit
        end do
        dpn = dp
    end subroutine

!>
!   Recurrence relation of the Jacobi polynomial
!
    subroutine recur(pn,dpn,pn1,x,nn,alf,bta,b,c)
        real(wp), intent(out) :: pn
        real(wp), intent(out) :: dpn
        real(wp), intent(out) :: pn1
        real(wp), intent(in) :: x
        integer, intent(in) :: nn
        real(wp), intent(in) :: alf
        real(wp), intent(in) :: bta
        real(wp), intent(in) :: b(*)
        real(wp), intent(in) :: c(*)

        integer :: j
        real(wp) :: p1, p, dp1, dp, q, dq

        p1  = 1.0_wp
        p   = x + (alf-bta)/(alf+bta+2.0_wp)
        dp1 = 0.0_wp
        dp  = 1.0_wp
        do j = 2, nn
            q = (x - b(j))*p - c(j)*p1
            dq = (x - b(j))*dp + p - c(j)*dp1
            p1 = p
            p = q
            dp1 = dp
            dp = dq
        end do
        pn = p
        dpn = dp
        pn1 = p1
    end subroutine


!>
!   Calculates the zeros \(x_i\) of the `nn`-th order
!   Laguerre polynomial \(L^{\alpha}_n(x)\) for the segment \((0,\infty)\).
!
!   The smallest zero will be stored in `x(1)`. Also
!   calculates the corresponding coefficients \(A_i\))
!   of the `nn`-th order Laguerre quadrature formula
!   of degree `2*nn-1`:
!   \[\int^{+\infty}_0 x^\alpha e^{-x} f(x) dx \simeq \sum_{i=1}^nn A_i f(x_i)\]
!
!   This subroutine must be given the coefficients
!   \[b_n = \alpha + 2n - 1\]
!   \[c_n = (n-1)(\slpha + n - 1\]
!           C(N) = (N-1)(ALF + N - 1)
!   in the recursion relation
!   \[P_n = (x-b_n)P_{n-1} - c_n P_{n-1}\]
!   for all \(n\) less than or equal to the highest degree `nn`.
!
    subroutine laguer(nn,x,a,alf,b,c,eps,csx,csa,tsx,tsa)
        integer, intent(in) :: nn           !! Degree of polynomial
        real(wp), intent(out) :: x(nn)      !! Zeros of the polynomial
        real(wp), intent(out) :: a(nn)      !! Quadrature weights
        real(wp), intent(in) :: alf         !! Exponent $\alfa$ of the Laguerre polynomial
        real(wp), intent(in) :: b(nn),c(nn) !! !! Coefficients of the recurrence relation
        real(wp), intent(in), optional :: eps         !! Desired accuracy
        real(wp), intent(out) :: csx        !! Calculated sum of zeros \(x_i\)
        real(wp), intent(out) :: csa        !! Calculated sum of quadrature weights \(A_i\)
        real(wp), intent(out) :: tsx        !! True sum of zeros \(x_i\)
        real(wp), intent(out) :: tsa        !! True sum of quadrature weights \(A_i\)

        real(wp) :: fn, cc, xt, fi, r1, r2, ratio, dpn, pn1
        integer :: j, i
        real(wp) :: eps_

        eps_ = epsm
        if (present(eps)) eps_ = eps

        fn = nn
        csx = 0
        csa = 0
        cc = gamma(alf + 1.0_wp)
        tsx = fn*(fn + alf)
        tsa = cc
        do j = 2, nn
          cc = cc*c(j)
        end do
        do i = 1, nn
            select case(i)
            case(1) ! smallest zero
                xt = (1.0_wp + alf)*(3.0_wp + 0.92_wp*alf)/(1.0_wp + 2.4_wp*fn + 1.8_wp*alf)
            case(2) ! second zero
                xt = xt + (15.0_wp + 6.25_wp*alf)/(1.0_wp + 0.9_wp*alf + 2.5_wp*fn)
            case default ! all other zeros
                fi = i - 2
                r1 = (1.0_wp + 2.55_wp*fi)/(1.9_wp*fi)
                r2 = 1.26_wp*fi*alf/(1.0_wp + 3.5_wp*fi)
                ratio = (r1+r2)/(1.0_wp + 0.3_wp*alf)
                xt = xt + ratio*(xt - x(i-2))
            end select

            call lgroot(xt,nn,alf,dpn,pn1,b,c,eps_)
            x(i) = xt
            a(i) = cc/dpn/pn1
            write(*,'(f6.2,2i3,2(1x,e14.8),2x,2(1x,e14.8))') alf,nn,i,xt,a(i)
            csx = csx + xt
            csa = csa + a(i)
        end do
        write(*,'(f6.2,2i3,2(1x,e14.8),2x,2(1x,e14.8))') alf,nn,i,csx,csa,tsx,tsa

    end subroutine

!>
!   Improves the approximate root `x` of the Laguerre polynomial
!
!   In addition we also obtain the derivative of $L(n)$ at $x$
!   and the value of $L(n-1)$ at x
!
    pure subroutine lgroot(x,nn,alf,dpn,pn1,b,c,eps)
        real(wp), intent(inout) :: x        !! Approximate root value
        integer, intent(in) :: nn           !! Degree of polynomial
        real(wp), intent(in) :: alf         !! Exponent $\alfa$ of the Laguerre polynomial
        real(wp), intent(out) :: dpn        !! Derivative of L(n) at x
        real(wp), intent(out) :: pn1        !! Value of L(n-1) at x
        real(wp), intent(in) :: b(*), c(*)  !! Coefficients of the recurrence relation
        real(wp), intent(in) :: eps         !! Desired accuracy

        integer :: iter
        real(wp) :: d, p, dp

        do iter = 1, max_iter
            call lgrecur(p,dp,pn1,x,nn,alf,b,c)
            d = p/dp
            x = x - d
            if (abs(d/x) < eps) exit
        end do
        dpn = dp

    end subroutine lgroot

!>
!   Recurrence relation of the Laguerre polynomial
!
    pure subroutine lgrecur(pn,dpn,pn1,x,nn,alf,b,c)
        real(wp), intent(out) :: pn
        real(wp), intent(out) :: dpn
        real(wp), intent(out) :: pn1
        real(wp), intent(in) :: x
        integer, intent(in) :: nn
        real(wp), intent(in) :: alf
        real(wp), intent(in) :: b(*)
        real(wp), intent(in) :: c(*)

        integer :: j
        real(wp) :: p1, p, dp1, dp, q, dq

        p1 = 1.0_wp
        p = x - alf - 1.0_wp
        dp1 = 0.0_wp
        dp = 1.0_wp
        do j = 2, nn
            q = (x - b(j))*p - c(j)*p1
            dq = (x - b(j))*dp + p - c(j)*dp1
            p1 = p
            p = q
            dp1 = dp
            dp = dq
        end do
        pn = p
        dpn = dp
        pn1 = p1
    end subroutine


!>
!   Calculates the zeros x(i) of the nn-th order Hermite polynomial
!
!   The largest zero will be stored in `x(1)`. Also calculates the
!   corresponding coefficients `a(i)` of the nn-th order Gauss-Hermite
!   quadrature formula of degree `2*nn-1`.
!
!   The initial estimate for the largest zero is one given by
!   Szegö (1959, p. 131):
!
!   \[x_n \simeq x_n^* = (2n+1)^{1/2} - 1.85575(2n+1)^{-1/6}\]
!
!### References
! * Szegö, G., *Orthogonal Polynomials*, Am. Math. Soc. Colloquim Publ., **23**
!   (1959)
!
    subroutine hermit(nn,x,a,eps)
        integer, intent(in) :: nn           !! Number of points
        real(wp), intent(out) :: x(nn)      !! Zeros of the Hermite polynomial
        real(wp), intent(out) :: a(nn)      !! Weights of the quadrature formula
        real(wp), intent(in), optional :: eps         !! Desired accuracy

        real(wp) :: fn, cc, s, xt, dpn, pn1
        integer :: n1, n2, i, ni
        real(wp) :: eps_

        eps_ = epsm
        if (present(eps)) eps_ = eps

        fn = real(nn,wp)
        n1 = nn - 1
        n2 = (nn + 1)/2
        cc = sqrt(pi)*gamma(fn)/(2**n1)
        s = (2.0_wp*fn + 1.0_wp)**(1.0_wp/6.0_wp)
        do i = 1, n2
            select case(i)
            case(1)
                xt = s**3 - 1.85575_wp/(5.0_wp)!**(1._wp/6._wp)       ! largest zero
            case(2)
                xt = xt - 1.14_wp*fn**0.426_wp/xt  ! second zero
            case(3)
                xt = 1.86_wp*xt - 0.86_wp*x(1)     ! third zero
            case(4)
                xt = 1.91_wp*xt - 0.91_wp*x(2)     ! fourth zero
            case default
                xt = 2.0_wp*xt - x(i-2)         ! all other zeros
            end select

            call hroot(xt,nn,dpn,pn1,eps_)
            x(i) = xt
            a(i) = cc/dpn/pn1
            write(*,'(2i4,2(2x,e14.8))') nn,i,xt,a(i)
            ni = nn-i+1
            x(ni) = -xt
            a(ni) = a(i)
        end do
    end subroutine


!>
!   Improves the approximate root `x` of the Hermite polynomial
!
!   In addition we also obtain the derivative of $H(n)$ at $x$
!   and the value of $H(n-1)$ at x
!
    subroutine hroot(x,nn,dpn,pn1,eps)
        real(wp), intent(inout) :: x        !! Approximate root value
        integer, intent(in) :: nn           !! Degree of polynomial
        real(wp), intent(out) :: dpn        !! Derivative of H(n) at x
        real(wp), intent(out) :: pn1        !! Value of H(n-1) at x
        real(wp), intent(in) :: eps         !! Desired accuracy

        integer :: iter
        real(wp) :: p, dp, d

        do iter = 1, max_iter
            call hrecur(p,dp,pn1,x,nn)
            d = p/dp
            x = x - d
            if (abs(d) < eps) exit
        end do
        dpn = dp
    end subroutine

!>
!   Recurrence relation of the Hermite polynomial
!
    pure subroutine hrecur(pn,dpn,pn1,x,nn)
        real(wp), intent(out) :: pn, dpn, pn1
        real(wp), intent(in) :: x
        integer, intent(in) :: nn

        integer :: j
        real(wp) :: fj2
        real(wp) :: p1, p, dp1, dp, q, dq

        p1 = 1.0_wp
        p = x
        dp1 = 0.0_wp
        dp = 1.0_wp
        do j = 2, nn
            fj2 = (real(j,wp) - 1.0_wp)/2.0_wp
            q = x*p - fj2*p1
            dq = x*dp + p - fj2*dp1
            p1 = p
            p = q
            dp1 = dp
            dp = dq
        end do
        pn = p
        dpn = dp
        pn1 = p1
    end subroutine

end module