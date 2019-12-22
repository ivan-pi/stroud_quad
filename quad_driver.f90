program test_quad_mod

    use quad_mod, only: jacobi, hermit, laguer, pi, wp
    implicit none

    ! Number of quadrature points
    integer, parameter :: n = 5

    real(wp), parameter :: eps = epsilon(1.0_wp)
    real(wp) :: x(n), a(n), b(n), c(n)
    real(wp) :: alf,bta,csx,csa,tsx,tsa

    write(*,*) "Gauss-Legendre Quadrature"
    alf = 0.0_wp
    bta = 0.0_wp
    call r_legendre(n,alf,bta,b,c)
    call jacobi(n,x,a,alf,bta,b,c,eps,csx,csa,tsx,tsa)

    write(*,*)
    write(*,*) "Generalized Gauss-Laguerre Quadrature"
    alf = 0.0_wp
    call r_laguerre(n,alf,b,c)
    call laguer(n,x,a,alf,b,c,eps,csx,csa,tsx,tsa)

    write(*,*)
    write(*,*) "Gauss-Hermite Quadrature"
    call hermit(n,x,a,eps)

    write(*,*)
    write(*,*) "Gauss-Chebyshev Quadrature (first kind)"
    alf = -0.5_wp
    bta = -0.5_wp
    b = 0.0_wp
    c(2) = 0.5_wp
    c(3:) = 0.25_wp
    call jacobi(n,x,a,alf,bta,b,c,eps,csx,csa,tsx,tsa)

    write(*,*)
    write(*,*) "Gauss-Chebyshev Quadrature (second kind)"
    alf = 0.5_wp
    bta = 0.5_wp
    b = 0.0_wp
    c(2:) = 0.25_wp
    call jacobi(n,x,a,alf,bta,b,c,eps,csx,csa,tsx,tsa)

contains

    ! Laguerre polynomial recurrence relation
    subroutine r_laguerre(nn,alf,b,c)
        integer, intent(in) :: nn
        real(wp), intent(in) :: alf
        real(wp), intent(out) :: b(nn), c(nn)
        integer :: i
        do i = 1, nn
            b(i) = alf + 2*i - 1
            c(i) = (i-1)*(alf+i-1)
        end do
    end subroutine

    ! Legendre polynomial recurrence relation
    subroutine r_legendre(nn,alf,bta,b,c)
        integer, intent(in) :: nn
        real(wp), intent(in) :: alf, bta
        real(wp), intent(out) :: b(nn), c(nn)
        integer :: n
        real(wp) :: num, denom

        do n = 1, nn
            b(n) = (alf + bta)*(bta - alf)/(2*n + alf + bta)/(2*n + alf + bta - 2)

            if (n == 2) then
                num = 4*(alf + 1)*(bta + 1)
                denom = (alf + bta + 3)*(alf + bta + 2)**2
                c(n) = num/denom
            else
                num = 4*(n-1)*(n+alf-1)*(n+bta-1)*(n+alf+bta-1)
                denom = (2*n+alf+bta-1)*(2*n + alf +bta - 2)**2*(2*n + alf + bta - 3)
                c(n) = num/denom
            end if
        end do
    end subroutine

    ! Second, third, and fourth Hermite quadrature rules
    ! subroutine get_hermite(nn,a)
    !     integer, intent(in) :: nn
    !     real(wp), intent(out) :: a(nn)

    !     select case(nn)
    !     case(2)
    !         a(1) = 0.5_wp*sqrt(pi)
    !         a(2) = a(1)
    !     case(3)
    !         a(1) = sqrt(pi)/6.0_wp
    !         a(2) = 2.0_wp*sqrt(pi)/3.0_wp
    !         a(3) = a(1)
    !     case(4)
    !         a(1) = sqrt(pi)/(4.0_wp*(3.0_wp + sqrt(6.0_wp)))
    !         a(2) = sqrt(pi)/(4.0_wp*(3.0_wp - sqrt(6.0_wp)))
    !         a(3) = a(2)
    !         a(4) = a(1)
    !     case default
    !         return
    !     end select
    ! end subroutine

end program