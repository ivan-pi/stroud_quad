module quad_mod

    implicit none
    private

    public :: jacobi,laguer,hermit

    public :: legendre

    interface
        subroutine jacobi(nn,x,a,alf,bta,b,c,eps,csx,csa,tsx,tsa)
            integer, intent(in) :: nn
            real, intent(out) :: x(nn), a(nn)
            real, intent(in) :: alf,bta
            real, intent(in) :: b(nn),c(nn)
            real, intent(in) :: eps
            real, intent(out) :: csx,csa,tsx,tsa
        end subroutine
        subroutine laguer(nn,x,a,alf,b,c,eps,csx,csa,tsx,tsa)
            integer, intent(in) :: nn
            real, intent(out) :: x(nn), a(nn)
            real, intent(in) :: alf
            real, intent(in) :: b(nn),c(nn)
            real, intent(in) :: eps
            real, intent(out) :: csx,csa,tsx,tsa
        end subroutine
        subroutine hermit(nn,x,a,eps)
            integer, intent(in) :: nn
            real, intent(out) :: x(nn), a(nn)
            real, intent(in) :: eps
        end subroutine
    end interface

contains


    subroutine legendre(nn,alf,bta,b,c)
        integer, intent(in) :: nn
        real, intent(in) :: alf, bta
        real, intent(out) :: b(nn), c(nn)
        integer :: n
        real :: num, denom

        do n = 1, nn
            b(n) = (alf + bta)*(bta - alf)/(2*n + alf + bta)/(2*n + alf + bta - 2)

            if (n == 2) then
                c(n) = 4.*(alf + 1)*(bta + 1)/(alf + bta + 3)/(alf + bta + 2)**2
            else
                num = 4.*(n-1)*(n+alf-1)*(n+bta-1)*(n+alf+bta-1)
                denom = (2*n+alf+bta-1)*(2*n + alf +bta - 2)**2*(2*n + alf + bta - 3)
                c(n) = num/denom
            end if
        end do
    end subroutine

end module

program main

    use quad_mod
    implicit none

    integer, parameter :: n = 5
    real :: x(n),a(n),eps
    real :: b(n), c(n)
    real :: alf, bta,csx,csa,tsx,tsa
    integer :: i

    eps = epsilon(eps)

    print *, "Hermite-Gauss quadrature"
    call hermit(n,x,a,eps)

    print *, "Legendre quadrature"
    alf = 0.
    bta = 0.

    call legendre(n,alf,bta,b,c)
    call jacobi(n,x,a,alf,bta,b,c,eps,csx,csa,tsx,tsa)


    print *, "Chebyshev first kind quadrature"

    alf = -0.5
    bta = -0.5
    b = 0
    c(2) = 0.5
    c(3:) = 0.25

    call jacobi(n,x,a,alf,bta,b,c,eps,csx,csa,tsx,tsa)


    print *, "Chebyshev second kind quadrature"

    alf = -0.5
    bta = -0.5
    b = 0
    c(2:) = 0.25

    call jacobi(n,x,a,alf,bta,b,c,eps,csx,csa,tsx,tsa)


    print *, "Generalized Laguerre-Gauss Quadrature"
    alf = 0.5
    bta = 0.5
    alf = 0.0
    b = [(alf + 2*i - 1,i=1,n)]
    c = [((i-1)*(alf+i-1),i=1,n)]
    call laguer(n,x,a,alf,b,c,eps,csx,csa,tsx,tsa)

end program