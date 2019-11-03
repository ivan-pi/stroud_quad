module quad_mod

    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    private

    public :: hermite, laguerre, wp, pi

    real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
    real(wp), parameter :: eps = epsilon(1.0_wp)
    integer, parameter :: max_iter = 10

contains

    !                  CALCULATES THE ZEROS X(I) OF THE NN-TH ORDER
    !                JACOBI POLYNOMIAL PN(ALF,BTA) FOR THE SEGMENT (-1,1).
    !                THE LARGEST ZERO WILL BE STORED IN X(l). ALSO
    !                CALCULATES THE CORRESPONDING COEFFICIENTS  A(I)
    !                OF THE NN-TH ORDER GAUSS-JACOBI QUADRATURE FORMULA
    !                OF DEGREE 2*NN-1.
    !                  THIS SUBROUTINE MUST BE GIVEN THE COEFFICIENTS

    !                                  (ALF+BTA)(BTA-ALF)
    !                      B(N) =  --------------------------
    !                              (ALF+BTA+2N)(ALF+BTA+2N-2)

    !                            4(N-1MALF+N-1) ( BTA+N-1) (ALF+BTA+N-1)
    !                C(N) =  ---------------------------------------------
    !                        (ALF+BTA+2N-1)(ALF+BTA+2N-2)**2(ALF+BTA+2N-3)

    !                IN THE RECURSION RELATION

    !                      P(N) = (X - B(N))*P(N-1) - C(N)*P(N-2)

    !                FOR ALL N LESS THAN OR EQUAL TO THE HIGHEST DEGREE NN.

    !                    CSX = CALC SUM X(I)    TSX = TRUE SUM X(I)
    !                    CSA = CALC SUM A(I)    TSA = TRUE SUM A(I)
    subroutine jacobi(nn,x,a,alf,bta,b,c,eps,csx,csa,tsx,tsa)
        integer, intent(in) :: nn
        real(wp), intent(out) :: x(50)
        real(wp), intent(out) :: a(50)
        real(wp), intent(in) :: alf
        real(wp), intent(in) :: bta
        real(wp), intent(in) :: b(50)
        real(wp), intent(in) :: c(50)
        real(wp), intent(in) :: eps
        real(wp), intent(out) :: csx
        real(wp), intent(out) :: csa
        real(wp), intent(out) :: tsx
        real(wp), intent(out) :: tsa

        integer :: i
        real(wp) :: fn, beta, cc
        real(wp) :: an, bn, r1, r2, r3, ratio, xt

        fn = nn
        csx = 0.0_wp
        csa = 0.0_wp
        beta = exp(log_gamma(alf+1.0_wp) + log_gamma(bta+1.0_wp) - log_gamma(alf+bta+2.0_wp))
        cc = 2.0_wp**(alf+bta+1.0_wp)*beta
        tsx = fn*(bta - alf)/(alf + bta + 2.0_wp*fn)
        tsa = cc
        do  j=2,nn
          cc = cc*c(j)
        end do
    do i = 1, nn

        select case(i)
        case(1)                  ! largest zero
            an = alf/fn
            bn = bta/fn
            r1 = (1. + alf)*(2.78/(4. + fn*fn) + .768*an/fn)
            r2 = 1. + 1.48*an + .96*bn + .452*an*an + .83*an*bn
            xt = 1. - r1/r2 
        case(2)                  ! second zero
            r1 = (4.1+alf)/((1.+alf)*(1.+.156*alf))
            r2 = 1. + .06*(fn-8.)*(1.+.12*alf)/fn
            r3 = 1. + .012*bta*(1.+.25*ABS(alf))/fn
            ratio = r1*r2*r3
            xt = xt - ratio*(1.-xt)
        case(3)
            r1 = (1.67+.28*alf)/(1.+.37*alf)
            r2 = 1. + .22*(fn-8.)/fn
            r3 = 1. + 8.*bta/((6.28+bta)*fn*fn)
            ratio = r1*r2*r3
            xt = xt - ratio*(x(1) - xt)
        case(nn-1) ! second last zero
            r1 = (1. + .235*bta)/(.766+.119*bta)
            r2 = 1./( 1. + .639*(fn-4.)/(1.+.71*(fn-4.)) )
            r3 = 1./( 1. + 20.*alf/((7.5+alf)*fn*fn) )
            ratio = r1*r2*r3
            xt = xt + ratio*(xt-x(i-2))
        case(nn) ! last zero
            r1 = (1.+.37*bta)/(1.67+.28*bta)
            r2 = 1./( 1. + .22*(fn-8.)/fn )
            r3 = 1./( 1. + 8.*alf/((6.28+alf)*fn*fn) )
            ratio = r1*r2*r3
            xt = xt + ratio*(xt-x(i-2))
        case default !middle zeros
            xt = 3.*x(i-1) - 3.*x(i-2) + x(i-3)
        end select
        call root(xt,nn,alf,bta,dpn,pn1,b,c,eps)
        x(i) = xt
        a(i) = cc/(dpn*pn1)
        WRITE(*,20) alf,bta,nn,i,xt,a(i)
        csx = csx + xt
        csa = csa + a(i)

      IF( i-1  < 0) THEN
        GO TO    12
      ELSE IF ( i-1  == 0) THEN
        GO TO     2
      ELSE
        GO TO     3
      END IF
    !                  LARGEST ZERO
      2 an = alf/fn
      bn = bta/fn
      r1 = (1.+alf)*(2.78/(4.+fn*fn) + .768*an/fn)
      r2 = 1. + 1.48*an + .96*bn + .452*an*an + .83*an*bn
      xt = 1. - r1/r2
      GO TO 11
      3 IF( i-2  < 0) THEN
        GO TO    12
      ELSE IF ( i-2  == 0) THEN
        GO TO     4
      ELSE
        GO TO     5
      END IF
    !                  SECOND ZERO
      4 r1 = (4.1+alf)/((1.+alf)*(1.+.156*alf))
      r2 = 1. + .06*(fn-8.)*(1.+.12*alf)/fn
      r3 = 1. + .012*bta*(1.+.25*ABS(alf))/fn
      ratio = r1*r2*r3
      xt = xt - ratio*(1.-xt)
      GO TO 11
      5 IF( i-3  < 0) THEN
        GO TO    12
      ELSE IF ( i-3  == 0) THEN
        GO TO     6
      ELSE
        GO TO     7
      END IF
    !                  THIRD ZERO
      6 r1 = (1.67+.28*alf)/(1.+.37*alf)
      r2 = 1. + .22*(fn-8.)/fn
      r3 = 1. + 8.*bta/((6.28+bta)*fn*fn)
      ratio = r1*r2*r3
      xt = xt - ratio*(x(1) - xt)
      GO TO 11
      7 IF(nn-i-1 < 0) THEN
        GO TO    10
      ELSE IF (nn-i-1 == 0) THEN
        GO TO     9
      END IF
    !                  MIDDLE ZEROS
      8 xt = 3.*x(i-1) - 3.*x(i-2) + x(i-3)
      GO TO 11
    !                  SECOND LAST ZERO
      9 r1 = (1. + .235*bta)/(.766+.119*bta)
      r2 = 1./( 1. + .639*(fn-4.)/(1.+.71*(fn-4.)) )
      r3 = 1./( 1. + 20.*alf/((7.5+alf)*fn*fn) )
      ratio = r1*r2*r3
      xt = xt + ratio*(xt-x(i-2))
      GO TO 11
    !                  LAST ZERO
      10 r1 = (1.+.37*bta)/(1.67+.28*bta)
      r2 = 1./( 1. + .22*(fn-8.)/fn )
      r3 = 1./( 1. + 8.*alf/((6.28+alf)*fn*fn) )
      ratio = r1*r2*r3
      xt = xt + ratio*(xt-x(i-2))
      
      11 CALL root(xt,nn,alf,bta,dpn,pn1,b,c,eps)
        x(i) = xt
        a(i) = cc/(dpn*pn1)
        WRITE(*,20) alf,bta,nn,i,xt,a(i)
        csx = csx + xt
        csa = csa + a(i)
    end do

    write(*,20) alf,bta,nn,i,csx,csa,tsx,tsa
20  format(2f6.2,2i3,2(1x,e14.8),2x,2(1x,e14.8))

    end subroutine

    !                  IMPROVES THE APPROXIMATE ROOT X
    !                IN ADDITION WE ALSO OBTAIN
    subroutine root(x,nn,alf,bta,dpn,pn1,b,c,eps)
        real(wp), intent(inout) :: x
        integer, intent(in) :: nn
        real(wp), intent(in) :: alf
        real(wp), intent(in) :: bta
        real(wp), intent(out) :: dpn !! DERIVATIVE OF P(N) AT X
        real(wp), intent(out) :: pn1 !! VALUE OF P(N-l) AT X
        real(wp), intent(in) :: b(*)
        real(wp), intent(in) :: c(*)
        real(wp), intent(in) :: eps

        integer :: iter
        real(wp) :: d, p, dp

        iter = 0
        do
            iter = iter + 1
            call recur(p,dp,pn1,x,nn,alf,bta,b,c)
            d = p/dp
            x = x - d
            if (abs(d) < eps) exit
            if (iter > max_iter) exit
        end do
        dpn = dp
    end subroutine

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


    !                  CALCULATES THE ZEROS  X(I)  OF THE NN-TH ORDER
    !                LAGUERRE POLYNOMIAL LN(ALF) FOR THE SEGMENT (0,INF)
    !                THE SMALLEST ZERO WILL BE STORED IN X(1). ALSO
    !                CALCULATES THE CORRESPONDING COEFFICIENTS  A(I)
    !                OF THE NN-TH ORDER LAGUERRE QUADRATURE FORMULA
    !                OF DEGREE 2*NN-1.
    !                  THIS SUBROUTINE MUST BE GIVEN THE COEFFICIENTS

    !                           B(N) =  (ALF + 2N - 1)

    !                        C(N) = (N-1)(ALF + N - 1)

    !                IN THE RECURSION RELATION

    !                      P(N) = (X - B(N))*P(N-1) - C(N)*P(N-2)

    !                FOR ALL N LESS THAN OR EQUAL TO THE HIGHEST DEGREE NN.

    !                    CSX = CALC SUM X(I)    TSX = TRUE SUM X(I)
    !                    CSA = CALC SUM A(I)    TSA = TRUE SUM A(I)
    subroutine laguerre(nn,x,a,alf,b,c,eps,csx,csa,tsx,tsa)
        integer, intent(in) :: nn
        real(wp), intent(out) :: x(nn)
        real(wp), intent(out) :: a(nn)
        real(wp), intent(in) :: alf
        real(wp), intent(in) :: b(nn)
        real(wp), intent(in) :: c(nn)
        real(wp), intent(in) :: eps
        real(wp), intent(out) :: csx
        real(wp), intent(out) :: csa
        real(wp), intent(out) :: tsx
        real(wp), intent(out) :: tsa

        real(wp) :: fn, cc, xt, fi, r1, r2, ratio, dpn, pn1
        integer :: j, i

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
            case(1)
                xt = (1.0_wp + alf)*(3.0_wp + 0.92_wp*alf)/(1.0_wp + 2.4_wp*fn + 1.8_wp*alf) ! smallest zero
            case(2)
                xt = xt + (15.0_wp + 6.25_wp*alf)/(1.0_wp + 9.0_wp*alf + 2.5_wp*fn)     ! second zero
            case default
                fi = i - 2                                      ! all other zeros                   
                r1 = (1.0_wp + 2.55_wp*fi)/(1.9_wp*fi)
                r2 = 1.26_wp*fi*alf/(1.0_wp + 3.5_wp*fi)
                ratio = (r1+r2)/(1.0_wp + 0.3_wp*alf)
                xt = xt + ratio*(xt - x(i-2))
            end select
          
            call lgroot(xt,nn,alf,dpn,pn1,b,c,eps)
            x(i) = xt
            a(i) = cc/dpn/pn1
            write(*,'(f6.2,2i3,2(1x,e14.8),2x,2(1x,e14.8))') alf,nn,i,xt,a(i)
            csx = csx + xt
            csa = csa + a(i)
        end do
        write(*,'(f6.2,2i3,2(1x,e14.8),2x,2(1x,e14.8))') alf,nn,i,csx,csa,tsx,tsa

    end subroutine


    !                  IMPROVES THE APPROXIMATE ROOT  C
    !                IN ADDITION WE ALSO OBTAIN
    pure subroutine lgroot(x,nn,alf,dpn,pn1,b,c,eps)
        real(wp), intent(inout) :: x
        integer, intent(in) :: nn
        real(wp), intent(in) :: alf
        real(wp), intent(out) :: dpn !! DERIVATIVE OF P(N) AT X
        real(wp), intent(out) :: pn1 !! VALUE OF P(N-1) AT X
        real(wp), intent(in) :: b(*)
        real(wp), intent(in) :: c(*)
        real(wp), intent(in) :: eps

        integer :: iter
        real(wp) :: d, p, dp

        iter = 0
        do
            iter = iter + 1
            call lgrecur(p,dp,pn1,x,nn,alf,b,c)
            d = p/dp
            x = x - d
            if (abs(d/x) < eps) exit
            if (iter > max_iter) exit
        end do
        dpn = dp

    end subroutine lgroot


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


    !                  CALCULATES THE ZEROS  X(I)  OF THE NN-TH ORDER
    !                HERMITE POLYNOMIAL.  THE LARGEST ZERO WILL BE
    !                STORED IN X(1).  ALSO CALCULATES THE CORRESPONDING
    !                COEFFICIENTS  A(I)  OF THE NN-TH ORDER GAUSS-HERMITE
    !                QUADRATURE FORMULA OF DEGREE 2*NN-1.
    subroutine hermite(nn,x,a,eps)
        integer, intent(in) :: nn
        real(wp), intent(out) :: x(nn)
        real(wp), intent(out) :: a(nn)
        real(wp), intent(in) :: eps

        real(wp) :: fn, cc, s, xt, dpn, pn1
        integer :: n1, n2, i, ni

        fn = real(nn,wp)
        n1 = nn - 1
        n2 = (nn + 1)/2
        cc = sqrt(pi)*gamma(fn)/(2**n1)
        s = (2.0_wp*fn + 1.0_wp)**(1.0_wp/6.0_wp)
        do i = 1, n2
            select case(i)
            case(1)
                xt = s**3 - 1.85575_wp/5.0_wp       ! largest zero
            case(2)
                xt = xt - 1.14_wp*fn**0.426_wp/xt  ! second zero
            case(3)
                xt = 1.86_wp*xt - 0.86_wp*x(1)     ! third zero
            case(4)
                xt = 1.91_wp*xt - 0.91_wp*x(2)     ! fourth zero
            case default
                xt = 2.0_wp*xt - x(i-2)         ! all other zeros
            end select
              
            call hroot(xt,nn,dpn,pn1,eps)
            x(i) = xt
            a(i) = cc/dpn/pn1
            write(*,'(2i4,2(2x,e14.8))') nn,i,xt,a(i)
            ni = nn-i+1
            x(ni) = -xt
            a(ni) = a(i)
        end do
    end subroutine hermite


    !>                  IMPROVES THE APPROXIMATE ROOT  X
    !                IN ADDITION WE ALSO OBTAIN
    !                    DPN = DERIVATIVE OF H(N) AT X
    !                    PN1 = VALUE OF H(N-1) AT X
    pure subroutine hroot(x,nn,dpn,pn1,eps)
        real(wp), intent(inout) :: x       !! the approximate root x
        integer, intent(in) :: nn
        real(wp), intent(out) :: dpn       !! derivative of H(N) at x
        real(wp), intent(out) :: pn1       !! value of H(n-1) at x
        real(wp), intent(in) :: eps

        integer :: iter
        real(wp) :: p, dp, d

        iter = 0
        do
            iter = iter + 1
            call hrecur(p,dp,pn1,x,nn)
            d = p/dp
            x = x - d
            if (abs(d) < eps) exit
            if (iter > max_iter) exit
        end do
        dpn = dp
    end subroutine


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

program test_stroud
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use quad_mod, only: hermite, laguerre, pi
    implicit none

    real(wp), parameter :: eps = epsilon(1.0_wp)
    integer, parameter :: n = 2
    real(wp) :: x(n), a(n), aa(n)
    real(wp) :: b(n), c(n)
    real(wp) :: alf, bta,csx,csa,tsx,tsa
    integer :: i


    call hermite(n,x,a,eps)

    call get_hermite(n,aa)
    print *, "Hermite quadrature"
    do i = 1, n
        print *, i, x(i), a(i), aa(i)
    end do

    print *, "Generalized Laguerre-Gauss Quadrature"

    alf = 0.0_wp
    b = [(alf + 2*i - 1,i=1,n)]
    c = [((i-1)*(alf+i-1),i=1,n)]
    call laguerre(n,x,a,alf,b,c,eps,csx,csa,tsx,tsa)

contains

    subroutine get_hermite(nn,a)
        integer, intent(in) :: nn
        real(wp), intent(out) :: a(nn)

        select case(nn)
        case(2)
            a(1) = 0.5_wp*sqrt(pi)
            a(2) = a(1)
        case(3)
            a(1) = sqrt(pi)/6.0_wp
            a(2) = 2.0_wp*sqrt(pi)/3.0_wp
            a(3) = a(1)
        case(4)
            a(1) = sqrt(pi)/(4.0_wp*(3.0_wp + sqrt(6.0_wp)))
            a(2) = sqrt(pi)/(4.0_wp*(3.0_wp - sqrt(6.0_wp)))
            a(3) = a(2)
            a(4) = a(1)
        case default
            return
        end select
    end subroutine
end program