      SUBROUTINE HERMIT(NN,X,A,EPS)
C                  CALCULATES THE ZEROS  X(I)  OF THE NN-TH ORDER
C                HERMITE POLYNOMIAL.  THE LARGEST ZERO WILL BE
C                STORED IN X(1).  ALSO CALCULATES THE CORRESPONDING
C                COEFFICIENTS  A(I)  OF THE NN-TH ORDER GAUSS-HERMITE
C                QUADRATURE FORMULA OF DEGREE 2*NN-1.
C
      DIMENSION X(50),A(50)
      FN = NN
      N1 = NN - 1
      N2 = (NN+1)/2
      CC = 1.7724538509*GAMMA(FN)/(2.**N1)
      S = (2.*FN+1.)**.1667
      DO 10 I=1,N2
      IF ( I-1 ) 10,1,2
C                  LARGEST ZERO
    1 XT = S**3 - 1.85575/5
      GO TO 9
    2 IF ( I-2 ) 10,3,4
C                  SECOND ZERO
    3 XT = XT - 1.14*FN**.426/XT
      GO TO 9
    4 IF ( I-3 ) 10,5,6
C                  THIRD ZERO
    5 XT = 1.86*XT - .86*X(1)
      GO TO 9
    6 IF ( I-4 ) 10,7,8
C                  FOURTH ZERO
    7 XT = 1.91*XT - .91*X(2)
      GO TO 9
C                  ALL OTHER ZEROS
    8 XT = 2.*XT - X(I-2)
C
    9 CALL HROOT(XT,NN,DPN,PN1,EPS)
      X(I) = XT
      A(I) = CC/DPN/PN1
      WRITE(*,20) NN,I,XT,A(I)
      NI = NN-I+1
      X(NI)=-XT
   10 A(NI) = A(I)
   20 FORMAT(2I4,2(2X,E14.8))
      RETURN
      END   