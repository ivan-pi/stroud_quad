      SUBROUTINE LAGUER(NN,X,A,ALF,B,C,EPS,CSX,CSA,TSX,TSA)
C                  CALCULATES THE ZEROS  X(I)  OF THE NN-TH ORDER
C                LAGUERRE POLYNOMIAL LN(ALF) FOR THE SEGMENT (0,INF)
C                THE SMALLEST ZERO WILL BE STORED IN X(1). ALSO
C                CALCULATES THE CORRESPONDING COEFFICIENTS  A(I)
C                OF THE NN-TH ORDER LAGUERRE QUADRATURE FORMULA
C                OF DEGREE 2*NN-1.
C                  THIS SUBROUTINE MUST BE GIVEN THE COEFFICIENTS
C
C                           B(N) =  (ALF + 2N - 1)
C
C                        C(N) = (N-1)(ALF + N - 1)
C
C                IN THE RECURSION RELATION
C
C                      P(N) = (X - B(N))*P(N-1) - C(N)*P(N-2)
C
C                FOR ALL N LESS THAN OR EQUAL TO THE HIGHEST DEGREE NN.
C
C                    CSX = CALC SUM X(I)    TSX = TRUE SUM X(I)
C                    CSA = CALC SUM A(I)    TSA = TRUE SUM A(I)
C
      DIMENSION X(50),A(50), B(50),C(50)
      EXTERNAL GAMMA
      FN = NN
      CSX = 0.
      CSA = 0.
      CC = GAMMA(ALF+1.)
      TSX = FN*(FN+ALF)
      TSA = CC
      DO 1 J=2,NN
    1 CC = CC*C(J)
      DO 7 I=1,NN
      IF( I-1 )  6,2,3
C                  SMALLEST ZERO
    2 XT = (1.+ALF)*(3.+0.92*ALF)/(1.+2.4*FN+1.8*ALF)
      GO TO 6
    3 IF( I-2 )  6,4,5
C                  SECOND ZERO
    4 XT = XT + (15.+6.25*ALF)/(1.+9.*ALF+2.5*FN)
      GO TO 6
C                  ALL OTHER ZEROS
    5 FI = I - 2
      R1 = (1.+2.55*FI)/(1.9*FI)
      R2 = 1.26*FI*ALF/(1.+3.5*FI)
      RATIO = (R1+R2)/(1.+.3*ALF)
      XT = XT + RATIO*(XT-X(I-2))
C
    6 CALL LGROOT(XT,NN,ALF,DPN,PN1,B,C,EPS)
      X(I) = XT
      A(I) = CC/DPN/PN1
      WRITE(*,20) ALF,NN,I,XT,A(I)
      CSX = CSX + XT
    7 CSA = CSA + A(I)
      WRITE(*,20) ALF,NN,I,CSX,CSA,TSX,TSA
   20 FORMAT( F6.2,2I3,2(1X,E14.8),2X,2(1X,E14.8))
      RETURN
      END