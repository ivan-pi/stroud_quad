      SUBROUTINE JACOBI(NN,X,A,ALF,BTA,B,C,EPS,CSX,CSA,TSX,TSA)
C                  CALCULATES THE ZEROS X(I) OF THE NN-TH ORDER
C                JACOBI POLYNOMIAL PN(ALF,BTA) FOR THE SEGMENT (-1,1).
C                THE LARGEST ZERO WILL BE STORED IN X(l). ALSO
C                CALCULATES THE CORRESPONDING COEFFICIENTS  A(I)
C                OF THE NN-TH ORDER GAUSS-JACOBI QUADRATURE FORMULA
C                OF DEGREE 2*NN-1.
C                  THIS SUBROUTINE MUST BE GIVEN THE COEFFICIENTS
C
C                                  (ALF+BTA)(BTA-ALF)
C                      B(N) =  --------------------------
C                              (ALF+BTA+2N)(ALF+BTA+2N-2)
C
C                            4(N-1MALF+N-1) ( BTA+N-1) (ALF+BTA+N-1)
C                C(N) =  ---------------------------------------------
C                        (ALF+BTA+2N-1)(ALF+BTA+2N-2)**2(ALF+BTA+2N-3)
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
      FN = NN
      CSX = 0.
      CSA = 0.
      BETA = EXP(FLGAMA(ALF+1.) + FLGAMA(BTA+1.) - FLGAMA(ALF+BTA+2.))
      CC = 2.**(ALF+BTA+1.)*BETA
      TSX = FN*(BTA-ALF)/(ALF+BTA+2.*FN)
      TSA = CC
      DO 1 J=2,NN
    1 CC = CC*C(J)
      DO 12 I=1,NN
      IF( I-1 ) 12,2,3
C                  LARGEST ZERO
    2 AN = ALF/FN
      BN = BTA/FN
      R1 = (1.+ALF)*(2.78/(4.+FN*FN) + .768*AN/FN)
      R2 = 1. + 1.48*AN + .96*BN + .452*AN*AN + .83*AN*BN
      XT = 1. - R1/R2
      GO TO 11
    3 IF( I-2 ) 12,4,5
C                  SECOND ZERO
    4 R1 = (4.1+ALF)/((1.+ALF)*(1.+.156*ALF))
      R2 = 1. + .06*(FN-8.)*(1.+.12*ALF)/FN
      R3 = 1. + .012*BTA*(1.+.25*ABS(ALF))/FN
      RATIO = R1*R2*R3
      XT = XT - RATIO*(1.-XT)
      GO TO 11
    5 IF( I-3 ) 12,6,7
C                  THIRD ZERO
    6 R1 = (1.67+.28*ALF)/(1.+.37*ALF)
      R2 = 1. + .22*(FN-8.)/FN
      R3 = 1. + 8.*BTA/((6.28+BTA)*FN*FN)
      RATIO = R1*R2*R3
      XT = XT - RATIO*(X(1) - XT)
      GO TO 11
    7 IF(NN-I-1) 10,9,8
C                  MIDDLE ZEROS
    8 XT = 3.*X(I-1) - 3.*X(I-2) + X(I-3)
      GO TO 11
C                  SECOND LAST ZERO
    9 R1 = (1. + .235*BTA)/(.766+.119*BTA)
      R2 = 1./( 1. + .639*(FN-4.)/(1.+.71*(FN-4.)) )
      R3 = 1./( 1. + 20.*ALF/((7.5+ALF)*FN*FN) )
      RATIO = R1*R2*R3
      XT = XT + RATIO*(XT-X(I-2))
      GO TO 11
C                  LAST ZERO
   10 R1 = (1.+.37*BTA)/(1.67+.28*BTA)
      R2 = 1./( 1. + .22*(FN-8.)/FN )
      R3 = 1./( 1. + 8.*ALF/((6.28+ALF)*FN*FN) )
      RATIO = R1*R2*R3
      XT = XT + RATIO*(XT-X(I-2))
C
   11 CALL ROOT(XT,NN,ALF,BTA,DPN,PN1,B,C,EPS)
      X(I) = XT
      A(I) = CC/(DPN*PN1)
      WRITE(*,20) ALF,BTA,NN,I,XT,A(I)
      CSX = CSX + XT
   12 CSA = CSA + A(I)
      WRITE(*,20) ALF,BTA,NN,I,CSX,CSA,TSX,TSA
   20 FORMAT(2F6.2,2I3,2(1X,E14.8),2X,2(1X,E14.8))
      RETURN
      END