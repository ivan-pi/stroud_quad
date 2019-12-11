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

      SUBROUTINE ROOT(X,NN,ALF,BTA,DPN,PNl,B,C,EPS)
C                  IMPROVES THE APPROXIMATE ROOT X
C                IN ADDITION WE ALSO OBTAIN
C                    DPN = DERIVATIVE OF P(N) AT X
C                    PN1 = VALUE OF P(N-l) AT X
C
      DIMENSION B(50),C(50)
      ITER = 0
    1 ITER = ITER + 1
      CALL RECUR(P,DP,PNl,X,NN,ALF,BTA,B,C)
      D = P/DP
      X = X - D
      IF( ABS(D)-EPS )  3,3,2
    2 IF( ITER - 10 )   1,3,3
    3 DPN = DP
      RETURN
      END

      SUBROUTINE RECUR(PN,DPN,PN1,X,NN,ALF,BTA,B,C)
      DIMENSION B(50),C(50)
      P1  = 1.
      P   = X + (ALF-BTA)/(ALF+BTA+2.)
      DP1 = 0.
      DP  = 1.
      DO 1 J=2,NN
      Q   = (X - B(J))*P - C(J)*P1
      DQ  = (X - B(J))*DP + P - C(J)*DP1
      P1  = P
      P   = Q
      DP1 = DP
    1 DP = DQ
      PN = P
      DPN = DP
      PN1 = P1
      RETURN
      END

      FUNCTION FLGAMA(W)
C                  CALCULATES LOG(BASE E)GAMMA(W) FOR W REAL AND
C                GAMMA(W) POSITIVE* USES STIRLlNGtS APPROXIMATION*
C                ACCURATE TO ABOUT 12 SIGNIFICANT PLACES*
C
      PI=3.141592653589793
      X=W
      M=0
      FK=-1.
      IF(X-.5)1,2,2
C          W LESS EQ .5
    1 M=1
      XPI=X*PI
      X=1.-X
    2 FK=FK+1.
      IF(X+FK-6.)2,2,3
    3 Z=X+FK
      ZZ=Z*Z
C          LOG GAMMA(Z). Z GREATER 6.
      Y=(Z-.5)*LOG(Z)-Z+ .9189385332047+(((((-4146./ZZ+1820.)/ZZ
     1         -1287.)/ZZ+1716.)/ZZ-6006.)/ZZ+180180.)/Z/2162160.
      IF(FK)6,6,4
    4 IK=FK
      DO 5 I=1,IK
      FK=FK-1.
    5 Y=Y-LOG(X+FK)
    6 IF(M)7,11,7
    7 P=PI/SIN(XPI )
      IF(P)8,8,10
    8 WRITE(*,9) W
    9 FORMAT(2X,6HGAMMA(E11.4,13H) IS NEGATIVE)
      Y=0.
      GO TO 11
   10 Y=LOG(P)-Y
   11 FLGAMA=Y
      RETURN
      END 

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


      SUBROUTINE LGROOT(X,NN,ALF,DPN,PN1,B,C,EPS)
C                  IMPROVES THE APPROXIMATE ROOT  C
C                IN ADDITION WE ALSO OBTAIN
C                    DPN = DERIVATIVE OF P(N) AT X
C                    PN1 0 VALUE OF P(N-1) AT X
C
      DIMENSION B(50),C(50)
      ITER = 0
    1 ITER = ITER + 1
      CALL LGRECR(P,DP,PN1,X,NN,ALF,B,C)
      D = P/DP
      X = X - D
      IF ( ABS(D/X)-EPS ) 3,3,2
    2 IF ( ITER - 10 )  1,3,3
    3 DPN = DP
      RETURN
      END


      SUBROUTINE LGRECR(PN,DPN,PN1,X,NN,ALF,B,C)
      DIMENSION  B(50),C(50)
      P1  = 1.
      P   = X - ALF - 1.
      DP1 = 0.
      DP  = 1.
      DO 1 J=2,NN
      Q   = (X - B(J))*P - C(J)*P1
      DQ  = (X - B(J))*DP + P - C(J)*DP1
      P1 = P
      P = Q
      DP1 = DP
    1 DP = DQ
      PN = P
      DPN = DP
      PN1 = P1
      RETURN
      END  


      FUNCTION GAMMAM(X)
C                  COMPUTES THE GAMMA FUNCTION OF X BY HASTINGS
C                APPROXIMATION.  X  MUST BE IN RANGE 0 TO 70
C
      GAM(Y) = (((((((.035868343*Y - .193527818)*Y + .482199394)*Y -
     1  .756704078)*Y + .918206857)*Y - .897056937)*Y + .988205891)*Y
     2 - .577191652)*Y + 1.0
      Z = X
      IF( Z ) 1,1,4
    1 GAMMAM = 0.
      WRITE(*,2) z
    2 FORMAT(2X,19HARG ERROR FOR GAMMA , E16.8)
      GO TO 14
    4 IF( Z - 70. ) 6,1,1
    6 IF( Z - 1. ) 8,7,9
    7 GAMMAM = 1.
      GO TO 14
    8 GAMMAM = GAM(Z)/Z
      GO TO 14
    9 ZA = 1.
   10 Z = Z-1.
      IF( Z-1. ) 13,11,12
   11 GAMMAM = ZA
      GO TO 14
   12 ZA = ZA*Z
      GO TO 10
   13 GAMMAM = ZA*GAM(Z)
   14 RETURN
      END   


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


      SUBROUTINE HROOT(X,NN,DPN,PN1,EPS)
C                  IMPROVES THE APPROXIMATE ROOT  X
C                IN ADDITION WE ALSO OBTAIN
C                    DPN = DERIVATIVE OF H(N) AT X
C                    PN1 = VALUE OF H(N-1) AT X
C
      ITER = 0
    1 ITER = ITER + 1
      CALL HRECUR(P,DP,PN1,X,NN)
      D = P/DP
      X = X - D
      IF ( ABS(D)-EPS ) 3,3,2
    2 IF ( ITER - 10 ) 1,3,3
    3 DPN = DP
      RETURN
      END


      SUBROUTINE HRECUR(PN,DPN,PN1,X,NN)
      P1 = 1.
      P = X
      DP1 = 0.
      DP = 1.
      DO 1 J=2,NN
      FJ = J
      FJ2 = (FJ-1.)/2.
      Q = X*P - FJ2*P1
      DQ = X*DP + P - FJ2*DP1
      P1 = P
      P = Q
      DP1 = DP
    1 DP = DQ
      PN = P
      DPN = DP   
      PN1 = P1
      RETURN
      END