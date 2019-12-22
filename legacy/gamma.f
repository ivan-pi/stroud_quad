      FUNCTION GAMMA(X)
C                  COMPUTES THE GAMMA FUNCTION OF  X  BY HASTINGS
C                APPROXIMATION.  X  MUST BE IN RANGE  0 TO 70.
C
      GAM(Y) = (((((((.035868343*Y - .193527818)*Y + .482199394)*Y -
     1  .756704078)*Y + .918206857)*Y - .897056937)*Y + .988205891)*Y
     2 - .577191652)*Y + 1.0
      Z = X
      IF( Z ) 1,1,4
    1 GAMMA = 0.
      WRITE(*,2) Z
    2 FORMAT(2X,19HARG ERROR FOR GAMMA , E16.8)
      GO TO 14
    4 IF ( Z - 70. ) 6,1,1
    6 IF ( Z - 1. ) 8,7,9
    7 GAMMA = 1.
      GO TO 14
    8 GAMMA = GAM(Z)/Z
      GO TO 14
    9 ZA = 1.
   10 Z = Z-1.
      IF( Z-1. ) 13,11,12
   11 GAMMA = ZA
      GO TO 14
   12 ZA=ZA*Z
      GO TO 10
   13 GAMMA = ZA*GAM(Z)
   14 RETURN
      END