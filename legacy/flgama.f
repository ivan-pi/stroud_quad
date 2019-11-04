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