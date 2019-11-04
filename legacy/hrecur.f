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