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