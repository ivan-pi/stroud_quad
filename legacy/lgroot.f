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