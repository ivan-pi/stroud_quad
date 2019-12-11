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