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