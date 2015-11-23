*********************************************************************
      SUBROUTINE EIGRS(A, N, N1, NE, NV, EPS, W, LW, E, V, IER)
************************************************************************
*  EIGENVALUES AND EIGENVECTORS FOR A SYMMETRIC MATRIX                 *
*     BY HOUSEHOLDER-BISECTION-INVERSE ITERATION METHOD.               *
*  PARAMETERS                                                          *  
*    (1) A: 2-DIM. ARRAY CONTAINING THE SYMMETRIC MATRIX               * 
*    (2) N: ORDER OF THE MATRIX (A)                                    *
*    (3) N1: ROW SIZE OF THE 2-DIM. ARRAY (A)                          *
*    (4) NE: NUMBER OF NEEDED EIGENVALUES                              *
*    (5) NV: NUMBER OF NEEDED EIGENVECTORS                             *
*    (6) EPS: TOLERANCE FOR CONVERGENCE                                *
*    (7) W: 2-DIM. WORKING ARRAY                                       *
*    (8) LW: 1-DIM. WORKING ARRAY                                      *
*    (9) E: 1-DIM. ARRAY CONTAINING COMPUTED EIGENVALUES               *
*   (10) V: 2-DIM. ARRAY CONTAINING COMPUTED EIGENVECTORS              *
*   (11) IER: ERROR CODE                                               *
*  COPYRIGHT   T. OGUNI   JUNE 30 1989    VERSION 1.0                  *
************************************************************************  
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 SW
      DIMENSION A(N1,N), W(N1,6), LW(N), E(1), V(N1,1)
C
      NEA = IABS(NE)
      NVA = IABS(NV)
      IF (NEA .EQ. 0 .OR. N1 .LT. N .OR. N .LT. NEA .OR. NEA .LT. 
     * NVA .OR. N .LT. 2) THEN
       WRITE(*,*) '(SUBR. EIGRS) INVALID ARGUMENT.',NV,NE,N,N1
       IER = 2
       RETURN
      ENDIF
      IF (EPS .LT. 0.0) EPS = 1.0D-6
      IF (N .EQ. 2) THEN
       W(1,1) = A(2,1)
       T = 0.5D0 * (A(1,1) + A(2,2))
       R = A(1,1) * A(2,2) - A(2,1) * A(2,1)
       D = T * T - R
       Q = DABS(T) + DSQRT(D)
       IF (T .LT. 0.0) Q = - Q
       T = T * DFLOAT(NE)
       IF (T .GE.0.0) THEN
        E(1) = Q
        IF (NEA .EQ. 2) E(2) = R / Q
       ELSE
        E(1) = R / Q
        IF (NEA .EQ. 2) E(2) = Q
       ENDIF
C  REDUCE TO TRIDIAGONAL FORM BY HOUSEHOLDER'S METHOD.
      ELSE
       DO 130 K=1,N-2
        S = 0.0D0
        DO 60 I=K+1,N
   60    S = S + A(I,K) * A(I,K)
        W(K,1) = 0.0D0
       IF (S .NE. 0) THEN
        SR = DSQRT(S)
        A1 = A(K+1,K)
        IF (A1 .LT. 0.0) SR = - SR
        W(K,1) = - SR
        R = 1.0D0 / (S + A1 * SR)
        A(K+1,K) = A1 + SR
        DO 90 I=K+1,N
         S = 0.0D0
         DO 70 J=K+1,I
   70     S = S + A(I,J) * A(J,K)
         IF (I .NE. N) THEN
          DO 80 J=I+1,N
   80      S = S + A(J,I) * A(J,K)
         ENDIF
   90    W(I,1) = S * R
         S = 0.0D0
         DO 100 I=K+1,N
  100     S = S + A(I,K) * W(I,1)
         T = 0.5D0 * S * R
         DO 110 I=K+1,N
  110     W(I,1) = W(I,1) - T * A(I,K)
         DO 120 J=K+1,N
          WJ1 = W(J,1)
          AJK = A(J,K)
          DO 120 I=J,N
  120      A(I,J) = A(I,J) - A(I,K)*WJ1 - W(I,1)*AJK
        ENDIF
  130 CONTINUE
        W(N-1,1) = A(N,N-1)
C  COMPUTE EIGENVALUES BY BISECTION METHOD.
         DO 135 I=1,N
  135     W(I,6) = A(I,I)
      R=DMAX1((DABS(W(1,6))+DABS(W(1,1))),(DABS(W(N-1,1))
     *   +DABS(W(N,6))))
      DO 140 I=2,N-1
       T=DABS(W(I-1,1))+DABS(W(I,6))+DABS(W(I,1))
       IF (T .GT. R) R = T
  140 CONTINUE
      EPS1 = R * 1.0D-6
      EPS2 = R * EPS
      DO 150 I=1,N-1
  150  W(I,2) = W(I,1) * W(I,1)
      IF (NE .LT. 0) R = - R
      F = R
      DO 160 I=1,NEA
  160  E(I) = - R
      DO 240 K=1,NEA
       D = E(K)
  170  T = 0.5D0 * (D + F)
      IF (DABS(D-F).GT.EPS2 .AND. T.NE.D .AND. T.NE.F) THEN
       J = 0
       I = 1
  180 Q = W(I,6) - T
  190 IF (Q .GE. 0.0) J = J + 1
      IF (Q .EQ. 0.0) GO TO 200
      I = I + 1
      IF (I .GT. N) GO TO 210
      Q = W(I,6) - T - W(I-1,2) / Q
      IF (L .NE. 1) GO TO 190
      J = J + 1
      I = I - 1
  200 I = I + 2
      IF (I .LE. N) GO TO 180
  210 IF (NE .LT. 0) J = N - J
      IF (J .LT. K) THEN
      F = T
      GO TO 170
      ENDIF
      D = T
      M = MIN0(J, NEA)
      DO 230 I=K,M
  230  E(I) = T
      GO TO 170
      ENDIF
  240 E(K) = T
      ENDIF
C  EIGENVECTORS BY INVERSE ITERATION.
        IF (NV .EQ. 0) RETURN
         IF (N .EQ. 2) THEN
          W(1,6) = A(1,1)
          W(2,6) = A(2,2)
         ENDIF
C
        W(N,1) = 0.0D0
        MM = 584287
        DO 410 I=1,NVA
         DO 260 J=1,N
          W(J,2) = W(J,6) - E(I)
          W(J,3) = W(J,1)
  260     V(J,I) = 1.0D0
         SW = 0
C  REDUCE TO TRIANGULAR FORM
         DO 280 J=1,N-1
          IF (DABS(W(J,2)) .GE. DABS(W(J,1))) THEN
           IF (W(J,2) .EQ. 0.0) W(J,2) = 1.0D-30
           W(J,5) = W(J,1) / W(J,2)
           LW(J) = 0
           W(J+1,2) = W(J+1,2) - W(J,5) * W(J,3)
           W(J,4) = 0.0D0
          ELSE
           W(J,5) = W(J,2) / W(J,1)
           LW(J) = 1
           W(J,2) = W(J,1)
           T = W(J,3)
           W(J,3) = W(J+1,2)
           W(J,4) = W(J+1,3)
           W(J+1,2) = T - W(J,5) * W(J,3)
           W(J+1,3) = - W(J,5) * W(J,4)
          ENDIF
  280 CONTINUE
          IF (W(N,2) .EQ. 0.0) W(N,2) = 1.0D-30
C  BEGIN BACK SUBSTITUTION
         IF (I .NE. 1) THEN
          IF (DABS(E(I) - E(I-1)) .LT. EPS1) THEN
C  GENERATE RANDOM NUMBERS
           DO 290 J=1,N
            MM = MM * 48828125
  290       V(J,I) = DFLOAT(MM) * 0.4656613D-9
          ENDIF
        ENDIF
  300 CONTINUE
         T = V(N,I)
         R = V(N-1,I)
  310    V(N,I) = T / W(N,2)
         V(N-1,I) = (R - W(N-1,3) * V(N,I)) / W(N-1,2)
         IF (L .EQ. 1) THEN
          DO 320 J=1,N-2
  320      V(J,I) = V(J,I) * 1.0D-5
          T = T * 1.0D-5
          R = R * 1.0D-5
          GO TO 310
         ENDIF 
         IF (N .NE. 2) THEN
         K = N - 2
  340    T = V(K,I)
  350    V(K,I) = (T - W(K,3)*V(K+1,I) - W(K,4)*V(K+2,I)) / W(K,2)
        IF (L .EQ. 1) THEN
         DO 360 J=1,N
  360     V(J,I) = V(J,I) * 1.0D-5
         T = T * 1.0D-5
         GO TO 350
        ENDIF
        K = K - 1
        IF (K .GT. 0) GO TO 340
       ENDIF
        IF (SW .NE. 0) THEN
         SW = 1
         DO 400 J=1,N-1
          IF (LW(J) .NE. 0) THEN
           V(J+1,I) = V(J+1,I) - W(J,5) * V(J,I)
          ELSE
           T = V(J,I)
           V(J,I) = V(J+1,I)
           V(J+1,I) = T - W(J,5) * V(J+1,I)
          ENDIF
  400    CONTINUE
        GO TO 300
       ENDIF
  410 CONTINUE
C  BEGIN BACK TRANSFORMATION
      IF (N .NE. 2) THEN
       DO 415 I=1,N-2
  415   W(I,1) = - W(I,1) * A(I+1,I)
       DO 460 I=1,NVA
        K = N - 2
  420   R = W(K,1)
        IF (R .NE. 0.0) THEN
         R = 1.0D0 / R
         S = 0.0D0
         DO 430 J=K+1,N
  430     S = S + A(J,K) * V(J,I)
         R = R * S
         DO 440 J=K+1,N
  440     V(J,I) = V(J,I) - R * A(J,K)
        ENDIF
        K = K - 1
        IF (K .GE. 1) GO TO 420
  460 CONTINUE
C  NORMALIZE EIGENVECTORS.  MAX ELEMENT IS 1.
      ENDIF
      DO 490 I=1,NVA
       T = DABS(V(1,I))
       K = 1
       DO 480 J=2,N
        R = DABS(V(J,I))
        IF (T .LT. R) THEN
         T = R
         K = J
        ENDIF
  480  CONTINUE
       T = 1.0D0 / V(K,I)
       DO 490 J=1,N
  490   V(J,I) = V(J,I) * T
       IF (NV .LT. 0) RETURN
C  ORTHONORMALIZE AS NORM IS 1.
      DO 550 I=1,NVA
       IF (I .NE. 1) THEN
        IF (DABS(E(I) - E(I-1)) .LT. EPS1) THEN
         DO 510 J=M,I-1
          S = 0.0D0
          DO 500 K=1,N
  500      S = S + V(K,J) * V(K,I)
          DO 510 K=1,N
  510      V(K,I) = V(K,I) - S * V(K,J)
        ELSE
         M = I
        ENDIF
       ENDIF
       S = 0.0D0
       DO 540 J=1,N
  540   S = S + V(J,I) * V(J,I)
       T = 0.0D0
       IF (S .NE. 0.0) T = DSQRT(1.0D0 / S)
       DO 550 J=1,N
  550   V(J,I) = V(J,I) * T
       RETURN
       END


C*********************************************************************
C*********************************************************************

      subroutine ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
c
      integer i,j,n,nm,ierr,matz
      double precision ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fm1(2,n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex hermitian matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex hermitian matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1, fv2, and  fm1  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
c
         do 30 j = 1, n
            zr(j,i) = 0.0d0
   30    continue
c
         zr(i,i) = 1.0d0
   40 continue
c
      call  tql2(nm,n,w,fv1,zr,ierr)
      if (ierr .ne. 0) go to 50
      call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
   50 return
      end

C*********************************************************************
C*********************************************************************
      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
      double precision f,g,h,fi,gi,hh,si,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix
c     to a real symmetric tridiagonal matrix using
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex hermitian input matrix.
c          only the lower triangle of the matrix need be supplied.
c
c     on output
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction in their full lower
c          triangles.  their strict upper triangles and the
c          diagonal of ar are unaltered.
c
c        d contains the diagonal elements of the the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        tau contains further information about the transformations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
c
      do 100 i = 1, n
  100 d(i) = ar(i,i)
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
c
         if (scale .ne. 0.0d0) go to 140
         tau(1,l) = 1.0d0
         tau(2,l) = 0.0d0
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
c
  140    do 150 k = 1, l
            ar(i,k) = ar(i,k) / scale
            ai(i,k) = ai(i,k) / scale
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
  150    continue
c
         e2(i) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = pythag(ar(i,l),ai(i,l))
c     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0d0) go to 160
         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         ar(i,l) = g * ar(i,l)
         ai(i,l) = g * ai(i,l)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(i,l) = g
  170    f = 0.0d0
c
         do 240 j = 1, l
            g = 0.0d0
            gi = 0.0d0
c     .......... form element of a*u ..........
            do 180 k = 1, j
               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
  180       continue
c
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
  240    continue
c
         hh = f / (h + h)
c     .......... form reduced a ..........
         do 260 j = 1, l
            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
c
            do 260 k = 1, j
               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k)
     x                           + fi * tau(2,k) + gi * ai(i,k)
               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k)
     x                           - fi * e(k) - gi * ar(i,k)
  260    continue
c
  270    do 280 k = 1, l
            ar(i,k) = scale * ar(i,k)
            ai(i,k) = scale * ai(i,k)
  280    continue
c
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * dsqrt(h)
  300 continue
c
      return
      end

C*********************************************************************
C*********************************************************************

C*********************************************************************
C*********************************************************************
C**** for old version, "send otqlrat from eispack"
C** From dana!moler Tue, 1 Sep 87 10:15:40 PDT
C** New TQLRAT
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      DOUBLE PRECISION D(N),E2(N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
C
C     This subroutine is a translation of the Algol procedure tqlrat,
C     Algorithm 464, Comm. ACM 16, 689(1973) by Reinsch.
C
C     This subroutine finds the eigenvalues of a symmetric
C     tridiagonal matrix by the rational QL method.
C
C     On input
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E2 contains the squares of the subdiagonal elements of the
C          input matrix in its last N-1 positions.  E2(1) is arbitrary.
C
C      On output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct and
C          ordered for indices 1,2,...IERR-1, but may not be
C          the smallest eigenvalues.
C
C        E2 has been destroyed.
C
C        IERR is set to
C          zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG for  DSQRT(A*A + B*B) .
C
C     Questions and comments should be directed to Burton S. Garbow,
C     Mathematics and Computer Science Div, Argonne National Laboratory
C
C     This version dated August 1987.
C     Modified by C. Moler to fix underflow/overflow difficulties,
C     especially on the VAX and other machines where epslon(1.0d0)**2
C     nearly underflows.  See the loop involving statement 102 and
C     the two statements just before statement 200.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0D0
      T = 0.0D0
      E2(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
         H = DABS(D(L)) + DSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
         if (c .ne. 0.0d0) go to 105
C        Spliting tolerance underflowed.  Look for larger value.
         do 102 i = l, n
            h = dabs(d(i)) + dsqrt(e2(i))
            if (h .gt. t) t = h
  102    continue
         b = epslon(t)
         c = b * b
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = PYTHAG(P,1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
C           Avoid division by zero on next pass
            if (g .eq. 0.0d0) g = epslon(d(i))
            h = g * (p / r)
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END

C*********************************************************************
C*********************************************************************

C*********************************************************************
C*********************************************************************
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end

C*********************************************************************
C*********************************************************************

C*********************************************************************
C*********************************************************************
      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
c
      integer i,j,k,l,m,n,nm
      double precision ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      double precision h,s,si
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  htridi.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction by  htridi  in their
c          full lower triangles except for the diagonal of ar.
c
c        tau contains further information about the transformations.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the hermitian
c                tridiagonal matrix. ..........
      do 50 k = 1, n
c
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
c
      if (n .eq. 1) go to 200
c     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0d0) go to 140
c
         do 130 j = 1, m
            s = 0.0d0
            si = 0.0d0
c
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
  110       continue
c     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
c
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end

C*********************************************************************
C*********************************************************************

C*********************************************************************
C*********************************************************************

      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
C*********************************************************************
C*********************************************************************

C*********************************************************************
C*********************************************************************
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end

C*********************************************************************
C*********************************************************************


