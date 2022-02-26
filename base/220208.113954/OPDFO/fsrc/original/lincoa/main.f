C     Calculate the tetrahedron of least volume that encloses the points
C       (XP(J),YP(J),ZP(J)), J=1,2,...,NP. Our method requires the origin
C       to be strictly inside the convex hull of these points. There are
C       twelve variables that define the four faces of each tetrahedron
C       that is considered. Each face has the form ALPHA*X+BETA*Y+GAMMA*Z=1,
C       the variables X(3K-2), X(3K-1) and X(3K) being the values of ALPHA,
C       BETA and GAMMA for the K-th face, K=1,2,3,4. Let the set T contain
C       all points in three dimensions that can be reached from the origin
C       without crossing a face. Because the volume of T may be infinite,
C       the objective function is the smaller of FMAX and the volume of T,
C       where FMAX is set to an upper bound on the final volume initially.
C       There are 4*NP linear constraints on the variables, namely that each
C       of the given points (XP(J),YP(J),ZP(J)) shall be in T. Let XS = min
C       XP(J), YS = min YP(J), ZS = min ZP(J) and SS = max XP(J)+YP(J)+ZP(J),
C       where J runs from 1 to NP. The initial values of the variables are
C       X(1)=1/XS, X(5)=1/YS, X(9)=1/ZS, X(2)=X(3)=X(4)=X(6)=X(7) =X(8)=0
C       and X(10)=X(11)=X(12)=1/SS, which satisfy the linear constraints,
C       and which provide the bound FMAX=(SS-XS-YS-ZS)**3/6. Other details
C       of the test calculation are given below, including the choice of
C       the data points (XP(J),YP(J),ZP(J)), J=1,2,...,NP. The smaller final
C       value of the objective function in the case NPT=35 shows that the
C       problem has local minima.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON FMAX
      DIMENSION XP(50),YP(50),ZP(50),A(12,200),B(200),X(12),W(500000)
C
C     Set some constants.
C
      ONE=1.0D0
      TWO=2.0D0
      ZERO=0.0D0
      PI=4.0D0*DATAN(ONE)
      IA=12
      N=12
C
C     Set the data points.
C
      NP=50
      SUMX=ZERO
      SUMY=ZERO
      SUMZ=ZERO
      DO 10 J=1,NP
      THETA=DFLOAT(J-1)*PI/DFLOAT(NP-1)
      XP(J)=DCOS(THETA)*DCOS(TWO*THETA)
      SUMX=SUMX+XP(J)
      YP(J)=DSIN(THETA)*DCOS(TWO*THETA)
      SUMY=SUMY+YP(J)
      ZP(J)=DSIN(TWO*THETA)
   10 SUMZ=SUMZ+ZP(J)
      SUMX=SUMX/DFLOAT(NP)
      SUMY=SUMY/DFLOAT(NP)
      SUMZ=SUMZ/DFLOAT(NP)
      DO 20 J=1,NP
      XP(J)=XP(J)-SUMX
      YP(J)=YP(J)-SUMY
   20 ZP(J)=ZP(J)-SUMZ
C
C     Set the linear constraints.
C
      M=4*NP
      DO 30 K=1,M
      B(K)=ONE
      DO 30 I=1,N
   30 A(I,K)=ZERO
      DO 40 J=1,NP
      DO 40 I=1,4
      K=4*J+I-4
      IW=3*I
      A(IW-2,K)=XP(J)
      A(IW-1,K)=YP(J)
   40 A(IW,K)=ZP(J)
C
C     Set the initial vector of variables. The JCASE=1,6 loop gives six
C       different choices of NPT when LINCOA is called.
C
      XS=ZERO
      YS=ZERO
      ZS=ZERO
      SS=ZERO
      DO 50 J=1,NP
      XS=DMIN1(XS,XP(J))
      YS=DMIN1(YS,YP(J))
      ZS=DMIN1(ZS,ZP(J))
   50 SS=DMAX1(SS,XP(J)+YP(J)+ZP(J))
      FMAX=(SS-XS-YS-ZS)**3/6.0D0
      DO 80 JCASE=1,6
      DO 60 I=2,8
   60 X(I)=ZERO
      X(1)=ONE/XS
      X(5)=ONE/YS
      X(9)=ONE/ZS
      X(10)=ONE/SS
      X(11)=ONE/SS
      X(12)=ONE/SS
C
C     Call of LINCOA, which provides the printing given at the end of this
C       note.
C
      NPT=5*JCASE+10
      RHOBEG=1.0D0
      RHOEND=1.0D-6
      IPRINT=1
      MAXFUN=10000
      PRINT 70, NPT,RHOEND
   70 FORMAT (//4X,'Output from LINCOA with  NPT =',I4,
     1  '  and  RHOEND =',1PD12.4)
      CALL LINCOA (N,NPT,M,A,IA,B,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
   80 CONTINUE
      STOP
      END
