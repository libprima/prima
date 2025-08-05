C------------------------------------------------------------------------------
C     Main program of test problems in Report DAMTP 1992/NA5.
C------------------------------------------------------------------------------
      COMMON NPROB
      DIMENSION X(10),XOPT(10),W(3000),IACT(51)
      DO 180 NPROB=1,10
      IF (NPROB .EQ. 1) THEN
C
C     Minimization of a simple quadratic function of two variables.
C
          PRINT 10
   10     FORMAT (/7X,'Output from test problem 1 (Simple quadratic)')
          N=2
          M=0
          XOPT(1)=-1.0
          XOPT(2)=0.0
      ELSE IF (NPROB .EQ. 2) THEN
C
C     Easy two dimensional minimization in unit circle.
C
          PRINT 20
   20     FORMAT (/7X,'Output from test problem 2 (2D unit circle ',
     1      'calculation)')
          N=2
          M=1
          XOPT(1)=SQRT(0.5)
          XOPT(2)=-XOPT(1)
      ELSE IF (NPROB .EQ. 3) THEN
C
C     Easy three dimensional minimization in ellipsoid.
C
          PRINT 30
   30     FORMAT (/7X,'Output from test problem 3 (3D ellipsoid ',
     1      'calculation)')
          N=3
          M=1
          XOPT(1)=1.0/SQRT(3.0)
          XOPT(2)=1.0/SQRT(6.0)
          XOPT(3)=-1.0/3.0
      ELSE IF (NPROB .EQ. 4) THEN
C
C     Weak version of Rosenbrock's problem.
C
          PRINT 40
   40     FORMAT (/7X,'Output from test problem 4 (Weak Rosenbrock)')
          N=2
          M=0
          XOPT(1)=-1.0
          XOPT(2)=1.0
      ELSE IF (NPROB .EQ. 5) THEN
C
C     Intermediate version of Rosenbrock's problem.
C
          PRINT 50
   50     FORMAT (/7X,'Output from test problem 5 (Intermediate ',
     1      'Rosenbrock)')
          N=2
          M=0
          XOPT(1)=-1.0
          XOPT(2)=1.0
      ELSE IF (NPROB .EQ. 6) THEN
C
C     This problem is taken from Fletcher's book Practical Methods of
C     Optimization and has the equation number (9.1.15).
C
          PRINT 60
   60     FORMAT (/7X,'Output from test problem 6 (Equation ',
     1      '(9.1.15) in Fletcher)')
          N=2
          M=2
          XOPT(1)=SQRT(0.5)
          XOPT(2)=XOPT(1)
      ELSE IF (NPROB .EQ. 7) THEN
C
C     This problem is taken from Fletcher's book Practical Methods of
C     Optimization and has the equation number (14.4.2).
C
          PRINT 70
   70     FORMAT (/7X,'Output from test problem 7 (Equation ',
     1      '(14.4.2) in Fletcher)')
          N=3
          M=3
          XOPT(1)=0.0
          XOPT(2)=-3.0
          XOPT(3)=-3.0
      ELSE IF (NPROB .EQ. 8) THEN
C
C     This problem is taken from page 66 of Hock and Schittkowski's book Test
C     Examples for Nonlinear Programming Codes. It is their test problem Number
C     43, and has the name Rosen-Suzuki.
C
          PRINT 80
   80     FORMAT (/7X,'Output from test problem 8 (Rosen-Suzuki)')
          N=4
          M=3
          XOPT(1)=0.0
          XOPT(2)=1.0
          XOPT(3)=2.0
          XOPT(4)=-1.0
      ELSE IF (NPROB .EQ. 9) THEN
C
C     This problem is taken from page 111 of Hock and Schittkowski's
C     book Test Examples for Nonlinear Programming Codes. It is their
C     test problem Number 100.
C
          PRINT 90
   90     FORMAT (/7X,'Output from test problem 9 (Hock and ',
     1      'Schittkowski 100)')
          N=7
          M=4
          XOPT(1)=2.330499
          XOPT(2)=1.951372
          XOPT(3)=-0.4775414
          XOPT(4)=4.365726
          XOPT(5)=-0.624487
          XOPT(6)=1.038131
          XOPT(7)=1.594227
      ELSE IF (NPROB .EQ. 10) THEN
C
C     This problem is taken from page 415 of Luenberger's book Applied
C     Nonlinear Programming. It is to maximize the area of a hexagon of
C     unit diameter.
C
          PRINT 100
  100     FORMAT (/7X,'Output from test problem 10 (Hexagon area)')
          N=9
          M=14
      END IF
      DO 160 ICASE=1,2
      DO 120 I=1,N
  120 X(I)=1.0
      RHOBEG=0.5
      RHOEND=0.001
      IF (ICASE .EQ. 2) RHOEND=0.0001
      IPRINT=1
      MAXFUN=2000
      CALL COBYLA (N,M,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,IACT)
      IF (NPROB .EQ. 10) THEN
          TEMPA=X(1)+X(3)+X(5)+X(7)
          TEMPB=X(2)+X(4)+X(6)+X(8)
          TEMPC=0.5/SQRT(TEMPA*TEMPA+TEMPB*TEMPB)
          TEMPD=TEMPC*SQRT(3.0)
          XOPT(1)=TEMPD*TEMPA+TEMPC*TEMPB
          XOPT(2)=TEMPD*TEMPB-TEMPC*TEMPA
          XOPT(3)=TEMPD*TEMPA-TEMPC*TEMPB
          XOPT(4)=TEMPD*TEMPB+TEMPC*TEMPA
          DO 130 I=1,4
  130     XOPT(I+4)=XOPT(I)
      END IF
      TEMP=0.0
      DO 140 I=1,N
  140 TEMP=TEMP+(X(I)-XOPT(I))**2
      PRINT 150, SQRT(TEMP)
  150 FORMAT (/5X,'Least squares error in variables =',1PE16.6)
  160 CONTINUE
      PRINT 170
  170 FORMAT (2X,'----------------------------------------------',
     1  '--------------------')
  180 CONTINUE
      STOP
      END
