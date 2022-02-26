      SUBROUTINE CALCFC (N,M,X,F,CON)
      COMMON NPROB
      DIMENSION X(*),CON(*)
      IF (NPROB .EQ. 1) THEN
C
C     Test problem 1 (Simple quadratic)
C     
          F=10.0*(X(1)+1.0)**2+X(2)**2
      ELSE IF (NPROB .EQ. 2) THEN
C
C    Test problem 2 (2D unit circle calculation)
C
          F=X(1)*X(2)
          CON(1)=1.0-X(1)**2-X(2)**2
      ELSE IF (NPROB .EQ. 3) THEN
C
C     Test problem 3 (3D ellipsoid calculation)
C
          F=X(1)*X(2)*X(3)
          CON(1)=1.0-X(1)**2-2.0*X(2)**2-3.0*X(3)**2
      ELSE IF (NPROB .EQ. 4) THEN
C
C     Test problem 4 (Weak Rosenbrock)
C
          F=(X(1)**2-X(2))**2+(1.0+X(1))**2
      ELSE IF (NPROB .EQ. 5) THEN
C
C     Test problem 5 (Intermediate Rosenbrock)
C
          F=10.0*(X(1)**2-X(2))**2+(1.0+X(1))**2
      ELSE IF (NPROB .EQ. 6) THEN
C
C     Test problem 6 (Equation (9.1.15) in Fletcher's book)
C
          F=-X(1)-X(2)
          CON(1)=X(2)-X(1)**2
          CON(2)=1.0-X(1)**2-X(2)**2
      ELSE IF (NPROB .EQ. 7) THEN
C
C     Test problem 7 (Equation (14.4.2) in Fletcher's book)
C
          F=X(3)
          CON(1)=5.0*X(1)-X(2)+X(3)
          CON(2)=X(3)-X(1)**2-X(2)**2-4.0*X(2)
          CON(3)=X(3)-5.0*X(1)-X(2)
      ELSE IF (NPROB .EQ. 8) THEN
C
C     Test problem 8 (Rosen-Suzuki)
C
          F=X(1)**2+X(2)**2+2.0*X(3)**2+X(4)**2-5.0*X(1)-5.0*X(2)
     1      -21.0*X(3)+7.0*X(4)
          CON(1)=8.0-X(1)**2-X(2)**2-X(3)**2-X(4)**2-X(1)+X(2)
     1      -X(3)+X(4)
          CON(2)=10.0-X(1)**2-2.0*X(2)**2-X(3)**2-2.0*X(4)**2+X(1)+X(4)
          CON(3)=5.0-2.0*X(1)**2-X(2)**2-X(3)**2-2.0*X(1)+X(2)+X(4)
      ELSE IF (NPROB .EQ. 9) THEN
C
C     Test problem 9 (Hock and Schittkowski 100)
C
          F=(X(1)-10.0)**2+5.0*(X(2)-12.0)**2+X(3)**4+3.0*(X(4)-11.0)**2
     1      +10.0*X(5)**6+7.0*X(6)**2+X(7)**4-4.0*X(6)*X(7)-10.0*X(6)
     2      -8.0*X(7)
          CON(1)=127.0-2.0*X(1)**2-3.0*X(2)**4-X(3)-4.0*X(4)**2-5.0*X(5)
          CON(2)=282.0-7.0*X(1)-3.0*X(2)-10.0*X(3)**2-X(4)+X(5)
          CON(3)=196.0-23.0*X(1)-X(2)**2-6.0*X(6)**2+8.0*X(7)
          CON(4)=-4.0*X(1)**2-X(2)**2+3.0*X(1)*X(2)-2.0*X(3)**2-5.0*X(6)
     1      +11.0*X(7)
      ELSE IF (NPROB .EQ. 10) THEN
C
C     Test problem 10 (Hexagon area)
C
          F=-0.5*(X(1)*X(4)-X(2)*X(3)+X(3)*X(9)-X(5)*X(9)+X(5)*X(8)
     1      -X(6)*X(7))
          CON(1)=1.0-X(3)**2-X(4)**2
          CON(2)=1.0-X(9)**2
          CON(3)=1.0-X(5)**2-X(6)**2
          CON(4)=1.0-X(1)**2-(X(2)-X(9))**2
          CON(5)=1.0-(X(1)-X(5))**2-(X(2)-X(6))**2
          CON(6)=1.0-(X(1)-X(7))**2-(X(2)-X(8))**2
          CON(7)=1.0-(X(3)-X(5))**2-(X(4)-X(6))**2
          CON(8)=1.0-(X(3)-X(7))**2-(X(4)-X(8))**2
          CON(9)=1.0-X(7)**2-(X(8)-X(9))**2
          CON(10)=X(1)*X(4)-X(2)*X(3)
          CON(11)=X(3)*X(9)
          CON(12)=-X(5)*X(9)
          CON(13)=X(5)*X(8)-X(6)*X(7)
          CON(14)=X(9)
      END IF
      RETURN
      END
