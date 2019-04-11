        SUBROUTINE FACT1(A1,B1,C1,F1,U1,p2,p3,i1,i2)
C*-------------------------------------------
C*-------------------------------------------
      REAL *8 ZNAM
      REAL *8 X1(1000),Y1(1000),YS1(1000)
      REAL *8 A1(1000),B1(1000),C1(1000),F1(1000),U1(1000)
C*-------------------------------------------
      X1(I1-1)=0.
      Y1(I1-1)=P2
      DO1 IPK=I1,I2
      ZNAM=1./(B1(IPK)+A1(IPK)*X1(IPK-1))
      X1(IPK)=-C1(IPK)*ZNAM
    1 Y1(IPK)=(F1(IPK)-A1(IPK)*Y1(IPK-1))*ZNAM
      U1(I2+1)=P3
      DO2 IPK=I1,I2+1
      KREAL=I2-IPK+I1
    2 U1(KREAL)=U1(KREAL+1)*X1(KREAL)+Y1(KREAL)
      RETURN
      END

        SUBROUTINE TrFACT1(A1,B1,C1,F1,U1,p2,p3,i1,i2)
C*-------------------------------------------
C*-------------------------------------------
      include 'model.par'
      REAL *8 ZNAM
      REAL *8 X1(1000),Y1(1000,ntrac)
      REAL *8 A1(1000),B1(1000),C1(1000),F1(1000,ntrac),U1(1000,ntrac)
C*-------------------------------------------
      X1(I1-1)=0.
	Do nt=1,ntrac
         Y1(I1-1,nt)=P2
      Enddo
      DO 1 IPK=I1,I2
         ZNAM=1./(B1(IPK)+A1(IPK)*X1(IPK-1))
         X1(IPK)=-C1(IPK)*ZNAM
	   Do 1 nt=1,ntrac
    1 Y1(IPK,nt)=(F1(IPK,nt)-A1(IPK)*Y1(IPK-1,nt))*ZNAM
      Do nt=1,ntrac
         U1(I2+1,nt)=P3
      Enddo
      DO 2 IPK=I1,I2+1
      KREAL=I2-IPK+I1
	Do 2 nt=1,ntrac
    2 U1(KREAL,nt)=U1(KREAL+1,nt)*X1(KREAL)+Y1(KREAL,nt)
      RETURN
      END
