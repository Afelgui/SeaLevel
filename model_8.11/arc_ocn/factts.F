      SUBROUTINE FACTTS(A1,B1,C1,F1,Fs1,U1,Us1,p2,ps2,p3,i1,i2)
C*-------------------------------------------
C*-------------------------------------------
      REAL *8 ZNAM
      REAL *8 X1(1000),Y1(1000),YS1(1000)
      REAL *8 A1(1000),B1(1000),C1(1000),F1(1000),U1(1000),
     *        FS1(1000),US1(1000)
C*-------------------------------------------
      X1(I1-1)=0.
      Y1(I1-1)=P2
      YS1(I1-1)=PS2
      DO 1 IPK=I1,I2
      ZNAM=1.D+0/(B1(IPK)+A1(IPK)*X1(IPK-1))
      X1(IPK)=-C1(IPK)*ZNAM
      YS1(IPK)=(FS1(IPK)-A1(IPK)*YS1(IPK-1))*ZNAM
      Y1(IPK)=(F1(IPK)-A1(IPK)*Y1(IPK-1))*ZNAM
c	if (i1.eq.2.and.i2.eq.6.and.ipk.gt.1.and.k.eq.2) then
c      print *, ipk, (34.-X1(IPK)*34.-Ys1(IPK))
c      end if
1     continue
      U1(I2+1)=P3
      US1(I2+1)=P3
      DO 2 IPK=I1,I2
      KREAL=I2-IPK+I1
      US1(KREAL)=US1(KREAL+1)*X1(KREAL)+YS1(KREAL)
      U1(KREAL)=U1(KREAL+1)*X1(KREAL)+Y1(KREAL)
c		if (i1.eq.2.and.i2.eq.6.and.kreal.gt.1.and.k.eq.2) then
c      print*,kreal,(us1(kreal)-X1(kreal)*us1(kreal+1)-Ys1(kreal))
c    *                               ,us1(kreal)
c      end if
    2  continue
      RETURN
      END
