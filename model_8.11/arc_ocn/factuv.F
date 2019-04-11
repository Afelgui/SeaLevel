      subroutine factuv(a11,a12,c11,c12,b11,b22,b12,fu1,fv1,
     *	u1,v1,i1,m,n1,n2)
      Include 'model.par'
          common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      REAL *8 A11(kh,nh),A12(kh),C11(kh,nh),C12(kh)
     *       ,B11(kh,nh),B12(kh,nh),FU1(kh,nh),FV1(kh,nh)
     *       ,u1(kh,nh),v1(kh,nh),B22(kh,nh)
     *       ,X11(kh,nh),X12(kh,nh),X21(kh,nh),X22(kh,nh)
     *       ,Y1(kh),Y2(kh)
     *       ,W11(kh,nh),W12(kh,nh),W21(kh,nh),W22(kh,nh)
      REAL *8 det
      DO 39 N=N1,N2
      X12(I1-1,N)=0.
      X11(I1-1,N)=0.
      X21(I1-1,N)=0.
      X22(I1-1,N)=0.
       i2=kbtv(m,n)
      DO 5 K=I1,i2
      W11(K,N)=(B11(K,n)-A12(K)*X12(K-1,N)+A11(K,N)*X22(K-1,N))
      W22(K,N)=(B22(K,n)+A11(K,N)*X11(K-1,N)+A12(K)*X21(K-1,N))
      W12(K,N)=-(B12(K,n)+A11(K,N)*X12(K-1,N)+A12(K)*X22(K-1,N))
      W21(K,N)=-(-B12(K,n)-A12(K)*X11(K-1,N)+A11(K,N)*X21(K-1,N))
      ddet=W11(K,N)*W22(K,N)-W21(K,N)*W12(K,N)
      if(ddet.eq.0)print *, 'ddet=0 in adapt',k,n
      if(ddet.eq.0)read *
      DET=1./(W11(K,N)*W22(K,N)-W21(K,N)*W12(K,N))
      W11(K,N)=W11(K,N)*DET
      W22(K,N)=W22(K,N)*DET
      W12(K,N)=W12(K,N)*DET
      W21(K,N)=W21(K,N)*DET
      X11(K,N)=-(W11(K,N)*C11(K,n)-W12(K,N)*C12(K))
      X21(K,N)=-(W21(K,N)*C11(K,n)-W22(K,N)*C12(K))
      X12(K,N)=-(W11(K,N)*C12(K)+W12(K,N)*C11(K,n))
      X22(K,N)=-(W21(K,N)*C12(K)+W22(K,N)*C11(K,n))
    5 CONTINUE
   39 CONTINUE
      DO 40 N=N1,N2
      Y1(I1-1)=0.
      Y2(I1-1)=0.
      i2=kbtv(m,n)
      DO 45 K=i1,i2
      F1=FU1(K,N)-A11(K,N)*Y1(K-1)-A12(K)*Y2(K-1)
      F2=FV1(K,N)+A12(K)*Y1(K-1)-A11(K,N)*Y2(K-1)

		if (n.eq.64) then
	continue
	end if


      Y1(K)=W11(K,N)*F1+W12(K,N)*F2
      Y2(K)=W21(K,N)*F1+W22(K,N)*F2
   45 CONTINUE
      U1(I2+1,N)=0.
      V1(I2+1,N)=0.
      DO 6 K=I1,I2
      I=I2+I1-K
      U1(I,N)=U1(I+1,N)*X11(I,N)+V1(I+1,N)*X12(I,N)+Y1(I)
    6 V1(I,N)=U1(I+1,N)*X21(I,N)+V1(I+1,N)*X22(I,N)+Y2(I)
   40 CONTINUE
      RETURN
      END
