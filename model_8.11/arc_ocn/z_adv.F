	


      SUBROUTINE z_adv1(UU,VV,WW,WWvel,ftok,lvl)
      Include 'model.par'
      COMMON/DAY/Nyear,month,nday
      COMMON/NH2/ NH2
      real *8 z,dz,dzk,dx,dy,h,h1,tax,tay,sl,slb
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      common/Wind/Tax(mh,nh),Tay(mh,nh)
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      real *8 rho(nodes,kgr)
     *       ,UU(nodes,kgr),VV(nodes,kgr),WW(nodes,kgr)
     *       ,WWvel(nodes,kgr)
     *       ,Ftok(mh,nh),uxy(mh,nh),vxy(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt1
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      COMMON/V1/ROZ(mh,nh)
      common/v5/grad_x(mh,nh),grad_y(mh,nh)
      common/v6/Fui(mh,nh),Fvi(mh,nh),Fui1(mh,nh),Fvi1(mh,nh)	
	common/brt/UBrT(nodes), VBrT(nodes)
      Dimension IR1(kh,nh),IR2(kh,nh)
      REAL *8 A11(kh,nh),A12(kh),C11(kh,nh),C12(kh)
     *       ,B11(kh,nh),B12(kh,nh),FU1(kh,nh),FV1(kh,nh)
     *       ,u1(kh,nh),v1(kh,nh),w(kh,nh),B22(kh,nh)
     *       ,F1(kh),F2(kh),FU(kh),FV(kh)
     *       ,V(kh,nh),U(kh,nh)
     *       ,UN(kh)
     *       ,UL(kh,nh),VL(kh,nh),FSL(kh,nh)
     *       ,CZ(kh,nh)
     *       ,CVZ(kh,nh),HGM(kh,nh),VI(nh,kh)
     *       ,FSU(48),FSV(48),Fui,Fvi,fui1,fvi1
     *       ,grad_x,grad_y
     *       ,u_int(mh,nh),v_int(mh,nh)
     *       ,ubrt,vbrt
     *        ,lvl(mh,nh),lvl0(mh,nh),lvl_uv
     *          ,u_sum(mh,nh) 
      REAL *8 FUH,FVH,FU0,FV0,QC,T1,T2,T3,T4,DET,roz
      REAL *8 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,hx,hy
      REAL *8 DX0,DX1,DY0,DY1,DZ0,DZ1,ZZ1,ZZ2,XMULT,CW
      real *8 rav(4,kh),den(2,nh,kh),dlam1,dlam2,dtet1,dtet2
     *       ,pc(kh),dcent(kh),df(4,kh),dlam(kh),dtet(kh)
     *       ,plam(kh),ptet(kh),ta,sa,h0,hh1,hh2,vmod
      integer *2 kc(4,kh)
      dimension h4(4)
      integer nu_diag
      common/iodiag/nu_diag
C*-------------------------------------------
C
c	goto 618
       
	 ll=0
                         N11=N1-1
                          N21=N2+1
                           KZ1=2
                           KZ2=KZ0-2
                           KKZ1=KKZ-1
                           KKZ2=KKZ-2
C***********************************************
C***********************************************
c        Write(nu_diag,*)'begin adapt'

      DO 10 K=1,kh

   10 UN(K)=0.1
      UN(1)=1. 

                  do n=1,nh
                 do m=1,mh
                roz(m,n)=0.
               end do
              end do
      ll=0

      DO 200 m=1,mh-1

  550 FORMAT(1X,'M=',I3)
         DO 11 K=KZ1,KKZ
         DO 11 N=N0,NN
            IZ1=IPH(M,N)
            IZ2=IPV(M,N)
            IZ1=IZ1*IZ2
            IR1(K,N)=0
            IR2(K,N)=0
            SQ1=(KBT(M,N)/K)*(KBT(M+1,N)/K)*(KBT(M+1,N+1)/K)*IPH(M,N)
            SQ2=(KBT(M,N)/K)*(KBT(M,N+1)/K)*(KBT(M+1,N+1)/K)*IPV(M,N)
            IF(SQ1.GT.0)IR1(K,N)=1
            IF(SQ2.GT.0)IR2(K,N)=1
            IF(IZ1.EQ.0)GO TO 89
            IR1(K,N)=IR1(K,N)*IR2(K,N)
            IR2(K,N)=IR1(K,N)
   89    CONTINUE
   11    CONTINUE



                 IF(M.EQ.1)GOTO 55
c        Write(nu_diag,*)'27'
      DO 27 N=N0,NN
            K2=KP(M,N)
             k3=kp(m,n+1)
              k1=kp(m,n-1)
               k5=kp(m+1,n)
                k4=kp(m-1,n)
      DO 27 K=1,kgr
                 u(k+1,n)=uu(k2,k)
                  v(k+1,n)=vv(k2,k)
                   u1(k+1,n)=0.
                    v1(k+1,n)=0.
	k2=kp(m,n)
	k3=kp(m,n+1)
	k5=kp(m+1,n)
	k7=kp(m+1,n+1)
   27 w(K+1,N)=WWvel(k2,k)

  333  format(10(1x,10f6.2,/))
  101 LL=LL+1
c       print *,iydu(1,ll),ll
      IF(ll.le.kyLu.and.IYDu(1,LL).EQ.M)GOTO 100
      LL=LL-1
      GOTO 50
  100 N1=IYDu(2,LL)
      N2=IYDu(3,LL)
      N11=N1-1
      N21=N2+1
      DY1=DY(1,1)


c        Write(nu_diag,*)'42'
      DO 42 N=N1,N2
      DT=dt1
            KZV=KBTV(M,N)
            IF(N.EQ.N11) KZV=KBTV(M,N+1)
c           svm=SMB(N)
            H0=H1(M,N)
            hm=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.
      dy1=dly(m,n)
      dy0=dly(m+1,n)
            DO 43 K=KZ1,KZV
               ZZ1=1.
               DZ0=DZ(K-1)
               DZ1=DZ(K)
			 DZ2=DZK(K)

		     lvl_uv=(lvl(m,n)+lvl(m+1,n)+lvl(m,n+1)+lvl(m+1,n+1))/4.

			if (k.eq.kz1) dz1=dz(k)-lvl_uv/2.d+0			
	    	if (k.eq.kz1) dz2=DZK(K)-lvl_uv
              

c               IF(K.EQ.KZV) DZ2=0.5*DZ0
               IF(K.EQ.2)ZZ1=0.
               ZZ2=1
                if(k.eq.kzv) zz2=0
                if(hm.eq.0) then
                print *,m,n
                stop
                end if
               CZ(K,N)=-UN(K)/DZ1*hm*ZZ2
   30          HGM(K,N)=hm*dz2

   	if (m.eq.87.and.n.eq.55) then
	continue
	end if

c               CVZ(K,N)=w(k,n)*hm/2.
c               CZ(K,N)=CZ(K,N)-0.125*DABS(CVZ(K,N))
   43       CONTINUE
   42    CONTINUE
C***************************************
C***************************************
      FU1M=0.
      FV1M=0.
c        Write(nu_diag,*)'39'
      DO 39 N=N1,N2
                FU0=TAX(m,n)
                 FV0=TAY(m,n)
            hm=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.
      DX1=(dlx(m,n)+dlx(m,n+1))/2.
      DY1=(dly(m,n)+dly(m+1,n))/2.
                   K2=KP(M,N)
                     kzv=kbtv(m,n)
	hh1=h1(m,n)
         Z1=iPH(M,N)
          Z2=iPV(M,N)
           Z12=Z1*Z2
						 IF(Z12.EQ.0.)then
								hh1=0.1
						         do k=kz1,kkz
								  C11(K,n)=0
						          A11(K,N)=0
						            B11(K,n)=HGM(K,N)
								   B22(K,n)=HGM(K,N)

									 B12(K,n)=0
									  Fu1(K,n)=0
									 Fv1(K,n)=0
								end do
							  GOTO 39
							 end if
             H0=H1(M,N)
c                 SVM=SMB(N)
c                  SL0=-SLB(N)*H0/SVM
c               FUH=-0.5*R*(Ftok(M+1,N+1)-Ftok(M+1,N)+
c     *                     Ftok(m,n+1)-Ftok(m,n))/DY1/SVM/h0
c               FVH=0.5*R*(Ftok(m+1,n)-Ftok(m,n)+
c     *                    Ftok(m+1,n+1)-Ftok(m,n+1))/DX1*RAD/h0
c


   62              CW=DX(M,n)*DY(m,n)
    
c((((
      DO 24 K=KZ1,kzv
               DZ1=DZK(K)
c               IF(K.EQ.KZV) DZ1=0.5*DZ(K-1)
          C12(K)=0.
           C11(K,n)=CZ(K,N)
            A11(K,N)=CZ(K-1,N)
             A12(K)=0.
              B11(K,n)=HGM(K,N)/dt-c11(k,n)-a11(k,n)
              B22(K,n)=B11(k,n)
               B12(K,n)=0
                Fu1(K,n)=Hgm(k,n)/dt*u(k,n)
                 Fv1(K,n)=Hgm(k,n)/dt*v(k,n)

   24    CONTINUE
       vmod=sqrt(u(kzv,n)*u(kzv,n)+v(kzv,n)*v(kzv,n))
       fuh=-2.5e-3*vmod
	B11(kzv,n)=B11(kzv,n)
     *	-FUH*hm
	B22(kzv,n)=B22(kzv,n)
     *   	-FUH*hm
c                ))))
                  Fu1(KZ1,n)=Fu1(KZ1,n)+FU0*hm
                Fv1(KZ1,n)=Fv1(KZ1,n)+FV0*hm

	if (m.eq.87.and.n.eq.55) then
	continue
	end if


   39 continue
      i1=kz1
c      i2=kzv
c        Write(nu_diag,*)'factuv'
      call factuv(a11,a12,c11,c12,b11,b22,b12,fu1,fv1,
     *         	u1,v1,i1,m,n1,n2)
      do n=n1,n2
	hh1=h1(m,n)+0.000000001
                     kzv=kbtv(m,n)
      DX1=(dlx(m,n)+dlx(m,n+1))/2.
      DY1=(dly(m,n)+dly(m+1,n))/2.
c       vmod=sqrt(u(kzv,n)*u(kzv,n)+v(kzv,n)*v(kzv,n))
c       fuh=-2.5e-3*vmod*u(kzv,n)
c       fvh=-2.5e-3*vmod*v(kzv,n)

       do k=kz1,kzv
	fui(m,n)=fui(m,n)-
     *     	(a11(k,n)*u1(k-1,n)+c11(k,n)*u1(k+1,n)
     *        -(c11(k,n)+a11(k,n))*u1(k,n))/hh1/dx1/dy1
	fvi(m,n)=fvi(m,n)-
     *    	(a11(k,n)*v1(k-1,n)+c11(k,n)*v1(k+1,n)
     *        -(c11(k,n)+a11(k,n))*v1(k,n))/hh1/dx1/dy1
	end do
	end do
      go to 101
   50 continue
  660 FORMAT(48(1X,13F6.1,/))
      DO  K=1,kgr
      DO  N=N0,NN
      k2=kp(m,n) 
c      uInteg(k2,k)=uu(k2,k)
c      vInteg(k2,k)=vv(k2,k)
      UU(K2,K)=U1(K+1,N)
      VV(K2,K)=V1(K+1,N)
      end do
      end do
   55 CONTINUE
  200 CONTINUE
c        Write(nu_diag,*)'after 200'

       DO  Kl=1,kxLu

	if (kl.eq.148) then
	continue
	end if

         n=IxDu(1,Kl)
         I1=IxDu(2,Kl)
         I2=IxDu(3,Kl)
           DO m=I1,I2
       
       kzv=kbtv(m,n)-1
       ubot=uu(kp(m,n),kzv)
       vbot=vv(kp(m,n),kzv)



       vmod=sqrt(ubot*ubot+vbot*vbot)
        if(h1(m,n).gt.0.) then
       fui(m,n)=fui(m,n)-2.5e-3*vmod*ubot/h1(m,n)
       fvi(m,n)=fvi(m,n)-2.5e-3*vmod*vbot/h1(m,n)
       end if
           end do
       end do


618	 continue

	u_int=0
	v_int=0
	fui1=0
	fvi1=0
	u_sum=0.

	DO  Kl=1,kxLu
            n=IxDu(1,Kl)
            I1=IxDu(2,Kl)
            I2=IxDu(3,Kl)
          do m=I1,I2
             kzv=kbtv(m,n)
	     k2=kp(m,n)
		lvl_uv=(lvl(m,n)+lvl(m+1,n)+lvl(m,n+1)+lvl(m+1,n+1))/4.
	   u_int(m,n)= uu(kp(m,n),1)* (dzk(2)-lvl_uv)/(h1(m,n)-lvl_uv)
         v_int(m,n)= vv(kp(m,n),1)*(dzk(2)-lvl_uv)/(h1(m,n)-lvl_uv)

             do k=3,kzv
                dz1=dzk(k) 
         u_int(m,n)=u_int(m,n)+uu(kp(m,n),k-1)*dz1/(h1(m,n)-lvl_uv)
         v_int(m,n)=v_int(m,n)+vv(kp(m,n),k-1)*dz1/(h1(m,n)-lvl_uv)	  
           end do

			if (kbt(m,n).eq.kbt(m,n+1).and.kbt(m,n).eq.kbt(m+1,n+1)
     *		.and.kbt(m,n).eq.kbt(m+1,n)) then
		    continue
c			print  *,'still' , m,n,kzv,u_sum(m,n)
c			read*
		end if


c	fui1(m,n)=(-Ubrt(K2)+U_int(m,N))/dt1
c     fvi1(m,n)=(-Vbrt(K2)+V_int(m,N))/dt1
	ubrt(k2)=u_int(m,n)
	vbrt(k2)=v_int(m,n)
             do k=2,kzv
                dz1=dzk(k)  
c				uu(kp(m,n),k-1)=uu(kp(m,n),k-1)-u_int(m,n)
c				vv(kp(m,n),k-1)=vv(kp(m,n),k-1)-v_int(m,n)
	  
              end do

		
		kzv=kbtv(m,n)
		u_sum(m,n)= uu(kp(m,n),1)* (dzk(2)-lvl_uv)/(h1(m,n)-lvl_uv)      


		do k=3,kzv
			lvl_uv=(lvl(m,n)+lvl(m+1,n)+lvl(m,n+1)+lvl(m+1,n+1))/4.
			dz1=dzk(k)
			u_sum(m,n)=u_sum(m,n)+uu(kp(m,n),k-1)*dz1/(h1(m,n)-lvl_uv)
		end do




		if (u_sum(m,n).gt.1.e-2) then
c			print  *,'achtung!',m,n,kzv,u_sum(m,n)
c			read*
		end if



       end do
      end do






c        Write(nu_diag,*)'end adapt'
 600  FORMAT(2X,11F7.2)
      RETURN
      END