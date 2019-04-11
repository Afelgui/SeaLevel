

      SUBROUTINE new_ad(rho,UU,VV,WW,WWvel,ftok,lvl)
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
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)

          COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
   
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
     *       ,FSU(48),FSV(48)
     *       ,grad_x,grad_y
     *        ,lvl(mh,nh),lvl_uv
      REAL *8 FUH,FVH,FU0,FV0,QC,T1,T2,T3,T4,DET,roz,cfi
      REAL *8 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,hx,hy
      REAL *8 DX0,DX1,DY0,DY1,DZ0,DZ1,ZZ1,ZZ2,XMULT,CW
      real *8 rav(4,kh),den(2,nh,kh),dlam1,dlam2,dtet1,dtet2
     *       ,pc(kh),dcent(kh),df(4,kh),dlam(kh),dtet(kh)
     *       ,plam(kh),ptet(kh),ta,sa,h0,hh1,hh2,fui,fvi,vmod

	
      integer *2 kc(4,kh)
      dimension h4(4)
      integer nu_diag
      common/iodiag/nu_diag
C*-------------------------------------------
C
       ll=0
	cfi=0.55
                         N11=N1-1
                          N21=N2+1
                           KZ1=2
                           KZ2=KZ0-2
                           KKZ1=KKZ-1
                           KKZ2=KKZ-2
C***********************************************
C***********************************************
      DO 10 K=1,kh
   10 UN(K)=1.
c      DO 33 N=1,nh
c   33 SU(N)=AYV(N)/SMB(N)/RAD/RAD
                  do n=1,nh
                 do m=1,mh
                roz(m,n)=0.
               end do
              end do
      ll=0

c	goto 201

      DO 200 m=1,mh-1

	if (m.eq.99) then
	continue
	end if

c      PRINT550,M
c      read *
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

      DO 12 K=KZ1,KKZ
       DO 12 N=1,nh
            den(1,n,k)=rho(kp(m,n),k-1)*0.98
             den(2,n,k)=rho(kp(m+1,n),k-1)*0.98
   12         continue


                 IF(M.EQ.1)GOTO 55
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

	if (ll.eq.178) then
	continue
	end if

c       print *,iydu(1,ll),ll
      IF(ll.le.kyLu.and.IYDu(1,LL).EQ.M)GOTO 100
      LL=LL-1
      GOTO 50
  100 N1=IYDu(2,LL)
      N2=IYDu(3,LL)
      N11=N1-1
      N21=N2+1
      DY1=DY(1,1)
		if (m.eq.60) then
	continue
	end if


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
c               IF(K.EQ.KZV) DZ2=0.5*DZ0
               IF(K.EQ.2)ZZ1=0.
               ZZ2=1
                if(k.eq.kzv) zz2=0
                if(hm.eq.0) then
                print *,m,n
                stop
                end if
   30          HGM(K,N)=hm*dz2
   43       CONTINUE
   42    CONTINUE
C***************************************
C***************************************
      FU1M=0.
      FV1M=0.


      DO 39 N=N1,N2
            hm=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.
      DX1=(dlx(m,n)+dlx(m,n+1))/2.
      DY1=(dly(m,n)+dly(m+1,n))/2.
                   K2=KP(M,N)
                     kzv=kbtv(m,n)
         Z1=iPH(M,N)
          Z2=iPV(M,N)
           Z12=Z1*Z2


            IF(Z12.EQ.0.)then
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
                   HX=H(m+1,n)-H(m,n)+H(m+1,n+1)-H(m,n+1)
                    HY=H(m+1,n+1)-H(m+1,n)+H(m,n+1)-H(m,n)

      fui=0
      fvi=0
        do k=1,kh
        dlam(k)=0
        dtet(k)=0
        end do

      do 57 k=2,kzv
               T1=den(1,N,K)
                T2=den(1,N+1,K)
                 T3=den(2,N,K)
                  T4=den(2,N+1,K)
                   SQ1=IR1(K,N)
                    SQ2=IR2(K,N)
              dlam(k)=(SQ1*(T3-T1)+SQ2*(T4-T2))/dx1/2.
              dtet(k)=(SQ2*(T2-T1)+SQ1*(T4-T3))/dy1/2.
   57  continue


c================================old variant
       plam(2)=0
       ptet(2)=0
        do k=3,kzv

		if (m.eq.42.and.n.eq.64.and.k.eq.3) then
	continue
	end if
                plam(k)=plam(k-1)+
     *                  (dlam(k-1)+dlam(k))*dz(k-1)/2.
                ptet(k)=ptet(k-1)+
     *                   (dtet(k-1)+dtet(k))*dz(k-1)/2.
      
         end do
	         fui=0
      fvi=0

      DO 38 K=KZ1,kzv
               DZ1=DZK(K)
               IF(K.EQ.KZV) DZ1=0.5*DZ(K-1)
       fui=fui-plam(k)*dz1
       fvi=fvi-ptet(k)*dz1
       fu(k)=-plam(k)
       fv(k)=-ptet(k)
   38 CONTINUE


	if (m.eq.60.and.n.eq.60) then
	continue
	end if
c================================ end old variant

c================================new variant
c       plam(1)=0
c       ptet(1)=0
c        do k=2,kzv
c			  plam(k)=plam(k-1)+
c     *                  dlam(k-1)*dzk(k-1)/2.+
c     *                  dlam(k)*dzk(k)/2.
c                ptet(k)=ptet(k-1)+
c     *                  dtet(k-1)*dzk(k-1)/2.+
c     *                  dtet(k)*dzk(k)/2.
      
c         end do
c	fui=0
c      fvi=0
c      DO 138 K=KZ1,kzv
c               DZ1=DZK(K)
c               IF(K.EQ.KZV) DZ1=0.5*DZ(K-1)
c       fui=fui-plam(k-1)*dzk(k-1)/2.-plam(k)*dzk(k)/2.
c       fvi=fvi-ptet(k-1)*dzk(k-1)/2.-ptet(k)*dzk(k)/2.
c       fu(k)=-plam(k)
c       fv(k)=-ptet(k)
c  138 CONTINUE
c================================ end new variant
	lvl_uv=(lvl(m,n)+lvl(m+1,n)+lvl(m,n+1)+lvl(m+1,n+1))/4.D+0
       grad_x(m,n)=fui/(h0-lvl_uv)
       grad_y(m,n)=fvi/(h0-lvl_uv)
c      do k=kz1,kkz
c      fu(k)=fu(k)-fui
c      fv(k)=fv(k)-fvi
c      end do

   62              CW=DX(M,n)*DY(m,n)

    
c((((


      DO 24 K=KZ1,kzv
               DZ1=DZK(K)
c               IF(K.EQ.KZV) DZ1=0.5*DZ(K-1)
		C12(K)=0.
           C11(K,n)=0
            A11(K,N)=0
             A12(K)=0.
              B11(K,n)=HGM(K,N)/dt
              B22(K,n)=B11(k,n)
               B12(K,n)=(-slb(m,n)*hgm(k,n))*cfi
            Fu1(K,n)=Hgm(k,n)/dt*u(k,n)
     *        +FU(K)*hgm(k,n)-(1.-cfi)*(-slb(m,n)*hgm(k,n))*v(k,n)
            Fv1(K,n)=Hgm(k,n)/dt*v(k,n)
     *             +Fv(K)*hgm(k,n)+(1.-cfi)*(-slb(m,n)*hgm(k,n))*u(k,n)

   24    CONTINUE

	if (m.eq.101.and.n.eq.93) then
	continue
	end if

   39 continue
      i1=kz1
c      i2=kzv
      call factuv(a11,a12,c11,c12,b11,b22,b12,fu1,fv1,
     *         	u1,v1,i1,m,n1,n2)
	continue
      go to 101
   50 continue
  660 FORMAT(48(1X,13F6.1,/))
      DO  K=1,kgr
      DO  N=N0,NN
      k2=kp(m,n) 
c===================modificirovAnnaya skorost==============
      UU(K2,K)=U1(K+1,N)        
      VV(K2,K)=V1(K+1,N)


	if (m.eq.99.and.n.eq.79) then
	continue
	end if

      end do
      end do
   55 CONTINUE
  200 CONTINUE
  201 CONTINUE

       do m=1,mh
       do n=1,nh
       uxy(m,n)=0
       vxy(m,n)=0
       end do
       end do

	mtt=141
		ntt=111
	write(2,*), 'after new_ad',uu(kp(mtt,ntt),1), vv(kp(mtt,ntt),1) 



           DO  Kl=1,kxLu
            n=IxDu(1,Kl)
             I1=IxDu(2,Kl)
              I2=IxDu(3,Kl)
           do m=I1,I2

		k2=kp(m,n)


					IF (m.eq.59.and.n.eq.2) then
c	print *,'uu_sum',m,n,k,uu_sum
	continue
c	read*
	end if

             kzv=kbtv(m,n)
	   lvl_uv=(lvl(m,n)+lvl(m+1,n)+lvl(m,n+1)+lvl(m+1,n+1))/4.D+0
	   uxy(m,n)= uu(kp(m,n),1)*(dzk(2)-lvl_uv)/(h1(m,n)-lvl_uv)
         vxy(m,n)= vv(kp(m,n),1)*(dzk(2)-lvl_uv)/(h1(m,n)-lvl_uv)

        
                    do k=3,kzv
                      dz1=dzk(k)
                 uxy(m,n)=uxy(m,n)+uu(kp(m,n),k-1)*dz1/(h1(m,n)-lvl_uv)
                 vxy(m,n)=vxy(m,n)+vv(kp(m,n),k-1)*dz1/(h1(m,n)-lvl_uv)
c	 uxy(m,n)=uxy(m,n)+uu(kp(m,n),k-1)*dz1/h1(m,n)
c                 vxy(m,n)=vxy(m,n)+vv(kp(m,n),k-1)*dz1/h1(m,n)
                    end do
		

					IF (m.eq.59.and.n.eq.2) then
c	print *,'uu_sum',m,n,k,uu_sum
	continue
c	read*
	end if

	      do k=2,kzv
                	uu(kp(m,n),k-1)=uu(kp(m,n),k-1)-uxy(m,n)
				vv(kp(m,n),k-1)=vv(kp(m,n),k-1)-vxy(m,n)				
              end do

	k2=kp(m,n)
				IF (m.eq.59.and.n.eq.2) then
c	print *,'uu_sum',m,n,k,uu_sum
	continue
c	read*
	end if


c==========================barotropnaya=================	
c	     ubrt(k2)=uxy(m,n)
c		 vbrt(k2)=vxy(m,n)  
           end do
           end do


  600 FORMAT(2X,11F7.2)
      RETURN
      END
