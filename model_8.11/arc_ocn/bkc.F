      SUBROUTINE BKC(T,S,UW,VW)
      Include 'model.par'
      COMMON/DAY/myear,month,nday
      COMMON/NH2/ NH2
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      real *8 z,dz,dzk,dx,dy,h,h1
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      COMMON/botV/h1(mh,nh)
      COMMON/botTS/h(mh,nh)
      real *8 tflux,sflux,tax,tay
      common/Wind/Tax(mh,nh),Tay(mh,nh)    
      real *8 sl,slb
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),slb(mh,nh)
#ifdef TRACERS
      real *8 Ptrs(nodes,kgr,ntrac)
      COMMON/Trac/Ptrs
      real *8 Ptr1(kh,ntrac),Ptr2(kh,nh,ntrac),Ptrx(111,ntrac)
      Common/TracBKC/Ptr1,Ptrx
#endif

      common/zbks/zzbks(kgr)

      real *8 T(nodes,kgr),S(nodes,kgr)
     *       ,UW(nodes,kgr),VW(nodes,kgr)
      COMMON/BL1/Ihbl(nodes)
      REAL*8 G,PI,DX1,DX2,DY1,DY2,A,A1,A2,A3,A0,HH,HF,ROSSBY,CORR
      REAL*8 U0(kh,nh),V0(kh,nh),U1(kh,nh),V1(kh,nh),
     *UU(kh),VV(kh),TT(kh),SS(kh),RR(kh),
     *UU0(kh,nh),VV0(kh,nh),T0(kh,nh),S0(kh,nh),
     *PART
      COMMON/MixB/G,ROSSBY,CORR,PART,KZT,KZU
      integer nu_diag
      common/iodiag/nu_diag
c      RS(AA,BB)=28.152-AA*(0.0735+0.00469*AA)+
c     *(0.802-0.002*AA)*(BB-35.)
      PART=0.9
C************************************************
      do k=1,kgr
      zzbks(k)=z(k+1)
      end do
      PI=3.14159265
      G=980.
      KZT=kgr
      KZU=kgr
      IH=1
      IK=mh
      IE=IK-1
      LL=0
       do k=1,kh
       do n=1,nh
         u0(k,n)=0
         v0(k,n)=0
         end do
         end do
C************************************************
      DO 1 I=IH,IE
            DO 10 J=1,nh
               K0=KP(I+1,J)
            DO 10 K=2,kgr+1
               T0(K,J)=T(K0,K-1)
               S0(K,J)=S(K0,K-1)
#ifdef TRACERS
               Do nt=1,ntrac
                  Ptr2(K,J,nt)=Ptrs(K0,K-1,nt)
               End Do
#endif
               U1(K,J)=UW(K0,K-1)
               V1(K,J)=VW(K0,K-1)
                  UU0(K,J)=U1(K,J)
   10                VV0(K,J)=V1(K,J)
         IT=I+1
 2003    LL=LL+1
         IF(IYD(1,LL).EQ.IT.and.ll.le.kyLine)GO TO 2004
         LL=LL-1
         GO TO 3
 2004 JH=IYD(2,LL)
         JE=IYD(3,LL)
C*    PRINT 4000,I,LL
 4000 FORMAT(1X,' ** ',I3,' ** ',I3)
         DO 2 J=JH,JE
            K0=KP(IT,J)
                        do k=1,kh
                        uu(k)=0
                        vv(k)=0
                        tt(k)=0
                        ss(k)=0
#ifdef TRACERS
                           do nt=1,ntrac
                              Ptr1(k,nt)=0
                           end do
#endif
                        end do
            FACTOR=-1.
  337       CONTINUE
            IF(FACTOR)338,339,2
  338       IF(K0.EQ.1)GO TO 340
            DX1=DX(it-1,j)
            DX2=DX(it,j)
            DY1=DY(it,j-1)
            DY2=DY(it,j)
            A=DX2*DY2+DX2*DY1+DX1*DY2+DX1*DY1
      A=1./A
            A1=A*DX1*DY1
            A2=A*DX1*DY2
            A3=A*DX2*DY1
            A0=A*DX2*DY2
            kf=kbt(it,j)-1
            DO 11 K=1,Kf
               KK=K+1
               UU(K)=A1*U0(KK,J-1)+A2*U0(KK,J)+A3*U1(KK,J-1)+
     *               A0*U1(KK,J)
               VV(K)=A1*V0(KK,J-1)+A2*V0(KK,J)+A3*V1(KK,J-1)+
     *               A0*V1(KK,J)
   11       CONTINUE
            CORR=SL(it,J)
            DO 12 K=1,Kf
               TT(K)=T(K0,K)
               SS(K)=S(K0,K)
#ifdef TRACERS
               do nt=1,ntrac
                  Ptr1(k,nt)=Ptrs(K0,K,nt)
               end do
#endif
   12       CONTINUE
            HH=H(IT,J)
            K1=KP(IT-1,J)
            K2=KP(IT,J-1)
            K3=KP(IT-1,J-1)
         TX=0.25*(TAX(IT,J)+TAX(IT-1,J)+TAX(IT,J-1)+TAX(IT-1,J-1))
         TY=0.25*(TAY(IT,J)+TAY(IT-1,J)+TAY(IT,J-1)+TAY(IT-1,J-1))
            GO TO 341
  339       IF((iPV(it,j)*iPH(it,j)).EQ.0.)GO TO 2
            K1=KP(IT,J+1)
            K2=KP(IT+1,J)
            K3=KP(IT+1,J+1)
            CORR=SLb(it,J)
            kf=kbtv(it,j)-1
            DO 13 K=1,kf
               UU(K)=U1(K+1,J)
               VV(K)=V1(K+1,J)
   13       CONTINUE
            DO 14 K=1,Kf
               te=0.
               se=0.
               pe=0.
               if(kbt(it,j)-1.ge.k)then
                  te=te+t(k0,k)
                  se=se+s(k0,k)
                  pe=pe+1.0
               endif
               if(kbt(it,j+1)-1.ge.k)then
                  te=te+t(k1,k)
                  se=se+s(k1,k)
                  pe=pe+1.0
               endif
               if(kbt(it+1,j)-1.ge.k)then
                  te=te+t(k2,k)
                  se=se+s(k2,k)
                  pe=pe+1.0
               endif
               if(kbt(it+1,j+1)-1.ge.k)then
                  te=te+t(k3,k)
                  se=se+s(k3,k)
                  pe=pe+1.0
               endif
               if(pe.eq.0)then
                  print *,'????????????????',it,j
               endif
               tt(k)=te/pe
               ss(k)=se/pe
c               TT(K)=0.25*(T(K0,K)+T(K1,K)+T(K2,K)+T(K3,K))
c               SS(K)=0.25*(S(K0,K)+S(K1,K)+S(K2,K)+S(K3,K))
#ifdef TRACERS
               do nt=1,ntrac
                  Ptr1(k,nt)=0.
               end do
#endif
   14       CONTINUE
            HH=H1(it,j)
            TX=TAX(IT,J)
            TY=TAY(IT,J)
  341       CONTINUE
         ROSSBY=(TX**2+TY**2)**0.25
C        IF(J.NE.28)GO TO 500
C        IF(I.NE.33)GO TO 500
C        DO 501 K=1,12
C 501    RR(K)=RS(TT(K),SS(K))
C        PRINT 502,(TT(K),SS(K),RR(K),K=1,12)
  502 FORMAT(3(3X,12(F5.2,2X),/))
  500 CONTINUE
      If(KF.GE.2) CALL ACOCLA(TT,SS,UU,VV,HH,KF,HF,KDEPTH,FACTOR)
  510 CONTINUE
            IF(FACTOR)30,31,2
   30    CONTINUE
C         IF(iHbl(K0).LT.KDEPTH)  iHbl(K0)=KDEPTH
         KDEPTH=0
         DO k=1,kgr
            if(0.01*ZZBKS(k).GT.HF)exit
            KDEPTH=k+1
         END DO
         IF(iHbl(K0).LT.KDEPTH)  iHbl(K0)=KDEPTH
         HF=0.
            DO 32 K=1,Kf
#ifdef TRACERS
               do nt=1,ntrac
                  Ptr2(k+1,j,nt)=Ptr1(k,nt)
               end do
#endif
               T0(K+1,J)=TT(K)
   32          S0(K+1,J)=SS(K)
  340       FACTOR=FACTOR+1.
            GO TO 337
   31    CONTINUE
         HF=0.
            DO 33 K=1,Kf
               UU0(K+1,J)=UU(K)
   33          VV0(K+1,J)=VV(K)
    2    CONTINUE
         GO TO 2003
C************************************************
C************************************************
    3    DO 3003 J=2,nh-1
            K0=KP(IT,J)
            DO 3003 K=1,kgr
               T(K0,K)=T0(K+1,J)
               S(K0,K)=S0(K+1,J)
              UW(K0,K)=UU0(K+1,J)
              VW(K0,K)=VV0(K+1,J)
#ifdef TRACERS
               do nt=1,ntrac
                  Ptrs(K0,K,nt)=Ptr2(K+1,J,nt)
               end do
#endif
 3003 CONTINUE
         DO 3004 J=2,nh-1
         DO 3004 K=2,kh-1
            U0(K,J)=U1(K,J)
 3004       V0(K,J)=V1(K,J)
    1 CONTINUE
C******************************************
C******************************************
 3005 CONTINUE
C*    PRINT 4001
 4001 FORMAT(1X,' ALL RIGHT')
      IF(NDAY/30*30.NE.NDAY)GOTO 1223
c      DO 1221 I=1,mh
cc      PRINT 1222,(iHbl(KP(I,J)),J=5,44)
 1222 FORMAT(1X,40I3)
 1221 CONTINUE
 1223 CONTINUE
C************************************************
      RETURN
      END
      SUBROUTINE ACOCLA(TM,SM,UU,VV,H,KF,HHH,KDEPTH,FACTOR)
      Include 'model.par'
      common/zbks/zzbks(kgr)
      REAL*8 TM(kh),SM(kh),UU(kh),VV(kh),HHZ(111),
     *HZU(kh),HZ(kh),T(111),S(111),U(111),V(111),R(111),
     *R1(111),UN(111),VN(111),F(111),DZ(111),HZN(kd),
     *G,RICR,H,HHH,AU,AV,DDZ,SAC,DU,DV,DR,DIZ,AAS,BBS,CCS,ALS,BTS,
     *AI1,AIS1,AI2,AI3,AI4,AI5,AIS5,AI6,WVS,WUS,UUS,VVS,SQARU,
     *SQARV,ROSSBY,CORR,PART,RS,aa,bb
      COMMON/MixB/ G,ROSSBY,CORR,PART,NZT,NZU
#ifdef TRACERS
      real *8 TAIS1(ntrac),TAIS5(ntrac)
      real *8 Ptr1(kh,ntrac),Ptrx(111,ntrac)
      Common/TracBKC/Ptr1,Ptrx
#endif
      integer nu_diag
      common/iodiag/nu_diag
      logical fine
c      DATA HZN/0.,30.,100.,300.,500.,700.,1000.,1500.,2500.,
c     *3500.,5000.,0./
c      RS(AA,BB)=1.+(28.152-AA*(0.0735+0.00469*AA)+
c     *(0.802-0.002*AA)*(BB-35.))*0.001
         fine=.false.
         RICR=1.
         HZ(1)=0.
       do k=1,kgr
       hz(k)=zzbks(k)
       end do
       hz(kh-1)=h
       hz(kh)=h
         NN=kf
         R(1)=(dens(TM(1),SM(1),0.)+1025.)/1000.
         AU=UU(1)
         AV=VV(1)
      DO 610 K=2,NN
         KK=K
         DDZ=HZ(K)-HZ(K-1)
         SAC=0.
         IF(K.NE.NN)SAC=0.5*(HZ(K+1)-HZ(K))
         DU=(UU(K)-AU)/DDZ
         DV=(VV(K)-AV)/DDZ
         DU=DU**2+DV**2
         R(K)=(1025.+dens(TM(K),SM(K),0.))/1000.
         DR=2*(R(K)-R(1))/(DDZ*(R(K)+R(1)))
         R1(K)=G*DR-RICR*DU
         IF(R1(K))609,609,611
  609    AU=(AU*HZ(K-1)+0.5*DDZ*(AU+UU(K))+SAC*UU(K))/(HZ(K)+SAC)
         AV=(AV*HZ(K-1)+0.5*DDZ*(AV+VV(K))+SAC*VV(K))/(HZ(K)+SAC)
  610 CONTINUE
         KK=NN
  611    KE=KK-1
         DIZ=DDZ/10.
c      PRINT 703,(TM(K),SM(K),UU(K),VV(K),K=1,12)
      if(.not.fine)ke=nn
      DO 606 K=1,KE
      U(K)=UU(K)
      V(K)=VV(K)
#ifdef TRACERS
         do nt=1,ntrac
            Ptrx(k,nt)=Ptr1(k,nt)
         end do
#endif
      T(K)=TM(K)
      S(K)=SM(K)
  606 HHZ(K)=HZ(K)
      if(.not.fine) goto 6019
      DO 607 KK=1,9
      K=KK+KE
      U(K)=UU(KE)+KK*DIZ*(UU(KE+1)-UU(KE))/DDZ
      V(K)=VV(KE)+KK*DIZ*(VV(KE+1)-VV(KE))/DDZ
#ifdef TRACERS
         do nt=1,ntrac
            Ptrx(k,nt)=Ptr1(ke,nt)+
     *         kk*diz*(Ptr1(ke+1,nt)-Ptr1(ke,nt))/DDZ
         end do
#endif
      T(K)=TM(KE)+KK*DIZ*(TM(KE+1)-TM(KE))/DDZ
      S(K)=SM(KE)+KK*DIZ*(SM(KE+1)-SM(KE))/DDZ
  607 HHZ(K)=HHZ(K-1)+DIZ
      KE1=KE+1
      DO 608 K=KE1,NN
      KI=K+9
      U(KI)=UU(K)
      V(KI)=VV(K)
#ifdef TRACERS
         do nt=1,ntrac
            Ptrx(ki,nt)=Ptr1(k,nt)
         end do
#endif
      T(KI)=TM(K)
      S(KI)=SM(K)
  608 HHZ(KI)=HZ(K)
      N=9+NN
 6019 if(.not.fine)n=nn
      N1=N-1
      DO 11 I=2,N
   11 DZ(I)=HHZ(I)-HHZ(I-1)
      DZ(1)=0.
  200 DO 51 I=1,N
   51 R(I)=(1025.+dens(T(I),S(I),0.))/1000.
      NUF=1
      CALL MIX(U,V,T,S,R,N,DZ,HHZ,NUF)
      HHH=HHZ(NUF)/100.
      KDEPTH=NUF
      IF(KDEPTH.EQ.1)KDEPTH=0
      if(.not.fine)then
      do k=1,n
         uu(k)=u(k)
         vv(k)=v(k)
#ifdef TRACERS
            do nt=1,ntrac
               ptr1(k,nt)=ptrx(k,nt)
            end do
#endif
         tm(k)=t(k)
         sm(k)=s(k)
      end do
      return
      endif
      UU(1)=U(1)
      VV(1)=V(1)
#ifdef TRACERS
      do nt=1,ntrac
         Ptr1(1,nt)=Ptrx(1,nt)
      end do
#endif
      TM(1)=T(1)
      SM(1)=S(1)
      IF(KE.LE.2)GO TO 707
      KE1=KE-1
      DO 612 K=2,KE1
      UU(K)=U(K)
      VV(K)=V(K)
#ifdef TRACERS
         do nt=1,ntrac
            Ptr1(k,nt)=Ptrx(k,nt)
         end do
#endif
      TM(K)=T(K)
  612 SM(K)=S(K)
  707 AAS=0.
      BBS=0.
      IF(KE.GE.2)AAS=0.5*10*(HZ(KE)-HZ(KE-1))/(HZ(KE+1)-HZ(KE))
      IF(KE.LE.(NN-2))BBS=0.5*10*(HZ(KE+2)-HZ(KE+1))/(HZ(KE+1)-HZ(KE))
      CCS=0.5*10
      ALS=0.5+AAS
      BTS=0.5+BBS
      AI1=ALS*T(KE)+BTS*T(KE+10)
#ifdef TRACERS
      do nt=1,ntrac
         TAIS1(nt)=ALS*Ptrx(KE,nt)+BTS*Ptrx(KE+10,nt)
         TAIS5(nt)=ALS*Ptrx(KE,nt)*HHZ(KE)+BTS*Ptrx(KE+10,nt)*HHZ(KE+10)
      end do
#endif
      AIS1=ALS*S(KE)+BTS*S(KE+10)
      AI2=ALS*U(KE)+BTS*U(KE+10)
      AI3=ALS*V(KE)+BTS*V(KE+10)
      AIU=ALS*U(KE)**2+BTS*U(KE+10)**2
      AIV=ALS*V(KE)**2+BTS*V(KE+10)**2
      AI5=ALS*T(KE)*HHZ(KE)+BTS*T(KE+10)*HHZ(KE+10)
      AIS5=ALS*S(KE)*HHZ(KE)+BTS*S(KE+10)*HHZ(KE+10)
      AI6=ALS*T(KE)*V(KE)+BTS*T(KE+10)*V(KE+10)
      KE1=KE+1
      KE9=KE+9
      DO 701 J=KE1,KE9
         AI1=AI1+T(J)
#ifdef TRACERS
         do nt=1,ntrac
            TAIS1(nt)=TAIS1(nt)+Ptrx(J,nt)
            TAIS5(nt)=TAIS5(nt)+Ptrx(J,nt)*HHZ(J)
         end do
#endif
         AIS1=AIS1+S(J)
         AI2=AI2+U(J)
         AI3=AI3+V(J)
         AIU=AIU+U(J)**2
         AIV=AIV+V(J)**2
         AI5=AI5+T(J)*HHZ(J)
         AIS5=AIS5+S(J)*HHZ(J)
  701    AI6=AI6+T(J)*V(J)
C      print *,701
C      read *
      IF(KE.EQ.1)GO TO 704
C*$$$ TM(KE)=(AI1*HZ(KE+1)-AI5)/((AAS+CCS)*(HZ(KE+1)-HZ(KE)))
C*$$$ SM(KE)=(AIS1*HZ(KE+1)-AIS5)/((AAS+CCS)*(HZ(KE+1)-HZ(KE)))
C*$$$ TM(KE+1)=(AI5-AI1*HZ(KE))/((BBS+CCS)*(HZ(KE+1)-HZ(KE)))
C*$$$ SM(KE+1)=(AIS5-AIS1*HZ(KE))/((BBS+CCS)*(HZ(KE+1)-HZ(KE)))
      TM(KE)=T(KE)
      TM(KE+1)=T(KE+10)
#ifdef TRACERS
      do nt=1,ntrac
         Ptr1(KE,nt)=Ptrx(KE,nt)
         Ptr1(KE+1,nt)=Ptrx(KE+10,nt)
      end do
#endif
      SM(KE)=S(KE)
      SM(KE+1)=S(KE+10)
C*    IF(DABS(TM(KE)-TM(KE+1)).LE.1.E-6)GO TO 702
      WVS=1.
      WUS=1.
      IF(U(KE).GE.U(KE+10))WUS=-1.
      IF(V(KE).GE.V(KE+10))WVS=-1.
      UUS=AAS+BBS+2.*CCS
      VVS=(BBS+CCS)/(AAS+CCS)
      SQARU=UUS*AIU-AI2**2
      SQARV=UUS*AIV-AI3**2
C*    PRINT 803,SQARU,SQARV
  803 FORMAT(25X,E12.6,E12.6)
      IF(SQARU.LT.0.)SQARU=0.
      IF(SQARV.LT.0.)SQARV=0.
 
c        The old way to calculate UU(KE), VV(KE)
c        UU(KE+1), VV(KE+1) preserving vertical integral
c        of u, v, u^2 and v^2 returned back because
c        otherwise vertical velosity at the bottom
c        becomes different from zero.

      UU(KE)=(AI2-WUS*DSQRT(SQARU*VVS))/UUS
      VV(KE)=(AI3-WVS*DSQRT(SQARV*VVS))/UUS
      UU(KE+1)=(AI2+WUS*DSQRT(SQARU/VVS))/UUS
      VV(KE+1)=(AI3+WVS*DSQRT(SQARV/VVS))/UUS
c      UU(KE)=U(KE)
c      UU(KE+1)=U(KE+10)
c      VV(KE)=V(KE)
c      VV(KE+1)=V(KE+10)
      GO TO 705
  702 VV(KE)=V(KE)
      UU(KE)=U(KE)
      UU(KE+1)=U(KE)
      VV(KE+1)=V(KE)
      GO TO 705
  704 CONTINUE
C*    TM(2)=(AI1-TM(1)*(AAS+CCS))/(BBS+CCS)
C*    SM(2)=(AIS1-SM(1)*(AAS+CCS))/(BBS+CCS)
C*    UU(2)=(AI2-UU(1)*(AAS+CCS))/(BBS+CCS)
C*    VV(2)=(AI3-VV(1)*(AAS+CCS))/(BBS+CCS)
      TM(2)=T(11)
      SM(2)=S(11)
#ifdef TRACERS
      do nt=1,ntrac
         Ptr1(2,nt)=Ptrx(11,nt)
      end do
#endif
      UU(2)=U(11)
      VV(2)=V(11)
  705 CONTINUE
      IF(KE.GT.(NN-2))GO TO 801
      KE2=KE+2
      DO 614 K=KE2,NN
      UU(K)=U(K+9)
      VV(K)=V(K+9)
#ifdef TRACERS
      do nt=1,ntrac
         Ptr1(K,nt)=Ptrx(K+9,nt)
      end do
#endif
      TM(K)=T(K+9)
  614 SM(K)=S(K+9)
  801 CONTINUE
c      PRINT 703,(TM(K),SM(K),UU(K),VV(K),K=1,12)
c      READ *
  703 FORMAT(1X,'*',/,12(4(E12.5,3X),/))
  616 CONTINUE
  620 CONTINUE
C      print *, 'al is ok'
C      read *
   33 RETURN
      END
c      SUBROUTINE MIX(U,V,T,S,R,N,DZ,HHZ,NUF)
      SUBROUTINE MIX(U,V,T,S,R,N,DZ,HHZ,NUF)
      Include 'model.par'

c      New way of mean U and V calculation is introduced
c      by means of UAR and VAR variables in order to
c      conserve zero value of vertical velocity component
c      at the bottom

      REAL*8
C      Real
     *  U(*),V(*),T(*),S(*),R(*),RI(251),RU(251)
     *,DZ(*),HHZ(*),PART,EP,QUOTE
     *,G,RA,TA,SA,UA,VA,UAR,VAR,ZA,E,AMS,R1,R2,ROSSBY,FC,EED,EG
      COMMON/MixB/G,ROSSBY,FC,PART,KZU,KZT
      integer nu_diag
      common/iodiag/nu_diag
#ifdef TRACERS
      real *8 PA(ntrac), PAo(ntrac)
      real *8 Ptr1(kh,ntrac),Ptrx(111,ntrac)
      Common/TracBKC/Ptr1,Ptrx
#endif
c      RS(AA,BB)=1.+(28.152-AA*(0.0735+0.00469*AA)+
c     *(0.802-0.002*AA)*(BB-35.))*0.001
      nuf=1
      jstop=0
      If(Rossby.eq.0.)Rossby=0.000001
      QUOTE=1.0-PART
      HHtot=0.
c      PRINT 703,(R(K)-1.025,R(K)-1.025,U(K),V(K),K=1,12)
c      READ *
  703 FORMAT(1X,'*',/,12(4(E12.5,3X),/))
      J=0
      II=0
      K=1
    1 K=K+1
      IF(K.GT.N)Then
         IF(J.EQ.0)GO TO 1111
         IM=N
         GOTO 4
      END IF
         RI(K)=G*DZ(K)*(R(K)-R(K-1))
         RU(K)=R(K)*R(K-1)*((U(K)-U(K-1))**2+(V(K)-V(K-1))**2)/
     *         (R(K)+R(K-1))
c         Print *, k, ri(k), ru(k), g
c         Read *
            IF(J.EQ.0)GO TO 3
         IF(RI(K).LE.RU(K))GO TO 1
         IM=K-1
         GO TO 4
    3    IF(RI(K).GT.RU(K))GO TO 1
         IN=K-1
         J=1
         GO TO 1
    4    NI=0
         MI=0
         IIo=IN
         JJo=IM
    8    RA=0.
         TA=0.
         SA=0.
         UA=0.
         VA=0.
         UAR=0.
         VAR=0.
#ifdef TRACERS
         do nt=1,ntrac
            PA(nt)=0.
         enddo
#endif
         ZA=0.
         II=IN-NI+1
         JJ=IM+MI
c         Print *,ii,jj
c         read *
         II1=II-1
         E=0.
         DO 2 L=II1,JJ
         RA=RA+R(L)*0.5*(DZ(L)+DZ(L+1))
         TA=TA+T(L)*0.5*(DZ(L)+DZ(L+1))
         SA=SA+S(L)*0.5*(DZ(L)+DZ(L+1))
#ifdef TRACERS
         do nt=1,ntrac
            PA(nt)=PA(nt)+Ptrx(L,nt)*0.5*(DZ(L)+DZ(L+1))
         enddo
#endif
         UA=UA+U(L)*R(L)*0.5*(DZ(L)+DZ(L+1))
         VA=VA+V(L)*R(L)*0.5*(DZ(L)+DZ(L+1))
         UAR=UAR+U(L)*0.5*(DZ(L)+DZ(L+1))
         VAR=VAR+V(L)*0.5*(DZ(L)+DZ(L+1))
         ZA=ZA+0.5*(DZ(L)+DZ(L+1))
         E=E+R(L)*(U(L)**2+V(L)**2)*0.5*(DZ(L)+DZ(L+1))
    2    CONTINUE
         AMS=1./ZA
         TA=TA*AMS
         SA=SA*AMS
         UAR=UAR*AMS
         VAR=VAR*AMS
#ifdef TRACERS
         do nt=1,ntrac
            PA(nt)=PA(nt)*AMS
         enddo
#endif
C-----------------------------------------------------
Comment out next two lines if you don't want Direchlet
cccc         IF(II.EQ.2)TA=T(1)
cccc         IF(II.EQ.2)SA=S(1)
C-----------------------------------------------------
         UA=UA/RA
         VA=VA/RA
         RA=(1025.+dens(TA,SA,0))/1000.
         E=E-(UA**2+VA**2)*RA*ZA
         EED=0.87*E*(1.+0.04*HHZ(JJ)*abs(FC)/ROSSBY)
         IF((E-EED).LT.0.
c     *  .OR.ii1.ne.1
     *   )EED=E
         E=E-EED
         EP=0.
         DO 5 L=II1,JJ
            QUOTE=1.
            EG=G*(DZ(L)+DZ(L+1))*(RA-R(L))*HHZ(L)
            IF(EP.LE.0.)GO TO 21
            QUOTE=exp(-7.631*DZ(L)/(HHZ(L-1)+HHZ(L)))
   21       EP=EP*QUOTE+EG
    5    CONTINUE
         E=E+EP
c sledujuchee dobavleno, chtoby iskljuchit' vsjakoe vovlechenie
         e=0.0
c konec dobavlenija
         IF(ni+mi.eq.0)Then
            RAo=RA
            TAo=TA
            SAo=SA
#ifdef TRACERS
            do nt=1,ntrac
               PAo(nt)=PA(nt)
            enddo
#endif
            UAo=UA
            VAo=VA
            UARo=UAR
            VARo=VAR
         End if
         R1=1.
         IF(II.NE.2)R1=RA-R(II-2)
         R2=1.
         IF(JJ.NE.N)R2=R(JJ+1)-RA
         IF(E)7,7,10
   10    IF((II.NE.2).AND.(((NI*R1).LE.(MI*R2)).OR.(JJ.EQ.N)))NI=NI+1
         IF((JJ.NE.N).AND.(((NI*R1).GE.(MI*R2)).OR.(II.EQ.2)))MI=MI+1
         IF((JJ.EQ.N).AND.(II.EQ.2))GO TO 11
         iio=ii1
         jjo=jj
         RAo=RA
         TAo=TA
         SAo=SA
         UAo=UA
         VAo=VA
         UARo=UAR
         VARo=VAR
         GO TO 8
   11    JSTOP=1
    7    DO 9 L=IIo,JJo
c            U(L)=UAo
c            V(L)=VAo
            U(L)=UARo
            V(L)=VARo
            T(L)=TAo
            S(L)=SAo
#ifdef TRACERS
            do nt=1,ntrac
               Ptrx(L,nt)=PAo(nt)
            enddo
#endif
            R(L)=RAo
    9    CONTINUE
         IF((II.EQ.2).AND.(JJ.GT.1))Then
                NUF=JJ
                hhtot=hhz(jj)*0.01
         Else
                hhtot=hhtot+(hhz(jj)-hhz(ii1))*0.01
         End If
         IF(JJ.EQ.N)GO TO 12
         J=0
         IF(JSTOP.EQ.1)GO TO 12
      K=JJ
      GO TO 1
 1111 CONTINUE
   12 IF((II.EQ.2).AND.(JJ.GT.1))NUF=JJ
C      print *,'mix'
C      read *
      RETURN
      END
