        SUBROUTINE tras_difXZ(Sal)
C**********************************************
C*                                            *
C*    CROSXZ:                                 *
C*                RU(T)=10*RU(S)              *
C**********************************************
      Include 'model.par'
      real *8 Sal(nodes,kgr)
      common/flux/tflux(mh,nh),sflux(mh,nh)
      real *8 tflux,sflux
      common/restor/Sstar(mh,nh),gsalt
      real *8 Sstar,gsalt
      COMMON/DAY/myear,month,nday
      COMMON/NH2/ NH2
      real *8 z,dz,dzk,dx,dy,h,h1
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/clame/hl1(mh,nh),hl2(mh,nh)
      common/Cfsmag/cfsmag(mh,nh,kgr)
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      COMMON/d1/ Rad,Gro,R,dt
c      Real *8 SL,SM,SMB,SLB,DX,DY,Z,DZ,DZK
      COMMON/RUNPARS/nmyear

      real *8 A1(1000),B1(1000),C1(1000),F1(1000),U1(1000),
     *        FS1(1000),US1(1000)

      REAL *8 TS(mh,kh),SS(mh,kh),umx1(mh,kh),rivx(mh,nh)
     *       ,ts1(mh,kh),SS1(mh,kh) ,ts2(mh,kh),SS2(mh,kh)
     *        ,gk1,gk2
     *        ,ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
      REAL *8 RZ1(mh,kh),RX(mh,kh)
      REAL *8 CZ(mh,kh),CX(mh,kh),ru_new(nh,kh)
      REAL *8 HGM(mh,kh),czs(mh,kh),hgm0(mh),ghan,cnr

      REAL *8 DZ0,DZ1,Axx,Bxx,Cxx,tadd,sadd,dzgam
      integer nu_diag
      common/iodiag/nu_diag
      COMMON/NU/RU1(nh)
      if(myear.le.nmyear)then
         tt=(myear-nmyear+1)*360+(month-1)*30+min(nday,30)-1
         tt=0.5*(1.0+cos(tt*3.14159265/720.0))
      else
         tt=0.0
      endif
            abh=0.5e+6*(1.-tt)+1.0e+6*tt
            ash=0.5e+6*(1.-tt)+1.0e+6*tt
       cnr=0.5
c	 if(nday/5*5.eq.nday) cnr=1
      KZ1=2
      NZV=0
      LL=0
      NZL=0
       astar=0.8
c       astar=0.1
	 arr=1.05/3.14159265
      alam=4.5e-5
      zstar=250000

      DO k=2,kkz
	do n=1,95
c         RU(N)=1.
      ru_new(n,k)=astar+arr*atan(alam*(z(k)-zstar))
	end do
	do n=1,nh
	ru_new(n,k)=0.05
	end do
      end do 
C*------------------------------------------------
      DO 73 N=N0,nh-1
c       Write(nu_diag,*)n
c       read *
       DO 222 K=1,kh
      DO 222 M=1,mh
             CZ(M,K)=0.
             CZS(M,K)=0.
             CX(M,K)=0.
             HGM(M,K)=0.
  222 continue       
      l=1
C*------------------------------------------------
      NZL=NZL+1
      DO 81 M=1,mh
      K2=KP(M,N)
      DY0=DY(m,N-1)
      DY1=DY(m,N)
      DO 81 K=2,KKZ
            SS(M,K)=SAL(K2,K-1)
                       ss1(m,k)=ss(m,k)
                         ss2(m,k)=ss(m,k)
               HGM(M,K)=2./DT*DlX(M,n)*DlY(m,n)*DZ(K)
   81 CONTINUE
      DO 74 K=2,kkz
      DO 74 M=M0,MM
            umx1(m,k)=abh+(ash-abh)*exp(-z(k)/50000.)
c            umx1(m,k)=1.e+9

   74 CONTINUE
C*************************************************
      LK=1
      NCROSS=0
    2 LL=LL+1
      X3=0.
      IF(ll.le.kxLine.and.IXD(1,LL).EQ.N)GOTO 3
c      Write(nu_diag,*)ll,n,ixd(1,ll)
c      read *
      LL=LL-1
      GOTO 201
    3 continue
      M1=IXD(2,LL)
      M2=IXD(3,LL)
         NCROSS=NCROSS+1
          LN=LK+1
           NRA=0
            MCc1=M1
             MCc2=M2
         lk=l
      DO 20 K=2,KKZ
      DO 5 M=MCc1,Mcc2

           RZ1(M,K)=ru_new(n,k)/2.
           RX(M,K)=UMX1(M,k)+cfsmag(m,n,k-1)*0
      if(n.lt.4)rx(m,k)=5.e+9    
    5 CONTINUE
   20 CONTINUE
       
C*****************************
C*****************************
      DO 40 M=M1,M2
       ghan=(dz(2)+dz(3)*0.5)/(30.*86400.0)
       KLAM=KBT(M,N)
	k2=kp(m,n)
              ds1=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.
           ds2=(dlx(m-1,n)+dlx(m-1,n+1))*(dly(m,n)+dly(m-1,n))/4.
           ds3=(dlx(m-1,n)+dlx(m-1,n-1))*(dly(m,n-1)+dly(m-1,n-1))/4.
           ds4=(dlx(m,n)+dlx(m,n-1))*(dly(m,n-1)+dly(m+1,n-1))/4.
         DO 41 K=KZ1,KLAM
             DZ0=DZ(K-1)
             DZ1=DZ(K)
               top1=ikbt(m,n,k)
               top2=ikbt(m-1,n,k)
               top3=ikbt(m-1,n-1,k) 
               top4=ikbt(m,n-1,k)
                     bot1=ikbt(m,n,k+1)
                     bot2=ikbt(m-1,n,k+1)
                     bot3=ikbt(m-1,n-1,k+1) 
                     bot4=ikbt(m,n-1,k+1)
	top1=1
	top2=1
	top3=1
	top4=1
	bot1=1
	bot2=1
	bot3=1
	bot4=1

           Gam=((ds1*top1+ds2*top2+ds3*top3+ds4*top4)*dz0
     *      +(ds1*bot1+ds2*bot2+ds3*bot3+ds4*bot4)*dz1)/8. 
       dzgam=(dz0+dz1)/2.
          if(k.eq.klam)dzgam=dz0/2.
          if(gam.eq.0) gam=ds1*dzgam
            zs1=1
          if(k.eq.2)then
               zs1=0
           hgm0(m)=(ds1*top1+ds2*top2+ds3*top3+ds4*top4)/4.*3.
             
          end if
      HGM(M,K)=6.*GAM*2.0



    

      CX(M,K)=-DZ0*(RX(M,K-1)+RX(M,K))*1.5
     *        *0.5*((dly(m,n)+dly(m+1,n))*0.5/dlx(m,n)*top1
     *             +(dly(m,n-1)+dly(m+1,n-1))*0.5/dlx(m,n)*top4)
     *        -DZ1*(RX(M,K+1)+RX(M,K))*1.5
     *        *0.5*((dly(m,n)+dly(m+1,n))*0.5/dlx(m,n)*bot1
     *             +(dly(m,n-1)+dly(m+1,n-1))*0.5/dlx(m,n)*bot4)

      cx(m,k)=cx(m,k)*inbt(m+1,n,k)


      CZs(M,K)=-((ds1*bot1+ds4*bot4)*0.5*(RZ1(M,K)+RZ1(M,K+1))*1.5
     *         +(ds2*bot2+ds3*bot3)*0.5*(RZ1(M,K)+RZ1(M,K+1))*1.5)
     *         /DZ1

   41 CONTINUE
       cz(m,klam)=0.
       czs(m,klam)=0.
   40 CONTINUE

          do k=1,kkz
            cx(m1-1,k)=0
            cx(m2,k)=0
            end do
C*****************************


   16 CONTINUE

      do 7 k=kz1,klam
        A1(K)=CZs(M,K-1)
        C1(K)=CZs(M,K)
        B1(K)=-(A1(K)+C1(K))

         Axx=CX(M-1,K)
         Cxx=CX(M,K)
         Bxx=-(Axx+Cxx)
         F1(K)=HGM(M,K)*SS(M,K)/dt
     *         -(A1(K)*SS(M,K-1)+B1(K)*SS(M,K)+C1(K)*SS(M,K+1))
     *         -2*(Axx*SS(m-1,K)+Bxx*SS(m,K)+Cxx*SS(m+1,K))
         B1(K)=B1(K)+HGM(M,K)/dt
    7 continue
c      kz0=3
c*    Rivers!!!!!!!!!!!!!!
c      if(mask(m,n).eq.1) then
c         P2=slev(kp(m,n))
c         p3=0.
c         I1=3
c         I2=Klam
c         Call Fact1(A1,B1,C1,F1,U1,p2,p3,i1,i2)  
c         do  K=KZ1,klam
c            SS1(M,K)=U1(K)
c         end do
c      else
         f1(kz0)=f1(kz0)-sflux(m,n)*hgm0(m)*(3-kz0)
     *	   +gsalt*Sstar(m,n)*hgm0(m)
	B1(kz0)=B1(kz0)+gsalt*hgm0(m)
         P2=0
         p3=0.
         I1=KZ0
         I2=Klam
         Call Fact1(A1,B1,C1,F1,U1,p2,p3,i1,i2)  
         do  K=KZ1,klam
            SS1(M,K)=U1(K)
         end do
    4 CONTINUE

C************************************

      DO 8 K=KZ0,KKZ
      DO 9 M=M1,M2
        A1(M)=CX(M-1,K)
        C1(M)=CX(M,K)
        B1(M)=-(A1(M)+C1(M))
        FS1(M)=HGM(M,K)*SS1(M,K)/dt
     *        +(A1(M)*SS(M-1,K)+B1(M)*SS(M,K)+C1(M)*SS(M+1,K))
c     *        -SS1(m,k)*hgm(m,k)/dlx(m,n)*rivx(m,n)
c     *          *inbt(m,n,k)
         B1(M)=B1(M)+HGM(M,K)/dt
  9   continue
         P2=0.
         ps2=0.
         p3=0.
         I1=M1
         I2=M2
         Call Fact1(A1,B1,C1,Fs1,Us1,p2,ps2,i1,i2)

      DO  M=M1,M2
            SS2(M,K)=US1(M)
      end do      
    8 CONTINUE

      goto 2
  201 CONTINUE
                 DO  M=1,mh
                     K2=KP(M,N)
                 DO  K=kz0,KKZ
                   SAL(K2,K-1)=SS2(M,K)
                 end do
	            end do
   73 CONTINUE
C*************************************************
      RETURN
      END

