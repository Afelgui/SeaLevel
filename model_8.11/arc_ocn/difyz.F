          SUBROUTINE difYZ(Temp,Sal,rivy,tbin,sbin)
c	implicit none
C*************************************************
C*                                               *
C*      CROSYZ:   T,S(N=47)= T,S(KLIM)           *
C*                RU(T)=10*RU(S)                 *
C*************************************************
      Include 'model.par'
      common/flux/tflux(mh,nh),sflux(mh,nh)
      real *8 tflux,sflux
      common/levitus/Tlev(nodes),Slev(nodes)
      real *8 tlev,slev
      common/restor/Sstar(mh,nh),gsalt
      real *8 Sstar,gsalt
      COMMON/DAY/myear,month,nday
      COMMON/NH2/ NH2
      real *8 z,dz,dzk,dx,dy,h,h1,sl,slb
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
      common/Cfsmag/cfsmag(mh,nh,kgr)
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      real *8 Temp(nodes,kgr),Sal(nodes,kgr)
      COMMON/RUNPARS/nmyear
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/north/tatl(16,4,23),satl(16,4,23)
c      REAL *8 SL,SM,SMB,SLB,DX,DY,Z,DZ,DZK
      real *8 A1(1000),B1(1000),C1(1000),F1(1000),U1(1000),
     *        FS1(1000),US1(1000)
      REAL *8 AS1(kh),BS1(kh),CS1(kh),PXS1(kh),CZS(nh,kh)
      REAL *8 TS(nh,kh),SS(nh,kh),hgm0(nh)
     *       ,ts1(nh,kh),SS1(nh,kh),ts2(nh,kh),SS2(nh,kh)
     *        ,rivy(mh,nh),tbin(3),sbin(3)
      REAL *8 HGM(nh,kh),CZ(nh,kh),CY(nh,kh)
     *,        tsec(nh,kh),ssec(nh,kh),ga(nh,kh)
      REAL *8 RZ1(nh,kh),RY(nh,kh),SU(nh,kh),ru_new(nh,kh)
      REAL *8 DZ0,DZ1,ZZ1,ZZ2,Axx,Bxx,Cxx,dzgam,gk1,gk2,ghan,gk3,cnr
     *        ,ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
	real *8 p2,p3,ps2
      integer nu_diag
      common/iodiag/nu_diag
     

C*********************************************
       astar=0.8
c       astar=0.1
	 arr=1.05/3.14159265
      alam=4.5e-5
      zstar=250000

      DO k=2,kkz
	do n=1,108
c         RU(N)=1.
      ru_new(n,k)=astar+arr*atan(alam*(z(k)-zstar))
	end do
	do n=1,nh
	ru_new(n,k)=0.05
	end do
      end do 
	
      
           ga_s=2.31e-7/5.
           ga_b=4.63e-8

      if(myear.le.nmyear)then
         tt=(myear-nmyear+1)*360+(month-1)*30+min(nday,30)-1
         tt=0.5*(1.0+cos(tt*3.14159265/720.0))
      else
         tt=0.0
      endif
            abh=0.5e+6*(1.-tt)+1.0e+6*tt
            ash=0.5e+6*(1.-tt)+1.0e+6*tt
      
      NZL=0
      LL=0
      KZ1=2
      
	DO 223 K=1,kh
      DO 223 N=1,nh
             RZ1(N,K)=0.
             RY(N,K)=0.
     
  223 continue  
C**********************************************
      DO 200 M=M0,MM
      do 1 n=1,nh
      do  1 k=1,kh
            umy1=abh+(ash-abh)*exp(-z(k)/50000.)
c            umy1=1.e+9
   1  continue
      NZL=NZL+1
          DO 110 N=1,Nh-1
      DX0=DX(M-1,n)
      DX1=DX(M,n)
             K2=KP(M,N)
           do k=1,kh
              CZS(N,K)=0.
             CY(N,K)=0.
            cy(n,k)=0.
            cz(n,k)=0.
          end do
           DO 110 K=2,KKZ
              tsec(n,k)=0
              ssec(n,k)=0
       ga(n,k)=0.
              TS(N,K)=TEMP(K2,K-1)
              SS(N,K)=SAL(K2,K-1)
                  ts1(n,k)=ts(n,k)
                    ts2(n,k)=ts(n,k)
                       ss1(n,k)=ss(n,k)
                         ss2(n,k)=ss(n,k)
      

               HGM(n,K)=2./DT*DlX(M,n)*DlY(m,n)*DZ(K)

  110        CONTINUE
C*********************************************
      lk=1
      x3=0.
      l=1
           jstep=1
           if(x3.eq.1)jstep=kkz-2
C*********************************************
      NCROSS=0
    2 LL=LL+1
      IF(ll.le.kyLine.and.IYD(1,LL).EQ.M)GOTO 3
      LL=LL-1
      GOTO 201
    3 N1=IYD(2,LL)
      N2=IYD(3,LL)
C      Write(nu_diag,*)M,N1,N2
C      READ *
      NCROSS=NCROSS+1
      LN=LK+1
      NRA=0
      IF(N1.EQ.N2)GOTO2
         lk=l
      DO 12 K=2,KKZ
      DO 12 N=N1,N2
            K2=KP(M,N)
            RZ1(N,K)=ru_new(n,k)/2.
            umy1=abh+(ash-abh)*exp(-z(k)/50000.)
            ry(n,k)=Umy1+cfsmag(m,n,k-1)*0

   12 CONTINUE
      DO  K=2,KKZ
       ry(2,k)=5.e+9
	ry(3,k)=5.e+9
	ry(4,k)=1.e+9
	ry(5,k)=5.e+8
	ry(6,k)=1.e+8
	ry(7,k)=5.e+7
	end do
      DO 40 N=N1,N2
       ghan=(dz(2)+dz(3)*0.5)/(30.*86400.0)
       KTET=KBT(M,N)
	k2=kp(m,n)
           ds1=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.
           ds2=(dlx(m-1,n)+dlx(m-1,n+1))*(dly(m,n)+dly(m-1,n))/4.
           ds3=(dlx(m-1,n)+dlx(m-1,n-1))*(dly(m,n-1)+dly(m-1,n-1))/4.
           ds4=(dlx(m,n)+dlx(m,n-1))*(dly(m,n-1)+dly(m+1,n-1))/4.
      DO 41 K=KZ1,KTET
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
      
           zs1=1
             dzgam=(dz0+dz1)/2.
          if(k.eq.ktet)dzgam=dz0/2.
          if(gam.eq.0) gam=ds1*dzgam
          if(k.eq.2)then
               zs1=0
             hgm0(n)=(ds1*top1+ds2*top2+ds3*top3+ds4*top4)/4.*3.*2.
             
          end if
                      
      HGM(n,K)=6.D+0*GAM*2.D+0
      Cy(n,K)=-DZ0*(Ry(n,K-1)+Ry(n,K))*1.5D+0
     *        *0.5D+0*((dlx(m,n)+dlx(m,n+1))*0.5D+0/dly(m,n)*top1
     *             +(dlx(m-1,n)+dlx(m-1,n+1))*0.5D+0/dly(m,n)*top2)
     *        -DZ1*(Ry(n,K+1)+Ry(n,K))*1.5D+0
     *        *0.5D+0*((dlx(m,n)+dlx(m,n+1))*0.5D+0/dlx(m,n)*bot1
     *             +(dlx(m-1,n)+dlx(m-1,n+1))*0.5D+0/dly(m,n)*bot2)

      cy(n,k)=cy(n,k)*inbt(m,n+1,k)

      CZ(n,K)=-((ds1*bot1+ds2*bot2)*0.5*(RZ1(n,K)+RZ1(n,K+1))*1.5
     *         +(ds3*bot3+ds4*bot4)*0.5*(RZ1(n,K)+RZ1(n,K+1))*1.5)
     *         /DZ1

      CZs(n,K)=-((ds1*bot1+ds2*bot2)*0.5*(RZ1(n,K)+RZ1(n,K+1))*1.5
     *          +(ds3*bot3+ds4*bot4)*0.5*(RZ1(n,K)+RZ1(n,K+1))*1.5)
     *         /DZ1

   41 CONTINUE

       cz(n,ktet)=0.
       czs(n,ktet)=0.


   40 CONTINUE

      do k=1,kkz
      cy(n2,k)=0.
      cy(n1-1,k)=0.
      end do

C*****************************
   14 CONTINUE
      NRA=NRA+1


c	goto 66

      DO 4 N=N1,N2
C------------------------------------------
            KTET=KBT(M,N)
            DO 6 K=KZ1,KTET
  
	if (m.eq.69.and.n.eq.181.and.k.eq.2) then
	continue
	end if

          A1(K)=CZ(N,K-1)
          C1(K)=CZ(N,K)
          B1(K)=-(A1(K)+C1(K))
      Axx=Cy(n-1,K)
      Cxx=Cy(n,K)
      Bxx=-(Axx+Cxx)
          F1(K)=HGM(N,K)/dt*TS(N,K)-
     *         (A1(K)*TS(N,K-1)+B1(K)*TS(N,K)+C1(K)*TS(N,K+1))
     *     -2*(Axx*TS(n-1,K)+Bxx*TS(n,K)+Cxx*TS(n+1,K))
 

    6      B1(K)=B1(K)+HGM(N,K)*(1./dt+ga(n,k))
         f1(kz0)=f1(kz0)-tflux(m,n)*hgm0(n)*(3-kz0)
              I1=KZ0
              I2=Ktet
              P2=tlev(kp(m,n))*(KZ0-2)
              P3=0.
           Call Fact1(A1,B1,C1,F1,U1,p2,p3,i1,i2)  
                     DO  K=KZ1,Ktet
                     Ts1(N,K)=U1(K)
                     end do
      DO 7 K=KZ1,Ktet
         A1(K)=CZs(N,K-1)
         C1(K)=CZs(N,K)
         B1(K)=-(A1(K)+C1(K))
      Axx=Cy(n-1,K)
      Cxx=Cy(n,K)
      Bxx=-(Axx+Cxx)
         F1(K)=HGM(N,K)/dt*SS(N,K)-
     *          (A1(K)*SS(N,K-1)+B1(K)*SS(N,K)+C1(K)*SS(N,K+1))
     *     -2*(Axx*SS(n-1,K)+Bxx*SS(n,K)+Cxx*SS(n+1,K))
     *        +ga(n,k)*hgm(n,k)*ssec(n,k)
         B1(K)=B1(K)+HGM(N,K)*(1./dt+ga(n,k))
    7 continue


         f1(kz0)=f1(kz0)-sflux(m,n)*hgm0(n)*(3-kz0)
     *	   +gsalt*Sstar(m,n)*hgm0(n)
	B1(kz0)=B1(kz0)+gsalt*hgm0(n)
         I1=KZ0
         I2=ktet
         P2=slev(kp(m,n))*(KZ0-2)
         P3=0.
         Call Fact1(A1,B1,C1,F1,U1,p2,p3,i1,i2)  
         DO  K=KZ1,ktet
            SS1(N,K)=U1(K)
         end do
    4 CONTINUE
C*****************************
c      IF(NRA.EQ.4)GOTO2
C*****************************
   15 CONTINUE
      NRA=NRA+1

66	continue
c	goto 67

      DO 8 K=KZ0,KKZ
      DO 9 N=N1,N2

		if (m.eq.60.and.n.eq.4.and.k.eq.2) then
		continue
	end if

		if (m.eq.68.and.n.eq.181) then
c		print *,k
	end if

      A1(N)=6.D+0*CY(N-1,K)
      C1(N)=6.D+0*CY(N,K)
      B1(N)=-6.D+0*CY(N-1,K)-6.D+0*CY(N,K)
  301 F1(N)=6.D+0*HGM(N,K)/dt*TS1(N,K)
     *      +(A1(N)*TS(N-1,K)+B1(N)*TS(N,K)+C1(N)*TS(N+1,K))
c     F1(N)=HGM(N,K)/dt

      FS1(N)=6.D+0*HGM(N,K)/dt*SS1(N,K)+
     *      (A1(N)*SS(N-1,K)+B1(N)*SS(N,K)+C1(N)*SS(N+1,K))

c     FS1(N)=HGM(N,K)/dt*SS1(N,K)
c     *          *inbt(m,n,k)

    9 B1(N)=B1(N)+6.D+0*HGM(N,K)/dt
c***   Bering Strait: n=181 , m=68,69,70

c      if(n2.eq.181) then

c         F1(N2)=F1(N2)
c     *	    +(tbin(m-67)-TS1(n2,k))*HGM(n2,k)*rivy(m,n2)/dly(m,n2)
c     *      *inbt(m,n2,k)
c
c         FS1(N2)=FS1(N2)+sbin(m-67)*HGM(n2,k)*rivy(m,n2)/dly(m,n2)
c     *      *inbt(m,n2,k)

c         FS1(N2)=FS1(N2)+SS1(n2,k)*hgm(n2,k)/dly(m,n2)*rivy(m,n2)

c      end if

            p2=0.
            ps2=0.
            p3=0.
             I1=N1
             I2=N2
      call FACTTS(A1,B1,C1,F1,Fs1,U1,Us1,p2,ps2,p3,i1,i2)
c	if (m.eq.60.and.k.eq.2) then
c	read*
c	end if
      DO 10 N=N1,N2
            SS2(N,K)=US1(N)
   10       TS2(N,K)=U1(N)
    8 CONTINUE

   67  continue

      goto 2
  201 CONTINUE
      DO  N=1,nh-1
      K2=KP(M,N)
      DO  K=kz0,KKZ
                  TEMP(K2,K-1)=TS2(N,K)
                  SAL(K2,K-1)=SS2(N,K)
	end do
	end do
  200 CONTINUE

C********************************************
      RETURN
      END

 
