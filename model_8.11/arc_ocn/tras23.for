        SUBROUTINE tras23(nt,Sal,uu,vv,ww,rivx,rivy,tbin,sbin,Ftok)
C**********************************************
C*                                            *
C*    CROSXZ:                                 *
C*                RU(T)=10*RU(S)              *
C**********************************************
      Include 'model.par'
      real *8 Sal(nodes,kgr),Ftok(mh,nh)
     *       ,uu(nodes,kgr),vv(nodes,kgr),ww(nodes,kgr)
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
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      COMMON/d1/ Rad,Gro,R,dt
c      Real *8 SL,SM,SMB,SLB,DX,DY,Z,DZ,DZK
      common/adv_fluxes/ tyflux(mh,nh,kh),syflux(mh,nh,kh)
	REAL *8 syflux,tyflux
      real *8 A1(1000),B1(1000),C1(1000),F1(1000),U1(1000),
     *        FS1(1000),US1(1000)

      REAL *8 SS(mh,kh),umx1(mh,kh)
     *       ,SS1(mh,kh) ,SS2(mh,kh)
     *        ,gk1,gk2
     *        ,ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
      REAL *8 RZ1(mh,kh),RX(mh,kh)
     *	,tzflux(mh,kh),txflux(mh,kh),szflux(mh,kh),sxflux(mh,kh)
      REAL *8 CVZ(mh,kh),CVX(mh,kh),CZ(mh,kh),CX(mh,kh),
     *	rivx(mh,nh),rivy(mh,nh),tbin(3),sbin(3)
      REAL *8 HGM(mh,kh),czs(mh,kh),hgm0(mh),ghan,cnr
     *    ,gradl,gradls,curvl,curvls,sfll,sfl,sflls,sfls
     *        ,curz,curm,wer,ut,ub
     *    ,gradup,gradups,curvz,curvzs,sflup,sfup,sflups,sfups
      integer tracx(mh,nh),tracxy(mh,nh)
	Common/Xtrac/tracx,tracxy 


      REAL *8 DZ0,DZ1,Axx,Bxx,Cxx,tadd,sadd,dzgam
      COMMON/NU/RU1(nh)
c	 if(nday/5*5.eq.nday) cnr=1
      KZ1=2
      NZV=0
      LL=0
      NZL=0
C*------------------------------------------------
      DO 73 N=N0,nh-1
c       print *,n
c       read *
       DO 222 K=1,kh
      DO 222 M=1,mh
           HGM(M,K)=0.
	sxflux(m,k)=0
	szflux(m,k)=0
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
C*************************************************
      LK=1
      NCROSS=0
    2 LL=LL+1
      X3=0.
      IF(ll.le.kxLine.and.IXD(1,LL).EQ.N)GOTO 3
c      print *,ll,n,ixd(1,ll)
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
c      end if
         lk=l
      DO 20 K=2,KKZ
C*    IF(X3.EQ.0)U(M1,K)=U(M1,K)*0.5
C*    IF(X3.EQ.0)U(M2-1,K)=U(M2-1,K)*0.5
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
c               top1=ikbt(m,n,k)
c               top2=ikbt(m-1,n,k)
c               top3=ikbt(m-1,n-1,k) 
c               top4=ikbt(m,n-1,k)
c                     bot1=ikbt(m,n,k+1)
c                     bot2=ikbt(m-1,n,k+1)
c                     bot3=ikbt(m-1,n-1,k+1) 
c                     bot4=ikbt(m,n-1,k+1)

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
CCCCCCCCCC+++++++++++++++++++++++++++++++++++++++++
       Gam=((ds1+ds2+ds3+ds4)*dz0
     *      +(ds1+ds2+ds3+ds4)*dz1)/8. 
       dzgam=(dz0+dz1)/2.
          if(k.eq.klam)dzgam=dz0/2.
          if(gam.eq.0) gam=ds1*dzgam
            zs1=1
          if(k.eq.2)then
               zs1=0
             hgm0(m)=(ds1*top1+ds2*top2+ds3*top3+ds4*top4)/4.*3.
          end if
      HGM(M,K)=6.*GAM*2.0

    

      gradls=(ss(m,k)-ss(m-1,k))/dlx(m-1,n)*inbt(m-1,n,k)*inbt(m,n,k)
      curm= (uu(kp(m-1,n-1),k-1)*(dly(m-1,n-1)+dly(m,n-1))
     *+uu(kp(m-1,n),k-1)*(dly(m-1,n)+dly(m,n)))
     */(dly(m,n)+dly(m,n-1)+dly(m-1,n)+dly(m-1,n-1))
     */dlx(m-1,n)*inbt(m-1,n,k)*inbt(m,n,k)*dt
      if(n.eq.2) curm=0.D+0


c--begin---
      if(curm.gt.0) then
	 icurm=inbt(m-1,n,k)*inbt(m-2,n,k)
         if(icurm.eq.0) then
            curvl=0
            curvls=0
         else
          curvls=2*((ss(m,k)-ss(m-1,k))/dlx(m-1,n)
     *          *inbt(m-1,n,k)*inbt(m,n,k)
     *             -(ss(m-1,k)-ss(m-2,k))/dlx(m-2,n)
     *          *inbt(m-1,n,k)*inbt(m-2,n,k))
     *              /(dlx(m-1,n)
     *             +dlx(m-2,n))
	    end if
      end if
c--end------

c===begin===
      if(curm.lt.0) then
	 icurm=inbt(m,n,k)*inbt(m+1,n,k)
         if(icurm.eq.0) then
	                     curvl=0
		curvls=0

          	else
          curvls=2*((ss(m+1,k)-ss(m,k))/dlx(m,n)
     *          *inbt(m+1,n,k)*inbt(m,n,k)
     *             -(ss(m,k)-ss(m-1,k))/dlx(m-1,n)
     *          *inbt(m-1,n,k)*inbt(m,n,k))
     *              /(dlx(m-1,n)+dlx(m,n))
	    end if
      end if
c===end=====
      sflls=0.5*(ss(m,k)+ss(m-1,k))
	sfls=0.5*(ss(m,k)+ss(m-1,k))
     *    -curm*gradls*dlx(m-1,n)/2
     *     -(1-curm**2)*curvls*(dlx(m-1,n)**2)/6.


c++++++++zflux
      if(k.gt.2) then
      gradups=(ss(m,k)-ss(m,k-1))/dz(k-1)
      curz=ww(kp(m,n),k-2)/dz(k-1)*dt
	else
      gradup=0
      gradups=0
      curz=0
         end if

c--begin---
      if(curz.gt.0) then
         if(k.le.3) then
	                     curvz=0
	                     curvzs=0
          	else
          curvzs=2*((ss(m,k)-ss(m,k-1))/dz(k-1)
     *             -(ss(m,k-1)-ss(m,k-2))/dz(k-2))
     *              /(dz(k-1)+dz(k-2))
	    end if
      end if
c--end------

c===begin===
      if(curz.lt.0) then
         if(k.eq.2.) then
	                     curvz=0
	                     curvzs=0
          	else
          curvzs=2*((ss(m,k+1)-ss(m,k))/dz(k)*inbt(m,n,k+1)
     *             -(ss(m,k)-ss(m,k-1))/dz(k-1)*inbt(m,n,k))
     *              /(dz(k)+dz(k-1))
	    end if
      end if
c===end=====
      sflups=0.5*(ss(m,k-1)+ss(m,k))*zs1
	sfups=sflups 
     *      -curz*gradups*dz(k-1)/2
     *     -(1-curz**2)*curvzs*(dz(k-1)**2)/6.


      sxflux(m-1,k)=sfls*curm
     *      *6*dzk(k)*dlx(m-1,n)*
     * 0.25*(dly(m,n)+dly(m,n-1)+dly(m-1,n)+dly(m-1,n-1))

      szflux(m,k-1)=6*curz*sfups*dz0*
     *           (ds1*top1+ds2*top2+ds3*top3+ds4*top4)/4.
c      if(n.eq.nh-1.and.m.le.70.and.m.ge.68)szflux(m,k)=szflux(m,k-1)
c      if(n.eq.nh-1.and.m.le.70.and.m.ge.68)tzflux(m,k)=tzflux(m,k-1)
   41 CONTINUE
   40 CONTINUE

c      goto 700
      m=m1
	if(rivx(m,n).ne.0) then
c      print *, 'rivx',m,m,k
	klam=KBT(M,N)
	                     strac=1.
                           if(nt.eq.tracx(m,n))strac=100.
	                    if(nt.eq.tracxy(m,n))strac=100.
c	if(nt.eq.tracx(m,n))strac=1.
         do K=KZ1,klam

	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m-1,n))
     *           *inbt(m,n,k)*inbt(m,n+1,k)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m-1,n-1))
     *           *inbt(m,n,k)*inbt(m,n-1,k)
	

         curm=(ut*(dly(m-1,n)+dly(m,n))
     *        +ub*(dly(m-1,n-1)+dly(m,n-1)))*dt
c     *        /(dly(m,n)+dly(m-1,n)+dly(m,n-1)+dly(m-1,n-1))
c     *        /dlx(m-1,n)


	sfls=0.5*(ss(m,k)+strac)
      sxflux(m-1,k)=sfls*curm
     *             *6*dzk(k)*0.25
c     * *dlx(m-1,n)*(dly(m,n)+dly(m-1,n)+dly(m,n-1)+dly(m-1,n-1))

      end do
	end if

      m=m2
	if(rivx(m,n).ne.0) then
	klam=KBT(M,N)
	                     strac=1.
c                           if(nt.eq.tracx(m,n))strac=1
                	if(nt.eq.tracx(m,n))strac=100
	            if(nt.eq.tracxy(m,n))strac=100.
         do K=KZ1,klam


	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m+1,n))
     *           *inbt(m,n,k)*inbt(m,n+1,k)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m+1,n-1))
     *           *inbt(m,n,k)*inbt(m,n-1,k)


      curm= (ut*(dly(m,n)+dly(m+1,n))
     *       +ub*(dly(m,n-1)+dly(m+1,n-1)))
c     *       /(dly(m,n)+dly(m,n-1)+dly(m+1,n)+dly(m+1,n-1))
c     *       /dlx(m,n)
     *       *dt
       
	sfls=0.5*(ss(m,k)+strac)
      sxflux(m,k)=sfls*curm
     *      *6*dzk(k)*0.25
c     *      *dlx(m,n)
c     *      *(dly(m,n)+dly(m,n-1)+dly(m+1,n)+dly(m+1,n-1))

      end do
	end if


  700 continue



   16 CONTINUE

      DO 8 K=KZ0,KKZ
      DO 9 M=M1,M2

        FS1(M)=HGM(M,K)*SS(M,K)
     *    -(sxflux(m,k)-sxflux(m-1,k))*2
     *    -(szflux(m,k)-szflux(m,k-1))*2
     *    -(syflux(m,n,k)-syflux(m,n-1,k))*2
c     *        -SS(m,k)*hgm(m,k)/dly(m,n)*rivy(m,n)*dt
c     *        -SS(m,k)*hgm(m,k)/dlx(m,n)*rivx(m,n)*dt

         B1(M)=HGM(M,K)

c       if(m.le.70.and.m.ge.68.and.n.eq.nh-1) then

c	F1(m)=F1(m)

c     *	+(tbin(m-67)-TS(m,k))*HGM(m,k)*rivy(m,n)/dly(m,n)

c     *         *inbt(m,n,k)*dt

c	  FS1(m)=FS1(m)+sbin(m-67)*HGM(m,k)*rivy(m,n)/dly(m,n)*dt

c     *         *inbt(m,n,k)

c	end if
  9   continue
         P2=0.
         ps2=0.
         p3=0.
         I1=M1
         I2=M2
c       if(n.ge.nc1.and.n.le.nc2.and.a1(m1).ne.0) then
c      Call FctCTS(A1,B1,C1,F1,Fs1,U1,Us1,i1,i2)
c      else  
c        Call FactTS(A1,B1,C1,F1,Fs1,U1,Us1,p2,ps2,p3,i1,i2)
c      end if
      DO  M=M1,M2
c            SS2(M,K)=US1(M)
c            TS2(M,K)=U1(M)
            SS2(M,K)=FS1(m)/b1(m)
      end do      
    8 CONTINUE

      goto 2
  201 CONTINUE
                 DO  M=2,mh
                     K2=KP(M,N)
                 DO  K=kz0,KKZ
                   SAL(K2,K-1)=SS2(M,K)
                 end do
	            end do
   73 CONTINUE
C*************************************************
      RETURN
      END
