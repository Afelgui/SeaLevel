        SUBROUTINE qstx(Temp,Sal,uu,vv,ww,ww1,rivx,rivy,tbin,sbin
     *	                                  ,Ftok,lvl,lvl0)
C*                                            *
C*    CROSXZ:                                 *
C*                RU(T)=10*RU(S)              *
C**********************************************
      Include 'model.par'
      real *8 Temp(nodes,kgr),Sal(nodes,kgr),Ftok(mh,nh)
     *       ,uu(nodes,kgr),vv(nodes,kgr),ww(nodes,kgr),ww1(nodes,kgr)
      common/flux/tflux(mh,nh),sflux(mh,nh)
      real *8 tflux,sflux
      common/levitus/Tlev(nodes),Slev(nodes)
      real *8 tlev,slev
      COMMON/DAY/myear,month,nday
      COMMON/NH2/ NH2
      real *8 z,dz,dzk,dx,dy,h,h1
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
	
	common/cwz/   cwz00(mh,nh),	cwzvel00(mh,nh)

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
	common/iodiag/nu_diag
c      Real *8 SL,SM,SMB,SLB,DX,DY,Z,DZ,DZK
         common/adv_fluxes/ tyflux(mh,nh,kh),syflux(mh,nh,kh),
     *						txflux(mh,nh,kh),sxflux(mh,nh,kh)
     *						,tzflux(mh,nh,kh),szflux(mh,nh,kh)

	REAL *8 syflux,tyflux
      real *8 A1(1000),B1(1000),C1(1000),F1(1000),U1(1000),
     *        FS1(1000),US1(1000)
       REAL *8  k_lt(mh,kh)
      REAL *8 TS(mh,kh),SS(mh,kh),umx1(mh,kh)
     *       ,ts1(mh,kh),SS1(mh,kh) ,ts2(mh,kh),SS2(mh,kh)
     *        ,gk1,gk2
     *        ,ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
      REAL *8 RZ1(mh,kh),RX(mh,kh),lvl(mh,nh),lvl0(mh,nh)
     *	,tzflux,szflux,txflux,sxflux
c     *	,tzflux(mh,kh),szflux(mh,kh),txflux,sxflux
c     * ,txflux(mh,kh),sxflux(mh,kh)
      REAL *8 CVZ(mh,kh),CVX(mh,kh),CZ(mh,kh),CX(mh,kh),
     *	rivx(mh,nh),rivy(mh,nh),tbin(3),sbin(3)
      REAL *8 HGM(mh,kh),czs(mh,kh),hgm0(mh),ghan,cnr
     *    ,gradl,gradls,curvl,curvls,sfll,sfl,sflls,sfls
     *        ,curz,curm,wer,ut,ub,wer1
     *    ,gradup,gradups,curvz,curvzs,sflup,sfup,sflups,sfups

      REAL *8 DZ0,DZ1,Axx,Bxx,Cxx,tadd,sadd,dzgam,cwz00,cwzvel00

	REAL *8 u1_,u2_,v1_,v2_,curm1
     *	,txflux1,tzflux1,txflux2,tzflux2,tyflux1,tyflux2
      COMMON/NU/RU1(nh)
			common/iso_slopes/ sx(mh,nh,kh),sy(mh,nh,kh)
	     real *8 sx,sy
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
	txflux(m,n,k)=0
	sxflux(m,n,k)=0
	tzflux(m,n,k)=0
	szflux(m,n,k)=0
  222 continue  
  
  
  
      do 4 m=1,mh
      do  4 k=1,kh
       k_lt(m,k)=1.e+6
   4  continue


     
      l=1
C*------------------------------------------------
      NZL=NZL+1
      DO 81 M=1,mh
      K2=KP(M,N)
      DY0=DY(m,N-1)
      DY1=DY(m,N)
      DO 81 K=2,KKZ
            TS(M,K)=TEMP(K2,K-1)
            SS(M,K)=SAL(K2,K-1)
c	 SS(M,K)=1.
                  ts1(m,k)=ts(m,k)
                    ts2(m,k)=ts(m,k)
                       ss1(m,k)=ss(m,k)
                         ss2(m,k)=ss(m,k)
               HGM(M,K)=2.D+0/DT*DlX(M,n)*DlY(m,n)*DZ(K)
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
		if (m.eq.87) then
	continue
	end if
       ghan=(dz(2)+dz(3)*0.5)/(30.*86400.0)
       KLAM=KBT(M,N)
	k2=kp(m,n)
              ds1=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.D+0
           ds2=(dlx(m-1,n)+dlx(m-1,n+1))*(dly(m,n)+dly(m-1,n))/4.D+0
           ds3=(dlx(m-1,n)+dlx(m-1,n-1))*(dly(m,n-1)+dly(m-1,n-1))/4.D+0
           ds4=(dlx(m,n)+dlx(m,n-1))*(dly(m,n-1)+dly(m+1,n-1))/4.D+0
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
     *      +(ds1*bot1+ds2*bot2+ds3*bot3+ds4*bot4)*dz1)/8.D+0 
CCCCCCCCCC+++++++++++++++++++++++++++++++++++++++++
       Gam=((ds1+ds2+ds3+ds4)*dz0
     *      +(ds1+ds2+ds3+ds4)*dz1)/8. D+0
       dzgam=(dz0+dz1)/2.D+0
          if(k.eq.klam)dzgam=dz0/2.D+0
          if(gam.eq.0) gam=ds1*dzgam
            zs1=1
          if(k.eq.2)then
               zs1=0
             hgm0(m)=(ds1*top1+ds2*top2+ds3*top3+ds4*top4)/4.D+0*3.D+0
          end if



      HGM(M,K)=6.*GAM*2.0

    

      gradl=(ts(m,k)-ts(m-1,k))/dlx(m-1,n)*inbt(m-1,n,k)*inbt(m,n,k)
      gradls=(ss(m,k)-ss(m-1,k))/dlx(m-1,n)*inbt(m-1,n,k)*inbt(m,n,k)
      curm= (uu(kp(m-1,n-1),k-1)*(dly(m-1,n-1)+dly(m,n-1))
     *+uu(kp(m-1,n),k-1)*(dly(m-1,n)+dly(m,n)))
     */(dly(m,n)+dly(m,n-1)+dly(m-1,n)+dly(m-1,n-1))
     */dlx(m-1,n)*inbt(m-1,n,k)*inbt(m,n,k)*dt
c      if(n.eq.2) curm=0.D+0


c--begin---
      if(curm.gt.0) then
	 icurm=inbt(m-1,n,k)*inbt(m-2,n,k)
         if(icurm.eq.0) then
            curvl=0.
            curvls=0.
         else
            curvl=2*((ts(m,k)-ts(m-1,k))/dlx(m-1,n)
     *          *inbt(m-1,n,k)*inbt(m,n,k)
     *             -(ts(m-1,k)-ts(m-2,k))/dlx(m-2,n)
     *          *inbt(m-1,n,k)*inbt(m-2,n,k))
     *              /(dlx(m-1,n)
     *             +dlx(m-2,n))
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
	                     curvl=0.
		curvls=0.

          	else
          curvl=2*((ts(m+1,k)-ts(m,k))/dlx(m,n)
     *          *inbt(m+1,n,k)*inbt(m,n,k)
     *             -(ts(m,k)-ts(m-1,k))/dlx(m-1,n)
     *          *inbt(m-1,n,k)*inbt(m,n,k))
     *              /(dlx(m-1,n)+dlx(m,n))
          curvls=2*((ss(m+1,k)-ss(m,k))/dlx(m,n)
     *          *inbt(m+1,n,k)*inbt(m,n,k)
     *             -(ss(m,k)-ss(m-1,k))/dlx(m-1,n)
     *          *inbt(m-1,n,k)*inbt(m,n,k))
     *              /(dlx(m-1,n)+dlx(m,n))
	    end if
      end if
c===end=====
      sfll=0.5D+0*(ts(m,k)+ts(m-1,k))
	sfl=0.5D+0*(ts(m,k)+ts(m-1,k))
     *    -curm*gradl*dlx(m-1,n)/2.D+0
     *     -(1-curm**2)*curvl*(dlx(m-1,n)**2)/6.D+0
      sflls=0.5D+0*(ss(m,k)+ss(m-1,k))
	sfls=0.5D+0*(ss(m,k)+ss(m-1,k))
     *    -curm*gradls*dlx(m-1,n)/2.D+0
     *     -(1-curm**2)*curvls*(dlx(m-1,n)**2)/6.D+0


	
		if (m.eq.61.and.n.eq.158.and.k.eq.4) then
	continue
	end if


		if (m.eq.60.and.n.eq.158.and.k.eq.4) then
	continue
	end if


    	goto 765
c     ULTIMATE
c------------temp-----------------


      if(curm.gt.0) then
         ii=max(m1+1,m-2)
	sups=inbt(ii,n,k)*ts(ii,k)+(1-inbt(ii,n,k))*ts(m,k)
	if(abs(ts(m,k)-sups).lt.1.e-5) then
	  sfl=ts(m-1,k)
	else
	   scnor=(ts(m-1,k)-sups)/(ts(m,k)-sups)
	   srnor=(sfl-sups)/(ts(m,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(scnor/curm,1.)
              clwlim=scnor
                 if(srnor.gt.cuplim) srnor=cuplim
                 if(srnor.lt.clwlim) srnor=clwlim
          else
               srnor=scnor
          end if    
      sfl=sups+srnor*(ts(m,k)-sups)
	end if
	else
c#####
         ii=min(m2,m+1)
	sups=inbt(ii,n,k)*ts(ii,k)+(1-inbt(ii,n,k))*ts(m-1,k)
	if (abs(ts(m-1,k)-sups).lt.1.e-5) then



	  sfl=ts(m,k)	  
	  	else
	   scnor=(ts(m,k)-sups)/(ts(m-1,k)-sups)
	   slnor=(sfl-sups)/(ts(m-1,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(-scnor/curm,1.)
              clwlim=scnor
                 if(slnor.gt.cuplim) slnor=cuplim
                 if(slnor.lt.clwlim) slnor=clwlim
          else
               slnor=scnor
          end if    
      sfl=sups+slnor*(ts(m-1,k)-sups)
	end if
	end if
c--------------sal---------------------
	if(curm.gt.0) then
         ii=max(m1+1,m-2)
	sups=inbt(ii,n,k)*ss(ii,k)+(1-inbt(ii,n,k))*ss(m,k)
	if (abs(ss(m,k)-sups).lt.1.e-5) then
	  sfls=ss(m-1,k)
	else
	   scnor=(ss(m-1,k)-sups)/(ss(m,k)-sups)
	   srnor=(sfls-sups)/(ss(m,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(scnor/curm,1.)
              clwlim=scnor
                 if(srnor.gt.cuplim) srnor=cuplim
                 if(srnor.lt.clwlim) srnor=clwlim
          else
               srnor=scnor
          end if    
      sfls=sups+srnor*(ss(m,k)-sups)
	end if
	else
c#####
         ii=min(m2-1,m+1)
	sups=inbt(ii,n,k)*ss(ii,k)+(1-inbt(ii,n,k))*ss(m-1,k)
	if (abs(ss(m-1,k)-sups).lt.1.e-5) then
	  sfls=ss(m,k)
	else
	   scnor=(ss(m,k)-sups)/(ss(m-1,k)-sups)
	   slnor=(sfls-sups)/(ss(m-1,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(-scnor/curm,1.)
              clwlim=scnor
                 if(slnor.gt.cuplim) slnor=cuplim
                 if(slnor.lt.clwlim) slnor=clwlim
          else
               slnor=scnor
          end if    
      sfls=sups+slnor*(ss(m-1,k)-sups)
	end if
	end if



 765  continue



c++++++++zflux
      if(k.gt.2) then

      gradup=(ts(m,k)-ts(m,k-1))/dz(k-1)
      gradups=(ss(m,k)-ss(m,k-1))/dz(k-1)
      curz=(ww(kp(m,n),k-2))/dz(k-1)*dt

	else
      gradup=0
      gradups=0
      curz=0
         end if

c--begin---
      if(curz.gt.0) then
         if(k.le.3) then
	                     curvz=0.
	                     curvzs=0.
          	else
          curvz=2*((ts(m,k)-ts(m,k-1))/dz(k-1)
     *             -(ts(m,k-1)-ts(m,k-2))/dz(k-2))
     *              /(dz(k-1)+dz(k-2))
          curvzs=2*((ss(m,k)-ss(m,k-1))/dz(k-1)
     *             -(ss(m,k-1)-ss(m,k-2))/dz(k-2))
     *              /(dz(k-1)+dz(k-2))
	    end if
      end if
c--end------

c===begin===
      if(curz.lt.0) then
         if(k.eq.2.) then
	                     curvz=0.
	                     curvzs=0.
          	else
          curvz=2.D+0*((ts(m,k+1)-ts(m,k))/dz(k)*inbt(m,n,k+1)
     *             -(ts(m,k)-ts(m,k-1))/dz(k-1)*inbt(m,n,k))
     *              /(dz(k)+dz(k-1))
          curvzs=2.D+0*((ss(m,k+1)-ss(m,k))/dz(k)*inbt(m,n,k+1)
     *             -(ss(m,k)-ss(m,k-1))/dz(k-1)*inbt(m,n,k))
     *              /(dz(k)+dz(k-1))
	    end if
      end if
c===end=====
	sflup=(ts(m,k)*dzk(k)+ts(m,k-1) *dzk(k-1))
     *		  /(dzk(k)+dzk(k-1))

c      sflup=0.5D+0*(ts(m,k-1)+ts(m,k))*zs1
	sfup=sflup 
     *      -curz*gradup*dz(k-1)/2.D+0
     *     -(1-curz**2)*curvz*(dz(k-1)**2)/6.D+0

	sflups=(ss(m,k)*dzk(k)+ss(m,k-1) *dzk(k-1))
     *		  /(dzk(k)+dzk(k-1))

c      sflups=0.5D+0*(ss(m,k-1)+ss(m,k))*zs1
	sfups=sflups 
     *      -curz*gradups*dz(k-1)/2.D+0
     *     -(1-curz**2)*curvzs*(dz(k-1)**2)/6.D+0

	
		if (m.eq.60.and.n.eq.158.and.k.eq.3) then
	continue
	end if

		if (m.eq.60.and.n.eq.158.and.k.eq.4) then
	continue
	end if

      goto 746
c     ULTIMATE
c=====================temp==================
	

      if(curz.gt.0) then
         ii=max(2,k-2)
	sups=inbt(m,n,ii)*ts(m,ii)+(1-inbt(m,n,ii))*ts(m,k)
	if (abs(ts(m,k)-sups).lt.1.e-5) then
	  sfup=ts(m,k-1)
	else
	   scnor=(ts(m,k-1)-sups)/(ts(m,k)-sups)
	   srnor=(sfup-sups)/(ts(m,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(scnor/curz,1.)
              clwlim=scnor
                 if(srnor.gt.cuplim) srnor=cuplim
                 if(srnor.lt.clwlim) srnor=clwlim
          else
               srnor=scnor
          end if    
      sfup=sups+srnor*(ts(m,k)-sups)
	end if
	
c#####
	else
c#####
         ii=min(klam-1,k+1)
	sups=inbt(m,n,ii)*ts(m,ii)+(1-inbt(m,n,ii))*ts(m,k-1)
	if (abs(ts(m,k-1)-sups).lt.1.e-5) then
	  sfup=ts(m,k)
	else
	   scnor=(ts(m,k)-sups)/(ts(m,k-1)-sups)
	   slnor=(sfup-sups)/(ts(m,k-1)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(-scnor/curz,1.)
              clwlim=scnor
                 if(slnor.gt.cuplim) slnor=cuplim
                 if(slnor.lt.clwlim) slnor=clwlim
          else
               slnor=scnor
          end if    
      sfup=sups+slnor*(ts(m,k-1)-sups)
	end if
	end if
c===============sal=========

	if(curz.gt.0) then
         ii=max(2,k-2)
	sups=inbt(m,n,ii)*ss(m,ii)+(1-inbt(m,n,ii))*ss(m,k)
	if (abs(ss(m,k)-sups).lt.1.e-5)  then
	  sfups=ss(m,k-1)
	else
	   scnor=(ss(m,k-1)-sups)/(ss(m,k)-sups)
	   srnor=(sfups-sups)/(ss(m,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(scnor/curz,1.)
              clwlim=scnor
                 if(srnor.gt.cuplim) srnor=cuplim
                 if(srnor.lt.clwlim) srnor=clwlim
          else
               srnor=scnor
          end if    
      sfups=sups+srnor*(ss(m,k)-sups)
	end if
	
c#####
	else
c#####
         ii=min(klam-1,k+1)
	sups=inbt(m,n,ii)*ss(m,ii)+(1-inbt(m,n,ii))*ss(m,k-1)
	if(abs(ss(m,k-1)-sups).lt.1.e-5) then
	  sfups=ss(m,k)
	else
	   scnor=(ss(m,k)-sups)/(ss(m,k-1)-sups)
	   slnor=(sfup-sups)/(ss(m,k-1)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(-scnor/curz,1.)
              clwlim=scnor
                 if(slnor.gt.cuplim) slnor=cuplim
                 if(slnor.lt.clwlim) slnor=clwlim
          else
               slnor=scnor
          end if    
      sfups=sups+slnor*(ss(m,k-1)-sups)
	end if
	end if

 746  continue


      txflux(m-1,n,k)=sfl*curm
     *      *6.D+0*dzk(k)
     *    *dlx(m-1,n)*
     * 0.25D+0*(dly(m,n)+dly(m,n-1)+dly(m-1,n)+dly(m-1,n-1))



      sxflux(m-1,n,k)=sfls*curm
     *      *6.D+0*dzk(k)*dlx(m-1,n)*
     * 0.25D+0*(dly(m,n)+dly(m,n-1)+dly(m-1,n)+dly(m-1,n-1))




      tzflux(m,n,k-1)=6.D+0*curz*sfup*dz0*
     *           (ds1*top1+ds2*top2+ds3*top3+ds4*top4)/4.D+0

	tzflux(m,n,k-1)=3.D+0*curz*sfup*dz0
     *         *(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))/2.D+0

	if (m.eq.140.and.n.eq.137.and.k.eq.1) then
	continue
	end if

      szflux(m,n,k-1)=6.D+0*curz*sfups*dz0
     *         *(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))/4.D+0

	continue

c      szflux(m,n,k-1)=6.D+0*curz*sfups*dz0*
c     *           (ds1*top1+ds2*top2+ds3*top3+ds4*top4)/4.D+0


   41 CONTINUE
   40 CONTINUE

c      goto 700
      m=m1
	if(rivx(m,n).ne.0) then
c      print *, 'rivx',m,m,k
	klam=KBT(M,N)
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
	sfl=0.5D+0*(ts(m,k)+ts(m,k))
	sfls=0.5D+0*(ss(m,k)+0.D+0)
       txflux(m-1,n,k)=sfl*curm
     *               *6.D+0*dzk(k)*0.25D+0
c     * *dlx(m-1,n)*(dly(m,n)+dly(m-1,n)+dly(m,n-1)+dly(m-1,n-1))

      sxflux(m-1,n,k)=sfls*curm
     *             *6.D+0*dzk(k)*0.25D+0
c     * *dlx(m-1,n)*(dly(m,n)+dly(m-1,n)+dly(m,n-1)+dly(m-1,n-1))

      end do
	end if

      m=m2
	if(rivx(m,n).ne.0) then
	klam=KBT(M,N)
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
       
	sfl=0.5*(ts(m,k)+ts(m,k))
	sfls=0.5*(ss(m,k)+0.)
       txflux(m,n,k)=sfl*curm
     *      *6*dzk(k)*0.25
c     *      *dlx(m,n)
c     *      *(dly(m,n)+dly(m,n-1)+dly(m+1,n)+dly(m+1,n-1))
      sxflux(m,n,k)=sfls*curm
     *      *6*dzk(k)*0.25
c     *      *dlx(m,n)
c     *      *(dly(m,n)+dly(m,n-1)+dly(m+1,n)+dly(m+1,n-1))

	if (m.eq.140.and.n.eq.137.and.k.eq.2) then
	continue
	end if



      end do
	end if


  700 continue



   16 CONTINUE
      if (n.eq.nh-1) then
	jhgffds=1
	end if

      DO 8 K=KZ0,KKZ
      DO 9 M=M1,M2
      

	if (m.eq.109.and.n.eq.2.and.k.eq.4) then
	continue
	end if


	if (k.eq.2) then
       wer=
     *    -(txflux(m,n,k)-txflux(m-1,n,k))*2.D+0
     *    -(tyflux(m,n,k)-tyflux(m,n-1,k))*2.D+0
     *    -(tzflux(m,n,k))*2.D+0*dzk(k)/(dzk(k)-lvl0(m,n))
	else
       wer=
     *    -(txflux(m,n,k)-txflux(m-1,n,k))*2.D+0
     *   -(tzflux(m,n,k)-tzflux(m,n,k-1))*2.D+0 !*dzk(2)/(dzk(2)-lvl(m,n)))
     *    -(tyflux(m,n,k)-tyflux(m,n-1,k))*2.D+0
      end if
    

			if (k.eq.2) then
				wer1=
     *			3.D+0*dt*TS(m,k)*cwz00(m,n)*dzk(k)/(dzk(k)-lvl0(m,n))
			else
			 	wer1=0.
			end if

		if (m.eq.60.and.n.eq.158.and.k.eq.4) then
	continue
	end if

			if (m.eq.70.and.n.eq.70.and.k.eq.3) then
	continue
	end if

c==============TEST TEST TEST TEST===========================================
c==============TEST TEST TEST TEST===========================================
c==============TEST TEST TEST TEST===========================================
c==============TEST TEST TEST TEST===========================================
	go to 165
	if (m.eq.57.and.n.eq.3.and.k.eq.3) then
	continue

       ds1=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.D+0
       ds2=(dlx(m-1,n)+dlx(m-1,n+1))*(dly(m,n)+dly(m-1,n))/4.D+0
       ds3=(dlx(m-1,n)+dlx(m-1,n-1))*(dly(m,n-1)+dly(m-1,n-1))/4.
       ds4=(dlx(m,n)+dlx(m,n-1))*(dly(m,n-1)+dly(m+1,n-1))/4.
         

	v2_=
     * 	(vv(kp(M,n),k-1)*(dlx(m,n)+dlx(m,n+1))
     *    +vv(kp(m-1,n),k-1)*(dlx(m-1,n)+dlx(m-1,n+1)))*dzk(k)
	
	v1_=
     * 	(vv(kp(M,n-1),k-1)*(dlx(m,n-1)+dlx(m,n))
     *    +vv(kp(m-1,n-1),k-1)*(dlx(m-1,n-1)+dlx(m-1,n)))*dzk(k)
	
	
      u2_=
     *	(uu(kp(M,n),k-1)*(dly(m,n)+dly(m+1,n))
     *    +uu(kp(m,n-1),k-1)*(dly(m,n-1)+dly(m+1,n-1)))*dzk(k)
	
      u1_=
     *	(uu(kp(M-1,n),k-1)*(dly(m-1,n)+dly(m,n))
     *    +uu(kp(m-1,n-1),k-1)*(dly(m-1,n-1)+dly(m,n-1)))*dzk(k)
	

	tzflux2=3.D+0*(ww(kp(m,n),k-1))*dt
     *	*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))/2.D+0
	tzflux1=3.D+0*(ww(kp(m,n),k-2))*dt
     *	*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))/2.D+0

	txflux1=u1_*dt*3.D+0/2.D+0
	txflux2=u2_*dt*3.D+0/2.D+0
	tyflux1=v1_*dt*3.D+0/2.D+0
	tyflux2=v2_*dt*3.D+0/2.D+0

	wer=
     *    -(u2_-u1_)
     *    -(v2_-v1_)
     * -(
     *  ww(kp(m,n),k-1)*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))
     * -ww(kp(m,n),k-2)*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))
c   *dzk(k)/(dzk(k)-lvl(m,n))
     *  )

c	wer=wer*dt*3.d+0
	
	wer=
     *    -(txflux2-txflux1)*2.D+0
     *    -(tyflux2-tyflux1)*2.D+0
c    *    -(tzflux2*dzk(k)/(dzk(k)-lvl(m,n)))*2.D+0
     *	-(tzflux2-tzflux1)*2.D+0 

c	wer1=
c     *			cwz00(m,n)*dzk(k)/(dzk(k)-lvl(m,n))
c	wer1=wer1*dt*3.d+0

c	wer=
c     *    -(u2_-u1_)/dzk(k)*(dzk(k)-lvl(m,n))
c     *    -(v2_-v1_)/dzk(k)*(dzk(k)-lvl(m,n))
c     * -(ww(kp(m,n),k-1)*(ds1+ds2+ds3+ds4))

c	wer1=
c     *			cwz00(m,n)

	continue
	end if

165	continue

c=============END TEST ===========================================
c==============END  TEST ==========================================
c==============END  TEST ==========================================
c==============END  TEST ==========================================



	F1(M)=HGM(M,K)*TS(M,K)+wer+wer1

 	if (k.eq.2) then
       wer=
     *    -(sxflux(m,n,k)-sxflux(m-1,n,k))*2.D+0
     *    -(syflux(m,n,k)-syflux(m,n-1,k))*2.D+0
     *    -szflux(m,n,k)*2.D+0 *dzk(k)/(dzk(k)-lvl0(m,n))

	else
         wer=
     *    -(sxflux(m,n,k)-sxflux(m-1,n,k))*2.D+0
     *    -(szflux(m,n,k)-szflux(m,n,k-1))*2.D+0
     *    -(syflux(m,n,k)-syflux(m,n-1,k))*2.D+0
	end if
			if (k.eq.2) then
				wer1=
     *			3.D+0*dt*SS(m,k)*cwz00(m,n)*dzk(k)/(dzk(k)-lvl0(m,n))
			else
			 	wer1=0.
			end if

	if (m.eq.140.and.n.eq.137.and.k.eq.2) then
	print *,'potoki'
	print *, sxflux(m,n,k),sxflux(m-1,n,k)
	print *, syflux(m,n,k),syflux(m,n-1,k)
	print *, szflux(m,n,k),szflux(m,n,k-1)
	print *

	continue
	end if

c==============TEST TEST TEST    SSSSSSSSSSSS===========================================
c==============TEST TEST TEST    SSSSSSSSSSSS===========================================
c==============TEST TEST TEST    SSSSSSSSSSSS==========================================
c==============TEST TEST TEST    SSSSSSSSSSSS==========================================
	go to 166
	if (m.eq.60.and.n.eq.158.and.k.eq.4) then
	continue

       ds1=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.D+0
       ds2=(dlx(m-1,n)+dlx(m-1,n+1))*(dly(m,n)+dly(m-1,n))/4.D+0
       ds3=(dlx(m-1,n)+dlx(m-1,n-1))*(dly(m,n-1)+dly(m-1,n-1))/4.
       ds4=(dlx(m,n)+dlx(m,n-1))*(dly(m,n-1)+dly(m+1,n-1))/4.
         

	v2_=
     * 	(vv(kp(M,n),k-1)*(dlx(m,n)+dlx(m,n+1))
     *    +vv(kp(m-1,n),k-1)*(dlx(m-1,n)+dlx(m-1,n+1)))*dzk(k)
	
	v1_=
     * 	(vv(kp(M,n-1),k-1)*(dlx(m,n-1)+dlx(m,n))
     *    +vv(kp(m-1,n-1),k-1)*(dlx(m-1,n-1)+dlx(m-1,n)))*dzk(k)
	
	
      u2_=
     *	(uu(kp(M,n),k-1)*(dly(m,n)+dly(m+1,n))
     *    +uu(kp(m,n-1),k-1)*(dly(m,n-1)+dly(m+1,n-1)))*dzk(k)
	
      u1_=
     *	(uu(kp(M-1,n),k-1)*(dly(m-1,n)+dly(m,n))
     *    +uu(kp(m-1,n-1),k-1)*(dly(m-1,n-1)+dly(m,n-1)))*dzk(k)
	


	tzflux2=3.D+0*(ww(kp(m,n),k-1))*dt
     *	*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))/2.D+0
	tzflux1=3.D+0*(ww(kp(m,n),k-2))*dt
     *	*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))/2.D+0

	txflux1=u1_*dt*3.D+0/2.D+0
	txflux2=u2_*dt*3.D+0/2.D+0


	tyflux1=v1_*dt*3.D+0/2.D+0
	tyflux2=v2_*dt*3.D+0/2.D+0

	tzflux2=30.5783446479082*3.D+0*(ww(kp(m,n),k-1))*dt
     *	*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))/2.D+0
	tzflux1=30.3713164075531*3.D+0*(ww(kp(m,n),k-2))*dt
     *	*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))/2.D+0

	txflux1=30.3936926506469*u1_*dt*3.D+0/2.D+0
	txflux2=30.3463382969417*u2_*dt*3.D+0/2.D+0


	tyflux1=30.3976243660402*v1_*dt*3.D+0/2.D+0
	tyflux2=30.3404841142016*v2_*dt*3.D+0/2.D+0



	wer=
     *    -(u2_-u1_)
     *    -(v2_-v1_)
     * -(
     *  ww(kp(m,n),k-1)*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))
     * -ww(kp(m,n),k-2)*(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))
c   *dzk(k)/(dzk(k)-lvl(m,n))
     *  )

c	wer=wer*dt*3.d+0
	
	wer=
     *    -(txflux2-txflux1)*2.D+0
     *    -(tyflux2-tyflux1)*2.D+0
c    *    -(tzflux2*dzk(k)/(dzk(k)-lvl(m,n)))*2.D+0
     *	-(tzflux2-tzflux1)*2.D+0 

c	wer1=
c     *			cwz00(m,n)*dzk(k)/(dzk(k)-lvl(m,n))
c	wer1=wer1*dt*3.d+0

c	wer=
c     *    -(u2_-u1_)/dzk(k)*(dzk(k)-lvl(m,n))
c     *    -(v2_-v1_)/dzk(k)*(dzk(k)-lvl(m,n))
c     * -(ww(kp(m,n),k-1)*(ds1+ds2+ds3+ds4))

c	wer1=
c     *			cwz00(m,n)

	continue
	end if

166	continue

c=============END TEST ===========================================
c==============END  TEST ==========================================
c==============END  TEST ==========================================
c==============END  TEST ==========================================



	FS1(M)=HGM(M,K)*SS(M,K)+wer+wer1

      B1(M)=HGM(M,K)

	
	if (m.eq.140.and.n.eq.137.and.k.eq.2) then
	continue
	end if

	if (m.eq.73.and.n.eq.2.and.k.eq.4) then
	continue
	end if

  9   continue
         P2=0.
         ps2=0.
         p3=0.
         I1=M1
         I2=M2

      DO  M=M1,M2
c            SS2(M,K)=US1(M)
c            TS2(M,K)=U1(M)
            SS2(M,K)=FS1(m)/b1(m)
            TS2(M,K)=F1(m)/b1(m)
      end do      
    8 CONTINUE

      goto 2
  201 CONTINUE
                 DO  M=2,mh
                     K2=KP(M,N)
                 DO  K=kz0,KKZ

	if (isNaN(TS2(M,K))) then
	print *,'NaN temp!!!!', m,n,k, TS2(M,K)
	read *
	end if


c	if(SS2(M,K).lt.1000.and.SS2(M,K).lt.0) then
c        Write(nu_diag,*)'salt fixing',m,n,k, SS2(M,K)
c	ss2(m,k)=0.
c      read *
c       end if





	if(SS2(M,K).lt.1000.and.(SS2(M,K).ge.55.or.SS2(M,K).lt.0)) then
        Write(nu_diag,*)'salt!!!>55',m,n,k, SS2(M,K)
c      read *
c	 Write(nu_diag,*)'uu',m,n,k,uu(kp(m,n),K), uu(kp(m-1,n),K)
c     * , uu(kp(m,n-1),K), uu(kp(m-1,n-1),K)
c		 Write(nu_diag,*)'vv',m,n,k,vv(kp(m,n),K), vv(kp(m-1,n),K)
c     * , vv(kp(m,n-1),K), vv(kp(m-1,n-1),K)
c     read *
       end if




                   TEMP(K2,K-1)=TS2(M,K)
                   SAL(K2,K-1)=SS2(M,K)
                 end do
	            end do
   73 CONTINUE
C*************************************************
      RETURN
      END