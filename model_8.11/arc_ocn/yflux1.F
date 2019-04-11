          SUBROUTINE yflux1(Temp,Sal,uu,vv,ww,ww1,tbin,sbin,rivy,ftok)
C*************************************************
C*                                               *
C*      CROSYZ:   T,S(N=47)= T,S(KLIM)           *
C*                RU(T)=10*RU(S)                 *
C*************************************************
      Include 'model.par'
      Real*8 tbin(3), sbin(3)
      common/flux/tflux(mh,nh),sflux(mh,nh)
      real *8 tflux,sflux
      common/levitus/Tlev(nodes),Slev(nodes)
      real *8 tlev,slev
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
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      real *8 Temp(nodes,kgr),Sal(nodes,kgr),Ftok(mh,nh)
     *       ,uu(nodes,kgr),vv(nodes,kgr),ww(nodes,kgr),ww1(nodes,kgr)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/north/tatl(16,4,23),satl(16,4,23)
      common/clim_ts/tclim(mh,6,kgr),sclim(mh,6,kgr)
      real *8 tclim,sclim,vp,vl,s_riv
      integer nu_diag
      common/iodiag/nu_diag
c      REAL *8 SL,SM,SMB,SLB,DX,DY,Z,DZ,DZK
      real *8 A1(1000),B1(1000),C1(1000),F1(1000),U1(1000),
     *        FS1(1000),US1(1000)
      REAL *8 AS1(kh),BS1(kh),CS1(kh),PXS1(kh),CZS(nh,kh)
      REAL *8 TS(nh,kh),SS(nh,kh),hgm0(nh)
     *       ,ts1(nh,kh),SS1(nh,kh),ts2(nh,kh),SS2(nh,kh)
c     *	,tzflux(nh,kh),szflux(nh,kh)
	         common/adv_fluxes/ tyflux(mh,nh,kh),syflux(mh,nh,kh),
     *						txflux(mh,nh,kh),sxflux(mh,nh,kh)
     *						,tzflux(mh,nh,kh),szflux(mh,nh,kh)

     
      REAL *8 tyflux,syflux,txflux,sxflux,tzflux,szflux
c     *	,tzflux(nh,kh),szflux(nh,kh)
      real *8 rivy(mh,nh)

 

      REAL *8 HGM(nh,kh),CZ(nh,kh),CY(nh,kh),CVY(nh,kh),CVZ(nh,kh)
     *,        tsec(nh,kh),ssec(nh,kh),ga(nh,kh),asg(6)
      REAL *8 RZ1(nh,kh),RY(nh,kh),SU(nh,kh)
      REAL *8 DZ0,DZ1,ZZ1,ZZ2,Axx,Bxx,Cxx,dzgam,gk1,gk2,ghan,gk3,cnr
     *        ,ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
     *    ,gradl,gradls,curvl,curvls,sfll,sfl,sflls,sfls
     *        ,curz,curm
     *    ,gradup,gradups,curvz,curvzs,sflup,sfup,sflups,sfups
	REAL *8 scnor,slnor,cuplim,clwlim,sups,srnor
      data asg /10.,20.,50.,100.,200.,300./
C*********************************************
           ga_s=2.31e-7/5.
      ga_b=4.63e-8
      NZL=0
      LL=0
      KZ1=2
C**********************************************
      DO 200 M=M0,MM
      do 1 n=1,nh
      do  1 k=1,kh
c            umy1=abh+(ash-abh)*exp(-z(k)/50000.)
	tyflux(m,n,k)=0
	syflux(m,n,k)=0
	tzflux(m,n,k)=0
	szflux(m,n,k)=0
   1  continue
      NZL=NZL+1
          DO 110 N=1,Nh
      DX0=DX(M-1,n)
      DX1=DX(M,n)
             K2=KP(M,N)
           do k=1,kh
              CZS(N,K)=0.
             CY(N,K)=0.
            cy(n,k)=0.
            cvy(n,k)=0.
            cz(n,k)=0.
            cvz(n,k)=0.
          end do
           DO 110 K=2,KKZ
              tsec(n,k)=0
              ssec(n,k)=0
       ga(n,k)=0.
              TS(N,K)=TEMP(K2,K-1)
              SS(N,K)=SAL(K2,K-1)
c	SS(N,K)=1.
                  ts1(n,k)=ts(n,k)
                    ts2(n,k)=ts(n,k)
                       ss1(n,k)=ss(n,k)
                         ss2(n,k)=ss(n,k)
      

               HGM(n,K)=2./DT*DlX(M,n)*DlY(m,n)*DZ(K)

  110        CONTINUE
       do k=2,kh-1
       do n=2,2
          ga(n,k)=1./(360*3600.*24.*5)*0
       
          tsec(n,k)=tclim(m,n,k-1)
          ssec(n,k)=sclim(m,n,k-1)
       end do
       end do
c       end if
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
C      PRINT *,M,N1,N2
C      READ *
      NCROSS=NCROSS+1
      LN=LK+1
      NRA=0
      IF(N1.EQ.N2)GOTO2
         lk=l

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
c       Gam=((ds1*top1+ds2*top2+ds3*top3+ds4*top4)*dz0
c     *      +(ds1*bot1+ds2*bot2+ds3*bot3+ds4*bot4)*dz1)/8. 
CCCCCCCCCCC++++++++++++++++++++++++++++++++++++
      Gam=((ds1+ds2+ds3+ds4)*dz0
     *      +(ds1+ds2+ds3+ds4)*dz1)/8. 
      
           zs1=1
             dzgam=(dz0+dz1)/2.
          if(k.eq.ktet)dzgam=dz0/2.
          if(gam.eq.0) gam=ds1*dzgam
          if(k.eq.2)then
               zs1=0
             hgm0(n)=(ds1*top1+ds2*top2+ds3*top3+ds4*top4)/4.*3.
          end if
                      
      HGM(n,K)=6.*GAM*2.
      gradl=(ts(n,k)-ts(n-1,k))/dly(m,n-1)*inbt(m,n-1,k)*inbt(m,n,k)
       gradls=(ss(n,k)-ss(n-1,k))/dly(m,n-1)*inbt(m,n-1,k)*inbt(m,n,k)
      curm= (vv(kp(m-1,n-1),k-1)*(dlx(m-1,n-1)+dlx(m-1,n))
     *+vv(kp(m,n-1),k-1)*(dlx(m,n-1)+dlx(m,n)))
     */(dlx(m,n)+dlx(m,n-1)+dlx(m-1,n)+dlx(m-1,n-1))
     */dly(m,n-1)*inbt(m,n-1,k)*inbt(m,n,k)*dt


c--begin---
      if(curm.gt.0) then
	 icurm=inbt(m,n-1,k)*inbt(m,n-2,k)
         if(icurm.eq.0) then
	                     curvl=0.
	curvls=0
          curvl=2.D+0*((ts(n,k)-ts(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k)
c     *             -(ts(n-1,k)-ts(n-2,k))/dly(m,n-2)
c     *          *inbt(m,n-1,k)*inbt(m,n-2,k)
     *              )/(dly(m,n-1)
     *             +dly(m,n-2))
          curvls=2.D+0*((ss(n,k)-ss(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k)
c     *             -(ss(n-1,k)-ss(n-2,k))/dly(m,n-2)
c     *          *inbt(m,n-1,k)*inbt(m,n-2,k)
     *              )/(dly(m,n-1)
     *             +dly(m,n-2))
          	else
          curvl=2.D+0*((ts(n,k)-ts(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k)
     *             -(ts(n-1,k)-ts(n-2,k))/dly(m,n-2)
     *          *inbt(m,n-1,k)*inbt(m,n-2,k))
     *              /(dly(m,n-1)
     *             +dly(m,n-2))
          curvls=2.D+0*((ss(n,k)-ss(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k)
     *             -(ss(n-1,k)-ss(n-2,k))/dly(m,n-2)
     *          *inbt(m,n-1,k)*inbt(m,n-2,k))
     *              /(dly(m,n-1)
     *             +dly(m,n-2))
	    end if
      end if
c--end------

c===begin===
      if(curm.lt.0) then
	 icurm=inbt(m,n,k)*inbt(m,n+1,k)
         if(icurm.eq.0) then
	                     curvl=0
		curvls=0

          	else
          curvl=2.D+0*((ts(n+1,k)-ts(n,k))/dly(m,n)
     *          *inbt(m,n+1,k)*inbt(m,n,k)
     *             -(ts(n,k)-ts(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k))
     *              /(dly(m,n-1)+dly(m,n))
          curvls=2.D+0*((ss(n+1,k)-ss(n,k))/dly(m,n)
     *          *inbt(m,n+1,k)*inbt(m,n,k)
     *             -(ss(n,k)-ss(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k))
     *              /(dly(m,n-1)+dly(m,n))
	    end if
      end if
c===end=====




      sfll=0.5*(ts(n,k)+ts(n-1,k))
	sfl=0.5*(ts(n,k)+ts(n-1,k))
     *    -curm*gradl*dly(m,n-1)/2
     *     -(1-curm**2)*curvl*(dly(m,n-1)**2)/6.
      sflls=0.5*(ss(n,k)+ss(n-1,k))
	sfls=0.5*(ss(n,k)+ss(n-1,k))
     *    -curm*gradls*dly(m,n-1)/2
     *     -(1-curm**2)*curvls*(dly(m,n-1)**2)/6.


	if (m.eq.60.and.n.eq.158.and.k.eq.4) then
	continue
	end if

		if (m.eq.60.and.n.eq.159.and.k.eq.4) then
	continue
	end if


   	goto 765
c     ULTIMATE
c------------temp-----------------
       
      if(curm.gt.0) then
         ii=max(n1+1,n-2)
	sups=inbt(m,ii,k)*ts(ii,k)+(1-inbt(m,ii,k))*ts(n,k)

		if(abs(ts(n,k)-sups).lt.1.e-5) then
	  sfl=ts(n-1,k)


	else
	   scnor=(ts(n-1,k)-sups)/(ts(n,k)-sups)
	   srnor=(sfl-sups)/(ts(n,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(scnor/curm,1.)
              clwlim=scnor
                 if(srnor.gt.cuplim) srnor=cuplim
                 if(srnor.lt.clwlim) srnor=clwlim
          else
               srnor=scnor
          end if    
      sfl=sups+srnor*(ts(n,k)-sups)
	end if
	else
c#####
         ii=min(n2,n+1)
	sups=inbt(m,ii,k)*ts(ii,k)+(1-inbt(m,ii,k))*ts(n-1,k)
	
		if(abs(ts(n-1,k)-sups).lt.1.e-5) then
	  sfl=ts(n,k)

	else
	   scnor=(ts(n,k)-sups)/(ts(n-1,k)-sups)
	   slnor=(sfl-sups)/(ts(n-1,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(-scnor/curm,1.)
              clwlim=scnor
                 if(slnor.gt.cuplim) slnor=cuplim
                 if(slnor.lt.clwlim) slnor=clwlim
          else
               slnor=scnor
          end if    
      sfl=sups+slnor*(ts(n-1,k)-sups)
	end if
	end if
c--------------sal---------------------
	if(curm.gt.0) then
         ii=max(n1+1,n-2)
	sups=inbt(m,ii,k)*ss(ii,k)+(1-inbt(m,ii,k))*ss(n,k)
	
		if(abs(ss(n,k)-sups).lt.1.e-5) then
	  sfls=ss(n-1,k)


	else
	   scnor=(ss(n-1,k)-sups)/(ss(n,k)-sups)
	   srnor=(sfls-sups)/(ss(n,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(scnor/curm,1.)
              clwlim=scnor
                 if(srnor.gt.cuplim) srnor=cuplim
                 if(srnor.lt.clwlim) srnor=clwlim
          else
               srnor=scnor
          end if    
      sfls=sups+srnor*(ss(n,k)-sups)
	end if
	else
c#####
         ii=min(n2,n+1)
	sups=inbt(m,ii,k)*ss(ii,k)+(1-inbt(m,ii,k))*ss(n-1,k)

		if(abs(ss(n-1,k)-sups).lt.1.e-5) then
	  sfls=ss(n,k)

	else
	   scnor=(ss(n,k)-sups)/(ss(n-1,k)-sups)
	   slnor=(sfls-sups)/(ss(n-1,k)-sups)
      
          if(scnor.le.1.and.scnor.ge.0) then
              cuplim=min(-scnor/curm,1.)
              clwlim=scnor
                 if(slnor.gt.cuplim) slnor=cuplim
                 if(slnor.lt.clwlim) slnor=clwlim
          else
               slnor=scnor
          end if    
      sfls=sups+slnor*(ss(n-1,k)-sups)
	end if
	end if



 765  continue


      tyflux(m,n-1,k)=sfl*curm
     *      *6*dzk(k)*dly(m,n-1)*
     * 0.25D+0*(dlx(m,n)+dlx(m,n-1)+dlx(m-1,n)+dlx(m-1,n-1))



      syflux(m,n-1,k)=sfls*curm
     *      *6*dzk(k)*dly(m,n-1)*
     * 0.25*(dlx(m,n)+dlx(m,n-1)+dlx(m-1,n)+dlx(m-1,n-1))

	if (m.eq.135.and.n.eq.141.and.k.eq.2) then
	continue
	end if

			if (m.eq.140.and.n.eq.133.and.k.eq.2) then
	continue
	end if



   41 CONTINUE
   40 CONTINUE


c       go to 700
c******* Rivers +Bering Strait inflow
      n=n1
	if(rivy(m,n).ne.0.or.n.eq.2) then
	KTET=KBT(M,N)
         do K=KZ1,KTET

      if(n.eq.2) then
      curm= (vv(kp(m-1,n),k-1)*(dlx(m-1,n)+dlx(m-1,n+1))
     *       +vv(kp(m,n),k-1)*(dlx(m,n)+dlx(m,n+1)))
c     *       /(dlx(m,n)+dlx(m,n+1)+dlx(m-1,n)+dlx(m-1,n+1))
c     *       /dly(m,n)
     *       *inbt(m,n,k)*inbt(m,n+1,k)*dt

	    	if(curm.gt.0) then
                    sfls=0.5*(ss(n,k)+sclim(m,3,k-1))
                    sfl=0.5*(ts(n,k)+tclim(m,3,k-1))


              else
                sfls=0.5*(ss(n,k)+ss(n,k))
                sfl=0.5*(ts(n,k)+ts(n,k))
	    	end if

c=========================SOUTH BOUND LIKE RIVER=============
             vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n-1))
     *           *inbt(m,n,k)*inbt(m+1,n,k)
             vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n-1))
     *           *inbt(m,n,k)*inbt(m-1,n,k)

             curm= (vl*(dlx(m-1,n)+dlx(m-1,n-1))
     *             +vp*(dlx(m,n)+dlx(m,n-1)))*dt
     *            /(dlx(m,n-1)+dlx(m,n)+dlx(m-1,n-1)+dlx(m-1,n))
     *            /dly(m,n-1)

            	if(curm.gt.0) then
                    sfls=0.5*(ss(n,k)+sclim(m,3,k-1))
                    sfl=0.5*(ts(n,k)+tclim(m,3,k-1))
               else
                sfls=0.5*(ss(n,k)+ss(n,k))
                sfl=0.5*(ts(n,k)+ts(n,k))



           	end if
c=================================================================

      else
             vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n-1))
     *           *inbt(m,n,k)*inbt(m+1,n,k)
             vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n-1))
     *           *inbt(m,n,k)*inbt(m-1,n,k)

             curm= (vl*(dlx(m-1,n)+dlx(m-1,n-1))
     *             +vp*(dlx(m,n)+dlx(m,n-1)))*dt
     *            /(dlx(m,n-1)+dlx(m,n)+dlx(m-1,n-1)+dlx(m-1,n))
     *            /dly(m,n-1)

            sfl=0.5*(ts(n,k)+ts(n,k))
c	    sfl=0.5*(ts(n,k)+0.)

	s_riv=0.0

      curvls=2.D+0*((ss(n,k)-s_riv)/dly(m,n-1)-(s_riv-s_riv)/dly(m,n-2))
     *              /(dly(m,n-1)
     *             +dly(m,n-2))

	curvl=0.

	gradls=(ss(n,k)-s_riv)/dly(m,n-1)
	gradl=0.0

	sfls=0.5*(ss(n,k)+s_riv)
     *    -curm*gradls*dly(m,n-1)/2.
     *     -(1-curm**2)*curvls*(dly(m,n-1)**2)/6.
      sfl=0.5*(ts(n,k)+ts(n,k))


      end if

       tyflux(m,n-1,k)=sfl*curm
     *                *6*dzk(k)*0.25
     *        *dly(m,n-1)*(dlx(m,n-1)+dlx(m,n)+dlx(m-1,n-1)+dlx(m-1,n))
       syflux(m,n-1,k)=sfls*curm
     *         *6*dzk(k)*0.25
     *        *dly(m,n-1)*(dlx(m,n-1)+dlx(m,n)+dlx(m-1,n-1)+dlx(m-1,n))

	if (m.eq.107.and.n.eq.2.and.k.eq.4) then
	continue
	end if

		if (m.eq.140.and.n.eq.132.and.k.eq.2) then
	continue
	end if


         end do
	end if



      n=n1
	if(rivy(m,n).ne.0) then
	KTET=KBT(M,N)
        do K=KZ1,KTET

	s_riv=0.0
      gradls=(ss(n+1,k)-ss(n,k))/dly(m,n)*inbt(m,n+1,k)*inbt(m,n,k)
      gradl=(ts(n+1,k)-ts(n,k))/dly(m,n)*inbt(m,n+1,k)*inbt(m,n,k)
      curm= (vv(kp(m-1,n),k-1)*(dlx(m-1,n)+dlx(m-1,n+1))
     *+vv(kp(m,n),k-1)*(dlx(m,n)+dlx(m,n+1)))
     */(dlx(m,n)+dlx(m,n+1)+dlx(m-1,n)+dlx(m-1,n+1))
     */dly(m,n)*inbt(m,n+1,k)*inbt(m,n,k)*dt


      if(curm.gt.0) then
        
          curvls=2.D+0*((ss(n+1,k)-ss(n,k))/dly(m,n)
     *          *inbt(m,n+1,k)*inbt(m,n,k)
     *             -(ss(n,k)-s_riv)/dly(m,n-1)
     *              )/(dly(m,n-1)
     *             +dly(m,n))

	     curvl=2.D+0*((ts(n+1,k)-ts(n,k))/dly(m,n)
     *          *inbt(m,n+1,k)*inbt(m,n,k)
     *             -(ts(n,k)-ts(n,k))/dly(m,n-1)
     *              )/(dly(m,n-1)
     *             +dly(m,n))

          else
       curvls=0.
	 curvl=0.0
      end if


      if(curm.lt.0) then
	 icurm=inbt(m,n+2,k)*inbt(m,n+1,k)
         if(icurm.eq.0) then	                 
		curvls=0.
	    curvl=0.0
          	else
         curvls=2.D+0*((ss(n+2,k)-ss(n+1,k))/dly(m,n+1)
     *          *inbt(m,n+2,k)*inbt(m,n+1,k)
     *             -(ss(n+1,k)-ss(n,k))/dly(m,n)
     *          *inbt(m,n,k)*inbt(m,n+1,k))
     *              /(dly(m,n+1)+dly(m,n))
	   curvl=2.D+0*((ts(n+2,k)-ts(n+1,k))/dly(m,n+1)
     *          *inbt(m,n+2,k)*inbt(m,n+1,k)
     *             -(ts(n+1,k)-ts(n,k))/dly(m,n)
     *          *inbt(m,n,k)*inbt(m,n+1,k))
     *              /(dly(m,n+1)+dly(m,n))

	    end if
      end if
 
            sfl=0.5*(ts(n,k)+ts(n+1,k))
     *    -curm*gradl*dly(m,n)/2.
     *     -(1-curm**2)*curvl*(dly(m,n)**2)/6.

	sfls=0.5*(ss(n,k)+ss(n+1,k))
     *    -curm*gradls*dly(m,n)/2.
     *     -(1-curm**2)*curvls*(dly(m,n)**2)/6.

	 tyflux(m,n,k)=sfl*curm
     *         *6*dzk(k)*dly(m,n)*
     *         0.25*(dlx(m,n)+dlx(m,n+1)+dlx(m-1,n)+dlx(m-1,n+1))


      syflux(m,n,k)=sfls*curm
     *         *6*dzk(k)*0.25
     *        *dly(m,n)*(dlx(m,n+1)+dlx(m,n)+dlx(m-1,n+1)+dlx(m-1,n))

		if (m.eq.140.and.n.eq.131.and.k.eq.2) then
	continue
	end if
        end do
	end if



      n=n2
		if (m.eq.69.and.n.eq.181) then
	continue
	end if

	if(rivy(m,n).ne.0) then
         KTET=KBT(M,N)
         DO K=KZ1,KTET

             vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *           /(dlx(m,n)+dlx(m,n+1))
     *           *inbt(m,n,k)*inbt(m+1,n,k)
             vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *           /(dlx(m-1,n)+dlx(m-1,n+1))
     *          *inbt(m,n,k)*inbt(m-1,n,k)

         curm= (vl*(dlx(m-1,n+1)+dlx(m-1,n))
     *         +vp*(dlx(m,n+1)+dlx(m,n)))
     *         /(dlx(m,n)+dlx(m,n+1)+dlx(m-1,n)+dlx(m-1,n+1))
     *         /dly(m,n)*dt



c*********** Bering Strait  or not
      if(m.le.70.and.m.ge.68.and.n2.eq.nh-1)then
	s_riv=ss(n,k)
	else
	s_riv=0.0
	end if

          
         
      curvls=2.D+0*((s_riv-s_riv)/dly(m,n+1)-(s_riv-ss(n,k))/dly(m,n))
     *              /(dly(m,n+1)+dly(m,n))
	gradls=(s_riv-ss(n,k))/dly(m,n)
	curvl=0.0
	gradl=0.0

	sfl=0.5*(ts(n,k)+ts(n,k))
	sfls=0.5*(ss(n,k)+s_riv)
     *    -curm*gradls*dly(m,n)/2.
     *     -(1-curm**2)*curvls*(dly(m,n)**2)/6.


            tyflux(m,n,k)=sfl*curm
     *         *6*dzk(k)*dly(m,n)*
     *         0.25*(dlx(m,n)+dlx(m,n+1)+dlx(m-1,n)+dlx(m-1,n+1))


            syflux(m,n,k)=sfls*curm
     *         *6*dzk(k)*dly(m,n)*
     *         0.25*(dlx(m,n)+dlx(m,n+1)+dlx(m-1,n)+dlx(m-1,n+1))

	if (m.eq.69.and.n.eq.181.and.k.eq.2) then
	continue
	end if


      fghjhk=0
         ENDDO
      endif


	n=n2		
	if(rivy(m,n).ne.0) then
         KTET=KBT(M,N)
         DO K=KZ1,KTET    
 
	      if(m.le.70.and.m.ge.68.and.n2.eq.nh-1)then
	s_riv=ss(n,k)
		  else 
	s_riv=0.0 
		  end if   


	gradl=(ts(n,k)-ts(n-1,k))/dly(m,n-1)*inbt(m,n-1,k)*inbt(m,n,k)	        
      gradls=(ss(n,k)-ss(n-1,k))/dly(m,n-1)*inbt(m,n-1,k)*inbt(m,n,k)

      curm= (vv(kp(m-1,n-1),k-1)*(dlx(m-1,n-1)+dlx(m-1,n))
     *+vv(kp(m,n-1),k-1)*(dlx(m,n-1)+dlx(m,n)))
     */(dlx(m,n)+dlx(m,n-1)+dlx(m-1,n)+dlx(m-1,n-1))
     */dly(m,n-1)*inbt(m,n-1,k)*inbt(m,n,k)*dt

		   if(curm.lt.0) then              
		        curvl=2.D+0*((ts(n,k)-ts(n,k))/dly(m,n)
     *		            -(ts(n,k)-ts(n-1,k))/dly(m,n-1)
     *		        *inbt(m,n-1,k)*inbt(m,n,k))
     *		            /(dly(m,n-1)+dly(m,n))
	           
		        curvls=2.D+0*((s_riv-ss(n,k))/dly(m,n)
     *		            -(ss(n,k)-ss(n-1,k))/dly(m,n-1)
     *		        *inbt(m,n-1,k)*inbt(m,n,k))
     *		            /(dly(m,n-1)+dly(m,n))

			else 

			curvls=0.
	        curvl=0.0

			end if


	       if(curm.gt.0) then
			 icurm=inbt(m,n-1,k)*inbt(m,n-2,k)
			if(icurm.eq.0) then					
			 curvls=0.				
               curvls=2.D+0*((ss(n,k)-ss(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k)
c     *             -(ss(n-1,k)-ss(n-2,k))/dly(m,n-2)
c     *          *inbt(m,n-1,k)*inbt(m,n-2,k)
     *              )/(dly(m,n-1)
     *             +dly(m,n-2))

	         curvl=2.D+0*((ts(n,k)-ts(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k)
c     *             -(ss(n-1,k)-ss(n-2,k))/dly(m,n-2)
c     *          *inbt(m,n-1,k)*inbt(m,n-2,k)
     *              )/(dly(m,n-1)
     *             +dly(m,n-2))

          	 else

			       curvl=2.D+0*((ts(n,k)-ts(n-1,k))/dly(m,n-1)
     *					 *inbt(m,n-1,k)*inbt(m,n,k)
     *						 -(ts(n-1,k)-ts(n-2,k))/dly(m,n-2)
     *					 *inbt(m,n-1,k)*inbt(m,n-2,k))
     *					  /(dly(m,n-1)
     *					+dly(m,n-2))

		           curvls=2.D+0*((ss(n,k)-ss(n-1,k))/dly(m,n-1)
     *					 *inbt(m,n-1,k)*inbt(m,n,k)
     *						 -(ss(n-1,k)-ss(n-2,k))/dly(m,n-2)
     *					 *inbt(m,n-1,k)*inbt(m,n-2,k))
     *					  /(dly(m,n-1)
     *					+dly(m,n-2))
	          end if
                end if

    
	sfls=0.5*(ss(n,k)+ss(n-1,k))   
     *    -curm*gradls*dly(m,n-1)/2.
     *     -(1-curm**2)*curvls*(dly(m,n-1)**2)/6.

	    
	sfl=0.5*(ts(n,k)+ts(n-1,k))   
     *    -curm*gradl*dly(m,n-1)/2.
     *     -(1-curm**2)*curvl*(dly(m,n-1)**2)/6.



          syflux(m,n-1,k)=sfls*curm
     *      *6*dzk(k)*dly(m,n-1)*
     * 0.25*(dlx(m,n)+dlx(m,n-1)+dlx(m-1,n)+dlx(m-1,n-1))

	          tyflux(m,n-1,k)=sfl*curm
     *      *6*dzk(k)*dly(m,n-1)*
     * 0.25*(dlx(m,n)+dlx(m,n-1)+dlx(m-1,n)+dlx(m-1,n-1))





         ENDDO
      endif



  700 continue


      goto 2
 201  CONTINUE
  200 CONTINUE

C********************************************

      RETURN
      END