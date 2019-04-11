          SUBROUTINE tras1(nt,Sal,uu,vv,ww,tbin,sbin,rivy,ftok)
C*************************************************
C*                                               *
C*      CROSYZ:   T,S(N=47)= T,S(KLIM)           *
C*                RU(T)=10*RU(S)                 *
C*************************************************
      Include 'model.par'
      Real*8 tbin(3), sbin(3)
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
      real *8 Sal(nodes,kgr),Ftok(mh,nh)
     *       ,uu(nodes,kgr),vv(nodes,kgr),ww(nodes,kgr)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/north/tatl(16,4,23),satl(16,4,23)
      common/clim_ts/tclim(mh,6,kgr),sclim(mh,6,kgr)
      real *8 tclim,sclim,vp,vl
      integer nu_diag
      common/iodiag/nu_diag
c      REAL *8 SL,SM,SMB,SLB,DX,DY,Z,DZ,DZK
      real *8 A1(1000),B1(1000),C1(1000),F1(1000),U1(1000),
     *        FS1(1000),US1(1000)
      REAL *8 AS1(kh),BS1(kh),CS1(kh),PXS1(kh),CZS(nh,kh)
      REAL *8 SS(nh,kh),hgm0(nh)
     *       ,ts1(nh,kh),SS1(nh,kh),SS2(nh,kh)
     *	,tzflux(nh,kh),szflux(nh,kh)
      common/adv_fluxes/ txflux(mh,nh,kh),sxflux(mh,nh,kh)
      REAL *8 txflux,sxflux
      real *8 rivy(mh,nh)
      integer tracx(mh,nh),tracxy(mh,nh)
	Common/Xtrac/tracx,tracxy

 

      REAL *8 HGM(nh,kh),CZ(nh,kh),CY(nh,kh),CVY(nh,kh),CVZ(nh,kh)
     *,        tsec(nh,kh),ssec(nh,kh),ga(nh,kh),asg(6)
      REAL *8 RZ1(nh,kh),RY(nh,kh),SU(nh,kh)
      REAL *8 DZ0,DZ1,ZZ1,ZZ2,Axx,Bxx,Cxx,dzgam,gk1,gk2,ghan,gk3,cnr
     *        ,ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
     *    ,gradl,gradls,curvl,curvls,sfll,sfl,sflls,sfls
     *        ,curz,curm
     *    ,gradup,gradups,curvz,curvzs,sflup,sfup,sflups,sfups
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
	sxflux(m,n,k)=0
	szflux(n,k)=0
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
              ssec(n,k)=0
       ga(n,k)=0.
              SS(N,K)=SAL(K2,K-1)
                       ss1(n,k)=ss(n,k)
                         ss2(n,k)=ss(n,k)
      

               HGM(n,K)=2./DT*DlX(M,n)*DlY(m,n)*DZ(K)

  110        CONTINUE
       do k=2,kh-1
       do n=2,2
          ga(n,k)=1./(360*3600.*24.*5)*0
       
          ssec(n,k)=0+0*sclim(m,n,k-1)
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
       gradls=(ss(n,k)-ss(n-1,k))/dly(m,n-1)*inbt(m,n-1,k)*inbt(m,n,k)
      curm= (vv(kp(m-1,n-1),k-1)*(dlx(m-1,n-1)+dlx(m-1,n))
     *+vv(kp(m,n-1),k-1)*(dlx(m,n-1)+dlx(m,n)))
     */(dlx(m,n)+dlx(m,n-1)+dlx(m-1,n)+dlx(m-1,n-1))
     */dly(m,n-1)*inbt(m,n-1,k)*inbt(m,n,k)*dt
c--begin---
      if(curm.gt.0) then
	 icurm=inbt(m,n-1,k)*inbt(m,n-2,k)
         if(icurm.eq.0) then
	                     curvl=0
	curvls=0
          curvls=2*((ss(n,k)-ss(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k)
c     *             -(ss(n-1,k)-ss(n-2,k))/dly(m,n-2)
c     *          *inbt(m,n-1,k)*inbt(m,n-2,k)
     *              )/(dly(m,n-1)
     *             +dly(m,n-2))
          	else
          curvls=2*((ss(n,k)-ss(n-1,k))/dly(m,n-1)
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
          curvls=2*((ss(n+1,k)-ss(n,k))/dly(m,n)
     *          *inbt(m,n+1,k)*inbt(m,n,k)
     *             -(ss(n,k)-ss(n-1,k))/dly(m,n-1)
     *          *inbt(m,n-1,k)*inbt(m,n,k))
     *              /(dly(m,n-1)+dly(m,n))
	    end if
      end if
c===end=====
      sflls=0.5*(ss(n,k)+ss(n-1,k))
	sfls=0.5*(ss(n,k)+ss(n-1,k))
     *    -curm*gradls*dly(m,n-1)/2
     *     -(1-curm**2)*curvls*(dly(m,n-1)**2)/6.


      sxflux(m,n-1,k)=sfls*curm
     *      *6*dzk(k)*dly(m,n-1)*
     * 0.25*(dlx(m,n)+dlx(m,n-1)+dlx(m-1,n)+dlx(m-1,n-1))

   41 CONTINUE
   40 CONTINUE


c       go to 700
c******* Rivers +Bering Strait inflow
      n=n1
	if(rivy(m,n).ne.0.or.n.eq.2) then
	                     strac=1.
                           if(nt.eq.tracx(m,n))strac=100.
	                     if(nt.eq.tracxy(m,n))strac=100.
c	                         if(nt.eq.tracx(m,n))strac=1.
	KTET=KBT(M,N)
         do K=KZ1,KTET

      if(n.eq.2) then
      curm= (vv(kp(m-1,n),k-1)*(dlx(m-1,n)+dlx(m-1,n+1))
     *       +vv(kp(m,n),k-1)*(dlx(m,n)+dlx(m,n+1)))
c     *       /(dlx(m,n)+dlx(m,n+1)+dlx(m-1,n)+dlx(m-1,n+1))
c     *       /dly(m,n)
     *       *inbt(m,n,k)*inbt(m,n+1,k)*dt

	if(curm.gt.0) then
                    sfls=0.5*(ss(n,k)+1.)
           else
                sfls=0.5*(ss(n,k)+ss(n,k))
	end if

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
c     *            /(dlx(m,n-1)+dlx(m,n)+dlx(m-1,n-1)+dlx(m-1,n))
c     *            /dly(m,n-1)

            sfls=0.5*(ss(n,k)+strac)
      end if

       sxflux(m,n-1,k)=sfls*curm
     *         *6*dzk(k)*0.25
c     *         *dly(m,n-1)*(dlx(m,n-1)+dlx(m,n)+dlx(m-1,n-1)+dlx(m-1,n))
         end do
	end if

      n=n2
	if(rivy(m,n).ne.0) then
	                     strac=1.
c	                    if(nt.eq.tracx(m,n))strac=1.
                           if(nt.eq.tracx(m,n))strac=100.
	                       if(nt.eq.tracxy(m,n))strac=100.
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



  


            sfls=0.5*(ss(n,k)+strac)





            sxflux(m,n,k)=sfls*curm
     *         *6*dzk(k)*dly(m,n)*
     *         0.25*(dlx(m,n)+dlx(m,n+1)+dlx(m-1,n)+dlx(m-1,n+1))
         ENDDO
      endif

  700 continue


      goto 2
 201  CONTINUE
  200 CONTINUE

C********************************************

      RETURN
      END
