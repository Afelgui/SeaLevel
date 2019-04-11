       SUBROUTINE V_full_new(UU,VV,WW,WWvel,ftok,rivx,rivy,uu_ts,vv_ts)
      Include 'model.par'
      COMMON/NH2/ NH2
      REAL *8 SL,SLB,DX,DY,Z,DZ,DZK,h,h1,hh
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      integer *2 ikbt,inbt
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      real *8 UU(nodes,kgr),VV(nodes,kgr)
     *       ,WWvel(nodes,kgr),WW(nodes,kgr)
     *       ,ftok(mh,nh),ftok2(mh,nh)
     *       ,ubrt,vbrt,ubrt_ts(nodes),vbrt_ts(nodes)
     *         ,uu_ts(mh,nh,kgr),vv_ts(mh,nh,kgr),bbrrtt
     *         ,rivx(mh,nh),rivy(mh,nh)
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
	common/brt/UBrT(nodes), VBrT(nodes)
      real *8 ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
     * dz0,dz1,v1,v2,u1,u2,u3,u4,v3,vb,vt,up,ul
     *,v4,cwz,vp,vl,ut,ub,uctrl 
       do k=1,kgr
       do n=1,nodes
       WWvel(n,k)=0
       end do
       end do
      DO 14 K=1,kxlu
       N=IXDu(1,K)
        M1=IXDu(2,K)
         M2=IXDu(3,K)
      DO 14 M=M1,M2
          K2=KP(M,N)
      IF(H1(M,N).EQ.0)GOTO 14
c      UBrT(K2)=-1.D+0/H1(M,N)*
c     *           ((ftok(M,N+1)-ftok(M,N))+
c     *            (ftok(M+1,N+1)-ftok(M+1,N)))
c     *             *IPV(M,N)*IPH(M,N)/(dly(m,n)+dly(m+1,n))
  
c      VBrT(K2)=1.D+0/H1(M,N)*
c     *          ((ftok(M+1,N)-ftok(M,N))+
c     *            (ftok(M+1,N+1)-ftok(M,N+1)))
c     *             *IPH(M,N)*IPV(M,N)/(dlx(m,n)+dlx(m,n+1))

   14 CONTINUE


c      UBrT(1)=0.
c      VBrT(1)=0.
c	   UBrT_ts(1)=0.
c      VBrT(1)=0.

      DO 15 KL=1,kxlu
        N=IXDu(1,KL)
         M1=IXDu(2,KL)
          M2=IXDu(3,KL)
      DO 15 M=M1,M2
          K2=KP(M,N)
           KZV=KBTV(M,N)-1
      DO K=1,KZV
		if (m.eq.45.and.n.eq.30) then
	continue
	end if
            uu(k2,k)=(uu(k2,k)+ubrt(k2))
            vv(k2,k)=(vv(k2,k)+vbrt(k2))
      end do
           if(kzv.lt.kgr.and.kzv.gt.0) then
           do k=kzv+1,kgr
           uu(k2,k)=0
           vv(k2,k)=0.
           end do
           end if
  15  continue

	u_ts=0.
	v_ts=0.
c**********************************************
c		*		*		*		*		*
c		u2		u1		*		*		*
c		*		*		*		*	T	*
c		*		*		*		*		*
c************u_ts*****************u************
c		*		*		*		*		*
c		*		u4		*		*		*
c		u3    v_ts		*		*	T	*
c		*		*		*		*		*
c**********************************************
      DO 16 KL=1,kyline
        m=IyD(1,KL)
         n1=IyD(2,KL)
          n2=IyD(3,KL)
	n11=max(n1-1,2)  !!!ranshe 2 nachinat nelzya
      DO 16 n=n11,n2
	
          K2=KP(M,N)
           KZV=max(KBT(M,N)-1,KBT(M,N+1)-1)
		iv2=ikbt(m,n,2)
		iv1=ikbt(m,n-1,2)
		iu2=ikbt(m,n,2)
		iu1=ikbt(m-1,n,2)


      DO K=1,KZV


c*********************u2*********u1*******************	
c*****************          *          ***************	
c*****************          *          ***************
c*****************         u_ts        ***************	
c*****************          *          ***************
c*****************          *          ****************			
c*********************u3****T****u4******************	
         
		u1=(uu(kp(m,n),k)*(dly(m,n)+dly(m+1,n))
     *				+uu(kp(m,n+1),k)*(dly(m,n+1)+dly(m+1,n+1)))/4.D+0
		u2=(uu(kp(m-1,n),k)*(dly(m-1,n)+dly(m,n))
     *			+uu(kp(m-1,n+1),k)*(dly(m-1,n+1)+dly(m,n+1)))/4.D+0
		u3=(uu(kp(m-1,n-1),k)*(dly(m-1,n-1)+dly(m,n-1))
     *				+uu(kp(m-1,n),k)*(dly(m-1,n)+dly(m,n)))/4.D+0
		u4=(uu(kp(m,n-1),k)*(dly(m,n-1)+dly(m+1,n-1))
     *					+uu(kp(m,n),k)*(dly(m,n)+dly(m+1,n)))/4.D+0

	if (n.eq.2) then
	u3=0.
	u4=0.
	end if


		if (m.eq.66.and.n.eq.2.and.k.eq.30) then
	continue
	end if



		uu_ts(m,n,k)=(u1+u2+u3+u4)/4.d+0/dly(m,n)	


		if(rivx(m,n).ne.0.and.rivx(m,n+1).ne.0.and.iu2.eq.0) then

		up=-2.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           *inbt(m,n,k+1)*inbt(m,n+1,k+1)
		 uu_ts(m,n,k)=(u2+u3+2.*up)/4.d+0/dly(m,n)
		end if	


		if(rivx(m,n).ne.0.and.rivx(m,n+1).ne.0.and.iu1.eq.0) then

		ul=-2.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           *inbt(m,n,k+1)*inbt(m,n+1,k+1)
		 uu_ts(m,n,k)=(u1+u4+2.*ul)/4.d+0/dly(m,n)
		end if	


      end do
           if(kzv.lt.kgr.and.kzv.gt.0) then
           do k=kzv+1,kgr
	uu_ts(m,n,k)=0.
           end do
           end if
  16  continue



      DO 17 KL=1,kxline
        N=IXD(1,KL)
         M1=IXD(2,KL)
          M2=IXD(3,KL)
	M11=max(M1-1,2)
      DO 17 M=M1-1,M2
          K2=KP(M,N)
           
      KZV=max(KBT(M,N)-1,KBT(M+1,N)-1)
      DO K=1,KZV

		v1=(vv(kp(m,n),k)*(dlx(m,n)+dlx(m,n+1))
     *				+vv(kp(m+1,n),k)*(dlx(m+1,n)+dlx(m+1,n+1)))/4.D+0
		v2=(vv(kp(m,n-1),k)*(dlx(m,n-1)+dlx(m,n))
     *			+vv(kp(m+1,n-1),k)*(dlx(m+1,n-1)+dlx(m+1,n)))/4.D+0

	if (m.ne.1) then	
		v3=(vv(kp(m-1,n-1),k)*(dlx(m-1,n-1)+dlx(m-1,n))
     *				+vv(kp(m,n-1),k)*(dlx(m,n-1)+dlx(m,n)))/4.D+0

		v4=(vv(kp(m-1,n),k)*(dlx(m-1,n)+dlx(m-1,n+1))
     *				+vv(kp(m,n),k)*(dlx(m,n)+dlx(m,n+1)))/4.D+0
	else
		v3=0.
		v4=0.
	end if


	if (n.eq.2) then
	v2=v1
	v3=v4
	end if

	if (m.eq.2.and.n.eq.42) then
	continue
	end if

		vv_ts(m,n,k)=(v1+v2+v3+v4)/4.D+0/dlx(m,n)

	
		if(rivy(m,n).ne.0.and.rivy(m+1,n).ne.0.and.iv2.eq.0) then

		vt=2.D+0/(h(m,n)+h(m+1,n))*
     *           (ftok(M+1,N)-ftok(M,N))
     *           *inbt(m,n,k+1)*inbt(m+1,n,k+1)
		 vv_ts(m,n,k)=(v3+v4+2.*vt)/4.d+0/dlx(m,n)
		end if	


		if(rivy(m,n).ne.0.and.rivy(m+1,n).ne.0.and.iv1.eq.0) then

		vb=2.D+0/(h(m,n)+h(m+1,n))*
     *           (ftok(M+1,N)-ftok(M,N))
     *           *inbt(m,n,k+1)*inbt(m+1,n,k+1)
		 vv_ts(m,n,k)=(v1+v2+2.*vb)/4.d+0/dlx(m,n)
		end if	



      end do
           if(kzv.lt.kgr.and.kzv.gt.0) then
           do k=kzv+1,kgr
		vv_ts(m,n,k)=0.
           end do
           end if
  17  continue



      RETURN
      END



	 SUBROUTINE V_brt(UU,VV,WW,WWvel,lvl,lvl0)
      Include 'model.par'
      COMMON/NH2/ NH2
      REAL *8 SL,SLB,DX,DY,Z,DZ,DZK,h,h1
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
	   common/v5/grad_x(mh,nh),grad_y(mh,nh)
      integer *2 ikbt,inbt
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt


      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)

      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
		common/sum/lvl_sum(mh,nh),ubrt_sum(nodes),vbrt_sum(nodes)
      real *8 UU(nodes,kgr),VV(nodes,kgr)
     *       ,WWvel(nodes,kgr),WW(nodes,kgr)
     *       ,ftok(mh,nh),grad_x,grad_y
     *       ,ubrt,vbrt
      real *8 dlx,dly,lvl(mh,nh),lvl0(mh,nh),lvlx,lvly,mt(mh,nh)
      common/d3/dlx(mh,nh),dly(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
	common/brt/UBrT(nodes), VBrT(nodes)
      real *8 ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
     * dz0,dz1,v1,v2,u1,u2,cwz,vp,vl,ut,ub
     *         ,lvl_sum,ubrt_sum,vbrt_sum  

	real *8 u4,u3,wer,dx2,v3, v4, dy2

		g=9.8*100

	ctime=0.55
c	ctime=1
c	cgam=0.55
c	cgam=0
	cgam=1
c      dt=3*3600.
	

       do k=1,kgr
       do n=1,nodes
       WWvel(n,k)=0
       end do
       end do

      DO 15 K=1,kxlu
       N=IXDu(1,K)
        M1=IXDu(2,K)
         M2=IXDu(3,K)
      DO 15 M=M1,M2
          K2=KP(M,N)
      IF(H1(M,N).EQ.0)GOTO 15
C**********************************************************
C************************2*****1***************************
C***************************X******************************
C************************3*****4***************************
		mt(m,n)=1/(1+dt**2*ctime**2*slb(m,n)**2)*ikbt(m,n,2)

		u1=UBrT(K2)
		v1=VBrT(K2)

	A=mt(m,n)
     *   *((1-ctime*(1-ctime)*slb(m,n)**2*dt**2)*u1
     *                           +slb(m,n)*dt*v1)
	B=mt(m,n)
     *   *(-slb(m,n)*dt*u1
     *    +(1-ctime*(1-ctime)*slb(m,n)**2*dt**2)*v1)	


		G1=mt(m,n)
     *   *(dt*grad_x(m,n)+ctime*slb(m,n)*dt**2*grad_y(m,n))

       	G2=mt(m,n)
     *   *(-ctime*slb(m,n)*dt**2*grad_x(m,n)+dt*grad_y(m,n))

	
     
      r=mt(m,n)*dt
      q=r*dt*slb(m,n)*ctime
	G1=r*grad_x(m,n)+q*grad_y(m,n)
	G2=-q*grad_x(m,n)+r*grad_y(m,n)
			
		lvlx=cgam*(lvl(m+1,n)+lvl(m+1,n+1)-lvl(m,n)-lvl(m,n+1))
     *								/(dlx(m,n)+dlx(m,n+1))
     *	+(1-cgam)*(lvl0(m+1,n)+lvl0(m+1,n+1)-lvl0(m,n)-lvl0(m,n+1))
     *								/(dlx(m,n)+dlx(m,n+1))
    
	    lvly=cgam*(lvl(m,n+1)+lvl(m+1,n+1)-lvl(m,n)-lvl(m+1,n))
     *	                              /(dly(m,n)+dly(m+1,n))
     *   +(1-cgam)*(lvl0(m,n+1)+lvl0(m+1,n+1)-lvl0(m,n)-lvl0(m+1,n))
     *	                              /(dly(m,n)+dly(m+1,n))
  
	
	UBrT(K2)=A+G1+r*g*lvlx	+q*g*lvly
	VBrT(K2)=B+G2-q*g*lvlx+r*g*lvly
	UBrT(K2)=UBrT(K2)*ikbt(m,n,2)
	VBrT(K2)=VBrT(K2)*ikbt(m,n,2)

	if (m.eq.64.and.n.eq.163) then
c		print *
c		print *
	write(2,*),'  ' 
	write(2,*),'grad x',grad_x(m,n),'    grad_y', grad_y(m,n)
	write(2,*),'lvlx',lvlx,'    lvly',lvly
	write(2,*),'A',A,'G1',G1
	write(2,*),'r*g*lvlx+q*g*lvly',r*g*lvlx+q*g*lvly
        write(2,*),'B',B,'G2',G2
	  write(2,*),'-q*g*lvlx+r*g*lvly',-q*g*lvlx+r*g*lvly
	   write(2,*),'  ' 
	continue
	end if


   15 CONTINUE
      UBrT(1)=0.
      VBrT(1)=0.




	      DO 16 KL=1,kxline
        N=IXD(1,KL)
         M1=IXD(2,KL)
          M2=IXD(3,KL)
      DO 16 M=M1,M2
          K2=KP(M,N)


		u4=(ubrt(kp(m,n-1))*(dly(m,n-1)+dly(m+1,n-1))
     *					+ubrt(kp(m,n))*(dly(m,n)+dly(m+1,n)))
     *              /(dly(m,n-1)+dly(m+1,n-1)+dly(m,n)+dly(m+1,n))
		u4=u4* ((h(m,n)+h(m+1,n))/2.-(lvl(m,n)+lvl(m+1,n))/2. )
	    u3=(ubrt(kp(m-1,n-1))*(dly(m-1,n-1)+dly(m,n-1))
     *					+ubrt(kp(m-1,n))*(dly(m-1,n)+dly(m,n)))
     *              /(dly(m-1,n-1)+dly(m,n-1)+dly(m-1,n)+dly(m,n))
	
		u3=u3* ((h(m,n)+h(m-1,n))/2.-(lvl(m,n)+lvl(m-1,n))/2. )


		v4=(vbrt(kp(m-1,n))*(dlx(m-1,n)+dlx(m-1,n+1))
     *					+vbrt(kp(m,n))*(dlx(m,n)+dlx(m,n+1)))
     *              /(dlx(m-1,n)+dlx(m-1,n+1)+dlx(m,n)+dlx(m,n+1))
		v4=v4* ((h(m,n)+h(m,n+1))/2.-(lvl(m,n)+lvl(m,n+1))/2. )

		v3=(vbrt(kp(m-1,n-1))*(dlx(m-1,n-1)+dlx(m-1,n))
     *					+vbrt(kp(m,n-1))*(dlx(m,n-1)+dlx(m,n)))
     *              /(dlx(m-1,n-1)+dlx(m-1,n)+dlx(m,n-1)+dlx(m,n))
		v3=v3* ((h(m,n-1)+h(m,n))/2.-(lvl(m,n-1)+lvl(m,n))/2. )


	    dx2=(dlx(m,n)+dlx(m-1,n))/2.
	  dy2=(dly(m,n)+dly(m,n-1))/2.

          wer=(lvl(m,n)-lvl0(m,n))/dt- (
     *           	(u4-u3)/dx2
     *            +(v4-v3)/dy2
     *           	)

	if (wer.gt.1.e-7) then
c	print *,'wer',wer,'m',m,'n',n
c	read *
	end if 

  16  continue




      DO 17 KL=1,kxlu
        N=IXDu(1,KL)
         M1=IXDu(2,KL)
          M2=IXDu(3,KL)
      DO 17 M=M1,M2
          K2=KP(M,N)
           KZV=KBTV(M,N)-1
      DO K=1,KZV
c           uu(k2,k)=(uu(k2,k)+ubrt(k2))
c            vv(k2,k)=(vv(k2,k)+vbrt(k2))
      end do
           if(kzv.lt.kgr.and.kzv.gt.0) then
           do k=kzv+1,kgr
           uu(k2,k)=0.
           vv(k2,k)=0.
           end do
           end if

  17  continue

  

      RETURN
      END



      SUBROUTINE 
     *	wT_wV(UU,VV,ww,WW1,wwvel,WWvel1,ftok,rivx,rivy,lvl,lvl0)
      Include 'model.par'
      COMMON/NH2/ NH2
      REAL *8 SL,SLB,DX,DY,Z,DZ,DZK,h,h1
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      integer *2 ikbt,inbt
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ

      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      real *8 rivx(mh,nh),rivy(mh,nh),cwz00,cwzvel00
      real *8 UU(nodes,kgr),VV(nodes,kgr),WWvel(nodes,kgr)
     *       ,WWvel1(nodes,kgr),WW1(nodes,kgr),WW(nodes,kgr)
     *       ,ftok(mh,nh)
     *       ,glvlx,glvly
     *       ,ubrt,vbrt
     *       ,ur1,ur2,ur3,ur4,vr1,vr2,vr3,vr4
     *       ,vp,vl,ub,ut
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
	common/brt/UBrT(nodes), VBrT(nodes)
      real *8 ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
     * dz0,dz1,v1,v2,u1,u2,u3,u4,v3,v4,wer,wer1
     * ,cwz,cwz0,lvl(mh,nh),lvl0(mh,nh)  
     *  ,wd(nodes,kgr),wd1(nodes,kgr)
     * ,glvlx1,glvlx2,glvlx3,glvlx4
     * ,glvly1,glvly2,glvly3,glvly4,etat,hh
	common/cwz/   cwz00(mh,nh),	cwzvel00(mh,nh)


c	dt=3*3600.
       do k=1,kgr
       do n=1,nodes
       WWvel(n,k)=0
       end do
       end do


      do 20 Kl=1,kxLine
        N=IXD(1,Kl)
         M1=IXd(2,Kl)
          M2=IXD(3,Kl)
           DO 20 M=M1,M2
           K2=KP(M,N)
           kzv=kbt(m,n)-1
	if (m.eq.65.and.n.eq.18) then
	continue
	end if

           cwz=0 
           ds1=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.D+0
           ds2=(dlx(m-1,n)+dlx(m-1,n+1))*(dly(m,n)+dly(m-1,n))/4.D+0
           ds3=(dlx(m-1,n)+dlx(m-1,n-1))*(dly(m,n-1)+dly(m-1,n-1))/4.D+0
           ds4=(dlx(m,n)+dlx(m,n-1))*(dly(m,n-1)+dly(m+1,n-1))/4.D+0

c******** type of boundary for the river inflow ************
	iv2=ikbt(m,n,2)+ikbt(m-1,n,2)
	iv1=ikbt(m,n-1,2)+ikbt(m-1,n-1,2)
	iu2=ikbt(m,n,2)+ikbt(m,n-1,2)
	iu1=ikbt(m-1,n,2)+ikbt(m-1,n-1,2)

c	Glvlx=(lvl(m+1,n)-lvl(m-1,n))/(dlx(m,n)+dlx(m-1,n))
c	Glvly=(lvl(m,n+1)-lvl(m,n-1))/(dly(m,n)+dly(m,n-1))

c	Glvlx1=(lvl(m+1,n+1)+lvl(m+1,n)-lvl(m,n+1)-lvl(m,n))
c     *                            /(dlx(m,n)+dlx(m,n+1))
c	Glvlx2=(lvl(m,n+1)+lvl(m,n)-lvl(m-1,n+1)-lvl(m-1,n))
c     *                            /(dlx(m-1,n)+dlx(m-1,n+1))
c	Glvlx3=(lvl(m,n)+lvl(m,n-1)-lvl(m-1,n)-lvl(m-1,n-1))
c     *                            /(dlx(m-1,n-1)+dlx(m-1,n))
cc	Glvlx4=(lvl(m+1,n)+lvl(m+1,n-1)-lvl(m,n)-lvl(m,n-1))
c     *                            /(dlx(m,n-1)+dlx(m,n))

c	Glvly1=(lvl(m+1,n+1)+lvl(m,n+1)-lvl(m+1,n)-lvl(m,n))
c     *                            /(dly(m,n)+dly(m+1,n))

c	Glvly2=(lvl(m,n+1)+lvl(m-1,n+1)-lvl(m,n)-lvl(m-1,n))
c     *                            /(dly(m-1,n)+dly(m,n))

c	Glvly3=(lvl(m,n)+lvl(m-1,n)-lvl(m,n-1)-lvl(m-1,n-1))
c     *                            /(dly(m-1,n-1)+dly(m,n-1))
c	Glvly4=(lvl(m+1,n)+lvl(m,n)-lvl(m+1,n-1)-lvl(m,n-1))
c     *                            /(dly(m,n-1)+dly(m+1,n-1))


c	Glvlx=Glvlx1*ds1+Glvlx2*ds2+Glvlx3*ds3+Glvlx4*ds4
c	Glvly=Glvly1*ds1+Glvly2*ds2+Glvly3*ds3+Glvly4*ds4

c*************************************************************
c				*		*		*
c				*	u2	*	u1	*
c				*		*		*
c*************************************************************
c				*		*		*
c				*	u3	*	u4	*
c				*		*		*
c*************************************************************
	
c	u1=uu(kp(m,n),1)
c	u2=uu(kp(m-1,n),1)
c	u3=uu(kp(m-1,n-1),1)
c	u4=uu(kp(m,n-1),1)
c	v1=vv(kp(m,n),1)
c	v2=vv(kp(m-1,n),1)
c	v3=vv(kp(m-1,n-1),1)
c	v4=vv(kp(m,n-1),1)
c	u=(u1*ds1+u2*ds2+u3*ds3+u4*ds4)/(ds1+ds2+ds3+ds4)
c	v=(v1*ds1+v2*ds2+v3*ds3+v4*ds4)/(ds1+ds2+ds3+ds4)


c	cwz0=+(lvl(m,n)-lvl0(m,n))/dt*(ds1+ds2+ds3+ds4)
c     *                             + u*Glvlx+v*Glvly



c*************************************************************
c				*	    	*		    *
c				*	u2 B	*	u1 A	*
c				*	    	*		    *
c*************************************************************
c				*		    *		    *
c				*	u3 C    *	u4 D	*
c				*		    *		    *
c*************************************************************
c	goto 67

      k=1


	dx0=(dlx(m,n)+dlx(m-1,n))/2.
	dy0=(dly(m,n)+dly(m,n-1))/2.

	etat=+(lvl(m,n)-lvl0(m,n))/dt
     *  *(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))


	u1=uu(kp(M,n),k)*(dly(m,n)+dly(m+1,n))
	u2=uu(kp(M-1,n),k)*(dly(m-1,n)+dly(m,n))
	u3=uu(kp(m-1,n-1),k)*(dly(m-1,n-1)+dly(m,n-1))
	u4=uu(kp(m,n-1),k)*(dly(m,n-1)+dly(m+1,n-1))


	v1=vv(kp(M,n),k)*(dlx(m,n)+dlx(m,n+1))
	v2=vv(kp(M-1,n),k)*(dlx(m-1,n)+dlx(m-1,n+1))
	v3=vv(kp(M-1,n-1),k)*(dlx(m-1,n-1)+dlx(m-1,n))
	v4=vv(kp(M,n-1),k)*(dlx(m,n-1)+dlx(m,n))


      u1=u1*(+lvl0(m+1,n)+lvl0(m+1,n+1)+lvl0(m,n+1)-3*lvl0(m,n))/4.D+0
      u4=u4*(+lvl0(m+1,n-1)+lvl0(m+1,n)-3*lvl0(m,n)+lvl0(m,n-1))/4.D+0
      u3=u3*(+lvl0(m,n-1)-3*lvl0(m,n)+lvl0(m-1,n)+lvl0(m-1,n-1))/4.D+0
      u2=u2*(-3*lvl0(m,n)+lvl0(m,n+1)+lvl0(m-1,n+1)+lvl0(m-1,n))/4.D+0
      v1=v1*(+lvl0(m+1,n)+lvl0(m+1,n+1)+lvl0(m,n+1)-3*lvl0(m,n))/4.D+0
      v2=v2*(-3*lvl0(m,n)+lvl0(m,n+1)+lvl0(m-1,n+1)+lvl0(m-1,n))/4.D+0
      v3=v3*(+lvl0(m,n-1)-3*lvl0(m,n)+lvl0(m-1,n)+lvl0(m-1,n-1))/4.D+0
      v4=v4*(+lvl0(m+1,n-1)+lvl0(m+1,n)-3*lvl0(m,n)+lvl0(m,n-1))/4.D+0


	ur1=(h1(m,n)-
     *     (lvl0(m,n)+lvl0(m+1,n)+lvl0(m,n+1)+lvl0(m+1,n+1))/4.D+0)
     *   *ubrt(kp(m,n))*(dly(m,n)+dly(m+1,n))


	ur2=(h1(m-1,n)-
     *     (lvl0(m-1,n)+lvl0(m,n)+lvl0(m-1,n+1)+lvl0(m,n+1))/4.D+0)
     *   *ubrt(kp(m-1,n))*(dly(m-1,n)+dly(m,n))


	ur3=(h1(m-1,n-1)-
     *     (lvl0(m-1,n-1)+lvl0(m,n-1)+lvl0(m-1,n)+lvl0(m,n))/4.D+0)
     *   *ubrt(kp(m-1,n-1))*(dly(m-1,n-1)+dly(m,n-1))

	ur4=(h1(m,n-1)-
     *     (lvl0(m,n-1)+lvl0(m+1,n-1)+lvl0(m,n)+lvl0(m+1,n))/4.D+0)
     *   *ubrt(kp(m,n-1))*(dly(m,n-1)+dly(m+1,n-1))

	vr1=(h1(m,n)-
     *     (lvl0(m,n)+lvl0(m+1,n)+lvl0(m,n+1)+lvl0(m+1,n+1))/4.D+0)
     *   *vbrt(kp(m,n))*(dlx(m,n)+dlx(m,n+1))

	vr2=(h1(m-1,n)-
     *     (lvl0(m-1,n)+lvl0(m,n)+lvl0(m-1,n+1)+lvl0(m,n+1))/4.D+0)
     *   *vbrt(kp(m-1,n))*(dlx(m-1,n)+dlx(m-1,n+1))

	vr3=(h1(m-1,n-1)-
     *     (lvl0(m-1,n-1)+lvl0(m,n-1)+lvl0(m-1,n)+lvl0(m,n))/4.D+0)
     *   *vbrt(kp(m-1,n-1))*(dlx(m-1,n-1)+dlx(m-1,n))

	vr4=(h1(m,n-1)-
     *     (lvl0(m,n-1)+lvl0(m+1,n-1)+lvl0(m,n)+lvl0(m+1,n))/4.D+0)
     *   *vbrt(kp(m,n-1))*(dlx(m,n)+dlx(m,n-1))

c	goto 26

	if (n.eq.2) then
c		vr3=vr2
c		vr4=vr1
c		ur1=0
c		ur2=0
c		ur3=0
c		ur4=0
c		v3=v2
c		v4=v1
c		u1=0
c		u2=0
c		u3=0
c		u4=0
	end if

		if (m.eq.73.and.n.eq.2) then
	continue
	end if




	if(rivy(m,n).ne.0.and.iv2.eq.0) then

      vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n+1))
     *           *inbt(m,n,k+1)*inbt(m+1,n,k+1)
      vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n+1))
     *           *inbt(m,n,k+1)*inbt(m-1,n,k+1)

 	hh=min(h(m,n),h(m+1,n))
	vr1=(hh-
     *     (lvl0(m,n)+lvl0(m+1,n)+lvl0(m,n)+lvl0(m+1,n))/4.D+0)
     *   *vp*(dlx(m,n)+dlx(m,n+1))

	hh=min(h(m,n),h(m-1,n))
	vr2=(hh-
     *     (lvl0(m-1,n)+lvl0(m,n)+lvl0(m-1,n)+lvl0(m,n))/4.D+0)
     *   *vl*(dlx(m-1,n)+dlx(m-1,n+1))

	v1=vp*(dlx(m,n)+dlx(m,n+1))
	v2=vl*(dlx(m-1,n)+dlx(m-1,n+1))

	v1=v1*(+lvl0(m+1,n)+lvl0(m+1,n)+lvl0(m,n)-3*lvl0(m,n))/4.D+0
      v2=v2*(-3*lvl0(m,n)+lvl0(m,n)+lvl0(m-1,n)+lvl0(m-1,n))/4.D+0


     	end if


	if(rivy(m,n).ne.0.and.iv1.eq.0.or.n.eq.2) then

      vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n-1))
     *           *inbt(m,n,k+1)*inbt(m+1,n,k+1)
      vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n-1))
     *           *inbt(m,n,k+1)*inbt(m-1,n,k+1)
	hh=min(h(m,n),h(m-1,n))
	vr3=(hh-
     *     (lvl0(m-1,n)+lvl0(m,n)+lvl0(m-1,n)+lvl0(m,n))/4.D+0)
     *   *vl*(dlx(m-1,n-1)+dlx(m-1,n))

	hh=min(h(m,n),h(m+1,n))
		
		vr4=(hh-
     *     (lvl0(m,n)+lvl0(m+1,n)+lvl0(m,n)+lvl0(m+1,n))/4.D+0)
     *   *vp*(dlx(m,n)+dlx(m,n-1))

	v3=vl*(dlx(m-1,n-1)+dlx(m-1,n))
	v4=vp*(dlx(m,n-1)+dlx(m,n))

	v3=v3*(+lvl0(m,n)-3*lvl0(m,n)+lvl0(m-1,n)+lvl0(m-1,n))/4.D+0
      v4=v4*(+lvl0(m+1,n)+lvl0(m+1,n)-3*lvl0(m,n)+lvl0(m,n))/4.D+0

	end if



	if(rivx(m,n).ne.0.and.iu2.eq.0) then

c	lvl(m+1,n)=lvl(m,n)
c		lvl(m+1,n+1)=lvl(m,n+1)
c		lvl(m+1,n-1)=lvl(m,n-1)

	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m+1,n))
     *           *inbt(m,n,k+1)*inbt(m,n+1,k+1)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m+1,n-1))
     *           *inbt(m,n,k+1)*inbt(m,n-1,k+1)
	hh=min(h(m,n),h(m,n+1))
	ur1=(hh-
     *     (lvl0(m,n)+lvl0(m,n)+lvl0(m,n+1)+lvl0(m,n+1))/4.D+0)
     *   *ut*(dly(m,n)+dly(m+1,n))

	hh=min(h(m,n),h(m,n-1))


		ur4=(hh-
     *     (lvl0(m,n-1)+lvl0(m,n)+lvl0(m,n-1)+lvl0(m,n))/4.D+0)
     *   *ub*(dly(m,n-1)+dly(m+1,n-1))
	


	u1=ut*(dly(m,n)+dly(m+1,n))
	u1=u1*(+lvl0(m,n)+lvl0(m,n+1)+lvl0(m,n+1)-3*lvl0(m,n))/4.D+0

	u4=ub*(dly(m,n-1)+dly(m+1,n-1))
	u4=u4*(+lvl0(m,n-1)+lvl0(m,n)-3*lvl0(m,n)+lvl0(m,n-1))/4.D+0

	if (m.eq.128.and.n.eq.147) then
	continue
	end if


      end if


	if(rivx(m,n).ne.0.and.iu1.eq.0) then
	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m-1,n))
     *           *inbt(m,n,k+1)*inbt(m,n+1,k+1)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m-1,n-1))
     *           *inbt(m,n-1,k+1)*inbt(m,n,k+1)
	hh=min(h(m,n),h(m,n+1))
 	ur2=(hh-
     *     (lvl0(m,n)+lvl0(m,n)+lvl0(m,n+1)+lvl0(m,n+1))/4.D+0)
     *   *ut*(dly(m-1,n)+dly(m,n))
	hh=min(h(m,n),h(m,n-1))

	ur3=(hh-
     *     (lvl0(m,n-1)+lvl0(m,n-1)+lvl0(m,n)+lvl0(m,n))/4.D+0)
     *   *ub*(dly(m-1,n-1)+dly(m,n-1))

	u2=ut*(dly(m-1,n)+dly(m,n))
	u2=u2*(-3*lvl0(m,n)+lvl0(m,n+1)+lvl0(m,n+1)+lvl0(m,n))/4.D+0

	u3=ub*(dly(m-1,n-1)+dly(m,n-1))
	u3=u3*(+lvl0(m,n-1)-3*lvl0(m,n)+lvl0(m,n)+lvl0(m,n-1))/4.D+0
	end if

26      continue

	etat=(ur1+ur4-ur2-ur3)+(vr1+vr2-vr3-vr4)

c		etat=+(lvl(m,n)-lvl0(m,n))/dt
c    *  *(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))


	


	
	
	cwz0=etat
     * +u1
     * +u4
     * -u3
     * -u2
     * +v1
     * +v2
     * -v3
     * -v4






67		continue


c	cwz=cwz0/(dzk(2)-lvl(m,n))*dzk(2)
	cwz=cwz0
	cwz00(m,n)=cwz0

			if (isNan(cwz0)) then
	print *,m,n,cwz0
      read *
	end if

	
	 do k=1,kzv
                         DZ0=DZ(K)
                         DZ1=DZ(K+1)


      v2=
     * 	(vv(kp(M,n),k)*(dlx(m,n)+dlx(m,n+1))
     *    +vv(kp(m-1,n),k)*(dlx(m-1,n)+dlx(m-1,n+1)))*dzk(k+1)
	
	v1=
     * 	(vv(kp(M,n-1),k)*(dlx(m,n-1)+dlx(m,n))
     *    +vv(kp(m-1,n-1),k)*(dlx(m-1,n-1)+dlx(m-1,n)))*dzk(k+1)
	
	
      u2=
     *	(uu(kp(M,n),k)*(dly(m,n)+dly(m+1,n))
     *    +uu(kp(m,n-1),k)*(dly(m,n-1)+dly(m+1,n-1)))*dzk(k+1)
	
      u1=
     *	(uu(kp(M-1,n),k)*(dly(m-1,n)+dly(m,n))
     *    +uu(kp(m-1,n-1),k)*(dly(m-1,n-1)+dly(m,n-1)))*dzk(k+1)




	 if (k.eq.1) then
	u2=u2*(dzk(2)-lvl0(m,n))/dzk(2)
	u1=u1*(dzk(2)-lvl0(m,n))/dzk(2)
	v1=v1*(dzk(2)-lvl0(m,n))/dzk(2)
	v2=v2*(dzk(2)-lvl0(m,n))/dzk(2)
	end if



c	goto 34

cc	if(n.eq.nh-1) then
	if(rivy(m,n).ne.0.and.iv2.eq.0) then

      vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n+1))
     *           *inbt(m,n,k+1)*inbt(m+1,n,k+1)
      vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n+1))
     *           *inbt(m,n,k+1)*inbt(m-1,n,k+1)
      v2=
     * 	(vp*(dlx(m,n)+dlx(m,n+1))
     *    +vl*(dlx(m-1,n)+dlx(m-1,n+1)))*dzk(k+1)

	if (k.eq.1) then
	v2=v2*(dzk(2)-lvl0(m,n))/dzk(2)
	end if


	end if


	if(rivy(m,n).ne.0.and.iv1.eq.0.or.n.eq.2) then

      vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n-1))
     *           *inbt(m,n,k+1)*inbt(m+1,n,k+1)
      vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n-1))
     *           *inbt(m,n,k+1)*inbt(m-1,n,k+1)
      v1=
     * 	(vp*(dlx(m,n)+dlx(m,n-1))
     *    +vl*(dlx(m-1,n)+dlx(m-1,n-1)))*dzk(k+1)


			if (k.eq.1) then
	v1=v1*(dzk(2)-lvl0(m,n))/dzk(2)
	end if


	end if

	if(n.eq.2)  then
c	v2=v1
c	u2=0.
c	u1=0.
	end if

	if(rivx(m,n).ne.0.and.iu2.eq.0) then
	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m+1,n))
     *           *inbt(m,n,k+1)*inbt(m,n+1,k+1)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m+1,n-1))
     *           *inbt(m,n,k+1)*inbt(m,n-1,k+1)
      u2=
     *	(ut*(dly(m,n)+dly(m+1,n))
     *    +ub*(dly(m,n-1)+dly(m+1,n-1)))*dzk(k+1)

		if (k.eq.1) then
	u2=u2*(dzk(2)-lvl0(m,n))/dzk(2)
	end if

      end if


	if(rivx(m,n).ne.0.and.iu1.eq.0) then
	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m-1,n))
     *           *inbt(m,n,k+1)*inbt(m,n+1,k+1)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m-1,n-1))
     *           *inbt(m,n-1,k+1)*inbt(m,n,k+1)
	
      u1=
     *	(ut*(dly(m-1,n)+dly(m,n))
     *    +ub*(dly(m-1,n-1)+dly(m,n-1)))*dzk(k+1)

	if (k.eq.1) then
	u1=u1*(dzk(2)-lvl0(m,n))/dzk(2)
	end if

	end if

		if (m.eq.90.and.n.eq.2) then
	continue
	end if

34	 continue
c================for rivers we must do it too=========================



c=====================  eto w   ===================
       	if (m.eq.73.and.n.eq.2) then
	continue
	end if
	 cwz=cwz-(v2-v1)-(u2-u1)




c	u1=uu(kp(M,n),1)*(dly(m,n)+dly(m+1,n))

c		ur1=(h1(m,n)-
c     *     (lvl0(m,n)+lvl0(m+1,n)+lvl0(m,n+1)+lvl0(m+1,n+1))/4.D+0)
c     *   *ubrt(kp(m,n))*(dly(m,n)+dly(m+1,n))

	
   
c	etat=ur1

c	cwz0=(h1(m,n)-
c     *  (lvl0(m,n)+lvl0(m+1,n)+lvl0(m,n+1)+lvl0(m+1,n+1))/4.D+0)*u1
c     * +u1*(lvl0(m+1,n)+lvl0(m+1,n+1)+lvl0(m,n+1)-3*lvl0(m,n))/4.D+0



c	wer=cwz0-(h1(m,n)-lvl0(m,n))*u1


      if (m.eq.109.and.n.eq.86) then
	continue
	end if


c     * +wer1
      	WD(k2,k)=cwz/4.D+0 !!vertikalni potok
	if (m.eq.70.and.n.eq.70.and.k.eq.2) then
	continue
	end if
c        ww(k2,k)=cwz/(ds1+ds2+ds3+ds4)
        ww(k2,k)=cwz/(dlx(m,n)+dlx(m-1,n))/(dly(m,n)+dly(m,n-1))

	 if (ww(k2,kzv).gt.1.e-14.or.ww(k2,kzv).lt.-1.e-14) then
	print *, 'ww>1.e-14', ww(k2,k),m,n,k,kzv
	read *
	end if
c===================== eto omega  ===================
	if (k.eq.1) then
		ww1(k2,k)=ww(k2,k)*dzk(2)/(dzk(2)-lvl(m,n))
		WD1(k2,k)=WD(k2,k)*dzk(2)/(dzk(2)-lvl(m,n))
	else
		ww1(k2,k)=ww(k2,k)
		WD1(k2,k)=WD(k2,k)
	end if

           end do



      k=1

c	u2=uu(kp(M-1,n),k)*(dly(m-1,n)+dly(m,n))

c
c	ur2=(h1(m-1,n)-
c     *     (lvl0(m-1,n)+lvl0(m,n)+lvl0(m-1,n+1)+lvl0(m,n+1))/4.D+0)
c     *   *ubrt(kp(m-1,n))*(dly(m-1,n)+dly(m,n))

	etat=-ur2

c	cwz0=etat
c    * -u2*(3*lvl0(m,n)-lvl0(m,n+1)-lvl0(m-1,n+1)-lvl0(m-1,n))/4.D+0

c	wer1=-(h1(m,n)-lvl0(m,n))*u2

c	wer=cwz0-wer1

	

   20 continue 
   21  continue 

               
	do 30 Kl=1,kxLu
        N=IXDu(1,Kl)
         M1=IXdu(2,Kl)
          M2=IXDu(3,Kl)

           DO 30 M=M1-1,M2+1
           K2=KP(M,N)
           kzv=kbtv(m,n)-1
	kzv=max(kbt(m,n),kbt(m+1,n),kbt(m,n+1),kbt(m+1,n+1))-1    
		
           do k=1,kzv 
	WWvel(k2,k)=
     *	(wd(k2,k)+wd(kp(m,n+1),k)+wd(kp(m+1,n),k)+wd(kp(m+1,n+1),k))   
     *      /(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))

	WWvel1(k2,k)=
     *(wd1(k2,k)+wd1(kp(m,n+1),k)+wd1(kp(m+1,n),k)+wd1(kp(m+1,n+1),k))   
     *      /(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))

		if (m.eq.63.and.n.eq.12) then
	continue
	end if

		 if (wwvel(k2,kzv).gt.1.e-14.or.wwvel(k2,kzv).lt.-1.e-14) then
	print *, 'wwvel>1.e-14', wwvel(k2,k),m,n,k,kzv
	read *
	end if


	    
45	continue
           end do
   30 continue 


	go to 888

	do 32 Kl=1,kyLu
        m=IyDu(1,Kl)
         n1=Iydu(2,Kl)
          n2=IyDu(3,Kl)

         n=n1-1
           K2=KP(M,N)
	     kzv=max(kbt(m,n),kbt(m+1,n),kbt(m,n+1),kbt(m+1,n+1))-1
           cwz=0. 	
           do k=1,kzv
                 DZ0=DZ(K)
                 DZ1=DZ(K+1)
				v2=vv_ts(m,n+1,k)
     *                 *dlx(m,n+1)
 
			    v1=vv_ts(m,n,k)
     *                 *dlx(m,n)

				u2=uu_ts(m+1,n,k)
     *                 *dly(m+1,n)
 
				u1=uu_ts(m,n,k)
     *                 *dly(m,n)
		    	cwz=cwz
     *     			-(v2-v1)*dzk(k+1)
     *				-(u2-u1)*dzk(k+1)

		WWvel(k2,k)=
     *      4.D+0*cwz/(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))

c		WWvel(k2,k)=(wd(k2,k)+wd(kp(m,n+1),k)+wd(kp(m+1,n),k)+wd(kp(m+1,n+1),k))   
c     *      /(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))
	    



	     end do
	n=n2+1

           K2=KP(M,N)
      	 kzv=max(kbt(m,n),kbt(m+1,n),kbt(m,n+1),kbt(m+1,n+1))-1
           cwz=0. 	
           do k=1,kzv
                         DZ0=DZ(K)
                         DZ1=DZ(K+1)
				v2=vv_ts(m,n+1,k)
     *                 *dlx(m,n+1)
 
			    v1=vv_ts(m,n,k)
     *                 *dlx(m,n)

				u2=uu_ts(m+1,n,k)
     *                 *dly(m+1,n)
 
				u1=uu_ts(m,n,k)
     *                 *dly(m,n)
		    	cwz=cwz
     *     			-(v2-v1)*dzk(k+1)
     *				-(u2-u1)*dzk(k+1)

	        WWvel(k2,k)=
     *      	  4.D+0*cwz/(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))
c	WWvel(k2,k)=(wd(k2,k)+wd(kp(m,n+1),k)+wd(kp(m+1,n),k)+wd(kp(m+1,n+1),k))   
c     *      /(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))
	     end do

   32 continue 
  888 continue




   31  continue             
       lkjh=0

c	print *,'wwwww',ww(kp(60,3),1)
      RETURN
      END





      SUBROUTINE 
     *	wT_wV_bot(UU,VV,ww,WW1,wwvel,WWvel1,ftok,rivx,rivy,lvl,lvl0)
      Include 'model.par'
      COMMON/NH2/ NH2
      REAL *8 SL,SLB,DX,DY,Z,DZ,DZK,h,h1
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      integer *2 ikbt,inbt
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ

      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      real *8 rivx(mh,nh),rivy(mh,nh),cwz00,cwzvel00
      real *8 UU(nodes,kgr),VV(nodes,kgr),WWvel(nodes,kgr)
     *       ,WWvel1(nodes,kgr),WW1(nodes,kgr),WW(nodes,kgr)
     *       ,ftok(mh,nh)
     *       ,glvlx,glvly
     *       ,ubrt,vbrt
     *       ,ur1,ur2,ur3,ur4,vr1,vr2,vr3,vr4
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
	common/brt/UBrT(nodes), VBrT(nodes)
      real *8 ds1,ds2,ds3,ds4,top1,top2,top3,top4,bot1,bot2,bot3,bot4
     * dz0,dz1,v1,v2,u1,u2,u3,u4,v3,v4,wer,wer1
     * ,cwz,cwz0,lvl(mh,nh),lvl0(mh,nh)  
     *  ,wd(nodes,kgr),wd1(nodes,kgr)
     * ,glvlx1,glvlx2,glvlx3,glvlx4
     * ,glvly1,glvly2,glvly3,glvly4,etat
	common/cwz/   cwz00(mh,nh),	cwzvel00(mh,nh)


c	dt=3*3600.
       do k=1,kgr
       do n=1,nodes
       WWvel(n,k)=0
       end do
       end do


      do 20 Kl=1,kxLine
        N=IXD(1,Kl)
         M1=IXd(2,Kl)
          M2=IXD(3,Kl)
           DO 20 M=M1,M2
           K2=KP(M,N)
           kzv=kbt(m,n)-1
	if (m.eq.65.and.n.eq.18) then
	continue
	end if

           cwz=0 
           ds1=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))/4.D+0
           ds2=(dlx(m-1,n)+dlx(m-1,n+1))*(dly(m,n)+dly(m-1,n))/4.D+0
           ds3=(dlx(m-1,n)+dlx(m-1,n-1))*(dly(m,n-1)+dly(m-1,n-1))/4.D+0
           ds4=(dlx(m,n)+dlx(m,n-1))*(dly(m,n-1)+dly(m+1,n-1))/4.D+0

c******** type of boundary for the river inflow ************
	iv2=ikbt(m,n,2)+ikbt(m-1,n,2)
	iv1=ikbt(m,n-1,2)+ikbt(m-1,n-1,2)
	iu2=ikbt(m,n,2)+ikbt(m,n-1,2)
	iu1=ikbt(m-1,n,2)+ikbt(m-1,n-1,2)


c*************************************************************
c				*		*		*
c				*	u2	*	u1	*
c				*		*		*
c*************************************************************
c				*		*		*
c				*	u3	*	u4	*
c				*		*		*
c*************************************************************
	

c*************************************************************
c				*	    	*		    *
c				*	u2 B	*	u1 A	*
c				*	    	*		    *
c*************************************************************
c				*		    *		    *
c				*	u3 C    *	u4 D	*
c				*		    *		    *
c*************************************************************


67		continue


	cwz00(m,n)=cwz


	
	 do k=kzv,1,-1
                         DZ0=DZ(K)
                         DZ1=DZ(K+1)


      v2=
     * 	(vv(kp(M,n),k)*(dlx(m,n)+dlx(m,n+1))
     *    +vv(kp(m-1,n),k)*(dlx(m-1,n)+dlx(m-1,n+1)))*dzk(k+1)
	
	v1=
     * 	(vv(kp(M,n-1),k)*(dlx(m,n-1)+dlx(m,n))
     *    +vv(kp(m-1,n-1),k)*(dlx(m-1,n-1)+dlx(m-1,n)))*dzk(k+1)
	
	
      u2=
     *	(uu(kp(M,n),k)*(dly(m,n)+dly(m+1,n))
     *    +uu(kp(m,n-1),k)*(dly(m,n-1)+dly(m+1,n-1)))*dzk(k+1)
	
      u1=
     *	(uu(kp(M-1,n),k)*(dly(m-1,n)+dly(m,n))
     *    +uu(kp(m-1,n-1),k)*(dly(m-1,n-1)+dly(m,n-1)))*dzk(k+1)


	 if (k.eq.1) then
	u2=u2*(dzk(2)-lvl0(m,n))/dzk(2)
	u1=u1*(dzk(2)-lvl0(m,n))/dzk(2)
	v1=v1*(dzk(2)-lvl0(m,n))/dzk(2)
	v2=v2*(dzk(2)-lvl0(m,n))/dzk(2)
	end if

	goto 34

cc	if(n.eq.nh-1) then
	if(rivy(m,n).ne.0.and.iv2.eq.0) then

      vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n+1))
     *           *inbt(m,n,k+1)*inbt(m+1,n,k+1)
      vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n+1))
     *           *inbt(m,n,k+1)*inbt(m-1,n,k+1)
      v2=
     * 	(vp*(dlx(m,n)+dlx(m,n+1))
     *    +vl*(dlx(m-1,n)+dlx(m-1,n+1)))*dzk(k+1)

	end if

	if(rivy(m,n).ne.0.and.iv1.eq.0) then

      vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n-1))
     *           *inbt(m,n,k+1)*inbt(m+1,n,k+1)
      vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n-1))
     *           *inbt(m,n,k+1)*inbt(m-1,n,k+1)
      v1=
     * 	(vp*(dlx(m,n)+dlx(m,n-1))
     *    +vl*(dlx(m-1,n)+dlx(m-1,n-1)))*dzk(k+1)

	end if

	if(n.eq.2)  then
	v2=v1
	u2=0.
	u1=0.
	end if

	if(rivx(m,n).ne.0.and.iu2.eq.0) then
	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m+1,n))
     *           *inbt(m,n,k+1)*inbt(m,n+1,k+1)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m+1,n-1))
     *           *inbt(m,n,k+1)*inbt(m,n-1,k+1)
      u2=
     *	(ut*(dly(m,n)+dly(m+1,n))
     *    +ub*(dly(m,n-1)+dly(m+1,n-1)))*dzk(k+1)
      end if


	if(rivx(m,n).ne.0.and.iu1.eq.0) then
	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m-1,n))
     *           *inbt(m,n,k+1)*inbt(m,n+1,k+1)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m-1,n-1))
     *           *inbt(m,n-1,k+1)*inbt(m,n,k+1)
	
      u1=
     *	(ut*(dly(m-1,n)+dly(m,n))
     *    +ub*(dly(m-1,n-1)+dly(m,n-1)))*dzk(k+1)
	end if

34	 continue

c=====================  eto w   ===================
       	if (m.eq.57.and.n.eq.3.and.k.eq.2) then
	continue
	end if
	 cwz=cwz+(v2-v1)+(u2-u1)

      	WD(k2,k)=cwz/4.D+0 !!vertikalni potok
	if (m.eq.70.and.n.eq.70.and.k.eq.2) then
	continue
	end if
c        ww(k2,k)=cwz/(ds1+ds2+ds3+ds4)
        ww(k2,k)=cwz/(dlx(m,n)+dlx(m-1,n))/(dly(m,n)+dly(m,n-1))
c===================== eto omega  ===================
	if (k.eq.1) then
		ww1(k2,k)=ww(k2,k)*dzk(2)/(dzk(2)-lvl(m,n))
		WD1(k2,k)=WD(k2,k)*dzk(2)/(dzk(2)-lvl(m,n))
	else
		ww1(k2,k)=ww(k2,k)
		WD1(k2,k)=WD(k2,k)
	end if

           end do


   20 continue 
   21  continue 

               
	do 30 Kl=1,kxLu
        N=IXDu(1,Kl)
         M1=IXdu(2,Kl)
          M2=IXDu(3,Kl)

           DO 30 M=M1-1,M2+1
           K2=KP(M,N)
           kzv=kbtv(m,n)-1
	kzv=max(kbt(m,n),kbt(m+1,n),kbt(m,n+1),kbt(m+1,n+1))-1    
		
           do k=1,kzv 
	WWvel(k2,k)=
     *	(wd(k2,k)+wd(kp(m,n+1),k)+wd(kp(m+1,n),k)+wd(kp(m+1,n+1),k))   
     *      /(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))

	WWvel1(k2,k)=
     *(wd1(k2,k)+wd1(kp(m,n+1),k)+wd1(kp(m+1,n),k)+wd1(kp(m+1,n+1),k))   
     *      /(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))
	    
45	continue
           end do
   30 continue 


	go to 888

	do 32 Kl=1,kyLu
        m=IyDu(1,Kl)
         n1=Iydu(2,Kl)
          n2=IyDu(3,Kl)

         n=n1-1
           K2=KP(M,N)
	     kzv=max(kbt(m,n),kbt(m+1,n),kbt(m,n+1),kbt(m+1,n+1))-1
           cwz=0. 	
           do k=1,kzv
                 DZ0=DZ(K)
                 DZ1=DZ(K+1)
				v2=vv_ts(m,n+1,k)
     *                 *dlx(m,n+1)
 
			    v1=vv_ts(m,n,k)
     *                 *dlx(m,n)

				u2=uu_ts(m+1,n,k)
     *                 *dly(m+1,n)
 
				u1=uu_ts(m,n,k)
     *                 *dly(m,n)
		    	cwz=cwz
     *     			-(v2-v1)*dzk(k+1)
     *				-(u2-u1)*dzk(k+1)

		WWvel(k2,k)=
     *      4.D+0*cwz/(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))

c		WWvel(k2,k)=(wd(k2,k)+wd(kp(m,n+1),k)+wd(kp(m+1,n),k)+wd(kp(m+1,n+1),k))   
c     *      /(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))
	    



	     end do
	n=n2+1

           K2=KP(M,N)
      	 kzv=max(kbt(m,n),kbt(m+1,n),kbt(m,n+1),kbt(m+1,n+1))-1
           cwz=0. 	
           do k=1,kzv
                         DZ0=DZ(K)
                         DZ1=DZ(K+1)
				v2=vv_ts(m,n+1,k)
     *                 *dlx(m,n+1)
 
			    v1=vv_ts(m,n,k)
     *                 *dlx(m,n)

				u2=uu_ts(m+1,n,k)
     *                 *dly(m+1,n)
 
				u1=uu_ts(m,n,k)
     *                 *dly(m,n)
		    	cwz=cwz
     *     			-(v2-v1)*dzk(k+1)
     *				-(u2-u1)*dzk(k+1)

	        WWvel(k2,k)=
     *      	  4.D+0*cwz/(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))
c	WWvel(k2,k)=(wd(k2,k)+wd(kp(m,n+1),k)+wd(kp(m+1,n),k)+wd(kp(m+1,n+1),k))   
c     *      /(dlx(m,n+1)+dlx(m,n))/(dly(m+1,n)+dly(m,n))
	     end do

   32 continue 
  888 continue




   31  continue             
       lkjh=0

c	print *,'wwwww',ww(kp(60,3),1)
      RETURN
      END
