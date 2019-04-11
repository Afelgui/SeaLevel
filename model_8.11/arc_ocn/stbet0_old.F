       SUBROUTINE STBET00(lvl,lvl0,Gam,rivx,rivy,ftok) 
      Include 'model.par'
      COMMON/DAY/myear,month,nday
      COMMON/NH2/ NH2
      real *8 z,Dz,Dzk,dx,dy,h,h1
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt

	    

     	Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      common/clame/hl1(mh,nh),hl2(mh,nh)
      real *8 Dzet(mh,nh),lvl(mh,nh),eta(mh,nh),lvl0(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      COMMON/V1/ROZ(mh,nh)
		common/sum/lvl_sum(mh,nh),ubrt_sum(nodes),vbrt_sum(nodes)
	common/brt/UBrT(nodes), VBrT(nodes)
	     common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
	common/v5/grad_x(mh,nh),grad_y(mh,nh)
	  common/d3/dlx(mh,nh),dly(mh,nh)
	 real *8 dlx,dly
      common/Wind/Tax(mh,nh),Tay(mh,nh)
     
      COMMON/V4/Risl2
      REAL *8 ROZ,FANT,Risl2,f00,f000,gg5,fg1,gg1,fg2,gg2(kk_isl)
       REAL *8 Q1,Q2,Q3,Q4,Q5,Q6,DTR
     *	 ,DX0,DX1,DY0,DY1,dx2,dy2,dx4,dy4,dx5,dy5,hl1,hl2,sl,slb
      Real *8 A1(1000),B1(1000),C1(1000),FF1(1000),U1(1000)
      REAL *8 ccx,acx,ccy,acy,cdx,adx,cdy,ady,cnx,cny,anx,any
     *         ha,ha1,hb,hb1,dly1,dlx1,f1,f2,f4,cgam,acd
     *         ,ct1,ct2,ct3,ct4,ct5
     * ,ct6,ct7,ct8,ct9,AB,GG
      REAL *8 GAM(mh,nh),tok_n(mh),Flh(mh,nh), u_br(mh,nh), v_br(mh,nh)
     *	,etan,rg(mh,nh),q(mh,nh),hmt(mh,nh),AA(mh,nh),BB(mh,nh)
     *,G1(mh,nh),G2(mh,nh),lvlx(mh,nh),lvly(mh,nh),cx(mh,nh),cy(mh,nh)
     * ,lvl_uv
	real *8 lvl1,lvl2,lvl3,lvl4,lvl5,lvl6,lvl7,lvl8,lvl9
     *	   ,lvlx1,lvlx2,lvlx4,lvlx5,lvly1,lvly2,lvly4,lvly5
     *	   ,r1,r2,r4,r5,a_1,a_2,a_4,a_5,b_1,b_2,b_4,b_5
      real *8 rivx(mh,nh),rivy(mh,nh),ftok(mh,nh),vbrt,ubrt
     *  ,grad_x,grad_y,lvl_sum,tax,tay,hgm,hh,Fus,vl,vp,ub,ut
	      integer nu_diag
	 integer *2 ikbt,inbt
      common/iodiag/nu_diag
      real *8 cnr,ctime
      SAVE aw1,aw2,aw3
      Data aw1,aw2,aw3/0.,0.,0./
      dimension NPsi(mh,nh)
      cnr=0.
	g=9.8*100

	ctime=0.55
c	ctime=1
	cgam=1
c	cgam=0
	raz=0
c      dt=3*3600/1.
	
c	goto 108

	do n=2,nh-1
      do m=2,mh-1
      k2=kp(m,n)
      u_br(m,n)=ubrt(k2)
      v_br(m,n)=vbrt(k2)
		ikbt(m,n,2)=inbt(m,n,2)*inbt(m+1,n,2)*inbt(m,n+1,2)
     *                                *inbt(m+1,n+1,2)

c******** type of boundary for the river inflow ************

	iv2=ikbt(m,n,2)+ikbt(m-1,n,2)
	iv1=ikbt(m,n-1,2)+ikbt(m-1,n-1,2)
	iu2=ikbt(m,n,2)+ikbt(m,n-1,2)
	iu1=ikbt(m-1,n,2)+ikbt(m-1,n-1,2)

	if(rivy(m,n).ne.0.and.rivy(m+1,n).ne.0.and.iv2.eq.0) then
	v_br(m,n)=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n+1))
     *           *inbt(m,n,2)*inbt(m+1,n,2)
	end if

	if(rivy(m,n).ne.0.and.rivy(m+1,n).ne.0.and.iv1.eq.0) then
	v_br(m,n-1)=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n-1)+dlx(m,n))
     *           *inbt(m,n,2)*inbt(m+1,n,2)
	end if

	if(rivx(m,n).ne.0.and.rivx(m,n+1).ne.0.and.iu2.eq.0) then
	u_br(m,n)=4.D+0/(h(m,n)+h(m,n+1))*
     *          (ftok(M,N+1)-ftok(M,N))
     *          /(dly(m,n)+dly(m+1,n))
     *           *inbt(m,n,2)*inbt(m,n+1,2)
	end if

	if(rivx(m,n).ne.0.and.rivx(m,n+1).ne.0.and.iu1.eq.0) then
	u_br(m-1,n)=4.D+0/(h(m,n)+h(m,n+1))*
     *          (ftok(M,N+1)-ftok(M,N))
     *          /(dlx(m-1,n)+dlx(m,n))
     *           *inbt(m,n,2)*inbt(m,n+1,2)
	end if





	end do
	end do




        do n=1,nh
        do m=1,mh

   
	if (m.ne.mh.and.n.ne.nh) then

c	lvl_uv=(lvl(m,n)+lvl(m+1,n)+lvl(m,n+1)+lvl(m+1,n+1))/4.D+0
c	hmt(m,n)=(H1(m,n)-lvl_uv)/(1+dt**2*ctime**2*slb(m,n)**2)
c     *	*ikbt(m,n,2)


c      eta(m,n)=lvl(m,n)
c	  lvl0(m,n)=lvl(m,n)
c	hmt(m,n)=H1(m,n)/(1+dt**2*ctime**2*slb(m,n)**2)*ikbt(m,n,1)

	end if
        
        end do
        end do


      DTR=1./DT
      ENERGY=0.





      DO 888 K=1,kxLu
      n=IxDu(1,K)
      I1=IxDu(2,K)
      I2=IxDu(3,K)

	if (n.eq.179) then
	continue
	end if


	m11=max(i1,2)
      DO 88 m=m11,I2    
c*****************************************************

c*****************************************************************************
c***********************7-----8-----9*****************************************
c***********************|**4**|**5**|*****************************************
c***********************4-----5-----6*****************************************
c***********************|**1**|**2**|*****************************************
c***********************1-----2-----3*****************************************	
cc*****************************************************************************
	lvl_uv=(lvl(m,n)+lvl(m+1,n)+lvl(m,n+1)+lvl(m+1,n+1))/4.D+0


c	grad1_x(m,n)=grad_x(m,n)+tax(m,n)/(h1(m,n)-lvl_uv)
c	grad1_y(m,n)=grad_y(m,n)+tay(m,n)/(h1(m,n)-lvl_uv)

c	grad_x(m,n)=tax(m,n)/(h1(m,n)-lvl_uv)
c	grad_y(m,n)=tay(m,n)/(h1(m,n)-lvl_uv)

	lvlx(m,n)=(lvl(m+1,n)+lvl(m+1,n+1)-lvl(m,n)-lvl(m,n+1))    			
	lvly(m,n)=(lvl(m,n+1)+lvl(m+1,n+1)-lvl(m,n)-lvl(m+1,n))

	if (m.eq.67.and.n.eq.180) then
	continue
	end if

c*****************************************************************************
c***********************7-----8-----9*****************************************
c***********************|**4**|**5**|*****************************************
c***********************4-----5-----6*****************************************
c***********************|**1**|**2**|*****************************************
c***********************1-----2-----3*****************************************	
c*****************************************************************************

   88 CONTINUE
 888  continue


      
	DO 988 K=1,kxLine
      n=IxD(1,K)
      I1=IxD(2,K)
      I2=IxD(3,K)
      DO 98 m=I1,I2

	 eta(m,n)=lvl(m,n)
	 lvl0(m,n)=lvl(m,n)
	
c-----------------------------

   98 CONTINUE
 988  continue

c
c**** main contur******************
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      f00=0
      fg2=0

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     gg2(1)=0
      DO4 N=1,1000
      A1(N)=0.
      B1(N)=0.
      C1(N)=0.
      U1(N)=0.
    4 FF1(N)=0.
      KK=0                                 
    2 KK=KK+1
c      print *,'iter', kk
      YU=0.


      DO 3 KL=1,kxLine
      N=IxD(1,KL)
      I1=IxD(2,KL)
      I2=IxD(3,KL)

      DO 5 M=I1,I2
	dx0=(dlx(m,n)+dlx(m-1,n))/2.
	dy0=(dly(m,n)+dly(m,n-1))/2.
	hgm=4*dx0*dy0
	dx5=dlx(m,n)+dlx(m,n+1)
	dy5=dly(m,n)+dly(m+1,n)
c
c	lvl_uv=(eta(m,n)+eta(m+1,n)+eta(m,n+1)+eta(m+1,n+1))/4.D+0
c	hmt(m,n)=(H1(m,n)-lvl_uv)/(1+dt**2*ctime**2*slb(m,n)**2)

	lvl_uv=(eta(m,n)+eta(m+1,n)+eta(m,n+1)+eta(m+1,n+1))/4.D+0
	hmt(m,n)=(H1(m,n)-lvl_uv)/(1+dt**2*ctime**2*slb(m,n)**2)
     *	*ikbt(m,n,2)

c******** type of boundary for the river inflow ************
	iv2=ikbt(m,n,2)+ikbt(m-1,n,2)
	iv1=ikbt(m,n-1,2)+ikbt(m-1,n-1,2)
	iu2=ikbt(m,n,2)+ikbt(m,n-1,2)
	iu1=ikbt(m-1,n,2)+ikbt(m-1,n-1,2)

	if(rivy(m,n).ne.0.and.rivy(m+1,n).ne.0.and.iv2.eq.0) then
	lvl_uv=(eta(m,n)+eta(m+1,n))/2.D+0
	hh=min(H(m,n),H(m+1,n))
	hmt(m,n)=(hh-lvl_uv)/(1+dt**2*ctime**2*slb(m,n)**2)
	hmt(m,n)=0.
	lvlx(m,n)=(lvl(m+1,n)-lvl(m,n)) *2.   !!!!eto perenesti  v nachalo
	lvly(m,n)=0.
	end if

	if(rivy(m,n).ne.0.and.rivy(m+1,n).ne.0.and.iv1.eq.0) then	
	lvl_uv=(eta(m,n)+eta(m+1,n))/2.D+0
	hh=min(H(m,n),H(m+1,n))
	hmt(m,n-1)=(hh-lvl_uv)/(1+dt**2*ctime**2*slb(m,n)**2)
	hmt(m,n-1)=0.
	lvlx(m,n-1)=(lvl(m+1,n)-lvl(m,n)) *2.
	lvly(m,n-1)=0.
	end if



	if(rivx(m,n).ne.0.and.rivx(m,n+1).ne.0.and.iu2.eq.0) then	
	lvl_uv=(eta(m,n)+eta(m,n+1))/2.D+0
	hh=min(H(m,n),H(m,n+1))
	hmt(m,n)=(hh-lvl_uv)/(1+dt**2*ctime**2*slb(m,n)**2)
	hmt(m,n)=0.
	lvlx(m,n)=0.
	lvly(m,n)=(lvl(m,n+1)-lvl(m,n)) *2.
	end if

	if(rivx(m,n).ne.0.and.rivx(m,n+1).ne.0.and.iu1.eq.0) then		
	lvl_uv=(eta(m,n)+eta(m,n+1))/2.D+0
	hh=min(H(m,n),H(m,n+1))
	hmt(m-1,n)=(hh-lvl_uv)/(1+dt**2*ctime**2*slb(m,n)**2)
	hmt(m-1,n)=0.
	lvlx(m-1,n)=0.
	lvly(m-1,n)=(lvl(m,n+1)-lvl(m,n)) *2.
	end if




	If (m.eq.140.and.n.eq.136) then

	continue
	end if



	rg(m,n)=hmt(m,n)*dt

	q(m,n)=rg(m,n)*ctime*slb(m,n)*dt

	cx(m,n)=rg(m,n)*g*dy5/dx5
	cy(m,n)=rg(m,n)*g*dx5/dy5

	G1(m,n)=(rg(m,n)*grad_x(m,n)+q(m,n)*grad_y(m,n))*dy5
	G2(m,n)=(-q(m,n)*grad_x(m,n)+rg(m,n)*grad_y(m,n))*dx5
	
		AA(m,n)=dy5
     *   *((1-ctime*(1-ctime)*slb(m,n)**2*dt**2)*u_br(m,n)
     *                           +slb(m,n)*dt*v_br(m,n))
	BB(m,n)=dx5
     *   *(-slb(m,n)*dt*u_br(m,n)
     *    +(1-ctime*(1-ctime)*slb(m,n)**2*dt**2)*v_br(m,n))	

	if (m.eq.135.and.n.eq.140) then
	continue
	end if

    5 CONTINUE
	
    3 CONTINUE



      DO 33 KL=1,kxLine
      N=IxD(1,KL)
      I1=IxD(2,KL)
      I2=IxD(3,KL)

	if (n.eq.150) then
	continue
	end if

      DO 55 M=I1,I2

	dx0=(dlx(m,n)+dlx(m-1,n))/2.
	dy0=(dly(m,n)+dly(m,n-1))/2.
	hgm=4*dx0*dy0

	AB=hmt(m,n)*AA(m,n)+hmt(m,n-1)*AA(m,n-1)
     *	 -(hmt(m-1,n-1)*AA(m-1,n-1)+hmt(m-1,n)*AA(m-1,n))+
     *        hmt(m,n)*BB(m,n)+hmt(m-1,n)*BB(m-1,n)
     *	 -(hmt(m-1,n-1)*BB(m-1,n-1)+hmt(m,n-1)*BB(m,n-1))

c	if (AA(m,n).ne.0.and.hmt(m,n).eq.0 ) then
c	print *,m,n
c	read*
c	end if


	GG=G1(m,n)+G1(m,n-1)-G1(m-1,n-1)-G1(m-1,n)+
     *  G2(m-1,n)+G2(m,n)-G2(m-1,n-1)-G2(m,n-1)

	etan=(1-cgam)*(
     *      cx(m,n)*lvlx(m,n)+cx(m,n-1)*lvlx(m,n-1)
     *     -cx(m-1,n-1)*lvlx(m-1,n-1)-cx(m-1,n)*lvlx(m-1,n)
     *     +cy(m-1,n)*lvly(m-1,n)+cy(m,n)*lvly(m,n)
     *     -cy(m-1,n-1)*lvly(m-1,n-1)-cy(m,n-1)*lvly(m-1,n)
     *     +q(m,n-1)*lvly(m,n-1)+q(m,n)*lvly(m,n)
     *	 -q(m-1,n-1)*lvly(m-1,n-1)-q(m-1,n)*lvly(m-1,n)
     *     -q(m-1,n)*lvlx(m-1,n)-q(m,n)*lvlx(m,n)
     *	 +q(m-1,n-1)*lvlx(m-1,n-1)+q(m,n-1)*lvlx(m,n-1))
     *         +hgm/dt*lvl(m,n)

	ct4=cgam*(
     *     cx(m-1,n-1)+cx(m-1,n)-2*g*q(m-1,n-1)+2*g*q(m-1,n)
     *     -cy(m-1,n)-cy(m-1,n-1))


	ct6=cgam*(
     *           cx(m,n-1)+cx(m,n)-cy(m,n)-cy(m,n-1)
     *	      -2*q(m,n)*g+2*q(m,n-1)*g)


	ct5=-cgam*(
     *     cx(m,n)+cx(m-1,n)+cx(m,n-1)+cx(m-1,n-1)
     *     +cy(m,n)+cy(m-1,n)+cy(m,n-1)+cy(m-1,n-1))
	
	
	ct1=cgam*(cx(m-1,n-1)+cy(m-1,n-1))

	ct2=cgam*(-cx(m,n-1)-cx(m-1,n-1)+cy(m,n-1)+cy(m-1,n-1)
     *			-2*q(m,n-1)*g+2*q(m-1,n-1)*g)
	ct3=cgam*(cx(m,n-1)+cy(m,n-1))
	ct7=cgam*(cx(m-1,n)+cy(m-1,n))
	ct8=cgam*(-cx(m,n)-cx(m-1,n)+cy(m,n)+cy(m-1,n)
     *			-2*q(m-1,n)*g+2*q(m,n)*g)
	ct9=cgam*(cx(m,n)+cy(m,n))

	A1(M)=-ct4
	C1(M)=-ct6
	B1(M)=-ct5+hgm/dt

	Fus=0.0 !dobavka vtekaiuschaia iz rek


c******** type of boundary for the river inflow ************
	iv2=ikbt(m,n,2)+ikbt(m-1,n,2)
	iv1=ikbt(m,n-1,2)+ikbt(m-1,n-1,2)
	iu2=ikbt(m,n,2)+ikbt(m,n-1,2)
	iu1=ikbt(m-1,n,2)+ikbt(m-1,n-1,2)
c-------------------------------------------------------------


	
	if (rivx(m,n).ne.0.and.iu2.eq.0) then

	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m+1,n))
     *           *inbt(m,n,2)*inbt(m,n+1,2)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m+1,n-1))
     *           *inbt(m,n,2)*inbt(m,n-1,2)


		hh=min(H(m,n),H(m,n+1))
     *		-0.5*(eta(m,n)+eta(m,n+1))

	Fus=2.* ut*0.5*(dly(m,n)+dly(m+1,n))*hh

		hh=min(H(m,n),H(m,n-1))
     *		-0.5*(eta(m,n)+eta(m,n-1))
	fus=fus
     *      +2.* ub*0.5*(dly(m,n-1)+dly(m+1,n-1))*hh
	end if


		if (rivx(m,n).ne.0.and.iu1.eq.0) then
	ut=-4.D+0/(h(m,n)+h(m,n+1))*
     *           (ftok(M,N+1)-ftok(M,N))
     *           /(dly(m,n)+dly(m-1,n))
     *           *inbt(m,n,2)*inbt(m,n+1,2)
  
	ub=-4.D+0/(h(m,n-1)+h(m,n))*
     *           (ftok(M,N)-ftok(M,N-1))
     *           /(dly(m,n-1)+dly(m-1,n-1))
     *           *inbt(m,n,2)*inbt(m,n-1,2)


		hh=min(H(m,n),H(m,n+1))
     *		-0.5*(eta(m,n)+eta(m,n+1))

	Fus=-2.* ut*0.5*(dly(m,n)+dly(m-1,n))*hh

		hh=min(H(m,n),H(m,n-1))
     *		-0.5*(eta(m,n)+eta(m,n-1))
	fus=fus
     *      -2.* ub*0.5*(dly(m,n-1)+dly(m-1,n-1))*hh
	end if



		if (rivy(m,n).ne.0.and.iv2.eq.0) then  
              vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n+1))
     *           *inbt(m,n,2)*inbt(m+1,n,2)
             vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n+1))
     *           *inbt(m,n,2)*inbt(m-1,n,2)

		hh=min(H(m,n),H(m+1,n))
     *		-0.5*(eta(m,n)+eta(m+1,n))

	Fus=2.* vp*0.5*(dlx(m,n)+dlx(m,n+1))*hh

		hh=min(H(m,n),H(m-1,n))
     *		-0.5*(eta(m,n)+eta(m-1,n))
	fus=fus
     *      +2.*vl*0.5*(dlx(m-1,n)+dlx(m-1,n+1))*hh
	end if


		if (rivy(m,n).ne.0.and.iv1.eq.0.or.n.eq.2) then

              vp=4.D+0/(h(m,n)+h(m+1,n))*
     *          (ftok(M+1,N)-ftok(M,N))
     *          /(dlx(m,n)+dlx(m,n-1))
     *           *inbt(m,n,2)*inbt(m+1,n,2)
             vl=4.D+0/(h(m,n)+h(m-1,n))*
     *          (ftok(M,N)-ftok(M-1,N))
     *          /(dlx(m-1,n)+dlx(m-1,n-1))
     *           *inbt(m,n,2)*inbt(m-1,n,2)

		hh=min(H(m,n),H(m+1,n))
     *		-0.5*(eta(m,n)+eta(m+1,n))

	Fus=-2.* vp*0.5*(dlx(m,n)+dlx(m,n-1))*hh

		hh=min(H(m,n),H(m-1,n))
     *		-0.5*(eta(m,n)+eta(m-1,n))
	fus=fus
     *      -2.*vl*0.5*(dlx(m-1,n)+dlx(m-1,n-1))*hh

	end if



      FF1(M)=etan+AB+GG
     *  +ct1*eta(m-1,n-1)+ct2*eta(m,n-1)+ct3*eta(m+1,n-1)
     *  +ct7*eta(m-1,n+1)+ct8*eta(m,n+1)+ct9*eta(m+1,n+1)
     *  +Fus


	If (m.eq.68.and.n.eq.181) then
	continue
c		print *,'bering', eta(68,181), eta(69,181),eta(70,181)
c		print *,'bering', eta(68,180), eta(69,180),eta(70,180)
c		print *,'bering', eta(68,179), eta(69,179),eta(70,179)
c	print *

	end if

	

	If (m.eq.69.and.n.eq.181) then
	continue
	end if


		If (m.eq.128.and.n.eq.147.and.mod(kk,10).eq.1) then
	continue
	end if




c	if (rivy(m,n).ne.0.and.rivy(m-1,n).ne.0.and.rivy(m+1,n).eq.0) then 
c	  a1(m)= -hgm
c	  b1(m)=hgm
c	  C1(M)=0
c        FF1(M)= -slb(m-1,n)*(ftok(m,n)-ftok(m-1,n))/H(m,n)/g *hgm
c	 end if


	
c	 if (n.eq.2.and.m.gt.i1) then 
c	  a1(m)= -hgm
c	  b1(m)=hgm
c	  C1(M)=0
c        FF1(M)= -slb(m-1,n)*(ftok(m,n)-ftok(m-1,n))/H(m,n)/g *hgm
c	 end if


c	 if (rivx(m,n).ne.0.and.rivx(m,n-1).ne.0) then 
c	  a1(m)=0
c	  b1(m)=hgm
c	  C1(M)=0
c      FF1(M)=eta(m,n-1)*hgm
c     *    	+slb(m,n-1)*(ftok(m,n)-ftok(m,n-1))/H(m,n)/g *hgm
c	 end if

	continue



   55 CONTINUE
      P2=0.
      P3=0.
	

	if (n.eq.181) then
	continue
	end if


         call FACT1(A1,B1,C1,FF1,U1,p2,p3,i1,i2)
caaa      end if  

	if (n.ne.2) then	
		DO 7 M=I1,I2
		If(dabs(u1(m)-eta(m,n)).gt.yu) yu=dabs(u1(m)-eta(m,n))	
    7		eta(M,N)=U1(M)
	end if

   33 CONTINUE

c      DO 9 M=I1,I2      	
c    9 eta(M,2)=eta(M,3)



      YU1=0.
      DO 8 M=M0,MM
      DO 8 N=N0,NN
    8 If(YU1.lt.dabs(eta(m,n)))YU1=DABS(eta(M,N))

      UR=YU/YU1

      IF(UR.GT.1.E-4)GOTO2

   20 CONTINUE
c** End of iteration
c	print *
c      print *,'kol-vo iter', kk

       do n=1,nh
         do m=1,mh
         lvl(m,n)=eta(m,n)
c     *	   *inbt(m,n,1)

         end do
       end do


         do m=1,mh
         lvl(m,1)= lvl(m,2)
		lvl0(m,1)= lvl0(m,2)
c     *	   *inbt(m,n,1)
         end do
   




108	 continue

	goto 166

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

          lvl(m,n)=lvl0(m,n)+dt* (
     *           	(u4-u3)/dx2
     *            +(v4-v3)/dy2
     *           	)

	if (wer.gt.1.e-2) then
c	print *,'wer',wer,'m',m,'n',n
	end if 

  16  continue

  166  continue


c		 print 116,lvl(60,45),lvl_sum(60,45)
c116    format(' lvl and lvl_sum 'f12.1,x,f12.1)

      RETURN



      END
