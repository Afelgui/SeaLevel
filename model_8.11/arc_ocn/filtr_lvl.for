       SUBROUTINE filtr_lvl(lvl,lvl0,Gam,rivx,rivy,ftok) 
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
     * ,ds1,ds2,ds3,ds4,ds5
      REAL *8 GAM(mh,nh),tok_n(mh),Flh(mh,nh), u_br(mh,nh), v_br(mh,nh)
     *	,etan,rg(mh,nh),q(mh,nh),hmt(mh,nh),AA(mh,nh),BB(mh,nh)
     *,G1(mh,nh),G2(mh,nh),lvlx(mh,nh),lvly(mh,nh),cx(mh,nh),cy(mh,nh)
     * ,lvl_uv,lvl1(mh,nh)
	real *8 lvl2,lvl3,lvl4,lvl5,lvl6,lvl7,lvl8,lvl9
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


	DO 133 KL=1,kxLine
      N=IxD(1,KL)
      I1=IxD(2,KL)
      I2=IxD(3,KL)

      DO 155 M=I1,I2

	if (m.ne.1.and.m.ne.mh.and.n.ne.1.and.n.ne.nh) then
c		print *,'ok'
	ds1=(dlx(m,n)+dlx(m-1,n))*(dly(m,n)+dly(m,n-1))

	if (m.ne.2) then
	ds2=(dlx(m-1,n)+dlx(m-2,n))*(dly(m-1,n)+dly(m-1,n-1))
	else
	ds2=(dlx(m-1,n)+dlx(m-1,n))*(dly(m-1,n)+dly(m-1,n-1))
	end if 

	ds3=(dlx(m+1,n)+dlx(m,n))*(dly(m+1,n)+dly(m+1,n-1))

	if (n.ne.2) then
	ds4=(dlx(m,n-1)+dlx(m-1,n-1))*(dly(m,n-1)+dly(m,n-2))
	else
	ds4=(dlx(m,n-1)+dlx(m-1,n-1))*(dly(m,n-1)+dly(m,n-1))
	end if

	ds5=(dlx(m,n+1)+dlx(m-1,n+1))*(dly(m,n+1)+dly(m,n))

	i1=inbt(m,n,2)
	i2=inbt(m-1,n,2)
	i3=inbt(m+1,n,2)

	i4=inbt(m,n-1,2)
	i5=inbt(m,n+1,2)
	lvl1(m,n)=(40*lvl(m,n)*i1*ds1+lvl(m-1,n)*i2*ds2+lvl(m+1,n)*i3*ds3
     *            	+lvl(m,n-1)*i4*ds4+lvl(m,n+1)*i5*ds5)
     *                 /(40*i1*ds1+i2*ds2+i3*ds3+i4*ds4+i5*ds5)

	end if
155	 continue
133	continue



	do ll=1,kxline
		m1=ixd(2,ll)
		m2=ixd(3,ll)
		n=ixd(1,ll)
		do m=m1,m2
		lvl(m,n)=lvl1(m,n)  
     		end do
	end do


         do m=1,mh
         lvl(m,1)= lvl(m,2)
c		lvl0(m,1)= lvl0(m,2)
c     *	   *inbt(m,n,1)
         end do
   
108	 continue


c		 print 116,lvl(60,45),lvl_sum(60,45)
c116    format(' lvl and lvl_sum 'f12.1,x,f12.1)

      RETURN



      END
