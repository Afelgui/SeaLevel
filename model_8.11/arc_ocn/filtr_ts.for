       SUBROUTINE filtr_ts(lvl,lvl0,temp,Sal) 
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
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)

      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      common/clame/hl1(mh,nh),hl2(mh,nh)
      real *8 Dzet(mh,nh),lvl(mh,nh),eta(mh,nh),lvl0(mh,nh)

	
	real *8 Temp(nodes,kgr),Sal(nodes,kgr),Temp1(nodes,kgr),Sal1(nodes,kgr)


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
   
       REAL *8 Q1,Q2,Q3,Q4,Q5,Q6,DTR
     *	 ,DX0,DX1,DY0,DY1,dx2,dy2,dx4,dy4,dx5,dy5,hl1,hl2,sl,slb
     
      REAL *8 ds1,ds2,ds3,ds4,ds5
      REAL *8 GAM(mh,nh)

      real *8 rivx(mh,nh),rivy(mh,nh),ftok(mh,nh),vbrt,ubrt
     *  ,grad_x,grad_y,lvl_sum,tax,tay,hgm,hh,Fus,vl,vp,ub,ut
	      integer nu_diag
	 integer *2 ikbt,inbt
      common/iodiag/nu_diag
      real *8 cnr,ctime

   

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

	k=2	
	i1=inbt(m,n,k)
	i2=inbt(m-1,n,k)
	i3=inbt(m+1,n,k)
	i4=inbt(m,n-1,k)
	i5=inbt(m,n+1,k)
	temp1(kp(m,n),k-1)=(8*temp(kp(m,n),k-1)*i1*ds1*(dz(2)-lvl(m,n))
     *                  +temp(kp(m-1,n),k-1)*i2*ds2*(dz(2)-lvl(m-1,n))
     *                  +temp(kp(m+1,n),k-1)*i3*ds3*(dz(2)-lvl(m+1,n))
     *                  +temp(kp(m,n-1),k-1)*i4*ds4*(dz(2)-lvl(m,n-1))
     *                  +temp(kp(m,n+1),k-1)*i5*ds5*(dz(2)-lvl(m,n+1)))
     *                 /(
     *                  8*i1*ds1*(dz(2)-lvl(m,n))
     *                  +i2*ds2*(dz(2)-lvl(m-1,n))
     *                  +i3*ds3*(dz(2)-lvl(m+1,n))
     *                  +i4*ds4*(dz(2)-lvl(m,n-1))
     *                  +i5*ds5*(dz(2)-lvl(m,n+1)) 
     *                    )

		sal1(kp(m,n),k-1)=(8*sal(kp(m,n),k-1)*i1*ds1*(dz(2)-lvl(m,n))
     *                  +sal(kp(m-1,n),k-1)*i2*ds2*(dz(2)-lvl(m-1,n))
     *                  +sal(kp(m+1,n),k-1)*i3*ds3*(dz(2)-lvl(m+1,n))
     *                  +sal(kp(m,n-1),k-1)*i4*ds4*(dz(2)-lvl(m,n-1))
     *                  +sal(kp(m,n+1),k-1)*i5*ds5*(dz(2)-lvl(m,n+1)))
     *                 /(
     *                  8*i1*ds1*(dz(2)-lvl(m,n))
     *                  +i2*ds2*(dz(2)-lvl(m-1,n))
     *                  +i3*ds3*(dz(2)-lvl(m+1,n))
     *                  +i4*ds4*(dz(2)-lvl(m,n-1))
     *                  +i5*ds5*(dz(2)-lvl(m,n+1)) 
     *                    )
	if (m.eq.37.and.n.eq.157.and.k.eq.2) then
	continue
	end if

	kzv=min(kbt(m,n),10)

	do k=3,kzv
	i1=inbt(m,n,k)
	i2=inbt(m-1,n,k)
	i3=inbt(m+1,n,k)
	i4=inbt(m,n-1,k)
	i5=inbt(m,n+1,k)
	temp1(kp(m,n),k-1)=(8*temp(kp(m,n),k-1)*i1*ds1
     *                	+temp(kp(m-1,n),k-1)*i2*ds2
     *                    +temp(kp(m+1,n),k-1)*i3*ds3
     *            	    +temp(kp(m,n-1),k-1)*i4*ds4
     *                    +temp(kp(m,n+1),k-1)*i5*ds5)
     *                 /(8*i1*ds1+i2*ds2+i3*ds3+i4*ds4+i5*ds5)

	sal1(kp(m,n),k-1)=(8*sal(kp(m,n),k-1)*i1*ds1
     *                +sal(kp(m-1,n),k-1)*i2*ds2
     *                +sal(kp(m+1,n),k-1)*i3*ds3
     *                +sal(kp(m,n-1),k-1)*i4*ds4
     *                +sal(kp(m,n+1),k-1)*i5*ds5)
     *               /(8*i1*ds1+i2*ds2+i3*ds3+i4*ds4+i5*ds5)

	continue
	end do

	end if
155	 continue
133	continue

	do ll=1,kxline
		m1=max(ixd(2,ll),2)
		m2=min(ixd(3,ll),mh-1)
		n=max(2,min(ixd(1,ll),nh-1))
		do m=m1,m2
	        kzv=min(kbt(m,n),10)
	        do k=2,kzv
			   temp(kp(m,n),k-1)=temp1(kp(m,n),k-1)
			   sal(kp(m,n),k-1)=sal1(kp(m,n),k-1)
	        end do
     		end do
	end do

   
108	 continue



      RETURN



      END
