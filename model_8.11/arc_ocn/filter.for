      SUBROUTINE filter(lvl)
	
	Include 'model.par'
      COMMON/NH2/ NH2
      REAL *8 SL,SLB,DX,DY,Z,DZ,DZK,h,h1

      integer *2 ikbt,inbt
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt

      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)     
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
	common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
	
      real *8 lvl(mh,nh),lvl1(mh,nh)
      common/d3/dlx(mh,nh),dly(mh,nh)


	do ll=1,kxline
		m1=ixd(2,ll)
		m2=ixd(3,ll)
		n=ixd(1,ll)
		do m=m1+1,m2-1

		lvl1(m,n)=(
     *		lvl(m,n)*inbt(m,n,2)
     *      +lvl(m-1,n)*inbt(m-1,n,2)+lvl(m+1,n)*inbt(m+1,n,2)
     *      +lvl(m,n-1)*inbt(m,n-1,2)+lvl(m,n+1)*inbt(m,n+1,2)
     *                  )
     */(inbt(m,n,2)
     *+inbt(m-1,n,2)+inbt(m+1,n,2)+inbt(m,n-1,2)+inbt(m,n+1,2))
     		end do
	end do


	do ll=1,kxline
		m1=ixd(2,ll)
		m2=ixd(3,ll)
		n=ixd(1,ll)
		do m=m1+1,m2-1
		lvl(m,n)=lvl1(m,n)  
     		end do
	end do

	return
	end 