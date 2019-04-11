      SUBROUTINE read_slp(slp)
	      Include 'model.par'

      real *8 slp(mh,nh,12)
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      real *8 Z,Dz,Dzk,dx,dy,h,h1
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly,sl,slb
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/d3/dlx(mh,nh),dly(mh,nh)
	real *8 p(mh,nh)
	real *8 tx(mh,nh),ty(mh,nh)
      real *8 xx(mh,nh),yy(mh,nh),zz(kh),xxu(mh,nh),yyu(mh,nh)

      pi=3.14159265
      ro_a=1.22                                                     ! kg/m3 
      Open(11,file='C:\dina\BG\sources\slp-2007-monthly-my.dat'
     *        ,form='binary',status='unknown')
c      Open(11,file='c:/elen/bg-exp/sl_pres_c.dat'
c     *        ,form='binary',status='unknown')
      do kmonth=1,12 
        read(11)p
          do n=1,nh
	     do m=1,mh
	       slp(m,n,kmonth)=p(m,n)
	     end do
	    end do
      end do
        close(11,STATUS='KEEP')


      return
      end 

c---------------------
      SUBROUTINE wind_monthly(slp)
	      Include 'model.par'
      COMMON/DAY/myear,month,nday,nsecs

      real *8 slp(mh,nh,12)
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      real *8 Z,Dz,Dzk,dx,dy,h,h1
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly,sl,slb
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/d3/dlx(mh,nh),dly(mh,nh)
	real *8 p(mh,nh)
	real *8 tx(mh,nh),ty(mh,nh)
      real *8 xx(mh,nh),yy(mh,nh),zz(kh),xxu(mh,nh),yyu(mh,nh)
      real *8 tax,tay
      common/Wind/Tax(mh,nh),Tay(mh,nh)    
      integer *4 monlen(12)
       data om,pi,monlen/0.729e-4,3.14159265,
     *31,28,31,30,31,30,31,31,30,31,30,31/
      nd30=monlen(month)
      ro_a=1.22                                                     ! kg/m3 
        mt1=month
        mt2=month+1
	        if(month.eq.12) mt2=1


      do m=1,mh
        do n=1,nh

         p(m,n)=(slp(m,n,mt1)+(slp(m,n,mt2)-slp(m,n,mt1))/nd30*(nday-1))
         end do
      end do

c     wind in m/sec
      do n=1,nh-1
	do m=1,mh-1
	tay(m,n)=1.e+4/ro_a/slb(m,n)
     *          *(p(m+1,n+1)+p(m+1,n)-p(m,n+1)-p(m,n))
     *         /(dlx(m,n)+dlx(m,n+1))
	tax(m,n)=-1.e+4/ro_a/slb(m,n)
     *          *(p(m+1,n+1)-p(m+1,n)+p(m,n+1)-p(m,n))
     *         /(dly(m,n)+dly(m+1,n))
	 end do
	end do


      do n=1,nh-1
	do m=1,mh-1
      wind=sqrt(tax(m,n)**2+tay(m,n)**2)
	if(wind.lt.15)then
	tx(m,n)=0.7*(cos(pi/6.)*tax(m,n)-sin(pi/6.)*tay(m,n))
	ty(m,n)=0.7*(sin(pi/6.)*tax(m,n)+cos(pi/6.)*tay(m,n))
      else
      tx(m,n)=0.8*(cos(pi/9.)*tax(m,n)-sin(pi/9.)*tay(m,n))
	ty(m,n)=0.8*(sin(pi/9.)*tax(m,n)+cos(pi/9.)*tay(m,n))
	end if
	 end do
	end do
c     wind stress in g/(sm*sec2)
      do n=1,nh-1
	do m=1,mh-1
      wind=sqrt(tx(m,n)**2+ty(m,n)**2)
      Dao=(1.1+0.04*wind)*1.e-3
	tax(m,n)=ro_a*wind*Dao*tx(m,n)*10
	tay(m,n)=ro_a*wind*Dao*ty(m,n)*10
	 end do
	end do

      return
      end 





      SUBROUTINE start_TS(Temp,Sal)
	Include 'model.par'
      real *8 Temp(nodes,kgr),Sal(nodes,kgr)

      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      real *8 Z,Dz,Dzk,dx,dy,h,h1
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
	real *4 p(mh,nh)
      real *8 z_bg(20),t_bg(20),s_bg(20),t_av(kgr),s_av(kgr)
      
	Open (12,file='C:\dina\BG\sources\TS-av.dat'
     *         ,form='formatted',status='unknown')
      read (12,1)
      do k=1,20 
    	   read(12,1)z1,t1,s1,d1
            z_bg(k)=z1
            t_bg(k)=t1
            s_bg(k)=s1
	end do
    1 format(f8.0,3x,f6.3,3x,f6.3,3x,f6.3)
      z_bg=z_bg*100.
      do k=2,kkz
	do m=1,19
	if (z(k).ge.z_bg(m).and.z(k).lt.z_bg(m+1)) then
	t_av(k-1)=t_bg(m)+(t_bg(m+1)-t_bg(m))/(z_bg(m+1)-z_bg(m))
     *          *(z(k)-z_bg(m))
	s_av(k-1)=s_bg(m)+(s_bg(m+1)-s_bg(m))/(z_bg(m+1)-z_bg(m))
     *          *(z(k)-z_bg(m))
	end if
	if (z(k).ge.z_bg(20)) then
	t_av(k-1)=t_bg(20)
	s_av(k-1)=s_bg(20)
      end if
	end do
	end do
	           do k=1,kgr  
	            do n=1,nh
	             do m=1,mh
	temp(kp(m,n),k)=1000.
	                 if(kbt(m,n).ge.k+1) then
	                    temp(kp(m,n),k)= 0*1+t_av(k)
						sal(kp(m,n),k)=0*34+s_av(k)
	                 end if
                   end do
	            end do
	           end do
      Print *,(Temp(kp(50,2),k),k=1,kgr)
	      Print *,(Sal(kp(50,2),k),k=1,kgr)
	k2=kp(101,93)
c	print *,temp(k2,5),sal(k2,5)
	read *


	return 
	end

c---------------------
      SUBROUTINE wind_bg
	      Include 'model.par'

      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      real *8 Z,Dz,Dzk,dx,dy,h,h1
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly,sl,slb
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/d3/dlx(mh,nh),dly(mh,nh)
	real *4 p(mh,nh)
	real *8 tx(mh,nh),ty(mh,nh)
      real *8 xx(mh,nh),yy(mh,nh),zz(kh),xxu(mh,nh),yyu(mh,nh)
      real *8 tax,tay
      common/Wind/Tax(mh,nh),Tay(mh,nh)    

      pi=3.14159265
      ro_a=1.22 
	                                                 ! kg/m3 
      Open(11,file='C:\dina\BG\sources\sl_pres_id.dat'
     *        ,form='binary',status='unknown')

c	     Open(11,file='C:\elena\BG-exp\sources\slp-2007-monthly-my.dat'
c     *        ,form='binary',status='unknown')


c      Open(11,file='C:\elena\BG-exp\sources\sl_pres_c.dat'
c     *        ,form='binary',status='unknown')

      read(11)p
      close(11,STATUS='KEEP')

c     wind in m/sec
      do n=1,nh-1
	do m=1,mh-1

	if (m.eq.74.and.n.eq.1.) then
	continue
	end if

	tay(m,n)=1.e+4/ro_a/slb(m,n)
     *          *(p(m+1,n+1)+p(m+1,n)-p(m,n+1)-p(m,n))
     *         /(dlx(m,n)+dlx(m,n+1))
	tax(m,n)=-1.e+4/ro_a/slb(m,n)
     *          *(p(m+1,n+1)-p(m+1,n)+p(m,n+1)-p(m,n))
     *         /(dly(m,n)+dly(m+1,n))

		if (isNaN(tay(m,n))) tay(m,n)=0
		if (isNaN(tax(m,n))) tax(m,n)=0

	 end do
	end do


      do n=1,nh-1
	do m=1,mh-1
      wind=sqrt(tax(m,n)**2+tay(m,n)**2)
c	print *,m,n,wind
	if(wind.lt.15)then
	tx(m,n)=0.7*(cos(pi/6.)*tax(m,n)-sin(pi/6.)*tay(m,n))
	ty(m,n)=0.7*(sin(pi/6.)*tax(m,n)+cos(pi/6.)*tay(m,n))
      else
      tx(m,n)=0.8*(cos(pi/9.)*tax(m,n)-sin(pi/9.)*tay(m,n))
	ty(m,n)=0.8*(sin(pi/9.)*tax(m,n)+cos(pi/9.)*tay(m,n))
	end if
	 end do
	end do
c     wind stress in g/(sm*sec2)


      do n=1,nh-1
	do m=1,mh-1
      wind=sqrt(tx(m,n)**2+ty(m,n)**2)
      Dao=(1.1+0.04*wind)*1.e-3
	tax(m,n)=ro_a*wind*Dao*tx(m,n)*10
	tay(m,n)=ro_a*wind*Dao*ty(m,n)*10
	 end do
	end do

      return
      end 

