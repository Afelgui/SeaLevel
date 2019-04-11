        SUBROUTINE rivers(triv,rivx,rivy,Ftok)
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
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      real *8 rivx(mh,nh),rivy(mh,nh),triv(14),Ftok(mh,nh),tb
      hmin=5000.
	do n=1,nh
	  do m=1,mh
	   rivx(m,n)=0
           rivy(m,n)=0
	  end do
	end do
	triv(14)=0.8*1.e+12
c	triv(14)=0.
c	triv(1)=0.
      triv=0.0

c-1     Pechora	
	rivy(141,111)=triv(1)/dlx(141,111)/hmin/2.
	rivy(142,111)=rivy(141,111)

c-2     Ob+Pur
	rivy(139,131)=triv(2)/dlx(139,131)/hmin/4.
	rivy(140,131)=rivy(139,131)*2
	rivy(141,131)=rivy(139,131)

c-3     Yenisey
	rivx(140,136)=triv(3)/dly(140,136)/hmin/2.
	rivx(140,137)=rivx(140,136)

c-4     Olenek
c	triv(4)=triv(4)*2000
	rivx(124,161)=triv(4)/dly(124,161)/hmin/2.
	rivx(124,162)=rivx(124,161)
c-5     Yana	
	rivy(114,169)=triv(5)/dlx(114,169)/hmin/2.
	rivy(115,169)=rivy(114,169)
c-6     Indigirka
c	rivx(103,171)=triv(6)/dly(103,171)/hmin/2.
c	rivx(103,172)=rivx(103,171)
	rivy(101,173)=triv(6)/dlx(101,173)/hmin/2.
	rivy(102,173)=rivy(101,173)

c-7     Kolyma
c	triv(7)=triv(7)*2000	
	rivy(94,175)=triv(7)/dlx(94,175)/hmin/2.
	rivy(95,175)=rivy(94,175)
c-8     Mackenzie	
	rivy(42,165)=triv(8)/dlx(42,165)/hmin/2.
	rivy(43,165)=rivy(42,165)
c-9     Sev.Dvina
	rivy(134,102)=triv(9)/dlx(134,102)/hmin/2.
	rivy(135,102)=rivy(134,102)
c-10     Lena

c	triv(10)=triv(10)*200	
c	rivy(118,166)=triv(10)/dlx(118,166)/hmin/4.
c	rivy(119,166)=rivy(118,166)
c	rivx(118,165)=triv(10)/dly(118,165)/hmin/4.
c	rivx(118,166)=rivx(118,165)
	rivx(118,165)=triv(10)/dly(118,165)/hmin/2.
	rivx(118,164)=rivx(118,165)
c-11     Khatanga	
	rivy(125,156)=triv(11)/dlx(125,156)/hmin/2.
	rivy(126,156)=rivy(125,156)
	
c-12     Taimyra
	rivx(128,146)=triv(12)/dly(128,146)/hmin/2.
	rivx(128,147)=rivx(128,146)
c-13     Pyasina	

	rivy(135,141)=triv(13)/dlx(135,141)/hmin/2.
	rivy(136,141)=rivy(135,141)
c-14   Bering Strait
	rivy(68,181)=triv(14)/dlx(68,181)/hmin/4.
	rivy(69,181)=rivy(68,181)*2
	rivy(70,181)=rivy(68,181)
    
      tb=-triv(8)
      do n=165,181
	  do m=43,68
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do
	Ftok(69,181)=tb-triv(14)/2

	tb=tb-triv(14)
      do n=175,181
	  do m=70,94
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do

	tb=tb-triv(7)
c      do n=172,175
c	  do m=95,103
      do n=173,175
	  do m=95,101
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do
      



	tb=tb-triv(6)
	Ftok(102,173)=tb
      do n=161,172
	  do m=101,114
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do

	tb=tb-triv(5)
      do n=165,169
	  do m=115,119
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do

c	Ftok(118,166)=tb-triv(10)/2

	tb=tb-triv(10)
      do n=162,164
	  do m=118,124
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do

	tb=tb-triv(4)
      do n=156,161
	  do m=124,126
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do
	tb=tb-triv(11)
       n=156
	  do m=122,125
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	
      do n=141,155
	  do m=114,128
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do

	tb=tb-triv(12)
      do n=141,146
	  do m=128,135
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do

	tb=tb-triv(13)
      do n=137,141
	  do m=136,140
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do

	tb=tb-triv(3)
      do n=131,136
	  do m=140,141
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do

	Ftok(140,131)=tb-triv(2)/2

	tb=tb-triv(2)
      n=131
	  do m=137,139
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	
      do n=111,130
	  do m=137,145
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do

	tb=tb-triv(1)
      do n=102,111
	  do m=134,141
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do
	tb=tb-triv(9)
	Ftok(134,102)=tb
      do n=2,104
	  do m=82,133
	  if(mask(m,n).eq.1) Ftok(m,n)=tb
	end do
	end do
      do m=60,111
	Ftok(m,2)=Ftok(111,2)/(111-60)*(m-60)
	end do
      kjhgf=0

	
c	print *,'triv',triv
c	read *

      RETURN

      END
