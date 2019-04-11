      Include 'model.par'

      integer*4 ierr,stop_now,ctrl,monlen(12)
      COMMON/DAY/myear,month,nday,nsecs
      COMMON/NH2/ NH2
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      common/Cfsmag/cfsmag(mh,nh,kgr)
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      real *8 z,dz,dzk,dx,dy,h,h1
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      COMMON/botV/h1(mh,nh)
      COMMON/botTS/h(mh,nh)
c#ifdef TRACERS
	real *8 Ptrs(nodes,kgr,ntrac),Pacif(nodes,kgr)
	COMMON/Trac/Ptrs
c#endif
      real *8 tflux,sflux,tax,tay
      common/Wind/Tax(mh,nh),Tay(mh,nh)    
      common/flux/tflux(mh,nh),sflux(mh,nh)
      common/restor/Sstar(mh,nh),gsalt
      real *8 Sstar,gsalt
      common/clim_ts/tclim(mh,6,kgr),sclim(mh,6,kgr)
      real *8 tclim,sclim
      real *8 slp(mh,nh,12)

	common/sum/lvl_sum(mh,nh),ubrt_sum(nodes),vbrt_sum(nodes)
      common/levitus/Tlev(nodes),Slev(nodes)
      real *8 tlev,slev,ubrt,vbrt
      real *8 Temp(nodes,kgr),Temp1(nodes,kgr),Sal(nodes,kgr)
     *        ,rho(nodes,kgr)
     *        ,UU(nodes,kgr),VV(nodes,kgr),WW1(nodes,kgr),WW(nodes,kgr)
     *        ,Dzet(mh,nh),Ftok(mh,nh),lvl(mh,nh),lvl0(mh,nh)
     *        ,uu_ts(mh,nh,kgr),vv_ts(mh,nh,kgr)
     *        ,el(mh,nh)
      real *8 Tyr(nodes,kgr),Syr(nodes,kgr)
     *         ,Uyr(nodes,kgr),Vyr(nodes,kgr)
     *         ,Wyr(nodes,kgr),Ftyr(mh,nh)
     *         ,Txyr(mh,nh),Tyyr(mh,nh),lvlyr(mh,nh)
     *    ,Hflyr(mh,nh),Sflyr(mh,nh),WWvel1(mh,nh,kgr),WWvel(mh,nh,kgr)
      real *8 Tmth(nodes,kgr),Smth(nodes,kgr)
     *         ,Umth(nodes,kgr),Vmth(nodes,kgr)
     *         ,Wmth(nodes,kgr),Fmth(mh,nh)
     *         ,Txmth(mh,nh),Tymth(mh,nh)
     *         ,Hflmth(mh,nh),Sflmth(mh,nh),elmth(mh,nh),lvlmth(mh,nh)
      REAL *8 Risl2,Ke
      real *8 ROZ,FANT,Fui,Fvi,sl,slb
      real *8 dz1, uu_sum,lvl_uv
 
      REAL *8 grad_x,grad_y,tmer(nh)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),slb(mh,nh)
      COMMON/V1/ROZ(mh,nh)
      COMMON/BL1/Ihl(nodes)
      COMMON/V4/Risl2
      common/v5/grad_x(mh,nh),grad_y(mh,nh)
       common/v6/Fui1(mh,nh),Fvi1(mh,nh)
c      common/north/tarc(16,4,23),sarc(16,4,23)
	common/brt/UBrT(nodes), VBrT(nodes)
      real *8 Tl1(mh,nh),Sl1(mh,nh),Tbstr(3,12),Sbstr(3,12)
     *          ,tl2(mh,nh),sl2(mh,nh),rtrans(13,12)
     *          ,tax1(mh,nh),tax2(mh,nh),tay1(mh,nh),tay2(mh,nh)
     *          ,TmeanZ(kgr),SmeanZ(kgr),kub(nodes) 
     *          ,txflux,tyflux,sxflux,syflux,cwz00  
      Real *8 Gam(mh,nh),cx(mh,nh),cy(mh,nh),
     *       ccx(mh,nh),ccy(mh,nh),ccxy(mh,nh),dcx(mh,nh),dcy(mh,nh)
      real *8 rivx(mh,nh),rivy(mh,nh),triv(14),tbin(3),sbin(3)
     *         ,lvl_sum,ubrt_sum,vbrt_sum,cwz
      COMMON/RUNPARS/nmyear
      common/adv_fluxes/ tyflux(mh,nh,kh),syflux(mh,nh,kh),
     *						txflux(mh,nh,kh),sxflux(mh,nh,kh)
     *						,tzflux(mh,nh,kh),szflux(mh,nh,kh)

	common/cwz/   cwz00(mh,nh),	cwzvel00(mh,nh)

      integer nu_diag
	integer matr1(3,3),matr2(3,3),matr3(3,3)
      common/iodiag/nu_diag
      Logical restart, renew
       data om,pi,monlen/0.729e-4,3.14159265,
     *31,28,31,30,31,30,31,31,30,31,30,31/
c**************
         nmyear=1948
	   myear=0


	   month=1
	   nday=1
	   nsecs=0
         nts=0
         nh00=81
         nu_diag=6
         	 restart=.true.   ! if true, initialize using restart file instead of defaults
         renew=.true.     ! if true, do not initialize date with one stored in restart file



      restart=.true.



       Call Org(restart,renew)

		
       dt=2*3600./1.

           call viscosity
	
                   call nulstep(Temp,Sal,uu,vv,ww,ww1,ftok,dzet)
                   call RDDISK1(Temp,Sal,UU,VV,WW1,Ftok,restart)


      call check(Temp,Sal)


                   month0=month
                   nday0=nday

      Open(15,file='c:\dina\input data\sourses\ttt00.dat.z'
     *    ,form='binary'
     *        ,access='sequential'
     *    ,status='unknown')

      Open(2,file='C:\dina/result/output_text_sdfg.txt'
c     *     ,access='sequential',recordtype='fixed'
     *     ,status='unknown')

          do k=1,kgr  
             read(15)tl2
	           do n=1,nh
                 	do m=1,mh
          	       if(isNaN(tl2(m,n))) tl2(m,n)=10000
Cc	               temp(kp(m,n),k)=tl2(m,n)
                  end do
	           end do
           do n=1,6
            do m=1,mh
	         tclim(m,n,k)=tl2(m,n)
            end do
	     end do
            do m=1,mh
	         temp(kp(m,2),k)=temp(kp(m,3),k)
            end do
           end do
	
	do k2=1,nodes
	ubrt(k2)=0
	vbrt(k2)=0
	end do


      CLOSE(unit=15,status='keep')

      Open(15,file='c:\dina\input data\sourses\sss00.dat.z'
     *    ,form='binary'
     *        ,access='sequential'
     *    ,status='unknown')

           do k=1,kgr  
              read(15)sl2
	            do n=1,nh
	             do m=1,mh
	                if(isNaN(sl2(m,n))) sl2(m,n)=10000
Cc	                sal(kp(m,n),k)=sl2(m,n)
                   end do
	            end do
           do n=1,6
            do m=1,mh
	         sclim(m,n,k)=sl2(m,n)
            end do
	     end do 
            do m=1,mh
	         sal(kp(m,2),k)=sal(kp(m,3),k)
            end do
            end do

	      CLOSE(unit=15,status='keep')

	
      Open(41,file='c:\dina\input data\sourses\tax.dat'
     *	,form='binary'
     *        ,access='sequential'
     *        ,status='unknown')
      Open(42,file='c:\dina\input data\sourses\tay.dat'
     * 	,form='binary'
     *        ,access='sequential'
     *        ,status='unknown')
      Open(43,file='c:\dina\input data\sourses\temp1.dat',
     *    	form='binary'
     *     ,access='sequential'
     *     ,status='unknown')
      Open(44,file='c:\dina\input data\sourses\salt1.dat',
     *    	form='binary'
     *     ,access='sequential'
     *      ,status='unknown')

c	do k=1,6 !chtob nachat s 6
       read(41)((tax2(m,n),m=1,mh),n=1,nh)
         read(42)((tay2(m,n),m=1,mh),n=1,nh)
         read(43)((tl2(m,n),m=1,mh),n=1,nh)
         read(44)((sl2(m,n),m=1,mh),n=1,nh)
c	end do
	   
	    
	do n=1,nh
	do m=1,mh
	if(isNaN(tax2(m,n))) tax2(m,n)=0
	if(isNaN(tay2(m,n))) tay2(m,n)=0
	if(isNaN(tl2(m,n))) tl2(m,n)=1000
	if(isNaN(sl2(m,n))) sl2(m,n)=1000
	end do
	end do

  
      Open(45,file='c:\dina\input data\sourses\river.dat',
     *    	form='binary',status='unknown')
      Open(47,file='c:\dina\input data\sourses\ts-bstr.dat',
     *    	form='binary'
     *     ,access='sequential'
     *     ,status='unknown')
      read(45)rtrans
      close(45)
	read(47)Tbstr,Sbstr
	close(47)
        Write(nu_diag,*)'begin'


      Tfrz=-1.8
      Do n=1,nh
         Do m=1,mh
	    k00=kp(m,n)
	    if(k00.gt.1)then
               Do k=1,kbt(m,n)-1
	          if(Temp(k00,k).lt.Tfrz)then
		     Temp(k00,k)=Tfrz
		  endif
	       Enddo
	    Endif
	 Enddo
      Enddo

c      call OCNCDF('ini',Temp,Sal,UU,VV,WW,Ftok,tax,tay,tflux,sflux)
      myear=myear-1





	goto 154

	     Open(15,file='c:\dina\output\TTT_19.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
            read(15)((temp(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')


	     Open(15,file='c:\dina\output\SSS_19.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
            read(15)((sal(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')


	     Open(15,file='c:\dina\output\UUU_19.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
            read(15)((uu(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      
      Open(15,file='c:\dina\output\vvv_19.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
            read(15)((vv(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')

      Open(15,file='c:\dina\output\www_19.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
            read(15)((ww(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
	
      Open(15,file='c:\dina\output\lvl_19.dat',
     *   form='binary',status='unknown')
            len=mh*nh
            read(15)((lvl(i,j),i=1,mh),j=1,nh)
      CLOSE(unit=15,status='keep')

      Open(15,file='c:\dina\output\lvl0_19.dat',
     *   form='binary',status='unknown')
            len=mh*nh
            read(15)((lvl0(i,j),i=1,mh),j=1,nh)
      CLOSE(unit=15,status='keep')

      Open(15,file='c:\dina\output\cwz00_19.dat',
     *   form='binary',status='unknown')
            len=mh*nh
            read(15)((cwz00(i,j),i=1,mh),j=1,nh)
      CLOSE(unit=15,status='keep')

154	 continue

	go to 155
	
	do k=1,kgr
         do i=2,nodes
           Temp(i,k)=1.
            end do
      end do


	do k=1,kgr
         do i=2,nodes
            sal(i,k)=34.
            end do
      end do


	      do k=1,kgr
            do n=1,nh
              do m=1,mh
c                 temp(kp(m,n),k)=1000.
c                 sal(kp(m,n),k)=1000.
                 if(kbt(m,n).ge.k+1) then
                    temp(kp(m,n),k)=1.
                    sal(kp(m,n),k)=34.
c                    if(n.ge.14) temp(kp(m,n),k)=0
                 end if
                 
              end do
           end do
         end do

155		continue

c	sal(kp(140,137),1)=1.
c	sal(kp(140,136),1)=1.
c	sal(kp(139,137),1)=1.
c	sal(kp(139,136),1)=1.
c	sal(kp(140,135),1)=1.


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   Main Do Loop
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      myear=0
	dt=3600*2
      DO 100 nyear=1,12

       call Wyear(Tyr,Syr,Uyr,Vyr,Wyr,Ftyr,lvlyr)
      if(nyear.gt.20) then
      do m=1,mh
         do n=1,nh
c	tax(m,n)=0
c	tay(m,n)=0
         end do
      end do
      end if
	myear=myear+1
      do m=1,mh
         do n=1,nh
             Ftyr(m,n)=0.
             Txyr(m,n)=0.
             Tyyr(m,n)=0.
             Hflyr(m,n)=0.
             Sflyr(m,n)=0.
             lvlyr(m,n)=0.
	     
         end do
      end do
       do k=1,kgr
         do i=2,nodes
            Tyr(i,k)=0.
            Syr(i,k)=0.
            Uyr(i,k)=0.
            Vyr(i,k)=0.
            Wyr(i,k)=0.
         end do
      end do
        
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   Start New Year
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
	 
	 mn0=max(month0,1)
       do 500 month=mn0,12
        month0=0
        Write(nu_diag,*)'year=',myear,'     month=',month
        write(2,*)'year=',myear,'     month=',month
       do k=1,kgr
         do i=2,nodes
            Tmth(i,k)=0.
            Smth(i,k)=0.
            Umth(i,k)=0.
            Vmth(i,k)=0.
            Wmth(i,k)=0.
         end do
      end do
        do m=1,mh
          do n=1,nh
             Fmth(m,n)=0.
             Txmth(m,n)=0.
             Tymth(m,n)=0.
             Hflmth(m,n)=0.
             Sflmth(m,n)=0.
              elmth(m,n)=0.
          	lvlmth(m,n)=0.
          end do
        end do
            do i=1,nodes
              ihl(i)=0
            end do

              do n=1,nh
                 do m=1,mh
                  tax1(m,n)=tax2(m,n)
                  tay1(m,n)=tay2(m,n)
                  tl1(m,n)=tl2(m,n)
                  sl1(m,n)=sl2(m,n)
                  end do
               end do

c	mnth=mod(month+5,12)  !!! chtoby nachat s 6 mesyaca
	mnth=month



			if(mnth.eq.12) then
!              if(month.eq.12) then
                 rewind(41)
                 rewind(42)
                 rewind(43)
                 rewind(44)
               end if
         read(41)((tax2(m,n),m=1,mh),n=1,nh)
         read(42)((tay2(m,n),m=1,mh),n=1,nh)
         read(43)((tl2(m,n),m=1,mh),n=1,nh)
         read(44)((sl2(m,n),m=1,mh),n=1,nh)

	do n=1,nh
	do m=1,mh
	if(isNaN(tax2(m,n))) tax2(m,n)=0
	if(isNaN(tay2(m,n))) tay2(m,n)=0
	if(isNaN(tl2(m,n))) tl2(m,n)=1000
	if(isNaN(sl2(m,n))) sl2(m,n)=1000
	if(isNaN(sl2(m,n))) sl2(m,n)=1000
	end do
	end do
 601  continue

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   Start New Month
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      nd1=1
 	nd1=max(1,nday0)
      nd30=monlen(month)

      do 600 nday=nd1,nd30
		print *
	print *,nday
		print *

	 Write(nu_diag,*)'year=',myear,'   month=',month
	 write(2,*)'nday=',nday
	      nday0=0



        mth=month+1
	  mth=mnth+1

	        if(mnth.eq.12) mth=1



        do k=1,13
	  triv(k)=(rtrans(k,mnth)+(rtrans(k,mth)-rtrans(k,mnth))
     *			/nd30*(nday-1))*1.e+6
c	  triv(k)=rtrans(k,2)
	  end do
        do k=1,3
	  tbin(k)=(tbstr(k,mnth)+(tbstr(k,mth)-tbstr(k,mnth))
     *			/nd30*(nday-1))
	  sbin(k)=(sbstr(k,mnth)+(sbstr(k,mth)-sbstr(k,mnth))
     *			/nd30*(nday-1))
	  end do
      Write(nu_diag,*)'                  Day ',nday

         do m=1,mh
            do n=1,nh
			i=kp(m,n)
              tax(m,n)=(tax1(m,n)+(tax2(m,n)-tax1(m,n))/nd30*(nday-1))
              tay(m,n)=(tay1(m,n)+(tay2(m,n)-tay1(m,n))/nd30*(nday-1))
              tlev(i)=tl1(m,n)+(tl2(m,n)-tl1(m,n))/nd30*(nday-1)
              slev(i)=sl1(m,n)+(sl2(m,n)-sl1(m,n))/nd30*(nday-1)
              if (tlev(i).le.-1.8) then
                 tax(m,n)=tax(m,n)*0.5
                 tay(m,n)=tay(m,n)*0.5
              end if

         end do
         end do
	
c	tax=-tax
c	tay=-tay

c		triv(13)=0.
       call rivers(triv,rivx,rivy,Ftok)
c	rivx=0.
c	rivy=0.


c	print *,kbt(mtt,ntt)


	mtt=62
		ntt=165



	
c	mtt=136
c		ntt=141

c		mtt=102
c		ntt=173


		k2=kp(mtt,ntt)
	kkk=2
c
c	print *,'river', rivx(mtt,ntt)


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   Start New Day
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          do lljjkk=1,12
	

c
c	tax=0.
c	tay=0.
c	print *,'lljjkk=',lljjkk

	if (lljjkk.eq.6.and.nday.eq.2) then
	continue
	end if

  !==close Pyasina=

             call density(rho,Temp,Sal)

c	  call wT_wV(UU,VV,ww,WW1,wwvel,WWvel1,ftok,rivx,rivy,lvl,lvl0)


  	write(2,*), 'before qst',uu(kp(mtt,ntt),1), vv(kp(mtt,ntt),1)

c	
c             call qst3_UV(uu,vv,uu_ts,vv_ts,wwvel1,ww1,lvl)
	      call qst_UV(uu,vv,uu_ts,vv_ts,wwvel1,ww1,lvl,lvl0)
	
		write(2,*), 'after qst' , uu(kp(mtt,ntt),1), vv(kp(mtt,ntt),1)

c	print *,'qst',uu(k2,1),vv(k2,1),ubrt(k2),vbrt(k2)
c	read *
		   call advx_V1(uu,vv)

c	print *,'xy vyaz',uu(k2,1),vv(k2,1),ubrt(k2),vbrt(k2)
c	read *
c	write(2,*), 'before z_adv',uu(kp(mtt,ntt),1), vv(kp(mtt,ntt),1)
             call z_adv1(UU,VV,WW1,WWvel1,Ftok,lvl)


c	print *,'vyaz',uu(k2,1),vv(k2,1),ubrt(k2),vbrt(k2)
c	read *
	write(2,*), 'after z_adv'
	write(2,*), uu(kp(mtt,ntt),1), vv(kp(mtt,ntt),1)
		    call new_ad1(rho,UU,VV,WW1,WWvel1,Ftok,lvl)  

c

c	print *,'nad',uu(k2,1),vv(k2,1),ubrt(k2),vbrt(k2)
c	read *	

c             call right1  

	  call  STBET01(lvl,lvl0,Gam,rivx,rivy,ftok)  !yug kak reka
c	  call  STBET0(lvl,lvl0,Gam,rivx,rivy,ftok)  !yug svobodny

        call v_brt(UU,VV,WW1,WWvel1,lvl,lvl0)

	
c	print *,'before',lvl(mtt-1,ntt),lvl(mtt,ntt),lvl(mtt+1,ntt)

c=============filtration is ready only for implicit scheme   cgam=1======================
c
c        call filtr_lvl(lvl,lvl0,Gam,rivx,rivy,ftok) 
c	   call v_brt(UU,VV,WW1,WWvel1,lvl,lvl0)
c=============end filtration=============================================================

c	print *,'after',lvl(mtt-1,ntt),lvl(mtt,ntt),lvl(mtt+1,ntt)
c	print *


c	read *
	call V_full_new(UU,VV,WW,WWvel,ftok,rivx,rivy,uu_ts,vv_ts)



c	print *,'vfull', uu(k2,1),vv(k2,1),ubrt(k2),vbrt(k2)
c	read *

	wwvel=0.

c      call wT_wV(UU,VV,ww,WW1,wwvel,WWvel1,ftok,rivx,rivy,lvl,lvl0)
	  call 	wT_wV1(UU,VV,ww,WW1,wwvel,WWvel1,uu_ts,vv_ts,
     * ftok,rivx,rivy,lvl,lvl0)
                     
c	 goto 666
        do n=1,nh
        do m=1,mh
         i=kp(m,n)
        tflux(m,n)=0.*(temp(i,1)-tlev(i))
     *              *2500./(3*3600*24.)*2
     *                            *(3-kz0)
        sflux(m,n)=0.*(sal(i,1)-slev(i))
     *              *2500./(5*3600*24.)*2
     *                              *(3-kz0)
        end do
        end do

csmag      call smag(uu,vv)

c	write(2,*) 'temp' 
c	write(2,*) , temp(kp(mtt,ntt),1:33)
c	write(2,*),'sal'
c	write(2,*),sal(kp(mtt,ntt),1:33)
        write(2,*),'  ' 
	write(2,*),'lljjkk=',lljjkk
	   write(2,*),'  ' 
	write(2,*),'uu='
	write(2,*),uu(kp(mtt,ntt),1:6)
	write(2,*),'vv='
	write(2,*),vv(kp(mtt,ntt),1:6)
	write(2,*),'ww='
	write(2,*),ww(kp(mtt,ntt),1:6)

      write(2,*),'  ' 
	
	write(2,*),'ubrt'
	write(2,*),ubrt(kp(mtt,ntt))
	write(2,*),'vbrt'
	write(2,*),vbrt(kp(mtt,ntt))

      write(2,*),'  ' 


	write(2,*),'level='
	write(2,*),lvl(mtt-1,ntt+1),lvl(mtt,ntt+1),lvl(mtt+1,ntt+1)
	write(2,*),lvl(mtt-1,ntt),lvl(mtt,ntt),lvl(mtt+1,ntt)
	write(2,*),lvl(mtt-1,ntt-1),lvl(mtt,ntt-1),lvl(mtt+1,ntt-1)

	  write(2,*),'  '  
	   write(2,*),'  ' 


      call  yflux1(Temp,Sal,uu,vv,ww,ww1,tbin,sbin,rivy,ftok)
      call qstx1(Temp,Sal,uu,vv,ww,ww1,rivx,rivy
     *  ,tbin,sbin,Ftok,lvl,lvl0)


c	read *

	kkk=1

c	read *	


c	print *
	
	do k=1,kgr
         do i=2,nodes
           if (Temp(i,k).lt.0.99999999) then 
c			print *,i,k,Temp(i,k)
c			read *
			end if
         end do
      end do

	
	if(nday.eq.15) then
c	print *,temp(kp(68,181),1),temp(kp(69,181),1),temp(kp(70,181),1)
	continue
	end if

		     call  difYZ(Temp,Sal,rivy,tbin,sbin)
              call  difXZ(Temp,Sal,rivx)
	
c	print *,temp(kp(mtt,ntt),kkk),temp(kp(mtt+1,ntt),kkk)
c	print *,lvl(mtt,ntt),lvl(mtt+1,ntt)

  666  continue

c	print *,temp(kp(68,181),1),temp(kp(69,181),1),temp(kp(70,181),1)
c		print *,sal(kp(60,4),1),sal(kp(61,4),1),sal(kp(62,4),1)

c      print *
c	print *,temp(kp(68,181),1),temp(kp(69,181),1),temp(kp(70,181),1)
c	print *
c	print *
c	print *

c           call bkc(temp,sal,uu,vv)


	  call 	wT_wV1(UU,VV,ww,WW1,wwvel,WWvel1,uu_ts,vv_ts,
     * ftok,rivx,rivy,lvl,lvl0)

c      call wT_wV(UU,VV,ww,WW1,wwvel,WWvel1,ftok,rivx,rivy,lvl,lvl0)


   	 call WTDISK(Temp,Sal,UU,VV,WW,Ftok,lvl,lvl0)

c	print *,'1 step'

	end do
c	call kin_en(UU,VV,Ke)
c      write(81)myear,month,nday,Ke
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   The End of the Day
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c	      Print *,(Sal(kp(50,2),k),k=1,kgr)	
c        print *,temp(kp(70,181),1)
   	 call WTDISK(Temp,Sal,UU,VV,WW,Ftok,lvl,lvl0)

	if (mod(nday,5).eq.0) call filtr_ts(lvl,lvl0,temp,Sal) 

        do k=1,kgr
          do i=2,nodes
            Tmth(i,k)=Tmth(i,k)+Temp(i,k)/float(nd30)
            Smth(i,k)=Smth(i,k)+Sal(i,k)/float(nd30)
            Umth(i,k)=Umth(i,k)+UU(i,k)/float(nd30)
            Vmth(i,k)=Vmth(i,k)+VV(i,k)/float(nd30)
            Wmth(i,k)=Wmth(i,k)+WW(i,k)/float(nd30)
          end do
        end do
        do m=1,mh
           do n=1,nh
              lvlmth(m,n)=lvlmth(m,n)+lvl(m,n)/float(nd30)
              Fmth(m,n)=Fmth(m,n)+Ftok(m,n)/float(nd30)
              Txmth(m,n)=Txmth(m,n)+tax(m,n)/float(nd30)
              Tymth(m,n)=Tymth(m,n)+tay(m,n)/float(nd30)
              Hflmth(m,n)=Hflmth(m,n)+tflux(m,n)/float(nd30)
              Sflmth(m,n)=Sflmth(m,n)+sflux(m,n)/float(nd30)
           end do
        end do




  600  continue
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   The End of the Month
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        call Wt_month(Temp,Sal,UU,VV,WW,Ftok,lvl,lvl0)
	print *,temp(kp(68,181),1),temp(kp(69,181),1),temp(kp(70,181),1)

        do k=1,kgr
          do i=2,nodes
            Tyr(i,k)=Tyr(i,k)+Tmth(i,k)/12.
            Syr(i,k)=Syr(i,k)+Smth(i,k)/12.
            Uyr(i,k)=Uyr(i,k)+Umth(i,k)/12.
            Vyr(i,k)=Vyr(i,k)+Vmth(i,k)/12.
            Wyr(i,k)=Wyr(i,k)+Wmth(i,k)/12.
          end do
        end do
        do m=1,mh
           do n=1,nh
              Ftyr(m,n)=Ftyr(m,n)+Fmth(m,n)/12.
              Txyr(m,n)=Txyr(m,n)+Txmth(m,n)/12.
              Tyyr(m,n)=Tyyr(m,n)+Tymth(m,n)/12.
              Hflyr(m,n)=Hflyr(m,n)+Hflmth(m,n)/12.
              Sflyr(m,n)=Sflyr(m,n)+Sflmth(m,n)/12.
              lvlyr(m,n)=lvlyr(m,n)+lvlmth(m,n)/12.
           end do
        end do

c        call OCNCDF('mth',Tmth,Smth,Umth,Vmth,Wmth,Fmth,
c     *     Txmth,Tymth,Hflmth,Sflmth)
c        call OCNCDF('mth',Temp,Sal,UU,VV,WW,Ftok,
c     *     Txmth,Tymth,Hflmth,Sflmth)
  500  continue

       call Wyear(Temp,Sal,UU,VV,WW,Ftok,lvl)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   The End of the Year
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c        call WTDISK(Temp,Sal,UU,VV,WW,Ftok,lvl,lvl0)
c        call OCNCDF('yr',Tyr,Syr,Uyr,Vyr,Wyr,Ftyr,
c     *     Txyr,Tyyr,Hflyr,Sflyr)
  100    continue
 236  continue


 235  continue

c 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   The End of the Main Loop
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  111 CONTINUE
       Write(nu_diag,*)'ok!'
       close(81)
    5 format(16f5.2)
  200 FORMAT(1X,'MODEL TIME=', i7,'YEARS  ',I6,'MONTHS   ', 'DAY=',
     *I3)
  140 FORMAT(1X,I10,E10.2)
  203 FORMAT(1X,'TAX:',I2,1X,I2,5X,'TAY:',I2,1X,I2)
  204 FORMAT(1X,11f6.2)
      if(nu_diag.ne.6)Close(nu_diag)
      STOP
      END

      subroutine adjust(temp,sal)
      include 'model.par'
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      real *8 z,dz,dzk,dx,dy
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh),dy(nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      real *8 temp(nodes,kgr),sal(nodes,kgr),
     *        ta,sa,tb,sb,ro1,ro2       
      COMMON/BL1/Ihl(nodes)
      integer nu_diag
      common/iodiag/nu_diag
c       dimension ihl(nodes)
      DO  Kl=1,kxLine
         n=IxD(1,Kl)
         I1=IxD(2,Kl)
         I2=IxD(3,Kl)
           DO m=I1,I2
             k2=kp(m,n)
             kzv=kbt(m,n)-1
             do k=2,kzv
      ta=temp(k2,k-1)
      tb=temp(k2,k)
      sa=sal(k2,k-1)
      sb=sal(k2,k)
      ro1=dens(ta,sa,0.)
      ro2=dens(tb,sb,0.)
c      if(k2.eq.256) then
c       Write(nu_diag,*)ro1,ta,sa
c      read *
c      end if
c      ro1=ro_ts(ta,sa)
c      ro2=ro_ts(tb,sb)
      if(ro1.gt.ro2) then
       ta=0.5*(temp(k2,k)+temp(k2,k-1))
       sa=0.5*(sal(k2,k)+sal(k2,k-1)) 
      temp(k2,k)=ta
      temp(k2,k-1)=ta
      sal(k2,k)=sa
      sal(k2,k-1)=sa
      if(ihl(k2).lt.k) ihl(k2)=k
      end if
      end do
      end do
      end do 
      return
      end


      
      subroutine nulstep(temp,sal,uu,vv,ww,ww1,ftok,dzet)
      Include 'model.par'
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      common /Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      common /Gr4/kbt(mh,nh),kbtv(mh,nh)
      real *8 Temp(nodes,kgr),Sal(nodes,kgr)
     *        ,UU(nodes,kgr),VV(nodes,kgr),WW(nodes,kgr),WW1(nodes,kgr)
     *        ,Dzet(mh,nh),Ftok(mh,nh)
c#ifdef TRACERS
	real *8 Ptrs(nodes,kgr,ntrac)
	COMMON/Trac/Ptrs
c#endif
      do k=1,kgr
        do n=1,nodes
           uu(n,k)=0
           vv(n,k)=0
           ww(n,k)=0
           ww1(n,k)=0
        end do
      end do
c        do k=2,kgr
c        temp(k2,k)=(tax(m,n)-1.)*exp(-z(k+1)*h(m,n)/50000.)+1.
c         end do

          do k=1,kgr
            do n=1,nh
              do m=1,mh
                 temp(kp(m,n),k)=1000.
                 sal(kp(m,n),k)=1000.
c#ifdef TRACERS
                 do nt=1,ntrac
                    Ptrs(kp(m,n),k,nt)=0.0
                 end do
c#endif
                 if(kbt(m,n).ge.k+1) then
                    temp(kp(m,n),k)=1.
                    sal(kp(m,n),k)=34.6
c                    if(n.ge.14) temp(kp(m,n),k)=0
                 end if
                 
              end do
           end do
         end do   
         do n=1,nh
           do m=1,mh
              dzet(m,n)=0
              ftok(m,n)=0.
           end do
         end do
      return
      end         

      subroutine viscosity
      Include 'model.par'
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      common/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      real *8 z,dz,dzk,dx,dy
      do  N=1,nh
          AYV(N)=0.5E+7
      end do
      DO  M=1,mh
          AXV(M)=0.5E+7
      end do
c       axv(103)=5.e+8
c       axv(104)=5.e+8
      abh=0.5e+7
      ash=1.e+7
      do k=2,kd
            umx(k)=abh+(ash-abh)*exp(-z(k)/50000.)
            umy(k)=abh+(ash-abh)*exp(-z(k)/50000.)
      end do
       astar=0.8
c       astar=0.1
	 arr=1.05/3.14159265
      alam=4.5e-5
      zstar=250000
      DO 10 N=2,kkz
c         RU(N)=1.
      ru(n)=astar+arr*atan(alam*(z(n)-zstar))
  10  continue
       ru(1)=0
       ru(kh)=0 
	do n=1,kkz
	ru(n)=0.1
	end do
	print *,ru
c       ru(2)=1.
       return
       end 


      subroutine filter2(ftok)
      Include 'model.par'
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      common /Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      common /Gr4/kbt(mh,nh),kbtv(mh,nh)
      real *8 Dzet(mh,nh),Ftok(mh,nh)
	alpha=0.9
        do n=1,nh
	do m=1,mh
	dzet(m,n)=ftok(m,n)
	end do
	end do
      DO  KL=1,kxLine
      N=IxD(1,KL)
      I1=IxD(2,KL)
      I2=IxD(3,KL)
      DO M=I1,I2
       if(mask(m,n).eq.1) then
	   Ftok(m,n)=(alpha*0.5*dzet(m-1,n)+(1.0-alpha)*dzet(m,n)
     *           +alpha*0.5*dzet(m+1,n))
	 end if
      end do
	end do
	do n=1,nh
	do m=1,mh
	dzet(m,n)=ftok(m,n)
	end do
	end do
      DO  KL=1,kyLine
      m=IyD(1,KL)
      I1=IyD(2,KL)
      I2=IyD(3,KL)
      DO n=I1,I2
       if(mask(m,n).eq.1) then
	   Ftok(m,n)=(alpha*0.5*dzet(m,n-1)+(1-alpha)*dzet(m,n)
     *           +alpha*0.5*dzet(m,n+1))
	 end if
      end do
	end do
       return
       end 

      subroutine filter3(a)
      Include 'model.par'
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      common /Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      common /Gr4/kbt(mh,nh),kbtv(mh,nh)
      real *8 a(nodes,kgr),dzet(mh,nh),arr(mh,nh)

      do 100 k=1,kgr

      	do n=1,nh
	     do m=1,mh
	       dzet(m,n)=a(kp(m,n),k)
	       arr(m,n)=dzet(m,n)
	     end do
	    end do
      DO  KL=1,kxLine
        N=IxD(1,KL)
        I1=IxD(2,KL)
        I2=IxD(3,KL)
         DO M=I1,I2
           ki1=kbt(m-1,n)-1
	     ki2=kbt(m,n)-1
	     ki3=kbt(m+1,n)-1 
          if(k.le.ki1.and.k.le.ki2.and.k.le.ki3) then
           arr(m,n)=(0.5*dzet(m-1,n)+3*dzet(m,n)+0.5*dzet(m+1,n))/4.
	    else
          if(k.le.ki1.and.k.le.ki2.and.k.gt.ki3) then
	       arr(m,n)=(0.5*dzet(m-1,n)+1.5*dzet(m,n))/2.
	     end if

          if(k.gt.ki1.and.k.le.ki2.and.k.le.ki3) then
	       arr(m,n)=(0.5*dzet(m+1,n)+1.5*dzet(m,n))/2.
	     end if

	    end if
         end do
	end do

      do n=1,nh
	do m=1,mh
	dzet(m,n)=arr(m,n)
	end do
	end do
      DO  KL=1,kyLine
        m=IyD(1,KL)
        I1=IyD(2,KL)
        I2=IyD(3,KL)
          DO n=I1,I2
            ki1=kbt(m,n-1)-1
	      ki2=kbt(m,n)-1
	      ki3=kbt(m,n+1)-1 
           if(k.le.ki1.and.k.le.ki2.and.k.le.ki3) then
	       arr(m,n)=(0.5*dzet(m,n-1)+3*dzet(m,n)+0.5*dzet(m,n+1))/4.
	    else
          if(k.le.ki1.and.k.le.ki2.and.k.gt.ki3) then
	       arr(m,n)=(1.5*dzet(m,n-1)+0.5*dzet(m,n))/2.
	     end if

          if(k.gt.ki1.and.k.le.ki2.and.k.le.ki3) then
	       arr(m,n)=(0.5*dzet(m,n+1)+1.5*dzet(m,n))/2.
	     end if

	     end if
          end do
	end do

	do n=1,nh
	do m=1,mh
	a(kp(m,n),k)=arr(m,n)
	end do
	end do

  100 continue

       return
       end 
      subroutine check(Temp,Sal)
      Include 'model.par'
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      common/Cfsmag/cfsmag(mh,nh,kgr)
      real *8 z,dz,dzk,dx,dy,h,h1
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      real *8 Temp(nodes,kgr),Sal(nodes,kgr)
	real *8 dx1,dy1,dudx,dudy,dvdx,dvdy,cf(kgr)
	do k=1,kgr
         DO  Ki=1,kxLine
            N=IXD(1,Ki)
            M1=IXD(2,Ki)
            M2=IXD(3,Ki)
            DO M=M1,M2
	         k2=kp(m,n)
	if (Temp(kp(m,n),k).gt.0.and.k.ge.30.and.n.gt.100) then
	temp(kp(m,n),k)=temp(kp(m,n),k-1)
	sal(kp(m,n),k)=sal(kp(m,n),k-1)
	end if
	         if(sal(k2,k).gt.1000.and.inbt(m,n,k+1).eq.1) then
	            i1=inbt(m,n-1,k+1)
	            i3=inbt(m,n+1,k+1)
	            i4=inbt(m-1,n,k+1)
	            i5=inbt(m+1,n,k+1)
	            i6=inbt(m-1,n-1,k+1)
	            i7=inbt(m+1,n+1,k+1)
	            i8=inbt(m+1,n-1,k+1)
	            i9=inbt(m-1,n+1,k+1)
	            if(sal(kp(m,n-1),k).gt.1000.) i1=0
	            if(sal(kp(m,n+1),k).gt.1000.) i3=0
	            if(sal(kp(m-1,n),k).gt.1000.) i4=0
	            if(sal(kp(m+1,n),k).gt.1000.) i5=0
	            if(sal(kp(m-1,n-1),k).gt.1000.) i6=0
	            if(sal(kp(m+1,n+1),k).gt.1000.) i7=0
	            if(sal(kp(m+1,n-1),k).gt.1000.) i8=0
	            if(sal(kp(m-1,n+1),k).gt.1000.) i9=0
	

	temp(k2,k)=	(i1*temp(kp(m,n-1),k)
     *        	+i3*temp(kp(m,n+1),k)
     *        	+i4*temp(kp(m-1,n),k)
     *        	+i5*temp(kp(m+1,n),k)
     *	        +i6*temp(kp(m-1,n-1),k)
     *	        +i7*temp(kp(m+1,n+1),k)
     *	        +i8*temp(kp(m+1,n-1),k)
     *	        +i9*temp(kp(m-1,n+1),k))
     *      /(i1+i3+i4+i5+i6+i7+i8+i9)

	sal(k2,k)=	(i1*sal(kp(m,n-1),k)
     *        	+i3*sal(kp(m,n+1),k)
     *        	+i4*sal(kp(m-1,n),k)
     *        	+i5*sal(kp(m+1,n),k)
     *	        +i6*sal(kp(m-1,n-1),k)
     *	        +i7*sal(kp(m+1,n+1),k)
     *	        +i8*sal(kp(m+1,n-1),k)
     *	        +i9*sal(kp(m-1,n+1),k))
     *      /(i1+i3+i4+i5+i6+i7+i8+i9)

      print *, m,n,k,temp(k2,k),sal(k2,k)

c	read *
	         end if

	 
	      end do
	   end do
	end do

	return 
	end
