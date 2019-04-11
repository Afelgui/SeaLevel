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
	   common/clame/hl1(mh,nh),hl2(mh,nh)
	common/sum/lvl_sum(mh,nh),ubrt_sum(nodes),vbrt_sum(nodes)
      real *8 tflux,sflux,tax,tay
      common/Wind/Tax(mh,nh),Tay(mh,nh)    
      common/flux/tflux(mh,nh),sflux(mh,nh)
      common/restor/Sstar(mh,nh),gsalt
      real *8 Sstar,gsalt
      common/clim_ts/tclim(mh,6,kgr),sclim(mh,6,kgr)
	common/brt/UBrT(nodes), VBrT(nodes)
	common/brkl/Uu_Brkl(nodes,kgr), Vv_brkl(nodes,kgr)
c	common/vel/WWvel(nodes,kgr)
	common/int/Uu_int(nodes), Vv_int(nodes)

      real *8 tclim,sclim


      common/levitus/Tlev(nodes),Slev(nodes)
      real *8 tlev,slev,ubrt,vbrt

      real *8 Temp(nodes,kgr),Temp1(nodes,kgr),Sal(nodes,kgr)
     *,rho(nodes,kgr)
     *        ,UU(nodes,kgr),VV(nodes,kgr),WW1(nodes,kgr),WW(nodes,kgr)
     *                 ,Dzet(mh,nh),Ftok(mh,nh),lvl(mh,nh),lvl0(mh,nh)
     *         ,uu_ts(mh,nh,kgr),vv_ts(mh,nh,kgr)
      real *8 Tyr(nodes,kgr),Syr(nodes,kgr)
     *         ,Uyr(nodes,kgr),Vyr(nodes,kgr)
     *         ,Wyr(nodes,kgr),Ftyr(mh,nh)
     *         ,Txyr(mh,nh),Tyyr(mh,nh)
     *    ,Hflyr(mh,nh),Sflyr(mh,nh),WWvel1(nodes,kgr),WWvel(nodes,kgr)
      real *8 Tmth(nodes,kgr),Smth(nodes,kgr)
     *         ,Umth(nodes,kgr),Vmth(nodes,kgr)
     *         ,Wmth(nodes,kgr),Fmth(mh,nh)
     *         ,Txmth(mh,nh),Tymth(mh,nh)
     *         ,Hflmth(mh,nh),Sflmth(mh,nh)
      REAL *8 Risl2
 
      REAL *8 grad_x,grad_y,tmer(nh)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),slb(mh,nh)
      COMMON/V1/ROZ(mh,nh)
      COMMON/BL1/Ihl(nodes)
      COMMON/V4/Risl2
      common/v5/grad_x(mh,nh),grad_y(mh,nh)

      common/v6/Fui1(mh,nh),Fvi1(mh,nh)
      real *8 ROZ,FANT,Fui1,Fvi1,sl,slb
c      common/north/tarc(16,4,23),sarc(16,4,23)
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
      common/iodiag/nu_diag
      Logical restart, renew
       data om,pi,monlen/0.729e-4,3.14159265,
     *31,28,31,30,31,30,31,31,30,31,30,31/
c**************
         nmyear=1947
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
c       dt=3*3600./1.
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
c                   call nulstep(Temp,Sal,uu,vv,ww,ftok,dzet)



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


	      Open(15,file='c:\dina\output\nUUU_17.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
c            read(15)((uu(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      
      Open(15,file='c:\dina\output\nvvv_17.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
c           read(15)((vv(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')

      Open(15,file='c:\dina\output\nwww_17.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
c            read(15)((ww(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')

      Open(15,file='c:\dina\output\nlvl_17.dat',
     *   form='binary',status='unknown')
            len=mh*nh*8
c            read(15)((lvl(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')


      Open(15,file='c:\dina\output\cwz00_17.dat',
     *   form='binary',status='unknown')
            len=mh*nh*8
c            read(15)((cwz00(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')



	do k=1,kgr
         do i=2,nodes
c           Temp(i,k)=1.
            end do
      end do


	do k=1,kgr
         do i=2,nodes
c            sal(i,k)=1.
            end do
      end do

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   Main Do Loop
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO 100 nyear=1,5
      myear=myear+1
      do m=1,mh
         do n=1,nh
             Ftyr(m,n)=0.
             Txyr(m,n)=0.
             Tyyr(m,n)=0.
             Hflyr(m,n)=0.
             Sflyr(m,n)=0.
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


              if(month.eq.12) then
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

      nd1=max(1,nday0)
      nd30=monlen(month)
c      if(month.eq.2.and.4*int(myear/4).eq.myear)nd30=nd30+1
      do 600 nday=nd1,30
      nday0=0
        mth=month+1
	        if(month.eq.12) mth=1
        do k=1,13
	  triv(k)=(rtrans(k,month)+(rtrans(k,mth)-rtrans(k,month))
     *			/nd30*(nday-1))*1.e+6
	  end do
        do k=1,3
	  tbin(k)=(tbstr(k,month)+(tbstr(k,mth)-tbstr(k,month))
     *			/nd30*(nday-1))
	  sbin(k)=(sbstr(k,month)+(sbstr(k,mth)-sbstr(k,month))
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

c       call rivers(triv,rivx,rivy,Ftok)






       do k=1,kgr
         do i=2,nodes
c            Temp(i,k)=4.
            end do
      end do


	do k=1,kgr
         do i=2,nodes
c          Temp(i,k)=1.
            end do
      end do

	  do k=1,kgr
         do i=2,nodes
c         sal(i,k)=1.
            end do
      end do


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   Start New Day
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


          do lljjkk=1,2
           do lkj=1,4


             call density(rho,Temp,Sal)
		call wT_wV(UU,VV,ww,WW1,wwvel,WWvel1,ftok,rivx,rivy,lvl,lvl0)


             call qst3_UV(uu,vv,uu_ts,vv_ts,wwvel1,ww1,lvl)
			call advx_V1(uu,vv)
             call z_adv1(UU,VV,WW1,WWvel1,Ftok,lvl)





		    call new_ad(rho,UU,VV,WW1,WWvel1,Ftok,lvl)   
             call right1      
	        call  STBET0(lvl,lvl0,Gam,rivx,rivy) 
	        call v_brt(UU,VV,WW1,WWvel1,lvl,lvl0)
			call V_full_new(UU,VV,WW,WWvel,ftok,uu_ts,vv_ts)
          call wT_wV(UU,VV,ww,WW1,wwvel,WWvel1,ftok,rivx,rivy,lvl,lvl0)
            


  486         continue
	
               

        do n=1,nh
        do m=1,mh
         i=kp(m,n)
        tflux(m,n)=0*(temp(i,1)-tlev(i))
     *              *2500./(3*3600*24.)*2
     *                            *(3-kz0)
        sflux(m,n)=0*(sal(i,1)-slev(i))
     *              *2500./(5*3600*24.)*2
     *                              *(3-kz0)
        end do
        end do


	if (nday.eq.13) then
	continue
	end if
	
	call  yflux(Temp,Sal,uu,vv,ww,ww1,tbin,sbin,rivy,ftok)
	call qstx(Temp,Sal,uu,vv,ww,ww1,rivx,rivy,tbin,sbin,Ftok,lvl,lvl0)


		print *,'temp'
	print *,temp(kp(75,151),19:31)
	print *,'sal'
	print *,sal(kp(75,151),19:31)
	print *,'uu'
	print *,uu(kp(75,151),19:31)
	print *,'vv'
	print *,vv(kp(75,151),19:31)
		print *,'ww'
	print *,ww(kp(75,151),19:31)


412    format(8(1x,9F6.2,/))
117    format('temp',f12.4,x,f12.4,x,f12.4)
116    format(i4,' lvl 'f12.1,x,f12.1,x,f7.3,x,f7.3)
167    format(' sal ' , f7.3,x,f7.3,x,f7.3,x,f7.3,x,f7.3,x,f7.3)
118    format(' temp ',f7.3,x,f7.3,x,f7.3,x,f7.3,x,f7.3,x,f7.3)
119   format(' ww ',f8.5,x,f8.5,x,f8.5,x,f8.5,x,f8.5,x,f8.5)
120   format(' uu ',f11.5,x,f11.5,x,f11.5,x,f11.5,x,f11.5,x,f11.5)
121   format(' vv ',f11.5,x,f11.5,x,f11.5,x,f11.5,x,f11.5,x,f11.5)



c		Print *,'   ww',ww(kp(101,47),16)
	
c		Print *,' sal    ' , sal(kp(101,47),16)

c	print *,kbt(101,47),kbt(82,37),kbt(83,37)
c	print *
c	print *,kbt(82,36),kbt(82,37),kbt(82,38)

c		print *
c	print *,kbt(83,38),kbt(81,38),kbt(83,36),kbt(81,36)
c	read *
c
	call  difYZ(Temp,Sal,rivy,tbin,sbin)
      call  difXZ(Temp,Sal,rivx)


c           call bkc(temp,sal,uu,vv)

c--------bkc
      end do
	end do
	 call WTDISK(Temp,Sal,UU,VV,WW,Ftok,lvl,lvl0)
		print *, 'lvl',lvl(70,70) 
	 

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   The End of the Day
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        do k=1,kgr
          do i=2,nodes
            Tmth(i,k)=Tmth(i,k)+Temp(i,k)/float(nd30)
            Smth(i,k)=Smth(i,k)+Sal(i,k)/float(nd30)
            Umth(i,k)=Umth(i,k)+UU(i,k)/float(nd30)
            Vmth(i,k)=Vmth(i,k)+VV(i,k)/float(nd30)
            Wmth(i,k)=Wmth(i,k)+WW1(i,k)/float(nd30)
          end do
        end do
        do m=1,mh
           do n=1,nh
              Fmth(m,n)=Fmth(m,n)+Ftok(m,n)/float(nd30)
              Txmth(m,n)=Txmth(m,n)+tax(m,n)/float(nd30)
              Tymth(m,n)=Tymth(m,n)+tay(m,n)/float(nd30)
              Hflmth(m,n)=Hflmth(m,n)+tflux(m,n)/float(nd30)
              Sflmth(m,n)=Sflmth(m,n)+sflux(m,n)/float(nd30)
           end do
        end do
c*    every day writing
c        call OCNCDF('day',Temp,Sal,UU,VV,WW,Ftok,
c     *     Tax,Tay,Tflux,Sflux)
c	 print 116,lvl(60,45),lvl_sum(60,45),ubrt(2526),vbrt(2526)
  600  continue

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   The End of the Month
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       call WTDISK(Temp,Sal,UU,VV,WW,Ftok,lvl,lvl0)
	print *,'another point'
	print *,'temp'
	print *,temp(kp(100,100),17:31)
	print *,'sal'
	print *,sal(kp(100,100),17:31)
	print *,'uu'
	print *,uu(kp(100,100),17:31)
	print *,'vv'
	print *,vv(kp(100,100),17:31)
	print *,'ww'
	print *,ww(kp(100,100),17:31)
	


c	 call WTDISK(Tmth,Smth,Umth,Vmth,Wmth,Fmth,lvl)
  700    format(20(13i4,/))         
  701  format(20(13f5.1,/))
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
           end do
        end do
c        Write(nu_diag,*)'wtdisk'
c        call WTDISK(Temp,Sal,UU,VV,WW,Ftok)
c        Write(nu_diag,*)'ocncdf'
c        call OCNCDF('mth',Tmth,Smth,Umth,Vmth,Wmth,Fmth,
c     *     Txmth,Tymth,Hflmth,Sflmth)
  500  continue
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   The End of the Year
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c        call OCNCDF('yr',Tyr,Syr,Uyr,Vyr,Wyr,Ftyr,
c     *     Txyr,Tyyr,Hflyr,Sflyr)
	call WTDISK(Temp,Sal,UU,VV,WW1,Ftok,lvl,lvl0)
  100    continue
 236  continue


 235  continue

c        call WTDISK(Temp,Sal,UU,VV,WW,Ftok)
c        call Wyear(Tyr,Syr,Uyr,Vyr,Ftyr)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                   The End of the Main Loop
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       close(41)
       close(42)
       close(43)
       close(44)
c   the end of the cycle
c-----------------------------
  111 CONTINUE
       Write(nu_diag,*)'ok!'
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
     *        ,UU(nodes,kgr),VV(nodes,kgr),WW(nodes,kgr)
     *        ,Dzet(mh,nh),Ftok(mh,nh),WW1(nodes,kgr)
      do k=1,kgr
        do n=1,nodes
           uu(n,k)=0
           vv(n,k)=0
           ww(n,k)=0
			ww1(n,k)=0

        end do
      end do

          do k=1,kgr
            do n=1,nh
              do m=1,mh
                 temp(kp(m,n),k)=1000.
                 sal(kp(m,n),k)=1000.
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
	ru(n)=0.07
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

c      print *, m,n,k,temp(k2,k),sal(k2,k)

c	read *
	         end if

	 
	      end do
	   end do
	end do

	return 
	end