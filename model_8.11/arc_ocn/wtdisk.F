       SUBROUTINE WTDISK(Temp,Sal,UU,VV,WW,Ftok,lvl,lvl0)
C*-------------------------------------------
      Include 'model.par'
      COMMON/DAY/myear,month,nday
      common/Wind/Tax(mh,nh),Tay(mh,nh)    
      common/flux/tflux(mh,nh),sflux(mh,nh)
      integer maxmon(12), myr, mmo, mdy
      Character*80 rpath
      character rname*21
	common/cwz/   cwz00(mh,nh),	cwzvel00(mh,nh)
      data maxmon/31,28,31,30,31,30,31,31,30,31,30,31/
      real *8 Temp(nodes,kgr),Sal(nodes,kgr)
     *       ,UU(nodes,kgr),VV(nodes,kgr),WW(nodes,kgr),Ftok(mh,nh)
     *    ,lvl(mh,nh),lvl0(mh,nh),cwz00

      rname='restartYYYY-MM-DD.ocn'

      rpath='c:\dina\output\'
      Open(15,file=rpath(1:len_trim(rpath))//'nTTT_8xmas.dat'
     *,form='binary'
     *,status='unknown')
            len=nodes*kgr*8
            Write(15)((Temp(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file=rpath(1:len_trim(rpath))//'nSSS_8xmas.dat'
     *,form='binary'


     *,status='unknown')
            Write(15)((Sal(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file=rpath(1:len_trim(rpath))//'nUUU_8xmas.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((UU(i,k),i=1,nodes),k=1,kgr)

      CLOSE(unit=15,status='keep')

      Open(15,file=rpath(1:len_trim(rpath))//'nVVV_8xmas.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((VV(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')

      Open(15,file=rpath(1:len_trim(rpath))//'nWWW_8xmas.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((WW(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')



       Open(15,file=rpath(1:len_trim(rpath))//'nlvl_8xmas.dat'
     * 	,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((lvl(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')

	
       Open(15,file=rpath(1:len_trim(rpath))//'nlvl0_8xmas.dat'
     * 	,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((lvl0(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')

	
       Open(15,file=rpath(1:len_trim(rpath))//'cwz00_8xmas.dat'
     * 	,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((cwz00(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')

      Open(15,file=rpath(1:len_trim(rpath))//'nSTREAM_8xmas.dat'
     * 	,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((Ftok(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')

      RETURN
      END


	 SUBROUTINE Wt_month(Temp,Sal,UU,VV,WW,Ftok,lvl,lvl0)
C*-------------------------------------------
      Include 'model.par'
      COMMON/DAY/myear,month,nday
      common/Wind/Tax(mh,nh),Tay(mh,nh)    
      common/flux/tflux(mh,nh),sflux(mh,nh)
      integer maxmon(12), myr, mmo, mdy
      Character*80 rpath
      character rname*21
	
	common/cwz/   cwz00(mh,nh),	cwzvel00(mh,nh)

      data maxmon/31,28,31,30,31,30,31,31,30,31,30,31/
      real *8 Temp(nodes,kgr),Sal(nodes,kgr)
     *       ,UU(nodes,kgr),VV(nodes,kgr),WW(nodes,kgr),
     *         Ftok(mh,nh),el(mh,nh),lvl(mh,nh),lvl0(mh,nh)
     *   ,cwz00
      integer yer
	character*4 cyer	
	character*4 kod

      yer=month*100+myear

      Write(cyer,'(i4.4)')yer
c	kod='wobr'
c	kod='xmas'
	kod='ts_f'
c	kod='lv_f'

      if (myear.le.20)then
       rpath='C:\dina\result'
	else
      rpath='C:\dina\result'
	end if


      Open(15,file=rpath(1:len_trim(rpath))
     * //'/TT/TTT'//kod//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            len=nodes*kgr*8
            Write(15)((Temp(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/SS/SSS'//kod//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((Sal(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/UU/UUU'//kod//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((UU(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/VV/VVV'//kod//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((VV(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/WW/WWWh'//kod//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((WW(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')

      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/ftok/ftok'//kod//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((Ftok(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')

      Open(15,file=rpath(1:len_trim(rpath))
     *//'/LVL/lvl'//kod//cyer//'.dat'
     * 	,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((lvl(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')

	      Open(15,file=rpath(1:len_trim(rpath))
     *//'/LVL/lvl0'//kod//cyer//'.dat'
     * 	,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((lvl0(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')

	      Open(15,file=rpath(1:len_trim(rpath))
     *//'/LVL/cwz00'//kod//cyer//'.dat'
     * 	,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((cwz00(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')


      RETURN
      END


        SUBROUTINE Wyear(Temp,Sal,UU,VV,WW,Ftok,lvl)
C*-------------------------------------------
      Include 'model.par'
      COMMON/DAY/myear,month,nday
      Character*80 rpath
      real *8 Temp(nodes,kgr),Sal(nodes,kgr)
     *       ,UU(nodes,kgr),VV(nodes,kgr),Ftok(mh,nh),lvl(mh,nh)
     *       ,WW(nodes,kgr)
	character*2 cyer	
	character*4 kod

      yer=myear

c      Write(cyer,'(i2.2)')yer
	write(cyer,'(i2)')myear
c	kod='wobr'
c	kod='xmas'
	kod='ts_f'
c	kod='lv_f'

      rpath='c:/dina/result/'
      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/TT/TTT'//kod//'-yr-'//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            len=nodes*kgr*8
            Write(15)((Temp(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/SS/SSS'//kod//'-yr-'//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((Sal(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/UU/UUU'//kod//'-yr-'//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((UU(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/VV/VVV'//kod//'-yr-'//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((VV(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')

	
      Open(15,file=rpath(1:len_trim(rpath))
     *//'/WW/WWW'//kod//'-yr-'//cyer//'.dat'
     *,form='binary'
     *,status='unknown')
            Write(15)((WW(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')


      Open(15,file=rpath(1:len_trim(rpath))
     *           	//'/ftok/FTOK'//kod//'-yr-'//cyer//'.dat'
     * 	,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((Ftok(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')

      Open(15,file=rpath(1:len_trim(rpath))
     *	//'/LVL/LVL'//kod//'-yr'//cyer//'.dat'
     * 	,form='binary'
     *,status='unknown')
            len=mh*nh*8
         Write(15)((lvl(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')



       OPEN(UNIT=15,FILE=rpath(1:len_trim(rpath))//'yday.dat'
     *	  ,STATUS= 'UNKNOWN')
       write (15,'(3(2x,i12))')myear
       CLOSE(UNIT=15,STATUS='KEEP')


      RETURN
      END
